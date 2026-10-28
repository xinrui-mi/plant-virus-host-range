suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(parallel)
})

source(file.path(
    Sys.getenv("ANALYSIS_ROOT"),
    "00_config",
    "00_config.R"
))
source(file.path(
    ANALYSIS_ROOT,
    "00_config",
    "00_helpers.R"
))

base <- file.path(
    ANALYSIS_ROOT,
    "03_arthropod_transmission_ecology"
)
make_dirs(base)

###############################################################################
# Configuration
###############################################################################

N_REP <- as.integer(Sys.getenv(
    "FIG4_INTERACTION_SENS_N",
    "100"
))

HOST_KEEP_FRACTION <- as.numeric(Sys.getenv(
    "FIG4_INTERACTION_HOST_KEEP",
    "0.70"
))

VIRUS_KEEP_FRACTION <- as.numeric(Sys.getenv(
    "FIG4_INTERACTION_VIRUS_KEEP",
    "0.70"
))

SENS_CORES <- as.integer(Sys.getenv(
    "FIG4_INTERACTION_SENS_CORES",
    as.character(min(NCORES, 8L))
))

if (!is.finite(N_REP) || N_REP < 1L) {
    stop("N_REP must be >= 1")
}
if (!is.finite(HOST_KEEP_FRACTION) ||
    HOST_KEEP_FRACTION <= 0 ||
    HOST_KEEP_FRACTION > 1) {
    stop("HOST_KEEP_FRACTION must be in (0, 1]")
}
if (!is.finite(VIRUS_KEEP_FRACTION) ||
    VIRUS_KEEP_FRACTION <= 0 ||
    VIRUS_KEEP_FRACTION > 1) {
    stop("VIRUS_KEEP_FRACTION must be in (0, 1]")
}

GROUP_LEVELS <- c(
    "no sequence support",
    "sequence-supported"
)

INTERACTION_TERM <-
    "taxonomy_distance:factor(group)sequence-supported"

###############################################################################
# Reconstruct the exact formal-virus input used in final 04B
###############################################################################

seq <- fread_any(
    SEQ_FINAL,
    sep = "\t",
    showProgress = FALSE
)

assoc <- fread_any(
    ASSOC_FINAL,
    sep = "\t",
    showProgress = FALSE
)

vc <- fread(
    file.path(
        ANALYSIS_ROOT,
        "01_dataset_summary",
        "tables",
        "01_VC80_internal_consistency_QC.tsv"
    ),
    sep = "\t",
    showProgress = FALSE
)

seq <- merge(
    seq,
    vc[, .(
        vc80_final_id,
        vc80_support_qc_pass
    )],
    by = "vc80_final_id",
    all.x = TRUE
)

seq[, arth_supported :=
        toupper(std_text(arthropod_sequence_supported)) == "TRUE"]

seq[, current_family :=
        fifelse(
            current_ictv_family_eligible(
                virus_ictv_match_rank,
                virus_family
            ),
            std_text(virus_family),
            ""
        )]

# Exact sequence-level aggregation used by 04B.
sm <- seq[
    present(virus_taxid),
    .(
        sequence_n =
            uniqueN(sequence_id),

        arthropod_supported_sequence_n =
            uniqueN(sequence_id[arth_supported]),

        arthropod_support_fraction =
            mean(arth_supported),

        any_arthropod_sequence_support =
            any(arth_supported),

        qc_pass_VC80_support =
            uniqueN(
                vc80_final_id[
                    vc80_support_qc_pass %in% TRUE
                ]
            ),

        virus_family =
            mode_value(current_family),

        virus_genus =
            mode_value(virus_genus)
    ),
    by = virus_taxid
]

# Exact association aggregation used by 04B.
am <- assoc[
    present(virus_taxid),
    .(
        association_n = .N,

        host_taxid_n =
            uniqueN(
                host_taxid[present(host_taxid)]
            ),

        host_family_n =
            uniqueN(
                host_family[present(host_family)]
            ),

        host_order_n =
            uniqueN(
                host_order[present(host_order)]
            ),

        country_n =
            uniqueN(
                country[present(country)]
            ),

        continent_n =
            uniqueN(
                continent[present(continent)]
            )
    ),
    by = virus_taxid
]

metric <- merge(
    sm,
    am,
    by = "virus_taxid",
    all.x = TRUE
)

metric[, seq_support_group :=
    fifelse(
        any_arthropod_sequence_support %in% TRUE,
        "sequence-supported",
        "no sequence support"
    )
]

###############################################################################
# Reconstruct exact host taxonomy and host-pair universe used by 04B
###############################################################################

hostcols <- c(
    "host_taxid",
    "host_species",
    "host_genus",
    "host_family",
    "host_order",
    "host_class",
    "host_phylum",
    "host_kingdom"
)

host <- unique(
    assoc[
        complete_host_taxonomy(assoc),
        ..hostcols
    ]
)

if (host[, uniqueN(host_taxid)] != nrow(host)) {
    stop("host taxonomy conflict")
}

host <- host[order(host_taxid)]
host[, hi := .I]

hm <- host[, .(
    host_taxid,
    hi
)]

ij <- combn(
    seq_len(nrow(host)),
    2L
)

pi <- ij[1, ]
pj <- ij[2, ]

dd <- host_distance_vectors(
    host,
    pi,
    pj
)

formal_groups <- metric[, .(
    virus_taxid,
    seq_support_group
)]

e <- unique(
    merge(
        assoc[
            present(virus_taxid) &
            complete_host_taxonomy(assoc),
            .(
                host_taxid,
                virus_taxid
            )
        ],
        formal_groups,
        by = "virus_taxid"
    )
)

e <- merge(
    e,
    hm,
    by = "host_taxid"
)

###############################################################################
# Core function: rebuild the same group-wise sharing data and GLM as final 04B
###############################################################################

fit_interaction <- function(
    edge_dt,
    retained_hosts = NULL,
    retained_viruses = NULL
) {

    interaction_rows <- vector(
        "list",
        length(GROUP_LEVELS)
    )

    diag_rows <- vector(
        "list",
        length(GROUP_LEVELS)
    )

    for (k in seq_along(GROUP_LEVELS)) {

        gp <- GROUP_LEVELS[k]

        ee <- edge_dt[
            seq_support_group == gp
        ]

        if (!is.null(retained_hosts)) {
            ee <- ee[
                hi %in% retained_hosts
            ]
        }

        if (!is.null(retained_viruses)) {
            ee <- ee[
                virus_taxid %in% retained_viruses
            ]
        }

        hh <- sort(unique(ee$hi))
        vl <- sort(unique(ee$virus_taxid))

        if (length(hh) < 2L) {
            stop(
                "Fewer than two hosts in group: ",
                gp
            )
        }

        if (length(vl) < 1L) {
            stop(
                "No virus TaxIDs in group: ",
                gp
            )
        }

        host_present <- logical(nrow(host))
        host_present[hh] <- TRUE

        mask <-
            host_present[pi] &
            host_present[pj]

        if (!any(mask)) {
            stop(
                "No host pairs in group: ",
                gp
            )
        }

        X <- sparseMatrix(
            i = ee$hi,
            j = match(
                ee$virus_taxid,
                vl
            ),
            x = 1,
            dims = c(
                nrow(host),
                length(vl)
            )
        )

        SS <- tcrossprod(X)

        sh <- as.integer(
            as.numeric(
                SS[
                    cbind(
                        pi[mask],
                        pj[mask]
                    )
                ]
            ) > 0
        )

        interaction_rows[[k]] <-
            data.table(
                group = gp,
                taxonomy_distance = dd[mask],
                share = sh
            )[
                ,
                .(
                    total_pairs = .N,
                    positive_pairs =
                        sum(share),
                    sharing_probability =
                        mean(share)
                ),
                by = .(
                    group,
                    taxonomy_distance
                )
            ]

        diag_rows[[k]] <- data.table(
            group = gp,
            host_n = length(hh),
            virus_taxid_n = length(vl),
            host_pair_n = sum(mask),
            positive_pair_n = sum(sh)
        )
    }

    inter <- rbindlist(
        interaction_rows,
        use.names = TRUE
    )

    # Force the same reference group every time.
    inter[, group :=
        factor(
            group,
            levels = GROUP_LEVELS
        )
    ]

    fit <- glm(
        cbind(
            positive_pairs,
            total_pairs - positive_pairs
        ) ~
            taxonomy_distance *
            factor(group),
        family = binomial(),
        data = inter
    )

    cf <- coef(fit)

    if (!all(
        c(
            "taxonomy_distance",
            INTERACTION_TERM
        ) %in% names(cf)
    )) {
        stop(
            "Required coefficients absent from fitted model"
        )
    }

    beta_unsupported <-
        unname(
            cf["taxonomy_distance"]
        )

    delta_beta <-
        unname(
            cf[INTERACTION_TERM]
        )

    beta_sequence_supported <-
        beta_unsupported +
        delta_beta

    d0 <- rbindlist(diag_rows)
    d_no <- d0[
        group == "no sequence support"
    ]
    d_yes <- d0[
        group == "sequence-supported"
    ]

    stats <- data.table(
        beta_no_sequence_support =
            beta_unsupported,

        beta_sequence_supported =
            beta_sequence_supported,

        delta_beta =
            delta_beta,

        host_n_no_sequence_support =
            d_no$host_n,

        host_n_sequence_supported =
            d_yes$host_n,

        virus_taxid_n_no_sequence_support =
            d_no$virus_taxid_n,

        virus_taxid_n_sequence_supported =
            d_yes$virus_taxid_n,

        host_pair_n_no_sequence_support =
            d_no$host_pair_n,

        host_pair_n_sequence_supported =
            d_yes$host_pair_n
    )

    list(
        inter = inter,
        fit = fit,
        stats = stats
    )
}

###############################################################################
# Mandatory audit:
# full-data reconstruction MUST reproduce final frozen Fig3C exactly.
###############################################################################

observed <- fit_interaction(e)

obs_inter <- copy(observed$inter)
obs_inter[, group := as.character(group)]

saved_inter <- fread(
    file.path(
        base,
        "tables",
        "13_host_distance_by_sequence_evidence_group.tsv"
    )
)

obs_inter <- obs_inter[
    order(
        group,
        taxonomy_distance
    )
]

saved_inter <- saved_inter[
    order(
        group,
        taxonomy_distance
    )
]

if (nrow(obs_inter) != nrow(saved_inter)) {
    stop(
        "Observed reconstruction does not match frozen Table 13 row count"
    )
}

key_cols <- c(
    "group",
    "taxonomy_distance",
    "total_pairs",
    "positive_pairs"
)

for (cc in key_cols) {
    if (!identical(
        as.character(obs_inter[[cc]]),
        as.character(saved_inter[[cc]])
    )) {
        stop(
            "Observed reconstruction mismatch in column: ",
            cc
        )
    }
}

prob_diff <- max(
    abs(
        obs_inter$sharing_probability -
        saved_inter$sharing_probability
    ),
    na.rm = TRUE
)

if (
    !is.finite(prob_diff) ||
    prob_diff > 1e-12
) {
    stop(
        "Observed sharing_probability mismatch; max diff = ",
        prob_diff
    )
}

saved_model <- fread(
    file.path(
        base,
        "tables",
        "14_host_distance_interaction_model.tsv"
    )
)

saved_delta <- saved_model[
    term == INTERACTION_TERM,
    estimate
]

if (length(saved_delta) != 1L) {
    stop(
        "Could not identify unique frozen interaction coefficient"
    )
}

observed_delta <-
    observed$stats$delta_beta

delta_diff <-
    abs(
        observed_delta -
        saved_delta
    )

if (
    !is.finite(delta_diff) ||
    delta_diff > 1e-10
) {
    stop(
        "Observed beta reconstruction mismatch; diff = ",
        delta_diff
    )
}

cat("\nFULL-DATA AUDIT PASS\n")
cat(
    "Observed beta, no sequence support = ",
    observed$stats$beta_no_sequence_support,
    "\n",
    sep = ""
)
cat(
    "Observed beta, sequence-supported = ",
    observed$stats$beta_sequence_supported,
    "\n",
    sep = ""
)
cat(
    "Observed delta beta = ",
    observed_delta,
    "\n",
    sep = ""
)
cat(
    "Frozen delta beta = ",
    saved_delta,
    "\n",
    sep = ""
)
cat(
    "Max sharing-probability difference = ",
    prob_diff,
    "\n\n",
    sep = ""
)

###############################################################################
# Active node universes for sensitivity analyses
#
# Match the frozen Fig3 sensitivity framework: sampling is performed from
# nodes that actually occur in the formal association network used here.
###############################################################################

ACTIVE_HOST_IDS <- sort(unique(e$hi))
ACTIVE_VIRUS_IDS <- sort(unique(e$virus_taxid))

###############################################################################
# One host-node subsampling iteration
###############################################################################

run_host_subsample <- function(it) {

    tryCatch({

        set.seed(
            SEED +
            410000L +
            it
        )

        n_keep <- max(
            2L,
            floor(
                HOST_KEEP_FRACTION *
                length(ACTIVE_HOST_IDS)
            )
        )

        retained_hosts <- sort(
            sample(
                ACTIVE_HOST_IDS,
                size = n_keep,
                replace = FALSE
            )
        )

        z <- fit_interaction(
            e,
            retained_hosts =
                retained_hosts
        )

        out <- copy(z$stats)

        out[, `:=`(
            sensitivity =
                "host_node_subsampling",

            iteration =
                it,

            requested_fraction =
                HOST_KEEP_FRACTION,

            retained_global_node_n =
                n_keep,

            status =
                "PASS",

            error_message =
                ""
        )]

        out

    }, error = function(err) {

        data.table(
            sensitivity =
                "host_node_subsampling",

            iteration =
                it,

            requested_fraction =
                HOST_KEEP_FRACTION,

            status =
                "FAIL",

            error_message =
                conditionMessage(err)
        )
    })
}

###############################################################################
# One virus-node subsampling iteration
#
# Match the frozen Fig3 sensitivity framework: randomly retain the requested
# fraction from the complete active formal-virus node set. Sequence-evidence
# group labels are preserved after sampling, and the same interaction model
# is rebuilt from the retained virus-host association network.
###############################################################################

run_virus_subsample <- function(it) {

    tryCatch({

        set.seed(
            SEED +
            420000L +
            it
        )

        n_keep <- max(
            2L,
            floor(
                VIRUS_KEEP_FRACTION *
                length(ACTIVE_VIRUS_IDS)
            )
        )

        retained_viruses <- sort(
            sample(
                ACTIVE_VIRUS_IDS,
                size = n_keep,
                replace = FALSE
            )
        )

        z <- fit_interaction(
            e,
            retained_viruses = retained_viruses
        )

        out <- copy(z$stats)

        out[, `:=`(
            sensitivity =
                "virus_taxid_subsampling",

            iteration =
                it,

            requested_fraction =
                VIRUS_KEEP_FRACTION,

            retained_global_node_n =
                n_keep,

            status =
                "PASS",

            error_message =
                ""
        )]

        out

    }, error = function(err) {

        data.table(
            sensitivity =
                "virus_taxid_subsampling",

            iteration =
                it,

            requested_fraction =
                VIRUS_KEEP_FRACTION,

            status =
                "FAIL",

            error_message =
                conditionMessage(err)
        )
    })
}

###############################################################################
# Reproducible parallel execution
###############################################################################

run_parallel <- function(FUN) {

    ids <- seq_len(N_REP)

    if (
        .Platform$OS.type == "unix" &&
        SENS_CORES > 1L
    ) {

        mclapply(
            ids,
            FUN,
            mc.cores =
                min(
                    SENS_CORES,
                    N_REP
                ),
            mc.preschedule = TRUE,
            mc.set.seed = FALSE
        )

    } else {

        lapply(
            ids,
            FUN
        )
    }
}

cat(
    "Running host-node subsampling: ",
    N_REP,
    " iterations\n",
    sep = ""
)

host_res <- rbindlist(
    run_parallel(
        run_host_subsample
    ),
    fill = TRUE
)

cat(
    "Running virus-TaxID subsampling: ",
    N_REP,
    " iterations\n",
    sep = ""
)

virus_res <- rbindlist(
    run_parallel(
        run_virus_subsample
    ),
    fill = TRUE
)

iterations <- rbindlist(
    list(
        host_res,
        virus_res
    ),
    fill = TRUE
)

setcolorder(
    iterations,
    c(
        "sensitivity",
        "iteration",
        "requested_fraction",
        "status",
        "beta_no_sequence_support",
        "beta_sequence_supported",
        "delta_beta",
        "host_n_no_sequence_support",
        "host_n_sequence_supported",
        "virus_taxid_n_no_sequence_support",
        "virus_taxid_n_sequence_supported",
        "host_pair_n_no_sequence_support",
        "host_pair_n_sequence_supported",
        "retained_global_node_n",
        "error_message"
    )
)

###############################################################################
# Sensitivity summaries
###############################################################################

ok <- iterations[
    status == "PASS"
]

qfun <- function(x, p) {
    as.numeric(
        quantile(
            x,
            probs = p,
            na.rm = TRUE,
            names = FALSE
        )
    )
}

summary_dt <- ok[
    ,
    .(
        requested_iterations =
            N_REP,

        successful_iterations =
            .N,

        observed_beta_no_sequence_support =
            observed$stats$beta_no_sequence_support,

        beta_no_sequence_support_median =
            median(
                beta_no_sequence_support,
                na.rm = TRUE
            ),

        beta_no_sequence_support_q025 =
            qfun(
                beta_no_sequence_support,
                0.025
            ),

        beta_no_sequence_support_q975 =
            qfun(
                beta_no_sequence_support,
                0.975
            ),

        observed_beta_sequence_supported =
            observed$stats$beta_sequence_supported,

        beta_sequence_supported_median =
            median(
                beta_sequence_supported,
                na.rm = TRUE
            ),

        beta_sequence_supported_q025 =
            qfun(
                beta_sequence_supported,
                0.025
            ),

        beta_sequence_supported_q975 =
            qfun(
                beta_sequence_supported,
                0.975
            ),

        observed_delta_beta =
            observed_delta,

        delta_beta_median =
            median(
                delta_beta,
                na.rm = TRUE
            ),

        delta_beta_q025 =
            qfun(
                delta_beta,
                0.025
            ),

        delta_beta_q975 =
            qfun(
                delta_beta,
                0.975
            ),

        delta_beta_positive_fraction =
            mean(
                delta_beta > 0,
                na.rm = TRUE
            )
    ),
    by = sensitivity
]

###############################################################################
# QC
###############################################################################

qc <- iterations[
    ,
    .(
        requested_iterations =
            N_REP,

        successful_iterations =
            sum(status == "PASS"),

        failed_iterations =
            sum(status != "PASS")
    ),
    by = sensitivity
]

qc[, `:=`(
    full_data_table13_row_match =
        TRUE,

    full_data_probability_max_abs_diff =
        prob_diff,

    full_data_interaction_beta =
        observed_delta,

    frozen_interaction_beta =
        saved_delta,

    full_data_interaction_beta_abs_diff =
        delta_diff
)]

###############################################################################
# Write NEW outputs only. Existing final 04B outputs are untouched.
###############################################################################

write_tsv(
    iterations,
    file.path(
        base,
        "tables",
        "17_host_distance_interaction_node_subsampling_iterations.tsv"
    )
)

write_tsv(
    summary_dt,
    file.path(
        base,
        "tables",
        "18_host_distance_interaction_node_subsampling_summary.tsv"
    )
)

write_tsv(
    qc,
    file.path(
        base,
        "qc",
        "02_Fig4_host_distance_interaction_sensitivity_QC.tsv"
    )
)

cat("\n========================================\n")
cat("Fig4 host-distance interaction sensitivity PASS\n")
cat("========================================\n\n")

print(summary_dt)

cat("\nQC\n")
print(qc)

if (any(qc$failed_iterations > 0L)) {
    warning(
        "One or more sensitivity iterations failed. ",
        "Inspect the iteration table before interpretation."
    )
}

flush.console()
invisible(gc())

if (!interactive()) {
    quit(
        save = "no",
        status = 0,
        runLast = FALSE
    )
}
