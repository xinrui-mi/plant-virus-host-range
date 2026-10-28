#!/usr/bin/env python3
"""Formal TSWV membership audit.

Identifies the formal current TSWV members present in the released sequence input
(virus TaxID 3052585 or the current species name) and, when the third-party TSWV
source GenBank flat file is available, cross-checks those accessions against it.

The TSWV source GenBank flat file is not redistributed with this release. When it
is unavailable the audit still reports the formal-TSWV membership found in the
released sequence input, and states explicitly that the GenBank cross-check and
the downstream gene extraction could not be performed.
"""
import os,gzip,csv,sys
from pathlib import Path
from Bio import SeqIO
ROOT=Path(os.environ['ANALYSIS_ROOT']); REPO=Path(os.environ.get('REPO_ROOT',ROOT.parent))
seqfile=Path(os.environ.get('SEQ_FINAL', REPO/'data'/'analysis_input'/'sequence_analysis_input.tsv.gz'))
source_gb=Path(os.environ.get('TSWV_SOURCE_GB','') or '')
base=ROOT/'06_TSWV_evolution'; (base/'tables').mkdir(parents=True,exist_ok=True);(base/'qc').mkdir(parents=True,exist_ok=True);(base/'work').mkdir(parents=True,exist_ok=True)
rows=[]
_opener = gzip.open if str(seqfile).endswith('.gz') else open
with _opener(seqfile,'rt',encoding='utf-8-sig') as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
def txt(x): return (x or '').strip()
formal=[r for r in rows if txt(r.get('virus_taxid'))=='3052585' or txt(r.get('virus_species'))=='Orthotospovirus tomatomaculae']
formal_ids={txt(r['sequence_id']) for r in formal if txt(r['sequence_id'])}
if not formal_ids: raise SystemExit('No formal current TSWV members in the released sequence input')
gb_available = bool(str(source_gb)) and source_gb.exists()
if not gb_available:
    print('TSWV source GenBank not available (TSWV_SOURCE_GB unset or missing): '
          'the GenBank cross-check and the formal-TSWV GenBank subset are skipped.')
gb_ids=set(); rec_by_id={}
if gb_available:
    for rec in SeqIO.parse(str(source_gb),'genbank'):
        gb_ids.add(rec.id);rec_by_id[rec.id]=rec
missing=sorted(formal_ids-gb_ids); extra=sorted(gb_ids-formal_ids)
shared=sorted(formal_ids&gb_ids)
# Write fresh formal-TSWV GenBank subset if all records are locally available.
subset=base/'work'/'TSWV_formal_current.gb'
if not missing:
    SeqIO.write([rec_by_id[x] for x in shared],str(subset),'genbank')
else:
    subset.write_text('')
with (base/'tables'/'03_TSWV_formal_membership.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['sequence_id','in_source_GenBank'])
    for x in sorted(formal_ids):w.writerow([x,str(x in gb_ids).upper()])
with (base/'tables'/'05_TSWV_membership_comparison.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['category','sequence_n','sequence_ids'])
    w.writerow(['formal_TSWV',len(formal_ids),';'.join(sorted(formal_ids))])
    w.writerow(['formal_and_sourceGB_shared',len(shared),';'.join(shared)])
    w.writerow(['formal_missing_from_sourceGB',len(missing),';'.join(missing)])
    w.writerow(['sourceGB_not_in_formal_TSWV',len(extra),';'.join(extra)])
gate=[]
if not gb_available: gate.append('TSWV_SOURCE_GB_NOT_AVAILABLE')
elif missing: gate.append('FORMAL_TSWV_ACCESSIONS_MISSING_FROM_LOCAL_SOURCE_GB')
status='GATE_REQUIRED' if gate else 'PASS'
with (base/'qc'/'02_TSWV_membership_GATE.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['status','formal_TSWV_n','sourceGB_shared_n','formal_missing_sourceGB_n','sourceGB_extra_n','gate_reasons','fresh_subset_GB'])
    w.writerow([status,len(formal_ids),len(shared),len(missing),len(extra),';'.join(gate),str(subset)])
print('TSWV membership audit',status,'formal',len(formal_ids),'missingGB',len(missing))
if gate: sys.exit(20)
