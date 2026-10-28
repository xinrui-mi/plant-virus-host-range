#!/usr/bin/env python3
from pathlib import Path
import csv, math, statistics
from collections import defaultdict, Counter

FAMS=["Geminiviridae","Caulimoviridae","Alphaflexiviridae","Betaflexiviridae","Closteroviridae","Mitoviridae","Partitiviridae","Potyviridae","Rhabdoviridae","Secoviridae","Solemoviridae","Tombusviridae"]
E=1e-5

def fval(x):
    try:return float(x)
    except:return float('inf')
def nlog(x):
    return 300.0 if x<=0 else min(300.0,-math.log10(x))
def qtile(v,q):
    v=sorted(v)
    if not v:return None
    if len(v)==1:return float(v[0])
    p=(len(v)-1)*q; a=int(math.floor(p)); b=int(math.ceil(p))
    return float(v[a]) if a==b else v[a]*(b-p)+v[b]*(p-a)

def read_fa(p):
    d={}; h=None; b=[]
    for line in open(p):
        line=line.strip()
        if not line:continue
        if line.startswith('>'):
            if h is not None:d[h]=''.join(b)
            h=line[1:].split()[0]; b=[]
        else:b.append(line)
    if h is not None:d[h]=''.join(b)
    return d

def write_fa(rows,p):
    with open(p,'w') as f:
        for h,s in rows:
            f.write('>'+h+'\n')
            for i in range(0,len(s),80):f.write(s[i:i+80]+'\n')

def read_tsv(p):
    with open(p) as f:
        yield from csv.DictReader(f,delimiter='\t')

def domtbl(p):
    d=defaultdict(list)
    for line in open(p):
        if not line.strip() or line.startswith('#'):continue
        x=line.split()
        if len(x)<22:continue
        r=dict(model=x[0],model_len=int(x[2]),query=x[3],query_len=int(x[5]),full_e=fval(x[6]),i_e=fval(x[12]),hmm_from=int(x[15]),hmm_to=int(x[16]),ali_from=int(x[17]),ali_to=int(x[18]))
        d[r['query']].append(r)
    return d

def best(h, exact=None, prefix=None):
    z=[]
    for x in h:
        if x['i_e']>E:continue
        ok=(exact is not None and x['model'] in exact) or (prefix is not None and any(x['model'].startswith(y) for y in prefix))
        if ok:z.append(x)
    return min(z,key=lambda x:x['i_e']) if z else None

def models(h):return {x['model'] for x in h if x['i_e']<=E}

def seco_motifs(seq,pro,pol):
    plo=max(1,pro['ali_from']-20); phi=min(len(seq),pro['ali_to']+20)
    rlo=max(1,pol['ali_from']-20); rhi=min(len(seq),pol['ali_to']+20)
    target=pro['ali_from']+0.75*(pro['ali_to']-pro['ali_from'])
    c=[]
    for i in range(plo-1,phi-3):
        if seq[i]=='G' and seq[i+2:i+4]=='CG':c.append((0,i+3))
    if not c:
        for i in range(plo-1,phi-1):
            if seq[i:i+2]=='CG':c.append((1,i+1))
    if not c:return None
    cpos=sorted(c,key=lambda x:(x[0],abs(x[1]-target)))[0][1]
    tg=pol['ali_from']+0.78*(pol['ali_to']-pol['ali_from'])
    g=[i+1 for i in range(rlo-1,rhi-2) if seq[i:i+3]=='GDD']
    if not g:return None
    gpos=min(g,key=lambda x:abs(x-tg)); epos=gpos+2
    return (cpos,epos) if epos>cpos else None

def eval_marker(fam,seq,h):
    L=len(seq); ms=models(h)
    if fam=='Geminiviridae':
        a=best(h,exact={'Gemini_AL1'})
        if a and 250<=L<=700:return (nlog(a['i_e']),1,L,'Rep_protein','Gemini_AL1',a['i_e'],'NA')
    elif fam=='Caulimoviridae':
        a=best(h,exact={'RVT_1'}); b=best(h,exact={'RT_RNaseH','RT_RNaseH_2'})
        if a and b and a['ali_from']<b['ali_to']:
            s,e=max(1,a['ali_from']),min(L,b['ali_to'])
            if 180<=e-s+1<=700:return (nlog(a['i_e'])+nlog(b['i_e']),s,e,'Pol_RT_to_RNaseH',a['model']+'+'+b['model'],min(a['i_e'],b['i_e']),'NA')
    elif fam in {'Alphaflexiviridae','Betaflexiviridae'}:
        a=best(h,exact={'Vmethyltransf'}); b=best(h,exact={'Viral_helicase1'}); c=best(h,exact={'RdRP_2'})
        if a and b and c and 700<=L<=3000 and a['ali_from']<b['ali_from']<c['ali_from']:
            return (nlog(a['i_e'])+nlog(b['i_e'])+nlog(c['i_e'])+L/1000,1,L,'entire_replication_protein','Vmethyltransf+Viral_helicase1+RdRP_2',min(a['i_e'],b['i_e'],c['i_e']),'NA')
    elif fam=='Closteroviridae':
        a=best(h,exact={'HSP70'})
        if a and 350<=L<=850:return (nlog(a['i_e']),1,L,'HSP70h','HSP70',a['i_e'],'NA')
    elif fam=='Mitoviridae':
        a=best(h,exact={'Mitovir_RNA_pol'})
        if a and 350<=L<=1200:return (nlog(a['i_e'])+L/2000,1,L,'full_or_near_full_RdRP','Mitovir_RNA_pol',a['i_e'],'NA')
    elif fam=='Partitiviridae':
        a=best(h,prefix={'RdRP_'})
        if a and 300<=L<=1100:return (nlog(a['i_e']),1,L,'complete_RdRP_ORF',a['model'],a['i_e'],'NA')
    elif fam=='Potyviridae':
        a=best(h,exact={'RdRP_1'}); hall={'Peptidase_C6','Potyvirid-P3','Poty_PP','Peptidase_S30','Peptidase_C4','Poty_coat'} & ms
        if a and len(hall)>=3 and 1600<=L<=5000:return (nlog(a['i_e'])+5*len(hall)+L/1500,1,L,'major_polyprotein','RdRP_1+'+','.join(sorted(hall)),a['i_e'],'NA')
    elif fam=='Rhabdoviridae':
        a=best(h,exact={'Mononeg_RNA_pol'}); b=best(h,exact={'Mononeg_mRNAcap'})
        if a and 1200<=L<=3200:return (nlog(a['i_e'])+(nlog(b['i_e']) if b else 0)+L/2000,1,L,'full_or_near_full_L','Mononeg_RNA_pol'+('+Mononeg_mRNAcap' if b else ''),min(a['i_e'],b['i_e']) if b else a['i_e'],'NA')
    elif fam=='Secoviridae':
        a=best(h,exact={'Peptidase_C3G'}); b=best(h,prefix={'RdRP_'})
        if a and b and a['ali_from']<b['ali_from']:
            m=seco_motifs(seq,a,b)
            if m:
                s,e=m
                if 150<=e-s+1<=1200:return (nlog(a['i_e'])+nlog(b['i_e']),s,e,'ProPol_CG_to_GDD','Peptidase_C3G+'+b['model'],min(a['i_e'],b['i_e']),'CG_and_GDD_found')
            return (-1e9,None,None,'ProPol_CG_to_GDD','Peptidase_C3G+'+b['model'],min(a['i_e'],b['i_e']),'MOTIF_FAIL')
    elif fam=='Solemoviridae':
        a=best(h,exact={'RdRP_4'})
        if a and 300<=L<=900:return (nlog(a['i_e']),1,L,'RdRP_ORF','RdRP_4',a['i_e'],'NA')
    elif fam=='Tombusviridae':
        a=best(h,exact={'RdRP_3'})
        if a and 250<=L<=900:return (nlog(a['i_e']),1,L,'polymerase_readthrough_or_frameshift_conserved_ORF','RdRP_3',a['i_e'],'NA')
    return None

def load_master(p):return {r['sequence_id']:r for r in read_tsv(p)}
def load_hosts(p):
    hf=defaultdict(set); ht=defaultdict(set)
    with open(p) as f:
        rd=csv.DictReader(f,delimiter='\t')
        if 'host_family' not in rd.fieldnames:raise RuntimeError('association master lacks host_family')
        for r in rd:
            sid=r['sequence_id']; a=(r.get('host_family') or '').strip(); b=(r.get('host_taxid') or '').strip()
            if a and a.lower() not in {'na','nan','none','unclassified'}:hf[sid].add(a)
            if b and b.lower() not in {'na','nan','none'}:ht[sid].add(b)
    return hf,ht

def main(root):
    root=Path(root).resolve(); base=root/'three_source_manuscript_analysis_v1_1/Fig3_host_phylogenetic_constraint'; inr=base/'pipeline_results/03K_family_specific_markers'; orfd=inr/'01_orf_recovery'; pfd=inr/'02_pfam'; out=inr/'03_final_marker_recovery'; out.mkdir(parents=True,exist_ok=True)
    sm=load_master(root/'final_three_source_pipeline_v1/results/20_three_source_sequence_master_FINAL.tsv'); hf,ht=load_hosts(root/'final_three_source_pipeline_v1/results/20_three_source_association_master_FINAL.tsv')
    plan={r['virus_family']:r for r in read_tsv(base/'pipeline_results/03J_family_specific_marker_audit/01_family_marker_plan.tsv')}
    sums=[]; rejs=[]
    for fam in FAMS:
        seqs=read_fa(orfd/f'{fam}.candidate_orfs.faa'); meta={r['orf_id']:r for r in read_tsv(orfd/f'{fam}.candidate_orfs.tsv')}; hits=domtbl(pfd/f'{fam}.Pfam.domtblout')
        by=defaultdict(list)
        for oid,m in meta.items():by[m['sequence_id']].append(oid)
        sel={}; reject=Counter(); motif_fail=0
        for sid,oids in by.items():
            cand=[]
            for oid in oids:
                if oid not in seqs:continue
                ev=eval_marker(fam,seqs[oid],hits.get(oid,[]))
                if not ev:continue
                if ev[-1]=='MOTIF_FAIL':motif_fail+=1;continue
                cand.append((ev[0],oid,ev))
            if not cand:reject['no_marker_candidate_passing_rule']+=1;continue
            cand.sort(reverse=True,key=lambda x:x[0]); _,oid,ev=cand[0]; s,e=ev[1],ev[2]; sel[sid]=(oid,ev,seqs[oid][s-1:e])
        if motif_fail:reject['secoviridae_expected_domains_but_CG_GDD_motif_fail']+=motif_fail
        rows=[]; fa=[]; lens=[]; insp=set(); recsp=set(); recg=set(); rect=set(); input_sids=set(by)
        for sid in input_sids:
            sp=(sm.get(sid,{}).get('virus_species') or '').strip()
            if sp:insp.add(sp)
        for sid,(oid,ev,mseq) in sel.items():
            m=meta[oid]; x=sm.get(sid,{}); sp=(x.get('virus_species') or '').strip(); ge=(x.get('virus_genus') or '').strip(); tx=(x.get('virus_taxid') or '').strip()
            if sp:recsp.add(sp)
            if ge:recg.add(ge)
            if tx:rect.add(tx)
            lens.append(len(mseq)); fa.append((f'{sid}|marker={ev[3]}',mseq))
            rows.append(dict(virus_family=fam,sequence_id=sid,orf_id=oid,translation_table=m['translation_table'],orf_aa_length=m['aa_length'],marker_type=ev[3],marker_start_aa=ev[1],marker_end_aa=ev[2],marker_length_aa=len(mseq),marker_domain_evidence=ev[4],best_marker_i_evalue=f'{ev[5]:.6g}',motif_status=ev[6],virus_taxid=tx,virus_genus=ge,virus_species=sp,host_annotated='YES' if hf.get(sid) else 'NO',host_family_n=len(hf.get(sid,set()))))
        write_fa(fa,out/f'{fam}.marker.raw.faa')
        fields=['virus_family','sequence_id','orf_id','translation_table','orf_aa_length','marker_type','marker_start_aa','marker_end_aa','marker_length_aa','marker_domain_evidence','best_marker_i_evalue','motif_status','virus_taxid','virus_genus','virus_species','host_annotated','host_family_n']
        with open(out/f'{fam}.marker_map.tsv','w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
        hostann=[sid for sid in sel if hf.get(sid)]; hfu=set(); htu=set()
        for sid in sel:hfu|=hf.get(sid,set());htu|=ht.get(sid,set())
        p=plan.get(fam,{})
        sums.append(dict(virus_family=fam,**{'03J_distinct_species_level_labels':p.get('distinct_virus_species_n','NA')},input_sequence_n=len(input_sids),marker_recovered_sequence_n=len(sel),marker_recovery_pct_sequence=f'{100*len(sel)/len(input_sids):.4f}' if input_sids else 'NA',input_nonmissing_species_label_n=len(insp),marker_recovered_species_label_n=len(recsp),marker_recovery_pct_species_label=f'{100*len(recsp)/len(insp):.4f}' if insp else 'NA',marker_recovered_genus_n=len(recg),marker_recovered_virus_taxid_n=len(rect),marker_length_min=min(lens) if lens else 'NA',marker_length_q05=f'{qtile(lens,0.05):.1f}' if lens else 'NA',marker_length_median=f'{statistics.median(lens):.1f}' if lens else 'NA',marker_length_q95=f'{qtile(lens,0.95):.1f}' if lens else 'NA',marker_length_max=max(lens) if lens else 'NA',marker_host_annotated_sequence_n=len(hostann),marker_host_annotation_pct=f'{100*len(hostann)/len(sel):.4f}' if sel else 'NA',marker_host_taxid_n=len(htu),marker_host_family_n=len(hfu),AA95_representative_n='PENDING_03K_C',alignment_QC_tree_tip_n='PENDING_03K_C',host_annotated_tree_tip_n='PENDING_03K_C',final_eligibility='PENDING_03K_C',final_exclusion_reason='PENDING_03K_C'))
        for reason,n in reject.items():rejs.append(dict(virus_family=fam,reason=reason,n=n))
    sf=list(sums[0].keys())
    with open(out/'03K_B_marker_recovery_summary.tsv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=sf,delimiter='\t');w.writeheader();w.writerows(sums)
    with open(out/'03K_B_rejection_reason_counts.tsv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=['virus_family','reason','n'],delimiter='\t');w.writeheader();w.writerows(rejs)
    rules=[('Geminiviridae','Rep protein'),('Caulimoviridae','Pol RT–RNase H'),('Alphaflexiviridae','entire replication protein'),('Betaflexiviridae','entire replication protein'),('Closteroviridae','HSP70h'),('Mitoviridae','full/near-full RdRP'),('Partitiviridae','complete RdRP'),('Potyviridae','major polyprotein'),('Rhabdoviridae','full/near-full L'),('Secoviridae','Pro–Pol CG→GDD'),('Solemoviridae','RdRP'),('Tombusviridae','conserved polymerase readthrough/frameshift ORF')]
    with open(out/'03K_B_marker_rules.tsv','w',newline='') as f:w=csv.writer(f,delimiter='\t');w.writerow(['virus_family','marker']);w.writerows(rules)
    print('03K-B family-specific marker extraction PASS');print('Output:',out)
    for r in sums:print(r['virus_family'],f"marker={r['marker_recovered_sequence_n']}/{r['input_sequence_n']}",f"species={r['marker_recovered_species_label_n']}/{r['input_nonmissing_species_label_n']}",f"host_families={r['marker_host_family_n']}",f"median_aa={r['marker_length_median']}")

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--project-root',required=True);a=ap.parse_args();main(a.project_root)
