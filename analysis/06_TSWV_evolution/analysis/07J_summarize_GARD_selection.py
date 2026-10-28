#!/usr/bin/env python3
from __future__ import annotations
import csv, json, os, re, sys
from pathlib import Path

GENES=['RdRp','GnGc','NSm','N','NSs']


def root_from_env_or_file():
    if os.environ.get('ANALYSIS_ROOT'): return Path(os.environ['ANALYSIS_ROOT']).resolve()
    return Path(__file__).resolve().parents[2]


def write_tsv(p,rows,fields):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def read_tsv(p):
    with p.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))


def walk(o,path=''):
    if isinstance(o,dict):
        for k,v in o.items(): yield from walk(v,path+'/'+str(k))
    elif isinstance(o,list):
        for i,v in enumerate(o): yield from walk(v,path+'/'+str(i))
    else: yield path,o


def header_name(h):
    if isinstance(h,list) and h:return str(h[0])
    return str(h)


def tables(o,path=''):
    if isinstance(o,dict):
        if isinstance(o.get('headers'),list) and isinstance(o.get('content'),dict):
            yield path,[header_name(x) for x in o['headers']],o['content']
        for k,v in o.items(): yield from tables(v,path+'/'+str(k))
    elif isinstance(o,list):
        for i,v in enumerate(o): yield from tables(v,path+'/'+str(i))


def num(x):
    try:return float(x)
    except:return None


def overall_p(js):
    hits=[]
    for p,v in walk(js):
        leaf=p.rsplit('/',1)[-1].lower()
        if 'p-value' in leaf or leaf in {'pvalue','p value'}:
            z=num(v)
            if z is not None and 0<=z<=1:hits.append((p,z))
    preferred=[x for x in hits if 'test results' in x[0].lower() or 'test_results' in x[0].lower()]
    pool=preferred or hits
    return min(pool,key=lambda x:(len(x[0]),x[0])) if pool else ('',None)


def partition_order(content):
    keys=list(content.keys())
    def kf(k):
        m=re.search(r'(\d+)',str(k))
        return (0,int(m.group(1))) if m else (1,str(k))
    return sorted(keys,key=kf)


def parse_meme(js,gene,partmap):
    out=[]
    for path,hdr,content in tables(js):
        pidx=[i for i,h in enumerate(hdr) if 'p-value' in h.lower() or 'p value' in h.lower()]
        if not pidx:continue
        pi=pidx[0]
        for rank,part in enumerate(partition_order(content),1):
            rows=content[part]
            if not isinstance(rows,list):continue
            codon_start=partmap.get(rank,{}).get('codon_start')
            for ri,row in enumerate(rows,1):
                if not isinstance(row,list) or len(row)<=pi:continue
                pv=num(row[pi])
                if pv is None or not 0<=pv<=1:continue
                global_codon=(int(codon_start)+ri-1) if codon_start not in (None,'') else ''
                out.append({'gene':gene,'method':'MEME','partition_rank':rank,'json_partition':part,
                            'site_index_within_partition':ri,'global_alignment_codon':global_codon,
                            'statistic':'p_value','value':pv,'selected_threshold':'p<=0.10',
                            'selected':str(pv<=0.10).upper(),'json_table_path':path})
    seen=set();z=[]
    for r in out:
        k=(r['partition_rank'],r['site_index_within_partition'],r['value'])
        if k not in seen:seen.add(k);z.append(r)
    return z


def parse_fubar(js,gene,partmap):
    out=[]
    for path,hdr,content in tables(js):
        low=[h.lower().replace(' ','') for h in hdr]
        pos=[];neg=[]
        for i,h in enumerate(low):
            if 'prob[alpha<beta]' in h or ('posterior' in h and ('positive' in h or 'diversifying' in h)):pos.append(i)
            if 'prob[alpha>beta]' in h or ('posterior' in h and ('negative' in h or 'purifying' in h)):neg.append(i)
        if not pos and not neg:continue
        for rank,part in enumerate(partition_order(content),1):
            rows=content[part]
            if not isinstance(rows,list):continue
            codon_start=partmap.get(rank,{}).get('codon_start')
            for ri,row in enumerate(rows,1):
                if not isinstance(row,list):continue
                global_codon=(int(codon_start)+ri-1) if codon_start not in (None,'') else ''
                for kind,idxs in [('positive',pos),('negative',neg)]:
                    if not idxs:continue
                    v=num(row[idxs[0]]) if len(row)>idxs[0] else None
                    if v is None or not 0<=v<=1:continue
                    out.append({'gene':gene,'method':'FUBAR','partition_rank':rank,'json_partition':part,
                                'site_index_within_partition':ri,'global_alignment_codon':global_codon,
                                'statistic':f'posterior_{kind}','value':v,'selected_threshold':'posterior>=0.90',
                                'selected':str(v>=0.90).upper(),'json_table_path':path})
    seen=set();z=[]
    for r in out:
        k=(r['partition_rank'],r['site_index_within_partition'],r['statistic'],r['value'])
        if k not in seen:seen.add(k);z.append(r)
    return z


def main():
    root=root_from_env_or_file();base=root/'06_TSWV_evolution'
    status_path=base/'tables'/'16_TSWV_GARD_final_status.tsv'
    audit_path=base/'tables'/'17_TSWV_selection_partition_audit.tsv'
    if not status_path.exists() or not audit_path.exists():
        raise SystemExit('Run 07O and 07P before 07J')
    status={r['gene']:r for r in read_tsv(status_path)}
    aud=read_tsv(audit_path)
    partmaps={g:{} for g in GENES}
    for r in aud:
        partmaps[r['gene']][int(r['partition_index'])]=r

    summary=[];sites=[];bad=[]
    for gene in GENES:
        sr=status[gene]
        row={'gene':gene,'GARD_status':sr['status'],'GARD_breakpoint_n':sr['final_breakpoint_n'],
             'GARD_baseline_cAIC':sr['baseline_cAIC'],'GARD_bestModelAICc':sr['bestModelAICc'],
             'selection_partition_n':len(partmaps[gene])}
        for method in ['BUSTED','MEME','FUBAR']:
            p=base/'pipeline_results'/'05_selection'/gene/f'{gene}.{method}.json'
            row[f'{method}_status']='missing'
            if not p.exists() or p.stat().st_size==0:
                bad.append(f'{gene}:{method}:missing');continue
            try:
                js=json.load(p.open())
            except Exception as e:
                bad.append(f'{gene}:{method}:invalid_json');continue
            row[f'{method}_status']='present_valid'
            if method=='BUSTED':
                path,pv=overall_p(js);row['BUSTED_pvalue']=pv if pv is not None else ''
                row['BUSTED_pvalue_path']=path;row['BUSTED_significant_p_lt_0.05']=str(pv is not None and pv<0.05).upper()
            elif method=='MEME':
                z=parse_meme(js,gene,partmaps[gene]);sites.extend(z)
                row['MEME_sites_parsed_n']=len({(x['partition_rank'],x['site_index_within_partition']) for x in z})
                row['MEME_sites_p_le_0.10']=sum(x['selected']=='TRUE' for x in z)
            else:
                z=parse_fubar(js,gene,partmaps[gene]);sites.extend(z)
                row['FUBAR_sites_parsed_n']=len({(x['partition_rank'],x['site_index_within_partition']) for x in z})
                row['FUBAR_positive_sites_posterior_ge_0.90']=sum(x['selected']=='TRUE' and x['statistic']=='posterior_positive' for x in z)
                row['FUBAR_negative_sites_posterior_ge_0.90']=sum(x['selected']=='TRUE' and x['statistic']=='posterior_negative' for x in z)
        summary.append(row)

    fields=[]
    for r in summary:
        for k in r:
            if k not in fields:fields.append(k)
    write_tsv(base/'tables'/'18_TSWV_GARD_selection_summary.tsv',summary,fields)
    site_fields=['gene','method','partition_rank','json_partition','site_index_within_partition','global_alignment_codon','statistic','value','selected_threshold','selected','json_table_path']
    write_tsv(base/'tables'/'19_TSWV_selection_sites_long.tsv',sites,site_fields)
    qc=[{'metric':'status','value':'PASS' if not bad else 'FAIL'},
        {'metric':'gene_n','value':len(GENES)},
        {'metric':'selection_json_expected_n','value':len(GENES)*3},
        {'metric':'selection_json_valid_n','value':len(GENES)*3-len(bad)},
        {'metric':'MEME_selected_row_n','value':sum(x['method']=='MEME' and x['selected']=='TRUE' for x in sites)},
        {'metric':'FUBAR_positive_selected_row_n','value':sum(x['method']=='FUBAR' and x['statistic']=='posterior_positive' and x['selected']=='TRUE' for x in sites)},
        {'metric':'FUBAR_negative_selected_row_n','value':sum(x['method']=='FUBAR' and x['statistic']=='posterior_negative' and x['selected']=='TRUE' for x in sites)}]
    write_tsv(base/'qc'/'09_TSWV_selection_QC.tsv',qc,['metric','value'])
    print('Fig7 07J selection summary', 'PASS' if not bad else 'FAIL')
    for r in summary:print(r)
    if bad:
        print('Missing/invalid:',','.join(bad),file=sys.stderr);return 2
    print('Raw HyPhy JSON remains source of truth. global_alignment_codon is alignment-coordinate mapping, not yet protein functional-domain annotation.')
    return 0

if __name__=='__main__':raise SystemExit(main())
