# 把原始 .Spectrum 檔整理成學生用的 tidy CSV（不更動原始檔）
import numpy as np, csv, os
root=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..')
def load(f):
    d=np.loadtxt(os.path.join(root,f),skiprows=1); o=np.argsort(d[:,0]); return d[o,0],d[o,1]
waterf={0.25:'water_0.25mm_Abs_1.Spectrum',0.5:'water_0.5mm_Abs_1.Spectrum',1:'water_1.0_abs__Abs_1.Spectrum',
        2:'water_2mm_1_Abs_1.Spectrum',5:'water_5mm_Abs_1_Abs_1.Spectrum'}
acid={'3%':3.32,'5%':4.87,'10%':7.45}
rows=[]; wide={}
for p,wf in waterf.items():
    ps='%g'%p
    items=[('water',0.0,0.0,wf,1)]+[('vinegar_'+k,float(k[:-1]),v,f'Acetic_Acid_{ps}mm_{k}_1_Abs_1.Spectrum',1) for k,v in acid.items()]
    if p==0.25: items.append(('water',0.0,0.0,'water_025mm_Abs_1.Spectrum',2))
    for s,lab,atago,f,rep in items:
        wl,A=load(f)
        for w,a in zip(wl,A): rows.append([p,s,lab,atago,rep,round(w,3),a])
        if rep==1: wide[f'{ps}mm_{atago:g}']=A
        WL=wl
out=os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(out,'spectra_long.csv'),'w',newline='',encoding='utf-8') as fh:
    w=csv.writer(fh); w.writerow(['path_mm','sample','label_pct','acid_pct','replicate','wavelength_nm','absorbance'])
    for r in rows: w.writerow(r[:6]+['%.6f'%r[6]])
print(len(rows),'rows;',len(WL),'wavelengths', WL[0],WL[-1])
