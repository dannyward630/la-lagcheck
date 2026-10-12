"""Reproduce the Figure 3 comparison from public numerical inputs.

The mapping is the LA-LagCheck default diagnosis of R01 stored in
reference_mapping.json. Run from any working directory. The output directory
must be new. The numerical audit is compared with expected_figure3.json.
"""
from pathlib import Path
import argparse
import re
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
import reconstruction as r
import feature_orientation as fo
local = fo.audit
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
OUT = args.output.resolve()
OUT.mkdir(parents=True, exist_ok=False)
HERE = Path(__file__).resolve().parent
ns = vars(r)
for k in ('qc','ql','ext','lo','hi','WIN','el_h','disp_um','blk_h','retained','tau_blk','tau_c','tmid_blk','model','V'):
    globals()[k] = ns[k]
main = local['comparisons'][1]
shift = local['constant_translation_mm']
window = local['display_window_mm']
data = HERE.parents[1] / 'study_summaries/Supplementary_Data_S1.csv'
M = pd.read_csv(data)
assert len(M)==29 and set(M.acquisition)=={f'R{i:02d}' for i in range(1,30)}
assert np.allclose(M.first_to_last_displacement_mm,np.abs(M.rate_ppm)*1e-6*M.duration_h*3600*M.speed_um_s/1000)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':8,'axes.labelsize':7.7,
                     'xtick.labelsize':7.5,'ytick.labelsize':7.5,'pdf.fonttype':42,'svg.fonttype':'none',
                     'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6})
cm=plt.get_cmap('viridis').copy();cm.set_bad('#cccccc')
blue='#0072B2';gray='.35'
SCALE=232.8/245  # Printed size 161.5 mm by 232.8 mm; all text is scaled by the same factor.
fig=plt.figure(figsize=(170*SCALE/25.4,245*SCALE/25.4))
gs=fig.add_gridspec(6,2,height_ratios=[1.05,.54,.54,.88,.065,.65],left=.075,right=.985,top=.955,bottom=.055,hspace=.65,wspace=.22)
whole_hours=float((ns['le'][-1]-ns['T0'])/3600)
for j,img in enumerate((qc,ql)):
    ax=fig.add_subplot(gs[0,j])
    im=ax.imshow(img,extent=ext,vmin=lo,vmax=hi,cmap=cm,aspect='equal',interpolation='nearest',rasterized=True)
    for key in ('B','C'):
        w=WIN[key]
        ax.add_patch(Rectangle((w['x'][0],w['y'][0]),np.ptp(w['x']),np.ptp(w['y']),fill=False,ec='white',lw=.9))
        ax.text(w['x'][0]+.15,w['y'][1]-.13,f'({key})',color='white',fontsize=7.5)
    xl=window[0]-(shift if j==0 else 0)
    ax.add_patch(Rectangle((xl,window[2]),window[1]-window[0],window[3]-window[2],fill=False,ec='white',lw=.9))
    ax.text(xl-.25,window[2]+.2,'(D)',color='white',ha='right',va='top',fontsize=7.5)
    ax.set_title(f'(A) Whole R01 scan: {whole_hours:.2f} h\nConstant offset' if j==0 else f'Whole R01 scan: {whole_hours:.2f} h\nAfter time mapping',loc='left')
    ax.set(xlabel='x (mm)',ylabel='y (mm)')
for grid_row,key in ((1,'B'),(2,'C')):
    w=WIN[key]
    assert len(w['blocks'])==2
    for j,img in enumerate((qc,ql)):
        ax=fig.add_subplot(gs[grid_row,j])
        ax.imshow(img,extent=ext,vmin=lo,vmax=hi,cmap=cm,aspect='equal',interpolation='nearest',rasterized=True)
        ax.set(xlim=w['x'],ylim=w['y'][::-1],xlabel='x (mm)',ylabel='y (mm)')
        for gx in np.linspace(w['x'][0],w['x'][1],5)[1:-1]:
            ax.axvline(gx,color='white',lw=.6,ls=':')
        span=f"{w['h'][0]:.1f} to {w['h'][1]:.1f} h"
        ax.set_title(f"({key}) {w['name']} displacement: {span}\nConstant offset" if j==0 else f"{span}\nAfter time mapping",loc='left')
        if j==0:
            xm=np.mean(w['x']);ya=w['y'][0]+.57*np.ptp(w['y']);delta=w['shift_mm']
            ax.annotate('',xy=(xm+delta,ya),xytext=(xm,ya),arrowprops=dict(arrowstyle='-|>',color='white',lw=1.4))
            ax.text(xm+delta/2,ya-.20,f'Δx {delta:+.2f} mm',color='white',ha='center',fontsize=7.5,
                    bbox=dict(fc='black',alpha=.55,ec='none',pad=1.2))
        rms=w['rms_const_mm'] if j==0 else w['rms_corr_mm']
        ax.text(.5,.03,f"Block-lag fit RMS (n = {len(w['blocks'])}): {rms:.4f} mm",transform=ax.transAxes,color='white',ha='center',fontsize=7.5,
                bbox=dict(fc='black',alpha=.55,ec='none',pad=1))
for j,img in enumerate((qc,ql)):
    ax=fig.add_subplot(gs[3,j]);e=list(ext)
    if j==0:e[:2]=[x+shift for x in e[:2]]
    ax.imshow(img,extent=e,vmin=lo,vmax=hi,cmap=cm,aspect='equal',interpolation='nearest',rasterized=True)
    p=main['constant'] if j==0 else main['mapped']
    center=np.array(p['center']);vec=np.array(p['axis']);line=np.array([center-.62*vec,center+.62*vec])
    effects=[pe.Stroke(linewidth=2.3,foreground='black'),pe.Normal()]
    ax.plot(line[:,0],line[:,1],color='white',lw=1.2,path_effects=effects)
    ax.plot([center[0]]*2,[center[1]-.65,center[1]+.65],color='white',lw=.8,ls=':',path_effects=effects)
    ax.set(xlim=window[:2],ylim=window[2:][::-1],xlabel='x (mm)',ylabel='y (mm)')
    ax.set_title('(D) Local feature: 3.0 to 3.8 h\nConstant offset' if j==0 else '3.0 to 3.8 h\nAfter time mapping',loc='left')
    ax.text(.04,.96,f"{p['angle_deg']:.1f}°",transform=ax.transAxes,color='white',fontsize=9,va='top',
            bbox=dict(fc='black',alpha=.6,ec='none',pad=1.5))
    if j==0: ax.text(.5,.025,f'Translated {shift:+.3f} mm',transform=ax.transAxes,color='white',fontsize=7.5,ha='center',
                    bbox=dict(fc='black',alpha=.6,ec='none',pad=.8))
cax=fig.add_subplot(gs[4,:]);pos=cax.get_position();cax.set_position([pos.x0+.20*pos.width,pos.y0,.60*pos.width,pos.height])
cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.set_label(r'Common colour scale: log$_{10}$ Fe56 (counts s$^{-1}$)',fontsize=7.5)
cb.ax.tick_params(labelsize=7.5)
ax=fig.add_subplot(gs[5,0])
for interval,label in [(WIN['B']['h'],'(B)'),(WIN['C']['h'],'(C)'),((3.0,3.8),'(D)')]:
    ax.axvspan(*interval,color='.5',alpha=.12,lw=0)
    ax.text(np.mean(interval),-1.20,label,ha='center',fontsize=7.5)
ax.axhline(0,color='.6',lw=.6)
ax.plot(el_h,disp_um/1000,color='black',lw=1)
ax.plot(blk_h[retained],(tau_blk[retained]-tau_c)*375/1000,'o',ms=3.5,color=blue)
ax.plot(blk_h[~retained],(tau_blk[~retained]-tau_c)*375/1000,'o',ms=3.5,mfc='white',mec=blue)
tz=float(el_h[np.argmin(abs(disp_um))]);ax.axvline(tz,color=gray,ls=':',lw=.7)
ax.text(tz+.1,-.67,f'Zero at {tz:.2f} h',fontsize=7.5,color=gray)
ax.text(.03,.9,f'Start {disp_um[0]/1000:+.3f} mm',transform=ax.transAxes,color=blue,fontsize=7.5)
ax.text(el_h[-1]-.1,1.13,f'End {disp_um[-1]/1000:+.3f} mm',color=blue,ha='right',fontsize=7.5)
ax.set(xlabel='Elapsed time (h)',ylabel='Δx (mm)',ylim=(-1.3,1.35))
ax.set_title('(E) Position difference over time',loc='left');ax.grid(alpha=.2,lw=.5)
ax=fig.add_subplot(gs[5,1]);markers={225:'o',281:'^',375:'s'}
for _,row in M.iterrows():
    highlight=row.acquisition=='R01'
    ax.plot(row.duration_h,row.first_to_last_displacement_mm,marker=markers[int(row.speed_um_s)],ls='none',
            color='black' if highlight else gray,ms=4.2 if highlight else 3.3,zorder=3)
for rid,position in [('R01',(12,.4)),('R29',(14.5,5.55))]:
    row=M[M.acquisition==rid].iloc[0]
    ax.annotate(rid,(row.duration_h,row.first_to_last_displacement_mm),xytext=position,fontsize=7.5,
                arrowprops=dict(arrowstyle='-',color=gray,lw=.5))
hh=np.linspace(0,22.6,100)
for speed,style in ((225,'-'),(281,'-.'),(375,'--')):
    ax.plot(hh,170e-6*hh*3600*speed/1000,style,color='.7',lw=.7,zorder=0)
ax.legend(handles=[Line2D([],[],marker=markers[s],ls=sty,color=gray,ms=3,label=f'{s} µm s$^{{-1}}$')
                   for s,sty in ((225,'-'),(281,'-.'),(375,'--'))],frameon=False,loc='upper left',fontsize=7.5)
ax.set(xlabel='Acquisition duration (h)',ylabel='First-to-last displacement\nfrom the fitted lag rate (mm)',xlim=(0,23),ylim=(0,6.2))
ax.set_title('(F) 29 routine acquisitions',loc='left');ax.grid(alpha=.2,lw=.5)
assert 29==len(M)
for panel in fig.axes:
    for artist in list(panel.texts) + [panel.title, panel._left_title, panel._right_title]:
        artist.set_text(re.sub(r'\(([A-F])\)', lambda m: '(' + m.group(1).lower() + ')', artist.get_text()))
for panel in fig.axes:
    panel.xaxis.labelpad=1.5
from matplotlib.text import Text
for artist in fig.findobj(Text):
    artist.set_fontsize(artist.get_fontsize()*SCALE)
fig.canvas.draw()
figbox=fig.bbox
outside=[]
check_artists=[]
for panel in fig.axes:
    check_artists += list(panel.texts) + [panel.title,panel._left_title,panel._right_title,panel.xaxis.label,panel.yaxis.label]
    if panel.get_legend() is not None: check_artists += list(panel.get_legend().get_texts())
for artist in check_artists:
    if not artist.get_visible() or not artist.get_text(): continue
    bb=artist.get_window_extent(fig.canvas.get_renderer())
    if bb.width and bb.height and (bb.x0 < -1 or bb.y0 < -1 or bb.x1 > figbox.x1+1 or bb.y1 > figbox.y1+1):
        outside.append(artist.get_text())
assert not outside, outside
for suffix in ('png','pdf','svg'):
    fig.savefig(OUT/f'Figure_3_R01_reproduced.{suffix}',dpi=600)
def feature_hours(img,xx):
    # Acquisition time of the scan lines that cross the extracted feature (threshold 5.0).
    _,mask=fo.measure(img,xx,5.0)
    rows=np.flatnonzero(mask.any(axis=1))
    lines=np.flatnonzero(np.isin(np.round(ns['ys'],1),ns['uy'][rows]))
    return [float(el_h[lines].min()),float(el_h[lines].max())]
audit=dict(rate_iter_ppm=ns['R_IT']*1e6,starting_alignment_s=ns['ALIGN'],iterated_fit_intercept_s=ns['RES'],
           constant_offset_s=tau_c,colour_scale_log10=[float(lo),float(hi)],array_shape=list(qc.shape),
           local_orientation=local,early_window=WIN['B'],
           late_window=WIN['C'],whole_duration_h=whole_hours,panel_order=['a','b','c','d','e','f'],n_acquisitions=len(M),displacement_source='study_summaries/Supplementary_Data_S1.csv',
           figure_width_mm=170*SCALE,figure_height_mm=245*SCALE,minimum_font_pt=7.5*SCALE,
           feature_h_constant=feature_hours(qc,fo.x+shift),feature_h_mapped=feature_hours(ql,fo.x),
           early_shift_spot_widths=abs(WIN['B']['shift_mm'])*1000/15,late_shift_spot_widths=abs(WIN['C']['shift_mm'])*1000/15,
           blocks=dict(block=[int(b) for b in ns['blk'].block],elapsed_h=[float(v) for v in blk_h],
                       lag_minus_constant_offset_mm=[float(v) for v in (tau_blk-tau_c)*V/1000],
                       retained_after_iteration=[bool(v) for v in retained]),
           n_blocks_drawn_filled=int(retained.sum()),n_blocks_drawn_open=int((~retained).sum()),
           mapping_zero_h=tz,first_mm=float(disp_um[0]/1000),last_mm=float(disp_um[-1]/1000),
           r01_summary_single_pass_rate_ppm=float(M[M.acquisition=='R01'].rate_ppm.iloc[0]),
           r01_summary_first_to_last_displacement_mm=float(M[M.acquisition=='R01'].first_to_last_displacement_mm.iloc[0]),
           input_representation='Anonymised saved-workspace signals and intervals; not a native vendor-export validation.')
(OUT/'Figure_3_R01_reproduced_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print('Saved Figure 3 reproduction from anonymised R01 inputs.')
# Compare every number of the audit with the stored expectation.
expected=json.loads((HERE/'expected_figure3.json').read_text(encoding='utf-8'))
def compare(a,b,path=''):
    if isinstance(a,dict):
        assert isinstance(b,dict) and set(a)==set(b),(path,sorted(set(a)^set(b)))
        return sum(compare(a[k],b[k],path+'/'+k) for k in a)
    if isinstance(a,(list,tuple)):
        assert len(a)==len(b),path
        return sum(compare(x,y,path+f'[{i}]') for i,(x,y) in enumerate(zip(a,b)))
    if isinstance(a,float) or isinstance(b,float):
        assert abs(a-b)<=1e-9*max(1.,abs(b)),(path,a,b)
        return 1
    assert a==b,(path,a,b)
    return 1
n_values=compare(json.loads(json.dumps(audit)),expected)
print(f'Audit matches expected_figure3.json: {n_values} values compared, relative tolerance 1e-9.')
