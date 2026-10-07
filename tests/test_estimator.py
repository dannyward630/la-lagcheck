import unittest
import numpy as np
from lagcheck.models import SignalData,LaserData,InputError
from lagcheck.estimator import diagnose


def synthetic(rate=170., offset=3., step=0.):
    rng=np.random.default_rng(148)
    lines=[]
    cursor=5.
    for i in range(240):
        length=float(rng.uniform(5,18))
        lines.append({'Start Time':cursor,'End Time':cursor+length})
        cursor+=length+float(rng.uniform(3,10))
    t=np.arange(0,cursor+30,.1)
    signal=np.ones(len(t))*2.
    origin=lines[0]['Start Time']
    for i,row in enumerate(lines):
        start=origin+(1+rate*1e-6)*(row['Start Time']-origin)+offset
        end=origin+(1+rate*1e-6)*(row['End Time']-origin)+offset
        if i>=120:
            start+=step;end+=step
        # Fractional exposure over each 100 ms sample provides a known physical timing model.
        fraction=np.maximum(0,np.minimum(t+.05,end)-np.maximum(t-.05,start))/.1
        signal+=fraction*1000
    signal+=rng.normal(0,.05,len(t))
    return SignalData(t,signal),LaserData(lines)


class EstimatorTests(unittest.TestCase):
    def test_positive_and_negative_known_rates_and_honest_qc(self):
        for rate in (170.,-170.,0.):
            with self.subTest(rate=rate):
                ms,laser=synthetic(rate)
                r,_=diagnose(ms,laser,{'tolerance_s':.2,'alignment_offset_s':3.})
                self.assertAlmostEqual(r['rate_ppm'],rate,delta=3.)
                self.assertIn(r['status'],('supported','qc_failed'))
                if not r['closure']['converged']:
                    self.assertEqual(r['status'],'qc_failed')
                    self.assertIn('closure_converged',r['failed_checks'])

    def test_offset_is_not_rate(self):
        results=[]
        for offset in (2.,8.):
            ms,laser=synthetic(100.,offset)
            r,_=diagnose(ms,laser,{'tolerance_s':.2,'alignment_offset_s':offset})
            self.assertEqual(r['status'],'supported',r)
            results.append(r['rate_ppm'])
        self.assertAlmostEqual(*results,delta=.1)

    def test_step_is_not_accepted_as_a_single_rate(self):
        ms,laser=synthetic(step=2.)
        r,_=diagnose(ms,laser,{'tolerance_s':.2,'alignment_offset_s':3.})
        self.assertNotEqual(r['status'],'supported')

    def test_temporary_middle_shift_cannot_be_hidden_by_outlier_rejection(self):
        ms,laser=synthetic(100.,3.)
        ms.signal=np.full(len(ms.time),2.)
        for i,row in enumerate(laser.lines):
            extra=8. if 96<=i<144 else 0.
            start=5+(row['Start Time']-5)*1.0001+3+extra
            end=5+(row['End Time']-5)*1.0001+3+extra
            fraction=np.maximum(0,np.minimum(ms.time+.05,end)-np.maximum(ms.time-.05,start))/.1
            ms.signal+=fraction*1000
        ms.signal+=np.random.default_rng(149).normal(0,.05,len(ms.time))
        r,blocks=diagnose(ms,laser,{'tolerance_s':.2,'alignment_offset_s':3.})
        self.assertEqual(r['status'],'qc_failed',r)
        self.assertIn('no_contradictory_high_quality_blocks',r['failed_checks'])
        self.assertEqual(r['contradictory_high_quality_blocks'],[4,5])
        self.assertTrue(all(b['gate_pass'] for b in blocks if b['block'] in (4,5)))

    def test_alternative_blocking_keeps_contradictory_evidence(self):
        # Actual waveform, not a mocked fit: the 16 shifted lines straddle two
        # default blocks, but fill block 7 exactly when 15 blocks are used.
        # Before the fix this returned supported, 100.779581 ppm, although
        # the alternative block had correlation .995 and an 8 s residual.
        ms,laser=synthetic(100.,3.)
        ms.signal=np.full(len(ms.time),2.)
        for i,row in enumerate(laser.lines):
            extra=8. if 112<=i<128 else 0.
            start=5+(row['Start Time']-5)*1.0001+3+extra
            end=5+(row['End Time']-5)*1.0001+3+extra
            fraction=np.maximum(0,np.minimum(ms.time+.05,end)-np.maximum(ms.time-.05,start))/.1
            ms.signal+=fraction*1000
        ms.signal+=np.random.default_rng(149).normal(0,.05,len(ms.time))
        r,_=diagnose(ms,laser,{'tolerance_s':.2,'alignment_offset_s':3.})
        self.assertEqual(r['contradictory_high_quality_blocks'],[])
        self.assertEqual(r['status'],'qc_failed')
        self.assertIn('all_sensitivity_evidence_supported',r['failed_checks'])
        variant=next(v for v in r['sensitivity'] if v['config']=='n_blocks=15')
        self.assertFalse(variant['evidence_pass'])
        self.assertIn(7,variant['contradictory_high_quality_blocks'])
        block=next(b for b in variant['blocks'] if b['block']==7)
        self.assertGreater(block['peak'],.99)
        self.assertGreater(block['margin'],.2)

    def test_clipped_last_block_cannot_claim_unobserved_coverage(self):
        # A uniform, perfectly shifted waveform remains highly correlated in
        # the observed part of the final block. Its last logged interval ends
        # 3000 s after the acquired signal, so nominal block endpoints cannot
        # establish coverage of the full log.
        rng=np.random.default_rng(148)
        lines=[]
        cursor=5.
        for _ in range(240):
            length=round(float(rng.uniform(5,18)),1)
            lines.append({'Start Time':cursor,'End Time':round(cursor+length,1)})
            cursor=round(cursor+length+float(rng.uniform(3,10)),1)
        lines[-1]['End Time']+=3000.
        t=np.arange(0,cursor+30,.1)
        signal=np.full(len(t),2.)
        for row in lines:
            signal[(t>=round(row['Start Time']-3,1)) &
                   (t<round(row['End Time']-3,1))]=1000.
        maximum_visible=(t[-1]-t[0])/(lines[-1]['End Time']-lines[0]['Start Time'])
        self.assertLess(maximum_visible,.6)
        r,blocks=diagnose(SignalData(t,signal),LaserData(lines),
                          {'tolerance_s':.2,'alignment_offset_s':0.})
        self.assertIn(r['status'],('qc_failed','insufficient_evidence'))
        if 'coverage_fraction' in r:
            self.assertLessEqual(r['coverage_fraction'],maximum_visible+1e-9)
        if blocks:
            last=next(b for b in blocks if b['block']==9)
            self.assertGreater(last['peak'],.99)
            self.assertFalse(last['fit_eligible'])
            self.assertLess(last['observed_n_lines'],last['n_lines'])

    def test_partial_tail_block_still_contributes_contradictory_evidence(self):
        # The 23 completely visible tail intervals clearly shift by -8 s.
        # Only the last interval extends outside the acquired signal. Excluding
        # this partial block from fitting must not suppress its visible evidence.
        # The unmodified pre-R3 estimator returned supported at 87.7% coverage.
        rng=np.random.default_rng(148)
        lines=[]
        cursor=5.
        for _ in range(240):
            length=round(float(rng.uniform(5,18)),1)
            lines.append({'Start Time':cursor,'End Time':round(cursor+length,1)})
            cursor=round(cursor+length+float(rng.uniform(3,10)),1)
        lines[-1]['End Time']+=50.
        t=np.arange(0,cursor+30,.1)
        signal=np.full(len(t),2.)
        for i,row in enumerate(lines):
            lag=-11. if i>=216 else -3.
            signal[(t>=round(row['Start Time']+lag,1)) &
                   (t<round(row['End Time']+lag,1))]=1000.
        r,blocks=diagnose(SignalData(t,signal),LaserData(lines),
                          {'tolerance_s':.2,'alignment_offset_s':0.})
        self.assertGreater(r['coverage_fraction'],.8)
        self.assertTrue(r['checks']['coverage'])
        last=next(b for b in blocks if b['block']==9)
        self.assertFalse(last['fit_eligible'])
        self.assertEqual(last['observed_n_lines'],23)
        self.assertGreater(last['peak'],.99)
        self.assertGreater(last['margin'],.2)
        self.assertAlmostEqual(last['lag_s'],-11.,delta=.01)
        self.assertEqual(r['status'],'qc_failed')
        self.assertIn(9,r['contradictory_high_quality_blocks'])
        self.assertIn('no_contradictory_high_quality_blocks',r['failed_checks'])

    def test_flat_signal_has_no_estimate(self):
        ms,laser=synthetic()
        ms.signal[:]=1
        r,_=diagnose(ms,laser,{'tolerance_s':.2})
        self.assertEqual(r['status'],'insufficient_evidence')
        self.assertIsNone(r['rate_ppm'])

    def test_too_few_lines(self):
        ms,laser=synthetic()
        laser.lines=laser.lines[:30]
        r,_=diagnose(ms,laser,{'tolerance_s':.2})
        self.assertEqual(r['reason'],'insufficient_laser_intervals')

    def test_gap_rejected(self):
        ms,laser=synthetic()
        ms.time[len(ms.time)//2:]+=3
        r,_=diagnose(ms,laser,{'tolerance_s':.2})
        self.assertEqual(r['reason'],'time_grid_not_supported_by_index_shift_estimator')

    def test_smooth_grid_warp_rejected(self):
        ms,laser=synthetic()
        ms.time+=.4*np.sin(np.linspace(0,2*np.pi,len(ms.time)))
        r,_=diagnose(ms,laser,{'tolerance_s':.2})
        self.assertEqual(r['reason'],'index_grid_approximation_exceeds_timing_budget')

    def test_explicit_tolerance_required(self):
        ms,laser=synthetic()
        with self.assertRaises(InputError): diagnose(ms,laser,{})

    def test_automatic_constant_alignment(self):
        ms,laser=synthetic(170.,3.)
        r,_=diagnose(ms,laser,{'tolerance_s':.2,'coarse_window_s':20.})
        self.assertTrue(r['alignment']['ok'],r)
        self.assertAlmostEqual(r['rate_ppm'],170.,delta=3.)

    def test_fixed_physical_lag_on_warped_grid_is_rejected(self):
        rng=np.random.default_rng(20260919)
        t=np.r_[0,np.cumsum(.05*(1+np.linspace(-.09,.09,24000)))]
        lines=[]; cursor=20.
        y=np.ones(len(t))
        for _ in range(200):
            duration=rng.uniform(2.5,3.8)
            lines.append({'Start Time':cursor,'End Time':cursor+duration})
            y[(t>=cursor+1)&(t<cursor+duration+1)]=1000
            cursor+=duration+rng.uniform(.6,2.5)
        r,_=diagnose(SignalData(t,y),LaserData(lines),{'tolerance_s':.15,'alignment_offset_s':0.})
        self.assertEqual(r['status'],'insufficient_evidence')
        self.assertIsNone(r['rate_ppm'])

if __name__=='__main__':unittest.main()
