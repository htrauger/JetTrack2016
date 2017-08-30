Implementation of this code for jet-track analysis:

(1) make_ntuples_run2 contains code used to skim the HiForest to produce mini-ntuples with only the quantities relevant for this analysis. JFF-dependent JEC are also applied to jtpt at this analysis step, and written into a branch called "corrpt".  Only used for pre-approval studies. 

(2) JetTrack_Inclusive/JetTrack_Inclusive3.C is a unified analyzer that produces correlations from PbPb and pp data as well as pythia and pythia+hydjet simulation.

The output from this analyzer is then (locally) taken through the analysis procedure described in AN-16-298:

(3) Mixed event correction (me_correct/me_correct3.cxx) applies mixed event correction, residual ME shape correction (if flag is set), AND background subtraction -- run once per MC or data scenario to produce a single file containing pTweighted and unweighted correlations as well as raw signals (per-jet normalized), mixed event correlations (normalized to 1), and 2D backgrounds

(4) Determination of the background fluctuation bias correction (spill_over/spill_over.cxx), runs on SubeNon0 correlations -- flag controls pTweighted or not (must be run separately for each)

(5) Determination of the residual JFF-JEC/swapping correction (jff_residual/jff_residual.cxx), runs on Sube0 correlations --  flag controls pTweighted or not (must be run separately for each)

(6) Projection in dEta/dPhi, correction, and application of systematic uncertainties (particle_yields/particle_yields.cxx) --  flag controls pTweighted or not (must be run separately for each)

(7) Yields by dR and jet shapes (jet_shapes_result/jet_shapes_result.cxx) --  flag controls pTweighted or not (must be run separately for each)

(8) Final plotting of all PAS plots (results_plotting/results_plotting.cxx) -- four scenarios (pTweighted or not, with QM reference or not), plot_all.sh runs all four

