#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLatex.h"
#include "THStack.h"

#include <iostream>
#include <vector>
#include <fstream>

#include "../JetTrack2016_functions.h"

//******************************************************
// VERSION TO DO JET-SHAPES FOR PAS
//******************************************************

Int_t jet_shapes_result(bool is_number=0, bool use_highpT_bin = kTRUE, bool do_ref=kFALSE){

  TFile *fMC[6];
  //  TFile *f_ref[7];

  // TFile *f_ref = new TFile("../FROZEN_PREAPPROVAL/jet_shapes_result/Jet_Shapes.root");
  TFile *f_ref = new TFile("../../HIN-15-011/dR_studies/Jet_Shapes.root");

  TCanvas *c_jetshape[5];
   
  TString jetetacut, etalabel,centlabel,pTlabel,ajlabel;
  float eta_ymax;

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins, limR;
  float deta, dphi, r, bc, bg_err,temp1,temp2, rbin, temperr, err, width_temp_x, width_temp_y, width_temp, norm_temp, zerobin, temp, norm_tot, err_temp, cont;
  
 
  float RBins[20] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.,1.2,1.4,1.6,1.8,2.0};
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 10;

  float PtBins[nPtBins+1] = {100, 300};
  
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt1000"};
  

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

   float TrkPtBins[nTrkPtBins] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 300};
   TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300","" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};

  float mean_pts[nTrkPtBins] = {0.844,1.35,2.35,3.37,5.07,9.72,13.8,17.9,22.};
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.25);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
 
 
  TH2D* result[6][nCBins][nPtBins][nTrkPtBins];

  TH2D* resultMC_gen[6][nCBins][nPtBins][nTrkPtBins];
  TH2D* resultMC_reco[6][nCBins][nPtBins][nTrkPtBins];

  TH2D* resultMC[6][nCBins][nPtBins][nTrkPtBins];

  TF1 *gaus_eta[12][nTrkPtBins][nCBins];
  TF1 *gaus_phi[12][nTrkPtBins][nCBins];

  TH1D *spill_over_dEta[12][nTrkPtBins][nCBins];
  TH1D *spill_over_dPhi[12][nTrkPtBins][nCBins];

 
  float par0, par1, par2, par3, par4;

  TH1D *background_syst_rebin[12][6][4][5];

   
  TH1D* JetShapeMC[6][nCBins][nTrkPtBins];
  TH1D* JetShapeMC_reco[6][nCBins][nTrkPtBins];

  TH1D* JetShape[6][nCBins][nTrkPtBins];
  TH1D* JetShape2[6][nCBins][nTrkPtBins];
  TH1D* JetShape_noerr_up[6][nCBins][nTrkPtBins];
  TH1D* JetShape_noerr_down[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ref[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ref_diff[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff[6][nCBins][nTrkPtBins];

  TH1D* JetShape_ref_ratio[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ratio[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ratio2[6][nCBins][nTrkPtBins];

  TH1D* JetShape_diff2[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff_noerr_up[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff_noerr_down[6][nCBins][nTrkPtBins];
  TH1D* JetShape2_geo[6][nCBins][nTrkPtBins];

  TH1D* JetShape_syst[6][nCBins][nTrkPtBins];
  TGraphAsymmErrors* JetShape_graph[6][nCBins][nTrkPtBins];



  THStack *JetShape_Stack_Up[6][nCBins];
  THStack *JetShape_Diff_Stack_Up[6][nCBins];

  THStack *JetShape_Stack_Down[6][nCBins];
  THStack *JetShape_Diff_Stack_Down[6][nCBins];

 
  double temp_cont, nextr, nextl, cos_weight,me00_range, mc_error,me_error, bg_error, me_error2, bg_error2;
  
  TString stem, datalabel,me00_range_string,stem_mc;

  TFile *f_deta_dphi = new TFile("../particle_yields/Particle_Yields.root","READ");

  if(!is_number)f_deta_dphi = new TFile("../particle_yields/Particle_Yields_pTweighted.root","READ");

  TH2D *bg_err_pp = (TH2D*)f_deta_dphi->Get("bg_err_pp")->Clone("bg_err_pp");
  TH2D *bg_err_pbpb = (TH2D*)f_deta_dphi->Get("bg_err_PbPb")->Clone("bg_err_PbPb");
  TH2D *me_err_pp = (TH2D*)f_deta_dphi->Get("me_err_pp")->Clone("me_err_pp");
  TH2D *me_err_pbpb = (TH2D*)f_deta_dphi->Get("me_err_PbPb")->Clone("me_err_PbPb");

 
  float norm;

  TFile *fin_pp  = new TFile("../me_correct/pp_Inclusive_Correlations.root", "READ");
  TFile *fin_PbPb = new TFile("../me_correct/PbPb_Inclusive_Correlations.root", "READ");
  TFile *fin_run1 = new TFile("../../JetTrack2015/bg_fit/Dijet_Correlations.root", "READ");
  TFile *fout;
  if(is_number) fout  = new TFile("Jet_Shapes.root","RECREATE"); 
  else 	 fout  = new TFile("Jet_Shapes_pTweighted.root","RECREATE"); 

  //fout  = new TFile("Jet_Shapes_NoCorrection.root","RECREATE"); 

  TFile *f_jff_pyth_gen = new TFile("../me_correct/Pythia_GenJet_GenTrack_Inclusive_Correlations.root", "READ");
  TFile *f_jff_pyth_reco = new TFile("../me_correct/Pythia_RecoJet_GenTrack_Inclusive_Correlations.root", "READ");
  TFile *f_jff_pyth_reco_reco = new TFile("../me_correct/Pythia_RecoJet_RecoTrack_Inclusive_Correlations.root", "READ");
  TFile *f_jff_hyd_gen = new TFile("../me_correct/HydJet_GenJet_GenTrack_Sube0_Inclusive_Correlations.root", "READ");
  TFile *f_jff_hyd_reco = new TFile("../me_correct/HydJet_RecoJet_GenTrack_Sube0_Inclusive_Correlations.root", "READ");
  TFile *f_jff_hyd_reco_reco = new TFile("../me_correct/HydJet_RecoJet_RecoTrack_Inclusive_Correlations.root", "READ");
  
  TFile *f_spillover = new TFile("../spill_over/Inclusive_Hydjet_SpillOvers.root");
  if(!is_number) f_spillover = new TFile("../spill_over/Inclusive_Hydjet_SpillOvers_pTweighted.root");

 
  TH1D *Integral_Pt[12][5];
  TH1D *Integral_syst_Pt[12][5];
  TH1D *Integral_diff_Pt[12][5];
  TH1D *Integral_diff_syst_Pt[12][5];


  double integral[12][nCBins][nTrkPtBins];
  double integral_err[12][nCBins][nTrkPtBins];
  double integral_syst_err[12][nCBins][nTrkPtBins];

  float rel_err_pbpb =TMath::Sqrt(0.05*0.05+0.01*0.01+0.04*0.04);
  float rel_err_pp =TMath::Sqrt(0.05*0.05+0.04*0.04+0.03*0.03);

  TPaveText *labels;

  TF2* ClosureFit = new TF2("ClosureFit", "[0]/2/TMath::Pi()/[1]/[2]*TMath::Exp(-1.*(x*x/[1]/[1]/2))*TMath::Exp(-1.*(y*y/[2]/[2]/2))",-5.,5.,-TMath::Pi()/2.,3*TMath::Pi()/2.);

  TCanvas *mc_canvas = new TCanvas("mc_canvas","",10,10,1500,4000);
  mc_canvas->Divide(4,9,0.,0.);

  //-----------------------
  // Start getting histos
  //-----------------------
  for(int g=0; g<2; g++){


    switch(g){
    case 0:
      if(is_number) stem = "Yield_BkgSub_";
      else  stem = "Yield_BkgSub_pTweighted";
      //   else  stem = "Yield_pTweighted_PbPb_";
      datalabel = "Inclusive";
      break;
    case 1:
      if(is_number) stem = "Yield_BkgSub_";
      //   else stem = "Yield_pTweighted_pp_";
      else stem = "Yield_BkgSub_pTweighted";
      //stem = "Yield_BkgSub_";
      datalabel = "Inclusive";
      break;
    }

    if(is_number) stem_mc = "Raw_Yield_";
    else  stem_mc = "Raw_Yield_pTweighted";

    for (int ibin=0;ibin<nCBins;ibin++){

      if(g==1&&ibin > 0) continue;

      for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
    
	for (int ibin3=0;ibin3<nTrkPtBins-1;ibin3++){

	  cout<<g<<" "<<ibin<<" "<<ibin3<<endl;

	  if(g==0){


	    result[g][ibin][ibin2][ibin3] = (TH2D*)fin_PbPb->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + datalabel+"_PbPb_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    cout<<"1"<<endl;
	    if(ibin3<4){
	 
	      resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_hyd_gen->Get((TString)(stem_mc+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Yield_BkgSub_GenGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      if(ibin3>6) resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_hyd_reco_reco->Get((TString)(stem_mc + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Yield_BkgSub_RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    	    
	      else resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_hyd_reco->Get((TString)(stem_mc + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Yield_BkgSub_RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	      TString	spillover_name =(TString)("Eta_SpillOver_Points_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]);
	      spill_over_dEta[g][ibin][ibin3] = (TH1D*)f_spillover->Get(spillover_name)->Clone(spillover_name);
	    
	 
	      spillover_name.ReplaceAll("Eta","Phi");
	  
	      spill_over_dPhi[g][ibin][ibin3] = (TH1D*)f_spillover->Get(spillover_name)->Clone(spillover_name);
	  

	      gaus_eta[g][ibin][ibin3] = (TF1*)f_spillover->Get((TString)("Eta_SpillOver_Fit_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	  
	      gaus_phi[g][ibin][ibin3] = (TF1*)f_spillover->Get((TString)("Phi_SpillOver_Fit_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	      cout<<"3"<<endl;
	   
	      par0= gaus_eta[g][ibin][ibin3]->GetParameter(1);
	      par1= gaus_eta[g][ibin][ibin3]->GetParameter(2);
	      par2= gaus_phi[g][ibin][ibin3]->GetParameter(2);

	      ClosureFit->SetParameter(0,par0);
	      ClosureFit->SetParameter(1,par1);
	      ClosureFit->SetParameter(2,par2);

	    
	      resultMC[g][ibin][0][ibin3] = new TH2D((TString)("Closure_2D_"+datalabel+"_"+CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]),"",500,-5,5,200,-TMath::Pi()/2,3*TMath::Pi()/2);

	      resultMC[g][ibin][0][ibin3]->Eval(ClosureFit);

	      width_temp_x = resultMC[g][ibin][0][ibin3]->GetXaxis()->GetBinWidth(1);
	      width_temp_y = resultMC[g][ibin][0][ibin3]->GetYaxis()->GetBinWidth(1);

	      resultMC[g][ibin][0][ibin3]->Scale(width_temp_x*width_temp_y);

	    
	      err = TMath::Sqrt(spill_over_dEta[g][ibin][ibin3]->GetBinError(1)*spill_over_dEta[g][ibin][ibin3]->GetBinError(1)+spill_over_dPhi[g][ibin][ibin3]->GetBinError(1)*spill_over_dPhi[g][ibin][ibin3]->GetBinError(1))*width_temp_x*width_temp_y/spill_over_dPhi[g][ibin][ibin3]->GetBinWidth(1)/spill_over_dEta[g][ibin][ibin3]->GetBinWidth(1);


	      

	      for(int k = 0; k< resultMC[g][ibin][0][ibin3]->GetNbinsX()+1; k++){
		for(int m = 0; m< resultMC[g][ibin][0][ibin3]->GetNbinsY()+1; m++){

		  resultMC[g][ibin][0][ibin3]->SetBinError(k,m,err);
		}
	      }
	    
	    }else{
 
	    resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_hyd_gen->Get((TString)(stem_mc+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Yield_BkgSub_GenGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    if(ibin3>6) resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_hyd_reco_reco->Get((TString)(stem_mc + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Yield_BkgSub_RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    	    
	    else resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_hyd_reco->Get((TString)(stem_mc + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Yield_BkgSub_RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));



	    }	  	    
	  }else{
	    
	    result[g][ibin][ibin2][ibin3] = (TH2D*)fin_pp->Get((TString)(stem +"Cent0_Cent10_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + datalabel+"_pp_"+ "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    
	    if(ibin3<4){
	      resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_pyth_gen->Get((TString)(stem_mc+ CBin_strs[0] + "_" + CBin_strs[1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"GenGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    
	      if(ibin3>6) resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_pyth_reco_reco->Get((TString)(stem_mc + CBin_strs[0] + "_" + CBin_strs[1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      else  resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_pyth_reco->Get((TString)(stem_mc + CBin_strs[0] + "_" + CBin_strs[1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    }else{
	      resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_pyth_gen->Get((TString)(stem_mc+ CBin_strs[0] + "_" + CBin_strs[1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"GenGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    
	      if(ibin3>6) resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_pyth_reco_reco->Get((TString)(stem_mc + CBin_strs[0] + "_" + CBin_strs[1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      else  resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_pyth_reco->Get((TString)(stem_mc + CBin_strs[0] + "_" + CBin_strs[1] + "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));



	    }
	  }
	  
	  resultMC_reco[g][ibin][0][ibin3]->Add(resultMC_gen[g][ibin][0][ibin3],-1.);

	  for(int k = 0; k< resultMC_reco[g][ibin][0][ibin3]->GetNbinsX()+1; k++){
	    for(int m = 0; m< resultMC_reco[g][ibin][0][ibin3]->GetNbinsY()+1; m++){


	      deta = resultMC_reco[g][ibin][0][ibin3]->GetXaxis()->GetBinCenter(k);
	      dphi = resultMC_reco[g][ibin][0][ibin3]->GetYaxis()->GetBinCenter(m);
	   
	      r = TMath::Sqrt(deta*deta+dphi*dphi);
	      
	      
	      if(r>0.5||(r>0.3&&ibin3>3)){
		resultMC_reco[g][ibin][0][ibin3]->SetBinContent(k,m,0.);
		resultMC_reco[g][ibin][0][ibin3]->SetBinError(k,m,0.);
	      }
	      
	    }
	   
	  }
	  
	
	
	  if(g==0&&ibin3<4){

	    if(ibin3==0&&is_number)resultMC[g][ibin][0][ibin3]->Scale(0.3); // The spill-overs are pre-scaled for # correlations, but nothing else is.

	 	    //	    if(!is_number)resultMC[g][ibin][0][ibin3]->Scale(mean_pts[ibin3]); //we use unweighted spillovers and weight on average
	    resultMC_reco[g][ibin][ibin2][ibin3]->Add(resultMC[g][ibin][0][ibin3]);
	  }
	 

	  result[g][ibin][ibin2][ibin3]->Add(resultMC_reco[g][ibin][0][ibin3],-1.); //applies correction



	  if(is_number){                                            //this is the dpt step
	    if(ibin3==0) result[g][ibin][ibin2][ibin3]->Scale(1./0.3);
	    if(ibin3> 3&& ibin3<nTrkPtBins-1) result[g][ibin][ibin2][ibin3]->Scale(1./4);

	    if(ibin3==0&&g==0)   resultMC[g][ibin][ibin2][ibin3]->Scale(1./0.3);
	    if(ibin3==0)   resultMC_reco[g][ibin][ibin2][ibin3]->Scale(1./0.3);

	    // if(ibin3> 3&& ibin3<nTrkPtBins-1) resultMC[g][ibin][ibin2][ibin3]->Scale(1./4);
	    if(ibin3> 3&& ibin3<nTrkPtBins-1) resultMC_reco[g][ibin][ibin2][ibin3]->Scale(1./4);
	    
	  }

	  if(ibin3==0){
	    result[g][ibin][ibin2][9] = (TH2D*) result[g][ibin][ibin2][ibin3]->Clone(((TString) ("Summed_result_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])));
	 
	    resultMC_reco[g][ibin][ibin2][9] = (TH2D*) resultMC_reco[g][ibin][ibin2][ibin3]->Clone(((TString) ("Summed_result_MC_reco_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+
1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])));


	    if(g==0)  resultMC[g][ibin][ibin2][9] = (TH2D*) resultMC[g][ibin][ibin2][ibin3]->Clone(((TString) ("Summed_result_MC_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+
1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])));

	  }else if(!is_number&&ibin3>0&&ibin3<nTrkPtBins-1){
	    
	    result[g][ibin][ibin2][9]->Add(result[g][ibin][ibin2][ibin3]);
	    resultMC_reco[g][ibin][ibin2][9]->Add(resultMC_reco[g][ibin][ibin2][ibin3]);
	    if(g==0&&ibin3<4) resultMC[g][ibin][ibin2][9]->Add(resultMC[g][ibin][ibin2][ibin3]);
	 
	  }else if(is_number&&ibin3>0&&ibin3<nTrkPtBins-2){
	    
	    result[g][ibin][ibin2][9]->Add(result[g][ibin][ibin2][ibin3]);
	    resultMC_reco[g][ibin][ibin2][9]->Add(resultMC_reco[g][ibin][ibin2][ibin3]);
	    if(g==0&&ibin3<4)resultMC[g][ibin][ibin2][9]->Add(resultMC[g][ibin][ibin2][ibin3]);
	 
	  }
	  
	}
    


	for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

	  if(g==0){	  
	  JetShape2[g][ibin][ibin3] = new TH1D((TString)("JetShape2_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  
 
	  JetShapeMC[g][ibin][ibin3] = new TH1D((TString)("JetShapeMC_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  

	  JetShapeMC_reco[g][ibin][ibin3] = new TH1D((TString)("JetShapeMC_reco_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  

	  JetShape2_geo[g][ibin][ibin3] = new TH1D((TString)("JetShape_Geometry_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  
	 
	  }else{
	    JetShape2[g][ibin][ibin3] = new TH1D((TString)("JetShape2_"+stem+ datalabel+"_pp_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  
	    JetShapeMC[g][ibin][ibin3] = new TH1D((TString)("JetShapeMC_"+stem+ datalabel+"_pp_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  
	    JetShapeMC_reco[g][ibin][ibin3] = new TH1D((TString)("JetShapeMC_reco_"+stem+ datalabel+"_pp_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  

	    JetShape2_geo[g][ibin][ibin3] = new TH1D((TString)("JetShape_Geometry_"+stem+ datalabel+"_pp_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  

	  }
	  

	  int kmin = result[g][ibin][0][ibin3]->GetXaxis()->FindBin(-TMath::Pi()/2.+.0001);
	  int kmax = result[g][ibin][0][ibin3]->GetXaxis()->FindBin(TMath::Pi()/2.-.0001);
      

	  int lmin = result[g][ibin][0][ibin3]->GetYaxis()->FindBin(-TMath::Pi()/2.+.0001);
	  int lmax = result[g][ibin][0][ibin3]->GetYaxis()->FindBin(TMath::Pi()/2.-.0001);
	
	
	  for(int k = kmin; k<kmax+1; k++){  //
	    for(int l = lmin; l<lmax+1; l++){

	      deta = result[g][ibin][0][ibin3]->GetXaxis()->GetBinCenter(k);
	      dphi = result[g][ibin][0][ibin3]->GetYaxis()->GetBinCenter(l);
	   
	      r = TMath::Sqrt(deta*deta+dphi*dphi);
	      
	      
	      rbin = JetShape2[g][ibin][ibin3]->FindBin(r);

 
	      //catch-all, if you want
	      /*
	      if(rbin > 14){
		rbin = 14;
		cout<<r<<" "<<rbin<<endl;
	      }
	      */
	      bc = result[g][ibin][0][ibin3]->GetBinContent(k,l);
	      err = result[g][ibin][0][ibin3]->GetBinError(k,l);

	      temp1 = JetShape2[g][ibin][ibin3]->GetBinContent(rbin);
	      temp2 = JetShape2_geo[g][ibin][ibin3]->GetBinContent(rbin);
	    
	      temperr = JetShape2[g][ibin][ibin3]->GetBinError(rbin);
	  
	      JetShape2[g][ibin][ibin3]->SetBinContent(rbin,temp1+bc);
	      JetShape2[g][ibin][ibin3]->SetBinError(rbin,TMath::Sqrt(temperr*temperr+err*err));
	      
	      JetShape2_geo[g][ibin][ibin3]->SetBinContent(rbin,temp2+1.);
	      JetShape2_geo[g][ibin][ibin3]->SetBinError(rbin,0.);

	   
	      if(g==0&&(ibin3<4||ibin3==9)){
   
		bc = resultMC[g][ibin][0][ibin3]->GetBinContent(k,l);
		err = resultMC[g][ibin][0][ibin3]->GetBinError(k,l);

		temp1 = JetShapeMC[g][ibin][ibin3]->GetBinContent(rbin);
		temperr = JetShapeMC[g][ibin][ibin3]->GetBinError(rbin);
	  
		JetShapeMC[g][ibin][ibin3]->SetBinContent(rbin,temp1+bc);
		JetShapeMC[g][ibin][ibin3]->SetBinError(rbin,TMath::Sqrt(temperr*temperr+err*err));
	      }



   
	      bc = resultMC_reco[g][ibin][0][ibin3]->GetBinContent(k,l);
	      err = resultMC_reco[g][ibin][0][ibin3]->GetBinError(k,l);

	      temp1 = JetShapeMC_reco[g][ibin][ibin3]->GetBinContent(rbin);
	      temperr = JetShapeMC_reco[g][ibin][ibin3]->GetBinError(rbin);
	  
	      JetShapeMC_reco[g][ibin][ibin3]->SetBinContent(rbin,temp1+bc);
	      JetShapeMC_reco[g][ibin][ibin3]->SetBinError(rbin,TMath::Sqrt(temperr*temperr+err*err));

	  	  
	    }
	  }

	  if(g==0&&(ibin3<4||ibin3==9))	  JetShapeMC_reco[g][ibin][ibin3]->Add(JetShapeMC[g][ibin][ibin3],-1.);  //this step isolates jff correction in jet shape form for errors.
		
	  JetShape[g][ibin][ibin3] = (TH1D*)JetShape2[g][ibin][ibin3]->Clone((TString)("JetShape_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	
	
	  for(int k = 1; k<JetShape2[g][ibin][ibin3]->GetNbinsX()+1; k++){
	    
	    float temp = JetShape2[g][ibin][ibin3]->GetBinContent(k)/JetShape2[g][ibin][ibin3]->GetBinWidth(k);
	    float err = JetShape2[g][ibin][ibin3]->GetBinError(k)/JetShape2[g][ibin][ibin3]->GetBinWidth(k);
	    
	    JetShape2[g][ibin][ibin3]->SetBinContent(k,temp);
	    JetShape2[g][ibin][ibin3]->SetBinError(k,err);

	    JetShape[g][ibin][ibin3]->SetBinContent(k,temp);
	    JetShape[g][ibin][ibin3]->SetBinError(k,err);

	    temp = JetShapeMC[g][ibin][ibin3]->GetBinContent(k)/JetShapeMC[g][ibin][ibin3]->GetBinWidth(k);
	    err = JetShapeMC[g][ibin][ibin3]->GetBinError(k)/JetShapeMC[g][ibin][ibin3]->GetBinWidth(k);

	    JetShapeMC[g][ibin][ibin3]->SetBinContent(k,temp);
	    JetShapeMC[g][ibin][ibin3]->SetBinError(k,err);


	    temp = JetShapeMC_reco[g][ibin][ibin3]->GetBinContent(k)/JetShapeMC[g][ibin][ibin3]->GetBinWidth(k);
	    err = JetShapeMC_reco[g][ibin][ibin3]->GetBinError(k)/JetShapeMC[g][ibin][ibin3]->GetBinWidth(k);

	    JetShapeMC_reco[g][ibin][ibin3]->SetBinContent(k,temp);
	    JetShapeMC_reco[g][ibin][ibin3]->SetBinError(k,err);	    
	   

	  }
	  /*
	  if(is_number){
	    norm = JetShape[g][ibin][ibin3]->GetBinWidth(1);

	    JetShape[g][ibin][ibin3]->Scale(1./norm);  
	    JetShape2[g][ibin][ibin3]->Scale(1./norm);     
	  }
	  */

	  /*
	
	  JetShapeMC[g][ibin][ibin3]->Divide(JetShape2_geo[g][ibin][ibin3]);

	  JetShapeMC[g][ibin][ibin3]->Scale(1./norm);
	  */
	
	 

	  if(g==0)	  JetShape_syst[g][ibin][ibin3] = (TH1D*) JetShape2[g][ibin][ibin3]->Clone((TString)("Jet_Shape_SystErr_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  else 	  JetShape_syst[g][ibin][ibin3] = (TH1D*) JetShape2[g][ibin][ibin3]->Clone((TString)("Jet_Shape_SystErr_pp_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	
	  if(g==0)    me_error = me_err_pbpb->GetBinContent(ibin3+1,ibin+1);
	  else   me_error = me_err_pp->GetBinContent(ibin3+1,ibin+1);

	  if(g==0)    bg_error = bg_err_pbpb->GetBinContent(ibin3+1,ibin+1);
	  else   bg_error = bg_err_pp->GetBinContent(ibin3+1,ibin+1);


	  // if(is_number){
	    if(ibin3==9){
	      if(g==0)    me_error = me_err_pbpb->GetBinContent(1,ibin+1);
	      else   me_error = me_err_pp->GetBinContent(1,ibin+1);

	      if(g==0)    bg_error = bg_err_pbpb->GetBinContent(1,ibin+1);
	      else   bg_error = bg_err_pp->GetBinContent(1,ibin+1);

	      for(int i = 2; i<8; i++){
		if(g==0)    me_error += me_err_pbpb->GetBinContent(i,ibin+1);
		else   me_error += me_err_pp->GetBinContent(i,ibin+1);

		if(g==0)    bg_error += bg_err_pbpb->GetBinContent(i,ibin+1);
		else   bg_error += bg_err_pp->GetBinContent(i,ibin+1);
	      }
	    }

	    //TURN ON TO TAKE BKG AND ME ERRORS FROM NON-PTWEIGHTED
	    /*
	  }else{

	    if(ibin3==9){
	      if(g==0)    me_error = me_err_pbpb->GetBinContent(1,ibin+1)*mean_pts[0];
	      else   me_error = me_err_pp->GetBinContent(1,ibin+1)*mean_pts[0];
		
	      if(g==0)    bg_error = bg_err_pbpb->GetBinContent(1,ibin+1)*mean_pts[0];
	      else   bg_error = bg_err_pp->GetBinContent(1,ibin+1)*mean_pts[0];

	      for(int i = 2; i<9; i++){
		if(g==0)    me_error += me_err_pbpb->GetBinContent(i,ibin+1)*mean_pts[i];
		else   me_error += me_err_pp->GetBinContent(i,ibin+1)*mean_pts[i];

		if(g==0)    bg_error += bg_err_pbpb->GetBinContent(i,ibin+1)*mean_pts[i];
		else   bg_error += bg_err_pp->GetBinContent(i,ibin+1)*mean_pts[i];
	      }
	    }

	  }
	    */

	  for(int k = 1; k<JetShape_syst[g][ibin][ibin3]->GetNbinsX()+1; k++){

	    
	    bc =JetShape_syst[g][ibin][ibin3]->GetBinContent(k);

	    if(g==0)  mc_error = TMath::Sqrt((JetShapeMC[g][ibin][ibin3]->GetBinContent(k)*JetShapeMC[g][ibin][ibin3]->GetBinContent(k)*.18*.18)+(JetShapeMC_reco[g][ibin][ibin3]->GetBinContent(k)*JetShapeMC_reco[g][ibin][ibin3]->GetBinContent(k)*.5*.5));	 
	    else   mc_error = JetShapeMC_reco[g][ibin][ibin3]->GetBinContent(k)*.5;

	    cout<<"MC ERRORS: "<<g<<" "<<ibin<<" "<<ibin3<<" "<<k<<" "<<JetShapeMC[g][ibin][ibin3]->GetBinContent(k)<<" "<<JetShapeMC_reco[g][ibin][ibin3]->GetBinContent(k)<<" "<<mc_error<<endl;
   
	  
	    //phase factor of the area of the dR ring is pi^r^2, divide by +/-1 integration range.
	    bg_error2 = TMath::Pi()*((JetShape_syst[g][ibin][ibin3]->GetBinLowEdge(k+1)*JetShape_syst[g][ibin][ibin3]->GetBinLowEdge(k+1))-JetShape_syst[g][ibin][ibin3]->GetBinLowEdge(k)*JetShape_syst[g][ibin][ibin3]->GetBinLowEdge(k))*bg_error/2.;

	    me_error2 = TMath::Pi()*((JetShape_syst[g][ibin][ibin3]->GetBinLowEdge(k+1)*JetShape_syst[g][ibin][ibin3]->GetBinLowEdge(k+1))-JetShape_syst[g][ibin][ibin3]->GetBinLowEdge(k)*JetShape_syst[g][ibin][ibin3]->GetBinLowEdge(k))*me_error/2.;
	    
	    //   cout<<"Error is: "<<bg_error<<" "<<me_error<<endl;

	    if(g==0)	JetShape_syst[g][ibin][ibin3]->SetBinError(k,TMath::Sqrt(bc*bc*rel_err_pbpb*rel_err_pbpb+mc_error*mc_error+bg_error2*bg_error2+me_error2*me_error2)); 
	    else   JetShape_syst[g][ibin][ibin3]->SetBinError(k,TMath::Sqrt(bc*bc*rel_err_pp*rel_err_pp+mc_error*mc_error+bg_error2*bg_error2+me_error2*me_error2)); 

	    //JetShape_syst[g][ibin][ibin3]->SetBinError(k,TMath::Sqrt(bc*bc*(.05*.05+0.04*0.04+0.03*0.03)+bg_error2*bg_error2+me_error2*me_error2)); 
	 	      
  }

	  cout<<"done"<<endl;


	  /*
	  if(is_number&&ibin3==9){

	    JetShape_syst[g][ibin][ibin3] = (TH1D*) JetShape_syst[g][ibin][0] ->Clone((TString)("Jet_Shape_Syst_Tot_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    
	    for(int k = 1; k<8; k++){
	      JetShape_syst[g][ibin][ibin3]->Add(JetShape_syst[g][ibin][k] );
	    } 
	  }else
	  */
	  fout->cd();

	  JetShape2[g][ibin][ibin3]->Write();
	  JetShape_syst[g][ibin][ibin3]->Write();


	  if(!is_number&&ibin3==9){
	  
	
	  
	  if(ibin<3){
	    JetShape2[g][ibin][ibin3]->GetYaxis()->SetLabelSize(0.);
	    JetShape2[g][ibin][ibin3]->GetYaxis()->SetTitleSize(0.);

	  }
	    int end_bin = JetShape2[g][ibin][ibin3]->FindBin(.9999);

	    norm =  JetShape2[g][ibin][ibin3]->Integral(1,end_bin,"width");
	    cout<<"norm is "<<norm<<" "<<g<<" "<<ibin<<" "<<ibin3<<endl;

	    //norm = 1.;

	    //	    if(!use_highpT_bin) norm = JetShape[g][ibin][ibin3]->GetBinContent(1)/JetShape_ref[g][ibin][ibin3]->GetBinContent(1);	 
	    for(int k = 0; k<10; k++){
	      JetShape[g][ibin][k]->Scale(1./norm);  
	      JetShape2[g][ibin][k]->Scale(1./norm);  
	      //  JetShapeMC[g][ibin][k]->Scale(1./norm);  
	      JetShape_syst[g][ibin][k]->Scale(1./norm);  
	  

	    }
	  }


	  //	  JetShape_ref[g][ibin][ibin3] = (TH1D*)f_in_ref->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + "_Ref_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  cout<<"and here"<<endl;
	
	
	  if(is_number){

	    int end_bin = JetShape2[g][ibin][ibin3]->FindBin(.999);


	    integral[g][ibin][ibin3] =  JetShape2[g][ibin][ibin3]->IntegralAndError(1,end_bin,integral_err[g][ibin][ibin3],"width");
	 

	    
	    float temp_err  = 0.;

	    for(int k = 1; k<JetShape_syst[g][ibin][ibin3]->GetNbinsX()+1; k++){
	      temp_err+=JetShape_syst[g][ibin][ibin3]->GetBinError(k)*JetShape_syst[g][ibin][ibin3]->GetBinWidth(k);

	    }

	    integral_syst_err[g][ibin][ibin3] = temp_err;
	    cout<<"norm is "<<integral[g][ibin][ibin3]<<" "<<integral_err[g][ibin][ibin3]<<" "<<integral_syst_err[g][ibin][ibin3]<<" "<<g<<" "<<ibin<<" "<<ibin3<<endl;
	  }



	  if(g==0)	  mc_canvas->cd(4*(ibin3+1)-ibin);

	
	  JetShapeMC[g][ibin][ibin3]->SetMinimum(-10.);
	  JetShapeMC[g][ibin][ibin3]->SetMaximum(15.);

	  JetShapeMC[g][ibin][ibin3]->SetMarkerStyle(10);
	  JetShapeMC[g][ibin][ibin3]->SetMarkerSize(1.);
	  JetShapeMC[g][ibin][ibin3]->GetYaxis()->SetLabelSize(0.08);

	  JetShapeMC[g][ibin][ibin3]->GetXaxis()->SetRangeUser(0.,1.);
	  JetShapeMC[g][ibin][ibin3]->Draw();
	
	  JetShapeMC_reco[g][ibin][ibin3]->SetMarkerStyle(10);
	  JetShapeMC_reco[g][ibin][ibin3]->SetMarkerSize(1);
	  JetShapeMC_reco[g][ibin][ibin3]->SetMarkerColor(kRed);
	  JetShapeMC_reco[g][ibin][ibin3]->SetLineColor(kRed);
	  JetShapeMC_reco[g][ibin][ibin3]->Draw("same");

	
	  JetShape_graph[g][ibin][ibin3] = new TGraphAsymmErrors(JetShape_syst[g][ibin][ibin3]);

	  JetShape_graph[g][ibin][ibin3]->SetName((TString)("Jet_Shape_SystErrGraph_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	  JetShape_graph[g][ibin][ibin3] = new TGraphAsymmErrors(JetShape_syst[g][ibin][ibin3]);

	  JetShape_graph[g][ibin][ibin3]->SetName((TString)("Jet_Shape_SystErrGraph_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));



	
	  JetShape_graph[g][ibin][ibin3]->SetMarkerStyle(20);
	  JetShape_graph[g][ibin][ibin3]->SetMarkerSize(1);
	  JetShape_graph[g][ibin][ibin3]->SetMarkerColor(kWhite);
	  JetShape_graph[g][ibin][ibin3]->SetFillColor(kBlack);
	  JetShape_graph[g][ibin][ibin3]->SetFillStyle(3004);

	  if(g==1){

	    for(int k = 0; k<nCBins; k++){
	      JetShape_syst[g+1][k][ibin3] = (TH1D*) JetShape_syst[g-1][k][ibin3]->Clone((TString)("Jet_Shape_Syst_Ratio_"+ CBin_strs[k] + "_" + CBin_strs[k+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      if(is_number)	      JetShape_syst[g+1][k][ibin3]->Add(JetShape_syst[g][0][ibin3],-1.);
	      else 	      JetShape_syst[g+1][k][ibin3]->Divide(JetShape_syst[g][0][ibin3]);


	      JetShape_graph[g+1][k][ibin3] = new TGraphAsymmErrors(JetShape_syst[g+1][k][ibin3]);

	      JetShape_graph[g+1][k][ibin3]->SetName((TString)("Jet_Shape_Ratio_Graph_"+ CBin_strs[k] + "_" + CBin_strs[k+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	     

	      JetShape_graph[g+1][k][ibin3]->SetMarkerStyle(20);
	      JetShape_graph[g+1][k][ibin3]->SetMarkerSize(1);
	      JetShape_graph[g+1][k][ibin3]->SetMarkerColor(kWhite);
	      JetShape_graph[g+1][k][ibin3]->SetFillColor(kBlack);
	      JetShape_graph[g+1][k][ibin3]->SetFillStyle(3004);

	    }
	  }
	  cout<<"and done with loop"<<endl;
	}

      }

    }
  }
 
  if(is_number) mc_canvas->SaveAs("MC_Corr_Test.png");
  else  mc_canvas->SaveAs("MC_Corr_Test_pTweighted.png");
  cout<<"ready to draw"<<endl;

  for(int k = 0; k<nTrkPtBins; k++){
   
    TString canvas_name = "JetShape_Comparison_";
  
    if(k<9){ 
      canvas_name+=TrkPtBin_strs[k]; canvas_name+="_"; canvas_name+=TrkPtBin_strs[k+1];  
    } else {
      canvas_name+=TrkPtBin_strs[1]; canvas_name+="_"; canvas_name+=TrkPtBin_strs[9];  
    }
   
    c_jetshape[k] = new TCanvas(canvas_name," ",10,10,1500,800);
    c_jetshape[k]->Divide(4,2,0.,0.);

  
    for(int ibin = 0; ibin<nCBins; ibin++){
      
      c_jetshape[k]->cd(4-ibin);

      JetShape2[0][ibin][k]->GetXaxis()->SetTitle("Radius(r)");
      JetShape2[0][ibin][k]->GetXaxis()->CenterTitle(true);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelFont(42);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelOffset(0.02);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelSize(0.08);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleSize(0.08);
      JetShape2[0][ibin][k]->GetXaxis()->SetTickLength(0.025);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleOffset(1.3);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleFont(42);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitle("#Rho(#Deltar)");
      JetShape2[0][ibin][k]->GetYaxis()->CenterTitle(true);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelFont(42);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelOffset(0.004);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.075);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitleSize(0.09);
      JetShape2[0][ibin][k]->GetYaxis()->SetTickLength(0.025);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitleOffset(1.);
     
      JetShape2[0][ibin][k]->SetMinimum(0.000008);
      JetShape2[0][ibin][k]->SetMaximum(14.5);

      // JetShape2[0][ibin][k]->SetMinimum(0.5);
      //JetShape2[0][ibin][k]->SetMaximum(1500.1);

      JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.07);
    

      gPad->SetLogy();

      JetShape2[0][ibin][k]->SetMarkerSize(1);
      JetShape2[0][ibin][k]->SetLineColor(kBlack);
      JetShape2[0][ibin][k]->SetMarkerColor(kBlack);
      JetShape2[0][ibin][k]->SetMarkerStyle(24);


      if(ibin<3){
	JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.);
	JetShape2[0][ibin][k]->GetYaxis()->SetTitleSize(0.);
      }

      gPad->SetLogy();

      JetShape2[0][ibin][k]->GetXaxis()->SetRangeUser(0.,0.99);
      JetShape2[0][ibin][k]->Draw();

      
  
      JetShape2[1][0][k]->SetMarkerSize(1);
      JetShape2[1][0][k]->SetLineColor(kBlue);
      JetShape2[1][0][k]->SetMarkerColor(kBlue);
      JetShape2[1][0][k]->SetMarkerStyle(24);
      

      JetShape2[1][0][k]->Draw("same");

  

      if(ibin==3){
	TLegend *legend = new TLegend(0.4,0.7,0.95,0.85);
	//	if(((k>0&&k<5)||k==6)&&!is_subleading)	legend->AddEntry(JetShape_ref[0][ibin][k],"PbPb JetShapes");
	legend->AddEntry(JetShape2[0][ibin][k],"PbPb");
	//if(((k>0&&k<5)||k==6)&&!is_subleading)	legend->AddEntry(JetShape_ref[1][0][k],"pp JetShapes");
	legend->AddEntry(JetShape2[1][0][k],"pp");
	legend->SetLineColor(kWhite);
	legend->SetFillColor(kWhite);
	legend->SetTextSize(0.065);
	legend->Draw();
      }else if(ibin==2){

	TLatex  *pttex = new TLatex(0.5,0.81,TrkPtBin_labels[k]);
	pttex->SetNDC();
	pttex->SetTextSize(0.075);
	pttex->Draw();

      }

    
      TLatex  *centtex ;
      centtex = new TLatex(0.55,0.9,CBin_labels[ibin]);
      centtex->SetNDC();
      centtex->SetTextSize(0.075);
      if(ibin==3){  
	centtex = new TLatex(0.6,0.9,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->SetTextSize(0.07);
      }
      
      centtex->Draw();

      c_jetshape[k]->cd(8-ibin);
  
      
    
      TString diff_name ="JetShape_diff_"; diff_name+=CBin_strs[ibin]; diff_name+= "_"; diff_name += CBin_strs[ibin+1]; diff_name+= k;

      TString diff2_name ="JetShape_diff2_"; diff2_name+=CBin_strs[ibin]; diff2_name+= "_"; diff2_name += CBin_strs[ibin+1]; diff2_name+= k;
      
     
      JetShape_diff[0][ibin][k] = (TH1D*)JetShape[0][ibin][k]->Clone(diff_name);
      JetShape_diff2[0][ibin][k] = (TH1D*)JetShape2[0][ibin][k]->Clone(diff2_name);

      JetShape_diff[0][ibin][k]->Add( JetShape[1][0][k],-1. );
      JetShape_diff2[0][ibin][k]->Add( JetShape2[1][0][k],-1. );
   

      JetShape_diff2[0][ibin][k]->SetMinimum(-2);
      JetShape_diff2[0][ibin][k]->SetMaximum(4.);
     

      JetShape_diff2[0][ibin][k]->SetMarkerStyle(24);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitle("#rho(r)_{PbPb} - #rho(r)_{pp}");
      JetShape_diff2[0][ibin][k]->GetXaxis()->SetNdivisions(408);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetNdivisions(612);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetLabelSize(0.07);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitleSize(0.09);

      if(ibin==3){
	JetShape_diff2[0][ibin][k]->GetXaxis()->SetTitleSize(0.9*JetShape_diff2[0][ibin][k]->GetXaxis()->GetTitleSize());
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitleSize(0.9*JetShape_diff2[0][ibin][k]->GetYaxis()->GetTitleSize());
	JetShape_diff2[0][ibin][k]->GetXaxis()->SetLabelSize(0.9*JetShape_diff2[0][ibin][k]->GetXaxis()->GetLabelSize());
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetLabelSize(0.9*JetShape_diff2[0][ibin][k]->GetYaxis()->GetLabelSize());
      }

      if(ibin<3){
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetLabelSize(0.);
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitleSize(0.);
      }

      JetShape_diff2[0][ibin][k]->SetAxisRange(0.,.99);
      JetShape_diff2[0][ibin][k]->Draw();


      TString ratio_name ="JetShape_ratio_"; ratio_name+=CBin_strs[ibin]; ratio_name+= "_"; ratio_name += CBin_strs[ibin+1]; ratio_name+= k;

      TString ratio2_name ="JetShape_ratio2_"; ratio2_name+=CBin_strs[ibin]; ratio2_name+= "_"; ratio2_name += CBin_strs[ibin+1]; ratio2_name+= k;
      
     
      JetShape_ratio[0][ibin][k] = (TH1D*)JetShape[0][ibin][k]->Clone(ratio_name);
      JetShape_ratio2[0][ibin][k] = (TH1D*)JetShape2[0][ibin][k]->Clone(ratio2_name);

      if(!is_number){ 
	JetShape_ratio[0][ibin][k]->Divide( JetShape[1][0][k]);
	JetShape_ratio2[0][ibin][k]->Divide( JetShape2[1][0][k] );
      
      }else{
	JetShape_ratio[0][ibin][k]->Add( JetShape[1][0][k],-1.);
	JetShape_ratio2[0][ibin][k]->Add( JetShape2[1][0][k],-1.);
      }

      JetShape_ratio2[0][ibin][k]->SetMinimum(-1.5);
      JetShape_ratio2[0][ibin][k]->SetMaximum(9.);
     

      JetShape_ratio2[0][ibin][k]->SetMarkerStyle(20);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitle("#Rho(r)_{PbPb} / #Rho(r)_{pp}");
      JetShape_ratio2[0][ibin][k]->GetXaxis()->SetNdivisions(408);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetNdivisions(612);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetLabelSize(0.07);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitleSize(0.08);

      if(ibin==3){
	JetShape_ratio2[0][ibin][k]->GetXaxis()->SetTitleSize(0.9*JetShape_ratio2[0][ibin][k]->GetXaxis()->GetTitleSize());
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitleSize(0.9*JetShape_ratio2[0][ibin][k]->GetYaxis()->GetTitleSize());
	JetShape_ratio2[0][ibin][k]->GetXaxis()->SetLabelSize(0.9*JetShape_ratio2[0][ibin][k]->GetXaxis()->GetLabelSize());
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetLabelSize(0.9*JetShape_ratio2[0][ibin][k]->GetYaxis()->GetLabelSize());
      }

      if(ibin<3){
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetLabelSize(0.);
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitleSize(0.);
      }

      JetShape_ratio2[0][ibin][k]->SetAxisRange(0.,.99);
      JetShape_ratio2[0][ibin][k]->Draw();


      TLine *l = new TLine(0.,1.,2.,1.);
      l->SetLineStyle(2);
      l->Draw("same");
      /*
	ratio_name ="JetShape_ref_ratio_"; ratio_name+=CBin_strs[ibin]; ratio_name+= "_"; ratio_name += CBin_strs[ibin+1]; ratio_name+= k;
      
     
	JetShape_ref_ratio[0][ibin][k] = (TH1D*)JetShape_ref[0][ibin][k]->Clone(ratio_name);

	JetShape_ref_ratio[0][ibin][k]->Divide( JetShape_ref[1][0][k] );

	JetShape_ref_ratio[0][ibin][k]->SetMarkerStyle(20);
	JetShape_ref_ratio[0][ibin][k]->SetMarkerColor(kRed);
	JetShape_ref_ratio[0][ibin][k]->Draw("same");
      */
         
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelSize(0.);

     
      if(ibin==3){
	TLegend *legend = new TLegend(0.2,0.7,0.4,0.85);
	//	if(((k>0&&k<5)||k==6)&&!is_subleading)	legend->AddEntry(JetShape_ref_ratio[0][ibin][k],"PbPb - pp JetShapes");
	legend->AddEntry(JetShape_ratio2[0][ibin][k],"PbPb / pp");
	legend->SetLineColor(kWhite);
	legend->SetFillColor(kWhite);
	legend->SetTextSize(0.06);
	legend->Draw();
      }      


    }


    canvas_name+=".pdf";
    
    c_jetshape[k]->SaveAs(canvas_name);
    
    canvas_name.ReplaceAll(".pdf",".png");
    c_jetshape[k]->SaveAs(canvas_name);


  }

  for(int ibin = 0; ibin<nCBins; ibin++){
    
    for(int k = 0; k<nTrkPtBins; k++){

      
 
      JetShape_noerr_up[0][ibin][k] = new TH1D((TString)("JetShape_PbPb_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
      JetShape_noerr_down[0][ibin][k] = new TH1D((TString)("JetShape_PbPb_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
   
      if(ibin==0){
	JetShape_noerr_up[1][ibin][k] = new TH1D((TString)("JetShape_pp_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
	JetShape_noerr_down[1][ibin][k] = new TH1D((TString)("JetShape_pp_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
      }
    

      JetShape_diff_noerr_up[0][ibin][k] = new TH1D((TString)("JetShape_Diff_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
      JetShape_diff_noerr_down[0][ibin][k] = new TH1D((TString)("JetShape_Diff_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
   
      

      for(int l = 0; l< JetShape[0][ibin][k]->GetNbinsX()+1; l++){
      
	bc = JetShape[0][ibin][k]->GetBinContent(l);
          
	if(bc>0){ 
	  JetShape_noerr_up[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_noerr_down[0][ibin][k]->SetBinContent(l,0.);	

	}else{
	  JetShape_noerr_down[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_noerr_up[0][ibin][k]->SetBinContent(l,0.);	
	}

  
	bc = JetShape_diff[0][ibin][k]->GetBinContent(l);

    
	if(bc>0){ 
	  JetShape_diff_noerr_up[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_diff_noerr_down[0][ibin][k]->SetBinContent(l,0.);	

	}else{
	  JetShape_diff_noerr_down[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_diff_noerr_up[0][ibin][k]->SetBinContent(l,0.);	
	} 

	if(ibin==0){

	  bc = JetShape[1][ibin][k]->GetBinContent(l);
      
	  if(bc>0){ 
	    JetShape_noerr_up[1][ibin][k]->SetBinContent(l,bc);	
	    JetShape_noerr_down[1][ibin][k]->SetBinContent(l,0.);	

	  }else{
	    JetShape_noerr_down[1][ibin][k]->SetBinContent(l,bc);	
	    JetShape_noerr_up[1][ibin][k]->SetBinContent(l,0.);	
	  }
	}
      
      }


      // cout<<"Have set for "<<ibin<<" "<<k<<" "<<JetShape_noerr_up[0][ibin][k]->GetBinContent(1)<<" "<<JetShape_noerr_up[0][ibin][k]->GetBinCenter(1)<<endl;

      JetShape_noerr_up[0][ibin][k]->SetLineColor(kBlack);;
      JetShape_noerr_up[1][0][k]->SetLineColor(kBlack);;
      JetShape_noerr_down[0][ibin][k]->SetLineColor(kBlack);;
      JetShape_noerr_down[1][0][k]->SetLineColor(kBlack);;
      JetShape_diff_noerr_up[0][ibin][k]->SetLineColor(kBlack);;
      JetShape_diff_noerr_down[0][ibin][k]->SetLineColor(kBlack);;

      switch(k){
      case 0:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kBlue-9);
	JetShape_noerr_up[1][0][k]->SetFillColor(kBlue-9);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kBlue-9);
	JetShape_noerr_down[1][0][k]->SetFillColor(kBlue-9);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kBlue-9);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kBlue-9);
	break;
      case 1:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_noerr_up[1][0][k]->SetFillColor(kYellow-9);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_noerr_down[1][0][k]->SetFillColor(kYellow-9);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kYellow-9);
	break;
      case 2:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_noerr_up[1][0][k]->SetFillColor(kOrange+1);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_noerr_down[1][0][k]->SetFillColor(kOrange+1);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kOrange+1);
	break;
      case 3:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_noerr_up[1][0][k]->SetFillColor(kViolet-5);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_noerr_down[1][0][k]->SetFillColor(kViolet-5);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kViolet-5);
	break;
      case 4:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_noerr_up[1][0][k]->SetFillColor(kGreen+3);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_noerr_down[1][0][k]->SetFillColor(kGreen+3);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kGreen+3);
	break;


      case 5:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kRed);
	JetShape_noerr_up[1][0][k]->SetFillColor(kRed);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kRed);
	JetShape_noerr_down[1][0][k]->SetFillColor(kRed);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kRed);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kRed);
	break;

    
      case 6:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_noerr_up[1][0][k]->SetFillColor(kRed+1);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_noerr_down[1][0][k]->SetFillColor(kRed+1);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kRed+1);
	break;

      case 7:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kRed+2);
	JetShape_noerr_up[1][0][k]->SetFillColor(kRed+2);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kRed+2);
	JetShape_noerr_down[1][0][k]->SetFillColor(kRed+2);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kRed+2);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kRed+2);
	break;

      case 8:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kRed+3);
	JetShape_noerr_up[1][0][k]->SetFillColor(kRed+3);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kRed+3);
	JetShape_noerr_down[1][0][k]->SetFillColor(kRed+3);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kRed+3);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kRed+3);
	break;

      default:
	break;
      }

      JetShape_noerr_up[0][ibin][k]->SetFillStyle(1001);
      JetShape_noerr_up[1][0][k]->SetFillStyle(1001);
      JetShape_diff_noerr_up[0][ibin][k]->SetFillStyle(1001);

     	  
      JetShape_noerr_down[0][ibin][k]->SetFillStyle(1001);
      JetShape_noerr_down[1][0][k]->SetFillStyle(1001);
      JetShape_diff_noerr_down[0][ibin][k]->SetFillStyle(1001);


      
    } //k
    

    JetShape_Stack_Up[0][ibin]=new THStack((TString)("JetShapeStack_PbPb_Up_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Stack_Up[1][0]=new THStack((TString)("JetShapeStack_pp_Up_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Diff_Stack_Up[0][ibin]=new THStack((TString)("JetShapeStack_Diff_Up_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");


    JetShape_Stack_Down[0][ibin]=new THStack((TString)("JetShapeStack_PbPb_Down_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Stack_Down[1][0]=new THStack((TString)("JetShapeStack_pp_Down_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Diff_Stack_Down[0][ibin]=new THStack((TString)("JetShapeStack_Diff_Down_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");

    if(is_number){
      for(int k = 0; k<nTrkPtBins-2; k++){
	JetShape_Stack_Up[0][ibin]->Add(JetShape_noerr_up[0][ibin][nTrkPtBins-3-k]);
	JetShape_Stack_Up[1][0]->Add(JetShape_noerr_up[1][0][nTrkPtBins-3-k]);
	JetShape_Diff_Stack_Up[0][ibin]->Add(JetShape_diff_noerr_up[0][ibin][nTrkPtBins-3-k]);

	JetShape_Stack_Down[0][ibin]->Add(JetShape_noerr_down[0][ibin][nTrkPtBins-3-k]);
	JetShape_Stack_Down[1][0]->Add(JetShape_noerr_down[1][0][nTrkPtBins-3-k]);
	JetShape_Diff_Stack_Down[0][ibin]->Add(JetShape_diff_noerr_down[0][ibin][nTrkPtBins-3-k]);
      }

    }else{

      for(int k = 0; k<nTrkPtBins-1; k++){
	JetShape_Stack_Up[0][ibin]->Add(JetShape_noerr_up[0][ibin][k]);
	JetShape_Stack_Up[1][0]->Add(JetShape_noerr_up[1][0][k]);
	JetShape_Diff_Stack_Up[0][ibin]->Add(JetShape_diff_noerr_up[0][ibin][k]);

	JetShape_Stack_Down[0][ibin]->Add(JetShape_noerr_down[0][ibin][k]);
	JetShape_Stack_Down[1][0]->Add(JetShape_noerr_down[1][0][k]);
	JetShape_Diff_Stack_Down[0][ibin]->Add(JetShape_diff_noerr_down[0][ibin][k]);
      }


    }
  
    if(!is_number){

      JetShape_Stack_Up[0][ibin]->SetMaximum(23.5);
      JetShape_Stack_Up[1][0]->SetMaximum(23.5);
      JetShape_Stack_Up[0][ibin]->SetMinimum(.005);
      JetShape_Stack_Up[1][0]->SetMinimum(.005);
  
    }else{

      JetShape_Stack_Up[0][ibin]->SetMaximum(40.);
      JetShape_Stack_Up[1][0]->SetMaximum(40.);
      JetShape_Stack_Up[0][ibin]->SetMinimum(-1.);
      JetShape_Stack_Up[1][0]->SetMinimum(-1.);
  
       


    }


  } //ibin
  
 


  TCanvas *PAS_plot = new TCanvas("JetShape_ForPAS","",2500,1000);
  PAS_plot->Divide(5,2,0.,0.);

  PAS_plot->cd(6);
  TLegend *legend = new TLegend(0.1,0.05,0.9,.9);
  legend->AddEntry(JetShape_noerr_up[0][0][0],"0.7 < p_{T}^{assoc.}< 1 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][1],"1 < p_{T}^{assoc.}< 2 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][2],"2 < p_{T}^{assoc.}< 3 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][3],"3 < p_{T}^{assoc.}< 4 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][4],"4 < p_{T}^{assoc.}< 8 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][5],"8 < p_{T}^{assoc.}< 12 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][6],"12 < p_{T}^{assoc.}< 16 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][7],"16 < p_{T}^{assoc.}< 20 GeV","f");
 

  if(!use_highpT_bin||is_number){
    legend->AddEntry(JetShape_graph[0][0][9],"Total 0.7 < p_{T}^{assoc.}< 20 GeV","lpfe");
  }else{ 
    legend->AddEntry(JetShape_noerr_up[0][0][8],"20 < p_{T}^{assoc.}< 300 GeV","f");
    legend->AddEntry(JetShape_graph[0][0][9],"Total 0.7 < p_{T}^{assoc.} < 300 GeV","lpfe");
  }
  legend->AddEntry(JetShape_noerr_up[0][0][4],"|#eta_{track}|< 2.4","");
  legend->SetTextSize(0.055);
  legend->SetLineColor(kWhite);
  legend->SetFillColor(kWhite);
  legend->Draw();

  
    

  PAS_plot->cd(0);
  
  TLatex *type_tex;
  if(is_number)type_tex = new TLatex(0.05,0.92,"Particle Yield by #Deltar");
  else  type_tex = new TLatex(0.05,0.92,"Inclusive Jet Shape");
 
  type_tex->SetTextSize(0.035);
  type_tex->SetLineColor(kWhite);
  type_tex->SetNDC();
  type_tex->Draw();
   
  TLatex   *luminosity_tex_pp = new TLatex(0.2,0.92,"pp 25 pb^{-1} (5.02 TeV)");
  luminosity_tex_pp->SetTextFont(43);
  luminosity_tex_pp->SetTextSizePixels(35);
  luminosity_tex_pp->SetLineColor(kWhite);
  luminosity_tex_pp->SetNDC();
  luminosity_tex_pp->Draw();
 
  TLatex   *luminosity_tex_PbPb = new TLatex(0.4,0.92,"PbPb 404 #mub^{-1} (5.02 TeV)");
  luminosity_tex_PbPb->SetTextFont(43);
  luminosity_tex_PbPb->SetTextSizePixels(35);
  luminosity_tex_PbPb->SetLineColor(kWhite);
  luminosity_tex_PbPb->SetNDC();
  luminosity_tex_PbPb->Draw();
 
  TLatex   *jet_reco_tex = new TLatex(0.6,0.92,"ak4CaloJets, p_{T}> 120, |#eta_{jet}| < 1.6");
  jet_reco_tex->SetTextFont(43);
  jet_reco_tex->SetTextSizePixels(35);
  jet_reco_tex->SetLineColor(kWhite);
  jet_reco_tex->SetNDC();
  jet_reco_tex->Draw();


 
    PAS_plot->cd(1);
  
 
  JetShape_Stack_Up[1][0]->Draw();
  JetShape_Stack_Up[1][0]->GetXaxis()->SetRangeUser(0.,0.99);
  JetShape_Stack_Up[1][0]->GetXaxis()->SetNdivisions(505);

   if(use_highpT_bin&&!is_number)  gPad->SetLogy();
  JetShape_Stack_Up[1][0]->GetYaxis()->SetLabelSize(0.08);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleSize(0.09);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleOffset(.7);
  if(is_number) JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("Y = #frac{1}{N_{jet}} #frac{dN}{d#Deltar}");
  else JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("     #rho(#Deltar)");
JetShape_Stack_Up[1][0]->GetYaxis()->CenterTitle();
  JetShape_Stack_Down[1][0]->Draw("same");
  JetShape_Stack_Up[1][0]->Draw("same");

 

JetShape_graph[1][0][9]->Draw("same e2 P");
  JetShape2[1][0][9]->Draw("same");
  if(do_ref){
    JetShape_ref[1][0][9]->SetMarkerColor(kBlue);
    JetShape_ref[1][0][9]->SetMarkerStyle(20);
    JetShape_ref[1][0][9]->SetMarkerSize(1);
    JetShape_ref[1][0][9]->SetLineColor(kBlue);
    JetShape_ref[1][0][9]->Draw("same");
  }  
  labels = new TPaveText(0.2,0.8,0.45,.99,"NDC");
    
  labels->SetName("labels");
  labels->SetFillColor(0);
  labels->SetLineColor(0);
  labels->SetTextAlign(11);
  labels->AddText("pp reference");
  labels->SetTextSize(0.08);
  labels->Draw("same");


  TLine *l_dr = new TLine(0.,0.,TMath::Pi()/2.,0.);
  l_dr->SetLineStyle(2);
  l_dr->Draw();


  TLatex *cms_tex_dphi = new TLatex(0.25,0.8,"CMS");
  cms_tex_dphi->SetTextFont(63);
  cms_tex_dphi->SetTextSizePixels(35);
  cms_tex_dphi->SetLineColor(kWhite);
  cms_tex_dphi->SetNDC();
  cms_tex_dphi->Draw(); 


TLatex *prelim_tex_dphi = new TLatex(0.4,0.8,"Preliminary");
  prelim_tex_dphi->SetTextFont(53);
  prelim_tex_dphi->SetTextSizePixels(35);
  prelim_tex_dphi->SetLineColor(kWhite);
  prelim_tex_dphi->SetNDC();
  prelim_tex_dphi->Draw(); 


  gPad->RedrawAxis();
  

  for(int ibin = 0; ibin<4; ibin++){
    PAS_plot->cd(5-ibin);


  
    JetShape_Stack_Up[0][ibin]->Draw();

 
    JetShape_Stack_Up[0][ibin]->GetXaxis()->SetRangeUser(0.,0.99);
    JetShape_Stack_Up[0][ibin]->GetXaxis()->SetNdivisions(505);

    if(!is_number) gPad->SetLogy();
    JetShape_Stack_Up[0][ibin]->GetYaxis()->SetLabelSize(0.);
    JetShape_Stack_Down[0][ibin]->Draw("same");

 

    JetShape_graph[0][ibin][9]->Draw("same e2 P");
    JetShape2[0][ibin][9]->Draw("same");

    if(do_ref){
      JetShape_ref[0][ibin][9]->SetMarkerColor(kBlue);
      JetShape_ref[0][ibin][9]->SetMarkerStyle(20);
      JetShape_ref[0][ibin][9]->SetMarkerSize(1);
      JetShape_ref[0][ibin][9]->SetLineColor(kBlue);
      JetShape_ref[0][ibin][9]->Draw("same");
    }  

    TPaveText  *labels = new TPaveText(0.05,0.8,0.45,.99,"NDC");
    
    labels->SetName("labels");
    labels->SetFillColor(0);
    labels->SetLineColor(0);
    labels->SetTextAlign(11);
    labels->AddText((TString)("PbPb "+CBin_labels[ibin]));
    labels->SetTextSize(0.08);
    labels->Draw("same");

 
    PAS_plot->cd(10-ibin);

    if(is_number){
      JetShape_ratio[0][ibin][9]->SetMinimum(-5.);
      JetShape_ratio[0][ibin][9]->SetMaximum(20.);

    }else{
      JetShape_ratio[0][ibin][9]->SetMinimum(0.);
      JetShape_ratio[0][ibin][9]->SetMaximum(3.5);
    }

    JetShape_ratio[0][ibin][9]->GetXaxis()->SetRangeUser(0.,0.99);
    JetShape_ratio[0][ibin][9]->Draw();
    JetShape_ratio[0][ibin][9]->GetYaxis()->SetLabelSize(0.08); 
    JetShape_ratio[0][ibin][9]->GetYaxis()->CenterTitle();
    if(is_number) JetShape_ratio[0][ibin][9]->GetYaxis()->SetTitle("Y_{PbPb} - Y_{pp}");
    else JetShape_ratio[0][ibin][9]->GetYaxis()->SetTitle("#rho(r)_{PbPb}/#rho(r)_{pp}");
    JetShape_ratio[0][ibin][9]->GetYaxis()->SetTitleSize(0.0);
    JetShape_ratio[0][ibin][9]->GetYaxis()->SetLabelSize(0.0);
    JetShape_ratio[0][ibin][9]->GetYaxis()->SetTitleOffset(1.);

    JetShape_ratio[0][ibin][9]->GetXaxis()->SetLabelSize(0.08); 
    JetShape_ratio[0][ibin][9]->GetXaxis()->CenterTitle();
    JetShape_ratio[0][ibin][9]->GetXaxis()->SetTitle("#Deltar");
    JetShape_ratio[0][ibin][9]->GetXaxis()->SetTitleSize(0.09);
    JetShape_ratio[0][ibin][9]->GetXaxis()->SetTitleOffset(0.6);
    JetShape_ratio[0][ibin][9]->GetXaxis()->CenterTitle();
    JetShape_ratio[0][ibin][9]->GetXaxis()->SetNdivisions(505);
 

    JetShape_ratio[0][ibin][9]->SetMarkerStyle(20);
    JetShape_ratio[0][ibin][9]->SetMarkerSize(1);
    JetShape_ratio[0][ibin][9]->SetMarkerColor(kBlack);
 
    if(is_number){
      JetShape_Diff_Stack_Up[0][ibin]->Draw("same");
      JetShape_Diff_Stack_Down[0][ibin]->Draw("same");
      JetShape_ratio[0][ibin][9]->Draw("same");
    }

    JetShape_graph[2][ibin][9]->Draw("same e2 P");
    if(is_number) JetShape_ratio[0][ibin][9]->SetMarkerStyle(24);
    JetShape_ratio[0][ibin][9]->Draw("same");
    if(do_ref){
      JetShape_ref_ratio[0][ibin][9]->SetMarkerColor(kBlue);
      JetShape_ref_ratio[0][ibin][9]->SetLineColor(kBlue);
      JetShape_ref_ratio[0][ibin][9]->SetMarkerStyle(20);
      JetShape_ref_ratio[0][ibin][9]->SetMarkerSize(1);
      JetShape_ref_ratio[0][ibin][9]->Draw("same");
    }
    /*

    cout<<"here"<<endl;

    PAS_plot->cd(15-ibin);

    if(is_number){
      JetShape_diff[0][ibin][9]->SetMinimum(-10.);
      JetShape_diff[0][ibin][9]->SetMaximum(25.);

    }else{
      JetShape_diff[0][ibin][9]->SetMinimum(0.);
      JetShape_diff[0][ibin][9]->SetMaximum(3.5);
    }

    JetShape_diff[0][ibin][9]->GetXaxis()->SetRangeUser(0.,0.99);
    JetShape_diff[0][ibin][9]->Draw();
    JetShape_diff[0][ibin][9]->GetYaxis()->SetLabelSize(0.08); 
    JetShape_diff[0][ibin][9]->GetYaxis()->CenterTitle();
    if(is_number) JetShape_diff[0][ibin][9]->GetYaxis()->SetTitle("Y_{PbPb} - Y_{pp}");
    else JetShape_diff[0][ibin][9]->GetYaxis()->SetTitle("#rho(r)_{PbPb}/#rho(r)_{pp}");
    JetShape_diff[0][ibin][9]->GetYaxis()->SetTitleSize(0.0);
    JetShape_diff[0][ibin][9]->GetYaxis()->SetLabelSize(0.0);
    JetShape_diff[0][ibin][9]->GetYaxis()->SetTitleOffset(1.);

    JetShape_diff[0][ibin][9]->GetXaxis()->SetLabelSize(0.08); 
    JetShape_diff[0][ibin][9]->GetXaxis()->CenterTitle();
    JetShape_diff[0][ibin][9]->GetXaxis()->SetTitle("#Deltar");
    JetShape_diff[0][ibin][9]->GetXaxis()->SetTitleSize(0.09);
    JetShape_diff[0][ibin][9]->GetXaxis()->SetTitleOffset(0.6);
    JetShape_diff[0][ibin][9]->GetXaxis()->CenterTitle();
    JetShape_diff[0][ibin][9]->GetXaxis()->SetNdivisions(505);
 

    JetShape_diff[0][ibin][9]->SetMarkerStyle(20);
    JetShape_diff[0][ibin][9]->SetMarkerSize(1);
    JetShape_diff[0][ibin][9]->SetMarkerColor(kBlack);
 

    JetShape_diff[0][ibin][9]->Draw("same");

    cout<<"and here"<<endl;
    
    if(do_ref){
      JetShape_ref_diff[0][ibin][9]->SetMarkerColor(kBlue);
      JetShape_ref_diff[0][ibin][9]->SetLineColor(kBlue);
      JetShape_ref_diff[0][ibin][9]->SetMarkerStyle(20);
      JetShape_ref_diff[0][ibin][9]->SetMarkerSize(1);
      JetShape_ref_diff[0][ibin][9]->Draw("same");
    }
    */

    if(ibin==3){
      TLegend *legend_ratio = new TLegend(0.02,0.7,0.9,0.9);
      if(is_number)legend_ratio->AddEntry(JetShape_ratio[0][ibin][9],"PbPb - pp");
      else legend_ratio->AddEntry(JetShape_ratio[0][ibin][9],"PbPb / pp");
      if(do_ref) legend_ratio->AddEntry(JetShape_ref_ratio[0][ibin][9],"2.76 TeV Leading Jets (50-100%, 0-30%)");
      legend_ratio->SetLineColor(kWhite);
      legend_ratio->SetFillColor(kWhite);
      legend_ratio->SetTextSize(0.08);
      legend_ratio->Draw("same");
    }


    if(ibin < 2&&do_ref){
      TLegend *legend_ratio = new TLegend(0.02,0.8,0.9,0.9);
    
      if(do_ref) legend_ratio->AddEntry(JetShape_ref_ratio[0][ibin][9],"2.76 TeV 0-30%");
      legend_ratio->SetLineColor(kWhite);
      legend_ratio->SetFillColor(kWhite);
      legend_ratio->SetTextSize(0.08);
      legend_ratio->Draw("same");
    }


    if(!is_number){
      TLine *line = new TLine(0.,1.,1.,1.);
      line->SetLineStyle(2);
      line->Draw();
    }else{
      TLine *line = new TLine(0.,0.,1.,0.);
      line->SetLineStyle(2);
      line->Draw();

    }
    /*
    TPave *cover_x_r = new TPave(0.9,0.,1.3,0.11);
    cover_x_r->SetFillColor(kWhite);
    cover_x_r->SetLineColor(kWhite);
    cover_x_r->SetOption("NDC NB");

    //   cover_x_r->Draw();

    TPave *cover_x_l = new TPave(-.03,0.,0.1,0.11);
    cover_x_l->SetFillColor(kWhite);
    cover_x_l->SetLineColor(kWhite);
    cover_x_l->SetOption("NDC NB");


    //  cover_x_l->Draw();
    */
  }

   PAS_plot->cd(6);
   if(is_number){
     TGaxis *dummy_axis_jetshape = new TGaxis(1.,0.13,1.0,.975,-5.,20.);

     dummy_axis_jetshape->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetYaxis());
     dummy_axis_jetshape->SetTitleOffset(1.);
     dummy_axis_jetshape->SetTickSize(0.);
     dummy_axis_jetshape->CenterTitle();
     dummy_axis_jetshape->SetTitleSize(0.08);
     dummy_axis_jetshape->SetLabelSize(0.08);
     dummy_axis_jetshape->SetTitle("Y_{PbPb} - Y_{pp}");
     dummy_axis_jetshape->Draw();
   }else{
     TGaxis *dummy_axis_jetshape = new TGaxis(1.,0.13,1.0,.975,0.,3.5);
     dummy_axis_jetshape->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetYaxis());
     dummy_axis_jetshape->SetTitleOffset(1.2);
     dummy_axis_jetshape->SetTickSize(0.);
     dummy_axis_jetshape->CenterTitle();
     dummy_axis_jetshape->SetTitleSize(0.08);
     dummy_axis_jetshape->SetLabelSize(0.06);
     dummy_axis_jetshape->SetTitle("#rho(#Deltar)_{PbPb}/#rho(#Deltar)_{pp}");
     dummy_axis_jetshape->Draw();
   }
 
  TGaxis *dummy_axis_r = new TGaxis(0.,1.,1.,1.,0.,1.);

  dummy_axis_r->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetXaxis());
  dummy_axis_r->SetTitleOffset(0.6);
  dummy_axis_r->SetTitleSize(0.09);
  dummy_axis_r->SetTickSize(0.);
  dummy_axis_r->SetNdivisions(505);
  //dummy_axis_r->Draw();
  /*
TPave *cover_x_b = new TPave(0.9,0.,0.995,0.17);
   cover_x_b->SetFillColor(kWhite);
   cover_x_b->SetLineColor(kWhite);
   cover_x_b->SetOption("NDC NB");
   //  cover_x_b->Draw();
   */
   PAS_plot->cd(4);
   //   dummy_axis_jetshape->Draw();
   // dummy_axis_r->Draw();
 

   if(!use_highpT_bin){
     PAS_plot->SaveAs("JetShapes_PAS.pdf");
     PAS_plot->SaveAs("JetShapes_PAS.png");
   }else{
     if(!is_number){
       if(do_ref){
	 PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted_WithRun1Ref.pdf");
	 PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted_WithRun1Ref.png");
       }else{
        PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted.pdf");
	PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted.png");


       }
     }else{
       PAS_plot->SaveAs("ParticleYield_by_dR.pdf");
       PAS_plot->SaveAs("ParticleYield_by_dR.png");

     }
   }

   if(is_number){


     for(int ibin = 0; ibin < nCBins; ibin++){

       for(int g = 0; g< 2; g++){

	 if(g==0){

	   Integral_Pt[g][ibin] = new TH1D((TString)("Integral_PbPb"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]),"",9,pTbins);
	   Integral_syst_Pt[g][ibin] = new TH1D((TString)("Integral_PbPb_Syst"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]),"",9,pTbins);
         
	 }else{

	   if(ibin > 0 )continue;

	   Integral_Pt[g][ibin] = new TH1D("Integral_pp","",9,pTbins);
	   Integral_syst_Pt[g][ibin] = new TH1D("Integral_pp_Syst","",9,pTbins);
	 }
     
	 cout<<"here"<<endl;
	 for(int ibin3 = 0; ibin3<nTrkPtBins-2; ibin3++){
      
	   Integral_Pt[g][ibin]->SetBinContent(ibin3+2,integral[g][ibin][ibin3]);
	   Integral_Pt[g][ibin]->SetBinError(ibin3+2,integral_err[g][ibin][ibin3]);
	
	   Integral_syst_Pt[g][ibin]->SetBinContent(ibin3+2,integral[g][ibin][ibin3]);
	   Integral_syst_Pt[g][ibin]->SetBinError(ibin3+2,integral_syst_err[g][ibin][ibin3]);
	   
	 }
	 cout<<"and here"<<endl;
	 if(g==0){
	   Integral_Pt[g][ibin]->SetMarkerStyle(10);
	   Integral_syst_Pt[g][ibin]->SetMarkerStyle(10);
	   Integral_syst_Pt[g][ibin]->SetFillColor(kMagenta-7);
	 }else{
	   Integral_Pt[g][ibin]->SetMarkerStyle(4);
	   Integral_syst_Pt[g][ibin]->SetMarkerStyle(4);
	   Integral_syst_Pt[g][ibin]->SetFillColor(kGray);
	 }
	 Integral_Pt[g][ibin]->SetMarkerSize(2);

	 Integral_Pt[g][ibin]->SetLineColor(kBlack);
	 Integral_Pt[g][ibin]->SetMarkerColor(kBlack);

	 Integral_syst_Pt[g][ibin]->SetMarkerSize(2);

	 Integral_syst_Pt[g][ibin]->SetLineColor(kBlack);
	 Integral_syst_Pt[g][ibin]->SetMarkerColor(kBlack);
  
       }
     }

     TCanvas *c_integral = new TCanvas("Integral_Canvas","",10,10,2000,1000);
     c_integral->Divide(4,2,0,0);

     for(int ibin = 0; ibin < nCBins; ibin++){
       c_integral->cd(4-ibin);


       Integral_Pt[0][ibin]->GetXaxis()->SetRangeUser(0.01,19.5);
  
       Integral_Pt[0][ibin]->SetMaximum(14.);
       Integral_Pt[0][ibin]->SetMinimum(-2.);

       Integral_Pt[0][ibin]->GetXaxis()->SetLabelSize(0.06);
       Integral_Pt[0][ibin]->GetXaxis()->SetTitleSize(0.06);
       Integral_Pt[0][ibin]->GetXaxis()->SetTitle("p_{T}^{assoc}");

       if(ibin==3){
	 Integral_Pt[0][ibin]->GetYaxis()->SetLabelSize(0.06);
	 Integral_Pt[0][ibin]->GetYaxis()->SetTitleSize(0.06);
	 Integral_Pt[0][ibin]->GetYaxis()->SetTitle("Tracks per jet");
       }else{
	 Integral_Pt[0][ibin]->GetYaxis()->SetLabelSize(0.);
       }
   
       Integral_Pt[0][ibin]->Draw();
       Integral_syst_Pt[0][ibin]->Draw("same e2");
       Integral_syst_Pt[1][0]->Draw("same e2");
       Integral_Pt[1][0]->Draw("same");
      
       Integral_Pt[0][ibin]->Draw("same");
       Integral_Pt[1][0]->Draw("same");
   
       if(ibin==3){
	 TLegend *int_legend = new TLegend(0.18,0.6,0.9,0.85);
	 int_legend->AddEntry( Integral_syst_Pt[0][ibin],"PbPb integral |#Deltar| < 1.0");
	 int_legend->AddEntry( Integral_syst_Pt[1][0],"pp integral |#Deltar| < 1.0");
	 int_legend->SetLineColor(kWhite);
	 int_legend->SetTextSize(0.06);
	 int_legend->Draw();

	 labels = new TPaveText(0.18,0.85,0.45,0.95,"NDC");
    
	 labels->SetName("labels");
	 labels->SetFillColor(0);
	 labels->SetLineColor(0);
	 labels->SetTextAlign(11);
	 labels->AddText(CBin_labels[ibin]);
	 labels->SetTextSize(0.06);
	 labels->Draw("same");
  
       }else{
	 labels = new TPaveText(0.05,0.85,0.45,0.95,"NDC");
    
	 labels->SetName("labels");
	 labels->SetFillColor(0);
	 labels->SetLineColor(0);
	 labels->SetTextAlign(11);
	 labels->AddText(CBin_labels[ibin]);
	 labels->SetTextSize(0.06);
	 labels->Draw("same");
  

       }

       TLine *int_zero = new TLine(0.,0.,20.,0.);
       int_zero->SetLineStyle(2);
       int_zero->SetLineColor(kBlack);
       int_zero->Draw();

       TPave *cover_x = new TPave(0.9,0.1,1.0,0.2);
       cover_x->SetFillColor(kBlack);
       cover_x->Draw();


       c_integral->cd(8-ibin);

       Integral_diff_Pt[0][ibin] = (TH1D*)  Integral_Pt[0][ibin]->Clone((TString)("Integral_Diff_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));
       Integral_diff_Pt[0][ibin]->Add(Integral_Pt[1][0],-1.);
  
       Integral_diff_syst_Pt[0][ibin] = (TH1D*)  Integral_syst_Pt[0][ibin]->Clone((TString)("Integral_Diff_Syst_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));
       Integral_diff_syst_Pt[0][ibin]->Add(Integral_syst_Pt[1][0],-1.);
       
    
       Integral_diff_Pt[0][ibin]->SetMaximum(8.5);
       Integral_diff_Pt[0][ibin]->SetMinimum(-2.);
       Integral_diff_Pt[0][ibin]->Draw();
       Integral_diff_syst_Pt[0][ibin]->Draw("same e2");
     Integral_diff_Pt[0][ibin]->Draw("same");

       if(ibin==3){
	 labels = new TPaveText(0.18,0.85,0.45,0.95,"NDC");
	 labels->SetTextSize(0.055);
       }else{
	 labels = new TPaveText(0.05,0.85,0.45,0.95,"NDC");
	 labels->SetTextSize(0.06);
       }
       labels->SetName("labels");
       labels->SetFillColor(0);
       labels->SetLineColor(0);
       labels->SetTextAlign(11);
       labels->AddText((TString)("PbPb (" + CBin_labels[ibin]+") minus pp"));
   
       labels->Draw("same");
  
       int_zero->Draw();

       c_integral->cd(0);


       type_tex = new TLatex(0.02,0.92,"Total Particle Yields");
       type_tex->SetTextSize(0.035);
       type_tex->SetLineColor(kWhite);
       type_tex->SetNDC();
       type_tex->Draw();
       
       luminosity_tex_pp->Draw();
       luminosity_tex_PbPb->Draw();
       jet_reco_tex->Draw();



       Integral_Pt[0][ibin]->Write();
       Integral_syst_Pt[0][ibin]->Write();
       Integral_diff_Pt[0][ibin]->Write();
       Integral_diff_syst_Pt[0][ibin]->Write();

       if(ibin==0){

	 Integral_Pt[1][0]->Write();
	 Integral_syst_Pt[1][0]->Write();
       }



     }

     c_integral->SaveAs("Integral_dR.png");
     c_integral->SaveAs("Integral_dR.pdf");


   }
  
  return 0;
}// main loop
