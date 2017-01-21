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

Int_t jet_shapes_closures3(bool is_number=0, bool use_highpT_bin = kTRUE, int type_code = 0, int do_closures = 0, int ref_is_reco = 0){

  TFile *fMC[6], *f_ref[7];

  TCanvas *c_jetshape;
   
  TString jetetacut, etalabel,centlabel,pTlabel,ajlabel;
  float eta_ymax;

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins, limR;
  float deta, dphi, r, bc, bg_err,temp1,temp2, rbin, temperr, err, width_temp_x, width_temp_y, width_temp, norm_temp, zerobin, temp, norm_tot, err_temp, cont, ratio_axis_min, ratio_axis_max, dr_axis_max, jetshape_axis_min, jetshape_axis_max;
  
 
  float RBins[20] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.,1.2,1.4,1.6,1.8,2.0};
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 10;

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

   float TrkPtBins[nTrkPtBins] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 300};
   TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300","ALL" };
   TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
 
  float mean_pts[nTrkPtBins] = {0.844,1.35,2.35,3.37,5.07,9.72,13.8,17.9,22.};

  enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RightGen, SpilledUnderGen, UnmatchedGen, RightReco, SpilledReco, UnmatchedReco, RecoGenSube0,RecoGenNoSube0,GenGenSube0,GenGenNoSube0,MatchedRecoGenSube0,MatchedRecoGenNoSube0,SwappedRecoGenSube0,SwappedRecoGenNoSube0, UnMatchedRecoGenSube0,UnMatchedRecoGenNoSube0,n_data_mc_types};


  TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RightGenJet_GenTrack","SpilledUnderJet_GenTrack","UnmatchedGenJet_GenTrack","RightRecoJet_GenTrack","SpilledReco_GenTrack","UnmatchedReco_GenTrack","RecoJet_GenTrack_Sube0","RecoJet_GenTrack_NoSube0","GenJet_GenTrack_Sube0","GenJet_GenTrack_NoSube0","MatchedRecoJet_GenTrack_Sube0","MatchedRecoJet_GenTrack_NoSube0","SwappedRecoJet_GenTrack_Sube0","SwappedRecoJet_GenTrack_NoSube0","UnmatchedRecoJet_GenTrack_Sube0","UnmatchedRecoJet_GenTrack_NoSube0",};
  
  TString data_mc_type_labels[9] = {"Data","RecoJet, RecoTrack","RecoJet, GenTrack","GenJet, RecoTrack", "GenJet, GenTrack","RecoReco over GenGen","RecoGen over GenGen", "GenReco over GenGen"};


  TString do_closures_string = "";


  TFile *f_in_gen;

  if(do_closures==1){
    if(type_code==4)cout<<"THIS WILL BE ZERO BY DEFINITION!!"<<endl;
    if(ref_is_reco){

      do_closures_string = "_over_RecoGen";
      f_in_gen = new TFile((TString)("Jet_Shapes_ClosuresRecoJet_GenTrack.root"));
   
    }else{
      do_closures_string = "_over_GenGen";
      f_in_gen = new TFile((TString)("Jet_Shapes_ClosuresGenJet_GenTrack.root"));
    }
  }


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

 
  double temp_cont, nextr, nextl, cos_weight,me00_range, mc_error;
  
  TString stem, datalabel,me00_range_string,stem_mc;

 
 
  float norm;

  TFile *fin_PbPb = new TFile((TString)("../me_correct/Hydjet_"+data_mc_type_strs[type_code]+"_Inclusive_Correlations.root"));
  TFile *fin_pp = new TFile((TString)("../me_correct/Pythia_"+data_mc_type_strs[type_code]+"_Inclusive_Correlations.root"));
  TFile *fout  = new TFile((TString)("Jet_Shapes_Closures"+data_mc_type_strs[type_code]+do_closures_string+".root"),"RECREATE"); 

  TFile *f_jff_pyth_gen = new TFile("../me_correct/Pythia_GenJet_GenTrack_Inclusive_Correlations.root", "READ");
  TFile *f_jff_pyth_reco = new TFile("../me_correct/Pythia_RecoJet_GenTrack_Inclusive_Correlations.root", "READ");
  TFile *f_jff_pyth_reco_reco = new TFile("../me_correct/Pythia_RecoJet_RecoTrack_Inclusive_Correlations.root", "READ");

  TFile *f_jff_hyd_gen2 = new TFile("../me_correct/HydJet_GenJet_GenTrack_Inclusive_Correlations.root", "READ");
  TFile *f_jff_hyd_gen = new TFile("../me_correct/HydJet_GenJet_GenTrack_Sube0_Inclusive_Correlations.root", "READ");
  TFile *f_jff_hyd_reco = new TFile("../me_correct/HydJet_RecoJet_GenTrack_Sube0_Inclusive_Correlations.root", "READ");
  
TFile *f_jff_hyd_reco_reco = new TFile("../me_correct/HydJet_RecoJet_RecoTrack_Inclusive_Correlations.root", "READ");



  TFile *f_spillover = new TFile("../spill_over/Inclusive_Hydjet_SpillOvers.root");
  if(!is_number) f_spillover = new TFile("../spill_over/Inclusive_Hydjet_SpillOvers_pTweighted.root");


 
  TH1D *Integral_Pt[12][5];
  TH1D *Integral_diff_Pt[12][5];


  double integral[12][nCBins][nTrkPtBins];
  double integral_err[12][nCBins][nTrkPtBins];
  

  TPaveText *labels;

  TF2* ClosureFit = new TF2("ClosureFit", "[0]/2/TMath::Pi()/[1]/[2]*TMath::Exp(-1.*(x*x/[1]/[1]/2))*TMath::Exp(-1.*(y*y/[2]/[2]/2))",-5.,5.,-TMath::Pi()/2.,3*TMath::Pi()/2.);

  TCanvas *c_dummy = new TCanvas("dummy");


  //-----------------------
  // Start getting histos
  //-----------------------
  for(int g=0; g<2; g++){

  
    switch(g){
    case 0:
      if(is_number) stem = "Yield_BkgSub_";
      else  stem = "Yield_BkgSub_pTweighted";
      datalabel = "Inclusive";
      break;
    case 1:
      if(is_number) stem = "Yield_BkgSub_";
      else stem = "Yield_BkgSub_pTweighted";
      //stem = "Yield_BkgSub_";
      datalabel = "Inclusive";
      break;
    }
  
    /*
  switch(g){
    case 0:
      if(is_number) stem = "Yield_Hydjet_";
      else  stem = "Yield_pTweighted_Hydjet_";
      datalabel = "Inclusive";
      break;
    case 1:
      if(is_number) stem = "Yield_Pythia_";
      else stem = "Yield_pTweighted_Pythia_";
      //stem = "Yield_BkgSub_";
      datalabel = "Inclusive";
      break;
    }
    */

    for (int ibin=0;ibin<nCBins;ibin++){

      if(g==1&&ibin > 0) continue;

      for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
    
	for (int ibin3=0;ibin3<nTrkPtBins-1;ibin3++){

	  if(g==0){

	    result[g][ibin][ibin2][ibin3] = (TH2D*)fin_PbPb->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + datalabel+"_PbPb_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  }else{
	    result[g][ibin][ibin2][ibin3] = (TH2D*)fin_pp->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + datalabel+"_PbPb_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  }
	  if((type_code==1||type_code==2)){
	    if(g==0){
	
	      if(ibin3>3){
		resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_hyd_gen2->Get((TString)(stem+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"GenGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

		resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_hyd_reco_reco->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    

	      }else{
	    

		resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_hyd_gen->Get((TString)(stem+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"GenGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

		resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_hyd_reco->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	 
	 
		TString	spillover_name =(TString)("Eta_SpillOver_Points_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]);
		spill_over_dEta[g][ibin][ibin3] = (TH1D*)f_spillover->Get(spillover_name)->Clone(spillover_name);
	    
	 
		spillover_name.ReplaceAll("Eta","Phi");
	  
		spill_over_dPhi[g][ibin][ibin3] = (TH1D*)f_spillover->Get(spillover_name)->Clone(spillover_name);
	  

		gaus_eta[g][ibin][ibin3] = (TF1*)f_spillover->Get((TString)("Eta_SpillOver_Fit_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	  
		gaus_phi[g][ibin][ibin3] = (TF1*)f_spillover->Get((TString)("Phi_SpillOver_Fit_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	   
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
	    
	      
		/*
		  if(!is_number){ //pT-weighted spill-overs aren't well controlled.  weight by hand. 
		

		  resultMC[g][ibin][0][ibin3]->Scale(mean_pts[ibin3]);
	

		  }
		*/

	      }
	 
	  
	    }else{
	    

	      resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_pyth_gen->Get((TString)(stem+ CBin_strs[0] + "_" + CBin_strs[1] + "_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"pp_GenGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    
	      if(ibin3>3)	    resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_pyth_reco_reco->Get((TString)(stem + CBin_strs[0] + "_" + CBin_strs[1] + "_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"pp_RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	      else   resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_pyth_reco->Get((TString)(stem + CBin_strs[0] + "_" + CBin_strs[1] + "_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"pp_RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	    }
	
		  	    
	    resultMC_reco[g][ibin][0][ibin3]->Add(   resultMC_gen[g][ibin][0][ibin3],-1.);

	    
	    for(int k = 0; k< resultMC_reco[g][ibin][0][ibin3]->GetNbinsX()+1; k++){
	      for(int m = 0; m< resultMC_reco[g][ibin][0][ibin3]->GetNbinsY()+1; m++){


		deta = resultMC_reco[g][ibin][0][ibin3]->GetXaxis()->GetBinCenter(k);
		dphi = resultMC_reco[g][ibin][0][ibin3]->GetYaxis()->GetBinCenter(m);
	   
		r = TMath::Sqrt(deta*deta+dphi*dphi);
	      
	
		if(r>1.0){
	      

		  resultMC_reco[g][ibin][0][ibin3]->SetBinContent(k,m,0.);
		  resultMC_reco[g][ibin][0][ibin3]->SetBinError(k,m,0.);
		}
	      }
	    }
	
	    if(g==0&&ibin3<4){
	      if(ibin3==0)resultMC[g][ibin][0][ibin3]->Scale(0.3);
	    	
	      resultMC_reco[g][ibin][ibin2][ibin3]->Add(resultMC[g][ibin][0][ibin3]);
	
	    }


	    result[g][ibin][ibin2][ibin3]->Add(resultMC_reco[g][ibin][0][ibin3],-1.);
	  

	  }

	  if(result[g][ibin][0][ibin3]->GetEntries()<1e-10){

	    cout<<"Corrections and/or references completely cancel correlations"<<endl;
	    //      return 0;
	  }

	  if(is_number){
	    if(ibin3==0) result[g][ibin][ibin2][ibin3]->Scale(1./0.3);
	    if(ibin3> 3&& ibin3<nTrkPtBins-1) result[g][ibin][ibin2][ibin3]->Scale(1./4);
	    if(type_code<3){
	    if(ibin3==0) resultMC_reco[g][ibin][ibin2][ibin3]->Scale(1./0.3);
	    if(ibin3> 3&& ibin3<nTrkPtBins-1) resultMC_reco[g][ibin][ibin2][ibin3]->Scale(1./4);
	    }
	  }

	

	  if(ibin3==0){
	    result[g][ibin][ibin2][9] = (TH2D*) result[g][ibin][ibin2][ibin3]->Clone(((TString) ("Summed_result_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])));

	  }else if(!is_number&&ibin3>0&&ibin3<nTrkPtBins-1){
	    
	    result[g][ibin][ibin2][9]->Add(result[g][ibin][ibin2][ibin3]);
	 
	  }else if(is_number&&ibin3>0&&ibin3<nTrkPtBins-2){
	    
	    result[g][ibin][ibin2][9]->Add(result[g][ibin][ibin2][ibin3]);
	 
	  }

	}


	for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

	  if(g==0){	  
	  JetShape2[g][ibin][ibin3] = new TH1D((TString)("JetShape2_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  

	  JetShape2_geo[g][ibin][ibin3] = new TH1D((TString)("JetShape_Geometry_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  
	  }else{
	    JetShape2[g][ibin][ibin3] = new TH1D((TString)("JetShape2_"+stem+ datalabel+"_pp_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  

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
	      	      
	  	  
	    }
	  }

	  JetShape[g][ibin][ibin3] = (TH1D*)JetShape2[g][ibin][ibin3]->Clone((TString)("JetShape_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	
	  for(int k = 1; k<JetShape2[g][ibin][ibin3]->GetNbinsX()+1; k++){
	    
	    float temp = JetShape2[g][ibin][ibin3]->GetBinContent(k)/JetShape2[g][ibin][ibin3]->GetBinWidth(k);
	    float err = JetShape2[g][ibin][ibin3]->GetBinError(k)/JetShape2[g][ibin][ibin3]->GetBinWidth(k);
	    
	    JetShape2[g][ibin][ibin3]->SetBinContent(k,temp);
	    JetShape2[g][ibin][ibin3]->SetBinError(k,err);

	    JetShape[g][ibin][ibin3]->SetBinContent(k,temp);
	    JetShape[g][ibin][ibin3]->SetBinError(k,err);

	  }

	  if(is_number){

	    int end_bin = JetShape2[g][ibin][ibin3]->FindBin(.999);
	    integral[g][ibin][ibin3] =  JetShape2[g][ibin][ibin3]->IntegralAndError(1,end_bin,integral_err[g][ibin][ibin3],"width");

	    // cout<<"norm is "<<integral[g][ibin][ibin3]<<" "<<g<<" "<<ibin<<" "<<ibin3<<endl;
	  
	  }

	  if(!is_number&&ibin3==9){



	    int end_bin = JetShape2[g][ibin][ibin3]->FindBin(.999);

	    norm =  JetShape2[g][ibin][ibin3]->Integral(1,end_bin,"width");

	    norm = 1.;

	     cout<<"norm is "<<norm<<" "<<g<<" "<<ibin<<" "<<ibin3<<endl;

	    //	    if(!use_highpT_bin) norm = JetShape[g][ibin][ibin3]->GetBinContent(1)/JetShape_ref[g][ibin][ibin3]->GetBinContent(1);	 
	    for(int k = 0; k<10; k++){
	      JetShape[g][ibin][k]->Scale(1./norm);  
	      JetShape2[g][ibin][k]->Scale(1./norm);  
	    }
	  }
	  
	  if(ibin<3){
	    JetShape2[g][ibin][ibin3]->GetYaxis()->SetLabelSize(0.);
	    JetShape2[g][ibin][ibin3]->GetYaxis()->SetTitleSize(0.);

	  }

	  /*
	  JetShapeMC[g][ibin][ibin3]->Divide(JetShape2_geo[g][ibin][ibin3]);

	  JetShapeMC[g][ibin][ibin3]->Scale(1./norm);
	  */

	  fout->cd();
	  
	  JetShape2[g][ibin][ibin3]->SetAxisRange(0.,0.99,"x");
	  JetShape[g][ibin][ibin3]->SetAxisRange(0.,0.99,"x");



	
	}

	for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

	  if(do_closures==1){
	 
	      if(g==0)	  JetShape2[g+2][ibin][ibin3] = (TH1D*)f_in_gen->Get((TString)("JetShape2_Yield_BkgSub_Inclusive_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString)("Yield_BkgSub_Ref_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_" + TrkPtBin_strs[ibin] + "_" + TrkPtBin_strs[ibin3+1]));
	    
	      else  JetShape2[g+2][ibin][ibin3] = (TH1D*)f_in_gen->Get((TString)("JetShape2_Yield_BkgSub_Inclusive_pp_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString)("Yield_BkgSub_Ref_pp_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	  	     
	      for(int k = 0; k<JetShape2[g][ibin][ibin3]->GetNbinsX()+1; k++){
		if(TMath::Abs( (JetShape2[g+2][ibin][ibin3]->GetBinContent(k)- JetShape2[g+2][ibin][ibin3]->GetBinError(k)))<1.5*JetShape2[g+2][ibin][ibin3]->GetBinError(k)){
		  JetShape2[g][ibin][ibin3]->SetBinContent(k,0.);
		  JetShape2[g][ibin][ibin3]->SetBinError(k,0.);
		}
	      }
		

	      cout<<"Integral: "<<g<<" "<<ibin<<" "<<ibin3<<" "<<JetShape2[g][ibin][ibin3]->Integral()<<" "<<JetShape2[g+2][ibin][ibin3]->Integral()<<endl;
	      JetShape2[g][ibin][ibin3]->Divide( JetShape2[g+2][ibin][ibin3]);
	      JetShape[g][ibin][ibin3]->Divide( JetShape2[g+2][ibin][ibin3]);

	 
		  
	  }
	  JetShape2[g][ibin][ibin3]->Write();
	}

      }

    }
  }
 

  cout<<"ready to draw"<<endl;
  TString canvas_name = "JetShape_Comparison_";

  if(is_number) canvas_name = "ParticleYield_dR_Comparison_";
  
  canvas_name+=data_mc_type_strs[type_code];
    
  canvas_name+=do_closures_string;
   
  c_jetshape= new TCanvas(canvas_name," ",10,10,1500,3600);
  c_jetshape->Divide(4,9,0.,0.);

 
  for(int k = 0; k<nTrkPtBins; k++){
   
    for(int ibin = 0; ibin<nCBins; ibin++){
      
      c_jetshape->cd(4*(k+1)-ibin);

      JetShape2[0][ibin][k]->GetXaxis()->SetTitle("Radius(r)");
      JetShape2[0][ibin][k]->GetXaxis()->CenterTitle(true);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelFont(42);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelOffset(0.02);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelSize(0.08);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleSize(0.08);
      JetShape2[0][ibin][k]->GetXaxis()->SetTickLength(0.025);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleOffset(1.3);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleFont(42);
     
      if(is_number)    JetShape2[0][ibin][k]->GetYaxis()->SetTitle("Y(#Delta r)");
      else  JetShape2[0][ibin][k]->GetYaxis()->SetTitle("#rho(#Deltar)");
      JetShape2[0][ibin][k]->GetYaxis()->CenterTitle(true);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelFont(42);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelOffset(0.004);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.075);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitleSize(0.09);
      JetShape2[0][ibin][k]->GetYaxis()->SetTickLength(0.025);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitleOffset(1.);
     

      if(!do_closures){
	JetShape2[0][ibin][k]->SetMinimum(-1.);
	JetShape2[0][ibin][k]->SetMaximum(1.);
      }else if(is_number){
		JetShape2[0][ibin][k]->SetMinimum(.5);
		JetShape2[0][ibin][k]->SetMaximum(1.5);
      }else{

	JetShape2[0][ibin][k]->SetMinimum(.5);
	JetShape2[0][ibin][k]->SetMaximum(1.5);

      }
      // JetShape2[0][ibin][k]->SetMinimum(0.5);
      //JetShape2[0][ibin][k]->SetMaximum(1500.1);

      JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.07);
    

      //    gPad->SetLogy();

      JetShape2[0][ibin][k]->SetMarkerSize(2);
      JetShape2[0][ibin][k]->SetLineColor(kBlack);
      JetShape2[0][ibin][k]->SetMarkerColor(kBlack);
      JetShape2[0][ibin][k]->SetMarkerStyle(20);


    
     
      JetShape2[1][0][k]->SetMarkerSize(2);
      JetShape2[1][0][k]->SetLineColor(kBlack);
      JetShape2[1][0][k]->SetMarkerColor(kBlack);
      JetShape2[1][0][k]->SetMarkerStyle(24);
      

      if(ibin<3){
	JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.);
	JetShape2[0][ibin][k]->GetYaxis()->SetTitleSize(0.);
      }

      //   gPad->SetLogy();

      JetShape2[0][ibin][k]->GetXaxis()->SetRangeUser(0.,.99);
      //    JetShape2[0][ibin][k]->Add( JetShape2[1][0][k],-1.);
      JetShape2[0][ibin][k]->Draw("p");

      
  
      JetShape2[1][0][k]->Draw("same p");


      TLine *l = new TLine(0.,1.,2.,1.);
      l->SetLineStyle(2);
      l->SetLineColor(kBlack);
      l->Draw();
  


      TLine *l_up = new TLine(0.,1.05,2.,1.05);
      l_up->SetLineStyle(2);
      l_up->SetLineColor(kBlack);
      l_up->Draw();
 

     TLine *l_down = new TLine(0.,.95,2.,.95);
      l_down->SetLineStyle(2);
      l_down->SetLineColor(kBlack);
      l_down->Draw();
  
      /*
      if(((k>0&&k<5)||k==6)&&!is_subleading){
	JetShape_ref[0][ibin][k]->Draw("same");  



	JetShape_ref[1][0][k]->SetMarkerStyle(20);

	JetShape_ref[1][0][k]->Draw("same");
      }
      */

      if(ibin==3&&k==0){
	TLegend *legend = new TLegend(0.2,0.65,0.5,0.75);
	legend->AddEntry(JetShape2[0][ibin][k],"PbPb");
	legend->AddEntry(JetShape2[1][0][k],"pp");
	legend->SetLineColor(kWhite);
	legend->SetFillColor(kWhite);
	legend->SetTextSize(0.065);
	legend->Draw();
      }
      TLatex  *pttex = new TLatex(0.05,0.81,TrkPtBin_labels[k]);
      pttex->SetNDC();
      pttex->SetTextSize(0.075);
     
    
      TLatex  *centtex ;
      centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
      centtex->SetNDC();
      centtex->SetTextSize(0.075);
  
      if(ibin==3){  
	centtex = new TLatex(0.2,0.9,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->SetTextSize(0.07);


	pttex = new TLatex(0.2,0.81,TrkPtBin_labels[k]);
	pttex->SetNDC();
	pttex->SetTextSize(0.075);

      }

      pttex->Draw();      
      centtex->Draw();

     
     

      c_dummy->cd();

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
      /*
      JetShape_diff2[0][ibin][k]->Draw();
      
      if(((k>0&&k<5)||k==6)&&!is_subleading){
	diff_name ="JetShape_ref_diff_"; diff_name+=CBin_strs[ibin]; diff_name+= "_"; diff_name += CBin_strs[ibin+1]; diff_name+= k;
      
     
	JetShape_ref_diff[0][ibin][k] = (TH1D*)JetShape_ref[0][ibin][k]->Clone(diff_name);

	JetShape_ref_diff[0][ibin][k]->Add( JetShape_ref[1][0][k],-1. );

	JetShape_ref_diff[0][ibin][k]->SetMarkerStyle(20);
	JetShape_ref_diff[0][ibin][k]->Draw("same");
      }

      */

      //    cout<<"about to make ratio "<<ibin<<" "<<k<<endl;

      TString ratio_name ="JetShape_ratio_"; ratio_name+=CBin_strs[ibin]; ratio_name+= "_"; ratio_name += CBin_strs[ibin+1]; ratio_name+= k;

      TString ratio2_name ="JetShape_ratio2_"; ratio2_name+=CBin_strs[ibin]; ratio2_name+= "_"; ratio2_name += CBin_strs[ibin+1]; ratio2_name+= k;
      
     
      JetShape_ratio[0][ibin][k] = (TH1D*)JetShape[0][ibin][k]->Clone(ratio_name);
      JetShape_ratio2[0][ibin][k] = (TH1D*)JetShape2[0][ibin][k]->Clone(ratio2_name);

      if(!is_number){ 
	JetShape_ratio[0][ibin][k]->Divide( JetShape[1][0][k]);
	JetShape_ratio2[0][ibin][k]->Divide( JetShape2[1][0][k] );
      
      }else{
	JetShape_ratio[0][ibin][k]->Divide( JetShape[1][0][k]);
	JetShape_ratio2[0][ibin][k]->Divide( JetShape2[1][0][k]);
      }

      JetShape_ratio2[0][ibin][k]->SetMinimum(-1.);
      JetShape_ratio2[0][ibin][k]->SetMaximum(9.);
     

      JetShape_ratio2[0][ibin][k]->SetMarkerStyle(20);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitle("#rho(r)_{PbPb} - #rho(r)_{pp}");
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

      //    JetShape_ratio2[0][ibin][k]->Draw();
      /*
      if(((k>0&&k<5)||k==6)&&!is_subleading){
	ratio_name ="JetShape_ref_ratio_"; ratio_name+=CBin_strs[ibin]; ratio_name+= "_"; ratio_name += CBin_strs[ibin+1]; ratio_name+= k;
      
     
	JetShape_ref_ratio[0][ibin][k] = (TH1D*)JetShape_ref[0][ibin][k]->Clone(ratio_name);

	JetShape_ref_ratio[0][ibin][k]->Divide( JetShape_ref[1][0][k] );

	JetShape_ref_ratio[0][ibin][k]->SetMarkerStyle(20);
	JetShape_ref_ratio[0][ibin][k]->Draw("same");
      }
      */
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelSize(0.);

      //   cout<<"and done"<<endl;

    }


  }
  cout<<"here"<<endl;
  
  canvas_name+=".pdf";
    
  c_jetshape->SaveAs(canvas_name);
    
  canvas_name.ReplaceAll(".pdf",".png");
  c_jetshape->SaveAs(canvas_name);
  
  cout<<"saved individual canvases"<<endl;

  for(int ibin = 0; ibin<nCBins; ibin++){
    
    for(int k = 0; k<nTrkPtBins; k++){

      
 
      JetShape_noerr_up[0][ibin][k] = new TH1D((TString)("JetShape_PbPb_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
      JetShape_noerr_down[0][ibin][k] = new TH1D((TString)("JetShape_PbPb_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
   
      if(ibin==0){
	JetShape_noerr_up[1][0][k] = new TH1D((TString)("JetShape_pp_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
	JetShape_noerr_down[1][0][k] = new TH1D((TString)("JetShape_pp_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
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

	  bc = JetShape[1][0][k]->GetBinContent(l);
      
	  if(bc>0){ 
	    JetShape_noerr_up[1][0][k]->SetBinContent(l,bc);	
	    JetShape_noerr_down[1][0][k]->SetBinContent(l,0.);	

	  }else{
	    JetShape_noerr_down[1][0][k]->SetBinContent(l,bc);	
	    JetShape_noerr_up[1][0][k]->SetBinContent(l,0.);	
	  }
	}
      
      }


      //   cout<<"Have set for "<<ibin<<" "<<k<<" "<<JetShape_noerr_up[0][ibin][k]->GetBinContent(1)<<" "<<JetShape_noerr_up[0][ibin][k]->GetBinCenter(1)<<endl;

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

    dr_axis_max = 0.99;


  
    if(!is_number){

      //jetshape_axis_max = 23.5;
      //    jetshape_axis_min = .005;

      jetshape_axis_max = 2000.;
      jetshape_axis_min = .1;
      ratio_axis_min = 0.;
      ratio_axis_max = 6.5;

    }else{

      jetshape_axis_max = 34.;
      jetshape_axis_min = -1.;
      ratio_axis_min = -10.;
      ratio_axis_max = 25.;
  
    }

    if(do_closures==1){

      jetshape_axis_max = 25.;
      jetshape_axis_min = -10.;
   
      ratio_axis_max = 25.;
      ratio_axis_min = -10.;
 
      if(!is_number){
	jetshape_axis_max = 25.;
	jetshape_axis_min = -5.;
    
      }
 
    }



    JetShape_Stack_Up[0][ibin]->SetMaximum(jetshape_axis_max);
    JetShape_Stack_Up[1][0]->SetMaximum(jetshape_axis_max);
    JetShape_Stack_Up[0][ibin]->SetMinimum(jetshape_axis_min);
    JetShape_Stack_Up[1][0]->SetMinimum(jetshape_axis_min);
    JetShape_ratio[0][ibin][9]->SetMinimum(ratio_axis_min);
    JetShape_ratio[0][ibin][9]->SetMaximum(ratio_axis_max);
 
  


  } //ibin

  cout<<"Starting stacked plot drawing"<<endl;


  

  TCanvas *PAS_plot = new TCanvas("JetShape_ForPAS","",1200,800);
  PAS_plot->Divide(3,2,0.,0.);


  PAS_plot->cd(4);
  TLegend *legend = new TLegend(0.1,0.05,0.9,.9);
  legend->AddEntry(JetShape_noerr_up[0][0][0],"0.5 < p_{T}^{assoc.}< 1 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][1],"1 < p_{T}^{assoc.}< 2 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][2],"2 < p_{T}^{assoc.}< 3 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][3],"3 < p_{T}^{assoc.}< 4 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][4],"4 < p_{T}^{assoc.}< 8 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][5],"8 < p_{T}^{assoc.}< 12 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][6],"12 < p_{T}^{assoc.}< 16 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][7],"16 < p_{T}^{assoc.}< 20 GeV","f");
 

  if(!use_highpT_bin||is_number){
    //   legend->AddEntry(JetShape_graph[0][0][9],"Total 0.5 < p_{T}^{assoc.}<  GeV","lpfe");
  }else{ 
    legend->AddEntry(JetShape_noerr_up[0][0][8],"p_{T}^{assoc.}> 20 GeV","f");
    //  legend->AddEntry(JetShape_graph[0][0][9],"Total p_{T}^{assoc.}> 0.5 GeV","lpfe");
  }
  legend->AddEntry(JetShape_noerr_up[0][0][4],"|#eta_{track}|< 2.4","");
  legend->SetTextSize(0.055);
  legend->SetLineColor(kWhite);
  legend->SetFillColor(kWhite);
  legend->Draw();

  
    

  PAS_plot->cd(0);
  
  TLatex *type_tex;
  if(is_number)type_tex = new TLatex(0.02,0.96,"Particle Yield by #Deltar");
  else  type_tex = new TLatex(0.02,0.96,"Inclusive Jet Shape");
 
  type_tex->SetTextSize(0.035);
  type_tex->SetLineColor(kWhite);
  type_tex->SetNDC();
  type_tex->Draw();
   
  TLatex   *mc_type_tex;

  if(do_closures){
    mc_type_tex = new TLatex(0.02,0.92,data_mc_type_labels[type_code+4]);
    mc_type_tex->SetTextFont(43);
    mc_type_tex->SetTextSizePixels(25);
    mc_type_tex->SetLineColor(kWhite);
    mc_type_tex->SetNDC();
    mc_type_tex->Draw();

  }else{

    mc_type_tex = new TLatex(0.02,0.92,data_mc_type_labels[type_code]);
    mc_type_tex->SetTextFont(43);
    mc_type_tex->SetTextSizePixels(25);
    mc_type_tex->SetLineColor(kWhite);
    mc_type_tex->SetNDC();
    mc_type_tex->Draw();

  }
 
    TLatex   *luminosity_tex_PbPb = new TLatex(0.3,0.92,"PYTHIA+HYDJET");
    luminosity_tex_PbPb->SetTextFont(43);
    luminosity_tex_PbPb->SetTextSizePixels(25);
    luminosity_tex_PbPb->SetLineColor(kWhite);
    luminosity_tex_PbPb->SetNDC();
    luminosity_tex_PbPb->Draw();
 
    TLatex   *jet_reco_tex = new TLatex(0.605,0.96,"ak4CaloJets, |#eta_{jet}| < 1.6");
    jet_reco_tex->SetTextFont(43);
    jet_reco_tex->SetTextSizePixels(25);
    jet_reco_tex->SetLineColor(kWhite);
    jet_reco_tex->SetNDC();
    jet_reco_tex->Draw();

    TLatex   *jet_cut_tex = new TLatex(0.605,0.92,"p_{T}> 120");
    jet_cut_tex->SetTextFont(43);
    jet_cut_tex->SetTextSizePixels(25);
    jet_cut_tex->SetLineColor(kWhite);
    jet_cut_tex->SetNDC();
    jet_cut_tex->Draw();
   
  PAS_plot->cd(1);
  


 
  JetShape_Stack_Up[1][0]->Draw();
  JetShape_Stack_Up[1][0]->GetXaxis()->SetRangeUser(0.,dr_axis_max);
  JetShape_Stack_Up[1][0]->GetXaxis()->SetNdivisions(505);

   if(use_highpT_bin&&!is_number&&!do_closures)  gPad->SetLogy();
  JetShape_Stack_Up[1][0]->GetYaxis()->SetLabelSize(0.06);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleSize(0.06);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleOffset(1.1);
  if(is_number) JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("Y = #frac{1}{N_{jet}} #frac{dN}{d#Deltar}");
  else JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("#rho(#Deltar)");
  //JetShape_Stack_Up[1][0]->GetYaxis()->CenterTitle();
  JetShape_Stack_Down[1][0]->Draw("same");
  JetShape_Stack_Up[1][0]->Draw("same");

 

  // JetShape_graph[1][0][9]->Draw("same e2 P");
  JetShape2[1][0][9]->Draw("same");
  // if(!is_subleading)JetShape_ref[1][0][9]->Draw("same");

  TLatex  *label_pp = new TLatex(0.2,0.9,"PYTHIA");
  label_pp->SetTextSize(0.09);
  label_pp->SetLineColor(kWhite);
  label_pp->SetNDC();
  label_pp->Draw();

  TLine *l_dr = new TLine(0.,0.,TMath::Pi()/2.,0.);
  l_dr->SetLineStyle(2);
  l_dr->Draw();


  TLatex *cms_tex_dphi = new TLatex(0.25,0.8,"CMS");
  cms_tex_dphi->SetTextFont(63);
  cms_tex_dphi->SetTextSizePixels(30);
  cms_tex_dphi->SetLineColor(kWhite);
  cms_tex_dphi->SetNDC();
  cms_tex_dphi->Draw(); 


TLatex *prelim_tex_dphi = new TLatex(0.45,0.8,"Preliminary");
  prelim_tex_dphi->SetTextFont(53);
  prelim_tex_dphi->SetTextSizePixels(30);
  prelim_tex_dphi->SetLineColor(kWhite);
  prelim_tex_dphi->SetNDC();
  prelim_tex_dphi->Draw(); 


  gPad->RedrawAxis();
  
  PAS_plot->cd(2);


  
  JetShape_Stack_Up[0][3]->Draw();

 
  JetShape_Stack_Up[0][3]->GetXaxis()->SetRangeUser(0.,0.99);
  JetShape_Stack_Up[0][3]->GetXaxis()->SetNdivisions(505);

   if(use_highpT_bin&&!is_number&&!do_closures)   gPad->SetLogy();
  JetShape_Stack_Up[0][3]->GetYaxis()->SetLabelSize(0.);
  JetShape_Stack_Down[0][3]->Draw("same");

 

  // JetShape_graph[0][3][9]->Draw("same e2 P");
  JetShape2[0][3][9]->Draw("same");

  // if(!is_subleading) JetShape_ref[0][3][9]->Draw("same");

  TLatex  *label_per = new TLatex(0.05,0.9,"HYD. Cent. 50-100%");
  label_per->SetTextSize(0.09);
  label_per->SetLineColor(kWhite);
  label_per->SetNDC();
  label_per->Draw();


  PAS_plot->cd(3);
 
  
  JetShape_Stack_Up[0][0]->Draw();

  JetShape_Stack_Up[0][0]->GetXaxis()->SetRangeUser(0.,0.99);
   if(use_highpT_bin&&!is_number&&!do_closures)   gPad->SetLogy();
  JetShape_Stack_Up[0][0]->GetYaxis()->SetLabelSize(0.);

  JetShape_Stack_Down[0][0]->Draw("same");

  JetShape_Stack_Up[0][0]->GetXaxis()->SetRangeUser(0.,0.99);
  JetShape_Stack_Up[0][0]->GetXaxis()->SetNdivisions(505);

   
  //  JetShape_graph[0][0][9]->Draw("same e2 P");
   JetShape2[0][0][9]->Draw("same");

   //if(!is_subleading) JetShape_ref[0][0][9]->Draw("same");
 
  TLatex  *label_cent = new TLatex(0.05,0.9,"HYD. Cent. 0-10%");
  label_cent->SetTextSize(0.09);
  label_cent->SetLineColor(kWhite);
  label_cent->SetNDC();
  label_cent->Draw();

 
  PAS_plot->cd(5);

  JetShape_ratio[0][3][9]->GetXaxis()->SetRangeUser(0.,0.99);
 JetShape_ratio[0][3][9]->Draw();
 JetShape_ratio[0][3][9]->GetYaxis()->SetLabelSize(0.08); 
 JetShape_ratio[0][3][9]->GetYaxis()->CenterTitle();
 if(is_number) JetShape_ratio[0][3][9]->GetYaxis()->SetTitle("Y_{HYD} - Y_{PYTH}");
 else JetShape_ratio[0][3][9]->GetYaxis()->SetTitle("#rho(r)_{HYD}/#rho(r)_{PYTH}");
 JetShape_ratio[0][3][9]->GetYaxis()->SetTitleSize(0.0);
 JetShape_ratio[0][3][9]->GetYaxis()->SetLabelSize(0.0);
 JetShape_ratio[0][3][9]->GetYaxis()->SetTitleOffset(1.);

 JetShape_ratio[0][3][9]->GetXaxis()->SetLabelSize(0.08); 
 JetShape_ratio[0][3][9]->GetXaxis()->CenterTitle();
 JetShape_ratio[0][3][9]->GetXaxis()->SetTitle("#Deltar");
 JetShape_ratio[0][3][9]->GetXaxis()->SetTitleSize(0.09);
 JetShape_ratio[0][3][9]->GetXaxis()->SetTitleOffset(0.6);
 JetShape_ratio[0][3][9]->GetXaxis()->CenterTitle();
 JetShape_ratio[0][3][9]->GetXaxis()->SetNdivisions(505);
 

 JetShape_ratio[0][3][9]->SetMarkerStyle(20);
 JetShape_ratio[0][3][9]->SetMarkerSize(1);
 JetShape_ratio[0][3][9]->SetMarkerColor(kBlack);
 
 if(is_number){
   JetShape_Diff_Stack_Up[0][3]->Draw("same");
   JetShape_Diff_Stack_Down[0][3]->Draw("same");
   JetShape_ratio[0][3][9]->Draw("same");
 }

  //  JetShape_graph[0][3][9]->Draw("same e2 P");
  JetShape_ratio[0][3][9]->Draw("same");
  // if(!is_subleading) JetShape_ref_ratio[0][3][9]->Draw("same");


  TLatex  *label_ratio = new TLatex(0.05,0.9,"");
  label_ratio->SetTextSize(0.09);
  label_ratio->SetLineColor(kWhite);
  label_ratio->SetNDC();
  label_ratio->Draw();

  TLegend *legend_ratio = new TLegend(0.02,0.7,0.9,0.9);
  if(is_number)legend_ratio->AddEntry(JetShape_ratio[0][3][9],"HYD - PYTH");
  else legend_ratio->AddEntry(JetShape_ratio[0][3][9],"HYD / PYTH");
  //  if(!is_subleading)  legend_ratio->AddEntry(JetShape_ref_ratio[0][3][9],"PLB 730 (2014)");
  legend_ratio->SetLineColor(kWhite);
  legend_ratio->SetFillColor(kWhite);
  legend_ratio->SetTextSize(0.08);
  legend_ratio->Draw("same");

  if(!is_number){
  TLine *line = new TLine(0.,1.,1.,1.);
  line->SetLineStyle(2);
  line->Draw();
  }else{
 TLine *line = new TLine(0.,0.,1.,0.);
  line->SetLineStyle(2);
  line->Draw();

  }
  
   TPave *cover_x_r = new TPave(0.9,0.,1.1,0.13);
   cover_x_r->SetFillColor(kWhite);
   cover_x_r->SetLineColor(kWhite);
   cover_x_r->SetOption("NDC NB");

   //   cover_x_r->Draw();

  TPave *cover_x_l = new TPave(-.0175,0.,0.1,0.125);
   cover_x_l->SetFillColor(kWhite);
   cover_x_l->SetLineColor(kWhite);
   cover_x_l->SetOption("NDC NB");


   //  cover_x_l->Draw();
 


  PAS_plot->cd(6);

 JetShape_ratio[0][0][9]->SetMarkerStyle(20);
 JetShape_ratio[0][0][9]->SetMarkerSize(1);
 JetShape_ratio[0][0][9]->SetMarkerColor(kBlack);
 

 JetShape_ratio[0][0][9]->GetXaxis()->SetRangeUser(0.,.99);
 JetShape_ratio[0][0][9]->Draw();
 JetShape_ratio[0][0][9]->GetYaxis()->SetLabelSize(0.0); 
  JetShape_ratio[0][0][9]->GetXaxis()->SetLabelSize(0.08); 
 JetShape_ratio[0][0][9]->GetXaxis()->CenterTitle();
 JetShape_ratio[0][0][9]->GetXaxis()->SetTitle("#Deltar");
 JetShape_ratio[0][0][9]->GetXaxis()->SetTitleSize(0.09);
 JetShape_ratio[0][0][9]->GetXaxis()->SetTitleOffset(.6);
 JetShape_ratio[0][0][9]->GetXaxis()->CenterTitle();
 JetShape_ratio[0][0][9]->GetXaxis()->SetNdivisions(505);

 
 // JetShape_ratio[0][0][9]->Draw();
  //  JetShape_graph[0][0][9]->Draw("same e2 P");
 
 if(is_number){
   JetShape_Diff_Stack_Up[0][0]->Draw("same");
   JetShape_Diff_Stack_Down[0][0]->Draw("same");
   JetShape_ratio[0][0][9]->Draw("same");
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
  


 cover_x_l->Draw();
   cover_x_r->Draw();


   PAS_plot->cd(4);
   if(is_number){
     TGaxis *dummy_axis_jetshape = new TGaxis(1.,0.13,1.0,.975,ratio_axis_min,ratio_axis_max);

     dummy_axis_jetshape->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetYaxis());
     dummy_axis_jetshape->SetTitleOffset(1.);
     dummy_axis_jetshape->SetTickSize(0.);
     dummy_axis_jetshape->CenterTitle();
     dummy_axis_jetshape->SetTitleSize(0.08);
     dummy_axis_jetshape->SetLabelSize(0.06);
     dummy_axis_jetshape->SetTitle("Y_{HYD} - Y_{PYTH}");
     dummy_axis_jetshape->Draw();
   }else{
     TGaxis *dummy_axis_jetshape = new TGaxis(1.,0.13,1.0,.975,ratio_axis_min,ratio_axis_max);
     dummy_axis_jetshape->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetYaxis());
     dummy_axis_jetshape->SetTitleOffset(1.);
     dummy_axis_jetshape->SetTickSize(0.);
     dummy_axis_jetshape->CenterTitle();
     dummy_axis_jetshape->SetTitleSize(0.08);
     dummy_axis_jetshape->SetLabelSize(0.08);
     dummy_axis_jetshape->SetTitle("#rho(#Deltar)_{HYD}/#rho(#Deltar)_{PYTH}");
     dummy_axis_jetshape->Draw();
   }
 
  TGaxis *dummy_axis_r = new TGaxis(0.,1.,1.,1.,0.,dr_axis_max);

  dummy_axis_r->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetXaxis());
  dummy_axis_r->SetTitleOffset(0.6);
  dummy_axis_r->SetTitleSize(0.09);
  dummy_axis_r->SetTickSize(0.);
  dummy_axis_r->SetNdivisions(505);
  //dummy_axis_r->Draw();
 
TPave *cover_x_b = new TPave(0.9,0.,0.995,0.17);
   cover_x_b->SetFillColor(kWhite);
   cover_x_b->SetLineColor(kWhite);
   cover_x_b->SetOption("NDC NB");
   cover_x_b->Draw();

   PAS_plot->cd(4);
   //   dummy_axis_jetshape->Draw();
   // dummy_axis_r->Draw();
 

   if(!use_highpT_bin){
     PAS_plot->SaveAs((TString)("JetShapes_PAS"+data_mc_type_strs[type_code]+do_closures_string+".pdf"));
     PAS_plot->SaveAs((TString)("JetShapes_PAS"+data_mc_type_strs[type_code]+do_closures_string+".png"));
   }else{
     if(!is_number){
       PAS_plot->SaveAs((TString)("JetShapes_WithHighpT_pTweighted"+data_mc_type_strs[type_code]+do_closures_string+".pdf"));
       PAS_plot->SaveAs((TString)("JetShapes_WithHighpT_pTweighted"+data_mc_type_strs[type_code]+do_closures_string+".png"));
     }else{
       PAS_plot->SaveAs((TString)("ParticleYield_by_dR"+data_mc_type_strs[type_code]+do_closures_string+".pdf"));
       PAS_plot->SaveAs((TString)("ParticleYield_by_dR"+data_mc_type_strs[type_code]+do_closures_string+".png"));

     }
   }

   if(is_number){


     for(int ibin = 0; ibin < nCBins; ibin++){

       for(int g = 0; g< 2; g++){

	 if(g==0){

	   Integral_Pt[g][ibin] = new TH1D((TString)("Integral_PYTHIA_HYDJET"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]),"",9,pTbins);
         
	 }else{

	   if(ibin > 0 )continue;

	   Integral_Pt[g][ibin] = new TH1D("Integral_Phi_pp","",9,pTbins);
	 }
     
	 for(int ibin3 = 0; ibin3<nTrkPtBins-1; ibin3++){
      
	   Integral_Pt[g][ibin]->SetBinContent(ibin3+2,integral[g][ibin][ibin3]);
	   Integral_Pt[g][ibin]->SetBinError(ibin3+2,integral_err[g][ibin][ibin3]);
	
	 }

	 if(g==0){
	   Integral_Pt[g][ibin]->SetMarkerStyle(10);
	 }else{
	   Integral_Pt[g][ibin]->SetMarkerStyle(4);

	 }
	 Integral_Pt[g][ibin]->SetMarkerSize(2);

	 Integral_Pt[g][ibin]->SetLineColor(kBlack);
	 Integral_Pt[g][ibin]->SetMarkerColor(kBlack);
  
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
       Integral_Pt[1][0]->Draw("same");
      
       Integral_Pt[0][ibin]->Draw();
       Integral_Pt[1][0]->Draw("same");
   
       if(ibin==3){
	 TLegend *int_legend = new TLegend(0.18,0.6,0.9,0.85);
	 int_legend->AddEntry( Integral_Pt[0][ibin],"HYDJET integral |#Delta#eta|<1.0, |#Delta#phi|<1.0");
	 int_legend->AddEntry( Integral_Pt[1][0],"PYTHIA integral |#Delta#eta|<1.0, |#Delta#phi|<1.0");
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
      
       Integral_diff_Pt[0][ibin]->SetMaximum(8.5);
       Integral_diff_Pt[0][ibin]->SetMinimum(-2.);
       Integral_diff_Pt[0][ibin]->Draw();

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
       labels->AddText((TString)("HYDJET (" + CBin_labels[ibin]+") minus PYTHIA"));
   
       labels->Draw("same");
  
       int_zero->Draw();

       Integral_Pt[0][ibin]->Write();
       Integral_Pt[1][0]->Write();
       Integral_diff_Pt[0][ibin]->Write();


     }

     c_integral->SaveAs((TString)("Integral_dR"+data_mc_type_strs[type_code]+do_closures_string+".png"));
     c_integral->SaveAs((TString)("Integral_dR"+data_mc_type_strs[type_code]+do_closures_string+".pdf"));


   }
  
  return 0;
}// main loop
