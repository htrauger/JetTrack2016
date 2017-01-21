#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include "TRandom.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLatex.h"
#include "THStack.h"
#include "TGaxis.h"


#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "../JetTrack2016_functions.h"


using namespace std;

Int_t particle_yields(bool is_number = kTRUE){

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.25);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);

  TCanvas *c_yields_phi;
  TCanvas *c_PbPb_pp_phi;

  TCanvas *c_yields_eta;
  TCanvas *c_PbPb_pp_eta;



  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 9;



  TH2D *result[12][nTrkPtBins][nCBins];

  THStack *signal_dEta_stack[12][nCBins];
  THStack *signal_dPhi_stack[12][nCBins];

  THStack *signal_dEta_diff_stack_up[12][nCBins];
  THStack *signal_dEta_diff_stack_down[12][nCBins];
  THStack *signal_dPhi_diff_stack_up[12][nCBins];
  THStack *signal_dPhi_diff_stack_down[12][nCBins];


  TH1D *signal_dPhi[12][nTrkPtBins][nCBins];
  TH1D *signal_dPhi_syst[12][nTrkPtBins][nCBins];
  TH1D *signal_dPhi_rebin[12][nTrkPtBins][nCBins];
  TH1D *signal_dPhi_PbPb_pp[12][nTrkPtBins][nCBins];  
  TH1D *signal_dPhi_PbPb_pp_syst[12][nTrkPtBins][nCBins];  
  TH1D *signal_dPhi_PbPb_pp_up[12][nTrkPtBins][nCBins];  
  TH1D *signal_dPhi_PbPb_pp_down[12][nTrkPtBins][nCBins];  

  TH1D *background_diff_rebin[12][nTrkPtBins][nCBins];
  TH1D *background_syst_rebin[12][nTrkPtBins][nCBins];
  TH1D *background_diff_PbPb_pp[12][nTrkPtBins][nCBins];

  TH1D *signal_dEta[12][nTrkPtBins][nCBins];
  TH1D *signal_dEta_syst[12][nTrkPtBins][nCBins];
  TH1D *signal_dEta_rebin[12][nTrkPtBins][nCBins];
  TH1D *signal_dEta_PbPb_pp[12][nTrkPtBins][nCBins];  
  TH1D *signal_dEta_PbPb_pp_syst[12][nTrkPtBins][nCBins];  
  TH1D *signal_dEta_PbPb_pp_up[12][nTrkPtBins][nCBins];  
  TH1D *signal_dEta_PbPb_pp_down[12][nTrkPtBins][nCBins];  

  TH1D *spill_over_dEta[12][nTrkPtBins][nCBins];
  TH1D *spill_over_dPhi[12][nTrkPtBins][nCBins];
 
  TH1D *jff_residual_dEta[12][nTrkPtBins][nCBins];
  TH1D *jff_residual_dPhi[12][nTrkPtBins][nCBins];

  TH1D *Integral_phi_Pt[12][5];
  TH1D *Integral_eta_Pt[12][5];
  TH1D *Integral_diff_Pt[12][5];

 Double_t integral_eta[12][nTrkPtBins][nCBins];
 Double_t integral_eta_err[12][nTrkPtBins][nCBins];
 
 Double_t integral_phi[12][nTrkPtBins][nCBins];
 Double_t integral_phi_err[12][nTrkPtBins][nCBins];

 
  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

   float TrkPtBins[nTrkPtBins+1] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
  


  TFile *f_in_pbpb = new TFile("../me_correct/PbPb_Inclusive_Correlations.root");

  TFile *f_in_ref = new TFile("../HIN_14_016_comparison/Inclusive_Data_AllPlots.root");
 
  if(!is_number){

    cout<<"We don't study weighted correlations by dPhi"<<endl;
    return -1;
  }
 
  TFile *f_in_pp = new TFile("../me_correct/pp_Inclusive_Correlations.root");
 
  TFile *f_spillover = new TFile("../spill_over/Inclusive_Hydjet_SpillOvers.root");
     
  TFile *f_jff_pyth = new TFile("../jff_residual/Inclusive_Pythia_JFFResiduals.root");
  TFile *f_jff_hyd = new TFile("../jff_residual/Inclusive_Hydjet_JFFResiduals.root");
 
  if(!is_number){
    f_jff_pyth = new TFile("../jff_residual/Inclusive_Pythia_JFFResiduals_pTweighted.root");
    f_jff_hyd = new TFile("../jff_residual/Inclusive_Hydjet_JFFResiduals_pTweighted.root");
 

  }

  TFile *f_out;
  
  if(is_number)   f_out = new TFile("Particle_Yields.root","RECREATE");
  else f_out = new TFile("Particle_Yields_pTweighted.root","RECREATE");

  TString in_name, pTlabel,centlabel,Ajlabel;

  int lbin, rbin;

  double bin_width_phi, bin_width_eta, diff_max, diff_min, signal_min, signal_max, stacked_max, stacked_max_diff, stacked_min_diff, stacked_min, bc, err;

  float rel_err_pbpb =TMath::Sqrt(0.05*0.05+0.01*0.01+0.04*0.04);
  float rel_err_pp =TMath::Sqrt(0.05*0.05+0.04*0.04+0.03*0.03);

  TF1 *fit_line_left = new TF1("fit_line_left","[0]");
  TF1 *fit_line_right = new TF1("fit_line_right","[0]");

  TCanvas *dummy = new TCanvas("dummy");

  TPaveText *labels;

  float me_err[12][6][5];
  float bg_err[12][6][5];

  TH2D *bg_err_pp = new TH2D("bg_err_pp","",9,0.5,8.5,4,0.5,3.5);

  TH2D *bg_err_pbpb = new TH2D("bg_err_PbPb","",9,0.5,8.5,4,0.5,3.5);

  TH2D *me_err_pp = new TH2D("me_err_pp","",9,0.5,8.5,4,0.5,3.5);

  TH2D *me_err_pbpb = new TH2D("me_err_PbPb","",9,0.5,8.5,4,0.5,3.5);



  c_yields_phi = new TCanvas("yields_phi","",10,10,1600,2400);
  c_yields_phi->Divide(4,8,0,0);

  c_PbPb_pp_phi = new TCanvas("yields_phi_diff","",10,10,1600,2400);
  c_PbPb_pp_phi->Divide(4,8,0,0);
    
  c_yields_eta = new TCanvas("yields_eta","",10,10,1600,2400);
  c_yields_eta->Divide(4,8,0,0);

  c_PbPb_pp_eta = new TCanvas("yields_eta_diff","",10,10,1600,2400);
  c_PbPb_pp_eta->Divide(4,8,0,0);
 

  
  for(int i = 0; i<nTrkPtBins-1; i++){

    for(int j = 0; j<4; j++){

      for(int g = 0; g<2; g++){


	if(g==0){

	  cout<<"Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]<<endl;

	  if(is_number)	  result[g][i][j] = (TH2D*)f_in_pbpb->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
	  else	  result[g][i][j] = (TH2D*)f_in_pbpb->Get((TString)("Yield_BkgSub_pTweighted" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));


	    
	}else{
	  
	  if(j>0)continue;

	  if(is_number) result[g][i][j] = (TH2D*)f_in_pp->Get((TString)("Yield_BkgSub_Cent0_Cent10_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_pp_Cent0_Cent10_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
	  else  result[g][i][j] = (TH2D*)f_in_pp->Get((TString)("Yield_BkgSub_pTweightedCent0_Cent10_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_pp_pTweightedCent0_Cent10_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
	}
	cout<<"got hist"<<endl;      
	  
	if(i==0) result[g][i][j]->Scale(1./.3);
	if(i> 3&& i<nTrkPtBins) result[g][i][j]->Scale(1./4);

	lbin = result[g][i][j]->GetXaxis()->FindBin(-etalim+.0001);
	rbin = result[g][i][j]->GetXaxis()->FindBin(etalim-.0001);

	  
	if(g==0)	signal_dPhi[g][i][j] = (TH1D*)result[g][i][j]->ProjectionY((TString)("Proj_dPhi_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);
	else	signal_dPhi[g][i][j] = (TH1D*)result[g][i][j]->ProjectionY((TString)("Proj_dPhi_pp_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);


	int nbins_half = result[g][i][j]->GetNbinsY()/2;

	for(int k = 1; k<nbins_half/2+1; k++){

	  bc = (signal_dPhi[g][i][j]->GetBinContent(k)+signal_dPhi[g][i][j]->GetBinContent(nbins_half+1-k))/2.;
	  err = TMath::Sqrt(signal_dPhi[g][i][j]->GetBinError(k)*signal_dPhi[g][i][j]->GetBinError(k)+signal_dPhi[g][i][j]->GetBinError(nbins_half+1-k)*signal_dPhi[g][i][j]->GetBinError(nbins_half+1-k))/2.;
	  signal_dPhi[g][i][j]->SetBinContent(k,bc);
	  signal_dPhi[g][i][j]->SetBinContent(nbins_half+1-k,bc);
	  signal_dPhi[g][i][j]->SetBinError(k,err);
	  signal_dPhi[g][i][j]->SetBinError(nbins_half+1-k,err);
	}


	if( signal_dPhi[g][i][j]->GetBinWidth(1)!= signal_dPhi[g][i][j]->GetBinWidth(1))cout<<"Widths do not match!"<<endl;
	  
	bin_width_phi =  signal_dPhi[g][i][j]->GetBinWidth(1);

	signal_dPhi[g][i][j]->Scale(1./bin_width_phi);
	
	signal_dPhi_rebin[g][i][j] = Rebin_dPhi(signal_dPhi[g][i][j]);

	lbin = result[g][i][j]->GetYaxis()->FindBin(-philim+.0001);
	rbin = result[g][i][j]->GetYaxis()->FindBin(philim-.0001);


	  
	if(g==0)	signal_dEta[g][i][j] = (TH1D*)result[g][i][j]->ProjectionX((TString)("Proj_dEta_PbPb_"+  CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);
	else 	signal_dEta[g][i][j] = (TH1D*)result[g][i][j]->ProjectionX((TString)("Proj_dEta_pp_"+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);

	bin_width_eta =  signal_dEta[g][i][j]->GetBinWidth(1);

	signal_dEta[g][i][j]->Scale(1./bin_width_eta);

	int nbins_eta = result[g][i][j]->GetNbinsX();

	  for(int k = 1; k<nbins_eta/2+1; k++){

	  bc = (signal_dEta[g][i][j]->GetBinContent(k)+signal_dEta[g][i][j]->GetBinContent(nbins_eta+1-k))/2.;
	  err = TMath::Sqrt(signal_dEta[g][i][j]->GetBinError(k)*signal_dEta[g][i][j]->GetBinError(k)+signal_dEta[g][i][j]->GetBinError(nbins_eta+1-k)*signal_dEta[g][i][j]->GetBinError(nbins_eta+1-k))/2.;
	  signal_dEta[g][i][j]->SetBinContent(k,bc);
	  signal_dEta[g][i][j]->SetBinContent(nbins_eta+1-k,bc);
	  signal_dEta[g][i][j]->SetBinError(k,err);
	  signal_dEta[g][i][j]->SetBinError(nbins_eta+1-k,err);

	  }

      
	signal_dEta_rebin[g][i][j] = Rebin_dEta(signal_dEta[g][i][j]);

	//------------------
	//Apply corrections
	//-------------------

	
	if(g==1){
	  
	  
	  TString	jff_name =(TString)("JFF_Residual_Eta_Cent0_Cent10_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]);
	  jff_residual_dEta[g][i][j] = (TH1D*)f_jff_pyth->Get(jff_name)->Clone(jff_name);

	  jff_name.ReplaceAll("Eta","Phi");
	  jff_residual_dPhi[g][i][j] = (TH1D*)f_jff_pyth->Get(jff_name)->Clone(jff_name);

	  signal_dEta_rebin[g][i][j]->Add(jff_residual_dEta[g][i][0],-1.);
	  signal_dPhi_rebin[g][i][j]->Add(jff_residual_dPhi[g][i][0],-1.);
	    
	
	}else{

	  cout<<"starting corrs"<<endl;
	  if(i>4){
	    /*
	    TString	jff_name =(TString)("JFF_Residual_Eta_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]);
	    jff_residual_dEta[g][i][j] = (TH1D*)f_jff_hyd->Get(jff_name)->Clone(jff_name);
	    */	    
	    TString	jff_name =(TString)("JFF_Residual_Eta_"+CBin_strs[0]+"_"+CBin_strs[1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]);
	    jff_residual_dEta[g][i][j] = (TH1D*)f_jff_hyd->Get(jff_name)->Clone(jff_name);
	    


	    jff_name.ReplaceAll("Eta","Phi");
	    jff_residual_dPhi[g][i][j] = (TH1D*)f_jff_hyd->Get(jff_name)->Clone(jff_name);

	    signal_dEta_rebin[g][i][j]->Add(jff_residual_dEta[g][i][j],-1.);
	    signal_dPhi_rebin[g][i][j]->Add(jff_residual_dPhi[g][i][j],-1.);
	  
	  }else{
	    /*
	    TString	jff_name =(TString)("JFF_Residual_Eta_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]);
	    jff_residual_dEta[g][i][j] = (TH1D*)f_jff_hyd->Get(jff_name)->Clone(jff_name);

	    */

	    TString	jff_name =(TString)("JFF_Residual_Eta_"+CBin_strs[0]+"_"+CBin_strs[1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]);
	    jff_residual_dEta[g][i][j] = (TH1D*)f_jff_hyd->Get(jff_name)->Clone(jff_name);
	    
	    
	    jff_name.ReplaceAll("Eta","Phi");

	    jff_residual_dPhi[g][i][j] = (TH1D*)f_jff_hyd->Get(jff_name)->Clone(jff_name);
	   
	    signal_dEta_rebin[g][i][j]->Add(jff_residual_dEta[g][i][j],-1.);
	    signal_dPhi_rebin[g][i][j]->Add(jff_residual_dPhi[g][i][j],-1.);
	    


	    TString	spillover_name =(TString)("Eta_SpillOver_Points_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]);
	    spill_over_dEta[g][i][j] = (TH1D*)f_spillover->Get(spillover_name)->Clone(spillover_name);
	    
	    cout<<"got one"<<endl;
	    spillover_name.ReplaceAll("Eta","Phi");
	    spill_over_dPhi[g][i][j] = (TH1D*)f_spillover->Get(spillover_name)->Clone(spillover_name);
	  
	    if(i<3){
	      signal_dEta_rebin[g][i][j]->Add(spill_over_dEta[g][i][j],-1.);
	      signal_dPhi_rebin[g][i][j]->Add(spill_over_dPhi[g][i][j],-1.);
	    }
	  }


	
	}
	//---------------------------------------------------
	//---------------------------------------------------
	cout<<"done corrs"<<endl;



	signal_dEta_rebin[g][i][j]->SetAxisRange(-2.49,2.49);
	signal_dPhi_rebin[g][i][j]->SetAxisRange(-1.49,1.49);

	signal_dEta_rebin[g][i][j]->SetMarkerSize(1);
	signal_dPhi_rebin[g][i][j]->SetMarkerSize(1);


	if(g==0){
	  signal_dPhi_rebin[g][i][j]->SetMarkerColor(kBlack);
	  signal_dPhi_rebin[g][i][j]->SetLineColor(kBlack);
	  signal_dEta_rebin[g][i][j]->SetMarkerColor(kBlack);
	  signal_dEta_rebin[g][i][j]->SetLineColor(kBlack);

	  signal_dEta_rebin[g][i][j]->SetMarkerStyle(10);
	  signal_dPhi_rebin[g][i][j]->SetMarkerStyle(10);

	}else{
	  signal_dPhi_rebin[g][i][j]->SetMarkerColor(kBlack);
	  signal_dPhi_rebin[g][i][j]->SetLineColor(kBlack);
	  signal_dEta_rebin[g][i][j]->SetMarkerColor(kBlack);
	  signal_dEta_rebin[g][i][j]->SetLineColor(kBlack);

	  signal_dEta_rebin[g][i][j]->SetMarkerStyle(24);
	  signal_dPhi_rebin[g][i][j]->SetMarkerStyle(24);
	  
	}
      
	cout<<"starting to get raw hists for dEta calc"<<endl;

	if(g==0){

	  if(is_number)	  result[g][i][j] = (TH2D*)f_in_pbpb->Get((TString)("Yield_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
	  else	  result[g][i][j] = (TH2D*)f_in_pbpb->Get((TString)("Yield_pTweighted_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
	    
	}else{
	  
	  if(is_number) result[g][i][j] = (TH2D*)f_in_pp->Get((TString)("Yield_pp_Cent0_Cent10_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_pp_Cent0_Cent10_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
	  else  result[g][i][j] = (TH2D*)f_in_pp->Get((TString)("Yield_pTweighted_pp_Cent0_Cent10_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_pp_pTweightedCent0_Cent10_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
	}
	  
	if(i==0) result[g][i][j]->Scale(1./.3);
	if(i> 3&& i<nTrkPtBins) result[g][i][j]->Scale(1./4);

	lbin = result[g][i][j]->GetXaxis()->FindBin(-etalim+.0001);
	rbin = result[g][i][j]->GetXaxis()->FindBin(etalim-.0001);

	  
	if(g==0)	signal_dPhi[g][i][j] = (TH1D*)result[g][i][j]->ProjectionY((TString)("Proj_dPhi_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);
	else	signal_dPhi[g][i][j] = (TH1D*)result[g][i][j]->ProjectionY((TString)("Proj_dPhi_pp_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);

	if( signal_dPhi[g][i][j]->GetBinWidth(1)!= signal_dPhi[g][i][j]->GetBinWidth(1))cout<<"Widths do not match!"<<endl;
	  
	bin_width_phi =  signal_dPhi[g][i][j]->GetBinWidth(1);

	signal_dPhi[g][i][j]->Scale(1./bin_width_phi);
	

	lbin = result[g][i][j]->GetYaxis()->FindBin(-philim+.0001);
	rbin = result[g][i][j]->GetYaxis()->FindBin(philim-.0001);


	  
	if(g==0)	signal_dEta[g][i][j] = (TH1D*)result[g][i][j]->ProjectionX((TString)("Proj_dEta_PbPb_"+  CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);
	else 	signal_dEta[g][i][j] = (TH1D*)result[g][i][j]->ProjectionX((TString)("Proj_dEta_pp_"+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);

	bin_width_eta =  signal_dEta[g][i][j]->GetBinWidth(1);

	signal_dEta[g][i][j]->Scale(1./bin_width_eta);
    

	float avg_dev = ( (  signal_dEta_rebin[g][i][j] ->GetBinContent(  signal_dEta_rebin[g][i][j]->FindBin(-1.51))- signal_dEta_rebin[g][i][j] ->GetBinContent(  signal_dEta_rebin[g][i][j]->FindBin(-2.01)))+(signal_dEta_rebin[g][i][j] ->GetBinContent(  signal_dEta_rebin[g][i][j]->FindBin(1.51))- signal_dEta_rebin[g][i][j] ->GetBinContent(  signal_dEta_rebin[g][i][j]->FindBin(2.01))))/2.;

	float max_bin =  max( TMath::Abs( signal_dEta_rebin[g][i][j] ->GetBinContent(  signal_dEta_rebin[g][i][j]->FindBin(-1.49))), TMath::Abs(signal_dEta_rebin[g][i][j] ->GetBinContent(  signal_dEta_rebin[g][i][j]->FindBin(1.49))));


	bg_err[g][i][j] = avg_dev;
	cout<<g<<" "<<i<<" "<<j<<" "<<max_bin<<" "<<avg_dev<<" "<<bg_err[g][i][j] <<endl;
	
	cout<<"here 0"<<endl;

	float line_left, line_right; 
	

	dummy->cd(0);
       	signal_dEta[g][i][j]->Fit("fit_line_left","","",-2.5,-1.5);
	
	line_left = fit_line_left->GetParameter(0);

	signal_dEta[g][i][j]->Fit("fit_line_right","","",1.5,2.5);
	
	line_right = fit_line_right->GetParameter(0);

	me_err[g][i][j] = TMath::Abs(line_left - line_right);

	if(g==0)	bg_err_pbpb->SetBinContent(i+1,j+1,bg_err[g][i][j]);
	else	bg_err_pp->SetBinContent(i+1,j+1,bg_err[g][i][j]);

	if(g==0)	me_err_pbpb->SetBinContent(i+1,j+1,me_err[g][i][j]);
	else	me_err_pp->SetBinContent(i+1,j+1,me_err[g][i][j]);

	
	cout<<bg_err[g][i][j]<<" "<<me_err[g][i][j]<<endl;


      }
	//---------------------------------------------------
	//---------------------------------------------------

	// CALCULATE SYST ERR

	//-------------------


	signal_dEta_syst[0][i][j] = (TH1D*)signal_dEta_rebin[0][i][j]->Clone((TString)("dEta_Syst_PbPb_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	signal_dEta_syst[1][i][j] = (TH1D*)signal_dEta_rebin[1][i][0]->Clone((TString)("dEta_Syst_pp_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	signal_dPhi_syst[0][i][j] = (TH1D*)signal_dPhi_rebin[0][i][j]->Clone((TString)("dPhi_Syst_PbPb_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	signal_dPhi_syst[1][i][j] = (TH1D*)signal_dPhi_rebin[1][i][0]->Clone((TString)("dPhi_Syst_pp_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

  
	for(int k = 1; k<signal_dPhi_syst[0][i][j]->GetNbinsX()+1; k++){

	  bc =  signal_dPhi_syst[0][i][j]->GetBinContent(k);

	  if(i<4) err = TMath::Sqrt(bc*rel_err_pbpb*bc*rel_err_pbpb+jff_residual_dPhi[0][i][j]->GetBinContent(k)*jff_residual_dPhi[0][i][j]->GetBinContent(k)/4.+spill_over_dPhi[0][i][j]->GetBinContent(k)*spill_over_dPhi[0][i][j]->GetBinContent(k)*.58*.58+me_err[0][i][j]*me_err[0][i][j]+bg_err[0][i][j]*bg_err[0][i][j]);

	  else   err = TMath::Sqrt(bc*rel_err_pbpb*bc*rel_err_pbpb+jff_residual_dPhi[0][i][j]->GetBinContent(k)*jff_residual_dPhi[0][i][j]->GetBinContent(k)/4.+me_err[0][i][j]*me_err[0][i][j]+bg_err[0][i][j]*bg_err[0][i][j]);

	  signal_dPhi_syst[0][i][j]->SetBinError(k,err);

	  cout<<i<<" "<<j<<" "<<k<<" "<<bc<<" "<<bc*rel_err_pbpb<<" "<<jff_residual_dPhi[0][i][j]->GetBinContent(k)/2.<<" "<<me_err[0][i][j]<<" "<<bg_err[0][i][j]<<endl;
	  
	  bc =  signal_dPhi_syst[1][i][0]->GetBinContent(k);

	  err =TMath::Sqrt(rel_err_pp*bc*rel_err_pp*bc+jff_residual_dPhi[1][i][0]->GetBinError(1)*jff_residual_dPhi[1][i][0]->GetBinError(1)/4.+me_err[1][i][0]*me_err[1][i][0]+bg_err[1][i][0]*bg_err[1][i][0]);

	  signal_dPhi_syst[1][i][0]->SetBinError(k,err);

	}


	for(int k = 1; k<signal_dEta_syst[0][i][j]->GetNbinsX()+1; k++){

	  bc =  signal_dEta_syst[0][i][j]->GetBinContent(k);
	  
	  if(i<4)	  err = TMath::Sqrt(bc*rel_err_pbpb*bc*rel_err_pbpb+jff_residual_dEta[0][i][j]->GetBinError(1)*jff_residual_dEta[0][i][j]->GetBinError(1)/4.+spill_over_dEta[0][i][j]->GetBinError(1)*spill_over_dEta[0][i][j]->GetBinError(1)*.58*.58+me_err[0][i][j]*me_err[0][i][j]+bg_err[0][i][j]*bg_err[0][i][j]);
	  else  err = TMath::Sqrt(bc*rel_err_pbpb*bc*rel_err_pbpb+jff_residual_dEta[0][i][j]->GetBinError(1)*jff_residual_dEta[0][i][j]->GetBinError(1)/4.+me_err[0][i][j]*me_err[0][i][j]+bg_err[0][i][j]*bg_err[0][i][j]);

	  signal_dEta_syst[0][i][j]->SetBinError(k,err);


	  bc =  signal_dEta_syst[1][i][0]->GetBinContent(k);

	  err =TMath::Sqrt(rel_err_pp*bc*rel_err_pp*bc+jff_residual_dEta[1][i][0]->GetBinError(1)*jff_residual_dEta[1][i][0]->GetBinError(1)+me_err[1][i][0]*me_err[1][i][0]+bg_err[1][i][0]*bg_err[1][i][0]);

	  signal_dEta_syst[1][i][0]->SetBinError(k,err);

	}
          
	cout<<"got syst err"<<endl;

      
      switch(i){
      case 0: 
	signal_max = 13.;
	signal_min = -1.;
	diff_max = 13.;
	diff_min = -1.;
	break;
      case 1: 
	signal_max = 8.;
	signal_min = -.5;
	diff_max = 8.;
	diff_min = -.5;
	break;
      case 2: 
	signal_max = 6.;
	signal_min = -.5;
	diff_max = 6.;
	diff_min = -.5;
	break;
      case 3:
	signal_max = 6.;
	signal_min = -.5; 
	diff_max = 6.;
	diff_min = -.5;
	break;
      case 4: 
	signal_max = 4.;
	signal_min = -.5;
	diff_max = 4.;
	diff_min = -.5;
	break;
      case 5: 
	signal_max = 2.;
	signal_min = -.5;
	diff_max = 2.;
	diff_min = -.5;
	break;
      default: 
	break;
      }


      stacked_min = -3.;
      stacked_min_diff = -2.;
      stacked_max = 35.;
      stacked_max_diff = 13.;
     
      c_yields_phi->cd(4*(i+1)-j);
     
      signal_dPhi_rebin[0][i][j]->SetMinimum(signal_min);
      signal_dPhi_rebin[0][i][j]->SetMaximum(signal_max);


      if(j==3){
	signal_dPhi_rebin[0][i][j]->GetYaxis()->SetLabelSize(0.08);
	signal_dPhi_rebin[0][i][j]->GetYaxis()->SetTitleSize(0.08);
	signal_dPhi_rebin[0][i][j]->GetYaxis()->SetTitle("1/N_{evt} 1/d#Delta#phi");
      }else{
	signal_dPhi_rebin[0][i][j]->GetYaxis()->SetLabelSize(0.0);
      }

      if(i==7){
	signal_dPhi_rebin[0][i][j]->GetXaxis()->SetLabelSize(0.08);
	signal_dPhi_rebin[0][i][j]->GetXaxis()->SetTitleSize(0.08);
	signal_dPhi_rebin[0][i][j]->GetXaxis()->SetTitle("#Delta#phi");
	signal_dPhi_rebin[0][i][j]->GetXaxis()->CenterTitle();

      }

      signal_dPhi_syst[0][i][j]->SetFillColor(kYellow);
      signal_dPhi_syst[1][i][0]->SetFillColor(kViolet);
    
      signal_dPhi_rebin[0][i][j]->Draw();
      signal_dPhi_syst[0][i][j]->Draw("same e2");
      signal_dPhi_syst[1][i][0]->Draw("same e2");
      signal_dPhi_rebin[1][i][0]->SetMarkerStyle(20);
      signal_dPhi_syst[1][i][0]->SetMarkerStyle(20);
      signal_dPhi_rebin[1][i][0]->Draw("same");
      signal_dPhi_rebin[0][i][j]->Draw("same");

      //------------
      // Drawing
      //------------

      TLine *zero = new TLine(-1.5,0.,1.5,0.);
      zero->SetLineStyle(2);
      zero->Draw();

      if(j==3){
	labels = new TPaveText(0.18,0.75,0.45,0.95,"NDC");
      }else{
	labels = new TPaveText(0.05,0.75,0.45,0.95,"NDC");
      }  
      labels->SetName("labels");
      labels->SetFillColor(0);
      labels->SetLineColor(0);
      labels->SetTextAlign(11);
      labels->AddText(CBin_labels[j]);
      labels->AddText(TrkPtBin_labels[i]);
      labels->SetTextSize(0.06);
      labels->Draw("same");

      if(i==0&&j==3){
	TLegend *legend = new TLegend(0.2,0.65,0.95,0.85);

	legend->AddEntry(signal_dPhi_rebin[0][i][j],"PbPb");
	legend->AddEntry(signal_dPhi_rebin[1][i][0],"pp");

	legend->SetTextSize(0.06);
	legend->SetLineColor(0);
	legend->Draw();
      }

      c_PbPb_pp_phi->cd(4*(i+1)-j);
      

      signal_dPhi_PbPb_pp[0][i][j] = (TH1D*)signal_dPhi_rebin[0][i][j]->Clone((TString)("Diff_dPhi_"+ CBin_strs[j] + "_" + CBin_strs[j+1] +"_"+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
      signal_dPhi_PbPb_pp[0][i][j]->Add(signal_dPhi_rebin[1][i][0],-1.);



      signal_dPhi_PbPb_pp_syst[0][i][j] = (TH1D*)signal_dPhi_syst[0][i][j]->Clone((TString)("Diff_dPhi_syst_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_"+TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
      signal_dPhi_PbPb_pp_syst[0][i][j]->Add(signal_dPhi_syst[1][i][0],-1.);


      signal_dPhi_PbPb_pp[0][i][j]->SetMaximum(diff_max);
      signal_dPhi_PbPb_pp[0][i][j]->SetMinimum(diff_min);

      signal_dPhi_PbPb_pp[0][i][j]->SetMarkerColor(kBlack);
      signal_dPhi_PbPb_pp[0][i][j]->SetLineColor(kBlack);


      if(j==3){
	signal_dPhi_PbPb_pp[0][i][j]->GetYaxis()->SetLabelSize(0.08);
	signal_dPhi_PbPb_pp[0][i][j]->GetYaxis()->SetTitleSize(0.08);
	signal_dPhi_PbPb_pp[0][i][j]->GetYaxis()->SetTitle("1/N_{evt} 1/d#Delta#phi");
      }else{
	signal_dPhi_PbPb_pp[0][i][j]->GetYaxis()->SetLabelSize(0.0);
      }

      if(i==7){
	signal_dPhi_PbPb_pp[0][i][j]->GetXaxis()->SetLabelSize(0.1);
	signal_dPhi_PbPb_pp[0][i][j]->GetXaxis()->SetTitleSize(0.1);
	signal_dPhi_PbPb_pp[0][i][j]->GetXaxis()->SetTitle("#Delta#phi");
	signal_dPhi_PbPb_pp[0][i][j]->GetXaxis()->CenterTitle();

      }
 

      signal_dPhi_PbPb_pp[0][i][j]->Draw();
      signal_dPhi_PbPb_pp_syst[0][i][j]->Draw("same e2");
      signal_dPhi_PbPb_pp[0][i][j]->Draw("same");

      zero->Draw();
      labels->SetLineColor(0);
      labels->Draw("same");

      if(j==3&&i==0){
	TLegend *legend_diff = new TLegend(0.2,0.65,0.95,0.85);

	legend_diff->AddEntry(signal_dPhi_PbPb_pp[0][i][j],"PbPb - pp");

	legend_diff->SetLineColor(0);
	legend_diff->SetTextSize(.06);
	legend_diff->Draw();

      }

      //_____________
      //  Draw dEta
      //-------------


      c_yields_eta->cd(4*(i+1)-j);

      signal_dEta_rebin[0][i][j]->SetMinimum(signal_min);
      signal_dEta_rebin[0][i][j]->SetMaximum(signal_max);

      if(j==3){
	signal_dEta_rebin[0][i][j]->GetYaxis()->SetLabelSize(0.08);
	signal_dEta_rebin[0][i][j]->GetYaxis()->SetTitleSize(0.08);
	signal_dEta_rebin[0][i][j]->GetYaxis()->SetTitle("1/N_{evt} 1/d#Delta#eta");
      }else{
	signal_dEta_rebin[0][i][j]->GetYaxis()->SetLabelSize(0.0);
      }

      if(i==7){
	signal_dEta_rebin[0][i][j]->GetXaxis()->SetLabelSize(0.08);
	signal_dEta_rebin[0][i][j]->GetXaxis()->SetTitleSize(0.08);
	signal_dEta_rebin[0][i][j]->GetXaxis()->SetTitle("#Delta#eta");
	signal_dEta_rebin[0][i][j]->GetXaxis()->CenterTitle();

      }
      signal_dEta_syst[0][i][j]->SetFillColor(kYellow);
      signal_dEta_syst[1][i][0]->SetFillColor(kViolet);
    
      signal_dEta_rebin[0][i][j]->Draw();
      signal_dEta_syst[0][i][j]->Draw("same e2");
      signal_dEta_syst[1][i][0]->Draw("same e2");
      signal_dEta_rebin[1][i][0]->SetMarkerStyle(20);
      signal_dEta_syst[1][i][0]->SetMarkerStyle(20);
      signal_dEta_rebin[1][i][0]->Draw("same");
      signal_dEta_rebin[0][i][j]->Draw("same");

      zero->Draw();
      labels->Draw("same");

      if(i==0&&j==3){
	TLegend *legend = new TLegend(0.2,0.65,0.95,0.85);

	legend->AddEntry(signal_dPhi_syst[0][i][j],"2.76 TeV ak3Calo pp");
	legend->AddEntry(signal_dPhi_syst[1][i][0],"5.02 TeV ak4Calo pp");
	legend->SetTextSize(0.06);
	legend->SetLineColor(0);
	legend->Draw();
      }

      c_PbPb_pp_eta->cd(4*(i+1)-j);
      //c_yields_eta->cd(4*i+j+3);



      signal_dEta_PbPb_pp[0][i][j] = (TH1D*)signal_dEta_rebin[0][i][j]->Clone((TString)("Diff_dEta_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_"+TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
      signal_dEta_PbPb_pp[0][i][j]->Add(signal_dEta_rebin[1][i][0],-1.);


      signal_dEta_PbPb_pp_syst[0][i][j] = (TH1D*)signal_dEta_syst[0][i][j]->Clone((TString)("Diff_dEta_syst_"+ CBin_strs[j] + "_" + CBin_strs[j+1] +"_"+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
      signal_dEta_PbPb_pp_syst[0][i][j]->Add(signal_dEta_syst[1][i][0],-1.);


      signal_dEta_PbPb_pp[0][i][j]->SetMaximum(diff_max);
      signal_dEta_PbPb_pp[0][i][j]->SetMinimum(diff_min);

   signal_dEta_PbPb_pp[0][i][j]->SetMarkerColor(kBlack);
      signal_dEta_PbPb_pp[0][i][j]->SetLineColor(kBlack);


      if(j==3){
	signal_dEta_PbPb_pp[0][i][j]->GetYaxis()->SetLabelSize(0.08);
	signal_dEta_PbPb_pp[0][i][j]->GetYaxis()->SetTitleSize(0.08);
	signal_dEta_PbPb_pp[0][i][j]->GetYaxis()->SetTitle("1/N_{evt} 1/d#Delta#eta");
      }else{
	signal_dEta_PbPb_pp[0][i][j]->GetYaxis()->SetLabelSize(0.0);
      }

      if(i==7){
	signal_dEta_PbPb_pp[0][i][j]->GetXaxis()->SetLabelSize(0.08);
	signal_dEta_PbPb_pp[0][i][j]->GetXaxis()->SetTitleSize(0.08);
	signal_dEta_PbPb_pp[0][i][j]->GetXaxis()->SetTitle("#Delta#eta");
	signal_dEta_PbPb_pp[0][i][j]->GetXaxis()->CenterTitle();

      }


      signal_dEta_PbPb_pp[0][i][j]->Draw();
      signal_dEta_PbPb_pp_syst[0][i][j]->Draw("same e2");
      signal_dEta_PbPb_pp[0][i][j]->Draw("same");

      zero->Draw();
      labels->SetLineColor(0);
      labels->Draw("same");

      if(j==3&&i==0){
	TLegend *legend_diff = new TLegend(0.2,0.55,0.95,0.75);

	legend_diff->AddEntry(signal_dEta_PbPb_pp[0][i][j],"PbPb - pp");
	legend_diff->SetLineColor(0);
	legend_diff->SetTextSize(0.06);
	legend_diff->Draw();

	dummy->cd();

      }

      //---------------------
      //  Write everything!
      //----------------------

      f_out->cd();

      signal_dPhi_rebin[0][i][j]->Write();
      signal_dPhi_rebin[1][i][0]->Write();
      signal_dEta_rebin[0][i][j]->Write();
      signal_dEta_rebin[1][i][0]->Write();


      llimiteta = signal_dEta_rebin[0][i][j]->FindBin(-etalim+.0001);
      rlimiteta = signal_dEta_rebin[0][i][j]->FindBin(etalim-.0001);


      integral_eta[0][i][j] = signal_dEta_rebin[0][i][j]->IntegralAndError(llimiteta, rlimiteta, integral_eta_err[0][i][j],"width");
      integral_eta[1][i][0] = signal_dEta_rebin[1][i][0]->IntegralAndError(llimiteta, rlimiteta, integral_eta_err[1][i][0],"width");

      float check = 0;
    
      llimitphi = signal_dPhi_rebin[0][i][j]->FindBin(-philim+.0001);
      rlimitphi = signal_dPhi_rebin[0][i][j]->FindBin(philim-.0001);

      integral_phi[0][i][j] = signal_dPhi_rebin[0][i][j]->IntegralAndError(llimitphi, rlimitphi, integral_phi_err[0][i][j],"width");
      integral_phi[1][i][0] = signal_dPhi_rebin[1][i][0]->IntegralAndError(llimitphi, rlimitphi, integral_phi_err[1][i][0],"width");
     

  
      signal_dEta_PbPb_pp_up[0][i][j] = (TH1D*)signal_dEta_PbPb_pp[0][i][j]->Clone((TString)("Diff_dEta_Up"+ CBin_strs[j] + "_" + CBin_strs[j+1] +"_"+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

      signal_dEta_PbPb_pp_down[0][i][j] = (TH1D*)signal_dEta_PbPb_pp[0][i][j]->Clone((TString)("Diff_dEta_Down"+ CBin_strs[j] + "_" + CBin_strs[j+1] +"_"+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));


      signal_dPhi_PbPb_pp_up[0][i][j] = (TH1D*)signal_dPhi_PbPb_pp[0][i][j]->Clone((TString)("Diff_dPhi_Up"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_"+TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

      signal_dPhi_PbPb_pp_down[0][i][j] = (TH1D*)signal_dPhi_PbPb_pp[0][i][j]->Clone((TString)("Diff_dPhi_Down"+ CBin_strs[j] + "_" + CBin_strs[j+1] +"_"+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
     
    }
  }


  c_yields_phi->SaveAs((TString)("Yields_dPhi_Inclusive.png"));
  c_PbPb_pp_phi->SaveAs((TString)("Yields_dPhi_PbPb_minus_pp_Inclusive.png"));
  c_yields_eta->SaveAs((TString)("Yields_dEta_Inclusive.png"));
  c_PbPb_pp_eta->SaveAs((TString)("Yields_dEta_PbPb_minus_pp_Inclusive.png"));


 c_yields_phi->SaveAs((TString)("Yields_dPhi_Inclusive.pdf"));
  c_PbPb_pp_phi->SaveAs((TString)("Yields_dPhi_PbPb_minus_pp_Inclusive.pdf"));
  c_yields_eta->SaveAs((TString)("Yields_dEta_Inclusive.pdf"));
  c_PbPb_pp_eta->SaveAs((TString)("Yields_dEta_PbPb_minus_pp_Inclusive.pdf"));


 

  // NOW MAKE AND DRAW STACKED PLOTS!





  for(int j = 0; j < nCBins; j++){
    

    signal_dEta_rebin[0][8][j] = (TH1D*) signal_dEta_rebin[0][0][j]->Clone((TString)("Combined_dEta_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
    signal_dPhi_rebin[0][8][j] = (TH1D*) signal_dPhi_rebin[0][0][j]->Clone((TString)("Combined_dPhi_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
 
   signal_dEta_syst[0][8][j] = (TH1D*) signal_dEta_syst[0][0][j]->Clone((TString)("Combined_dEta_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
   signal_dPhi_syst[0][8][j] = (TH1D*) signal_dPhi_syst[0][0][j]->Clone((TString)("Combined_dPhi_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
 
   if(j==0){
     signal_dEta_rebin[1][8][j] = (TH1D*) signal_dEta_rebin[1][0][j]->Clone((TString)("Combined_dEta_pp_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
     signal_dPhi_rebin[1][8][j] = (TH1D*) signal_dPhi_rebin[1][0][j]->Clone((TString)("Combined_dPhi_pp_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
 
     signal_dEta_syst[1][8][j] = (TH1D*) signal_dEta_syst[1][0][j]->Clone((TString)("Combined_dEta_pp_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
     signal_dPhi_syst[1][8][j] = (TH1D*) signal_dPhi_syst[1][0][j]->Clone((TString)("Combined_dPhi_pp_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
   }


    signal_dEta_PbPb_pp[0][8][j] = (TH1D*) signal_dEta_PbPb_pp[0][0][j]->Clone((TString)("Combined_dEta_PbPb_pp_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
    signal_dPhi_PbPb_pp[0][8][j] = (TH1D*) signal_dPhi_PbPb_pp[0][0][j]->Clone((TString)("Combined_dPhi_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
 
   signal_dEta_PbPb_pp_syst[0][8][j] = (TH1D*) signal_dEta_PbPb_pp_syst[0][0][j]->Clone((TString)("Combined_dEta_PbPb_pp_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
   signal_dPhi_PbPb_pp_syst[0][8][j] = (TH1D*) signal_dPhi_PbPb_pp_syst[0][0][j]->Clone((TString)("Combined_dPhi_PbPb_pp_" + CBin_strs[j] + "_" + CBin_strs[j+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
 

    for(int k = 1; k<8; k++){
 

      signal_dEta_rebin[0][8][j]->Add(signal_dEta_rebin[0][k][j]);
      signal_dPhi_rebin[0][8][j]->Add(signal_dPhi_rebin[0][k][j]);

      signal_dEta_syst[0][8][j]->Add(signal_dEta_syst[0][k][j]);
      signal_dPhi_syst[0][8][j]->Add(signal_dPhi_syst[0][k][j]);


      signal_dEta_PbPb_pp[0][8][j]->Add(signal_dEta_PbPb_pp[0][k][j]);
      signal_dPhi_PbPb_pp[0][8][j]->Add(signal_dPhi_PbPb_pp[0][k][j]);

      signal_dEta_PbPb_pp_syst[0][8][j]->Add(signal_dEta_PbPb_pp_syst[0][k][j]);
      signal_dPhi_PbPb_pp_syst[0][8][j]->Add(signal_dPhi_PbPb_pp_syst[0][k][j]);
     

     
      if(j==0){
	signal_dEta_syst[1][8][j]->Add(signal_dEta_syst[1][k][j]);
	signal_dPhi_syst[1][8][j]->Add(signal_dPhi_syst[1][k][j]);

	signal_dEta_rebin[1][8][j]->Add(signal_dEta_rebin[1][k][j]);
	signal_dPhi_rebin[1][8][j]->Add(signal_dPhi_rebin[1][k][j]);

      }


    }
  

 
    signal_dEta_syst[0][8][j]->SetFillColor(kBlack);
    signal_dPhi_syst[0][8][j]->SetFillColor(kBlack);

    signal_dEta_syst[0][8][j]->SetMarkerColor(kWhite);
    signal_dPhi_syst[0][8][j]->SetMarkerColor(kWhite);

    

    signal_dEta_syst[0][8][j]->SetFillStyle(3004);
    signal_dPhi_syst[0][8][j]->SetFillStyle(3004);

    signal_dEta_rebin[0][8][j]->SetMarkerStyle(24);
    signal_dPhi_rebin[0][8][j]->SetMarkerStyle(24);



 
    signal_dEta_PbPb_pp_syst[0][8][j]->SetFillColor(kBlack);
    signal_dPhi_PbPb_pp_syst[0][8][j]->SetFillColor(kBlack);

    signal_dEta_PbPb_pp_syst[0][8][j]->SetMarkerColor(kWhite);
    signal_dPhi_PbPb_pp_syst[0][8][j]->SetMarkerColor(kWhite);

    

    signal_dEta_PbPb_pp_syst[0][8][j]->SetFillStyle(3004);
    signal_dPhi_PbPb_pp_syst[0][8][j]->SetFillStyle(3004);

    signal_dEta_PbPb_pp[0][8][j]->SetMarkerStyle(24);
    signal_dPhi_PbPb_pp[0][8][j]->SetMarkerStyle(24);




     
    if(j==0){
      signal_dEta_syst[1][8][j]->SetFillColor(kBlack);
      signal_dPhi_syst[1][8][j]->SetFillColor(kBlack);
 
      signal_dEta_syst[1][8][j]->SetMarkerColor(kWhite);
      signal_dPhi_syst[1][8][j]->SetMarkerColor(kWhite);

 

      signal_dEta_syst[1][8][j]->SetFillStyle(3004);
      signal_dPhi_syst[1][8][j]->SetFillStyle(3004);

      signal_dEta_rebin[1][8][j]->SetMarkerStyle(24);
      signal_dPhi_rebin[1][8][j]->SetMarkerStyle(24);


    }




    signal_dEta_stack[0][j] = new THStack((TString)("Signal_dEta_Diff_Stack_"+CBin_strs[j]+CBin_strs[j+1]),"");
    signal_dPhi_stack[0][j] = new THStack((TString)("Signal_dPhi_Diff_Stack_"+CBin_strs[j]+CBin_strs[j+1]),"");

    signal_dEta_diff_stack_up[0][j] = new THStack((TString)("Signal_dEta_Diff_Stack_Up"+CBin_strs[j]+CBin_strs[j+1]),"");
    signal_dPhi_diff_stack_up[0][j] = new THStack((TString)("Signal_dPhi_Diff_Stack_Up"+CBin_strs[j]+CBin_strs[j+1]),"");

    signal_dEta_diff_stack_down[0][j] = new THStack((TString)("Signal_dEta_Diff_Stack_Down"+CBin_strs[j]+CBin_strs[j+1]),"");
    signal_dPhi_diff_stack_down[0][j] = new THStack((TString)("Signal_dPhi_Diff_Stack_Down"+CBin_strs[j]+CBin_strs[j+1]),"");


    if(j==0){
      signal_dEta_stack[1][j] = new THStack((TString)("Signal_dEta_Stack_pp"),"");
      signal_dPhi_stack[1][j] = new THStack((TString)("Signal_dPhi_Stack_pp"),"");
      
    }
 
    


    for(int i = 0; i < nTrkPtBins - 1; i++){
      signal_dEta_rebin[0][i][j]->SetMarkerSize(0.);
      signal_dPhi_rebin[0][i][j]->SetMarkerSize(0.);
      
      signal_dEta_rebin[1][i][0]->SetMarkerSize(0.);
      signal_dPhi_rebin[1][i][0]->SetMarkerSize(0.);
   
      signal_dEta_PbPb_pp_up[0][i][j]->SetMarkerSize(0.);
      signal_dPhi_PbPb_pp_up[0][i][j]->SetMarkerSize(0.);
    
      switch(i){
      case 0:
	signal_dEta_rebin[0][i][j]->SetFillColor(kBlue-9);
	signal_dPhi_rebin[0][i][j]->SetFillColor(kBlue-9);
      	signal_dEta_rebin[1][i][0]->SetFillColor(kBlue-9);
	signal_dPhi_rebin[1][i][0]->SetFillColor(kBlue-9);

	signal_dEta_PbPb_pp_up[0][i][j]->SetFillColor(kBlue-9);
	signal_dPhi_PbPb_pp_up[0][i][j]->SetFillColor(kBlue-9);
   
	signal_dEta_PbPb_pp_down[0][i][j]->SetFillColor(kBlue-9);
	signal_dPhi_PbPb_pp_down[0][i][j]->SetFillColor(kBlue-9);
    
	break;
      case 1:
	signal_dEta_rebin[0][i][j]->SetFillColor(kYellow-9);
	signal_dPhi_rebin[0][i][j]->SetFillColor(kYellow-9);
     	signal_dEta_rebin[1][i][0]->SetFillColor(kYellow-9);
	signal_dPhi_rebin[1][i][0]->SetFillColor(kYellow-9);
   
	signal_dEta_PbPb_pp_up[0][i][j]->SetFillColor(kYellow-9);
	signal_dPhi_PbPb_pp_up[0][i][j]->SetFillColor(kYellow-9);
      
	signal_dEta_PbPb_pp_down[0][i][j]->SetFillColor(kYellow-9);
	signal_dPhi_PbPb_pp_down[0][i][j]->SetFillColor(kYellow-9);
   
	break;
      case 2:
	signal_dEta_rebin[0][i][j]->SetFillColor(kOrange+1);
	signal_dPhi_rebin[0][i][j]->SetFillColor(kOrange+1);
      	signal_dEta_rebin[1][i][0]->SetFillColor(kOrange+1);
	signal_dPhi_rebin[1][i][0]->SetFillColor(kOrange+1);
   
	signal_dEta_PbPb_pp_up[0][i][j]->SetFillColor(kOrange+1);
	signal_dPhi_PbPb_pp_up[0][i][j]->SetFillColor(kOrange+1);

	signal_dEta_PbPb_pp_down[0][i][j]->SetFillColor(kOrange+1);
	signal_dPhi_PbPb_pp_down[0][i][j]->SetFillColor(kOrange+1);

      	
	break;
      case 3:
	signal_dEta_rebin[0][i][j]->SetFillColor(kViolet-5);
	signal_dPhi_rebin[0][i][j]->SetFillColor(kViolet-5);
	signal_dEta_rebin[1][i][0]->SetFillColor(kViolet-5);
	signal_dPhi_rebin[1][i][0]->SetFillColor(kViolet-5);

	signal_dEta_PbPb_pp_up[0][i][j]->SetFillColor(kViolet-5);
	signal_dPhi_PbPb_pp_up[0][i][j]->SetFillColor(kViolet-5);

	signal_dEta_PbPb_pp_down[0][i][j]->SetFillColor(kViolet-5);
	signal_dPhi_PbPb_pp_down[0][i][j]->SetFillColor(kViolet-5);

   	break;
      case 4:
	signal_dEta_rebin[0][i][j]->SetFillColor(kGreen+3);
	signal_dPhi_rebin[0][i][j]->SetFillColor(kGreen+3);
	signal_dEta_rebin[1][i][0]->SetFillColor(kGreen+3);
	signal_dPhi_rebin[1][i][0]->SetFillColor(kGreen+3);

	signal_dEta_PbPb_pp_up[0][i][j]->SetFillColor(kGreen+3);
	signal_dPhi_PbPb_pp_up[0][i][j]->SetFillColor(kGreen+3);
    
	signal_dEta_PbPb_pp_down[0][i][j]->SetFillColor(kGreen+3);
	signal_dPhi_PbPb_pp_down[0][i][j]->SetFillColor(kGreen+3);
    
 	break;
  
      case 5:
	signal_dEta_rebin[0][i][j]->SetFillColor(kRed);
	signal_dPhi_rebin[0][i][j]->SetFillColor(kRed);
      	signal_dEta_rebin[1][i][0]->SetFillColor(kRed);
	signal_dPhi_rebin[1][i][0]->SetFillColor(kRed);
 
	signal_dEta_PbPb_pp_up[0][i][j]->SetFillColor(kRed);
	signal_dPhi_PbPb_pp_up[0][i][j]->SetFillColor(kRed);

 	signal_dEta_PbPb_pp_down[0][i][j]->SetFillColor(kRed);
	signal_dPhi_PbPb_pp_down[0][i][j]->SetFillColor(kRed);
    
  	break;
      case 6:
	signal_dEta_rebin[0][i][j]->SetFillColor(kRed+1);
	signal_dPhi_rebin[0][i][j]->SetFillColor(kRed+1);
   	signal_dEta_rebin[1][i][0]->SetFillColor(kRed+1);
	signal_dPhi_rebin[1][i][0]->SetFillColor(kRed+1);
   
	signal_dEta_PbPb_pp_up[0][i][j]->SetFillColor(kRed+1);
	signal_dPhi_PbPb_pp_up[0][i][j]->SetFillColor(kRed+1);

	signal_dEta_PbPb_pp_down[0][i][j]->SetFillColor(kRed+1);
	signal_dPhi_PbPb_pp_down[0][i][j]->SetFillColor(kRed+1);

	break;
	
      case 7:
	signal_dEta_rebin[0][i][j]->SetFillColor(kRed+2);
	signal_dPhi_rebin[0][i][j]->SetFillColor(kRed+2);
	signal_dEta_rebin[1][i][0]->SetFillColor(kRed+2);
	signal_dPhi_rebin[1][i][0]->SetFillColor(kRed+2);
   
	signal_dEta_PbPb_pp_up[0][i][j]->SetFillColor(kRed+2);
	signal_dPhi_PbPb_pp_up[0][i][j]->SetFillColor(kRed+2);

	signal_dEta_PbPb_pp_down[0][i][j]->SetFillColor(kRed+2);
	signal_dPhi_PbPb_pp_down[0][i][j]->SetFillColor(kRed+2);

	break;

      default:
	break;
      }
    }
    

    cout<<"here"<<endl;
      


    for(int i = 0; i < nTrkPtBins- 1; i++){

      for(int k = 0 ; k< signal_dEta_rebin[0][i][j]->GetNbinsX()+1; k++){
	signal_dEta_rebin[0][nTrkPtBins- 2 - i][j]->SetBinError(k,0.);
	signal_dEta_PbPb_pp_up[0][nTrkPtBins- 2 - i][j]->SetBinError(k,0.);
	signal_dEta_PbPb_pp_down[0][nTrkPtBins- 2 - i][j]->SetBinError(k,0.);
	if(	signal_dEta_PbPb_pp_up[0][nTrkPtBins- 2 - i][j]->GetBinContent(k)<0)	signal_dEta_PbPb_pp_up[0][nTrkPtBins- 2 - i][j]->SetBinContent(k,0.);
	if(	signal_dEta_PbPb_pp_down[0][nTrkPtBins- 2 - i][j]->GetBinContent(k)>0)	signal_dEta_PbPb_pp_down[0][nTrkPtBins- 2 - i][j]->SetBinContent(k,0.);
      }

      signal_dEta_stack[0][j]->Add(signal_dEta_rebin[0][nTrkPtBins- 2 - i][j]);
      signal_dEta_diff_stack_up[0][j]->Add(signal_dEta_PbPb_pp_up[0][nTrkPtBins- 2 - i][j]);
      signal_dEta_diff_stack_down[0][j]->Add(signal_dEta_PbPb_pp_down[0][nTrkPtBins- 2 - i][j]);

      for(int k = 0 ; k< signal_dPhi_rebin[0][i][j]->GetNbinsX()+1; k++){
	signal_dPhi_rebin[0][nTrkPtBins- 2 - i][j]->SetBinError(k,0.);
	signal_dPhi_PbPb_pp_up[0][nTrkPtBins- 2 - i][j]->SetBinError(k,0.);
	signal_dPhi_PbPb_pp_down[0][nTrkPtBins- 2 - i][j]->SetBinError(k,0.);
	if(	signal_dPhi_PbPb_pp_up[0][nTrkPtBins- 2 - i][j]->GetBinContent(k)<0)	signal_dPhi_PbPb_pp_up[0][nTrkPtBins- 2 - i][j]->SetBinContent(k,0.);
	if(	signal_dPhi_PbPb_pp_down[0][nTrkPtBins- 2 - i][j]->GetBinContent(k)>0)	signal_dPhi_PbPb_pp_down[0][nTrkPtBins- 2 - i][j]->SetBinContent(k,0.);

      }

      signal_dPhi_stack[0][j]->Add(signal_dPhi_rebin[0][nTrkPtBins- 2 - i][j]);
      signal_dPhi_diff_stack_up[0][j]->Add(signal_dPhi_PbPb_pp_up[0][nTrkPtBins- 2 - i][j]);
      signal_dPhi_diff_stack_down[0][j]->Add(signal_dPhi_PbPb_pp_down[0][nTrkPtBins- 2 - i][j]);
      
      if(j==0){
	for(int k = 0 ; k< signal_dEta_rebin[1][i][j]->GetNbinsX()+1; k++){
	  signal_dEta_rebin[1][nTrkPtBins- 2 - i][j]->SetBinError(k,0.);
	}

	signal_dEta_stack[1][j]->Add(signal_dEta_rebin[1][nTrkPtBins- 2 - i][j]);

	for(int k = 0 ; k< signal_dPhi_rebin[1][i][j]->GetNbinsX()+1; k++){
	  signal_dPhi_rebin[1][nTrkPtBins- 2 - i][j]->SetBinError(k,0.);
	}

	signal_dPhi_stack[1][j]->Add(signal_dPhi_rebin[1][nTrkPtBins- 2 - i][j]);

      }
    }
  }

  cout<<"here"<<endl;

  TCanvas *c_stacked_eta = new TCanvas("c_stacked_eta","",10,10,2500,1000);
  c_stacked_eta->Divide(5,2,0.,0.);

  c_stacked_eta->cd(1);


  signal_dEta_stack[1][0]->SetMinimum(stacked_min);
    signal_dEta_stack[1][0]->SetMaximum(stacked_max);

    signal_dEta_stack[1][0]->Draw();
    signal_dEta_stack[1][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
    signal_dEta_stack[1][0]->GetXaxis()->SetTitleSize(0.06);
    signal_dEta_stack[1][0]->GetXaxis()->SetLabelSize(0.06);
    signal_dEta_stack[1][0]->GetXaxis()->CenterTitle();
    signal_dEta_stack[1][0]->GetXaxis()->SetTitle("#Delta#eta");
    
    signal_dEta_stack[1][0]->GetYaxis()->SetTitleSize(0.06);
    signal_dEta_stack[1][0]->GetYaxis()->SetLabelSize(0.06);
    signal_dEta_stack[1][0]->GetYaxis()->SetTitleOffset(1.2);
    signal_dEta_stack[1][0]->GetYaxis()->SetTitle("Y = #frac{1}{N_{jets}} #frac{dN}{d#Delta#eta}");

       labels = new TPaveText(0.18,0.75,0.45,0.95,"NDC");
    
    labels->SetName("labels");
    labels->SetFillColor(0);
    labels->SetLineColor(0);
    labels->SetTextAlign(11);
    labels->AddText("pp reference");
    labels->SetTextSize(0.06);
    labels->Draw("same");

    signal_dEta_syst[1][8][0]->Draw("same e2");
    signal_dEta_rebin[1][8][0]->Draw("same");
  
  for(int j = 0; j<nCBins; j++){
    c_stacked_eta->cd(5-j);
    signal_dEta_stack[0][j]->SetMinimum(stacked_min);
    signal_dEta_stack[0][j]->SetMaximum(stacked_max);


    signal_dEta_stack[0][j]->Draw();
    signal_dEta_stack[0][j]->GetXaxis()->SetRangeUser(-1.5,1.5);   
   signal_dEta_stack[0][j]->GetXaxis()->SetTitle("#Delta#eta");
   
    signal_dEta_stack[0][j]->GetYaxis()->SetLabelSize(0.0);
  
   
    labels = new TPaveText(0.05,0.75,0.45,0.95,"NDC");
    
    labels->SetName("labels");
    labels->SetFillColor(0);
    labels->SetLineColor(0);
    labels->SetTextAlign(11);
    labels->AddText((TString)("PbPb "+CBin_labels[j]));
    labels->SetTextSize(0.06);
    labels->Draw("same");

    signal_dEta_syst[0][8][j]->Draw("same e2");
    signal_dEta_rebin[0][8][j]->Draw("same");

    
    c_stacked_eta->cd(10-j);

    signal_dEta_diff_stack_up[0][j]->SetMinimum(stacked_min_diff);
    signal_dEta_diff_stack_up[0][j]->SetMaximum(stacked_max_diff);

    signal_dEta_diff_stack_up[0][j]->Draw();
    signal_dEta_diff_stack_up[0][j]->GetXaxis()->SetRangeUser(-1.5,1.5);
    signal_dEta_diff_stack_up[0][j]->GetXaxis()->SetTitleSize(0.06);
    signal_dEta_diff_stack_up[0][j]->GetXaxis()->SetLabelSize(0.06);
    signal_dEta_diff_stack_up[0][j]->GetXaxis()->CenterTitle();
    signal_dEta_diff_stack_up[0][j]->GetXaxis()->SetTitle("#Delta#eta");
    
    if(j==3){
      signal_dEta_diff_stack_up[0][j]->GetYaxis()->SetTitleSize(0.06);
      signal_dEta_diff_stack_up[0][j]->GetYaxis()->SetLabelSize(0.06);
      signal_dEta_diff_stack_up[0][j]->GetYaxis()->SetTitleOffset(1.5);
      signal_dEta_diff_stack_up[0][j]->GetYaxis()->SetTitle("Y_{PbPb} - Y_{pp}");
    }else{
      signal_dEta_diff_stack_up[0][j]->GetYaxis()->SetLabelSize(0.0);
    }
    signal_dEta_diff_stack_up[0][j]->Draw();
    signal_dEta_diff_stack_down[0][j]->Draw("same");

  labels = new TPaveText(0.05,0.75,0.45,0.95,"NDC");
    
    labels->SetName("labels");
    labels->SetFillColor(0);
    labels->SetLineColor(0);
    labels->SetTextAlign(11);
    labels->AddText((TString)("PbPb ("+CBin_labels[j]+") minus pp"));
    labels->SetTextSize(0.06);
    labels->Draw("same");

    cout<<"about to draw syst"<<endl;


    signal_dEta_PbPb_pp_syst[0][8][j]->Draw("same e2");
    signal_dEta_PbPb_pp[0][8][j]->Draw("same");


  }

  c_stacked_eta->cd(6);
  TLegend *legend = new TLegend(0.1,0.05,0.9,.9);
  legend->AddEntry(signal_dEta_rebin[1][0][0],"0.7 < p_{T}^{assoc.}< 1 GeV","f");
  legend->AddEntry(signal_dEta_rebin[1][1][0],"1 < p_{T}^{assoc.}< 2 GeV","f");
  legend->AddEntry(signal_dEta_rebin[1][2][0],"2 < p_{T}^{assoc.}< 3 GeV","f");
  legend->AddEntry(signal_dEta_rebin[1][3][0],"3 < p_{T}^{assoc.}< 4 GeV","f");
  legend->AddEntry(signal_dEta_rebin[1][4][0],"4 < p_{T}^{assoc.}< 8 GeV","f");
  legend->AddEntry(signal_dEta_rebin[1][5][0],"8 < p_{T}^{assoc.}< 12 GeV","f");
  legend->AddEntry(signal_dEta_rebin[1][6][0],"12 < p_{T}^{assoc.}< 16 GeV","f");
  legend->AddEntry(signal_dEta_rebin[1][7][0],"16 < p_{T}^{assoc.}< 20 GeV","f");
  
  legend->SetTextSize(0.06);
  legend->SetLineColor(kWhite);
  legend->Draw();

  TGaxis *dummy_axis_diff = new TGaxis(1.,0.13,1.0,.975,stacked_min_diff, stacked_max_diff);

  dummy_axis_diff->ImportAxisAttributes( signal_dEta_stack[1][0]->GetYaxis());
  dummy_axis_diff->SetTitleOffset(1.);
  dummy_axis_diff->SetTickSize(0.);
  dummy_axis_diff->CenterTitle();
  dummy_axis_diff->SetTitleSize(0.08);
  dummy_axis_diff->SetLabelSize(0.06);
  dummy_axis_diff->SetTitle("Y_{PbPb} - Y_{pp}");
  dummy_axis_diff->Draw();
  
  c_stacked_eta->cd(0);

  TLatex *type_tex;
  type_tex = new TLatex(0.05,0.92,"Particle Yield by #Delta#eta");
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


  
  c_stacked_eta->SaveAs("Yield_dEta_Stacked.png");
  c_stacked_eta->SaveAs("Yield_dEta_Stacked.pdf");
 



  TCanvas *c_stacked_phi = new TCanvas("c_stacked_phi","",10,10,2500,1000);
  c_stacked_phi->Divide(5,2,0.,0.);

  c_stacked_phi->cd(1);



  signal_dPhi_stack[1][0]->SetMinimum(stacked_min);
    signal_dPhi_stack[1][0]->SetMaximum(stacked_max);

    signal_dPhi_stack[1][0]->Draw();
    signal_dPhi_stack[1][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
    signal_dPhi_stack[1][0]->GetXaxis()->SetTitleSize(0.06);
    signal_dPhi_stack[1][0]->GetXaxis()->SetLabelSize(0.06);
    signal_dPhi_stack[1][0]->GetXaxis()->CenterTitle();
    signal_dPhi_stack[1][0]->GetXaxis()->SetTitle("#Delta#phi");
    
    signal_dPhi_stack[1][0]->GetYaxis()->SetTitleSize(0.06);
    signal_dPhi_stack[1][0]->GetYaxis()->SetLabelSize(0.06);
    signal_dPhi_stack[1][0]->GetYaxis()->SetTitleOffset(1.2);
    signal_dPhi_stack[1][0]->GetYaxis()->SetTitle("Y = #frac{1}{N_{jets}} #frac{dN}{d#Delta#phi}");


    labels = new TPaveText(0.18,0.75,0.45,0.95,"NDC");
    
    labels->SetName("labels");
    labels->SetFillColor(0);
    labels->SetLineColor(0);
    labels->SetTextAlign(11);
    labels->AddText("pp reference");
    labels->SetTextSize(0.06);
    labels->Draw("same");

    signal_dPhi_syst[1][8][0]->Draw("same e2");
    signal_dPhi_rebin[1][8][0]->Draw("same");
  

  for(int j = 0; j<nCBins; j++){
    c_stacked_phi->cd(5-j);
    signal_dPhi_stack[0][j]->SetMinimum(stacked_min);
    signal_dPhi_stack[0][j]->SetMaximum(stacked_max);

  
    signal_dPhi_stack[0][j]->Draw();
    signal_dPhi_stack[0][j]->GetXaxis()->SetRangeUser(-1.5,1.5); 
   signal_dPhi_stack[0][j]->GetXaxis()->SetTitle("#Delta#phi");
   
    signal_dPhi_stack[0][j]->GetYaxis()->SetLabelSize(0.0);

    labels = new TPaveText(0.05,0.75,0.45,0.95,"NDC");
    
    labels->SetName("labels");
    labels->SetFillColor(0);
    labels->SetLineColor(0);
    labels->SetTextAlign(11);
    labels->AddText((TString)("PbPb "+CBin_labels[j]));
    labels->SetTextSize(0.06);
    labels->Draw("same");

    signal_dPhi_syst[0][8][j]->Draw("same e2");
    signal_dPhi_rebin[0][8][j]->Draw("same");
    
    c_stacked_phi->cd(10-j);

    signal_dPhi_diff_stack_up[0][j]->SetMinimum(stacked_min_diff);
    signal_dPhi_diff_stack_up[0][j]->SetMaximum(stacked_max_diff);

    signal_dPhi_diff_stack_up[0][j]->Draw();
    signal_dPhi_diff_stack_up[0][j]->GetXaxis()->SetRangeUser(-1.5,1.5);
    signal_dPhi_diff_stack_up[0][j]->GetXaxis()->SetTitleSize(0.06);
    signal_dPhi_diff_stack_up[0][j]->GetXaxis()->SetLabelSize(0.06);
    signal_dPhi_diff_stack_up[0][j]->GetXaxis()->CenterTitle();
    signal_dPhi_diff_stack_up[0][j]->GetXaxis()->SetTitle("#Delta#phi");
    
    if(j==3){
      signal_dPhi_diff_stack_up[0][j]->GetYaxis()->SetTitleSize(0.06);
      signal_dPhi_diff_stack_up[0][j]->GetYaxis()->SetLabelSize(0.06);
      signal_dPhi_diff_stack_up[0][j]->GetYaxis()->SetTitleOffset(1.5);
      signal_dPhi_diff_stack_up[0][j]->GetYaxis()->SetTitle("Y_{PbPb} - Y_{pp}");
    }else{
      signal_dPhi_diff_stack_up[0][j]->GetYaxis()->SetLabelSize(0.0);
    }
    signal_dPhi_diff_stack_up[0][j]->Draw();
    signal_dPhi_diff_stack_down[0][j]->Draw("same");

    labels = new TPaveText(0.05,0.75,0.45,0.95,"NDC");
    
    labels->SetName("labels");
    labels->SetFillColor(0);
    labels->SetLineColor(0);
    labels->SetTextAlign(11);
    labels->AddText((TString)("PbPb ("+CBin_labels[j]+") minus pp"));
    labels->SetTextSize(0.06);
    labels->Draw("same");


    signal_dPhi_PbPb_pp_syst[0][8][j]->Draw("same e2");
    signal_dPhi_PbPb_pp[0][8][j]->Draw("same");


  

  }

  c_stacked_phi->cd(6);


  legend->Draw();

  dummy_axis_diff->Draw();

  c_stacked_phi->cd(0);

  type_tex = new TLatex(0.05,0.92,"Particle Yield by #Delta#phi");
  type_tex->SetTextSize(0.035);
  type_tex->SetLineColor(kWhite);
  type_tex->SetNDC();
  type_tex->Draw();
   
  luminosity_tex_pp->Draw();
  luminosity_tex_PbPb->Draw();
  jet_reco_tex->Draw();
    
  c_stacked_phi->SaveAs("Yield_dPhi_Stacked.png");
  c_stacked_phi->SaveAs("Yield_dPhi_Stacked.pdf");



  for(int j = 0; j < nCBins; j++){

    for(int g = 0; g< 2; g++){

      if(g==0){

    	Integral_phi_Pt[g][j] = new TH1D((TString)("Integral_Phi_PbPb"+CBin_strs[j]+"_"+CBin_strs[j+1]),"",9,pTbins);
   	Integral_eta_Pt[g][j] = new TH1D((TString)("Integral_Eta_PbPb"+CBin_strs[j]+"_"+CBin_strs[j+1]),"",9,pTbins);
         
      }else{

	if(j > 0 )continue;

	Integral_phi_Pt[g][j] = new TH1D("Integral_Phi_pp","",9,pTbins);
	Integral_eta_Pt[g][j] = new TH1D("Integral_Eta_pp","",9,pTbins);
      }

     
      for(int i = 0; i<nTrkPtBins; i++){
      
	cout<< integral_phi[g][i][j]<<" "<<integral_phi_err[g][i][j]<<" "<<integral_eta[g][i][j]<<" "<<integral_eta_err[g][i][j]<<endl;

	Integral_phi_Pt[g][j]->SetBinContent(i+2,integral_phi[g][i][j]);
	Integral_phi_Pt[g][j]->SetBinError(i+2,integral_phi_err[g][i][j]);

	Integral_eta_Pt[g][j]->SetBinContent(i+2,integral_eta[g][i][j]);
	Integral_eta_Pt[g][j]->SetBinError(i+2,integral_eta_err[g][i][j]);

      
      }

      if(g==0){
	Integral_phi_Pt[g][j]->SetMarkerStyle(10);
	Integral_eta_Pt[g][j]->SetMarkerStyle(10);
      }else{
	Integral_phi_Pt[g][j]->SetMarkerStyle(4);
	Integral_eta_Pt[g][j]->SetMarkerStyle(4);

      }
      Integral_phi_Pt[g][j]->SetMarkerSize(2);
      Integral_eta_Pt[g][j]->SetMarkerSize(2);

      Integral_phi_Pt[g][j]->SetLineColor(kBlack);
      Integral_phi_Pt[g][j]->SetMarkerColor(kBlack);
  
      Integral_eta_Pt[g][j]->SetLineColor(kBlack);
      Integral_eta_Pt[g][j]->SetMarkerColor(kBlack);


    }
  }

  TCanvas *c_integral = new TCanvas("Integral_Canvas","",10,10,2000,1000);
  c_integral->Divide(4,2,0,0);

  for(int j = 0; j < nCBins; j++){
    c_integral->cd(4-j);


    Integral_phi_Pt[0][j]->GetXaxis()->SetRangeUser(0.01,19.5);
  
    Integral_phi_Pt[0][j]->SetMaximum(14.);
    Integral_phi_Pt[0][j]->SetMinimum(-3.);

    Integral_phi_Pt[0][j]->GetXaxis()->SetLabelSize(0.06);
    Integral_phi_Pt[0][j]->GetXaxis()->SetTitleSize(0.06);
    Integral_phi_Pt[0][j]->GetXaxis()->SetTitle("p_{T}^{assoc}");

    if(j==3){
    Integral_phi_Pt[0][j]->GetYaxis()->SetLabelSize(0.06);
    Integral_phi_Pt[0][j]->GetYaxis()->SetTitleSize(0.06);
    Integral_phi_Pt[0][j]->GetYaxis()->SetTitle("Tracks per Jet");
    }else{
      Integral_phi_Pt[0][j]->GetYaxis()->SetLabelSize(0.);
    }
   
    Integral_phi_Pt[0][j]->Draw();
    Integral_phi_Pt[1][0]->Draw("same");
    Integral_eta_Pt[0][j]->Draw("same");
    Integral_eta_Pt[1][0]->Draw("same");

   
    Integral_phi_Pt[0][j]->Draw();
    Integral_phi_Pt[1][0]->Draw("same");
    Integral_eta_Pt[0][j]->Draw("same");
    Integral_eta_Pt[1][0]->Draw("same");

    if(j==3){
      TLegend *int_legend = new TLegend(0.18,0.6,0.9,0.85);
      int_legend->AddEntry( Integral_phi_Pt[0][j],"PbPb integral |#Delta#eta|<1.0, |#Delta#phi|<1.0");
      int_legend->AddEntry( Integral_phi_Pt[1][0],"pp integral |#Delta#eta|<1.0, |#Delta#phi|<1.0");
      int_legend->SetLineColor(kWhite);
      int_legend->SetTextSize(0.06);
      int_legend->Draw();

      labels = new TPaveText(0.18,0.85,0.45,0.95,"NDC");
    
      labels->SetName("labels");
      labels->SetFillColor(0);
      labels->SetLineColor(0);
      labels->SetTextAlign(11);
      labels->AddText(CBin_labels[j]);
      labels->SetTextSize(0.06);
      labels->Draw("same");
  
    }else{
      labels = new TPaveText(0.05,0.85,0.45,0.95,"NDC");
    
      labels->SetName("labels");
      labels->SetFillColor(0);
      labels->SetLineColor(0);
      labels->SetTextAlign(11);
      labels->AddText(CBin_labels[j]);
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


    c_integral->cd(8-j);

    Integral_diff_Pt[0][j] = (TH1D*)  Integral_phi_Pt[0][j]->Clone((TString)("Integral_Diff_"+CBin_strs[j]+"_"+CBin_strs[j+1]));
    Integral_diff_Pt[0][j]->Add(Integral_phi_Pt[1][0],-1.);

    Integral_diff_Pt[2][j] = (TH1D*) f_in_ref->Get((TString)("Integrated_Yield_Eta_"+CBin_strs[j]+"_"+CBin_strs[j+1]))->Clone((TString)("Integrated_Yield_Eta_"+CBin_strs[j]+"_"+CBin_strs[j+1]));

    Integral_diff_Pt[2][j]->SetLineColor(kBlue);
    Integral_diff_Pt[2][j]->SetMarkerColor(kBlue);
    Integral_diff_Pt[2][j]->SetMarkerStyle(20);
    Integral_diff_Pt[2][j]->SetMarkerSize(1.5);
  
    Integral_diff_Pt[0][j]->SetMaximum(8.5);
    Integral_diff_Pt[0][j]->SetMinimum(-1.5);
    Integral_diff_Pt[0][j]->Draw();

    Integral_diff_Pt[2][j]->Draw("same");
    
    if(j==3){
      TLegend *int_legend = new TLegend(0.18,0.6,0.9,0.85);
      int_legend->AddEntry( Integral_diff_Pt[0][j],"PbPb - pp 5.02 TeV");
      int_legend->AddEntry( Integral_diff_Pt[2][0],"PbPb - pp 2.76 TeV)");
      int_legend->SetLineColor(kWhite);
      int_legend->SetTextSize(0.06);
      int_legend->Draw();

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
    labels->AddText((TString)("PbPb (" + CBin_labels[j]+") minus pp"));
   
    labels->Draw("same");
  
    int_zero->Draw();

    Integral_phi_Pt[0][j]->Write();
    Integral_phi_Pt[1][0]->Write();
    Integral_eta_Pt[0][j]->Write();
    Integral_eta_Pt[1][0]->Write();
    Integral_diff_Pt[0][j]->Write();


    c_integral->cd(0);
    type_tex = new TLatex(0.02,0.96,"Integrated Particle Yields");
    type_tex->SetTextSize(0.035);
    type_tex->SetLineColor(kWhite);
    type_tex->SetNDC();
    type_tex->Draw();
   
    luminosity_tex_pp->Draw();
    luminosity_tex_PbPb->Draw();
    jet_reco_tex->Draw();

  }


  bg_err_pbpb->Write();
  bg_err_pp->Write();
 
  me_err_pbpb->Write();
  me_err_pp->Write();



  c_integral->SaveAs("Integral_Yield.png");
  c_integral->SaveAs("Integral_Yield.pdf");




  
  return 0;
}
