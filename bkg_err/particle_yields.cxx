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

  TCanvas *c_yields_eta[2];
   
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 9;



  TH2D *result[12][nTrkPtBins][nCBins];

  TH1D *signal_dPhi[12][nTrkPtBins][nCBins];
  TH1D *signal_dPhi_syst[12][nTrkPtBins][nCBins];
  TH1D *signal_dPhi_rebin[12][nTrkPtBins][nCBins];
  TH1D *background_diff_rebin[12][nTrkPtBins][nCBins];
  TH1D *background_syst_rebin[12][nTrkPtBins][nCBins];
  TH1D *background_diff_PbPb_pp[12][nTrkPtBins][nCBins];

  TH1D *signal_dEta[12][nTrkPtBins][nCBins];
  TH1D *signal_dEta_syst[12][nTrkPtBins][nCBins];
  TH1D *signal_dEta_rebin[12][nTrkPtBins][nCBins];

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

  TF1 *right_fit = new TF1("right_fit","x+[0]");
  TF1 *left_fit = new TF1("left_fit","x+[0]");

 
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

  double bin_width_phi, bin_width_eta, diff_max, diff_min, signal_min, signal_max, stacked_max, stacked_min, bc, err;

  float rel_err_pbpb =TMath::Sqrt(0.05*0.05+0.05*0.05+0.04*0.04);
  float rel_err_pp =TMath::Sqrt(0.05*0.05+0.04*0.04+0.03*0.03);

  TF1 *fit_line_left = new TF1("fit_line_left","[0]",-2.5,-1.5);
  TF1 *fit_line_right = new TF1("fit_line_right","[0]",1.5,2.5);

  TCanvas *dummy = new TCanvas("dummy");

  TPaveText *labels;

  float me_err[12][6][5];
  float bg_err[12][6][5];
  float bkg_err[12][6][5];
  TH2D *bg_err_pp = new TH2D("bg_err_pp","",9,0.5,8.5,4,0.5,3.5);

  TH2D *bg_err_pbpb = new TH2D("bg_err_PbPb","",9,0.5,8.5,4,0.5,3.5);

  TH2D *me_err_pp = new TH2D("me_err_pp","",9,0.5,8.5,4,0.5,3.5);

  TH2D *me_err_pbpb = new TH2D("me_err_PbPb","",9,0.5,8.5,4,0.5,3.5);

  for(int g = 0; g<2; g++){
    c_yields_eta[g] = new TCanvas(Form("yields_eta_%d",g),"",10,10,2000,2400);
    c_yields_eta[g]->Divide(5,8,0,0);
  }
  
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

      
	  //	signal_dEta_rebin[g][i][j] = Rebin_dEta(signal_dEta[g][i][j]);
	  signal_dEta_rebin[g][i][j] = (TH1D*)signal_dEta[g][i][j]->Clone(Form("SignalEta%d%d%d",g,i,j));
	  signal_dEta_rebin[g][i][j]->Rebin(10);
	  signal_dEta_rebin[g][i][j]->Scale(1./10);
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
	  
	    TString	jff_name =(TString)("JFF_Residual_Eta_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]);
	    jff_residual_dEta[g][i][j] = (TH1D*)f_jff_hyd->Get(jff_name)->Clone(jff_name);
	    
	    jff_name.ReplaceAll("Eta","Phi");
	    jff_residual_dPhi[g][i][j] = (TH1D*)f_jff_hyd->Get(jff_name)->Clone(jff_name);

	    signal_dEta_rebin[g][i][j]->Add(jff_residual_dEta[g][i][j],-1.);
	    signal_dPhi_rebin[g][i][j]->Add(jff_residual_dPhi[g][i][j],-1.);
	  
	  }else{

	    TString	jff_name =(TString)("JFF_Residual_Eta_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]);
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

	signal_dEta_rebin[g][i][j]->SetAxisRange(-3.001,3.001);

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
	bkg_err[g][i][j] = TMath::Max(fit_line_left->GetParError(0),fit_line_right->GetParError(0));
	  
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

    
	for(int g = 0; g<2; g++){
	  if(g==1&&j>0) continue;

	  float lbin = signal_dEta_rebin[g][i][j]->FindBin(-1.5+.001);
	  float rbin = signal_dEta_rebin[g][i][j]->FindBin(1.5-.001);

	  for(int k = lbin; k<rbin+1; k++){
	    signal_dEta_rebin[g][i][j]->SetBinContent(k,0.);
	    signal_dEta_rebin[g][i][j]->SetBinError(k,0.);
	    signal_dEta_syst[g][i][j]->SetBinContent(k,0.);
	    signal_dEta_syst[g][i][j]->SetBinError(k,0.);
	  }

	
	  switch(i){
	  case 0: 
	    signal_max = .5;
	    signal_min = -.5;
	    break;
	  case 1: 
	    signal_max = .5;
	    signal_min = -.5;
	    break;
	  case 2: 
	    signal_max = .2;
	    signal_min = -.2;
	    break;
	  case 3:
	    signal_max = .2;
	    signal_min = -.2; 
	    break;
	  case 4: 
	    signal_max = .1;
	    signal_min = -.1;
	    break;
	  case 5: 
	    signal_max = .1;
	    signal_min = -.1;
	    break;
	  default: 
	    break;
	  }
            
	  //_____________
	  //  Draw dEta
	  //-------------

	  fit_line_left->SetLineColor(kViolet);
	  fit_line_right->SetLineColor(kRed);


	  if(g==0)	  c_yields_eta[g]->cd(5*(i+1)-j);
	  else 	  c_yields_eta[0]->cd(5*(i)+1);

	  signal_dEta_rebin[g][i][j]->SetMinimum(signal_min);
	  signal_dEta_rebin[g][i][j]->SetMaximum(signal_max);

	  if(g==1){
	    signal_dEta_rebin[g][i][j]->GetYaxis()->SetLabelSize(0.08);
	    signal_dEta_rebin[g][i][j]->GetYaxis()->SetTitleSize(0.08);
	    signal_dEta_rebin[g][i][j]->GetYaxis()->SetTitle("1/N_{evt} 1/d#Delta#eta");
	  }else{
	    signal_dEta_rebin[g][i][j]->GetYaxis()->SetLabelSize(0.0);
	  }

	  if(i==7){
	    signal_dEta_rebin[g][i][j]->GetXaxis()->SetLabelSize(0.08);
	    signal_dEta_rebin[g][i][j]->GetXaxis()->SetTitleSize(0.08);
	    signal_dEta_rebin[g][i][j]->GetXaxis()->SetTitle("#Delta#eta");
	    signal_dEta_rebin[g][i][j]->GetXaxis()->CenterTitle();

	  }
	  // signal_dEta_syst[g][i][j]->SetFillColor(kYellow);
	  signal_dEta_rebin[g][i][j]->Draw("p");
	  //	  signal_dEta_syst[g][i][j]->Draw("same e5");
	
	  TLine *zero = new TLine(-3.,0.,3.,0.);
	  zero->SetLineStyle(2);
	  zero->Draw();
	  



	  if(g==1){
	    labels = new TPaveText(0.18,0.65,0.45,0.95,"NDC");
	    labels->SetName("labels");
	    labels->SetFillColor(0);
	    labels->SetLineColor(0);
	    labels->SetTextAlign(11);
	    labels->AddText("pp reference");
	    labels->AddText(TrkPtBin_labels[i]);
	    labels->AddText(Form("Bkg error: %f",bkg_err[g][i][j]));
	    labels->AddText(Form("ME error: %f",	me_err[g][i][j]));
	    labels->SetTextSize(0.06);
	    labels->Draw("same");

	  }else{
	    labels = new TPaveText(0.05,0.65,0.45,0.95,"NDC");
	    labels->SetName("labels");
	    labels->SetFillColor(0);
	    labels->SetLineColor(0);
	    labels->SetTextAlign(11);
	    labels->AddText(CBin_labels[j]);
	    labels->AddText(TrkPtBin_labels[i]);
	    labels->AddText(Form("Bkg error: %f",bkg_err[g][i][j]));
	    labels->AddText(Form("ME error: %f",	me_err[g][i][j]));
	    labels->SetTextSize(0.06);
	    labels->Draw("same");

	  }  
	
	  if(i==0&&g==1){
	    TLegend *legend = new TLegend(0.5,0.85,0.95,0.95);

	    legend->AddEntry(signal_dEta_rebin[g][i][j],"PbPb");
	    legend->AddEntry(signal_dEta_rebin[g][i][j],"pp");
	
	    legend->SetTextSize(0.06);
	    legend->SetLineColor(0);
	    legend->Draw();
	  }


	  zero->Draw();
	  labels->Draw("same");

	  fit_line_left->Draw("same");
	  fit_line_right->Draw("same");


	}
    
    }
  }

  
  c_yields_eta[0]->SaveAs((TString)("Yields_dEta_Inclusive.png"));
  c_yields_eta[0]->SaveAs((TString)("Yields_dEta_Inclusive.pdf"));
 

  return 0;
}
