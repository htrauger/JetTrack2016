
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
#include "TLegend.h"

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TExec.h"
#include "TLatex.h"

#include "../../HIN-14-016/HIN_14_016_functions.h"

#include <iostream>
#include <vector>
#include <fstream>


using namespace std;

Int_t PAS_plots8(bool draw_ref = kFALSE){

#include "../../HIN-14-016/HIN_14_016_universals.h"

  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.15);
  gStyle->SetPadLeftMargin  (0.16);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);
  gStyle->SetTextFont(43);
  gStyle->SetCanvasBorderMode(0);
 
 

  TFile *fin[12];
  TFile *fin_ref[12];
  TFile *fin2[12];

  TFile *fbkgfit = new TFile("../../HIN-14-016/bg_fits/PbPb_Leading_Yield_and_Bkg.root","READ");
  TFile *fbkgsummed = new TFile("../../HIN-14-016/me_correct/PbPb_Leading_Correlations.root","READ");
 
  TFile *fbkgfit_sub = new TFile("../../HIN-14-016/bg_fits/PbPb_SubLeading_Yield_and_Bkg.root","READ");
  TFile *fbkgsummed_sub = new TFile("../../HIN-14-016/me_correct/PbPb_SubLeading_Correlations.root","READ");
  


  Double_t xAxis[5] = {-100,-50,-30,-10,0}; 
  TH1D* int_cent[12][4];
  TH1D* blank[4];
  TH1D* blank2[4];

  TH1D *check_new_phi_rebin[12][5][4];
  TH1D *check_new_phi_ref[12][5][4];
  TH1D *check_new_phi_syst[12][5][4];
  
  TH1D *check_new_eta_rebin[12][5][4];
  TH1D *check_new_eta_ref[12][5][4];
  TH1D *check_new_eta_syst[12][5][4];

  TH1D *PbPb_pp_phi[12][5][4];
  TH1D *PbPb_pp_phi_syst[12][5][4];

  TH1D *PbPb_pp_eta[12][5][4];
  TH1D *PbPb_pp_eta_syst[12][5][4];

  TH1D *PbPb_pp_phi_ref[12][5][4];
  TH1D *PbPb_pp_eta_ref[12][5][4];

  TH1D *HYDJET_PYTHIA_eta[12][5][4];
  TH1D *HYDJET_PYTHIA_phi[12][5][4];

 
  TF1 *gaus_phi[12][5][4];
  TF1 *gaus_eta[12][5][4];

  
  TH1D *Integral_phi_syst[12][4];
  TH1D *Integral_phi_Pt[12][4];

  TH1D *dummy_pTaxis[4];
 
  TH1D *Integral_eta_Pt[12][4];
  TH1D *Integral_eta_ref_Pt[12][4];
  TGraphAsymmErrors *Integral_eta_syst[12][4];
  TGraphAsymmErrors *Integral_eta_cent[12][4];

  TGraph *Error_up_eta_pT[12][4];
  TGraph *Error_down_eta_pT[12][4];
  TH1D *Error_up_eta_cent[12][4];
  TH1D *Error_down_eta_cent[12][4];

  TLegend *l40,*l41,*l42;

  TString datalabel, in_name, centlabel, pTlabel;
  TString   leadingrefetaname[5][4];

  TCanvas *cFigure1[4],*cFigure2[4],*cFigure3[4],*cFigure4[4],*cFigure5[4],*cFigure6[4],*cFigure7;

  float raw_min,raw_max,mixed_min,mixed_max,yield_min,yield_max,result_min,result_max,  val_l,val_r,err_l,err_r;

 
  TCanvas *cintegral_eta_pT[12];
  TCanvas *cintegral_phi_pT[12];
  TCanvas *cError_up_eta_pT[12];
   
  TCanvas *cintegral_eta_cent[12];
  TCanvas *cintegral_phi_cent[12];
  TCanvas *cError_up_eta_cent[12];
  
  TLine *linePhi,*lineEta,*linePt, *lineCent;

  TLatex *tex24phi,*tex24eta;


  float value, error, evalpt, value2, dx_eta, dx_phi;

  vector<float> pTbin_centers;
  pTbin_centers.push_back(1.5);
  pTbin_centers.push_back(2.5);
  pTbin_centers.push_back(3.5);
  pTbin_centers.push_back(6.0);
  vector<float> pTbin_errors;
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(2.);

  float  dummy_pTbins[] = {0.8,1.0,2.0,3.0,4.0,8.0,8.2 };


  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
 
  for(int g = 0; g<2; g++){
    
   
    TString integral_eta_pT_name = "integral_eta_pT";
    integral_eta_pT_name+=g;

    cintegral_eta_pT[g] = new TCanvas(integral_eta_pT_name,"",10,10,1500,500);
    cintegral_eta_pT[g]->Divide(4,1,0.,0.);


    TString integral_phi_pT_name = "integral_phi_pT";
    integral_phi_pT_name+=g;

    cintegral_phi_pT[g] = new TCanvas(integral_phi_pT_name,"",10,10,1500,500);
    cintegral_phi_pT[g]->Divide(4,1,0.,0.);

    TString integral_eta_cent_name = "integral_eta_cent";
    integral_eta_cent_name+=g;

    cintegral_eta_cent[g] = new TCanvas(integral_eta_cent_name,"",10,10,1600,450);
    cintegral_eta_cent[g]->Divide(4,1,0.0001,0.001);

    TString integral_phi_cent_name = "integral_phi_cent";
    integral_phi_cent_name+=g;

    cintegral_phi_cent[g] = new TCanvas(integral_phi_cent_name,"",10,10,1500,500);
    cintegral_phi_cent[g]->Divide(4,1,0.,0.);


    TString Error_up_eta_pT_name = "Error_up_eta_pT";
    Error_up_eta_pT_name+=g;

    cError_up_eta_pT[g] = new TCanvas(Error_up_eta_pT_name,"",10,10,1500,400);
    cError_up_eta_pT[g]->Divide(4,1,0.,0.);

     
    TString Error_up_eta_cent_name = "Error_up_eta_cent";
    Error_up_eta_cent_name+=g;

    cError_up_eta_cent[g] = new TCanvas(Error_up_eta_cent_name,"",10,10,1500,400);
    cError_up_eta_cent[g]->Divide(4,1,0.,0.);



    //Open files
  
    switch(g){
    case 0:
      fin[g] = new TFile("../../HIN-14-016/study_yield/Inclusive_Data_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../particle_yields/Particle_Yields.root", "READ");
      //  fin[g] = new TFile("../study_yield/Inclusive_Data_NoSpillOver_AllPlots.root", "READ");
      datalabel = "Inclusive";     break;
    case 1:
      fin[g] = new TFile("../../HIN-14-016/study_yield/Inclusive_Data_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../particle_yields/Particle_Yields.root", "READ");
      datalabel = "Inclusive";     break;
     
    default:
      break;
    }

    cout<<g<<endl;

    for(int i=0; i<4; i++){
  
      for (int j=0; j<4; j++){

	

	in_name =	make_name("Result_",g,i,j,0,centlabel,pTlabel);
	

	//-------------------------------------
	//  Get and assign all histograms 
	//-------------------------------------


	TString newchecknamePhi_rebin = in_name;
	TString newchecknameEta_rebin = in_name;


	newchecknameEta_rebin.ReplaceAll("Result", "New_check_Eta");
	newchecknameEta_rebin += "rebin";
	newchecknamePhi_rebin.ReplaceAll("Result", "New_check_Phi");
	newchecknamePhi_rebin += "rebin";
	
	check_new_eta_rebin[g][i][j] = (TH1D*)fin[g]->Get(newchecknameEta_rebin)->Clone((TString)(newchecknameEta_rebin+"new"));
	check_new_eta_rebin[g][i][j]->SetAxisRange(-1.35,1.35,"x");

	
	check_new_phi_rebin[g][i][j] = (TH1D*)fin[g]->Get(newchecknamePhi_rebin)->Clone((TString)(newchecknamePhi_rebin+"new"));
	check_new_phi_rebin[g][i][j]->SetAxisRange(-1.35,1.35,"x");

	if(g==0){
	  TString newchecknameEtasyst = newchecknameEta_rebin;
	  newchecknameEtasyst.ReplaceAll("New_check","Syst_error");
	  newchecknameEtasyst.ReplaceAll("rebin","");

	  check_new_eta_syst[g][i][j] = (TH1D*)fin[g]->Get(newchecknameEtasyst)->Clone((TString)(newchecknameEtasyst+"new"));
	  check_new_eta_syst[g][i][j]->SetAxisRange(-1.4,1.41,"x");

	  TString newchecknamePhisyst = newchecknamePhi_rebin;
	  newchecknamePhisyst.ReplaceAll("New_check","Syst_error");
	  newchecknamePhisyst.ReplaceAll("rebin","");
	
	  check_new_phi_syst[g][i][j] = (TH1D*)fin[g]->Get(newchecknamePhisyst)->Clone((TString)(newchecknamePhisyst+"new"));
	  check_new_phi_syst[g][i][j]->SetAxisRange(-1.4,1.41,"x");
	}
	
	TString newchecknameEta_ref = in_name;
	TString newchecknamePhi_ref = in_name;
	newchecknameEta_ref.ReplaceAll("Result", "Proj_dEta");
	newchecknamePhi_ref.ReplaceAll("Result", "Proj_dPhi");

	newchecknameEta_ref.ReplaceAll("Pt100_Pt300_", "");
	newchecknamePhi_ref.ReplaceAll("Pt100_Pt300_", "");

	newchecknameEta_ref+="_Rebin";
	newchecknamePhi_ref+="_Rebin";
	  

	if(g==1){
	  newchecknameEta_ref.ReplaceAll("Cent0_Cent10_", "");
	  newchecknamePhi_ref.ReplaceAll("Cent0_Cent10_", "");
	  newchecknameEta_ref.ReplaceAll("Cent10_Cent30_", "");
	  newchecknamePhi_ref.ReplaceAll("Cent10_Cent30_", "");
	  newchecknameEta_ref.ReplaceAll("Cent30_Cent50_", "");
	  newchecknamePhi_ref.ReplaceAll("Cent30_Cent50_", "");
	  newchecknameEta_ref.ReplaceAll("Cent50_Cent100_", "");
	  newchecknamePhi_ref.ReplaceAll("Cent50_Cent100_", "");
	  
	}

	cout<<newchecknameEta_ref<<endl;
	check_new_eta_ref[g][i][j] = (TH1D*)fin_ref[g]->Get(newchecknameEta_ref)->Clone(newchecknameEta_ref);
	check_new_phi_ref[g][i][j] = (TH1D*)fin_ref[g]->Get(newchecknamePhi_ref)->Clone(newchecknamePhi_ref);

	check_new_eta_ref[g][i][j]->SetMarkerSize(1);
	check_new_phi_ref[g][i][j]->SetMarkerSize(1);


	if(g==1){
	check_new_phi_ref[g][i][j]->SetMarkerStyle(4);
	check_new_eta_ref[g][i][j]->SetMarkerStyle(4);

	}else{
	check_new_phi_ref[g][i][j]->SetMarkerStyle(20);
	check_new_eta_ref[g][i][j]->SetMarkerStyle(20);

	}

	check_new_eta_ref[g][i][j]->SetMarkerColor(kRed);
	check_new_eta_ref[g][i][j]->SetLineColor(kRed);


	check_new_phi_ref[g][i][j]->SetMarkerColor(kRed);
	check_new_phi_ref[g][i][j]->SetLineColor(kRed);


	cout<<"here"<<endl;


	nbins = check_new_eta_rebin[g][i][j]->GetNbinsX()+1;
	  
	for(int k = 1; k<nbins/2+1;k++){
	    
	  val_l = check_new_eta_rebin[g][i][j]->GetBinContent(k);
	  err_l  = check_new_eta_rebin[g][i][j]->GetBinError(k);
	    
	  val_r = check_new_eta_rebin[g][i][j]->GetBinContent(nbins-k);
	  err_r = check_new_eta_rebin[g][i][j]->GetBinError(nbins-k);
	    
	  cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<TMath::Sqrt(err_l*err_l+err_r*err_r)/2<<endl;

	  check_new_eta_rebin[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	  check_new_eta_rebin[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	  if(g==0||g==2||g==4){
	    check_new_eta_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_eta_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	  }
	  if(k<nbins/2){
	    check_new_eta_rebin[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    check_new_eta_rebin[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	  }
	}




	nbins = check_new_phi_rebin[g][i][j]->GetNbinsX()+1;
	  
	for(int k = 1; k<nbins/2+1;k++){
	    
	  val_l = check_new_phi_rebin[g][i][j]->GetBinContent(k);
	  err_l  = check_new_phi_rebin[g][i][j]->GetBinError(k);
	    
	  val_r = check_new_phi_rebin[g][i][j]->GetBinContent(nbins-k);
	  err_r = check_new_phi_rebin[g][i][j]->GetBinError(nbins-k);
	    
	  cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<TMath::Sqrt(err_l*err_l+err_r*err_r)/2<<endl;

	  check_new_phi_rebin[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	  check_new_phi_rebin[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	  if(g==0||g==2||g==4){
	    check_new_phi_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_phi_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	  }
	  if(k<nbins/2){
	    check_new_phi_rebin[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    check_new_phi_rebin[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	  }
	}


	check_new_phi_rebin[g][i][j]->GetYaxis()->SetTitle("    Y #equiv #frac{1}{N_{jet}} #frac{d^{2}N}{d#Delta#phi dp_{T}}   (GeV^{-1})");
	check_new_eta_rebin[g][i][j]->GetYaxis()->SetTitle("    Y #equiv #frac{1}{N_{jet}} #frac{d^{2}N}{d#Delta#eta dp_{T}}   (GeV^{-1})");


	check_new_phi_rebin[g][i][j]->GetYaxis()->SetTitleOffset(0.8);
	check_new_eta_rebin[g][i][j]->GetYaxis()->SetTitleOffset(0.8);



	
	if(g==1||g==3||g==5){
	  
	  TString PbPbppname_eta =in_name;
	  PbPbppname_eta.ReplaceAll("Result_pp","PbPb_minus_pp_eta");

	  TString PbPbppname_eta_syst =in_name;
	  PbPbppname_eta_syst.ReplaceAll("Result_pp","PbPb_minus_pp_eta_syst");



	  PbPb_pp_eta[g][i][j] = (TH1D*)fin[g-1]->Get(PbPbppname_eta)->Clone(PbPbppname_eta);
	  PbPb_pp_eta_syst[g][i][j] = (TH1D*)fin[g-1]->Get(PbPbppname_eta_syst)->Clone(PbPbppname_eta_syst);
	

	  TString PbPbppname_phi =in_name;
	  PbPbppname_phi.ReplaceAll("Result_pp","PbPb_minus_pp_phi");

	  TString PbPbppname_phi_syst =in_name;
	  PbPbppname_phi_syst.ReplaceAll("Result_pp","PbPb_minus_pp_phi_syst");

	  PbPb_pp_phi[g][i][j] = (TH1D*)fin[g-1]->Get(PbPbppname_phi)->Clone(PbPbppname_phi);
	  PbPb_pp_phi_syst[g][i][j] = (TH1D*)fin[g-1]->Get(PbPbppname_phi_syst)->Clone(PbPbppname_phi_syst);

	  cout<<"here"<<endl;

	  PbPb_pp_eta_ref[g][i][j] = (TH1D*)check_new_eta_ref[g-1][i][j]->Clone((TString)(PbPbppname_eta+"Ref"));
	  PbPb_pp_eta_ref[g][i][j]->Add(check_new_eta_ref[g][i][j],-1.);
	  PbPb_pp_phi_ref[g][i][j] = (TH1D*)check_new_phi_ref[g-1][i][j]->Clone((TString)(PbPbppname_phi+"Ref"));
	  PbPb_pp_phi_ref[g][i][j]->Add(check_new_phi_ref[g][i][j],-1.);
	  PbPb_pp_eta_ref[g][i][j]->SetLineColor(kRed);
	  PbPb_pp_eta_ref[g][i][j]->SetMarkerColor(kRed);

	  PbPb_pp_phi_ref[g][i][j]->SetLineColor(kRed);
	  PbPb_pp_phi_ref[g][i][j]->SetMarkerColor(kRed);

	 
	  nbins = PbPb_pp_eta[g][i][j]->GetNbinsX()+1;
	  
	  for(int k = 1; k<nbins/2+1;k++){
	    val_l = PbPb_pp_eta[g][i][j]->GetBinContent(k);
	 
	    err_l  = PbPb_pp_eta[g][i][j]->GetBinError(k);
	    
	    val_r = PbPb_pp_eta[g][i][j]->GetBinContent(nbins-k);
	    err_r = PbPb_pp_eta[g][i][j]->GetBinError(nbins-k);
	    
	    cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<TMath::Sqrt(err_l*err_l+err_r*err_r)/2<<endl;

	    PbPb_pp_eta[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    PbPb_pp_eta[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	    PbPb_pp_eta_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    PbPb_pp_eta_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);

	    if(k<nbins/2){
	      PbPb_pp_eta[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      PbPb_pp_eta[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    }
	   
	  }

	  nbins = PbPb_pp_phi[g][i][j]->GetNbinsX()+1;
	
	  for(int k = 1; k<nbins/2+1;k++){
	    val_l = PbPb_pp_phi[g][i][j]->GetBinContent(k);
	 
	    err_l  = PbPb_pp_phi[g][i][j]->GetBinError(k);
	    
	    val_r = PbPb_pp_phi[g][i][j]->GetBinContent(nbins-k);
	    err_r = PbPb_pp_phi[g][i][j]->GetBinError(nbins-k);
	    
	    cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<TMath::Sqrt(err_l*err_l+err_r*err_r)/2<<endl;

	    PbPb_pp_phi[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    PbPb_pp_phi[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	    PbPb_pp_phi_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    PbPb_pp_phi_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);

	    if(k<nbins/2){
	      PbPb_pp_phi[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      PbPb_pp_phi[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    }
	   
	  }	

	  PbPb_pp_eta[g][i][j]->SetAxisRange(-1.4,1.41,"x");
	  PbPb_pp_eta_syst[g][i][j]->SetAxisRange(-1.4,1.41,"x");
	  PbPb_pp_phi[g][i][j]->SetAxisRange(-1.4,1.41,"x");
	  PbPb_pp_phi_syst[g][i][j]->SetAxisRange(-1.4,1.41,"x");

	  cout<<"time for integrals"<<endl;
	  /*
	  if(i==0){ //OBVIOUSLY only one set of integrals per centrality class!
	
	
	    TString IntegratedYieldname_eta = in_name;
	    IntegratedYieldname_eta.ReplaceAll("Result_pp","Integrated_Yield_Eta");
	    IntegratedYieldname_eta.ReplaceAll("Pt100_Pt300_","");
	    IntegratedYieldname_eta.ReplaceAll("_TrkPt1_TrkPt2","");
 

	    TString IntegratedYieldname_eta_ref = in_name;
	    IntegratedYieldname_eta_ref.ReplaceAll("Result_pp","Integrated_Yield_Eta_Ref");
	    IntegratedYieldname_eta_ref.ReplaceAll("Pt100_Pt300_","");
	    IntegratedYieldname_eta_ref.ReplaceAll("_TrkPt1_TrkPt2","");
 
	  
	 
	    Integral_eta_ref_Pt[g][j] = (TH1D*)fin_ref[g-1]->Get(IntegratedYieldname_eta)->Clone(IntegratedYieldname_eta_ref);
	    Integral_eta_ref_Pt[g][j]->SetName(IntegratedYieldname_eta_ref);
	    Integral_eta_ref_Pt[g][j]->SetAxisRange(1.01,7.8,"x");



	  
	 
	    Integral_eta_Pt[g][j] = (TH1D*)fin[g-1]->Get(IntegratedYieldname_eta)->Clone(IntegratedYieldname_eta);
	    Integral_eta_Pt[g][j]->SetAxisRange(1.0,7.8,"x");

 
	  
	    TString IntegratedYieldname_eta_syst = IntegratedYieldname_eta;
	    IntegratedYieldname_eta_syst.ReplaceAll("Eta","Eta_Syst");
	    Integral_eta_syst[g][j] = (TGraphAsymmErrors*)fin[g-1]->Get(IntegratedYieldname_eta_syst)->Clone(IntegratedYieldname_eta_syst);
	 
	    TString IntegratedYieldname_phi = IntegratedYieldname_eta;
	    IntegratedYieldname_phi.ReplaceAll("Eta","Phi");
	    Integral_phi_Pt[g][j] = (TH1D*)fin[g-1]->Get(IntegratedYieldname_phi)->Clone(IntegratedYieldname_phi);
	    Integral_phi_Pt[g][j]->SetAxisRange(1.01,7.8,"x");	  
 


	    TString IntegratedYieldname_phi_syst = IntegratedYieldname_eta_syst;
	    IntegratedYieldname_phi_syst.ReplaceAll("Eta","Phi");
	    Integral_phi_syst[g][j] = (TH1D*)fin[g-1]->Get(IntegratedYieldname_phi_syst)->Clone(IntegratedYieldname_phi_syst);

	
	
	
	    TString ErrorUpEtaPt_name = IntegratedYieldname_eta;
	    ErrorUpEtaPt_name.ReplaceAll("Integrated_Yield_Eta","Integral_Error_Up");
	    Error_up_eta_pT[g][j] = (TGraph*)fin[g-1]->Get(ErrorUpEtaPt_name)->Clone(ErrorUpEtaPt_name);

	    TString ErrorDownEtaPt_name = ErrorUpEtaPt_name;
	    ErrorDownEtaPt_name.ReplaceAll("Up","Down");
	    Error_down_eta_pT[g][j] = (TGraph*)fin[g-1]->Get(ErrorDownEtaPt_name)->Clone(ErrorDownEtaPt_name);
		
	  } //end no pT classes for integrals
	 
	} //end integrals only for  g==1,3,5

	  */
      } //end j
 
    } //end i
  
      /*
    //*****************************
    //   Plot integrals
    //****************************
     
    if(g==1||g==3||g==5){
    
      //-----------------------------------
      //    pT integral plotting loop
      //----------------------------------


     

      for(int j=0;j<4;j++){

	cintegral_eta_pT[g]->cd(j+1);
      

	Integral_eta_Pt[g][j]->SetMarkerColor(1);
	Integral_eta_Pt[g][j]->SetLineColor(1);
	Integral_eta_syst[g][j]->SetMarkerSize(0);
	Integral_eta_Pt[g][j]->SetMarkerSize(2);
	Integral_eta_ref_Pt[g][j]->SetMarkerSize(2);

	switch(g){
	case 1: 
	  Integral_eta_Pt[g][j]->SetMarkerStyle(10);
	  Integral_eta_ref_Pt[g][j]->SetMarkerStyle(10);
	  Integral_eta_syst[g][j]->SetFillColor(90);
	  Error_up_eta_pT[g][j]->SetMarkerStyle(10);
	  Error_up_eta_pT[g][j]->SetFillColor(90);
	  Error_down_eta_pT[g][j]->SetMarkerStyle(10);
	  Error_down_eta_pT[g][j]->SetFillColor(90);
	  break;
	case 3:
	  Integral_eta_Pt[g][j]->SetMarkerStyle(34);
	  Integral_eta_ref_Pt[g][j]->SetMarkerStyle(34);
	  Integral_eta_syst[g][j]->SetFillColor(30);
	  Error_up_eta_pT[g][j]->SetMarkerStyle(34);
	  Error_up_eta_pT[g][j]->SetFillColor(30);
	  Error_down_eta_pT[g][j]->SetMarkerStyle(34);
	  Error_down_eta_pT[g][j]->SetFillColor(30);
	  break;
	case 5: 
	  Integral_eta_Pt[g][j]->SetMarkerStyle(21);
	  Integral_eta_ref_Pt[g][j]->SetMarkerStyle(21);
	  Integral_eta_syst[g][j]->SetFillColor(kOrange-2);
	  Error_up_eta_pT[g][j]->SetMarkerStyle(21);
	  Error_up_eta_pT[g][j]->SetFillColor(kOrange-2);
	  Error_down_eta_pT[g][j]->SetMarkerStyle(21);
	  Error_down_eta_pT[g][j]->SetFillColor(kOrange-2);
	  break;
	default:
	  Integral_eta_Pt[g][j]->SetMarkerStyle(10);
	  break;
	}

	Error_down_eta_pT[g][j]->SetMarkerColor(kRed);


	if(g==1){Integral_eta_Pt[g][j]->SetMaximum(4.2);
	  Integral_eta_Pt[g][j]->SetMinimum(-.45);}
	if(g==3){Integral_eta_Pt[g][j]->SetMaximum(10.2);
	  Integral_eta_Pt[g][j]->SetMinimum(-1.);}
	if(g==5){Integral_eta_Pt[g][j]->SetMaximum(4.2);
	  Integral_eta_Pt[g][j]->SetMinimum(-.45);}


	    
	Integral_eta_Pt[g][j]->Draw("p");
	 
	
	Integral_eta_Pt[g][j]->GetYaxis()->SetLabelSize(tstitle);
	   
	Integral_eta_Pt[g][j]->GetXaxis()->SetLabelSize(tstitle);
	Integral_eta_Pt[g][j]->GetXaxis()->SetTitle("Track p_{T} (GeV)");
	Integral_eta_Pt[g][j]->GetXaxis()->SetTitleSize(tstitle);
	Integral_eta_Pt[g][j]->GetXaxis()->SetTitleOffset(xoffset);
	Integral_eta_Pt[g][j]->GetYaxis()->SetTitle("(dN/dp_{T})_{PbPb}- (dN/dp_{T})_{pp}   (GeV^{-1})");
	Integral_eta_Pt[g][j]->GetYaxis()->SetTitleOffset(yoffset);
	Integral_eta_Pt[g][j]->GetYaxis()->SetTitleSize(tstitle2);

 
	Integral_eta_Pt[g][j]->GetXaxis()->SetRangeUser(1.000001,7.7);
	Integral_eta_Pt[g][j]->GetXaxis()->CenterTitle();
	Integral_eta_Pt[g][j]->GetYaxis()->CenterTitle();
	   
	if(j>0){
	  Integral_eta_Pt[g][j]->GetYaxis()->SetTitleSize(0.0);
	  Integral_eta_Pt[g][j]->GetYaxis()->SetLabelSize(0.0);
	}



	if(g==5){
	   
	  Integral_eta_syst[g][j]->Draw("same 2");
	  Integral_eta_Pt[g][j]->Draw("same p");
	}
	    

	if(g==3){

	  Integral_eta_syst[g][j]->Draw("same 2");
	  Integral_eta_Pt[g][j]->Draw("same p");
	   
	}

	if(g==1){
	  Integral_eta_syst[g][j]->Draw("same 2");
	  Integral_eta_Pt[g][j]->Draw("same p");

	}
	 
	if(g>5){
	  Integral_eta_syst[g][j]->Draw("same2");

	}
	if(g==1||g==3||g==5){
	  Integral_eta_ref_Pt[g][j]->SetLineColor(kRed);
	  Integral_eta_ref_Pt[g][j]->SetMarkerColor(kRed);
	  Integral_eta_ref_Pt[g][j]->Draw("same p");
	}
	if(j==0){ 
	  l40 = new TLegend(textalign2-0.05,texty2-.1,0.6,texty4-0.05);
	  l40->SetName("l40");
	  l40->SetTextFont(63);
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){l40->AddEntry(Integral_eta_Pt[g][j],"Inclusive Jets","lp");
	  }
	  if(g==3){l40->AddEntry(Integral_eta_Pt[g][j],"SubLeading Jets","lp");
	  }
	  if(g==5){l40->AddEntry(Integral_eta_Pt[g][j],"Leading Jets","lp");
	  }
	  //	  l40->AddEntry( Integral_eta_ref_Pt[g][j],"PAS Reference", "lp");
	  l40->Draw("same");
	}
	
	
	    
	linePt = new TLine(1.,0,8.,0);
	linePt->SetLineStyle(2);
	linePt->SetLineWidth(2);
	linePt->Draw("same");

	drawlabels_int_pt(g,j);



	cintegral_phi_pT[g]->cd(j+1);
    
	Integral_phi_Pt[g][j]->SetMarkerColor(1);
	Integral_phi_Pt[g][j]->SetLineColor(1);
	Integral_phi_syst[g][j]->SetMarkerSize(0);
	Integral_phi_Pt[g][j]->SetMarkerSize(2);
  
       
	  
	    
	Integral_phi_Pt[g][j]->SetMarkerSize(2);
	Integral_phi_Pt[g][j]->SetMarkerColor(1);
	Integral_phi_Pt[g][j]->SetLineColor(1);
	   

	switch(g){
	case 1: 
	  Integral_phi_Pt[g][j]->SetMarkerStyle(10);
	  Integral_phi_syst[g][j]->SetFillColor(90);
	  break;
	case 3:
	  Integral_phi_Pt[g][j]->SetMarkerStyle(34);
	  Integral_phi_syst[g][j]->SetFillColor(30);
	  break;
	case 5: 
	  Integral_phi_Pt[g][j]->SetMarkerStyle(21);
	  Integral_phi_syst[g][j]->SetFillColor(kOrange-2);
	  break;
	default:
	  Integral_phi_Pt[g][j]->SetMarkerStyle(10);
	  break;
	}


	Integral_phi_Pt[g][j]->SetMinimum(-1.);
	Integral_phi_Pt[g][j]->SetMaximum(10.2);
	 

	Integral_phi_Pt[g][j]->SetMarkerColor(kRed);
	Integral_phi_Pt[g][j]->SetLineColor(kRed);
	Integral_phi_Pt[g][j]->SetMarkerStyle(29);
	Integral_phi_Pt[g][j]->SetMarkerSize(5);
	Integral_phi_Pt[g][j]->Draw("p");
	  


	Integral_phi_Pt[g][j]->GetYaxis()->SetLabelSize(tstitle2);
	   
	Integral_phi_Pt[g][j]->GetXaxis()->SetLabelSize(tstitle);
	Integral_phi_Pt[g][j]->GetXaxis()->SetTitle("Track p_{T} (GeV)");
	Integral_phi_Pt[g][j]->GetXaxis()->SetTitleSize(tstitle);
	Integral_phi_Pt[g][j]->GetXaxis()->SetTitleOffset(xoffset);
	Integral_phi_Pt[g][j]->GetYaxis()->SetTitle("Y_{PbPb}-Y_{pp}   (GeV^{-1})");
	Integral_phi_Pt[g][j]->GetYaxis()->SetTitleOffset(yoffset);
	Integral_phi_Pt[g][j]->GetYaxis()->SetTitleSize(tstitle2);

 
	Integral_phi_Pt[g][j]->GetXaxis()->SetRangeUser(1.000001,7.7);
	Integral_phi_Pt[g][j]->GetXaxis()->CenterTitle();
	Integral_phi_Pt[g][j]->GetYaxis()->CenterTitle();
	   
	if(j>0){
	  Integral_phi_Pt[g][j]->GetYaxis()->SetTitleSize(0.0);
	  Integral_phi_Pt[g][j]->GetYaxis()->SetLabelSize(0.0);
	}
	  
	Integral_eta_Pt[g][j]->Draw("same");


	if(j==0){
	  l40 = new TLegend(textalign2-0.05,texty2-.1,0.6,texty4-0.05);
	  l40->SetName("l40");
	  l40->SetTextSizePixels(tspixels);
	  l40->SetTextFont(43);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){l40->AddEntry(Integral_eta_Pt[g][j],"Inclusive Jets Integrated #Delta#eta","lp");
	    l40->AddEntry(Integral_phi_Pt[g][j],"Inclusive Jets Integrated #Delta#phi","lp");}
	  if(g==3){l40->AddEntry(Integral_eta_Pt[g][j],"SubLeading Jets Int. #Delta#eta","lp");
	    l40->AddEntry(Integral_phi_Pt[g][j],"SubLeading Jets Int. #Delta#phi","lp");}
	  if(g==5){l40->AddEntry(Integral_eta_Pt[g][j],"Leading Jets Integrated #Delta#eta","lp");
	    l40->AddEntry(Integral_phi_Pt[g][j],"Leading Jets Integrated #Delta#phi","lp");}
	  l40->Draw("same");
	}
	



	//	if(j==0){l40->Draw("same");}
	    
	linePt->Draw("same");

	  
	drawlabels_int_pt(g,j);



	// Draw error up evolution plot

	cError_up_eta_pT[g]->cd(j+1);
	
	gStyle->SetOptTitle(0);

	Error_up_eta_pT[g][j]->SetMinimum(-0.1);
	Error_up_eta_pT[g][j]->SetMaximum(4.);
	
	
	    
	Error_up_eta_pT[g][j]->Draw("p");
	


	Error_up_eta_pT[g][j]->GetYaxis()->SetLabelSize(tstitle);
	   
	Error_up_eta_pT[g][j]->GetXaxis()->SetLabelSize(tstitle);
	Error_up_eta_pT[g][j]->GetXaxis()->SetTitle("Track p_{T} (GeV)");
	Error_up_eta_pT[g][j]->GetXaxis()->SetTitleSize(tstitle);
	Error_up_eta_pT[g][j]->GetXaxis()->SetTitleOffset(xoffset);
	Error_up_eta_pT[g][j]->GetYaxis()->SetTitle("Y_{PbPb}-Y_{pp}   (GeV^{-1})");
	Error_up_eta_pT[g][j]->GetYaxis()->SetTitleOffset(yoffset);
	Error_up_eta_pT[g][j]->GetYaxis()->SetTitleSize(tstitle2);

 
	Error_up_eta_pT[g][j]->GetXaxis()->SetRangeUser(1.000001,7.7);
	Error_up_eta_pT[g][j]->GetXaxis()->CenterTitle();
	Error_up_eta_pT[g][j]->GetYaxis()->CenterTitle();
	   
	if(j>0){
	  Error_up_eta_pT[g][j]->GetYaxis()->SetTitleSize(0.0);
	  Error_up_eta_pT[g][j]->GetYaxis()->SetLabelSize(0.0);
	}
	
	Error_up_eta_pT[g][j]->SetMarkerSize(2);
	Error_up_eta_pT[g][j]->Draw();
		
	drawlabels_int_pt(g,j);

	if(j==0){ 
	  l40 = new TLegend(textalign2-0.05,texty2-.1,0.6,texty4-0.05);
	  l40->SetName("l40");
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){l40->AddEntry(Error_up_eta_pT[g][j],"Inclusive ErrorUp","lp");
	    l40->AddEntry(Error_down_eta_pT[g][j],"Inclusive ErrorDown","lp");}
	  if(g==3){l40->AddEntry(Error_up_eta_pT[g][j],"SubLeading ErrorDown","lp");
	    l40->AddEntry(Error_down_eta_pT[g][j],"SubLeading ErrorDown","lp");}
	  if(g==5){l40->AddEntry(Error_up_eta_pT[g][j],"Leading ErrorUp","lp");
	    l40->AddEntry(Error_down_eta_pT[g][j],"Leading ErrorDown","lp");}
	  l40->Draw("same");
	}
	
	linePt->Draw("same");
	
	Error_down_eta_pT[g][j]->SetLineColor(kRed);
	Error_down_eta_pT[g][j]->SetMarkerSize(2);
	Error_down_eta_pT[g][j]->Draw("same pl");
	


      } //close j


      TString IntegralSaveName_eta = in_name;
      IntegralSaveName_eta.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      IntegralSaveName_eta.ReplaceAll("Result","Integral_Yield");
      IntegralSaveName_eta.ReplaceAll("pp","PbPb_minus_pp_eta");
      IntegralSaveName_eta +=".pdf";
      cintegral_eta_pT[g]->SaveAs(IntegralSaveName_eta);


      TString IntegralSaveName_phi = in_name;
      IntegralSaveName_phi.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      IntegralSaveName_phi.ReplaceAll("Result","Integral_Yield");
      IntegralSaveName_phi.ReplaceAll("pp","PbPb_minus_pp_phi");
      IntegralSaveName_phi +=".pdf";
      cintegral_phi_pT[g]->SaveAs(IntegralSaveName_phi);
    

      TString ErrorUpSaveNamePt = in_name;
      ErrorUpSaveNamePt.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      ErrorUpSaveNamePt.ReplaceAll("Result","Error_pT");
      ErrorUpSaveNamePt.ReplaceAll("pp","PbPb_minus_pp_phi");
      ErrorUpSaveNamePt +=".pdf";
      cError_up_eta_pT[g]->SaveAs(ErrorUpSaveNamePt);

      
    

      //---------------------------------------
      //Now for integral by centrality plots
      //--------------------------------------

      for(int i = 0; i<4; i++){

	cintegral_eta_cent[g]->cd(i+1);
	

	TString histnamecent = "integral_cent";
	histnamecent+=g;
	histnamecent+=i;

	int_cent[g][i] = new TH1D(histnamecent,"",4,xAxis);

	TString histnameblank = "blank_hist";
	histnameblank+=g;
	histnameblank+=i;

	blank[i] = new TH1D(histnameblank,"",4,xAxis);

	
	TString histnameblank2 = "blank_hist2";
	histnameblank2+=g;
	histnameblank2+=i;

	blank2[i] = new TH1D(histnameblank2,"",4,xAxis);

	
	for(int k=0; k<4; k++){

	  value = Integral_eta_Pt[g][k]->GetBinContent(i+1);
	  error = Integral_eta_syst[g][k]->GetErrorYhigh(i);

	  int_cent[g][i]->SetBinContent(k+1,value);
	  int_cent[g][i]->SetBinError(k+1,error);

	}

	

	TString ErrorUpEtaCent_name = "ErrorUpEtaCent";
	ErrorUpEtaCent_name+=g;
	ErrorUpEtaCent_name+=i;
	Error_up_eta_cent[g][i] = new TH1D(ErrorUpEtaCent_name,"",4,xAxis);

	TString ErrorDownEtaCent_name = "ErrorDownEtaCent";
	ErrorDownEtaCent_name+=g;
	ErrorDownEtaCent_name+=i;
	Error_down_eta_cent[g][i] = new TH1D(ErrorDownEtaCent_name,"",4,xAxis);

		
	for(int k=0; k<4; k++){
	  
	  evalpt = pTbin_centers.at(i);

	  value = Error_up_eta_pT[g][k]->Eval(evalpt);
	 
	  Error_up_eta_cent[g][i]->SetBinContent(k+1,value);
	  
	  value = Error_down_eta_pT[g][k]->Eval(evalpt);
	  
	  Error_down_eta_cent[g][i]->SetBinContent(k+1,value);
	  
	}


	//Markers etc. for all histograms at once
	//----------------------------------------


	int_cent[g][i]->SetMarkerSize(2);
	int_cent[g][i]->SetMarkerColor(1);
	int_cent[g][i]->SetLineColor(1);

	switch(g){
	case 1:
	  int_cent[g][i]->SetFillColor(90);
	  int_cent[g][i]->SetMarkerStyle(10);
	  Error_up_eta_cent[g][i]->SetMarkerStyle(10);
	  Error_down_eta_cent[g][i]->SetMarkerStyle(10);
	  break;     
	case 3:
	  int_cent[g][i]->SetFillColor(30);
	  int_cent[g][i]->SetMarkerStyle(34);
	  Error_up_eta_cent[g][i]->SetMarkerStyle(34);
	  Error_down_eta_cent[g][i]->SetMarkerStyle(34);
	  break;  
   	case 5:
	  int_cent[g][i]->SetFillColor(kOrange-2);
	  int_cent[g][i]->SetMarkerStyle(21);
	  Error_up_eta_cent[g][i]->SetMarkerStyle(21);
	  Error_down_eta_cent[g][i]->SetMarkerStyle(21);
	  break;     
	}

	Error_down_eta_cent[g][i]->SetMarkerColor(kRed);
	
	//Plot aesthetics for every canvas.
	//----------------------------------
	
	blank[i]->SetMinimum(-1.);
	blank[i]->SetMaximum(10.2);
	blank[i]->GetXaxis()->SetTitle("Centrality (%)");
	blank[i]->GetXaxis()->SetTitleOffset(1.1);
	blank[i]->GetXaxis()->CenterTitle(true);
	blank[i]->GetXaxis()->SetTitleSize(ts);

	blank[i]->GetYaxis()->SetTitle("(dN/dp_{T})_{PbPb}- (dN/dp_{T})_{pp} (GeV^{-1})");
	blank[i]->GetYaxis()->SetTitleSize(0.);
	blank[i]->GetYaxis()->CenterTitle(true);
	//	blank[i]->GetYaxis()->SetLabelOffset(yoffset);
	//	blank[i]->GetYaxis()->SetLabelSize(0.);
   
	blank[i]->GetYaxis()->SetTickLength(0.025);

	blank[i]->GetXaxis()->SetBinLabel(1,"50-100");
	blank[i]->GetXaxis()->SetBinLabel(2,"30-50");
	blank[i]->GetXaxis()->SetBinLabel(3,"10-30");
	blank[i]->GetXaxis()->SetBinLabel(4," 0-10");
    
	blank[i]->GetXaxis()->SetLabelSize(0.08);
	blank[i]->GetXaxis()->SetLabelOffset(0.015);
	
	blank[i]->GetXaxis()->LabelsOption("h");
	blank[i]->GetXaxis()->SetTickLength(0.0);

	switch(i){
	case 0: 
	  //  gPad->SetLeftMargin(0.2);
	  blank[i]->GetYaxis()->SetTitleSize(ts);
	  blank[i]->GetXaxis()->SetTitleOffset(1.1);
	  blank[i]->GetXaxis()->SetTitleSize(0.07);
	  blank[i]->SetLabelSize(0.95*blank[i]->GetXaxis()->GetLabelSize());
	  blank[i]->GetYaxis()->SetLabelSize(ts2);
	  break;
	case 3:
	  // gPad->SetRightMargin(0.02);
	  break;
	default:
	  break;
	}


	TLine *line1, *line2, *line3;


	switch(i){
	case 0:
	  blank[i]->SetMaximum(6.1);
	  line1 = new TLine(-50,-0.5,-50.,-0.1);
	  line2 = new TLine(-30,-0.5,-30.,-0.1);
	  line3 = new TLine(-10,-0.5,-10.,-0.1);
	  break;
	case 1:	 
	  blank[i]->SetMaximum(1.8);
	  line1 = new TLine(-50,-0.38,-50.,-0.22);
	  line2 = new TLine(-30,-0.38,-30.,-0.22);
	  line3 = new TLine(-10,-0.38,-10.,-0.22);	 
	  break;
	case 2:
	  blank[i]->SetMaximum(0.7);
	  line1 = new TLine(-50,-0.34,-50.,-0.26);
	  line2 = new TLine(-30,-0.34,-30.,-0.26);
	  line3 = new TLine(-10,-0.34,-10.,-0.26);
	  break;
	case 3:	
	  blank[i]->SetMaximum(0.4);
	  line1 = new TLine(-50,-0.33,-50.,-0.27);
	  line2 = new TLine(-30,-0.33,-30.,-0.27);
	  line3 = new TLine(-10,-0.33,-10.,-0.27);
	  break;
	}

	blank[i]->GetYaxis()->SetLabelSize(ts);

	//	blank[i]->SetMaximum(1.1);
	blank[i]->SetMinimum(-.3);

	blank[i]->GetYaxis()->SetNdivisions(408);
	blank[i]->GetYaxis()->SetTitleOffset(0.8);
	blank[i]->Draw();


	
	line1->Draw();
	line2->Draw();
	line3->Draw();
	


	
	if(i==0){ 
	  l40 = new TLegend(.17,texty2-.05,.8,texty3-0.05);
	  l40->SetName("l40");
	  l40->SetTextFont(tspixels);
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){l40->AddEntry(int_cent[g][i],"Inclusive Jets","lpfe");}
	  if(g==3){l40->AddEntry(int_cent[g][i],"Subleading Jets","lpfe");}
	  if(g==5){l40->AddEntry(int_cent[g][i],"Leading Jets","lpfe");}
	  l40->Draw("same");
	}
	drawlabels_int_cent2(g,i);
	 
	lineCent = new TLine(-100.,0.,0.,0.);
	lineCent->SetLineStyle(2);
	lineCent->SetLineWidth(2);
	lineCent->Draw("same");
	
	int_cent[g][i]->Draw("same bp0 e2");

	lineCent->Draw("same");

	gPad->RedrawAxis();

	//----------------------------
	// Draw error up evolution plot
	//----------------------------
	cError_up_eta_cent[g]->cd(i+1);
	

	gStyle->SetOptTitle(0);


	blank2[i]->GetXaxis()->SetTitle("Centrality (%)");
	blank2[i]->GetXaxis()->SetTitleOffset(1.1);
	blank2[i]->GetXaxis()->CenterTitle(true);
	blank2[i]->GetXaxis()->SetTitleSize(ts);

	blank2[i]->GetYaxis()->SetTitle("Y_{PbPb}-Y_{pp}   (GeV^{-1})");
	blank2[i]->GetYaxis()->SetTitleSize(0.);
	blank2[i]->GetYaxis()->CenterTitle(true);
	//	blank2[i]->GetYaxis()->SetLabelOffset(yoffset);
	blank2[i]->GetYaxis()->SetLabelSize(0.);
   
	blank2[i]->GetYaxis()->SetTickLength(0.025);

	blank2[i]->GetXaxis()->SetBinLabel(1,"50-100");
	blank2[i]->GetXaxis()->SetBinLabel(2,"30-50");
	blank2[i]->GetXaxis()->SetBinLabel(3,"10-30");
	blank2[i]->GetXaxis()->SetBinLabel(4," 0-10");
    
	blank2[i]->GetXaxis()->SetLabelSize(0.09);
	blank2[i]->GetXaxis()->SetLabelOffset(0.015);
	blank2[i]->GetXaxis()->SetTicks("+-");
	blank2[i]->GetXaxis()->LabelsOption("h");
	blank2[i]->GetXaxis()->SetTickLength(0.025);

	switch(i){
	case 0: 
	  gPad->SetLeftMargin(0.2);
	  blank2[i]->GetYaxis()->SetTitleSize(0.09);
	  blank2[i]->GetXaxis()->SetTitleOffset(1.2);
	  blank2[i]->GetXaxis()->SetTitleSize(0.06);
	  blank2[i]->SetLabelSize(0.9*blank2[i]->GetXaxis()->GetLabelSize());
	  blank2[i]->GetYaxis()->SetLabelSize(ts);
	  break;
	case 3:
	  gPad->SetRightMargin(0.02);
	  break;
	default:
	  break;
	}
	
	
	blank2[i]->SetMinimum(-0.1);
	blank2[i]->SetMaximum(4.);


	blank2[i]->Draw();





	Error_up_eta_cent[g][i]->SetMarkerSize(2);
	Error_up_eta_cent[g][i]->Draw("same p");
	Error_up_eta_cent[g][i]->Draw("same pl");
		
	
	Error_down_eta_cent[g][i]->SetLineColor(kRed);
	Error_down_eta_cent[g][i]->SetMarkerSize(2);
	Error_down_eta_cent[g][i]->Draw("same pl");
	


	drawlabels_int_cent2(g,i);

	if(i==0){ 
	  l40 = new TLegend(textalign2,texty2-.1,0.6,texty4-0.05);
	  l40->SetName("l40");
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){l40->AddEntry(Error_up_eta_pT[g][i],"Inclusive ErrorUp","lp");
	    l40->AddEntry(Error_down_eta_pT[g][i],"Inclusive ErrorDown","lp");}
	  if(g==3){l40->AddEntry(Error_up_eta_pT[g][i],"SubLeading ErrorDown","lp");
	    l40->AddEntry(Error_down_eta_pT[g][i],"SubLeading ErrorDown","lp");}
	  if(g==5){l40->AddEntry(Error_up_eta_pT[g][i],"Leading ErrorUp","lp");
	    l40->AddEntry(Error_down_eta_pT[g][i],"Leading ErrorDown","lp");}
	  l40->Draw("same");
	}	
	
	lineCent->Draw("same");
	
	
	
      } //close i



      cintegral_eta_cent[g]->cd(0);

      TLatex *canvas_title = new TLatex(0.05,0.9,"CMS Preliminary");
      canvas_title->SetTextSizePixels(tspixels);
      canvas_title->SetTextFont(63);
      canvas_title->Draw();

      TLatex *canvas_title2 = new TLatex(0.295,0.9,"PbPb 166 #mub^{-1} (2.76 TeV)                     pp 5.3 pb^{-1} (2.76 TeV)");
      canvas_title2->SetTextSizePixels(tspixels);
      canvas_title2->Draw();


      TString IntegralSaveName_cent = in_name;
      IntegralSaveName_cent.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      IntegralSaveName_cent.ReplaceAll("Result","Integral_Yield");
      IntegralSaveName_cent.ReplaceAll("pp","PbPb_minus_pp_cent");
      IntegralSaveName_cent +=".pdf";
      cintegral_eta_cent[g]->SaveAs(IntegralSaveName_cent);
      IntegralSaveName_cent.ReplaceAll(".pdf",".png");
      cintegral_eta_cent[g]->SaveAs(IntegralSaveName_cent);



      TString ErrorUpSaveName_cent = in_name;
      ErrorUpSaveName_cent.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      ErrorUpSaveName_cent.ReplaceAll("Result","Error_Cent");
      ErrorUpSaveName_cent.ReplaceAll("pp","PbPb_minus_pp_cent");
      ErrorUpSaveName_cent +=".pdf";
      cError_up_eta_cent[g]->SaveAs(ErrorUpSaveName_cent);


    }// close "only integrate for pp"

  } // close g

  cout<<"Now starting on PAS Plots..."<<endl;

      */
    }
  }



  //*************************************
  //   Now draw PAS plots
  //***********************************

  for(int i = 0;i <4; i++){
   
    
    TString pTrange;

    switch(i){
    case 0: 
      pTrange = "TrkPt1_TrkPt2"; 
      raw_min = 10.;
      raw_max = 50.;
      mixed_min = 10./34.;
      mixed_max = 50/34.;
      yield_min = 35.;
      yield_max = 43.;
      result_min = -.5;
      result_max = 8.5;
      break;
    case 1: 
      pTrange = "TrkPt2_TrkPt3"; 
      raw_min = 1.;
      raw_max = 10.;
      mixed_min = 1./4.5;
      mixed_max = 10./4.5;
      yield_min = 4.;
      yield_max = 10.;
      result_min = -.45;
      result_max = 8.5;
      break;
    case 2: 
      pTrange = "TrkPt3_TrkPt4"; 
      raw_min = 0.;
      raw_max = 5.;
      mixed_min = 0.;
      mixed_max = 6.;
      yield_min = 0.;
      yield_max = 5.;
      result_min = -.45;
      result_max = 8.5;
      break;
    case 3: 
      pTrange = "TrkPt4_TrkPt8"; 
      raw_min = 0.;
      raw_max = 2.;
      mixed_min = 0.;
      mixed_max = 10.;
      yield_min = -.1;
      yield_max = 1.9;
      result_min = -.1;
      result_max = 1.9;
      break;
    }
       

    TF1 *gen_gaus = new TF1("gen_gaus",
			  " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))  \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);
   
    TF1 *gen_gaus_up = new TF1("gen_gaus_up",
			       " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6])) \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6]))) ",-TMath::Pi()/2,3*TMath::Pi()/2);
   
  

    TF1 *gen_gaus_down = new TF1("gen_gaus_down",
				 " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6])) \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);
   

    TF1 *gen_gaus_level = new TF1("gen_gaus_level",
			  " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))  \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);



    TF1 *gen_gaus_bg = new TF1("gen_gaus_bg",
			  " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))  \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);

     TF1 *sub_gen_gaus = new TF1("sub_gen_gaus",
			    " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6])) \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);
   
    TF1 *sub_gen_gaus_up = new TF1("sub_gen_gaus_up",
			       " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6])) \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);
   
  

    TF1 *sub_gen_gaus_down = new TF1("sub_gen_gaus_down",
				 " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6])) \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);
   
    

    TF1 *sub_gen_gaus_bg = new TF1("sub_gen_gaus_bg",
			  " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))  \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);

     
    TString figure3_name = "PAS_Figure_3_";
    figure3_name+=pTrange;
    cFigure3[i] = new TCanvas(figure3_name," ",10,10,1550,850);
    cFigure3[i]->Divide(4,2,0.,0.);

    TString figure4_name = "PAS_Figure_4_";
    figure4_name+=pTrange;
    cFigure4[i] = new TCanvas(figure4_name," ",10,10,1550,850);
    cFigure4[i]->Divide(4,2,0.,0.);

    int ndivisions;
    float diff_min, diff_max;

    switch(i){
    case 0: 
      result_min = -.8;
      result_max = 10.7;
      diff_min = -.8;
      diff_max = 10.7;
      ndivisions = 408;
      break;
    case 1: 
      result_min = -.8;
      result_max = 7.9;
      diff_min = -.8;
      diff_max = 7.9;
      ndivisions = 612;
      break;
    case 2: 
      result_min = -.45;
      result_max = 6.6;
      diff_min = -.45-.4;
      diff_max = 6.6-.4;
      ndivisions = 408;
      break;
    case 3: 
      result_min = -.25;
      result_max = 3.9;
      diff_min = -.25-.6;
      diff_max = 3.9-.6;
      ndivisions = 306;
      break;
    }


    for (int j=0; j<4; j++){
		   
	
      //---------------------
      //  PAS FIGURE 3
      //----------------------
	
      cFigure3[i]->cd(j+1);
      check_new_eta_rebin[0][i][j]->SetMaximum(result_max);
      check_new_eta_rebin[0][i][j]->SetMinimum(result_min);
      check_new_eta_rebin[0][i][j]->GetXaxis()->SetLabelSize(0.);
      check_new_eta_rebin[0][i][j]->GetYaxis()->SetNdivisions(ndivisions);
          
      check_new_eta_rebin[0][i][j]->GetYaxis()->CenterTitle(1);
      check_new_eta_rebin[0][i][j]->Draw();
   
      check_new_eta_syst[0][i][j]->SetMarkerSize(1.2);
     
      check_new_eta_syst[0][i][j]->SetMarkerStyle(20);
      check_new_eta_syst[0][i][j]->Draw("same e2");
     
 
      check_new_eta_rebin[0][i][j]->Draw("same");
      check_new_eta_rebin[1][i][j]->SetMarkerStyle(24);
      check_new_eta_rebin[1][i][j]->SetMarkerSize(1.2);
      check_new_eta_rebin[1][i][j]->Draw("same");
     
      if(draw_ref) check_new_eta_ref[0][i][j]->Draw("same");
      if(draw_ref) check_new_eta_ref[1][i][j]->Draw("same");
      drawlabels(1,i,j);

       

      if(j==3){
	tex24eta = new TLatex(textalign,texty2,phirangelabel);
	tex24eta->SetName("tex24eta");
	tex24eta->SetNDC();
	tex24eta->SetTextSizePixels(tspixels);
	tex24eta->Draw();
      }

      lineEta = new TLine(-1.5,0,1.5,0); 
      lineEta->Draw("same");
	
      if(j==0){ 
	l40 = new TLegend(textalign2-0.03,texty3-0.1,.95,texty2-0.05);
	l40->SetName("l40");
	l40->SetTextSizePixels(tspixels);
	l40->SetFillColor(kWhite);
	l40->SetLineColor(kWhite);
	l40->AddEntry(check_new_eta_syst[0][i][j],"PbPb Inclusive Jets","lpfe");
	l40->AddEntry(check_new_eta_rebin[1][i][j],"pp Inclusive Jets","lpfe");
	l40->Draw("same");
      }

   if(j==1){ 
	l40 = new TLegend(textalign-0.03,texty3-0.1,.95,texty2-0.05);
	l40->SetName("l40");
	l40->SetTextSizePixels(tspixels);
	l40->SetFillColor(kWhite);
	l40->SetLineColor(kWhite);
	l40->AddEntry(check_new_eta_ref[0][i][j],"PbPb Inclusive Jets","lpfe");
	l40->AddEntry(check_new_eta_ref[1][i][j],"pp Inclusive Jets","lpfe");
	l40->Draw("same");
      }


      cFigure3[i]->cd(j+5);
      PbPb_pp_eta_syst[1][i][j]->GetYaxis()->SetTitle("Y_{PbPb}-Y_{pp}   (GeV^{-1})");


      PbPb_pp_eta_syst[1][i][j]->SetMinimum(diff_min);
      PbPb_pp_eta_syst[1][i][j]->SetMaximum(diff_max);
      PbPb_pp_eta_syst[1][i][j]->SetMarkerSize(1.2);
      PbPb_pp_eta_syst[1][i][j]->GetYaxis()->SetNdivisions(ndivisions);

      PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetLabelSize(ts2);
      PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetTitleSize(tstitle);
      if(j==0){ PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetLabelSize(ts2-0.01);
	PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetTitleSize(ts);
	PbPb_pp_eta_syst[1][i][j]->GetYaxis()->SetTitleSize(ts);
	PbPb_pp_eta_syst[1][i][j]->GetYaxis()->SetLabelSize(ts2);
	PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetLabelOffset(0.013);
	PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetTitleOffset(xoffset2);
      }
      //PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetLabelOffset(xoffset);
      PbPb_pp_eta_syst[1][i][j]->SetFillColor(90);
      //  PbPb_pp_eta_syst[1][i][j]->SetMarkerSize(0);
      PbPb_pp_eta_syst[1][i][j]->Draw("e2");
      PbPb_pp_eta[1][i][j]->Draw("same");

      if(draw_ref) PbPb_pp_eta_ref[1][i][j]->Draw("same");
	 
      lineEta->Draw("same");

      TPave *cover_x = new TPave(-1.9,-1.8,-1.2,-.88);
      
      cover_x->SetFillColor(kWhite);
      cover_x->SetShadowColor(kWhite);
      cover_x->SetLineColor(kWhite);
      cover_x->Draw();
  
      TPave *cover_x2 = new TPave(1.1,-1.8,1.8,-.88);
      
      cover_x2->SetFillColor(kWhite);
      cover_x2->SetShadowColor(kWhite);
      cover_x2->SetLineColor(kWhite);
      cover_x2->Draw();

    	
 


      TPave *labelcover = new TPave(1.3,diff_min-0.8,1.6,diff_min-0.05);
      labelcover->SetLineColor(kWhite);
      labelcover->SetOption("nb");
      labelcover->SetFillColor(kWhite);
      if(j<3){ labelcover->Draw(); }

      if(j==0){ 
	l40 = new TLegend(textalign2-0.03,texty2-0.15,.95,texty1-0.05);
	l40->SetName("l40");
	l40->SetTextSizePixels(tspixels);
	l40->SetFillColor(kWhite);
	l40->SetLineColor(kWhite);
	l40->AddEntry(PbPb_pp_eta_syst[1][i][j],"PbPb - pp Inclusive Jets","lpfe");
	l40->AddEntry(PbPb_pp_eta_ref[1][i][j],"PbPb - pp 5.02 TeV","lpfe");
	l40->Draw("same");
      }

	    
      
      //-----------
      // Figure 4
      //-----------


      cFigure4[i]->cd(j+1);
      check_new_phi_rebin[0][i][j]->SetMaximum(result_max);
      check_new_phi_rebin[0][i][j]->SetMinimum(result_min);
      check_new_phi_rebin[0][i][j]->GetYaxis()->SetNdivisions(ndivisions);
      check_new_phi_rebin[0][i][j]->GetXaxis()->SetLabelSize(0.);
      check_new_phi_rebin[0][i][j]->GetYaxis()->CenterTitle(1);
      check_new_phi_rebin[0][i][j]->Draw();
      check_new_phi_syst[0][i][j]->SetMarkerSize(1.2);
      check_new_phi_rebin[1][i][j]->SetMarkerSize(1.2);
      check_new_phi_syst[0][i][j]->SetMarkerStyle(20);
      check_new_phi_syst[0][i][j]->Draw("same e2");
      check_new_phi_rebin[0][i][j]->Draw("same");

      if(draw_ref) check_new_phi_ref[0][i][j]->Draw("same");
      if(draw_ref) check_new_phi_ref[1][i][j]->Draw("same");

      check_new_phi_rebin[1][i][j]->Draw("same");
	 
      drawlabels(1,i,j);
	  
	
      if(j==0){ 
	l40 = new TLegend(textalign2-0.03,texty3-0.1,.95,texty2-0.05);
	l40->SetName("l40");
	l40->SetTextSizePixels(tspixels);
	l40->SetFillColor(kWhite);
	l40->SetLineColor(kWhite);
	l40->AddEntry(check_new_phi_syst[0][i][j],"PbPb Inclusive Jets","lpfe");
	l40->AddEntry(check_new_phi_rebin[1][i][j],"pp Inclusive Jets","lpfe");
	l40->Draw("same");
      }

	
      if(j==1){ 
	l40 = new TLegend(textalign-0.03,texty3-0.1,.95,texty2-0.05);
	l40->SetName("l40");
	l40->SetTextSizePixels(tspixels);
	l40->SetFillColor(kWhite);
	l40->SetLineColor(kWhite);
	l40->AddEntry(check_new_phi_ref[0][i][j],"PbPb 5.02 TeV","lpfe");
	l40->AddEntry(check_new_phi_ref[1][i][j],"pp 5.02 TeV","lpfe");
	l40->Draw("same");
      }

	  
      if(j==3){
	tex24phi = new TLatex(textalign,texty2,etarangelabel);
	tex24phi->SetName("tex24phi");
	tex24phi->SetNDC();
	tex24phi->SetTextSizePixels(tspixels);
	tex24phi->Draw();
      }


      linePhi = new TLine(-TMath::Pi()/2,0.,TMath::Pi()/2,0.);
      linePhi->Draw("same");

  
      cFigure4[i]->cd(j+5);
      PbPb_pp_phi_syst[1][i][j]->GetYaxis()->SetTitle("Y_{PbPb}-Y_{pp}   (GeV^{-1})");
      PbPb_pp_phi_syst[1][i][j]->SetMinimum(diff_min);
      PbPb_pp_phi_syst[1][i][j]->SetMaximum(diff_max);
      PbPb_pp_phi_syst[1][i][j]->SetMarkerSize(1.2);
      PbPb_pp_phi_syst[1][i][j]->GetYaxis()->SetNdivisions(ndivisions);
      PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetRangeUser(-1.4999,1.5);
      PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetLabelSize(ts2);
      PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetTitleSize(tstitle);
     
      if(j==0){ PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetLabelSize(ts2-0.01);
	PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetTitleSize(ts);
	PbPb_pp_phi_syst[1][i][j]->GetYaxis()->SetTitleSize(ts);
	PbPb_pp_phi_syst[1][i][j]->GetYaxis()->SetLabelSize(ts2);
	PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetTitleOffset(xoffset2);
	PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetLabelOffset(0.013);
      }
      PbPb_pp_phi_syst[1][i][j]->SetFillColor(90);
      //   PbPb_pp_phi_syst[1][i][j]->SetMarkerSize(0);
      PbPb_pp_phi_syst[1][i][j]->Draw("e2");
      PbPb_pp_phi[1][i][j]->Draw("same p");

      if(draw_ref) PbPb_pp_phi_ref[1][i][j]->Draw("same p");
      linePhi->Draw("same");


      labelcover = new TPave(1.3,diff_min-0.8,1.6,diff_min-0.05);
      labelcover->SetLineColor(kWhite);
      labelcover->SetOption("nb");
      labelcover->SetFillColor(kWhite);
      if(j<3){ labelcover->Draw(); }

      if(j==0){ 
	l40 = new TLegend(textalign2-0.03,texty2-0.15,.95,texty1-0.05);
	l40->SetName("l40");
	l40->SetTextSizePixels(tspixels);
	l40->SetFillColor(kWhite);
	l40->SetLineColor(kWhite);
	l40->AddEntry(PbPb_pp_phi_syst[1][i][j],"PbPb - pp Inclusive Jets","lpfe");
	l40->AddEntry(PbPb_pp_phi_ref[1][i][j],"PbPb - pp 5.02 TeV","lpfe");
	l40->Draw("same");
      }
  
        cover_x->Draw();
      cover_x2->Draw();

    } //close j
  
  
      //--------------------------------------------
      //   Save all PAS plots except for Figure 7
      //-----------------------------------------
  

    cFigure3[i]->cd(0);
    TLatex *canvas_title = new TLatex(0.06,0.95,"CMS ");
    canvas_title->SetTextSizePixels(tspixels);
    canvas_title->SetTextFont(63);
    canvas_title->Draw();
  
    TLatex *canvas_title2 = new TLatex(0.292,0.95,"PbPb 166 #mub^{-1} (2.76 TeV)                pp 5.3 pb^{-1} (2.76 TeV)");
    canvas_title2->SetTextSizePixels(tspixels);
    canvas_title2->Draw();


    cFigure4[i]->cd(0);
    canvas_title->Draw();
    canvas_title2->Draw();
 
    if(draw_ref){
      figure3_name+="_WithRef";
      figure4_name+="_WithRef";
    }

    figure3_name+=".pdf";
    cFigure3[i]->SaveAs(figure3_name);
    figure3_name.ReplaceAll("pdf","png");
    cFigure3[i]->SaveAs(figure3_name);   
    figure3_name.ReplaceAll("png","C");
    cFigure3[i]->SaveAs(figure3_name);   


    figure4_name+=".pdf";
    cFigure4[i]->SaveAs(figure4_name);
    figure4_name.ReplaceAll("pdf","png");
    cFigure4[i]->SaveAs(figure4_name);   
   figure4_name.ReplaceAll("png","C");
    cFigure4[i]->SaveAs(figure4_name);   
  } //closes i loop for PAS plots

  return 0;
  
} //Close main loop

