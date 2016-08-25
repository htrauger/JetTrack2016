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
#include "TStyle.h"
#include "TLatex.h"
#include "TGraphErrors.h"


#include <iostream>
#include <vector>
#include <fstream>

#include "../../HIN-14-016/HIN_14_016_functions.h"


Int_t jff_residual(bool is_number = kTRUE){


  gROOT->ForceStyle();
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(1);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.15);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);



  const int nCBins= 4;
  const int nPtBins=1;
  const int nTrkPtBins=9;


  enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RightGen, SpilledUnderGen, UnmatchedGen, RightReco, SpilledReco, UnmatchedReco, RecoGenSube0,RecoGenNoSube0,GenGenSube0,GenGenNoSube0,MatchedRecoGenSube0,MatchedRecoGenNoSube0,SwappedRecoGenSube0,SwappedRecoGenNoSube0, UnMatchedRecoGenSube0,UnMatchedRecoGenNoSube0,n_data_mc_types};


  TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RightGenJet_GenTrack","SpilledUnderJet_GenTrack","UnmatchedGenJet_GenTrack","RightRecoJet_GenTrack","SpilledReco_GenTrack","UnmatchedReco_GenTrack","RecoJet_GenTrack_Sube0","RecoJet_GenTrack_NoSube0","GenJet_GenTrack_Sube0","GenJet_GenTrack_NoSube0","MatchedRecoJet_GenTrack_Sube0","MatchedRecoJet_GenTrack_NoSube0","SwappedRecoJet_GenTrack_Sube0","SwappedRecoJet_GenTrack_NoSube0","UnmatchedRecoJet_GenTrack_Sube0","UnmatchedRecoJet_GenTrack_NoSube0",};

  int data_mc_type_code = -999;

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
 TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins+1] = {0.7, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8","TrkPt12","TrkPt16","TrkPt20","TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12","12<pT<16","16<pT<20","pT>20"};

  float x, offset, value;
  TF1 *do_offset = new TF1("do_offset","-1.*[0]+x-x",-3.,3.);

  TPaveText *labels;

  Double_t xAxis[nTrkPtBins] = {-100,-50,-30,-10,0}; 
  TH1D* int_cent[12][nTrkPtBins];
  TH1D* blank[nTrkPtBins];
  TH1D* blank2[nTrkPtBins];

   
  TFile *fin[12];
  TFile *fin_ref[12];
  TFile *fout[12];
  TFile *fclosures[12];

  TH2D *result[12][nTrkPtBins][4];
  TH2D *result2[12][nTrkPtBins][4];


  TH2D* background[12][nTrkPtBins][4];
  TH1D* background_left[12][nTrkPtBins][4];
  TH1D* background_right[12][nTrkPtBins][4];
  TH1D* background_proj[12][nTrkPtBins][4];


  TH1D *phi_proj[12][nTrkPtBins][4];
  TH1D *phi_proj_rebin[12][nTrkPtBins][4];
  TH1D *phi_proj_rebin2[12][nTrkPtBins][4];
  TH1D *eta_proj[12][nTrkPtBins][4];
  TH1D *eta_proj_rebin[12][nTrkPtBins][4];
  TH1D *eta_proj_rebin2[12][nTrkPtBins][4];

  TH1D *eta_proj_ref[12][nTrkPtBins][4];
  TH1D *phi_proj_ref[12][nTrkPtBins][4];


  TH1D *jff_residual_eta[12][nTrkPtBins][4];
  TH1D *jff_residual_phi[12][nTrkPtBins][4];

  TCanvas *corr_canvas_eta[12][nTrkPtBins];
  TCanvas *corr_canvas_phi[12][nTrkPtBins];



  vector<float> pTbin_centers;
  pTbin_centers.push_back(0.85);
  
  pTbin_centers.push_back(1.5);
  pTbin_centers.push_back(2.5);
  pTbin_centers.push_back(3.5);
 
  pTbin_centers.push_back(6.0);
  pTbin_centers.push_back(10.0);
  pTbin_centers.push_back(14.0);
  pTbin_centers.push_back(18.0);
  pTbin_centers.push_back(22.0);
  
  vector<float> pTbin_errors;
  pTbin_errors.push_back(0.075);
 
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  
  pTbin_errors.push_back(2.);
  pTbin_errors.push_back(2.);
  pTbin_errors.push_back(2.);
  pTbin_errors.push_back(2.);
  pTbin_errors.push_back(2.);

 
  

  vector<float> Closure_integral_eta0;
  vector<float> Closure_integral_phi0;
  vector<float> Closure_integral_eta1;
  vector<float> Closure_integral_phi1;
  vector<float> Closure_integral_eta2;
  vector<float> Closure_integral_phi2;
  vector<float> Closure_integral_eta3;
  vector<float> Closure_integral_phi3;
 
 
  TGraphErrors *Closure_integral_eta_pT[12][4];
  TGraphErrors *Closure_integral_phi_pT[12][4];


  TGraphErrors *Closure_integral_eta_pT2[12][4];
 

  vector<float> closure_integral_values, closure_integral_errors;

  TH1D *Closure_integral_eta_cent[12][nTrkPtBins];
  TH1D *Closure_integral_eta_cent2[12][nTrkPtBins];

  TLine *lineCent, *linePt;


  TCanvas *cintegral_eta_pT[12];
  TCanvas *cintegral_phi_pT[12];
   
  TCanvas *cintegral_eta_cent[12];
  TCanvas *cintegral_phi_cent[12];


  
  TString in_name, plotname, outname, funcname, centlabel, datalabel, jettype,jettype2, pTlabel;

 

  TF1 *gaus1d = new TF1("gaus1d","[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))");

  TF1 *gaus_phi[12][nTrkPtBins][4];
  TF1 *gaus_eta[12][nTrkPtBins][4];
 
  TLegend *lcheck, *leta, *lHminusP;
  
  TLine *linePhi, *lineEta;
 
  TLegend *l40,*l41,*l42;

  int mc_type_code;

  int llimiteta1, rlimiteta1,llimiteta2, rlimiteta2; 
  Double_t check_ymax, check_ymin, dx_eta, dx_phi, bc, err, evalpt, temp1, err1;

  /////////////////////////

  etalim = 1.;
  philim = 1.;

  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
 
  int gstart = 0; 
  int gend = 8;


  for(int g=gstart; g<gend; g++){


    cout<<"starting "<<g<<endl;
    //  There will only be one big "g-loop".
    
    if(g>1&&g<6)continue;
    switch(g){

    case 0:

      fin[g] = new TFile("../me_correct/Pythia_RecoJet_GenTrack_Inclusive_Correlations.root","READ");
      mc_type_code = 2;
      if(is_number)  jettype = "";
      else jettype = "pTweighted";
      jettype2 = "Inclusive_";
      break;
    case 1:
      fin[g] = new TFile("../me_correct/Pythia_GenJet_GenTrack_Inclusive_Correlations.root","READ");
      mc_type_code = 4;
      if(is_number)  jettype = "";
      else jettype = "pTweighted";

      jettype2 = "Inclusive_";
      if(is_number)  fout[g] = new TFile("Inclusive_Pythia_JFFResiduals.root", "RECREATE");
      else   fout[g] = new TFile("Inclusive_Pythia_JFFResiduals_pTweighted.root", "RECREATE");
      break; 

    case 6:
      fin[g] = new TFile("../me_correct/HydJet_RecoJet_GenTrack_Sube0_Inclusive_Correlations.root","READ");
      
      mc_type_code = 2;
      if(is_number)  jettype = "";
      else jettype = "pTweighted";
      jettype2 = "Inclusive";
      break;
    case 7:
      fin[g] = new TFile("../me_correct/HydJet_GenJet_GenTrack_Sube0_Inclusive_Correlations.root","READ");
      
      mc_type_code = 4;
      if(is_number)  jettype = "";
      else jettype = "pTweighted";

      jettype2 = "Inclusive_";
      if(is_number)   fout[g] = new TFile("Inclusive_Hydjet_JFFResiduals.root", "RECREATE");
      else    fout[g] = new TFile("Inclusive_Hydjet_JFFResiduals_pTweighted.root", "RECREATE");
      break;
 
    default: 
      cout<<"Invalid input code for inclusive studies<<endl"<<endl;
      return -1;
      break;
    }


    cout<<"got input files"<<endl;
    //----------------------------------------------------
    //  Start of main i & j loops 
    //-----------------------------------------------------

    for(int i=0; i<nTrkPtBins; i++){
	
      if(g==7){

	corr_canvas_eta[g][i] = new TCanvas(Form("CorrCanvasEta%d%d",g,i)," ",10,10,1500,800);
	corr_canvas_eta[g][i]->Divide(4,2,0.,0.);
  
	corr_canvas_phi[g][i] = new TCanvas(Form("CorrCanvasPhi%d%d",g,i)," ",10,10,1500,800);
	corr_canvas_phi[g][i]->Divide(4,2,0.,0.);
      }
    


      for (int j=0; j<4; j++){

	if(g<6&&j>0)continue;
       

	TString in_name = "Yield_BkgSub_"; in_name+=jettype;in_name+= CBin_strs[j]; in_name+="_"; in_name+= CBin_strs[j+1]; in_name+= "_Pt100_Pt300_"; in_name+=TrkPtBin_strs[i]; in_name+="_"; in_name+=TrkPtBin_strs[i+1];

	cout<<in_name<<endl;
	/*
	if(g==7){
	  in_name = "Raw_Yield_"; in_name+= CBin_strs[j]; in_name+="_"; in_name+= CBin_strs[j+1]; in_name+= "_Pt100_Pt300_"; in_name+=TrkPtBin_strs[i]; in_name+="_"; in_name+=TrkPtBin_strs[i+1];
	}
	if(g==6){	
	  in_name = "Raw_Yield_"; in_name+= CBin_strs[j]; in_name+="_"; in_name+= CBin_strs[j+1]; in_name+= "_Pt100_Pt300_"; in_name+=TrkPtBin_strs[i]; in_name+="_"; in_name+=TrkPtBin_strs[i+1];
	}
	*/

	cout<<g<<" "<<i<<" "<<j<<" "<<in_name<<endl;
	result[g][i][j] = (TH2D*)fin[g]->Get(in_name)->Clone(in_name);


	if(i>3){
	  result[g][i][j]->Scale(1./4.);
	}else if(i==0){
	  result[g][i][j]->Scale(1/.3);
	}

	//-------------------------------
	//dEta projection
	//------------------------

	TString eta_proj_name= in_name;
	eta_proj_name.ReplaceAll("Yield_BkgSub","Eta_Proj");
	eta_proj_name.ReplaceAll("hJetTrackSignalBackground","Eta_Proj");
	eta_proj_name.ReplaceAll("Yield","Eta_Proj");
	    
	llimiteta = result[g][i][j]->GetXaxis()->FindBin(-etalim+.001);
	rlimiteta = result[g][i][j]->GetXaxis()->FindBin(etalim-.001);

	llimitphi = result[g][i][j]->GetYaxis()->FindBin(-philim+.001);
	rlimitphi = result[g][i][j]->GetYaxis()->FindBin(philim-.001);
	    

	eta_proj[g][i][j] = result[g][i][j]->ProjectionX(eta_proj_name,llimitphi,rlimitphi);
	dx_eta = eta_proj[g][i][j]->GetBinWidth(1);
	eta_proj[g][i][j]->Scale(1/dx_eta);

	//	  eta_proj[g][i][j]->Scale(0.8);
	  

	TString eta_proj_name_rebin = eta_proj_name;
	eta_proj_name_rebin.ReplaceAll("Eta_Proj","Eta_Proj_Rebin");
	 
	eta_proj_rebin[g][i][j] = (TH1D*)Rebin_dEta(eta_proj[g][i][j]);
	eta_proj_rebin[g][i][j]->SetName(eta_proj_name_rebin);


	//-------------------------------
	//dPhi projection
	//------------------------

	TString phi_proj_name= in_name;
	phi_proj_name.ReplaceAll("Yield_BkgSub","Phi_Proj");
	phi_proj_name.ReplaceAll("hJetTrackSignalBackground","Phi_Proj");
	phi_proj_name.ReplaceAll("Yield","Phi_Proj");

	phi_proj[g][i][j] = result[g][i][j]->ProjectionY(phi_proj_name,llimiteta,rlimiteta);
	dx_phi = phi_proj[g][i][j]->GetBinWidth(1);
	phi_proj[g][i][j]->Scale(1/dx_phi);

	TString phi_proj_name_rebin = phi_proj_name;
	phi_proj_name_rebin.ReplaceAll("Phi_Proj","Phi_Proj_Rebin");
	 
	phi_proj_rebin[g][i][j] = (TH1D*)Rebin_dPhi(phi_proj[g][i][j]);
	phi_proj_rebin[g][i][j]->SetName(phi_proj_name_rebin);

	float totbins = eta_proj_rebin[g][i][j]->GetNbinsX();

	//	  offset= (eta_proj_rebin[g][i][j]->GetBinContent(1)+eta_proj_rebin[g][i][j]->GetBinContent(2)+eta_proj_rebin[g][i][j]->GetBinContent(3)+eta_proj_rebin[g][i][j]->GetBinContent(totbins-2)+eta_proj_rebin[g][i][j]->GetBinContent(totbins-1)+eta_proj_rebin[g][i][j]->GetBinContent(totbins))/6.;

	offset= (eta_proj_rebin[g][i][j]->GetBinContent(eta_proj_rebin[g][i][j]->FindBin(1.001))+eta_proj_rebin[g][i][j]->GetBinContent(eta_proj_rebin[g][i][j]->FindBin(-1.01)))/2.;

	do_offset->SetParameter(0, offset);

	//	eta_proj_rebin[g][i][j]->Add(do_offset);
	//	phi_proj_rebin[g][i][j]->Add(do_offset);

	check_ymax  = 4.;
	check_ymin = -1.;

	switch(i){
	case 0: 
	  check_ymax = 11.;
	  check_ymin = -1.75; 
	  break;
	case 1: 
	  check_ymax = 6.5; 
	  check_ymin = -1.25; 
	  break;
	case 2: 
	  check_ymax = 4.2; 
	  check_ymin = -1.05; 
	  break;
	case 3: 
	  check_ymax = 3.2; 
	  check_ymin = -.75; 
	   
	  break;
	case 4: 
	  check_ymax = 21;
	  check_ymin = -1.; 
	  break;
	default: 
	  check_ymax = 10.;
	  check_ymin = -1.; 
	  break;
	}
	
	eta_proj_rebin[g][i][j]->SetLineColor(kBlack);
	eta_proj_rebin[g][i][j]->SetMarkerColor(kBlack);
	eta_proj_rebin[g][i][j]->SetMarkerStyle(20);
	eta_proj_rebin[g][i][j]->SetMarkerSize(1);
	eta_proj_rebin[g][i][j]->SetMinimum(check_ymin);
	eta_proj_rebin[g][i][j]->SetMaximum(check_ymax);

	phi_proj_rebin[g][i][j]->SetLineColor(kBlack);
	phi_proj_rebin[g][i][j]->SetMarkerColor(kBlack);
	phi_proj_rebin[g][i][j]->SetMarkerStyle(20);
	phi_proj_rebin[g][i][j]->SetMarkerSize(1);
	phi_proj_rebin[g][i][j]->SetMinimum(check_ymin);
	phi_proj_rebin[g][i][j]->SetMaximum(check_ymax);

	if(g<6){
	  phi_proj_rebin[g][i][j]->SetLineColor(kRed);
	  phi_proj_rebin[g][i][j]->SetMarkerColor(kRed);
	  eta_proj_rebin[g][i][j]->SetLineColor(kRed);
	  eta_proj_rebin[g][i][j]->SetMarkerColor(kRed);
	}


	//-------------------

	//   Saving & Plotting!

	//--------------------

	if(g%2==0)continue;

	TString residual_name_eta = in_name; 
	residual_name_eta.ReplaceAll("Yield_BkgSub","JFF_Residual_Eta");
	residual_name_eta.ReplaceAll("Yield","JFF_Residual_Eta");

	jff_residual_eta[g][i][j] = (TH1D*)eta_proj_rebin[g-1][i][j]->Clone(residual_name_eta);
	jff_residual_eta[g][i][j]->Add(eta_proj_rebin[g][i][j],-1.);
	jff_residual_eta[g][i][j]->SetMinimum(check_ymin-1.);
	jff_residual_eta[g][i][j]->SetMaximum(check_ymax-1.);
	if(i==4){
	  jff_residual_eta[g][i][j]->SetMinimum(-10.);
	  jff_residual_eta[g][i][j]->SetMaximum(11.);
	}
	
	eta_proj_rebin[g-1][i][j]->Write();
	eta_proj_rebin[g][i][j]->Write();
	jff_residual_eta[g][i][j]->Write();



	TString residual_name_phi = in_name; 
	residual_name_phi.ReplaceAll("Yield_BkgSub","JFF_Residual_Phi");
	residual_name_phi.ReplaceAll("Yield","JFF_Residual_Phi");

	jff_residual_phi[g][i][j] = (TH1D*)phi_proj_rebin[g-1][i][j]->Clone(residual_name_phi);
	jff_residual_phi[g][i][j]->Add(phi_proj_rebin[g][i][j],-1.);
	jff_residual_phi[g][i][j]->SetMinimum(check_ymin-1.);
	jff_residual_phi[g][i][j]->SetMaximum(check_ymax-1.);
	
	if(i==4){
	  jff_residual_phi[g][i][j]->SetMinimum(-10.);
	  jff_residual_phi[g][i][j]->SetMaximum(11.);
	}


	phi_proj_rebin[g-1][i][j]->Write();
	phi_proj_rebin[g][i][j]->Write();
	jff_residual_phi[g][i][j]->Write();

      	      
	if((i==0&&j==0)||(g%2!=0&&i==0)){
	  Closure_integral_eta0.clear();
	  Closure_integral_phi0.clear();
	  Closure_integral_eta1.clear();
	  Closure_integral_phi1.clear();
	  Closure_integral_eta2.clear();
	  Closure_integral_phi2.clear();
	  Closure_integral_eta3.clear();
	  Closure_integral_phi3.clear();
	}

	llimiteta = jff_residual_eta[g][i][j]->GetXaxis()->FindBin(-1.0+.0001);
	rlimiteta = jff_residual_eta[g][i][j]->GetXaxis()->FindBin(1.0-.0001);
	     
	double Yield_eta = jff_residual_eta[g][i][j]->Integral(llimiteta,rlimiteta,"width");	      

	//	    if(Yield_eta<0.){Yield_eta=0.;}

	switch(j){
	case 0:
	  Closure_integral_eta0.push_back(Yield_eta);
	  Closure_integral_phi0.push_back(Yield_eta);
	  break;
	case 1:
	  Closure_integral_eta1.push_back(Yield_eta);
	  Closure_integral_phi1.push_back(Yield_eta);
	  break;
	case 2:
	  Closure_integral_eta2.push_back(Yield_eta);
	  Closure_integral_phi2.push_back(Yield_eta);
	  break;
	case 3:
	  Closure_integral_eta3.push_back(Yield_eta);
	  Closure_integral_phi3.push_back(Yield_eta);
	  break;
	}
      

	cout<<g<<" "<<i<<" "<<j<<" "<<Yield_eta<<endl;

	if(g==7){
	
	  corr_canvas_eta[g][i]->cd(4-j);


	  eta_proj_rebin[g][i][j]->SetMarkerStyle(4);
	  if(j==3) eta_proj_rebin[g][i][j]->GetYaxis()->SetLabelSize(0.06);
	  else eta_proj_rebin[g][i][j]->GetYaxis()->SetLabelSize(0.0);
	  eta_proj_rebin[g][i][j]->Draw();



	  
	  eta_proj_rebin[g-1][i][j]->Draw("same");

	  eta_proj_rebin[g-7][i][0]->Draw("same");
	  eta_proj_rebin[g-6][i][0]->SetMarkerStyle(4);
	  eta_proj_rebin[g-6][i][0]->Draw("same");

	  //	  drawlabels(g,i,j);


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

	  TLine *l_eta = new TLine(-1.5,0.,1.5,0.);
	  l_eta->SetLineStyle(2);
	  l_eta->Draw();
	  
	  corr_canvas_eta[g][i]->cd(8-j);

	  if(j==3) jff_residual_eta[g][i][j]->GetYaxis()->SetLabelSize(0.06);
	  else jff_residual_eta[g][i][j]->GetYaxis()->SetLabelSize(0.0);
	  jff_residual_eta[g][i][j]->GetXaxis()->SetLabelSize(0.06);
	  jff_residual_eta[g][i][j]->GetXaxis()->SetTitleSize(0.06);
	  jff_residual_eta[g][i][j]->GetXaxis()->SetTitle("#Delta#eta");
	  jff_residual_eta[g][i][j]->GetXaxis()->CenterTitle();
	  
	  jff_residual_eta[g][i][j]->Draw();
	  jff_residual_eta[g-6][i][0]->Draw("same");
	  l_eta->Draw();
	  
	  TLegend *legend = new TLegend(0.2,0.75,0.9,0.95);
	  legend->AddEntry( eta_proj_rebin[g-1][i][j],"P+H RecoGen");
	  legend->AddEntry( eta_proj_rebin[g][i][j],"P+H GenGen");
	  legend->AddEntry( eta_proj_rebin[g-7][i][0],"Pythia RecoGen");
	  legend->AddEntry( eta_proj_rebin[g-6][i][0],"Pythia GenGen");

	  legend->SetTextSize(0.05);
	  legend->SetLineColor(kWhite);

	  if(j==3)legend->Draw();

					


	  corr_canvas_phi[g][i]->cd(4-j);

	  phi_proj_rebin[g][i][j]->SetMarkerStyle(4);
	  if(j==3) phi_proj_rebin[g][i][j]->GetYaxis()->SetLabelSize(0.06);
	  else phi_proj_rebin[g][i][j]->GetYaxis()->SetLabelSize(0.0);
	  phi_proj_rebin[g][i][j]->Draw();


	  
	  phi_proj_rebin[g-1][i][j]->Draw("same");

	  phi_proj_rebin[g-7][i][0]->Draw("same");
	  phi_proj_rebin[g-6][i][0]->SetMarkerStyle(4);
	  phi_proj_rebin[g-6][i][0]->Draw("same");

	  //	  drawlabels(g,i,j);

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

	  TLine *l_phi = new TLine(-1.5,0.,1.5,0.);
	  l_phi->SetLineStyle(2);
	  l_phi->Draw();
	  
	  corr_canvas_phi[g][i]->cd(8-j);

	  if(j==3) jff_residual_phi[g][i][j]->GetYaxis()->SetLabelSize(0.06);
	  else jff_residual_phi[g][i][j]->GetYaxis()->SetLabelSize(0.0);
	  jff_residual_phi[g][i][j]->GetXaxis()->SetLabelSize(0.06);
	  jff_residual_phi[g][i][j]->GetXaxis()->SetTitleSize(0.06);
	  jff_residual_phi[g][i][j]->GetXaxis()->SetTitle("#Delta#phi");
	  jff_residual_phi[g][i][j]->GetXaxis()->CenterTitle();
	  
	  jff_residual_phi[g][i][j]->Draw();
	  jff_residual_phi[g-6][i][0]->Draw("same");
	  l_phi->Draw();
	
	  if(j==3)legend->Draw();

	}
     
      }
 
      if(g==7){

	TString save_name_eta = "JFF_Residual_Corrections_Eta_";
	save_name_eta+=jettype2;
	save_name_eta+=jettype;
	save_name_eta+=TrkPtBin_strs[i]; save_name_eta+="_"; save_name_eta+=TrkPtBin_strs[i+1];
	save_name_eta+=".png";
	corr_canvas_eta[g][i]->SaveAs(save_name_eta);
	save_name_eta.ReplaceAll(".png",".pdf");
	corr_canvas_eta[g][i]->SaveAs(save_name_eta);

	TString save_name_phi = "JFF_Residual_Corrections_Phi_";
	save_name_phi+=jettype2;
	save_name_phi+=jettype;
	save_name_phi+=TrkPtBin_strs[i]; save_name_phi+="_"; save_name_phi+=TrkPtBin_strs[i+1];
	save_name_phi+=".png";
	corr_canvas_phi[g][i]->SaveAs(save_name_phi);
	save_name_phi.ReplaceAll(".png",".pdf");
	corr_canvas_phi[g][i]->SaveAs(save_name_phi);

      }

    }//i
  

    if(is_number){

      cout<<"starting integrals"<<endl;
      TString integral_eta_pT_name = "integral_eta_pT";
      integral_eta_pT_name+=g;

      cintegral_eta_pT[g] = new TCanvas(integral_eta_pT_name,"",10,10,1500,500);
      cintegral_eta_pT[g]->Divide(4,1,0.,0.);


      TString integral_eta_cent_name = "integral_eta_cent";
      integral_eta_cent_name+=g;

      cintegral_eta_cent[g] = new TCanvas(integral_eta_cent_name,"",10,10,3000,500);
      cintegral_eta_cent[g]->Divide(8,1,0.,0.);

  
  
      for(int j = 0; j<4; j++){

	cout<<"here "<<j<<endl;

	if(g<6&&j>0)continue;

	//    in_name = make_name("Result_",g,3,j,0,centlabel,pTlabel);


	cintegral_eta_pT[g]->cd(4-j);

	TString ClosureIntegralEtaPt_name = in_name;
	ClosureIntegralEtaPt_name.ReplaceAll("Yield_BkgSub","Closure_Integral_Eta");
	ClosureIntegralEtaPt_name.ReplaceAll("_TrkPt4_TrkPt8","");
	ClosureIntegralEtaPt_name.ReplaceAll("_Pt100_Pt300",jettype);
	ClosureIntegralEtaPt_name+=g;
 

	cout<<"Sizes: "<<pTbin_centers.size()<<" "<<Closure_integral_eta3.size()<<endl;

	switch(j){
	case 0:
	  Closure_integral_eta_pT[g][j] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta0[0],&pTbin_errors[0],&closure_integral_errors[0]);
	  break;
	case 1:
	  Closure_integral_eta_pT[g][j] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta1[0],&pTbin_errors[0],&closure_integral_errors[0]);
	  break;
	case 2:
	  Closure_integral_eta_pT[g][j] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta2[0],&pTbin_errors[0],&closure_integral_errors[0]);
	  break;
	case 3:
	  Closure_integral_eta_pT[g][j] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta3[0],&pTbin_errors[0],&closure_integral_errors[0]);
	  break;

	}

	Closure_integral_eta_pT[g][j]->SetName(ClosureIntegralEtaPt_name);
	cout<<g<<ClosureIntegralEtaPt_name<<endl;
      

	Closure_integral_eta_pT[g][j]->SetMarkerColor(1);
	Closure_integral_eta_pT[g][j]->SetMarkerSize(1);
	Closure_integral_eta_pT[g][j]->SetLineColor(1);
	Closure_integral_eta_pT[g][j]->SetMarkerStyle(10);


	Closure_integral_eta_pT[g][j]->SetMinimum(-2.);
	Closure_integral_eta_pT[g][j]->SetMaximum(3.);
	    
	Closure_integral_eta_pT[g][j]->GetXaxis()->SetRangeUser(.5,20.);
	Closure_integral_eta_pT[g][j]->GetYaxis()->SetNdivisions(306);
	Closure_integral_eta_pT[g][j]->Draw("p X A");
	 

	Closure_integral_eta_pT[g][j]->GetYaxis()->SetLabelSize(0.07);
	   


	Closure_integral_eta_pT[g][j]->GetXaxis()->SetTitle("Track p_{T} (GeV/c)");
	Closure_integral_eta_pT[g][j]->GetXaxis()->SetTitleSize(0.06);
	Closure_integral_eta_pT[g][j]->GetXaxis()->SetTitleOffset(xoffset+0.2);
	Closure_integral_eta_pT[g][j]->GetYaxis()->SetTitle("(dN/dp_{T})_{P+H} - (dN/dp_{T})_{PYTH} (GeV/c)^{-1}");

	Closure_integral_eta_pT[g][j]->GetXaxis()->SetNdivisions(8);
   
	
	Closure_integral_eta_pT[g][j]->GetXaxis()->CenterTitle();
	Closure_integral_eta_pT[g][j]->GetYaxis()->CenterTitle();
	   
	if(j<3){
	  Closure_integral_eta_pT[g][j]->GetYaxis()->SetTitleSize(0.0);
	  Closure_integral_eta_pT[g][j]->GetYaxis()->SetLabelSize(0.0);
	  Closure_integral_eta_pT[g][j]->GetXaxis()->SetTitleSize(0.07);
	  Closure_integral_eta_pT[g][j]->GetXaxis()->SetLabelSize(0.07);
	  Closure_integral_eta_pT[g][j]->GetXaxis()->SetTitleOffset(xoffset+0.15);
	}else{
	  Closure_integral_eta_pT[g][j]->GetXaxis()->SetLabelSize(ts3);
	  Closure_integral_eta_pT[g][j]->GetXaxis()->SetLabelOffset(0.015);
	  Closure_integral_eta_pT[g][j]->GetYaxis()->SetTitleOffset(1.);
	  Closure_integral_eta_pT[g][j]->GetYaxis()->SetTitleSize(0.06);
	  Closure_integral_eta_pT[g][j]->GetYaxis()->SetLabelSize(0.06);
	}



	Closure_integral_eta_pT[g][j]->SetMarkerSize(2);



	linePt = new TLine(.5,0,20.,0);
	linePt->SetLineStyle(2);
	linePt->SetLineWidth(1);
	linePt->Draw("same");

	cout<<"here "<<g<<endl;
    
	if(g!=7)continue;

	
	closure_integral_values.clear();
	closure_integral_errors.clear();
	
	for(int k = 0; k<4; k++){
	  double pt_val, x_val;
	  
	  Closure_integral_eta_pT[g][j]->GetPoint(k,x_val,pt_val);
	  closure_integral_values.push_back(pt_val);
	  closure_integral_errors.push_back(pt_val/2.);

	  cout<<pt_val<<endl;
	  
	}


      
	 	  
	Closure_integral_eta_pT2[g][j] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&closure_integral_values[0],&pTbin_errors[0],&closure_integral_errors[0]);

  
	Closure_integral_eta_pT[1][0]->SetMarkerColor(kRed);
	Closure_integral_eta_pT[1][0]->SetLineColor(kRed);
	Closure_integral_eta_pT[1][0]->Draw("same p X");


	if(g==7&&j==3){ 
	  l40 = new TLegend(0.2,0.7,0.8,0.8);
	  l40->SetName("l40");
	  l40->SetTextFont(43);
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);


	  l40->AddEntry(Closure_integral_eta_pT[7][j],"Inclusive P+H","p");

	  l40->AddEntry(Closure_integral_eta_pT[1][0],"Inclusive PYTHIA","p");

	  l40->Draw("same");

	}
      
	//     drawlabels_int_pt2(g,j);
      
	if(j==0){
	  labels = new TPaveText(0.18,0.85,0.45,0.95,"NDC");
	}else{
	  labels = new TPaveText(0.05,0.85,0.45,0.95,"NDC");
	}  
	labels->SetName("labels");
	labels->SetFillColor(0);
	labels->SetLineColor(0);
	labels->SetTextAlign(11);
	labels->AddText(CBin_labels[j]);
	labels->SetTextSize(0.06);
	labels->Draw("same");

      
      }
      cout<<"and here"<<endl;

      if(g%2==0)continue;
  
      for(int i = 0; i<nTrkPtBins; i++){
	
	TString ClosureIntegralEtaCent_name = "ClosureIntegralEtaCent";
	ClosureIntegralEtaCent_name+=g;
	ClosureIntegralEtaCent_name+=i;
	Closure_integral_eta_cent[g][i] = new TH1D(ClosureIntegralEtaCent_name,"",4,xAxis);

	for(int k=0; k<4; k++){
	  evalpt = pTbin_centers.at(i);
	  if(g<6) value = Closure_integral_eta_pT[g][0]->Eval(evalpt);
	  else	value = Closure_integral_eta_pT[g][k]->Eval(evalpt);
	  Closure_integral_eta_cent[g][i]->SetBinContent(4-k,value);
	}
  
	switch(g){
	case 1: 
	  Closure_integral_eta_cent[g][i]->SetMarkerStyle(10);
	  break;
 
	case 7: 
	  Closure_integral_eta_cent[g][i]->SetMarkerStyle(10);
	  break;
  
	default:
	  Closure_integral_eta_cent[g][i]->SetMarkerStyle(10);
	  break;
	}


	Closure_integral_eta_cent[g][i]->SetMarkerSize(2);
	Closure_integral_eta_cent[g][i]->SetLineColor(kBlack);

	Closure_integral_eta_cent[g][i]->SetLineColor(kBlack);
	Closure_integral_eta_cent[g][i]->SetLineColor(kBlack);
	Closure_integral_eta_cent[g][i]->GetYaxis()->SetNdivisions(306);

   
	if(g==7){


	  cintegral_eta_cent[g]->cd(i+1);


	  TString histnameblank = "blank_hist";
	  histnameblank+=g;
	  histnameblank+=i;

	  blank[i] = new TH1D(histnameblank,"",4,xAxis);

	
	  TString histnameblank2 = "blank_hist2";
	  histnameblank2+=g;
	  histnameblank2+=i;

	
	  blank[i]->SetMinimum(-2.);
	  blank[i]->SetMaximum(3.);
	  blank[i]->GetXaxis()->SetTitle("Centrality (%)");
	  blank[i]->GetXaxis()->SetTitleOffset(1.1);
	  blank[i]->GetXaxis()->CenterTitle(true);
	  blank[i]->GetXaxis()->SetTitleSize(0.07);

	  blank[i]->GetYaxis()->SetTitle("(dN/dp_{T})_{PbPb}- (dN/dp_{T})_{pp} (GeV/c)^{-1}");
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
	    blank[i]->GetYaxis()->SetTitleSize(0.07);
	    blank[i]->GetXaxis()->SetTitleOffset(1.1);
	    blank[i]->GetXaxis()->SetTitleSize(0.07);
	    blank[i]->SetLabelSize(0.95*blank[i]->GetXaxis()->GetLabelSize());
	    blank[i]->GetYaxis()->SetLabelSize(0.06);
	    break;
	  case 3:
	    // gPad->SetRightMargin(0.02);
	    break;
	  default:
	    break;
	  }

	  //----------------------------------

	  blank[i]->GetXaxis()->SetTitle("Centrality (%)");
	  blank[i]->GetXaxis()->SetTitleOffset(1.1);
	  blank[i]->GetXaxis()->CenterTitle(true);
	  blank[i]->GetXaxis()->SetTitleSize(0.07);

	  blank[i]->GetYaxis()->SetTitle("(dN/dp_{T})_{PbPb}- (dN/dp_{T})_{pp} (GeV/c)^{-1}");
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
	    blank[i]->GetYaxis()->SetTitleSize(0.07);
	    blank[i]->GetXaxis()->SetTitleOffset(1.1);
	    blank[i]->GetXaxis()->SetTitleSize(0.07);
	    blank[i]->SetLabelSize(0.95*blank[i]->GetXaxis()->GetLabelSize());
	    blank[i]->GetYaxis()->SetLabelSize(0.06);
	    break;
	  default:
	    blank[i]->GetYaxis()->SetLabelSize(0.);
	    break;
	  }



	  blank[i]->Draw();

	  Closure_integral_eta_cent[1][i]->SetMarkerColor(kRed);
	  Closure_integral_eta_cent[1][i]->Draw("same p X");

	  Closure_integral_eta_cent[7][i]->SetMarkerColor(kBlack);
	  Closure_integral_eta_cent[7][i]->Draw("same p X");

	  Closure_integral_eta_cent[1][i]->Write();

	  TLatex *pt_label = new TLatex(0.2,0.95, TrkPtBin_labels[i]);
	  pt_label->SetNDC();
	  pt_label->SetTextSizePixels(tspixels);
	  pt_label->Draw();

	  if(i==0)	l40->Draw("same");

	}

      }
  
   			       
      cintegral_eta_pT[g]->cd(0);
								      
      TLatex *canvas_title = new TLatex(0.06,0.9,"CMS Preliminary Simulation");
      canvas_title->SetTextSizePixels(tspixels);
      canvas_title->SetTextFont(63);
      canvas_title->Draw();

      TLatex *canvas_title2 = new TLatex(0.295,0.9,"PYTHIA+HYDJET");
      canvas_title2->SetTextSizePixels(tspixels);
      canvas_title2->Draw();

      cintegral_eta_pT[g]->SaveAs("Integral_JFFResidual_pT.pdf");
      cintegral_eta_pT[g]->SaveAs("Integral_JFFResidual_pT.png");


      cintegral_eta_cent[g]->cd(0);

      canvas_title->Draw();

      canvas_title2->Draw();

   

      cintegral_eta_cent[g]->SaveAs("Integral_JFFResidual_Cent.png");
      cintegral_eta_cent[g]->SaveAs("Integral_JFFResidual_Cent.pdf");
    }

  }//g
   
  return 0;
}
