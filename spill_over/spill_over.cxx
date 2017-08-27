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


Int_t spill_over(bool is_number = kTRUE, bool quark_gluon = kFALSE, bool rx_plane = kFALSE, bool show_errs = kFALSE){



  if(quark_gluon&&rx_plane){

    cerr<<"this code is only good for one comparison at a time"<<endl;
    return -1;

  }


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

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins+1] = {0.7, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8","TrkPt12","TrkPt16","TrkPt20","TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12","12<pT<16","16<pT<20","pT>20"};

  float x, offset, value;
  TF1 *do_offset = new TF1("do_offset","1.*[0]+x-x",-3.,3.);

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
  TH1D *phi_proj_syst[12][nTrkPtBins][4];

  TH1D *eta_proj[12][nTrkPtBins][4];
  TH1D *eta_proj_rebin[12][nTrkPtBins][4];
  TH1D *eta_proj_syst[12][nTrkPtBins][4];


  TH1D *eta_proj_ref[12][nTrkPtBins][4];
  TH1D *phi_proj_ref[12][nTrkPtBins][4];

  TH1D *eta_spill_over[12][nTrkPtBins][4];
  TH1D *phi_spill_over[12][nTrkPtBins][4];



  TCanvas *corr_canvas_eta[12];
  TCanvas *corr_canvas_phi[12];


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
  pTbin_errors.push_back(0.15);
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
 

  vector<float> closure_integral_values, closure_integral_errors_all, closure_integral_errors0,  closure_integral_errors1,  closure_integral_errors2,  closure_integral_errors3;

  TH1D *Closure_integral_eta_cent[12][nTrkPtBins];
  TH1D *Closure_integral_eta_cent2[12][nTrkPtBins];

  TLine *lineCent, *linePt;


  TCanvas *cintegral_eta_pT[12];
  TCanvas *cintegral_phi_pT[12];
   
  TCanvas *cintegral_eta_cent[12];
  TCanvas *cintegral_phi_cent[12];


  float me_err[12][6][5];
  float bg_err[12][6][5];

  
  TString in_name, plotname, outname, funcname, centlabel, datalabel, jettype,jettype2, pTlabel;

 

  TF1 *gaus1d = new TF1("gaus1d","[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))");

  TF1 *gaus_phi[12][nTrkPtBins][4];
  TF1 *gaus_eta[12][nTrkPtBins][4];

  TF1 *gaus_phi_up[12][nTrkPtBins][4];
  TF1 *gaus_eta_up[12][nTrkPtBins][4];
  
  TF1 *gaus_phi_down[12][nTrkPtBins][4];
  TF1 *gaus_eta_down[12][nTrkPtBins][4];
 


  TF1 *fit_line_left = new TF1("fit_line_left","[0]");
  TF1 *fit_line_right = new TF1("fit_line_right","[0]");

  TLegend *lcheck, *leta, *lHminusP;
  
  TLine *linePhi, *lineEta;
 
  TLegend *l40,*l41,*l42;

  int llimiteta1, rlimiteta1,llimiteta2, rlimiteta2; 
  Double_t check_ymax, check_ymin, dx_eta, dx_phi, bc, err, evalpt, temp1, err1, err_temp;

  /////////////////////////

  etalim = 1.;
  philim = 1.;

  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
 
  int gstart = 5; 
  int gend = 8;

  TCanvas *dummy = new TCanvas("toy_canvas");


 
  cout<<"got input files"<<endl;
  //----------------------------------------------------
  //  Start of main i & j loops 
  //-----------------------------------------------------

  corr_canvas_eta[7] = new TCanvas(Form("CorrCanvasEta%d",7)," ",10,10,1500,2400);
  corr_canvas_eta[7]->Divide(4,9,0.,0.);

  corr_canvas_phi[7] = new TCanvas(Form("CorrCanvasPhi%d",7)," ",10,10,1500,2400);
  corr_canvas_phi[7]->Divide(4,9,0.,0.);

  if(is_number)  fout[7] = new TFile("Inclusive_Hydjet_SpillOvers.root", "RECREATE");
  else   fout[7] = new TFile("Inclusive_Hydjet_SpillOvers_pTweighted.root", "RECREATE");

  for(int i=0; i<nTrkPtBins; i++){

    for (int j=0; j<4; j++){


      for(int g=gstart; g<gend; g++){

	cout<<"starting "<<g<<endl;
    
	switch(g){
 

	case 5:
	  if(rx_plane)	  fin[g] = new TFile("../me_correct/HydJet_RecoJet_GenTrack_NoSube0_Inclusive_Correlations_InRxPlaneOnly.root","READ");
	  else	  fin[g] = new TFile("../me_correct/HydJet_RecoJet_GenTrack_NoSube0_Inclusive_Correlations_QuarkOnly.root","READ");
	  jettype = "";
	  if(!is_number) jettype = "pTweighted";
	  jettype2 = "Inclusive";
	  break;

	case 6:
	  if(rx_plane)  fin[g] = new TFile("../me_correct/HydJet_RecoJet_GenTrack_NoSube0_Inclusive_Correlations_OutRxPlaneOnly.root","READ");
	  fin[g] = new TFile("../me_correct/HydJet_RecoJet_GenTrack_NoSube0_Inclusive_Correlations_GluonOnly.root","READ");
	  jettype = "";
	  if(!is_number) jettype = "pTweighted";
	  jettype2 = "Inclusive";
	  break;

	case 7:
	  fin[g] = new TFile("../me_correct/HydJet_RecoJet_GenTrack_NoSube0_Inclusive_Correlations.root","READ");
	  jettype = "";
	  if(!is_number) jettype = "pTweighted";
	  jettype2 = "Inclusive";
	 
	  break;
 
	default: 
	  cout<<"Invalid input code for inclusive studies<<endl"<<endl;
	  return -1;
	  break;
	}



	TString in_name = "Yield_BkgSub_"; in_name+=jettype;in_name+= CBin_strs[j]; in_name+="_"; in_name+= CBin_strs[j+1]; in_name+= "_Pt100_Pt1000_"; in_name+=TrkPtBin_strs[i]; in_name+="_"; in_name+=TrkPtBin_strs[i+1];

	cout<<g<<" "<<i<<" "<<j<<" "<<in_name<<endl;
	result[g][i][j] = (TH2D*)fin[g]->Get(in_name)->Clone((TString)(in_name+"_All"));

	if(is_number){
	  if(i>3){
	    result[g][i][j]->Scale(1./4.);
	  }else if(i==0){
	    result[g][i][j]->Scale(1/.3);
	  }

	}
	cout<<"got one"<<endl;


	//-------------------------------
	//dEta projection
	//------------------------

	TString eta_proj_name= in_name;
	eta_proj_name.ReplaceAll("Yield_BkgSub","Eta_Proj");
	eta_proj_name.ReplaceAll("hJetTrackSignalBackground","Eta_Proj");
	eta_proj_name.ReplaceAll("Yield","Eta_Proj");
	    
	llimiteta = result[g][i][j]->GetXaxis()->FindBin(-1.5+.001);
	rlimiteta = result[g][i][j]->GetXaxis()->FindBin(1.5-.001);

	llimitphi = result[g][i][j]->GetYaxis()->FindBin(-TMath::Pi()/2.+.001);
	rlimitphi = result[g][i][j]->GetYaxis()->FindBin(TMath::Pi()/2.-.001);
	    

	eta_proj[g][i][j] = result[g][i][j]->ProjectionX(eta_proj_name,llimitphi,rlimitphi);
	dx_eta = eta_proj[g][i][j]->GetBinWidth(1);
	eta_proj[g][i][j]->Scale(1/dx_eta);
	  

	TString eta_proj_name_rebin = eta_proj_name;
	eta_proj_name_rebin.ReplaceAll("Eta_Proj","Eta_Proj_Rebin");

	 
	eta_proj_rebin[g][i][j] = (TH1D*)Rebin_dEta(eta_proj[g][i][j]);
	eta_proj_rebin[g][i][j]->SetName(eta_proj_name_rebin);


	eta_proj_name_rebin.ReplaceAll("Eta_Proj","Eta_Proj_Syst");
	eta_proj_syst[g][i][j] = (TH1D*)eta_proj_rebin[g][i][j]->Clone(eta_proj_name_rebin);


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

	phi_proj_name_rebin.ReplaceAll("Phi_Proj","Phi_Proj_Syst");
	phi_proj_syst[g][i][j] = (TH1D*)phi_proj_rebin[g][i][j]->Clone(phi_proj_name_rebin);


	float totbins = eta_proj_rebin[g][i][j]->GetNbinsX();
		
	//	offset= -(phi_proj_rebin[g][i][j]->GetBinContent(phi_proj_rebin[g][i][j]->FindBin(1.001))+phi_proj_rebin[g][i][j]->GetBinContent(phi_proj_rebin[g][i][j]->FindBin(-1.01)))/2.;
     

	//	do_offset->SetParameter(0, offset);
	//	eta_proj_rebin[g][i][j]->Add(do_offset);
	//	phi_proj_rebin[g][i][j]->Add(do_offset);
	
	check_ymax  = 4.;
	check_ymin = -1.;


	eta_proj_rebin[g][i][j]->SetMarkerStyle(20);
	eta_proj_rebin[g][i][j]->SetMarkerSize(1);
	phi_proj_rebin[g][i][j]->SetMarkerStyle(20);
	phi_proj_rebin[g][i][j]->SetMarkerSize(1);

	//SYSTEMATIC UNCERTAINTY SET-UP // 

	bg_err[g][i][j] =TMath::Abs( (  eta_proj_rebin[g][i][j] ->GetBinContent(  eta_proj_rebin[g][i][j]->FindBin(-1.51))- eta_proj_rebin[g][i][j] ->GetBinContent(  eta_proj_rebin[g][i][j]->FindBin(-2.01)))+(eta_proj_rebin[g][i][j] ->GetBinContent(  eta_proj_rebin[g][i][j]->FindBin(1.51))- eta_proj_rebin[g][i][j] ->GetBinContent(  eta_proj_rebin[g][i][j]->FindBin(2.01))))/2.;

	float line_left, line_right; 
	/*

	dummy->cd(0);
       	eta_proj_rebin[g][i][j]->Fit("fit_line_left","q","",-2.5,-1.5);

      	line_left = fit_line_left->GetParameter(0);

	eta_proj_rebin[g][i][j]->Fit("fit_line_right","q","",1.5,2.5);
	
	line_right = fit_line_right->GetParameter(0);

	me_err[g][i][j] = TMath::Abs(line_left - line_right);

	*/
	me_err[g][i][j] = 0.;




      } // g loop
      

      cout<<"got everything we need"<<endl;
      
      //-------------------

      //   Saving & Plotting!

      //--------------------
    
      dummy->cd();
      /*
      eta_proj_rebin[1][i][j] = (TH1D*)fin2[6]->Get( (TString)("JFF_Residual_Eta_"+jettype+ CBin_strs[j]+"_"+CBin_strs[j+1] +"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]))->Clone((TString)("JFF_Residual_Eta_"+jettype+ CBin_strs[j]+"_"+CBin_strs[j+1] +"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]+"_Reco"));

      phi_proj_rebin[1][i][j] = (TH1D*)fin2[6]->Get( (TString)("JFF_Residual_Phi_"+jettype+ CBin_strs[j]+"_"+CBin_strs[j+1] +"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]))->Clone((TString)("JFF_Residual_Phi_"+jettype+ CBin_strs[j]+"_"+CBin_strs[j+1] +"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]+"_Reco"));
  
   
      //    eta_proj_rebin[7][i][j]->Add(eta_proj_rebin[6][i][0],-1.);
      // phi_proj_rebin[7][i][j]->Add(phi_proj_rebin[6][i][0],-1.);
 
      eta_proj_rebin[7][i][j]->Add(eta_proj_rebin[1][i][j],-1.);
      phi_proj_rebin[7][i][j]->Add(phi_proj_rebin[1][i][j],-1.);
 
      */
    
      eta_proj_rebin[5][i][j]->SetLineColor(kBlue);
      eta_proj_rebin[5][i][j]->SetMarkerColor(kBlue);

      eta_proj_rebin[6][i][j]->SetLineColor(kRed);
      eta_proj_rebin[6][i][j]->SetMarkerColor(kRed);
      
      eta_proj_rebin[7][i][j]->SetLineColor(kBlack);
      eta_proj_rebin[7][i][j]->SetMarkerColor(kBlack);
   
      eta_proj_rebin[7][i][j]->SetMinimum(check_ymin);
      eta_proj_rebin[7][i][j]->SetMaximum(check_ymax);

      phi_proj_rebin[7][i][j]->SetLineColor(kBlack);
      phi_proj_rebin[7][i][j]->SetMarkerColor(kBlack);

      phi_proj_rebin[5][i][j]->SetLineColor(kBlue);
      phi_proj_rebin[5][i][j]->SetMarkerColor(kBlue);

      phi_proj_rebin[6][i][j]->SetLineColor(kRed);
      phi_proj_rebin[6][i][j]->SetMarkerColor(kRed);

     
      phi_proj_rebin[7][i][j]->SetMinimum(check_ymin);
      phi_proj_rebin[7][i][j]->SetMaximum(check_ymax);
	     

      cout<<eta_proj_rebin[7][i][j]->GetName()<<endl;

      if(j<2&&i<3){
	llimiteta = eta_proj_rebin[7][i][j]->GetXaxis()->FindBin(-1.5+.0001);
	rlimiteta = eta_proj_rebin[7][i][j]->GetXaxis()->FindBin(1.5-.0001);
      }else{
	llimiteta = eta_proj_rebin[7][i][j]->GetXaxis()->FindBin(-1.+.0001);
	rlimiteta = eta_proj_rebin[7][i][j]->GetXaxis()->FindBin(1.-.0001);

      }

      //FIRST VARY AND FIT FOR SYST ERR

      
      do_offset->SetParameter(0,bg_err[7][i][j]+me_err[7][i][j]);

      eta_proj_syst[7][i][j]->Add(do_offset);

      double Yield_eta_up = eta_proj_syst[7][i][j]->Integral(llimiteta,rlimiteta,"width");	      
     	      
      if(Yield_eta_up < 0.) Yield_eta_up = 0.;

      gaus1d->FixParameter(0,0.);
      gaus1d->FixParameter(1,Yield_eta_up);
      //      gaus1d->ReleaseParameter(2);
      gaus1d->SetParameter(2,0.3);


      eta_proj_syst[7][i][j]->Fit("gaus1d","","",-1.5,1.5);

   TString   gaus_eta_name = "Eta_SpillOver_Fit_Up_"; gaus_eta_name+= CBin_strs[j]; gaus_eta_name+="_"; gaus_eta_name+= CBin_strs[j+1]; gaus_eta_name+= "_Pt100_Pt300_"; gaus_eta_name+=TrkPtBin_strs[i]; gaus_eta_name+="_"; gaus_eta_name+=TrkPtBin_strs[i+1];

      gaus_eta_up[7][i][j] = new TF1(gaus_eta_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

      for(int k= 0; k<3; k++){
	
	double temp = gaus1d->GetParameter(k);
	gaus_eta_up[7][i][j]->SetParameter(k, temp);

      }



   
      
      do_offset->SetParameter(0,-2*bg_err[7][i][j]-2*me_err[7][i][j]);

      eta_proj_syst[7][i][j]->Add(do_offset);

      double Yield_eta_down = eta_proj_syst[7][i][j]->Integral(llimiteta,rlimiteta,"width");	      
     	      
      if(Yield_eta_down < 0.) Yield_eta_down = 0.;

      gaus1d->FixParameter(1,Yield_eta_down);
      //      gaus1d->ReleaseParameter(2);
      gaus1d->SetParameter(2,0.3);


      eta_proj_syst[7][i][j]->Fit("gaus1d","","",-1.5,1.5);

      gaus_eta_name = "Eta_SpillOver_Fit_Down_"; gaus_eta_name+= CBin_strs[j]; gaus_eta_name+="_"; gaus_eta_name+= CBin_strs[j+1]; gaus_eta_name+= "_Pt100_Pt300_"; gaus_eta_name+=TrkPtBin_strs[i]; gaus_eta_name+="_"; gaus_eta_name+=TrkPtBin_strs[i+1];

      gaus_eta_down[7][i][j] = new TF1(gaus_eta_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

      for(int k= 0; k<3; k++){
	
	double temp = gaus1d->GetParameter(k);
	gaus_eta_down[7][i][j]->SetParameter(k, temp);

      }

      for(int k = 0; k<eta_proj_syst[7][i][j]->GetNbinsX()+1; k++){

	err_temp = (gaus_eta_up[7][i][j]->Eval(eta_proj_syst[7][i][j]->GetBinCenter(k))-gaus_eta_down[7][i][j]->Eval(eta_proj_syst[7][i][j]->GetBinCenter(k)))/2.;
	eta_proj_syst[7][i][j]->SetBinContent(k,err_temp);
	eta_proj_syst[7][i][j]->SetBinError(k,0.);
      }


      double Yield_eta = eta_proj_rebin[7][i][j]->Integral(llimiteta,rlimiteta,"width");	      
      
      cout<<Yield_eta<<endl;
	     	     
	    
      if(Yield_eta<0.){Yield_eta=0.;}

	      
      gaus1d->FixParameter(1,Yield_eta);
	      
      //      gaus1d->ReleaseParameter(2);

      gaus1d->SetParameter(2,0.3);
 

      if(i<5)  eta_proj_rebin[7][i][j]->Fit("gaus1d","","",-1.5,1.5);
	      
      gaus_eta_name = "Eta_SpillOver_Fit_"; gaus_eta_name+= CBin_strs[j]; gaus_eta_name+="_"; gaus_eta_name+= CBin_strs[j+1]; gaus_eta_name+= "_Pt100_Pt300_"; gaus_eta_name+=TrkPtBin_strs[i]; gaus_eta_name+="_"; gaus_eta_name+=TrkPtBin_strs[i+1];

      gaus_eta[7][i][j] = new TF1(gaus_eta_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

      for(int k= 0; k<3; k++){
	
	double temp = gaus1d->GetParameter(k);
	gaus_eta[7][i][j]->SetParameter(k, temp);

      }

   
      //THEN DO NOMINAL

      err_temp = gaus1d->GetParError(2);

      gaus_eta_name = "Eta_SpillOver_Points_"; gaus_eta_name+= CBin_strs[j]; gaus_eta_name+="_"; gaus_eta_name+= CBin_strs[j+1]; gaus_eta_name+= "_Pt100_Pt300_"; gaus_eta_name+=TrkPtBin_strs[i]; gaus_eta_name+="_"; gaus_eta_name+=TrkPtBin_strs[i+1];

      eta_spill_over[7][i][j] = (TH1D*)eta_proj_rebin[7][i][j]->Clone(gaus_eta_name);
      
      float temp_diff = 0;

      for(int k=0; k< eta_spill_over[7][i][j]->GetNbinsX()+1; k++){
	evalpt = eta_spill_over[7][i][j]->GetBinCenter(k);
	bc = gaus1d->Eval(evalpt);

	temp_diff += TMath::Abs(bc - eta_spill_over[7][i][j]->GetBinContent(k));

	eta_spill_over[7][i][j]->SetBinContent(k,bc);
	eta_spill_over[7][i][j]->SetBinError(k,err_temp);


      }

      cout<<"TEMP DIFF = "<<temp_diff<<" as percent = "<<temp_diff/eta_spill_over[7][i][j]->Integral(1,eta_spill_over[7][i][j]->GetNbinsX()+1);
	
      double err_temp= gaus1d->GetParError(2);

      gaus_eta[7][i][j]->SetLineColor(kBlue);
	


      gaus1d->FixParameter(1,Yield_eta);
	      
      //      gaus1d->ReleaseParameter(2);

      gaus1d->SetParameter(2,0.4);
      gaus1d->SetParLimits(2,0.2,1.);

        //FIRST VARY AND FIT FOR SYST ERR

      
      do_offset->SetParameter(0,bg_err[7][i][j]+me_err[7][i][j]);

      phi_proj_syst[7][i][j]->Add(do_offset);

      double Yield_phi_up = phi_proj_syst[7][i][j]->Integral(llimitphi,rlimitphi,"width");	      
     	      
        gaus1d->FixParameter(1,Yield_phi_up);
       gaus1d->SetParameter(2,0.3);


      phi_proj_syst[7][i][j]->Fit("gaus1d","","",-1.5,1.5);

      TString   gaus_phi_name = "Phi_SpillOver_Fit_Up_"; gaus_phi_name+= CBin_strs[j]; gaus_phi_name+="_"; gaus_phi_name+= CBin_strs[j+1]; gaus_phi_name+= "_Pt100_Pt300_"; gaus_phi_name+=TrkPtBin_strs[i]; gaus_phi_name+="_"; gaus_phi_name+=TrkPtBin_strs[i+1];

      gaus_phi_up[7][i][j] = new TF1(gaus_phi_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

      for(int k= 0; k<3; k++){
	
	double temp = gaus1d->GetParameter(k);
	gaus_phi_up[7][i][j]->SetParameter(k, temp);

      }



   
      
      do_offset->SetParameter(0,-2*bg_err[7][i][j]-2*me_err[7][i][j]);

      phi_proj_syst[7][i][j]->Add(do_offset);

      double Yield_phi_down = phi_proj_syst[7][i][j]->Integral(llimitphi,rlimitphi,"width");	      
     	      
      gaus1d->FixParameter(1,Yield_phi_down);
      //      gaus1d->ReleaseParameter(2);
      gaus1d->SetParameter(2,0.3);


      phi_proj_syst[7][i][j]->Fit("gaus1d","","",-1.5,1.5);

      gaus_phi_name = "Phi_SpillOver_Fit_Down_"; gaus_phi_name+= CBin_strs[j]; gaus_phi_name+="_"; gaus_phi_name+= CBin_strs[j+1]; gaus_phi_name+= "_Pt100_Pt300_"; gaus_phi_name+=TrkPtBin_strs[i]; gaus_phi_name+="_"; gaus_phi_name+=TrkPtBin_strs[i+1];

      gaus_phi_down[7][i][j] = new TF1(gaus_phi_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

      for(int k= 0; k<3; k++){
	
	double temp = gaus1d->GetParameter(k);
	gaus_phi_down[7][i][j]->SetParameter(k, temp);

      }

      for(int k = 0; k<phi_proj_syst[7][i][j]->GetNbinsX()+1; k++){

	err_temp = (gaus_phi_up[7][i][j]->Eval(phi_proj_syst[7][i][j]->GetBinCenter(k))-gaus_phi_down[7][i][j]->Eval(phi_proj_syst[7][i][j]->GetBinCenter(k)))/2.;
	phi_proj_syst[7][i][j]->SetBinContent(k,err_temp);
	phi_proj_syst[7][i][j]->SetBinError(k,0.);
      }


      double Yield_phi = phi_proj_rebin[7][i][j]->Integral(llimitphi,rlimitphi,"width");	      
      
      cout<<Yield_phi<<endl;
	     	     
	    
      if(Yield_phi<0.){Yield_phi=0.;}

	      
      gaus1d->FixParameter(1,Yield_phi);
	      
      //      gaus1d->ReleaseParameter(2);

      gaus1d->SetParameter(2,0.3);
 

      if(i<5)  phi_proj_rebin[7][i][j]->Fit("gaus1d","","",-1.5,1.5);
	      
      gaus_phi_name = "Phi_SpillOver_Fit_"; gaus_phi_name+= CBin_strs[j]; gaus_phi_name+="_"; gaus_phi_name+= CBin_strs[j+1]; gaus_phi_name+= "_Pt100_Pt300_"; gaus_phi_name+=TrkPtBin_strs[i]; gaus_phi_name+="_"; gaus_phi_name+=TrkPtBin_strs[i+1];

      gaus_phi[7][i][j] = new TF1(gaus_phi_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

      for(int k= 0; k<3; k++){
	
	double temp = gaus1d->GetParameter(k);
	gaus_phi[7][i][j]->SetParameter(k, temp);

      }


      if(i<5)  phi_proj_rebin[7][i][j]->Fit("gaus1d","","",-1.5,1.5);



      gaus_phi_name = "Phi_SpillOver_Fit_"; gaus_phi_name+= CBin_strs[j]; gaus_phi_name+="_"; gaus_phi_name+= CBin_strs[j+1]; gaus_phi_name+= "_Pt100_Pt300_"; gaus_phi_name+=TrkPtBin_strs[i]; gaus_phi_name+="_"; gaus_phi_name+=TrkPtBin_strs[i+1];


	    
      gaus_phi[7][i][j] = new TF1(gaus_phi_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);


      for(int k= 0; k<3; k++){
	double temp = gaus1d->GetParameter(k);
	gaus_phi[7][i][j]->SetParameter(k, temp);
      }

      err_temp = gaus1d->GetParError(2);
	 
      gaus_phi[7][i][j]->SetLineColor(kBlue);
      gaus_phi[7][i][j]->Draw("same"); 


      err_temp = gaus1d->GetParError(2);

      gaus_phi_name.ReplaceAll("Fit","Points");

      phi_spill_over[7][i][j] = (TH1D*)phi_proj_rebin[7][i][j]->Clone(gaus_phi_name);
      
      for(int k=0; k< phi_spill_over[7][i][j]->GetNbinsX()+1; k++){
	evalpt = phi_spill_over[7][i][j]->GetBinCenter(k);
	bc = gaus1d->Eval(evalpt);
	phi_spill_over[7][i][j]->SetBinContent(k,bc);
	phi_spill_over[7][i][j]->SetBinError(k,err_temp);
      }

      //WRITE STUFF OUT

      fout[7]->cd();

      eta_spill_over[7][i][j]->Write();
      eta_proj_syst[7][i][j]->Write();

      gaus_eta[7][i][j]->Write();
      gaus_eta_up[7][i][j]->Write();
      gaus_eta_down[7][i][j]->Write();

      phi_spill_over[7][i][j]->Write();
      phi_proj_syst[7][i][j]->Write();
      gaus_phi[7][i][j]->Write();
      gaus_phi_up[7][i][j]->Write();
      gaus_phi_down[7][i][j]->Write();

      /*
      	      
      if(i==0){
	Closure_integral_eta0.clear();
	Closure_integral_phi0.clear();
	Closure_integral_eta1.clear();
	Closure_integral_phi1.clear();
	Closure_integral_eta2.clear();
	Closure_integral_phi2.clear();
	Closure_integral_eta3.clear();
	Closure_integral_phi3.clear();

      }
*/
      //	    if(Yield_eta<0.){Yield_eta=0.;}

      switch(j){
      case 0:
	Closure_integral_eta0.push_back(Yield_eta);
	//	closure_integral_errors0.push_back((Yield_eta_up - Yield_eta_down)/2.);
	closure_integral_errors0.push_back(Yield_eta*.18);
	break;
      case 1:
	Closure_integral_eta1.push_back(Yield_eta);
	closure_integral_errors1.push_back(Yield_eta*.18);
	//closure_integral_errors1.push_back((Yield_eta_up - Yield_eta_down)/2.);
	break;
      case 2:
	Closure_integral_eta2.push_back(Yield_eta);
	closure_integral_errors2.push_back(Yield_eta*.18);
	//	closure_integral_errors2.push_back((Yield_eta_up - Yield_eta_down)/2.);
	break;
      case 3:
	Closure_integral_eta3.push_back(Yield_eta);
	closure_integral_errors3.push_back(Yield_eta*.18);
	//	closure_integral_errors3.push_back((Yield_eta_up - Yield_eta_down)/2.);
	break;
      }
      

      cout<<i<<" "<<j<<" "<<Yield_eta<<" "<<Closure_integral_eta0.size()<<endl;

      corr_canvas_eta[7]->cd(4*(i+1)-j);


      if(j==3) eta_proj_rebin[7][i][j]->GetYaxis()->SetLabelSize(0.06);
      else eta_proj_rebin[7][i][j]->GetYaxis()->SetLabelSize(0.0);
	

      eta_proj_rebin[7][i][j]->GetXaxis()->SetTitleSize(0.06);
      eta_proj_rebin[7][i][j]->GetXaxis()->SetTitle("#Delta#eta");
      eta_proj_rebin[7][i][j]->GetXaxis()->SetLabelSize(0.06);


      eta_proj_rebin[7][i][j]->GetXaxis()->SetRangeUser(-2.49,2.49);
      eta_proj_rebin[7][i][j]->Draw();
  
      if(quark_gluon || rx_plane)   eta_proj_rebin[5][i][j]->Draw("same");
      if(quark_gluon || rx_plane)  eta_proj_rebin[6][i][j]->Draw("same");

      if(show_errs){
	  eta_proj_syst[7][i][j]->SetLineColor(kYellow);
	  eta_proj_syst[7][i][j]->SetLineWidth(4);

	  eta_proj_syst[7][i][j]->Draw("C same");
	  gaus_eta_up[7][i][j]->Draw("same");
	  gaus_eta_down[7][i][j]->Draw("same");
	}
      eta_proj_rebin[7][i][j]->Draw("same");



	  
      //  eta_proj_rebin[6][i][0]->Draw("same");

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

      gaus_eta[7][i][j]->Draw("same");
	  
      TLegend *legend = new TLegend(0.2,0.5,0.9,0.7);
      //      legend->AddEntry( eta_proj_rebin[6][i][0],"Pythia+Hydjet RecoGen - GenGen");
      legend->AddEntry( eta_proj_rebin[7][i][j],"Nominal RecoGen Sube>0");
      if(quark_gluon)   legend->AddEntry( eta_proj_rebin[5][i][j],"Quark Only");
      if(quark_gluon)      legend->AddEntry( eta_proj_rebin[6][i][j],"Gluon Only");

      if(rx_plane)   legend->AddEntry( eta_proj_rebin[5][i][j],"In Rx-Plane Jets");
      if(rx_plane)      legend->AddEntry( eta_proj_rebin[6][i][j],"Out Rx-Plane Jets");

      legend->SetTextSize(0.06);
      legend->SetLineColor(kWhite);

      if(j==3&&i==0)legend->Draw();


      corr_canvas_phi[7]->cd(4*(i+1)-j);

      if(j==3) phi_proj_rebin[7][i][j]->GetYaxis()->SetLabelSize(0.06);
      else phi_proj_rebin[7][i][j]->GetYaxis()->SetLabelSize(0.0);



      phi_proj_rebin[7][i][j]->GetXaxis()->SetTitleSize(0.06);
      phi_proj_rebin[7][i][j]->GetXaxis()->SetTitle("#Delta#eta");
      phi_proj_rebin[7][i][j]->GetXaxis()->SetLabelSize(0.06);

      phi_proj_rebin[7][i][j]->Draw();

      if(quark_gluon || rx_plane)   phi_proj_rebin[5][i][j]->Draw("same");
      if(quark_gluon || rx_plane)   phi_proj_rebin[6][i][j]->Draw("same");



	if(show_errs){
	  phi_proj_syst[7][i][j]->SetLineColor(kYellow);
	  phi_proj_syst[7][i][j]->SetLineWidth(4);

	  phi_proj_syst[7][i][j]->Draw("C same");
	  gaus_phi_up[7][i][j]->Draw("same");
	  gaus_phi_down[7][i][j]->Draw("same");
	}

      phi_proj_rebin[7][i][j]->Draw("same");


	  
      //   phi_proj_rebin[6][i][0]->Draw("same");

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
	  
      gaus_phi[7][i][j]->Draw("same");
	
      if(j==3&&i==0)legend->Draw();

    }
    dummy->cd();   
       


    
  }//i
  
  TString save_name_eta = "SpillOver_Corrections_Eta";
  if(!is_number) save_name_eta+="_";
  save_name_eta+=jettype;

  if(quark_gluon)   save_name_eta+="_QG_Comparison";
  if(rx_plane)   save_name_eta+="_RXplane_Comparison";
  //   save_name_eta+=TrkPtBin_strs[i]; save_name_eta+="_"; save_name_eta+=TrkPtBin_strs[i+1];
  save_name_eta+=".png";
  corr_canvas_eta[7]->SaveAs(save_name_eta);
  save_name_eta.ReplaceAll(".png",".pdf");
  corr_canvas_eta[7]->SaveAs(save_name_eta);

  TString save_name_phi = "SpillOver_Corrections_Phi";
  if(!is_number) save_name_phi+="_";
  save_name_phi+=jettype;

  if(quark_gluon)   save_name_phi+="_QG_Comparison";
  if(rx_plane)   save_name_phi+="_RXplane_Comparison";

  //    save_name_phi+=TrkPtBin_strs[i]; save_name_phi+="_"; save_name_phi+=TrkPtBin_strs[i+1];
  save_name_phi+=".png";
  corr_canvas_phi[7]->SaveAs(save_name_phi);
  save_name_phi.ReplaceAll(".png",".pdf");
  corr_canvas_phi[7]->SaveAs(save_name_phi);

  if(!is_number) return 0;

  TString integral_eta_pT_name = "integral_eta_pT";
  
  cintegral_eta_pT[7] = new TCanvas(integral_eta_pT_name,"",10,10,1500,500);
  cintegral_eta_pT[7]->Divide(4,1,0.,0.);


  TString integral_eta_cent_name = "integral_eta_cent";
 
  cintegral_eta_cent[7] = new TCanvas(integral_eta_cent_name,"",10,10,3000,500);
  cintegral_eta_cent[7]->Divide(8,1,0.,0.);

  
  
  for(int j = 0; j<4; j++){

    
    //    in_name = make_name("Result_",g,3,j,0,centlabel,pTlabel);


    cintegral_eta_pT[7]->cd(4-j);

    TString ClosureIntegralEtaPt_name = in_name;
     
    ClosureIntegralEtaPt_name +=j;

    cout<<pTbin_centers.size()<<" "<<Closure_integral_eta0.size()<<endl;


    switch(j){
    case 0:
      Closure_integral_eta_pT[7][j] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta0[0],&pTbin_errors[0],&closure_integral_errors0[0]);
      break;
    case 1:
      Closure_integral_eta_pT[7][j] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta1[0],&pTbin_errors[0],&closure_integral_errors1[0]);
      break;
    case 2:
      Closure_integral_eta_pT[7][j] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta2[0],&pTbin_errors[0],&closure_integral_errors2[0]);
      break;
    case 3:
      Closure_integral_eta_pT[7][j] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta3[0],&pTbin_errors[0],&closure_integral_errors3[0]);
      break;

    }

    Closure_integral_eta_pT[7][j]->SetName(ClosureIntegralEtaPt_name);
        

    Closure_integral_eta_pT[7][j]->SetMarkerColor(1);
    Closure_integral_eta_pT[7][j]->SetMarkerSize(1);
    Closure_integral_eta_pT[7][j]->SetLineColor(1);
    Closure_integral_eta_pT[7][j]->SetMarkerStyle(10);
    Closure_integral_eta_pT[7][j]->SetFillColor(kGreen);


    Closure_integral_eta_pT[7][j]->SetMinimum(-2.);
    Closure_integral_eta_pT[7][j]->SetMaximum(5.);
	    
    Closure_integral_eta_pT[7][j]->GetXaxis()->SetRangeUser(.5,20.);
    Closure_integral_eta_pT[7][j]->GetYaxis()->SetNdivisions(306);
    Closure_integral_eta_pT[7][j]->Draw("p e2 A");

    Closure_integral_eta_pT[7][j]->Draw("same p X");
	 

    Closure_integral_eta_pT[7][j]->GetYaxis()->SetLabelSize(0.07);
	   


    Closure_integral_eta_pT[7][j]->GetXaxis()->SetTitle("Track p_{T} (GeV/c)");
    Closure_integral_eta_pT[7][j]->GetXaxis()->SetTitleSize(0.06);
    Closure_integral_eta_pT[7][j]->GetXaxis()->SetTitleOffset(xoffset+0.2);
    Closure_integral_eta_pT[7][j]->GetYaxis()->SetTitle("(dN/dp_{T})_{P+H} - (dN/dp_{T})_{PYTH} (GeV/c)^{-1}");

    Closure_integral_eta_pT[7][j]->GetXaxis()->SetNdivisions(8);
   
	
    Closure_integral_eta_pT[7][j]->GetXaxis()->CenterTitle();
    Closure_integral_eta_pT[7][j]->GetYaxis()->CenterTitle();
	   
    if(j<3){
      Closure_integral_eta_pT[7][j]->GetYaxis()->SetTitleSize(0.0);
      Closure_integral_eta_pT[7][j]->GetYaxis()->SetLabelSize(0.0);
      Closure_integral_eta_pT[7][j]->GetXaxis()->SetTitleSize(0.07);
      Closure_integral_eta_pT[7][j]->GetXaxis()->SetLabelSize(0.07);
      Closure_integral_eta_pT[7][j]->GetXaxis()->SetTitleOffset(xoffset+0.15);
    }else{
      Closure_integral_eta_pT[7][j]->GetXaxis()->SetLabelSize(ts3);
      Closure_integral_eta_pT[7][j]->GetXaxis()->SetLabelOffset(0.015);
      Closure_integral_eta_pT[7][j]->GetYaxis()->SetTitleOffset(1.);
      Closure_integral_eta_pT[7][j]->GetYaxis()->SetTitleSize(0.06);
      Closure_integral_eta_pT[7][j]->GetYaxis()->SetLabelSize(0.06);
    }



    Closure_integral_eta_pT[7][j]->SetMarkerSize(1.5);



    linePt = new TLine(.5,0,20.,0);
    linePt->SetLineStyle(2);
    linePt->SetLineWidth(1);
    linePt->Draw("same");
	
    closure_integral_values.clear();
    closure_integral_errors_all.clear();
  	
    for(int k = 0; k<4; k++){
      double pt_val, x_val;
	  
      Closure_integral_eta_pT[7][j]->GetPoint(k,x_val,pt_val);
      Closure_integral_eta_pT[7][j]->GetPoint(k,x_val,pt_val);
      closure_integral_values.push_back(pt_val);
      closure_integral_errors_all.push_back(pt_val/2.);

	  
    }


      
	 	  
    Closure_integral_eta_pT2[7][j] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&closure_integral_values[0],&pTbin_errors[0],&closure_integral_errors_all[0]);

  
    if(j==3){ 
      l40 = new TLegend(0.2,0.7,0.8,0.8);
      l40->SetName("l40");
      l40->SetTextFont(43);
      l40->SetTextSizePixels(tspixels);
      l40->SetFillColor(kWhite);
      l40->SetLineColor(kWhite);


      l40->AddEntry(Closure_integral_eta_pT[7][j],"Inclusive P+H","p");


      l40->Draw("same");

    }
      
    //     drawlabels_int_pt2(g,j);
      
    if(j==3){
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
    
  for(int i = 0; i<nTrkPtBins; i++){
	
    TString ClosureIntegralEtaCent_name = "ClosureIntegralEtaCent";
    ClosureIntegralEtaCent_name+=i;
    
       Closure_integral_eta_cent[7][i] = new TH1D(ClosureIntegralEtaCent_name,"",4,xAxis);

    for(int k=0; k<4; k++){
      evalpt = pTbin_centers.at(i);
      value = Closure_integral_eta_pT[7][k]->Eval(evalpt);
      Closure_integral_eta_cent[7][i]->SetBinContent(4-k,value);
    }
  
    Closure_integral_eta_cent[7][i]->SetMarkerStyle(10);
    

    Closure_integral_eta_cent[7][i]->SetMarkerSize(1.5);
    Closure_integral_eta_cent[7][i]->SetLineColor(kBlack);

    Closure_integral_eta_cent[7][i]->SetLineColor(kBlack);
    Closure_integral_eta_cent[7][i]->SetLineColor(kBlack);
    Closure_integral_eta_cent[7][i]->GetYaxis()->SetNdivisions(306);

    cintegral_eta_cent[7]->cd(i+1);


    TString histnameblank = "blank_hist";
    histnameblank+=i;

    blank[i] = new TH1D(histnameblank,"",4,xAxis);

	
    TString histnameblank2 = "blank_hist2";
      histnameblank2+=i;

	
    blank[i]->SetMinimum(-2.);
    blank[i]->SetMaximum(5.);
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


    Closure_integral_eta_cent[7][i]->SetMarkerColor(kBlack);
    Closure_integral_eta_cent[7][i]->Draw("same p X");

    Closure_integral_eta_cent[7][i]->Write();
  

    TLatex *pt_label = new TLatex(0.2,0.95, TrkPtBin_labels[i]);
    pt_label->SetNDC();
    pt_label->SetTextSizePixels(tspixels);
    pt_label->Draw();

    if(i==0)	l40->Draw("same");


  
    cout<<TrkPtBin_strs[i]<<" "<<Closure_integral_eta0.at(i)<<" "<<closure_integral_errors0.at(i)<<" "<<Closure_integral_eta1.at(i)<<" "<<closure_integral_errors1.at(i)<<" "<<Closure_integral_eta2.at(i)<<" "<<closure_integral_errors2.at(i)<<" "<<Closure_integral_eta3.at(i)<<" "<<closure_integral_errors3.at(i)<<" "<<endl;
  }

			       
  cintegral_eta_pT[7]->cd(0);
								      
  TLatex *canvas_title = new TLatex(0.06,0.9,"CMS Preliminary Simulation");
  canvas_title->SetTextSizePixels(tspixels);
  canvas_title->SetTextFont(63);
  canvas_title->Draw();

  TLatex *canvas_title2 = new TLatex(0.05,0.9,"Pythia+Hydjet SpillOver");
  canvas_title2->SetTextSizePixels(tspixels);
  canvas_title2->Draw();

  cintegral_eta_pT[7]->SaveAs("Integral_SpillOver_pT.pdf");
  cintegral_eta_pT[7]->SaveAs("Integral_SpillOver_pT.png");


  cintegral_eta_cent[7]->cd(0);

  canvas_title->Draw();

  canvas_title2->Draw();

   

  cintegral_eta_cent[7]->SaveAs("Integral_SpillOver_Cent.png");
  cintegral_eta_cent[7] ->SaveAs("Integral_SpillOver_Cent.pdf");
    

  return 0;
}
