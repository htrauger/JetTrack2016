
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

#include "../JetTrack2016_functions.h"

#include <iostream>
#include <vector>
#include <fstream>


using namespace std;

Int_t width_determination(bool with_ref = 0){

  
#include "../HIN_14_016_universals.h"

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.15);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);
  gStyle->SetTextFont(43);


  
  TFile *fin[12];
  TFile *fin_ref[12];
  Double_t xAxis[5] = {-100,-50,-30,-10,0}; 
  
  TH1D* int_cent[12][4];
  TH1D* blank[4];
  TH1D* blank2[4];
  TH1D* blank3[4];
  
  

  TH1D *check_new_phi_rebin[12][5][4];
  TH1D *check_new_phi_syst[12][5][4];
  
  TH1D *check_new_eta_rebin[12][5][4];
  TH1D *check_new_eta_syst[12][5][4];

  TH1D *PbPb_pp_eta[12][5][4];
  TH1D *PbPb_pp_eta_syst[12][5][4];

  TH1D *check_new_eta_up[12][5][4];
  TH1D *check_new_eta_down[12][5][4];


  TH1D *PbPb_pp_eta_up[12][5][4];
  TH1D *PbPb_pp_eta_down[12][5][4];



  TH1D *check_new_phi_up[12][5][4];
  TH1D *check_new_phi_down[12][5][4];

  TH1D *PbPb_pp_phi_up[12][5][4];
  TH1D *PbPb_pp_phi_down[12][5][4];



  float ref_bc, ref_err, new_bc, new_err, Width_val_eta_temp, Width_err_eta_temp, Width_val_phi_temp, Width_err_phi_temp, Width_diff_val_eta_temp, Width_diff_err_eta_temp, Width_diff_val_phi_temp, Width_diff_err_phi_temp, val_l,val_r,err_l,err_r, spill;

  TGraphAsymmErrors *Width_eta[12][4];
  TGraphAsymmErrors *Width_diff_eta[12][4];
  TGraphAsymmErrors *Width_phi[12][4];
  TGraphAsymmErrors *Width_diff_phi[12][4];
 TGraphAsymmErrors *Width_diff_phi_ref[12][4];
 TGraphAsymmErrors *Width_diff_eta_ref[12][4];
  
  TGraphAsymmErrors *Width_eta_cent[12][4];
  TGraphAsymmErrors *Width_diff_eta_cent[12][4];
  
  TGraphAsymmErrors *Width_eta_ref[12][4];
  TGraphAsymmErrors *Width_phi_ref[12][4];
  TGraphAsymmErrors *Width_eta_cent_ref[12][4];


  TGraphAsymmErrors *Width_PbPb_minus_pp_eta[12][4];
  TGraphAsymmErrors *Width_PbPb_minus_pp_phi[12][4];

  TCanvas *Width_eta_canvas[12];
  TCanvas *Width_phi_canvas[12];
  TCanvas *cWidth_eta_cent[12];
  TCanvas *check_fits_diff[12]; 


  vector<float> Width_val_eta;
  vector<float> Width_err_eta_up;
  vector<float> Width_err_eta_down;


  vector<float> Width_diff_val_eta;
  vector<float> Width_diff_err_eta_up;
  vector<float> Width_diff_err_eta_down;

 
  vector<float> Width_val_phi;
  vector<float> Width_err_phi_up;
  vector<float> Width_err_phi_down;


  vector<float> Width_err_eta_up_combined;
  vector<float> Width_err_eta_down_combined;

  vector<float> Width_err_phi_up_combined;
  vector<float> Width_err_phi_down_combined;


  vector<float> Width_diff_val_phi;
  vector<float> Width_diff_err_phi_up;
  vector<float> Width_diff_err_phi_down;

   
  vector<float> Width_val_cent_eta;
  vector<float> Width_err_cent_eta_up;
  vector<float> Width_err_cent_eta_down;

   
  vector<float> Width_val_cent_phi;
  vector<float> Width_err_cent_phi_up;
  vector<float> Width_err_cent_phi_down;

    
  vector<float> Width_diff_val_cent_eta;
  vector<float> Width_diff_err_cent_eta_up;
  vector<float> Width_diff_err_cent_eta_down;

   
  vector<float> Width_diff_val_cent_phi;
  vector<float> Width_diff_err_cent_phi_up;
  vector<float> Width_diff_err_cent_phi_down;

  float evalpt, value, exh,exl,eyh,eyl;
  


  TLegend *l40,*l41,*l42;

  TString datalabel, in_name, centlabel, pTlabel;
  TString   leadingrefetaname[5][4];
  float raw_min,raw_max,mixed_min,mixed_max,yield_min,yield_max,result_min,result_max, temp_int, temp_sigma, full_int;
  int counter, lbin, rbin;

  TLine *linePhi,*lineEta, *lineCent;

  TLatex *tex24phi,*tex24eta;



  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 9;


  vector<float> pTbin_centers;
  pTbin_centers.push_back(1.5);
  pTbin_centers.push_back(2.5);
  pTbin_centers.push_back(3.5);
  pTbin_centers.push_back(6.0);
  pTbin_centers.push_back(10.0);
  pTbin_centers.push_back(14.0);
  pTbin_centers.push_back(18.0);
  pTbin_centers.push_back(22.0);
  
  vector<float> pTbin_errors;
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(2.);
  pTbin_errors.push_back(2.);
  pTbin_errors.push_back(2.);
  pTbin_errors.push_back(2.);
  pTbin_errors.push_back(2.);
 

 
  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

   float TrkPtBins[nTrkPtBins+1] = {05, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.5<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
  


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


  vector<float> cent_centers;
  cent_centers.push_back(-75);
  cent_centers.push_back(-40);
  cent_centers.push_back(-20);
  cent_centers.push_back(-5);
  vector<float> cent_errors;
  cent_errors.push_back(25);
  cent_errors.push_back(10);
  cent_errors.push_back(10);
  cent_errors.push_back(5);

  TFile *fout[6];
  
  TCanvas *check_fits_eta[6];
  TCanvas *check_fits_phi[6];


 
  
  float par0, par1, par2, par3, par4, special_err, error_up, error_down,   val_2, err_syst, dx_eta;

  //-------------------------
  //   PARAMETER SET

  float range = 1.;
 
  //----------------------------



  TF1 *double_gaus = new TF1("double_gaus","[0]*TMath::Exp(-pow(TMath::Abs(x)/[1],2))+[2]*TMath::Exp(-pow(TMath::Abs(x)/[3],2))+[4]",-range,range);

  TF1 *gaus_1_eta[12][5][4];
  TF1 *gaus_2_eta[12][5][4];
  TF1 *gaus_tot_eta[12][5][4];

 TF1 *gaus_1_phi[12][5][4];
  TF1 *gaus_2_phi[12][5][4];
  TF1 *gaus_tot_phi[12][5][4];

  TF1 *flat_line = new TF1("flat_line","-[0]+x-x",-3.5,3.5);


  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------

  fin[0] =  new TFile("../particle_yields/Particle_Yields.root", "READ");

  for(int g = 0; g<2; g++){
      
    cout<<"starting getting files for g = "<<g<<endl;
  
    if(g==1){
      fout[g] = new TFile((TString)("Width_Output_"+datalabel+".root"),"RECREATE");
    }

    cout<<"got files"<<endl;

    TString Width_eta_canvas_name = "Width_eta";  Width_eta_canvas_name+=g;
    Width_eta_canvas[g] = new TCanvas(Width_eta_canvas_name,"",10,10,1500,850);
    Width_eta_canvas[g]->Divide(4,2,0.,0.);


    TString Width_phi_canvas_name = "Width_phi";  Width_phi_canvas_name+=g;
    Width_phi_canvas[g] = new TCanvas(Width_phi_canvas_name,"",10,10,1500,850);
    Width_phi_canvas[g]->Divide(4,2,0.,0.);


    TString cWidth_eta_cent_name = "Width_eta_cent";  cWidth_eta_cent_name+=g;
    if(g==5){  cWidth_eta_cent[g] = new TCanvas(cWidth_eta_cent_name,"",10,10,1500,800);
      cWidth_eta_cent[g]->Divide(2,2,0.,0.);
    }else{
      cWidth_eta_cent[g] = new TCanvas(cWidth_eta_cent_name,"",10,10,1500,400);
      cWidth_eta_cent[g]->Divide(2,1,0.,0.);
    }

    TString checkfitname_eta = "CheckFitEta";
    checkfitname_eta+=g;
    check_fits_eta[g] = new TCanvas(checkfitname_eta," ",10,10,1500,1600);
    check_fits_eta[g]->Divide(4,4,0.,0.);
    
    TString checkfitname_phi = "CheckFitPhi";
    checkfitname_phi+=g;
    check_fits_phi[g] = new TCanvas(checkfitname_phi," ",10,10,1500,1600);
    check_fits_phi[g]->Divide(4,4,0.,0.);
    


    TString checkfitdiffname = "CheckFitDiff";
    checkfitdiffname+=g;
    check_fits_diff[g] = new TCanvas(checkfitdiffname," ",10,10,1500,800);
    check_fits_diff[g]->Divide(4,2,0.,0.);

    cout<<"did canvases"<<endl;

    

    TCanvas *dummy = new TCanvas("dummy");
  
    for(int j = 0; j<4; j++){

      Width_val_eta.clear();
      Width_err_eta_up.clear();
      Width_err_eta_down.clear();
 

      Width_diff_val_eta.clear();
      Width_diff_err_eta_up.clear();
      Width_diff_err_eta_down.clear();
 
      Width_val_phi.clear();
      Width_err_phi_up.clear();
      Width_err_phi_down.clear();
 

      Width_diff_val_phi.clear();
      Width_diff_err_phi_up.clear();
      Width_diff_err_phi_down.clear();
 



      for(int i=0; i<4; i++){

	dummy->cd();

	//-------------------------------------
	//  Get and assign all histograms 
	//-------------------------------------

	if(g==0) check_new_eta_rbin[g][i][j] = (TH1D*)result[g][i][j]->ProjectionX((TString)("Proj_dEta_PbPb_"+  CBin_strs[j] + "_" + CBin_strs[j+1]+" "+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]+"_Rebin"));
	else 	check_new_eta_rbin[g][i][j] = (TH1D*)result[g][i][j]->ProjectionX((TString)("Proj_dEta_pp_"+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]+"_Rebin"));

	if(g==0) check_new_phi_rbin[g][i][j] = (TH1D*)result[g][i][j]->ProjectionX((TString)("Proj_dPhi_PbPb_"+  CBin_strs[j] + "_" + CBin_strs[j+1]+" "+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]+"_Rebin"));
	else 	check_new_phi_rbin[g][i][j] = (TH1D*)result[g][i][j]->ProjectionX((TString)("Proj_dPhi_pp_"+ TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]+"_Rebin"));


       	cout<<"got histos"<<endl;
	/*	
	if(g==0){
 

	  TString newchecknameEtasyst = newchecknameEta_rebin;
	  newchecknameEtasyst.ReplaceAll("New_check","Syst_error");
	  newchecknameEtasyst.ReplaceAll("rebin","");
	  check_new_eta_syst[g][i][j] = (TH1D*)fin[g]->Get(newchecknameEtasyst)->Clone(newchecknameEtasyst);
	

	  TString newchecknamePhisyst = newchecknamePhi_rebin;
	  newchecknamePhisyst.ReplaceAll("New_check","Syst_error");
	  newchecknamePhisyst.ReplaceAll("rebin","");
	  check_new_phi_syst[g][i][j] = (TH1D*)fin[g]->Get(newchecknamePhisyst)->Clone(newchecknamePhisyst);




	  newchecknameEta_rebin.ReplaceAll("rebin","up");
	  check_new_eta_up[g][i][j] = (TH1D*)check_new_eta_rebin[g][i][j]->Clone(newchecknameEta_rebin);
	
	  newchecknameEta_rebin.ReplaceAll("up","down");
	  check_new_eta_down[g][i][j] = (TH1D*)check_new_eta_rebin[g][i][j]->Clone(newchecknameEta_rebin);
	  
	  nbins = check_new_eta_rebin[g][i][j]->GetNbinsX()+1;
	  
	  for(int k = 1; k<nbins/2+1;k++){
	    val_l = check_new_eta_rebin[g][i][j]->GetBinContent(k);
	 
	    err_l  = check_new_eta_rebin[g][i][j]->GetBinError(k);
	    
	    val_r = check_new_eta_rebin[g][i][j]->GetBinContent(nbins-k);
	    err_r = check_new_eta_rebin[g][i][j]->GetBinError(nbins-k);
	    
	    //cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<endl;

	    check_new_eta_rebin[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_eta_rebin[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	    check_new_eta_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_eta_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	   
	    if(k<nbins/2){
	      check_new_eta_rebin[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      check_new_eta_rebin[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    }


	    if(g==0||g==2||g==4){
	      val_2 = check_new_eta_rebin[g][i][j]->GetBinContent(k);
	      err_syst  = check_new_eta_syst[g][i][j]->GetBinError(k);
	      
	      check_new_eta_up[g][i][j]->SetBinContent(k,val_2+err_syst);
	      check_new_eta_up[g][i][j]->SetBinContent(nbins-k,val_2+err_syst);
	      check_new_eta_down[g][i][j]->SetBinContent(k,val_2-err_syst);
	      check_new_eta_down[g][i][j]->SetBinContent(nbins-k,val_2-err_syst);

	      if(k<nbins/2){
		check_new_eta_up[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_eta_up[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_eta_down[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_eta_down[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      }
	    }
	
	  }
	
	*/
      

	  newchecknamePhi_rebin.ReplaceAll("rebin","up");
	  check_new_phi_up[g][i][j] = (TH1D*)check_new_phi_rebin[g][i][j]->Clone(newchecknamePhi_rebin);
	
	  newchecknamePhi_rebin.ReplaceAll("up","down");
	  check_new_phi_down[g][i][j] = (TH1D*)check_new_phi_rebin[g][i][j]->Clone(newchecknamePhi_rebin);


	  
	  nbins = check_new_phi_rebin[g][i][j]->GetNbinsX()+1;
	  
	  for(int k = 1; k<nbins/2+1;k++){
	    val_l = check_new_phi_rebin[g][i][j]->GetBinContent(k);
	 
	    err_l  = check_new_phi_rebin[g][i][j]->GetBinError(k);
	    
	    val_r = check_new_phi_rebin[g][i][j]->GetBinContent(nbins-k);
	    err_r = check_new_phi_rebin[g][i][j]->GetBinError(nbins-k);
	    
	    // cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<endl;

	    check_new_phi_rebin[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_phi_rebin[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	    check_new_phi_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_phi_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	   
	    if(k<nbins/2){
	      check_new_phi_rebin[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      check_new_phi_rebin[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    }


	    if(g==0||g==2||g==4){
	      val_2 = check_new_phi_rebin[g][i][j]->GetBinContent(k);
	      err_syst  = check_new_phi_syst[g][i][j]->GetBinError(k);
	    
	      check_new_phi_up[g][i][j]->SetBinContent(k,val_2+err_syst);
	      check_new_phi_up[g][i][j]->SetBinContent(nbins-k,val_2+err_syst);
	      check_new_phi_down[g][i][j]->SetBinContent(k,val_2-err_syst);
	      check_new_phi_down[g][i][j]->SetBinContent(nbins-k,val_2-err_syst);

	      if(k<nbins/2){
		check_new_phi_up[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_phi_up[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_phi_down[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_phi_down[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      }
	    }
	
	  }
	}
	
	//Difference histos

	if(g==1){
	 
	  PbPb_pp_eta[g][i][j] = (TH1D*)check_new_eta_rebin[g-1][i][j]->Clone(Difference_name);

	  PbPb_pp_eta_up[g][i][j] = (TH1D*)check_new_eta_up[g-1][i][j]->Clone(Difference_name);

	  PbPb_pp_eta_down[g][i][j] = (TH1D*)check_new_eta_down[g-1][i][j]->Clone(Difference_name);


	  PbPb_pp_eta[g][i][j]->Add(check_new_eta_rebin[g][i][j],-1.);

	  PbPb_pp_eta_up[g][i][j]->Add(check_new_eta_rebin[g][i][j],-1.);

	  PbPb_pp_eta_down[g][i][j]->Add(check_new_eta_rebin[g][i][j],-1.);
	
	}

	//--------------------------------------
	//    Get Width and error
	//--------------------------------------

	cout<<"*********************************************"<<endl;

	cout<<in_name<<endl;

	cout<<g<<" "<<i<<" "<<j<<" "<<endl;
	
	double_gaus->ReleaseParameter(0);
	double_gaus->ReleaseParameter(2);
	double_gaus->ReleaseParameter(1);
	double_gaus->ReleaseParameter(3);
	//	double_gaus->ReleaseParameter(4);

	double_gaus->FixParameter(4,0.);

	double_gaus->SetParLimits(0,0.0,10.);  
	double_gaus->SetParLimits(2,0.0,10.);  
	

	double_gaus->SetParameter(1,0.3);
	//	double_gaus->SetParLimits(1,0.1,0.67*range);
	double_gaus->SetParLimits(1,0.1,range);
	double_gaus->SetParameter(3,0.6);
	//	double_gaus->SetParLimits(3,0.1,0.67*range);
	double_gaus->SetParLimits(3,0.1,range);
	double_gaus->SetParameter(0,1.);
	double_gaus->SetParameter(2,1.);

	check_new_eta_rebin[g][i][j]->Fit(double_gaus,"","",-range,range);

	
	gaus_1_eta[g][i][j] = new TF1("gaus_1_eta","[0]*TMath::Exp(-pow(TMath::Abs(x)/[1],2))",-range,range);
	gaus_2_eta[g][i][j] = new TF1("gaus_2_eta","[0]*TMath::Exp(-pow(TMath::Abs(x)/[1],2))",-range,range);
	gaus_tot_eta[g][i][j] =  new TF1("gaus_tot_eta","[0]*TMath::Exp(-pow(TMath::Abs(x)/[1],2))+[2]*TMath::Exp(-pow(TMath::Abs(x)/[3],2))",-range,range);


	par0 = double_gaus->GetParameter(0);
	par1 = double_gaus->GetParameter(1);
	par2 = double_gaus->GetParameter(2);
	par3 = double_gaus->GetParameter(3);


	//	temp_sigma = TMath::Sqrt((par0*par1*par1+par2*par3*par3)/(par0+par2));

	gaus_1_eta[g][i][j]->SetParameter(0,par0);
	gaus_1_eta[g][i][j]->SetParameter(1,par1);
	gaus_2_eta[g][i][j]->SetParameter(0,par2);
	gaus_2_eta[g][i][j]->SetParameter(1,par3);

	gaus_tot_eta[g][i][j]->SetParameter(0,par0);
	gaus_tot_eta[g][i][j]->SetParameter(1,par1);
	gaus_tot_eta[g][i][j]->SetParameter(2,par2);
	gaus_tot_eta[g][i][j]->SetParameter(3,par3);


	full_int = gaus_tot_eta[g][i][j]->Integral(0,range);
	temp_int = 0;
	temp_sigma = 0;

	counter = 0;

	while(temp_int <0.67*full_int){
	  temp_sigma = range*counter/50000;
	  temp_int = gaus_tot_eta[g][i][j]->Integral(0, temp_sigma);
	  counter++;
	}

	Width_val_eta.push_back(temp_sigma);

	if(g==0||g==2||g==4){
	  check_new_eta_up[g][i][j]->Fit(double_gaus,"","",-range,range);

	  //	  double_gaus->SetParameter(4,0.);

	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }
	
	  error_up = temp_sigma-Width_val_eta.at(i);


	  check_new_eta_down[g][i][j]->Fit(double_gaus,"","",-range,range);

	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }
	  


		  
	  error_down = temp_sigma-Width_val_eta.at(i);

	  //	  error_down = error_up;
	
	  if(error_up>0&&error_down<0){
	    Width_err_eta_up.push_back(TMath::Abs(error_up));
	    Width_err_eta_down.push_back(TMath::Abs(error_down));
	  }else if(error_up<0&&error_down>0){
	    Width_err_eta_down.push_back(TMath::Abs(error_up));
	    Width_err_eta_up.push_back(TMath::Abs(error_down));
	  } else {
	    Width_err_eta_up.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    Width_err_eta_down.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    cerr<<"Both errors go same way "<<error_up<<" "<<error_down<<endl;
	    //return -1;
	  }

	  
	} //close error for PbPb only
      

	///------------------
	//   Phi

	//-----------


	cout<<in_name<<endl;

	cout<<g<<" "<<i<<" "<<j<<" "<<endl;

	check_new_phi_rebin[g][i][j]->Fit(double_gaus,"","",-range,range);

	
	gaus_1_phi[g][i][j] = new TF1("gaus_1_phi","[0]*TMath::Exp(-pow(TMath::Abs(x)/[1],2))",-range,range);
	gaus_2_phi[g][i][j] = new TF1("gaus_2_phi","[0]*TMath::Exp(-pow(TMath::Abs(x)/[1],2))",-range,range);
	gaus_tot_phi[g][i][j] =  new TF1("gaus_tot_phi","[0]*TMath::Exp(-pow(TMath::Abs(x)/[1],2))+[2]*TMath::Exp(-pow(TMath::Abs(x)/[3],2))",-range,range);


	par0 = double_gaus->GetParameter(0);
	par1 = double_gaus->GetParameter(1);
	par2 = double_gaus->GetParameter(2);
	par3 = double_gaus->GetParameter(3);


	//	temp_sigma = TMath::Sqrt((par0*par1*par1+par2*par3*par3)/(par0+par2));

	gaus_1_phi[g][i][j]->SetParameter(0,par0);
	gaus_1_phi[g][i][j]->SetParameter(1,par1);
	gaus_2_phi[g][i][j]->SetParameter(0,par2);
	gaus_2_phi[g][i][j]->SetParameter(1,par3);

	gaus_tot_phi[g][i][j]->SetParameter(0,par0);
	gaus_tot_phi[g][i][j]->SetParameter(1,par1);
	gaus_tot_phi[g][i][j]->SetParameter(2,par2);
	gaus_tot_phi[g][i][j]->SetParameter(3,par3);


	full_int = gaus_tot_phi[g][i][j]->Integral(0,range);
	temp_int = 0;
	temp_sigma = 0;

	counter = 0;
	while(temp_int <0.67*full_int){
	  temp_sigma = range*counter/50000;
	  temp_int = gaus_tot_phi[g][i][j]->Integral(0, temp_sigma);
	  counter++;
	}

	Width_val_phi.push_back(temp_sigma);

	if(g==0||g==2||g==4){
	  check_new_phi_up[g][i][j]->Fit(double_gaus,"","",-range,range);

	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }
	
	  error_up = temp_sigma-Width_val_phi.at(i);

	  check_new_phi_down[g][i][j]->Fit(double_gaus,"","",-range,range);

	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }
	  
		  
	  error_down = temp_sigma-Width_val_phi.at(i);

	  //	  error_down = error_up;
	
	  if(error_up>0&&error_down<0){
	    Width_err_phi_up.push_back(TMath::Abs(error_up));
	    Width_err_phi_down.push_back(TMath::Abs(error_down));
	  }else if(error_up<0&&error_down>0){
	    Width_err_phi_down.push_back(TMath::Abs(error_up));
	    Width_err_phi_up.push_back(TMath::Abs(error_down));
	  } else {
	    Width_err_phi_up.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    Width_err_phi_down.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    cerr<<"Both errors go same way "<<error_up<<" "<<error_down<<endl;
	    //return -1;
	  }

	  
	} //close error for PbPb only
      
      }// close i


      cout<<"ready to set width plots"<<endl;

    
      Width_eta[g][j] = new TGraphAsymmErrors(pTbin_centers.size(),&pTbin_centers[0],&Width_val_eta[0],&pTbin_errors[0],&pTbin_errors[0],&Width_err_eta_down[0],&Width_err_eta_up[0]);
      TString Width_eta_name = "RMS_eta"; Width_eta_name+= g; Width_eta_name+= j;
      Width_eta[g][j]->SetName(Width_eta_name);
   
     
   
      Width_phi[g][j] = new TGraphAsymmErrors(pTbin_centers.size(),&pTbin_centers[0],&Width_val_phi[0],&pTbin_errors[0],&pTbin_errors[0],&Width_err_phi_down[0],&Width_err_phi_up[0]);


      TString Width_phi_name = "RMS_phi"; Width_phi_name+= g; Width_phi_name+= j;


      Width_phi[g][j]->SetName(Width_phi_name);

    
      if(with_ref){
	Width_eta_ref[g][j] = (TGraphAsymmErrors*)fin_ref[g]->Get(Width_eta_name)->Clone((TString)(Width_eta_name+="_ref"));
	Width_phi_ref[g][j] = (TGraphAsymmErrors*)fin_ref[g]->Get(Width_phi_name)->Clone((TString)(Width_phi_name+="_ref"));


	Width_eta_ref[g][j]->SetMarkerColor(kRed);      
	Width_eta_ref[g][j]->SetLineColor(kRed);      

	Width_phi_ref[g][j]->SetMarkerColor(kRed);      
	Width_phi_ref[g][j]->SetLineColor(kRed);      
      }
      if(g%2!=0){
	
	//------------------------
	//   Width of difference
	//-----------------------

	cout<<"going to try"<<endl;
	if(with_ref){
	  Width_diff_eta_ref[g][j] = (TGraphAsymmErrors*)Width_eta_ref[g-1][j]->Clone(Form("WidthDiffEtaRef%d%d",g,j));

	  Width_diff_phi_ref[g][j] = (TGraphAsymmErrors*)Width_phi_ref[g-1][j]->Clone(Form("WidthDiffPhiRef%d%d",g,j));
	  cout<<"here"<<endl;

	  double x,y1,y2,y;
	  for(int k = 0; k<4; k++){
	    Width_eta_ref[g-1][j]->GetPoint(k,x,y1);
	    Width_eta_ref[g][j]->GetPoint(k,x,y2);
	    y = y1-y2;

	    Width_diff_eta_ref[g][j]->SetPoint(k,x,y);
      

	    Width_phi_ref[g-1][j]->GetPoint(k,x,y1);
	    Width_phi_ref[g][j]->GetPoint(k,x,y2);
	    y = y1-y2;

	    Width_diff_phi_ref[g][j]->SetPoint(k,x,y);
      

	  }
	}
	cout<<"got it"<<endl;

   	for(int i = 0; i<2; i++){
	  ///--------------
	  //   Difference
	  //---------------

	  //	  double_gaus->ReleaseParameter(2);  //we're going to let there be a negative contribution, since it's a difference

	  PbPb_pp_eta[g][i][j]->Fit(double_gaus,"","",-range,range);

	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }

	
	  
	  Width_diff_val_eta.push_back(temp_sigma);
	

	  PbPb_pp_eta_up[g][i][j]->Fit(double_gaus,"","",-range-1.,range+1.);
	
	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }
	   

	
	  error_up = temp_sigma-Width_diff_val_eta.at(i);

	 
	  
	  PbPb_pp_eta_down[g][i][j]->Fit(double_gaus,"","",-range-1.,range+1.);

	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }
	 
	  error_down = temp_sigma-Width_val_eta.at(i);
 
	  if(error_up>0.01&&error_down<-0.01){
	    Width_diff_err_eta_up.push_back(TMath::Abs(error_up));
	    Width_diff_err_eta_down.push_back(TMath::Abs(error_down));
	  }else if(error_up<-0.01&&error_down>0.01){
	    Width_diff_err_eta_down.push_back(TMath::Abs(error_up));
	    Width_diff_err_eta_up.push_back(TMath::Abs(error_down));
	  } else {
	    cerr<<"Both errors go same way"<<endl;
	    Width_diff_err_eta_up.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    Width_diff_err_eta_down.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    //return -1;
	  }
	 
	}
      
      
	Width_diff_eta[g][j] = new TGraphAsymmErrors(pTbin_centers.size(),&pTbin_centers[0],&Width_diff_val_eta[0],&pTbin_errors[0],&pTbin_errors[0],&Width_diff_err_eta_down[0],&Width_diff_err_eta_up[0]);
      

	TString Width_diff_eta_name = "Width_diff_eta"; Width_diff_eta_name+= g; Width_diff_eta_name+= j;

	Width_diff_eta[g][j]->SetName(Width_diff_eta_name);

      }
  
    } //j

    cout<<"ready to start drawing"<<endl;

    for(int k = 0; k<4; k++){
      float temp = (Width_eta[g][0]->GetErrorYhigh(k)+Width_eta[g][1]->GetErrorYhigh(k)+Width_eta[g][2]->GetErrorYhigh(k)+Width_eta[g][3]->GetErrorYhigh(k))/4.;
      for(int j = 0; j<4; j++){
	Width_eta[g][j]->SetPointEYhigh(k,temp);
      }
      temp = (Width_phi[g][0]->GetErrorYhigh(k)+Width_phi[g][1]->GetErrorYhigh(k)+Width_phi[g][2]->GetErrorYhigh(k)+Width_phi[g][3]->GetErrorYhigh(k))/4.;
      for(int j = 0; j<4; j++){
	Width_phi[g][j]->SetPointEYhigh(k,temp);
      }
 
      temp = (Width_eta[g][0]->GetErrorYlow(k)+Width_eta[g][1]->GetErrorYlow(k)+Width_eta[g][2]->GetErrorYlow(k)+Width_eta[g][3]->GetErrorYlow(k))/4.;
      for(int j = 0; j<4; j++){
	Width_eta[g][j]->SetPointEYlow(k,temp);
      }
      temp = (Width_phi[g][0]->GetErrorYlow(k)+Width_phi[g][1]->GetErrorYlow(k)+Width_phi[g][2]->GetErrorYlow(k)+Width_phi[g][3]->GetErrorYlow(k))/4.;
      for(int j = 0; j<4; j++){
	Width_phi[g][j]->SetPointEYlow(k,temp);
      }

    }
        
    for(int j = 0; j<4; j++){

      if(g==1||g==3||g==5){
	TString Width_PbPb_minus_pp_eta_name = "Width_PbPb_minus_pp_eta"; Width_PbPb_minus_pp_eta_name+= g; Width_PbPb_minus_pp_eta_name+= j;
	Width_PbPb_minus_pp_eta[g][j]  = (TGraphAsymmErrors*) Width_eta[g-1][j]->Clone(Width_PbPb_minus_pp_eta_name);
    
	for(int k = 0; k<4; k++){
	  
	  float temp_x = pTbin_centers.at(k);
	  
	  float PbPb_val = Width_eta[g-1][j]->Eval(temp_x);
	  float pp_val = Width_eta[g][j]->Eval(temp_x);

	  Width_PbPb_minus_pp_eta[g][j]->SetPoint(k,temp_x,PbPb_val-pp_val);

	}
           
	Width_eta_canvas[g]->cd(j+1);
	Width_eta[g-1][j]->Draw("AP");

	Width_eta[g-1][j]->GetXaxis()->SetRangeUser(.8,8.2);


	switch(g){
	case 1: 
	  Width_eta[g-1][j]->SetMarkerStyle(10);
	  Width_eta[g][j]->SetMarkerStyle(4);
	  if(with_ref) 	  Width_eta_ref[g-1][j]->SetMarkerStyle(10);
	  if(with_ref) 	  Width_diff_eta_ref[g][j]->SetMarkerStyle(10);
	  if(with_ref) 	  Width_eta_ref[g][j]->SetMarkerStyle(4);
	  Width_eta[g-1][j]->SetFillColor(90);
	  Width_eta[g-1][j]->SetLineColor(90);

	  break;
	case 3: 
	  Width_eta[g-1][j]->SetMarkerStyle(34);
	  Width_eta[g-1][j]->SetMarkerSize(2);
	  if(with_ref) 	  Width_eta_ref[g-1][j]->SetMarkerStyle(34);
	  if(with_ref) 	  Width_diff_eta_ref[g][j]->SetMarkerStyle(34);
	  if(with_ref) 	  Width_eta_ref[g-1][j]->SetMarkerSize(2);
	  if(with_ref) 	  Width_diff_eta_ref[g][j]->SetMarkerSize(2);
	  Width_eta[g-1][j]->SetFillColor(30);
	  Width_eta[g-1][j]->SetLineColor(30);
	  Width_eta[g][j]->SetMarkerStyle(28);
	  break;  
	case 5: 
	  Width_eta[g-1][j]->SetMarkerStyle(21);
	  if(with_ref) 	  Width_eta_ref[g-1][j]->SetMarkerStyle(21);
	  if(with_ref) 	  Width_diff_eta_ref[g][j]->SetMarkerStyle(21);
	  Width_eta[g-1][j]->SetFillColor(kOrange-2);
	  Width_eta[g-1][j]->SetLineColor(kOrange-2);
	  Width_eta[g][j]->SetMarkerStyle(25);
	  break;
	}

	Width_eta[g-1][j]->SetTitle("");

	Width_eta[g][j]->SetMarkerSize(3);
	Width_eta[g-1][j]->SetMarkerSize(2);

	if(with_ref){
	  Width_eta_ref[g-1][j]->SetTitle("");

	  Width_eta_ref[g][j]->SetMarkerSize(2);
	  Width_eta_ref[g-1][j]->SetMarkerSize(2);
	}
	
	if(j==0){
	  Width_eta[g-1][j]->GetYaxis()->SetTitle("#Delta#eta Width");
	  Width_eta[g-1][j]->GetYaxis()->SetTitleSize(tstitle);
	  Width_eta[g-1][j]->GetYaxis()->SetLabelSize(ts2);
	  Width_eta[g-1][j]->GetYaxis()->CenterTitle();

	}else{
	  Width_eta[g-1][j]->GetYaxis()->SetLabelSize(0);
	}

	Width_eta[g-1][j]->SetMarkerColor(kBlack);
	Width_eta[g-1][j]->SetMinimum(-0.02);
	Width_eta[g-1][j]->SetMaximum(0.75);

	Width_eta[g-1][j]->SetLineWidth(0);


	Width_eta[g-1][j]->GetXaxis()->SetNdivisions(8);
	Width_eta[g-1][j]->GetYaxis()->SetNdivisions(306);

	Width_eta[g-1][j]->Draw("p Z e2");

	Width_eta[g][j]->SetMarkerSize(2);

	Width_eta[g-1][j]->Draw("p Z same");
	Width_eta[g][j]->Draw("p X same");

	if(with_ref)	Width_eta_ref[g-1][j]->Draw("same P");
	if(with_ref)	Width_eta_ref[g][j]->Draw("same P X");
	
	if(j==0){

	  l40 = new TLegend(textalign2-.05,texty3-0.05,0.6,texty1-0.05);
	  if(g==1){  
	    l40->AddEntry(Width_eta[g-1][j],"PbPb Inclusive Jets","lpfe"); 
	    l40->AddEntry(Width_eta[g][j],"pp Inclusive Jets","p"); 
	  }
	  if(g==3){  
	    l40->AddEntry(Width_eta[g-1][j],"PbPb Subleading Jets","lpfe"); 
	    l40->AddEntry(Width_eta[g][j],"pp Subleading Jets","p"); 
	  }
	  if(g==5){  
	    l40->AddEntry(Width_eta[g-1][j],"PbPb Leading Jets","lpfe"); 
	    l40->AddEntry(Width_eta[g][j],"pp Leading Jets","p"); 
	  }
	  l40->SetTextFont(43);
	  l40->SetTextSizePixels(tspixels);
	  l40->SetLineColor(kWhite);
	  l40->Draw();
	}

	TLine *line = new TLine(.8,0.,8.2,0.);
	//	line->SetLineStyle(2);
	//	line->Draw("same");

	drawlabels_int_pt2(g,j);

	  
	if(j==3){
	  tex24eta = new TLatex(textalign,texty2,phirangelabel);
	  tex24eta->SetName("tex24phi");
	  tex24eta->SetNDC();
	  tex24eta->SetTextSizePixels(tspixels);
	  tex24eta->Draw();
	}


	//-------------------
	//  Draw  WidthPbPb - Widthpp
	//----------------
      


	Width_eta_canvas[g]->cd(j+5);
	Width_PbPb_minus_pp_eta[g][j]->Draw("AP");
	Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetRangeUser(0.8,8.2);
	Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetNdivisions(8);
	Width_PbPb_minus_pp_eta[g][j]->GetYaxis()->SetNdivisions(204);

	switch(g){
	case 1: 
	  Width_PbPb_minus_pp_eta[g][j]->SetMarkerStyle(10);
	  Width_PbPb_minus_pp_eta[g][j]->SetFillColor(90);
	  Width_PbPb_minus_pp_eta[g][j]->SetLineColor(90);
	  break;
	case 3: 
	  Width_PbPb_minus_pp_eta[g][j]->SetMarkerStyle(34);
	  Width_PbPb_minus_pp_eta[g][j]->SetMarkerSize(2);
	  Width_PbPb_minus_pp_eta[g][j]->SetFillColor(30);
	  Width_PbPb_minus_pp_eta[g][j]->SetLineColor(30);
	  break;  
	case 5: 
	  Width_PbPb_minus_pp_eta[g][j]->SetMarkerStyle(21);
	  Width_PbPb_minus_pp_eta[g][j]->SetFillColor(kOrange-2);
	  Width_PbPb_minus_pp_eta[g][j]->SetLineColor(kOrange-2);
	  break;
	}


	Width_PbPb_minus_pp_eta[g][j]->SetMarkerSize(2);
	Width_PbPb_minus_pp_eta[g][j]->SetTitle("");
	Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetLabelSize(ts);
	Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetTitleSize(tstitle-.002);
	Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->CenterTitle();
	Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetTitleOffset(0.83);
	Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV)");



	if(j==0){
	  Width_PbPb_minus_pp_eta[g][j]->GetYaxis()->SetTitle("(#Delta#eta Width)_{PbPb} - (#Delta#eta Width)_{pp}");
	  Width_PbPb_minus_pp_eta[g][j]->GetYaxis()->SetTitleSize(ts);
	  Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetTitleSize(ts-.004);
	  Width_PbPb_minus_pp_eta[g][j]->GetYaxis()->SetLabelSize(ts3);
	  Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetLabelSize(ts3);
	  Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetLabelOffset(Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->GetLabelOffset()*3.);

	  Width_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetTitleOffset(1.);
	  Width_PbPb_minus_pp_eta[g][j]->GetYaxis()->CenterTitle();

	}else{
	  Width_PbPb_minus_pp_eta[g][j]->GetYaxis()->SetLabelSize(0);
	}


	Width_PbPb_minus_pp_eta[g][j]->SetMarkerColor(kBlack);
	Width_PbPb_minus_pp_eta[g][j]->SetMinimum(-0.05);
	Width_PbPb_minus_pp_eta[g][j]->SetMaximum(0.26);

	Width_PbPb_minus_pp_eta[g][j]->SetLineWidth(0);

	Width_PbPb_minus_pp_eta[g][j]->Draw("p Z e2");

	Width_PbPb_minus_pp_eta[g][j]->Draw("p Z same");

	cout<<"here"<<endl;
	if(with_ref)	Width_diff_eta_ref[g][j]->Draw("same P");

	
	if(j==0){

	  l40 = new TLegend(textalign2-.05,texty2-0.05,0.6,texty1-0.05);
	  l40->AddEntry(Width_PbPb_minus_pp_eta[g][j],"PbPb - pp","lpfe"); 
	  l40->SetTextFont(43);	    
	  l40->SetTextSizePixels(tspixels);
	  l40->SetLineColor(kWhite);
	  l40->Draw();
	}



	line->SetLineColor(kBlack);
	line->Draw("same");
	
	//-------------------------------------
	//    Start dPhi
	//----------------------------


	TString Width_PbPb_minus_pp_phi_name = "Width_PbPb_minus_pp_phi"; Width_PbPb_minus_pp_phi_name+= g; Width_PbPb_minus_pp_phi_name+= j;

	Width_PbPb_minus_pp_phi[g][j]  = (TGraphAsymmErrors*) Width_phi[g-1][j]->Clone(Width_PbPb_minus_pp_phi_name);

	for(int k = 0; k<4; k++){
	  
	  float temp_x = pTbin_centers.at(k);
	  
	  float PbPb_val = Width_phi[g-1][j]->Eval(temp_x);
	  float pp_val = Width_phi[g][j]->Eval(temp_x);

	  Width_PbPb_minus_pp_phi[g][j]->SetPoint(k,temp_x,PbPb_val-pp_val);


	}



	//-------------------
	//  Draw  WidthPbPb - Widthpp
	//----------------
   
	Width_phi_canvas[g]->cd(j+5);
	Width_PbPb_minus_pp_phi[g][j]->Draw("AP");
	Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetRangeUser(0.8,8.2);
	Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetNdivisions(8);
	Width_PbPb_minus_pp_phi[g][j]->GetYaxis()->SetNdivisions(204);

	switch(g){
	case 1: 
	  Width_PbPb_minus_pp_phi[g][j]->SetMarkerStyle(10);
	  Width_PbPb_minus_pp_phi[g][j]->SetFillColor(90);
	  Width_PbPb_minus_pp_phi[g][j]->SetLineColor(90);
	  break;
	case 3: 
	  Width_PbPb_minus_pp_phi[g][j]->SetMarkerStyle(34);
	  Width_PbPb_minus_pp_phi[g][j]->SetMarkerSize(2);
	  Width_PbPb_minus_pp_phi[g][j]->SetFillColor(30);
	  Width_PbPb_minus_pp_phi[g][j]->SetLineColor(30);
	  break;  
	case 5: 
	  Width_PbPb_minus_pp_phi[g][j]->SetMarkerStyle(21);
	  Width_PbPb_minus_pp_phi[g][j]->SetFillColor(kOrange-2);
	  Width_PbPb_minus_pp_phi[g][j]->SetLineColor(kOrange-2);
	  break;
	}


	Width_PbPb_minus_pp_phi[g][j]->SetMarkerSize(2);
	Width_PbPb_minus_pp_phi[g][j]->SetTitle("");
	Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetLabelSize(ts);
	Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetTitleSize(tstitle-.002);
	Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->CenterTitle();
	Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetTitleOffset(0.82);
	Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV)");


	if(j==0){
	  Width_PbPb_minus_pp_phi[g][j]->GetYaxis()->SetTitle("(#Delta#phi Width)_{PbPb} - (#Delta#phi Width)_{pp}");
	  Width_PbPb_minus_pp_phi[g][j]->GetYaxis()->SetTitleSize(ts);
	  Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetTitleSize(ts-.004);
	  Width_PbPb_minus_pp_phi[g][j]->GetYaxis()->SetLabelSize(ts3);
	  Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetLabelSize(ts3);
	  Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetTitleOffset(1.);
	  Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetLabelOffset( Width_PbPb_minus_pp_phi[g][j]->GetXaxis()->GetLabelOffset()*3.);
	  Width_PbPb_minus_pp_phi[g][j]->GetYaxis()->CenterTitle();

	}else{
	  Width_PbPb_minus_pp_phi[g][j]->GetYaxis()->SetLabelSize(0);
	}

	Width_PbPb_minus_pp_phi[g][j]->SetMarkerColor(kBlack);
	Width_PbPb_minus_pp_phi[g][j]->SetMinimum(-0.05);
	Width_PbPb_minus_pp_phi[g][j]->SetMaximum(0.26);

	Width_PbPb_minus_pp_phi[g][j]->SetLineWidth(0);

	Width_PbPb_minus_pp_phi[g][j]->Draw("p Z e2");

	Width_PbPb_minus_pp_phi[g][j]->Draw("p Z same");

	if(with_ref)	Width_diff_phi_ref[g][j]->Draw("same P");

	line->Draw("same");

	
	if(j==0){

	  l40 = new TLegend(textalign2-.05,texty2-0.05,0.6,texty1-0.05);
	  l40->AddEntry(Width_PbPb_minus_pp_phi[g][j],"PbPb - pp","lpfe"); 
	  l40->SetTextFont(43);	    
	  l40->SetTextSizePixels(tspixels);
	  l40->SetLineColor(kWhite);
	  l40->Draw();
	}

      

	Width_phi_canvas[g]->cd(j+1);
	Width_phi[g-1][j]->Draw("AP");
	Width_phi[g-1][j]->GetXaxis()->SetRangeUser(0.8,8.2);
	Width_phi[g-1][j]->GetXaxis()->SetNdivisions(8);
	Width_phi[g-1][j]->GetYaxis()->SetNdivisions(306);

	switch(g){
	case 1: 
	  Width_phi[g-1][j]->SetMarkerStyle(10);
	  Width_phi[g][j]->SetMarkerStyle(4);
	  if(with_ref) 	  Width_phi_ref[g-1][j]->SetMarkerStyle(10);
	  if(with_ref) 	  Width_diff_phi_ref[g][j]->SetMarkerStyle(10);
	  if(with_ref) 	  Width_phi_ref[g][j]->SetMarkerStyle(4);
	  Width_phi[g-1][j]->SetFillColor(90);
	  Width_phi[g-1][j]->SetLineColor(90);

	  break;
	case 3: 
	  Width_phi[g-1][j]->SetMarkerStyle(34);
	  Width_phi[g-1][j]->SetMarkerSize(2);
	  if(with_ref) 	  Width_phi_ref[g-1][j]->SetMarkerStyle(34);
	  if(with_ref) 	  Width_phi_ref[g-1][j]->SetMarkerSize(2);
	  if(with_ref) 	  Width_diff_phi_ref[g][j]->SetMarkerStyle(34);
	  if(with_ref) 	  Width_diff_phi_ref[g][j]->SetMarkerSize(2);

	  Width_phi[g-1][j]->SetFillColor(30);
	  Width_phi[g-1][j]->SetLineColor(30);
	  Width_phi[g][j]->SetMarkerStyle(28);
	  break;  
	case 5: 
	  Width_phi[g-1][j]->SetMarkerStyle(21);
	  if(with_ref) 	  Width_phi_ref[g-1][j]->SetMarkerStyle(21);
	  if(with_ref) 	  Width_phi_ref[g-1][j]->SetMarkerSize(2);
	  if(with_ref) 	  Width_diff_phi_ref[g][j]->SetMarkerStyle(21);
	  if(with_ref) 	  Width_diff_phi_ref[g][j]->SetMarkerSize(2);
	  Width_phi[g-1][j]->SetFillColor(kOrange-2);
	  Width_phi[g-1][j]->SetLineColor(kOrange-2);
	  Width_phi[g][j]->SetMarkerStyle(25);
	  break;
	}

	Width_phi[g-1][j]->SetTitle("");

	Width_phi[g][j]->SetMarkerSize(3);
	Width_phi[g-1][j]->SetMarkerSize(2);

	if(with_ref) 	Width_phi_ref[g][j]->SetMarkerSize(2);
	if(with_ref) 	Width_phi_ref[g-1][j]->SetMarkerSize(2);

	if(j==0){
	  Width_phi[g-1][j]->GetYaxis()->SetTitle("#Delta#phi Width");
	  Width_phi[g-1][j]->GetYaxis()->SetTitleSize(tstitle);
	  Width_phi[g-1][j]->GetYaxis()->SetLabelSize(ts2);
	  Width_phi[g-1][j]->GetYaxis()->CenterTitle();

	}else{
	  Width_phi[g-1][j]->GetYaxis()->SetLabelSize(0);
	}

	Width_phi[g-1][j]->SetMarkerColor(kBlack);
	Width_phi[g-1][j]->SetMinimum(-0.02);
	Width_phi[g-1][j]->SetMaximum(0.75);

	Width_phi[g-1][j]->SetLineWidth(0);

	Width_phi[g-1][j]->Draw("p Z e2");
	Width_phi[g][j]->SetMarkerSize(2);

	Width_phi[g-1][j]->Draw("p Z same");
	Width_phi[g][j]->Draw("p X same");


	if(with_ref)	Width_phi_ref[g-1][j]->Draw("same P");
	if(with_ref)	Width_phi_ref[g][j]->Draw("same P X");
	
	if(j==0){

	  l40 = new TLegend(textalign2-.05,texty3-0.05,0.6,texty1-0.05);
	  if(g==1){  
	    l40->AddEntry(Width_phi[g-1][j],"PbPb Inclusive Jets","lpfe"); 
	    l40->AddEntry(Width_phi[g][j],"pp Inclusive Jets","p"); 
	  }
	  if(g==3){  
	    l40->AddEntry(Width_phi[g-1][j],"PbPb Subleading Jets","lpfe"); 
	    l40->AddEntry(Width_phi[g][j],"pp Subleading Jets","p"); 
	  }
	  if(g==5){  
	    l40->AddEntry(Width_phi[g-1][j],"PbPb Leading Jets","lpfe"); 
	    l40->AddEntry(Width_phi[g][j],"pp Leading Jets","p"); 

	  }
	  l40->SetTextFont(43);	    
	  l40->SetTextSizePixels(tspixels);
	  l40->SetLineColor(kWhite);
	  l40->Draw();
	}

	//	line->Draw("same");

	drawlabels_int_pt2(g,j);

	  
	if(j==3){
	  tex24phi = new TLatex(textalign,texty2,etarangelabel);
	  tex24phi->SetName("tex24phi");
	  tex24phi->SetNDC();
	  tex24phi->SetTextSizePixels(tspixels);
	  tex24phi->Draw();
	}



	Width_eta[g-1][j]->Write();
	Width_eta[g][j]->Write();
	Width_PbPb_minus_pp_eta[g][j]->Write();


	Width_phi[g-1][j]->Write();
	Width_phi[g][j]->Write();
	Width_PbPb_minus_pp_phi[g][j]->Write();

    

	//*****************************************
	//*****************************************
	
	// Width of difference as function of cent.
	
	//*****************************************
	//*****************************************

      }  //close only pp
      
      cout<<"made it here "<<g<<" "<<" "<<j<<endl;

      for(int i = 0; i<4; i++){
	
	check_fits_eta[g]->cd(4*i+j+1);

	check_new_eta_rebin[g][i][j]->SetAxisRange(-2.499,2.4999);

	check_new_eta_rebin[g][i][j]->SetMarkerColor(kBlack);
	check_new_eta_rebin[g][i][j]->SetLineColor(kBlack);
	check_new_eta_rebin[g][i][j]->SetMarkerStyle(20);
	check_new_eta_rebin[g][i][j]->SetMarkerSize(1);
	check_new_eta_rebin[g][i][j]->Draw();

	if(g==0||g==2||g==4){	
	  check_new_eta_syst[g][i][j]->Draw("same e2");
	  check_new_eta_up[g][i][j]->SetMarkerSize(0.);
	  check_new_eta_up[g][i][j]->Draw("same");
	  check_new_eta_down[g][i][j]->SetMarkerSize(0.);
	  check_new_eta_down[g][i][j]->Draw("same");}
	
	gaus_tot_eta[g][i][j]->SetLineWidth(2);
	gaus_tot_eta[g][i][j]->SetLineColor(kBlue);


	gaus_1_eta[g][i][j]->SetLineWidth(2);
	gaus_1_eta[g][i][j]->SetLineStyle(2);
	gaus_1_eta[g][i][j]->SetLineColor(kBlack);
	gaus_1_eta[g][i][j]->Draw("same");
	
	gaus_2_eta[g][i][j]->SetLineWidth(2);
	gaus_2_eta[g][i][j]->SetLineStyle(2);
	gaus_2_eta[g][i][j]->SetLineColor(kBlack);
	gaus_2_eta[g][i][j]->Draw("same");
	
	
	
	check_new_eta_rebin[g][i][j]->Draw("same");
	gaus_tot_eta[g][i][j]->Draw("same");

      
	if(i==0){	drawlabels(g,i,j); 
	}else if(j==0){

	  TLatex *tex01;
	  switch(i){
	  case 0:
	    tex01 = new TLatex(textalign2,texty1,"1 < p_{T}^{trk}< 2 GeV");
	    break;
	  case 1:
	    tex01 = new TLatex(textalign2,texty1,"2 < p_{T}^{trk}< 3 GeV");
	    break;
	  case 2:
	    tex01 = new TLatex(textalign2,texty1,"3 < p_{T}^{trk}< 4 GeV");
	    break;
	  case 3:
	    tex01 = new TLatex(textalign2,texty1,"4 < p_{T}^{trk}< 8 GeV");
	    break;
	  }
	    
	  tex01->SetName("tex01");
	  tex01->SetNDC();
	  tex01->SetTextSize(ts2);
	  tex01->SetLineWidth(2);
	  tex01->Draw("same");
	}



	check_fits_phi[g]->cd(4*i+j+1);

	check_new_phi_rebin[g][i][j]->SetAxisRange(-1.499,1.4999);



	check_new_phi_rebin[g][i][j]->Draw();

	if(g==0||g==2||g==4){	
	  check_new_phi_syst[g][i][j]->Draw("same e2");
	  check_new_phi_up[g][i][j]->SetMarkerSize(0.);
	  check_new_phi_up[g][i][j]->Draw("same");
	  check_new_phi_down[g][i][j]->SetMarkerSize(0.);
	  check_new_phi_down[g][i][j]->Draw("same");
	}

	gaus_tot_phi[g][i][j]->SetLineWidth(2);
	gaus_tot_phi[g][i][j]->SetLineColor(kBlue);


	gaus_1_phi[g][i][j]->SetLineWidth(2);
	gaus_1_phi[g][i][j]->SetLineStyle(2);
	gaus_1_phi[g][i][j]->SetLineColor(kBlack);
	gaus_1_phi[g][i][j]->Draw("same");
	
	gaus_2_phi[g][i][j]->SetLineWidth(2);
	gaus_2_phi[g][i][j]->SetLineStyle(2);
	gaus_2_phi[g][i][j]->SetLineColor(kBlack);
	gaus_2_phi[g][i][j]->Draw("same");

	check_new_phi_rebin[g][i][j]->GetFunction("double_gaus")->SetLineColor(kBlue);	

	check_new_phi_rebin[g][i][j]->Draw("same");
      
	if(i==0){	drawlabels(g,i,j); 
	}else if(j==0){

	  TLatex *tex01;
	  switch(i){
	  case 0:
	    tex01 = new TLatex(textalign2,texty1,"1 < p_{T}^{trk}< 2 GeV");
	    break;
	  case 1:
	    tex01 = new TLatex(textalign2,texty1,"2 < p_{T}^{trk}< 3 GeV");
	    break;
	  case 2:
	    tex01 = new TLatex(textalign2,texty1,"3 < p_{T}^{trk}< 4 GeV");
	    break;
	  case 3:
	    tex01 = new TLatex(textalign2,texty1,"4 < p_{T}^{trk}< 8 GeV");
	    break;
	  }
	    
	  tex01->SetName("tex01");
	  tex01->SetNDC();
	  tex01->SetTextSize(ts2);
	  tex01->SetLineWidth(2);
	  tex01->Draw("same");
	}

	dummy->cd();

      	
	if((g==1||g==3||g==5)&&i<2){

	  check_fits_diff[g]->cd(4*i+j+1);
	  
	  
	  PbPb_pp_eta[g][i][j]->SetMinimum(-1.);
	  if(i==0){
	    PbPb_pp_eta[g][i][j]->SetMaximum(7.);
	  }  else{
	    PbPb_pp_eta[g][i][j]->SetMaximum(4.);
	  }
	  
	  PbPb_pp_eta[g][i][j]->SetMarkerColor(kBlack);
	  PbPb_pp_eta[g][i][j]->SetLineColor(kBlack);
	  PbPb_pp_eta[g][i][j]->SetMarkerStyle(20);
	  PbPb_pp_eta[g][i][j]->SetMarkerSize(1);

	  PbPb_pp_eta[g][i][j]->Draw();
	  //	  PbPb_pp_eta_syst[g][i][j]->Draw("same e2");


	  PbPb_pp_eta_up[g][i][j]->SetMarkerSize(0);
	  PbPb_pp_eta_down[g][i][j]->SetMarkerSize(0);
	  PbPb_pp_eta_up[g][i][j]->Draw("same");
	  PbPb_pp_eta_down[g][i][j]->Draw("same");

	  PbPb_pp_eta[g][i][j]->GetFunction("double_gaus")->SetLineColor(kBlue);
	  PbPb_pp_eta[g][i][j]->Draw("same");

	  if(i==0){  drawlabels(g,i,j); 
	  }else if(j==0){

	    TLatex *tex01;
	    switch(i){
	    case 0:
	      tex01 = new TLatex(textalign2,texty1,"1 < p_{T}^{trk}< 2 GeV");
	      break;
	    case 1:
	      tex01 = new TLatex(textalign2,texty1,"2 < p_{T}^{trk}< 3 GeV");
	      break;
	    case 2:
	      tex01 = new TLatex(textalign2,texty1,"3 < p_{T}^{trk}<4 GeV");
	      break;
	    case 3:
	      tex01 = new TLatex(textalign2,texty1,"4 < p_{T}^{trk}<8 GeV");
	      break;
	    }
	    
	    tex01->SetName("tex01");
	    tex01->SetNDC();
	    tex01->SetTextSize(ts2);
	    tex01->SetLineWidth(2);
	    tex01->Draw("same");
	  }

	  lineEta = new TLine(-2.,0.,2.,0.);
	  lineEta->SetLineColor(kBlack);
	  lineEta->Draw();

	}

	dummy->cd();
      

      }  //end i

    } // end j

    Width_eta_canvas[g]->cd(0);
    TLatex *canvas_title = new TLatex(0.06,0.95,"CMS ");
    canvas_title->SetTextSizePixels(tspixels);
    canvas_title->SetTextFont(63);
    canvas_title->Draw();

    TLatex *canvas_title2 = new TLatex(0.292,0.95,"PbPb 166 #mub^{-1} (2.76 TeV)               pp 5.3 pb^{-1} (2.76 TeV)");
    canvas_title2->SetTextSizePixels(tspixels);
    canvas_title2->Draw();

    Width_phi_canvas[g]->cd(0);
    canvas_title->Draw();
    canvas_title2->Draw();

    check_fits_eta[g]->cd(0);
    check_fits_phi[g]->cd(0);

    if(g%2==0){
      check_fits_eta[g]->SaveAs((TString)"Width_Check_Fits_Eta_PbPb_"+datalabel+".png");
      check_fits_eta[g]->SaveAs((TString)"Width_Check_Fits_Eta_PbPb_"+datalabel+".eps");
      check_fits_phi[g]->SaveAs((TString)"Width_Check_Fits_Phi_PbPb_"+datalabel+".png");
      check_fits_phi[g]->SaveAs((TString)"Width_Check_Fits_Phi_PbPb_"+datalabel+".eps");
    }else{
      check_fits_eta[g]->SaveAs((TString)"Width_Check_Fits_Eta_pp_"+datalabel+".eps");
      check_fits_diff[g]->SaveAs((TString)"Width_Check_Fits_Diff_"+datalabel+".eps");
      check_fits_eta[g]->SaveAs((TString)"Width_Check_Fits_Eta_pp_"+datalabel+".png");
      check_fits_diff[g]->SaveAs((TString)"Width_Check_Fits_Diff_"+datalabel+".png");
      check_fits_phi[g]->SaveAs((TString)"Width_Check_Fits_Phi_pp"+datalabel+".eps");
      check_fits_diff[g]->SaveAs((TString)"Width_Check_Fits_Diff_"+datalabel+".eps");
      check_fits_phi[g]->SaveAs((TString)"Width_Check_Fits_Phi_pp_"+datalabel+".png");
      check_fits_diff[g]->SaveAs((TString)"Width_Check_Fits_Diff_"+datalabel+".png");
     
      Width_eta_canvas[g]->SaveAs((TString)("Width_Eta_"+datalabel+".pdf"));
      Width_phi_canvas[g]->SaveAs((TString)("Width_Phi_"+datalabel+".pdf"));
      Width_eta_canvas[g]->SaveAs((TString)("Width_Eta_"+datalabel+".png"));
      Width_phi_canvas[g]->SaveAs((TString)("Width_Phi_"+datalabel+".png"));
    }
   
  
  
    //******************************************************

      
    //---------------------------------------
    //Now for width by centrality plots
    //--------------------------------------
      
    cout<<"ready to start widths by centrality"<<endl;

    for(int i = 0; i<2; i++){
      
      Width_val_cent_eta.clear();
      Width_err_cent_eta_up.clear();
      Width_err_cent_eta_down.clear();


      for(int k=0; k<4; k++){
	  
	evalpt = pTbin_centers.at(i);
	Width_val_cent_eta.push_back( Width_eta[g][k]->Eval(evalpt));
	
	eyl = 0.;
	eyh = 0.;
		
	Width_err_cent_eta_up.push_back(eyh);
	Width_err_cent_eta_down.push_back(eyl);

      }
    
      Width_eta_cent[g][i] = new TGraphAsymmErrors(4,&cent_centers[0], &Width_val_cent_eta[0], &cent_errors[0], &cent_errors[0],  &Width_err_cent_eta_down[0], &Width_err_cent_eta_up[0]);

      TString	Width_eta_cent_name = "Width_eta_cent"; 	Width_eta_cent_name+= g; Width_eta_cent_name+= i; 
      Width_eta_cent[g][i]->SetName(Width_eta_cent_name );
    
    
      Width_diff_val_cent_eta.clear();
      Width_diff_err_cent_eta_up.clear();
      Width_diff_err_cent_eta_down.clear();

      cout<<"made it here starting rms eta cent"<<endl;

      if(g%2!=0){
	

	cout<<"Here"<<endl;

	for(int k=0; k<4; k++){
	  
	  evalpt = pTbin_centers.at(i);
	  Width_diff_val_cent_eta.push_back(Width_diff_eta[g][k]->Eval(evalpt));

	  eyl = Width_diff_eta[g][k]->GetErrorYlow(i);
	  eyh = Width_diff_eta[g][k]->GetErrorYhigh(i);
	
	  if(k==0&&i==1&&(g==1||g==5)&&eyh<0.67) eyh = 0.67-Width_diff_eta[g][k]->Eval(evalpt);  //our nominal procedure doesn't 
	  if(k==0&&i==1&&g==5&&eyl<Width_diff_eta[g][k+1]->GetErrorYlow(i)) eyl = Width_diff_eta[g][k+1]->GetErrorYlow(i);  //our nominal procedure doesn't 
		
	  Width_diff_err_cent_eta_up.push_back(eyh);
	  Width_diff_err_cent_eta_down.push_back(eyl);

	  cout<<g<<" "<<i<<" "<<k<<" "<<evalpt<<" "<<Width_diff_eta[g][k]->Eval(evalpt)<<" "<<eyl<<" "<<eyh<<endl;

	}
	Width_diff_eta_cent[g][i] = new TGraphAsymmErrors(4,&cent_centers[0], &Width_diff_val_cent_eta[0], &cent_errors[0], &cent_errors[0],  &Width_diff_err_cent_eta_down[0], &Width_diff_err_cent_eta_up[0]);


	TString	Width_diff_eta_cent_name = "Width_diff_eta_cent"; 	Width_diff_eta_cent_name+= g; Width_diff_eta_cent_name+= i; 
	Width_diff_eta_cent[g][i]->SetName(Width_diff_eta_cent_name );


    
	cWidth_eta_cent[g]->cd(i+1);

	
	TString histnameblank = "blank_hist";
	histnameblank+=g;
	histnameblank+=i;

	blank[i] = new TH1D(histnameblank,"",4,xAxis);
	
	//Markers etc. for all histograms at once
	//----------------------------------------

     

	Width_eta_cent[g-1][i]->SetMarkerSize(2);
	Width_eta_cent[g][i]->SetMarkerSize(2);
	Width_eta_cent[g][i]->SetMarkerColor(1);
	Width_eta_cent[g][i]->SetLineColor(1);

	switch(g){
	case 1:
	  Width_eta_cent[g-1][i]->SetFillColor(90);
	  Width_eta_cent[g-1][i]->SetMarkerStyle(10);
	  Width_eta_cent[g][i]->SetMarkerStyle(4);
	  break;     
	case 3:
	  Width_eta_cent[g-1][i]->SetFillColor(30);
	  Width_eta_cent[g-1][i]->SetMarkerStyle(34);
	  Width_eta_cent[g][i]->SetMarkerStyle(28);
	  break;  
	case 5:
	  Width_eta_cent[g-1][i]->SetFillColor(kOrange-2);
	  Width_eta_cent[g-1][i]->SetMarkerStyle(21);
	  Width_eta_cent[g][i]->SetMarkerStyle(25);

	  break;     
	}

	//Plot aesthetics for every canvas.
	//----------------------------------
	
	blank[i]->SetMinimum(0.);
	blank[i]->SetMaximum(0.9);
	blank[i]->GetXaxis()->SetTitle("Centrality (%)");
	blank[i]->GetXaxis()->SetTitleOffset(1.1);
	blank[i]->GetXaxis()->CenterTitle(true);
	blank[i]->GetXaxis()->SetTitleSize(ts);

	blank[i]->GetYaxis()->SetTitle("#Delta#eta Width");
	blank[i]->GetYaxis()->SetTitleSize(0.);
	blank[i]->GetYaxis()->CenterTitle(true);
	//	blank[i]->GetYaxis()->SetLabelOffset(yoffset);
	blank[i]->GetYaxis()->SetLabelSize(0.);
   
	blank[i]->GetYaxis()->SetTickLength(0.025);

	blank[i]->GetXaxis()->SetBinLabel(1,"50-100");
	blank[i]->GetXaxis()->SetBinLabel(2,"30-50");
	blank[i]->GetXaxis()->SetBinLabel(3,"10-30");
	blank[i]->GetXaxis()->SetBinLabel(4," 0-10");
    
	blank[i]->GetXaxis()->SetLabelSize(0.09);
	blank[i]->GetXaxis()->SetLabelOffset(0.015);
	//	blank[i]->GetXaxis()->SetTicks("+-");
	blank[i]->GetXaxis()->LabelsOption("h");
	blank[i]->GetXaxis()->SetTickLength(0.025);

	switch(i){
	case 0: 
	  gPad->SetLeftMargin(0.1);
	  blank[i]->GetYaxis()->SetTitleSize(0.09);
	  blank[i]->GetXaxis()->SetTitleOffset(1.2);
	  blank[i]->GetYaxis()->SetTitleOffset(0.5);
	  blank[i]->GetXaxis()->SetTitleSize(0.06);
	  blank[i]->SetLabelSize(0.95*blank[i]->GetXaxis()->GetLabelSize());
	  blank[i]->GetYaxis()->SetLabelSize(ts);
	  break;
	case 3:
	  gPad->SetRightMargin(0.02);
	  break;
	default:
	  break;
	}
	blank[i]->Draw();
	Width_eta_cent[g-1][i]->Draw("same p e2");
	Width_eta_cent[g][i]->SetLineWidth(0);
	Width_eta_cent[g][i]->Draw("same p X");
	
	if(i==0){ 
	  l40 = new TLegend(textalign2,texty2-.05,0.6,texty4-0.05);
	  l40->SetName("l40");
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){
	    l40->AddEntry(Width_eta_cent[g-1][i],"Inclusive PbPb","lpfe");
	    l40->AddEntry(Width_eta_cent[g][i],"Inclusive pp","lpe");
	  }
	  if(g==3){
	    l40->AddEntry(Width_eta_cent[g-1][i],"Subleading PbPb","lpfe");
	    l40->AddEntry(Width_eta_cent[g][i],"Subleading pp","lpe");
	  }
	  if(g==5){
	    l40->AddEntry(Width_eta_cent[g-1][i],"Leading PbPb","lpfe");
	    l40->AddEntry(Width_eta_cent[g][i],"Leading pp","lpe");
	  }
	  l40->Draw("same");
	}
	drawlabels_int_cent(g,i);
	 
	lineCent = new TLine(-100.,0.,0.,0.);
	lineCent->SetLineStyle(2);
	lineCent->SetLineWidth(2);
	lineCent->Draw("same");
	
	//---------------------
	//  Width of diff
	//--------------------
      

	cout<<"here starting rms of diff"<<endl;      

	  cWidth_eta_cent[g]->cd(i+1);


	  TString histnameblank2 = "blank2_hist";
	  histnameblank2+=g;
	  histnameblank2+=i;

	  blank2[i] = new TH1D(histnameblank2,"",4,xAxis);
	
	  //Markers etc. for all histograms at once
	  //----------------------------------------


	  Width_diff_eta_cent[g][i]->SetMarkerSize(2);

	  switch(g){
	  case 1:
	    Width_diff_eta_cent[g][i]->SetFillColor(90);
	    Width_diff_eta_cent[g][i]->SetMarkerStyle(10);
	    break;     
	  case 3:
	    Width_diff_eta_cent[g][i]->SetFillColor(30);
	    Width_diff_eta_cent[g][i]->SetMarkerStyle(34);
	    break;  
	  case 5:
	    Width_diff_eta_cent[g][i]->SetFillColor(kOrange-2);
	    Width_diff_eta_cent[g][i]->SetMarkerStyle(21);
	    break;     
	  }

	  //Plot aesthetics for every canvas.
	  //----------------------------------
	
	  blank2[i]->SetMinimum(-0.01);
	  blank2[i]->SetMaximum(1.35);
	  blank2[i]->GetXaxis()->SetTitle("Centrality (%)");
	  blank2[i]->GetXaxis()->CenterTitle(true);
	  blank2[i]->GetXaxis()->SetTitleSize(ts);

	  blank2[i]->GetYaxis()->SetTitle("#Delta#eta Width");
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
	  //	blank2[i]->GetXaxis()->SetTicks("+-");
	  blank2[i]->GetXaxis()->LabelsOption("h");

	  blank2[i]->GetXaxis()->SetTickLength(0.0);
	  blank2[i]->GetYaxis()->SetNdivisions(8);
	

	  if(i==0){
	    gPad->SetLeftMargin(0.1);
	    blank2[i]->GetYaxis()->SetTitleSize(tstitle+0.01);
	    blank2[i]->GetXaxis()->SetTitleOffset(1.);
	    blank2[i]->GetYaxis()->SetTitleOffset(0.5);
	    blank2[i]->GetXaxis()->SetTitleSize(ts);
	    blank2[i]->GetYaxis()->SetLabelSize(tstitle);
	  
	  }

	  blank2[i]->Draw();

	  Width_eta_cent[g-1][i]->SetLineWidth(2);
	  Width_eta_cent[g-1][i]->SetMarkerStyle(1);
	  Width_eta_cent[g-1][i]->SetFillColor(kBlack);
	  Width_diff_eta_cent[g][i]->Draw("same p e2");
	  Width_eta_cent[g-1][i]->Draw("same  ");

	  
	  if(i==0){ 
	    l40 = new TLegend(0.4,texty1+0.05,0.9,texty2);
	    l40->SetName("l40");
	    l40->SetTextSize(ts+0.01);
	    l40->SetFillColor(kWhite);
	    l40->SetLineColor(kWhite);
	    l40->AddEntry(Width_diff_eta_cent[g][i],"Width of PbPb - pp","pfe");
	    l40->AddEntry(Width_eta_cent[g-1][i],"Width of PbPb correlation","l");
	    l40->Draw("same");



	    TLatex *pTtex = new TLatex(0.15,texty2,"1 < p_{T}^{trk}<2 GeV");
	    pTtex->SetTextSizePixels(tspixels);
	    pTtex->SetNDC();
	    pTtex->SetTextFont(43);
	    pTtex->Draw();

	    TLatex *datatex = new TLatex(0.15,texty1,(TString)(datalabel+" Jets"));
	    datatex->SetTextSizePixels(tspixels);
	    datatex->SetTextFont(43);
	    datatex->SetNDC();
	    datatex->Draw();

	  }else{
	  
	    TLatex *pTtex = new TLatex(0.05,texty3,"2 < p_{T}^{trk}< 3");
	    pTtex->SetTextSizePixels(tspixels);
	    pTtex->SetNDC();
	    pTtex->SetTextFont(43);
	    pTtex->Draw();

	    TLatex *tex2 = new TLatex(0.05,texty2,"anti-k_{T} jets,  R=0.3,  |#eta_{jet}|<1.6,  Projected |#Delta#phi|<1");
	    tex2->SetTextSizePixels(tspixels);
	    tex2->SetNDC();
	    tex2->SetTextFont(43);
	    tex2->Draw();

	    if(g<2){
	      TLatex *jettex = new TLatex(0.05,texty1,"p_{T,jet}> 120 GeV");
	      jettex->SetTextSizePixels(tspixels);
	      jettex->SetNDC();
	      jettex->SetTextFont(43);
	      jettex->Draw();
	    }else{
	      TLatex *jettex = new TLatex(0.05,texty1,"p_{T,jet1}> 120 GeV,  p_{T,jet2}> 50 GeV,  |#Delta#phi_{1,2}| > 5#pi/6");
	      jettex->SetTextSizePixels(tspixels);
	      jettex->SetNDC();
	      jettex->SetTextFont(43);
	      jettex->Draw();
	    }



	  }

       
		

	  TLine *line1 = new TLine(-50,-0.05,-50.,0.05);
	  TLine *line2 = new TLine(-30,-0.05,-30.,0.05);
	  TLine *line3 = new TLine(-10,0.-0.05,-10.,0.05);

	  if(g!=5){
	    line1->Draw();
	    line2->Draw();
	    line3->Draw();
	  }


	  gPad->RedrawAxis();	


	  if(g==5){
	    cWidth_eta_cent[g]->cd(i+3);

	    TString histnameblank3 = "blank3_hist";
	    histnameblank3+=g;
	    histnameblank3+=i;
	  
	    blank3[i] = (TH1D*)blank2[i]->Clone(histnameblank3);

	    if(i==0){ 
	      gPad->SetLeftMargin(0.1);
	      blank3[i]->GetYaxis()->SetTitleSize(tstitle);
	      blank3[i]->GetYaxis()->SetTitleOffset(0.55);
	      blank3[i]->GetYaxis()->SetLabelSize(ts);
	      blank3[i]->Draw();
	    }
	    blank3[i]->Draw();
	    Width_diff_eta_cent[3][i]->Draw("same p e2");
	    Width_eta_cent[2][i]->Draw("same  ");

	    line1->Draw();
	    line2->Draw();
	    line3->Draw();

	    // lineCent->Draw("same");

	  
	    if(i==0){ 
	      l40 = new TLegend(0.4,texty1+0.05,0.9,texty2);
	      l40->SetName("l40");
	      l40->SetTextSize(ts2);
	      l40->SetFillColor(kWhite);
	      l40->SetLineColor(kWhite);
	      l40->AddEntry(Width_diff_eta_cent[3][i],"Width of PbPb - pp","pfe");
	      l40->AddEntry(Width_eta_cent[2][i],"Width of PbPb correlation","l");
	      l40->Draw("same");



	      TLatex *datatex2 = new TLatex(0.15,texty1,"Subleading Jets");
	      datatex2->SetTextSizePixels(tspixels);
	      datatex2->SetNDC();
	      datatex2->SetTextFont(43);
	      datatex2->Draw();


	      TLatex *zero = new TLatex(0.0805,.98,"0");
	      zero->SetNDC();
	      zero->SetTextFont(43);
	      zero->SetTextSizePixels(tspixels);
	      zero->Draw();



	    }
	  

	    gPad->RedrawAxis();
	  
	    cWidth_eta_cent[g]->cd(i+1);
	    //	  blank2[i]->Draw("same");

	  gPad->RedrawAxis();

	}

      }
           		
    } //close i


    if(g==1||g==3||g==5){
      cWidth_eta_cent[g]->cd(0);
      
      if(g==5){
	TLatex *canvas_title = new TLatex(0.06,0.95,"CMS ");
	canvas_title->SetTextSizePixels(tspixels);
	canvas_title->SetTextFont(63);
	canvas_title->Draw();

	TLatex *canvas_title2 = new TLatex(0.28,0.95,"PbPb 166 #mub^{-1} (2.76 TeV)                       pp 5.3 pb^{-1} (2.76 TeV)");
	canvas_title2->SetTextSizePixels(tspixels);
	canvas_title2->SetTextFont(43);
	canvas_title2->Draw();
	
      }else{
	TLatex *canvas_title = new TLatex(0.06,0.9,"CMS ");
	canvas_title->SetTextSizePixels(tspixels);
	canvas_title->SetTextFont(63);
	canvas_title->Draw();

	TLatex *canvas_title2 = new TLatex(0.28,0.9,"PbPb 166 #mub^{-1} (2.76 TeV)                       pp 5.3 pb^{-1} (2.76 TeV)");
	canvas_title2->SetTextSizePixels(tspixels);
	canvas_title2->SetTextFont(43);
	canvas_title2->Draw();


      }

      TString WidthSaveName_cent = in_name;
      WidthSaveName_cent.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      WidthSaveName_cent.ReplaceAll("Yield_pp","Width_Eta_Cent");
      WidthSaveName_cent +=".pdf";
      cWidth_eta_cent[g]->SaveAs(WidthSaveName_cent);
    
      WidthSaveName_cent.ReplaceAll("pdf","png");
      WidthSaveName_cent.ReplaceAll("Width","WidthOfDifference");
      cWidth_eta_cent[g]->SaveAs(WidthSaveName_cent);
    
      //fout[g]->cd();
      //   fout[g]->Write();
    
    } 

  } //close g 

  cout<<" this should be the last thing...."<<endl;

  return 0;

} //Close main loop

