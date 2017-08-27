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
//PLOTTING ONLY CODE FOR THE TEAM
//******************************************************

/*
 NOTES AND INSTRUCTIONS: 

 This code is designed to produce standard PAS plots, with one optional reference controlled by the do_ref flag.  

 Jet shapes are not normalized in the files, so that it is possible to make non-normalized versions.  This code by default includes normalization.

*/

Int_t results_plotting(bool is_number=0,bool do_ref=kFALSE){

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins, limR;
  float deta, dphi, r, bc, bg_err,temp1,temp2, rbin, temperr, err, width_temp_x, width_temp_y, width_temp, norm_temp, zerobin, temp, norm_tot, err_temp, cont;
    
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 10;

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300","" };
  TString TrkPtBin_strs2[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300","" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};

  float mean_pts[nTrkPtBins] = {0.844,1.35,2.35,3.37,5.07,9.72,13.8,17.9,22.};

  float RBins[20] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.,1.2,1.4,1.6,1.8,2.0};
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadTopMargin   (0.25);
  gStyle->SetPadLeftMargin  (0.2);
  gStyle->SetPadRightMargin (0.05);
 
  TH1D* JetShape[6][nCBins][nTrkPtBins];
  TH1D* JetShape_noerr_up[6][nCBins][nTrkPtBins];
  TH1D* JetShape_noerr_down[6][nCBins][nTrkPtBins];
  
  TH1D* JetShape_diff_noerr_up[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff_noerr_down[6][nCBins][nTrkPtBins];

  TH1D* JetShape_ref[6][nCBins][nTrkPtBins];
  
  TH1D* JetShape_ref_diff[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff[6][nCBins][nTrkPtBins];

  TH1D* JetShape_ref_ratio[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ratio[6][nCBins][nTrkPtBins];
  
  TH1D* JetShape_syst[6][nCBins][nTrkPtBins];
  TH1D* JetShape_syst_ref[6][nCBins][nTrkPtBins];

  THStack *JetShape_Stack_Up[6][nCBins];
  THStack *JetShape_Diff_Stack_Up[6][nCBins];

  THStack *JetShape_Stack_Down[6][nCBins];
  THStack *JetShape_Diff_Stack_Down[6][nCBins];


  THStack *signal_dEta_stack[12][nCBins];
  THStack *signal_dPhi_stack[12][nCBins];

  THStack *signal_dEta_diff_stack_up[12][nCBins];
  THStack *signal_dEta_diff_stack_down[12][nCBins];
  THStack *signal_dPhi_diff_stack_up[12][nCBins];
  THStack *signal_dPhi_diff_stack_down[12][nCBins];


  TH1D *signal_dPhi[12][nCBins][nTrkPtBins];
  TH1D *signal_dPhi_syst[12][nCBins][nTrkPtBins];
  TH1D *signal_dPhi_rebin[12][nCBins][nTrkPtBins];
  TH1D *signal_dPhi_PbPb_pp[12][nCBins][nTrkPtBins];  
  TH1D *signal_dPhi_PbPb_pp_syst[12][nCBins][nTrkPtBins];  
  TH1D *signal_dPhi_PbPb_pp_up[12][nCBins][nTrkPtBins];  
  TH1D *signal_dPhi_PbPb_pp_down[12][nCBins][nTrkPtBins];  

  TH1D *signal_dEta[12][nCBins][nTrkPtBins];
  TH1D *signal_dEta_syst[12][nCBins][nTrkPtBins];
  TH1D *signal_dEta_rebin[12][nCBins][nTrkPtBins];
  TH1D *signal_dEta_PbPb_pp[12][nCBins][nTrkPtBins];  
  TH1D *signal_dEta_PbPb_pp_syst[12][nCBins][nTrkPtBins];  
  TH1D *signal_dEta_PbPb_pp_up[12][nCBins][nTrkPtBins];  
  TH1D *signal_dEta_PbPb_pp_down[12][nCBins][nTrkPtBins];  



  double temp_cont, nextr, nextl, cos_weight,me00_range, mc_error,me_error, bg_error, me_error2, bg_error2;
  
  TString stem, ref_stem, datalabel,me00_range_string,stem_mc;
 
  float norm;
 
  TH1D *Integral_Pt[12][5];
  TH1D *Integral_syst_Pt[12][5];
  TH1D *Integral_diff_Pt[12][5];
  TH1D *Integral_diff_syst_Pt[12][5];
  TGraphAsymmErrors *Integral_diff_graph_ref[12][5];

  double integral[12][nCBins][nTrkPtBins];
  double integral_err[12][nCBins][nTrkPtBins];
  double integral_syst_err[12][nCBins][nTrkPtBins];

  TFile *f_in_ref_int = new TFile("Inclusive_Data_AllPlots.root");

  TPaveText *labels;

  //-----------------------------------
  // Pick input and reference files here
  //----------------------------------
  
 
  TFile *f_in, *f_in_ref, *f_in_eta_phi;

  //HALLIE'S FILE PATHS FOR EASY ITERATION
  

  if(is_number){
    f_in = new TFile("../jet_shapes_result/Jet_Shapes.root");
    f_in_eta_phi = new TFile("../particle_yields/Particle_Yields.root");
    f_in_ref = new TFile("Jet_Shapes_QM.root");
    //  f_in_ref = new TFile("Jet_Shapes_Preapproval.root");

  }else{

    f_in = new TFile("../jet_shapes_result/Jet_Shapes_pTweighted.root");
    f_in_ref = new TFile("Jet_Shapes_pTweighted_QM.root");
    //   f_in_ref = new TFile("Jet_Shapes_pTweighted_Preapproval.root");
    //  f_in_ref = new TFile("Jet_Shapes_Run1.root");

  }

  //FILE PATHS IN GITHUB CODE:
  /*
 if(is_number){
    f_in = new TFile("Jet_Shapes.root");
    f_in_eta_phi = new TFile("Particle_Yields.root");
    f_in_ref = new TFile("Jet_Shapes_Preapproval.root");

  }else{

    f_in = new TFile("Jet_Shapes_pTweighted.root");
    //  f_in_ref = new TFile("Jet_Shapes_pTweighted_Preapproval.root");
    f_in_ref = new TFile("Jet_Shapes_Run1.root");

  }
  */
  TF1 *temp_func = new TF1("temp_func","[0]"); // this exists because you can't scale a graph normally
  //-----------------------------------
  // Consolidated style settings
  //----------------------------------
  

  float label_font = 42;
  float x_label_offset = 0.02;
  float x_label_size = 0.08;
  float x_title_size = 0.08;
  float x_title_offset = .9;
  float x_tick_length = 0.025;

  float y_label_size = 0.08;
  float y_title_size = 0.1;
  float y_label_offset = 0.004;
  float y_title_offset = 0.9;
  float y_tick_length = 0.025;


  float r_max = 0.999;

  float number_y_min = -5.;
  float number_y_max = 43.;
  
  float jetshape_y_min_normed = .005;
  float jetshape_y_max_normed = 23.5;
  
 
  float jetshape_y_min = .5;
  float jetshape_y_max = 2000.;
 
  float jetshape_ratio_y_min = 0.;
  float jetshape_ratio_y_max = 3.2;

  //  if(do_ref)jetshape_ratio_y_max = 6.5;

  float number_diff_y_min = -1.5;
  float number_diff_y_max = 15.;

  float number_diff_y_min_r = -5.5;
  float number_diff_y_max_r = 19.5;

  float integral_max = 13.;
  float integral_min = -1.;

  float integral_diff_max = 8.2;
  float integral_diff_min = -1.;


 
  //***********************************
  //***********************************

  //-----------------------
  // Start getting histos
  //-----------------------

    
  if(is_number) stem = "JetShape2_Yield_BkgSub_Inclusive_";
  else stem = "JetShape2_Yield_BkgSub_pTweightedInclusive_";

  ref_stem = "Jet_Shape_SystErr_";


  for(int g=0; g<2; g++){
  
    for (int ibin=0;ibin<nCBins;ibin++){

      if(g==1&&ibin > 0) continue;
    
	for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){


	  cout<<(TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1])<<endl;
	  if(g==0){
	    JetShape[g][ibin][ibin3] = (TH1D*)f_in->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + "_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    cout<<"0"<<ibin<<" "<<ibin3<<endl;
	    if(do_ref)    JetShape_ref[g][ibin][ibin3] = (TH1D*)f_in_ref->Get((TString)(ref_stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + "_Ref_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
		  
	    JetShape_syst[g][ibin][ibin3] = (TH1D*)f_in->Get((TString)("Jet_Shape_SystErr_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Jet_Shape_SystErr_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    cout<<"1"<<endl;
	 
	    if(is_number&&ibin3<8){

	      signal_dPhi_rebin[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("Proj_dPhi_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]+"_Rebin"))->Clone((TString) ("Proj_dPhi_PbPb_"  + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 
	      
	      signal_dEta_rebin[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("Proj_dEta_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]+"_Rebin"))->Clone((TString) ("Proj_dEta_PbPb_"  + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 

	      signal_dPhi_syst[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("dPhi_Syst_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("dPhi_Syst_PbPb_"  + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 
	      signal_dEta_syst[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("dEta_Syst_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("dEta_Syst_PbPb_"  + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 


	    }

	  
	  }else{
	  

	    JetShape[g][ibin][ibin3] = (TH1D*)f_in->Get((TString)(stem +"pp_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + "_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    
	    cout<<"0 for pp"<<endl;
	    

	    cout<<(TString)(stem +"pp_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1])<<endl;
	    if(do_ref)	    JetShape_ref[g][ibin][ibin3] = (TH1D*)f_in_ref->Get((TString)("Jet_Shape_SystErr_pp_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + "_Ref_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    
	    
	    JetShape_syst[g][ibin][ibin3] = (TH1D*)f_in->Get((TString)("Jet_Shape_SystErr_pp_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Jet_Shape_SystErr_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    cout<<"1"<<endl;

	    if(is_number&&ibin3<8){
	      signal_dPhi_rebin[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("Proj_dPhi_pp_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]+"_Rebin"))->Clone((TString) ("Proj_dPhi_pp_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 
	      signal_dEta_rebin[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("Proj_dEta_pp_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]+"_Rebin"))->Clone((TString) ("Proj_dEta_pp_"  +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 

	      signal_dPhi_syst[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("dPhi_Syst_pp_Cent0_Cent10_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("dPhi_Syst_pp_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 

	      signal_dEta_syst[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("dEta_Syst_pp_Cent0_Cent10_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("dEta_Syst_pp_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 



	    }
	  
	    cout<<"2"<<endl;
	  }
	

	  JetShape_syst[g][ibin][ibin3]->SetMarkerStyle(20);
	  JetShape_syst[g][ibin][ibin3]->SetMarkerSize(1);
	  JetShape_syst[g][ibin][ibin3]->SetMarkerColor(kWhite);
	  JetShape_syst[g][ibin][ibin3]->SetFillColor(kBlack);
	  JetShape_syst[g][ibin][ibin3]->SetFillStyle(3004);
	
	
	  JetShape[g][ibin][ibin3]->GetXaxis()->SetLabelFont(label_font);
	  JetShape[g][ibin][ibin3]->GetXaxis()->SetLabelOffset(x_label_offset);
	  JetShape[g][ibin][ibin3]->GetXaxis()->SetLabelSize(x_label_size);
	  JetShape[g][ibin][ibin3]->GetXaxis()->SetTitleSize(x_title_size);
	  JetShape[g][ibin][ibin3]->GetXaxis()->SetTickLength(x_tick_length);
	  JetShape[g][ibin][ibin3]->GetXaxis()->SetTitleOffset(x_title_offset);
	  JetShape[g][ibin][ibin3]->GetXaxis()->SetTitleFont(label_font);
	  JetShape[g][ibin][ibin3]->GetYaxis()->SetTitle("#Rho(#Deltar) (GeV)");
	  JetShape[g][ibin][ibin3]->GetYaxis()->CenterTitle(true);
	  JetShape[g][ibin][ibin3]->GetYaxis()->SetLabelFont(label_font);
	  JetShape[g][ibin][ibin3]->GetYaxis()->SetLabelOffset(y_label_offset);
	  JetShape[g][ibin][ibin3]->GetYaxis()->SetLabelSize(y_label_size);
	  JetShape[g][ibin][ibin3]->GetYaxis()->SetTitleSize(y_title_size);
	  JetShape[g][ibin][ibin3]->GetYaxis()->SetTickLength(y_tick_length);
	  JetShape[g][ibin][ibin3]->GetYaxis()->SetTitleOffset(y_title_offset);
	  JetShape[g][ibin][ibin3]->GetYaxis()->SetLabelSize(y_label_size);

	  JetShape[g][ibin][ibin3]->SetAxisRange(0.,r_max);
	  JetShape[g][ibin][ibin3]->SetMarkerSize(1);
	  JetShape[g][ibin][ibin3]->SetLineColor(kBlack);
	  JetShape[g][ibin][ibin3]->SetMarkerColor(kBlack);
	  JetShape[g][ibin][ibin3]->SetMarkerStyle(24);

	  cout<<"3"<<endl;
	  if(do_ref){
	    JetShape_ref[g][ibin][ibin3]->SetMarkerSize(2);
	    if(is_number)   JetShape_ref[g][ibin][ibin3]->SetLineColor(kRed);
	    else    JetShape_ref[g][ibin][ibin3]->SetLineColor(kSpring);
	    
	    if(is_number)	    JetShape_ref[g][ibin][ibin3]->SetMarkerColor(kRed);
	    else 	    JetShape_ref[g][ibin][ibin3]->SetMarkerColor(kSpring);
       
	    if(is_number)	    JetShape_ref[g][ibin][ibin3]->SetFillColor(kRed);
	    else 	    JetShape_ref[g][ibin][ibin3]->SetFillColor(kSpring-1);
	
	    JetShape_ref[g][ibin][ibin3]->SetMarkerStyle(34);
	    JetShape_ref[g][ibin][ibin3]->SetFillStyle(3005);
	  
	  }

	  JetShape[g+3][ibin][ibin3] = (TH1D*) JetShape[g][ibin][ibin3]->Clone((TString) (stem + "_NormedJetShape_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  JetShape_syst[g+3][ibin][ibin3] = (TH1D*) JetShape_syst[g][ibin][ibin3]->Clone((TString) (stem + "_NormedJetShape_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  

	  if(!is_number&&ibin3==9){

	    cout<<"4"<<endl;
	    if(do_ref){
	      JetShape_ref[g+3][ibin][ibin3] = (TH1D*) JetShape_ref[g][ibin][ibin3]->Clone((TString) (stem + "_NormedJetShape_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	 
	      cout<<g<<endl;
	      /*
	      if(g==0)	      JetShape_syst_ref[g+3][ibin][ibin3] = (TH1D*)f_in_ref->Get((TString)("Jet_Shape_SystErr_"+ CBin_strs[ibin]  +"_"+ CBin_strs[ibin+1] + "_TrkPt300_"))->Clone((TString) (stem + "_Syst_PbPb_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs2[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	      else 	      JetShape_syst_ref[g+3][ibin][ibin3] = (TH1D*)f_in_ref->Get((TString)("Jet_Shape_SystErr_pp_TrkPt300_"))->Clone((TString) (stem + "_Syst_pp_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs2[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	      */
	    }

	    cout<<"5"<<endl;
	    norm = JetShape[g][ibin][ibin3]->Integral("width");
	  
	    for(int k = 0; k<10; k++){
	      JetShape[g+3][ibin][k]->Scale(1./norm);
	      JetShape_syst[g+3][ibin][k]->Scale(1./norm);
	    }

	    if(do_ref){
	  
	      norm = JetShape_ref[g][ibin][ibin3]->Integral("width");
	
	      JetShape_ref[g+3][ibin][ibin3]->Scale(1./norm);
	      // JetShape_syst_ref[g+3][ibin][ibin3]->Scale(1./norm);
	    }
	    
	  }
	}
    }
  }

  cout<<"got all histograms"<<endl;


  for (int ibin=0;ibin<nCBins;ibin++){ 
    
    for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

      for(int g = 0; g<4; g++){

	if(g==1||g==2)continue;
	if(is_number&&g>0)continue;
	TString norm_string = "PerJet_";
	if(g==3) norm_string="NormalizedToOne_";
      
	if(is_number){
	
	  JetShape_syst[g+2][ibin][ibin3] = (TH1D*) JetShape_syst[g][ibin][ibin3]->Clone((TString) ("Jet_Shape_Syst_Diff_" +norm_string+ CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  JetShape_syst[g+2][ibin][ibin3]->Add(JetShape_syst[g+1][0][ibin3],-1.);
	}else{
	  JetShape_syst[g+2][ibin][ibin3] = (TH1D*) JetShape_syst[g][ibin][ibin3]->Clone((TString) ("Jet_Shape_Syst_Ratio_" +norm_string+ CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  JetShape_syst[g+2][ibin][ibin3]->Divide(JetShape_syst[g+1][0][ibin3]);

	}

	JetShape_syst[g+2][ibin][ibin3]->SetMarkerStyle(20);
	JetShape_syst[g+2][ibin][ibin3]->SetMarkerSize(1);
	  JetShape_syst[g+2][ibin][ibin3]->SetMarkerColor(kWhite);
	  JetShape_syst[g+2][ibin][ibin3]->SetFillColor(kBlack);
	  JetShape_syst[g+2][ibin][ibin3]->SetFillStyle(3004);
	  TString diff_name ="JetShape_diff_"; diff_name+=norm_string; diff_name+=CBin_strs[ibin]; diff_name+= "_"; diff_name += CBin_strs[ibin+1]; diff_name+= ibin3;
	  JetShape_diff[g][ibin][ibin3] = (TH1D*)JetShape[g][ibin][ibin3]->Clone(diff_name);
	  JetShape_diff[g][ibin][ibin3]->Add( JetShape[g+1][0][ibin3],-1. );
    
	  TString ratio_name ="JetShape_ratio_"; ratio_name+=norm_string; ratio_name+=CBin_strs[ibin]; ratio_name+= "_"; ratio_name += CBin_strs[ibin+1]; ratio_name+= ibin3;

      
	  JetShape_ratio[g][ibin][ibin3] = (TH1D*)JetShape[g][ibin][ibin3]->Clone(ratio_name);
	  JetShape_ratio[g][ibin][ibin3]->Divide(JetShape[g+1][0][ibin3]);
	  JetShape_ratio[g][ibin][ibin3]->GetXaxis()->SetNdivisions(505);
  
	  if(!is_number&&do_ref&&ibin3==9){
	    JetShape_ref_ratio[g][ibin][ibin3] = (TH1D*)JetShape_ref[g][ibin][ibin3]->Clone((TString)(ratio_name+"_Ref"));
	    JetShape_ref_ratio[g][ibin][ibin3]->Divide( JetShape_ref[g+1][0][ibin3]);
	    /*
	    if(g==3){
	      JetShape_syst_ref[g+2][ibin][ibin3] = (TH1D*)JetShape_syst_ref[g][ibin][ibin3]->Clone((TString)(ratio_name+"_SystRef"));

	      JetShape_syst_ref[g+2][ibin][ibin3]->Divide( JetShape_ref[g+1][0][ibin3]);

       
	      JetShape_syst_ref[g+2][ibin][ibin3]->SetMarkerStyle(34);
	      JetShape_syst_ref[g+2][ibin][ibin3]->SetMarkerSize(2);
	      JetShape_syst_ref[g+2][ibin][ibin3]->SetMarkerColor(kSpring);
	      JetShape_syst_ref[g+2][ibin][ibin3]->SetFillColor(kYellow-9);
	      JetShape_syst_ref[g+2][ibin][ibin3]->SetFillStyle(3001);

       
	      JetShape_syst_ref[g+1][0][ibin3]->SetMarkerStyle(34);
	      JetShape_syst_ref[g+1][0][ibin3]->SetMarkerSize(2);
	      JetShape_syst_ref[g+1][0][ibin3]->SetMarkerColor(kSpring);
	      JetShape_syst_ref[g+1][0][ibin3]->SetFillColor(kYellow-9);
	      JetShape_syst_ref[g+1][0][ibin3]->SetFillStyle(3001);


       
	      JetShape_syst_ref[g][ibin][ibin3]->SetMarkerStyle(34);
	      JetShape_syst_ref[g][ibin][ibin3]->SetMarkerSize(2);
	      JetShape_syst_ref[g][ibin][ibin3]->SetMarkerColor(kSpring);
	      JetShape_syst_ref[g][ibin][ibin3]->SetFillColor(kYellow-4);
	      JetShape_syst_ref[g][ibin][ibin3]->SetFillStyle(3001);

	    }
	    */
	  }
      }
    }
  }

  cout<<"starting stack plots"<<endl;
  //------------------------
  // STACK PLOTS
  //------------------------


  for(int ibin = 0; ibin<nCBins; ibin++){
    
    for(int g = 0; g<5; g++){
      if(is_number&&g>1)continue;
      if(g==2)continue;
      if((g==1||g==4)&&ibin>0)continue;

 
      TString type_string = "Jet_Shape_";
      switch(g){
      case 0:
	type_string+="PbPb_PerJet_";
	break;
      case 1:
	type_string+="pp_PerJet_";
	break;
    
      case 3:
	type_string+="PbPb_NormedToOne_";
	break;
      case 4:
	type_string+="pp_NormedToOne_";
	break;
    
      }

      for(int k = 0; k<nTrkPtBins-1; k++){
	JetShape_noerr_up[g][ibin][k] = new TH1D((TString)(type_string+"NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
	JetShape_noerr_down[g][ibin][k] = new TH1D((TString)(type_string+"NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
	if(g==0||g==3){   
	  JetShape_diff_noerr_up[g][ibin][k] = new TH1D((TString)(type_string+"Diff_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
	  JetShape_diff_noerr_down[g][ibin][k] = new TH1D((TString)(type_string+"Diff_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
	}

	for(int l = 0; l< JetShape[g][ibin][k]->GetNbinsX()+1; l++){
      
	  bc = JetShape[g][ibin][k]->GetBinContent(l);
          
	  if(bc>0){ 
	    JetShape_noerr_up[g][ibin][k]->SetBinContent(l,bc);	
	    JetShape_noerr_down[g][ibin][k]->SetBinContent(l,0.);	

	  }else{
	    JetShape_noerr_down[g][ibin][k]->SetBinContent(l,bc);	
	    JetShape_noerr_up[g][ibin][k]->SetBinContent(l,0.);	
	  }

	  if(g==0||g==3){  
	    bc = JetShape_diff[g][ibin][k]->GetBinContent(l);
    
	    if(bc>0){ 
	      JetShape_diff_noerr_up[g][ibin][k]->SetBinContent(l,bc);	
	      JetShape_diff_noerr_down[g][ibin][k]->SetBinContent(l,0.);	

	    }else{
	      JetShape_diff_noerr_down[g][ibin][k]->SetBinContent(l,bc);	
	      JetShape_diff_noerr_up[g][ibin][k]->SetBinContent(l,0.);	
	    } 
	  }
	}
	JetShape_noerr_up[g][ibin][k]->SetLineColor(kBlack);;
	JetShape_noerr_down[g][ibin][k]->SetLineColor(kBlack);;
	if(g==0||g==3)	JetShape_diff_noerr_up[g][ibin][k]->SetLineColor(kBlack);;
	if(g==0||g==3)	JetShape_diff_noerr_down[g][ibin][k]->SetLineColor(kBlack);;

	switch(k){
	case 0:
	  JetShape_noerr_up[g][ibin][k]->SetFillColor(kBlue-9);
	  JetShape_noerr_down[g][ibin][k]->SetFillColor(kBlue-9);
	  if(g==0||g==3) JetShape_diff_noerr_up[g][ibin][k]->SetFillColor(kBlue-9);
	  if(g==0||g==3) JetShape_diff_noerr_down[g][ibin][k]->SetFillColor(kBlue-9);
	  break;
	case 1:
	  JetShape_noerr_up[g][ibin][k]->SetFillColor(kYellow-9);
	  JetShape_noerr_down[g][ibin][k]->SetFillColor(kYellow-9);
	  if(g==0||g==3) JetShape_diff_noerr_up[g][ibin][k]->SetFillColor(kYellow-9);
	  if(g==0||g==3) JetShape_diff_noerr_down[g][ibin][k]->SetFillColor(kYellow-9);
	  break;
	case 2:
	  JetShape_noerr_up[g][ibin][k]->SetFillColor(kOrange+1);
	  JetShape_noerr_down[g][ibin][k]->SetFillColor(kOrange+1);
	  if(g==0||g==3) JetShape_diff_noerr_up[g][ibin][k]->SetFillColor(kOrange+1);
	  if(g==0||g==3) JetShape_diff_noerr_down[g][ibin][k]->SetFillColor(kOrange+1);
	  break;
	case 3:
	  JetShape_noerr_up[g][ibin][k]->SetFillColor(kViolet-5);
	  JetShape_noerr_down[g][ibin][k]->SetFillColor(kViolet-5);
	  if(g==0||g==3) JetShape_diff_noerr_up[g][ibin][k]->SetFillColor(kViolet-5);
	  if(g==0||g==3) JetShape_diff_noerr_down[g][ibin][k]->SetFillColor(kViolet-5);
	  break;
	case 4:
	  JetShape_noerr_up[g][ibin][k]->SetFillColor(kGreen+3);
	  JetShape_noerr_down[g][ibin][k]->SetFillColor(kGreen+3);
	  if(g==0||g==3) JetShape_diff_noerr_up[g][ibin][k]->SetFillColor(kGreen+3);
	  if(g==0||g==3) JetShape_diff_noerr_down[g][ibin][k]->SetFillColor(kGreen+3);
	  break;


	case 5:
	  JetShape_noerr_up[g][ibin][k]->SetFillColor(kRed);
	  JetShape_noerr_down[g][ibin][k]->SetFillColor(kRed);
	  if(g==0||g==3) JetShape_diff_noerr_up[g][ibin][k]->SetFillColor(kRed);
	  if(g==0||g==3) JetShape_diff_noerr_down[g][ibin][k]->SetFillColor(kRed);
	  break;

    
	case 6:
	  JetShape_noerr_up[g][ibin][k]->SetFillColor(kRed+1);
	  JetShape_noerr_down[g][ibin][k]->SetFillColor(kRed+1);
	  if(g==0||g==3) JetShape_diff_noerr_up[g][ibin][k]->SetFillColor(kRed+1);
	  if(g==0||g==3) JetShape_diff_noerr_down[g][ibin][k]->SetFillColor(kRed+1);
	  break;

	case 7:
	  JetShape_noerr_up[g][ibin][k]->SetFillColor(kRed+2);
	  JetShape_noerr_down[g][ibin][k]->SetFillColor(kRed+2);
	  if(g==0||g==3) JetShape_diff_noerr_up[g][ibin][k]->SetFillColor(kRed+2);
	  if(g==0||g==3) JetShape_diff_noerr_down[g][ibin][k]->SetFillColor(kRed+2);
	  break;

	case 8:
	  JetShape_noerr_up[g][ibin][k]->SetFillColor(kRed+3);
	  JetShape_noerr_down[g][ibin][k]->SetFillColor(kRed+3);
	  if(g==0||g==3) JetShape_diff_noerr_up[g][ibin][k]->SetFillColor(kRed+3);
	  if(g==0||g==3) JetShape_diff_noerr_down[g][ibin][k]->SetFillColor(kRed+3);
	  break;

	default:
	  break;
	}

	JetShape_noerr_up[g][ibin][k]->SetFillStyle(1001);
	if(g==0||g==3) JetShape_diff_noerr_up[g][ibin][k]->SetFillStyle(1001);
	JetShape_noerr_down[g][ibin][k]->SetFillStyle(1001);
	if(g==0||g==3) JetShape_diff_noerr_down[g][ibin][k]->SetFillStyle(1001);


      
      } //k
    
 
      JetShape_Stack_Up[g][ibin]=new THStack((TString)("JetShapeStack_PbPb_Up_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
      JetShape_Diff_Stack_Up[g][ibin]=new THStack((TString)("JetShapeStack_Diff_Up_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");

      JetShape_Stack_Down[g][ibin]=new THStack((TString)("JetShapeStack_PbPb_Down_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
      JetShape_Diff_Stack_Down[g][ibin]=new THStack((TString)("JetShapeStack_Diff_Down_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");

      if(is_number){
	for(int k = 0; k<nTrkPtBins-2; k++){
	  JetShape_Stack_Up[g][ibin]->Add(JetShape_noerr_up[g][ibin][nTrkPtBins-3-k]);
	  if(g==0) JetShape_Diff_Stack_Up[g][ibin]->Add(JetShape_diff_noerr_up[g][ibin][nTrkPtBins-3-k]);

	  JetShape_Stack_Down[g][ibin]->Add(JetShape_noerr_down[g][ibin][nTrkPtBins-3-k]);
	  if(g==0)	  JetShape_Diff_Stack_Down[g][ibin]->Add(JetShape_diff_noerr_down[g][ibin][nTrkPtBins-3-k]);
	}

      }else{

	for(int k = 0; k<nTrkPtBins-1; k++){
	  JetShape_Stack_Up[g][ibin]->Add(JetShape_noerr_up[g][ibin][k]);
	  if(g==0||g==3) JetShape_Diff_Stack_Up[g][ibin]->Add(JetShape_diff_noerr_up[g][ibin][k]);

	  JetShape_Stack_Down[g][ibin]->Add(JetShape_noerr_down[g][ibin][k]);
	   if(g==0||g==3) JetShape_Diff_Stack_Down[g][ibin]->Add(JetShape_diff_noerr_down[g][ibin][k]);
	}

      }
  
      if(is_number){

	JetShape_Stack_Up[g][ibin]->SetMaximum(number_y_max);
	JetShape_Stack_Up[g][ibin]->SetMinimum(number_y_min);
  
      }else if(g<3){

	JetShape_Stack_Up[g][ibin]->SetMaximum(jetshape_y_max);
	JetShape_Stack_Up[g][ibin]->SetMinimum(jetshape_y_min);
      }else{
      	JetShape_Stack_Up[g][ibin]->SetMaximum(jetshape_y_max_normed);
	JetShape_Stack_Up[g][ibin]->SetMinimum(jetshape_y_min_normed);

      }
    }
  
  
  } //ibin
  
  
  cout<<"ready to draw"<<endl;
 
  TCanvas *PAS_plot = new TCanvas("JetShape_ForPAS","",2500,1000);
  PAS_plot->Divide(5,2,0.,0.);

  PAS_plot->cd(6);
  TLegend *legend = new TLegend(0.1,0.05,0.9,.85);
  legend->AddEntry(JetShape_noerr_up[0][0][0],"0.7 < p_{T}^{assoc.}< 1 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][1],"1 < p_{T}^{assoc.}< 2 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][2],"2 < p_{T}^{assoc.}< 3 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][3],"3 < p_{T}^{assoc.}< 4 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][4],"4 < p_{T}^{assoc.}< 8 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][5],"8 < p_{T}^{assoc.}< 12 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][6],"12 < p_{T}^{assoc.}< 16 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][7],"16 < p_{T}^{assoc.}< 20 GeV","f");
 

  if(is_number){
    legend->AddEntry(JetShape_syst[0][0][9],"Total 0.7 < p_{T}^{assoc.}< 20 GeV","lpfe");
  }else{ 
    legend->AddEntry(JetShape_noerr_up[0][0][8],"20 < p_{T}^{assoc.}< 300 GeV","f");
    legend->AddEntry(JetShape_syst[0][0][9],"Total 0.7 < p_{T}^{assoc.} < 300 GeV","lpfe");
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
   
  TLatex   *luminosity_tex_pp = new TLatex(0.2,0.92,"pp 27.4 pb^{-1} (5.02 TeV)");
  luminosity_tex_pp->SetTextFont(43);
  luminosity_tex_pp->SetTextSizePixels(30);
  luminosity_tex_pp->SetLineColor(kWhite);
  luminosity_tex_pp->SetNDC();
  luminosity_tex_pp->Draw();
 
  TLatex   *luminosity_tex_PbPb = new TLatex(0.4,0.92,"PbPb 404 #mub^{-1} (5.02 TeV)");
  luminosity_tex_PbPb->SetTextFont(43);
  luminosity_tex_PbPb->SetTextSizePixels(30);
  luminosity_tex_PbPb->SetLineColor(kWhite);
  luminosity_tex_PbPb->SetNDC();
  luminosity_tex_PbPb->Draw();
 
  TLatex   *jet_reco_tex = new TLatex(0.65,0.92,"ak4CaloJets, p_{T}> 120 GeV, |#eta_{jet}| < 1.6");
  jet_reco_tex->SetTextFont(43);
  jet_reco_tex->SetTextSizePixels(30);
  jet_reco_tex->SetLineColor(kWhite);
  jet_reco_tex->SetNDC();
  jet_reco_tex->Draw();


 
  PAS_plot->cd(1);
 
  JetShape_Stack_Up[1][0]->Draw();
  
  JetShape_Stack_Up[1][0]->GetXaxis()->SetRangeUser(0.,r_max);
    JetShape_Stack_Up[1][0]->GetXaxis()->SetNdivisions(505);

  if(!is_number)  gPad->SetLogy();
  JetShape_Stack_Up[1][0]->GetYaxis()->SetLabelSize(y_label_size);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleSize(y_title_size);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleOffset(y_title_offset);
  if(is_number) JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("Y = #frac{1}{N_{jets}} #frac{dN}{d#Deltar}");
  else JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("     #Rho(#Deltar) (GeV)");
  JetShape_Stack_Up[1][0]->GetYaxis()->CenterTitle();

  JetShape_Stack_Down[1][0]->Draw("same");
  JetShape_Stack_Up[1][0]->Draw("same");
  
 
  JetShape_syst[1][0][9]->Draw("same e2 P");

  JetShape[1][0][9]->Draw("same");

  cout<<"here"<<endl;

  if(do_ref){
    JetShape_ref[1][0][9]->Draw("same e2");
  } 
  cout<<"and here"<<endl;
  labels = new TPaveText(0.28,0.8,0.45,.99,"NDC");
    
  labels->SetName("labels");
  labels->SetFillColor(0);
  labels->SetLineColor(0);
  labels->SetTextAlign(11);
  labels->AddText("pp reference");
  labels->SetTextSize(x_label_size);
  labels->Draw("same");


  TLine *l_dr = new TLine(0.,0.,TMath::Pi()/2.,0.);
  l_dr->SetLineStyle(2);
  l_dr->Draw();


  TLatex *cms_tex = new TLatex(0.3,0.8,"CMS");
  cms_tex->SetTextFont(63);
  cms_tex->SetTextSizePixels(35);
  cms_tex->SetLineColor(kWhite);
  cms_tex->SetNDC();
  cms_tex->Draw(); 


  TLatex *prelim_tex = new TLatex(0.45,0.8,"Preliminary");
  prelim_tex->SetTextFont(53);
  prelim_tex->SetTextSizePixels(35);
  prelim_tex->SetLineColor(kWhite);
  prelim_tex->SetNDC();
  prelim_tex->Draw(); 


  gPad->RedrawAxis();

  for(int ibin = 0; ibin<4; ibin++){
    PAS_plot->cd(5-ibin);


  
    JetShape_Stack_Up[0][ibin]->Draw();

 
    JetShape_Stack_Up[0][ibin]->GetXaxis()->SetRangeUser(0.,r_max);
    JetShape_Stack_Up[0][ibin]->GetXaxis()->SetNdivisions(505);

    if(!is_number) gPad->SetLogy();
    JetShape_Stack_Up[0][ibin]->GetYaxis()->SetLabelSize(0.);
    JetShape_Stack_Down[0][ibin]->Draw("same");

 

    JetShape_syst[0][ibin][9]->Draw("same e2 P");
    JetShape[0][ibin][9]->Draw("same");

  

    TPaveText  *labels = new TPaveText(0.05,0.8,0.45,.99,"NDC");
    
    labels->SetName("labels");
    labels->SetFillColor(0);
    labels->SetLineColor(0);
    labels->SetTextAlign(11);
    labels->AddText((TString)("PbPb "+CBin_labels[ibin]));
    labels->SetTextSize(x_label_size);
    labels->Draw("same");


    if(do_ref){
      JetShape_ref[0][ibin][9]->Draw("same e2 P");
    }

   
    PAS_plot->cd(10-ibin);

    if(is_number){

    
      JetShape_diff[0][ibin][9]->GetYaxis()->SetTitleSize(0.); 
      JetShape_diff[0][ibin][9]->GetYaxis()->SetLabelSize(0.); 
      JetShape_diff[0][ibin][9]->GetXaxis()->SetNdivisions(505);
      JetShape_diff[0][ibin][9]->GetXaxis()->SetTitle("#Deltar");
      JetShape_diff[0][ibin][9]->GetXaxis()->CenterTitle();
      JetShape_diff[0][ibin][9]->GetXaxis()->SetTitleOffset(x_title_offset);
      JetShape_diff[0][ibin][9]->GetXaxis()->SetTitleSize(x_title_size);
      JetShape_diff[0][ibin][9]->GetXaxis()->SetLabelSize(x_label_size);

      JetShape_diff[0][ibin][9]->SetMinimum(number_diff_y_min_r);
      JetShape_diff[0][ibin][9]->SetMaximum(number_diff_y_max_r);
      JetShape_diff[0][ibin][9]->GetXaxis()->SetRangeUser(0.,r_max);
      JetShape_diff[0][ibin][9]->Draw();
      JetShape_diff[0][ibin][9]->GetYaxis()->SetTitle("Y_{PbPb} - Y_{pp}");
      JetShape_Diff_Stack_Up[0][ibin]->Draw("same");
      JetShape_Diff_Stack_Down[0][ibin]->Draw("same");
      JetShape_diff[0][ibin][9]->Draw("same");
      JetShape_syst[2][ibin][9]->Draw("same e2 P");
      JetShape_diff[0][ibin][9]->SetMarkerStyle(24);
      JetShape_diff[0][ibin][9]->Draw("same");

      if(do_ref){
	JetShape_ref[3][ibin][9] = (TH1D*)	JetShape_ref[0][ibin][9]->Clone(Form("JetShapeRefDiff%d",ibin));
	JetShape_ref[3][ibin][9]->Add(JetShape_ref[1][0][9],-1.);
	JetShape_ref[3][ibin][9]->Draw("same e2 P");
      }

    }else{

      JetShape_ratio[0][ibin][9]->SetMarkerStyle(34);
      
      JetShape_ratio[0][ibin][9]->GetYaxis()->SetTitleSize(0.); 
      JetShape_ratio[0][ibin][9]->GetYaxis()->SetLabelSize(0.); 
      JetShape_ratio[0][ibin][9]->GetXaxis()->SetTitle("#Deltar");
      JetShape_ratio[0][ibin][9]->GetXaxis()->CenterTitle();
      JetShape_ratio[0][ibin][9]->GetXaxis()->SetTitleSize(x_title_size);
      JetShape_ratio[0][ibin][9]->SetMinimum(jetshape_ratio_y_min);
      JetShape_ratio[0][ibin][9]->SetMaximum(jetshape_ratio_y_max);
      JetShape_ratio[0][ibin][9]->GetYaxis()->SetTitleOffset(1.3);
      JetShape_ratio[0][ibin][9]->GetXaxis()->SetRangeUser(0.,r_max);
      JetShape_ratio[0][ibin][9]->Draw();
      JetShape_ratio[0][ibin][9]->GetYaxis()->SetTitle("#Rho(r)_{PbPb}/#Rho(r)_{pp}");

      JetShape_syst[2][ibin][9]->Draw("same e2 P");
      JetShape_ratio[0][ibin][9]->Draw("same");
 
      if(!is_number&&do_ref)    JetShape_ref_ratio[0][ibin][9]->Draw("same");
  
    }

    TLatex  *label_ratio = new TLatex(0.05,0.9,"");
    label_ratio->SetTextSize(0.09);
    label_ratio->SetLineColor(kWhite);
    label_ratio->SetNDC();
    label_ratio->Draw();

    if(ibin==3){
      TLegend *legend_ratio = new TLegend(0.02,0.7,0.9,0.9);
      if(is_number)legend_ratio->AddEntry(JetShape_diff[0][ibin][9],"PbPb - pp");
      else legend_ratio->AddEntry(JetShape_ratio[0][ibin][9],"PbPb / pp");
      if(!is_number&&do_ref) legend_ratio->AddEntry(JetShape_ref_ratio[0][ibin][9],"Results at QM","lpfe");
      legend_ratio->SetLineColor(kWhite);
      legend_ratio->SetFillColor(kWhite);
      legend_ratio->SetTextSize(0.07);
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
   
  }


  PAS_plot->cd(6);
  if(is_number){
    TGaxis *dummy_axis_jetshape = new TGaxis(1.,0.18,1.0,.975,number_diff_y_min_r, number_diff_y_max_r);

    dummy_axis_jetshape->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetYaxis());
    dummy_axis_jetshape->SetTitleOffset(1.1);
    dummy_axis_jetshape->SetTickSize(0.);
    dummy_axis_jetshape->CenterTitle();
    dummy_axis_jetshape->SetTitleSize(y_title_size*.8);
    dummy_axis_jetshape->SetLabelSize(y_label_size*.8);
    dummy_axis_jetshape->SetTitle("Y_{PbPb} - Y_{pp}");
    dummy_axis_jetshape->Draw();
  }else{
    TGaxis *dummy_axis_jetshape = new TGaxis(1.,0.18,1.0,.975,jetshape_ratio_y_min,jetshape_ratio_y_max);
    dummy_axis_jetshape->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetYaxis());
    dummy_axis_jetshape->SetTitleOffset(1.1);
    dummy_axis_jetshape->SetTickSize(0.);
    dummy_axis_jetshape->CenterTitle();
    dummy_axis_jetshape->SetTitleSize(y_title_size*.8);
    dummy_axis_jetshape->SetLabelSize(y_label_size*.8);
    dummy_axis_jetshape->SetTitle("#Rho(#Deltar)_{PbPb}/#Rho(#Deltar)_{pp}");
    dummy_axis_jetshape->Draw();
  }
 
  TGaxis *dummy_axis_r = new TGaxis(.235,1.,1.,1.,0.,r_max);

  dummy_axis_r->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetXaxis());
  dummy_axis_r->SetTitleOffset(x_title_offset);
  dummy_axis_r->SetTitleSize(x_title_size);
  dummy_axis_r->SetTickSize(0.);
  dummy_axis_r->SetNdivisions(505);
  dummy_axis_r->Draw();

  PAS_plot->cd(4);
 
  if(do_ref){
    if(!is_number){
        
      PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted_QMRef.pdf");
      PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted_QMRef.png");

    }else{
      PAS_plot->SaveAs("ParticleYield_by_dR_QMRef.pdf");
      PAS_plot->SaveAs("ParticleYield_by_dR_QMRef.png");

    }
  }else{

   if(!is_number){
        
      PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted.pdf");
      PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted.png");

    }else{
      PAS_plot->SaveAs("ParticleYield_by_dR.pdf");
      PAS_plot->SaveAs("ParticleYield_by_dR.png");

    }




  }

  if(is_number){


    TCanvas *c_integral = new TCanvas("Integral_Canvas","",10,10,2000,1000);
    c_integral->Divide(4,2,0,0);

    for(int ibin = 0; ibin < nCBins; ibin++){
      c_integral->cd(4-ibin);


      Integral_Pt[0][ibin] = (TH1D*)f_in->Get((TString)("Integral_PbPb"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]))->Clone((TString)("Integral_PbPb"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));

      Integral_syst_Pt[0][ibin] = (TH1D*)f_in->Get((TString)("Integral_PbPb_Syst"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]))->Clone((TString)("Integral_PbPb"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));

      cout<<"got integrals"<<endl;
    
      if(ibin==0){
	Integral_Pt[1][0] = (TH1D*)f_in->Get((TString)("Integral_pp"));
	Integral_syst_Pt[1][0] = (TH1D*)f_in->Get((TString)("Integral_pp_Syst"));

	Integral_Pt[1][0]->SetMarkerStyle(4);
	Integral_syst_Pt[1][0]->SetMarkerStyle(4);
	Integral_syst_Pt[1][0]->SetFillColor(kGray);

	Integral_Pt[1][0]->SetMarkerSize(2);

	Integral_Pt[1][0]->SetLineColor(kBlack);
	Integral_Pt[1][0]->SetMarkerColor(kBlack);

	Integral_syst_Pt[1][0]->SetMarkerSize(2);

	Integral_syst_Pt[1][0]->SetLineColor(kBlack);
	Integral_syst_Pt[1][0]->SetMarkerColor(kBlack);
  
      }

      Integral_Pt[0][ibin]->SetMarkerSize(2);

      Integral_Pt[0][ibin]->SetLineColor(kBlack);
      Integral_Pt[0][ibin]->SetMarkerColor(kBlack);

      Integral_syst_Pt[0][ibin]->SetMarkerSize(2);

      Integral_syst_Pt[0][ibin]->SetLineColor(kBlack);
      Integral_syst_Pt[0][ibin]->SetMarkerColor(kBlack);
      
      Integral_Pt[0][ibin]->SetMarkerStyle(20);
      Integral_syst_Pt[0][ibin]->SetMarkerStyle(20);
      
      Integral_syst_Pt[0][ibin]->SetFillColor(kCyan-6);
 
      Integral_Pt[0][ibin]->GetXaxis()->SetRangeUser(0.01,19.5);
  
      Integral_Pt[0][ibin]->SetMaximum(integral_max);
      Integral_Pt[0][ibin]->SetMinimum(integral_min);

    
      if(ibin==3){
	Integral_Pt[0][ibin]->GetYaxis()->SetLabelSize(y_title_size);
	Integral_Pt[0][ibin]->GetYaxis()->SetTitleSize(y_title_size);
	Integral_Pt[0][ibin]->GetYaxis()->CenterTitle();
	Integral_Pt[0][ibin]->GetYaxis()->SetTitle("dN/dp_{T}");
      }else{
	Integral_Pt[0][ibin]->GetYaxis()->SetLabelSize(0.);
      }
      Integral_Pt[0][ibin]->GetXaxis()->SetLabelSize(0.);
   
      Integral_Pt[0][ibin]->Draw();
      Integral_syst_Pt[0][ibin]->Draw("same e2");
      Integral_syst_Pt[1][0]->Draw("same e2");
      Integral_Pt[1][0]->Draw("same");
      
      Integral_Pt[0][ibin]->Draw("same");
      Integral_Pt[1][0]->Draw("same");
   
      if(ibin==3){
	TLegend *int_legend = new TLegend(0.28,0.45,0.9,0.7);
	int_legend->AddEntry( Integral_syst_Pt[0][ibin],"PbPb integral |#Deltar| < 1.0");
	int_legend->AddEntry( Integral_syst_Pt[1][0],"pp integral |#Deltar| < 1.0");
	int_legend->SetLineColor(kWhite);
	int_legend->SetTextSize(x_label_size*.9);
	int_legend->Draw();

	labels = new TPaveText(0.28,0.85,0.45,0.95,"NDC");
    
	labels->SetName("labels");
	labels->SetFillColor(0);
	labels->SetLineColor(0);
	labels->SetTextAlign(11);
	labels->AddText(CBin_labels[ibin]);
	labels->SetTextSize(x_label_size);
	labels->Draw("same");


	TLatex *cms_tex = new TLatex(0.28,0.75,"CMS");
	cms_tex->SetTextFont(63);
	cms_tex->SetTextSizePixels(35);
	cms_tex->SetLineColor(kWhite);
	cms_tex->SetNDC();
	cms_tex->Draw(); 


	TLatex *prelim_tex = new TLatex(0.43,0.75,"Preliminary");
	prelim_tex->SetTextFont(53);
	prelim_tex->SetTextSizePixels(35);
	prelim_tex->SetLineColor(kWhite);
	prelim_tex->SetNDC();
	prelim_tex->Draw(); 
  
      }else{
	labels = new TPaveText(0.05,0.85,0.45,0.95,"NDC");
    
	labels->SetName("labels");
	labels->SetFillColor(0);
	labels->SetLineColor(0);
	labels->SetTextAlign(11);
	labels->AddText(CBin_labels[ibin]);
	labels->SetTextSize(x_label_size);
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
      Integral_diff_Pt[2][ibin] = (TH1D*) f_in_ref_int->Get((TString)("Integrated_Yield_Eta_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]))->Clone((TString)("Integrated_Yield_Eta_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));

      Integral_diff_graph_ref[2][ibin] = (TGraphAsymmErrors*) f_in_ref_int->Get((TString)("Integrated_Yield_Eta_Syst_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]))->Clone((TString)("Integral_Error_Up_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));
   
      Integral_diff_graph_ref[2][ibin]->SetFillStyle(3001);
      Integral_diff_graph_ref[2][ibin]->SetFillColor(kYellow-4);

      Integral_diff_Pt[2][ibin]->SetLineColor(kRed);
      Integral_diff_Pt[2][ibin]->SetMarkerColor(kRed);
      Integral_diff_Pt[2][ibin]->SetMarkerStyle(34);
      Integral_diff_Pt[2][ibin]->SetMarkerSize(2);

      Integral_diff_graph_ref[2][ibin]->SetLineColor(kRed);
      Integral_diff_graph_ref[2][ibin]->SetMarkerColor(kRed);
      Integral_diff_graph_ref[2][ibin]->SetMarkerStyle(34);
      Integral_diff_graph_ref[2][ibin]->SetMarkerSize(2);
    
      Integral_diff_Pt[0][ibin]->SetMaximum(integral_diff_max);
      Integral_diff_Pt[0][ibin]->SetMinimum(integral_diff_min);

      if(ibin==3){ 
	Integral_diff_Pt[0][ibin]->GetYaxis()->SetTitleSize(y_title_size*0.8);
	Integral_diff_Pt[0][ibin]->GetYaxis()->SetLabelSize(y_label_size*0.8);
	Integral_diff_Pt[0][ibin]->GetYaxis()->CenterTitle();
	Integral_diff_Pt[0][ibin]->GetXaxis()->SetTitleSize(x_title_size*.8);
	Integral_diff_Pt[0][ibin]->GetXaxis()->SetLabelSize(x_label_size*.8);
	Integral_diff_Pt[0][ibin]->GetXaxis()->CenterTitle();
	Integral_diff_Pt[0][ibin]->GetXaxis()->SetTitle("p_{T}^{assoc}");
	
      }else{
	Integral_diff_Pt[0][ibin]->GetXaxis()->SetTitleSize(x_title_size*.9);
	Integral_diff_Pt[0][ibin]->GetXaxis()->SetLabelSize(x_label_size*.9);
	Integral_diff_Pt[0][ibin]->GetXaxis()->CenterTitle();

      }
      Integral_diff_Pt[0][ibin]->Draw();
      Integral_diff_syst_Pt[0][ibin]->Draw("same e2");
        
      Integral_diff_graph_ref[2][ibin]->Draw("same 2");
       

      Integral_diff_Pt[0][ibin]->Draw("same");
      Integral_diff_Pt[2][ibin]->Draw("same");

      if(ibin==3){
	labels = new TPaveText(0.28,0.85,0.45,0.95,"NDC");
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

      if(ibin==3){
	TLegend *int_legend = new TLegend(0.28,0.6,0.9,0.85);
	int_legend->AddEntry( Integral_diff_syst_Pt[0][ibin],"5.02 TeV");
	int_legend->AddEntry( Integral_diff_graph_ref[2][0],"2.76 TeV");
	int_legend->SetLineColor(kWhite);
	int_legend->SetTextSize(0.06);
	int_legend->Draw();

	labels = new TPaveText(0.28,0.85,0.45,0.95,"NDC");
	labels->SetTextSize(0.055);
      }else{
	labels = new TPaveText(0.05,0.85,0.45,0.95,"NDC");
	labels->SetTextSize(0.06);
      }

    
  
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

   

    }

    c_integral->SaveAs("Integral_dR.png");
    c_integral->SaveAs("Integral_dR.pdf");

    if(is_number){



      for(int ibin = 0; ibin < nCBins; ibin++){

	for(int ibin3 = 0; ibin3<8; ibin3++){




	  signal_dEta_PbPb_pp[0][ibin][ibin3] = (TH1D*)signal_dEta_rebin[0][ibin][ibin3]->Clone((TString)("Diff_dEta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	  signal_dEta_PbPb_pp[0][ibin][ibin3]->Add(signal_dEta_rebin[1][0][ibin3],-1.);


	  signal_dEta_PbPb_pp_syst[0][ibin][ibin3] = (TH1D*)signal_dEta_syst[0][ibin][ibin3]->Clone((TString)("Diff_dEta_syst_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	  signal_dEta_PbPb_pp_syst[0][ibin][ibin3]->Add(signal_dEta_syst[1][0][ibin3],-1.);

	  signal_dPhi_PbPb_pp[0][ibin][ibin3] = (TH1D*)signal_dPhi_rebin[0][ibin][ibin3]->Clone((TString)("Diff_dPhi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	  signal_dPhi_PbPb_pp[0][ibin][ibin3]->Add(signal_dPhi_rebin[1][0][ibin3],-1.);


	  signal_dPhi_PbPb_pp_syst[0][ibin][ibin3] = (TH1D*)signal_dPhi_syst[0][ibin][ibin3]->Clone((TString)("Diff_dPhi_syst_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	  signal_dPhi_PbPb_pp_syst[0][ibin][ibin3]->Add(signal_dPhi_syst[1][0][ibin3],-1.);



  
	  signal_dEta_PbPb_pp_up[0][ibin][ibin3] = (TH1D*)signal_dEta_PbPb_pp[0][ibin][ibin3]->Clone((TString)("Diff_dEta_Up"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	  signal_dEta_PbPb_pp_down[0][ibin][ibin3] = (TH1D*)signal_dEta_PbPb_pp[0][ibin][ibin3]->Clone((TString)("Diff_dEta_Down"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));


	  signal_dPhi_PbPb_pp_up[0][ibin][ibin3] = (TH1D*)signal_dPhi_PbPb_pp[0][ibin][ibin3]->Clone((TString)("Diff_dPhi_Up"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	  signal_dPhi_PbPb_pp_down[0][ibin][ibin3] = (TH1D*)signal_dPhi_PbPb_pp[0][ibin][ibin3]->Clone((TString)("Diff_dPhi_Down"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));


	}


	signal_dEta_rebin[0][ibin][8] = (TH1D*) signal_dEta_rebin[0][ibin][0]->Clone((TString)("Combined_dEta_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
	signal_dPhi_rebin[0][ibin][8] = (TH1D*) signal_dPhi_rebin[0][ibin][0]->Clone((TString)("Combined_dPhi_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
 
	signal_dEta_syst[0][ibin][8] = (TH1D*) signal_dEta_syst[0][ibin][0]->Clone((TString)("Combined_dEta_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
	signal_dPhi_syst[0][ibin][8] = (TH1D*) signal_dPhi_syst[0][ibin][0]->Clone((TString)("Combined_dPhi_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
 
	if(ibin==0){
	  signal_dEta_rebin[1][ibin][8] = (TH1D*) signal_dEta_rebin[1][0][ibin]->Clone((TString)("Combined_dEta_pp_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
	  signal_dPhi_rebin[1][ibin][8] = (TH1D*) signal_dPhi_rebin[1][0][ibin]->Clone((TString)("Combined_dPhi_pp_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
 
	  signal_dEta_syst[1][ibin][8] = (TH1D*) signal_dEta_syst[1][0][ibin]->Clone((TString)("Combined_dEta_pp_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
	  signal_dPhi_syst[1][ibin][8] = (TH1D*) signal_dPhi_syst[1][0][ibin]->Clone((TString)("Combined_dPhi_pp_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
	}


	signal_dEta_PbPb_pp[0][ibin][8] = (TH1D*) signal_dEta_PbPb_pp[0][ibin][0]->Clone((TString)("Combined_dEta_PbPb_pp_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
	signal_dPhi_PbPb_pp[0][ibin][8] = (TH1D*) signal_dPhi_PbPb_pp[0][ibin][0]->Clone((TString)("Combined_dPhi_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
 
	signal_dEta_PbPb_pp_syst[0][ibin][8] = (TH1D*) signal_dEta_PbPb_pp_syst[0][ibin][0]->Clone((TString)("Combined_dEta_PbPb_pp_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
	signal_dPhi_PbPb_pp_syst[0][ibin][8] = (TH1D*) signal_dPhi_PbPb_pp_syst[0][ibin][0]->Clone((TString)("Combined_dPhi_PbPb_pp_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[0] + "_" + TrkPtBin_strs[8]));
 


	for(int k = 1; k<8; k++){
 

	  signal_dEta_rebin[0][ibin][8]->Add(signal_dEta_rebin[0][ibin][k]);
	  signal_dPhi_rebin[0][ibin][8]->Add(signal_dPhi_rebin[0][ibin][k]);

	  signal_dEta_syst[0][ibin][8]->Add(signal_dEta_syst[0][ibin][k]);
	  signal_dPhi_syst[0][ibin][8]->Add(signal_dPhi_syst[0][ibin][k]);


	  signal_dEta_PbPb_pp[0][ibin][8]->Add(signal_dEta_PbPb_pp[0][ibin][k]);
	  signal_dPhi_PbPb_pp[0][ibin][8]->Add(signal_dPhi_PbPb_pp[0][ibin][k]);

	  signal_dEta_PbPb_pp_syst[0][ibin][8]->Add(signal_dEta_PbPb_pp_syst[0][ibin][k]);
	  signal_dPhi_PbPb_pp_syst[0][ibin][8]->Add(signal_dPhi_PbPb_pp_syst[0][ibin][k]);
     

     
	  if(ibin==0){
	    signal_dEta_syst[1][ibin][8]->Add(signal_dEta_syst[1][ibin][k]);
	    signal_dPhi_syst[1][ibin][8]->Add(signal_dPhi_syst[1][ibin][k]);

	    signal_dEta_rebin[1][ibin][8]->Add(signal_dEta_rebin[1][ibin][k]);
	    signal_dPhi_rebin[1][ibin][8]->Add(signal_dPhi_rebin[1][ibin][k]);

	  }


	}
  


	signal_dEta_syst[0][ibin][8]->SetFillColor(kBlack);
	signal_dPhi_syst[0][ibin][8]->SetFillColor(kBlack);

	signal_dEta_syst[0][ibin][8]->SetMarkerColor(kWhite);
	signal_dPhi_syst[0][ibin][8]->SetMarkerColor(kWhite);

    

	signal_dEta_syst[0][ibin][8]->SetFillStyle(3004);
	signal_dPhi_syst[0][ibin][8]->SetFillStyle(3004);

	signal_dEta_rebin[0][ibin][8]->SetMarkerStyle(24);
	signal_dPhi_rebin[0][ibin][8]->SetMarkerStyle(24);



 
	signal_dEta_PbPb_pp_syst[0][ibin][8]->SetFillColor(kBlack);
	signal_dPhi_PbPb_pp_syst[0][ibin][8]->SetFillColor(kBlack);

	signal_dEta_PbPb_pp_syst[0][ibin][8]->SetMarkerColor(kWhite);
	signal_dPhi_PbPb_pp_syst[0][ibin][8]->SetMarkerColor(kWhite);

    

	signal_dEta_PbPb_pp_syst[0][ibin][8]->SetFillStyle(3004);
	signal_dPhi_PbPb_pp_syst[0][ibin][8]->SetFillStyle(3004);

	signal_dEta_PbPb_pp[0][ibin][8]->SetMarkerStyle(24);
	signal_dPhi_PbPb_pp[0][ibin][8]->SetMarkerStyle(24);



     
	if(ibin==0){
	  signal_dEta_syst[1][0][8]->SetFillColor(kBlack);
	  signal_dPhi_syst[1][0][8]->SetFillColor(kBlack);
 
	  signal_dEta_syst[1][0][8]->SetMarkerColor(kWhite);
	  signal_dPhi_syst[1][0][8]->SetMarkerColor(kWhite);

 

	  signal_dEta_syst[1][0][8]->SetFillStyle(3004);
	  signal_dPhi_syst[1][0][8]->SetFillStyle(3004);

	  signal_dEta_rebin[1][0][8]->SetMarkerStyle(24);
	  signal_dPhi_rebin[1][0][8]->SetMarkerStyle(24);


	}




	signal_dEta_stack[0][ibin] = new THStack((TString)("Signal_dEta_Stack_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
	signal_dPhi_stack[0][ibin] = new THStack((TString)("Signal_dPhi_Stack_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");

	signal_dEta_diff_stack_up[0][ibin] = new THStack((TString)("Signal_dEta_Diff_Stack_Up"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
	signal_dPhi_diff_stack_up[0][ibin] = new THStack((TString)("Signal_dPhi_Diff_Stack_Up"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");

	signal_dEta_diff_stack_down[0][ibin] = new THStack((TString)("Signal_dEta_Diff_Stack_Down"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
	signal_dPhi_diff_stack_down[0][ibin] = new THStack((TString)("Signal_dPhi_Diff_Stack_Down"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");


	if(ibin==0){
	  signal_dEta_stack[1][ibin] = new THStack((TString)("Signal_dEta_Stack_pp"),"");
	  signal_dPhi_stack[1][ibin] = new THStack((TString)("Signal_dPhi_Stack_pp"),"");
      
	}
 
	for(int ibin3 = 0; ibin3 < nTrkPtBins - 2; ibin3++){
	  signal_dEta_rebin[0][ibin][ibin3]->SetMarkerSize(0.);
	  signal_dPhi_rebin[0][ibin][ibin3]->SetMarkerSize(0.);
      
	  signal_dEta_rebin[1][0][ibin3]->SetMarkerSize(0.);
	  signal_dPhi_rebin[1][0][ibin3]->SetMarkerSize(0.);
   
	  signal_dEta_PbPb_pp_up[0][ibin][ibin3]->SetMarkerSize(0.);
	  signal_dPhi_PbPb_pp_up[0][ibin][ibin3]->SetMarkerSize(0.);
    
	  switch(ibin3){
	  case 0:
	    signal_dEta_rebin[0][ibin][ibin3]->SetFillColor(kBlue-9);
	    signal_dPhi_rebin[0][ibin][ibin3]->SetFillColor(kBlue-9);
	    signal_dEta_rebin[1][0][ibin3]->SetFillColor(kBlue-9);
	    signal_dPhi_rebin[1][0][ibin3]->SetFillColor(kBlue-9);

	    signal_dEta_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kBlue-9);
	    signal_dPhi_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kBlue-9);
   
	    signal_dEta_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kBlue-9);
	    signal_dPhi_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kBlue-9);
    
	    break;
	  case 1:
	    signal_dEta_rebin[0][ibin][ibin3]->SetFillColor(kYellow-9);
	    signal_dPhi_rebin[0][ibin][ibin3]->SetFillColor(kYellow-9);
	    signal_dEta_rebin[1][0][ibin3]->SetFillColor(kYellow-9);
	    signal_dPhi_rebin[1][0][ibin3]->SetFillColor(kYellow-9);
   
	    signal_dEta_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kYellow-9);
	    signal_dPhi_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kYellow-9);
      
	    signal_dEta_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kYellow-9);
	    signal_dPhi_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kYellow-9);
   
	    break;
	  case 2:
	    signal_dEta_rebin[0][ibin][ibin3]->SetFillColor(kOrange+1);
	    signal_dPhi_rebin[0][ibin][ibin3]->SetFillColor(kOrange+1);
	    signal_dEta_rebin[1][0][ibin3]->SetFillColor(kOrange+1);
	    signal_dPhi_rebin[1][0][ibin3]->SetFillColor(kOrange+1);
   
	    signal_dEta_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kOrange+1);
	    signal_dPhi_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kOrange+1);

	    signal_dEta_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kOrange+1);
	    signal_dPhi_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kOrange+1);

      	
	    break;
	  case 3:
	    signal_dEta_rebin[0][ibin][ibin3]->SetFillColor(kViolet-5);
	    signal_dPhi_rebin[0][ibin][ibin3]->SetFillColor(kViolet-5);
	    signal_dEta_rebin[1][0][ibin3]->SetFillColor(kViolet-5);
	    signal_dPhi_rebin[1][0][ibin3]->SetFillColor(kViolet-5);

	    signal_dEta_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kViolet-5);
	    signal_dPhi_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kViolet-5);

	    signal_dEta_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kViolet-5);
	    signal_dPhi_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kViolet-5);

	    break;
	  case 4:
	    signal_dEta_rebin[0][ibin][ibin3]->SetFillColor(kGreen+3);
	    signal_dPhi_rebin[0][ibin][ibin3]->SetFillColor(kGreen+3);
	    signal_dEta_rebin[1][0][ibin3]->SetFillColor(kGreen+3);
	    signal_dPhi_rebin[1][0][ibin3]->SetFillColor(kGreen+3);

	    signal_dEta_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kGreen+3);
	    signal_dPhi_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kGreen+3);
    
	    signal_dEta_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kGreen+3);
	    signal_dPhi_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kGreen+3);
    
	    break;
  
	  case 5:
	    signal_dEta_rebin[0][ibin][ibin3]->SetFillColor(kRed);
	    signal_dPhi_rebin[0][ibin][ibin3]->SetFillColor(kRed);
	    signal_dEta_rebin[1][0][ibin3]->SetFillColor(kRed);
	    signal_dPhi_rebin[1][0][ibin3]->SetFillColor(kRed);
 
	    signal_dEta_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kRed);
	    signal_dPhi_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kRed);

	    signal_dEta_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kRed);
	    signal_dPhi_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kRed);
    
	    break;
	  case 6:
	    signal_dEta_rebin[0][ibin][ibin3]->SetFillColor(kRed+1);
	    signal_dPhi_rebin[0][ibin][ibin3]->SetFillColor(kRed+1);
	    signal_dEta_rebin[1][0][ibin3]->SetFillColor(kRed+1);
	    signal_dPhi_rebin[1][0][ibin3]->SetFillColor(kRed+1);
   
	    signal_dEta_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kRed+1);
	    signal_dPhi_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kRed+1);

	    signal_dEta_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kRed+1);
	    signal_dPhi_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kRed+1);

	    break;
	
	  case 7:
	    signal_dEta_rebin[0][ibin][ibin3]->SetFillColor(kRed+2);
	    signal_dPhi_rebin[0][ibin][ibin3]->SetFillColor(kRed+2);
	    signal_dEta_rebin[1][0][ibin3]->SetFillColor(kRed+2);
	    signal_dPhi_rebin[1][0][ibin3]->SetFillColor(kRed+2);
   
	    signal_dEta_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kRed+2);
	    signal_dPhi_PbPb_pp_up[0][ibin][ibin3]->SetFillColor(kRed+2);

	    signal_dEta_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kRed+2);
	    signal_dPhi_PbPb_pp_down[0][ibin][ibin3]->SetFillColor(kRed+2);

	    break;

	  default:
	    break;
	  }
	}
    

	for(int ibin3 = 0; ibin3 < nTrkPtBins- 2; ibin3++){

	  for(int k = 0 ; k< signal_dEta_rebin[0][ibin][ibin3]->GetNbinsX()+1; k++){

	    signal_dEta_rebin[0][ibin][nTrkPtBins- 3 - ibin3]->SetBinError(k,0.);


	    signal_dEta_PbPb_pp_up[0][ibin][nTrkPtBins- 3 - ibin3]->SetBinError(k,0.);

	    signal_dEta_PbPb_pp_down[0][ibin][nTrkPtBins- 3 - ibin3]->SetBinError(k,0.);
	
	    if(	signal_dEta_PbPb_pp_up[0][ibin][nTrkPtBins- 3 - ibin3]->GetBinContent(k)<0)	signal_dEta_PbPb_pp_up[0][ibin][nTrkPtBins- 3 - ibin3]->SetBinContent(k,0.);
	    if(	signal_dEta_PbPb_pp_down[0][ibin][nTrkPtBins- 3 - ibin3]->GetBinContent(k)>0)	signal_dEta_PbPb_pp_down[0][ibin][nTrkPtBins- 3 - ibin3]->SetBinContent(k,0.);
	  }
	
	  signal_dEta_stack[0][ibin]->Add(signal_dEta_rebin[0][ibin][nTrkPtBins- 3 - ibin3]);
	  signal_dEta_diff_stack_up[0][ibin]->Add(signal_dEta_PbPb_pp_up[0][ibin][nTrkPtBins- 3 - ibin3]);
	  signal_dEta_diff_stack_down[0][ibin]->Add(signal_dEta_PbPb_pp_down[0][ibin][nTrkPtBins- 3 - ibin3]);

	  for(int k = 0 ; k< signal_dPhi_rebin[0][ibin][ibin3]->GetNbinsX()+1; k++){
	    signal_dPhi_rebin[0][ibin][nTrkPtBins- 3 - ibin3]->SetBinError(k,0.);
	    signal_dPhi_PbPb_pp_up[0][ibin][nTrkPtBins- 3 - ibin3]->SetBinError(k,0.);
	    signal_dPhi_PbPb_pp_down[0][ibin][nTrkPtBins- 3 - ibin3]->SetBinError(k,0.);
	    if(	signal_dPhi_PbPb_pp_up[0][ibin][nTrkPtBins- 3 - ibin3]->GetBinContent(k)<0)	signal_dPhi_PbPb_pp_up[0][ibin][nTrkPtBins- 3 - ibin3]->SetBinContent(k,0.);
	    if(	signal_dPhi_PbPb_pp_down[0][ibin][nTrkPtBins- 3 - ibin3]->GetBinContent(k)>0)	signal_dPhi_PbPb_pp_down[0][ibin][nTrkPtBins- 3 - ibin3]->SetBinContent(k,0.);

	  }

	  signal_dPhi_stack[0][ibin]->Add(signal_dPhi_rebin[0][ibin][nTrkPtBins- 3 - ibin3]);
	  signal_dPhi_diff_stack_up[0][ibin]->Add(signal_dPhi_PbPb_pp_up[0][ibin][nTrkPtBins- 3 - ibin3]);
	  signal_dPhi_diff_stack_down[0][ibin]->Add(signal_dPhi_PbPb_pp_down[0][ibin][nTrkPtBins- 3 - ibin3]);
      
	  if(ibin==0){
	    for(int k = 0 ; k< signal_dEta_rebin[1][ibin][ibin3]->GetNbinsX()+1; k++){
	      signal_dEta_rebin[1][ibin][nTrkPtBins- 3 - ibin3]->SetBinError(k,0.);
	    }

	    signal_dEta_stack[1][ibin]->Add(signal_dEta_rebin[1][ibin][nTrkPtBins- 3 - ibin3]);

	    for(int k = 0 ; k< signal_dPhi_rebin[1][ibin][ibin3]->GetNbinsX()+1; k++){
	      signal_dPhi_rebin[1][ibin][nTrkPtBins- 3 - ibin3]->SetBinError(k,0.);
	    }

	    signal_dPhi_stack[1][ibin]->Add(signal_dPhi_rebin[1][ibin][nTrkPtBins- 3 - ibin3]);

	  }
	}
      }
      cout<<"ready to draw"<<endl;

      TCanvas *c_stacked_eta = new TCanvas("c_stacked_eta","",10,10,2500,1000);
      c_stacked_eta->Divide(5,2,0.,0.);

      c_stacked_eta->cd(1);

   
      signal_dEta_stack[1][0]->SetMinimum(number_y_min);
      signal_dEta_stack[1][0]->SetMaximum(number_y_max);

      signal_dEta_stack[1][0]->Draw();
      signal_dEta_stack[1][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
      signal_dEta_stack[1][0]->GetXaxis()->SetTitleSize(x_title_size);
      signal_dEta_stack[1][0]->GetXaxis()->SetLabelSize(y_title_size);
      signal_dEta_stack[1][0]->GetXaxis()->CenterTitle();
      signal_dEta_stack[1][0]->GetXaxis()->SetTitle("#Delta#eta");
    
      signal_dEta_stack[1][0]->GetYaxis()->SetTitleSize(y_title_size);
      signal_dEta_stack[1][0]->GetYaxis()->CenterTitle();
      signal_dEta_stack[1][0]->GetYaxis()->SetLabelSize(y_label_size);
      signal_dEta_stack[1][0]->GetYaxis()->SetTitleOffset(y_title_offset);
      signal_dEta_stack[1][0]->GetYaxis()->SetTitle("Y = #frac{1}{N_{jets}} #frac{dN}{d#Delta#eta}");

 
      TLatex *cms_tex = new TLatex(0.28,0.78,"CMS");
      cms_tex->SetTextFont(63);
      cms_tex->SetTextSizePixels(35);
      cms_tex->SetLineColor(kWhite);
      cms_tex->SetNDC();
      cms_tex->Draw(); 


      TLatex *prelim_tex = new TLatex(0.43,0.78,"Preliminary");
      prelim_tex->SetTextFont(53);
      prelim_tex->SetTextSizePixels(35);
      prelim_tex->SetLineColor(kWhite);
      prelim_tex->SetNDC();
      prelim_tex->Draw(); 

      labels = new TPaveText(0.28,0.85,0.45,0.95,"NDC");
    
      labels->SetName("labels");
      labels->SetFillColor(0);
      labels->SetLineColor(0);
      labels->SetTextAlign(11);
      labels->AddText("pp reference");
      labels->SetTextSize(x_label_size);
      labels->Draw("same");

      signal_dEta_syst[1][0][8]->Draw("same e2");
      signal_dEta_rebin[1][0][8]->Draw("same");
  
      for(int ibin = 0; ibin<nCBins; ibin++){
	c_stacked_eta->cd(5-ibin);
	signal_dEta_stack[0][ibin]->SetMinimum(number_y_min);
	signal_dEta_stack[0][ibin]->SetMaximum(number_y_max);


	signal_dEta_stack[0][ibin]->Draw();
	signal_dEta_stack[0][ibin]->GetXaxis()->SetRangeUser(-1.5,1.5);   
	signal_dEta_stack[0][ibin]->GetXaxis()->SetTitle("#Delta#eta");
   
	signal_dEta_stack[0][ibin]->GetYaxis()->SetLabelSize(0.0);
  
   
	labels = new TPaveText(0.05,0.85,0.45,0.95,"NDC");
    
	labels->SetName("labels");
	labels->SetFillColor(0);
	labels->SetLineColor(0);
	labels->SetTextAlign(11);
	labels->AddText((TString)("PbPb "+CBin_labels[ibin]));
	labels->SetTextSize(x_label_size);
	labels->Draw("same");

	signal_dEta_syst[0][ibin][8]->Draw("same e2");
	signal_dEta_rebin[0][ibin][8]->Draw("same");

    
	c_stacked_eta->cd(10-ibin);

	signal_dEta_diff_stack_up[0][ibin]->SetMinimum(number_diff_y_min);
	signal_dEta_diff_stack_up[0][ibin]->SetMaximum(number_diff_y_max);

	signal_dEta_diff_stack_up[0][ibin]->Draw();
	signal_dEta_diff_stack_up[0][ibin]->GetXaxis()->SetRangeUser(-1.5,1.5);
	signal_dEta_diff_stack_up[0][ibin]->GetXaxis()->SetTitleSize(x_title_size);

	signal_dEta_diff_stack_up[0][ibin]->GetXaxis()->SetLabelSize(x_label_size);
	signal_dEta_diff_stack_up[0][ibin]->GetXaxis()->CenterTitle();
	signal_dEta_diff_stack_up[0][ibin]->GetXaxis()->SetTitle("#Delta#eta");
    
	signal_dEta_diff_stack_up[0][ibin]->GetYaxis()->SetTitleSize(0.);
	signal_dEta_diff_stack_up[0][ibin]->GetYaxis()->SetLabelSize(0.);
	signal_dEta_diff_stack_up[0][ibin]->GetYaxis()->CenterTitle();
	signal_dEta_diff_stack_up[0][ibin]->GetYaxis()->SetTitle("Y_{PbPb} - Y_{pp}");
	
	signal_dEta_diff_stack_up[0][ibin]->Draw();
	signal_dEta_diff_stack_down[0][ibin]->Draw("same");

	labels = new TPaveText(0.05,0.85,0.45,0.95,"NDC");
    
	labels->SetName("labels");
	labels->SetFillColor(0);
	labels->SetLineColor(0);
	labels->SetTextAlign(11);
	labels->AddText((TString)("PbPb ("+CBin_labels[ibin]+") minus pp"));
	labels->SetTextSize(x_label_size*0.8);
	labels->Draw("same");


	signal_dEta_PbPb_pp_syst[0][ibin][8]->Draw("same e2");
	signal_dEta_PbPb_pp[0][ibin][8]->Draw("same");


      }

      c_stacked_eta->cd(6);
      legend->Draw();

      TGaxis *dummy_axis_diff = new TGaxis(1.,0.18,1.0,.975,number_diff_y_min, number_diff_y_max);

      dummy_axis_diff->ImportAxisAttributes( signal_dEta_stack[1][0]->GetYaxis());
      dummy_axis_diff->SetTitleOffset(y_title_offset);
  
      dummy_axis_diff->CenterTitle();
      dummy_axis_diff->SetTitleSize(y_title_size*.8);
      dummy_axis_diff->SetLabelSize(y_label_size*.8);
      dummy_axis_diff->SetTitle("Y_{PbPb} - Y_{pp}");
      dummy_axis_diff->SetTickSize(0.);
      dummy_axis_diff->Draw();


      TGaxis *dummy_axis_eta = new TGaxis(.235,1.,1.,1., -1.5, 1.5);

      dummy_axis_eta->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetXaxis());
      dummy_axis_eta->SetTitleOffset(x_title_offset);
      dummy_axis_eta->SetTitleSize(x_title_size);
    
      dummy_axis_eta->SetNdivisions(505);
      dummy_axis_eta->SetTickSize(0.);
      dummy_axis_eta->Draw();

  
      c_stacked_eta->cd(0);

      TLatex *type_tex;
      type_tex = new TLatex(0.05,0.92,"Particle Yield by #Delta#eta");
      type_tex->SetTextSize(0.035);
      type_tex->SetLineColor(kWhite);
      type_tex->SetNDC();
      type_tex->Draw();
   
      TLatex   *luminosity_tex_pp = new TLatex(0.2,0.92,"pp 27.4 pb^{-1} (5.02 TeV)");
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
 
      TLatex   *jet_reco_tex = new TLatex(0.6,0.92,"ak4CaloJets, p_{T }> 120 GeV, | #eta_{jet}| < 1.6");
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

      // gPad->SetLeftMargin(0.25);

      signal_dPhi_stack[1][0]->SetMinimum(number_y_min);
      signal_dPhi_stack[1][0]->SetMaximum(number_y_max);

      signal_dPhi_stack[1][0]->Draw();
      signal_dPhi_stack[1][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
      signal_dPhi_stack[1][0]->GetXaxis()->SetTitleSize(x_title_size);
      signal_dPhi_stack[1][0]->GetXaxis()->SetLabelSize(x_title_size);
      signal_dPhi_stack[1][0]->GetXaxis()->CenterTitle();
      signal_dPhi_stack[1][0]->GetXaxis()->SetTitle("#Delta#phi");
    
      signal_dPhi_stack[1][0]->GetYaxis()->SetTitleSize(y_title_size);
      signal_dPhi_stack[1][0]->GetYaxis()->SetLabelSize(y_label_size);
      signal_dPhi_stack[1][0]->GetYaxis()->CenterTitle();
      signal_dPhi_stack[1][0]->GetYaxis()->SetTitleOffset(y_title_offset);
      signal_dPhi_stack[1][0]->GetYaxis()->SetTitle("Y = #frac{1}{N_{jets}} #frac{dN}{d#Delta#phi}");

      cms_tex->Draw(); 

      prelim_tex->Draw(); 


      labels = new TPaveText(0.28,0.85,0.45,0.95,"NDC");
    
      labels->SetName("labels");
      labels->SetFillColor(0);
      labels->SetLineColor(0);
      labels->SetTextAlign(11);
      labels->AddText("pp reference");
      labels->SetTextSize(x_label_size);
      labels->Draw("same");

      signal_dPhi_syst[1][0][8]->Draw("same e2");
      signal_dPhi_rebin[1][0][8]->Draw("same");
  

      for(int ibin = 0; ibin<nCBins; ibin++){
	c_stacked_phi->cd(5-ibin);
	signal_dPhi_stack[0][ibin]->SetMinimum(number_y_min);
	signal_dPhi_stack[0][ibin]->SetMaximum(number_y_max);

  
	signal_dPhi_stack[0][ibin]->Draw();
	signal_dPhi_stack[0][ibin]->GetXaxis()->SetRangeUser(-1.5,1.5); 
	signal_dPhi_stack[0][ibin]->GetXaxis()->SetTitle("#Delta#phi");
   
	signal_dPhi_stack[0][ibin]->GetYaxis()->SetLabelSize(0.0);

	labels = new TPaveText(0.05,0.85,0.45,0.95,"NDC");
    
	labels->SetName("labels");
	labels->SetFillColor(0);
	labels->SetLineColor(0);
	labels->SetTextAlign(11);
	labels->AddText((TString)("PbPb "+CBin_labels[ibin]));
	labels->SetTextSize(x_label_size);
	labels->Draw("same");

	signal_dPhi_syst[0][ibin][8]->Draw("same e2");
	signal_dPhi_rebin[0][ibin][8]->Draw("same");
    
	c_stacked_phi->cd(10-ibin);

	signal_dPhi_diff_stack_up[0][ibin]->SetMinimum(number_diff_y_min);
	signal_dPhi_diff_stack_up[0][ibin]->SetMaximum(number_diff_y_max);

	signal_dPhi_diff_stack_up[0][ibin]->Draw();
	signal_dPhi_diff_stack_up[0][ibin]->GetXaxis()->SetRangeUser(-1.5,1.5);
	signal_dPhi_diff_stack_up[0][ibin]->GetXaxis()->SetTitleSize(x_title_size);
	signal_dPhi_diff_stack_up[0][ibin]->GetXaxis()->SetLabelSize(x_label_size);
	signal_dPhi_diff_stack_up[0][ibin]->GetXaxis()->CenterTitle();
	signal_dPhi_diff_stack_up[0][ibin]->GetXaxis()->SetTitle("#Delta#phi");
    
	signal_dPhi_diff_stack_up[0][ibin]->GetYaxis()->SetTitle("Y_{PbPb} - Y_{pp}");
	signal_dPhi_diff_stack_up[0][ibin]->GetYaxis()->SetLabelSize(0.0);
	signal_dPhi_diff_stack_up[0][ibin]->GetYaxis()->SetTitleSize(0.0);

	signal_dPhi_diff_stack_up[0][ibin]->Draw();
	signal_dPhi_diff_stack_down[0][ibin]->Draw("same");

	labels = new TPaveText(0.05,0.85,0.45,0.95,"NDC");
    
	labels->SetName("labels");
	labels->SetFillColor(0);
	labels->SetLineColor(0);
	labels->SetTextAlign(11);
	labels->AddText((TString)("PbPb ("+CBin_labels[ibin]+") minus pp"));
	labels->SetTextSize(x_label_size*.8);
	labels->Draw("same");


	signal_dPhi_PbPb_pp_syst[0][ibin][8]->Draw("same e2");
	signal_dPhi_PbPb_pp[0][ibin][8]->Draw("same");


  

      }

      c_stacked_phi->cd(6);


      legend->Draw();

      dummy_axis_diff->Draw();

      TGaxis *dummy_axis_phi = new TGaxis(.235,1.,1.,1., -1.5, 1.5);

      dummy_axis_phi->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetXaxis());
      dummy_axis_phi->SetTitleOffset(x_title_offset);
      dummy_axis_phi->SetTitleSize(x_title_size);
 
      dummy_axis_phi->SetNdivisions(505);
      dummy_axis_phi->SetTickSize(0.);
      dummy_axis_phi->Draw();


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


    }
  }else{


  TCanvas *PAS_plot_3row = new TCanvas("JetShape_With3panels","",2500,1500);
  PAS_plot_3row->Divide(5,3,0.,0.);

  PAS_plot_3row->cd(6);


  PAS_plot_3row->cd(0);
  
  TLatex *type_tex;
  if(is_number)type_tex = new TLatex(0.05,0.945,"Particle Yield by #Deltar");
  else  type_tex = new TLatex(0.05,0.945,"Inclusive Jet Shape");
  type_tex->SetTextSize(0.025);
  type_tex->SetLineColor(kWhite);
  type_tex->SetNDC();
  type_tex->Draw();
   
  TLatex   *luminosity_tex_pp = new TLatex(0.2,0.945,"pp 27.4 pb^{-1} (5.02 TeV)");
  luminosity_tex_pp->SetTextFont(43);
  luminosity_tex_pp->SetTextSizePixels(30);
  luminosity_tex_pp->SetLineColor(kWhite);
  luminosity_tex_pp->SetNDC();
  luminosity_tex_pp->Draw();
 
  TLatex   *luminosity_tex_PbPb = new TLatex(0.4,0.945,"PbPb 404 #mub^{-1} (5.02 TeV)");
  luminosity_tex_PbPb->SetTextFont(43);
  luminosity_tex_PbPb->SetTextSizePixels(30);
  luminosity_tex_PbPb->SetLineColor(kWhite);
  luminosity_tex_PbPb->SetNDC();
  luminosity_tex_PbPb->Draw();
 
  TLatex   *jet_reco_tex = new TLatex(0.6,0.945,"ak4CaloJets, p_{T}> 120 GeV, |#eta_{jet}| < 1.6");
  jet_reco_tex->SetTextFont(43);
  jet_reco_tex->SetTextSizePixels(30);
  jet_reco_tex->SetLineColor(kWhite);
  jet_reco_tex->SetNDC();
  jet_reco_tex->Draw();


 
  PAS_plot_3row->cd(1);
 
  JetShape_Stack_Up[1][0]->Draw();
  
  JetShape_Stack_Up[1][0]->GetXaxis()->SetRangeUser(0.,r_max);
    JetShape_Stack_Up[1][0]->GetXaxis()->SetNdivisions(505);

  if(!is_number)  gPad->SetLogy();
  JetShape_Stack_Up[1][0]->GetYaxis()->SetLabelSize(y_label_size);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleSize(y_title_size);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleOffset(y_title_offset);
  if(is_number) JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("Y = #frac{1}{N_{jets}} #frac{dN}{d#Deltar}");
  else JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("   #Rho(#Deltar)  (GeV)");
  JetShape_Stack_Up[1][0]->GetYaxis()->CenterTitle();

  JetShape_Stack_Down[1][0]->Draw("same");
  JetShape_Stack_Up[1][0]->Draw("same");
  
 
  JetShape_syst[1][0][9]->Draw("same e2 P");

  JetShape[1][0][9]->Draw("same");

  if(do_ref){
    JetShape_ref[1][0][9]->Draw("same e2 P");
  } 
 
  labels = new TPaveText(0.28,0.8,0.45,.99,"NDC");
    
  labels->SetName("labels");
  labels->SetFillColor(0);
  labels->SetLineColor(0);
  labels->SetTextAlign(11);
  labels->AddText("pp reference");
  labels->SetTextSize(x_label_size);
  labels->Draw("same");


  TLine *l_dr = new TLine(0.,0.,TMath::Pi()/2.,0.);
  l_dr->SetLineStyle(2);
  l_dr->Draw();


  TLatex *cms_tex = new TLatex(0.3,0.8,"CMS");
  cms_tex->SetTextFont(63);
  cms_tex->SetTextSizePixels(35);
  cms_tex->SetLineColor(kWhite);
  cms_tex->SetNDC();
  cms_tex->Draw(); 


  TLatex *prelim_tex = new TLatex(0.45,0.8,"Preliminary");
  prelim_tex->SetTextFont(53);
  prelim_tex->SetTextSizePixels(35);
  prelim_tex->SetLineColor(kWhite);
  prelim_tex->SetNDC();
  prelim_tex->Draw(); 
  

  gPad->RedrawAxis();

  for(int ibin = 0; ibin<4; ibin++){
    PAS_plot_3row->cd(5-ibin);


  
    JetShape_Stack_Up[0][ibin]->Draw();

 
    JetShape_Stack_Up[0][ibin]->GetXaxis()->SetRangeUser(0.,r_max);
    JetShape_Stack_Up[0][ibin]->GetXaxis()->SetNdivisions(505);

    if(!is_number) gPad->SetLogy();
    JetShape_Stack_Up[0][ibin]->GetYaxis()->SetLabelSize(0.);
    JetShape_Stack_Down[0][ibin]->Draw("same");

 

    JetShape_syst[0][ibin][9]->Draw("same e2 P");
    JetShape[0][ibin][9]->Draw("same");

    if(do_ref){
      JetShape_ref[0][ibin][9]->Draw("same e2");
    }  



    TPaveText  *labels = new TPaveText(0.05,0.8,0.45,.99,"NDC");
    
    labels->SetName("labels");
    labels->SetFillColor(0);
    labels->SetLineColor(0);
    labels->SetTextAlign(11);
    labels->AddText((TString)("PbPb "+CBin_labels[ibin]));
    labels->SetTextSize(x_label_size);
    labels->Draw("same");

  }
  PAS_plot_3row->cd(6);
 
  JetShape_Stack_Up[4][0]->Draw();
  
  JetShape_Stack_Up[4][0]->GetXaxis()->SetRangeUser(0.,r_max);
  JetShape_Stack_Up[4][0]->GetXaxis()->SetNdivisions(505);

  if(!is_number)  gPad->SetLogy();
  JetShape_Stack_Up[4][0]->GetYaxis()->SetLabelSize(y_label_size);
  JetShape_Stack_Up[4][0]->GetYaxis()->SetTitleSize(y_title_size);
  JetShape_Stack_Up[4][0]->GetYaxis()->SetTitleOffset(y_title_offset);
  if(is_number) JetShape_Stack_Up[4][0]->GetYaxis()->SetTitle("Y = #frac{1}{N_{jets}} #frac{dN}{d#Deltar}");
  else JetShape_Stack_Up[4][0]->GetYaxis()->SetTitle("     #rho(#Deltar)");
  JetShape_Stack_Up[4][0]->GetYaxis()->CenterTitle();

  JetShape_Stack_Down[4][0]->Draw("same");
  JetShape_Stack_Up[4][0]->Draw("same");
  
 
  JetShape_syst[4][0][9]->Draw("same e2 P");

  JetShape[4][0][9]->Draw("same");

  if(do_ref){
    JetShape_ref[4][0][9]->Draw("same e2");
    // JetShape_syst_ref[4][0][9]->Draw("same p e2");
  } 
 
  labels = new TPaveText(0.28,0.8,0.45,.99,"NDC");
    
  labels->SetName("labels");
  labels->SetFillColor(0);
  labels->SetLineColor(0);
  labels->SetTextAlign(11);
  labels->AddText("pp reference");
  labels->SetTextSize(x_label_size);
  labels->Draw("same");

  l_dr->Draw();
  

  for(int ibin = 0; ibin<4; ibin++){
    PAS_plot_3row->cd(10-ibin);


  
    JetShape_Stack_Up[3][ibin]->Draw();

 
    JetShape_Stack_Up[3][ibin]->GetXaxis()->SetRangeUser(0.,r_max);
    JetShape_Stack_Up[3][ibin]->GetXaxis()->SetNdivisions(505);

    if(!is_number) gPad->SetLogy();
    JetShape_Stack_Up[3][ibin]->GetYaxis()->SetLabelSize(0.);
    JetShape_Stack_Down[3][ibin]->Draw("same");

 

    JetShape_syst[3][ibin][9]->Draw("same e2 P");
    JetShape[3][ibin][9]->Draw("same");


    JetShape_syst[3][ibin][9]->Draw("same e2 P");
    JetShape[3][ibin][9]->Draw("same");

    if(do_ref){
      JetShape_ref[3][ibin][9]->SetMarkerColor(kSpring);
      JetShape_ref[3][ibin][9]->SetFillColor(kSpring-1);
      JetShape_ref[3][ibin][9]->SetFillStyle(3005);
      //  JetShape_syst_ref[3][ibin][9]->SetMarkerColor(kSpring);
      JetShape_ref[3][ibin][9]->Draw("same e2");
      //   JetShape_syst_ref[3][ibin][9]->SetFillColor(kWhite);
      //JetShape_syst_ref[3][ibin][9]->Draw("same p e2");
    }  



    TPaveText  *labels = new TPaveText(0.05,0.8,0.45,.99,"NDC");
    
    labels->SetName("labels");
    labels->SetFillColor(0);
    labels->SetLineColor(0);
    labels->SetTextAlign(11);
    labels->AddText((TString)("PbPb "+CBin_labels[ibin]));
    labels->SetTextSize(x_label_size);
    labels->Draw("same");

    
   
    PAS_plot_3row->cd(15-ibin);

    JetShape_ratio[3][ibin][9]->SetMarkerStyle(20);
      
    JetShape_ratio[3][ibin][9]->GetYaxis()->SetTitleSize(0.); 
    JetShape_ratio[3][ibin][9]->GetYaxis()->SetLabelSize(0.); 
    JetShape_ratio[3][ibin][9]->GetXaxis()->SetTitle("#Deltar");
    JetShape_ratio[3][ibin][9]->GetXaxis()->CenterTitle();
    JetShape_ratio[3][ibin][9]->GetXaxis()->SetTitleSize(x_title_size);
    JetShape_ratio[3][ibin][9]->SetMinimum(jetshape_ratio_y_min);
    JetShape_ratio[3][ibin][9]->SetMaximum(jetshape_ratio_y_max);
    JetShape_ratio[3][ibin][9]->GetYaxis()->SetTitleOffset(1.3);
    JetShape_ratio[3][ibin][9]->GetXaxis()->SetRangeUser(0.,r_max);
    JetShape_ratio[3][ibin][9]->Draw();
    JetShape_ratio[3][ibin][9]->GetYaxis()->SetTitle("#Rho(r)_{PbPb}/#Rho(r)_{pp}");

    JetShape_syst[5][ibin][9]->Draw("same e2 P");
    JetShape_ratio[3][ibin][9]->Draw("same");
 
    if(do_ref){
      //  JetShape_syst_ref[5][ibin][9]->Draw("same e2 P");
      JetShape_ref_ratio[3][ibin][9]->Draw("same e2");
    }
    JetShape_syst[5][ibin][9]->Draw("same e2 P");
    JetShape_ratio[3][ibin][9]->Draw("same");

    TLatex  *label_ratio = new TLatex(0.05,0.9,"");
    label_ratio->SetTextSize(0.09);
    label_ratio->SetLineColor(kWhite);
    label_ratio->SetNDC();
    label_ratio->Draw();

    if(ibin==3){
      TLegend *legend_ratio = new TLegend(0.02,0.7,0.9,0.9);
      legend_ratio->AddEntry(JetShape_ratio[0][ibin][9],"#rho(r)_{PbPb} / #rho(r)_{pp}");
      //     if(do_ref) legend_ratio->AddEntry(JetShape_syst_ref[5][ibin][9],"QM Results");
      if(do_ref) legend_ratio->AddEntry(JetShape_ref[3][ibin][9],"QM Results");
      legend_ratio->SetLineColor(kWhite);
      legend_ratio->SetFillColor(kWhite);
      legend_ratio->SetTextSize(0.07);
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
  }

  PAS_plot_3row->cd(11);
  legend->Draw();
 
  TGaxis *dummy_axis_jetshape = new TGaxis(1.,0.18,1.0,.975,jetshape_ratio_y_min,jetshape_ratio_y_max);
  dummy_axis_jetshape->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetYaxis());
  dummy_axis_jetshape->SetTitleOffset(1.1);
  dummy_axis_jetshape->SetTickSize(0.);
  dummy_axis_jetshape->CenterTitle();
  dummy_axis_jetshape->SetTitleSize(y_title_size*.8);
  dummy_axis_jetshape->SetLabelSize(y_label_size*.8);
  dummy_axis_jetshape->SetTitle("#rho(#Deltar)_{PbPb}/#rho(#Deltar)_{pp}");
  dummy_axis_jetshape->Draw();

 
  TGaxis *dummy_axis_r = new TGaxis(.235,1.,1.,1.,0.,r_max);

  dummy_axis_r->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetXaxis());
  dummy_axis_r->SetTitleOffset(x_title_offset);
  dummy_axis_r->SetTitleSize(x_title_size);
  dummy_axis_r->SetTickSize(0.);
  dummy_axis_r->SetNdivisions(505);
  dummy_axis_r->Draw();

 
 
  if(do_ref){
          
    PAS_plot_3row->SaveAs("JetShapes_WithHighpT_pTweighted_3Row_QMRef.pdf");
    PAS_plot_3row->SaveAs("JetShapes_WithHighpT_pTweighted_3Row_QMRef.png");

  }else{

        
      PAS_plot_3row->SaveAs("JetShapes_WithHighpT_pTweighted_3Row.pdf");
      PAS_plot_3row->SaveAs("JetShapes_WithHighpT_pTweighted_3Row.png");
  
  }

  }
  
  return 0;
}// main loop
