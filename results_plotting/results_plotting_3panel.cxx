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

Int_t results_plotting_3panel(bool is_number=0,bool do_ref=kFALSE){

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
  TGraphAsymmErrors* JetShape_graph[6][nCBins][nTrkPtBins];

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
  
  TString stem, datalabel,me00_range_string,stem_mc;
 
  float norm;
 
  TH1D *Integral_Pt[12][5];
  TH1D *Integral_syst_Pt[12][5];
  TH1D *Integral_diff_Pt[12][5];
  TH1D *Integral_diff_syst_Pt[12][5];

  double integral[12][nCBins][nTrkPtBins];
  double integral_err[12][nCBins][nTrkPtBins];
  double integral_syst_err[12][nCBins][nTrkPtBins];

  TPaveText *labels;

  //-----------------------------------
  // Pick input and reference files here
  //----------------------------------
  
 
  TFile *f_in, *f_in_ref, *f_in_eta_phi;

  if(is_number){
    f_in = new TFile("../jet_shapes_result/Jet_Shapes.root");
    f_in_eta_phi = new TFile("../particle_yields/Particle_Yields.root");
    f_in_ref = new TFile("Jet_Shapes_Preapproval.root");

  }else{

    f_in = new TFile("../jet_shapes_result/Jet_Shapes_pTweighted.root");
    f_in_ref = new TFile("Jet_Shapes_pTweighted_Preapproval.root");

  }

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
  float y_title_offset = 1.1;
  float y_tick_length = 0.025;


  float r_max = 0.999;

  float number_y_min = -5.;
  float number_y_max = 40.;
  
  float jetshape_y_min_normed = .005;
  float jetshape_y_max_normed = 23.5;
  
 
  float jetshape_y_min = .5;
  float jetshape_y_max = 2000.;
 
  float jetshape_ratio_y_min = 0.;
  float jetshape_ratio_y_max = 3.2;

  float number_diff_y_min = -3.5;
  float number_diff_y_max = 13.;

  float number_diff_y_min_r = -5.5;
  float number_diff_y_max_r = 15.5;

  float integral_max = 13.;
  float integral_min = -1.;

  float integral_diff_max = 6.2;
  float integral_diff_min = -1.;


 
  //***********************************
  //***********************************

  //-----------------------
  // Start getting histos
  //-----------------------

    
  if(is_number) stem = "JetShape2_Yield_BkgSub_Inclusive_";
  else stem = "JetShape2_Yield_BkgSub_pTweightedInclusive_";


  for(int g=0; g<2; g++){
  
    for (int ibin=0;ibin<nCBins;ibin++){

      if(g==1&&ibin > 0) continue;
    
      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

	if(g==0){
	  JetShape[g][ibin][ibin3] = (TH1D*)f_in->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + "_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	 
	  if(do_ref)    JetShape_ref[g][ibin][ibin3] = (TH1D*)f_in_ref->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + "_Ref_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  JetShape_syst[g][ibin][ibin3] = (TH1D*)f_in->Get((TString)("Jet_Shape_SystErr_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Jet_Shape_SystErr_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	 
	  if(is_number&&ibin3<8){

	    signal_dPhi_rebin[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("Proj_dPhi_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]+"_Rebin"))->Clone((TString) ("Proj_dPhi_PbPb_"  + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 
	      
	    signal_dEta_rebin[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("Proj_dEta_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]+"_Rebin"))->Clone((TString) ("Proj_dEta_PbPb_"  + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 

	    signal_dPhi_syst[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("dPhi_Syst_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("dPhi_Syst_PbPb_"  + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 
	    signal_dEta_syst[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("dEta_Syst_PbPb_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("dEta_Syst_PbPb_"  + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 


	  }

	}else{
	  JetShape[g][ibin][ibin3] = (TH1D*)f_in->Get((TString)(stem +"pp_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + "_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	 
	  if(do_ref)   JetShape_ref[g][ibin][ibin3] = (TH1D*)f_in_ref->Get((TString)(stem +"pp_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + "_Ref_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  JetShape_syst[g][ibin][ibin3] = (TH1D*)f_in->Get((TString)("Jet_Shape_SystErr_pp_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Jet_Shape_SystErr_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	


	  if(is_number&&ibin3<8){
	    signal_dPhi_rebin[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("Proj_dPhi_pp_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]+"_Rebin"))->Clone((TString) ("Proj_dPhi_pp_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 
	    signal_dEta_rebin[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("Proj_dEta_pp_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]+"_Rebin"))->Clone((TString) ("Proj_dEta_pp_"  +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 

	    signal_dPhi_syst[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("dPhi_Syst_pp_Cent0_Cent10_Pt100_Pt300_"+ TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("dPhi_Syst_pp_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 

	    signal_dEta_syst[g][ibin][ibin3] = (TH1D*)f_in_eta_phi->Get((TString)("dEta_Syst_pp_Cent0_Cent10_Pt100_Pt300_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("dEta_Syst_pp_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])); 



	  }
	}
	   

	JetShape_syst[g][ibin][ibin3]->SetMarkerStyle(20);
	JetShape_syst[g][ibin][ibin3]->SetMarkerSize(1);
	JetShape_syst[g][ibin][ibin3]->SetMarkerColor(kWhite);
	JetShape_syst[g][ibin][ibin3]->SetFillColor(kBlack);
	JetShape_syst[g][ibin][ibin3]->SetFillStyle(3004);
	
	
	JetShape[g][ibin][ibin3]->GetXaxis()->SetLabelFont(label_font);
	JetShape[g][ibin][ibin3]->GetXaxis()->SetLabelOffset(x_label_offset);
	JetShape[g][ibin][ibin3]->GetXaxis()->SetLabelSize(x_label_size);
	JetShape[g][ibin][ibin3]->GetXaxis()->SetTitleSize(x_label_size);
	JetShape[g][ibin][ibin3]->GetXaxis()->SetTickLength(x_tick_length);
	JetShape[g][ibin][ibin3]->GetXaxis()->SetTitleOffset(x_title_offset);
	JetShape[g][ibin][ibin3]->GetXaxis()->SetTitleFont(label_font);
	JetShape[g][ibin][ibin3]->GetYaxis()->SetTitle("#Rho(#Deltar)");
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
	if(do_ref){
	  JetShape_ref[g][ibin][ibin3]->SetMarkerSize(1);
	  JetShape_ref[g][ibin][ibin3]->SetLineColor(kCyan);
	  JetShape_ref[g][ibin][ibin3]->SetMarkerColor(kCyan);
	  JetShape_ref[g][ibin][ibin3]->SetMarkerStyle(20);
	}

		
      }

    }
  }

  cout<<"got all histograms"<<endl;


  for (int ibin=0;ibin<nCBins;ibin++){ 
    
    for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

      if(is_number){
	
	JetShape_syst[2][ibin][ibin3] = (TH1D*) JetShape_syst[0][ibin][ibin3]->Clone((TString) ("Jet_Shape_Syst_Diff_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	JetShape_syst[2][ibin][ibin3]->Add(JetShape_syst[1][0][ibin3],-1.);
      }else{
	JetShape_syst[2][ibin][ibin3] = (TH1D*) JetShape_syst[0][ibin][ibin3]->Clone((TString) ("Jet_Shape_Syst_Ratio_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	JetShape_syst[2][ibin][ibin3]->Divide(JetShape_syst[1][0][ibin3]);
      }

      JetShape_syst[2][ibin][ibin3]->SetMarkerStyle(20);
      JetShape_syst[2][ibin][ibin3]->SetMarkerSize(1);
      JetShape_syst[2][ibin][ibin3]->SetMarkerColor(kWhite);
      JetShape_syst[2][ibin][ibin3]->SetFillColor(kBlack);
      JetShape_syst[2][ibin][ibin3]->SetFillStyle(3004);
	

      TString diff_name ="JetShape_diff_"; diff_name+=CBin_strs[ibin]; diff_name+= "_"; diff_name += CBin_strs[ibin+1]; diff_name+= ibin3;
      JetShape_diff[0][ibin][ibin3] = (TH1D*)JetShape[0][ibin][ibin3]->Clone(diff_name);
      JetShape_diff[0][ibin][ibin3]->Add( JetShape[1][0][ibin3],-1. );

      if(do_ref){    
	JetShape_ref_diff[0][ibin][ibin3] = (TH1D*)JetShape_ref[0][ibin][ibin3]->Clone((TString)(diff_name+"_Ref"));
	JetShape_ref_diff[0][ibin][ibin3]->Add( JetShape_ref[1][0][ibin3],-1. );
      }
    
      TString ratio_name ="JetShape_ratio_"; ratio_name+=CBin_strs[ibin]; ratio_name+= "_"; ratio_name += CBin_strs[ibin+1]; ratio_name+= ibin3;

      
      JetShape_ratio[0][ibin][ibin3] = (TH1D*)JetShape[0][ibin][ibin3]->Clone(ratio_name);
      JetShape_ratio[0][ibin][ibin3]->Divide(JetShape[1][0][ibin3]);
      JetShape_ratio[0][ibin][ibin3]->GetXaxis()->SetNdivisions(505);
  
      if(do_ref){
	JetShape_ref_ratio[0][ibin][ibin3] = (TH1D*)JetShape_ref[0][ibin][ibin3]->Clone((TString)(ratio_name+"_Ref"));
	JetShape_ref_ratio[0][ibin][ibin3]->Divide( JetShape_ref[1][0][ibin3]);
      }
    }
  }
         
 
  //------------------------
  // STACK PLOTS
  //------------------------


  for(int ibin = 0; ibin<nCBins; ibin++){
    
    for(int k = 0; k<nTrkPtBins-1; k++){

      
 
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

      JetShape_Stack_Up[0][ibin]->SetMaximum(jetshape_y_max);
      JetShape_Stack_Up[1][0]->SetMaximum(jetshape_y_max);
      JetShape_Stack_Up[0][ibin]->SetMinimum(jetshape_y_min);
      JetShape_Stack_Up[1][0]->SetMinimum(jetshape_y_min);
  
    }else{

      JetShape_Stack_Up[0][ibin]->SetMaximum(number_y_max);
      JetShape_Stack_Up[1][0]->SetMaximum(number_y_max);
      JetShape_Stack_Up[0][ibin]->SetMinimum(number_y_min);
      JetShape_Stack_Up[1][0]->SetMinimum(number_y_min);

  
    }
  
  } //ibin
  

  cout<<"ready to draw"<<endl;
 
  TCanvas *PAS_plot = new TCanvas("JetShape_ForPAS","",2500,1500);
  PAS_plot->Divide(5,3,0.,0.);

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
 
  TLatex   *jet_reco_tex = new TLatex(0.6,0.92,"ak4CaloJets, p_{T}> 120 GeV, |#eta_{jet}| < 1.6");
  jet_reco_tex->SetTextFont(43);
  jet_reco_tex->SetTextSizePixels(35);
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
  else JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("     #Rho(#Deltar)");
  JetShape_Stack_Up[1][0]->GetYaxis()->CenterTitle();

  JetShape_Stack_Down[1][0]->Draw("same");
  JetShape_Stack_Up[1][0]->Draw("same");
  
 
  JetShape_syst[1][0][9]->Draw("same e2 P");

  JetShape[1][0][9]->Draw("same");

  if(do_ref){
    JetShape_ref[1][0][9]->Draw("same");
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
    PAS_plot->cd(5-ibin);


  
    JetShape_Stack_Up[0][ibin]->Draw();

 
    JetShape_Stack_Up[0][ibin]->GetXaxis()->SetRangeUser(0.,r_max);
    JetShape_Stack_Up[0][ibin]->GetXaxis()->SetNdivisions(505);

    if(!is_number) gPad->SetLogy();
    JetShape_Stack_Up[0][ibin]->GetYaxis()->SetLabelSize(0.);
    JetShape_Stack_Down[0][ibin]->Draw("same");

 

    JetShape_syst[0][ibin][9]->Draw("same e2 P");
    JetShape[0][ibin][9]->Draw("same");

    if(do_ref){
      JetShape_ref[0][ibin][9]->SetMarkerColor(kCyan);
      JetShape_ref[0][ibin][9]->SetMarkerStyle(20);
      JetShape_ref[0][ibin][9]->SetMarkerSize(1);
      JetShape_ref[0][ibin][9]->SetLineColor(kCyan);
      JetShape_ref[0][ibin][9]->Draw("same");
    }  

    TPaveText  *labels = new TPaveText(0.05,0.8,0.45,.99,"NDC");
    
    labels->SetName("labels");
    labels->SetFillColor(0);
    labels->SetLineColor(0);
    labels->SetTextAlign(11);
    labels->AddText((TString)("PbPb "+CBin_labels[ibin]));
    labels->SetTextSize(x_label_size);
    labels->Draw("same");


    cout<<"start new norm code"<<endl;

    for (int ibin=0;ibin<nCBins;ibin++){ 
    
      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){


	if(!is_number&&ibin3==9){
	  for(int g=0; g<2; g++){

	    norm = JetShape[g][ibin][ibin3]->Integral("width");
	    for(int k = 0; k<10; k++){
	      JetShape[g][ibin][k]->Scale(1./norm);
	      JetShape_syst[g][ibin][k]->Scale(1./norm);
	    }
	    if(do_ref){
	      norm = JetShape_ref[g][ibin][ibin3]->Integral("width");
	      for(int k = 0; k<10; k++){
		JetShape_ref[g][ibin][k]->Scale(1./norm);
	      }
	    }
	  }
	}

	if(is_number){
	
	  JetShape_syst[2][ibin][ibin3] = (TH1D*) JetShape_syst[0][ibin][ibin3]->Clone((TString) ("Jet_Shape_Syst_Diff_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  JetShape_syst[2][ibin][ibin3]->Add(JetShape_syst[1][0][ibin3],-1.);
	}else{
	  JetShape_syst[2][ibin][ibin3] = (TH1D*) JetShape_syst[0][ibin][ibin3]->Clone((TString) ("Jet_Shape_Syst_Ratio_" + CBin_strs[ibin+1] + "_" +TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  JetShape_syst[2][ibin][ibin3]->Divide(JetShape_syst[1][0][ibin3]);
	}

	JetShape_syst[2][ibin][ibin3]->SetMarkerStyle(20);
	JetShape_syst[2][ibin][ibin3]->SetMarkerSize(1);
	JetShape_syst[2][ibin][ibin3]->SetMarkerColor(kWhite);
	JetShape_syst[2][ibin][ibin3]->SetFillColor(kBlack);
	JetShape_syst[2][ibin][ibin3]->SetFillStyle(3004);

	TString ratio_name ="JetShape_ratio_"; ratio_name+=CBin_strs[ibin]; ratio_name+= "_"; ratio_name += CBin_strs[ibin+1]; ratio_name+= ibin3;

      
	JetShape_ratio[0][ibin][ibin3] = (TH1D*)JetShape[0][ibin][ibin3]->Clone(ratio_name);
	JetShape_ratio[0][ibin][ibin3]->Divide(JetShape[1][0][ibin3]);
	JetShape_ratio[0][ibin][ibin3]->GetXaxis()->SetNdivisions(505);
  
	if(do_ref){
	  JetShape_ref_ratio[0][ibin][ibin3] = (TH1D*)JetShape_ref[0][ibin][ibin3]->Clone((TString)(ratio_name+"_Ref"));
	  JetShape_ref_ratio[0][ibin][ibin3]->Divide( JetShape_ref[1][0][ibin3]);
	}
      }
      
    }
         
    cout<<"here"<<endl;
    //------------------------
    // STACK PLOTS
    //------------------------


    for(int ibin = 0; ibin<nCBins; ibin++){
    
      for(int k = 0; k<nTrkPtBins-1; k++){

      
 
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

	JetShape_Stack_Up[0][ibin]->SetMaximum(jetshape_y_max);
	JetShape_Stack_Up[1][0]->SetMaximum(jetshape_y_max);
	JetShape_Stack_Up[0][ibin]->SetMinimum(jetshape_y_min);
	JetShape_Stack_Up[1][0]->SetMinimum(jetshape_y_min);
  
      }else{

	JetShape_Stack_Up[0][ibin]->SetMaximum(number_y_max);
	JetShape_Stack_Up[1][0]->SetMaximum(number_y_max);
	JetShape_Stack_Up[0][ibin]->SetMinimum(number_y_min);
	JetShape_Stack_Up[1][0]->SetMinimum(number_y_min);

  
      }
  
    } //ibin

    cout<<"drawing ratios"<<endl;
  
   
    PAS_plot->cd(15-ibin);

    JetShape_ratio[0][ibin][9]->SetMarkerStyle(20);
      
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
 
    if(do_ref)    JetShape_ref_ratio[0][ibin][9]->Draw("same");
   

    TLatex  *label_ratio = new TLatex(0.05,0.9,"");
    label_ratio->SetTextSize(0.09);
    label_ratio->SetLineColor(kWhite);
    label_ratio->SetNDC();
    label_ratio->Draw();

    if(ibin==3){
      TLegend *legend_ratio = new TLegend(0.02,0.7,0.9,0.9);
      if(is_number)legend_ratio->AddEntry(JetShape_diff[0][ibin][9],"PbPb - pp");
      else legend_ratio->AddEntry(JetShape_ratio[0][ibin][9],"PbPb / pp");
      if(do_ref) legend_ratio->AddEntry(JetShape_ref_ratio[0][ibin][9],"At Preapproval");
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
   
  }


  PAS_plot->cd(6);
  if(is_number){
    TGaxis *dummy_axis_jetshape = new TGaxis(1.,0.18,1.0,.975,number_diff_y_min_r, number_diff_y_max_r);

    dummy_axis_jetshape->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetYaxis());
    dummy_axis_jetshape->SetTitleOffset(1.1);
    dummy_axis_jetshape->SetTickSize(0.);
    dummy_axis_jetshape->CenterTitle();
    dummy_axis_jetshape->SetTitleSize(y_title_size*.8);
    dummy_axis_jetshape->SetLabelSize(y_label_size);
    dummy_axis_jetshape->SetTitle("Y_{PbPb} - Y_{pp}");
    dummy_axis_jetshape->Draw();
  }else{
    TGaxis *dummy_axis_jetshape = new TGaxis(1.,0.18,1.0,.975,jetshape_ratio_y_min,jetshape_ratio_y_max);
    dummy_axis_jetshape->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetYaxis());
    dummy_axis_jetshape->SetTitleOffset(1.1);
    dummy_axis_jetshape->SetTickSize(0.);
    dummy_axis_jetshape->CenterTitle();
    dummy_axis_jetshape->SetTitleSize(y_title_size);
    dummy_axis_jetshape->SetLabelSize(y_label_size);
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
          
    PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted_PreapprovalRef.pdf");
    PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted_PreapprovalRef.png");
  }else{
        
    PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted.pdf");
    PAS_plot->SaveAs("JetShapes_WithHighpT_pTweighted.png");
  }

  return 0;
}// main loop
