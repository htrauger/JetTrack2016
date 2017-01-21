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

Int_t jet_shapes_run1_run2(){


  //THIS VERSION DOES DATA VS. DATA
  //-------------------------
  // Display scale parameters
  //--------------------------

  float rho_max =3000.;
  float rho_min = .5;

  float ratio_min = 0.;
  float ratio_max = 6.5;


  float diff_min = -70.;
  float diff_max = 90.;



  TFile *fMC[9], *f_ref[7];

  TCanvas *c_jetshape[5];
   
  TString jetetacut, etalabel,centlabel,pTlabel,Ajlabel;
  float eta_ymax;

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins, limR;
  float deta, dphi, r, bc, bg_err,temp1,temp2, rbin, temperr, err, width_temp_x, width_temp_y, width_temp, norm_temp, zerobin, temp, norm_tot, err_temp, cont;
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 10;

  float RBins[20] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.,1.2,1.4,1.6,1.8,2.0};
  //float RBins[20] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.,1.2,1.4,1.6,2.0};

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  
  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins+1] = {0.5, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300","" };
  TString TrkPtBin_strs2[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" ,""};
  TString TrkPtBin_labels[nTrkPtBins] = {"0.5<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.25);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
 
  TH1D* JetShape[6][nCBins][nTrkPtBins];
  TH1D* JetShape2[12][nCBins][nTrkPtBins];

  TH1D* JetShape_ratio[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff[6][nCBins][nTrkPtBins];

  double temp_cont, nextr, nextl, cos_weight,me00_range, mc_error;
  
  TString stem, datalabel,me00_range_string,stem_mc;

 
 
  float norm, err1;


  TFile *fin, *fin2;

  TFile *fin_shape_run1 =new TFile ("../../JetTrack2015/dR_studies/Jet_Shapes.root");
  TFile *fin_shape_run1_inc =new TFile ("../Run1_Studies/jet_shapes_result/Jet_Shapes.root");
  TFile *fin_shape_run1_inc_gen =new TFile ("../Run1_Studies/jet_shapes_closures/Jet_Shapes_ClosuresGenJet_GenTrack.root");
  TFile *fin_shape_run2 = new TFile ("../jet_shapes_result/Jet_Shapes.root");
  TFile *fin_shape_run2_lead = new TFile ("../Leading_Studies/jet_shapes_results/Jet_Shapes.root");
  TFile *fin_shape_gen = new TFile ("../jet_shapes_closures/Jet_Shapes_Closures.root");
 
  TFile *fin_shape_gen_run1 = new TFile ("../jet_shapes_closures/Jet_Shapes_Closures.root");

 TFile *fin_shape_gen_lead = new TFile ("../Leading_Studies/jet_shapes_closures/Jet_Shapes_ClosuresGenJet_GenTrack.root");

  //-----------------------
  // Start getting histos
  //-----------------------
  for(int g=0; g<4; g++){

    if(g==1||g==2)continue;

    cout<<"here!"<<endl;

    stem = "Yield_BkgSub_pTweighted";

    for (int ibin=0;ibin<nCBins;ibin++){
  
     

      if(g==0){
    
	 JetShape2[g][ibin][0] = (TH1D*)fin_shape_run2->Get((TString)("JetShape2_Yield_BkgSub_pTweightedInclusive_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_TrkPt300_"))->Clone((TString)("JetShape2_Yield_BkgSub_pTweightedInclusive_PbPb_TrkPt300"));
	  
	 JetShape2[g+1][ibin][0] = (TH1D*)fin_shape_run1->Get((TString)("JetShape2_Yield_BkgSub_PbPb_Leading_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_TrkPt300_"))->Clone((TString)("JetShape2_Run1_Yield_BkgSub_PbPb_TrkPt300_"));

	 JetShape2[g+2][ibin][0] = (TH1D*)fin_shape_run1_inc->Get((TString)("JetShape2_Yield_BkgSub_pTweightedInclusive_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_TrkPt300_"))->Clone((TString)("JetShape2_Run1_Yield_BkgSub_PbPb_Inclusive_TrkPt300_"));
	 

	 }else{


	if(ibin==0){
	  JetShape2[g][ibin][0] = (TH1D*)fin_shape_run2->Get((TString)("JetShape2_Yield_BkgSub_pTweightedInclusive_pp_TrkPt300_"))->Clone((TString)("JetShape2_Yield_BkgSub_pTweightedInclusive_PbPb_TrkPt300")); //black
	  
	  JetShape2[g+1][ibin][0] = (TH1D*)fin_shape_run1->Get((TString)("JetShape2_Yield_BkgSub_pp_Leading_Cent0_Cent10_TrkPt300_"))->Clone((TString)("JetShape2_Run1_Yield_BkgSub_pp_TrkPt300_")); //blue

	  JetShape2[g+2][ibin][0] = (TH1D*)fin_shape_run1_inc->Get((TString)("JetShape2_Yield_BkgSub_pTweightedInclusive_pp_TrkPt300_"))->Clone((TString)("JetShape2_Run1_Yield_BkgSub_pp_Inclusive_TrkPt300_")); //red


	  JetShape2[g+3][ibin][0] = (TH1D*)fin_shape_gen->Get((TString)("JetShape2_Yield_BkgSub_pTweightedInclusive_Cent0_Cent10_TrkPt300_"))->Clone((TString)("JetShape2_GenYield_BkgSub_pTweightedInclusive_Cent0_Cent10_TrkPt300_")); //orange

	  JetShape2[g+4][ibin][0] = (TH1D*)fin_shape_gen_run1->Get((TString)("JetShape2_Yield_BkgSub_pTweightedInclusive2_Cent0_Cent10_TrkPt300_"))->Clone((TString)("JetShape2_GenYield_BkgSub_pTweightedLeading_Cent0_Cent10_TrkPt300_")); //green

	  JetShape2[g+5][ibin][0] = (TH1D*)fin_shape_run1_inc_gen->Get((TString)("JetShape2_Yield_BkgSub_pTweightedInclusive_pp_TrkPt300_"))->Clone((TString)("JetShape2_GenYield_Run1_Inclusive_BkgSub_pTweightedInclusive2_Cent0_Cent10_TrkPt300_")); //violet

	  JetShape2[g+6][ibin][0] = (TH1D*)fin_shape_gen_lead->Get((TString)("JetShape2_Yield_BkgSub_pTweightedInclusive_pp_TrkPt300_"))->Clone((TString)("JetShape2_GenYield_Run2_Leading_BkgSub_pTweightedLeading2_Cent0_Cent10_TrkPt300_")); //cyan

	  JetShape2[g+7][ibin][0] = (TH1D*)fin_shape_run2_lead->Get((TString)("JetShape2_Yield_BkgSub_pTweightedInclusive_pp_TrkPt300_"))->Clone((TString)("JetShape2_pp_Yield_Run2_Leading_BkgSub_pTweightedLeading2_Cent0_Cent10_TrkPt300_")); //dark green


	}
      }
     

      if(g==3){
	JetShape_ratio[g][ibin][0] = (TH1D*)JetShape2[g-3][ibin][0]->Clone((TString)("JetShape2_Ratio_Inclusive_Cent0_Cent10_TrkPt300_"));
	JetShape_ratio[g+1][ibin][0] = (TH1D*)JetShape2[g-2][ibin][0]->Clone((TString)("JetShape2_Ratio_Leading_Cent0_Cent10_TrkPt300_"));
	JetShape_ratio[g+2][ibin][0] = (TH1D*)JetShape2[g-1][ibin][0]->Clone((TString)("JetShape2_Ratio_Inclusive_Cent0_Cent10_TrkPt300_"));

	JetShape_ratio[g][ibin][0]->Divide(JetShape2[g][0][0]);
	JetShape_ratio[g+1][ibin][0]->Divide(JetShape2[g+1][0][0]);
	JetShape_ratio[g+2][ibin][0]->Divide(JetShape2[g+2][0][0]);


	JetShape_diff[g][ibin][0] = (TH1D*)JetShape2[g-3][ibin][0]->Clone((TString)("JetShape2_Diff_Inclusive_Cent0_Cent10_TrkPt300_"));
	JetShape_diff[g+1][ibin][0] = (TH1D*)JetShape2[g-2][ibin][0]->Clone((TString)("JetShape2_Diff_Leading_Cent0_Cent10_TrkPt300_"));
	JetShape_diff[g+2][ibin][0] = (TH1D*)JetShape2[g-1][ibin][0]->Clone((TString)("JetShape2_Diff_Inclusive_Cent0_Cent10_TrkPt300_"));

	JetShape_diff[g][ibin][0]->Add(JetShape2[g][0][0],-1.);
	JetShape_diff[g+1][ibin][0]->Add(JetShape2[g+1][0][0],-1.);
	JetShape_diff[g+2][ibin][0]->Add(JetShape2[g+2][0][0],-1.);



      }
    }
  }

  cout<<"drawing"<<endl;

  TCanvas *PAS_plot = new TCanvas("JetShape_ForPAS","",1200,1200);
  PAS_plot->Divide(3,3,0.,0.);
 
  
  PAS_plot->cd(1);


  JetShape2[3][0][0]->SetMarkerStyle(4);
  JetShape2[3][0][0]->SetMarkerSize(1);
  JetShape2[3][0][0]->SetMarkerColor(kBlack);
  JetShape2[3][0][0]->SetLineColor(kBlack);
  JetShape2[3][0][0]->SetMaximum(rho_max);
  JetShape2[3][0][0]->SetMinimum(rho_min);
  JetShape2[3][0][0]->GetXaxis()->SetRangeUser(0.,.99);
  JetShape2[3][0][0]->GetYaxis()->SetLabelSize(0.06);
  JetShape2[3][0][0]->Draw();

  JetShape2[4][0][0]->SetMarkerStyle(4);
  JetShape2[4][0][0]->SetMarkerSize(1);
  JetShape2[4][0][0]->SetMarkerColor(kBlue);
  JetShape2[4][0][0]->SetLineColor(kBlue);
  JetShape2[4][0][0]->Draw("same");

  JetShape2[5][0][0]->SetMarkerStyle(10);
  JetShape2[5][0][0]->SetMarkerSize(1);
  JetShape2[5][0][0]->SetMarkerColor(kRed+1);
  JetShape2[5][0][0]->SetLineColor(kRed+1);
  //  JetShape2[5][0][0]->Draw("same");



 JetShape2[6][0][0]->SetMarkerStyle(4);
  JetShape2[6][0][0]->SetMarkerSize(1);
  JetShape2[6][0][0]->SetMarkerColor(kViolet);
  JetShape2[6][0][0]->SetLineColor(kViolet);
  JetShape2[6][0][0]->SetLineWidth(3);
  // JetShape2[6][0][0]->SetLineStyle(2);
  // JetShape2[6][0][0]->SetFillColor(kViolet);
  // JetShape2[6][0][0]->Draw("hist C same");
 JetShape2[6][0][0]->Draw("same");


  JetShape2[7][0][0]->SetMarkerStyle(4);
  JetShape2[7][0][0]->SetMarkerSize(1);
  JetShape2[7][0][0]->SetMarkerColor(kGreen);
  JetShape2[7][0][0]->SetLineColor(kGreen);
  JetShape2[7][0][0]->SetLineWidth(3);
  // JetShape2[7][0][0]->SetLineStyle(2);
  //JetShape2[7][0][0]->Draw("same hist C");
  JetShape2[7][0][0]->Draw("same ");


  JetShape2[3][0][0]->Draw("same");
  JetShape2[4][0][0]->Draw("same");

  JetShape2[9][0][0]->SetMarkerStyle(10);
  JetShape2[9][0][0]->SetMarkerSize(1);
  JetShape2[9][0][0]->SetMarkerColor(kCyan);
  JetShape2[9][0][0]->SetLineColor(kCyan);
  // JetShape2[9][0][0]->Draw("same");


  JetShape2[8][0][0]->SetMarkerStyle(10);
  JetShape2[8][0][0]->SetMarkerSize(1);
  JetShape2[8][0][0]->SetMarkerColor(kOrange);
  JetShape2[8][0][0]->SetLineColor(kOrange);
  //  JetShape2[8][0][0]->Draw("same");



  JetShape2[10][0][0]->SetMarkerStyle(10);
  JetShape2[10][0][0]->SetMarkerSize(1);
  JetShape2[10][0][0]->SetMarkerColor(kPink-9);
  JetShape2[10][0][0]->SetLineColor(kPink-9);
  //  JetShape2[10][0][0]->Draw("same");


  gPad->SetLogy();


  TLatex *label_pp = new TLatex(0.3,.9,"pp");
  label_pp->SetNDC();
  label_pp->Draw();


  TLegend *legend = new TLegend(0.35,0.55,0.9,0.95);
  legend->AddEntry(JetShape2[3][0][0],"5.02 TeV Inclusive");
  // legend->AddEntry(JetShape2[5][0][0],"2.76 TeV Inclusive");
  // legend->AddEntry(JetShape2[10][0][0],"5.02 TeV Leading");
  legend->AddEntry(JetShape2[4][0][0],"2.76 TeV Leading");

  
  legend->AddEntry(JetShape2[6][0][0],"5.02 TeV Inclusive Gen Pythia","lp");
  // legend->AddEntry(JetShape2[9][0][0],"5.02 TeV Leading Gen Pythia");
  // legend->AddEntry(JetShape2[8][0][0],"2.76 TeV Inclusive Gen Pythia");
  legend->AddEntry(JetShape2[7][0][0],"2.76 TeV Leading Gen Pythia","lp");
 
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kWhite);
  legend->SetTextSize(0.05);
  legend->Draw();

  PAS_plot->cd(2);

  JetShape2[0][3][0]->SetMarkerStyle(10);
  JetShape2[0][3][0]->SetMarkerSize(1);
  JetShape2[0][3][0]->SetMarkerColor(kBlack);
  JetShape2[0][3][0]->SetLineColor(kBlack);
  JetShape2[0][3][0]->SetMaximum(rho_max);
  JetShape2[0][3][0]->SetMinimum(rho_min);
  JetShape2[0][3][0]->GetXaxis()->SetRangeUser(0.,.99);
  JetShape2[0][3][0]->GetYaxis()->SetLabelSize(0);
  JetShape2[0][3][0]->Draw();

  JetShape2[1][3][0]->SetMarkerStyle(10);
  JetShape2[1][3][0]->SetMarkerSize(1);
  JetShape2[1][3][0]->SetMarkerColor(kBlue);
  JetShape2[1][3][0]->SetLineColor(kBlue);
  JetShape2[1][3][0]->Draw("same");


  JetShape2[2][3][0]->SetMarkerStyle(10);
  JetShape2[2][3][0]->SetMarkerSize(1);
  JetShape2[2][3][0]->SetMarkerColor(kRed);
  JetShape2[2][3][0]->SetLineColor(kRed);
  // JetShape2[2][3][0]->Draw("same");

  gPad->SetLogy();

  TLatex *label_per = new TLatex(0.1,.9,"PbPb Cent. 50-100%");
  label_per->SetNDC();
  label_per->Draw();



  PAS_plot->cd(3);

  JetShape2[0][0][0]->SetMarkerStyle(10);
  JetShape2[0][0][0]->SetMarkerSize(1);
  JetShape2[0][0][0]->SetMarkerColor(kBlack);
  JetShape2[0][0][0]->SetLineColor(kBlack);
  JetShape2[0][0][0]->SetMaximum(rho_max);
  JetShape2[0][0][0]->SetMinimum(rho_min);
  JetShape2[0][0][0]->GetXaxis()->SetRangeUser(0.,.99);
  JetShape2[0][0][0]->GetYaxis()->SetLabelSize(0);
  JetShape2[0][0][0]->Draw();

  JetShape2[1][0][0]->SetMarkerStyle(10);
  JetShape2[1][0][0]->SetMarkerSize(1);
  JetShape2[1][0][0]->SetMarkerColor(kBlue);
  JetShape2[1][0][0]->SetLineColor(kBlue);
  JetShape2[1][0][0]->Draw("same");


  JetShape2[2][0][0]->SetMarkerStyle(10);
  JetShape2[2][0][0]->SetMarkerSize(1);
  JetShape2[2][0][0]->SetMarkerColor(kRed);
  JetShape2[2][0][0]->SetLineColor(kRed);
  // JetShape2[2][0][0]->Draw("same");

  gPad->SetLogy();

TLatex *label_cent = new TLatex(0.1,.9,"PbPb Cent. 0-10%");
  label_cent->SetNDC();
  label_cent->Draw();




  PAS_plot->cd(5);

  JetShape_ratio[3][3][0]->SetMarkerStyle(10);
  JetShape_ratio[3][3][0]->SetMarkerSize(1);
  JetShape_ratio[3][3][0]->SetMarkerColor(kBlack);
  JetShape_ratio[3][3][0]->SetLineColor(kBlack);
 JetShape_ratio[3][3][0]->SetMaximum(ratio_max);
  JetShape_ratio[3][3][0]->SetMinimum(ratio_min);
   JetShape_ratio[3][3][0]->GetXaxis()->SetRangeUser(0.,.99);
   JetShape_ratio[3][3][0]->GetYaxis()->SetLabelSize(0.06);
  JetShape_ratio[3][3][0]->Draw();

  JetShape_ratio[4][3][0]->SetMarkerStyle(10);
  JetShape_ratio[4][3][0]->SetMarkerSize(1);
  JetShape_ratio[4][3][0]->SetMarkerColor(kBlue);
  JetShape_ratio[4][3][0]->SetLineColor(kBlue);
  JetShape_ratio[4][3][0]->Draw("same");


  JetShape_ratio[5][3][0]->SetMarkerStyle(10);
  JetShape_ratio[5][3][0]->SetMarkerSize(1);
  JetShape_ratio[5][3][0]->SetMarkerColor(kRed);
  JetShape_ratio[5][3][0]->SetLineColor(kRed);
  //JetShape_ratio[5][3][0]->Draw("same");



  TLatex *label_ratio = new TLatex(0.1,.9,"PbPb / pp");
  label_ratio->SetNDC();
  label_ratio->Draw();

  TLine *ratio_line = new TLine(0.,1.,1.,1.);
  ratio_line->SetLineStyle(2);
  ratio_line->Draw();


  PAS_plot->cd(6);

  JetShape_ratio[3][0][0]->SetMarkerStyle(10);
  JetShape_ratio[3][0][0]->SetMarkerSize(1);
  JetShape_ratio[3][0][0]->SetMarkerColor(kBlack);
  JetShape_ratio[3][0][0]->SetLineColor(kBlack);
  JetShape_ratio[3][0][0]->SetMaximum(ratio_max);
  JetShape_ratio[3][0][0]->SetMinimum(ratio_min);
   JetShape_ratio[3][0][0]->GetXaxis()->SetRangeUser(0.,.99);
  JetShape_ratio[3][0][0]->Draw();

  JetShape_ratio[4][0][0]->SetMarkerStyle(10);
  JetShape_ratio[4][0][0]->SetMarkerSize(1);
  JetShape_ratio[4][0][0]->SetMarkerColor(kBlue);
  JetShape_ratio[4][0][0]->SetLineColor(kBlue);
  JetShape_ratio[4][0][0]->Draw("same");


  JetShape_ratio[5][0][0]->SetMarkerStyle(10);
  JetShape_ratio[5][0][0]->SetMarkerSize(1);
  JetShape_ratio[5][0][0]->SetMarkerColor(kRed);
  JetShape_ratio[5][0][0]->SetLineColor(kRed);
  // JetShape_ratio[5][0][0]->Draw("same");

  //label_ratio->Draw();


  ratio_line->Draw();

  PAS_plot->cd(8);

  JetShape_diff[3][3][0]->SetMarkerStyle(10);
  JetShape_diff[3][3][0]->SetMarkerSize(1);
  JetShape_diff[3][3][0]->SetMarkerColor(kBlack);
  JetShape_diff[3][3][0]->SetLineColor(kBlack);
 JetShape_diff[3][3][0]->SetMaximum(diff_max);
  JetShape_diff[3][3][0]->SetMinimum(diff_min);
  JetShape_diff[3][3][0]->GetXaxis()->SetRangeUser(0.,.99);
  JetShape_diff[3][3][0]->GetYaxis()->SetLabelSize(0.06);
  JetShape_diff[3][3][0]->GetXaxis()->SetLabelSize(0.06);
  JetShape_diff[3][3][0]->GetXaxis()->SetTitleSize(0.06);
  JetShape_diff[3][3][0]->GetXaxis()->CenterTitle();
  JetShape_diff[3][3][0]->GetXaxis()->SetTitle("#Delta r");
  JetShape_diff[3][3][0]->Draw();

  JetShape_diff[4][3][0]->SetMarkerStyle(10);
  JetShape_diff[4][3][0]->SetMarkerSize(1);
  JetShape_diff[4][3][0]->SetMarkerColor(kBlue);
  JetShape_diff[4][3][0]->SetLineColor(kBlue);
  JetShape_diff[4][3][0]->Draw("same");


  JetShape_diff[5][3][0]->SetMarkerStyle(10);
  JetShape_diff[5][3][0]->SetMarkerSize(1);
  JetShape_diff[5][3][0]->SetMarkerColor(kRed);
  JetShape_diff[5][3][0]->SetLineColor(kRed);
  //  JetShape_diff[5][3][0]->Draw("same");



  TLatex *label_diff = new TLatex(0.1,.9,"PbPb - pp");
  label_diff->SetNDC();
  label_diff->Draw();

 TLine *diff_line = new TLine(0.,0.,1.,0.);
  diff_line->SetLineStyle(2);
  diff_line->Draw();





  PAS_plot->cd(9);

  JetShape_diff[3][0][0]->SetMarkerStyle(10);
  JetShape_diff[3][0][0]->SetMarkerSize(1);
  JetShape_diff[3][0][0]->SetMarkerColor(kBlack);
  JetShape_diff[3][0][0]->SetLineColor(kBlack);
  JetShape_diff[3][0][0]->SetMaximum(diff_max);
  JetShape_diff[3][0][0]->SetMinimum(diff_min);
   JetShape_diff[3][0][0]->GetXaxis()->SetRangeUser(0.,.99);
   JetShape_diff[3][0][0]->GetYaxis()->SetLabelSize(0.0);
   JetShape_diff[3][0][0]->GetXaxis()->SetLabelSize(0.06);
 JetShape_diff[3][0][0]->GetXaxis()->SetTitleSize(0.06);
 JetShape_diff[3][0][0]->GetXaxis()->CenterTitle();
 JetShape_diff[3][0][0]->GetXaxis()->SetTitle("#Delta r");
  JetShape_diff[3][0][0]->Draw();

  JetShape_diff[4][0][0]->SetMarkerStyle(10);
  JetShape_diff[4][0][0]->SetMarkerSize(1);
  JetShape_diff[4][0][0]->SetMarkerColor(kBlue);
  JetShape_diff[4][0][0]->SetLineColor(kBlue);
  JetShape_diff[4][0][0]->Draw("same");


  JetShape_diff[5][0][0]->SetMarkerStyle(10);
  JetShape_diff[5][0][0]->SetMarkerSize(1);
  JetShape_diff[5][0][0]->SetMarkerColor(kRed);
  JetShape_diff[5][0][0]->SetLineColor(kRed);
  //  JetShape_diff[5][0][0]->Draw("same");

  diff_line->Draw();

  //label_diff->Draw();

   PAS_plot->SaveAs("Comparison_Run1_Run2_PbPb_pp.png");


  TCanvas *ratio_plot = new TCanvas("RatioPlot","",1200,900);
  ratio_plot->Divide(3,2,0.,0.);
  
  ratio_plot->cd(1);

  JetShape2[3][0][0]->Draw();
  JetShape2[4][0][0]->Draw("same");
  JetShape2[6][0][0]->Draw("same");
  JetShape2[7][0][0]->Draw("same");
  label_pp->Draw("same");
  gPad->SetLogy();

  legend->Draw();

  ratio_plot->cd(2);

  JetShape2[0][3][0]->Draw();
  JetShape2[1][3][0]->Draw("same");
 
label_per->Draw("same");

 gPad->SetLogy();

  ratio_plot->cd(3);

  JetShape2[0][0][0]->Draw();
  JetShape2[1][0][0]->Draw("same");
gPad->SetLogy();
 
 label_cent->Draw("same");

 


 TH1D *Ratio2_pp = (TH1D*)JetShape2[3][0][0]->Clone("ratio2_pp");

 TH1D *Ratio2_per = (TH1D*)JetShape2[0][3][0]->Clone("ratio2_pp");

 TH1D *Ratio2_cent = (TH1D*)JetShape2[0][0][0]->Clone("ratio2_pp");
  
 for(int k = 0; k<JetShape2[1][0][0]->GetNbinsX()+1; k++){
   
   Ratio2_pp->SetBinContent(k,JetShape2[3][0][0]->GetBinContent(k)/JetShape2[4][0][0]->GetBinContent(k));
   Ratio2_pp->SetBinError(k,TMath::Sqrt(JetShape2[3][0][0]->GetBinError(k)*JetShape2[3][0][0]->GetBinError(k)/JetShape2[3][0][0]->GetBinContent(k)/JetShape2[3][0][0]->GetBinContent(k)+JetShape2[4][0][0]->GetBinError(k)*JetShape2[4][0][0]->GetBinError(k)/JetShape2[4][0][0]->GetBinContent(k)/JetShape2[4][0][0]->GetBinContent(k)));

   Ratio2_per->SetBinContent(k,JetShape2[0][3][0]->GetBinContent(k)/JetShape2[1][3][0]->GetBinContent(k));
   Ratio2_per->SetBinError(k,TMath::Sqrt(JetShape2[0][3][0]->GetBinError(k)*JetShape2[0][3][0]->GetBinError(k)/JetShape2[0][3][0]->GetBinContent(k)/JetShape2[0][3][0]->GetBinContent(k)+JetShape2[1][3][0]->GetBinError(k)*JetShape2[1][3][0]->GetBinError(k)/JetShape2[1][3][0]->GetBinContent(k)/JetShape2[1][3][0]->GetBinContent(k)));

   Ratio2_cent->SetBinContent(k,JetShape2[0][0][0]->GetBinContent(k)/JetShape2[1][0][0]->GetBinContent(k));
   Ratio2_cent->SetBinError(k,TMath::Sqrt(JetShape2[0][0][0]->GetBinError(k)*JetShape2[0][0][0]->GetBinError(k)/JetShape2[0][0][0]->GetBinContent(k)/JetShape2[0][0][0]->GetBinContent(k)+JetShape2[1][0][0]->GetBinError(k)*JetShape2[1][0][0]->GetBinError(k)/JetShape2[1][0][0]->GetBinContent(k)/JetShape2[1][0][0]->GetBinContent(k)));


 }



 ratio_plot->cd(4);



 Ratio2_pp->SetMinimum(0.);
 Ratio2_pp->SetMaximum(6.5);

 Ratio2_pp->GetXaxis()->SetRangeUser(0.,0.99);
 Ratio2_pp->GetXaxis()->SetRangeUser(0.,.99);
  Ratio2_pp->GetYaxis()->SetLabelSize(0.0);
  Ratio2_pp->GetXaxis()->SetLabelSize(0.06);
Ratio2_pp->GetXaxis()->SetTitleSize(0.06);
Ratio2_pp->GetXaxis()->CenterTitle();
Ratio2_pp->GetXaxis()->SetTitle("#Delta r");

 Ratio2_pp->Draw();

 ratio_line->Draw();

 TLatex *ratio_tex = new TLatex (0.2,0.9,"5.02 Inclusive / 2.76 Leading");
 ratio_tex->SetLineColor(kWhite);
 ratio_tex->SetNDC();
 ratio_tex->Draw();


 ratio_plot->cd(5);


 Ratio2_per->SetMinimum(0.);
 Ratio2_per->SetMaximum(6.5);

 Ratio2_per->GetXaxis()->SetRangeUser(0.,0.99);
 Ratio2_per->GetXaxis()->SetRangeUser(0.,.99);
  Ratio2_per->GetYaxis()->SetLabelSize(0.0);
  Ratio2_per->GetXaxis()->SetLabelSize(0.06);
Ratio2_per->GetXaxis()->SetTitleSize(0.06);
Ratio2_per->GetXaxis()->CenterTitle();
Ratio2_per->GetXaxis()->SetTitle("#Delta r");


 Ratio2_per->GetXaxis()->SetRangeUser(0.,0.99);
 Ratio2_per->Draw();

 ratio_line->Draw();


  ratio_plot->cd(6);


 Ratio2_cent->SetMinimum(0.);
 Ratio2_cent->SetMaximum(6.5);

 Ratio2_cent->GetXaxis()->SetRangeUser(0.,0.99);
 Ratio2_cent->GetXaxis()->SetRangeUser(0.,.99);
  Ratio2_cent->GetYaxis()->SetLabelSize(0.0);
  Ratio2_cent->GetXaxis()->SetLabelSize(0.06);
Ratio2_cent->GetXaxis()->SetTitleSize(0.06);
Ratio2_cent->GetXaxis()->CenterTitle();
Ratio2_cent->GetXaxis()->SetTitle("#Delta r");


 Ratio2_cent->GetXaxis()->SetRangeUser(0.,0.99);
 Ratio2_cent->Draw();

 ratio_line->Draw();





  ratio_plot->SaveAs("Ratios_Run1_Run2_PbPb_pp.png");


  TCanvas *comparison_plot = new TCanvas("ComparisonPlot","",2000,500);
  comparison_plot->Divide(5,1,0.,0.);


  comparison_plot->cd(1);

  legend = new TLegend(0.05,0.55,0.9,0.95);
  legend->AddEntry(JetShape2[3][0][0],"5.02 TeV Inclusive pp Data");
  legend->AddEntry(JetShape2[5][0][0],"2.76 TeV Inclusive pp Data");
  legend->AddEntry(JetShape2[10][0][0],"5.02 TeV Leading pp Data");
  legend->AddEntry(JetShape2[4][0][0],"2.76 TeV Leading pp Data");

  legend->AddEntry(JetShape2[6][0][0],"5.02 TeV Inclusive Gen Pythia");
  legend->AddEntry(JetShape2[9][0][0],"5.02 TeV Leading Gen Pythia");
  legend->AddEntry(JetShape2[8][0][0],"2.76 TeV Inclusive Gen Pythia");
  legend->AddEntry(JetShape2[7][0][0],"2.76 TeV Leading Gen Pythia");
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kWhite);
  legend->SetTextSize(0.05);
  legend->Draw();


  comparison_plot->cd(2);


  JetShape2[4][0][0]->SetMarkerStyle(10);
  JetShape2[4][0][0]->SetMarkerSize(1);
  JetShape2[4][0][0]->SetMarkerColor(kBlue);
  JetShape2[4][0][0]->SetLineColor(kBlue);
  JetShape2[4][0][0]->SetMaximum(rho_max);
  JetShape2[4][0][0]->SetMinimum(rho_min);
  JetShape2[4][0][0]->GetXaxis()->SetRangeUser(0.,.99);
  JetShape2[4][0][0]->Draw();

  JetShape2[5][0][0]->Draw("same");

  JetShape2[7][0][0]->Draw("same");

  JetShape2[8][0][0]->Draw("same");

  gPad->SetLogy();

  TLatex *label1 = new TLatex(0.3,0.9,"2.76 TeV");
  label1->SetNDC();
  label1->Draw();


  comparison_plot->cd(3);


  JetShape2[3][0][0]->Draw("same");

  JetShape2[6][0][0]->Draw("same");

  JetShape2[9][0][0]->Draw("same");
JetShape2[10][0][0]->Draw("same");

  gPad->SetLogy();

 TLatex *label2 = new TLatex(0.3,0.9,"5.02 TeV");
  label2->SetNDC();
  label2->Draw();



 comparison_plot->cd(4);

 JetShape2[3][0][0]->Draw();
 JetShape2[5][0][0]->Draw("same");
 JetShape2[6][0][0]->Draw("same");
 JetShape2[8][0][0]->Draw("same");

 

 gPad->SetLogy();

 TLatex *label3 = new TLatex(0.3,0.9,"All Inclusive");
  label3->SetNDC();
  label3->Draw();




 comparison_plot->cd(5);

 JetShape2[4][0][0]->Draw();
 JetShape2[7][0][0]->Draw("same");
 JetShape2[9][0][0]->Draw("same");
JetShape2[10][0][0]->Draw("same");

 TLatex *label4 = new TLatex(0.3,0.9,"All Leading");
  label4->SetNDC();
  label4->Draw();



 gPad->SetLogy();

  comparison_plot->SaveAs("Pythia_pp_Studies.png");


  return 0;
}// main loop
