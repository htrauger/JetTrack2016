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

Int_t spill_jet_studies(bool do_corr = kTRUE){

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.25);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);

  TCanvas *c_yields_phi, *c_yields_eta, *c_yields_eta_diff, *c_double_diff;

  const int nCBins = 5;
  const int nPtBins = 1;
  const int nTrkPtBins = 9;



  TH2D *result[12][nTrkPtBins][nCBins];
  TH2D *result2[12][nTrkPtBins][nCBins];

  TH1D *signal_dPhi[12][nTrkPtBins][nCBins];
  
  TH1D *signal_dEta[12][nTrkPtBins][nCBins];
  TH1D *spill_over_dPhi[12][nTrkPtBins][nCBins];
  TH1D *spill_over_dEta[12][nTrkPtBins][nCBins];
  TH1D *signal_dEta_diff[12][nTrkPtBins][nCBins];

  TH1D *jff_residual_dEta[12][nTrkPtBins][nCBins];

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%","pp"};

   float TrkPtBins[nTrkPtBins+1] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
  
  TString Type_strs[4] = {"Nominal","In","Out","Out_E5"};

  TFile *f_in_pbpb = new TFile("../me_correct/PbPb_Inclusive_Correlations.root");
  TFile *f_in_pbpb_in = new TFile("../me_correct/PbPb_Inclusive_Correlations_SpilledIn.root");
  TFile *f_in_pbpb_out_e5 = new TFile("../me_correct/PbPb_Inclusive_Correlations_SpilledOut_Eta0p5.root");
  TFile *f_in_pbpb_out = new TFile("../me_correct/PbPb_Inclusive_Correlations_SpilledOut.root");
 
  TFile *f_in_sube0 = new TFile("../me_correct/Hydjet_RecoJet_GenTrack_NoSube0_Inclusive_Correlations.root");
  TFile *f_in_sube0_in = new TFile("../me_correct/Hydjet_RecoJet_GenTrack_NoSube0_Inclusive_Correlations_SpillInOnly.root");
  TFile *f_in_sube0_out = new TFile("../me_correct/Hydjet_RecoJet_GenTrack_NoSube0_Inclusive_Correlations_SpillOutOnly.root");


  TFile *f_spill_over = new TFile("../spill_over/Inclusive_Hydjet_SpillOvers.root");
 
  /*
  TFile *f_spill_over_preapp = new TFile("../FROZEN_PREAPPROVAL/spill_over/Inclusive_Hydjet_SpillOvers.root");
  TFile *f_spill_over_qm = new TFile("../FROZEN_QM/spill_over/Inclusive_Hydjet_SpillOvers.root");

  */
  TString in_name, pTlabel,centlabel,Ajlabel;

  int lbin, rbin;

  double bin_width_phi, bin_width_eta, diff_max, diff_min, signal_min, signal_max, stacked_max, stacked_max_diff, stacked_min_diff, stacked_min, bc, err;
    
  c_yields_phi = new TCanvas("yields_phi","",10,10,3200,3200);
  c_yields_phi->Divide(4,8,0,0);

   
  c_yields_eta = new TCanvas("yields_eta","",10,10,3200,3200);
  c_yields_eta->Divide(4,8,0,0);


  for(int i = 0; i<nTrkPtBins-1; i++){

    for(int j = 0; j<4; j++){

      for(int g = 0; g<3; g++){

	
	cout<<g<<" "<<i<<" "<<j<<endl;

	if(g==0){
	  result[g][i][j] = (TH2D*)f_in_pbpb->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	  result2[g][i][j] = (TH2D*)f_in_sube0->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_Sube0_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
	
	}else if(g==1){
	  
	  result[g][i][j] = (TH2D*)f_in_pbpb_in->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_SpillIn_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	  result2[g][i][j] = (TH2D*)f_in_sube0_in->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_SpillInSube0_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	}else if(g==2){
	  
	  result[g][i][j] = (TH2D*)f_in_pbpb_out->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_SpillOut_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	  result2[g][i][j] = (TH2D*)f_in_sube0_out->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_SpillOutSube0_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	}else if(g==3){
	    
	    result[g][i][j] = (TH2D*)f_in_pbpb_out_e5->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_SpillIn_Etap5_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
	}else{

	  cout<<"not a sample!"<<endl;

	  return -1;
	}

	if(i>3){
	  result[g][i][j]->Scale(1./4.);
	  result2[g][i][j]->Scale(1./4.);
	}else if(i==0){
	  result[g][i][j]->Scale(1/.3);
	  result2[g][i][j]->Scale(1/.3);
	}

	int lbin = result[g][i][j]->GetXaxis()->FindBin(-.99);
	int rbin = result[g][i][j]->GetXaxis()->FindBin(.99);

	signal_dPhi[g][i][j] = (TH1D*)result[g][i][j]->ProjectionY((TString)("Proj_dPhi_" +Type_strs[g]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);

	signal_dPhi[g][i][j]->SetMarkerStyle(10);
	signal_dPhi[g][i][j]->SetMarkerSize(1);

	signal_dPhi[g][i][j]->Rebin(5);

	spill_over_dPhi[g][i][j] = (TH1D*)result2[g][i][j]->ProjectionY((TString)("SpillOver_dPhi_" +Type_strs[g]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);

	spill_over_dPhi[g][i][j]->SetMarkerStyle(4);
	spill_over_dPhi[g][i][j]->SetMarkerSize(1);

	spill_over_dPhi[g][i][j]->Rebin(5);


	float norm = 	signal_dPhi[g][i][j]->GetBinWidth(1);

	signal_dPhi[g][i][j]->Scale(1./norm);
	spill_over_dPhi[g][i][j]->Scale(1./norm);


	lbin = result[g][i][j]->GetYaxis()->FindBin(-.99);
	rbin = result[g][i][j]->GetYaxis()->FindBin(.99);

	signal_dEta[g][i][j] = (TH1D*)result[g][i][j]->ProjectionX((TString)("Proj_dEta_" +Type_strs[g]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);

	spill_over_dEta[g][i][j] = (TH1D*)result2[g][i][j]->ProjectionX((TString)("SpillOver_dEta_" +Type_strs[g]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);

	signal_dEta[g][i][j]->SetMarkerStyle(10);
	signal_dEta[g][i][j]->SetMarkerSize(1);

	signal_dEta[g][i][j]->Rebin(5);

	spill_over_dEta[g][i][j]->SetMarkerStyle(4);
	spill_over_dEta[g][i][j]->SetMarkerSize(1);

	spill_over_dEta[g][i][j]->Rebin(5);


	norm = 	signal_dEta[g][i][j]->GetBinWidth(1);

	signal_dEta[g][i][j]->Scale(1./norm);
	spill_over_dEta[g][i][j]->Scale(1./norm);


	if(do_corr){

	  signal_dEta[g][i][j]->Add(spill_over_dEta[g][i][j],-1.);
	  signal_dPhi[g][i][j]->Add(spill_over_dPhi[g][i][j],-1.);
	}

      }
      
    }
  
    cout<<"ready to draw"<<endl;
    for(int j = 0; j<4; j++){
      c_yields_phi->cd(4*(i+1)-j);

      signal_dPhi[0][i][j]->SetLineColor(kViolet);
      signal_dPhi[0][i][j]->SetMarkerColor(kViolet);

      signal_dPhi[1][i][j]->SetLineColor(kRed);
      signal_dPhi[1][i][j]->SetMarkerColor(kRed);

      signal_dPhi[2][i][j]->SetLineColor(kBlue);
      signal_dPhi[2][i][j]->SetMarkerColor(kBlue);

      //    signal_dPhi[3][i][j]->SetLineColor(kGreen);
      // signal_dPhi[3][i][j]->SetMarkerColor(kGreen);

      spill_over_dPhi[0][i][j]->SetLineColor(kViolet);
      spill_over_dPhi[0][i][j]->SetMarkerColor(kViolet);

      spill_over_dPhi[1][i][j]->SetLineColor(kRed);
      spill_over_dPhi[1][i][j]->SetMarkerColor(kRed);

      spill_over_dPhi[2][i][j]->SetLineColor(kBlue);
      spill_over_dPhi[2][i][j]->SetMarkerColor(kBlue);



      double ymax = 23.;
      double ymin = -2.;

      if(i>2){
	ymax = 10.;
	ymin = -1.;
      }
      if(i>4){
	ymax = 4.;
	ymin = -1.;
      }

      //  ymax = 3.;


      signal_dPhi[0][i][j]->SetMaximum(ymax);
      signal_dPhi[0][i][j]->SetMinimum(ymin);

      signal_dPhi[0][i][j]->GetXaxis()->SetRangeUser(-TMath::Pi()/2., 3.*TMath::Pi()/2.);
     
      if(j==3)  signal_dPhi[0][i][j]->GetYaxis()->SetLabelSize(0.1);
      else   signal_dPhi[0][i][j]->GetYaxis()->SetLabelSize(0.);

      signal_dPhi[0][i][j]->GetXaxis()->SetLabelSize(0.1);

      
      signal_dPhi[0][i][j]->Draw();
      signal_dPhi[1][i][j]->Draw("same");
      signal_dPhi[2][i][j]->Draw("same");
      /*
      spill_over_dPhi[0][i][j]->Draw("same");
      spill_over_dPhi[1][i][j]->Draw("same");
      spill_over_dPhi[2][i][j]->Draw("same");
      //signal_dPhi[3][i][j]->Draw("same");
      */
      TLine *zero = new TLine(-TMath::Pi()/2.,0.,3.*TMath::Pi()/2.,0.);
      zero->SetLineStyle(2);
      zero->Draw("same");

      if(i==0&&j==3){
	TLegend *l = new TLegend(0.5,0.5,0.9,0.9);
	l->AddEntry(signal_dPhi[0][i][j],"Nominal");
	l->AddEntry(signal_dPhi[1][i][j],"Spilled In");
	l->AddEntry(signal_dPhi[2][i][j],"Spilled Out");
	//	l->AddEntry(signal_dPhi[3][i][j],"Spilled Out |#eta_{jet}|<0.5");
	l->SetLineColor(kWhite);
	l->SetTextSize(0.08);
	l->Draw();
      }

      if(j<3){
	  
	TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[j]);
	centtex->SetNDC();
	centtex->Draw();
	 
	TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[i]);
	pttex->SetNDC();
	pttex->Draw();
     
      }else{

	TLatex *centtex = new TLatex(0.2,0.9,CBin_labels[j]);
	centtex->SetNDC();
	centtex->Draw();
	
	TLatex *pttex = new TLatex(0.2,0.85,TrkPtBin_labels[i]);
	pttex->SetNDC();
	pttex->Draw();
      }


  c_yields_eta->cd(4*(i+1)-j);

      signal_dEta[0][i][j]->SetLineColor(kViolet);
      signal_dEta[0][i][j]->SetMarkerColor(kViolet);

      signal_dEta[1][i][j]->SetLineColor(kRed);
      signal_dEta[1][i][j]->SetMarkerColor(kRed);

      signal_dEta[2][i][j]->SetLineColor(kBlue);
      signal_dEta[2][i][j]->SetMarkerColor(kBlue);


      spill_over_dEta[0][i][j]->SetLineColor(kViolet);
      spill_over_dEta[0][i][j]->SetMarkerColor(kViolet);

      spill_over_dEta[1][i][j]->SetLineColor(kRed);
      spill_over_dEta[1][i][j]->SetMarkerColor(kRed);

      spill_over_dEta[2][i][j]->SetLineColor(kBlue);
      spill_over_dEta[2][i][j]->SetMarkerColor(kBlue);

      //   signal_dEta[3][i][j]->SetLineColor(kGreen);
      // signal_dEta[3][i][j]->SetMarkerColor(kGreen);

      signal_dEta[0][i][j]->SetMaximum(ymax);
      signal_dEta[0][i][j]->SetMinimum(ymin);

      signal_dEta[0][i][j]->GetXaxis()->SetRangeUser(-2.5,2.5);
     
      if(j==3)  signal_dEta[0][i][j]->GetYaxis()->SetLabelSize(0.1);
      else   signal_dEta[0][i][j]->GetYaxis()->SetLabelSize(0.);


      signal_dEta[0][i][j]->GetXaxis()->SetLabelSize(0.1);
      
      signal_dEta[0][i][j]->Draw();
      signal_dEta[1][i][j]->Draw("same");
      signal_dEta[2][i][j]->Draw("same");
      /*
      spill_over_dEta[0][i][j]->Draw("same");
      spill_over_dEta[1][i][j]->Draw("same");
      spill_over_dEta[2][i][j]->Draw("same");
      //   signal_dEta[3][i][j]->Draw("same");
      */
      TLine *zero2 = new TLine(-2.5,0.,2.5,0.);
      zero2->SetLineStyle(2);
      zero2->Draw("same");


      if(i==0&&j==3){
	TLegend *l = new TLegend(0.6,0.5,0.9,0.9);
	l->AddEntry(signal_dEta[0][i][j],"Nominal");
	l->AddEntry(signal_dEta[1][i][j],"Spilled In");
	l->AddEntry(signal_dEta[2][i][j],"Spilled Out");
	//	l->AddEntry(signal_dEta[3][i][j],"Spilled Out |#eta_{jet}|<0.5");
	l->SetLineColor(kWhite);
	l->SetTextSize(0.08);
	l->Draw();
      }

      if(j<3){
	  
	TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[j]);
	centtex->SetNDC();
	centtex->Draw();
	 
	TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[i]);
	pttex->SetNDC();
	pttex->Draw();
     
      }else{

	TLatex *centtex = new TLatex(0.2,0.9,CBin_labels[j]);
	centtex->SetNDC();
	centtex->Draw();
	
	TLatex *pttex = new TLatex(0.2,0.85,TrkPtBin_labels[i]);
	pttex->SetNDC();
	pttex->Draw();
      }
    }
  }

  if(do_corr){
    c_yields_eta->SaveAs("Spill_Over_Overlay_dEta_Corrected.png");
    c_yields_phi->SaveAs("Spill_Over_Overlay_dPhi_Corrected.png");
  }else{
    c_yields_eta->SaveAs("Spill_Over_Overlay_dEta.png");
    c_yields_phi->SaveAs("Spill_Over_Overlay_dPhi.png");
  }
  return 0;
}
