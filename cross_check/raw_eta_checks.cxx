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

Int_t raw_eta_checks(bool spill_in_out = kFALSE){

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.25);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);

  TCanvas *c_yields_eta, *c_yields_eta_diff, *c_double_diff;

  const int nCBins = 5;
  const int nPtBins = 1;
  const int nTrkPtBins = 9;



  TH2D *result[12][nTrkPtBins][nCBins];
  TH2D *result2[12][nTrkPtBins][nCBins];

  TH1D *signal_dPhi[12][nTrkPtBins][nCBins];
  TH1D *signal_dEta[12][nTrkPtBins][nCBins];
  TH1D *spill_over_dEta[12][nTrkPtBins][nCBins];
  TH1D *signal_dEta_diff[12][nTrkPtBins][nCBins];

  TH1D *jff_residual_dEta[12][nTrkPtBins][nCBins];

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%","pp"};

   float TrkPtBins[nTrkPtBins+1] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
  
  TString Type_strs[4] = {"Latest","Preapp","QM","Aug. Re-run"};

  TFile *f_in_pbpb = new TFile("../me_correct/PbPb_Inclusive_Correlations.root");
  TFile *f_in_pbpb_old = new TFile("../me_correct/PbPb_Inclusive_Correlations_Old.root");
  TFile *f_in_pbpb_preapp = new TFile("../FROZEN_PREAPPROVAL/me_correct/PbPb_Inclusive_Correlations.root");
  TFile *f_in_pbpb_qm = new TFile("../FROZEN_QM/me_correct/PbPb_Inclusive_Correlations.root");


  TFile *f_in_pp = new TFile("../me_correct/pp_Inclusive_Correlations.root");
  TFile *f_in_pp_preapp = new TFile("../FROZEN_PREAPPROVAL/me_correct/pp_Inclusive_Correlations.root");
  TFile *f_in_pp_qm = new TFile("../FROZEN_QM/me_correct/pp_Inclusive_Correlations.root");


  TFile *f_spill_over = new TFile("../spill_over/Inclusive_Hydjet_SpillOvers.root");
  TFile *f_spill_over_preapp = new TFile("../FROZEN_PREAPPROVAL/spill_over/Inclusive_Hydjet_SpillOvers.root");
  TFile *f_spill_over_qm = new TFile("../FROZEN_QM/spill_over/Inclusive_Hydjet_SpillOvers.root");


  TString in_name, pTlabel,centlabel,Ajlabel;

  int lbin, rbin;

  TFile *f_jff_hyd = new TFile("../jff_residual/Inclusive_Hydjet_JFFResiduals.root");
  TFile *f_jff_preapp = new TFile("../mc_raw_correlations/JFFcorrs_sube0_runWithoutJFFs_drumTune.root");


  TFile *f_jff_pyth = new TFile("../jff_residual/Inclusive_Pythia_JFFResiduals.root");
  TFile *f_jff_pyth_preapp = new TFile("../FROZEN_PREAPPROVAL/jff_residual/Inclusive_Pythia_JFFResiduals.root");


  double bin_width_phi, bin_width_eta, diff_max, diff_min, signal_min, signal_max, stacked_max, stacked_max_diff, stacked_min_diff, stacked_min, bc, err;
    
  c_yields_eta = new TCanvas("yields_eta","",10,10,2000,2400);
  c_yields_eta->Divide(5,8,0,0);

  c_yields_eta_diff = new TCanvas("yields_eta_diff","",10,10,2000,2400);
  c_yields_eta_diff->Divide(5,8,0,0);



  c_double_diff = new TCanvas("double_diff","",10,10,2000,2400);
  c_double_diff->Divide(5,8,0,0);

  for(int i = 0; i<nTrkPtBins-1; i++){

    for(int j = 0; j<4; j++){

      for(int g = 0; g<4; g++){

	cout<<"starting"<<endl;

	
	if(g==0) result[g][i][j] = (TH2D*)f_in_pbpb->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	else if(g==1) result[g][i][j] = (TH2D*)f_in_pbpb_preapp->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_PreApp_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	else if(g==2) result[g][i][j] = (TH2D*)f_in_pbpb_qm->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_QM_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	else result[g][i][j] = (TH2D*)f_in_pbpb_old->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_Aug_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	if(i>3){
	  result[g][i][j]->Scale(1./4.);
	}else if(i==0){
	  result[g][i][j]->Scale(1/.3);
	}

	int lbin = result[g][i][j]->GetYaxis()->FindBin(-.99);
	int rbin = result[g][i][j]->GetYaxis()->FindBin(.99);
	signal_dEta[g][i][j] = (TH1D*)result[g][i][j]->ProjectionX((TString)("Proj_" +Type_strs[g]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);

	signal_dEta[g][i][j]->SetMarkerStyle(10);
	signal_dEta[g][i][j]->SetMarkerSize(1);

	float norm = 	signal_dEta[g][i][j]->GetBinWidth(1);

	signal_dEta[g][i][j]->Scale(1./norm);
	signal_dEta[g][i][j] = Rebin_dEta(signal_dEta[g][i][j]);


	
	TString	spillover_name =(TString)("Eta_SpillOver_Points_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]);
	cout<<spillover_name<<endl;
	
	if(g==1) 	spill_over_dEta[g][i][j] = (TH1D*)f_spill_over_preapp->Get(spillover_name)->Clone(spillover_name);
	else if(g==2) 	spill_over_dEta[g][i][j] = (TH1D*)f_spill_over_qm->Get(spillover_name)->Clone(spillover_name);
	else 	spill_over_dEta[g][i][j] = (TH1D*)f_spill_over->Get(spillover_name)->Clone(spillover_name);

	if(!spill_in_out)	signal_dEta[g][i][j]->Add(spill_over_dEta[g][i][j],-1.);

      }
    
      cout<<"here"<<endl;

      TString	jff_name =(TString)("JFF_Residual_Eta_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]);
      jff_residual_dEta[0][i][j] = (TH1D*)f_jff_hyd->Get(jff_name)->Clone(jff_name);

      jff_name.ReplaceAll("Pt1000","Pt300");

      cout<<"here"<<endl;
      result2[1][i][j] = (TH2D*)f_jff_preapp->Get(Form("JFFcorrs_cent%d_pt%d",j,i+1))->Clone(Form("JFFcorrs_cent%d_pt%d",j,i+1));


      if(i>3){
	result2[1][i][j]->Scale(1./4.);
      }else if(i==0){
	result2[1][i][j]->Scale(1/.3);
      }


      cout<<"and here"<<endl;
      lbin = result2[1][i][j]->GetYaxis()->FindBin(-.99);
      rbin = result2[1][i][j]->GetYaxis()->FindBin(.99);
      
      jff_residual_dEta[1][i][j] = (TH1D*)result2[1][i][j]->ProjectionX((TString)("JFF_" +Type_strs[1]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);

      jff_residual_dEta[1][i][j]->SetMarkerStyle(10);
      jff_residual_dEta[1][i][j]->SetMarkerSize(1);

      float norm = jff_residual_dEta[1][i][j]->GetBinWidth(1);

      jff_residual_dEta[1][i][j]->Scale(1./norm);
      jff_residual_dEta[1][i][j] = Rebin_dEta( jff_residual_dEta[1][i][j]);

      jff_name+="Double_Diff";
      jff_residual_dEta[3][i][j]= (TH1D*)jff_residual_dEta[0][i][j]->Clone((TString)(jff_name+"Diff"));
      jff_residual_dEta[3][i][j]->Add(jff_residual_dEta[1][i][j],-1.);

    }
  
    cout<<"starting pp"<<endl;

    for(int g = 0; g<3; g++){
      if(g==0) result[g][i][4] = (TH2D*)f_in_pp->Get((TString)("Yield_BkgSub_"+CBin_strs[0]+"_"+CBin_strs[1]+"_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
      cout<<"0"<<endl;
      if(g==1) result[g][i][4] = (TH2D*)f_in_pp_preapp->Get((TString)("Yield_BkgSub_"+CBin_strs[0]+"_"+CBin_strs[1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_PreApp_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
      cout<<"1"<<endl;
      if(g==2) result[g][i][4] = (TH2D*)f_in_pp_qm->Get((TString)("Yield_BkgSub_"+CBin_strs[0]+"_"+CBin_strs[1]+"_Pt100_Pt300_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Yield_BkgSub_QM_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
      cout<<"2"<<endl;

      if(i>3){
	result[g][i][4]->Scale(1./4.);
      }else if(i==0){
	result[g][i][4]->Scale(1/.3);
      }
  
      int lbin = result[g][i][4]->GetYaxis()->FindBin(-.99);
      int rbin = result[g][i][4]->GetYaxis()->FindBin(.99);

      signal_dEta[g][i][4] = (TH1D*)result[g][i][4]->ProjectionX((TString)("Proj_pp_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);
  
      signal_dEta[g][i][4]->SetMarkerStyle(10);
      signal_dEta[g][i][4]->SetMarkerSize(1);
  
      float norm = 	signal_dEta[g][i][4]->GetBinWidth(1);

      signal_dEta[g][i][4]->Scale(1./norm);
      signal_dEta[g][i][4] = Rebin_dEta(signal_dEta[g][i][4]);


      cout<<"here"<<endl;
      TString	jff_name =(TString)("JFF_Residual_Eta_"+CBin_strs[0]+"_"+CBin_strs[1]+"_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]);
      jff_residual_dEta[0][i][4] = (TH1D*)f_jff_pyth->Get(jff_name)->Clone((TString)(jff_name+"pp"));

      jff_name.ReplaceAll("Pt1000","Pt300");
      jff_residual_dEta[1][i][4] = (TH1D*)f_jff_pyth_preapp->Get(jff_name)->Clone((TString)(jff_name+"pp_Preapp"));

      /*
      jff_name+="Double_Diff";
      jff_residual_dEta[3][i][4]= (TH1D*)jff_residual_dEta[0][i][4]->Clone((TString)(jff_name+"Diff"));
      jff_residual_dEta[3][i][4]->Add(jff_residual_dEta[1][i][4],-1.);
      */

    }
  

    cout<<"ready to draw"<<endl;
    for(int j = 0; j<5; j++){
      c_yields_eta->cd(5*(i+1)-j);

      signal_dEta[0][i][j]->SetLineColor(kViolet);
      signal_dEta[0][i][j]->SetMarkerColor(kViolet);

      signal_dEta[1][i][j]->SetLineColor(kRed);
      signal_dEta[1][i][j]->SetMarkerColor(kRed);

      signal_dEta[2][i][j]->SetLineColor(kBlue);
      signal_dEta[2][i][j]->SetMarkerColor(kBlue);

      if(j<4){
	signal_dEta[3][i][j]->SetLineColor(kOrange);
	signal_dEta[3][i][j]->SetMarkerColor(kOrange);
      }

      signal_dEta[0][i][j]->SetMaximum( signal_dEta[0][i][j]->GetMaximum()+1.);
      signal_dEta[0][i][j]->SetMinimum(-1.5);

      signal_dEta[0][i][j]->SetAxisRange(-2.49,2.49);
     
      if(j==4)  signal_dEta[0][i][j]->GetYaxis()->SetLabelSize(0.06);
      else   signal_dEta[0][i][j]->GetYaxis()->SetLabelSize(0.);

      
      signal_dEta[0][i][j]->Draw();
      signal_dEta[1][i][j]->Draw("same");
      signal_dEta[2][i][j]->Draw("same");
      // if(j<4)  signal_dEta[3][i][j]->Draw("same");

      if(!spill_in_out){
	jff_residual_dEta[0][i][j]->SetMarkerStyle(10);
	jff_residual_dEta[0][i][j]->SetMarkerSize(1);
	jff_residual_dEta[0][i][j]->SetMarkerColor(kBlack);
	jff_residual_dEta[0][i][j]->SetLineColor(kBlack);
	jff_residual_dEta[0][i][j]->Draw("same");


	jff_residual_dEta[1][i][j]->SetMarkerStyle(10);
	jff_residual_dEta[1][i][j]->SetMarkerSize(1);
	jff_residual_dEta[1][i][j]->SetMarkerColor(kGreen);
	jff_residual_dEta[1][i][j]->SetLineColor(kGreen);
	jff_residual_dEta[1][i][j]->Draw("same");
      }

      TLine *zero = new TLine(-2.5,0.,2.5,0.);
      zero->SetLineStyle(2);
      zero->Draw("same");


      if(i==0&&j==4){
	TLegend *l = new TLegend(0.7,0.45,0.9,0.95);
	l->AddEntry(signal_dEta[0][i][j],"Latest Data");
	l->AddEntry(signal_dEta[1][i][j],"Preapproval Data");
	l->AddEntry(signal_dEta[2][i][j],"QM Data");
	//  	l->AddEntry(signal_dEta[3][i][3],"Aug. Data, New MC");
	if(!spill_in_out)	l->AddEntry(jff_residual_dEta[0][i][j],"Latest JFF");
	if(!spill_in_out)	l->AddEntry(jff_residual_dEta[1][i][j],"JFF No JEC");

	l->SetLineColor(kWhite);
	l->SetTextSize(0.05);
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


      cout<<"ready for diff "<<i<<" "<<j<<endl;
      c_yields_eta_diff->cd(5*(i+1)-j);

      signal_dEta_diff[1][i][j] = (TH1D*)signal_dEta[0][i][j]->Clone((TString)("Diff_Preapp_" +Type_strs[1]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

      signal_dEta_diff[2][i][j] = (TH1D*)signal_dEta[2][i][j]->Clone((TString)("Diff_QM_" +Type_strs[2]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

      if(j<4)    signal_dEta_diff[3][i][j] = (TH1D*)signal_dEta[3][i][j]->Clone((TString)("Diff_Aug_" +Type_strs[3]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));


      signal_dEta_diff[1][i][j]->Add(signal_dEta[1][i][j],-1.);
      signal_dEta_diff[2][i][j]->Add(signal_dEta[1][i][j],-1.);

      signal_dEta_diff[1][i][j]->SetLineColor(kRed);
      signal_dEta_diff[1][i][j]->SetMarkerColor(kRed);

      signal_dEta_diff[2][i][j]->SetLineColor(kBlue);
      signal_dEta_diff[2][i][j]->SetMarkerColor(kBlue);

      signal_dEta_diff[1][i][j]->SetMaximum(2.5);
      signal_dEta_diff[1][i][j]->SetMinimum(-0.5);
      signal_dEta_diff[1][i][j]->Draw();
      signal_dEta_diff[2][i][j]->Draw("same");
      if(j<4)  signal_dEta_diff[3][i][j]->Draw("same");

      if(!spill_in_out){
      jff_residual_dEta[3][i][j]->SetMarkerStyle(10);
      jff_residual_dEta[3][i][j]->SetMarkerSize(1);
      jff_residual_dEta[3][i][j]->SetMarkerColor(kGreen);
      jff_residual_dEta[3][i][j]->SetLineColor(kGreen);
      jff_residual_dEta[3][i][j]->Draw("same");
      }
      zero->Draw("same");
     
      if(i==0&&j==4){
	TLegend *l = new TLegend(0.5,0.6,0.9,0.9);

	l->AddEntry(signal_dEta_diff[1][i][j],"Latest - Preapproval");
	l->AddEntry(signal_dEta_diff[2][i][j],"QM - Preapproval");
	if(!spill_in_out)	l->AddEntry( jff_residual_dEta[3][i][j],"JFF Change");
	l->SetLineColor(kWhite);
	l->SetTextSize(0.06);
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

  if(spill_in_out){
    c_yields_eta->SaveAs("Raw_Eta_Overlay_SpillIn_SpillOut.png");
    c_yields_eta_diff->SaveAs("Raw_Eta_Diffs_SpillIn_SpillOut.png");
  }else{
    c_yields_eta->SaveAs("Raw_Eta_Overlay.png");
    c_yields_eta_diff->SaveAs("Raw_Eta_Diffs.png");
  }
  
  return 0;
}
