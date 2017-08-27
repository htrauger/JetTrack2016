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

Int_t jff_qg_eval(bool quark_gluon = kTRUE, bool rx_plane = kFALSE, bool spilled_in_out = kFALSE, bool is_jff = kTRUE){

  cout<<"quark gluon: "<<quark_gluon<<" rx plane: "<<rx_plane<<" spill in/out: "<<spilled_in_out<<" JFF? "<<is_jff<<endl;

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.25);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);

  TCanvas *c_yields_phi, *c_yields_phi_sub;

  const int nCBins = 5;
  const int nPtBins = 1;
  const int nTrkPtBins = 9;



  TH2D *result[12][nTrkPtBins][nCBins];
 
  TH1D *signal_dPhi[12][nTrkPtBins][nCBins];
  TH1D *signal_dPhi_sub[12][nTrkPtBins][nCBins];
 
  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%","pp"};

   float TrkPtBins[nTrkPtBins+1] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
  
  TString Type_strs[5] = {"Reco","Gen","Corr","Quark","Gluon"};

  TFile *f_in[5];

  f_in[0] = new TFile("../me_correct/HydJet_RecoJet_GenTrack_Sube0_Inclusive_Correlations.root");
  f_in[1] = new TFile("../me_correct/HydJet_GenJet_GenTrack_Sube0_Inclusive_Correlations.root");
  f_in[2] = new TFile("../me_correct/HydJet_GenJet_GenTrack_Sube0_Inclusive_Correlations.root");

  f_in[3] = new TFile("../me_correct/HydJet_RecoJet_GenTrack_Sube0_Inclusive_Correlations_QuarkOnly.root");
  f_in[4] = new TFile("../me_correct/HydJet_RecoJet_GenTrack_Sube0_Inclusive_Correlations_GluonOnly.root");
  
 
  TString in_name, pTlabel,centlabel,Ajlabel;

  int lbin, rbin;

  double bin_width_phi, bin_width_eta, diff_max, diff_min, signal_min, signal_max, stacked_max, stacked_max_diff, stacked_min_diff, stacked_min, bc, err, ymax, ymin;
    
  c_yields_phi = new TCanvas("yields_phi","",10,10,3200,3200);
  c_yields_phi->Divide(4,8,0,0);

   
  c_yields_phi_sub = new TCanvas("yields_phi_sub","",10,10,3200,3200);
  c_yields_phi_sub->Divide(4,8,0,0);

  TF1 *zyam = new TF1("zyam","[0]+x-x",-TMath::Pi()/2.,3.*TMath::Pi()/2.);


  for(int i = 0; i<nTrkPtBins-1; i++){

    for(int j = 0; j<4; j++){

      for(int g = 0; g<5; g++){

	result[g][i][j] = (TH2D*)f_in[g]->Get((TString)("Yield_BkgSub_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Raw_Yield_" + Type_strs[g]+"_"+CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	if(i>3){
	  result[g][i][j]->Scale(1./4.);
	}else if(i==0){
	  result[g][i][j]->Scale(1/.3);
	}

	int lbin = result[g][i][j]->GetXaxis()->FindBin(-0.99);
	int rbin = result[g][i][j]->GetXaxis()->FindBin(.99);

	signal_dPhi[g][i][j] = (TH1D*)result[g][i][j]->ProjectionY((TString)("Proj_dPhi_" +Type_strs[g]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);

	signal_dPhi[g][i][j]->SetMarkerStyle(10);
	signal_dPhi[g][i][j]->SetMarkerSize(1);

	signal_dPhi[g][i][j]->Rebin(5);

	float norm = 	signal_dPhi[g][i][j]->GetBinWidth(1);

	signal_dPhi[g][i][j]->Scale(1./norm);

	float bkg_level =( signal_dPhi[g][i][j]->GetBinContent(signal_dPhi[g][i][j]->FindBin(-TMath::Pi()/2.))+signal_dPhi[g][i][j]->GetBinContent(signal_dPhi[g][i][j]->FindBin(TMath::Pi()/2.)))/2.;

	zyam->SetParameter(0,bkg_level);

	
	signal_dPhi[g][i][j]->Add(zyam,-1.);

	signal_dPhi_sub[g][i][j] = (TH1D*)signal_dPhi[g][i][j]->Clone((TString)("Proj_dPhi_SubTracked" +Type_strs[g]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));
	

	for(int k = 0; k< signal_dPhi_sub[g][i][j]->GetNbinsX()/2+1; k++){
	  bc = 	signal_dPhi_sub[g][i][j]->GetBinContent(k)-signal_dPhi[g][i][j]->GetBinContent(k+signal_dPhi_sub[g][i][j]->GetNbinsX()/2);
	  err = TMath::Sqrt(signal_dPhi[g][i][j]->GetBinError(k+signal_dPhi_sub[g][i][j]->GetNbinsX()/2)*signal_dPhi[g][i][j]->GetBinError(k+signal_dPhi_sub[g][i][j]->GetNbinsX()/2)+signal_dPhi_sub[g][i][j]->GetBinError(k)*signal_dPhi_sub[g][i][j]->GetBinError(k));
	  
	  signal_dPhi_sub[g][i][j]->SetBinContent(k,bc);
	  signal_dPhi_sub[g][i][j]->SetBinError(k,err);

	}

	/*

	lbin = result[g][i][j]->GetYaxis()->FindBin(-.99);
	rbin = result[g][i][j]->GetYaxis()->FindBin(.99);

	signal_dEta[g][i][j] = (TH1D*)result[g][i][j]->ProjectionX((TString)("Proj_dEta_" +Type_strs[g]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);

	signal_dEta[g][i][j]->SetMarkerStyle(10);
	signal_dEta[g][i][j]->SetMarkerSize(1);

	signal_dEta[g][i][j]->Rebin(5);

	norm = 	signal_dEta[g][i][j]->GetBinWidth(1);

	signal_dEta[g][i][j]->Scale(1./norm);

	*/
      }
      
    }

    
      for(int j = 0; j<4; j++){

	signal_dPhi_sub[2][i][j] = (TH1D*)signal_dPhi_sub[0][i][j]->Clone(Form("Corr%d%d",i,j));
	signal_dPhi_sub[2][i][j]->Add(signal_dPhi_sub[1][i][j],-1.);

      c_yields_phi->cd(4*(i+1)-j);

      signal_dPhi[0][i][j]->SetLineColor(kBlack);
      signal_dPhi[0][i][j]->SetMarkerColor(kBlack);

      signal_dPhi[1][i][j]->SetLineColor(kBlue);
      signal_dPhi[1][i][j]->SetMarkerColor(kBlue);

      signal_dPhi[2][i][j]->SetLineColor(kRed);
      signal_dPhi[2][i][j]->SetMarkerColor(kRed);

      signal_dPhi[3][i][j]->SetLineColor(kBlue);
      signal_dPhi[3][i][j]->SetMarkerColor(kBlue);

      signal_dPhi[4][i][j]->SetLineColor(kRed);
      signal_dPhi[4][i][j]->SetMarkerColor(kRed);

      signal_dPhi_sub[0][i][j]->SetLineColor(kBlack);
      signal_dPhi_sub[0][i][j]->SetMarkerColor(kBlack);

      signal_dPhi_sub[1][i][j]->SetLineColor(kBlue);
      signal_dPhi_sub[1][i][j]->SetMarkerColor(kBlue);

      signal_dPhi_sub[2][i][j]->SetLineColor(kRed);
      signal_dPhi_sub[2][i][j]->SetMarkerColor(kRed);

      signal_dPhi_sub[3][i][j]->SetLineColor(kBlue);
      signal_dPhi_sub[3][i][j]->SetMarkerColor(kBlue);

      signal_dPhi_sub[4][i][j]->SetLineColor(kRed);
      signal_dPhi_sub[4][i][j]->SetMarkerColor(kRed);

   
      switch(i){
      case 0:
	ymax = 15.;
	ymin = -1.;
       	break;
      case 1:
	ymax = 6.;
	ymin = -1.;
	break; 
      case 2:
	ymax = 4.;
	ymin = -1.;
	break; 
      case 3:
	ymax = 4.;
	ymin = -.5;
	break;
      default:
	ymax = 0.4;
	ymin = -0.4;
	break;
      }

      signal_dPhi[0][i][j]->SetMaximum(ymax);
      signal_dPhi[0][i][j]->SetMinimum(ymin);
     
      signal_dPhi[0][i][j]->GetXaxis()->SetRangeUser(-TMath::Pi()/2., 3.*TMath::Pi()/2.);
     
      if(j==3)  signal_dPhi[0][i][j]->GetYaxis()->SetLabelSize(0.1);
      else   signal_dPhi[0][i][j]->GetYaxis()->SetLabelSize(0.);

      signal_dPhi[0][i][j]->GetXaxis()->SetLabelSize(0.1);

      
      signal_dPhi[0][i][j]->Draw();
      if(rx_plane||spilled_in_out)     signal_dPhi[1][i][j]->Draw("same");
      if(rx_plane||spilled_in_out)  signal_dPhi[2][i][j]->Draw("same");
      if(quark_gluon)     signal_dPhi[3][i][j]->Draw("same");
      if(quark_gluon)  signal_dPhi[4][i][j]->Draw("same");
      //signal_dPhi[3][i][j]->Draw("same");

      TLine *zero = new TLine(-TMath::Pi()/2.,0.,3.*TMath::Pi()/2.,0.);
      zero->SetLineStyle(2);
      zero->Draw("same");

      if(i==0&&j==3){
	TLegend *l = new TLegend(0.5,0.5,0.9,0.9);
	l->AddEntry(signal_dPhi[0][i][j],"Nominal");
	if(rx_plane)	l->AddEntry(signal_dPhi[1][i][j],"In Rx-Plane");
	if(rx_plane)	l->AddEntry(signal_dPhi[2][i][j],"Out Rx-Plane");
	if(spilled_in_out)	l->AddEntry(signal_dPhi[1][i][j],"Spilled Out");
	if(spilled_in_out)	l->AddEntry(signal_dPhi[2][i][j],"Spilled In");
	if(quark_gluon)	l->AddEntry(signal_dPhi[3][i][j],"Quark Jets");
	if(quark_gluon)	l->AddEntry(signal_dPhi[4][i][j],"Gluon Jets");
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


      c_yields_phi_sub->cd(4*(i+1)-j);

      if(is_jff){
        signal_dPhi_sub[0][i][j]->SetMaximum(ymax);
	signal_dPhi_sub[0][i][j]->SetMinimum(ymin);
      }else{
	signal_dPhi_sub[0][i][j]->SetMaximum(2.);
	signal_dPhi_sub[0][i][j]->SetMinimum(-.5);

      }
      signal_dPhi_sub[0][i][j]->GetXaxis()->SetRangeUser(-TMath::Pi()/2., TMath::Pi()/2.);
     
      if(j==3)  signal_dPhi_sub[0][i][j]->GetYaxis()->SetLabelSize(0.1);
      else   signal_dPhi_sub[0][i][j]->GetYaxis()->SetLabelSize(0.);

      signal_dPhi_sub[0][i][j]->GetXaxis()->SetLabelSize(0.1);

      
      signal_dPhi_sub[0][i][j]->Draw();
      if(rx_plane||spilled_in_out)     signal_dPhi_sub[1][i][j]->Draw("same");
      if(rx_plane||spilled_in_out)  signal_dPhi_sub[2][i][j]->Draw("same");
      if(quark_gluon)     signal_dPhi_sub[3][i][j]->Draw("same");
      if(quark_gluon)  signal_dPhi_sub[4][i][j]->Draw("same");
      //signal_dPhi_sub[3][i][j]->Draw("same");

      zero->Draw("same");

      if(i==0&&j==3){
	TLegend *l = new TLegend(0.5,0.5,0.9,0.9);
	l->AddEntry(signal_dPhi_sub[0][i][j],"Nominal");
	if(rx_plane)	l->AddEntry(signal_dPhi_sub[1][i][j],"In Rx-Plane");
	if(rx_plane)	l->AddEntry(signal_dPhi_sub[2][i][j],"Out Rx-Plane");
	if(quark_gluon)	l->AddEntry(signal_dPhi_sub[3][i][j],"Quark Jets");
	if(quark_gluon)	l->AddEntry(signal_dPhi_sub[4][i][j],"Gluon Jets");
	if(spilled_in_out)	l->AddEntry(signal_dPhi[1][i][j],"Spilled Out");
	if(spilled_in_out)	l->AddEntry(signal_dPhi[2][i][j],"Spilled In");
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

    int lbin =    signal_dPhi_sub[1][i][0]->FindBin(-1.);
    int rbin =    signal_dPhi_sub[1][i][0]->FindBin(1.);

    //    if(rx_plane) cout<<TrkPtBin_strs[i]<<" "<<signal_dPhi_sub[0][i][0]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[1][i][0]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[2][i][0]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[0][i][1]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[1][i][1]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[2][i][1]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[0][i][2]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[1][i][2]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[2][i][2]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[0][i][3]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[1][i][3]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[2][i][3]->Integral(lbin,rbin)<<" "<<endl;

    if(quark_gluon) cout<<TrkPtBin_strs[i]<<" "<<signal_dPhi_sub[2][i][0]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[3][i][0]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[4][i][0]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[2][i][1]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[3][i][1]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[4][i][1]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[2][i][2]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[3][i][2]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[4][i][2]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[2][i][3]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[3][i][3]->Integral(lbin,rbin)<<" "<<signal_dPhi_sub[4][i][3]->Integral(lbin,rbin)<<" "<<endl;


  }

  if(rx_plane)  c_yields_phi->SaveAs("Raw_Spill_Over_By_RxPlane.png");
  if(spilled_in_out)  c_yields_phi->SaveAs("Raw_Spill_Over_By_SpilledJets.png"); 
  if(quark_gluon&&is_jff)  c_yields_phi->SaveAs("Raw_Sube0_QuarkGluon.png");
  if(quark_gluon&&is_jff)  c_yields_phi->SaveAs("Raw_Sube0_QuarkGluon.png");
 

  if(rx_plane)  c_yields_phi_sub->SaveAs("Raw_Spill_Over_By_RxPlane_BkgSub.png");
  if(spilled_in_out)  c_yields_phi_sub->SaveAs("Raw_Spill_Over_By_SpilledJets_BkgSub.png");
  if(quark_gluon&&!is_jff)  c_yields_phi_sub->SaveAs("Raw_Spill_Over_QuarkGluon_BkgSub.png");
  if(quark_gluon&&is_jff)  c_yields_phi_sub->SaveAs("Raw_Sube0_QuarkGluon_BkgSub.png");
  

  return 0;
}
