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

Int_t data_mc_v2(){

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
 
  TH1D *sideband[12][nTrkPtBins][nCBins];
  TH1D *sideband_right[12][nTrkPtBins][nCBins];
 
  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%","pp"};

  float TrkPtBins[nTrkPtBins+1] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
  
  TString Type_strs[5] = {"Data","SubeNon0"};

  TFile *f_in[5];

  f_in[0] = new TFile("../me_correct/PbPb_Inclusive_Correlations.root");
  f_in[1] = new TFile("../me_correct/HydJet_RecoJet_GenTrack_NoSube0_Inclusive_Correlations.root");

   
  TString in_name, pTlabel,centlabel,Ajlabel;

  int lbin, rbin;

  double bin_width_phi, bin_width_eta, diff_max, diff_min, signal_min, signal_max, stacked_max, stacked_max_diff, stacked_min_diff, stacked_min, bc, err, ymax, ymin;
    
  c_yields_phi = new TCanvas("yields_phi","",10,10,3200,3200);
  c_yields_phi->Divide(4,8,0,0);

   
  c_yields_phi_sub = new TCanvas("yields_phi_sub","",10,10,3200,3200);
  c_yields_phi_sub->Divide(4,8,0,0);

  TF1 *fourier = new TF1("fourier","[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))");
 

  double v2_val[12][6][4][3];


  fourier->SetParName(0,"Bkg level");
  fourier->SetParName(1,"V_{1}");
  fourier->SetParName(2,"V_{2}");
  fourier->SetParName(3,"V_{3}");
 

  for(int i = 0; i<nTrkPtBins-1; i++){

    for(int j = 0; j<4; j++){

      for(int g = 0; g<2; g++){

	if(g==0)	result[g][i][j] = (TH2D*)f_in[g]->Get((TString)("Yield_PbPb_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Raw_Yield_" + Type_strs[g]+"_"+CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	else	result[g][i][j] = (TH2D*)f_in[g]->Get((TString)("Yield_Hydjet_" + CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]))->Clone((TString)("Raw_Yield_" + Type_strs[g]+"_"+CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]));

	if(i>3){
	  result[g][i][j]->Scale(1./4.);
	}else if(i==0){
	  result[g][i][j]->Scale(1/.3);
	}

	int lbin = result[g][i][j]->GetXaxis()->FindBin(-2.499);
	int rbin = result[g][i][j]->GetXaxis()->FindBin(-1.499);

	sideband[g][i][j] = (TH1D*)result[g][i][j]->ProjectionY((TString)("Proj_dPhi_" +Type_strs[g]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);


	lbin = result[g][i][j]->GetXaxis()->FindBin(1.499);
	rbin = result[g][i][j]->GetXaxis()->FindBin(2.499);

	sideband_right[g][i][j] = (TH1D*)result[g][i][j]->ProjectionY((TString)("Proj_dPhi_" +Type_strs[g]+"_"+ CBin_strs[j] + "_" + CBin_strs[j+1] + "_Pt100_Pt1000_" + TrkPtBin_strs[i] + "_" + TrkPtBin_strs[i+1]),lbin,rbin);


	sideband[g][i][j]->Add(sideband_right[g][i][j]);
	sideband[g][i][j]->SetMarkerStyle(10);
	sideband[g][i][j]->SetMarkerSize(1);

	sideband[g][i][j]->Rebin(5);

	float norm = 	sideband[g][i][j]->GetBinWidth(1);

	sideband[g][i][j]->Scale(1./norm);


      }
      
    }
    for(int j = 0; j<4; j++){
      c_yields_phi->cd(4*(i+1)-j);

      sideband[0][i][j]->SetLineColor(kBlack);
      sideband[0][i][j]->SetMarkerColor(kBlack);

      sideband[1][i][j]->SetLineColor(kBlue);
      sideband[1][i][j]->SetMarkerColor(kBlue);
  

      ymax = max( sideband[0][i][j]->GetMaximum(), sideband[1][i][j]->GetMaximum())+.1;
      ymin = min( sideband[0][i][j]->GetMinimum(), sideband[1][i][j]->GetMinimum())-.1;

      sideband[0][i][j]->SetMaximum(ymax);
      sideband[0][i][j]->SetMinimum(ymin);
     
      sideband[0][i][j]->GetXaxis()->SetRangeUser(-TMath::Pi()/2., TMath::Pi()/2.);
     
      if(j==3)  sideband[0][i][j]->GetYaxis()->SetLabelSize(0.1);
      else   sideband[0][i][j]->GetYaxis()->SetLabelSize(0.);

      sideband[0][i][j]->GetXaxis()->SetLabelSize(0.1);

      
      sideband[0][i][j]->Draw();
      sideband[0][i][j]->Fit("fourier","q");

      v2_val[0][i][j] = fourier->GetParameter(2);
      
      sideband[1][i][j]->Draw("same");

      sideband[1][i][j]->Fit("fourier","q");

      v2_val[1][i][j] = fourier->GetParameter(2);
    
      TLine *zero = new TLine(-TMath::Pi()/2.,0.,3.*TMath::Pi()/2.,0.);
      zero->SetLineStyle(2);
      zero->Draw("same");

      if(i==0&&j==3){
	TLegend *l = new TLegend(0.5,0.5,0.9,0.9);
	l->AddEntry(sideband[0][i][j],"Data");
	l->AddEntry(sideband[1][i][j],"Hydjet");
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
   
  c_yields_phi->SaveAs("Data_MC_v2.png");

  return 0;
}
