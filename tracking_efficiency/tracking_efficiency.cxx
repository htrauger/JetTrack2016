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


#include <iostream>
#include <vector>
#include <fstream>

#include "../JetTrack2016_functions.h"


Int_t tracking_efficiency(bool is_pythia = kFALSE){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

 
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 9;

  float CBins[nCBins+1] = {0, 10, 30, 50, 100};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%", "Cent. 30-50%","Cent. 50-100%"};

   float TrkPtBins[nTrkPtBins+1] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt0p7","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt999" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
 

  TH1D *pt_gen[nCBins][nTrkPtBins];
  TH1D *pt_reco[nCBins][nTrkPtBins];
  TH1D *pt_corr[nCBins][nTrkPtBins];

  TH1D *pt_reco_ratio[nCBins][nTrkPtBins];
  TH1D *pt_corr_ratio[nCBins][nTrkPtBins];


  TH1D *eta_gen[nCBins][nTrkPtBins];
  TH1D *eta_reco[nCBins][nTrkPtBins];
  TH1D *eta_corr[nCBins][nTrkPtBins];

  TH1D *eta_reco_ratio[nCBins][nTrkPtBins];
  TH1D *eta_corr_ratio[nCBins][nTrkPtBins];



  TH1D *phi_gen[nCBins][nTrkPtBins];
  TH1D *phi_reco[nCBins][nTrkPtBins];
  TH1D *phi_corr[nCBins][nTrkPtBins];

  TH1D *phi_reco_ratio[nCBins][nTrkPtBins];
  TH1D *phi_corr_ratio[nCBins][nTrkPtBins];

  TCanvas *pt_canvas;

  if(is_pythia){
    pt_canvas= new TCanvas("pt_canvas","",10,10,500,800);
    pt_canvas->Divide(1,2,0,0);

  }else{
    pt_canvas= new TCanvas("pt_canvas","",10,10,1500,800);
    pt_canvas->Divide(4,2,0,0);


  }

  TCanvas *eta_canvas[nTrkPtBins];
  TCanvas *phi_canvas[nTrkPtBins];

  for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){
    TString etaname = "eta_canvas"; etaname+=ibin3;
   
    TString phiname = "phi_canvas"; phiname+=ibin3;
    
    if(is_pythia){
      eta_canvas[ibin3] = new TCanvas(etaname,"",10,10,500,800);
      eta_canvas[ibin3]->Divide(1,2,0,0);
      phi_canvas[ibin3] = new TCanvas(phiname,"",10,10,500,800);
      phi_canvas[ibin3]->Divide(1,2,0,0);
   

    }else{
      eta_canvas[ibin3] = new TCanvas(etaname,"",10,10,1500,800);
      eta_canvas[ibin3]->Divide(4,2,0,0);
      phi_canvas[ibin3] = new TCanvas(phiname,"",10,10,1500,800);
      phi_canvas[ibin3]->Divide(4,2,0,0);
   

    }

  }


  TFile *f_Gen, *f_Reco;
  
  if(is_pythia){
    f_Gen  = new TFile("../mc_raw_correlations/Pythia_GenJet_GenTrack_Aug23.root");
     f_Reco = new TFile("../mc_raw_correlations/Pythia_GenJet_RecoTrack_Aug23.root");

  }else{
    f_Gen  = new TFile("../mc_raw_correlations/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_GenGenReduced_fineBin.root");
    f_Reco = new TFile("../mc_raw_correlations/PbPb_5TeV_MC_PythiaHydjet_MixHistos_GenRecoReduced_refpt_newJetTrackCorrections_Merged_fineBin.root");
  }
  

  cout<<"got file"<<endl;
  float norm;
  
  for(int ibin = 0; ibin<nCBins; ibin++){

    if(is_pythia&&ibin>0)continue;

    for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){

         
      pt_gen[ibin][ibin3]=(TH1D*)f_Gen->Get((TString)("GenJet_GenTrack_TrkPt"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone((TString)("GenJet_GenTrack_TrkPt"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

      cout<<"got gen spectrum"<<endl;

      pt_reco[ibin][ibin3]=(TH1D*)f_Reco->Get((TString)("GenJet_RecoTrack_TrkPt"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone((TString)("GenJet_RecoTrack_TrkPt"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));
      
      
      pt_corr[ibin][ibin3]=(TH1D*)f_Reco->Get((TString)("GenJet_RecoTrack_TrkPt_weighted"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone((TString)("GenJet_RecoTrack_TrkPt_weighted"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));


      cout<<"got all"<<endl;

      pt_gen[ibin][ibin3]->Rebin(25);
      pt_reco[ibin][ibin3]->Rebin(25);
      pt_corr[ibin][ibin3]->Rebin(25);
            /*
      norm = pt_gen[ibin][ibin3]->Integral();
      pt_gen[ibin][ibin3]->Scale(1./norm);
     
      norm = pt_reco[ibin][ibin3]->Integral();
      pt_reco[ibin][ibin3]->Scale(1./norm);

      norm = pt_corr[ibin][ibin3]->Integral();
      pt_corr[ibin][ibin3]->Scale(1./norm);
     
      */

      pt_gen[ibin][ibin3]->SetLineColor(kBlack);
      pt_gen[ibin][ibin3]->SetMarkerColor(kBlack);

      pt_gen[ibin][ibin3]->SetMarkerSize(2);
      pt_gen[ibin][ibin3]->SetMarkerStyle(21);
   
      pt_reco[ibin][ibin3]->SetMarkerSize(1);
      pt_reco[ibin][ibin3]->SetMarkerStyle(20);
   
      pt_corr[ibin][ibin3]->SetMarkerSize(1);
      pt_corr[ibin][ibin3]->SetMarkerStyle(20);
    
      if(is_pythia)pt_canvas->cd(1);
      else pt_canvas->cd(4-ibin);
	
      if(ibin3==0){


	gPad->SetLogy();
	//gPad->SetLogx();


	pt_gen[ibin][ibin3]->SetAxisRange(.7,19.5);
	if(is_pythia){
	  pt_gen[ibin][ibin3]->SetMinimum(5000.);
	  pt_gen[ibin][ibin3]->SetMaximum(1e8);
	}else{
	  pt_gen[ibin][ibin3]->SetMinimum(1.);
	  pt_gen[ibin][ibin3]->SetMaximum(1e7);

	}
	pt_gen[ibin][ibin3]->GetYaxis()->SetLabelSize(0.);
	if(ibin==3||is_pythia){ pt_gen[ibin][ibin3]->GetYaxis()->SetLabelSize(0.05); }


	pt_gen[ibin][ibin3]->Draw();
      }else{
	pt_gen[ibin][ibin3]->Draw("same");
      }	

      pt_reco[ibin][ibin3]->SetLineColor(kBlue);
      pt_reco[ibin][ibin3]->SetMarkerColor(kBlue);
      pt_reco[ibin][ibin3]->Draw("same");
    
      pt_corr[ibin][ibin3]->SetLineColor(kRed);
      pt_corr[ibin][ibin3]->SetMarkerColor(kRed);
      pt_corr[ibin][ibin3]->Draw("same");

      if(ibin==3||is_pythia){
	TLegend *lpt = new TLegend(0.15,0.75,0.6,0.95);
	lpt->AddEntry(pt_gen[ibin][ibin3],"Gen Tracks");
	lpt->AddEntry(pt_reco[ibin][ibin3], "Reco Tracks");
	lpt->AddEntry(pt_corr[ibin][ibin3], "Corr. Tracks");
	lpt->SetLineColor(kWhite);
	lpt->SetTextSize(0.05);
	lpt->Draw();
      }

     
      TLatex *cent = new TLatex(0.6,0.9,CBin_labels[ibin]);
      cent->SetNDC();
      if(!is_pythia) cent->Draw();
         
      pt_reco_ratio[ibin][ibin3]=(TH1D*)pt_reco[ibin][ibin3]->Clone((TString)("GenJet_RecoTrack_TrkPt_Ratio_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

      pt_corr_ratio[ibin][ibin3]=(TH1D*)pt_corr[ibin][ibin3]->Clone((TString)("GenJet_RecoTrack_TrkPt_Weighted_Ratio_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

      pt_reco_ratio[ibin][ibin3]->Divide(pt_gen[ibin][ibin3]);
      pt_corr_ratio[ibin][ibin3]->Divide(pt_gen[ibin][ibin3]);
      
      if(is_pythia)pt_canvas->cd(2);
      else pt_canvas->cd(8-ibin);
      if(ibin3==0){
	
	pt_reco_ratio[ibin][ibin3]->GetYaxis()->SetLabelSize(0.);

	if(ibin==3||is_pythia){	pt_reco_ratio[ibin][ibin3]->GetYaxis()->SetLabelSize(0.05); }
	pt_reco_ratio[ibin][ibin3]->GetXaxis()->SetLabelSize(0.05);
	pt_reco_ratio[ibin][ibin3]->GetXaxis()->SetTitleSize(0.05);
	pt_reco_ratio[ibin][ibin3]->GetXaxis()->SetTitle("p_{T}^{assoc.} (GeV/c)"); 
	pt_reco_ratio[ibin][ibin3]->GetXaxis()->SetTitleOffset(0.8);
	pt_reco_ratio[ibin][ibin3]->GetXaxis()->CenterTitle();
	

	//	gPad->SetLogx();

	pt_reco_ratio[ibin][ibin3]->SetMinimum(0.);
	pt_reco_ratio[ibin][ibin3]->SetMaximum(2.);
	pt_reco_ratio[ibin][ibin3]->SetAxisRange(.5,19.5);

	pt_reco_ratio[ibin][ibin3]->Draw();
	pt_corr_ratio[ibin][ibin3]->Draw("same");

	TLine *line_pT = new TLine(1.,1.,19.5,1.);
	line_pT->SetLineStyle(2);
	line_pT->Draw();

	TLine *line_pT_down = new TLine(.7,.95,19.5,.95);
	line_pT_down->SetLineStyle(2);
	line_pT_down->Draw();

	TLine *line_pT_up = new TLine(.7,1.05,19.5,1.05);
	line_pT_up->SetLineStyle(2);
	line_pT_up->Draw();


      }else{
	pt_reco_ratio[ibin][ibin3]->Draw("same");
	pt_corr_ratio[ibin][ibin3]->Draw("same");

	if(ibin3==8)	pt_corr_ratio[ibin][ibin3]->Fit("pol0");
      }
  
      if(ibin==3||is_pythia){
	TLegend *lpt_ratio = new TLegend(0.15,0.8,0.6,0.95);
	lpt_ratio->AddEntry(pt_reco_ratio[ibin][ibin3], "Reco/Gen");
	lpt_ratio->AddEntry(pt_corr_ratio[ibin][ibin3], "Corr/Gen");
	lpt_ratio->SetLineColor(kWhite);
	lpt_ratio->SetTextSize(0.05);
	lpt_ratio->Draw();
      }

      cout<<"done pt"<<endl;

      //------
      // Eta
      //------

   
      eta_gen[ibin][ibin3]=(TH1D*)f_Gen->Get((TString)("GenJet_GenTrack_TrkEta"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone((TString)("GenJet_GenTrack_TrkPt"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

      eta_reco[ibin][ibin3]=(TH1D*)f_Reco->Get((TString)("GenJet_RecoTrack_TrkEta"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone((TString)("GenJet_RecoTrack_TrkPt"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));
      
      
      eta_corr[ibin][ibin3]=(TH1D*)f_Reco->Get((TString)("GenJet_RecoTrack_TrkEta_weighted"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone((TString)("GenJet_RecoTrack_TrkEta_weighted"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));



      eta_gen[ibin][ibin3]->Rebin(2);
      eta_reco[ibin][ibin3]->Rebin(2);
      eta_corr[ibin][ibin3]->Rebin(2);
  
      eta_gen[ibin][ibin3]->SetLineColor(kBlack);
      eta_gen[ibin][ibin3]->SetMarkerColor(kBlack);

      eta_gen[ibin][ibin3]->SetMarkerSize(2);
      eta_gen[ibin][ibin3]->SetMarkerStyle(21);
   
      eta_reco[ibin][ibin3]->SetMarkerSize(1);
      eta_reco[ibin][ibin3]->SetMarkerStyle(20);
   
      eta_corr[ibin][ibin3]->SetMarkerSize(1);
      eta_corr[ibin][ibin3]->SetMarkerStyle(20);
    
      if(is_pythia) eta_canvas[ibin3]->cd(1);
      else eta_canvas[ibin3]->cd(4-ibin);
	
    
      eta_gen[ibin][ibin3]->SetAxisRange(-2.4,2.4);

      eta_gen[ibin][ibin3]->SetMinimum(0.);
      eta_gen[ibin][ibin3]->SetMaximum(1.5* eta_gen[ibin][ibin3]->GetMaximum());

      eta_gen[ibin][ibin3]->GetYaxis()->SetLabelSize(0.);
      eta_gen[ibin][ibin3]->GetYaxis()->SetLabelSize(0.05); 


      eta_gen[ibin][ibin3]->Draw();
    

      eta_reco[ibin][ibin3]->SetLineColor(kBlue);
      eta_reco[ibin][ibin3]->SetMarkerColor(kBlue);
      eta_reco[ibin][ibin3]->Draw("same");
    
      eta_corr[ibin][ibin3]->SetLineColor(kRed);
      eta_corr[ibin][ibin3]->SetMarkerColor(kRed);
      eta_corr[ibin][ibin3]->Draw("same");

      if(ibin==3||is_pythia){
	TLegend *lpt = new TLegend(0.15,0.75,0.6,0.95);
	lpt->AddEntry(eta_gen[ibin][ibin3],"Gen Tracks");
	lpt->AddEntry(eta_reco[ibin][ibin3], "Reco Tracks");
	lpt->AddEntry(eta_corr[ibin][ibin3], "Corr. Tracks");
	lpt->SetLineColor(kWhite);
	lpt->SetTextSize(0.05);
	lpt->Draw();
      }

      
      if(!is_pythia) cent->Draw();

      TLatex *trkpt = new TLatex(0.6,0.85,TrkPtBin_labels[ibin3]);
      trkpt->SetNDC();
      trkpt->Draw();
  
  
      eta_reco_ratio[ibin][ibin3]=(TH1D*)eta_reco[ibin][ibin3]->Clone((TString)("GenJet_RecoTrack_TrkEta_Ratio_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

      eta_corr_ratio[ibin][ibin3]=(TH1D*)eta_corr[ibin][ibin3]->Clone((TString)("GenJet_RecoTrack_TrkEta_Weighted_Ratio_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

      eta_reco_ratio[ibin][ibin3]->Divide(eta_gen[ibin][ibin3]);
      eta_corr_ratio[ibin][ibin3]->Divide(eta_gen[ibin][ibin3]);
      
      if(is_pythia)  eta_canvas[ibin3]->cd(2);
      else eta_canvas[ibin3]->cd(8-ibin);
    	
      //   eta_reco_ratio[ibin][ibin3]->GetYaxis()->SetLabelSize(0.);

      eta_reco_ratio[ibin][ibin3]->GetYaxis()->SetLabelSize(0.05); 
      eta_reco_ratio[ibin][ibin3]->GetXaxis()->SetLabelSize(0.05);
      eta_reco_ratio[ibin][ibin3]->GetXaxis()->SetTitleSize(0.05);
      eta_reco_ratio[ibin][ibin3]->GetXaxis()->SetTitle("#eta_{track}"); 
      eta_reco_ratio[ibin][ibin3]->GetXaxis()->SetTitleOffset(0.8);
      eta_reco_ratio[ibin][ibin3]->GetXaxis()->CenterTitle();
	

      //	gPad->SetLogx();

      eta_reco_ratio[ibin][ibin3]->SetMinimum(0.);
      eta_reco_ratio[ibin][ibin3]->SetMaximum(2.);
      eta_reco_ratio[ibin][ibin3]->SetAxisRange(-2.4,2.4);

      eta_reco_ratio[ibin][ibin3]->Draw();
      eta_corr_ratio[ibin][ibin3]->Draw("same");
      eta_corr_ratio[ibin][ibin3]->Fit("pol0");

      TLine *line_eta = new TLine(-2.5,1.,2.5,1.);
      line_eta->SetLineStyle(2);
      line_eta->Draw();

      TLine *line_eta_up = new TLine(-2.5,.95,2.5,.95);
      line_eta_up->SetLineStyle(2);
      line_eta_up->Draw();

      TLine *line_eta_down = new TLine(-2.5,1.05,2.5,1.05);
      line_eta_down->SetLineStyle(2);
      line_eta_down->Draw();

      if(ibin==3||is_pythia){
	TLegend *leta_ratio = new TLegend(0.15,0.8,0.6,0.95);
	leta_ratio->AddEntry(eta_reco_ratio[ibin][ibin3], "Reco/Gen");
	leta_ratio->AddEntry(eta_corr_ratio[ibin][ibin3], "Corr/Gen");
	leta_ratio->SetLineColor(kWhite);
	leta_ratio->SetTextSize(0.05);
	leta_ratio->Draw();
      }

      cout<<"done eta"<<endl;


      //------
      // Phi
      //------

   
      phi_gen[ibin][ibin3]=(TH1D*)f_Gen->Get((TString)("GenJet_GenTrack_TrkPhi"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone((TString)("GenJet_GenTrack_TrkPt"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

      phi_reco[ibin][ibin3]=(TH1D*)f_Reco->Get((TString)("GenJet_RecoTrack_TrkPhi"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone((TString)("GenJet_RecoTrack_TrkPt"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));
      
      
      phi_corr[ibin][ibin3]=(TH1D*)f_Reco->Get((TString)("GenJet_RecoTrack_TrkPhi_weighted"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone((TString)("GenJet_RecoTrack_TrkPhi_weighted"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));



      phi_gen[ibin][ibin3]->Rebin(4);
      phi_reco[ibin][ibin3]->Rebin(4);
      phi_corr[ibin][ibin3]->Rebin(4);
  
      phi_gen[ibin][ibin3]->SetLineColor(kBlack);
      phi_gen[ibin][ibin3]->SetMarkerColor(kBlack);

      phi_gen[ibin][ibin3]->SetMarkerSize(2);
      phi_gen[ibin][ibin3]->SetMarkerStyle(21);
   
      phi_reco[ibin][ibin3]->SetMarkerSize(1);
      phi_reco[ibin][ibin3]->SetMarkerStyle(20);
   
      phi_corr[ibin][ibin3]->SetMarkerSize(1);
      phi_corr[ibin][ibin3]->SetMarkerStyle(20);
    
      if(is_pythia)    phi_canvas[ibin3]->cd(1);
      else phi_canvas[ibin3]->cd(4-ibin);
	
    
      phi_gen[ibin][ibin3]->SetAxisRange(-TMath::Pi()/2.,TMath::Pi()/2.);

      phi_gen[ibin][ibin3]->SetMinimum(0.);
      phi_gen[ibin][ibin3]->SetMaximum(1.5* phi_gen[ibin][ibin3]->GetMaximum());

      //  phi_gen[ibin][ibin3]->GetYaxis()->SetLabelSize(0.);
      phi_gen[ibin][ibin3]->GetYaxis()->SetLabelSize(0.05); 


      phi_gen[ibin][ibin3]->Draw();
    

      phi_reco[ibin][ibin3]->SetLineColor(kBlue);
      phi_reco[ibin][ibin3]->SetMarkerColor(kBlue);
      phi_reco[ibin][ibin3]->Draw("same");
    
      phi_corr[ibin][ibin3]->SetLineColor(kRed);
      phi_corr[ibin][ibin3]->SetMarkerColor(kRed);
      phi_corr[ibin][ibin3]->Draw("same");

      if(ibin==3||is_pythia){
	TLegend *lpt = new TLegend(0.15,0.75,0.6,0.95);
	lpt->AddEntry(phi_gen[ibin][ibin3],"Gen Tracks");
	lpt->AddEntry(phi_reco[ibin][ibin3], "Reco Tracks");
	lpt->AddEntry(phi_corr[ibin][ibin3], "Corr. Tracks");
	lpt->SetLineColor(kWhite);
	lpt->SetTextSize(0.05);
	lpt->Draw();
      }

      
      if(!is_pythia) cent->Draw();
      trkpt->Draw();
  
      phi_reco_ratio[ibin][ibin3]=(TH1D*)phi_reco[ibin][ibin3]->Clone((TString)("GenJet_RecoTrack_TrkPhi_Ratio_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

      phi_corr_ratio[ibin][ibin3]=(TH1D*)phi_corr[ibin][ibin3]->Clone((TString)("GenJet_RecoTrack_TrkPhi_Weighted_Ratio_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

      phi_reco_ratio[ibin][ibin3]->Divide(phi_gen[ibin][ibin3]);
      phi_corr_ratio[ibin][ibin3]->Divide(phi_gen[ibin][ibin3]);
     
      if(is_pythia)   phi_canvas[ibin3]->cd(2);
      else  phi_canvas[ibin3]->cd(8-ibin);
    	
      phi_reco_ratio[ibin][ibin3]->GetYaxis()->SetLabelSize(0.);

      phi_reco_ratio[ibin][ibin3]->GetYaxis()->SetLabelSize(0.05); 
      phi_reco_ratio[ibin][ibin3]->GetXaxis()->SetLabelSize(0.05);
      phi_reco_ratio[ibin][ibin3]->GetXaxis()->SetTitleSize(0.05);
      phi_reco_ratio[ibin][ibin3]->GetXaxis()->SetTitleOffset(0.8);
      phi_reco_ratio[ibin][ibin3]->GetXaxis()->SetTitle("#phi_{track}"); 
      phi_reco_ratio[ibin][ibin3]->GetXaxis()->CenterTitle();
	

      //	gPad->SetLogx();

      phi_reco_ratio[ibin][ibin3]->SetMinimum(0.);
      phi_reco_ratio[ibin][ibin3]->SetMaximum(2.);
      phi_reco_ratio[ibin][ibin3]->SetAxisRange(-TMath::Pi()/2.,TMath::Pi()/2.);

      phi_reco_ratio[ibin][ibin3]->Draw();
      phi_corr_ratio[ibin][ibin3]->Draw("same");
      phi_corr_ratio[ibin][ibin3]->Fit("pol0");

      TLine *line_phi = new TLine(-TMath::Pi()/2.,1.,TMath::Pi()/2.,1.);
      line_phi->SetLineStyle(2);
      line_phi->Draw();


      TLine *line_phi_up = new TLine(-TMath::Pi()/2.,1.05,TMath::Pi()/2.,1.05);
      line_phi_up->SetLineStyle(2);
      line_phi_up->Draw();

   TLine *line_phi_down = new TLine(-TMath::Pi()/2.,.95,TMath::Pi()/2.,.95);
      line_phi_down->SetLineStyle(2);
      line_phi_down->Draw();

      if(ibin==3||is_pythia){
	TLegend *lphi_ratio = new TLegend(0.15,0.8,0.6,0.95);
	lphi_ratio->AddEntry(phi_reco_ratio[ibin][ibin3], "Reco/Gen");
	lphi_ratio->AddEntry(phi_corr_ratio[ibin][ibin3], "Corr/Gen");
	lphi_ratio->SetLineColor(kWhite);
	lphi_ratio->SetTextSize(0.05);
	lphi_ratio->Draw();
      }
      

      cout<<"done phi"<<endl;
    }
  }

  cout<<"done all"<<endl;

  if(is_pythia){

    pt_canvas->SaveAs("TrackingEfficiencyPtPythia.png");
    pt_canvas->SaveAs("TrackingEfficiencyPtPythia.pdf");


  }else{
    pt_canvas->SaveAs("TrackingEfficiencyPt.png");
    pt_canvas->SaveAs("TrackingEfficiencyPt.pdf");
  }

  for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){
    

    if(is_pythia){
      eta_canvas[ibin3]->SaveAs((TString)("TrackingEfficiencyEtaPythia_"+TrkPtBin_strs[ibin3]+TrkPtBin_strs[ibin3+1]+".png"));
      eta_canvas[ibin3]->SaveAs((TString)("TrackingEfficiencyEtaPythia_"+TrkPtBin_strs[ibin3]+TrkPtBin_strs[ibin3+1]+".pdf"));
 
	phi_canvas[ibin3]->SaveAs((TString)("TrackingEfficiencyPhiPythia_"+TrkPtBin_strs[ibin3]+TrkPtBin_strs[ibin3+1]+".png"));
	phi_canvas[ibin3]->SaveAs((TString)("TrackingEfficiencyPhiPythia_"+TrkPtBin_strs[ibin3]+TrkPtBin_strs[ibin3+1]+".pdf"));

      }else{

	eta_canvas[ibin3]->SaveAs((TString)("TrackingEfficiencyEta_"+TrkPtBin_strs[ibin3]+TrkPtBin_strs[ibin3+1]+".png"));
	eta_canvas[ibin3]->SaveAs((TString)("TrackingEfficiencyEta_"+TrkPtBin_strs[ibin3]+TrkPtBin_strs[ibin3+1]+".pdf"));
 
	phi_canvas[ibin3]->SaveAs((TString)("TrackingEfficiencyPhi_"+TrkPtBin_strs[ibin3]+TrkPtBin_strs[ibin3+1]+".png"));
	phi_canvas[ibin3]->SaveAs((TString)("TrackingEfficiencyPhi_"+TrkPtBin_strs[ibin3]+TrkPtBin_strs[ibin3+1]+".pdf"));



    }
   
  }
  return 0; 
  
}     
