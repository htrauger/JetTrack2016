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


Int_t make_jetspectra(float etacut = 1.6, bool is_pp = kFALSE, bool jff_jec = kFALSE){

  //-------------Input files & set eta-----------------

  TFile *fin, *fout_inc,*fout_lead,*fout_sub, *fspectra, *f_pbpb_spectra;
  TString jetetacut, etalabel, datalabel;
  float eta_ymax;
  
  if( fabs(etacut-1.6)<1e-4 && is_pp == kFALSE && !jff_jec){
    fin = new TFile("Data_PbPb_Aug18.root","READ");
    fspectra = new TFile("PbPb_JetSpectra.root","RECREATE");
    jetetacut = "JetEtaCut1.6";
    etalabel = "|#eta_{jet}|<1.6";
    datalabel = "PbPb_";
    eta_ymax = 0.5;
  
  } else if( fabs(etacut-1.6)<1e-4 && is_pp == kFALSE && jff_jec){
    fin = new TFile("../data_raw_correlations/PbPb_5TeV_Data_MixHistos_newJetTrackCorrections_FullDataMerged_fineBin.root","READ");
    fspectra = new TFile("PbPb_JetSpectra_JFFCorr.root","RECREATE");
    jetetacut = "JetEtaCut1.6";
    etalabel = "|#eta_{jet}|<1.6";
    datalabel = "PbPb_";
    eta_ymax = 0.5;
  
  } else if( fabs(etacut-1.6)<1e-4 && is_pp == kTRUE && jff_jec){
    fin = new TFile("","READ");
    fspectra = new TFile("pp_JetSpectra_JFFCorr.root","RECREATE");
    fout_inc = new TFile("pp_Inclusive_Correlations.root","RECREATE");
    f_pbpb_spectra = new TFile("PbPb_JetSpectra_JFFCorr.root","READ");
    jetetacut = "JetEtaCut1.6";
    etalabel = "|#eta_{jet}|<1.6";
    datalabel = "pp_";
    eta_ymax = 0.5;
    
  }  else if( fabs(etacut-1.6)<1e-4 && is_pp == kTRUE && !jff_jec){
    fin = new TFile("Data_pp_July10.root","READ");
    fspectra = new TFile("pp_JetSpectra.root","RECREATE");
    fout_inc = new TFile("pp_Inclusive_Correlations.root","RECREATE");
    f_pbpb_spectra = new TFile("PbPb_JetSpectra.root","READ");
    jetetacut = "JetEtaCut1.6";
    etalabel = "|#eta_{jet}|<1.6";
    datalabel = "pp_";
    eta_ymax = 0.5;
  } else{
    cout<<etacut<<endl;
    cerr<<"No data exists for that jet cut range."<<endl;
    return -1;
  }
  


  //----------------------------------------------------

  
  //-----------------------
  //  Set projection limits
  //-----------------------
  
  float etalim = 1.0;
  float philim = 1.0;

  //-----------------------
  
  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins;
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 9;

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};
 
 float TrkPtBins[nTrkPtBins+1] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
   
  TString ref_name;
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
  

  TH1F* NEvents;
  TH1F* NEvents_test;
  TH1F* NEvents_after_noise;
  TH1F* NEvents_dijets;
  TH1F* all_jets_corrpT[nCBins][nPtBins];
  TH1F* all_jets_phi[nCBins][nPtBins];
  TH1F* all_jets_eta[nCBins][nPtBins];

  TH1F* pbpb_inclusive_corrpT[nCBins][nPtBins];

  TH1F* ratio_inclusive_corrpT[nCBins][nPtBins];

  TH2D* hJetTrackSignalBackground[nCBins][nPtBins][nTrkPtBins];
 
  TH1D* CorrelationInclusiveEtaProj[nCBins][nPtBins][nTrkPtBins];
  TH1D* CorrelationInclusivePhiProj[nCBins][nPtBins][nTrkPtBins];

  TH1D *raw_eta_rebin[nCBins][nPtBins];
  TH1D *raw_phi_rebin[nCBins][nPtBins];

  TString desc = "Data";
    
  TCanvas *JetCanvas = new TCanvas("JetCanvas"," ",10,10,1500,800);
  JetCanvas->Divide(4,2,0.000000,0.0000000);

  TLegend *JetLegend, CorrelationLegend,*lproj;
 

  TString legendentry;
 

  int ibin2 =0;
  //-----------------------
  // Start getting histos
  //-----------------------

  cout<<"start cent loop"<<endl;
  for (int ibin=0;ibin<nCBins;ibin++){

    if(!is_pp||ibin==0){
      all_jets_corrpT[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) (datalabel+"all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));
    
      all_jets_phi[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_all_jets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) (datalabel+"all_jets_phi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));

      all_jets_eta[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_all_jets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) (datalabel+"all_jets_eta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));
    }
       
    if(is_pp){

      cout<<ibin<<" "<<"all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]<<endl;
	
      pbpb_inclusive_corrpT[ibin][ibin2] = (TH1F*)f_pbpb_spectra->Get((TString)("PbPb_all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("pbpb_inclusive_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));
	
      if(ibin>0)  continue;

    }      
 
    cout<<"here"<<endl;
   
    //-------------------------
    // Normalize all!
    //------------------------

    float norm_temp;
    float width_temp;
    float width_temp_x;
    float width_temp_y;
    float width_temp_ref;
   

    norm_temp =  all_jets_corrpT[ibin][ibin2]->GetEntries();
    width_temp =  all_jets_corrpT[ibin][ibin2]->GetBinWidth(1);
    all_jets_corrpT[ibin][ibin2]->Scale(1/norm_temp/width_temp);
     
    norm_temp =  all_jets_phi[ibin][ibin2]->GetEntries();
    width_temp =  all_jets_phi[ibin][ibin2]->GetBinWidth(1);
    all_jets_phi[ibin][ibin2]->Scale(1/norm_temp/width_temp);
      
    norm_temp =  all_jets_eta[ibin][ibin2]->GetEntries();
    width_temp =  all_jets_eta[ibin][ibin2]->GetBinWidth(1);
    all_jets_eta[ibin][ibin2]->Scale(1/norm_temp/width_temp);
  


    cout<<"made it here "<<ibin<<endl;

    //****************************
    // Plotting starts here!
    //***************************
 
    /*
      TLatex *centtex = new TLatex(0.18,0.9,CBin_labels[ibin]);
      centtex->SetNDC();
      centtex->Draw();
    
      TLatex *etatex = new TLatex(0.18,0.85,etalabel);
      etatex->SetNDC();
      etatex->Draw();
    */

  
    fspectra->cd();

    all_jets_corrpT[ibin][0]->Write();
    all_jets_eta[ibin][0]->Write();
    all_jets_phi[ibin][0]->Write();
    
  }

  if(is_pp){
    for(int ibin = 0; ibin<4; ibin++){
    

   JetCanvas->cd(4-ibin);

   
      all_jets_corrpT[0][0]->SetAxisRange(120.,300.,"x");   
      all_jets_corrpT[0][0]->GetXaxis()->SetTitle("p_{T}");
      all_jets_corrpT[0][0]->SetMaximum(0.9);
      all_jets_corrpT[0][0]->SetMinimum(0.00001);
      if(ibin==3){ 
	all_jets_corrpT[0][0]->GetYaxis()->SetTitle("Fraction of jets by pT");
	all_jets_corrpT[0][0]->GetYaxis()->SetTitleSize(0.08);
	all_jets_corrpT[0][0]->GetYaxis()->SetLabelSize(0.06);
	all_jets_corrpT[0][0]->GetYaxis()->SetTitleOffset(0.8);
      }else{
	all_jets_corrpT[0][0]->GetYaxis()->SetLabelSize(0.0);
      }

      all_jets_corrpT[0][0]->SetMarkerStyle(20);
      all_jets_corrpT[0][0]->SetMarkerSize(1);

      all_jets_corrpT[0][0]->SetMarkerColor(kBlack);
      all_jets_corrpT[0][0]->SetLineColor(kBlack);


      pbpb_inclusive_corrpT[ibin][0]->SetAxisRange(120.,300.,"x");   
      pbpb_inclusive_corrpT[ibin][0]->GetXaxis()->SetTitle("p_{T}");
      pbpb_inclusive_corrpT[ibin][0]->SetMaximum(0.9);
      pbpb_inclusive_corrpT[ibin][0]->SetMinimum(0.00001);
      if(ibin==3){ 
	pbpb_inclusive_corrpT[ibin][0]->GetYaxis()->SetTitle("Fraction of jets by pT");
	pbpb_inclusive_corrpT[ibin][0]->GetYaxis()->SetTitleSize(0.08);
	pbpb_inclusive_corrpT[ibin][0]->GetYaxis()->SetLabelSize(0.06);
	pbpb_inclusive_corrpT[ibin][0]->GetYaxis()->SetTitleOffset(0.8);
      }else{
	pbpb_inclusive_corrpT[ibin][0]->GetYaxis()->SetLabelSize(0.0);
      }

 
      pbpb_inclusive_corrpT[ibin][0]->SetMarkerStyle(20);
      pbpb_inclusive_corrpT[ibin][0]->SetMarkerSize(1);

      pbpb_inclusive_corrpT[ibin][0]->SetMarkerColor(kRed);
      pbpb_inclusive_corrpT[ibin][0]->SetLineColor(kRed);


      pbpb_inclusive_corrpT[ibin][0]->Draw();


      
      if(ibin==3){
      TLatex *cent = new TLatex(0.2,0.85, CBin_labels[ibin]);
      cent->SetNDC();
      cent->Draw();
      }else{
TLatex *cent = new TLatex(0.05,0.85, CBin_labels[ibin]);
      cent->SetNDC();
      cent->Draw();
  
      }

    all_jets_corrpT[0][0]->Draw("same");
 
  
    gPad->SetLogy();
  
    JetLegend = new TLegend(0.2,0.2,0.6,0.4);
    JetLegend->SetTextSize(0.06);
    JetLegend->SetLineColor(kWhite);
    JetLegend->AddEntry(all_jets_corrpT[ibin][0],"PbPb ak4PuCalo","lpfe");
    JetLegend->AddEntry(pbpb_inclusive_corrpT[ibin][0],"pp ak4Calo","lpfe");
    if(ibin==3) JetLegend->Draw("same");
  

      ratio_inclusive_corrpT[ibin][0] = (TH1F*) all_jets_corrpT[0][0]->Clone((TString) ("ratio_inclusive_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[0] + "_" + PtBin_strs[0+1]));

   
      ratio_inclusive_corrpT[ibin][0]->Divide(pbpb_inclusive_corrpT[ibin][0]);
   

      JetCanvas->cd(8-ibin);

      ratio_inclusive_corrpT[ibin][0]->SetAxisRange(120.,300.);

      ratio_inclusive_corrpT[ibin][0]->GetXaxis()->SetLabelSize(0.06);
      ratio_inclusive_corrpT[ibin][0]->GetXaxis()->SetNdivisions(505);
      ratio_inclusive_corrpT[ibin][0]->GetXaxis()->SetTitleOffset(0.8);
      ratio_inclusive_corrpT[ibin][0]->GetXaxis()->SetTitleSize(0.08);
  

      if(ibin==3){

	ratio_inclusive_corrpT[ibin][0]->GetYaxis()->SetLabelSize(0.06);
	ratio_inclusive_corrpT[ibin][0]->GetYaxis()->CenterTitle();
   
	ratio_inclusive_corrpT[ibin][0]->GetYaxis()->SetTitleSize(0.08);
	ratio_inclusive_corrpT[ibin][0]->GetYaxis()->SetTitle("PbPb/pp");
  
	ratio_inclusive_corrpT[ibin][0]->GetYaxis()->SetTitleOffset(0.8);
      }

      ratio_inclusive_corrpT[ibin][0]->SetMinimum(0.);
      ratio_inclusive_corrpT[ibin][0]->SetMaximum(2.);
    
      ratio_inclusive_corrpT[ibin][0]->Draw();

      TLine *l = new TLine(120.,1.,300.,1.);
      l->SetLineStyle(2);
      l->Draw();
      if(ibin==3){
      TLatex *inc_tex = new TLatex(0.2,0.9, "PbPb / pp");
      inc_tex->SetNDC();
      inc_tex->Draw();
      }     
    } // centrality bin loop
}

  cout<<"made it here, too"<<endl;

  cout<<datalabel<<endl;
  cout<<jetetacut<<endl;
  
  JetCanvas->cd(0);

  if(jff_jec){
    JetCanvas->SaveAs((TString)("JetSummary_"+datalabel+"JFFCorr.pdf"));
    JetCanvas->SaveAs((TString)("JetSummary_"+datalabel+"JFFCorr.png"));
   }else{
    JetCanvas->SaveAs((TString)("JetSummary_"+datalabel+"NoCorr.pdf"));
    JetCanvas->SaveAs((TString)("JetSummary_"+datalabel+"NoCorr.png"));
  }

  return 0;

} // main loop






  
     
