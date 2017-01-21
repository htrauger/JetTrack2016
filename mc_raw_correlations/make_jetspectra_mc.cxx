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


Int_t make_jetspectra_mc(float etacut = 1.6, bool is_pp = kFALSE){

  //-------------Input files & set eta-----------------

  TFile *fin, *fout_inc,*fout_lead,*fout_sub, *fspectra, *data_spectra;
  TString jetetacut, etalabel, datalabel;
  float eta_ymax;

  
  
  if( fabs(etacut-1.6)<1e-4 && is_pp == kFALSE){
    fin = new TFile("../mc_raw_correlations/PbPb_5TeV_MC_PythiaHydjet_MixHistos_RecoGenReduced_Merged_refpt_newJetTrackCorrections_fineBin.root","READ");
    fspectra = new TFile("HydJet_JetSpectra.root","RECREATE");
    data_spectra = new TFile("../data_raw_correlations/PbPb_JetSpectra.root");
    jetetacut = "JetEtaCut16";
    etalabel = "|#eta_{jet}|<1.6";
    datalabel = "HydJet_RecoJet_";
    eta_ymax = 0.5;
  
  } else if( fabs(etacut-1.6)<1e-4 && is_pp == kTRUE){
    fin = new TFile("../FROZEN_PREAPPROVAL/mc_raw_correlations/Pythia_RecoJet_GenTrack_Aug23.root","READ");
    fspectra = new TFile("Pythia_Reco_JetSpectra.root","RECREATE");
    data_spectra = new TFile("../data_raw_correlations/pp_JetSpectra.root");
    jetetacut = "JetEtaCut16";
    etalabel = "|#eta_{jet}|<1.6";
    datalabel = "Pythia_RecoJet_";
    eta_ymax = 0.5;
    cout<<"opened files"<<endl;

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

 float TrkPtBins[nTrkPtBins+1] = {0.5, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.5<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};

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
  TH1F* dPhi_hist[nCBins];
   
  TH1F* all_jets_corrpT[nCBins][nPtBins];
  TH1F* all_jets_phi[nCBins][nPtBins];
  TH1F* all_jets_eta[nCBins][nPtBins];
 
  TH1F* data_corrpT[nCBins][nPtBins];
   TH1F* data_eta[nCBins][nPtBins];
   TH1F* data_phi[nCBins][nPtBins];
 
  TH2D* hJetTrackSignalBackground[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundLeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundSubLeading[nCBins][nPtBins][nTrkPtBins];

  TH1D* CorrelationInclusiveEtaProj[nCBins][nPtBins][nTrkPtBins];
  TH1D* CorrelationInclusivePhiProj[nCBins][nPtBins][nTrkPtBins];

  TH1D* CorrelationLeadingMinusSubEtaProj[nCBins][nPtBins][nTrkPtBins];
  TH1D* CorrelationLeadingMinusSubPhiProj[nCBins][nPtBins][nTrkPtBins];

  TH1D *raw_eta_rebin[nCBins][nPtBins];
  TH1D *raw_phi_rebin[nCBins][nPtBins];

  TString desc = "RecoJet_GenTrack";
    
  
  TCanvas *SpectraCanvas_pp = new TCanvas("SpectraCanvas_pp"," ",10,10,1200,400);
  SpectraCanvas_pp->Divide(3,1,0.00000001,0.000000001);
   
  TCanvas *SpectraCanvas = new TCanvas("SpectraCanvas"," ",10,10,1500,400);
  SpectraCanvas->Divide(4,1,0.00000001,0.000000001);

  TCanvas *EtaCanvas = new TCanvas("EtaCanvas"," ",10,10,1500,400);
  EtaCanvas->Divide(4,1,0.00000001,0.000000001);
  
  TCanvas *PhiCanvas = new TCanvas("PhiCanvas"," ",10,10,1500,400);
  PhiCanvas->Divide(4,1,0.00000001,0.000000001);

  TLegend *JetLegend, CorrelationLegend,*lproj;
 
  TString legendentry;
 
  TLatex *jet_tex;
  TLatex *cent_tex;
 

  //-----------------------
  // Start getting histos
  //-----------------------

 
  for (int ibin=0;ibin<nCBins;ibin++){

    if(is_pp&&ibin>0){continue;}

    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 

      cout<<desc + "_all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]<<endl;
      
      all_jets_corrpT[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone((TString) ("all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ));
    
      cout<<"got it"<<endl;

      all_jets_phi[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_all_jets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300" ))->Clone((TString) ("all_jets_phi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ));

      all_jets_eta[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_all_jets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300" ))->Clone((TString) ("all_jets_eta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ));
      cout<<"here"<<endl;

      cout<<"all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] <<endl;

      if(!is_pp){
	data_corrpT[ibin][ibin2] = (TH1F*)data_spectra->Get((TString) ("PbPb_all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"))->Clone((TString)("Data_inclusive_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ));

	data_eta[ibin][ibin2] = (TH1F*)data_spectra->Get((TString) ("PbPb_all_jets_eta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"))->Clone((TString)("Data_inclusive_eta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ));


	data_phi[ibin][ibin2] = (TH1F*)data_spectra->Get((TString) ("PbPb_all_jets_phi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"))->Clone((TString)("Data_inclusive_phi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ));


      }else{

	data_corrpT[ibin][ibin2] = (TH1F*)data_spectra->Get((TString) ("pp_all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"))->Clone((TString)("Data_pp_inclusive_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ));


	data_eta[ibin][ibin2] = (TH1F*)data_spectra->Get((TString) ("pp_all_jets_eta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"))->Clone((TString)("Data_pp_inclusive_eta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ));

	data_phi[ibin][ibin2] = (TH1F*)data_spectra->Get((TString) ("pp_all_jets_phi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"))->Clone((TString)("Data_pp_inclusive_phi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ));


      }
      cout<<"and here"<<endl;
    } // ibin2


   
    //-------------------------
    // Normalize all!
    //------------------------

    float norm_temp;
    float width_temp;
    float width_temp_x;
    float width_temp_y;
    float width_temp_ref;
    
    
    for (int ibin2=0;ibin2<nPtBins;ibin2++){

      all_jets_eta[ibin][ibin2]->Rebin(2);
      data_eta[ibin][ibin2]->Rebin(2);
      data_eta[ibin][ibin2]->Scale(1./2);

      all_jets_phi[ibin][ibin2]->Rebin(4);

      data_phi[ibin][ibin2]->Rebin(4);

      data_phi[ibin][ibin2]->Scale(1./4);

      
      norm_temp =  all_jets_corrpT[ibin][ibin2]->Integral();
      width_temp =  all_jets_corrpT[ibin][ibin2]->GetBinWidth(1);
      all_jets_corrpT[ibin][ibin2]->Scale(1/norm_temp/width_temp);
     
      norm_temp =  all_jets_phi[ibin][ibin2]->Integral();
      width_temp =  all_jets_phi[ibin][ibin2]->GetBinWidth(1);
      all_jets_phi[ibin][ibin2]->Scale(1/norm_temp/width_temp);
      
      norm_temp =  all_jets_eta[ibin][ibin2]->Integral();
      width_temp =  all_jets_eta[ibin][ibin2]->GetBinWidth(1);
      all_jets_eta[ibin][ibin2]->Scale(1/norm_temp/width_temp);

    } // ibin2

    cout<<"ready to plot"<<endl;
    //****************************
    // Plotting starts here!
    //***************************


    int inc_entries =   data_corrpT[ibin][0]->GetEntries();
      
    if(is_pp)SpectraCanvas_pp->cd(1);
    else SpectraCanvas->cd(4-ibin);
    gPad->SetLogy();



    all_jets_corrpT[ibin][0]->GetXaxis()->SetTitle("p_{T}");
    all_jets_corrpT[ibin][0]->GetXaxis()->SetLabelSize(0.06);
    all_jets_corrpT[ibin][0]->GetYaxis()->SetLabelSize(0.06);
    all_jets_corrpT[ibin][0]->GetXaxis()->SetTitleSize(0.08);
    all_jets_corrpT[ibin][0]->GetXaxis()->SetTitleOffset(0.8);

 
    cout<<"now here"<<endl;  


  
    data_corrpT[ibin][0]->SetMarkerStyle(20);
    data_corrpT[ibin][0]->SetMarkerSize(1);
    data_corrpT[ibin][0]->SetMarkerColor(kBlue+2);
    data_corrpT[ibin][0]->SetLineColor(kBlue+2);


    cout<<"and here"<<endl;
    all_jets_corrpT[ibin][0]->SetLineColor(kCyan);
    all_jets_corrpT[ibin][0]->SetMarkerColor(kCyan);
    all_jets_corrpT[ibin][0]->SetMarkerStyle(34);
    all_jets_corrpT[ibin][0]->SetMarkerSize(2);
    

    all_jets_corrpT[ibin][0]->SetAxisRange(50.,499.);
    all_jets_corrpT[ibin][0]->SetMaximum(0.9);
    all_jets_corrpT[ibin][0]->SetMinimum(0.000001);
    all_jets_corrpT[ibin][0]->Draw();
   
    data_corrpT[ibin][0]->Draw("same");
  

    cent_tex = new TLatex(0.2,0.9,CBin_labels[ibin]);
    cent_tex->SetNDC();
    cent_tex->SetTextSize(0.06);
    if(!is_pp) cent_tex->Draw();

    jet_tex = new TLatex(0.2,0.85,"Inclusive");
    jet_tex->SetNDC();
    jet_tex->SetTextSize(0.06);
    jet_tex->Draw();
 
   

    TLegend *l = new TLegend(0.5,0.65,0.9,0.85);
    if(is_pp){
      l->AddEntry( all_jets_corrpT[ibin][0],"Pythia");
      l->AddEntry(data_corrpT[ibin][0],"pp Data"); 
    }else{
      l->AddEntry( all_jets_corrpT[ibin][0],"Pythia+Hydjet");
      l->AddEntry(data_corrpT[ibin][0],"PbPb Data"); 
    }
    TString jetstring = "(";
    jetstring+=inc_entries;
    jetstring+=" Jets)";
    l->AddEntry((TObject*)0,jetstring,"");
    l->SetTextSize(0.06);
    l->SetLineColor(kWhite);
    l->Draw();

    if(is_pp)SpectraCanvas_pp->cd(4);
    else SpectraCanvas->cd(8-ibin);
    gPad->SetLogy();
   
    if(is_pp)SpectraCanvas_pp->cd(7);
    else SpectraCanvas->cd(12-ibin);
    gPad->SetLogy();



    ///ETA///
   
    if(is_pp)SpectraCanvas_pp->cd(2);
    else EtaCanvas->cd(4-ibin);
     
    data_eta[ibin][0]->SetMarkerStyle(20);
    data_eta[ibin][0]->SetMarkerSize(1);
    data_eta[ibin][0]->SetMarkerColor(kBlue+2);
    data_eta[ibin][0]->SetLineColor(kBlue+2);

    all_jets_eta[ibin][0]->SetLineColor(kCyan);
    all_jets_eta[ibin][0]->SetMarkerColor(kCyan);
    all_jets_eta[ibin][0]->SetMarkerStyle(34);
    all_jets_eta[ibin][0]->SetMarkerSize(2);
    

    all_jets_eta[ibin][0]->GetXaxis()->SetTitle("#eta");
    all_jets_eta[ibin][0]->GetXaxis()->SetLabelSize(0.06);
    all_jets_eta[ibin][0]->GetYaxis()->SetLabelSize(0.06);
    all_jets_eta[ibin][0]->GetXaxis()->SetTitleSize(0.08);
    all_jets_eta[ibin][0]->GetXaxis()->SetTitleOffset(0.8);

 

    all_jets_eta[ibin][0]->SetAxisRange(-3.,3.);
    all_jets_eta[ibin][0]->SetMaximum(0.8);
    all_jets_eta[ibin][0]->SetMinimum(0.0);
    all_jets_eta[ibin][0]->Draw();
   
    data_eta[ibin][0]->Draw("same");
  
  
    cent_tex = new TLatex(0.2,0.9,CBin_labels[ibin]);
    cent_tex->SetNDC();
    cent_tex->SetTextSize(0.06);
    if(!is_pp) cent_tex->Draw();

    jet_tex = new TLatex(0.2,0.85,"Inclusive");
    jet_tex->SetNDC();
    jet_tex->SetTextSize(0.06);
    jet_tex->Draw();
 
   

    l = new TLegend(0.5,0.65,0.9,0.85);
    if(is_pp){
      l->AddEntry( all_jets_eta[ibin][0],"Pythia");
      l->AddEntry(data_eta[ibin][0],"pp Data"); 
    }else{
      l->AddEntry( all_jets_eta[ibin][0],"Pythia+Hydjet");
      l->AddEntry(data_eta[ibin][0],"PbPb Data"); 
    }
    jetstring = "(";
    jetstring+=inc_entries;
    jetstring+=" Jets)";
    l->AddEntry((TObject*)0,jetstring,"");
    l->SetTextSize(0.06);
    l->SetLineColor(kWhite);
    l->Draw();

    //PHI//

    if(is_pp)SpectraCanvas_pp->cd(3);
    else PhiCanvas->cd(4-ibin);


    all_jets_phi[ibin][0]->GetXaxis()->SetTitle("#phi");
    all_jets_phi[ibin][0]->GetXaxis()->SetLabelSize(0.06);
    all_jets_phi[ibin][0]->GetYaxis()->SetLabelSize(0.06);
    all_jets_phi[ibin][0]->GetXaxis()->SetTitleSize(0.08);
    all_jets_phi[ibin][0]->GetXaxis()->SetTitleOffset(0.8);

 
   
  
    data_phi[ibin][0]->SetMarkerStyle(20);
    data_phi[ibin][0]->SetMarkerSize(1);
    data_phi[ibin][0]->SetMarkerColor(kBlue+2);
    data_phi[ibin][0]->SetLineColor(kBlue+2);

    all_jets_phi[ibin][0]->SetLineColor(kCyan);
    all_jets_phi[ibin][0]->SetMarkerColor(kCyan);
    all_jets_phi[ibin][0]->SetMarkerStyle(34);
    all_jets_phi[ibin][0]->SetMarkerSize(2);
    

    all_jets_phi[ibin][0]->SetAxisRange(-TMath::Pi(),TMath::Pi()-.0001);
    all_jets_phi[ibin][0]->SetMaximum(0.5);
    all_jets_phi[ibin][0]->SetMinimum(0.0);
    all_jets_phi[ibin][0]->Draw();
   
    data_phi[ibin][0]->Draw("same");
  
    cent_tex = new TLatex(0.2,0.9,CBin_labels[ibin]);
    cent_tex->SetNDC();
    cent_tex->SetTextSize(0.06);
    if(!is_pp)  cent_tex->Draw();

    jet_tex = new TLatex(0.2,0.85,"Inclusive");
    jet_tex->SetNDC();
    jet_tex->SetTextSize(0.06);
    jet_tex->Draw();
 
   

    l = new TLegend(0.5,0.65,0.9,0.85);
    if(is_pp){
      l->AddEntry( all_jets_phi[ibin][0],"Pythia");
      l->AddEntry(data_phi[ibin][0],"pp Data"); 
    }else{
      l->AddEntry( all_jets_phi[ibin][0],"Pythia+Hydjet");
      l->AddEntry(data_phi[ibin][0],"PbPb Data"); 
    }
   jetstring = "(";
    jetstring+=inc_entries;
    jetstring+=" Jets)";
    l->AddEntry((TObject*)0,jetstring,"");
    l->SetTextSize(0.06);
    l->SetLineColor(kWhite);
    l->Draw();


    fspectra->cd();

    all_jets_corrpT[ibin][0]->Write();

  }

  cout<<"ready to save"<<endl;

  if(is_pp){
    cout<<"here"<<endl;
    //   cout<<SpectraCanvas->GetName()<<endl;
    
    SpectraCanvas_pp->SaveAs("JetKinematics_pp_Pythia.png");
    SpectraCanvas_pp->SaveAs("JetKinematics_pp_Pythia.pdf");
 
   cout<<"saved"<<endl;
  }else{
    SpectraCanvas->SaveAs("JetKinematics_PbPb_Hydjet.png");
    SpectraCanvas->SaveAs("JetKinematics_PbPb_Hydjet.pdf");

    EtaCanvas->SaveAs("JetEta_PbPb_Hydjet.png");
    EtaCanvas->SaveAs("JetEta_PbPb_Hydjet.pdf");
    
    PhiCanvas->SaveAs("JetPhi_PbPb_Hydjet.png");
    PhiCanvas->SaveAs("JetPhi_PbPb_Hydjet.pdf");

  }
  return 0;

} // main loop






  
     
