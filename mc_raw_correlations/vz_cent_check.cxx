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


Int_t vz_cent_check(){

  gStyle->SetOptStat(0);

  float norm;
  


  TFile *f_pp = new TFile("../EventAndJet_Output/VertexCentJetInfo_Data_pp_Merged.root");
  TFile *f_pyth = new TFile("Pythia_GenJet_GenTrack_Aug23.root");

  TFile *f_hyd = new TFile("Hydjet_GenJet_GenTrack_Aug23.root");
   TFile *f_PbPb = new TFile("../EventAndJet_Output/VertexCentJetInfo_Data_PbPb_part0.root");


   cout<<"here"<<endl;

  TH1D *hyd_cent_old = (TH1D*)f_hyd->Get("GenJet_GenTrack_Centrality")->Clone("hyd_cent_old");
  norm = hyd_cent_old->Integral("width");
  hyd_cent_old->Scale(1./norm);

  cout<<"got cent"<<endl;

  TH1D *hyd_cent_new = (TH1D*)f_hyd->Get("GenJet_GenTrack_Centrality_Reweighted")->Clone("hyd_cent_new");
  norm = hyd_cent_new->Integral("width");
  hyd_cent_new->Scale(1./norm);

  cout<<"got hyd"<<endl;

  TH1D *PbPb_cent= (TH1D*)f_PbPb->Get("CentDist")->Clone("PbPb_cent");
  PbPb_cent->Rebin(4);
  norm = PbPb_cent->Integral("width");
  PbPb_cent->Scale(1./norm);

  TCanvas *c_hyd_cent = new TCanvas("c_hyd_cent");
  
  hyd_cent_new->SetLineColor(kBlue);
  hyd_cent_new->SetMarkerColor(kBlue);
  hyd_cent_new->SetMarkerSize(1);
  hyd_cent_new->SetMarkerStyle(20);
  hyd_cent_new->GetXaxis()->SetRangeUser(0.,195.);
  hyd_cent_new->GetXaxis()->SetTitle("CMS Centrality Bin (hiBin)");

  hyd_cent_new->Draw();

  hyd_cent_old->SetLineColor(kMagenta+2);
  hyd_cent_old->SetMarkerColor(kMagenta+2);
  hyd_cent_old->SetMarkerSize(1);
  hyd_cent_old->SetMarkerStyle(20);
  hyd_cent_old->Draw("same");
  
  PbPb_cent->SetLineColor(kRed);
  PbPb_cent->SetMarkerColor(kRed);
  PbPb_cent->SetMarkerSize(1);
  PbPb_cent->SetMarkerStyle(20);
  PbPb_cent->Draw("same");

  hyd_cent_new->Draw("same");

  TLegend *legend = new TLegend(0.6,0.8,0.95,0.95);
  legend->AddEntry(PbPb_cent,"Data");
  legend->AddEntry(hyd_cent_old,"MC Before Reweighting");
  legend->AddEntry(hyd_cent_new,"MC After Reweighting");
  legend->Draw();
  
  c_hyd_cent->SaveAs("HydjetCentralityReweighting.png");
  c_hyd_cent->SaveAs("HydjetCentralityReweighting.pdf");


  TH1D *hyd_vz_old = (TH1D*)f_hyd->Get("GenJet_GenTrack_Vz")->Clone("hyd_vz_old");
  norm = hyd_vz_old->Integral("width");
  hyd_vz_old->Scale(1./norm);

  TH1D *hyd_vz_new = (TH1D*)f_hyd->Get("GenJet_GenTrack_Vz_Reweighted")->Clone("hyd_vz_new");
  norm = hyd_vz_new->Integral("width");
  hyd_vz_new->Scale(1./norm);

  TH1D *PbPb_vz= (TH1D*)f_PbPb->Get("VertexDist")->Clone("PbPb_vz");
  PbPb_vz->Rebin(2);
  norm = PbPb_vz->Integral("width");
  PbPb_vz->Scale(1./norm);

  TCanvas *c_hyd_vz = new TCanvas("c_hyd_vz");
  
  hyd_vz_new->SetLineColor(kBlue);
  hyd_vz_new->SetMarkerColor(kBlue);
  hyd_vz_new->SetMarkerSize(1);
  hyd_vz_new->SetMarkerStyle(20);
  hyd_vz_new->GetXaxis()->SetRangeUser(-20.,20.);
  hyd_vz_new->GetXaxis()->SetTitle("Vertex z");

  hyd_vz_new->Draw();

  hyd_vz_old->SetLineColor(kMagenta+2);
  hyd_vz_old->SetMarkerColor(kMagenta+2);
  hyd_vz_old->SetMarkerSize(1);
  hyd_vz_old->SetMarkerStyle(20);
  hyd_vz_old->Draw("same");
  
  PbPb_vz->SetLineColor(kRed);
  PbPb_vz->SetMarkerColor(kRed);
  PbPb_vz->SetMarkerSize(1);
  PbPb_vz->SetMarkerStyle(20);
  PbPb_vz->Draw("same");

  hyd_vz_new->Draw("same");

  legend->Draw();
  
  c_hyd_vz->SaveAs("HydjetVzReweighting.png");
  c_hyd_vz->SaveAs("HydjetVzReweighting.pdf");

 
  cout<<"here"<<endl;

  //  TH1D *pyth_vz_old = (TH1D*)f_pyth->Get("GenJet_GenTrack_Vz_Merged")->Clone("pyth_vz_old");
  TH1D *pyth_vz_old = (TH1D*)f_pyth->Get("GenJet_GenTrack_Vz")->Clone("pyth_vz_old");
  norm = pyth_vz_old->Integral("width");
  pyth_vz_old->Scale(1./norm);


  cout<<"and here"<<endl;
  // TH1D *pyth_vz_new = (TH1D*)f_pyth->Get("GenJet_GenTrack_Vz_Reweighted_Merged")->Clone("pyth_vz_new");
  TH1D *pyth_vz_new = (TH1D*)f_pyth->Get("GenJet_GenTrack_Vz_Reweighted")->Clone("pyth_vz_new");
  norm = pyth_vz_new->Integral("width");
  pyth_vz_new->Scale(1./norm);

  cout<<"got pbpb"<<endl;
  TH1D *pp_vz= (TH1D*)f_pp->Get("VertexDist")->Clone("pp_vz");
  pp_vz->Rebin(2);
  norm = pp_vz->Integral("width");
  pp_vz->Scale(1./norm);

  TCanvas *c_pyth_vz = new TCanvas("c_pyth_vz");
  
  pyth_vz_new->SetLineColor(kBlue);
  pyth_vz_new->SetMarkerColor(kBlue);
  pyth_vz_new->SetMarkerSize(1);
  pyth_vz_new->SetMarkerStyle(20);
  pyth_vz_new->GetXaxis()->SetRangeUser(-20.,20.);
  pyth_vz_new->GetXaxis()->SetTitle("Vertex z");

 

  pyth_vz_old->SetLineColor(kMagenta+2);
  pyth_vz_old->SetMarkerColor(kMagenta+2);
  pyth_vz_old->SetMarkerSize(1);
  pyth_vz_old->SetMarkerStyle(20);
  pyth_vz_old->Draw();

  pyth_vz_new->Draw("same");
  
  pp_vz->SetLineColor(kRed);
  pp_vz->SetMarkerColor(kRed);
  pp_vz->SetMarkerSize(1);
  pp_vz->SetMarkerStyle(20);
  pp_vz->Draw("same");

  pyth_vz_new->Draw("same");
 

legend = new TLegend(0.6,0.8,0.95,0.95);
  legend->AddEntry(pp_vz,"Data");
  legend->AddEntry(pyth_vz_old,"MC Before Reweighting");
  legend->AddEntry(pyth_vz_new,"MC After Reweighting");
  
  
  legend->Draw();
  
  c_pyth_vz->SaveAs("PythiaVzReweighting.png");
  c_pyth_vz->SaveAs("PythiaVzReweighting.pdf");

  



  return 0; 

}     
