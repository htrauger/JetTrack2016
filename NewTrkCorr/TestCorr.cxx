#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include <vector>
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "trkTree.C"


using namespace std;

TFile *f_cymbal = new TFile("/data/kurtjung/JetTrackCorr_skims/miniSkims/miniTrkTreeCymbal.root");
TTree *t_cymbal = (TTree*)f_cymbal->Get("trkTree");

trkTree *trkTree_cymbal = new trkTree(t_cymbal);

TFile *f_drum = new TFile("/data/kurtjung/JetTrackCorr_skims/miniSkims/miniTrkTreeDrum.root");
TTree *t_drum = (TTree*)f_drum->Get("trkTree");

trkTree *trkTree_drum = new trkTree(t_drum);

TFile *f_data = new TFile("/data/kurtjung/JetTrackCorr_skims/miniSkims/miniTrackTreeData_full.root");
TTree *t_data = (TTree*)f_data->Get("trkTree");

trkTree *trkTree_data = new trkTree(t_data);



int TestCorr(bool do_cuts = 1){


  gStyle->SetOptStat(0);

  TH1D *h_eta = new TH1D("h_eta","h_eta",100,-2.5,2.5);
  TH1D *h_eta_good = new TH1D("h_eta_good","h_eta_good",100,-2.5,2.5);

  TH1D *h_eta_cent = new TH1D("h_eta_cent","h_eta_cent",100,-2.5,2.5);


  TH1D *h_eta_good_dz = new TH1D("h_eta_good_dz","h_eta_good_dz",100,-2.5,2.5);
  TH1D *h_eta_good_dz_cent = new TH1D("h_eta_good_dz_cent","h_eta_good_dz_cent",100,-2.5,2.5);

  TH2D *h_dz_dzerr[5][2];
  TH2D *h_bad_dz = new TH2D("h_bad_dz","",50,-2.5,2.5,50,16.,100.);
  TH2D *h_dz_dxy = new TH2D("h_dz_dxy","",50,0.,0.03,50,0.,0.03);

  TH1D *h_drum_eta = new TH1D("h_drum_eta","h_drum_eta",100,-2.5,2.5);
  TH1D *h_drum_eta_cent = new TH1D("h_drum_eta_cent","h_drum_eta_cent",100,-2.5,2.5);

  TH1D *h_drum_eta_good_dz = new TH1D("h_drum_eta_good_dz","h_drum_eta_good_dz",100,-2.5,2.5);
  TH1D *h_drum_eta_good_dz_cent = new TH1D("h_drum_eta_good_dz_cent","h_drum_eta_good_dz_cent",100,-2.5,2.5);
  TH2D *h_drum_dz_dzerr[5][2];
  TH2D *h_drum_bad_dz = new TH2D("h_drum_bad_dz","",50,-2.5,2.5,50,16.,100.);
  TH2D *h_drum_dz_dxy = new TH2D("h_drum_dz_dxy","",50,0.,0.015,50,0.,0.015);

  TH1D *h_data_eta = new TH1D("h_data_eta","h_data_eta",100,-2.5,2.5);
  TH1D *h_data_eta_cent = new TH1D("h_data_eta_cent","h_data_eta_cent",100,-2.5,2.5);

  TH1D *h_data_eta_good_dz = new TH1D("h_data_eta_good_dz","h_data_eta_good_dz",100,-2.5,2.5);
  TH1D *h_data_eta_good_dz_cent = new TH1D("h_data_eta_good_dz_cent","h_data_eta_good_dz_cent",100,-2.5,2.5);
  TH2D *h_data_dz_dzerr[5][2];
  TH2D *h_data_bad_dz = new TH2D("h_data_bad_dz","",50,-2.5,2.5,50,16.,100.);
  TH2D *h_data_dz_dxy = new TH2D("h_data_dz_dxy","",50,0.,0.03,50,0.,0.03);


  for(int i = 0; i<5; i++){
    h_dz_dzerr[i][0] = new TH2D(Form("h_dz_dzerr_cent%d_midrapidity",i),"",50,0.,.03,50,0.,.01);
    h_dz_dzerr[i][1] = new TH2D(Form("h_dz_dzerr_cent%d_largeta",i),"",50,0.,.03,50,0.,.01);

    h_drum_dz_dzerr[i][0] = new TH2D(Form("h_drum_dz_dzerr_cent%d_midrapidity",i),"",50,0.,.03,50,0.,.01);
    h_drum_dz_dzerr[i][1] = new TH2D(Form("h_drum_dz_dzerr_cent%d_largeta",i),"",50,0.,.03,50,0.,.01);

    h_data_dz_dzerr[i][0] = new TH2D(Form("h_data_dz_dzerr_cent%d_midrapidity",i),"",50,0.,.03,50,0.,.01);
    h_data_dz_dzerr[i][1] = new TH2D(Form("h_data_dz_dzerr_cent%d_largeta",i),"",50,0.,.03,50,0.,.01);
  

    h_dz_dzerr[i][0]->Sumw2();
    h_dz_dzerr[i][1]->Sumw2();
  
    h_drum_dz_dzerr[i][0]->Sumw2();
    h_drum_dz_dzerr[i][1]->Sumw2();

    h_data_dz_dzerr[i][0]->Sumw2();
    h_data_dz_dzerr[i][1]->Sumw2();

  }

  TH1D *h_vz = new TH1D("h_vz","",50,-20.,20.);
  TH1D *h_cent = new TH1D("h_cent","",200,0.,200.);

  TH1D *h_drum_vz = new TH1D("h_drum_vz","",50,-20.,20.);
  TH1D *h_drum_cent = new TH1D("h_drum_cent","",200,0.,200.);

  TH1D *h_data_vz = new TH1D("h_data_vz","",50,-20.,20.);
  TH1D *h_data_cent = new TH1D("h_data_cent","",200,0.,200.);



  h_eta->Sumw2();
  h_eta_good_dz->Sumw2();

  h_eta_cent->Sumw2();
  h_eta_good_dz_cent->Sumw2();


  h_drum_eta->Sumw2();
  h_drum_eta_good_dz->Sumw2();

  h_drum_eta_cent->Sumw2();
  h_drum_eta_good_dz_cent->Sumw2();


  h_data_eta->Sumw2();
  h_data_eta_good_dz->Sumw2();

  h_data_eta_cent->Sumw2();
  h_data_eta_good_dz_cent->Sumw2();


  float vz_w = 1.;

  float cent_count = 0.;
  float per_count = 0.;

  Long64_t n_events = trkTree_cymbal->fChain->GetEntries();

  cout<<n_events<<endl;

  for(int i = 0; i<n_events; i++){

    if(i%100000==0)cout<<"On event: "<<i<<endl;

    trkTree_cymbal->GetEntry(i);
    if((float)abs(trkTree_cymbal->vz) > 15.){
      //  cout<<trkTree_cymbal->vz<<endl;
      continue;
    }

    if(!trkTree_cymbal->evtSel)continue;
    
    vz_w = trkTree_cymbal->weight;

    if((int) trkTree_cymbal->hiBin > 145){
    
      h_cent->Fill(trkTree_cymbal->hiBin,vz_w);
      

      for(int k = 0; k< (int)trkTree_cymbal->trkPt->size(); k++){
  
	if(do_cuts&&(!trkTree_cymbal->passEtCut->at(k)||!trkTree_cymbal->passChi2Cuts->at(k)||!trkTree_cymbal->passNhitCuts->at(k))) continue;
	if(do_cuts&&(trkTree_cymbal->trkDxy->at(k)/trkTree_cymbal->trkDxyErr->at(k) > 3.0)) continue;
	

	if(trkTree_cymbal->trkPt->at(k) < 16.) continue;

	if(trkTree_cymbal->trkDz->at(k)/trkTree_cymbal->trkDzErr->at(k) > 3.0)	h_dz_dxy->Fill(trkTree_cymbal->trkDz->at(k),trkTree_cymbal->trkDxy->at(k),vz_w);
      
	//   std::cout<<i<<" "<<k<<" "<<trkTree_cymbal->trkPt->at(k)<<" "<<endl;

	h_eta->Fill(trkTree_cymbal->trkEta->at(k),vz_w);

	if(trkTree_cymbal->trkDz->at(k)/trkTree_cymbal->trkDzErr->at(k) < 3.0) h_eta_good_dz->Fill(trkTree_cymbal->trkEta->at(k),vz_w);
	else h_bad_dz->Fill(trkTree_cymbal->trkEta->at(k),trkTree_cymbal->trkPt->at(k),vz_w);
     

	//	if(trkTree_cymbal->goodTrack->at(k)==1&&trkTree_cymbal->passEtCut->at(k)==1) h_eta_good->Fill(trkTree_cymbal->trkEta->at(k));

	if(abs(trkTree_cymbal->trkEta->at(k))<1.6)	h_dz_dzerr[4][0]->Fill(trkTree_cymbal->trkDz->at(k),trkTree_cymbal->trkDzErr->at(k),vz_w);
	else	h_dz_dzerr[4][1]->Fill(trkTree_cymbal->trkDz->at(k),trkTree_cymbal->trkDzErr->at(k),vz_w);
      }

    }else if((int) trkTree_cymbal->hiBin < 20){

      h_cent->Fill(trkTree_cymbal->hiBin,vz_w);


      for(int k = 0; k< (int)trkTree_cymbal->trkPt->size(); k++){
	
	  if(do_cuts&&(!trkTree_cymbal->passEtCut->at(k)||!trkTree_cymbal->passChi2Cuts->at(k)||!trkTree_cymbal->passNhitCuts->at(k))) continue;
	
	  if(do_cuts&&(trkTree_cymbal->trkDxy->at(k)/trkTree_cymbal->trkDxyErr->at(k) > 3.0)) continue;


	if(trkTree_cymbal->trkPt->at(k) < 16.) continue;
      
	//   std::cout<<i<<" "<<k<<" "<<trkTree_cymbal->trkPt->at(k)<<" "<<endl;

	h_eta_cent->Fill(trkTree_cymbal->trkEta->at(k),vz_w);

	if(trkTree_cymbal->trkDz->at(k)/trkTree_cymbal->trkDzErr->at(k) < 3.0) h_eta_good_dz_cent->Fill(trkTree_cymbal->trkEta->at(k),vz_w);
     
	//if(trkTree_cymbal->goodTrack->at(k)==1&&trkTree_cymbal->passEtCut->at(k)==1) h_eta_good_cent->Fill(trkTree_cymbal->trkEta->at(k));

	if(abs(trkTree_cymbal->trkEta->at(k))<1.6)	h_dz_dzerr[0][0]->Fill(trkTree_cymbal->trkDz->at(k),trkTree_cymbal->trkDzErr->at(k),vz_w);
	else	h_dz_dzerr[0][1]->Fill(trkTree_cymbal->trkDz->at(k),trkTree_cymbal->trkDzErr->at(k),vz_w);

      
      }
    }
  }
  int bin_edge = h_cent->FindBin(145);

  per_count = h_cent->Integral(bin_edge,200);
  h_eta->Scale(1./per_count);
  h_eta_good_dz->Scale(1./per_count);

  bin_edge = h_cent->FindBin(20);
  cent_count = h_cent->Integral(1,bin_edge);
  h_eta_cent->Scale(1./cent_count);
  h_eta_good_dz_cent->Scale(1./cent_count);


  ///////DRUM LOOP//////////
  cent_count = 0;
  per_count = 0;


 n_events = trkTree_drum->fChain->GetEntries();

  cout<<n_events<<endl;

  for(int i = 0; i<n_events; i++){

    if(i%100000==0)cout<<"On event: "<<i<<endl;

    trkTree_drum->GetEntry(i);

    if(!trkTree_drum->evtSel)continue;

    if((float) abs(trkTree_drum->vz) > 15.) continue;

    vz_w = trkTree_drum->weight;

    if((int) trkTree_drum->hiBin > 160){

      h_drum_cent->Fill(trkTree_drum->hiBin,vz_w);

      for(int k = 0; k< (int)trkTree_drum->trkPt->size(); k++){
	
	if(do_cuts&&(!trkTree_drum->passEtCut->at(k)||!trkTree_drum->passChi2Cuts->at(k)||!trkTree_drum->passNhitCuts->at(k))) continue;

	if(do_cuts&&(trkTree_drum->trkDxy->at(k)/trkTree_drum->trkDxyErr->at(k) > 3.0)) continue;
	if(trkTree_drum->trkPt->at(k) < 16.) continue;

	if(trkTree_drum->trkDz->at(k)/trkTree_drum->trkDzErr->at(k) > 3.0) 	h_drum_dz_dxy->Fill(trkTree_drum->trkDz->at(k),trkTree_drum->trkDxy->at(k),vz_w);
      
	//   std::cout<<i<<" "<<k<<" "<<trkTree_drum->trkPt->at(k)<<" "<<endl;

	h_drum_eta->Fill(trkTree_drum->trkEta->at(k),vz_w);

	//	if(trkTree_drum->goodTrack->at(k)==1&&trkTree_drum->passEtCut->at(k)==1) h_drum_eta_good->Fill(trkTree_drum->trkEta->at(k));

	if(abs(trkTree_drum->trkEta->at(k))<1.6)	h_drum_dz_dzerr[4][0]->Fill(trkTree_drum->trkDz->at(k),trkTree_drum->trkDzErr->at(k),vz_w);
	else	h_drum_dz_dzerr[4][1]->Fill(trkTree_drum->trkDz->at(k),trkTree_drum->trkDzErr->at(k),vz_w);
	

	if(trkTree_drum->trkDz->at(k)/trkTree_drum->trkDzErr->at(k) < 3.0) h_drum_eta_good_dz->Fill(trkTree_drum->trkEta->at(k),vz_w);
	else h_drum_bad_dz->Fill(trkTree_drum->trkEta->at(k),trkTree_drum->trkPt->at(k),vz_w);
     


      }



    }else if((int) trkTree_drum->hiBin < 20){
      h_drum_cent->Fill(trkTree_drum->hiBin,vz_w);

      for(int k = 0; k< (int)trkTree_drum->trkPt->size(); k++){

	if(do_cuts&&(!trkTree_drum->passEtCut->at(k)||!trkTree_drum->passChi2Cuts->at(k)||!trkTree_drum->passNhitCuts->at(k))) continue;
	
	if(do_cuts&&(trkTree_drum->trkDxy->at(k)/trkTree_drum->trkDxyErr->at(k) > 3.0) )continue;

	if(trkTree_drum->trkPt->at(k) < 16.) continue;
      
	//   std::cout<<i<<" "<<k<<" "<<trkTree_drum->trkPt->at(k)<<" "<<endl;

	h_drum_eta_cent->Fill(trkTree_drum->trkEta->at(k),vz_w);
	if(trkTree_drum->trkDz->at(k)/trkTree_drum->trkDzErr->at(k) < 3.0) h_drum_eta_good_dz_cent->Fill(trkTree_drum->trkEta->at(k),vz_w);
     
	//if(trkTree_drum->goodTrack->at(k)==1&&trkTree_drum->passEtCut->at(k)==1) h_drum_eta_good_cent->Fill(trkTree_drum->trkEta->at(k));

	if(abs(trkTree_drum->trkEta->at(k))<1.6)	h_drum_dz_dzerr[0][0]->Fill(trkTree_drum->trkDz->at(k),trkTree_drum->trkDzErr->at(k),vz_w);
	else	h_drum_dz_dzerr[0][1]->Fill(trkTree_drum->trkDz->at(k),trkTree_drum->trkDzErr->at(k),vz_w);

      
      }
    }
  }
  
  bin_edge = h_drum_cent->FindBin(160);

  per_count = h_drum_cent->Integral(bin_edge,200);
  h_drum_eta->Scale(1./per_count);
  h_drum_eta_good_dz->Scale(1./per_count);

  bin_edge = h_drum_cent->FindBin(20);
  cent_count = h_drum_cent->Integral(1,bin_edge);
  h_drum_eta_cent->Scale(1./cent_count);
  h_drum_eta_good_dz_cent->Scale(1./cent_count);

  ///////DATA LOOP//////////



  n_events = trkTree_data->fChain->GetEntries();

  cout<<n_events<<endl;

  for(int i = 0; i<n_events; i++){

    if(i%100000==0)cout<<"On event: "<<i<<endl;

    trkTree_data->GetEntry(i);

    if(!trkTree_data->evtSel)continue;

    if((float) abs(trkTree_data->vz) > 15.) continue;

    vz_w = 1.;

    if((int) trkTree_data->hiBin > 145){
      h_data_cent->Fill(trkTree_data->hiBin,vz_w);
      
      for(int k = 0; k< (int)trkTree_data->trkPt->size(); k++){
	if(do_cuts&&(!trkTree_data->passEtCut->at(k)||!trkTree_data->passChi2Cuts->at(k)||!trkTree_data->passNhitCuts->at(k))) continue;
	if(do_cuts&&(trkTree_data->trkDxy->at(k)/trkTree_data->trkDxyErr->at(k) > 3.0)) continue;

	if(trkTree_data->trkPt->at(k) < 16.) continue;

	h_data_dz_dxy->Fill(trkTree_data->trkDz->at(k),trkTree_data->trkDxy->at(k),vz_w);
      
	//   std::cout<<i<<" "<<k<<" "<<trkTree_data->trkPt->at(k)<<" "<<endl;

	h_data_eta->Fill(trkTree_data->trkEta->at(k),vz_w);

	//	if(trkTree_data->goodTrack->at(k)==1&&trkTree_data->passEtCut->at(k)==1) h_data_eta_good->Fill(trkTree_data->trkEta->at(k));

	if(abs(trkTree_data->trkEta->at(k))<1.6)	h_data_dz_dzerr[4][0]->Fill(trkTree_data->trkDz->at(k),trkTree_data->trkDzErr->at(k),vz_w);
	else	h_data_dz_dzerr[4][1]->Fill(trkTree_data->trkDz->at(k),trkTree_data->trkDzErr->at(k),vz_w);
	

	if(trkTree_data->trkDz->at(k)/trkTree_data->trkDzErr->at(k) < 3.0) h_data_eta_good_dz->Fill(trkTree_data->trkEta->at(k),vz_w);
	else h_data_bad_dz->Fill(trkTree_data->trkEta->at(k),trkTree_data->trkPt->at(k),vz_w);
     


      }



    }else if((int) trkTree_data->hiBin < 20){
        h_data_cent->Fill(trkTree_data->hiBin);

      for(int k = 0; k< (int)trkTree_data->trkPt->size(); k++){
	if(do_cuts&&(!trkTree_data->passEtCut->at(k)||!trkTree_data->passChi2Cuts->at(k)||!trkTree_data->passNhitCuts->at(k))) continue;

	if(do_cuts&&(trkTree_data->trkDxy->at(k)/trkTree_data->trkDxyErr->at(k) > 3.0)) continue;

	if(trkTree_data->trkPt->at(k) < 16.) continue;
      
	//   std::cout<<i<<" "<<k<<" "<<trkTree_data->trkPt->at(k)<<" "<<endl;

	h_data_eta_cent->Fill(trkTree_data->trkEta->at(k),vz_w);

	if(trkTree_data->trkDz->at(k)/trkTree_data->trkDzErr->at(k) < 3.0) h_data_eta_good_dz_cent->Fill(trkTree_data->trkEta->at(k),vz_w);
     

	//if(trkTree_data->goodTrack->at(k)==1&&trkTree_data->passEtCut->at(k)==1) h_data_eta_good_cent->Fill(trkTree_data->trkEta->at(k));

	if(abs(trkTree_data->trkEta->at(k))<1.6)	h_data_dz_dzerr[0][0]->Fill(trkTree_data->trkDz->at(k),trkTree_data->trkDzErr->at(k),vz_w);
	else	h_data_dz_dzerr[0][1]->Fill(trkTree_data->trkDz->at(k),trkTree_data->trkDzErr->at(k),vz_w);

      }
    }
  }
  
  bin_edge = h_data_cent->FindBin(145);

  per_count = h_data_cent->Integral(bin_edge,200);
  h_data_eta->Scale(1./per_count);
  h_data_eta_good_dz->Scale(1./per_count);

  bin_edge = h_data_cent->FindBin(20);
  cent_count = h_data_cent->Integral(1,bin_edge);
  h_data_eta_cent->Scale(1./cent_count);
  h_data_eta_good_dz_cent->Scale(1./cent_count);


  TCanvas *c_eta = new TCanvas("c_eta");

 
  h_eta->SetLineColor(kRed);
  h_eta->SetMarkerColor(kRed);
  h_eta->SetMarkerSize(1);
  h_eta->SetMarkerStyle(10);

  h_eta->SetMinimum(0.);
  h_eta->SetMaximum( h_eta->GetMaximum()*2.);
  h_eta->Draw();

  h_drum_eta->SetLineColor(kBlue);
  h_drum_eta->SetMarkerColor(kBlue);
  h_drum_eta->SetMarkerSize(1);
  h_drum_eta->SetMarkerStyle(10);

  h_drum_eta->Draw("same");

  h_data_eta->SetLineColor(kGreen-2);
  h_data_eta->SetMarkerColor(kGreen-2);
  h_data_eta->SetMarkerSize(1);
  h_data_eta->SetMarkerStyle(10);

  h_data_eta->Draw("same");


  TLegend *l = new TLegend(0.2,0.8,0.5,0.9); 
 l->AddEntry(h_eta,"Cymbal");
 l->AddEntry(h_drum_eta,"Drum");
 l->AddEntry(h_data_eta,"Data");

 l->Draw("same");




  cout<<h_eta->Integral()<<" "<<h_eta_good->Integral()<<" "<<h_eta_good->Integral()/h_eta->Integral()<<endl;
  if(do_cuts) c_eta->SaveAs("Eta_Cent73_Cent100_Raw.png");
  else c_eta->SaveAs("Eta_Cent73_Cent100_NoCuts_Raw.png");


  h_eta_good_dz->SetLineColor(kRed);
  h_eta_good_dz->SetMarkerColor(kRed);
  h_eta_good_dz->SetMarkerSize(1);
  h_eta_good_dz->SetMarkerStyle(10);


  h_eta_good_dz->Divide( h_eta_good_dz,h_eta,1,1,"b");
  h_eta_good_dz->SetMinimum(0.5);
  h_eta_good_dz->SetMaximum(1.5);
  h_eta_good_dz->Draw();

  h_drum_eta_good_dz->SetLineColor(kBlue);
  h_drum_eta_good_dz->SetMarkerColor(kBlue);
  h_drum_eta_good_dz->SetMarkerSize(1);
  h_drum_eta_good_dz->SetMarkerStyle(10);


  h_drum_eta_good_dz->Divide( h_drum_eta_good_dz,h_drum_eta,1,1,"b");
  h_drum_eta_good_dz->Draw("same");

  h_data_eta_good_dz->SetLineColor(kGreen-2);
  h_data_eta_good_dz->SetMarkerColor(kGreen-2);
  h_data_eta_good_dz->SetMarkerSize(1);
  h_data_eta_good_dz->SetMarkerStyle(10);


  h_data_eta_good_dz->Divide( h_data_eta_good_dz,h_data_eta,1,1,"b");
  h_data_eta_good_dz->Draw("same");

 cout<<h_eta->Integral()<<" "<<h_eta_good->Integral()<<" "<<h_eta_good->Integral()/h_eta->Integral()<<endl;
  if(do_cuts) c_eta->SaveAs("Eta_Cent73_Cent100.png");
  else c_eta->SaveAs("Eta_Cent73_Cent100_NoCuts.png");



 h_eta_cent->SetLineColor(kRed);
  h_eta_cent->SetMarkerColor(kRed);
  h_eta_cent->SetMarkerSize(1);
  h_eta_cent->SetMarkerStyle(10);

  h_eta_cent->SetMinimum(0.);
  h_eta_cent->SetMaximum( h_eta_cent->GetMaximum()*2.);
  h_eta_cent->Draw();

  h_drum_eta_cent->SetLineColor(kBlue);
  h_drum_eta_cent->SetMarkerColor(kBlue);
  h_drum_eta_cent->SetMarkerSize(1);
  h_drum_eta_cent->SetMarkerStyle(10);

  h_drum_eta_cent->Draw("same");

  h_data_eta_cent->SetLineColor(kGreen-2);
  h_data_eta_cent->SetMarkerColor(kGreen-2);
  h_data_eta_cent->SetMarkerSize(1);
  h_data_eta_cent->SetMarkerStyle(10);


  h_data_eta_cent->Draw("same");

  l->Draw("same");

  if(do_cuts)  c_eta->SaveAs("Eta_Cent0_Cent10_Raw.png");  
  else c_eta->SaveAs("Eta_Cent0_Cent10_NoCuts_Raw.png");  




 h_eta_good_dz_cent->Divide( h_eta_good_dz_cent,h_eta_cent,1,1,"b");
 h_eta_good_dz_cent->SetLineColor(kRed);
  h_eta_good_dz_cent->SetMarkerColor(kRed);
  h_eta_good_dz_cent->SetMarkerSize(1);
  h_eta_good_dz_cent->SetMarkerStyle(10);

  h_eta_good_dz_cent->SetMinimum(0.5);
  h_eta_good_dz_cent->SetMaximum(1.5);
  h_eta_good_dz_cent->Draw();

  h_drum_eta_good_dz_cent->SetLineColor(kBlue);
  h_drum_eta_good_dz_cent->SetMarkerColor(kBlue);
  h_drum_eta_good_dz_cent->SetMarkerSize(1);
  h_drum_eta_good_dz_cent->SetMarkerStyle(10);


  h_drum_eta_good_dz_cent->Divide( h_drum_eta_good_dz_cent,h_drum_eta_cent,1,1,"b");
  h_drum_eta_good_dz_cent->Draw("same");

  h_data_eta_good_dz_cent->SetLineColor(kGreen-2);
  h_data_eta_good_dz_cent->SetMarkerColor(kGreen-2);
  h_data_eta_good_dz_cent->SetMarkerSize(1);
  h_data_eta_good_dz_cent->SetMarkerStyle(10);


  h_data_eta_good_dz_cent->Divide( h_data_eta_good_dz_cent,h_data_eta_cent,1,1,"b");
  h_data_eta_good_dz_cent->Draw("same");

  l->Draw("same");

  if(do_cuts)  c_eta->SaveAs("Eta_Cent0_Cent10.png");  
  else c_eta->SaveAs("Eta_Cent0_Cent10_NoCuts.png");  




  cout<<h_eta->Integral()<<" "<<h_eta_good->Integral()<<" "<<h_eta_good->Integral()/h_eta->Integral()<<endl;

  TCanvas *c_dz = new TCanvas("c_dz","",10,10,800,800);
  c_dz->Divide(2,2,0.,0.);

  c_dz->cd(1);
  h_dz_dzerr[4][0]->Draw("colz");

  TLine *three = new TLine(0.,0.,0.3,0.1);
  three->Draw("same");


  TLatex *l_1 = new TLatex(0.3,0.15,"|#eta|<1.6, hiBin>145");
  l_1->SetNDC();
  l_1->Draw("same");

  c_dz->cd(2);
 h_dz_dzerr[0][0]->Draw("colz");
  three->Draw("same");

  TLatex *l_2 = new TLatex(0.3,0.15,"|#eta|<1.6, hiBin<20");
  l_2->SetNDC();
  l_2->Draw("same");

 
  c_dz->cd(3);
 
 h_dz_dzerr[4][1]->Draw("colz");
  three->Draw("same");
  TLatex *l_3 = new TLatex(0.3,0.15,"|#eta|>1.6, hiBin>145");
  l_3->SetNDC();
  l_3->Draw("same");

  c_dz->cd(4);
  h_dz_dzerr[0][1]->Draw("colz");
  three->Draw("same");

  TLatex *l_4 = new TLatex(0.3,0.15, "|#eta|>1.6, hiBin<20");
  l_4->SetNDC();
  l_4->Draw("same");

  c_dz->SaveAs("Dz_Error_Corr.png");



  c_dz->cd(1);
  h_data_dz_dzerr[4][0]->Draw("colz");

  three->Draw("same");


  l_1->Draw("same");

  c_dz->cd(2);
  h_data_dz_dzerr[0][0]->Draw("colz");
  three->Draw("same");
  l_2->Draw("same");

 
  c_dz->cd(3);
 
  h_data_dz_dzerr[4][1]->Draw("colz");
  three->Draw("same");
   l_3->Draw("same");

  c_dz->cd(4);
  h_data_dz_dzerr[0][1]->Draw("colz");
  three->Draw("same");
  l_4->Draw("same");



  c_dz->SaveAs("Dz_Data_Error_Corr.png");

  c_dz->cd(1);
  h_drum_dz_dzerr[4][0]->Draw("colz");

  three->Draw("same");
  l_1 = new TLatex(0.3,0.15,"|#eta|<1.6, hiBin>160");
  l_1->SetNDC();
  l_1->Draw("same");

  c_dz->cd(2);
  h_drum_dz_dzerr[0][0]->Draw("colz");
  three->Draw("same");
  l_2->Draw("same");

 
  c_dz->cd(3);
 
  h_drum_dz_dzerr[4][1]->Draw("colz");
  three->Draw("same");
  l_3 = new TLatex(0.3,0.15,"|#eta|>1.6, hiBin>160");
  l_3->SetNDC();
  l_3->Draw("same");
  
  c_dz->cd(4);
  h_drum_dz_dzerr[0][1]->Draw("colz");
  three->Draw("same");
  l_4->Draw("same");

  c_dz->SaveAs("Dz_Drum_Error_Corr.png");


  TFile *f_out = new TFile("TrkCorrStudies.root","RECREATE");

  h_bad_dz->Write();
  h_drum_bad_dz->Write();
  h_data_bad_dz->Write();
  h_dz_dxy->Write();
  h_drum_dz_dxy->Write();
 
  h_dz_dzerr[0][0]->Write();
  h_dz_dzerr[0][1]->Write();
  h_dz_dzerr[4][0]->Write();
  h_dz_dzerr[4][1]->Write();

  h_drum_dz_dzerr[0][0]->Write();
  h_drum_dz_dzerr[0][1]->Write();
  h_drum_dz_dzerr[4][0]->Write();
  h_drum_dz_dzerr[4][1]->Write();

 h_data_dz_dzerr[0][0]->Write();
  h_data_dz_dzerr[0][1]->Write();
  h_data_dz_dzerr[4][0]->Write();
  h_data_dz_dzerr[4][1]->Write();

  h_cent->Write();
  h_drum_cent->Write();
  h_data_cent->Write();

  return 0;
}

