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

int TestCorr2(){


  gStyle->SetOptStat(0);

  const int n_types = 4;
  const int n_cents = 5;
  const int n_eta = 2;
  const int n_cut_scenarios = 4;

  TString types[n_types] = {"Cymbal","Drum","Data","DataOverDrum"};
  float cents[n_cents+1] = {0,10,30,50,70,100};
  TString cent_labels[n_cents] = {"Cent0_Cent10", "Cent10_Cent30", "Cent30_Cent50", "Cent50_Cent70", "Cent70_Cent100"};
  TString cent_tex[n_cents] = {"0-10%", "10-30%", "30-50%", "50-70%", "70-100%"};
  TString eta_bins[n_eta] = {"MidEta","LargeEta"};
  TString cuts[n_cut_scenarios] = {"NoCuts","AllCuts","AllCutsButDz","NoCutsButDz"};

  TH1D *h_eta[n_types][n_cents][n_eta][n_cut_scenarios];
  TH2D *h_dz_dzerr[n_types][n_cents][n_eta][n_cut_scenarios];
 
  TH1D *h_vz[n_types];
  TH1D *h_cent[n_types];

  for(int i_type = 0; i_type < n_types; i_type++){
    for(int i_cent = 0; i_cent < n_cents; i_cent++){
      for(int i_eta = 0; i_eta < n_eta; i_eta++){
	for(int i_cut = 0; i_cut < n_cut_scenarios; i_cut++){
	  
	  h_eta[i_type][i_cent][i_eta][i_cut] = new TH1D((TString)("h_eta_"+types[i_type]+"_"+cent_labels[i_cent]+"_"+eta_bins[i_eta]+"_"+cuts[i_cut]), "",100, -2.5, 2.5);

	  h_dz_dzerr[i_type][i_cent][i_eta][i_cut] = new TH2D((TString)("h_dz_dzerr_"+types[i_type]+"_"+cent_labels[i_cent]+"_"+eta_bins[i_eta]+"_"+cuts[i_cut]), "",50,0.,.03,50,0.,.01);

	  h_eta[i_type][i_cent][i_eta][i_cut]->Sumw2();
	  h_dz_dzerr[i_type][i_cent][i_eta][i_cut]->Sumw2();

	}	
      }    
    }
    h_vz[i_type] = new TH1D((TString)("h_vz_"+types[i_type]),"",50,-20.,20.);
    h_cent[i_type] = new TH1D((TString)("h_cent_"+types[i_type]),"",200,0.,200.);
   
    h_vz[i_type]->Sumw2();
    h_cent[i_type]->Sumw2();
  }


  float vz_w = 1.;

  int i_cent = 0;
  int i_eta = 0;
  int i_type = 0;
 
  TFile *f_cymbal = new TFile("/data/kurtjung/JetTrackCorr_skims/miniSkims/miniTrkTreeCymbal.root");
  TTree *t_cymbal = (TTree*)f_cymbal->Get("trkTree");

  TFile *f_drum = new TFile("/data/kurtjung/JetTrackCorr_skims/miniSkims/miniTrkTreeDrum.root");
  TTree *t_drum = (TTree*)f_drum->Get("trkTree");

  TFile *f_data = new TFile("/data/kurtjung/JetTrackCorr_skims/miniSkims/miniTrackTreeData_full.root");
  TTree *t_data = (TTree*)f_data->Get("trkTree");

  trkTree *trk_tree;

  for(int i_type = 0; i_type < n_types; i_type++){

    if(i_type==0)  trk_tree = new trkTree(t_cymbal);
    else if(i_type==1) trk_tree = new trkTree(t_drum);
    else  trk_tree = new trkTree(t_data);
    
    Long64_t n_events = trk_tree->fChain->GetEntries();


    n_events = 50000;
    cout<<n_events<<endl;

    for(int i = 0; i<n_events; i++){

      if(i%100000==0)cout<<"On file: "+types[i_type]<<"  event: "<<i<<endl;

      trk_tree->GetEntry(i);
    
      if(TMath::Abs(trk_tree->vz) > 15.){
	//  cout<<trk_tree->vz<<endl;
	continue;
      }

      if(!trk_tree->evtSel)continue;
    
      if(i_type <2)  vz_w = trk_tree->weight;
      else vz_w = 1.;

      //  cout<<vz_w<<endl;

      h_cent[i_type]->Fill(trk_tree->hiBin,vz_w);
      h_vz[i_type]->Fill(trk_tree->vz,vz_w);

      for(int j = 0; j< 5; j++){
	if(trk_tree->hiBin/2 >= cents[j] && trk_tree->hiBin/2 < cents[j+1]){
	  i_cent = j;
	}
      }
 
      //   cout<<trk_tree->hiBin/2<<" "<<cents[i_cent]<<" "<<cents[i_cent+1]<<endl;
     
      for(int k = 0; k< (int)trk_tree->trkPt->size(); k++){

	if(trk_tree->trkPt->at(k) < 16.) continue;
  
	
	if(abs(trk_tree->trkEta->at(k)) < 1.6) i_eta = 0;
	else i_eta = 1;

	h_eta[i_type][i_cent][i_eta][0]->Fill(trk_tree->trkEta->at(k),vz_w);

	if(trk_tree->passEtCut->at(k)&&trk_tree->passChi2Cuts->at(k)&&trk_tree->passNhitCuts->at(k)&&trk_tree->trkDxy->at(k)/trk_tree->trkDxyErr->at(k) < 3.0&&trk_tree->trkDz->at(k)/trk_tree->trkDzErr->at(k) < 3.0)   h_eta[i_type][i_cent][i_eta][1]->Fill(trk_tree->trkEta->at(k),vz_w);

	if(trk_tree->passEtCut->at(k)&&trk_tree->passChi2Cuts->at(k)&&trk_tree->passNhitCuts->at(k)&&(trk_tree->trkDxy->at(k)/trk_tree->trkDxyErr->at(k) < 3.0))   h_eta[i_type][i_cent][i_eta][2]->Fill(trk_tree->trkEta->at(k),vz_w);
	
	if(trk_tree->trkDz->at(k)/trk_tree->trkDzErr->at(k) < 3.0)   h_eta[i_type][i_cent][i_eta][3]->Fill(trk_tree->trkEta->at(k),vz_w);

	h_dz_dzerr[i_type][i_cent][i_eta][0]->Fill(trk_tree->trkDz->at(k),trk_tree->trkDzErr->at(k),vz_w);
    
      }
    }
      

    //Normalization
      
    int lbin = h_cent[i_type]->FindBin(cents[i_cent]+.0001);
    int rbin = h_cent[i_type]->FindBin(cents[i_cent+1]-.0001);
    
    double norm = 1./h_cent[i_type]->Integral(lbin,rbin);

    cout<<"norm: "<<norm<<endl;
    
    for(int i_cent = 0; i_cent < n_cents; i_cent++){
      for(int i_eta = 0; i_eta < n_eta; i_eta++){

	h_dz_dzerr[i_type][i_cent][i_eta][0]->Scale(norm);
	
	for(int i_cut = 0; i_cut < n_cut_scenarios; i_cut++){
	
	  //	  cout<<"before: "<<  h_eta[i_type][i_cent][i_eta][i_cut]->Integral()<<endl;
	  h_eta[i_type][i_cent][i_eta][i_cut]->Scale(norm);  

	  //	  cout<<"after: "<<norm<<" "<< h_eta[i_type][i_cent][i_eta][i_cut]->Integral()<<endl;
	}
	
      }
      
    }   
      
  }
  //DRAWING AND OUTPUT
  
  TCanvas *c_eta[n_cut_scenarios];

  TFile *f_out = new TFile("TrkCorrStudies2.root","RECREATE");

    for(int i_cut = 0; i_cut < n_cut_scenarios; i_cut++){
     c_eta[i_cut] = new TCanvas(Form("c_eta%d",i_cut),"",10,10,2500,1000);
      c_eta[i_cut]->Divide(5,2,0.,0.);

      for(int i_eta = 0; i_eta < 1; i_eta++){
	for(int i_cent = 0; i_cent < n_cents; i_cent++){
	  c_eta[i_cut]->cd(5*(i_type+1) - i_cent);

	  h_eta[0][i_cent][i_eta][i_cut]->SetLineColor(kRed);
	  h_eta[0][i_cent][i_eta][i_cut]->SetMarkerColor(kRed);
	  h_eta[0][i_cent][i_eta][i_cut]->SetMarkerSize(1);
	  h_eta[0][i_cent][i_eta][i_cut]->SetMarkerStyle(10);

	  h_eta[0][i_cent][i_eta][i_cut]->SetMinimum(0.);
	  if(i_cent == 0) h_eta[0][i_cent][i_eta][i_cut]->SetMaximum(0.8);
	  if(i_cent == 1) h_eta[0][i_cent][i_eta][i_cut]->SetMaximum(0.2);
	  if(i_cent == 2) h_eta[0][i_cent][i_eta][i_cut]->SetMaximum(0.08);
	  if(i_cent == 3) h_eta[0][i_cent][i_eta][i_cut]->SetMaximum(0.02);
	  if(i_cent == 4) h_eta[0][i_cent][i_eta][i_cut]->SetMaximum(0.007);

	  if(i_eta==0)	  h_eta[0][i_cent][i_eta][i_cut]->Draw();
	  else 	  h_eta[0][i_cent][i_eta][i_cut]->Draw("same");

	  h_eta[1][i_cent][i_eta][i_cut]->SetLineColor(kBlue);
	  h_eta[1][i_cent][i_eta][i_cut]->SetMarkerColor(kBlue);
	  h_eta[1][i_cent][i_eta][i_cut]->SetMarkerSize(1);
	  h_eta[1][i_cent][i_eta][i_cut]->SetMarkerStyle(10);
	  h_eta[1][i_cent][i_eta][i_cut]->Draw("same");

	  h_eta[2][i_cent][i_eta][i_cut]->SetLineColor(kGreen-2);
	  h_eta[2][i_cent][i_eta][i_cut]->SetMarkerColor(kGreen-2);
	  h_eta[2][i_cent][i_eta][i_cut]->SetMarkerSize(1);
	  h_eta[2][i_cent][i_eta][i_cut]->SetMarkerStyle(10);
	  h_eta[2][i_cent][i_eta][i_cut]->Draw("same");

	  h_eta[3][i_cent][i_eta][i_cut] = (TH1D*)h_eta[0][i_cent][i_eta][i_cut]->Clone((TString)("h_eta_"+types[i_type]+"_"+cent_labels[i_cent]+"_"+eta_bins[i_eta]+"_"+cuts[i_cut]));
	  h_eta[3][i_cent][i_eta][i_cut]->Divide(h_eta[2][i_cent][i_eta][i_cut]);

	  h_eta[3][i_cent][i_eta][i_cut]->SetLineColor(kBlack);
	  h_eta[3][i_cent][i_eta][i_cut]->SetMarkerColor(kBlack);


	  if(i_cent==4&&i_eta==0){
	    TLegend *l = new TLegend(0.2,0.6,0.8,0.9); 
	    l->AddEntry(h_eta[0][i_cent][i_eta][i_cut],"Cymbal");
	    l->AddEntry(h_eta[1][i_cent][i_eta][i_cut],"Drum");
	    l->AddEntry(h_eta[2][i_cent][i_eta][i_cut],"Data");
	    l->AddEntry(h_eta[3][i_cent][i_eta][i_cut],"Cymbal/Data");
	    l->SetLineColor(kWhite);
	    l->Draw("same");
	  }
	 
	  TLatex *cent_label = new TLatex(0.6,0.9,cent_tex[i_cent]);
	  cent_label->SetLineColor(kWhite);
	  cent_label->SetNDC();
	  cent_label->Draw();



	  c_eta[i_cut]->cd(10 - i_cent);
	  if(i_eta==0) h_eta[3][i_cent][i_eta][i_cut]->Draw();
	  else  h_eta[3][i_cent][i_eta][i_cut]->Draw("same");
	  
	  h_eta[0][i_cent][i_eta][i_cut]->Write();
	  h_eta[1][i_cent][i_eta][i_cut]->Write();
	  h_eta[2][i_cent][i_eta][i_cut]->Write();
	  h_dz_dzerr[0][i_cent][i_eta][i_cut]->Write();
	  h_dz_dzerr[1][i_cent][i_eta][i_cut]->Write();
	  h_dz_dzerr[2][i_cent][i_eta][i_cut]->Write();


	}
      }
    
     
      c_eta[i_cut]->SaveAs((TString)("Raw_Eta_"+cuts[i_cut]+".png"));
    }
    
    for(int i_type = 0; i_type < n_types; i_type++){
      h_cent[i_type]->Write();
      h_vz[i_type]->Write();
    }

    return 0;
  
}

