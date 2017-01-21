#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TF1.h"

#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>

#include <fstream>
#include "TMath.h"
#include <vector>

using namespace std;

int trigger_check(){




  Int_t hiBin;
  Int_t HLT_L1MinimumBiasHF2AND_v1_Prescl;
  Int_t HLT_AK4CaloJet80_Eta5p1ForPPRef_v1;

  Float_t jtpt[206];
  

  TBranch *b_hiBin;
  TBranch *b_HLT_L1MinimumBiasHF2AND_v1_Prescl; 
  TBranch *b_HLT_AK4CaloJet80_Eta5p1ForPPRef_v1;
  TBranch *b_jtpt;  


  TH1D *trig_eff_all = new TH1D("trig_eff_all","",30,0,200);
  TH1D *trig_eff_cent0_cent10 = new TH1D("trig_eff_cent0_cent10","",30,0,200);
  TH1D *trig_eff_cent10_cent30 = new TH1D("trig_eff_cent10_cent30","",30,0,200);
  TH1D *trig_eff_cent30_cent50 = new TH1D("trig_eff_cent30_cent50","",30,0,200);
  TH1D *trig_eff_cent50_cent100 = new TH1D("trig_eff_cent50_cent100","",30,0,200);


  TH1D *trig_denom_all = new TH1D("trig_denom_all","",30,0,200);
  TH1D *trig_denom_cent0_cent10 = new TH1D("trig_denom_cent0_cent10","",30,0,200);
  TH1D *trig_denom_cent10_cent30 = new TH1D("trig_denom_cent10_cent30","",30,0,200);
  TH1D *trig_denom_cent30_cent50 = new TH1D("trig_denom_cent30_cent50","",30,0,200);
  TH1D *trig_denom_cent50_cent100 = new TH1D("trig_denom_cent50_cent100","",30,0,200);

  trig_eff_all->Sumw2();
  trig_eff_cent0_cent10->Sumw2();
  trig_eff_cent10_cent30->Sumw2();
  trig_eff_cent30_cent50->Sumw2();
  trig_eff_cent50_cent100->Sumw2();
  
  trig_denom_all->Sumw2();
  trig_denom_cent0_cent10->Sumw2();
  trig_denom_cent10_cent30->Sumw2();
  trig_denom_cent30_cent50->Sumw2();
  trig_denom_cent50_cent100->Sumw2();
  

  for(int file_n = 0; file_n < 1; file_n++){

    TFile *f_test = TFile::Open(Form( "root://cms-xrd-global.cern.ch///store/user/rbi/merged/MinBias_TuneCUETP8M1_5p02TeV-pythia8-HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2/%d.root",file_n));

    cout<<"0"<<endl;
    TTree *hlt_tree2 = (TTree*)f_test->Get("hltanalysis/HltTree");
     
    hlt_tree2->SetBranchAddress("HLT_L1MinimumBiasHF2AND_v1_Prescl", &HLT_L1MinimumBiasHF2AND_v1_Prescl, &b_HLT_L1MinimumBiasHF2AND_v1_Prescl);
    hlt_tree2->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1ForPPRef_v1", &HLT_AK4CaloJet80_Eta5p1ForPPRef_v1, &b_HLT_AK4CaloJet80_Eta5p1ForPPRef_v1);
    cout<<"1"<<endl;
    TTree *evt_tree = (TTree*) f_test->Get("hiEvtAnalyzer/HiTree");
     int hiBin = 0;
    //evt_tree->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
    //hiBin = evt_tree->hiBin;
     cout<<"2"<<endl;

    TTree *jet_tree = (TTree*)f_test->Get("ak4CaloJetAnalyzer/t");
   
    jet_tree->SetBranchAddress("jtpt",&jtpt, &b_jtpt);

    
    Long64_t n_evt = hlt_tree2->GetEntriesFast();

    cout<<"File #: "<<file_n<<"  Total events: "<<n_evt<<endl;

    n_evt = 100000;
 
    for(int evi = 0; evi < n_evt; evi++){
   
      for(int j4i = 0; j4i< 50 ; j4i++){

	jet_tree->GetEvent(evi);
	if(jtpt[j4i]>500||jtpt[j4i]<25)continue;
      
	hlt_tree2->GetEvent(evi);


	if(HLT_L1MinimumBiasHF2AND_v1_Prescl==0) continue;
     
	trig_denom_all->Fill(jtpt[j4i]);

	evt_tree->GetEvent(evi);
	if(hiBin < 20)   trig_denom_cent0_cent10->Fill(jtpt[j4i]);
	else  if(hiBin >= 20 && hiBin < 60)   trig_denom_cent10_cent30->Fill(jtpt[j4i]);
	else if(hiBin >= 60 && hiBin < 100)   trig_denom_cent30_cent50->Fill(jtpt[j4i]);
	else if(hiBin >= 100 && hiBin < 200)   trig_denom_cent50_cent100->Fill(jtpt[j4i]);


	if(HLT_AK4CaloJet80_Eta5p1ForPPRef_v1==0) continue;

    
	trig_eff_all->Fill(jtpt[j4i]);

	if(hiBin < 20)   trig_eff_cent0_cent10->Fill(jtpt[j4i]);
	else  if(hiBin >= 20 && hiBin < 60)   trig_eff_cent10_cent30->Fill(jtpt[j4i]);
	else if(hiBin >= 60 && hiBin < 100)   trig_eff_cent30_cent50->Fill(jtpt[j4i]);
	else if(hiBin >= 100 && hiBin < 200)   trig_eff_cent50_cent100->Fill(jtpt[j4i]);

      }

    }
 
  TFile *f_out = new TFile(Form("TriggerEfficiency_File%d.root",file_n),"RECREATE");
  trig_eff_all->Write();
  trig_eff_cent0_cent10->Write();
  trig_eff_cent10_cent30->Write();
  trig_eff_cent30_cent50->Write();
  trig_eff_cent50_cent100->Write();
  
  trig_denom_all->Write();
  trig_denom_cent0_cent10->Write();
  trig_denom_cent10_cent30->Write();
  trig_denom_cent30_cent50->Write();
  trig_denom_cent50_cent100->Write();
 

  f_out->Close();

  f_test->Close();

  }


  return 1;
}
