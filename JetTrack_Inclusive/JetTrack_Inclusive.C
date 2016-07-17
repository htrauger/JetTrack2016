/// reco PbPb
#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "mixing_tree.h"
#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>
#include <vector>
#include "TCanvas.h"

#include "Jan_24_pp_Iterative/getTrkCorr.h"

using namespace std;

#define nCBins 4
#define nPtBins 1 
#define nTrkPtBins 9

float trkPtCut=0.5;

int parti = -999;
bool is_data = false;

enum enum_dataset_types {e_Data_PbPb,e_Data_pp,e_HydJet15,e_HydJet30,e_HydJet50, e_HydJet80, e_HydJet120,e_HydJet170,e_HydJet220,e_HydJet280, e_HydJet370,e_Pythia15,e_Pythia30,e_Pythia50, e_Pythia80, e_Pythia120,e_Pythia170,e_Pythia220,e_Pythia280, e_Pythia370, e_n_dataset_types};
int dataset_type_code = -999;

TString dataset_type_strs[e_n_dataset_types] = {"Data_PbPb","Data_pp","HydJet15","HydJet30","HydJet50","HydJet80", "HydJet120", "HydJet170","HydJet220","HydJet280","HydJet370","Pythia15","Pythia30","Pythia50","Pythia80", "Pythia120", "Pythia170","Pythia220","Pythia280","Pythia370"};
//Hydjet80 = 5
//Pythia80 = 14

TString dataset_type_file_names[e_n_dataset_types] = {"ClusterData_PbPb.txt","ClusterData_pp.txt","Hydjet15.txt","Hydjet30.txt","Hydjet50.txt","Hydjet80.txt", "Hydjet120.txt", "Hydjet170.txt","Hydjet220.txt","Hydjet280.txt","Hydjet370.txt","Pythia15.txt","Pythia30.txt","Pythia50.txt","Pythia80.txt", "Pythia120.txt", "Pythia170.txt","Pythia220.txt","Pythia280.txt","Pythia370.txt"};

int dataset_pthats[e_n_dataset_types+1] = {0,0,15,30,50,80,120,170,220,280,370,15,30,50,80,120,170,220,280,370,999};


enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RightGen, SpilledUnderGen, UnmatchedGen, RightReco, SpilledReco, UnmatchedReco, RecoGenSube0,RecoGenNoSube0,GenGenSube0,GenGenNoSube0,MatchedRecoGenSube0,MatchedRecoGenNoSube0,SwappedRecoGenSube0,SwappedRecoGenNoSube0, UnMatchedRecoGenSube0,UnMatchedRecoGenNoSube0,n_data_mc_types};


TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RightGenJet_GenTrack","SpilledUnderJet_GenTrack","UnmatchedGenJet_GenTrack","RightRecoJet_GenTrack","SpilledReco_GenTrack","UnmatchedReco_GenTrack","RecoJet_GenTrack_Sube0","RecoJet_GenTrack_NoSube0","GenJet_GenTrack_Sube0","GenJet_GenTrack_NoSube0","MatchedRecoJet_GenTrack_Sube0","MatchedRecoJet_GenTrack_NoSube0","SwappedRecoJet_GenTrack_Sube0","SwappedRecoJet_GenTrack_NoSube0","UnmatchedRecoJet_GenTrack_Sube0","UnmatchedRecoJet_GenTrack_NoSube0",};
int data_mc_type_code = -999;


float PtBins[nPtBins+1] = {100, 300};
TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};

float CBins[nCBins+1] = {0, 20, 60, 100, 200};
TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};

float TrkPtBins[nTrkPtBins+1] = {0.5, 1, 2, 3, 4, 8, 12, 16, 20, 300};
TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8","TrkPt12","TrkPt16","TrkPt20","TrkPt300" };

float track_pt, trk_corr;

#include "hist_class_def_HT.h"

//Auxillary functions defined below
void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug=false);


///***************************************************************
//     MAIN LOOP STARTS HERE!  
//*************************************************************


int main(int argc, char *argv[]){
 
  assert(argc == 3);
  dataset_type_code = atoi(argv[1]);    //// pick datasets you want to run over
  
  parti = atoi(argv[2]);

  if(dataset_type_code == e_Data_PbPb || dataset_type_code == e_Data_pp){
    is_data = true;
    data_mc_type_code = 0;
  } else{
    is_data = false;
    data_mc_type_code = 2;
  }
    
  bool do_mixing = kFALSE;

  std::cout<<"dataset_type_code is " <<dataset_type_code<<" "<<dataset_type_strs[dataset_type_code]<<endl;
  std::cout << "Running with trkPtCut " << trkPtCut << std::endl;
    
  std::vector<TString> file_names;   file_names.clear();

  ReadFileList( file_names, dataset_type_file_names[dataset_type_code], true);
  
  cout<<"got file"<<endl;

  bool is_pp = kFALSE;
    
  if(dataset_type_code == e_Data_pp||dataset_type_code>10){is_pp = kTRUE;}

  int n_data_mc_types_used = 5;

  hist_class *my_hists[n_data_mc_types];
 
  if(is_data){
    my_hists[data_mc_type_code] = new hist_class(data_mc_type_strs[data_mc_type_code], is_data);
  }else{
  
    for(int mc_type_i = 1; mc_type_i < n_data_mc_types_used; mc_type_i++){
      
      cout<<data_mc_type_strs[mc_type_i]<<endl;

      my_hists[mc_type_i] = new hist_class(data_mc_type_strs[mc_type_i], is_data);
    }

    data_mc_type_code = 2;
  }
   
  cout<<"made hist classes"<<endl;
  
  
  //****************************************
  //        ALL CUTS ARE HERE!
  //****************************************
 
  const double etacut = 1.6;
  // const double pTmaxcut = 300.; // not applied
  const double pTmincut = 120.;
  const double trketamaxcut = 2.4;
  const double max_trkPt = 300;
 
  
  //****************************************


  double cent, eta, pt, phi, rmin, r_reco, jeteta, jetphi, vz,deta, dphi, reco_eta, gen_eta, reco_phi, gen_phi, dr, closest_dr, jet_dir_eta, jet_dir_phi;
  bool foundjet,   is_inclusive;
  int closest_j4i;
	

  double wvz = 1.;
  double wcen = 1.;
 
  TF1 *fit_cen, *fit_vz;

  int unmatched_counter = 0;

  //////////###### PTHAT SAMPLES ###########///////////////

  if(!is_data){
    TFile *f_vertex_cent = new TFile("VertexCentReweightingFits.root","READ");

    if(!is_pp)    fit_cen = (TF1*)f_vertex_cent->Get((TString)("Fit_Cent_"+dataset_type_strs[dataset_type_code]))->Clone((TString)("Fit_Cent_"+dataset_type_strs[dataset_type_code]));

    fit_vz = (TF1*)f_vertex_cent->Get((TString)("Fit_Vz_"+dataset_type_strs[dataset_type_code]))->Clone((TString)("Fit_Vz_"+dataset_type_strs[dataset_type_code]));
 
  }

  //----------------------------------------------------------------
  //   Set tracking efficiency calculation
  //-------------------------------------------------------------

  TrkCorr* trkCorr;

  if(is_pp){
    trkCorr = new TrkCorr("Jan_24_pp_Iterative/");
  
  }else{
    trkCorr = new TrkCorr("Jan18_PbPb_Iterative/");
    
  }

  cout<<"TrkCorr has been set for: "<<dataset_type_strs[dataset_type_code]<<endl;

  //-----------------------------------------------------------------------
  //  ** START ** READING ** THE ** PRIMARY ** TREE **
  //-----------------------------------------------------------------------


  cout<<"Am I pp? "<<is_pp<<endl;

  assert(parti <= (int) file_names.size() );

  for(int fi = 0; fi < (int) file_names.size(); fi++) {
    if( parti >= 0 && parti != fi ) continue;
    TFile *my_file = TFile::Open(file_names.at(fi));
    std::cout << "Current file: " << ", file_name: " << file_names.at(fi) << ", number " << fi << " of " << file_names.size() << std::endl;
    if(my_file->IsZombie()) {
      std::cout << "Is zombie" << std::endl;
    }
    
  
    TTree *inp_tree = (TTree*)my_file->Get("mixing_tree");
    mixing_tree *my_primary = new mixing_tree(inp_tree);
    std::cout << "Successfully retrieved tree from input file!" << std::endl;
    Long64_t n_evt = my_primary->fChain->GetEntriesFast();


    TString me_file_name;
    //  if(is_data&&!is_pp){
    //    me_file_name = "/data/htrauger/PbPb_MinimumBias_7_16/Data_PbPb_MinimumBias_Combined.root";
    // }else{
    if( parti< (int) file_names.size()-1){
      me_file_name = file_names.at(parti+1);
    }else{
      me_file_name = file_names.at(0);
    }


    TFile  *me_file= new TFile(me_file_name,"READ");
    TTree *inp_tree2 = (TTree*)me_file->Get("mixing_tree");
    mixing_tree *me_tree = new mixing_tree(inp_tree2);
  
 

    Long64_t nme = me_tree->fChain->GetEntriesFast();

  
      
    int meptrig = 50;  //note:  this is mixing events per event in which at least one jet (reco or gen in MC) passes selection.
    
    gRandom->SetSeed(0);
    Long64_t  me = gRandom->Rndm()*nme;

   
    TH1D * centbins = new TH1D("centbins","centbins. JUST A DUMMY REALLY", 40, 0.0, 200.0);
    TH1D * vzbins = new TH1D("vzbins","vzbins. JUST A DUMMY REALLY", 30, -15., 15.);
    int jet_cent, jet_vzbin, me_passes;
    vector <Int_t> me_cent;
    vector <Int_t> me_vzbin;
    vector <Int_t> me_evt_sel;
  
   
    if(do_mixing){ 

      cout<<me_file_name<<endl;
      cout<<"There are "<<nme<<" events in the MB file. First, we run through them all once."<<endl;
 
      for(int mei = 0; mei <nme; mei++){
	me_tree->fChain->GetEntry(mei);

	me_passes = 0;

	if(!is_pp){  
	  me_cent.push_back((centbins->FindBin(me_tree->hiBin)));
	}
	if((is_pp&&is_data && me_tree->pBeamScrapingFilter==1 && me_tree->pPAprimaryVertexFilter==1&&me_tree->HBHENoiseFilterResult==1) || 
	   (!is_pp&&is_data  && me_tree->pprimaryVertexFilter==1 && me_tree->pcollisionEventSelection==1&&me_tree->HBHENoiseFilterResult==1) ||
	   !is_data){
	  me_passes = 1;
	}

	me_evt_sel.push_back(me_passes);

	me_vzbin.push_back((vzbins->FindBin(me_tree->vz->at(0))));
	if(vzbins->FindBin(me_tree->vz->at(0))>31){cout<<"THIS IS A PROBLEM!!"<<endl;}
      }
   
    }

    // n_evt = 10000;
    
    ///==========================   Event Loop starts ===================================
    ///==========================   Event Loop starts ===================================
 
    for(int evi = 0; evi < n_evt; evi++) {
    
      my_primary->fChain->GetEntry(evi);


      if (evi%1000==0) std::cout << " I am running on file " << fi+1 << " of " << ((int) file_names.size()) << ", evi: " << evi << " of " << n_evt << std::endl;
      
      Int_t hiBin = 0;
      if(is_data) hiBin = my_primary->hiBin;
   
      vz = my_primary->vz->at(0);

      int ibin2 = 0;  int ibin3=0;

      my_hists[data_mc_type_code]->NEvents->Fill(hiBin/2.0);
   
      if(is_data) {

     	int noise_event_selection = my_primary->HBHENoiseFilterResult;
	if(noise_event_selection==0){ 
	  //cout<<"failed noise"<<endl;
	  continue;      
	}

	if(is_pp){
	  int event_selection = my_primary->pPAprimaryVertexFilter; 
	  if(event_selection==0) {
	    //cout<<"failed vertex filter"<<endl;
	    continue; 
	  }

	  int beam_scraping_selection = my_primary->pBeamScrapingFilter;
	  if(beam_scraping_selection==0){
	    //cout<<"failed beam-scraping"<<endl;
	  continue;
	}

	
	}else{
	
	  int pbpb_event_selection = my_primary->pcollisionEventSelection;

	  if(pbpb_event_selection==0){
	    //cout<<"failed PbPb event selection "<<evi<<endl;
	    continue;
	  }

	  int event_selection = my_primary->pprimaryVertexFilter;
	  if(event_selection==0){
	    cout<<"failed vertex filter "<<event_selection<<" "<<evi<<endl;
	    continue;
	  
	  }
	
	}
      }else{
	
	double evt_pthat = my_primary->pthat;
	if(evt_pthat > dataset_pthats[dataset_type_code+1]){continue; }


      }
      
      if(fabs(vz) > 15.) continue;      

    
      if(!is_data){data_mc_type_code = 4;}

      my_hists[data_mc_type_code]->NEvents_after_noise->Fill(hiBin/2.0);
      my_hists[data_mc_type_code]->Centrality->Fill(hiBin);
      my_hists[data_mc_type_code]->Vz->Fill(vz);
 
      wvz=1;
      wcen=1;

     
      if(!is_data){
	
     	wvz = fit_vz->Eval(vz);
	my_hists[data_mc_type_code]->Vz_new->Fill(vz,wvz);
	
	if(!is_pp){
	  wcen = fit_cen->Eval(1.*hiBin);
	  my_hists[data_mc_type_code]->Centrality_new->Fill(hiBin, wcen);
	}		       
      }
          
      if(!is_data){ data_mc_type_code = 2; } //General event info we put in RecoGen, since we use this for nominal...Setting this here is actually redundant.
    
      for (int ibin=0;ibin<nCBins; ibin ++){


	if (!is_pp&&(my_primary->hiBin<CBins[ibin] || my_primary->hiBin >=CBins[ibin+1])){ continue; }
    
	if(is_pp&&ibin > 0)continue; // no pp reweighting

	foundjet = kFALSE;

	
	for(int j4i = 0; j4i < (int) my_primary->calo_jtpt->size(); j4i++) {

	  if(!is_data){ data_mc_type_code = 2; } //General event info we put in RecoGen, since we use this for nominal...Setting this here is actually redundant.


	  is_inclusive = kFALSE;
	  if( my_primary->calo_trackMax->at(j4i)/my_primary->calo_rawpt->at(j4i) > 0.98 ||my_primary->calo_trackMax->at(j4i)/my_primary->calo_rawpt->at(j4i) < 0.01) continue;
	  if( fabs(my_primary->calo_jteta->at(j4i)) > etacut ) continue;
	  //  if( my_primary->calo_jtpt->at(j4i) > pTmaxcut ) continue;
	  if( my_primary->calo_jtpt->at(j4i) > pTmincut ){
	    is_inclusive = kTRUE;  foundjet = kTRUE;
	  } 
	  	  
	  ibin2 = 0;  ibin3=0;
        
	  for(int pti = 0; pti < nPtBins; pti++) {
	    if (my_primary->calo_jtpt->at(j4i) >=PtBins[pti] && my_primary->calo_jtpt->at(j4i) < PtBins[pti+1])  ibin2 = pti ;
	  }

	  //Determine gen-jet direction once, for use later in sube0 only residual jff-jec scans.

	  jet_dir_eta = my_primary->calo_jteta->at(j4i);
	  jet_dir_phi = my_primary->calo_jtphi->at(j4i);

	  if(!is_data){
	    closest_dr = 999.;
	    closest_j4i = -1;

	    for(int j4i_gen = 0; j4i_gen < (int) my_primary->genpt->size(); j4i_gen++) {

	      gen_phi = my_primary->genphi->at(j4i_gen);

	      if(jet_dir_phi - gen_phi > 2.*TMath::Pi()) gen_phi += 2.*TMath::Pi();
	      if(jet_dir_phi - gen_phi < -2.*TMath::Pi()) gen_phi -= 2.*TMath::Pi();

	      gen_eta = my_primary->geneta->at(j4i_gen);
	  
	      dr = TMath::Sqrt((jet_dir_eta-gen_eta)*(jet_dir_eta-gen_eta)+(jet_dir_phi-gen_phi)*(jet_dir_phi-gen_phi));
	      
	      if(dr<closest_dr){
		closest_j4i = j4i_gen;
		closest_dr = dr;
	      }
	    }// j4i_gen;
	
	    if(is_inclusive){
	      if(closest_dr<0.3){
		jet_dir_eta = my_primary->geneta->at(closest_j4i);	
		jet_dir_phi = my_primary->genphi->at(closest_j4i);
	      }else{
		cout<<"No gen jet found: ("<<jet_dir_eta<<","<<jet_dir_phi<<") pT = "<<my_primary->calo_jtpt->at(j4i)<<endl;
		unmatched_counter++;

	      }
	    }

	  }//!is_data



	  if(is_inclusive == kTRUE){
	    my_hists[data_mc_type_code]->all_jets_corrpT[ibin][ibin2]->Fill(my_primary->calo_jtpt->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->all_jets_phi[ibin][ibin2]->Fill(my_primary->calo_jtphi->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->all_jets_eta[ibin][ibin2]->Fill(my_primary->calo_jteta->at(j4i), wvz*wcen); 
	  }
	
	  if(!is_pp) {cent = my_primary->hiBin; }
	  if(!is_data){data_mc_type_code = 1; }


	  for(int tracks =0; tracks < (int) my_primary->trkPt->size(); tracks++){
	    if(fabs(my_primary->trkEta->at(tracks))>=trketamaxcut) continue;
	    if (my_primary->highPurity->at(tracks)!=1) continue;
	    if(my_primary->trkPt->at(tracks)<=trkPtCut) continue;
	    if(my_primary->trkPt->at(tracks) > max_trkPt) continue;


	    for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
	      if (my_primary->trkPt->at(tracks) >=TrkPtBins[trkpti] && my_primary->trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti;
	    } /// trkpti loop
	  

	    //  Prepare for and call efficiency calculation
	 
	    eta= my_primary->trkEta->at(tracks);
	    pt= my_primary->trkPt->at(tracks);
	    phi= my_primary->trkPhi->at(tracks);
	    rmin = 99;
	    /*
	    for(int ijet=0;ijet<(int) my_primary->calo_jtpt->size();ijet++){
	      jeteta = my_primary->calo_jteta->at(ijet);
	      jetphi = my_primary->calo_jtphi->at(ijet);
	    
	      if(fabs(jeteta)>2 || my_primary->calo_jtpt->at(ijet)<50) continue;
	   
	      r_reco=sqrt(pow(jeteta-eta,2)+pow(acos(cos(jetphi-phi)),2));
	      if(r_reco<rmin)rmin=r_reco;
	    }
	    */
	
	    if(is_pp)  trk_corr = trkCorr->getTrkCorr(pt,eta,phi,0,rmin);
	    else   trk_corr = trkCorr->getTrkCorr(pt,eta,phi,hiBin,rmin);

	    track_pt = my_primary->trkPt->at(tracks);

	    //---------------------------
	    // Now we are ready to fill!
	    //---------------------------
	
	    
	    my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->trkPt->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->trkEta->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->trkPhi->at(tracks),wvz*wcen);
	    
	    my_hists[data_mc_type_code]->TrkPt_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkPt->at(tracks),trk_corr*wvz*wcen);
	    my_hists[data_mc_type_code]->TrkEta_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkEta->at(tracks),trk_corr*wvz*wcen);
	    my_hists[data_mc_type_code]->TrkPhi_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkPhi->at(tracks),trk_corr*wvz*wcen);

	    if(is_inclusive == kTRUE){
	   
	    
	      deta = my_primary->calo_jteta->at(j4i) - my_primary->trkEta->at(tracks);
	      dphi = my_primary->calo_jtphi->at(j4i) - my_primary->trkPhi->at(tracks);
	 
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta,dphi, trk_corr*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, track_pt*trk_corr*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);

	    }
	  	    	    
	  } // Track loop
	  	  
	  if(!is_data){

	    data_mc_type_code = 2;
	    //-------------------------------
	    //   These jets, but gen tracks
	    //-------------------------------

	    for(int tracks =0; tracks < (int) my_primary->pt->size(); tracks++){
	      if(fabs(my_primary->eta->at(tracks))>=trketamaxcut) continue;
	      if(my_primary->pt->at(tracks)<=trkPtCut) continue;
	      if(my_primary->pt->at(tracks) > max_trkPt) continue;
	      if(my_primary->chg->at(tracks)==0) continue;

	  
	      for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
		if (my_primary->pt->at(tracks) >=TrkPtBins[trkpti] && my_primary->pt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
	      } /// trkpti loop
	  
	      my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->pt->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->eta->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->phi->at(tracks),wvz*wcen);

	      track_pt = my_primary->pt->at(tracks);
	    
	      if(is_inclusive == kTRUE){
	   
	    
		deta = my_primary->calo_jteta->at(j4i) - my_primary->eta->at(tracks);
		dphi = my_primary->calo_jtphi->at(j4i) - my_primary->phi->at(tracks);
	 
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
		my_hists[data_mc_type_code]->hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, track_pt*wvz*wcen);

	    
	      }
	    } // Gen Particles
	      /*
	      //------------------------------------------------
	      //Repeat, now subdividing into sube==0 and sube>0.   **duplicate both sube0 scan left for back-compatibility and checks**
	      //-----------------------------------------------
	    
	      if(is_pp||my_primary->sube->at(tracks)==0) data_mc_type_code = 11;
	      else data_mc_type_code = 12;

	    
	      my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->pt->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->eta->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->phi->at(tracks),wvz*wcen);
	    
	      if(is_inclusive == kTRUE){
	   
	    
		deta = jet_dir_eta - my_primary->eta->at(tracks);
		dphi = jet_dir_phi - my_primary->phi->at(tracks);
	 
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
		my_hists[data_mc_type_code]->hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, track_pt*wvz*wcen);
	    
	      }

	      data_mc_type_code =2;
	
	    } // Gen particle loop
	      
	      //---------------
	      //Reco Rightjets
	      //---------------

	  
	      //------- Begin RightJets --------
	
	    reco_eta = my_primary->calo_jteta->at(j4i);
	    reco_phi = my_primary->calo_jtphi->at(j4i);

	    closest_dr = 999.;
	    closest_j4i = -1;

	    for(int j4i_gen = 0; j4i_gen < (int) my_primary->genpt->size(); j4i_gen++) {

	      gen_phi = my_primary->genphi->at(j4i_gen);
	      gen_eta = my_primary->geneta->at(j4i_gen);
	  
	      dr = TMath::Sqrt((reco_eta-gen_eta)*(reco_eta-gen_eta)+(reco_phi-gen_phi)*(reco_phi-gen_phi));
	      
	      if(dr<closest_dr){
		closest_j4i = j4i_gen;
		closest_dr = dr;
	      }
	    }// j4i_gen;
       
	
	    //------- End RightJets -------

	
	    if( closest_dr<0.3&&(my_primary->genpt->at(closest_j4i)>120.)){

	      data_mc_type_code = 8;
	      
	    }else if(closest_dr<0.3&&(my_primary->genpt->at(closest_j4i)<=120.)){

	      data_mc_type_code = 9;
	    }else{
	      data_mc_type_code = 10;
	    }

	    if(is_inclusive == kTRUE){
	      my_hists[data_mc_type_code]->all_jets_corrpT[ibin][ibin2]->Fill(my_primary->calo_jtpt->at(j4i), wvz*wcen); 
	      my_hists[data_mc_type_code]->all_jets_phi[ibin][ibin2]->Fill(my_primary->calo_jtphi->at(j4i), wvz*wcen); 
	      my_hists[data_mc_type_code]->all_jets_eta[ibin][ibin2]->Fill(my_primary->calo_jteta->at(j4i), wvz*wcen); 


	      for(int tracks =0; tracks < (int) my_primary->pt->size(); tracks++){
		if(fabs(my_primary->eta->at(tracks))>=trketamaxcut) continue;
		if(my_primary->pt->at(tracks)<=trkPtCut) continue;
		if(my_primary->chg->at(tracks)==0) continue;
		//	if(my_primary->sube->at(tracks)!=0) continue;  //only pythi for these closures

		track_pt = my_primary->pt->at(tracks);
	  
		for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
		  if (my_primary->pt->at(tracks) >=TrkPtBins[trkpti] && my_primary->pt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
		} /// trkpti loop
	  
		my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->pt->at(tracks),wvz*wcen);
		my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->eta->at(tracks),wvz*wcen);
		my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->phi->at(tracks),wvz*wcen);
	    
		if(is_inclusive == kTRUE){
	   
	    
		  deta = my_primary->calo_jteta->at(j4i) - my_primary->eta->at(tracks);
		  dphi = my_primary->calo_jtphi->at(j4i) - my_primary->phi->at(tracks);
	 
		  while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		  while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		  my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
		  my_hists[data_mc_type_code]->hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, track_pt*wvz*wcen);
	    
		}
	      
		   
	      } // Gen particle loop
	    }

	    */
	  	  
	  }//!is_data 
	
	}  /// Closes jpti loop.  THIS MEANS THAT WE TAKE ALL JETS IN AN EVENT >120 GeV, not just the hardest jets.
      
	if(foundjet==kTRUE){my_hists[data_mc_type_code]->NEvents_test->Fill(hiBin/2.);}


	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////////////
  
      	  
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////   We do ALL of this again, but this time with generated jets.     //////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////////////////////


	if(!is_data){
	data_mc_type_code = 4;

    
	//----------------------------------------------------------------------------
	// Have dijet information.  Time to start filling bins.
	//----------------------------------------------------------------------------

	//Loop over cent bins, but we pick only the right one to fill for PbPb.  We fill all cent bins (properly weighted each time) for pp.

     	for(int j4i = 0; j4i < (int) my_primary->genpt->size(); j4i++) {

	  data_mc_type_code = 4;

	  is_inclusive = kFALSE;
	  if( fabs(my_primary->geneta->at(j4i)) > etacut ) continue;
	  //	  if( my_primary->genpt->at(j4i) > pTmaxcut ) continue;
	  if(my_primary->genpt->at(j4i) > pTmincut ){ is_inclusive = kTRUE; 
	    foundjet = kTRUE;  	
	  } 


	  ibin2 = 0;  ibin3=0;
        
	  for(int pti = 0; pti < nPtBins; pti++) {
	    if (my_primary->genpt->at(j4i) >=PtBins[pti] && my_primary->genpt->at(j4i) < PtBins[pti+1])  ibin2 = pti ;
	  }
      
	  if(is_inclusive == kTRUE){
	 
	    my_hists[data_mc_type_code]->all_jets_corrpT[ibin][ibin2]->Fill(my_primary->genpt->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->all_jets_phi[ibin][ibin2]->Fill(my_primary->genphi->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->all_jets_eta[ibin][ibin2]->Fill(my_primary->geneta->at(j4i), wvz*wcen); 
	  }

	  if(!is_data){	  data_mc_type_code = 3; }
       
	  if(!is_pp) {cent = my_primary->hiBin; }

	  for(int tracks =0; tracks < (int) my_primary->trkPt->size(); tracks++){
	    if(fabs(my_primary->trkEta->at(tracks))>=trketamaxcut) continue;
	    if (my_primary->highPurity->at(tracks)!=1) continue;
	    if(my_primary->trkPt->at(tracks)<=trkPtCut) continue;
	    if(my_primary->trkPt->at(tracks)> max_trkPt) continue;

	  
	    for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
	      if (my_primary->trkPt->at(tracks) >=TrkPtBins[trkpti] && my_primary->trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
	    } /// trkpti loop
	  

	      //  Prepare for and call efficiency calculation

	    eta= my_primary->trkEta->at(tracks);
	    pt= my_primary->trkPt->at(tracks);
	    phi= my_primary->trkPhi->at(tracks);
	    rmin = 99;
	    /*
	    for(int ijet=0;ijet<(int) my_primary->genpt->size();ijet++){
	      jeteta = my_primary->geneta->at(ijet);
	      jetphi = my_primary->genphi->at(ijet);
	    
	      if(fabs(jeteta)>2 || my_primary->genpt->at(ijet)<50) continue;
	   
	      r_reco=sqrt(pow(jeteta-eta,2)+pow(acos(cos(jetphi-phi)),2));
	      if(r_reco<rmin)rmin=r_reco;
	    }
	    */
	 
	    if(is_pp)	    trk_corr = trkCorr->getTrkCorr(pt,eta,phi,0,rmin);
	    else 	    trk_corr = trkCorr->getTrkCorr(pt,eta,phi,cent,rmin);

	    track_pt = my_primary->trkPt->at(tracks);
	 	 	
	    //---------------------------
	    // Now we are ready to fill!
	    //---------------------------
	
	    my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->trkPt->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->trkEta->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->trkPhi->at(tracks),wvz*wcen);
	    
	    my_hists[data_mc_type_code]->TrkPt_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkPt->at(tracks),trk_corr*wvz*wcen);
	    my_hists[data_mc_type_code]->TrkEta_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkEta->at(tracks),trk_corr*wvz*wcen);
	    my_hists[data_mc_type_code]->TrkPhi_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkPhi->at(tracks),trk_corr*wvz*wcen);



	    if(is_inclusive == kTRUE){
	   
	    
	      deta = my_primary->geneta->at(j4i) - my_primary->trkEta->at(tracks);
	      dphi = my_primary->genphi->at(j4i) - my_primary->trkPhi->at(tracks);
	 
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta,dphi, trk_corr*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, track_pt*trk_corr*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    
	    }
	  

	  } // Track loop
	
	    //-------------------------------
	    //   These jets, but gen tracks
	    //-------------------------------

	  data_mc_type_code = 4;
	    
	  for(int tracks =0; tracks < (int) my_primary->pt->size(); tracks++){
	    if(fabs(my_primary->eta->at(tracks))>=trketamaxcut) continue;
	    if(my_primary->pt->at(tracks)<=trkPtCut) continue;
	    if(my_primary->pt->at(tracks)> max_trkPt) continue;
	    if(my_primary->chg->at(tracks)==0) continue;
	 
	    //	    if(my_primary->sube->at(tracks)!=0) continue;

	  
	    for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
	      if (my_primary->pt->at(tracks) >=TrkPtBins[trkpti] && my_primary->pt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
	    } /// trkpti loop
		
	    my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->pt->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->eta->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->phi->at(tracks),wvz*wcen);
	    
	
	    track_pt = my_primary->pt->at(tracks);

	    if(is_inclusive == kTRUE){
	   
	    
	      deta = my_primary->geneta->at(j4i) - my_primary->eta->at(tracks);
	      dphi = my_primary->genphi->at(j4i) - my_primary->phi->at(tracks);
	 
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, track_pt*wvz*wcen);
	    
	    }
	    /*
	    //-------------------------------
	    //  Same but split by sube==0
	    //-----------------------------

	    if(!is_pp){

	      if(my_primary->sube->at(tracks)==0) data_mc_type_code = 13;
	      else data_mc_type_code = 14;
	 
	  
	      my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->pt->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->eta->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->phi->at(tracks),wvz*wcen);
	    
	      if(is_inclusive == kTRUE){
	   
	    
		deta = my_primary->geneta->at(j4i) - my_primary->eta->at(tracks);
		dphi = my_primary->genphi->at(j4i) - my_primary->phi->at(tracks);
	 
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
		my_hists[data_mc_type_code]->hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, track_pt*wvz*wcen);
	    
	      }
	   
	    } //!is_pp
	    */

	    data_mc_type_code =4;
	
	  
	  } // Gen particle loop
	} // Gen jet loop
	} //!is_data
    	  //----------------------------------------------------
	  //      EVENT MIXING STARTS HERE!  (ALL JET TYPES MATCHED AT ONCE)
	  //-----------------------------------------------------
       
	if(do_mixing&&foundjet){ // new scheme:  match event once if we have found a jet
	  jet_cent = 0;
	  if(!is_pp){ jet_cent = centbins->FindBin(my_primary->hiBin);}
	  jet_vzbin = vzbins->FindBin(my_primary->vz->at(0));
      
	  //	  cout << "mixing now " << me <<" "<<nme<<endl;

	  int startovercheck = 0;
	  int mevi = 0;
	  while(mevi< meptrig &&startovercheck <2){  //
	    me++;
	    if(me>=nme){
	      me=0;
	      cout<<"starting over, startovercheck = "<<startovercheck<<" evi= "<<evi<<" jet_cent = "<<jet_cent<<" jet_vzbin = "<<jet_vzbin<<" mevi = "<<mevi<<" "<<me<<" "<<nme<<endl;  
	      assert(startovercheck<20);
	      startovercheck++; 
	    }
	  
	    if(me_evt_sel.at(me)==0) { continue;     }
	  
	    //  Centrality matching
	    if (!is_pp&&(me_cent.at(me)!=jet_cent)){ continue; }

	    // Vz matching
	    if(me_vzbin.at(me)==0||me_vzbin.at(me)==31){ continue; }
	   
	    if(jet_vzbin!= me_vzbin.at(me)){ continue; }
	    
	    me_tree->fChain->GetEntry(me);
	    mevi++;

	
	    for(int j4i = 0; j4i < (int) my_primary->calo_jtpt->size(); j4i++) {

	      if(!is_data){ data_mc_type_code = 2; } //General event info we put in RecoGen, since we use this for nominal...Setting this here is actually redundant.
	      if( my_primary->calo_trackMax->at(j4i)/my_primary->calo_rawpt->at(j4i) > 0.98 ||my_primary->calo_trackMax->at(j4i)/my_primary->calo_rawpt->at(j4i) < 0.01) continue;
	      if( fabs(my_primary->calo_jteta->at(j4i)) > etacut ) continue;
	      //  if( my_primary->calo_jtpt->at(j4i) > pTmaxcut ) continue;
	      if( my_primary->calo_jtpt->at(j4i) < pTmincut) continue;
	     
	      if(!is_data) data_mc_type_code = 1;

	      for(int tracks =0; tracks < (int) me_tree->trkPt->size(); tracks++){
		if(fabs(me_tree->trkEta->at(tracks))>=trketamaxcut) continue;
		if (me_tree->highPurity->at(tracks)!=1) continue;
		if(me_tree->trkPt->at(tracks)<=trkPtCut) continue;
		if(me_tree->trkPt->at(tracks)> max_trkPt) continue;
	  
	  
		for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
		  if (me_tree->trkPt->at(tracks) >=TrkPtBins[trkpti] && me_tree->trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
		} /// trkpti loop

	
		eta= me_tree->trkEta->at(tracks);
		pt= me_tree->trkPt->at(tracks);
		phi= me_tree->trkPhi->at(tracks);
		if(!is_pp){cent = me_cent.at(me);}
	  
		rmin = 99;
		/*
		for(int ijet=0;ijet<(int) me_tree->calo_jtpt->size();ijet++){
		  jeteta = me_tree->calo_jteta->at(ijet);
		  jetphi = me_tree->calo_jtphi->at(ijet);
		 	    
		  if(fabs(jeteta)>2 || me_tree->calo_jtpt->at(ijet)<50) continue;
		  r_reco=sqrt(pow(jeteta-eta,2)+pow(acos(cos(jetphi-phi)),2));
		  if(r_reco<rmin)rmin=r_reco;
		}
		*/

		if(is_pp)	    trk_corr = trkCorr->getTrkCorr(pt,eta,phi,0,rmin);
		else 	    trk_corr = trkCorr->getTrkCorr(pt,eta,phi,cent,rmin);


		//---------------------------
		// Now we are ready to fill!
		//---------------------------
		 
		my_hists[data_mc_type_code]->ME_TrkPt[ibin][ibin2][ibin3]->Fill(me_tree->trkPt->at(tracks),wvz*wcen);
		my_hists[data_mc_type_code]->ME_TrkEta[ibin][ibin2][ibin3]->Fill(me_tree->trkEta->at(tracks),wvz*wcen);
		my_hists[data_mc_type_code]->ME_TrkPhi[ibin][ibin2][ibin3]->Fill(me_tree->trkPhi->at(tracks),wvz*wcen);
	    
		my_hists[data_mc_type_code]->ME_TrkPt_weighted[ibin][ibin2][ibin3]->Fill(me_tree->trkPt->at(tracks),trk_corr*wvz*wcen);
		my_hists[data_mc_type_code]->ME_TrkEta_weighted[ibin][ibin2][ibin3]->Fill(me_tree->trkEta->at(tracks),trk_corr*wvz*wcen);
		my_hists[data_mc_type_code]->ME_TrkPhi_weighted[ibin][ibin2][ibin3]->Fill(me_tree->trkPhi->at(tracks),trk_corr*wvz*wcen);

	    
		   	    
		deta = my_primary->calo_jteta->at(j4i) - me_tree->trkEta->at(tracks);
		dphi = my_primary->calo_jtphi->at(j4i) - me_tree->trkPhi->at(tracks);
	 
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackME[ibin][ibin2][ibin3]->Fill(deta,dphi, trk_corr*wvz*wcen);
		my_hists[data_mc_type_code]->hJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	
	  
	      }  //track mixed event for data


	      if(!is_data){

		data_mc_type_code = 2;
		  
		  
		for(int tracks =0; tracks < (int) me_tree->pt->size(); tracks++){
		  if(fabs(me_tree->eta->at(tracks))>=trketamaxcut) continue;
		  if (me_tree->chg->at(tracks)==0) continue;
		  if(me_tree->pt->at(tracks)<=trkPtCut) continue;
		  if(me_tree->pt->at(tracks)> max_trkPt) continue;
	  
	  
		  for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
		    if (me_tree->pt->at(tracks) >=TrkPtBins[trkpti] && me_tree->pt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
		  } /// trkpti loop
	 
		  //---------------------------
		  // Now we are ready to fill!
		  //---------------------------
		 
		  my_hists[data_mc_type_code]->ME_TrkPt[ibin][ibin2][ibin3]->Fill(me_tree->pt->at(tracks),wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkEta[ibin][ibin2][ibin3]->Fill(me_tree->eta->at(tracks),wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkPhi[ibin][ibin2][ibin3]->Fill(me_tree->phi->at(tracks),wvz*wcen);
	    
		  my_hists[data_mc_type_code]->ME_TrkPt_weighted[ibin][ibin2][ibin3]->Fill(me_tree->pt->at(tracks),trk_corr*wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkEta_weighted[ibin][ibin2][ibin3]->Fill(me_tree->eta->at(tracks),trk_corr*wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkPhi_weighted[ibin][ibin2][ibin3]->Fill(me_tree->phi->at(tracks),trk_corr*wvz*wcen);

	    	   	    
		  deta = my_primary->calo_jteta->at(j4i) - me_tree->eta->at(tracks);
		  dphi = my_primary->calo_jtphi->at(j4i) - me_tree->phi->at(tracks);
	 
		  while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		  while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		  my_hists[data_mc_type_code]->hJetTrackME[ibin][ibin2][ibin3]->Fill(deta,dphi, trk_corr*wvz*wcen);
		  my_hists[data_mc_type_code]->hJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	  
	  
		}  //track mixed event for RecoGen MC
		
	      }// !is_data for gen tracks

	    } //mixing recojet loop

	    

	      // NOW MIXING GENJETS
	    if(!is_data){

	      for(int j4i = 0; j4i < (int) my_primary->genpt->size(); j4i++) {

		if( fabs(my_primary->geneta->at(j4i)) > etacut ) continue;
		if(my_primary->genpt->at(j4i) < pTmincut ) continue;

		if(!is_data) data_mc_type_code = 3;

		for(int tracks =0; tracks < (int) me_tree->trkPt->size(); tracks++){
		  if(fabs(me_tree->trkEta->at(tracks))>=trketamaxcut) continue;
		  if (me_tree->highPurity->at(tracks)!=1) continue;
		  if(me_tree->trkPt->at(tracks)<=trkPtCut) continue;
		  if(me_tree->trkPt->at(tracks)> max_trkPt) continue;
	  
	  
		  for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
		    if (me_tree->trkPt->at(tracks) >=TrkPtBins[trkpti] && me_tree->trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
		  } /// trkpti loop

	
		  eta= me_tree->trkEta->at(tracks);
		  pt= me_tree->trkPt->at(tracks);
		  phi= me_tree->trkPhi->at(tracks);
		  if(!is_pp){cent = me_cent.at(me);}
	  
		  rmin = 99;
		  /*
		  for(int ijet=0;ijet<(int) me_tree->genpt->size();ijet++){
		    jeteta = me_tree->geneta->at(ijet);
		    jetphi = me_tree->genphi->at(ijet);
		 	    
		    if(fabs(jeteta)>2 || me_tree->genpt->at(ijet)<50) continue;
		    r_reco=sqrt(pow(jeteta-eta,2)+pow(acos(cos(jetphi-phi)),2));
		    if(r_reco<rmin)rmin=r_reco;
		  }
		  */

		  if(is_pp)	    trk_corr = trkCorr->getTrkCorr(pt,eta,phi,0,rmin);
		  else 	    trk_corr = trkCorr->getTrkCorr(pt,eta,phi,cent,rmin);
			    
		  //---------------------------
		  // Now we are ready to fill!
		  //---------------------------
		 
		  my_hists[data_mc_type_code]->ME_TrkPt[ibin][ibin2][ibin3]->Fill(me_tree->trkPt->at(tracks),wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkEta[ibin][ibin2][ibin3]->Fill(me_tree->trkEta->at(tracks),wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkPhi[ibin][ibin2][ibin3]->Fill(me_tree->trkPhi->at(tracks),wvz*wcen);
	    
		  my_hists[data_mc_type_code]->ME_TrkPt_weighted[ibin][ibin2][ibin3]->Fill(me_tree->trkPt->at(tracks),trk_corr*wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkEta_weighted[ibin][ibin2][ibin3]->Fill(me_tree->trkEta->at(tracks),trk_corr*wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkPhi_weighted[ibin][ibin2][ibin3]->Fill(me_tree->trkPhi->at(tracks),trk_corr*wvz*wcen);


		  
		  deta = my_primary->geneta->at(j4i) - me_tree->trkEta->at(tracks);
		  dphi = my_primary->genphi->at(j4i) - me_tree->trkPhi->at(tracks);
	 
		  while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		  while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		  my_hists[data_mc_type_code]->hJetTrackME[ibin][ibin2][ibin3]->Fill(deta,dphi, trk_corr*wvz*wcen);
		  my_hists[data_mc_type_code]->hJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	  
	  
		}  //track mixed event for data


		data_mc_type_code = 4;
		  
		  
		for(int tracks =0; tracks < (int) me_tree->pt->size(); tracks++){
		  if(fabs(me_tree->eta->at(tracks))>=trketamaxcut) continue;
		  if (me_tree->chg->at(tracks)==0) continue;
		  if(me_tree->pt->at(tracks)<=trkPtCut) continue;
		  if(me_tree->pt->at(tracks)> max_trkPt) continue;
	  
	  
		  for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
		    if (me_tree->pt->at(tracks) >=TrkPtBins[trkpti] && me_tree->pt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
		  } /// trkpti loop

		  //---------------------------
		  // Now we are ready to fill!
		  //---------------------------
		 
		  my_hists[data_mc_type_code]->ME_TrkPt[ibin][ibin2][ibin3]->Fill(me_tree->pt->at(tracks),wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkEta[ibin][ibin2][ibin3]->Fill(me_tree->eta->at(tracks),wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkPhi[ibin][ibin2][ibin3]->Fill(me_tree->phi->at(tracks),wvz*wcen);
	    
		  deta = my_primary->geneta->at(j4i) - me_tree->eta->at(tracks);
		  dphi = my_primary->genphi->at(j4i) - me_tree->phi->at(tracks);
	 
		  while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		  while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		  my_hists[data_mc_type_code]->hJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	  
			  
		}  //track mixed event for GenGen MC
	
	      } // !is data
	  	  
	    } // gen jet mixing loop

	  } //meptrig events per triggered event
	  
	} //only mix if do_mixing and foundjet
	      
      } //cent is a big loop
          
      	//    cout<<"here at end"<<endl;

    } ///we do EVERYTHING one event at a time.
    
  }//FILE LOOP  (sort of a dummy, since we run on one file at a time).  
 
  cout<<"There were a total of "<<unmatched_counter<<" jets we could not match over all selections."<<endl;

  cout<<"Ready to write"<<endl;
 
  if(is_data){
    my_hists[data_mc_type_code]->Write(0);
    
  }else{
  
    for(int mc_type_i =  1; mc_type_i < n_data_mc_types_used; mc_type_i++){
      cout<<data_mc_type_strs[mc_type_i]<<endl;
      cout<<mc_type_i<<endl;
      my_hists[mc_type_i]->Write(mc_type_i);
    }
  }
  std::cout << "I am FINALLY done!!!" << std::endl;
  
} // end main


void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug)
{
  ifstream file_stream(file_of_names);
  std::string line;
  my_file_names.clear();
  if( debug ) std::cout << "Open file " << file_of_names << " to extract files to run over" << std::endl;
  if( file_stream.is_open() ) {
    if( debug ) std::cout << "Opened " << file_of_names << " for reading" << std::endl;
    int line_num = 0;
    while( !file_stream.eof() ) {
      getline(file_stream, line);
      if( debug ) std::cout << line_num << ": " << line << std::endl;
      TString tstring_line(line);
      if( tstring_line.CompareTo("", TString::kExact) != 0 ) my_file_names.push_back(tstring_line);
      line_num++;
    }
  } else {
    std::cout << "Error, could not open " << file_of_names << " for reading" << std::endl;
    assert(0);
  }
}

