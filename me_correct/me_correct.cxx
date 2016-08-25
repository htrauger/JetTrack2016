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


Int_t me_correct(int data_mc_type_code=0, bool is_pp=kFALSE){

  //-------------Input files & set eta-----------------

  TString datalabel;
  float eta_ymax;

  TLatex *centtex, *pttex;

  TF1 *pol0 = new TF1("pol0","[0]+x-x",-.3,.3);


  TFile *fin, *fin_me, *fin_me2, *fin_me3, *fin_jets, *fout_inc, *fout_inc_pTweighted;


  enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RightGen, SpilledUnderGen, UnmatchedGen, RightReco, SpilledReco, UnmatchedReco, RecoGenSube0,RecoGenNoSube0,GenGenSube0,GenGenNoSube0,MatchedRecoGenSube0,MatchedRecoGenNoSube0,SwappedRecoGenSube0,SwappedRecoGenNoSube0, UnMatchedRecoGenSube0,UnMatchedRecoGenNoSube0,n_data_mc_types};


  TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RightGenJet_GenTrack","SpilledUnderJet_GenTrack","UnmatchedGenJet_GenTrack","RightRecoJet_GenTrack","SpilledReco_GenTrack","UnmatchedReco_GenTrack","RecoJet_GenTrack_Sube0","RecoJet_GenTrack_NoSube0","GenJet_GenTrack_Sube0","GenJet_GenTrack_NoSube0","MatchedRecoJet_GenTrack_Sube0","MatchedRecoJet_GenTrack_NoSube0","SwappedRecoJet_GenTrack_Sube0","SwappedRecoJet_GenTrack_NoSube0","UnmatchedRecoJet_GenTrack_Sube0","UnmatchedRecoJet_GenTrack_NoSube0",};


  TString desc = data_mc_type_strs[data_mc_type_code];


  if(data_mc_type_code==0){
    if(is_pp){
      datalabel = "pp";
      fin = new TFile("../data_raw_correlations/Data_pp_Aug18.root","READ");
      fin_me = new TFile("../data_raw_correlations/Data_pp_July10.root","READ");
      fin_me2 = new TFile("../data_raw_correlations/Data_pp_July10.root","READ");
      fout_inc = new TFile("pp_Inclusive_Correlations.root","RECREATE");
    }else{
      datalabel = "PbPb";
      fin = new TFile("../data_raw_correlations/Data_PbPb_Aug18.root","READ");
      fin_me = new TFile("../data_raw_correlations/Data_PbPb_Partial_Aug14.root","READ");
      fin_me3 = new TFile("../data_raw_correlations/Data_PbPb_ME_Aug23.root","READ");
      fin_me2 = new TFile("../data_raw_correlations/Data_pp_July10.root","READ");
      fout_inc = new TFile("PbPb_Inclusive_Correlations.root","RECREATE");
    }
  }else{
    if(is_pp){
      datalabel = "Pythia";
      if(data_mc_type_code==1||data_mc_type_code==2||data_mc_type_code==11||data_mc_type_code==12) fin = new TFile((TString)("../mc_raw_correlations/Pythia_"+data_mc_type_strs[data_mc_type_code]+"_Aug23.root"),"READ");
      else fin = new TFile((TString)("../mc_raw_correlations/Pythia_"+data_mc_type_strs[data_mc_type_code]+"_Aug23.root"),"READ");
    
      fin_me = new TFile((TString)("../mc_raw_correlations/Pythia_"+data_mc_type_strs[data_mc_type_code]+"_Aug10.root"),"READ");
      if(data_mc_type_code==1)    fin_jets = new TFile((TString)("../mc_raw_correlations/Pythia_RecoJet_GenTrack_Aug23.root"),"READ");
      if(data_mc_type_code==3)    fin_jets = new TFile((TString)("../mc_raw_correlations/Pythia_GenJet_GenTrack_Aug23.root"),"READ");
      fin_me2 = new TFile("../data_raw_correlations/Data_pp_July10.root","READ");
      fout_inc = new TFile((TString)("Pythia_"+data_mc_type_strs[data_mc_type_code]+"_Inclusive_Correlations.root"),"RECREATE");

    }else{
      datalabel = "Hydjet";
      if(data_mc_type_code==1||data_mc_type_code==2||data_mc_type_code==11||data_mc_type_code==12)  fin = new TFile((TString)("../mc_raw_correlations/HydJet_"+data_mc_type_strs[data_mc_type_code]+"_Aug23.root"),"READ");
      else fin = new TFile((TString)("../mc_raw_correlations/HydJet_"+data_mc_type_strs[data_mc_type_code]+"_Aug23.root"),"READ"); 
      if(data_mc_type_code==1||data_mc_type_code==11||data_mc_type_code==12)    fin_jets = new TFile((TString)("../mc_raw_correlations/HydJet_RecoJet_GenTrack_Aug23.root"),"READ");
      if(data_mc_type_code==3||data_mc_type_code==13||data_mc_type_code==14)    fin_jets = new TFile((TString)("../mc_raw_correlations/HydJet_GenJet_GenTrack_Aug23.root"),"READ");

      if(data_mc_type_code==2)  fin_me = new TFile((TString)("../mc_raw_correlations/HydJet_RecoJet_GenTrack_NoSube0_Aug23.root"),"READ");
      if(data_mc_type_code==4)  fin_me = new TFile((TString)("../mc_raw_correlations/HydJet_GenJet_GenTrack_NoSube0_Aug23.root"),"READ");
      fin_me2 = new TFile("../data_raw_correlations/Data_pp_July10.root","READ");
      fout_inc = new TFile((TString)("Hydjet_"+data_mc_type_strs[data_mc_type_code]+"_Inclusive_Correlations.root"),"RECREATE");
      }
    }
    //----------------------------------------------------



  TCanvas *me_proj_canvas, *result_proj_canvas, *result_proj_eta_canvas;

  if(is_pp){

    me_proj_canvas = new TCanvas("me_proj_canvas","",0,0,400,2400);
    me_proj_canvas->Divide(1,9,0.0000,0.0000);


    result_proj_canvas = new TCanvas("result_proj_canvas","",0,0,400,2400);
    result_proj_canvas->Divide(1,9,0.0000,0.0000);

    result_proj_eta_canvas = new TCanvas("result_proj_canvas_eta","",0,0,400,2400);
    result_proj_eta_canvas->Divide(1,9,0.0000,0.0000);

  }else{

    me_proj_canvas = new TCanvas("me_proj_canvas","",0,0,1600,2400);
    me_proj_canvas->Divide(4,9,0.0000,0.0000);


    result_proj_canvas = new TCanvas("result_proj_canvas","",0,0,1600,2400);
    result_proj_canvas->Divide(4,9,0.0000,0.0000);

    result_proj_eta_canvas = new TCanvas("result_proj_canvas_eta","",0,0,1600,2400);
    result_proj_eta_canvas->Divide(4,9,0.0000,0.0000);

  }

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins;
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 9;

  float me00,me00_lead,me00_sub;
  int me00binl, me00binr, n_me00bins;

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins+1] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs2[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
   TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
  
  TString ref_name;
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
  
  TH1F* all_jets_corrpT[nCBins][nPtBins];
 
  TH2D* hJetTrackSignalBackground[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackground_pTweighted[nCBins][nPtBins][nTrkPtBins];
  
  TH2D* hJetTrackME[nCBins][nPtBins][nTrkPtBins];
 
  TH2D *yield_inc[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_inc_pTweighted[nCBins][nPtBins][nTrkPtBins];

  TH2D *yield_bgsub[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_bgsub_pTweighted[nCBins][nPtBins][nTrkPtBins];
 
 
  TH1D *yield_inc_proj[nCBins][nPtBins][nTrkPtBins];
  TH1D *yield_inc_proj_eta[nCBins][nPtBins][nTrkPtBins];
  TH1D *test_me_inc_proj[nCBins][nPtBins][nTrkPtBins];
 
  TH1D *me_proj[nCBins][nPtBins][nTrkPtBins];

  TH2D* background[nCBins][nPtBins][nTrkPtBins];
  TH2D* background_pTweighted[nCBins][nPtBins][nTrkPtBins];
  TH1D* background_left[nCBins][nPtBins][nTrkPtBins];
  TH1D* background_right[nCBins][nPtBins][nTrkPtBins];
  TH1D* background_proj[nCBins][nPtBins][nTrkPtBins];

  float norm_temp, width_temp, width_temp_x, width_temp_y, width_temp_ref, max_bin, max_cont,bc,err,err1, dx_phi, dx_eta, temp1;
    
  //-----------------------
  // Start getting histos
  //-----------------------

  
  for (int ibin=0;ibin<nCBins;ibin++){

    if(is_pp && ibin > 0 ) continue;
  
    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 

      cout<<desc + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]<<endl;

         
      if(data_mc_type_code ==1||data_mc_type_code ==11||data_mc_type_code ==12){  

	all_jets_corrpT[ibin][ibin2] = (TH1F*)fin_jets->Get((TString) ("RecoJet_GenTrack_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));

      }else  if(data_mc_type_code ==3||data_mc_type_code ==13||data_mc_type_code ==14){  

	all_jets_corrpT[ibin][ibin2] = (TH1F*)fin_jets->Get((TString) ("GenJet_GenTrack_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));

	}else{
	all_jets_corrpT[ibin][ibin2] = (TH1F*)fin->Get((TString) (data_mc_type_strs[data_mc_type_code] + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));


      }
         
      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){
	if(data_mc_type_code > 0 && (data_mc_type_code%2==0||data_mc_type_code > 10)){ //if gen

	  hJetTrackSignalBackground[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackSignalBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
		  
	  hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackSignalBackground_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	 
	  if(ibin3>5||is_pp){
	    
	    hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin_me2->Get((TString)("Data_hJetTrackME"+ CBin_strs[0] + "_" + CBin_strs[1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	
  	  }else {

	      hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackSignalBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    cout<<"here"<<endl;

	    if(!is_pp&&data_mc_type_code == 2&&!(ibin3==0&&ibin==0))    hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)("RecoJet_GenTrack_NoSube0_hJetTrackSignalBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    //if(!is_pp&&data_mc_type_code == 4&&!(ibin3==0&&ibin==0))	      hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)("GenJet_GenTrack_NoSube0_hJetTrackSignalBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    if(!is_pp&&data_mc_type_code == 4&&!(ibin==0&&ibin3==1))	      hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)("GenJet_GenTrack_NoSube0_hJetTrackSignalBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    if(!is_pp&&data_mc_type_code == 2&&(ibin3<1&&ibin==0))	      hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)("RecoJet_GenTrack_NoSube0_hJetTrackSignalBackground_notrkcorr"+ CBin_strs[1] + "_" + CBin_strs[2] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    //	    if(!is_pp&&data_mc_type_code == 4&&(ibin3==0&&ibin==0))	      hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)("GenJet_GenTrack_NoSube0_hJetTrackSignalBackground_notrkcorr"+ CBin_strs[1] + "_" + CBin_strs[2] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  }

	}else{
	  cout<<desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]<<endl;

	  hJetTrackSignalBackground[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackSignalBackground_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	 
	  if(ibin3>5||is_pp){  //high-pT, just use pp for geometry
 	    hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin_me2->Get((TString)("Data_hJetTrackME"+ CBin_strs[0] + "_" + CBin_strs[1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	 
	    //	  }else if(data_mc_type_code!=0||(ibin3==0&&is_pp)){
	  }else if(data_mc_type_code!=0||ibin3==0){
	    cout<<"using this one"<<endl;
	    hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc+"_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    /* 
	  }else if(data_mc_type_code==0&is_pp==0&&ibin3==0){

	    hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin_me3->Get((TString)(desc+"_hJetTrackME"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    */  
	  }else{
	    hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)(desc+"_hJetTrackME"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  }

	}
      } /// ibin3

    } // ibin2

    //-------------------------
    // Normalize and subtract bkg
    //------------------------


    

    cout<<"got hists for centbin "<<ibin<<endl;
   
    for (int ibin2=0;ibin2<nPtBins;ibin2++){

      for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){

	width_temp_x = 1.;
	width_temp_y = 1.;

	norm_temp = all_jets_corrpT[ibin][ibin2]->Integral();
	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	hJetTrackME[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y/50.);
   
	//     INCLUSIVE

	fout_inc->cd();
	
	if(is_pp)  	me_proj_canvas->cd(ibin3+1);
	else me_proj_canvas->cd(4*(ibin3+1)-ibin);

	yield_inc[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackground[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	yield_inc_pTweighted[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Clone((TString)("Yield_pTweighted_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	int lbin = 1;
	int rbin = 200;

	if((data_mc_type_code!=0 && ibin3<6)||ibin3==0){
	  lbin = 82;
	  rbin = 119;
	}

	if((data_mc_type_code==12||data_mc_type_code==14||data_mc_type_code==2||data_mc_type_code==4) && (ibin3<6&&!is_pp)){
	  lbin = 82;
	  rbin = 200;
	}

	      
	me_proj[ibin][ibin2][ibin3] = hJetTrackME[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),lbin,rbin);
	
	cout<<ibin<<" "<<ibin3<<" "<< hJetTrackME[ibin][ibin2][ibin3]->GetEntries()<<endl;;

	me_proj[ibin][ibin2][ibin3]->Scale (1./(rbin-lbin+1));
	me_proj[ibin][ibin2][ibin3] -> Fit("pol0","","",-.4,.4);

	me00 = 	me_proj[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);



	for(int k = 1;  k<yield_inc[ibin][ibin2][ibin3]->GetNbinsX(); k++){
	  temp1 = me_proj[ibin][ibin2][ibin3]->GetBinContent(k);
	  err1 = me_proj[ibin][ibin2][ibin3]->GetBinError(k);
	    
	  for(int m = 1;  m<yield_inc[ibin][ibin2][ibin3]->GetNbinsY(); m++){
	    hJetTrackME[ibin][ibin2][ibin3]->SetBinContent(k,m,temp1);
	    hJetTrackME[ibin][ibin2][ibin3]->SetBinError(k,m,err1);
	  }
	}
	
      hJetTrackME[ibin][ibin2][ibin3]->Scale(1./me00);

	//	me_proj[ibin][ibin2][ibin3]->Scale (1./me00);

	me_proj[ibin][ibin2][ibin3]->SetAxisRange(-2.99,2.99);


	
	me_proj[ibin][ibin2][ibin3]->Draw();


	if(ibin<3&&!is_pp){
	  
	  TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	 
	  TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}
	if(ibin==3||is_pp){
	  if(!is_pp){
	    TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	    centtex->SetNDC();
	    centtex->Draw();
	  }
	  TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}
      }
    }
  }
 
  for (int ibin=0;ibin<nCBins;ibin++){

    if(is_pp && ibin > 0 ) continue;
  
    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){ 
 
	yield_inc[ibin][ibin2][ibin3]->Divide(hJetTrackME[ibin][ibin2][ibin3]);
	yield_inc_pTweighted[ibin][ibin2][ibin3]->Divide(hJetTrackME[ibin][ibin2][ibin3]);

     
	if(is_pp) result_proj_eta_canvas->cd(ibin3+1);
	result_proj_eta_canvas->cd(4*(ibin3+1)-ibin);


	llimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(-TMath::Pi()/2+0.0001);
	rlimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(TMath::Pi()/2-0.0001);

	yield_inc_proj_eta[ibin][ibin2][ibin3] = (TH1D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->ProjectionX((TString) ("Yield_Proj_Eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),llimitphi,rlimitphi);

	yield_inc_proj_eta[ibin][ibin2][ibin3]->Scale(1./(rlimiteta-llimiteta+1));
	
	yield_inc_proj_eta[ibin][ibin2][ibin3]->Rebin(10);

	yield_inc_proj_eta[ibin][ibin2][ibin3]->SetAxisRange(-2.999,2.999);
	yield_inc_proj_eta[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	yield_inc_proj_eta[ibin][ibin2][ibin3]->SetMarkerColor(kViolet-2);
	yield_inc_proj_eta[ibin][ibin2][ibin3]->SetLineColor(kViolet-2);
	yield_inc_proj_eta[ibin][ibin2][ibin3]->Draw();



	if(ibin<3&&!is_pp){
	  
	  TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	 
	  TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}
	if(ibin==3||is_pp){
	  if(!is_pp){
	    TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	    centtex->SetNDC();
	    centtex->Draw();
	  }
	  TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}
	//	cout<<"start background code"<<endl;





	if(is_pp) result_proj_canvas->cd(ibin3+1);
	else 	result_proj_canvas->cd(4*(ibin3+1)-ibin);      

	if((data_mc_type_code==2||data_mc_type_code == 12)&&is_pp&&ibin3==0){
	  llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	  rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5-0.0001);

	}else{
	  llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(1.5+0.0001);
	  rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(2.5-0.0001);


	}

	background_proj[ibin][ibin2][ibin3] = (TH1D*)yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("ProjectedBackground%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);

	if(data_mc_type_code == 0&&is_pp&&ibin==0&ibin3==0){
	  llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	  rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5-0.0001);
	}else{
	  llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	  rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5-0.0001);
	}

	background_left[ibin][ibin2][ibin3] = (TH1D*)yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("LeftSideBackground%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);
	      
	background_proj[ibin][ibin2][ibin3]->Add(background_left[ibin][ibin2][ibin3]);

	dx_phi =  background_proj[ibin][ibin2][ibin3]->GetBinWidth(1);

	background_proj[ibin][ibin2][ibin3]->Scale(1./2/(rlimiteta-llimiteta+1));

	  
	background[ibin][ibin2][ibin3] = (TH2D*)yield_inc[ibin][ibin2][ibin3]->Clone((TString) ("SummedBkg_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));



	for(int k = 1;  k<yield_inc[ibin][ibin2][ibin3]->GetNbinsY(); k++){
	  temp1 = background_proj[ibin][ibin2][ibin3]->GetBinContent(k);
	  err1 = background_proj[ibin][ibin2][ibin3]->GetBinError(k);
	    
	  for(int m = 1;  m<yield_inc[ibin][ibin2][ibin3]->GetNbinsX(); m++){
	    background[ibin][ibin2][ibin3]->SetBinContent(m,k,temp1);
	    background[ibin][ibin2][ibin3]->SetBinError(m,k,err1);
	  }
	}

      

	yield_bgsub[ibin][ibin2][ibin3] = (TH2D*)yield_inc[ibin][ibin2][ibin3]->Clone((TString) ("Yield_BkgSub_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	yield_bgsub[ibin][ibin2][ibin3]->Add(background[ibin][ibin2][ibin3],-1.);



	if((data_mc_type_code==2||data_mc_type_code == 12)&&is_pp&&ibin3==0){
	  llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	  rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5-0.0001);

	}else{
	  llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(1.5+0.0001);
	  rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(2.5-0.0001);


	}

	background_proj[ibin][ibin2][ibin3] = (TH1D*)yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("ProjectedBackground%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);

	if(data_mc_type_code == 0&&is_pp&&ibin==0&ibin3==0){
	  llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	  rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5-0.0001);
	}else{
	  llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	  rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5-0.0001);
	}

	background_left[ibin][ibin2][ibin3] = (TH1D*)yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("LeftSideBackground%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);




	if((data_mc_type_code==2||data_mc_type_code == 12)&&ibin<2&&is_pp&&ibin3==0){
	  llimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	  rlimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5-0.0001);

	}else{
	  llimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(1.5+0.0001);
	  rlimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(2.5-0.0001);


	}

	background_proj[ibin][ibin2][ibin3] = (TH1D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->ProjectionY(Form("ProjectedBackground%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);

	if(data_mc_type_code == 0&&is_pp&&ibin==0&ibin3==0){
	  llimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	  rlimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5-0.0001);
	}else{
	  llimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	  rlimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-1.5-0.0001);
	}

	background_left[ibin][ibin2][ibin3] = (TH1D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->ProjectionY(Form("LeftSideBackground%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);







	background_left[ibin][ibin2][ibin3] = (TH1D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->ProjectionY(Form("LeftSideBackground%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);
	      
	background_proj[ibin][ibin2][ibin3]->Add(background_left[ibin][ibin2][ibin3]);

	dx_phi =  background_proj[ibin][ibin2][ibin3]->GetBinWidth(1);

	background_proj[ibin][ibin2][ibin3]->Scale(1./2/(rlimiteta-llimiteta+1));

	  
	background_pTweighted[ibin][ibin2][ibin3] = (TH2D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->Clone((TString) ("SummedBkg_pTweighted_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));



	for(int k = 1;  k<yield_inc[ibin][ibin2][ibin3]->GetNbinsY(); k++){
	  temp1 = background_proj[ibin][ibin2][ibin3]->GetBinContent(k);
	  err1 = background_proj[ibin][ibin2][ibin3]->GetBinError(k);
	    
	  for(int m = 1;  m<yield_inc[ibin][ibin2][ibin3]->GetNbinsX(); m++){
	    background_pTweighted[ibin][ibin2][ibin3]->SetBinContent(m,k,temp1);
	    background_pTweighted[ibin][ibin2][ibin3]->SetBinError(m,k,err1);
	  }
	}

      

	yield_bgsub_pTweighted[ibin][ibin2][ibin3] = (TH2D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->Clone((TString) ("Yield_BkgSub_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	yield_bgsub_pTweighted[ibin][ibin2][ibin3]->Add(background_pTweighted[ibin][ibin2][ibin3],-1.);

	

	if(data_mc_type_code!=0){
	  
	  //	  cout<<"symmetrizing X"<<endl;
  
	  int nbins = yield_bgsub[ibin][ibin2][ibin3]->GetNbinsX();
	  
	  for(int m = 1; m< yield_bgsub[ibin][ibin2][ibin3]->GetNbinsY()+1; m++){
	    for(int k = 1; k<nbins/2+1; k ++){
	      bc = (yield_bgsub[ibin][ibin2][ibin3]->GetBinContent(k,m)+yield_bgsub[ibin][ibin2][ibin3]->GetBinContent(nbins+1-k,m))/2.;
	      err = TMath::Sqrt(yield_bgsub[ibin][ibin2][ibin3]->GetBinError(k,m)*yield_bgsub[ibin][ibin2][ibin3]->GetBinError(k,m)+yield_bgsub[ibin][ibin2][ibin3]->GetBinError(nbins+1-k,m)*yield_bgsub[ibin][ibin2][ibin3]->GetBinError(nbins+1-k,m))/2.;
	      
	      yield_bgsub[ibin][ibin2][ibin3]->SetBinContent(k,m,bc);
	      yield_bgsub[ibin][ibin2][ibin3]->SetBinError(k,m,err);
	      
	      yield_bgsub[ibin][ibin2][ibin3]->SetBinContent(nbins+1-k,m,bc);
	      yield_bgsub[ibin][ibin2][ibin3]->SetBinError(nbins+1-k,m,err);
	      //	      cout<<k<<" "<<nbins+1-k<<" "<<bc<<endl;
	    }
	    
	  }

	  //	  cout<<"symmetrizing Y"<<endl;
	  nbins = yield_bgsub[ibin][ibin2][ibin3]->GetNbinsY()/2;
	  
	  for(int m = 1; m< yield_bgsub[ibin][ibin2][ibin3]->GetNbinsX()+1; m++){
	    for(int k = 1; k<nbins/2+1; k ++){
	      
	      bc = (yield_bgsub[ibin][ibin2][ibin3]->GetBinContent(m,k)+yield_bgsub[ibin][ibin2][ibin3]->GetBinContent(m,nbins+1-k))/2.;
	      err = TMath::Sqrt(yield_bgsub[ibin][ibin2][ibin3]->GetBinError(m,k)*yield_bgsub[ibin][ibin2][ibin3]->GetBinError(m,k)+yield_bgsub[ibin][ibin2][ibin3]->GetBinError(m,nbins+1-k)*yield_bgsub[ibin][ibin2][ibin3]->GetBinError(m,nbins+1-k))/2.;
	      
	      yield_bgsub[ibin][ibin2][ibin3]->SetBinContent(m,k,bc);
	      yield_bgsub[ibin][ibin2][ibin3]->SetBinError(m,k,err);
	      
	      yield_bgsub[ibin][ibin2][ibin3]->SetBinContent(m,nbins+1-k,bc);
	      yield_bgsub[ibin][ibin2][ibin3]->SetBinError(m,nbins+1-k,err);

	      //     cout<<k<<" "<<nbins+1-k<<" "<<bc<<endl;
	      
	    }
	    
	  }

	  //	  cout<<"symmetrizing X"<<endl;
	  nbins = yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetNbinsX();
	  
	  for(int m = 1; m< yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetNbinsY()+1; m++){
	    for(int k = 1; k<nbins/2+1; k ++){
	      bc = (yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinContent(k,m)+yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinContent(nbins+1-k,m))/2.;
	      err = TMath::Sqrt(yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinError(k,m)*yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinError(k,m)+yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinError(nbins+1-k,m)*yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinError(nbins+1-k,m))/2.;
	      
	      yield_bgsub_pTweighted[ibin][ibin2][ibin3]->SetBinContent(k,m,bc);
	      yield_bgsub_pTweighted[ibin][ibin2][ibin3]->SetBinError(k,m,err);
	      
	      yield_bgsub_pTweighted[ibin][ibin2][ibin3]->SetBinContent(nbins+1-k,m,bc);
	      yield_bgsub_pTweighted[ibin][ibin2][ibin3]->SetBinError(nbins+1-k,m,err);
	      //	      cout<<k<<" "<<nbins+1-k<<" "<<bc<<endl;
	    }
	    
	  }

	  //	  cout<<"symmetrizing Y"<<endl;
	  nbins = yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetNbinsY()/2.;
	  
	  for(int m = 1; m< yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetNbinsX()+1; m++){
	    for(int k = 1; k<nbins/2+1; k ++){
	      
	      bc = (yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinContent(m,k)+yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinContent(m,nbins+1-k))/2.;
	      err = TMath::Sqrt(yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinError(m,k)*yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinError(m,k)+yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinError(m,nbins+1-k)*yield_bgsub_pTweighted[ibin][ibin2][ibin3]->GetBinError(m,nbins+1-k))/2.;
	      
	      yield_bgsub_pTweighted[ibin][ibin2][ibin3]->SetBinContent(m,k,bc);
	      yield_bgsub_pTweighted[ibin][ibin2][ibin3]->SetBinError(m,k,err);
	      
	      yield_bgsub_pTweighted[ibin][ibin2][ibin3]->SetBinContent(m,nbins+1-k,bc);
	      yield_bgsub_pTweighted[ibin][ibin2][ibin3]->SetBinError(m,nbins+1-k,err);

	      //     cout<<k<<" "<<nbins+1-k<<" "<<bc<<endl;
	      
	    }
	    
	  }

	}



	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Write();
	hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Write();
	yield_inc[ibin][ibin2][ibin3]->Write();
	yield_inc_pTweighted[ibin][ibin2][ibin3]->Write();
	yield_bgsub[ibin][ibin2][ibin3]->Write();
	yield_bgsub_pTweighted[ibin][ibin2][ibin3]->Write();
	background[ibin][ibin2][ibin3]->Write();
	background_pTweighted[ibin][ibin2][ibin3]->Write();
	hJetTrackME[ibin][ibin2][ibin3]->Write();


	llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(2.5-0.0001);

	yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->ProjectionY((TString) ("Yield_Proj_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),llimiteta,rlimiteta);

	yield_inc_proj[ibin][ibin2][ibin3]->Scale(1./(rlimiteta-llimiteta+1));
	
	yield_inc_proj[ibin][ibin2][ibin3]->Rebin(8);
	yield_inc_proj[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	yield_inc_proj[ibin][ibin2][ibin3]->SetMarkerColor(kViolet);
	yield_inc_proj[ibin][ibin2][ibin3]->SetLineColor(kViolet);
	background_proj[ibin][ibin2][ibin3]->Rebin(8);
	background_proj[ibin][ibin2][ibin3]->SetLineColor(kBlack);
	background_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
	background_proj[ibin][ibin2][ibin3]->SetMarkerStyle(10);

	yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(background_proj[ibin][ibin2][ibin3]->GetMinimum()-0.001);
	yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(background_proj[ibin][ibin2][ibin3]->GetMinimum()+0.01);
	yield_inc_proj[ibin][ibin2][ibin3]->Draw();

	if(ibin<3&&!is_pp){
	  
	  TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	 
	  TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}
	if(ibin==3||is_pp){
	  if(!is_pp){
	    TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	    centtex->SetNDC();
	    centtex->Draw();
	  }
	  TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}	
	background_proj[ibin][ibin2][ibin3]->Draw("same");

  

      }//ibin3;
      
    } // ibin2

  } // ibin ( centrality ) loop

  me_proj_canvas->SaveAs((TString)("All_ME_Projections_Inclusive_"+datalabel+"_"+data_mc_type_strs[data_mc_type_code]+".png"));
  result_proj_canvas->SaveAs((TString)("All_Background_Projections_Inclusive_"+datalabel+"_"+data_mc_type_strs[data_mc_type_code]+".png"));

  result_proj_eta_canvas->SaveAs((TString)("All_dEta_Projections_Inclusive_"+datalabel+"_"+data_mc_type_strs[data_mc_type_code]+".png"));

  if(data_mc_type_code==0){
 
    result_proj_canvas->SaveAs((TString)("All_Background_Projections_Inclusive_"+datalabel+"_"+data_mc_type_strs[data_mc_type_code]+".pdf"));

    result_proj_eta_canvas->SaveAs((TString)("All_dEta_Projections_Inclusive_"+datalabel+"_"+data_mc_type_strs[data_mc_type_code]+".pdf"));

  }
    
  return 0;
} // main loop






  
     
