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


Int_t me_correct2(int data_mc_type_code=0, bool is_pp=kFALSE, bool do_residual = kFALSE, bool signal_only_me = kFALSE){

  //-------------Input files & set eta-----------------

  TString datalabel;
  float eta_ymax;

  TLatex *centtex, *pttex;

  TF1 *pol0 = new TF1("pol0","[0]+x-x",-.3,.3);

  TF1 *abs_fit = new TF1("abs_fit","[0]*x^2+[1]",-3.,3.);
 

  TFile *fin, *fin2, *fin_me, *fout_inc;
  enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RightGen, SpilledUnderGen, UnmatchedGen, RightReco, SpilledReco, UnmatchedReco, RecoGenSube0,RecoGenNoSube0,GenGenSube0,GenGenNoSube0,MatchedRecoGenSube0,MatchedRecoGenNoSube0,SwappedRecoGenSube0,SwappedRecoGenNoSube0, UnMatchedRecoGenSube0,UnMatchedRecoGenNoSube0,n_data_mc_types};


  TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RightGenJet_GenTrack","SpilledUnderJet_GenTrack","UnmatchedGenJet_GenTrack","RightRecoJet_GenTrack","SpilledReco_GenTrack","UnmatchedReco_GenTrack","RecoJet_GenTrack_Sube0","RecoJet_GenTrack_NoSube0","GenJet_GenTrack_Sube0","GenJet_GenTrack_NoSube0","MatchedRecoJet_GenTrack_Sube0","MatchedRecoJet_GenTrack_NoSube0","SwappedRecoJet_GenTrack_Sube0","SwappedRecoJet_GenTrack_NoSube0","UnmatchedRecoJet_GenTrack_Sube0","UnmatchedRecoJet_GenTrack_NoSube0",};

  TString desc = data_mc_type_strs[data_mc_type_code];

  desc.ReplaceAll("_Sube0","");
  desc.ReplaceAll("_NoSube0","");


  if(data_mc_type_code==0){
    if(is_pp){
      datalabel = "pp";
      //fin = new TFile("../data_raw_correlations/pp_Data_withMix_inclJetBinning_finalJFFs.root","READ");
      //fin_me = new TFile("../data_raw_correlations/pp_Data_withMix_inclJetBinning_finalJFFs.root","READ");
      fin = new TFile("../data_raw_correlations/pp_Data_withMix_inclJetBinning_finalJFFs_60vzMixBins.root","READ");
      fin_me = new TFile("../data_raw_correlations/pp_Data_withMix_inclJetBinning_finalJFFs_60vzMixBins.root","READ");
      fout_inc = new TFile("pp_Inclusive_Correlations.root","RECREATE");
    }else{
      datalabel = "PbPb"; 
      //fin = new TFile("../data_raw_correlations/PbPbData_noMix_inclJetBinning_finalJFFs_spillOutJets_eta0p5.root","READ");
      fin = new TFile("../data_raw_correlations/PbPb_Data_withMix_CymbalTune_fullCymbalCorrs_inclJetBinning_withPtWeightME.root","READ");
      //      fin_me = new TFile("../data_raw_correlations/PbPb_Data_withMix_CymbalTune_fullCymbalCorrs_mixEvtWingCorrected_absSinUnweight_inclJetBinning.root","READ");

      fin2 = new TFile(" ../data_raw_correlations/PbPbData_noMix_inclJetBinning_finalJFFs_officialTrkCorrs.root","READ");
      fin_me = new TFile("../data_raw_correlations/PbPb_Data_withMix_CymbalTune_fullCymbalCorrs_inclJetBinning_withPtWeightME.root","READ");
      fout_inc = new TFile("PbPb_Inclusive_Correlations.root","RECREATE");
    }


  }else{
    if(is_pp){
      datalabel = "Pythia";
      switch(data_mc_type_code){
  
	 case 1:
	fin = new TFile("../mc_raw_correlations/pp_Pythia6MC_RecoReco_withMix_inclJetBinning_finalJFFs.root","READ");
	fin_me = new TFile("../mc_raw_correlations/pp_Pythia6MC_RecoReco_withMix_inclJetBinning_finalJFFs.root","READ");
	break;
      case 2:
	fin = new TFile("../mc_raw_correlations/pp_Pythia6MC_RecoGen_withMix_inclJetBinning_finalJFFs_pfID1.root","READ");
	//	fin_me = new TFile("../mc_raw_correlations/pp_5TeV_MC_Pythia_MixHistos_Merged_RecoGenReduced_fixVz_fixTrkCorr_fineBin.root","READ");
	fin_me = new TFile("../mc_raw_correlations/pp_Pythia6MC_RecoGen_withMix_inclJetBinning_finalJFFs_pfID1.root","READ");
	
	break;
      case 4:
	fin = new TFile("../mc_raw_correlations/pp_Pythia6MC_GenGen_withMix_inclJetBinning_finalJFFs.root","READ");
	fin_me = new TFile("../mc_raw_correlations/pp_Pythia6MC_GenGen_withMix_inclJetBinning_finalJFFs.root","READ");
	break;
      }
	fout_inc = new TFile((TString)("Pythia_"+data_mc_type_strs[data_mc_type_code]+"_Inclusive_Correlations.root"),"RECREATE");      
    }else{
      datalabel = "Hydjet";
      switch(data_mc_type_code){
      case 1: 
	fin = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_RecoReco_CymbalTune_inclJetBinning.root","READ");
	fin_me = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_RecoReco_CymbalTune_inclJetBinning.root","READ");
	break;
      case 2: 
	fin = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_RecoGen_CymbalTune_inclJetBinning.root","READ");
	fin_me = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_RecoGen_CymbalTune_inclJetBinning.root","READ");
	break;
      case 3: 
	fin = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_GenReco_CymbalTune_inclJetBinning.root","READ");
	fin_me = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_GenReco_CymbalTune_inclJetBinning.root","READ");
	break;
      case 4: 
	fin = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_GenGen_CymbalTune_inclJetBinning.root","READ");
	fin_me = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_GenGen_CymbalTune_inclJetBinning.root","READ");

	break;

      case 11: 
	fin = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_RecoGen_sube0_CymbalTune_fullCymbalCorrs_inclJetBinning.root","READ");
	//	fin_me = new TFile("../mc_raw_correlations/PbPb_5TeV_MC_PythiaHydjet_MixHistos_RecoGenReduced_Merged_refpt_newJetTrackCorrections_fineBin.root","READ");
	fin_me = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_RecoGen_sube0_CymbalTune_fullCymbalCorrs_inclJetBinning.root","READ");
	break;

      case 12: 
	fin = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_RecoGen_subeNon0_CymbalTune_fullCymbalCorrs_inclJetBinning.root","READ");
	fin_me = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_RecoGen_subeNon0_CymbalTune_fullCymbalCorrs_inclJetBinning.root","READ");
	break;

      case 13: 
	fin = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_GenGen_sube0_CymbalTune_fullCymbalCorrs_inclJetBinning.root","READ");
	//	fin_me = new TFile("../mc_raw_correlations/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_GenGenReduced_fineBin.root","READ");
	fin_me = new TFile("../mc_raw_correlations/PbPb_5TeVMC_withMix_GenGen_sube0_CymbalTune_fullCymbalCorrs_inclJetBinning.root","READ");
	break;

      case 14: 
	fin = new TFile("../mc_raw_correlations/PbPb_5TeV_MC_PythiaHydjet_MixHistos_GenGenReduced_subeNon0_Merged_withMix_newJetTrackCorrections_fineBin.root","READ");
	fin_me = new TFile("../mc_raw_correlations/PbPb_5TeV_MC_PythiaHydjet_MixHistos_GenGenReduced_subeNon0_Merged_withMix_newJetTrackCorrections_fineBin.root","READ");
	break;

	
      default: 
	cout<<"Sorry, you are mistaken, this is not a sample we have"<<endl;
	return -1;
      }

      fout_inc = new TFile((TString)("Hydjet_"+data_mc_type_strs[data_mc_type_code]+"_Inclusive_Correlations.root"),"RECREATE");      
       }
    }
    //----------------------------------------------------



  TCanvas *me_proj_canvas, *result_proj_canvas, *result_proj_eta_canvas;

  TCanvas *me_proj_canvas2, *result_proj_canvas2, *result_proj_eta_canvas2;

  if(is_pp){

    me_proj_canvas = new TCanvas("me_proj_canvas","",0,0,400,2400);
    me_proj_canvas->Divide(1,9,0.0000,0.0000);


    result_proj_canvas = new TCanvas("result_proj_canvas","",0,0,400,2400);
    result_proj_canvas->Divide(1,9,0.0000,0.0000);

    result_proj_eta_canvas = new TCanvas("result_proj_canvas_eta","",0,0,400,2400);
    result_proj_eta_canvas->Divide(1,9,0.0000,0.0000);

    me_proj_canvas2 = new TCanvas("me_proj_canvas2","",0,0,400,2400);
    me_proj_canvas2->Divide(1,9,0.0000,0.0000);


    result_proj_canvas2 = new TCanvas("result_proj_canvas2","",0,0,400,2400);
    result_proj_canvas2->Divide(1,9,0.0000,0.0000);

    result_proj_eta_canvas2 = new TCanvas("result_proj_eta_canvas2","",0,0,400,2400);
    result_proj_eta_canvas2->Divide(1,9,0.0000,0.0000);


  }else{

    me_proj_canvas = new TCanvas("me_proj_canvas","",0,0,1600,2400);
    me_proj_canvas->Divide(4,9,0.0000,0.0000);


    result_proj_canvas = new TCanvas("result_proj_canvas","",0,0,1600,2400);
    result_proj_canvas->Divide(4,9,0.0000,0.0000);

    result_proj_eta_canvas = new TCanvas("result_proj_canvas_eta","",0,0,1600,2400);
    result_proj_eta_canvas->Divide(4,9,0.0000,0.0000);

    me_proj_canvas2 = new TCanvas("me_proj_canvas2","",0,0,1600,2400);
    me_proj_canvas2->Divide(4,9,0.0000,0.0000);


    result_proj_canvas2 = new TCanvas("result_proj_canvas2","",0,0,1600,2400);
    result_proj_canvas2->Divide(4,9,0.0000,0.0000);

    result_proj_eta_canvas2 = new TCanvas("result_proj_canvas_eta2","",0,0,1600,2400);
    result_proj_eta_canvas2->Divide(4,9,0.0000,0.0000);

  }

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins;
  
  const int nCBins = 5;
  const int nPtBins = 1;
  const int nTrkPtBins = 9;

  float me00,me00_lead,me00_sub;
  int me00binl, me00binr, n_me00bins;

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt1000"};
  
  /*
  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};
  */

  float CBins[nCBins+1] = {0, 20, 60, 100, 140, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent70", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent.50-70%","Cent. 70-100%"};


  float TrkPtBins[nTrkPtBins+1] = {07, 1, 2, 3, 4, 8, 12, 16, 20, 999};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_strs3[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_strs2[nTrkPtBins+1] = {"TrkPt0p7","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt999" };
   TString TrkPtBin_labels[nTrkPtBins] = {"0.7<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};

  TString ref_name;
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
  
  TH1F* all_jets_corrpT[nCBins][nPtBins];
 
  TH2D* hJetTrackSignalBackground[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackground2[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackground_pTweighted[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackground_pTweighted2[nCBins][nPtBins][nTrkPtBins];
  
  TH2D* hJetTrackME[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackME2[nCBins][nPtBins][nTrkPtBins];
 
  TH2D *yield_inc[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_inc_pTweighted[nCBins][nPtBins][nTrkPtBins];

  TH2D *yield_bgsub[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_bgsub_pTweighted[nCBins][nPtBins][nTrkPtBins];
 
 
  TH1D *yield_inc_proj[nCBins][nPtBins][nTrkPtBins];
  TH1D *yield_inc_proj2[nCBins][nPtBins][nTrkPtBins];
  TH1D *yield_inc_proj_eta[nCBins][nPtBins][nTrkPtBins];
  TH1D *yield_inc_proj_eta2[nCBins][nPtBins][nTrkPtBins];
  TH1D *test_me_inc_proj[nCBins][nPtBins][nTrkPtBins];
 
  TH1D *me_proj[nCBins][nPtBins][nTrkPtBins];
  TH1D *me_proj2[nCBins][nPtBins][nTrkPtBins];
  TH1D *me_proj3[nCBins][nPtBins][nTrkPtBins];

  TH2D* background[nCBins][nPtBins][nTrkPtBins];
  TH2D* background_pTweighted[nCBins][nPtBins][nTrkPtBins];
  TH1D* background_left[nCBins][nPtBins][nTrkPtBins];
  TH1D* background_right[nCBins][nPtBins][nTrkPtBins];
  TH1D* background_proj[nCBins][nPtBins][nTrkPtBins];
  TH1D* background_proj2[nCBins][nPtBins][nTrkPtBins];

  TH1D* eta_side_proj[nCBins][nPtBins][nTrkPtBins];
  TH1D* eta_side_proj2[nCBins][nPtBins][nTrkPtBins];

  float norm_temp, width_temp, width_temp_x, width_temp_y, width_temp_ref, max_bin, max_cont,bc,err,err1, dx_phi, dx_eta, temp1, pt_x;
    
  //-----------------------
  // Start getting histos
  //-----------------------

  
  for (int ibin=0;ibin<nCBins;ibin++){

    if(is_pp && ibin > 0 ) continue;
  
    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 

      cout<<desc + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]<<endl;

      all_jets_corrpT[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc+ "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));
         
      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

	if(data_mc_type_code > 0 && (data_mc_type_code%2==0||data_mc_type_code > 10)){ //if gen
	  cout<<"1"<<endl;	
	  hJetTrackSignalBackground[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackSignalBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  cout<<"2"<<endl;
	  hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackSignalBackground_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Raw_Yield_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  cout<<"3"<<endl;
	  hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)(desc+"_hJetTrackME_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  hJetTrackME2[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)(desc+"_hJetTrackME_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Mixed_Event2_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  cout<<"4"<<endl;
	  
	}else{ //this is Reco
	  cout<<"1"<<endl;
	  if(data_mc_type_code==0&&!is_pp&&ibin3<3&&ibin<2)  hJetTrackSignalBackground[ibin][ibin2][ibin3] = (TH2D*) fin2->Get((TString)(desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  else  hJetTrackSignalBackground[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  cout<<"2"<<endl;
	  if(data_mc_type_code==0&&!is_pp&&ibin3<3&&ibin<2) hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3] = (TH2D*) fin2->Get((TString)(desc + "_hJetTrackSignalBackground_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Raw_Yield_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  else	  hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackSignalBackground_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Raw_Yield_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  cout<<"3"<<endl;
	  //	  hJetTrackME2[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc+"_hJetTrackSignalBackground_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Mixed_Event2_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  //	  hJetTrackME3[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc+"_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Mixed_Event3_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  if(signal_only_me||(data_mc_type_code==0&&!is_pp&&ibin3<3&&ibin<2)){

	    hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin2->Get((TString)(desc+"_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    hJetTrackME2[ibin][ibin2][ibin3] = (TH2D*) fin2->Get((TString)(desc+"_hJetTrackSignalBackground_pTweighted"+ CBin_strs[ibin] +"_"+CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Mixed_Event2_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  }else{
	    // if(data_mc_type_code==0&&!is_pp&&ibin3<3&&ibin<2)	hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc+"_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    //else
	    hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)(desc+"_hJetTrackME"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    cout<<"4"<<endl;

	    //	    if(data_mc_type_code==0&&!is_pp&&ibin3<3&&ibin<2) hJetTrackME2[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)(desc+"_hJetTrackSignalBackground_pTweighted"+ CBin_strs[ibin] +"_"+CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Mixed_Event2_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
 
	    //else 
	    if(data_mc_type_code==0&&!is_pp)	  hJetTrackME2[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)(desc+"_hJetTrackME_pTweighted"+ CBin_strs[ibin] +"_"+ CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Mixed_Event2_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    else hJetTrackME2[ibin][ibin2][ibin3] = (TH2D*) fin_me->Get((TString)(desc+"_hJetTrackME"+ CBin_strs[ibin] +"_"+ CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs2[ibin3+1]))->Clone((TString) ("Mixed_Event2_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    cout<<"5"<<endl;

	  }
	}
      } /// ibin3

    } // ibin2

    //-------------------------
    // Normalize and subtract bkg
    //------------------------

  } //ibin
    

  for (int ibin=0;ibin<nCBins-1;ibin++){

    if(ibin > 0 && is_pp) continue;

    cout<<"got hists for centbin "<<ibin<<endl;
   
    for (int ibin2=0;ibin2<nPtBins;ibin2++){

      if(ibin==3){
	 
	all_jets_corrpT[ibin][ibin2]->Add(all_jets_corrpT[ibin+1][ibin2]);

	TString temp_name = all_jets_corrpT[ibin][ibin2]->GetName();
	temp_name.ReplaceAll("Cent70","Cent100");
	all_jets_corrpT[ibin][ibin2]->SetName( temp_name);
      }


      for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){

	if(ibin==3){

	  cout<<"Combining bins: "<<all_jets_corrpT[ibin][ibin2]->GetEntries()<<" "<<all_jets_corrpT[ibin+1][ibin2]->GetEntries()<<" "<< hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetEntries()<<" "<< hJetTrackSignalBackground[ibin+1][ibin2][ibin3]->GetEntries()<<endl;

	  hJetTrackSignalBackground[ibin][ibin2][ibin3]->Add(hJetTrackSignalBackground[ibin+1][ibin2][ibin3]);
	  hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Add(hJetTrackSignalBackground_pTweighted[ibin+1][ibin2][ibin3]);
	  // hJetTrackME[ibin][ibin2][ibin3]->Add(hJetTrackME[ibin+1][ibin2][ibin3]);
	  //	  hJetTrackME2[ibin][ibin2][ibin3]->Add(hJetTrackME2[ibin+1][ibin2][ibin3]);
	 

	 TString temp_name =  hJetTrackSignalBackground[ibin][ibin2][ibin3]->GetName();
	  temp_name.ReplaceAll("Cent70","Cent100");
	  hJetTrackSignalBackground[ibin][ibin2][ibin3]->SetName(temp_name);

	  temp_name =  hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->GetName();
	  temp_name.ReplaceAll("Cent70","Cent100");
	  hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->SetName(temp_name);

	  temp_name =  hJetTrackME[ibin][ibin2][ibin3]->GetName();
	  temp_name.ReplaceAll("Cent70","Cent100");
	  hJetTrackME[ibin][ibin2][ibin3]->SetName(temp_name);


	  temp_name =  hJetTrackME2[ibin][ibin2][ibin3]->GetName();
	  temp_name.ReplaceAll("Cent70","Cent100");
	  hJetTrackME2[ibin][ibin2][ibin3]->SetName(temp_name);

	}
	
	
	width_temp_x = 1.;
	width_temp_y = 1.;
	
	norm_temp = all_jets_corrpT[ibin][ibin2]->Integral();
	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);

	hJetTrackME[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y/50.);

	hJetTrackME2[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y/50.);
   
	//     INCLUSIVE

	fout_inc->cd();
	
	if(is_pp)  	me_proj_canvas->cd(ibin3+1);
	else me_proj_canvas->cd(4*(ibin3+1)-ibin);


	if(ibin==3)	yield_inc[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackground[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+2]+ "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	else 	yield_inc[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackground[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+  "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	if(ibin==3)	yield_inc_pTweighted[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Clone((TString)("Yield_pTweighted_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+2]+  "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	else 	yield_inc_pTweighted[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Clone((TString)("Yield_pTweighted_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+  "_Pt100_Pt1000_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	int lbin = 1;
	int rbin = 200;


	//	if(data_mc_type_code==0&&!is_pp&&ibin3<3&&ibin<2){
	if((data_mc_type_code==0&&!is_pp&&ibin3<3&&ibin<2)||signal_only_me){
	  lbin = 90;
	  //  rbin = 11;
	}

	me_proj[ibin][ibin2][ibin3] = hJetTrackME[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),lbin,rbin);
	
	me_proj[ibin][ibin2][ibin3]->Scale (1./(rbin-lbin+1));
	me_proj[ibin][ibin2][ibin3] -> Fit("pol0","q","",-.3,.3);

	me00 = 	me_proj[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);


	
	for(int k = 1;  k<yield_inc[ibin][ibin2][ibin3]->GetNbinsX()+1; k++){
	  temp1 = me_proj[ibin][ibin2][ibin3]->GetBinContent(k);
	  err1 = me_proj[ibin][ibin2][ibin3]->GetBinError(k);
	    
	  for(int m = 1;  m<yield_inc[ibin][ibin2][ibin3]->GetNbinsY()+1; m++){
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


	if(is_pp)  	me_proj_canvas2->cd(ibin3+1);
	else me_proj_canvas2->cd(4*(ibin3+1)-ibin);

	me_proj2[ibin][ibin2][ibin3] = hJetTrackME2[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_temp2_%d%d%d",ibin,ibin2,ibin3),lbin,rbin);
	
	me_proj2[ibin][ibin2][ibin3]->Scale (1./(rbin-lbin+1));
	me_proj2[ibin][ibin2][ibin3] -> Fit("pol0","q","",-.3,.3);

	me00 = 	me_proj2[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);


	
	for(int k = 1;  k<yield_inc[ibin][ibin2][ibin3]->GetNbinsX()+1; k++){
	  temp1 = me_proj2[ibin][ibin2][ibin3]->GetBinContent(k);
	  err1 = me_proj2[ibin][ibin2][ibin3]->GetBinError(k);
	    
	  for(int m = 1;  m<yield_inc[ibin][ibin2][ibin3]->GetNbinsY()+1; m++){
	    hJetTrackME2[ibin][ibin2][ibin3]->SetBinContent(k,m,temp1);
	    hJetTrackME2[ibin][ibin2][ibin3]->SetBinError(k,m,err1);
	  }
	}
	
      hJetTrackME2[ibin][ibin2][ibin3]->Scale(1./me00);

	//	me_proj[ibin][ibin2][ibin3]->Scale (1./me00);

	me_proj2[ibin][ibin2][ibin3]->SetAxisRange(-2.99,2.99);


	
	me_proj2[ibin][ibin2][ibin3]->Draw();

  
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

 
  for (int ibin=0;ibin<nCBins-1;ibin++){

    if(is_pp && ibin > 0 ) continue;
  
    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){ 
 
	  yield_inc_pTweighted[ibin][ibin2][ibin3]->Divide(hJetTrackME2[ibin][ibin2][ibin3]);

	  yield_inc[ibin][ibin2][ibin3]->Divide(hJetTrackME[ibin][ibin2][ibin3]);
	
	
     
	if(is_pp) result_proj_eta_canvas->cd(ibin3+1);
	result_proj_eta_canvas->cd(4*(ibin3+1)-ibin);


	if(do_residual){
	  llimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(1.5+.001);
	  rlimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(2.-.001);


	  eta_side_proj[ibin][ibin2][ibin3] = (TH1D*)yield_inc[ibin][ibin2][ibin3]->ProjectionX((TString) ("Eta_Sideband"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),llimitphi,rlimitphi);
      
	  eta_side_proj[ibin][ibin2][ibin3]->Scale(1./(rlimitphi-llimitphi+1));

	  eta_side_proj[ibin][ibin2][ibin3]->Rebin(10);
	  eta_side_proj[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	  eta_side_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
	  eta_side_proj[ibin][ibin2][ibin3]->SetLineColor(kBlack);
	  eta_side_proj[ibin][ibin2][ibin3]->SetAxisRange(-3.,3.);
	  if(ibin3 < 2)	eta_side_proj[ibin][ibin2][ibin3]->SetMaximum(	eta_side_proj[ibin][ibin2][ibin3]->GetMaximum()+.005);
	  else 	eta_side_proj[ibin][ibin2][ibin3]->SetMaximum(	eta_side_proj[ibin][ibin2][ibin3]->GetMaximum()+.002);
	  eta_side_proj[ibin][ibin2][ibin3]->SetMinimum(	eta_side_proj[ibin][ibin2][ibin3]->GetMinimum()-.001);
	  eta_side_proj[ibin][ibin2][ibin3]->Draw();

	  eta_side_proj[ibin][ibin2][ibin3]->Fit("abs_fit","","",-3.,3.);


	  if(ibin3<6){
	    for(int m = 1; m< yield_inc[ibin][ibin2][ibin3]->GetNbinsY()+1; m++){
	      for(int k = 1; k< yield_inc[ibin][ibin2][ibin3]->GetNbinsX()+1; k ++){

		pt_x = 	yield_inc[ibin][ibin2][ibin3]->GetXaxis()->GetBinCenter(k);
	    
		norm_temp = abs_fit->Eval(pt_x)/eta_side_proj[ibin][ibin2][ibin3]->GetBinContent(eta_side_proj[ibin][ibin2][ibin3]->FindBin(0.));

		bc =  yield_inc[ibin][ibin2][ibin3]->GetBinContent(k,m);
		err =  yield_inc[ibin][ibin2][ibin3]->GetBinError(k,m);
	      
		yield_inc[ibin][ibin2][ibin3]->SetBinContent(k,m,bc/norm_temp);
		yield_inc[ibin][ibin2][ibin3]->SetBinError(k,m,err/norm_temp);	 

	      }
	    
	    }
	  }
	}

	llimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(1.5+.001);
	rlimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(2.-.001);


	eta_side_proj[ibin][ibin2][ibin3] = (TH1D*)yield_inc[ibin][ibin2][ibin3]->ProjectionX((TString) ("Eta_Sideband"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),llimitphi,rlimitphi);
      
	eta_side_proj[ibin][ibin2][ibin3]->Scale(1./(rlimitphi-llimitphi+1));

	eta_side_proj[ibin][ibin2][ibin3]->Rebin(10);
	eta_side_proj[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	eta_side_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
	eta_side_proj[ibin][ibin2][ibin3]->SetLineColor(kBlack);
	eta_side_proj[ibin][ibin2][ibin3]->SetAxisRange(-3.,3.);
	if(ibin3 < 2)	eta_side_proj[ibin][ibin2][ibin3]->SetMaximum(	eta_side_proj[ibin][ibin2][ibin3]->GetMaximum()+.005);
	   else 	eta_side_proj[ibin][ibin2][ibin3]->SetMaximum(	eta_side_proj[ibin][ibin2][ibin3]->GetMaximum()+.002);
	eta_side_proj[ibin][ibin2][ibin3]->SetMinimum(	eta_side_proj[ibin][ibin2][ibin3]->GetMinimum()-.001);
	eta_side_proj[ibin][ibin2][ibin3]->Draw("same");


	eta_side_proj[ibin][ibin2][ibin3]->Fit("abs_fit","","",-3.,3.);


	llimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(-TMath::Pi()/2+0.0001);
	rlimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(TMath::Pi()/2-0.0001);

	yield_inc_proj_eta[ibin][ibin2][ibin3] = (TH1D*)yield_inc[ibin][ibin2][ibin3]->ProjectionX((TString) ("Yield_Proj_Eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),llimitphi,rlimitphi);

	yield_inc_proj_eta[ibin][ibin2][ibin3]->Scale(1./(rlimitphi-llimitphi+1));
	
	yield_inc_proj_eta[ibin][ibin2][ibin3]->Rebin(10);

	yield_inc_proj_eta[ibin][ibin2][ibin3]->SetAxisRange(-2.999,2.999);
	yield_inc_proj_eta[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	yield_inc_proj_eta[ibin][ibin2][ibin3]->SetMarkerColor(kViolet-2);
	yield_inc_proj_eta[ibin][ibin2][ibin3]->SetLineColor(kViolet-2);
	yield_inc_proj_eta[ibin][ibin2][ibin3]->SetMaximum(yield_inc_proj_eta[ibin][ibin2][ibin3]->GetMaximum()+0.003);
	yield_inc_proj_eta[ibin][ibin2][ibin3]->SetMinimum(yield_inc_proj_eta[ibin][ibin2][ibin3]->GetMinimum()-0.001);
	if(ibin3<3)	yield_inc_proj_eta[ibin][ibin2][ibin3]->SetMinimum(yield_inc_proj_eta[ibin][ibin2][ibin3]->GetMinimum()-0.005);
	yield_inc_proj_eta[ibin][ibin2][ibin3]->Draw("same");





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


	if(is_pp) result_proj_eta_canvas2->cd(ibin3+1);
	result_proj_eta_canvas2->cd(4*(ibin3+1)-ibin);


	if(do_residual){

	  llimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(1.5+0.0001);
	  rlimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(2.-0.0001);


	  eta_side_proj2[ibin][ibin2][ibin3] = (TH1D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->ProjectionX((TString) ("Eta_Sideband2"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),llimitphi,rlimitphi);
      
	  eta_side_proj2[ibin][ibin2][ibin3]->Scale(1./(rlimitphi-llimitphi+1));

	  eta_side_proj2[ibin][ibin2][ibin3]->Rebin(10);


	  eta_side_proj2[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	  eta_side_proj2[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
	  eta_side_proj2[ibin][ibin2][ibin3]->SetLineColor(kBlack);
	  eta_side_proj2[ibin][ibin2][ibin3]->SetAxisRange(-3.,3.);
	  if(ibin3 < 2) eta_side_proj2[ibin][ibin2][ibin3]->SetMaximum(	eta_side_proj2[ibin][ibin2][ibin3]->GetMaximum()+.007);
	  else	eta_side_proj2[ibin][ibin2][ibin3]->SetMaximum(	eta_side_proj2[ibin][ibin2][ibin3]->GetMaximum()+.002);
	  eta_side_proj2[ibin][ibin2][ibin3]->SetMinimum(	eta_side_proj2[ibin][ibin2][ibin3]->GetMinimum()-.001);
	  eta_side_proj2[ibin][ibin2][ibin3]->Draw();

	  //	abs_fit->SetParameter(0,0.003);
	  //	abs_fit->SetParameter(1,0.03);
	  //	abs_fit->SetParameter(2,	eta_side_proj2[ibin][ibin2][ibin3]->GetBinContent(	eta_side_proj2[ibin][ibin2][ibin3]->FindBin(0.)));
       
	  eta_side_proj2[ibin][ibin2][ibin3]->Fit("abs_fit","","",-3.,3.);

	  if(ibin3<6){
	    for(int m = 1; m< yield_inc_pTweighted[ibin][ibin2][ibin3]->GetNbinsY()+1; m++){
	      for(int k = 1; k < yield_inc_pTweighted[ibin][ibin2][ibin3]->GetNbinsX()+1; k ++){


		pt_x = 	yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->GetBinCenter(k);

		norm_temp = abs_fit->Eval(pt_x)/eta_side_proj2[ibin][ibin2][ibin3]->GetBinContent(eta_side_proj2[ibin][ibin2][ibin3]->FindBin(0.));
		bc =  yield_inc_pTweighted[ibin][ibin2][ibin3]->GetBinContent(k,m);
		err =  yield_inc_pTweighted[ibin][ibin2][ibin3]->GetBinError(k,m);
	      
		yield_inc_pTweighted[ibin][ibin2][ibin3]->SetBinContent(k,m,bc/norm_temp);
		yield_inc_pTweighted[ibin][ibin2][ibin3]->SetBinError(k,m,err/norm_temp);	 

	      }
	    
	    }
	  }
	}
	
	llimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(1.5+0.0001);
	rlimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(2.-0.0001);


	eta_side_proj2[ibin][ibin2][ibin3] = (TH1D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->ProjectionX((TString) ("Eta_Sideband2"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),llimitphi,rlimitphi);
      
	eta_side_proj2[ibin][ibin2][ibin3]->Scale(1./(rlimitphi-llimitphi+1));

	eta_side_proj2[ibin][ibin2][ibin3]->Rebin(10);


	eta_side_proj2[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	eta_side_proj2[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
	eta_side_proj2[ibin][ibin2][ibin3]->SetLineColor(kBlack);
	eta_side_proj2[ibin][ibin2][ibin3]->SetAxisRange(-3.,3.);
	if(ibin3 < 2) eta_side_proj2[ibin][ibin2][ibin3]->SetMaximum(	eta_side_proj2[ibin][ibin2][ibin3]->GetMaximum()+.007);
	else	eta_side_proj2[ibin][ibin2][ibin3]->SetMaximum(	eta_side_proj2[ibin][ibin2][ibin3]->GetMaximum()+.002);
	eta_side_proj2[ibin][ibin2][ibin3]->SetMinimum(	eta_side_proj2[ibin][ibin2][ibin3]->GetMinimum()-.001);
	eta_side_proj2[ibin][ibin2][ibin3]->Draw();

	eta_side_proj2[ibin][ibin2][ibin3]->Fit("abs_fit","","",-3.,3.);


	llimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(-TMath::Pi()/2+0.0001);
	rlimitphi = yield_inc[ibin][ibin2][ibin3]->GetYaxis()->FindBin(TMath::Pi()/2-0.0001);

	yield_inc_proj_eta2[ibin][ibin2][ibin3] = (TH1D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->ProjectionX((TString) ("Yield_Proj_Eta2"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),llimitphi,rlimitphi);

	yield_inc_proj_eta2[ibin][ibin2][ibin3]->Scale(1./(rlimitphi-llimitphi+1));
	
	yield_inc_proj_eta2[ibin][ibin2][ibin3]->Rebin(10);

	yield_inc_proj_eta2[ibin][ibin2][ibin3]->SetAxisRange(-2.999,2.999);
	yield_inc_proj_eta2[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	yield_inc_proj_eta2[ibin][ibin2][ibin3]->SetMarkerColor(kViolet-2);
	yield_inc_proj_eta2[ibin][ibin2][ibin3]->SetLineColor(kViolet-2);

	yield_inc_proj_eta2[ibin][ibin2][ibin3]->SetMaximum(yield_inc_proj_eta2[ibin][ibin2][ibin3]->GetMinimum()+0.002);
	yield_inc_proj_eta2[ibin][ibin2][ibin3]->SetMinimum(yield_inc_proj_eta2[ibin][ibin2][ibin3]->GetMinimum()-0.001);
	if(ibin3<3)	yield_inc_proj_eta2[ibin][ibin2][ibin3]->SetMinimum(yield_inc_proj_eta2[ibin][ibin2][ibin3]->GetMinimum()-0.005);
	yield_inc_proj_eta2[ibin][ibin2][ibin3]->Draw("same");



      
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


	float bg_inner = 1.5;
	float bg_outer = 2.5;

	if(is_pp){
	  bg_inner = 1.5;
	  bg_outer = 2.5;
	  
	}
	


	llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(bg_inner+0.0001);
	rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(bg_outer-0.0001);


	if(data_mc_type_code==0&&!is_pp&&ibin==0&&ibin3==1){

	  bg_outer = 3.0;

	llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-bg_outer+0.0001);
	rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-bg_inner-0.0001);

	}


	background_proj[ibin][ibin2][ibin3] = (TH1D*)yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("ProjectedBackground%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);


	llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-bg_outer+0.0001);
	rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-bg_inner-0.0001);

	
	background_left[ibin][ibin2][ibin3] = (TH1D*)yield_inc[ibin][ibin2][ibin3]->ProjectionY(Form("LeftSideBackground%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);

	      
	background_proj[ibin][ibin2][ibin3]->Add(background_left[ibin][ibin2][ibin3]);

	dx_phi =  background_proj[ibin][ibin2][ibin3]->GetBinWidth(1);

	background_proj[ibin][ibin2][ibin3]->Scale(1./2/(rlimiteta-llimiteta+1));

	  
	if(ibin==3)	background[ibin][ibin2][ibin3] = (TH2D*)yield_inc[ibin][ibin2][ibin3]->Clone((TString) ("SummedBkg_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+2] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	else 	background[ibin][ibin2][ibin3] = (TH2D*)yield_inc[ibin][ibin2][ibin3]->Clone((TString) ("SummedBkg_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));



	for(int k = 1;  k<yield_inc[ibin][ibin2][ibin3]->GetNbinsY(); k++){
	  temp1 = background_proj[ibin][ibin2][ibin3]->GetBinContent(k);
	  err1 = background_proj[ibin][ibin2][ibin3]->GetBinError(k);
	    
	  for(int m = 1;  m<yield_inc[ibin][ibin2][ibin3]->GetNbinsX(); m++){
	    background[ibin][ibin2][ibin3]->SetBinContent(m,k,temp1);
	    background[ibin][ibin2][ibin3]->SetBinError(m,k,err1);
	  }
	}

      

	if(ibin==3)	yield_bgsub[ibin][ibin2][ibin3] = (TH2D*)yield_inc[ibin][ibin2][ibin3]->Clone((TString) ("Yield_BkgSub_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+2] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	else 	yield_bgsub[ibin][ibin2][ibin3] = (TH2D*)yield_inc[ibin][ibin2][ibin3]->Clone((TString) ("Yield_BkgSub_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	yield_bgsub[ibin][ibin2][ibin3]->Add(background[ibin][ibin2][ibin3],-1.);



	llimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(bg_inner+0.0001);
	rlimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(bg_outer-0.0001);


	background_proj2[ibin][ibin2][ibin3] = (TH1D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->ProjectionY(Form("ProjectedBackground2%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);

	llimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-bg_outer+0.0001);
	rlimiteta = yield_inc_pTweighted[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-bg_inner-0.0001);

	background_left[ibin][ibin2][ibin3] = (TH1D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->ProjectionY(Form("LeftSideBackground%d%d%d",ibin,ibin2,ibin3),llimiteta,rlimiteta);

	      
	background_proj2[ibin][ibin2][ibin3]->Add(background_left[ibin][ibin2][ibin3]);

	dx_phi =  background_proj2[ibin][ibin2][ibin3]->GetBinWidth(1);

	background_proj2[ibin][ibin2][ibin3]->Scale(1./2/(rlimiteta-llimiteta+1));

	  
	if(ibin==3)	background_pTweighted[ibin][ibin2][ibin3] = (TH2D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->Clone((TString) ("SummedBkg_pTweighted_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+2] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	else	background_pTweighted[ibin][ibin2][ibin3] = (TH2D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->Clone((TString) ("SummedBkg_pTweighted_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));



	for(int k = 1;  k<yield_inc[ibin][ibin2][ibin3]->GetNbinsY(); k++){
	  temp1 = background_proj2[ibin][ibin2][ibin3]->GetBinContent(k);
	  err1 = background_proj2[ibin][ibin2][ibin3]->GetBinError(k);
	    
	  for(int m = 1;  m<yield_inc[ibin][ibin2][ibin3]->GetNbinsX(); m++){
	    background_pTweighted[ibin][ibin2][ibin3]->SetBinContent(m,k,temp1);
	    background_pTweighted[ibin][ibin2][ibin3]->SetBinError(m,k,err1);
	  }
	}

      

	if(ibin==3)	yield_bgsub_pTweighted[ibin][ibin2][ibin3] = (TH2D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->Clone((TString) ("Yield_BkgSub_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+2] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	else 	yield_bgsub_pTweighted[ibin][ibin2][ibin3] = (TH2D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->Clone((TString) ("Yield_BkgSub_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	yield_bgsub_pTweighted[ibin][ibin2][ibin3]->Add(background_pTweighted[ibin][ibin2][ibin3],-1.);

	
	//	if(data_mc_type_code!=0){
	  
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
	
      
	  //	}
	

	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Write();
	hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Write();
	yield_inc[ibin][ibin2][ibin3]->Write();
	yield_inc_pTweighted[ibin][ibin2][ibin3]->Write();
	yield_bgsub[ibin][ibin2][ibin3]->Write();
	yield_bgsub_pTweighted[ibin][ibin2][ibin3]->Write();
	background[ibin][ibin2][ibin3]->Write();
	background_pTweighted[ibin][ibin2][ibin3]->Write();
	hJetTrackME[ibin][ibin2][ibin3]->Write();



	if(is_pp) result_proj_canvas->cd(ibin3+1);
	else 	result_proj_canvas->cd(4*(ibin3+1)-ibin);      



	llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(2.5-0.0001);

	yield_inc_proj[ibin][ibin2][ibin3] = (TH1D*)yield_inc[ibin][ibin2][ibin3]->ProjectionY((TString) ("Yield_Proj_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),llimiteta,rlimiteta);

	yield_inc_proj[ibin][ibin2][ibin3]->Scale(1./(rlimiteta-llimiteta+1));
	
	//	yield_inc_proj[ibin][ibin2][ibin3]->Rebin(5);
	yield_inc_proj[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	yield_inc_proj[ibin][ibin2][ibin3]->SetMarkerColor(kViolet);
	yield_inc_proj[ibin][ibin2][ibin3]->SetLineColor(kViolet);
	//	background_proj[ibin][ibin2][ibin3]->Rebin(5);
	background_proj[ibin][ibin2][ibin3]->SetLineColor(kBlack);
	background_proj[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
	background_proj[ibin][ibin2][ibin3]->SetMarkerStyle(10);

	yield_inc_proj[ibin][ibin2][ibin3]->SetMinimum(background_proj[ibin][ibin2][ibin3]->GetMinimum()-0.001);
	yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(background_proj[ibin][ibin2][ibin3]->GetMaximum()+0.001);

	//	yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(.002);
	yield_inc_proj[ibin][ibin2][ibin3]->SetAxisRange(-1.5,1.5);

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
      


	if(is_pp) result_proj_canvas2->cd(ibin3+1);
	else 	result_proj_canvas2->cd(4*(ibin3+1)-ibin);      
  


	llimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(-2.5+0.0001);
	rlimiteta = yield_inc[ibin][ibin2][ibin3]->GetXaxis()->FindBin(2.5-0.0001);

	yield_inc_proj2[ibin][ibin2][ibin3] = (TH1D*)yield_inc_pTweighted[ibin][ibin2][ibin3]->ProjectionY((TString) ("Yield_Proj2_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),llimiteta,rlimiteta);

	yield_inc_proj2[ibin][ibin2][ibin3]->Scale(1./(rlimiteta-llimiteta+1));
	
	//	yield_inc_proj2[ibin][ibin2][ibin3]->Rebin(5);
	yield_inc_proj2[ibin][ibin2][ibin3]->SetMarkerStyle(10);
	yield_inc_proj2[ibin][ibin2][ibin3]->SetMarkerColor(kViolet);
	yield_inc_proj2[ibin][ibin2][ibin3]->SetLineColor(kViolet);
	//	background_proj2[ibin][ibin2][ibin3]->Rebin(5);
	background_proj2[ibin][ibin2][ibin3]->SetLineColor(kBlack);
	background_proj2[ibin][ibin2][ibin3]->SetMarkerColor(kBlack);
	background_proj2[ibin][ibin2][ibin3]->SetMarkerStyle(10);

	yield_inc_proj2[ibin][ibin2][ibin3]->SetMinimum(background_proj2[ibin][ibin2][ibin3]->GetMinimum()-0.001);
	yield_inc_proj2[ibin][ibin2][ibin3]->SetMaximum(background_proj2[ibin][ibin2][ibin3]->GetMaximum()+0.001);

	//	yield_inc_proj[ibin][ibin2][ibin3]->SetMaximum(.002);
	yield_inc_proj2[ibin][ibin2][ibin3]->SetAxisRange(-1.5,1.5);

	yield_inc_proj2[ibin][ibin2][ibin3]->Draw();

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
	background_proj2[ibin][ibin2][ibin3]->Draw("same");
      





      }//ibin3;
      
    } // ibin2

  } // ibin ( centrality ) loop

  me_proj_canvas->SaveAs((TString)("All_ME_Projections_Inclusive_"+datalabel+"_"+data_mc_type_strs[data_mc_type_code]+".png"));
  result_proj_canvas->SaveAs((TString)("All_Background_Projections_Inclusive_"+datalabel+"_"+data_mc_type_strs[data_mc_type_code]+".png"));

  result_proj_eta_canvas->SaveAs((TString)("All_dEta_Projections_Inclusive_"+datalabel+"_"+data_mc_type_strs[data_mc_type_code]+".png"));



  me_proj_canvas2->SaveAs((TString)("All_ME_Projections_Inclusive_pTweighted_"+datalabel+"_"+data_mc_type_strs[data_mc_type_code]+".png"));
  result_proj_canvas2->SaveAs((TString)("All_Background_Projections_Inclusive_pTweighted_"+datalabel+"_"+data_mc_type_strs[data_mc_type_code]+".png"));

  result_proj_eta_canvas2->SaveAs((TString)("All_dEta_Projections_Inclusive_pTweighted_"+datalabel+"_"+data_mc_type_strs[data_mc_type_code]+".png"));

  return 0;
} // main loop






  
     
