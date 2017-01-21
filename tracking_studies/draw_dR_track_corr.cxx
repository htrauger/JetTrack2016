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


Int_t draw_dR_track_corr(){

  //-------------Input files & set eta-----------------

  TString datalabel;
  float eta_ymax;

  TLatex *centtex, *pttex;

  TF1 *pol0 = new TF1("pol0","[0]+x-x",-.3,.3);


  TFile *fin_nocorr, *fin_corr;

  enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RightGen, SpilledUnderGen, UnmatchedGen, RightReco, SpilledReco, UnmatchedReco, RecoGenSube0,RecoGenNoSube0,GenGenSube0,GenGenNoSube0,MatchedRecoGenSube0,MatchedRecoGenNoSube0,SwappedRecoGenSube0,SwappedRecoGenNoSube0, UnMatchedRecoGenSube0,UnMatchedRecoGenNoSube0,n_data_mc_types};


  TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RightGenJet_GenTrack","SpilledUnderJet_GenTrack","UnmatchedGenJet_GenTrack","RightRecoJet_GenTrack","SpilledReco_GenTrack","UnmatchedReco_GenTrack","RecoJet_GenTrack_Sube0","RecoJet_GenTrack_NoSube0","GenJet_GenTrack_Sube0","GenJet_GenTrack_NoSube0","MatchedRecoJet_GenTrack_Sube0","MatchedRecoJet_GenTrack_NoSube0","SwappedRecoJet_GenTrack_Sube0","SwappedRecoJet_GenTrack_NoSube0","UnmatchedRecoJet_GenTrack_Sube0","UnmatchedRecoJet_GenTrack_NoSube0",};

 int data_mc_type_code = 0;

 bool is_pp = kFALSE;

  TString desc = data_mc_type_strs[data_mc_type_code];

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins;
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 9;

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins+1] = {05, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.5<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
  
  TString ref_name;
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
  
  TH2D* hJetTrackSignalBackground_dRcorr[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackground_no_dRcorr[nCBins][nPtBins][nTrkPtBins];
  TH2D* ratio[nCBins][nPtBins][nTrkPtBins];
  
  float norm_temp, width_temp, width_temp_x, width_temp_y, width_temp_ref, max_bin, max_cont,bc,err,err1, dx_phi, dx_eta, temp1;
    

  TH1F* all_jets_corrpT_dRcorr[nCBins][nPtBins];
  TH1F* all_jets_corrpT_no_dRcorr[nCBins][nPtBins];


  //-----------------------
  // Start getting histos
  //-----------------------

 TCanvas *c = new TCanvas("canvas","",0,0,2000,2400);
 c->Divide(4,9,0.00001,0.00001);


 for(int g = 0; g<1; g++){

   if(g==1) is_pp = kTRUE;


   if(data_mc_type_code==0){
     if(is_pp){
       datalabel = "pp";
       fin_nocorr = new TFile("../data_raw_correlations/Data_pp_No_dR_Aug15.root","READ");
       fin_corr = new TFile("../data_raw_correlations/Data_pp_July10.root","READ");
     }else{
       datalabel = "PbPb";
       fin_nocorr = new TFile("../data_raw_correlations/Data_PbPb_No_dR_Aug15.root","READ");
       fin_corr = new TFile("../data_raw_correlations/Data_PbPb_Partial_Aug14.root","READ");
     }
   }else{
     cout<<"not set up for MC yet"<<endl;
   }


  
   for (int ibin=0;ibin<nCBins;ibin++){

     if(is_pp && ibin > 0 ) continue;
  
     for (int ibin2=0;ibin2<nPtBins;ibin2++){ 

       cout<<desc + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]<<endl;
      
       all_jets_corrpT_dRcorr[ibin][ibin2] = (TH1F*)fin_corr->Get((TString) (data_mc_type_strs[data_mc_type_code] + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("dR_corr_all_jets_corrpT_"+ datalabel+CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));

       all_jets_corrpT_no_dRcorr[ibin][ibin2] = (TH1F*)fin_nocorr->Get((TString) (data_mc_type_strs[data_mc_type_code] + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("No_dR_corr_all_jets_corrpT_"+datalabel+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));

       for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

	 hJetTrackSignalBackground_dRcorr[ibin][ibin2][ibin3] = (TH2D*) fin_corr->Get((TString)(desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_dRcorr"+datalabel+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	 hJetTrackSignalBackground_no_dRcorr[ibin][ibin2][ibin3] = (TH2D*) fin_nocorr->Get((TString)(desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_No_dRcorr"+datalabel+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
      	  
	 norm_temp = all_jets_corrpT_dRcorr[ibin][ibin2]->Integral();

	 hJetTrackSignalBackground_dRcorr[ibin][ibin2][ibin3]->Scale(1./norm_temp);

	 if(ibin3>4)	 hJetTrackSignalBackground_dRcorr[ibin][ibin2][ibin3]->Rebin2D(10,8);
	 else 	 hJetTrackSignalBackground_dRcorr[ibin][ibin2][ibin3]->Rebin2D(10,4);

	
	 norm_temp = all_jets_corrpT_no_dRcorr[ibin][ibin2]->Integral();

	 hJetTrackSignalBackground_no_dRcorr[ibin][ibin2][ibin3]->Scale(1./norm_temp);

	 if(ibin3>4) hJetTrackSignalBackground_no_dRcorr[ibin][ibin2][ibin3]->Rebin2D(10,8);
	 else	 hJetTrackSignalBackground_no_dRcorr[ibin][ibin2][ibin3]->Rebin2D(10,4);

	 ratio[ibin][ibin2][ibin3] = (TH2D*)	hJetTrackSignalBackground_dRcorr[ibin][ibin2][ibin3]->Clone((TString)("Ratio_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3]));

	 ratio[ibin][ibin2][ibin3]->Divide(hJetTrackSignalBackground_no_dRcorr[ibin][ibin2][ibin3]);

	 //	 c->cd(5*(ibin3+1)-ibin-4*g);
	 c->cd(4*(ibin3+1)-ibin);



	 ratio[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-1.,1.);
	 ratio[ibin][ibin2][ibin3]->GetYaxis()->SetRangeUser(-1.,1.);
	 ratio[ibin][ibin2][ibin3]->SetMinimum(0.9);
	 ratio[ibin][ibin2][ibin3]->SetMaximum(1.1);

	 ratio[ibin][ibin2][ibin3]->GetXaxis()->SetLabelSize(0.08);
	 ratio[ibin][ibin2][ibin3]->GetYaxis()->SetLabelSize(0.08);
	 ratio[ibin][ibin2][ibin3]->GetZaxis()->SetLabelSize(0.08);


	 ratio[ibin][ibin2][ibin3]->GetXaxis()->SetTitleSize(0.08);
	 ratio[ibin][ibin2][ibin3]->GetYaxis()->SetTitleSize(0.08);
	 ratio[ibin][ibin2][ibin3]->GetZaxis()->SetTitleSize(0.08);

	 ratio[ibin][ibin2][ibin3]->GetXaxis()->SetNdivisions(505);
	 ratio[ibin][ibin2][ibin3]->GetYaxis()->SetNdivisions(505);
	 ratio[ibin][ibin2][ibin3]->GetZaxis()->SetNdivisions(505);
	 

	 ratio[ibin][ibin2][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
	 ratio[ibin][ibin2][ibin3]->GetYaxis()->SetTitle("#Delta#phi");
	 ratio[ibin][ibin2][ibin3]->GetZaxis()->SetTitle("dR_corr/No_dR_corr");

	 


	 ratio[ibin][ibin2][ibin3]->Draw("surf1");


	 gPad->SetTheta(55);

	  
	 TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	 centtex->SetNDC();
	 centtex->Draw();
	 
	 TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	 pttex->SetNDC();
	 pttex->Draw();
      
    
       }//ibin3;
      
     } // ibin2

   } // ibin ( centrality ) loop
 }
   c->SaveAs("Ratio_dRcorr_to_no_dRcorr.png");
   
  return 0;
} // main loop






  
     
