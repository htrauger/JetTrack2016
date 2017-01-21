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
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLatex.h"
#include "THStack.h"

#include <iostream>
#include <vector>
#include <fstream>

#include "../../JetTrack2016_functions.h"

//******************************************************
// VERSION TO DO JET-SHAPES FOR PAS
//******************************************************

Int_t jet_shapes_result(bool use_highpT_bin = kTRUE){

  TFile *fMC[6], *f_ref[7];

  TCanvas *c_jetshape[5];
   
  TString jetetacut, etalabel,centlabel,pTlabel,ajlabel;
  float eta_ymax;

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins, limR;
  float deta, dphi, r, bc, bg_err,temp1,temp2, rbin, temperr, err, width_temp_x, width_temp_y, width_temp, norm_temp, zerobin, temp, norm_tot, err_temp, cont;
  
 
  float RBins[20] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.,1.2,1.4,1.6,1.8,2.0};
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 7;

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  

  float CBins[nCBins+1] = {0,20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins+1] = {05, 1, 2, 3, 4, 8, 300};
  TString TrkPtBin_strs2[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt300" };
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.5<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","pT>8"};
  
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.25);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
 
 
  TH2D* result[6][nCBins][nPtBins][nTrkPtBins];

  TH2D* resultMC_gen[6][nCBins][nPtBins][nTrkPtBins];
  TH2D* resultMC_reco[6][nCBins][nPtBins][nTrkPtBins];

  TH2D* resultMC2[6][nCBins][nPtBins][nTrkPtBins];

  TH2D* resultMC[6][nCBins][nPtBins][nTrkPtBins];

  TF1 *gaus_eta[4][1][4];
  TF1 *gaus_phi[4][1][4];
 
  float par0, par1, par2, par3, par4;

  TH1D *background_syst_rebin[12][6][4][5];

   
  TH1D* JetShapeMC[6][nCBins][nTrkPtBins];

  TH1D* JetShape[6][nCBins][nTrkPtBins];
  TH1D* JetShape2[6][nCBins][nTrkPtBins];
  TH1D* JetShape_noerr_up[6][nCBins][nTrkPtBins];
  TH1D* JetShape_noerr_down[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ref[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ref_diff[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff[6][nCBins][nTrkPtBins];

  TH1D* JetShape_ref_ratio[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ratio[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ratio2[6][nCBins][nTrkPtBins];

  TH1D* JetShape_diff2[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff_noerr_up[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff_noerr_down[6][nCBins][nTrkPtBins];
  TH1D* JetShape2_geo[6][nCBins][nTrkPtBins];

  TH1D* JetShape_syst[6][nCBins][nTrkPtBins];
  TGraphAsymmErrors* JetShape_graph[6][nCBins][nTrkPtBins];



  THStack *JetShape_Stack_Up[6][nCBins];
  THStack *JetShape_Diff_Stack_Up[6][nCBins];

  THStack *JetShape_Stack_Down[6][nCBins];
  THStack *JetShape_Diff_Stack_Down[6][nCBins];

 
  double temp_cont, nextr, nextl, cos_weight,me00_range, mc_error;
  
  TString stem, datalabel,me00_range_string,stem_mc;

 
 
  float norm;


  TFile *fin = new TFile("../../me_correct/pp_Inclusive_Correlations.root", "READ");
   TFile *fin_run1 = new TFile("../../../JetTrack2015/bg_fit/Dijet_Correlations.root", "READ");
   //TFile *fin_run1 = new TFile("../../../JetTrack2015/bg_fit/Dijet_Correlations_NoPtWeight.root", "READ");
  TFile *fout  = new TFile("Jet_Shapes.root","RECREATE"); 
  /*
  TFile *f_spillover = new TFile("../spill_over/Dijet_SpillOvers.root");

  TFile *f_jff_gen2 = new TFile("../me_correct_mc/Pythia_GenJet_GenTrack_Dijet_Correlations.root");  
  TFile *f_jff_reco2 = new TFile("../me_correct_mc/Pythia_RecoJet_RecoTrack_Dijet_Correlations.root"); 
 
  TFile *f_jff_reco_reco = new TFile("../me_correct_mc/HydJet_RecoJet_RecoTrack_Dijet_Correlations.root"); 
  TFile *f_jff_reco = new TFile("../me_correct_mc/HydJet_RecoJet_GenTrack_Sube0_Dijet_Correlations.root"); 
  TFile *f_jff_gen = new TFile("../me_correct_mc/HydJet_GenJet_GenTrack_Sube0_Dijet_Correlations.root");  
  */


  // TFile *f_jff_reco = new TFile("../me_correct_mc/Pythia_RecoJet_GenTrack_Dijet_Correlations.root"); 
 

  

  TF2* ClosureFit = new TF2("ClosureFit", "[0]/2/TMath::Pi()/[1]/[2]*TMath::Exp(-1.*(x*x/[1]/[1]/2))*TMath::Exp(-1.*(y*y/[2]/[2]/2))",-5.,5.,-TMath::Pi()/2.,3*TMath::Pi()/2.);


  //-----------------------
  // Start getting histos
  //-----------------------
  for(int g=0; g<2; g++){

    cout<<"here!"<<endl;

    switch(g){
    case 0:
      stem = "Yield_BkgSub_pp_";
      datalabel = "Inclusive";
      break;
    case 1:
      stem = "Yield_BkgSub_pTweighted";
      //stem = "Yield_BkgSub_";
      datalabel = "Inclusive";
      break;
    }

    for (int ibin=0;ibin<nCBins;ibin++){

      if(ibin > 0) continue; //temporarily I only consider 
  
      for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
    
	for (int ibin3=0;ibin3<nTrkPtBins-1;ibin3++){

	  if(g==0){
	    result[g][ibin][ibin2][ibin3] = (TH2D*)fin_run1->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + datalabel+"_Run1_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  }else{
	    if(ibin3<5){
   cout<<stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]<<endl;
	    
	    result[g][ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    }else if( ibin3==5){
	      result[g][ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_TrkPt8_TrkPt12"))->Clone((TString) (stem + datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    }
	    cout<<"got hist"<<endl;

	  }

	  width_temp_x = result[g][ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1);
	  width_temp_y = result[g][ibin][ibin2][ibin3]->GetYaxis()->GetBinWidth(1);
	  
	  //	  result[g][ibin][ibin2][ibin3]->Scale(1./width_temp_x/width_temp_y);


	  /*
	  if(g==0)	  resultMC2[g][ibin][0][ibin3] = new TH2D((TString)("Closure_2D_PbPb_"+CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]),"",500,-5,5,200,-TMath::Pi()/2,3*TMath::Pi()/2);
	  else	  resultMC2[g][ibin][0][ibin3] = new TH2D((TString)("Closure_2D_pp_"+CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]),"",500,-5,5,200,-TMath::Pi()/2,3*TMath::Pi()/2);
	

	  if(g==0){
	   
	  gaus_eta[ibin][0][ibin3] = (TF1*)f_spillover->Get((TString)("GausFit_Eta_pTweighted_Leading_Hydjet_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));



	   gaus_phi[ibin][0][ibin3] = (TF1*)f_spillover->Get((TString)("GausFit_Phi_pTweighted_Leading_Hydjet_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));


	    par0= gaus_eta[ibin][0][ibin3]->GetParameter(1);
	    par1= gaus_eta[ibin][0][ibin3]->GetParameter(2);
	    par2= gaus_phi[ibin][0][ibin3]->GetParameter(2);

	    ClosureFit->SetParameter(0,par0);
	    ClosureFit->SetParameter(1,par1);
	    ClosureFit->SetParameter(2,par2);
	    

	    resultMC2[g][ibin][0][ibin3]->Eval(ClosureFit);

	    resultMC2[g][ibin][ibin2][ibin3]->Scale(width_temp_x*width_temp_y);
	        
	
	    result[g][ibin][ibin2][ibin3]->Add(resultMC2[g][ibin][0][ibin3],-1.);
	 
	  }

	  if(ibin==0&&g==1){

	   

	      resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_gen2->Get((TString)("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_Leading_PythiaGen"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_reco2->Get((TString)("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_Leading_PythiaReco"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	
	    resultMC_reco[g][ibin][0][ibin3]->Add(resultMC_gen[g][ibin][0][ibin3],-1.);
	  
	    resultMC_reco[g][ibin][0][ibin3]->Scale(1./width_temp_x/width_temp_y);

	  }else if(g==0){

	   

	      resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_gen->Get((TString)("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_Leading_PythiaForPbPbGen"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      if(ibin3>3)      resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_reco_reco->Get((TString)("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_Leading_PythiaReco"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	      else      resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_reco->Get((TString)("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_Leading_PythiaReco"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    resultMC_reco[g][ibin][0][ibin3]->Add(resultMC_gen[g][ibin][0][ibin3],-1.);
	  
	    resultMC_reco[g][ibin][0][ibin3]->Scale(1./width_temp_x/width_temp_y);

	  }

	  if(g%2==0||ibin==0){
	    
	    for(int k = 1; k< resultMC_reco[g][ibin][0][ibin3]->GetNbinsX(); k++){
	      for(int m = 1; m< resultMC_reco[g][ibin][0][ibin3]->GetNbinsY(); m++){

		deta = resultMC_reco[g][ibin][0][ibin3]->GetXaxis()->GetBinCenter(k);
		dphi = resultMC_reco[g][ibin][0][ibin3]->GetYaxis()->GetBinCenter(m);
		
		r = TMath::Sqrt(deta*deta+dphi*dphi);
	  
		if(r>0.3){
		  resultMC_reco[g][ibin][0][ibin3]->SetBinContent(k,m,0.);
		  resultMC_reco[g][ibin][0][ibin3]->SetBinError(k,m,0.);

	      }
	    }
	    }
	  }


	  if(g==0) result[g][ibin][ibin2][ibin3]->Add(resultMC_reco[g][ibin][0][ibin3],-1.);
	  else 	  result[g][ibin][ibin2][ibin3]->Add(resultMC_reco[g][0][0][ibin3],-1.);


	  resultMC[g][ibin][ibin2][ibin3] = (TH2D*)resultMC2[g][ibin][ibin2][ibin3]->Clone((TString)("Combined_MC_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  resultMC[g][ibin][ibin2][ibin3]->Add(resultMC_reco[g][0][0][ibin3],0.02/0.07);
	 	 

	  if(ibin3==0){

	    resultMC[g][ibin][ibin2][6] = (TH2D*) result[g][ibin][ibin2][ibin3]->Clone(((TString) ("Summed_result_MC"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])));
	  }else if(ibin3>0&&(use_highpT_bin||ibin3<5)){

	    resultMC[g][ibin][ibin2][6]->Add(resultMC[g][ibin][ibin2][ibin3]);
	  }
	}
	  */

	
	  if(ibin3==0){
	    result[g][ibin][ibin2][6] = (TH2D*) result[g][ibin][ibin2][ibin3]->Clone(((TString) ("Summed_result_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])));

	  }else if(ibin3>0&&(use_highpT_bin||ibin3<5)){

	    result[g][ibin][ibin2][6]->Add(result[g][ibin][ibin2][ibin3]);
	 
	  }
	}
    
	for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

	  
	  JetShape2[g][ibin][ibin3] = new TH1D((TString)("JetShape2_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  

	  JetShape2_geo[g][ibin][ibin3] = new TH1D((TString)("JetShape_Geometry_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",19,RBins);  

	  

	  int kmin = result[g][ibin][0][ibin3]->GetXaxis()->FindBin(-TMath::Pi()/2.+.0001);
	  int kmax = result[g][ibin][0][ibin3]->GetXaxis()->FindBin(TMath::Pi()/2.-.0001);
      

	  int lmin = result[g][ibin][0][ibin3]->GetYaxis()->FindBin(-TMath::Pi()/2.+.0001);
	  int lmax = result[g][ibin][0][ibin3]->GetYaxis()->FindBin(TMath::Pi()/2.-.0001);
	
	
	  for(int k = kmin; k<kmax+1; k++){  //
	    for(int l = lmin; l<lmax+1; l++){

	      deta = result[g][ibin][0][ibin3]->GetXaxis()->GetBinCenter(k);
	      dphi = result[g][ibin][0][ibin3]->GetYaxis()->GetBinCenter(l);
	   
	      r = TMath::Sqrt(deta*deta+dphi*dphi);
	      
	      
	      rbin = JetShape2[g][ibin][ibin3]->FindBin(r);
	   
	      //catch-all, if you want
	      /*
	      if(rbin > 14){
		rbin = 14;
		cout<<r<<" "<<rbin<<endl;
	      }
	      */
	      bc = result[g][ibin][0][ibin3]->GetBinContent(k,l);
	      err = result[g][ibin][0][ibin3]->GetBinError(k,l);

	      temp1 = JetShape2[g][ibin][ibin3]->GetBinContent(rbin);
	      temp2 = JetShape2_geo[g][ibin][ibin3]->GetBinContent(rbin);
	    
	      temperr = JetShape2[g][ibin][ibin3]->GetBinError(rbin);
	  
	      JetShape2[g][ibin][ibin3]->SetBinContent(rbin,temp1+bc);
	      JetShape2[g][ibin][ibin3]->SetBinError(rbin,TMath::Sqrt(temperr*temperr+err*err));
	      
	      JetShape2_geo[g][ibin][ibin3]->SetBinContent(rbin,temp2+1.);
	      JetShape2_geo[g][ibin][ibin3]->SetBinError(rbin,0.);
	      
	      
	  	  
	    }
	  }

	  /*
	  JetShapeMC[g][ibin][ibin3] = new TH1D((TString)("JetShape_MC_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",16,RBins);  

	  for(int k = kmin; k<kmax+1; k++){  //
	    for(int l = lmin; l<lmax+1; l++){

	      deta = resultMC[g][ibin][0][ibin3]->GetXaxis()->GetBinCenter(k);
	      dphi = resultMC[g][ibin][0][ibin3]->GetYaxis()->GetBinCenter(l);
	   
	      r = TMath::Sqrt(deta*deta+dphi*dphi);
	  
	      rbin = JetShapeMC[g][ibin][ibin3]->FindBin(r);

	      bc = resultMC[g][ibin][0][ibin3]->GetBinContent(k,l);
	      err = resultMC[g][ibin][0][ibin3]->GetBinError(k,l);

	      temp1 = JetShapeMC[g][ibin][ibin3]->GetBinContent(rbin);
	      temperr = JetShapeMC[g][ibin][ibin3]->GetBinError(rbin);
	 	  
	      JetShapeMC[g][ibin][ibin3]->SetBinContent(rbin,temp1+bc);
	      JetShapeMC[g][ibin][ibin3]->SetBinError(rbin,TMath::Sqrt(temperr*temperr+err*err));
	    }

	  
	}
	*/

	  
	  //  JetShape2[g][ibin][ibin3]->Divide(JetShape2_geo[g][ibin][ibin3]);
	
	  JetShape[g][ibin][ibin3] = (TH1D*)JetShape2[g][ibin][ibin3]->Clone((TString)("JetShape_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	
	  /*

	    if(g==0&&ibin==0){
	      if(ibin3>0&&ibin3<5){
		f_ref[ibin3] = new TFile((TString)("JetShapesRef_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]+".root"),"READ");
	      }else if(ibin3==6){
		f_ref[ibin3] = new TFile((TString)("JetShapesRef_"+TrkPtBin_strs[1]+"_"+TrkPtBin_strs[5]+".root"),"READ");
	      }
	    }
	    if(ibin3==1){

	      JetShape_ref[0][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pbpb_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone("JetShape_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]);
	
	      JetShape_ref[1][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pp_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone("JetShape_pp_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]);
     
	    }else if(ibin3>0&&ibin3<5){

	      JetShape_ref[0][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pbpb_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone("JetShape_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]);
	
	      JetShape_ref[1][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pp_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone("JetShape_pp_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]);
     

	    }else if(ibin3==6){

	      JetShape_ref[0][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pbpb_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone("JetShape_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[1]+"_"+TrkPtBin_strs[5]);
	      JetShape_ref[1][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pp_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone("JetShape_pp_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[1]+"_"+TrkPtBin_strs[5]);
     

	    }

	    if(ibin3==6){
	   
	      if(ibin==0){
		JetShape_ref[1][ibin][ibin3]->SetBinContent(1,12.485);
		JetShape_ref[1][ibin][ibin3]->SetBinContent(2,4.683);
		JetShape_ref[1][ibin][ibin3]->SetBinContent(3,1.5);
		JetShape_ref[1][ibin][ibin3]->SetBinContent(4,0.722);
		JetShape_ref[1][ibin][ibin3]->SetBinContent(5,0.396);
		JetShape_ref[1][ibin][ibin3]->SetBinContent(6,0.215);

		JetShape_ref[1][ibin][ibin3]->SetBinError(1,0.024);
		JetShape_ref[1][ibin][ibin3]->SetBinError(2,0.014);
		JetShape_ref[1][ibin][ibin3]->SetBinError(3,0.006);
		JetShape_ref[1][ibin][ibin3]->SetBinError(4,0.004);
		JetShape_ref[1][ibin][ibin3]->SetBinError(5,0.002);
		JetShape_ref[1][ibin][ibin3]->SetBinError(6,0.002);



		JetShape_ref[0][ibin][ibin3]->SetBinContent(1,12.636);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(2,4.560);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(3,1.318);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(4,0.707);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(5,0.486);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(6,0.294);

		JetShape_ref[0][ibin][ibin3]->SetBinError(1,0.05);
		JetShape_ref[0][ibin][ibin3]->SetBinError(2,0.028);
		JetShape_ref[0][ibin][ibin3]->SetBinError(3,0.012);
		JetShape_ref[0][ibin][ibin3]->SetBinError(4,0.009);
		JetShape_ref[0][ibin][ibin3]->SetBinError(5,0.009);
		JetShape_ref[0][ibin][ibin3]->SetBinError(6,0.009);

	      }


	      if(ibin==1){
		JetShape_ref[0][ibin][ibin3]->SetBinContent(1,12.93);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(2,4.36);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(3,1.28);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(4,0.685);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(5,0.455);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(6,0.273);

		JetShape_ref[0][ibin][ibin3]->SetBinError(1,0.046);
		JetShape_ref[0][ibin][ibin3]->SetBinError(2,0.025);
		JetShape_ref[0][ibin][ibin3]->SetBinError(3,0.011);
		JetShape_ref[0][ibin][ibin3]->SetBinError(4,0.007);
		JetShape_ref[0][ibin][ibin3]->SetBinError(5,0.007);
		JetShape_ref[0][ibin][ibin3]->SetBinError(6,0.007);

	      }
   
	      if(ibin==2){
		JetShape_ref[0][ibin][ibin3]->SetBinContent(1,13.168);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(2,4.281);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(3,1.278);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(4,0.640);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(5,0.399);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(6,0.233);

		JetShape_ref[0][ibin][ibin3]->SetBinError(1,0.071);
		JetShape_ref[0][ibin][ibin3]->SetBinError(2,0.038);
		JetShape_ref[0][ibin][ibin3]->SetBinError(3,0.016);
		JetShape_ref[0][ibin][ibin3]->SetBinError(4,0.010);
		JetShape_ref[0][ibin][ibin3]->SetBinError(5,0.008);
		JetShape_ref[0][ibin][ibin3]->SetBinError(6,0.007);

	      }


	      if(ibin==3){
		JetShape_ref[0][ibin][ibin3]->SetBinContent(1,(13.086+12.98)/2.);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(2,(4.356+4.322)/2.);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(3,(1.320+1.396)/2.);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(4,(0.637+0.707)/2.);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(5,(0.395+0.371)/2.);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(6,(0.230+0.215)/2.);

		JetShape_ref[0][ibin][ibin3]->SetBinError(1,0.13);
		JetShape_ref[0][ibin][ibin3]->SetBinError(2,0.071);
		JetShape_ref[0][ibin][ibin3]->SetBinError(3,0.030);
		JetShape_ref[0][ibin][ibin3]->SetBinError(4,0.019);
		JetShape_ref[0][ibin][ibin3]->SetBinError(5,0.012);
		JetShape_ref[0][ibin][ibin3]->SetBinError(6,0.010);

	      }

	      
	      if(ibin==1){

		JetShape_ref[0][ibin][ibin3]->Scale(0.53);
		JetShape_ref[0][ibin][ibin3]->Add(JetShape_ref[1][ibin-1][ibin3],0.47);
	

	      }

	    }

	    
	  

	    if(((ibin3>0&&ibin3<5)||ibin3==6)){

		JetShape_ref[0][ibin][ibin3]->SetMarkerStyle(20);
		JetShape_ref[0][ibin][ibin3]->SetMarkerSize(1);
		JetShape_ref[0][ibin][ibin3]->SetMarkerColor(kCyan);
		JetShape_ref[0][ibin][ibin3]->SetLineColor(kCyan);


		JetShape_ref[1][ibin][ibin3]->SetMarkerStyle(20);
		JetShape_ref[1][ibin][ibin3]->SetMarkerSize(1);
		JetShape_ref[1][ibin][ibin3]->SetMarkerColor(kCyan);
		JetShape_ref[1][ibin][ibin3]->SetLineColor(kCyan);
	
	    }

	  */
	
	  for(int k = 1; k<JetShape2[g][ibin][ibin3]->GetNbinsX()+1; k++){
	    
	    float temp = JetShape2[g][ibin][ibin3]->GetBinContent(k)/JetShape2[g][ibin][ibin3]->GetBinWidth(k);
	    float err = JetShape2[g][ibin][ibin3]->GetBinError(k)/JetShape2[g][ibin][ibin3]->GetBinWidth(k);
	    
	    JetShape2[g][ibin][ibin3]->SetBinContent(k,temp);
	    JetShape2[g][ibin][ibin3]->SetBinError(k,err);

	    JetShape[g][ibin][ibin3]->SetBinContent(k,temp);
	    JetShape[g][ibin][ibin3]->SetBinError(k,err);

	  }	    
	  cout<<"here"<<endl;

	  //	  norm =  JetShape2[g][ibin][ibin3]->Integral(1,6,"width");

	  norm = 1.;
      
	  if(ibin3==6){

	    cout<<"norm is "<<norm<<endl;

	    //	    if(!use_highpT_bin) norm = JetShape[g][ibin][ibin3]->GetBinContent(1)/JetShape_ref[g][ibin][ibin3]->GetBinContent(1);	 

	    JetShape[g][ibin][ibin3]->Scale(1./norm);  
	    JetShape[g][ibin][0]->Scale(1./norm);  
	    JetShape[g][ibin][1]->Scale(1./norm);  
	    JetShape[g][ibin][2]->Scale(1./norm);  
	    JetShape[g][ibin][3]->Scale(1./norm);  
	    JetShape[g][ibin][4]->Scale(1./norm);  
	    JetShape[g][ibin][5]->Scale(1./norm);  

	 
	  }
	  
	  
	  if(ibin<3){
	    JetShape2[g][ibin][ibin3]->GetYaxis()->SetLabelSize(0.);
	    JetShape2[g][ibin][ibin3]->GetYaxis()->SetTitleSize(0.);

	  }
	  /*
	
	  JetShapeMC[g][ibin][ibin3]->Divide(JetShape2_geo[g][ibin][ibin3]);

	  JetShapeMC[g][ibin][ibin3]->Scale(1./norm);
	  */
	
	  fout->cd();
	  
	  //	  JetShape2[g][ibin][ibin3]->SetAxisRange(0.,0.99,"x");
	  //JetShape[g][ibin][ibin3]->SetAxisRange(0.,0.99,"x");
	  JetShape2[g][ibin][ibin3]->Write();


	  cout<<"here "<<ibin<<" "<<ibin3<<endl;

	  if(ibin3<6){
	    background_syst_rebin[g][ibin][ibin3][0] = (TH1D*)fin_run1->Get((TString)("Syst_pp_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Syst_" + datalabel+"_Run1_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    
	  }else{
	  
	    background_syst_rebin[g][ibin][ibin3][0] = (TH1D*) background_syst_rebin[g][ibin][0][0] ->Clone((TString) ("Syst_" + datalabel+"_Run1_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[0]+"_" +TrkPtBin_strs[6]));
	    for(int k = 0; k<ibin3; k++){
	      background_syst_rebin[g][ibin][ibin3][0]->Add(background_syst_rebin[g][ibin][k][0]);
	    }
	  }

	  cout<<"and here"<<endl;

	  JetShape_syst[g][ibin][ibin3] = (TH1D*) JetShape2[g][ibin][ibin3]->Clone((TString)("Jet_Shape_SystErr_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	

	  for(int k = 1; k<JetShape_syst[g][ibin][ibin3]->GetNbinsX()+1; k++){

	    bg_err = background_syst_rebin[g][ibin][ibin3][0]->GetBinError(1)/background_syst_rebin[g][ibin][ibin3][0]->GetBinWidth(1)/5.*JetShape2_geo[g][ibin][ibin3]->GetBinContent(k)/JetShape2_geo[g][ibin][ibin3]->GetBinWidth(k)*width_temp_x*width_temp_y;

	    //	    bg_err = 0.;

	    bc = TMath::Abs(JetShape_syst[g][ibin][ibin3]->GetBinContent(k))+TMath::Abs(JetShape_syst[g][ibin][ibin3]->GetBinError(k));

	    
	    
	    JetShape_syst[g][ibin][ibin3]->SetBinError(k,TMath::Sqrt(bc*bc*(.05*.05+0.04*0.04+0.03*0.03+0.07*0.07)+bg_err*bg_err)); 
	 	      
	  }
	
	  JetShape_graph[g][ibin][ibin3] = new TGraphAsymmErrors(JetShape_syst[g][ibin][ibin3]);

	  JetShape_graph[g][ibin][ibin3]->SetName((TString)("Jet_Shape_SystErrGraph_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	  JetShape_graph[g][ibin][ibin3] = new TGraphAsymmErrors(JetShape_syst[g][ibin][ibin3]);

	  JetShape_graph[g][ibin][ibin3]->SetName((TString)("Jet_Shape_SystErrGraph_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	  cout<<"and even here"<<endl;

	  for(int k = 1; k<JetShape_syst[g][ibin][ibin3]->GetNbinsX()+1; k++){


	    float evalpt = JetShape_syst[g][ibin][ibin3]->GetBinCenter(k);
 
	    //    mc_error = JetShapeMC[g][ibin][ibin3]->GetBinContent(k);
	    mc_error = 0.;

	    JetShape_graph[g][ibin][ibin3]->SetPointEYhigh(k-1,TMath::Sqrt(JetShape_syst[g][ibin][ibin3]->GetBinError(k)*JetShape_syst[g][ibin][ibin3]->GetBinError(k)+mc_error*mc_error/4.));

	    JetShape_graph[g][ibin][ibin3]->SetPointEYlow(k-1,TMath::Sqrt(JetShape_syst[g][ibin][ibin3]->GetBinError(k)*JetShape_syst[g][ibin][ibin3]->GetBinError(k)+mc_error*mc_error/4.));


	  }
	  JetShape_graph[g][ibin][ibin3]->SetMarkerStyle(20);
	  JetShape_graph[g][ibin][ibin3]->SetMarkerSize(1);
	  JetShape_graph[g][ibin][ibin3]->SetMarkerColor(kWhite);
	  JetShape_graph[g][ibin][ibin3]->SetFillColor(kBlack);
	  JetShape_graph[g][ibin][ibin3]->SetFillStyle(3004);
	  
	  if(g==1){

	    JetShape_syst[g+1][ibin][ibin3] = (TH1D*) JetShape_syst[g-1][ibin][ibin3]->Clone((TString)("Jet_Shape_Syst_Diff_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    JetShape_syst[g+1][ibin][ibin3]->Divide(JetShape_syst[g][ibin][ibin3]);


	    JetShape_graph[g+1][ibin][ibin3] = new TGraphAsymmErrors(JetShape_syst[g+1][ibin][ibin3]);

	    JetShape_graph[g+1][ibin][ibin3]->SetName((TString)("Jet_Shape_Diff_Graph_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	    JetShape_graph[g+1][ibin][ibin3] = new TGraphAsymmErrors(JetShape_syst[g+1][ibin][ibin3]);

	    JetShape_graph[g+1][ibin][ibin3]->SetName((TString)("Jet_Shape_Diff_Graph_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  	  
	    for(int k = 1; k<JetShape_syst[g+1][ibin][ibin3]->GetNbinsX()+1; k++){

	      //   mc_error = JetShapeMC[g][ibin][ibin3]->GetBinContent(k);
	      mc_error = 0.;

	      float evalpt = JetShape_syst[g+1][ibin][ibin3]->GetBinCenter(k);
	
	      JetShape_graph[g+1][ibin][ibin3]->SetPointEYhigh(k-1,TMath::Sqrt(JetShape_syst[g+1][ibin][ibin3]->GetBinError(k)*JetShape_syst[g+1][ibin][ibin3]->GetBinError(k)+mc_error*mc_error/4.));

	      JetShape_graph[g+1][ibin][ibin3]->SetPointEYlow(k-1,TMath::Sqrt(JetShape_syst[g+1][ibin][ibin3]->GetBinError(k)*JetShape_syst[g+1][ibin][ibin3]->GetBinError(k)+mc_error*mc_error/4.));


	    }
	 
	    JetShape_graph[g+1][ibin][ibin3]->SetMarkerStyle(20);
	    JetShape_graph[g+1][ibin][ibin3]->SetMarkerSize(1);
	    JetShape_graph[g+1][ibin][ibin3]->SetMarkerColor(kWhite);
	    JetShape_graph[g+1][ibin][ibin3]->SetFillColor(kBlack);
	    JetShape_graph[g+1][ibin][ibin3]->SetFillStyle(3004);
	  }

	}
	  
      }

    }

    cout<<"done "<<g<<endl;
  }
 

  cout<<"ready to draw"<<endl;

  for(int k = 0; k<7; k++){
   
    TString canvas_name = "JetShape_Comparison_";
  
    if(k<6){ 
      canvas_name+=TrkPtBin_strs[k]; canvas_name+="_"; canvas_name+=TrkPtBin_strs[k+1];  
    } else {
      canvas_name+=TrkPtBin_strs[1]; canvas_name+="_"; canvas_name+=TrkPtBin_strs[6];  
    }
   
    c_jetshape[k] = new TCanvas(canvas_name," ",10,10,1500,800);
    c_jetshape[k]->Divide(4,2,0.,0.);


    for(int ibin = 0; ibin<4; ibin++){
  
      if(ibin>0) continue;
       
      c_jetshape[k]->cd(4-ibin);

      JetShape2[0][ibin][k]->GetXaxis()->SetTitle("Radius(r)");
      JetShape2[0][ibin][k]->GetXaxis()->CenterTitle(true);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelFont(42);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelOffset(0.02);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelSize(0.08);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleSize(0.08);
      JetShape2[0][ibin][k]->GetXaxis()->SetTickLength(0.025);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleOffset(1.3);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleFont(42);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitle("#rho(#Deltar)");
      JetShape2[0][ibin][k]->GetYaxis()->CenterTitle(true);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelFont(42);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelOffset(0.004);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.075);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitleSize(0.09);
      JetShape2[0][ibin][k]->GetYaxis()->SetTickLength(0.025);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitleOffset(1.);
     
      // JetShape2[0][ibin][k]->SetMinimum(0.008);
      //JetShape2[0][ibin][k]->SetMaximum(14.5);

      JetShape2[0][ibin][k]->SetMinimum(0.01);
      JetShape2[0][ibin][k]->SetMaximum(1500.1);

      JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.07);
    

      gPad->SetLogy();

      JetShape2[0][ibin][k]->SetMarkerSize(1);
      JetShape2[0][ibin][k]->SetLineColor(kBlack);
      JetShape2[0][ibin][k]->SetMarkerColor(kBlack);
      JetShape2[0][ibin][k]->SetMarkerStyle(10);
     
      JetShape2[1][ibin][k]->SetMarkerSize(1);
      JetShape2[1][ibin][k]->SetLineColor(kBlack);
      JetShape2[1][ibin][k]->SetMarkerColor(kBlack);
      JetShape2[1][ibin][k]->SetMarkerStyle(24);
      

      if(ibin<3){
	JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.);
	JetShape2[0][ibin][k]->GetYaxis()->SetTitleSize(0.);
      }

      gPad->SetLogy();

      JetShape2[0][ibin][k]->Draw();

      
  
      JetShape2[1][0][k]->Draw("same");
      /*
      if(((k>0&&k<5)||k==6)&&!is_subleading){
	JetShape_ref[0][ibin][k]->Draw("same");  



	JetShape_ref[1][0][k]->SetMarkerStyle(20);

	JetShape_ref[1][0][k]->Draw("same");
      }
      */

      if(ibin==3){
	TLegend *legend = new TLegend(0.4,0.55,0.95,0.85);
	//	if(((k>0&&k<5)||k==6)&&!is_subleading)	legend->AddEntry(JetShape_ref[0][ibin][k],"PbPb JetShapes");
	legend->AddEntry(JetShape2[0][ibin][k],"PbPb Present Study");
	//if(((k>0&&k<5)||k==6)&&!is_subleading)	legend->AddEntry(JetShape_ref[1][ibin][k],"pp JetShapes");
	legend->AddEntry(JetShape2[1][ibin][k],"pp Present Study");
	legend->SetLineColor(kWhite);
	legend->SetFillColor(kWhite);
	legend->SetTextSize(0.065);
	legend->Draw();
      }else if(ibin==2){
	if(k<4){
	TLatex  *pttex = new TLatex(0.5,0.81,TrkPtBin_labels[k]);
	pttex->SetNDC();
	pttex->SetTextSize(0.075);
	pttex->Draw();
	}else{
	  TLatex  *pttex = new TLatex(0.5,0.81,"1<p_{T}^{assoc.}<8 GeV");
	  pttex->SetNDC();
	  pttex->SetTextSize(0.075);
	  pttex->Draw();
	}
      }

    
      TLatex  *centtex ;
      centtex = new TLatex(0.55,0.9,CBin_labels[ibin]);
      centtex->SetNDC();
      centtex->SetTextSize(0.075);
      if(ibin==3){  
	centtex = new TLatex(0.6,0.9,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->SetTextSize(0.07);
      }
      
      centtex->Draw();

      c_jetshape[k]->cd(8-ibin);

   
    
      TString diff_name ="JetShape_diff_"; diff_name+=CBin_strs[ibin]; diff_name+= "_"; diff_name += CBin_strs[ibin+1]; diff_name+= k;

      TString diff2_name ="JetShape_diff2_"; diff2_name+=CBin_strs[ibin]; diff2_name+= "_"; diff2_name += CBin_strs[ibin+1]; diff2_name+= k;
      
     
      JetShape_diff[0][ibin][k] = (TH1D*)JetShape[0][ibin][k]->Clone(diff_name);
      JetShape_diff2[0][ibin][k] = (TH1D*)JetShape2[0][ibin][k]->Clone(diff2_name);

      JetShape_diff[0][ibin][k]->Add( JetShape[1][0][k],-1. );
      JetShape_diff2[0][ibin][k]->Add( JetShape2[1][0][k],-1. );
   

      JetShape_diff2[0][ibin][k]->SetMinimum(-2);
      JetShape_diff2[0][ibin][k]->SetMaximum(4.);
     

      JetShape_diff2[0][ibin][k]->SetMarkerStyle(24);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitle("#rho(r)_{PbPb} - #rho(r)_{pp}");
      JetShape_diff2[0][ibin][k]->GetXaxis()->SetNdivisions(408);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetNdivisions(612);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetLabelSize(0.07);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitleSize(0.09);

      if(ibin==3){
	JetShape_diff2[0][ibin][k]->GetXaxis()->SetTitleSize(0.9*JetShape_diff2[0][ibin][k]->GetXaxis()->GetTitleSize());
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitleSize(0.9*JetShape_diff2[0][ibin][k]->GetYaxis()->GetTitleSize());
	JetShape_diff2[0][ibin][k]->GetXaxis()->SetLabelSize(0.9*JetShape_diff2[0][ibin][k]->GetXaxis()->GetLabelSize());
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetLabelSize(0.9*JetShape_diff2[0][ibin][k]->GetYaxis()->GetLabelSize());
      }

      if(ibin<3){
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetLabelSize(0.);
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitleSize(0.);
      }

      JetShape_diff2[0][ibin][k]->Draw();
      /*
      if(((k>0&&k<5)||k==6)&&!is_subleading){
	diff_name ="JetShape_ref_diff_"; diff_name+=CBin_strs[ibin]; diff_name+= "_"; diff_name += CBin_strs[ibin+1]; diff_name+= k;
      
     
	JetShape_ref_diff[0][ibin][k] = (TH1D*)JetShape_ref[0][ibin][k]->Clone(diff_name);

	JetShape_ref_diff[0][ibin][k]->Add( JetShape_ref[1][0][k],-1. );

	JetShape_ref_diff[0][ibin][k]->SetMarkerStyle(20);
	JetShape_ref_diff[0][ibin][k]->Draw("same");
      }

      */

      TString ratio_name ="JetShape_ratio_"; ratio_name+=CBin_strs[ibin]; ratio_name+= "_"; ratio_name += CBin_strs[ibin+1]; ratio_name+= k;

      TString ratio2_name ="JetShape_ratio2_"; ratio2_name+=CBin_strs[ibin]; ratio2_name+= "_"; ratio2_name += CBin_strs[ibin+1]; ratio2_name+= k;
      
     
      JetShape_ratio[0][ibin][k] = (TH1D*)JetShape[0][ibin][k]->Clone(ratio_name);
      JetShape_ratio2[0][ibin][k] = (TH1D*)JetShape2[0][ibin][k]->Clone(ratio2_name);

      JetShape_ratio[0][ibin][k]->Divide( JetShape[1][0][k]);
      JetShape_ratio2[0][ibin][k]->Divide( JetShape2[1][0][k] );
   

      JetShape_ratio2[0][ibin][k]->SetMinimum(-1.);
      JetShape_ratio2[0][ibin][k]->SetMaximum(9.);
     

      JetShape_ratio2[0][ibin][k]->SetMarkerStyle(20);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitle("#rho(r)_{PbPb} - #rho(r)_{pp}");
      JetShape_ratio2[0][ibin][k]->GetXaxis()->SetNdivisions(408);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetNdivisions(612);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetLabelSize(0.07);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitleSize(0.08);

      if(ibin==3){
	JetShape_ratio2[0][ibin][k]->GetXaxis()->SetTitleSize(0.9*JetShape_ratio2[0][ibin][k]->GetXaxis()->GetTitleSize());
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitleSize(0.9*JetShape_ratio2[0][ibin][k]->GetYaxis()->GetTitleSize());
	JetShape_ratio2[0][ibin][k]->GetXaxis()->SetLabelSize(0.9*JetShape_ratio2[0][ibin][k]->GetXaxis()->GetLabelSize());
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetLabelSize(0.9*JetShape_ratio2[0][ibin][k]->GetYaxis()->GetLabelSize());
      }

      if(ibin<3){
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetLabelSize(0.);
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitleSize(0.);
      }

      JetShape_ratio2[0][ibin][k]->Draw();
      /*
      if(((k>0&&k<5)||k==6)&&!is_subleading){
	ratio_name ="JetShape_ref_ratio_"; ratio_name+=CBin_strs[ibin]; ratio_name+= "_"; ratio_name += CBin_strs[ibin+1]; ratio_name+= k;
      
     
	JetShape_ref_ratio[0][ibin][k] = (TH1D*)JetShape_ref[0][ibin][k]->Clone(ratio_name);

	JetShape_ref_ratio[0][ibin][k]->Divide( JetShape_ref[1][0][k] );

	JetShape_ref_ratio[0][ibin][k]->SetMarkerStyle(20);
	JetShape_ref_ratio[0][ibin][k]->Draw("same");
      }
      */
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelSize(0.);

     
      if(ibin==3){
	TLegend *legend = new TLegend(0.2,0.7,0.4,0.85);
	//	if(((k>0&&k<5)||k==6)&&!is_subleading)	legend->AddEntry(JetShape_ref_ratio[0][ibin][k],"PbPb - pp JetShapes");
	legend->AddEntry(JetShape_ratio2[0][ibin][k],"PbPb - pp Present Study");
	legend->SetLineColor(kWhite);
	legend->SetFillColor(kWhite);
	legend->SetTextSize(0.06);
	legend->Draw();
      }      

    }
 

    canvas_name+=".pdf";
    
    c_jetshape[k]->SaveAs(canvas_name);
    
    canvas_name.ReplaceAll(".pdf",".png");
    c_jetshape[k]->SaveAs(canvas_name);

  }

  for(int ibin = 0; ibin<4; ibin++){
 
    if(ibin > 0 )continue;
    

    for(int k = 0; k<7; k++){

 
      JetShape_noerr_up[0][ibin][k] = new TH1D((TString)("JetShape_PbPb_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
      JetShape_noerr_down[0][ibin][k] = new TH1D((TString)("JetShape_PbPb_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
   
      JetShape_noerr_up[1][ibin][k] = new TH1D((TString)("JetShape_pp_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
      JetShape_noerr_down[1][ibin][k] = new TH1D((TString)("JetShape_pp_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
    
      JetShape_diff_noerr_up[0][ibin][k] = new TH1D((TString)("JetShape_Diff_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
      JetShape_diff_noerr_down[0][ibin][k] = new TH1D((TString)("JetShape_Diff_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",19,RBins);  
   
      

      for(int l = 0; l< JetShape[0][ibin][k]->GetNbinsX()+1; l++){
      
	bc = JetShape[0][ibin][k]->GetBinContent(l);
          
	if(bc>0){ 
	  JetShape_noerr_up[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_noerr_down[0][ibin][k]->SetBinContent(l,0.);	

	}else{
	  JetShape_noerr_down[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_noerr_up[0][ibin][k]->SetBinContent(l,0.);	
	}

  
	bc = JetShape_diff[0][ibin][k]->GetBinContent(l);

    
	if(bc>0){ 
	  JetShape_diff_noerr_up[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_diff_noerr_down[0][ibin][k]->SetBinContent(l,0.);	

	}else{
	  JetShape_diff_noerr_down[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_diff_noerr_up[0][ibin][k]->SetBinContent(l,0.);	
	} 

	bc = JetShape[1][ibin][k]->GetBinContent(l);
      
	if(bc>0){ 
	  JetShape_noerr_up[1][ibin][k]->SetBinContent(l,bc);	
	  JetShape_noerr_down[1][ibin][k]->SetBinContent(l,0.);	

	}else{
	  JetShape_noerr_down[1][ibin][k]->SetBinContent(l,bc);	
	  JetShape_noerr_up[1][ibin][k]->SetBinContent(l,0.);	
	}
      }
     
      switch(k){
      case 0:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kBlue-9);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kBlue-9);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kBlue-9);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kBlue-9);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kBlue-9);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kBlue-9);
	break;
      case 1:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kYellow-9);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kYellow-9);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kYellow-9);
	break;
      case 2:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kOrange+1);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kOrange+1);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kOrange+1);
	break;
      case 3:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kViolet-5);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kViolet-5);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kViolet-5);
	break;
      case 4:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kGreen+3);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kGreen+3);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kGreen+3);
	break;
      case 5:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kRed+1);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kRed+1);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kRed+1);
	break;
      default:
	break;
      }

      JetShape_noerr_up[0][ibin][k]->SetFillStyle(1001);
      JetShape_noerr_up[1][ibin][k]->SetFillStyle(1001);
      JetShape_diff_noerr_up[0][ibin][k]->SetFillStyle(1001);

     	  
      JetShape_noerr_down[0][ibin][k]->SetFillStyle(1001);
      JetShape_noerr_down[1][ibin][k]->SetFillStyle(1001);
      JetShape_diff_noerr_down[0][ibin][k]->SetFillStyle(1001);


      
    } //k
    
    cout<<"ready to make stack for "<<ibin<<endl;

    JetShape_Stack_Up[0][ibin]=new THStack((TString)("JetShapeStack_PbPb_Up_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Stack_Up[1][ibin]=new THStack((TString)("JetShapeStack_pp_Up_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Diff_Stack_Up[0][ibin]=new THStack((TString)("JetShapeStack_Diff_Up_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");


    JetShape_Stack_Down[0][ibin]=new THStack((TString)("JetShapeStack_PbPb_Down_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Stack_Down[1][ibin]=new THStack((TString)("JetShapeStack_pp_Down_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Diff_Stack_Down[0][ibin]=new THStack((TString)("JetShapeStack_Diff_Down_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");

    cout<<"made stack"<<endl;

    if(!use_highpT_bin){
      for(int k = 0; k<5; k++){
	JetShape_Stack_Up[0][ibin]->Add(JetShape_noerr_up[0][ibin][4-k]);
	JetShape_Stack_Up[1][ibin]->Add(JetShape_noerr_up[1][ibin][4-k]);
	JetShape_Diff_Stack_Up[0][ibin]->Add(JetShape_diff_noerr_up[0][ibin][4-k]);

	JetShape_Stack_Down[0][ibin]->Add(JetShape_noerr_down[0][ibin][4-k]);
	JetShape_Stack_Down[1][ibin]->Add(JetShape_noerr_down[1][ibin][4-k]);
	JetShape_Diff_Stack_Down[0][ibin]->Add(JetShape_diff_noerr_down[0][ibin][4-k]);
      }

    }else{

      for(int k = 0; k<6; k++){
	JetShape_Stack_Up[0][ibin]->Add(JetShape_noerr_up[0][ibin][k]);
	JetShape_Stack_Up[1][ibin]->Add(JetShape_noerr_up[1][ibin][k]);
	JetShape_Diff_Stack_Up[0][ibin]->Add(JetShape_diff_noerr_up[0][ibin][k]);

	JetShape_Stack_Down[0][ibin]->Add(JetShape_noerr_down[0][ibin][k]);
	JetShape_Stack_Down[1][ibin]->Add(JetShape_noerr_down[1][ibin][k]);
	JetShape_Diff_Stack_Down[0][ibin]->Add(JetShape_diff_noerr_down[0][ibin][k]);
      }


    }

    cout<<"here"<<endl;
  
    JetShape_Stack_Up[0][ibin]->SetMaximum(170.1);
    JetShape_Stack_Up[1][ibin]->SetMaximum(170.1);
    JetShape_Stack_Up[0][ibin]->SetMinimum(-.95);
    JetShape_Stack_Up[1][ibin]->SetMinimum(-.95);
  
  
    
    JetShape_Diff_Stack_Up[0][ibin]->SetMinimum(-1.41);
    JetShape_Diff_Stack_Up[0][ibin]->SetMaximum(1.41);
    
    

    if(use_highpT_bin){

      JetShape_Stack_Up[0][ibin]->SetMaximum(2000.);
      JetShape_Stack_Up[1][ibin]->SetMaximum(2000.);
      JetShape_Stack_Up[0][ibin]->SetMinimum(.05);
      JetShape_Stack_Up[1][ibin]->SetMinimum(.05);
  

      JetShape_Diff_Stack_Up[0][ibin]->SetMinimum(-1.1);
      JetShape_Diff_Stack_Up[0][ibin]->SetMaximum(1.45);
   


    }
    
  


  } //ibin
  

  cout<<"drawing"<<endl;

  TCanvas *PAS_plot = new TCanvas("JetShape_ForPAS","",1200,800);
  PAS_plot->Divide(3,2,0.,0.);


  PAS_plot->cd(5);
  TLegend *legend = new TLegend(0.1,0.05,0.9,.9);
  legend->AddEntry(JetShape_noerr_up[0][0][0],"0.5 < p_{T}^{assoc.}< 1 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][1],"1 < p_{T}^{assoc.}< 2 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][2],"2 < p_{T}^{assoc.}< 3 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][3],"3 < p_{T}^{assoc.}< 4 GeV","f");
  legend->AddEntry(JetShape_noerr_up[0][0][4],"4 < p_{T}^{assoc.}< 8 GeV","f");
 

  if(!use_highpT_bin){
    legend->AddEntry(JetShape_graph[0][0][6],"Total 0.5 < p_{T}^{assoc.}< 8 GeV","lpfe");
  }else{ 
    legend->AddEntry(JetShape_noerr_up[0][0][5],"p_{T}^{assoc.}> 8 GeV","f");
    legend->AddEntry(JetShape_graph[0][0][6],"Total p_{T}^{assoc.}> 0.5 GeV","lpfe");
  }
  legend->AddEntry(JetShape_noerr_up[0][0][4],"|#eta_{track}|< 2.4","");
  legend->SetTextSize(0.055);
  legend->SetLineColor(kWhite);
  legend->SetFillColor(kWhite);
  legend->Draw();

  
    

  PAS_plot->cd(0);
  
  TLatex *type_tex;
 type_tex = new TLatex(0.3,0.96,"Inclusive Jet Shape");
 
  type_tex->SetTextSize(0.035);
  type_tex->SetLineColor(kWhite);
  type_tex->SetNDC();
  type_tex->Draw();
   
  TLatex   *luminosity_tex_pp = new TLatex(0.3,0.92,"pp 5.3 pb^{-1} (2.76 TeV)");
  luminosity_tex_pp->SetTextFont(43);
  luminosity_tex_pp->SetTextSizePixels(25);
    luminosity_tex_pp->SetLineColor(kWhite);
    luminosity_tex_pp->SetNDC();
    luminosity_tex_pp->Draw();
 
    /*  
    TLatex   *luminosity_tex_PbPb = new TLatex(0.3,0.92,"PbPb 166 #mub^{-1} (2.76 TeV)");
    luminosity_tex_PbPb->SetTextFont(43);
    luminosity_tex_PbPb->SetTextSizePixels(25);
    luminosity_tex_PbPb->SetLineColor(kWhite);
    luminosity_tex_PbPb->SetNDC();
    luminosity_tex_PbPb->Draw();
    */

    TLatex   *jet_reco_tex = new TLatex(0.605,0.96,"anti-k_{T}, |#eta_{jet}| < 1.6");
    jet_reco_tex->SetTextFont(43);
    jet_reco_tex->SetTextSizePixels(25);
    jet_reco_tex->SetLineColor(kWhite);
    jet_reco_tex->SetNDC();
    jet_reco_tex->Draw();

    TLatex   *jet_cut_tex = new TLatex(0.605,0.92,"p_{T,1}> 120, p_{T,2}> 50 GeV, #Delta#phi_{1,2}> 5#pi/6");
    jet_cut_tex->SetTextFont(43);
    jet_cut_tex->SetTextSizePixels(25);
    jet_cut_tex->SetLineColor(kWhite);
    jet_cut_tex->SetNDC();
    jet_cut_tex->Draw();
    


  
  PAS_plot->cd(2);
  
 
  JetShape_Stack_Up[1][0]->Draw();
  //  JetShape_Stack_Up[1][0]->GetXaxis()->SetRangeUser(0.,0.99);
  JetShape_Stack_Up[1][0]->GetXaxis()->SetNdivisions(505);

   if(use_highpT_bin)  gPad->SetLogy();
  JetShape_Stack_Up[1][0]->GetYaxis()->SetLabelSize(0.08);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleSize(0.09);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleOffset(0.8);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("#rho(#Deltar)");
  //JetShape_Stack_Up[1][0]->GetYaxis()->CenterTitle();
  JetShape_Stack_Down[1][0]->Draw("same");
  JetShape_Stack_Up[1][0]->Draw("same");

 

  JetShape_graph[1][0][6]->Draw("same e2 P");
  JetShape2[1][0][6]->Draw("same");
  // if(!is_subleading)JetShape_ref[1][0][6]->Draw("same");

  TLatex  *label_pp = new TLatex(0.2,0.9,"pp 5.02 TeV ak4Calo");
  label_pp->SetTextSize(0.09);
  label_pp->SetLineColor(kWhite);
  label_pp->SetNDC();
  label_pp->Draw();

  TLine *l_dr = new TLine(0.,0.,TMath::Pi()/2.,0.);
  l_dr->SetLineStyle(2);
  l_dr->Draw();


  TLatex *cms_tex_dphi = new TLatex(0.25,0.8,"CMS");
  cms_tex_dphi->SetTextFont(63);
  cms_tex_dphi->SetTextSizePixels(30);
  cms_tex_dphi->SetLineColor(kWhite);
  cms_tex_dphi->SetNDC();
  cms_tex_dphi->Draw(); 


TLatex *prelim_tex_dphi = new TLatex(0.45,0.8,"Preliminary");
  prelim_tex_dphi->SetTextFont(53);
  prelim_tex_dphi->SetTextSizePixels(30);
  prelim_tex_dphi->SetLineColor(kWhite);
  prelim_tex_dphi->SetNDC();
  prelim_tex_dphi->Draw(); 


  gPad->RedrawAxis();

  cout<<"and here"<<endl;
  /*
  PAS_plot->cd(2);


  
  JetShape_Stack_Up[0][3]->Draw();

  JetShape_Stack_Up[0][3]->GetXaxis()->SetRangeUser(0.,0.99);
  JetShape_Stack_Up[0][3]->GetXaxis()->SetNdivisions(505);

  gPad->SetLogy();
  JetShape_Stack_Up[0][3]->GetYaxis()->SetLabelSize(0.);
  JetShape_Stack_Down[0][3]->Draw("same");

 

  JetShape_graph[0][3][6]->Draw("same e2 P");
  JetShape2[0][3][6]->Draw("same");

  // if(!is_subleading) JetShape_ref[0][3][6]->Draw("same");

  TLatex  *label_per = new TLatex(0.05,0.9,"PbPb Cent. 50-100%");
  label_per->SetTextSize(0.09);
  label_per->SetLineColor(kWhite);
  label_per->SetNDC();
  label_per->Draw()
;

  */
  PAS_plot->cd(3);
  

  JetShape_Stack_Up[0][0]->GetXaxis()->SetRangeUser(0.,1.);
  JetShape_Stack_Up[0][0]->Draw();
  if(use_highpT_bin) gPad->SetLogy();
  JetShape_Stack_Up[0][0]->GetYaxis()->SetLabelSize(0.);
  JetShape_Stack_Down[0][0]->Draw("same");

  // JetShape_Stack_Up[0][0]->GetXaxis()->SetRangeUser(0.,0.99);
  JetShape_Stack_Up[0][0]->GetXaxis()->SetNdivisions(505);

   
   JetShape_graph[0][0][6]->Draw("same e2 P");
   JetShape2[0][0][6]->Draw("same");

   //if(!is_subleading) JetShape_ref[0][0][6]->Draw("same");
 
  TLatex  *label_cent = new TLatex(0.05,0.9,"pp 2.76 TeV ak3Calo");
  label_cent->SetTextSize(0.09);
  label_cent->SetLineColor(kWhite);
  label_cent->SetNDC();
  label_cent->Draw();

  /*
  PAS_plot->cd(5);


  JetShape_ratio[0][3][6]->SetMinimum(0.);
  JetShape_ratio[0][3][6]->SetMaximum(6.5);


 JetShape_ratio[0][3][6]->Draw();
 JetShape_ratio[0][3][6]->GetYaxis()->SetLabelSize(0.08); 
 JetShape_ratio[0][3][6]->GetYaxis()->CenterTitle();
 JetShape_ratio[0][3][6]->GetYaxis()->SetTitle("#rho(r)_{PbPb}/#rho(r)_{pp}");
 JetShape_ratio[0][3][6]->GetYaxis()->SetTitleSize(0.0);
 JetShape_ratio[0][3][6]->GetYaxis()->SetLabelSize(0.0);
 JetShape_ratio[0][3][6]->GetYaxis()->SetTitleOffset(1.);

 JetShape_ratio[0][3][6]->GetXaxis()->SetLabelSize(0.08); 
 JetShape_ratio[0][3][6]->GetXaxis()->CenterTitle();
 JetShape_ratio[0][3][6]->GetXaxis()->SetTitle("#Deltar");
 JetShape_ratio[0][3][6]->GetXaxis()->SetTitleSize(0.09);
 JetShape_ratio[0][3][6]->GetXaxis()->SetTitleOffset(0.6);
 JetShape_ratio[0][3][6]->GetXaxis()->CenterTitle();
 JetShape_ratio[0][3][6]->GetXaxis()->SetNdivisions(505);
 

 JetShape_ratio[0][3][6]->SetMarkerStyle(20);
 JetShape_ratio[0][3][6]->SetMarkerSize(1);
 JetShape_ratio[0][3][6]->SetMarkerColor(kBlack);
  JetShape_ratio[0][3][6]->Draw();

  JetShape_graph[2][3][6]->Draw("same e2 P");
  JetShape_ratio[0][3][6]->Draw("same");
  // if(!is_subleading) JetShape_ref_ratio[0][3][6]->Draw("same");


  TLatex  *label_ratio = new TLatex(0.05,0.9,"");
  label_ratio->SetTextSize(0.09);
  label_ratio->SetLineColor(kWhite);
  label_ratio->SetNDC();
  label_ratio->Draw();

  TLegend *legend_ratio = new TLegend(0.02,0.7,0.9,0.9);
  legend_ratio->AddEntry(JetShape_ratio[0][3][6],"PbPb/pp");
  //  if(!is_subleading)  legend_ratio->AddEntry(JetShape_ref_ratio[0][3][6],"PLB 730 (2014)");
  legend_ratio->SetLineColor(kWhite);
  legend_ratio->SetFillColor(kWhite);
  legend_ratio->SetTextSize(0.08);
  legend_ratio->Draw("same");


  */
  TLine *line = new TLine(0.,1.,1.,1.);
  line->SetLineStyle(2);
  // line->Draw();

  
   TPave *cover_x_r = new TPave(0.9,0.,1.1,0.13);
   cover_x_r->SetFillColor(kWhite);
   cover_x_r->SetLineColor(kWhite);
   cover_x_r->SetOption("NDC NB");

   //   cover_x_r->Draw();

  TPave *cover_x_l = new TPave(-.0175,0.,0.1,0.125);
   cover_x_l->SetFillColor(kWhite);
   cover_x_l->SetLineColor(kWhite);
   cover_x_l->SetOption("NDC NB");


   //  cover_x_l->Draw();
 


   cout<<"here too"<<endl;
  PAS_plot->cd(6);


  JetShape_ratio[0][0][6]->SetMinimum(0.);
  JetShape_ratio[0][0][6]->SetMaximum(6.5);


 JetShape_ratio[0][0][6]->SetMarkerStyle(20);
 JetShape_ratio[0][0][6]->SetMarkerSize(1);
 JetShape_ratio[0][0][6]->SetMarkerColor(kBlack);
 

 JetShape_ratio[0][0][6]->GetXaxis()->SetRangeUser(0.,1.);
 JetShape_ratio[0][0][6]->Draw();
 JetShape_ratio[0][0][6]->GetYaxis()->SetLabelSize(0.0); 
  JetShape_ratio[0][0][6]->GetXaxis()->SetLabelSize(0.08); 
 JetShape_ratio[0][0][6]->GetXaxis()->CenterTitle();
 JetShape_ratio[0][0][6]->GetXaxis()->SetTitle("#Deltar");
 JetShape_ratio[0][0][6]->GetXaxis()->SetTitleSize(0.09);
 JetShape_ratio[0][0][6]->GetXaxis()->SetTitleOffset(.6);
 JetShape_ratio[0][0][6]->GetXaxis()->CenterTitle();
 JetShape_ratio[0][0][6]->GetXaxis()->SetNdivisions(505);

 
  JetShape_ratio[0][0][6]->Draw();
  JetShape_graph[2][0][6]->Draw("same e2 P");
  JetShape_ratio[0][0][6]->Draw("same");
  // if(!is_subleading) JetShape_ref_ratio[0][0][6]->Draw("same");
 
 TLegend *legend_ratio = new TLegend(0.02,0.7,0.9,0.9);
  legend_ratio->AddEntry(JetShape_ratio[0][0][6],"2.76 TeV / 5.02 TeV");
  legend_ratio->SetLineColor(kWhite);
  legend_ratio->SetFillColor(kWhite);
  legend_ratio->SetTextSize(0.08);
  legend_ratio->Draw("same");



  line->Draw();
 cover_x_l->Draw();
   cover_x_r->Draw();


  PAS_plot->cd(4);
 
  TGaxis *dummy_axis_jetshape = new TGaxis(1.,0.13,1.0,.975,0.,6.5);

  dummy_axis_jetshape->ImportAxisAttributes( JetShape_ratio[0][0][6]->GetYaxis());
  dummy_axis_jetshape->SetTitleOffset(1.5);
  dummy_axis_jetshape->SetTickSize(0.);
  dummy_axis_jetshape->CenterTitle();
  dummy_axis_jetshape->SetTitleSize(0.08);
 dummy_axis_jetshape->SetLabelSize(0.08);

  dummy_axis_jetshape->SetTitle("#rho(#Deltar)_{2.76}/#rho(#Deltar)_{5.02}");
  
 
  TGaxis *dummy_axis_r = new TGaxis(0.,1.,1.,1.,0.,2.);

  dummy_axis_r->ImportAxisAttributes( JetShape_ratio[0][0][6]->GetXaxis());
  dummy_axis_r->SetTitleOffset(0.6);
  dummy_axis_r->SetTitleSize(0.09);
  dummy_axis_r->SetTickSize(0.);
  dummy_axis_r->SetNdivisions(505);
  //dummy_axis_r->Draw();
 
TPave *cover_x_b = new TPave(0.9,0.,0.995,0.17);
   cover_x_b->SetFillColor(kWhite);
   cover_x_b->SetLineColor(kWhite);
   cover_x_b->SetOption("NDC NB");
   cover_x_b->Draw();

   PAS_plot->cd(5);
   dummy_axis_jetshape->Draw();
   dummy_axis_r->Draw();
 

   if(!use_highpT_bin){
     PAS_plot->SaveAs("JetShapes_PAS.pdf");
     PAS_plot->SaveAs("JetShapes_PAS.png");
   }else{
     PAS_plot->SaveAs("JetShapes_WithHighpT_PAS.pdf");
     PAS_plot->SaveAs("JetShapes_WithHighpT_PAS.png");
   }

  
  return 0;
}// main loop
