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

#include "../JetTrack2016_functions.h"

//******************************************************
// VERSION TO DO JET-SHAPES FOR PAS
//******************************************************

Int_t jet_shapes_run1_run2(bool use_highpT_bin = kTRUE, bool is_number = kFALSE){


  //-------------------------
  // Display scale parameters
  //--------------------------

  float rho_max = 8000.;
  float rho_min = 0.1;

  float ratio_min = 0.;
  float ratio_max = 2.;

  




  TFile *fMC[9], *f_ref[7];

  TCanvas *c_jetshape[5];
   
  TString jetetacut, etalabel,centlabel,pTlabel,Ajlabel;
  float eta_ymax;

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins, limR;
  float deta, dphi, r, bc, bg_err,temp1,temp2, rbin, temperr, err, width_temp_x, width_temp_y, width_temp, norm_temp, zerobin, temp, norm_tot, err_temp, cont;
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 10;

  float RBins[20] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.,1.2,1.4,1.6,1.8,2.0};

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  
  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

   float TrkPtBins[nTrkPtBins+1] = {0.5, 1, 2, 3, 4, 8, 12, 16, 20, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_strs2[nTrkPtBins+1] = {"TrkPt07","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.5<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","8<pT<12", "12<pT<16","16<pT<20","pT>20"};
 
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
   
  TH1D* JetShapeMC[6][nCBins][nTrkPtBins];


  TH2D* background[6][nCBins][nPtBins][nTrkPtBins];
  TH1D* background_left[6][nCBins][nPtBins][nTrkPtBins];
  TH1D* background_right[6][nCBins][nPtBins][nTrkPtBins];
  TH1D* background_proj[6][nCBins][nPtBins][nTrkPtBins];


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

  TFile *f_jff_gen = new TFile("../me_correct/Pythia_GenJet_GenTrack_Inclusive_Correlations.root", "READ");
  TFile *f_jff_reco = new TFile("../me_correct/Pythia_RecoJet_GenTrack_Inclusive_Correlations.root", "READ");
 

  THStack *JetShape_Stack_Up[6][nCBins];
  THStack *JetShape_Diff_Stack_Up[6][nCBins];

  THStack *JetShape_Stack_Down[6][nCBins];
  THStack *JetShape_Diff_Stack_Down[6][nCBins];

 
  double temp_cont, nextr, nextl, cos_weight,me00_range, mc_error;
  
  TString stem, datalabel,me00_range_string,stem_mc;

 
 
  float norm, err1;


  TFile *fin, *fin2;

  fin  = new TFile("../../JetTrack2015/bg_fit/Dijet_Correlations.root", "READ");
  
  fin2  = new TFile("../me_correct/pp_Inclusive_Correlations.root", "READ");
 
  TFile *fout  = new TFile("Jet_Shapes_Closures.root","RECREATE"); 

 
 
  //-----------------------
  // Start getting histos
  //-----------------------
  for(int g=0; g<2; g++){

    cout<<"here!"<<endl;

    switch(g){
    case 0:
      if(is_number) stem = "Yield_BkgSub_";
      else  stem = "Yield_BkgSub_pTweighted";
      datalabel = "Inclusive";
      break;
    case 1:
      if(is_number) stem = "Yield_BkgSub_";
      else stem = "Yield_BkgSub_pTweighted";
      //stem = "Yield_BkgSub_";
      datalabel = "Inclusive";
      break;
    }

    for (int ibin=0;ibin<nCBins;ibin++){
  
      if(ibin > 0)continue;

      for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
    
	for (int ibin3=0;ibin3<nTrkPtBins-1;ibin3++){

	  if(g==0){

	    cout<<stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]<<endl;


	    if(ibin3<5)     result[g][ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)("Yield_BkgSub_pp_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + datalabel+"_Run1_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    else  result[g][ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)("Yield_BkgSub_pp_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_" + TrkPtBin_strs2[5] + "_" + TrkPtBin_strs[9]))->Clone((TString) (stem + datalabel+"_Run1_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    cout<<"here"<<endl;
	      
	    if(ibin3>5){

	      for(int k = 0; k< result[g][ibin][ibin2][ibin3]->GetNbinsX()+1; k++){
		for(int m = 0; m< result[g][ibin][ibin2][ibin3]->GetNbinsY()+1; m++){
		  result[g][ibin][ibin2][ibin3]->SetBinContent(k,m,0.);
		  result[g][ibin][ibin2][ibin3]->SetBinError(k,m,0.);
		}
	
	      }
	    }

	  }else{

	    result[g][ibin][ibin2][ibin3] = (TH2D*)fin2->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_" + TrkPtBin_strs2[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + datalabel+"_Run2_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    
	  }

	  width_temp_x = result[g][ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1);
	  width_temp_y = result[g][ibin][ibin2][ibin3]->GetYaxis()->GetBinWidth(1);
	
	  /*
	    if(ibin==0&&g==0){ //recoreco pythia

	      cout<<"here 0"<<endl;
	      resultMC_gen[g][ibin][0][ibin3] = (TH2D*)fin2->Get((TString)(stem+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"GenGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      cout<<"here"<<endl;
	      resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_reco->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)(stem+"RecoGen_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	      resultMC_reco[g][ibin][0][ibin3]->Add(   resultMC_gen[g][ibin][0][ibin3],-1.);
	      result[g][ibin][ibin2][ibin3]->Add(resultMC_reco[g][0][0][ibin3],-1.);
	  
	      cout<<"and here"<<endl;
	    }
	    cout<<"dealt with corrections"<<endl;
	  */
	  
	
	
	  if(ibin3==0){
	    result[g][ibin][ibin2][9] = (TH2D*) result[g][ibin][ibin2][ibin3]->Clone(((TString) ("Summed_result_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+TrkPtBin_strs[0]+"_" +TrkPtBin_strs[10])));

	  }else if(ibin3>0&&(g==1||ibin3<6)){

	    result[g][ibin][ibin2][9]->Add(result[g][ibin][ibin2][ibin3]);
	 
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


	JetShape[g][ibin][ibin3] = (TH1D*)JetShape2[g][ibin][ibin3]->Clone((TString)("JetShape_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	

	
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
      
	if(ibin3==9){

	  cout<<"norm is "<<norm<<endl;

	  //	    if(!is_subleading&&!use_highpT_bin) norm = JetShape[g][ibin][ibin3]->GetBinContent(1)/JetShape_ref[g][ibin][ibin3]->GetBinContent(1);	 

	  JetShape[g][ibin][ibin3]->Scale(1./norm);  
	  JetShape[g][ibin][0]->Scale(1./norm);  
	  JetShape[g][ibin][1]->Scale(1./norm);  
	  JetShape[g][ibin][2]->Scale(1./norm);  
	  JetShape[g][ibin][3]->Scale(1./norm);  
	  JetShape[g][ibin][4]->Scale(1./norm);  
	  JetShape[g][ibin][5]->Scale(1./norm);  
	  JetShape[g][ibin][6]->Scale(1./norm);  
	  JetShape[g][ibin][7]->Scale(1./norm);  
	  JetShape[g][ibin][8]->Scale(1./norm);  

	 
	}
	  
	  
	if(ibin<3){
	  JetShape2[g][ibin][ibin3]->GetYaxis()->SetLabelSize(0.);
	  JetShape2[g][ibin][ibin3]->GetYaxis()->SetTitleSize(0.);

	}
	 	
	fout->cd();
	  

	JetShape2[g][ibin][ibin3]->Write();


	cout<<"here"<<endl;
      }
	  
    }

  }

  cout<<"done "<<g<<endl;
}
 

cout<<"ready to draw"<<endl;

for(int k = 0; k<10; k++){
   
  TString canvas_name = "JetShape_Closure_";
  
  if(k<9){ 
    canvas_name+=TrkPtBin_strs[k]; canvas_name+="_"; canvas_name+=TrkPtBin_strs[k+1];  
  } else {
    canvas_name+=TrkPtBin_strs[1]; canvas_name+="_"; canvas_name+=TrkPtBin_strs[10];  
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
    JetShape2[0][ibin][k]->GetYaxis()->SetTitle("#rho(r)");
    JetShape2[0][ibin][k]->GetYaxis()->CenterTitle(true);
    JetShape2[0][ibin][k]->GetYaxis()->SetLabelFont(42);
    JetShape2[0][ibin][k]->GetYaxis()->SetLabelOffset(0.004);
    JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.075);
    JetShape2[0][ibin][k]->GetYaxis()->SetTitleSize(0.09);
    JetShape2[0][ibin][k]->GetYaxis()->SetTickLength(0.025);
    JetShape2[0][ibin][k]->GetYaxis()->SetTitleOffset(1.);
     
    // JetShape2[0][ibin][k]->SetMinimum(0.008);
    //JetShape2[0][ibin][k]->SetMaximum(14.5);

    JetShape2[0][ibin][k]->SetMinimum(0.1);
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
   

    //  JetShape_ratio[0][ibin][k]->Add( JetShape[1][0][k]);
    // JetShape_ratio2[0][ibin][k]->Add( JetShape2[1][0][k] );
   

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
    

  for(int k = 0; k<9; k++){

 
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
	JetShape_noerr_up[1][0][k]->SetFillColor(kBlue-9);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kBlue-9);
	JetShape_noerr_down[1][0][k]->SetFillColor(kBlue-9);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kBlue-9);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kBlue-9);
	break;
      case 1:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_noerr_up[1][0][k]->SetFillColor(kYellow-9);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_noerr_down[1][0][k]->SetFillColor(kYellow-9);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kYellow-9);
	break;
      case 2:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_noerr_up[1][0][k]->SetFillColor(kOrange+1);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_noerr_down[1][0][k]->SetFillColor(kOrange+1);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kOrange+1);
	break;
      case 3:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_noerr_up[1][0][k]->SetFillColor(kViolet-5);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_noerr_down[1][0][k]->SetFillColor(kViolet-5);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kViolet-5);
	break;
      case 4:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_noerr_up[1][0][k]->SetFillColor(kGreen+3);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_noerr_down[1][0][k]->SetFillColor(kGreen+3);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kGreen+3);
	break;


      case 5:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kRed);
	JetShape_noerr_up[1][0][k]->SetFillColor(kRed);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kRed);
	JetShape_noerr_down[1][0][k]->SetFillColor(kRed);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kRed);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kRed);
	break;

    
      case 6:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_noerr_up[1][0][k]->SetFillColor(kRed+1);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_noerr_down[1][0][k]->SetFillColor(kRed+1);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kRed+1);
	break;

      case 7:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kRed+2);
	JetShape_noerr_up[1][0][k]->SetFillColor(kRed+2);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kRed+2);
	JetShape_noerr_down[1][0][k]->SetFillColor(kRed+2);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kRed+2);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kRed+2);
	break;

      case 8:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kRed+3);
	JetShape_noerr_up[1][0][k]->SetFillColor(kRed+3);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kRed+3);
	JetShape_noerr_down[1][0][k]->SetFillColor(kRed+3);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kRed+3);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kRed+3);
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
    for(int k = 0; k<8; k++){
      JetShape_Stack_Up[0][ibin]->Add(JetShape_noerr_up[0][ibin][4-k]);
      JetShape_Stack_Up[1][ibin]->Add(JetShape_noerr_up[1][ibin][4-k]);
      JetShape_Diff_Stack_Up[0][ibin]->Add(JetShape_diff_noerr_up[0][ibin][4-k]);

      JetShape_Stack_Down[0][ibin]->Add(JetShape_noerr_down[0][ibin][4-k]);
      JetShape_Stack_Down[1][ibin]->Add(JetShape_noerr_down[1][ibin][4-k]);
      JetShape_Diff_Stack_Down[0][ibin]->Add(JetShape_diff_noerr_down[0][ibin][4-k]);
    }

  }else{

    for(int k = 0; k<9; k++){
      JetShape_Stack_Up[0][ibin]->Add(JetShape_noerr_up[0][ibin][k]);
      JetShape_Stack_Up[1][ibin]->Add(JetShape_noerr_up[1][ibin][k]);
      JetShape_Diff_Stack_Up[0][ibin]->Add(JetShape_diff_noerr_up[0][ibin][k]);

      JetShape_Stack_Down[0][ibin]->Add(JetShape_noerr_down[0][ibin][k]);
      JetShape_Stack_Down[1][ibin]->Add(JetShape_noerr_down[1][ibin][k]);
      JetShape_Diff_Stack_Down[0][ibin]->Add(JetShape_diff_noerr_down[0][ibin][k]);
    }


  }

  cout<<"here"<<endl;
    
  if(use_highpT_bin){

    JetShape_Stack_Up[0][ibin]->SetMaximum(rho_max);
    JetShape_Stack_Up[1][ibin]->SetMaximum(rho_max);
    JetShape_Stack_Up[0][ibin]->SetMinimum(rho_min);
    JetShape_Stack_Up[1][ibin]->SetMinimum(rho_min);
      
  }
    
  


 } //ibin
  

cout<<"drawing"<<endl;

TCanvas *PAS_plot = new TCanvas("JetShape_ForPAS","",1200,800);
PAS_plot->Divide(3,2,0.,0.);


PAS_plot->cd(4);
TLegend *legend = new TLegend(0.1,0.05,0.9,.9);
legend->AddEntry(JetShape_noerr_up[0][0][0],"0.5 < p_{T}^{assoc.}< 1 GeV","f");
legend->AddEntry(JetShape_noerr_up[0][0][1],"1 < p_{T}^{assoc.}< 2 GeV","f");
legend->AddEntry(JetShape_noerr_up[0][0][2],"2 < p_{T}^{assoc.}< 3 GeV","f");
legend->AddEntry(JetShape_noerr_up[0][0][3],"3 < p_{T}^{assoc.}< 4 GeV","f");
legend->AddEntry(JetShape_noerr_up[0][0][4],"4 < p_{T}^{assoc.}< 8 GeV","f");
 
if(use_highpT_bin){
  legend->AddEntry(JetShape_noerr_up[0][0][5],"p_{T}^{assoc.}> 8 GeV","f");
    
 }
 
legend->AddEntry(JetShape_noerr_up[0][0][4],"|#eta_{track}|< 2.4","");
legend->SetTextSize(0.055);
legend->SetLineColor(kWhite);
legend->SetFillColor(kWhite);
legend->Draw();

  
    

PAS_plot->cd(0);

TLatex *type_tex = new TLatex(0.3,0.96,"Inclusive Jet Shape");
 
type_tex->SetTextSize(0.035);
type_tex->SetLineColor(kWhite);
type_tex->SetNDC();
type_tex->Draw();
   
 TLatex   *luminosity_tex_pp;

 luminosity_tex_pp = new TLatex(0.3,0.92,"pp 5.02 TeV");
 luminosity_tex_pp->SetTextFont(43);
 luminosity_tex_pp->SetTextSizePixels(25);
 luminosity_tex_pp->SetLineColor(kWhite);
 luminosity_tex_pp->SetNDC();
 luminosity_tex_pp->Draw();

  
TLatex   *jet_reco_tex = new TLatex(0.605,0.96,"anti-k_{T}, |#eta_{jet}| < 1.6");
jet_reco_tex->SetTextFont(43);
jet_reco_tex->SetTextSizePixels(25);
jet_reco_tex->SetLineColor(kWhite);
jet_reco_tex->SetNDC();
jet_reco_tex->Draw();

 TLatex   *jet_cut_tex = new TLatex(0.605,0.92,"p_{T} > 120");
jet_cut_tex->SetTextFont(43);
jet_cut_tex->SetTextSizePixels(25);
jet_cut_tex->SetLineColor(kWhite);
jet_cut_tex->SetNDC();
jet_cut_tex->Draw();
    


  
PAS_plot->cd(2);
  
 
JetShape_Stack_Up[1][0]->Draw();
 JetShape_Stack_Up[1][0]->GetXaxis()->SetRangeUser(0.,01.499);
JetShape_Stack_Up[1][0]->GetXaxis()->SetNdivisions(505);


if(use_highpT_bin)  gPad->SetLogy();
JetShape_Stack_Up[1][0]->GetYaxis()->SetLabelSize(0.08);
JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleSize(0.09);
JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleOffset(0.8);
JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("#rho(#Deltar)  (GeV)");
//JetShape_Stack_Up[1][0]->GetYaxis()->CenterTitle();
JetShape_Stack_Down[1][0]->Draw("same");
JetShape_Stack_Up[1][0]->Draw("same");

 
JetShape2[1][0][9]->Draw("same");
// if(!is_subleading)JetShape_ref[1][0][9]->Draw("same");

TLatex  *label_pp = new TLatex(0.1,0.9,"pp 5.02 TeV");
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
 
PAS_plot->cd(3);
  
JetShape_Stack_Up[0][0]->Draw();
if(use_highpT_bin) gPad->SetLogy();
JetShape_Stack_Up[0][0]->GetYaxis()->SetLabelSize(0.);
JetShape_Stack_Down[0][0]->Draw("same");

 JetShape_Stack_Up[0][0]->GetXaxis()->SetRangeUser(0.,01.499);
JetShape_Stack_Up[0][0]->GetXaxis()->SetNdivisions(505);

   
JetShape2[0][0][9]->Draw("same");

//if(!is_subleading) JetShape_ref[0][0][9]->Draw("same");
 
 TLatex  *label_cent;
 label_cent = new TLatex(0.05,0.9,"pp 2.76 TeV");
 label_cent->SetTextSize(0.09);
label_cent->SetLineColor(kWhite);
label_cent->SetNDC();
label_cent->Draw();


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


JetShape_ratio[0][0][9]->SetMinimum(0.);
JetShape_ratio[0][0][9]->SetMaximum(ratio_max);


JetShape_ratio[0][0][9]->SetMarkerStyle(20);
JetShape_ratio[0][0][9]->SetMarkerSize(1);
JetShape_ratio[0][0][9]->SetMarkerColor(kBlack);

 JetShape_ratio[0][0][9]->GetXaxis()->SetRangeUser(0.,01.499);
 

JetShape_ratio[0][0][9]->Draw();
JetShape_ratio[0][0][9]->GetYaxis()->SetLabelSize(0.0); 
JetShape_ratio[0][0][9]->GetXaxis()->SetLabelSize(0.08); 
JetShape_ratio[0][0][9]->GetXaxis()->CenterTitle();
JetShape_ratio[0][0][9]->GetXaxis()->SetTitle("#Deltar");
JetShape_ratio[0][0][9]->GetXaxis()->SetTitleSize(0.09);
JetShape_ratio[0][0][9]->GetXaxis()->SetTitleOffset(.6);
JetShape_ratio[0][0][9]->GetXaxis()->CenterTitle();
JetShape_ratio[0][0][9]->GetXaxis()->SetNdivisions(505);

 
JetShape_ratio[0][0][9]->Draw();

JetShape_ratio[0][0][9]->Draw("same");
// if(!is_subleading) JetShape_ref_ratio[0][0][9]->Draw("same");
 
TLegend *legend_ratio = new TLegend(0.02,0.7,0.9,0.9);
legend_ratio->AddEntry(JetShape_ratio[0][0][9],"2.76 TeV / 5.02 TeV");
legend_ratio->SetLineColor(kWhite);
legend_ratio->SetFillColor(kWhite);
legend_ratio->SetTextSize(0.08);
legend_ratio->Draw("same");



line->Draw();
cover_x_l->Draw();
cover_x_r->Draw();


PAS_plot->cd(4);
 
TGaxis *dummy_axis_jetshape = new TGaxis(1.,0.13,1.0,.975,0.,ratio_max);

dummy_axis_jetshape->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetYaxis());
dummy_axis_jetshape->SetTitleOffset(1.5);
dummy_axis_jetshape->SetTickSize(0.);
dummy_axis_jetshape->CenterTitle();
dummy_axis_jetshape->SetTitleSize(0.08);
dummy_axis_jetshape->SetLabelSize(0.08);

 dummy_axis_jetshape->SetTitle("#rho(#Deltar)_{2.76}/#rho(#Deltar)_{5.02}");
   
 
TGaxis *dummy_axis_r = new TGaxis(0.,1.,1.,1.,0.,1.5);

dummy_axis_r->ImportAxisAttributes( JetShape_ratio[0][0][9]->GetXaxis());
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
 

 PAS_plot->SaveAs("JetShapes_WithHighpT_Run1_Run2.pdf");
 PAS_plot->SaveAs("JetShapes_WithHighpT_Run1_Run2.png");
 
  
return 0;
}// main loop
