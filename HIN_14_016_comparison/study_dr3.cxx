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
#include "TStyle.h"
#include "TLatex.h"


#include <iostream>
#include <vector>
#include <fstream>

#include "../HIN_14_016_functions.h"

Int_t study_dr3(){

  TFile *fin[6], *fMC[6];
   
  TString jetetacut, etalabel;
  float eta_ymax;

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins, limR;
  float deta, dphi, r, bc, temp1,temp2, rbin, temperr, err, width_temp_x, width_temp_y, width_temp, norm_temp, zerobin, temp, norm_tot, err_temp, cont;
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 5;

  float RBins[13] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.6,0.8,1.,1.4};

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt999" };
  TString TrkPtBin_labels[nTrkPtBins] = {"1<pT<2","2<pT<3","3<pT<4","4<pT<8","pT>8"};
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
 
 
  TH2D* result[6][nCBins][nPtBins][nTrkPtBins];
  TH2D* resultMC[6][nCBins][nPtBins][nTrkPtBins];
  TH1D* R1[6][nCBins][nPtBins];
  TH1D* R2[6][nCBins][nPtBins];
  TH1D* EtaProj[6][nCBins][nPtBins];
  TH1D* PhiProj[6][nCBins][nPtBins];
  TH1D* JetShape[6][nCBins][nPtBins];
  TGraphErrors* JetShape_graph[6][nCBins][nPtBins];
  TH1D* Ratio[6][nCBins][nPtBins];
  TGraphErrors* Running_int_graph[6][nCBins][nPtBins];
  TH1D* Running_int[6][nCBins][nPtBins];

  TH2D* Raw_Yield_TrkPt1_TrkPt8[6][nCBins][nPtBins];
  TH2D* Closure_2D_TrkPt1_TrkPt8[6][nCBins][nPtBins];
  TF1* FitR2 = new TF1("FitR2","[0]+[1]*x",0.,1.);

  TF2* ClosureFit = new TF2("ClosureFit", "[0]/2/TMath::Pi()/[1]/[2]*TMath::Exp(-1.*(x*x/[1]/[1]/2))*TMath::Exp(-1.*(y*y/[2]/[2]/2))",-1.,1.,-1.,1.);

  TLegend *JetLegend, CorrelationLegend,*lproj;
  TString legendentry;
  TLatex *centtex;
   
  TCanvas *AllJetShapesCanvas = new TCanvas((TString)("AllJetShapesCanvas")," ",10,10,1500,1200);
  AllJetShapesCanvas->Divide(4,2,0.00000,0.0000000);
 
  TCanvas *ClosureTestInclusive = new TCanvas((TString)("ClosureTestInclusive")," ",10,10,1500,1200);
  ClosureTestInclusive->Divide(4,3,0.0001,0.00001);
  
  TCanvas *AllJetShapesDisplay = new TCanvas((TString)("AllJetShapesDisplay")," ",10,10,1500,600);
  AllJetShapesDisplay->Divide(2,1,0.00000,0.0000000);
 

 
  TString stem, datalabel;
  fin[0] = new TFile("../bg_fits/OfficialRecoPbPb_Inclusive120_2Dyield_and_NewBkg_files.root", "READ");
  fMC[0] = new TFile("../bg_fits/HYDJET_Inc120_2Dyield_and_NewBkg_files.root","READ");
  fin[1] = new TFile("../bg_fits/pp_Inclusive120_2Dyield_and_NewBkg_files.root", "READ");
  fMC[1] = new TFile("../bg_fits/PYTHIA_Inc120_2Dyield_and_NewBkg_files.root","READ");
  fin[2] = new TFile("../bg_fits/OfficialRecoPbPb_SubLeading50_2Dyield_and_NewBkg_files.root", "READ");
  fMC[2] = new TFile("../bg_fits/HYDJET_Subleading50_2Dyield_and_NewBkg_files.root","READ");
  fin[3] = new TFile("../bg_fits/pp_SubLeading50_2Dyield_and_NewBkg_files.root", "READ");
  fMC[3] = new TFile("../bg_fits/PYTHIA_Subleading50_2Dyield_and_NewBkg_files.root","READ");
  fin[4] = new TFile("../bg_fits/OfficialRecoPbPb_Leading120_2Dyield_and_NewBkg_files.root", "READ");
  fMC[4] = new TFile("../bg_fits/HYDJET_Leading120_2Dyield_and_NewBkg_files.root","READ");
  fin[5] = new TFile("../bg_fits/pp_Leading120_2Dyield_and_NewBkg_files.root", "READ");
  fMC[5] = new TFile("../bg_fits/PYTHIA_Leading120_2Dyield_and_NewBkg_files.root","READ");
  
 
   TFile *fout = new TFile("PbPb_dR_Correlations.root","RECREATE"); 

    //-----------------------
    // Start getting histos
    //-----------------------
 for(int g=0; g<6; g++){

   switch(g){
   case 0:
     stem = "Result_PbPb_";
     datalabel = "Inclusive";
     break;
   case 1:
     stem = "Result_pp_";
     datalabel = "Inclusive";
     break;
   case 2:
     stem = "Result_PbPb_";
     datalabel = "SubLeading";
     break;
   case 3:
     stem = "Result_pp_";
     datalabel = "SubLeading";
     break;
   case 4:
     stem = "Result_PbPb_";
     datalabel = "Leading";
     break;
   case 5:
     stem = "Result_pp_";
     datalabel = "Leading";
     break;
   }
   
    for (int ibin=0;ibin<nCBins;ibin++){
  
      for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
    
	for (int ibin3=0;ibin3<nTrkPtBins-1;ibin3++){


	 
	  result[g][ibin][ibin2][ibin3] = (TH2D*) fin[g]->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	
	 
	  width_temp_x = result[g][ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1);
	  width_temp_y = result[g][ibin][ibin2][ibin3]->GetYaxis()->GetBinWidth(1);

	  result[g][ibin][ibin2][ibin3]->Scale(1./width_temp_x/width_temp_y);

	
	  if((g==0||g==2||g==4)) {

	    resultMC[g][ibin][ibin2][ibin3] = (TH2D*) fMC[g]->Get((TString)("SummedResult_PYTHIAHYDJET_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("SummedResult_PYTHIAHYDJET_"+ stem + datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    resultMC[g][ibin][ibin2][ibin3]->Scale(1./width_temp_x/width_temp_y);
	  }

	
	
	  if(g==1||g==3||g==5) {
	    
	    resultMC[g][ibin][ibin2][ibin3] = (TH2D*) fMC[g]->Get((TString)("SummedResult_PYTHIA_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("SummedResult_PYTHIA_"+ stem + datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	    resultMC[g][ibin][ibin2][ibin3]->Scale(1./width_temp_x/width_temp_y);
	  }

	} 
	
	R1[g][ibin][ibin2] = new TH1D((TString)("R1_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]),"",12,RBins);  
	R2[g][ibin][ibin2] = new TH1D((TString)("R2_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]),"",12,RBins);  

	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2] = (TH2D*) result[g][ibin][ibin2][0]->Clone((TString) ("Summed_" + stem + datalabel + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[0]+"_" +TrkPtBin_strs[4]));     
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->Add( result[g][ibin][ibin2][1]);
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->Add( result[g][ibin][ibin2][2]);
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->Add( result[g][ibin][ibin2][3]);
     
	if((g==0||g==2||g==4)||ibin==0){

	  Closure_2D_TrkPt1_TrkPt8[g][ibin][ibin2] = (TH2D*) resultMC[g][ibin][ibin2][0]->Clone((TString) ("Closure_2D_"+stem+ datalabel + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[0]+"_" +TrkPtBin_strs[4]));     
	
	  Closure_2D_TrkPt1_TrkPt8[g][ibin][ibin2]->Add( resultMC[g][ibin][ibin2][1]);
	  Closure_2D_TrkPt1_TrkPt8[g][ibin][ibin2]->Add( resultMC[g][ibin][ibin2][2]);
	  Closure_2D_TrkPt1_TrkPt8[g][ibin][ibin2]->Add( resultMC[g][ibin][ibin2][3]);
	}
      }
    }
 }

 for (int ibin=0;ibin<nCBins;ibin++){
   for (int ibin2=0;ibin2<nPtBins;ibin2++){ 

     Closure_2D_TrkPt1_TrkPt8[0][ibin][ibin2]->Add(Closure_2D_TrkPt1_TrkPt8[1][0][ibin2],-1.);
     Closure_2D_TrkPt1_TrkPt8[2][ibin][ibin2]->Add(Closure_2D_TrkPt1_TrkPt8[3][0][ibin2],-1.);
     Closure_2D_TrkPt1_TrkPt8[4][ibin][ibin2]->Add(Closure_2D_TrkPt1_TrkPt8[5][0][ibin2],-1.);
     
     //Constrain and fit
     int kmin = Closure_2D_TrkPt1_TrkPt8[0][ibin][ibin2]->GetXaxis()->FindBin(-1+.0001);
     int kmax = Closure_2D_TrkPt1_TrkPt8[0][ibin][ibin2]->GetXaxis()->FindBin(1.-.0001);

     int lmin = Closure_2D_TrkPt1_TrkPt8[0][ibin][ibin2]->GetYaxis()->FindBin(-1.+.0001);
     int lmax = Closure_2D_TrkPt1_TrkPt8[0][ibin][ibin2]->GetYaxis()->FindBin(1.-.0001);
     
     ClosureFit->SetParLimits(1,0.,1.);
     ClosureFit->SetParLimits(2,0.,1.);
      

     ClosureTestInclusive->cd(4-ibin);
     temp = 0.;
     for(int k = kmin; k<kmax+1; k++){  
       for(int l = lmin; l<lmax+1; l++){
	 temp+= Closure_2D_TrkPt1_TrkPt8[0][ibin][ibin2]->GetBinContent(k,l);
       }
     }
     ClosureFit->FixParameter(0,temp);
     Closure_2D_TrkPt1_TrkPt8[0][ibin][ibin2]->Fit(ClosureFit);
     //  Raw_Yield_TrkPt1_TrkPt8[0][ibin][ibin2]->Add(ClosureFit,-1.);

     Closure_2D_TrkPt1_TrkPt8[0][ibin][ibin2]->GetXaxis()->SetRangeUser(-1.0,1.0);
     Closure_2D_TrkPt1_TrkPt8[0][ibin][ibin2]->GetYaxis()->SetRangeUser(-1.0,1.0);
     Closure_2D_TrkPt1_TrkPt8[0][ibin][ibin2]->Draw("surf1");

     ClosureTestInclusive->cd(12-ibin);

     temp = 0.;
     for(int k = kmin; k<kmax+1; k++){  
       for(int l = lmin; l<lmax; l++){
	 temp+= Closure_2D_TrkPt1_TrkPt8[2][ibin][ibin2]->GetBinContent(k,l);
       }
     }
          ClosureFit->FixParameter(0,temp);
     Closure_2D_TrkPt1_TrkPt8[2][ibin][ibin2]->Fit(ClosureFit);
     //  Raw_Yield_TrkPt1_TrkPt8[2][ibin][ibin2]->Add(ClosureFit,-1.);


     Closure_2D_TrkPt1_TrkPt8[2][ibin][ibin2]->GetXaxis()->SetRangeUser(-1.0,1.0);
     Closure_2D_TrkPt1_TrkPt8[2][ibin][ibin2]->GetYaxis()->SetRangeUser(-1.0,1.0);
     Closure_2D_TrkPt1_TrkPt8[2][ibin][ibin2]->Draw("surf1");     

   ClosureTestInclusive->cd(8-ibin);

     temp = 0.;
     for(int k = kmin; k<kmax+1; k++){  
       for(int l = lmin; l<lmax; l++){
	 temp+= Closure_2D_TrkPt1_TrkPt8[4][ibin][ibin2]->GetBinContent(k,l);
       }
     }
     ClosureFit->FixParameter(0,temp);
     Closure_2D_TrkPt1_TrkPt8[4][ibin][ibin2]->Fit(ClosureFit);
     //   Raw_Yield_TrkPt1_TrkPt8[4][ibin][ibin2]->Add(ClosureFit,-1.);


     Closure_2D_TrkPt1_TrkPt8[4][ibin][ibin2]->GetXaxis()->SetRangeUser(-1.0,1.0);
     Closure_2D_TrkPt1_TrkPt8[4][ibin][ibin2]->GetYaxis()->SetRangeUser(-1.0,1.0);
     Closure_2D_TrkPt1_TrkPt8[4][ibin][ibin2]->Draw("surf1");
        
   }
 }

 //------------------------------------------------------
 //  dR calculation on summed, spill-over corrected 2Ds
 //------------------------------------------------------

for(int g=0; g<6; g++){

   switch(g){
   case 0:
     stem = "Result_PbPb_";
     datalabel = "Inclusive";
     break;
   case 1:
     stem = "Result_pp_";
     datalabel = "Inclusive";
     break;
   case 2:
     stem = "Result_PbPb_";
     datalabel = "SubLeading";
     break;
   case 3:
     stem = "Result_pp_";
     datalabel = "SubLeading";
     break;
   case 4:
     stem = "Result_PbPb_";
     datalabel = "Leading";
     break;
   case 5:
     stem = "Result_pp_";
     datalabel = "Leading";
     break;
   }

  for (int ibin=0;ibin<nCBins;ibin++){
    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
   
	int kmin = Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetXaxis()->FindBin(-4+.0001);
	int kmax = Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetXaxis()->FindBin(4.-.0001);


	int lmin = Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetYaxis()->FindBin(-TMath::Pi()+.0001);
	int lmax = Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetYaxis()->FindBin(TMath::Pi()-.0001);
      
	for(int k = kmin; k<kmax+1; k++){  //
	  for(int l = lmin; l<lmax+1; l++){

	    deta = Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetXaxis()->GetBinCenter(k);
	    dphi =  Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetYaxis()->GetBinCenter(l);

	   
	    r = TMath::Sqrt(deta*deta+dphi*dphi);
	  
	    rbin = R1[g][ibin][ibin2]->FindBin(r);

	    bc = Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetBinContent(k,l);
	    err = Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetBinError(k,l);
	  
	  
	    temp1 = R1[g][ibin][ibin2]->GetBinContent(rbin);
	    temp2 = R2[g][ibin][ibin2]->GetBinContent(rbin);

	    temperr = R1[g][ibin][ibin2]->GetBinError(rbin);
	  
	    R1[g][ibin][ibin2]->SetBinContent(rbin,temp1+bc);
	    R1[g][ibin][ibin2]->SetBinError(rbin,TMath::Sqrt(temperr*temperr+err*err));
	
	    R2[g][ibin][ibin2]->SetBinContent(rbin,temp2+1);
	  
	
	  }

	}

	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetXaxis()->SetTitle("#eta");
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetYaxis()->SetTitle("#phi");
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetXaxis()->SetTitleSize(0.08);
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetYaxis()->SetTitleSize(0.08);
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetXaxis()->SetTitleOffset(1.);
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetYaxis()->SetTitleOffset(1.);
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetXaxis()->SetRangeUser(-1.999,1.999);
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->GetYaxis()->SetRangeUser(-TMath::Pi()/2+0.0001,3*TMath::Pi()/2-0.0001);
	Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->SetMaximum(15.);


	JetShape[g][ibin][ibin2]= (TH1D*)R1[g][ibin][ibin2]->Clone((TString)("Jet_Shape_"+stem+datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]));

	Running_int[g][ibin][ibin2] = (TH1D*)R1[g][ibin][ibin2]->Clone((TString)("Running_int_"+ stem+ datalabel+"_"+CBin_strs[ibin] + "_" + CBin_strs[ibin+1])); 


	for(int k = 0; k<24; k++){

	  if(k<2){
	    JetShape[g][ibin][ibin2]->SetBinContent(k,0);
	    JetShape[g][ibin][ibin2]->SetBinError(k,0);
	    continue;
	  }
	  temp =  R2[g][ibin][ibin2]->GetBinContent(k);
	  if(temp==0){continue;}
	  bc = JetShape[g][ibin][ibin2]->GetBinContent(k)/temp;
	  err = JetShape[g][ibin][ibin2]->GetBinError(k)/temp;
	  JetShape[g][ibin][ibin2]->SetBinContent(k,bc);
	  JetShape[g][ibin][ibin2]->SetBinError(k,err);
	
	  if(g==1||g==3||g==5){

	    JetShape[g][ibin][ibin2]->SetBinError(k,0.05*bc);
	    
	  }
	}

	JetShape[g][ibin][ibin2]->SetMarkerStyle(10);
	JetShape[g][ibin][ibin2]->SetMarkerSize(2);
	JetShape[g][ibin][ibin2]->SetLineColor(kBlack);
	JetShape[g][ibin][ibin2]->SetMinimum(0.01);
	JetShape[g][ibin][ibin2]->SetAxisRange(0.,0.999,"x");
	/*

	if(g==1||g==3||g==5){
	  JetShape[g][ibin][ibin2]->SetMarkerStyle(4);
	}
	*/
      } // ibin2
   		       
    } // centrality bin loop
  
  } // end g  


 for (int ibin=0;ibin<nCBins;ibin++){
    
   for (int ibin2=0;ibin2<nPtBins;ibin2++){ 

     AllJetShapesCanvas->cd(8-ibin);

     for(int g = 0; g<6; g++){
	  
       norm_tot = JetShape[g][ibin][ibin2]->Integral(1,12,"width");

       for(int k = 1; k<13; k++){

	 
	 norm_temp = JetShape[g][ibin][ibin2]->Integral(1,k,"width");

	 err_temp = 0.;

	 if(k>1){
	   err_temp = Running_int[g][ibin][ibin2]->GetBinError(k-1);
	 }
	 cont = JetShape[g][ibin][ibin2]->GetBinContent(k);
	 if(cont > 0.){
	   err = JetShape[g][ibin][ibin2]->GetBinError(k)/cont*norm_temp/norm_tot;
	 } else { err = 0.;}
	    
	 Running_int[g][ibin][ibin2]->SetBinContent(k,norm_temp/norm_tot);
	 Running_int[g][ibin][ibin2]->SetBinError(k,TMath::Sqrt(err*err+err_temp*err_temp)/norm_tot);
	  

	 //cout<<g<<" "<<k<<" "<<norm_temp/norm_tot<<" "<<err<<" "<<err_temp<<" "<<TMath::Sqrt(err*err+err_temp*err_temp)<<endl;

	 if(g==1||g==3||g==5){

	   Running_int[g][ibin][ibin2]->SetBinError(k,0.005);
	   if(k<3){Running_int[g][ibin][ibin2]->SetBinError(k,0.03);}  //this is purely for cosmetic reasons.  delete if you will display error bars!!!
	    
	 }	    

	 //	  gPad->SetLogx();

	 Running_int[g][ibin][ibin2]->SetMarkerSize(2);
	   
	 switch(g){
	 case 0:
	   Running_int[g][ibin][ibin2]->SetMarkerStyle(10);
	   Running_int[g][ibin][ibin2]->SetMarkerColor(kBlack);
	   Running_int[g][ibin][ibin2]->SetLineColor(kBlack);
	   break;
	 case 1:
	   Running_int[g][ibin][ibin2]->SetMarkerStyle(10);
	   Running_int[g][ibin][ibin2]->SetMarkerColor(kBlack);
	   Running_int[g][ibin][ibin2]->SetLineColor(kBlack);
	   break;
	 case 2:
	   Running_int[g][ibin][ibin2]->SetMarkerStyle(34);
	   Running_int[g][ibin][ibin2]->SetMarkerColor(30);
	   Running_int[g][ibin][ibin2]->SetLineColor(30);
	   break;
	 case 3:
	   Running_int[g][ibin][ibin2]->SetMarkerStyle(3);
	   Running_int[g][ibin][ibin2]->SetMarkerSize(2);
	   Running_int[g][ibin][ibin2]->SetMarkerColor(30);
	   Running_int[g][ibin][ibin2]->SetLineColor(30);
	   break;
	 case 4:
	   Running_int[g][ibin][ibin2]->SetMarkerStyle(21);
	   Running_int[g][ibin][ibin2]->SetMarkerColor(kOrange-2);
	   Running_int[g][ibin][ibin2]->SetLineColor(kOrange-2);
	   break;
	 case 5:
	   Running_int[g][ibin][ibin2]->SetMarkerStyle(3);
	   Running_int[g][ibin][ibin2]->SetMarkerColor(kOrange-2);
	   Running_int[g][ibin][ibin2]->SetLineColor(kOrange-2);
	   break;
	 }


	 Running_int[g][ibin][ibin2]->SetBinContent(3,0.);
	 Running_int[g][ibin][ibin2]->SetBinError(3,0.);

	 Running_int[g][ibin][ibin2]->GetXaxis()->SetTitle("r");
	 Running_int[g][ibin][ibin2]->GetXaxis()->CenterTitle();

	 Running_int[g][ibin][ibin2]->GetXaxis()->SetTitleSize(0.07);
	 Running_int[g][ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
	   
	 Running_int[g][ibin][ibin2]->GetXaxis()->SetTitleOffset(0.8);

	 Running_int[g][ibin][ibin2]->SetMaximum(1.1);
	 Running_int[g][ibin][ibin2]->SetAxisRange(0.,0.99,"x");
	 Running_int[g][ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
	 Running_int[g][ibin][ibin2]->GetYaxis()->SetTitleSize(0.07);
	 Running_int[g][ibin][ibin2]->GetYaxis()->SetTitleOffset(0.9);
	 Running_int[g][ibin][ibin2]->GetYaxis()->CenterTitle();
	 Running_int[g][ibin][ibin2]->SetMinimum(0.3);
	 Running_int[g][ibin][ibin2]->GetYaxis()->SetTitle("#int_{0}^{r}( dN/d#Delta#etad#Delta#phi )dr");

	 if(ibin<3){
	      
	   Running_int[g][ibin][ibin2]->GetXaxis()->SetTitleSize(0.08);
	   Running_int[g][ibin][ibin2]->GetXaxis()->SetLabelSize(0.06);
	   Running_int[g][ibin][ibin2]->GetYaxis()->SetLabelSize(0.);
	   Running_int[g][ibin][ibin2]->GetYaxis()->SetTitleSize(0.);

	 }
	 if(g==2||g==4){   Running_int[g][ibin][ibin2]->Draw("same");	}    
	 if(g==3||g==5){  
	    
	   Running_int[g][ibin][ibin2]->Draw("same");	}    

       }
     }
   	
     if(ibin==0){

       TLegend *intlegend = new TLegend(0.4,0.4,0.95,0.7);
       //  intlegend->AddEntry(Running_int[0][ibin][ibin2],"Inclusive PbPb/pp","lpfe");
       intlegend->AddEntry(Running_int[4][ibin][ibin2],"Leading PbPb","lpfe");
       intlegend->AddEntry(Running_int[5][ibin][ibin2],"Leading pp","lpfe");
       intlegend->AddEntry(Running_int[2][ibin][ibin2],"SubLeading PbPb","lpfe");
       intlegend->AddEntry(Running_int[3][ibin][ibin2],"SubLeading pp","lpfe");
       intlegend->SetTextSize(0.055);
       intlegend->SetLineColor(kWhite);
       intlegend->Draw("same");
     }

     TLine *line = new TLine(0.,1.,1.,1.);
     line->SetLineStyle(2);
     line->Draw();


     AllJetShapesCanvas->cd(4-ibin);
	

     gPad->SetLogy();
     JetShape[0][ibin][ibin2]->Scale(100);
     JetShape[1][ibin][ibin2]->Scale(100);
     JetShape[1][ibin][ibin2]->SetMinimum(0.1);
     JetShape[1][ibin][ibin2]->SetMaximum(50000);
     JetShape[1][ibin][ibin2]->SetAxisRange(0.,0.999,"x");
     JetShape[1][ibin][ibin2]->GetYaxis()->SetLabelSize(0.06);
     JetShape[1][ibin][ibin2]->GetYaxis()->SetTitleSize(0.08);
     if(ibin<3){	JetShape[1][ibin][ibin2]->GetYaxis()->SetTitleSize(0.); 	JetShape[1][ibin][ibin2]->GetYaxis()->SetLabelSize(0.);}
     JetShape[1][ibin][ibin2]->GetYaxis()->SetTitle("dN/d#Delta#eta d#Delta#phi");
     JetShape[1][ibin][ibin2]->GetYaxis()->CenterTitle();
     JetShape[1][ibin][ibin2]->SetMarkerStyle(3);
     JetShape[1][ibin][ibin2]->Draw();
     JetShape[0][ibin][ibin2]->Draw("same");

     JetShape[3][ibin][ibin2]->SetMarkerColor(30);
     JetShape[2][ibin][ibin2]->SetMarkerColor(30);
     JetShape[3][ibin][ibin2]->SetLineColor(30);

     JetShape[2][ibin][ibin2]->SetLineColor(30);
     JetShape[3][ibin][ibin2]->SetMarkerStyle(3);
     JetShape[2][ibin][ibin2]->SetMarkerStyle(34);

     JetShape[3][ibin][ibin2]->SetAxisRange(0.,0.999,"x");
     JetShape[3][ibin][ibin2]->GetYaxis()->SetLabelSize(0.07);
     JetShape[3][ibin][ibin2]->GetYaxis()->SetTitleSize(0.08);
     JetShape[3][ibin][ibin2]->GetYaxis()->SetTitleOffset(0.9);
     if(ibin<3){	JetShape[3][ibin][ibin2]->GetYaxis()->SetTitleSize(0.); 	JetShape[3][ibin][ibin2]->GetYaxis()->SetLabelSize(0.);}
     JetShape[3][ibin][ibin2]->GetYaxis()->SetTitle("dN/d#Delta#eta d#Delta#phi");
     JetShape[3][ibin][ibin2]->GetYaxis()->CenterTitle();
     JetShape[3][ibin][ibin2]->SetMinimum(0.1);
     JetShape[3][ibin][ibin2]->SetMaximum(5000);
     JetShape[3][ibin][ibin2]->Draw("same");
     JetShape[2][ibin][ibin2]->Draw("same");

     JetShape[3][ibin][ibin2]->Draw("same");

     JetShape[5][ibin][ibin2]->Scale(10);
     JetShape[4][ibin][ibin2]->Scale(10);

     JetShape[5][ibin][ibin2]->SetMarkerColor(kOrange-2);
     JetShape[4][ibin][ibin2]->SetMarkerColor(kOrange-2);
     JetShape[5][ibin][ibin2]->SetLineColor(kOrange-2);
     JetShape[4][ibin][ibin2]->SetLineColor(kOrange-2);
	
	
     JetShape[5][ibin][ibin2]->SetMarkerStyle(3);
     JetShape[4][ibin][ibin2]->SetMarkerStyle(21);


	
     JetShape[4][ibin][ibin2]->Draw("same");
     JetShape[5][ibin][ibin2]->Draw("same");
     if(ibin ==0){
       TLegend *alllegend = new TLegend(0.4,0.65,0.99,0.98);
       alllegend->AddEntry(JetShape[0][ibin][ibin2],"Inc. PbPb x 100","lpfe");
       alllegend->AddEntry(JetShape[1][ibin][ibin2],"Inclusive pp x 100","lpfe");
       alllegend->AddEntry(JetShape[4][ibin][ibin2],"Leading PbPb x 10","lpfe");
       alllegend->AddEntry(JetShape[5][ibin][ibin2],"Leading pp x 10","lpfe");
       alllegend->AddEntry(JetShape[2][ibin][ibin2],"SubLeading PbPb","lpfe");
       alllegend->AddEntry(JetShape[3][ibin][ibin2],"SubLeading pp","lpfe");
       alllegend->SetTextSize(0.055);
       alllegend->SetLineColor(kWhite);
       alllegend->Draw();
     }


     if(ibin==3){	centtex = new TLatex(0.2,0.9,CBin_labels[ibin]);
       centtex->SetTextSize(0.06);
       centtex->SetNDC();
       centtex->Draw();
     }

     if(ibin<3){centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
       centtex->SetTextSize(0.07);
       centtex->SetNDC();
       centtex->Draw();
     }




     for(int g = 0; g<6; g++){


       Raw_Yield_TrkPt1_TrkPt8[g][ibin][ibin2]->Write();
       JetShape[g][ibin][ibin2]->Write();
       Running_int[g][ibin][ibin2]->Write();
     }
	
   } 
 }
     
 AllJetShapesCanvas->SaveAs("AllJetShapes.pdf");

 AllJetShapesDisplay->cd(1);
	
   Running_int[4][3][0]->SetMarkerStyle(25);
   Running_int[5][0][0]->SetMarkerStyle(25);
   Running_int[5][0][0]->SetLineColor(kBlack);
   Running_int[5][0][0]->SetMarkerColor(kBlack);
   Running_int[5][0][0]->GetYaxis()->SetTitleSize(0.07);
   Running_int[5][0][0]->GetYaxis()->SetLabelSize(0.05);
   Running_int[5][0][0]->GetXaxis()->SetLabelSize(0.05);
   Running_int[5][0][0]->GetYaxis()->SetTitleOffset(0.8);
   Running_int[5][0][0]->Draw();
   Running_int[4][3][0]->Draw("same");
   Running_int[4][0][0]->Draw("same");


   TLegend *intlegend1 = new TLegend(0.4,0.3,0.95,0.6);
   intlegend1->AddEntry(Running_int[5][0][0],"Leading pp","lpfe");
   intlegend1->AddEntry(Running_int[4][3][0],"Leading PbPb Cent. 50-100%","lpfe");
   intlegend1->AddEntry(Running_int[4][0][0],"Leading PbPb Cent. 0-10%","lpfe");
   intlegend1->SetTextSize(0.04);
   intlegend1->SetLineColor(kWhite);
   intlegend1->Draw("same");


   TLine *line_disp = new TLine(0.,1.,1.,1.);
   line_disp->SetLineStyle(2);
   line_disp->Draw("same");



   AllJetShapesDisplay->cd(2);
	
   Running_int[2][3][0]->SetMarkerStyle(28);
   Running_int[3][0][0]->SetMarkerStyle(28);
   Running_int[3][0][0]->SetLineColor(kBlack);
   Running_int[3][0][0]->SetMarkerColor(kBlack);
   Running_int[3][0][0]->GetXaxis()->SetLabelSize(0.05);
   Running_int[3][0][0]->Draw();
   Running_int[2][3][0]->Draw("same");
   Running_int[2][0][0]->Draw("same");


   TLegend *intlegend2 = new TLegend(0.3,0.3,0.95,0.6);
   intlegend2->AddEntry(Running_int[3][0][0],"SubLeading pp","lpfe");
   intlegend2->AddEntry(Running_int[2][3][0],"SubLeading PbPb Cent. 50-100%","lpfe");
   intlegend2->AddEntry(Running_int[2][0][0],"SubLeading PbPb Cent. 0-10%","lpfe");
   intlegend2->SetTextSize(0.04);
   intlegend2->SetLineColor(kWhite);
   intlegend2->Draw("same");

   line_disp->Draw("same");

 
   AllJetShapesDisplay->SaveAs("AllJetShapesDisplay.pdf");
   

 // AllJetShapesCanvas->SaveAs("AllJetShapes.C");
 
 // ClosureTestInclusive->SaveAs("ClosureTest.pdf");

  

     return 0;
} // main loop
