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
#include "TProfile.h"


#include <iostream>
#include <vector>
#include <fstream>

#include "../JetTrack2016_functions.h"


Int_t jet_resolution(bool is_pythia = 0){


  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);

  gStyle->SetOptTitle(0);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);


  TFile *f_in;
  if(is_pythia){
    f_in = new TFile("VertexCentJetInfo_Pythia80_Merged.root","READ");
  }else{
    f_in = new TFile("VertexCentJetInfo_Hydjet_NewOfficialJEC.root","READ");
  }
  TCanvas *jet_res_canvas;
  if(is_pythia){
    jet_res_canvas= new TCanvas("jet_res_canvas","",500,800);
    jet_res_canvas->Divide(1,2,0,0);
  }else{
    jet_res_canvas= new TCanvas("jet_res_canvas","",1500,800);
    jet_res_canvas->Divide(4,2,0,0);

  }

  const int nCBins = 4;

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

  TString corr_strs[2] = {"","JFFCorr_"};

  TH2D *jet_res[2][nCBins];
  TH1D *jet_res_proj[2][nCBins][100];
  TH1D *sigma[2][nCBins];

  TH1D *mu[2][nCBins];
  TProfile *mu2[2][nCBins];
  /*
  TF1 *gaus = new TF1("gaus","[0]*TMath::Exp(-.5*pow(((x-[1])/[2]),2))");
 
  gaus->SetParameter(0,10.);
  gaus->SetParLimits(0,0.,300.);
  gaus->SetParameter(1,1.);
  gaus->SetParLimits(1,0.5,1.5);
  gaus->SetParameter(2,0.15);
  gaus->SetParLimits(2,0.,1.);
  */


  for(int g = 0; g<1; g++){

    for(int j = 0; j<4; j++){

      if(is_pythia&&j>0){continue;}

      //  jet_res[g][j] = (TH2D*)f_in->Get((TString)(corr_strs[g]+"Jet_Resolution_Merged_"+CBin_strs[j]+"_"+CBin_strs[j+1]))->Clone((TString)(corr_strs[g]+"Jet_Resolution_Merged_"+CBin_strs[j]+"_"+CBin_strs[j+1]));

      jet_res[g][j] = (TH2D*)f_in->Get((TString)(corr_strs[g]+"Jet_Resolution_"+CBin_strs[j]+"_"+CBin_strs[j+1]))->Clone((TString)(corr_strs[g]+"Jet_Resolution_Merged_"+CBin_strs[j]+"_"+CBin_strs[j+1]));

      cout<<"here"<<endl;

      jet_res[g][j]->Rebin2D(4,1);
      
      mu2[g][j] = (TProfile*)jet_res[g][j]->ProfileX((TString)(corr_strs[g]+"Jet_Resolution_Profile_"+CBin_strs[j]+"_"+CBin_strs[j+1]),1,200);
      sigma[g][j] = new TH1D((TString)(corr_strs[g]+"Jet_Resolution_Sigma_"+CBin_strs[j]+"_"+CBin_strs[j+1]),"",jet_res[g][j]->GetNbinsX(),jet_res[g][j]->GetXaxis()->GetBinLowEdge(1),jet_res[g][j]->GetXaxis()->GetBinLowEdge(jet_res[g][j]->GetNbinsX()+1));

      sigma[g][j]->Sumw2();

      for(int k = 1; k<jet_res[g][j]->GetNbinsX(); k++){

	TString kstring = "Bin"; kstring+=k;
	jet_res_proj[g][j][k] = (TH1D*)jet_res[g][j]->ProjectionY((TString)(corr_strs[g]+"Jet_Resolution_Merged_"+kstring+"_"+CBin_strs[j]+"_"+CBin_strs[j+1]),k,k);
	float rms = jet_res_proj[g][j][k]->GetRMS();
	float mean = jet_res_proj[g][j][k]->GetMean();

	cout<<k<<" "<<sigma[g][j]->GetBinCenter(k)<<" "<<jet_res[g][j]->GetXaxis()->GetBinCenter(k)<<" "<<mean<<" "<<rms<<endl;
	
	sigma[g][j]->SetBinContent(k,rms);
	sigma[g][j]->SetBinError(k,0.001);

	//	mu2[g][j]->SetBinContent(k,mean);
	//mu2[g][j]->SetBinError(k,rms);
      }
      
    
      /*    
      mu[g][j] = (TH1D*)gDirectory->Get((TString)(corr_strs[g]+"Jet_Resolution_Merged_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_1"))->Clone((TString)(corr_strs[g]+"Jet_Resolution_Mu_"+CBin_strs[j]+"_"+CBin_strs[j+1]));
   
      sigma[g][j] = (TH1D*)gDirectory->Get((TString)(corr_strs[g]+"Jet_Resolution_Merged_"+CBin_strs[j]+"_"+CBin_strs[j+1]+"_2"))->Clone((TString)(corr_strs[g]+"Jet_Resolution_Sigma_"+CBin_strs[j]+"_"+CBin_strs[j+1]));

      */
   
      cout<<"ready to draw"<<endl;

      if(is_pythia) jet_res_canvas->cd(1);
      else jet_res_canvas->cd(4-j);

      //   sigma[g][j]->Divide(mu2[g][j]);
    
      sigma[g][j]->SetMinimum(-0.001);
      sigma[g][j]->SetMaximum(.3);
      sigma[g][j]->SetAxisRange(101.,300.);
      
      sigma[g][j]->SetLineColor(kBlack);
      sigma[g][j]->SetMarkerColor(kBlack);
      sigma[g][j]->SetMarkerStyle(20);
      sigma[g][j]->SetMarkerSize(1);
      //sigma[g][j]->GetYaxis()->SetTitle("#sigma/#mu");
      sigma[g][j]->GetYaxis()->SetTitle("#sigma(p_{T}^{Reco}/p_{T}^{Gen})");
      sigma[g][j]->GetYaxis()->SetTitleSize(0.09);
      sigma[g][j]->GetYaxis()->CenterTitle();
      sigma[g][j]->GetYaxis()->SetLabelSize(0.06);
      sigma[g][j]->GetYaxis()->SetTitleOffset(0.8);

      if(j!=3&&!is_pythia){
	sigma[g][j]->GetYaxis()->SetTitleSize(0.0);
	sigma[g][j]->GetYaxis()->SetLabelSize(0.0);
      }
      
      //gPad->SetLogx();
      
      if(g==0){
	sigma[g][j]->Draw();
      }else{
	sigma[g][j]->SetLineColor(kRed);
	sigma[g][j]->SetMarkerColor(kRed);
	sigma[g][j]->Draw("same");
      
	if(j==3||is_pythia){
	  TLegend *legend = new TLegend(0.2,0.7,0.7,0.85);
	  legend->AddEntry(	sigma[0][j],"Reco Jets (db JEC)");
	  legend->AddEntry(	sigma[1][j],"Reco with JFF-JEC");
	  legend->SetLineColor(kWhite);
	  legend->SetTextSize(0.06);
	  legend->Draw();
	 
	}
      }

      if(!is_pythia){
	TLatex *cent = new TLatex(0.65,0.9,CBin_labels[j]);
	cent->SetNDC();
	cent->SetTextSize(0.06);
	cent->Draw();
      }

      if(j==3&&!is_pythia){
	TLatex *label = new TLatex(0.2,0.9,"PYTHIA+HYDJET");
	label->SetNDC();
	label->SetTextSize(0.06);
	label->Draw();
      }else if(is_pythia){
	TLatex *label = new TLatex(0.2,0.9,"PYTHIA");
	label->SetNDC();
	label->SetTextSize(0.06);
	label->Draw();
      }
      if(j==2){
	TLatex *label = new TLatex(0.05,0.9,"akVs3Calo");
	label->SetNDC();
	label->SetTextSize(0.06);
	label->Draw();
      }
      
      if(j==1){
	TLatex *label = new TLatex(0.05,0.9,"HLT_Jet80");
	label->SetNDC();
	label->SetTextSize(0.06);
	label->Draw();
      }
      
      if(is_pythia) jet_res_canvas->cd(2);
      else jet_res_canvas->cd(8-j);
      /*
      mu[g][j]->SetMinimum(0.9);
      mu[g][j]->SetMaximum(1.1);
      mu[g][j]->SetAxisRange(50.,329.);

      
      mu[g][j]->SetLineColor(kBlack);
      mu[g][j]->SetMarkerColor(kBlack);
      mu[g][j]->SetMarkerStyle(20);
      mu[g][j]->SetMarkerSize(1);
      mu[g][j]->GetYaxis()->SetTitle("#mu(p_{T}^{Reco}/p_{T}^{Gen})");
      mu[g][j]->GetYaxis()->SetTitleSize(0.08);
      mu[g][j]->GetYaxis()->CenterTitle();
      mu[g][j]->GetYaxis()->SetLabelSize(0.05);
      mu[g][j]->GetYaxis()->SetTitleOffset(0.8);

      if(j!=3){  
	mu[g][j]->GetYaxis()->SetTitleSize(0.0); 
	mu[g][j]->GetYaxis()->SetLabelSize(0.0); 
      }
      
      mu[g][j]->GetXaxis()->CenterTitle();
      mu[g][j]->GetXaxis()->SetTitleSize(0.07);
      mu[g][j]->GetXaxis()->SetLabelSize(0.06);
      mu[g][j]->GetXaxis()->SetTitle("p_{T}^{Gen}");
      mu[g][j]->GetXaxis()->SetTitleOffset(0.8);
   
      //     mu[g][j]->Draw();

      */

      mu2[g][j]->GetYaxis()->SetTitle("#mu(p_{T}^{Reco}/p_{T}^{Gen})");
      mu2[g][j]->GetYaxis()->SetTitleSize(0.08);
      mu2[g][j]->GetYaxis()->CenterTitle();
      mu2[g][j]->GetYaxis()->SetLabelSize(0.05);
      mu2[g][j]->GetYaxis()->SetTitleOffset(0.8);

      mu2[g][j]->GetXaxis()->CenterTitle();
      mu2[g][j]->GetXaxis()->SetTitleSize(0.07);
      mu2[g][j]->GetXaxis()->SetLabelSize(0.06);
    
      mu2[g][j]->GetXaxis()->SetTitle("p_{T}^{Gen}");
      mu2[g][j]->GetXaxis()->SetTitleOffset(0.8);
      if(j==3){ mu2[g][j]->GetXaxis()->SetLabelSize(0.05); }
      if(j==3){ mu2[g][j]->GetXaxis()->SetTitleSize(0.06); }
      if(j==3){ mu2[g][j]->GetXaxis()->SetTitleOffset(0.9); }
      mu2[g][j]->SetMinimum(0.9);
      mu2[g][j]->SetMaximum(1.15);
      mu2[g][j]->SetAxisRange(101.,300.);


      if(j!=3&&!is_pythia){  
	mu2[g][j]->GetYaxis()->SetTitleSize(0.0); 
	mu2[g][j]->GetYaxis()->SetLabelSize(0.0); 
      }
      


      mu2[g][j]->SetLineColor(kBlack);
      mu2[g][j]->SetMarkerColor(kBlack);
      mu2[g][j]->SetMarkerStyle(20);
      mu2[g][j]->SetMarkerSize(1);

      //  gPad->SetLogx();
 
      if(g==0){
	mu2[g][j]->Draw();
      }else{
	mu2[g][j]->SetLineColor(kRed);
	mu2[g][j]->SetMarkerColor(kRed);
	mu2[g][j]->Draw("same");
      }

      TLine *l_mu = new TLine(90.,1.,329.,1.);
      l_mu->SetLineColor(kBlack);
      l_mu->Draw();

      TLine *l_mu_up = new TLine(90.,1.02,329.,1.02);
      l_mu_up->SetLineStyle(2);
      l_mu_up->SetLineColor(kBlack);
      l_mu_up->Draw();

      TLine *l_mu_down = new TLine(90.,.98,329.,.98);
      l_mu_down->SetLineStyle(2);
      l_mu_down->SetLineColor(kBlack);
      l_mu_down->Draw();
  
      /*
      jet_res_canvas->cd(12-j);

      //  gPad->SetLogx();
      if(g==1){
	//	jet_res[1][j]->Divide(jet_res[0][j]);
	jet_res[1][j]->Draw("colz");
      } 
      */

    }
    //  if(g==1){
      if(is_pythia){
	jet_res_canvas->SaveAs("JetResolution_Pythia.pdf");
	jet_res_canvas->SaveAs("JetResolution_Pythia.png");
      }else{
	jet_res_canvas->SaveAs("JetResolution_PythiaHydjet.pdf");
	jet_res_canvas->SaveAs("JetResolution_PythiaHydjet.png");
      }
      //  }

  }
    return 0;
}
