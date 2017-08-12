
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
#include "TLegend.h"

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TExec.h"
#include "TLatex.h"

#include <iostream>
#include <vector>
#include <fstream>

#include "../JetTrack2016_functions.h"

using namespace std;

Int_t supplementary_figures(){

  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.15);
  gStyle->SetPadLeftMargin  (0.16);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);
  gStyle->SetTextFont(43);
  gStyle->SetCanvasBorderMode(0);
 
 

  TFile *fin[12];
  TFile *fin_ref[12];
  TFile *fin2[12];


  Double_t xAxis[5] = {-100,-50,-30,-10,0}; 
  TH1D* int_cent[12][4];
  TH1D* blank[4];
  TH1D* blank2[4];

  TH1D *check_new_phi_rebin[12][5][4];
  TH1D *check_new_phi_ref[12][5][4];
  TH1D *check_new_phi_syst[12][5][4];
  
  TH1D *check_new_eta_rebin[12][5][4];
  TH1D *check_new_eta_ref[12][5][4];
  TH1D *check_new_eta_syst[12][5][4];

  TH1D *PbPb_pp_phi[12][5][4];
  TH1D *PbPb_pp_phi_syst[12][5][4];

  TH1D *PbPb_pp_eta[12][5][4];
  TH1D *PbPb_pp_eta_syst[12][5][4];

  TH1D *PbPb_pp_phi_ref[12][5][4];
  TH1D *PbPb_pp_eta_ref[12][5][4];

  TH1D *HYDJET_PYTHIA_eta[12][5][4];
  TH1D *HYDJET_PYTHIA_phi[12][5][4];

 
  TF1 *gaus_phi[12][5][4];
  TF1 *gaus_eta[12][5][4];

  
  TH1D *Integral_phi_syst[12][4];
  TH1D *Integral_phi_Pt[12][4];

  TH1D *dummy_pTaxis[4];
 
  TH1D *Integral_eta_Pt[12][4];
  TH1D *Integral_eta_ref_Pt[12][4];
  TGraphAsymmErrors *Integral_eta_syst[12][4];
  TGraphAsymmErrors *Integral_eta_cent[12][4];

  TGraph *Error_up_eta_pT[12][4];
  TGraph *Error_down_eta_pT[12][4];
  TH1D *Error_up_eta_cent[12][4];
  TH1D *Error_down_eta_cent[12][4];

  TLegend *l40,*l41,*l42;

  TString datalabel, in_name, centlabel, pTlabel;
  TString   leadingrefetaname[5][4];

  TCanvas *cFigure1,*cFigure2, *cFigure3, *cFigure4;

  float raw_min,raw_max,mixed_min,mixed_max,yield_min,yield_max,result_min,result_max,  val_l,val_r,err_l,err_r;
 
  TLine *linePhi,*lineEta,*linePt, *lineCent;

  TLatex *tex24phi,*tex24eta;


  float value, error, evalpt, value2, dx_eta, dx_phi;

float tspixels = 25;

  vector<float> pTbin_centers;
  pTbin_centers.push_back(1.5);
  pTbin_centers.push_back(2.5);
  pTbin_centers.push_back(3.5);
  pTbin_centers.push_back(6.0);
  vector<float> pTbin_errors;
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(2.);

  float  dummy_pTbins[] = {0.8,1.0,2.0,3.0,4.0,8.0,8.2 };


  //*************************************
  //   Now draw PAS plots
  //***********************************

 
    TString figure1_name = "ME_Correction";
    cFigure1 = new TCanvas(figure1_name," ",10,10,1500,500);
    cFigure1->Divide(3,1,0.001,0.001);
  
    TString figure2_name = "Background_Subtraction";
    cFigure2 = new TCanvas(figure2_name," ",10,10,1500,500);
    cFigure2->Divide(3,1,0.001,0.001);

    TString figure3_name = "Results_Presentation";
    cFigure3 = new TCanvas(figure3_name," ",10,10,1500,500);
    cFigure3->Divide(3,1,0.001,0.001);

    TString figure4_name = "Results_2D";
    cFigure4 = new TCanvas(figure4_name," ",10,10,500,500);
     
    //--------------------
    //   PAS FIGURE 1
    //--------------------------
    cFigure1->cd(1);

    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);


    TFile *fbkgsummed = TFile::Open("../me_correct/PbPb_Inclusive_Correlations.root");
       
    TString temp_name = ("Raw_Yield_Cent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2");

    cout<<temp_name<<endl;

    TH2D *raw_corr = (TH2D*)fbkgsummed->Get(temp_name)->Clone(temp_name);
    raw_min = 10.;
    raw_max = 60.;
    yield_min = 45.;
    yield_max = 55.;
    mixed_min = 0.;
    mixed_max = 1.4;

    result_min = -.5;
    result_max = 9.5;
    
    raw_corr->Rebin2D(10,4);
       
    dx_eta = raw_corr->GetXaxis()->GetBinWidth(1);
    dx_phi = raw_corr->GetYaxis()->GetBinWidth(1);
     
    raw_corr->Scale(1/dx_eta/dx_phi);
    raw_corr->SetMinimum(raw_min);
    raw_corr->SetMaximum(raw_max);
    raw_corr->SetLineColor(kBlack);
    raw_corr->GetXaxis()->SetRangeUser(-3.,3.);
    raw_corr->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    raw_corr->GetXaxis()->SetTitle("#Delta#eta");
    raw_corr->GetYaxis()->SetTitle("#Delta#phi");
    raw_corr->GetXaxis()->CenterTitle();
    raw_corr->GetYaxis()->CenterTitle();
    raw_corr->GetZaxis()->CenterTitle();
    raw_corr->GetYaxis()->SetTitleSize(0.06);
    raw_corr->GetXaxis()->SetTitleSize(0.06);
    raw_corr->GetXaxis()->SetTitleOffset(1.);
    raw_corr->GetYaxis()->SetTitleOffset(1.);
    raw_corr->GetZaxis()->SetTitle("S(#Delta#eta,#Delta#phi)_{Raw}");
    raw_corr->GetZaxis()->SetTitleSize(0.06);
    raw_corr->GetZaxis()->SetTitleOffset(1.);
    raw_corr->GetZaxis()->SetNdivisions(205);
     
      
    raw_corr->Draw("surf1");

    //  drawlabels_PAS_bkg(i);

   

    cFigure1->cd(2);

    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);

    temp_name.ReplaceAll("Raw_Yield","Mixed_Event");

    TH2D *mixed_bkg = (TH2D*) fbkgsummed->Get(temp_name)->Clone(temp_name);
      
    mixed_bkg->Rebin2D(10,4);
    mixed_bkg->Scale(1/40.);


    mixed_bkg->GetXaxis()->SetRangeUser(-3.,3.);
    mixed_bkg->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.0001);
    mixed_bkg->SetLineColor(kBlack);
    mixed_bkg->GetXaxis()->SetTitle("#Delta#eta");
    mixed_bkg->GetYaxis()->SetTitle("#Delta#phi");
    mixed_bkg->GetXaxis()->CenterTitle();
    mixed_bkg->GetYaxis()->CenterTitle();
    mixed_bkg->GetZaxis()->CenterTitle();
    mixed_bkg->GetYaxis()->SetTitleSize(0.06);
    mixed_bkg->GetXaxis()->SetTitleSize(0.06);
    mixed_bkg->GetZaxis()->SetTitle("ME(#Delta#eta,#Delta#phi)");
    mixed_bkg->GetXaxis()->SetTitleOffset(1.);
    mixed_bkg->GetYaxis()->SetTitleOffset(1.);
    mixed_bkg->GetZaxis()->SetTitleOffset(1.);
    mixed_bkg->SetMinimum(mixed_min);
    mixed_bkg->SetMaximum(mixed_max);
    
    mixed_bkg->GetXaxis()->SetTitle("#Delta#eta");
    mixed_bkg->GetYaxis()->SetTitle("#Delta#phi");
    mixed_bkg->GetZaxis()->SetTitle("ME(#Delta#eta,#Delta#phi)");
    mixed_bkg->GetZaxis()->SetTitleSize(0.06);
    mixed_bkg->GetZaxis()->SetNdivisions(205);
     
    
    mixed_bkg->Draw("surf1");

    //   drawlabels_PAS_bkg(i);

    // leading_tex->Draw();

    cFigure1->cd(3);
 
    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);

    temp_name.ReplaceAll("Mixed_Event","Yield_PbPb");
    TH2D *raw_yield = (TH2D*) fbkgsummed->Get(temp_name)->Clone(temp_name);

    raw_yield->Rebin2D(10,4);
    raw_yield->Scale(1/dx_eta/dx_phi);
      
    raw_yield->GetXaxis()->SetRangeUser(-3.,3.);
    raw_yield->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    raw_yield->GetXaxis()->SetTitle("#Delta#eta");
    raw_yield->GetYaxis()->SetTitle("#Delta#phi");
    raw_yield->GetXaxis()->CenterTitle();
    raw_yield->GetYaxis()->CenterTitle();
    raw_yield->GetZaxis()->CenterTitle();
    raw_yield->GetYaxis()->SetTitleSize(0.06);
    raw_yield->GetXaxis()->SetTitleSize(0.06);
    raw_yield->GetZaxis()->SetTitle("S(#Delta#eta,#Delta#phi)");
    raw_yield->GetZaxis()->SetTitleSize(0.06);
    raw_yield->GetXaxis()->SetTitleOffset(1.);
    raw_yield->GetYaxis()->SetTitleOffset(1.);
    raw_yield->GetZaxis()->SetTitleOffset(1.);
     
    raw_yield->SetMinimum(yield_min);
    raw_yield->SetMaximum(yield_max);
    raw_yield->GetZaxis()->SetNdivisions(205);
    raw_yield->SetLineColor(kBlack);
    raw_yield->Draw("surf1");
      

    TLatex *tex00, *tex01, *tex02,*tex03,*tex04,*tex05, *tex11;



    gStyle->SetTextFont(43);

  

    cFigure1->cd(0);

    TLatex *cms_tex = new TLatex(0.04,0.9,"CMS");
      cms_tex->SetTextFont(63);
      cms_tex->SetTextSizePixels(tspixels);
      cms_tex->SetLineColor(kWhite);
      cms_tex->SetNDC();
      cms_tex->Draw(); 


      TLatex *prelim_tex = new TLatex(0.08,0.9,"Preliminary");
      prelim_tex->SetTextFont(53);
      prelim_tex->SetTextSizePixels(tspixels);
      prelim_tex->SetLineColor(kWhite);
      prelim_tex->SetNDC();
      prelim_tex->Draw(); 
  
      TLatex   *luminosity_tex_pp = new TLatex(0.2,0.9,"pp 27.4 pb^{-1} (5.02 TeV)");
      luminosity_tex_pp->SetTextFont(43);
      luminosity_tex_pp->SetTextSizePixels(25);
      luminosity_tex_pp->SetLineColor(kWhite);
      luminosity_tex_pp->SetNDC();
      //  luminosity_tex_pp->Draw();
 
      TLatex   *luminosity_tex_PbPb = new TLatex(0.2,0.9,"PbPb 404 #mub^{-1} (5.02 TeV)");
      luminosity_tex_PbPb->SetTextFont(43);
      luminosity_tex_PbPb->SetTextSizePixels(tspixels);
      luminosity_tex_PbPb->SetLineColor(kWhite);
      luminosity_tex_PbPb->SetNDC();
      luminosity_tex_PbPb->Draw();
 
      TLatex   *jet_reco_tex = new TLatex(0.4,0.9,"anti-k_{T} calorimeter jets, R=0.4, p_{T}> 120 GeV, |#eta_{jet}| < 1.6");
      jet_reco_tex->SetTextFont(43);
      jet_reco_tex->SetTextSizePixels(tspixels);
      jet_reco_tex->SetLineColor(kWhite);
      jet_reco_tex->SetNDC();
      jet_reco_tex->Draw();

      TLatex *track_tex = new TLatex(0.8,0.9,"1 < p_{T}^{trk} < 2 GeV");
      track_tex->SetName("tex04");
      track_tex->SetNDC();
      track_tex->SetTextFont(43);
      track_tex->SetTextSizePixels(tspixels);
      track_tex->Draw();




      cout<<"done fig 1"<<endl;
    
      //---------------------------
      //   PAS FIGURE 2
      //--------------------------



      cFigure2->cd(1);

      gPad->SetTheta(60.839);
      gPad->SetPhi(38.0172);

      raw_yield->Draw("surf1");
    //   drawlabels_PAS_bkg(i);


      

    cFigure2->cd(2);
  gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);


    temp_name.ReplaceAll("Yield_PbPb","SummedBkg");

    TH2D *background = (TH2D*)fbkgsummed->Get(temp_name)->Clone(temp_name);
     
    background->Rebin2D(10,4);
    background->Scale(1/dx_eta/dx_phi);
      
    background->GetXaxis()->SetRangeUser(-3.,3.);
    background->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    background->GetXaxis()->SetTitle("#Delta#eta");
    background->GetYaxis()->SetTitle("#Delta#phi");
    background->GetZaxis()->SetTitle("B(#Delta#eta,#Delta#phi)");
    background->GetXaxis()->CenterTitle();
    background->GetYaxis()->CenterTitle();
    background->GetZaxis()->CenterTitle();
    background->GetXaxis()->SetLabelSize(0.045);
    background->GetXaxis()->SetTitleSize(0.06);
    background->GetXaxis()->SetTitleOffset(1.);
    background->GetYaxis()->SetLabelSize(0.045);
    background->GetYaxis()->SetTitleSize(0.06);
    background->GetYaxis()->SetTitleOffset(1.);
    background->GetZaxis()->SetLabelSize(0.04);
    background->GetZaxis()->SetTitleSize(0.06);
    background->GetZaxis()->SetTitleOffset(1.);
    background->GetZaxis()->SetNdivisions(205);
    background->SetMinimum(yield_min);
    background->SetMaximum(yield_max);
    background->SetLineColor(kBlack);
    background->Draw("surf1");
      


    cFigure2->cd(3);

    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);

    temp_name.ReplaceAll("SummedBkg","Yield_BkgSub");

    TH2D *result = (TH2D*)fbkgsummed->Get(temp_name)->Clone(temp_name);
     
    result->Rebin2D(10,4);
    result->Scale(1/dx_eta/dx_phi);
      
    result->GetXaxis()->SetRangeUser(-3.,3.);
    result->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    result->GetXaxis()->SetTitle("#Delta#eta");
    result->GetYaxis()->SetTitle("#Delta#phi");
    result->GetZaxis()->SetTitle("S(#Delta#eta,#Delta#phi) - B(#Delta#eta,#Delta#phi)");
    result->GetXaxis()->CenterTitle();
    result->GetYaxis()->CenterTitle();
    result->GetZaxis()->CenterTitle();
    result->GetXaxis()->SetLabelSize(0.045);
    result->GetXaxis()->SetTitleSize(0.06);
    result->GetXaxis()->SetTitleOffset(1.);
    result->GetYaxis()->SetLabelSize(0.045);
    result->GetYaxis()->SetTitleSize(0.06);
    result->GetYaxis()->SetTitleOffset(1.);
    result->GetZaxis()->SetLabelSize(0.04);
    result->GetZaxis()->SetTitleSize(0.06);
    result->GetZaxis()->SetTitleOffset(1.);
    result->GetZaxis()->SetNdivisions(205);
    result->SetMinimum(result_min);
    result->SetMaximum(result_max);
    result->SetLineColor(kBlack);
    result->Draw("surf1");
  gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);

      
    //   drawlabels_PAS_bkg(i);

    //   track_tex->Draw();
 
    cFigure2->cd(0);

    cms_tex->Draw();
    prelim_tex->Draw();
    //  luminosity_tex_pp->Draw();
    luminosity_tex_PbPb->Draw();
    jet_reco_tex->Draw();


    track_tex->Draw();

    cout<<"done figure 2"<<endl;
    //Other plots go by centrality class



    TFile *f_eta_phi = TFile::Open("../particle_yields/Particle_Yields.root ");

    TH1D *eta_proj = (TH1D*)f_eta_phi->Get("Proj_dEta_PbPb_Cent0_Cent10_TrkPt1_TrkPt2_Rebin")->Clone("eta_proj");

    TH1D *eta_syst = (TH1D*)f_eta_phi->Get("dEta_Syst_PbPb_Cent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2")->Clone("eta_syst");


    TH1D *phi_proj = (TH1D*)f_eta_phi->Get("Proj_dPhi_PbPb_Cent0_Cent10_TrkPt1_TrkPt2_Rebin")->Clone("phi_proj");
    TH1D *phi_syst = (TH1D*)f_eta_phi->Get("dPhi_Syst_PbPb_Cent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2")->Clone("phi_syst");

    TFile *f_dr = TFile::Open("../jet_shapes_result/Jet_Shapes.root");
    cout<<"here"<<endl;

    TH1D *dr_proj = (TH1D*)f_dr->Get("JetShape2_Yield_BkgSub_Inclusive_Cent0_Cent10_TrkPt1_TrkPt2")->Clone("dr_proj");
    cout<<"and here"<<endl;
    TH1D *dr_syst = (TH1D*)f_dr->Get("Jet_Shape_SystErr_Cent0_Cent10_TrkPt1_TrkPt2")->Clone("dr_syst");

    cout<<"got syst"<<endl;

    cFigure3->cd(1);

    eta_proj->GetYaxis()->SetTitleSize(0.06);
    eta_proj->GetYaxis()->SetLabelSize(0.06);
    eta_proj->GetXaxis()->SetLabelSize(0.06);
    eta_proj->GetXaxis()->SetTitleSize(0.06);
    eta_proj->GetYaxis()->SetTitle("Y = #frac{1}{N_{jets}} #frac{dN}{d#Delta#eta}");
    eta_proj->GetXaxis()->SetTitle("#Delta#eta");
    eta_proj->GetXaxis()->SetRangeUser(-2.,2.);
    eta_proj->GetXaxis()->CenterTitle();
    eta_proj->SetMaximum(10.);
    eta_proj->SetMinimum(-.5);
    eta_proj->Draw();
    eta_syst->SetFillColor(kCyan-6);
    eta_syst->Draw("e2 same");
    eta_proj->Draw("same");

    TLine *line = new TLine(-1.5,0.,1.5,0.);
    line->SetLineStyle(2);
    line->Draw("same");

    TLegend *legend = new TLegend(0.2,0.75,0.8,0.8);
    legend->SetLineColor(kWhite);
    legend->AddEntry(eta_syst,"PbPb 1 < p_{T}^{trk} < 2 GeV");
    legend->Draw();


    cFigure3->cd(2);

    phi_proj->GetYaxis()->SetTitleSize(0.06);
    phi_proj->GetYaxis()->SetLabelSize(0.06);
    phi_proj->GetXaxis()->SetLabelSize(0.06);
    phi_proj->GetXaxis()->SetTitleSize(0.06);
    phi_proj->GetXaxis()->SetTitle("#Delta#phi");
    phi_proj->GetYaxis()->SetTitle("Y = #frac{1}{N_{jets}} #frac{dN}{d#Delta#phi}");
    phi_proj->GetXaxis()->CenterTitle();
    phi_proj->SetMaximum(10.);
    phi_proj->SetMinimum(-.5);
    phi_proj->Draw();
    phi_syst->SetFillColor(kCyan-6);
    phi_syst->Draw("e2 same");
    phi_proj->Draw("same");
    line->Draw("same");

    cFigure3->cd(3);

  dr_proj->GetYaxis()->SetTitleSize(0.06);
    dr_proj->GetYaxis()->SetLabelSize(0.06);
    dr_proj->GetXaxis()->SetLabelSize(0.06);
    dr_proj->GetXaxis()->SetTitleSize(0.06);
    dr_proj->GetXaxis()->SetTitle("#Delta r");
    dr_proj->GetYaxis()->SetTitle("Y = #frac{1}{N_{jets}} #frac{dN}{d#Deltar}");
    dr_proj->GetXaxis()->CenterTitle();
    dr_proj->SetMaximum(10.);
    dr_proj->SetMinimum(-.5);
    dr_proj->GetXaxis()->SetRangeUser(0.,1.);
    dr_proj->GetXaxis()->SetNdivisions(505);
    dr_proj->Draw(); 
    dr_proj->SetMarkerColor(kBlack);
    dr_proj->SetLineColor(kBlack);
    dr_proj->SetMarkerStyle(20);
    dr_proj->SetMarkerSize(1);
    dr_syst->SetFillColor(kCyan-6);
    dr_syst->Draw("e2 same");
    dr_proj->Draw("same");

    TLine *line_r = new TLine(0.,0.,1.,0.);
    line_r->SetLineStyle(2);
    line_r->Draw("same");
 

    cFigure3->cd(0);

    cms_tex->Draw();
    prelim_tex->Draw();
    // luminosity_tex_pp->Draw();
    luminosity_tex_PbPb->Draw();
    jet_reco_tex->Draw();
  


    int ndivisions;
    float diff_min, diff_max;
  
    //--------------------------------------------
    //   Save all PAS plots except for Figure 7
    //-----------------------------------------
  
   
    cFigure1->Update();
    cFigure2->Update();
   
    figure1_name+=".pdf";
    cFigure1->SaveAs(figure1_name);
    figure1_name.ReplaceAll("pdf","png");
    cFigure1->SaveAs(figure1_name);     
  
    figure2_name+=".pdf";
    cFigure2->SaveAs(figure2_name);
    figure2_name.ReplaceAll("pdf","png");
    cFigure2->SaveAs(figure2_name);   
     

    figure3_name+=".pdf";
    cFigure3->SaveAs(figure3_name);
    figure3_name.ReplaceAll("pdf","png");
    cFigure3->SaveAs(figure3_name);   




    cFigure4->cd();

    result->GetXaxis()->SetRangeUser(-1.5,1.5);
    result->GetYaxis()->SetRangeUser(-1.5,1.5);
    result->GetYaxis()->SetTitleOffset(1.5);
    result->SetMinimum(-1.);
    result->Draw("surf1");

    cms_tex = new TLatex(0.04,0.9,"CMS");
      cms_tex->SetTextFont(63);
      cms_tex->SetTextSizePixels(18);
      cms_tex->SetLineColor(kWhite);
      cms_tex->SetNDC();
      cms_tex->Draw(); 


      prelim_tex = new TLatex(0.14,0.9,"Preliminary");
      prelim_tex->SetTextFont(53);
      prelim_tex->SetTextSizePixels(18);
      prelim_tex->SetLineColor(kWhite);
      prelim_tex->SetNDC();
      prelim_tex->Draw(); 
 
      luminosity_tex_pp = new TLatex(0.4,0.9,"pp 27.4 pb^{-1}");
      luminosity_tex_pp->SetTextFont(43);
      luminosity_tex_pp->SetTextSizePixels(18);
      luminosity_tex_pp->SetLineColor(kWhite);
      luminosity_tex_pp->SetNDC();
      //   luminosity_tex_pp->Draw();
 
      luminosity_tex_PbPb = new TLatex(0.5,0.9,"PbPb 404 #mub^{-1} (5.02 TeV)");
      luminosity_tex_PbPb->SetTextFont(43);
      luminosity_tex_PbPb->SetTextSizePixels(18);
      luminosity_tex_PbPb->SetLineColor(kWhite);
      luminosity_tex_PbPb->SetNDC();
      luminosity_tex_PbPb->Draw();
 
      jet_reco_tex = new TLatex(0.04,0.86,"anti-k_{T} calorimeter jets, R=0.4, p_{T}> 120 GeV, |#eta_{jet}| < 1.6");
      jet_reco_tex->SetTextFont(43);
      jet_reco_tex->SetTextSizePixels(18);
      jet_reco_tex->SetLineColor(kWhite);
      jet_reco_tex->SetNDC();
      jet_reco_tex->Draw();

      track_tex = new TLatex(0.04,0.8,"1 < p_{T}^{trk} < 2 GeV");
      track_tex->SetTextFont(43);
      track_tex->SetTextSizePixels(18);
      track_tex->SetLineColor(kWhite);
      track_tex->SetNDC();
      track_tex->Draw();



  figure4_name+=".pdf";
    cFigure4->SaveAs(figure4_name);
    figure4_name.ReplaceAll("pdf","png");
    cFigure4->SaveAs(figure4_name);   
      
    return 0;
  
} //Close main loop

