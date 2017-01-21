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
#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>
#include <vector>
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;

#define nCBins 4
#define nPtBins 1
#define nTrkPtBins 9

float trkPtCut=1;

int parti = -999;
bool is_data = false;


enum enum_dataset_types {e_Data_PbPb,e_Data_pp,e_HydJet, e_Pythia, e_n_dataset_types};
int dataset_type_code = -999;

TString dataset_type_strs[e_n_dataset_types] = {"Data_PbPb","Data_pp","HydJet","Pythia"};

TString dataset_type_file_names[e_n_dataset_types] = {"ClusterData_PbPb.txt","ClusterData_pp.txt","Hydjet.txt","Pythia.txt"};

int pthat_cuts[6] = {50,80,120,170,220,280};

float pthat_weight[6] = {0.88403, 0.101112, 0.010379, 0.00422609, 0.000252783, 5.04625e-05};


float PtBins[nPtBins+1] = {100, 300};
TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};

float CBins[nCBins+1] = {0, 20, 60, 100, 200};
TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

float TrkPtBins[nTrkPtBins+1] = {0.5, 1, 2, 3, 4, 8, 12, 16, 20, 300};
TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8","TrkPt12","TrkPt16","TrkPt20","TrkPt300" };


Int_t make_reweighting(int dataset_min = 2, int dataset_max = 3){

  TH1D *cent[e_n_dataset_types];
  TH1D *vertex[e_n_dataset_types];

  double wvz = 1.;
  double wcen = 1.;
 
  //////////###### PTHAT SAMPLES ###########///////////////
  TFile * wtfile_vtx_mc;

  TH1F* hWeight_vtx;
  TH1F* hWeight_MC_vtx;

  TH1F* hWeight_cent;
  TH1F* hWeight_MC_cent;

  int pthat =0;
  int pthatmax =0;
  

  TF1 *fit_vz[e_n_dataset_types];
  TF1 *fit_cen[e_n_dataset_types];

  TFile * wtfile_vtx;
  
  TFile *out_file = new TFile("VertexCentReweightingFits.root","RECREATE");

  for(int dataset_type_code = dataset_min; dataset_type_code< dataset_max; dataset_type_code++){

    fit_vz[dataset_type_code] = new TF1((TString)("Fit_Vz_"+dataset_type_strs[dataset_type_code]),"[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",-15.,15.);
    fit_cen[dataset_type_code] = new TF1((TString)("Fit_Cent_"+dataset_type_strs[dataset_type_code]),"[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.,200.);


    cout<<"about to try reweighting"<<endl;  

    TString vertex_cent_name;

    pthat = pthat_cuts[dataset_type_code];
    pthatmax=pthat_cuts[dataset_type_code+1];

    if(dataset_type_code==3){
  
      vertex_cent_name = "VertexCentJetInfo_Pythia.root"; 
      wtfile_vtx_mc = TFile::Open(vertex_cent_name,"READ");
      wtfile_vtx = TFile::Open("VertexCentJetInfo_Data_pp_Merged.root", "readonly");
    }else if(dataset_type_code == 2){
   
      vertex_cent_name = "VertexCentJetInfo_HydJet.root"; 
      wtfile_vtx_mc = TFile::Open(vertex_cent_name,"READ");
      wtfile_vtx = TFile::Open("VertexCentJetInfo_Data_PbPb_part0.root", "readonly");
    }

    hWeight_MC_vtx = (TH1F*) ((TH1F*)wtfile_vtx_mc->Get("VertexDist"))->Clone("VertexDistMCNormalized");
    hWeight_MC_vtx->Scale(1./hWeight_MC_vtx->Integral());

  
    hWeight_vtx = (TH1F*) ((TH1F*)wtfile_vtx->Get("VertexDist"))->Clone("VertexDistNormalized");

    hWeight_vtx->Scale(1./hWeight_vtx->Integral());
    hWeight_vtx->Divide(hWeight_MC_vtx);
 
    TCanvas *vz_canvas = new TCanvas("vz_canvas");
    hWeight_vtx->Draw();
    
    hWeight_vtx->Fit((TString)("Fit_Vz_"+dataset_type_strs[dataset_type_code]),"","",-15.,15.);
   
   
    TString vz_canvas_name = "TestVzFit_"; vz_canvas_name+= dataset_type_strs[dataset_type_code]; vz_canvas_name+=".png";
    vz_canvas->SaveAs(vz_canvas_name);
  
    hWeight_MC_cent = (TH1F*) ((TH1F*)wtfile_vtx_mc->Get("CentDist"))->Clone("CentDistMCNormalized");
    hWeight_MC_cent->Scale(1./hWeight_MC_cent->Integral());

   
    hWeight_cent = (TH1F*) ((TH1F*)wtfile_vtx->Get("CentDist"))->Clone("CentDistNormalized");
    hWeight_cent->Divide(hWeight_MC_cent);

    hWeight_cent->Scale(1./hWeight_cent->Integral());
    
    fit_cen[dataset_type_code]->SetParameter(1,.03);
    fit_cen[dataset_type_code]->SetParameter(2,-.03);
    hWeight_cent->Fit((TString)("Fit_Cent_"+dataset_type_strs[dataset_type_code]),"","",0.,190.);

    TString cent_canvas_name = vz_canvas_name; cent_canvas_name.ReplaceAll("Vz","Cent");
    TCanvas *cent_canvas = new TCanvas("cent_canvas");
    hWeight_cent->Draw();
    cent_canvas->SaveAs(cent_canvas_name);

  
  

    out_file->cd();
    
    fit_vz[dataset_type_code]->Write();
    if(dataset_type_code<11)fit_cen[dataset_type_code]->Write();
  }
  return 0;

}
