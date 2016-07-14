#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TF1.h"
#include <string>
#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>

#include <fstream>
#include "TMath.h"

#include <utility>
#include <vector>


Int_t HLT_HIPuAK4CaloJet80_Eta5p1_v1;
TBranch *b_HLT_HIPuAK4CaloJet80_Eta5p1_v1;



//JETS// 



ULong64_t           evt;
Float_t         b;
Int_t           nref;
Int_t           ngen;
Float_t         pthat;  

Float_t         genpt[206];   //[nref]
Float_t         geneta[206];   //[nref]
Float_t         genphi[206];   //[nref]

Float_t         trackMax[206];   //[nref] 
Float_t         rawpt[206];   //[nref]
Float_t         jtpt[206];   //[nref]
Float_t         jteta[206];   //[nref]
Float_t         jty[206];   //[nref]
Float_t         jtphi[206];   //[nref]


TBranch        *b_evt;   //!
TBranch        *b_b;   //!
TBranch        *b_nref;   //!
TBranch        *b_ngen;   //!
TBranch        *b_pthat;   //!
TBranch        *b_genpt;   //!
TBranch        *b_geneta;   //!
TBranch        *b_genphi;   //!

TBranch        *b_trackMax;   //!
TBranch        *b_rawpt;   //!
TBranch        *b_jtpt;   //!
TBranch        *b_jteta;   //!
TBranch        *b_jty;   //!
TBranch        *b_jtphi;   //!


///TRACKS//


Int_t           nTrk;
Float_t         trkPt[12739];   //[nTrk]
Float_t         trkPtError[12739];   //[nTrk]
Int_t           trkNHit[12739];   //[nTrk]
Int_t           trkNlayer[12739];   //[nTrk]
Float_t         trkEta[12739];   //[nTrk]
Float_t         trkPhi[12739];   //[nTrk]
Float_t         pfEcal[12739];   //[nTrk]
Float_t         pfHcal[12739];   //[nTrk]

Bool_t         trkMVALoose[12739];   //[nTrk]
Bool_t         trkMVATight[12739];   //[nTrk]
Bool_t          highPurity[12739];   //[nTrk]
Float_t         trkDxy1[12739];   //[nTrk]
Float_t         trkDxyError1[12739];   //[nTrk]
Float_t         trkDz1[12739];   //[nTrk]
Float_t         trkDzError1[12739];   //[nTrk]
 


TBranch        *b_nTrk;   //!
TBranch        *b_trkPt;   //!
TBranch        *b_trkPtError;   //!
TBranch        *b_trkNHit;   //!
TBranch        *b_trkNlayer;   //!
TBranch        *b_trkEta;   //!
TBranch        *b_trkPhi;   //!
TBranch        *b_pfHcal;   //!
TBranch        *b_pfEcal;   //!
TBranch       *b_trkMVALoose;
TBranch         *b_trkMVATight;
TBranch *b_highPurity;
TBranch        *b_trkDxy1;   //!
TBranch        *b_trkDxyError1;   //!
TBranch        *b_trkDz1;   //!
TBranch        *b_trkDzError1;   //!


//GEN PARTICLES//
Int_t           mult;


   std::vector<Float_t>  *pt;   
   std::vector<Float_t>  *eta;   
   std::vector<Float_t>  *phi;   
   std::vector<Int_t>  *pdg;   
   std::vector<Int_t> *chg;   
   std::vector<Int_t>  *sube;

TBranch        *b_mult;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_chg;   //!
   TBranch        *b_sube;   //!

//EVENT//

UInt_t         lumi;
Float_t         vx;
Float_t         vy;
Float_t         vz;
Int_t           hiBin;

TBranch        *b_lumi;   //!
TBranch        *b_vx;   //!
TBranch        *b_vy;   //!
TBranch        *b_vz;   //!
TBranch        *b_hiBin;   //!


//FILTERS//

Int_t ana_step;

Int_t pHBHENoiseFilterResultProducer;
 Int_t HBHENoiseFilterResult;
 Int_t  HBHENoiseFilterResultRun1;
 Int_t  HBHENoiseFilterResultRun2Loose;
 Int_t  HBHENoiseFilterResultRun2Tight;
 Int_t  HBHEIsoNoiseFilterResult;
  Int_t pPAprimaryVertexFilter;
  Int_t pBeamScrapingFilter;
  Int_t pVertexFilterCutG;
  Int_t pVertexFilterCutGloose;
  Int_t pVertexFilterCutGtight;
  Int_t pVertexFilterCutGplus;
  Int_t pVertexFilterCutE;
  Int_t pVertexFilterCutEandG;


Int_t pcollisionEventSelection;
Int_t pprimaryVertexFilter;
Int_t pclusterCompatibilityFilter;
Int_t superFilterPath; 
Int_t phfCoincFilter1;
Int_t phfCoincFilter2;
Int_t phfCoincFilter3;
Int_t phfCoincFilter4;
Int_t phfCoincFilter5;


TBranch *b_ana_step;
TBranch *b_pHBHENoiseFilterResultProducer;
 TBranch *b_HBHENoiseFilterResult;
 TBranch *b_HBHENoiseFilterResultRun1;
 TBranch *b_HBHENoiseFilterResultRun2Loose;
 TBranch *b_HBHENoiseFilterResultRun2Tight;
 TBranch *b_HBHEIsoNoiseFilterResult;
  TBranch *b_pPAprimaryVertexFilter;
  TBranch *b_pBeamScrapingFilter;
  TBranch *b_pVertexFilterCutG;
  TBranch *b_pVertexFilterCutGloose;
  TBranch *b_pVertexFilterCutGtight;
  TBranch *b_pVertexFilterCutGplus;
  TBranch *b_pVertexFilterCutE;
  TBranch *b_pVertexFilterCutEandG;


TBranch *b_pcollisionEventSelection;
TBranch *b_pprimaryVertexFilter;
TBranch *b_pclusterCompatibilityFilter;
TBranch *b_superFilterPath; 
TBranch *b_phfCoincFilter1;
TBranch *b_phfCoincFilter2;
TBranch *b_phfCoincFilter3;
TBranch *b_phfCoincFilter4;
TBranch *b_phfCoincFilter5;



//PFCAND//

  Int_t         nPFpart;   //[nref]

std::vector<Int_t> *pfId;
std::vector<Float_t> *pfPt;
std::vector<Float_t> *pfVsPt;
std::vector<Float_t> *pfEta;
std::vector<Float_t> *pfPhi;
Float_t sumpt;


   // List of branches
   TBranch        *b_nPFpart;   //!
   TBranch        *b_pfId;   //!
   TBranch        *b_pfPt;   //!
   TBranch        *b_pfVsPt;   //!
   TBranch        *b_pfEta;   //!
   TBranch        *b_pfPhi;   //!
   TBranch        *b_sumpt;   //!

