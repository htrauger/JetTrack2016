#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TF1.h"

#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>

#include <fstream>
#include "TMath.h"
#include <vector>



#include "HiForest_Branches_MC.h"
//#include "fragmentation_JEC_fits.h"



using namespace std;



//arg 1 = which data set, arg 2 = output file number

enum enum_dataset_types {e_Data2015,e_Data_pp,e_HydJet15,e_HydJet30,e_HydJet50, e_HydJet80, e_HydJet120,e_HydJet170,e_HydJet220,e_HydJet280, e_HydJet370,e_Pythia15,e_Pythia30,e_Pythia50, e_Pythia80, e_Pythia120,e_Pythia170,e_Pythia220,e_Pythia280, e_Pythia370, e_n_dataset_types};

TString dataset_type_strs[e_n_dataset_types] = {"Data2015","Data_pp","HydJet15","HydJet30","HydJet50","HydJet80", "HydJet120", "HydJet170","HydJet220","HydJet280","HydJet370","Pythia15","Pythia30","Pythia50","Pythia80", "Pythia120", "Pythia170","Pythia220","Pythia280","Pythia370"};

int dataset_pthats[e_n_dataset_types+1] = {0,0,   15,30,50,80,120,170,220,280,370,    15,30,50,80,120,170,220,280,370,999};


int make_ntuples_run2_MC(bool doCrab=0, int jobID=0, int endfile = 10, int dataset_type_code = 0){

  bool is_data = false;
  bool do_PbPb = 1;

  
  if(dataset_type_code == 0 || dataset_type_code == 1) is_data = true;
 
  cout<<dataset_type_code<<endl;
  if(dataset_type_code==1 ||dataset_type_code > 10){do_PbPb = 0;}
  
  if(is_data){
    cout<<"This version is only for Monte Carlo"<<endl;
    return -1;
    
  }


  std::vector<Float_t> *v_trkEta = new std::vector<Float_t>();   v_trkEta->clear();
  std::vector<Float_t> *v_trkPhi = new std::vector<Float_t>();   v_trkPhi->clear();
  std::vector<Float_t> *v_trkPt = new std::vector<Float_t>();   v_trkPt->clear();
  std::vector<Float_t> *v_trkMVALoose = new std::vector<Float_t>();   v_trkMVALoose->clear();
  std::vector<Float_t> *v_trkMVATight = new std::vector<Float_t>();   v_trkMVATight->clear();
  std::vector<Float_t> *v_highPurity = new std::vector<Float_t>();   v_highPurity->clear();
 
  std::vector<Float_t> *v_vz = new std::vector<Float_t>();   v_vz->clear();


  std::vector<Float_t> *v_calo_jteta = new std::vector<Float_t>();   v_calo_jteta->clear();
  std::vector<Float_t> *v_calo_jtphi = new std::vector<Float_t>();   v_calo_jtphi->clear();
  std::vector<Float_t> *v_calo_jtpt = new std::vector<Float_t>();   v_calo_jtpt->clear();
  std::vector<Float_t> *v_calo_rawpt = new std::vector<Float_t>();   v_calo_rawpt->clear();
  std::vector<Float_t> *v_calo_trackMax = new std::vector<Float_t>();   v_calo_trackMax->clear();  
   
  std::vector<Float_t> *v_pf_jteta = new std::vector<Float_t>();   v_pf_jteta->clear();
  std::vector<Float_t> *v_pf_jtphi = new std::vector<Float_t>();   v_pf_jtphi->clear();
  std::vector<Float_t> *v_pf_jtpt = new std::vector<Float_t>();   v_pf_jtpt->clear();
  std::vector<Float_t> *v_pf_rawpt = new std::vector<Float_t>();   v_pf_rawpt->clear();
  std::vector<Float_t> *v_pf_trackMax = new std::vector<Float_t>();   v_pf_trackMax->clear();  
   
  std::vector<Float_t> *v_trkDxy1 = new std::vector<Float_t>();   v_trkDxy1->clear();
  std::vector<Float_t> *v_trkDxyError1 = new std::vector<Float_t>();   v_trkDxyError1->clear();
  std::vector<Float_t> *v_trkDz1 = new std::vector<Float_t>();   v_trkDz1->clear();
  std::vector<Float_t> *v_trkDzError1 = new std::vector<Float_t>();   v_trkDzError1->clear();
  std::vector<Float_t> *v_trkPtError = new std::vector<Float_t>();   v_trkPtError->clear();


  std::vector<Int_t> *v_pfId = new std::vector<Int_t>(); v_pfId->clear();
  std::vector<Float_t> *v_pfPt = new std::vector<Float_t>(); v_pfPt->clear();
  std::vector<Float_t> *v_pfVsPt = new std::vector<Float_t>(); v_pfVsPt->clear();
  std::vector<Float_t> *v_pfEta = new std::vector<Float_t>(); v_pfEta->clear();
  std::vector<Float_t> *v_pfPhi = new std::vector<Float_t>(); v_pfPhi->clear();
  std::vector<Float_t> *v_sumpt = new std::vector<Float_t>(); v_sumpt->clear();
  std::vector<Float_t> *v_nPFpart = new std::vector<Float_t>(); v_nPFpart->clear();



  std::vector<Int_t> *v_sube= new std::vector<Int_t>();  v_sube->clear();
  std::vector<Float_t> *v_pt = new std::vector<Float_t>();  v_pt->clear();
  std::vector<Float_t> *v_phi = new std::vector<Float_t>();  v_phi->clear();
  std::vector<Float_t> *v_eta = new std::vector<Float_t>();  v_eta->clear();
  std::vector<Int_t> *v_chg= new std::vector<Int_t>();  v_chg->clear();
  std::vector<Int_t> *v_pdg= new std::vector<Int_t>();  v_pdg->clear();
  std::vector<Float_t> *v_pPt = new std::vector<Float_t>();  v_pPt->clear();
  std::vector<Float_t> *v_pPhi = new std::vector<Float_t>();  v_pPhi->clear();
  std::vector<Float_t> *v_pEta = new std::vector<Float_t>();  v_pEta->clear();
  std::vector<Float_t> *v_geneta = new std::vector<Float_t>();   v_geneta->clear();
  std::vector<Float_t> *v_genphi = new std::vector<Float_t>();   v_genphi->clear();
  std::vector<Float_t> *v_genpt = new std::vector<Float_t>();   v_genpt->clear();
 
  std::vector<Float_t> *v_calo_refpt = new std::vector<Float_t>();   v_calo_refpt->clear();
  std::vector<Float_t> *v_calo_refeta = new std::vector<Float_t>();   v_calo_refeta->clear();
  std::vector<Float_t> *v_calo_refphi = new std::vector<Float_t>();   v_calo_refphi->clear();
  std::vector<Float_t> *v_calo_refdrjt = new std::vector<Float_t>();   v_calo_refdrjt->clear();
  std::vector<Float_t> *v_calo_refparton_pt = new std::vector<Float_t>();   v_calo_refparton_pt->clear();
  std::vector<Float_t> *v_calo_refparton_flavor = new std::vector<Float_t>();   v_calo_refparton_flavor->clear();

  std::vector<Float_t> *v_pf_refpt = new std::vector<Float_t>();   v_pf_refpt->clear();
  std::vector<Float_t> *v_pf_refeta = new std::vector<Float_t>();   v_pf_refeta->clear();
  std::vector<Float_t> *v_pf_refphi = new std::vector<Float_t>();   v_pf_refphi->clear();
  std::vector<Float_t> *v_pf_refdrjt = new std::vector<Float_t>();   v_pf_refdrjt->clear();
  std::vector<Float_t> *v_pf_refparton_pt = new std::vector<Float_t>();   v_pf_refparton_pt->clear();
  std::vector<Float_t> *v_pf_refparton_flavor = new std::vector<Float_t>();   v_pf_refparton_flavor->clear();
  
  TTree *calo_jet_tree;
  TTree *pf_jet_tree;
  TTree *track_tree;
  TTree *hlt_tree;
  TTree *hlt_tree2;
  TTree *evt_tree;
  TTree *gen_tree;
  TTree *pf_tree;

  std::string in_file_name;
  TString output_file_base;


  TFile *output_file;
  TFile *my_file;

  if(doCrab){
    in_file_name = Form("job_input_file_list_%d.txt",jobID);
  }else{

   
  
    if(is_data&&!do_PbPb){

      in_file_name = "root://cms-xrd-global.cern.ch///store/group/phys_heavyions/kjung/pp5TeV_Jet80PD_AOD_cmssw758_24Feb2016/HighPtJet80/crab_pp5Tev_Jet80PD_24Feb16/160224_163329/0000/HiForestAOD_1.root";
    }else if(is_data&&do_PbPb){
      in_file_name = "root://cms-xrd-global.cern.ch///store/user/velicanu/HIHardProbes/HIHardProbes-HIRun2015-PromptReco-v1-FOREST-1-v22/160126_202951/0000/HiForestAOD_1.root";
    }else if(!is_data&&!do_PbPb){

      if(dataset_type_code!=14)cout<<"WATCH OUT, YOU ARE RUNNING ON PTHAT80"<<endl;

      in_file_name = "Pythia8_HiForest.txt";

    }else{
      cerr<<"need to set up to run on that sample..."<<endl;
    }
  
    cout<<"File name is "<<in_file_name<<endl;
  }
 
  std::ifstream instr(in_file_name.c_str(), std::ifstream::in);
  if(!instr.is_open()) cout << "filelist not found!! Exiting..." << endl;
  std::string filename;
  int ifile=0;
  
  while(instr>>filename && ifile<endfile){
    filename.erase(std::remove(filename.begin(), filename.end(), '"'), filename.end());
    filename.erase(std::remove(filename.begin(), filename.end(), ','), filename.end());
    filename.erase(std::remove(filename.begin(), filename.end(), '['), filename.end());
    filename.erase(std::remove(filename.begin(), filename.end(), ']'), filename.end());
    cout<<"File name is "<< filename <<endl;
    ifile++;
 
   
    int pos = filename.find_first_of('s');
    string reducedfn = filename.substr(pos);
    string xrdPrefix = "root://xrootd.unl.edu//";
    TFile *my_file = TFile::Open((xrdPrefix+reducedfn).c_str());
  
    
    if(!my_file){ cout << "File cannot be found!!" << endl; exit(1); }	

    if(my_file->IsZombie()) { 
      std::cout << "Is zombie" << std::endl;
    }    
	
    cout<<"opened file"<<endl;

    if(do_PbPb){
      hlt_tree2 = (TTree*)my_file->Get("hltanalysis/HltTree");
      hlt_tree2->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1", &HLT_HIPuAK4CaloJet80_Eta5p1_v1, &b_HLT_HIPuAK4CaloJet80_Eta5p1_v1);
    }
  
  
    if(do_PbPb){
      calo_jet_tree = (TTree*)  my_file->Get("akVs4CaloJetAnalyzer/t");
      pf_jet_tree = (TTree*)  my_file->Get("akVs4PFJetAnalyzer/t");
    }else{
      calo_jet_tree = (TTree*)  my_file->Get("ak4CaloJetAnalyzer/t");
      pf_jet_tree = (TTree*)  my_file->Get("ak4PFJetAnalyzer/t");
    }

  
    calo_jet_tree->SetBranchAddress("nref", &nref, &b_nref);
    calo_jet_tree->SetBranchAddress("trackMax", trackMax, &b_trackMax);
    calo_jet_tree->SetBranchAddress("rawpt", rawpt, &b_rawpt);
    calo_jet_tree->SetBranchAddress("jtpt", jtpt, &b_jtpt);
    calo_jet_tree->SetBranchAddress("jteta", jteta, &b_jteta);
    calo_jet_tree->SetBranchAddress("jtphi", jtphi, &b_jtphi);
 
    pf_jet_tree->SetBranchAddress("nref", &nref, &b_nref);
    pf_jet_tree->SetBranchAddress("trackMax", trackMax, &b_trackMax);
    pf_jet_tree->SetBranchAddress("rawpt", rawpt, &b_rawpt);
    pf_jet_tree->SetBranchAddress("jtpt", jtpt, &b_jtpt);
    pf_jet_tree->SetBranchAddress("jteta", jteta, &b_jteta);
    pf_jet_tree->SetBranchAddress("jtphi", jtphi, &b_jtphi);
 
    if(!is_data){

      calo_jet_tree->SetBranchAddress("ngen",&ngen, &b_ngen);
      calo_jet_tree->SetBranchAddress("pthat",&pthat, &b_pthat);
      calo_jet_tree->SetBranchAddress("genpt", genpt, &b_genpt);
      calo_jet_tree->SetBranchAddress("geneta", geneta, &b_geneta);
      calo_jet_tree->SetBranchAddress("genphi", genphi, &b_genphi);


      calo_jet_tree->SetBranchAddress("refpt", &refpt, &b_refpt);
      calo_jet_tree->SetBranchAddress("refeta", &refeta, &b_refeta);
      calo_jet_tree->SetBranchAddress("refphi",&refphi, &b_refphi);
      calo_jet_tree->SetBranchAddress("refdrjt", &refdrjt, &b_refdrjt);
      calo_jet_tree->SetBranchAddress("refparton_pt", &refparton_pt, &b_refparton_pt);
      calo_jet_tree->SetBranchAddress("refparton_flavor", &refparton_flavor, &b_refparton_flavor);

      

      pf_jet_tree->SetBranchAddress("refpt", &refpt, &b_refpt);
      pf_jet_tree->SetBranchAddress("refeta", &refeta, &b_refeta);
      pf_jet_tree->SetBranchAddress("refphi",&refphi, &b_refphi);
      pf_jet_tree->SetBranchAddress("refdrjt", &refdrjt, &b_refdrjt);
      pf_jet_tree->SetBranchAddress("refparton_pt", &refparton_pt, &b_refparton_pt);
      pf_jet_tree->SetBranchAddress("refparton_flavor", &refparton_flavor, &b_refparton_flavor);



      gen_tree = (TTree*) my_file->Get("HiGenParticleAna/hi");

      gen_tree->SetBranchAddress("mult",&mult, &b_mult);
      gen_tree->SetBranchAddress("pt",&pt, &b_pt);
      gen_tree->SetBranchAddress("eta",&eta, &b_eta);
      gen_tree->SetBranchAddress("phi",&phi, &b_phi);
      gen_tree->SetBranchAddress("pdg",&pdg, &b_pdg);
      gen_tree->SetBranchAddress("chg",&chg, &b_chg);
      gen_tree->SetBranchAddress("sube",&sube, &b_sube);
    }
 

    evt_tree = (TTree*) my_file->Get("hiEvtAnalyzer/HiTree");

    evt_tree->SetBranchAddress("evt", &evt, &b_evt);
    evt_tree->SetBranchAddress("vx", &vx, &b_vx);
    evt_tree->SetBranchAddress("vy", &vy, &b_vy);
    evt_tree->SetBranchAddress("vz", &vz, &b_vz);
    evt_tree->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
 
    hlt_tree = (TTree*) my_file->Get("skimanalysis/HltTree");



    hlt_tree->SetBranchAddress("ana_step", &ana_step, &b_ana_step);
    hlt_tree->SetBranchAddress("pHBHENoiseFilterResultProducer", &pHBHENoiseFilterResultProducer, &b_pHBHENoiseFilterResultProducer);
    hlt_tree->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult, &b_HBHENoiseFilterResult);
    hlt_tree->SetBranchAddress("HBHENoiseFilterResultRun1", &HBHENoiseFilterResultRun1, &b_HBHENoiseFilterResultRun1);
    hlt_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose, &b_HBHENoiseFilterResultRun2Loose);
    hlt_tree->SetBranchAddress("HBHENoiseFilterResultRun2Tight",&HBHENoiseFilterResultRun2Tight, &b_HBHENoiseFilterResultRun2Tight);
    hlt_tree->SetBranchAddress("HBHEIsoNoiseFilterResult", &HBHEIsoNoiseFilterResult, &b_HBHEIsoNoiseFilterResult);


    if(do_PbPb){

      hlt_tree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection,&b_pcollisionEventSelection);
      hlt_tree->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter,&b_pprimaryVertexFilter);
      hlt_tree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter,&b_pclusterCompatibilityFilter);
      hlt_tree->SetBranchAddress("superFilterPath",&superFilterPath,&b_superFilterPath);
      hlt_tree->SetBranchAddress("phfCoincFilter1",&phfCoincFilter1,&b_phfCoincFilter1);
      hlt_tree->SetBranchAddress("phfCoincFilter2",&phfCoincFilter2,&b_phfCoincFilter2);
      hlt_tree->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3,&b_phfCoincFilter3);
      hlt_tree->SetBranchAddress("phfCoincFilter4",&phfCoincFilter4,&b_phfCoincFilter4);
      hlt_tree->SetBranchAddress("phfCoincFilter5",&phfCoincFilter5,&b_phfCoincFilter5);
 

    }else{
      hlt_tree->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter,&b_pPAprimaryVertexFilter);
      hlt_tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter,&b_pBeamScrapingFilter);
      hlt_tree->SetBranchAddress("pVertexFilterCutG",&pVertexFilterCutG, &b_pVertexFilterCutG);
      hlt_tree->SetBranchAddress("pVertexFilterCutGloose",&pVertexFilterCutGloose, &b_pVertexFilterCutGloose);
      hlt_tree->SetBranchAddress("pVertexFilterCutGtight",&pVertexFilterCutGtight,&b_pVertexFilterCutGtight);
      hlt_tree->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplus, &b_pVertexFilterCutGplus);
      hlt_tree->SetBranchAddress("pVertexFilterCutE",&pVertexFilterCutE, &b_pVertexFilterCutE);
      hlt_tree->SetBranchAddress("pVertexFilterCutEandG", &pVertexFilterCutEandG, &b_pVertexFilterCutEandG);

    }
 
    if(do_PbPb){
      track_tree = (TTree*) my_file->Get("anaTrack/trackTree");
    }else{
      track_tree = (TTree*) my_file->Get("ppTrack/trackTree");
    }
 
    track_tree->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
    track_tree->SetBranchAddress("trkPt", trkPt, &b_trkPt);
    track_tree->SetBranchAddress("trkPtError", trkPtError, &b_trkPtError);
    track_tree->SetBranchAddress("trkEta", trkEta, &b_trkEta);
    track_tree->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
    track_tree->SetBranchAddress("pfHcal", pfHcal, &b_pfHcal);
    track_tree->SetBranchAddress("pfEcal", pfEcal, &b_pfEcal);
    track_tree->SetBranchAddress("highPurity", highPurity, &b_highPurity);
    track_tree->SetBranchAddress("trkMVALoose", trkMVALoose, &b_trkMVALoose);
    track_tree->SetBranchAddress("trkMVATight", trkMVATight, &b_trkMVATight);
    track_tree->SetBranchAddress("highPurity", highPurity, &b_highPurity);
   

    pf_tree = (TTree*)my_file->Get("pfcandAnalyzer/pfTree");
  

    pf_tree = (TTree*)my_file->Get("pfcandAnalyzer/pfTree");
  
    pf_tree->SetBranchAddress("nPFpart", &nPFpart, &b_nPFpart);
    pf_tree->SetBranchAddress("pfId", &pfId, &b_pfId);
    pf_tree->SetBranchAddress("pfPt", &pfPt, &b_pfPt);
    if(do_PbPb)  pf_tree->SetBranchAddress("pfPuPt", &pfVsPt, &b_pfVsPt);
    pf_tree->SetBranchAddress("pfEta", &pfEta, &b_pfEta);
    pf_tree->SetBranchAddress("pfPhi", &pfPhi, &b_pfPhi);
    pf_tree->SetBranchAddress("sumpt", &sumpt, &b_sumpt);

    int n_evt = evt_tree->GetEntriesFast();
  
  
    /*
    if(!is_data && do_PbPb){
      output_file_base= "/afs/cern.ch/work/h/htrauger/public/JetShapes2016Skims/Hydjet/";
    }else if (is_data&&do_PbPb){
      output_file_base= "skim_output/PbPb_Demo/";
    }else if (is_data&&!do_PbPb){
      output_file_base= "/afs/cern.ch/work/h/htrauger/public/JetShapes2016Skims/pp/";
    }else if (!is_data&&!do_PbPb){
      output_file_base= "/afs/cern.ch/work/h/htrauger/public/JetShapes2016Skims/Pythia/";
   
    }else{
      cerr<<"nope, we can't handle that data set"<<endl;
      return -1;
    }
    */
    if(doCrab){
    output_file_base = "./";
    }else{
      output_file_base= "/afs/cern.ch/work/h/htrauger/public/JetShapes2016Skims/";
    }
    output_file_base = dataset_type_strs[dataset_type_code];

    output_file = new TFile((TString) (output_file_base + ".root"), "RECREATE");


    output_file->cd();

    TTree *mixing_tree = new TTree("mixing_tree", "");

    mixing_tree->Branch("nTrk", &nTrk, "nTrk/I");
    mixing_tree->Branch("trkEta", "vector<Float_t>", &v_trkEta);
    mixing_tree->Branch("trkPhi", "vector<Float_t>", &v_trkPhi);
    mixing_tree->Branch("trkPt", "vector<Float_t>", &v_trkPt);
    mixing_tree->Branch("highPurity", "vector<Float_t>", &v_highPurity);
    mixing_tree->Branch("trkMVALoose", "vector<Float_t>", &v_trkMVALoose);
    mixing_tree->Branch("trkMVATight", "vector<Float_t>", &v_trkMVATight);
    mixing_tree->Branch("vz", "vector<Float_t>", &v_vz);
 

    mixing_tree->Branch("ana_step", &ana_step, "anastep/I");
    mixing_tree->Branch("pHBHENoiseFilterResultProducer", &pHBHENoiseFilterResultProducer, "pHBHENoiseFilterResultProducer/I");
    mixing_tree->Branch("HBHENoiseFilterResult",&HBHENoiseFilterResult, "HBHENoiseFilterResult/I");
    mixing_tree->Branch("HBHENoiseFilterResultRun1", &HBHENoiseFilterResultRun1, "HBHENoiseFilterResultRun1/I");
    mixing_tree->Branch("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose, "HBHENoiseFilterResultRun2Loose/I");
    mixing_tree->Branch("HBHENoiseFilterResultRun2Tight",&HBHENoiseFilterResultRun2Tight, "HBHENoiseFilterResultRun2Tight/I");
    mixing_tree->Branch("HBHEIsoNoiseFilterResult", &HBHEIsoNoiseFilterResult, "HBHEIsoNoiseFilterResult/I");


    if(do_PbPb){

      mixing_tree->Branch("pcollisionEventSelection",&pcollisionEventSelection,"pcollisionEventSelection/I");
      mixing_tree->Branch("pprimaryVertexFilter",&pprimaryVertexFilter,"pprimaryVertexFilter/I");
      mixing_tree->Branch("pclusterCompatibilityFilter",&pclusterCompatibilityFilter,"pclusterCompatibilityFilter/I");
      mixing_tree->Branch("superFilterPath",&superFilterPath,"superFilterPath/I");
      mixing_tree->Branch("phfCoincFilter1",&phfCoincFilter1,"phfCoincFilter1/I");
      mixing_tree->Branch("phfCoincFilter2",&phfCoincFilter2,"phfCoincFilter2/I");
      mixing_tree->Branch("phfCoincFilter3",&phfCoincFilter3,"phfCoincFilter3/I");
      mixing_tree->Branch("phfCoincFilter4",&phfCoincFilter4,"phfCoincFilter4/I");
      mixing_tree->Branch("phfCoincFilter5",&phfCoincFilter5,"phfCoincFilter5/I");

    }else{
      mixing_tree->Branch("pPAprimaryVertexFilter",&pPAprimaryVertexFilter, "pPAprimaryVertexFilter/I");
      mixing_tree->Branch("pBeamScrapingFilter",&pBeamScrapingFilter, "pBeamScrapingFilter/I");
      mixing_tree->Branch("pVertexFilterCutG",&pVertexFilterCutG, "pVertexFilterCutG/I");
      mixing_tree->Branch("pVertexFilterCutGloose",&pVertexFilterCutGloose, "pVertexFilterCutGloose/I");
      mixing_tree->Branch("pVertexFilterCutGtight",&pVertexFilterCutGtight, "pVertexFilterCutGtight/I");
      mixing_tree->Branch("pVertexFilterCutGplus", &pVertexFilterCutGplus, "pVertexFilterCutGplus/I");
      mixing_tree->Branch("pVertexFilterCutE",&pVertexFilterCutE, "pVertexFilterCutE/I");
      mixing_tree->Branch("pVertexFilterCutEandG", &pVertexFilterCutEandG, "VertexFilterCutEandG/I");
    }
    mixing_tree->Branch("hiBin", &hiBin, "hiBin/I");
    
    mixing_tree->Branch("calo_jteta", "vector<Float_t>", &v_calo_jteta);
    mixing_tree->Branch("calo_jtphi", "vector<Float_t>", &v_calo_jtphi);
    mixing_tree->Branch("calo_jtpt", "vector<Float_t>", &v_calo_jtpt);
    mixing_tree->Branch("calo_rawpt", "vector<Float_t>", &v_calo_rawpt);
    mixing_tree->Branch("calo_trackMax", "vector<Float_t>", &v_calo_trackMax);
    
    mixing_tree->Branch("pf_jteta", "vector<Float_t>", &v_pf_jteta);
    mixing_tree->Branch("pf_jtphi", "vector<Float_t>", &v_pf_jtphi);
    mixing_tree->Branch("pf_jtpt", "vector<Float_t>", &v_pf_jtpt);
    mixing_tree->Branch("pf_rawpt", "vector<Float_t>", &v_pf_rawpt);
    mixing_tree->Branch("pf_trackMax", "vector<Float_t>", &v_pf_trackMax);
    
    mixing_tree->Branch("pthat", &pthat, "pthat/F");
    mixing_tree->Branch("trkDxy1", "vector<Float_t>", &v_trkDxy1);
    mixing_tree->Branch("trkDxyError1", "vector<Float_t>", &v_trkDxyError1);
    mixing_tree->Branch("trkDz1", "vector<Float_t>", &v_trkDz1);
    mixing_tree->Branch("trkDzError1", "vector<Float_t>", &v_trkDzError1);
    mixing_tree->Branch("trkPtError", "vector<Float_t>", &v_trkPtError);
  
    mixing_tree->Branch("nPFpart", &nPFpart, "nPFpart/I");
    mixing_tree->Branch("pfId", "vector<Int_t>", &v_pfId);
    mixing_tree->Branch("pfPt", "vector<Float_t>", &v_pfPt);
    mixing_tree->Branch("pfVsPt", "vector<Float_t>", &v_pfVsPt);
    mixing_tree->Branch("pfEta", "vector<Float_t>", &v_pfEta);
    mixing_tree->Branch("pfPhi", "vector<Float_t>", &v_pfPhi);
  
    if(!is_data){
      mixing_tree->Branch("mult", &mult, "mult/I");
      mixing_tree->Branch("pt", "vector<Float_t>", &v_pt);
      mixing_tree->Branch("phi", "vector<Float_t>", &v_phi);
      mixing_tree->Branch("eta", "vector<Float_t>", &v_eta);
      mixing_tree->Branch("chg", "vector<Int_t>", &v_chg);
      mixing_tree->Branch("pdg", "vector<Int_t>", &v_pdg);
      mixing_tree->Branch("sube", "vector<Int_t>", &v_sube);
      mixing_tree->Branch("pPt", "vector<Float_t>", &v_pPt);
      mixing_tree->Branch("pPhi", "vector<Float_t>", &v_pPhi);
      mixing_tree->Branch("pEta", "vector<Float_t>", &v_pEta);

      mixing_tree->Branch("geneta", "vector<Float_t>", &v_geneta);
      mixing_tree->Branch("genphi", "vector<Float_t>", &v_genphi);
      mixing_tree->Branch("genpt", "vector<Float_t>", &v_genpt);

      mixing_tree->Branch("calo_refpt", "vector<Float_t>", &v_calo_refpt);
      mixing_tree->Branch("calo_refeta", "vector<Float_t>", &v_calo_refeta);
      mixing_tree->Branch("calo_refphi", "vector<Float_t>", &v_calo_refphi);
      mixing_tree->Branch("calo_refdrjt", "vector<Float_t>", &v_calo_refdrjt);
      mixing_tree->Branch("calo_refparton_pt", "vector<Float_t>", &v_calo_refparton_pt);
      mixing_tree->Branch("calo_refparton_flavor", "vector<Float_t>", &v_calo_refparton_flavor);

      mixing_tree->Branch("pf_refpt", "vector<Float_t>", &v_pf_refpt);
      mixing_tree->Branch("pf_refeta", "vector<Float_t>", &v_pf_refeta);
      mixing_tree->Branch("pf_refphi", "vector<Float_t>", &v_pf_refphi);
      mixing_tree->Branch("pf_refdrjt", "vector<Float_t>", &v_pf_refdrjt);
      mixing_tree->Branch("pf_refparton_pt", "vector<Float_t>", &v_pf_refparton_pt);
      mixing_tree->Branch("pf_refparton_flavor", "vector<Float_t>", &v_pf_refparton_flavor);



    }

    int ev_min = 0;
    Int_t ev_max = n_evt+1;
     
    if(ev_max > 1000000) ev_max = 100000;
    

    cout << "ev_min: " << ev_min << ", Entries: " << n_evt << std::endl;

    std::cout << "Will run from event number " << ev_min << " to " << ev_max -1 << "\n";
    for(int evi = ev_min; evi < ev_max; evi++) {
      
      if( evi % 1000 == 0 )  std::cout << "evi: " << evi <<  " of " << n_evt << "\n";
      //if( evi > 1000 ) break;


      if(do_PbPb){
	hlt_tree2->GetEntry(evi);
	if(HLT_HIPuAK4CaloJet80_Eta5p1_v1==0){
	  continue;
	}     
      }

   
   
      evt_tree->GetEntry(evi);
   
      hlt_tree->GetEntry(evi);
  
      track_tree->GetEntry(evi);
  
     
      calo_jet_tree->GetEntry(evi);

      for(int j4i = 0; j4i < nref ; j4i++) {

	if( fabs(jteta[j4i]) > 2. ) continue;

	if(jtpt[j4i] < 25) continue;

	v_calo_jteta->push_back(jteta[j4i]);
	v_calo_jtphi->push_back(jtphi[j4i]);
	v_calo_jtpt->push_back(jtpt[j4i]);
	v_calo_rawpt->push_back(rawpt[j4i]);
	v_calo_trackMax->push_back(trackMax[j4i]);

	if(!is_data){
	  v_calo_refpt->push_back(refpt[j4i]);
	  v_calo_refeta->push_back(refeta[j4i]);
	  v_calo_refphi->push_back(refphi[j4i]);
	  v_calo_refdrjt->push_back(refdrjt[j4i]);
	  v_calo_refparton_pt->push_back(refparton_pt[j4i]);
	  v_calo_refparton_flavor->push_back(refparton_flavor[j4i]);
	  
	}

      } /// jet loop
     
      if(!is_data){

	for(int j4i_gen = 0; j4i_gen < ngen ; j4i_gen++) {

	  if( fabs(geneta[j4i_gen]) > 2 ) continue;
	  if( genpt[j4i_gen] < 30 ) continue;
	
	  v_geneta->push_back(geneta[j4i_gen]);
	  v_genphi->push_back(genphi[j4i_gen]);
	  v_genpt->push_back(genpt[j4i_gen]);

	} /// genjet loop

      }
      
      pf_jet_tree->GetEntry(evi);

      for(int j4i = 0; j4i < nref ; j4i++) {

	if( fabs(jteta[j4i]) > 2. ) continue;
  
	if(jtpt[j4i] < 25) continue;

	v_pf_jteta->push_back(jteta[j4i]);
	v_pf_jtphi->push_back(jtphi[j4i]);
	v_pf_jtpt->push_back(jtpt[j4i]);
	v_pf_rawpt->push_back(rawpt[j4i]);
	v_pf_trackMax->push_back(trackMax[j4i]);


	if(!is_data){
	  v_pf_refpt->push_back(refpt[j4i]);
	  v_pf_refeta->push_back(refeta[j4i]);
	  v_pf_refphi->push_back(refphi[j4i]);
	  v_pf_refdrjt->push_back(refdrjt[j4i]);
	  v_pf_refparton_pt->push_back(refparton_pt[j4i]);
	  v_pf_refparton_flavor->push_back(refparton_flavor[j4i]);
	  
	}

      } /// jet loop
      
      // cout<<"ran jet loop"<<endl;
   
      pf_tree->GetEntry(evi);

      if(is_data){
	for(int pfi = 0; pfi< nPFpart ; pfi++) {
    
	  v_pfId->push_back(pfId->at(pfi));
	  v_pfPt->push_back(pfPt->at(pfi));
	  if(do_PbPb) v_pfVsPt->push_back(pfVsPt->at(pfi));
	  v_pfEta->push_back(pfEta->at(pfi));
	  v_pfPhi->push_back(pfPhi->at(pfi));
	  v_sumpt->push_back(sumpt);
     
	} /// particle flow candidate loop
      }
  
      //// reco track loop
      for(int itrk=0;itrk<nTrk;itrk++){

	//very basic cuts

	if(highPurity[itrk]!=1) continue;
     
	if(trkPtError[itrk]/trkPt[itrk]>=0.3 || TMath::Abs(trkDz1[itrk]/trkDzError1[itrk])>=3.0 ||TMath::Abs(trkDxy1[itrk]/trkDxyError1[itrk])>=3.0) continue ;

	float Et = (pfHcal[itrk]+pfEcal[itrk])/TMath::CosH(trkEta[itrk]);
	if(!(trkPt[itrk]<20 || (Et>0.2*trkPt[itrk] && Et>trkPt[itrk]-80))) continue;


	float eta=trkEta[itrk];
	if(fabs(eta)>=2.4) continue; //acceptance of the tracker   

	float pt=trkPt[itrk];
	if(pt <= 0.5) continue; //pt min
	if(pt >= 300 ) continue; //pt max

	// reco track quantities

	
	v_trkEta->push_back(trkEta[itrk]);
	v_trkPhi->push_back(trkPhi[itrk]);
	v_trkPt->push_back(trkPt[itrk]);
	v_highPurity->push_back(highPurity[itrk]);
	v_trkMVALoose->push_back(trkMVALoose[itrk]);
	v_trkMVATight->push_back(trkMVATight[itrk]);

	
      }    

      //  cout<<"ran track loop"<<endl;

      if(!is_data){


	gen_tree->GetEvent(evi);
	//// gen track loop
	for(int itrk=0;itrk<mult;itrk++){

	
	  v_trkEta->push_back(trkEta[itrk]);
	  v_trkPhi->push_back(trkPhi[itrk]);
	  v_trkPt->push_back(trkPt[itrk]);
	  v_highPurity->push_back(highPurity[itrk]);
	  v_trkMVALoose->push_back(trkMVALoose[itrk]);
	  v_trkMVATight->push_back(trkMVATight[itrk]);

	  v_pt->push_back(pt->at(itrk));
	  v_eta->push_back(eta->at(itrk));
	  v_phi->push_back(phi->at(itrk));
	  v_pdg->push_back(pdg->at(itrk));
	  v_chg->push_back(chg->at(itrk));
	  v_sube->push_back(sube->at(itrk));
	  
  	}   
      } // !is_data
      //  cout<<"ran track loop"<<endl;


      v_vz->push_back(vz);
   

      ///// Fill it


      mixing_tree->Fill();

   	 
      v_trkEta->clear();
      v_trkPhi->clear();
      v_trkPt->clear();
      v_trkMVALoose->clear();
      v_trkMVATight->clear();
    
      v_highPurity->clear();
  
      v_vz->clear();
  
      v_calo_jteta->clear();
      v_calo_jtphi->clear();
      v_calo_jtpt->clear();
      v_calo_rawpt->clear();
      v_calo_trackMax->clear();

      v_calo_refpt->clear();
      v_calo_refeta->clear();
      v_calo_refphi->clear();
      v_calo_refdrjt->clear();
      v_calo_refparton_pt->clear();
      v_calo_refparton_flavor->clear();
      
      v_pf_jteta->clear();
      v_pf_jtphi->clear();
      v_pf_jtpt->clear();
      v_pf_rawpt->clear();
      v_pf_trackMax->clear();
      

      v_pf_refpt->clear();
      v_pf_refeta->clear();
      v_pf_refphi->clear();
      v_pf_refdrjt->clear();
      v_pf_refparton_pt->clear();
      v_pf_refparton_flavor->clear();

      v_trkDxy1->clear();
      v_trkDxyError1->clear();
      v_trkDz1->clear();
      v_trkDzError1->clear();
      v_trkPtError->clear();

      v_pfId->clear();
      v_pfPt->clear();
      v_pfVsPt->clear();
      v_pfEta->clear();
      v_pfPhi->clear();
      v_sumpt->clear();
      v_nPFpart->clear();

      if(!is_data){
	v_pt->clear();
	v_phi->clear();
	v_eta->clear();
	v_chg->clear();
	v_pdg->clear();
	v_sube->clear();

	v_pPt->clear();
	v_pPhi->clear();
	v_pEta->clear();


	v_geneta->clear();
	v_genphi->clear();
	v_genpt->clear();
      }

      if( evi % 1000 == 0 ) std::cout << "Filled successfully" << std::endl;

    }  ///event loop
   
  
    cout<<"writing"<<endl;

    output_file->Write();

    output_file->Close();

  }
  


  return 0;
  
}




