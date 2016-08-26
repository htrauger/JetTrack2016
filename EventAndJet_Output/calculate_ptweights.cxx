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


Int_t calculate_ptweights(int min_sample = 3, int max_sample = 19, bool is_pythia = 0){


enum enum_dataset_types {e_Data2011,e_Data_pp,e_HydJet15,e_HydJet30,e_HydJet50, e_HydJet80, e_HydJet120,e_HydJet170,e_HydJet220,e_HydJet280, e_HydJet370,e_Pythia15,e_Pythia30,e_Pythia50, e_Pythia80, e_Pythia120,e_Pythia170, e_n_dataset_types};
int dataset_type_code = -999;

 double weight[e_n_dataset_types];
 int total_events[e_n_dataset_types+1] = {0, 0, 0, 0,0,0,0, 0 ,  0 ,  0 , 0};
 float cross_section[e_n_dataset_types+1] = {0., 0.,5.335E-01 , 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06 ,6.160E-07, 1.088E-07};


 int dataset_pthats[e_n_dataset_types+1] = {0,0,15,30,50,80,120,170,220,280,370,500};

 int this_i= 0;
 TFile *f_unweighted;
 if(is_pythia) f_unweighted = new TFile("VertexCentJetInfo_Pythia.root");
 else  f_unweighted = new TFile("VertexCentJetInfo_HydJet_NoWeights.root");

 int min_bin, max_bin;

 TH1D *h_unweighted = (TH1D*)f_unweighted->Get("PthatDist")->Clone("PthatDist");

 TH1D *h_weighted = (TH1D*)f_unweighted->Get("PthatDist")->Clone("PthatDistWeighted");

 for(int sample_i = min_sample; sample_i<max_sample+1; sample_i++){

   min_bin = h_unweighted->GetXaxis()->FindBin(dataset_pthats[sample_i]+0.001);
   max_bin = h_unweighted->GetXaxis()->FindBin(dataset_pthats[sample_i+1]-0.001);

   if(sample_i == max_sample) max_bin = h_unweighted->GetNbinsX();

   cout<<min_bin<<" "<<max_bin<<endl;

   total_events[sample_i] = h_unweighted->Integral(min_bin,max_bin);

   cout<<sample_i<<" "<<dataset_pthats[sample_i]<<" to "<<dataset_pthats[sample_i+1]<<": "<<min_bin<<" "<<max_bin<<" "<<   total_events[sample_i] <<endl;

 } 


  double summed_cross_section = cross_section[min_sample]-cross_section[min_sample+1];
  double summed_total_events = total_events[min_sample];
  double summed_weight = (cross_section[min_sample]-cross_section[min_sample+1])/total_events[min_sample];



  for(int sample_i = min_sample+1; sample_i<max_sample+1; sample_i++){
    summed_cross_section += (cross_section[sample_i]-cross_section[sample_i+1]);
    summed_total_events += total_events[sample_i];
    summed_weight +=  (cross_section[sample_i]-cross_section[sample_i+1])/total_events[sample_i];

  }

    for(int sample_i = min_sample; sample_i<max_sample+1; sample_i++){
     
      weight[sample_i] = (cross_section[sample_i]-cross_section[sample_i+1])/total_events[sample_i]/summed_weight;
      
      cout<<sample_i<<" "<<dataset_pthats[sample_i]<<" "<<   weight[sample_i]<<endl;
      
    }

    for(int sample_i = min_sample; sample_i<max_sample+1; sample_i++){
      cout<< weight[sample_i]<<", ";
    }
    
    

    for(int k = 1; k< h_unweighted->GetNbinsX()+1; k++){
            
      for(int i = min_sample; i<max_sample + 1; i++){


	if(h_unweighted->GetBinLowEdge(k) >= dataset_pthats[i] && h_unweighted->GetBinLowEdge(k+1) < dataset_pthats[i+1] ){

	  this_i = i;
	}
      }
  
      float bc = h_unweighted->GetBinContent(k)*weight[this_i];
     
      h_weighted->SetBinContent(k,bc);
    
    }


    TCanvas *c = new TCanvas("my_canvas");
    
    h_weighted->Draw();

    h_unweighted->SetLineColor(kRed);
    h_unweighted->Draw("same");

    gPad->SetLogy();
  
    c->SaveAs("PthatCheckDist.png");
  
  return 0;

}
