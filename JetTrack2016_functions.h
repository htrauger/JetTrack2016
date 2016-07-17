
#include "JetTrack2016_universals.h"


TH1D *Rebin_dPhi(TH1D* hold){

  Double_t dx_new;

  Double_t bin_bounds_phi[]   =  {-1.50796, -1.00531,-0.879646, -.75398, -0.628319,-0.502655, -0.376991, -0.251327, -0.125664, 0.125664, 0.251327, 0.376991, 0.502655, 0.628319,.75398, 0.879646, 1.00531,1.50796};

  //  Double_t bin_bounds_phi[]=  {-1.507964,-1.382300,-1.25664,-1.130973,-1.005310,-0.879646,-0.753982,-0.62832,-0.502656,-0.376991,-0.251327,-0.125664,0.125664,0.251327,0.376991,0.502656,0.62832,0.753982,0.879646,1.005310,1.130973,1.25664,1.382300,1.507964};
  Int_t new_bins = sizeof(bin_bounds_phi) / sizeof(bin_bounds_phi[0])-1; // Should be 35 or so now
  
  // important cross check. These need to be binboundaries of the old one!
  for (int i=0; i<new_bins+1 ; ++i ){
    bool binok=false;
    for (int j=1; j<=hold->GetNbinsX()+1 ; ++j ){
      if ( fabs( bin_bounds_phi[i] - hold->GetBinLowEdge(j) )<=1e-3 ) {
	binok=true;
	break;
      }
    }
    if (!binok ){
      cerr << "Bin edge " << bin_bounds_phi[i] << " is NOT a bin boundary" << endl;
      return 0;
    }
  }
    
  TString new_name = TString( hold->GetName() ) + "_Rebin" ;
  TH1D* hnew = (TH1D*) hold->Rebin ( new_bins, new_name, bin_bounds_phi );
  
  // Now renormalize per dPhi
  Double_t dx = hold->GetXaxis()->GetBinWidth(1);
  for (int j=1; j<=hnew->GetNbinsX()+1 ; ++j ){
    dx_new = hnew->GetXaxis()->GetBinWidth(j);
    hnew->SetBinContent( j, hnew->GetBinContent( j ) * dx / dx_new );
    hnew->SetBinError(   j, hnew->GetBinError( j ) * dx / dx_new );
  }
      

  return hnew;
}



TH1D *Rebin_dEta(TH1D* hold){
  // New bin boundaries
  Double_t bin_bounds_eta[]
    //  =   {-2.5,-2.,-1.5,  -1,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1, 1.5,2.,2.5};
    =   {-2.5,-2.,-1.5,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.5,2.,2.5};
 
  Int_t new_bins = sizeof(bin_bounds_eta) / sizeof(bin_bounds_eta[0])-1; // Should be 35 or so now
  
  Double_t dx_new;
  // important cross check. These need to be binboundaries of the old one!
  for (int i=0; i<new_bins+1 ; ++i ){
    bool binok=false;
    for (int j=1; j<=hold->GetNbinsX()+1 ; ++j ){
      if ( fabs( bin_bounds_eta[i] - hold->GetBinLowEdge(j) )<1e-4 ) {
	binok=true;
	break;
      }
    }
    if (!binok ){
      cerr << "Bin edge " << bin_bounds_eta[i] << " is NOT a bin boundary" << endl;
      return 0;
    }
  }

  TString new_name = TString( hold->GetName() ) + "_Rebin" ;
  TH1D* hnew = (TH1D*) hold->Rebin ( new_bins, new_name, bin_bounds_eta );
  
  // Now renormalize per dPhi
  Double_t dx = hold->GetXaxis()->GetBinWidth(1);

  for (int j=1; j<=hnew->GetNbinsX()+1 ; ++j ){
    dx_new = hnew->GetXaxis()->GetBinWidth(j);
    hnew->SetBinContent( j, hnew->GetBinContent( j ) * dx / dx_new );
    hnew->SetBinError(   j, hnew->GetBinError( j ) * dx / dx_new);
  }
      

  return hnew;
}




TH1D *Rebin_dEta2(TH1D* hold){
  // New bin boundaries
  Double_t bin_bounds_eta[]
    //  =   {-2.5,-2.,-1.5,  -1,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1, 1.5,2.,2.5};
    =   {-2.,-1.5,  -1,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1, 1.5,2.};
 
  Int_t new_bins = sizeof(bin_bounds_eta) / sizeof(bin_bounds_eta[0])-1; // Should be 35 or so now
  

  // important cross check. These need to be binboundaries of the old one!
  for (int i=0; i<new_bins+1 ; ++i ){
    bool binok=false;
    for (int j=1; j<=hold->GetNbinsX()+1 ; ++j ){
      if ( fabs( bin_bounds_eta[i] - hold->GetBinLowEdge(j) )<1e-4 ) {
	binok=true;
	break;
      }
    }
    if (!binok ){
      cerr << "Bin edge " << bin_bounds_eta[i] << " is NOT a bin boundary" << endl;
      return 0;
    }
  }

  TString new_name = TString( hold->GetName() ) + "_Rebin" ;
  TH1D* hnew = (TH1D*) hold->Rebin ( new_bins, new_name, bin_bounds_eta );
  
  // Now renormalize per dPhi
  Double_t dx = hold->GetXaxis()->GetBinWidth(1);
  for (int j=1; j<=hnew->GetNbinsX()+1 ; ++j ){
    hnew->SetBinContent( j, hnew->GetBinContent( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
    hnew->SetBinError(   j, hnew->GetBinError( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
  }
      

  return hnew;
}






TH1D *Rebin_dEta3(TH1D* hold){
  // New bin boundaries
  Double_t bin_bounds_eta[]
    =   {-3.,-2.5,-2.,-1.5,  -1,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1, 1.5,2.,2.5,3.};
 
  Int_t new_bins = sizeof(bin_bounds_eta) / sizeof(bin_bounds_eta[0])-1; // Should be 35 or so now
  

  // important cross check. These need to be binboundaries of the old one!
  for (int i=0; i<new_bins+1 ; ++i ){
    bool binok=false;
    for (int j=1; j<=hold->GetNbinsX()+1 ; ++j ){
      if ( fabs( bin_bounds_eta[i] - hold->GetBinLowEdge(j) )<1e-4 ) {
	binok=true;
	break;
      }
    }
    if (!binok ){
      cerr << "Bin edge " << bin_bounds_eta[i] << " is NOT a bin boundary" << endl;
      return 0;
    }
  }

  TString new_name = TString( hold->GetName() ) + "_Rebin" ;
  TH1D* hnew = (TH1D*) hold->Rebin ( new_bins, new_name, bin_bounds_eta );
  
  // Now renormalize per dPhi
  Double_t dx = hold->GetXaxis()->GetBinWidth(1);
  for (int j=1; j<=hnew->GetNbinsX()+1 ; ++j ){
    hnew->SetBinContent( j, hnew->GetBinContent( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
    hnew->SetBinError(   j, hnew->GetBinError( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
  }
      

  return hnew;
}


 



TH1D *Rebin_dPhi_full(TH1D* hold){
  // New bin boundaries
  Double_t bin_bounds_phi[]=  {-1.507968,-1.382304,-1.25664,-1.130976,-1.005312,-0.879648,-0.753984,-0.62832,-0.502656,-0.376992,-0.251328,-0.125664,0.125664,0.251328,0.376992,0.502656,0.62832,0.753984,0.879648,1.005312,1.130976,1.25664,1.382304,1.507968,1.633632,1.82212,2.01062,2.19911,2.38761,2.57611,2.7646,2.9531,3.14159,3.33009,3.51858,3.70708,3.89557,4.08407,4.27257,4.46106,4.58673};

  Int_t new_bins = sizeof(bin_bounds_phi) / sizeof(bin_bounds_phi[0])-1; // Should be 35 or so now
  

  // important cross check. These need to be binboundaries of the old one!
  for (int i=0; i<new_bins+1 ; ++i ){
    bool binok=false;
    for (int j=1; j<=hold->GetNbinsX()+1 ; ++j ){
      if ( fabs( bin_bounds_phi[i] - hold->GetBinLowEdge(j) )<=1e-3 ) {
	binok=true;
	break;
      }
    }
    if (!binok ){
      cerr << "Bin edge " << bin_bounds_phi[i] << " is NOT a bin boundary" << endl;
      return 0;
    }
  }
    
  TString new_name = TString( hold->GetName() ) + "_Rebin" ;
  TH1D* hnew = (TH1D*) hold->Rebin ( new_bins, new_name, bin_bounds_phi );
  
  // Now renormalize per dPhi
  Double_t dx = hold->GetXaxis()->GetBinWidth(1);
  for (int j=1; j<=hnew->GetNbinsX()+1 ; ++j ){
    hnew->SetBinContent( j, hnew->GetBinContent( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
    hnew->SetBinError(   j, hnew->GetBinError( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
  }
      

  return hnew;
}
