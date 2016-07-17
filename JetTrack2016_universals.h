

  //*********************************************************
  // SET INTEGRAL LIMITS AND LABELS
  //*************************************************

  double etalim = 1.0;
  double philim = 1.0;

  TString etarangelabel = "Projected |#Delta#eta| < 1.0";
  TString phirangelabel = "Projected |#Delta#phi| < 1.0";

  float llimitphi,rlimitphi,llimiteta,rlimiteta,nbins;
  //*********************************************************

Int_t nbounds_phi = 18;
Int_t nbounds_eta = 20;
Int_t nbounds_eta2 = 18;
 

Double_t bin_bounds_phi[18]   =  {-1.50796, -1.00531,-0.879646, -.75398, -0.628319,-0.502655, -0.376991, -0.251327, -0.125664, 0.125664, 0.251327, 0.376991, 0.502655, 0.628319,.75398, 0.879646, 1.00531,1.50796};

Double_t bin_bounds_eta[20]
    =   {-2.5,-2.,-1.5,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.5,2.,2.5};


 Double_t bin_bounds_eta2[18]
    =   {-2.,-1.5,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.5,2.};

Double_t pTbins[10] = {0.,0.5,1.,2.,3.,4.,8.,12.,16.,20.};

 

 

