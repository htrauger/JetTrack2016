To use:

place all files in this folder into a subdirectory
include getTrkCorr.h header file (which should also grab the other header file)
make a trkCorr object at the start of your program somewhere:

TrkCorr* trkCorr = new TrkCorr("*Name of subdirectory*/"); //must have the '/' after the subdirectory name

To call get the correction, which is applied multiplicatively, use:
Note that you should use hiBin, NOT centrality (i.e. 0-199 range)

trkCorr->getTrkCorr(pt,eta,0,hiBin,0);


Tracking selection to be used:

|eta|<2.4, 0.5<pt<300 GeV

if(highPurity[j]!=1) continue;
if(trkPtError[j]/trkPt[j]>0.3 || TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;

float Et = (pfHcal[j]+pfEcal[j])/TMath::CosH(trkEta[j]);
if(!(trkPt[j]<20 || (Et>0.2*trkPt[j] && Et>trkPt[j]-80))) continue;


