//********************************************************************
//     Example (very naive for the moment) of the data analysis 
//                    using the ESD classes
//     Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//     Modified by Basanta K. Nandi for PMD analysis
//********************************************************************

#if !defined( __CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include "TTree.h"
  #include "TFile.h"
  #include "TH1F.h"
  #include "TCanvas.h"
  #include "TStyle.h"
  #include "TStopwatch.h"

  #include "AliESD.h"
#endif

extern TStyle *gStyle;

Int_t AliPMDesdanal() { 
   TStopwatch timer;

   gStyle->SetOptStat(111110);
   gStyle->SetOptFit(1);

//****** File with the ESD
   TFile *ef=TFile::Open("AliESDcheck.root");
   if (!ef || !ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}
   AliESD* event = new AliESD;
   TTree* tree = (TTree*) ef->Get("esdTree");
   if (!tree) {cerr<<"no ESD tree found\n"; return 1;};
   tree->SetBranchAddress("ESD", &event);

   Int_t n=0;

//******* The loop over events
   while (tree->GetEvent(n)) {
     cout<<endl<<"Processing event number : "<<n++<<endl;


     Int_t npmdcl=event->GetNumberOfPmdTracks();
     cout<<"Number of PMD tracks : "<<npmdcl<<endl; 

     //****** The loop over PMD clusters
     while (npmdcl--) {
       AliESDPmdTrack *pmdtr = event->GetPmdTrack(npmdcl);
       
       Int_t   det   = pmdtr->GetDetector();
       Float_t theta = pmdtr->GetTheta();
       Float_t phi   = pmdtr->GetPhi();
       Float_t adc   = pmdtr->GetClusterADC();
       Float_t pid   = pmdtr->GetClusterPID();
       
     }
   }
   timer.Stop();
   timer.Print();

   return 0;
}
