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
   TFile *ef=TFile::Open("AliESDs.root");
   if (!ef || !ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}
   AliESDEvent * event = new AliESDEvent;
   TTree* tree = (TTree*) ef->Get("esdTree");
   if (!tree) {cerr<<"no ESD tree found\n"; return 1;};
   event->ReadFromTree(tree);
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
       Float_t clsX  = pmdtr->GetClusterX();
       Float_t clsY  = pmdtr->GetClusterY();
       Float_t clsZ  = pmdtr->GetClusterZ();
       Float_t ncell = pmdtr->GetClusterCells();
       Float_t adc   = pmdtr->GetClusterADC();
       Float_t pid   = pmdtr->GetClusterPID();
       
     }
   }

   delete event;

   timer.Stop();
   timer.Print();

   return 0;
}