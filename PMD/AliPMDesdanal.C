//********************************************************************
//     Example (very naive for the moment) of the data analysis 
//                    using the ESD classes
//     Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//     Modified by Basanta K. Nandi for PMD analysis
//********************************************************************

#if !defined( __CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include "TKey.h"
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

   Int_t n=0;
   TKey *key=0;
   TIter next(ef->GetListOfKeys());

//******* The loop over events
   while ((key=(TKey*)next())!=0) {
     cout<<endl<<"Processing event number : "<<n++<<endl;

     AliESD *event=(AliESD*)key->ReadObj();


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
