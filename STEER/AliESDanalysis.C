//********************************************************************
//     Example (very basic for the moment) of the data analysis 
//                    using the ESD classes
//
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#ifndef __CINT__
  #include <Riostream.h>
  #include "TKey.h"
  #include "TFile.h"
  #include "TH2F.h"
  #include "TCanvas.h"
  #include "TStopwatch.h"

  #include "AliESD.h"
#endif

Int_t AliESDanalysis(Int_t nev=1) { 
   TH2F *tpcHist=new TH2F("tpcHist","dE/dX vs momentum",50,0.,2.,50,0.,400.);
   TH2F *itsHist=new TH2F("itsHits","dE/dX vs momentum",50,0.,2.,50,0.,200.);
   
   TFile *ef=TFile::Open("AliESDs.root");
   if (!ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}

   TStopwatch timer;
   Int_t rc=0,n=0;
   TKey *key=0;
   TIter next(ef->GetListOfKeys());

   //******* The loop over events
   while ((key=(TKey*)next())!=0) {
     cerr<<"Processing event number : "<<n++<<endl;
     AliESD *event=(AliESD*)key->ReadObj();

     Int_t ntrk=event->GetNumberOfTracks();
     cerr<<"Number of ESD tracks : "<<ntrk<<endl; 
     //****** The loop over tracks
     while (ntrk--) {
       AliESDtrack *t=event->GetTrack(ntrk);
       Double_t p=t->GetP();
       if (t->GetStatus()&AliESDtrack::kTPCin) {
	 Double_t dedx=t->GetTPCsignal();
         tpcHist->Fill(p,dedx,1);
       }
       if (t->GetStatus()&AliESDtrack::kITSin) {
	 Double_t dedx=t->GetITSsignal();
         itsHist->Fill(p,dedx,1);
       }
     } 
     delete event;
   }
   timer.Stop(); timer.Print();

   TCanvas *c1=new TCanvas("c1","",0,0,600,1200);
   c1->Divide(1,2);

   c1->cd(1);
   tpcHist->Draw();
   c1->cd(2);
   itsHist->Draw();

   ef->Close();

   return rc;
}
