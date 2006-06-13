/* $Id$ */

#include "AlidNdEtaCorrection.h"

#include <TCanvas.h>

//____________________________________________________________________
ClassImp(AlidNdEtaCorrection)

//____________________________________________________________________
AlidNdEtaCorrection::AlidNdEtaCorrection(Char_t* name) 
  : TNamed(name, name)
{  
  // constructor
  //

  fNtrackToNparticleCorrection = new CorrectionMatrix2D("nTrackToNPart", "nTrackToNPart",80,-20,20,120,-6,6);

  Float_t binLimitsN[]   = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 
			    10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 25.5, 30.5, 40.5, 50.5, 100.5, 300.5};
  Float_t binLimitsVtx[] = {-20,-15,-10,-6,-3,0,3,6,10,15,20};
  
  fVertexRecoCorrection        = new CorrectionMatrix2D("vtxReco",       "vtxReco",10,binLimitsVtx ,22,binLimitsN);

  fTriggerBiasCorrection       = new CorrectionMatrix2D("triggerBias",   "triggerBias",120,-6,6,100, 0, 10);

  fNtrackToNparticleCorrection ->SetAxisTitles("vtx z [cm]", "#eta");
  fVertexRecoCorrection        ->SetAxisTitles("vtx z [cm]", "n particles/tracks/tracklets?");
  
  fTriggerBiasCorrection       ->SetAxisTitles("#eta", "p_{T} [GeV/c]");
}

//____________________________________________________________________
void
AlidNdEtaCorrection::Finish(Int_t nEventsAll, Int_t nEventsTriggered) {  
  //
  // finish method
  //
  // divide the histograms in the CorrectionMatrix2D objects to get the corrections

  
  fNtrackToNparticleCorrection->Divide();

  fVertexRecoCorrection->Divide();

  fTriggerBiasCorrection->GetMeasuredHistogram()->Scale(Double_t(nEventsTriggered)/Double_t(nEventsAll));
  fTriggerBiasCorrection->Divide();

}

//____________________________________________________________________
Long64_t 
AlidNdEtaCorrection::Merge(TCollection* list) {
  // Merge a list of dNdEtaCorrection objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of measured and generated histograms
  TList* collectionNtrackToNparticle = new TList;
  TList* collectionVertexReco        = new TList;
  TList* collectionTriggerBias       = new TList;

  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AlidNdEtaCorrection* entry = dynamic_cast<AlidNdEtaCorrection*> (obj);
    if (entry == 0) 
      continue;

    collectionNtrackToNparticle ->Add(entry->GetNtrackToNpraticleCorrection());
    collectionVertexReco        ->Add(entry->GetVertexRecoCorrection());
    collectionTriggerBias        ->Add(entry->GetTriggerBiasCorrection());

    count++;
  }
  fNtrackToNparticleCorrection ->Merge(collectionNtrackToNparticle);
  fVertexRecoCorrection        ->Merge(collectionVertexReco);
  fTriggerBiasCorrection        ->Merge(collectionTriggerBias);
  
  delete collectionNtrackToNparticle;
  delete collectionVertexReco;
  delete collectionTriggerBias;

  return count+1;
}


//____________________________________________________________________
Bool_t
AlidNdEtaCorrection::LoadHistograms(Char_t* fileName, Char_t* dir) {
  //
  // loads the histograms
  //

  fNtrackToNparticleCorrection ->LoadHistograms(fileName, dir);
  fVertexRecoCorrection        ->LoadHistograms(fileName, dir);
  fTriggerBiasCorrection       ->LoadHistograms(fileName, dir);
  
  return kTRUE;
}


//____________________________________________________________________
void
AlidNdEtaCorrection::SaveHistograms() {
  //
  // save the histograms
  //

  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());

  fNtrackToNparticleCorrection ->SaveHistograms();
  fVertexRecoCorrection        ->SaveHistograms();
  fTriggerBiasCorrection       ->SaveHistograms();

  gDirectory->cd("../");
}

//____________________________________________________________________
void AlidNdEtaCorrection::DrawHistograms()
{
  //
  // call the draw histogram method of the two CorrectionMatrix2D objects

  fNtrackToNparticleCorrection ->DrawHistograms();
  fVertexRecoCorrection        ->DrawHistograms();
  fTriggerBiasCorrection       ->DrawHistograms();

  
}
