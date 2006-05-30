/* $Id$ */

#include "dNdEtaCorrection.h"

#include <TCanvas.h>

//____________________________________________________________________
ClassImp(dNdEtaCorrection)

//____________________________________________________________________
dNdEtaCorrection::dNdEtaCorrection(Char_t* name) 
  : TNamed(name, name)
{  
  // constructor
  //

  fNtrackToNparticleCorrection = new CorrectionMatrix2D("nTrackToNPart","",80,-20,20,120,-6,6);
  fEventBiasCorrection         = new CorrectionMatrix2D("EventBias",    "",80,-20,20,120,-6,6);

  fNtrackToNparticleCorrection ->SetAxisTitles("vtx z [cm]", "#eta");
  fEventBiasCorrection         ->SetAxisTitles("vtx z [cm]", "#eta");

}

//____________________________________________________________________
void
dNdEtaCorrection::Finish() {  
  //
  // finish method
  //
  // divide the histograms in the CorrectionMatrix2D objects to get the corrections


  fNtrackToNparticleCorrection->Divide();
  fEventBiasCorrection        ->Divide();

}

//____________________________________________________________________
Long64_t 
dNdEtaCorrection::Merge(TCollection* list) {
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
  TList* collectionEventBias         = new TList;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    dNdEtaCorrection* entry = dynamic_cast<dNdEtaCorrection*> (obj);
    if (entry == 0) 
      continue;

    collectionNtrackToNparticle ->Add(entry->GetNtrackToNpraticleCorrection());
    collectionEventBias         ->Add(entry->GetEventBiasCorrection());
    
    count++;
  }
  fNtrackToNparticleCorrection ->Merge(collectionNtrackToNparticle);
  fEventBiasCorrection         ->Merge(collectionEventBias);
  
  delete collectionNtrackToNparticle;
  delete collectionEventBias;

  return count+1;
}




//____________________________________________________________________
void
dNdEtaCorrection::RemoveEdges(Float_t cut, Int_t nBinsVtx, Int_t nBinsEta) {
  //
  // removes the edges of the correction maps
  //

  fNtrackToNparticleCorrection ->RemoveEdges(cut, nBinsVtx, nBinsEta);
  fEventBiasCorrection         ->RemoveEdges(cut, nBinsVtx, nBinsEta);
}

//____________________________________________________________________
Bool_t
dNdEtaCorrection::LoadHistograms(Char_t* fileName, Char_t* dir) {
  //
  // loads the histograms
  //

  fNtrackToNparticleCorrection ->LoadHistograms(fileName, dir);
  fEventBiasCorrection         ->LoadHistograms(fileName, dir);
  
  return kTRUE;
}


//____________________________________________________________________
void
dNdEtaCorrection::SaveHistograms() {
  //
  // save the histograms
  //

  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());

  fNtrackToNparticleCorrection ->SaveHistograms();
  fEventBiasCorrection         ->SaveHistograms();

  gDirectory->cd("../");
}

//____________________________________________________________________
void dNdEtaCorrection::DrawHistograms()
{
  //
  // call the draw histogram method of the two CorrectionMatrix2D objects

  fNtrackToNparticleCorrection ->DrawHistograms();
  fEventBiasCorrection         ->DrawHistograms();

  fEventBiasCorrection->Get1DCorrection("x")->Draw();
  
}
