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

  fNtrackToNparticleCorrection = new CorrectionMatrix2D("nTrackToNPart", "nTrackToNPart",80,-20,20,120,-6,6);
  fEventBiasCorrection         = new CorrectionMatrix2D("EventBias",    "EventBias ",80,-20,20,120,-6,6);

  fNtrackToNparticleCorrection ->SetAxisTitles("vtx z [cm]", "#eta");
  fEventBiasCorrection         ->SetAxisTitles("vtx z [cm]", "#eta");

  fhVtxZAllEvents  = new TH1F("vtx_z_all_events", "vtx_z_all_events",80,-20,20);
  fhVtxZUsedEvents = new TH1F("vtx_z_used_events","vtx_z_used_events",80,-20,20);

  fhVtxZAllEvents ->Sumw2();
  fhVtxZUsedEvents->Sumw2();
}

//____________________________________________________________________
void
dNdEtaCorrection::Finish() {  
  //
  // finish method
  //
  // divide the histograms in the CorrectionMatrix2D objects to get the corrections


  fNtrackToNparticleCorrection->Divide();

  // normalize event bias histograms to the number of events

  TH2F* meas = fEventBiasCorrection->GetMeasuredHistogram();
  TH2F* gene = fEventBiasCorrection->GetGeneratedHistogram();
  for (Int_t i_vtx=0; i_vtx<=meas->GetNbinsX(); i_vtx++) {
    Int_t nEventsAll  = (Int_t)fhVtxZAllEvents->GetBinContent(i_vtx);
    Int_t nEventsUsed = (Int_t)fhVtxZUsedEvents->GetBinContent(i_vtx);

    if (nEventsAll<10)  nEventsAll=0;
    if (nEventsUsed<10) nEventsUsed=0;

    for (Int_t i_eta=0; i_eta<=meas->GetNbinsY(); i_eta++) {
      Float_t valueMeas=0;
      Float_t errorMeas=0;

      Float_t valueGene=0;
      Float_t errorGene=0;
      
      if (nEventsUsed!=0) {
	valueMeas = meas->GetBinContent(i_vtx, i_eta)/Float_t(nEventsUsed);
	errorMeas = meas->GetBinError(i_vtx, i_eta)/Float_t(nEventsUsed);
      }
      meas->SetBinContent(i_vtx, i_eta, valueMeas);
      meas->SetBinError(i_vtx, i_eta, errorMeas);

      if (nEventsAll!=0) {
	valueGene = gene->GetBinContent(i_vtx, i_eta)/Float_t(nEventsAll);
	errorGene = gene->GetBinError(i_vtx, i_eta)/Float_t(nEventsAll);
      }
      gene->SetBinContent(i_vtx, i_eta, valueGene);
      gene->SetBinError(i_vtx, i_eta, errorGene);

    }
  }
  fEventBiasCorrection->SetMeasuredHistogram(meas);
  fEventBiasCorrection->SetGeneratedHistogram(gene);
  
  fEventBiasCorrection->Divide();
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
  TList* collectionVtxHistAllEvents  = new TList;
  TList* collectionVtxHistUsedEvents = new TList;

  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    dNdEtaCorrection* entry = dynamic_cast<dNdEtaCorrection*> (obj);
    if (entry == 0) 
      continue;

    collectionNtrackToNparticle ->Add(entry->GetNtrackToNpraticleCorrection());
    collectionEventBias         ->Add(entry->GetEventBiasCorrection());
    collectionVtxHistAllEvents  ->Add(entry->GetVertexZHistogramAllEvents());
    collectionVtxHistUsedEvents ->Add(entry->GetVertexZHistogramUsedEvents());

    count++;
  }
  fNtrackToNparticleCorrection ->Merge(collectionNtrackToNparticle);
  fEventBiasCorrection         ->Merge(collectionEventBias);

  fhVtxZAllEvents ->Merge(collectionVtxHistAllEvents );
  fhVtxZUsedEvents->Merge(collectionVtxHistUsedEvents);

  
  delete collectionNtrackToNparticle;
  delete collectionEventBias;
  delete collectionVtxHistAllEvents;
  delete collectionVtxHistUsedEvents;

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

  fhVtxZAllEvents ->Write();
  fhVtxZUsedEvents->Write();

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
