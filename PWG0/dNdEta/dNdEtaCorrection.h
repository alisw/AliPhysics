/* $Id$ */

#ifndef DNDETACORRECTION_H
#define DNDETACORRECTION_H


// ------------------------------------------------------
//
// Class to handle corrections for dN/dEta measurements
//
// ------------------------------------------------------
//
// TODO:
// - add documentation
// - add status: generate or use maps
// - add functionality to set the bin sizes
// 

#include <TNamed.h>
#include <TFile.h>

#include <CorrectionMatrix2D.h>


class dNdEtaCorrection : public TNamed
{
protected:  
  
  CorrectionMatrix2D* fNtrackToNparticleCorrection; // handles the track-to-vertex correction
  CorrectionMatrix2D* fEventBiasCorrection;         // handles the event bias correction

  TH1F* fhVtxZAllEvents;
  TH1F* fhVtxZUsedEvents;

public:
  dNdEtaCorrection(Char_t* name="dndeta_correction");

  void FillEvent(Float_t vtx)     {fhVtxZAllEvents->Fill(vtx);}
  void FillUsedEvent(Float_t vtx) {fhVtxZUsedEvents->Fill(vtx);}
  void FillParticleAllEvents(Float_t vtx, Float_t eta)         {fEventBiasCorrection->FillGene(vtx, eta);}
  void FillParticleWhenUsedEvent(Float_t vtx, Float_t eta)     {fEventBiasCorrection->FillMeas(vtx, eta);
                                                                fNtrackToNparticleCorrection->FillGene(vtx, eta);}
  void FillParticleWhenMeasuredTrack(Float_t vtx, Float_t eta) {fNtrackToNparticleCorrection->FillMeas(vtx, eta);}

  void Finish();

  CorrectionMatrix2D* GetNtrackToNpraticleCorrection() {return fNtrackToNparticleCorrection;}
  CorrectionMatrix2D* GetEventBiasCorrection()         {return fEventBiasCorrection;}

  TH1F* GetVertexZHistogramAllEvents()  {return fhVtxZAllEvents;}
  TH1F* GetVertexZHistogramUsedEvents() {return fhVtxZUsedEvents;}

  virtual Long64_t Merge(TCollection* list);

  void    SaveHistograms();
  Bool_t  LoadHistograms(Char_t* fileName, Char_t* dir = "dndeta_correction");
  Bool_t  LoadCorrection(Char_t* fileName, Char_t* dir = "dndeta_correction") 
    {return LoadHistograms(fileName, dir);}
  
  void DrawHistograms();
  
  void RemoveEdges(Float_t cut=2, Int_t nBinsVtx=0, Int_t nBinsEta=0);
  
  Float_t GetCorrection(Float_t vtx, Float_t eta) 
    {return fNtrackToNparticleCorrection->GetCorrection(vtx, eta) * fEventBiasCorrection->GetCorrection(vtx, eta);}
  

  ClassDef(dNdEtaCorrection,0)
};

#endif
