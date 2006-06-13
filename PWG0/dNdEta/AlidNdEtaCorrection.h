/* $Id$ */

#ifndef ALIDNDETACORRECTION_H
#define ALIDNDETACORRECTION_H


// ------------------------------------------------------
//
// Class to handle corrections for dN/dEta measurements
//
// ------------------------------------------------------
//
// TODO:
// - make the ntrack to npart correction 3D
// - add documentation
// - add status: generate or use maps
// - add functionality to set the bin sizes
// 

#include <TNamed.h>
#include <TFile.h>

#include <CorrectionMatrix2D.h>


class AlidNdEtaCorrection : public TNamed
{
protected:  
  
  CorrectionMatrix2D* fNtrackToNparticleCorrection; // handles the track-to-vertex correction
  CorrectionMatrix2D* fVertexRecoCorrection;        // handles the vertex reco (n tracks vs vtx)

  CorrectionMatrix2D* fTriggerBiasCorrection;          // MB to desired sample

public:
  AlidNdEtaCorrection(Char_t* name="dndeta_correction");

  void FillEvent(Float_t vtx, Float_t n)                        {fVertexRecoCorrection->FillGene(vtx, n);}
  void FillEventWithReconstructedVertex(Float_t vtx, Float_t n) {fVertexRecoCorrection->FillMeas(vtx, n);}
  
  void FillParticle(Float_t vtx, Float_t eta, Float_t pt=0)                  {fNtrackToNparticleCorrection->FillGene(vtx, eta);}
  void FillParticleWhenMeasuredTrack(Float_t vtx, Float_t eta, Float_t pt=0) {fNtrackToNparticleCorrection->FillMeas(vtx, eta);}
  
  void FillParticleAllEvents(Float_t eta, Float_t pt=0)          {fTriggerBiasCorrection->FillGene(eta, pt);}
  void FillParticleWhenEventTriggered(Float_t eta, Float_t pt=0) {fTriggerBiasCorrection->FillMeas(eta, pt);}

  void Finish(Int_t nEventsAll = 1, Int_t nEventsTriggered = 1);

  CorrectionMatrix2D* GetNtrackToNpraticleCorrection() {return fNtrackToNparticleCorrection;}
  CorrectionMatrix2D* GetVertexRecoCorrection()        {return fVertexRecoCorrection;}
  CorrectionMatrix2D* GetTriggerBiasCorrection()       {return fTriggerBiasCorrection;}

  virtual Long64_t Merge(TCollection* list);

  void    SaveHistograms();
  Bool_t  LoadHistograms(Char_t* fileName, Char_t* dir = "dndeta_correction");
  Bool_t  LoadCorrection(Char_t* fileName, Char_t* dir = "dndeta_correction") 
    {return LoadHistograms(fileName, dir);}
  
  void DrawHistograms();
  
  //  void RemoveEdges(Float_t cut=2, Int_t nBinsVtx=0, Int_t nBinsEta=0);
  
  Float_t GetNtracksToNpartCorrection(Float_t vtx, Float_t eta, Float_t pt) 
    {return fNtrackToNparticleCorrection->GetCorrection(vtx, eta);}  
  
  Float_t GetVertexRecoCorrection(Float_t vtx, Float_t n) {return fVertexRecoCorrection->GetCorrection(vtx, n);}

  Float_t GetTriggerBiasCorrection(Float_t eta, Float_t pt=0) {return fTriggerBiasCorrection->GetCorrection(eta, pt);}
  

  ClassDef(AlidNdEtaCorrection,0)
};

#endif
