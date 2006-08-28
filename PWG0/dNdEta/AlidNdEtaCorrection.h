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
// - update MERge function
//

#include <TNamed.h>

#include <AliCorrectionMatrix2D.h>
#include <AliCorrectionMatrix3D.h>

class AlidNdEtaCorrection : public TNamed
{
public:
  AlidNdEtaCorrection();
  AlidNdEtaCorrection(const Char_t* name, const Char_t* title);

  ~AlidNdEtaCorrection();

  // fVertexRecoCorrection, fTriggerCorrection
  void FillEvent(Float_t vtx, Float_t n)                        {fTriggerCorrection->FillGene(vtx, n);}
  void FillEventWithTrigger(Float_t vtx, Float_t n)             {fVertexRecoCorrection->FillGene(vtx, n); fTriggerCorrection->FillMeas(vtx, n);}
  void FillEventWithTriggerWithReconstructedVertex(Float_t vtx, Float_t n) {fVertexRecoCorrection->FillMeas(vtx, n);}

  // fTrack2ParticleCorrection
  void FillParticle(Float_t vtx, Float_t eta, Float_t pt)                  {fTrack2ParticleCorrection->FillGene(vtx, eta, pt);}
  void FillParticleWhenMeasuredTrack(Float_t vtx, Float_t eta, Float_t pt) {fTrack2ParticleCorrection->FillMeas(vtx, eta, pt);}

  // fTriggerBiasCorrection
  void FillParticleAllEvents(Float_t eta, Float_t pt)          {fTriggerBiasCorrection->FillGene(eta, pt);}
  void FillParticleWhenEventTriggered(Float_t eta, Float_t pt) {fTriggerBiasCorrection->FillMeas(eta, pt);}

  //void IncreaseEventCount() { fNEvents++; }
  //void IncreaseTriggeredEventCount() { fNTriggeredEvents++; }

  void Finish();

  AliCorrectionMatrix3D* GetTrack2ParticleCorrection()    {return fTrack2ParticleCorrection;}
  AliCorrectionMatrix2D* GetVertexRecoCorrection()        {return fVertexRecoCorrection;}
  AliCorrectionMatrix2D* GetTriggerBiasCorrection()       {return fTriggerBiasCorrection;}
  AliCorrectionMatrix2D* GetTriggerCorrection()       {return fTriggerCorrection;}

  virtual Long64_t Merge(TCollection* list);

  void    SaveHistograms();
  Bool_t  LoadHistograms(const Char_t* fileName, const Char_t* dir = "dndeta_correction");
  Bool_t  LoadCorrection(const Char_t* fileName, const Char_t* dir = "dndeta_correction")
    {return LoadHistograms(fileName, dir);}
  
  void DrawHistograms();
  
  //  void RemoveEdges(Float_t cut=2, Int_t nBinsVtx=0, Int_t nBinsEta=0);
  
  Float_t GetTrack2ParticleCorrection(Float_t vtx, Float_t eta, Float_t pt)
    {return fTrack2ParticleCorrection->GetCorrection(vtx, eta, pt);}

  Float_t GetVertexRecoCorrection(Float_t vtx, Float_t n) {return fVertexRecoCorrection->GetCorrection(vtx, n);}

  Float_t GetTriggerCorrection(Float_t vtx, Float_t n) {return fTriggerCorrection->GetCorrection(vtx, n);}

  Float_t GetTriggerBiasCorrection(Float_t eta, Float_t pt=0) {return fTriggerBiasCorrection->GetCorrection(eta, pt);}

  Float_t GetMeasuredFraction(Float_t ptCutOff, Float_t eta = -1, Bool_t debug = kFALSE);

  //void SetNEvents(Long64_t events) { fNEvents = events; }

  Int_t GetNevents();
  Int_t GetNtriggeredEvents();

  void ReduceInformation();

protected:
  AliCorrectionMatrix3D* fTrack2ParticleCorrection; //-> handles the track-to-particle correction, function of vtx_z, eta, pt
  AliCorrectionMatrix2D* fVertexRecoCorrection;     //-> handles the vertex reconstruction efficiency, function of n_clustersITS and vtx_z
  AliCorrectionMatrix2D* fTriggerCorrection;       //-> handles the trigger efficiency, function of n_clustersITS and vtx_z

  AliCorrectionMatrix2D* fTriggerBiasCorrection;          //-> MB to desired sample

  //Long64_t fNEvents;
  //Long64_t fNTriggeredEvents;

private:
  AlidNdEtaCorrection(const AlidNdEtaCorrection&);
  AlidNdEtaCorrection& operator=(const AlidNdEtaCorrection&);

  ClassDef(AlidNdEtaCorrection, 1)
};

#endif
