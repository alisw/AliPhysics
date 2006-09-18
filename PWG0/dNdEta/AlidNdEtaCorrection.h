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
#include <AliLog.h>

class AlidNdEtaCorrection : public TNamed
{
public:
  AlidNdEtaCorrection();
  AlidNdEtaCorrection(const Char_t* name, const Char_t* title);

  ~AlidNdEtaCorrection();

  // fTrack2ParticleCorrection
  void FillParticle(Float_t vtx, Float_t eta, Float_t pt)                  {fTrack2ParticleCorrection->FillGene(vtx, eta, pt);}
  void FillParticleWhenMeasuredTrack(Float_t vtx, Float_t eta, Float_t pt) {fTrack2ParticleCorrection->FillMeas(vtx, eta, pt);}

  // fVertexRecoCorrection & fTriggerBiasCorrection
  void FillEventWithTriggerWithReconstructedVertex(Float_t vtx, Float_t n) {fVertexRecoCorrection->FillMeas(vtx, n);}
  void FillEventWithTrigger(Float_t vtx, Float_t n);                         
  void FillEventAll(Float_t vtx, Float_t n, Char_t* opt="INEL"); 
    
  void Finish();

  AliCorrectionMatrix3D* GetTrack2ParticleCorrection()  {return fTrack2ParticleCorrection;}
  AliCorrectionMatrix2D* GetVertexRecoCorrection()      {return fVertexRecoCorrection;}
  AliCorrectionMatrix2D* GetTriggerBiasCorrectionINEL() {return fTriggerBiasCorrectionMBToINEL;}
  AliCorrectionMatrix2D* GetTriggerBiasCorrectionNSD()  {return fTriggerBiasCorrectionMBToNSD;}
  AliCorrectionMatrix2D* GetTriggerBiasCorrectionND()   {return fTriggerBiasCorrectionMBToND;}
  AliCorrectionMatrix2D* GetTriggerBiasCorrection(Char_t* opt="INEL");

  virtual Long64_t Merge(TCollection* list);

  void    SaveHistograms();
  Bool_t  LoadHistograms(const Char_t* fileName, const Char_t* dir = "dndeta_correction");
  Bool_t  LoadCorrection(const Char_t* fileName, const Char_t* dir = "dndeta_correction") {return LoadHistograms(fileName, dir);}
  
  void DrawHistograms();
  
  //  void RemoveEdges(Float_t cut=2, Int_t nBinsVtx=0, Int_t nBinsEta=0);
  
  Float_t GetTrack2ParticleCorrection(Float_t vtx, Float_t eta, Float_t pt)
    {return fTrack2ParticleCorrection->GetCorrection(vtx, eta, pt);}
  
  Float_t GetVertexRecoCorrection     (Float_t vtx, Float_t n) {return fVertexRecoCorrection->GetCorrection(vtx, n);}
  Float_t GetTriggerBiasCorrectionINEL(Float_t vtx, Float_t n) {return fTriggerBiasCorrectionMBToINEL->GetCorrection(vtx, n);} 
  Float_t GetTriggerBiasCorrectionNSD (Float_t vtx, Float_t n) {return fTriggerBiasCorrectionMBToNSD->GetCorrection(vtx, n);} 
  Float_t GetTriggerBiasCorrectionND  (Float_t vtx, Float_t n) {return fTriggerBiasCorrectionMBToND->GetCorrection(vtx, n);} 
  Float_t GetTriggerBiasCorrection    (Float_t vtx, Float_t n, Char_t* opt="INEL");

  Float_t GetMeasuredFraction(Float_t ptCutOff, Float_t eta = -100, Bool_t debug = kFALSE);


  void ReduceInformation();

protected:
  AliCorrectionMatrix3D* fTrack2ParticleCorrection;       //-> handles the track-to-particle correction, function of vtx_z, eta, pt
  AliCorrectionMatrix2D* fVertexRecoCorrection;           //-> handles the vertex reconstruction efficiency, function of n_clustersITS and vtx_z
  AliCorrectionMatrix2D* fTriggerBiasCorrectionMBToINEL;  //-> handles the trigger bias MB->INEL, function of n and vtx_z
  AliCorrectionMatrix2D* fTriggerBiasCorrectionMBToNSD;   //-> handles the trigger bias MB->NSD,  function of n and vtx_z
  AliCorrectionMatrix2D* fTriggerBiasCorrectionMBToND;    //-> handles the trigger bias MB->ND,   function of n and vtx_z

  //Long64_t fNEvents;
  //Long64_t fNTriggeredEvents;

private:
  AlidNdEtaCorrection(const AlidNdEtaCorrection&);
  AlidNdEtaCorrection& operator=(const AlidNdEtaCorrection&);

  ClassDef(AlidNdEtaCorrection, 1)
};

#endif
