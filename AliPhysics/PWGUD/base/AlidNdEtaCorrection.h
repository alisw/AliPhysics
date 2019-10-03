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
// - add functionality to set the bin sizes
// - update MERge function
//

#include <TCollection.h>
#include <TNamed.h>
#include "AliPWG0Helper.h"

class AliCorrection;
class TH1;

class AlidNdEtaCorrection : public TNamed
{
public:
  enum CorrectionType {
    kNone = 0,
    kTrack2Particle,  // measured events
    kVertexReco,      // MB sample
    kINEL,
    kNSD,
    kND,
    kOnePart
  };

  AlidNdEtaCorrection();
  AlidNdEtaCorrection(const Char_t* name, const Char_t* title, AliPWG0Helper::AnalysisMode analysis = (AliPWG0Helper::AnalysisMode) (AliPWG0Helper::kTPC | AliPWG0Helper::kFieldOn));

  virtual Long64_t Merge(TCollection* list);

  ~AlidNdEtaCorrection();

  void FillMCParticle(Float_t vtx, Float_t eta, Float_t pt, Bool_t trigger, Bool_t vertex, Int_t processType);
  void FillTrackedParticle(Float_t vtx, Float_t eta, Float_t pt, Double_t weight=1.);
  void FillEvent(Float_t vtx, Float_t n, Bool_t trigger, Bool_t vertex, Int_t processType);

  void Finish();

  AliCorrection* GetTrack2ParticleCorrection()  {return fTrack2ParticleCorrection;}
  AliCorrection* GetVertexRecoCorrection()      {return fVertexRecoCorrection;}
  AliCorrection* GetTriggerBiasCorrectionINEL() {return fTriggerBiasCorrectionMBToINEL;}
  AliCorrection* GetTriggerBiasCorrectionNSD()  {return fTriggerBiasCorrectionMBToNSD;}
  AliCorrection* GetTriggerBiasCorrectionND()   {return fTriggerBiasCorrectionMBToND;}
  AliCorrection* GetTriggerBiasCorrectionOnePart()   {return fTriggerBiasCorrectionMBToOnePart;}
  AliCorrection* GetCorrection(CorrectionType correctionType);

  void    Reset(void);
  void    Add(AlidNdEtaCorrection* aCorrectionsToAdd, Float_t c=1);
  void    Scale(Float_t c);

  void    SaveHistograms();
  Bool_t  LoadHistograms(const Char_t* dir = 0);
  void    DrawHistograms();
  void    DrawOverview(const char* canvasName = 0);

  Float_t GetMeasuredFraction(CorrectionType correctionType, Float_t ptCutOff, Float_t eta = -100, Int_t vertexBegin = -1, Int_t vertexEnd = -1, Bool_t debug = kFALSE);
  TH1*    GetMeasuredEventFraction(CorrectionType correctionType, Int_t multCut);

  void    ReduceInformation();

protected:
  AliCorrection* fTrack2ParticleCorrection;       //-> handles the track-to-particle correction (only track level (vtx_z, eta, pt))
  AliCorrection* fVertexRecoCorrection;           //-> handles the vertex reconstruction efficiency, (n, vtx_z)
  AliCorrection* fTriggerBiasCorrectionMBToINEL;  //-> handles the trigger bias MB->INEL, function of n and vtx_z
  AliCorrection* fTriggerBiasCorrectionMBToNSD;   //-> handles the trigger bias MB->NSD,  function of n and vtx_z
  AliCorrection* fTriggerBiasCorrectionMBToND;    //-> handles the trigger bias MB->ND,   function of n and vtx_z
  AliCorrection* fTriggerBiasCorrectionMBToOnePart;    //-> handles the trigger bias MB->OnePart,   function of n and vtx_z

private:
  AlidNdEtaCorrection(const AlidNdEtaCorrection&);
  AlidNdEtaCorrection& operator=(const AlidNdEtaCorrection&);

  ClassDef(AlidNdEtaCorrection, 2)
};

#endif
