#ifndef ALIEMCALTRACKINGQATASK_H
#define ALIEMCALTRACKINGQATASK_H

#include "AliAnalysisTaskEmcal.h"

class AliParticleContainer;
class THnSparse;
class TH3;

class AliEmcalTrackingQATask : public AliAnalysisTaskEmcal {

 public:
  AliEmcalTrackingQATask();
  AliEmcalTrackingQATask(const char *name); 
  virtual ~AliEmcalTrackingQATask();

  void                   UserCreateOutputObjects();
  void                   SetGeneratorLevelName(const char* name);
  void                   SetDetectorLevelName(const char* name);

 protected:
  Bool_t                 FillHistograms()                               ;
  void                   ExecOnce()                                     ;
  void                   AllocateFindableParticlesTHnSparse()           ;
  void                   AllocateMatchedParticlesTHnSparse()            ;
  void                   FillFindableParticlesTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt);
  void                   FillMatchedParticlesTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt,
						       Double_t trackEta, Double_t trackPhi, Double_t trackPt, Byte_t trackType);

  // Task configuration

  // Service fields (non-streamed)
  AliParticleContainer* fGeneratorLevel        ; //! generator level container
  AliParticleContainer* fDetectorLevel         ; //! detector level container

  // Histograms
  TH3***                fTracksAll             ; //! all tracks
  TH3***                fTracksSelected        ; //! selected tracks (e.g. remove injected signal in HIJING productions)
  TH3**                 fParticlesAllPhysPrim  ; //! all physical primary particles
  TH3**                 fParticlesSelected     ; //! selected physical primary particles (e.g. remove injected signal in HIJING productions)
  THnSparse*            fParticlesFindable     ; //! findable physical primary particles (use PDG and charge selection)
  THnSparse*            fParticlesMatched      ; //! primary particles matched to detector level tracks
  
 private:
  AliEmcalTrackingQATask(const AliEmcalTrackingQATask&);            // not implemented
  AliEmcalTrackingQATask &operator=(const AliEmcalTrackingQATask&); // not implemented

  ClassDef(AliEmcalTrackingQATask, 1) // Track QA task (efficiency and pt resolution)
};
#endif
