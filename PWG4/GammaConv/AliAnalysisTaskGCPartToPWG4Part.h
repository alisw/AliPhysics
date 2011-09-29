#ifndef AliAnalysisTaskGCPartToPWG4Part_cxx
#define AliAnalysisTaskGCPartToPWG4Part_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Svein Lindal

class TH1F;
class AliESDEvent;
class AliGammaConversionAODObject;
class AliAODConversionPhoton;
class AliAODPWG4ParticleCorrelation;
class AliAODPWG4Particle;
class TClonesArray;
class TString;
class AliMCAnalysisUtils;
class AliAODMCHeader;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskGCPartToPWG4Part : public AliAnalysisTaskSE {

public:

  AliAnalysisTaskGCPartToPWG4Part(); 
  AliAnalysisTaskGCPartToPWG4Part(const char *name);
  virtual ~AliAnalysisTaskGCPartToPWG4Part();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetDeltaAODFileName(TString string) { fDeltaAODFileName = string;}
  void SetGammaBranchName(TString string) { fAODBranchName = string; }
 

  void SetDebugLevel(Int_t debugLevel) { fDebugLevel = debugLevel; }
  Int_t GetDebugLevel() const { return fDebugLevel; }

  void SetGammaCutId(TString cut) { fGammaCutString = Form("GammaConv_%s", cut.Data());}
  void SetPionCutId(TString cut) { fPionCutString = Form("GammaConv_%s", cut.Data());}

  
 private:

  //Clean up
  void CleanUp();

  //Get the AOD event from whereever it might be accessible
  AliAODEvent * GetAODEvent();

  Bool_t BothTracksPresent(const AliAODConversionPhoton * const photon, const TClonesArray * const tracks) const;
  Bool_t BothGammaPresent(const AliAODConversionPhoton * const pion, const TClonesArray * const photons, const TClonesArray * const tracks) const;

  //Get Conversion gammas branch
  TClonesArray * GetConversionGammas(const AliAODEvent * aodEvent) const;
  TClonesArray * GetPions(const AliAODEvent * aodEvent) const;
  TClonesArray * GetAODBranch(const AliAODEvent * aodEvent, TString branchName) const;

  //Fill AOD tree with PWG4 particles
  AliAODPWG4ParticleCorrelation * AddToAOD(AliGammaConversionAODObject * aodO, TClonesArray * branch, TString detector);
  AliAODPWG4ParticleCorrelation * AddToAOD(AliAODConversionPhoton * aodO, TClonesArray * branch, TString detector);
  AliAODPWG4ParticleCorrelation * AddPionToAOD(AliAODConversionPhoton * pion, TClonesArray * branch, TString detector, TClonesArray * photons);  
  //Process conv gamma
  void ProcessConvGamma( const AliAODEvent * const aodEvent );

  TString     fDeltaAODFileName;//! File where Gamma Conv AOD is located, if not in default AOD
  TString     fGammaCutString;   //! The cut string of the conversion analysis used to produce input AOD
  TString     fPionCutString;   //! The cut string of the conversion analysis used to produce input AOD
  TString     fAODBranchName;
  TClonesArray * fAODPWG4Photons;
  TClonesArray * fAODPWG4Pi0;


  Int_t fDebugLevel;


  AliAnalysisTaskGCPartToPWG4Part(const AliAnalysisTaskGCPartToPWG4Part&); // not implemented
  AliAnalysisTaskGCPartToPWG4Part& operator=(const AliAnalysisTaskGCPartToPWG4Part&); // not implemented

  //Int_t CheckTag(AliAODPWG4ParticleCorrelation * particle, TClonesArray * tracks, TClonesArray * arrayMC, AliAODMCHeader * mcHeader);
  
  ClassDef(AliAnalysisTaskGCPartToPWG4Part, 1); // example of analysis
};

#endif
