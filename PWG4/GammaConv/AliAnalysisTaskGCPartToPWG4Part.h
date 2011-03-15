#ifndef AliAnalysisTaskGCPartToPWG4Part_cxx
#define AliAnalysisTaskGCPartToPWG4Part_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Svein Lindal

class TH1F;
class AliESDEvent;
class AliGammaConversionAODObject;
class AliAODConversionParticle;
class AliAODPWG4ParticleCorrelation;
class AliAODPWG4Particle;
class TClonesArray;
class TString;
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
  
 private:

  //Clean up
  void CleanUp();

  //Get the AOD event from whereever it might be accessible
  AliAODEvent * GetAODEvent();

  //Get Conversion gammas branch
  TClonesArray * GetConversionGammas(const AliAODEvent * aodEvent);

  //Fill AOD tree with PWG4 particles
  AliAODPWG4ParticleCorrelation * AddToAOD(AliGammaConversionAODObject * aodO, TClonesArray * branch, TString detector);
  AliAODPWG4ParticleCorrelation * AddToAOD(AliAODConversionParticle * aodO, TClonesArray * branch, TString detector);
  
  //Process conv gamma
  void ProcessConvGamma( const AliAODEvent * const aodEvent );

  TString     fDeltaAODFileName;//! File where Gamma Conv AOD is located, if not in default AOD
  TString     fAODBranchName;
  TClonesArray * fAODPWG4Particles;

  Int_t fDebugLevel;


  AliAnalysisTaskGCPartToPWG4Part(const AliAnalysisTaskGCPartToPWG4Part&); // not implemented
  AliAnalysisTaskGCPartToPWG4Part& operator=(const AliAnalysisTaskGCPartToPWG4Part&); // not implemented

  //Int_t CheckTag(AliAODPWG4ParticleCorrelation * particle, TClonesArray * tracks, TClonesArray * arrayMC, AliAODMCHeader * mcHeader);
  
  ClassDef(AliAnalysisTaskGCPartToPWG4Part, 1); // example of analysis
};

#endif
