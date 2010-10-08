#ifndef AliAnalysisTaskGammaJet_cxx
#define AliAnalysisTaskGammaJet_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class AliESDEvent;
class AliGammaConversionAODObject;
class AliAODPWG4ParticleCorrelation;
class TClonesArray;
class TString;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskGammaJet : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskGammaJet(); 
  AliAnalysisTaskGammaJet(const char *name);
  virtual ~AliAnalysisTaskGammaJet() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetDeltaAODFileName(TString string) { fDeltaAODFileName = string;}
  //Move all of these somewhere
  inline Float_t GetMinPt() const { return fMinPt;} 
  void SetMinPt( const Float_t pt ) { fMinPt = pt; }
  inline Float_t GetConeSize () const { return fConeSize; }
  void SetConeSize ( const Float_t cs ) { fConeSize = cs; }
  inline Float_t GetPtThreshold () const { return fPtThreshold; }
  void SetPtThreshold ( Float_t ptt ) { fPtThreshold = ptt; }
 



  
 private:

  //Clean up
  void CleanUp();

  //Get the AOD event from whereever it might be accessible
  AliAODEvent * GetAODEvent();

  //Get Conversion gammas branch
  TClonesArray * GetConversionGammas(const AliAODEvent * aodEvent);

  //Create PWG4 particles from the aod objects
  AliAODPWG4ParticleCorrelation * PWG4PartFromGammaConvAODObject(AliGammaConversionAODObject * gcObject, TString detector);

  //Fill AOD tree with PWG4 particles
  void FillPWG4PartCorrBranch( TClonesArray * gcBranch, TClonesArray * partCorrBranch, TString detector);
  
  //Is particle isolated
  Bool_t IsIsolated( AliAODPWG4ParticleCorrelation * particle, TClonesArray * tracks, Float_t coneSize, Float_t ptThreshold);

  //Process conv gamma
  void ProcessConvGamma( const AliAODEvent * const aodEvent );

  //Process calorimeters
  void ProcessCalorimeters( const AliAODEvent * const aodEvent );
  
  //Correlate particle with jets
  void CorrelateWithJets(AliAODPWG4ParticleCorrelation * photon, const TClonesArray * const jets);

  //Is eta - phi distance smaller than conesize ?
  inline Bool_t IsInCone(Float_t dEta, Float_t dPhi, Float_t coneSize) {   
    return ( (dEta*dEta + dPhi*dPhi) < coneSize*coneSize);
  }

  TList       *fOutputList; //! Output list
  TH1F        *fHistPt; //! Pt spectrum
  TH1F        *fHistPtPhos; //! Pt spectrum
  TH1F        *fHistPtEmcal; //! Pt spectrum
  TH1F        *fHistPtJets; //! Pt spectrum
  TH1F        *fHistGammaJets; //!Phi correlations
  TH1F        *fHistGammaJetsIso; //!Phi correlations
  
  
  Float_t fMinPt; //Minimum pt for correlation
  Float_t fConeSize; //cone size for isolation
  Float_t fPtThreshold; //Threshold pt for isolation

  TString     fDeltaAODFileName;//! File where Gamma Conv AOD is located, if not in default AOD

  TClonesArray * fPhotons;

  AliAnalysisTaskGammaJet(const AliAnalysisTaskGammaJet&); // not implemented
  AliAnalysisTaskGammaJet& operator=(const AliAnalysisTaskGammaJet&); // not implemented
  
  ClassDef(AliAnalysisTaskGammaJet, 2); // example of analysis
};

#endif
