#ifndef ALIANALYSISTASKSIGNEDBF_H
#define ALIANALYSISTASKSIGNEDBF_H

//Analysis task for signed BF studies
//Author: Panos.Christakoglou@cern.ch

class TList;
class TH1F;
class TH2F;
class TObjArray;
class TLorentzVector;
class TVector2;
class TVector3;

class AliVEvent;
class AliAnalysisUtils;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSignedBF : public AliAnalysisTaskSE {
 public:  
  AliAnalysisTaskSignedBF(const char *name = "AliAnalysisTaskptspectra");
  virtual ~AliAnalysisTaskSignedBF(); 
   
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);

  //===============Getters================//
  Double_t GetRefMultiOrCentrality(AliVEvent *event);
  TObjArray *GetAcceptedTracks(AliVEvent *event);
  //======================================//

  //===============Setters================//
  void SetAnalysisLevel(const char* gAnalysisLevel = "AOD") {
    fAnalysisLevel = gAnalysisLevel;
  }
  
  void SetMultiplicityEstimator(const char* gMultiplicityEstimator = "V0M") {
    fMultiplicityEstimator = gMultiplicityEstimator;
  }

  void SetAODtrackCutBit(Int_t bit){
    fnAODtrackCutBit = bit;
  }  

  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;
    fVyMax = vy;
    fVzMax = vz;
  }

  void SetUseAdditionalVtxCuts() {
  fUseAdditionalVtxCuts=kTRUE;}

  //percentile
  void SetPercentileRange(Double_t min, Double_t max) { 
    fCentralityPercentileMin = min;
    fCentralityPercentileMax = max;
  }

  void SetKinematicsCutsAOD(Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax){
    fPtMin  = ptmin;  fPtMax  = ptmax;
    fEtaMin = etamin; fEtaMax = etamax;
  }

  void UseOfflineTrigger() {fUseOfflineTrigger = kTRUE;}
  void CheckPileUp() {fCheckPileUp = kTRUE;}
  //======================================//
  
 private:
  Double_t    IsEventAccepted(AliVEvent* event);
  Double_t    GetEventPlane(AliVEvent* event);
  void CalculateSignedBFEbyE(TObjArray *cObjAcceptedParticles, 
			     Double_t gCentrality,
			     Double_t gReactionPlane);

  TList *fListQA; //fList object
  TList *fListBF; //fList object

  //=============Event level=============//
  TH1F *fHistEventStats; //event stats
  TH1F *fHistVx; //vx
  TH1F *fHistVy; //vy
  TH1F *fHistVz; //vz
  TH1F *fHistMultiplicityPercentile; //multiplicity/centrality percentile
  TH1F *fHistMultiplicity;//multiplicity distribution
  TH2F *fHistMultiplicityvsPercentile;//multiplicity vs multiplicity percentile correlation
  TH2F *fHistEventPlane;//event plane angle

  TH1F *fHistNumberOfAcceptedParticles; //accepted particles,
  TH1F *fHistNumberOfAcceptedPositiveParticles; //acepted positive particles
  TH1F *fHistNumberOfAcceptedNegativeParticles; //acepted negative particles
  //=====================================//

  //=============Track level=============//
  TH1F *fHistPt; //pT spectrum
  //=====================================//

  //=============BF=============//
  TH1F *fHistP; //number of positive particle per event
  TH1F *fHistN; //number of negative particle per event

  TH1F *fHistPNLabOut; // +- (lab frame) out-of-plane
  TH1F *fHistNPLabOut; // -+ (lab frame) out-of-plane 
  TH1F *fHistPPLabOut; // ++ (lab frame) out-of-plane
  TH1F *fHistNNLabOut;// -- (lab frame) out-of-plane
  TH1F *fHistPNLabIn; // +- (lab frame) in-plane
  TH1F *fHistNPLabIn; // -+ (lab frame) in-plane
  TH1F *fHistPPLabIn; // ++ (lab frame) in-plane
  TH1F *fHistNNLabIn; // -- (lab frame) in-plane
  TH2F *fHistDeltaBLabOut; //DB (lab frame) out-of-plane 
  TH2F *fHistDeltaBLabIn; //DB (lab frame) in-plane 

  TH1F *fHistPNRestOut; // +- (rest frame) out-of-plane
  TH1F *fHistNPRestOut; // -+ (rest frame) out-of-plane 
  TH1F *fHistPPRestOut; // ++ (rest frame) out-of-plane
  TH1F *fHistNNRestOut;// -- (rest frame) out-of-plane
  TH1F *fHistPNRestIn; // +- (rest frame) in-plane
  TH1F *fHistNPRestIn; // -+ (rest frame) in-plane
  TH1F *fHistPPRestIn; // ++ (rest frame) in-plane
  TH1F *fHistNNRestIn; // -- (rest frame) in-plane
  TH2F *fHistDeltaBRestOut; //DB (rest frame) out-of-plane 
  TH2F *fHistDeltaBRestIn; //DB (rest frame) in-plane 

  TLorentzVector *fLVParticle1; //TLorentzVector particle 1
  TLorentzVector *fLVParticle2; //TLorentzVector particle 2
  TLorentzVector *fLVParticlePair; //TLorentzVector particle pair
  TVector2 *fV2Particle1; //TVector2 particle 1
  TVector2 *fV2Particle2; //TVector2 particle 2
  TVector3 *fV3Particle1; //TVector3 particle 1
  TVector3 *fV3Particle2; //TVector3 particle 2
  //============================//

  //============Other variables==========//
  AliAnalysisUtils *fUtils;//AliAnalysisUtils
  TString fAnalysisLevel; //"AOD" or "MC"
  TString fMultiplicityEstimator; //multiplicity estimator ("V0M")
  Bool_t fUseOfflineTrigger;//use offline trigger
  Double_t fVxMax, fVyMax, fVzMax;//vx vy vz max values
  Bool_t fUseAdditionalVtxCuts;// additional vertex cuts
  Bool_t fCheckPileUp;//pile up checks
  Double_t fCentralityPercentileMin, fCentralityPercentileMax;//centrality ranges
  Double_t fEventPlane; //event plane
  Int_t fnAODtrackCutBit;// AOD filter bit

  Double_t fPtMin;//only used for AODs
  Double_t fPtMax;//only used for AODs
  Double_t fEtaMin;//only used for AODs
  Double_t fEtaMax;//only used for AODs

  AliAnalysisTaskSignedBF(const AliAnalysisTaskSignedBF&); // not implemented
  AliAnalysisTaskSignedBF& operator=(const AliAnalysisTaskSignedBF&); // not implemented
  
  ClassDef(AliAnalysisTaskSignedBF, 1); // example of analysis
};



#endif
