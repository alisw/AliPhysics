
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
 * See cxx source for full Copyright notice */ 

/*
          AliAnalysisTaskLambdaOverK0sJets class

          This program obtains the production of K0s and Lambdas and calculates 
	  the correlation (in the variables phi and eta) with respect to a
	  high-pt charged particle.
	  It works with MC info and AOD tree.
	  WARNING: The Mixed Event part is under construction.
	  Origin: X. Sanchez Castro August2012, xsanchez@cern.ch
*/


#ifndef ALIANALYSISTASKLAMBDAOVERK0SJETS_H
#define ALIANALYSISTASKLAMBDAOVERK0SJETS_H
 
#include "AliAnalysisTaskSE.h"

class AliAODEvent;
class AliPIDResponse;
class AliAODTrack;
class AliAODVertex;
class AliAODv0;

class TH1F;
class TH2F;
class TH3F;
class TList;
class TString;

const int    kN1 = 4; 
const double kPtBinV0[kN1+1] = {2.,2.5,3.,4.,5.};

class AliAnalysisTaskLambdaOverK0sJets : public AliAnalysisTaskSE {

 public:

  enum V0LoopStep_t { kTriggerCheck=1, kCorrelation=2, kMixedEvent=3 };

  AliAnalysisTaskLambdaOverK0sJets(const char *name = "AliAnalysisTaskLambdaOverK0sJets");
  virtual ~AliAnalysisTaskLambdaOverK0sJets() {}

  // Setter for global variables in the event
  void SetMC(Bool_t isMC=kTRUE) {fIsMC=isMC;} 
  void SetPID(Bool_t usePID=kTRUE) {fUsePID=usePID;} 
  void SetCentrality(Double_t min=0., Double_t max=90.) {fCentMin=min;fCentMax=max;} 
  void SetQA(Bool_t doQA=kFALSE){fDoQA=doQA;}
  void SetTriggerPt(Double_t ptMinTrig=8., Double_t ptMaxTrig=50.) {fTrigPtMin=ptMinTrig;fTrigPtMax=ptMaxTrig;} 
  void SetTriggerEta(Double_t etaMaxTrig=0.8){fTrigEtaMax=etaMaxTrig;} 
  void SetCheckIDTrig(Bool_t checkIDTrig=kFALSE){fCheckIDTrig=checkIDTrig;}
  void SetSeparateInjectedPart(Bool_t doSep=kTRUE) {fSeparateInjPart=doSep;} 

  // Setters for V0 candidate selection
  // TO BE FIXED!!!
  void SetV0Cuts(Double_t *cutsV0){
    //   1.  Daughter cuts
    fMinPtDaughter=cutsV0[0];
    fMaxEtaDaughter=cutsV0[1];
    fMaxDCADaughter=cutsV0[2];
    //   2.  V0 candidate
    fYMax=cutsV0[3];
    fDCAToPrimVtx=cutsV0[4];
    fMinCPA=cutsV0[5];
    fNSigma=cutsV0[6];
    fMinCtau=cutsV0[7];
    fMaxCtau=cutsV0[8];
  }

  //   1.  Daughter cuts
  void SetMinPtDaughter(Double_t minPtDaughter=0.160) {fMinPtDaughter=minPtDaughter;} 
  void SetMaxEtaDaughter(Double_t maxEta=0.8) {fMaxEtaDaughter=maxEta;} 
  void SetMaxDCADaughter(Double_t maxDCA=1.0) {fMaxDCADaughter=maxDCA;} 
  //   2.  V0 candidate
  void SetMaxY(Double_t yMax=0.5) {fYMax=yMax;} 
  void SetDCAToPrimVtx(Double_t dcaToPrimVtx=0.1) {fDCAToPrimVtx=dcaToPrimVtx;}
  void SetMinCPA(Double_t minCPA=0.998) {fMinCPA=minCPA;} 
  void SetNSigmaPID(Double_t nSigma=3) {fNSigma=nSigma;} 
  void SetCtau(Double_t minCtau = 0., Double_t maxCtau = 3.) {fMinCtau=minCtau;fMaxCtau=maxCtau;} 

  // Getters
  Double_t GetMinCentr() { return fCentMin; }
  Double_t GetMaxCentr() { return fCentMax; }

  // Main functions
  virtual void     UserCreateOutputObjects();
  virtual Bool_t   AcceptV0(AliAODVertex *vtx, const AliAODv0 *v0);
  virtual void     RecCascade(AliAODTrack *trk1,const AliAODTrack *trk2,const AliAODTrack *trkBch,TString histo);
  virtual TArrayD* V0Loop(AliAODTrack *trkLP, V0LoopStep_t step, Bool_t isTriggered);
  virtual TArrayD* TriggerParticle();
    
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);  

 private: 

  AliAnalysisTaskLambdaOverK0sJets(const AliAnalysisTaskLambdaOverK0sJets&);           //not implemented
  AliAnalysisTaskLambdaOverK0sJets& operator=(const AliAnalysisTaskLambdaOverK0sJets&);//not implemented 

  AliAODEvent *fAOD;
  Bool_t   fIsMC;                        //  Use MC data 
  Bool_t   fUsePID;                      //  Use PID for tracks
  Double_t fCentMin;                     //  Minimum centrality
  Double_t fCentMax;                     //  Maximum centrality
  Bool_t   fDoQA;                        //  Do Auality Assurance?
  Double_t fTrigPtMin;                   //  Minimum pt for trigger particle
  Double_t fTrigPtMax;                   //  Maximum pt for trigger particle
  Double_t fTrigEtaMax;                  //  Maximum eta for trigger particle
  Bool_t   fCheckIDTrig;
  Bool_t   fSeparateInjPart;             //  Separate MC injected particles in case of correlation 
  Int_t    fEndOfHijingEvent;            //  Limit natural-injected MC  particles 
  AliPIDResponse *fPIDResponse;          //  PID Response


  Double_t fMinPtDaughter;               //  Minimum transverse momentum for V0's daughters
  Double_t fMaxEtaDaughter;              //  Maximum pseudo-rapidity for V0's daughters  
  Double_t fMaxDCADaughter;              //  Maximum Distance of Closest Approach between daughters (given in sigmas)
  Double_t fYMax;                        //  Maximum rapidity for V0
  Double_t fDCAToPrimVtx;                //  Mimimum distance of closest approach of daughters to the vertex           
  Double_t fMinCPA;                      //  Minimum Cosine of the Pointing Angle to the vertex for V0  
  Double_t fNSigma;                      //  Number of sigmas for PID wi dE/dx
  Double_t fMinCtau;                     //  Minimum ctau
  Double_t fMaxCtau;                     //  Maximum ctau

  TList*  fOutput;                       //! List of histograms
  TList*  fOutputQA;                     //! List of histograms

  TH1F*   fEvents;                       //! Counter for the number of events in each step

  TH1F*   fCentrality;                   //! Event centrality per centil
  TH1F*   fPrimaryVertexX;               //! Primary vertex position in X
  TH1F*   fPrimaryVertexY;               //! Primary vertex position in Y
  TH1F*   fPrimaryVertexZ;               //! Primary vertex position in Z
  TH2F*   fNumberPileUp;                 //! Number of pile up: SPD vs Tracks 
  TH2F*   fCentMult;                     //! Event centrality vs Track multiplicity
  TH2F*   fdEdx;                         //! dEdx
  TH2F*   fdEdxPid;                      //! dEdx with PID

  TH3F*   fTriggerMCPtCent;              //! Trigger particle MC: pt vs centrality
  TH2F*   fTriggerPtCent;                //! Trigger particle: pt vs centrality
  TH2F*   fTriggerEtaPhi;                //! Trigger particle: eta vs phi
  TH1F*   fCheckTriggerFromV0Daug;       //! Trigger particle: it is a daughter from a V0-candidate
  TH1F*   fTriggerComingFromDaug;        //! Trigger particle: pt when LP is a daughter from a V0-candidate
  TH1F*   fTriggerIsV0;                  //! Trigger particle: the V0 is the highest-pt particle
  TH3F*   fCheckIDTrigPtK0s;
  TH3F*   fCheckIDTrigPhiK0s;
  TH3F*   fCheckIDTrigPtLambda;
  TH3F*   fCheckIDTrigPhiLambda;

  TH1F*   fInjectedParticles;            //! Number of injected particles

  TH1F*   fK0sMCPt;                      //! K0s MC: pt
  TH2F*   fK0sMCPtRap;                   //! K0s MC: pt vs rapidity
  TH2F*   fK0sMCPtEta;                   //! K0s MC: pt vs pseudo-rapidity
  TH3F*   fK0sMCPtLt;                    //! K0s MC: pt vs decay lenght vs centrality
  TH1F*   fK0sAssocPt;                   //! K0s Assoc: pt
  TH3F*   fK0sAssocPtLt;                 //! K0s Assoc: pt vs decay lenght vs centrality
  TH3F*   fK0sAssocPtLtArm;              //! K0s Assoc: pt vs decay lenght vs centrality
  TH2F*   fK0sAssocPtRap;                //! K0s Assoc: pt vs rapidity
  TH2F*   fK0sAssocPtEta;                //! K0s Assoc: pt vs pseudo-rapidity

  TH1F*   fLambdaMCPt;                   //! Lambda MC: pt
  TH2F*   fLambdaMCPtRap;                //! Lambda MC: pt vs rapidity
  TH2F*   fLambdaMCPtEta;                //! Lambda MC: pt vs pseudo-rapidity
  TH3F*   fLambdaMCPtLt;                 //! Lambda MC: pt vs decay lenght vs centrality
  TH1F*   fLambdaAssocPt;                //! Lambda Assoc: pt
  TH3F*   fLambdaAssocPtLt;              //! Lambda Assoc: pt vs decay lenght vs centrality
  TH3F*   fLambdaAssocPtLtArm;           //! Lambda Assoc: pt vs decay lenght vs centrality
  TH2F*   fLambdaAssocPtRap;             //! Lambda Assoc: pt vs rapidity
  TH2F*   fLambdaAssocPtEta;             //! Lambda Assoc: pt vs pseudo-rapidity

  TH3F*   fHistArmenterosPodolanski;     //! Armenteros-Podolanski plot inside 3 sigma of the signal
  TH3F*   fHistArmPodBckg;               //! Armenteros-Podolanski plot outside 3 sigma of the signal      

  TH3F*   fK0sMass;                      //! Mass for K0s
  TH2F*   fK0sPtLtSB;                    //! K0s: Side-band subtracted lt vs pt
  TH3F*   fK0sPtvsEta;                   //! K0s: pt vs eta
  TH3F*   fK0sPtvsRap;                   //! K0s: pt vs rap
  TH2F*   fK0sEtaPhi;                    //! K0s: eta vs phi
  TH3F*   fK0sMassPtPhi;                 //! K0s: mass vs phi

  TH3F*   fK0sMassPtvsPtL;               //! K0s: mass, pt vs pt of leading particle
  TH3F*   fK0sSiPtL;                     //! K0s: mass, vs leading particle
  TH2F*   fK0sDaughtersPt;               //! K0s: pt of daughters
  TH3F*   fK0sdPhiPtAssocPtL;            //! K0s: Delta phi,pt vs pt of the leading particle
  TH3F*   fK0sDCADaugToPrimVtx;          //! K0s: DCA to primary vertex of daughters vs leading particle's pt inside a radio wrt the near-side peak
   
  TH3F*   fK0sdPhidEtaMC[kN1];           //! K0s MC: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fK0sdPhidEtaMCCent[kN1];       //! K0s MC in central events: Delta phi,Delta eta vs pt of the leading particle

  TH3F*   fK0sdPhidEtaPtL[kN1];          //! K0s: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fK0sdPhidEtaPtLCent[kN1];      //! K0s in central events: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fK0sdPhidEtaPtLBckg[kN1];      //! K0s background: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fK0sdPhidEtaPtLCentBckg[kN1];  //! K0s background in central events: Delta phi,Delta eta vs pt of the leading particle

  TH3F*   fK0sdPhidEtaPtL2[kN1];         //! K0s: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fK0sdPhidEtaPtLCent2[kN1];     //! K0s in central events: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fK0sdPhidEtaPtLBckg2[kN1];     //! K0s background: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fK0sdPhidEtaPtLCentBckg2[kN1]; //! K0s background in central events: Delta phi,Delta eta vs pt of the leading particle

  TH2F*   fK0sBckgDecLength;             //! K0s background: Decay lenght vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fK0sBckgDCADaugToPrimVtx;      //! K0s background: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH2F*   fK0sdEdxPosDaug;               //! K0s background: dE/dx of the positive daughter particle inside a radio wrt the near-side peak
  TH2F*   fK0sdEdxNegDaug;               //! K0s background: dE/dx of the negative daughter particle inside a radio wrt the near-side peak
  TH2F*   fK0sBckgEtaPhi;                //! K0s background: Phi vs Eta inside a radio wrt the near-side peak
  TH2F*   fK0sBckgPhiRadio;              //! K0s background: Phi vs radio inside a radio wrt the near-side peak
  TH2F*   fK0sBckgDCANegDaugToPrimVtx;   //! K0s background: DCA of Negative daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  TH2F*   fK0sBckgDCAPosDaugToPrimVtx;   //! K0s background: DCA of Positive daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  TH2F*   fK0sMassCascade;               //! K0s background: Poddible mismatching of tracks due to cascades decays

  TH3F*   fLambdaMass;                   //! Mass for Lambda
  TH2F*   fLambdaPtLtSB;                 //! Lambda: l vs p with side-band subtraction
  TH3F*   fLambdaPtvsEta;                //! Lambda: pt vs eta
  TH3F*   fLambdaPtvsRap;                //! Lambda: pt vs rap
  TH2F*   fLambdaEtaPhi;                 //! Lambda: eta vs phi
  TH3F*   fLambdaMassPtPhi;              //! Lambda: mass vs phi 

  TH2F*   fLambdadEdx;                   //! Lambda: dE/dx for proton 
  TH1F*   fCPA;                          //! Lambda: Cosine of the pointing angle
  TH1F*   fDCA;                          //! Lambda: DCA between daughters

  TH3F*   fLambdaMassPtvsPtL;            //! Lambda: mass, pt vs pt of leading particle
  TH3F*   fLambdaSiPtL;                  //! Lambda: mass, vs leading particle
  TH2F*   fLambdaDaughtersPt;            //! Lambda: pt of daughters
  TH3F*   fLambdadPhiPtAssocPtL;         //! Lambda: Delta phi,pt vs pt of the leading particle
  TH3F*   fLambdaDCADaugToPrimVtx;       //! Lambda: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak

  TH3F*   fLambdadPhidEtaMC[kN1];          //! Lambda MC: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fLambdadPhidEtaMCCent[kN1];      //! Lambda MC in central events: Delta phi,Delta eta vs pt of the leading particle

  TH3F*   fLambdadPhidEtaPtL[kN1];         //! Lambda: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fLambdadPhidEtaPtLCent[kN1];     //! Lambda in central events: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fLambdadPhidEtaPtLBckg[kN1];     //! Lambda background: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fLambdadPhidEtaPtLCentBckg[kN1]; //! Lambda background in central events: Delta phi,Delta eta vs pt of the leading particle

  TH3F*   fLambdadPhidEtaPtL2[kN1];        //! Lambda: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fLambdadPhidEtaPtLCent2[kN1];    //! Lambda in central events: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fLambdadPhidEtaPtLBckg2[kN1];    //! Lambda background: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fLambdadPhidEtaPtLCentBckg2[kN1];//! Lambda background in central events: Delta phi,Delta eta vs pt of the leading particle

  TH2F*   fLambdaBckgDecLength;            //! Lambda background: Decay lenght vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fLambdaBckgDCADaugToPrimVtx;     //! Lambda background: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH2F*   fLambdadEdxPosDaug;              //! Lambda background: dE/dx of the positive daughter particle inside a radio wrt the near-side peak
  TH2F*   fLambdadEdxNegDaug;              //! Lambda background: dE/dx of the negative daughter particle inside a radio wrt the near-side peak
  TH2F*   fLambdaBckgEtaPhi;               //! Lambda background: Phi vs Eta inside a radio wrt the near-side peak
  TH2F*   fLambdaBckgPhiRadio ;            //! Lambda background: Phi vs radio inside a radio wrt the near-side peak
  TH2F*   fLambdaBckgDCANegDaugToPrimVtx;  //! Lambda background: DCA of Negative daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  TH2F*   fLambdaBckgDCAPosDaugToPrimVtx;  //! Lambda background: DCA of Positive daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  TH2F*   fLambdaMassCascade;              //! Lambda background: Poddible mismatching of tracks due to cascades decays
        

  TH3F*   fK0sPIDPosDaug;
  TH3F*   fK0sPIDNegDaug;
  TH3F*   fK0sBckgPIDPosDaug;
  TH3F*   fK0sBckgPIDNegDaug;

  TH3F*   fK0sPhiEtaPosDaug;
  TH3F*   fK0sPhiEtaNegDaug;
  TH3F*   fK0sBckgPhiEtaPosDaug;
  TH3F*   fK0sBckgPhiEtaNegDaug;

  TH2F*   fK0sDCAPosDaug;
  TH2F*   fK0sDCANegDaug;
  TH2F*   fK0sBckgDCAPosDaug;
  TH2F*   fK0sBckgDCANegDaug;

  TH2F*   fK0sDifPtPosDaug;
  TH2F*   fK0sDifPtNegDaug;
  TH2F*   fK0sBckgDifPtPosDaug;
  TH2F*   fK0sBckgDifPtNegDaug;

  TH3F*   fK0sDecayPos;
  TH3F*   fK0sBckgDecayPos;
  TH2F*   fK0sDecayVertex;
  TH2F*   fK0sBckgDecayVertex;
  TH2F*   fK0sDecayVertexZoom;
  TH2F*   fK0sBckgDecayVertexZoom;

  TH2F*   fK0sCPA;
  TH2F*   fK0sBckgCPA;
  TH2F*   fK0sDCAV0Daug;
  TH2F*   fK0sBckgDCAV0Daug;

  TH3F*   fLambdaPIDPosDaug;
  TH3F*   fLambdaPIDNegDaug;
  TH3F*   fLambdaBckgPIDPosDaug;
  TH3F*   fLambdaBckgPIDNegDaug;

  TH3F*   fLambdaPhiEtaPosDaug;
  TH3F*   fLambdaPhiEtaNegDaug;
  TH3F*   fLambdaBckgPhiEtaPosDaug;
  TH3F*   fLambdaBckgPhiEtaNegDaug;

  TH2F*   fLambdaDCAPosDaug;
  TH2F*   fLambdaDCANegDaug;
  TH2F*   fLambdaBckgDCAPosDaug;
  TH2F*   fLambdaBckgDCANegDaug;

  TH2F*   fLambdaDifPtPosDaug;
  TH2F*   fLambdaDifPtNegDaug;
  TH2F*   fLambdaBckgDifPtPosDaug;
  TH2F*   fLambdaBckgDifPtNegDaug;

  TH3F*   fLambdaDecayPos;
  TH3F*   fLambdaBckgDecayPos;
  TH2F*   fLambdaDecayVertex;
  TH2F*   fLambdaBckgDecayVertex;
  TH2F*   fLambdaDecayVertexZoom;
  TH2F*   fLambdaBckgDecayVertexZoom;

  TH2F*   fLambdaCPA;
  TH2F*   fLambdaBckgCPA;
  TH2F*   fLambdaDCAV0Daug;
  TH2F*   fLambdaBckgDCAV0Daug;

  ClassDef(AliAnalysisTaskLambdaOverK0sJets,1);

};

#endif
