
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
const float  kPtBinV0[kN1+1] = {2.,2.5,3.,4.,5.};

class AliAnalysisTaskLambdaOverK0sJets : public AliAnalysisTaskSE {

 public:

  enum V0LoopStep_t { kTriggerCheck=1, kCorrelation=2, kMixedEvent=3 };

  AliAnalysisTaskLambdaOverK0sJets(const char *name = "AliAnalysisTaskLambdaOverK0sJets");
  virtual ~AliAnalysisTaskLambdaOverK0sJets() {}

  // Setter for global variables in the event
  void SetMC(Bool_t isMC=kTRUE) {fIsMC=isMC;} 
  void SetPID(Bool_t usePID=kTRUE) {fUsePID=usePID;} 
  void SetCentrality(Float_t min=0., Float_t max=90.) {fCentMin=min;fCentMax=max;} 
  void SetQA(Bool_t doQA=kFALSE){fDoQA=doQA;}
  void SetTriggerPt(Float_t ptMinTrig=8., Float_t ptMaxTrig=50.) {fTrigPtMin=ptMinTrig;fTrigPtMax=ptMaxTrig;} 
  void SetTriggerEta(Float_t etaMaxTrig=0.8){fTrigEtaMax=etaMaxTrig;} 
  void SetCheckIDTrig(Bool_t checkIDTrig=kFALSE){fCheckIDTrig=checkIDTrig;}
  void SetSeparateInjectedPart(Bool_t doSep=kTRUE) {fSeparateInjPart=doSep;} 

  // Setters for V0 candidate selection
  // TO BE FIXED!!!
  void SetV0Cuts(Float_t *cutsV0){
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
  void SetMinPtDaughter(Float_t minPtDaughter=0.160) {fMinPtDaughter=minPtDaughter;} 
  void SetMaxEtaDaughter(Float_t maxEta=0.8) {fMaxEtaDaughter=maxEta;} 
  void SetMaxDCADaughter(Float_t maxDCA=1.0) {fMaxDCADaughter=maxDCA;} 
  //   2.  V0 candidate
  void SetMaxY(Float_t yMax=0.5) {fYMax=yMax;} 
  void SetDCAToPrimVtx(Float_t dcaToPrimVtx=0.1) {fDCAToPrimVtx=dcaToPrimVtx;}
  void SetMinCPA(Float_t minCPA=0.998) {fMinCPA=minCPA;} 
  void SetNSigmaPID(Float_t nSigma=3) {fNSigma=nSigma;} 
  void SetCtau(Float_t minCtau = 0., Float_t maxCtau = 3.) {fMinCtau=minCtau;fMaxCtau=maxCtau;} 

  // Getters
  Float_t GetMinCentr() { return fCentMin; }
  Float_t GetMaxCentr() { return fCentMax; }

  // Main functions
  virtual void     UserCreateOutputObjects();
  virtual Bool_t   AcceptV0(AliAODVertex *vtx, const AliAODv0 *v0);
  virtual void     RecCascade(AliAODTrack *trk1,const AliAODTrack *trk2,const AliAODTrack *trkBch,TString histo);
  virtual void     V0Loop(V0LoopStep_t step, Bool_t isTriggered);
  virtual void     TriggerParticle();
    
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);  

 private: 

  AliAnalysisTaskLambdaOverK0sJets(const AliAnalysisTaskLambdaOverK0sJets&);           //not implemented
  AliAnalysisTaskLambdaOverK0sJets& operator=(const AliAnalysisTaskLambdaOverK0sJets&);//not implemented 

  AliAODEvent *fAOD;
  Bool_t   fIsMC;                        //  Use MC data 
  Bool_t   fUsePID;                      //  Use PID for tracks
  Float_t  fCentMin;                     //  Minimum centrality
  Float_t  fCentMax;                     //  Maximum centrality
  Bool_t   fDoQA;                        //  Do Auality Assurance?
  Float_t  fTrigPtMin;                   //  Minimum pt for trigger particle
  Float_t  fTrigPtMax;                   //  Maximum pt for trigger particle
  Float_t  fTrigEtaMax;                  //  Maximum eta for trigger particle
  Bool_t   fCheckIDTrig;                 //  Do comparison with V0's daughter tracks?
  Bool_t   fSeparateInjPart;             //  Separate MC injected particles in case of correlation 
  Int_t    fEndOfHijingEvent;            //  Limit natural-injected MC  particles 
  AliPIDResponse *fPIDResponse;          //  PID Response


  Float_t fMinPtDaughter;                //  Minimum transverse momentum for V0's daughters
  Float_t fMaxEtaDaughter;               //  Maximum pseudo-rapidity for V0's daughters  
  Float_t fMaxDCADaughter;               //  Maximum Distance of Closest Approach between daughters (given in sigmas)
  Float_t fYMax;                         //  Maximum rapidity for V0
  Float_t fDCAToPrimVtx;                 //  Mimimum distance of closest approach of daughters to the vertex            
  Float_t fMinCPA;                       //  Minimum Cosine of the Pointing Angle to the vertex for V0  
  Float_t fNSigma;                       //  Number of sigmas for PID wi dE/dx
  Float_t fMinCtau;                      //  Minimum ctau
  Float_t fMaxCtau;                      //  Maximum ctau

  Int_t   fIdTrigger;                    //  ID track of the trigger particle
  Int_t   fIsTrigFromV0daug;             //  Flag: trigger particle is a V0's daughter
  Int_t   fIsV0LP;                       //  Flag: V0 has the highest pt in the event
  Float_t fPtV0LP;                       //  Pt of the leading V0
  Int_t   fIsSndCheck;                   //  Flag: trigger particle is the second leaidng particle
  

  TList*  fOutput;                       //! List of histograms for main analysis
  TList*  fOutputQA;                     //! List of histograms for Quality Assurance

  TH1F*   fEvents;                       //! Counter for the number of events in each step
  TH1F*   fCentrality;                   //! Event centrality per centil
  TH1F*   fPrimaryVertexX;               //! Primary vertex position in X
  TH1F*   fPrimaryVertexY;               //! Primary vertex position in Y
  TH1F*   fPrimaryVertexZ;               //! Primary vertex position in Z
  TH2F*   fCentMult;                     //! Event centrality vs Track multiplicity
  TH2F*   fdEdx;                         //! dEdx
  TH2F*   fdEdxPid;                      //! dEdx with PID

  TH3F*   fTriggerMCPtCent;              //! Trigger particle MC: pt vs centrality
  TH3F*   fTriggerPtCent;                //! Trigger particle: pt vs centrality vs Z vertex
  TH2F*   fTriggerEtaPhi;                //! Trigger particle: eta vs phi
  TH1F*   fCheckTriggerFromV0Daug;       //! Trigger particle: it is a daughter from a V0-candidate
  TH1F*   fTriggerComingFromDaug;        //! Trigger particle: pt when LP is a daughter from a V0-candidate
  TH1F*   fTriggerIsV0;                  //! Trigger particle: the V0 is the highest-pt particle
  TH3F*   fCheckIDTrigPtK0s;             //! Trigger particle: pt comparison between trigger track and K0s daughter track
  TH3F*   fCheckIDTrigPhiK0s;            //! Trigger particle: phi comparison between trigger track and K0s daughter track
  TH3F*   fCheckIDTrigPtLambda;          //! Trigger particle: pt comparison between trigger track and Lambda daughter track
  TH3F*   fCheckIDTrigPhiLambda;         //! Trigger particle: phi comparison between trigger track and Lambda daughter track

  TH1F*   fInjectedParticles;            //! Number of injected particles

  TH1F*   fK0sMCPt;                      //! K0s MC: pt
  TH3F*   fK0sMCPtRap;                   //! K0s MC: pt vs rapidity
  TH3F*   fK0sMCPtPhiEta;                //! K0s MC: pt vs pseudo-rapidity
  TH1F*   fK0sAssocPt;                   //! K0s Assoc: pt
  TH3F*   fK0sAssocPtArm;                //! K0s Assoc: pt vs decay lenght vs centrality
  TH3F*   fK0sAssocPtRap;                //! K0s Assoc: pt vs rapidity
  TH3F*   fK0sAssocPtPhiEta;             //! K0s Assoc: pt vs pseudo-rapidity
  TH3F*   fK0sMCResPhi;                  //! K0s Assoc: phi resolution

  TH1F*   fLambdaMCPt;                   //! Lambda MC: pt
  TH3F*   fLambdaMCPtRap;                //! Lambda MC: pt vs rapidity
  TH3F*   fLambdaMCPtPhiEta;             //! Lambda MC: pt vs pseudo-rapidity
  TH1F*   fLambdaAssocPt;                //! Lambda Assoc: pt
  TH3F*   fLambdaAssocPtArm;             //! Lambda Assoc: pt vs decay lenght vs centrality
  TH3F*   fLambdaAssocPtRap;             //! Lambda Assoc: pt vs rapidity
  TH3F*   fLambdaAssocPtPhiEta;          //! Lambda Assoc: pt vs pseudo-rapidity
  TH3F*   fLambdaMCResPhi;               //! Lambda Assoc: phi resolution

  TH1F*   fAntiLambdaMCPt;               //! AntiLambda MC: pt
  TH3F*   fAntiLambdaMCPtRap;            //! AntiLambda MC: pt vs rapidity
  TH3F*   fAntiLambdaMCPtPhiEta;         //! AntiLambda MC: pt vs pseudo-rapidity
  TH1F*   fAntiLambdaAssocPt;            //! AntiLambda Assoc: pt
  TH3F*   fAntiLambdaAssocPtArm;         //! AntiLambda Assoc: pt vs decay lenght vs centrality
  TH3F*   fAntiLambdaAssocPtRap;         //! AntiLambda Assoc: pt vs rapidity
  TH3F*   fAntiLambdaAssocPtPhiEta;      //! AntiLambda Assoc: pt vs pseudo-rapidity
  TH3F*   fAntiLambdaMCResPhi;           //! AntiLambda Assoc: phi resolution

  /// ====== Histohgrmas for Correlations ====== ///

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

  TH3F*   fAntiLambdaMass;                     //! Mass for AntiLambda
  TH2F*   fAntiLambdaPtLtSB;                   //! AntiLambda: l vs p with side-band subtraction
  TH3F*   fAntiLambdaPtvsEta;                  //! AntiLambda: pt vs eta
  TH3F*   fAntiLambdaPtvsRap;                  //! AntiLambda: pt vs rap
  TH2F*   fAntiLambdaEtaPhi;                   //! AntiLambda: eta vs phi
  TH3F*   fAntiLambdaMassPtPhi;                //! Lambda: mass vs phi 

  TH3F*   fAntiLambdaMassPtvsPtL;              //! AntiLambda: mass, pt vs pt of leading particle
  TH3F*   fAntiLambdaSiPtL;                    //! AntiLambda: mass, vs leading particle
  TH2F*   fAntiLambdaDaughtersPt;              //! AntiLambda: pt of daughters
  TH3F*   fAntiLambdadPhiPtAssocPtL;           //! AntiLambda: Delta phi,pt vs pt of the leading particle
  TH3F*   fAntiLambdaDCADaugToPrimVtx;         //! AntiLambda: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak

  TH3F*   fAntiLambdadPhidEtaMC[kN1];          //! AntiLambda MC: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fAntiLambdadPhidEtaMCCent[kN1];      //! AntiLambda MC in central events: Delta phi,Delta eta vs pt of the leading particle

  TH3F*   fAntiLambdadPhidEtaPtL[kN1];         //! AntiLambda: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fAntiLambdadPhidEtaPtLCent[kN1];     //! AntiLambda in central events: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fAntiLambdadPhidEtaPtLBckg[kN1];     //! AntiLambda background: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fAntiLambdadPhidEtaPtLCentBckg[kN1]; //! AntiLambda background in central events: Delta phi,Delta eta vs pt of the leading particle

  TH3F*   fAntiLambdadPhidEtaPtL2[kN1];        //! AntiLambda: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fAntiLambdadPhidEtaPtLCent2[kN1];    //! AntiLambda in central events: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fAntiLambdadPhidEtaPtLBckg2[kN1];    //! AntiLambda background: Delta phi,Delta eta vs pt of the leading particle
  TH3F*   fAntiLambdadPhidEtaPtLCentBckg2[kN1];//! AntiLambda background in central events: Delta phi,Delta eta vs pt of the leading particle

  TH2F*   fAntiLambdaBckgDecLength;            //! AntiLambda background: Decay lenght vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fAntiLambdaBckgDCADaugToPrimVtx;     //! AntiLambda background: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH2F*   fAntiLambdadEdxPosDaug;              //! AntiLambda background: dE/dx of the positive daughter particle inside a radio wrt the near-side peak
  TH2F*   fAntiLambdadEdxNegDaug;              //! AntiLambda background: dE/dx of the negative daughter particle inside a radio wrt the near-side peak
  TH2F*   fAntiLambdaBckgEtaPhi;               //! AntiLambda background: Phi vs Eta inside a radio wrt the near-side peak
  TH2F*   fAntiLambdaBckgPhiRadio ;            //! AntiLambda background: Phi vs radio inside a radio wrt the near-side peak
  TH2F*   fAntiLambdaBckgDCANegDaugToPrimVtx;  //! AntiLambda background: DCA of Negative daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  TH2F*   fAntiLambdaBckgDCAPosDaugToPrimVtx;  //! AntiLambda background: DCA of Positive daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  TH2F*   fAntiLambdaMassCascade;              //! AntiLambda background: Poddible mismatching of tracks due to cascades decays

    
  ///  ==== Quality Assurance plots === ///

  //           K0s            //

  TH3F*   fK0sPIDPosDaug;                    //! K0s: Pos. track PID
  TH3F*   fK0sPIDNegDaug;                    //! K0s: Neg. track PID
  TH3F*   fK0sBckgPIDPosDaug;                //! K0s Bckg: Pos. track PID
  TH3F*   fK0sBckgPIDNegDaug;                //! K0s Bckg: Neg. track PID

  TH3F*   fK0sPhiEtaPosDaug;                 //! K0s: Pos. track phi vs eta 
  TH3F*   fK0sPhiEtaNegDaug;                 //! K0s: Neg. track phi vs eta
  TH3F*   fK0sBckgPhiEtaPosDaug;             //! K0s Bckg: Pos. track phi vs eta  
  TH3F*   fK0sBckgPhiEtaNegDaug;             //! K0s Bckg: Neg. track phi vs eta

  TH2F*   fK0sDCAPosDaug;                    //! K0s: Pos. track DCA to primary vertex
  TH2F*   fK0sDCANegDaug;                    //! K0s: Neg. track DCA to primary vertex
  TH2F*   fK0sBckgDCAPosDaug;                //! K0s Bckg: Pos. track DCA to primary vertex
  TH2F*   fK0sBckgDCANegDaug;                //! K0s Bckg: Neg. track DCA to primary vertex

  TH2F*   fK0sDifPtPosDaug;                  //! K0s: Pos. track diference berween pt of the daughter and the V0
  TH2F*   fK0sDifPtNegDaug;                  //! K0s: Neg. track diference berween pt of the daughter and the V0
  TH2F*   fK0sBckgDifPtPosDaug;              //! K0s Bckg: Pos. track diference berween pt of the daughter and the V0
  TH2F*   fK0sBckgDifPtNegDaug;              //! K0s Bckg: Neg. track diference berween pt of the daughter and the V0

  TH3F*   fK0sDecayPos;                      //! K0s: 2D decay position  
  TH3F*   fK0sBckgDecayPos;                  //! K0s Bckg: 2D decay position   
  TH2F*   fK0sDecayVertex;                   //! K0s: decay lenght
  TH2F*   fK0sBckgDecayVertex;               //! K0s Bckg: decay lenght 
  TH2F*   fK0sDecayVertexZoom;               //! K0s: decay lenght Zoom
  TH2F*   fK0sBckgDecayVertexZoom;           //! K0s Bckg: decay lenght Zoom

  TH2F*   fK0sCPA;                           //! K0s: cosine of the pointing angle
  TH2F*   fK0sBckgCPA;                       //! K0s Bckg: cosine of the pointing angle
  TH2F*   fK0sDCAV0Daug;                     //! K0s: distance of the closest approach to the primary vertex
  TH2F*   fK0sBckgDCAV0Daug;                 //! K0s Bckg: distance of the closest approach to the primary vertex 

  TH3F*   fK0sNClustersTPC;                  //! K0s: Numbers of TPC clusters of the daughter tracks 
  TH3F*   fK0sBckgNClustersTPC;              //! K0s Bckg: Numbers of TPC clusters of the daughter tracks 
  TH3F*   fK0sNClustersITSPos;               //! K0s: Pos. Daug. Numbers of ITS clusters of the daughter tracks 
  TH3F*   fK0sNClustersITSNeg;               //! K0s: Neg. Daug. Numbers of ITS clusters of the daughter tracks 
  TH3F*   fK0sBckgNClustersITSPos;           //! K0s Bckg: Pos. Daug. Numbers of ITS clusters of the daughter tracks 
  TH3F*   fK0sBckgNClustersITSNeg;           //! K0s Bckg: Neg. Daug. Numbers of ITS clusters of the daughter tracks 

  //          Lambda          //

  TH3F*   fLambdaPIDPosDaug;                    //! Lambda: Pos. track PID
  TH3F*   fLambdaPIDNegDaug;                    //! Lambda: Neg. track PID
  TH3F*   fLambdaBckgPIDPosDaug;                //! Lambda Bckg: Pos. track PID
  TH3F*   fLambdaBckgPIDNegDaug;                //! Lambda Bckg: Neg. track PID

  TH3F*   fLambdaPhiEtaPosDaug;                 //! Lambda: Pos. track phi vs eta 
  TH3F*   fLambdaPhiEtaNegDaug;                 //! Lambda: Neg. track phi vs eta
  TH3F*   fLambdaBckgPhiEtaPosDaug;             //! Lambda Bckg: Pos. track phi vs eta  
  TH3F*   fLambdaBckgPhiEtaNegDaug;             //! Lambda Bckg: Neg. track phi vs eta

  TH2F*   fLambdaDCAPosDaug;                    //! Lambda: Pos. track DCA to primary vertex
  TH2F*   fLambdaDCANegDaug;                    //! Lambda: Neg. track DCA to primary vertex
  TH2F*   fLambdaBckgDCAPosDaug;                //! Lambda Bckg: Pos. track DCA to primary vertex
  TH2F*   fLambdaBckgDCANegDaug;                //! Lambda Bckg: Neg. track DCA to primary vertex

  TH2F*   fLambdaDifPtPosDaug;                  //! Lambda: Pos. track diference berween pt of the daughter and the V0
  TH2F*   fLambdaDifPtNegDaug;                  //! Lambda: Neg. track diference berween pt of the daughter and the V0
  TH2F*   fLambdaBckgDifPtPosDaug;              //! Lambda Bckg: Pos. track diference berween pt of the daughter and the V0
  TH2F*   fLambdaBckgDifPtNegDaug;              //! Lambda Bckg: Neg. track diference berween pt of the daughter and the V0

  TH3F*   fLambdaDecayPos;                      //! Lambda: 2D decay position  
  TH3F*   fLambdaBckgDecayPos;                  //! Lambda Bckg: 2D decay position   
  TH2F*   fLambdaDecayVertex;                   //! Lambda: decay lenght
  TH2F*   fLambdaBckgDecayVertex;               //! Lambda Bckg: decay lenght 
  TH2F*   fLambdaDecayVertexZoom;               //! Lambda: decay lenght Zoom
  TH2F*   fLambdaBckgDecayVertexZoom;           //! Lambda Bckg: decay lenght Zoom

  TH2F*   fLambdaCPA;                           //! Lambda: cosine of the pointing angle
  TH2F*   fLambdaBckgCPA;                       //! Lambda Bckg: cosine of the pointing angle
  TH2F*   fLambdaDCAV0Daug;                     //! Lambda: distance of the closest approach to the primary vertex
  TH2F*   fLambdaBckgDCAV0Daug;                 //! Lambda Bckg: distance of the closest approach to the primary vertex 

  TH3F*   fLambdaNClustersTPC;                  //! Lambda: Numbers of TPC clusters of the daughter tracks 
  TH3F*   fLambdaBckgNClustersTPC;              //! Lambda Bckg: Numbers of TPC clusters of the daughter tracks 
  TH3F*   fLambdaNClustersITSPos;               //! Lambda: Pos. Daug. Numbers of ITS clusters of the daughter tracks 
  TH3F*   fLambdaNClustersITSNeg;               //! Lambda: Neg. Daug. Numbers of ITS clusters of the daughter tracks 
  TH3F*   fLambdaBckgNClustersITSPos;           //! Lambda Bckg: Pos. Daug. Numbers of ITS clusters of the daughter tracks 
  TH3F*   fLambdaBckgNClustersITSNeg;           //! Lambda Bckg: Neg. Daug. Numbers of ITS clusters of the daughter tracks 

  //        AntiLambda        //

  TH3F*   fAntiLambdaPIDPosDaug;                    //! AntiLambda: Pos. track PID
  TH3F*   fAntiLambdaPIDNegDaug;                    //! AntiLambda: Neg. track PID
  TH3F*   fAntiLambdaBckgPIDPosDaug;                //! AntiLambda Bckg: Pos. track PID
  TH3F*   fAntiLambdaBckgPIDNegDaug;                //! AntiLambda Bckg: Neg. track PID

  TH3F*   fAntiLambdaPhiEtaPosDaug;                 //! AntiLambda: Pos. track phi vs eta 
  TH3F*   fAntiLambdaPhiEtaNegDaug;                 //! AntiLambda: Neg. track phi vs eta
  TH3F*   fAntiLambdaBckgPhiEtaPosDaug;             //! AntiLambda Bckg: Pos. track phi vs eta  
  TH3F*   fAntiLambdaBckgPhiEtaNegDaug;             //! AntiLambda Bckg: Neg. track phi vs eta

  TH2F*   fAntiLambdaDCAPosDaug;                    //! AntiLambda: Pos. track DCA to primary vertex
  TH2F*   fAntiLambdaDCANegDaug;                    //! AntiLambda: Neg. track DCA to primary vertex
  TH2F*   fAntiLambdaBckgDCAPosDaug;                //! AntiLambda Bckg: Pos. track DCA to primary vertex
  TH2F*   fAntiLambdaBckgDCANegDaug;                //! AntiLambda Bckg: Neg. track DCA to primary vertex

  TH2F*   fAntiLambdaDifPtPosDaug;                  //! AntiLambda: Pos. track diference berween pt of the daughter and the V0
  TH2F*   fAntiLambdaDifPtNegDaug;                  //! AntiLambda: Neg. track diference berween pt of the daughter and the V0
  TH2F*   fAntiLambdaBckgDifPtPosDaug;              //! AntiLambda Bckg: Pos. track diference berween pt of the daughter and the V0
  TH2F*   fAntiLambdaBckgDifPtNegDaug;              //! AntiLambda Bckg: Neg. track diference berween pt of the daughter and the V0

  TH3F*   fAntiLambdaDecayPos;                      //! AntiLambda: 2D decay position  
  TH3F*   fAntiLambdaBckgDecayPos;                  //! AntiLambda Bckg: 2D decay position   
  TH2F*   fAntiLambdaDecayVertex;                   //! AntiLambda: decay lenght
  TH2F*   fAntiLambdaBckgDecayVertex;               //! AntiLambda Bckg: decay lenght 
  TH2F*   fAntiLambdaDecayVertexZoom;               //! AntiLambda: decay lenght Zoom
  TH2F*   fAntiLambdaBckgDecayVertexZoom;           //! AntiLambda Bckg: decay lenght Zoom

  TH2F*   fAntiLambdaCPA;                           //! AntiLambda: cosine of the pointing angle
  TH2F*   fAntiLambdaBckgCPA;                       //! AntiLambda Bckg: cosine of the pointing angle
  TH2F*   fAntiLambdaDCAV0Daug;                     //! AntiLambda: distance of the closest approach to the primary vertex
  TH2F*   fAntiLambdaBckgDCAV0Daug;                 //! AntiLambda Bckg: distance of the closest approach to the primary vertex 

  TH3F*   fAntiLambdaNClustersTPC;                  //! AntiLambda: Numbers of TPC clusters of the daughter tracks 
  TH3F*   fAntiLambdaBckgNClustersTPC;              //! AntiLambda Bckg: Numbers of TPC clusters of the daughter tracks 
  TH3F*   fAntiLambdaNClustersITSPos;               //! AntiLambda: Pos. Daug. Numbers of ITS clusters of the daughter tracks 
  TH3F*   fAntiLambdaNClustersITSNeg;               //! AntiLambda: Neg. Daug. Numbers of ITS clusters of the daughter tracks 
  TH3F*   fAntiLambdaBckgNClustersITSPos;           //! AntiLambda Bckg: Pos. Daug. Numbers of ITS clusters of the daughter tracks 
  TH3F*   fAntiLambdaBckgNClustersITSNeg;           //! AntiLambda Bckg: Neg. Daug. Numbers of ITS clusters of the daughter tracks 


  ClassDef(AliAnalysisTaskLambdaOverK0sJets,1);

};

#endif
