
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
 * See cxx source for full Copyright notice */ 

/*
                AliAnalysisTaskLambdaOverK0sJets class

                This program obtains the production of K0s and Lambdas and calculates 
                the correlation (in the variables phi and eta) with respect to the
                triggers particles (high-pt charged particles).
                It works with MC information and AOD tree.
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

class TObjArray;

// pt fpr V0
const int    kN1 = 8; 
const float  kPtBinV0[kN1+1] = {2.0,2.25,2.5,2.75,3.0,3.5,4.0,5.0,7.0};

// pt for charged particles
const int    kNc = 3; 
const float  kPtBinCharged[kN1+1] = {2.0,2.25,3.0,4.0};

// pt bins for Xi minus
const int    kN2 = 12; 
const float  kPtBinV02[kN2+1] = {0.0,2.0,2.25,2.5,2.75,3.0,3.5,4.0,5.0,7.0,10.0,15.0,1000.};

const int    kN3 = 3; 
const float  kPtBinV03[kN3+1] = {0.0,2.0,7.0,1000.};
// -------

const int    kNVtxZ = 10; 
const double kBinVtxZ[kNVtxZ+1] = {-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.};

const int    kNCent  = 4;
const double kBinCent[kNCent+1] = {0.0,5.0,10.0,20.0,40.0};

//  ------------------------------------
//  Inv. Mass width as function of the centrality
//  Linear polimomial dependence: sigma(pt) = a0 * a1*pt

const double kCteK0s2010[kNCent] = {0.00367, 0.00363, 0.00358, 0.00348};
const double kLinearK0s2010[kNCent] = {6.148E-4, 5.937E-4, 5.741E-4, 5.693E-4};

const double kCteK0s2011[kNCent] = {0.00354, 0.00348, 0.00360, 0.00352};
const double kLinearK0s2011[kNCent] = {6.526E-4, 6.497E-4, 5.853E-4, 5.808E-4};

const double kCteLambda2010[kNCent] = {0.00113, 0.00114, 0.00119, 0.00119};
const double kLinearLambda2010[kNCent] = {3.062E-4, 2.900E-4, 2.629E-4, 2.440E-4};

const double kCteLambda2011[kNCent] = {9.81E-4, 9.212E-4, 9.876E-4, 0.00106};
const double kLinearLambda2011[kNCent] = {3.878E-4, 3.965E-4, 3.611E-4 , 3.351E-4};

const double kCteAntiLambda2010[kNCent] = {0.00109, 0.00134, 0.00117, 0.00116};
const double kLinearAntiLambda2010[kNCent] = {3.245E-4, 2.308E-4, 2.707E-4, 2.562E-4};

const double kCteAntiLambda2011[kNCent] = {9.859E-4, 0.00111, 0.00104, 0.00110};
const double kLinearAntiLambda2011[kNCent] = {3.881E-4, 3.379E-4, 3.490E-4, 3.166E-4};

// -------------------------------------

class AliAnalysisTaskLambdaOverK0sJets : public AliAnalysisTaskSE {

 public:
  
  enum V0LoopStep_t { kTriggerCheck=1, kReconstruction=2 };

  AliAnalysisTaskLambdaOverK0sJets(const char *name = "AliAnalysisTaskLambdaOverK0sJets");
  virtual ~AliAnalysisTaskLambdaOverK0sJets();

  // Setter for global variables in the event
  void SetCollisionType(TString data="PbPb2010") {fCollision=data;}
  void SetMC(Bool_t isMC=kTRUE) {fIsMC=isMC;} 
  void SetPID(Bool_t usePID=kTRUE) {fUsePID=usePID;} 
  void SetCentrality(Float_t min=0., Float_t max=90.) {fCentMin=min;fCentMax=max;} 
  void SetQA(Bool_t doQA=kFALSE){fDoQA=doQA;}
  void SetDoMix(Bool_t doMixEvt=kTRUE) {fDoMixEvt=doMixEvt;} 
  void SetTriggerPt(Float_t ptMinTrig=8., Float_t ptMaxTrig=50.) {fTrigPtMin=ptMinTrig;fTrigPtMax=ptMaxTrig;} 
  void SetTriggerEta(Float_t etaMaxTrig=0.8){fTrigEtaMax=etaMaxTrig;} 
  void SetCheckIDTrig(Bool_t checkIDTrig=kFALSE){fCheckIDTrig=checkIDTrig;}
  void SetSeparateInjectedPart(Bool_t doSep=kTRUE) {fSeparateInjPart=doSep;} 

  //   1.  Daughter cuts
  void SetMinPtDaughter(Float_t minPtDaughter=0.160) {fMinPtDaughter=minPtDaughter;} 
  void SetMaxEtaDaughter(Float_t maxEta=0.8) {fMaxEtaDaughter=maxEta;} 
  void SetMaxDCADaughter(Float_t maxDCA=1.0) {fMaxDCADaughter=maxDCA;} 
  void SetDCAToPrimVtx(Float_t dcaToPrimVtx=0.1) {fDCAToPrimVtx=dcaToPrimVtx;}
  void SetNSigmaPID(Float_t nSigma=3) {fNSigma=nSigma;} 
  void SetNClsTPC(Float_t nClsTPC=70.) {fDaugNClsTPC=nClsTPC;}
  //   2.  V0 candidate
  void SetMaxY(Float_t yMax=0.5) {fYMax=yMax;} 
  void SetMinCPA(Float_t minCPA=0.998) {fMinCPA=minCPA;} 
  void SetCtau(Float_t minCtau = 0., Float_t maxCtau = 3.) {fMinCtau=minCtau;fMaxCtau=maxCtau;} 

  // Getters
  Float_t GetMinCentr() { return fCentMin; }
  Float_t GetMaxCentr() { return fCentMax; }

  // Main functions
  virtual void     UserCreateOutputObjects();
  virtual Bool_t   AcceptTrack(AliAODTrack *t); 
  virtual Bool_t   AcceptTrackV0(const AliAODTrack *t);
  virtual Bool_t   AcceptV0(AliAODVertex *vtx, const AliAODv0 *v0);
  virtual void     RecCascade(AliAODTrack *trk1,const AliAODTrack *trk2,const AliAODTrack *trkBch,TString histo);
  virtual void     V0Loop(V0LoopStep_t step, Bool_t isTriggered, Int_t iArray, Int_t idTrig);
  virtual void     TriggerParticle();
    
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);  

 private: 

  AliAnalysisTaskLambdaOverK0sJets(const AliAnalysisTaskLambdaOverK0sJets&);           //not implemented
  AliAnalysisTaskLambdaOverK0sJets& operator=(const AliAnalysisTaskLambdaOverK0sJets&);//not implemented 

  AliAODEvent *fAOD;
  TString  fCollision;                   //  Data: PbPb2010 / PbPb2011
  Bool_t   fIsMC;                        //  Use MC data 
  Bool_t   fUsePID;                      //  Use PID for tracks
  Float_t  fCentMin;                     //  Minimum centrality
  Float_t  fCentMax;                     //  Maximum centrality
  Bool_t   fDoQA;                        //  Do Auality Assurance?
  Bool_t   fDoMixEvt;                    //  Do Mixed Events
  Float_t  fTrigPtMin;                   //  Minimum pt for trigger particle
  Float_t  fTrigPtMax;                   //  Maximum pt for trigger particle
  Float_t  fTrigPtMCMin;                 //  Minimum pt for trigger particle in MC
  Float_t  fTrigPtMCMax;                 //  Maximum pt for trigger particle in MC
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
  Float_t fDaugNClsTPC;                  //  Number of TPC clusters for daughters
  Float_t fMinCtau;                      //  Minimum ctau
  Float_t fMaxCtau;                      //  Maximum ctau

  Int_t   fIdTrigger;                    //  ID track of the trigger particle
  Int_t   fIsV0LP;                       //  Flag: V0 has the highest pt in the event
  Float_t fPtV0LP;                       //  Pt of the leading V0
  Int_t   fIsSndCheck;                   //  Flag: trigger particle is the second leaidng particle

  TList*  fOutput;                       //! List of histograms for main analysis
  TList*  fOutputQA;                     //! List of histograms for Quality Assurance
  TList*  fOutputME;                     //! List of histograms for Mixed Events
  TList** fMEList;                       //![] List of Mixed Events

  TObjArray* fTriggerParticles;          // Trigger particle array
  TObjArray* fChargedAssocParticles;      // Trigger particle array
  TObjArray* fTriggerPartMC;             // MC Trigger particle array
  TObjArray* fAssocParticles;            // Associated particle array
  TObjArray* fAssocPartMC;               // MC Associated particle array
  TObjArray* fXiTriggerPartMC;           // Xi leading particle: MC Trigger particle array

  TH1F*   fEvents;                       //! Counter for the number of events in each step
  TH1F*   fCentrality;                   //! Event centrality per centil
  TH1F*   fCentrality2;                  //! Event centrality per centil with |VtxZ|<10cm
  TH2F*   fCentralityTrig;               //! Event centrality per trigger
  TH1F*   fPrimaryVertexX;               //! Primary vertex position in X
  TH1F*   fPrimaryVertexY;               //! Primary vertex position in Y
  TH1F*   fPrimaryVertexZ;               //! Primary vertex position in Z

  TH1F*   fTriggerEventPlane;            //! Distance between the trigger particle direction and the event plane angle

  TH2F*   fTriggerMCPtCent;              //! Trigger particle MC: pt vs centrality
  TH3F*   fTriggerMCResPt;               //! Trigger particle MC: pt resolution
  TH3F*   fTriggerMCResEta;              //! Trigger particle MC: eta resolution
  TH3F*   fTriggerMCResPhi;              //! Trigger particle MC: phi resolution
  TH3F*   fTriggerPtCent;                //! Trigger particle: pt vs centrality vs Z vertex
  TH2F*   fNTrigPerEvt;                  //! Trigger particle: Number of particle triggers per event
  TH1F*   fTriggerWiSPDHit;              //! Trigger particle: Has Hits in the SPD?
  TH2F*   fTriggerEtaPhi;                //! Trigger particle: eta vs phi
  TH1F*   fCheckTriggerFromV0Daug;       //! Trigger particle: it is a daughter from a V0-candidate
  TH1F*   fTriggerComingFromDaug;        //! Trigger particle: pt when LP is a daughter from a V0-candidate
  TH1F*   fTriggerIsV0;                  //! Trigger particle: the V0 is the highest-pt particle
  TH3F*   fCheckIDTrigPtK0s;             //! Trigger particle: pt comparison between trigger track and K0s daughter track
  TH3F*   fCheckIDTrigPhiK0s;            //! Trigger particle: phi comparison between trigger track and K0s daughter track
  TH3F*   fCheckIDTrigEtaK0s;            //! Trigger particle: eta comparison between trigger track and K0s daughter track
  TH3F*   fCheckIDTrigNclsK0s;           //! Trigger particle: number of cluster of the daughter particle 
  TH3F*   fCheckIDTrigPtLambda;          //! Trigger particle: pt comparison between trigger track and Lambda daughter track
  TH3F*   fCheckIDTrigPhiLambda;         //! Trigger particle: phi comparison between trigger track and Lambda daughter track
  TH3F*   fCheckIDTrigEtaLambda;         //! Trigger particle: eta comparison between trigger track and Lambda daughter track
  TH3F*   fCheckIDTrigNclsLambda;        //! Trigger particle: number of cluster of the daughter particle  
  TH3F*   fCheckIDTrigPtAntiLambda;      //! Trigger particle: pt comparison between trigger track and AntiLambda daughter track
  TH3F*   fCheckIDTrigPhiAntiLambda;     //! Trigger particle: phi comparison between trigger track and AntiLambda daughter track
  TH3F*   fCheckIDTrigEtaAntiLambda;     //! Trigger particle: eta comparison between trigger track and AntiLambda daughter track
  TH3F*   fCheckIDTrigNclsAntiLambda;    //! Trigger particle: number of cluster of the daughter particle 

  // ==============  Monte Carlo  ================= //
  TH1F*   fInjectedParticles;            //! Number of injected particles

  //           K0s            //
  TH1F*   fK0sMCPt;                      //! K0s MC: pt
  TH3F*   fK0sMCPtRap;                   //! K0s MC: pt vs rapidity
  TH3F*   fK0sMCPtRap2;                  //! K0s MC: pt vs rapidity (is Natural)
  TH3F*   fK0sMCPtRapVtx;                //! K0s MC: pt vs Z vtx position  vs centrality
  TH3F*   fK0sMCPtRapEmbeded;            //! K0s MC: pt vs rapidity  (embeded particles)
  TH3F*   fK0sMCPtRapVtxEmbeded;         //! K0s MC: pt vs Z vtx position rapidity  vs centrality (embeded particles)
  TH3F*   fK0sMCPtPhiEta[kNCent];        //! K0s MC: pt vs pseudo-rapidity

  TH1F*   fK0sAssocPt;                         //! K0s Assoc: pt
  TH3F*   fK0sAssocPtArm;                      //! K0s Assoc: pt vs rapidity vs centrality (arm. pod. cut)
  TH3F*   fK0sAssocPtRap;                      //! K0s Assoc: pt vs rapidity vs centrality
  TH3F*   fK0sAssocPtPhiEta[kNCent];           //! K0s Assoc: pt vs pseudo-rapidity

  TH3F*   fK0sAssocPtMassArm[kNCent];          //! K0s Assoc: mass vs pt vs centrality
  TH3F*   fK0sAssocMassPtVtx[kNCent];          //! K0s Assoc: mass vs pt vs Z vertex position
  TH3F*   fK0sAssocMassPtDCADaug[kNCent];      //! K0s Assoc: mass vs pt vs dca between daughters
  TH3F*   fK0sAssocMassPtCPA[kNCent];          //! K0s Assoc: mass vs pt vs cpa
  TH3F*   fK0sAssocMassPtDCAPV[kNCent];        //! K0s Assoc: mass vs pt vs dca to prim. vtx
  TH3F*   fK0sAssocMassPtDaugNClsTPC[kNCent];  //! K0s Assoc: mass vs pt vs num. of tpc clusters

  TH3F*   fK0sAssocPtRapEmbeded;                      //! K0s Assoc: pt vs rapidity vs centrality  (embeded particles)
  TH3F*   fK0sAssocPtMassArmEmbeded[kNCent];          //! K0s Assoc: mass vs pt vs rapidity  (embeded particles)
  TH3F*   fK0sAssocMassPtVtxEmbeded[kNCent];          //! K0s Assoc: mass vs pt vs Z vertex position  (embeded particles)
  TH3F*   fK0sAssocMassPtDCADaugEmbeded[kNCent];      //! K0s Assoc: mass vs pt vs dca between daughters  (embeded particles)
  TH3F*   fK0sAssocMassPtCPAEmbeded[kNCent];          //! K0s Assoc: mass vs pt vs cpa  (embeded particles)
  TH3F*   fK0sAssocMassPtDCAPVEmbeded[kNCent];        //! K0s Assoc: mass vs pt vs dca to prim. vtx  (embeded particles)
  TH3F*   fK0sAssocMassPtDaugNClsTPCEmbeded[kNCent];  //! K0s Assoc: mass vs pt vs num. o ftpc clusters (embeded particles)

  TH3F*   fK0sMCResEta;                  //! K0s Assoc: eta resolution
  TH3F*   fK0sMCResPhi;                  //! K0s Assoc: phi resolution


  //           Lambda            //
  TH1F*   fLambdaMCPt;                   //! Lambda MC: pt
  TH3F*   fLambdaMCPtRap;                //! Lambda MC: pt vs rapidity
  TH3F*   fLambdaMCPtRap2;               //! Lambda MC: pt vs rapidity  (is Natural)
  TH3F*   fLambdaMCPtRapVtx;             //! Lambda MC: pt vs Z vtx position rapidity vs centrality
  TH3F*   fLambdaMCPtRapEmbeded;         //! Lambda MC: pt vs rapidity (embeded particles)
  TH3F*   fLambdaMCPtRapVtxEmbeded;      //! Lambda MC: pt vs Z vtx position vs centrality  (embeded particles)
  TH2F*   fLambdaMCFromXi;               //! Lambda MC: coming from Xi
  TH3F*   fLambdaMCPtPhiEta[kNCent];     //! Lambda MC: pt vs pseudo-rapidity

  TH1F*   fLambdaAssocPt;                //! Lambda Assoc: pt
  TH3F*   fLambdaAssocPtRap;             //! Lambda Assoc: pt vs rapidity
  TH2F*   fLambdaAssocFromXi;            //! Lambda Assoc: coming from Xi
  TH3F*   fLambdaAssocPtPhiEta[kNCent];  //! Lambda Assoc: pt vs pseudo-rapidity

  TH3F*   fLambdaAssocMassPtRap[kNCent];          //! Lambda Assoc: pt vs rapidity vs mass
  TH3F*   fLambdaAssocMassPtRap2[kNCent];         //! Lambda Assoc: pt vs rapidity vs mass (wo Cross contamination)
  TH3F*   fLambdaAssocMassPtVtx[kNCent];          //! Lambda Assoc: mass vs pt vs Z vertex position
  TH3F*   fLambdaAssocMassPtDCADaug[kNCent];      //! Lambda Assoc: mass vs pt vs dca btween daughters
  TH3F*   fLambdaAssocMassPtCPA[kNCent];          //! Lambda Assoc: mass vs pt vs cpa
  TH3F*   fLambdaAssocMassPtDCAPV[kNCent];        //! Lambda Assoc: mass vs pt vs dca to prim vtx
  TH3F*   fLambdaAssocMassPtDaugNClsTPC[kNCent];  //! Lambda Assoc: mass vs pt vs num.of tpc clusters

  TH3F*   fLambdaAssocMassPtRapEmbeded[kNCent];          //! Lambda Assoc: pt vs rapidity vs mass (embeded)
  TH3F*   fLambdaAssocMassPtRapEmbeded2[kNCent];         //! Lambda Assoc: pt vs rapidity vs mass  (wo Cross contamination) (embeded)
  TH3F*   fLambdaAssocMassPtVtxEmbeded[kNCent];          //! Lambda Assoc: mass vs pt vs Z vertex position  (embeded particles)
  TH3F*   fLambdaAssocMassPtDCADaugEmbeded[kNCent];      //! Lambda Assoc: mass vs pt vs dca between daughters  (embeded particles)
  TH3F*   fLambdaAssocMassPtCPAEmbeded[kNCent];          //! Lambda Assoc: mass vs pt vs cpa  (embeded particles)
  TH3F*   fLambdaAssocMassPtDCAPVEmbeded[kNCent];        //! Lambda Assoc: mass vs pt vs dca to prim vtx  (embeded particles)
  TH3F*   fLambdaAssocMassPtDaugNClsTPCEmbeded[kNCent];  //! Lambda Assoc: mass vs pt vs num. of tpc clusters  (embeded particles)

  TH3F*   fLambdaMCResEta;               //! Lambda Assoc: eta resolution
  TH3F*   fLambdaMCResPhi;               //! Lambda Assoc: phi resolution


  //           AntiLambda            //
  TH1F*   fAntiLambdaMCPt;                  //! AntiLambda MC: pt
  TH3F*   fAntiLambdaMCPtRap;               //! AntiLambda MC: pt vs rapidity
  TH3F*   fAntiLambdaMCPtRap2;              //! AntiLambda MC: pt vs rapidity  (is Natural)
  TH3F*   fAntiLambdaMCPtRapVtx;            //! AntiLambda MC: pt vs rapidity vs Z vtx position
  TH3F*   fAntiLambdaMCPtRapEmbeded;        //! AntiLambda MC: pt vs rapidity (embeded particles)
  TH3F*   fAntiLambdaMCPtRapVtxEmbeded;     //! AntiLambda MC: pt vs rapidity  vs Z vtx position 
  TH2F*   fAntiLambdaMCFromXi;              //! AntiLambda MC: coming from Xi
  TH3F*   fAntiLambdaMCPtPhiEta[kNCent];    //! AntiLambda MC: pt vs pseudo-rapidity

  TH1F*   fAntiLambdaAssocPt;                 //! AntiLambda Assoc: pt
  TH3F*   fAntiLambdaAssocPtRap;              //! AntiLambda Assoc: pt vs rapidity vscentrality
  TH2F*   fAntiLambdaAssocFromXi;             //! AntiLambda Assoc: coming from Xi
  TH3F*   fAntiLambdaAssocPtPhiEta[kNCent];   //! AntiLambda Assoc: pt vs pseudo-rapidity

  TH3F*   fAntiLambdaAssocMassPtRap[kNCent];          //! AntiLambda Assoc: mass vs pt vs rapidity
  TH3F*   fAntiLambdaAssocMassPtRap2[kNCent];         //! AntiLambda Assoc: mass vs pt vs rapidity  (wo Cross contamination)
  TH3F*   fAntiLambdaAssocMassPtVtx[kNCent];          //! AntiLambda Assoc: mass vs pt vs Z vtx position
  TH3F*   fAntiLambdaAssocMassPtDCADaug[kNCent];      //! AntiLambda Assoc: mass vs pt vs Dca between daughters
  TH3F*   fAntiLambdaAssocMassPtCPA[kNCent];          //! AntiLambda Assoc: mass vs pt vs cpa
  TH3F*   fAntiLambdaAssocMassPtDCAPV[kNCent];        //! AntiLambda Assoc: mass vs pt vs dca to prim. vtx
  TH3F*   fAntiLambdaAssocMassPtDaugNClsTPC[kNCent];  //! AntiLambda Assoc: mass vs pt vs num. of tpc clusters

  TH3F*   fAntiLambdaAssocMassPtRapEmbeded[kNCent];          //! AntiLambda Assoc: mass vs pt vs rapidity  (embeded)
  TH3F*   fAntiLambdaAssocMassPtRapEmbeded2[kNCent];         //! AntiLambda Assoc: mass vs pt vs rapidity  (wo Cross contamination) (embeded)
  TH3F*   fAntiLambdaAssocMassPtVtxEmbeded[kNCent];          //! AntiLambda Assoc: mass vs pt vs Z vtx. position  (embeded particles)
  TH3F*   fAntiLambdaAssocMassPtDCADaugEmbeded[kNCent];      //! AntiLambda Assoc: mass vs pt vs dca between daughters  (embeded particles)
  TH3F*   fAntiLambdaAssocMassPtCPAEmbeded[kNCent];          //! AntiLambda Assoc: mass vs pt vs cpa  (embeded particles)
  TH3F*   fAntiLambdaAssocMassPtDCAPVEmbeded[kNCent];        //! AntiLambda Assoc: mass vs pt vs dca to prim. vtx  (embeded particles)
  TH3F*   fAntiLambdaAssocMassPtDaugNClsTPCEmbeded[kNCent];  //! AntiLambda Assoc: mass vs pt vs num. of tpc clusters  (embeded particles)

  TH3F*   fAntiLambdaMCResEta;             //! AntiLambda Assoc: eta resolution
  TH3F*   fAntiLambdaMCResPhi;             //! AntiLambda Assoc: phi resolution


  /// ====== Histograms for Correlations ====== ///

  TH3F*   fHistArmenterosPodolanski;     //! Armenteros-Podolanski plot inside 3 sigma of the signal
  TH3F*   fHistArmPodBckg;               //! Armenteros-Podolanski plot outside 3 sigma of the signal      



  //           K0s            //
  TH3F*   fK0sMass;                      //! Mass for K0s
  TH3F*   fK0sMassEmbeded;               //! Mass for K0s embeded
  TH3F*   fK0sMassPtEta;                 //! K0s: mass vs pt vs eta
  TH3F*   fK0sMassPtRap[kNCent];         //! K0s: mass vs pt vs rap vs centrality
  TH3F*   fK0sMassPtPhi;                 //! K0s: mass vs pt vs phi

  TH2F*   fK0sDaughtersPt;               //! K0s: pt of daughters
  TH3F*   fK0sDCADaugToPrimVtx;          //! K0s: DCA to primary vertex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fK0sSpatialRes;                //! K0s: Spatial resolution  
   
  TH3F*   fK0sdPhidEtaMC[kNCent*kN1];             //! K0s MC: Delta phi,Delta eta vs Z vertex position
  TH3F*   fK0sdPhidEtaPtL[kNVtxZ*kNCent*kN1];     //! K0s: Delta phi,Delta eta vs Z vertex position
  //TH3F*   fK0sdPhidEtaPtLBckg[kNCent*kN1];      //! K0s background: Delta phi,Delta eta vs Z vertex position
 
  TH2F*   fK0sBckgDecLength;             //! K0s background: Decay lenght vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fK0sBckgDCADaugToPrimVtx;      //! K0s background: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH2F*   fK0sBckgEtaPhi;                //! K0s background: Phi vs Eta inside a radio wrt the near-side peak
  TH2F*   fK0sBckgPhiRadio;              //! K0s background: Phi vs radio inside a radio wrt the near-side peak
  TH2F*   fK0sBckgDCANegDaugToPrimVtx;   //! K0s background: DCA of Negative daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  TH2F*   fK0sBckgDCAPosDaugToPrimVtx;   //! K0s background: DCA of Positive daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  TH2F*   fV0MassCascade;                //! V0s candiates: Possible mismatching of tracks due to cascades decays


  //           Lambda            //
  TH3F*   fLambdaMass;                   //! Mass for Lambda
  TH3F*   fLambdaMassEmbeded;            //! Mass for Lambda embeded
  TH3F*   fLambdaMass2;                  //! Mass for Lambda (rejecting crosscontamination)
  TH3F*   fLambdaMass2Embeded;           //! Mass for Lambda embded (rejecting crosscontamination)
  TH3F*   fLambdaMassPtEta;              //! Lambda: mass vs pt vs eta
  TH3F*   fLambdaMassPtRap[kNCent];      //! Lambda: mass vs pt vs rap
  TH3F*   fLambdaMassPtPhi;              //! Lambda: mass vs pt vs phi 

  TH2F*   fLambdaDaughtersPt;            //! Lambda: pt of daughters
  TH3F*   fLambdaDCADaugToPrimVtx;       //! Lambda: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fLambdaSpatialRes;             //! Lambda: Spatial resolution  

  TH3F*   fLambdadPhidEtaMC[kNCent*kN1];            //! Lambda MC: Delta phi,Delta eta vs Z vertex position
  TH3F*   fLambdadPhidEtaPtL[kNVtxZ*kNCent*kN1];    //! Lambda: Delta phi,Delta eta vs Z vertex position
  //TH3F*   fLambdadPhidEtaPtLBckg[kNCent*kN1];     //! Lambda background: Delta phi,Delta eta vs Z vertex position
 

  TH2F*   fLambdaBckgDecLength;            //! Lambda background: Decay lenght vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fLambdaBckgDCADaugToPrimVtx;     //! Lambda background: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH2F*   fLambdaBckgEtaPhi;               //! Lambda background: Phi vs Eta inside a radio wrt the near-side peak
  TH2F*   fLambdaBckgPhiRadio ;            //! Lambda background: Phi vs radio inside a radio wrt the near-side peak
  TH2F*   fLambdaBckgDCANegDaugToPrimVtx;  //! Lambda background: DCA of Negative daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  TH2F*   fLambdaBckgDCAPosDaugToPrimVtx;  //! Lambda background: DCA of Positive daughter to the primary vertex inside the radio 0.4 wrt the near-side peak


  //           AntiLambda            //
  TH3F*   fAntiLambdaMass;                     //! Mass for AntiLambda
  TH3F*   fAntiLambdaMassEmbeded;              //! Mass for AntiLambda embeded
  TH3F*   fAntiLambdaMass2;                    //! Mass for AntiLambda (rejecting crosscontamination)
  TH3F*   fAntiLambdaMass2Embeded;             //! Mass for AntiLambda embded (rejecting crosscontamination)

  TH3F*   fAntiLambdaMassPtEta;                //! AntiLambda: pt vs eta
  TH3F*   fAntiLambdaMassPtRap[kNCent];        //! AntiLambda: pt vs rap
  TH3F*   fAntiLambdaMassPtPhi;                //! Lambda: mass vs phi 

  TH2F*   fAntiLambdaDaughtersPt;              //! AntiLambda: pt of daughters
  TH3F*   fAntiLambdaDCADaugToPrimVtx;         //! AntiLambda: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fAntiLambdaSpatialRes;               //! AntiLambda: Spatial resolution  

  TH3F*   fAntiLambdadPhidEtaMC[kNCent*kN1];            //! AntiLambda MC: Delta phi,Delta eta vs Z vertex position
  TH3F*   fAntiLambdadPhidEtaPtL[kNVtxZ*kNCent*kN1];    //! AntiLambda: Delta phi,Delta eta vs pt of the leading particle
  //TH3F*   fAntiLambdadPhidEtaPtLBckg[kNCent*kN1];     //! AntiLambda background: Delta phi,Delta eta vs Z vertex position

  TH2F*   fAntiLambdaBckgDecLength;            //! AntiLambda background: Decay lenght vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fAntiLambdaBckgDCADaugToPrimVtx;     //! AntiLambda background: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH2F*   fAntiLambdaBckgEtaPhi;               //! AntiLambda background: Phi vs Eta inside a radio wrt the near-side peak
  TH2F*   fAntiLambdaBckgPhiRadio ;            //! AntiLambda background: Phi vs radio inside a radio wrt the near-side peak
  TH2F*   fAntiLambdaBckgDCANegDaugToPrimVtx;  //! AntiLambda background: DCA of Negative daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  TH2F*   fAntiLambdaBckgDCAPosDaugToPrimVtx;  //! AntiLambda background: DCA of Positive daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  
  //           Xi Minus          //
  TH2F*   fXiMinusPtMCAssoc;                         //! Xi Minus MC: Pt vs Centrality when they are associated particles
  TH2F*   fXiMinusPtMCTrigger;                       //! Xi Minus MC: Pt vs Centrality when they are trigger particles
  TH3F*   fXiMinusdPhidEtaMC[kNCent*kN2];            //! Xi Minus MC: Delta phi,Delta eta vs Z vertex position
  TH3F*   fXiMinusdPhidEtaMC2[kNCent*kN2];           //! Xi Minus MC: Delta phi,Delta eta vs Z vertex position
  TH3F*   fXiMinusdPhidEtaMC3[kNCent*kN3];           //! Xi Minus MC: Delta phi,Delta eta vs Z vertex position

  //        Gamma converison      //
  TH3F*  fGammaConversiondPhidEta[kNCent];        //! Gamma conversion: Delta phi,Delta eta vs Z vertex position

  //        Charged particles     //
  TH2F*  fChargeddPhidEta[kNVtxZ*kNCent*kNc];     //! Charged particles: Delta phi,Delta eta 
    
  ///  ==== Quality Assurance plots === ///

  //           K0s            //
  TH2F*   fK0sPtPosDaug;                     //! K0s: Pos. pt
  TH2F*   fK0sPtNegDaug;                     //! K0s: Neg. pt
  TH2F*   fK0sBckgPtPosDaug;                 //! K0s Bckg: Pos. pt
  TH2F*   fK0sBckgPtNegDaug;                 //! K0s Bckg: Neg. pt

  TH3F*   fK0sPhiEtaPosDaug;                 //! K0s: Pos. track phi vs eta 
  TH3F*   fK0sPhiEtaNegDaug;                 //! K0s: Neg. track phi vs eta
  TH3F*   fK0sBckgPhiEtaPosDaug;             //! K0s Bckg: Pos. track phi vs eta  
  TH3F*   fK0sBckgPhiEtaNegDaug;             //! K0s Bckg: Neg. track phi vs eta

  TH2F*   fK0sDCAPosDaug;                    //! K0s: Pos. track DCA to primary vertex
  TH2F*   fK0sDCANegDaug;                    //! K0s: Neg. track DCA to primary vertex
  TH2F*   fK0sBckgDCAPosDaug;                //! K0s Bckg: Pos. track DCA to primary vertex
  TH2F*   fK0sBckgDCANegDaug;                //! K0s Bckg: Neg. track DCA to primary vertex

  TH3F*   fK0sDecayPos;                      //! K0s: 2D decay position  
  TH3F*   fK0sBckgDecayPos;                  //! K0s Bckg: 2D decay position   
  TH2F*   fK0sDecayVertex;                   //! K0s: decay lenght
  TH2F*   fK0sBckgDecayVertex;               //! K0s Bckg: decay lenght 

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
  TH2F*   fLambdaPtPosDaug;                     //! Lambda: Pos. pt
  TH2F*   fLambdaPtNegDaug;                     //! Lambda: Neg. pt
  TH2F*   fLambdaBckgPtPosDaug;                 //! Lambda Bckg: Pos. pt
  TH2F*   fLambdaBckgPtNegDaug;                 //! Lambda Bckg: Neg. pt

  TH3F*   fLambdaPhiEtaPosDaug;                 //! Lambda: Pos. track phi vs eta 
  TH3F*   fLambdaPhiEtaNegDaug;                 //! Lambda: Neg. track phi vs eta
  TH3F*   fLambdaBckgPhiEtaPosDaug;             //! Lambda Bckg: Pos. track phi vs eta  
  TH3F*   fLambdaBckgPhiEtaNegDaug;             //! Lambda Bckg: Neg. track phi vs eta

  TH2F*   fLambdaDCAPosDaug;                    //! Lambda: Pos. track DCA to primary vertex
  TH2F*   fLambdaDCANegDaug;                    //! Lambda: Neg. track DCA to primary vertex
  TH2F*   fLambdaBckgDCAPosDaug;                //! Lambda Bckg: Pos. track DCA to primary vertex
  TH2F*   fLambdaBckgDCANegDaug;                //! Lambda Bckg: Neg. track DCA to primary vertex

  TH3F*   fLambdaDecayPos;                      //! Lambda: 2D decay position  
  TH3F*   fLambdaBckgDecayPos;                  //! Lambda Bckg: 2D decay position   
  TH2F*   fLambdaDecayVertex;                   //! Lambda: decay lenght
  TH2F*   fLambdaBckgDecayVertex;               //! Lambda Bckg: decay lenght 

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
  TH2F*   fAntiLambdaPtPosDaug;                     //! AntiLambda: Pos. pt
  TH2F*   fAntiLambdaPtNegDaug;                     //! AntiLambda: Neg. pt
  TH2F*   fAntiLambdaBckgPtPosDaug;                 //! AntiLambda Bckg: Pos. pt
  TH2F*   fAntiLambdaBckgPtNegDaug;                 //! AntiLambda Bckg: Neg. pt

  TH3F*   fAntiLambdaPhiEtaPosDaug;                 //! AntiLambda: Pos. track phi vs eta 
  TH3F*   fAntiLambdaPhiEtaNegDaug;                 //! AntiLambda: Neg. track phi vs eta
  TH3F*   fAntiLambdaBckgPhiEtaPosDaug;             //! AntiLambda Bckg: Pos. track phi vs eta  
  TH3F*   fAntiLambdaBckgPhiEtaNegDaug;             //! AntiLambda Bckg: Neg. track phi vs eta

  TH2F*   fAntiLambdaDCAPosDaug;                    //! AntiLambda: Pos. track DCA to primary vertex
  TH2F*   fAntiLambdaDCANegDaug;                    //! AntiLambda: Neg. track DCA to primary vertex
  TH2F*   fAntiLambdaBckgDCAPosDaug;                //! AntiLambda Bckg: Pos. track DCA to primary vertex
  TH2F*   fAntiLambdaBckgDCANegDaug;                //! AntiLambda Bckg: Neg. track DCA to primary vertex

  TH3F*   fAntiLambdaDecayPos;                      //! AntiLambda: 2D decay position  
  TH3F*   fAntiLambdaBckgDecayPos;                  //! AntiLambda Bckg: 2D decay position   
  TH2F*   fAntiLambdaDecayVertex;                   //! AntiLambda: decay lenght
  TH2F*   fAntiLambdaBckgDecayVertex;               //! AntiLambda Bckg: decay lenght

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


  ///  ==== Mixed Events plots === ///
  TH2F*  fChargeddPhidEtaME[kNVtxZ*kNCent*kNc+1];             //! K0s Mixed Events
  TH2F*  fK0sdPhidEtaME[kNVtxZ*kNCent*kN1+1];             //! K0s Mixed Events
  TH2F*  fLambdadPhidEtaME[kNVtxZ*kNCent*kN1+1];          //! Lambda Mixed Events
  TH2F*  fAntiLambdadPhidEtaME[kNVtxZ*kNCent*kN1+1];      //! AntiLambda Mixed Events

  ClassDef(AliAnalysisTaskLambdaOverK0sJets,1);

};


/*  
    Based on AliV0ChBasicParticle class of AliAnalysisTaskV0ChCorrelations.
    Keeps basic information to reduce memory consumption for event mixing.
*/
class AliMiniParticle : public AliVParticle
{
 public:
 AliMiniParticle(Float_t centrality, Float_t vtxZ, Int_t id,Double_t pt, Double_t phi, 
		 Double_t eta, Int_t negDaugMC, Int_t posDaugMC, Short_t candidate)
   :fCentrality(centrality), fVtxZ(vtxZ),  fId(id), fPt(pt),
    fPhi(phi), fEta(eta), fNegDaugMC(negDaugMC), fPosDaugMC(posDaugMC), fCandidate(candidate)
  {
  }
  
  virtual ~AliMiniParticle() {}
  
  // event
  virtual Float_t Centrality() const { return fCentrality; }
  virtual Float_t VtxZ() const { return fVtxZ; }

  virtual Int_t   ID()  const { return fId; }  
  // kinematics
  virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }

  virtual Double_t Pt() const { return fPt; }
  virtual Double_t P()  const { AliFatal("Not implemented"); return 0; }
  virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

  virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
  virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

  virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
    
  virtual Double_t Phi()        const { return fPhi; }
  virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }
  virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
  virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }
    
  virtual Double_t Eta()        const { return fEta; }
  virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }

  virtual Short_t Charge()      const { AliFatal("Not implemented"); return 0; }
  virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
  // PID
  virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
  virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
  virtual Int_t   NegDaugMCLabel() const { return fNegDaugMC; } 
  virtual Int_t   PosDaugMCLabel() const { return fPosDaugMC; }
  virtual Short_t WhichCandidate() const { return fCandidate; }
  
 private:
  Float_t fCentrality; // centrality of the event
  Float_t fVtxZ;       // vertex postition in the event
  Int_t   fId;         // ID related either to AliAODtrack or AliAODv0
  Float_t fPt;         // pt 
  Float_t fPhi;        // phi
  Float_t fEta;        // eta 
  Int_t   fNegDaugMC;  // MC origin of negative daughter
  Int_t   fPosDaugMC;  // MC origin of positive daughter
  Short_t fCandidate;  // Candidate: 0-Not trigger, 1-Trigger, 2-Gamma Conversion, 3-K0s candidates, 4-Lambda candidates, 5-AntiLambda candidates
  
  ClassDef( AliMiniParticle, 1); // class required for event mixing
};

#endif
