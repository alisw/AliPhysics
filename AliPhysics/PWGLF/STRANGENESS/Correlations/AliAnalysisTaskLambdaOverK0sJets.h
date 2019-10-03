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

//try

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
class THnSparse;
class TList;
class TString;

class TObjArray;

// pt for V0
const int    kN1 = 6;
const float  kPtBinV0[kN1+1] = {1.5, 2.0, 2.6, 3.2, 4.0, 5.5, 8.0};

const int    kNVtxZ = 10; 
const double kBinVtxZ[kNVtxZ+1] = {-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.};

const int    kNCent  = 4;
const double kBinCent[kNCent+1] = {0.0,5.0,10.0,20.0,40.0};

//  ------------------------------------
//  Inv. Mass width as function of the centrality
//  Linear polimomial dependence: sigma(pt) = a0 * a1*pt

const double kCteK0s2010[kNCent] = {0.00348, 0.00351, 0.00346, 0.00318};
const double kLinearK0s2010[kNCent] = {8.024E-4, 7.403E-4, 7.250E-4, 7.665E-4};

const double kCteK0s2011[kNCent] = {0.00338, 0.00328, 0.00333, 0.00326};
const double kLinearK0s2011[kNCent] = {8.336E-4, 8.385E-4, 7.891E-4, 7.851E-4};

const double kCteLambda2010[kNCent] = {0.00145, 0.00122, 0.00140, 0.00135};
const double kLinearLambda2010[kNCent] = {2.233E-4, 2.836-4, 2.105-4, 2.076E-4};

const double kCteLambda2011[kNCent] = {0.00130, 0.00123, 0.00114, 0.00121};
const double kLinearLambda2011[kNCent] = {3.002E-4, 3.067E-4, 3.207E-4, 2.813E-4};

const double kCteAntiLambda2010[kNCent] = {0.00109, 0.00134, 0.00117, 0.00116};
const double kLinearAntiLambda2010[kNCent] = {3.245E-4, 2.308E-4, 2.707E-4, 2.562E-4};

const double kCteAntiLambda2011[kNCent] = {9.859E-4, 0.00111, 0.00104, 0.00110};
const double kLinearAntiLambda2011[kNCent] = {3.881E-4, 3.379E-4, 3.490E-4, 3.166E-4};

// -------------------------------------

class AliAnalysisTaskLambdaOverK0sJets : public AliAnalysisTaskSE {

 public:
  
  enum V0LoopStep_t  { kTriggerCheck=1, kReconstruction=2 };

  AliAnalysisTaskLambdaOverK0sJets(const char *name = "AliAnalysisTaskLambdaOverK0sJets");
  virtual ~AliAnalysisTaskLambdaOverK0sJets();

  // Setter for global variables in the event
  void SetCollisionType(TString data="PbPb2010") {fCollision=data;}
  void SetMC(Bool_t isMC=kTRUE) {fIsMC=isMC;} 
  void SetPID(Bool_t usePID=kTRUE) {fUsePID=usePID;} 
  void SetCentrality(Float_t min=0., Float_t max=90.) {fCentMin=min;fCentMax=max;} 
  void SetQA(Bool_t doQA=kFALSE){fDoQA=doQA;}
  void SetDoMix(Bool_t doMixEvt=kTRUE) {fDoMixEvt=doMixEvt;} 
  void SetTriggerFilterBit(Int_t triggerFB=768){fTriggerFB=triggerFB;}
  void SetTriggerPt(Float_t ptMinTrig=8., Float_t ptMaxTrig=50.) {fTrigPtMin=ptMinTrig;fTrigPtMax=ptMaxTrig;} 
  void SetTriggerEta(Float_t etaMaxTrig=0.8){fTrigEtaMax=etaMaxTrig;}
  void SetTriggerNCls(Float_t nclsTrig=10){fTriggerNCls=nclsTrig;}
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
  void SetEtaCut(Bool_t etaCut=kFALSE) {fUseEtaCut=etaCut;} 
  void SetMaxY(Float_t yMax=0.5) {fYMax=yMax;} 
  void SetMinCPA(Float_t minCPA=0.998) {fMinCPA=minCPA;} 
  void SetCtau(Float_t minCtau = 0., Float_t maxCtau = 3.) {fMinCtau=minCtau;fMaxCtau=maxCtau;} 

  // Setting variables for splitting cut
  void SetTPCRadius(Float_t tpcRadius=125.) {fTPCRadius=tpcRadius;}    
  void SetFracSharedTPCcls(Float_t fracSharedTPCcls=0.4) {fFracTPCcls=fracSharedTPCcls;}                     
  void SetDiffSharedTPCcls(Float_t diffSharedTPCcls=0.06) {fDiffTrigDaugFracTPCSharedCls=diffSharedTPCcls;}

  // Getters
  Float_t GetMinCentr() { return fCentMin; }
  Float_t GetMaxCentr() { return fCentMax; }

  // Main functions
  virtual void     UserCreateOutputObjects();
  virtual Bool_t   AcceptTrack(const AliAODTrack *t); 
  virtual Bool_t   AcceptTrackV0(const AliAODTrack *t);
  virtual Bool_t   AcceptV0(AliAODVertex *vtx, const AliAODv0 *v0);
  virtual Float_t  GetMultiplicity();
  virtual Double_t ThetaS(TString part);
  virtual Double_t EtaS(TString part);
  virtual Float_t  dEtaS();
  virtual Float_t  dPhiSAtR125();
  virtual void     SetSftPosR125(const AliAODTrack *track,const Float_t bfield,const Float_t priVtx[3], TString part);
  virtual void     RecCascade(const AliAODTrack *trk1,const AliAODTrack *trk2,const AliAODTrack *trkBch,TString histo);
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
  Int_t    fTriggerFB;                   //  Trigger track filter bit
  Float_t  fTrigPtMin;                   //  Minimum pt for trigger particle
  Float_t  fTrigPtMax;                   //  Maximum pt for trigger particle
  Float_t  fTrigPtMCMin;                 //  Minimum pt for trigger particle in MC
  Float_t  fTrigPtMCMax;                 //  Maximum pt for trigger particle in MC
  Float_t  fTrigEtaMax;                  //  Maximum eta for trigger particle
  Float_t  fTriggerNCls;                 //  Number of Crossed pas rows in the TPC
  Bool_t   fCheckIDTrig;                 //  Do comparison with V0's daughter tracks?
  Bool_t   fSeparateInjPart;             //  Separate MC injected particles in case of correlation 
  Int_t    fEndOfHijingEvent;            //  Limit natural-injected MC  particles 
  AliPIDResponse *fPIDResponse;          //  PID Response

  Float_t fMinPtDaughter;                //  Minimum transverse momentum for V0's daughters
  Float_t fMaxEtaDaughter;               //  Maximum pseudo-rapidity for V0's daughters  
  Float_t fMaxDCADaughter;               //  Maximum Distance of Closest Approach between daughters (given in sigmas)
  Bool_t  fUseEtaCut;                    //  Swicth between rapidity or pseudo-rapidity cut
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

  Float_t fTPCRadius;                    // Radial position of TPC to obtain the separation between the trigger particle and the daughter particle
  Float_t fTrigSftR125[3];               // Shifted position of the daughter track to the Primary verterx
  Float_t fDaugSftR125[3];               // Shifted position of the trigger track to the Primary verterx
  Float_t fFracTPCcls;                   // Threshold for the fraction of TPC shared clusters for single track
  Float_t fDiffTrigDaugFracTPCSharedCls; // Allowed dispertion in the fraction of TPC shared clusters between trigger particle and the daughter track


  TList*  fOutput;                       //! List of histograms for main analysis
  TList*  fOutputQA;                     //! List of histograms for Quality Assurance
  TList*  fOutputME;                     //! List of histograms for Mixed Events
  TList** fMEList;                       //![] List of Mixed Events

  TObjArray* fTriggerParticles;          // Trigger particle array
  TObjArray* fTriggerPartMC;             // MC Trigger particle array
  TObjArray* fAssocParticles;            // Associated particle array
  TObjArray* fAssocPartMC;               // MC Associated particle array
  
  TH1F*   fEvents;                       //! Counter for the number of events in each step
  TH2F*   fEvtPerCent;                   //! Counter for the number of events in each step per centrality bin
  TH1F*   fCentrality;                   //! Event centrality per centil
  TH1F*   fCentrality2;                  //! Event centrality per centil with |VtxZ|<10cm
  TH2F*   fCentralityTrig;               //! Event centrality per trigger
  TH2F*   fPrimayVtxGlobalvsSPD;         //! Zvtx tracking vs Zvtx SPD
  TH1F*   fPrimaryVertexX;               //! Primary vertex position in X
  TH1F*   fPrimaryVertexY;               //! Primary vertex position in Y
  TH1F*   fPrimaryVertexZ;               //! Primary vertex position in Z

  TH2F*   fChargedMultiplicity;          //! Charged multiplicity vs centrality bin  
  TH1F*   fTriggerEventPlane;            //! Distance between the trigger particle direction and the event plane angle

  TH2F*   fTriggerMCPtCent;              //! Trigger particle MC: pt vs centrality
  TH3F*   fTriggerMCResPt;               //! Trigger particle MC: pt resolution
  TH3F*   fTriggerMCResEta;              //! Trigger particle MC: eta resolution
  TH3F*   fTriggerMCResPhi;              //! Trigger particle MC: phi resolution
  TH3F*   fTriggerPtCent;                //! Trigger particle: pt vs centrality vs Z vertex
  TH3F*   fTriggerPtCentCh;              //! Trigger particle: pt vs centrality vs Z vertex for hh correlations
  TH2F*   fNTrigPerEvt;                  //! Trigger particle: Number of particle triggers per event
  TH1F*   fTriggerWiSPDHit;              //! Trigger particle: Has Hits in the SPD?
  TH2F*   fTriggerEtaPhi;                //! Trigger particle: eta vs phi
  TH2F*   fTrigFracShTPCcls;             //! Trigger particle: pt vs fraction of shared TPC cls
  TH2F*   fTriggerDCA;                   //! Trigger particle: dca to primary vertex
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
  TH3F*   fK0sMCPtRapVtx[kNCent];        //! K0s MC: pt vs Z vtx position  vs centrality
  TH3F*   fK0sMCPtRapEmbeded;            //! K0s MC: pt vs rapidity  (embeded particles)
  TH3F*   fK0sMCPtRapVtxEmbeded[kNCent]; //! K0s MC: pt vs Z vtx position rapidity  vs centrality (embeded particles)
  THnSparse* fK0sMCPtRapPtDaugPt[kNCent]; //! K0s MC: pt vs rapidity vs transverse momemtum of daughters
  THnSparse* fK0sMCPtRapPtDaugPtEmbeded[kNCent]; //! K0s MC: pt vs rapidity vs transverse momemtum of daughters
  TH3F*   fK0sMCPtPhiEta[kNCent];        //! K0s MC: pt vs pseudo-rapidity

  TH1F*   fK0sAssocPt;                         //! K0s Assoc: pt
  TH3F*   fK0sAssocPtArm;                      //! K0s Assoc: pt vs rapidity vs centrality (arm. pod. cut)
  TH3F*   fK0sAssocPtRap;                      //! K0s Assoc: pt vs rapidity vs centrality
  TH3F*   fK0sAssocPtRapEmbeded;               //! K0s Assoc: pt vs rapidity vs centrality  (embeded particles)
  TH3F*   fK0sAssocPtPhiEta[kNCent];           //! K0s Assoc: pt vs pseudo-rapidity

  THnSparse*  fK0sAssocPtMassArm[kNCent];          //! K0s Assoc: mass vs pt vs centrality
  THnSparse*  fK0sAssocMassPtVtx[kNCent];          //! K0s Assoc: mass vs pt vs Z vertex position
  THnSparse*  fK0sAssocMassPtDCADaug[kNCent];      //! K0s Assoc: mass vs pt vs dca between daughters
  THnSparse*  fK0sAssocMassPtCPA[kNCent];          //! K0s Assoc: mass vs pt vs cpa
  THnSparse*  fK0sAssocMassPtDCAPV[kNCent];        //! K0s Assoc: mass vs pt vs dca to prim. vtx
  THnSparse*  fK0sAssocMassPtDaugNClsTPC[kNCent];  //! K0s Assoc: mass vs pt vs num. of tpc clusters
  THnSparse*  fK0sAssocMassPtShTPCcls[kNCent];     //! K0s Assoc: mass vs pt vs fraction of shared TPC cls
  THnSparse*  fK0sAssocMassPtDaugPt[kNCent];       //! K0s Assoc: mass vs pt vs transverse momemtum of daughters
  THnSparse*  fK0sAssocMassPtCtau[kNCent];         //! K0s Assoc: mass vs pt vs proper time of life
  THnSparse*  fK0sAssocMassPtFidVolume[kNCent];    //! K0s Assoc: mass vs pt vs fiducial volume

  THnSparse*  fK0sAssocPtMassArmEmbeded[kNCent];          //! K0s Assoc: mass vs pt vs rapidity  (embeded particles)
  THnSparse*  fK0sAssocMassPtVtxEmbeded[kNCent];          //! K0s Assoc: mass vs pt vs Z vertex position  (embeded particles)
  THnSparse*  fK0sAssocMassPtDCADaugEmbeded[kNCent];      //! K0s Assoc: mass vs pt vs dca between daughters  (embeded particles)
  THnSparse*  fK0sAssocMassPtCPAEmbeded[kNCent];          //! K0s Assoc: mass vs pt vs cpa  (embeded particles)
  THnSparse*  fK0sAssocMassPtDCAPVEmbeded[kNCent];        //! K0s Assoc: mass vs pt vs dca to prim. vtx  (embeded particles)
  THnSparse*  fK0sAssocMassPtDaugNClsTPCEmbeded[kNCent];  //! K0s Assoc: mass vs pt vs num. o ftpc clusters (embeded particles)
  THnSparse*  fK0sAssocMassPtShTPCclsEmbeded[kNCent];     //! K0s Assoc: mass vs pt vs fraction of shared TPC cls
  THnSparse*  fK0sAssocMassPtDaugPtEmbeded[kNCent];       //! K0s Assoc: mass vs pt vs transverse momemtum of daughters
  THnSparse*  fK0sAssocMassPtCtauEmbeded[kNCent];         //! K0s Assoc: mass vs pt vs proper time of life
  THnSparse*  fK0sAssocMassPtFidVolumeEmbeded[kNCent];    //! K0s Assoc: mass vs pt vs fiducial volume


  TH3F*   fK0sMCResEta;                 //! K0s Assoc: eta resolution
  TH3F*   fK0sMCResPhi;                 //! K0s Assoc: phi resolution
  TH3F*   fK0sMCResPt;                  //! K0s Assoc: pt resolution
  TH3F*   fK0sPosMCResEta;              //! K0s Pos. Daughter: eta resolution
  TH3F*   fK0sPosMCResPhi;              //! K0s Pos. Daughter: phi resolution
  TH3F*   fK0sPosMCResPt;               //! K0s Pos. Daughter: pt resolution
  TH3F*   fK0sNegMCResEta;              //! K0s Neg. Daughter: eta resolution
  TH3F*   fK0sNegMCResPhi;              //! K0s Neg. Daughter: phi resolution
  TH3F*   fK0sNegMCResPt;               //! K0s Neg. Daughter: pt resolution

  //           Lambda            //
  TH1F*   fLambdaMCPt;                        //! Lambda MC: pt
  TH3F*   fLambdaMCPtRap;                     //! Lambda MC: pt vs rapidity
  TH3F*   fLambdaMCPtRap2;                    //! Lambda MC: pt vs rapidity  (is Natural)
  TH3F*   fLambdaMCPtRapVtx[kNCent];          //! Lambda MC: pt vs Z vtx position rapidity vs centrality
  TH3F*   fLambdaMCPtRapEmbeded;              //! Lambda MC: pt vs rapidity (embeded particles)
  TH3F*   fLambdaMCPtRapVtxEmbeded[kNCent];   //! Lambda MC: pt vs Z vtx position vs centrality  (embeded particles)
  THnSparse*  fLambdaMCPtRapPtDaugPt[kNCent]; //! Lambda MC: pt vs rapidity vs transverse momemtum of daughters
  THnSparse*  fLambdaMCPtRapPtDaugPtEmbeded[kNCent]; //! Lambda MC: pt vs rapidity vs transverse momemtum of daughters
  TH2F*   fLambdaMCFromXi;                    //! Lambda MC: coming from Xi
  TH3F*   fLambdaMCPtPhiEta[kNCent];          //! Lambda MC: pt vs pseudo-rapidity

  TH1F*   fLambdaAssocPt;                //! Lambda Assoc: pt
  TH3F*   fLambdaAssocPtRap;             //! Lambda Assoc: pt vs rapidity
  TH2F*   fLambdaAssocFromXi;            //! Lambda Assoc: coming from Xi
  TH3F*   fLambdaAssocPtPhiEta[kNCent];  //! Lambda Assoc: pt vs pseudo-rapidity

  THnSparse*  fLambdaAssocMassPtRap[kNCent];          //! Lambda Assoc: pt vs rapidity vs mass
  THnSparse*  fLambdaAssocMassPtRap2[kNCent];         //! Lambda Assoc: pt vs rapidity vs mass (wo Cross contamination)
  THnSparse*  fLambdaAssocMassPtVtx[kNCent];          //! Lambda Assoc: mass vs pt vs Z vertex position
  THnSparse*  fLambdaAssocMassPtDCADaug[kNCent];      //! Lambda Assoc: mass vs pt vs dca btween daughters
  THnSparse*  fLambdaAssocMassPtCPA[kNCent];          //! Lambda Assoc: mass vs pt vs cpa
  THnSparse*  fLambdaAssocMassPtDCAPV[kNCent];        //! Lambda Assoc: mass vs pt vs dca to prim vtx
  THnSparse*  fLambdaAssocMassPtDaugNClsTPC[kNCent];  //! Lambda Assoc: mass vs pt vs num.of tpc clusters
  THnSparse*  fLambdaAssocMassPtShTPCcls[kNCent];     //! Lambda Assoc: mass vs pt vs fraction of shared TPC cls
  THnSparse*  fLambdaAssocMassPtDaugPt[kNCent];       //! Lambda Assoc: mass vs pt vs transverse momemtum of daughters
  THnSparse*  fLambdaAssocMassPtCtau[kNCent];         //! Lambda Assoc: mass vs pt vs proper time of life
  THnSparse*  fLambdaAssocMassPtFidVolume[kNCent];    //! Lambda Assoc: mass vs pt vs fiducial volume

  THnSparse*  fLambdaAssocMassPtRapEmbeded[kNCent];          //! Lambda Assoc: pt vs rapidity vs mass (embeded)
  THnSparse*  fLambdaAssocMassPtRapEmbeded2[kNCent];         //! Lambda Assoc: pt vs rapidity vs mass  (wo Cross contamination) (embeded)
  THnSparse*  fLambdaAssocMassPtVtxEmbeded[kNCent];          //! Lambda Assoc: mass vs pt vs Z vertex position  (embeded particles)
  THnSparse*  fLambdaAssocMassPtDCADaugEmbeded[kNCent];      //! Lambda Assoc: mass vs pt vs dca between daughters  (embeded particles)
  THnSparse*  fLambdaAssocMassPtCPAEmbeded[kNCent];          //! Lambda Assoc: mass vs pt vs cpa  (embeded particles)
  THnSparse*  fLambdaAssocMassPtDCAPVEmbeded[kNCent];        //! Lambda Assoc: mass vs pt vs dca to prim vtx  (embeded particles)
  THnSparse*  fLambdaAssocMassPtDaugNClsTPCEmbeded[kNCent];  //! Lambda Assoc: mass vs pt vs num. of tpc clusters  (embeded particles)
  THnSparse*  fLambdaAssocMassPtShTPCclsEmbeded[kNCent];     //! Lambda Assoc: mass vs pt vs fraction of shared TPC cls
  THnSparse*  fLambdaAssocMassPtDaugPtEmbeded[kNCent];       //! Lambda Assoc: mass vs pt vs transverse momemtum of daughters
  THnSparse*  fLambdaAssocMassPtCtauEmbeded[kNCent];         //! Lambda Assoc: mass vs pt vs proper time of life
  THnSparse*  fLambdaAssocMassPtFidVolumeEmbeded[kNCent];    //! Lambda Assoc: mass vs pt vs fiducial volume


  TH3F*   fLambdaMCResEta;                 //! Lambda Assoc: eta resolution
  TH3F*   fLambdaMCResPhi;                 //! Lambda Assoc: phi resolution
  TH3F*   fLambdaMCResPt;                  //! Lambda Assoc: pt resolution
  TH3F*   fLambdaPosMCResEta;              //! Lambda Pos. Daughter: eta resolution
  TH3F*   fLambdaPosMCResPhi;              //! Lambda Pos. Daughter: phi resolution
  TH3F*   fLambdaPosMCResPt;               //! Lambda Pos. Daughter: pt resolution
  TH3F*   fLambdaNegMCResEta;              //! Lambda Neg. Daughter: eta resolution
  TH3F*   fLambdaNegMCResPhi;              //! Lambda Neg. Daughter: phi resolution
  TH3F*   fLambdaNegMCResPt;               //! Lambda Neg. Daughter: pt resolution

  //           AntiLambda            //
  TH1F*   fAntiLambdaMCPt;                  //! AntiLambda MC: pt
  TH3F*   fAntiLambdaMCPtRap;               //! AntiLambda MC: pt vs rapidity
  TH3F*   fAntiLambdaMCPtRap2;              //! AntiLambda MC: pt vs rapidity  (is Natural)
  TH3F*   fAntiLambdaMCPtRapVtx[kNCent];          //! AntiLambda MC: pt vs rapidity vs Z vtx position
  TH3F*   fAntiLambdaMCPtRapEmbeded;              //! AntiLambda MC: pt vs rapidity (embeded particles)
  TH3F*   fAntiLambdaMCPtRapVtxEmbeded[kNCent];   //! AntiLambda MC: pt vs rapidity  vs Z vtx position 
  THnSparse*  fAntiLambdaMCPtRapPtDaugPt[kNCent]; //! AntiLambda MC: pt vs rapidity vs transverse momemtum of daughters
  THnSparse*  fAntiLambdaMCPtRapPtDaugPtEmbeded[kNCent]; //! AntiLambda MC: pt vs rapidity vs transverse momemtum of daughters
  TH2F*   fAntiLambdaMCFromXi;                    //! AntiLambda MC: coming from Xi
  TH3F*   fAntiLambdaMCPtPhiEta[kNCent];          //! AntiLambda MC: pt vs pseudo-rapidity

  TH1F*   fAntiLambdaAssocPt;                 //! AntiLambda Assoc: pt
  TH3F*   fAntiLambdaAssocPtRap;              //! AntiLambda Assoc: pt vs rapidity vscentrality
  TH2F*   fAntiLambdaAssocFromXi;             //! AntiLambda Assoc: coming from Xi
  TH3F*   fAntiLambdaAssocPtPhiEta[kNCent];   //! AntiLambda Assoc: pt vs pseudo-rapidity

  THnSparse*  fAntiLambdaAssocMassPtRap[kNCent];          //! AntiLambda Assoc: mass vs pt vs rapidity
  THnSparse*  fAntiLambdaAssocMassPtRap2[kNCent];         //! AntiLambda Assoc: mass vs pt vs rapidity  (wo Cross contamination)
  THnSparse*  fAntiLambdaAssocMassPtVtx[kNCent];          //! AntiLambda Assoc: mass vs pt vs Z vtx position
  THnSparse*  fAntiLambdaAssocMassPtDCADaug[kNCent];      //! AntiLambda Assoc: mass vs pt vs Dca between daughters
  THnSparse*  fAntiLambdaAssocMassPtCPA[kNCent];          //! AntiLambda Assoc: mass vs pt vs cpa
  THnSparse*  fAntiLambdaAssocMassPtDCAPV[kNCent];        //! AntiLambda Assoc: mass vs pt vs dca to prim. vtx
  THnSparse*  fAntiLambdaAssocMassPtDaugNClsTPC[kNCent];  //! AntiLambda Assoc: mass vs pt vs num. of tpc clusters
  THnSparse*  fAntiLambdaAssocMassPtShTPCcls[kNCent];     //! AntiLambda Assoc: mass vs pt vs fraction of shared TPC cls
  THnSparse*  fAntiLambdaAssocMassPtDaugPt[kNCent];       //! AntiLambda Assoc: mass vs pt vs transverse momemtum of daughters
  THnSparse*  fAntiLambdaAssocMassPtCtau[kNCent];         //! AntiLambda Assoc: mass vs pt vs proper time of life
  THnSparse*  fAntiLambdaAssocMassPtFidVolume[kNCent];    //! AntiLambda Assoc: mass vs pt vs fiducial volume

  THnSparse*  fAntiLambdaAssocMassPtRapEmbeded[kNCent];          //! AntiLambda Assoc: mass vs pt vs rapidity  (embeded)
  THnSparse*  fAntiLambdaAssocMassPtRapEmbeded2[kNCent];         //! AntiLambda Assoc: mass vs pt vs rapidity  (wo Cross contamination) (embeded)
  THnSparse*  fAntiLambdaAssocMassPtVtxEmbeded[kNCent];          //! AntiLambda Assoc: mass vs pt vs Z vtx. position  (embeded particles)
  THnSparse*  fAntiLambdaAssocMassPtDCADaugEmbeded[kNCent];      //! AntiLambda Assoc: mass vs pt vs dca between daughters  (embeded particles)
  THnSparse*  fAntiLambdaAssocMassPtCPAEmbeded[kNCent];          //! AntiLambda Assoc: mass vs pt vs cpa  (embeded particles)
  THnSparse*  fAntiLambdaAssocMassPtDCAPVEmbeded[kNCent];        //! AntiLambda Assoc: mass vs pt vs dca to prim. vtx  (embeded particles)
  THnSparse*  fAntiLambdaAssocMassPtDaugNClsTPCEmbeded[kNCent];  //! AntiLambda Assoc: mass vs pt vs num. of tpc clusters  (embeded particles)
  THnSparse*  fAntiLambdaAssocMassPtShTPCclsEmbeded[kNCent];    //! AntiLambda Assoc: mass vs pt vs fraction of shared TPC cls
  THnSparse*  fAntiLambdaAssocMassPtDaugPtEmbeded[kNCent];       //! AntiLambda Assoc: mass vs pt vs transverse momemtum of daughters
  THnSparse*  fAntiLambdaAssocMassPtCtauEmbeded[kNCent];         //! AntiLambda Assoc: mass vs pt vs proper time of life
  THnSparse*  fAntiLambdaAssocMassPtFidVolumeEmbeded[kNCent];    //! AntiLambda Assoc: mass vs pt vs fiducial volume


  TH3F*   fAntiLambdaMCResEta;              //! AntiLambda Assoc: eta resolution
  TH3F*   fAntiLambdaMCResPhi;              //! AntiLambda Assoc: phi resolution
  TH3F*   fAntiLambdaMCResPt;               //! AntiLambda Assoc: pt resolution
  TH3F*   fAntiLambdaPosMCResEta;           //! AntiLambda Pos. Daughter: eta resolution
  TH3F*   fAntiLambdaPosMCResPhi;           //! AntiLambda Pos. Daughter: phi resolution
  TH3F*   fAntiLambdaPosMCResPt;            //! AntiLambda Pos. Daughter: pt resolution
  TH3F*   fAntiLambdaNegMCResEta;           //! AntiLambda Neg. Daughter: eta resolution
  TH3F*   fAntiLambdaNegMCResPhi;           //! AntiLambda Neg. Daughter: phi resolution
  TH3F*   fAntiLambdaNegMCResPt;            //! AntiLambda Neg. Daughter: pt resolution



  /// ====== Histograms for Correlations ====== ///

  TH3F*   fHistArmenterosPodolanski;     //! Armenteros-Podolanski plot inside 3 sigma of the signal
  TH3F*   fHistArmPodBckg;               //! Armenteros-Podolanski plot outside 3 sigma of the signal      

  //           K0s            //
  TH3F*   fK0sMass;                      //! Mass for K0s
  TH3F*   fK0sMassEmbeded;               //! Mass for K0s embeded
  TH3F*   fK0sMassPtEta;                 //! K0s: mass vs pt vs eta
  TH3F*   fK0sMassPtRap[kNCent];         //! K0s: mass vs pt vs rap vs centrality
  TH3F*   fK0sMassPtPhi;                 //! K0s: mass vs pt vs phi
  TH3F*   fK0sPosDaugFracShTPCcls;       //! K0s Pos. Daug: mass vs pt vs fraction of shared TPC cls
  TH3F*   fK0sNegDaugFracShTPCcls;       //! K0s Neg Daug: mass vs pt vs fraction of shared TPC cls

  TH2F*   fK0sDaughtersPt;                       //! K0s: pt of daughters
  THnSparse* fK0sPosDaugdPhiSdEtaS[kNCent];      //! K0s: Positive daughter: delta(phi)* delta(eta)*    
  THnSparse* fK0sNegDaugdPhiSdEtaS[kNCent];      //! K0s: Negative daughter: delta(phi)* delta(eta)* 
  THnSparse* fK0sPosDaugSplCheckCovMat[kNCent];  //! K0s: Check Covariance Matrix elemenets between trigger trcak and daughter track
  THnSparse* fK0sNegDaugSplCheckCovMat[kNCent];  //! K0s: Check Covariance Matrix elemenets between trigger trcak and daughter track
  THnSparse* fK0sPosMCResdEtaSdPhiS[kNCent];     //! K0s: Positive daughter: resolution for  delta(phi)* delta(eta)*
  THnSparse* fK0sNegMCResdEtaSdPhiS[kNCent];     //! K0s: Negative daughter: resolution for  delta(phi)* delta(eta)* 
  TH3F*   fK0sPosDaugFracShTPCclsTrig;           //! K0s Pos. Daug: mass vs pt vs fraction of shared TPC cls
  TH3F*   fK0sNegDaugFracShTPCclsTrig;           //! K0s Neg Daug: mass vs pt vs fraction of shared TPC cls
  TH3F*   fK0sDCADaugToPrimVtx;                  //! K0s: DCA to primary vertex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fK0sSpatialRes;                        //! K0s: Spatial resolution  
   
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
  TH3F*   fLambdaPosDaugFracShTPCcls;    //! Lambda Pos. Daug: mass vs pt vs fraction of shared TPC cls
  TH3F*   fLambdaNegDaugFracShTPCcls;    //! Lambda Neg Daug: mass vs pt vs fraction of shared TPC cls

  TH2F*   fLambdaDaughtersPt;                       //! Lambda: pt of daughters
  THnSparse* fLambdaPosDaugdPhiSdEtaS[kNCent];      //! Lambda:Positive daughter: delta(phi)* delta(eta)*    
  THnSparse* fLambdaNegDaugdPhiSdEtaS[kNCent];      //! Lambda: Negative daughter: delta(phi)* delta(eta)* 
  THnSparse* fLambdaPosDaugSplCheckCovMat[kNCent];  //! Lambda: Check Covariance Matrix elemenets between trigger trcak and daughter track
  THnSparse* fLambdaNegDaugSplCheckCovMat[kNCent];  //! Lambda: Check Covariance Matrix elemenets between trigger trcak and daughter track
  THnSparse* fLambdaPosMCResdEtaSdPhiS[kNCent];     //! Lambda: Positive daughter: resolution for  delta(phi)* delta(eta)*
  THnSparse* fLambdaNegMCResdEtaSdPhiS[kNCent];     //! Lambda: Negative daughter: resolution for  delta(phi)* delta(eta)* 
  TH3F*   fLambdaPosDaugFracShTPCclsTrig;           //! Lambda Pos. Daug: mass vs pt vs fraction of shared TPC cls
  TH3F*   fLambdaNegDaugFracShTPCclsTrig;           //! Lambda Neg Daug: mass vs pt vs fraction of shared TPC cls
  TH3F*   fLambdaDCADaugToPrimVtx;                  //! Lambda: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fLambdaSpatialRes;                        //! Lambda: Spatial resolution  

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
  TH3F*   fAntiLambdaMassPtPhi;                //! AntiLambda: mass vs phi 
  TH3F*   fAntiLambdaPosDaugFracShTPCcls;      //! AntiLambda Pos. Daug: mass vs pt vs fraction of shared TPC cls
  TH3F*   fAntiLambdaNegDaugFracShTPCcls;      //! AntiLambda Neg Daug: mass vs pt vs fraction of shared TPC cls

  TH2F*   fAntiLambdaDaughtersPt;                       //! AntiLambda: pt of daughters
  THnSparse* fAntiLambdaPosDaugdPhiSdEtaS[kNCent];      //! AntiLambda: Positive daughter: delta(phi)* delta(eta)*    
  THnSparse* fAntiLambdaNegDaugdPhiSdEtaS[kNCent];      //! AntiLambda: Negative daughter: delta(phi)* delta(eta)* 
  THnSparse* fAntiLambdaPosDaugSplCheckCovMat[kNCent];  //! AntiLambda: Check Covariance Matrix elemenets between trigger trcak and daughter track
  THnSparse* fAntiLambdaNegDaugSplCheckCovMat[kNCent];  //! AntiLambda: Check Covariance Matrix elemenets between trigger trcak and daughter track
  THnSparse* fAntiLambdaPosMCResdEtaSdPhiS[kNCent];     //! AntiLambda: Positive daughter: resolution for  delta(phi)* delta(eta)*
  THnSparse* fAntiLambdaNegMCResdEtaSdPhiS[kNCent];     //! AntiLambda: Negative daughter: resolution for  delta(phi)* delta(eta)* 
  TH3F*   fAntiLambdaPosDaugFracShTPCclsTrig;           //! AntiLambda Pos. Daug: mass vs pt vs fraction of shared TPC cls
  TH3F*   fAntiLambdaNegDaugFracShTPCclsTrig;           //! AntiLambda Neg Daug: mass vs pt vs fraction of shared TPC cls
  TH3F*   fAntiLambdaDCADaugToPrimVtx;                  //! AntiLambda: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fAntiLambdaSpatialRes;                        //! AntiLambda: Spatial resolution  

  TH3F*   fAntiLambdadPhidEtaMC[kNCent*kN1];            //! AntiLambda MC: Delta phi,Delta eta vs Z vertex position
  TH3F*   fAntiLambdadPhidEtaPtL[kNVtxZ*kNCent*kN1];    //! AntiLambda: Delta phi,Delta eta vs pt of the leading particle
  //TH3F*   fAntiLambdadPhidEtaPtLBckg[kNCent*kN1];     //! AntiLambda background: Delta phi,Delta eta vs Z vertex position

  TH2F*   fAntiLambdaBckgDecLength;            //! AntiLambda background: Decay lenght vs leading particle's pt inside a radio wrt the near-side peak
  TH3F*   fAntiLambdaBckgDCADaugToPrimVtx;     //! AntiLambda background: DCA to primary vrtex of daughters vs leading particle's pt inside a radio wrt the near-side peak
  TH2F*   fAntiLambdaBckgEtaPhi;               //! AntiLambda background: Phi vs Eta inside a radio wrt the near-side peak
  TH2F*   fAntiLambdaBckgPhiRadio ;            //! AntiLambda background: Phi vs radio inside a radio wrt the near-side peak
  TH2F*   fAntiLambdaBckgDCANegDaugToPrimVtx;  //! AntiLambda background: DCA of Negative daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  TH2F*   fAntiLambdaBckgDCAPosDaugToPrimVtx;  //! AntiLambda background: DCA of Positive daughter to the primary vertex inside the radio 0.4 wrt the near-side peak
  
    
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

  TH2F*   fK0sCTau;                          //! K0s: ctau 
  TH2F*   fK0sBckgCTau;                      //! K0s Bckg: ctau

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

  TH2F*   fLambdaCTau;                          //! Lambda: ctau 
  TH2F*   fLambdaBckgCTau;                      //! Lambda Bckg: ctau

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

  TH2F*   fAntiLambdaCTau;                          //! AntiLambda: ctau 
  TH2F*   fAntiLambdaBckgCTau;                      //! AntiLambda Bckg: ctau

  ///  ==== Mixed Events plots === ///
  TH3F*  fK0sdPhidEtaME[kNVtxZ*kNCent*kN1+1];             //! K0s Mixed Events
  TH3F*  fLambdadPhidEtaME[kNVtxZ*kNCent*kN1+1];          //! Lambda Mixed Events
  TH3F*  fAntiLambdadPhidEtaME[kNVtxZ*kNCent*kN1+1];      //! AntiLambda Mixed Events

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
