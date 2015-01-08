#ifndef ALIANALYSISTASKPROTONLAMBDA_H
#define ALIANALYSISTASKPROTONLAMBDA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Task to study femtoscopic proton-lambda correlations
// Author: Hans Beck, Hans.Beck@cern.ch

// Includes. All classes used in the .h file need to be included.
// Your task compiles with forward declaration and including in the .cxx
// file, but forward declaration will fail when loading this task from
// a library and just including this .h file.
/* #include <TH1F.h> */
/* #include <TH2F.h> */
/* #include <TH3F.h> */
/* #include <TList.h> */
/* #include <TAxis.h> */
class TH1F;
class TH2F;
class TH3F;
#include <THnSparse.h>
class TList;
class TAxis;

#include <AliAnalysisTaskSE.h>
/* #include <AliTPCPIDResponse.h> */
/* #include <AliPIDResponse.h> */
class AliTPCPIDResponse;
class AliPIDResponse;

#include <AliAODEvent.h>
#include <AliAODVertex.h>
#include <AliAODv0.h>
#include <AliAODTrack.h>

class AliAnalysisTaskProtonLambda : public AliAnalysisTaskSE {
  // Forward declaration for nested classes
 private:
  class FemtoBufferTrack;
  class FemtoBufferV0;
  class FemtoBufferEvent;
  class FemtoBuffer;
 public:
  AliAnalysisTaskProtonLambda();                 // Dummy constructor
  AliAnalysisTaskProtonLambda(const char *name); // Constructor
  virtual ~AliAnalysisTaskProtonLambda();                // Destructor
  AliAnalysisTaskProtonLambda(const AliAnalysisTaskProtonLambda& atpl); // Not implemented
  AliAnalysisTaskProtonLambda& operator=(const AliAnalysisTaskProtonLambda& atpl); //Not implemented
  
  void   UserCreateOutputObjects();  // Called once at the beginning
  void   UserExec(Option_t *option); // Main loop
  void   Terminate(Option_t *);      // Called once at the end
  
 private:
  void   ProcessOffline(AliAODv0 *v0, AliAODTrack *pTrack, AliAODTrack *nTrack);    //! Use offline V0 finder
  void   ProcessOnTheFly(AliAODv0 *v0, AliAODTrack *pTrack, AliAODTrack *nTrack);   //! Use on-the-fly V0 finder
  void   ProcessTOF(AliAODTrack *track);         // Use Tof for 1.0 GeV/c < p < 3.25 GeV/c
  void   ProcessTPC(AliAODTrack *track);         // Use TPC for p < 0.75 GeV/c
  void   ProcessHybrid(AliAODTrack *track);      // Preselection via TPC, identification w/ Tof 0.75 GeV/c < p < 1.0 GeV/c
  void   CleaningProcedure();         // Assure uniqueness of all daughters and primaries
  void   ProcessReal();               // TTR cut & qinv calc for real events
  void   ProcessMixed();              // TTR cut & qinv calc for mixed events
  void   ProcessRealBackground();     // Check procedure with lambda background from sideband, real events
  void   ProcessMixedBackground();    // Check procedure with lambda background from sideband, real events

  Float_t Qinv(FemtoBufferV0 v01, FemtoBufferV0 v02);           // Calcs Qinv for (anti-)Lam(anti-)Lam
  Float_t Qinv(FemtoBufferV0 v0, FemtoBufferTrack track);       // Calcs Qinv for (anti-)Lam(anti-)Pro
  Float_t Qinv(FemtoBufferTrack track, FemtoBufferV0 v0);       // Calls Qinv(v0, track)
  Float_t QinvProPro(FemtoBufferTrack proTrack1,
		     FemtoBufferTrack proTrack2);               // Calcs Qinv for proton proton
  Float_t QinvPioPro(FemtoBufferTrack pioTrack,
		     FemtoBufferTrack proTrack);                // Calcs Qinv for pion proton
  Float_t QinvConstr(FemtoBufferV0 v0, FemtoBufferTrack track); // w/ vtx constraint for pri pro
  Float_t QinvConstr(FemtoBufferTrack track, FemtoBufferV0 v0); // Calls QinvConstr(v0, track)
  Float_t Minv(FemtoBufferV0 v01, FemtoBufferV0 v02);           // Calcs Minv for (anti-)Lam(anti-)Lam
  Float_t Minv(FemtoBufferV0 v0, FemtoBufferTrack track);       // Calcs Minv for (anti-)Lam(anti-)Pro
  Float_t Minv(FemtoBufferTrack track, FemtoBufferV0 v0);       // Calls Minv(v0, track)
  Float_t mt(FemtoBufferTrack track, FemtoBufferV0 v0);         // mt of pair
  Float_t mt(FemtoBufferV0 v0, FemtoBufferTrack track);         //    ""
  Float_t mt(FemtoBufferV0 v01, FemtoBufferV0 v02);             //    ""
  static Float_t ktSquared(FemtoBufferV0 v01, FemtoBufferV0 v02);      // kt squared
  static Float_t ktSquared(FemtoBufferTrack track, FemtoBufferV0 v0);  //    ""
  static Float_t ktSquared(FemtoBufferV0 v0, FemtoBufferTrack track);  //    ""
  static Float_t calcDist(const Float_t r1[3], const Float_t r2[3]);   // Calcs the distance between vector r1 and r2
  //  static Float_t calcDistSft(const Float_t r1[3], const Float_t r2[3], 
  //			     const FemtoBufferEvent *evt1,
  //			     const FemtoBufferEvent *evt2);            // Calcs the distance between vector r1 and r2 shifting them first to originate both in (0,0,0)
  Bool_t  goodDCA(AliAODTrack *track);                  //! Does a cut on the DCAxy value and fills a histogram
  Float_t RapidityProton(AliAODTrack *track);           //! Rapidity assuming proton mass
  //  void    getTPConlyV0Info(const AliAODTrack *posDaughter,const AliAODTrack *negDaughter, Double_t tpcV0Mom[3], Double_t TPConlyV0MinvLam, Double_t TPConlyV0MinvALam);   //! Get V0 momentum and invariant mass using TPC only tracks

  // For AOD tracks, need a function to get the dca
  static Float_t DCAxy(const AliAODTrack *track, const AliVEvent *evt); 
  void StoreGlobalTrackReference(AliAODTrack *track);      // Store pointer to global track
  void ResetGlobalTrackReference();                        // Reset all pointers to global tracks
  Bool_t acceptTrack(const AliAODTrack *track);            // Famous crossed rows / findable clusters
  Bool_t GoodTPCFitMapSharedMap(const AliAODTrack *track); // Rejects shared clusters, primaries
  Bool_t GoodTPCFitMapSharedMap(const AliAODTrack *pTrack,
				const AliAODTrack *nTrack);// Rejects shared clusters, two V0 daughters
  Float_t GetCorrectedTOFSignal(const AliVTrack *track);   // Return correct TOF signal for old and new AODs
  void FillDedxHist(const AliVTrack *track);               // Fill dE/dx histograms

  const Bool_t   fkUseOnTheFly;            //! Use on-the-fly or offline V0 finder
  const Float_t  fkAbsZvertexCut;          //! Values of Zvertex and centrality cut
  const Float_t  fkCentCut;                //! are needed for event mixing
  const Float_t  fkLamMass;                //! PDG-Mass of the lambda particle
  const Float_t  fkProMass;                //! PDG-Mass of the lambda particle
  const Float_t  fkPioMass;                //! PDG-Mass of the lambda particle

  Double_t       fPrimaryVtxPosition[3];   //! Position of prim. vertex

  AliPIDResponse  *fPIDResponse;           //! PID response object
  AliTPCPIDResponse *fTpcResponse;         //! Bethe-Bloch parametrisation


  FemtoBuffer     *fFemtoBuffer;           //! Event mixing: event collection
  
  AliAODEvent     *fAOD;                   //! AOD event
  AliAODVertex    *fPrimaryVtx;            //! AOD vertex

  TList           *fOutputList;            //! V0 output list 
  TList           *fOutputPrimaries;       //! Primaries output list
  TList           *fOutput2Part;           //! Two-particle output list

  // Store pointers to global tracks for pid and dca
  AliAODTrack     **fGTI;                  //! Array of pointers, just nicely sorted according to the id
  const UShort_t  fTrackBuffSize;          //! Size fo the above array, ~12000 for PbPb
  //
  //! ---------        Histograms        ---------
  //
  TH1F        *fHistGoodEvent;                  //! Control hist for event cuts

  // Low memory version: just get rid of a lot of histograms
  /* TH2F        *fHistPrimaryVertexPosXY;         //! Primary vertex position in XY */
  /* TH1F        *fHistPrimaryVertexPosZ;          //! Primary vertex position in Z */
  /* TH1F        *fHistTrackMultiplicity;          //! Track multiplicity distribution */

  /* TH1F        *fHistShareV0pos;                 //! Number of shared clusters, V0 pos daughter */
  /* TH1F        *fHistShareV0neg;                 //! Number of shared clusters, V0 neg daughter */

  /* TH2F        *fHistPosTofBeforeCut;            //! Tof signal for pos. daughter w/o cut */
  /* TH2F        *fHistPosTofAfterCut;             //! Tof signal for pos. daughter with cut */
  /* TH2F        *fHistNegTofBeforeCut;            //! Tof signal for pos. daughter w/o cut */
  /* TH2F        *fHistNegTofAfterCut;             //! Tof signal for pos. daughter with cut */
  /* TH2F        *fHistPosTpcBeforeCut;            //! Tpc signal for pos. daughter w/o cut */
  /* TH2F        *fHistPosTpcAfterCut;             //! Tpc signal for pos. daughter with cut */
  /* TH2F        *fHistNegTpcBeforeCut;            //! Tpc signal for pos. daughter w/o cut */
  /* TH2F        *fHistNegTpcAfterCut;             //! Tpc signal for pos. daughter with cut */
  
  /* TH2F        *fHistGoodV0;                     //! Control histogram for cuts */
  /* TH2F        *fHistCorrectSigns;               //! Check of sign of tracks */
  /* TH2F        *fHistDcaPosToPrimVertex;         //! DCA of pos daughter to primary vertex */
  /* TH2F        *fHistDcaNegToPrimVertex;         //! DCA of neg daughter to primary vertex */
  /* TH2F        *fHistDcaPosToPrimVertexZoom;     //! Zoom */
  /* TH2F        *fHistDcaNegToPrimVertexZoom;     //! Zoom */
  /* TH2F        *fHistRadiusV0;                   //! Radial position of V0 vertex */
  /* TH2F        *fHistDecayLengthV0;              //! Dist prim. vertex - V0 vertex */
  /* TH2F        *fHistDcaV0Daughters;             //! DCA of daughters */
  /* TH2F        *fHistChi2;                       //! Chi2 for V0 vertex */
  /* TH2F        *fHistCosPointAngle;              //! Cosine of pointing angle */
  /* TH2F        *fHistCosPointAngleZoom;          //! Zoom */

  TH1F        *fHistSideBandOffLam;             //! Off: side band background lambda
  TH1F        *fHistSideBandOffALam;            //! Off: side band background anti-lambda

  TH2F        *fHistTPCNclsPosOffLam;           //! Off: lam-mass vs TPCNcls(pos daughter)
  TH2F        *fHistTPCNclsNegOffLam;           //! Off: lam-mass vs TPCNcls(neg daughter)
  TH2F        *fHistTPCNclsPosOffALam;          //! Off: alam-mass vs TPCNcls(pos daughter)
  TH2F        *fHistTPCNclsNegOffALam;          //! Off: alam-mass vs TPCNcls(neg daughter)
  /* TH2F        *fHistPosNsigmaTpcOffLam;         //! */
  /* TH2F        *fHistPosNsigmaTpcOffALam;        //! */
  /* TH2F        *fHistNegNsigmaTpcOffLam;         //! */
  /* TH2F        *fHistNegNsigmaTpcOffALam;        //! */
  /* TH2F        *fHistUseTofOffLam;               //! */
  /* TH2F        *fHistUseTofOffALam;              //! */
  /* TH2F        *fHistDcaPosOffLam;               //! Off: lam-mass vs DCA(pos - prim. vertex) */
  /* TH2F        *fHistDcaPosOffALam;              //! Off: alam-mass vs DCA(pos - prim. vertex) */
  /* TH2F        *fHistDcaNegOffLam;               //! Off: lam-mass vs DCA(neg - prim. vertex) */
  /* TH2F        *fHistDcaNegOffALam;              //! Off: alam-mass vs DCA(neg - prim. vertex) */
  /* TH2F        *fHistDcaV0DaughtersOffLam;       //! Off: lam-mass vs DCA(V0-daughters) */
  /* TH2F        *fHistDcaV0DaughtersOffALam;      //! Off: alam-mass vs DCA(V0-daughters) */
  /* TH2F        *fHistCosPointLamOff;             //! Off: lam-mass vs cosine(pointing angle) */
  /* TH2F        *fHistCosPointALamOff;            //! Off: Alam-mass vs cosine(pointing angle) */
  /* TH2F        *fHistCosPointLamZoomOff;         //! Zoom */
  /* TH2F        *fHistCosPointALamZoomOff;        //! Zoom */
  /* TH2F        *fHistV0RadiusLamOff;             //! Off: lam-mass vs radial position(V0) */
  /* TH2F        *fHistV0RadiusALamOff;            //! Off: alam-mass vs radial position(V0) */
  /* TH2F        *fHistV0DecayLengthLamOff;        //! Off: lam-mass vs decay length(V0) */
  /* TH2F        *fHistV0DecayLengthALamOff;       //! Off: alam-mass vs decay length(V0) */
  /* TH2F        *fHistDcaV0PriVertexLamOff;       //! Off: lam-mass vs dca (prim. vertex - V0) */
  /* TH2F        *fHistDcaV0PriVertexALamOff;      //! Off: alam-mass vs dca (prim. vertex - V0) */
  TH1F        *fHistMassLambdaOff;              //! Off: invariant mass assuming lambda
  TH1F        *fHistMassAntiLambdaOff;          //! Off: invariant mass assuming anti-lambda
  /* TH2F        *fHistPtVsMassLambdaOff;          //! Off: lam-mass vs pt(V0) */
  /* TH2F        *fHistPtVsMassAntiLambdaOff;      //! Off: alam-mass vs pt(V0) */
  TH3F        *fHistYPtMassLamOff;              //! Off: lam-mass vs y-pt
  TH3F        *fHistYPtMassALamOff;             //! Off: alam-mass vs y-pt
  /* TH2F        *fHistPtVsYLambdaOff;             //! Off: y-pt (lambda) */
  /* TH2F        *fHistPtVsYAntiLambdaOff;         //! Off: y-pt (anti-lambda) */
  
  TH1F        *fHistSideBandOnLam;              //! On-the-fly: side band background lambda
  TH1F        *fHistSideBandOnALam;             //! On-the-fly: side band background anti-lambda

  /* TH2F        *fHistLikeSignOnLam;              //! On-the-fly: lambda mass vs signs of pairs */
  /* TH2F        *fHistLikeSignOnALam;             //! On-the-fly: antilambda mass vs signs of pairs */
  TH2F        *fHistTPCNclsPosOnLam;            //! On: lam-mass vs TPCNcls(pos daughter)
  TH2F        *fHistTPCNclsNegOnLam;            //! On: lam-mass vs TPCNcls(neg daughter)
  TH2F        *fHistTPCNclsPosOnALam;           //! On: alam-mass vs TPCNcls(pos daughter)
  TH2F        *fHistTPCNclsNegOnALam;           //! On: alam-mass vs TPCNcls(neg daughter)
  /* TH2F        *fHistPosNsigmaTpcOnLam;          //! */
  /* TH2F        *fHistPosNsigmaTpcOnALam;         //! */
  /* TH2F        *fHistNegNsigmaTpcOnLam;          //! */
  /* TH2F        *fHistNegNsigmaTpcOnALam;         //!  */
  /* TH2F        *fHistUseTofOnLam;                //! */
  /* TH2F        *fHistUseTofOnALam;               //! */
  /* TH2F        *fHistDcaPosOnLam;                //! On: lam-mass vs DCA(pos - prim. vertex) */
  /* TH2F        *fHistDcaPosOnALam;               //! On: alam-mass vs DCA(pos - prim. vertex) */
  /* TH2F        *fHistDcaNegOnLam;                //! On: lam-mass vs DCA(neg - prim. vertex) */
  /* TH2F        *fHistDcaNegOnALam;               //! On: alam-mass vs DCA(neg - prim. vertex) */
  /* TH2F        *fHistDcaV0DaughtersOnLam;        //! On: lam-mass vs DCA(V0-daughters) */
  /* TH2F        *fHistDcaV0DaughtersOnALam;       //! On: alam-mass vs DCA(V0-daughters) */
  /* TH2F        *fHistCosPointLamOn;              //! On: lam-mass vs cosine(pointing angle) */
  /* TH2F        *fHistCosPointALamOn;             //! On: Alam-mass vs cosine(pointing angle) */
  /* TH2F        *fHistCosPointLamZoomOn;          //! Zoom */
  /* TH2F        *fHistCosPointALamZoomOn;         //! Zoom */
  /* TH2F        *fHistV0RadiusLamOn;              //! On: lam-mass vs radial position(V0) */
  /* TH2F        *fHistV0RadiusALamOn;             //! On: alam-mass vs radial position(V0) */
  /* TH2F        *fHistV0DecayLengthLamOn;         //! On: lam-mass vs decay length(V0) */
  /* TH2F        *fHistV0DecayLengthALamOn;        //! On: alam-mass vs decay length(V0) */
  /* TH2F        *fHistDcaV0PriVertexLamOn;        //! On: lam-mass vs dca (prim. vertex - V0) */
  /* TH2F        *fHistDcaV0PriVertexALamOn;       //! On: alam-mass vs dca (prim. vertex - V0) */
  /* TH2F        *fHistChi2TPCPosLamOn;            //! On: Chi2 TPC clusters */
  /* TH2F        *fHistChi2TPCPosALamOn;           //! */
  /* TH2F        *fHistChi2TPCNegLamOn;            //! */
  /* TH2F        *fHistChi2TPCNegALamOn;           //! */
  /* TH1F        *fHistMinvTPConlyLamOn;           //! On: minv w/ TPC only daughters, lambda */
  /* TH1F        *fHistMinvTPConlyALamOn;          //! On:        --""--             , anti-lambda */
  TH1F        *fHistMassLambdaOn;               //! On: invariant mass assuming lambda
  TH1F        *fHistMassAntiLambdaOn;           //! On: invariant mass assuming anti-lambda
  /* TH2F        *fHistPtVsMassLambdaOn;           //! On: lam-mass vs pt(V0) */
  /* TH2F        *fHistPtVsMassAntiLambdaOn;       //! On: alam-mass vs pt(V0) */
  TH3F        *fHistYPtMassLamOn;               //! On: lam-mass vs y-pt
  TH3F        *fHistYPtMassALamOn;              //! On: alam-mass vs y-pt
  /* TH2F        *fHistPtVsYLambdaOn;              //! On: y-pt (lambda) */
  /* TH2F        *fHistPtVsYAntiLambdaOn;          //! On: y-pt (anti-lambda) */

  /* TH3F        *fHistMomDiffLam;                 //! Momentum difference standard V0 / */
  /* TH3F        *fHistMomDiffALam;                //! TPC only V0 */
  /* TH3F        *fHistMomDiffBgLam;               //! */
  /* TH3F        *fHistMomDiffBgALam;              //! */
  /* TH3F        *fHistMomDiffWoSPDLam;            //! Momentum difference standard V0 / */
  /* TH3F        *fHistMomDiffWoSPDALam;           //! TPC only V0, only tracks without  */
  /* TH3F        *fHistMomDiffWoSPDBgLam;          //! SPD hits */
  /* TH3F        *fHistMomDiffWoSPDBgALam;         //! */

  // Primary particles
  TH1F        *fPriHistShare;                   //! Primaries: number of shared clusters
  // ToF histograms
  /* TH1F        *fPriHistPosNsigmaTof;            //! Primaries: Nsigma TOF distribution (pos) */
  TH2F        *fPriHistPosNsigmaTofVsP;         //! Pri: Nsigma TOF vs total momentum (pos)
  TH2F        *fPriHistPosNsigmaTofVsPt;        //! Pri: Nsigma TOF vs transverse momentum (pos)
  /* TH1F        *fPriHistNegNsigmaTof;            //! Primaries: Nsigma TOF distribution (neg) */
  TH2F        *fPriHistNegNsigmaTofVsP;         //! Pri: Nsigma TOF vs total momentum (neg)
  TH2F        *fPriHistNegNsigmaTofVsPt;        //! Pri: Nsigma TOF vs transverse momentum (neg)
  TH2F        *fPriHistTOFsignalPosVsP;         //! Pri: TOF signal vs p (pos)
  TH2F        *fPriHistTOFsignalPosVsPt;        //! Pri: TOF signal vs pt (pos)
  TH2F        *fPriHistTOFsignalNegVsP;         //! Pri: TOF signal vs p (neg)
  TH2F        *fPriHistTOFsignalNegVsPt;        //! Pri: TOF signal vs pt (neg)
  // Hybrid analysis: dEdx & ToF histograms
  TH1F        *fPriHistHybridTOFsigPosWoTPC;    //! Pri: TOF signal without dEdx selection (pos)
  TH1F        *fPriHistHybridTOFsigPosTPCok;    //! Pri: TOF signal with dEdx selection (pos)
  TH1F        *fPriHistHybridTOFsigNegWoTPC;    //! Pri: TOF signal without dEdx selection (neg)
  TH1F        *fPriHistHybridTOFsigNegTPCok;    //! Pri: TOF signal with dEdx selection (neg)

  /* TH1F        *fPriHistHasTofPos;               //! Pri: 1 = TOFpid bit is there, 0 else */
  TH2F        *fPriHistTPCsignalPos;            //! Pri: TPC dE/dx signal vs p
  /* TH2F        *fPriHistNsigmaTPCPos;            //! Pri: Nsigma TPC distribution vs p */
  /* TH2F        *fPriHistTPCsignalTOFcutPos;      //! Pri: dE/dx vs p + TOF sigma(proton) > -10 */
  /* TH2F        *fPriHistNsigmaTPCTOFcutPos;      //! Pri: Nsigma TPC vs p + TOF sigma(proton) > -10 */
  TH2F        *fPriHistTPCsignalLowPPos;        //! Pri: dEdx for 0.1 < p < 0.3
  TH2F        *fPriHistTPCsignalMedPPos;        //! Pri: dEdx for 0.3 < p < 0.9
  TH2F        *fPriHistTPCsignalHigPPos;        //! Pri: dEdx for 0.9 < p < 1.9

  /* TH1F        *fPriHistHasTofNeg;               //! Pri: 1 = TOFpid bit is there, 0 else */
  TH2F        *fPriHistTPCsignalNeg;            //! Pri: TPC dE/dx signal vs p
  /* TH2F        *fPriHistNsigmaTPCNeg;            //! Pri: Nsigma TPC distribution vs p */
  /* TH2F        *fPriHistTPCsignalTOFcutNeg;      //! Pri:   " + TOF sigma(proton) > -10 */
  /* TH2F        *fPriHistNsigmaTPCTOFcutNeg;      //! Pri: Nsigma TPC vs p + TOF sigma(proton) > -10 */
  TH2F        *fPriHistTPCsignalLowPNeg;        //! Pri: dEdx for 0.1 < p < 0.3
  TH2F        *fPriHistTPCsignalMedPNeg;        //! Pri: dEdx for 0.3 < p < 0.9
  TH2F        *fPriHistTPCsignalHigPNeg;        //! Pri: dEdx for 0.9 < p < 1.9

  TH3F        *fPriHistDCAxyYPtPro;             //! Pri: DCAxy vs (y,pt) for protons
  TH3F        *fPriHistDCAxyYPtAPro;            //! Pri: DCAxy vs (y,pt) for anti-protons

  /*        // Histograms for two particle observables */
  /* TH2F        *f2HistLamLamMeanMinDistProReal;     //! TTR: Real events */
  /* TH2F        *f2HistLamLamMeanMinDistPioReal;     //!  */
  /* TH2F        *f2HistLamProMeanMinDistProReal;     //!  */
  /* TH2F        *f2HistALamALamMeanMinDistAProReal;  //!  */
  /* TH2F        *f2HistALamALamMeanMinDistPioReal;   //!  */
  /* TH2F        *f2HistALamAProMeanMinDistAProReal;  //!  */
  /* TH2F        *f2HistSftLamLamMeanMinDistProReal;     //! TTR: Shifted pri vtx */
  /* TH2F        *f2HistSftLamLamMeanMinDistPioReal;     //!  */
  /* TH2F        *f2HistSftLamProMeanMinDistProReal;     //!  */
  /* TH2F        *f2HistSftALamALamMeanMinDistAProReal;  //!  */
  /* TH2F        *f2HistSftALamALamMeanMinDistPioReal;   //!  */
  /* TH2F        *f2HistSftALamAProMeanMinDistAProReal;  //!  */
  /* TH2F        *f2HistSftIrocLamLamMeanMinDistProReal;     //! TTR: Shifted pri vtx: IROCs */
  /* TH2F        *f2HistSftIrocLamLamMeanMinDistPioReal;     //!  */
  /* TH2F        *f2HistSftIrocLamProMeanMinDistProReal;     //!  */
  /* TH2F        *f2HistSftIrocALamALamMeanMinDistAProReal;  //!  */
  /* TH2F        *f2HistSftIrocALamALamMeanMinDistPioReal;   //!  */
  /* TH2F        *f2HistSftIrocALamAProMeanMinDistAProReal;  //!  */
  /* TH2F        *f2HistSftOrocLamLamMeanMinDistProReal;     //! TTR: Shifted pri vtx: OROCs */
  /* TH2F        *f2HistSftOrocLamLamMeanMinDistPioReal;     //!  */
  /* TH2F        *f2HistSftOrocLamProMeanMinDistProReal;     //!  */
  /* TH2F        *f2HistSftOrocALamALamMeanMinDistAProReal;  //!  */
  /* TH2F        *f2HistSftOrocALamALamMeanMinDistPioReal;   //!  */
  /* TH2F        *f2HistSftOrocALamAProMeanMinDistAProReal;  //!  */
    
  // transverse mass mt of the pair
  /* TH1F        *f2HistMtLamLamReal;     //! */
  TH1F        *f2HistMtLamProReal;     //!
  /* TH1F        *f2HistMtALamALamReal;   //! */
  TH1F        *f2HistMtALamAProReal;   //!
  // .. for low q pairs only
  /* TH1F        *f2HistMtLowQLamLamReal;     //! */
  TH1F        *f2HistMtLowQLamProReal;     //!
  /* TH1F        *f2HistMtLowQALamALamReal;   //! */
  TH1F        *f2HistMtLowQALamAProReal;   //!

  // September '12: use a THnSparseF with 4 dimensions
  // for pro lam: qinv propro, mean d propro
  // min d propro, qinv lampro
  THnSparseF       *LamProReal;  //!
  THnSparseF       *ALamAProReal;//!

  // Changed on march 4th: all correlation functions vs 
  // distances to try out diff. cuts
  /* TH3F        *f3HistLamLamQinvReal;                //! Qinv distribution real events */
  /* TH3F        *f3HistALamALamQinvReal;              //!                 */
  /* TH3F        *f3HistLamLamMinvReal;                //! Minv distribution real events */
  /* TH3F        *f3HistLamProMinvReal;                //! */
  /* TH3F        *f3HistALamALamMinvReal;              //!                 */
  /* TH3F        *f3HistALamAProMinvReal;              //! */

 // Background lambdas: real events
  /* TH2F        *f2HistBgLamBgLamMeanMinDistProReal;   //! TTR: Real events */
  /* TH2F        *f2HistBgLamBgLamMeanMinDistPioReal;   //!  */
  /* TH2F        *f2HistBgLamProMeanMinDistProReal;     //!  */
  /* TH2F        *f2HistBgALamBgALamMeanMinDistAProReal;//!  */
  /* TH2F        *f2HistBgALamBgALamMeanMinDistPioReal; //!  */
  /* TH2F        *f2HistBgALamAProMeanMinDistAProReal;  //!  */
  /* TH2F        *f2HistSftBgLamBgLamMeanMinDistProReal;   //! TTR: Shifted pri vtx */
  /* TH2F        *f2HistSftBgLamBgLamMeanMinDistPioReal;   //!  */
  /* TH2F        *f2HistSftBgLamProMeanMinDistProReal;     //!  */
  /* TH2F        *f2HistSftBgALamBgALamMeanMinDistAProReal;//!  */
  /* TH2F        *f2HistSftBgALamBgALamMeanMinDistPioReal; //!  */
  /* TH2F        *f2HistSftBgALamAProMeanMinDistAProReal;  //!  */
  /* TH2F        *f2HistSftIrocBgLamBgLamMeanMinDistProReal;   //! TTR: Shifted pri vtx: IROCs */
  /* TH2F        *f2HistSftIrocBgLamBgLamMeanMinDistPioReal;   //!  */
  /* TH2F        *f2HistSftIrocBgLamProMeanMinDistProReal;     //!  */
  /* TH2F        *f2HistSftIrocBgALamBgALamMeanMinDistAProReal;//!  */
  /* TH2F        *f2HistSftIrocBgALamBgALamMeanMinDistPioReal; //!  */
  /* TH2F        *f2HistSftIrocBgALamAProMeanMinDistAProReal;  //!  */
  /* TH2F        *f2HistSftOrocBgLamBgLamMeanMinDistProReal;   //! TTR: Shifted pri vtx: OROCs */
  /* TH2F        *f2HistSftOrocBgLamBgLamMeanMinDistPioReal;   //!  */
  /* TH2F        *f2HistSftOrocBgLamProMeanMinDistProReal;     //!  */
  /* TH2F        *f2HistSftOrocBgALamBgALamMeanMinDistAProReal;//!  */
  /* TH2F        *f2HistSftOrocBgALamBgALamMeanMinDistPioReal; //!  */
  /* TH2F        *f2HistSftOrocBgALamAProMeanMinDistAProReal;  //!  */

  // September '12: use a THnSparseF with 4 dimensions
  // for pro lam: qinv propro, mean d propro
  // min d propro, qinv lampro
  THnSparseF       *BgLamProReal;  //!
  THnSparseF       *BgALamAProReal;//!

  /* TH3F        *f3HistBgLamBgLamQinvReal;             //! Qinv distribution real events */
  /* TH3F        *f3HistBgALamBgALamQinvReal;           //!                 */
 
  // Mixed events
  /* TH2F        *f2HistLamLamMeanMinDistProMixed;      //! TTR: Mixed events */
  /* TH2F        *f2HistLamLamMeanMinDistPioMixed;      //!  */
  /* TH2F        *f2HistLamProMeanMinDistProMixed;      //!  */
  /* TH2F        *f2HistALamALamMeanMinDistAProMixed;   //!  */
  /* TH2F        *f2HistALamALamMeanMinDistPioMixed;    //!  */
  /* TH2F        *f2HistALamAProMeanMinDistAProMixed;   //!  */
  /* TH2F        *f2HistSftLamLamMeanMinDistProMixed;      //! TTR: Mixed events, shifted pri vtx */
  /* TH2F        *f2HistSftLamLamMeanMinDistPioMixed;      //!  */
  /* TH2F        *f2HistSftLamProMeanMinDistProMixed;      //!  */
  /* TH2F        *f2HistSftALamALamMeanMinDistAProMixed;   //!  */
  /* TH2F        *f2HistSftALamALamMeanMinDistPioMixed;    //!  */
  /* TH2F        *f2HistSftALamAProMeanMinDistAProMixed;   //!  */
  /* TH2F        *f2HistSftIrocLamLamMeanMinDistProMixed;      //! TTR: Mixed events, shifted pri vtx, IROCs */
  /* TH2F        *f2HistSftIrocLamLamMeanMinDistPioMixed;      //!  */
  /* TH2F        *f2HistSftIrocLamProMeanMinDistProMixed;      //!  */
  /* TH2F        *f2HistSftIrocALamALamMeanMinDistAProMixed;   //!  */
  /* TH2F        *f2HistSftIrocALamALamMeanMinDistPioMixed;    //!  */
  /* TH2F        *f2HistSftIrocALamAProMeanMinDistAProMixed;   //!  */
  /* TH2F        *f2HistSftOrocLamLamMeanMinDistProMixed;      //! TTR: Mixed events, shifted pri vtx, OROCs */
  /* TH2F        *f2HistSftOrocLamLamMeanMinDistPioMixed;      //!  */
  /* TH2F        *f2HistSftOrocLamProMeanMinDistProMixed;      //!  */
  /* TH2F        *f2HistSftOrocALamALamMeanMinDistAProMixed;   //!  */
  /* TH2F        *f2HistSftOrocALamALamMeanMinDistPioMixed;    //!  */
  /* TH2F        *f2HistSftOrocALamAProMeanMinDistAProMixed;   //!  */

  // September '12: use a THnSparseF with 4 dimensions
  // for pro lam: qinv propro, mean d propro
  // min d propro, qinv lampro
  THnSparseF       *LamProMixed; //!
  THnSparseF       *ALamAProMixed;//!

  /* TH3F        *f3HistLamLamQinvMixed;                //! Qinv distribution mixed events */
  /* TH3F        *f3HistALamALamQinvMixed;              //!                 */
  /* TH3F        *f3HistLamLamMinvMixed;                //! Minv distribution mixed events */
  /* TH3F        *f3HistLamProMinvMixed;                //! */
  /* TH3F        *f3HistALamALamMinvMixed;              //!                 */
  /* TH3F        *f3HistALamAProMinvMixed;              //! */
   // Background lambdas: mixed events
  /* TH2F        *f2HistBgLamBgLamMeanMinDistProMixed;   //! TTR: Mixed events */
  /* TH2F        *f2HistBgLamBgLamMeanMinDistPioMixed;   //!  */
  /* TH2F        *f2HistBgLamProMeanMinDistProMixed;     //!  */
  /* TH2F        *f2HistBgALamBgALamMeanMinDistAProMixed;//!  */
  /* TH2F        *f2HistBgALamBgALamMeanMinDistPioMixed; //!  */
  /* TH2F        *f2HistBgALamAProMeanMinDistAProMixed;  //!  */
  /* TH2F        *f2HistSftBgLamBgLamMeanMinDistProMixed;   //! TTR: Mixed events, shifted pri vtx */
  /* TH2F        *f2HistSftBgLamBgLamMeanMinDistPioMixed;   //!  */
  /* TH2F        *f2HistSftBgLamProMeanMinDistProMixed;     //!  */
  /* TH2F        *f2HistSftBgALamBgALamMeanMinDistAProMixed;//!  */
  /* TH2F        *f2HistSftBgALamBgALamMeanMinDistPioMixed; //!  */
  /* TH2F        *f2HistSftBgALamAProMeanMinDistAProMixed;  //!  */
  /* TH2F        *f2HistSftIrocBgLamBgLamMeanMinDistProMixed;   //! TTR: Mixed events, shifted pri vtx, IROCs */
  /* TH2F        *f2HistSftIrocBgLamBgLamMeanMinDistPioMixed;   //!  */
  /* TH2F        *f2HistSftIrocBgLamProMeanMinDistProMixed;     //!  */
  /* TH2F        *f2HistSftIrocBgALamBgALamMeanMinDistAProMixed;//!  */
  /* TH2F        *f2HistSftIrocBgALamBgALamMeanMinDistPioMixed; //!  */
  /* TH2F        *f2HistSftIrocBgALamAProMeanMinDistAProMixed;  //!  */
  /* TH2F        *f2HistSftOrocBgLamBgLamMeanMinDistProMixed;   //! TTR: Mixed events, shifted pri vtx, OROCs */
  /* TH2F        *f2HistSftOrocBgLamBgLamMeanMinDistPioMixed;   //!  */
  /* TH2F        *f2HistSftOrocBgLamProMeanMinDistProMixed;     //!  */
  /* TH2F        *f2HistSftOrocBgALamBgALamMeanMinDistAProMixed;//!  */
  /* TH2F        *f2HistSftOrocBgALamBgALamMeanMinDistPioMixed; //!  */
  /* TH2F        *f2HistSftOrocBgALamAProMeanMinDistAProMixed;  //!  */

  // September '12: use a THnSparseF with 4 dimensions
  // for pro lam: qinv propro, mean d propro
  // min d propro, qinv lampro
  THnSparseF       *BgLamProMixed;   //!
  THnSparseF       *BgALamAProMixed; //!

  /* TH3F        *f3HistBgLamBgLamQinvMixed;             //! Qinv distribution mixed events */
  /* TH3F        *f3HistBgALamBgALamQinvMixed;           //!                 */

 private:
  //  ------------------------------------
  //
  //    Nested classes for event mixing
  //    Features are 
  //    * low memory usage as only needed variables
  //      are stored and not, e.g., the full AliAODTrack
  //    * it's fast, as memory allocation occours
  //      only once and no deep copying is done
  //      when doing the fifo shift and reseting
  //      the event. Resetting and shifting could
  //      be reduced to only a few integer assignments.
  //      
  //  ------------------------------------


  class FemtoBufferTrack{// Charged tracks
  public:
    FemtoBufferTrack();                               // Constructor
    FemtoBufferTrack(const AliAODTrack *track,
		     const Float_t bfield,
		     const Float_t priVtx[3]);           // Constructor
    ~FemtoBufferTrack(){;}                            // Destructor, nothing to do.
    FemtoBufferTrack(const FemtoBufferTrack& fbt);    // Copy constructor
    FemtoBufferTrack& operator=(const FemtoBufferTrack& fbt); // Assignment operator
    void Set(const AliAODTrack *track,                // Set parameters of the FemtoBufferTrack
	     const Float_t bfield,                    // to these of the given AliAODtrack
	     const Float_t priVtx[3]);    
    void Set(const AliAODTrack *track,                // Overloaded fct that just
	     const Float_t bfield,	              // calls Set(..,Float)
	     const Double_t priVtx[3]);
    // Two functions to get the global and shifted positions. Both is not trivial as
    // tracks in ALICE get propagated in the local system
    // Shifted for studying all events shifted to (0,0,0)
    void GetShiftedPositionAtShiftedRadii(const AliAODTrack *track,
					  const Float_t bfield, const Float_t priVtx[3]); 
    // Global if not shifted
    void GetGlobalPositionAtGlobalRadii(const AliAODTrack *track,
					const Float_t bfield);
    Bool_t UseIt() const {return fID!=65535;}          // Use entry? eg. set to false by cleaning procedure
    void   SetBadFlag(){fID=65535;}                   // Can only set 'bad track' flag
    // Interpretation of ALICE coding rule RC14 is that FemtoBufferTrack is a private member
    // of AliAnalysisTaskProtonLambda, so the following data member are private
    UShort_t fID;               //! Unique track id (->AliAODTrack.h), UShort_t goes to 65000
    Double_t fP[3];             //! Momentum of track
    /* Float_t  fXglobal[9][3];    //! Global positions at different global radii */
    Float_t  fXshifted[9][3];   //! Shifted positions at different shifted radii
  };
  
  class FemtoBufferV0{// V0 vertices
  public:
    // Constructors n stuff
    FemtoBufferV0();                          // Constructor
    FemtoBufferV0(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter,
		  const Double_t bfield, Double_t priVtxPos[3]);   // Constructor
    ~FemtoBufferV0(){;}                       // Destructor, nothing to do.
    FemtoBufferV0(const FemtoBufferV0 &fbv);  // Copy constructor
    FemtoBufferV0& operator=(const FemtoBufferV0 &fbv);// Assignment operator
    // Functions
    void Set(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter
	     ,const Double_t bfield, Double_t priVtxPos[3]); // Set properties of this to these of v0
    Bool_t UseIt() const {return fCosPoint>0;} // Use entry? eg. set to false by cleaning procedure
    void SetBadFlag(){fCosPoint=-9999.;fPosDaughter.SetBadFlag();
      fNegDaughter.SetBadFlag();}               // Can only set 'bad track' flag
    // Data member
    Double_t fP[3];                 //! Momentum (x,y,z)
    Float_t fCosPoint;              //! Cosine of pointing angle 
    FemtoBufferTrack fPosDaughter;  //! Positive daughter of the V0
    FemtoBufferTrack fNegDaughter;  //! Negative daughter fo the V0
  };
  
  class FemtoBufferEvent{// Event 
  public:
    FemtoBufferEvent();                       // Constructor
    FemtoBufferEvent(const UShort_t priTrackBuff,const UShort_t V0Buff, 
		     const Double_t bfield,const Double_t priVtxPos[3]); // Constructor
    FemtoBufferEvent(const UShort_t priTrackBuff,const UShort_t V0Buff); // Constructor
    FemtoBufferEvent(const FemtoBufferEvent &fbe);           // Copy constructor
    // Assignment operator won't change the size of the arrays!
    FemtoBufferEvent& operator=(const FemtoBufferEvent &fbe); // Assignment operator
    ~FemtoBufferEvent();                       // Destructor
    void Reset(const Double_t bfield, const Double_t priVtxPos[3]);// Resets the event with new variables given
    
    // Functions to add particles to the event
    void AddPro(const AliAODTrack *track);      // Add a proton to this event  
    void AddAPro(const AliAODTrack *track);     // Add a anti-proton this actual event  
    void AddLam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter); 
    // Add a lamba to this event
    void AddALam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter);
    // Add a anti-lambda to this event  
    void AddBgLam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter); 
    // Add lambda background to this event
    void AddBgALam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter); 
    // Add anti-lambda background to this event
    
    // Getters for the event properties, no setters as you're not supposed to change them :)
    // The fBfield and fPriVtxPos get set on each 'reset' for a new event and the limits get set 
    // on creation of the event
    Double_t GetBfield()const{return fBfield;}            // Getter for magnetic field of event
    void GetVtxPos(Double_t xyz[3])const{for(Int_t i=0;i<3;i++)xyz[i]=fPriVtxPos[i];} // Get the xyz of the vertex.
    UShort_t GetPriTrackLim()const{return fPriTrackLim;}  // Get the size of the array for primary tracks
    UShort_t GetV0Lim()const{return fV0Lim;}              // Get the size of the array for V0s
    // The number of tracks stored in the event
    UShort_t GetNPro()const{return fNProTracks;}       // Number of stored protons
    UShort_t GetNAPro()const{return fNAProTracks;}     // .. anti-protons
    UShort_t GetNLam()const{return fNLamTracks;}       // .. lambda 
    UShort_t GetNALam()const{return fNALamTracks;}     // .. anti-lambda
    UShort_t GetNBgLam()const{return fNBgLamTracks;}   // .. background lambda
    UShort_t GetNBgALam()const{return fNBgALamTracks;} // .. background anti-lambda

    // Data member, sorry for this private public private mixture,
    // but what's this reorder warning? sometimes it's there, st not..
  private:
    // Size of the arrays for tracks. Constant and private to reflect the idea that
    // the arrays get allocated only once to avoid excessive memory allocation
    const UShort_t fPriTrackLim;   // Limit for primary tracks
    const UShort_t fV0Lim;         // Limit for V0s
  public:
    // Pointer to array of ...
    FemtoBufferTrack *fProTracks;          //! Proton tracks
    FemtoBufferTrack *fAProTracks;         //! Anti-proton tracks
    
    FemtoBufferV0    *fLamTracks;          //! Lambda tracks
    FemtoBufferV0    *fALamTracks;         //! Anti-lambda tracks
    
    FemtoBufferV0    *fBgLamTracks;        //! Background lambda tracks
    FemtoBufferV0    *fBgALamTracks;       //! Background anti-lambda tracks
  private:
    // Number of stored tracks in the event
    UShort_t fNProTracks;    // Number of stored protons
    UShort_t fNAProTracks;   // Number of stored anti-protons
    UShort_t fNLamTracks;    // Number of stored lambdas
    UShort_t fNALamTracks;   // Number of stored anti-lambdas
    UShort_t fNBgLamTracks;  // Number of stored lambdas
    UShort_t fNBgALamTracks; // Number of stored anti-lambdas

    // Double_t needed??? magnetic field probably is like 5.0 and not 5.0000120047
    Double_t fBfield;               // Magnetic field in ALICE unit [kG]
    Double_t fPriVtxPos[3];         // Primary vtx position
  };
  
  class FemtoBuffer { // Holds the events
  public:
    FemtoBuffer();           // Dummy constructor
    FemtoBuffer(const UChar_t ZvertexBins,const UChar_t CentBins,const UChar_t MixBuff,
		const UShort_t PriTrackLim,const UShort_t V0Lim,const Float_t AbsZvertexCut,const Float_t CentCut); // Constructor
    FemtoBuffer(const FemtoBuffer &fb); //Ctor
    FemtoBuffer& operator=(const AliAnalysisTaskProtonLambda::FemtoBuffer&); // Assignment
    ~FemtoBuffer();          // Destructor
    void ShiftAndAdd(AliAODEvent *evt); // Discard last event, shift all, set first one
    void ShiftAndAdd(const Double_t bfield,const Double_t priVtxPos[3],const Float_t centrality); // Discard last event, shift all, set first one
    FemtoBufferEvent *GetEvt(const UChar_t i)const{return fCurEvt[i];}; // Returns a pointer to the i'th event of the current event mixing class
    UChar_t GetMixBuffSize()const{return fkMixBuffSize;}// Returns the number of events held in every mixing bin
                                 
  private:
    const UChar_t fkZvertexBins;           // Number of bins in Zvertex
    const UChar_t fkCentBins;              // Number of bins in centrality
    const UChar_t fkMixBuffSize;           // Number of stored events
    const UShort_t fkPriTrackLim;          // Buffer size protons per event
    const UShort_t fkV0Lim;                // Buffer size lambdas per event

    const TAxis *fZvertexAxis;          //! To find Zvertex bin
    const TAxis *fCentAxis;             //! To find centrality bin

    FemtoBufferEvent **fCurEvt;//! Array of pointer to the set of the current events
                                     //  Note that the pointer won't be constant
    FemtoBufferEvent ****fEC;  //! The internal thing where the events
                                     //  are stored. 
                                     //  fEC stands for event collection.
  };
  
  ClassDef(AliAnalysisTaskProtonLambda, 1);
};

#endif
