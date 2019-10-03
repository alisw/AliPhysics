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
class TH1F;
#include <TH2F.h>
class TH3F;
#include <THn.h>
class TList;
class TAxis;

#include <AliAnalysisTaskSE.h>
class AliTPCPIDResponse;
class AliPIDResponse;

#include <AliAODEvent.h>
#include <AliAODVertex.h>
#include <AliAODv0.h>
#include <AliAODTrack.h>

class AliAnalysisTaskProtonLambda2d : public AliAnalysisTaskSE {
  // enum for cut variation, for now DCA
 public:
  enum {kStd,kWide,kStrict}; //0,1,2
  
  // Forward declaration for nested classes
 private:
  class FemtoBufferTrack;
  class FemtoBufferV0;
  class FemtoBufferEvent;
  class FemtoBuffer;
 public:
  AliAnalysisTaskProtonLambda2d();                 // Dummy constructor
  AliAnalysisTaskProtonLambda2d(const char *name
				,const Float_t centCutLo=0.
				,const Float_t centCutHi=10.
				,const Float_t minvCut=0.004
				,const Bool_t doBgLam=kFALSE
				,const Bool_t useOnTheFly=kTRUE
				,const TString centEst="V0M"
				,const Int_t dcaFlag=kStd); // Constructor
  virtual ~AliAnalysisTaskProtonLambda2d();                // Destructor
  
  void   UserCreateOutputObjects();  // Called once at the beginning
  void   UserExec(Option_t *option); // Main loop
  void   Terminate(Option_t *);      // Called once at the end

  void SetLamPurHist(const TH2F &h){fLamPur = h;fLamPur.SetName("hLamPurT");}   // Setters for the corrections
  void SetALamPurHist(const TH2F &h){fALamPur = h;fALamPur.SetName("hALamPurT");}//             "

  void SetLamFdPur010Hist(const TH2F &h){fLamFdPur010=h;fLamFdPur010.SetName("hLamFdPur010T");}    // lam fd pur
  void SetLamFdPur1020Hist(const TH2F &h){fLamFdPur1020=h;fLamFdPur1020.SetName("hLamFdPur1020T");}
  void SetLamFdPur2040Hist(const TH2F &h){fLamFdPur2040=h;fLamFdPur2040.SetName("hLamFdPur2040T");}
  void SetLamFdPur4060Hist(const TH2F &h){fLamFdPur4060=h;fLamFdPur4060.SetName("hLamFdPur4060T");} 
  void SetALamFdPur010Hist(const TH2F &h){fALamFdPur010=h;fALamFdPur010.SetName("hALamFdPur010T");}  // alam fd pur
  void SetALamFdPur1020Hist(const TH2F &h){fALamFdPur1020=h;fALamFdPur1020.SetName("hALamFdPur1020T");}
  void SetALamFdPur2040Hist(const TH2F &h){fALamFdPur2040=h;fALamFdPur2040.SetName("hALamFdPur2040T");}
  void SetALamFdPur4060Hist(const TH2F &h){fALamFdPur4060=h;fALamFdPur4060.SetName("hALamFdPur4060T");}
  
  void SetProFdPurHist(const TH2F &h){fProFdPur=h;fProFdPur.SetName("hProFdPurT");}    // pro fd pur
  void SetAProFdPurHist(const TH2F &h){fAProFdPur=h;fAProFdPur.SetName("hAProFdPurT");}

  /* Bool_t SetCentRange(const Float_t CentCutLo,const Float_t CentCutHi){ */
  /*   fkCentCutLo=CentCutLo;fkCentCutHi=CentCutHi;return kTRUE;} // Setter for the centrality cut */

  // Setters to skip some unneccessary stuff
  //  void   SetDoLamOnly(Bool_t doLamOnly=kTRUE){fdoLamOnly=doLamOnly;} 
  //void   SetDoDCAHists(Bool_t b=kTRUE){fDoDCAHists=b;}
  //  void   SetDoMtHists(Bool_t b=kTRUE){fDoMtHists=b;}

 private:
  AliAnalysisTaskProtonLambda2d(const AliAnalysisTaskProtonLambda2d& atpl); // Not implemented
  AliAnalysisTaskProtonLambda2d& operator=(const AliAnalysisTaskProtonLambda2d& atpl); //Not implemented

  void   ProcessOffline(const AliAODv0 *v0,const AliAODTrack *pTrack,const AliAODTrack *nTrack);    //! Use offline V0 finder
  void   ProcessOnTheFly(const AliAODv0 *v0,const AliAODTrack *pTrack,const AliAODTrack *nTrack);   //! Use on-the-fly V0 finder
  void   ProcessTOF(const AliAODTrack *track);         // Use Tof for 1.0 GeV/c < p < 3.25 GeV/c
  void   ProcessTPC(const AliAODTrack *track);         // Use TPC for p < 0.75 GeV/c
  void   ProcessHybrid(const AliAODTrack *track);      // Preselection via TPC, identification w/ Tof 0.75 GeV/c < p < 1.0 GeV/c
  void   CleaningProcedure();         // Assure uniqueness of all daughters and primaries
  void   ProcessReal();               // TTR cut & qinv calc for real events
  void   ProcessMixed();              // TTR cut & qinv calc for mixed events
  void   ProcessRealBackground();     // Check procedure with lambda background from sideband, real events
  void   ProcessMixedBackground();    // Check procedure with lambda background from sideband, real events

  // Pair observables
  Float_t dEtaS(const FemtoBufferTrack &track1,
	       const FemtoBufferTrack &track2); // delta eta star(R=1.25m)
  Float_t dThetaS(const FemtoBufferTrack &track1,
		  const FemtoBufferTrack &track2); // delta theta star(R=1.25m)
  Float_t dPhiSAtR12(const FemtoBufferTrack &track1,
		     const FemtoBufferTrack &track2); // delta phi star (R=1.25m)
  Float_t Qinv(const FemtoBufferV0 &v01,
	       const FemtoBufferV0 &v02);           // Calcs Qinv for (anti-)Lam(anti-)Lam
  Float_t Qinv(const FemtoBufferV0 &v0,
	       const FemtoBufferTrack &track);       // Calcs Qinv for (anti-)Lam(anti-)Pro
  Float_t Qinv(const FemtoBufferTrack &track,
	       const FemtoBufferV0 &v0);       // Calls Qinv(v0, track)
  Float_t QinvProPro(const FemtoBufferTrack &proTrack1,
		     const FemtoBufferTrack &proTrack2);               // Calcs Qinv for proton proton
  Float_t QinvPioPro(const FemtoBufferTrack &pioTrack,
		     const FemtoBufferTrack &proTrack);                // Calcs Qinv for pion proton
  Float_t QinvConstr(const FemtoBufferV0 &v0,
		     const FemtoBufferTrack &track); // w/ vtx constraint for pri pro
  Float_t QinvConstr(const FemtoBufferTrack &track,
		     const FemtoBufferV0 &v0); // Calls QinvConstr(v0, track)
  Float_t Minv(const FemtoBufferV0 &v01,
	       const FemtoBufferV0 &v02);           // Calcs Minv for (anti-)Lam(anti-)Lam
  Float_t Minv(const FemtoBufferV0 &v0,
	       const FemtoBufferTrack &track);       // Calcs Minv for (anti-)Lam(anti-)Pro
  Float_t Minv(const FemtoBufferTrack &track,
	       const FemtoBufferV0 &v0);       // Calls Minv(v0, track)
  Float_t mt(const FemtoBufferTrack &track,
	     const FemtoBufferV0 &v0);         // mt of pair
  Float_t mt(const FemtoBufferV0 &v0,
	     const FemtoBufferTrack &track);         //    ""
  Float_t mt(const FemtoBufferV0 &v01,
	     const FemtoBufferV0 &v02);             //    ""
  static Float_t ktSquared(const FemtoBufferV0 &v01,
			   const FemtoBufferV0 &v02);      // kt squared
  static Float_t ktSquared(const FemtoBufferTrack &track,
			   const FemtoBufferV0 &v0);  //    ""
  static Float_t ktSquared(const FemtoBufferV0 &v0,
			   const FemtoBufferTrack &track);  //    ""
  static Float_t yPair(const FemtoBufferTrack &track,
		       const FemtoBufferV0 &v0);      // pair rapidity
  static Float_t yPair(const FemtoBufferV0 &v0,
		       const FemtoBufferTrack &track);      // pair rapidity

  Bool_t  goodDCA(const AliAODTrack *track);          // Does a cut on the DCAxy value and fills a histogram
  Float_t RapidityProton(const AliAODTrack *track);   // Rapidity assuming proton mass
  // For AOD tracks, need a function to get the dca
  static Bool_t DCAxyz(const AliAODTrack *track, const AliVEvent *evt,Double_t dca[2]); 
  static Bool_t getRapAndPtBin(const Float_t rap,const Float_t pt,
			UChar_t &rapBin,UChar_t &ptBin);
  void StoreGlobalTrackReference(const AliAODTrack *track);      // Store pointer to global track
  void ResetGlobalTrackReference();                        // Reset all pointers to global tracks
  Bool_t acceptTrack(const AliAODTrack *track);            // Famous crossed rows / findable clusters
  Bool_t GoodTPCFitMapSharedMap(const AliAODTrack *track); // Rejects shared clusters, primaries
  Bool_t GoodTPCFitMapSharedMap(const AliAODTrack *pTrack,
				const AliAODTrack *nTrack);// Rejects shared clusters, two V0 daughters
  Float_t GetCorrectedTOFSignal(const AliAODTrack *track);   // Return correct TOF signal for old and new AODs
  static Float_t GetTOFSignal(const AliAODTrack *track);     // Wrapper for AliAODTrack::GetTOFSignal for crash
  static Float_t GetProtonTime(const AliAODTrack *track);    // Wrapper for AliAODTrack::GetIntegratedTimes 

  void FillDedxHist(const AliAODTrack *track);               // Fill dE/dx histograms

  Float_t GetLamPur(Float_t y,Float_t pt); // Returns the purity for given y,pt 
  Float_t GetALamPur(Float_t y,Float_t pt); // Returns the purity for given y,pt 
  Float_t GetLamFdPur(Float_t y,Float_t pt,const Float_t cent); // Returns the purity for given y,pt,cent
  Float_t GetALamFdPur(Float_t y,Float_t pt,const Float_t cent); // Returns the purity for given y,pt,cent 
  Float_t GetProFdPur(Float_t y,Float_t pt); // Returns the purity for given y,pt 
  Float_t GetAProFdPur(Float_t y,Float_t pt); // Returns the purity for given y,pt
  Float_t GetPur(Float_t y, Float_t pt,const TH2 *h); // Returns bin cont. of h at y,pt
 
  //
  // ----  Data member  ----
  //
  const Bool_t   fkDoStdHists;             // Do some std histograms
  const Bool_t   fkDoDCAHists;             // Create the 2d DCA histograms?
  const Bool_t   fkDoBgLamALam;            // Proton side-band lambda correlation?
  const Bool_t   fkUseOnTheFly;            // Use on-the-fly or offline V0 finder
  const Float_t  fkAbsZvertexCut;          // Values of Zvertex and centrality cut
  const Float_t  fkCentCutLo;              // are needed for event mixing
  const Float_t  fkCentCutHi;              //
  const Float_t  fkMinvCut;                // Minv cut on lambda and anti-lambda
  const TString  fkCentEst;                // Centrality estimator used
  //  Bool_t   fdoLamOnly;                     // Get's raw lambdas only, disabling the rest
  const Float_t  fdEtaSCut;                // TTR cut in dEtaS
  const Float_t  fdPhiSCut;                // TTR cut in dPhiS
  const Float_t  fkLamMass;                // PDG-Mass of the lambda
  const Float_t  fkProMass;                // PDG-Mass of the proton
  const Float_t  fkPioMass;                // PDG-Mass of the pion
  const Int_t    fkDCAFlag;                // enum std wide strict
  
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
  const AliAODTrack **fGTI;                //! Array of pointers, just nicely sorted according to the id
  const UShort_t  fTrackBuffSize;          //! Size of the above array, ~12000 for PbPb

  // Correction histograms, to be set by the AddTaskMacro
  TH2F fLamPur,fALamPur;                   // Correction for the (anti-)lambda purity in y,pt
  TH2F fLamFdPur010,fLamFdPur1020,fLamFdPur2040,fLamFdPur4060,
    fALamFdPur010,fALamFdPur1020,fALamFdPur2040,fALamFdPur4060; // Feed-down correction for lambda
  TH2F fProFdPur,fAProFdPur;               // Material and weak decays for (anti-)protons
  
  //
  //  ---------        Histograms        ---------
  //
  TH1F        *fHistGoodEvent;                  //! Control hist for event cuts
  TH2F        *fHistCentComp;                   //! V0M (std) vs ZEMvsZDC centrality

  TH1F        *fHistSideBandOffLam;             //! Off: side band background lambda
  TH1F        *fHistSideBandOffALam;            //! Off: side band background anti-lambda
  TH1F        *fHistMassLambdaOff;              //! Off: invariant mass assuming lambda
  TH1F        *fHistMassAntiLambdaOff;          //! Off: invariant mass assuming anti-lambda
  TH3F        *fHistYPtMassLamOff;              //! Off: lam-mass vs y-pt
  TH3F        *fHistYPtMassALamOff;             //! Off: alam-mass vs y-pt
  
  TH1F        *fHistSideBandOnLam;              //! On-the-fly: side band background lambda
  TH1F        *fHistSideBandOnALam;             //! On-the-fly: side band background anti-lambda
  TH1F        *fHistMassLambdaOn;               //! On: invariant mass assuming lambda
  TH1F        *fHistMassAntiLambdaOn;           //! On: invariant mass assuming anti-lambda
  TH3F        *fHistYPtMassLamOn;               //! On: lam-mass vs y-pt
  TH3F        *fHistYPtMassALamOn;              //! On: alam-mass vs y-pt

  // Primary particles
  TH1F        *fPriHistShare;                   //! Primaries: number of shared clusters
  // ToF histograms
  TH2F        *fPriHistTOFsignalPosVsP;         //! Pri: TOF signal vs p (pos)
  TH2F        *fPriHistTOFsignalNegVsP;         //! Pri: TOF signal vs p (neg)
  // Hybrid analysis: dEdx & ToF histograms
  TH1F        *fPriHistHybridTOFsigPosWoTPC;    //! Pri: TOF signal without dEdx selection (pos)
  TH1F        *fPriHistHybridTOFsigPosTPCok;    //! Pri: TOF signal with dEdx selection (pos)
  TH1F        *fPriHistHybridTOFsigNegWoTPC;    //! Pri: TOF signal without dEdx selection (neg)
  TH1F        *fPriHistHybridTOFsigNegTPCok;    //! Pri: TOF signal with dEdx selection (neg)
  // dEdx
  TH2F        *fPriHistTPCsignalPos;            //! Pri: TPC dE/dx signal vs p
  TH2F        *fPriHistTPCsignalLowPPos;        //! Pri: dEdx for 0.1 < p < 0.3
  TH2F        *fPriHistTPCsignalMedPPos;        //! Pri: dEdx for 0.3 < p < 0.9
  TH2F        *fPriHistTPCsignalHigPPos;        //! Pri: dEdx for 0.9 < p < 1.9
  TH2F        *fPriHistTPCsignalNeg;            //! Pri: TPC dE/dx signal vs p
  TH2F        *fPriHistTPCsignalLowPNeg;        //! Pri: dEdx for 0.1 < p < 0.3
  TH2F        *fPriHistTPCsignalMedPNeg;        //! Pri: dEdx for 0.3 < p < 0.9
  TH2F        *fPriHistTPCsignalHigPNeg;        //! Pri: dEdx for 0.9 < p < 1.9
  // All
  TH2F        ***fPriHistDCAxyzYPtPro;          //! Pri: DCAxy vs DCAz (y,pt) for protons
  TH2F        ***fPriHistDCAxyzYPtAPro;         //! Pri: DCAxy vs DCAz (y,pt) for anti-protons

  // Monitor needed buffer size
  TH1F        *f2HistNLamBefClean;              //! Number of lambdas in buffer before cleaning procedure
  TH1F        *f2HistNProBefClean;              //! Number of protons in buffer before cleaning procedure
  TH1F        *f2HistNALamBefClean;             //! Number of anti-lambdas in buffer before cleaning procedure
  TH1F        *f2HistNAProBefClean;             //! Number of anti-protons in buffer before cleaning procedure
  
  // Histograms for two particle observables..
  // ..2d distance at R=1.2m (with shifted vertices)
  TH2F        *f2HistLamProAngDistSft2dAtR12Real;     //!
  TH2F        *f2HistALamAProAngDistSft2dAtR12Real;   //!
  TH2F        *f2HistBgLamProAngDistSft2dAtR12Real;   //!
  TH2F        *f2HistBgALamAProAngDistSft2dAtR12Real; //!

  TH2F        *f2HistLamProAngDistSft2dAtR12Mixed;    //!
  TH2F        *f2HistALamAProAngDistSft2dAtR12Mixed;  //!
  TH2F        *f2HistBgLamProAngDistSft2dAtR12Mixed;  //!
  TH2F        *f2HistBgALamAProAngDistSft2dAtR12Mixed;//!
    
  // transverse mass mt of the pair
  TH1F        *f2HistMtLamProReal;     //!
  TH1F        *f2HistMtALamAProReal;   //!
  // .. for low q pairs only
  TH1F        *f2HistMtLowQLamProReal;     //!
  TH1F        *f2HistMtLowQALamAProReal;   //!
  // mT vs phi*
  TH2F        *f2HistMtVsPhiSLowQ; //!

  // Corr fct vs distances
  THnF       *fLamProReal;  //!
  THnF       *fALamAProReal;//!
  // Pairs vs qinv, w and w/o purity correction
  TH1F *f2HistLamProWLamPurReal,*f2HistLamProWLamFdPurReal,
    *f2HistLamProWProFdPurReal,*f2HistLamProWoPurReal,
    *f2HistALamAProWALamPurReal,*f2HistALamAProWALamFdPurReal,
    *f2HistALamAProWAProFdPurReal,*f2HistALamAProWoPurReal;
  // Pairs vs qinv vs mT, w and w/o purity correction
  TH2F *f2HistLamProWLamPurmTReal,*f2HistLamProWLamFdPurmTReal,
    *f2HistLamProWProFdPurmTReal,*f2HistLamProWoPurmTReal,
    *f2HistALamAProWALamPurmTReal,*f2HistALamAProWALamFdPurmTReal,
    *f2HistALamAProWAProFdPurmTReal,*f2HistALamAProWoPurmTReal;
  
  // Purity fluctuations for qinv .0 - .2
  TH1F *f2HistPurFlucLamPro,*f2HistPurFlucALamAPro;
  // Phase space of lambdas, protons and pairs for qinv < .2,
  // also for lambdas of pairs with .4 < qinv < .5
  TH2F *f2HistLamYPtPair, *f2HistALamYPtPair,
    *f2HistProYPtPair, *f2HistAProYPtPair,
    *f2HistYKtPair, *f2HistYKtAPair,
    *f2HistLamYPtHiPair, *f2HistALamYPtHiPair;

  // Background lambdas: real events
  THnF       *fBgLamProReal;  //!
  THnF       *fBgALamAProReal;//!
  // Mixed events
  THnF       *fLamProMixed; //!
  THnF       *fALamAProMixed;//!
   // Background lambdas: mixed events
  THnF       *fBgLamProMixed;   //!
  THnF       *fBgALamAProMixed; //!

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
    /* void GetShiftedPositionAtShiftedRadii(const AliAODTrack *track, */
    /* 					  const Float_t bfield, const Float_t priVtx[3]);  */
    /* // Global if not shifted */
    /* void GetGlobalPositionAtGlobalRadii(const AliAODTrack *track, */
    /* 					const Float_t bfield); */
    // Position at R=1.25m in shifted coordinate system
    void SetSftPosR125(const AliAODTrack *track
		       ,const Float_t bfield
		       ,const Float_t priVtx[3]); 

    // Retrieve kinematics
    Double_t RapPro() const;        // Calcs rapidity assuming proton mass
    Double_t Pt() const;            // Calcs pt
    // Cleaning procedure flag
    Bool_t UseIt() const {return fID!=65535;}          // Use entry? eg. set to false by cleaning procedure
    void   SetBadFlag(){fID=65535;}                   // Can only set 'bad track' flag
    // Check for good propagation
    /* Bool_t GoodPropR12(){return fXshifted[2][0]>-9998.;} // Checks for good propagation to R=1.2(5)m */
    Bool_t GoodPropR12()const{return fXSftR125[0]>-9998.;} // Checks for good propagation to R=1.2(5)m
    // Getter for dervied TTR variables
    Double_t ThetaS()const; // Returns the longitudinal angle of the particles propagated position at R=1.25m
    Double_t EtaS()const;   // Returns the corresponding eta of a pri. part. with pos at R=1.25m

    // Interpretation of ALICE coding rule RC14 is that FemtoBufferTrack is a private member
    // of AliAnalysisTaskProtonLambda2d, so the following data member are private
    UShort_t fID;          // Unique track id (->AliAODTrack.h), UShort_t goes to 65000
    Double_t fP[3];        // Momentum of track
    /* Float_t  fXglobal[9][3];    //! Global positions at different global radii */
    /* Float_t  fXshifted[9][3];   //! Shifted positions at different shifted radii */
    Float_t fXSftR125[3];  // Spatial position at R=1.25m in shifted coordinate system
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
    Double_t RapLam() const;        // Calcs rapidity assuming lambda mass
    Double_t Pt() const;            // Calcs pt
    // Data member
    Double_t fP[3];                 // Momentum (x,y,z)
    Float_t fCosPoint;              // Cosine of pointing angle 
    FemtoBufferTrack fPosDaughter;  // Positive daughter of the V0
    FemtoBufferTrack fNegDaughter;  // Negative daughter fo the V0
  };
  
  class FemtoBufferEvent{// Event 
  public:
    FemtoBufferEvent();                       // Constructor
    FemtoBufferEvent(const UShort_t priTrackBuff,const UShort_t V0Buff, 
		     const Bool_t DoBgLamALam,
		     const Double_t bfield,const Double_t priVtxPos[3]); // Constructor
    FemtoBufferEvent(const UShort_t priTrackBuff,const UShort_t V0Buff,
		     const Bool_t DoBgLamALam); // Constructor
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
    Bool_t   GetDoBgLamALam() const{return fDoBgLamALam;} // Store bg lam?
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
    const Bool_t   fDoBgLamALam;   // Store bg lam?
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
		const UShort_t PriTrackLim,const UShort_t V0Lim,const Float_t AbsZvertexCut,
		const Float_t CentCutLo,const Float_t CentCutHi,const Bool_t DoBgLamALam); // Constructor
    FemtoBuffer(const FemtoBuffer &fb); //Ctor
    FemtoBuffer& operator=(const AliAnalysisTaskProtonLambda2d::FemtoBuffer&); // Assignment
    ~FemtoBuffer();          // Destructor
    void ShiftAndAdd(AliAODEvent *evt,const TString centEst); // Discard last event, shift all, set first one
    void ShiftAndAdd(const Double_t bfield,const Double_t priVtxPos[3],const Float_t centrality); // Discard last event, shift all, set first one
    FemtoBufferEvent *GetEvt(const UChar_t i)const{return fCurEvt[i];}; // Returns a pointer to the i'th event of the current event mixing class
    UChar_t GetMixBuffSize()const{return fkMixBuffSize;}// Returns the number of events held in every mixing bin
                                 
  private:
    const UChar_t fkZvertexBins;           // Number of bins in Zvertex
    const UChar_t fkCentBins;              // Number of bins in centrality
    const UChar_t fkMixBuffSize;           // Number of stored events
    const Bool_t  fDoBgLamALam;            // Enable background lam/alam
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
  
  ClassDef(AliAnalysisTaskProtonLambda2d, 1);
};

#endif
