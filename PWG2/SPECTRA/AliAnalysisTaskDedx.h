#ifndef AliAnalysisTaskDedx_cxx
#define AliAnalysisTaskDedx_cxx

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskDedx class
//            This task is for QAing the dE/dx from the ESD
//              Origin: B.H. Nov2007, hippolyt@in2p3.fr
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskDedx : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDedx();
  AliAnalysisTaskDedx(const char    *rName,
		      const Bool_t   rAllConstrainedFlag,
		      const Bool_t   rMidPseudoRapidityFlag,
		      const Bool_t   rSelTrackRemoveKink,
		      const Bool_t   rSelTrackWithOnTheFlyV0,
		      const Int_t    rSelTrackMinClustersTPC,
		      const Int_t    rSelTrackMinClustersITS,
		      const Float_t  rSelTrackMaxChi2PerClusterTPC,
		      const Float_t  rSelTrackMaxChi2PerClusterITS,
		      const Double_t rSelTrackMaxCov11,
		      const Double_t rSelTrackMaxCov22,
		      const Double_t rSelTrackMaxCov33,
		      const Double_t rSelTrackMaxCov44,
		      const Double_t rSelTrackMaxCov55,
		      const Double_t rSelV0MaxDcaDaughters,
		      const Double_t rSelV0MinDecayLength);

  virtual ~AliAnalysisTaskDedx() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void   SetCollidingSystems(Int_t collidingSystems = 0) {fCollidingSystems = collidingSystems;}
  void   SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}

  void   SetAllConstrainedFlag(Bool_t allConstrainedFlag = 0)                      {fAllConstrainedFlag = allConstrainedFlag;}
  void   SetMidPseudoRapidityFlag(Bool_t midPseudoRapidityFlag = 0)                {fMidPseudoRapidityFlag = midPseudoRapidityFlag;}

  void   SetSelTrackRemoveKink(Bool_t selTrackRemoveKink = 0)                      {fSelTrackRemoveKink = selTrackRemoveKink;}
  void   SetSelTrackWithOnTheFlyV0(Bool_t selTrackWithOnTheFlyV0 = 0)              {fSelTrackWithOnTheFlyV0 = selTrackWithOnTheFlyV0;}
  void   SetSelTrackMinClustersTPC(Int_t selTrackMinClustersTPC = 0)               {fSelTrackMinClustersTPC = selTrackMinClustersTPC;}
  void   SetSelTrackMinClustersITS(Int_t selTrackMinClustersITS = 0)               {fSelTrackMinClustersITS = selTrackMinClustersITS;}
  void   SetSelTrackMaxChi2PerClusterTPC(Float_t selTrackMaxChi2PerClusterTPC = 0) {fSelTrackMaxChi2PerClusterTPC = selTrackMaxChi2PerClusterTPC;}
  void   SetSelTrackMaxChi2PerClusterITS(Float_t selTrackMaxChi2PerClusterITS = 0) {fSelTrackMaxChi2PerClusterITS = selTrackMaxChi2PerClusterITS;}
  void   SetSelTrackMaxCov11(Double_t selTrackMaxCov11 = 0)                        {fSelTrackMaxCov11 = selTrackMaxCov11;}
  void   SetSelTrackMaxCov22(Double_t selTrackMaxCov22 = 0)                        {fSelTrackMaxCov22 = selTrackMaxCov22;}
  void   SetSelTrackMaxCov33(Double_t selTrackMaxCov33 = 0)                        {fSelTrackMaxCov33 = selTrackMaxCov33;}
  void   SetSelTrackMaxCov44(Double_t selTrackMaxCov44 = 0)                        {fSelTrackMaxCov44 = selTrackMaxCov44;}
  void   SetSelTrackMaxCov55(Double_t selTrackMaxCov55 = 0)                        {fSelTrackMaxCov55 = selTrackMaxCov55;}
  void   SetSelV0MaxDcaDaughters(Double_t selV0MaxDcaDaughters = 0)                {fSelV0MaxDcaDaughters = selV0MaxDcaDaughters;}
  void   SetSelV0MinDecayLength(Double_t selV0MinDecayLength = 0)                  {fSelV0MinDecayLength = selV0MinDecayLength;}
  
 private:

  Bool_t         IsAccepted(AliESDtrack *track);
  //  Float_t        GetSigmaToVertex(AliESDtrack* esdTrack);
  TString      fAnalysisType;                       //  ESD or AOD
  Int_t        fCollidingSystems;                   //  Colliding systems 0/1 for pp/PbPb
  TList       *fListHist;                           //! List of histograms
  TH1F        *fHistPtot;                           //! Ptot spectrum
  TH1F        *fHistMultiplicity;                   //! Multiplicity distribution
  TH2F        *fHistTPCDedxVsMomentum;              //! TPC dE/dx vs momemtum
  TH2F        *fHistITSDedxVsMomentum;              //! ITS dE/dx vs momemtum
  TH1F        *fHistMassK0;                         //! Invariant Mass of K0s
  TH1F        *fHistMassLambda;                     //! Invariant Mass of Lambda
  TH1F        *fHistMassAntiLambda;                 //! Invariant Mass of Anti-Lambda
  TH2F        *fHistTPCDedxVsMomPosK0;              //! TPC dE/dx vs momemtum for K0 positive daughter
  TH2F        *fHistTPCDedxVsMomNegK0;              //! TPC dE/dx vs momemtum for K0 negative daughter
  TH2F        *fHistTPCDedxVsMomPosLambda;          //! TPC dE/dx vs momemtum for Lambda positive daughter
  TH2F        *fHistTPCDedxVsMomNegLambda;          //! TPC dE/dx vs momemtum for Lambda negative daughter
  TH2F        *fHistTPCDedxVsMomPosAntiLambda;      //! TPC dE/dx vs momemtum for Anti-Lambda positive daughter
  TH2F        *fHistTPCDedxVsMomNegAntiLambda;      //! TPC dE/dx vs momemtum for Anti-Lambda negative daughter
  TH1F        *fHistDiffInOutMomentum;              //! Difference between inner and outer TPC momenta
  TH1F        *fHistDiffPrimOutMomentum;            //! Difference between primary and outer TPC momenta
  TH1F        *fHistDiffPrimMeanMomentum;           //! Difference between primary and (inner+outer)/2 TPC momenta
  TH1F        *fHistPercPrimMeanMomentum;           //! Percentage between primary and (inner+outer)/2 TPC momenta
  TH1F        *fHistPrimEta;                        //! Pseudorapidity distribution
  TH2F        *fHistPercPrimMeanMomentumVsEta;      //! Same as fHistDiffPrimMeanMomentum but vs pseudorapidity
  TH2F        *fHistPercPrimMeanMomentumVsPrim;     //! Same as fHistDiffPrimMeanMomentum but vs primary momentum
	  
  TH1F        *fHistMultiplicityCuts;               //! Same as above but once primary track cuts applied
  TH2F        *fHistTPCDedxVsMomentumCuts;          //!
  TH2F        *fHistITSDedxVsMomentumCuts;          //!
  TH1F        *fHistMassK0Cuts;                     //!
  TH1F        *fHistMassLambdaCuts;                 //!
  TH1F        *fHistMassAntiLambdaCuts;             //!
  TH2F        *fHistTPCDedxVsMomPosK0Cuts;          //!
  TH2F        *fHistTPCDedxVsMomNegK0Cuts;          //!
  TH2F        *fHistTPCDedxVsMomPosLambdaCuts;      //!
  TH2F        *fHistTPCDedxVsMomNegLambdaCuts;      //!
  TH2F        *fHistTPCDedxVsMomPosAntiLambdaCuts;  //!
  TH2F        *fHistTPCDedxVsMomNegAntiLambdaCuts;  //!
  TH1F        *fHistDiffInOutMomentumCuts;          //!
  TH1F        *fHistDiffPrimOutMomentumCuts;        //!
  TH1F        *fHistDiffPrimMeanMomentumCuts;       //!
  TH1F        *fHistPercPrimMeanMomentumCuts;       //!
  TH1F        *fHistPrimEtaCuts;                    //!
  TH2F        *fHistPercPrimMeanMomentumVsEtaCuts;  //!
  TH2F        *fHistPercPrimMeanMomentumVsPrimCuts; //!

  Bool_t       fAllConstrainedFlag;                 // Primary vertex constrain requirement
  Bool_t       fMidPseudoRapidityFlag;              // Mid-eta requirement

                                                    // Track selections: streaming allowed and needed !
  Bool_t       fSelTrackRemoveKink;                 // Remove kink candidates
  Bool_t       fSelTrackWithOnTheFlyV0;             // Select daughter tracks from on-the-fly V0s
  Int_t        fSelTrackMinClustersTPC;             // Minimum number of cluster in the TPC
  Int_t        fSelTrackMinClustersITS;             // Minimum number of cluster in the ITS
  Float_t      fSelTrackMaxChi2PerClusterTPC;       // Maximum chisq per cluster in the TPC
  Float_t      fSelTrackMaxChi2PerClusterITS;       // Maximum chisq per cluster in the ITS
  Double_t     fSelTrackMaxCov11;                   // Maximum value for cov.mat. diag. element
  Double_t     fSelTrackMaxCov22;                   // Maximum value for cov.mat. diag. element
  Double_t     fSelTrackMaxCov33;                   // Maximum value for cov.mat. diag. element
  Double_t     fSelTrackMaxCov44;                   // Maximum value for cov.mat. diag. element
  Double_t     fSelTrackMaxCov55;                   // Maximum value for cov.mat. diag. element
  Double_t     fSelV0MaxDcaDaughters;               // Maximum value for DCA between V0 daughter tracks
  Double_t     fSelV0MinDecayLength;                // Minimum value for V0 decay length

  AliAnalysisTaskDedx(const AliAnalysisTaskDedx&);            // not implemented
  AliAnalysisTaskDedx& operator=(const AliAnalysisTaskDedx&); // not implemented
  
  ClassDef(AliAnalysisTaskDedx, 1);
};

#endif
