#ifndef AliAnalysisTaskESDDedx_cxx
#define AliAnalysisTaskESDDedx_cxx

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskESDDedx class
//            This task is for QAing the dE/dx from the ESD
//              Origin: B.H. Nov2007, hippolyt@in2p3.fr
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class AliESDEvent;

#include "AliAnalysisTask.h"

class AliAnalysisTaskESDDedx : public AliAnalysisTask {
 public:
  AliAnalysisTaskESDDedx(const char    *rName = "AliAnalysisTaskESDDedx",
			 const Bool_t   rAllConstrainedFlag     = kFALSE,
			 const Bool_t   rMidPseudoRapidityFlag  = kFALSE,
			 const Bool_t   rSelTrackRemoveKink     = kTRUE,
			 const Bool_t   rSelTrackWithOnTheFlyV0 = kFALSE,
			 const Int_t    rSelTrackMinClustersTPC = 50,
			 const Int_t    rSelTrackMinClustersITS =  0,
			 const Float_t  rSelTrackMaxChi2PerClusterTPC = 3.5,
			 const Float_t  rSelTrackMaxChi2PerClusterITS = 10,
			 const Double_t rSelTrackMaxCov11 = 2.0,
			 const Double_t rSelTrackMaxCov22 = 2.0,
			 const Double_t rSelTrackMaxCov33 = 0.5,
			 const Double_t rSelTrackMaxCov44 = 0.5,
			 const Double_t rSelTrackMaxCov55 = 2.0,
			 const Double_t rSelV0MaxDcaDaughters = 0.5,
			 const Double_t rSelV0MinDecayLength  = 0.0);

  virtual ~AliAnalysisTaskESDDedx() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:

  Bool_t         IsAccepted(AliESDtrack *track);
  //  Float_t        GetSigmaToVertex(AliESDtrack* esdTrack);

  AliESDEvent *fESD;                                //! ESD object
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

  AliAnalysisTaskESDDedx(const AliAnalysisTaskESDDedx&);            // not implemented
  AliAnalysisTaskESDDedx& operator=(const AliAnalysisTaskESDDedx&); // not implemented
  
  ClassDef(AliAnalysisTaskESDDedx, 1);
};

#endif
