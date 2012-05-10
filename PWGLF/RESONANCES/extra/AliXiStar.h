#ifndef ALIXISTAR_H
#define ALIXISTAR_H
//
// Class AliXiStar
//
// AliXiStar
// author:
//        Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//



class TH1F;
class TH1D;
class TH2D;
class TH3D;
class TProfile;

class AliESDEvent;
class AliAODEvent;
class AliESDtrackCuts;
class AliESDpid;

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODPid.h"
#include "AliESDpid.h"
#include "AliXiStarEventCollection.h"

class AliXiStar : public AliAnalysisTaskSE {
 public:
 
  AliXiStar();
  AliXiStar(const char *name, Bool_t AODdecision, Bool_t MCdecision, Int_t CutListOption=0);
  virtual ~AliXiStar();
  AliXiStar(const AliXiStar &obj); 
  AliXiStar &operator=(const AliXiStar &obj);

 private:
 
  virtual void   UserCreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  

  void XiStarInit();// initialization of fixed values
  Double_t LinearPropagateToDCA(AliESDtrack*, AliESDtrack*, Double_t);// for linear propagation
  Double_t Det(Double_t, Double_t, Double_t, Double_t) const;// for linear propagation
  Double_t Det(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t) const;// for linear propagation


 
  enum {
    kNbinsM              = 300,// mult bins for certain histograms
    kXiStarCode          = 3324,// Xi(1530)^0 MC code
    kXiCode              = 3312,// Xi- MC code
    kLambdaCode          = 3122,// Lambda MC code
    kProtonCode          = 2212,// Proton+ MC code
    kPionCode            = 211// Pion+ MC code
  };
  
  const char* fname;// name of class
  AliAODEvent            *fAOD; //!    // AOD object
  AliESDEvent            *fESD; //!    // ESD object
  TList                  *fOutputList; //! Compact Output list
  AliESDtrackCuts        *fTrackCut; //! ESD track cuts
  AliPIDResponse         *fPIDResponse; //! PID object
  
   
  AliXiStarEventCollection ***fEC; //! The event collection 
  AliXiStarEventStruct *fEvt; //! The current event type
  AliXiStarTrackStruct *fTempStruct; //! A temporary track storage.  Eventually put into fEvt
  
  //

  Int_t fZvertexBins;// number of Z-vertex bins for event-mixing
  Int_t fEventsToMix;// number of maximum events to mix
  Int_t fMultBins;// number of multiplicity bins for event-mixing
  Int_t fMultLimits[11+1];// the multiplicity edges of the mult bins
  Bool_t fMCcase;// switch for MC data or real data
  Bool_t fAODcase;// switch for AODs or ESDs
  Int_t fEventCounter;// The event counter
 
  // cut list data members
  Float_t fSigmaCutProton;// Nsigma cut proton
  Float_t fSigmaCutPionFirst;// Nsigma cut pion from lambda
  Float_t fSigmaCutPionSecond;// Nsigma cut pion from Xi
  Float_t fSigmaCutPionThird;// Nsigma cut pion from Xi(1530)
  Float_t fDCAVtxProton;// dca to Primary Vertex of proton
  Float_t fDCAVtxPionFirst;// dca to Primary Vertex of pion from lambda
  Float_t fDCAVtxPionSecond;// dca to Primary Vertex of pion from Xi
  Float_t fDCAVtxLambda;// dca to Primary Vertex of Lambda
  Float_t fDCAProtonPion;// dca of lambda daughters to each other
  Float_t fDCALambdaPion;// dca of Xi daughters to each other
  Float_t fLambdaDecayLengthXY;// decay length of Lambda in xy
  Float_t fXiDecayLengthXY;// decay length of Xi in xy
  Float_t fMaxDecayLength;// max decay length
  Float_t fLamCosTheta;// cosine of pointing angle for Lambda
  Float_t fXiCosTheta;// cosine of pointing angle for Xi
  Float_t fXiStarCosTheta;// cosine of pointing angle for XiStar
  Float_t fMassWindow;// Mass window of acceptance for Lambda and Xi candidates
  
  Double_t fCovMatrix[21];// Covarience matrix of track
  Double_t fTrueMassPr, fTrueMassPi, fTrueMassK, fTrueMassLam, fTrueMassXi;// The PDG mass values
  
  
  AliESDtrack* fESDTrack4; //! esdtrack for XiStar's daughter pion
  AliESDtrack* fXiTrack; //! esdtrack for XiStar's daughter Xi
  
  Int_t fCutList;// Cut List option (mean values or systematic variations)
  
  ClassDef(AliXiStar, 1); 
};

#endif
