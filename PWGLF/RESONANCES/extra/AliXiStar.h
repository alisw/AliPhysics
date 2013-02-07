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
    kPionCode            = 211,// Pion+ MC code
    kNCutVariations      = 13,// number of cut variations
    kNCuts               = 15// number of cut types
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
  Float_t fMaxDecayLength;// max decay length
  Float_t fMassWindow;// Mass window of acceptance for Lambda and Xi candidates
  
  Double_t fCovMatrix[21];// Covarience matrix of track
  Double_t fTrueMassPr, fTrueMassPi, fTrueMassK, fTrueMassLam, fTrueMassXi;// The PDG mass values
  
  /////////////////////////////////////

  struct St_CutType {
    TH3F *fXi; //!
    TH3F *fXibar; //!
    //    
    TH3F *fXiMinusPiPlus; //!
    TH3F *fXiMinusPiMinus; //!
    TH3F *fXiPlusPiPlus; //!
    TH3F *fXiPlusPiMinus; //!
    
    TH3F *fXiMinusPiPlusbkg; //!
    TH3F *fXiMinusPiMinusbkg; //!
    TH3F *fXiPlusPiPlusbkg; //!
    TH3F *fXiPlusPiMinusbkg; //!
    //
    TH3F *fMCrecXi; //!
    TH3F *fMCrecXibar; //!
        
    TH3F *fMCrecXiMinusPiPlus; //!
    TH3F *fMCrecXiPlusPiMinus; //!
    //    
    /*TH2F *fMCrecXiStarxiy; //!
    TH2F *fMCrecXiStarpiony; //!
    TH2F *fMCrecXilambday; //!
    TH2F *fMCrecXipiony; //!
    TH2F *fMCrecLamprotony; //!
    TH2F *fMCrecLampiony; //!*/
  };
  struct St_CutType CutVar[kNCutVariations]; //!
  
  
  AliESDtrack* fESDTrack4; //! esdtrack for XiStar's daughter pion
  AliESDtrack* fXiTrack; //! esdtrack for XiStar's daughter Xi
  
  Int_t fCutList;// Cut List option (mean values or systematic variations)

  Float_t fDecayParameters[kNCuts];// array of reconstruction kinematics
  Float_t fCutValues[kNCutVariations][kNCuts];// array of reconstruction kinematics

  ClassDef(AliXiStar, 1); 
};

#endif
