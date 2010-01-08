#ifndef ALIANALYSISTASKSDDRP
#define ALIANALYSISTASKSDDRP

class TList;
class TH1F;
class TTree;
class TString;
class AliESDEvent;
class AliESDfriend;
class AliITSresponseSDD;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSDDRP : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskSDDRP();
  virtual ~AliAnalysisTaskSDDRP();
  virtual void   UserExec(Option_t *option);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t *option);

  void SetRunNumber(Int_t nrun){
    fRunNumber=nrun;
  }
  void SetMinITSPoints(Int_t minp=3){
    fMinITSpts=minp;
  }
  void SetUseOnlyCINT1BTriggers(Bool_t use=kTRUE){
    fOnlyCINT1BTrig=use;
  }
  void SetMinPfordEdx(Float_t minp=1.5){
    fMinPfordEdx=minp;
  }
  Bool_t CheckModule(Int_t lay, Int_t lad, Int_t det) const;


 private:
  AliAnalysisTaskSDDRP(const AliAnalysisTaskSDDRP &source);
  AliAnalysisTaskSDDRP& operator=(const AliAnalysisTaskSDDRP &source);
  
  TList*  fOutput;          //! ntuple with output of vertexers
  TH1F*   fHistNEvents;     //! histo with N of events  
  TH1F*   fRecPMod;         //! histo with module occupancy (RecP) 
  TH1F*   fTrackPMod;       //! histo with module occupancy (TrP)
  TH1F*   fGoodAnMod;       //! histo good anodes per module 
  TH1F*   fRecPLadLay3;     //! histo with ladder occupancy on layer3 (RecP) 
  TH1F*   fRecPLadLay4;     //! histo with ladder occupancy on layer4 (RecP)
  TH1F*   fTrackPLadLay3;   //! histo with ladder occupancy on layer3 (TrP)
  TH1F*   fTrackPLadLay4;   //! histo with ladder occupancy on layer4 (TrP)
  TH1F*   fGoodAnLadLay3;   //! histo good anodes per ladder on layer3 
  TH1F*   fGoodAnLadLay4;   //! histo good anodes per ladder on layer4 
  TH1F*   fDriftTimeRP;     //! histo with drift time distribution (RecP)
  TH1F*   fDriftTimeTP;     //! histo with drift time distribution (TrP)
  TH1F*   fSignalTime[8];   //! histos of dE/dx in time windows
  AliESDEvent  *fESD;       // ESD object
  AliESDfriend *fESDfriend; // ESD friend object
  AliITSresponseSDD* fResp; // ResponseSDD object
  Int_t   fRunNumber;       // Run number
  Int_t   fMinITSpts;       // Minimum number of points per track
  Float_t fMinPfordEdx;     // Minimum momentum for dE/dx
  Bool_t  fOnlyCINT1BTrig;  // Flag for using all events or only intections
  Bool_t  fInitialised;     // True if initialised
  ClassDef(AliAnalysisTaskSDDRP,1);  
};


#endif
