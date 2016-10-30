#ifndef AliAnalysisTaskSigma0Run2_H
#define AliAnalysisTaskSigma0Run2_H

// ROOT includes
#include <TList.h>
#include <TClonesArray.h>

// AliROOT includes
#include "AliV0ReaderV1.h"
#include "AliV0ReaderStrange.h"
#include <AliAnalysisTaskSE.h>

// forward delcarations
class AliESDtrackCuts;
class AliVParticle;
class AliVTrack;
class AliPIDResponse;

class AliAnalysisTaskSigma0Run2: public AliAnalysisTaskSE {
public:
  AliAnalysisTaskSigma0Run2();
  AliAnalysisTaskSigma0Run2(const char* name);
  virtual ~AliAnalysisTaskSigma0Run2();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  
  void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
  void SetV0ReaderStrangeName(TString name){fV0ReaderStrangeName=name; return;}
  void SetIsHeavyIon(Bool_t isHeavyIon){fIsHeavyIon=isHeavyIon;}  
  void SetIsMC(Bool_t isMC){fIsMC=isMC;}  
  void SetIsQA(Bool_t isQA){fIsQA=isQA;}  
  
  void ProcessMCParticles();
  
private:
  
  AliV0ReaderV1       *fV0Reader;             // V0Reader for basic conversion photon 
  AliV0ReaderStrange  *fV0ReaderStrange;      // V0Reader for basic conversion photon selection
  TString             fV0ReaderName;     
  TString             fV0ReaderStrangeName;          
  TClonesArray        *fReaderGammas;         // array with photons from fV0Reader
  TClonesArray        *fReaderV0s;
  AliESDEvent         *fESDEvent;
  
  Bool_t              fIsMC;
  Bool_t              fIsHeavyIon;
  Bool_t              fIsQA;
  
  AliMCEvent          *fMCEvent;    // current MC event
  AliStack            *fMCStack;    // current MC stack
  
  TList               *fOutputContainer;
  TList               *fReco;
  TList               *fQA;
  TList               *fMC;
  TList               *fMCtruth;
  
  //Reco Histos
  TH1F                *fHistLambdaInvMass;
  TH2F                *fHistLambdaInvMassPt;
  TH2F                *fHistLambdaInvMassEta;
  TH2F                *fHistLambdaAngle;
  TH2F                *fHistLambdaR;
  TH1F                *fHistPhotonPt;
  TH2F                *fHistPhotonInvMassPt;
  TH2F                *fHistPhotonInvMassEta;
  TH2F                *fHistPhotonAngle;
  TH2F                *fHistPhotonR;
  TH1F                *fHistSigmaInvMass;
  TH2F                *fHistSigmaInvMassPt;
  TH2F                *fHistSigmaInvMassEta;
  TH2F                *fHistSigmaAngle;
  TH2F                *fHistSigmaR;
    
  //MC histos
  TH1F                *fHistPtMC;
  
  
  //MC truth histos
  
  
  //QA histos
  TH1F                *fHistNevents;
  TH1F                *fHistNTracks;
  TH2F                *fHistNTracksPt;
  TH1F                *fHistNV0;
  TH1F                *fHistPt;
  TH1F                *fHistZvertex;
  
  AliESDtrackCuts     *fESDCuts;              //! track cuts for ESD analysis
  AliPIDResponse      *fPIDResponse;          //! pid response
  
  AliAnalysisTaskSigma0Run2(const AliAnalysisTaskSigma0Run2 &task);
  AliAnalysisTaskSigma0Run2& operator= (const AliAnalysisTaskSigma0Run2 &task);
    
  ClassDef(AliAnalysisTaskSigma0Run2, 0);
};
#endif
