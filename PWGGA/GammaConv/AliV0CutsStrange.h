#ifndef ALIV0CUTSSTRANGE_H
#define ALIV0CUTSSTRANGE_H

#include "AliAODpidUtil.h"
#include "AliConversionPhotonBase.h"
#include "AliAODConversionMother.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliMCEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisCuts.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisManager.h"


class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class AliPIDResponse;
class AliKFVertex;
class TH1F;
class TH2F;
class TF1;
class TProfile;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;
class AliAODMCParticle;


class AliV0CutsStrange : public AliAnalysisCuts {
  
public: 
  enum cutIds {
    kv0FinderType,
    kclsTPCCut,
    kProtonPIDcut,
    kPionPIDcut,
    kNCuts
  };
  
   enum V0Cuts {
        kV0In=0,
        kOnFly,
        kNoTracks,
        kNoV0,
        kSameSign,
        kdEdxCuts,
        kConvPointFail,
        kPhotonCuts,
        kEventPlane,
        kV0Out
    };
  
  AliV0CutsStrange(const char *name="V0ReaderCutsStrange", const char * title="V0ReaderCutsStrange");
  AliV0CutsStrange(const AliV0CutsStrange&);
  AliV0CutsStrange& operator=(const AliV0CutsStrange&);

  virtual ~AliV0CutsStrange();

  Bool_t SetCutIds(TString cutString); 
  Int_t fCuts[kNCuts];
  Bool_t SetCut(cutIds cutID, Int_t cut);
  Bool_t UpdateCutString();
  
  static const char * fgkCutNames[kNCuts];
  
  Bool_t InitializeCutsFromCutString(const TString analysisCutSelection);   
    
  Bool_t InitPIDResponse();
  void SetPIDResponse(AliPIDResponse * pidResponse) {fPIDResponse = pidResponse;}
  AliPIDResponse * GetPIDResponse() { return fPIDResponse;}
  
  TString GetCutNumber();
  
  void PrintCuts();
  void PrintCutsWithValues();
  
  void InitCutHistograms(TString name="",Bool_t preCut = kTRUE);
  void SetFillCutHistograms(TString name="",Bool_t preCut = kTRUE){if(!fHistograms){InitCutHistograms(name,preCut);};}
  TList *GetCutHistograms(){return fHistograms;}
  void FillV0CutIndex(Int_t v0cut){if(fHistoCutIndex)fHistoCutIndex->Fill(v0cut);}  
  
  void SetV0ReaderName(TString name){fV0ReaderStrangeName = name; return;}
  
  Bool_t PhotonIsSelectedMC(TParticle *particle,AliMCEvent *mcEvent,Bool_t checkForConvertedGamma);
  
  AliVTrack * GetTrack(AliVEvent * event, Int_t label);
  AliESDtrack *GetESDTrack(AliESDEvent * event, Int_t label);
  
  Bool_t GetPIDpion(AliVTrack *fCurrentTrack);
  Bool_t GetPIDproton(AliVTrack *fCurrentTrack);
  
  // Set Individual Cuts
  Bool_t SetV0Finder(Int_t v0FinderType);
  Bool_t SetTPCClusterCut(Int_t clsTPCCut);
  Bool_t SetProtonPIDCut(Int_t pPIDcut);
  Bool_t SetPionPIDCut(Int_t pPIDcut);
  
  Bool_t SelectV0Finder(Bool_t onfly){
    if(onfly == fUseOnFlyV0Finder) return kTRUE;
    else return kFALSE;
  }
  Int_t GetV0FinderSameSign(){return fUseOnFlyV0FinderSameSign;}

  void SetIsQA(Bool_t isQA){fIsQA=isQA;} 
  
  virtual   Bool_t IsSelected(TObject* /*obj*/) {return kTRUE;}
  virtual   Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
  
protected:
  TList*            fHistograms;                          //
  AliPIDResponse*   fPIDResponse;                         //
  
  TString           fV0ReaderStrangeName;
  TObjString*       fCutString;                           // cut number used for analysis
  TString           fCutStringRead;
  
  Bool_t            fIsQA;
  
  //cuts
  
  Bool_t            fUseOnFlyV0Finder;                    // flag
  Int_t             fUseOnFlyV0FinderSameSign;            // int to set same sign pairing
  Double_t          fUseCorrectedTPCClsInfo;
  Double_t          fMinClsTPCToF;
  Double_t          fMinClsTPC;
  Double_t          fPIDTPCnSigmaProtonLow;
  Double_t          fPIDTPCnSigmaProtonUp;
  Double_t          fPIDTOFnSigmaProtonLow;
  Double_t          fPIDTOFnSigmaProtonUp;
  Double_t          fPIDTPCnSigmaPionLow;
  Double_t          fPIDTPCnSigmaPionUp;
  Double_t          fPIDTOFnSigmaPionLow;
  Double_t          fPIDTOFnSigmaPionUp;
  
  
  
  // Histograms
  TH1F*             fHistoCutIndex;                       // bookkeeping for cuts
  TH2F*             fHistodEdxCutsProton;                 // bookkeeping proton ID
  TH2F*             fHistoTPCdEdxProtonBefore;
  TH2F*             fHistoTPCdEdxSigmaProtonBefore;
  TH2F*             fHistoTPCdEdxProtonAfter;
  TH2F*             fHistoTPCdEdxSigmaProtonAfter;
  TH2F*             fHistoTOFdEdxProtonBefore;
  TH2F*             fHistoTOFdEdxSigmaProtonBefore;
  TH2F*             fHistoTOFdEdxProtonAfter;
  TH2F*             fHistoTOFdEdxSigmaProtonAfter;
  TH2F*             fHistodEdxCutsPion;                   // bookkeeping pion ID
  TH2F*             fHistoTPCdEdxPionBefore;
  TH2F*             fHistoTPCdEdxSigmaPionBefore;
  TH2F*             fHistoTPCdEdxPionAfter;
  TH2F*             fHistoTPCdEdxSigmaPionAfter;    
  TH2F*             fHistoTOFdEdxPionBefore;
  TH2F*             fHistoTOFdEdxSigmaPionBefore;
  TH2F*             fHistoTOFdEdxPionAfter;
  TH2F*             fHistoTOFdEdxSigmaPionAfter;  
  
  Bool_t            fPreSelCut;                           // Flag for preselection cut used in V0Reader
  
private:
  
  ClassDef(AliV0CutsStrange,2)
};

#endif
