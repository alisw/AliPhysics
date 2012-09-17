#ifndef ALICONVERSIONMESONCUTS_H
#define ALICONVERSIONMESONCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: Svein Lindal, Daniel Lohner												*

#include "AliAODpidUtil.h"
#include "AliConversionPhotonBase.h"
#include "AliAODConversionMother.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliStack.h"
#include "AliAnalysisCuts.h"
#include "TH1F.h"

class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class AliKFVertex;
class TH1F;
class TH2F;
class AliPIDResponse;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;


using namespace std;

class AliConversionMesonCuts : public AliAnalysisCuts {
	
 public: 


  enum cutIds {
	kMesonKind,
	kchi2MesonCut,
	kalphaMesonCut,               
	kRapidityMesonCut,    
	kRCut,
	kBackgroundScheme,            
	kDegreesForRotationMethod,    
	kNumberOfRotations,           
	kElecShare,
	kToCloseV0s,
	kNCuts
  };

  Bool_t SetCutIds(TString cutString); 
  Int_t fCuts[kNCuts];
  Bool_t SetCut(cutIds cutID, Int_t cut);
  Bool_t UpdateCutString(cutIds cutID, Int_t value);

  static const char * fgkCutNames[kNCuts];

  Bool_t InitializeCutsFromCutString(const TString analysisCutSelection);
  void FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0);
   
  AliConversionMesonCuts(const char *name="MesonCuts", const char * title="Meson Cuts");
  virtual ~AliConversionMesonCuts();                            //virtual destructor

  virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
  virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

  TString GetCutNumber();

  // Cut Selection
  Bool_t MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal=kTRUE);
  Bool_t MesonIsSelectedMC(TParticle *fMCMother,AliStack *fMCStack, Bool_t bMCDaughtersInAcceptance=kFALSE);
	Bool_t MesonIsSelectedMCDalitz(TParticle *fMCMother,AliStack *fMCStack);
  void PrintCuts();

  void InitCutHistograms(TString name="",Bool_t preCut = kTRUE);
  void SetFillCutHistograms(TString name="",Bool_t preCut = kTRUE){if(!fHistograms){InitCutHistograms(name,preCut);};}
  TList *GetCutHistograms(){return fHistograms;}

  ///Cut functions
  Bool_t RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s);
  Bool_t RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0);
  // Set Individual Cuts
  Bool_t SetRCut(Int_t RCut);
  Bool_t SetMesonKind(Int_t mesonKind);
  Bool_t SetChi2MesonCut(Int_t chi2MesonCut);
  Bool_t SetAlphaMesonCut(Int_t alphaMesonCut);
  Bool_t SetRapidityMesonCut(Int_t RapidityMesonCut);
  Bool_t SetBackgroundScheme(Int_t BackgroundScheme);
  Bool_t SetNDegreesForRotationMethod(Int_t DegreesForRotationMethod);
  Bool_t SetNumberOfRotations(Int_t NumberOfRotations);
  Bool_t SetSharedElectronCut(Int_t sharedElec);
  Bool_t SetToCloseV0sCut(Int_t toClose);


  // Request Flags
  Bool_t UseRotationMethod(){return fUseRotationMethodInBG;}
  Int_t NumberOfRotationEvents(){return fnumberOfRotationEventsForBG;}
  Int_t NDegreesRotation(){return fnDegreeRotationPMForBG;}
  Bool_t DoBGProbability(){return fdoBGProbability;}
  Bool_t UseElecSharingCut(){return fDoSharedElecCut;}
  Bool_t UseToCloseV0sCut(){return fDoToCloseV0sCut;}

  protected:
  TList *fHistograms;
  //cuts
  Int_t fMesonKind;
  Double_t fMaxR; //r cut  
  Double_t fChi2CutMeson; //chicut meson
  Double_t fAlphaMinCutMeson; // min value for meson alpha cut
  Double_t fAlphaCutMeson; // max value for meson alpha cut
  Double_t fRapidityCutMeson; // max value for meson rapidity
  Bool_t fUseRotationMethodInBG; // flag to apply rotation method for meson bg estimation
  Bool_t fdoBGProbability; // flag to use probability method for meson bg estimation
  Bool_t fUseTrackMultiplicityForBG; // flag to use track multiplicity for meson bg estimation (else V0 mult)
  Int_t fnDegreeRotationPMForBG; //
  Int_t fnumberOfRotationEventsForBG; //
  Float_t fOpeningAngle; // min opening angle for meson
  Bool_t fDoToCloseV0sCut; //
  Double_t fminV0Dist; //
  Bool_t fDoSharedElecCut; //
  TRandom3 fRandom; //
  Int_t *fElectronLabelArray; // Array with elec/pos v0 label
  
  // Histograms
  TObjString *fCutString; // cut number used for analysis
  TH1F *hMesonCuts; // bookkeeping for meson cuts
  TH1F *hMesonBGCuts; // bookkeeping for meson bg cuts


private:

  AliConversionMesonCuts(const AliConversionMesonCuts&); // not implemented
  AliConversionMesonCuts& operator=(const AliConversionMesonCuts&); // not implemented


  ClassDef(AliConversionMesonCuts,2)
};


#endif
