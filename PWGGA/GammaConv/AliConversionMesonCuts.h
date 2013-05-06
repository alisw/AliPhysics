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
	kBackgroundScheme,            
	kNumberOfBGEvents,           
	kDegreesForRotationMethod,    
	kRapidityMesonCut,    
	kRCut,
	kalphaMesonCut,               
	kchi2MesonCut,
	kElecShare,
	kToCloseV0s,
	kuseMCPSmearing,
	kNCuts
  };

  Bool_t SetCutIds(TString cutString); 
  Int_t fCuts[kNCuts];
  Bool_t SetCut(cutIds cutID, Int_t cut);
  Bool_t UpdateCutString();

  static const char * fgkCutNames[kNCuts];

  Bool_t InitializeCutsFromCutString(const TString analysisCutSelection);
  void FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0);
   
  AliConversionMesonCuts(const char *name="MesonCuts", const char * title="Meson Cuts");
  AliConversionMesonCuts(const AliConversionMesonCuts&);
  AliConversionMesonCuts& operator=(const AliConversionMesonCuts&);

  virtual ~AliConversionMesonCuts();                            //virtual destructor

  virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
  virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

  TString GetCutNumber();

  // Cut Selection
  Bool_t MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal=kTRUE, Double_t fRapidityShift=0.);
  Bool_t MesonIsSelectedMC(TParticle *fMCMother,AliStack *fMCStack, Double_t fRapidityShift=0.);
  Bool_t MesonIsSelectedMCDalitz(TParticle *fMCMother,AliStack *fMCStack, Double_t fRapidityShift=0.);
  Bool_t MesonIsSelectedMCChiC(TParticle *fMCMother,AliStack *fMCStack, Int_t &, Int_t &, Int_t &, Double_t fRapidityShift=0. );
  void PrintCuts();

  void InitCutHistograms(TString name="");
  void SetFillCutHistograms(TString name=""){if(!fHistograms){InitCutHistograms(name);};}
  TList *GetCutHistograms(){return fHistograms;}
  void SmearParticle(AliAODConversionPhoton * photon);
  
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
  Bool_t SetNumberOfBGEvents(Int_t NumberOfBGEvents);
  Bool_t SetMCPSmearing(Int_t useMCPSmearing);
  Bool_t SetSharedElectronCut(Int_t sharedElec);
  Bool_t SetToCloseV0sCut(Int_t toClose);


  // Request Flags
  Bool_t UseRotationMethod(){return fUseRotationMethodInBG;}
  Bool_t UseTrackMultiplicity(){return fUseTrackMultiplicityForBG;}
  Int_t GetNumberOfBGEvents(){return fNumberOfBGEvents;}
  Int_t NDegreesRotation(){return fnDegreeRotationPMForBG;}
  Bool_t DoBGCalculation(){return fDoBG;}
  Bool_t DoBGProbability(){return fdoBGProbability;}
  Bool_t UseElecSharingCut(){return fDoSharedElecCut;}
  Bool_t UseToCloseV0sCut(){return fDoToCloseV0sCut;}
  Bool_t UseMCPSmearing(){return fUseMCPSmearing;}
  Int_t BackgroundHandlerType(){return fBackgroundHandler;}
  
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
  Bool_t fDoBG; // flag to intialize BG
  Bool_t fdoBGProbability; // flag to use probability method for meson bg estimation
  Bool_t fUseTrackMultiplicityForBG; // flag to use track multiplicity for meson bg estimation (else V0 mult)
  Int_t fnDegreeRotationPMForBG; //
  Int_t fNumberOfBGEvents; //
  Float_t fOpeningAngle; // min opening angle for meson
  Bool_t fDoToCloseV0sCut; //
  Double_t fminV0Dist; //
  Bool_t fDoSharedElecCut; //
  Bool_t fUseMCPSmearing; // flag
  Double_t fPBremSmearing;//
  Double_t fPSigSmearing; //
  Double_t fPSigSmearingCte; //
  TF1 *fBrem; //
  TRandom3 fRandom; //
  Int_t fElectronLabelArraySize;
  Int_t *fElectronLabelArray; //[fElectronLabelArraySize] Array with elec/pos v0 label
  Int_t fBackgroundHandler; //
  
  // Histograms
  TObjString *fCutString; // cut number used for analysis
  TH1F *hMesonCuts; // bookkeeping for meson cuts
  TH1F *hMesonBGCuts; // bookkeeping for meson bg cuts


private:


  ClassDef(AliConversionMesonCuts,3)
};


#endif
