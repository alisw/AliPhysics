#ifndef ALIVERTEXINGHFUTILS_H
#define ALIVERTEXINGHFUTILS_H


/* $Id$ */

////////////////////////////////////////////////////////////////////
/// \class AliVertexingHFUtils
///                                                               //
/// \brief Class with functions useful for different D2H analyses //
/// - event plane resolution                                      //
/// - <pt> calculation with side band subtraction                 //
/// - tracklet multiplicity calculation                            //
/// \author Origin: F.Prino, Torino, prino@to.infn.it              //
///                                                               //
////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "TVirtualPad.h"
#include "AliHFInvMassFitter.h"

#include <vector>

#include "Fit/Fitter.h"
#include "Fit/Chi2FCN.h"
#include "Math/WrappedMultiTF1.h"
#include "Fit/BinData.h"
#include "HFitInterface.h"

class AliMCEvent;
class AliMCParticle;
class AliAODMCParticle;
class AliAODMCHeader;
class AliGenEventHeader;
class AliAODEvent;
class TProfile;
class TClonesArray;
class TH1F;
class TH2F;
class TF1;

using std::vector;

class AliVertexingHFUtils : public TObject{
 public:
  AliVertexingHFUtils();
  AliVertexingHFUtils(Int_t k);
  virtual ~AliVertexingHFUtils() {};

  /// Significance calculator
  static void ComputeSignificance(Double_t signal, Double_t  errsignal, Double_t  background, Double_t  errbackground, Double_t &significance,Double_t &errsignificance);

  /// Functions for Event plane resolution
  void SetK(Int_t k){fK=k;}
  void SetSubEvResol(Double_t res){fSubRes=res;}
  void SetSubEventHisto(const TH1F* hSub){
    fSubRes=GetSubEvResol(hSub);
  }
  Int_t GetK() const  {return fK;}
  Double_t GetSubEvResol() const  {return fSubRes;}
  Double_t Pol(Double_t x) const {return Pol(x,fK);}
  Double_t FindChi() const {return FindChi(fSubRes,fK);}
  Double_t GetFullEvResol() const {return GetFullEvResol(fSubRes,fK);}
  static Double_t FindChi(Double_t res,  Int_t k=1);
  static Double_t Pol(Double_t x, Int_t k);
  static Double_t ResolK1(Double_t x);
  static Double_t GetSubEvResol(const TH1F* hSubEvCorr){
    if(hSubEvCorr) return TMath::Sqrt(hSubEvCorr->GetMean());
    else return 1.;
  }
  static Double_t GetSubEvResolLowLim(const TH1F* hSubEvCorr){
    if(hSubEvCorr) return TMath::Sqrt(hSubEvCorr->GetMean()-hSubEvCorr->GetMeanError());
    else return 1.;
  }
  static Double_t GetSubEvResolHighLim(const TH1F* hSubEvCorr){
    if(hSubEvCorr) return TMath::Sqrt(hSubEvCorr->GetMean()+hSubEvCorr->GetMeanError());
    else return 1.;
  }
  static Double_t GetFullEvResol(Double_t resSub, Int_t k=1);
  static Double_t GetFullEvResol(const TH1F* hSubEvCorr, Int_t k=1);
  static Double_t GetFullEvResolLowLim(const TH1F* hSubEvCorr, Int_t k=1);
  static Double_t GetFullEvResolHighLim(const TH1F* hSubEvCorr, Int_t k=1);
  static TString  GetGenerator(Int_t label, AliAODMCHeader* header);
  static Bool_t IsTrackInjected(Int_t label,AliAODMCHeader *header,TClonesArray *arrayMC);
  static Bool_t IsTrackInjected(AliAODTrack *track,AliAODMCHeader *header,TClonesArray *arrayMC);
  static void GetTrackPrimaryGenerator(AliAODTrack *track,AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen);
  static void GetTrackPrimaryGenerator(Int_t label,AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen);
  static Bool_t IsCandidateInjected(AliAODRecoDecayHF *cand, AliAODMCHeader *header,TClonesArray *arrayMC);
  static Bool_t IsCandidateInjected(AliAODRecoDecayHF *cand, AliAODEvent* aod, AliAODMCHeader *header,TClonesArray *arrayMC);
  static Bool_t HasCascadeCandidateAnyDaughInjected(AliAODRecoCascadeHF *cand, AliAODMCHeader *header,TClonesArray *arrayMC);
  static Int_t PreSelectITSUpgrade(TClonesArray* arrayMC, AliAODMCHeader *header, TObjArray aodTracks, Int_t nDaug, Int_t pdgabs, const Int_t *pdgDg);

  /// Functions for tracklet multiplcity calculation
  void SetEtaRangeForTracklets(Double_t mineta, Double_t maxeta){
    fMinEtaForTracklets=mineta;
    fMaxEtaForTracklets=maxeta;
  }
  static Int_t GetNumberOfTrackletsInEtaRange(AliAODEvent* ev, Double_t mineta, Double_t maxeta);
  Int_t GetNumberOfTrackletsInEtaRange(AliAODEvent* ev) const {
    return GetNumberOfTrackletsInEtaRange(ev,fMinEtaForTracklets,fMaxEtaForTracklets);
  }
  static Int_t GetGeneratedMultiplicityInEtaRange(TClonesArray* arrayMC, Double_t mineta, Double_t maxeta);
  static Int_t GetGeneratedPrimariesInEtaRange(TClonesArray* arrayMC, Double_t mineta, Double_t maxeta);
  static Int_t GetGeneratedPhysicalPrimariesInEtaRange(TClonesArray* arrayMC, Double_t mineta, Double_t maxeta);

  /// Functions for event shape variables
  static void GetSpherocity(AliAODEvent* aod,
                            Double_t &spherocity, Double_t &phiRef,
                            Double_t etaMin=-0.8, Double_t etaMax=0.8,
                            Double_t ptMin=0.15, Double_t ptMax=10.,
                            Int_t filtbit1=256, Int_t filtbit2=512,
                            Int_t minMult=3, Double_t phiStepSizeDeg=0.1,
                            Int_t nTrksToSkip=0, Int_t* idToSkip=0x0);

  static void GetGeneratedSpherocity(TClonesArray *arrayMC,
                                     Double_t &spherocity, Double_t &phiRef,
                                     Double_t etaMin=-0.8, Double_t etaMax=0.8,
                                     Double_t ptMin=0.15, Double_t ptMax=10.,
                                     Int_t minMult=3, Double_t phiStepSizeDeg=0.1);

  static Double_t GetSphericity(AliAODEvent* aod,
                                Double_t etaMin=-0.8, Double_t etaMax=0.8,
                                Double_t ptMin=0.15, Double_t ptMax=10.,
                                Int_t filtbit1=256, Int_t filtbit2=512,
                                Int_t minMult=3);

  /// Utilities for V0 multiplicity checks
  static Double_t GetVZEROAEqualizedMultiplicity(AliAODEvent* ev);
  static Double_t GetVZEROCEqualizedMultiplicity(AliAODEvent* ev);

  /// Functions for computing average pt
  static void AveragePt(Float_t& averagePt, Float_t& errorPt, Float_t ptmin, Float_t ptmax, TH2F* hMassD, Float_t massFromFit, Float_t sigmaFromFit, 
                        TF1* funcB2, Float_t sigmaRangeForSig=2.5, Float_t sigmaRangeForBkg=4.5, Float_t minMass=0., Float_t maxMass=3., Int_t rebin=1);

  /// Functions for processing trigger information
  static Bool_t CheckT0TriggerFired(AliAODEvent* aodEv);

  /// Rebinning of invariant mass histograms
  static TH1D* RebinHisto(TH1* hOrig, Int_t reb, Int_t firstUse=-1);
  static TH1* AdaptTemplateRangeAndBinning(const TH1 *hRef,TH1 *hData, Double_t minFit, Double_t maxFit);

  /// Functions for computing true impact parameter of D meson
  static Double_t GetTrueImpactParameterDzero(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partDp);
  static Double_t GetTrueImpactParameterDplus(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partDp);

  static Double_t GetCorrectedNtracklets(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult);

  /// Functions to check the decay tree
  static Int_t CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Bool_t searchUpToQuark=kTRUE);
  static Int_t CheckOrigin(AliMCEvent* mcEvent, AliMCParticle *mcPart, Bool_t searchUpToQuark=kTRUE);
  static Bool_t IsTrackFromCharm(AliAODTrack* tr, TClonesArray* arrayMC);
  static Bool_t IsTrackFromBeauty(AliAODTrack* tr, TClonesArray* arrayMC);
  static Bool_t IsTrackFromHadronDecay(Int_t pdgMoth, AliAODTrack* tr, TClonesArray* arrayMC);
  static Double_t GetBeautyMotherPt(TClonesArray* arrayMC, AliAODMCParticle *mcPart);
  static Double_t GetBeautyMotherPtAndPDG(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t &pdgGranma);
  static Int_t CheckD0Decay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckD0Decay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);
  static Int_t CheckDplusDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckDplusDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);
  static Int_t CheckDplusKKpiDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckDplusKKpiDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);
  static Int_t CheckDplusK0spiDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckDsDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckDsDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);
  static Int_t CheckDsK0sKDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckDstarDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckDstarDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);
  static Int_t CheckLcpKpiDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckLcpKpiDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);
  static Int_t CheckLcV0bachelorDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckLcV0bachelorDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);
  static Int_t CheckXicXipipiDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckBplusDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckBplusDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);
  static Int_t CheckB0toDminuspiDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckB0toDminuspiDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);
  static Int_t CheckBsDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab, Bool_t ITS2UpgradeProd=kFALSE);
  static Int_t CheckBsDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab, Bool_t ITS2UpgradeProd=kFALSE);
  static Int_t CheckLbDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab);
  static Int_t CheckLbDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);


  /// Simultaneus fit
  /// GlobalChi2 structure for simultaneus in-plane - out-of-plane fit
  struct GlobalInOutOfPlaneChi2 {
    GlobalInOutOfPlaneChi2(ROOT::Math::IMultiGenFunction & fInPlane, ROOT::Math::IMultiGenFunction & fOutOfPlane, Int_t npars, vector<UInt_t> commonpars) :
    fChi2_InPlane(&fInPlane),
    fChi2_OutOfPlane(&fOutOfPlane),
    fNpars(npars),
    fCommonPars() {
      fCommonPars.clear();
      fCommonPars = commonpars;
    }

    Double_t operator() (const Double_t *par) const {
      const UInt_t npars = fNpars;
      Double_t pInPlane[npars];
      for(UInt_t iPar=0; iPar<npars; iPar++)
        pInPlane[iPar]=par[iPar];

      UInt_t iParOutOfPlane = fNpars;
      Double_t pOutOfPlane[npars];
      vector<UInt_t> veccopy = fCommonPars;
      vector<UInt_t>::iterator iter;
      for(UInt_t iPar=0; iPar<npars; iPar++) {
        iter = find(veccopy.begin(),veccopy.end(),iPar);
        if(iter!=veccopy.end()) { //is common
          pOutOfPlane[iPar] = par[iPar];
        }
        else {
          pOutOfPlane[iPar] = par[iParOutOfPlane];
          iParOutOfPlane++;
        }
      }

      return (*fChi2_InPlane)(pInPlane) + (*fChi2_OutOfPlane)(pOutOfPlane);
    }

    const ROOT::Math::IMultiGenFunction *fChi2_InPlane;
    const ROOT::Math::IMultiGenFunction *fChi2_OutOfPlane;
    UInt_t fNpars;
    vector<UInt_t> fCommonPars;
  };

  static ROOT::Fit::FitResult DoInPlaneOutOfPlaneSimultaneusFit(AliHFInvMassFitter *&massfitterInPlane, AliHFInvMassFitter *&massfitterOutOfPlane, 
                                                                TH1F* hMassInPlane, TH1F* hMassOutOfPlane, Double_t MinMass, Double_t MaxMass, 
                                                                Double_t massD, vector<UInt_t> commonpars);
  
  /// Helper functions for D-meson analyses
  static Double_t ComputeMaxd0MeasMinusExp(AliAODRecoDecayHF *cand, Double_t bfield);
  static Double_t CombineNsigmaTPCTOF(Double_t nsigmaTPC, Double_t nsigmaTOF);

 private:

  Int_t fK;             /// ratio of measured harmonic to event plane harmonic
  Double_t fSubRes;     /// sub-event resolution = sqrt(<cos[n(phiA-phiB)] >)
  Double_t fMinEtaForTracklets; /// min eta for counting tracklets
  Double_t fMaxEtaForTracklets; /// min eta for counting tracklets

  /// \cond CLASSIMP
  ClassDef(AliVertexingHFUtils,0);
  /// \endcond
};
#endif
