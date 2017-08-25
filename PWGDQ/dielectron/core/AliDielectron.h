#ifndef ALIDIELECTRON_H
#define ALIDIELECTRON_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#####################################################
//#                                                   #
//#             Class AliDielectron                   #
//#         Main Class for e+e- analysis              #
//#                                                   #
//#  by WooJin J. Park, GSI / W.J.Park@gsi.de         #
//#     Ionut C. Arsene, GSI / I.C.Arsene@gsi.de      #
//#     Magnus Mager, CERN / Magnus.Mager@cern.ch     #
//#     Jens Wiechula, Uni HD / Jens.Wiechula@cern.ch #
//#                                                   #
//#####################################################


#include <TNamed.h>
#include <TObjArray.h>
#include <THnBase.h>
#include <TSpline.h>

#include <AliAnalysisFilter.h>
#include <AliKFParticle.h>

#include "AliDielectronHistos.h"
#include "AliDielectronHF.h"
#include "AliDielectronCutQA.h"

class AliEventplane;
class AliVEvent;
class AliMCEvent;
class THashList;
class AliDielectronCF;
class AliDielectronDebugTree;
class AliDielectronTrackRotator;
class AliDielectronPair;
class AliDielectronSignalMC;
class AliDielectronMixingHandler;

//________________________________________________________________
class AliDielectron : public TNamed {

  friend class AliDielectronMixingHandler; //mixing as friend class

public:
  enum EPairType { kEv1PP=0, kEv1PM, kEv1MM,
      kEv1PEv2P, kEv1MEv2P, kEv2PP,
      kEv1PEv2M, kEv1MEv2M, kEv2PM,
      kEv2MM, kEv1PMRot };
  enum ELegType  { kEv1P, kEv1M, kEv2P, kEv2M };

  AliDielectron();
  AliDielectron(const char* name, const char* title);
  virtual ~AliDielectron();

  void Init();

  void Process(/*AliVEvent *ev1, */TObjArray *arr);
  Bool_t Process(AliVEvent *ev1, AliVEvent *ev2=0);

  const AliAnalysisFilter& GetEventFilter() const { return fEventFilter; }
  const AliAnalysisFilter& GetTrackFilter() const { return fTrackFilter; }
  const AliAnalysisFilter& GetPairFilter()  const { return fPairFilter;  }

  AliAnalysisFilter& GetEventFilter()       { return fEventFilter;       }
  AliAnalysisFilter& GetTrackFilter()       { return fTrackFilter;       }
  AliAnalysisFilter& GetPairFilter()        { return fPairFilter;        }
  AliAnalysisFilter& GetPairPreFilter()     { return fPairPreFilter1;     }
  AliAnalysisFilter& GetPairPreFilter2()     { return fPairPreFilter2;     }
  AliAnalysisFilter& GetPairPreFilterLegs() { return fPairPreFilterLegs1; }
  AliAnalysisFilter& GetPairPreFilterLegs2() { return fPairPreFilterLegs2; }
  AliAnalysisFilter& GetEventPlanePreFilter(){ return fEventPlanePreFilter; }
  AliAnalysisFilter& GetEventPlanePOIPreFilter(){ return fEventPlanePOIPreFilter; }
  AliDielectronQnEPcorrection* GetQnTPCACcuts(){ return fQnTPCACcuts; }

  void  SetMotherPdg( Int_t pdgMother ) { fPdgMother=pdgMother; }
  void  SetLegPdg(Int_t pdgLeg1, Int_t pdgLeg2) { fPdgLeg1=pdgLeg1; fPdgLeg2=pdgLeg2; }
  Int_t GetMotherPdg() const { return fPdgMother; }
  Int_t GetLeg1Pdg()   const { return fPdgLeg1;   }
  Int_t GetLeg2Pdg()   const { return fPdgLeg2;   }

  void SetCutQA(Bool_t qa=kTRUE) { fCutQA=qa; }
  void SetNoPairing(Bool_t noPairing=kTRUE) { fNoPairing=noPairing; }
  void SetProcessLS(Bool_t doLS=kTRUE) { fProcessLS=doLS; }
  void SetUseKF(Bool_t useKF=kTRUE) { fUseKF=useKF; }
  const TObjArray* GetTrackArray(Int_t i) const {return (i>=0&&i<4)?&fTracks[i]:0;}
  const TObjArray* GetPairArray(Int_t i)  const {return (i>=0&&i<11)?
      static_cast<TObjArray*>(fPairCandidates->UncheckedAt(i)):0;}

  TObjArray** GetPairArraysPointer() { return &fPairCandidates; }
  void SetPairArraysPointer( TObjArray *arr) { fPairCandidates=arr; }
  void SetHistogramArray(AliDielectronHF * const histoarray) { fHistoArray=histoarray; }
  const TObjArray * GetHistogramArray() const { return fHistoArray?fHistoArray->GetHistArray():0x0; }
  const TObjArray * GetQAHistArray() const { return fQAmonitor?fQAmonitor->GetQAHistArray():0x0; }

  void SetHistogramManager(AliDielectronHistos * const histos) { fHistos=histos; }
  AliDielectronHistos* GetHistoManager() const { return fHistos; }
  const THashList * GetHistogramList() const { return fHistos?fHistos->GetHistogramList():0x0; }

  Bool_t HasCandidates() const { return GetPairArray(1)?GetPairArray(1)->GetEntriesFast()>0:0; }
  Bool_t HasCandidatesLikeSign() const {
    return (GetPairArray(0)&&GetPairArray(2)) ? (GetPairArray(0)->GetEntriesFast()>0 || GetPairArray(2)->GetEntriesFast()>0) : 0;
  }

  Bool_t HasCandidatesTR() const {return GetPairArray(10)?GetPairArray(10)->GetEntriesFast()>0:0;}
  void SetCFManagerPair(AliDielectronCF * const cf) { fCfManagerPair=cf; }
  AliDielectronCF* GetCFManagerPair() const { return fCfManagerPair; }

  void SetPreFilterEventPlane(Bool_t setValue=kTRUE){fPreFilterEventPlane=setValue;};
  void SetLikeSignSubEvents(Bool_t setValue=kTRUE){fLikeSignSubEvents=setValue;};

  void SetPreFilterUnlikeOnly(Bool_t setValue=kTRUE){fPreFilterUnlikeOnly1=setValue;};
  void SetPreFilterLikeOnly(Bool_t setValue=kTRUE){fPreFilterLikeOnly1=setValue;};
  void SetPreFilterAllSigns(Bool_t setValue=kTRUE){fPreFilterAllSigns1=setValue;};
  void SetPreFilterPhotons(Bool_t setValue=kTRUE){fPreFilterPhotons1=setValue;};
  void SetPreFilterOnlyOnePair(Bool_t setValue=kTRUE){fPreFilterOnlyOnePair1=setValue;};

  void SetPreFilterUnlikeOnly2(Bool_t setValue=kTRUE){fPreFilterUnlikeOnly2=setValue;};
  void SetPreFilterLikeOnly2(Bool_t setValue=kTRUE){fPreFilterLikeOnly2=setValue;};
  void SetPreFilterAllSigns2(Bool_t setValue=kTRUE){fPreFilterAllSigns2=setValue;};
  void SetPreFilterPhotons2(Bool_t setValue=kTRUE){fPreFilterPhotons2=setValue;};
  void SetPreFilterOnlyOnePair2(Bool_t setValue=kTRUE){fPreFilterOnlyOnePair2=setValue;};

  void SetTrackRotator(AliDielectronTrackRotator * const rot) { fTrackRotator=rot; }
  AliDielectronTrackRotator* GetTrackRotator() const { return fTrackRotator; }

  void SetRotatePP(Bool_t const rotate){ fRotatePP = rotate;}
  void SetRotateMM(Bool_t const rotate){ fRotateMM = rotate;}


  Bool_t GetRotatePP() const { return fRotatePP; }
  Bool_t GetRotateMM() const { return fRotateMM; }


  void SetMixingHandler(AliDielectronMixingHandler *mix) { fMixing=mix; }
  AliDielectronMixingHandler* GetMixingHandler() const { return fMixing; }

  void SetHasMC(Bool_t hasMC) { fHasMC = hasMC; }
  Bool_t GetHasMC() const     { return fHasMC;  }

  void SetStoreRotatedPairs(Bool_t storeTR) {fStoreRotatedPairs = storeTR;}
  void SetDontClearArrays(Bool_t dontClearArrays=kTRUE) { fDontClearArrays=dontClearArrays; }
  Bool_t DontClearArrays() const { return fDontClearArrays; }

  void AddSignalMC(AliDielectronSignalMC* signal);

  void SetDebugTree(AliDielectronDebugTree * const tree) { fDebugTree=tree; }

  const TObjArray* GetMCSignals() const { return fSignalsMC; }
  static const char* TrackClassName(Int_t i) { return (i>=0&&i<4)?fgkTrackClassNames[i]:""; }
  static const char* PairClassName(Int_t i)  { return (i>=0&&i<11)?fgkPairClassNames[i]:""; }

  void SetEstimatorFilename(const Char_t* filename) {fEstimatorFilename = filename;} // Does not work on the Grid
  void SetEstimatorObjArray(TObjArray* array) {fEstimatorObjArray = array;} // Should work on the Grid
  void SetTRDcorrectionFilename(const Char_t* filename) {fTRDpidCorrectionFilename = filename;}
  void SetVZEROCalibrationFilename(const Char_t* filename) {fVZEROCalibrationFilename = filename;}
  void SetVZERORecenteringFilename(const Char_t* filename) {fVZERORecenteringFilename = filename;}
  void SetZDCRecenteringFilename(const Char_t* filename) {fZDCRecenteringFilename = filename;}
  void InitLegEffMap(TString filename, TString generatedname="hGenerated", TString foundname="hFound")  { fLegEffMap=InitEffMap(filename,generatedname,foundname)  ;}
  void InitPairEffMap(TString filename, TString generatedname="hGenerated", TString foundname="hFound") { fPairEffMap=InitEffMap(filename,generatedname,foundname) ;}

  void SetCentroidCorrArr(TObjArray *arrFun, Bool_t bHisto, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetWidthCorrArr(TObjArray *arrFun, Bool_t bHisto, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetCentroidCorrFunction(TF1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetWidthCorrFunction(TF1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetCentroidCorrFunction(TH1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetWidthCorrFunction(TH1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  /// @TODO: implement smart way to use generic functions for all detectors.
  void SetCentroidCorrFunctionITS(TF1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetWidthCorrFunctionITS(TF1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetCentroidCorrFunctionITS(TH1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetWidthCorrFunctionITS(TH1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetCentroidCorrFunctionTOF(TF1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetWidthCorrFunctionTOF(TF1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetCentroidCorrFunctionTOF(TH1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);
  void SetWidthCorrFunctionTOF(TH1 *fun, UInt_t varx, UInt_t vary=0, UInt_t varz=0);

  void SetQnTPCACcuts(AliDielectronQnEPcorrection *acCuts){ fQnTPCACcuts = acCuts; fACremovalIsSetted = kTRUE;}
  void SetQnVectorNormalisation(TString norm) {fQnVectorNorm = norm;}
  void SaveDebugTree();
  Bool_t DoEventProcess() const { return fEventProcess; }
  void SetEventProcess(Bool_t setValue=kTRUE) { fEventProcess=setValue; }
  Bool_t GammaTracksUsed() const { return fUseGammaTracks; }
  void SetUseGammaTracks(Bool_t setValue=kTRUE) { fUseGammaTracks=setValue; }
  void  FillHistogramsFromPairArray(Bool_t pairInfoOnly=kFALSE);

private:

  Bool_t fCutQA;                    // monitor cuts
  AliDielectronCutQA *fQAmonitor;   // monitoring of cuts

  TObjArray *fPostPIDCntrdCorrArr;   // post pid correction array containing the run-wise objects for electron sigma centroids in TPC
  TObjArray *fPostPIDWdthCorrArr;    // post pid correction array containing the run-wise objects for electron sigma widths in TPC
  TH1 *fPostPIDCntrdCorr;   // post pid correction object for electron sigma centroids in TPC
  TH1 *fPostPIDWdthCorr;    // post pid correction object for electron sigma widths in TPC
  TH1 *fPostPIDCntrdCorrITS;// post pid correction object for electron sigma centroids in ITS
  TH1 *fPostPIDWdthCorrITS; // post pid correction object for electron sigma widths in ITS
  TH1 *fPostPIDCntrdCorrTOF;// post pid correction object for electron sigma centroids in TOF
  TH1 *fPostPIDWdthCorrTOF; // post pid correction object for electron sigma widths in TOF
                            /// @TODO: smart way to store corrections, potentially for all detectors.
  TObject *fLegEffMap;      // single electron efficiency map
  TObject *fPairEffMap;      // pair efficiency map
  AliAnalysisFilter fEventFilter;    // Event cuts
  AliAnalysisFilter fTrackFilter;    // leg cuts
  AliAnalysisFilter fPairPreFilter1;  // pair prefilter cuts
  AliAnalysisFilter fPairPreFilter2;  // pair prefilter cuts
  AliAnalysisFilter fPairPreFilterLegs1; // Leg filter after the pair prefilter cuts
  AliAnalysisFilter fPairPreFilterLegs2; // Leg filter after the pair prefilter cuts
  AliAnalysisFilter fPairFilter;     // pair cuts
  AliAnalysisFilter fEventPlanePreFilter;  // event plane prefilter cuts
  AliAnalysisFilter fEventPlanePOIPreFilter;  // PoI cuts in the event plane prefilter
  AliDielectronQnEPcorrection *fQnTPCACcuts; // QnFramework est. 2016 ac removal
  TString fQnVectorNorm;

  Int_t fPdgMother;     // pdg code of mother tracks
  Int_t fPdgLeg1;       // pdg code leg1
  Int_t fPdgLeg2;       // pdg code leg2

  TObjArray* fSignalsMC;      // array of AliDielectronSignalMC

  Bool_t fNoPairing;    // if to skip pairing, can be used for track QA only
  Bool_t fProcessLS; // do the like-sign pairing (default kTRUE)
  Bool_t fUseKF;    // if to skip pairing, can be used for track QA only

  AliDielectronHF *fHistoArray;   // Histogram framework
  AliDielectronHistos *fHistos;   // Histogram manager
                                  //  Streaming and merging should be handled
                                  //  by the analysis framework
  TBits *fUsedVars;               // used variables

  TObjArray fTracks[4];           //! Selected track candidates
                                  //  0: Event1, positive particles
                                  //  1: Event1, negative particles
                                  //  2: Event2, positive particles
                                  //  3: Event2, negative particles

  TObjArray *fPairCandidates;     //! Pair candidate arrays
                                  //TODO: better way to store it? TClonesArray?

  AliDielectronCF *fCfManagerPair;//Correction Framework Manager for the Pair
  AliDielectronTrackRotator *fTrackRotator; //Track rotator
  Bool_t fRotatePP; // combine rotated positive tracks
  Bool_t fRotateMM; // combine rotated negative tracks
  AliDielectronDebugTree *fDebugTree;  // Debug tree output
  AliDielectronMixingHandler *fMixing; // handler for event mixing

  Bool_t fPreFilterEventPlane;  //Filter for the Eventplane determination in TPC
  Bool_t fACremovalIsSetted;
  Bool_t fLikeSignSubEvents;    //Option for dividing into subevents, sub1 ++ sub2 --
  Bool_t fPreFilterUnlikeOnly1;  //Apply PreFilter either in +- or to ++/--/+- individually
  Bool_t fPreFilterLikeOnly1;    //Apply PreFilter either in -- and ++ or to ++/--/+- individually
  Bool_t fPreFilterAllSigns1;    //Apply PreFilter find in  ++/--/+- and remove from all
  Bool_t fPreFilterPhotons1;    //Apply PreFilter to search for photons
  Bool_t fPreFilterOnlyOnePair1; // each track can only belong to one pair in the prefilter
  Bool_t fPreFilterUnlikeOnly2;  //Apply PreFilter either in +- or to ++/--/+- individually
  Bool_t fPreFilterLikeOnly2;    //Apply PreFilter either in -- and ++ or to ++/--/+- individually
  Bool_t fPreFilterAllSigns2;    //Apply PreFilter find in  ++/--/+- and remove from all
  Bool_t fPreFilterPhotons2;    //Apply PreFilter to search for photons
  Bool_t fPreFilterOnlyOnePair2; // each track can only belong to one pair in the prefilter
  Bool_t fHasMC;                //If we run with MC, at the moment only needed in AOD
  Bool_t fStoreRotatedPairs;    //It the rotated pairs should be stored in the pair array
  Bool_t fDontClearArrays;      //Don't clear the arrays at the end of the Process function, needed for external use of pair and tracks
  Bool_t fEventProcess;         //Process event (or pair array)
  Bool_t fUseGammaTracks;       // use function SetGammaTracks for MCtruth photons

  void FillTrackArrays(AliVEvent * const ev, Int_t eventNr=0);
  void EventPlanePreFilter(Int_t arr1, Int_t arr2, TObjArray arrTracks1, TObjArray arrTracks2, const AliVEvent *ev);
  void PairPreFilter(Int_t arr1, Int_t arr2, TObjArray &arrTracks1, TObjArray &arrTracks2, const AliVEvent *ev, Int_t prefilterN);
  void FillPairArrays(Int_t arr1, Int_t arr2, const AliVEvent *ev = 0x0);
  void FillPairArrayTR();

  Int_t GetPairIndex(Int_t arr1, Int_t arr2) const {return arr1>=arr2?arr1*(arr1+1)/2+arr2:arr2*(arr2+1)/2+arr1;}

  void InitPairCandidateArrays();
  void ClearArrays();

  TObjArray* PairArray(Int_t i);
  TObject* InitEffMap(TString filename, TString generatedname, TString foundname);

  static const char* fgkTrackClassNames[4];   //Names for track arrays
  static const char* fgkPairClassNames[11];   //Names for pair arrays

  TString fEstimatorFilename;                // name for the pp multiplicity estimators filename
  TObjArray* fEstimatorObjArray;		     // Grid compatible version of the above
  TString fTRDpidCorrectionFilename;         // name for the file containing the single particle TRD pid corrections
  TString fVZEROCalibrationFilename;         // file containing VZERO channel-by-channel calibration
  TString fVZERORecenteringFilename;         // file containing VZERO Q-vector recentering averages
  TString fZDCRecenteringFilename;         // file containing ZDCQ-vector recentering averages

  void ProcessMC(AliVEvent *ev1);

  void  FillHistograms(const AliVEvent *ev, Bool_t pairInfoOnly=kFALSE);
  void  FillMCHistograms(const AliVEvent *ev);
  void  FillMCHistograms(Int_t label1, Int_t label2, Int_t nSignal);
  void  FillHistogramsMC(const AliMCEvent *ev,  AliVEvent *ev1);
  void  FillHistogramsPair(AliDielectronPair *pair,Bool_t fromPreFilter=kFALSE);
  void  FillHistogramsTracks(TObjArray **tracks);

  void  FillDebugTree();

  AliDielectron(const AliDielectron &c);
  AliDielectron &operator=(const AliDielectron &c);

  ClassDef(AliDielectron,17);
};

inline void AliDielectron::InitPairCandidateArrays()
{
  //
  // initialise all pair candidate arrays
  //
  fPairCandidates->SetOwner();
  for (Int_t i=0;i<11;++i){
    TObjArray *arr=new TObjArray;
    fPairCandidates->AddAt(arr,i);
    arr->SetOwner();
  }
}

inline TObjArray* AliDielectron::PairArray(Int_t i)
{
  //
  // for internal use only: unchecked return of track array for fast access
  //
  return static_cast<TObjArray*>(fPairCandidates->UncheckedAt(i));
}

inline void AliDielectron::ClearArrays()
{
  //
  // Reset the Arrays
  //
  for (Int_t i=0;i<4;++i){
    fTracks[i].Clear();
  }
  for (Int_t i=0;i<11;++i){
    if (PairArray(i)) PairArray(i)->Delete();
  }
}

#endif
