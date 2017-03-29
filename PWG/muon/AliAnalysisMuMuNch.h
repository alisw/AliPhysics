#ifndef ALIANALYSISMUMUNCH_H
#define ALIANALYSISMUMUNCH_H

/**
 *
 * \class AliAnalysisMuMuNch
 * \brief Charged particle multiplicity analysis (for plotting muon variables against it)
 * \author L. Aphecetche and J. Martin Blanco (Subatech)
 */

#include "AliAnalysisMuMuBase.h"
#include "TRandom3.h"
#include "TMap.h"

class TH2F;
class TH2;
class AliVEvent;
class AliAODTracklets;
class TAxis;
class TF1;
class AliAODVertex;
class TObjArray;

class AliAnalysisMuMuNch : public AliAnalysisMuMuBase
{
public:

  AliAnalysisMuMuNch(TRootIOCtor* ioCtor);

  AliAnalysisMuMuNch(TH2F* spdCorrection=0x0, TProfile* spdMeanCorrection = 0x0, Double_t meanTrRef=-1., Double_t etaMin=-0.5, Double_t etaMax=0.5
                     , Double_t zmin=-10., Double_t zmax=10., Bool_t disableHistos=kFALSE ,Bool_t computeResolution=kFALSE);

  AliAnalysisMuMuNch(TH2F* spdCorrection, TH2F* spdFluctuations, Double_t etaMin=-0.5, Double_t etaMax=0.5
                                         , Double_t zMin=-10., Double_t zMax=10., Bool_t disableHistos=kFALSE, Bool_t computeResolution=kFALSE);

  AliAnalysisMuMuNch(TProfile* spdMeanCorrection, TProfile* spdMeanCorrectionToCompare, Double_t meanTrRef=-1., Double_t etaMin=-1.0,
                     Double_t etaMax=1.0, Double_t zMin=-10., Double_t zMax=10., Double_t etaMinToCompare=-0.5, Double_t etaMaxToCompare=0.5,Bool_t disableHistos=kFALSE, Bool_t computeResolution=kFALSE);

  AliAnalysisMuMuNch(TObjArray* spdCorrectionList, Double_t meanTrRef=-1., TObjArray* runWeightsList=0x0, Double_t etaMin=-0.5, Double_t etaMax=0.5
                                         , Double_t zmin=-10, Double_t zmax=10, Bool_t disableHistos=kFALSE ,Bool_t computeResolution=kFALSE);

  AliAnalysisMuMuNch(const char* V0side/*="V0A*/, TProfile* V0MeanCorrection=0x0, Double_t meanV0mRef=-1., TProfile* spdMeanCorrection=0x0,
                     Double_t meanTrRef=-1., Double_t etaMin=-1., Double_t etaMax=1., Double_t zmin=-10, Double_t zmax=10, Bool_t disableHistos=kFALSE);

  virtual ~AliAnalysisMuMuNch();

  void DefineGeneratorName(const char* genName);

  Bool_t HasAtLeastNTrackletsInEtaRange(const AliVEvent& event, Int_t n,
                                        Double_t& etaMin, Double_t& etaMax) const;

  void NameOfHasAtLeastNTrackletsInEtaRange(TString& name, Int_t n,
                                            Double_t& etaMin, Double_t& etaMax) const;

  virtual void SetEvent(AliVEvent* event, AliMCEvent* mcEvent=0x0);

  void SetRun(const AliInputEventHandler* eventHandler);

  virtual void Terminate(Option_t *);

protected:

  void AddHisto(const char* eventSelection,
                const char* triggerClassName,
                const char* centrality,
                const char* histoname,
                Double_t z,
                TH1* h,
                Bool_t isMC=kFALSE);

  void AttachSPDAcceptance(UInt_t dataType,
                           const char* eventSelection,
                           const char* triggerClassName,
                           const char* centrality,const char* histoname);

  Double_t ApplyFluctuations(Double_t ntrCorr);

  void FillHistosForMCEvent(const char* eventSelection, const char* triggerClassName,
              const char* centrality);


  void FillHistosForEvent(const char* eventSelection, const char* triggerClassName,
                          const char* centrality);

  virtual void DefineHistogramCollection(const char* eventSelection, const char* triggerClassName,
                                         const char* centrality, Bool_t mix = kFALSE);

  void DefineSPDAcceptance();

  void DefineSPDFluctuationsMap(TH2F* spdFluctuations);

  Bool_t GetEtaRangeSPD(Double_t spdZVertex, Double_t etaRange[]);

  Double_t GetSPDCorrection(Double_t zvert, Double_t eta) const;

  Double_t GetTrackletsMeanCorrection(Double_t zvert, Int_t nTracklets,Bool_t corrToCompare=kFALSE) const;

  Double_t GetV0MeanCorrection(Double_t zvert, Int_t mV0) const;

  AliAODTracklets* GetTracklets(AliVEvent* event);

  Bool_t IsMCtrackFromGenerator(Int_t indexMC) const;

private:

  Double_t NumberOfTrackletsInEtaRange(const AliVEvent& event, Double_t& etamin,
                                       Double_t& etamax, Bool_t corrected=kFALSE) const;

  void DefineSPDCorrectionMap(TObjArray* spdCorrectionList);

private:
  TH2F* fSPDOneOverAccxEff; // Nch/Tracklets_SPD (eta vs z). SPD AccxEffCorrection for tracklets
  TObjArray* fSPDFluctuationsList; // Array for the fluctuations distributions in Ntr corrected by SPD AccxEff slices
  TMap* fSPDFluctuationsMap; // Map for the fluctuations in Ntr corrected by AccxEff
  TMap* fSPDCorrectionMap; // Map for the SPD AccxEffCorrections or Mean Tracklets corrections for different subperiods
  TObjArray* fSPDCorrectionList; // List of SPD AccxEffCorrections or Mean Tracklets corrections for different subperiods
  TProfile* fSPDMeanTracklets; // <Tracklets_SPD> vs zvtx. 'mean' correction for tracklets
  TProfile* fSPDMeanTrackletsCorrToCompare; // <Tracklets_SPD> vs zvtx. 'mean' correction for tracklets in the range for comparison with the main range
  TProfile* fV0MeanMult; // <V0Multiplicty> vs zvtx. 'mean' correction for V0 multiplicity
  TAxis* fEtaAxis; // Eta axis used for the histos
  TAxis* fZAxis;  // Z vertex axis used for the histos
  AliVEvent* fCurrentEvent; //! cache of the current event
  Double_t fMeanTrRef; // Mean reference number of tracklets for mean tracklets correction
  Double_t fMeanV0Ref;  // Mean reference V0 multiplicity for mean V0 multiplicity correction
  Double_t fEtaMin; // Minimum tracklet eta value
  Double_t fEtaMax; // Maximum tracklet eta value
  Double_t fEtaMinToCompare; // Minimum tracklet eta value for the comparison
  Double_t fEtaMaxToCompare; // Maximum tracklet eta value for the comparison
  Double_t fetaRange[2];
  Double_t fZMin; // Minimum z vertex value
  Double_t fZMax; // Maximum z vertex value
  Bool_t fResolution; // Flag to set the resolution computation
  TRandom3* frand;
  TString* fGeneratorHeaderClass; // Class of the header MC generator

  TF1* fSPD1LR; // SPD acceptance shape
  TF1* fSPD1LL; // SPD acceptance shape
  TF1* fSPD2LR; // SPD acceptance shape
  TF1* fSPD2LL; // SPD acceptance shape

  TObjArray* fMCWeightList; // List with the run weights
  Double_t fMCWeight; // Weight of current MC run
  TString* fV0side; // Which V0 side will be use to estimate multiplicicty

  ClassDef(AliAnalysisMuMuNch,7) // implementation of AliAnalysisMuMuBase for Nch analysis
};

#endif
