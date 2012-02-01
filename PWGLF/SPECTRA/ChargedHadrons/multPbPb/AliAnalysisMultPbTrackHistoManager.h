#ifndef ALIANALYSISMULTPBTRACKHISTOMANAGER_H
#define ALIANALYSISMULTPBTRACKHISTOMANAGER_H

#include "AliHistoListWrapper.h"
#include "TH3D.h"

class TH3D;
class TH1D;
class TH1I;
class AliMCParticle;
class TMap;

//#define WEIGHTED_DCA
#define TRANSVERSE_DCA


//-------------------------------------------------------------------------
//                      AliAnalysisMultPbTrackHistoManager
// 
// 
//
//
// Author: Michele Floris, CERN
//-------------------------------------------------------------------------


class AliAnalysisMultPbTrackHistoManager : public AliHistoListWrapper {

public:

  typedef enum {kHistoGen, kHistoRec, kHistoRecPrim, kHistoRecSecWeak, kHistoRecSecMat, kHistoRecFake, kHistoRecHighestMeanPt, kHistoRecLowestMeanPt, kNHistos} Histo_t;
  typedef enum {kStatAll,  kStatCentr, kStatPhysSel, kStatVtx, kStatZDCCut, kStatVtxRangeCut, kNStatBins} Stat_t;
  typedef enum {kPartPiPlus, kPartKPlus, kPartP, kPartLPlus, kPartPiMinus, kPartKMinus, kPartPBar, kPartLMinus, kPartOther, kNPart} Part_t;

  // these bits define the behaviour at merging
  // if kKeepMax(Mib)Mean is set, the histogram with the highest (lowest) mean is kept
  enum {kKeepMaxMean = BIT(22), kKeepMinMean = BIT(23)};



  AliAnalysisMultPbTrackHistoManager();
  AliAnalysisMultPbTrackHistoManager(const char * name,const char * title);
  AliAnalysisMultPbTrackHistoManager(const AliAnalysisMultPbTrackHistoManager& obj) ;
  ~AliAnalysisMultPbTrackHistoManager();
  
  // Setters
  void SetSuffix(const char * suffix) { fHNameSuffix = suffix;}

  // Histo getters
  TH3D * GetHistoPtEtaVz(Histo_t id, Int_t particle = -1);
  TH1D * GetHistoPt (Histo_t id, Float_t minEta = -22222, Float_t maxEta = -22222, Float_t minVz  = -22222, Float_t maxVz  = -22222, Bool_t scaleWidth = kTRUE);
  TH1D * GetHistoEta(Histo_t id, Float_t minPt  = -22222, Float_t maxPt  = -22222, Float_t minVz  = -22222, Float_t maxVz  = -22222, Bool_t scaleWidth = kTRUE);
  TH1D * GetHistoVz (Histo_t id, Float_t minPt  = -22222, Float_t maxPt  = -22222, Float_t minEta = -22222, Float_t maxEta = -22222, Bool_t scaleWidth = kTRUE);
  TH2D * GetHistoPtVz (Histo_t id, Float_t minEta = -22222, Float_t maxEta = -22222, Bool_t scaleWidth = kFALSE);

  TH1I * GetHistoStats();
  TH2D * GetHistoDCA(Histo_t id);
  TH1D * GetHistoMult(Histo_t id);
  TH2D * GetHistoV0vsNtracks(Histo_t id);


  TH1D * GetHistoSpecies(Histo_t id);
  TH1D * GetHistoProcess(Histo_t id);
  TH1D * GetHistoVzEvent(Histo_t id);
  TH1D * GetHistoPtEvent(Histo_t id);
  TH1D * GetHistoMeanPt (Histo_t id);

  TH2D * GetHistoElectronCutQA();

  TMap * GetParticleSpecies(){return fParticleSpecies;}

  // Misch utils
  void ScaleHistos (Double_t nev, Option_t * option="");
  Int_t GetLocalParticleID(AliMCParticle * part);
  void FillParticleID(Histo_t id, AliMCParticle * part) { GetHistoSpecies(id)->Fill(GetLocalParticleID(part));}
  void FillSpeciesMap(Int_t pdgCode);


  // Histo bookers
  TH3D * BookHistoPtEtaVz(const char * name, const char * title);
  TH2D * BookHistoDCA(const char * name, const char * title);
  TH1I * BookHistoStats();
  TH1D * BookHistoMult(const char * name, const char * title);

  // 
  TH1 * GetHisto(const char * name);

  virtual void Print(Option_t* option = "") const {fList->Print(option);}


  Long64_t Merge(TCollection* list); // need to override this because of the mean pts

private:

  static const char * kStatStepNames[];       // names of the step hist
  static const char * kHistoPtEtaVzNames[];   // names of the 3D histograms pt/eta/vz
  static const char * kHistoDCANames[];   // names of the DCA histograms 
  static const char * kHistoPrefix[];   // prefix for histo names // FIXME: remove the others and keep only this 
  static const char * kSpeciesName[];   // Particle species
  TString fHNameSuffix; // Suffix added to all histo names. Useful if you have in the same session e.g. MC and data.
  
  TMap * fParticleSpecies;// Map assiciating PDG code to number of particles in the physical primary sample.

  AliAnalysisMultPbTrackHistoManager& operator=(const AliAnalysisMultPbTrackHistoManager& task);
  
  ClassDef(AliAnalysisMultPbTrackHistoManager, 4)


};

#endif
