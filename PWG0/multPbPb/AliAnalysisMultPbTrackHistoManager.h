#ifndef ALIANALYSISMULTPBTRACKHISTOMANAGER_H
#define ALIANALYSISMULTPBTRACKHISTOMANAGER_H

#include "AliHistoListWrapper.h"
#include "TH3D.h"

class TH3D;
class TH1D;
class TH1I;

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

  typedef enum {kHistoGen, kHistoRec, kHistoRecPrim, kHistoRecSecWeak, kHistoRecSecMat, kHistoRecFake, kNHistos} Histo_t;
  typedef enum {kStatAll, kStatCentr, kStatVtx, kNStatBins} Stat_t;

  AliAnalysisMultPbTrackHistoManager();
  AliAnalysisMultPbTrackHistoManager(const char * name,const char * title);
  AliAnalysisMultPbTrackHistoManager(const AliAnalysisMultPbTrackHistoManager& obj) ;
  ~AliAnalysisMultPbTrackHistoManager();
  
  // Setters
  void SetSuffix(const char * suffix) { fHNameSuffix = suffix;}

  // Histo getters
  TH3D * GetHistoPtEtaVz(Histo_t id);
  TH1D * GetHistoPt (Histo_t id, Float_t minEta = -22222, Float_t maxEta = -22222, Float_t minVz  = -22222, Float_t maxVz  = -22222, Bool_t scaleWidth = kTRUE);
  TH1D * GetHistoEta(Histo_t id, Float_t minPt  = -22222, Float_t maxPt  = -22222, Float_t minVz  = -22222, Float_t maxVz  = -22222, Bool_t scaleWidth = kTRUE);
  TH1D * GetHistoVz (Histo_t id, Float_t minPt  = -22222, Float_t maxPt  = -22222, Float_t minEta = -22222, Float_t maxEta = -22222, Bool_t scaleWidth = kTRUE);

  TH1I * GetHistoStats();
  TH1D * GetHistoDCA(Histo_t id);

  // Misch utils
  void ScaleHistos (Double_t nev, Option_t * option="");
  


  // Histo bookers
  TH3D * BookHistoPtEtaVz(const char * name, const char * title);
  TH1D * BookHistoDCA(const char * name, const char * title);
  TH1I * BookHistoStats();

  // 
  TH1 * GetHisto(const char * name);

private:

  static const char * kStatStepNames[];       // names of the step hist
  static const char * kHistoPtEtaVzNames[];   // names of the 3D histograms pt/eta/vz
  static const char * kHistoDCANames[];   // names of the DCA histograms 
  TString fHNameSuffix; // Suffix added to all histo names. Useful if you have in the same session e.g. MC and data.

  AliAnalysisMultPbTrackHistoManager& operator=(const AliAnalysisMultPbTrackHistoManager& task);
  
  ClassDef(AliAnalysisMultPbTrackHistoManager, 2)


};

#endif
