#ifndef ALISPECTRAAODEVENTCUTS_H
#define ALISPECTRAAODEVENTCUTS_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliSpectraAODEventCuts
//
//
//
//
// Authors: Michele Floris, CERN, Philip Versteeg, UU, Redmer Bertens, UU
//-------------------------------------------------------------------------

class AliAODEvent;
class AliSpectraAODTrackCuts;
class AliSpectraAODHistoManager;

#include "TNamed.h"

class AliSpectraAODEventCuts : public TNamed
{
 public:
  enum {  kProcessedEvents = 0,kPhysSelEvents,kAcceptedEvents, kVtxRange, kVtxCentral, kVtxNoEvent, kNVtxCuts};

  // Constructors
 AliSpectraAODEventCuts() : TNamed(), fAOD(0), fIsMC(0), fUseCentPatchAOD049(0), fIsSelected(0), fCentralityCutMin(0), fCentralityCutMax(0), fHistoCuts(0),fHistoVtxBefSel(0),fHistoVtxAftSel(0),fHistoEtaBefSel(0),fHistoEtaAftSel(0),fHistoNChAftSel(0) {}
  AliSpectraAODEventCuts(const char *name);
  virtual  ~AliSpectraAODEventCuts() {}
  
  void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
  Bool_t GetIsMC()           const           { return fIsMC;};
  
  void SetUseCentPatchAOD049(Bool_t useCentPatchAOD049 = kFALSE)    {fUseCentPatchAOD049 = useCentPatchAOD049; };
  Bool_t GetUseCentPatchAOD049()           const           { return fUseCentPatchAOD049;};
  
  // Methods
  Bool_t IsSelected(AliAODEvent * aod,AliSpectraAODTrackCuts     *trackcuts);
  Bool_t CheckVtxRange();
  Bool_t CheckCentralityCut();
  void  SetCentralityCutMin(Float_t cut)  { fCentralityCutMin = cut; }
  void  SetCentralityCutMax(Float_t cut)  { fCentralityCutMax = cut; }

   
  TH1I * GetHistoCuts()         {  return fHistoCuts; }
  TH1F * GetHistoVtxBefSel()         {  return fHistoVtxBefSel; }
  TH1F * GetHistoVtxAftSel()         {  return fHistoVtxAftSel; }
  TH1F * GetHistoEtaBefSel()         {  return fHistoEtaBefSel; }
  TH1F * GetHistoEtaAftSel()         {  return fHistoEtaAftSel; }
  TH1F * GetHistoNChAftSel()         {  return fHistoNChAftSel; }
  Float_t  GetCentralityMin()  const {  return fCentralityCutMin; }
  Float_t  GetCentralityMax()  const {  return fCentralityCutMax; }
  void   PrintCuts();
  Double_t ApplyCentralityPatchAOD049();

  Float_t  NumberOfEvents()     { return fHistoCuts->GetBinContent(kAcceptedEvents+1); }
  Float_t  NumberOfProcessedEvents()     { return fHistoCuts->GetBinContent(kProcessedEvents+1); }
  Float_t  NumberOfPhysSelEvents()     { return fHistoCuts->GetBinContent(kPhysSelEvents+1); }

  Long64_t Merge(TCollection* list);


 private:
  
  AliAODEvent     *fAOD;              //! AOD event
  Bool_t          fIsMC;// true if processing MC
  Bool_t          fUseCentPatchAOD049;// Patch for centrality selection on AOD049
  AliSpectraAODTrackCuts     *fTrackCuts;              //! track cuts
  Bool_t          fIsSelected;        // True if cuts are selected
  Float_t         fCentralityCutMin;     // minimum centrality percentile
  Float_t         fCentralityCutMax;     // maximum centrality percentile
  TH1I            *fHistoCuts;        // Cuts statistics
  TH1F            *fHistoVtxBefSel;        // Vtx distr before event selection
  TH1F            *fHistoVtxAftSel;        // Vtx distr after event selection
  TH1F            *fHistoEtaBefSel;        // Eta distr before event selection
  TH1F            *fHistoEtaAftSel;        // Eta distr after event selection
  TH1F            *fHistoNChAftSel;        // NCh distr after event selection
  AliSpectraAODEventCuts(const AliSpectraAODEventCuts&);
  AliSpectraAODEventCuts& operator=(const AliSpectraAODEventCuts&);
  
  ClassDef(AliSpectraAODEventCuts, 2);
  
};
#endif

