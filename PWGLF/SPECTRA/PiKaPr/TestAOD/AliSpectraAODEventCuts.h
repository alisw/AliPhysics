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
//class AliSpectraAODHistoManager;

#include "TH1I.h"
#include "TNamed.h"

class AliSpectraAODEventCuts : public TNamed
{
 public:
  enum {  kProcessedEvents = 0,kPhysSelEvents,kAcceptedEvents, kVtxRange, kVtxCentral, kVtxNoEvent, kQVector, kNVtxCuts};

  // Constructors
 AliSpectraAODEventCuts() : TNamed(), fAOD(0), fTrackBits(0), fIsMC(0), fCentFromV0(0), fUseCentPatchAOD049(0),fTrackCuts(0), fIsSelected(0), fCentralityCutMin(0), fCentralityCutMax(0), fQVectorCutMin(0), fQVectorCutMax(0), fVertexCutMin(0), fVertexCutMax(0), fMultiplicityCutMin(0), fMultiplicityCutMax(0), fHistoCuts(0),fHistoVtxBefSel(0),fHistoVtxAftSel(0),fHistoEtaBefSel(0),fHistoEtaAftSel(0),fHistoNChAftSel(0),fHistoQVector(0),fHistoEP(0) {}
  AliSpectraAODEventCuts(const char *name);
  virtual  ~AliSpectraAODEventCuts() {}
  
  void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
  Bool_t GetIsMC()           const           { return fIsMC;};
  void SetCentFromV0(Bool_t centFromV0 = kFALSE)    {fCentFromV0 = centFromV0; };
  Bool_t GetCentFromV0()           const           { return fCentFromV0;};
  
  void SetUseCentPatchAOD049(Bool_t useCentPatchAOD049 = kFALSE)    {fUseCentPatchAOD049 = useCentPatchAOD049; };
  Bool_t GetUseCentPatchAOD049()           const           { return fUseCentPatchAOD049;};
  
  // Methods
  Bool_t IsSelected(AliAODEvent * aod,AliSpectraAODTrackCuts     *trackcuts);
  Bool_t CheckVtxRange();
  Bool_t CheckCentralityCut();
  Bool_t CheckMultiplicityCut();
  Bool_t CheckQVectorCut();
  void  SetCentralityCutMin(Float_t cut)  { fCentralityCutMin = cut; }
  void  SetCentralityCutMax(Float_t cut)  { fCentralityCutMax = cut; }
  //void  SetQVectorPosCut(Float_t min,Float_t max)  { fQVectorPosCutMin = min; fQVectorPosCutMax = max; }
  //void  SetQVectorNegCut(Float_t min,Float_t max)  { fQVectorNegCutMin = min; fQVectorNegCutMax = max; }
  void  SetQVectorCut(Float_t min,Float_t max)  { fQVectorCutMin = min; fQVectorCutMax = max; }
  void  SetVertexCut(Float_t min,Float_t max)  { fVertexCutMin = min; fVertexCutMax = max; }
  void  SetMultiplicityCut(Float_t min,Float_t max)  { fMultiplicityCutMin = min; fMultiplicityCutMax = max; }
  void SetTrackBits(UInt_t TrackBits) {fTrackBits=TrackBits;}

  UInt_t GetTrackType()  const    { return fTrackBits;}
  TH1I * GetHistoCuts()         {  return fHistoCuts; }
  TH1F * GetHistoVtxBefSel()         {  return fHistoVtxBefSel; }
  TH1F * GetHistoVtxAftSel()         {  return fHistoVtxAftSel; }
  TH1F * GetHistoEtaBefSel()         {  return fHistoEtaBefSel; }
  TH1F * GetHistoEtaAftSel()         {  return fHistoEtaAftSel; }
  TH1F * GetHistoNChAftSel()         {  return fHistoNChAftSel; }
  /* TH1F * GetHistoQVectorPos()         {  return fHistoQVectorPos; } */
  /* TH1F * GetHistoQVectorNeg()         {  return fHistoQVectorNeg; } */
  TH1F * GetHistoQVector()         {  return fHistoQVector; }
  TH1F * GetHistoEP()         {  return fHistoEP; }
  Float_t  GetCentralityMin()  const {  return fCentralityCutMin; }
  Float_t  GetCentralityMax()  const {  return fCentralityCutMax; }
  //Float_t  GetQVectorPosCutMin()  const {  return fQVectorPosCutMin; }
  //Float_t  GetQVectorPosCutMax()  const {  return fQVectorPosCutMax; }
  //Float_t  GetQVectorNegCutMin()  const {  return fQVectorNegCutMin; }
  //Float_t  GetQVectorNegCutMax()  const {  return fQVectorNegCutMax; }
  Float_t  GetQVectorCutMin()  const {  return fQVectorCutMin; }
  Float_t  GetQVectorCutMax()  const {  return fQVectorCutMax; }
  Float_t  GetVertexCutMin()  const {  return fVertexCutMin; }
  Float_t  GetVertexCutMax()  const {  return fVertexCutMax; }
  Float_t  GetMultiplicityCutMin()  const {  return fMultiplicityCutMin; }
  Float_t  GetMultiplicityCutMax()  const {  return fMultiplicityCutMax; }
  void   PrintCuts();
  Double_t ApplyCentralityPatchAOD049();

  Float_t  NumberOfEvents()     { return fHistoCuts->GetBinContent(kAcceptedEvents+1); }
  Float_t  NumberOfProcessedEvents()     { return fHistoCuts->GetBinContent(kProcessedEvents+1); }
  Float_t  NumberOfPhysSelEvents()     { return fHistoCuts->GetBinContent(kPhysSelEvents+1); }

  Long64_t Merge(TCollection* list);


 private:
  
  AliAODEvent     *fAOD;              //! AOD event
  UInt_t           fTrackBits;       // Type of track to be used in the Qvector calculation
  Bool_t          fIsMC;// true if processing MC
  Bool_t          fCentFromV0;// default centrality with tracks
  Bool_t          fUseCentPatchAOD049;// Patch for centrality selection on AOD049
  AliSpectraAODTrackCuts     *fTrackCuts;              //! track cuts
  Bool_t          fIsSelected;        // True if cuts are selected
  Float_t         fCentralityCutMin;     // minimum centrality percentile
  Float_t         fCentralityCutMax;     // maximum centrality percentile
  /* Float_t         fQVectorPosCutMin;     // minimum qvecPos */
  /* Float_t         fQVectorPosCutMax;     // maximum qvecPos */
  /* Float_t         fQVectorNegCutMin;     // minimum qvecNeg */
  /* Float_t         fQVectorNegCutMax;     // maximum qvecNeg */
  Float_t         fQVectorCutMin;     // minimum qvecPos
  Float_t         fQVectorCutMax;     // maximum qvecPos
  Float_t         fVertexCutMin;     // minimum vertex position
  Float_t         fVertexCutMax;     // maximum vertex position
  Float_t         fMultiplicityCutMin;     // minimum multiplicity position
  Float_t         fMultiplicityCutMax;     // maximum multiplicity position
  TH1I            *fHistoCuts;        // Cuts statistics
  TH1F            *fHistoVtxBefSel;        // Vtx distr before event selection
  TH1F            *fHistoVtxAftSel;        // Vtx distr after event selection
  TH1F            *fHistoEtaBefSel;        // Eta distr before event selection
  TH1F            *fHistoEtaAftSel;        // Eta distr after event selection
  TH1F            *fHistoNChAftSel;        // NCh distr after event selection
  //TH1F            *fHistoQVectorPos;        // QVectorPos
  //TH1F            *fHistoQVectorNeg;        // QVectorNeg
  TH1F            *fHistoQVector;        // QVector
  TH1F            *fHistoEP;        // EP
  AliSpectraAODEventCuts(const AliSpectraAODEventCuts&);
  AliSpectraAODEventCuts& operator=(const AliSpectraAODEventCuts&);
  
  ClassDef(AliSpectraAODEventCuts, 2);
  
};
#endif

