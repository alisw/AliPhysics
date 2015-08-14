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
class TProfile;

#include "TNamed.h"
#include "TFile.h"
#include "TKey.h"
#include "AliOADBContainer.h"
#include "AliVEvent.h"

class AliSpectraAODEventCuts : public TNamed
{
public:
  enum {  kProcessedEvents = 0,kPhysSelEvents,kAcceptedEvents, kVtxRange, kVtxCentral, kVtxNoEvent, kQVector, kNVtxCuts};

  // Constructors
  AliSpectraAODEventCuts() :
    TNamed(),
    fAOD(0),
    fSelectBit(AliVEvent::kMB),
    fCentralityMethod(),
    fTrackBits(1),
    fIsMC(0),
    fIsLHC10h(0),
    fTrackCuts(0),
    fIsSelected(0),
    fCentralityCutMin(0.),
    fCentralityCutMax(999),
    fQVectorCutMin(-999.),
    fQVectorCutMax(999.),
    fVertexCutMin(-10.),
    fVertexCutMax(10.),
    fMultiplicityCutMin(-999.),
    fMultiplicityCutMax(99999.),
    fRejectionFractionTPC(-999.),
    fEtaTPCmin(-0.4),
    fEtaTPCmax(0.4),
    fqTPC(-999.),
    fqV0C(-999.),
    fqV0A(-999.),
    fqV0Cx(-999.),
    fqV0Ax(-999.),
    fqV0Cy(-999.),
    fqV0Ay(-999.),
    fPsiV0C(-999.),
    fPsiV0A(-999.),
    fPsiTPC(-999.),
    fCent(-999.),
    fOutput(0),
    fCalib(0),
    fRun(-1),
    fMultV0(0),
    fV0Cpol1(-1),
    fV0Cpol2(-1),
    fV0Cpol3(-1),
    fV0Cpol4(-1),
    fV0Apol1(-1),
    fV0Apol2(-1),
    fV0Apol3(-1),
    fV0Apol4(-1),
    fQvecIntList(0),
    fQvecIntegral(0), 
    fSplineArrayV0A(0),
    fSplineArrayV0C(0),
    fSplineArrayTPC(0),
    fQgenIntegral(0), 
    fSplineArrayV0Agen(0),
    fSplineArrayV0Cgen(0),
    fSplineArrayTPCgen(0),
    fQvecMC(0),
    fQtrkbit(128),
    fNch(0),
    fQvecCalibType(0),
    fV0Aeff(0)
  {
    for (Int_t i = 0; i<10; i++){
      fMeanQxa2[i] = -1;
      fMeanQya2[i] = -1;
      fMeanQxc2[i] = -1;
      fMeanQyc2[i] = -1;
    }
  }
  AliSpectraAODEventCuts(const char *name);
  virtual  ~AliSpectraAODEventCuts() {}

  void  SetEventSelectionBit( UInt_t val )        { fSelectBit = val;  }
  void  SetCentralityMethod(const char* method) { fCentralityMethod = method; }
  void  SetTrackBits(UInt_t TrackBits) {fTrackBits=TrackBits;}
  void  SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
  void  SetIsLHC10h(Bool_t isLHC10h = kFALSE)    {fIsLHC10h = isLHC10h; };
  void  SetCentralityCutMin(Float_t cut)  { fCentralityCutMin = cut; }
  void  SetCentralityCutMax(Float_t cut)  { fCentralityCutMax = cut; }
  void  SetQVectorCut(Float_t min,Float_t max)  { fQVectorCutMin = min; fQVectorCutMax = max; }
  void  SetVertexCut(Float_t min,Float_t max)  { fVertexCutMin = min; fVertexCutMax = max; }
  void  SetMultiplicityCut(Float_t min,Float_t max)  { fMultiplicityCutMin = min; fMultiplicityCutMax = max; }
  void  SetRejectionFractionTPC(Float_t rej)  { fRejectionFractionTPC = rej; }
  void  SetEtaTPCmin(Float_t min)  { fEtaTPCmin = min; }
  void  SetEtaTPCmax(Float_t max)  { fEtaTPCmax = max; }
  void  SetQtrkbit (UInt_t val)  { fQtrkbit = val; }


  UInt_t GetEventSelectionBit()   const           { return fSelectBit;}
  TString GetCentralityMethod()   const           { return fCentralityMethod;}
  UInt_t GetTrackType()          const    { return fTrackBits;}
  Bool_t GetIsMC()                const           { return fIsMC;}
  Bool_t GetIsLHC10h()                const           { return fIsLHC10h;}
  Float_t  GetCentralityMin()      const {  return fCentralityCutMin; }
  Float_t  GetCentralityMax()     const {  return fCentralityCutMax; }
  Float_t  GetQVectorCutMin()    const {  return fQVectorCutMin; }
  Float_t  GetQVectorCutMax()   const {  return fQVectorCutMax; }
  Float_t  GetVertexCutMin()     const {  return fVertexCutMin; }
  Float_t  GetVertexCutMax()    const {  return fVertexCutMax; }
  Float_t  GetMultiplicityCutMin()  const {  return fMultiplicityCutMin; }
  Float_t  GetMultiplicityCutMax()  const {  return fMultiplicityCutMax; }
  Double_t  GetqTPC()  const {  return fqTPC; }
  Double_t  GetqV0C()  const {  return fqV0C; }
  Double_t  GetqV0A()  const {  return fqV0A; }
  Double_t  GetqV0Cx()  const {  return fqV0Cx; }
  Double_t  GetqV0Ax()  const {  return fqV0Ax; }
  Double_t  GetqV0Cy()  const {  return fqV0Cy; }
  Double_t  GetqV0Ay()  const {  return fqV0Ay; }
  Double_t  GetPsiV0C()  const {  return fPsiV0C; }
  Double_t  GetPsiV0A()  const {  return fPsiV0A; }
  Double_t  GetPsiTPC()  const {  return fPsiTPC; }
  Double_t  GetCent()  const {  return fCent; }
  TList *GetOutputList()       {return fOutput;};
  TList *GetCalibList()       {return fCalib;};
  void SetCalibFile(TFile *f)    {
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
      AliOADBContainer * obj=(AliOADBContainer*)key->ReadObj();
      fCalib->Add(obj);
    }
  };

  void SetQvecIntegralFile(TFile *f)    {
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
      TObject * h=(TObject*)key->ReadObj();
      fQvecIntList->Add(h);
    }
  };

  // Methods
  Bool_t IsSelected(AliAODEvent * aod,AliSpectraAODTrackCuts     *trackcuts);
  Bool_t CheckVtxRange();
  Bool_t CheckCentralityCut();
  Bool_t CheckMultiplicityCut();
  Bool_t CheckQVectorCut();
  Double_t CalculateQVectorLHC10h();  //procedure to calculate the q vector for PbPb 2010
  Double_t CalculateQVectorTPC();  //procedure to calculate the q vector from TPC
  void   PrintCuts();
  Bool_t OpenInfoCalbration(Int_t run);
  Short_t  GetCentrCode(AliVEvent* ev);

  Double_t CalculateQVector(); //q vector calculation using Event plane task

  Float_t  NumberOfEvents()     { return ((TH1I*)fOutput->FindObject("fHistoCuts"))->GetBinContent(kAcceptedEvents+1); }
  Float_t  NumberOfProcessedEvents()     { return ((TH1I*)fOutput->FindObject("fHistoCuts"))->GetBinContent(kProcessedEvents+1); }
  Float_t  NumberOfPhysSelEvents()     { return ((TH1I*)fOutput->FindObject("fHistoCuts"))->GetBinContent(kPhysSelEvents+1); }

  Long64_t Merge(TCollection* list);

  Double_t GetQvecPercentile(Int_t v0side);  
  Bool_t CheckSplineArray(TObjArray * splarr, Int_t n);
  TObjArray *GetSplineArrayV0A() { return fSplineArrayV0A; }
  TObjArray *GetSplineArrayV0C() { return fSplineArrayV0C; }

  Double_t GetQvecMC() {return fQvecMC;}

  Int_t GetNch() { return fNch; }

  void SetQVecCalibType(Int_t val) { fQvecCalibType=val; }  //0. centrality - 1. Nch
  Int_t GetNchBin(TH2D * h);

  Double_t CalculateQVectorMC(Int_t v0side, Int_t type);
  Double_t GetQvecPercentileMC(Int_t v0side, Int_t type);

  Int_t CheckVZEROchannel(Int_t vzeroside, Double_t eta, Double_t phi);
  Int_t CheckVZEROacceptance(Double_t eta);

private:

  AliAODEvent     *fAOD;              //! AOD event
  UInt_t           fSelectBit;            // Select events according to AliAnalysisTaskJetServices bit maps 
  TString          fCentralityMethod;     // Method to determine centrality
  UInt_t           fTrackBits;       // Type of track to be used in the Qvector calculation
  Bool_t          fIsMC; // true if processing MC
  Bool_t          fIsLHC10h; // for LHC10h
  AliSpectraAODTrackCuts     *fTrackCuts;              //! track cuts
  Bool_t          fIsSelected;        // True if cuts are selected
  Float_t         fCentralityCutMin;     // minimum centrality percentile
  Float_t         fCentralityCutMax;     // maximum centrality percentile
  Float_t         fQVectorCutMin;     // minimum qvecPos
  Float_t         fQVectorCutMax;     // maximum qvecPos
  Float_t         fVertexCutMin;     // minimum vertex position
  Float_t         fVertexCutMax;     // maximum vertex position
  Float_t         fMultiplicityCutMin;     // minimum multiplicity position
  Float_t         fMultiplicityCutMax;     // maximum multiplicity position
  Double_t       fRejectionFractionTPC;            //rejection fraction in TPC calculation, (e.g. 0.1 means that 10 % of the tracks will be rejected in the calculation)
  Double_t       fEtaTPCmin;            //q vector in the TPC
  Double_t       fEtaTPCmax;            //q vector in the TPC
  Double_t       fqTPC;            //q vector in the TPC
  Double_t       fqV0C;            //q vector in the VZERO-C
  Double_t       fqV0A;            //q vector in the VZERO-A
  Double_t       fqV0Cx;            //q vector in the VZERO-C, x-comp
  Double_t       fqV0Ax;            //q vector in the VZERO-A, x-comp
  Double_t       fqV0Cy;            //q vector in the VZERO-C, y-comp
  Double_t       fqV0Ay;            //q vector in the VZERO-A, y-comp
  Double_t       fPsiV0C;            //EP from VZERO-C
  Double_t       fPsiV0A;            //EP from VZERO-A
  Double_t       fPsiTPC;            //EP from TPC
  Double_t       fCent;            //centrality according to fCentralityMethod
  TList            *fOutput;        // output list 
  TList            *fCalib;        // output list 
  Int_t fRun;                       // run number - for calibration
  TProfile* fMultV0;                //! profile from V0 multiplicity
  Float_t fV0Cpol1;                 // mean V0C multiplicity - from fit profile - ring 1
  Float_t fV0Cpol2;                 // mean V0C multiplicity - from fit profile - ring 2
  Float_t fV0Cpol3;                 // mean V0C multiplicity - from fit profile - ring 3
  Float_t fV0Cpol4;                 // mean V0C multiplicity - from fit profile - ring 4
  Float_t fV0Apol1;                 // mean V0A multiplicity - from fit profile - ring 1
  Float_t fV0Apol2;                 // mean V0A multiplicity - from fit profile - ring 2
  Float_t fV0Apol3;                 // mean V0A multiplicity - from fit profile - ring 3
  Float_t fV0Apol4;                 // mean V0A multiplicity - from fit profile - ring 4
  Float_t fMeanQxa2[10];             // mean Qxa values for centr - recentering
  Float_t fMeanQya2[10];             // mean Qya values for centr - recentering
  Float_t fMeanQxc2[10];             // mean Qxc values for centr - recentering
  Float_t fMeanQyc2[10];             // mean Qyc values for centr - recentering

  TList *fQvecIntList;            // List with Qvec Integrated vs centrality distribution
  TH2D * fQvecIntegral;           // ! Integrated Qvec distribution
  TObjArray * fSplineArrayV0A;    // TSpline array for VZERO-A
  TObjArray * fSplineArrayV0C;    // TSpline array for VZERO-C
  TObjArray * fSplineArrayTPC;    // TSpline array for TPC
  TH2D * fQgenIntegral;           // ! Integrated Qvec distribution for generated tracks
  TObjArray * fSplineArrayV0Agen;    // TSpline array for VZERO-A for generated tracks
  TObjArray * fSplineArrayV0Cgen;    // TSpline array for VZERO-C for generated tracks
  TObjArray * fSplineArrayTPCgen;    // TSpline array for TPC for generated tracks
  Double_t fQvecMC; //q-vector value from MC
  UInt_t fQtrkbit; // filterbit for q2_TPC distribution

  Int_t fNch;
  Int_t fQvecCalibType; //0. centrality - 1. Nch
  TH1F * fV0Aeff; // VZEROA efficiency prim+sec / gen.

  AliSpectraAODEventCuts(const AliSpectraAODEventCuts&);
  AliSpectraAODEventCuts& operator=(const AliSpectraAODEventCuts&);

  ClassDef(AliSpectraAODEventCuts, 14);

};
#endif

