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
    fqV0C(-999.),
    fqV0A(-999.),
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
    fV0Apol4(-1)
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
  Double_t  GetqV0C()  const {  return fqV0C; }
  Double_t  GetqV0A()  const {  return fqV0A; }
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
  
  // Methods
  Bool_t IsSelected(AliAODEvent * aod,AliSpectraAODTrackCuts     *trackcuts);
  Bool_t CheckVtxRange();
  Bool_t CheckCentralityCut();
  Bool_t CheckMultiplicityCut();
  Bool_t CheckQVectorCut();
  Double_t CalculateQVectorLHC10h();
  void   PrintCuts();
  Bool_t OpenInfoCalbration(Int_t run);
  Short_t  GetCentrCode(AliVEvent* ev);
  
  Float_t  NumberOfEvents()     { return ((TH1I*)fOutput->FindObject("fHistoCuts"))->GetBinContent(kAcceptedEvents+1); }
  Float_t  NumberOfProcessedEvents()     { return ((TH1I*)fOutput->FindObject("fHistoCuts"))->GetBinContent(kProcessedEvents+1); }
  Float_t  NumberOfPhysSelEvents()     { return ((TH1I*)fOutput->FindObject("fHistoCuts"))->GetBinContent(kPhysSelEvents+1); }

  Long64_t Merge(TCollection* list);
  

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
  Double_t       fqV0C;            //q vector in the VZERO-C
  Double_t       fqV0A;            //q vector in the VZERO-A
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

  AliSpectraAODEventCuts(const AliSpectraAODEventCuts&);
  AliSpectraAODEventCuts& operator=(const AliSpectraAODEventCuts&);
  
  ClassDef(AliSpectraAODEventCuts, 3);
  
};
#endif

