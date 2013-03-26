#ifndef ALISPECTRABOTHTRACKCUTS_H
#define ALISPECTRABOTHTRACKCUTS_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliSpectraBothTrackCuts
//
//
//
//
// Authors: Michele Floris, CERN, Philip Versteeg, UU, Redmer Bertens, UU
//-------------------------------------------------------------------------

class AliAODEvent;
//class AliSpectraBothHistoManager;
class TH1I;
class AliAODMCParticle;
class AliAODTrack;
class  AliESDtrackCuts;
#include "AliSpectraBothHistoManager.h"
#include "TNamed.h"
#include "AliESDtrackCuts.h"
 
using namespace AliSpectraNameSpaceBoth;

class AliSpectraBothTrackCuts : public TNamed
{
 public:
  
  enum { kTrkBit = 0, kTrkCuts, kTrkEta, kTrkDCA, kTrkP, kTrkPt,kTrkPtTOF,kTOFMatching,kTrTOFout,kTrTIME,kTrTOFpid,kAccepted,kNTrkCuts};
  enum {kAODobject=0,kESDobject,kotherobject};
  
 AliSpectraBothTrackCuts() : TNamed(), fIsSelected(0), fTrackBits(0), fMinTPCcls(0), fEtaCutMin(0), fEtaCutMax(0),fDCACut(0) ,fPCut(0), fPtCut(0),fYCutMax(0),fYCutMin(0), fPtCutTOFMatching(0),fAODtrack(0), fHashitinSPD1(0),fusedadditionalcuts(kTRUE),
fHistoCuts(0), fHistoNSelectedPos(0), fHistoNSelectedNeg(0), fHistoNMatchedPos(0), fHistoNMatchedNeg(0), fHistoEtaPhiHighPt(0), fHistoNclustersITS(0),fTrack(0),fCuts(0) {}
  
  AliSpectraBothTrackCuts(const char *name);
  virtual  ~AliSpectraBothTrackCuts() {} // To be implemented
  
  Bool_t IsSelected(AliVTrack * track,Bool_t FillHistStat);
  
  void SetTrackType(UInt_t bit);
  Bool_t CheckTrackType();
  Bool_t CheckTrackCuts();
  Bool_t CheckEtaCut();
  Bool_t CheckYCut(BothParticleSpecies_t specie); // not included in standard cuts
  Bool_t CheckDCACut();
  Bool_t CheckPCut();
  Bool_t CheckPtCut();
  Bool_t CheckTOFMatching(Bool_t FillHistStat);
  void PrintCuts() const;
  
   UInt_t GetTrackType()  const    { return fTrackBits;}
   TH1I * GetHistoCuts()      { return fHistoCuts; }
   TH1F * GetHistoNSelectedPos()      { return fHistoNSelectedPos; } 
   TH1F * GetHistoNSelectedNeg()      { return fHistoNSelectedNeg; }
   TH1F * GetHistoNMatchedPos()      { return fHistoNMatchedPos; }
   TH1F * GetHistoNMatchedNeg()      { return fHistoNMatchedNeg; }
   TH2F * GetHistoEtaPhiHighPt()      { return fHistoEtaPhiHighPt; }
   TH1F * GetHistoNclustersITS()  {return fHistoNclustersITS;}  
   void SetEta(Float_t etamin,Float_t etamax)   { fEtaCutMin = etamin;fEtaCutMax = etamax; }
   void SetDCA(Float_t dca)   { fDCACut = dca; }
   void SetP(Float_t p)       { fPCut = p; }
   void SetPt(Float_t pt)     { fPtCut = pt; }
   void SetY(Float_t ymax,Float_t ymin) { fYCutMax = ymax;fYCutMin=ymin;}
   void SetPtTOFMatching(Float_t pt)     { fPtCutTOFMatching = pt; }
   void SetTrackBits(UInt_t TrackBits) {fTrackBits=TrackBits;}
   void SetMinTPCcls(UInt_t MinTPCcls) {fMinTPCcls=MinTPCcls;}
   void SetHashitinSPD1 (Bool_t value) {fHashitinSPD1=value;}	
   void SetUsedAdditionalCuts (Bool_t value) {fusedadditionalcuts=value;}	
   Float_t GetEtaMin()       const    { return fEtaCutMin; }
   Float_t GetEtaMax()       const    { return fEtaCutMax; }
   Float_t GetYMax()         const    { return fYCutMax; }
   Float_t GetYMin()         const    { return fYCutMin; }
   Float_t GetY()         const    { return 0.5*(fYCutMax-fYCutMin); }
   Float_t GetDCA()       const    { return fDCACut; }
   Float_t GetP()         const    { return fPCut; }
   Float_t GetPt()        const    { return fPtCut; }
   Float_t GetPtTOFMatching()        const    { return fPtCutTOFMatching; } 
   Long64_t Merge(TCollection* list);
   void SetAliESDtrackCuts(AliESDtrackCuts*  cuts ){fCuts=cuts;}
   
 private:
   
   Bool_t           fIsSelected;      // True if cuts are selected
   UInt_t           fTrackBits;       // Type of track to be used
   UInt_t           fMinTPCcls;       // min number of clusters in the TPC
   Float_t          fEtaCutMin;          // Allowed absolute maximum value of Eta
   Float_t          fEtaCutMax;          // Allowed absolute maximum value of Eta
   Float_t          fDCACut;          // Maximum value of DCA
   Float_t          fPCut;            // Maximum value of P
   Float_t          fPtCut;           // Maximum value of Pt
   Float_t          fYCutMax;           // Maximum value of Y 
   Float_t          fYCutMin;           // Minimum value of Y
   Float_t          fPtCutTOFMatching;           // TOF Matching
   Int_t 	     fAODtrack; // 0 ESD track connected , 1 AOD track conected , else nothing
   Bool_t 	    fHashitinSPD1; // Check if SPD1 has a hit 	
   Bool_t          fusedadditionalcuts;          //If set to true the TPCrefit, ITSrefit, SPDany and Ncluster cut is check 		
   TH1I             *fHistoCuts;       // Cuts statistics
   TH1F             *fHistoNSelectedPos;       // Selected positive tracks
   TH1F             *fHistoNSelectedNeg;       // Selected negative tracks
   TH1F             *fHistoNMatchedPos;       // Matched positive tracks
   TH1F             *fHistoNMatchedNeg;       // Matched negative tracks
   TH2F             *fHistoEtaPhiHighPt;       // EtaPhi distr at high pt (>1.5 GeV/c)
   TH1F            *fHistoNclustersITS;      // Number of clusters in ITS 
       	 
 		
   AliVTrack      *fTrack;           //! Track pointer
   AliESDtrackCuts *fCuts;      //! cuts  
   static const char * kBinLabel[]; // labels of stat histo

   
   AliSpectraBothTrackCuts(const AliSpectraBothTrackCuts&);
   AliSpectraBothTrackCuts& operator=(const AliSpectraBothTrackCuts&);
   
   ClassDef(AliSpectraBothTrackCuts, 5);
};
#endif

