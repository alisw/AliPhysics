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

#include <AliAnalysisFilter.h>
#include <AliKFParticle.h>

#include "AliDielectronHistos.h"
#include "AliDielectronPair.h"

class AliVEvent;
class THashList;
class AliDielectronCF;
class AliDielectronDebugTree;

//________________________________________________________________
class AliDielectron : public TNamed {
  
public:
  enum ParticleValues { kPx=0, kPy, kPz, kPt, kP, kXv, kYv, kZv, kOneOverPt,
      kPhi, kTheta, kE, kM, kEta, kY, kCharge, kNParticleValues };
  enum PairValues { kChi2NDF=kNParticleValues, kDecayLength, kR, kOpeningAngle, kMerr, kNPairValues };
  
  AliDielectron();
  AliDielectron(const char* name, const char* title);
  virtual ~AliDielectron();

  void Init();
  
  void Process(AliVEvent *ev1, AliVEvent *ev2=0);

  const AliAnalysisFilter& GetEventFilter() const { return fEventFilter; }
  const AliAnalysisFilter& GetTrackFilter() const { return fTrackFilter; }
  const AliAnalysisFilter& GetPairFilter()  const { return fPairFilter;  }

  AliAnalysisFilter& GetEventFilter() { return fEventFilter; }
  AliAnalysisFilter& GetTrackFilter() { return fTrackFilter; }
  AliAnalysisFilter& GetPairFilter()  { return fPairFilter;  }

  void  SetMotherPdg( Int_t pdgMother ) { fPdgMother=pdgMother; }
  void  SetLegPdg(Int_t pdgLeg1, Int_t pdgLeg2) { fPdgLeg1=pdgLeg1; fPdgLeg2=pdgLeg2; }
  Int_t GetMotherPdg() const { return fPdgMother; }
  Int_t GetLeg1Pdg()   const { return fPdgLeg1;   }
  Int_t GetLeg2Pdg()   const { return fPdgLeg2;   }
  
  const TObjArray* GetTrackArray(Int_t i) const {return (i>=0&&i<4)?&fTracks[i]:0;}
  const TObjArray* GetPairArray(Int_t i)  const {return (i>=0&&i<10)?
      static_cast<TObjArray*>(fPairCandidates->UncheckedAt(i)):0;}

  TObjArray** GetPairArraysPointer() { return &fPairCandidates; }
  
  void SetHistogramManager(AliDielectronHistos * const histos) { fHistos=histos; }
  const THashList * GetHistogramList() const { return fHistos?fHistos->GetHistogramList():0x0; }

  Bool_t HasCandidates() const { return GetPairArray(1)?GetPairArray(1)->GetEntriesFast()>0:0; }

  void SetCFManagerPair(AliDielectronCF * const cf) { fCfManagerPair=cf; }
  AliDielectronCF* GetCFManagerPair() const { return fCfManagerPair; }

  void SetDebugTree(AliDielectronDebugTree * const tree) { fDebugTree=tree; }
  
  static const char* TrackClassName(Int_t i) { return (i>=0&&i<4)?fgkTrackClassNames[i]:""; }
  static const char* PairClassName(Int_t i)  { return (i>=0&&i<10)?fgkPairClassNames[i]:""; }

  void SaveDebugTree();
  
private:

  
  AliAnalysisFilter fEventFilter; // Event cuts
  AliAnalysisFilter fTrackFilter; // leg cuts
  AliAnalysisFilter fPairFilter;  // pair cuts
  
  Int_t fPdgMother;     // pdg code of mother tracks
  Int_t fPdgLeg1;       // pdg code leg1
  Int_t fPdgLeg2;       // pdg code leg2
  
  AliDielectronHistos *fHistos;   // Histogram manager
                                  //  Streaming and merging should be handled
                                  //  by the analysis framework
  
  TObjArray fTracks[4];           //! Selected track candidates
                                  //  0: Event1, positive particles
                                  //  1: Event1, negative particles
                                  //  2: Event2, positive particles
                                  //  3: Event2, negative particles

  TObjArray *fPairCandidates;     //! Pair candidate arrays
                                  //TODO: better way to store it? TClonesArray?

  AliDielectronCF *fCfManagerPair;//Correction Framework Manager for the Pair

  AliDielectronDebugTree *fDebugTree;  // Debug tree output
  
  void FillTrackArrays(AliVEvent * const ev, Int_t eventNr=0);
  void FillPairArrays(Int_t arr1, Int_t arr2);
  
  Int_t GetPairIndex(Int_t arr1, Int_t arr2) const {return arr1>=arr2?arr1*(arr1+1)/2+arr2:arr2*(arr2+1)/2+arr1;}

  void InitPairCandidateArrays();
  void ClearArrays();
  
  TObjArray* PairArray(Int_t i);
  
  static const char* fgkTrackClassNames[4];   //Names for track arrays
  static const char* fgkPairClassNames[10];   //Names for pair arrays

  void ProcessMC();
  
  void  FillHistograms(const AliVEvent *ev);
  void  FillDebugTree();
  
  AliDielectron(const AliDielectron &c);
  AliDielectron &operator=(const AliDielectron &c);
  
  ClassDef(AliDielectron,3);
};

inline void AliDielectron::InitPairCandidateArrays()
{
  //
  // initialise all pair candidate arrays
  //
  fPairCandidates->SetOwner();
  for (Int_t i=0;i<10;++i){
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
  for (Int_t i=0;i<10;++i){
    PairArray(i)->Delete();
  }
}

#endif
