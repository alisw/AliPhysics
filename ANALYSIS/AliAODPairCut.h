#ifndef ALIAODPAIRCUT_H
#define ALIAODPAIRCUT_H

/* $Id$ */

//Piotr Skowronski@cern.ch
//Class implements cut on the pair of particles
//
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
#include <TNamed.h> 
#include "AliAODPairBaseCut.h"

class AliAODParticleCut;
class AliAODPairBaseCut;

/******************************************************************/

class AliAODPairCut: public TNamed
{
 public:
  AliAODPairCut();
  AliAODPairCut(const AliAODPairCut& in);
  AliAODPairCut& operator = (const AliAODPairCut& in);
  
  virtual ~AliAODPairCut();
  virtual Bool_t Rejected(AliAODPair* pair) const;
  virtual Bool_t PassPairProp(AliAODPair* pair) const;
     
  virtual Bool_t IsEmpty() const {return kFALSE;}
  void SetFirstPartCut(AliAODParticleCut* cut);  //sets the cut on the first particle
  void SetSecondPartCut(AliAODParticleCut* cut); //sets the cut on the second particle
  
  void SetPartCut(AliAODParticleCut* cut);//sets the the same cut on both particles
  
  virtual void AddBasePairCut(AliAODPairBaseCut* cut);
  
  virtual void Print();

  void SetDeltaERange(Double_t min, Double_t max);
  void SetDeltaPRange(Double_t min, Double_t max);
  
  void SetQInvRange(Double_t min, Double_t max);
  void SetKtRange(Double_t min, Double_t max);
  void SetKStarRange(Double_t min, Double_t max);
  void SetKStarOutRange(Double_t min, Double_t max);
  void SetKStarSideRange(Double_t min, Double_t max);
  void SetKStarLongRange(Double_t min, Double_t max);
  void SetQOutLCMSRange(Double_t min, Double_t max);
  void SetQSideLCMSRange(Double_t min, Double_t max);
  void SetQLongLCMSRange(Double_t min, Double_t max);
  void SetAvSeparationRange(Double_t min,Double_t max = 10e5);//Anti-Merging Cut
  void SetITSSeparation(Int_t layer, Double_t drphi=0.01,Double_t dz = 0.08);//Anti-Merging Cut for first pixel layer
  void SetClusterOverlapRange(Double_t min,Double_t max);//Anti-Splitting Max range -0.5 1.0
      
  AliAODParticleCut* GetFirstPartCut() const {return fFirstPartCut;}
  AliAODParticleCut* GetSecondPartCut() const {return fSecondPartCut;}
  
 protected:
  AliAODParticleCut*      fFirstPartCut;//cut on first particle in pair
  AliAODParticleCut*      fSecondPartCut;//cut on second particle in pair
  
  AliAODPairBaseCut** fCuts; //! array of poiters to base cuts
  Int_t fNCuts;//Number of cuts in fCuts array
  
  
  AliAODPairBaseCut* FindCut(AliAODPairBaseCut::EAODPairCutProperty cut);
 private:
  static const Int_t fgkMaxCuts; // Max number of cuts
  ClassDef(AliAODPairCut,2)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliAODPairEmptyCut:  public AliAODPairCut
{
  //Empty - it passes possitively all particles - it means returns always False
  //Class describing cut on pairs of particles
 public:
  AliAODPairEmptyCut(){};
  AliAODPairEmptyCut(const AliAODPairEmptyCut& in):AliAODPairCut(in){};
  virtual ~AliAODPairEmptyCut(){};
  
  Bool_t Rejected(AliAODPair*) const {return kFALSE;} //accpept everything
  Bool_t IsEmpty() const {return kTRUE;}
  
  ClassDef(AliAODPairEmptyCut,1)
};


#endif
