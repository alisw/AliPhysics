#ifndef ALIHBTPAIRCUT_H
#define ALIHBTPAIRCUT_H

/* $Id$ */

//Piotr Skowronski@cern.ch
//Class implements cut on the pair of particles
//
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
 
#include "AliHBTPair.h"

class AliHBTParticleCut;
class AliHbtBasePairCut;

enum AliHBTPairCutProperty
{
  kHbtPairCutPropQInv, //Q invariant
  kHbtPairCutPropKt,
  kHbtPairCutPropKStar,
  kHbtPairCutPropQSideLCMS,
  kHbtPairCutPropQOutLCMS,
  kHbtPairCutPropQLongLCMS,
  kHbtPairCutPropDeltaPhi,
  kHbtPairCutPropDeltaTheta,
  kHbtPairCutPropDeltaP,
  kHbtPairCutPropDeltaPt,
  kHbtPairCutPropAvSepar,
  kHbtPairCutPropSepar,
  kHbtPairCutPropClOverlap,
  kHbtPairCutPropPixelSepar,
  kHbtPairCutPropNone
};
/******************************************************************/

class AliHBTPairCut: public TNamed
{
 public:
  AliHBTPairCut();
  AliHBTPairCut(const AliHBTPairCut& in);
  AliHBTPairCut& operator = (const AliHBTPairCut& in);
  
  virtual ~AliHBTPairCut();
  virtual Bool_t Pass(AliHBTPair* pair) const;
  virtual Bool_t PassPairProp(AliHBTPair* pair) const;
     
  virtual Bool_t IsEmpty() const {return kFALSE;}
  void SetFirstPartCut(AliHBTParticleCut* cut);  //sets the cut on the first particle
  void SetSecondPartCut(AliHBTParticleCut* cut); //sets the cut on the second particle
  
  void SetPartCut(AliHBTParticleCut* cut);//sets the the same cut on both particles
  
  virtual void AddBasePairCut(AliHbtBasePairCut* cut);
  
  virtual void Print();
  
  void SetQInvRange(Double_t min, Double_t max);
  void SetKtRange(Double_t min, Double_t max);
  void SetKStarRange(Double_t min, Double_t max);
  void SetQOutCMSLRange(Double_t min, Double_t max);
  void SetQSideCMSLRange(Double_t min, Double_t max);
  void SetQLongCMSLRange(Double_t min, Double_t max);
  void SetAvSeparationRange(Double_t min,Double_t max = 10e5);//Anti-Merging Cut
  void SetITSSeparation(Int_t layer, Double_t drphi=0.01,Double_t dz = 0.08);//Anti-Merging Cut for first pixel layer
  void SetClusterOverlapRange(Double_t min,Double_t max);//Anti-Splitting Max range -0.5 1.0
      
  AliHBTParticleCut* GetFirstPartCut() const {return fFirstPartCut;}
  AliHBTParticleCut* GetSecondPartCut() const {return fSecondPartCut;}
  
 protected:
  AliHBTParticleCut*      fFirstPartCut;//cut on first particle in pair
  AliHBTParticleCut*      fSecondPartCut;//cut on second particle in pair
  
  AliHbtBasePairCut** fCuts; //! array of poiters to base cuts
  Int_t fNCuts;//Number of cuts in fCuts array
  
  
  AliHbtBasePairCut* FindCut(AliHBTPairCutProperty cut);
 private:
  static const Int_t fgkMaxCuts; // Max number of cuts
  ClassDef(AliHBTPairCut,3)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/
#endif
