#ifndef AliAODPairCUT_H
#define AliAODPairCUT_H

/* $Id$ */

//Piotr Skowronski@cern.ch
//Class implements cut on the pair of particles
//
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
 
#include "AliAODPair.h"

class AliAODParticleCut;
class AliAODBasePairCut;

enum AliAODPairCutProperty
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

class AliAODPairCut: public TNamed
{
 public:
  AliAODPairCut();
  AliAODPairCut(const AliAODPairCut& in);
  AliAODPairCut& operator = (const AliAODPairCut& in);
  
  virtual ~AliAODPairCut();
  virtual Bool_t Pass(AliAODPair* pair) const;
  virtual Bool_t PassPairProp(AliAODPair* pair) const;
     
  virtual Bool_t IsEmpty() const {return kFALSE;}
  void SetFirstPartCut(AliAODParticleCut* cut);  //sets the cut on the first particle
  void SetSecondPartCut(AliAODParticleCut* cut); //sets the cut on the second particle
  
  void SetPartCut(AliAODParticleCut* cut);//sets the the same cut on both particles
  
  virtual void AddBasePairCut(AliAODBasePairCut* cut);
  
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
      
  AliAODParticleCut* GetFirstPartCut() const {return fFirstPartCut;}
  AliAODParticleCut* GetSecondPartCut() const {return fSecondPartCut;}
  
 protected:
  AliAODParticleCut*      fFirstPartCut;//cut on first particle in pair
  AliAODParticleCut*      fSecondPartCut;//cut on second particle in pair
  
  AliAODBasePairCut** fCuts; //! array of poiters to base cuts
  Int_t fNCuts;//Number of cuts in fCuts array
  
  
  AliAODBasePairCut* FindCut(AliAODPairCutProperty cut);
 private:
  static const Int_t fgkMaxCuts; // Max number of cuts
  ClassDef(AliAODPairCut,2)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliAODEmptyPairCut:  public AliAODPairCut
{
  //Empty - it passes possitively all particles - it means returns always False
  //Class describing cut on pairs of particles
 public:
  AliAODEmptyPairCut(){};
  AliAODEmptyPairCut(const AliAODEmptyPairCut& in):AliAODPairCut(in){};
  virtual ~AliAODEmptyPairCut(){};
  
  Bool_t Pass(AliAODPair*) const {return kFALSE;} //accpept everything
  Bool_t IsEmpty() const {return kTRUE;}
  
  ClassDef(AliAODEmptyPairCut,1)
};



/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliAODBasePairCut: public TObject
{
  //This class defines the range of some property - pure virtual
  //Property is coded by AliAODCutTypes type
   
 public:
     
  AliAODBasePairCut(Double_t min = 0.0, Double_t max = 0.0, AliAODPairCutProperty prop= kHbtPairCutPropNone):
    fMin(min),fMax(max),fProperty(prop){}
  
  virtual   ~AliAODBasePairCut(){}
     
  virtual Bool_t    Pass(AliAODPair* pair) const;
  
  void      SetRange(Double_t min, Double_t max){fMin = min; fMax = max;}
  
  void      SetMinimum(Double_t min){fMin = min;}
  void      SetMaximum(Double_t max){fMax = max;}
  
  Double_t  GetMinimum() const {return fMin;}
  Double_t  GetMaximum() const {return fMax;}
  
  AliAODPairCutProperty GetProperty() const {return fProperty;}
  
 protected:
  virtual Double_t  GetValue(AliAODPair* pair) const = 0;
  
  Double_t fMin; // Lower boundary of the range
  Double_t fMax; // Upper boundary of the range
  
  AliAODPairCutProperty fProperty; // The property itself
  
  ClassDef(AliAODBasePairCut,1)
 
 };
/******************************************************************/

inline Bool_t AliAODBasePairCut::Pass(AliAODPair* pair) const
{
  //checks if pair proprty is in range
  //null pointer check is made by AliAODPairCut, so here is unnecesary
  
  Double_t value = GetValue(pair);
  if ( (value > fMin) && (value <fMax ) ) return kFALSE; //accepted
  else return kTRUE; //rejected
}
/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliAODQInvCut: public AliAODBasePairCut
{
 public:
  AliAODQInvCut(Double_t min = 0.0, Double_t max = 0.0):AliAODBasePairCut(min,max,kHbtPairCutPropQInv){}
  virtual ~AliAODQInvCut(){}
 protected:
  virtual Double_t  GetValue(AliAODPair* pair) const {return pair->GetQInv();}
  
  ClassDef(AliAODQInvCut,1)
 };
/******************************************************************/

class AliAODKtCut: public AliAODBasePairCut {
 public:
  AliAODKtCut(Double_t min = 0.0, Double_t max = 0.0):AliAODBasePairCut(min,max,kHbtPairCutPropKt){}
  virtual ~AliAODKtCut(){}
 protected:
  virtual Double_t  GetValue(AliAODPair* pair) const {return pair->GetKt();}

  ClassDef(AliAODKtCut,1)
 };
/******************************************************************/

class AliAODKStarCut: public AliAODBasePairCut
{
 public:
  AliAODKStarCut(Double_t min = 0.0, Double_t max = 0.0):AliAODBasePairCut(min,max,kHbtPairCutPropKStar){}
  virtual ~AliAODKStarCut(){}
 protected:
  virtual Double_t  GetValue(AliAODPair* pair) const {return pair->GetKStar();}

  ClassDef(AliAODKStarCut,1)
};
/******************************************************************/

class AliAODQSideLCMSCut: public AliAODBasePairCut
{
 public:
  AliAODQSideLCMSCut(Double_t min = 0.0, Double_t max = 0.0):
    AliAODBasePairCut(min,max,kHbtPairCutPropQSideLCMS){}
  virtual ~AliAODQSideLCMSCut(){}
 protected:
  virtual Double_t  GetValue(AliAODPair* pair) const 
    {return pair->GetQSideLCMS();}

  ClassDef(AliAODQSideLCMSCut,1)
};
/******************************************************************/


class AliAODQOutLCMSCut: public AliAODBasePairCut
{
 public:
  AliAODQOutLCMSCut(Double_t min = 0.0, Double_t max = 0.0):
    AliAODBasePairCut(min,max,kHbtPairCutPropQOutLCMS){}
  virtual ~AliAODQOutLCMSCut(){}
 protected:
  virtual Double_t  GetValue(AliAODPair* pair) const 
    {return pair->GetQOutLCMS();}
  
  ClassDef(AliAODQOutLCMSCut,1)
};
/******************************************************************/

class AliAODQLongLCMSCut: public AliAODBasePairCut
{
 public:
  AliAODQLongLCMSCut(Double_t min = 0.0, Double_t max = 0.0):
    AliAODBasePairCut(min,max,kHbtPairCutPropQLongLCMS){}
  virtual ~AliAODQLongLCMSCut(){}
 protected:
  virtual Double_t  GetValue(AliAODPair* pair) const 
    {return pair->GetQLongLCMS();}

  ClassDef(AliAODQLongLCMSCut,1)
};
/******************************************************************/

class AliAODDeltaPhiCut: public AliAODBasePairCut
{
 public:
  AliAODDeltaPhiCut(Double_t min = 0.0, Double_t max = 0.0):
    AliAODBasePairCut(min,max,kHbtPairCutPropDeltaPhi){}
  virtual ~AliAODDeltaPhiCut(){}
 protected:
  virtual Double_t  GetValue(AliAODPair* pair) const 
    {return TMath::Abs(pair->GetDeltaPhi());}

  ClassDef(AliAODDeltaPhiCut,1)
};
/******************************************************************/

class AliAODDeltaThetaCut: public AliAODBasePairCut
{
 public:
  AliAODDeltaThetaCut(Double_t min = 0.0, Double_t max = 0.0):
    AliAODBasePairCut(min,max,kHbtPairCutPropDeltaTheta){}
  virtual ~AliAODDeltaThetaCut(){}
 protected:
  virtual Double_t  GetValue(AliAODPair* pair) const 
    {return TMath::Abs(pair->GetDeltaTheta());}

  ClassDef(AliAODDeltaThetaCut,1)
};
/******************************************************************/

class AliAODCluterOverlapCut: public AliAODBasePairCut
{
 public:
  AliAODCluterOverlapCut(Double_t min = 0.0, Double_t max = 1e5):
    AliAODBasePairCut(min,max,kHbtPairCutPropClOverlap){}
  virtual ~AliAODCluterOverlapCut(){}

 protected:
  virtual Double_t  GetValue(AliAODPair* pair) const;
  ClassDef(AliAODCluterOverlapCut,1)
};
/******************************************************************/
  
class AliAODAvSeparationCut: public AliAODBasePairCut
{
 public:
  AliAODAvSeparationCut(Double_t min = 0.0, Double_t max = 1e5):
    AliAODBasePairCut(min,max,kHbtPairCutPropAvSepar){}
  virtual ~AliAODAvSeparationCut(){}
  
 protected:
  virtual Double_t  GetValue(AliAODPair* pair) const;
  ClassDef(AliAODAvSeparationCut,1)
};
/******************************************************************/
  
class AliAODSeparationCut: public AliAODBasePairCut
{
 public:
  AliAODSeparationCut(Double_t min = 0.0, Double_t max = 1e5, Int_t point = 0):
    AliAODBasePairCut(min,max,kHbtPairCutPropSepar),fPoint(point){}
  virtual ~AliAODSeparationCut(){}
  
 protected:
  Int_t fPoint;//index of the point that distance should be measured
  virtual Double_t  GetValue(AliAODPair* pair) const;
  ClassDef(AliAODSeparationCut,1)
};
/******************************************************************/
  
class AliAODITSSeparationCut: public AliAODBasePairCut
{
//Anti merging cut for the first layer of pixels
 public:
  AliAODITSSeparationCut(Int_t layer = 0, Double_t deltarphi = 0.01, Double_t deltaz = 0.08):
    AliAODBasePairCut(deltarphi,deltaz,kHbtPairCutPropPixelSepar),fLayer(layer){}
  virtual ~AliAODITSSeparationCut(){}
  Bool_t   Pass(AliAODPair* pair) const;
  Int_t    GetLayer() const {return fLayer;}
 protected:
  Int_t fLayer;//index of the layer that distance should be measured 0: 1st pixels
  virtual Double_t  GetValue(AliAODPair* /*pair*/) const {return 0.0;}//not used
  ClassDef(AliAODITSSeparationCut,1)
};
/******************************************************************/

class AliAODOutSideSameSignCut: public AliAODBasePairCut
{
 public:
  AliAODOutSideSameSignCut(){}
  virtual ~AliAODOutSideSameSignCut(){}
  virtual Bool_t Pass(AliAODPair *p) const;
 protected:
  virtual Double_t  GetValue(AliAODPair* /*pair*/) const {return 0.0;}
  ClassDef(AliAODOutSideSameSignCut,1)
};
/******************************************************************/

class AliAODOutSideDiffSignCut: public AliAODBasePairCut
{
 public:
  AliAODOutSideDiffSignCut(){}
  virtual ~AliAODOutSideDiffSignCut(){}
  virtual Bool_t Pass(AliAODPair *p) const;
 protected:
  virtual Double_t  GetValue(AliAODPair* /*pair*/) const {return 0.0;}
  ClassDef(AliAODOutSideDiffSignCut,1)
};
/******************************************************************/

class AliAODLogicalOperPairCut:  public AliAODBasePairCut
 {
   public:
     AliAODLogicalOperPairCut();
     AliAODLogicalOperPairCut(AliAODBasePairCut* first, AliAODBasePairCut* second);
     virtual   ~AliAODLogicalOperPairCut();
   protected:
     Double_t  GetValue(AliAODPair * /*pair*/) const {MayNotUse("GetValue");return 0.0;}

     AliAODBasePairCut* fFirst;   //second cut
     AliAODBasePairCut* fSecond;  //first cut
   private:
    class  AliAODDummyBasePairCut: public AliAODBasePairCut
     {
       Double_t  GetValue(AliAODPair* /*pair*/) const {return 0.0;}
       Bool_t    Pass(AliAODPair* /*pair*/) const;
     };

    ClassDef(AliAODLogicalOperPairCut,1)
 };
/******************************************************************/

class AliAODOrPairCut: public AliAODLogicalOperPairCut
{
   public:
     AliAODOrPairCut(){}
     AliAODOrPairCut(AliAODBasePairCut* first, AliAODBasePairCut* second):AliAODLogicalOperPairCut(first,second){}
     virtual   ~AliAODOrPairCut(){}
     Bool_t    Pass(AliAODPair *p) const;
     ClassDef(AliAODOrPairCut,1)
};
/******************************************************************/

class AliAODAndPairCut: public AliAODLogicalOperPairCut
{
   public:
     AliAODAndPairCut(){}
     AliAODAndPairCut(AliAODBasePairCut* first, AliAODBasePairCut* second):AliAODLogicalOperPairCut(first,second){}
     virtual   ~AliAODAndPairCut(){}
     Bool_t    Pass(AliAODPair *p) const;
     ClassDef(AliAODAndPairCut,1)
};

#endif
