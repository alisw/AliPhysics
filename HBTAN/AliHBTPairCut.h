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
  ClassDef(AliHBTPairCut,2)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliHBTEmptyPairCut:  public AliHBTPairCut
{
  //Empty - it passes possitively all particles - it means returns always False
  //Class describing cut on pairs of particles
 public:
  AliHBTEmptyPairCut(){};
  AliHBTEmptyPairCut(const AliHBTEmptyPairCut& in):AliHBTPairCut(in){};
  virtual ~AliHBTEmptyPairCut(){};
  
  Bool_t Pass(AliHBTPair*) const {return kFALSE;} //accpept everything
  Bool_t IsEmpty() const {return kTRUE;}
  
  ClassDef(AliHBTEmptyPairCut,1)
};



/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliHbtBasePairCut: public TObject
{
  //This class defines the range of some property - pure virtual
  //Property is coded by AliHBTCutTypes type
   
 public:
     
  AliHbtBasePairCut(Double_t min = 0.0, Double_t max = 0.0, AliHBTPairCutProperty prop= kHbtPairCutPropNone):
    fMin(min),fMax(max),fProperty(prop){}
  
  virtual   ~AliHbtBasePairCut(){}
     
  virtual Bool_t    Pass(AliHBTPair* pair) const;
  
  void      SetRange(Double_t min, Double_t max){fMin = min; fMax = max;}
  
  void      SetMinimum(Double_t min){fMin = min;}
  void      SetMaximum(Double_t max){fMax = max;}
  
  Double_t  GetMinimum() const {return fMin;}
  Double_t  GetMaximum() const {return fMax;}
  
  AliHBTPairCutProperty GetProperty() const {return fProperty;}
  
 protected:
  virtual Double_t  GetValue(AliHBTPair* pair) const = 0;
  
  Double_t fMin; // Lower boundary of the range
  Double_t fMax; // Upper boundary of the range
  
  AliHBTPairCutProperty fProperty; // The property itself
  
  ClassDef(AliHbtBasePairCut,1)
 
 };
/******************************************************************/

inline Bool_t AliHbtBasePairCut::Pass(AliHBTPair* pair) const
{
  //checks if pair proprty is in range
  //null pointer check is made by AliHBTPairCut, so here is unnecesary
  
  Double_t value = GetValue(pair);
  if ( (value > fMin) && (value <fMax ) ) return kFALSE; //accepted
  else return kTRUE; //rejected
}
/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliHBTQInvCut: public AliHbtBasePairCut
{
 public:
  AliHBTQInvCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBasePairCut(min,max,kHbtPairCutPropQInv){}
  virtual ~AliHBTQInvCut(){}
 protected:
  virtual Double_t  GetValue(AliHBTPair* pair) const {return pair->GetQInv();}
  
  ClassDef(AliHBTQInvCut,1)
 };
/******************************************************************/

class AliHBTKtCut: public AliHbtBasePairCut {
 public:
  AliHBTKtCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBasePairCut(min,max,kHbtPairCutPropKt){}
  virtual ~AliHBTKtCut(){}
 protected:
  virtual Double_t  GetValue(AliHBTPair* pair) const {return pair->GetKt();}

  ClassDef(AliHBTKtCut,1)
 };
/******************************************************************/

class AliHBTKStarCut: public AliHbtBasePairCut
{
 public:
  AliHBTKStarCut(Double_t min = 0.0, Double_t max = 0.0):AliHbtBasePairCut(min,max,kHbtPairCutPropKStar){}
  virtual ~AliHBTKStarCut(){}
 protected:
  virtual Double_t  GetValue(AliHBTPair* pair) const {return pair->GetKStar();}

  ClassDef(AliHBTKStarCut,1)
};
/******************************************************************/

class AliHBTQSideLCMSCut: public AliHbtBasePairCut
{
 public:
  AliHBTQSideLCMSCut(Double_t min = 0.0, Double_t max = 0.0):
    AliHbtBasePairCut(min,max,kHbtPairCutPropQSideLCMS){}
  virtual ~AliHBTQSideLCMSCut(){}
 protected:
  virtual Double_t  GetValue(AliHBTPair* pair) const 
    {return pair->GetQSideLCMS();}

  ClassDef(AliHBTQSideLCMSCut,1)
};
/******************************************************************/


class AliHBTQOutLCMSCut: public AliHbtBasePairCut
{
 public:
  AliHBTQOutLCMSCut(Double_t min = 0.0, Double_t max = 0.0):
    AliHbtBasePairCut(min,max,kHbtPairCutPropQOutLCMS){}
  virtual ~AliHBTQOutLCMSCut(){}
 protected:
  virtual Double_t  GetValue(AliHBTPair* pair) const 
    {return pair->GetQOutLCMS();}
  
  ClassDef(AliHBTQOutLCMSCut,1)
};
/******************************************************************/

class AliHBTQLongLCMSCut: public AliHbtBasePairCut
{
 public:
  AliHBTQLongLCMSCut(Double_t min = 0.0, Double_t max = 0.0):
    AliHbtBasePairCut(min,max,kHbtPairCutPropQLongLCMS){}
  virtual ~AliHBTQLongLCMSCut(){}
 protected:
  virtual Double_t  GetValue(AliHBTPair* pair) const 
    {return pair->GetQLongLCMS();}

  ClassDef(AliHBTQLongLCMSCut,1)
};
/******************************************************************/

class AliHBTDeltaPhiCut: public AliHbtBasePairCut
{
 public:
  AliHBTDeltaPhiCut(Double_t min = 0.0, Double_t max = 0.0):
    AliHbtBasePairCut(min,max,kHbtPairCutPropDeltaPhi){}
  virtual ~AliHBTDeltaPhiCut(){}
 protected:
  virtual Double_t  GetValue(AliHBTPair* pair) const 
    {return TMath::Abs(pair->GetDeltaPhi());}

  ClassDef(AliHBTDeltaPhiCut,1)
};
/******************************************************************/

class AliHBTDeltaThetaCut: public AliHbtBasePairCut
{
 public:
  AliHBTDeltaThetaCut(Double_t min = 0.0, Double_t max = 0.0):
    AliHbtBasePairCut(min,max,kHbtPairCutPropDeltaTheta){}
  virtual ~AliHBTDeltaThetaCut(){}
 protected:
  virtual Double_t  GetValue(AliHBTPair* pair) const 
    {return TMath::Abs(pair->GetDeltaTheta());}

  ClassDef(AliHBTDeltaThetaCut,1)
};
/******************************************************************/

class AliHBTCluterOverlapCut: public AliHbtBasePairCut
{
 public:
  AliHBTCluterOverlapCut(Double_t min = 0.0, Double_t max = 1e5):
    AliHbtBasePairCut(min,max,kHbtPairCutPropClOverlap){}
  virtual ~AliHBTCluterOverlapCut(){}

 protected:
  virtual Double_t  GetValue(AliHBTPair* pair) const;
  ClassDef(AliHBTCluterOverlapCut,1)
};
/******************************************************************/
  
class AliHBTAvSeparationCut: public AliHbtBasePairCut
{
 public:
  AliHBTAvSeparationCut(Double_t min = 0.0, Double_t max = 1e5):
    AliHbtBasePairCut(min,max,kHbtPairCutPropAvSepar){}
  virtual ~AliHBTAvSeparationCut(){}
  
 protected:
  virtual Double_t  GetValue(AliHBTPair* pair) const;
  ClassDef(AliHBTAvSeparationCut,1)
};
/******************************************************************/
  
class AliHBTSeparationCut: public AliHbtBasePairCut
{
 public:
  AliHBTSeparationCut(Double_t min = 0.0, Double_t max = 1e5, Int_t point = 0):
    AliHbtBasePairCut(min,max,kHbtPairCutPropSepar),fPoint(point){}
  virtual ~AliHBTSeparationCut(){}
  
 protected:
  Int_t fPoint;//index of the point that distance should be measured
  virtual Double_t  GetValue(AliHBTPair* pair) const;
  ClassDef(AliHBTSeparationCut,1)
};
/******************************************************************/
  
class AliHBTITSSeparationCut: public AliHbtBasePairCut
{
//Anti merging cut for the first layer of pixels
 public:
  AliHBTITSSeparationCut(Int_t layer = 0, Double_t deltarphi = 0.01, Double_t deltaz = 0.08):
    AliHbtBasePairCut(deltarphi,deltaz,kHbtPairCutPropPixelSepar),fLayer(layer){}
  virtual ~AliHBTITSSeparationCut(){}
  Bool_t   Pass(AliHBTPair* pair) const;
  Int_t    GetLayer() const {return fLayer;}
 protected:
  Int_t fLayer;//index of the layer that distance should be measured 0: 1st pixels
  virtual Double_t  GetValue(AliHBTPair* /*pair*/) const {return 0.0;}//not used
  ClassDef(AliHBTITSSeparationCut,1)
};
/******************************************************************/

class AliHBTOutSideSameSignCut: public AliHbtBasePairCut
{
 public:
  AliHBTOutSideSameSignCut(){}
  virtual ~AliHBTOutSideSameSignCut(){}
  virtual Bool_t Pass(AliHBTPair *p) const;
 protected:
  virtual Double_t  GetValue(AliHBTPair* /*pair*/) const {return 0.0;}
  ClassDef(AliHBTOutSideSameSignCut,1)
};
/******************************************************************/

class AliHBTOutSideDiffSignCut: public AliHbtBasePairCut
{
 public:
  AliHBTOutSideDiffSignCut(){}
  virtual ~AliHBTOutSideDiffSignCut(){}
  virtual Bool_t Pass(AliHBTPair *p) const;
 protected:
  virtual Double_t  GetValue(AliHBTPair* /*pair*/) const {return 0.0;}
  ClassDef(AliHBTOutSideDiffSignCut,1)
};
/******************************************************************/

class AliHBTLogicalOperPairCut:  public AliHbtBasePairCut
 {
   public:
     AliHBTLogicalOperPairCut();
     AliHBTLogicalOperPairCut(AliHbtBasePairCut* first, AliHbtBasePairCut* second);
     virtual   ~AliHBTLogicalOperPairCut();
   protected:
     Double_t  GetValue(AliHBTPair * /*pair*/) const {MayNotUse("GetValue");return 0.0;}

     AliHbtBasePairCut* fFirst;   //second cut
     AliHbtBasePairCut* fSecond;  //first cut
   private:
    class  AliHBTDummyBasePairCut: public AliHbtBasePairCut
     {
       Double_t  GetValue(AliHBTPair* /*pair*/) const {return 0.0;}
       Bool_t    Pass(AliHBTPair* /*pair*/) const;
     };

    ClassDef(AliHBTLogicalOperPairCut,1)
 };
/******************************************************************/

class AliHBTOrPairCut: public AliHBTLogicalOperPairCut
{
   public:
     AliHBTOrPairCut(){}
     AliHBTOrPairCut(AliHbtBasePairCut* first, AliHbtBasePairCut* second):AliHBTLogicalOperPairCut(first,second){}
     virtual   ~AliHBTOrPairCut(){}
     Bool_t    Pass(AliHBTPair *p) const;
     ClassDef(AliHBTOrPairCut,1)
};
/******************************************************************/

class AliHBTAndPairCut: public AliHBTLogicalOperPairCut
{
   public:
     AliHBTAndPairCut(){}
     AliHBTAndPairCut(AliHbtBasePairCut* first, AliHbtBasePairCut* second):AliHBTLogicalOperPairCut(first,second){}
     virtual   ~AliHBTAndPairCut(){}
     Bool_t    Pass(AliHBTPair *p) const;
     ClassDef(AliHBTAndPairCut,1)
};

#endif
