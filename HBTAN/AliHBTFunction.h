//Piotr Skowronski@cern.ch

#ifndef ALIHBTFUNCTION_H
#define ALIHBTFUNCTION_H

#include "AliHBTParticleCut.h"
#include "AliHBTPairCut.h"
#include "AliHBTPair.h"

#include <TH2.h>
#include <TH3.h>


class AliHBTAnalysis;

class AliHBTFunction: public TNamed
//Abstract base class for HBT functions
{
  public:
    AliHBTFunction();
    virtual ~AliHBTFunction(){}
    
    virtual TH1* GetNumerator() =0;
    virtual TH1* GetDenominator() =0;
    virtual TH1* GetResult() = 0;
    
    virtual void Write();
    
    TH1* GetRatio(Double_t normfactor = 1.0);
    void Rename(const Char_t * name); //renames the function and histograms ==title is the same that name
    void Rename(const Char_t * name, const Char_t * title); //renames and retitle the function and histograms
    
    void SetPairCut(AliHBTPairCut*);
    
    virtual AliHBTPair* CheckPair(AliHBTPair* pair);
    
  protected:
    
    AliHBTPairCut*      fPairCut;
    
  public:  
   ClassDef(AliHBTFunction,1)
};
inline AliHBTPair* AliHBTFunction::CheckPair(AliHBTPair* pair)
{
  //check if pair and both particles meets the cut criteria
  if(fPairCut->Pass(pair)) //if the pair is BAD
   {//it is BAD 
    pair = pair->GetSwapedPair();
    if(fPairCut->Pass(pair)) //so try reverse combination
     { 
       return 0x0;//it is BAD as well - so return
     }
   }
  return pair; 
}



/******************************************************************/
/******************************************************************/
/******************************************************************/
class AliHBTTwoPartFctn: public AliHBTFunction
{
  public:
    AliHBTTwoPartFctn(){}
    virtual ~AliHBTTwoPartFctn(){}
    
    virtual void ProcessSameEventParticles(AliHBTPair* pair) = 0;
    virtual void ProcessDiffEventParticles(AliHBTPair* pair) = 0;
    
    
  protected:
  public:  
   ClassDef(AliHBTTwoPartFctn,1)
  
};
/******************************************************************/
/******************************************************************/
/******************************************************************/
class AliHBTFourPartFctn: public AliHBTFunction
{
  public:
    AliHBTFourPartFctn(){};
    virtual ~AliHBTFourPartFctn(){};
    
    virtual void 
    ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair) = 0;
    virtual void 
    ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair) = 0;
	     
  protected:
  public:  
   ClassDef(AliHBTFourPartFctn,1)
  
};
/******************************************************************/
/******************************************************************/
/******************************************************************/


class AliHBTTwoPartFctn1D: public AliHBTTwoPartFctn
{
 public:
  AliHBTTwoPartFctn1D(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual ~AliHBTTwoPartFctn1D();
  
  
  TH1* GetNumerator(){return fNumerator;}
  TH1* GetDenominator(){return fDenominator;}

  void ProcessSameEventParticles(AliHBTPair* pair);
  void ProcessDiffEventParticles(AliHBTPair* pair);
  
 protected:
  //retruns velue to be histogrammed
  virtual Double_t GetValue(AliHBTPair* pair) = 0; 
  
  TH1D* fNumerator;
  TH1D* fDenominator;
  
 public:
  ClassDef(AliHBTTwoPartFctn1D,1)
};

inline void 
AliHBTTwoPartFctn1D::ProcessSameEventParticles(AliHBTPair* pair)
{
 //Fills the numerator
   pair = CheckPair(pair);
   if(pair) fNumerator->Fill(GetValue(pair));
}
  
  
inline void
AliHBTTwoPartFctn1D::ProcessDiffEventParticles(AliHBTPair* pair)
 {
  //fills denumerator
   pair = CheckPair(pair);
   if(pair) fDenominator->Fill(GetValue(pair));

  }
/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliHBTTwoPartFctn2D: public AliHBTTwoPartFctn
{
 public:
  AliHBTTwoPartFctn2D(Int_t nXbins = 200, Double_t maxXval = 1.5, Double_t minXval = 0.0, 
                      Int_t nYbins = 200, Double_t maxYval = .15, Double_t minYval =-0.15);
  ~AliHBTTwoPartFctn2D();
  
  TH1* GetNumerator(){return fNumerator;}
  TH1* GetDenominator(){return fDenominator;}
  
  void ProcessSameEventParticles(AliHBTPair* pair);
  void ProcessDiffEventParticles(AliHBTPair* pair);
 
  virtual void GetValues(AliHBTPair* pair, Double_t&, Double_t&) = 0;
  
 protected:
  TH2D* fNumerator;
  TH2D* fDenominator;
  
 public:
  ClassDef(AliHBTTwoPartFctn2D,1)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

class AliHBTTwoPartFctn3D: public AliHBTTwoPartFctn
{
 public:
  AliHBTTwoPartFctn3D(Int_t nXbins = 200, Double_t maxXval = 1.5, Double_t minXval = 0.0, 
                      Int_t nYbins = 200, Double_t maxYval = .15, Double_t minYval =-0.15, 
                      Int_t nZbins = 200, Double_t maxZval = .15, Double_t minZval =-0.15);
	    
  virtual ~AliHBTTwoPartFctn3D();

  TH1* GetNumerator(){return fNumerator;}
  TH1* GetDenominator(){return fDenominator;}

 protected:
  TH3D* fNumerator;
  TH3D* fDenominator;
 public:
  ClassDef(AliHBTTwoPartFctn3D,1)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/



/******************************************************************/
/******************************************************************/
/******************************************************************/
class AliHBTFourPartFctn2D: public AliHBTFourPartFctn
{
 public:
  AliHBTFourPartFctn2D(Int_t nXbins = 200, Double_t maxXval = 1.5, Double_t minXval = 0.0, 
                       Int_t nYbins = 200, Double_t maxYval = .15, Double_t minYval =-0.15);
  ~AliHBTFourPartFctn2D();
  
  TH1* GetNumerator(){return fNumerator;}
  TH1* GetDenominator(){return fDenominator;}
  
  void ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
 
  virtual void GetValues(AliHBTPair*,AliHBTPair*, Double_t&, Double_t&) = 0;
  
 protected:
  TH2D* fNumerator;
  TH2D* fDenominator;
  
 public:
  ClassDef(AliHBTFourPartFctn2D,1)
};


/******************************************************************/
/******************************************************************/
/******************************************************************/

/******************************************************************/
/******************************************************************/
/******************************************************************/



#endif
