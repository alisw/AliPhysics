#ifndef ALIHBTFUNCTION_H
#define ALIHBTFUNCTION_H

/* Id: $ */

///////////////////////////////////////////////////////
//                                                   //
// AliHBTFunction                                    //
//                                                   //
// Abstract Base Calss for all the function classes  //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

#include <TH1.h>
#include <TH2D.h>
#include <TH3D.h>

#include "AliHBTPairCut.h"
#include "AliHBTPair.h"


class AliHBTAnalysis;
class AliHBTParticleCut;

class AliHBTFunction: public TNamed
{
  public:
    AliHBTFunction();
    AliHBTFunction(const char* name, const char* title);
    AliHBTFunction(const AliHBTFunction & source);
    
    virtual ~AliHBTFunction();
    
    AliHBTFunction & operator= (const AliHBTFunction & source);

    virtual TH1* GetNumerator() const = 0;
    virtual TH1* GetDenominator() const = 0;
    virtual TH1* GetResult() = 0;

    virtual void WriteFunction();
    virtual void InitFunction();
    
    TH1* GetRatio(Double_t normfactor = 1.0);
    void Rename(const Char_t * name); //renames the function and histograms ==title is the same that name
    void Rename(const Char_t * name, const Char_t * title); //renames and retitle the function and histograms
    
    void SetPairCut(AliHBTPairCut* cut);
    
    virtual AliHBTPair* CheckPair(AliHBTPair* pair);
    void  SetWriteNumAndDen(Bool_t flag = kFALSE){fWriteNumAndDen = flag;}
  protected:
    virtual void BuildHistos() = 0;//builds default histograms
    AliHBTPairCut*   fPairCut;     //pair cut
    Bool_t           fWriteNumAndDen; //flag indicating whether numerator and denominator should be writted together with a result
    ClassDef(AliHBTFunction,3)
};
/******************************************************************/
inline AliHBTPair* AliHBTFunction::CheckPair(AliHBTPair* pair)
{
  //check if pair and both particles meets the cut criteria
  if(fPairCut->Pass(pair)) //if the pair is BAD
   {//it is BAD 
    pair = pair->GetSwapedPair();
    if(pair)
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
class AliHBTCorrelFunction
{
  public:
    AliHBTCorrelFunction():fRatio(0x0){}
    AliHBTCorrelFunction(const AliHBTCorrelFunction& in):fRatio((in.fRatio)?(TH1*)in.fRatio->Clone():0x0){}
    virtual ~AliHBTCorrelFunction(){delete fRatio;}
    
    AliHBTCorrelFunction& operator=(const AliHBTCorrelFunction& in);
   
  protected:
    TH1* fRatio;//!pointer to the ratio(result)
    
  ClassDef(AliHBTCorrelFunction,1)
};

/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTOnePairFctn                                 //
//                                                   //
// Abstract Base Calss for Functions that need       //
// one pair to fill function                         //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

class AliHBTOnePairFctn
{
  public:
    AliHBTOnePairFctn(){}
    virtual ~AliHBTOnePairFctn(){}
    
    virtual void ProcessSameEventParticles(AliHBTPair* pair) = 0;
    virtual void ProcessDiffEventParticles(AliHBTPair* pair) = 0;

    virtual void Init() = 0;
    virtual void Write() = 0;
    virtual const char* Name() = 0;
    
   ClassDef(AliHBTOnePairFctn,2)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTTwoPairFctn                                 //
//                                                   //
// Abstract Base Calss for Functions that need       //
// two pairs to fill function,                       //
// one reconstructed track and corresponding         //
// simulated pair                                    //
// Basically resolution functions                    //
// Lednicky's algorithm uses that as well            //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

class AliHBTTwoPairFctn
{
  public:
    AliHBTTwoPairFctn(){};
    virtual ~AliHBTTwoPairFctn(){};
    
    virtual void 
    ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair) = 0;
    virtual void 
    ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair) = 0;
    
    virtual void Init() = 0;
    virtual void Write() = 0;
    virtual const char* Name() = 0;
	     
   ClassDef(AliHBTTwoPairFctn,2)
  
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTFunction1D                                  //
//                                                   //
// Base Calss for 1-dimensinal Functions             //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////


class AliHBTFunction1D: public AliHBTFunction
{
 public:
  AliHBTFunction1D();//default conmstructor
  AliHBTFunction1D(Int_t nbins, Float_t maxXval, Float_t minXval);
  AliHBTFunction1D(const Char_t *name, const Char_t *title);
  AliHBTFunction1D(const Char_t *name, const Char_t *title,
                   Int_t nbins, Float_t maxXval, Float_t minXval);

  AliHBTFunction1D(const AliHBTFunction1D & source);
  AliHBTFunction1D & operator= (const AliHBTFunction1D& /*source*/);

  virtual ~AliHBTFunction1D();
  
  TH1* GetNumerator() const {return fNumerator;}//returns numerator histogram
  TH1* GetDenominator() const {return fDenominator;}//returns denominator histogram

  Double_t Scale();
  void SetNumberOfBinsToScale(Int_t n = fgkDefaultNBinsToScale){fNBinsToScale = n;}

 protected:
  //returns value to be histogrammed
  virtual void BuildHistos(Int_t nbins, Float_t max, Float_t min);
  virtual void BuildHistos();
  Double_t Scale(TH1D* num,TH1D* den);
  
  TH1D* fNumerator; // Numerator histogram
  TH1D* fDenominator; // Denumerator histogram
  UInt_t fNBinsToScale; // Number of bins to scale

 private:
  //this must be declared before constructors because they are used as a default arguments
  static const Int_t fgkDefaultNBins;//default number of Bins in histograms
  static const Float_t fgkDefaultMin;//Default min value of histograms
  static const Float_t fgkDefaultMax;//Default max value of histograms
  static const UInt_t fgkDefaultNBinsToScale;//Default number of bins used for scaling to tale

  ClassDef(AliHBTFunction1D,2)
};

/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTFunction2D                                  //
//                                                   //
// Base Calss for 2-dimensinal Functions             //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////
 
class AliHBTFunction2D: public AliHBTFunction
{
 public:
  AliHBTFunction2D();

  AliHBTFunction2D(const Char_t *name, const Char_t *title);

  AliHBTFunction2D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval);

  AliHBTFunction2D(const Char_t *name, const Char_t *title,
                      Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval);
	  
  AliHBTFunction2D(const AliHBTFunction2D & source);

  AliHBTFunction2D & operator= (const AliHBTFunction2D & /*source*/);

  virtual ~AliHBTFunction2D();
  
  TH1* GetNumerator() const {return fNumerator;}
  TH1* GetDenominator() const {return fDenominator;}
  
  void SetNumberOfBinsToScale(UInt_t xn = fgkDefaultNBinsToScaleX, 
                              UInt_t yn = fgkDefaultNBinsToScaleY);
  
  Double_t Scale();
 protected:
  virtual void BuildHistos(Int_t nxbins, Float_t xmax, Float_t xmin,
                           Int_t nybins, Float_t ymax, Float_t ymin);
  virtual void BuildHistos();
  
  TH2D* fNumerator; // Numerator histogram
  TH2D* fDenominator; // Denominator histogram
  
  //definition of area used for scaling -> Scale is calculated this 
  //way that after division tale is on 1
  UInt_t fNBinsToScaleX;//number of bins on X axis
  UInt_t fNBinsToScaleY;//number of bins on Y axis

 private:
  //this must be declared before constructors because they are used as a default arguments
  static const Int_t fgkDefaultNBinsX;//default number of Bins in X axis in histograms
  static const Float_t fgkDefaultMinX;//Default min value of X axis in histograms
  static const Float_t fgkDefaultMaxX;//Default max value of X axis inhistograms
  static const Int_t fgkDefaultNBinsY;//default number of Bins in histograms
  static const Float_t fgkDefaultMinY;//Default min value of histograms
  static const Float_t fgkDefaultMaxY;//Default max value of histograms

  static const UInt_t fgkDefaultNBinsToScaleX;//Default number of X bins used for scaling to tale
  static const UInt_t fgkDefaultNBinsToScaleY;//Default number of bins used for scaling to tale

  ClassDef(AliHBTFunction2D,2)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTFunction3D                               //
//                                                   //
// Base Calss for 3-dimensinal Functions that need   //
// one pair to fill function                         //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

class AliHBTFunction3D: public AliHBTFunction
{
 public:
  AliHBTFunction3D();

  AliHBTFunction3D(const Char_t *name, const Char_t *title);

  AliHBTFunction3D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                   Int_t nYbins, Double_t maxYval, Double_t minYval, 
                   Int_t nZbins, Double_t maxZval, Double_t minZval);

  AliHBTFunction3D(const Char_t *name, const Char_t *title,
                   Int_t nXbins, Double_t maxXval, Double_t minXval, 
                   Int_t nYbins, Double_t maxYval, Double_t minYval, 
                   Int_t nZbins, Double_t maxZval, Double_t minZval);

  AliHBTFunction3D(const AliHBTFunction3D & source);
  AliHBTFunction3D & operator= (const AliHBTFunction3D & /*source*/);

  virtual ~AliHBTFunction3D();//destructor

  TH1* GetNumerator() const {return fNumerator;}
  TH1* GetDenominator() const {return fDenominator;}


  void SetNumberOfBinsToScale(UInt_t xn = fgkDefaultNBinsToScaleX, 
                              UInt_t yn = fgkDefaultNBinsToScaleY,
                              UInt_t zn = fgkDefaultNBinsToScaleZ);

  Double_t Scale();

 protected:
  virtual void BuildHistos(Int_t nxbins, Float_t xmax, Float_t xmin,
                           Int_t nybins, Float_t ymax, Float_t ymin,
	       Int_t nzbins, Float_t zmax, Float_t zmin);
  virtual void BuildHistos();
  
  TH3F* fNumerator; // Numerator histogram
  TH3F* fDenominator; // Denominator histogram
  
  //definition of area used for scaling -> Scale is calculated this 
  //way that after division tale is on 1
  UInt_t fNBinsToScaleX;//number of bins on X axis
  UInt_t fNBinsToScaleY;//number of bins on Y axis
  UInt_t fNBinsToScaleZ;//number of bins on Z axis
  
 private:
  //this must be declared before constructors because they are used as a default arguments
  static const Int_t fgkDefaultNBinsX;//default number of Bins in X axis in histograms
  static const Float_t fgkDefaultMinX;//Default min value of X axis in histograms
  static const Float_t fgkDefaultMaxX;//Default max value of X axis inhistograms
  static const Int_t fgkDefaultNBinsY;//default number of Bins in Y axis in histograms
  static const Float_t fgkDefaultMinY;//Default min value of Y axis in histograms
  static const Float_t fgkDefaultMaxY;//Default max value of Y axis inhistograms
  static const Int_t fgkDefaultNBinsZ;//default number of Bins in Z axis in histograms
  static const Float_t fgkDefaultMinZ;//Default min value of Z axis in histograms
  static const Float_t fgkDefaultMaxZ;//Default max value of Z axis inhistograms

  static const UInt_t fgkDefaultNBinsToScaleX;//Default number of X bins used for scaling to tale
  static const UInt_t fgkDefaultNBinsToScaleY;//Default number of Y bins used for scaling to tale
  static const UInt_t fgkDefaultNBinsToScaleZ;//Default number of Z bins used for scaling to tale
  
  ClassDef(AliHBTFunction3D,2)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTOnePairFctn1D                               //
//                                                   //
// Base Calss for 1-dimensinal Functions that need   //
// one pair to fill function                         //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

class AliHBTOnePairFctn1D: public AliHBTOnePairFctn, public AliHBTFunction1D
{
 public:
  AliHBTOnePairFctn1D(){}//default conmstructor
  AliHBTOnePairFctn1D(Int_t nbins, Float_t maxXval, Float_t minXval);
  AliHBTOnePairFctn1D(const Char_t *name, const Char_t *title);
  AliHBTOnePairFctn1D(const Char_t *name, const Char_t *title,
                      Int_t nbins, Float_t maxXval, Float_t minXval);
  virtual ~AliHBTOnePairFctn1D(){}

  void ProcessSameEventParticles(AliHBTPair* pair);
  void ProcessDiffEventParticles(AliHBTPair* pair);
  void Write(){WriteFunction();}
  void Init(){InitFunction();}
  const char* Name(){return GetName();}
  
 protected:
  //retruns velue to be histogrammed
  virtual Double_t GetValue(AliHBTPair* pair) const = 0; 
  ClassDef(AliHBTOnePairFctn1D,2)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTOnePairFctn2D                               //
//                                                   //
// Base Calss for 2-dimensinal Functions that need   //
// one pair to fill function                         //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

class AliHBTOnePairFctn2D: public AliHBTOnePairFctn, public AliHBTFunction2D
{
 public:
  AliHBTOnePairFctn2D(){}

  AliHBTOnePairFctn2D(const Char_t *name, const Char_t *title);

  AliHBTOnePairFctn2D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval);

  AliHBTOnePairFctn2D(const Char_t *name, const Char_t *title,
                      Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval);
	  
  virtual ~AliHBTOnePairFctn2D(){}
  
  void ProcessSameEventParticles(AliHBTPair* pair);
  void ProcessDiffEventParticles(AliHBTPair* pair);
  void Write(){WriteFunction();}
  void Init(){InitFunction();}
  const char* Name(){return GetName();}
 protected:
  virtual void GetValues(AliHBTPair* pair, Double_t& x, Double_t& y) const = 0;
  ClassDef(AliHBTOnePairFctn2D,2)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTOnePairFctn3D                               //
//                                                   //
// Base Calss for 3-dimensinal Functions that need   //
// one pair to fill function                         //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

class AliHBTOnePairFctn3D: public AliHBTOnePairFctn, public AliHBTFunction3D
{
 public:
  AliHBTOnePairFctn3D(){}

  AliHBTOnePairFctn3D(const Char_t *name, const Char_t *title);

  AliHBTOnePairFctn3D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval, 
                      Int_t nZbins, Double_t maxZval, Double_t minZval);

  AliHBTOnePairFctn3D(const Char_t *name, const Char_t *title,
                      Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval, 
                      Int_t nZbins, Double_t maxZval, Double_t minZval);
   
  virtual ~AliHBTOnePairFctn3D(){}//destructor

  void ProcessSameEventParticles(AliHBTPair* pair);
  void ProcessDiffEventParticles(AliHBTPair* pair);
  void Write(){WriteFunction();}
  void Init(){InitFunction();}
  const char* Name(){return GetName();}
 protected:
  virtual void GetValues(AliHBTPair* pair, Double_t& x, Double_t& y, Double_t& z) const = 0;
 ClassDef(AliHBTOnePairFctn3D,2)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTTwoPairFctn1D                               //
//                                                   //
// Base Calss for 1-dimensinal Functions that need   //
// two pair (simulated and reconstructed)            //
// to fill function                                  //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

class AliHBTTwoPairFctn1D: public AliHBTTwoPairFctn, public AliHBTFunction1D
{
 public:
  AliHBTTwoPairFctn1D(){}//default conmstructor
  AliHBTTwoPairFctn1D(Int_t nbins, Float_t maxXval, Float_t minXval);
  AliHBTTwoPairFctn1D(const Char_t *name, const Char_t *title);
  AliHBTTwoPairFctn1D(const Char_t *name, const Char_t *title,
                      Int_t nbins, Float_t maxXval, Float_t minXval);
  virtual ~AliHBTTwoPairFctn1D(){}

  void ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void Write(){WriteFunction();}
  void Init(){InitFunction();}
  const char* Name(){return GetName();}
  
 protected:
  virtual Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair) const = 0;

  ClassDef(AliHBTTwoPairFctn1D,2)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTTwoPairFctn2D                               //
//                                                   //
// Base Calss for 2-dimensinal Functions that need   //
// two pair (simulated and reconstructed)            //
// to fill function                                  //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

class AliHBTTwoPairFctn2D: public AliHBTTwoPairFctn, public AliHBTFunction2D
{
 public:
  AliHBTTwoPairFctn2D(){}

  AliHBTTwoPairFctn2D(const Char_t *name, const Char_t *title);

  AliHBTTwoPairFctn2D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval);

  AliHBTTwoPairFctn2D(const Char_t *name, const Char_t *title,
                      Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval);
	  
  virtual ~AliHBTTwoPairFctn2D(){}
  
  void ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void Write(){WriteFunction();}
  void Init(){InitFunction();}
  const char* Name(){return GetName();}

 protected:
  virtual void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y) const = 0;

  ClassDef(AliHBTTwoPairFctn2D,2)
};
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTTwoPairFctn3D                               //
//                                                   //
// Base Calss for 3-dimensinal Functions that need   //
// two pair (simulated and reconstructed)            //
// to fill function                                  //
//                                                   //
// Piotr.Skowronski@cern.ch                          //
// http://alisoft.cern.ch/people/skowron/analyzer    //
//                                                   //
///////////////////////////////////////////////////////

class AliHBTTwoPairFctn3D: public AliHBTTwoPairFctn, public AliHBTFunction3D
{
 public:
  AliHBTTwoPairFctn3D(){}
  
  AliHBTTwoPairFctn3D(const Char_t *name, const Char_t *title);

  AliHBTTwoPairFctn3D(Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval, 
                      Int_t nZbins, Double_t maxZval, Double_t minZval);

  AliHBTTwoPairFctn3D(const Char_t *name, const Char_t *title,
                      Int_t nXbins, Double_t maxXval, Double_t minXval, 
                      Int_t nYbins, Double_t maxYval, Double_t minYval, 
                      Int_t nZbins, Double_t maxZval, Double_t minZval);
   
  virtual ~AliHBTTwoPairFctn3D(){}//destructor

  void ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void Write(){WriteFunction();}
  void Init(){InitFunction();}
  const char* Name(){return GetName();}

 protected:
  virtual void GetValues(AliHBTPair* trackpair, AliHBTPair* partpair, Double_t& x, Double_t& y, Double_t& z) const = 0;

  ClassDef(AliHBTTwoPairFctn3D,2)
};

#endif
