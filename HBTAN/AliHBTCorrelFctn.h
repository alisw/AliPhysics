#ifndef ALIHBTCORRELFUNCTION_H
#define ALIHBTCORRELFUNCTION_H
//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//
// class AliHBTQInvCorrelFctn
// class AliHBTQOutLCMSCorrelFctn
// class AliHBTQLongLCMSCorrelFctn
// class AliHBTQSideLCMSCorrelFctn
// class AliHBTInvMassCorrelFctn
// class AliHBTTwoKStarCorrelFctn
//
// Set of functions:
//   Q Invaraint Correlation Function
//   Invariant Mass Function
//
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
// Piotr.Skowronski@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#include "AliHBTFunction.h"
#include <Riostream.h>

/*************************************************************************************/ 
class AliHBTQInvCorrelFctn: public AliHBTOnePairFctn1D, public AliHBTCorrelFunction
{
//Q Invaraint Correlation Function
//1D two particle function 
 public:
   AliHBTQInvCorrelFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
   virtual ~AliHBTQInvCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair) const {return pair->GetQInv();}
 private:  
   ClassDef(AliHBTQInvCorrelFctn,2)
 
};
/*************************************************************/

class AliHBTOutSideLongFctn: public AliHBTOnePairFctn3D, public AliHBTCorrelFunction
{

  public:
    AliHBTOutSideLongFctn(Int_t nXbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0,
                          Int_t nYbins = 100, Double_t maxYval = 0.15, Double_t minYval = 0.0,
	      Int_t nZbins = 100, Double_t maxZval = 0.15, Double_t minZval = 0.0);
    virtual  ~AliHBTOutSideLongFctn(){}

    TH1* GetResult();
    void UseAbsoluteValues(Bool_t flag){fAbs = flag;}
 
  protected:
    void GetValues(AliHBTPair* pair, Double_t& x, Double_t& y, Double_t& z) const;
    
    Bool_t fAbs;//flag indicating if absolute values of qout, qside and qlong should be histogrammed
  ClassDef(AliHBTOutSideLongFctn,1)
};

/*************************************************************************************/ 

class AliHBTQOutLCMSCorrelFctn: public AliHBTOnePairFctn1D, public AliHBTCorrelFunction
{
//Q OutLCMSaraint Correlation Function
//1D two particle function 
 public:
   AliHBTQOutLCMSCorrelFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
   virtual ~AliHBTQOutLCMSCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair) const {return pair->GetQOutLCMS();}
 private:  
    ClassDef(AliHBTQOutLCMSCorrelFctn,2)
};
/*************************************************************************************/ 

class AliHBTQLongLCMSCorrelFctn: public AliHBTOnePairFctn1D, public AliHBTCorrelFunction
{
//Q LongLCMSaraint Correlation Function
//1D two particle function 
 public:
   AliHBTQLongLCMSCorrelFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
   virtual ~AliHBTQLongLCMSCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair) const {return pair->GetQLongLCMS();}
 private:  
    ClassDef(AliHBTQLongLCMSCorrelFctn,2)
};
/*************************************************************************************/ 

class AliHBTQSideLCMSCorrelFctn: public AliHBTOnePairFctn1D, public AliHBTCorrelFunction
{
//Q SideLCMSaraint Correlation Function
//1D two particle function 
 public:
   AliHBTQSideLCMSCorrelFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
   virtual ~AliHBTQSideLCMSCorrelFctn(){}
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair) const {return pair->GetQSideLCMS();}
 private:  
    ClassDef(AliHBTQSideLCMSCorrelFctn,2)
};
/*************************************************************************************/ 

class AliHBTInvMassCorrelFctn: public AliHBTOnePairFctn1D
{
//   Invariant Mass Function 
 public:
   AliHBTInvMassCorrelFctn(Int_t nbins = 2000, Double_t maxXval = 2., Double_t minXval = 0.0);
   virtual ~AliHBTInvMassCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair) const { return pair->GetInvMass();}
 private:  
    ClassDef(AliHBTInvMassCorrelFctn,1)
};

/*************************************************************************************/ 

class AliHBTTwoKStarCorrelFctn: public AliHBTOnePairFctn1D, public AliHBTCorrelFunction
{
//   Correlation Function of 2*KStar
 public:
   AliHBTTwoKStarCorrelFctn(Int_t nbins = 200, Double_t maxXval = 0.15, Double_t minXval = 0.0);
   virtual ~AliHBTTwoKStarCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair) const { return 2.0*pair->GetKStar();}
 private:  
    ClassDef(AliHBTTwoKStarCorrelFctn,2)
};

/*************************************************************************************/ 

class AliHBTAvSeparCorrelFctn: public AliHBTOnePairFctn1D, public AliHBTCorrelFunction
{
//   Correlation Function of 2*KStar
 public:
   AliHBTAvSeparCorrelFctn(Int_t nbins = 200, Double_t maxXval = 30, Double_t minXval = 0.0);
   virtual ~AliHBTAvSeparCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair) const { return pair->GetAvarageDistance();}
 private:  
    ClassDef(AliHBTAvSeparCorrelFctn,2)
};

#endif
