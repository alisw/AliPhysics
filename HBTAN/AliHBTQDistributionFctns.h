#ifndef ALIHBTQOUTDISTRIBUTIONVSKTLFCTN_H
#define ALIHBTQOUTDISTRIBUTIONVSKTLFCTN_H

/////////////////////////////////////////////////////////////////////////////
// 
// class AliHBTQInvDistributionVsKtFctn;    //QInvCMSLC   Distribution Vs   Kt
// class AliHBTQOutDistributionVsKtFctn;    //QOutCMSLC   Distribution Vs   Kt
// class AliHBTQSideDistributionVsKtFctn;   //QSideCMSLC  Distribution Vs   Kt
// class AliHBTQLongDistributionVsKtFctn;   //QLongCMSLC  Distribution Vs   Kt

// class AliHBTQOutDistributionVsQInvFctn;    //QOutCMSLC   Distribution Vs   QInv
// class AliHBTQSideDistributionVsQInvFctn;   //QSideCMSLC  Distribution Vs   QInv
// class AliHBTQLongDistributionVsQInvFctn;   //QLongCMSLC  Distribution Vs   QInv
// class AliHBTPtDiffDistributionVsQInvFctn;
//
// added by Zbigniew.Chajecki@cern.ch
// this classes create distribution functions of pair momentum 
//
/////////////////////////////////////////////////////////////////////////////

class AliHBTQInvDistributionVsKtFctn;    //QInvCMSLC   Distribution Vs   Kt
class AliHBTQOutDistributionVsKtFctn;    //QOutCMSLC   Distribution Vs   Kt
class AliHBTQSideDistributionVsKtFctn;   //QSideCMSLC  Distribution Vs   Kt
class AliHBTQLongDistributionVsKtFctn;   //QLongCMSLC  Distribution Vs   Kt

class AliHBTQOutDistributionVsQInvFctn;    //QOutCMSLC   Distribution Vs   QInv
class AliHBTQSideDistributionVsQInvFctn;   //QSideCMSLC  Distribution Vs   QInv
class AliHBTQLongDistributionVsQInvFctn;   //QLongCMSLC  Distribution Vs   QInv
class AliHBTPtDiffDistributionVsQInvFctn;

#include "AliHBTFunction.h"

/***********************************************************************/
/***********************************************************************/
class AliHBTQOutDistributionVsKtFctn: public AliHBTOnePairFctn2D
 {
  public: 
   AliHBTQOutDistributionVsKtFctn(Int_t nXbins = 200, Double_t maxXval = 1., Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   virtual ~AliHBTQOutDistributionVsKtFctn(){}
   TH1* GetResult(){return this->GetNumerator();}
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQOutCMSLC();
     x = partpair->GetKt();
    }
  protected:
    ClassDef(AliHBTQOutDistributionVsKtFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
class AliHBTQSideDistributionVsKtFctn: public AliHBTOnePairFctn2D
 {
  public: 
   AliHBTQSideDistributionVsKtFctn(Int_t nXbins = 200, Double_t maxXval = 1.2, Double_t minXval = -0.1, 
                             Int_t nYbins = 500, Double_t maxYval = 1.2, Double_t minYval =-1.2);
   virtual ~AliHBTQSideDistributionVsKtFctn(){}
   TH1* GetResult(){return this->GetNumerator();}
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQSideCMSLC();
     x = partpair->GetKt();
    }
  protected:
    ClassDef(AliHBTQSideDistributionVsKtFctn,1)
 };
/***********************************************************************/
/***********************************************************************/

class AliHBTQLongDistributionVsKtFctn: public AliHBTOnePairFctn2D
 {
  public: 
   AliHBTQLongDistributionVsKtFctn(Int_t nXbins = 200, Double_t maxXval = 1.2, Double_t minXval = -0.1, 
                             Int_t nYbins = 500, Double_t maxYval = 1.2, Double_t minYval =-1.2);
   virtual ~AliHBTQLongDistributionVsKtFctn(){}
   TH1* GetResult(){return this->GetNumerator();}
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQLongCMSLC();
     x = partpair->GetKt();
    }
  protected:
    ClassDef(AliHBTQLongDistributionVsKtFctn,1)
 };
/***********************************************************************/
/***********************************************************************/

class AliHBTQInvDistributionVsKtFctn: public AliHBTOnePairFctn2D
 {
  public: 
   AliHBTQInvDistributionVsKtFctn(Int_t nXbins = 200, Double_t maxXval = 1.2, Double_t minXval = -0.1, 
                             Int_t nYbins = 500, Double_t maxYval = 1.2, Double_t minYval =-1.2);
   virtual ~AliHBTQInvDistributionVsKtFctn(){}
   TH1* GetResult(){return this->GetNumerator();}
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQInv();
     x = partpair->GetKt();
    }
  protected:
    ClassDef(AliHBTQInvDistributionVsKtFctn,1)
 };

/***********************************************************************/
/***********************************************************************/
class AliHBTQOutDistributionVsQInvFctn: public AliHBTOnePairFctn2D
 {
  public: 
   AliHBTQOutDistributionVsQInvFctn(Int_t nXbins = 200, Double_t maxXval = 1., Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   virtual ~AliHBTQOutDistributionVsQInvFctn(){}
   TH1* GetResult(){return this->GetNumerator();}
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQOutCMSLC();
     x = partpair->GetQInv();
    }
  protected:
    ClassDef(AliHBTQOutDistributionVsQInvFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
class AliHBTQSideDistributionVsQInvFctn: public AliHBTOnePairFctn2D
 {
  public: 
   AliHBTQSideDistributionVsQInvFctn(Int_t nXbins = 200, Double_t maxXval = 1.2, Double_t minXval = -0.1, 
                             Int_t nYbins = 500, Double_t maxYval = 1.2, Double_t minYval =-1.2);
   virtual ~AliHBTQSideDistributionVsQInvFctn(){}
   TH1* GetResult(){return this->GetNumerator();}
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQSideCMSLC();
     x = partpair->GetQInv();
    }
  protected:
    ClassDef(AliHBTQSideDistributionVsQInvFctn,1)
 };
/***********************************************************************/
/***********************************************************************/

class AliHBTQLongDistributionVsQInvFctn: public AliHBTOnePairFctn2D
 {
  public: 
   AliHBTQLongDistributionVsQInvFctn(Int_t nXbins = 200, Double_t maxXval = 1.2, Double_t minXval = -0.1, 
                             Int_t nYbins = 500, Double_t maxYval = 1.2, Double_t minYval =-1.2);
   virtual ~AliHBTQLongDistributionVsQInvFctn(){}
   TH1* GetResult(){return this->GetNumerator();}
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQLongCMSLC();
     x = partpair->GetQInv();
    }
  protected:
    ClassDef(AliHBTQLongDistributionVsQInvFctn,1)
 };
/***********************************************************************/
/***********************************************************************/
class AliHBTPtDiffDistributionVsQInvFctn: public AliHBTOnePairFctn2D
 {
  public: 
   AliHBTPtDiffDistributionVsQInvFctn(Int_t nXbins = 800, Double_t maxXval = 4.0, Double_t minXval = 0., 
                             Int_t nYbins = 500, Double_t maxYval = 0.1, Double_t minYval =-0.1);
   virtual ~AliHBTPtDiffDistributionVsQInvFctn(){}
   TH1* GetResult(){return this->GetNumerator();}
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->Particle1()->Pt() - partpair->Particle2()->Pt();
     x = partpair->GetQInv();
    }
  protected:
    ClassDef(AliHBTPtDiffDistributionVsQInvFctn,1)
 };
/***********************************************************************/
/***********************************************************************/

#endif
