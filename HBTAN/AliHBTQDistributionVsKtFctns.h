#ifndef ALIHBTQDISTRIBUTIONVSKTFCTN_H
#define ALIHBTQDISTRIBUTIONVSKTFCTN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////////////////////
// 
// class AliHBTQInvDistributionVsKtFctn;    //QInvLCMS   Distribution Vs   Kt
// class AliHBTQOutDistributionVsKtFctn;    //QOutLCMS   Distribution Vs   Kt
// class AliHBTQSideDistributionVsKtFctn;   //QSideLCMS  Distribution Vs   Kt
// class AliHBTQLongDistributionVsKtFctn;   //QLongLCMS  Distribution Vs   Kt
//
// added by Zbigniew.Chajecki@cern.ch
// this classes create distribution functions of pair momentum 
//
/////////////////////////////////////////////////////////////////////////////

#include "AliHBTFunction.h"

/***********************************************************************/
/***********************************************************************/
class AliHBTQOutDistributionVsKtFctn: public AliHBTOnePairFctn2D
 {
  public: 
   AliHBTQOutDistributionVsKtFctn(Int_t nXbins = 200, Double_t maxXval = 1., Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   virtual ~AliHBTQOutDistributionVsKtFctn(){}
   TH1* GetResult();
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQOutLCMS();
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
   TH1* GetResult();
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQSideLCMS();
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
   TH1* GetResult();
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQLongLCMS();
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
   TH1* GetResult();
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQInv();
     x = partpair->GetKt();
    }
  protected:
    ClassDef(AliHBTQInvDistributionVsKtFctn,1)
 };


#endif

