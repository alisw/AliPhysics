#ifndef ALIHBTQDISTRIBUTIONVSQINVFCTN_H
#define ALIHBTQDISTRIBUTIONVSQINVFCTN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////////////////////
// 
// class AliHBTQOutDistributionVsQInvFctn;    //QOutLCMS   Distribution Vs   QInv
// class AliHBTQSideDistributionVsQInvFctn;   //QSideLCMS  Distribution Vs   QInv
// class AliHBTQLongDistributionVsQInvFctn;   //QLongLCMS  Distribution Vs   QInv
// class AliHBTPtDiffDistributionVsQInvFctn;
//
// added by Zbigniew.Chajecki@cern.ch
// this classes create distribution functions of pair momentum 
//
/////////////////////////////////////////////////////////////////////////////

#include "AliHBTFunction.h"

class AliHBTQOutDistributionVsQInvFctn: public AliHBTOnePairFctn2D
 {
  public: 
   AliHBTQOutDistributionVsQInvFctn(Int_t nXbins = 200, Double_t maxXval = 1., Double_t minXval = 0.0, 
                             Int_t nYbins = 500, Double_t maxYval = .15, Double_t minYval =-0.15);
   virtual ~AliHBTQOutDistributionVsQInvFctn(){}
   TH1* GetResult();
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQOutLCMS();
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
   TH1* GetResult();
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQSideLCMS();
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
   TH1* GetResult();
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->GetQLongLCMS();
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
   TH1* GetResult();
   void GetValues(AliHBTPair* partpair, Double_t& x, Double_t& y) const
    {
     y = partpair->Particle1()->Pt() - partpair->Particle2()->Pt();
     x = partpair->GetQInv();
    }
  protected:
    ClassDef(AliHBTPtDiffDistributionVsQInvFctn,1)
 };
 
#endif

