#ifndef ALIHBTRDISTRIBUTIONS_H
#define ALIHBTRDISTRIBUTIONS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////////////////////
// AliHBTRStarDistribution
// AliHBTRDistribution
// added by Zbigniew.Chajecki@cern.ch
// this classes create distribution functions of pair momentum 
//
/////////////////////////////////////////////////////////////////////////////

#include "AliHBTFunction.h"

class AliHBTRStarDistribution: public AliHBTOnePairFctn1D
{
  public:
    AliHBTRStarDistribution(Int_t nXbins = 500, Double_t maxXval = 5e-11, Double_t minXval = 0.);
    virtual ~AliHBTRStarDistribution(){}
    TH1* GetResult();
  protected:
    Double_t GetValue(AliHBTPair* partpair) const
    {
      return partpair->GetRStar();
    }
   
  private:
   ClassDef(AliHBTRStarDistribution,1)
};

/***********************************************************************/
/***********************************************************************/

class AliHBTRDistribution: public AliHBTOnePairFctn1D
{
  public:
    AliHBTRDistribution(Int_t nXbins = 500, Double_t maxXval = 5e-11, Double_t minXval = 0.);
    virtual ~AliHBTRDistribution(){}
    TH1* GetResult();
  protected:
    Double_t GetValue(AliHBTPair* partpair) const
    {
      return partpair->GetR();
    }
   
  private:
   ClassDef(AliHBTRDistribution,1)
};
#endif

