#ifndef ALIPADMAX_H
#define ALIPADMAX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCclusterKr.h,v 1.8 2007/12/31 16:07:15 matyja Exp $ */

//-------------------------------------------------------
//                    TPC Kr Cluster Class
//
//   Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl
//-------------------------------------------------------

#include "AliTPCvtpr.h"
//_____________________________________________________________________________
class AliPadMax: public AliTPCvtpr{
public:
  AliPadMax();
  ~AliPadMax();
  AliPadMax(AliTPCvtpr vtpr,
	    Short_t beg,Short_t end,Short_t sum);

  AliPadMax(Short_t max,Short_t nt,Short_t np,Short_t nr,
	    Double_t x,Double_t y,Double_t t,
	    Short_t beg,Short_t end,Short_t sum);


  void SetBegin(Short_t q){fBegin=q;}
  void SetEnd(Short_t q){fEnd=q;}
  void SetSum(Short_t q){fSumAdc=q;}

  Short_t GetBegin(){return fBegin;}
  Short_t GetEnd(){return fEnd;}
  Short_t GetSum(){return fSumAdc;}

private:
  Short_t fBegin;
  Short_t fEnd;//nt-1;//end of decreasing
  Short_t fSumAdc;

  ClassDef(AliPadMax,1)  // Time Projection Chamber Kr clusters
};

#endif


