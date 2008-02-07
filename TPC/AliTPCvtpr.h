#ifndef ALITPCVTPR_H
#define ALITPCVTPR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCclusterKr.h,v 1.8 2007/12/31 16:07:15 matyja Exp $ */

//-------------------------------------------------------
//                    TPC cordinate Class
//
//   Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl
//-------------------------------------------------------

#include "TObject.h"

//_____________________________________________________________________________
class AliTPCvtpr: public TObject{
public:
  AliTPCvtpr();
  ~AliTPCvtpr();

  AliTPCvtpr(Short_t max,Short_t nt,Short_t np,Short_t nr);
  AliTPCvtpr(const AliTPCvtpr & param);
  AliTPCvtpr & operator = (const AliTPCvtpr & param);


  void SetAdc(Short_t q){fAdc=q;}
  void SetTime(Short_t q){fTime=q;}
  void SetPad(Short_t q){fPad=q;}
  void SetRow(Short_t q){fRow=q;}


  Short_t GetAdc(){return fAdc;}
  Short_t GetTime(){return fTime;}
  Short_t GetPad(){return fPad;}
  Short_t GetRow(){return fRow;}


protected:
  Short_t fAdc;
  Short_t fTime;
  Short_t fPad;
  Short_t fRow;

private:

  ClassDef(AliTPCvtpr,1)  // TPC coordinate class
};

#endif


