#ifndef ALILEGO_H
#define ALILEGO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
//    Utility class to compute and draw Radiation Length Map                 //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TH2.h>

class AliLegoGenerator;

class AliLego : public TNamed  {

public:
  AliLego();
  AliLego(const char *title, Int_t ntheta,Float_t themin,
	  Float_t themax, Int_t nphi, Float_t phimin,
	  Float_t phimax,Float_t rmin,Float_t rmax,Float_t zmax);
  AliLego(const AliLego &lego) {lego.Copy(*this);}
  virtual ~AliLego();
  void  Copy(AliLego &lego) const;
  virtual void  StepManager();
  virtual void  BeginEvent();
  virtual void  FinishEvent();
  virtual void  FinishRun();
  virtual AliLego &operator=(const AliLego &lego) 
  {lego.Copy(*this);return(*this);}
  
private:
  AliLegoGenerator *fGener;     //Lego generator
   Float_t    fTotRadl;         //Total Radiation length
   Float_t    fTotAbso;         //Total absorption length
   Float_t    fTotGcm2;         //Total G/CM2 traversed
   TH2F      *fHistRadl;        //Radiation length map 
   TH2F      *fHistAbso;        //Interaction length map
   TH2F      *fHistGcm2;        //g/cm2 length map
   TH2F      *fHistReta;        //Radiation length map as a function of eta
   
  ClassDef(AliLego,1) //Utility class to compute and draw Radiation Length Map

};


#endif
