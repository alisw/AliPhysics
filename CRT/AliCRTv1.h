#ifndef ALICRTV1_H
#define ALICRTV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: CRTv1         //
////////////////////////////////////////////////

#include "AliCRTv0.h"

class AliCRTv1 : public AliCRTv0
{
 
public:
                   AliCRTv1();
                   AliCRTv1(const char *name, const char *title);
                   AliCRTv1(const AliCRTv1& crt);
                   AliCRTv1& operator= (const AliCRTv1& crt);
   virtual         ~AliCRTv1() {}

   virtual void    CreateGeometry();
   virtual void    Init();
   virtual Int_t   IsVersion() const {return 1;}
   virtual void    DrawDetector();
   virtual TString Version(void) {return TString("v1");}
   virtual void    StepManager();
   
   void IncludeRICH(Bool_t status = kTRUE) {fRICHStatus=status;}
   void IncludeMagnet(Bool_t status = kTRUE) {fMagnetStatus=status;}
   void IncludeTPC(Bool_t status = kTRUE) {fTPCStatus=status;}

protected:
   virtual void CreateMolasse();
   virtual void CreateShafts();

   void    CreateRICHGeometry();
   void    CreateTPCGeometry();
   void    CreateMagnetGeometry();

   Bool_t fCRTModule;
   Bool_t fRICHStatus;
   Bool_t fTPCStatus;
   Bool_t fMagnetStatus;
   Bool_t fCRTStatus;

private: 
  ClassDef(AliCRTv1,1)  //Class for CRT, version 1, Shafts outside of AliHALL

};

#endif
