#ifndef AliRICHReconstructor_h
#define AliRICHReconstructor_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: */

#include <TTask.h>
#include <TROOT.h>

class AliRICHReconstructor: public TTask 
{
public:
// ctor & dtor:      
   AliRICHReconstructor(){} 
   AliRICHReconstructor(const char *sName,const char *sTitle):TTask(sName,sTitle)
	 {gROOT->GetListOfTasks()->Add(this);gROOT->GetListOfBrowsables()->Add(this);}  
   virtual ~AliRICHReconstructor(){}
// The following staff is defined in .cxx:   
   virtual void  Exec(Option_t *option);
protected:
   ClassDef(AliRICHReconstructor,1)     // Digits or Clusters to RecoClonesArray 
};
#endif // AliRICHReconstructor_h
