// -*- mode: C++ -*- 
#ifndef ALIESDCALOTRIGGER_H
#define ALIESDCALOTRIGGER_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//                          Class AliESDCaloTrigger
//   This is a class that summarizes the Trigger Data of EMCal and Phos
//   for the ESD   
//   Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------



#include <TNamed.h>
#include <TArrayF.h>




class AliESDCaloTrigger : public TNamed {
public:
  AliESDCaloTrigger();
  AliESDCaloTrigger(const  AliESDCaloTrigger& ctrig);
  AliESDCaloTrigger& operator=(const  AliESDCaloTrigger& ctrig);
  virtual ~AliESDCaloTrigger();

  // does this create mem leak? CKB use new with placement?
  void AddTriggerPosition(const TArrayF & array)  { 
    if(fTriggerPosition) delete fTriggerPosition;
    fTriggerPosition =  new TArrayF(array);
  }

  void AddTriggerAmplitudes(const TArrayF & array) { 
    if(fTriggerAmplitudes)delete fTriggerAmplitudes;
    fTriggerAmplitudes  = new TArrayF(array); 
  }
  
  void Reset(); 

  TArrayF* GetTriggerPosition()    {return fTriggerPosition;}
  TArrayF* GetTriggerAmplitudes()  {return fTriggerPosition;}
  

private:

  TArrayF *fTriggerAmplitudes; // Amplitude of PHOS or EMCal Trigger
  TArrayF *fTriggerPosition;   // Position of PHOS or EMCal Trigger

  ClassDef(AliESDCaloTrigger,1)
};


#endif

