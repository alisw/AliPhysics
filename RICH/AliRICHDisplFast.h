#ifndef AliRICHDisplFast_h
#define AliRICHDisplFast_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TTask.h>
#include "AliRICH.h"
#include <AliRun.h>
#include "TH2.h"
class AliRICH;

class AliRICHDisplFast : public TTask 
{
public :
               AliRICHDisplFast() {;}
  virtual     ~AliRICHDisplFast() {;}      
  static  void DrawSectors();                               //Draw sectors in plot 
  void ShowEvent(Int_t iEvtNmin,Int_t iEvtNmax);            // Display looping on events
  void ShowEvent(Int_t iEvent) {ShowEvent(iEvent,iEvent);}  // Display only one event
  virtual void Exec(Option_t *opt=0);                       //virtual do the main job
protected:  
  ClassDef(AliRICHDisplFast,0)                              //Utility class to draw the current event topology
};
    
#endif // #ifdef AliRICHDisplFast_cxx

