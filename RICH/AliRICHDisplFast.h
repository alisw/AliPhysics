#ifndef AliRICHDisplFast_h
#define AliRICHDisplFast_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TTask.h>

class AliRunLoader;

class AliRICHDisplFast : public TTask 
{
public :
               AliRICHDisplFast() {;}
  virtual     ~AliRICHDisplFast() {;}      
  void ShowEvent(Int_t iEvtNmin,Int_t iEvtNmax);                                 //Display looping on events
  void ShowEvent(Int_t iEvent)                      {ShowEvent(iEvent,iEvent);}  //Display only one event
  static Int_t    Nparticles   (Int_t iPID,Int_t iEvent=0,AliRunLoader *p=0);    //returns a number of particles of given type
  virtual void Exec(Option_t *opt=0);                                            //virtual do the main job
protected:  
  ClassDef(AliRICHDisplFast,0)                              //Utility class to draw the current event topology
};
    
#endif //AliRICHDisplFast_cxx

