#ifndef AliGenRadioactive_h
#define AliGenRadioactive_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Simple radioactive source generator.
// Currently Sr90 only

#include <AliGenerator.h>
#include <TH1.h>

typedef enum {kSr90} RadioNuclei_t;     

class AliGenRadioactive : public AliGenerator
{
 public:
           AliGenRadioactive() : AliGenerator(),fPartId(0),fGenH1(0) {;}//default ctor
           AliGenRadioactive(Int_t iSource,Int_t iNparts);                  //main ctor
  virtual ~AliGenRadioactive()                                           {if(fGenH1) delete fGenH1;}//dtor
  virtual void Generate();                                                              //interface from AliGenerator, generates current event 


protected:
  Int_t    fPartId;          //ID of secondary from radioactive nucleus
  TH1F    *fGenH1;           //Histogram containg exp spectrum to be used to sample the energy
  ClassDef(AliGenRadioactive,1) // Radioactive source generator
};
#endif
