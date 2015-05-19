#ifndef ALIALGSENSTOF_H
#define ALIALGSENSTOF_H

#include "AliAlgSens.h"


class TObjArray;


/*--------------------------------------------------------
  TOF sensor
  -------------------------------------------------------*/

// Author: ruben.shahoyan@cern.ch


class AliAlgSensTOF : public AliAlgSens
{
 public:
  AliAlgSensTOF(const char* name=0, Int_t vid=0, Int_t iid=0, Int_t isec=0);
  virtual ~AliAlgSensTOF();
  //
  virtual void   SetTrackingFrame();
  //
  Int_t GetSector()                      const {return fSector;}
  void  SetSector(UInt_t sc)                   {fSector = (UChar_t)sc;}
  //
 protected:
  //
  UChar_t fSector;                      // sector ID
  //
  ClassDef(AliAlgSensTOF,1)
};


#endif
