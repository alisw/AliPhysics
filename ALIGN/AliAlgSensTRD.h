#ifndef ALIALGSENSTRD_H
#define ALIALGSENSTRD_H

#include "AliAlgSens.h"


class TObjArray;


/*--------------------------------------------------------
  TRD sensor
  -------------------------------------------------------*/

// Author: ruben.shahoyan@cern.ch


class AliAlgSensTRD : public AliAlgSens
{
 public:
  AliAlgSensTRD(const char* name=0, Int_t vid=0, Int_t iid=0, Int_t isec=0);
  virtual ~AliAlgSensTRD();
  //
  Int_t GetSector()                      const {return fSector;}
  void  SetSector(UInt_t sc)                   {fSector = (UChar_t)sc;}
  //
  virtual void   SetTrackingFrame();
  //
 protected:
  //
  UChar_t fSector;                      // sector ID

  ClassDef(AliAlgSensTRD,1)
};


#endif
