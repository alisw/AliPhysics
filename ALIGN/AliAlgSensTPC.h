#ifndef ALIALGSENSTPC_H
#define ALIALGSENSTPC_H

#include "AliAlgSens.h"


class TObjArray;


/*--------------------------------------------------------
  TPC sensor (chamber)
  -------------------------------------------------------*/

// Author: ruben.shahoyan@cern.ch


class AliAlgSensTPC : public AliAlgSens
{
 public:
  AliAlgSensTPC(const char* name=0, Int_t vid=0, Int_t iid=0, Int_t isec=0);
  virtual ~AliAlgSensTPC();
  //
  Int_t GetSector()                      const {return fSector;}
  void  SetSector(UInt_t sc)                   {fSector = (UChar_t)sc;}
  //
  virtual void   SetTrackingFrame();
  //
 protected:
  //
  UChar_t fSector;                      // sector ID

  ClassDef(AliAlgSensTPC,1)
};


#endif
