#ifndef AliT0CalibLaserData_H
#define AliT0CalibLaserData_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TGNumberEntry.h"
#include "TObject.h"

class AliT0CalibLaserData : public TObject
{
 public:
  AliT0CalibLaserData ();
  AliT0CalibLaserData(const AliT0CalibLaserData &calibda) : TObject(calibda),
  	           fRunNumber(0) {}
  AliT0CalibLaserData & operator= (const AliT0CalibLaserData  &) {return *this;}
  virtual ~AliT0CalibLaserData() {}
  void ReadHistSize(Int_t rNumber=905);
  void DoOk();
  void            ReadData();

 private:
 // void            ReadData();
  Int_t           fRunNumber;
  TGNumberEntry * fEntries[30];
  double          fHistLimits[30];

  ClassDef(AliT0CalibLaserData,1)
};

#endif
