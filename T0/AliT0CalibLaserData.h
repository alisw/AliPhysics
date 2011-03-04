#ifndef AliT0CalibLaserData_H
#define AliT0CalibLaserData_H
/***************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Alla Maevskaya INR RAS alla@inr.ru
 *
 * See cxx source for full Copyright notice                               
 ***************************************************************************/
#include "TGNumberEntry.h"
#include "TGTextEntry.h"
#include "TObject.h"

class AliT0CalibLaserData : public TObject
{
 public:
  AliT0CalibLaserData ();
  AliT0CalibLaserData(const AliT0CalibLaserData &calibda) : TObject(calibda),
    fTEntry(0),
    fFileName(" ")
    {
      for ( Int_t i=0; i<30; i++ )
      {
         fEntries[i] = NULL;
         fHistLimits[i] = 0;
      }
    }
  AliT0CalibLaserData & operator= (const AliT0CalibLaserData  &) {return *this;}
  virtual ~AliT0CalibLaserData() {}
  void ReadHistSize();
  void DoOk();
  void            ReadData();
  void OpenFile();

 private:
  TGNumberEntry * fEntries[30];     //for GUI histogram limits
  TGTextEntry * fTEntry;     //for GUI file name
  double          fHistLimits[30];  // histogram limits
  const char *fFileName;     // file name


  ClassDef(AliT0CalibLaserData,2)
};

#endif
