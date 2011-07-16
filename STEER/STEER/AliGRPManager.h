#ifndef ALIGRPMANAGER_H
#define ALIGRPMANAGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliGRPManager class                                                    //
// The class can be used in order to access and read the Global Run       //
// Parameters entry from OCDB.                                            //
// It has a methods to set the magnetic field instanton and return        //
// the run and event info objects.                                        //
//                                                                        //
// cvetan.cheshkov@cern.ch 15/06/2009                                     //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRunInfo;
class AliGRPObject;

class AliGRPManager: public TObject {
public:
  AliGRPManager();
  virtual ~AliGRPManager();

  const AliGRPObject* GetGRPData() const { return fGRPData; }

  void        SetGRPEntry(AliGRPObject* source);
  Bool_t      ReadGRPEntry();
  Bool_t      SetMagField();

  AliRunInfo* GetRunInfo();

private:
  
  AliGRPObject*  fGRPData;        // Data from the GRP/GRP/Data CDB folder

  AliGRPManager(const AliGRPManager& man);
  AliGRPManager& operator = (const AliGRPManager& man);

  ClassDef(AliGRPManager, 0)      // class for accessing GRP
};

#endif
