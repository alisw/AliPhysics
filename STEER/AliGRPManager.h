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

  Bool_t      ReadGRPEntry();
  Bool_t      SetMagField();

  AliRunInfo* GetRunInfo();

private:
  Bool_t SetFieldMap(Float_t l3Current=30000., Float_t diCurrent=6000., 
		     Float_t l3Pol=-1., Float_t dipPol=-1.,
		     Int_t convention=0, Bool_t uniform = kFALSE, 
		     Float_t benergy=7000., const Char_t* btype="pp",
		     const Char_t* path="$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root");
  
  AliGRPObject*  fGRPData;        // Data from the GRP/GRP/Data CDB folder

  AliGRPManager(const AliGRPManager& man);
  AliGRPManager& operator = (const AliGRPManager& man);

  ClassDef(AliGRPManager, 0)      // class for accessing GRP
};

#endif
