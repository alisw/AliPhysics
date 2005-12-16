/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calibration parameters by accessing the CDB           //
//                                                                           //
// Request an instance with AliTRDcalibDB::Instance()                 //
// If a new event is processed set the event number with SetRun              //
// Then request the calibration data                                         // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TRandom.h>

#include <AliCDBManager.h>

#include "AliTRDcalibDB.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDCommonParam.h"

#include "AliTRDCalROC.h"
#include "AliTRDCalChamber.h"
#include "AliTRDCalStack.h"
#include "AliTRDCalPad.h"
#include "AliTRDCalDet.h"
#include "AliTRDCalGlobals.h"

ClassImp(AliTRDcalibDB)

AliTRDcalibDB* AliTRDcalibDB::fgInstance = 0;
Bool_t AliTRDcalibDB::fgTerminated = kFALSE;

//_ singleton implementation __________________________________________________
AliTRDcalibDB* AliTRDcalibDB::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  // 
  
  if (fgTerminated != kFALSE)
    return 0;

  if (fgInstance == 0)
    fgInstance = new AliTRDcalibDB();
  
  return fgInstance;
}

void AliTRDcalibDB::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class and sets the terminated flag, instances cannot be requested anymore
  // This function can be called several times.
  //
  
  fgTerminated = kTRUE;
  
  if (fgInstance != 0)
  {
    delete fgInstance;
    fgInstance = 0;
  }
}

//_____________________________________________________________________________
AliTRDcalibDB::AliTRDcalibDB()
{
  //
  // constructor
  //
  
  fLocator = 0;
  // TODO Default runnumber is set to 0, this should be changed later to an invalid value (e.g. -1) to prevent
  // TODO invalid calibration data to be used.
  fRun = 0;
  
  AliCDBManager* manager = AliCDBManager::Instance();
  if (!manager)
  {
    std::cout << "Failed to get instance of AliCDBManager." << std::endl;
    return;
  }
  
  fLocator = manager->GetStorage("local://$ALICE_ROOT");
  
  for (Int_t i=0; i<kCDBCacheSize; ++i)
  {
    fCDBCache[i] = 0;
    fCDBEntries[i] = 0;
  }
};

//_____________________________________________________________________________
AliTRDcalibDB::~AliTRDcalibDB() 
{
  //
  // destructor
  //
  
  Invalidate();
};

//_____________________________________________________________________________
void AliTRDcalibDB::SetRun(Long64_t run)
{
  //
  // Sets current run number. Calibration data is read from the corresponding file. 
  // When the run number changes the caching is invalidated.
  //
  
  if (fRun == run)
    return;
  
  fRun = run;
  Invalidate();
}
  
//_____________________________________________________________________________
void AliTRDcalibDB::Invalidate()
{
  //
  // Invalidates cache (when run number is changed).
  //
  
  for (Int_t i=0; i<kCDBCacheSize; ++i)
  {
    if (fCDBEntries[i])
    {
      if (fCDBEntries[i]->IsOwner() != kFALSE && fCDBCache[i])
        delete fCDBCache[i];
      delete fCDBEntries[i];
      fCDBEntries[i] = 0;
      fCDBCache[i] = 0;
    }
  }
}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::GetChamberPos(Int_t det, Float_t* xyz)
{
  //
  // Returns the deviation of the chamber position from the nominal position.
  //
  
  AliTRDCalChamber* chamber = dynamic_cast<AliTRDCalChamber*>(GetCachedCDBObject(kIDChamber));
  if (!chamber)
    return kFALSE;
  
  const Float_t* kvalues = chamber->GetChamberPos(det);
  if (!kvalues)
    return kFALSE;
  
  xyz[0] = kvalues[0];
  xyz[1] = kvalues[1];
  xyz[2] = kvalues[2];
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::GetChamberRot(Int_t det, Float_t* xyz)
{
  //
  // Returns the rotation of the chamber from the nominal position.
  //
  
  AliTRDCalChamber* chamber = dynamic_cast<AliTRDCalChamber*>(GetCachedCDBObject(kIDChamber));
  if (!chamber)
    return kFALSE;
  
  const Float_t* kvalues = chamber->GetChamberRot(det);
  if (!kvalues)
    return kFALSE;
  
  xyz[0] = kvalues[0];
  xyz[1] = kvalues[1];
  xyz[2] = kvalues[2];
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::GetStackPos(Int_t chamber, Int_t sector, Float_t* xyz)
{
  //
  // Returns the deviation of the stack position from the nominal position.
  //
  
  AliTRDCalStack* stack = dynamic_cast<AliTRDCalStack*>(GetCachedCDBObject(kIDStack));
  if (!stack)
    return kFALSE;
  
  const Float_t* kvalues = stack->GetStackPos(chamber, sector);
  if (!kvalues)
    return kFALSE;
  
  xyz[0] = kvalues[0];
  xyz[1] = kvalues[1];
  xyz[2] = kvalues[2];
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTRDcalibDB::GetStackRot(Int_t chamber, Int_t sector, Float_t* xyz)
{
  //
  // Returns the rotation of the stack from the nominal position.
  //
  
  AliTRDCalStack* stack = dynamic_cast<AliTRDCalStack*>(GetCachedCDBObject(kIDStack));
  if (!stack)
    return kFALSE;
  
  const Float_t* kvalues = stack->GetStackRot(chamber, sector);
  if (!kvalues)
    return kFALSE;
  
  xyz[0] = kvalues[0];
  xyz[1] = kvalues[1];
  xyz[2] = kvalues[2];
  
  return kTRUE;
}

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetVdrift(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns the drift velocity for the given pad.
  //
  
  AliTRDCalPad* calPad = dynamic_cast<AliTRDCalPad*> (GetCachedCDBObject(kIDVdrift));
  if (!calPad)
    return -1;

  AliTRDCalROC* roc = calPad->GetCalROC(det);
  if (!roc)
    return -1;

  return roc->GetValue(col, row);
};

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetT0(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns t0 for the given pad.
  //
  
  AliTRDCalPad* calPad = dynamic_cast<AliTRDCalPad*> (GetCachedCDBObject(kIDT0));
  if (!calPad)
    return -1;

  AliTRDCalROC* roc = calPad->GetCalROC(det);
  if (!roc)
    return -1;

  return roc->GetValue(col, row);
};

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetGainFactor(Int_t det, Int_t col, Int_t row)
{
  //
  // Returns the gain factor for the given pad.
  //
  
  AliTRDCalPad* calPad = dynamic_cast<AliTRDCalPad*> (GetCachedCDBObject(kIDGainFactor));
  if (!calPad)
    return -1;

  AliTRDCalROC* roc = calPad->GetCalROC(det);
  if (!roc)
    return -1;

  return roc->GetValue(col, row);
};

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetSamplingFrequency()
{
  //
  // Returns the sampling frequency of the TRD read-out.
  //
  
  AliTRDCalGlobals* calGlobal = dynamic_cast<AliTRDCalGlobals*> (GetCachedCDBObject(kIDGlobals));
  if (!calGlobal)
    return -1;  
  
  return calGlobal->GetSamplingFrequency();
}
  
//_____________________________________________________________________________
Int_t AliTRDcalibDB::GetNumberOfTimeBins()
{
  //
  // Returns the number of time bins which are read-out.
  //
  
  AliTRDCalGlobals* calGlobal = dynamic_cast<AliTRDCalGlobals*> (GetCachedCDBObject(kIDGlobals));
  if (!calGlobal)
    return -1;  
  
  return calGlobal->GetNumberOfTimeBins();
}

//_____________________________________________________________________________
Float_t AliTRDcalibDB::GetOmegaTau(Float_t vdrift)
{
  //
  // Returns omega*tau (tan(Lorentz-angle)) for a given drift velocity <vd> 
  // and a B-field <b> for Xe/CO2 (15%).
  // The values are according to a GARFIELD simulation.
  //
  // This function basically does not belong to the calibration class. It should be moved somewhere else. 
  // However, currently it is in use by simulation and reconstruction.
  //
  
  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam)
    return -1;
  Float_t field = commonParam->GetField();

  const Int_t kNb = 5;
  Float_t p0[kNb] = {  0.004810,  0.007412,  0.010252,  0.013409,  0.016888 };
  Float_t p1[kNb] = {  0.054875,  0.081534,  0.107333,  0.131983,  0.155455 };
  Float_t p2[kNb] = { -0.008682, -0.012896, -0.016987, -0.020880, -0.024623 };
  Float_t p3[kNb] = {  0.000155,  0.000238,  0.000330,  0.000428,  0.000541 };

  Int_t ib = ((Int_t) (10 * (field - 0.15)));
  ib       = TMath::Max(  0,ib);
  ib       = TMath::Min(kNb,ib);

  Float_t alphaL = p0[ib] 
      + p1[ib] * vdrift
      + p2[ib] * vdrift*vdrift
      + p3[ib] * vdrift*vdrift*vdrift;

  return TMath::Tan(alphaL);
}

