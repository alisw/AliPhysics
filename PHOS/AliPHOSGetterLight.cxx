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

//_________________________________________________________________________
// 
//
//*-- Author :  D.Peressounko (RRC KI) 
//////////////////////////////////////////////////////////////////////////////
//Class designed to perform fast on-flight reconstruction of raw data
//without writing it to file and without creation of all aliroot folder 
//structure. This class keeps lists of (raw) digits, RecPoints TrackSegments and 
//RecParticles for one (current) event. In addition, for convinience, it returns 
//pointers to PHOS geometry, Clusterizser, TSMaker and PID maker.

// --- ROOT system ---
// --- Standard library ---

// --- AliRoot header files ---
//#include "AliPHOSCalibrationDB.h"
#include "AliPHOSGetterLight.h"

ClassImp(AliPHOSGetterLight)


//____________________________________________________________________________ 
AliPHOSGetterLight::AliPHOSGetterLight() : 
  AliPHOSGetter(0),
  fDigits(0),
  fEmcRecPoints(0),
  fCpvRecPoints(0),
  fTS(0),
  fRP(0),
  fClusterizer(0), 
  fTSM(0),
  fPID(0)
{
  // ctor
  SetRawDigits(kTRUE) ;
  //  fRawDigits =kTRUE;
  fgObjGetter = this ;
}
//____________________________________________________________________________ 
AliPHOSGetterLight::AliPHOSGetterLight(const char* /*alirunFileName*/, const char* /*version*/, Option_t * /*openingOption*/) :
  AliPHOSGetter(0),
  fDigits(new TClonesArray("AliPHOSDigit",256)),
  fEmcRecPoints(new TObjArray(50)),
  fCpvRecPoints(new TObjArray(0)),
  fTS(new TClonesArray("AliPHOSTrackSegment",50)),
  fRP(new TClonesArray("AliPHOSRecParticle",50)),
  fClusterizer (0),
  fTSM(0),
  fPID(0)
{
  //Create containers of reconstructed objects for one event
  fEmcRecPoints->SetOwner(kTRUE) ;
  fCpvRecPoints->SetOwner(kTRUE) ;

  SetRawDigits(kTRUE) ;
  //  fRawDigits = kTRUE ;
  fgObjGetter = this ;
}

//____________________________________________________________________________ 
AliPHOSGetterLight::AliPHOSGetterLight(const AliPHOSGetterLight & obj) : 
  AliPHOSGetter(obj),
  fDigits(0),
  fEmcRecPoints(0),
  fCpvRecPoints(0),
  fTS(0),
  fRP(0),
  fClusterizer(0), 
  fTSM(0),
  fPID(0)
{
  // cpy ctor requested by Coding Convention 
  Fatal("cpy ctor", "not implemented") ;
} 

//____________________________________________________________________________ 
  AliPHOSGetterLight::~AliPHOSGetterLight()
{
  //Delete containers owned by Getter and do not touch other pointers
  if(fDigits){ delete fDigits ; fDigits = 0 ;}
  if(fEmcRecPoints){ delete fEmcRecPoints; fEmcRecPoints = 0 ;}
  if(fCpvRecPoints){ delete fCpvRecPoints; fCpvRecPoints = 0 ;}
  if(fTS){ delete fTS; fTS = 0 ;}
  if(fRP){ delete fRP; fRP = 0 ;}
}
//____________________________________________________________________________ 
AliPHOSGetterLight * AliPHOSGetterLight::Instance(const char* alirunFileName, const char* version, Option_t * openingOption) 
{
  // Creates and returns the pointer of the unique instance
  // Must be called only when the environment has changed
  
  if(!fgObjGetter){ // first time the getter is called 
    fgObjGetter = (AliPHOSGetter*) new AliPHOSGetterLight(alirunFileName, version, openingOption) ;
  }
  return (AliPHOSGetterLight*) fgObjGetter ;
}
//____________________________________________________________________________ 
AliPHOSGetterLight * AliPHOSGetterLight::Instance(void) 
{
  return (AliPHOSGetterLight*) fgObjGetter ;
}
