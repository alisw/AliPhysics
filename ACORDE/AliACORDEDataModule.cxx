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




///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class of  Calibration for ACORDE  Data Modules                            //
// Pedro Gonzalez Zamora    pedro.gonzalez@fcfm.buap.mx                      //
// Irais Bautista Guzman    irais@fcfm.buap.mx                               //
// Arturo Fernandez Tellez afernan@cern.ch                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "AliACORDEDataModule.h"

ClassImp(AliACORDEDataModule)


AliACORDEDataModule::AliACORDEDataModule():
TNamed(),
fRate(0x0),
fStatus(kTRUE)
{
 // Default constructor
}

//________________________________________________________________
AliACORDEDataModule::AliACORDEDataModule(Float_t value,Bool_t status,const char* name):
TNamed(),
fRate(value),
fStatus(status)
{
 
SetName(name);

}

//________________________________________________________________
AliACORDEDataModule::~AliACORDEDataModule()
{
  // Destructor
}
//__________________________________________________________________
Float_t AliACORDEDataModule::GetRate()
{
 return fRate;
}
//_________________________________________________________________
Bool_t AliACORDEDataModule::GetStatus()
{
return fStatus;
}
void AliACORDEDataModule::SetRate(Float_t value)
{
    fRate = value;
}
 
