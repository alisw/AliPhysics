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
// Short description  
//
//*-- Author :  D.Peressounko (RRC KI & SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"  
#include "AliPHOSCalibrManager.h"
#include "AliPHOSCalibrationDB.h"
ClassImp(AliPHOSCalibrationDB)


//____________________________________________________________________________ 
  AliPHOSCalibrationDB::AliPHOSCalibrationDB():fPedestals(), fGains()
{
}
//____________________________________________________________________________ 
AliPHOSCalibrationDB::AliPHOSCalibrationDB(const char * database):
  TNamed("AliPHOSCalibrationDB",database),fPedestals("Pedestals"),fGains("Gains"){

} 
//____________________________________________________________________________ 
  AliPHOSCalibrationDB::~AliPHOSCalibrationDB()
{
}
//____________________________________________________________________________ 
Float_t AliPHOSCalibrationDB::Calibrate(Int_t amp, Int_t absId)const 
{
  //if absID is known, return calibrated energy, else - zero
  Float_t ret = (amp - fPedestals.Data(absId))*fGains.Data(absId); 
  if(ret > 0)
    return ret ;
  else
    return 0.0000001 ; //Should not be zero - to avoid FPE
}

//____________________________________________________________________________
void AliPHOSCalibrationDB::GetParameters(void){
  //In this method we read calibration parameters using AliPHOSCalibrManager
  //This manager should be configured to read from correct source.

  AliPHOSCalibrManager * man = AliPHOSCalibrManager::GetInstance() ;
  if(!man){
    Error("GetParameters","AliPHOSCalibrManager not instanciated") ;
    return ;
  }


  man->GetParameters(fPedestals) ;
  man->GetParameters(fGains) ;
}

//____________________________________________________________________________
AliPHOSCalibrationDB& AliPHOSCalibrationDB::operator=(AliPHOSCalibrationDB const & cdb)
{
  //
  fPedestals=cdb.fPedestals ;
  fGains = cdb.fGains;
  return *this ;
}
