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

#include "AliTRDCalOnlineGainTableMCM.h"
#include "AliTRDCalOnlineGainTable.h"

//////////////////////////////////////////////////////////////////////////////////////////////
//
// Data structure to store gaintables of the online calibration in the OCDB
// consisting of three classes:
// AliTRDCalOnlineGainTable 
// AliTRDCalOnlineGainTableROC 
// AliTRDCalOnlineGainTableMCM
//
// AliTRDCalOnlineGainTable is the main class from which all stored data can be accessed.
// The two sub-classes AliTRDCalOnlineGainTableROC and AliTRDCalOnlineGainTableMCM
// contain the gaintables on ROC level and on the MCM level respectively.
//
// The online calibration is used to compensate gain deviations on the pad level.
// For the offline reconstruction the online calibration has to be undone. 
// The corresponding gain correction factor that was used by the online gain filter can be accessed 
// via the functions AliTRDCalOnlineGainTable::GetGainCorrectionFactor(Int_t det, Int_t row, Int_t col) 
// and AliTRDCalOnlineGainTable::GetGainCorrectionFactor(Int_t sector, Int_t stack, Int_t layer, Int_t row, Int_t col).
//
// With the class AliTRDCalOnlineGainTablesMCM all values used for the 
// online calibration can be set and accessed on the MCM/channel level
//
//////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliTRDCalOnlineGainTableMCM);

//_____________________________________________________________________________
AliTRDCalOnlineGainTableMCM::AliTRDCalOnlineGainTableMCM()
  :TObject()
  ,fAdcdac(0)
  ,fMCMGain(0.0)
{
  //
  // constructor
  //

  fAdcdac = -1;
  fMCMGain = AliTRDCalOnlineGainTable::UnDef;
  
  for (int i=0; i<21;i++) {

    fFGFN[i] = -1;
    fFGAN[i] = -1;
  }

}

//_____________________________________________________________________________
AliTRDCalOnlineGainTableMCM::~AliTRDCalOnlineGainTableMCM()
{
  //
  // destructor
  //

}

//_____________________________________________________________________________
Float_t AliTRDCalOnlineGainTableMCM::GetGainCorrectionFactor(Int_t channel)
{
  //
  // returns the Gain Correction Factor of the given channel that was used by the online gain filter
  // 1.0 means no correction
  // 0.9 means a correction of -10%
  // 1.1 means a correction of +10% etc.
  //

  if (fAdcdac == 0){
    if (fFGFN[channel] < 0){
      return -1.;
    }
    else if(fFGFN[channel] > 511){
      return -999.;
    }
    else{
      return (fFGFN[channel]/2048.)+0.875;
 
    }
  }
  else{
    // if the reference voltage of the ADC is not the default value 
    // this taken into account for the Gain Correction Factor
    Float_t fAdcdac_Correction = ( 1./(1.+((Float_t)fAdcdac/31.)*(0.4/1.05))); 
    return (fAdcdac_Correction*((fFGFN[channel]/2048.)+0.875));
  }

}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTableMCM::GetAdcdac()
{
  //
  // returns an integer between 0 and 31 which corresponds to an ADC reference voltage between 1.05V and 1.45V
  // U_Ref =  (1.05V + (fAdcdac/31)*0.4V
  // fAdcdac is the same value for all ADCs within one MCM
  //

  return fAdcdac;

}

//_____________________________________________________________________________
Float_t AliTRDCalOnlineGainTableMCM::GetMCMGain()
{
  //
  // returns the Gain Factor which would lead to a Gain Correction Factor of 1.0
  // this value is the same for all channels within one MCM
  // this value is used for the online PID
  //

  return fMCMGain;

}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTableMCM::GetFGAN(Int_t channel)
{
  //
  // returns the Gain Correction Filter Additive as an interger between 0 and 15 
  // as it is loaded into the TRAP
  //

  if (fFGAN[channel] < 0){
    return -1;
  }
  else if(fFGAN[channel] > 511){
    return -999;
  }
  else{
    return fFGAN[channel];
      
  }

}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTableMCM::GetFGFN(Int_t channel)
{
  //
  // returns the Gain Correction Filter Factors as an interger between 0 and 511 
  // as it is loaded into the TRAP 
  //

  if (fFGFN[channel] < 0){
    return -1;
  }
  else if(fFGFN[channel] > 511){
    return -999;
  }
  else{
    return fFGFN[channel];
      
  }

}

