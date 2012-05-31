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

#include <AliTRDfeeParam.h>
#include "AliTRDCalOnlineGainTableROC.h"
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
// AliTRDCalOnlineGainTableROC is a class to allocate MCM Gain Tables 
// and to access all stored calibration values from the ROC level by indicating row and col
//
//////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliTRDCalOnlineGainTableROC);

//_____________________________________________________________________________
AliTRDCalOnlineGainTableROC::AliTRDCalOnlineGainTableROC()
  :TObject()
{
  //
  // constructor
  //

  for (int i=0; i<128; i++) {
    fMCMGainTables[i] = 0;
  }
}

//_____________________________________________________________________________
AliTRDCalOnlineGainTableROC::AliTRDCalOnlineGainTableROC(const AliTRDCalOnlineGainTableROC& other)
  :TObject(other)
{
  //
  // copy constructor
  //

  for (int i=0; i<128; i++) {
    if (other.GetGainTableMCM(i)) {
      fMCMGainTables[i] = new AliTRDCalOnlineGainTableMCM( *(other.GetGainTableMCM(i)) );
    } else {
      fMCMGainTables[i] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTRDCalOnlineGainTableROC& AliTRDCalOnlineGainTableROC::operator=(const AliTRDCalOnlineGainTableROC& other)
{
  //
  // assignment operator
  //

  for (int i=0; i<128; i++) {

    if (fMCMGainTables[i]) {
      delete fMCMGainTables[i];
    }

    if (other.GetGainTableMCM(i)) {
      fMCMGainTables[i] = new AliTRDCalOnlineGainTableMCM( *(other.GetGainTableMCM(i)) );
    } else {
      fMCMGainTables[i] = 0;
    }
  }

  return *this;

}

//_____________________________________________________________________________
AliTRDCalOnlineGainTableROC::~AliTRDCalOnlineGainTableROC()
{
  //
  // destructor
  //

  for (int i=0; i<128; i++) {
    if (fMCMGainTables[i]) {
      delete fMCMGainTables[i];
    }
  }

}

//_____________________________________________________________________________
Float_t AliTRDCalOnlineGainTableROC::GetGainCorrectionFactor(Int_t row, Int_t col)
{
  //
  // chooses ROB/MCM/channel from row/col
  // returns the Gain Correction Factor of the given channel
  //

  AliTRDfeeParam * para = AliTRDfeeParam::Instance(); 

  Int_t rob = para->GetROBfromPad(row,col);
  Int_t mcm = para->GetMCMfromPad(row,col);
  Int_t channel =19-(col%18);

  AliTRDCalOnlineGainTableMCM* gtbl=GetGainTableMCM(rob,mcm);

  if (gtbl) {
    return gtbl->GetGainCorrectionFactor(channel);
  } else {
    return AliTRDCalOnlineGainTable::UnDef;
  }

}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTableROC::GetAdcdac(Int_t row, Int_t col)
{
  //
  // chooses ROB/MCM/channel from row/col
  // returns the ADC's reference voltage of the given MCM
  //

  AliTRDfeeParam * para = AliTRDfeeParam::Instance(); 

  Int_t rob = para->GetROBfromPad(row,col);
  Int_t mcm = para->GetMCMfromPad(row,col);

  AliTRDCalOnlineGainTableMCM* gtbl=GetGainTableMCM(rob,mcm);

  if (gtbl) {
    return gtbl->GetAdcdac();
  } else {
    return -999;
  }

}

//_____________________________________________________________________________
Float_t AliTRDCalOnlineGainTableROC::GetMCMGain(Int_t row, Int_t col)
{
  //
  // chooses ROB/MCM/channel from row/col
  // returns the Gain Factor which would lead to a Correction Factor of 1.0 within the given MCM
  //

  AliTRDfeeParam * para = AliTRDfeeParam::Instance(); 

  Int_t rob = para->GetROBfromPad(row,col);
  Int_t mcm = para->GetMCMfromPad(row,col);

  AliTRDCalOnlineGainTableMCM* gtbl=GetGainTableMCM(rob,mcm);

  if (gtbl) {
    return gtbl->GetMCMGain();
  } else {
    return AliTRDCalOnlineGainTable::UnDef;
  }

}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTableROC::GetFGAN(Int_t row, Int_t col)
{
  //
  // chooses ROB/MCM/channel from row/col
  // returns the Gain Correction Filter Factor of the given channel
  //

  AliTRDfeeParam * para = AliTRDfeeParam::Instance(); 

  Int_t rob = para->GetROBfromPad(row,col);
  Int_t mcm = para->GetMCMfromPad(row,col);
  Int_t channel =19-(col%18);

  AliTRDCalOnlineGainTableMCM* gtbl=GetGainTableMCM(rob,mcm);

  if (gtbl) {
    return gtbl->GetFGAN(channel);
  } else {
    return -999;
  }

}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTableROC::GetFGFN(Int_t row, Int_t col)
{
  //
  // chooses ROB/MCM/channel from row/col
  // returns the Gain Correction Filter Additive of the given channel
  //

  AliTRDfeeParam * para = AliTRDfeeParam::Instance(); 

  Int_t rob = para->GetROBfromPad(row,col);
  Int_t mcm = para->GetMCMfromPad(row,col);
  Int_t channel =19-(col%18);

  AliTRDCalOnlineGainTableMCM* gtbl=GetGainTableMCM(rob,mcm);

  if (gtbl) {
    return gtbl->GetFGFN(channel);
  } else {
    return -999;
  }

}  

//_____________________________________________________________________________
void AliTRDCalOnlineGainTableROC::AllocateGainTableMCM(Int_t rob, Int_t mcm)
{
  //
  // allocates a Gain Table for the given MCM
  //

  Int_t index = rob*16 + mcm;

  if (fMCMGainTables[index]) {
    delete fMCMGainTables[index];
  }

  fMCMGainTables[index] = new AliTRDCalOnlineGainTableMCM;
}

