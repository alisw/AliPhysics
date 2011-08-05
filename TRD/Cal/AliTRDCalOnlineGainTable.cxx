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
//////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliTRDCalOnlineGainTable);

const Float_t AliTRDCalOnlineGainTable::UnDef=-999.;

//_____________________________________________________________________________
AliTRDCalOnlineGainTable::AliTRDCalOnlineGainTable()
  :TObject()
{
  //
  // constructor
  //

  for (int i=0; i<540; i++) {
    fROCGainTables[i] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalOnlineGainTable::AliTRDCalOnlineGainTable(const AliTRDCalOnlineGainTable& other)
  :TObject(other)
{
  //
  // copy constructor
  //

  for (int i=0; i<540; i++) {
    if (other.GetGainTableROC(i)) {
      fROCGainTables[i] = new AliTRDCalOnlineGainTableROC( *(other.GetGainTableROC(i)) );
    } else {
      fROCGainTables[i] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTRDCalOnlineGainTable& AliTRDCalOnlineGainTable::operator=(const AliTRDCalOnlineGainTable& other)
{
  //
  // assignment operator
  //

  for (int i=0; i<540; i++) {

    if (fROCGainTables[i]) {
      delete fROCGainTables[i];
    }

    if (other.GetGainTableROC(i)) {
      fROCGainTables[i] = new AliTRDCalOnlineGainTableROC( *(other.GetGainTableROC(i)) );
    } else {
      fROCGainTables[i] = 0;
    }
  }

  return *this;

}

//_____________________________________________________________________________
AliTRDCalOnlineGainTable::~AliTRDCalOnlineGainTable()
{
  //
  // destructor
  //

  for (int i=0; i<540; i++) {
    if (fROCGainTables[i]) {
      delete fROCGainTables[i];
    }
  }

}

//_____________________________________________________________________________
Float_t AliTRDCalOnlineGainTable::GetGainCorrectionFactor(Int_t det
                                                        , Int_t row
                                                        , Int_t col) const
{
  //
  // returns the Gain Correction Factor of the channel
  // given by det, row, col
  //

  AliTRDCalOnlineGainTableROC* gtbl = GetGainTableROC(det);

  if (gtbl) {
    return gtbl->GetGainCorrectionFactor(row,col);
  } else {
    return AliTRDCalOnlineGainTable::UnDef;
  }

}

//_____________________________________________________________________________
Float_t AliTRDCalOnlineGainTable::GetGainCorrectionFactor(Int_t sector
                                                        , Int_t stack
                                                        , Int_t layer
							, Int_t row
                                                        , Int_t col) const
{ 
  //
  // returns the Gain Correction Factor of the channel
  // given by sector, stack, layer, row, col
  //

  return GetGainCorrectionFactor(30*sector + 6*stack + layer, row, col);
}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTable::GetAdcdac(Int_t det, Int_t row, Int_t col)
{
  //
  // returns the ADC's reference voltage of the channel 
  // given by det, row, col
  //

  AliTRDCalOnlineGainTableROC* gtbl = GetGainTableROC(det);

  if (gtbl) {
    return gtbl->GetAdcdac(row,col);
  } else {
    return -999;
  }

}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTable::GetAdcdac(Int_t sector, Int_t stack, Int_t layer, 
						Int_t row, Int_t col)
{ 
  //
  // returns the ADC's reference voltage of the channel 
  // given by sector, stack, layer, row, col
  //

  return GetAdcdac(30*sector + 6*stack + layer, row, col);
}

//_____________________________________________________________________________
Float_t AliTRDCalOnlineGainTable::GetMCMGain(Int_t det, Int_t row, Int_t col)
{ 
  //
  // returns the Gain Factor which would lead to a Correction Factor of 1.0  
  // within the MCM given by det, row, col
  //

  AliTRDCalOnlineGainTableROC* gtbl = GetGainTableROC(det);

  if (gtbl) {
    return gtbl->GetMCMGain(row,col);
  } else {
    return AliTRDCalOnlineGainTable::UnDef;
  }

}

//_____________________________________________________________________________
Float_t AliTRDCalOnlineGainTable::GetMCMGain(Int_t sector, Int_t stack, Int_t layer, 
						Int_t row, Int_t col)
{  
  //
  // returns the Gain Factor which would lead to a Correction Factor of 1.0  
  // within the MCM given by sector, stack, layer, row, col
  //

  return GetMCMGain(30*sector + 6*stack + layer, row, col);
}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTable::GetFGAN(Int_t det, Int_t row, Int_t col)
{
  //
  // returns the Gain Correction Filter Additive of the channel
  // given by det, row, col
  //

  AliTRDCalOnlineGainTableROC* gtbl = GetGainTableROC(det);

  if (gtbl) {
    return gtbl->GetFGAN(row,col);
  } else {
    return -999;
  }

}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTable::GetFGAN(Int_t sector, Int_t stack, Int_t layer, 
						Int_t row, Int_t col)
{ 
  //
  // returns the Gain Correction Filter Additive of the channel
  // given by sector, stack, layer, row, col
  //

  return GetFGAN(30*sector + 6*stack + layer, row, col);
}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTable::GetFGFN(Int_t det, Int_t row, Int_t col)
{  
  //
  // returns the Gain Correction Filter Factor of the channel
  // given by det, row, col
  //

  AliTRDCalOnlineGainTableROC* gtbl = GetGainTableROC(det);

  if (gtbl) {
    return gtbl->GetFGFN(row,col);
  } else {
    return -999;
  }
}

//_____________________________________________________________________________
Short_t AliTRDCalOnlineGainTable::GetFGFN(Int_t sector, Int_t stack, Int_t layer, 
						Int_t row, Int_t col)
{  
  //
  // returns the Gain Correction Filter Factor of the channel
  // given by sector, stack, layer, row, col
  //

  return GetFGFN(30*sector + 6*stack + layer, row, col);
}

//_____________________________________________________________________________
void AliTRDCalOnlineGainTable::AllocateGainTableROC(Int_t det)
{ 
  //
  // allocates a Gain Table for the given detector
  //

  if (fROCGainTables[det]) {
    delete fROCGainTables[det];
  }

  fROCGainTables[det] = new AliTRDCalOnlineGainTableROC;
}

