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

//
// *** Class AliRsnPIDIndex ***
//
// It sorts the indexes of all tracks in an AliRsnEvent
// for a fast retrieval of them according to charge and PID.
//
// author: M. Vala (email: martin.vala@cern.ch)
//

#include <TObject.h>

#include "AliLog.h"

#include "AliRsnPIDIndex.h"

ClassImp(AliRsnPIDIndex)

//_____________________________________________________________________________
AliRsnPIDIndex::AliRsnPIDIndex(Int_t num)
{
//
// Default constructor
//
  Int_t i, j, k;
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < AliRsnPID::kSpecies+1; j++)
    {
      fNumOfIndex[i][j] = 0;
      fIndex[i][j].Set(num);
      for (k = 0; k < num; k++) fIndex[i][j].AddAt(-1, k);
    }
  }
}

//_____________________________________________________________________________
AliRsnPIDIndex::AliRsnPIDIndex(const AliRsnPIDIndex & copy)
    : TObject(copy)
{
//
// Copy constructor.
// Creates new instances of all collections
// to store a copy of all objects.
//

  Int_t i, j, k, size;
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < AliRsnPID::kSpecies+1; j++)
    {
      fNumOfIndex[i][j] = copy.fNumOfIndex[i][j];
      size = copy.fIndex[i][j].GetSize();
      fIndex[i][j].Set(size);
      for (k = 0; k < size; k++)
      {
        fIndex[i][j].AddAt(copy.fIndex[i][j].At(k), k);
      }
    }
  }
}

//_____________________________________________________________________________
AliRsnPIDIndex& AliRsnPIDIndex::operator= (const AliRsnPIDIndex & copy)
{
//
// Assignment operator.
// Creates new instances of all collections
// to store a copy of all objects.
//

  Int_t i, j, k, size;
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < AliRsnPID::kSpecies+1; j++)
    {
      fNumOfIndex[i][j] = copy.fNumOfIndex[i][j];
      size = copy.fIndex[i][j].GetSize();
      fIndex[i][j].Set(size);
      for (k = 0; k < size; k++)
      {
        fIndex[i][j].AddAt(copy.fIndex[i][j].At(k), k);
      }
    }
  }

  // return this object
  return (*this);
}

AliRsnPIDIndex::~AliRsnPIDIndex()
{
//
// Destructor.
// Does nothing.
//
}

//_____________________________________________________________________________
void AliRsnPIDIndex::Print(Option_t* /*option*/) const
{
//
// Prints AliRsnPIDIndex info
//
  Int_t i, j;
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < AliRsnPID::kSpecies + 1; j++)
    {
      AliInfo(Form(" [%d][%d] %d %d", i, j, fIndex[i][j].GetSize(), fNumOfIndex[i][j]));
    }
  }
}

//_____________________________________________________________________________
void AliRsnPIDIndex::AddIndex(const Int_t index, Char_t sign, AliRsnPID::EType type)
{
//
// Adds index to corresponding TArrayI
//
  Int_t iCharge = ChargeIndex(sign);
  Int_t iType = (Int_t)type;
  fIndex[iCharge][iType].AddAt(index, fNumOfIndex[iCharge][iType]);
  fNumOfIndex[iCharge][iType]++;
}

//_____________________________________________________________________________
void AliRsnPIDIndex::AddIndex(const Int_t index, Short_t sign, Int_t type)
{
//
// Adds index to corresponding TArrayI
//

  fIndex[sign][type].AddAt(index, fNumOfIndex[sign][type]);
  fNumOfIndex[sign][type]++;
}

//_____________________________________________________________________________
void AliRsnPIDIndex::SetCorrectIndexSize()
{
//
// Sets Correct sizes to all TArrayI
//

  Int_t i, j;
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < AliRsnPID::kSpecies + 1; j++)
    {
      fIndex[i][j].Set(fNumOfIndex[i][j]);
    }
  }
}

//_____________________________________________________________________________
TArrayI* AliRsnPIDIndex::GetTracksArray(Char_t sign, AliRsnPID::EType type)
{
//
// Returns the array of indexes of tracks whose charge
// and PID correspond to the passed arguments:
//   1) sign of particle ('+' or '-')
//   2) PID of particle (from AliRsnPID::EType)
// Otherwise returns null pointer.
//

  Int_t icharge = ChargeIndex(sign);
  if (icharge < 0) return (TArrayI *) 0x0;
  if (type < AliRsnPID::kElectron || type > AliRsnPID::kSpecies)
  {
    AliError(Form("Index %d out of range", type));
    return (TArrayI *) 0x0;
  }

  return &fIndex[icharge][type];
}

//_____________________________________________________________________________
TArrayI* AliRsnPIDIndex::GetCharged(Char_t sign)
{
//
// Returns the array of indexes of tracks whose charge
// corresponds to the passed argument
// Otherwise returns a null pointer.
//

  // check that argument is meaningful
  Int_t icharge = ChargeIndex(sign);
  if (icharge < 0) return (TArrayI *)0x0;

  // count total number of tracks with that charge
  // and create output object of appropriate size
  Int_t i, total = 0;
  for (i = 0; i <= AliRsnPID::kSpecies; i++) total += fIndex[icharge][i].GetSize();
  TArrayI *output = new TArrayI(total);

  // add all indexes
  Int_t j, counter = 0;
  for (i = 0; i <= AliRsnPID::kSpecies; i++)
  {
    for (j = 0; j < fIndex[icharge][i].GetSize(); j++)
    {
      output->AddAt(fIndex[icharge][i].At(j), counter++);
    }
  }

  return output;
}

//_____________________________________________________________________________
Int_t AliRsnPIDIndex::ChargeIndex(Char_t sign) const
{
//
// Returns the array index corresponding to charge
//

  if (sign == '+') return 0;
  else if (sign == '-') return 1;
  else
  {
    AliError(Form("Character '%c' not recognized as charge sign", sign));
    return -1;
  }
}
