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

/*
$Log$
Revision 1.2  2005/11/03 13:09:19  hristov
Removing meaningless const declarations (linuxicc)

Revision 1.1  2005/10/11 12:31:50  masera
Preprocessor classes for SPD (Paul Nilsson)

*/

///////////////////////////////////////////////////////////////////////////
// AliITSBadChannelsSPD implementation by P. Nilsson 2005
// AUTHOR/CONTACT: Paul.Nilsson@cern.ch
//
// The class is used by the AliITSPreprocessorSPD class to store the
// final noisy and dead channel objects in the calibration database for
// the SPD.
//
// An instance of this container class contains all the "bad" channels,
// i.e. the noisy or dead channels of the SPD. It contains TObjArrays
// for each module of the SPD (240 in total for ALICE, and 6 for the 2004
// joint ITS beam test. The instance object should, once filled with data,
// be stored in the calibration database. This is done by the SPD
// preprocessor.
///////////////////////////////////////////////////////////////////////////

#include "AliITSBadChannelsSPD.h"

ClassImp(AliITSBadChannelsSPD)

//__________________________________________________________________________
AliITSBadChannelsSPD::AliITSBadChannelsSPD(void) :
  fIndexArraySize(0),
  fBadChannelsArraySize(0),
  fIndexArray(0x0),
  fBadChannelsArray(0x0)
{
  // Default constructor
}

//__________________________________________________________________________
AliITSBadChannelsSPD::AliITSBadChannelsSPD(const AliITSBadChannelsSPD &bc) :
  TObject(bc),
fIndexArraySize(bc.fIndexArraySize),
fBadChannelsArraySize(bc.fBadChannelsArraySize),
fIndexArray(0),
fBadChannelsArray(0){
  // Default copy constructor

  // Create new arrays
  fIndexArray = new Int_t[fIndexArraySize];
  fBadChannelsArray = new Int_t[fBadChannelsArraySize];

  // Copy arrays
  for (Int_t i = 0; i < fIndexArraySize; i++)
    {
      fIndexArray[i] = bc.fIndexArray[i];
    }
  for (Int_t i = 0; i < fBadChannelsArraySize; i++)
    {
      fBadChannelsArray[i] = bc.fBadChannelsArray[i];
    }
}

//__________________________________________________________________________
AliITSBadChannelsSPD::~AliITSBadChannelsSPD(void)
{
  // Default destructor

  delete [] fIndexArray;
  fIndexArray = 0;
  delete [] fBadChannelsArray;
  fBadChannelsArray = 0;
}

//__________________________________________________________________________
AliITSBadChannelsSPD& AliITSBadChannelsSPD::operator=(const AliITSBadChannelsSPD &bc)
{
  // Assignment operator
  
  // Guard against self-assignment
  if (this != &bc)
    {
      // Copy array sizes
      fIndexArraySize = bc.fIndexArraySize;
      fBadChannelsArraySize = bc.fBadChannelsArraySize;

      delete [] fIndexArray;
      fIndexArray = new Int_t[fIndexArraySize];

      delete [] fBadChannelsArray;
      fBadChannelsArray = new Int_t[fBadChannelsArraySize];

      // Copy arrays
      for (Int_t i = 0; i < fIndexArraySize; i++)
	{
	  fIndexArray[i] = bc.fIndexArray[i];
	}
      for (Int_t i = 0; i < fBadChannelsArraySize; i++)
	{
	  fBadChannelsArray[i] = bc.fBadChannelsArray[i];
	}
    }

  return *this;
}


//__________________________________________________________________________
void AliITSBadChannelsSPD::Put(Int_t* &badChannelsArray, const Int_t &badChannelsArraySize,
			       Int_t* &indexArray, const Int_t &indexArraySize)
{
  // Add the bad channels and index arrays

  fIndexArraySize = indexArraySize;
  fBadChannelsArraySize = badChannelsArraySize;

  fIndexArray = new Int_t[fIndexArraySize];
  fBadChannelsArray = new Int_t[fBadChannelsArraySize];

  // Copy the arrays
  for (Int_t i = 0; i < fIndexArraySize; i++)
    {
      fIndexArray[i] = indexArray[i];
    }
  for (Int_t i = 0; i < fBadChannelsArraySize; i++)
    {
      fBadChannelsArray[i] = badChannelsArray[i];
    }

}

//__________________________________________________________________________
Bool_t AliITSBadChannelsSPD::Get(Int_t* &badChannelsArray, Int_t* &indexArray) const
{
  // Get the bad channels and the index arrays

  Bool_t status = kTRUE;

  // Set the array pointers
  if (fIndexArraySize > 0)
    {
      badChannelsArray = fBadChannelsArray;
      indexArray = fIndexArray;
    }
  else
    {
      status = kFALSE;
    }

  return status;
}

//__________________________________________________________________________
Int_t* AliITSBadChannelsSPD::CreateModuleArray(Int_t module) const
{
  // Create an Int_t array for a given module

  Int_t *moduleArray = 0;

  const Int_t kSize = AliITSBadChannelsSPD::GetModuleArraySize(module);
  if (kSize > 0)
    {
      // Create a new array
      moduleArray = new Int_t[kSize];

      // Copy the module data into the module array from the bad channels array
      const Int_t kPosition = fIndexArray[module];
      for (Int_t index = 0; index < kSize; index++)
	{
	  moduleArray[index] = fBadChannelsArray[kPosition + index];
	}
    }

  return moduleArray;
}

//__________________________________________________________________________
TObjArray* AliITSBadChannelsSPD::CreateModuleObjArray(Int_t module) const
{
  // Create a TObjArray for a given module

  TObjArray *moduleArray = 0;

  const Int_t kSize = AliITSBadChannelsSPD::GetModuleObjArraySize(module);
  if (kSize > 0)
    {
      // Create a new array
      moduleArray = new TObjArray(kSize);

      // Copy the module data into the module array from the bad channels array

      // Get the start position of the data (skip the number of bad channels, i.e. the first stored number)
      Int_t position = fIndexArray[module] + 1;

      Int_t i = 0;
      while (i < kSize)
	{
	  // Create and add the current channel
	  AliITSChannelSPD *channel = new AliITSChannelSPD(fBadChannelsArray[position++], fBadChannelsArray[position++]);
	  moduleArray->Add(channel);

	  // Go to next bad channel
	  i++;
	}
    }

  return moduleArray;
}
