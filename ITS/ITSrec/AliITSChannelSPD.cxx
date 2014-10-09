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
// AliITSChannelSPD implementation by P. Nilsson 2005
// AUTHOR/CONTACT: Paul.Nilsson@cern.ch
//
// Objects of this class are stored in TObjArrays and should be
// interpreted as "bad" channels, i.e. either noisy or dead channels
// depending on where they are stored.
//
// A channel has the structure:
//
//                           Int_t fColumn:    Column in the SPD module
//                           Int_t fRow:       Row in the SPD module
//
// The class is used by the AliITSPreprocessorSPD class to store noisy
// and dead channels in the calibration database for the SPD.
//
// A channel can be compared with other channels (the equality operator
// is defined). This is e.g. useful for the clustering algorithm. A noisy
// channel should not be used in the clustering
///////////////////////////////////////////////////////////////////////////

#include "AliITSChannelSPD.h"

ClassImp(AliITSChannelSPD)

//__________________________________________________________________________
AliITSChannelSPD::AliITSChannelSPD(void) :
fColumn(-1),
fRow(-1)
{
  // Default constructor
}


//__________________________________________________________________________
AliITSChannelSPD::AliITSChannelSPD(const AliITSChannelSPD &ch) :
  TObject(ch),
fColumn(ch.fColumn),
fRow(ch.fRow){
  // Copy constructor

}

//__________________________________________________________________________
AliITSChannelSPD& AliITSChannelSPD::operator=(const AliITSChannelSPD &ch)
{
  // Assignment operator
  
  // Guard against self-assignment
  if (this != &ch)
    {
      // Copy the data members
      fColumn = ch.fColumn;
      fRow = ch.fRow;
    }
  return *this;
}

//__________________________________________________________________________
Bool_t AliITSChannelSPD::operator==(const AliITSChannelSPD &channel) const
{
  // Equality operator
  // For comparisons between AliITSChannelSPD objects

  return ( ((fColumn == channel.fColumn) && (fRow == channel.fRow)) );
}

//__________________________________________________________________________
AliITSChannelSPD::AliITSChannelSPD(Int_t column, Int_t row):
fColumn(column),
fRow(row){
  // Constructor for already existing channel


}
