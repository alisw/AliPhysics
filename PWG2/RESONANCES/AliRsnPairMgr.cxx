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
// *** Class AliRsnPairMgr ***
//
// TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#include "AliLog.h"

#include "AliRsnPairMgr.h"

ClassImp ( AliRsnPairMgr )

//_____________________________________________________________________________
AliRsnPairMgr::AliRsnPairMgr ( const char*name )
    : TNamed ( name,name),fPairs(0)
{
//=========================================================
// Default constructor
//=========================================================

}

//_____________________________________________________________________________
AliRsnPairMgr::~AliRsnPairMgr()
{
//=========================================================
// Destructor
//=========================================================

}

//_____________________________________________________________________________
void AliRsnPairMgr::AddPair ( AliRsnPair * pair )
{
//=========================================================
// Adds pair
//=========================================================

  fPairs.Add ( ( AliRsnPair * ) pair );
}

//_____________________________________________________________________________
void AliRsnPairMgr::PrintPairs()
{
//=========================================================
// Prints all pairs
//=========================================================
  AliRsnPair * pair=0;
  for ( Int_t i=0;i<fPairs.GetEntriesFast() ;i++ )
  {
    pair = ( AliRsnPair * ) fPairs.At ( i );
    pair->Print();
  }
}
