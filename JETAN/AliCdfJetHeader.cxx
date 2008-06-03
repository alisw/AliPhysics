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

//---------------------------------------------------------------------
// PxCone CDF Algorithm Jet finder header class
// Stores parameters of CDF Jet finder
//---------------------------------------------------------------------

#include "AliCdfJetHeader.h"
ClassImp ( AliCdfJetHeader )

////////////////////////////////////////////////////////////////////////

AliCdfJetHeader::AliCdfJetHeader() :
    fRadius (0.7),
    fPtMin  (0.),
    fPtMax (0.),
    fEtaMin (0.),
    fEtaMax (0.),
    fPhiMin (0.),
    fPhiMax (0.)
  {
  // Constructor
  }

////////////////////////////////////////////////////////////////////////

// void AliCdfJetHeader::PrintParameters() const
//   {
//   // prints out parameters of jet algorithm
//
//   cout << " CDF simplified jet algorithm " << endl;
//   cout << " * Jet parameters: " << endl;
//   cout << " Cone size: " << fRadius << endl;
//
//
//   }
