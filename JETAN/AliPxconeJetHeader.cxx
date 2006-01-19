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
// Jet header base class 
// Stores a comment which describes the jet analysis
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include <Riostream.h> 
#include "AliPxconeJetHeader.h"
ClassImp(AliPxconeJetHeader)
 
 
////////////////////////////////////////////////////////////////////////

AliPxconeJetHeader::AliPxconeJetHeader():
  AliJetHeader("AliPxconeJetHeader")
{
  // Constructor
  SetMode();
  SetRadius();
  SetMinPt();
  SetOverlap();
}
 

////////////////////////////////////////////////////////////////////////
 
void AliPxconeJetHeader::PrintParameters() const

{
  // prints out parameters of jet algorithm
  cout << " PXCONE jet algorithm " << endl;
  cout << "  Running mode: " << fMode << endl;
  cout << "  Cone size: " << fRadius << endl;
  cout << "  Minimum jet energy: " << fMinPt << endl;
  cout << "  Overlap fraction: " << fOverlap << endl;
}
