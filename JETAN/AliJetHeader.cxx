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

#include "AliJetHeader.h"
ClassImp(AliJetHeader)
 
 
////////////////////////////////////////////////////////////////////////

AliJetHeader::AliJetHeader():
  TNamed("AliJetHeader", "Jet Header"),
  fComment("No comment"),
  fJetEtaMax(0.5),
  fJetEtaMin(-0.5),
  fJetPhiMax(2*TMath::Pi()),
  fJetPhiMin(0.),
  fDebug(0)
{
  //
  // Default constructor
  //

}
 
////////////////////////////////////////////////////////////////////////

AliJetHeader::AliJetHeader(const char * name):
  TNamed(name, "Jet Header"),
  fComment("No comment"),
  fJetEtaMax(0.5),
  fJetEtaMin(-0.5),
  fJetPhiMax(2*TMath::Pi()),
  fJetPhiMin(0.),
  fDebug(0)
{
  //
  // Constructor
  //

}

////////////////////////////////////////////////////////////////////////

 
