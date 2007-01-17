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
 
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// base class for Jet Reader Header 
// Author: jgcn@mda.cinvestav.mx
//-------------------------------------------------------------------------

#include "AliJetReaderHeader.h"

ClassImp(AliJetReaderHeader)

////////////////////////////////////////////////////////////////////////

AliJetReaderHeader::AliJetReaderHeader():  
 TNamed("AliJetReaderHeader", "Jet Reader Header"),
 fFirst(0),
 fLast(-1),
 fOption(0),
 fSignalPerBg(0),
 fFiducialEtaMin(-0.9),
 fFiducialEtaMax(0.9),
 fPtCut(2.0),
 fComment("No comment"),
 fDir(""),
 fBgDir(""),
 fPattern("")
{
  // Default constructor
}

////////////////////////////////////////////////////////////////////////

AliJetReaderHeader::AliJetReaderHeader(const char * name):  
 TNamed(name, "Jet Reader Header"),
 fFirst(0),
 fLast(-1),
 fOption(0),
 fSignalPerBg(0),
 fFiducialEtaMin(-0.9),
 fFiducialEtaMax(0.9),
 fPtCut(2.0),
 fComment("No comment"),
 fDir(""),
 fBgDir(""),
 fPattern("")
{
  // Constructor
}

////////////////////////////////////////////////////////////////////////

AliJetReaderHeader::~AliJetReaderHeader()
{
  // destructor

}
