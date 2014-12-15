/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class to store the HLT status           //
// (when HLT is in mode C SDD data are compressed,               // 
// see AliITSCompressRawDataSDD)                                 //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "TString.h"
#include "AliITSHLTforSDD.h"

ClassImp(AliITSHLTforSDD)
//______________________________________________________________________
AliITSHLTforSDD::AliITSHLTforSDD():TObject(),fHLTmodeC(0)
{
  // default constructor
}
//______________________________________________________________________
AliITSHLTforSDD::AliITSHLTforSDD(TString hltMode):TObject(),fHLTmodeC(0)
{
  // standard constructor
  TSubString firstChar = hltMode(0,1);
  if (firstChar == "C") fHLTmodeC=kTRUE;
}

