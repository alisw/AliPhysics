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

/* $Id$ */

//_________________________________________________________________________
//  Base class for the clusterization algorithm (pure abstract)
//*--
//*-- Author: Yves Schutz  SUBATECH 
// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TGeometry.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h" 

// --- Standard library ---
#include <Riostream.h>
#include <stdlib.h>   

// --- AliRoot header files ---
#include "AliRun.h" 
#include "AliEMCALClusterizer.h"
#include "AliHeader.h" 
#include "AliEMCALGetter.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALDigitizer.h"

ClassImp(AliEMCALClusterizer)

//____________________________________________________________________________
  AliEMCALClusterizer::AliEMCALClusterizer():TTask("","")
{
  // ctor
  fSplitFile = 0 ;  
  fToSplit  = kFALSE ;

}

//____________________________________________________________________________
AliEMCALClusterizer::AliEMCALClusterizer(const char* headerFile, const char* name, const Bool_t toSplit):
TTask(name, headerFile)
{
  // ctor
  fToSplit  = toSplit ;
  fSplitFile = 0 ;  

}

//____________________________________________________________________________
AliEMCALClusterizer::~AliEMCALClusterizer()
{
  // dtor
  
  fSplitFile = 0 ;
 }

