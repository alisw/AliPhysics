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
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TGeometry.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"

// --- Standard library ---
#include <iostream.h>
#include <stdlib.h>   

// --- AliRoot header files ---
#include "AliRun.h" 
#include "AliPHOSClusterizer.h"
#include "AliHeader.h" 
#include "AliPHOSGetter.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSDigitizer.h"

ClassImp(AliPHOSClusterizer)

//____________________________________________________________________________
  AliPHOSClusterizer::AliPHOSClusterizer():TTask("","")
{
  // ctor
  fSplitFile= 0 ; 
  fToSplit  = kFALSE ;

}

//____________________________________________________________________________
AliPHOSClusterizer::AliPHOSClusterizer(const char* headerFile, const char* name, const Bool_t toSplit):
TTask(name, headerFile)
{
  // ctor
  fToSplit  = toSplit ;
  fSplitFile= 0 ; 

}

//____________________________________________________________________________
AliPHOSClusterizer::~AliPHOSClusterizer()
{
  // dtor
         
  fSplitFile = 0 ; 
}


