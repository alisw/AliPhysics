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
//  Algorithm class for the identification of particles detected in EMCAL        
//  base  class  of identified particle  
//  Why should I put meaningless comments
//  just to satisfy
//  the code checker                
                         
//                  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko


// --- ROOT system ---
#include "TGeometry.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
 
// --- Standard library ---
#include <stdlib.h>


// --- AliRoot header files ---
#include "AliRun.h" 
#include "AliEMCALPID.h"
#include "AliHeader.h" 

ClassImp(AliEMCALPID)

//____________________________________________________________________________
  AliEMCALPID::AliEMCALPID():TTask("","")
{
  // ctor
  fSplitFile= 0 ; 

}


//____________________________________________________________________________
AliEMCALPID::AliEMCALPID(const char* headerFile, const char * name, const Bool_t toSplit):TTask(name, headerFile)
{
  // ctor

  fToSplit = toSplit ;
  fSplitFile= 0 ; 
}

//____________________________________________________________________________
AliEMCALPID::~AliEMCALPID()
{
  // dtor
        
  fSplitFile = 0 ;
}
