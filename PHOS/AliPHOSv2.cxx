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
// Version of AliPHOSv1 which keeps all hits in TreeH
// AddHit, StepManager,and FinishEvent are redefined 
//                  
//*-- Author: Gines MARTINEZ (SUBATECH)
//*-- Modified Nov. 22 2000 by Dmitri Peressounko
// All hits are stored.
// Note, that primaries will not be assigned to digits:
// becouse of tiny energy deposition at each step.
//  

// --- ROOT system ---

#include "TBRIK.h"
#include "TNode.h"
#include "TRandom.h"

// --- Standard library ---

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strstream.h>

// --- AliRoot header files ---

#include "AliPHOSv2.h"
#include "AliPHOSHit.h"
#include "AliPHOSDigit.h"
#include "AliPHOSReconstructioner.h"
#include "AliRun.h"
#include "AliConst.h"

ClassImp(AliPHOSv2)

//____________________________________________________________________________
AliPHOSv2::AliPHOSv2()
{
  // default ctor

}

//____________________________________________________________________________
AliPHOSv2::AliPHOSv2(const char *name, const char *title):
AliPHOSv1(name,title)
{
  // ctor
}

//____________________________________________________________________________
AliPHOSv2::~AliPHOSv2()
{
  // dtor
}

//____________________________________________________________________________
void AliPHOSv2::AddHit(Int_t shunt, Int_t primary, Int_t tracknumber, Int_t Id, Float_t * hits, Int_t pid)
{
  // Add a hit to the hit list.

  AliPHOSHit *newHit ;

  newHit = new AliPHOSHit(shunt, primary, tracknumber, Id, hits, pid) ;

  new((*fHits)[fNhits]) AliPHOSHit(*newHit) ;    
  fNhits++ ;

  delete newHit;

}


