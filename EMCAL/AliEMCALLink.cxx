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
//  Algorithm class used only by AliEMCALTrackSegmentMaker 
//  Links recpoints into tracksegments                               
//*-- Author: Dmitri Peressounko (SUBATECH)
//*-- Author: Adapted from PHOS by Y. Schutz (SUBATECH)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliEMCALLink.h"

ClassImp(AliEMCALLink)
//____________________________________________________________________________
  AliEMCALLink::AliEMCALLink(Float_t prod, Int_t eca, Int_t rp)
{
  // ctor

  if (gDebug == 2 ) 
    printf("ctor: prod = %f, ec=%d , rp=%d", prod, eca, rp) ;  
  fProd   = prod ;  
  fECAN   = eca ;
  fOtherN = rp ;
}

//____________________________________________________________________________
Int_t AliEMCALLink::Compare(const TObject * obj) const
{
  // Compare according to the distance between two recpoints in a track segment 

  Int_t rv ;

  AliEMCALLink * link = (AliEMCALLink *)obj ;

  if(GetProd() < link->GetProd() ) 
    rv = -1 ;
  else 
    rv = 1 ;

  return rv ;
}
