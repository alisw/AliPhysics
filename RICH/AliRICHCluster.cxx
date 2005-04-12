//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************

#include "AliRICHCluster.h"
#include <AliLog.h>

 
ClassImp(AliRICHCluster)
//__________________________________________________________________________________________________
void AliRICHCluster::Print(Option_t*)const
{
//Print current cluster  
  const char *status=0;
  switch(fStatus){
    case      kRaw: status="raw"     ;break;
    case kResolved: status="resolved";break;
    case    kEmpty: status="empty"   ;break;
  }
  if(fDigits)    
    ::Info("cluster","cfm=%10i, cs=%2i, SiMa=%6i, Shape=%5i, x=%7.3f, y=%7.3f, Q=%6i, %s with %i digits",
                             fCFM,fChamber,fSize,fShape,fX,fY,fQdc,status,fDigits->GetEntriesFast());
  else
    AliInfo(Form("cfm=%10i, cs=%2i, SiMa=%6i, Shape=%5i, x=%7.3f, y=%7.3f, Q=%6i, %s with %i digits",
                             fCFM,fChamber,fSize,fShape,fX,fY,fQdc,status,0));
    
}//Print()
//__________________________________________________________________________________________________
