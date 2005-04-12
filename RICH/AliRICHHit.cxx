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

#include "AliRICHHit.h"
#include <AliLog.h>
 
ClassImp(AliRICHHit)
//__________________________________________________________________________________________________
void AliRICHHit::Print(Option_t*)const
{
//Print hit
  AliInfo(Form("Ch=%1i,TID=%6i,Elos=%9.3f eV,IN(%6.2f,%6.2f,%6.2f)-OUT(%6.2f,%6.2f,%6.2f)=%9.4f"
      ,fChamber,fTrack,fEloss*1e9,fInX3.X() ,fInX3.Y() ,fInX3.Z(),
                                  fOutX3.X(),fOutX3.Y(),fOutX3.Z(),Length()));
}
//__________________________________________________________________________________________________
