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

#include "AliRICHHit.h" //class header
#include <TPDGCode.h>   //Print() 
 
ClassImp(AliRICHHit)
//__________________________________________________________________________________________________
void AliRICHHit::Print(Option_t*)const
{
//Print hit
  char *sPart=Form("pid=%i",fPid);
  switch(fPid){
    case kProton:      sPart="p+  ";break;
    case kProtonBar:   sPart="a-  ";break;
    case kKPlus:       sPart="K+  ";break;
    case kKMinus:      sPart="K-  ";break;
    case kPiPlus:      sPart="pi+ ";break;
    case kPiMinus:     sPart="pi- ";break;
    case kMuonPlus:    sPart="mu+ ";break;
    case kMuonMinus:   sPart="mu- ";break;
    case kElectron:    sPart="e-  ";break;
    case kPositron:    sPart="e+  ";break;
    case 50000050:     sPart="ckov";break;
    case 50000051:     sPart="feed";break;
  }

  Printf("%s TID=%6i,Ch=(%2i),E=%9.3f eV, pos=(%7.2f,%7.2f,%7.2f)cm",sPart,Track(),C(),fE*1e9,fX,fY,fZ);
}
//__________________________________________________________________________________________________
