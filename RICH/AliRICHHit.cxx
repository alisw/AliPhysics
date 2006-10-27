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
#include <TString.h>
 
ClassImp(AliRICHHit)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHHit::Print(Option_t*)const
{
//Print hit
  char *sPart=Form("pid=%i",Pid());
  switch(Pid()){
    case kProton:      sPart="p+  ";break;
    case kProtonBar:   sPart="p-  ";break;
    case kKPlus:       sPart="K+  ";break;
    case kKMinus:      sPart="K-  ";break;
    case kPiPlus:      sPart="Pi+ ";break;
    case kPiMinus:     sPart="Pi- ";break;
    case kMuonPlus:    sPart="Mu+ ";break;
    case kMuonMinus:   sPart="Mu- ";break;
    case kElectron:    sPart="e-  ";break;
    case kPositron:    sPart="e+  ";break;
    case 50000050:     sPart="ckov";break;
    case 50000051:     sPart="feed";break;
  }

  Printf("%s Ch:%i TID:%6i,E:%9.3f eV, LORS:(%7.2f,%7.2f) MARS:(%7.2f,%7.2f,%7.2f)cm",sPart,Ch(),Tid(),E()*1e9,LorsX(),LorsY(),X(),Y(),Z());
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
