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

#include "AliHMPIDHit.h"  //class header
#include <TPDGCode.h>     //Draw() Print()
#include <TMarker.h>      //Draw()
#include <TClonesArray.h> //Hit2Sdi()
#include "AliHMPIDParam.h" 
ClassImp(AliHMPIDHit)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDHit::Draw(Option_t*)
{
// Draw option of the hits in the display
  Int_t iMark;
  switch(Pid()){
    case 50000050:   iMark=4;  break;
    case 50000051:   iMark=27; break;
    default:         iMark=26; break;
  }    
  TMarker *pMark=new TMarker(fLx,fLy,iMark); pMark->SetMarkerColor(kRed); pMark->Draw();
}//Draw
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDHit::Hit2Sdi(TClonesArray *pSdiLst,Int_t iHow)const
{
// Adds  sdigits of this hit to the list
// Arguments: pSdiLst- sigits list where to add new sdgits
//            iHow- how many pads to check 
//   Returns: none
  Int_t pc,px,py;
  AliHMPIDParam::Lors2Pad(fLx,fLy,pc,px,py); if(py<0) return; //check if the hit in dead zone. Should never happen during trasport!

  AliHMPIDDigit dig;
  Int_t iSdiCnt=pSdiLst->GetEntries();                       //list of sdigits contains sdigits from previous ivocations of Hit2Sdi, do not override them

  for(Int_t i=-iHow;i<=iHow;i++){                            //horizontal loop
    for(Int_t j=-iHow;j<=iHow;j++){                          //vertical loop
      if(dig.Set(fCh,pc,px+i,py+j,fTrack)) continue;
      dig.SetQ(fQ*dig.IntMathieson(fLx,fLy));
      new((*pSdiLst)[iSdiCnt++]) AliHMPIDDigit(dig);
    }
  }
}//Hit2Sdi
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDHit::Print(Option_t *opt)const
{
//Print hit
  const char *sPart=Form("pid=%i",Pid());
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

  Printf("%sHIT: ch=%i                 (%7.6f,%7.6f) Q=%8.3f TID= %5i, MARS=(%7.2f,%7.2f,%7.2f) %s  %s",
         opt,  Ch(),                    fLx,fLy,  fQ,     fTrack,         X(),  Y(),  Z(),   sPart, 
                        (AliHMPIDParam::IsInDead(LorsX(),LorsY()))? "IN DEAD ZONE":"");
}//Print
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
