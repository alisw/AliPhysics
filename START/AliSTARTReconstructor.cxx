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

#include <Riostream.h>

#include <TDirectory.h>

#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliSTARTLoader.h"
#include "AliSTARTdigit.h"
#include "AliSTARTReconstructor.h"
#include <AliESD.h>
#include "AliSTARTRecPoint.h"

ClassImp(AliSTARTReconstructor)

  void AliSTARTReconstructor::Reconstruct(/*AliRunLoader* runLoader*/) 
{
// nothing to be done
}

void AliSTARTReconstructor::FillESD(AliRunLoader* rl, AliESD *pESD) const
{
  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/
  
  Float_t c = 0.3;  //speed of light mm/ps
  Int_t channelWigth=25; //ps
  if (!rl) {
    Error("Reconstruct", "No run loader");
    return;
  }

  if (rl->GetDebug()>1) Info("Reconstruct","START!!!");

  AliSTARTLoader* pStartLoader = (AliSTARTLoader*) rl->GetLoader("STARTLoader");
 
  pStartLoader->LoadDigits();
  AliSTARTdigit* pDigits=pStartLoader->Digits();
  if (!pDigits) {
    Error("Reconstruct", "no digits found");
    return;
  }

  if (rl->GetDebug()>1) pDigits->Dump();
  if(pDigits) {
    Int_t   besttimeright = pDigits->GetBestTimeRight();
    Int_t besttimeleft  = pDigits->GetBestTimeLeft();
    Float_t besttimerightPs = Float_t (besttimeright*channelWigth);
    Float_t besttimeleftPs  = Float_t (besttimeleft*channelWigth);
   Float_t Zposit=(c*(besttimerightPs-besttimeleftPs)-(3500.-697))/2;
  
    
    pESD->SetT0zVertex(Zposit);
    
    if (rl->GetDebug()>1) {
      cout<<" vertex in ESD "<< pESD->GetT0zVertex()<<endl;
    }
    
  } // vertex in 3 sigma
  pStartLoader->UnloadDigits();
 }






