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
#include <stdlib.h>

#include <TDirectory.h>
#include <TVirtualMC.h>

#include <AliRun.h>
#include <AliRunLoader.h>
#include "AliSTART.h"
#include "AliSTARTLoader.h"
#include "AliSTARTdigit.h"
#include "AliSTARThit.h"
#include "AliSTARTvertex.h"
#include <AliESD.h>

ClassImp(AliSTARTvertex)

void AliSTARTvertex::Reconstruct(AliRunLoader* rl, AliESD *pESD) 
{
  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/

  Int_t timediff;
  Float_t timePs;
  
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

  if(pDigits->GetTimeDiff()<TMath::Abs(1000)) {
    timediff=pDigits->GetTimeDiff();     //time in number of channels
    timePs=(512-timediff)*2.5;       // time in Ps channel_width =10ps
    // Float_t c = 299792458/1.e9;  //speed of light cm/ps
    Float_t c = 0.3;  //speed of light mm/ps
    Float_t Zposit=timePs*c;// for 0 vertex
    
    if (rl->GetDebug()>1) {
      cout << "timediff " << timediff
	   <<" timePs " << timePs 
	   <<" Zposit "<<Zposit<<endl;
    }
    
    fZposition = Zposit;
    pESD->SetT0zVertex(Zposit);
    
    if (rl->GetDebug()>1) {
      cout<<" vertex in ESD "<< pESD->GetT0zVertex()<<endl;
    }
    
  } // vertex in 3 sigma
  pStartLoader->UnloadDigits();
 }






