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
#include "AliSTARTdigit.h"
#include "AliSTARThit.h"
#include "AliSTARTvertex.h"
#include <AliESD.h>

ClassImp(AliSTARTvertex)

AliSTARTvertex::AliSTARTvertex( Int_t * Zposit)
{
  //
  //     The creator for the AliSTARTvertex class. This routine fills the
  // AliSTARTvertex data members from the array vertex.
  // The order of the elements in the vertex array are
  //  fZposition = vertex[0],
  //
 
  Zposit = &fZposition ;
}

void AliSTARTvertex::Reconstruct(AliESD *pESD) 
{
  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/

  Int_t timediff;
  Float_t timePs;
  
  
  AliRunLoader* rl = AliRunLoader::Open("galice.root");
    //,AliConfig::fgkDefaultEventFolderName,"read");
  if (rl == 0x0)
   {
     cerr<<"Can not open session for file galice.root\n";
     return;
   }
  
#ifdef DEBUG
 Info("Reconstruct","START!!!");
#endif
  AliLoader* pStartLoader = rl->GetLoader("STARTLoader");
 
 // Event ------------------------- LOOP  
   

  Int_t iNevents=rl->GetNumberOfEvents();
  
  for (Int_t evNumber=0; evNumber<iNevents; evNumber++)
    {
    rl->GetEvent(evNumber);
  pStartLoader ->LoadDigits("READ");

#ifdef DEBUG
  gDirectory->ls();
#endif
    AliSTARTdigit* pDigits=(AliSTARTdigit*)gDirectory->Get("START_D");

#ifdef DEBUG
    pDigits->Dump();
#endif  
     if(pDigits->GetTimeDiff()<TMath::Abs(1000))
      {
	timediff=pDigits->GetTimeDiff();     //time in number of channels
	timePs=(512-timediff)*2.5;       // time in Ps channel_width =10ps
	// Float_t c = 299792458/1.e9;  //speed of light cm/ps
	Float_t c = 0.3;  //speed of light mm/ps
	Float_t Zposit=timePs*c;// for 0 vertex
#ifdef DEBUG
	cout<<"timediff "<< timediff<<" timePs "<<timePs<<" Zposit "<<Zposit<<endl;
#endif 
	pESD->SetT0zVertex(Zposit);

#ifdef DEBUG
	cout<<" vertex in ESD "<< pESD->GetT0zVertex()<<endl;
#endif
      } // vertex in 3 sigma
     
    } //event loop
 }






