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

#include "AliRun.h"
#include "AliSTART.h"
#include "AliSTARTdigit.h"
#include "AliSTARThit.h"
#include "AliSTARTvertex.h"

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

void AliSTARTvertex::Reconstruct() 
{
  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/

  Int_t timediff;
  Float_t timePs;
  char nameTD[8],nameTR[8];
 char filename[100];
  sprintf(filename,"galice.root");
  AliRunLoader* rl = AliRunLoader::Open("galice.root",AliConfig::fgkDefaultEventFolderName,"read");
  if (rl == 0x0)
   {
     cerr<<"Can not open session for file galice.root\n";
     return;
   }

  rl->LoadgAlice();
  gAlice = rl->GetAliRun();
  
  AliSTART* START  = (AliSTART *)gAlice->GetDetector("START");
  
  rl->LoadHeader();
  rl->LoadKinematics("READ");
  Int_t retval;
  AliLoader* lstart = rl->GetLoader("STARTLoader");
  lstart->LoadDigits("READ");
 
  AliSTARTdigit *digits;
  AliSTARTvertex *fvertex;
 
  digits = new AliSTARTdigit();
  fvertex = new AliSTARTvertex();

 // Event ------------------------- LOOP  
   
  // gAlice->GetEvent(evNumber);
  Int_t iNevents=rl->GetNumberOfEvents();
  cout<<"  nevents   "<<iNevents<<endl;
  
  for (Int_t evNumber=0; evNumber<iNevents; evNumber++){
    rl->GetEvent(evNumber);
    lstart->LoadDigits("READ");
    gDirectory->ls();

    sprintf(nameTD,"START_D_%d",evNumber);
    TObject *td = (TObject*)gDirectory->Get(nameTD);
    printf("%s\n",nameTD);
    //   td->Dump();
    if (!td) {
      cerr<<"something wrong with input...."<<endl;
      exit(111);
    }
    td->Read(nameTD);
    digits->Read(nameTD);
    if(digits->GetTimeDiff()<TMath::Abs(1000))
      {
	timediff=digits->GetTimeDiff();     //time in number of channels
	timePs=(512-timediff)*2.5;       // time in Ps channel_width =10ps
	cout<<"timediff "<< timediff<<" timePs "<<timePs<<endl;
	// Float_t c = 299792458/1.e9;  //speed of light cm/ps
	Float_t c = 0.3;  //speed of light mm/ps
	Float_t Zposit=timePs*c;// for 0 vertex
	cout<<" Zposit "<<Zposit<<endl;
	fvertex->Set((Int_t) Zposit);
      }
    lstart->LoadRecPoints("UPDATE");
    sprintf(nameTR,"START_R_%d",evNumber);
    printf("%s\n",nameTR);
    fvertex->Write(nameTR);


  }
}






