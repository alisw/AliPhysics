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
/*
$Log$
Revision 1.1  2000/03/24 17:46:58  alla
Vertex reconstruction

*/ 
#include <TObject.h>
#include "AliSTARTvertex.h"
#include "AliSTARTdigit.h"
#include "AliSTARThit.h"
#include "AliSTART.h"
#include "AliRun.h"
#include "AliMC.h"

ClassImp(AliSTARTvertex)

AliSTARTvertex::AliSTARTvertex( Int_t * Zposit)
{
  //
  //     The creator for the AliSTARTvertex class. This routine fills the
  // AliSTARTvertex data members from the array vertex.
  // The order of the elements in the vertex array are
  // fEvent = digits[0], fZposition = vertex[1],
  // fTime_diff = Vertex[2]
  // Therefore the array digits is expected to be at least 3 elements long.
  //

  Zposit = &fZposition ;
}

void AliSTARTvertex::Reconstruct(Int_t evNumber=1) 
{
 
  Int_t timediff;
  Float_t timePs;
  char nameTD[8],nameTR[8];

  TBranch *bRec=0;
  TBranch *bd;
  AliSTARTdigit *digits;
  AliSTARTvertex *vertex;
 
  Int_t buffersize=256;
  Int_t split=1;
 
  // TParticle *particle;
  digits = new AliSTARTdigit();
  vertex = new AliSTARTvertex();

 // Event ------------------------- LOOP  
   
  sprintf(nameTD,"TreeD%d",evNumber);
  printf("%s\n",nameTD);
  TTree *td = (TTree*)gDirectory->Get(nameTD);
  bd = td->GetBranch("START");
  bd->SetAddress(&digits);
  bd->GetEvent(0);
  printf(" Digits: "); digits->MyDump();    
  sprintf(nameTR,"TreeR%d",evNumber);
  TTree *tr = new TTree(nameTR,"START");
  bRec = tr->Branch("START","AliSTARTvertex",&vertex,buffersize,split);

  //  td->Print(); td->Show(0); td->GetBranch("START")->Dump();
        digits->MyDump();
       printf("digits-> %d \n",digits->GetTime());
  
  if(digits->GetTime()!=999999)
    {
      timediff=digits->GetTime();     //time in number of channels
      timePs=(timediff-128)*10.;       // time in Ps channel_width =10ps
      printf(" timediff %d in PS %f\n",timediff,timePs);
      Float_t c = 299792458/1.e9;  //speed of light cm/ps
      //Float_t c = 0.3;  //speed of light mm/ps
      Float_t Zposit=timePs*c;// for 0 vertex
      //    Float_t Zposit=timePs*c/2.;// for spread vertex
      //      printf(" Z position %f\n",Zposit);
      //      vertex->GetVertex();
      vertex->Set(Zposit);
      tr->Fill();
      tr->Write();
      //hTimediff->Fill(timePs);
      //hVertex->Fill(Zposit);
      }

}




