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
  // Create START digit
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
  // AliSTART *START  = (AliSTART*) gAlice->GetDetector("START");

 // Event ------------------------- LOOP  
   
    sprintf(nameTD,"TreeD%d",evNumber);
    printf("%s\n",nameTD);
    TTree *TD = (TTree*)gDirectory->Get(nameTD);
    bd = TD->GetBranch("START");
    bd->SetAddress(&digits);
    bd->GetEvent(0);
    printf(" Digits: "); digits->MyDump();
    
    sprintf(nameTR,"TreeR%d",evNumber);
    TTree *TR = new TTree(nameTR,"START");
    bRec = TR->Branch("START","AliSTARTvertex",&vertex,buffersize,split);

    //  TD->Print(); TD->Show(0); TD->GetBranch("START")->Dump();
  
    if(digits->fTime_average!=0)
      {
      timediff=digits->fTime_diff;     //time in number of channels
      timePs=(timediff-128)*10.;       // time in Ps
      printf(" timediff %d in PS %f\n",timediff,timePs);
      Float_t c = 299792458/1.e9;  //speed of light cm/ps
      Float_t Zposit=timePs*c;
      printf(" Z position %f\n",Zposit);
      //      vertex->GetVertex();
      vertex->Set(Zposit);
      TR->Fill();
      TR->Write();
      //hTimediff->Fill(timePs);
      //hVertex->Fill(Zposit);
      }

}




