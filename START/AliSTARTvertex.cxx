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
Revision 1.6.6.2  2002/07/24 09:50:10  alibrary
Updating VirtualMC

Revision 1.9  2002/07/23 11:48:05  alla
new Digits structure

Revision 1.8  2002/04/16 10:52:41  hristov
Wrong usage of exit() corrected (Sun)

Revision 1.7  2002/04/15 08:04:01  alla
Digits and reconstruction with TObject

Revision 1.6  2001/10/19 05:29:38  alla
bug in meduim fixed

Revision 1.5  2001/07/27 13:03:12  hristov
Default Branch split level set to 99

Revision 1.4  2000/12/22 16:17:15  hristov
Updated  START code from Alla

Revision 1.3  2000/10/02 21:28:13  fca
Removal of useless dependecies via forward declarations
 
Revision 1.2  2000/07/13 16:41:29  fca
New START corrected for coding conventions

Revision 1.1  2000/03/24 17:46:58  alla
Vertex reconstruction

*/ 
#include "TObject.h"
#include "AliSTARTvertex.h"
#include "AliSTARTdigit.h"
#include "AliSTARThit.h"
#include "AliSTART.h"
#include "AliRun.h"
#include "AliMC.h"

//#include "TTree.h"
#include "TDirectory.h"
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

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

void AliSTARTvertex::Reconstruct(Int_t evNumber=1) 
{
  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/

  Int_t timediff;
  Float_t timePs;
  char nameTD[8],nameTR[8];

  AliSTARTdigit *digits;
  AliSTARTvertex *fvertex;
 
  digits = new AliSTARTdigit();
  fvertex = new AliSTARTvertex();

 // Event ------------------------- LOOP  
   
  // gAlice->GetEvent(evNumber);

  sprintf(nameTD,"START_D_%d",evNumber);
  TObject *td = (TObject*)gDirectory->Get(nameTD);
  printf("%s\n",nameTD);
  
  if (!td) {
    cerr<<"something wrong with output...."<<endl;
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
     /*
        TTree *outTreeR = gAlice->TreeR();
    if (!outTreeR) {
      cerr<<"something wrong with output...."<<endl;
      exit(111);
    }
    */
  sprintf(nameTR,"START_R_%d",evNumber);
  printf("%s\n",nameTR);
  //  TDirectory *wd = gDirectory;
  //  outTreeR->GetDirectory()->cd();
    fvertex->Write(nameTR);
    //  wd->cd();
}






