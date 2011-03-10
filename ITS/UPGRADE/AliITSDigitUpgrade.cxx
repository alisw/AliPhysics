/**************************************************************************
 * Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliITSDigitUpgrade.h"
#include "AliITSsegmentationUpgrade.h"
#include "AliLog.h"

///////////////////////////////////////////////////////////////////
//                                                               //
// Class defining the digit object
// for ITS upgrade
// Inherits from AliDigit
//                                                               //
///////////////////////////////////////////////////////////////////

ClassImp(AliITSDigitUpgrade)
//______________________________________________________________________
  AliITSDigitUpgrade::AliITSDigitUpgrade():AliDigit(),
					   fPixId(9999),
					   fSignal(0),
					   fNLayer(0),
						fModule(0),
					   fNelectrons(0)
{for(Int_t i=0; i<3 ; i++) fSignalID[i]=-1;} //default creator
//_______________________________________________________________________
AliITSDigitUpgrade::AliITSDigitUpgrade(Int_t *digits): AliDigit(digits),
						       fPixId(9999), 
						       fSignal(0),
						       fNLayer(0),
							fModule(0), 
						       fNelectrons(0)
{for(Int_t i=0; i<3 ; i++) fSignalID[i]=-1;} //standard creator digits only
//____________________________________________________________________________________________________
AliITSDigitUpgrade::AliITSDigitUpgrade(ULong_t pixid, Float_t eloss): AliDigit(),
								      fPixId(pixid), 
								      fSignal(eloss),
								      fNLayer(0), 
									fModule(0),
								      fNelectrons(0)
{for(Int_t i=0; i<3 ; i++) fSignalID[i]=-1;} //standard creator digits only
//____________________________________________________________________________________________________
AliITSDigitUpgrade::AliITSDigitUpgrade(const AliITSDigitUpgrade &d):AliDigit(d),
								    fPixId(d.fPixId),
								    fSignal(d.fSignal),
								    fNLayer(d.fNLayer), 
								    fModule(d.fModule), 
								    fNelectrons(d.fNelectrons)
{for(Int_t i=0; i<3 ; i++) fSignalID[i]=d.fSignalID[i];} //copy constructor
//____________________________________________________________________________________________________
void  AliITSDigitUpgrade::GetPosition(Int_t ilayer, Int_t nx, Int_t nz, Double_t &xloc, Double_t &zloc){
  AliITSsegmentationUpgrade *s =new AliITSsegmentationUpgrade();
  if(s->GetCellSizeX(ilayer)!=0) xloc= (nx)*(s->GetCellSizeX(ilayer))+0.5*(s->GetCellSizeX(ilayer));
  else AliError("Upgrade segmentation not initalized");

  if(s->GetCellSizeZ(ilayer)!=0)
    zloc=(nz)*(s->GetCellSizeZ(ilayer))+0.5*(s->GetCellSizeZ(ilayer))-(s->GetHalfLength(ilayer));
  else AliError("Upgrade segmentation not initalized");
  delete s;
}
//____________________________________________________________________________________________________
void AliITSDigitUpgrade::PrintInfo(){
  //Standard output format for this class
  Double_t xz[2]={-1,-1};
  GetPosition(fNLayer,GetxPixelNumber(),GetzPixelNumber(),xz[0],xz[1]);
   
  printf("pixid  %10.0i (%6.3f,%6.3f) in layer %i \n",(Int_t)fPixId,xz[0],xz[1],fNLayer);
  printf("pixid  %u ",(UInt_t)fPixId);
  printf(" (xloc, zloc)= (%6.3f, %6.3f) in layer %i and module %i \n",xz[0],xz[1],fNLayer, fModule);
  printf(" Eloss %f  Nel %f   track ID %i   %i  %i ", fSignal, fNelectrons,fTracks[0],fTracks[1],fTracks[2]);
  printf(" ElossID %f  %f %f  \n", fSignalID[0],fSignalID[1],fSignalID[2]);
}

