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
// Class defining the digit object                               //
// for ITS upgrade                                               //
// Inherits from AliDigit                                        //
//                                                               //
///////////////////////////////////////////////////////////////////

ClassImp(AliITSDigitUpgrade)
//______________________________________________________________________
  AliITSDigitUpgrade::AliITSDigitUpgrade():
    AliDigit(),
    fPixId(9999),
    fSignal(0),
    fNLayer(0),
    fModule(0),
    fNelectrons(0),
    fNTracksIdMC(0)
{
 //
 // default constructor
 //
  for(Int_t i=0; i<kMaxLab ; i++) {
    fTrackIdMC[i]=-1;
    fSignalID[i]=-1;
  }
} 
//_______________________________________________________________________
AliITSDigitUpgrade::AliITSDigitUpgrade(Int_t *digits): 
  AliDigit(digits),
  fPixId(9999), 
  fSignal(0),
  fNLayer(0),
  fModule(0), 
  fNelectrons(0),
  fNTracksIdMC(0)
{
 //
 // constructor
 // 
  for(Int_t i=0; i<kMaxLab ; i++) {
    fTrackIdMC[i]=-1;
    fSignalID[i]=-1;
  }
} 
//____________________________________________________________________________________________________
AliITSDigitUpgrade::AliITSDigitUpgrade(ULong_t pixid, Float_t eloss): 
  AliDigit(),
  fPixId(pixid), 
  fSignal(eloss),
  fNLayer(0), 
  fModule(0),
  fNelectrons(0),
  fNTracksIdMC(0)
{
 //
 // Used constructor in simulation
 //
  for(Int_t i=0; i<kMaxLab ; i++) {
    fTrackIdMC[i]=-1;
    fSignalID[i]=-1;
  }
} 
//____________________________________________________________________________________________________
AliITSDigitUpgrade::AliITSDigitUpgrade(const AliITSDigitUpgrade &d):
  AliDigit(d),
  fPixId(d.fPixId),
  fSignal(d.fSignal),
  fNLayer(d.fNLayer), 
  fModule(d.fModule), 
  fNelectrons(d.fNelectrons),
  fNTracksIdMC(d.fNTracksIdMC)
{
 // 
 // copy constructor
 //
  for(Int_t i=0; i<kMaxLab ; i++) {
    fTrackIdMC[i]=d.fTrackIdMC[i];
    fSignalID[i]=d.fSignalID[i];
  }
  for(Int_t i=0; i<3 ; i++) fSignalID[i]=d.fSignalID[i];
} 

//____________________________________________________________________________________________________
void AliITSDigitUpgrade::AddTrackID(Int_t tid) { 
  // 
  // Add an MC label (track ID) to the "expanded list"
  //
  if (fNTracksIdMC==kMaxLab) {
    AliWarning("Max. numbers of labels reached!"); 
  } else {
    fTrackIdMC[fNTracksIdMC]=tid; 
    fNTracksIdMC++;
  } 
}

//____________________________________________________________________________________________________
void  AliITSDigitUpgrade::GetPosition(Int_t ilayer, Int_t nx, Int_t nz, Double_t &xloc, Double_t &zloc){
  //
  // Determines the local coordinates of the center of the pixel (a digit is a pixel)
  //

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
  //
  //Standard output format for this class
  //

  Double_t xz[2]={-1,-1};
  GetPosition(fNLayer,GetxPixelNumber(),GetzPixelNumber(),xz[0],xz[1]);
   
  printf("pixid  %10.0i (%6.3f,%6.3f) in layer %i \n",(Int_t)fPixId,xz[0],xz[1],fNLayer);
  printf("pixid  %u ",(UInt_t)fPixId);
  printf(" (xloc, zloc)= (%6.3f, %6.3f) in layer %i and module %i \n",xz[0],xz[1],fNLayer, fModule);
  printf(" Eloss %f  Nel %f \n ",fSignal, fNelectrons);
  printf(" MC Track Ids =(");
  if (kMaxLab<=fNTracksIdMC) {
    for (Int_t i=0; i<kMaxLab; i++) { printf("%d,",fTrackIdMC[i]); }
    printf(")\n ElossID = (");
    for (Int_t i=0; i<kMaxLab; i++) { printf("%lf,",fSignalID[i]); }
  } else {
    for (Int_t i=0; i<fNTracksIdMC; i++) { printf("%d,",fTrackIdMC[i]); }
    printf(")\n ElossID = (");
    for (Int_t i=0; i<fNTracksIdMC; i++) { printf("%lf,",fSignalID[i]); }
  }
  printf(")\n");
}

