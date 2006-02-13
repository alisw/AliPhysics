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

/*$Log$
author: Chiara Zampolli, zampolli@bo.infn.it
 */  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF calibration : array of AliTOFChannels                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TClass.h"
#include "AliTOFGeometryV4.h"
#include "AliTOFCalSector.h"
#include "AliTOFCal.h"
#include "AliTOFChannel.h"

extern TROOT *gROOT;

ClassImp(AliTOFCal)

//________________________________________________________________

AliTOFCal::AliTOFCal():TObject(){
  // fCalp = 0;
  fNSector = AliTOFGeometryV4::NSectors();
  fNPlate = AliTOFGeometryV4::NPlates();
  fNStripA = AliTOFGeometryV4::NStripA();
  fNStripB = AliTOFGeometryV4::NStripB();
  //  fNStripC = AliTOFGeometryV4::NStripC();
  fNStripC = 20;
  fNpadZ = AliTOFGeometryV4::NpadZ();
  fNpadX = AliTOFGeometryV4::NpadX();
  fnpad = 0;
  fPads = 0x0;
  gROOT->GetListOfBrowsables()->Add(this);
}
//________________________________________________________________

AliTOFCal::AliTOFCal(const AliTOFCal& cal):
  TObject(cal)
  {
    fNSector = cal.fNSector;
    fNPlate = cal.fNPlate;
    fNStripA = cal.fNStripA;
    fNStripB = cal.fNStripB;
    fNStripC = cal.fNStripC;
    fNpadZ = cal.fNpadZ;
    fNpadX = cal.fNpadX;
    fnpad = cal.fnpad;
    for (Int_t i = 0; i<fnpad; i++){
      fPads[i]=cal.fPads[i];
    }
    gROOT->GetListOfBrowsables()->Add(this);
  }
//____________________________________________________________________________ 
AliTOFCal::~AliTOFCal()
{
  delete [] fPads;
}
//________________________________________________________________

void AliTOFCal::Browse(TBrowser *b)
{
  char name[10];
  for(Int_t i=0; i<fNSector; ++i) {
    snprintf(name,sizeof(name),"Sector %2.2d",i);
    b->Add(new AliTOFCalSector(&fPads[i*fnpad/fNSector]),name);
  }
}
//________________________________________________________________

void AliTOFCal::CreateArray(){
  fnpad = AliTOFGeometryV4::NSectors()*(2*(20+AliTOFGeometryV4::NStripB())+AliTOFGeometryV4::NStripA())*AliTOFGeometryV4::NpadZ()*AliTOFGeometryV4::NpadX();
  fPads= new AliTOFChannel[fnpad];
}
