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
Revision 1.3  2006/03/28 14:56:48  arcelli
updates to handle new V5 geometry & some re-arrangements

Revision 1.2  2006/02/13 17:22:26  arcelli
just Fixing Log info

Revision 1.1  2006/02/13 16:10:48  arcelli
Add classes for TOF Calibration (C.Zampolli)

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
#include "AliLog.h"
#include "AliTOFGeometryV5.h"
#include "AliTOFCalSector.h"
#include "AliTOFCal.h"
#include "AliTOFChannel.h"

extern TROOT *gROOT;

ClassImp(AliTOFCal)

//________________________________________________________________

AliTOFCal::AliTOFCal():TObject(){
  //main ctor
  fGeom = 0x0;
  fNSector = 0;
  fNPlate  = 0;
  fNStripA = 0;
  fNStripB = 0;
  fNStripC = 0;
  fNpadZ = 0;
  fNpadX = 0;
  fnpad = 0;
  fPads = 0x0;
  gROOT->GetListOfBrowsables()->Add(this);
 }
//________________________________________________________________

AliTOFCal::AliTOFCal(AliTOFGeometry *geom):TObject(){
  //ctor with geom
  fGeom = geom;
  fNSector = fGeom->NSectors();
  fNPlate  = fGeom->NPlates();
  fNStripA = fGeom->NStripA();
  fNStripB = fGeom->NStripB();
  fNStripC = fGeom->NStripC();
  fNpadZ = fGeom->NpadZ();
  fNpadX = fGeom->NpadX();
  fnpad = fNSector*(2*(fNStripC+fNStripB)+fNStripA)*fNpadZ*fNpadX;
  fPads = 0x0;
  gROOT->GetListOfBrowsables()->Add(this);
}
//________________________________________________________________

AliTOFCal::AliTOFCal(const AliTOFCal& cal):
  TObject(cal)
  {
    //copy ctor 
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
  //dtor
  gROOT->GetListOfBrowsables()->Remove(this);
  delete [] fPads;
}
//________________________________________________________________

void AliTOFCal::Browse(TBrowser *b)
{
  //add cal obj to list of browsables
  char name[10];
  for(Int_t i=0; i<fNSector; ++i) {
    snprintf(name,sizeof(name),"Sector %2.2d",i);
    b->Add(new AliTOFCalSector(&fPads[i*fnpad/fNSector]),name);
  }
}
//________________________________________________________________

void AliTOFCal::CreateArray(){
  //create cal channel array
  if(fGeom==0x0){
    AliInfo("V5 TOF Geometry is taken as a default");
    AliTOFGeometry *geom= new AliTOFGeometryV5();
    fNSector = geom->NSectors();
    fNPlate  = geom->NPlates();
    fNStripA = geom->NStripA();
    fNStripB = geom->NStripB();
    fNStripC = geom->NStripC();
    fNpadZ = geom->NpadZ();
    fNpadX = geom->NpadX();
    fnpad = fNSector*(2*(fNStripC+fNStripB)+fNStripA)*fNpadZ*fNpadX;
    delete geom;
  }
  fPads= new AliTOFChannel[fnpad];
}
