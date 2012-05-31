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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF online calibration : array of AliTOFChannelOnline           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"

#include "AliLog.h"

#include "AliTOFCalOnline.h"
#include "AliTOFGeometryV5.h"

extern TROOT *gROOT;

ClassImp(AliTOFCalOnline)

//________________________________________________________________

AliTOFCalOnline::AliTOFCalOnline():
  TObject(),
  fNSector(0),
  fNPlate(0),
  fNStripA(0),
  fNStripB(0),
  fNStripC(0),
  fNpadZ(0),
  fNpadX(0),
  fnpad(0),
  fGeom(0x0),
  fPads(0x0)
{
  //main ctor
 }
//________________________________________________________________

AliTOFCalOnline::AliTOFCalOnline(AliTOFGeometry *geom):
  TObject(),
  fNSector(0),
  fNPlate(0),
  fNStripA(0),
  fNStripB(0),
  fNStripC(0),
  fNpadZ(0),
  fNpadX(0),
  fnpad(0),
  fGeom(geom),
  fPads(0x0)
{
  //ctor with geom
  fNSector = fGeom->NSectors();
  fNPlate  = fGeom->NPlates();
  fNStripA = fGeom->NStripA();
  fNStripB = fGeom->NStripB();
  fNStripC = fGeom->NStripC();
  fNpadZ = fGeom->NpadZ();
  fNpadX = fGeom->NpadX();
  fnpad = fNSector*(2*(fNStripC+fNStripB)+fNStripA)*fNpadZ*fNpadX;
}
//________________________________________________________________

AliTOFCalOnline::AliTOFCalOnline(const AliTOFCalOnline& cal):
  TObject(cal),
  fNSector(0),
  fNPlate(0),
  fNStripA(0),
  fNStripB(0),
  fNStripC(0),
  fNpadZ(0),
  fNpadX(0),
  fnpad(0),
  fGeom(0x0),
  fPads(0x0)  
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
  }
//____________________________________________________________________________ 
AliTOFCalOnline& AliTOFCalOnline::operator=(const AliTOFCalOnline& cal)
  {
    //assignment operator
    this->fNSector = cal.fNSector;
    this->fNPlate = cal.fNPlate;
    this->fNStripA = cal.fNStripA;
    this->fNStripB = cal.fNStripB;
    this->fNStripC = cal.fNStripC;
    this->fNpadZ = cal.fNpadZ;
    this->fNpadX = cal.fNpadX;
    this->fnpad = cal.fnpad;
    for (Int_t i = 0; i<fnpad; i++){
      this->fPads[i]=cal.fPads[i];
    }
    return *this;

  }
//____________________________________________________________________________ 
AliTOFCalOnline::~AliTOFCalOnline()
{
  //dtor
  delete [] fPads;
}
//________________________________________________________________

void AliTOFCalOnline::CreateArray(){
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
  fPads= new AliTOFChannelOnline[fnpad];
}
