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
Revision 1.3  2006/03/28 14:58:16  arcelli
updates to handle new V5 geometry & some re-arrangements

Revision 1.2  2006/02/13 16:53:00  decaro
just Fixing Log info

Revision 1.1  2006/02/13 16:10:48  arcelli
Add classes for TOF Calibration (C.Zampolli)

author: Chiara Zampolli, zampolli@bo.infn.it
*/  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF calibration : Sectors                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TClass.h"
#include "AliLog.h"
#include "AliTOFGeometryV5.h"
#include "AliTOFCalPlateA.h"
#include "AliTOFCalPlateB.h"
#include "AliTOFCalPlateC.h"
#include "AliTOFCalSector.h"
#include "AliTOFChannel.h"

extern TROOT *gROOT;

ClassImp(AliTOFCalSector)

//________________________________________________________________

AliTOFCalSector::AliTOFCalSector(){
  //main ctor
  fCh = 0;
  fGeom=0x0;
  fNPlate=0;
  fNStripA=0;
  fNStripB=0;
  fNStripC=0;
  fNpadZ=0;
  fNpadX=0;
  gROOT->GetListOfBrowsables()->Add(this);

}
//________________________________________________________________

AliTOFCalSector::AliTOFCalSector(AliTOFChannel *ch):
  fCh(ch)
{
  //ctor with channel
  fGeom=0x0;
  fNPlate=0;
  fNStripA=0;
  fNStripB=0;
  fNStripC=0;
  fNpadZ=0;
  fNpadX=0;
  gROOT->GetListOfBrowsables()->Add(this);
}
//________________________________________________________________

AliTOFCalSector::AliTOFCalSector(AliTOFGeometry *geom){
  //ctor with geom
  fCh = 0;
  fGeom= geom; 
  fNPlate  = fGeom->NPlates();
  fNStripA = fGeom->NStripA();
  fNStripB = fGeom->NStripB();
  fNStripC = fGeom->NStripC();
  fNpadZ = fGeom->NpadZ();
  fNpadX = fGeom->NpadX();
  gROOT->GetListOfBrowsables()->Add(this);

}
//________________________________________________________________

AliTOFCalSector::AliTOFCalSector(AliTOFGeometry *geom,AliTOFChannel *ch):
  fCh(ch)
{
  // ctor with channel and geom
  fGeom= geom; 
  fNPlate  = fGeom->NPlates();
  fNStripA = fGeom->NStripA();
  fNStripB = fGeom->NStripB();
  fNStripC = fGeom->NStripC();
  fNpadZ = fGeom->NpadZ();
  fNpadX = fGeom->NpadX();
  gROOT->GetListOfBrowsables()->Add(this);
}
//________________________________________________________________

AliTOFCalSector::AliTOFCalSector(const AliTOFCalSector& sec):
  TObject(sec)
  {
    //copy ctor
    fCh = sec.fCh;
    fNPlate = sec.fNPlate;
    fNStripA = sec.fNStripA;
    fNStripB = sec.fNStripB;
    fNStripC = sec.fNStripC;
    fNpadZ = sec.fNpadZ;
    fNpadX = sec.fNpadX;
    gROOT->GetListOfBrowsables()->Add(this);
  }
//________________________________________________________________

AliTOFCalSector::~AliTOFCalSector()
{
  //dtor
  gROOT->GetListOfBrowsables()->Remove(this);
  delete[] fCh;
}

//________________________________________________________________

void AliTOFCalSector::Browse(TBrowser *b){
  //add cal obj to list of browsables
  if(fGeom==0x0){
    AliTOFGeometry *geom= new AliTOFGeometryV5(); 
    AliInfo("V5 TOF Geometry is taken as the default");
    fNPlate  = geom->NPlates();
    fNStripA = geom->NStripA();
    fNStripB = geom->NStripB();
    fNStripC = geom->NStripC();
    fNpadZ = geom->NpadZ();
    fNpadX = geom->NpadX();
    delete geom;
  }
  b->Add(new AliTOFCalPlateC(fCh),        "Plate0");
  b->Add(new AliTOFCalPlateB(&fCh[fNStripC*96]),"Plate1");
  b->Add(new AliTOFCalPlateA(&fCh[(fNStripC+fNStripB)*fNpadZ*fNpadX]),"Plate2");
  b->Add(new AliTOFCalPlateB(&fCh[(fNStripC+2*fNStripB)*fNpadZ*fNpadX]),"Plate3");
  b->Add(new AliTOFCalPlateC(&fCh[2*(fNStripC+fNStripB)*fNpadZ*fNpadX]),"Plate4");
}
 
