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
Revision 1.2  2006/02/13 17:22:26  arcelli
just Fixing Log info

Revision 1.1  2006/02/13 16:10:48  arcelli
Add classes for TOF Calibration (C.Zampolli)

author: Chiara Zampolli, zampolli@bo.infn.it
*/  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF calibration : PlateB                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TClass.h"
#include "AliLog.h"
#include "AliTOFGeometryV5.h"
#include "AliTOFCalStrip.h"
#include "AliTOFCalPlateB.h"
#include "AliTOFChannel.h"

ClassImp(AliTOFCalPlateB)

//________________________________________________________________

AliTOFCalPlateB::AliTOFCalPlateB(){
  fCh = 0;
  fGeom= 0x0; 
  fNStripB = 0;
  fNpadZ = 0;
  fNpadX = 0;
}
//________________________________________________________________

AliTOFCalPlateB::AliTOFCalPlateB(AliTOFChannel *ch) : fCh(ch)
{
  fGeom= 0x0; 
  fNStripB = 0;
  fNpadZ = 0;
  fNpadX = 0;
}
//________________________________________________________________

AliTOFCalPlateB::AliTOFCalPlateB(AliTOFGeometry *geom){
  fCh = 0;
  fGeom = geom;  
  fNStripB = fGeom->NStripB();
  fNpadZ = fGeom->NpadZ();
  fNpadX = fGeom->NpadX();
}
//________________________________________________________________

AliTOFCalPlateB::AliTOFCalPlateB(AliTOFGeometry *geom, AliTOFChannel *ch): fCh(ch)
{
  fGeom = geom;  
  fNStripB = fGeom->NStripB();
  fNpadZ = fGeom->NpadZ();
  fNpadX = fGeom->NpadX();
}

//________________________________________________________________

AliTOFCalPlateB::~AliTOFCalPlateB()
{
  delete[] fCh;
}

//________________________________________________________________

AliTOFCalPlateB::AliTOFCalPlateB(const AliTOFCalPlateB& pl):
  TObject(pl)
  {
    fCh = pl.fCh;
    fNStripB = pl.fNStripB;
    fNpadZ = pl.fNpadZ;
    fNpadX = pl.fNpadX;
    fGeom = pl.fGeom;

  }
//________________________________________________________________

void AliTOFCalPlateB::Browse(TBrowser *b){

  if(fGeom==0x0){
    AliTOFGeometry *geom = new AliTOFGeometryV5(); 
    AliInfo("V5 TOF Geometry is taken as the default");
    fNStripB = geom->NStripB();
    fNpadZ = geom->NpadZ();
    fNpadX = geom->NpadX();
    delete geom;
  }
  char name[10];
  for(Int_t i=0; i<fNStripB; ++i) {
    snprintf(name,sizeof(name),"Strip %2.2d",i);
    b->Add(new AliTOFCalStrip(&fCh[i*fNpadZ*fNpadX]),name);
  }
}
