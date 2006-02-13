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
#include "AliTOFGeometryV4.h"
#include "AliTOFCalPlateA.h"
#include "AliTOFCalPlateB.h"
#include "AliTOFCalPlateC.h"
#include "AliTOFCalSector.h"
#include "AliTOFChannel.h"

extern TROOT *gROOT;

ClassImp(AliTOFCalSector)

//________________________________________________________________

AliTOFCalSector::AliTOFCalSector(){
  fCh = 0;
  fNSector = AliTOFGeometryV4::NSectors();
  fNPlate = AliTOFGeometryV4::NPlates();
  fNStripA = AliTOFGeometryV4::NStripA();
  fNStripB = AliTOFGeometryV4::NStripB();
  fNStripC = 20;
  //  fNStripC = AliTOFGeometryV4::NStripC();
  fNpadZ = AliTOFGeometryV4::NpadZ();
  fNpadX = AliTOFGeometryV4::NpadX();
  gROOT->GetListOfBrowsables()->Add(this);

}
//________________________________________________________________

AliTOFCalSector::AliTOFCalSector(AliTOFChannel *ch):
  fCh(ch)
{
  fNSector = AliTOFGeometryV4::NSectors();
  fNPlate = AliTOFGeometryV4::NPlates();
  fNStripA = AliTOFGeometryV4::NStripA();
  fNStripB = AliTOFGeometryV4::NStripB();
  fNStripC = 20;
  //  fNStripC = AliTOFGeometryV4::NStripC();
  fNpadZ = AliTOFGeometryV4::NpadZ();
  fNpadX = AliTOFGeometryV4::NpadX();
}
//________________________________________________________________

AliTOFCalSector::AliTOFCalSector(const AliTOFCalSector& sec):
  TObject(sec)
  {
    fCh = sec.fCh;
    fNSector = sec.fNSector;
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
  delete[] fCh;
}

//________________________________________________________________

void AliTOFCalSector::Browse(TBrowser *b){

  b->Add(new AliTOFCalPlateC(fCh),        "Plate0");
  b->Add(new AliTOFCalPlateB(&fCh[fNStripC*96]),"Plate1");
  b->Add(new AliTOFCalPlateA(&fCh[(fNStripC+fNStripB)*fNpadZ*fNpadX]),"Plate2");
  b->Add(new AliTOFCalPlateB(&fCh[(fNStripC+2*fNStripB)*fNpadZ*fNpadX]),"Plate3");
  b->Add(new AliTOFCalPlateC(&fCh[2*(fNStripC+fNStripB)*fNpadZ*fNpadX]),"Plate4");
}
 
