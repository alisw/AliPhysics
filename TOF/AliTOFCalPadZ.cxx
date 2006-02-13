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
// class for TOF calibration : PadZ                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TClass.h"
#include "AliTOFGeometryV4.h"
#include "AliTOFChannel.h"
#include "AliTOFCalPadZ.h"

ClassImp(AliTOFCalPadZ)

//________________________________________________________________

AliTOFCalPadZ::AliTOFCalPadZ(){
  fCh = 0;
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

AliTOFCalPadZ::AliTOFCalPadZ(AliTOFChannel *ch):
  fCh(ch)
{  
  fNSector = AliTOFGeometryV4::NSectors();
  fNPlate = AliTOFGeometryV4::NPlates();
  fNStripA = AliTOFGeometryV4::NStripA();
  fNStripB = AliTOFGeometryV4::NStripB();
  //  fNStripC = AliTOFGeometryV4::NStripC();
  fNStripC = 20;
  fNpadZ = AliTOFGeometryV4::NpadZ();
  fNpadX = AliTOFGeometryV4::NpadX();

}
//________________________________________________________________

AliTOFCalPadZ::~AliTOFCalPadZ()
  {
    delete[] fCh;
  }
//________________________________________________________________

void AliTOFCalPadZ::Browse(TBrowser *b)
{
  char name[10];
  for(Int_t i=0; i<fNpadX; ++i) {
    snprintf(name,sizeof(name),"PadX %2.2d",i);
    b->Add(new AliTOFChannel(fCh[i]),name);
  }
}
