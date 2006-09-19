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
***************************************************************************/

/*
$Log$
Revision 1.10  2006/08/22 13:26:05  arcelli
removal of effective c++ warnings (C.Zampolli)

Revision 1.9  2006/08/10 14:46:54  decaro
TOF raw data format: updated version

Revision 1.8  2006/05/04 19:41:42  hristov
Possibility for partial TOF geometry (S.Arcelli)

Revision 1.7  2006/04/27 13:13:29  hristov
Moving the destructor to the implementation file

Revision 1.6  2006/04/20 22:30:49  hristov
Coding conventions (Annalisa)

Revision 1.5  2006/04/16 22:29:05  hristov
Coding conventions (Annalisa)

Revision 1.4  2006/04/05 08:35:38  hristov
Coding conventions (S.Arcelli, C.Zampolli)

Revision 1.3  2006/03/31 13:49:07  arcelli
Removing some junk printout

Revision 1.2  2006/03/31 11:26:30  arcelli
 changing CDB Ids according to standard convention

Revision 1.1  2006/03/28 14:54:48  arcelli
class for TOF alignment

author: Silvia Arcelli, arcelli@bo.infn.it
*/  

/////////////////////////////////////////////////////////
//                                                     //
//            Class for alignment procedure            //
//                                                     //
//                                                     //
//                                                     //
/////////////////////////////////////////////////////////

#include <Rtypes.h>

#include "TRandom.h"

#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliAlignObjAngles.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliTOFAlignment.h"

ClassImp(AliTOFAlignment)

//_____________________________________________________________________________
AliTOFAlignment::AliTOFAlignment():
  TTask("AliTOFAlignment",""),
  fNTOFAlignObj(0),
  fTOFAlignObjArray(0x0)
 { 
  //AliTOFalignment main Ctor

}
//_____________________________________________________________________________
AliTOFAlignment::AliTOFAlignment(const AliTOFAlignment &t):
  TTask("AliTOFAlignment",""),
  fNTOFAlignObj(0),
  fTOFAlignObjArray(0x0)
{ 
  //AliTOFAlignment copy Ctor

  fNTOFAlignObj=t.fNTOFAlignObj;
  fTOFAlignObjArray=t.fTOFAlignObjArray;

}

//_____________________________________________________________________________
AliTOFAlignment& AliTOFAlignment::operator=(const AliTOFAlignment &t){ 
  //AliTOFAlignment assignment operator

  this->fNTOFAlignObj=t.fNTOFAlignObj;
  this->fTOFAlignObjArray=t.fTOFAlignObjArray;
  return *this;

}

//_____________________________________________________________________________
AliTOFAlignment::~AliTOFAlignment() {delete fTOFAlignObjArray;}

//_____________________________________________________________________________
void AliTOFAlignment::Smear( Float_t *tr, Float_t *rot)
{
  //Introduce Random Offset/Tilts
  fTOFAlignObjArray = new TObjArray(kMaxAlignObj);
  Float_t dx, dy, dz;  // shifts
  Float_t dpsi, dtheta, dphi; // angular displacements
  TRandom *rnd   = new TRandom(1567);
 
  Int_t nSMTOF = 18;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t iIndex=0; //dummy volume index
  //  AliAlignObj::ELayerID iLayer = AliAlignObj::kTOF;
  //  Int_t iIndex=1; //dummy volume index
  UShort_t dvoluid = AliAlignObj::LayerToVolUID(iLayer,iIndex); //dummy volume identity 
  Int_t i;
  for (i = 0; i<nSMTOF ; i++) {
    Char_t  path[100];
    sprintf(path,"/ALIC_1/B077_1/BSEGMO%i_1/BTOF%i_1",i,i);

    dx = (rnd->Gaus(0.,1.))*tr[0];
    dy = (rnd->Gaus(0.,1.))*tr[1];
    dz = (rnd->Gaus(0.,1.))*tr[2];
    dpsi   = rot[0];
    dtheta = rot[1];
    dphi   = rot[2];
    AliAlignObjAngles *o =new AliAlignObjAngles(path, dvoluid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    fTOFAlignObjArray->Add(o);
  }

  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  delete rnd;
}

//_____________________________________________________________________________
void AliTOFAlignment::Align( Float_t *tr, Float_t *rot)
{
  //Introduce Offset/Tilts

  fTOFAlignObjArray = new TObjArray(kMaxAlignObj);
  Float_t dx, dy, dz;  // shifts
  Float_t dpsi, dtheta, dphi; // angular displacements


  Int_t nSMTOF = 18;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t iIndex=0; //dummy volume index
  //  AliAlignObj::ELayerID iLayer = AliAlignObj::kTOF;
  //  Int_t iIndex=1; //dummy volume index
  UShort_t dvoluid = AliAlignObj::LayerToVolUID(iLayer,iIndex); //dummy volume identity 
  Int_t i;
  for (i = 0; i<nSMTOF ; i++) {

    Char_t  path[100];
    sprintf(path,"/ALIC_1/B077_1/BSEGMO%i_1/BTOF%i_1",i,i);
    dx = tr[0];
    dy = tr[1];
    dz = tr[2];
    dpsi   = rot[0];
    dtheta = rot[1];
    dphi   = rot[2];
    
    AliAlignObjAngles *o =new AliAlignObjAngles(path, dvoluid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    fTOFAlignObjArray->Add(o);
  }
  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
}
//_____________________________________________________________________________
void AliTOFAlignment::WriteParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write Align Par on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  Char_t *sel1 = "AlignPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBId idTOFAlign(out,minrun,maxrun);
  AliCDBMetaData *mdTOFAlign = new AliCDBMetaData();
  mdTOFAlign->SetResponsible("TOF");
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  man->Put(fTOFAlignObjArray,idTOFAlign,mdTOFAlign);
}
//_____________________________________________________________________________
void AliTOFAlignment::ReadParFromCDB(Char_t *sel, Int_t nrun)
{
  //Read Align Par from CDB
  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  Char_t *sel1 = "AlignPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBEntry *entry = man->Get(out,nrun);
  fTOFAlignObjArray=(TObjArray*)entry->GetObject();
  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes from CDB: %d",fNTOFAlignObj));

}
//_____________________________________________________________________________
void AliTOFAlignment::WriteSimParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write Sim Align Par on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  Char_t *sel1 = "AlignSimPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBId idTOFAlign(out,minrun,maxrun);
  AliCDBMetaData *mdTOFAlign = new AliCDBMetaData();
  mdTOFAlign->SetResponsible("TOF");
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  man->Put(fTOFAlignObjArray,idTOFAlign,mdTOFAlign);
}
//_____________________________________________________________________________
void AliTOFAlignment::ReadSimParFromCDB(Char_t *sel, Int_t nrun){
  //Read Sim Align Par from CDB
  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  Char_t *sel1 = "AlignSimPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBEntry *entry = man->Get(out,nrun);
  fTOFAlignObjArray=(TObjArray*)entry->GetObject();
  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes from CDB: %d",fNTOFAlignObj));

}
//_____________________________________________________________________________
void AliTOFAlignment::WriteOnCDBforDC()
{
  //Write Align Par on CDB for DC06
  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBId idTOFAlign("TOF/Align/Data",0,0);
  AliCDBMetaData *mdTOFAlign = new AliCDBMetaData();
  mdTOFAlign->SetComment("Alignment objects for ideal geometry, i.e. applying them to TGeo has to leave geometry unchanged");
  mdTOFAlign->SetResponsible("TOF");
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  man->Put(fTOFAlignObjArray,idTOFAlign,mdTOFAlign);
}
//_____________________________________________________________________________
void AliTOFAlignment::ReadFromCDBforDC()
{
  //Read Sim Align Par from CDB for DC06
  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBEntry *entry = man->Get("TOF/Align/Data",0);
  fTOFAlignObjArray=(TObjArray*)entry->GetObject();
  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes from CDB: %d",fNTOFAlignObj));

}
