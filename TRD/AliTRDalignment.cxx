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

/* $Id$ */

//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// An AliTRDalignment object contains the alignment data (3 shifts and 3 tilts) //
// for all the alignable volumes of the TRD, i.e. for 18 supermodules and 540   //
// chambers. The class provides simple tools for reading and writing these data //
// in different formats, and for generating fake data that can be used to       //
// simulate misalignment.                                                       //
// The six alignment variables have the following meaning:                      //
// shift in rphi                                                                //
// shift in z                                                                   //
// shift in r                                                                   //
// tilt around rphi                                                             //
// tilt around z                                                                //
// tilt around r                                                                //
// The shifts are in cm and the tilts are in degrees.                           //
// The currently supported formats are:                                         //
// - ascii                                                                      //
// - root file containing a TClonesArray of alignment objects                   //
// - offline conditions database                                                //
// - OCDB-like root file                                                        //
// - geometry file (like misaligned_geometry.root)                              //
//                                                                              //
// D.Miskowiec, November 2006                                                   //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <string>

#include "TMath.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoPhysicalNode.h"
#include "TClonesArray.h"

#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliAlignObjAngles.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "AliCDBId.h"

#include "AliTRDalignment.h"

ClassImp(AliTRDalignment)

//_____________________________________________________________________________
AliTRDalignment::AliTRDalignment() 
  :TObject()
  ,fRan(0)
{
  //
  // constructor
  //

  SetZero();

}

//_____________________________________________________________________________
AliTRDalignment::AliTRDalignment(const AliTRDalignment& source) 
  :TObject(source)
  ,fRan(source.fRan) 
{
  //
  // copy constructor
  //

  for (int i = 0; i <  18; i++) {
    SetSm(i,source.fSm[i]);
  }
  for (int i = 0; i < 540; i++) {
    SetCh(i,source.fCh[i]);
  }

}

//_____________________________________________________________________________
AliTRDalignment& AliTRDalignment::operator=(const AliTRDalignment &source) 
{
  //
  // assignment operator
  //

  if (this != &source) {
    for (int i = 0; i <  18; i++) SetSm(i,source.fSm[i]);
    for (int i = 0; i < 540; i++) SetCh(i,source.fCh[i]);
  }

  return *this;

}

//_____________________________________________________________________________
AliTRDalignment& AliTRDalignment::operator+=(const AliTRDalignment &source) 
{
  //
  // addition operator
  //

  for (int i = 0; i <  18; i++) {
    for (int j = 0; j < 6; j++) {
      this->fSm[i][j] =+ source.fSm[i][j];
    }
  }
  for (int i = 0; i < 540; i++) {
    for (int j = 0; j < 6; j++) {
      this->fCh[i][j] =+ source.fCh[i][j];
    }
  }

  return *this;

}

//_____________________________________________________________________________
AliTRDalignment& AliTRDalignment::operator-=(const AliTRDalignment &source) 
{
  //
  // subtraction operator
  //

  for (int i = 0; i <  18; i++) {
    for (int j = 0; j < 6; j++) {
      fSm[i][j] -= source.fSm[i][j];
    }
  }
  for (int i = 0; i < 540; i++) {
    for (int j = 0; j < 6; j++) {
      fCh[i][j] -= source.fCh[i][j];
    }
  }

  return *this;

}

//_____________________________________________________________________________
Bool_t AliTRDalignment::operator==(const AliTRDalignment &source) const
{
  //
  // comparison operator
  //

  Bool_t areEqual = 1;

  for (int i = 0; i <  18; i++) {
    for (int j = 0; j < 6; j++) {
      areEqual &= (fSm[i][j] == source.fSm[i][j]);
    }
  }
  for (int i = 0; i < 540; i++) {
    for (int j = 0; j < 6; j++) {
      areEqual &= (fCh[i][j] == source.fCh[i][j]);
    }
  }

  return areEqual;

}

//_____________________________________________________________________________
void AliTRDalignment::SetSmZero() 
{
  //
  // reset to zero supermodule data
  //

  memset(&fSm[0][0],0,sizeof(fSm));

}

//_____________________________________________________________________________
void AliTRDalignment::SetChZero() 
{
  //
  // reset to zero chamber data
  //

  memset(&fCh[0][0],0,sizeof(fCh));

}

//_____________________________________________________________________________
void AliTRDalignment::SetSmRandom(Double_t a[6]) 
{
  //
  // generate random gaussian supermodule data with sigmas a
  //

  double x[6];

  for (Int_t i = 0; i < 18; i++) {
    fRan.Rannor(x[0],x[1]);
    fRan.Rannor(x[2],x[3]);
    fRan.Rannor(x[4],x[5]);
    for (Int_t j = 0; j < 6; j++) {
      x[j] *= a[j];
    }
    SetSm(i,x);
    //PrintSm(i);
  }

}

//_____________________________________________________________________________
void AliTRDalignment::SetChRandom(Double_t a[6]) 
{
  //
  // generate random gaussian chamber data with sigmas a
  //

  double x[6];

  for (Int_t i = 0; i < 540; i++) {
    fRan.Rannor(x[0],x[1]);
    fRan.Rannor(x[2],x[3]);
    fRan.Rannor(x[4],x[5]);
    for (Int_t j = 0; j < 6; j++) {
      x[j] *= a[j];
    }
    SetCh(i,x);
    //PrintCh(i);
  }

}

//_____________________________________________________________________________
void AliTRDalignment::SetSmFull() 
{
  //
  // generate random gaussian supermodule data similar to the misalignment 
  // expected from the mechanical precision 
  //

  Double_t a[6];

  a[0] = 0.3; // phi
  a[1] = 0.3; // z
  a[2] = 0.3; // r
  a[3] = 0.4/1000.0 / TMath::Pi()*180.0; // phi
  a[4] = 2.0/1000.0 / TMath::Pi()*180.0; // z
  a[5] = 0.4/1000.0 / TMath::Pi()*180.0; // r

  SetSmRandom(a);

}

//_____________________________________________________________________________
void AliTRDalignment::SetChFull() 
{
  //
  // generate random gaussian chamber data similar to the misalignment 
  // expected from the mechanical precision 
  //

  Double_t a[6];

  a[0] = 0.1; // phi
  a[1] = 0.1; // z
  a[2] = 0.1; // r
  a[3] = 1.0/1000.0 / TMath::Pi()*180.0; // phi
  a[4] = 1.0/1000.0 / TMath::Pi()*180.0; // z
  a[5] = 0.7/1000.0 / TMath::Pi()*180.0; // r

  SetChRandom(a);

}

//_____________________________________________________________________________
void AliTRDalignment::SetSmResidual() 
{
  //
  // generate random gaussian supermodule data similar to the misalignment 
  // remaining after full calibration
  // I assume that it will be negligible
  //

  SetSmZero();

}

//_____________________________________________________________________________
void AliTRDalignment::SetChResidual() 
{
  //
  // generate random gaussian chamber data similar to the misalignment 
  // remaining after full calibration
  //

  Double_t a[6];

  a[0] = 0.002; // phi
  a[1] = 0.003; // z
  a[2] = 0.007; // r
  a[3] = 0.3/1000.0 / TMath::Pi()*180.0; // phi
  a[4] = 0.3/1000.0 / TMath::Pi()*180.0; // z
  a[5] = 0.1/1000.0 / TMath::Pi()*180.0; // r

  SetChRandom(a);

}

//_____________________________________________________________________________
void AliTRDalignment::PrintSm(Int_t i, FILE *fp) const 
{
  //
  // print the supermodule data
  //

  fprintf(fp,"%4d   %11.4f %11.4f  %11.4f      %11.5f  %11.5f  %11.5f   %6d  %s\n"
	 ,i,fSm[i][0],fSm[i][1],fSm[i][2],fSm[i][3],fSm[i][4],fSm[i][5]
	 ,0,GetSmName(i));

}

//_____________________________________________________________________________
void AliTRDalignment::PrintCh(Int_t i, FILE *fp) const 
{
  //
  // print the chamber data
  //

  fprintf(fp,"%4d   %11.4f %11.4f  %11.4f      %11.5f  %11.5f  %11.5f   %6d  %s\n"
	 ,i,fCh[i][0],fCh[i][1],fCh[i][2],fCh[i][3],fCh[i][4],fCh[i][5]
	 ,GetVoi(i),GetChName(i));

}

//_____________________________________________________________________________
void AliTRDalignment::ReadAscii(char *filename) 
{
  //
  // read the alignment data from ascii file
  //

  double x[6];      // alignment data
  int volid;        // volume id
  std::string syna; // symbolic name
  int j;            // dummy index

  fstream fi(filename,fstream::in);
  if (!fi) {
    AliFatal(Form("cannot open input file %s",filename));
  }

  // supermodules

  for (int i = 0; i < 18; i++) {
    fi>>j>>x[0]>>x[1]>>x[2]>>x[3]>>x[4]>>x[5]>>volid>>syna;
    if (j != i) {
      AliError(Form("sm %d expected, %d found",i,j));
    }
    if (volid != 0) {
      AliError(Form("sm %d volume id %d expected, %d found",i,0,volid));
    }
    std::string symnam = GetSmName(i);
    if (syna != symnam) {
      AliError(Form("sm %d name %s expected, %s found",i,symnam.data(),syna.data()));
    }
    SetSm(i,x);
  }

  // chambers

  for (int i = 0; i < 540; i++) {
    fi>>j>>x[0]>>x[1]>>x[2]>>x[3]>>x[4]>>x[5]>>volid>>syna;
    if (j != i) {
      AliError(Form("ch %d expected, %d found",i,j));
    }
    if (volid != GetVoi(i)) {
      AliError(Form("ch %d volume id %d expected, %d found",i,GetVoi(i),volid));
    }
    std::string symnam = GetChName(i);
    if (syna != symnam) {
      AliError(Form("ch %d name %s expected, %s found",i,symnam.data(),syna.data()));
    }
    SetCh(i,x);
  }

  fi.close();

}

//_____________________________________________________________________________
void AliTRDalignment::ReadRoot(char *filename) 
{
  //
  // read the alignment data from root file
  // here I expect a fixed order and number of elements
  // it would be much better to identify the alignment objects 
  // one by one and set the parameters of the corresponding sm or ch
  //

  TFile fi(filename,"READ");

  if (fi.IsOpen()) {
    TClonesArray *ar = (TClonesArray*) fi.Get("TRDAlignObjs");
    ArToNumbers(ar);
    fi.Close();
  } 
  else {
    AliError(Form("cannot open input file %s",filename));
  }

  return;

}

//_____________________________________________________________________________
void AliTRDalignment::ReadDB(char *filename) 
{
  //
  // read the alignment data from database file
  //

  TFile fi(filename,"READ");

  if (fi.IsOpen()) {
    AliCDBEntry  *e  = (AliCDBEntry *) fi.Get("AliCDBEntry");
    e->PrintMetaData();
    TClonesArray *ar = (TClonesArray *) e->GetObject();
    ArToNumbers(ar);
    fi.Close();
  } 
  else {
    AliError(Form("cannot open input file %s",filename));
  }

  return;

}

//_____________________________________________________________________________
void AliTRDalignment::ReadDB(char *db, char *path, Int_t run
                           , Int_t version, Int_t subversion)
{
  //
  // read the alignment data from database
  //

  AliCDBManager *cdb     = AliCDBManager::Instance();
  AliCDBStorage *storLoc = cdb->GetStorage(db);
  AliCDBEntry   *e       = storLoc->Get(path,run,version,subversion);
  e->PrintMetaData();
  TClonesArray  *ar      =  (TClonesArray *) e->GetObject();
  ArToNumbers(ar);

}

//_____________________________________________________________________________
void AliTRDalignment::ReadGeo(char *misaligned) 
{
  //
  // determine misalignment by comparing original and misaligned matrix 
  // of the last node on the misaligned_geometry file 
  // an alternative longer way is in attic.C
  //

  TGeoHMatrix *ideSm[18];  // ideal
  TGeoHMatrix *ideCh[540];
  TGeoHMatrix *misSm[18];  // misaligned
  TGeoHMatrix *misCh[540];

  // read misaligned and original matrices

  TGeoManager::Import(misaligned);
  for (int i = 0; i < 18; i++) {
    TGeoPNEntry      *pne  = gGeoManager->GetAlignableEntry(GetSmName(i));
    if (!pne) {
      AliError(Form("no such physical node entry: %s",GetSmName(i))); 
      return;
    }
    TGeoPhysicalNode *node = pne->GetPhysicalNode();
    if (!node) {
      AliError(Form("physical node entry %s has no physical node",GetSmName(i))); 
      return;
    }
    misSm[i] = new TGeoHMatrix(*node->GetNode(node->GetLevel())->GetMatrix());
    ideSm[i] = new TGeoHMatrix(*node->GetOriginalMatrix());
  }
  for (int i = 0; i < 540; i++) {
    TGeoPNEntry      *pne  = gGeoManager->GetAlignableEntry(GetChName(i));
    if (!pne) {
      AliError(Form("no such physical node entry: %s",GetChName(i))); 
      return;
    }
    TGeoPhysicalNode *node = pne->GetPhysicalNode();
    if (!node) {
      AliError(Form("physical node entry %s has no physical node",GetChName(i))); 
      return;
    }
    misCh[i] = new TGeoHMatrix(*node->GetNode(node->GetLevel())->GetMatrix());
    ideCh[i] = new TGeoHMatrix(*node->GetOriginalMatrix());
  }

  // calculate the local misalignment matrices as inverse misaligned times ideal

  for (int i = 0; i < 18; i++) {
    TGeoHMatrix mat(ideSm[i]->Inverse()); 
    mat.Multiply(misSm[i]);
    double *tra = mat.GetTranslation();
    double *rot = mat.GetRotationMatrix();
    double pars[6];
    pars[0] = tra[0];
    pars[1] = tra[1];
    pars[2] = tra[2];
    if (TMath::Abs(rot[0])<1e-7 || TMath::Abs(rot[8])<1e-7) {
      AliError("Failed to extract roll-pitch-yall angles!");
      return;
    }
    double raddeg = TMath::RadToDeg();
    pars[3] = raddeg * TMath::ATan2(-rot[5],rot[8]);
    pars[4] = raddeg * TMath::ASin(rot[2]);
    pars[5] = raddeg * TMath::ATan2(-rot[1],rot[0]);
    SetSm(i,pars);
  }

  for (int i = 0; i < 540; i++) {
    TGeoHMatrix mat(ideCh[i]->Inverse()); 
    mat.Multiply(misCh[i]);
    double *tra = mat.GetTranslation();
    double *rot = mat.GetRotationMatrix();
    double pars[6];
    pars[0] = tra[0];
    pars[1] = tra[1];
    pars[2] = tra[2];
    if(TMath::Abs(rot[0])<1e-7 || TMath::Abs(rot[8])<1e-7) {
      AliError("Failed to extract roll-pitch-yall angles!");
      return;
    }
    double raddeg = TMath::RadToDeg();
    pars[3] = raddeg * TMath::ATan2(-rot[5],rot[8]);
    pars[4] = raddeg * TMath::ASin(rot[2]);
    pars[5] = raddeg * TMath::ATan2(-rot[1],rot[0]);
    SetCh(i,pars);
  }

  // cleanup
  for (int i = 0; i <  18; i++) delete ideSm[i];
  for (int i = 0; i <  18; i++) delete misSm[i];
  for (int i = 0; i < 540; i++) delete ideCh[i];
  for (int i = 0; i < 540; i++) delete misCh[i];

  return;

}

//_____________________________________________________________________________
void AliTRDalignment::ReadSurveyReport(char *filename) 
{
  // read survey report and set the supermodule parameters correspondingly

  fstream fi(filename,fstream::in);
  if (!fi) {
    AliFatal(Form("cannot open input file %s",filename));
  }

  // to be continued...

}

//_____________________________________________________________________________
void AliTRDalignment::ReadAny(char *filename) 
{
  //
  // read the alignment data from any kind of file
  //

  TString fist(filename);
  if (fist.EndsWith(".txt")) {
    ReadAscii(filename);
  }
  if (fist.EndsWith(".dat")) {
    ReadAscii(filename);
  }
  if (fist.EndsWith(".root")) {
    if (fist.Contains("Run")) {
      ReadDB(filename);
    }
    else {
      ReadRoot(filename);
    }
  }

}

//_____________________________________________________________________________
void AliTRDalignment::WriteAscii(char *filename) const
{
  //
  // store the alignment data on ascii file
  //

  FILE *fp = fopen(filename, "w");
  if (!fp) {
    AliError(Form("cannot open output file %s",filename));
    return;
  }

  PrintSm(fp);
  PrintCh(fp);
  
  fclose(fp);

}

//_____________________________________________________________________________
void AliTRDalignment::WriteRoot(char *filename) 
{
  //
  // store the alignment data on root file
  //

  TClonesArray *ar = new TClonesArray("AliAlignObjAngles",10000);
  NumbersToAr(ar);
  TFile fo(filename,"RECREATE");
  if (fo.IsOpen()) {
    fo.cd();
    fo.WriteObject(ar,"TRDAlignObjs","kSingleKey");
    fo.Close();
  } 
  else {
    AliError(Form("cannot open output file %s",filename));
  }

  delete ar;

}

//_____________________________________________________________________________
void AliTRDalignment::WriteDB(char *filename, char *comment, Int_t run0, Int_t run1) 
{
  //
  // dumping on a DB-like file
  //

  TClonesArray   *ar = new TClonesArray("AliAlignObjAngles",10000);
  NumbersToAr(ar);
  char *path = "di1/di2/di3";
  AliCDBId id(path,run0,run1);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Dariusz Miskowiec");
  md->SetComment(comment);
  AliCDBEntry    *e  = new AliCDBEntry(ar, id, md);
  TFile fi(filename,"RECREATE");
  if (fi.IsOpen()) {
    e->Write();
    fi.Close();
  } 
  else {
    AliError(Form("cannot open input file %s",filename));
  }

  delete e;
  delete md;
  delete ar;

  return;

}

//_____________________________________________________________________________
void AliTRDalignment::WriteDB(char *db, char *path, char *comment, Int_t run0, Int_t run1) 
{
  //
  // store the alignment data in database
  //

  TClonesArray   *ar      = new TClonesArray("AliAlignObjAngles",10000);
  NumbersToAr(ar);
  AliCDBManager  *cdb     = AliCDBManager::Instance();
  AliCDBStorage  *storLoc = cdb->GetStorage(db);
  AliCDBMetaData *md      = new AliCDBMetaData();
  md->SetResponsible("Dariusz Miskowiec");
  md->SetComment(comment);
  AliCDBId id(path,run0,run1);
  storLoc->Put(ar,id,md);
  md->Delete();
  delete ar;

}

//_____________________________________________________________________________
void AliTRDalignment::WriteGeo(char *filename) 
{
  //
  // apply misalignment to (currently loaded ideal) geometry and store the 
  // resulting geometry on a root file
  //

  TClonesArray *ar = new TClonesArray("AliAlignObjAngles",10000);
  NumbersToAr(ar);
  for (int i = 0; i < ar->GetEntriesFast(); i++) {
    AliAlignObj *alobj = (AliAlignObj *) ar->UncheckedAt(i);
    alobj->ApplyToGeometry();
  }
  delete ar;
  gGeoManager->Export(filename);

}

//_____________________________________________________________________________
Double_t AliTRDalignment::GetSmRMS(Int_t xyz) const 
{
  //
  // rms fSm[][xyz]
  //

  Double_t s1 = 0.0;
  Double_t s2 = 0.0;
  for (int i = 0; i < 18; i++) {
    s1 += fSm[i][xyz];
    s2 += fSm[i][xyz]*fSm[i][xyz];
  }
  Double_t rms2 = s2/18.0 - s1*s1/18.0/18.0;

  return rms2>0 ? sqrt(rms2) : 0.0;

}

//_____________________________________________________________________________
Double_t AliTRDalignment::GetChRMS(Int_t xyz) const
{
  //
  // rms fCh[][xyz]
  //

  Double_t s1 =0.0;
  Double_t s2 =0.0;
  for (int i = 0; i < 540; i++) {
    s1 += fCh[i][xyz];
    s2 += fCh[i][xyz]*fCh[i][xyz];
  }
  Double_t rms2 = s2/540.0 - s1*s1/540.0/540.0;

  return rms2>0 ? sqrt(rms2) : 0.0;

}

//_____________________________________________________________________________
void AliTRDalignment::PrintSmRMS() const
{
  //
  // dump rms of fSm
  //

  printf("       %11.4f %11.4f  %11.4f      %11.5f  %11.5f  %11.5f  supermodule rms\n"
	,GetSmRMS(0),GetSmRMS(1),GetSmRMS(2),GetSmRMS(3),GetSmRMS(4),GetSmRMS(5));

}

//_____________________________________________________________________________
void AliTRDalignment::PrintChRMS() const
{
  //
  // dump rms of fCh
  //

  printf("       %11.4f %11.4f  %11.4f      %11.5f  %11.5f  %11.5f  chamber rms\n"
	,GetChRMS(0),GetChRMS(1),GetChRMS(2),GetChRMS(3),GetChRMS(4),GetChRMS(5));

}

//_____________________________________________________________________________
void AliTRDalignment::ArToNumbers(TClonesArray *ar) 
{
  //
  // read numbers from the array of AliAlignObj objects and fill fSm and fCh
  //

  LoadIdealGeometry();
  SetZero();
  double pa[6];
  for (int i = 0; i <  18; i++) {
    AliAlignObj *aao = (AliAlignObj *) ar->At(i);
    aao->GetLocalPars(pa,pa+3);
    SetSm(i,pa);
  }
  for (int i = 0; i < 540; i++) {
    AliAlignObj *aao = (AliAlignObj *) ar->At(18+i);
    aao->GetLocalPars(pa,pa+3);
    SetCh(i,pa);
  }

}

//_____________________________________________________________________________
void AliTRDalignment::NumbersToAr(TClonesArray *ar) 
{
  //
  // build array of AliAlignObj objects based on fSm and fCh data
  //

  LoadIdealGeometry();
  TClonesArray &alobj = *ar;
  int nobj = 0;
  for (int i = 0; i <  18; i++) {      
      new(alobj[nobj]) AliAlignObjAngles(GetSmName(i)
                                        ,0 
					,fSm[i][0],fSm[i][1],fSm[i][2]
					,fSm[i][3],fSm[i][4],fSm[i][5]
					,0);
    nobj++;
  }

  for (int i = 0; i < 540; i++) {
    new(alobj[nobj]) AliAlignObjAngles(GetChName(i)
                                      ,GetVoi(i)
				      ,fCh[i][0],fCh[i][1],fCh[i][2]
				      ,fCh[i][3],fCh[i][4],fCh[i][5]
				      ,0);
    nobj++;
  }

}

//_____________________________________________________________________________
void AliTRDalignment::LoadIdealGeometry(char *filename) 
{
  //
  // load ideal geometry from filename
  // it is needed for operations on AliAlignObj objects
  // this needs to be straightened out
  // particularly, sequences LoadIdealGeometry("file1"); LoadIdealGeometry("file2"); 
  // do not work as one would naturally expect
  //

  static int attempt = 0; // which reload attempt is it? just to avoid endless loops

  if (!gGeoManager) {
    TGeoManager::Import(filename);
  }
  if (!gGeoManager) {
    AliFatal(Form("cannot open geometry file %s",filename));
  }
  if (gGeoManager->GetListOfPhysicalNodes()->GetEntries()) {
    if (attempt) {
      AliFatal(Form("geometry on file %s is not ideal",filename));
    }
    AliWarning("current geometry is not ideal - it contains physical nodes");
    AliWarning(Form("reloading geometry from %s - attempt nr %d",filename,attempt));
    gGeoManager = 0;
    attempt++;
    LoadIdealGeometry(filename);
  }

  attempt = 0;

}
