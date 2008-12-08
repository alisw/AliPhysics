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
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// An AliTRDalignment object contains the alignment data (3 shifts and 3     //
// tilts) for all the alignable volumes of the TRD, i.e. for 18 supermodules //
// and 540 chambers. The class provides simple tools for reading and writing //
// these data in different formats, and for generating fake data that can be //
// used to simulate misalignment.                                            //
// The six alignment variables have the following meaning:                   //
// shift in rphi                                                             //
// shift in z                                                                //
// shift in r                                                                //
// tilt around rphi                                                          //
// tilt around z                                                             //
// tilt around r                                                             //
// The shifts are in cm and the tilts are in degrees.                        //
// The currently supported formats are:                                      //
// - ascii                                                                   //
// - root file containing a TClonesArray of alignment objects                //
// - offline conditions database                                             //
// - OCDB-like root file                                                     //
// - geometry file (like misaligned_geometry.root)                           //
//                                                                           //
// Some examples of usage (in an aliroot session):                           //
// AliTRDalignment a,b,c,d,e;                                                //
// double xsm[]={0,0,0,-70,0,0};                                             //
// double xch[]={0,0,-50,0,0,0};                                             //
// a.SetSm(4,xsm);                                                           // 
// a.SetCh(120,xch);                                                         //
// a.WriteAscii("kuku.dat");                                                 //
// TGeoManager::Import("geometry.root"); a.WriteRoot("kuku.root");           //
// TGeoManager::Import("geometry.root"); a.WriteDB("kukudb.root",0,0);       //
// TGeoManager::Import("geometry.root");                                     //
// a.WriteDB("local://$ALICE_ROOT", "TRD/Align/Data", 0,0);                  //
// TGeoManager::Import("geometry.root"); a.WriteGeo("kukugeometry.root");    //
//                                                                           //
// b.ReadAscii("kuku.dat");                                                  //
// TGeoManager::Import("geometry.root"); c.ReadRoot("kuku.root");            //
// TGeoManager::Import("geometry.root"); d.ReadDB("kukudb.root");            //
// TGeoManager::Import("kukugeometry.root"); e.ReadCurrentGeo();             //
//                                                                           //
// e.PrintSm(4);                                                             //
// e.PrintCh(120);                                                           // 
// a.PrintRMS();                                                             //
// b.PrintRMS();                                                             //
// e.PrintRMS();                                                             //
//                                                                           //
//                                                                           //
// D.Miskowiec, November 2006                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
//#include <string>

#include "TMath.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoPhysicalNode.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TFitter.h"
#include "TMinuit.h"

#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliAlignObjParams.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "AliCDBId.h"
#include "AliSurveyObj.h"
#include "AliSurveyPoint.h"

#include "AliTRDalignment.h"

void trdAlignmentFcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *x, Int_t iflag);

ClassImp(AliTRDalignment)

//_____________________________________________________________________________
AliTRDalignment::AliTRDalignment() 
  :TObject()
  ,fComment()
  ,fRan(0)
{
  //
  // constructor
  //

  SetZero();

  for (int i=0; i<18; i++) for (int j=0; j<2; j++) for (int k=0; k<2; k++) for (int l=0; l<2; l++) {
    fSurveyX[i][j][k][l] = 0.0;
    fSurveyY[i][j][k][l] = 0.0;
    fSurveyZ[i][j][k][l] = 0.0;
    fSurveyEX[i][j][k][l] = 0.0;
    fSurveyEY[i][j][k][l] = 0.0;
    fSurveyEZ[i][j][k][l] = 0.0;
  }

  // Initialize the nominal positions of the survey points 
  // in the local frame of supermodule (where y is the long side, 
  // z corresponds to the radius in lab, and x to the phi in lab).
  // Four survey marks are on each z-side of the supermodule. 
  //               A           B
  //           ----o-----------o----        x |
  //           \                   /          |
  //            \                 /           |
  //             \               /            |
  //              \             /             |
  //               ---o-----o---              -------------->
  //                  C     D                              y
  // 
  // For the purpose of this explanation lets define the origin such that 
  // the supermodule occupies 0 < x < 77.9 cm. Then the coordinates (x,y) 
  // are (in cm) 
  // A (76.2,-30.25)
  // B (76.2,+30.25)
  // C ( 2.2,-22.5 )
  // D ( 2.2,+22.5 )
  // 

  double x[2] = {22.5,30.25};                   // lab phi, or tracking-y
  double y[2] = {353.0, -353.0};                // lab z; inc. 2 cm survey target offset
  double z[2] = {-(77.9/2.0-2.0),77.9/2.0-1.5}; // lab r, or better tracking-x

  for (int j=0; j<2; j++) for (int k=0; k<2; k++) for (int l=0; l<2; l++) {
    fSurveyX0[j][k][l] = -TMath::Power(-1,l) * x[k];
    fSurveyY0[j][k][l] = y[j];
    fSurveyZ0[j][k][l] = z[k];
  }

}

//_____________________________________________________________________________
AliTRDalignment::AliTRDalignment(const AliTRDalignment& source) 
  :TObject(source)
  ,fComment(source.fComment)
  ,fRan(source.fRan)
{
  //
  // copy constructor
  //

  for (int i=0; i<18; i++) SetSm(i,source.fSm[i]);
  for (int i=0; i<540; i++) SetCh(i,source.fCh[i]);
  for (int i=0; i<18; i++) for (int j=0; j<2; j++) for (int k=0; k<2; k++) for (int l=0; l<2; l++) {
    fSurveyX[i][j][k][l] = source.fSurveyX[i][j][k][l];
    fSurveyY[i][j][k][l] = source.fSurveyY[i][j][k][l];
    fSurveyZ[i][j][k][l] = source.fSurveyZ[i][j][k][l];
    fSurveyEX[i][j][k][l] = source.fSurveyEX[i][j][k][l];
    fSurveyEY[i][j][k][l] = source.fSurveyEY[i][j][k][l];
    fSurveyEZ[i][j][k][l] = source.fSurveyEZ[i][j][k][l];
  }
  for (int j=0; j<2; j++) for (int k=0; k<2; k++) for (int l=0; l<2; l++) {
    fSurveyX0[j][k][l] = source.fSurveyX0[j][k][l];
    fSurveyY0[j][k][l] = source.fSurveyY0[j][k][l];
    fSurveyZ0[j][k][l] = source.fSurveyZ0[j][k][l];
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
    for (int i=0; i<18; i++) for (int j=0; j<2; j++) for (int k=0; k<2; k++) for (int l=0; l<2; l++) {
      fSurveyX[i][j][k][l] = source.fSurveyX[i][j][k][l];
      fSurveyY[i][j][k][l] = source.fSurveyY[i][j][k][l];
      fSurveyZ[i][j][k][l] = source.fSurveyZ[i][j][k][l];
      fSurveyEX[i][j][k][l] = source.fSurveyEX[i][j][k][l];
      fSurveyEY[i][j][k][l] = source.fSurveyEY[i][j][k][l];
      fSurveyEZ[i][j][k][l] = source.fSurveyEZ[i][j][k][l];
    }
    for (int j=0; j<2; j++) for (int k=0; k<2; k++) for (int l=0; l<2; l++) {
      fSurveyX0[j][k][l] = source.fSurveyX0[j][k][l];
      fSurveyY0[j][k][l] = source.fSurveyY0[j][k][l];
      fSurveyZ0[j][k][l] = source.fSurveyZ0[j][k][l];
    }
    fComment = source.fComment;
  }

  return *this;

}

//_____________________________________________________________________________
AliTRDalignment& AliTRDalignment::operator*=(double fac) 
{
  //
  // multiplication operator
  //

  for (int i = 0; i <  18; i++) for (int j = 0; j < 6; j++) this->fSm[i][j] *= fac;
  for (int i = 0; i < 540; i++) for (int j = 0; j < 6; j++) this->fCh[i][j] *= fac;

  return *this;

}

//_____________________________________________________________________________
AliTRDalignment& AliTRDalignment::operator+=(const AliTRDalignment &source) 
{
  //
  // addition operator
  //

  for (int i = 0; i <  18; i++) for (int j = 0; j < 6; j++) this->fSm[i][j] += source.fSm[i][j];
  for (int i = 0; i < 540; i++) for (int j = 0; j < 6; j++) this->fCh[i][j] += source.fCh[i][j];

  return *this;

}

//_____________________________________________________________________________
AliTRDalignment& AliTRDalignment::operator-=(const AliTRDalignment &source) 
{
  //
  // subtraction operator
  //

  for (int i = 0; i <  18; i++) for (int j = 0; j < 6; j++) fSm[i][j] -= source.fSm[i][j];
  for (int i = 0; i < 540; i++) for (int j = 0; j < 6; j++) fCh[i][j] -= source.fCh[i][j];

  return *this;

}

//_____________________________________________________________________________
Bool_t AliTRDalignment::operator==(const AliTRDalignment &source) const
{
  //
  // comparison operator
  //

  Bool_t areEqual = 1;

  for (int i = 0; i <  18; i++) for (int j = 0; j < 6; j++) areEqual &= (fSm[i][j] == source.fSm[i][j]);
  for (int i = 0; i < 540; i++) for (int j = 0; j < 6; j++) areEqual &= (fCh[i][j] == source.fCh[i][j]);

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
void AliTRDalignment::SetSmRandom(double a[6]) 
{
  //
  // generate random gaussian supermodule data with sigmas a
  //

  double x[6];
  double xmax[6]={999, 0.6, 999, 999, 999, 999};

  for (int i = 0; i < 18; i++) {
    for (int j = 0; j < 6; j++) {
      do {x[j] = fRan.Gaus(0,a[j]);} while (TMath::Abs(x[j]) > xmax[j]);
    }
    SetSm(i,x);
    //PrintSm(i);
  }

}

//_____________________________________________________________________________
void AliTRDalignment::SetChRandom(double a[6]) 
{
  //
  // generate random gaussian chamber data with sigmas a
  //

  double x[6];

  for (int i = 0; i < 540; i++) {
    fRan.Rannor(x[0],x[1]);
    fRan.Rannor(x[2],x[3]);
    fRan.Rannor(x[4],x[5]);
    for (int j = 0; j < 6; j++) x[j] *= a[j];
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

  double a[6];

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

  double a[6];

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

  double a[6];

  a[0] = 0.002; // phi
  a[1] = 0.003; // z
  a[2] = 0.007; // r
  a[3] = 0.3/1000.0 / TMath::Pi()*180.0; // phi
  a[4] = 0.3/1000.0 / TMath::Pi()*180.0; // z
  a[5] = 0.1/1000.0 / TMath::Pi()*180.0; // r

  SetChRandom(a);

}

//_____________________________________________________________________________
void AliTRDalignment::PrintSm(int i, FILE *fp) const 
{
  //
  // print the supermodule data
  //

  fprintf(fp,"%4d   %11.4f %11.4f  %11.4f      %11.5f  %11.5f  %11.5f   %6d  %s\n"
	 ,i,fSm[i][0],fSm[i][1],fSm[i][2],fSm[i][3],fSm[i][4],fSm[i][5]
	 ,0,GetSmName(i));

}

//_____________________________________________________________________________
void AliTRDalignment::PrintCh(int i, FILE *fp) const 
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
    AliError(Form("cannot open input file %s",filename));
    return;
  }

  // supermodules

  for (int i = 0; i < 18; i++) {
    fi>>j>>x[0]>>x[1]>>x[2]>>x[3]>>x[4]>>x[5]>>volid>>syna;
    if (j != i) AliError(Form("sm %d expected, %d found",i,j));
    if (volid != 0) AliError(Form("sm %d volume id %d expected, %d found",i,0,volid));
    std::string symnam = GetSmName(i);
    if (syna != symnam) AliError(Form("sm %d name %s expected, %s found",i,symnam.data(),syna.data()));
    SetSm(i,x);
  }

  // chambers

  for (int i = 0; i < 540; i++) {
    fi>>j>>x[0]>>x[1]>>x[2]>>x[3]>>x[4]>>x[5]>>volid>>syna;
    if (j != i) AliError(Form("ch %d expected, %d found",i,j));
    if (volid != GetVoi(i)) AliError(Form("ch %d volume id %d expected, %d found",i,GetVoi(i),volid));
    std::string symnam = GetChName(i);
    if (syna != symnam) AliError(Form("ch %d name %s expected, %s found",i,symnam.data(),syna.data()));
    SetCh(i,x);
  }

  fi.close();

}

//_____________________________________________________________________________
void AliTRDalignment::ReadCurrentGeo() 
{
  //
  // use currently loaded geometry to determine misalignment by comparing 
  // original and misaligned matrix of the last node
  // Now, original, does not mean "ideal". It is the matrix before the alignment. 
  // So, if alignment was applied more than once, the numbers extracted will 
  // represent just the last alignment. -- check this!
  //

  TGeoPNEntry *pne;
  TGeoHMatrix *ideSm[18];  // ideal
  TGeoHMatrix *misSm[18];  // misaligned
  for (int i = 0; i < 18; i++) if ((pne = gGeoManager->GetAlignableEntry(GetSmName(i)))) {

    // read misaligned and original matrices

    TGeoPhysicalNode *node = pne->GetPhysicalNode();
    if (!node) AliError(Form("physical node entry %s has no physical node",GetSmName(i)));
    if (!node) continue;
    misSm[i] = new TGeoHMatrix(*node->GetNode(node->GetLevel())->GetMatrix());
    ideSm[i] = new TGeoHMatrix(*node->GetOriginalMatrix());

    // calculate the local misalignment matrices as inverse misaligned times ideal

    TGeoHMatrix mat(ideSm[i]->Inverse()); 
    mat.Multiply(misSm[i]);
    double *tra = mat.GetTranslation();
    double *rot = mat.GetRotationMatrix();
    double pars[6];
    pars[0] = tra[0];
    pars[1] = tra[1];
    pars[2] = tra[2];
    if (TMath::Abs(rot[0])<1e-7 || TMath::Abs(rot[8])<1e-7) AliError("Failed to extract roll-pitch-yall angles!");
    double raddeg = TMath::RadToDeg();
    pars[3] = raddeg * TMath::ATan2(-rot[5],rot[8]);
    pars[4] = raddeg * TMath::ASin(rot[2]);
    pars[5] = raddeg * TMath::ATan2(-rot[1],rot[0]);
    SetSm(i,pars);

    // cleanup

    delete ideSm[i];
    delete misSm[i];
  }

  TGeoHMatrix *ideCh[540]; // ideal
  TGeoHMatrix *misCh[540]; // misaligned
  for (int i = 0; i < 540; i++) if ((pne = gGeoManager->GetAlignableEntry(GetChName(i)))) {

    // read misaligned and original matrices

    TGeoPhysicalNode *node = pne->GetPhysicalNode();
    if (!node) AliError(Form("physical node entry %s has no physical node",GetChName(i)));
    if (!node) continue;
    misCh[i] = new TGeoHMatrix(*node->GetNode(node->GetLevel())->GetMatrix());
    ideCh[i] = new TGeoHMatrix(*node->GetOriginalMatrix());

    // calculate the local misalignment matrices as inverse misaligned times ideal

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

    // cleanup
    delete ideCh[i];
    delete misCh[i];
  }

  return;

}

//_____________________________________________________________________________
void AliTRDalignment::ReadRoot(char *filename) 
{
  //
  // read the alignment data from root file
  //

  TFile fi(filename,"READ");

  if (fi.IsOpen()) {
    TClonesArray *ar = (TClonesArray*) fi.Get("TRDAlignObjs");
    ArToNumbers(ar);
    fi.Close();
  } 
  else AliError(Form("cannot open input file %s",filename));

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
    fComment.SetString(e->GetMetaData()->GetComment());
    TClonesArray *ar = (TClonesArray *) e->GetObject();
    ArToNumbers(ar);
    fi.Close();
  } 
  else AliError(Form("cannot open input file %s",filename));

  return;

}

//_____________________________________________________________________________
void AliTRDalignment::ReadDB(char *db, char *path, int run
                           , int version, int subversion)
{
  //
  // read the alignment data from database
  //

  AliCDBManager *cdb     = AliCDBManager::Instance();
  AliCDBStorage *storLoc = cdb->GetStorage(db);
  AliCDBEntry   *e       = storLoc->Get(path,run,version,subversion);
  if (e) {
    e->PrintMetaData();
    fComment.SetString(e->GetMetaData()->GetComment());
    TClonesArray  *ar      =  (TClonesArray *) e->GetObject();
    ArToNumbers(ar);
  }
}

//_____________________________________________________________________________
Bool_t AliTRDalignment::DecodeSurveyPointName(TString pna, Int_t &sm, Int_t &iz, 
					      Int_t &ir, Int_t &iphi) {
  // decode the survey point name and extract the sm, z, r and phi indices
  
  if (pna(0,6)!="TRD_sm") {
    AliError(Form("unexpected point name: %s",pna.Data()));
    return kFALSE;
  }
  sm = atoi(pna(6,2).Data()); // supermodule number
  iz = -1;
  if (pna(8) == 'a') iz=0; // anticlockwise, positive z
  if (pna(8) == 'c') iz=1; // clockwise, negative z
  ir = -1;
  if (pna(9) == 'l') ir=0; // low radius
  if (pna(9) == 'h') ir=1; // high radius
  iphi = -1;
  if (pna(10) == '0') iphi = 0; // low phi within supermodule
  if (pna(10) == '1') iphi = 1; // high phi within supermodule
  if (sm>=0 && sm<18 && iz>=0 && iz<2 && ir>=0 && ir<2 && iphi>=0 && iphi<2) return kTRUE;
  AliError(Form("cannot decode point name: %s",pna.Data()));
  return kFALSE;
}

//_____________________________________________________________________________
void AliTRDalignment::ReadSurveyReport(char *filename) 
{
  //
  // Read survey report and store the numbers in fSurveyX, fSurveyY, fSurveyZ, 
  // and fSurveyE.  Store the survey info in the fComment.
  // Each supermodule has 8 survey points. The point names look like 
  // TRD_sm08ah0 and have the following meaning. 
  //
  // sm00..17 mean supermodule 0 through 17, following the phi.
  // Supermodule 00 is between phi=0 and phi=20 degrees.
  //
  // a or c denotes the anticlockwise and clockwise end of the supermodule
  // in z. Clockwise end is where z is negative and where the muon arm sits.
  //
  // l or h denote low radius and high radius holes
  //
  // 0 or 1 denote the hole at smaller and at larger phi, respectively.
  //

  // read the survey file

  fstream in(filename,fstream::in);
  if (!in) {
    AliError(Form("cannot open input file %s",filename));
    return;
  }

  // loop through the lines of the file until the beginning of data

  TString title,date,subdetector,url,version,observations,system,units;
  while (1) {
    char pee=in.peek();
    if (pee==EOF) break; 
    TString line;
    line.ReadLine(in);
    if (line.Contains("Title:"))        title.ReadLine(in);
    if (line.Contains("Date:"))         date.ReadLine(in);
    if (line.Contains("Subdetector:"))  subdetector.ReadLine(in);
    if (line.Contains("URL:"))          url.ReadLine(in);
    if (line.Contains("Version:"))      version.ReadLine(in);
    if (line.Contains("Observations:")) observations.ReadLine(in);
    if (line.Contains("System:"))       system.ReadLine(in);
    if (line.Contains("Units:"))        units.ReadLine(in);
    if (line.Contains("Data:"))         break;
  }

  // check what we found so far (watch out, they have \r at the end)

  std::cout<<"title .........."<<title<<std::endl;
  std::cout<<"date ..........."<<date<<std::endl;
  std::cout<<"subdetector ...."<<subdetector<<std::endl;
  std::cout<<"url ............"<<url<<std::endl;
  std::cout<<"version ........"<<version<<std::endl;
  std::cout<<"observations ..."<<observations<<std::endl;
  std::cout<<"system ........."<<system<<std::endl;
  std::cout<<"units .........."<<units<<std::endl;

  if (!subdetector.Contains("TRD")) {
    AliWarning(Form("Not a TRD survey file, subdetector = %s",subdetector.Data()));
    return;
  }
  double tocm = 0; // we want to have it in cm
  if (units.Contains("mm"))      tocm = 0.1;
  else if (units.Contains("cm")) tocm = 1.0;
  else if (units.Contains("m"))  tocm = 100.0;
  else if (units.Contains("pc")) tocm = 3.24078e-15;
  else {
    AliError(Form("unexpected units: %s",units.Data()));
    return;
  }
  if (!system.Contains("ALICEPH")) {
    AliError(Form("wrong system: %s, should be ALICEPH",system.Data()));
    return;
  }

  // scan the rest of the file which should contain list of surveyed points
  // for every point, decode the point name and store the numbers in the right 
  // place in the arrays fSurveyX etc.

  while (1) {
    TString pna; // point name
    char type, target;
    double x,y,z,precision;
    
    in >> pna >> x >> y >> z >> type >> target >> precision;  
    if (in.fail()) break;
    Int_t i,j,k,l;
    if (DecodeSurveyPointName(pna,i,j,k,l)) {
      fSurveyX[i][j][k][l] = tocm*x;
      fSurveyY[i][j][k][l] = tocm*y;
      fSurveyZ[i][j][k][l] = tocm*z;
      fSurveyEX[i][j][k][l] = precision/10; // "precision" is supposed to be in mm
      fSurveyEY[i][j][k][l] = precision/10; // "precision" is supposed to be in mm
      fSurveyEZ[i][j][k][l] = precision/10; // "precision" is supposed to be in mm
      // if, at some point, separate precision numbers for x,y,z show up in the 
      // survey reports the function will fail here
      printf("decoded %s %02d %d %d %d  %8.2f %8.2f %8.2f %6.2f %6.2f %6.2f\n", 
	     pna.Data(), i, j, k, l,
	     fSurveyX[i][j][k][l], fSurveyY[i][j][k][l], fSurveyZ[i][j][k][l],
	     fSurveyEX[i][j][k][l], fSurveyEY[i][j][k][l], fSurveyEZ[i][j][k][l]);
    } else AliError(Form("cannot decode point name: %s",pna.Data()));
  }
  in.close();
  TString info = "Survey "+title+" "+date+" "+url+" "+version+" "+observations;
  info.ReplaceAll("\r","");
  fComment.SetString(info.Data());
 
}

//_____________________________________________________________________________
void AliTRDalignment::ReadSurveyReport(AliSurveyObj *so) 
{
  //
  // Read survey report and store the numbers in fSurveyX, fSurveyY, fSurveyZ, 
  // and fSurveyE.  Store the survey info in the fComment.
  // Each supermodule has 8 survey points. The point names look like 
  // TRD_sm08ah0 and have the following meaning. 
  //
  // sm00..17 mean supermodule 0 through 17, following the phi.
  // Supermodule 00 is between phi=0 and phi=20 degrees.
  //
  // a or c denotes the anticlockwise and clockwise end of the supermodule
  // in z. Clockwise end is where z is negative and where the muon arm sits.
  //
  // l or h denote low radius and high radius holes
  //
  // 0 or 1 denote the hole at smaller and at larger phi, respectively.
  //

  // read and process the data from the survey object

  Int_t size = so->GetEntries();
  printf("-> %d\n", size);

  TString title        = so->GetReportTitle();
  TString date         = so->GetReportDate();
  TString subdetector  = so->GetDetector();
  TString url          = so->GetURL();
  TString report       = so->GetReportNumber();
  TString version      = so->GetReportVersion();
  TString observations = so->GetObservations();
  TString system       = so->GetCoordSys();
  TString units        = so->GetUnits();

  // check what we found so far (watch out, they have \r at the end)

  std::cout<<"title .........."<<title<<std::endl;
  std::cout<<"date ..........."<<date<<std::endl;
  std::cout<<"subdetector ...."<<subdetector<<std::endl;
  std::cout<<"url ............"<<url<<std::endl;
  std::cout<<"version ........"<<version<<std::endl;
  std::cout<<"observations ..."<<observations<<std::endl;
  std::cout<<"system ........."<<system<<std::endl;
  std::cout<<"units .........."<<units<<std::endl;

  if (!subdetector.Contains("TRD")) {
    AliWarning(Form("Not a TRD survey file, subdetector = %s",subdetector.Data()));
    return;
  }
  double tocm = 0; // we want to have it in cm
  if (units.Contains("mm"))      tocm = 0.1;
  else if (units.Contains("cm")) tocm = 1.0;
  else if (units.Contains("m"))  tocm = 100.0;
  else if (units.Contains("pc")) tocm = 3.24078e-15;
  else {
    AliError(Form("unexpected units: %s",units.Data()));
    return;
  }
  if (!system.Contains("ALICEPH")) {
    AliError(Form("wrong system: %s, should be ALICEPH",system.Data()));
    return;
  }

  // for every survey point, decode the point name and store the numbers in 
  // the right place in the arrays fSurveyX etc.

  TObjArray *points = so->GetData();
  for (int ip = 0; ip<points->GetEntries(); ++ip) {
    AliSurveyPoint *po = (AliSurveyPoint *) points->At(ip);
    TString pna = po->GetPointName();
    Int_t i,j,k,l;
    if (DecodeSurveyPointName(pna,i,j,k,l)) {
      fSurveyX[i][j][k][l] = tocm*po->GetX();
      fSurveyY[i][j][k][l] = tocm*po->GetY();
      fSurveyZ[i][j][k][l] = tocm*po->GetZ();
      fSurveyEX[i][j][k][l] = po->GetPrecisionX()/10; // "precision" is supposed to be in mm
      fSurveyEY[i][j][k][l] = po->GetPrecisionY()/10;
      fSurveyEZ[i][j][k][l] = po->GetPrecisionZ()/10;
      printf("decoded %s %02d %d %d %d  %8.2f %8.2f %8.2f %6.2f %6.2f %6.2f\n", 
	     pna.Data(), i, j, k, l,
	     fSurveyX[i][j][k][l], fSurveyY[i][j][k][l], fSurveyZ[i][j][k][l],
	     fSurveyEX[i][j][k][l], fSurveyEY[i][j][k][l], fSurveyEZ[i][j][k][l]);
    } else AliError(Form("cannot decode point name: %s",pna.Data()));
  }

  TString info = "Survey "+title+" "+date+" "+url+" "+report+" "+version+" "+observations;
  info.ReplaceAll("\r","");
  fComment.SetString(info.Data());			 
}

//_____________________________________________________________________________
double AliTRDalignment::SurveyChi2(int i, double *a) {

  //
  // Compare the survey results to the ideal positions of the survey marks
  // in the local frame of supermodule. When transforming, use the alignment 
  // parameters a[6]. Return chi-squared.
  //

  if (!IsGeoLoaded()) return 0;
  printf("Survey of supermodule %d\n",i);
  AliAlignObjParams al(GetSmName(i),0,a[0],a[1],a[2],a[3],a[4],a[5],0);
  TGeoPNEntry      *pne  = gGeoManager->GetAlignableEntry(GetSmName(i));
  if (!pne) AliError(Form("no such physical node entry: %s",GetSmName(i)));
  TGeoPhysicalNode *node = pne->GetPhysicalNode();
  if (!node) AliError(Form("physical node entry %s has no physical node",GetSmName(i))); 

  //  al.ApplyToGeometry();    
  //  node = pne->GetPhysicalNode(); // changed in the meantime
  //  TGeoHMatrix *ma = node->GetMatrix();

  // a less destructive method (it does not modify geometry), gives the same result:

  TGeoHMatrix *ma = new TGeoHMatrix();
  al.GetLocalMatrix(*ma);
  ma->MultiplyLeft(node->GetMatrix()); // global trafo, modified by a[]

  double chi2=0;
  printf("              sm   z   r  phi    x (lab phi)  y (lab z)   z (lab r)   all in cm\n");
  for (int j=0; j<2; j++) for (int k=0; k<2; k++) for (int l=0; l<2; l++) {
    if (fSurveyEX[i][j][k][l] == 0.0 
	&& fSurveyEY[i][j][k][l] == 0.0 
	&& fSurveyEZ[i][j][k][l] == 0.0) continue; // no data for this survey point
    double master[3] = {fSurveyX[i][j][k][l],fSurveyY[i][j][k][l],fSurveyZ[i][j][k][l]};
    double local[3];
    ma->MasterToLocal(master,local);
    double dx = local[0]-fSurveyX0[j][k][l];
    double dy = local[1]-fSurveyY0[j][k][l];
    double dz = local[2]-fSurveyZ0[j][k][l];
    chi2 += dx*dx/fSurveyEX[i][j][k][l]/fSurveyEX[i][j][k][l];
    chi2 += dy*dy/fSurveyEY[i][j][k][l]/fSurveyEY[i][j][k][l];
    chi2 += dz*dz/fSurveyEZ[i][j][k][l]/fSurveyEZ[i][j][k][l];
    printf("local survey %3d %3d %3d %3d %12.3f %12.3f %12.3f\n",i,j,k,l,local[0],local[1],local[2]);
    printf("local ideal                  %12.3f %12.3f %12.3f\n",fSurveyX0[j][k][l],
	   fSurveyY0[j][k][l],fSurveyZ0[j][k][l]);
    printf("difference                   %12.3f %12.3f %12.3f\n",dx,dy,dz);
  }
  printf("chi2 = %.2f\n",chi2);
  return chi2;
}

//_____________________________________________________________________________
void trdAlignmentFcn(int &npar, double *g, double &f, double *par, int iflag) {

  // 
  // Standard function as needed by Minuit-like minimization procedures. 
  // For the set of parameters par calculates and returns chi-squared.
  //

  // smuggle a C++ object into a C function
  AliTRDalignment *alignment = (AliTRDalignment*) gMinuit->GetObjectFit(); 

  f = alignment->SurveyChi2(par);
  if (iflag==3) {}
  if (npar) {} 
  if (g) {} // no warnings about unused stuff...

}

//_____________________________________________________________________________
void AliTRDalignment::SurveyToAlignment(int i,char *flag) {

  //
  // Find the supermodule alignment parameters needed to make the survey 
  // results coincide with the ideal positions of the survey marks.
  // The string flag should look like "101000"; the six characters corresponds 
  // to the six alignment parameters and 0/1 mean that the parameter should 
  // be fixed/released in the fit. 

  if (strlen(flag)!=6) {
    AliError(Form("unexpected flag: %s",flag));
    return;
  }

  printf("Finding alignment matrix for supermodule %d\n",i);
  fIbuffer[0] = i; // store the sm number in the buffer so minuit can see it

  TFitter fitter(100);
  gMinuit->SetObjectFit(this);
  fitter.SetFCN(trdAlignmentFcn);
  fitter.SetParameter(0,"dx",0,0.5,0,0);
  fitter.SetParameter(1,"dy",0,0.5,0,0);
  fitter.SetParameter(2,"dz",0,0.5,0,0);
  fitter.SetParameter(3,"rx",0,0.1,0,0);
  fitter.SetParameter(4,"ry",0,0.1,0,0);
  fitter.SetParameter(5,"rz",0,0.1,0,0);

  for (int j=0; j<6; j++) if (flag[j]=='0') fitter.FixParameter(j);

  double arglist[100];
  arglist[0] = 2;
  fitter.ExecuteCommand("SET PRINT", arglist, 1);
  fitter.ExecuteCommand("SET ERR", arglist, 1);
  arglist[0]=50;
  //fitter.ExecuteCommand("SIMPLEX", arglist, 1);
  fitter.ExecuteCommand("MINIMIZE", arglist, 1);
  fitter.ExecuteCommand("CALL 3", arglist,0);
  double a[6];
  for (int j=0; j<6; j++) a[j] = fitter.GetParameter(j);
  SetSm(i,a);
  for (int j=0; j<6; j++) printf("%10.3f ",fitter.GetParameter(j));   
  printf("\n");
  for (int j=0; j<6; j++) printf("%10.3f ",fitter.GetParError(j));
  printf("\n");

}

//_____________________________________________________________________________
void AliTRDalignment::ReadAny(char *filename) 
{
  //
  // read the alignment data from any kind of file
  //

  TString fist(filename);
  if (fist.EndsWith(".txt")) ReadAscii(filename);
  if (fist.EndsWith(".dat")) ReadAscii(filename);
  if (fist.EndsWith(".root")) {
    if (fist.Contains("Run")) ReadDB(filename);
    else ReadRoot(filename);
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

  TClonesArray *ar = new TClonesArray("AliAlignObjParams",10000);
  NumbersToAr(ar);
  TFile fo(filename,"RECREATE");
  if (fo.IsOpen()) {
    fo.cd();
    fo.WriteObject(ar,"TRDAlignObjs","kSingleKey");
    fo.Close();
  } 
  else AliError(Form("cannot open output file %s",filename));

  delete ar;

}

//_____________________________________________________________________________
void AliTRDalignment::WriteDB(char *filename, int run0, int run1) 
{
  //
  // dumping on a DB-like file
  //

  TClonesArray   *ar = new TClonesArray("AliAlignObjParams",10000);
  NumbersToAr(ar);
  const Char_t *path = "TRD/Align/Data";
  AliCDBId id(path,run0,run1);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Dariusz Miskowiec");
  md->SetComment(fComment.GetString().Data());
  AliCDBEntry    *e  = new AliCDBEntry(ar, id, md);
  TFile fi(filename,"RECREATE");
  if (fi.IsOpen()) {
    e->Write();
    fi.Close();
  } 
  else AliError(Form("cannot open input file %s",filename));

  delete e;
  delete md;
  delete ar;

  return;

}

//_____________________________________________________________________________
void AliTRDalignment::WriteDB(char *db, char *path, int run0, int run1) 
{
  //
  // store the alignment data in database
  //

  TClonesArray   *ar      = new TClonesArray("AliAlignObjParams",10000);
  NumbersToAr(ar);
  AliCDBManager  *cdb     = AliCDBManager::Instance();
  AliCDBStorage  *storLoc = cdb->GetStorage(db);
  AliCDBMetaData *md      = new AliCDBMetaData();
  md->SetResponsible("Dariusz Miskowiec");
  md->SetComment(fComment.GetString().Data());
  AliCDBId id(path,run0,run1);
  storLoc->Put(ar,id,md);
  md->Delete();
  delete ar;

}

//_____________________________________________________________________________
void AliTRDalignment::WriteGeo(char *filename) 
{
  //
  // apply misalignment to current geometry and store the 
  // resulting geometry on a root file
  //

  TClonesArray *ar = new TClonesArray("AliAlignObjParams",10000);
  NumbersToAr(ar);
  delete ar;
  gGeoManager->Export(filename);

}

//_____________________________________________________________________________
double AliTRDalignment::GetSmRMS(int xyz) const 
{
  //
  // rms fSm[][xyz]
  //

  double s1 = 0.0;
  double s2 = 0.0;
  for (int i = 0; i < 18; i++) {
    s1 += fSm[i][xyz];
    s2 += fSm[i][xyz]*fSm[i][xyz];
  }
  double rms2 = s2/18.0 - s1*s1/18.0/18.0;

  return rms2>0 ? sqrt(rms2) : 0.0;

}

//_____________________________________________________________________________
double AliTRDalignment::GetChRMS(int xyz) const
{
  //
  // rms fCh[][xyz]
  //

  double s1 =0.0;
  double s2 =0.0;
  for (int i = 0; i < 540; i++) {
    s1 += fCh[i][xyz];
    s2 += fCh[i][xyz]*fCh[i][xyz];
  }
  double rms2 = s2/540.0 - s1*s1/540.0/540.0;

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
  // for each of the alignment objects in array ar extract the six local 
  // alignment parameters; recognize by name to which supermodule or chamber 
  // the alignment object pertains; set the respective fSm or fCh
  //

  ar->Sort();
  if (!IsGeoLoaded()) return;
  for (int i = 0; i < ar->GetEntries(); i++) {
    AliAlignObj *aao = (AliAlignObj *) ar->At(i);
    aao->ApplyToGeometry();
  }
  SetZero();
  ReadCurrentGeo();

}

//_____________________________________________________________________________
void AliTRDalignment::NumbersToAr(TClonesArray *ar) 
{
  //
  // build array of AliAlignObj objects based on fSm and fCh data
  // at the same time, apply misalignment to the currently loaded geometry
  // it is important to apply misalignment of supermodules before creating 
  // alignment objects for chambers
  //

  if (!IsGeoLoaded()) return;
  TClonesArray &alobj = *ar;
  int nobj = 0;
  for (int i = 0; i <  18; i++) {      
      new(alobj[nobj]) AliAlignObjParams(GetSmName(i)
                                        ,0 
					,fSm[i][0],fSm[i][1],fSm[i][2]
					,fSm[i][3],fSm[i][4],fSm[i][5]
					,0);
    ((AliAlignObj *) alobj[nobj])->ApplyToGeometry();
    nobj++;
  }

  for (int i = 0; i < 540; i++) {
    if (gGeoManager->GetAlignableEntry(GetChName(i))) {
      new(alobj[nobj]) AliAlignObjParams(GetChName(i)
					 ,GetVoi(i)
					 ,fCh[i][0],fCh[i][1],fCh[i][2]
					 ,fCh[i][3],fCh[i][4],fCh[i][5]
					 ,0);
      ((AliAlignObj *) alobj[nobj])->ApplyToGeometry();
      nobj++;
    }
  }
  AliInfo("current geometry modified");

}

//_____________________________________________________________________________
int AliTRDalignment::IsGeoLoaded() 
{
  //
  // check whether a geometry is loaded
  // issue a warning if geometry is not ideal
  //

  if (gGeoManager) {
    if (gGeoManager->GetListOfPhysicalNodes()->GetEntries()) AliWarning("current geometry is not ideal");
    return 1;
  } else {
    AliError("first load geometry by calling TGeoManager::Import(filename)");
    return 0;
  }

}

//_____________________________________________________________________________
