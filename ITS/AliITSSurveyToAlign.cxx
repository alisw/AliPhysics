/**************************************************************************
 * Copyright(c) 2008-2010, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////
//   Class to convert survey tables in alignment objects
//   for SSD and SDD
//   origin: Marco Van Leeuwen (m.vanleeuwen1@uu.nl)
//           Panos.Christakoglou (Panos.Christakoglou@cern.ch)
//           Martin Poghosyan (Martin.Poghosyan@to.infn.it)
//////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TGeoPhysicalNode.h"
#include "TMatrixD.h"
#include "TMath.h"

#include "AliITSSurveyToAlign.h"
#include "AliSurveyPoint.h"
#include "AliAlignObjParams.h"
#include "AliGeomManager.h"

#include "AliLog.h"

#include "AliCDBManager.h"

#include "AliITSgeomTGeo.h"

ClassImp(AliITSSurveyToAlign)

const Double_t AliITSSurveyToAlign::fgkLocR[6][3]={{ 3.24,0.21905,-2.4},
                                                   { 3.58,0.21905, 0. },
						   { 3.24,0.21905,+2.4},
						   {-3.24,0.21905,+2.4},
						   {-3.58,0.21905, 0. },
						   {-3.24,0.21905,-2.4}};

const Double_t AliITSSurveyToAlign::fgkLocL[6][3]={{-3.24,0.21905, 2.4},
						   {-3.58,0.21905, 0. },
						   {-3.24,0.21905,-2.4},
						   { 3.24,0.21905,-2.4},
						   { 3.58,0.21905, 0. },
						   { 3.24,0.21905, 2.4}};

const Double_t kRadToDeg = 180./TMath::Pi();

//________________________________________________________________________
AliITSSurveyToAlign::AliITSSurveyToAlign(Int_t run, Int_t repSDD, Int_t repVerSDD, Int_t repModSSD, Int_t repModVerSSD, Int_t repLaddSSD, Int_t repLaddVerSSD) :
  AliSurveyToAlignObjs(),
  fRun(run),
  fSDDrepNumber(repSDD),
  fSDDrepVersion(repVerSDD),
  fSSDModuleRepNumber(repModSSD),
  fSSDModuleRepVersion(repModVerSSD),
  fSSDLadderRepNumber(repLaddSSD),
  fSSDLadderRepVersion(repLaddVerSSD)
 {
  //
  //  default constructor
  //  Arguments are report numbers for survey data. 
  //  The defaults point to reports from detector construction
  // 
}

//_________________________________________________________________________
AliITSSurveyToAlign::AliITSSurveyToAlign(const AliITSSurveyToAlign &align) :
  AliSurveyToAlignObjs(align),
  fRun(align.fRun),
  fSDDrepNumber(align.fSDDrepNumber),
  fSDDrepVersion(align.fSDDrepVersion),
  fSSDModuleRepNumber(align.fSSDModuleRepNumber),
  fSSDModuleRepVersion(align.fSSDModuleRepVersion),
  fSSDLadderRepNumber(align.fSSDLadderRepNumber),
  fSSDLadderRepVersion(align.fSSDLadderRepVersion)
{
  //
  //  copy constructor 
  //
}

//__________________________________________________________________________
AliITSSurveyToAlign & AliITSSurveyToAlign::operator =(const AliITSSurveyToAlign& /* align */) {
  //
  // assignment operator - dummy
  //

  return (*this);
}

//__________________________________________________________________________
AliITSSurveyToAlign::~AliITSSurveyToAlign() {
  //
  // destructor
  //
}

//______________________________________________________________________
void AliITSSurveyToAlign::Run() { 
  //
  // Runs the full chain
  // User should call StoreAlignObjToFile or StoreAlignObjToCDB afterwards to 
  // store output (not included here to leave the choice between the two)
  //

  // Load ideal geometry from the OCDB
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(fRun);
  AliGeomManager::LoadGeometry();

  if(!CreateAlignObjs()) AliError("Construction of alignment objects from survey failed!");
} 

//______________________________________________________________________
Bool_t AliITSSurveyToAlign::CreateAlignObjs() { 
  // Fill the array of alignment objects with alignment objects
  // from survey for all three subdetectors
  //

  //for SPD
  CreateAlignObjDummySPD();

  // for SDD
  if(!LoadSurveyFromAlienFile("ITS", fSDDrepNumber, fSDDrepVersion)){
      AliError("Loading of alignment objects from survey for SDD failed!");
      return kFALSE;
  }
  CreateAlignObjSDD();

  // for SSD ladders
  if(!LoadSurveyFromAlienFile("ITS", fSSDLadderRepNumber, fSSDLadderRepVersion)){
      AliError("Loading of alignment objects from survey for SSD ladders failed!");
      return kFALSE;
  }
  CreateAlignObjSSDLadders();

  // for SSD modules
  if(!ApplyAlignObjSSDLadders()) return kFALSE; // needed to build correctly the objects for SSD modules
  if(!LoadSurveyFromAlienFile("ITS", fSSDModuleRepNumber, fSSDModuleRepVersion)){
      AliError("Loading of alignment objects from survey for SSD modules failed!");
      return kFALSE;
  }
  CreateAlignObjSSDModules();

  return kTRUE;
}

//______________________________________________________________________
void AliITSSurveyToAlign::CreateAlignObjDummySPD(){
  // 
  // Create alignObjs for SPD
  //    For the moment, uses 0,0,0,0,0,0
  //
  for(Int_t imod = 0; imod < 240; imod++) {
    Int_t ilayer = (imod < 80) ? AliGeomManager::kSPD1 : AliGeomManager::kSPD2;
    Int_t imodule = (imod < 80) ? imod : imod - 80;

    Int_t uid = AliGeomManager::LayerToVolUID(ilayer,imodule);
    const Char_t *symname = AliGeomManager::SymName(uid);

    new((*fAlignObjArray)[imod]) AliAlignObjParams(symname, uid, 0., 0., 0., 0., 0., 0., kTRUE);
  }//module loop

}

//______________________________________________________________________
void AliITSSurveyToAlign::CreateAlignObjSDD(){
  //
  // Create alignment objects for SDD
  // Called by Run()
  //
  Int_t uid = 0;
  const char* symname = 0;
  AliSurveyPoint* pt = 0;
 
  Int_t iModuleIndex=240;
  Int_t iModule0=0;
  Int_t iLadder0=0;
  Int_t iLayer0=3;
  Int_t nModules=0;

  if (fSurveyPoints == 0 || fSurveyPoints->GetEntries() == 0) {
    AliWarning("SDD survey data are not available, using zero values");
    CreateAlignObjDummySDD();
    return;
  }

  for(Int_t imod = 1; imod < fSurveyPoints->GetEntries(); imod++) {
    pt = (AliSurveyPoint*) fSurveyPoints->At(imod);
    if(!pt) continue;

    Int_t iLayer, iLadder, iModule, iPoint;
    ReadPointNameSDD(pt->GetName(),iLayer, iLadder, iModule, iPoint);

    if(iModule==iModule0)
    {
      fSDDmeP[iPoint][0]=pt->GetX();
      fSDDmeP[iPoint][1]=pt->GetY();
      fSDDmeP[iPoint][2]=pt->GetZ();
      fSDDmeP[iPoint][3]=pt->GetPrecisionX();
      fSDDmeP[iPoint][4]=pt->GetPrecisionY();
      fSDDmeP[iPoint][5]=pt->GetPrecisionZ();
      fSDDisMe[iPoint]=kTRUE;

      if(iLayer==3) uid = AliGeomManager::LayerToVolUID(iLayer0,iModuleIndex-240);
      if(iLayer==4) uid = AliGeomManager::LayerToVolUID(iLayer0,iModuleIndex-324);
      symname = AliGeomManager::SymName(uid);
      GetIdPosSDD(uid,iLayer0, iModule0, iPoint);
      nModules++;
    }
    //    cout << "Points red module " << imod << endl;
    if((iModule!=iModule0)||(imod==(fSurveyPoints->GetEntries()-1)))
    {
      ConvertToRSofModulesAndRotSDD(iLayer0, iModule0);

      Double_t tet = 0.;
      Double_t psi =0.;
      Double_t phi = 0.;
      Double_t x0  = 0.;
      Double_t y0  =0.;
      Double_t z0  = 0.;

      if(nModules==2) CalcShiftSDD(x0,y0,z0);
      if(nModules>2)   CalcShiftRotSDD(tet, psi, phi, x0, y0, z0);
      tet*=kRadToDeg;
      psi*=kRadToDeg;
      phi*=kRadToDeg;

//    printf("%s  %d  %f  %f  %f  %f  %f  %f\n",symname, uid, x0/10., y0/10., z0/10., psi, tet, phi);
//      cout << "Allocate alignobjparams " << imod << endl;
      new((*fAlignObjArray)[iModuleIndex]) AliAlignObjParams(symname, uid, x0/10., y0/10., z0/10., psi, tet, phi, kFALSE);

      iModule0=iModule;
      iLayer0=iLayer;
      iLadder0=iLadder;
      nModules=0;
      iModuleIndex = AliITSgeomTGeo::GetModuleIndex(iLayer,iLadder+1,iModule+1);
      for(Int_t i=0; i<6;i++) fSDDisMe[i]=kFALSE;
      if(imod!=(fSurveyPoints->GetEntries()-1)) imod--;
    }
  }//module loop
}

//______________________________________________________________________
void AliITSSurveyToAlign::CreateAlignObjDummySDD(){
  // 
  // Create empty alignment objects
  // Used when fSurveySDD == 0
  //
  for(Int_t imod = 0; imod < 260; imod++) {

    Int_t ilayer = (imod < 84) ? AliGeomManager::kSDD1 : AliGeomManager::kSDD2;
    Int_t imodule = (imod < 84) ? imod : imod - 84;

    Int_t uid = AliGeomManager::LayerToVolUID(ilayer,imodule);
    const Char_t *symname = AliGeomManager::SymName(uid);

    new((*fAlignObjArray)[imod+240]) AliAlignObjParams(symname, uid, 0., 0., 0., 0., 0., 0., kTRUE);
  }//module loop
}

//______________________________________________________________________
void AliITSSurveyToAlign::CreateAlignObjSSDModules(){
  //
  // Create alignment objects for SSD modules
  // Objects for SSD ladders must be applied to geometry first
  //
  Double_t sx, sz;
  const Float_t kMu2Cm = 1e-4;
  const Float_t kSensLength = 7.464;
  const Int_t kSSDMODULES = 1698;

  if (fSurveyPoints == 0 || fSurveyPoints->GetEntries() == 0) {
    AliWarning("SSD module survey data not available; using dummy values");
    CreateAlignObjDummySSDModules();
    return;
  }

  // First do module-by-module

  for(Int_t imod = 500; imod < kSSDMODULES + 500; imod++) {
    Int_t iLayer, iLadder, iLaddMod;
    AliITSgeomTGeo::GetModuleId(imod,iLayer,iLadder,iLaddMod);  // returns 1-based numbers
 
    TString pname="ITS/SSD";
    pname += iLayer-1;
    pname += "/Ladder";
    pname += iLadder-1;
    pname += "/Sensor";
    pname += iLaddMod-1;
    AliSurveyPoint *pt1 = (AliSurveyPoint*) fSurveyPoints->FindObject(pname+"/Point0");
    AliSurveyPoint *pt2 = (AliSurveyPoint*) fSurveyPoints->FindObject(pname+"/Point1");
    if(!pt1 || !pt2) {
      AliWarning(Form("No Survey points for iladd %d imod %d",iLadder,iLaddMod));
      continue;
    }

    sx = 0.5*(pt1->GetX() + pt2->GetX()) * kMu2Cm;
    sz = 0.5*(pt1->GetZ() + pt2->GetZ()) * kMu2Cm;

    // Minus sign to change local coordinate convention 
    Float_t theta = -(pt2->GetZ() - pt1->GetZ())*kMu2Cm/kSensLength;

    theta *= kRadToDeg;
    Int_t iLayMod = imod - 500;
    if (iLayer == 6)
      iLayMod -= 748;
    Int_t uid = AliGeomManager::LayerToVolUID(iLayer,iLayMod);

    const Char_t *symname = AliGeomManager::SymName(uid);
    if (pname.CompareTo(symname) != 0)
      AliWarning(Form("Mapping mismatch survey point %s volume name %s",pname.Data(),symname));
    /*
    if (imod >= 676 && imod <= 697) {
      cout << "ilayer " << iLayer << " imod " << imod 
	   << " uid " << uid << " name " << symname 
	   << " survey shift " << sx << " " << 0 << " " << sz << endl
	   << " theta " << theta << endl;
    }
    */
    new((*fAlignObjArray)[imod]) AliAlignObjParams(symname, uid, sx, 0, sz, 0., theta, 0., kFALSE);
  } //module loop
}

//______________________________________________________________________
Bool_t AliITSSurveyToAlign::ApplyAlignObjSSDLadders(){
  //
  //   Apply alignment objects for SSD ladders to geometry, needed to correctly
  //   build alignment objects for SSD modules
  // 
  TClonesArray* tobeApplied = new TClonesArray("AliAlignObjParams",72);
  Int_t ii=0;
  for(Int_t jj=0; jj<fAlignObjArray->GetEntriesFast(); jj++)
  {
      AliAlignObjParams* ap = dynamic_cast<AliAlignObjParams*> (fAlignObjArray->UncheckedAt(jj));
      if(ap) 
      {
	  TString sName(ap->GetSymName());
	  if(sName.Contains("SSD") && sName.Contains("Ladder"))
	      (*tobeApplied)[ii++] = (AliAlignObjParams*) fAlignObjArray->UncheckedAt(jj);
      }
  }
  AliInfo(Form(" %d alignment objects for SSD ladders applied to geometry.",tobeApplied->GetEntriesFast()));

  return(AliGeomManager::ApplyAlignObjsToGeom(*tobeApplied));
}

//______________________________________________________________________
void AliITSSurveyToAlign::CreateAlignObjSSDLadders(){
  //
  //   Alignment objects from survey for SSD ladders (Torino data)
  // 
  const Float_t kLaddLen5 = 90.27;  // Layer 5: distance between mouting points
  const Float_t kLaddLen6 = 102.0;  // Layer 6: distance between mouting points
  const Float_t zLag = 2.927;         // Distance between V mounting point and Zloc = 0
                                    // = half ladder length - nom z-position of ladder from gGeoManager
  const Float_t kMu2Cm = 1e-4;

  TString ssdName = "ITS/SSD";

  TObjArray *ladderPoints = fSurveyPoints;  
  if (ladderPoints == 0 || ladderPoints->GetEntries() == 0) {
    AliWarning("No SSD Ladder alignment points found. Skipping");
    return;
  }
  if (ladderPoints->GetEntries()!= 2*(34+38)) {
    AliWarning(Form("Unexpected number of survey points %d, should be 144",ladderPoints->GetEntries())); 
  }
  Int_t iLadd = 0;
  for (Int_t ilayer =  4; ilayer <=  5; ilayer ++) {
    Int_t nLadder = 34; // layer 5
    if (ilayer == 5)
      nLadder = 38;     // layer 6

    for (Int_t iLadder = 0; iLadder < nLadder; iLadder++) {
      TString ladName = ssdName;
      ladName += ilayer;
      ladName += "/Ladder";
      ladName += iLadder;

      AliSurveyPoint *vPoint =  (AliSurveyPoint*) ladderPoints->FindObject(ladName+"/V");
      AliSurveyPoint *qPoint =  (AliSurveyPoint*) ladderPoints->FindObject(ladName+"/Q");
      if (vPoint == 0) {
	AliWarning(Form("Cannot find V side point for ladder %s",ladName.Data()));
	continue;
      }
      if (qPoint == 0) {
	AliWarning(Form("Cannot find Q side point for ladder %s",ladName.Data()));
	continue;
      }

      TString tmpStr;
      tmpStr.Insert(0,vPoint->GetName(),3);
      Int_t ladder = tmpStr.Atoi();
      tmpStr="";
      tmpStr.Insert(0,qPoint->GetName(),3);
      if (tmpStr.Atoi() != ladder) 
	AliError(Form("Survey data file error. Expect pairs of V,Q points. Got ladders %d %d",ladder,tmpStr.Atoi()));

      // Note: file gives meas-nom in local offline coordinates, 
      // ie. local z = - global z and local x = - global x (for ladder 508, i.e. top ladder)
      Double_t dxLoc = vPoint->GetX() * kMu2Cm;
      Double_t dyLoc = vPoint->GetY() * kMu2Cm;
      Double_t dzLoc = vPoint->GetZ() * kMu2Cm;

      // rot around z-axis
      Double_t phi = 0;  // Not measured
      // rot around y-axis
      Double_t theta = 0;
      Double_t psi = 0;

      // Note: local psi = -global psi, psi = atan(-(y(z1) - y(z0)) / (z1-z0))  
      // local theta = global theta = atan(dx/dz) 
      // V side is A side is large global z 
      // Q side is C side is large local z

      if (ladder >= 600) {
	theta = TMath::ATan((qPoint->GetX() - vPoint->GetX())*kMu2Cm/kLaddLen6);
	psi = TMath::ATan((vPoint->GetY() - qPoint->GetY())*kMu2Cm/kLaddLen6);
      }
      else {
	theta = TMath::ATan((qPoint->GetX() - vPoint->GetX())*kMu2Cm/kLaddLen5);
	psi = TMath::ATan((vPoint->GetY() - qPoint->GetY())*kMu2Cm/kLaddLen5);
      } 

      // Move along ladder to local Z = 0 point
      dxLoc += zLag*theta;
      dyLoc -= zLag*psi;

      // Convert to degrees
      theta *= kRadToDeg;
      psi *= kRadToDeg;
      AliDebug(1,Form("ladname %f %f %f %f %f %f ",dxLoc,dyLoc,dzLoc,psi,theta,phi));  
      
      new((*fAlignObjArray)[500+1698+iLadd]) AliAlignObjParams(ladName,0,dxLoc,dyLoc,dzLoc,psi,theta,phi,kFALSE);

      iLadd++;
    }  // Ladder loop
  }  // Layer loop
}

//______________________________________________________________________
void AliITSSurveyToAlign::CreateAlignObjDummySSDModules(){
  // 
  // Create empty alignment objects
  // Used when fSurveySSD == 0
  //
  for(Int_t imod = 0; imod < 1698; imod++) {
    Int_t ilayer = (imod < 748) ? AliGeomManager::kSSD1 : AliGeomManager::kSSD2;
    Int_t imodule = (imod < 748) ? imod : imod - 748;

    Int_t uid = AliGeomManager::LayerToVolUID(ilayer,imodule);
    const Char_t *symname = AliGeomManager::SymName(uid);

    new((*fAlignObjArray)[500+imod]) AliAlignObjParams(symname, uid, 0., 0., 0., 0., 0., 0., kTRUE);
  }//module loop
}


//______________________________________________________________________
void AliITSSurveyToAlign::GetIdPosSDD(Int_t uid, Int_t layer, Int_t module, Int_t iPoint)
{
  // 
  //    Utility function used by CreateAlignObjSDD
  // 
  TGeoHMatrix gMod = *AliGeomManager::GetMatrix(uid); //global matrix of sensor
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntryByUID(uid);
  // TString ladderPath = AliGeomManager::SymName(uid);
  TString ladderPath(pne->GetTitle());
  if(ladderPath.EndsWith("/")) ladderPath.Remove(TString::kTrailing,'/');
  ladderPath.Remove(ladderPath.Last('/'));
  ladderPath.Remove(ladderPath.Last('/'));
  gGeoManager->cd(ladderPath.Data());
  TGeoHMatrix gLad = *gGeoManager->GetCurrentMatrix(); // global matrix of ladder
  TGeoHMatrix rel = gMod; // to equal relative matrix ladder to sensor.
  TGeoHMatrix invgLad = gLad.Inverse();
  rel.MultiplyLeft(&invgLad);
  TGeoRotation* rr = new TGeoRotation("rr",90,90,0,0,90,180);
  TGeoCombiTrans* ct = 0;
  if(layer==3) ct= new TGeoCombiTrans(25.,0.,0.,rr);
  if(layer==4) ct= new TGeoCombiTrans(25.+7.5,0.,0.,rr);

  rel.MultiplyLeft(ct);
  
  if((layer==3)&&(module<3)) rel.LocalToMaster(fgkLocR[iPoint],fSDDidP[iPoint]);
  if((layer==3)&&(module>2)) rel.LocalToMaster(fgkLocL[iPoint],fSDDidP[iPoint]);
  if((layer==4)&&(module<4)) rel.LocalToMaster(fgkLocR[iPoint],fSDDidP[iPoint]);
  if((layer==4)&&(module>3)) rel.LocalToMaster(fgkLocL[iPoint],fSDDidP[iPoint]);

  for(Int_t i=0; i<3; i++) fSDDidP[iPoint][i]*=10;

}

//______________________________________________________________________
void AliITSSurveyToAlign::ReadPointNameSDD(const char str[], Int_t &iLayer, Int_t &iLader, Int_t &iModul, Int_t &iPoint) const
{
  // 
  //    Utility function used by CreateAlignObjSDD
  // 
  iLayer=-1;
  iLader=-1;
  iModul=-1;
  iPoint=-1;

  if(str[7]=='2') iLayer=3;
  if(str[7]=='3') iLayer=4;

  if(str[15]=='0') iLader=0;
  if(str[15]=='1') iLader=1;
  if(str[15]=='2') iLader=2;
  if(str[15]=='3') iLader=3;
  if(str[15]=='4') iLader=4;
  if(str[15]=='5') iLader=5;
  if(str[15]=='6') iLader=6;
  if(str[15]=='7') iLader=7;
  if(str[15]=='8') iLader=8;
  if(str[15]=='9') iLader=9;

  Int_t ord=0;
  if(str[16]=='0') {iLader=10*iLader+0; ord=1;}
  if(str[16]=='1') {iLader=10*iLader+1; ord=1;}
  if(str[16]=='2') {iLader=10*iLader+2; ord=1;}
  if(str[16]=='3') {iLader=10*iLader+3; ord=1;}
  if(str[16]=='4') {iLader=10*iLader+4; ord=1;}
  if(str[16]=='5') {iLader=10*iLader+5; ord=1;}
  if(str[16]=='6') {iLader=10*iLader+6; ord=1;}
  if(str[16]=='7') {iLader=10*iLader+7; ord=1;}
  if(str[16]=='8') {iLader=10*iLader+8; ord=1;}
  if(str[16]=='9') {iLader=10*iLader+9; ord=1;}

  if(str[23+ord]=='0') iModul=0;
  if(str[23+ord]=='1') iModul=1;
  if(str[23+ord]=='2') iModul=2;
  if(str[23+ord]=='3') iModul=3;
  if(str[23+ord]=='4') iModul=4;
  if(str[23+ord]=='5') iModul=5;
  if(str[23+ord]=='6') iModul=6;
  if(str[23+ord]=='7') iModul=7;
  if(str[23+ord]=='8') iModul=8;
  if(str[23+ord]=='9') iModul=9;

  if((str[25+ord]=='R')&&(str[26+ord]=='D')) iPoint=0;
  if((str[25+ord]=='R')&&(str[26+ord]=='C')) iPoint=1;
  if((str[25+ord]=='R')&&(str[26+ord]=='U')) iPoint=2;
  if((str[25+ord]=='L')&&(str[26+ord]=='U')) iPoint=3;
  if((str[25+ord]=='L')&&(str[26+ord]=='C')) iPoint=4;
  if((str[25+ord]=='L')&&(str[26+ord]=='D')) iPoint=5;
  return;
}


//______________________________________________________________________
void AliITSSurveyToAlign::ConvertToRSofModulesAndRotSDD(Int_t Layer, Int_t Module)
{
  // 
  //    Utility function used by CreateAlignObjSDD
  // 

  Double_t ymId;
  Double_t zmId;

  Double_t ymMe;
  Double_t zmMe;
  Double_t ymMeE;
  Double_t zmMeE;

  Double_t x0=fSDDidP[1][0];
  Double_t z0=fSDDidP[1][2]-0.52;
  for(Int_t i=0; i<6; i++)
    {
      fSDDidP[i][2]-=0.52;

      if(!fSDDisMe[i]) continue; 

      fSDDidP[i][0]-=x0;
      fSDDidP[i][2]-=z0;
      fSDDmeP[i][0]-=x0;
      fSDDmeP[i][2]-=z0;
				
      ymId=fSDDidP[i][1];
      zmId=fSDDidP[i][2];
			
      fSDDidP[i][2]=fSDDidP[i][0];
      fSDDidP[i][0]=ymId;
      fSDDidP[i][1]=zmId;
			
      ymMe=fSDDmeP[i][1];
      zmMe=fSDDmeP[i][2];
			
      ymMeE=fSDDmeP[i][4];
      zmMeE=fSDDmeP[i][5];
			
      fSDDmeP[i][2]=fSDDmeP[i][0];
      fSDDmeP[i][0]=ymMe;
      fSDDmeP[i][1]=zmMe;
      fSDDmeP[i][5]=fSDDmeP[i][3];
      fSDDmeP[i][3]=ymMeE;
      fSDDmeP[i][4]=zmMeE;
			

      if(((Layer==3)&&(Module>2))||((Layer==4)&&(Module>3)))
	{
	  fSDDidP[i][0]*=(-1);
	  fSDDidP[i][2]*=(-1);
	  fSDDmeP[i][0]*=(-1);
	  fSDDmeP[i][2]*=(-1);
	}
    }	
}


//______________________________________________________________________
void AliITSSurveyToAlign::CalcShiftSDD(Double_t &x0,Double_t &y0,Double_t &z0) const
{
    // Calculates the 3 shifts for the present SDD module
    // and sets the three reference arguments
    //
  Double_t xId, yId, zId;
  Double_t xMe, yMe, zMe, sX2, sY2, sZ2;
  Double_t aX=0., bX=0.;
  Double_t aY=0., bY=0.;
  Double_t aZ=0., bZ=0.;
  for(Int_t iP1=0; iP1<6; iP1++)
    {
      if(!fSDDisMe[iP1]) continue;
      xId=fSDDidP[iP1][0];
      yId=fSDDidP[iP1][1];
      zId=fSDDidP[iP1][2];
      xMe=fSDDmeP[iP1][0];
      yMe=fSDDmeP[iP1][1];
      zMe=fSDDmeP[iP1][2];
      sX2 =fSDDmeP[iP1][3]*fSDDmeP[iP1][3];
      sY2 =fSDDmeP[iP1][4]*fSDDmeP[iP1][4];
      sZ2 =fSDDmeP[iP1][5]*fSDDmeP[iP1][5];
      aX+=(1./sX2);
      bX+=((xMe-xId)/sX2); 
      aY+=(1./sY2);
      bY+=((yMe-yId)/sY2); 
      aZ+=(1./sZ2);
      bZ+=((zMe-zId)/sZ2); 
    }
  Double_t x1 = bX/aX;
  Double_t x2 = bY/aY;
  Double_t x3 = bZ/aZ;
  x0=x1;
  y0=x2;
  z0=x3;
  return;
}


//______________________________________________________________________
void AliITSSurveyToAlign::CalcShiftRotSDD(Double_t &tet,Double_t &psi,Double_t &phi,Double_t &x0,Double_t &y0,Double_t &z0)
{
    // Calculates the 3 shifts and 3 euler angles for the present SDD module
    // and sets the six reference arguments
    //
  TMatrixD pC(6,6);

  Double_t a[6][6];
  for(Int_t ii=0; ii<6; ii++){
      for(Int_t jj=0; jj<6; jj++){
	  a[ii][jj]=0.;
      }
  }

  Double_t c[6];
  for(Int_t ii=0; ii<6; ii++)
      c[ii]=0.;

  Double_t xId, yId, zId;
  Double_t xMe, yMe, zMe, sX2, sY2, sZ2;

  for(Int_t iP1=0; iP1<=6; iP1++)
    {
      if(!fSDDisMe[iP1]) continue;

      //ideal x,y,z for fiducial mark iP1
      xId= fSDDidP[iP1][0];
      yId= fSDDidP[iP1][1];
      zId= fSDDidP[iP1][2];

      //measured x,y,z for fiducial mark iP1
      xMe= fSDDmeP[iP1][0];
      yMe= fSDDmeP[iP1][1];
      zMe= fSDDmeP[iP1][2];

      //squared precisions of measured x,y,z for fiducial mark iP1
      sX2 = fSDDmeP[iP1][3]* fSDDmeP[iP1][3];
      sY2 = fSDDmeP[iP1][4]* fSDDmeP[iP1][4];
      sZ2 = fSDDmeP[iP1][5]* fSDDmeP[iP1][5];

      a[0][0]+=(zId*zId/sX2+xId*xId/sZ2);
      a[0][1]-=(zId*yId/sX2);
      a[0][2]-=(xId*yId/sZ2);
      a[0][3]-=(zId/sX2);
      a[0][4] =0.;
      a[0][5]+=(xId/sZ2);
      c[0]+=(xId*(zMe-zId)/sZ2-zId*(xMe-xId)/sX2); 

      a[1][0]-=(yId*zId/sX2);
      a[1][1]+=(xId*xId/sY2+yId*yId/sX2);
      a[1][2]-=(xId*zId/sY2);
      a[1][3]+=(yId/sX2);
      a[1][4]-=(xId/sY2);
      a[1][5] =0.;
      c[1]+=(yId*(xMe-xId)/sX2-xId*(yMe-yId)/sY2); 

      a[2][0]-=(yId*xId/sZ2);
      a[2][1]-=(xId*zId/sY2);
      a[2][2]+=(zId*zId/sY2+yId*yId/sZ2);
      a[2][3] =0.;
      a[2][4]+=(zId/sY2);
      a[2][5]-=(yId/sZ2);
      c[2]+=(zId*(yMe-yId)/sY2-yId*(zMe-zId)/sZ2); 

      a[3][0]-=(zId/sX2);
      a[3][1]+=(yId/sX2);
      a[3][2] =0.;
      a[3][3]+=(1./sX2);
      a[3][4] =0.;
      a[3][5] =0.;
      c[3]+=((xMe-xId)/sX2); 

      a[4][0] =0.;
      a[4][1]-=(xId/sY2);
      a[4][2]+=(zId/sY2);
      a[4][3] =0.;
      a[4][4]+=(1./sY2);
      a[4][5] =0.;
      c[4]+=((yMe-yId)/sY2); 

      a[5][0]+=(xId/sZ2);
      a[5][1] =0.;
      a[5][2]-=(yId/sZ2);
      a[5][3] =0.;
      a[5][4] =0.;
      a[5][5]+=(1./sZ2);
      c[5]+=((zMe-zId)/sZ2); 
    }

  ///////////////////////////////////////////////////////////////

  pC.SetMatrixArray(&(a[0][0]));
  TMatrixD p1(pC);
  TMatrixD p2(pC);
  TMatrixD p3(pC);
  TMatrixD p4(pC);
  TMatrixD p5(pC);
  TMatrixD p6(pC);

  for(Int_t raw=0; raw<6; raw++)
      p1[raw][0]=c[raw];
  for(Int_t raw=0; raw<6; raw++)
      p2[raw][1]=c[raw];
  for(Int_t raw=0; raw<6; raw++)
      p3[raw][2]=c[raw];
  for(Int_t raw=0; raw<6; raw++)
      p4[raw][3]=c[raw];
  for(Int_t raw=0; raw<6; raw++)
      p5[raw][4]=c[raw];
  for(Int_t raw=0; raw<6; raw++)
      p6[raw][5]=c[raw];

  // cout << "calculating determinants" << endl;
  Double_t det0=pC.Determinant();
  Double_t x1 = p1.Determinant()/det0;
  Double_t x2 = p2.Determinant()/det0;
  Double_t x3 = p3.Determinant()/det0;
  Double_t x4 = p4.Determinant()/det0;
  Double_t x5 = p5.Determinant()/det0;
  Double_t x6 = p6.Determinant()/det0;
  //cout << "calculating determinants done" << endl;
  if (x1 == 0) {
    AliInfo("p1 singular ");
    p1.Print();
  }
  if (x2 == 0) {
    AliInfo("p2 singular ");
    p2.Print();
  }
  if (x3 == 0) {
    AliInfo("p3 singular ");
    p3.Print();
  }
  if (x4 == 0) {
    AliInfo("p4 singular ");
    p4.Print();
  }
  if (x5 == 0) {
    AliInfo("p5 singular ");
    p5.Print();
  }
  if (x6 == 0) {
    AliInfo("p6 singular ");
    p6.Print();
  }


  tet=x1;
  psi=x2;
  phi=x3;
  x0=x4;
  y0=x5;
  z0=x6;
  return;
}

