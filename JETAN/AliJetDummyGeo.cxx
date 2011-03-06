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

//
// Temporarily added to define part of the EMCal geometry
// necessary for the jet finder
//
// M. Estienne
// Magali.Estienne@cern.ch
//

#include <Riostream.h>
#include <assert.h>
#include <TList.h>

// --- Root header files ---
#include <TVector3.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliJetDummyGeo.h"
#include "AliJetDummyShishKebabTrd1Module.h"

ClassImp(AliJetDummyGeo)

AliJetDummyGeo::AliJetDummyGeo():
  fArm1EtaMin(-0.7), 
  fArm1EtaMax(+0.7), 
  fArm1PhiMin(80.0), 
  fArm1PhiMax(200.0),
  fNumberOfSuperModules(12), 
  fSteelFrontThick(0.0), 
  fLateralSteelStrip(0.01),
  fIPDistance(428.0), 
  fPhiGapForSM(2.), 
  fNPhi(12),       
  fNZ(24),      
  fPhiModuleSize(12.26 - fPhiGapForSM / Float_t(fNPhi)), 
  fEtaModuleSize(fPhiModuleSize),
  fNPHIdiv(2), 
  fNETAdiv(2), 
  fPhiTileSize(fPhiModuleSize/Double_t(fNPHIdiv) - fLateralSteelStrip),
  fEtaTileSize(fEtaModuleSize/Double_t(fNETAdiv) - fLateralSteelStrip),
  fNECLayers(77), 
  fECScintThick(0.16), 
  fECPbRadThickness(0.16), 
  fSampling(12.327),
  fTrd1Angle(1.5), 
  fNCellsInModule(fNPHIdiv*fNETAdiv), 
  fNCellsInSupMod(fNCellsInModule*fNPhi*fNZ), 
  fNCells(fNCellsInSupMod*fNumberOfSuperModules-fNCellsInSupMod), 
  fLongModuleSize(fNECLayers*(fECScintThick + fECPbRadThickness)), 
  f2Trd1Dx2(fEtaModuleSize + 2.*fLongModuleSize*TMath::Tan(fTrd1Angle*TMath::DegToRad()/2.)), 
  fShellThickness(TMath::Sqrt(fLongModuleSize*fLongModuleSize + f2Trd1Dx2*f2Trd1Dx2)+fSteelFrontThick),
  fZLength(2.*ZFromEtaR(fIPDistance+fShellThickness,fArm1EtaMax)),
  fEtaMaxOfTRD1(0.), 
  fPhiBoundariesOfSM(0),
  fPhiCentersOfSM(0),
  fCentersOfCellsEtaDir(0),
  fCentersOfCellsXDir(0),
  fCentersOfCellsPhiDir(0),
  fEtaCentersOfCells(0),
  fPhiCentersOfCells(0),
  fShishKebabTrd1Modules(0),
  fDebug(0)
{
  // Constructor

  // Local coordinates
  fParSM[0] = GetShellThickness()/2.;        
  fParSM[1] = GetPhiModuleSize() * GetNPhi()/2.;
  fParSM[2] = 350./2.;

  fEnvelop[0] = fIPDistance; // mother volume inner radius
  fEnvelop[1] = fIPDistance + fShellThickness; // mother volume outer r.
  fEnvelop[2] = 1.00001*fZLength; // add some padding for mother volume. 

  // SM phi boundaries - (0,1),(2,3) .. (10,11) - has the same boundaries; Nov 7, 2006 
  fPhiBoundariesOfSM.Set(fNumberOfSuperModules);
  fPhiCentersOfSM.Set(fNumberOfSuperModules/2);
  fPhiBoundariesOfSM[0] = TMath::PiOver2() - TMath::ATan2(fParSM[1] , fIPDistance); // 1th and 2th modules)
  fPhiBoundariesOfSM[1] = TMath::PiOver2() + TMath::ATan2(fParSM[1] , fIPDistance);
  fPhiCentersOfSM[0]    = TMath::PiOver2();
  for(Int_t i=1; i<=4; i++) { // from 2th ro 9th
    fPhiBoundariesOfSM[2*i]   = fPhiBoundariesOfSM[0] + 20.*TMath::DegToRad()*i;
    fPhiBoundariesOfSM[2*i+1] = fPhiBoundariesOfSM[1] + 20.*TMath::DegToRad()*i;
    fPhiCentersOfSM[i]        = fPhiCentersOfSM[0]    + 20.*TMath::DegToRad()*i;
  }
  fPhiBoundariesOfSM[11] = 190.*TMath::DegToRad();
  fPhiBoundariesOfSM[10] = fPhiBoundariesOfSM[11] - TMath::ATan2((fParSM[1]) , fIPDistance);
  fPhiCentersOfSM[5]     = (fPhiBoundariesOfSM[10]+fPhiBoundariesOfSM[11])/2.; 

  Int_t nphism  = GetNumberOfSuperModules()/2;
  Double_t dphi = (GetArm1PhiMax() - GetArm1PhiMin())/nphism;
  Double_t rpos = (GetEnvelop(0) + GetEnvelop(1))/2.;
  Double_t phi, phiRad, xpos, ypos, zpos;
  for(Int_t i=0; i<nphism; i++){
    phi    = GetArm1PhiMin() + dphi*(2*i+1)/2.; // phi= 90, 110, 130, 150, 170, 190
    phiRad = phi*TMath::Pi()/180.;
    xpos = rpos * TMath::Cos(phiRad);
    ypos = rpos * TMath::Sin(phiRad);
    zpos = fParSM[2];
    if(i==5) {
      xpos += (fParSM[1]/2. * TMath::Sin(phiRad)); 
      ypos -= (fParSM[1]/2. * TMath::Cos(phiRad));
    }
    // pozitive z
    Int_t ind = 2*i;
    TGeoRotation *geoRot0 = new TGeoRotation("geoRot0", 90.0, phi, 90.0, 90.0+phi, 0.0, 0.0);
    fMatrixOfSM[ind] = new TGeoCombiTrans(Form("EmcalSM%2.2i",ind),
					  xpos,ypos, zpos, geoRot0);
    // negaive z
    ind++;
    Double_t phiy = 90. + phi + 180.;
    if(phiy>=360.) phiy -= 360.;
    TGeoRotation *geoRot1 = new TGeoRotation("geoRot1", 90.0, phi, 90.0, phiy, 180.0, 0.0);
    fMatrixOfSM[ind] = new TGeoCombiTrans(Form("EmcalSM%2.2i",ind),
					  xpos,ypos,-zpos, geoRot1);
  } // for

  CreateListOfTrd1Modules();

  if(fDebug > 0){
    for(Int_t i=0; i<6; i++){
      cout << "fMatrixOfSM[" << i << "]: " << fMatrixOfSM[i] << endl;
    }
    cout << "fArm1EtaMin: " << fArm1EtaMin << endl;
    cout << "fArm1EtaMax: " << fArm1EtaMax << endl;
    cout << "fArm1PhiMin: " << fArm1PhiMin << endl;
    cout << "fArm1PhiMax: " << fArm1PhiMax << endl;
    cout << "fNumberOfSuperModules: " << fNumberOfSuperModules << endl;
    cout << "fSteelFrontThick: " << fSteelFrontThick << endl;
    cout << "fIPDistance: " << fIPDistance  << endl;
    cout << "fZLength: " << fZLength << endl;
    cout << "fPhiGapForSM: " << fPhiGapForSM << endl;
    cout << "fNPhi: " << fNPhi << endl;
    cout << "fNZ: " << fNZ << endl;
    cout << "fPhiModuleSize: " << fPhiModuleSize << endl;
    cout << "fEtaModuleSize: " << fEtaModuleSize << endl;
    cout << "fNPHIdiv: " << fNPHIdiv << endl;
    cout << "fNETAdiv: " << fNETAdiv << endl;
    cout << "fNECLayers: " << fNECLayers << endl;
    cout << "fECScintThick: " <<  fECScintThick << endl;
    cout << "fECPbRadThickness: " << fECPbRadThickness << endl;
    cout << "fSampling: " << fSampling << endl;
    cout << "fTrd1Angle: " << fTrd1Angle << endl;
    cout << "fNCellsInModule: " << fNCellsInModule << endl;
    cout << "fNCellsInSupMod: " << fNCellsInSupMod << endl;
    cout << "fNCells: " << fNCells << endl;
    cout << "fLongModuleSize: " <<  fLongModuleSize << endl;
    cout << "f2Trd1Dx2: " << f2Trd1Dx2 << endl;
    cout << "fShellThickness: " << fShellThickness << endl;
    cout << "fEtaMaxOfTRD1: " << fEtaMaxOfTRD1 << endl;
  }
}

AliJetDummyGeo::AliJetDummyGeo(const AliJetDummyGeo& geom):
  TObject(geom),  
  fArm1EtaMin(geom.fArm1EtaMin), 
  fArm1EtaMax(geom.fArm1EtaMax), 
  fArm1PhiMin(geom.fArm1PhiMin), 
  fArm1PhiMax(geom.fArm1PhiMax),
  fNumberOfSuperModules(geom.fNumberOfSuperModules), 
  fSteelFrontThick(geom.fSteelFrontThick), 
  fLateralSteelStrip(geom.fLateralSteelStrip),
  fIPDistance(geom.fIPDistance), 
  fPhiGapForSM(geom.fPhiGapForSM), 
  fNPhi(geom.fNPhi),       
  fNZ(geom.fNZ),      
  fPhiModuleSize(geom.fPhiModuleSize), 
  fEtaModuleSize(geom.fEtaModuleSize),
  fNPHIdiv(geom.fNPHIdiv), 
  fNETAdiv(geom.fNETAdiv), 
  fPhiTileSize(geom.fPhiTileSize),
  fEtaTileSize(geom.fEtaTileSize),
  fNECLayers(geom.fNECLayers), 
  fECScintThick(geom.fECScintThick), 
  fECPbRadThickness(geom.fECPbRadThickness), 
  fSampling(geom.fSampling),
  fTrd1Angle(geom.fTrd1Angle), 
  fNCellsInModule(geom.fNCellsInModule), 
  fNCellsInSupMod(geom.fNCellsInSupMod), 
  fNCells(geom.fNCells),
  fLongModuleSize(geom.fLongModuleSize),
  f2Trd1Dx2(geom.f2Trd1Dx2),
  fShellThickness(geom.fShellThickness),
  fZLength(geom.fZLength),
  fEtaMaxOfTRD1(geom.fEtaMaxOfTRD1),
  fPhiBoundariesOfSM(geom.fPhiBoundariesOfSM),
  fPhiCentersOfSM(geom.fPhiCentersOfSM),
  fCentersOfCellsEtaDir(geom.fCentersOfCellsEtaDir),
  fCentersOfCellsXDir(geom.fCentersOfCellsXDir),
  fCentersOfCellsPhiDir(geom.fCentersOfCellsPhiDir),
  fEtaCentersOfCells(geom.fEtaCentersOfCells),
  fPhiCentersOfCells(geom.fPhiCentersOfCells),
  fShishKebabTrd1Modules(geom.fShishKebabTrd1Modules),
  fDebug(geom.fDebug)
{
  // Constructor
  // Local coordinates
  fParSM[0] = GetShellThickness()/2.;        
  fParSM[1] = GetPhiModuleSize() * GetNPhi()/2.;
  fParSM[2] = 350./2.;

  fEnvelop[0] = fIPDistance; // mother volume inner radius
  fEnvelop[1] = fIPDistance + fShellThickness; // mother volume outer r.
  fEnvelop[2] = 1.00001*fZLength; // add some padding for mother volume. 

  // SM phi boundaries - (0,1),(2,3) .. (10,11) - has the same boundaries; Nov 7, 2006 
  fPhiBoundariesOfSM.Set(fNumberOfSuperModules);
  fPhiCentersOfSM.Set(fNumberOfSuperModules/2);
  fPhiBoundariesOfSM[0] = TMath::PiOver2() - TMath::ATan2(fParSM[1] , fIPDistance); // 1th and 2th modules)
  fPhiBoundariesOfSM[1] = TMath::PiOver2() + TMath::ATan2(fParSM[1] , fIPDistance);
  fPhiCentersOfSM[0]    = TMath::PiOver2();
  for(Int_t i=1; i<=4; i++) { // from 2th ro 9th
    fPhiBoundariesOfSM[2*i]   = fPhiBoundariesOfSM[0] + 20.*TMath::DegToRad()*i;
    fPhiBoundariesOfSM[2*i+1] = fPhiBoundariesOfSM[1] + 20.*TMath::DegToRad()*i;
    fPhiCentersOfSM[i]        = fPhiCentersOfSM[0]    + 20.*TMath::DegToRad()*i;
  }
  fPhiBoundariesOfSM[11] = 190.*TMath::DegToRad();
  fPhiBoundariesOfSM[10] = fPhiBoundariesOfSM[11] - TMath::ATan2((fParSM[1]) , fIPDistance);
  fPhiCentersOfSM[5]     = (fPhiBoundariesOfSM[10]+fPhiBoundariesOfSM[11])/2.; 

  Int_t nphism  = GetNumberOfSuperModules()/2;
  Double_t dphi = (GetArm1PhiMax() - GetArm1PhiMin())/nphism;
  Double_t rpos = (GetEnvelop(0) + GetEnvelop(1))/2.;
  Double_t phi, phiRad, xpos, ypos, zpos;
  for(Int_t i=0; i<nphism; i++){
    phi    = GetArm1PhiMin() + dphi*(2*i+1)/2.; // phi= 90, 110, 130, 150, 170, 190
    phiRad = phi*TMath::Pi()/180.;
    xpos = rpos * TMath::Cos(phiRad);
    ypos = rpos * TMath::Sin(phiRad);
    zpos = fParSM[2];
    if(i==5) {
      xpos += (fParSM[1]/2. * TMath::Sin(phiRad)); 
      ypos -= (fParSM[1]/2. * TMath::Cos(phiRad));
    }
    // pozitive z
    Int_t ind = 2*i;
    TGeoRotation *geoRot0 = new TGeoRotation("geoRot0", 90.0, phi, 90.0, 90.0+phi, 0.0, 0.0);
    fMatrixOfSM[ind] = new TGeoCombiTrans(Form("EmcalSM%2.2i",ind),
					  xpos,ypos, zpos, geoRot0);
    // negaive z
    ind++;
    Double_t phiy = 90. + phi + 180.;
    if(phiy>=360.) phiy -= 360.;
    TGeoRotation *geoRot1 = new TGeoRotation("geoRot1", 90.0, phi, 90.0, phiy, 180.0, 0.0);
    fMatrixOfSM[ind] = new TGeoCombiTrans(Form("EmcalSM%2.2i",ind),
					  xpos,ypos,-zpos, geoRot1);

    delete geoRot0;
    delete geoRot1;

  } // for

  CreateListOfTrd1Modules();

}

//------------------------------------------------------------------------------------
AliJetDummyGeo::~AliJetDummyGeo()
{
  // Destructor
  // delete [] fMatrixOfSM;
}

//------------------------------------------------------------------------------------
void AliJetDummyGeo::EtaPhiFromIndex(Int_t absId, Float_t& eta, Float_t& phi) const
{
  // Nov 16, 2006- float to double
  // version for TRD1 only
  static TVector3 vglob;
  GetGlobal(absId, vglob);
  eta = vglob.Eta();
  phi = vglob.Phi();

}

//------------------------------------------------------------------------------------
void AliJetDummyGeo::GetGlobal(const Double_t *loc, Double_t *glob, Int_t ind) const
{
  // Figure out the global numbering of a given supermodule from the
  // local numbering
  // Alice numbering - Jun 03,2006
  //  if(fMatrixOfSM[0] == 0) GetTransformationForSM();

  if(ind>=0 && ind < GetNumberOfSuperModules()) {
    fMatrixOfSM[ind]->LocalToMaster(loc, glob);
  }
}

//------------------------------------------------------------------------------------
void AliJetDummyGeo::GetGlobal(Int_t absId , Double_t glob[3]) const
{ 
  // Alice numbering scheme - Jun 03, 2006
  static Int_t nSupMod, nModule, nIphi, nIeta;
  static Double_t loc[3];

  glob[0]=glob[1]=glob[2]=0.0; // bad case
  if(RelPosCellInSModule(absId, loc)) {
    GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
    fMatrixOfSM[nSupMod]->LocalToMaster(loc, glob);
  }
}

//------------------------------------------------------------------------------------
void AliJetDummyGeo::GetGlobal(Int_t absId , TVector3 &vglob) const
{ 
  // Alice numbering scheme - Jun 03, 2006
  static Double_t glob[3];

  GetGlobal(absId, glob);
  vglob.SetXYZ(glob[0], glob[1], glob[2]);

}

//------------------------------------------------------------------------------------
Bool_t AliJetDummyGeo::RelPosCellInSModule(Int_t absId, Double_t loc[3]) const
{
  // Alice numbering scheme - Jun 03, 2006
  loc[0] = loc[1] = loc[2]=0.0;
  if(RelPosCellInSModule(absId, loc[0],loc[1],loc[2])) {
    return kTRUE;
  }
  return kFALSE;
}

//------------------------------------------------------------------------------------
Bool_t AliJetDummyGeo::RelPosCellInSModule(Int_t absId, Double_t &xr, Double_t &yr, Double_t &zr) const
{
  // Look to see what the relative position inside a given cell is
  // for a recpoint.
  // Alice numbering scheme - Jun 08, 2006
  // In:
  // absId   - cell is as in Geant,     0<= absId   < fNCells;
  // OUT:
  // xr,yr,zr - x,y,z coordinates of cell with absId inside SM 

  // Shift index taking into account the difference between standard SM 
  // and SM of half size in phi direction
  const Int_t phiIndexShift = fCentersOfCellsPhiDir.GetSize()/4; // Nov 22, 2006; was 6 for cas 2X2
  static Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta;
  if(!CheckAbsCellId(absId)) return kFALSE;

  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
  GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi, ieta); 
 
  xr = fCentersOfCellsXDir.At(ieta);
  zr = fCentersOfCellsEtaDir.At(ieta);

  if(nSupMod<10) {
    yr = fCentersOfCellsPhiDir.At(iphi);
  } else {
    yr = fCentersOfCellsPhiDir.At(iphi + phiIndexShift);
  }

  return kTRUE;
}

//------------------------------------------------------------------------------------
Bool_t  AliJetDummyGeo::CheckAbsCellId(Int_t absId) const
{ 
  // May 31, 2006; only trd1 now
  if(absId<0 || absId >= fNCells) return kFALSE;
  else                            return kTRUE;
}

//------------------------------------------------------------------------------------
Bool_t AliJetDummyGeo::GetCellIndex(Int_t absId,Int_t &nSupMod,Int_t &nModule,Int_t &nIphi,Int_t &nIeta) const
{ 
  // 21-sep-04; 19-oct-05;
  // May 31, 2006; ALICE numbering scheme:
  // 
  // In:
  // absId   - cell is as in Geant,     0<= absId   < fNCells;
  // Out:
  // nSupMod - super module(SM) number, 0<= nSupMod < fNumberOfSuperModules;
  // nModule  - module number in SM,    0<= nModule  < fNCellsInSupMod/fNCellsInSupMod or(/2) for tow last SM (10th and 11th);
  // nIphi   - cell number in phi driection inside module; 0<= nIphi < fNPHIdiv; 
  // nIeta   - cell number in eta driection inside module; 0<= nIeta < fNETAdiv; 
  // 
  static Int_t tmp=0, sm10=0;
  if(!CheckAbsCellId(absId)) return kFALSE;

  sm10 = fNCellsInSupMod*10;
  if(absId >= sm10) { // 110 degree case; last two supermodules  
    nSupMod = (absId-sm10) / (fNCellsInSupMod/2) + 10;
    tmp     = (absId-sm10) % (fNCellsInSupMod/2);
  } else {
    nSupMod = absId / fNCellsInSupMod;
    tmp     = absId % fNCellsInSupMod;
  }

  nModule  = tmp / fNCellsInModule;
  tmp     = tmp % fNCellsInModule;
  nIphi   = tmp / fNPHIdiv;
  nIeta   = tmp % fNPHIdiv;

  return kTRUE;
}

//------------------------------------------------------------------------------------
void AliJetDummyGeo::GetCellPhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta, Int_t &iphi, Int_t &ieta) const
{ 
  // 
  // Added nSupMod; Nov 25, 05
  // Alice numbering scheme  - Jun 01,2006 
  // IN:
  // nSupMod - super module(SM) number, 0<= nSupMod < fNumberOfSuperModules;
  // nModule  - module number in SM,     0<= nModule  < fNCellsInSupMod/fNCellsInSupMod or(/2) for tow last SM (10th and 11th);
  // nIphi   - cell number in phi driection inside module; 0<= nIphi < fNPHIdiv; 
  // nIeta   - cell number in eta driection inside module; 0<= nIeta < fNETAdiv; 
  // 
 // OUT:
  // ieta, iphi - indexes of cell(tower) in two dimensional grid of SM
  // ieta - have to change from 0 to (fNZ*fNETAdiv-1)
  // iphi - have to change from 0 to (fNPhi*fNPHIdiv-1 or fNPhi*fNPHIdiv/2-1)
  //
  static Int_t iphim, ietam;

  GetModulePhiEtaIndexInSModule(nSupMod,nModule, iphim, ietam); 
  //  ieta  = ietam*fNETAdiv + (1-nIeta); // x(module) = -z(SM) 
  ieta  = ietam*fNETAdiv + (fNETAdiv - 1 - nIeta); // x(module) = -z(SM) 
  iphi  = iphim*fNPHIdiv + nIphi;     // y(module) =  y(SM) 

}


//------------------------------------------------------------------------------------
void AliJetDummyGeo::GetModulePhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule,  Int_t &iphim, Int_t &ietam) const
{ 
  // added nSupMod; - 19-oct-05 !
  // Alice numbering scheme        - Jun 01,2006 
  // ietam, iphi - indexes of module in two dimensional grid of SM
  // ietam - have to change from 0 to fNZ-1
  // iphim - have to change from 0 to nphi-1 (fNPhi-1 or fNPhi/2-1)
  static Int_t nphi;

  if(nSupMod>=10) nphi = fNPhi/2;
  else            nphi = fNPhi;

  ietam = nModule/nphi;
  iphim = nModule%nphi;
}

//------------------------------------------------------------------------------------
Bool_t AliJetDummyGeo::GetAbsCellIdFromEtaPhi(Double_t eta, Double_t phi, Int_t &absId) const
{
  // Nov 17,2006
  // stay here - phi problem as usual 
  static Int_t nSupMod, i, ieta, iphi, etaShift, nphi;
  static Double_t absEta=0.0, d=0.0, dmin=0.0, phiLoc;
  absId = nSupMod = - 1;
  if(SuperModuleNumberFromEtaPhi(eta, phi, nSupMod)) {
    // phi index first
    phi    = TVector2::Phi_0_2pi(phi);
    phiLoc = phi - fPhiCentersOfSM[nSupMod/2];
    nphi   = fPhiCentersOfCells.GetSize();
    if(nSupMod>=10) {
      phiLoc = phi - 190.*TMath::DegToRad();
      nphi /= 2;
    }

    dmin   = TMath::Abs(fPhiCentersOfCells[0]-phiLoc);
    iphi   = 0;
    for(i=1; i<nphi; i++) {
      d = TMath::Abs(fPhiCentersOfCells[i] - phiLoc);
      if(d < dmin) {
        dmin = d;
        iphi = i;
      }
    }
    // odd SM are turned with respect of even SM - reverse indexes

    // eta index
    absEta   = TMath::Abs(eta);
    etaShift = iphi*fCentersOfCellsEtaDir.GetSize();
    dmin     = TMath::Abs(fEtaCentersOfCells[etaShift]-absEta);
    ieta     = 0;
    for(i=1; i<fCentersOfCellsEtaDir.GetSize(); i++) {
      d = TMath::Abs(fEtaCentersOfCells[i+etaShift] - absEta);
      if(d < dmin) {
        dmin = d;
        ieta = i;
      }
    }

    if(eta<0) iphi = (nphi-1) - iphi;
    absId = GetAbsCellIdFromCellIndexes(nSupMod, iphi, ieta);

    return kTRUE;
  }
  return kFALSE;
}

//------------------------------------------------------------------------------------
Bool_t AliJetDummyGeo::SuperModuleNumberFromEtaPhi(Double_t eta, Double_t phi, Int_t &nSupMod) const
{ 
  // Return false if phi belongs a phi cracks between SM
 
  static Int_t i;

  if(TMath::Abs(eta) > fEtaMaxOfTRD1) return kFALSE;

  phi = TVector2::Phi_0_2pi(phi); // move phi to (0,2pi) boundaries
  for(i=0; i<6; i++) {
    if(phi>=fPhiBoundariesOfSM[2*i] && phi<=fPhiBoundariesOfSM[2*i+1]) {
      nSupMod = 2*i;
      if(eta < 0.0) nSupMod++;
      return kTRUE;
    }
  }
  return kFALSE;
}

//------------------------------------------------------------------------------------
Int_t  AliJetDummyGeo::GetAbsCellIdFromCellIndexes(Int_t nSupMod, Int_t iphi, Int_t ieta) const
{
  // Transition from super module number(nSupMod) and cell indexes (ieta,iphi) to absId
  static Int_t ietam, iphim, nModule;
  static Int_t nIeta, nIphi; // cell indexes in module

  GetModuleIndexesFromCellIndexesInSModule(nSupMod, iphi, ieta, ietam, iphim, nModule);

  nIeta = ieta%fNETAdiv;
  nIeta = fNETAdiv - 1 - nIeta;
  nIphi = iphi%fNPHIdiv;

  return GetAbsCellId(nSupMod, nModule, nIphi, nIeta);
}

//------------------------------------------------------------------------------------
void  AliJetDummyGeo::GetModuleIndexesFromCellIndexesInSModule(Int_t nSupMod, Int_t iphi, Int_t ieta, 
			Int_t &iphim, Int_t &ietam, Int_t &nModule) const
{
  // Transition from cell indexes (ieta,iphi) to module indexes (ietam,iphim, nModule)
  static Int_t nphi;
  nphi  = GetNumberOfModuleInPhiDirection(nSupMod);  

  ietam  = ieta/fNETAdiv;
  iphim  = iphi/fNPHIdiv;
  nModule = ietam * nphi + iphim; 
}

//------------------------------------------------------------------------------------
Int_t AliJetDummyGeo::GetAbsCellId(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta) const
{ 
  // 27-aug-04; 
  // corr. 21-sep-04; 
  //       13-oct-05; 110 degree case
  // May 31, 2006; ALICE numbering scheme:
  // 0 <= nSupMod < fNumberOfSuperModules
  // 0 <= nModule  < fNPHI * fNZ ( fNPHI * fNZ/2 for fKey110DEG=1)
  // 0 <= nIphi   < fNPHIdiv
  // 0 <= nIeta   < fNETAdiv
  // 0 <= absid   < fNCells
  static Int_t id=0; // have to change from 0 to fNCells-1
  if(nSupMod >= 10) { // 110 degree case; last two supermodules
    id  = fNCellsInSupMod*10 + (fNCellsInSupMod/2)*(nSupMod-10);
  } else {
    id  = fNCellsInSupMod*nSupMod;
  }
  id += fNCellsInModule *nModule;
  id += fNPHIdiv *nIphi;
  id += nIeta;
  if(id<0 || id >= fNCells) {
    id = -TMath::Abs(id); // if negative something wrong
  }
  return id;
}

//------------------------------------------------------------------------------------
Bool_t AliJetDummyGeo::GetPhiBoundariesOfSMGap(Int_t nPhiSec, Double_t &phiMin, Double_t &phiMax) const
{
  // 0<= nPhiSec <=4; phi in rad
  // 0;  gap boundaries between  0th&2th  | 1th&3th SM
  // 1;  gap boundaries between  2th&4th  | 3th&5th SM
  // 2;  gap boundaries between  4th&6th  | 5th&7th SM
  // 3;  gap boundaries between  6th&8th  | 7th&9th SM
  // 4;  gap boundaries between  8th&10th | 9th&11th SM
  if(nPhiSec<0 || nPhiSec >4) return kFALSE; 
  phiMin = fPhiBoundariesOfSM[2*nPhiSec+1];
  phiMax = fPhiBoundariesOfSM[2*nPhiSec+2];
  return kTRUE; 
}

//------------------------------------------------------------------------------------
void AliJetDummyGeo::CreateListOfTrd1Modules()
{
  // Generate the list of Trd1 modules
  // which will make up the EMCAL
  // geometry
    printf("CreateListOfTrd1Modules() \n");
    
  AliJetDummyShishKebabTrd1Module *mod=0, *mTmp=0; // current module
  if(fShishKebabTrd1Modules == 0) {
    fShishKebabTrd1Modules = new TList;
    fShishKebabTrd1Modules->SetName("ListOfTRD1");
    for(Int_t iz=0; iz< GetNZ(); iz++) { 
      if(iz==0) { 
        mod  = new AliJetDummyShishKebabTrd1Module(TMath::Pi()/2.,this);
      } else {
        mTmp  = new AliJetDummyShishKebabTrd1Module(*mod);
        mod   = mTmp;
      }
      fShishKebabTrd1Modules->Add(mod);
    }
  } else {
    printf(" Already exits : ");
  }
  mod = (AliJetDummyShishKebabTrd1Module*)fShishKebabTrd1Modules->At(fShishKebabTrd1Modules->GetSize()-1);
  fEtaMaxOfTRD1 = mod->GetMaxEtaOfModule(0);

  //  printf(" fShishKebabTrd1Modules has %i modules : max eta %5.4f \n", 
  //		  fShishKebabTrd1Modules->GetSize(),fEtaMaxOfTRD1);

  // Feb 20,2006;
  // Jun 01, 2006 - ALICE numbering scheme
  // define grid for cells in eta(z) and x directions in local coordinates system of SM
  // Works just for 2x2 case only -- ?? start here
  // 
  //
  // Define grid for cells in phi(y) direction in local coordinates system of SM
  // as for 2X2 as for 3X3 - Nov 8,2006
  // 

  //  printf(" Cells grid in phi directions : size %i\n", fCentersOfCellsPhiDir.GetSize());

  Int_t ind=0; // this is phi index
  Int_t ieta=0, nModule=0, iphiTemp;
  Double_t xr, zr, theta, phi, eta, r, x,y;
  xr = 0.;
  zr = 0.;
  
  TVector3 vglob;
  Double_t ytCenterModule=0.0, ytCenterCell=0.0;

  fCentersOfCellsPhiDir.Set(fNPhi*fNPHIdiv);
  fPhiCentersOfCells.Set(fNPhi*fNPHIdiv);

  Double_t R0 = GetIPDistance() + GetLongModuleSize()/2.;
  for(Int_t it=0; it<fNPhi; it++) { // cycle on modules
    ytCenterModule = -fParSM[1] + fPhiModuleSize*(2*it+1)/2;  // center of module
    for(Int_t ic=0; ic<fNPHIdiv; ic++) { // cycle on cells in module
      if(fNPHIdiv==2) {
        ytCenterCell = ytCenterModule + fPhiTileSize *(2*ic-1)/2.;
      } else if(fNPHIdiv==3){
        ytCenterCell = ytCenterModule + fPhiTileSize *(ic-1);
      } else if(fNPHIdiv==1){
        ytCenterCell = ytCenterModule;
      }
      fCentersOfCellsPhiDir.AddAt(ytCenterCell,ind);
      // Define grid on phi direction
      // Grid is not the same for different eta bin;
      // Effect is small but is still here
      phi = TMath::ATan2(ytCenterCell, R0);
      fPhiCentersOfCells.AddAt(phi, ind);

      //      printf(" ind %2.2i : y %8.3f ", ind, fCentersOfCellsPhiDir.At(ind)); 
      ind++;
    }
  }

  fCentersOfCellsEtaDir.Set(fNZ *fNETAdiv);
  fCentersOfCellsXDir.Set(fNZ *fNETAdiv);
  fEtaCentersOfCells.Set(fNZ *fNETAdiv * fNPhi*fNPHIdiv);
  if(fDebug>1) AliInfo(Form(" Cells grid in eta directions : size %i\n", fCentersOfCellsEtaDir.GetSize()));
  for(Int_t it=0; it<fNZ; it++) {
    AliJetDummyShishKebabTrd1Module *trd1 = GetShishKebabModule(it);
    nModule = fNPhi*it;
    for(Int_t ic=0; ic<fNETAdiv; ic++) {
      if(fNPHIdiv==2) {
        trd1->GetCenterOfCellInLocalCoordinateofSM(ic, xr, zr);      // case of 2X2
        GetCellPhiEtaIndexInSModule(0, nModule, 0, ic, iphiTemp, ieta); 
      } if(fNPHIdiv==3) {
        trd1->GetCenterOfCellInLocalCoordinateofSM3X3(ic, xr, zr);  // case of 3X3
        GetCellPhiEtaIndexInSModule(0, nModule, 0, ic, iphiTemp, ieta); 
      } if(fNPHIdiv==1) {
        trd1->GetCenterOfCellInLocalCoordinateofSM1X1(xr, zr);      // case of 1X1
        GetCellPhiEtaIndexInSModule(0, nModule, 0, ic, iphiTemp, ieta); 
      }
      fCentersOfCellsXDir.AddAt(Float_t(xr) - fParSM[0],ieta);
      fCentersOfCellsEtaDir.AddAt(Float_t(zr) - fParSM[2],ieta);
      // Define grid on eta direction for each bin in phi
      for(Int_t iphi=0; iphi<fCentersOfCellsPhiDir.GetSize(); iphi++) {
        x = xr + trd1->GetRadius();
        y = fCentersOfCellsPhiDir[iphi];
        r = TMath::Sqrt(x*x + y*y + zr*zr);
        theta = TMath::ACos(zr/r);
        eta   = AliJetDummyShishKebabTrd1Module::ThetaToEta(theta);
	//        ind   = ieta*fCentersOfCellsPhiDir.GetSize() + iphi;
        ind   = iphi*fCentersOfCellsEtaDir.GetSize() + ieta;
        fEtaCentersOfCells.AddAt(eta, ind);
      }
      //printf(" ieta %i : xr + trd1->GetRadius() %f : zr %f : eta %f \n", ieta, xr + trd1->GetRadius(), zr, eta);
    }
  }

  if(fDebug>1){
    for(Int_t i=0; i<fCentersOfCellsEtaDir.GetSize(); i++) {
      AliInfo(Form(" ind %2.2i : z %8.3f : x %8.3f", i+1, 
		   fCentersOfCellsEtaDir.At(i),fCentersOfCellsXDir.At(i)));
    }
  }
}

//------------------------------------------------------------------------------------
AliJetDummyShishKebabTrd1Module* AliJetDummyGeo::GetShishKebabModule(Int_t neta)
{
  // This method was too long to be included in the header file - the
  // rule checker complained about it's length, so we move it here. It
  // returns the shishkebabmodule at a given eta index point.

  static AliJetDummyShishKebabTrd1Module* trd1=0;
  if(fShishKebabTrd1Modules && neta>=0 && neta<fShishKebabTrd1Modules->GetSize()) {
    trd1 = (AliJetDummyShishKebabTrd1Module*)fShishKebabTrd1Modules->At(neta);
  } else trd1 = 0;
  return trd1;
}

//------------------------------------------------------------------------------------
void  AliJetDummyGeo::GetTransformationForSM()
{
  // Uses the geometry manager to load the transformation matrix
  // for the supermodules
  // Unused after 19 Jan, 2007 - keep for compatibility; 

  return;
  static Bool_t transInit=kFALSE;
  if(transInit) return;

  Int_t i=0;
  if(gGeoManager == 0) {
    Info("CreateTransformationForSM() "," Load geometry : TGeoManager::Import()");
    assert(0);
  }

  TGeoNode *tn = gGeoManager->GetTopNode();
  TGeoNode *node=0, *xen1 = 0;
  for(i=0; i<tn->GetNdaughters(); i++) {
    node = tn->GetDaughter(i);
    TString ns(node->GetName());
    if(ns.Contains(GetNameOfEMCALEnvelope())) {
      xen1 = node;
      break;
    }
  }

  if(!xen1) {
    Info("CreateTransformationForSM() "," geometry has not EMCAL envelope with name %s", 
    GetNameOfEMCALEnvelope());
    assert(0);
  }
  AliInfo(Form(" i %i : EMCAL Envelope is %s : #SM %i \n", i, xen1->GetName(), xen1->GetNdaughters()));
  for(i=0; i<xen1->GetNdaughters(); i++) {
    TGeoNodeMatrix *sm = (TGeoNodeMatrix*)xen1->GetDaughter(i);
    fMatrixOfSM[i] = sm->GetMatrix();
  }
  AliInfo(Form("transInit %d: ", transInit));
  transInit = kTRUE;
}

