/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$  */

//========================================================================
//
// This class is a helper, producing ITS aligmnent objects.
// It provides also some useful functions
// See the parameters of the misalignment at the end of this script.
//
// Main author: L. Gaudichet
// Contact: andrea.dainese@lnl.infn.it
//
//========================================================================

#include <TRandom3.h>
#include <TClonesArray.h>
#include <TClass.h>


#include "AliLog.h"
#include "AliAlignObjParams.h"
#include "AliITSMisalignMaker.h"


ClassImp(AliITSMisalignMaker)
  
const Int_t AliITSMisalignMaker::fgkNLadders[AliITSMisalignMaker::kNLayers] = {20,40,14,22,34,38};
const Int_t AliITSMisalignMaker::fgkNDetectors[AliITSMisalignMaker::kNLayers] = {4,4,6,8,22,25};


//________________________________________________________________________
AliITSMisalignMaker::AliITSMisalignMaker():
  fRnd(),
  fInd(0),
  fAlobj(TClonesArray("AliAlignObjParams",4000)),
  fStrSPD("ITS/SPD"),
  fStrSDD("ITS/SDD"),
  fStrSSD("ITS/SSD"),
  fStrStave("/Stave"),
  fStrHalfStave("/HalfStave"),
  fStrLadder("/Ladder"),
  fStrSector("/Sector"),
  fStrSensor("/Sensor")
{
  //
  // defaul constructor
  //
  fRnd.SetSeed(38217945);
}
//________________________________________________________________________
Int_t AliITSMisalignMaker::AddAlignObj(char* name,Double_t dx,Double_t dy,Double_t dz,
				Double_t dpsi,Double_t dtheta,Double_t dphi,const char* distrib) {
  //
  // misalignment by symname
  //
  Double_t vx=0.,vy=0.,vz=0.,vpsi=0.,vtheta=0.,vphi=0.;

  TString sdistrib(distrib);

  if(sdistrib==TString("gaussian")) {
    vx = GaussCut(0,dx/3.,dx); // mean, sigma, max absolute value 
    vy = GaussCut(0,dy/3.,dy);
    vz = GaussCut(0,dz/3.,dz);
    vpsi   = GaussCut(0,dpsi/3.,  dpsi );
    vtheta = GaussCut(0,dtheta/3.,dtheta);
    vphi   = GaussCut(0,dphi/3.,  dphi);
  }else if(sdistrib==TString("uniform")){ 
    vx = fRnd.Uniform(-dx,dx);
    vy = fRnd.Uniform(-dy,dy);
    vz = fRnd.Uniform(-dz,dz);
    vpsi = fRnd.Uniform(-dpsi,dpsi);
    vtheta = fRnd.Uniform(-dtheta,dtheta);
    vphi = fRnd.Uniform(-dphi,dphi);
  }else if(sdistrib==TString("fixed")){
    vx=dx;
    vy=dy;
    vz=dz;
    vpsi=dpsi;
    vtheta=dtheta;
    vphi=dphi;
  }else{
    AliFatal(Form("Invalid string \"%s\" specifying the misalignment type for the volume \"%s\""));
  }

  new(fAlobj[fInd]) AliAlignObjParams(name,0,vx,vy,vz,vpsi,vtheta,vphi,kFALSE);

  AliAlignObjParams* itsalobj = (AliAlignObjParams*) fAlobj.UncheckedAt(fInd);
  itsalobj->ApplyToGeometry();

  fInd++;

  return kTRUE;
}


//________________________________________________________________________
Int_t AliITSMisalignMaker::AddAlignObj(Int_t lay,Double_t dx,Double_t dy,Double_t dz,
				 Double_t dpsi,Double_t dtheta,Double_t dphi,Bool_t unif) {
  //
  // misalignment at the level of ladders/modules
  //
  lay+=1; // layers are numbered from 1 to 6 in AliGeomManager

  printf("LAYER %d  MODULES %d\n",lay,AliGeomManager::LayerSize(lay));

  for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(lay); iModule++) {

    Double_t vx,vy,vz,vpsi,vtheta,vphi;
    
    if(!unif) {
      vx = GaussCut(0,dx/3.,dx); // mean, sigma, max absolute value 
      vy = GaussCut(0,dy/3.,dy);
      vz = GaussCut(0,dz/3.,dz);
      vpsi   = GaussCut(0,dpsi/3.,  dpsi );
      vtheta = GaussCut(0,dtheta/3.,dtheta);
      vphi   = GaussCut(0,dphi/3.,  dphi);
    } else {
      vx = fRnd.Uniform(-dx,dx);
      vy = fRnd.Uniform(-dy,dy);
      vz = fRnd.Uniform(-dz,dz);
      vpsi = fRnd.Uniform(-dpsi,dpsi);
      vtheta = fRnd.Uniform(-dtheta,dtheta);
      vphi = fRnd.Uniform(-dphi,dphi);
    }
    
    UShort_t volid = AliGeomManager::LayerToVolUID(lay,iModule);
    const char *symname = AliGeomManager::SymName(volid);
    
    new(fAlobj[fInd]) AliAlignObjParams(symname,volid,vx,vy,vz,vpsi,vtheta,vphi,kFALSE);
    AliAlignObjParams* itsalobj = (AliAlignObjParams*) fAlobj.UncheckedAt(fInd);
    itsalobj->ApplyToGeometry();
    fInd++; 
  }

  return kTRUE;
}

//________________________________________________________________________
Int_t AliITSMisalignMaker::AddAlignObj(Int_t lay,Int_t ladd,Double_t dx,Double_t dy,Double_t dz,
				       Double_t dpsi,Double_t dtheta,Double_t dphi,
				       Double_t xShift,Double_t yShift,Double_t zShift,
				       Double_t psiShift,Double_t thetaShift,Double_t phiShift,
				       Bool_t unif) {
  //
  // misalignment at the level of half-staves/ladders (ladd=-1 means that all ladders are scanned)
  //
  Double_t vx,vy,vz,vpsi,vtheta,vphi;
  Double_t tr[3],rot[3];  
  
  Int_t laddMin = ladd;
  Int_t laddMax = laddMin+1;
  if (ladd<0) {
    laddMin = 0;
    laddMax = fgkNLadders[lay];
  }

  for (Int_t iLadd=laddMin; iLadd<laddMax; iLadd++) {

    Int_t nHS = 1; 
    if (lay<2) nHS = 2;
    for (Int_t iHalfStave=0; iHalfStave<nHS; iHalfStave++) {
      
      if(!unif) {
	vx = GaussCut(0,dx/3.,dx); // mean, sigma, max absolute value 
	vy = GaussCut(0,dy/3.,dy);
	vz = GaussCut(0,dz/3.,dz);
	vpsi   = GaussCut(0,dpsi/3.,  dpsi );
	vtheta = GaussCut(0,dtheta/3.,dtheta);
	vphi   = GaussCut(0,dphi/3.,  dphi);
      } else {
	vx = fRnd.Uniform(-dx,dx);
	vy = fRnd.Uniform(-dy,dy);
	vz = fRnd.Uniform(-dz,dz);
	vpsi = fRnd.Uniform(-dpsi,dpsi);
	vtheta = fRnd.Uniform(-dtheta,dtheta);
	vphi = fRnd.Uniform(-dphi,dphi);
      }

      TString name(GetHalfStaveLadderSymbName(lay,iLadd,iHalfStave));

      // first apply half-stave / ladder level misalignment
      AliAlignObjParams aaop(name.Data(),0,vx,vy,vz,vpsi,vtheta,vphi,kFALSE); // set them as local
      aaop.GetPars(tr,rot); // global

      // then, apply layer-level misalignment (only for SDD and SSD)
      if(lay>1) {
	tr[0] += xShift;
	tr[1] += yShift;
	tr[2] += zShift;
	rot[0] += psiShift;
	rot[1] += thetaShift;
	rot[2] += phiShift;
      }
      new(fAlobj[fInd]) AliAlignObjParams(name.Data(),0,tr[0],tr[1],tr[2],rot[0],rot[1],rot[2],kTRUE); // set them as global

      AliAlignObjParams* itsalobj = (AliAlignObjParams*) fAlobj.UncheckedAt(fInd);
      itsalobj->ApplyToGeometry();
      fInd++;
    }
  }
  
  return kTRUE;
}


//________________________________________________________________________
Int_t AliITSMisalignMaker::AddSectorAlignObj(Int_t sectMin,Int_t sectMax,
				       Double_t dx,Double_t dy,Double_t dz,
				       Double_t dpsi,Double_t dtheta,Double_t dphi,
				       Double_t xShift,Double_t yShift,Double_t zShift,
				       Double_t psiShift,Double_t thetaShift,Double_t phiShift,Bool_t unif) {
  //
  // misalignment at the level of SPD sectors and half-barrels
  // 

  if ((sectMin<1) || (sectMax>10)) return kFALSE;
  Double_t vx,vy,vz,vpsi,vtheta,vphi;
  Double_t tr[3],rot[3];  

  for (Int_t iSect = sectMin-1; iSect<sectMax; iSect++) {

    // first, apply sector level misalignment    
    if(!unif) {
      vx = GaussCut(0,dx/3.,dx); // mean, sigma, max absolute value 
      vy = GaussCut(0,dy/3.,dy);
      vz = GaussCut(0,dz/3.,dz);
      vpsi   = GaussCut(0,dpsi/3.,  dpsi );
      vtheta = GaussCut(0,dtheta/3.,dtheta);
      vphi   = GaussCut(0,dphi/3.,  dphi);
    } else {
      vx = fRnd.Uniform(-dx,dx);
      vy = fRnd.Uniform(-dy,dy);
      vz = fRnd.Uniform(-dz,dz);
      vpsi = fRnd.Uniform(-dpsi,dpsi);
      vtheta = fRnd.Uniform(-dtheta,dtheta);
      vphi = fRnd.Uniform(-dphi,dphi);
    }

    TString name(GetSymbName(0));
    name += fStrSector;
    name += iSect;


    AliAlignObjParams aaop(name.Data(),0,vx,vy,vz,vpsi,vtheta,vphi,kFALSE); // set them as local
    aaop.GetPars(tr,rot); // global

    // then, apply half-barrel level misalignment
    tr[0] += xShift;
    tr[1] += yShift;
    tr[2] += zShift;
    rot[0] += psiShift;
    rot[1] += thetaShift;
    rot[2] += phiShift;

    new(fAlobj[fInd]) AliAlignObjParams(name.Data(),0,tr[0],tr[1],tr[2],rot[0],rot[1],rot[2],kTRUE); // set them as global

    AliAlignObjParams* itsalobj = (AliAlignObjParams*) fAlobj.UncheckedAt(fInd);
    itsalobj->ApplyToGeometry();
    fInd++;
  }
  return kTRUE;
}

//________________________________________________________________________
Double_t AliITSMisalignMaker::GaussCut(Double_t mean,Double_t sigma,Double_t absMax) {
  //
  // random from gaussian with cut on tails
  //
  Double_t val = fRnd.Gaus(mean,sigma);
  while (TMath::Abs(val-mean)>absMax)
    val = fRnd.Gaus(mean,sigma);
  return val;
}

//________________________________________________________________________
const char* AliITSMisalignMaker::GetSymbName(Int_t layer) const {
  //
  // be careful : SPD0 and SPD1 are not physically separated 
  //
  TString name;
  switch (layer) {
  case 0:
  case 1: name = fStrSPD; name += layer; break;
  case 2:
  case 3: name = fStrSDD; name += layer; break;
  case 4:
  case 5: name = fStrSSD; name += layer; break;
  default: AliFatal("Wrong layer index");
  }
  return name.Data();
}

//________________________________________________________________________
const char* AliITSMisalignMaker::GetSymbName(Int_t layer, Int_t ladder, Int_t det) const {
  //
  // symname from layer, ladder, detector
  //
  TString symname(GetHalfStaveLadderSymbName(layer,ladder,det));
  if(layer<=2){
    symname+="Ladder";
  }else if(layer<=6){
    symname+="Sensor";
  }else{
    AliError("Invalid layer!");
    return 0;
  }
  symname+=det;
  return symname.Data();
}

//________________________________________________________________________
const char* AliITSMisalignMaker::GetSymbName(Int_t layer,Int_t ladd) const {
  //
  // Get logical names at the level of staves / ladders
  //
  TString name(GetSymbName(layer));
  if (layer==0) { // SPD1

    int sector = ladd/2;
    name += fStrSector;
    name += sector;
    int stave = ladd-sector*2;
    name += fStrStave;
    name += stave;
  }
  else if (layer==1) { // SPD2

    int sector = ladd/4;
    name += fStrSector;
    name += sector;
    int stave = ladd-sector*4;
    name += fStrStave;
    name += stave;
  }
  else if (layer>=2 && layer<=5) { // SDD and SSD
    name += fStrLadder;
    name += ladd;
  }
  else {
    AliFatal("Wrong layer index");
  }
  return name.Data();
}

//________________________________________________________________________
const char* AliITSMisalignMaker::GetHalfStaveLadderSymbName(Int_t layer,Int_t ladd,Int_t halfStave) const {
  //
  // Get logical names at the level of half-staves (SPD) or ladders (SDD and SSD)
  //
  TString name(GetSymbName(layer));
  if (layer==0) { // SPD1

    int sector = ladd/2;
    name += fStrSector;
    name += sector;
    int stave = ladd-sector*2;
    name += fStrStave;
    name += stave;
    name += fStrHalfStave;
    name += halfStave;
  }
  else if (layer==1) { // SPD2

    int sector = ladd/4;
    name += fStrSector;
    name += sector;
    int stave = ladd-sector*4;
    name += fStrStave;
    name += stave;
    name += fStrHalfStave;
    name += halfStave;
  }
  else if (layer>=2 && layer<=5) { // SDD and SSD
    name += fStrLadder;
    name += ladd;
  } 
  else {
    AliFatal("Wrong layer index");
  }
  return name.Data();
}

//________________________________________________________________________
const char* AliITSMisalignMaker::GetParentSymName(const char* symname) {
  //
  // symnane of parent volume
  //
  TString parent(symname);
  // Give the symname of 
  if(parent.BeginsWith('/')) parent.Remove(TString::kLeading,'/');
  if(parent.EndsWith("/")) parent.Remove(TString::kTrailing,'/');
  
  if(!parent.CountChar('/')) AliErrorClass("Not a valid symbolic name");

  Int_t layer,level;
  GetLayerAndLevel(symname,layer,level);
  if(level==1) return "ITS";
  
  parent.Remove(parent.Last('/'));
  
  if((layer==0 || layer==1) && level==2){
    parent.Remove(parent.Last('/'));
    parent[7]='0';
  }
    
  return parent.Data(); 
}

//________________________________________________________________________
Bool_t AliITSMisalignMaker::GetLayerAndLevel(const char* symname, Int_t &layer, Int_t &level) {
  //
  // given the symbolic name set layer and level
  //
  const char* basename[6] = {"ITS/SPD0/Sector","ITS/SPD1/Sector","ITS/SDD2/Ladder","ITS/SDD3/Ladder","ITS/SSD4/Ladder","ITS/SSD5/Ladder"};
  TString strSym(symname);
  if(strSym=="ITS"){
    level=0;
    layer=-1;
    return kTRUE;
  }
  Int_t i;
  for(i=0; i<6; i++){
    if(strSym.BeginsWith(basename[i])) break;
  }

  if(i>=6){
    AliErrorClass(Form("%s is not a valid symbolic name for an ITS alignable volume",strSym.Data()));
    return kFALSE;
  }
  
  layer=i;
  //The part above could be replaced by just
  // TString seventh = strSym[7];
  // layer = seventh.Atoi();
  // if we don't need to check the validity of the symname
  
  level=1;
  switch(layer){
    case 0:
    case 1:
      if(strSym.Contains("Stave")) level=2;
      if(strSym.Contains("Ladder")) level=3;
      break;
    case 2:
    case 3:
    case 4:
    case 5:
      if(strSym.Contains("Sensor")) level=2;
  }
  
  return kTRUE;
}

//________________________________________________________________________
Int_t AliITSMisalignMaker::GetNSisters(const char* symname) {
  //
  // number of volumes on same level
  //
  Int_t layer,level;
  if(!GetLayerAndLevel(symname,layer,level)) return -1;
  if(level==0) return -1;
  if(level==1) return GetNLadders(layer);
  if(level==2) return GetNDetectors(layer);
  AliErrorClass(Form("Invalid layer and level"));
  return -1;
}

//________________________________________________________________________
Int_t AliITSMisalignMaker::GetNDaughters(const char* symname) {
  //
  // number of daughter volumes
  // 
  Int_t layer,level;
  if(!GetLayerAndLevel(symname,layer,level)) return -1;
  if(level==0) {
    Int_t nLadders = 0;
    for(Int_t lay=0; lay<6; lay++) nLadders += GetNLadders(lay);
    return nLadders;
  }
  if(level==1) return GetNDetectors(layer);
  if(level==2){
    AliWarningClass(Form("Volume %s is a sensitive volume and has no alignable dauthers",symname));
    return -1;
  }
  AliErrorClass(Form("Invalid layer and level"));
  return -1;
}

/*
//________________________________________________________________________
TString AliITSMisalignMaker::GetSymbName(Int_t layer,Int_t ladd,Int_t mod) const {

  // Get logical names at the level of SPD ladders / SDD and SSD modules

  Int_t halfStave = mod/2;
  TString name = GetHalfStaveLadderSymbName(layer,ladd,halfStave);

  if (layer<2) { // SPD
    name += fStrLadder;
    name += mod;
  } 
  else if (layer>=2 && layer<=5) { // SDD and SSD
    name += fStrSensor;
    name += mod;
  }
  else {
    AliFatal("Wrong layer index");
  }
  return name;
}
*/

