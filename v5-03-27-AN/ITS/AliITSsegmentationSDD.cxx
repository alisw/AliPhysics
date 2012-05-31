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
#include <Riostream.h>
#include <TMath.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include "AliITSsegmentationSDD.h"
#include "AliITSDriftSpeedSDD.h"

/////////////////////////////////////////////////////////////////////////////
// Segmentation class for drift detectors                                  //
//                                                                         //
//     microcables                                  microcables            //
//       /\                                           /\                   //
//       ||                                           ||                   //
//       ||                                           ||                   //
//       ||                                           ||                   //
//        0                     256                    0                   //
//      0 |----------------------|---------------------| 511               //
//        |       time-bins      |     time-bins       |                   //
//        | a                    |                   a |                   //
//        | n                    |                   n |                   //
//    X <_|_o____________________|___________________o_|__                 //
//        | d                    |                   d |                   //
//        | e     LEFT SIDE      |   RIGHT SIDE      e |                   //
//        | s     CHANNEL 0      |   CHANNEL 1       s |                   //
//        |       Xlocal > 0     |   Xlocal < 0        |                   //
//    255 |----------------------|---------------------| 256               //
//                               |                                         //
//                               |                                         //
//                               V                                         //
//                               Z                                         //
/////////////////////////////////////////////////////////////////////////////

/* $Id$ */

const Float_t AliITSsegmentationSDD::fgkDxDefault = 35085.;
const Float_t AliITSsegmentationSDD::fgkDzDefault = 75264.;
const Float_t AliITSsegmentationSDD::fgkDyDefault = 300.;
const Float_t AliITSsegmentationSDD::fgkPitchDefault = 294.;
const Float_t AliITSsegmentationSDD::fgkClockDefault = 40.;
const Int_t AliITSsegmentationSDD::fgkHalfNanodesDefault = 256; 
const Int_t AliITSsegmentationSDD::fgkNsamplesDefault = 256;
const Int_t AliITSsegmentationSDD::fgkNchipsPerHybrid = 4;
const Int_t AliITSsegmentationSDD::fgkNanodesPerChip = 64;
const Float_t AliITSsegmentationSDD::fgkCm2Micron = 10000.;
const Float_t AliITSsegmentationSDD::fgkMicron2Cm = 1.0E-04;
ClassImp(AliITSsegmentationSDD)

//______________________________________________________________________
AliITSsegmentationSDD::AliITSsegmentationSDD(Option_t *opt) : AliITSsegmentation(),
fNsamples(0),
fNanodes(0),
fPitch(0),
fTimeStep(0),
fDriftSpeed(0),
fSetDriftSpeed(0){
  // Default constructor
  Init();
  if(strstr(opt,"TGeo")){
    if(!gGeoManager){
      AliError("Geometry is not initialized\n");
      return;
    }
    TGeoVolume *v=NULL;
    v = gGeoManager->GetVolume("ITSsddSensitivL3");
    if(!v){
      AliWarning("TGeo volume ITSsddSensitivL3 not found (hint: use v11Hybrid geometry)\n Using hardwired default values"); 
    }
    else {
      TGeoBBox *s=(TGeoBBox*)v->GetShape();
      SetDetSize(s->GetDX()*10000.,s->GetDZ()*20000.,s->GetDY()*20000.);
    }
  }
}

//______________________________________________________________________
void AliITSsegmentationSDD::Copy(TObject &obj) const {
  // protected method. copy this to obj
  AliITSsegmentation::Copy(obj);
  ((AliITSsegmentationSDD& ) obj).fNsamples = fNsamples;
  ((AliITSsegmentationSDD& ) obj).fNanodes = fNanodes;
  ((AliITSsegmentationSDD& ) obj).fPitch = fPitch;
  ((AliITSsegmentationSDD& ) obj).fTimeStep = fTimeStep;
  ((AliITSsegmentationSDD& ) obj).fDriftSpeed = fDriftSpeed;
  ((AliITSsegmentationSDD& ) obj).fSetDriftSpeed = fSetDriftSpeed;
}

//______________________________________________________________________
AliITSsegmentationSDD& AliITSsegmentationSDD::operator=(const AliITSsegmentationSDD &source){
   // = operator
   if(this==&source) return *this;
   source.Copy(*this);
   return *this;
}

//____________________________________________________________________________
AliITSsegmentationSDD::AliITSsegmentationSDD(const AliITSsegmentationSDD &source) :
    AliITSsegmentation(source),
fNsamples(0),
fNanodes(0),
fPitch(0),
fTimeStep(0),
fDriftSpeed(0),
fSetDriftSpeed(0){
  // copy constructor
  source.Copy(*this);
}

//----------------------------------------------------------------------
void AliITSsegmentationSDD::Init(){
// Standard initilisation routine
   fDriftSpeed=AliITSDriftSpeedSDD::DefaultDriftSpeed();
   fCorr=0;
   SetDetSize(fgkDxDefault,fgkDzDefault,fgkDyDefault);
   SetPadSize(fgkPitchDefault,fgkClockDefault);
   SetNPads(fgkHalfNanodesDefault,fgkNsamplesDefault);
}

//----------------------------------------------------------------------
void AliITSsegmentationSDD::
Neighbours(Int_t iX, Int_t iZ, Int_t* Nlist, Int_t Xlist[8], Int_t Zlist[8]) const {
  // returns neighbours for use in Cluster Finder routines and the like

    if(iX >= fNanodes) printf("iX > fNanodes %d %d\n",iX,fNanodes);
    if(iZ >= fNsamples) printf("iZ > fNsamples %d %d\n",iZ,fNsamples);
    *Nlist=4;
    Xlist[0]=Xlist[1]=iX;
    if(iX && (iX != fNanodes/2)) Xlist[2]=iX-1;
    else Xlist[2]=iX;
    if ((iX !=fNanodes/2 -1) && (iX != fNanodes)) Xlist[3]=iX+1;
    else Xlist[3]=iX;
    if(iZ) Zlist[0]=iZ-1;
    else Zlist[0]=iZ;
    if (iZ < fNsamples) Zlist[1]=iZ+1;
    else Zlist[1]=iZ;
    Zlist[2]=Zlist[3]=iZ;
}
//----------------------------------------------------------------------
Float_t AliITSsegmentationSDD::GetAnodeFromLocal(Float_t xloc,Float_t zloc) const {
  // returns anode coordinate (as float) starting from local coordinates
  Float_t xAnode=zloc*fgkCm2Micron/fPitch;
  if(xloc>0){   // left side (anodes 0-255, anode 0 at zloc<0)
    xAnode+=(Float_t)fNanodes/4;
  }else{ // right side (anodes 256-511, anode 0 at zloc>0)
    xAnode=3*fNanodes/4-xAnode;
  }
  return xAnode;
}

//----------------------------------------------------------------------
Float_t AliITSsegmentationSDD::GetLocalZFromAnode(Int_t nAnode) const{
  // returns local Z coordinate from anode number (integer)
  Float_t zAnode=(Float_t)nAnode+0.5;
  return GetLocalZFromAnode(zAnode);
}

//----------------------------------------------------------------------
Float_t AliITSsegmentationSDD::GetLocalZFromAnode(Float_t zAnode) const{
  // returns local Z coordinate from anode number (float)
  Float_t zloc=0.;
  if(zAnode<fNanodes/2){ // left side
    zloc=(zAnode*fPitch-fDz/2)*fgkMicron2Cm;
  }else{  // right side
    zAnode-=fNanodes/2;
    zloc=-(zAnode*fPitch-fDz/2)*fgkMicron2Cm;
  }
  return zloc;
}
//----------------------------------------------------------------------
Int_t AliITSsegmentationSDD::GetChipFromChannel(Int_t ix, Int_t iz) const {
  // returns chip number (in range 0-7) starting from channel number
  if(iz>=fNanodes  || iz<0 || ix>fNsamples){
    AliError("Bad cell number");
    return -1;
  }
  Int_t theChip=iz/fgkNanodesPerChip;
  return theChip;
}
//----------------------------------------------------------------------
Int_t AliITSsegmentationSDD::GetChipFromLocal(Float_t xloc, Float_t zloc) const {  
  // returns chip number (in range 0-7) starting from local coordinates
  Float_t detsize=fDz*fgkMicron2Cm;
  Float_t chipsize=detsize/(Float_t)fgkNchipsPerHybrid;
  zloc+=detsize/2.;
  if(zloc<-0.01 || zloc>detsize+0.01){ // 100 micron tolerance around edges
    AliError("Z local value out of sensitive SDD area");
    return -1;
  }
  Int_t iChip=int(zloc/chipsize);
  if(zloc<0.) iChip=0;  
  if(zloc>=detsize) iChip=fgkNchipsPerHybrid-1;
  if(iChip>=fgkNchipsPerHybrid || iChip<0){ 
    AliError(Form("Bad chip number %d",iChip));
    return -1;
  }
  Int_t iSide=GetSideFromLocalX(xloc);
  if(iSide==1) iChip=fgkNchipsPerHybrid-iChip+3;   // i.e. 7-iChip
  return iChip;
}
//----------------------------------------------------------------------
Int_t AliITSsegmentationSDD::GetChipsInLocalWindow(Int_t* array, Float_t zmin, Float_t zmax, Float_t xmin, Float_t xmax) const {
  // returns the numbers of the chips that read channels in a given region
  // of the module defined in local coordinates by zmin-zmax, xmin-max

  Int_t nChipInW = 0;
  Float_t zminDet=-fDz*fgkMicron2Cm/2.;
  Float_t zmaxDet=fDz*fgkMicron2Cm/2.;
  if(zmin<zminDet) zmin=zminDet;
  if(zmax>zmaxDet) zmax=zmaxDet;
  Float_t xminDet=-fDx*fgkMicron2Cm;
  Float_t xmaxDet=fDx*fgkMicron2Cm;
  if(xmin<xminDet) xmin=xminDet;
  if(xmax>xmaxDet) xmax=xmaxDet;
  Int_t n1=GetChipFromLocal(xmin,zmin);
  array[nChipInW]=n1;
  nChipInW++;
  Int_t n2=GetChipFromLocal(xmin,zmax);
  if(n2!=n1){
    Int_t imin=TMath::Min(n1,n2);
    Int_t imax=TMath::Max(n1,n2);
    for(Int_t ichip=imin; ichip<=imax; ichip++){
      if(ichip==n1) continue;
      array[nChipInW]=ichip;
      nChipInW++;    
    }
  }
  Int_t n3=GetChipFromLocal(xmax,zmin);
  if(n3!=n1){
    array[nChipInW]=n3;
    nChipInW++;
    Int_t n4=GetChipFromLocal(xmax,zmax);
    if(n4!=n3){
      Int_t imin=TMath::Min(n3,n4);
      Int_t imax=TMath::Max(n3,n4);
      for(Int_t ichip=imin; ichip<=imax; ichip++){
	if(ichip==n3) continue;
	array[nChipInW]=ichip;
	nChipInW++;    
      }
    }
  }
  return nChipInW;
}
//----------------------------------------------------------------------
void AliITSsegmentationSDD::GetPadIxz(Float_t x,Float_t z,
				      Int_t &timebin,Int_t &anode) const {
// Returns cell coordinates (time sample,anode)
// for given real local coordinates (x,z)

  // expects x, z in cm

  Float_t driftpath=fDx-TMath::Abs(x*fgkCm2Micron);
  timebin=(Int_t)(driftpath/fDriftSpeed/fTimeStep);
  anode=(Int_t)GetAnodeFromLocal(x,z);
  if(!fSetDriftSpeed){
    timebin=-999;
    AliWarning("Drift speed not set: timebin is dummy");
  }

}
//----------------------------------------------------------------------
void AliITSsegmentationSDD::GetPadCxz(Int_t timebin,Int_t anode,
				      Float_t &x ,Float_t &z) const{
  // Transform from cell to real local coordinates
  // returns x, z in cm

  // the +0.5 means that an # and time bin # should start from 0 !!! 
  // the +0.5 means that an # and time bin # should start from 0 !!! 

  Float_t driftpath=GetDriftTimeFromTb(timebin)*fDriftSpeed;
  if (anode < fNanodes/2){ // left side, positive x
    x=(fDx-driftpath)*fgkMicron2Cm;
  }else{ // right side, negative x
    x = -(fDx-driftpath)*fgkMicron2Cm;
  }
  z=GetLocalZFromAnode(anode);
  if(!fSetDriftSpeed){
    x=-9999.;
    AliWarning("Drift speed not set: x coord. is dummy");
  }
}
//----------------------------------------------------------------------
void AliITSsegmentationSDD::GetPadTxz(Float_t &x,Float_t &z) const{
    // Get anode and time bucket as floats - numbering from 0

    // expects x, z in cm

    Float_t xloc=x;
    Float_t zloc=z;
    Float_t driftpath=fDx-TMath::Abs(fgkCm2Micron*xloc);
    x=driftpath/fDriftSpeed/fTimeStep;
    if (xloc < 0) x = -x;
    z=GetAnodeFromLocal(xloc,zloc);
    if(!fSetDriftSpeed){
      x=-9999.;
      AliWarning("Drift speed not set: x coord. is dummy");
    }
}
//----------------------------------------------------------------------
void AliITSsegmentationSDD::Print(Option_t *opt) const {
  // Print SDD segmentation Parameters

   cout << "**************************************************" << endl;
   cout << "  Silicon Drift Detector Segmentation Parameters  " << opt << endl;
   cout << "**************************************************" << endl;
   cout << "Number of Time Samples: " << fNsamples << endl;
   cout << "Number of Anodes: " << fNanodes << endl;
   cout << "Time Step (ns): " << fTimeStep << endl;
   cout << "Anode Pitch (um): " << fPitch << endl;
   cout << "Full Detector Width     (x): " << fDx;
   cout<<" (Default is "<<fgkDxDefault<<") "<<endl;
   cout << "Half Detector Length    (z): " << fDz;
   cout<<" (Default is "<<fgkDzDefault<<") "<<endl;
   cout << "Full Detector Thickness (y): " << fDy;
   cout<<" (Default is "<<fgkDyDefault<<") "<<endl;
   cout << "**************************************************" << endl;

}

//----------------------------------------------------------------------
void AliITSsegmentationSDD::PrintDefaultParameters() const {
  // print SDD parameters defined as static const data members

  cout << "**************************************************" << endl;
  cout << "  Silicon Drift Detector Segmentation Parameters  " << endl;
  cout << " Default values defined as static const data members"<< endl;
  cout << " Actual values can be set with the relevant setters"<< endl; 
  cout << "**************************************************" << endl;
  cout<<"fgkDxDefault = "<<fgkDxDefault<<endl;
  cout<<"fgkDzDefault = "<<fgkDzDefault<<endl;
  cout<<"fgkDyDefault = "<<fgkDyDefault<<endl;
  cout<<"fgkPitchDefault = "<<fgkPitchDefault<<endl;
  cout<<"fgkClockDefault = "<<fgkClockDefault<<endl;
  cout<<"fgkHalfNanodesDefault = "<<fgkHalfNanodesDefault<<endl;
  cout<<"fgkNsamplesDefault = "<<fgkNsamplesDefault<<endl;
  cout << "**************************************************" << endl;

}

//______________________________________________________________________
Bool_t AliITSsegmentationSDD::LocalToDet(Float_t x,Float_t z,
                                         Int_t &ix,Int_t &iz) const {
// Transformation from Geant detector centered local coordinates (cm) to
// time bucket numbers ix and anode number iz.
// Input:
// Float_t   x      detector local coordinate x in cm with respect to the
//                  center of the sensitive volume.
// Float_t   z      detector local coordinate z in cm with respect to the
//                  center of the sensitive volulme.
// Output:
// Int_t    ix      detector x time coordinate. Has the range 0<=ix<fNsamples.
// Int_t    iz      detector z anode coordinate. Has the range 0<=iz<fNandoes.
//   A value of -1 for ix or iz indecates that this point is outside of the
// detector segmentation as defined.

    Float_t dx,dz;

    ix = -1; // default values
    iz = -1; // default values
    dx = -fgkMicron2Cm*Dx(); // lower left edge in cm.
    dz = -0.5*fgkMicron2Cm*Dz(); // lower left edge in cm.
    if(x<dx || x>-dx) {
      AliWarning(Form("Input argument %f out of range (%f, %f)",x,dx,-dx));
      return kFALSE; // outside of defined volume.
    }
    if(z<dz || z>-dz) {
      AliWarning(Form("Input argument %f out of range (%f, %f)",z,dz,-dz));
      return kFALSE; // outside of defined volume.
    }
    GetPadIxz(x,z,ix,iz);
    if(!fSetDriftSpeed){
      ix=-999;
      AliWarning("Drift speed not set: timebin is dummy");
      return kFALSE;
    }
    return kTRUE; // Found ix and iz, return.
}
//______________________________________________________________________
void AliITSsegmentationSDD::DetToLocal(Int_t ix,Int_t iz,Float_t &x,Float_t &z) const
{
// Transformation from Detector time bucket and anode coordiantes to Geant 
// detector centerd local coordinates (cm).
// Input:
// Int_t    ix      detector x time coordinate. Has the range 0<=ix<fNsamples.
// Int_t    iz      detector z anode coordinate. Has the range 0<=iz<fNandoes.
// Output:
// Float_t   x      detector local coordinate x in cm with respect to the
//                  center of the sensitive volume.
// Float_t   z      detector local coordinate z in cm with respect to the
//                  center of the sensitive volulme.
// If ix and or iz is outside of the segmentation range a value of -Dx()
// or -0.5*Dz() is returned.

  x=-Dx();
  z=-0.5*Dz();
  if(ix<0 || ix>=Npx()) {
    AliWarning(Form("Input argument %d out of range (0, %d)",ix,Npx()));
    return; // outside of detector
  }
  if(iz<0 || iz>=Npz()) {
    AliWarning(Form("Input argument %d out of range (0, %d)",iz,Npz()));
    return; // outside of detctor
  }
  GetPadCxz(ix,iz,x,z);
  if(!fSetDriftSpeed){
    x=-9999.;
    AliWarning("Drift speed not set: x coord. is dummy");
  }
  return; // Found x and z, return.
}
