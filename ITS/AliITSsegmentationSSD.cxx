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

#include <Riostream.h>
#include <TMath.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include "AliITSsegmentationSSD.h"

//////////////////////////////////////////////////////
// Segmentation class for                           //
// silicon strips                                   //
//                                                  //
//////////////////////////////////////////////////////
const Float_t AliITSsegmentationSSD::fgkDxDefault = 73000.;
const Float_t AliITSsegmentationSSD::fgkDzDefault = 40000.;
const Float_t AliITSsegmentationSSD::fgkDyDefault = 300.;
const Float_t AliITSsegmentationSSD::fgkPitchDefault = 95.;
const Int_t AliITSsegmentationSSD::fgkNstripsDefault = 768;
const Int_t AliITSsegmentationSSD::fgkNchipsPerSide = 6;
const Int_t AliITSsegmentationSSD::fgkNstripsPerChip = 128;

ClassImp(AliITSsegmentationSSD)

AliITSsegmentationSSD::AliITSsegmentationSSD(Option_t *opt): AliITSsegmentation(),
fNstrips(0),
fStereoP(0),
fStereoN(0),
fPitch(0),
fLayer(0){
    // default constructor
  SetDetSize(fgkDxDefault,fgkDzDefault,fgkDyDefault);
  SetPadSize(fgkPitchDefault,0.);
  SetNPads(fgkNstripsDefault,0);
  SetAngles();
  if(strstr(opt,"TGeo")){
    if(!gGeoManager){
      AliError("Geometry is not initialized\n");
      return;
    }
    TGeoVolume *v=NULL;
    v = gGeoManager->GetVolume("ITSssdSensitivL5");
    if(!v){
      AliWarning("TGeo volumeITSssdSensitivL5  not found (hint: use v11Hybrid geometry)\n Using hardwired default values"); 
    }
    else {
      TGeoBBox *s=(TGeoBBox*)v->GetShape();
      SetDetSize(s->GetDX()*20000.,s->GetDZ()*20000.,s->GetDY()*20000.);
    }
  }
}

//______________________________________________________________________
void AliITSsegmentationSSD::Copy(TObject &obj) const {
  // protected method. copy this to obj
  AliITSsegmentation::Copy(obj);
  ((AliITSsegmentationSSD& ) obj).Clear();
  ((AliITSsegmentationSSD& ) obj).fNstrips = fNstrips;
  ((AliITSsegmentationSSD& ) obj).fStereoP = fStereoP;
  ((AliITSsegmentationSSD& ) obj).fStereoN = fStereoN;
  ((AliITSsegmentationSSD& ) obj).fLayer   = fLayer;
  ((AliITSsegmentationSSD& ) obj).fPitch   = fPitch;
  ((AliITSsegmentationSSD& ) obj).fLayer   = fLayer;
 
}

//______________________________________________________________________
AliITSsegmentationSSD& AliITSsegmentationSSD::operator=(
                        const AliITSsegmentationSSD &source){
// Operator =
  if(this != &source){
    source.Copy(*this);
  }
  return *this;
}
//______________________________________________________________________
AliITSsegmentationSSD::AliITSsegmentationSSD(const AliITSsegmentationSSD &source):
    AliITSsegmentation(source),
fNstrips(0),
fStereoP(0),
fStereoN(0),
fPitch(0),
fLayer(0){
    // copy constructor
  source.Copy(*this);
}
//----------------------------------------------------------------------
void AliITSsegmentationSSD::Init(){
    // standard initalizer

    SetPadSize(fgkPitchDefault,0.);
    SetNPads(fgkNstripsDefault,0);
    SetAngles();
}
//----------------------------------------------------------------------
void AliITSsegmentationSSD::Angles(Float_t &aP,Float_t &aN) const{
  // P and N side stereo angles
  aP = fStereoP;
  aN = fStereoN;
}
//----------------------------------------------------------------------
void AliITSsegmentationSSD::SetLayer(Int_t l){
  //set fLayer data member (only 5 or 6 are allowed)
    if (l==5) fLayer =5;
    else if (l==6) fLayer =6;
    else AliError(Form("Layer can be 5 or 6, not %d",l));
}
//----------------------------------------------------------------------
void AliITSsegmentationSSD::GetPadTxz(Float_t &x,Float_t &z) const{
    // returns P and N sided strip numbers for a given location.
    // Transformation from microns detector center local coordinates
    // to detector P and N side strip numbers..
    /*                       _-  Z
                    + angle /    ^
        fNstrips           v     |   N-Side        ...0
            \-------/------------|-----------\--------\
            |\\\\\\/////////////.|\\\\\\\\\\\\\\\\\\\\|
            |0\\\\/////////////..|.\\\\\\\\\\\\\\\\\\\|
            |00\\/////////////...|..\\\\\\\\\\\\\\\\\\|
       X <--|000/////////////... |...\\\\\\\\\\\\\\\\\|
            |00/////////////...  | ...\\\\\\\\\\\\\\\\|
            |0/////////////...   |  ...\\\\\\\\\\\\\\\|
            |//////////////...   |  ...\\\\\\\\\\\\\\\|
            /-----\--------------|--------------------/
        fNstrips-1             P-Side              ...0
                     |0\
                     |00\
	Dead region: |000/
                     |00/
                     |0/
    // expects x, z in cm
    */

  /*
    Float_t stereoP, stereoN;
    Angles(stereoP,stereoN);
    Float_t tanP = TMath::Tan(stereoP);
    Float_t tanN = TMath::Tan(-stereoN);
    Float_t x1 = x;
    Float_t z1 = z;
    x1 += fDx/2;
    z1 += fDz/2;
    x   = (x1 - z1*tanP)/fPitch;
    z   = (x1 - tanN*(z1 - fDz))/fPitch;
  */

  Float_t P=0;
  Float_t N=0;
  if(fLayer==5) {
    P = 105.26*x - 0.7895*z + 382.000; //- 0.79*z + 381.89;
    N = P + 3.684*z - 4; 
  }
  else if(fLayer==6) {
    P = -105.26*x - 0.7895*z + 385.000; //- 0.79*z + 384.66;
    N = P + 3.684*z + 4;
  }
  else AliError("Layer can be 5 or 6");

  x=P;
  z=N;

}
//----------------------------------------------------------------------
void AliITSsegmentationSSD::GetPadIxz(Float_t x,Float_t z,Int_t &iP,Int_t &iN) const {
  // returns P and N sided strip numbers for a given location.
  // expects x, z in cm

  GetPadTxz(x,z);
  iP = Int_t(x+0.5);
  iN = Int_t(z+0.5);
}
//-------------------------------------------------------
void AliITSsegmentationSSD::GetPadCxz(Float_t iP,Float_t iN,Float_t &x,Float_t &z) const {
    // actually this is the GetCrossing(Float_t &,Float_t &)
    // returns local x, z  in cm 
  const Float_t kStartXzero=3.64325;
  const Float_t kDeltaXzero5or6=0.02239;
  const Float_t kDeltaZ5to6=7.6/7.0;

  z = 1.9*(iN-iP)/7.0;
  x = kStartXzero-(285*iN + 1045*iP)/140000.0;

  if (fLayer==5){
    z += kDeltaZ5to6;
    x = -x + kDeltaXzero5or6;
  }
  else if (fLayer==6) {
    z -= kDeltaZ5to6;
    x += kDeltaXzero5or6;
  }
  else {
    AliWarning("Layer shoudl be 5 or 6");
    x = -99999;
    z = -99999;
  }
}
//______________________________________________________________________
Bool_t AliITSsegmentationSSD::LocalToDet(Float_t x,Float_t z,
				       Int_t &iP,Int_t &iN) const {
    // Transformation from Geant cm detector center local coordinates
    // to detector P and N side strip numbers..
    /*                       _-  Z
                    + angle /    ^
        fNstrips           v     |   N-Side        ...0
            \-------/------------|-----------\--------\
            |\\\\\\/////////////.|\\\\\\\\\\\\\\\\\\\\|
            |0\\\\/////////////..|.\\\\\\\\\\\\\\\\\\\|
            |00\\/////////////...|..\\\\\\\\\\\\\\\\\\|
       X <--|000/////////////... |...\\\\\\\\\\\\\\\\\|
            |00/////////////...  | ...\\\\\\\\\\\\\\\\|
            |0/////////////...   |  ...\\\\\\\\\\\\\\\|
            |//////////////...   |  ...\\\\\\\\\\\\\\\|
            /-----\--------------|--------------------/
        fNstrips-1             P-Side              ...0
                     |0\
                     |00\
	Dead region: |000/
                     |00/
                     |0/
    */
  Float_t dx,dz;
  const Double_t kconst = 1.0E-04; // convert microns to cm.
  dx = 0.5*kconst*Dx();
  dz = 0.5*kconst*Dz();
  if( (x<-dx) || (x>dx) ) { 
    iP=-1; 
    AliWarning(Form("Input argument %f out of range (%f, %f)",x,dx,-dx));
    return kFALSE; // outside of defined volume.
  } // outside x range.
  if( (z<-dz) || (z>dz) ) { 
    iN=-1; 
    AliWarning(Form("Input argument %f out of range (%f, %f)",z,dz,-dz));
    return kFALSE; // outside of defined volume.
  }
  
  //x /= kconst;  // convert to microns
  //z /= kconst;  // convert to microns
  this->GetPadTxz(x,z);
  
  // first for P side
  iP = (Int_t) x;
  if(iP<0 || iP>=fNstrips) iP=-1; // strip number must be in range.
  // Now for N side)
  iN = (Int_t) z;
  if(iN<0 || iN>=fNstrips) iN=-1; // strip number must be in range.
  return kTRUE;
}
//----------------------------------------------------------------------
void AliITSsegmentationSSD::DetToLocal(Int_t ix,Int_t iPN,
				       Float_t &x,Float_t &z) const{
    // Transformation from detector segmentation/cell coordiantes starting
    // from 0. iPN=0 for P side and 1 for N side strip. Returned is z=0.0
    // and the corresponding x value..
    /*                       _-  Z
                    + angle /    ^
        fNstrips           v     |   N-Side        ...0
            \-------/------------|-----------\--------\
            |\\\\\\/////////////.|\\\\\\\\\\\\\\\\\\\\|
            |0\\\\/////////////..|.\\\\\\\\\\\\\\\\\\\|
            |00\\/////////////...|..\\\\\\\\\\\\\\\\\\|
       X <--|000/////////////... |...\\\\\\\\\\\\\\\\\|
            |00/////////////...  | ...\\\\\\\\\\\\\\\\|
            |0/////////////...   |  ...\\\\\\\\\\\\\\\|
            |//////////////...   |  ...\\\\\\\\\\\\\\\|
            /-----\--------------|--------------------/
        fNstrips-1             P-Side              ...0
                     |0\
                     |00\
	Dead region: |000/
                     |00/
                     |0/
    */
    // for strips p-side
    // x = a + b + z*tan(fStereoP); a = Dpx(iP)*(iP+0.5)-dx; b = dz*th;
    // for strips n-side
    // x = a + b + z*tan(fStereoP); a = Dpx(iN)*(iN+0.5)-dx; b = -dz*th;
    AliWarning("This function has not been verified. Should probably use GetPadCxz");
    const Double_t kconst = 1.0E-04; // convert microns to cm.
    Float_t flag=kconst*Dx(); // error value
    Double_t th=0.0,dx,dz,i,a,b=0.0,xb[4],zb[4];
    Float_t stereoP, stereoN;
    Angles(stereoP,stereoN);

    z = 0.0;  // Strip center in z.
    if(iPN<0 || iPN>1){// if error return full detector size in x.
	x = z = flag; 
	return;
    } // end if
    if(ix<0 || ix>=fNstrips) { // if error return full detector size in x.
	x = z = flag;
	return;
    } // end if
    i  = (Double_t) ix;      // convert to double
    dx = 0.5*kconst*Dx();    // half distance in x in cm
    dz = 0.5*kconst*Dz();    // half distance in z in cm
    a  = kconst*Dpx(ix)*(i+0.5)-dx; // Min x value.
    if(iPN==0){ //P-side angle defined backwards.
	th = TMath::Tan(stereoP); 
	b  = dz*th;
    }else if(iPN==1){ // N-side
	 th = TMath::Tan(-stereoN);
	 b  = -dz*th;
    } // end if
    // compute average/center position of the strip.
    xb[0] = +dx; if(th!=0.0) zb[0] = (+dx-a-b)/th; else zb[0] = 0.0;
    xb[1] = -dx; if(th!=0.0) zb[1] = (-dx-a-b)/th; else zb[1] = 0.0;
    xb[2] = a+b+dz*th; zb[2] = +dz;
    xb[3] = a+b-dz*th; zb[3] = -dz;
    x = 0.0; z = 0.0;
    for(Int_t j=0;j<4;j++){
	if(xb[j]>=-dx && xb[j]<=dx && zb[j]>=-dz && zb[j]<=dz){
	    x += xb[j];
	    z += zb[j];
	} // end if
    } // end for
    x *= 0.5;
    z *= 0.5;
    return;
}
//----------------------------------------------------------------------
Int_t AliITSsegmentationSSD::GetChipFromChannel(Int_t ix, Int_t iz) const {
  // returns chip number (in range 0-11) starting from channel number

  if( (iz>=fgkNstripsDefault) || (iz<0) || (ix<0) || (ix>1) ) {
    AliError("Bad cell number");
    return -1;
  }
  
  if(ix==1) iz = 1535-iz;
  Int_t theChip =iz/fgkNstripsPerChip;
  return theChip;

}
//----------------------------------------------------------------------
Int_t AliITSsegmentationSSD::GetChipFromLocal(Float_t xloc, Float_t zloc) const
{
  // returns chip numbers starting from local coordinates
  // The two Nside chip number and Pside chip number are 
  // coded as chip=Nchip*10+Pchip

  Int_t iP=0;
  Int_t iN=0;
  if (!LocalToDet(xloc,zloc,iP,iN) || 
      (iP<0) || (iP>=fNstrips) || (iN<0) || (iN>=fNstrips) ) {
    //AliWarning("Bad local coordinate");
    return -1;
  }

  Int_t iChip = GetChipFromChannel(0,iP);
  iChip += 10*GetChipFromChannel(1,iN); // add Nside

  return iChip;

}
//

Int_t AliITSsegmentationSSD::GetChipsInLocalWindow(Int_t* array, Float_t zmin, Float_t zmax, 
						   Float_t xmin, Float_t xmax) const {
  // returns chip number in a given xz window

  Int_t nChipInW = 0;

  Float_t zminDet=-fDz*1.0E-04/2.;
  Float_t zmaxDet=fDz*1.0E-04/2.;
  if(zmin<zminDet) zmin=zminDet;
  if(zmax>zmaxDet) zmax=zmaxDet;

  Float_t xminDet=-fDx*1.0E-04/2;
  Float_t xmaxDet=fDx*1.0E-04/2;
  if(xmin<xminDet) xmin=xminDet;
  if(xmax>xmaxDet) xmax=xmaxDet;

  Int_t n1N=-1;
  Int_t n1P=-1;
  Int_t n1=GetChipFromLocal(xmin,zmin); 
  if(n1!=-1) { // Note! Recpoint can be on the sensor but in the dead area not covered by strips!
    n1N = (Int_t) (n1/10); // N-side chip coded as 10*chip_index
    n1P = n1 - 10 * n1N; // P-side chip coded 0-5
    array[nChipInW]=n1P;
    nChipInW++;
    array[nChipInW]=n1N;
    nChipInW++;
  }
  
  Int_t n2N=-1;
  Int_t n2P=-1;
  Int_t n2=GetChipFromLocal(xmin,zmax);
  if(n2!=-1) { // Note! Recpoint can be on the sensor but in the dead area not covered by strips!
    n2N = (Int_t) (n2/10); // N-side chip coded as 10*chip_index
    n2P = n2 - 10 * n2N; // P-side chip coded 0-5
    if(n2P!=n1P) { array[nChipInW]=n2P; nChipInW++;}
    if(n2N!=n1N) { array[nChipInW]=n2N; nChipInW++;}
  }

  Int_t n3N=-1;
  Int_t n3P=-1;
  Int_t n3=GetChipFromLocal(xmax,zmin);
  if(n3!=-1) {
    n3N=(Int_t) (n3/10); // N-side chip coded as 10*chip_index
    n3P=n3 - 10 * n3N; // P-side chip coded 0-5
    if((n3P!=n1P)&&(n3P!=n2P)) { array[nChipInW]=n3P; nChipInW++;}
    if((n3N!=n1N)&&(n3N!=n2N)) { array[nChipInW]=n3N; nChipInW++;}
  }
  
  Int_t n4N=-1;
  Int_t n4P=-1;
  Int_t n4=GetChipFromLocal(xmax,zmax);
  if(n4!=-1) {
    n4N=(Int_t) (n4/10); // N-side chip coded as 10*chip_index
    n4P=n4 - 10 * n4N; // P-side chip coded 0-5
    if((n4P!=n1P)&&(n4P!=n2P)&&(n4P!=n3P)) { array[nChipInW]=n4P; nChipInW++;}
    if((n4N!=n1N)&&(n4N!=n2N)&&(n4N!=n3N)) { array[nChipInW]=n4N; nChipInW++;}
  }
  
  return nChipInW;

}

//----------------------------------------------------------------------
void AliITSsegmentationSSD::PrintDefaultParameters() const {
// Print default values for parameters. 
// Values specified as static const data members are shown

  cout<<"fgkDxDefault = "<<fgkDxDefault<<endl;
  cout<<"fgkDzDefault = "<<fgkDzDefault<<endl;
  cout<<"fgkDyDefault = "<<fgkDyDefault<<endl;
  cout<<"fgkPitchDefault = "<<fgkPitchDefault<<endl;
  cout<<"fgkNstripsDefault = "<<fgkNstripsDefault<<endl;
}
