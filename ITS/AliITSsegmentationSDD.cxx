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

#include "AliITSsegmentationSDD.h"
// #include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSgeomSDD.h"
#include "AliITSCalibration.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSresponseSDD.h"
//////////////////////////////////////////////////////
// Segmentation class for                           //
// drift detectors                                  //
//                                                  //
//////////////////////////////////////////////////////

const Float_t AliITSsegmentationSDD::fgkDxDefault = 35085.;
const Float_t AliITSsegmentationSDD::fgkDzDefault = 75264.;
const Float_t AliITSsegmentationSDD::fgkDyDefault = 300.;
const Float_t AliITSsegmentationSDD::fgkPitchDefault = 294.;
const Float_t AliITSsegmentationSDD::fgkClockDefault = 40.;
const Int_t AliITSsegmentationSDD::fgkHalfNanodesDefault = 256; 
const Int_t AliITSsegmentationSDD::fgkNsamplesDefault = 256;

ClassImp(AliITSsegmentationSDD)
//----------------------------------------------------------------------
AliITSsegmentationSDD::AliITSsegmentationSDD(AliITSgeom* geom,
					     AliITSCalibration *resp):
fNsamples(0),
fNanodes(0),
fPitch(0),
fTimeStep(0),
fDriftSpeed(0){
  // constructor
   fGeom=geom;
   AliITSCalibrationSDD* sp = (AliITSCalibrationSDD*)resp;
   fDriftSpeed=sp->GetDriftSpeed();
   fCorr=0;
   SetDetSize(fgkDxDefault,fgkDzDefault,fgkDyDefault);
   SetPadSize(fgkPitchDefault,fgkClockDefault);
   SetNPads(fgkHalfNanodesDefault,fgkNsamplesDefault);

}
//______________________________________________________________________
AliITSsegmentationSDD::AliITSsegmentationSDD() : AliITSsegmentation(),
fNsamples(0),
fNanodes(0),
fPitch(0),
fTimeStep(0),
fDriftSpeed(0){
  // Default constructor
   fDriftSpeed=0;  
   SetDetSize(fgkDxDefault,fgkDzDefault,fgkDyDefault);
   SetPadSize(fgkPitchDefault,fgkClockDefault);
   SetNPads(fgkHalfNanodesDefault,fgkNsamplesDefault);

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
fDriftSpeed(0){
  // copy constructor
  source.Copy(*this);
}

//----------------------------------------------------------------------
void AliITSsegmentationSDD::Init(){
  // Standard initilisation routine

   if(!fGeom) {
     Fatal("Init","the pointer to the ITS geometry class (AliITSgeom) is null\n");
     return;
   }
   AliITSgeomSDD *gsdd = (AliITSgeomSDD *) (fGeom->GetShape(3,1,1));

   const Float_t kconv=10000.;
   fDz = 2.*kconv*gsdd->GetDz();
   fDx = kconv*gsdd->GetDx();
   fDy = 2.*kconv*gsdd->GetDy();
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
void AliITSsegmentationSDD::GetPadIxz(Float_t x,Float_t z,
				      Int_t &timebin,Int_t &anode) const {
// Returns cell coordinates (time sample,anode) incremented by 1 !!!!! 
// for given real local coordinates (x,z)

    // expects x, z in cm

    const Float_t kconv=10000;  // cm->um

    x *= kconv; // Convert to microns
    z *= kconv; // Convert to microns
    Int_t na = fNanodes/2;
    Float_t driftpath=fDx-TMath::Abs(x);
    timebin=(Int_t)(driftpath/fDriftSpeed/fTimeStep);
    anode=(Int_t)(z/fPitch + na/2);
    if (x > 0) anode += na;

    timebin+=1;
    anode+=1;

}
//----------------------------------------------------------------------
void AliITSsegmentationSDD::GetPadCxz(Int_t timebin,Int_t anode,
				      Float_t &x ,Float_t &z) const{
    // Transform from cell to real local coordinates
    // returns x, z in cm

  // the +0.5 means that an # and time bin # should start from 0 !!! 
    const Float_t kconv=10000;  // um->cm
  // the +0.5 means that an # and time bin # should start from 0 !!! 

    Int_t na = fNanodes/2;
    Float_t driftpath=(timebin+0.5)*fTimeStep*fDriftSpeed;
    if (anode >= na) x=(fDx-driftpath)/kconv;
    else x = -(fDx-driftpath)/kconv;
    if (anode >= na) anode-=na;
    z=((anode+0.5)*fPitch-fDz/2)/kconv;

}
//----------------------------------------------------------------------
void AliITSsegmentationSDD::GetPadTxz(Float_t &x,Float_t &z) const{
    // Get anode and time bucket as floats - numbering from 0

    // expects x, z in cm

    const Float_t kconv=10000;  // cm->um

    Float_t x0=x;
    Int_t na = fNanodes/2;
    Float_t driftpath=fDx-TMath::Abs(kconv*x);
    x=driftpath/fDriftSpeed/fTimeStep;
    z=kconv*z/fPitch + (float)na/2;
    if (x0 < 0) x = -x;

}
//----------------------------------------------------------------------
void AliITSsegmentationSDD::GetLocal(Int_t module,Float_t *g ,Float_t *l) const {
  // returns local coordinates from global
    if(!fGeom) {
      Fatal("GetLocal","the pointer to the ITS geometry class (AliITSgeom) is null\n");
        return;
    }
    fGeom->GtoL(module,g,l);
}
//----------------------------------------------------------------------
void AliITSsegmentationSDD::GetGlobal(Int_t module,Float_t *l ,Float_t *g) const {
  // return global coordinates from local
    if(!fGeom) {
      Fatal("GetGlobal","the pointer to the ITS geometry class (AliITSgeom) is null\n");
    }

    fGeom->LtoG(module,l,g);

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
//     This segmentation geometry can be discribed as the following:
// {assumes 2*Dx()=7.0cm Dz()=7.5264cm, Dpx()=25ns,
//  res->DeriftSpeed()=7.3mic/ns, Dpz()=512. For other values a only the 
//  specific numbers will change not their layout.}
//
//        0                     191                    0
//      0 |----------------------|---------------------| 256
//        | a     time-bins      |     time-bins     a |
//        | n                    |                   n |
//        | o                    |___________________o_|__> X
//        | d                    |                   d |
//        | e                    |                   e |
//        | s                    |                   s |
//    255 |----------------------|---------------------| 511
//                               |
//                               V
//                               Z
    Float_t dx,dz,tb;
    const Float_t kconv = 1.0E-04; // converts microns to cm.

    ix = -1; // default values
    iz = -1; // default values
    dx = -kconv*Dx(); // lower left edge in cm.
    dz = -0.5*kconv*Dz(); // lower left edge in cm.
    if(x<dx || x>-dx) {
      Warning("LocalToDet","input argument %f out of range (%f, %f)",x,dx,-dx);
      return kFALSE; // outside of defined volume.
    }
    if(z<dz || z>-dz) {
      Warning("LocalToDet","input argument %f out of range (%f, %f)",z,dz,-dz);
      return kFALSE; // outside of defined volume.
    }
    tb = fDriftSpeed*fTimeStep*kconv; // compute size of time bin.
    if(x>0) dx = -(dx + x)/tb; // distance from + side in time bin units
    else dx = (x - dx)/tb;     // distance from - side in time bin units
    dz = (z - dz)/(kconv*fPitch); // distance in z in anode pitch units
    ix = (Int_t) dx;   // time bin
    iz = (Int_t) dz;   // anode
    if(x>0) iz += Npz()/2; // if x>0 then + side anodes values.
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
//     This segmentation geometry can be discribed as the following:
// {assumes 2*Dx()=7.0cm Dz()=7.5264cm, Dpx()=25ns,
//  res->DeriftSpeed()=7.3mic/ns, Dpz()=512. For other values a only the 
//  specific numbers will change not their layout.}
//
//        0                     191                    0
//      0 |----------------------|---------------------| 256
//        | a     time-bins      |     time-bins     a |
//        | n                    |                   n |
//        | o                    |___________________o_|__> X
//        | d                    |                   d |
//        | e                    |                   e |
//        | s                    |                   s |
//    255 |----------------------|---------------------| 511
//                               |
//                               V
//                               Z
    Int_t i,j;
    Float_t tb;
    const Float_t kconv = 1.0E-04; // converts microns to cm.

    if(iz>=Npz()/2) x = kconv*Dx(); // default value for +x side.
    else x = -kconv*Dx(); // default value for -x side.
    z = -0.5*kconv*Dz(); // default value.
    if(ix<0 || ix>=Npx()) {
      Warning("DetToLocal","input argument %d out of range (0, %d)",ix,Npx());
      return; // outside of detector
    }
    if(iz<0 || iz>=Npz()) {
      Warning("DetToLocal","input argument %d out of range (0, %d)",iz,Npz());
     return; // outside of detctor
    }
    tb = fDriftSpeed*fTimeStep*kconv; // compute size of time bin.
    if(iz>=Npz()/2) tb *= -1.0; // for +x side decrement frmo Dx().
    for(i=0;i<ix;i++) x += tb; // sum up to cell ix-1
    x += 0.5*tb; // add 1/2 of cell ix for center location.
    if(iz>=Npz()/2) iz -=Npz()/2;// If +x side don't count anodes from -x side.
    for(j=0;j<iz;j++) z += kconv*fPitch; // sum up cell iz-1
    z += 0.5*kconv*fPitch; // add 1/2 of cell iz for center location.
    return; // Found x and z, return.
}
