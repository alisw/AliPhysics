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

/*
$Log$
*/

///////////////////////////////////////////////////////////////////////
//  Manager and of geomety  classes for set: TPC                     //
//                                                                   //
//  !sectors are numbered from  0                                     //
//  !pad rows are numbered from 0                                     //
//  
//  12.6.   changed z relative 
//  Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk // 
//                                                                   //  
///////////////////////////////////////////////////////////////////////


#include <iostream.h>
#include <TMath.h>
#include <TObject.h>
#include "AliTPCSecGeo.h"
#include <AliTPCParam.h>


ClassImp(AliTPCParam)

const static  Int_t kMaxRows=600;


// default values  
const static   Int_t kMaxTBin =512; 


const static  Float_t kInnerRadiusLow = 89.45;
const static  Float_t kOuterRadiusLow = 143.725;
const static  Float_t kInnerRadiusUp  = 134.55;
const static  Float_t kOuterRadiusUp  = 248.275;

const static  Float_t kInnerAngle = 0.523598775; // 30 degrees
const static  Float_t kInnerAngleShift = 0;
const static  Float_t kOuterAngle = 0.261799387; //  15 degrees
const static  Float_t kOuterAngleShift = 0;

const static Float_t kPadPitchLength = 2.05;
const static Float_t kPadPitchWidth = 0.35;
const static Float_t kPadLength = 2.05;
const static Float_t kPadWidth = 0.35;

//  Number of wires per pad and wire-wire pitch
const static Int_t knWires = 5;
const static  Float_t  kDiffT = 2.2e-2; 
const static  Float_t  kDiffL = 2.2e-2; 
const static  Float_t  kDriftV  =2.85e6;

const static  Float_t  kOmegaTau = 0.145;
const static  Float_t  kAttCoef = 250.;
const static  Float_t  kOxyCont = 5.e-6;


const static  Float_t  kChipGain = 24;
const static  Float_t  kGasGain = 1e4;
const static  Float_t  kTSample = 2.e-7; //TSAMPLE
const static  Float_t  kTFWHM   = 2.5e-7;  //fwhm of charge distribution
 
const static  Float_t  kNoise = 500;  //default noise = 1000 el 
const static  Int_t    kZeroSup=5;
const static  Float_t  kPadCoupling=0.5;
// 
const static  Float_t  kEdgeSectorSpace = 1.15;
const static  Float_t  kDegtoRad = 0.01745329251994;
const static  Float_t  kRadtoDeg = 57.29577951309;




//___________________________________________
AliTPCParam::AliTPCParam()
{   
  //constructor set the default parameters
  SetDefault();  
}


void  AliTPCParam::SetSectorAngles(Float_t innerangle, Float_t innershift, Float_t outerangle,
			Float_t outershift, Bool_t inDegree)
{
  //
  // set opening angles  
  fInnerAngle = innerangle;       //opening angle of Inner sector
  fInnerAngleShift = innershift;  //shift of first inner sector center to the 0
  fOuterAngle = outerangle;       //opening angle of outer sector
  fOuterAngleShift = outershift;  //shift of first sector center to the 0  
  if (inDegree==kTRUE){
    fInnerAngle *=kDegtoRad;
    fInnerAngleShift *=kDegtoRad;
    fOuterAngle *=kDegtoRad;
    fOuterAngleShift *=kDegtoRad;
  }    
}


void AliTPCParam::CRXYZtoXYZ(Float_t *xyz,
	       const Int_t &sector, const Int_t & padrow, Int_t option) const  
{  
  //transform relative coordinates to absolute
  Bool_t rel = ( (option&2)!=0);
  Float_t row_first; 
  row_first = (sector<=fNInnerSector) ? fPadRowLow[0] : fPadRowUp[0]; 
  if (rel==kTRUE)  //if the position is relative to pad row  
    {
      xyz[0]+=row_first;
      xyz[0]+=(Int_t) padrow*fPadPitchLength;
    }  

  xyz[2]=z_end-xyz[2];
  if (sector<fNInnerSector) {
    if ( sector>=(fNInnerSector>>1))	xyz[2]*=-1.;
  } else {
    if ( (sector-fNInnerSector) >= (fNOuterSector>>1) )    xyz[2]*=-1;
  }

  Float_t x1=xyz[0];
  Float_t y1=xyz[1];
  Float_t cos,sin;
  AdjustAngles(sector,cos,sin);
  xyz[0]= x1*cos - y1*sin;
  xyz[1]= x1*sin + y1*cos;
}

void AliTPCParam::XYZtoCRXYZ(Float_t *xyz,
			     Int_t &sector, Int_t & padrow, Int_t option)
{
   //transform global position to the position relative to the sector padrow
  //if option=0  X calculate absolute            calculate sector
  //if option=1  X           absolute            use input sector
  //if option=2  X           relative to pad row calculate sector
  //if option=3  X           relative            use input sector
  //!!!!!!!!! WE start to calculate rows from row = 0
  
  Bool_t rel = ( (option&2)!=0);  
  //option 0 and 2  means that we don't have information about sector
  //we calculate sector
  if ((option&1)==0){
    Float_t angle;
    Float_t r = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
    if ((xyz[0]==0)&&(xyz[1]==0)) angle = 0;
    else
      {
	angle =TMath::ASin(xyz[1]/r);
	if   (xyz[0]<0)   angle=TMath::Pi()-angle;
	if ( (xyz[0]>0) && (xyz[1]<0) ) angle=2*TMath::Pi()+angle;
      }
    //transform global position to the position relative to the sector padrow
    //fistly calculate xyz[0] radius  for lover sector
    //bacause in this moment we dont know in which sector we are
    sector=Int_t((angle-fInnerAngleShift)/fInnerAngle);      
    Float_t x1;
    Float_t y1;
    //firstly we suppose that we are in inner sector
    Float_t cos,sin;
    AdjustAngles(sector,cos,sin);

    x1=xyz[0]*cos + xyz[1]*sin;
    y1=-xyz[0]*sin + xyz[1]*cos;
    if (x1>fOuterRadiusLow)
      {
	sector=Int_t((angle-fOuterAngleShift)/fOuterAngle)+fNInnerSector;
	AdjustAngles(sector,cos,sin);        
	x1=xyz[0]*cos + xyz[1]*sin;
	y1=-xyz[0]*sin + xyz[1]*cos;      
	if (xyz[2]<0) 	sector+=(fNOuterSector>>1);            
      }
    else   
      if (xyz[2]<0) sector+=(fNInnerSector>>1);    

  if  (x1<fOuterRadiusLow)   
    padrow =Int_t( (x1-fPadRowLow[0])/fPadPitchLength+0.5);
  else
    padrow = Int_t( (x1-fPadRowUp[0])/fPadPitchLength+0.5);
  if (rel==kTRUE)
      if (x1<fOuterRadiusLow)   x1-=padrow*fPadPitchLength+fPadRowLow[0];
      else
	x1-=padrow*fPadPitchLength+fPadRowUp[0];  
   xyz[0]=x1;
   xyz[1]=y1;    
   xyz[2]=z_end-TMath::Abs(xyz[2]);  
  }   //endif we don't have information about sector
  else{
    //if we have information about sector
    Float_t cos,sin;
    AdjustAngles(sector,cos,sin);   
    Float_t x1;
    Float_t y1;
    //rotate to given sector
    x1=xyz[0]*cos + xyz[1]*sin;
    y1=-xyz[0]*sin + xyz[1]*cos; 
    //calculate pad row number
    if (sector<fNInnerSector) {
      padrow =Int_t( (x1-fPadRowLow[0])/fPadPitchLength+1.5)-1;
    }
    else {
      padrow =Int_t( (x1-fPadRowUp[0])/fPadPitchLength+1.5)-1;
    }
    //if we store relative position calculate position relative to pad row
    if (rel==kTRUE){
      if (sector<fNInnerSector)
	x1-=padrow*fPadPitchLength+fPadRowLow[0];
      else 
	x1-=padrow*fPadPitchLength+fPadRowUp[0];
    }      
    xyz[0]=x1;
    xyz[1]=y1;
    xyz[2]=z_end-TMath::Abs(xyz[2]);  
  }
}

void AliTPCParam::CRYZtoTimePad(const Float_t &y, const Float_t &z,
				Float_t &time, Float_t &pad,
				Int_t sector, Int_t padrow)
{
  //transform position in cm to position in time slices and pads
  Float_t  nofpads = GetNPads(sector,padrow);
  Float_t padc=(nofpads+1)/2; // this is the "central" pad for a row
  pad = y/(fPadPitchWidth)+padc;
  time=z/fZWidth;  
}
void AliTPCParam::CRTimePadtoYZ(Float_t &y, Float_t &z,
				const Float_t &time, const Float_t &pad,
				Int_t sector, Int_t padrow)
{
  //transform position in time slices and pads  to cm 
   Float_t  nofpads = GetNPads(sector,padrow);
   Float_t padc=(nofpads+1)/2; // this is the "central" pad for a row
   y=(pad-padc)*fPadPitchWidth;
   z=time*fZWidth;
}

Int_t AliTPCParam::GetWire(Float_t & x)
{
  //
  //return wire number of pad for electron at relative position x
  //to the center of the pad
  //and adjust x to the wire position
  //we suppose that if the wire number is even the center wire
  //is at center of pad
  //
  Float_t xrel= x/fWWPitch;
  if ((fnWires>>1)==0) xrel+=1;
  else  xrel+=0.5;
  Int_t nw=Int_t(xrel);
  if (xrel<0) nw-=1;
  
  x=(nw*fWWPitch);
  if ((fnWires>>1)==0) x-=fWWPitch/2.;
  return nw;
}

Int_t AliTPCParam::GetIndex(Int_t sector, Int_t row)
{
  //
  //give index of the given sector and pad row 
  //no control if the sectors and rows  are reasonable !!!
  //
  if (sector<fNInnerSector) return sector*fnRowLow+row;
  return (fNInnerSector*fnRowLow)+(sector-fNInnerSector)*fnRowUp+row;  
}

Bool_t   AliTPCParam::AdjustSectorRow(Int_t index, Int_t & sector, Int_t &row)
{
  //
  //return sector and padrow for given index
  //if index is reasonable return true else return false
  //
  if ( (index<0) || (index>fNtRows))  return kFALSE;
  Int_t outindex = fNInnerSector*fnRowLow;
  if (index<outindex) {
    sector = index/fnRowLow;
    row    = index - sector*fnRowLow;
    return kTRUE;
  }
  index-= outindex;
  sector = index/fnRowUp;
  row    = index - sector*fnRowUp;
  return kTRUE;         
} 



Int_t AliTPCParam::GetPadRow(Int_t isec, Float_t  &x)
{
  //
  //return the pad row for given x (transformed) 
  //
  Float_t row_first=GetPadRowRadii(isec,0);
  Int_t row = Int_t(( x-row_first+1.5*fPadPitchLength)/fPadPitchLength)-1;
  //Int_t will make from -0.5 0 but we want to make -1 so we add and after substract 1
  x -=row* fPadPitchLength+row_first;
  if (  (row<0)||(row>=GetNRow(isec))) return -1;
  else return row;  
}

void AliTPCParam::SetDefault()
{
  //set default TPC param   
  fbStatus = kFALSE;
  //set sector  parameters
  fInnerRadiusLow = kInnerRadiusLow;
  fOuterRadiusLow = kOuterRadiusLow;
  fInnerRadiusUp  = kInnerRadiusUp;
  fOuterRadiusUp  = kOuterRadiusUp;   
  SetSectorAngles(kInnerAngle,kInnerAngleShift, kOuterAngle, kOuterAngleShift); 
  // set default pad size and shape
  fPadPitchLength  = kPadPitchLength;
  fPadPitchWidth   = kPadPitchWidth;
  fPadLength  = kPadLength;
  fPadWidth   = kPadWidth;   
  //
  fnWires = knWires;
  fWWPitch= kPadPitchLength/Float_t(knWires);
  fDiffT  = kDiffT;
  fDiffL  = kDiffL;
  fOmegaTau = kOmegaTau;
  fOxyCont  = kOxyCont;
  fAttCoef  = kAttCoef;
  fNoise  = kNoise;
  fChipGain = kChipGain;
  fGasGain = kGasGain;
  fZeroSup= kZeroSup;
  fPadCoupling= kPadCoupling;
  fTSample =kTSample;
  fTSigma  =kTFWHM/2.35; 
  fDriftV=kDriftV;  
  fMaxTBin = kMaxTBin;
  fbStatus = Update();
}

void  AliTPCParam::AdjustAngles(Int_t isec, Float_t &cos, Float_t &sin) const
{
  //
  //set cosinus and sinus of rotation angles for sector isec
  //
  cos=fRotAngle[isec*2];
  sin=fRotAngle[isec*2+1];
}
          
Bool_t AliTPCParam::Update()
{
  //
  // update some calculated parameter which must be updated after changing "base"
  // parameters 
  // for example we can change size of pads and according this recalculate number
  // of pad rows, number of of pads in given row ....
  //
  fbStatus = kFALSE;

  Int_t i,j;  //loop variables because HP 
  //-----------------Sector section------------------------------------------
  //calclulate number of sectors
  fNInnerSector = Int_t(4*TMath::Pi()/fInnerAngle+0.2); // number of inner sectors - factor 0.2 to don't
  //be influnced by inprecision
  if (fNInnerSector%2) return kFALSE;
  fNOuterSector = Int_t(4*TMath::Pi()/fOuterAngle+0.2); 
  if (fNOuterSector%2) return kFALSE;
  fNSector  = fNInnerSector+fNOuterSector;
  //calculate sin and cosine of rotations angle     
  //sectors angles numbering from 0
  j=fNInnerSector;
  Float_t angle = fInnerAngleShift; 
  for (i=0; i<fNInnerSector*2; i+=2, j+=2 , angle +=fInnerAngle){
    fRotAngle[i]=TMath::Cos(angle);
    fRotAngle[i+1]=TMath::Sin(angle);
    fRotAngle[j] =  fRotAngle[i];
    fRotAngle[j+1] =  fRotAngle[i+1];
  }
  angle = fOuterAngleShift; 
  j=(fNInnerSector+fNOuterSector/2)*2;
  for (i=fNInnerSector*2; i<fNSector*2; i+=2,j+=2, angle +=fOuterAngle){
    fRotAngle[i]=TMath::Cos(angle);
    fRotAngle[i+1]=TMath::Sin(angle);
    fRotAngle[j] =  fRotAngle[i];
    fRotAngle[j+1] =  fRotAngle[i+1];
  }

  
  //----------------PAD section------------------------------------
  //recalculate and check some geometric parameters 
  if (0.001>fPadPitchLength){
    cout<<"ERROR !!! Small pad pitch length \n"<<flush;
    return kFALSE;
  }
  if (fPadPitchLength<fPadLength) {
    cout<<"ERROR !!! Pitch length  smaller then length of pad \n"<<flush;
    return kFALSE;
  } 
  fnRowUp   = Int_t((0.01+fOuterRadiusUp-fOuterRadiusLow)/fPadPitchLength)+1; 
  if ( kMaxRows<fnRowUp) fnRowUp = kMaxRows;
  if (1>fnRowUp) return kFALSE;

  fnRowLow   = Int_t((0.01+fInnerRadiusUp-fInnerRadiusLow)/fPadPitchLength)+1;
  if ( kMaxRows<fnRowLow) fnRowUp = kMaxRows;
  if (1>fnRowLow) return kFALSE;
  // adjust upper sectors pad row positions and pad numbers
  for (i = 0;i<fnRowUp;i++) 
    {
       Float_t x  = fOuterRadiusLow +fPadPitchLength*(Float_t)i;
       //Float_t y =  x*2*tan(alpha_up/2)-kEdgeSectorSpace;
       Float_t y = (x-0.5*fPadPitchLength)*tan(fOuterAngle/2)-kEdgeSectorSpace
       -fPadPitchWidth/2.;
       fPadRowUp[i] = x;
       fnPadsUp[i] = 1+2*(Int_t)(y/fPadPitchWidth) ;        
       
    }
  // adjust lower sectors pad row positions and pad numbers 
  for (i = 0;i<fnRowLow;i++) 
    {
       Float_t x  = fInnerRadiusLow +fPadPitchLength*(Float_t)i;
       //  Float_t y =  x*2*tan(alpha_low/2)-kEdgeSectorSpace;
       Float_t y = (x-0.5*fPadPitchLength)*tan(fInnerAngle/2)-kEdgeSectorSpace
       -fPadPitchWidth/2.;
       fPadRowLow[i] = x;
       fnPadsLow[i] = 1+2*(Int_t)(y/fPadPitchWidth) ;
         
    }

  //that variable are not writen to the file there are calculated
  //
  fWWPitch= fPadPitchLength/Float_t(fnWires);
  fZWidth = fTSample*fDriftV;  
  fNtRows = fNInnerSector*fnRowLow+fNOuterSector*fnRowUp;
  fbStatus = kTRUE;
  return kTRUE;
}



Bool_t AliTPCParam::GetStatus()
{
  //get information about object consistency
  return fbStatus;
}

Int_t AliTPCParam::GetNRowLow() const
{
  //get the number of pad rows in low sector
  return fnRowLow;
}
Int_t AliTPCParam::GetNRowUp() const
{
  //get the number of pad rows in up sector
  return fnRowUp;
}
Float_t AliTPCParam::GetPadRowRadiiLow(Int_t irow) const
{
  //get the pad row (irow) radii
  if ( !(irow<0) && (irow<fnRowLow) ) 
    return  fPadRowLow[irow];
  else
    return 0;
}

Float_t AliTPCParam::GetPadRowRadiiUp(Int_t irow) const
{
  //get the pad row (irow) radii
 if ( !(irow<0) && (irow<fnRowUp) ) 
    return  fPadRowUp[irow];
  else
    return 0;
}

Int_t AliTPCParam::GetNPadsLow(Int_t irow) const
{
  //get the number of pads in row irow
  if ( !(irow<0) && (irow<fnRowLow) ) 
    return  fnPadsLow[irow];
  else
    return 0;
}


Int_t AliTPCParam::GetNPadsUp(Int_t irow) const
{
  //get the number of pads in row irow
  if ( !(irow<0) && (irow<fnRowUp) ) 
    return  fnPadsUp[irow];
  else
    return 0;
}


void AliTPCParam::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliTPC.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);
      if (R__v < 2) return;
      //sector parameters      
      R__b >> fInnerRadiusLow;
      R__b >> fInnerRadiusUp;
      R__b >> fOuterRadiusLow;
      R__b >> fOuterRadiusUp;
      R__b >> fInnerAngle;
      R__b >> fInnerAngleShift;
      R__b >> fOuterAngle;
      R__b >> fOuterAngleShift;
      //pad parameters
      R__b >> fPadPitchLength;
      R__b >> fPadPitchWidth;
      R__b >> fPadLength;
      R__b >> fPadWidth;

      R__b >> fnWires;
      //gas parameters
      R__b >>fDiffT;
      R__b >>fDiffL;
      R__b >>fGasGain;
      R__b >>fDriftV;
      R__b >>fOmegaTau;
      R__b >>fOxyCont;
      R__b >>fAttCoef;
      
      R__b >>fPadCoupling;
      R__b >>fZeroSup;
      R__b >>fNoise;
      R__b >>fChipGain;
      
      R__b >>fTSample;
      R__b >>fTSigma;     
      //
      Update();
   } else {
      R__b.WriteVersion(AliTPCParam::IsA());
      TObject::Streamer(R__b);      
      R__b << fInnerRadiusLow;
      R__b << fInnerRadiusUp;
      R__b << fOuterRadiusLow;
      R__b << fOuterRadiusUp;
      R__b << fInnerAngle;
      R__b << fInnerAngleShift;
      R__b << fOuterAngle;
      R__b << fOuterAngleShift;

      R__b << fPadPitchLength;
      R__b << fPadPitchWidth;
      R__b << fPadLength;
      R__b << fPadWidth;

      R__b << fnWires;
      
      R__b <<fDiffT;
      R__b <<fDiffL;
      R__b <<fGasGain;
      R__b <<fDriftV;
      R__b <<fOmegaTau;
      R__b <<fOxyCont;
      R__b <<fAttCoef;


      R__b <<fPadCoupling;
      R__b <<fZeroSup;
      R__b <<fNoise;
      R__b <<fChipGain;
      
      R__b <<fTSample;
      R__b <<fTSigma;                              
   }
}

