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
Revision 1.10  2000/07/10 20:57:39  hristov
Update of TPC code and macros by M.Kowalski

Revision 1.9  2000/06/30 12:07:50  kowal2
Updated from the TPC-PreRelease branch

Revision 1.8.4.4  2000/06/26 07:39:42  kowal2
Changes to obey the coding rules

Revision 1.8.4.3  2000/06/25 08:38:41  kowal2
Splitted from AliTPCtracking
  
Revision 1.8.4.2  2000/06/14 16:48:24  kowal2
Parameter setting improved. Removed compiler warnings

Revision 1.8.4.1  2000/06/09 07:12:21  kowal2  

Updated defaults

Revision 1.8  2000/04/17 09:37:33  kowal2
removed obsolete AliTPCDigitsDisplay.C

Revision 1.7.8.2  2000/04/10 08:44:51  kowal2

New transformations added
Different pad and pad-rows geometries for different sectors

Revision 1.7.8.1  2000/04/10 07:56:53  kowal2
Not used anymore - removed

Revision 1.7  1999/10/08 13:10:35  fca
Values in SetDefault are in radiants

Revision 1.6  1999/10/08 06:27:59  fca
Defaults updated

Revision 1.5  1999/10/05 17:18:27  fca
Correct GetWire check on even/odd fnWires

Revision 1.4  1999/09/29 09:24:34  fca
Introduction of the Copyright and cvs Log

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

//

#include <iostream.h>
#include <TMath.h>
#include <TObject.h>
#include <TRandom.h>
#include <AliTPCParam.h>




ClassImp(AliTPCParam)


//___________________________________________
AliTPCParam::AliTPCParam()
{   
  //
  //constructor sets the default parameters
  //

  fResponseBin = 0;
  fResponseWeight = 0;
  fRotAngle = 0;
  SetTitle("75x40_100x60");
  SetDefault();  
}

AliTPCParam::~AliTPCParam()
{
  //
  //destructor deletes some dynamicaly alocated variables
  //

  if (fResponseBin!=0)    delete [] fResponseBin;
  if (fResponseWeight!=0) delete [] fResponseWeight;
  if (fRotAngle      !=0) delete [] fRotAngle;

}




Int_t  AliTPCParam::Transform0to1(Float_t *xyz, Int_t * index)  const
{
  //
  // calculates sector number (index[1], undefined on input)
  // xyz intact
  //

  Float_t angle,x1;
  Int_t sector;
  Float_t r = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  if ((xyz[0]==0)&&(xyz[1]==0)) angle = 0.;
  else
    {
      angle =TMath::ASin(xyz[1]/r);
      if   (xyz[0]<0)   angle=TMath::Pi()-angle;
      if ( (xyz[0]>0) && (xyz[1]<0) ) angle=2*TMath::Pi()+angle;
    }

  sector=Int_t((angle-fInnerAngleShift)/fInnerAngle);      
 
  Float_t cos,sin;
  AdjustCosSin(sector,cos,sin);
  x1=xyz[0]*cos + xyz[1]*sin;

  if (x1>fOuterRadiusLow)
    {
      sector=Int_t((angle-fOuterAngleShift)/fOuterAngle)+fNInnerSector;      
      if (xyz[2]<0) 	sector+=(fNOuterSector>>1);            
    }
    else   
      if (xyz[2]<0) sector+=(fNInnerSector>>1);    
  index[1]=sector; // calculated sector number
  index[0]=1; // indicates system after transformation
  return sector;
}

Bool_t  AliTPCParam::Transform(Float_t *xyz, Int_t *index, Int_t* oindex)
{
  //transformation from input coodination system to output coordination system
  switch (index[0]){
  case 0:
    break;
  };

  return kFALSE;

}

Int_t AliTPCParam::GetPadRow(Float_t *xyz, Int_t *index) const 
{
  //
  //calculates pad row of point xyz - transformation to system 8 (digit system)
  //
  Int_t system = index[0];
  if (0==system) {
    Transform0to1(xyz,index); 
    system=1;
  }
  if (1==system) {
    Transform1to2(xyz,index); 
    system=2;
  }
    
  if (fGeometryType==0){ //straight row    
    if (2==system) {
      Transform2to3(xyz,index);       
      system=3;
    } 
    if (3==system) {
      Transform3to4(xyz,index);
      system=4; 
    }
    if (4==system) {
      Transform4to8(xyz,index);
      system=8;     
    }
    if (8==system) {
      index[0]=8;
      return index[2];
    } 
  }

  if (fGeometryType==1){ //cylindrical geometry    
    if (2==system) {
      Transform2to5(xyz,index);       
      system=5;
    } 
    if (5==system) {
      Transform2to3(xyz,index);
      system=6;
    }
    if (6==system) {
      Transform3to4(xyz,index); 
      system=7;
    }
    if (8==system) {
      index[0]=8;
      return index[2];
    }
  } 
  index[0]=system;
  return -1; //if no reasonable system     
}

void  AliTPCParam::SetSectorAngles(Float_t innerangle, Float_t innershift, Float_t outerangle,
			Float_t outershift)
{
  //
  // set opening angles  
  const static  Float_t  kDegtoRad = 0.01745329251994;
  fInnerAngle = innerangle;       //opening angle of Inner sector
  fInnerAngleShift = innershift;  //shift of first inner sector center to the 0
  fOuterAngle = outerangle;       //opening angle of outer sector
  fOuterAngleShift = outershift;  //shift of first sector center to the 0  
  fInnerAngle *=kDegtoRad;
  fInnerAngleShift *=kDegtoRad;
  fOuterAngle *=kDegtoRad;
  fOuterAngleShift *=kDegtoRad;
}

Float_t  AliTPCParam::GetInnerAngle() const
{
  //return angle 
  return fInnerAngle;

}

Float_t  AliTPCParam::GetInnerAngleShift() const
{  
  //return angle   
  return fInnerAngleShift;  
}
Float_t  AliTPCParam::GetOuterAngle() const
{ 
  //return angle 
  return fOuterAngle;
} 
Float_t  AliTPCParam::GetOuterAngleShift() const
{ 
  //return angle 

     return fOuterAngleShift;
} 


Int_t AliTPCParam::GetIndex(Int_t sector, Int_t row)
{
  //
  //give index of the given sector and pad row 
  //no control if the sectors and rows  are reasonable !!!
  //
  if (sector<fNInnerSector) return sector*fNRowLow+row;
  return (fNInnerSector*fNRowLow)+(sector-fNInnerSector)*fNRowUp+row;  
}

Bool_t   AliTPCParam::AdjustSectorRow(Int_t index, Int_t & sector, Int_t &row) const
{
  //
  //return sector and padrow for given index
  //if index is reasonable returns true else return false
  //
  if ( (index<0) || (index>fNtRows))  return kFALSE;
  Int_t outindex = fNInnerSector*fNRowLow;
  if (index<outindex) {
    sector = index/fNRowLow;
    row    = index - sector*fNRowLow;
    return kTRUE;
  }
  index-= outindex;
  sector = index/fNRowUp;
  row    = index - sector*fNRowUp;
  sector += fNInnerSector;
  return kTRUE;         
} 

void AliTPCParam::SetDefault()
{
  //
  //set default parameters
  //
  //const static  Int_t kMaxRows=600; 
  //
  //sector default parameters
  //
  const static  Float_t kInnerRadiusLow = 82.97;
  const static  Float_t kOuterRadiusLow = 133.58;
  const static  Float_t kInnerRadiusUp  = 133.17;
  const static  Float_t kOuterRadiusUp  = 247.78;
  const static  Float_t kInnerAngle = 20; // 20 degrees
  const static  Float_t kInnerAngleShift = 10;
  const static  Float_t kOuterAngle = 20; //  20 degrees
  const static  Float_t kOuterAngleShift = 10;
  const static  Float_t kInnerFrameSpace = 1.5;
  const static  Float_t kOuterFrameSpace = 1.5;
  const static  Float_t kInnerWireMount = 1.370825926;
  const static  Float_t kOuterWireMount = 1.370825926;
  const static  Float_t kZLength =250.;
  const static  Int_t   kGeometryType = 0; //straight rows 
  //
  //wires default parameters
  //
  const static Int_t    kNInnerWiresPerPad = 3;
  const static Int_t    kInnerDummyWire = 2;
  const static Float_t  kInnerOffWire = 0.5;
  const static Int_t    kNOuterWiresPerPad = 4;
  const static Int_t    kOuterDummyWire = 2;
  const static Float_t  kOuterOffWire = 0.5;
  //
  //pad default parameters
  // 
  const static Float_t  kInnerPadPitchLength = 0.75;
  const static Float_t  kInnerPadPitchWidth = 0.40;
  const static Float_t  kInnerPadLength = 0.75;
  const static Float_t  kInnerPadWidth = 0.40;
  const static Float_t  kOuterPadPitchLength = 1.0;
  const static Float_t  kOuterPadPitchWidth = 0.6;
  const static Float_t  kOuterPadLength = 1.0;
  const static Float_t  kOuterPadWidth = 0.6;
  const static Bool_t   kBMWPCReadout = kTRUE; //MWPC readout - another possibility GEM 
  const static Int_t    kNCrossRows = 1; //number of rows to cross-talk
  
  //
  //gas default parameters
  //
  const static  Float_t  kDiffT = 2.2e-2; 
  const static  Float_t  kDiffL = 2.2e-2;
  const static  Float_t  kGasGain = 2.e4;
  const static  Float_t  kDriftV  =2.83e6;
  const static  Float_t  kOmegaTau = 0.145;
  const static  Float_t  kAttCoef = 250.;
  const static  Float_t  kOxyCont = 5.e-6;
  //
  //electronic default parameters
  //
  const static  Float_t  kPadCoupling=0.5;
  const static  Int_t    kZeroSup=2;
  const static  Float_t  kNoise = 1000;                            
  const static  Float_t  kChipGain = 12;
  const static  Float_t  kChipNorm = 0.4;
  const static  Float_t  kTSample = 2.e-7; 
  const static  Float_t  kTFWHM   = 1.9e-7;  //fwhm of charge distribution
  const static  Int_t    kMaxTBin =445;  
  const static  Int_t    kADCSat  =1024;  
  const static  Float_t  kADCDynRange =2000.;  
  //
  //
  //
  const static  Float_t kBField =0.2; 
  const static  Float_t kNPrimLoss =10.9;
  const static  Float_t kNTotalLoss =39.9;
  // 
  //response constants
  //
  const static Int_t     kNResponseMax=100;
  const static Float_t   kResponseThreshold=0.01;     
  fbStatus = kFALSE;
  //
  //set sector parameters
  //
  SetInnerRadiusLow(kInnerRadiusLow);
  SetOuterRadiusLow(kOuterRadiusLow);
  SetInnerRadiusUp(kInnerRadiusUp);
  SetOuterRadiusUp(kOuterRadiusUp);
  SetInnerFrameSpace(kInnerFrameSpace);
  SetOuterFrameSpace(kOuterFrameSpace);
  SetInnerWireMount(kInnerWireMount);
  SetOuterWireMount(kOuterWireMount);
  SetSectorAngles(kInnerAngle,kInnerAngleShift,kOuterAngle,kOuterAngleShift);
  SetZLength(kZLength);
  SetGeometryType(kGeometryType);
  //
  //set wire parameters
  //
  SetInnerNWires(kNInnerWiresPerPad);
  SetInnerDummyWire(kInnerDummyWire);
  SetInnerOffWire(kInnerOffWire);
  SetOuterNWires(kNOuterWiresPerPad);
  SetOuterDummyWire(kOuterDummyWire);
  SetOuterOffWire(kOuterOffWire);
  //
  //set pad parameter
  //
  SetInnerPadPitchLength(kInnerPadPitchLength);
  SetInnerPadPitchWidth(kInnerPadPitchWidth);
  SetInnerPadLength(kInnerPadLength);
  SetInnerPadWidth(kInnerPadWidth);
  SetOuterPadPitchLength(kOuterPadPitchLength);
  SetOuterPadPitchWidth(kOuterPadPitchWidth);
  SetOuterPadLength(kOuterPadLength);
  SetOuterPadWidth(kOuterPadWidth); 
  SetMWPCReadout(kBMWPCReadout);
  SetNCrossRows(kNCrossRows);
  //
  //set gas paremeters
  //
  SetDiffT(kDiffT);
  SetDiffL(kDiffL);
  SetGasGain(kGasGain);
  SetDriftV(kDriftV);
  SetOmegaTau(kOmegaTau);
  SetAttCoef(kAttCoef);
  SetOxyCont(kOxyCont);
  //
  //set electronivc parameters  
  //
  SetPadCoupling(kPadCoupling);
  SetZeroSup(kZeroSup);
  SetNoise(kNoise);
  SetChipGain(kChipGain);
  SetChipNorm(kChipNorm);   
  SetTSample(kTSample);
  SetTFWHM(kTFWHM);
  SetMaxTBin(kMaxTBin);
  SetADCSat(kADCSat);
  SetADCDynRange(kADCDynRange);
  //set magnetic field
  SetBField(kBField);
  SetNPrimLoss(kNPrimLoss);
  SetNTotalLoss(kNTotalLoss);
  //
  //set response  parameters  
  //
  SetNResponseMax(kNResponseMax); 
  SetResponseThreshold(kResponseThreshold);
}

          
Bool_t AliTPCParam::Update()
{
  //
  // update some calculated parameter which must be updated after changing "base"
  // parameters 
  // for example we can change size of pads and according this recalculate number
  // of pad rows, number of of pads in given row ....
  //
  const Float_t kQel = 1.602e-19; // elementary charge
  fbStatus = kFALSE;

  Int_t i,j;  //loop variables because HP 
  //-----------------Sector section------------------------------------------
  //calclulate number of sectors
  fNInnerSector = Int_t(4*TMath::Pi()/fInnerAngle+0.2); 
       // number of inner sectors - factor 0.2 to don't be influnced by inprecision
  if (fNInnerSector%2) return kFALSE;
  fNOuterSector = Int_t(4*TMath::Pi()/fOuterAngle+0.2); 
  if (fNOuterSector%2) return kFALSE;
  fNSector  = fNInnerSector+fNOuterSector;

  if (fRotAngle!=0) delete [] fRotAngle;
  fRotAngle = new Float_t[4*fNSector];
  //calculate sin and cosine of rotations angle     
  //sectors angles numbering from 0

  j=fNInnerSector*2;
  Float_t angle = fInnerAngleShift; 
  for (i=0; j<fNInnerSector*4; i+=4, j+=4 , angle +=fInnerAngle){
    fRotAngle[i]=TMath::Cos(angle);
    fRotAngle[i+1]=TMath::Sin(angle);
    fRotAngle[j] =  fRotAngle[i];
    fRotAngle[j+1] =  fRotAngle[i+1];
    fRotAngle[i+2] =angle;
    fRotAngle[j+2] =angle;    
  }
  angle = fOuterAngleShift; 
  j=(fNInnerSector+fNOuterSector/2)*4;
  for (i=fNInnerSector*4; j<fNSector*4; i+=4,j+=4, angle +=fOuterAngle){
    fRotAngle[i]=TMath::Cos(angle);
    fRotAngle[i+1]=TMath::Sin(angle);
    fRotAngle[j] =  fRotAngle[i];
    fRotAngle[j+1] =  fRotAngle[i+1];
    fRotAngle[i+2] =angle;
    fRotAngle[j+2] =angle;    
  }
  fZWidth = fTSample*fDriftV;  
  fTotalNormFac = fPadCoupling*fChipNorm*kQel*1.e15*fChipGain*fADCSat/fADCDynRange;
  fNoiseNormFac = kQel*1.e15*fChipGain*fADCSat/fADCDynRange;
  //wire section 
  Int_t nwire;
  Float_t wspace; //available space for wire
  Float_t dummyspace; //dummyspace for wire

  fInnerWWPitch = Float_t((Double_t)fInnerPadPitchLength/(Double_t)fNInnerWiresPerPad);  
  wspace =fInnerRadiusUp-fInnerRadiusLow-2*fInnerOffWire;
  nwire = Int_t(wspace/fInnerWWPitch);
  wspace = Float_t(nwire)*fInnerWWPitch;
  dummyspace =(fInnerRadiusUp-fInnerRadiusLow-wspace)/2.; 
  fRInnerFirstWire = fInnerRadiusLow+dummyspace;
  fRInnerLastWire = fRInnerFirstWire+fInnerWWPitch*(Float_t)(nwire);

  fOuterWWPitch = Float_t((Double_t)fOuterPadPitchLength/(Double_t)fNOuterWiresPerPad);  
  wspace =fOuterRadiusUp-fOuterRadiusLow-2*fOuterOffWire;
  nwire = Int_t(wspace/fOuterWWPitch);
  wspace = Float_t(nwire)*fOuterWWPitch;
  dummyspace =(fOuterRadiusUp-fOuterRadiusLow-wspace)/2.; 
  fROuterFirstWire = fOuterRadiusLow+dummyspace;
  fROuterLastWire = fROuterFirstWire+fOuterWWPitch*(Float_t)(nwire);

  
  //
  //response data
  //
  if (fResponseBin==0) delete [] fResponseBin;
  if (fResponseWeight==0) delete [] fResponseBin;
  fResponseBin    = new Int_t[3*fNResponseMax];
  fResponseWeight = new Float_t[fNResponseMax];
  
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
  return fNRowLow;
}
Int_t AliTPCParam::GetNRowUp() const
{
  //get the number of pad rows in up sector
  return fNRowUp;
}
Float_t AliTPCParam::GetPadRowRadiiLow(Int_t irow) const
{
  //get the pad row (irow) radii
  if ( !(irow<0) && (irow<fNRowLow) ) 
    return  fPadRowLow[irow];
  else
    return 0;
}

Float_t AliTPCParam::GetPadRowRadiiUp(Int_t irow) const
{
  //get the pad row (irow) radii
 if ( !(irow<0) && (irow<fNRowUp) ) 
    return  fPadRowUp[irow];
  else
    return 0;
}

Int_t AliTPCParam::GetNPadsLow(Int_t irow) const
{
  //get the number of pads in row irow
  if ( !(irow<0) && (irow<fNRowLow) ) 
    return  fNPadsLow[irow];
  else
    return 0;
}


Int_t AliTPCParam::GetNPadsUp(Int_t irow) const
{
  //get the number of pads in row irow
  if ( !(irow<0) && (irow<fNRowUp) ) 
    return  fNPadsUp[irow];
  else
    return 0;
}

