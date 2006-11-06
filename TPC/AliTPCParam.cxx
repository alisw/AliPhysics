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

#include <AliTPCParam.h>

#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include "AliAlignObj.h"
#include "AliAlignObjAngles.h"
#include "AliLog.h"

ClassImp(AliTPCParam)


//___________________________________________
AliTPCParam::AliTPCParam()
            :AliDetectorParam(),
	     fbStatus(kFALSE),
             fInnerRadiusLow(0.),
             fInnerRadiusUp(0.),
             fOuterRadiusUp(0.),
             fOuterRadiusLow(0.),
	     fInnerAngle(0.),
	     fInnerAngleShift(0.),
	     fOuterAngle(0.),
	     fOuterAngleShift(0.),
	     fInnerFrameSpace(0.),
	     fOuterFrameSpace(0.),
	     fInnerWireMount(0.),
	     fOuterWireMount(0.),
	     fNInnerSector(0),
	     fNOuterSector(0),
	     fNSector(0),
	     fZLength(0),
	     fRotAngle(),
	     fGeometryType(0),
	     fTrackingMatrix(0),
	     fClusterMatrix(0), 
	     fGlobalMatrix(0),
	     fNInnerWiresPerPad(0),
	     fInnerWWPitch(0),
	     fInnerDummyWire(0),
	     fInnerOffWire(0.),
	     fRInnerFirstWire(0.),
	     fRInnerLastWire(0.),
	     fLastWireUp1(0.),
	     fNOuter1WiresPerPad(0),
	     fNOuter2WiresPerPad(0),
	     fOuterWWPitch(0.),
	     fOuterDummyWire(0),
	     fOuterOffWire(0.),
	     fROuterFirstWire(0.),
	     fROuterLastWire(0.),
	     fInnerPadPitchLength(0.),
	     fInnerPadPitchWidth(0.),
	     fInnerPadLength(0.),
	     fInnerPadWidth(0.),
	     fOuter1PadPitchLength(0.),
	     fOuter2PadPitchLength(0.),
	     fOuterPadPitchWidth(0.),
	     fOuter1PadLength(0.),
	     fOuter2PadLength(0.),
	     fOuterPadWidth(0.),
	     fBMWPCReadout(kFALSE),
	     fNCrossRows(0),
	     fNRowLow(0),
	     fNRowUp1(0),
	     fNRowUp2(0),
	     fNRowUp(0),
	     fNtRows(0),
	     fDiffT(0.),
	     fDiffL(0.),
	     fGasGain(0.),
	     fDriftV(0.),
	     fOmegaTau(0.),
	     fAttCoef(0.),
	     fOxyCont(0.),
	     fPadCoupling(0.),
	     fZeroSup(0),
	     fNoise(0.),
	     fChipGain(0.),
	     fChipNorm(0.),
	     fTSample(0.),
	     fZWidth(0.),
	     fTSigma(0.),
	     fMaxTBin(0),
	     fADCSat(0),
	     fADCDynRange(0.),
	     fTotalNormFac(0.),
	     fNoiseNormFac(0.),
	     fNResponseMax(0),
	     fResponseThreshold(0.),
	     fCurrentMax(0),
	     fResponseBin(0),
	     fResponseWeight(0),
	     fGateDelay(0.),
	     fL1Delay(0.),
	     fNTBinsBeforeL1(0),
	     fNTBinsL1(0.)   
{   
  //
  //constructor sets the default parameters
  //

  SetTitle("75x40_100x60_150x60");
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

  if (fTrackingMatrix) {
    for(Int_t i = 0; i < fNSector; i++)
      delete fTrackingMatrix[i];
    delete [] fTrackingMatrix;
  }

  if (fClusterMatrix) {
    for(Int_t i = 0; i < fNSector; i++)
      delete fClusterMatrix[i];
    delete [] fClusterMatrix;
  }

  if (fGlobalMatrix) {
    for(Int_t i = 0; i < fNSector; i++)
      delete fGlobalMatrix[i];
    delete [] fGlobalMatrix;
  }

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

  sector=Int_t(TMath::Nint((angle-fInnerAngleShift)/fInnerAngle));      
 
  Float_t cos,sin;
  AdjustCosSin(sector,cos,sin);
  x1=xyz[0]*cos + xyz[1]*sin;

  if (x1>fOuterRadiusLow)
    {
      sector=Int_t(TMath::Nint((angle-fOuterAngleShift)/fOuterAngle))+fNInnerSector;      
      if (xyz[2]<0) 	sector+=(fNOuterSector>>1);            
    }
    else   
      if (xyz[2]<0) sector+=(fNInnerSector>>1);    
  index[1]=sector; // calculated sector number
  index[0]=1; // indicates system after transformation
  return sector;
}

Bool_t  AliTPCParam::Transform(Float_t */*xyz*/, Int_t *index, Int_t* /*oindex*/)
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
  static const  Float_t  kDegtoRad = 0.01745329251994;
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


Int_t AliTPCParam::GetIndex(Int_t sector, Int_t row) const
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
  static const  Float_t kInnerRadiusLow = 83.65;
  static const  Float_t kInnerRadiusUp  = 133.3;
  static const  Float_t kOuterRadiusLow = 133.5;
  static const  Float_t kOuterRadiusUp  = 247.7;
  static const  Float_t kInnerAngle = 20; // 20 degrees
  static const  Float_t kInnerAngleShift = 10;
  static const  Float_t kOuterAngle = 20; //  20 degrees
  static const  Float_t kOuterAngleShift = 10;
  static const  Float_t kInnerFrameSpace = 1.5;
  static const  Float_t kOuterFrameSpace = 1.5;
  static const  Float_t kInnerWireMount = 1.2;
  static const  Float_t kOuterWireMount = 1.4;
  static const  Float_t kZLength =250.;
  static const  Int_t   kGeometryType = 0; //straight rows 
  static const Int_t kNRowLow = 63;
  static const Int_t kNRowUp1 = 64;
  static const Int_t kNRowUp2 = 32;
  static const Int_t  kNRowUp = 96;
  //
  //wires default parameters
  //
  static const Int_t    kNInnerWiresPerPad = 3;
  static const Int_t    kInnerDummyWire = 2;
  static const Float_t  kInnerWWPitch = 0.25;
  static const Float_t  kRInnerFirstWire = 84.475;
  static const Float_t  kRInnerLastWire = 132.475;
  static const Float_t  kInnerOffWire = 0.5;
  static const Int_t    kNOuter1WiresPerPad = 4;
  static const Int_t    kNOuter2WiresPerPad = 6;
  static const Float_t  kOuterWWPitch = 0.25;  
  static const Float_t  kROuterFirstWire = 134.225;
  static const Float_t  kROuterLastWire = 246.975;
  static const Int_t    kOuterDummyWire = 2;
  static const Float_t  kOuterOffWire = 0.5;
  //
  //pad default parameters
  // 
  static const Float_t  kInnerPadPitchLength = 0.75;
  static const Float_t  kInnerPadPitchWidth = 0.40;
  static const Float_t  kInnerPadLength = 0.75;
  static const Float_t  kInnerPadWidth = 0.40;
  static const Float_t  kOuter1PadPitchLength = 1.0;
  static const Float_t  kOuterPadPitchWidth = 0.6;
  static const Float_t  kOuter1PadLength = 1.0;
  static const Float_t  kOuterPadWidth = 0.6;
  static const Float_t  kOuter2PadPitchLength = 1.5;
  static const Float_t  kOuter2PadLength = 1.5;

  static const Bool_t   kBMWPCReadout = kTRUE; //MWPC readout - another possibility GEM 
  static const Int_t    kNCrossRows = 1; //number of rows to cross-talk
  
  //
  //gas default parameters
  //
  static const  Float_t  kDiffT = 2.2e-2; 
  static const  Float_t  kDiffL = 2.2e-2;
  static const  Float_t  kGasGain = 2.e4;
  static const  Float_t  kDriftV  =2.83e6;
  static const  Float_t  kOmegaTau = 0.145;
  static const  Float_t  kAttCoef = 250.;
  static const  Float_t  kOxyCont = 5.e-6;
  //
  //electronic default parameters
  //
  static const  Float_t  kPadCoupling=0.5;
  static const  Int_t    kZeroSup=2;
  static const  Float_t  kNoise = 1000;                            
  static const  Float_t  kChipGain = 12;
  static const  Float_t  kChipNorm = 0.4;
  static const  Float_t  kTSample = 2.e-7; 
  static const  Float_t  kTFWHM   = 1.9e-7;  //fwhm of charge distribution
  static const  Int_t    kMaxTBin =445;  
  static const  Int_t    kADCSat  =1024;  
  static const  Float_t  kADCDynRange =2000.;  
  // 
  //response constants
  //
  static const Int_t     kNResponseMax=100;
  static const Float_t   kResponseThreshold=0.01;     
  //L1 constants
  //  static const Float_t   kGateDelay=6.1e-6; //In s
  static const Float_t   kGateDelay=0.; //For the moment no gating
  //  static const Float_t   kL1Delay=6.5e-6; //In s
  static const Float_t   kL1Delay=0.; //For the moment no delay
  //  static const UShort_t  kNTBinsBeforeL1=14;
  static const UShort_t  kNTBinsBeforeL1=0; //For the moment no shift
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
  SetRowNLow(kNRowLow);
  SetRowNUp1 (kNRowUp1);
  SetRowNUp2(kNRowUp2);
  SetRowNUp(kNRowUp);
  //
  //set wire parameters
  //
  SetInnerNWires(kNInnerWiresPerPad);
  SetInnerDummyWire(kInnerDummyWire);
  SetInnerOffWire(kInnerOffWire);
  SetOuter1NWires(kNOuter1WiresPerPad);
  SetOuter2NWire(kNOuter2WiresPerPad);
  SetOuterDummyWire(kOuterDummyWire);
  SetOuterOffWire(kOuterOffWire);
  SetInnerWWPitch(kInnerWWPitch);
  SetRInnerFirstWire(kRInnerFirstWire);
  SetRInnerLastWire(kRInnerLastWire);
  SetOuterWWPitch(kOuterWWPitch);
  SetROuterFirstWire(kROuterFirstWire);
  SetROuterLastWire(kROuterLastWire);  
  //
  //set pad parameter
  //
  SetInnerPadPitchLength(kInnerPadPitchLength);
  SetInnerPadPitchWidth(kInnerPadPitchWidth);
  SetInnerPadLength(kInnerPadLength);
  SetInnerPadWidth(kInnerPadWidth);
  SetOuter1PadPitchLength(kOuter1PadPitchLength); 
  SetOuter2PadPitchLength(kOuter2PadPitchLength);
  SetOuterPadPitchWidth(kOuterPadPitchWidth);
  SetOuter1PadLength(kOuter1PadLength);
  SetOuter2PadLength(kOuter2PadLength);
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
//   //set magnetic field
//   SetBField(kBField);
//   SetNPrimLoss(kNPrimLoss);
//   SetNTotalLoss(kNTotalLoss);
  //
  //set response  parameters  
  //
  SetNResponseMax(kNResponseMax); 
  SetResponseThreshold(static_cast<int>(kResponseThreshold));
  //L1 data
  SetGateDelay(kGateDelay);
  SetL1Delay(kL1Delay);
  SetNTBinsBeforeL1(kNTBinsBeforeL1);
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
  /*  Int_t nwire;
  Float_t wspace; //available space for wire
  Float_t dummyspace; //dummyspace for wire
 
  wspace =fInnerRadiusUp-fInnerRadiusLow-2*fInnerOffWire;
  nwire = Int_t(wspace/fInnerWWPitch);
  wspace = Float_t(nwire)*fInnerWWPitch;
  dummyspace =(fInnerRadiusUp-fInnerRadiusLow-wspace)/2.;  
  wspace =fOuterRadiusUp-fOuterRadiusLow-2*fOuterOffWire;
  nwire = Int_t(wspace/fOuterWWPitch);
  wspace = Float_t(nwire)*fOuterWWPitch;
  dummyspace =(fOuterRadiusUp-fOuterRadiusLow-wspace)/2.; 
  fROuterFirstWire = fOuterRadiusLow+dummyspace;
  fROuterLastWire = fROuterFirstWire+fOuterWWPitch*(Float_t)(nwire);
  */
  
  //
  //response data
  //
  if (fResponseBin==0) delete [] fResponseBin;
  if (fResponseWeight==0) delete [] fResponseBin;
  fResponseBin    = new Int_t[3*fNResponseMax];
  fResponseWeight = new Float_t[fNResponseMax];

  //L1 data
  fNTBinsL1 = fL1Delay/fTSample - (Float_t)fNTBinsBeforeL1;
  fbStatus = kTRUE;
  return kTRUE;
}



Bool_t AliTPCParam::ReadGeoMatrices(){
  //
  // read geo matrixes
  //
  if (!gGeoManager){
    AliFatal("Geo manager not initialized\n");
  }
  AliAlignObjAngles o;
  //
  if (fTrackingMatrix) delete [] fTrackingMatrix;
  fTrackingMatrix = new TGeoHMatrix*[fNSector];
  if (fClusterMatrix) delete [] fClusterMatrix;
  fClusterMatrix = new TGeoHMatrix*[fNSector];
  if (fGlobalMatrix) delete [] fGlobalMatrix;
  fGlobalMatrix = new TGeoHMatrix*[fNSector];
  //
  for (Int_t isec=0; isec<fNSector; isec++) {
    fGlobalMatrix[isec] = 0;
    fClusterMatrix[isec]= 0;
    fTrackingMatrix[isec]=0;   
    AliAlignObj::ELayerID iLayer;
    Int_t iModule;

    if(isec<fNInnerSector) {
      iLayer = AliAlignObj::kTPC1;
      iModule = isec;
    }
    else {
      iLayer = AliAlignObj::kTPC2;
      iModule = isec - fNInnerSector;
    }

    UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
    const char *symname = AliAlignObj::SymName(volid);
    TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
    const char *path = symname;
    if(pne) path=pne->GetTitle();
    if (!gGeoManager->cd(path)) return kFALSE;      
    TGeoHMatrix *m = gGeoManager->GetCurrentMatrix();
 
    //
    TGeoRotation mchange; 
    mchange.RotateY(90); mchange.RotateX(90);
    Float_t ROCcenter[3]; 
    GetChamberCenter(isec,ROCcenter);
    //
    // Convert to global coordinate system
    //
    fGlobalMatrix[isec] = new TGeoHMatrix(*m);
    fGlobalMatrix[isec]->Multiply(&(mchange.Inverse()));
    TGeoTranslation center("center",-ROCcenter[0],-ROCcenter[1],-ROCcenter[2]);
    fGlobalMatrix[isec]->Multiply(&center);
    //
    //  cluster correction matrix
    //
    fClusterMatrix[isec] = new TGeoHMatrix;
    Double_t sectorAngle = 20.*(isec%18)+10;
    TGeoHMatrix  rotMatrix;
    rotMatrix.RotateZ(sectorAngle);
    if (GetGlobalMatrix(isec)->GetTranslation()[2]>0){
      //
      // mirrored system 
      //
      TGeoRotation mirrorZ;
      mirrorZ.SetAngles(90,0,90,90,180,0);
      fClusterMatrix[isec]->Multiply(&mirrorZ);
    }
    TGeoTranslation trans(0,0,GetZLength());
    fClusterMatrix[isec]->MultiplyLeft(&trans);
    fClusterMatrix[isec]->MultiplyLeft((GetGlobalMatrix(isec)));	
    fClusterMatrix[isec]->MultiplyLeft(&(rotMatrix.Inverse()));
  }
  return kTRUE;
}


Bool_t AliTPCParam::GetStatus() const
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
Int_t AliTPCParam::GetNRowUp1() const
{
  //get the number of pad rows in up1 sector
  return fNRowUp1;
}
Int_t AliTPCParam::GetNRowUp2() const
{
  //get the number of pad rows in up2 sector
  return fNRowUp2;
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
Float_t AliTPCParam::GetYInner(Int_t irow) const
{
  return fYInner[irow];
}


Float_t AliTPCParam::GetYOuter(Int_t irow) const
{
  return fYOuter[irow];
}

Int_t AliTPCParam::GetSectorIndex(Float_t angle, Int_t row, Float_t z) const
{
  // returns the sector index
  // takes as input the angle, index of the pad row and z position
  if(row<0) return -1;

  if (angle > 2.*TMath::Pi()) angle -= 2.*TMath::Pi();
  if (angle < 0.            ) angle += 2.*TMath::Pi();
 
  Int_t sector;
  if(row<fNRowLow) {
    sector=Int_t(TMath::Nint((angle-fInnerAngleShift)/fInnerAngle));
    if (z<0) sector += (fNInnerSector>>1);
  }
  else {
    sector=Int_t(TMath::Nint((angle-fOuterAngleShift)/fOuterAngle))+fNInnerSector;      
    if (z<0) sector += (fNOuterSector>>1);
  }    
  
  return sector;
}

Float_t AliTPCParam::GetChamberCenter(Int_t isec, Float_t * center) const
{
  // returns the default radial position
  // of the readout chambers

  const Float_t kROCcenterIn = 110.2;
  const Float_t kROCcenterOut = 188.45;

  if (isec<fNInnerSector){
    if (center){
      center[0] = kROCcenterIn;
      center[1] = 0; 
      center[2] = -5.51; 
    }
    return kROCcenterIn;
  }
  else{
    if (center){
      center[0] = kROCcenterOut;
      center[1] = 0; 
      center[2] = -5.61; 
    }
    return kROCcenterOut;
  }
}



