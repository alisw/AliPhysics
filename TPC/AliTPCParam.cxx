///////////////////////////////////////////////////////////////////////
//  Manager and of geomety  classes for set: TPC                     //
//                                                                   //
//  Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk // 
//                                                                   //  
///////////////////////////////////////////////////////////////////////

// MI change global variables for geometry 
// declaration of the global static variable 
// of AliTPCParam objects

#include <iostream.h>
#include <TMath.h>
//#include <TObject.h>
#include "AliTPCParam.h"
//some old TPC parameters in AliTPCSecGeo.h
#include "AliTPCSecGeo.h"


ClassImp(AliTPCParam)

// default values  

const  Float_t kInnerRadiusLow = 89.45;
const  Float_t kOuterRadiusLow = 143.725;
const  Float_t kInnerRadiusUp  = 134.55;
const  Float_t kOuterRadiusUp  = 248.275;

const Float_t kPadPitchLength = 2.05;
const Float_t kPadPitchWidth = 0.35;
const Float_t kPadLength = 2.05;
const Float_t kPadWidth = 0.35;
//  Number of wires per pad and wire-wire pitch
const Int_t knWires = 5;
const  Float_t  kDiffT = 2.2e-2; 
const  Float_t  kDiffL = 2.2e-2; 
const  Float_t  kDriftV  =2.85e6;

const  Float_t  kOmegaTau = 0.145;
const  Float_t  kAttCoef = 250.;
const  Float_t  kOxyCont = 5.e-6;



const  Float_t  kChipGain = 24;
const  Float_t  kGasGain = 1e4;
const  Float_t  kTSample = 2.e-7; //TSAMPLE
const  Float_t  kTFWHM   = 2.5e-7;  //fwhm of charge distribution
 
const  Float_t  kNoise = 500;  //default noise = 1000 el 
const  Int_t  kZeroSup=5;
const  Float_t  kPadCoupling=0.5;
// 
const  Float_t  kEdgeSectorSpace = 5.26;




//___________________________________________
AliTPCParam::AliTPCParam()
{   
  //constructor set the default parameters
  SetDefault();  
}


void AliTPCParam::CRXYZtoXYZ(Float_t *xyz,
	       const Int_t &sector, const Int_t & padrow, Int_t option) const  
{  
  //transform relative coordinates to absolute
  Bool_t rel = ( (option&2)!=0);
  Float_t row_first; 
  row_first = (sector<25) ? fPadRowLow[0] : fPadRowUp[0]; 
  if (rel==kTRUE)  //if we have 
    {
      xyz[0]+=row_first;
      xyz[0]+=(Int_t) padrow*fPadPitchLength;
    }  
  if (sector<25) 
    if ( sector>12)	xyz[2]*=-1.;
  else 
     if (sector>48)    xyz[2]*=-1;       
  Float_t x1=xyz[0];
  Float_t y1=xyz[1];
  Float_t cos,sin;
  AdjustAngles(sector,cos,sin);
  xyz[0]=x1*cos - y1*sin;
  xyz[1]=x1*sin + y1*cos;
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
  //fistly calculate xyz[0] "polomer for lover sector
    sector=Int_t(angle/alpha_low)+1;      
    Float_t x1;
    Float_t y1;
    //firstly we suppose that we are in inner sector
    Float_t cos,sin;
    AdjustAngles(sector,cos,sin);

    x1=xyz[0]*cos + xyz[1]*sin;
    y1=-xyz[0]*sin + xyz[1]*cos;
    if (x1>fOuterRadiusLow)
      {
	sector=Int_t(angle/alpha_up)+25;
	AdjustAngles(sector,cos,sin);        
	x1=xyz[0]*cos + xyz[1]*sin;
	y1=-xyz[0]*sin + xyz[1]*cos;      
	if (xyz[2]<0) 	sector+=24;            
      }
    else   
      if (xyz[2]<0) sector+=12;    
    if (xyz[2]<0) xyz[2]=-xyz[2];  
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
  }
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
    if (sector<25) {
      padrow =Int_t( (x1-fPadRowLow[0])/fPadPitchLength+1.5)-1;
      if ( sector>12)	xyz[2]=-xyz[2];
    }
    else {
      padrow =Int_t( (x1-fPadRowUp[0])/fPadPitchLength+1.5)-1;
      if (sector>48)    xyz[2]=-xyz[2];      
    }
    //if we store relative position calculate position relative to pad row
    if (rel==kTRUE){
      if (sector<25)
	x1-=padrow*fPadPitchLength+fPadRowLow[0];
      else 
	x1-=padrow*fPadPitchLength+fPadRowUp[0];
    }      
    xyz[0]=x1;
    xyz[1]=y1;
  }
}

void AliTPCParam::CRYZtoTimePad(const Float_t &y, const Float_t &z,
				Float_t &time, Float_t &pad,
				Int_t sector, Int_t padrow)
{
  //transform position in cm to position in time slices and pads
  Float_t  nofpads = (sector < 25) ? fnPadsLow[padrow] : fnPadsUp[padrow];
  Float_t padc=(nofpads+1)/2; // this is the "central" pad for a row
  pad = y/(fPadPitchWidth)+padc;
  time=(z_end-z)/(fDriftV*fTSample);  
  //  cout<<y<<"  "<<z<<"   "<<time<<"    "<<pad<<"  "<<
  //    sector<<"   "<<padrow<<"\n";	
}
void AliTPCParam::CRTimePadtoYZ(Float_t &y, Float_t &z,
				const Float_t &time, const Float_t &pad,
				Int_t sector, Int_t padrow)
{
  //transform position in time slices and pads  to cm 
   Float_t  nofpads = (sector < 25) ? fnPadsLow[padrow] : fnPadsUp[padrow];
   Float_t padc=(nofpads+1)/2; // this is the "central" pad for a row
   y=(pad-padc)*fPadPitchWidth;
   z=z_end-time*(fDriftV*fTSample);
   //   cout<<y<<"  "<<z<<"   "<<time<<"    "<<pad<<"  "<<
   //    sector<<"   "<<padrow<<"\n";	
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
  if (sector<25) return (sector-1)*fnRowLow+row;
  return (24*fnRowLow)+(sector-25)*fnRowUp+row;  
}

Bool_t   AliTPCParam::AdjustSectorRow(Int_t index, Int_t & sector, Int_t &row)
{
  //
  //return sector and padrow for given index
  //if index is reasonable return true else return false
  //
  if ( (index<0) || (index>fNtRows))  return kFALSE;
  Int_t outindex = 24*fnRowLow;
  if (index<outindex) {
    sector = index/fnRowLow;
    row    = index - sector*fnRowLow;
    sector++;
    return kTRUE;
  }
  index-= outindex;
  sector = index/fnRowUp;
  row    = index - sector*fnRowUp;
  sector++;
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
  //set radius parameters
  fInnerRadiusLow = kInnerRadiusLow;
  fOuterRadiusLow = kOuterRadiusLow;
  fInnerRadiusUp  = kInnerRadiusUp;
  fOuterRadiusUp  = kOuterRadiusUp; 
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
  //calculate sin and cosine of rotations angle   
  for (Int_t i=1; i<80; i++)
    {
      Float_t angle;
      if(i < 25){
	angle = (i < 13) ? (i-1)*alpha_low : (i-13)*alpha_low;
      }
      else {
	angle = (i < 49) ? (i-25)*alpha_up : (i-49)*alpha_up;
      }
      fRotAngle[i]=TMath::Cos(angle);
      fRotAngle[100+i]=TMath::Sin(angle);
    }
  fbStatus = Update();
}

void  AliTPCParam::AdjustAngles(Int_t isec, Float_t &cos, Float_t &sin) const
{
  //set cosinus and sinus of rotation angles for sector isec
  cos=fRotAngle[isec];
  sin=fRotAngle[100+isec];
}
          
Bool_t AliTPCParam::Update()
{
  fbStatus = kFALSE;
  Int_t i;
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
       Float_t y = (x-0.5*fPadPitchLength)*2.*tan(alpha_up/2)-kEdgeSectorSpace;
       fPadRowUp[i] = x;
       fnPadsUp[i] = (Int_t)(y/fPadPitchWidth) ;        
       if ((fnPadsUp[i]%2) == 0) fnPadsUp[i]-=1;        
    }
  // adjust lower sectors pad row positions and pad numbers 
  for (i = 0;i<fnRowLow;i++) 
    {
       Float_t x  = fInnerRadiusLow +fPadPitchLength*(Float_t)i;
       Float_t y = (x-0.5*fPadPitchLength)*2.*tan(alpha_low/2)-kEdgeSectorSpace;
       fPadRowLow[i] = x;
       fnPadsLow[i] = (Int_t)(y/fPadPitchWidth) ;
       if ((fnPadsLow[i]%2) == 0) fnPadsLow[i]-=1;        
    }

  //that variable are not writen to the file there are calculated
  //
  fWWPitch= fPadPitchLength/Float_t(fnWires);
  fZWidth = fTSample*fDriftV;  
  fNtRows = 24*fnRowLow+48*fnRowUp;

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
      
      R__b >> fInnerRadiusLow;
      R__b >> fInnerRadiusUp;
      R__b >> fOuterRadiusLow;
      R__b >> fOuterRadiusUp;

      R__b >> fPadPitchLength;
      R__b >> fPadPitchWidth;
      R__b >> fPadLength;
      R__b >> fPadWidth;

      R__b >> fnWires;
      
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
      fWWPitch= fPadPitchLength/Float_t(fnWires);
      fZWidth = fTSample*fDriftV;  
      fNtRows = 24*fnRowLow+48*fnRowUp;
      Update();
   } else {
      R__b.WriteVersion(AliTPCParam::IsA());
      TObject::Streamer(R__b);      
      R__b << fInnerRadiusLow;
      R__b << fInnerRadiusUp;
      R__b << fOuterRadiusLow;
      R__b << fOuterRadiusUp;

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

