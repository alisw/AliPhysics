#ifndef ALITPCPARAM_H
#define ALITPCPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC parameters          //
////////////////////////////////////////////////

#include "AliDetectorParam.h"
#include "TMath.h"



class AliTPCParam : public AliDetectorParam {
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //ALITPCParam object to be possible change 
  //geometry and some other parameters of TPC   
  //used by AliTPC and AliTPCSector 
 
public:
  AliTPCParam(); 
  virtual ~AliTPCParam();
  
  virtual Bool_t  Transform(Float_t *xyz, Int_t *index, Int_t* oindex);
  //transformation from input coodination system to output coordination system  
  Int_t  Transform0to1(Float_t *xyz, Int_t *index) const;
  //trasforamtion from global to global - adjust index[0] sector 
  //return value is equal to sector corresponding to global position
  void Transform1to2(Float_t *xyz, Int_t *index) const;
  //transformation to rotated coordinata 
  void Transform2to1(Float_t *xyz, Int_t *index) const;
  //transformation from rotated coordinata to global coordinata
  void Transform2to2(Float_t *xyz, Int_t *index, Int_t *oindex) const;
  //transform rotated coordinata of one sector to rotated
  //coordinata relative to another sector
  Float_t  Transform2to2NearestWire(Float_t *xyz, Int_t *index) const;
  //round x position to nearest wire
  Int_t   Transform2to3(Float_t *xyz, Int_t *index) const;
  //calulate coresponding index[2] -pad row for straight rows
  //does not change xyz[] 
  //return pad - row 
  void   Transform3to4(Float_t *xyz, Int_t *index) const;
  //valid only for straight rows straight rows
  //calculate xyz[0] position relative to given index
  //return pad - row 
  void   Transform4to3(Float_t *xyz, Int_t *index) const;
  //valid only for straight rows straight rows
  //transform  xyz[0] position relative to given index
  void   Transform2to5( Float_t *xyz, Int_t *index) const;
  //transform [x,y,z] to [r,rphi,z]
  void   Transform5to2(Float_t *xyz, Int_t *index) const;
  //transform [r,rphi,z] coordinata to [x,y,z] 
  void  Transform4to8(Float_t *xyz, Int_t *index) const;
  //transform xyz coordinata to 'digit' coordinata
  void  Transform8to4(Float_t *xyz, Int_t *index) const;
  //transform  'digit' coordinata to xyz coordinata   
  void  Transform6to8(Float_t *xyz, Int_t *index) const;
  //transform dr,f coordinata to 'digit' coordinata
  void  Transform8to6(Float_t *xyz, Int_t *index) const;
  //transform 'digit' coordinata to dr,f coordinata 

  virtual Int_t  Transform2toPadRow(Float_t *xyz, Int_t *index) const{return 0;}
  //transform rotated to

  virtual  Int_t GetPadRow(Float_t *xyz, Int_t *index) const ;
  //return pad row of point xyz - xyz is given in coordinate system -(given by index)
  //output system is 3 for straight row and 7 for cylindrical row
  virtual void XYZtoCRXYZ(Float_t *xyz, 
		      Int_t &sector, Int_t &padrow, Int_t option=3) const {;}
  //transform global position to the position relative to the sector padrow
  //if option=0  X calculate absolute            calculate sector
  //if option=1  X           absolute            use input sector
  //if option=2  X           relative to pad row calculate sector
  //if option=3  X           relative            use input sector

  virtual void CRXYZtoXYZ(Float_t *xyz,
			  const Int_t &sector, const Int_t & padrow, Int_t option=3) const {;}  
  //transform relative position  to the gloabal position

  virtual void CRTimePadtoYZ(Float_t &y, Float_t &z, 
			     const Float_t &time, const Float_t &pad,
			     Int_t sector, Int_t padrow ){;}
  //transform position in digit  units (time slices and pads)  to "normal" 
  //units (cm)   
  virtual void CRYZtoTimePad(const Float_t &y, const Float_t &z, 
			     Float_t &time, Float_t &pad,
			     Int_t sector, Int_t padrow){;}
  //transform position in cm to position in digit unit 
  virtual Int_t   CalcResponse(Float_t* x, Int_t * index, Int_t row){return 0;}
  //calculate bin response as function of the input position -x and the weight 
  //if row -pad row is equal -1 calculate response for each pad row 
  //otherwise it calculate only in given pad row
  //return number of valid response bin
  virtual void SetDefault();          //set defaut TPCparam 
  virtual Bool_t Update();            //recalculate and check geometric parameters 
  Bool_t GetStatus();         //get information about object consistency  
  Int_t GetIndex(Int_t sector, Int_t row);  //give index of the given sector and pad row 
  Int_t GetNSegmentsTotal() const {return fNtRows;} 
  Double_t GetLowMaxY(Int_t irow) const {return irow*0.;}
  Double_t GetUpMaxY(Int_t irow) const {return irow*0;}
  //additional geometrical function - for Belikov
 
  Bool_t   AdjustSectorRow(Int_t index, Int_t & sector, Int_t &row) const; //return sector and padrow
  //for given index

  void  AdjustCosSin(Int_t isec, Float_t &cos, Float_t &sin) const;
  //set cosinus and sinus of rotation angles for sector isec 
  Float_t GetAngle(Int_t isec) const;
  //
  //set sector parameters
  //
  void  SetInnerRadiusLow(Float_t InnerRadiusLow )  { fInnerRadiusLow=InnerRadiusLow;}
  void  SetOuterRadiusLow(Float_t OuterRadiusLow )  { fOuterRadiusLow=OuterRadiusLow;} 
  void  SetInnerRadiusUp(Float_t InnerRadiusUp)  {  fInnerRadiusUp= InnerRadiusUp;} 
  void  SetOuterRadiusUp(Float_t OuterRadiusUp) {  fOuterRadiusUp= OuterRadiusUp;}   
  void  SetSectorAngles(Float_t innerangle, Float_t innershift, Float_t outerangle,
			Float_t outershift);
  void  SetInnerFrameSpace(Float_t frspace) {fInnerFrameSpace = frspace;}
  void  SetOuterFrameSpace(Float_t frspace) {fOuterFrameSpace = frspace;}
  void  SetInnerWireMount(Float_t fmount) {fInnerWireMount = fmount;}
  void  SetOuterWireMount(Float_t fmount) {fOuterWireMount = fmount;}
  void  SetZLength(Float_t zlength) {fZLength = zlength;} 
  void  SetGeometryType(Int_t type) {fGeometryType = type;}
  //
  //set wire parameters
  //
  void  SetInnerNWires(Int_t nWires){  fNInnerWiresPerPad=nWires;}
  void  SetInnerDummyWire(Int_t dummy) {fInnerDummyWire  = dummy;}
  void  SetInnerOffWire(Float_t offset) {fInnerOffWire =offset;}    
  void  SetOuterNWires(Int_t nWires){  fNOuterWiresPerPad=nWires;}
  void  SetOuterDummyWire(Int_t dummy) {fOuterDummyWire  = dummy;}
  void  SetOuterOffWire(Float_t offset) {fOuterOffWire =offset;}    
  //
  //set pad parameter
  //
  void  SetInnerPadPitchLength(Float_t PadPitchLength){  fInnerPadPitchLength=PadPitchLength;}
  void  SetInnerPadPitchWidth(Float_t PadPitchWidth){  fInnerPadPitchWidth = PadPitchWidth;}
  void  SetInnerPadLength(Float_t PadLength){  fInnerPadLength=PadLength;}
  void  SetInnerPadWidth(Float_t PadWidth) {  fInnerPadWidth=PadWidth;}  
  void  SetOuterPadPitchLength(Float_t PadPitchLength){  fOuterPadPitchLength=PadPitchLength;}
  void  SetOuterPadPitchWidth(Float_t PadPitchWidth){  fOuterPadPitchWidth = PadPitchWidth;}
  void  SetOuterPadLength(Float_t PadLength){  fOuterPadLength=PadLength;}
  void  SetOuterPadWidth(Float_t PadWidth) {  fOuterPadWidth=PadWidth;} 
  void  SetMWPCReadout(Bool_t type) {fBMWPCReadout = type;}
  void  SetNCrossRows(Int_t rows){fNCrossRows = rows;}
  //
  //set gas paremeters
  //
  void  SetDiffT(Float_t DiffT){  fDiffT= DiffT;}
  void  SetDiffL(Float_t DiffL){  fDiffL=DiffL;}
  void  SetGasGain(Float_t GasGain){  fGasGain=GasGain;}
  void  SetDriftV(Float_t DriftV){  fDriftV= DriftV;}
  void  SetOmegaTau(Float_t OmegaTau){  fOmegaTau=OmegaTau;}
  void  SetAttCoef(Float_t AttCoef){  fAttCoef=AttCoef;}
  void  SetOxyCont(Float_t OxyCont){  fOxyCont=OxyCont;}
  //
  //set electronivc parameters  
  //
  void  SetPadCoupling(Float_t PadCoupling){  fPadCoupling=PadCoupling;}
  void  SetZeroSup(Int_t ZeroSup)    {  fZeroSup=ZeroSup;}
  void  SetNoise(Float_t Noise )     {  fNoise= Noise;}
  void  SetChipGain(Float_t ChipGain){  fChipGain= ChipGain;}
  void  SetChipNorm(Float_t ChipNorm){  fChipNorm= ChipNorm;}
  void  SetTSample(Float_t TSample)  {  fTSample=TSample;}
  void  SetTFWHM(Float_t fwhm)     {  fTSigma=fwhm/2.35;}
  void  SetMaxTBin(Int_t maxtbin)  {  fMaxTBin = maxtbin;}
  void  SetADCSat(Int_t adcsat)    {  fADCSat  = adcsat;}
  void  SetADCDynRange(Float_t adcdynrange) {fADCDynRange = adcdynrange;}
  //
  //set response  parameters  
  //
  void  SetNResponseMax(Int_t max) { fNResponseMax = max;} 
  void  SetResponseThreshold(Int_t threshold) {fResponseThreshold = threshold;}
  //
  //get sector parameters
  //
  Float_t  GetInnerRadiusLow() const {return fInnerRadiusLow;}
  Float_t  GetInnerRadiusUp() const {return fInnerRadiusUp;} 
  Float_t  GetOuterRadiusLow() const {return fOuterRadiusLow;} 
  Float_t  GetOuterRadiusUp() const {return fOuterRadiusUp;} 
  Float_t  GetInnerFrameSpace() const {return fInnerFrameSpace;}
  Float_t  GetOuterFrameSpace() const {return fOuterFrameSpace;}
  Float_t  GetInnerWireMount() const {return fInnerWireMount;}
  Float_t  GetOuterWireMount() const {return fOuterWireMount;}
  Float_t  GetInnerAngle() const ;
  Float_t  GetInnerAngleShift() const ;
  Float_t  GetOuterAngle() const ;
  Float_t  GetOuterAngleShift() const ; 
  Int_t    GetNInnerSector() const {return fNInnerSector;}
  Int_t    GetNOuterSector() const {return fNOuterSector;}
  Int_t    GetNSector() const {return fNSector;}
  Float_t  GetZLength() const {return fZLength;}
  Int_t    GetGeometryType() const {return fGeometryType;}

  //
  //get wires parameter
  //
  Int_t    GetInnerNWires() const {return fNInnerWiresPerPad;}
  Float_t  GetInnerWWPitch() const {return fInnerWWPitch;}  
  Int_t    GetInnerDummyWire() const {return fInnerDummyWire;}
  Float_t  GetInnerOffWire() const {return fInnerOffWire;}
  Float_t  GetRInnerFirstWire() const {return fRInnerFirstWire;}
  Float_t  GetRInnerLastWire() const {return fRInnerLastWire;}
  Int_t    GetOuterNWires() const {return fNOuterWiresPerPad;}
  Float_t  GetOuterWWPitch() const {return fOuterWWPitch;}  
  Int_t    GetOuterDummyWire() const {return fOuterDummyWire;}
  Float_t  GetOuterOffWire() const {return fOuterOffWire;}
  Float_t  GetROuterFirstWire() const {return fROuterFirstWire;}
  Float_t  GetROuterLastWire() const {return fROuterLastWire;}  
  Float_t  GetWWPitch(Int_t isector = 0) const  {
    return ( (isector < fNInnerSector) ? fInnerWWPitch :fOuterWWPitch);} 
  //
  //get pad  parameters
  //
  Float_t  GetInnerPadPitchLength() const {return fInnerPadPitchLength;}
  Float_t  GetInnerPadPitchWidth() const {return fInnerPadPitchWidth;}
  Float_t  GetInnerPadLength() const {return fInnerPadLength;}
  Float_t  GetInnerPadWidth() const  {return fInnerPadWidth;}  
  Float_t  GetOuterPadPitchLength() const {return fOuterPadPitchLength;}
  Float_t  GetOuterPadPitchWidth() const {return fOuterPadPitchWidth;}
  Float_t  GetOuterPadLength() const {return fOuterPadLength;}
  Float_t  GetOuterPadWidth()  const {return fOuterPadWidth;}  
  Bool_t   GetMWPCReadout() const {return fBMWPCReadout;}
  Int_t    GetNCrossRows() const {return fNCrossRows;}
  Float_t  GetPadPitchWidth(Int_t isector = 0) const  {
    return ( (isector < fNInnerSector) ? fInnerPadPitchWidth :fOuterPadPitchWidth);}
  Float_t  GetPadPitchLength(Int_t isector = 0)  const {
    return ( (isector < fNInnerSector) ? fInnerPadPitchLength :fOuterPadPitchLength);} 

  Int_t GetNRowLow() const;   //get the number of pad rows in low sector
  Int_t GetNRowUp() const;    //get the number of pad rows in up sector
  Int_t GetNRow(Int_t isec) const {return  ((isec<fNInnerSector) ?  fNRowLow:fNRowUp);}
  Int_t GetNRowsTotal(){return fNtRows;}  //get total nuber of rows
  Float_t GetPadRowRadiiLow(Int_t irow) const; //get the pad row (irow) radii
  Float_t GetPadRowRadiiUp(Int_t irow) const;  //get the pad row (irow) radii
  Float_t GetPadRowRadii(Int_t isec,Int_t irow) const {
    return ( (isec < fNInnerSector) ?GetPadRowRadiiLow(irow):GetPadRowRadiiUp(irow));}
    //retrun radii of the pad row irow in sector i
  Int_t GetNPadsLow(Int_t irow) const;    //get the number of pads in row irow 
  Int_t GetNPadsUp(Int_t irow) const;     //get the number of pads in row irow
  Int_t GetNPads(Int_t isector,Int_t irow) const{
     return ( (isector < fNInnerSector) ?GetNPadsLow(irow) : GetNPadsUp(irow));}    
  //
  //get GAS parameters 
  //
  Float_t  GetDiffT() const {return fDiffT;}
  Float_t  GetDiffL() const {return fDiffL;}
  Float_t  GetGasGain() const {return fGasGain;}
  Float_t  GetDriftV() const {return fDriftV;}
  Float_t  GetOmegaTau() const {return fOmegaTau;}
  Float_t  GetAttCoef() const {return fAttCoef;}
  Float_t  GetOxyCont() const {return fOxyCont;}
  //
  //get Electronic parameters
  //
  Float_t  GetPadCoupling() const {return fPadCoupling;}
  Int_t    GetZeroSup() const {return fZeroSup;}
  Float_t  GetNoise() const {return fNoise;}
  Float_t  GetChipGain() const {return fChipGain;}
  Float_t  GetChipNorm() const {return fChipNorm;}
  Float_t  GetTSample() const {return fTSample;}
  Float_t  GetZWidth() const {return fZWidth;}
  Float_t  GetTFWHM() const {return fTSigma*2.35;}
  Float_t  GetZSigma() const {return fTSigma*fDriftV;}  
  virtual  Float_t  GetZOffset() {return 3*fTSigma*fDriftV;}
  Int_t    GetMaxTBin() const {return fMaxTBin;}
  Int_t    GetADCSat() const {return fADCSat;}
  Float_t  GetADCDynRange() const {return fADCDynRange;}
  Float_t  GetTotalNormFac() const {return fTotalNormFac;}
  Float_t  GetNoiseNormFac() const {return fNoiseNormFac;}
  //
  // get response data
  //  
  Int_t * GetResBin(Int_t i);  
  //return response bin i  - bin given by  padrow [0] pad[1] timebin[2] 
  Float_t GetResWeight(Int_t i);
  //return  weight of response bin i
protected :

  Bool_t fbStatus;  //indicates consistency of the data
  //---------------------------------------------------------------------
  //   ALICE TPC sector geometry
  //--------------------------------------------------------------------  
  Float_t fInnerRadiusLow;    // lower radius of inner sector-IP
  Float_t fInnerRadiusUp;     // upper radius of inner  sector-IP
  Float_t fOuterRadiusUp;     // upper radius of outer  sector-IP
  Float_t fOuterRadiusLow;    // lower radius of outer sector-IP
  Float_t fInnerAngle;        //opening angle of Inner sector
  Float_t fInnerAngleShift;   //shift of first inner sector center to the 0
  Float_t fOuterAngle;        //opening angle of outer sector
  Float_t fOuterAngleShift;   //shift of first sector center to the 0  
  Float_t fInnerFrameSpace;   //space for inner frame in the phi direction 
  Float_t fOuterFrameSpace;   //space for outer frame in the phi direction 
  Float_t fInnerWireMount;    //space for wire mount, inner sector
  Float_t fOuterWireMount;    //space for wire mount, outer sector
  Int_t   fNInnerSector;      //!number of inner sectors             -calculated
  Int_t   fNOuterSector;      //!number of outer sectors             -calculated
  Int_t   fNSector;           //! total number of sectors            -calculated
  Float_t fZLength;           //length of the drift region of the TPC
  Float_t *fRotAngle;         //  sin and cos of rotation angles for 
                              //  diferent sectors -calculated
  Int_t   fGeometryType;      //type of geometry -0 straight rows
  //1-cylindrical
  //---------------------------------------------------------------------
  //   ALICE TPC wires  geometry - for GEM we can consider that it is gating  
  //--------------------------------------------------------------------
  Int_t   fNInnerWiresPerPad; //  Number of wires per pad
  Float_t fInnerWWPitch;      // pitch between wires  in inner sector     - calculated
  Int_t   fInnerDummyWire;    //number of wires without pad readout
  Float_t fInnerOffWire;      //oofset of first wire to the begining of the sector
  Float_t fRInnerFirstWire;   //position of the first wire                -calculated
  Float_t fRInnerLastWire;    //position of the last wire                 -calculated
  Int_t   fNOuterWiresPerPad; //  Number of wires per pad
  Float_t fOuterWWPitch;      // pitch between wires in outer sector      -calculated
  Int_t   fOuterDummyWire;    //number of wires without pad readout
  Float_t fOuterOffWire;      //oofset of first wire to the begining of the sector
  Float_t fROuterFirstWire;   //position of the first wire                -calulated
  Float_t fROuterLastWire;    //position of the last wire                 -calculated 
  //---------------------------------------------------------------------
  //   ALICE TPC pad parameters
  //--------------------------------------------------------------------
  Float_t   fInnerPadPitchLength;    //Inner pad pitch length
  Float_t   fInnerPadPitchWidth;     //Inner pad pitch width
  Float_t   fInnerPadLength;         //Inner pad  length
  Float_t   fInnerPadWidth;          //Inner pad  width
  Float_t   fOuterPadPitchLength;    //Outer pad pitch length
  Float_t   fOuterPadPitchWidth;     //Outer pad pitch width
  Float_t   fOuterPadLength;         //Outer pad  length
  Float_t   fOuterPadWidth;          //Outer pad  width
  Bool_t    fBMWPCReadout;           //indicate wire readout - kTRUE or GEM readout -kFALSE
  Int_t     fNCrossRows;             //number of rows to crostalk calculation
      
  Int_t fNRowLow;           //number of pad rows per low sector    -calculated
  Int_t fNRowUp;            //number of pad rows per sector up     -calculated
  Int_t fNtRows;            //total number of rows in TPC          -calculated
  Float_t  fPadRowLow[600]; //Lower sector, pad row radii          -calculated
  Float_t  fPadRowUp[600];  //Upper sector, pad row radii          -calculated 
  Int_t    fNPadsLow[600];  //Lower sector, number of pads per row -calculated
  Int_t    fNPadsUp[600];   //Upper sector, number of pads per row -calculated  
  //---------------------------------------------------------------------
  //   ALICE TPC Gas Parameters
  //--------------------------------------------------------------------
  Float_t  fDiffT;          //tangencial diffusion constant
  Float_t  fDiffL;          //longutudinal diffusion constant
  Float_t  fGasGain;        //gas gain constant
  Float_t  fDriftV;         //drift velocity constant
  Float_t  fOmegaTau;       //omega tau ExB coeficient
  Float_t  fAttCoef;        //attachment coefitients
  Float_t  fOxyCont;        //oxygen content
  //---------------------------------------------------------------------
  //   ALICE TPC  Electronics Parameters
  //--------------------------------------------------------------------
  Float_t fPadCoupling;     //coupling factor ration of  anode signal 
                            //and total pads signal  
  Int_t fZeroSup;           //zero suppresion constant
  Float_t fNoise;           //noise sigma constant
  Float_t fChipGain;        //preamp shaper constant
  Float_t fChipNorm;         //preamp shaper normalisation 
  Float_t fTSample;         //sampling time
  Float_t fZWidth;          //derived value calculated using TSample and driftw  -computed
  Float_t fTSigma;          //width of the Preamp/Shaper function
  Int_t   fMaxTBin;         //maximum time bin number
  Int_t   fADCSat;          //saturation value of ADC (10 bits)
  Float_t fADCDynRange;     //input dynamic range (mV)
  Float_t fTotalNormFac;    //full normalisation factor - calculated
  Float_t fNoiseNormFac;    //normalisation factor to transform noise in electron to ADC channel   
  
 
protected:
  //---------------------------------------------------------------------
  // ALICE TPC response data 
  //---------------------------------------------------------------------
  Int_t   fNResponseMax;   //maximal dimension of response        
  Float_t fResponseThreshold; //threshold for accepted response   
  Int_t   fCurrentMax;     //current maximal dimension            -calulated 
  Int_t   *fResponseBin;    //array with bins                     -calulated
  Float_t *fResponseWeight; //array with response                 -calulated

  ClassDef(AliTPCParam,2)  //parameter  object for set:TPC
};

 
inline Int_t * AliTPCParam::GetResBin(Int_t i)
{
  //return response bin i  - bin given by  padrow [0] pad[1] timebin[2] 
  if (i<fCurrentMax) return &fResponseBin[i*3];
  else return 0;
};
  
inline Float_t AliTPCParam::GetResWeight(Int_t i)
{
  //return  weight of response bin i
  if (i<fCurrentMax) return fResponseWeight[i];
  else return 0;
}


inline void  AliTPCParam::AdjustCosSin(Int_t isec, Float_t &cos, Float_t &sin) const
{
  //
  //set cosinus and sinus of rotation angles for sector isec
  //
  cos=fRotAngle[isec*4];
  sin=fRotAngle[isec*4+1];
}

inline Float_t   AliTPCParam::GetAngle(Int_t isec) const
{
  //
  //return rotation angle of given sector
  //
  return fRotAngle[isec*4+2];
}


inline void AliTPCParam::Transform1to2(Float_t *xyz, Int_t *index) const
{
  //transformation to rotated coordinates 
  //we must have information about sector!

  //rotate to given sector
  Float_t cos,sin;
  AdjustCosSin(index[1],cos,sin);   
  Float_t x1=xyz[0]*cos + xyz[1]*sin;
  Float_t y1=-xyz[0]*sin + xyz[1]*cos; 
  xyz[0]=x1;
  xyz[1]=y1;
  xyz[2]=fZLength-TMath::Abs(xyz[2]); 
  index[0]=2;
}

inline void AliTPCParam::Transform2to1(Float_t *xyz, Int_t *index) const
{
  //
  //transformation from  rotated coordinates to global coordinates
  //
  Float_t cos,sin;
  AdjustCosSin(index[1],cos,sin);   
  Float_t x1=xyz[0]*cos - xyz[1]*sin;
  Float_t y1=xyz[0]*sin + xyz[1]*cos; 
  xyz[0]=x1;
  xyz[1]=y1;
  xyz[2]=fZLength-xyz[2]; 
  if (index[1]<fNInnerSector)
    if ( index[1]>=(fNInnerSector>>1))	xyz[2]*=-1.;
  else 
    if ( (index[1]-fNInnerSector) > (fNOuterSector>>1) )    xyz[2]*=-1;      
  index[0]=1;
}

inline void AliTPCParam::Transform2to2(Float_t *xyz, Int_t *index, Int_t *oindex) const
{
  //transform rotated coordinats of one sector to rotated
  //coordinates relative to another sector
  Transform2to1(xyz,index);
  Transform1to2(xyz,oindex);
  index[0]=2;
  index[1]=oindex[1];  
}

inline Float_t  AliTPCParam::Transform2to2NearestWire(Float_t *xyz, Int_t *index)  const
{
  //
  // asigns the x-position of the closest wire to xyz[0], return the
  // electron to closest wire distance
  //
  Float_t xnew,dx;
  if (index[1]<fNInnerSector) {
     xnew = fRInnerFirstWire+TMath::Nint((xyz[0]-fRInnerFirstWire)/fInnerWWPitch)*fInnerWWPitch;
    }
    else {
     xnew = fROuterFirstWire+TMath::Nint((xyz[0]-fROuterFirstWire)/fOuterWWPitch)*fOuterWWPitch;
    }
  dx = xnew-xyz[0];
  xyz[0]=xnew;
  return  dx;
}

inline Int_t   AliTPCParam::Transform2to3(Float_t *xyz, Int_t *index)  const
{
  //
  //calulates coresponding pad row number, sets index[2] for straight rows
  //does not change xyz[] information
  //valid only for straight row
  //
  if  (index[1]<fNInnerSector)   
    index[2] =TMath::Nint((xyz[0]-fPadRowLow[0])/fInnerPadPitchLength);
  else
    index[2] = TMath::Nint((xyz[0]-fPadRowUp[0])/fOuterPadPitchLength);
  index[0]=3;
  return index[2];
}

inline void   AliTPCParam::Transform3to4(Float_t *xyz, Int_t *index)  const
{
  //
  //valid only for straight rows straight rows
  //calculate xyz[0] position relative to given index
  //
  if  (index[1]<fNInnerSector)   
    xyz[0] -=index[2]*fInnerPadPitchLength+fPadRowLow[0];
  else
    xyz[0] -=index[2]*fOuterPadPitchLength+fPadRowUp[0];
  index[0]  =4;
}

inline void   AliTPCParam::Transform4to3(Float_t *xyz, Int_t *index) const
{
  //
  //valid only for straight rows 
  //transforms  relative xyz[0] to the global one within given sector
  //
  if  (index[1]<fNInnerSector)   
    xyz[0] +=index[2]*fInnerPadPitchLength+fPadRowLow[0];
  else
    xyz[0] +=index[2]*fOuterPadPitchLength+fPadRowUp[0];
  index[0]  =3;
}


inline void   AliTPCParam::Transform2to5( Float_t *xyz, Int_t *index) const
{
  //
  //transform [x,y,z] to [r,phi,z]
  //
  Float_t angle;
  Float_t r = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  if ((xyz[0]==0)&&(xyz[1]==0)) angle = 0;
  else
    {
      angle =TMath::ASin(xyz[1]/r);
      if   (xyz[0]<0)   angle=TMath::Pi()-angle;
      if ( (xyz[0]>0) && (xyz[1]<0) ) angle=2*TMath::Pi()+angle;
    }
  xyz[0]=r;
  xyz[1]=angle;
  index[0]=5;
}

inline void   AliTPCParam::Transform5to2( Float_t *xyz, Int_t *index)  const
{
  //
  //transform [r,rphi,z] to [x,y,z] 
  //
  Float_t r = xyz[0];
  Float_t angle= xyz[1];
  xyz[0]=r*TMath::Cos(angle);
  xyz[1]=r*TMath::Sin(angle);
  index[0]=2;
}

inline void AliTPCParam::Transform4to8(Float_t *xyz, Int_t *index) const
{
  //
  //transform xyz coordinates to 'digit' coordinates
  //
  if  (index[1]<fNInnerSector) {    
    xyz[0]/=fInnerPadPitchLength;
    xyz[1]/=fInnerPadPitchWidth;
    xyz[2]/=fZWidth;
  }
  else{    
    xyz[0]/=fOuterPadPitchLength;
    xyz[1]/=fOuterPadPitchWidth;
    xyz[2]/=fZWidth;
  }        
  index[0]=8;
}

inline void AliTPCParam::Transform8to4(Float_t *xyz, Int_t *index) const
{
  //
  //transforms 'digit' coordinates to xyz coordinates
  //
  if  (index[1]<fNInnerSector) {    
    xyz[0]*=fInnerPadPitchLength;
    xyz[1]*=fInnerPadPitchWidth;
    xyz[2]*=fZWidth;
  }
  else{    
    xyz[0]*=fOuterPadPitchLength;
    xyz[1]*=fOuterPadPitchWidth;
    xyz[2]*=fZWidth;
  }        
  index[0]=4;
}

inline void  AliTPCParam::Transform6to8(Float_t *xyz, Int_t *index) const
{
  //
  //transforms cylindrical xyz coordinates to 'digit' coordinates
  //
  
  if  (index[1]<fNInnerSector) {    
    xyz[0]/=fInnerPadPitchLength;
    xyz[1]*=xyz[0]/fInnerPadPitchWidth;
    xyz[2]/=fZWidth;
  }
  else{   
    xyz[0]/=fOuterPadPitchLength;
    xyz[1]*=xyz[0]/fOuterPadPitchWidth;
    xyz[2]/=fZWidth;
  }        
  index[0]=8;
}

inline void  AliTPCParam::Transform8to6(Float_t *xyz, Int_t *index) const
{
  //
  //transforms 'digit' coordinates to cylindrical xyz coordinates 
  //
  
  if  (index[1]<fNInnerSector) {    
    xyz[0]*=fInnerPadPitchLength;
    xyz[1]/=xyz[0]/fInnerPadPitchWidth;
    xyz[2]*=fZWidth;
  }
  else{   
    xyz[0]*=fOuterPadPitchLength;
    xyz[1]/=xyz[0]/fOuterPadPitchWidth;
    xyz[2]*=fZWidth;
  }        
  index[0]=6;
}



#endif  
