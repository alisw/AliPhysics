#ifndef TPCParam_H
#define TPCParam_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC parameters          //
////////////////////////////////////////////////

#include "TObject.h"

// the last things from AliTPCSecGeo
//const Float_t z_end = 250.; 
//const Float_t alpha_low=0.523598775; // 30 degrees
//const Float_t alpha_up=0.261799387; //  15 degrees
//const Float_t q_el = 1.602e-19; // elementary charge
//const Float_t adc_sat = 1023; // dynamic range (10 bits)
//const Float_t dyn_range = 2000.; // output dynamic range (mV)



class AliTPCParam : public TObject {
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //ALITPCParam object to be possible change 
  //geometry and some other parameters of TPC   
  //used by AliTPC and AliTPCSector  
public:
  AliTPCParam(); 
  virtual ~AliTPCParam() {;}  //dummy destructor
  
  void XYZtoCRXYZ(Float_t *xyz, 
		  Int_t &sector, Int_t &padrow, Int_t option=3);
  //transform global position to the position relative to the sector padrow
  //if option=0  X calculate absolute            calculate sector
  //if option=1  X           absolute            use input sector
  //if option=2  X           relative to pad row calculate sector
  //if option=3  X           relative            use input sector

  void CRXYZtoXYZ(Float_t *xyz,
	    const Int_t &sector, const Int_t & padrow, Int_t option=3) const;  
  //transform relative position  to the gloabal position

  void CRTimePadtoYZ(Float_t &y, Float_t &z, 
		     const Float_t &time, const Float_t &pad,
		     Int_t sector, Int_t padrow );
  //transform position in digit  units (time slices and pads)  to "normal" 
  //units (cm)   
  void CRYZtoTimePad(const Float_t &y, const Float_t &z, 
		     Float_t &time, Float_t &pad,
  		     Int_t sector, Int_t padrow);
  //transform position in cm to position in digit unit  
  Double_t GetLowMaxY(Int_t irow) const {return irow*0.;}
  Double_t GetUpMaxY(Int_t irow) const {return irow*0;}
  //additional geometrical function
  Int_t GetPadRow(Int_t isec, Float_t &x);
  //return pad row for given sector and position x
  //if res=-1 it is out of sector
  
  Int_t GetWire(Float_t &x);      
  Int_t GetIndex(Int_t sector, Int_t row);  //give index of the given sector and pad row 
  Bool_t   AdjustSectorRow(Int_t index, Int_t & sector, Int_t &row); //return sector and padrow
  //for given index
  Int_t GetNRowsTotal(){return fNtRows;}  //get total nuber of rows
  void SetDefault();          //set defaut TPCparam
  Bool_t Update();            //recalculate and check geometric parameters 
  Bool_t GetStatus();         //get information about object consistency  


  void  AdjustAngles(Int_t isec, Float_t &cos, Float_t &sin) const;
  //set cosinus and sinus of rotation angles for sector isec
  Int_t GetNRowLow() const;   //get the number of pad rows in low sector
  Int_t GetNRowUp() const;    //get the number of pad rows in up sector
  Int_t GetNRow(Int_t isec) {return  ((isec<fNInnerSector) ?  fnRowLow:fnRowUp);}
  //get the nuber of pad row in given sector
  Float_t GetPadRowRadiiLow(Int_t irow) const; //get the pad row (irow) radii
  Float_t GetPadRowRadiiUp(Int_t irow) const;  //get the pad row (irow) radii
  Float_t GetPadRowRadii(Int_t isec,Int_t irow) const {
    return ( (isec < fNInnerSector) ?GetPadRowRadiiLow(irow):GetPadRowRadiiUp(irow));}
    //retrun raii of the pad row irow in sector i
  Int_t GetNPadsLow(Int_t irow) const;    //get the number of pads in row irow 
  Int_t GetNPadsUp(Int_t irow) const;     //get the number of pads in row irow
  Int_t GetNPads(Int_t isector,Int_t irow){
     return ( (isector < fNInnerSector) ?GetNPadsLow(irow) : GetNPadsUp(irow));}
    //get the number of pads  in given sector and row
  //  Int_t GetNPads(Int_t isector, Int_t irow) const;         
   //get the number of pads in sector isector and row irow

  void  SetInnerRadiusLow(Float_t InnerRadiusLow ) { fInnerRadiusLow=InnerRadiusLow;}
  void  SetOuterRadiusLow(Float_t OuterRadiusLow ){  fOuterRadiusLow=OuterRadiusLow;} 
  void  SetInnerRadiusUp(Float_t InnerRadiusUp){  fInnerRadiusUp= InnerRadiusUp;} 
  void  SetOuterRadiusUp(Float_t OuterRadiusUp){  fOuterRadiusUp= OuterRadiusUp;} 
  
  void  SetSectorAngles(Float_t innerangle, Float_t innershift, Float_t outerangle,
			Float_t outershift,Bool_t inDegree=kTRUE);

  void    SetInSecLowEdge(Float_t isle){fInSecLowEdge=isle;}
  void    SetInSecUpEdge(Float_t isue){fInSecUpEdge=isue;}
  void    SetOuSecLowEdge(Float_t osle){fOuSecLowEdge=osle;}
  void    SetOuSecUpEdge(Float_t osue){fOuSecUpEdge=osue;}

  void    SetEdge(Float_t edge){fEdge = edge;}
  void    SetDeadZone(Float_t zone){fDeadZone = zone;} 


  void  SetPadPitchLength(Float_t PadPitchLength){  fPadPitchLength=PadPitchLength;}
  void  SetPadPitchWidth(Float_t PadPitchWidth){  fPadPitchWidth = PadPitchWidth;}
  void  SetPadLength(Float_t PadLength){  fPadLength=PadLength;}
  void  SetPadWidth(Float_t PadWidth) {  fPadWidth=PadWidth;}  
  void  SetDiffT(Float_t DiffT){  fDiffT= DiffT;}
  void  SetDiffL(Float_t DiffL){  fDiffL=DiffL;}
  void  SetDriftV(Float_t DriftV){  fDriftV= DriftV;}
  void  SetOmegaTau(Float_t OmegaTau){  fOmegaTau=OmegaTau;}
  void  SetAttCoef(Float_t AttCoef){  fAttCoef=AttCoef;}
  void  SetOxyCont(Float_t OxyCont){  fOxyCont=OxyCont;}

  void  SetNoise(Float_t Noise ){  fNoise= Noise;}
  void  SetChipGain(Float_t ChipGain){  fChipGain= ChipGain;}
  void  SetGasGain(Float_t GasGain){  fGasGain=GasGain;}
  void  SetTSample(Float_t TSample){  fTSample=TSample;}
  void  SetTSigma(Float_t Sigma){  fTSigma=Sigma;}
  void  SetPadCoupling(Float_t PadCoupling){  fPadCoupling=PadCoupling;}
  void  SetNWires(Int_t nWires){  fnWires=nWires;}
  void  SetWWPitch(Float_t WWPitch){  fWWPitch=WWPitch;}
  void  SetZeroSup(Int_t ZeroSup){  fZeroSup=ZeroSup;}

  Float_t  GetInnerRadiusLow(){return fInnerRadiusLow;}
  Float_t  GetOuterRadiusLow(){return fOuterRadiusLow;} 
  Float_t  GetInnerRadiusUp(){return fInnerRadiusUp;} 
  Float_t  GetOuterRadiusUp(){return fOuterRadiusUp;} 

  Float_t  GetInnerAngle(){return fInnerAngle;}
  Float_t  GetInnerAngleShift(){return fInnerAngleShift;}
  Float_t  GetOuterAngle(){return fOuterAngle;} 
  Float_t  GetOuterAngleShift(){return fOuterAngleShift;} 
  Int_t    GetNInnerSector(){return fNInnerSector;}
  Int_t    GetNOuterSector(){return fNOuterSector;}
  Int_t    GetNSector(){return fNSector;}

  Float_t    GetInSecLowEdge(){return fInSecLowEdge;}
  Float_t    GetInSecUpEdge(){return fInSecUpEdge;}
  Float_t    GetOuSecLowEdge(){return fOuSecLowEdge;}
  Float_t    GetOuSecUpEdge(){return fOuSecUpEdge;}  

  Float_t    GetEdge(){return fEdge;}
  Float_t    GetDeadZone(){return fDeadZone;}

  Float_t  GetPadPitchLength(){return fPadPitchLength;}
  Float_t  GetPadPitchWidth(){return fPadPitchWidth;}
  Float_t  GetPadLength(){return fPadLength;}
  Float_t  GetPadWidth() {return fPadWidth;}  
  Float_t  GetDiffT(){return fDiffT;}
  Float_t  GetDiffL(){return fDiffL;}
  Float_t  GetDriftV(){return fDriftV;}
  Float_t  GetOmegaTau(){return fOmegaTau;}
  Float_t  GetAttCoef(){return fAttCoef;}
  Float_t  GetOxyCont(){return fOxyCont;}

  Float_t  GetNoise(){return fNoise;}
  Float_t  GetChipGain(){return fChipGain;}
  Float_t  GetGasGain(){return fGasGain;}
  Float_t  GetTSample(){return fTSample;}
  Float_t  GetTSigma(){return fTSigma;}
  Float_t  GetZWidth(){return fZWidth;}
  Float_t  GetZSigma(){return fTSigma*fDriftV;}  
  Float_t  GetPadCoupling(){return fPadCoupling;}
  Int_t    GetNWires(){return fnWires;}
  Float_t  GetWWPitch(){return fWWPitch;}
  Int_t    GetZeroSup(){return fZeroSup;}
  Int_t    GetMaxTBin(){return fMaxTBin;}

private :
  Bool_t fbStatus;  //indicates consistency of the data
  //---------------------------------------------------------------------
  //   ALICE TPC sector geometry
  //--------------------------------------------------------------------
  
  Float_t fInnerRadiusLow;  // lower radius of inner sector
  Float_t fOuterRadiusLow;  // lower radius of outer sector
  Float_t fInnerRadiusUp;   // upper radius of inner  sector
  Float_t fOuterRadiusUp;   // upper radius of outer  sector

  Float_t fInnerAngle;       //opening angle of Inner sector
  Float_t fInnerAngleShift;  //shift of first inner sector center to the 0
  Float_t fOuterAngle;       //opening angle of outer sector
  Float_t fOuterAngleShift;  //shift of first sector center to the 0
 
  Int_t   fNInnerSector;      //!number of inner sectors
  Int_t   fNOuterSector;      //!number of outer sectors
  Int_t   fNSector;           //! total number of sectors

  Float_t   fInSecLowEdge;      // inner sector lower edge
  Float_t   fInSecUpEdge;       // inner sector upper edge
  Float_t   fOuSecLowEdge;      // outer sector lower edge
  Float_t   fOuSecUpEdge;       // outer sector upper edge
  
  Float_t   fEdge;              // thickness of the sector edge
  Float_t   fDeadZone;          // dead zone due to the sector mounting etc.
 
  //---------------------------------------------------------------------
  //   ALICE TPC pad parameters
  //--------------------------------------------------------------------
  Float_t   fPadPitchLength;    //pad pitch length
  Float_t   fPadPitchWidth;     //pad pitch width
  Float_t   fPadLength;         //pad  length
  Float_t   fPadWidth;          //pad  width
  
  
  Int_t fnRowLow;           //  number of pad rows per low sector 
  Int_t fnRowUp;            //  number of pad rows per sector up 
  Float_t  fPadRowLow[600]; // Lower sector, pad row radii
  Float_t  fPadRowUp[600];  // Upper sector, pad row radii
  Int_t    fnPadsLow[600];     // Lower sector, number of pads per row
  Int_t    fnPadsUp[600];      //  Upper sector, number of pads per row
  Float_t fRotAngle[200];      //  sin and cos of rotation angles for 
                                 //  diferent sectors

  Int_t fnWires;            //  Number of wires per pad
  Float_t fWWPitch;         // pitch between wires   
  //---------------------------------------------------------------------
  //   ALICE TPC Gas Parameters
  //--------------------------------------------------------------------
  Float_t  fDiffT;          //tangencial diffusion constant
  Float_t  fDiffL;          //longutudinal diffusion constant
  Float_t  fGasGain;        //gas gain constant
  Float_t  fDriftV;          //drift velocity constant
  Float_t  fOmegaTau;       //omega tau ExB coeficient
  Float_t  fAttCoef;        //attachment coefitients
  Float_t  fOxyCont;        //oxygen content
  //---------------------------------------------------------------------
  //   ALICE TPC  Electronics Parameters
  //--------------------------------------------------------------------
  Float_t fPadCoupling;     //coupling factor ration of  anode signal 
                            //and total pads signal  
  Int_t fZeroSup;         //zero suppresion constant
  Float_t fNoise;         //noise sigma constant
  Float_t fChipGain;      //preamp shaper constant
  Float_t fTSample; // sampling time
  Float_t fZWidth;  //derived value calculated using TSample and driftw 
  Float_t fTSigma;  // width of the Preamp/Shaper function
  Int_t   fMaxTBin; //maximum time bin number
  //--------------------------------------------------------
  //
  Int_t fNtRows;  //total number of rows in TPC  
  ClassDef(AliTPCParam,2)  //parameter  object for set:TPC
};


/////////////////////////////////////////////////////////////////////////////
//
//---------------------------------------------------------------------
//   ALICE TPC Cluster Parameters
//--------------------------------------------------------------------
//
//
// Sigma rphi
/*const Float_t a_rphi=0.41818e-2;
const Float_t b_rphi=0.17460e-4;
const Float_t c_rphi=0.30993e-2;
const Float_t d_rphi=0.41061e-3;
// Sigma z
const Float_t a_z=0.39614e-2;
const Float_t b_z=0.22443e-4;
const Float_t c_z=0.51504e-1;
// Cluster width in rphi
const Float_t ac_rphi=0.18322;
const Float_t bc_rphi=0.59551e-3;
const Float_t cc_rphi=0.60952e-1;
// Cluster width in z
const Float_t ac_z=0.19081;
const Float_t bc_z=0.55938e-3;
const Float_t cc_z=0.30428;
*/



#endif  
