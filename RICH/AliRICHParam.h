#ifndef AliRICHParam_h
#define AliRICHParam_h

#include <TObject.h>
#include "AliRICHConst.h"

class AliRICHParam :public TObject  
{
public:
           AliRICHParam();  
  virtual ~AliRICHParam()                    {;}  
  
  inline  Int_t Neighbours(Int_t iPadX,Int_t iPadY,Int_t aListX[4],Int_t aListY[4])const;                      //pad->neibours
  inline  void   SigGenInit(Float_t x,Float_t y);
  inline  Bool_t SigGenCond(Float_t x,Float_t y);
  Int_t   Local2Pad(Float_t x,Float_t y,Int_t &padx,Int_t &pady)const;               //(x,y)->(padx,pady), returns sector code 
  Int_t   Local2PadX(Float_t x,Float_t y)    const {Int_t padx,pady;Local2Pad(x,y,padx,pady);return padx;}//(x,y)->padx
  Int_t   Local2PadY(Float_t x,Float_t y)    const {Int_t padx,pady;Local2Pad(x,y,padx,pady);return pady;}//(x,y)->pady
  void    Pad2Local(Int_t padx,Int_t pady,Float_t &x,Float_t &y);                                         //(padx,pady)->(x,y)
  Int_t   LocalX2Wire(Float_t x)             const {return  Int_t((x+PcSizeX()/2)/fWirePitch)+1;}         //x->wire number
  Float_t Wire2LocalX(Int_t iWireN)          const {return iWireN*fWirePitch-PcSizeX()/2;}                //wire number->x
  
  Float_t Gain(Float_t y);                                 //Returns total charge induced by single photon
  Float_t TotalCharge(Int_t iPID,Float_t eloss,Float_t y); //Returns total charge induced by particle lost eloss GeV
  Float_t PadCharge(Int_t /* iPadX */,Int_t /* iPadY */) {return 0;}   //Returns charge for a given pad
  void    FirstPad(Float_t x,Float_t y);
          
  static  Int_t   NpadsX()                   {return 144;}
  static  Int_t   NpadsY()                   {return 160;}   
  static  Int_t   NpadsXsec()                {return NpadsX()/3;}   
  static  Int_t   NpadsYsec()                {return NpadsY()/2;}   
  void    DeadZone(Float_t a)                {       fDeadZone=a;}
  Float_t DeadZone()                    const{return fDeadZone;}
  void    PadSize(Float_t x,Float_t y)       {       fPadSizeX=x;fPadSizeY=y;} 
  Float_t PadSizeX()                    const{return fPadSizeX;}
  Float_t PadSizeY()                    const{return fPadSizeY;}
  Float_t SectorSizeX()                 const{return NpadsX()*PadSizeX()/3;}
  Float_t SectorSizeY()                 const{return NpadsY()*PadSizeY()/2;}  
  Float_t PcSizeX()                     const{return NpadsX()*PadSizeX()+2*DeadZone();}
  Float_t PcSizeY()                     const{return NpadsY()*PadSizeY()+DeadZone();}
  Float_t WirePitch()                   const{return PadSizeX()/2;}
            
  void    Size(Float_t x,Float_t y,Float_t z){fSizeX=x;fSizeY=y;fSizeZ=z;}
  void    GeantSize(Float_t *pArr)      const{pArr[0]=fSizeX/2;pArr[1]=fSizeY/2;pArr[2]=fSizeZ/2;}  
  Float_t SizeX()                       const{return fSizeX;}
  Float_t SizeY()                       const{return fSizeY;}
  Float_t SizeZ()                       const{return fSizeZ;}   
  void    Offset(Float_t offset)             {       fOffset=offset;}  
  Float_t Offset()                      const{return fOffset;}  
  void    Angles(Float_t xy,Float_t yz)      {       fAngleXY=xy;fAngleYZ=yz;} 
  Float_t AngleYZ()                     const{return fAngleYZ*kD2r;} 
  Float_t AngleXY()                     const{return fAngleXY*kD2r;} 
  void    AngleRot(Float_t angle)            {       fAngleRot=angle;}
  Float_t AngleRot()                    const{return fAngleRot*kD2r;}                
  void    GapThickness(Float_t a)            {       fGapThickness=a;}    
  Float_t GapThickness()                const{return fGapThickness;}      
  void    ProximityGapThickness(Float_t a)   {       fProximityGapThickness=a;}
  Float_t ProximityGapThickness()       const{return fProximityGapThickness;}    
  void    QuartzLength(Float_t a)            {       fQuartzLength=a;}
  Float_t QuartzLength()                const{return fQuartzLength;}   
  void    QuartzWidth(Float_t a)             {       fQuartzWidth=a;}
  Float_t QuartzWidth()                 const{return fQuartzWidth;}
  void    QuartzThickness(Float_t a)         {       fQuartzThickness=a;}
  Float_t QuartzThickness()             const{return fQuartzThickness;}   
  void    OuterFreonLength(Float_t a)        {       fOuterFreonLength=a;}
  Float_t OuterFreonLength()            const{return fOuterFreonLength;}   
  void    OuterFreonWidth(Float_t a)         {       fOuterFreonWidth=a;}
  Float_t OuterFreonWidth()             const{return fOuterFreonWidth;}   
  void    InnerFreonLength(Float_t a)        {       fInnerFreonLength=a;}
  Float_t InnerFreonLength()            const{return fInnerFreonLength;}   
  void    InnerFreonWidth(Float_t a)         {       fInnerFreonWidth=a;}
  Float_t InnerFreonWidth()             const{return fInnerFreonWidth;}   
  void    FreonThickness(Float_t a)          {       fFreonThickness=a;}
  Float_t FreonThickness()              const{return fFreonThickness;}   
  void    RadiatorToPads(Float_t a)          {       fRadiatorToPads=a;}
  Float_t RadiatorToPads()              const{return fRadiatorToPads;}   
        
  void    SigmaIntegration(Float_t a)        {       fSigmaIntegration=a;}    
  Float_t SigmaIntegration()            const{return fSigmaIntegration;}    
  void    ChargeSpreadX(Float_t a)           {       fChargeSpreadX=a;}
  Float_t ChargeSpreadX()               const{return fChargeSpreadX;}    
  void    ChargeSpreadY(Float_t a)           {       fChargeSpreadY=a;}  
  Float_t ChargeSpreadY()               const{return fChargeSpreadY;}  
  Float_t AreaX()                       const{return fSigmaIntegration*fChargeSpreadX;} 
  Float_t AreaY()                       const{return fSigmaIntegration*fChargeSpreadY;} 
  void    ChargeSlope(Float_t a)             {       fChargeSlope=a;}
  Float_t ChargeSlope()                      {return fChargeSlope;}
  void    MaxAdc(Int_t a)                    {       fMaxAdc=a;}
  Int_t   MaxAdc()                      const{return fMaxAdc;}
  void    Pitch(Float_t a)                   {       fPitch=a;}
  Float_t Pitch()                       const{return fPitch;}
  void    AlphaFeedback(Float_t a)           {       fAlphaFeedback=a;}
  Float_t AlphaFeedback()               const{return fAlphaFeedback;}
  void    EIonisation(Float_t a)             {       fEIonisation=a;}
  Float_t EIonisation()                 const{return fEIonisation;}                            
  void    SqrtKx3(Float_t a)                 {       fSqrtKx3=a;};
  void    Kx2(Float_t a)                     {       fKx2=a;}
  void    Kx4(Float_t a)                     {       fKx4=a;}
  void    SqrtKy3(Float_t a)                 {       fSqrtKy3=a;}
  void    Ky2(Float_t a)                     {       fKy2=a;}
  void    Ky4(Float_t a)                     {       fKy4=a;}
  void    WireSag(Int_t a)                   {       fWireSag=a;}
  void    Voltage(Int_t a)                   {       fVoltage=a;}       
  Float_t Voltage()                     const{return fVoltage;}       
protected:
  Int_t   Local2Sector(Float_t &x,Float_t &y)const; //(x,y)->sector
  Int_t   Pad2Sector(Int_t &padx,Int_t &pady)const; //(padx,pady)->sector
  
  Float_t fDeadZone;                                //space between PC sectors, cm     
  Float_t fPadSizeX,fPadSizeY;                      //pad size, cm
  Float_t fWirePitch;                               //distance between wires along x
  
  Int_t   fCurrentPadX,fCurrentPadY;              //???
  Int_t   fCurrentWire;                           //???
    
  Float_t fSizeX;  Float_t fSizeY; Float_t fSizeZ;                                //chamber outer size, cm
  Float_t fAngleRot;                                                              //azimuthal rotation XY plane, deg  
  Float_t fAngleYZ;                                                               //angle between chambers YZ plane, deg
  Float_t fAngleXY;                                                               //angle between chambers XY plane, deg
  Float_t fOffset;                                                                //chambers offset from IP, cm   
  Float_t fGapThickness;                                                          //gap thickness, cm
  Float_t fProximityGapThickness;                                                 //proximity gap thickness, cm
  Float_t fQuartzLength;     Float_t fQuartzWidth;     Float_t fQuartzThickness;  //quartz window size, cm
  Float_t fOuterFreonLength; Float_t fOuterFreonWidth;                            //freon box outer size, cm
  Float_t fInnerFreonLength; Float_t fInnerFreonWidth;                            //freon box inner size, cm
  Float_t fFreonThickness;                                                        //freon thickness
  Float_t fRadiatorToPads;                                                        //distance from radiator to pads, cm
  
  Float_t fChargeSlope;              //Slope of the charge distribution
  Float_t fChargeSpreadX;            //Width of the charge distribution in x
  Float_t fChargeSpreadY;            //Width of the charge distribution in y
  Float_t fSigmaIntegration;         //Number of sigma's used for charge distribution
  Float_t fAlphaFeedback;            //Feedback photons coefficient
  Float_t fEIonisation;              //Mean ionisation energy
  Int_t   fMaxAdc;                   //Maximum ADC channel
  Float_t fSqrtKx3;                  //Mathieson parameters for x
  Float_t fKx2;                      //Mathieson parameters for x
  Float_t fKx4;                      //Mathieson parameters for x
  Float_t fSqrtKy3;                  //Mathieson parameters for y
  Float_t fKy2;                      //Mathieson parameters for y 
  Float_t fKy4;                      //Mathieson parameters for y
  Float_t fPitch;                    //Anode-cathode pitch
  Int_t   fWireSag;                  //Flag to turn on/off (0/1) wire sag
  Int_t   fVoltage;                  //Working voltage (2000, 2050, 2100, 2150)

  ClassDef(AliRICHParam,1)    //RICH main parameters
};
//__________________________________________________________________________________________________
void AliRICHParam::SigGenInit(Float_t x,Float_t y)
{//Initialises pad and wire position during stepping
  Local2Pad(x,y,fCurrentPadX,fCurrentPadY);
  fCurrentWire= (x>0) ? Int_t(x/fWirePitch)+1 : Int_t(x/fWirePitch)-1 ;
}
//__________________________________________________________________________________________________
Bool_t AliRICHParam::SigGenCond(Float_t x,Float_t y)
{//Signal will be generated if particle crosses pad boundary or boundary between two wires.
  Int_t curPadX,curPadY;
  Local2Pad(x,y,curPadX,curPadY);
  Int_t currentWire=(x>0) ? Int_t(x/fWirePitch)+1 : Int_t(x/fWirePitch)-1;
  if((curPadX != fCurrentPadX) || (curPadY != fCurrentPadY) || (currentWire!=fCurrentWire)) 
    return kTRUE;
  else
    return kFALSE;
}//Bool_t AliRICHParam::SigGenCond(Float_t x,Float_t y)
//__________________________________________________________________________________________________
Int_t AliRICHParam::Neighbours(Int_t iPadX,Int_t iPadY,Int_t listX[4],Int_t listY[4])const
{
  listX[0]=iPadX;   listY[0]=iPadY-1;       
  listX[1]=iPadX+1; listY[1]=iPadY;       
  listX[2]=iPadX;   listY[2]=iPadY+1;       
  listX[3]=iPadX-1; listY[3]=iPadY;       
  return 4;
}
//__________________________________________________________________________________________________
#endif //AliRICHParam_h
