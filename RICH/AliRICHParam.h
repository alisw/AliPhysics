#ifndef AliRICHParam_h
#define AliRICHParam_h

#include <TObject.h>
#include "AliRICHConst.h"

class AliRICHParam :public TObject  
{
public:
           AliRICHParam();  
  virtual ~AliRICHParam()                    {;}  
  
  void    Recalc();//Recalculate dependent parameters after changes applied  
  void    Segmentation(Int_t Nx, Int_t Ny)   {fNx=Nx;fNy=Ny;Recalc();}
  Int_t   Nx()                          const{return fNx;}
  Int_t   Ny()                          const{return fNy;}   
  void    DeadZone(Float_t a)                {       fDeadZone=a;Recalc();}
  Float_t DeadZone()                    const{return fDeadZone;}
  void    PadSize(Float_t x,Float_t y)       {       fPadX=x;fPadY=y;Recalc();} 
  Float_t PadX()                        const{return fPadX;}
  Float_t PadY()                        const{return fPadY;}
  Float_t PadPlaneWidth()               const{return fPadPlaneWidth;}
  Float_t PadPlaneLength()              const{return fPadPlaneLength;}  

  void    Size(Float_t x,Float_t y,Float_t z){fSizeX=x;fSizeY=y;fSizeZ=z;}
  void    GeantSize(Float_t *pArr)      const{pArr[0]=fSizeX/2;pArr[1]=fSizeY/2;pArr[2]=fSizeZ/2;}  
  Float_t SizeX()                       const{return fSizeX;}
  Float_t SizeY()                       const{return fSizeY;}
  Float_t SizeZ()                       const{return fSizeZ;}   
  void    Offset(Float_t offset)             {       fOffset=offset;}  
  Float_t Offset()                      const{return fOffset;}  
  void    Angles(Float_t xy,Float_t yz)      {       fAngleXY=xy;fAngleYZ=yz;} 
  Float_t AngleYZ()                     const{return fAngleYZ*d2r;} 
  Float_t AngleXY()                     const{return fAngleXY*d2r;} 
  void    AngleRot(Float_t angle)            {       fAngleRot=angle;}
  Float_t AngleRot()                    const{return fAngleRot*d2r;}                
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
  void    ChargeSlope(Float_t a)             {       fChargeSlope=a;}
  Float_t ChargeSlope()                      {return fChargeSlope;}
  void    MaxAdc(Float_t a)                  {       fMaxAdc=a;}
  Float_t MaxAdc()                      const{return fMaxAdc;}
  void    Pitch(Float_t a)                   {       fPitch=a;};
  Float_t Pitch()                       const{return fPitch;}
  void    AlphaFeedback(Float_t a)           {       fAlphaFeedback=a;}
  Float_t AlphaFeedback()               const{return fAlphaFeedback;}
  void    EIonisation(Float_t a)             {       fEIonisation=a;}
  Float_t EIonisation()                 const{return fEIonisation;}                            
  void    SqrtKx3(Float_t a)                 {       fSqrtKx3=a;};
  void    Kx2(Float_t a)                     {       fKx2=a;};
  void    Kx4(Float_t a)                     {       fKx4=a;};
  void    SqrtKy3(Float_t a)                 {       fSqrtKy3=a;};
  void    Ky2(Float_t a)                     {       fKy2=a;};
  void    Ky4(Float_t a)                     {       fKy4=a;};
  void    WireSag(Int_t a)                   {       fWireSag=a;};
  void    Voltage(Int_t a)                   {       fVoltage=a;};       
protected:
  Int_t   fNx;                //number of pads along X
  Int_t   fNy;                //number of pads along Y
  Float_t fDeadZone;          //spacer between PC planes, cm     
  Float_t fPadX;              //pad width, cm
  Float_t fPadY;              //pad lenght, cm
  Float_t fPadPlaneWidth;     //pad plane width, cm
  Float_t fPadPlaneLength;    //pad plane length, cm
  
  Float_t fSizeX;             //chamber length, cm
  Float_t fSizeY;             //chamber thickness, cm
  Float_t fSizeZ;             //chamber width, cm
  Float_t fAngleRot;          //azimuthal rotation angle in X-Y plane, deg  
  Float_t fAngleYZ;           //angle between RICH chambers in YZ plane, deg
  Float_t fAngleXY;           //angle between RICH chambers in XY plane, deg
  Float_t fOffset;            //chambers offset from IP, cm   
  Float_t fGapThickness;            //gap thickness, cm
  Float_t fProximityGapThickness;   //proximity gap thickness, cm
  Float_t fQuartzLength;            //quartz length, cm
  Float_t fQuartzWidth;             //quartz width, cm
  Float_t fQuartzThickness;         //quartz thickness, cm
  Float_t fOuterFreonLength;        //outer freon length, cm
  Float_t fOuterFreonWidth;         //outer freon width, cm
  Float_t fInnerFreonLength;        //inner freon length, cm
  Float_t fInnerFreonWidth;         //inner freon width, cm
  Float_t fFreonThickness;          //freon thickness
  Float_t fRadiatorToPads;          //distance from radiator to pads, cm

  Float_t fChargeSlope;              //Slope of the charge distribution
  Float_t fChargeSpreadX;            //Width of the charge distribution in x
  Float_t fChargeSpreadY;            //Width of the charge distribution in y
  Float_t fSigmaIntegration;         //Number of sigma's used for charge distribution
  Float_t fAlphaFeedback;            //Feedback photons coefficient
  Float_t fEIonisation;              //Mean ionisation energy
  Float_t fMaxAdc;                   //Maximum ADC channel
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

#endif //AliRICHParam_h
