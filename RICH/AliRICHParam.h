#ifndef AliRICHParam_h
#define AliRICHParam_h

#include <TObject.h>
#include "AliRICHConst.h"

class AliRICHParam :public TObject  
{
public:
                 AliRICHParam();  
         void    Recalc();//Recalculate dependent parameters after changes 
  inline void    Segmentation(Int_t Nx, Int_t Ny)   {fNx=Nx;fNy=Ny;Recalc();}
  inline Int_t   Nx()                          const{return fNx;}
  inline Int_t   Ny()                          const{return fNy;}   
  inline void    DeadZone(Float_t a)                {       fDeadZone=a;Recalc();}
  inline Float_t DeadZone()                    const{return fDeadZone;}
  inline void    PadSize(Float_t x,Float_t y)       {       fPadX=x;fPadY=y;Recalc();} 
  inline Float_t PadX()                        const{return fPadX;}
  inline Float_t PadY()                        const{return fPadY;}
  inline Float_t PadPlaneWidth()               const{return fPadPlaneWidth;}
  inline Float_t PadPlaneLength()              const{return fPadPlaneLength;}  
  inline void    Size(Float_t x,Float_t y,Float_t z){fSizeX=x;fSizeY=y;fSizeZ=z;}
  inline void    GeantSize(Float_t *pParam)    const{pParam[0]=fSizeX/2;pParam[1]=fSizeY/2;pParam[2]=fSizeZ/2;}  
  inline Float_t SizeX()                       const{return fSizeX;}
  inline Float_t SizeY()                       const{return fSizeY;}
  inline Float_t SizeZ()                       const{return fSizeZ;}   
  inline void    Offset(Float_t offset)             {       fOffset=offset;}  
  inline Float_t Offset()                      const{return fOffset;}  
  inline void    AnglesDeg(Float_t xy,Float_t yz)   {       fAngleXY=xy;fAngleYZ=yz;} 
  inline Float_t AngleYZ()                     const{return fAngleYZ*d2r;} 
  inline Float_t AngleXY()                     const{return fAngleXY*d2r;} 
  inline void    AngleRot(Float_t angle)            {       fAngleRot=angle;}
  inline Float_t AngleRot()                    const{return fAngleRot*d2r;}                
  inline void    GapThickness(Float_t a)            {       fGapThickness=a;}    
  inline Float_t GapThickness()                const{return fGapThickness;}      
  inline void    ProximityGapThickness(Float_t a)   {       fProximityGapThickness=a;}
  inline Float_t ProximityGapThickness()       const{return fProximityGapThickness;}    
  inline void    QuartzLength(Float_t a)            {       fQuartzLength=a;}
  inline Float_t QuartzLength()                const{return fQuartzLength;}   
  inline void    QuartzWidth(Float_t a)             {       fQuartzWidth=a;}
  inline Float_t QuartzWidth()                 const{return fQuartzWidth;}
  inline void    QuartzThickness(Float_t a)         {       fQuartzThickness=a;}
  inline Float_t QuartzThickness()             const{return fQuartzThickness;}   
  inline void    OuterFreonLength(Float_t a)        {       fOuterFreonLength=a;}
  inline Float_t OuterFreonLength()            const{return fOuterFreonLength;}   
  inline void    OuterFreonWidth(Float_t a)         {       fOuterFreonWidth=a;}
  inline Float_t OuterFreonWidth()             const{return fOuterFreonWidth;}   
  inline void    InnerFreonLength(Float_t a)        {       fInnerFreonLength=a;}
  inline Float_t InnerFreonLength()            const{return fInnerFreonLength;}   
  inline void    InnerFreonWidth(Float_t a)         {       fInnerFreonWidth=a;}
  inline Float_t InnerFreonWidth()             const{return fInnerFreonWidth;}   
  inline void    FreonThickness(Float_t a)          {       fFreonThickness=a;}
  inline Float_t FreonThickness()              const{return fFreonThickness;}   
  inline void    RadiatorToPads(Float_t a)          {       fRadiatorToPads=a;}
  inline Float_t RadiatorToPads()              const{return fRadiatorToPads;}   
               
  inline void    SigmaIntegration(Float_t a)        {       fSigmaIntegration=a;}    
  inline Float_t SigmaIntegration()            const{return fSigmaIntegration;}    
  inline void    ChargeSpreadX(Float_t a)           {       fChargeSpreadX=a;}
  inline Float_t ChargeSpreadX()               const{return fChargeSpreadX;}    
  inline void    ChargeSpreadY(Float_t a)           {       fChargeSpreadY=a;}
  inline Float_t ChargeSpreadY()               const{return fChargeSpreadY;}   
  inline void    ChargeSlope(Float_t a)             {       fChargeSlope=a;}
  inline Float_t ChargeSlope()                      {return fChargeSlope;}
  inline void    MaxAdc(Float_t a)                  {       fMaxAdc=a;}
  inline Float_t MaxAdc()                      const{return fMaxAdc;}
  inline void    Pitch(Float_t a)                   {       fPitch=a;};
  inline Float_t Pitch()                       const{return fPitch;}
  inline void    AlphaFeedback(Float_t a)           {       fAlphaFeedback=a;}
  inline Float_t AlphaFeedback()               const{return fAlphaFeedback;}
  inline void    EIonisation(Float_t a)             {       fEIonisation=a;}
  inline Float_t EIonisation()                 const{return fEIonisation;}                            
  inline void    SqrtKx3(Float_t a)                 {       fSqrtKx3=a;};
  inline void    Kx2(Float_t a)                     {       fKx2=a;};
  inline void    Kx4(Float_t a)                     {       fKx4=a;};
  inline void    SqrtKy3(Float_t a)                 {       fSqrtKy3=a;};
  inline void    Ky2(Float_t a)                     {       fKy2=a;};
  inline void    Ky4(Float_t a)                     {       fKy4=a;};
  inline void    WireSag(Int_t a)                   {       fWireSag=a;};
  inline void    Voltage(Int_t a)                   {       fVoltage=a;};       
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
  Float_t fAngleRot;          //azimuthal rotation angle in X-Y plane, grad  
  Float_t fAngleYZ;           //angle between RICH chambers in YZ plane, grad
  Float_t fAngleXY;           //angle between RICH chambers in XY plane, grad
  Float_t fOffset;            //chambers offset from IP, cm   
  Float_t fGapThickness;            //gap thickness, cm
  Float_t fProximityGapThickness;   //proximity gap thickness, cm
  Float_t fQuartzLength;            //quartz length
  Float_t fQuartzWidth;             //quartz width
  Float_t fQuartzThickness;         //quartz thickness
  Float_t fOuterFreonLength;        //outer freon length
  Float_t fOuterFreonWidth;         //outer freon width
  Float_t fInnerFreonLength;        //inner freon length
  Float_t fInnerFreonWidth;         //inner freon width
  Float_t fFreonThickness;          //freon thickness
  Float_t fRadiatorToPads;          //distance from radiator to pads

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
