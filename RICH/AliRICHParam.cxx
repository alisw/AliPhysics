#include "AliRICHParam.h"
#include "AliRICHConst.h"
#include <TMath.h>
#include <TRandom.h>
 
ClassImp(AliRICHParam)

// RICH main parameters manipulator
//__________________________________________________________________________________________________
AliRICHParam::AliRICHParam():
fNpadsX(0),fNpadsY(0),fNpadsXsector(0),fNpadsYsector(0),
fDeadZone(0),
fPadSizeX(0),fPadSizeY(0),
fSectorSizeX(0),fSectorSizeY(0),
fWirePitch(0),
fCurrentPadX(0),fCurrentPadY(0),fCurrentWire(0),
fSizeZ(0),
fAngleRot(0),fAngleYZ(0),fAngleXY(0),
fOffset(0),
fGapThickness(0),
fProximityGapThickness(0),
fQuartzLength(0),
fQuartzWidth(0),
fQuartzThickness(0),
fOuterFreonLength(0),
fOuterFreonWidth(0),
fInnerFreonLength(0),
fInnerFreonWidth(0),
fFreonThickness(0),
fRadiatorToPads(0),
fPcSizeX(0),
fPcSizeY(0),
fChargeSlope(0),
fChargeSpreadX(0),
fChargeSpreadY(0),
fSigmaIntegration(0),
fAlphaFeedback(0),
fEIonisation(0),
fMaxAdc(0),
fSqrtKx3(0),
fKx2(0),
fKx4(0),
fSqrtKy3(0),
fKy2(0),
fKy4(0),
fPitch(0),
fWireSag(0),
fVoltage(0)
{//defines the default parameters
  Segmentation         (144,160);           //nx,ny  for the whole chamber
  DeadZone             (3*kcm);             //spacer between PC planes
  PadSize              (8.4*kmm,8.0*kmm);     
  fWirePitch=PadSizeX()/2;
  
  Size                 (132.6*kcm,26*kcm,136.7*kcm);  //full length, not GEANT half notation
  AngleRot             (-60);                         //rotation of the whole RICH around Z, deg
  Angles               (20,19.5);                     //XY angle, YZ angle  deg  
  Offset               (490*kcm+1.267*kcm);           //1.267???????cm distance from IP to the center of module 
  GapThickness         (8*kcm);              
  ProximityGapThickness(0.4*kcm);            
  QuartzLength         (133*kcm);            
  QuartzWidth          (127.9*kcm);          
  QuartzThickness      (0.5*kcm);            
  OuterFreonLength     (133*kcm);            
  OuterFreonWidth      (41.3*kcm);           
  InnerFreonLength     (133*kcm);            
  InnerFreonWidth      (41.3*kcm);           
  FreonThickness       (1.5*kcm);            
  RadiatorToPads       (80*kmm);                 
  
  ChargeSlope(27.);
  ChargeSpreadX(0.18);ChargeSpreadY(0.18);
  SigmaIntegration(5.);
  MaxAdc(4096);
  AlphaFeedback(0.036);
  EIonisation(26.e-9);
  SqrtKx3(0.77459667);
  Kx2(0.962);
  Kx4(0.379);
  SqrtKy3(0.77459667);
  Ky2(0.962);
  Ky4(0.379);
  Pitch(0.25);
  WireSag(1);		      // 1->On, 0->Off
  Voltage(2150);	      // Should only be 2000, 2050, 2100 or 2150  
  
  Recalc();
}//AliRICHParam::named ctor 
//__________________________________________________________________________________________________
void AliRICHParam::Recalc()
{//recalculate  
  fNpadsXsector=NpadsX()/3;  fNpadsYsector=NpadsY()/2;
  fPcSizeX=NpadsX()*fPadSizeX+2*fDeadZone;
  fPcSizeY=NpadsY()*fPadSizeY+fDeadZone;
  fSectorSizeX=(fPcSizeX-2*fDeadZone)/3;
  fSectorSizeY=(fPcSizeY-fDeadZone)/2;  
}//void AliRICHParam::Recalc()
//__________________________________________________________________________________________________
Int_t AliRICHParam::Sector(Float_t &x, Float_t &y)const
{//Determines sector for a given hit (x,y) and trasform this point to the local system of that sector.
  
  Int_t sector=kBad;  
  if(x<=-fSectorSizeX/2-fDeadZone&&x>=-fPcSizeX/2)     {sector=1;x+=fPcSizeX/2;}
  else if(x>=-fSectorSizeX/2 && x<=fSectorSizeX/2)     {sector=2;x+=fSectorSizeX/2;}
  else if(x>= fSectorSizeX/2+fDeadZone&&x<=fPcSizeX/2) {sector=3;x-=fSectorSizeX/2+fDeadZone;}
  else if(x<-fPcSizeX/2||x>fPcSizeX/2)                 {Error("Sector","given x position is out of active PC area");return kBad;}
  else                                                 {return kBad;} //in dead zone

  if(y>=-fPcSizeY/2&&y<= -fDeadZone/2)                {y+=fPcSizeY/2;  return -sector;}
  else if(y>-fDeadZone/2&&y<fDeadZone/2)              {return kBad;} //in dead zone
  else if(y>=fDeadZone/2&&y<=fPcSizeY/2)              {y-=fDeadZone/2; return  sector;}
  else                                                {Error("Sector","given y position is out of active PC area");return kBad;}
}//Int_t AliRICHParam::Sector(Float_t x, Float_t y)
//__________________________________________________________________________________________________
Int_t AliRICHParam::L2P(Float_t x, Float_t y, Int_t &padx, Int_t &pady)const
{//returns pad numbers (iPadX,iPadY) for given point in local coordinates (x,y) 
 //count starts in lower left corner from 1,1 to 144,180
  
  padx=pady=kBad;
  Int_t sector=Sector(x,y);
  if(sector==kBad) return sector;
  
  padx=Int_t(x/fPadSizeX)+1; 
  if(padx>fNpadsXsector)          padx=fNpadsXsector;
  if(sector==2||sector==-2)       padx+=fNpadsXsector;
  else if(sector==3||sector==-3)  padx+=fNpadsXsector*2;
  
  pady=Int_t(y/fPadSizeY)+1;
  if(pady>fNpadsYsector)          padx=fNpadsYsector;
  if(sector>0)                    pady+=fNpadsYsector;    

  return sector;
}//void AliRICHParam::L2P(Float_t x, Float_t y, Int_t &padx, Int_t &pady)
//__________________________________________________________________________________________________
Float_t AliRICHParam::Gain(Float_t y)
{//Calculates the gain
  if(fWireSag){
    Float_t gainK=9e-6*TMath::Power(y,4)+2e-7*TMath::Power(y,3)-0.0316*TMath::Power(y,2)-3e-4*y+25.367;
    Float_t gain = (ChargeSlope()+ChargeSlope()*gainK/100)*0.9;
    return -gain*TMath::Log(gRandom->Rndm());
  }else     
    return -ChargeSlope()*TMath::Log(gRandom->Rndm());
}//Float_t AliRICHParam::IntPH(Float_t yhit)
//__________________________________________________________________________________________________
Float_t AliRICHParam::TotalCharge(Int_t iPID,Float_t eloss,Float_t y)
{//Get number of electrons and return charge
    
  if(iPID>50000)//it's photon no more then 1 electron after photoelectron conversion
    return Gain(y);
  else{  
    Int_t iNelectrons=Int_t(eloss/fEIonisation);if(iNelectrons==0) iNelectrons=1;
    Float_t charge=0;
    for(Int_t i=1;i<=iNelectrons;i++)
      charge+=Gain(y);
    return charge;
  }
}//Float_t AliRICHParam::TotalCharge(Int_t iPID,Float_t eloss, Float_t y)
//__________________________________________________________________________________________________
void AliRICHParam::FirstPad(Float_t x,Float_t y)
{
  Int_t padx,pady;
  L2P(x,y,padx,pady);
}//void AliRICHParam::FirstPad(Float_t x,Float_t y)
