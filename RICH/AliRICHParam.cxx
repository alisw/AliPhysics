#include "AliRICHParam.h"
#include "AliRICHConst.h"  //for units ????
 
ClassImp(AliRICHParam)

//______________________________________________________________________________
// RICH main parameters manipulator
AliRICHParam::AliRICHParam()
{//defines the default parameters
  Segmentation         (144,160);               //nx,ny   
  DeadZone             (3*cm);                  //spacer between PC planes
  PadSize              (8.4*mm,8.0*mm);   
  Size                 (80*cm,7*cm,60*cm);      //full length, not GEANT half notation
  AngleRot             (0*deg);                 //rotation of the whole RICH around Z
  AnglesDeg            (20,19.5);               //XY angle, YZ angle  deg  
  Offset               (490*cm+1.267*cm);       //1.267???????cm distance from IP to the center of module 
  GapThickness         (8*cm);               //Gap Thickness
  ProximityGapThickness(0.4*cm);             //Proximity Gap Thickness
  QuartzLength         (133*cm);             //Quartz Length
  QuartzWidth          (127.9*cm);           //Quartz Width
  QuartzThickness      (0.5*cm);             //Quartz Thickness
  OuterFreonLength     (133*cm);             //Outer Freon Length
  OuterFreonWidth      (41.3*cm);            //Outer Freon Width
  InnerFreonLength     (133*cm);             //Inner Freon Length
  InnerFreonWidth      (41.3*cm);            //Inner Freon Width
  FreonThickness       (1.5*cm);             //Freon Thickness
  RadiatorToPads       (0);               //Distance from radiator to pads
  
  SigmaIntegration(5.);
  ChargeSlope(27.);
  ChargeSpreadX(0.18);ChargeSpreadY(0.18);
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
//______________________________________________________________________________
void AliRICHParam::Recalc()
{//recalculate  
  Float_t csi_length=fNy*fPadY+fDeadZone;
  Float_t csi_width =fNx*fPadX+2*fDeadZone;
  fPadPlaneWidth = (csi_width - 2*fDeadZone)/3;
  fPadPlaneLength = (csi_length - fDeadZone)/2;  
}//void AliRICHParam::Recalc()
//______________________________________________________________________________
