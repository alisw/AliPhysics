#ifndef AliRICHParam_h
#define AliRICHParam_h

#include <TObject.h>
#include <TMath.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TError.h>
#include <TObjArray.h>


static const int kNCH=7;           //number of RICH chambers ???
static const int kNchambers=7;     //number of RICH chambers 
static const int kNpadsX = 160;    //number of pads along X in single chamber
static const int kNpadsY = 144;    //number of pads along Y in single chamber
static const int kBad=-101;        //useful static const to mark initial (uninitalised) values
static const int kNsectors=6;      //number of sectors per chamber

static const int kadc_satm  =  4096;  //dynamic range (10 bits)
static const int kCerenkov=50000050;  //??? go to something more general like TPDGCode
static const int kFeedback=50000051;  //??? go to something more general like TPDGCode

class AliRICHChamber;

class AliRICHParam :public TObject  
{
public:
                  AliRICHParam():TObject(),fpChambers(0)  {CreateChambers();}
  virtual        ~AliRICHParam()                          {delete fpChambers;}
         void     CreateChambers();
         AliRICHChamber* C(Int_t i)          {return (AliRICHChamber*)fpChambers->UncheckedAt(i-1);}      //returns pointer to chamber i
  static Int_t    NpadsX()                   {return kNpadsX;}                           //pads along X in chamber
  static Int_t    NpadsY()                   {return kNpadsY;}                           //pads along Y in chamber
  static Int_t    NpadsXsec()                {return NpadsX()/2;}                        //pads along X in sector
  static Int_t    NpadsYsec()                {return NpadsY()/3;}                        //pads along Y in sector
  static Double_t DeadZone()                 {return 2.6;}                               //dead zone size in cm  
  static Double_t PadSizeX()                 {return 0.8;}                               //pad size x in cm 
  static Double_t PadSizeY()                 {return 0.84;}                              //pad size y in cm   
  static Double_t SectorSizeX()              {return NpadsX()*PadSizeX()/2;}             //sector size x in cm
  static Double_t SectorSizeY()              {return NpadsY()*PadSizeY()/3;}             //sector size y in cm 
  static Double_t PcSizeX()                  {return NpadsX()*PadSizeX()+DeadZone();}    //photocathode size x in cm
  static Double_t PcSizeY()                  {return NpadsY()*PadSizeY()+2*DeadZone();}  //photocathode size y in cm 
  static Double_t SizeX()                    {return 132.6;}
  static Double_t SizeY()                    {return 26;}
  static Double_t SizeZ()                    {return 136.7;}                             
  static Double_t Offset()                   {return 490+1.267;}                         //distance from IP to center of chamber in cm 
  static Double_t AngleYZ()                  {return 19.5*TMath::DegToRad();}            //angle between chambers in YZ plane, rad
  static Double_t AngleXY()                  {return 20*TMath::DegToRad();}              //angle between chambers in XY plane, rad
  static Double_t AngleRot()                 {return fgAngleRot*TMath::DegToRad();}      //RICH rotation around Z, rad
  static Double_t FreonThickness()           {return 1.5;}   
  static Double_t QuartzThickness()          {return 0.5;}   
  
  static Double_t GapProx()                  {return 8.0;}                               //cm between CsI PC and radiator quartz window
  static Double_t GapColl()                  {return 7.0;}                               //cm between CsI PC and third wire grid (collection wires)     
  static Double_t GapAnod()                  {return 0.204;}                             //cm between CsI PC and first wire grid (anod wires)     
  static Double_t GapAmp()                   {return 0.445;}                             //cm between CsI PC and second wire grid (cathode wires)
  static Double_t PitchAnod()                {return PadSizeY()/2;}                      //cm between anode wires
  static Double_t PitchCath()                {return PadSizeY()/4;}                      //cm between cathode wires
  static Double_t PitchColl()                {return 0.5;}                               //cm between collect wires
  
  static Double_t GapThickness()             {return 8.0;}      
  static Double_t RadiatorToPads()           {return FreonThickness()+QuartzThickness()+GapThickness();}   
  static Double_t AnodeCathodeGap()          {return 0.2;}                               //between CsI PC and first wire grid     
  static Double_t QuartzLength()             {return 133;}   
  static Double_t QuartzWidth()              {return 127.9;}
  static Double_t OuterFreonLength()         {return 133;}   
  static Double_t OuterFreonWidth()          {return 41.3;}   
  static Double_t InnerFreonLength()         {return 133;}   
  static Double_t InnerFreonWidth()          {return 41.3;}   
  static Double_t IonisationPotential()      {return 26.0e-9;}                            //for CH4 in GeV taken from ????
  static TVector2 MathiesonDelta()           {return TVector2(5*0.18,5*0.18);}            //area of 5 sigmas of Mathieson distribution (cm)
  static Int_t    MaxQdc()                   {return 4095;}                               //QDC number of channels          
  static Double_t AlphaFeedback(Int_t )      {return 0.030;}                              //determines number of feedback photons
    
  static Bool_t   IsResolveClusters()         {return fgIsResolveClusters;}  //go after resolved clusters?
  static Bool_t   IsWireSag()                 {return fgIsWireSag;}          //take wire sagita in account?
  static Bool_t   IsRadioSrc()                {return fgIsRadioSrc;}         //add radioactive source inside CH4?
  static Int_t    HV(Int_t sector)            {
    if (sector>=1 && sector <=6)
      return fgHV[sector-1];
    else {
      ::Error("HV","Wrong sector %d",sector);
      return kBad;
    } 
  }       //high voltage for this sector
  static void     SetDeclustering(Bool_t a)   {fgIsResolveClusters=a;}  
  static void     SetRadioSrc(Bool_t a)       {fgIsRadioSrc=a;}  
  static void     SetWireSag(Bool_t status)   {fgIsWireSag=status;}  
  static void     SetHV(Int_t sector,Int_t hv){fgHV[sector-1]=hv;}  
  static void     SetAngleRot(Double_t rot)   {fgAngleRot =rot;}

  inline static TVector  Loc2Area(TVector2 x2);                                                    //return area affected by hit x2
  inline static TVector  Loc2Pad(TVector2 x2);                                                     //return pad containing given position
  inline static TVector2 Pad2Loc(TVector pad);                                                     //return center of the pad
         static TVector2 Pad2Loc(Int_t x,Int_t y) {TVector pad(2);pad[0]=x;pad[1]=y;return Pad2Loc(pad);}
  inline static Int_t    PadNeighbours(Int_t iPadX,Int_t iPadY,Int_t aListX[4],Int_t aListY[4]);   //number of neighbours for this pad
  
  inline static Double_t Mathieson(Double_t x1,Double_t x2,Double_t y1,Double_t y2);               //Mathienson integral over these limits
  inline static Double_t GainSag(Double_t x,Int_t sector);                                         //gain variations in %
         static Double_t QdcSlope(Int_t sec){switch(sec){case kBad: return 0;  default:   return 33;}} //weight of electon in QDC channels
         static Double_t Gain(TVector2 x2){if(IsWireSag()) return QdcSlope(Loc2Sec(x2))*(1+GainSag(x2.X(),Loc2Sec(x2))/100);else return QdcSlope(Loc2Sec(x2));}//gain for point in chamber RS 
  inline static Double_t FracQdc(TVector2 x2,TVector pad);                                         //charge fraction to pad from hit
  inline static Int_t    TotQdc(TVector2 x2,Double_t eloss);                                       //total charge for hit eloss=0 for photons
  inline        Bool_t   IsOverTh(Int_t c,TVector pad,Double_t q);                                 //is QDC of the pad registered by FEE  
         static Int_t    NsigmaTh()                    {return fgNsigmaTh;}                        //
         static Float_t  SigmaThMean()                 {return fgSigmaThMean;}                     //QDC electronic noise mean
         static Float_t  SigmaThSpread()               {return fgSigmaThSpread;}                   //QDC electronic noise width
                void     Print(const Option_t *opt="");                                            //virtual
  inline static void     PropogateHelix(TVector3 x0,TVector3 p0,Double_t s,TVector3 *x,TVector3 *p);                
                
  inline static Int_t    Loc2Sec(TVector2 &x2);             //return sector, x2->Sector RS
  inline static Int_t    Pad2Sec(TVector pad);              //return sector
protected:
         TObjArray *fpChambers;                             //list of chambers    
  static Bool_t     fgIsWireSag;                            //wire sagitta ON/OFF flag
  static Bool_t     fgIsResolveClusters;                    //declustering ON/OFF flag
  static Bool_t     fgIsRadioSrc;                           //radioactive source ON/OFF flag
  static Int_t      fgHV[6];                                //HV applied to anod wires
  static Double_t   fgAngleRot;                             //module rotation from up postion (0,0,490)cm
  static Int_t      fgNsigmaTh;                             //n. of sigmas to cut for zero suppression
  static Float_t    fgSigmaThMean;                          //sigma threshold value
  static Float_t    fgSigmaThSpread;                        //spread of sigma
  ClassDef(AliRICHParam,5)                                  //RICH main parameters class
};
//__________________________________________________________________________________________________
Int_t AliRICHParam::PadNeighbours(Int_t iPadX,Int_t iPadY,Int_t listX[4],Int_t listY[4])
{
// Determines all the neighbouring pads for the given one (iPadX,iPadY). Returns total number of these pads.
// Dead zones are taken into account.    
//   1  
// 2   3
//   4     
  Int_t nPads=0;
  if(iPadY!=NpadsY()&&iPadY!=2*NpadsYsec()&&iPadY!=NpadsYsec()){listX[nPads]=iPadX;   listY[nPads]=iPadY+1; nPads++;}       //1
  if(iPadX!=1&&iPadX!=NpadsXsec()+1)                           {listX[nPads]=iPadX-1; listY[nPads]=iPadY;   nPads++;}       //2
  if(iPadX!=NpadsXsec()&&iPadX!=NpadsX())                      {listX[nPads]=iPadX+1; listY[nPads]=iPadY;   nPads++;}       //3
  if(iPadY!=1&&iPadY!=NpadsYsec()+1&&2*NpadsYsec()+1)          {listX[nPads]=iPadX;   listY[nPads]=iPadY-1; nPads++;}       //4

  return nPads;
}//Pad2ClosePads()
//__________________________________________________________________________________________________
Int_t AliRICHParam::Loc2Sec(TVector2 &v2)
{
// Determines sector containing the given point and trasform this point to the local system of that sector.
// Returns sector code:                       
//y ^  5 6
//  |  3 4
//  |  1 2
//   -------> x  
  Double_t x0=0; Double_t x1=SectorSizeX(); Double_t x2=SectorSizeX()+DeadZone(); Double_t x3=PcSizeX();
  Double_t y0=0; Double_t y1=SectorSizeY(); Double_t y2=SectorSizeY()+DeadZone(); Double_t y3=2*SectorSizeY()+DeadZone(); 
  Double_t y4=PcSizeY()-SectorSizeY();      Double_t y5=PcSizeY();
  
  Int_t sector=kBad;  
  Double_t x=v2.X(),y=v2.Y();  
  if     (v2.X() >= x0 && v2.X() <= x1 )  {sector=1;}
  else if(v2.X() >= x2 && v2.X() <= x3 )  {sector=2; x=v2.X()-x2;}
  else                                    {sector=kBad; ::Error("Loc2Sec","Position %6.2f,%6.2f is out of chamber in X",v2.X(),v2.Y());return kBad;}
  
  if     (v2.Y() >= y0 && v2.Y() <= y1 )  {}                                  //sectors 1 or 2 
  else if(v2.Y() >= y2 && v2.Y() <= y3 )  {sector+=2; y=v2.Y()-y2;}           //sectors 3 or 4
  else if(v2.Y() >= y4 && v2.Y() <= y5 )  {sector+=4; y=v2.Y()-y4;}           //sectors 5 or 6
  else                                    {sector=kBad; ::Error("Loc2Sec","Position %6.2f,%6.2f is out of chamber in Y",v2.X(),v2.Y());return kBad;}
  v2.Set(x,y);
  return sector;
}//Loc2Sec(Double_t x, Double_t y)
//__________________________________________________________________________________________________
TVector AliRICHParam::Loc2Pad(TVector2 x2)
{
// Determines pad number TVector(padx,pady) containing the given point x2 defined the chamber RS.
// Pad count starts in lower left corner from 1,1 to 144,160 in upper right corner of a chamber.
// Returns sector number of the determined pad.      
//y ^  5 6
//  |  3 4
//  |  1 2
//   -------> x  
  TVector pad(2);
  Int_t sector=Loc2Sec(x2);//trasforms x2 to sector reference system
  if(sector==kBad) {pad[0]=pad[1]=kBad; return pad;}
  
  pad[0]=Int_t(x2.X()/PadSizeX())+1; if(pad[0]>NpadsXsec()) pad[0]= NpadsXsec();       
  if(sector==2||sector==4||sector==6)   pad[0]+=  NpadsXsec();     

  pad[1]=Int_t(x2.Y()/PadSizeY())+1; if(pad[1]>NpadsYsec()) pad[1]= NpadsYsec();
  if(sector==3||sector==4)   pad[1]+=NpadsYsec();    
  if(sector==5||sector==6)   pad[1]+=2*NpadsYsec();     
  return pad;
}
//__________________________________________________________________________________________________
Int_t AliRICHParam::Pad2Sec(TVector pad)
{
// Determines sector containing the given pad.
  Int_t sector=kBad;      
  if     (pad[0] >= 1           && pad[0] <=   NpadsXsec() )    {sector=1;}
  else if(pad[0] >  NpadsXsec() && pad[0] <=   NpadsX()    )    {sector=2;} 
  else                                                         ::Error("Pad2Sec","Wrong pad (%3.0f,%3.0f)",pad[0],pad[1]);
    
  if     (pad[1] >= 1             && pad[1] <=   NpadsYsec() )    {}
  else if(pad[1] >  NpadsYsec()   && pad[1] <= 2*NpadsYsec() )    {sector+=2;}
  else if(pad[1] >  2*NpadsYsec() && pad[1] <=   NpadsY()    )    {sector+=4;}
  else                                                         ::Error("Pad2Sec","Wrong pad (%3.0f,%3.0f)",pad[0],pad[1]);

  return sector;
}//Pad2Sec()
//__________________________________________________________________________________________________
TVector2 AliRICHParam::Pad2Loc(TVector pad)
{
// Returns position of the center of the given pad in local system of the chamber    
// y ^  5 6
//   |  3 4        chamber structure
//   |  1 2
//    -------> x  
  Double_t x=kBad,y=kBad;
  if(pad[0] > 0 && pad[0] <= NpadsXsec())//it's 1 or 3 or 5
    x=(pad[0]-0.5)*PadSizeX();
  else if(pad[0] > NpadsXsec() && pad[0] <= NpadsX())//it's 2 or 4 or 6
    x=(pad[0]-0.5)*PadSizeX()+DeadZone();
  else
    ::Error("Pad2Loc","Wrong pad (%3.0f,%3.0f)",pad[0],pad[1]);
  
  if(pad[1] > 0 && pad[1] <= NpadsYsec())//it's 1 or 2
    y=(pad[1]-0.5)*PadSizeY();
  else if(pad[1] > NpadsYsec() && pad[1] <= 2*NpadsYsec())//it's 3 or 4
    y=(pad[1]-0.5)*PadSizeY()+DeadZone();
  else if(pad[1] > 2*NpadsYsec() && pad[1]<= NpadsY())//it's 5 or 6
    y=(pad[1]-0.5)*PadSizeY()+2*DeadZone();
  else
    ::Error("Pad2Loc","Wrong pad (%3.0f,%3.0f)",pad[0],pad[1]);
    
  return TVector2(x,y);
}
//__________________________________________________________________________________________________
Double_t AliRICHParam::GainSag(Double_t x,Int_t sector)
{
// Returns % of gain variation due to wire sagita.
// All curves are parametrized as per sector basis, so y must be apriory transformed to the Sector RS.    
  x-=SectorSizeX()/2;
  if(x>SectorSizeX()) x-=SectorSizeX(); 
  switch(HV(sector)){
    case 2150: return 9e-6*TMath::Power(x,4)+2e-7*TMath::Power(x,3)-0.0316*TMath::Power(x,2)-3e-4*x+25.367;//%
    case 2100: return 8e-6*TMath::Power(x,4)+2e-7*TMath::Power(x,3)-0.0283*TMath::Power(x,2)-2e-4*x+23.015;
    case 2050: return 7e-6*TMath::Power(x,4)+1e-7*TMath::Power(x,3)-0.0254*TMath::Power(x,2)-2e-4*x+20.888;
    case 2000: return 6e-6*TMath::Power(x,4)+8e-8*TMath::Power(x,3)-0.0227*TMath::Power(x,2)-1e-4*x+18.961;
    default:   return 0;
  }
}
//__________________________________________________________________________________________________
Int_t AliRICHParam::TotQdc(TVector2 x2,Double_t eloss)
{
// Calculates the total charge produced by the eloss in point x2 (Chamber RS).
// Returns this change parametrised in QDC channels, or 0 if the hit in the dead zone.
// eloss=0 means photons which provided for only 1 electron
// eloss > 0 for Mip
  if(Loc2Sec(x2)==kBad) return 0; //hit in the dead zone     
  Int_t iNelectrons=Int_t(eloss/IonisationPotential()); if(iNelectrons==0) iNelectrons=1;
  Double_t qdc=0;
  for(Int_t i=1;i<=iNelectrons;i++) qdc+=-Gain(x2)*TMath::Log(gRandom->Rndm());
  return Int_t(qdc);
}
//__________________________________________________________________________________________________
Double_t AliRICHParam::FracQdc(TVector2 x2,TVector pad)
{
// Calculates the charge fraction induced to given pad by the hit from the given point.
// Integrated Mathieson distribution is used.  
  TVector2 center2=Pad2Loc(pad);//gives center of requested pad
  Double_t normXmin=(x2.X()-center2.X()-PadSizeX()/2)  /AnodeCathodeGap();
  Double_t normXmax=(x2.X()-center2.X()+PadSizeX()/2)  /AnodeCathodeGap();
  Double_t normYmin=(x2.Y()-center2.Y()-PadSizeY()/2)  /AnodeCathodeGap();
  Double_t normYmax=(x2.Y()-center2.Y()+PadSizeY()/2)  /AnodeCathodeGap();
  
  if(Loc2Sec(x2)!=Pad2Sec(pad)) return 0;//requested pad does not belong to the sector of the given hit position  
  else                          return Mathieson(normXmin, normYmin, normXmax, normYmax);
}
//__________________________________________________________________________________________________
Double_t AliRICHParam::Mathieson(Double_t xMin,Double_t yMin,Double_t xMax,Double_t yMax)
{
// All arguments are parametrised according to NIM A370(1988)602-603
// Returns a charge fraction.   
  const Double_t kSqrtKx3=0.77459667;const Double_t kX2=0.962;const Double_t kX4=0.379;
  const Double_t kSqrtKy3=0.77459667;const Double_t kY2=0.962;const Double_t kY4=0.379;

  Double_t ux1=kSqrtKx3*TMath::TanH(kX2*xMin);
  Double_t ux2=kSqrtKx3*TMath::TanH(kX2*xMax);    
  Double_t uy1=kSqrtKy3*TMath::TanH(kY2*yMin);
  Double_t uy2=kSqrtKy3*TMath::TanH(kY2*yMax);
  return 4*kX4*(TMath::ATan(ux2)-TMath::ATan(ux1))*kY4*(TMath::ATan(uy2)-TMath::ATan(uy1));
}  
//__________________________________________________________________________________________________
TVector AliRICHParam::Loc2Area(TVector2 x2)
{
// Calculates the area of disintegration for a given point. It's assumed here that this points lays on anode wire.
// Area is a rectangulare set of pads defined by its left-down and right-up coners.
  TVector area(4);
  TVector pad=Loc2Pad(x2); 
  area[0]=area[2]=pad[0]; area[1]=area[3]=pad[1];//area is just a pad fired  
  if(pad[0]!=1           && pad[0]!= NpadsXsec()+1                            ) area[0]--; //left down coner X
  if(pad[1]!=1           && pad[1]!= NpadsYsec()+1 && pad[1]!= 2*NpadsYsec()+1) area[1]--; //left down coner Y 
  if(pad[0]!=NpadsXsec() && pad[0]!= NpadsX()                                 ) area[2]++; //right up coner X
  if(pad[1]!=NpadsYsec() && pad[1]!= 2*NpadsYsec() && pad[1]!= NpadsY()       ) area[3]++; //right up coner Y
  return area;          
}
//__________________________________________________________________________________________________
Bool_t AliRICHParam::IsOverTh(Int_t ,TVector ,Double_t q)
{
// Checks if the current q is over threshold and FEE will save this value to data concentrator.
  return (q>NsigmaTh()*(SigmaThMean()+(1.-2*gRandom->Rndm())*SigmaThSpread()));
}
//__________________________________________________________________________________________________
void AliRICHParam::PropogateHelix(TVector3 x0,TVector3 p0,Double_t s,TVector3 *x,TVector3 *p)
{
// Propogates the helix given by (x0,p0) in MRS to the position of interest defined by helix length s  
  const Double_t c = 0.00299792458;
  const Double_t Bz = 0.5;       //field in Tesla
  const Double_t q = 1;          //charge in electron units
  Double_t a = -c*Bz*q;
 
  Double_t rho = a/p0.Mag();
  p->SetX(p0.X()*TMath::Cos(rho*s)-p0.Y()*TMath::Sin(rho*s));
  p->SetY(p0.Y()*TMath::Cos(rho*s)+p0.X()*TMath::Sin(rho*s));
  p->SetZ(p0.Z());
  x->SetX(x0.X()+p0.X()*TMath::Sin(rho*s)/a-p0.Y()*(1-TMath::Cos(rho*s))/a);
  x->SetY(x0.Y()+p0.Y()*TMath::Sin(rho*s)/a+p0.X()*(1-TMath::Cos(rho*s))/a);
  x->SetZ(x0.Z()+p0.Z()*s/p->Mag());
}
#endif //AliRICHParam_h
