#ifndef AliRICHParam_h
#define AliRICHParam_h

#include <TError.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TError.h>
#include <TObjArray.h>
#include <AliLog.h>
#include <TClass.h>


static const int kNchambers=7;     //number of RICH chambers 
static const int kNpadsX = 160;    //number of pads along X in single chamber
static const int kNpadsY = 144;    //number of pads along Y in single chamber
static const int kBad=-101;        //useful static const to mark initial (uninitalised) values
static const int kNsectors=6;      //number of sectors per chamber

static const int kCerenkov=50000050;  //??? go to something more general like TPDGCode
static const int kFeedback=50000051;  //??? go to something more general like TPDGCode

class AliRICHChamber;

class AliRICHParam :public TObject  
{
public:
                  AliRICHParam():TObject(),fpChambers(0)  {CreateChambers();}
  virtual        ~AliRICHParam()                          {delete fpChambers;}
         void     CreateChambers();
  AliRICHChamber* C(Int_t i)                 {return (AliRICHChamber*)fpChambers->UncheckedAt(i-1);}      //returns pointer to chamber i
  static Int_t    NpadsX()                   {return kNpadsX;}                           //pads along X in chamber
  static Int_t    NpadsY()                   {return kNpadsY;}                           //pads along Y in chamber
  static Int_t    NpadsXsec()                {return NpadsX()/2;}                        //pads along X in sector
  static Int_t    NpadsYsec()                {return NpadsY()/3;}                        //pads along Y in sector
  static Double_t DeadZone()                 {return 2.6;}                               //dead zone size in cm  
  static Double_t PadSizeX()                 {return 0.8;}                               //pad size x in cm 
  static Double_t PadSizeY()                 {return 0.84;}                              //pad size y in cm   
  
  static Double_t SectorSizeX()              {return NpadsX()*PadSizeX()/2;}             //sector size x in cm
  static Double_t SectorSizeY()              {return NpadsY()*PadSizeY()/3;}             //sector size y in cm 
  static Double_t PcSizeX()                  {return NpadsX()*PadSizeX()+DeadZone();}    //PC size x, cm
  static Double_t PcSizeY()                  {return NpadsY()*PadSizeY()+2*DeadZone();}  //PC size y, cm
   
  static Double_t Zfreon()                   {return 1.5;}                               //freon thinkness, cm
  static Double_t Zwin()                     {return 0.5;}                               //radiator quartz window, cm   
  static Double_t Pc2Win()                   {return 8.0;}                               //cm between CsI PC and radiator quartz window
  static Double_t Pc2Coll()                  {return 7.0;}                               //cm between CsI PC and third wire grid (collection wires)     
  static Double_t Pc2Anod()                  {return 0.204;}                             //cm between CsI PC and first wire grid (anod wires)     
  static Double_t Pc2Cath()                  {return 0.445;}                             //cm between CsI PC and second wire grid (cathode wires)
  static Double_t Freon2Pc()                 {return Zfreon()+Zwin()+Pc2Win();}          //cm between CsI PC and entrance to freon
  static Double_t PitchAnod()                {return PadSizeY()/2;}                      //cm between anode wires
  static Double_t PitchCath()                {return PadSizeY()/4;}                      //cm between cathode wires
  static Double_t PitchColl()                {return 0.5;}                               //cm between collection wires
  
  static Double_t IonisationPotential()      {return 26.0e-9;}                            //for CH4 in GeV taken from ????
  static TVector2 MathiesonDelta()           {return TVector2(5*0.18,5*0.18);}            //area of 5 sigmas of Mathieson distribution (cm)
  static Int_t    MaxQdc()                   {return 4095;}                               //QDC number of channels          
    
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
  static Double_t IndOfRefC6F14(Double_t eV)  {return eV*0.0172+1.177;}          // eV = photon energy in eV
  static Double_t IndOfRefSiO2(Double_t eV)   {Double_t e1=10.666,e2=18.125,f1=46.411,f2= 228.71;
                                     return TMath::Sqrt(1.+f1/(e1*e1-TMath::Power(eV,2))+f2/(e2*e2-TMath::Power(eV,2)));}//TDR p.35
  static Double_t IndOfRefCH4()               {return 1.000444;}

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
                
  inline static Int_t    Loc2Sec(TVector2 &x2);             //return sector, x2->Sector RS
  inline static Int_t    Pad2Sec(const TVector &pad);              //return sector
         static Bool_t   IsAccepted(const TVector2 &x2) {return ( x2.X()>=0 && x2.X()<=PcSizeX() && x2.Y()>=0 && x2.Y()<=PcSizeY() ) 
? kTRUE:kFALSE;}
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
  else                                    {return kBad;}
  
  if     (v2.Y() >= y0 && v2.Y() <= y1 )  {}                                  //sectors 1 or 2 
  else if(v2.Y() >= y2 && v2.Y() <= y3 )  {sector+=2; y=v2.Y()-y2;}           //sectors 3 or 4
  else if(v2.Y() >= y4 && v2.Y() <= y5 )  {sector+=4; y=v2.Y()-y4;}           //sectors 5 or 6
  else                                    {return kBad;}
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
Int_t AliRICHParam::Pad2Sec(const TVector &pad)
{
// Determines sector containing the given pad.
  Int_t sector=kBad;      
  if     (pad[0] >= 1           && pad[0] <=   NpadsXsec() )    {sector=1;}
  else if(pad[0] >  NpadsXsec() && pad[0] <=   NpadsX()    )    {sector=2;} 
  else                                                         AliDebugClass(1,Form("Wrong pad (%3.0f,%3.0f)",pad[0],pad[1]));
    
  if     (pad[1] >= 1             && pad[1] <=   NpadsYsec() )    {}
  else if(pad[1] >  NpadsYsec()   && pad[1] <= 2*NpadsYsec() )    {sector+=2;}
  else if(pad[1] >  2*NpadsYsec() && pad[1] <=   NpadsY()    )    {sector+=4;}
  else                                                         AliDebugClass(1,Form("Wrong pad (%3.0f,%3.0f)",pad[0],pad[1]));

  return sector;
}//Pad2Sec()
//__________________________________________________________________________________________________
TVector2 AliRICHParam::Pad2Loc(TVector pad)
{
// Returns position of the center of the given pad in local system of the chamber (cm)    
// y ^  5 6
//   |  3 4        sector numbers
//   |  1 2
//    -------> x  
  Double_t x=kBad,y=kBad;
  if(pad[0] > 0 && pad[0] <= NpadsXsec())//it's 1 or 3 or 5
    x=(pad[0]-0.5)*PadSizeX();
  else if(pad[0] > NpadsXsec() && pad[0] <= NpadsX())//it's 2 or 4 or 6
    x=(pad[0]-0.5)*PadSizeX()+DeadZone();
  else
    AliDebugClass(1,Form("Wrong pad (%3.0f,%3.0f)",pad[0],pad[1]));
  
  if(pad[1] > 0 && pad[1] <= NpadsYsec())//it's 1 or 2
    y=(pad[1]-0.5)*PadSizeY();
  else if(pad[1] > NpadsYsec() && pad[1] <= 2*NpadsYsec())//it's 3 or 4
    y=(pad[1]-0.5)*PadSizeY()+DeadZone();
  else if(pad[1] > 2*NpadsYsec() && pad[1]<= NpadsY())//it's 5 or 6
    y=(pad[1]-0.5)*PadSizeY()+2*DeadZone();
  else
    AliDebugClass(1,Form("Wrong pad (%3.0f,%3.0f)",pad[0],pad[1]));
    
  return TVector2(x,y);
}
//__________________________________________________________________________________________________
Double_t AliRICHParam::GainSag(Double_t x,Int_t sector)
{
// Returns % of gain variation due to wire sagita.
// All curves are parametrized as per sector basis, so x must be apriory transformed to the Sector RS.    
// Here x is a distance along wires.  
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
// eloss=0 means photon which produces 1 electron only eloss > 0 for Mip
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
  Double_t normXmin=(x2.X()-center2.X()-PadSizeX()/2)  /Pc2Cath();//parametrise for Mathienson
  Double_t normXmax=(x2.X()-center2.X()+PadSizeX()/2)  /Pc2Cath();
  Double_t normYmin=(x2.Y()-center2.Y()-PadSizeY()/2)  /Pc2Cath();
  Double_t normYmax=(x2.Y()-center2.Y()+PadSizeY()/2)  /Pc2Cath();
 
//requested pad might not belong to the sector of the given hit position, hence the check:
  return (Loc2Sec(x2)!=Pad2Sec(pad)) ? 0:Mathieson(normXmin, normYmin, normXmax, normYmax);
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
#endif //AliRICHParam_h
