#ifndef AliRICHParam_h
#define AliRICHParam_h

#include <TObject.h>
#include <TMath.h>
#include <TVector2.h>
#include <TRandom.h>
#include <TError.h>


static const int kNCH=7;           //number of RICH chambers 
static const int kNpadsX = 144;    //number of pads along X in single chamber
static const int kNpadsY = 160;    //number of pads along Y in single chamber
static const int kBad=-101;        //useful static const to mark initial (uninitalised) values
static const int kNsectors=6;      // nb. of sectors per chamber

static const int kadc_satm  = 4096;  //dynamic range (10 bits)
static const int kCerenkov=50000050;  //??? go to something more general like TPDGCode
static const int kFeedback=50000051;  //??? go to something more general like TPDGCode


class AliRICHParam :public TObject  
{
public:
           AliRICHParam()                    {;}
  virtual ~AliRICHParam()                    {;}
  static const Int_t   NpadsX()              {return kNpadsX;}                           //pads along X in chamber
  static const Int_t   NpadsY()              {return kNpadsY;}                           //pads along Y in chamber
  static Int_t    NpadsXsec()                {return NpadsX()/3;}                        //pads along X in sector
  static Int_t    NpadsYsec()                {return NpadsY()/2;}                        //pads along Y in sector
  static Double_t DeadZone()                 {return 2.6;}                               //dead zone size in cm  
  static Double_t PadSizeX()                 {return 0.84;}                              //pad size x in cm 
  static Double_t PadSizeY()                 {return 0.8;}                               //pad size y in cm   
  static Double_t SectorSizeX()              {return NpadsX()*PadSizeX()/3;}             //sector size x in cm
  static Double_t SectorSizeY()              {return NpadsY()*PadSizeY()/2;}             //sector size y in cm 
  static Double_t PcSizeX()                  {return NpadsX()*PadSizeX()+2*DeadZone();}  //photocathode size x in cm
  static Double_t PcSizeY()                  {return NpadsY()*PadSizeY()+DeadZone();}    //photocathode size y in cm 
  static Double_t WirePitch()                {return PadSizeX()/2;}                      //distance between anode wires
  static Double_t SizeX()                    {return 132.6;}
  static Double_t SizeY()                    {return 26;}
  static Double_t SizeZ()                    {return 136.7;}                             
  static Double_t Offset()                   {return 490+1.267;}                         //distance from IP to center of chamber in cm 
  static Double_t AngleYZ()                  {return 19.5*TMath::DegToRad();}            //angle between chambers in YZ plane, rad
  static Double_t AngleXY()                  {return 20*TMath::DegToRad();}              //angle between chambers in XY plane, rad
  static Double_t AngleRot()                 {return fgAngleRot*TMath::DegToRad();}      //RICH rotation around Z, rad
  static Double_t FreonThickness()           {return 1.5;}   
  static Double_t QuartzThickness()          {return 0.5;}   
  static Double_t GapThickness()             {return 8.0;}      
  static Double_t RadiatorToPads()           {return FreonThickness()+QuartzThickness()+GapThickness();}   
  static Double_t ProximityGap()             {return 0.445;}    
  static Double_t AnodeCathodeGap()          {return 0.2;}
  static Double_t QuartzLength()             {return 133;}   
  static Double_t QuartzWidth()              {return 127.9;}
  static Double_t OuterFreonLength()         {return 133;}   
  static Double_t OuterFreonWidth()          {return 41.3;}   
  static Double_t InnerFreonLength()         {return 133;}   
  static Double_t InnerFreonWidth()          {return 41.3;}   
  static Double_t IonisationPotential()      {return 26.0e-9;}                            
  static TVector2 MathiesonDelta()           {return TVector2(5*0.18,5*0.18);}    
  static Int_t    MaxQdc()                   {return 4095;}          
  static Double_t AlphaFeedback(Int_t sec)   {HV(sec);return 0.036;}
  
  static Bool_t   IsResolveClusters()         {return fgIsResolveClusters;}  //go after resolved clusters?
  static Bool_t   IsWireSag()                 {return fgIsWireSag;}          //take wire sagita in account?
  static Int_t    HV(Int_t sector)            {
    if (sector>=1 && sector <=6)
      return fgHV[sector-1];
    else {
      ::Error("HV","Wrong sector %d",sector);
      return kBad;
    } 
  }       //high voltage for this sector
  static void     IsResolveClusters(Bool_t a) {fgIsResolveClusters=a;}  
  static void     SetWireSag(Bool_t status)   {fgIsWireSag=status;}  
  static void     SetHV(Int_t sector,Int_t hv){fgHV[sector-1]=hv;}  
  static void     SetAngleRot(Double_t rot)   {fgAngleRot =rot;}

  inline static void     Loc2Area(TVector2 x2,Int_t &padxMin,Int_t &padyMin,Int_t &padxMax,Int_t &padyMax); //
  inline static Int_t    Loc2Pad(TVector2 x2,Int_t &padx,Int_t &pady);                             //return sector and pad
  inline static TVector2 Pad2Loc(Int_t padx,Int_t pady);                                           //return center of the pad
         static Int_t    Sector(Int_t padx,Int_t pady)          {return Pad2Sec(padx,pady);}       //sector of this pad
         static Int_t    Sector(TVector2 x2)                    {int x,y;return Loc2Pad(x2,x,y);}  //sector of this point
  inline static Int_t    PadNeighbours(Int_t iPadX,Int_t iPadY,Int_t aListX[4],Int_t aListY[4]);   //number of neighbours for this pad
  inline static TVector2 ShiftToWirePos(TVector2 x2);                                              //shift to the nearest wire
  
  inline static Double_t Mathieson(Double_t lx1,Double_t lx2,Double_t ly1,Double_t ly2);           //Mathienson integral over these limits
  inline static Double_t GainSag(Double_t y,Int_t sector);                                         //gain variations in %
  inline static Double_t QdcSlope(Int_t sec);                                                      //weight of electon in QDC channels
  inline static Double_t Gain(TVector2 x2);                                                        //gain for point in ChRS 
  inline static Double_t FracQdc(TVector2 x2,Int_t padx,Int_t pady);                               //charge fraction to pad from hit
  inline static Int_t    TotQdc(TVector2 x2,Double_t eloss);                                       //total charge for hit eloss=0 for photons
  inline        Bool_t   IsOverTh(Int_t iChamber, Int_t x, Int_t y, Double_t q);                   //  
         static Int_t   NsigmaTh()                    {return fgNsigmaTh;}                         //
         static Float_t SigmaThMean()                 {return fgSigmaThMean;}      //
         static Float_t SigmaThSpread()               {return fgSigmaThSpread;}    //
                void    GenSigmaThMap();                                           //generate pedestal map
         static void    Print();                
protected:
  inline static Int_t    Loc2Sec(TVector2 &x2);             //return sector, x2->Sector RS
  inline static Int_t    Pad2Sec(Int_t &padx,Int_t &pady);  //return sector, (padx,pady)->Sector RS
  static Bool_t  fgIsWireSag;                               //is wire sagitta taken into account
  static Bool_t  fgIsResolveClusters;                       //performs declustering or not
  static Int_t   fgHV[6];                                   //HV applied to anod wires
  static Double_t fgAngleRot;                               //rotation of RICH from up postion (0,0,490)cm
  static Float_t fSigmaThMap[kNCH][kNpadsX][kNpadsY];       //sigma of the pedestal distributions for all pads
  static Int_t fgNsigmaTh;                                  //n. of sigmas to cut for zero suppression
  static Float_t fgSigmaThMean;                             //sigma threshold value
  static Float_t fgSigmaThSpread;                           //spread of sigma
  ClassDef(AliRICHParam,4)    //RICH main parameters
};
//__________________________________________________________________________________________________
Int_t AliRICHParam::PadNeighbours(Int_t iPadX,Int_t iPadY,Int_t listX[4],Int_t listY[4])
{
// Determines all the neighbouring pads for the given one. Returns total amount of these pads.
// Dead zones are taken into account.    
  Int_t nPads=0;
  if(iPadY!=NpadsY()&&iPadY!=NpadsYsec())                      {listX[nPads]=iPadX;   listY[nPads]=iPadY+1; nPads++;}       
  if(iPadX!=NpadsXsec()&&iPadX!=2*NpadsXsec()&&iPadX!=NpadsX()){listX[nPads]=iPadX+1; listY[nPads]=iPadY;   nPads++;}       
  if(iPadY!=1&&iPadY!=NpadsYsec()+1)                           {listX[nPads]=iPadX;   listY[nPads]=iPadY-1; nPads++;}      
  if(iPadX!=1&&iPadX!=NpadsXsec()+1&&iPadX!=2*NpadsXsec()+1)   {listX[nPads]=iPadX-1; listY[nPads]=iPadY;   nPads++;}

  return nPads;
}//Pad2ClosePads()
//__________________________________________________________________________________________________
Int_t AliRICHParam::Loc2Sec(TVector2 &x2)
{
// Determines sector containing the given point and trasform this point to the local system of that sector.
// Returns sector code: 1 2 3
//                      4 5 6
  Int_t sector=kBad;  
  Double_t p1=-0.5*PcSizeX();      Double_t p2=-0.5*SectorSizeX()-DeadZone();  Double_t p3=-0.5*SectorSizeX();
  Double_t p4= 0.5*SectorSizeX();  Double_t p5= 0.5*SectorSizeX()+DeadZone();  Double_t p6= 0.5*PcSizeX();
  Double_t x,y;  
  if     (x2.X()>=p1&&x2.X()<=p2)    {sector=1;x=x2.X()+0.5*PcSizeX();}
  else if(x2.X()>=p3&&x2.X()<=p4)    {sector=2;x=x2.X()+0.5*SectorSizeX();}
  else if(x2.X()>=p5&&x2.X()<=p6)    {sector=3;x=x2.X()-0.5*SectorSizeX()-DeadZone();}
  else                               {return kBad;} //in dead zone or out of chamber
  
  if     (x2.Y()>=-0.5*PcSizeY() &&x2.Y()<=-0.5*DeadZone())  {y=x2.Y()+0.5*PcSizeY();sector+=3;} //sectors 4,5,6 
  else if(x2.Y()> -0.5*DeadZone()&&x2.Y()<  0.5*DeadZone())  {return kBad;}                      //in dead zone
  else if(x2.Y()>= 0.5*DeadZone()&&x2.Y()<= 0.5*PcSizeY())   {y=x2.Y()-0.5*DeadZone();}          //sectors 1,2,3
  else                                                       {return kBad;}                      //out of chamber    
  x2.Set(x,y);
  return sector;
}//Loc2Sec(Double_t x, Double_t y)
//__________________________________________________________________________________________________
Int_t AliRICHParam::Loc2Pad(TVector2 x2,Int_t &padx,Int_t &pady)
{
// Determines pad number (padx,pady) containing the given point x2 defined the chamber RS.
// Pad count starts in lower left corner from 1,1 to 144,160 in upper right corner of a chamber.
// Returns sector number of the determined pad.      
  Int_t sector=Loc2Sec(x2);//trasforms x2 to sector reference system
  if(sector==kBad) {padx=pady=kBad; return sector;}
  
  padx=Int_t(x2.X()/PadSizeX())+1; if(padx>NpadsXsec()) padx= NpadsXsec();       
  if(sector==2||sector==5)   padx+=  NpadsXsec();     // 1 2 3
  if(sector==3||sector==6)   padx+=2*NpadsXsec();     // 4 5 6

  pady=Int_t(x2.Y()/PadSizeY())+1; if(pady>NpadsYsec()) pady= NpadsYsec();
  if(sector<4)               pady+=NpadsYsec();    
  return sector;
}
//__________________________________________________________________________________________________
Int_t AliRICHParam::Pad2Sec(Int_t &padx, Int_t &pady)
{
// Determines sector containing the given pad (padx,pady) and trasform it to the local RS of that sector.
  Int_t sector=kBad;      
  if     (padx>=1            &&padx<=NpadsXsec())      {sector=1;}
  else if(padx> NpadsXsec()  &&padx<=NpadsXsec()*2)    {sector=2;padx-=NpadsXsec();}
  else if(padx> NpadsXsec()*2&&padx<=NpadsX())         {sector=3;padx-=NpadsXsec()*2;}
  else                                                 {return kBad;}

  if     (pady>=1            &&pady<=NpadsYsec())      {return sector+3;}
  else if(pady>NpadsYsec()   &&pady<=NpadsY())         {pady-=NpadsYsec();return sector;} 
  else                                                 {return kBad;}
}//Pad2Sec()
//__________________________________________________________________________________________________
TVector2 AliRICHParam::Pad2Loc(Int_t padx,Int_t pady)
{
// Returns position of the center of the given pad (padx,pady) in local RS of the chamber    
  Int_t sector=Pad2Sec(padx,pady);//shifts to sector RS
  if(sector==kBad) return TVector2(-101,-101);  
  Double_t x,y;
  if(sector<=3)
    y=0.5*DeadZone()+pady*PadSizeY()-0.5*PadSizeY();   // 1 2 3
  else{                                                // 4 5 6
    y=-0.5*PcSizeY()+pady*PadSizeY()-0.5*PadSizeY();
  }
  if(sector==1||sector==4)
    x=-0.5*PcSizeX()+padx*PadSizeX()-0.5*PadSizeX();
  else if(sector==2||sector==5)
    x=-0.5*SectorSizeX()+padx*PadSizeX()-0.5*PadSizeX();
  else
    x= 0.5*SectorSizeX()+DeadZone()+padx*PadSizeX()-0.5*PadSizeX();
  return TVector2(x,y);
}
//__________________________________________________________________________________________________
Double_t AliRICHParam::GainSag(Double_t y,Int_t sector)
{
// Returns % of gain variation due to wire sagita.
// All cureves are parametrized per sector basis, so y must be scaled to the Sector RS.    
  if(y>0) y-=SectorSizeY()/2; else  y+=SectorSizeY()/2; 
  switch(HV(sector)){
    case 2150: return 9e-6*TMath::Power(y,4)+2e-7*TMath::Power(y,3)-0.0316*TMath::Power(y,2)-3e-4*y+25.367;//%
    case 2100: return 8e-6*TMath::Power(y,4)+2e-7*TMath::Power(y,3)-0.0283*TMath::Power(y,2)-2e-4*y+23.015;
    case 2050: return 7e-6*TMath::Power(y,4)+1e-7*TMath::Power(y,3)-0.0254*TMath::Power(y,2)-2e-4*y+20.888;
    case 2000: return 6e-6*TMath::Power(y,4)+8e-8*TMath::Power(y,3)-0.0227*TMath::Power(y,2)-1e-4*y+18.961;
    default:   return 0;
  }
}
//__________________________________________________________________________________________________
Double_t AliRICHParam::QdcSlope(Int_t sec)
{
// Returns number of QDC channels per single electron at the unknown yet ???? point for a given sector 
  switch(sec){
    case kBad: return 0;
    default:   return 27;
  }
}
//__________________________________________________________________________________________________
Double_t AliRICHParam::Gain(TVector2 x2)
{ 
//   
  if(IsWireSag()) 
    return QdcSlope(Sector(x2))*(1+GainSag(x2.Y(),Sector(x2))/100);
  else
    return QdcSlope(Sector(x2));
}
//__________________________________________________________________________________________________
Int_t AliRICHParam::TotQdc(TVector2 x2,Double_t eloss)
{
// Calculates the total charge produced by the eloss in point x2 (Chamber RS).
// Returns this change parametrised in QDC channels.
// eloss=0 means photons which provided for only 1 electron
// eloss > 0 for Mip
  if(Sector(x2)==kBad) return 0; //hit in the dead zone     
  Int_t iNelectrons=Int_t(eloss/IonisationPotential()); if(iNelectrons==0) iNelectrons=1;
  Double_t qdc=0;
  for(Int_t i=1;i<=iNelectrons;i++) qdc+=-Gain(x2)*TMath::Log(gRandom->Rndm());
  return Int_t(qdc);
}
//__________________________________________________________________________________________________
Double_t AliRICHParam::FracQdc(TVector2 x2,Int_t padx,Int_t pady)
{
// Calculates the charge fraction for a given pad (padx,pady) from the given hit point.
// Mathieson distribution integrated is used.  
  TVector2 center2=Pad2Loc(padx,pady);//gives center of requested pad
  Double_t normXmin=(x2.X()-center2.X()-PadSizeX()/2)  /AnodeCathodeGap();
  Double_t normXmax=(x2.X()-center2.X()+PadSizeX()/2)  /AnodeCathodeGap();
  Double_t normYmin=(x2.Y()-center2.Y()-PadSizeY()/2)  /AnodeCathodeGap();
  Double_t normYmax=(x2.Y()-center2.Y()+PadSizeY()/2)  /AnodeCathodeGap();
  
  if(Sector(x2)!=Sector(padx,pady)) return 0;//requested pad does not belong to the sector of given point  
  else                              return Mathieson(normXmin, normYmin, normXmax, normYmax);
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
void AliRICHParam::Loc2Area(TVector2 x2,Int_t &iPadXmin,Int_t &iPadYmin,Int_t &iPadXmax,Int_t &iPadYmax)
{
// Calculates the area of disintegration for a given point. It's assumed here that this points lays on anode wire.
// Area is a rectangulare set of pads defined by its left-down and right-up coners.
  Loc2Pad(x2-MathiesonDelta(),iPadXmin,iPadYmin);
  Loc2Pad(x2+MathiesonDelta(),iPadXmax,iPadYmax);    
}
//__________________________________________________________________________________________________
Bool_t AliRICHParam::IsOverTh(Int_t c,Int_t x,Int_t y,Double_t q)
{
// Calculate the new charge subtracting pedestal and if the current digit is over threshold
  if (c>0 && x>0 && y>0 && c<kNCH && x<kNpadsX && y<kNpadsY)
    if(q>NsigmaTh()*fSigmaThMap[c-1][x-1][y-1]) return kTRUE;
  return kFALSE;
}
//__________________________________________________________________________________________________
TVector2 AliRICHParam::ShiftToWirePos(TVector2 x2)
{
// Calculate the position of the wire nearest to the hit
  Int_t padx,pady;
  Loc2Pad(x2,padx,pady);
  Double_t x;
  TVector2 center2=Pad2Loc(padx,pady);
  if(x2.X()>center2.X()) x=center2.X()+0.5*WirePitch();
  else                   x=center2.X()-0.5*WirePitch();
  x2.Set(x,x2.Y());
  return x2;
}
#endif //AliRICHParam_h
