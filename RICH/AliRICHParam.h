#ifndef AliRICHParam_h
#define AliRICHParam_h

#include "AliRICHConst.h"
#include <TObject.h>
#include <TMath.h>
#include <TVector3.h>
#include <TRandom.h>

class AliRICHParam :public TObject  
{
public:
           AliRICHParam()                    {;}
  virtual ~AliRICHParam()                    {;}
  static const Int_t   NpadsX()             {return kNpadsX;}
  static const Int_t   NpadsY()             {return kNpadsY;}   
  static Int_t   NpadsXsec()                {return NpadsX()/3;}   
  static Int_t   NpadsYsec()                {return NpadsY()/2;}   
  static Double_t DeadZone()                 {return 2.6;}
  static Double_t PadSizeX()                 {return 0.84;}
  static Double_t PadSizeY()                 {return 0.8;}
  static Double_t SectorSizeX()              {return NpadsX()*PadSizeX()/3;}
  static Double_t SectorSizeY()              {return NpadsY()*PadSizeY()/2;}  
  static Double_t PcSizeX()                  {return NpadsX()*PadSizeX()+2*DeadZone();}
  static Double_t PcSizeY()                  {return NpadsY()*PadSizeY()+DeadZone();}
  static Double_t WirePitch()                {return PadSizeX()/2;}
  static Double_t SizeX()                    {return 132.6;}
  static Double_t SizeY()                    {return 26;}
  static Double_t SizeZ()                    {return 136.7;}   
  static Double_t Offset()                   {return 490+1.267;}  
  static Double_t AngleYZ()                  {return 19.5*TMath::DegToRad();} 
  static Double_t AngleXY()                  {return 20*TMath::DegToRad();} 
  static Double_t FreonThickness()           {return 1.5;}   
  static Double_t QuartzThickness()          {return 0.5;}   
  static Double_t GapThickness()             {return 8.0;}      
  static Double_t RadiatorToPads()           {return FreonThickness()+QuartzThickness()+GapThickness();}   
  static Double_t ProximityGapThickness()    {return 0.4;}    
  static Double_t AnodeCathodeGap()          {return 0.2;}
  static Double_t QuartzLength()             {return 133;}   
  static Double_t QuartzWidth()              {return 127.9;}
  static Double_t OuterFreonLength()         {return 133;}   
  static Double_t OuterFreonWidth()          {return 41.3;}   
  static Double_t InnerFreonLength()         {return 133;}   
  static Double_t InnerFreonWidth()          {return 41.3;}   
  static Double_t IonisationPotential()      {return 26.0e-9;}                            
  static Double_t MathiesonDeltaX()          {return 5*0.18;}    
  static Double_t MathiesonDeltaY()          {return 5*0.18;}    
  static Int_t    MaxQdc()                   {return 4095;}          
  static Double_t QdcSlope(Int_t sec)        {HV(sec);return 27;}
  static Double_t AlphaFeedback(Int_t sec)   {HV(sec);return 0.036;}
  
  static Bool_t  IsResolveClusters()         {return fgIsResolveClusters;}  
  static Bool_t   IsWireSag()                {return fgIsWireSag;}
  static Int_t    HV(Int_t)                  {return fgHV;}
  static Double_t AngleRot()                 {return fgAngleRot*TMath::DegToRad();} 
    static void  SetResolveClusters(Bool_t a){fgIsResolveClusters=a;}  
    static void  SetWireSag(Bool_t status)   {fgIsWireSag=status;}  
    static void  SetHV(Int_t hv)             {fgHV       =hv;}  
    static void  SetAngleRot(Double_t rot)   {fgAngleRot =rot;}

  inline static Double_t Mathieson(Double_t lx1,Double_t lx2,Double_t ly1,Double_t ly2);   
  inline static void    Loc2Area(TVector3 hitX3,Int_t &padxMin,Int_t &padyMin,Int_t &padxMax,Int_t &padyMax);
  inline static Int_t   PadNeighbours(Int_t iPadX,Int_t iPadY,Int_t aListX[4],Int_t aListY[4]);
  inline static Int_t   Loc2Pad(Double_t x,Double_t y,Int_t &padx,Int_t &pady); 
  inline static void    Pad2Loc(Int_t padx,Int_t pady,Double_t &x,Double_t &y);  
  inline static Double_t GainVariation(Double_t y,Int_t sector);       
  inline static Int_t   Loc2TotQdc(TVector3 locX3,Double_t eloss,Int_t iPid, Int_t &sector);
  inline static Double_t Loc2PadFrac(TVector3 locX3,Int_t padx,Int_t pady);
  
  inline static Int_t   Loc2Sec(Double_t &x,Double_t &y); 
  inline static Int_t   Pad2Sec(Int_t &padx,Int_t &pady); 
  inline Bool_t IsOverTh(Int_t iChamber, Int_t x, Int_t y, Double_t q);
  static Int_t NsigmaTh() {return fgNsigmaTh;}
  static Float_t SigmaThMean() {return fgSigmaThMean;}
  static Float_t SigmaThSpread() {return fgSigmaThSpread;}
  void GenSigmaThMap();
protected:
  static Bool_t  fgIsWireSag;                           //is wire sagitta taken into account
  static Bool_t  fgIsResolveClusters;                   //performs declustering or not
  static Int_t   fgHV;                                  //HV applied to anod wires
  static Double_t fgAngleRot;                           //rotation of RICH from up postion (0,0,490)cm
  Float_t fSigmaThMap[kNCH][kNpadsX][kNpadsY];          // sigma of the pedestal distributions for all pads
  static Int_t fgNsigmaTh;                              // n. of sigmas to cut for zero suppression
  static Float_t fgSigmaThMean;                         // sigma threshold value
  static Float_t fgSigmaThSpread;                       // spread of sigma
  ClassDef(AliRICHParam,4)    //RICH main parameters
};
//__________________________________________________________________________________________________
Int_t AliRICHParam::PadNeighbours(Int_t iPadX,Int_t iPadY,Int_t listX[4],Int_t listY[4])
{
  Int_t nPads=0;
  if(iPadY<NpadsY()){listX[nPads]=iPadX;   listY[nPads]=iPadY+1; nPads++;}       
  if(iPadX<NpadsX()){listX[nPads]=iPadX+1; listY[nPads]=iPadY;   nPads++;}       
  if(iPadY>1)       {listX[nPads]=iPadX;   listY[nPads]=iPadY-1; nPads++;}      
  if(iPadX>1)       {listX[nPads]=iPadX-1; listY[nPads]=iPadY;   nPads++;}       
  return nPads;
}//Pad2ClosePads()
//__________________________________________________________________________________________________
Int_t AliRICHParam::Loc2Sec(Double_t &x,Double_t &y)
{//Determines sector for a given hit (x,y) and trasform this point to the local system of that sector.
  Int_t sector=kBad;  
  Double_t x1=-PcSizeX()/2;      Double_t x2=-SectorSizeX()/2-DeadZone();  Double_t x3=-SectorSizeX()/2;
  Double_t x4= SectorSizeX()/2;  Double_t x5= SectorSizeX()/2+DeadZone();  Double_t x6= PcSizeX()/2;

  if     (x>=x1&&x<=x2)    {sector=1;x+=PcSizeX()/2;}
  else if(x>=x3&&x<=x4)    {sector=2;x+=SectorSizeX()/2;}
  else if(x>=x5&&x<=x6)    {sector=3;x-=SectorSizeX()/2+DeadZone();}
  else                     {return kBad;} //in dead zone

  if     (y>=-PcSizeY()/2   &&y<=-DeadZone()/2)  {y+=PcSizeY()/2;  return sector;}
  else if(y> -DeadZone()/2  &&y<  DeadZone()/2)  {return kBad;} //in dead zone
  else if(y>= DeadZone()/2  &&y<= PcSizeY()/2)   {y-=DeadZone()/2; return sector+3;}
  else                                           {return kBad;}
}//Loc2Sec(Double_t x, Double_t y)
//__________________________________________________________________________________________________
Int_t AliRICHParam::Pad2Sec(Int_t &padx, Int_t &pady)
{//Determines sector for a given pad (padx,pady) and trasform this point to the local system of that sector.
  Int_t sector=kBad;      
  if     (padx>=1            &&padx<=NpadsXsec())      {sector=1;}
  else if(padx> NpadsXsec()  &&padx<=NpadsXsec()*2)    {sector=2;padx-=NpadsXsec();}
  else if(padx> NpadsXsec()*2&&padx<=NpadsX())         {sector=3;padx-=NpadsXsec()*2;}
  else                                                 {return kBad;}

  if     (pady>=1         &&pady<= NpadsYsec())     {return sector;}
  else if(pady>NpadsYsec()&&pady<= NpadsY())        {pady-=NpadsYsec();return sector+3;} 
  else                                              {return kBad;}
}//Pad2Sec()
//__________________________________________________________________________________________________
Int_t AliRICHParam::Loc2Pad(Double_t x, Double_t y, Int_t &padx, Int_t &pady)
{//returns pad numbers (iPadX,iPadY) for given point in local coordinates (x,y) 
 //count starts in lower left corner from 1,1 to 144,180
  
  padx=pady=kBad;
  Int_t sector=Loc2Sec(x,y);
  if(sector==kBad) return sector;
  
  padx=Int_t(x/PadSizeX())+1; 
  if(padx>NpadsXsec())            padx= NpadsXsec();
  if(sector==2||sector==5)       padx+=NpadsXsec();
  else if(sector==3||sector==6)  padx+=NpadsXsec()*2;
  
  pady=Int_t(y/PadSizeY())+1;
  if(pady>NpadsYsec())            padx= NpadsYsec();
  if(sector>0)                    pady+=NpadsYsec();    

  return sector;
}//Loc2Pad()
//__________________________________________________________________________________________________
void AliRICHParam::Pad2Loc(Int_t padx,Int_t pady,Double_t &x,Double_t &y)
{
  Int_t sector=Pad2Sec(padx,pady);  
  if(sector>3)
    y=0.5*DeadZone()+pady*PadSizeY()-0.5*PadSizeY();
  else{
    y=-0.5*PcSizeY()+pady*PadSizeY()-0.5*PadSizeY();
  }
  if(sector==1||sector==4)
    x=-0.5*PcSizeX()+padx*PadSizeX()-0.5*PadSizeX();
  else if(sector==2||sector==5)
    x=-0.5*SectorSizeX()+padx*PadSizeX()-0.5*PadSizeX();
  else
    x= 0.5*SectorSizeX()+DeadZone()+padx*PadSizeX()-0.5*PadSizeX();
  return;
}//Pad2Loc()
//__________________________________________________________________________________________________
Double_t AliRICHParam::GainVariation(Double_t y,Int_t sector)
{
  if(IsWireSag()){
    if(y>0) y-=SectorSizeY()/2; else  y+=SectorSizeY()/2; 
    switch(HV(sector)){
      case 2150:
      default:  
        return 9e-6*TMath::Power(y,4)+2e-7*TMath::Power(y,3)-0.0316*TMath::Power(y,2)-3e-4*y+25.367;//%
    }
  }else
    return 0;
}
//__________________________________________________________________________________________________
Int_t AliRICHParam::Loc2TotQdc(TVector3 x3,Double_t eloss,Int_t iPid,Int_t &sector)
{//calculates the total charge produced by the hit given in local refenrence system
  Double_t x=x3.X(),y=x3.Y();
  
  sector=Loc2Sec(x,y);
  
  Double_t gain=QdcSlope(sector)*(1+GainVariation(x3.Y(),sector)/100);

  
  if(iPid>50000){//it's photon => 1 electron
    return Int_t(gain*-TMath::Log(gRandom->Rndm()));
  }else{//it's MIP  
    Int_t iNelectrons=Int_t(eloss/IonisationPotential());
    if(iNelectrons==0) return 0;
    Double_t qdc=0;
    for(Int_t i=1;i<=iNelectrons;i++) qdc+=gain*-TMath::Log(gRandom->Rndm());
    return Int_t(qdc);
  }
}
//__________________________________________________________________________________________________
Double_t AliRICHParam::Loc2PadFrac(TVector3 hitX3,Int_t padx,Int_t pady)
{//
  Double_t padXcenter=0,padYcenter=0;  Pad2Loc(padx,pady,padXcenter,padYcenter);  

  //correction to the position of the nearest wire
  
  Double_t normXmin=(hitX3.X()-padXcenter-PadSizeX()/2)  /AnodeCathodeGap();
  Double_t normXmax=(hitX3.X()-padXcenter+PadSizeX()/2)  /AnodeCathodeGap();
  Double_t normYmin=(hitX3.Y()-padYcenter-PadSizeY()/2)  /AnodeCathodeGap();
  Double_t normYmax=(hitX3.Y()-padYcenter+PadSizeY()/2)  /AnodeCathodeGap();
  
  return Mathieson(normXmin,normYmin,normXmax,normYmax);
}//Loc2PadQdc()
//__________________________________________________________________________________________________
Double_t AliRICHParam::Mathieson(Double_t xMin,Double_t yMin,Double_t xMax,Double_t yMax)
{//see NIM A370(1988)602-603 
  const Double_t kSqrtKx3=0.77459667;const Double_t kX2=0.962;const Double_t kX4=0.379;
  const Double_t kSqrtKy3=0.77459667;const Double_t kY2=0.962;const Double_t kY4=0.379;

  Double_t ux1=kSqrtKx3*TMath::TanH(kX2*xMin);
  Double_t ux2=kSqrtKx3*TMath::TanH(kX2*xMax);    
  Double_t uy1=kSqrtKy3*TMath::TanH(kY2*yMin);
  Double_t uy2=kSqrtKy3*TMath::TanH(kY2*yMax);
  return 4*kX4*(TMath::ATan(ux2)-TMath::ATan(ux1))*kY4*(TMath::ATan(uy2)-TMath::ATan(uy1));
}  
//__________________________________________________________________________________________________
void AliRICHParam::Loc2Area(TVector3 hitX3,Int_t &iPadXmin,Int_t &iPadYmin,Int_t &iPadXmax,Int_t &iPadYmax)
{//calculates the area of disintegration for a given hit. Area is a rectangulare set pf pads
 //defined by its left-down and right-up coners
  //  hitX3.SetX(Shift2NearestWire(hitX3.X());
  Loc2Pad(hitX3.X()-MathiesonDeltaX(),hitX3.Y()-MathiesonDeltaY(),iPadXmin,iPadYmin);   
  Loc2Pad(hitX3.X()+MathiesonDeltaX(),hitX3.Y()+MathiesonDeltaY(),iPadXmax,iPadYmax);     
}//
//__________________________________________________________________________________________________
Bool_t AliRICHParam::IsOverTh(Int_t iChamber, Int_t x, Int_t y, Double_t q)
{// Calculate the new charge subtracting pedestal and if the current digit is over threshold
  if(q>NsigmaTh()*fSigmaThMap[iChamber-1][x-1][y-1]) return kTRUE; else return kFALSE;
}//
#endif //AliRICHParam_h
