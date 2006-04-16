#ifndef AliRICHParam_h
#define AliRICHParam_h

#include <TNamed.h>       //base class
#include <TMath.h>        //QdcTot()  
#include <TVector.h>      //old style
#include <TVector2.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TClonesArray.h>  //Hit2SDigs()
#include <AliLog.h>
#include <TGeoMatrix.h>    //Mars2Lors() Lors2Mars() 
#include <TF1.h>           //fields
#include <TF2.h>           //fields
#include "AliRICHDigit.h"  //Hit2Sdigs()
#include <TGeoManager.h>   //Instance()

static const int kNchambers=7;     //number of RICH chambers 
static const int kNpadsX = 160;    //number of pads along X in single chamber
static const int kNpadsY = 144;    //number of pads along Y in single chamber
static const int kNsectors=6;      //number of sectors per chamber

static const int kCerenkov=50000050;  //??? go to something more general like TPDGCode
static const int kFeedback=50000051;  //??? go to something more general like TPDGCode

// Class providing all the needed parametrised information
// to construct the geometry, to define segmentation and to provide response model
// In future will also provide all the staff needed for alignment and calibration


class AliRICHParam :public TNamed  
{
public:
//ctor&dtor    
  virtual        ~AliRICHParam()                                    {delete fIdxC6F14;fgInstance=0;}
//test methodes  
         void     Print(Option_t *opt="") const;                                         //print current parametrization
  static void     DrawAxis();
  static void     DrawSectors();
//flags staff         
  static inline AliRICHParam* Instance();                                //pointer to AliRICHParam singleton
  static        Int_t      Stack(Int_t evt=-1,Int_t tid=-1);              //Print stack info for event and tid
  static        Int_t      StackCount(Int_t pid,Int_t evt);               //Counts stack particles of given sort in given event  
  static inline Double_t   ErrLoc                  (Double_t thetaC,Double_t phiC,Double_t thetaT,Double_t phiT,Double_t beta);
  static inline Double_t   ErrGeom                 (Double_t thetaC,Double_t phiC,Double_t thetaT,Double_t phiT,Double_t beta);
  static inline Double_t   ErrCrom                 (Double_t thetaC,Double_t phiC,Double_t thetaT,Double_t phiT,Double_t beta);
  static inline Double_t   SigmaSinglePhotonFormula(Double_t thetaC,Double_t phiC,Double_t thetaT,Double_t phiT,Double_t beta);
//Geometrical properties  
  static        Int_t      NpadsX      ()   {return kNpadsX;}                           //number of pads along X in chamber
  static        Int_t      NpadsY      ()   {return kNpadsY;}                           //number of pads along Y in chamber
  static        Int_t      NpadsXsec   ()   {return NpadsX()/2;}                        //number of pads along X in sector
  static        Int_t      NpadsYsec   ()   {return NpadsY()/3;}                        //number of pads along Y in sector
  
  static        Double_t   AnodPitch   ()   {return PadSizeY()/2;}                      //cm between anode wires
  static        Double_t   AnodZ       ()   {return 7.806;}                             //Z positon of anod plane in LORS of the chamber, [cm]
  static        Double_t   CathPitch   ()   {return PadSizeY()/4;}                      //dist between cathode wires [cm]
  static        Double_t   CollPitch   ()   {return 0.5;}                               //dist between collection wires [cm]
  static        Double_t   DeadZone    ()   {return 2.6;}                               //dead zone thickness [cm]  
  static        Double_t   PadSizeX    ()   {return 0.8;}                               //pad size x [cm] 
  static        Double_t   PadSizeY    ()   {return 0.84;}                              //pad size y [cm]   
  static        Double_t   PcSizeX     ()   {return NpadsX()*PadSizeX()+DeadZone();}    //PC size x [cm]
  static        Double_t   PcSizeY     ()   {return NpadsY()*PadSizeY()+2*DeadZone();}  //PC size y [cm]
  static        Double_t   Pc2Cath     ()   {return 0.445;}                             //dist between PC entrance plane and cathode wires plane [cm]
  static        Double_t   Pc2Win      ()   {return PcZ();}                             //dist between PC entrance plane and window exit plane [cm]
  static        Double_t   PcZ         ()   {return 8.0;  }                             //Z positon of PC entrance plane in LORS of the chamber [cm]
  static        Double_t   RadThick    ()   {return 1.5;}                               //radiator thickness [cm]
  static        Double_t   RadZ        ()   {return -2.0; }                             //Z positon of radiator entrance plane in LORS of the chamber [cm]
  static        Double_t   SecSizeX    ()   {return NpadsX()*PadSizeX()/2;}             //sector size x [cm]
  static        Double_t   SecSizeY    ()   {return NpadsY()*PadSizeY()/3;}             //sector size y [cm ]
  static        Double_t   WinThick    ()   {return 0.5;}                               //radiator window thickness [cm]   
  
  
//trasformation methodes
         inline TVector3   Lors2Mars     (Int_t c,Double_t x,Double_t y,Int_t p=kPc); //LORS->MARS transform of point [cm] for chamber c and plane p
         inline TVector3   Lors2MarsVec  (Int_t c,const TVector3 &p                ); //LORS->MARS transform of vector for chamber c
         inline TVector2   Mars2Lors     (Int_t c,const TVector3 &x    ,Int_t p=kPc); //MARS->LORS transform of point [cm] for chamber c and plane p    
         inline TVector3   Mars2LorsVec  (Int_t c,const TVector3 &p                ); //MARS->LORS transform of vector for chamber c
  
  static inline TVector3   Lors2MarsOld  (Int_t c,Double_t x,Double_t y,Int_t p); //LORS->MARS transform of position (cm) for chamber c and plane p
  static inline TVector2   Mars2LorsOld  (Int_t c,const TVector3 &x,Int_t p    ); //MARS->LORS transform of position (cm) for chamber c and plane p    
  static inline TVector3   Center        (Int_t c,Int_t p                      ); //Center of plane p of chamber c in MARS (cm)
  static inline TVector3   Norm          (Int_t c                              ); //Norm vector to the chamber c in MARS (cm)
  static inline TGeoMatrix*Matrix        (Int_t iCh, Int_t iPlane              ); //TGeoMatrix for the given chamber plain
  
  static        Int_t      Pad2Cha       (Int_t pad                      ){return pad/100000000;                     }//abs pad -> chamber
  static        Int_t      Pad2Sec       (Int_t pad                      ){return pad%100000000/1000000;             }//abs pad -> sector
  static        Int_t      Pad2PadX      (Int_t pad                      ){return pad%1000000/1000;                  }//abs pad -> pad x 
  static        Int_t      Pad2PadY      (Int_t pad                      ){return pad%1000000%100;                   }//abs pad -> pad y
  static        Int_t      PadAbs        (Int_t c,Int_t s,Int_t x,Int_t y){return 100000000*c+1000000*s+1000*x+y;    }//(c,s,x,y) -> abs pad
  static inline TVector2   Pad2Loc       (Int_t pad                      );                                           //abs pad ->LORS
  static inline TVector2   Pad2Loc       (TVector pad                    );                                           //pad  -> LORS returns center of the pad
  static        TVector2   Pad2Loc       (Int_t x,Int_t y                ){TVector pad(2);pad[0]=x;pad[1]=y;return Pad2Loc(pad);}//return center of the pad (x,y)
  static inline TVector    Loc2Area      (const TVector2 &x2             );                                           //pads area affected by hit x2. area is LeftDown-RightUp pad numbers
  static inline Int_t      Loc2Sec       (const TVector2 &x2             );                                           //LORS -> sector
  static        Int_t      Loc2Sec       (Double_t x,Double_t y          ){return Loc2Sec(TVector2(x,y));}            //LORS -> sector
  static inline TVector    Loc2Pad       (const TVector2 &x2             );                                           //LORS -> pad
  static        TVector    Loc2Pad       (Double_t x,Double_t y          ){return Loc2Pad(TVector2(x,y));}            //LORS -> pad
  static inline Int_t      Pad2Sec       (const TVector &pad             );                                           //pad  -> sector
  static inline Int_t      PadNeighbours (Int_t iPadX,Int_t iPadY,Int_t aListX[4],Int_t aListY[4]);                   //pad -> list of it neighbours
  static        Bool_t     IsAccepted    (const TVector2 &x2             ){return ( x2.X()>=0 && x2.X()<=PcSizeX() && x2.Y()>=0 && x2.Y()<=PcSizeY() ) ? kTRUE:kFALSE;}
//optical properties methodes  
  static        Float_t    EckovMean     (                               ){return 6.766e-9;}                          //mean Ckov energy according to the total trasmission curve
  static        Float_t    EckovMin      (                               ){return 5.5e-9;}                            //min photon energy [GeV] defined in optical curves
  static        Float_t    EckovMax      (                               ){return 8.5e-9;}                            //min photon energy [GeV] defined in optical curves
  
                Float_t    IdxC6F14      (Float_t gev                    ){return fIdxC6F14->Eval(gev,fIdxC6F14->GetUniqueID());}  //n=f(Eckov) [GeV] for C6H14 used as radiator
  static        Float_t    IdxSiO2       (Float_t gev                    ){return TMath::Sqrt(1+46.411/(10.666*10.666-gev*gev*1e18)+228.71/(18.125*18.125-gev*gev*1e18));} //n=f(Eckov) [GeV] for SiO2 used as window TDR p.35
  static        Float_t    IdxCH4        (Float_t gev                    ){return 1+0.12489e-6/(2.62e-4 - TMath::Power(1239.84e-9/gev,-2));}              //n=f(Eckov) [GeV] for CF4 
  static        Float_t    AbsCH4        (Float_t gev                    );                                                                          //abs len=f(Eckov) [GeV] for CF4 
  
                void       CdbRead   (Int_t run,Int_t version           );                                           //read all calibration information for requested run
  
  static Double_t IonisationPotential()      {return 26.0e-9;}                            //for CH4 in GeV taken from ????
  static TVector2 MathiesonDelta()           {return TVector2(5*0.18,5*0.18);}            //area of 5 sigmas of Mathieson distribution (cm)
  static Int_t    MaxQdc()                   {return 4095;}                               //QDC number of channels          

  
  static Int_t    QthMIP()                   {return 100;}
  static Double_t DmatchMIP()                {return 1.;}
  static Double_t PmodMax()                  {return 6.5;}
  static Int_t    HV(Int_t sector)           {if (sector>=1 && sector <=6) return 2050;  else return -1;} //high voltage for this sector
//charge response methodes  
  inline static Double_t Mathieson(Double_t x1,Double_t x2,Double_t y1,Double_t y2);               //Mathienson integral over given limits
  
  inline static Double_t GainSag(Double_t x,Int_t sector);                                         //gain variations in %
         static Double_t Gain(const TVector2 &x2,Bool_t isSag=kTRUE){//gives chamber gain in terms of QDC channels for given point in local ref system
                          if(isSag) return QdcSlope(Loc2Sec(x2))*(1+GainSag(x2.X(),Loc2Sec(x2))/100);
                          else      return QdcSlope(Loc2Sec(x2));}
  inline static Double_t FracQdc(const TVector2 &x2,const TVector &pad);                           //charge fraction to pad from hit
  inline static Int_t    TotQdc(TVector2 x2,Double_t e    );                                       //total charge for Eloss (GeV) 0 for photons
         static Double_t QdcSlope(Int_t sec){switch(sec){case -1: return 0;  default:   return 33;}} //weight of electon in QDC channels
  
  static inline Int_t    Lors2Pad      (Double_t x,Double_t y                         );                                            //LORS (x,y) [cm] -> abs pad number        
  static        Double_t IonPot        (                                              ){return 26.0e-9;}                            //for CH4 in GeV taken from ????
  static inline Int_t    QdcTot        (Int_t iPad,Double_t e                         );                                            //total QDC generated by Eloss or Etot [GeV]
  static inline Double_t QdcSag        (Int_t iPad                                    );                                            //mean QDC variation due to sagita [0,1]
  static        Double_t QdcEle        (Int_t iPad,Bool_t isSag=kTRUE                 ){return isSag?33*(1+QdcSag(iPad)):33;}       //mean QDC per electron
  static inline Int_t    Hit2SDigs     (Int_t iPad,  Double_t e,TClonesArray* pSDigLst);                                            //hit->sdigits, returns Qtot
  static inline Int_t    Hit2SDigs     (TVector2 hit,Double_t e,TClonesArray* pSDigLst);                                            //hit->sdigits, returns Qtot, old style
  static        void     TestHit2SDigs (Double_t x,Double_t y,Double_t e,Bool_t isNew=kFALSE);                                      //test hit->sdigits
  
  inline static Bool_t   IsOverTh(Int_t c,TVector pad,Double_t q);                                 //is QDC of the pad registered by FEE  
         static Int_t    NsigmaTh()                    {return 4;}                        //
         static Float_t  SigmaThMean()                 {return 1.132;}                    //QDC electronic noise mean
         static Float_t  SigmaThSpread()               {return 0.035;}                    //QDC electronic noise width
                
         static Double_t CogCorr(Double_t x) {return 3.31267e-2*TMath::Sin(2*TMath::Pi()/PadSizeX()*x) //correction of cluster CoG due to sinoidal
                                                    -2.66575e-3*TMath::Sin(4*TMath::Pi()/PadSizeX()*x)
                                                    +2.80553e-3*TMath::Sin(6*TMath::Pi()/PadSizeX()*x)+0.0070;}
       
         TVector3 ForwardTracing(TVector3 entranceTrackPoint,TVector3 vectorTrack, Double_t thetaC, Double_t phiC); //it traces foward a photon from Emission Point to PC
  static TVector3 PlaneIntersect(const TVector3 &lineDir,const TVector3 &linePoint,const TVector3 &planeNorm,const TVector3 &planePoint); //intersection between line and plane
  static Double_t SnellAngle(Float_t n1, Float_t n2, Float_t theta1);                              // Snell law
  static void     AnglesInDRS(Double_t trackTheta,Double_t trackPhi,Double_t thetaCerenkov,Double_t phiCerenkov,Double_t &tout,Double_t &pout);//It finds photon angles in 
  static        Double_t AlphaFeedback(Int_t c,Int_t s) {c++;s++; return 0.02;} //for sector s of chamber c
//test part  
  static void     Test()                            {TestSeg();TestTrans();TestResp();}  //test all groups of methodes
  static void     TestResp();                                                            //test the response group of methodes
  static void     TestSeg();                                                             //test the segmentation group of methodes
  static void     TestTrans();                                                           //test the transform group of methodes

  static Double_t fgMass[5];                                // mass array
  enum EPlaneId {kCenter,kPc,kRad,kAnod,kNch=7};            //4 planes in chamber and total number of chambers
protected:
         AliRICHParam();             //default ctor is protected to enforce it to be singleton
  static AliRICHParam *fgInstance;   //static pointer  to instance of AliRICHParam singleton
  TF2         *fIdxC6F14;            //n=f(Ephot,T) [GeV] for radiator freon   C6F14
  TGeoHMatrix *fMatrix[kNchambers];  //poiners to matrices defining RICH chambers rotations-translations
  ClassDef(AliRICHParam,0)           //RICH main parameters class
};

AliRICHParam* AliRICHParam::Instance()
{
// Return pointer to the AliRICHParam singleton. 
// Arguments: none
//   Returns: pointer to the instance of AliRICHParam or 0 if no geometry       
  if(!fgInstance&&gGeoManager) new AliRICHParam; 
  else if(!gGeoManager)                        Printf("No geometry imported");
  return fgInstance;  
}//Instance()    
//__________________________________________________________________________________________________
Int_t AliRICHParam::PadNeighbours(Int_t iPadX,Int_t iPadY,Int_t listX[4],Int_t listY[4])
{
//Determines all the neighbouring pads for the given one (iPadX,iPadY). Returns total number of these pads.
//Dead zones are taken into account, meaning pads from different sector are not taken. 
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
Int_t AliRICHParam::Loc2Sec(const TVector2 &v2)
{
// Determines sector containing the given point. y ^  5 6
//                                                 |  3 4
//                                                 |  1 2
//                                                  -------> x 
// Arguments: v2- LORS position [cm]
//   Returns: sector code
  Double_t x0=0; Double_t x1=SecSizeX(); Double_t x2=SecSizeX()+DeadZone(); Double_t x3=PcSizeX();
  Double_t y0=0; Double_t y1=SecSizeY(); Double_t y2=SecSizeY()+DeadZone(); Double_t y3=2*SecSizeY()+DeadZone(); 
  Double_t y4=PcSizeY()-SecSizeY();      Double_t y5=PcSizeY();
  
  Int_t sector=-1;  
  if     (v2.X() >= x0 && v2.X() <= x1 )  sector=1;
  else if(v2.X() >= x2 && v2.X() <= x3 )  sector=2;
  else                                    return -1;
  
  if     (v2.Y() >= y0 && v2.Y() <= y1 )  ;                    //sectors 1 or 2 
  else if(v2.Y() >= y2 && v2.Y() <= y3 )  sector+=2;           //sectors 3 or 4
  else if(v2.Y() >= y4 && v2.Y() <= y5 )  sector+=4;           //sectors 5 or 6
  else                                    return -1;
  return sector;
}//Loc2Sec(Double_t x, Double_t y)
//__________________________________________________________________________________________________
TVector AliRICHParam::Loc2Pad(const TVector2 &loc)
{
//Determines pad number TVector(padx,pady) containing the given point x2 defined in the chamber RS.
//Pad count starts in lower left corner from 1,1 to 144,160 in upper right corner of a chamber.
//y ^  5 6
//  |  3 4
//  |  1 2
//   -------> x  
  TVector pad(2);
  Int_t sec=Loc2Sec(loc);//trasforms x2 to sector reference system
  if(sec==-1) {pad[0]=pad[1]=-1; return pad;}
//first we deal with x  
  if(sec==1||sec==3||sec==5)    pad[0]=           Int_t(            loc.X()   / PadSizeX() )+1; //sector 1 or 3 or 5
  else                          pad[0]=NpadsX() - Int_t( (PcSizeX()-loc.X())  / PadSizeX() )  ; //sector 2 or 4 or 6
//second deal with y
       if(sec==1||sec==2)       pad[1]=Int_t(             loc.Y()                / PadSizeY())+1;               //sector 1 or 2 
  else if(sec==3||sec==4)       pad[1]=Int_t( (loc.Y()-SecSizeY()-DeadZone()) / PadSizeY())+NpadsYsec()+1;  //sector 3 or 4    
  else                          pad[1]=NpadsY() - Int_t( (PcSizeY()-loc.Y())     / PadSizeY());                //sector 5 or 6        
  return pad;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHParam::Lors2Pad(Double_t x,Double_t y)
{
// Determines abs pad number containing the given point (x,y) defined in the chamber RS.
// Pad count starts in lower left corner from 1,1 to 144,160 in upper right corner of a chamber.
// y ^  5 6
//   |  3 4
//   |  1 2
//    -------> x
  Int_t padx,pady;
  if     (x>=          0           && x<=  SecSizeX()            )  padx=  1              + Int_t(            x /PadSizeX() ); //sector 1 or 3 or 5
  else if(x>=SecSizeX()+DeadZone() && x<=  PcSizeX()             )  padx=     NpadsX()    - Int_t( (PcSizeX()-x)/PadSizeX() ); //sector 2 or 4 or 6
  else                                                              return -1;                                             //dead zone or out of chamber


  if     (y>=          0           && y<=  SecSizeY()            )  pady= 1 +               Int_t(            y /PadSizeY()              ); //sector 1 or 2
  else if(y>=SecSizeY()+DeadZone() && y<=2*SecSizeY()+DeadZone() )  pady= 1 + NpadsYsec() + Int_t( (y-SecSizeY()-DeadZone()) / PadSizeY()); //sector 3 or 4
  else if(y>= PcSizeY()-SecSizeY() && y<=  PcSizeY()             )  pady=     NpadsY()    - Int_t( (PcSizeY()-y)/PadSizeY()              ); //sector 5 or 6
  else                                                              return -1;                                               //dead zone or out of chamber

  return AliRICHDigit::P2A(0,padx,pady);
}//Lors2Pad()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHParam::Pad2Sec(const TVector &pad)
{
//Determines sector containing the given pad.
  Int_t sector=-1;      
  if     (pad[0] >= 1           && pad[0] <=   NpadsXsec() )    {sector=1;}
  else if(pad[0] >  NpadsXsec() && pad[0] <=   NpadsX()    )    {sector=2;} 
  else                                                         AliDebugClass(1,Form("Wrong pad (%3.0f,%3.0f)",pad[0],pad[1]));
    
  if     (pad[1] >= 1             && pad[1] <=   NpadsYsec() )    {}
  else if(pad[1] >  NpadsYsec()   && pad[1] <= 2*NpadsYsec() )    {sector+=2;}
  else if(pad[1] >  2*NpadsYsec() && pad[1] <=   NpadsY()    )    {sector+=4;}
  else                                                         AliDebugClass(1,Form("Wrong pad (%3.0f,%3.0f)",pad[0],pad[1]));

  return sector;
}//Pad2Sec()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector2 AliRICHParam::Pad2Loc(TVector pad)
{
//Returns position of the center of the given pad in local system of the chamber (cm)    
// y ^  5 6
//   |  3 4        sector numbers
//   |  1 2
//    -------> x  
  Double_t x=-1,y=-1;
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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector2 AliRICHParam::Pad2Loc(Int_t pad)
{
// Converts absolute pad number to local position in LORS
// LORS is a chamber  reference system with origin in left-down coner looking from IP
// Arguments: pad- absolute pad number
//   Returns: pad center position as TVector2 in PCRS  
  TVector2 pos;
  pos.Set((Pad2PadX(pad)-0.5)*PadSizeX() , (Pad2PadY(pad)-0.5)*PadSizeY());//set to sector LORS
  return pos;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliRICHParam::GainSag(Double_t x,Int_t sector)
{
//Returns % of gain variation due to wire sagita.
//All curves are parametrized as per sector basis, so x must be apriory transformed to the Sector RS.    
//Here x is a distance along wires.  
  x-=SecSizeX()/2;
  if(x>SecSizeX()) x-=SecSizeX(); 
  switch(HV(sector)){
    case 2150: return 9e-6*TMath::Power(x,4)+2e-7*TMath::Power(x,3)-0.0316*TMath::Power(x,2)-3e-4*x+25.367;//%
    case 2100: return 8e-6*TMath::Power(x,4)+2e-7*TMath::Power(x,3)-0.0283*TMath::Power(x,2)-2e-4*x+23.015;
    case 2050: return 7e-6*TMath::Power(x,4)+1e-7*TMath::Power(x,3)-0.0254*TMath::Power(x,2)-2e-4*x+20.888;
    case 2000: return 6e-6*TMath::Power(x,4)+8e-8*TMath::Power(x,3)-0.0227*TMath::Power(x,2)-1e-4*x+18.961;
    default:   return 0;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliRICHParam::QdcSag(Int_t iPad)
{
// It was observed at BNL that wires are affected by gravitation field providing a significant sagita leading to the local electric field variation
// which means that different pads produce different signals.
// Arguments:  iPad- absolute pad number
//   Returns:  gain variation due to wire sagita 0 < QdcSag < 1.
// Curves are parametrised in terms of distance x (cm) along wires having 0 on the left edge of the photocathode
  Double_t x=AliRICHDigit::P2X(iPad)*PadSizeX()-0.5*PadSizeX(); //center of the padx (count from 1)
  switch(HV(iPad)){
    case 2150: return 0.01*(9e-6*TMath::Power(x,4)+2e-7*TMath::Power(x,3)-0.0316*TMath::Power(x,2)-3e-4*x+25.367);//function is a fit in % so multiply by 0.01
    case 2100: return 0.01*(8e-6*TMath::Power(x,4)+2e-7*TMath::Power(x,3)-0.0283*TMath::Power(x,2)-2e-4*x+23.015);
    case 2050: return 0.01*(7e-6*TMath::Power(x,4)+1e-7*TMath::Power(x,3)-0.0254*TMath::Power(x,2)-2e-4*x+20.888);
    case 2000: return 0.01*(6e-6*TMath::Power(x,4)+8e-8*TMath::Power(x,3)-0.0227*TMath::Power(x,2)-1e-4*x+18.961);
    default:   return 0;
  }
}//QdcSag()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHParam::QdcTot(Int_t iPad,Double_t e)
{
// Calculates the total charge produced by the hit. Method: 
// 1. number of electrons is calculated as energy lost in amp gas divided by ionisation potential (for photon only one electron as Etot is always less then ionization potential)
// 2. each electron imposes a charge distributed as Poisson with QdcEle() mean. Different pads produce different means. See QdcEle().
// Arguments:  iPad- absolute pad number contaning the hit;
//                e- Eloss for mip in amplification gas or Etot for photon
//   Returns:  charge parametrised in QDC channels.
  Int_t iNele=Int_t(e/IonPot()); if(iNele==0) iNele=1;//e < ion. pot. means it's photoelectron
  Double_t dQdc=0;
  for(Int_t i=1;i<=iNele;i++) dQdc+=-QdcEle(iPad)*TMath::Log(gRandom->Rndm());
  return Int_t(dQdc);
}//QdcTot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHParam::TotQdc(TVector2 x2,Double_t eloss)
{
//Calculates the total charge produced by the eloss in point x2 (Chamber RS).
//Returns this change parametrised in QDC channels, or 0 if the hit in the dead zone.
//eloss=0 means photon which produces 1 electron only eloss > 0 for Mip
  if(Loc2Sec(x2)==-1) return 0; //hit in the dead zone     
  Int_t iNelectrons=Int_t(eloss/IonisationPotential()); if(iNelectrons==0) iNelectrons=1;
  Double_t qdc=0;
  for(Int_t i=1;i<=iNelectrons;i++) qdc+=-Gain(x2)*TMath::Log(gRandom->Rndm());
  return Int_t(qdc);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliRICHParam::FracQdc(const TVector2 &x2,const TVector &pad)
{
//Calculates the charge fraction induced to given pad by the hit from the given point.
//Integrated Mathieson distribution is used.  
  TVector2 center2=Pad2Loc(pad);//gives center of requested pad
  Double_t normXmin=(x2.X()-center2.X()-PadSizeX()/2)  /Pc2Cath();//parametrise for Mathienson
  Double_t normXmax=(x2.X()-center2.X()+PadSizeX()/2)  /Pc2Cath();
  Double_t normYmin=(x2.Y()-center2.Y()-PadSizeY()/2)  /Pc2Cath();
  Double_t normYmax=(x2.Y()-center2.Y()+PadSizeY()/2)  /Pc2Cath();
 
//requested pad might not belong to the sector of the given hit position, hence the check:
  return (Loc2Sec(x2)!=Pad2Sec(pad)) ? 0:Mathieson(normXmin, normYmin, normXmax, normYmax);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliRICHParam::Mathieson(Double_t x1,Double_t y1,Double_t x2,Double_t y2)
{
// This is the answer to electrostatic problem of charge distrubution in MWPC described elsewhere. (NIM A370(1988)602-603)
// Arguments: x1- diff between center of distribution which is a hit position and left edge of interested pad divided by anod-cathode distance
//            x2- right edge of the pad
//            y1- up edge of the pad
//            y2- bottom edge of the pad
//  Returns: a charge fraction [0-1] imposed into the pad
  const Double_t kSqrtKx3=0.77459667;const Double_t kX2=0.962;const Double_t kX4=0.379;
  const Double_t kSqrtKy3=0.77459667;const Double_t kY2=0.962;const Double_t kY4=0.379;

  Double_t ux1=kSqrtKx3*TMath::TanH(kX2*x1);
  Double_t ux2=kSqrtKx3*TMath::TanH(kX2*x2);
  Double_t uy1=kSqrtKy3*TMath::TanH(kY2*y1);
  Double_t uy2=kSqrtKy3*TMath::TanH(kY2*y2);
  return 4*kX4*(TMath::ATan(ux2)-TMath::ATan(ux1))*kY4*(TMath::ATan(uy2)-TMath::ATan(uy1));
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector AliRICHParam::Loc2Area(const TVector2 &x2)
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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliRICHParam::IsOverTh(Int_t ,TVector ,Double_t qdc)
{
// Checks if the current QDC is over threshold and FEE will save this value to data concentrator.
// This is done on pad by pad level, so the pad pedestal map is to be used. ??????????????  
// Arguments: 
//   Returns:  true if QDC over treshold 
  return (qdc>NsigmaTh()*(SigmaThMean()+(1.-2*gRandom->Rndm())*SigmaThSpread())); //??????????? to be change to real values
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGeoMatrix* AliRICHParam::Matrix(Int_t iChamN,Int_t iPlane)
{
  TGeoHMatrix *pMatrix=new TGeoHMatrix;
        
  const Double_t kAngHor=19.5; //  horizontal angle between chambers  19.5 grad
  const Double_t kAngVer=20;   //  vertical angle between chambers    20   grad     
  const Double_t kAngCom=30;   //  common RICH rotation with respect to x axis  30   grad     
  
  pMatrix->RotateY(90);        //rotate around y since initial position is in XY plane -> now in YZ plane
  Double_t trans[3]={490,0,0}; //center of the chamber is on window-gap surface
    
  switch(iPlane){
    case kCenter:                          break; 
    case kPc    :   trans[0]+=PcZ();       break; 
    case kRad   :   trans[0]+=RadZ();      break; 
    case kAnod  :   trans[0]+=AnodZ();     break;      
    default:               return 0;       break; 
  }
  pMatrix->SetTranslation(trans); //now plane in YZ is shifted along x 
  
  switch(iChamN){
    case 1:                pMatrix->RotateY(kAngHor);  pMatrix->RotateZ(-kAngVer);  break; //right and down 
    case 2:                                            pMatrix->RotateZ(-kAngVer);  break; //down              
    case 3:                pMatrix->RotateY(kAngHor);                               break; //right 
    case 4:                                                                         break; //no rotation
    case 5:                pMatrix->RotateY(-kAngHor);                              break; //left   
    case 6:                                            pMatrix->RotateZ(kAngVer);   break; //up
    case 7:                pMatrix->RotateY(-kAngHor); pMatrix->RotateZ(kAngVer);   break; //left and up 
    default:               return 0;                                  break;
  }//switch(iChamber)
  pMatrix->RotateZ(kAngCom);     //apply common rotation  in XY plane
  return pMatrix;
}//Matrix()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector3  AliRICHParam::Lors2Mars(Int_t iChId,Double_t x,Double_t y,Int_t iPlnId)
{
// Trasform from LORS to MARS
// Arguments: iChId  - chamber code 1..7
//            x,y    - point in LORS 
//            iPlnN  - chamber plane code  might be kPc kRad kCenter kAnod    
  Double_t z=0;
  switch(iPlnId){
    case kPc    : z=PcZ()  ; break;
    case kAnod  : z=AnodZ(); break;
    case kCenter: z=0      ; break;
    case kRad   : z=RadZ() ; break;
  }
  Double_t lors[3]={x-0.5*PcSizeX(),y-0.5*PcSizeY(),z},  mars[3]; 
  fMatrix[iChId-1]->LocalToMaster(lors,mars);
  return TVector3(mars);
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector3  AliRICHParam::Lors2MarsVec(Int_t iCh,const TVector3 &p)
{
  Double_t mars[3], lors[3]; p.GetXYZ(lors);  
  fMatrix[iCh-1]->LocalToMasterVect(lors,mars);
  return TVector3(mars);
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector2  AliRICHParam::Mars2Lors(Int_t iChId,const TVector3 &x,Int_t iPlnId)
{
// Trasform from MARS to LORS
// Arguments: iChId  - chamber code 1..7
//            mars   - point in MARS 
//            iPlnN  - chamber plane code  might be kPc kRad kCenter kAnod    
  Double_t z=0;
  switch(iPlnId){
    case kPc    : z=PcZ()  ; break;
    case kAnod  : z=AnodZ(); break;
    case kCenter: z=0      ; break;
    case kRad   : z=RadZ() ; break;
  }
  Double_t lors[3],mars[3]; 
  x.GetXYZ(mars);
  fMatrix[iChId-1]->MasterToLocal(mars,lors);
  return TVector2(lors[0]+0.5*PcSizeX(),lors[1]+0.5*PcSizeY());
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector3  AliRICHParam::Mars2LorsVec(Int_t iCh,const TVector3 &p)
{
  Double_t mars[3], lors[3]; p.GetXYZ(mars);  
  fMatrix[iCh-1]->MasterToLocalVect(mars,lors);
  return TVector3(lors);
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector3  AliRICHParam::Lors2MarsOld(Int_t iChId,Double_t x,Double_t y,Int_t iPlnId)
{
// Trasform from LORS to MARS
// Arguments: iChId  - chamber code 0..6
//            x,y    - point in LORS 
//            iPlnN  - chamber plane code  might be kPc kRad kCenter kAnod    
  TGeoMatrix *pMatrix=Matrix(iChId,iPlnId);
  Double_t lors[3]={x-0.5*PcSizeX(),y-0.5*PcSizeY(),0},  mars[3];  pMatrix->LocalToMaster(lors,mars);  delete pMatrix;
  return TVector3(mars);
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector2  AliRICHParam::Mars2LorsOld(Int_t iChamN,const TVector3 &x,Int_t iPlaneN)
{
  TGeoMatrix *pMatrix=Matrix(iChamN,iPlaneN);
  Double_t mars[3]={x.X(),x.Y(),x.Z()}  , lors[3];  pMatrix->MasterToLocal(mars,lors);  delete pMatrix;
  return TVector2(lors[0]+0.5*PcSizeX(),lors[1]+0.5*PcSizeY());
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector3  AliRICHParam::Center(Int_t iChamN,Int_t iPlaneN)
{
  TGeoMatrix *pMatrix=Matrix(iChamN,iPlaneN);
  Double_t mars[3]  , lors[3]={0,0,0};  pMatrix->LocalToMaster(lors,mars);  delete pMatrix;
  return TVector3(mars);
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TVector3  AliRICHParam::Norm(Int_t iChamN)
{
  TGeoMatrix *pMatrix=Matrix(iChamN,kPc);
  Double_t mars[3] , lors[3]={0,0,1};  pMatrix->LocalToMasterVect(lors,mars);  delete pMatrix;
  return TVector3(mars);
}    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHParam::Hit2SDigs(Int_t iHitPad,Double_t e,TClonesArray *pSDigLst)
{
// Determines a number of pads affected by the hit and calculates the charge induced to each pad.
// Integrated Mathieson distribution is used. Invoked from AliRICHvX::Hits2SDigits()
// Arguments: iHitPad  - hit pad absolute number
//            e        - energy (GeV) of this hit (Eloss for mip or Etot for photon)
//            pSDigLst - pointer to clones array to store in calculated sdigits
//   Returns:          total QDC for this hit
  Int_t iQtot=QdcTot(iHitPad,e);                    //total QDC value collected for this hit
  Int_t a=1;                                        //analise current pad +- a pads in both directions
  Int_t iLeftX=0,iBotY=0,iRightX=0,iTopY=0;         //area of disintegration for cluster formation, shifts to hit pad, not pad numbers
  if(AliRICHDigit::P2X(iHitPad) >                         a) iLeftX =-a;//determine area of disintegration as hit pad +- parametrised number
  if(AliRICHDigit::P2X(iHitPad) < AliRICHDigit::kPadsSecX-a) iRightX= a;//of pads. this number is determined by5 sigmas of Mathieson shape
  if(AliRICHDigit::P2Y(iHitPad) >                         a) iBotY  =-a;//see RICH TDR page 29
  if(AliRICHDigit::P2Y(iHitPad) < AliRICHDigit::kPadsSecY-a) iTopY  = a;//also boundary conditions are checked (edge of sector aka PC)
  Int_t iPadsCnt=0;
  for(Int_t iShiftX=iLeftX;iShiftX<=iRightX;iShiftX++){//affected pads loop iShiftX is a distance (in pads) between hit pad and pad under analisys
    for(Int_t iShiftY=iBotY;iShiftY<=iTopY;iShiftY++){//affected pads loop
      iHitPad+=AliRICHDigit::kPadAbsX*iShiftX+iShiftY;
      Double_t x1=PadSizeX()/Pc2Cath()*(iShiftX-0.5);//parametrise for Mathienson
      Double_t x2=PadSizeX()/Pc2Cath()*(iShiftX+0.5);//parametrise for Mathienson
      Double_t y1=PadSizeY()/Pc2Cath()*(iShiftY-0.5);//parametrise for Mathienson
      Double_t y2=PadSizeY()/Pc2Cath()*(iShiftY+0.5);//parametrise for Mathienson
      (*pSDigLst)[iPadsCnt++]= new AliRICHDigit(iHitPad,iQtot*Mathieson(x1,y1,x2,y2));
    }//Y loop
  }//X loop
  return iQtot;
}//Hit2SDigs() for abs pad 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHParam::Hit2SDigs(TVector2 hitX2,Double_t e,TClonesArray *pSDigLst)
{
// Determines a number of pads affected by the hit and calculates the charge induced to each pad.
// Integrated Mathieson distribution is used. Invoked from AliRICHvX::Hits2SDigits()
// Arguments: hitX2    - hit position in LORS, cm
//            e        - energy (GeV) of this hit (Eloss for mip or Etot for photon)
//            pSDigLst - pointer to clones array to store in calculated sdigits
//   Returns:          total QDC for this hit
  Int_t iQtot=TotQdc(hitX2,e);//total charge produced by hit, 0 if hit in dead zone
  if(iQtot==0) return 0;

  TVector  hitPad=Loc2Pad(hitX2);  TVector2 padCenterX2=Pad2Loc(hitPad);     //shift the hit position to the nearest anod wire
  TVector2 anod;
  if((hitX2.Y()-padCenterX2.Y())>0) anod.Set(hitX2.X(),padCenterX2.Y()+AnodPitch()/2); //upper part of the pad: shift to upper anod wire
  else                              anod.Set(hitX2.X(),padCenterX2.Y()-AnodPitch()/2); //lower part of the pad: shift to lower anod wire

  TVector area=Loc2Area(anod);//determine affected pads, dead zones analysed inside
  TVector pad(2); //current pad
  Int_t iPadsCnt=0;
  for(pad[1]=area[1];pad[1]<=area[3];pad[1]++){//affected pads loop
    for(pad[0]=area[0];pad[0]<=area[2];pad[0]++){
      Double_t dQpad=iQtot*FracQdc(anod,pad);
      if(dQpad>0.1) (*pSDigLst)[iPadsCnt++]= new AliRICHDigit(pad,dQpad);//make sdigit if Qpad is large enough, meaning after merging there is a chance to go above threshold
    }//X loop
  }//Y loop
  return iQtot;
}//Hit2SDigs() for TVector2
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliRICHParam::SigmaSinglePhotonFormula(Double_t thetaCer, Double_t phiCer, Double_t theta, Double_t phi, Double_t beta)
{
  TVector3 v(-999,-999,-999);

  v.SetX(ErrLoc(thetaCer,phiCer,theta,phi,beta));
  v.SetY(ErrGeom(thetaCer,phiCer,theta,phi,beta));
  v.SetZ(ErrCrom(thetaCer,phiCer,theta,phi,beta));

  return v.Mag2();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliRICHParam::ErrLoc(Double_t thetaC, Double_t phiC, Double_t Ptheta, Double_t Pphi, Double_t beta)
{
//par->RefIdxC6F14(par->MeanCkovEnergy())
//(Float_t)1.29337525367736816e+00
Double_t RefC6F14m = 1.29337;
  Double_t Hgap = Pc2Win();
  Double_t dphi = phiC - Pphi;

  Double_t alpha =TMath::Cos(Ptheta)-TMath::Tan(thetaC)*TMath::Cos(dphi)*TMath::Sin(Ptheta);
  Double_t k = 1.-RefC6F14m*RefC6F14m+alpha*alpha/(beta*beta);
  if (k<0) k=0; //PH more investigation needed...

  Double_t mu = TMath::Sin(Ptheta)*TMath::Sin(Pphi) + TMath::Tan(thetaC)*(TMath::Cos(Ptheta)*TMath::Cos(dphi)*TMath::Sin(Pphi)
+ TMath::Sin(dphi)*TMath::Cos(Pphi));

  Double_t e = TMath::Sin(Ptheta)*TMath::Cos(Pphi)+TMath::Tan(thetaC)*(TMath::Cos(Ptheta)*TMath::Cos(dphi)*TMath::Cos(Pphi) -TMath::Sin(dphi)*TMath::Sin(Pphi));

  Double_t kk = beta*TMath::Sqrt(k)/(Hgap*alpha);
  Double_t dtdxc = kk*(k*(TMath::Cos(dphi)*TMath::Cos(Pphi) - TMath::Cos(Ptheta)*TMath::Sin(dphi)*TMath::Sin(Pphi)) - ( alpha*
  mu/(beta*beta) )*TMath::Sin(Ptheta)*TMath::Sin(dphi));

  Double_t dtdyc = kk*(k*(TMath::Cos(dphi)*TMath::Sin(Pphi) + TMath::Cos(Ptheta)*TMath::Sin(dphi)*TMath::Cos(Pphi)) + ( alpha*
  e/(beta*beta) )* TMath::Sin(Ptheta)*TMath::Sin(dphi));

  return  TMath::Sqrt(0.2*0.2*dtdxc*dtdxc + 0.25*0.25*dtdyc*dtdyc);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliRICHParam::ErrCrom(Double_t thetaC, Double_t phiC, Double_t Ptheta, Double_t Pphi, Double_t beta)
{
  Double_t dphi = phiC - Pphi;
  Double_t RefC6F14m = 1.29337;
  Double_t alpha =TMath::Cos(Ptheta)-TMath::Tan(thetaC)*TMath::Cos(dphi)*TMath::Sin(Ptheta);

  Double_t dtdn = TMath::Cos(Ptheta)*RefC6F14m*beta*beta/(alpha*TMath::Tan(thetaC));
            
  Double_t f = 0.00928*(7.75-5.635)/TMath::Sqrt(12.);

  return f*dtdn;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliRICHParam::ErrGeom(Double_t thetaC, Double_t phiC, Double_t Ptheta, Double_t Pphi, Double_t beta )
{

  Double_t Tr = RadThick();
  Double_t Xep = 0.5*Tr;

  Double_t dphi = phiC - Pphi;
  Double_t RefC6F14m = 1.29337;
  Double_t alpha =TMath::Cos(Ptheta)-TMath::Tan(thetaC)*TMath::Cos(dphi)*TMath::Sin(Ptheta);

  Double_t k = 1.-RefC6F14m*RefC6F14m+alpha*alpha/(beta*beta);

  Double_t Hgap = Pc2Win();


  Double_t eTr = (Tr - Xep)*beta*TMath::Sqrt(k)/(Hgap*alpha);
  Double_t lambda = 1.-TMath::Sin(Ptheta)*TMath::Sin(Ptheta)*TMath::Sin(phiC)*TMath::Sin(phiC);

  Double_t c = 1./(1.+ eTr*k/(alpha*alpha*TMath::Cos(thetaC)*TMath::Cos(thetaC)));
  Double_t I = beta*TMath::Tan(thetaC)*lambda*TMath::Power(k,1.5);
  Double_t II = 1.+eTr*beta*I;

  Double_t err = c * (I/(alpha*alpha*Hgap) +  II* (1.-lambda) / ( alpha*alpha*Hgap*beta*(1.+eTr)) );
  Double_t TrErr = Tr/(TMath::Sqrt(12.)*TMath::Cos(Ptheta));

  return TrErr*err;
}//ErrGeom()

#endif //AliRICHParam_h
