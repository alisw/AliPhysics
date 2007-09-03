#ifndef AliHMPIDParam_h
#define AliHMPIDParam_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TMath.h>
#include <TNamed.h>        //base class
#include <TGeoManager.h>   //Instance()
#include <TVector3.h>      //Lors2Mars() Mars2Lors()
 
// Class providing all the needed parametrised information
// to construct the geometry, to define segmentation and to provide response model
// In future will also provide all the staff needed for alignment and calibration

class AliHMPIDParam :public TNamed  
{
public:
//ctor&dtor    
  virtual        ~AliHMPIDParam()                                    {for(Int_t i=0;i<7;i++) delete fM[i]; delete fgInstance; fgInstance=0;}
         void     Print(Option_t *opt="") const;                                         //print current parametrization
  static inline AliHMPIDParam* Instance();                                //pointer to AliHMPIDParam singleton
  static inline AliHMPIDParam* InstanceNoGeo();                           //pointer to AliHMPIDParam singleton without geometry.root for MOOD, displays, ...
//geo info
  enum EChamberData{kMinCh=0,kMaxCh=6,kMinPc=0,kMaxPc=5};      //Segmenation
  enum EPadxData{kPadPcX=80,kMinPx=0,kMaxPx=79,kMaxPcx=159};   //Segmentation structure along x
  enum EPadyData{kPadPcY=48,kMinPy=0,kMaxPy=47,kMaxPcy=143};   //Segmentation structure along y 

  static Float_t SizePadX    (                               )     {return fgCellX; /*return 0.804;*/}                    //pad size x, [cm]  
  static Float_t SizePadY    (                               )     {return fgCellY; /*0.84*/}                           //pad size y, [cm]  

  static Float_t SizePcX    (                                )     {return fgPcX;}                                    // PC size x
  static Float_t SizePcY    (                                )     {return fgPcY;}                                    // PC size y
  static Float_t MaxPcX      (Int_t iPc                      )     {return fgkMaxPcX[iPc];}                           // PC limits
  static Float_t MaxPcY      (Int_t iPc                      )     {return fgkMaxPcY[iPc];}                           // PC limits
  static Float_t MinPcX      (Int_t iPc                      )     {return fgkMinPcX[iPc];}                           // PC limits
  static Float_t MinPcY      (Int_t iPc                      )     {return fgkMinPcY[iPc];}                           // PC limits
  static Int_t   Nsig        (                               )     {return fgSigmas;}                                 //Getter n. sigmas for noise
  static Float_t SizeAllX    (                               )     {return fgAllX/*fgkMaxPcX[5]*/;}                             //all PCs size x, [cm]        
  static Float_t SizeAllY    (                               )     {return fgAllY/*fgkMaxPcY[5]*/;}                             //all PCs size y, [cm]    

  static Float_t LorsX       (Int_t pc,Int_t padx             )     {return (padx    +0.5)*SizePadX()+fgkMinPcX[pc]; } //center of the pad x, [cm]

  static Float_t LorsY       (Int_t pc,Int_t pady            )     {return (pady    +0.5)*SizePadY()+fgkMinPcY[pc];   } //center of the pad y, [cm]

  inline static void   Lors2Pad(Float_t x,Float_t y,Int_t &pc,Int_t &px,Int_t &py);                                       //(x,y)->(pc,px,py) 

  static Int_t   Abs         (Int_t ch,Int_t pc,Int_t x,Int_t y)   {return ch*100000000+pc*1000000+x*1000+y;         } //(ch,pc,padx,pady)-> abs pad
  static Int_t   A2C         (Int_t pad                      )     {return pad/100000000;                             } //abs pad -> chamber
  static Int_t   A2P         (Int_t pad                      )     {return pad%100000000/1000000;                    } //abs pad -> pc 
  static Int_t   A2X         (Int_t pad                      )     {return pad%1000000/1000;                          } //abs pad -> pad X 
  static Int_t   A2Y         (Int_t pad                      )     {return pad%1000;                                  } //abs pad -> pad Y 

  static Bool_t  IsOverTh    (Float_t q                      )     {return q >= fgSigmas;                            } //is digit over threshold?

  inline static Bool_t IsInDead(Float_t x,Float_t y        );                                                          //is point in dead area?
  static Bool_t  IsInside    (Float_t x,Float_t y,Float_t d=0)     {return  x>-d&&y>-d&&x<fgkMaxPcX[kMaxPc]+d&&y<fgkMaxPcY[kMaxPc]+d; } //is point inside chamber boundary?


            Double_t   MeanIdxRad              ()const {return 1.29204;}   //<--TEMPORAR--> to be removed in future  Mean ref index C6F14
            Double_t   MeanIdxWin              ()const {return 1.57819;}   //<--TEMPORAR--> to be removed in future. Mean ref index quartz
            Float_t    DistCut                 ()const {return 1.0;}       //<--TEMPORAR--> to be removed in future. Cut for MIP-TRACK residual 
            Float_t    QCut                    ()const {return 100;}       //<--TEMPORAR--> to be removed in future. Separation PHOTON-MIP charge 
            Float_t    MultCut                 ()const {return 200;}       //<--TEMPORAR--> to be removed in future. Multiplicity cut to activate WEIGHT procedure 



  static        Int_t      Stack(Int_t evt=-1,Int_t tid=-1);              //Print stack info for event and tid
  static        Int_t      StackCount(Int_t pid,Int_t evt);               //Counts stack particles of given sort in given event  
  static        void       IdealPosition(Int_t iCh,TGeoHMatrix *m);       //ideal position of given chamber 
  //trasformation methodes
  void     Lors2Mars   (Int_t c,Float_t x,Float_t y,Double_t *m,Int_t pl=kPc)const{Double_t z=0; switch(pl){case kPc:z=8.0;break; case kAnod:z=7.806;break; case kRad:z=-1.25; break;}   Double_t l[3]={x-fX,y-fY,z};  fM[c]->LocalToMaster(l,m); }    
  TVector3 Lors2Mars   (Int_t c,Float_t x,Float_t y,            Int_t pl=kPc)const{Double_t m[3];Lors2Mars(c,x,y,m,pl); return TVector3(m);    }//MRS->LRS  
  void     Mars2Lors   (Int_t c,Double_t *m,Float_t &x ,Float_t &y          )const{Double_t l[3];fM[c]->MasterToLocal(m,l);x=l[0]+fX;y=l[1]+fY;}//MRS->LRS
  void     Mars2LorsVec(Int_t c,Double_t *m,Float_t &th,Float_t &ph         )const{Double_t l[3]; fM[c]->MasterToLocalVect(m,l); 
                                                                                   Float_t pt=TMath::Sqrt(l[0]*l[0]+l[1]*l[1]); 
                                                                                           th=TMath::ATan(pt/l[2]); 
                                                                                           ph=TMath::ATan2(l[1],l[0]);}    
  TVector3 Norm        (Int_t c                                             )const{Double_t n[3]; Norm(c,n); return TVector3(n);               }//norm 
  void     Norm        (Int_t c,Double_t *n                                 )const{Double_t l[3]={0,0,1};fM[c]->LocalToMasterVect(l,n);        }//norm
  void     Point       (Int_t c,Double_t *p,Int_t plane                     )const{Lors2Mars(c,0,0,p,plane);}      //point of given chamber plane

  enum EPlaneId {kPc,kRad,kAnod};            //3 planes in chamber 

  static Int_t fgSigmas;   //sigma Cut

  
protected:
  static /*const*/ Float_t fgkMinPcX[6];                                                           //limits PC
  static /*const*/ Float_t fgkMinPcY[6];                                                           //limits PC
  static /*const*/ Float_t fgkMaxPcX[6];                                                           //limits PC
  static /*const*/ Float_t fgkMaxPcY[6]; 

  static Float_t fgCellX, fgCellY, fgPcX, fgPcY, fgAllX, fgAllY;
         AliHMPIDParam(Bool_t noGeo);             //default ctor is protected to enforce it to be singleton

  static AliHMPIDParam *fgInstance;   //static pointer  to instance of AliHMPIDParam singleton

  TGeoHMatrix *fM[7];                 //pointers to matrices defining HMPID chambers rotations-translations
  Float_t fX;                         //x shift of LORS with respect to rotated MARS 
  Float_t fY;                         //y shift of LORS with respect to rotated MARS   

  
  
  ClassDef(AliHMPIDParam,0)           //HMPID main parameters class
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDParam* AliHMPIDParam::Instance()
{
// Return pointer to the AliHMPIDParam singleton. 
// Arguments: none
//   Returns: pointer to the instance of AliHMPIDParam or 0 if no geometry       
  if(!fgInstance) new AliHMPIDParam(kFALSE);                                //default setting for reconstruction, if no geometry.root -> AliFatal
  return fgInstance;  
}//Instance()    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDParam* AliHMPIDParam::InstanceNoGeo()
{
// Return pointer to the AliHMPIDParam singleton without the geometry.root. 
// Arguments: none
//   Returns: pointer to the instance of AliHMPIDParam or 0 if no geometry       
  if(!fgInstance) new AliHMPIDParam(kTRUE);                               //to avoid AliFatal, for MOOD and displays, use ideal geometry parameters
  return fgInstance;  
}//Instance()    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDParam::IsInDead(Float_t x,Float_t y)
{
// Check is the current point is outside of sensitive area or in dead zones
// Arguments: x,y -position
//   Returns: 1 if not in sensitive zone           
  for(Int_t iPc=0;iPc<6;iPc++)
    if(x>=fgkMinPcX[iPc] && x<=fgkMaxPcX[iPc] && y>=fgkMinPcY[iPc] && y<=fgkMaxPcY [iPc]) return kFALSE; //in current pc
  
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDParam::Lors2Pad(Float_t x,Float_t y,Int_t &pc,Int_t &px,Int_t &py)
{
// Check the pad of given position
// Arguments: x,y- position [cm] in LORS; pc,px,py- pad where to store the result
//   Returns: none
  pc=px=py=-1;
  if     (x>fgkMinPcX[0] && x<fgkMaxPcX[0]) {pc=0; px=Int_t( x               / SizePadX());}//PC 0 or 2 or 4
  else if(x>fgkMinPcX[1] && x<fgkMaxPcX[1]) {pc=1; px=Int_t((x-fgkMinPcX[1]) / SizePadX());}//PC 1 or 3 or 5
  else return;
  if     (y>fgkMinPcY[0] && y<fgkMaxPcY[0]) {      py=Int_t( y               / SizePadY());}//PC 0 or 1
  else if(y>fgkMinPcY[2] && y<fgkMaxPcY[2]) {pc+=2;py=Int_t((y-fgkMinPcY[2]) / SizePadY());}//PC 2 or 3
  else if(y>fgkMinPcY[4] && y<fgkMaxPcY[4]) {pc+=4;py=Int_t((y-fgkMinPcY[4]) / SizePadY());}//PC 4 or 5
  else return;
}
#endif
