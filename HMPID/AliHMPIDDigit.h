#ifndef AliHMPIDDigit_h
#define AliHMPIDDigit_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class of HMPID to manage digits ---> pads
//.
//.
//.

#include <AliDigit.h>      //base class  
#include "TMath.h"         //Mathieson()
#include <AliBitPacking.h> //Raw()


class TClonesArray;        //Hit2Sdi()
  
class AliHMPIDDigit :public AliDigit //TObject-AliDigit-AliHMPIDDigit
{
public:
  enum EChamberData{kMinCh=0,kMaxCh=6,kMinPc=0,kMaxPc=5};      //Segmenation     
  enum EPadxData{kPadPcX=80,kMinPx=0,kMaxPx=79,kMaxPcx=159};   //Segmentation structure along x
  enum EPadyData{kPadPcY=48,kMinPy=0,kMaxPy=47,kMaxPcy=143};   //Segmentation structure along y
//ctor&dtor    
  AliHMPIDDigit(                          ):AliDigit( ),fPad(Abs(-1,-1,-1,-1)),fQ(-1)  {}                         //default ctor
  AliHMPIDDigit(Int_t pad,Int_t q,Int_t *t):AliDigit(t),fPad(pad             ),fQ(q )  {}                         //digit ctor
  AliHMPIDDigit(const AliHMPIDDigit &d    ):AliDigit(d),fPad(d.fPad),fQ(d.fQ)          {}                         //copy ctor
  virtual ~AliHMPIDDigit()                                                             {}                         //dtor   
//framework part    
         Bool_t  IsSortable  (                               )const{return kTRUE;}                                                     //provision to use TObject::Sort() 
  inline Int_t   Compare     (const TObject *pObj            )const;                                                                   //provision to use TObject::Sort()
         void    Draw        (Option_t *opt=""               );                                                                        //TObject::Draw() overloaded
         void    Print       (Option_t *opt=""               )const;                                                                   //TObject::Print() overloaded
//private part  
  static Int_t   Abs         (Int_t ch,Int_t pc,Int_t x,Int_t y)   {return ch*100000000+pc*1000000+x*1000+y;                         } //(ch,pc,padx,pady)-> abs pad
  static Int_t   A2C         (Int_t pad                      )     {return pad/100000000;                                            } //abs pad -> chamber
  static Int_t   A2P         (Int_t pad                      )     {return pad%100000000/1000000;                                    } //abs pad -> pc
  static Int_t   A2X         (Int_t pad                      )     {return pad%1000000/1000;                                         } //abs pad -> pad X 
  static Int_t   A2Y         (Int_t pad                      )     {return pad%1000;                                                 } //abs pad -> pad Y 
         void    AddTidOffset(Int_t offset                   )     {for (Int_t i=0; i<3; i++) if (fTracks[i]>0) fTracks[i]+=offset;  } //needed for merging
         Int_t   Ch          (                               )const{return A2C(fPad);                                                } //chamber number
  static Bool_t  IsOverTh    (Float_t q                      )     {return q >= fgSigmas;                                            } //is digit over threshold?
  static Bool_t  IsInside    (Float_t x,Float_t y,Float_t margin=0){return x>-margin&&y>-margin&&x<SizeAllX()+margin&&y<SizeAllY()+margin;} //is point inside chamber boundary?
         Float_t LorsX       (                               )const{return LorsX(A2P(fPad),A2X(fPad));                               } //center of the pad x, [cm]
  static Float_t LorsX       (Int_t pc,Int_t padx            )     {return (padx    +0.5)*SizePadX()+(pc  %2)*(SizePcX()+SizeDead());} //center of the pad x, [cm]
         Float_t LorsY       (                               )const{return LorsY(A2P(fPad),A2Y(fPad));                               } //center of the pad y, [cm]
  static Float_t LorsY       (Int_t pc,Int_t pady            )     {return (pady    +0.5)*SizePadY()+(pc  /2)*(SizePcY()+SizeDead());} //center of the pad y, [cm]
  inline Float_t IntMathieson(Float_t x,Float_t y            )const;                                                                   //Mathieson distribution 
         Int_t   PadPcX      (                               )const{return A2X(fPad);}                                                 //pad pc x # 0..79
         Int_t   PadPcY      (                               )const{return A2Y(fPad);}                                                 //pad pc y # 0..47
         Int_t   PadChX      (                               )const{return (Pc()%2)*kPadPcX+PadPcX();}                                 //pad ch x # 0..159
         Int_t   PadChY      (                               )const{return (Pc()/2)*kPadPcY+PadPcY();}                                 //pad ch y # 0..143
         Int_t   Pad         (                               )const{return fPad;}                                                      //absolute id of this pad
         Int_t   Pc          (                               )const{return A2P(fPad);}                                                 //PC position number
         Float_t Q           (                               )const{return fQ;}                                                        //charge, [QDC]
  inline void    Raw         (UInt_t &w32,Int_t &ddl,Int_t &r,Int_t &d,Int_t &a)const;                                                 //digit->(w32,ddl,r,d,a)
  inline void    Raw         (UInt_t  w32,Int_t  ddl         );                                                                        //(w32,ddl)->digit
  inline void    Raw         (Int_t ddl,Int_t r,Int_t d,Int_t a);                                                                      //raw->abs pad number
  inline Bool_t  Set         (Int_t c,Int_t p,Int_t x,Int_t y,Int_t tid=0);                                                            //manual creation 
         void    SetQ        (Float_t q                      )     {fQ=q;}                                                             //manual creation 
         void    SetNsig     (Int_t sigmas                   )     {fgSigmas=sigmas;}                                                  //set n sigmas 
  static void    WriteRaw    (TObjArray *pDigLst             );                                                                        //write as raw stream     
  
  static Float_t CathAnoCath (                               )     {return 0.445;}                                                     //Cathode-Anode-cathode pitch
  static Float_t MaxPcX      (Int_t iPc                      )     {return fgkMaxPcX[iPc];}                                            // PC limits
  static Float_t MaxPcY      (Int_t iPc                      )     {return fgkMaxPcY[iPc];}                                            // PC limits
  static Float_t MinPcX      (Int_t iPc                      )     {return fgkMinPcX[iPc];}                                            // PC limits
  static Float_t MinPcY      (Int_t iPc                      )     {return fgkMinPcY[iPc];}                                            // PC limits
  static Int_t   Nsig        (                               )     {return fgSigmas;}                                                  //Getter n. sigmas for noise
  static Float_t SizeAllX    (                               )     {return fgkMaxPcX[5];}                                              //all PCs size x, [cm]        
  static Float_t SizeAllY    (                               )     {return fgkMaxPcY[5];}                                              //all PCs size y, [cm]    
  static Float_t SizeArea    (                               )     {return SizePcX()*SizePcY()*(kMaxPc-kMinPc+1);}                     //sence area, [cm^2]  
  static Float_t SizeDead    (                               )     {return 2.6;}                                                       //dead zone size x, [cm]         
  static Float_t SizeGap     (                               )     {return 8;  }
  static Float_t SizePadX    (                               )     {return 0.8;}                                                       //pad size x, [cm]  
  static Float_t SizePadY    (                               )     {return 0.84;}                                                      //pad size y, [cm]  
  static Float_t SizePcX     (                               )     {return fgkMaxPcX[0];}                                              //PC size x, [cm]        
  static Float_t SizePcY     (                               )     {return fgkMaxPcY[0];}                                              //PC size y, [cm]    
  static Float_t SizeWin     (                               )     {return 0.5;}                                                       //Quartz window width
  static Float_t SizeRad     (                               )     {return 1.5;}                                                       //Rad width   
  inline static Bool_t IsInDead(Float_t x,Float_t y        );                                                                          //is point in dead area?
  inline static void   Lors2Pad(Float_t x,Float_t y,Int_t &pc,Int_t &px,Int_t &py);                                                    //(x,y)->(pc,px,py) 
protected:                                                                   //AliDigit has fTracks[3]
  static Int_t fgSigmas;                                                                                                               //n. sigma to cut on charge 
  static const Float_t fgkMinPcX[6];                                                                                                   //limits PC
  static const Float_t fgkMinPcY[6];                                                                                                   //limits PC
  static const Float_t fgkMaxPcX[6];                                                                                                   //limits PC
  static const Float_t fgkMaxPcY[6];                                                                                                   //limits PC
  static const Float_t fgk1;                                                                                                           //Mathieson parameters
  static const Float_t fgk2;                                                                                                           //...
  static const Float_t fgkSqrtK3;                                                                                                      //...
  static const Float_t fgk4;                                                                                                           //...
  Int_t    fPad;                                                                                                                       //absolute pad number
  Float_t  fQ;                                                               //QDC value, fractions are permitted for summable procedure  
  ClassDef(AliHMPIDDigit,4)                                                  //HMPID digit class       
};//class AliHMPIDDigit

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigit::Lors2Pad(Float_t x,Float_t y,Int_t &pc,Int_t &px,Int_t &py)
{
// Check the pad of given position
// Arguments: x,y- position [cm] in LORS; pc,px,py- pad where to store the result
//   Returns: none
  pc=px=py=-1;
  if     (x>=          0          && x<=  SizePcX()            ) {pc=0; px=Int_t( x                           / SizePadX());}//PC 0 or 2 or 4
  else if(x>=SizePcX()+SizeDead() && x<=  SizeAllX()           ) {pc=1; px=Int_t((x-  SizePcX()-  SizeDead()) / SizePadX());}//PC 2 or 4 or 6
  else return;
  if     (y>=          0          && y<=  SizePcY()            ) {      py=Int_t( y                           / SizePadY());}//PC 0 or 1
  else if(y>=SizePcY()+SizeDead() && y<=2*SizePcY()+SizeDead() ) {pc+=2;py=Int_t((y-  SizePcY()-  SizeDead()) / SizePadY());}//PC 2 or 3
  else if(y>=SizeAllY()-SizePcY() && y<=  SizeAllY()           ) {pc+=4;py=Int_t((y-2*SizePcY()-2*SizeDead()) / SizePadY());}//PC 4 or 5
  else return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDDigit::Compare(const TObject *pObj) const
{
// Used in Sort() method to compare to objects. Note that abs pad structure is first x then y, hence will be sorted on column basis.
// This feature is used in digitizer to facilitate finding of sdigits for the same pad since they all will come together after sorting.
// Arguments: pObj - pointer to object to compare with
//   Retunrs: -1 if AbsPad less then in pObj, 1 if more and 0 if they are the same      
  if     (fPad==((AliHMPIDDigit*)pObj)->Pad()) return  0;
  else if(fPad >((AliHMPIDDigit*)pObj)->Pad()) return  1;
  else                                         return -1;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDDigit::IsInDead(Float_t x,Float_t y)
{
// Check is the current point is outside of sensitive area or in dead zones
// Arguments: x,y -position
//   Returns: 1 if not in sensitive zone           
  if(x<0 || x>SizeAllX() || y<0 || y>SizeAllY()) return kTRUE; //out of pc 
  
  if(x>SizePcX()  && x<SizePcX()+SizeDead())   return kTRUE; //in dead zone along x  
  
  if(y>SizePcY()                       && y<SizePcY()+SizeDead())      return kTRUE; //in first dead zone along y   
  if(y>SizeAllY()-SizePcY()-SizeDead() && y<SizeAllY()-SizePcY())      return kTRUE; //in second dead zone along y   
  return kFALSE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Float_t AliHMPIDDigit::IntMathieson(Float_t x,Float_t y)const
{
// Integration of Mathieson.
// This is the answer to electrostatic problem of charge distrubution in MWPC described elsewhere. (NIM A370(1988)602-603)
// Arguments: x,y- position of the center of Mathieson distribution
//  Returns: a charge fraction [0-1] imposed into the pad
//  K1    =0.28278796
//  K2    =0.96242952
//  SqrtK3=0.77459667
//  K4    =0.37932926

  Float_t ux1=fgkSqrtK3*TMath::TanH(fgk2*(x-LorsX()+0.5*SizePadX())/CathAnoCath());
  Float_t ux2=fgkSqrtK3*TMath::TanH(fgk2*(x-LorsX()-0.5*SizePadX())/CathAnoCath());
  Float_t uy1=fgkSqrtK3*TMath::TanH(fgk2*(y-LorsY()+0.5*SizePadY())/CathAnoCath());
  Float_t uy2=fgkSqrtK3*TMath::TanH(fgk2*(y-LorsY()-0.5*SizePadY())/CathAnoCath());
  return 4*fgk4*(TMath::ATan(ux2)-TMath::ATan(ux1))*fgk4*(TMath::ATan(uy2)-TMath::ATan(uy1));
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigit::Raw(UInt_t &w32,Int_t &ddl,Int_t &r,Int_t &d,Int_t &a)const
{
// Convert digit structure to raw word format
// Arguments: w32,ddl,r,d,a where to write the results
//   Returns: none
  Int_t y2a[6]={5,3,1,0,2,4};

                                  ddl=2*Ch()+Pc()%2;                     //DDL# 0..13
  Int_t tmp=1+Pc()/2*8+PadPcY()/6;  r=(Pc()%2)? 25-tmp:tmp;              //row r=1..24
                                    d=1+PadPcX()/8;                      //DILOGIC# 1..10
                                    a=y2a[PadPcY()%6]+6*(PadPcX()%8);    //ADDRESS 0..47        
      
  w32=0;    
  AliBitPacking::PackWord((UInt_t)fQ,w32, 0,11);  // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq        Qdc               bits (00..11) counts (0..4095)
  AliBitPacking::PackWord(        a ,w32,12,17);  // 3322 2222 2222 1111 1111 1000 0000 0000        DILOGIC address   bits (12..17) counts (0..47)
  AliBitPacking::PackWord(        d ,w32,18,21);  // 1098 7654 3210 9876 5432 1098 7654 3210        DILOGIC number    bits (18..21) counts (1..10)
  AliBitPacking::PackWord(        r ,w32,22,26);  //                                                Row number        bits (22..26) counts (1..24)  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigit::Raw(UInt_t w32,Int_t ddl)
{
// Converts a given raw data word to a digit
// Arguments: w32 - 32 bits raw data word
//            ddl - DDL idx  0 1 2 3 4 ... 13
//   Returns: none
  Int_t r = AliBitPacking::UnpackWord(w32,22,26); assert(1<=r&&r<=24);   //                                         Row number      (1..24)    
  Int_t d = AliBitPacking::UnpackWord(w32,18,21); assert(1<=d&&d<=10);   // 3322 2222 2222 1111 1111 1000 0000 0000 DILOGIC number  (1..10)
  Int_t a = AliBitPacking::UnpackWord(w32,12,17); assert(0<=a&&a<=47);   // 1098 7654 3210 9876 5432 1098 7654 3210 DILOGIC address (0..47)  
  Int_t q = AliBitPacking::UnpackWord(w32, 0,11); assert(0<=q&&q<=4095); // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq Qdc             (0..4095) 
  Raw(ddl,r,d,a);
  fQ=q;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigit::Raw(Int_t ddl,Int_t r,Int_t d,Int_t a)
{
  assert(0<=ddl&&ddl<=13); assert(1<=r&&r<=24); assert(1<=d&&d<=10);   assert(0<=a&&a<=47);  
  Int_t a2y[6]={3,2,4,1,5,0};//pady for a given address (for single DILOGIC chip)
                                  Int_t ch=ddl/2;
  Int_t tmp=(r-1)/8;              Int_t pc=(ddl%2)? 5-2*tmp:2*tmp; 
                                  Int_t px=(d-1)*8+a/6;
        tmp=(ddl%2)?(24-r):r-1;   Int_t py=6*(tmp%8)+a2y[a%6];
  fPad=Abs(ch,pc,px,py);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDDigit::Set(Int_t ch,Int_t pc,Int_t px,Int_t py,Int_t tid)
{
// Manual creation of digit
// Arguments: ch,pc,px,py,qdc,tid  
//   Returns: kTRUE if wrong digit
  if(px<kMinPx || px>kMaxPx) return kTRUE;
  if(py<kMinPy || py>kMaxPy) return kTRUE;

  fPad=Abs(ch,pc,px,py);fTracks[0]=tid;
  fQ=0;
  return kFALSE;
}
#endif
