#ifndef AliHMPIDDigit_h
#define AliHMPIDDigit_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliDigit.h>      //base class  
#include <TMath.h>         //Mathieson()
#include <TRandom.h>       //IsOverTh()  
#include <AliBitPacking.h> //Raw()

class AliHMPIDHit;          //Hit2Sdi()
class TClonesArray;        //Hit2Sdi()
 
class AliHMPIDDigit :public AliDigit //TObject-AliDigit-AliHMPIDDigit
{
public:
  enum EAbsPad {kChAbs=100000000,kPcAbs=1000000,kPadAbsX=1000,kPadAbsY=1};       //absolute pad number structure
  enum ERawData{kDilX=8,kDilY=6,kNdil=10,kNrow=24,kNddls=14};                    //RAW data structure
  enum EPadData{kPcX=2,kPcY=3,kPad1=0,kPadPcX=80,kPadPcY=48,kPadAllX=kPadPcX*kPcX,kPadAllY=kPadPcY*kPcY,kPcAll=kPcX*kPcY,kPadAll=kPadAllX*kPadAllY};   //Segmentation structure 
  enum EPadShif{kC=0,kU=1,kUR=2,kR=3,kDR=4,kD=5,kDL=6,kL=7,kUL=8};                 
//ctor&dtor    
           AliHMPIDDigit(                                                     ):AliDigit( ),fPad(Abs(-1,-1,-1,-1)),fQ(-1)  {}                //default ctor
           AliHMPIDDigit(Int_t pad,Int_t q,Int_t *t                           ):AliDigit(t),fPad(pad             ),fQ(q )  {}                //digit ctor 
  inline   AliHMPIDDigit(Int_t c,Float_t q,Int_t t,Float_t x,Float_t y,Int_t f=0);                                                               //sdigit ctor 
  virtual ~AliHMPIDDigit()                                {} //dtor
//framework part    
         Bool_t  IsSortable  (                               )const{return kTRUE;}                                                     //provision to use TObject::Sort() 
  inline Int_t   Compare     (const TObject *pObj            )const;                                                                   //provision to use TObject::Sort()
         void    Print       (Option_t *opt=""               )const;                                                                   //TObject::Print() overloaded
//private part  
  static Int_t   Abs         (Int_t c,Int_t s,Int_t x,Int_t y)     {return c*kChAbs+s*kPcAbs+x*kPadAbsX+y*kPadAbsY; }                  //(ch,pc,padx,pady)-> abs pad
  static Int_t   A2C         (Int_t pad                      )     {return pad/kChAbs;                              }                  //abs pad -> chamber
  static Int_t   A2P         (Int_t pad                      )     {return pad%kChAbs/kPcAbs;                       }                  //abs pad -> pc
  static Int_t   A2X         (Int_t pad                      )     {return pad%kPcAbs/kPadAbsX;                     }                  //abs pad -> pad X 
  static Int_t   A2Y         (Int_t pad                      )     {return pad%kPadAbsX;                            }                  //abs pad -> pad Y 
         Int_t   Addr        (                               )const{Int_t mapY2A[kDilY]={5,3,1,0,2,4}; return mapY2A[A2Y(fPad)%kDilY]+kDilY*(A2X(fPad)%kDilX);}//raw a=0..47
         void    AddTidOffset(Int_t offset                   )     {for (Int_t i=0; i<3; i++) if (fTracks[i]>0) fTracks[i]+=offset;};  //needed for merging
         Int_t   Ch          (                               )const{return A2C(fPad);                               }                  //chamber number
         Int_t   Dilogic     (                               )const{return 10-PadX()/kDilX;                         }                  //raw d=1..10
  static void    DrawPc      (Bool_t isFill=kTRUE            );                                                                        //draw PCs
         Int_t   Ddl         (                               )const{return (PadX()<kPadPcX) ? 2*Ch() : 2*Ch()+1;    }                  //DDL number 0..13
  static Float_t Hit2Sdi     (AliHMPIDHit *pHit,TClonesArray *);                                                                        //hit -> 9 sdigits, returns total QDC   
  static Bool_t  IsOverTh    (Float_t q                      )     {return q > 6;                                                  }   //is digit over threshold????
  static Bool_t  IsInside    (Float_t x,Float_t y            )     {return x>0&&y>0&&x<SizeAllX()&&y<SizeAllY();                   }   //is point inside pc boundary?
  inline static Bool_t  IsInDead    (Float_t x,Float_t y            );                                                                 //is point in dead area?
         Float_t LorsX       (                               )const{return (PadX()+0.5)*SizePadX()+(Pc()%2)*(SizePcX()+SizeDead());}   //center of the pad x, [cm]
         Float_t LorsY       (                               )const{return (PadY()+0.5)*SizePadY()+(Pc()/2)*(SizePcY()+SizeDead());}   //center of the pad y, [cm]
  inline Float_t Mathieson   (Float_t x,Float_t y            )const;                                                                   //Mathieson distribution 
         Int_t   PadX        (                               )const{return A2X(fPad);}                                                 //x position of the pad
         Int_t   PadY        (                               )const{return A2Y(fPad);}                                                 //y postion of the pad     
         Int_t   Pad         (                               )const{return fPad;}                                                      //absolute id of this pad
         Float_t Q           (                               )const{return fQ;}                                                        //charge, [QDC]
         Int_t   Pc          (                               )const{return A2P(fPad);}                                                 //PC position number
  static void    PrintSize   (                               );                                                                        //print all segmentation sizes      
  inline Int_t   Raw         (UInt_t &w32                    )const;                                                                   //raw
         Int_t   Row         (                               )const{Int_t r=1+Pc()/2*8+PadY()/kDilY; return (Pc()%2)?kNrow-r+1:r;}     //row r=1..24
         void    Set         (Int_t c,Int_t s,Int_t x,Int_t y)     {fPad=Abs(c,s,x,y);}                                                //set new digit
         void    ReadRaw     (Int_t ddl,Int_t r,Int_t d,Int_t a){Int_t mapA2Y[kDilY]={3,2,4,1,5,0};fPad=Abs(ddl/2,ddl%7,d*kDilX+a/kDilY,r*kDilY+mapA2Y[a%kDilY]);} //from raw
  inline void    ReadRaw     (Int_t ddl,UInt_t w32           );                                                                        //read raw word
  
  static Float_t SizeAllX    (                               )     {return SizePadX()*kPadAllX+SizeDead();}                            //all PCs size x, [cm]        
  static Float_t SizeAllY    (                               )     {return SizePadY()*kPadAllY+2*SizeDead();}                          //all PCs size y, [cm]    
  static Float_t SizeArea    (                               )     {return SizePcX()*SizePcY()*kPcAll;}                                //sence area, [cm^2]  
  static Float_t SizeDead    (                               )     {return 2.6;}                                                       //dead zone size x, [cm]         
  static Float_t SizeGap     (                               )     {return 8;  }
  static Float_t SizePadX    (                               )     {return 0.8;}                                                       //pad size x, [cm]  
  static Float_t SizePadY    (                               )     {return 0.84;}                                                      //pad size y, [cm]  
  static Float_t SizePcX     (                               )     {return SizePadX()*kPadPcX;}                                        //PC size x, [cm]        
  static Float_t SizePcY     (                               )     {return SizePadY()*kPadPcY;}                                        //PC size y, [cm]    
  static Float_t SizeWin     (                               )     {return 0.5;}
  static Float_t SizeRad     (                               )     {return 1.5;}     
  static void    TestSeg     (                               );                                                                        //test segmentation
         void    Zoom        (                               ); 
protected:                  //AliDigit has fTracks[3]
  Int_t   fPad;             //absolute pad number is chamber*kCham
  Float_t fQ;               //QDC value, fractions are permitted for summable procedure  
  ClassDef(AliHMPIDDigit,4) //HMPID digit class       
};//class AliHMPIDDigitN

typedef AliHMPIDDigit AliRICHDigit; // for backward compatibility

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDDigit::AliHMPIDDigit(Int_t c,Float_t q,Int_t t,Float_t x,Float_t y,Int_t flag):AliDigit(),fPad(Abs(-1,-1,-1,-1)),fQ(-1)  
{
// Creation of sdigit  
// Arguments: c- chamber 
//            q- total QDC
//            t -TID
//            x,y - hit position, LORS        
//   Returns: none    
  Int_t pc,padx,pady;
  if     (x>=          0          && x<=  SizePcX()            ) {pc=0; padx=Int_t( x                           / SizePadX());}//PC 0 or 2 or 4
  else if(x>=SizePcX()+SizeDead() && x<=  SizeAllX()           ) {pc=1; padx=Int_t((x-  SizePcX()-  SizeDead()) / SizePadX());}//PC 2 or 4 or 6
  else return;
  if     (y>=          0          && y<=  SizePcY()            ) {      pady=Int_t( y                           / SizePadY());}//PC 0 or 1
  else if(y>=SizePcY()+SizeDead() && y<=2*SizePcY()+SizeDead() ) {pc+=2;pady=Int_t((y-  SizePcY()-  SizeDead()) / SizePadY());}//PC 2 or 3
  else if(y>=SizeAllY()-SizePcY() && y<=  SizeAllY()           ) {pc+=4;pady=Int_t((y-2*SizePcY()-2*SizeDead()) / SizePadY());}//PC 4 or 5
  else return;
  
  switch(flag){
    case kUL:padx--;pady++;break;    case kU:pady++;break;     case kUR:padx++; pady++;break;
                                              
    case kL: padx--;       break;    case kC:       break;     case kR:padx++;         break;
                                                 
    case kDL:padx--;pady--;break;    case kD:pady--;break;     case kDR:padx++; pady--;break;                                            
  }
  if(padx<0 || padx>=kPadPcX) return;
  if(pady<0 || pady>=kPadPcY) return;
  fPad=Abs(c,pc,padx,pady);
  fQ=q*Mathieson(x,y);
  fTracks[0]=t; 
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
  else                                        return -1;
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
Float_t AliHMPIDDigit::Mathieson(Float_t x,Float_t y)const
{
// This is the answer to electrostatic problem of charge distrubution in MWPC described elsewhere. (NIM A370(1988)602-603)
// Arguments: x,y- position of the center of Mathieson distribution
//  Returns: a charge fraction [0-1] imposed into the pad
  const Float_t kSqrtK3=0.77459667,k2=0.962,k4=0.379;

  Float_t ux1=kSqrtK3*TMath::TanH(k2*(x-LorsX()+0.5*SizePadX())/0.425);
  Float_t ux2=kSqrtK3*TMath::TanH(k2*(x-LorsX()-0.5*SizePadX())/0.425);
  Float_t uy1=kSqrtK3*TMath::TanH(k2*(y-LorsY()+0.5*SizePadY())/0.425);
  Float_t uy2=kSqrtK3*TMath::TanH(k2*(y-LorsY()-0.5*SizePadY())/0.425);
  return 4*k4*(TMath::ATan(ux2)-TMath::ATan(ux1))*k4*(TMath::ATan(uy2)-TMath::ATan(uy1));
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDDigit::Raw(UInt_t &w32)const
{
// Convert digit structure to raw word format
// Arguments: 32 bits raw word to fill
//   Returns: DDL ID where to write this digit
  w32=0;
  AliBitPacking::PackWord((UInt_t)fQ        ,w32, 0,11);  // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq        Qdc               bits (00..11) counts (0..4095)
  AliBitPacking::PackWord(         Addr()   ,w32,12,17);  // 3322 2222 2222 1111 1111 1000 0000 0000        DILOGIC address   bits (12..17) counts (0..47)
  AliBitPacking::PackWord(         Dilogic(),w32,18,21);  // 1098 7654 3210 9876 5432 1098 7654 3210        DILOGIC number    bits (18..21) counts (1..10)
  AliBitPacking::PackWord(         Row()    ,w32,22,26);  //                                                Row number        bits (22..26) counts (1..24)  
  return Ddl(); //ddl 0..13 where to write this digit 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigit::ReadRaw(Int_t ddl,UInt_t w32)
{
// Converts a given raw data word to a digit
// Arguments: w32 - 32 bits raw data word
//            ddl - DDL idx  0 1 2 3 4 ... 13
//   Returns: none
        fQ = AliBitPacking::UnpackWord(w32, 0,11);  // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq         Qdc               bits (00..11) counts (0..4095)   
  UInt_t a = AliBitPacking::UnpackWord(w32,12,17);  // 3322 2222 2222 1111 1111 1000 0000 0000         DILOGIC address   bits (12..17) counts (0..47)
  UInt_t d = AliBitPacking::UnpackWord(w32,18,21);  // 1098 7654 3210 9876 5432 1098 7654 3210         DILOGIC number    bits (18..21) counts (1..10)
  UInt_t r = AliBitPacking::UnpackWord(w32,22,26);  //                                                 Row number        bits (22..26) counts (1..24)    
  ReadRaw(ddl,r,d,a);    
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
