#ifndef AliHMPIDDigitN_h
#define AliHMPIDDigitN_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliDigit.h>      //base class  

class AliHMPIDDigitN :public AliDigit //TObject-AliDigit-AliHMPIDDigit
{
public:
  enum EAbsPad {kChAbs=100000000,kPcAbs=1000000,kPadAbsX=1000,kPadAbsY=1};       //absolute pad number structure
  enum ERawData{kDilX=8,kDilY=6,kNdil=10,kNrow=24,kNddls=14};                    //RAW data structure
  enum EPadData{kPadsPcX=80,kPadsPcY=48,kPadsChamX=160,kPadsChamY=144,kNpc=6};   //Segmentation structure 
//ctor&dtor    
           AliHMPIDDigitN(                                     ):AliDigit(),fPad(Abs(-1,-1,-1,-1)) ,fQdc(-1)  {} //default ctor
  inline   AliHMPIDDigitN(Int_t c,Float_t x,Float_t y,Float_t q);                                                //ctor 
  virtual ~AliHMPIDDigitN()                                {} //dtor
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
  static void    DrawPc      (                               );                                                                        //draw PCs
         Int_t   Ddl         (                               )const{return (PadX()<kPadsPcX) ? 2*Ch()-2 : 2*Ch()-1;}                  //DDL number 0..13
         Float_t LorsX       (                               )const{return (PadX()+0.5)*SizePadX()+(Pc()%2)*(SizePcX()+SizeDead());}   //center of the pad x, [cm]
         Float_t LorsY       (                               )const{return (PadY()+0.5)*SizePadY()+(Pc()/2)*(SizePcY()+SizeDead());}   //center of the pad y, [cm]
         Int_t   PadX        (                               )const{return A2X(fPad);}                                                 //x position of the pad
         Int_t   PadY        (                               )const{return A2Y(fPad);}                                                 //y postion of the pad     
         Int_t   Pad         (                               )const{return fPad;}                                                      //absolute id of this pad
         Float_t Qdc         (                               )const{return fQdc;}                                                      //charge, [QDC]
         Int_t   Pc          (                               )const{return A2P(fPad);}                                                 //PC position number
  static void    PrintSize   (                               );                                                                        //print all segmentation sizes      
         Int_t   Row         (                               )const{Int_t r=1+Pc()/2*8+PadY()/kDilY; return (Pc()%2)?kNrow-r+1:r;}     //row r=1..24
         void    Set         (Int_t c,Int_t s,Int_t x,Int_t y)     {fPad=Abs(c,s,x,y);}                                                //set new digit
         void    ReadRaw     (Int_t ddl,Int_t r,Int_t d,Int_t a){Int_t mapA2Y[kDilY]={3,2,4,1,5,0};fPad=Abs(ddl/2,ddl%7,d*kDilX+a/kDilY,r*kDilY+mapA2Y[a%kDilY]);} //from raw
  static Float_t SizePadX    (                               )     {return 0.8;}                                                       //pad size x, [cm]        
  static Float_t SizePadY    (                               )     {return 0.84;}                                                      //pad size y, [cm]  
  static Float_t SizePcX     (                               )     {return SizePadX()*kPadsPcX;}                                      //PC size x, [cm]        
  static Float_t SizePcY     (                               )     {return SizePadY()*kPadsPcY;}                                      //PC size y, [cm]    
  static Float_t SizeAllX    (                               )     {return SizePadX()*kPadsChamX+SizeDead();}                          //all PCs size x, [cm]        
  static Float_t SizeAllY    (                               )     {return SizePadY()*kPadsChamY+2*SizeDead();}                        //all PCs size y, [cm]    
  static Float_t SizeDead    (                               )     {return 2.6;}                                                       //dead zone size x, [cm]         
  static void    TestSeg     (                               );                                                                        //test segmentation
         void    Zoom        (                               ); 
protected:
  Int_t   fPad;             //absolute pad number is chamber*kCham
  Float_t fQdc;             //QDC value, fractions are permitted for summable procedure  
  ClassDef(AliHMPIDDigitN,4) //HMPID digit class       
};//class AliHMPIDDigitN
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDDigitN::Compare(const TObject *pObj) const
{
// Used in Sort() method to compare to objects. Note that abs pad structure is first x then y, hence will be sorted on column basis.
// This feature is used in digitizer to facilitate finding of sdigits for the same pad since they all will come together after sorting.
// Arguments: pObj - pointer to object to compare with
//   Retunrs: -1 if AbsPad less then in pObj, 1 if more and 0 if they are the same      
  if     (fPad==((AliHMPIDDigitN*)pObj)->Pad()) return  0;
  else if(fPad >((AliHMPIDDigitN*)pObj)->Pad()) return  1;
  else                                         return -1;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDDigitN::AliHMPIDDigitN(Int_t c,Float_t x,Float_t y,Float_t q):AliDigit(),fPad(Abs(-1,-1,-1,-1)),fQdc(-1)
{
// Determines abs pad number containing the given point (x,y) defined in the chamber RS.
// Pad count starts in lower left corner from 1,1 to 144,160 in upper right corner of a chamber.
// y ^  4 5
//   |  2 3
//   |  0 1
//    -------> x
  Int_t pc,padx,pady;
  if     (x>=          0          && x<=  SizePcX()            ) {pc=0; padx=Int_t( x                           / SizePadX());}//PC 0 or 2 or 4
  else if(x>=SizePcX()+SizeDead() && x<=  SizeAllX()           ) {pc=1; padx=Int_t((x-  SizePcX()-  SizeDead()) / SizePadX());}//PC 2 or 4 or 6
  else return;
  if     (y>=          0          && y<=  SizePcY()            ) {      pady=Int_t( y                           / SizePadY());}//PC 0 or 1
  else if(y>=SizePcY()+SizeDead() && y<=2*SizePcY()+SizeDead() ) {pc+=2;pady=Int_t((y-  SizePcY()-  SizeDead()) / SizePadY());}//PC 2 or 3
  else if(y>=SizeAllY()-SizePcY() && y<=  SizeAllY()           ) {pc+=4;pady=Int_t((y-2*SizePcY()-2*SizeDead()) / SizePadY());}//PC 4 or 5
  else return;
  fPad=Abs(c,pc,padx,pady);
  fQdc=q;
}
#endif
