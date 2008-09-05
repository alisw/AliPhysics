#ifndef AliHMPIDDigit_h
#define AliHMPIDDigit_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class of HMPID to manage digits ---> pads
//.
//.
//.

#include <AliDigit.h>      //base class  
#include <AliRawReader.h>
#include <AliLog.h>
#include "TMath.h"         //Mathieson()
#include <AliBitPacking.h> //Raw()
#include "AliHMPIDParam.h"
//#include "AliHMPIDRawStream.h"

class TClonesArray;        //Hit2Sdi()
  
class AliHMPIDDigit :public AliDigit //TObject-AliDigit-AliHMPIDDigit
{
public:
    
//ctor&dtor    
  AliHMPIDDigit(                          ):AliDigit( ),fPad(AliHMPIDParam::Abs(-1,-1,-1,-1)),fQ(-1)  {}                         //default ctor
  AliHMPIDDigit(Int_t pad,Int_t q,Int_t *t):AliDigit(t),fPad(pad             ),fQ(q )  {if(fQ>4095)fQ=4095;}                     //digit ctor
  AliHMPIDDigit(Int_t pad,Int_t q         ):AliDigit( ),fPad(pad             ),fQ(q )  {if(fQ>4095)fQ=4095;}                     //digit ctor
  AliHMPIDDigit(const AliHMPIDDigit &d    ):AliDigit(d),fPad(d.fPad),fQ(d.fQ)          {}                                        //copy ctor
  virtual ~AliHMPIDDigit()                                                             {}                         //dtor   
//framework part    
         Bool_t  IsSortable  (                               )const{return kTRUE;}                                                     //provision to use TObject::Sort() 
  inline Int_t   Compare     (const TObject *pObj            )const;                                                                   //provision to use TObject::Sort()
         void    Draw        (Option_t *opt=""               );                                                                        //TObject::Draw() overloaded
         void    Print       (Option_t *opt=""               )const;                                                                   //TObject::Print() overloaded
//private part  

         void    AddTidOffset(Int_t offset                   )     {for (Int_t i=0; i<3; i++) if (fTracks[i]>0) fTracks[i]+=offset;  } //needed for merging
         Int_t   Ch          (                               )const{return AliHMPIDParam::A2C(fPad);                                                } //chamber number
 
         Float_t LorsX       (                               )const{return AliHMPIDParam::LorsX(AliHMPIDParam::A2P(fPad),AliHMPIDParam::A2X(fPad));                               } //center of the pad x, [cm]

         Float_t LorsY       (                               )const{return AliHMPIDParam::LorsY(AliHMPIDParam::A2P(fPad),AliHMPIDParam::A2Y(fPad));                               } //center of the pad y, [cm]
//  
  inline Double_t MathiesonX   (Double_t x                   )const;                                                                   //Mathieson distribution along wires X 
  inline Double_t MathiesonY   (Double_t x                   )const;                                                                   //Mathieson distribution perp to wires Y
  inline Double_t IntPartMathiX(Double_t z                   )const;                                                                   //integral in 1-dim of Mathieson X
  inline Double_t IntPartMathiY(Double_t z                   )const;                                                                   //integral in 1-dim of Mathieson Y
  inline Double_t IntMathieson (Double_t x,Double_t y        )const;                                                                   //integral in 2-dim of Mathieson  
         Int_t   PadPcX      (                               )const{return AliHMPIDParam::A2X(fPad);}                                                 //pad pc x # 0..79
         Int_t   PadPcY      (                               )const{return AliHMPIDParam::A2Y(fPad);}                                                 //pad pc y # 0..47
         Int_t   PadChX      (                               )const{return (Pc()%2)*AliHMPIDParam::kPadPcX+PadPcX();}                                 //pad ch x # 0..159
         Int_t   PadChY      (                               )const{return (Pc()/2)*AliHMPIDParam::kPadPcY+PadPcY();}                                 //pad ch y # 0..143
         Int_t   Pad         (                               )const{return fPad;}                                                      //absolute id of this pad
         Int_t   Pc          (                               )const{return AliHMPIDParam::A2P(fPad);}                                                 //PC position number
         Float_t Q           (                               )const{return fQ;}                                                        //charge, [QDC]
  inline void    Raw(UInt_t &w32,Int_t &ddl,Int_t &r,Int_t &d,Int_t &a)const;
  inline Bool_t  Set         (Int_t c,Int_t p,Int_t x,Int_t y,Int_t tid=0);                                                            //manual creation 
         void    SetQ        (Float_t q                      )     {fQ=q;if(fQ>4095)fQ=4095;}                                          //setter for charge 
         void    SetPad      (Int_t pad                      )     {fPad=pad;}                                                         //setter for pad
 
protected:                                                                   //AliDigit has fTracks[3]
                                                                               

  Int_t    fPad;                                                                                                                       //absolute pad number
  Float_t  fQ;                                                               //QDC value, fractions are permitted for summable procedure  
  ClassDef(AliHMPIDDigit,4)                                                  //HMPID digit class       
};//class AliHMPIDDigit

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

Double_t AliHMPIDDigit::MathiesonX(Double_t x)const
{
// Mathieson function.
// This is the answer to electrostatic problem of charge distrubution in MWPC described elsewhere. (NIM A370(1988)602-603)
// Arguments: x- position of the center of Mathieson distribution
//  Returns: value of the Mathieson function
  
  Double_t lambda = x/AliHMPIDParam::PitchAnodeCathode();
  Double_t tanh = TMath::TanH(AliHMPIDParam::K2x()*lambda);
  Double_t a=1-tanh*tanh;
  Double_t b=1+AliHMPIDParam::SqrtK3x()*AliHMPIDParam::SqrtK3x()*tanh*tanh;
  Double_t mathi = AliHMPIDParam::K1x()*a/b;
  return mathi;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t AliHMPIDDigit::MathiesonY(Double_t y)const
{
// Mathieson function.
// This is the answer to electrostatic problem of charge distrubution in MWPC described elsewhere. (NIM A370(1988)602-603)
// Arguments: x- position of the center of Mathieson distribution
//  Returns: value of the Mathieson function
  
  Double_t lambda = y/AliHMPIDParam::PitchAnodeCathode();
  Double_t tanh = TMath::TanH(AliHMPIDParam::K2y()*lambda);
  Double_t a=1-tanh*tanh;
  Double_t b=1+AliHMPIDParam::SqrtK3y()*AliHMPIDParam::SqrtK3y()*tanh*tanh;
  Double_t mathi = AliHMPIDParam::K1y()*a/b;
  return mathi;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t AliHMPIDDigit::IntPartMathiX(Double_t x)const
{
// Integration of Mathieson.
// This is the answer to electrostatic problem of charge distrubution in MWPC described elsewhere. (NIM A370(1988)602-603)
// Arguments: x,y- position of the center of Mathieson distribution
//  Returns: a charge fraction [0-1] imposed into the pad
  Double_t shift1 = -LorsX()+0.5*AliHMPIDParam::SizePadX();
  Double_t shift2 = -LorsX()-0.5*AliHMPIDParam::SizePadX();
    
  Double_t ux1=AliHMPIDParam::SqrtK3x()*TMath::TanH(AliHMPIDParam::K2x()*(x+shift1)/AliHMPIDParam::PitchAnodeCathode());
  Double_t ux2=AliHMPIDParam::SqrtK3x()*TMath::TanH(AliHMPIDParam::K2x()*(x+shift2)/AliHMPIDParam::PitchAnodeCathode());
  
  return AliHMPIDParam::K4x()*(TMath::ATan(ux2)-TMath::ATan(ux1));
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t AliHMPIDDigit::IntPartMathiY(Double_t y)const
{
// Integration of Mathieson.
// This is the answer to electrostatic problem of charge distrubution in MWPC described elsewhere. (NIM A370(1988)602-603)
// Arguments: x,y- position of the center of Mathieson distribution
//  Returns: a charge fraction [0-1] imposed into the pad
  Double_t shift1 = -LorsY()+0.5*AliHMPIDParam::SizePadY();
  Double_t shift2 = -LorsY()-0.5*AliHMPIDParam::SizePadY();
    
  Double_t uy1=AliHMPIDParam::SqrtK3y()*TMath::TanH(AliHMPIDParam::K2y()*(y+shift1)/AliHMPIDParam::PitchAnodeCathode());
  Double_t uy2=AliHMPIDParam::SqrtK3y()*TMath::TanH(AliHMPIDParam::K2y()*(y+shift2)/AliHMPIDParam::PitchAnodeCathode());
  
  return AliHMPIDParam::K4y()*(TMath::ATan(uy2)-TMath::ATan(uy1));
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t AliHMPIDDigit::IntMathieson(Double_t x,Double_t y)const
{
// Integration of Mathieson.
// This is the answer to electrostatic problem of charge distrubution in MWPC described elsewhere. (NIM A370(1988)602-603)
// Arguments: x,y- position of the center of Mathieson distribution
//  Returns: a charge fraction [0-1] imposed into the pad

  Double_t xm = IntPartMathiX(x);
  Double_t ym = IntPartMathiY(y);
  return 4*xm*ym;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDDigit::Raw(UInt_t &w32,Int_t &ddl,Int_t &r,Int_t &d,Int_t &a)const
{
// Convert digit structure to raw word format
// Arguments: w32,ddl,r,d,a where to write the results
//   Returns: none
  Int_t y2a[6]={5,3,1,0,2,4};

  ddl=2*Ch()+Pc()%2;                    //DDL# 0..13
  Int_t tmp=1+Pc()/2*8+PadPcY()/6;  r=(Pc()%2)? tmp:25-tmp;               //row r=1..24
  d=1+PadPcX()/8;                       //DILOGIC# 1..10
//  d=AliHMPIDRawStream::kNDILOGICAdd+1-d;                   ////flip according to Paolo (2-9-2008)
  d=10+1-d;                                                  ////flip according to Paolo (2-9-2008)
  a=y2a[PadPcY()%6]+6*(7-PadPcX()%8);   //ADDRESS 0..47        
    
  w32=0;   
  //Printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
  //
  //Printf("AliHMPIDDigit::Raw ddl: %d r: %d d: %d a: %d",ddl,r,d,a);
  //Printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
//  Bool_t isOK=kTRUE; isOK=
  AliBitPacking::PackWord((UInt_t)fQ,w32, 0,11);                      // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq        Qdc               bits (00..11) counts (0..4095)
  assert(0<=a&&a<=47);AliBitPacking::PackWord(        a ,w32,12,17);  // 3322 2222 2222 1111 1111 1000 0000 0000        DILOGIC address   bits (12..17) counts (0..47)
  assert(1<=d&&d<=10);AliBitPacking::PackWord(        d ,w32,18,21);  // 1098 7654 3210 9876 5432 1098 7654 3210        DILOGIC number    bits (18..21) counts (1..10)
  assert(1<=r&&r<=24);AliBitPacking::PackWord(        r ,w32,22,26);  //                                                Row number        bits (22..26) counts (1..24)  
                      AliBitPacking::PackWord((UInt_t)0, w32,27,27);  //To make sure set the 27th bit to Zero so we can distinguis it from the EoE
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDDigit::Set(Int_t ch,Int_t pc,Int_t px,Int_t py,Int_t tid)
{
// Manual creation of digit
// Arguments: ch,pc,px,py,qdc,tid  
//   Returns: kTRUE if wrong digit
  if(px<AliHMPIDParam::kMinPx || px>AliHMPIDParam::kMaxPx) return kTRUE;
  if(py<AliHMPIDParam::kMinPy || py>AliHMPIDParam::kMaxPy) return kTRUE;

  fPad=AliHMPIDParam::Abs(ch,pc,px,py);fTracks[0]=tid;
  fQ=0;
  return kFALSE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#endif
