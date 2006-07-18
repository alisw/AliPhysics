#ifndef AliRICHDigit_h
#define AliRICHDigit_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliDigit.h>      //base class  
#include <AliBitPacking.h> //ToRaw(), FromRaw()
#include <TVector.h>       //ctor

class AliRICHDigit :public AliDigit //TObject-AliDigit-AliRICHDigit
{
public:
  enum EAbsPad {kChamAbs=10000000,kSecAbs=1000000,kPadAbsX=1000,kPadAbsY=1};                            //absolute pad number structure
  enum ERawData{kDiloX=8,kDiloY=6,kNdilo=10};                                                           //DILOGIC structure
  enum EPadData{kFirstPad=1,kPadsSecX=80,kPadsSecY=48,kPadsChamX=160,kPadsChamY=144,kSecX=2,kSecY=3};   //Segmentation structure 
  enum EDdlData{kNddls=14};                                                                             //DDL structure
//ctor&dtor    
  AliRICHDigit()                                                :AliDigit(),fCFM(-1) ,fChamber(-1  )     ,fPadX(-1)      ,fPadY(-1)      ,fQdc(-1)  {}
  AliRICHDigit(Int_t pad,Double_t qdc,Int_t cfm=-1,Int_t tid=-1):AliDigit(),fCFM(cfm),fChamber(P2C(pad)) ,fPadX(P2X(pad)),fPadY(P2Y(pad)),fQdc(qdc) {fTracks[0]=tid;}
  AliRICHDigit(TVector pad,Double_t q                          ):AliDigit(),fCFM(-1) ,fChamber(-1)       ,fPadX((Int_t)pad[0])  ,fPadY((Int_t)pad[1])  ,fQdc(q)   {}
  AliRICHDigit(Int_t c,TVector pad,Double_t q,Int_t cfm,Int_t tid0,Int_t tid1,Int_t tid2):fCFM(cfm),fChamber(c)  
       {fPadX=(Int_t)pad[0];fPadY=(Int_t)pad[1];fQdc=q;fTracks[0]=tid0;fTracks[1]=tid1;fTracks[2]=tid2;}
  virtual ~AliRICHDigit() {;}
//framework part    
                Bool_t   IsSortable  (                               )const{return kTRUE;}                                                    //provision to use TObject::Sort()
         inline Int_t    Compare     (const TObject *pObj            )const;                                                                  //provision to use TObject::Sort()
                void     Print       (Option_t *opt=""               )const;                                                                  //TObject::Print() overloaded
//private part  
                Int_t    A           (                               )const{return   (PadY()-1)%kDiloY*kDiloX+(PadX()-1)%kDiloX;}             //DILOGIC address 0..47 invented ?????
                void     AddTidOffset(Int_t offset                   )     {for (Int_t i=0; i<3; i++) if (fTracks[i]>0) fTracks[i]+=offset;}; //needed for merging
                Int_t    Cfm         (                               )const{return fCFM;}                                                     //ckov-feed-mip mixture
                Int_t    Chamber     (                               )const{return fChamber;}                                                 //chamber number
                Int_t    C           (                               )const{return fChamber;}                                                 //chamber number 
                Int_t    D           (                               )const{return 1+(PadX()-1)/kDiloX;}                                      //DILOGIC chip number 1..10
                Int_t    Ddl         (                               )const{return (PadX()<=kPadsSecX) ? 2*C()-2 : 2*C()-1;}                  //DDL number 0..13
         inline Int_t    Dig2Raw     (        UInt_t &w              )const;                                                                  //returns DDL ID and fill raw word
                Int_t    PadX        (                               )const{return fPadX;}                                                    //x position of the pad
                Int_t    PadY        (                               )const{return fPadY;}                                                    //y postion of the pad     
                TVector  Pad         (                               )const{Float_t v[2]={fPadX,fPadY}; return TVector(2,v);}
                Int_t    PadAbs      (                               )const{return fChamber*kChamAbs+fPadX*kPadAbsX+fPadY;}                   //absolute id of this pad
  static inline Int_t    Pad2Sec     (Int_t x,Int_t y                );                                   
  static        Int_t    P2A         (Int_t c,        Int_t x,Int_t y)     {Int_t s=Pad2Sec(x,y);return c*kChamAbs+s*kSecAbs+x*kPadAbsX+y;}   //(cham,padx,pady)-> abs pad
  static        Int_t    P2A         (Int_t c,Int_t s,Int_t x,Int_t y)     {return Pad2Sec(x,y)==s?c*kChamAbs+s*kSecAbs+x*kPadAbsX+y:-1;}     //(cham,sec,padx,pady)-> abs pad
  static        Int_t    P2C         (Int_t pad                      )     {return pad/kChamAbs;}                                             //abs pad -> chamber
  static        Int_t    P2S         (Int_t pad                      )     {return pad%kChamAbs/kSecAbs;}                                     //abs pad -> sector 
  static        Int_t    P2X         (Int_t pad                      )     {return pad%kSecAbs/kPadAbsX;}                                     //abs pad -> pad X 
  static        Int_t    P2Y         (Int_t pad                      )     {return pad%kPadAbsX;}                                             //abs pad number-> pad Y 
                Double_t Qdc         (                               )const{return fQdc;}                                                     //charge, QDC
                Int_t    R           (                               )const{return 1+(PadY()-1)/kDiloY;}                                      //DILOGIC row number 1..24
         inline void     Raw2Dig     (Int_t ddl,UInt_t  w32          );                                                                       //(DDL,w32)->(ch,sec,x,y,QDC)
                Int_t    S           (                               )const{return -1;}                                                       //sector number  ?????
  static        void     Test        (                               );                                                                       //test raw-digit manipulations
protected:
  Int_t    fCFM;      //1000000*Ncerenkovs+1000*Nfeedbacks+Nmips  
  Int_t    fChamber;  //chamber number
  Int_t    fPadX;     //pad along X counts from kFirstPad to kPadsChamX 
  Int_t    fPadY;     //pad along Y counts from kFirstPad to kPadsChamY 
  Double_t fQdc;      //QDC value, fractions are permitted for summable procedure  
  ClassDef(AliRICHDigit,3) //RICH digit class       
};//class AliRICHDigit
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHDigit::Compare(const TObject *pObj) const
{
// Used in Sort() method to compare to objects. Note that abs pad structure is first x then y, hence will be sorted on column basis.
// This feature is used in digitizer to facilitate finding of sdigits for the same pad since they all will come together after sorting.
// Arguments: pObj - pointer to object to compare with
//   Retunrs: -1 if AbsPad less then in pObj, 1 if more and 0 if they are the same      
  if     (PadAbs()==((AliRICHDigit*)pObj)->PadAbs()) return  0;
  else if(PadAbs() >((AliRICHDigit*)pObj)->PadAbs()) return  1;
  else                                               return -1;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHDigit::Raw2Dig(Int_t ddl,UInt_t w32)
{
//Converts a given raw data word to a digit
//Arguments: w32 - 32 bits raw data word
//           ddl - DDL file number  0 1 2 3 4 ... 13
//  Returns: none
      fQdc = AliBitPacking::UnpackWord(w32, 0,11);  // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq         Qdc               bits (00..11) counts (0..4095)   
  UInt_t a = AliBitPacking::UnpackWord(w32,12,17);  // 3322 2222 2222 1111 1111 1000 0000 0000         DILOGIC address   bits (12..17) counts (0..47)
  UInt_t d = AliBitPacking::UnpackWord(w32,18,21);  // 1098 7654 3210 9876 5432 1098 7654 3210         DILOGIC number    bits (18..21) counts (1..10)
  UInt_t r = AliBitPacking::UnpackWord(w32,22,26);  //                                                 Row number        bits (22..26) counts (1..24)    
  
  fPadY    = (r-1)*kDiloY+a/kDiloX+1;
  fPadX    = (d-1)*kDiloX+a%kDiloX+1;      fPadX+=(ddl%2)*kDiloX*kNdilo;//if ddl is odd then right half of the chamber
  fChamber = (ddl+2)/2;  // ddl 0..13 to chamber 1..7
}  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHDigit::Dig2Raw(UInt_t &w32)const
{
// Convert digit structure to raw word format
// Arguments: 32 bits raw word to fill
//   Returns: DDL ID where to write this digit
  w32=0;
  AliBitPacking::PackWord((UInt_t)fQdc,w32, 0,11);  // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq        Qdc               bits (00..11) counts (0..4095)
  AliBitPacking::PackWord(         A(),w32,12,17);  // 3322 2222 2222 1111 1111 1000 0000 0000        DILOGIC address   bits (12..17) counts (0..47)
  AliBitPacking::PackWord(         D(),w32,18,21);  // 1098 7654 3210 9876 5432 1098 7654 3210        DILOGIC number    bits (18..21) counts (1..10)
  AliBitPacking::PackWord(         R(),w32,22,26);  //                                                Row number        bits (22..26) counts (1..24)  
  return Ddl(); //ddl 0..13 where to write this digit 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHDigit::Pad2Sec(Int_t padx,Int_t pady)
{
// Determines sector containing the given pad.
// Arguments: padx,pady - pad number
//   Returns: sector number    
// y ^  5 6         sectors map as seen from outside (IP is behind the page)
//   |  3 4
//   |  1 2
//   -------> x  
  Int_t sector=-1;      
  if     (padx >= kFirstPad && padx <=   kPadsSecX )    sector=1;
  else if(padx >  kPadsSecX && padx <=   kPadsChamX)    sector=2; 
  else                                                  return -1;//padx out of range     
  if     (pady >= kFirstPad && pady <=   kPadsSecY )    return sector;
  else if(pady >  kPadsSecY && pady <= 2*kPadsSecY )    return sector+2;
  else if(pady >2*kPadsSecY && pady <=   kPadsChamY)    return sector+4;
  else                                                  return -1; //pady out of range
}//Pad2Sec()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
