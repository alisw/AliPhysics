#ifndef AliRICHDigit_h
#define AliRICHDigit_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliDigit.h>      //base class  
#include <AliBitPacking.h> //Dig2Raw()
#include "AliRICHParam.h"

//RICH DDL ID allowed range is [0x700,0x714] or in decimal notation [1792,1812]. 20 DDL files are reserved. (kDdlOffset)
//RICH actually uses 14 DDLs (kNddls), 2 per chamber, even number for left part(1-3-5) odd number for right part(2-4-6) 
//So the chamber-DDL map is:
//N  0 1L=0x700 1792          N  1 1R=0x701 1793
//N  2 2L=0x702 1794          N  3 2R=0x703 1795
//N  4 3L=0x704 1796          N  5 3R=0x705 1797
//N  6 4L=0x706 1798          N  7 4R=0x707 1799
//N  8 5L=0x708 1800          N  9 5R=0x709 1801
//N 10 6L=0x70a 1802          N 11 6R=0x70b 1803
//N 12 7L=0x70c 1804          N 13 7R=0x70d 1805
//RICH has no any propriate header just uses the common one
//RICH chamber is divide on 2 halves vertically 
//Half chamber is divided by 24 rows counted from 1 to 24 (8 raws per sector) from top to bottom for left half chamber (sectors 1-3-5) 
//                                                                         and from bottom to top for right half chamber (sectors 2-4-6) as seen from MARS (0,0,0)
//Raw is composed from 10 DILOGIC chips (kNchips) counted from left to right from 1 to 10  as seen from MARS (0,0,0)
//So each DILOGIC chip serves 48 channels for the 8x6 pads box (kChipX,kChipY). Channels counted from 0 to 47.
//??????? Currently the exact mapping of DILOGIC addresses to pads is not known. So we invented horizontal zig-zag ???????  
//So RICH raw word is  32 bits word with structure:   
// 00000        rrrrr                      dddd                               aaaaaa                          qqqqqqqqqqqq 
// 5 bits zero  5 bits raw number (1..24)  4 bits DILOGIC chip number (1..10) 6 bits DILOGIC address (0..47)  12 bits QDC value (0..4095)

class AliRICHDigit :public AliDigit
{
public:
  enum EAbsPad {kChamber=10000000,kPadX=1000};                           //absolute pad number structure
  enum ERawProp{kChipX=8,kChipY=6,kNchips=10,kNddls=14,kRichRawId=7,kDdlOffset=0x700};//DILOGIC is 8x6 pads
  AliRICHDigit()                                  :AliDigit(),fCFM(-1),fChamber(-1  )  ,fPadX(-1)    ,fPadY(-1)    ,fQdc(-1) {}
  AliRICHDigit(Int_t c,Int_t x,Int_t y,Double_t q):AliDigit(),fCFM(-1),fChamber(10*c)  ,fPadX(x )    ,fPadY(y )    ,fQdc(q ) {}
  AliRICHDigit(Int_t c,TVector pad,Double_t q,Int_t cfm,Int_t tid0,Int_t tid1,Int_t tid2):fCFM(cfm)  
       {fPadX=(Int_t)pad[0];fPadY=(Int_t)pad[1];fQdc=q;fChamber=10*c+AliRICHParam::Pad2Sec(pad);fTracks[0]=tid0;fTracks[1]=tid1;fTracks[2]=tid2;}
  virtual ~AliRICHDigit() {;}
//framework part    
         Bool_t   IsSortable  (                   )const{return kTRUE;}                                                    //provision to use TObject::Sort()
  inline Int_t    Compare     (const TObject *pObj)const;                                                                  //provision to use TObject::Sort()
         void     Print       (Option_t *option="")const;                                                                  //TObject Print() overload
//private part  
         void     AddTidOffset(Int_t offset     )     {for (Int_t i=0; i<3; i++) if (fTracks[i]>0) fTracks[i]+=offset;}; //needed for merging
         Int_t    Cfm         (                 )const{return fCFM;}                                                     //particle mixture for this digit
         Int_t    Chamber     (                 )const{return fChamber/10;}                                              //chamber number
         Int_t    Sector      (                 )const{return fChamber%10;}                                              //sector number
         Int_t    PadX        (                 )const{return fPadX;}                                                    //x position of the pad
         Int_t    PadY        (                 )const{return fPadY;}                                                    //y postion of the pad     
         TVector  Pad         (                 )const{Float_t v[2]={fPadX,fPadY}; return TVector(2,v);}
         Int_t    PadAbs      (                 )const{return fChamber*kChamber+fPadX*kPadX+fPadY;}                      //absolute id of this pad
         Double_t Qdc         (                 )const{return fQdc;}                                                     //charge in terms of ADC channels
  inline Int_t    Dig2Raw     (        UInt_t &w)const;                                                                  //returns DDL ID and fill raw 32 bits word
  inline void     Raw2Dig     (Int_t d,UInt_t  w);                                                                       //(DDL,word32)->(ch,sec,padx,pady,QDC)
  static Int_t    P2C         (Int_t pad        )     {return pad/kChamber;}                                             //abs pad number-> chamber number
  static Int_t    P2X         (Int_t pad        )     {return pad%kChamber/kPadX;}                                       //abs pad number-> pad X number
  static Int_t    P2Y         (Int_t pad        )     {return pad%kChamber%kPadX;}                                       //abs pad number-> pad Y number
         void     Test        (                 );                                                                       //used to test all possible digit manipulations
protected:
  Int_t    fCFM;  //1000000*Ncerenkovs+1000*Nfeedbacks+Nmips  
  Int_t    fChamber;  //10*chamber number+ sector number 
  Int_t    fPadX;     //pad number along X
  Int_t    fPadY;     //pad number along Y
  Double_t fQdc;      //QDC value, fractions are permitted for summable procedure  
  ClassDef(AliRICHDigit,3) //RICH digit class       
};//class AliRICHDigit
//__________________________________________________________________________________________________
Int_t AliRICHDigit::Compare(const TObject *pObj) const
{
//Used in Sort() method to compare to objects. Note that abs pad structure is first x then y, hence will be sorted on column basis.
//This feature is used in digitizer to facilitate finding of sdigits for the same pad as they will be together after sorting.
//Arguments: pObj - pointer to object to compare with
//  Retunrs: -1 if AbsPad less then in pObj, 1 if more and 0 if they are the same      
  if(PadAbs()==((AliRICHDigit*)pObj)->PadAbs())
    return 0;
  else if(PadAbs()>((AliRICHDigit*)pObj)->PadAbs())
    return 1;
  else 
    return -1;
}
//__________________________________________________________________________________________________
void AliRICHDigit::Raw2Dig(Int_t ddl,UInt_t w32)
{
//Reads next raw word from raw data stream and convert     
//Arguments: w32 - 32 bits word as in raw data stream
//           ddl - DDL file number  0 1 2 3 4 ... 13
//  Returns: none
      fQdc = AliBitPacking::UnpackWord(w32, 0,11);  // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq 
  UInt_t a = AliBitPacking::UnpackWord(w32,12,17);  // 3322 2222 2222 1111 1111 1000 0000 0000 
  UInt_t d = AliBitPacking::UnpackWord(w32,18,21);  // 1098 7654 3210 9876 5432 1098 7654 3210 
  UInt_t r = AliBitPacking::UnpackWord(w32,22,26);  // r- iRawN d- iChiN a- iChiC              
  
  fPadY    = (r-1)*kChipY+a/kChipX+1;
  fPadX    = (d-1)*kChipX+a%kChipX+1;      fPadX+=(ddl%2)*kChipX*kNchips;//if ddl is odd then right half of the chamber
  TVector pad(2); pad[0]=fPadX;pad[1]=fPadY;
  fChamber = ((ddl+2)/2)*10+AliRICHParam::Pad2Sec(pad);  // ddl 0..13 to chamber 1..7
}  
//__________________________________________________________________________________________________
Int_t AliRICHDigit::Dig2Raw(UInt_t &w32)const
{
//Convert digit structure to raw word format
//Arguments: 32 bits raw word to fill
//  Returns: DDL ID where to write this digit
  Int_t ddl=2*Chamber()-1;              //chamber 1..7 -> DDL 0..13, this idDdl is for right half (sectors 2 4 6), to be decremented if d < kNchips
  UInt_t a =  (PadY()-1)%kChipY*kChipX+(PadX()-1)%kChipX;          //invented to be horizontal zig-zag
  UInt_t r =1+(PadY()-1)/kChipY;      
  UInt_t d =1+(PadX()-1)/kChipX;    
  if(d>kNchips)
    d-=kNchips;              //chip number more then kNchips means right half of chamber, goes to this ddl
  else
    ddl--;                   //chip number less then kNchips means left half of the chamber, goes to ddl-1 
  
  w32=0;
  AliBitPacking::PackWord((UInt_t)fQdc,w32, 0,11);  // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq
  AliBitPacking::PackWord(           a,w32,12,17);  // 3322 2222 2222 1111 1111 1000 0000 0000
  AliBitPacking::PackWord(           d,w32,18,21);  // 1098 7654 3210 9876 5432 1098 7654 3210 
  AliBitPacking::PackWord(           r,w32,22,26);  
  return ddl; //ddl 0..13 where to write this digit 
}

#endif
