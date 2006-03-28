#ifndef AliRICHDigit_h
#define AliRICHDigit_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliDigit.h>      //base class  
#include <AliBitPacking.h> //ToRaw(), FromRaw()
#include <TVector.h>       //ctor
/*
Any given LDC collects data from all FEE connected to the LDC by DDL. This data is stored in name.ddl file.
Name of this file is composed by detector name plus some value for example RICH1793.ddl
That value is calculated as sequensial number of detector LDC plus some predifined offset, unique for a given detector.
The value is expected to be within a given range assigned to detector.
For RICH, the offset number is 0x700 hex or 1792 decimal (reffered in the code as kDdlOffset). 
The range assigned for RICH is 0x700-0x714 hex or 1792-1812 decimal. It is 20 ddl files or 20 LDCs.
RICH actually uses 14 LDCs hence DAQ writes for RICH 14 ddl files (reffered in the code as kNddls).
RICH FEE is connected to LDC in the following way:
Single LDC serves one half of a chamber i.e. 3 photocathodes aka sectors,  even LDC for left part( sectors 1-3-5) and odd LDC for right part(2-4-6) 
So the LDC -chamber-ddl file name map is:
LDC  0 -> ch 1L -> file name value 0x700 1792          LDC  1 -> ch 1R -> file name value 0x701 1793
LDC  2 -> ch 2L -> file name value 0x702 1794          LDC  3 -> ch 2R -> file name value 0x703 1795
LDC  4 -> ch 3L -> file name value 0x704 1796          LDC  5 -> ch 3R -> file name value 0x705 1797
LDC  6 -> ch 4L -> file name value 0x706 1798          LDC  7 -> ch 4R -> file name value 0x707 1799
LDC  8 -> ch 5L -> file name value 0x708 1800          LDC  9 -> ch 5R -> file name value 0x709 1801
LDC 10 -> ch 6L -> file name value 0x70a 1802          LDC 11 -> ch 6R -> file name value 0x70b 1803
LDC 12 -> ch 7L -> file name value 0x70c 1804          LDC 13 -> ch 7R -> file name value 0x70d 1805

Programmatically, operations with ddl files are interfaced by class AliRawReader. In order to select some ddl files for detector,
one needs to provide a reserved id number of detector. For RICH, this number is 7 (reffered in the code as kRichRawId).

DDL file starts with common header described in AliRawDataHeader, which can be followed by private header.
RICH has no any private header, just uses the common one.

RICH FEE as seen by single LDC is composed from a number of DILOGIC chips organized in vertical stack of rows. 
Each DILOGIC chip serves 48 channels for the 8x6 pads (reffered in the code as kDiloX,kDiloY). Channels counted from 0 to 47.

??????? Currently the exact mapping of DILOGIC addresses to pads is not known. So we invented horizontal zig-zag ???????  

10 DILOGIC chips composes so called "row" in horizontal direction (reffered in the code as kNdilo), so the row is 80x6 pads structure. 
DILOGIC chips in the row are counted from left to right as seen from MARS (0,0,0), from 1 to 10.
24 rows are piled up forming the whole FEE served by single LDC, so one LDC sees 80x144 pads separated in 3 photocathodes aka sectors.
Rows are counted from 1 to 24 from top    to bottom for left  half chamber (sectors 1-3-5) as seen from MARS (0,0,0), meaning even LDC number
                          and from bottom to top    for right half chamber (sectors 2-4-6) as seen from MARS (0,0,0), meaning odd LDC number.

So RICH raw word is 32 bits with the structure:   
   00000             rrrrr                      dddd                               aaaaaa                          qqqqqqqqqqqq 
 5 bits zero  5 bits row number (1..24)  4 bits DILOGIC chip number (1..10) 6 bits DILOGIC address (0..47)  12 bits QDC value (0..4095)
*/

class AliRICHDigit :public AliDigit
{
public:
  enum EAbsPad {kChamAbs=10000000,kSecAbs=1000000,kPadAbsX=1000,kPadAbsY=1};                //absolute pad number structure
  enum ERawData{kDiloX=8,kDiloY=6,kNdilo=10};                                               //DILOGIC structure, see description above 
  enum EPadData{kFirstPad=1,kPadsSecX=80,kPadsSecY=48,kPadsChamX=160,kPadsChamY=144,kSecX=2,kSecY=3};   //Segmentation structure 
  enum EDdlData{kNddls=14,kDdlOffset=0x700,kRichRawId=7};                                   //Common DDL structure, see description above
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
                void     Print       (Option_t *opt=""               )const;                                                                  //TObject Print() overload
//private part  
                void     AddTidOffset(Int_t offset                   )     {for (Int_t i=0; i<3; i++) if (fTracks[i]>0) fTracks[i]+=offset;}; //needed for merging
                Int_t    Cfm         (                               )const{return fCFM;}                                                     //ckov-feed-mip mixture
                Int_t    Chamber     (                               )const{return fChamber;}                                                 //chamber number
                Int_t    C           (                               )const{return fChamber;}                                                 //chamber number 
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
         inline void     Raw2Dig     (Int_t ddl,UInt_t  w32          );                                                                       //(DDL,w32)->(ch,sec,x,y,QDC)
                Int_t    S           (                               )const{return -1;}                                                       //sector number  ?????
  static        void     Test        (                               );                                                                       //test raw-digit manipulations
protected:
  Int_t    fCFM;      //1000000*Ncerenkovs+1000*Nfeedbacks+Nmips  
  Int_t    fChamber;  //chamber
  Int_t    fPadX;     //pad along X
  Int_t    fPadY;     //pad along Y
  Double_t fQdc;      //QDC value, fractions are permitted for summable procedure  
  ClassDef(AliRICHDigit,3) //RICH digit class       
};//class AliRICHDigit
//__________________________________________________________________________________________________
Int_t AliRICHDigit::Compare(const TObject *pObj) const
{
//Used in Sort() method to compare to objects. Note that abs pad structure is first x then y, hence will be sorted on column basis.
//This feature is used in digitizer to facilitate finding of sdigits for the same pad since they all will come together after sorting.
//Arguments: pObj - pointer to object to compare with
//  Retunrs: -1 if AbsPad less then in pObj, 1 if more and 0 if they are the same      
  if     (PadAbs()==((AliRICHDigit*)pObj)->PadAbs()) return  0;
  else if(PadAbs() >((AliRICHDigit*)pObj)->PadAbs()) return  1;
  else                                               return -1;
}
//__________________________________________________________________________________________________
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
//__________________________________________________________________________________________________
Int_t AliRICHDigit::Dig2Raw(UInt_t &w32)const
{
//Convert digit structure to raw word format
//Arguments: 32 bits raw word to fill
//  Returns: DDL ID where to write this digit
  Int_t ddl=2*C()-1;                                         //chamber 1..7 -> DDL 0..13, this idDdl is for right half (sectors 2 4 6), to be decremented if d < kNchips
  UInt_t a =  (PadY()-1)%kDiloY*kDiloX+(PadX()-1)%kDiloX;    //invented to be horizontal zig-zag
  UInt_t r =1+(PadY()-1)/kDiloY;                             //Row number depends only on y and we have (1..24) rows per (1..144) pads
  UInt_t d =1+(PadX()-1)/kDiloX;                             //DILOGIC number depends only on x we have (1..10) chips per (1..80) pads
  if(d>kNdilo)
    d-=kNdilo;              //chip number more then kNdilo means right half of chamber, goes to this ddl
  else
    ddl--;                  //chip number less then kNdilo means left half of the chamber, goes to ddl-1 
  
  w32=0;
  AliBitPacking::PackWord((UInt_t)fQdc,w32, 0,11);  // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq        Qdc               bits (00..11) counts (0..4095)
  AliBitPacking::PackWord(           a,w32,12,17);  // 3322 2222 2222 1111 1111 1000 0000 0000        DILOGIC address   bits (12..17) counts (0..47)
  AliBitPacking::PackWord(           d,w32,18,21);  // 1098 7654 3210 9876 5432 1098 7654 3210        DILOGIC number    bits (18..21) counts (1..10)
  AliBitPacking::PackWord(           r,w32,22,26);  //                                                Row number        bits (22..26) counts (1..24)  
  return ddl; //ddl 0..13 where to write this digit 
}
//__________________________________________________________________________________________________
Int_t AliRICHDigit::Pad2Sec(Int_t padx,Int_t pady)
{
//Determines sector containing the given pad.
//Arguments: padx,pady - pad number
//  Returns: sector number    
//y ^  5 6         sectors map as seen from IP
//  |  3 4
//  |  1 2
//   -------> x  
  Int_t sector;      
  if     (padx >=     1      && padx <=   kPadsSecX )    sector=1;
  else if(padx >  kPadsSecX  && padx <=   kPadsChamX)    sector=2; 
  else                                                   return -1;//padx out of range     
  if     (pady >=     1      && pady <=   kPadsSecY )    return sector;
  else if(pady >   kPadsSecY && pady <= 2*kPadsSecY )    return sector+2;
  else if(pady > 2*kPadsSecY && pady <=   kPadsChamY)    return sector+4;
  else                                                   return -1; //pady out of range
}//Pad2Sec()
//__________________________________________________________________________________________________
#endif
