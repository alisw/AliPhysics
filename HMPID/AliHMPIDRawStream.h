#ifndef ALIHMPIDRAWSTREAM_H
#define ALIHMPIDRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading raw data digits for HMPID.
/// The data format is taken from the document provided by Paolo Martinengo.
///
/// cvetan.cheshkov@cern.ch 19/07/2007
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include "AliHMPIDParam.h"
#include <AliBitPacking.h>
#include <AliFstream.h>
#include "AliHMPIDDigit.h"
#include "AliDAQ.h"
#include "AliRawDataHeaderSim.h"

class AliRawReader;

class AliHMPIDRawStream: public TObject {
  public :
    AliHMPIDRawStream(AliRawReader* rawReader);
    AliHMPIDRawStream();
    
    virtual ~AliHMPIDRawStream();

    virtual void     Reset();
    virtual Bool_t   Next();
            void     InitVars(Int_t n);
            void     DelVars();
    
   static  inline Int_t GetPad(Int_t ddl,Int_t row,Int_t dil,Int_t pad);               //get absolute pad number
   static  inline Int_t GetFee(Int_t ddl,Int_t row);                                   //get the FEE number
   static  Int_t GetNDDL()     { return kNDDL;}                                        //return the number of max # of DDLs
   static  Int_t GetNErrors()  { return kSumErr;}                                      //return the number of max # of Error Types

  	    Int_t   GetNPads()         const { return fNPads;}                         //Get number of pads present in the stream
            Int_t*  GetPadArray()      const { return fPad;}                           //Get pad array from stream decoded
            Int_t*  GetChargeArray()   const { return fCharge;}                        //Get the charge of the pads from dedcoded stream 
            Int_t*  GetnDDLInStream()  const { return fnDDLInStream;}                  //Get the DDL input check array
            Int_t*  GetnDDLOutStream() const { return fnDDLOutStream;}                 //Get the DDL output check array
	    Int_t   Pc          ( Int_t ddl,Int_t row,Int_t dil,Int_t pad                            ) {return AliHMPIDParam::A2P(GetPad(ddl,row,dil,pad));}                                                 //PC position number
	    Int_t   PadPcX      ( Int_t ddl,Int_t row,Int_t dil,Int_t pad                            ) {return AliHMPIDParam::A2X(GetPad(ddl,row,dil,pad));}                                                 //pad pc x # 0..79
	    Int_t   PadPcY      ( Int_t ddl,Int_t row,Int_t dil,Int_t pad                            ) {return AliHMPIDParam::A2Y(GetPad(ddl,row,dil,pad));}                                                 //pad pc y # 0..47

    static inline const Char_t* GetErrName(Int_t eType);      
    inline  Bool_t SetZeroSup (Bool_t isSup);
    inline  Bool_t GetZeroSup()const; 
    inline  Int_t  GetErrors(Int_t ddl,Int_t eType)const;                                                                                   //Get errors and occurance
    Int_t  GetDDLNumber() const{ return fDDLNumber;} //return the number of DDL actually being decoded
    UInt_t GetLDCNumber() const{ return fLDCNumber;} //return the number of LDC actually being decoded
    UInt_t GetTimeStamp() const{ return fTimeStamp;} //return the time stamp of the event actually being decoded

    void   SetTurbo(Bool_t isTurbo){fTurbo=isTurbo;}     // Enable/Disable Turbo
    Bool_t GetTurbo(){ return fTurbo;}                   // Enable Turbo
    Bool_t Turbo();                                      // Read HMPID Raw data without error checks
    Bool_t ReadHMPIDRawData();                           // Read HMPID Raw data
    Bool_t ReadSegment(Int_t &cntSegment);               // Read Segment
    Bool_t ReadRow(Int_t &cntRow);                       // Read Row
    Bool_t ReadDilogic(Int_t &cntDilogic);               // Read Dilogic

    Bool_t CheckRow(UInt_t row);                         // Check Row
    Bool_t CheckDilogic(UInt_t dilogic);                 // Check Dilogic
    Bool_t CheckPad(UInt_t pad);                         // Check pad
    Bool_t CheckEoE(Int_t &nDil);                        // Check EoE
    Bool_t CheckRowMarker();                             // Check RowMarker
    Bool_t CheckSegment();                               // Check Segment
    void   DumpData(Int_t nw);                           // Dump Data
    void   StorePosition();                              //Debug purpose
    
    Double_t GetDdlDataSize()      { return 4.0*fRawDataSize;} //returns the data size for the DDL which is decoded in Next(); fRawDataSize = Bytes/4  
    
//    inline void    Raw            (UInt_t &w32,Int_t &ddl,Int_t &r,Int_t &d,Int_t &a);                                              //digit->(w32,ddl,r,d,a)
//    inline void    Raw            (Int_t ddl,Int_t r,Int_t d,Int_t a);                                                              //raw->abs pad number
//    inline Bool_t  Raw            (UInt_t  w32,Int_t  ddl,AliRawReader *pRR);                                                       //(w32,ddl)->digit

    inline void   WriteRaw       (TObjArray *pDigLst                             );                 //write as raw stream     
    inline void   WriteRowMarker  (AliFstream *ddl,UInt_t size);                                    //write row marker in simulation
    inline void   WriteEoE        (AliFstream *ddl,UInt_t row,UInt_t dil,UInt_t wordCnt);           //write Enf Of Event word in simulation
    inline void   WriteSegMarker  (AliFstream *ddl,UInt_t row, Int_t nwInSeg);                      //write Segment Marker word in simulation
    inline void   Write5FirmwareWords(AliFstream *ddl);                                             //write the firmware control words in simulation
    
//    inline TClonesArray  ReMap(TClonesArray *pDigIn);
enum EDirection {kFwd,kBwd};

enum Ebits {kbit0,kbit1 , kbit2, kbit3, kbit4, kbit5, kbit6, kbit7, kbit8,
                  kbit9 ,kbit10,kbit11,kbit12,kbit13,kbit14,kbit15,kbit16,
                  kbit17,kbit18,kbit19,kbit20,kbit21,kbit22,kbit23,kbit24,
                  kbit25,kbit26,kbit27,kbit28,kbit29,kbit30,kbit31,kbit32};
    
  enum EHMPIDRawStreamError { kRawDataSizeErr   = 0,  kRowMarkerErr     = 1,  kWrongRowErr      = 2,  kWrongDilogicErr  = 3,
                              kWrongPadErr      = 4,  kEoEFlagErr       = 5,  kEoESizeErr       = 6,  kEoEDILOGICErr    = 7,
                              kEoERowErr        = 8,  kBadSegWordErr    = 9,  kWrongSegErr      = 10, kRowMarkerSizeErr = 11,
                              kPedQZero         =12,  kSumErr           = 13  //This is always the last one, to retreive the number of errors
                            };                        //Always check the updated list of names in the .cxx file for print-out!
   
  enum {
      kNRows       = 24,                                    // Number of rows (starting from 1 !)//was25
      kNDILOGICAdd = 10,                                    // Number of DILOGIC addresses in a row (starting from 1 !) //was11
      kNPadAdd     = 48,                                    // Number of pad row
      kNRowsPerSegment = 8,                                 // Number of rows per segment
      kNDDL = 14
    };
 enum EHMPIDRawError {
    kInvalidRawDataWord = 1
  };

    
  private :

    AliHMPIDRawStream& operator = (const AliHMPIDRawStream& stream);
    AliHMPIDRawStream(const AliHMPIDRawStream& stream);

    Bool_t           GetWord(Int_t n=1,EDirection dir=kFwd);             // Get n-th word
    UInt_t           GetNextWord();                                      // Get next word
    Int_t            fNPads;                                             // counter of pads in one DDL
    Int_t           *fCharge;                                            // Array for charge values for all channels in one DDL
    Int_t           *fPad;                                               // Array for abs pad values for all channels in one DDL
    Int_t            fDDLNumber;                                         // index of current DDL number
    Int_t           *fnDDLInStream;                                      // if the DDL is in the raw data
    Int_t           *fnDDLOutStream;                                     // if the DDL is in the raw data
    UInt_t           fLDCNumber;                                         // index of current LDC number
    UInt_t           fTimeStamp;                                         // TimeStamp
    AliRawReader    *fRawReader;                                         // object for reading the raw data
    UChar_t         *fData;                                              // raw data
    Int_t          **fNumOfErr;                                          // Store the numner of errors for a given error type and a given DDL
    Int_t            fPosition;                                          // current word
    UInt_t           fWord;                                              // current position in fData
    Bool_t           fZeroSup;                                           // set if zero suppression is applied
    Int_t           *fPos;                                               // for debug purposes
    Int_t            fiPos;                                              // counter for debug
    Bool_t           fTurbo;                                             // kTRUE = Turbo decoding is called. DEFAULT: kFALSE = normal decoding is called
    Int_t            fRawDataSize;
    ClassDef(AliHMPIDRawStream, 4)                                       // base class for reading HMPID raw digits
};
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /*
void AliHMPIDRawStream::Raw(UInt_t &w32,Int_t &ddl,Int_t &r,Int_t &d,Int_t &a)
{
// Convert raw stream word to raw word format
// Arguments: w32,ddl,r,d,a where to write the results
// Returns: none
  Int_t y2a[6]={5,3,1,0,2,4};
  
  ddl=2*Ch(ddl,r,d,a)+Pc(ddl,r,d,a)%2;                                                          //DDL# 0..13
  Int_t tmp=1+Pc(ddl,r,d,a)/2*8+PadPcY(ddl,r,d,a)/6;  r=(Pc(ddl,r,d,a)%2)? 25-tmp:tmp;              //row r=1..24
  d=1+PadPcX(ddl,r,d,a)/8;                                                                  //DILOGIC# 1..10
  a=y2a[PadPcY(ddl,r,d,a)%6]+6*(PadPcX(ddl,r,d,a)%8);                                           //ADDRESS 0..47        
  
  w32=0;    
  AliBitPacking::PackWord((fCharge[fNPads]>4095)?4095:(UInt_t)fCharge[fNPads],w32, 0,11);       // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq        Qdc               bits (00..11) counts (0..4095)
  //molnarl: Since in simulation the the charge can be > than 4095 but not in real life we need to protect. If fQ>4095 after packing we will get 0 for the charge! 
  assert(0<=a&&a<=47);AliBitPacking::PackWord(        a ,w32,12,17);  // 3322 2222 2222 1111 1111 1000 0000 0000        DILOGIC address   bits (12..17) counts (0..47)
  assert(1<=d&&d<=10);AliBitPacking::PackWord(        d ,w32,18,21);  // 1098 7654 3210 9876 5432 1098 7654 3210        DILOGIC number    bits (18..21) counts (1..10)
  assert(1<=r&&r<=24);AliBitPacking::PackWord(        r ,w32,22,26);  //                                                Row number        bits (22..26) counts (1..24)  
}
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /*
Int_t AliHMPIDRawStream::Raw(Int_t ddl,Int_t r,Int_t d,Int_t a)
{
  //Assign absolute pad ID based on ddl,row,dil,pad
  //Arguments: DDL, row number, dilogic number, dilogic address(pad)
  //Returns  : nothing

  assert(0<=ddl&&ddl<=13); assert(1<=r&&r<=24); assert(1<=d&&d<=10);   assert(0<=a&&a<=47);  
  Int_t a2y[6]={3,2,4,1,5,0};//pady for a given address (for single DILOGIC chip)
                                  Int_t ch=ddl/2;
  Int_t tmp=(r-1)/8;              Int_t pc=(ddl%2)? 5-2*tmp:2*tmp; 
                                  Int_t px=(d-1)*8+a/6;
        tmp=(ddl%2)?(24-r):r-1;   Int_t py=6*(tmp%8)+a2y[a%6];
  return AliHMPIDParam::Abs(ch,pc,px,py);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDRawStream::Raw(UInt_t w32,Int_t ddl, AliRawReader *pRR)
{
// Converts a given raw data word to a digit
// Arguments: w32 - 32 bits raw data word
//            ddl - DDL idx  0 1 2 3 4 ... 13
//   Returns: none
  Int_t r = AliBitPacking::UnpackWord(w32,22,26); assert(1<=r&&r<=24);   //                                         Row number      (1..24)    
  Int_t d = AliBitPacking::UnpackWord(w32,18,21); assert(1<=d&&d<=10);   // 3322 2222 2222 1111 1111 1000 0000 0000 DILOGIC number  (1..10)
  Int_t a = AliBitPacking::UnpackWord(w32,12,17); assert(0<=a&&a<=47);   // 1098 7654 3210 9876 5432 1098 7654 3210 DILOGIC address (0..47)  
  Int_t q = AliBitPacking::UnpackWord(w32, 0,11); assert(0<=q&&q<=4095); // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq Qdc             (0..4095) 
  if (r<1 || r>24 || d<1 || d>10 || a<0 || a>47 || q<0 || q>4095) {
    AliWarning(Form("Invalid raw data word %x",w32));
    pRR->AddMajorErrorLog(kInvalidRawDataWord,Form("w=%x",w32));
    return kFALSE;
  }
  Raw(ddl,r,d,a);
  fCharge[ddl][r][d][a]=q;
  return kTRUE;
}
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDRawStream::GetPad(Int_t ddl,Int_t row,Int_t dil,Int_t pad)
{
  // The method returns the absolute pad number or -1 
  // in case the charge from the channels
  // has not been read or invalid arguments
 
 
  if(ddl<0 || ddl >13 || row<1 || row >25 || dil<1 || dil >10 || pad<0 || pad >47 ) return -1;
  Int_t a2y[6]={3,2,4,1,5,0};     //pady for a given padress (for single DILOGIC chip)
  Int_t ch=ddl/2;
  Int_t tmp=(24-row)/8;
  Int_t pc=(ddl%2)?5-2*tmp:2*tmp;
  Int_t px=(kNDILOGICAdd+1 - dil)*8-pad/6-1;  //flip according to Paolo (2-9-2008)

  tmp=(ddl%2)?row-1:(24-row);
  Int_t py=6*(tmp%8)+a2y[pad%6];

  return AliHMPIDParam::Abs(ch,pc,px,py);
}//GetPad()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDRawStream::GetFee(Int_t ddl,Int_t row)
{
   // The method returns the FEE number or -1 
   // in case of invalid argument(s)
 
   if(ddl<0 || ddl >13 || row<1 || row >25  ) return -1;
   Int_t fee=-1, kLeft=0, kRight=0;
   if(ddl%2==0) {kLeft=1;kRight=0;}
   if(ddl%2!=0) {kLeft=0;kRight=1;}
   if((kLeft==1 && row<=24 && row>=21) || (kRight==1 && row<=4  && row>=1)) fee=0;    //calculation of FEE values based on Giacomos pedestal macro
   if((kLeft==1 && row<=20 && row>=17) || (kRight==1 && row<=8  && row>=5)) fee=1;
   if((kLeft==1 && row<=16 && row>=13) || (kRight==1 && row<=12 && row>=9)) fee=2;
   if((kLeft==1 && row<=12 && row>=9)  || (kRight==1 && row<=16 && row>=13))fee=3;
   if((kLeft==1 && row<=8  && row>=5)  || (kRight==1 && row<=20 && row>=17))fee=4;
   if((kLeft==1 && row<=4  && row>=1)  || (kRight==1 && row<=24 && row>=21))fee=5;
   return fee;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::WriteRowMarker(AliFstream *ddl,UInt_t size)
{
  //Writes the row marker for real data and pedestal into the ddl stream
  //Arguments: ddl stream and the size of the block of the given row, the siye is at least the 10 EoE words!
  //Returns:   nothing
  UInt_t w32=0;
  UInt_t marker=13992;                                   //for pedestal=12968  ==  32a8 for zero suppressed 36a8
  AliBitPacking::PackWord(size,  w32, 16,31);            //number of roaw written after row marker (digits and EoE)
  AliBitPacking::PackWord(marker,w32,0,15);              //the marker word
  ddl->WriteBuffer((char*)&w32,sizeof(w32));              
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::WriteEoE(AliFstream *ddl,UInt_t row,UInt_t dil,UInt_t wordCnt  )
{
  //Writes the EoE word from real data and pedestals into the ddl stream
  //Arguments:  ddl stream, row number, dilogic number and the number of words before the EoE
  //Retursns:   nothing
  UInt_t e=1;
  UInt_t w32=0;
  assert(1<=row&&row<=24);      AliBitPacking::PackWord((UInt_t)row     ,w32,22,26);    // row number (1...24)
  assert(1<=dil&&dil<=10);      AliBitPacking::PackWord((UInt_t)dil     ,w32,18,21);    // DILOGIC number (1...10)
	                        AliBitPacking::PackWord(          e     ,w32, 7,17);   // event number -- not used
                        	AliBitPacking::PackWord((UInt_t)wordCnt ,w32, 0, 6);  // word counter (0...47)                                                           	AliBitPacking::PackWord((UInt_t)1       ,w32,27,27);  // bit 27 is always 1 by definition of EoE
                                AliBitPacking::PackWord((UInt_t)1       ,w32,27,27);  // bit 27 is always 1 by definition of EoE    
  ddl->WriteBuffer((char*)&w32,sizeof(w32));      
} 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
void AliHMPIDRawStream::WriteSegMarker(AliFstream *ddl,UInt_t row, Int_t nwInSeg)
{
  //Writes the segment marker (after 8 rows) into the ddl stream
  //Arguments: ddl stream and the segment: row 8 -> 0x5800, row 16 -> 5801, row 24 -> 5802 for pedestal
  //Retruns:   nothing
    UInt_t w32=0;

      //Segment marker: 2736 == ab0
      //AliBitPacking::PackWord((UInt_t)0   ,w32,27,31);          //zero out the rest of the bits, since they are not needed
      AliBitPacking::PackWord((UInt_t)2736   ,w32,20,31);       //ab0 the segment marker word
      AliBitPacking::PackWord((UInt_t)nwInSeg,w32, 8,19);       //number of words in the segment
      AliBitPacking::PackWord((UInt_t)(row/8),w32, 0, 7);       //segment 0,1,2    
      ddl->WriteBuffer((char*)&w32,sizeof(w32)); 
      //Printf("Segment word created is: %x",w32);
}      
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
void AliHMPIDRawStream::Write5FirmwareWords(AliFstream *ddl)
{
  //Before each DDL payload 5 words are written: 
  // 1.) Firmware version,              for sim = 999
  // 2.) Status and error bits from CD, for sim = 0
  // 3.) # FEE RESET received         , for sim = 0
  // 4.) # TTC READY                  , for sim = 0  
  // 5.) Spare/Reserved               , for sim = 0
  //Returns:   nothing
  UInt_t w32=0;
  AliBitPacking::PackWord((UInt_t)999,w32,0,31); ddl->WriteBuffer((char*)&w32,sizeof(w32));              
  AliBitPacking::PackWord((UInt_t) 10,w32,0,31); ddl->WriteBuffer((char*)&w32,sizeof(w32));              
  AliBitPacking::PackWord((UInt_t) 11,w32,0,31); ddl->WriteBuffer((char*)&w32,sizeof(w32));              
  AliBitPacking::PackWord((UInt_t) 12,w32,0,31); ddl->WriteBuffer((char*)&w32,sizeof(w32));              
  AliBitPacking::PackWord((UInt_t) 13,w32,0,31); ddl->WriteBuffer((char*)&w32,sizeof(w32));              
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
Bool_t AliHMPIDRawStream::SetZeroSup (Bool_t isSup)
{
  //Prevision to turn OFF zero suppression
  //Arguments: setter
  //Returns:   switch
  fZeroSup=isSup;
  return fZeroSup;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
Bool_t AliHMPIDRawStream::GetZeroSup()const
{
  if(fZeroSup==kTRUE) return kTRUE;
  else                return kFALSE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
void AliHMPIDRawStream::WriteRaw(TObjArray *pDigAll)
{
// Write a list of digits for a given chamber in raw data stream
// Arguments: pDigAll- list of digits 
//   Returns: none      
  Int_t  ddl,r,d,a;            //32 bits data word 
  Int_t  cntLpad,cntRpad;
  Int_t  cntLrow,cntRrow;
  Int_t  cntL=0,cntR=0;                           //data words counters for DDLs
  Int_t  cntLeoe,cntReoe;
  UInt_t posL,posR;
  UInt_t cntLseg,cntRseg;
  UInt_t cntwInLseg=0,cntwInRseg=0;
  Int_t  cntRdig=0,cntLdig=0;
  
  UInt_t posLmarker,posRmarker;
  Int_t digcnt=0;

  Int_t isDigThere[14][25][11][48];
  
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){//chambers loop
    cntL=0;cntR=0;   
    for(Int_t iddl=0;iddl<14;iddl++){
      for(Int_t irow=1;irow<=24;irow++){
	for(Int_t idil=1;idil<=10;idil++){
	  for(Int_t ipad=0;ipad<48;ipad++){
	    isDigThere[iddl][irow][idil][ipad]=-1;
	  }
	}
      }
    }
    
    AliFstream* ddlL;                                 //output streams, 2 per chamber
    AliFstream* ddlR;                          
    
    AliRawDataHeaderSim header; header.SetAttribute(0);  //empty DDL header
    
    ddlL = new AliFstream(AliDAQ::DdlFileName("HMPID",2*iCh+1)); //left and right looking at the IP
    ddlR = new AliFstream(AliDAQ::DdlFileName("HMPID",2*iCh));   //open both DDL of this chamber in parallel
    
    ddlL->WriteBuffer((char*)&header,sizeof(header));            //write dummy header as place holder, actual 
    ddlR->WriteBuffer((char*)&header,sizeof(header));            //will be rewritten later when total size of DDL is known
    
    UInt_t w32=0;                 //32 bits data word 
    digcnt=0;
    
    //added frimware control words
    Write5FirmwareWords(ddlL);  cntL+=5;
    Write5FirmwareWords(ddlR);  cntR+=5;
   
    
    TClonesArray *pDigCh=(TClonesArray *)pDigAll->At(iCh); //list of digits for current chamber 
   
    for(Int_t iDig=0;iDig<pDigCh->GetEntriesFast();iDig++){//digits loop
      AliHMPIDDigit *pDig1=(AliHMPIDDigit*)pDigCh->At(iDig);
      pDig1->Raw(w32,ddl,r,d,a);  //??????????
      isDigThere[ddl][r][d][a]=iDig;
    }  
    
    for(Int_t row = 1; row <= AliHMPIDRawStream::kNRows; row++){ //AliHMPIDRawStream::kNRows=25!
      cntRrow=0;cntLrow=0;cntLseg=0;cntRseg=0;// 
      cntLeoe=0;cntReoe=0;
      posLmarker=ddlL->Tellp(); WriteRowMarker(ddlL,(UInt_t)1);   cntL++; cntRrow++; cntwInRseg++;
      posRmarker=ddlR->Tellp(); WriteRowMarker(ddlR,(UInt_t)1);   cntR++; cntLrow++; cntwInLseg++;
      for(Int_t dil = 1; dil <= AliHMPIDRawStream::kNDILOGICAdd; dil++){ //AliHMPIDRawStream::kNDILOGICAdd = 11!
	cntLpad=0;cntRpad=0;
        for(Int_t pad = 0; pad < AliHMPIDRawStream::kNPadAdd; pad++){   //AliHMPIDRawStream::kNPadAdd     = 48
	  for ( Int_t iddl=2*iCh; iddl<=2*iCh+1;iddl++){
	    if (isDigThere[iddl][row][dil][pad]!=-1) {
	      AliHMPIDDigit *pDig=(AliHMPIDDigit*)pDigCh->At(isDigThere[iddl][row][dil][pad]);             
	      pDig->Raw(w32,ddl,r,d,a);  
	      if(pDig->Q() < 0 ) continue;                                                 //We can turn of the zero sup for pedestal simulation
                if(ddl%2){                                                                               //write raw digit selecting on DDL
		ddlL->WriteBuffer((char*)&w32,sizeof(w32));   cntL++; cntLpad++; cntLrow++;  cntLdig++; cntwInLseg++;//Printf(" WL: %x isDig: %d",w32,isDigThere[iddl][row][dil][pad]);
              }else{
		ddlR->WriteBuffer((char*)&w32,sizeof(w32));   cntR++; cntRpad++; cntRrow++;   cntRdig++;cntwInRseg++;//Printf(" WR: %x isDig: %d",w32,isDigThere[iddl][row][dil][pad]);
	      }
            }//ddl 
          }//isDig
	}//pad
        WriteEoE(ddlL,row,dil,cntLpad); cntL++;  cntLrow++;    cntLeoe++;   cntwInLseg++;                              //molnarl: write EoE markers
        WriteEoE(ddlR,row,dil,cntRpad); cntR++;  cntRrow++;    cntReoe++;   cntwInRseg++;
      }//dil
      if(row%8==0){                                               
        WriteSegMarker(ddlL,row,cntwInLseg); cntL++;  cntLseg++; cntwInLseg=0;
        WriteSegMarker(ddlR,row,cntwInRseg); cntR++;  cntRseg++;  cntwInRseg=0; 
      }
      posL=ddlL->Tellp();   ddlL->Seekp(posLmarker);    WriteRowMarker(ddlL,(UInt_t)(cntLrow-1)); ddlL->Seekp(posL);      //find the marker position write and  go back to the actual position to continue writing                    
      posR=ddlR->Tellp();   ddlR->Seekp(posRmarker);    WriteRowMarker(ddlR,(UInt_t)(cntRrow-1)); ddlR->Seekp(posR);                           
    }//row
    header.fSize=sizeof(header)+cntL*sizeof(w32); ddlL->Seekp(0); ddlL->WriteBuffer((char*)&header,sizeof(header)); delete ddlL; //rewrite header with size set to
    header.fSize=sizeof(header)+cntR*sizeof(w32); ddlR->Seekp(0); ddlR->WriteBuffer((char*)&header,sizeof(header)); delete ddlR; //number of bytes and close file
    
    //Printf("In Ch %d # digits written to LDD %d RDDL %d",iCh,cntLdig,cntRdig);
    
  }//chambers loop
}//WriteRaw()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDRawStream::GetErrors(Int_t ddl,Int_t eType)const
{
// Return the number of errors for a given error tye during raw data reading
// Arguments: errorType
// Returns: error or -999 if error Type does not exist
  
  if(eType < 0 || eType> kSumErr ||  ddl < 0 || ddl > kNDDL-1 ) return -999;
  else return fNumOfErr[ddl][eType];
} //GetErrors()     
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const Char_t* AliHMPIDRawStream::GetErrName(Int_t eType)
{
  // Return the name of the error for a given error tye during raw data reading
  // Arguments: errorType
  // Returns: error or -999 if error Type does not exist
  const Char_t *eName[]={ "kRawDataSizeErr",  "kRowMarkerErr" , "kWrongRowErr" , "kWrongDilogicErr",
                    "kWrongPadErr"   ,  "kEoEFlagErr"   , "kEoESizeErr"  , "kEoEDILOGICErr",
                    "kEoERowErr"     ,  "kBadSegWordErr", "kWrongSegErr" , "kRowMarkerSizeErr",
                    "kPedQZero"      ,  "kSumErr"        };                       
  const Char_t *eNoErr="NotDefinedErrorType";
  if(eType<0 || eType>kSumErr) return eNoErr;
  else                         return eName[eType];
}//GetErrName()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    
#endif
