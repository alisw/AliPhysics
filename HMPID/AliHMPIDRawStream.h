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
#include <TRandom.h>
#include "AliHMPIDParam.h"
#include <AliBitPacking.h>
#include <AliFstream.h>
#include "AliHMPIDDigit.h"
#include "AliDAQ.h"
class AliRawReader;

class AliHMPIDRawStream: public TObject {
  public :
    AliHMPIDRawStream(AliRawReader* rawReader);
    AliHMPIDRawStream();
    
    virtual ~AliHMPIDRawStream();

    virtual void     Reset();
    virtual Bool_t   Next();
            void     Init();
    
	    Int_t    Ch( Int_t ddl,Int_t row,Int_t dil,Int_t pad    ) {return AliHMPIDParam::A2C(fPad[ddl][row][dil][pad]); }            //chamber number
	    
	    Int_t GetDDLNumber()  const { return fDDLNumber; }                                                                            // Provide current DDL number
    inline Int_t GetCharge(Int_t ddl,Int_t row, Int_t dilogic, Int_t pad);                                                                         // Provide the charge observed in certain row,dilogic,pad channel
	    inline Int_t GetPad(Int_t ddl,Int_t row,Int_t dil,Int_t pad);                                                                        //
	    
	    Int_t   Pc          ( Int_t ddl,Int_t row,Int_t dil,Int_t pad                            ) {return AliHMPIDParam::A2P(fPad[ddl][row][dil][pad]);}                                                 //PC position number
	    Int_t   PadPcX      ( Int_t ddl,Int_t row,Int_t dil,Int_t pad                            ) {return AliHMPIDParam::A2X(fPad[ddl][row][dil][pad]);}                                                 //pad pc x # 0..79
	    Int_t   PadPcY      ( Int_t ddl,Int_t row,Int_t dil,Int_t pad                            ) {return AliHMPIDParam::A2Y(fPad[ddl][row][dil][pad]);}                                                 //pad pc y # 0..47
   
	    inline  Bool_t SetZeroSup (Bool_t isSup);
    inline  Bool_t GetZeroSup(); 
                 
    inline void    Raw            (UInt_t &w32,Int_t &ddl,Int_t &r,Int_t &d,Int_t &a);                                               //digit->(w32,ddl,r,d,a)
    inline void    Raw            (Int_t ddl,Int_t r,Int_t d,Int_t a);                                                               //raw->abs pad number
    inline Bool_t  Raw            (UInt_t  w32,Int_t  ddl,AliRawReader *pRR);                                                      //(w32,ddl)->digit
    inline void   SetCharge      (Int_t ddl,Int_t row,Int_t dil,Int_t pad,Int_t q);
    inline void    WriteRaw       (TObjArray *pDigLst                             );                                                      //write as raw stream     
    inline void   WriteRowMarker  (AliFstream *ddl,UInt_t size);
    inline void   WriteEoE        (AliFstream *ddl,UInt_t row,UInt_t dil,UInt_t wordCnt);  
    inline void   WriteSegMarker  (AliFstream *ddl,UInt_t row);   
    
//    inline TClonesArray  ReMap(TClonesArray *pDigIn);
    
    enum EHMPIDRawStreamError {
      kRawDataSizeErr = 1,
      kRowMarkerErr = 2,
      kWrongRowErr = 3,
      kWrongDilogicErr = 4,
      kWrongPadErr = 5,
      kEoEFlagErr = 6,
      kEoESizeErr = 7,
      kEoEDILOGICErr = 8,
      kEoERowErr = 9,
      kBadSegWordErr = 10,
      kWrongSegErr = 11
    };
    
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

    UInt_t           GetNextWord();

    Int_t            fCharge[kNDDL][kNRows+1][kNDILOGICAdd+1][kNPadAdd]; // Array for charge values for all channels in one DDL
    
    Int_t            fPad[kNDDL][kNRows+1][kNDILOGICAdd+1][kNPadAdd]; // Array for abs pad values for all channels in one DDL

    UInt_t           fRawWord[kNDDL][kNRows+1][kNDILOGICAdd+1][kNPadAdd];
        
    Int_t            fDDLNumber;    // index of current DDL number

    AliRawReader*    fRawReader;    // object for reading the raw data

    UChar_t*         fData;         // raw data

    Int_t            fPosition;     // current position in fData
    
    Bool_t           fZeroSup;

    ClassDef(AliHMPIDRawStream, 0)  // base class for reading HMPID raw digits
};
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
  AliBitPacking::PackWord((fCharge[ddl][r][d][a]>4095)?4095:(UInt_t)fCharge[ddl][r][d][a],w32, 0,11);       // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq        Qdc               bits (00..11) counts (0..4095)
  //molnarl: Since in simulation the the charge can be > than 4095 but not in real life we need to protect. If fQ>4095 after packing we will get 0 for the charge! 
  assert(0<=a&&a<=47);AliBitPacking::PackWord(        a ,w32,12,17);  // 3322 2222 2222 1111 1111 1000 0000 0000        DILOGIC address   bits (12..17) counts (0..47)
  assert(1<=d&&d<=10);AliBitPacking::PackWord(        d ,w32,18,21);  // 1098 7654 3210 9876 5432 1098 7654 3210        DILOGIC number    bits (18..21) counts (1..10)
  assert(1<=r&&r<=24);AliBitPacking::PackWord(        r ,w32,22,26);  //                                                Row number        bits (22..26) counts (1..24)  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::Raw(Int_t ddl,Int_t r,Int_t d,Int_t a)
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
  fPad[ddl][r][d][a]=AliHMPIDParam::Abs(ch,pc,px,py);
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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::SetCharge(Int_t ddl,Int_t row,Int_t dil,Int_t pad,Int_t q)
{ 
  //Setter for the charger in the raw stream
  //Arguments: DDL, row number, dilogic number, dilogic address(pad), charge
  //Returns:   Charge from the raw stream
     fCharge[ddl][row][dil][pad]=q;
     //  return fCharge[ddl][row][dil][pad];
      
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDRawStream::GetPad(Int_t ddl,Int_t row,Int_t dil,Int_t pad)
{
  // The method returns the absolute pad number or -1 
  // in case the charge from the channels
  // has not been read or invalid arguments
 
  assert(0<=ddl&&ddl<=13);  
  assert(1<=row&&row<=24); 
  assert(1<=dil&&dil<=10);   
  assert(0<=pad&&pad<=47);  
  
  Int_t a2y[6]={3,2,4,1,5,0};     //pady for a given padress (for single DILOGIC chip)
                                  Int_t ch=ddl/2;
  Int_t tmp=(row-1)/8;            Int_t pc=(ddl%2)? 5-2*tmp:2*tmp; 
                                  Int_t px=(dil-1)*8+pad/6;
        tmp=(ddl%2)?(24-row):row-1;   
                                  Int_t py=6*(tmp%8)+a2y[pad%6];
                                  
  return fPad[ddl][row][dil][pad]=AliHMPIDParam::Abs(ch,pc,px,py);
}//GetPad()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDRawStream::GetCharge(Int_t ddl,Int_t row, Int_t dilogic, Int_t pad)
{
  // The method returns the charge collected
  // in a particular channel
  // Return -1 in case the charge from the channels
  // has not been read or invalid arguments
  if (ddl < 0 || ddl > kNDDL) {
    AliError(Form("Wrong DDL index %d!",ddl));
    return 0;
  }  
  if (row < 1 || row > kNRows) {
    AliError(Form("Wrong row index %d!",row));
    return 0;
  }

  if (dilogic < 1 || dilogic > kNDILOGICAdd) {
    AliError(Form("Wrong DILOGIC address %d!",dilogic));
    return 0;
  }
  
  if (pad >= kNPadAdd) {
    AliError(Form("Wrong pad index %d!",pad));
    return 0;
  }

  return fCharge[ddl][row][dilogic][pad];
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDRawStream::WriteRowMarker(AliFstream *ddl,UInt_t size)
{
  //Writes the row marker for real data and pedestal into the ddl stream
  //Arguments: ddl stream and the size of the block of the given row, the siye is at least the 10 EoE words!
  //Returns:   nothing
  UInt_t w32=0;
  UInt_t marker=12968;                                    //ror marker: 32a8 in hexa; 12968 in decimal
  AliBitPacking::PackWord(size,  w32, 16,31);            //number of roaw written after row marker (digits and EoE)
  AliBitPacking::PackWord(marker,w32,0,15);              //32a8=12968
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
void AliHMPIDRawStream::WriteSegMarker(AliFstream *ddl,UInt_t row)
{
  //Writes the segment marker (after 8 rows) into the ddl stream
  //Arguments: ddl stream and the segment: row 8 -> 0x5800, row 16 -> 5801, row 24 -> 5802
  //Retruns:   nothing
    UInt_t w32=0;
                                AliBitPacking::PackWord((UInt_t)43791,        w32,16,31);       //43791==AB0F
                                AliBitPacking::PackWord((UInt_t)(22528+row/8),w32, 0,15);       //22528==5800
    
    ddl->WriteBuffer((char*)&w32,sizeof(w32)); 
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
Bool_t AliHMPIDRawStream::GetZeroSup()
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
    
    AliRawDataHeader header; header.SetAttribute(0);  //empty DDL header
    
    ddlL = new AliFstream(AliDAQ::DdlFileName("HMPID",2*iCh+1)); //left and right looking at the IP
    ddlR = new AliFstream(AliDAQ::DdlFileName("HMPID",2*iCh));   //open both DDL of this chamber in parallel
    
    ddlL->WriteBuffer((char*)&header,sizeof(header));            //write dummy header as place holder, actual 
    ddlR->WriteBuffer((char*)&header,sizeof(header));            //will be rewritten later when total size of DDL is known
    
    UInt_t w32=0;                 //32 bits data word 
    digcnt=0;
    
    TClonesArray *pDigCh=(TClonesArray *)pDigAll->At(iCh); //list of digits for current chamber 
    //Printf("::::::::::::::::::; Number of Digits to write: %d == %d",pDigCh->GetEntriesFast(),pDigCh->GetEntries());
    for(Int_t iDig=0;iDig<pDigCh->GetEntriesFast();iDig++){//digits loop
      AliHMPIDDigit *pDig1=(AliHMPIDDigit*)pDigCh->At(iDig);
      pDig1->Raw(w32,ddl,r,d,a);
      isDigThere[ddl][r][d][a]=iDig;
    }  
    
    for(Int_t row = 1; row <= AliHMPIDRawStream::kNRows; row++){ //AliHMPIDRawStream::kNRows=25!
      cntRrow=0;cntLrow=0;cntLseg=0;cntRseg=0;// 
      cntLeoe=0;cntReoe=0;
      posLmarker=ddlL->Tellp(); WriteRowMarker(ddlL,(UInt_t)1);   cntL++; cntRrow++; 
      posRmarker=ddlR->Tellp(); WriteRowMarker(ddlR,(UInt_t)1);   cntR++; cntLrow++; 
      for(Int_t dil = 1; dil <= AliHMPIDRawStream::kNDILOGICAdd; dil++){ //AliHMPIDRawStream::kNDILOGICAdd = 11!
	cntLpad=0;cntRpad=0;
        for(Int_t pad = 0; pad < AliHMPIDRawStream::kNPadAdd; pad++){   //AliHMPIDRawStream::kNPadAdd     = 48
	  for ( Int_t iddl=2*iCh; iddl<=2*iCh+1;iddl++){
	    if (isDigThere[iddl][row][dil][pad]!=-1) {
	      AliHMPIDDigit *pDig=(AliHMPIDDigit*)pDigCh->At(isDigThere[iddl][row][dil][pad]);             
	      pDig->Raw(w32,ddl,r,d,a);  
	      if(pDig->Q() < 0 ) continue;                                                 //We can turn of the zero sup for pedestal simulation
              //Printf("::::::::::::::: ddl from Digit : %d",ddl);
	      if(ddl%2){                                                                               //write raw digit selecting on DDL
		ddlL->WriteBuffer((char*)&w32,sizeof(w32));   cntL++; cntLpad++; cntLrow++;  cntLdig++;//Printf(" WL: %x isDig: %d",w32,isDigThere[iddl][row][dil][pad]);
              }else{
		ddlR->WriteBuffer((char*)&w32,sizeof(w32));   cntR++; cntRpad++; cntRrow++;   cntRdig++;//Printf(" WR: %x isDig: %d",w32,isDigThere[iddl][row][dil][pad]);
	      }
            }//ddl 
          }//isDig
	}//pad
        WriteEoE(ddlL,row,dil,cntLpad); cntL++;  cntLrow++;    cntLeoe++;                                 //molnarl: write EoE markers
        WriteEoE(ddlR,row,dil,cntRpad); cntR++;  cntRrow++;    cntReoe++;
      }//dil
      if(row%8==0){                                               
        WriteSegMarker(ddlL,row); cntL++;  cntLseg++;
        WriteSegMarker(ddlR,row); cntR++;  cntRseg++;
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

#endif
