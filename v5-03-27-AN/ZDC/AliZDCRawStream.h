#ifndef ALIZDCRAWSTREAM_H
#define ALIZDCRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////
//						//
//  Class to provide access to ZDC raw data	//
//  Author: Chiara Oppedisano			//
//						//
//////////////////////////////////////////////////

#include <TObject.h>
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliZDCChMap.h"

class AliRawReader;
struct AliRawDataHeader;


class AliZDCRawStream: public TObject {
  public :
    
    // Module type codes
    enum ZDCModules{kV965=1, kV830=2, kTRG=3, kTRGI=4, kPU=5, KV1290=6, kV775N=7}; 
    
    // Module type codes
    enum ZDCGeoAddr{kFirstADCGeo=0, kLastADCGeo=3, kADDADCGeo=5,
  	 kTDCFakeGeo=8, kZDCTDCGeo=4, kADDTDCGeo=6,
	 kScalerGeo=16, kPUGeo=29, kTrigScales=30, kTrigHistory=31};
    
    // Signal codes for ZDC 
    // Same codes used in DAQ configuration file
    // To be changed ONLY IF this file is changed!!! 
    // **** DO NOT CHANGE THE FOLLOWING LINES!!! ****
    enum ZDCSignal{
         kNotConnected=0, kVoid=1,
	 kZNAC=2, kZNA1=3, kZNA2=4, kZNA3=5, kZNA4=6,
	 kZPAC=7, kZPA1=8, kZPA2=9, kZPA3=10, kZPA4=11,
	 kZNCC=12, kZNC1=13, kZNC2=14, kZNC3=15, kZNC4=16,
	 kZPCC=17, kZPC1=18, kZPC2=19, kZPC3=20, kZPC4=21,
	 kZEM1=22, kZEM2=23,
	 kZDCAMon=24, kZDCCMon=25,
	 kZNACoot=26, kZNA1oot=27, kZNA2oot=28, kZNA3oot=29, kZNA4oot=30,
	 kZPACoot=31, kZPA1oot=32, kZPA2oot=33, kZPA3oot=34, kZPA4oot=35,
	 kZNCCoot=36, kZNC1oot=37, kZNC2oot=38, kZNC3oot=39, kZNC4oot=40,
	 kZPCCoot=41, kZPC1oot=42, kZPC2oot=43, kZPC3oot=44, kZPC4oot=45,
	 kZEM1oot=46, kZEM2oot=47,
	 kZDCAMonoot=48, kZDCCMonoot=49,
	 kL1MBI=50, kL1CNI=51, kL1SCI=52, kL1EMDI=53, kL0I=54, 
	 kL1MBO=55, kL1CNO=56, kL1SCO=57, kL1EMDO=58, 
	 kHMBCN=59, kHSCEMD=60,
	 kZNACD=61, kZNA1D=62, kZNA2D=63, kZNA3D=64, kZNA4D=65,
	 kZPACD=66, kZPA1D=67, kZPA2D=68, kZPA3D=69, kZPA4D=70,
	 kZNCCD=71, kZNC1D=72, kZNC2D=73, kZNC3D=74, kZNC4D=75,
	 kZPCCD=76, kZPC1D=77, kZPC2D=78, kZPC3D=79, kZPC4D=80,
	 kZEM1D=81, kZEM2D=82,
	 kZDCAMonD=83, kZDCCMonD=84,
	 kZNAD=85, kZPAD=86, kZNCD=87, kZPCD=88, kZEMD=89,
         kZNA0D=90, kZPA0D=91, kZNC0D=92, kZPC0D=93, k1kHzD=94, kGate=95, kAD=96, kCD=97, 
	 kAorCD=98, kAandCD=99, kZEMORD=100, kAorCorZEMORD=101, kAorCorZEMD=102, kAD0=103, kAD1=104, kAD2=105, 
	 kAD3=106, kAD4=107, kAD5=108, kAD6=109, kAD7=110, kAD8=111, kAD9=112, kAD10=113, 
	 kAD11=114, kAD12=115, kAD13=116, kAD14=117, kAD15=118, kAD0D=119, kAD1D=120, kAD2D=121,
	 kAD3D=122, kAD4D=123, kAD5D=124, kAD6D=125, kAD7D=126, kAD8D=127, kAD9D=128, kAD10D=129,
	 kAD11D=130, kAD12D=131, kAD13D=132, kAD14D=133, kAD15D=134
	 };
    
    // Error codes in raw data streaming
    enum EZDCRawStreamError{
       kCDHError = 1,
       kDARCError = 2,
       kZDCDataError = 3,
       kInvalidADCModule = 4,
       kInvalidSector = 5};
	 
    AliZDCRawStream(AliRawReader* rawReader); 
    virtual ~AliZDCRawStream();
    virtual Bool_t   Next();
    
    virtual void ReadChMap();

    virtual void ReadCDHHeader();

    UInt_t GetRawBuffer()      const {return fBuffer;}
    Int_t  GetReadOutCard()    const {return fReadOutCard;}
    
    Int_t  GetDeadfaceOffset() const {return fDeadfaceOffset;}
    Int_t  GetDeadbeefOffset() const {return fDeadbeefOffset;}
    Int_t  GetDataOffset()     const {return fDataOffset;}

    Int_t  GetSector(Int_t i) const {return fSector[i];}
    Int_t  GetModType()       const {return fModType;}
    Int_t  GetADCModule()     const {return fADCModule;}
    Int_t  GetADCNChannels()  const {return fADCNChannels;}
    Int_t  GetADCChannel()    const {return fADCChannel;}
    Int_t  GetADCValue()      const {return fADCValue;}
    Int_t  GetADCGain()       const {return fADCGain;}

    Int_t  GetModuleGEO()     const {return fADCModule;}
    Int_t  GetChannel()       const {return fADCChannel;}
    
    Bool_t IsZDCTDCHeader()   const {return fIsZDCTDCHeader;}
    Bool_t IsZDCTDCDatum()    const {return fIsZDCTDCdatum;}
    Int_t  GetZDCTDCDatum()   const {return fZDCTDCdatum;}
    
    Bool_t IsADDADCChannel()  const {return fIsADDChannel;}
    Int_t  GetADDADCDatum()   const {return fADDADCdatum;}
    Bool_t IsADDTDCdatum()    const {return fIsADDTDCdatum;}
    Int_t  GetADDTDCDatum()   const {return fADDTDCdatum;}
    
    AliCDBStorage *SetStorage(const char* uri);

    // Map from OCDB
    AliZDCChMap   *GetChMap() const;
    //  ADC map
    Int_t  GetNChannelsOn()	       const {return fNChannelsOn;}
    Int_t  GetCabledSignal()           const {return fCabledSignal;}
    Int_t  GetADCModFromMap(Int_t i)   const {return fMapADC[i][0];}
    Int_t  GetADCChFromMap(Int_t i)    const {return fMapADC[i][1];}
    Int_t  GetADCSignFromMap(Int_t i)  const {return fMapADC[i][2];}
    Int_t  GetDetectorFromMap(Int_t i) const {return fMapADC[i][3];}
    Int_t  GetTowerFromMap(Int_t i)    const {return fMapADC[i][4];}
    //  Scaler map
    Int_t  GetScalerModFromMap(Int_t i)  const {return fScalerMap[i][0];}
    Int_t  GetScalerChFromMap(Int_t i)   const {return fScalerMap[i][1];}
    Int_t  GetScalerSignFromMap(Int_t i) const {return fScalerMap[i][2];}
    Int_t  GetScDetectorFromMap(Int_t i) const {return fScalerMap[i][3];}
    Int_t  GetScTowerFromMap(Int_t i)    const {return fScalerMap[i][4];}
    //  TDC map
    Int_t  GetTDCModFromMap(Int_t i)  const {return fTDCMap[i][0];}
    Int_t  GetTDCChFromMap(Int_t i)   const {return fTDCMap[i][1];}
    Int_t  GetTDCSignFromMap(Int_t i) const {return fTDCMap[i][2];}
    
    Bool_t IsCalibration()   const {return fIsCalib;}
    Bool_t IsDARCHeader()    const {return fIsDARCHeader;}
    Bool_t IsHeaderMapping() const {return fIsHeaderMapping;}
    Bool_t IsChMapping()     const {return fIsChMapping;}
    Bool_t IsADCDataWord()   const {return fIsADCDataWord;}
    Bool_t IsADCHeader()     const {return fIsADCHeader;}
    Bool_t IsADCEOB()	     const {return fIsADCEOB;}
    Bool_t IsUnderflow()     const {return fIsUnderflow;}
    Bool_t IsOverflow()      const {return fIsOverflow;}
    
    UInt_t GetScGeo()           const {return fScGeo;}	    
    UInt_t GetScNWords()        const {return fScNWords;}	    
    UInt_t GetScTriggerSource() const {return fScTriggerSource;}	    
    UInt_t GetTriggerNumber()   const {return fScTriggerNumber;}
    UInt_t GetTriggerCount()    const {return fScEvCounter;}
    Bool_t IsScHeaderRead()     const {return fIsScHeaderRead;}
    Bool_t IsScEventGood()      const {return fIsScEventGood;}
    Bool_t IsScalerWord()       const {return fIsScalerWord;}
    
    UInt_t GetDetectorPattern() const {return fDetPattern;}
    
    Bool_t IstriggerHistoryWord() const {return fIsTriggerHistory;}
    Int_t  GetTriggerInput2CTP() const {return *fCPTInput;}
    Bool_t IsCPTInputMBTrigger() 
    	{if(fCPTInput[0]==1) return kTRUE; else return kFALSE;}
    Bool_t IsCPTInputCentralTrigger()
    	{if(fCPTInput[1]==1) return kTRUE; else return kFALSE;}
    Bool_t IsCPTInputSemiCentralTrigger()
    	{if(fCPTInput[2]==1) return kTRUE; else return kFALSE;}
    Bool_t IsCPTInputEMDTrigger()
    	{if(fCPTInput[3]==1) return kTRUE; else return kFALSE;}
    
    Bool_t IsADCEventGood() const {return fIsADCEventGood;} 
    Bool_t IsL0BitSet()     const {return fIsL0BitSet;}  
    Bool_t IsPileUpEvent()  const {return fIsPileUpEvent;} 
    
    void SetReadOutCard(Int_t icard)  {fReadOutCard=icard;}
    void SetDataOffset(Int_t iOffset) {fDataOffset=iOffset;}
    void SetNChannelsOn(Int_t val) {fNChannelsOn = val;}
    void SetSector(Int_t i, Int_t val) {fSector[i] = val;}
    void SetMapRead(Bool_t value) {fIsMapRead=value;}
    void SetMapADCMod(Int_t iraw, Int_t imod) {fMapADC[iraw][0]=imod;}
    void SetMapADCCh(Int_t iraw, Int_t ich)   {fMapADC[iraw][1]=ich;}
    void SetMapADCSig(Int_t iraw, Int_t isig) {fMapADC[iraw][2]=isig;}
    void SetMapDet(Int_t iraw, Int_t idet)    {fMapADC[iraw][3]=idet;}
    void SetMapTow(Int_t iraw, Int_t itow)    {fMapADC[iraw][4]=itow;}
    void SetReadCDH(Bool_t value) {fReadCDH=value;}
    void SetSODReading(Bool_t iset) {fSODReading = iset;}
     
  private :
    AliZDCRawStream(const AliZDCRawStream& stream);
    AliZDCRawStream& operator = (const AliZDCRawStream& stream);

    AliRawReader* fRawReader;    // object for reading the raw data
    
    // Data for buffer decoding
    UInt_t fBuffer;	      // [DARC header +] ADC buffer
    Int_t  fReadOutCard;      // 0 for DARC, 1 for ZRC
    UInt_t fEvType;	      // Event type
    Int_t  fPosition;	      // bit position in buffer data word
    
    // Boolean variables indicating data type
    Bool_t fIsCalib;	      // True when calibration run
    Bool_t fIsDARCHeader;     // True when DARC header
    Bool_t fIsHeaderMapping;  // True when reading header mapping
    Bool_t fIsChMapping;      // True when reading ch. mapping
    Bool_t fIsADCDataWord;    // True when data word
    Bool_t fIsADCHeader;      // True when ADC header
    Bool_t fIsADCEOB;	      // True when EOB
    Bool_t fSODReading;	      // True when reading SOD (DA)
    Bool_t fIsMapRead; 	      // True if map is already read
    Bool_t fReadCDH;          // False for sim raw data (uncorrect CDH!)
    
    Int_t  fDeadfaceOffset;   // deadface offset
    Int_t  fDeadbeefOffset;   // deadbeef offset
    Int_t  fDataOffset;       // data offset
        
    // ADC signal
    Int_t  fSector[2];    // [detector, sector]
    Int_t  fModType;	  // Module type
    Int_t  fADCModule;    // ADC module = GEO address for scaler, trigger card, P.U.
    Int_t  fADCNChannels; // number of ADC ch.
    Int_t  fADCChannel;   // ADC channel = ch. for scaler, trigger card, P.U.
    Int_t  fADCValue;	  // ADC channel
    Int_t  fADCGain;	  // ADC gain (0=high range; 1=low range)
    Bool_t fIsUnderflow;  // ADC underflow
    Bool_t fIsOverflow;   // ADC overflow
    
    // Scaler
    UInt_t fScGeo;           // scaler GEO address
    UInt_t fScNWords;        // no. of words in scaler event
    UInt_t fScTriggerSource; // Trigger source 
    UInt_t fScTriggerNumber; // no. of triggers
    Bool_t fIsScEventGood;   // true if scaler event is good
    Bool_t fIsScHeaderRead;  // true if scaler header is read
    Int_t  fScStartCounter;  // position in the buffer where scaler data begins
    UInt_t fScEvCounter;     // event counter
    Bool_t fIsScalerWord;    // is scaler word (not header)
    
    // Pattern Unit
    UInt_t fDetPattern;  // word from the pattern unit
    
    // Trigger card
    // (1) trigger counts
    Int_t  fTrigCountNWords;  // no. of words to read from trigger card scalers
    Bool_t fIsTriggerScaler;  // Trigger card scalers - 1st word read
    Int_t  fTrigCountStart;   // Trigger card scalers - counter
    Int_t  fMBTrigInput;      // MB          trigger input to trigger card
    Int_t  fCentralTrigInput; // CENTRAL     trigger input to trigger card
    Int_t  fSCentralTrigInput;// SEMICENTRAL trigger input to trigger card
    Int_t  fEMDTrigInput;     // EMD	     trigger input to trigger card
    Int_t  fL0Received;       // L0 received by the trigger card
    Int_t  fMBtrig2CTP;       // trigger input to the CTP for MB
    Int_t  fCentralTrig2CTP;  // trigger input to the CTP for CENTRAL
    Int_t  fSCentralTrig2CTP; // trigger input to the CTP for SEMICENTRAL
    Int_t  fEMDTrig2CTP;      // trigger input to the CTP for EMD
    // (2) trigger history
    Int_t  fTrigHistNWords;   // no. of words to read from trigger history data
    Bool_t fIsTriggerHistory; // Trigger history - 1st word read
    Int_t  fTrigHistStart;    // Trigger card history - counter
    Int_t  fPileUpBit1stWord; // Pile up bit from 1st word
    Int_t  fL0Bit1stWord;     // L0 bit from 1st word
    UInt_t fCentralTrigHist;  // history for CENTRAL trigger
    UInt_t fMBTrigHist;       // history for CENTRAL trigger
    Int_t  fPileUpBit2ndWord; // Pile up bit from 2nd word
    Int_t  fL0Bit2ndWord;     // L0 bit from 2nd word
    UInt_t fSCentralTrigHist; // history for SEMICENTRAL trigger
    UInt_t fEMDTrigHist;      // history for EMD trigger
    Int_t  fCPTInput[4];      // Trigger sent to the CTP
    
    // Channel mapping 
    Int_t fNChannelsOn;      // No. of signals/ADC ch. used
    Int_t fCurrentCh;	     // current mapped ADC ch.
    Int_t fCabledSignal;     // physics signal (from enum)
    Int_t fMapADC[48][5];    // ADC map {ADC mod., ch., signal, det., sec.}
    Int_t fCurrScCh;	     // current mapped scaler ch.
    Int_t fScalerMap[32][5]; // Scaler map {Scaler mod., ch., signal, det., sec.}
    Int_t fCurrTDCCh;	     // current mapped TDC ch.
    Int_t fTDCMap[32][3];    // TDC map {Scaler mod., ch., signal}
    
    // Checks over raw data event quality
    Bool_t fIsADCEventGood; // true if not valid datum not corrupted
    Bool_t fIsL0BitSet;     // true if L0 bit in history words = 1 
    Bool_t fIsPileUpEvent;  // true if pile up bits in history words = 0
    
    // ADD ADC
    Bool_t fIsADDChannel;   // true if datum is an ADD ADC channel
    Int_t  fADDADCdatum;    // ADD ADC datum
    
    // TDCs
    Bool_t fIsTDCHeaderRead; // true when streaming one of the 2 TDCs
    Int_t  fTDCStartCounter; // counts after a TDC header
    //
    Bool_t fIsZDCTDCHeader;  // true if datum is a ZDC TDC header
    Bool_t fIsZDCTDCdatum;   // true if the current is a TDC datum
    Int_t  fZDCTDCdatum;     // datum for ZDC TDC
    //
    Bool_t fIsADDTDCHeader;  // true if datum is an ADD TDC channel
    Bool_t fIsADDTDCdatum;   // true when streaming ADD TDC data
    Int_t  fADDTDCdatum;     // datum for ADD TDC
   
    ClassDef(AliZDCRawStream, 20)    // class for reading ZDC raw data
};

#endif
