#ifndef ALITPCMONITORALTRO_H
#define ALITPCMONITORALTRO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorAltro class
////
//// Class for decoding raw TPC data in the ALTRO format
//// 
//// Authors: Roland Bramm, 
////          Stefan Kniege, IKF, Frankfurt
////       
/////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TNamed.h>

class AliTPCMonitorAltro : public TNamed {
 public:
    AliTPCMonitorAltro(UInt_t* memory, Int_t size, Int_t fformat);
    AliTPCMonitorAltro(const  AliTPCMonitorAltro &altro);
    AliTPCMonitorAltro& operator= (const AliTPCMonitorAltro& altro);
    ~AliTPCMonitorAltro();
    
    void         Allocate40BitArray(); 
    void         Allocate10BitArray();
    
    void         Decodeto40Bit();
    void         Decodeto10Bit(Int_t equipment = -1); 

    Int_t        DecodeTrailer(Int_t pos);
    Int_t        DecodeTrailerVbb(Int_t pos);
    
    Long64_t*    Get40BitArray();
    Short_t     *Get10BitArray();     
    Int_t        Get40BitArraySize()     const { return fmemory[fsize-GetRCUTrailerSize()];}   //hier �ndern
    Int_t        Get10BitArraySize()     const { return fmemory[fsize-GetRCUTrailerSize()]*4;}  //number of 10 bit words from trailer
    const Char_t* GetActFileName()        const { return ffilename.Data();}
    
    static Int_t GetHwMaskFEC()                { return fgkHwMaskFEC;}
    static Int_t GetHwMaskBranch()             { return fgkHwMaskBranch;}
    static Int_t GetHwMaskFECChannel()         { return fgkHwMaskFECChannel;}
    static Int_t GetHwMaskAltroChannel()       { return fgkHwMaskAltroChannel;}
    static Int_t GetHwMaskAltroChip()          { return fgkHwMaskAltroChip;}

    static Int_t GetHwMaskRCU()                { return fgkHwMaskRCU;}
    
    Int_t        GetNextTrailerPos()     const { return fNextPos;}

    Int_t        GetTrailerNWords()      const { return fTrailerNWords   ;}
    Int_t        GetTrailerHwAddress()   const { return fTrailerHwAddress;}
    Int_t        GetTrailerDataPos()     const { return fTrailerDataPos  ;}
    Int_t        GetTrailerBlockPos()    const { return fTrailerBlockPos ;}
    Int_t        GetTrailerPos()         const { return fTrailerPos      ;} 

  Int_t        GetRCUTrailerSize()     const { Int_t ts=(GetAltroVersion()==0xaaaa||GetAltroVersion()==0xaabb)*
      (fmemory[fsize-1]&0x3F); return (ts>0)?ts:1;}
  UInt_t        GetAltroVersion()       const { return fmemory[fsize-1]>>16; }

    void         SetDataOffset(Int_t val){ foffset     =val ;} 
    void         SetWrite10Bit(Int_t wr) { fwrite10bit =wr  ;}
    
  
    void         SetActFilename(const Char_t* name){ ffilename=name; }
    void         SetVerbose(Int_t val)   { fverb=val;}
     
 private:
    
    Int_t                    fverb;                                                     // verbose flag           
    UInt_t*                  fmemory;                                                   // memory pointer fo payload               
    Int_t                    fsize;                                                     // size of fmemory 
    Long64_t*                f40BitArray;                                               // array to store 40 bit words 
    Short_t*                 f10BitArray;                                               // array to store 10 bit words   
    Int_t                    fdecoderPos;                                               // start position for decoding 40 bit words
         
    Bool_t                   fallocate40BitArray;                                       // flag for decoding to 40 bit words
    Bool_t                   fallocate10BitArray;                                       // flag for decoding to 10 bit words
    Int_t                    foffset ;                                                   // data offset (CDH length) 
    Int_t                    fwrite10bit;                                               // flag for writing 10 bit words to file 
    
    Int_t                    fTrailerNWords ;                                           // from Trailer: number of 40 bit words for channel  
    Int_t                    fTrailerHwAddress;                                         // from Trailer: hardware address for current channel
    Int_t                    fTrailerDataPos;                                           // from Trailer: position of first adc value 
    Int_t                    fTrailerBlockPos;                                          // from Trailer: number of 40 bit words for channel
    Int_t                    fTrailerPos;                                               // trailer position

    Int_t                    fNextPos;                                                  // position of next trailer
    TString                  ffilename;                                                 // name of processed file
    
    static const Int_t       fgk24BitOn                = 16777215;                        // bit masks for first 24 bits of 32  for decoding 32 bit words
    static const Int_t       fgk16BitOn                = 65535;                           // bit masks for first 24 bits of 24
    static const Int_t       fgk08BitOn                = 255;                             // bit masks for first 24 bits of 8
    
    
    static const Long64_t    fgkmask10                 = (Long64_t)0x00000000000003FFULL; // mask first   10 bit out of 4o0 bit word 
    static const Long64_t    fgkmask20                 = (Long64_t)0x00000000000FFC00ULL; // mask second  10 bit out of 4o0 bit word 
    static const Long64_t    fgkmask30                 = (Long64_t)0x000000003FF00000ULL; // mask third   10 bit out of 4o0 bit word 
    static const Long64_t    fgkmask40                 = (Long64_t)0x000000FFC0000000ULL; // mask fourth  10 bit out of 4o0 bit word 
    
    static const Long64_t    fgkTrailerTail            = (Long64_t)0x0000000000002AAAULL; // Tail of the Trailer set to 2AAA 
    static const Long64_t    fgkTrailerTailErr         = (Long64_t)0x0000000000002AEEULL; // Tail of the Trailer set to 2AEE if an error occured
    static const Long64_t    fgkTrailerMaskTail        = (Long64_t)0x000000fffC000000ULL; // mask for trailer
    static const Long64_t    fgkTrailerMaskHardw       = (Long64_t)0x0000000000000FFFULL; // mask for hardware address
    static const Long64_t    fgkTrailerMaskNWords      = (Long64_t)0x0000000003FF0000ULL; // mask for nwords  (number of 40 bit data words)
    
    static const Int_t       fgkHwMaskFEC              = 0x0780;                          // mask for fec in hardware address
    static const Int_t       fgkHwMaskBranch           = 0x0800;                          // mask for branch in hardware address
    static const Int_t       fgkHwMaskFECChannel       = 0x007f;                          // mask for fec channel  in hardware address
    static const Int_t       fgkHwMaskAltroChannel     = 0x000f;                          // mask for altro channel in hardware address
    static const Int_t       fgkHwMaskAltroChip        = 0x0070;                          // mask for altro chip  in hardware address
    static const Int_t       fgkHwMaskRCU              = 0x7000;                          // not part of the trailer added afterwards
    
    ClassDef(AliTPCMonitorAltro,1);
};
#endif
