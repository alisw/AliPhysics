#ifndef ALITRDRAWTPSTREAM_H
#define ALITRDRAWTPSTREAM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDrawTPStream.h 27696 2008-07-31 09:18:53Z cblume $ */

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
// This class provides access to pattern generated TRD raw data including            //
// configuration data.                                                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#define NREGS       433      // number of conf. registers
#define NCMD         16      // number of command registers
#define NRO          25      // number of command registers
#define N_BLOCKS      38      // number of blocks
#define N_BLOCKS_n    64      // number of blocks
#define N_PACKD_DAT 0xE0      // the max size of the packed conf., the absolute max is 256!
#define N_NO_BCST      3      // number of regs without broadcast
#define NDMEM     0x400      // number of DMEM words
#define NDBANK    0x100      // number of DBANK words
#define N_IMEM    0x1000      // number of IMEM words/CPU
#define IMEM_EMPTY 0x80000000 // mark for empty IMEM
#define DBANK_ADDR 0xF000     // start address of DBANK in GIO
#define DMEM_ADDR  0xC000     // start address of DMEM in GIO
#define ENDM_CONF  0x7FFF00FE // end marker for the packed configuration

#include "TObject.h"
#include "TString.h"
#include "AliTRDrawStreamBase.h"



class AliTRDrawTPStream : public AliTRDrawStreamBase
{ // class def begin

 public:

  struct SimpleRegs {
     const Char_t     * fName;  //! Name of the register 
     UInt_t	fAddr;    // Address in GIO of TRAP
     UInt_t	fNbits;   // Number of bits, from 1 to 32
     UInt_t	fResVal;  // reset value [mj]
  };

  struct CmdRegs{
     const Char_t     * fName;  //! Name of the command register
     UInt_t   	fAddr;    // Address in GIO of TRAP
  };

  AliTRDrawTPStream(Int_t rawVMajorOpt, UInt_t * pPos);
  AliTRDrawTPStream(const AliTRDrawTPStream& st);
  AliTRDrawTPStream &operator=(const AliTRDrawTPStream &);
  virtual ~AliTRDrawTPStream();

  Bool_t	DecodeTPdata();
  Bool_t	DecodeConfigdata();
  Bool_t	FillConfig();
  Int_t		ReadPacked(UInt_t *word, UInt_t *pdata, Int_t * const len);
  Int_t		UnPackConfN(const UInt_t *pData, Int_t maxLength);
  Int_t		SetU(UInt_t addr, UInt_t newVal);
  Int_t		AddrIsDmem(UInt_t addr) const;
  Int_t		AddrIsDbank(UInt_t addr) const;
  UInt_t	Addr2Idx(UInt_t addr) const;
  const Char_t	* Addr2Name(UInt_t addr) const; //!
  Char_t	CnfStat(UInt_t prop) const;
  void		PowerUp();
  void		DumpCnf(Int_t slv);
  


 protected:

  enum DbankProp {kDbankEmpty=0, kDbankHeader, kDbankData, kDbankNoB, kDbankCrc32, kDbankEheader, kScsnDat}; 

  SimpleRegs fTrapReg[NREGS];       // all TRAP configuration registers 
  CmdRegs    fCmdReg[NCMD];         // all TRAP command registers
  CmdRegs    fRoReg[NRO];           // all TRAP command registers

  UInt_t     	fCnfPro[NREGS];     // something ...
  UInt_t     	fDmemValid[NDMEM];  // 0- empty, 1- valid
  UInt_t    	fRegs[NREGS];       // the actual content of all conf. registers
  UInt_t    	fDmem[NDMEM];       // content of the DMEM, in GIO from 0xC000 to 0xC3FF
  UInt_t    	fDbank[NDBANK];     // 32 bit data, to be send to DBANK
  DbankProp   	fDbankPro[NDBANK];  // property: 0-empty, 1- header, 2- data, 3- data no broadcast, 4- crc-32 checksum, 5- empty header

  UInt_t	*fpPos;              //! current position in the buffer
  Int_t		fRawVMajorOpt;       // Raw data version


  ClassDef(AliTRDrawTPStream, 0)     // Pattern generated TRD raw data

}; 

#endif
