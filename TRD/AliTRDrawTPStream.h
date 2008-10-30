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


#define N_REGS       433      // number of conf. registers
#define N_CMD         16      // number of command registers
#define N_RO          25      // number of command registers
#define N_BLOCKS      38      // number of blocks
#define N_BLOCKS_n    64      // number of blocks
#define N_PACKD_DAT 0xE0      // the max size of the packed conf., the absolute max is 256!
#define N_NO_BCST      3      // number of regs without broadcast
#define N_DMEM     0x400      // number of DMEM words
#define N_DBANK    0x100      // number of DBANK words
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
     Char_t     * fName;  //! Name of the register 
     UInt_t	fAddr;    // Address in GIO of TRAP
     UInt_t	fNbits;   // Number of bits, from 1 to 32
     UInt_t	fResVal;  // reset value [mj]
  };

  struct CmdRegs{
     Char_t     * fName;  //! Name of the command register
     UInt_t   	fAddr;    // Address in GIO of TRAP
  };


  Bool_t	DecodeTPdata();
  Bool_t	DecodeConfigdata();
  Bool_t	FillConfig();
  Int_t		ReadPacked(UInt_t *word, UInt_t *pdata, Int_t *len);
  Int_t		UnPackConfN(UInt_t *pData, Int_t maxLength);
  Int_t		SetU(UInt_t addr, UInt_t newVal);
  Int_t		AddrIsDmem(UInt_t addr);
  Int_t		AddrIsDbank(UInt_t addr);
  UInt_t	Addr2Idx(UInt_t addr);
  Char_t	* Addr2Name(UInt_t addr); //!
  Char_t	CnfStat(UInt_t prop);
  void		PowerUp();
  void		DumpCnf(Int_t slv);
  
  enum DbankProp {kDbankEmpty=0, kDbankHeader, kDbankData, kDbankNoB, kDbankCrc32, kDbankEheader, kScsnDat}; 

  SimpleRegs fTrapReg[N_REGS];       // all TRAP configuration registers 
  CmdRegs    fCmdReg[N_CMD];         // all TRAP command registers
  CmdRegs    fRoReg[N_RO];           // all TRAP command registers


  //--------------------------------------------------------
  AliTRDrawTPStream(Int_t rawVMajorOpt, UInt_t * pPos);
  AliTRDrawTPStream(const AliTRDrawTPStream& st);
  AliTRDrawTPStream &operator=(const AliTRDrawTPStream &);
  virtual ~AliTRDrawTPStream();
  //--------------------------------------------------------


 protected:

  UInt_t     	fCnfPro[N_REGS];
  UInt_t     	fDmemValid[N_DMEM];  // 0- empty, 1- valid
  UInt_t    	fRegs[N_REGS];       // the actual content of all conf. registers
  UInt_t    	fDmem[N_DMEM];       // content of the DMEM, in GIO from 0xC000 to 0xC3FF
  UInt_t    	fDbank[N_DBANK];     // 32 bit data, to be send to DBANK
  DbankProp   	fDbankPro[N_DBANK];  // property: 0-empty, 1- header, 2- data, 3- data no broadcast, 4- crc-32 checksum, 5- empty header

  UInt_t	*fpPos;              //! current position in the buffer
  Int_t		fRawVMajorOpt;       // Raw data version


  ClassDef(AliTRDrawTPStream, 0)     // Pattern generated TRD raw data

}; 

#endif
