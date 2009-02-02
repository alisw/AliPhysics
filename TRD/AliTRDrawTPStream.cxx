#/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/* $Id: AliTRDrawTPStream.cxx 27797 2008-08-05 14:37:22Z cblume $ */

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
// This class provides access to pattern generated TRD raw data including            //
// configuration data.                                                               //
//                                                                                   //
// It is based on Venelin Angelov's c++ code decoding standalone                     // 
// configuration data                                                                //  
// http://alice.physi.uni-heidelberg.de/svn/trd/wconfigurations/trunk/C/trap_cnf.cpp //
// http://alice.physi.uni-heidelberg.de/svn/trd/wconfigurations/trunk/C/trap_cnf.h   //
//                                                                                   //
// Author: MinJung Kweon(minjung@physi.uni-heidelberg.de)                            // 
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "TString.h"
#include "AliLog.h"

#include "AliTRDrawStream.h"
#include "AliTRDrawTPStream.h"


#define GET_VALUE_AT(w,m,s) (( (w) >> (s)) & (m) )
#define MCM_HEADER_MASK_ERR(w) ( ((w) & (0xf)) == (0xc) ? 0 : 1) 
#define MCM_ROB_NUMBER(w) GET_VALUE_AT(w,0x7,28)
#define MCM_MCM_NUMBER(w) GET_VALUE_AT(w,0x0f,24)
#define MCM_EVENT_COUNTER(w) GET_VALUE_AT(w,0x00fffff,4)



ClassImp(AliTRDrawTPStream)

//---------------------------------------------------------------------
AliTRDrawTPStream::AliTRDrawTPStream(Int_t rawVMajorOpt, UInt_t * pPos)
  : AliTRDrawStreamBase()
  , fTrapReg() 
  , fCmdReg() 
  , fRoReg() 
  , fCnfPro()
  , fDmemValid()
  , fRegs()
  , fDmem()
  , fDbank()
  , fDbankPro()
  , fpPos(pPos)
  , fRawVMajorOpt(rawVMajorOpt) 
{
  //
  // default constructor
  //

  if (FillConfig() == kFALSE)
    AliError("Reading reset value failed.");

}

//---------------------------------------------------------------------
AliTRDrawTPStream::AliTRDrawTPStream(const AliTRDrawTPStream& /*st*/)
  : AliTRDrawStreamBase()
  , fTrapReg() 
  , fCmdReg() 
  , fRoReg() 
  , fCnfPro()
  , fDmemValid()
  , fRegs()
  , fDmem()
  , fDbank()
  , fDbankPro()
  , fpPos()
  , fRawVMajorOpt() 
{
  //
  // copy constructor
  //

  AliError("Not implemeneted.");

}

//---------------------------------------------------------------------
AliTRDrawTPStream &
AliTRDrawTPStream::operator=(const AliTRDrawTPStream &)
{
  //
  // we are not using this functionality
  //
  AliFatal("May not use.");
  return *this;
}

//---------------------------------------------------------------------
AliTRDrawTPStream::~AliTRDrawTPStream()
{
  //
  // destructor
  //
}

//---------------------------------------------------------------------
Bool_t AliTRDrawTPStream::DecodeTPdata()
{

  if (fRawVMajorOpt == 7)
    {
     AliInfo("This is configuration data event read by first trigger.");
     if(!AliTRDrawStream::fgEnableDecodeConfigData) return kTRUE;
     if (DecodeConfigdata() == kFALSE) // configuration data 
       {
        AliError("failed to to decode configuration data");
        return kFALSE;
       }
     else 
       return kTRUE;
    }
  else
    AliError("These are different type of test pattern data. You need other reader");

  return kFALSE;
}

//---------------------------------------------------------------------
Bool_t AliTRDrawTPStream::DecodeConfigdata()
{

    UInt_t packedConf[256];
    Int_t mcmPos, mcmsRead, lengthPacked;

    mcmsRead = 0;
    do
    {
        mcmPos = ReadPacked(fpPos, packedConf, &lengthPacked);
        if (mcmPos >= 0)
        {
            PowerUp();
            UnPackConfN(packedConf, lengthPacked);
            DumpCnf(mcmPos);
            mcmsRead++;
            AliInfo(Form("%d MCMs read up to now, last was MCM%02d\n",mcmsRead, mcmPos));
        }
    } while ((mcmsRead < 84) && (mcmPos >= 0)); // [mj] have to think about # of mcmsRead
    AliInfo("Done\n");

    return kTRUE;
}

//---------------------------------------------------------------------
Int_t AliTRDrawTPStream::ReadPacked(UInt_t *word, UInt_t *pData, Int_t *nWords)
{

    UInt_t vword = *word;

    Int_t  iLength;
    UInt_t err, robNum, mcmNum, chipId, NoEndMarker;

    iLength = 0;
    err = 0;

    // decode mcm header
    if(MCM_HEADER_MASK_ERR(vword)) err++;

    robNum = MCM_ROB_NUMBER(vword);
    mcmNum = MCM_MCM_NUMBER(vword);
    chipId = MCM_EVENT_COUNTER(vword);

    if (err == 0) {
      AliInfo(Form("MCM header ROB %d, MCM %02d, ChipId %d 0x%05x\n", robNum, mcmNum, chipId, chipId));
    }
    else 
      return -1;

    // read MCM data and store into array
    NoEndMarker = 1;
    do
    {
        word++;
        vword = *word;

        NoEndMarker = ((vword != ENDM_CONF) && (vword != (ENDM_CONF | 1)) && (vword != 0x10001000));
        *pData = vword;
        pData++;
        iLength++;
    } while (NoEndMarker && (iLength < 256));

    word++;       
    fpPos = word;

    *nWords = iLength;
    if (iLength == 0) 
      return -1;
    else
      return mcmNum;
}

//---------------------------------------------------------------------
void AliTRDrawTPStream::PowerUp() // power up
{
    // copy the reset values 
    for (Int_t i=0; i< N_REGS; i++)
    {
        fRegs[i] = fTrapReg[i].fResVal;
        fCnfPro[i] = 0;
    }
    // mark all DMEM cells as invalid
    for (Int_t i=0; i< N_DMEM; i++) fDmemValid[i] = 0;
    // mark all DBANK cells as empty
    for (Int_t i=0; i< N_DBANK; i++) fDbankPro[i] = kDbankEmpty;
}


//---------------------------------------------------------------------
Int_t AliTRDrawTPStream::UnPackConfN(UInt_t *pData, Int_t maxLength)
{
    Int_t debug = 0; // the debug mode not completely ready
    Int_t step, bwidth, nwords, idx, err, exitFlag, bitcnt, werr;
    UInt_t caddr;
    UInt_t dat, msk, header, dataHi;

    idx = 0; // index in PackedConf
    err = 0;
    while (idx < maxLength)
    {
        header = *pData;
        if (debug) printf("read 0x%08x  ",header);
        pData++;
        idx++;
        if (header & 0x01) // single data
          {
            dat = (header >> 2) & 0xFFFF;       // 16 bit data
            caddr = (header >> 18) & 0x3FFF;    // 14 bit address
            if (caddr != 0x1FFF)                // temp!!! because the end marker was wrong
            {
             if (header & 0x02)                 // check if > 16 bits
               {
                dataHi = *pData;
                if (debug) printf("read 0x%08x  ",dataHi);
                pData++;
                idx++;
                err += ((dataHi ^ (dat | 1)) & 0xFFFF) != 0;
                dat = (dataHi & 0xFFFF0000) | dat;
               }
               if (debug) printf("addr=0x%04x (%s) data=0x%08x\n",caddr, Addr2Name(caddr), dat);
               werr = SetU(caddr, dat);
               if (werr < 0)
                 {
                  printf("(single-write): non-existing address 0x%04x containing 0x%08x\n", caddr, header);
                 }
               if (idx > maxLength)
                 {
                  printf("(single-write): no more data, missing end marker\n");
                  return -err;
                 }
            }
            else
            {
             printf("(single-write): address 0x%04x => old endmarker?\n",caddr);
                return err;
            }
          }
        else               // block of data
          {
            step   =  (header >>  1) & 0x0003;
            bwidth = ((header >>  3) & 0x001F) + 1;
            nwords =  (header >>  8) & 0x00FF;
            caddr  =  (header >> 16) & 0xFFFF;
            exitFlag = (step == 0) || (step == 3) || (nwords == 0);
            if (exitFlag) return err;
            switch (bwidth)
            {
                case    15:
                case    10:
                case     7:
                case     6:
                case     5:
                {
                    msk = (1 << bwidth) - 1;
                    bitcnt = 0;
                    while (nwords > 0)
                    {
                        nwords--;
                        bitcnt -= bwidth;
                        if (bitcnt < 0)
                        {
                            header = *pData;
                            if (debug) printf("read 0x%08x  ",header);
                            pData++;
                            idx++;
                            err += (header & 1);
                            header = header >> 1;
                            bitcnt = 31 - bwidth;
                        }
                        if (debug) printf("addr=0x%04x (%s) data=0x%08x\n",caddr, Addr2Name(caddr), header & msk);
                        werr = SetU(caddr, header & msk);
                        if (werr < 0)
                        {
                          printf("(single-write): non-existing address 0x%04x containing 0x%08x\n", caddr, header);
                        }
                        caddr += step;
                        header = header >> bwidth;
                        if (idx >= maxLength)
                        {
                          printf("(block-write): no end marker! %d words read\n",idx);
                          return -err;
                        }
                    }
                    break;
                } // end case 5-15
                case 31:
                {
                    while (nwords > 0)
                    {
                        header = *pData;
                        if (debug) printf("read 0x%08x  ",header);
                        pData++;
                        idx++;
                        nwords--;
                        err += (header & 1);
                        if (debug) printf("addr=0x%04x (%s) data=0x%08x\n",caddr, Addr2Name(caddr), header >> 1);
                        werr = SetU(caddr, header >> 1);
                        if (werr < 0)
                        {
                            printf("(single-write): non-existing address 0x%04x containing 0x%08x\n", caddr, header);
                        }
                        caddr += step;
                        if (idx >= maxLength)
                        {
                            printf("no end marker! %d words read\n",idx);
                            return -err;
                        }
                    }
                    break;
                }
                default: return err;
            } // end switch
        } // end block case
    } // end while
    printf("no end marker! %d words read\n",idx);
    return -err; // only if the max length of the block reached!
}

//---------------------------------------------------------------------
void AliTRDrawTPStream::DumpCnf(Int_t slv)
{
    UInt_t idx;
    for (idx = 0; idx < N_REGS; idx++) // config. reg
       {
        if (slv >= 0)
          printf("%s\t0x%08x\t%3d %c\n", fTrapReg[idx].fName, (Int_t) fRegs[idx], slv, CnfStat(fCnfPro[idx]));
        else
          printf("%s\t0x%08x %c\n", fTrapReg[idx].fName, (Int_t) fRegs[idx], CnfStat(fCnfPro[idx]));
       }
}


//---------------------------------------------------------------------
const Char_t * AliTRDrawTPStream::Addr2Name(UInt_t addr)
{
    Int_t idx;
    idx = 0;
    if ( ( ( (addr >> 4) & 0xFFE) == 0x0C0) && ( ( (addr >> 2) & 1) == 1) )
    {
        addr = addr & 0x0C07;
    }
    while ((idx < N_REGS) && (fTrapReg[idx].fAddr != addr) ) idx++;
    if (idx < N_REGS)
        return fTrapReg[idx].fName;
    idx = 0;
    while ((idx < N_CMD) && (fCmdReg[idx].fAddr != addr)) idx++;
    if (idx < N_CMD)
        return fCmdReg[idx].fName;
    idx = 0;
    while ((idx < N_RO) && (fRoReg[idx].fAddr != addr)) idx++;
    if (idx < N_RO)
        return fRoReg[idx].fName;
    else
        return 0;
}

//---------------------------------------------------------------------
Char_t AliTRDrawTPStream::CnfStat(UInt_t prop)
{
    if (prop == 0) return 'U';
    else
    if (prop == 1) return 'R';
    else
    if (prop == 2) return 'I';
    else
                   return prop;
}

//---------------------------------------------------------------------
Int_t AliTRDrawTPStream::SetU(UInt_t addr, UInt_t newVal)
{
    Int_t i;
    UInt_t maxVal = 0;

    if (AddrIsDmem(addr))
    {
        fDmem[addr & 0x3FF] = newVal;
        fDmemValid[addr & 0x3FF] = 1;
        return 0;
    }
    else
    if (AddrIsDbank(addr))
    {
        fDbank[addr & 0xFF] = newVal;
        fDbankPro[addr & 0xFF] = kScsnDat;
        return 0;
    }
    else
    {
        i = Addr2Idx(addr);
        if (i < N_REGS) // found
        {
            fCnfPro[i] = 2;
            if (fTrapReg[i].fNbits < 32) // create the max value from the number of bits
            {
                maxVal = 1;
                maxVal = (maxVal << fTrapReg[i].fNbits) - 1;
            }
            if ( (fTrapReg[i].fNbits == 32) || (newVal <= maxVal) ) // in range
            {
                fRegs[i] = newVal;
                return 0;
            }
            else
            {   // out of range
                fRegs[i] = newVal & maxVal;
                printf("Out of range, writing 0x%08x to %d bits at addr = 0x%04x\n",newVal, fTrapReg[i].fNbits, addr);
                return -2;
            }
        }
        else    // not found
            {
                printf("(SetU): No such address, writing 0x%08x to addr = 0x%04x\n",newVal, addr);
                return -1; // no such address
            }
    }
}

//---------------------------------------------------------------------
Int_t AliTRDrawTPStream::AddrIsDmem(UInt_t addr)
{
    addr = (addr >> 10);
    return (addr == 0x30);
}

//---------------------------------------------------------------------
Int_t AliTRDrawTPStream::AddrIsDbank(UInt_t addr)
{
    addr = (addr >> 8);
    return (addr == 0xF0);
}

//---------------------------------------------------------------------
UInt_t AliTRDrawTPStream::Addr2Idx(UInt_t addr)
{
    Int_t idx;
    idx = 0;
    // check if global const
    if ( ( ( (addr >> 4) & 0xFFE) == 0x0C0) && ( ( (addr >> 2) & 1) == 1) )
    {
        addr = addr & 0x0C07;
    }
    // searching
    while ((idx < N_REGS) && (fTrapReg[idx].fAddr != addr)) idx++;
       // printf("Addr = 0x%04x; Idx = %d\n",addr, idx); // debugging
    return idx;
}

//---------------------------------------------------------------------
Bool_t AliTRDrawTPStream::FillConfig()
{

  const SimpleRegs trapReg[N_REGS] = {
    // Name         Address Nbits   Reset Value
    // Global state machine
    {"SML0",        0x0A00, 15,     0x4050},
    {"SML1",        0x0A01, 15,     0x4200},
    {"SML2",        0x0A02, 15,     0x4384},
    {"SMMODE",      0x0A03, 16,     0xF0E2},
    {"NITM0",       0x0A08, 14,     0x3FFF},
    {"NITM1",       0x0A09, 14,     0x3FFF},
    {"NITM2",       0x0A0A, 14,     0x3FFF},
    {"NIP4D",       0x0A0B, 7,      0x7F},
    {"CPU0CLK",     0x0A20, 5,      0x07},
    {"CPU1CLK",     0x0A22, 5,      0x07},
    {"CPU2CLK",     0x0A24, 5,      0x07},
    {"CPU3CLK",     0x0A26, 5,      0x07},
    {"NICLK",       0x0A28, 5,      0x07},
    {"FILCLK",      0x0A2A, 5,      0x07},
    {"PRECLK",      0x0A2C, 5,      0x07},
    {"ADCEN",       0x0A2E, 5,      0x07},
    {"NIODE",       0x0A30, 5,      0x07},
    {"NIOCE",       0x0A32, 5,      0x21}, // bit 5 is status bit (read-only)!
    {"NIIDE",       0x0A34, 5,      0x07},
    {"NIICE",       0x0A36, 5,      0x07},
    // Arbiter
    {"ARBTIM",      0x0A3F, 4,      0x0},
    // IVT of CPU0
    {"IA0IRQ0",     0x0B00, 12,     0x000},
    {"IA0IRQ1",     0x0B01, 12,     0x000},
    {"IA0IRQ2",     0x0B02, 12,     0x000},
    {"IA0IRQ3",     0x0B03, 12,     0x000},
    {"IA0IRQ4",     0x0B04, 12,     0x000},
    {"IA0IRQ5",     0x0B05, 12,     0x000},
    {"IA0IRQ6",     0x0B06, 12,     0x000},
    {"IA0IRQ7",     0x0B07, 12,     0x000},
    {"IA0IRQ8",     0x0B08, 12,     0x000},
    {"IA0IRQ9",     0x0B09, 12,     0x000},
    {"IA0IRQA",     0x0B0A, 12,     0x000},
    {"IA0IRQB",     0x0B0B, 12,     0x000},
    {"IA0IRQC",     0x0B0C, 12,     0x000},
    {"IRQSW0",      0x0B0D, 13,     0x1FFF},
    {"IRQHW0",      0x0B0E, 13,     0x0000},
    {"IRQHL0",      0x0B0F, 13,     0x0000},
    // IVT of CPU1
    {"IA1IRQ0",     0x0B20, 12,     0x000},
    {"IA1IRQ1",     0x0B21, 12,     0x000},
    {"IA1IRQ2",     0x0B22, 12,     0x000},
    {"IA1IRQ3",     0x0B23, 12,     0x000},
    {"IA1IRQ4",     0x0B24, 12,     0x000},
    {"IA1IRQ5",     0x0B25, 12,     0x000},
    {"IA1IRQ6",     0x0B26, 12,     0x000},
    {"IA1IRQ7",     0x0B27, 12,     0x000},
    {"IA1IRQ8",     0x0B28, 12,     0x000},
    {"IA1IRQ9",     0x0B29, 12,     0x000},
    {"IA1IRQA",     0x0B2A, 12,     0x000},
    {"IA1IRQB",     0x0B2B, 12,     0x000},
    {"IA1IRQC",     0x0B2C, 12,     0x000},
    {"IRQSW1",      0x0B2D, 13,     0x1FFF},
    {"IRQHW1",      0x0B2E, 13,     0x0000},
    {"IRQHL1",      0x0B2F, 13,     0x0000},
    // IVT of CPU2
    {"IA2IRQ0",     0x0B40, 12,     0x000},
    {"IA2IRQ1",     0x0B41, 12,     0x000},
    {"IA2IRQ2",     0x0B42, 12,     0x000},
    {"IA2IRQ3",     0x0B43, 12,     0x000},
    {"IA2IRQ4",     0x0B44, 12,     0x000},
    {"IA2IRQ5",     0x0B45, 12,     0x000},
    {"IA2IRQ6",     0x0B46, 12,     0x000},
    {"IA2IRQ7",     0x0B47, 12,     0x000},
    {"IA2IRQ8",     0x0B48, 12,     0x000},
    {"IA2IRQ9",     0x0B49, 12,     0x000},
    {"IA2IRQA",     0x0B4A, 12,     0x000},
    {"IA2IRQB",     0x0B4B, 12,     0x000},
    {"IA2IRQC",     0x0B4C, 12,     0x000},
    {"IRQSW2",      0x0B4D, 13,     0x1FFF},
    {"IRQHW2",      0x0B4E, 13,     0x0000},
    {"IRQHL2",      0x0B4F, 13,     0x0000},
    // IVT of CPU3
    {"IA3IRQ0",     0x0B60, 12,     0x000},
    {"IA3IRQ1",     0x0B61, 12,     0x000},
    {"IA3IRQ2",     0x0B62, 12,     0x000},
    {"IA3IRQ3",     0x0B63, 12,     0x000},
    {"IA3IRQ4",     0x0B64, 12,     0x000},
    {"IA3IRQ5",     0x0B65, 12,     0x000},
    {"IA3IRQ6",     0x0B66, 12,     0x000},
    {"IA3IRQ7",     0x0B67, 12,     0x000},
    {"IA3IRQ8",     0x0B68, 12,     0x000},
    {"IA3IRQ9",     0x0B69, 12,     0x000},
    {"IA3IRQA",     0x0B6A, 12,     0x000},
    {"IA3IRQB",     0x0B6B, 12,     0x000},
    {"IA3IRQC",     0x0B6C, 12,     0x000},
    {"IRQSW3",      0x0B6D, 13,     0x1FFF},
    {"IRQHW3",      0x0B6E, 13,     0x0000},
    {"IRQHL3",      0x0B6F, 13,     0x0000},
    // Global Counter/Timer
    {"CTGDINI",     0x0B80, 32,     0x00000000},
    {"CTGCTRL",     0x0B81, 12,     0xE3F},
    // CPU constants
    {"C08CPU0",     0x0C00, 32,     0x00000000},
    {"C09CPU0",     0x0C01, 32,     0x00000000},
    {"C10CPU0",     0x0C02, 32,     0x00000000},
    {"C11CPU0",     0x0C03, 32,     0x00000000},
    {"C12CPUA",     0x0C04, 32,     0x00000000},
    {"C13CPUA",     0x0C05, 32,     0x00000000},
    {"C14CPUA",     0x0C06, 32,     0x00000000},
    {"C15CPUA",     0x0C07, 32,     0x00000000},
    {"C08CPU1",     0x0C08, 32,     0x00000000},
    {"C09CPU1",     0x0C09, 32,     0x00000000},
    {"C10CPU1",     0x0C0A, 32,     0x00000000},
    {"C11CPU1",     0x0C0B, 32,     0x00000000},
    {"C08CPU2",     0x0C10, 32,     0x00000000},
    {"C09CPU2",     0x0C11, 32,     0x00000000},
    {"C10CPU2",     0x0C12, 32,     0x00000000},
    {"C11CPU2",     0x0C13, 32,     0x00000000},
    {"C08CPU3",     0x0C18, 32,     0x00000000},
    {"C09CPU3",     0x0C19, 32,     0x00000000},
    {"C10CPU3",     0x0C1A, 32,     0x00000000},
    {"C11CPU3",     0x0C1B, 32,     0x00000000},
    // NI interface
    {"NMOD",        0x0D40, 6,      0x08},
    {"NDLY",        0x0D41, 30,     0x24924924},
    {"NED",         0x0D42, 16,     0xA240},
    {"NTRO",        0x0D43, 18,     0x3FFFC},
    {"NRRO",        0x0D44, 18,     0x3FFFC},

    {"NES",         0x0D45, 32,     0x00000000},
    {"NTP",         0x0D46, 32,     0x0000FFFF},
    {"NBND",        0x0D47, 16,     0x6020},
    {"NP0",         0x0D48, 11,     0x44C},
    {"NP1",         0x0D49, 11,     0x44C},
    {"NP2",         0x0D4A, 11,     0x44C},
    {"NP3",         0x0D4B, 11,     0x44C},
    {"NCUT",        0x0D4C, 32,     0xFFFFFFFF},
    // Filter and Preprocessor
    {"TPPT0",       0x3000, 7,      0x01},
    {"TPFS",        0x3001, 7,      0x05},
    {"TPFE",        0x3002, 7,      0x14},
    {"TPPGR",       0x3003, 7,      0x15},
    {"TPPAE",       0x3004, 7,      0x1E},
    {"TPQS0",       0x3005, 7,      0x00},
    {"TPQE0",       0x3006, 7,      0x0A},
    {"TPQS1",       0x3007, 7,      0x0B},
    {"TPQE1",       0x3008, 7,      0x14},
    {"EBD",         0x3009, 3,      0x0},
    {"EBAQA",       0x300A, 7,      0x00},
    {"EBSIA",       0x300B, 7,      0x20},
    {"EBSF",        0x300C, 1,      0x1},
    {"EBSIM",       0x300D, 1,      0x1},
    {"EBPP",        0x300E, 1,      0x1},
    {"EBPC",        0x300F, 1,      0x1},

    {"EBIS",        0x3014, 10,     0x005},
    {"EBIT",        0x3015, 12,     0x028},
    {"EBIL",        0x3016, 8,      0xF0},
    {"EBIN",        0x3017, 1,      0x1},
    {"FLBY",        0x3018, 1,      0x0},
    {"FPBY",        0x3019, 1,      0x0},
    {"FGBY",        0x301A, 1,      0x0},
    {"FTBY",        0x301B, 1,      0x0},
    {"FCBY",        0x301C, 1,      0x0},
    {"FPTC",        0x3020, 2,      0x3},
    {"FPNP",        0x3021, 9,      0x078},
    {"FPCL",        0x3022, 1,      0x1},
    {"FGTA",        0x3028, 12,     0x014},
    {"FGTB",        0x3029, 12,     0x80C},
    {"FGCL",        0x302A, 1,      0x1},
    {"FTAL",        0x3030, 10,     0x0F6},
    {"FTLL",        0x3031, 9,      0x11D},
    {"FTLS",        0x3032, 9,      0x0D3},
    {"FCW1",        0x3038, 8,      0x1E},
    {"FCW2",        0x3039, 8,      0xD4},
    {"FCW3",        0x303A, 8,      0xE6},
    {"FCW4",        0x303B, 8,      0x4A},
    {"FCW5",        0x303C, 8,      0xEF},
    {"TPFP",        0x3040, 9,      0x037},
    {"TPHT",        0x3041, 14,     0x00A0},

    {"TPVT",        0x3042, 6,      0x00},
    {"TPVBY",       0x3043, 1,      0x0},
    {"TPCT",        0x3044, 5,      0x08},
    {"TPCL",        0x3045, 5,      0x01},
    {"TPCBY",       0x3046, 1,      0x1},
    {"TPD",         0x3047, 4,      0xF},
    {"TPCI0",       0x3048, 5,      0x00},
    {"TPCI1",       0x3049, 5,      0x00},
    {"TPCI2",       0x304A, 5,      0x00},
    {"TPCI3",       0x304B, 5,      0x00},

    {"ADCMSK",      0x3050, 21,     0x1FFFFF},
    {"ADCINB",      0x3051, 2,      0x2},
    {"ADCDAC",      0x3052, 5,      0x10},
    {"ADCPAR",      0x3053, 18,     0x195EF},
    {"ADCTST",      0x3054, 2,      0x0},
    {"SADCAZ",      0x3055, 1,      0x1},

    {"FGF0",        0x3080, 9,      0x000},
    {"FGF1",        0x3081, 9,      0x000},
    {"FGF2",        0x3082, 9,      0x000},
    {"FGF3",        0x3083, 9,      0x000},
    {"FGF4",        0x3084, 9,      0x000},
    {"FGF5",        0x3085, 9,      0x000},
    {"FGF6",        0x3086, 9,      0x000},
    {"FGF7",        0x3087, 9,      0x000},
    {"FGF8",        0x3088, 9,      0x000},
    {"FGF9",        0x3089, 9,      0x000},
    {"FGF10",       0x308A, 9,      0x000},
    {"FGF11",       0x308B, 9,      0x000},
    {"FGF12",       0x308C, 9,      0x000},
    {"FGF13",       0x308D, 9,      0x000},
    {"FGF14",       0x308E, 9,      0x000},
    {"FGF15",       0x308F, 9,      0x000},
    {"FGF16",       0x3090, 9,      0x000},
    {"FGF17",       0x3091, 9,      0x000},
    {"FGF18",       0x3092, 9,      0x000},
    {"FGF19",       0x3093, 9,      0x000},
    {"FGF20",       0x3094, 9,      0x000},

    {"FGA0",        0x30A0, 6,      0x00},
    {"FGA1",        0x30A1, 6,      0x00},
    {"FGA2",        0x30A2, 6,      0x00},
    {"FGA3",        0x30A3, 6,      0x00},
    {"FGA4",        0x30A4, 6,      0x00},
    {"FGA5",        0x30A5, 6,      0x00},
    {"FGA6",        0x30A6, 6,      0x00},
    {"FGA7",        0x30A7, 6,      0x00},
    {"FGA8",        0x30A8, 6,      0x00},
    {"FGA9",        0x30A9, 6,      0x00},
    {"FGA10",       0x30AA, 6,      0x00},
    {"FGA11",       0x30AB, 6,      0x00},
    {"FGA12",       0x30AC, 6,      0x00},
    {"FGA13",       0x30AD, 6,      0x00},
    {"FGA14",       0x30AE, 6,      0x00},
    {"FGA15",       0x30AF, 6,      0x00},
    {"FGA16",       0x30B0, 6,      0x00},
    {"FGA17",       0x30B1, 6,      0x00},
    {"FGA18",       0x30B2, 6,      0x00},
    {"FGA19",       0x30B3, 6,      0x00},
    {"FGA20",       0x30B4, 6,      0x00},
    // non-linearity table, 64 x 6 bits
    {"FLL00",       0x3100, 6,      0x00},
    {"FLL01",       0x3101, 6,      0x00},
    {"FLL02",       0x3102, 6,      0x00},
    {"FLL03",       0x3103, 6,      0x00},
    {"FLL04",       0x3104, 6,      0x00},
    {"FLL05",       0x3105, 6,      0x00},
    {"FLL06",       0x3106, 6,      0x00},
    {"FLL07",       0x3107, 6,      0x00},
    {"FLL08",       0x3108, 6,      0x00},
    {"FLL09",       0x3109, 6,      0x00},
    {"FLL0A",       0x310A, 6,      0x00},
    {"FLL0B",       0x310B, 6,      0x00},
    {"FLL0C",       0x310C, 6,      0x00},
    {"FLL0D",       0x310D, 6,      0x00},
    {"FLL0E",       0x310E, 6,      0x00},
    {"FLL0F",       0x310F, 6,      0x00},
    {"FLL10",       0x3110, 6,      0x00},
    {"FLL11",       0x3111, 6,      0x00},
    {"FLL12",       0x3112, 6,      0x00},
    {"FLL13",       0x3113, 6,      0x00},
    {"FLL14",       0x3114, 6,      0x00},
    {"FLL15",       0x3115, 6,      0x00},
    {"FLL16",       0x3116, 6,      0x00},
    {"FLL17",       0x3117, 6,      0x00},
    {"FLL18",       0x3118, 6,      0x00},
    {"FLL19",       0x3119, 6,      0x00},
    {"FLL1A",       0x311A, 6,      0x00},
    {"FLL1B",       0x311B, 6,      0x00},
    {"FLL1C",       0x311C, 6,      0x00},
    {"FLL1D",       0x311D, 6,      0x00},
    {"FLL1E",       0x311E, 6,      0x00},
    {"FLL1F",       0x311F, 6,      0x00},
    {"FLL20",       0x3120, 6,      0x00},
    {"FLL21",       0x3121, 6,      0x00},
    {"FLL22",       0x3122, 6,      0x00},
    {"FLL23",       0x3123, 6,      0x00},
    {"FLL24",       0x3124, 6,      0x00},
    {"FLL25",       0x3125, 6,      0x00},
    {"FLL26",       0x3126, 6,      0x00},
    {"FLL27",       0x3127, 6,      0x00},
    {"FLL28",       0x3128, 6,      0x00},
    {"FLL29",       0x3129, 6,      0x00},
    {"FLL2A",       0x312A, 6,      0x00},
    {"FLL2B",       0x312B, 6,      0x00},
    {"FLL2C",       0x312C, 6,      0x00},
    {"FLL2D",       0x312D, 6,      0x00},
    {"FLL2E",       0x312E, 6,      0x00},
    {"FLL2F",       0x312F, 6,      0x00},
    {"FLL30",       0x3130, 6,      0x00},
    {"FLL31",       0x3131, 6,      0x00},
    {"FLL32",       0x3132, 6,      0x00},
    {"FLL33",       0x3133, 6,      0x00},
    {"FLL34",       0x3134, 6,      0x00},
    {"FLL35",       0x3135, 6,      0x00},
    {"FLL36",       0x3136, 6,      0x00},
    {"FLL37",       0x3137, 6,      0x00},
    {"FLL38",       0x3138, 6,      0x00},
    {"FLL39",       0x3139, 6,      0x00},
    {"FLL3A",       0x313A, 6,      0x00},
    {"FLL3B",       0x313B, 6,      0x00},
    {"FLL3C",       0x313C, 6,      0x00},
    {"FLL3D",       0x313D, 6,      0x00},
    {"FLL3E",       0x313E, 6,      0x00},
    {"FLL3F",       0x313F, 6,      0x00},
    // end of non-lin table
    {"PASADEL",     0x3158, 8,      0xFF},
    {"PASAPHA",     0x3159, 6,      0x3F},
    {"PASAPRA",     0x315A, 6,      0x0F},
    {"PASADAC",     0x315B, 8,      0x80},
    {"PASACHM",     0x315C, 19,     0x7FFFF},
    {"PASASTL",     0x315D, 8,      0xFF},
    {"PASAPR1",     0x315E, 1,      0x0},
    {"PASAPR0",     0x315F, 1,      0x0},
    {"SADCTRG",     0x3161, 1,      0x0},
    {"SADCRUN",     0x3162, 1,      0x0},
    {"SADCPWR",     0x3163, 3,      0x7},
    {"L0TSIM",      0x3165, 14,     0x0050},
    {"SADCEC",      0x3166, 7,      0x00},
    {"SADCMC",      0x3170, 8,      0xC0},
    {"SADCOC",      0x3171, 8,      0x19},
    {"SADCGTB",     0x3172, 32,     0x37737700},
    {"SEBDEN",      0x3178, 3,      0x0},
    {"SEBDOU",      0x3179, 3,      0x0},
    // pos table, 128 x 5 bits
    {"TPL00",       0x3180, 5,      0x00},
    {"TPL01",       0x3181, 5,      0x00},
    {"TPL02",       0x3182, 5,      0x00},
    {"TPL03",       0x3183, 5,      0x00},
    {"TPL04",       0x3184, 5,      0x00},
    {"TPL05",       0x3185, 5,      0x00},
    {"TPL06",       0x3186, 5,      0x00},
    {"TPL07",       0x3187, 5,      0x00},
    {"TPL08",       0x3188, 5,      0x00},
    {"TPL09",       0x3189, 5,      0x00},
    {"TPL0A",       0x318A, 5,      0x00},
    {"TPL0B",       0x318B, 5,      0x00},
    {"TPL0C",       0x318C, 5,      0x00},
    {"TPL0D",       0x318D, 5,      0x00},
    {"TPL0E",       0x318E, 5,      0x00},
    {"TPL0F",       0x318F, 5,      0x00},
    {"TPL10",       0x3190, 5,      0x00},
    {"TPL11",       0x3191, 5,      0x00},
    {"TPL12",       0x3192, 5,      0x00},
    {"TPL13",       0x3193, 5,      0x00},
    {"TPL14",       0x3194, 5,      0x00},
    {"TPL15",       0x3195, 5,      0x00},
    {"TPL16",       0x3196, 5,      0x00},
    {"TPL17",       0x3197, 5,      0x00},
    {"TPL18",       0x3198, 5,      0x00},
    {"TPL19",       0x3199, 5,      0x00},
    {"TPL1A",       0x319A, 5,      0x00},
    {"TPL1B",       0x319B, 5,      0x00},
    {"TPL1C",       0x319C, 5,      0x00},
    {"TPL1D",       0x319D, 5,      0x00},
    {"TPL1E",       0x319E, 5,      0x00},
    {"TPL1F",       0x319F, 5,      0x00},
    {"TPL20",       0x31A0, 5,      0x00},
    {"TPL21",       0x31A1, 5,      0x00},
    {"TPL22",       0x31A2, 5,      0x00},
    {"TPL23",       0x31A3, 5,      0x00},
    {"TPL24",       0x31A4, 5,      0x00},
    {"TPL25",       0x31A5, 5,      0x00},
    {"TPL26",       0x31A6, 5,      0x00},
    {"TPL27",       0x31A7, 5,      0x00},
    {"TPL28",       0x31A8, 5,      0x00},
    {"TPL29",       0x31A9, 5,      0x00},
    {"TPL2A",       0x31AA, 5,      0x00},
    {"TPL2B",       0x31AB, 5,      0x00},
    {"TPL2C",       0x31AC, 5,      0x00},
    {"TPL2D",       0x31AD, 5,      0x00},
    {"TPL2E",       0x31AE, 5,      0x00},
    {"TPL2F",       0x31AF, 5,      0x00},
    {"TPL30",       0x31B0, 5,      0x00},
    {"TPL31",       0x31B1, 5,      0x00},
    {"TPL32",       0x31B2, 5,      0x00},
    {"TPL33",       0x31B3, 5,      0x00},
    {"TPL34",       0x31B4, 5,      0x00},
    {"TPL35",       0x31B5, 5,      0x00},
    {"TPL36",       0x31B6, 5,      0x00},
    {"TPL37",       0x31B7, 5,      0x00},
    {"TPL38",       0x31B8, 5,      0x00},
    {"TPL39",       0x31B9, 5,      0x00},
    {"TPL3A",       0x31BA, 5,      0x00},
    {"TPL3B",       0x31BB, 5,      0x00},
    {"TPL3C",       0x31BC, 5,      0x00},
    {"TPL3D",       0x31BD, 5,      0x00},
    {"TPL3E",       0x31BE, 5,      0x00},
    {"TPL3F",       0x31BF, 5,      0x00},
    {"TPL40",       0x31C0, 5,      0x00},
    {"TPL41",       0x31C1, 5,      0x00},
    {"TPL42",       0x31C2, 5,      0x00},
    {"TPL43",       0x31C3, 5,      0x00},
    {"TPL44",       0x31C4, 5,      0x00},
    {"TPL45",       0x31C5, 5,      0x00},
    {"TPL46",       0x31C6, 5,      0x00},
    {"TPL47",       0x31C7, 5,      0x00},
    {"TPL48",       0x31C8, 5,      0x00},
    {"TPL49",       0x31C9, 5,      0x00},
    {"TPL4A",       0x31CA, 5,      0x00},
    {"TPL4B",       0x31CB, 5,      0x00},
    {"TPL4C",       0x31CC, 5,      0x00},
    {"TPL4D",       0x31CD, 5,      0x00},
    {"TPL4E",       0x31CE, 5,      0x00},
    {"TPL4F",       0x31CF, 5,      0x00},
    {"TPL50",       0x31D0, 5,      0x00},
    {"TPL51",       0x31D1, 5,      0x00},
    {"TPL52",       0x31D2, 5,      0x00},
    {"TPL53",       0x31D3, 5,      0x00},
    {"TPL54",       0x31D4, 5,      0x00},
    {"TPL55",       0x31D5, 5,      0x00},
    {"TPL56",       0x31D6, 5,      0x00},
    {"TPL57",       0x31D7, 5,      0x00},
    {"TPL58",       0x31D8, 5,      0x00},
    {"TPL59",       0x31D9, 5,      0x00},
    {"TPL5A",       0x31DA, 5,      0x00},
    {"TPL5B",       0x31DB, 5,      0x00},
    {"TPL5C",       0x31DC, 5,      0x00},
    {"TPL5D",       0x31DD, 5,      0x00},
    {"TPL5E",       0x31DE, 5,      0x00},
    {"TPL5F",       0x31DF, 5,      0x00},
    {"TPL60",       0x31E0, 5,      0x00},
    {"TPL61",       0x31E1, 5,      0x00},
    {"TPL62",       0x31E2, 5,      0x00},
    {"TPL63",       0x31E3, 5,      0x00},
    {"TPL64",       0x31E4, 5,      0x00},
    {"TPL65",       0x31E5, 5,      0x00},
    {"TPL66",       0x31E6, 5,      0x00},
    {"TPL67",       0x31E7, 5,      0x00},
    {"TPL68",       0x31E8, 5,      0x00},
    {"TPL69",       0x31E9, 5,      0x00},
    {"TPL6A",       0x31EA, 5,      0x00},
    {"TPL6B",       0x31EB, 5,      0x00},
    {"TPL6C",       0x31EC, 5,      0x00},
    {"TPL6D",       0x31ED, 5,      0x00},
    {"TPL6E",       0x31EE, 5,      0x00},
    {"TPL6F",       0x31EF, 5,      0x00},
    {"TPL70",       0x31F0, 5,      0x00},
    {"TPL71",       0x31F1, 5,      0x00},
    {"TPL72",       0x31F2, 5,      0x00},
    {"TPL73",       0x31F3, 5,      0x00},
    {"TPL74",       0x31F4, 5,      0x00},
    {"TPL75",       0x31F5, 5,      0x00},
    {"TPL76",       0x31F6, 5,      0x00},
    {"TPL77",       0x31F7, 5,      0x00},
    {"TPL78",       0x31F8, 5,      0x00},
    {"TPL79",       0x31F9, 5,      0x00},
    {"TPL7A",       0x31FA, 5,      0x00},
    {"TPL7B",       0x31FB, 5,      0x00},
    {"TPL7C",       0x31FC, 5,      0x00},
    {"TPL7D",       0x31FD, 5,      0x00},
    {"TPL7E",       0x31FE, 5,      0x00},
    {"TPL7F",       0x31FF, 5,      0x00},
    // end of pos table
    {"MEMRW",       0xD000, 7,      0x79},
    {"MEMCOR",      0xD001, 9,      0x000},
    {"DMDELA",      0xD002, 4,      0x8},
    {"DMDELS",      0xD003, 4,      0x8}
  };

  const CmdRegs cmdReg[N_CMD] = {
    // Name      Address
    {"SMCMD"   , 0x0A04},
    {"SMOFFON" , 0x0A05},
    {"SMON"    , 0x0A06},
    {"SMOFF"   , 0x0A07},
    {"CPU0SS"  , 0x0A21},
    {"CPU1SS"  , 0x0A23},
    {"CPU2SS"  , 0x0A25},
    {"CPU3SS"  , 0x0A27},
    {"NICLKSS" , 0x0A29},
    {"FILCLKSS", 0x0A2B},
    {"PRECLKSS", 0x0A2D},
    {"ADCENSS" , 0x0A2F},
    {"NIODESS" , 0x0A31},
    {"NIOCESS" , 0x0A33},
    {"NIIDESS" , 0x0A35},
    {"NIICESS" , 0x0A37}
  };

  const CmdRegs roReg[N_RO] = {
    // NI
    {"NCTRL"  , 0x0DC0},
    {"NFE"    , 0x0DC1},
    {"NFSM"   , 0x0DC2},
    // event buffer parity violation counters
    {"EBP0"   , 0x3010},
    {"EBP1"   , 0x3011},
    {"EBP2"   , 0x3012},
    {"EBP3"   , 0x3013},
    // slow ADC
    {"SADCC0" , 0x3168},
    {"SADCC1" , 0x3169},
    {"SADCC2" , 0x316A},
    {"SADCC3" , 0x316B},
    {"SADCC4" , 0x316C},
    {"SADCC5" , 0x316D},
    {"SADCC6" , 0x316E},
    {"SADCC7" , 0x316F},
    // hamming counters
    {"HCNTI0" , 0xD010},
    {"HCNTI1" , 0xD011},
    {"HCNTI2" , 0xD012},
    {"HCNTI3" , 0xD013},
    {"HCNTD0" , 0xD014},
    {"HCNTD1" , 0xD015},
    {"HCNTD2" , 0xD016},
    {"HCNTD3" , 0xD017},

    {"CHIPID" , 0x3160},

    {"SEBDIN" , 0x317A}
  };


  for (Int_t i = 0; i < N_REGS; i++) {
     fTrapReg[i] = trapReg[i];
  }
  for (Int_t i = 0; i < N_CMD; i++) {
     fCmdReg[i] = cmdReg[i];
  }
  for (Int_t i = 0; i < N_RO; i++) {
     fRoReg[i] = roReg[i];
  }

  return kTRUE;
}
