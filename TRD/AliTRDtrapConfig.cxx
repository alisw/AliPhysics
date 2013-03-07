/**************************************************************************
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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRAP config                                                           //
//                                                                        //
//  Author: J. Klein (Jochen.Klein@cern.ch)                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"

#include "AliTRDgeometry.h"
#include "AliTRDfeeParam.h"
#include "AliTRDtrapConfig.h"

#include <fstream>
#include <iostream>
#include <iomanip>

const Int_t AliTRDtrapConfig::AliTRDtrapValue::fgkSize[] = {
  0,
  1,
  540,
  1080,
  8*18*540,
  4,
  6,
  8*18*30
};
Bool_t AliTRDtrapConfig::fgRegAddressMapInitialized = kFALSE;
AliTRDtrapConfig::TrapReg_t AliTRDtrapConfig::fgRegAddressMap[] = { };
const Int_t AliTRDtrapConfig::fgkRegisterAddressBlockStart[] = { 0x0a00, 0x3000, 0xd000 };
const Int_t AliTRDtrapConfig::fgkRegisterAddressBlockSize[]  = { 0x0400, 0x0200, 0x0004 };

AliTRDtrapConfig::AliTRDtrapConfig(const TString &name, const TString &title) :
  TNamed(name, title)
{
  // default constructor

  // initialize and reset the TRAP registers
  InitRegs();
  ResetRegs();

  for (Int_t iWord = 0; iWord < fgkDmemWords; ++iWord) {
    fDmem[iWord].SetAddress(iWord + fgkDmemStartAddress);
  }

  // initialize the map from address to register
  if (!fgRegAddressMapInitialized) {
    for (Int_t iReg = 0; iReg < kLastReg; iReg++) {
      Int_t addr = fRegisterValue[iReg].GetAddr();
      if (addr < fgkRegisterAddressBlockStart[0]) {
	AliError(Form("Register address 0x%04x not handled in register map", addr));
      }
      else if (addr < fgkRegisterAddressBlockStart[0] + fgkRegisterAddressBlockSize[0]) {
	fgRegAddressMap[addr - fgkRegisterAddressBlockStart[0]] = (TrapReg_t) iReg;
      }
      else if (addr < fgkRegisterAddressBlockStart[1]) {
	AliError(Form("Register address 0x%04x not handled in register map", addr));
      }
      else if (addr < fgkRegisterAddressBlockStart[1] + fgkRegisterAddressBlockSize[1]) {
	fgRegAddressMap[addr - fgkRegisterAddressBlockStart[1] + fgkRegisterAddressBlockSize[0]] = (TrapReg_t) iReg;
      }
      else if (addr < fgkRegisterAddressBlockStart[2]) {
	AliError(Form("Register address 0x%04x not handled in register map", addr));
      }
      else if (addr < fgkRegisterAddressBlockStart[2] + fgkRegisterAddressBlockSize[2]) {
	Int_t ind = addr - fgkRegisterAddressBlockStart[2] + fgkRegisterAddressBlockSize[1] + fgkRegisterAddressBlockSize[0];
	fgRegAddressMap[ind] = (TrapReg_t) iReg;
      }
      else {
	AliError(Form("Register address 0x%04x not handled in register map", addr));
      }
    }
    fgRegAddressMapInitialized = kTRUE;
  }
}


AliTRDtrapConfig::~AliTRDtrapConfig()
{
  // destructor
}


void AliTRDtrapConfig::InitRegs()
{
  // initialize all TRAP registers

  //                              Name          Address  Nbits   Reset Value
  fRegisterValue[kSML0]    .Init("SML0",        0x0A00, 15,     0x4050     );  // Global state machine
  fRegisterValue[kSML1]    .Init("SML1",        0x0A01, 15,     0x4200     );
  fRegisterValue[kSML2]    .Init("SML2",        0x0A02, 15,     0x4384     );
  fRegisterValue[kSMMODE]  .Init("SMMODE",      0x0A03, 16,     0xF0E2     );
  fRegisterValue[kSMCMD]   .Init("SMCMD",       0x0A04, 16,     0x0000     );
  fRegisterValue[kNITM0]   .Init("NITM0",       0x0A08, 14,     0x3FFF     );
  fRegisterValue[kNITM1]   .Init("NITM1",       0x0A09, 14,     0x3FFF     );
  fRegisterValue[kNITM2]   .Init("NITM2",       0x0A0A, 14,     0x3FFF     );
  fRegisterValue[kNIP4D]   .Init("NIP4D",       0x0A0B, 7,      0x7F       );
  fRegisterValue[kCPU0CLK] .Init("CPU0CLK",     0x0A20, 5,      0x07       );
  fRegisterValue[kCPU1CLK] .Init("CPU1CLK",     0x0A22, 5,      0x07       );
  fRegisterValue[kCPU2CLK] .Init("CPU2CLK",     0x0A24, 5,      0x07       );
  fRegisterValue[kCPU3CLK] .Init("CPU3CLK",     0x0A26, 5,      0x07       );
  fRegisterValue[kNICLK]   .Init("NICLK",       0x0A28, 5,      0x07       );
  fRegisterValue[kFILCLK]  .Init("FILCLK",      0x0A2A, 5,      0x07       );
  fRegisterValue[kPRECLK]  .Init("PRECLK",      0x0A2C, 5,      0x07       );
  fRegisterValue[kADCEN]   .Init("ADCEN",       0x0A2E, 5,      0x07       );
  fRegisterValue[kNIODE]   .Init("NIODE",       0x0A30, 5,      0x07       );
  fRegisterValue[kNIOCE]   .Init("NIOCE",       0x0A32, 6,      0x21       );  // bit 5 is status bit (read-only)!
  fRegisterValue[kNIIDE]   .Init("NIIDE",       0x0A34, 5,      0x07       );
  fRegisterValue[kNIICE]   .Init("NIICE",       0x0A36, 5,      0x07       );
  fRegisterValue[kARBTIM]  .Init("ARBTIM",      0x0A3F, 4,      0x0        );  // Arbiter
  fRegisterValue[kIA0IRQ0] .Init("IA0IRQ0",     0x0B00, 12,     0x000      );  // IVT of CPU0
  fRegisterValue[kIA0IRQ1] .Init("IA0IRQ1",     0x0B01, 12,     0x000      );
  fRegisterValue[kIA0IRQ2] .Init("IA0IRQ2",     0x0B02, 12,     0x000      );
  fRegisterValue[kIA0IRQ3] .Init("IA0IRQ3",     0x0B03, 12,     0x000      );
  fRegisterValue[kIA0IRQ4] .Init("IA0IRQ4",     0x0B04, 12,     0x000      );
  fRegisterValue[kIA0IRQ5] .Init("IA0IRQ5",     0x0B05, 12,     0x000      );
  fRegisterValue[kIA0IRQ6] .Init("IA0IRQ6",     0x0B06, 12,     0x000      );
  fRegisterValue[kIA0IRQ7] .Init("IA0IRQ7",     0x0B07, 12,     0x000      );
  fRegisterValue[kIA0IRQ8] .Init("IA0IRQ8",     0x0B08, 12,     0x000      );
  fRegisterValue[kIA0IRQ9] .Init("IA0IRQ9",     0x0B09, 12,     0x000      );
  fRegisterValue[kIA0IRQA] .Init("IA0IRQA",     0x0B0A, 12,     0x000      );
  fRegisterValue[kIA0IRQB] .Init("IA0IRQB",     0x0B0B, 12,     0x000      );
  fRegisterValue[kIA0IRQC] .Init("IA0IRQC",     0x0B0C, 12,     0x000      );
  fRegisterValue[kIRQSW0]  .Init("IRQSW0",      0x0B0D, 13,     0x1FFF     );
  fRegisterValue[kIRQHW0]  .Init("IRQHW0",      0x0B0E, 13,     0x0000     );
  fRegisterValue[kIRQHL0]  .Init("IRQHL0",      0x0B0F, 13,     0x0000     );
  fRegisterValue[kIA1IRQ0] .Init("IA1IRQ0",     0x0B20, 12,     0x000      );  // IVT of CPU1
  fRegisterValue[kIA1IRQ1] .Init("IA1IRQ1",     0x0B21, 12,     0x000      );
  fRegisterValue[kIA1IRQ2] .Init("IA1IRQ2",     0x0B22, 12,     0x000      );
  fRegisterValue[kIA1IRQ3] .Init("IA1IRQ3",     0x0B23, 12,     0x000      );
  fRegisterValue[kIA1IRQ4] .Init("IA1IRQ4",     0x0B24, 12,     0x000      );
  fRegisterValue[kIA1IRQ5] .Init("IA1IRQ5",     0x0B25, 12,     0x000      );
  fRegisterValue[kIA1IRQ6] .Init("IA1IRQ6",     0x0B26, 12,     0x000      );
  fRegisterValue[kIA1IRQ7] .Init("IA1IRQ7",     0x0B27, 12,     0x000      );
  fRegisterValue[kIA1IRQ8] .Init("IA1IRQ8",     0x0B28, 12,     0x000      );
  fRegisterValue[kIA1IRQ9] .Init("IA1IRQ9",     0x0B29, 12,     0x000      );
  fRegisterValue[kIA1IRQA] .Init("IA1IRQA",     0x0B2A, 12,     0x000      );
  fRegisterValue[kIA1IRQB] .Init("IA1IRQB",     0x0B2B, 12,     0x000      );
  fRegisterValue[kIA1IRQC] .Init("IA1IRQC",     0x0B2C, 12,     0x000      );
  fRegisterValue[kIRQSW1]  .Init("IRQSW1",      0x0B2D, 13,     0x1FFF     );
  fRegisterValue[kIRQHW1]  .Init("IRQHW1",      0x0B2E, 13,     0x0000     );
  fRegisterValue[kIRQHL1]  .Init("IRQHL1",      0x0B2F, 13,     0x0000     );
  fRegisterValue[kIA2IRQ0] .Init("IA2IRQ0",     0x0B40, 12,     0x000      );  // IVT of CPU2
  fRegisterValue[kIA2IRQ1] .Init("IA2IRQ1",     0x0B41, 12,     0x000      );
  fRegisterValue[kIA2IRQ2] .Init("IA2IRQ2",     0x0B42, 12,     0x000      );
  fRegisterValue[kIA2IRQ3] .Init("IA2IRQ3",     0x0B43, 12,     0x000      );
  fRegisterValue[kIA2IRQ4] .Init("IA2IRQ4",     0x0B44, 12,     0x000      );
  fRegisterValue[kIA2IRQ5] .Init("IA2IRQ5",     0x0B45, 12,     0x000      );
  fRegisterValue[kIA2IRQ6] .Init("IA2IRQ6",     0x0B46, 12,     0x000      );
  fRegisterValue[kIA2IRQ7] .Init("IA2IRQ7",     0x0B47, 12,     0x000      );
  fRegisterValue[kIA2IRQ8] .Init("IA2IRQ8",     0x0B48, 12,     0x000      );
  fRegisterValue[kIA2IRQ9] .Init("IA2IRQ9",     0x0B49, 12,     0x000      );
  fRegisterValue[kIA2IRQA] .Init("IA2IRQA",     0x0B4A, 12,     0x000      );
  fRegisterValue[kIA2IRQB] .Init("IA2IRQB",     0x0B4B, 12,     0x000      );
  fRegisterValue[kIA2IRQC] .Init("IA2IRQC",     0x0B4C, 12,     0x000      );
  fRegisterValue[kIRQSW2]  .Init("IRQSW2",      0x0B4D, 13,     0x1FFF     );
  fRegisterValue[kIRQHW2]  .Init("IRQHW2",      0x0B4E, 13,     0x0000     );
  fRegisterValue[kIRQHL2]  .Init("IRQHL2",      0x0B4F, 13,     0x0000     );
  fRegisterValue[kIA3IRQ0] .Init("IA3IRQ0",     0x0B60, 12,     0x000      );  // IVT of CPU3
  fRegisterValue[kIA3IRQ1] .Init("IA3IRQ1",     0x0B61, 12,     0x000      );
  fRegisterValue[kIA3IRQ2] .Init("IA3IRQ2",     0x0B62, 12,     0x000      );
  fRegisterValue[kIA3IRQ3] .Init("IA3IRQ3",     0x0B63, 12,     0x000      );
  fRegisterValue[kIA3IRQ4] .Init("IA3IRQ4",     0x0B64, 12,     0x000      );
  fRegisterValue[kIA3IRQ5] .Init("IA3IRQ5",     0x0B65, 12,     0x000      );
  fRegisterValue[kIA3IRQ6] .Init("IA3IRQ6",     0x0B66, 12,     0x000      );
  fRegisterValue[kIA3IRQ7] .Init("IA3IRQ7",     0x0B67, 12,     0x000      );
  fRegisterValue[kIA3IRQ8] .Init("IA3IRQ8",     0x0B68, 12,     0x000      );
  fRegisterValue[kIA3IRQ9] .Init("IA3IRQ9",     0x0B69, 12,     0x000      );
  fRegisterValue[kIA3IRQA] .Init("IA3IRQA",     0x0B6A, 12,     0x000      );
  fRegisterValue[kIA3IRQB] .Init("IA3IRQB",     0x0B6B, 12,     0x000      );
  fRegisterValue[kIA3IRQC] .Init("IA3IRQC",     0x0B6C, 12,     0x000      );
  fRegisterValue[kIRQSW3]  .Init("IRQSW3",      0x0B6D, 13,     0x1FFF     );
  fRegisterValue[kIRQHW3]  .Init("IRQHW3",      0x0B6E, 13,     0x0000     );
  fRegisterValue[kIRQHL3]  .Init("IRQHL3",      0x0B6F, 13,     0x0000     );
  fRegisterValue[kCTGDINI] .Init("CTGDINI",     0x0B80, 32,     0x00000000 );  // Global Counter/Timer
  fRegisterValue[kCTGCTRL] .Init("CTGCTRL",     0x0B81, 12,     0xE3F      );
  fRegisterValue[kC08CPU0] .Init("C08CPU0",     0x0C00, 32,     0x00000000 );  // CPU constants
  fRegisterValue[kC09CPU0] .Init("C09CPU0",     0x0C01, 32,     0x00000000 );
  fRegisterValue[kC10CPU0] .Init("C10CPU0",     0x0C02, 32,     0x00000000 );
  fRegisterValue[kC11CPU0] .Init("C11CPU0",     0x0C03, 32,     0x00000000 );
  fRegisterValue[kC12CPUA] .Init("C12CPUA",     0x0C04, 32,     0x00000000 );
  fRegisterValue[kC13CPUA] .Init("C13CPUA",     0x0C05, 32,     0x00000000 );
  fRegisterValue[kC14CPUA] .Init("C14CPUA",     0x0C06, 32,     0x00000000 );
  fRegisterValue[kC15CPUA] .Init("C15CPUA",     0x0C07, 32,     0x00000000 );
  fRegisterValue[kC08CPU1] .Init("C08CPU1",     0x0C08, 32,     0x00000000 );
  fRegisterValue[kC09CPU1] .Init("C09CPU1",     0x0C09, 32,     0x00000000 );
  fRegisterValue[kC10CPU1] .Init("C10CPU1",     0x0C0A, 32,     0x00000000 );
  fRegisterValue[kC11CPU1] .Init("C11CPU1",     0x0C0B, 32,     0x00000000 );
  fRegisterValue[kC08CPU2] .Init("C08CPU2",     0x0C10, 32,     0x00000000 );
  fRegisterValue[kC09CPU2] .Init("C09CPU2",     0x0C11, 32,     0x00000000 );
  fRegisterValue[kC10CPU2] .Init("C10CPU2",     0x0C12, 32,     0x00000000 );
  fRegisterValue[kC11CPU2] .Init("C11CPU2",     0x0C13, 32,     0x00000000 );
  fRegisterValue[kC08CPU3] .Init("C08CPU3",     0x0C18, 32,     0x00000000 );
  fRegisterValue[kC09CPU3] .Init("C09CPU3",     0x0C19, 32,     0x00000000 );
  fRegisterValue[kC10CPU3] .Init("C10CPU3",     0x0C1A, 32,     0x00000000 );
  fRegisterValue[kC11CPU3] .Init("C11CPU3",     0x0C1B, 32,     0x00000000 );
  fRegisterValue[kNMOD]    .Init("NMOD",        0x0D40, 6,      0x08       );  // NI interface
  fRegisterValue[kNDLY]    .Init("NDLY",        0x0D41, 30,     0x24924924 );
  fRegisterValue[kNED]     .Init("NED",         0x0D42, 16,     0xA240     );
  fRegisterValue[kNTRO]    .Init("NTRO",        0x0D43, 18,     0x3FFFC    );
  fRegisterValue[kNRRO]    .Init("NRRO",        0x0D44, 18,     0x3FFFC    );
  fRegisterValue[kNES]     .Init("NES",         0x0D45, 32,     0x00000000 );
  fRegisterValue[kNTP]     .Init("NTP",         0x0D46, 32,     0x0000FFFF );
  fRegisterValue[kNBND]    .Init("NBND",        0x0D47, 16,     0x6020     );
  fRegisterValue[kNP0]     .Init("NP0",         0x0D48, 11,     0x44C      );
  fRegisterValue[kNP1]     .Init("NP1",         0x0D49, 11,     0x44C      );
  fRegisterValue[kNP2]     .Init("NP2",         0x0D4A, 11,     0x44C      );
  fRegisterValue[kNP3]     .Init("NP3",         0x0D4B, 11,     0x44C      );
  fRegisterValue[kNCUT]    .Init("NCUT",        0x0D4C, 32,     0xFFFFFFFF );
  fRegisterValue[kTPPT0]   .Init("TPPT0",       0x3000, 7,      0x01       );  // Filter and Preprocessor
  fRegisterValue[kTPFS]    .Init("TPFS",        0x3001, 7,      0x05       );
  fRegisterValue[kTPFE]    .Init("TPFE",        0x3002, 7,      0x14       );
  fRegisterValue[kTPPGR]   .Init("TPPGR",       0x3003, 7,      0x15       );
  fRegisterValue[kTPPAE]   .Init("TPPAE",       0x3004, 7,      0x1E       );
  fRegisterValue[kTPQS0]   .Init("TPQS0",       0x3005, 7,      0x00       );
  fRegisterValue[kTPQE0]   .Init("TPQE0",       0x3006, 7,      0x0A       );
  fRegisterValue[kTPQS1]   .Init("TPQS1",       0x3007, 7,      0x0B       );
  fRegisterValue[kTPQE1]   .Init("TPQE1",       0x3008, 7,      0x14       );
  fRegisterValue[kEBD]     .Init("EBD",         0x3009, 3,      0x0        );
  fRegisterValue[kEBAQA]   .Init("EBAQA",       0x300A, 7,      0x00       );
  fRegisterValue[kEBSIA]   .Init("EBSIA",       0x300B, 7,      0x20       );
  fRegisterValue[kEBSF]    .Init("EBSF",        0x300C, 1,      0x1        );
  fRegisterValue[kEBSIM]   .Init("EBSIM",       0x300D, 1,      0x1        );
  fRegisterValue[kEBPP]    .Init("EBPP",        0x300E, 1,      0x1        );
  fRegisterValue[kEBPC]    .Init("EBPC",        0x300F, 1,      0x1        );
  fRegisterValue[kEBIS]    .Init("EBIS",        0x3014, 10,     0x005      );
  fRegisterValue[kEBIT]    .Init("EBIT",        0x3015, 12,     0x028      );
  fRegisterValue[kEBIL]    .Init("EBIL",        0x3016, 8,      0xF0       );
  fRegisterValue[kEBIN]    .Init("EBIN",        0x3017, 1,      0x1        );
  fRegisterValue[kFLBY]    .Init("FLBY",        0x3018, 1,      0x0        );
  fRegisterValue[kFPBY]    .Init("FPBY",        0x3019, 1,      0x0        );
  fRegisterValue[kFGBY]    .Init("FGBY",        0x301A, 1,      0x0        );
  fRegisterValue[kFTBY]    .Init("FTBY",        0x301B, 1,      0x0        );
  fRegisterValue[kFCBY]    .Init("FCBY",        0x301C, 1,      0x0        );
  fRegisterValue[kFPTC]    .Init("FPTC",        0x3020, 2,      0x3        );
  fRegisterValue[kFPNP]    .Init("FPNP",        0x3021, 9,      0x078      );
  fRegisterValue[kFPCL]    .Init("FPCL",        0x3022, 1,      0x1        );
  fRegisterValue[kFGTA]    .Init("FGTA",        0x3028, 12,     0x014      );
  fRegisterValue[kFGTB]    .Init("FGTB",        0x3029, 12,     0x80C      );
  fRegisterValue[kFGCL]    .Init("FGCL",        0x302A, 1,      0x1        );
  fRegisterValue[kFTAL]    .Init("FTAL",        0x3030, 10,     0x0F6      );
  fRegisterValue[kFTLL]    .Init("FTLL",        0x3031, 9,      0x11D      );
  fRegisterValue[kFTLS]    .Init("FTLS",        0x3032, 9,      0x0D3      );
  fRegisterValue[kFCW1]    .Init("FCW1",        0x3038, 8,      0x1E       );
  fRegisterValue[kFCW2]    .Init("FCW2",        0x3039, 8,      0xD4       );
  fRegisterValue[kFCW3]    .Init("FCW3",        0x303A, 8,      0xE6       );
  fRegisterValue[kFCW4]    .Init("FCW4",        0x303B, 8,      0x4A       );
  fRegisterValue[kFCW5]    .Init("FCW5",        0x303C, 8,      0xEF       );
  fRegisterValue[kTPFP]    .Init("TPFP",        0x3040, 9,      0x037      );
  fRegisterValue[kTPHT]    .Init("TPHT",        0x3041, 14,     0x00A0     );
  fRegisterValue[kTPVT]    .Init("TPVT",        0x3042, 6,      0x00       );
  fRegisterValue[kTPVBY]   .Init("TPVBY",       0x3043, 1,      0x0        );
  fRegisterValue[kTPCT]    .Init("TPCT",        0x3044, 5,      0x08       );
  fRegisterValue[kTPCL]    .Init("TPCL",        0x3045, 5,      0x01       );
  fRegisterValue[kTPCBY]   .Init("TPCBY",       0x3046, 1,      0x1        );
  fRegisterValue[kTPD]     .Init("TPD",         0x3047, 4,      0xF        );
  fRegisterValue[kTPCI0]   .Init("TPCI0",       0x3048, 5,      0x00       );
  fRegisterValue[kTPCI1]   .Init("TPCI1",       0x3049, 5,      0x00       );
  fRegisterValue[kTPCI2]   .Init("TPCI2",       0x304A, 5,      0x00       );
  fRegisterValue[kTPCI3]   .Init("TPCI3",       0x304B, 5,      0x00       );
  fRegisterValue[kADCMSK]  .Init("ADCMSK",      0x3050, 21,     0x1FFFFF   );
  fRegisterValue[kADCINB]  .Init("ADCINB",      0x3051, 2,      0x2        );
  fRegisterValue[kADCDAC]  .Init("ADCDAC",      0x3052, 5,      0x10       );
  fRegisterValue[kADCPAR]  .Init("ADCPAR",      0x3053, 18,     0x195EF    );
  fRegisterValue[kADCTST]  .Init("ADCTST",      0x3054, 2,      0x0        );
  fRegisterValue[kSADCAZ]  .Init("SADCAZ",      0x3055, 1,      0x1        );
  fRegisterValue[kFGF0]    .Init("FGF0",        0x3080, 9,      0x000      );
  fRegisterValue[kFGF1]    .Init("FGF1",        0x3081, 9,      0x000      );
  fRegisterValue[kFGF2]    .Init("FGF2",        0x3082, 9,      0x000      );
  fRegisterValue[kFGF3]    .Init("FGF3",        0x3083, 9,      0x000      );
  fRegisterValue[kFGF4]    .Init("FGF4",        0x3084, 9,      0x000      );
  fRegisterValue[kFGF5]    .Init("FGF5",        0x3085, 9,      0x000      );
  fRegisterValue[kFGF6]    .Init("FGF6",        0x3086, 9,      0x000      );
  fRegisterValue[kFGF7]    .Init("FGF7",        0x3087, 9,      0x000      );
  fRegisterValue[kFGF8]    .Init("FGF8",        0x3088, 9,      0x000      );
  fRegisterValue[kFGF9]    .Init("FGF9",        0x3089, 9,      0x000      );
  fRegisterValue[kFGF10]   .Init("FGF10",       0x308A, 9,      0x000      );
  fRegisterValue[kFGF11]   .Init("FGF11",       0x308B, 9,      0x000      );
  fRegisterValue[kFGF12]   .Init("FGF12",       0x308C, 9,      0x000      );
  fRegisterValue[kFGF13]   .Init("FGF13",       0x308D, 9,      0x000      );
  fRegisterValue[kFGF14]   .Init("FGF14",       0x308E, 9,      0x000      );
  fRegisterValue[kFGF15]   .Init("FGF15",       0x308F, 9,      0x000      );
  fRegisterValue[kFGF16]   .Init("FGF16",       0x3090, 9,      0x000      );
  fRegisterValue[kFGF17]   .Init("FGF17",       0x3091, 9,      0x000      );
  fRegisterValue[kFGF18]   .Init("FGF18",       0x3092, 9,      0x000      );
  fRegisterValue[kFGF19]   .Init("FGF19",       0x3093, 9,      0x000      );
  fRegisterValue[kFGF20]   .Init("FGF20",       0x3094, 9,      0x000      );
  fRegisterValue[kFGA0]    .Init("FGA0",        0x30A0, 6,      0x00       );
  fRegisterValue[kFGA1]    .Init("FGA1",        0x30A1, 6,      0x00       );
  fRegisterValue[kFGA2]    .Init("FGA2",        0x30A2, 6,      0x00       );
  fRegisterValue[kFGA3]    .Init("FGA3",        0x30A3, 6,      0x00       );
  fRegisterValue[kFGA4]    .Init("FGA4",        0x30A4, 6,      0x00       );
  fRegisterValue[kFGA5]    .Init("FGA5",        0x30A5, 6,      0x00       );
  fRegisterValue[kFGA6]    .Init("FGA6",        0x30A6, 6,      0x00       );
  fRegisterValue[kFGA7]    .Init("FGA7",        0x30A7, 6,      0x00       );
  fRegisterValue[kFGA8]    .Init("FGA8",        0x30A8, 6,      0x00       );
  fRegisterValue[kFGA9]    .Init("FGA9",        0x30A9, 6,      0x00       );
  fRegisterValue[kFGA10]   .Init("FGA10",       0x30AA, 6,      0x00       );
  fRegisterValue[kFGA11]   .Init("FGA11",       0x30AB, 6,      0x00       );
  fRegisterValue[kFGA12]   .Init("FGA12",       0x30AC, 6,      0x00       );
  fRegisterValue[kFGA13]   .Init("FGA13",       0x30AD, 6,      0x00       );
  fRegisterValue[kFGA14]   .Init("FGA14",       0x30AE, 6,      0x00       );
  fRegisterValue[kFGA15]   .Init("FGA15",       0x30AF, 6,      0x00       );
  fRegisterValue[kFGA16]   .Init("FGA16",       0x30B0, 6,      0x00       );
  fRegisterValue[kFGA17]   .Init("FGA17",       0x30B1, 6,      0x00       );
  fRegisterValue[kFGA18]   .Init("FGA18",       0x30B2, 6,      0x00       );
  fRegisterValue[kFGA19]   .Init("FGA19",       0x30B3, 6,      0x00       );
  fRegisterValue[kFGA20]   .Init("FGA20",       0x30B4, 6,      0x00       );
  fRegisterValue[kFLL00]   .Init("FLL00",       0x3100, 6,      0x00       );  // non-linearity table, 64 x 6 bits
  fRegisterValue[kFLL01]   .Init("FLL01",       0x3101, 6,      0x00       );
  fRegisterValue[kFLL02]   .Init("FLL02",       0x3102, 6,      0x00       );
  fRegisterValue[kFLL03]   .Init("FLL03",       0x3103, 6,      0x00       );
  fRegisterValue[kFLL04]   .Init("FLL04",       0x3104, 6,      0x00       );
  fRegisterValue[kFLL05]   .Init("FLL05",       0x3105, 6,      0x00       );
  fRegisterValue[kFLL06]   .Init("FLL06",       0x3106, 6,      0x00       );
  fRegisterValue[kFLL07]   .Init("FLL07",       0x3107, 6,      0x00       );
  fRegisterValue[kFLL08]   .Init("FLL08",       0x3108, 6,      0x00       );
  fRegisterValue[kFLL09]   .Init("FLL09",       0x3109, 6,      0x00       );
  fRegisterValue[kFLL0A]   .Init("FLL0A",       0x310A, 6,      0x00       );
  fRegisterValue[kFLL0B]   .Init("FLL0B",       0x310B, 6,      0x00       );
  fRegisterValue[kFLL0C]   .Init("FLL0C",       0x310C, 6,      0x00       );
  fRegisterValue[kFLL0D]   .Init("FLL0D",       0x310D, 6,      0x00       );
  fRegisterValue[kFLL0E]   .Init("FLL0E",       0x310E, 6,      0x00       );
  fRegisterValue[kFLL0F]   .Init("FLL0F",       0x310F, 6,      0x00       );
  fRegisterValue[kFLL10]   .Init("FLL10",       0x3110, 6,      0x00       );
  fRegisterValue[kFLL11]   .Init("FLL11",       0x3111, 6,      0x00       );
  fRegisterValue[kFLL12]   .Init("FLL12",       0x3112, 6,      0x00       );
  fRegisterValue[kFLL13]   .Init("FLL13",       0x3113, 6,      0x00       );
  fRegisterValue[kFLL14]   .Init("FLL14",       0x3114, 6,      0x00       );
  fRegisterValue[kFLL15]   .Init("FLL15",       0x3115, 6,      0x00       );
  fRegisterValue[kFLL16]   .Init("FLL16",       0x3116, 6,      0x00       );
  fRegisterValue[kFLL17]   .Init("FLL17",       0x3117, 6,      0x00       );
  fRegisterValue[kFLL18]   .Init("FLL18",       0x3118, 6,      0x00       );
  fRegisterValue[kFLL19]   .Init("FLL19",       0x3119, 6,      0x00       );
  fRegisterValue[kFLL1A]   .Init("FLL1A",       0x311A, 6,      0x00       );
  fRegisterValue[kFLL1B]   .Init("FLL1B",       0x311B, 6,      0x00       );
  fRegisterValue[kFLL1C]   .Init("FLL1C",       0x311C, 6,      0x00       );
  fRegisterValue[kFLL1D]   .Init("FLL1D",       0x311D, 6,      0x00       );
  fRegisterValue[kFLL1E]   .Init("FLL1E",       0x311E, 6,      0x00       );
  fRegisterValue[kFLL1F]   .Init("FLL1F",       0x311F, 6,      0x00       );
  fRegisterValue[kFLL20]   .Init("FLL20",       0x3120, 6,      0x00       );
  fRegisterValue[kFLL21]   .Init("FLL21",       0x3121, 6,      0x00       );
  fRegisterValue[kFLL22]   .Init("FLL22",       0x3122, 6,      0x00       );
  fRegisterValue[kFLL23]   .Init("FLL23",       0x3123, 6,      0x00       );
  fRegisterValue[kFLL24]   .Init("FLL24",       0x3124, 6,      0x00       );
  fRegisterValue[kFLL25]   .Init("FLL25",       0x3125, 6,      0x00       );
  fRegisterValue[kFLL26]   .Init("FLL26",       0x3126, 6,      0x00       );
  fRegisterValue[kFLL27]   .Init("FLL27",       0x3127, 6,      0x00       );
  fRegisterValue[kFLL28]   .Init("FLL28",       0x3128, 6,      0x00       );
  fRegisterValue[kFLL29]   .Init("FLL29",       0x3129, 6,      0x00       );
  fRegisterValue[kFLL2A]   .Init("FLL2A",       0x312A, 6,      0x00       );
  fRegisterValue[kFLL2B]   .Init("FLL2B",       0x312B, 6,      0x00       );
  fRegisterValue[kFLL2C]   .Init("FLL2C",       0x312C, 6,      0x00       );
  fRegisterValue[kFLL2D]   .Init("FLL2D",       0x312D, 6,      0x00       );
  fRegisterValue[kFLL2E]   .Init("FLL2E",       0x312E, 6,      0x00       );
  fRegisterValue[kFLL2F]   .Init("FLL2F",       0x312F, 6,      0x00       );
  fRegisterValue[kFLL30]   .Init("FLL30",       0x3130, 6,      0x00       );
  fRegisterValue[kFLL31]   .Init("FLL31",       0x3131, 6,      0x00       );
  fRegisterValue[kFLL32]   .Init("FLL32",       0x3132, 6,      0x00       );
  fRegisterValue[kFLL33]   .Init("FLL33",       0x3133, 6,      0x00       );
  fRegisterValue[kFLL34]   .Init("FLL34",       0x3134, 6,      0x00       );
  fRegisterValue[kFLL35]   .Init("FLL35",       0x3135, 6,      0x00       );
  fRegisterValue[kFLL36]   .Init("FLL36",       0x3136, 6,      0x00       );
  fRegisterValue[kFLL37]   .Init("FLL37",       0x3137, 6,      0x00       );
  fRegisterValue[kFLL38]   .Init("FLL38",       0x3138, 6,      0x00       );
  fRegisterValue[kFLL39]   .Init("FLL39",       0x3139, 6,      0x00       );
  fRegisterValue[kFLL3A]   .Init("FLL3A",       0x313A, 6,      0x00       );
  fRegisterValue[kFLL3B]   .Init("FLL3B",       0x313B, 6,      0x00       );
  fRegisterValue[kFLL3C]   .Init("FLL3C",       0x313C, 6,      0x00       );
  fRegisterValue[kFLL3D]   .Init("FLL3D",       0x313D, 6,      0x00       );
  fRegisterValue[kFLL3E]   .Init("FLL3E",       0x313E, 6,      0x00       );
  fRegisterValue[kFLL3F]   .Init("FLL3F",       0x313F, 6,      0x00       );
  fRegisterValue[kPASADEL] .Init("PASADEL",     0x3158, 8,      0xFF       );  // end of non-lin table
  fRegisterValue[kPASAPHA] .Init("PASAPHA",     0x3159, 6,      0x3F       );
  fRegisterValue[kPASAPRA] .Init("PASAPRA",     0x315A, 6,      0x0F       );
  fRegisterValue[kPASADAC] .Init("PASADAC",     0x315B, 8,      0x80       );
  fRegisterValue[kPASACHM] .Init("PASACHM",     0x315C, 19,     0x7FFFF    );
  fRegisterValue[kPASASTL] .Init("PASASTL",     0x315D, 8,      0xFF       );
  fRegisterValue[kPASAPR1] .Init("PASAPR1",     0x315E, 1,      0x0        );
  fRegisterValue[kPASAPR0] .Init("PASAPR0",     0x315F, 1,      0x0        );
  fRegisterValue[kSADCTRG] .Init("SADCTRG",     0x3161, 1,      0x0        );
  fRegisterValue[kSADCRUN] .Init("SADCRUN",     0x3162, 1,      0x0        );
  fRegisterValue[kSADCPWR] .Init("SADCPWR",     0x3163, 3,      0x7        );
  fRegisterValue[kL0TSIM]  .Init("L0TSIM",      0x3165, 14,     0x0050     );
  fRegisterValue[kSADCEC]  .Init("SADCEC",      0x3166, 7,      0x00       );
  fRegisterValue[kSADCMC]  .Init("SADCMC",      0x3170, 8,      0xC0       );
  fRegisterValue[kSADCOC]  .Init("SADCOC",      0x3171, 8,      0x19       );
  fRegisterValue[kSADCGTB] .Init("SADCGTB",     0x3172, 32,     0x37737700 );
  fRegisterValue[kSEBDEN]  .Init("SEBDEN",      0x3178, 3,      0x0        );
  fRegisterValue[kSEBDOU]  .Init("SEBDOU",      0x3179, 3,      0x0        );
  fRegisterValue[kTPL00]   .Init("TPL00",       0x3180, 5,      0x00       );  // pos table, 128 x 5 bits
  fRegisterValue[kTPL01]   .Init("TPL01",       0x3181, 5,      0x00       );
  fRegisterValue[kTPL02]   .Init("TPL02",       0x3182, 5,      0x00       );
  fRegisterValue[kTPL03]   .Init("TPL03",       0x3183, 5,      0x00       );
  fRegisterValue[kTPL04]   .Init("TPL04",       0x3184, 5,      0x00       );
  fRegisterValue[kTPL05]   .Init("TPL05",       0x3185, 5,      0x00       );
  fRegisterValue[kTPL06]   .Init("TPL06",       0x3186, 5,      0x00       );
  fRegisterValue[kTPL07]   .Init("TPL07",       0x3187, 5,      0x00       );
  fRegisterValue[kTPL08]   .Init("TPL08",       0x3188, 5,      0x00       );
  fRegisterValue[kTPL09]   .Init("TPL09",       0x3189, 5,      0x00       );
  fRegisterValue[kTPL0A]   .Init("TPL0A",       0x318A, 5,      0x00       );
  fRegisterValue[kTPL0B]   .Init("TPL0B",       0x318B, 5,      0x00       );
  fRegisterValue[kTPL0C]   .Init("TPL0C",       0x318C, 5,      0x00       );
  fRegisterValue[kTPL0D]   .Init("TPL0D",       0x318D, 5,      0x00       );
  fRegisterValue[kTPL0E]   .Init("TPL0E",       0x318E, 5,      0x00       );
  fRegisterValue[kTPL0F]   .Init("TPL0F",       0x318F, 5,      0x00       );
  fRegisterValue[kTPL10]   .Init("TPL10",       0x3190, 5,      0x00       );
  fRegisterValue[kTPL11]   .Init("TPL11",       0x3191, 5,      0x00       );
  fRegisterValue[kTPL12]   .Init("TPL12",       0x3192, 5,      0x00       );
  fRegisterValue[kTPL13]   .Init("TPL13",       0x3193, 5,      0x00       );
  fRegisterValue[kTPL14]   .Init("TPL14",       0x3194, 5,      0x00       );
  fRegisterValue[kTPL15]   .Init("TPL15",       0x3195, 5,      0x00       );
  fRegisterValue[kTPL16]   .Init("TPL16",       0x3196, 5,      0x00       );
  fRegisterValue[kTPL17]   .Init("TPL17",       0x3197, 5,      0x00       );
  fRegisterValue[kTPL18]   .Init("TPL18",       0x3198, 5,      0x00       );
  fRegisterValue[kTPL19]   .Init("TPL19",       0x3199, 5,      0x00       );
  fRegisterValue[kTPL1A]   .Init("TPL1A",       0x319A, 5,      0x00       );
  fRegisterValue[kTPL1B]   .Init("TPL1B",       0x319B, 5,      0x00       );
  fRegisterValue[kTPL1C]   .Init("TPL1C",       0x319C, 5,      0x00       );
  fRegisterValue[kTPL1D]   .Init("TPL1D",       0x319D, 5,      0x00       );
  fRegisterValue[kTPL1E]   .Init("TPL1E",       0x319E, 5,      0x00       );
  fRegisterValue[kTPL1F]   .Init("TPL1F",       0x319F, 5,      0x00       );
  fRegisterValue[kTPL20]   .Init("TPL20",       0x31A0, 5,      0x00       );
  fRegisterValue[kTPL21]   .Init("TPL21",       0x31A1, 5,      0x00       );
  fRegisterValue[kTPL22]   .Init("TPL22",       0x31A2, 5,      0x00       );
  fRegisterValue[kTPL23]   .Init("TPL23",       0x31A3, 5,      0x00       );
  fRegisterValue[kTPL24]   .Init("TPL24",       0x31A4, 5,      0x00       );
  fRegisterValue[kTPL25]   .Init("TPL25",       0x31A5, 5,      0x00       );
  fRegisterValue[kTPL26]   .Init("TPL26",       0x31A6, 5,      0x00       );
  fRegisterValue[kTPL27]   .Init("TPL27",       0x31A7, 5,      0x00       );
  fRegisterValue[kTPL28]   .Init("TPL28",       0x31A8, 5,      0x00       );
  fRegisterValue[kTPL29]   .Init("TPL29",       0x31A9, 5,      0x00       );
  fRegisterValue[kTPL2A]   .Init("TPL2A",       0x31AA, 5,      0x00       );
  fRegisterValue[kTPL2B]   .Init("TPL2B",       0x31AB, 5,      0x00       );
  fRegisterValue[kTPL2C]   .Init("TPL2C",       0x31AC, 5,      0x00       );
  fRegisterValue[kTPL2D]   .Init("TPL2D",       0x31AD, 5,      0x00       );
  fRegisterValue[kTPL2E]   .Init("TPL2E",       0x31AE, 5,      0x00       );
  fRegisterValue[kTPL2F]   .Init("TPL2F",       0x31AF, 5,      0x00       );
  fRegisterValue[kTPL30]   .Init("TPL30",       0x31B0, 5,      0x00       );
  fRegisterValue[kTPL31]   .Init("TPL31",       0x31B1, 5,      0x00       );
  fRegisterValue[kTPL32]   .Init("TPL32",       0x31B2, 5,      0x00       );
  fRegisterValue[kTPL33]   .Init("TPL33",       0x31B3, 5,      0x00       );
  fRegisterValue[kTPL34]   .Init("TPL34",       0x31B4, 5,      0x00       );
  fRegisterValue[kTPL35]   .Init("TPL35",       0x31B5, 5,      0x00       );
  fRegisterValue[kTPL36]   .Init("TPL36",       0x31B6, 5,      0x00       );
  fRegisterValue[kTPL37]   .Init("TPL37",       0x31B7, 5,      0x00       );
  fRegisterValue[kTPL38]   .Init("TPL38",       0x31B8, 5,      0x00       );
  fRegisterValue[kTPL39]   .Init("TPL39",       0x31B9, 5,      0x00       );
  fRegisterValue[kTPL3A]   .Init("TPL3A",       0x31BA, 5,      0x00       );
  fRegisterValue[kTPL3B]   .Init("TPL3B",       0x31BB, 5,      0x00       );
  fRegisterValue[kTPL3C]   .Init("TPL3C",       0x31BC, 5,      0x00       );
  fRegisterValue[kTPL3D]   .Init("TPL3D",       0x31BD, 5,      0x00       );
  fRegisterValue[kTPL3E]   .Init("TPL3E",       0x31BE, 5,      0x00       );
  fRegisterValue[kTPL3F]   .Init("TPL3F",       0x31BF, 5,      0x00       );
  fRegisterValue[kTPL40]   .Init("TPL40",       0x31C0, 5,      0x00       );
  fRegisterValue[kTPL41]   .Init("TPL41",       0x31C1, 5,      0x00       );
  fRegisterValue[kTPL42]   .Init("TPL42",       0x31C2, 5,      0x00       );
  fRegisterValue[kTPL43]   .Init("TPL43",       0x31C3, 5,      0x00       );
  fRegisterValue[kTPL44]   .Init("TPL44",       0x31C4, 5,      0x00       );
  fRegisterValue[kTPL45]   .Init("TPL45",       0x31C5, 5,      0x00       );
  fRegisterValue[kTPL46]   .Init("TPL46",       0x31C6, 5,      0x00       );
  fRegisterValue[kTPL47]   .Init("TPL47",       0x31C7, 5,      0x00       );
  fRegisterValue[kTPL48]   .Init("TPL48",       0x31C8, 5,      0x00       );
  fRegisterValue[kTPL49]   .Init("TPL49",       0x31C9, 5,      0x00       );
  fRegisterValue[kTPL4A]   .Init("TPL4A",       0x31CA, 5,      0x00       );
  fRegisterValue[kTPL4B]   .Init("TPL4B",       0x31CB, 5,      0x00       );
  fRegisterValue[kTPL4C]   .Init("TPL4C",       0x31CC, 5,      0x00       );
  fRegisterValue[kTPL4D]   .Init("TPL4D",       0x31CD, 5,      0x00       );
  fRegisterValue[kTPL4E]   .Init("TPL4E",       0x31CE, 5,      0x00       );
  fRegisterValue[kTPL4F]   .Init("TPL4F",       0x31CF, 5,      0x00       );
  fRegisterValue[kTPL50]   .Init("TPL50",       0x31D0, 5,      0x00       );
  fRegisterValue[kTPL51]   .Init("TPL51",       0x31D1, 5,      0x00       );
  fRegisterValue[kTPL52]   .Init("TPL52",       0x31D2, 5,      0x00       );
  fRegisterValue[kTPL53]   .Init("TPL53",       0x31D3, 5,      0x00       );
  fRegisterValue[kTPL54]   .Init("TPL54",       0x31D4, 5,      0x00       );
  fRegisterValue[kTPL55]   .Init("TPL55",       0x31D5, 5,      0x00       );
  fRegisterValue[kTPL56]   .Init("TPL56",       0x31D6, 5,      0x00       );
  fRegisterValue[kTPL57]   .Init("TPL57",       0x31D7, 5,      0x00       );
  fRegisterValue[kTPL58]   .Init("TPL58",       0x31D8, 5,      0x00       );
  fRegisterValue[kTPL59]   .Init("TPL59",       0x31D9, 5,      0x00       );
  fRegisterValue[kTPL5A]   .Init("TPL5A",       0x31DA, 5,      0x00       );
  fRegisterValue[kTPL5B]   .Init("TPL5B",       0x31DB, 5,      0x00       );
  fRegisterValue[kTPL5C]   .Init("TPL5C",       0x31DC, 5,      0x00       );
  fRegisterValue[kTPL5D]   .Init("TPL5D",       0x31DD, 5,      0x00       );
  fRegisterValue[kTPL5E]   .Init("TPL5E",       0x31DE, 5,      0x00       );
  fRegisterValue[kTPL5F]   .Init("TPL5F",       0x31DF, 5,      0x00       );
  fRegisterValue[kTPL60]   .Init("TPL60",       0x31E0, 5,      0x00       );
  fRegisterValue[kTPL61]   .Init("TPL61",       0x31E1, 5,      0x00       );
  fRegisterValue[kTPL62]   .Init("TPL62",       0x31E2, 5,      0x00       );
  fRegisterValue[kTPL63]   .Init("TPL63",       0x31E3, 5,      0x00       );
  fRegisterValue[kTPL64]   .Init("TPL64",       0x31E4, 5,      0x00       );
  fRegisterValue[kTPL65]   .Init("TPL65",       0x31E5, 5,      0x00       );
  fRegisterValue[kTPL66]   .Init("TPL66",       0x31E6, 5,      0x00       );
  fRegisterValue[kTPL67]   .Init("TPL67",       0x31E7, 5,      0x00       );
  fRegisterValue[kTPL68]   .Init("TPL68",       0x31E8, 5,      0x00       );
  fRegisterValue[kTPL69]   .Init("TPL69",       0x31E9, 5,      0x00       );
  fRegisterValue[kTPL6A]   .Init("TPL6A",       0x31EA, 5,      0x00       );
  fRegisterValue[kTPL6B]   .Init("TPL6B",       0x31EB, 5,      0x00       );
  fRegisterValue[kTPL6C]   .Init("TPL6C",       0x31EC, 5,      0x00       );
  fRegisterValue[kTPL6D]   .Init("TPL6D",       0x31ED, 5,      0x00       );
  fRegisterValue[kTPL6E]   .Init("TPL6E",       0x31EE, 5,      0x00       );
  fRegisterValue[kTPL6F]   .Init("TPL6F",       0x31EF, 5,      0x00       );
  fRegisterValue[kTPL70]   .Init("TPL70",       0x31F0, 5,      0x00       );
  fRegisterValue[kTPL71]   .Init("TPL71",       0x31F1, 5,      0x00       );
  fRegisterValue[kTPL72]   .Init("TPL72",       0x31F2, 5,      0x00       );
  fRegisterValue[kTPL73]   .Init("TPL73",       0x31F3, 5,      0x00       );
  fRegisterValue[kTPL74]   .Init("TPL74",       0x31F4, 5,      0x00       );
  fRegisterValue[kTPL75]   .Init("TPL75",       0x31F5, 5,      0x00       );
  fRegisterValue[kTPL76]   .Init("TPL76",       0x31F6, 5,      0x00       );
  fRegisterValue[kTPL77]   .Init("TPL77",       0x31F7, 5,      0x00       );
  fRegisterValue[kTPL78]   .Init("TPL78",       0x31F8, 5,      0x00       );
  fRegisterValue[kTPL79]   .Init("TPL79",       0x31F9, 5,      0x00       );
  fRegisterValue[kTPL7A]   .Init("TPL7A",       0x31FA, 5,      0x00       );
  fRegisterValue[kTPL7B]   .Init("TPL7B",       0x31FB, 5,      0x00       );
  fRegisterValue[kTPL7C]   .Init("TPL7C",       0x31FC, 5,      0x00       );
  fRegisterValue[kTPL7D]   .Init("TPL7D",       0x31FD, 5,      0x00       );
  fRegisterValue[kTPL7E]   .Init("TPL7E",       0x31FE, 5,      0x00       );
  fRegisterValue[kTPL7F]   .Init("TPL7F",       0x31FF, 5,      0x00       );
  fRegisterValue[kMEMRW]   .Init("MEMRW",       0xD000, 7,      0x79       );  // end of pos table
  fRegisterValue[kMEMCOR]  .Init("MEMCOR",      0xD001, 9,      0x000      );
  fRegisterValue[kDMDELA]  .Init("DMDELA",      0xD002, 4,      0x8        );
  fRegisterValue[kDMDELS]  .Init("DMDELS",      0xD003, 4,      0x8        );
}


void AliTRDtrapConfig::ResetRegs()
{
  // Reset the content of all TRAP registers to the reset values (see TRAP User Manual)

  for (Int_t iReg = 0; iReg < kLastReg; iReg++) {
    fRegisterValue[iReg].Reset();
  }
}


void AliTRDtrapConfig::ResetDmem()
{
  // reset the data memory

  for(Int_t iAddr = 0; iAddr < fgkDmemWords; iAddr++)
    fDmem[iAddr].Reset();
}


Int_t AliTRDtrapConfig::GetTrapReg(TrapReg_t reg, Int_t det, Int_t rob, Int_t mcm) const
{
  // get the value of an individual TRAP register
  // if it is individual for TRAPs a valid TRAP has to be specified

  if ((reg < 0) || (reg >= kLastReg)) {
    AliError("Non-existing register requested");
    return 0;
  }
  else {
    return fRegisterValue[reg].GetValue(det, rob, mcm);
  }
}


Bool_t AliTRDtrapConfig::SetTrapReg(TrapReg_t reg, Int_t value, Int_t det)
{
  // set a value for the given TRAP register on all chambers,

  return fRegisterValue[reg].SetValue(value, det);
}


Bool_t AliTRDtrapConfig::SetTrapReg(TrapReg_t reg, Int_t value, Int_t det, Int_t rob, Int_t mcm)
{
  // set the value for the given TRAP register of an individual MCM

  return fRegisterValue[reg].SetValue(value, det, rob, mcm);
}


UInt_t AliTRDtrapConfig::Peek(Int_t addr, Int_t det, Int_t rob, Int_t mcm) const
{
  // reading from given address

  if ( (addr >= fgkDmemStartAddress) &&
       (addr < (fgkDmemStartAddress + fgkDmemWords)) ) {
    return GetDmemUnsigned(addr, det, rob, mcm);
  }
  else {
    TrapReg_t mcmReg = GetRegByAddress(addr);
    if ( mcmReg >= 0 && mcmReg < kLastReg) {
      return (UInt_t) GetTrapReg(mcmReg, det, rob, mcm);
    }
  }

  AliError(Form("peek for invalid addr: 0x%04x", addr));
  return 0;
}


Bool_t AliTRDtrapConfig::Poke(Int_t addr, UInt_t value, Int_t det, Int_t rob, Int_t mcm)
{
  // writing to given address

  if ( (addr >= fgkDmemStartAddress) &&
       (addr < (fgkDmemStartAddress + fgkDmemWords)) ) {
    AliDebug(2, Form("DMEM 0x%08x : %i", addr, value));
    return SetDmem(addr, value, det, rob, mcm);
  }
  else {
    TrapReg_t mcmReg = GetRegByAddress(addr);
    if ( mcmReg >= 0 && mcmReg < kLastReg) {
      AliDebug(2, Form("Register: %s : %i\n", GetRegName(mcmReg), value));
      return SetTrapReg(mcmReg, (UInt_t) value, det, rob, mcm);
    }
  }

  AliError(Form("poke for invalid address: 0x%04x", addr));
  return kFALSE;
}


Bool_t AliTRDtrapConfig::SetDmemAlloc(Int_t addr, Alloc_t mode)
{
  addr = addr - fgkDmemStartAddress;

  if(addr < 0 || addr >=  fgkDmemWords) {
    AliError(Form("Invalid DMEM address: 0x%04x", addr+fgkDmemStartAddress));
    return kFALSE;
  }
  else {
    fDmem[addr].Allocate(mode);
    return kTRUE;
  }
}


Bool_t AliTRDtrapConfig::SetDmem(Int_t addr, UInt_t value, Int_t det)
{
  // Set the content of the given DMEM address

  addr = addr - fgkDmemStartAddress;

  if(addr < 0 || addr >=  fgkDmemWords) {
    AliError(Form("No DMEM address: 0x%08x", addr+fgkDmemStartAddress));
    return kFALSE;
  }

  if (!fDmem[addr].SetValue(value, det)) {
    AliError(Form("Problem writing to DMEM address 0x%04x", addr));
    return kFALSE;
  }
  else
    return kTRUE;
}


Bool_t AliTRDtrapConfig::SetDmem(Int_t addr, UInt_t value, Int_t det, Int_t rob, Int_t mcm)
{
  // Set the content of the given DMEM address
  addr = addr - fgkDmemStartAddress;

  if(addr < 0 || addr >=  fgkDmemWords) {
      AliError(Form("No DMEM address: 0x%08x", addr+fgkDmemStartAddress));
      return kFALSE;
   }

  if (!fDmem[addr].SetValue(value, det, rob, mcm)) {
    AliError(Form("Problem writing to DMEM address 0x%04x", addr));
    return kFALSE;
  }
  else
    return kTRUE;
}


UInt_t AliTRDtrapConfig::GetDmemUnsigned(Int_t addr, Int_t det, Int_t rob, Int_t mcm) const
{
  // get the content of the data memory at the given address
  // (only if the value is the same for all MCMs)

   addr = addr - fgkDmemStartAddress;

  if(addr < 0 || addr >=  fgkDmemWords) {
    AliError(Form("No DMEM address: 0x%08x", addr+fgkDmemStartAddress));
    return 0;
  }

  return fDmem[addr].GetValue(det, rob, mcm);
}


Bool_t AliTRDtrapConfig::PrintTrapReg(TrapReg_t reg, Int_t det, Int_t rob, Int_t mcm) const
{
  // print the value stored in the given register
  // if it is individual a valid MCM has to be specified

  if((det >= 0 && det < AliTRDgeometry::Ndet()) &&
     (rob >= 0 && rob < AliTRDfeeParam::GetNrobC1()) &&
     (mcm >= 0 && mcm < AliTRDfeeParam::GetNmcmRob() + 2)) {
    printf("%10s (%2i bits) at 0x%04x is 0x%08x and resets to: 0x%08x (currently individual mode)\n",
	   GetRegName((TrapReg_t) reg),
	   GetRegNBits((TrapReg_t) reg),
	   GetRegAddress((TrapReg_t) reg),
	   fRegisterValue[reg].GetValue(det, rob, mcm),
	   GetRegResetValue((TrapReg_t) reg));
  }
  else {
    AliError("Register value is MCM-specific: Invalid detector, ROB or MCM requested");
    return kFALSE;
  }

  return kTRUE;
}


Bool_t AliTRDtrapConfig::PrintTrapAddr(Int_t addr, Int_t det, Int_t rob, Int_t mcm) const
{
  // print the value stored at the given address in the MCM chip
  TrapReg_t reg = GetRegByAddress(addr);
  if (reg >= 0 && reg < kLastReg) {
    return PrintTrapReg(reg, det, rob, mcm);
  }
  else {
    AliError(Form("There is no register at address 0x%08x in the simulator", addr));
    return kFALSE;
  }
}


AliTRDtrapConfig::TrapReg_t AliTRDtrapConfig::GetRegByAddress(Int_t address) const
{
  // get register by its address
  // used for reading of configuration data as sent to real FEE

  if (address < fgkRegisterAddressBlockStart[0])
    return kLastReg;
  else if (address < fgkRegisterAddressBlockStart[0] + fgkRegisterAddressBlockSize[0])
    return fgRegAddressMap[address - fgkRegisterAddressBlockStart[0]];
  else if (address < fgkRegisterAddressBlockStart[1])
    return kLastReg;
  else if (address < fgkRegisterAddressBlockStart[1] + fgkRegisterAddressBlockSize[1])
    return fgRegAddressMap[address - fgkRegisterAddressBlockStart[1] + fgkRegisterAddressBlockSize[0]];
  else if (address < fgkRegisterAddressBlockStart[2])
    return kLastReg;
  else if (address < fgkRegisterAddressBlockStart[2] + fgkRegisterAddressBlockSize[2])
    return fgRegAddressMap[address - fgkRegisterAddressBlockStart[2] + fgkRegisterAddressBlockSize[1] + fgkRegisterAddressBlockSize[0]];
  else
    return kLastReg;
}


void AliTRDtrapConfig::PrintMemDatx(ostream &os, Int_t addr) const
{
  // print the content of the data memory as datx

   PrintMemDatx(os, addr, 0, 0, 127);
}

void AliTRDtrapConfig::PrintMemDatx(ostream &os, Int_t addr, Int_t det, Int_t rob, Int_t mcm) const
{
  // print the content of the data memory as datx

   if(addr < fgkDmemStartAddress || addr >= fgkDmemStartAddress+fgkDmemWords) {
      AliError(Form("Invalid DMEM address 0x%08x!", addr));
      return;
   }
   PrintDatx(os, addr, GetDmemUnsigned(addr, det, rob, mcm), rob, mcm);
}


void AliTRDtrapConfig::PrintMemDatx(ostream &os, TrapReg_t reg) const
{
  // print the content of the data memory as datx

   PrintMemDatx(os, reg, 0, 0, 127);
}


void AliTRDtrapConfig::PrintMemDatx(ostream &os, TrapReg_t reg, Int_t det, Int_t rob, Int_t mcm) const
{
  // print the content of the data memory as datx

   if(reg>= kLastReg) {
      AliError(Form("Invalid register %i!", reg));
      return;
   }
   PrintDatx(os, GetRegAddress(reg), GetTrapReg(reg, det, rob, mcm), rob, mcm);
}


void AliTRDtrapConfig::PrintDatx(ostream &os, UInt_t addr, UInt_t data, Int_t rob, Int_t mcm) const
{
  // print the value at the given address as datx

   os << std::setw(5) << 10
      << std::setw(8) << addr
      << std::setw(12) << data;
   if(mcm==127)
      os << std::setw(8) << 127;
   else
      os << std::setw(8) << AliTRDfeeParam::AliToExtAli(rob, mcm);

   os << std::endl;
}


AliTRDtrapConfig::AliTRDtrapValue::AliTRDtrapValue() :
  TObject(),
  fAllocMode(kAllocGlobal),
  fSize(1),
  fData(new UInt_t[1]),
  fValid(new Bool_t[1])
{
  fData[0] = 0;
  fValid[0] = kTRUE;
}


Bool_t AliTRDtrapConfig::AliTRDtrapValue::Allocate(Alloc_t alloc)
{
  // allocate memory for the specified granularity

  delete [] fData;
  delete [] fValid;

  fAllocMode = alloc;
  fSize = fgkSize[fAllocMode];

  if (fSize > 0) {
    fData = new UInt_t[fSize];
    fValid = new Bool_t[fSize];
    for (Int_t i = 0; i < fSize; ++i) {
      fData[i] = 0;
      fValid[i] = kFALSE;
    }
  }
  else {
    fData = 0x0;
    fValid = 0x0;
  }

  return kTRUE;
}


Int_t AliTRDtrapConfig::AliTRDtrapValue::GetIdx(Int_t det, Int_t rob, Int_t mcm) const
{
  // return Idx to access the data for the given position

  Int_t idx = -1;

  switch (fAllocMode) {
  case kAllocNone:
    idx = -1;
    break;
  case kAllocGlobal:
    idx = 0;
    break;
  case kAllocByDetector:
    idx = det;
    break;
  case kAllocByHC:
    idx = det + (rob % 2);
    break;
  case kAllocByMCM:
    idx = 18*8*det + 18*rob + mcm;
    break;
  case kAllocByLayer:
    idx = det % 6;
    break;
  case kAllocByMCMinSM:
    idx = 18*8*(det%30) + 18*rob + mcm;
    break;
  default:
    idx = -1;
    AliError("Invalid allocation mode");
  }
  if (idx < fSize) {
    return idx;
  }
  else {
    AliError(Form("Index too large %i (size %i) for %s", idx, fSize, this->GetName()));
    return  -1;
  }
}


Bool_t AliTRDtrapConfig::AliTRDtrapValue::SetData(UInt_t value)
{
  // set the given value everywhere

  for (Int_t i = 0; i < fSize; ++i) {
    fData[i] = value;
    fValid[i] = kFALSE;
  }

  return kTRUE;
}


Bool_t AliTRDtrapConfig::AliTRDtrapValue::SetData(UInt_t value, Int_t det)
{
  // set the data for a given detector

  Int_t idx = GetIdx(det, 0, 0);

  if (idx >= 0) {
    // short cut for detector-wise allocation
    if (fAllocMode == kAllocByDetector) {
      if (fValid[idx] && (fData[idx] != value)) {
	AliDebug(1, Form("Overwriting previous value %i of %s with %i for %i!",
			 fData[idx], this->GetName(), value, det));
      }
      fData[idx] = value;
      fValid[idx] = kTRUE;
      return kTRUE;
    }
    else {
      for (Int_t rob = 0; rob < 8; ++rob) {
	for (Int_t mcm = 0; mcm < 18; ++mcm) {
	  idx = GetIdx(det, rob, mcm);
	  if (fValid[idx] && (fData[idx] != value)) {
	    AliDebug(1, Form("Overwriting previous value %i of %s with %i for %i %i:%02i!",
			     fData[idx], this->GetName(), value, det, rob, mcm));
	  }
	  fData[idx] = value;
	  fValid[idx] = kTRUE;
	}
      }
      return kTRUE;
    }
  }

  if (fAllocMode == kAllocNone) {
    // assume nobody cares
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliTRDtrapConfig::AliTRDtrapValue::SetData(UInt_t value, Int_t det, Int_t rob, Int_t mcm)
{
  // set data for an individual MCM

  Int_t idx = GetIdx(det, rob, mcm);

  if (idx >= 0) {
    if (fValid[idx] && (fData[idx] != value)) {
      AliDebug(1, Form("Overwriting previous value %i of %s with %i for %i %i:%02i (idx: %i)!",
		       fData[idx], this->GetName(), value, det, rob, mcm, idx));
    }
    fData[idx] = value;
    fValid[idx] = kTRUE;
    return kTRUE;
  }
  else if (fAllocMode == kAllocNone) {
    return kTRUE;
  }
  else {
    AliError(Form("setting failed"));
    return kFALSE;
  }
}

UInt_t AliTRDtrapConfig::AliTRDtrapValue::GetData(Int_t det, Int_t rob, Int_t mcm) const
{
  // read data for the given MCM

  Int_t idx = GetIdx(det, rob, mcm);
  if (idx >= 0) {
    if (!fValid[idx])
      AliDebug(1,Form("reading from unwritten address: %s at idx %i: %i", this->GetName(), idx, fValid[idx]));
    return fData[idx];
  }
  else {
    AliError("read from invalid address");
    return 0;
  }
}

AliTRDtrapConfig::AliTRDtrapRegister::AliTRDtrapRegister() :
  AliTRDtrapValue(),
  fName("invalid"),
  fAddr(0),
  fNbits(0),
  fResetValue(0)
{
  // default constructor

}

AliTRDtrapConfig::AliTRDtrapRegister::~AliTRDtrapRegister()
{
  // destructor

}

void AliTRDtrapConfig::AliTRDtrapRegister::Init(const char* name, Int_t addr, Int_t nBits, Int_t resetValue)
{
  // init the TRAP register

  if (fAddr == 0) {
    fName = name;
    fAddr = addr;
    fNbits = nBits;
    fResetValue = resetValue;
  }
  else
    AliFatal("Re-initialising an existing TRAP register");
}

