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

ClassImp(AliTRDtrapConfig)

AliTRDtrapConfig* AliTRDtrapConfig::fgInstance = 0x0;
const Int_t AliTRDtrapConfig::fgkMaxMcm = AliTRDfeeParam::GetNmcmRob() + 2;
const Int_t AliTRDtrapConfig::fgkDmemStartAddress = 0xc000;

AliTRDtrapConfig::AliTRDtrapConfig() : 
  TObject(), fScaleQ0(0), fScaleQ1(0)
{
  // default constructor, initializing array of TRAP registers

  //                              Name          Address  Nbits   Reset Value
  fRegs[kSML0]    =   SimpleReg_t("SML0",        0x0A00, 15,     0x4050     );  // Global state machine
  fRegs[kSML1]    =   SimpleReg_t("SML1",        0x0A01, 15,     0x4200     );
  fRegs[kSML2]    =   SimpleReg_t("SML2",        0x0A02, 15,     0x4384     );
  fRegs[kSMMODE]  =   SimpleReg_t("SMMODE",      0x0A03, 16,     0xF0E2     );
  fRegs[kNITM0]   =   SimpleReg_t("NITM0",       0x0A08, 14,     0x3FFF     );
  fRegs[kNITM1]   =   SimpleReg_t("NITM1",       0x0A09, 14,     0x3FFF     );
  fRegs[kNITM2]   =   SimpleReg_t("NITM2",       0x0A0A, 14,     0x3FFF     );
  fRegs[kNIP4D]   =   SimpleReg_t("NIP4D",       0x0A0B, 7,      0x7F       );
  fRegs[kCPU0CLK] =   SimpleReg_t("CPU0CLK",     0x0A20, 5,      0x07       );
  fRegs[kCPU1CLK] =   SimpleReg_t("CPU1CLK",     0x0A22, 5,      0x07       );
  fRegs[kCPU2CLK] =   SimpleReg_t("CPU2CLK",     0x0A24, 5,      0x07       );
  fRegs[kCPU3CLK] =   SimpleReg_t("CPU3CLK",     0x0A26, 5,      0x07       );
  fRegs[kNICLK]   =   SimpleReg_t("NICLK",       0x0A28, 5,      0x07       );
  fRegs[kFILCLK]  =   SimpleReg_t("FILCLK",      0x0A2A, 5,      0x07       );
  fRegs[kPRECLK]  =   SimpleReg_t("PRECLK",      0x0A2C, 5,      0x07       );
  fRegs[kADCEN]   =   SimpleReg_t("ADCEN",       0x0A2E, 5,      0x07       );
  fRegs[kNIODE]   =   SimpleReg_t("NIODE",       0x0A30, 5,      0x07       );
  fRegs[kNIOCE]   =   SimpleReg_t("NIOCE",       0x0A32, 6,      0x21       );  // bit 5 is status bit (read-only)!
  fRegs[kNIIDE]   =   SimpleReg_t("NIIDE",       0x0A34, 5,      0x07       );
  fRegs[kNIICE]   =   SimpleReg_t("NIICE",       0x0A36, 5,      0x07       );
  fRegs[kARBTIM]  =   SimpleReg_t("ARBTIM",      0x0A3F, 4,      0x0        );  // Arbiter
  fRegs[kIA0IRQ0] =   SimpleReg_t("IA0IRQ0",     0x0B00, 12,     0x000      );  // IVT of CPU0
  fRegs[kIA0IRQ1] =   SimpleReg_t("IA0IRQ1",     0x0B01, 12,     0x000      );
  fRegs[kIA0IRQ2] =   SimpleReg_t("IA0IRQ2",     0x0B02, 12,     0x000      );
  fRegs[kIA0IRQ3] =   SimpleReg_t("IA0IRQ3",     0x0B03, 12,     0x000      );
  fRegs[kIA0IRQ4] =   SimpleReg_t("IA0IRQ4",     0x0B04, 12,     0x000      );
  fRegs[kIA0IRQ5] =   SimpleReg_t("IA0IRQ5",     0x0B05, 12,     0x000      );
  fRegs[kIA0IRQ6] =   SimpleReg_t("IA0IRQ6",     0x0B06, 12,     0x000      );
  fRegs[kIA0IRQ7] =   SimpleReg_t("IA0IRQ7",     0x0B07, 12,     0x000      );
  fRegs[kIA0IRQ8] =   SimpleReg_t("IA0IRQ8",     0x0B08, 12,     0x000      );
  fRegs[kIA0IRQ9] =   SimpleReg_t("IA0IRQ9",     0x0B09, 12,     0x000      );
  fRegs[kIA0IRQA] =   SimpleReg_t("IA0IRQA",     0x0B0A, 12,     0x000      );
  fRegs[kIA0IRQB] =   SimpleReg_t("IA0IRQB",     0x0B0B, 12,     0x000      );
  fRegs[kIA0IRQC] =   SimpleReg_t("IA0IRQC",     0x0B0C, 12,     0x000      );
  fRegs[kIRQSW0]  =   SimpleReg_t("IRQSW0",      0x0B0D, 13,     0x1FFF     );
  fRegs[kIRQHW0]  =   SimpleReg_t("IRQHW0",      0x0B0E, 13,     0x0000     );
  fRegs[kIRQHL0]  =   SimpleReg_t("IRQHL0",      0x0B0F, 13,     0x0000     );
  fRegs[kIA1IRQ0] =   SimpleReg_t("IA1IRQ0",     0x0B20, 12,     0x000      );  // IVT of CPU1
  fRegs[kIA1IRQ1] =   SimpleReg_t("IA1IRQ1",     0x0B21, 12,     0x000      );
  fRegs[kIA1IRQ2] =   SimpleReg_t("IA1IRQ2",     0x0B22, 12,     0x000      );
  fRegs[kIA1IRQ3] =   SimpleReg_t("IA1IRQ3",     0x0B23, 12,     0x000      );
  fRegs[kIA1IRQ4] =   SimpleReg_t("IA1IRQ4",     0x0B24, 12,     0x000      );
  fRegs[kIA1IRQ5] =   SimpleReg_t("IA1IRQ5",     0x0B25, 12,     0x000      );
  fRegs[kIA1IRQ6] =   SimpleReg_t("IA1IRQ6",     0x0B26, 12,     0x000      );
  fRegs[kIA1IRQ7] =   SimpleReg_t("IA1IRQ7",     0x0B27, 12,     0x000      );
  fRegs[kIA1IRQ8] =   SimpleReg_t("IA1IRQ8",     0x0B28, 12,     0x000      );
  fRegs[kIA1IRQ9] =   SimpleReg_t("IA1IRQ9",     0x0B29, 12,     0x000      );
  fRegs[kIA1IRQA] =   SimpleReg_t("IA1IRQA",     0x0B2A, 12,     0x000      );
  fRegs[kIA1IRQB] =   SimpleReg_t("IA1IRQB",     0x0B2B, 12,     0x000      );
  fRegs[kIA1IRQC] =   SimpleReg_t("IA1IRQC",     0x0B2C, 12,     0x000      );
  fRegs[kIRQSW1]  =   SimpleReg_t("IRQSW1",      0x0B2D, 13,     0x1FFF     );
  fRegs[kIRQHW1]  =   SimpleReg_t("IRQHW1",      0x0B2E, 13,     0x0000     );
  fRegs[kIRQHL1]  =   SimpleReg_t("IRQHL1",      0x0B2F, 13,     0x0000     );
  fRegs[kIA2IRQ0] =   SimpleReg_t("IA2IRQ0",     0x0B40, 12,     0x000      );  // IVT of CPU2
  fRegs[kIA2IRQ1] =   SimpleReg_t("IA2IRQ1",     0x0B41, 12,     0x000      );
  fRegs[kIA2IRQ2] =   SimpleReg_t("IA2IRQ2",     0x0B42, 12,     0x000      );
  fRegs[kIA2IRQ3] =   SimpleReg_t("IA2IRQ3",     0x0B43, 12,     0x000      );
  fRegs[kIA2IRQ4] =   SimpleReg_t("IA2IRQ4",     0x0B44, 12,     0x000      );
  fRegs[kIA2IRQ5] =   SimpleReg_t("IA2IRQ5",     0x0B45, 12,     0x000      );
  fRegs[kIA2IRQ6] =   SimpleReg_t("IA2IRQ6",     0x0B46, 12,     0x000      );
  fRegs[kIA2IRQ7] =   SimpleReg_t("IA2IRQ7",     0x0B47, 12,     0x000      );
  fRegs[kIA2IRQ8] =   SimpleReg_t("IA2IRQ8",     0x0B48, 12,     0x000      );
  fRegs[kIA2IRQ9] =   SimpleReg_t("IA2IRQ9",     0x0B49, 12,     0x000      );
  fRegs[kIA2IRQA] =   SimpleReg_t("IA2IRQA",     0x0B4A, 12,     0x000      );
  fRegs[kIA2IRQB] =   SimpleReg_t("IA2IRQB",     0x0B4B, 12,     0x000      );
  fRegs[kIA2IRQC] =   SimpleReg_t("IA2IRQC",     0x0B4C, 12,     0x000      );
  fRegs[kIRQSW2]  =   SimpleReg_t("IRQSW2",      0x0B4D, 13,     0x1FFF     );
  fRegs[kIRQHW2]  =   SimpleReg_t("IRQHW2",      0x0B4E, 13,     0x0000     );
  fRegs[kIRQHL2]  =   SimpleReg_t("IRQHL2",      0x0B4F, 13,     0x0000     );
  fRegs[kIA3IRQ0] =   SimpleReg_t("IA3IRQ0",     0x0B60, 12,     0x000      );  // IVT of CPU3
  fRegs[kIA3IRQ1] =   SimpleReg_t("IA3IRQ1",     0x0B61, 12,     0x000      );
  fRegs[kIA3IRQ2] =   SimpleReg_t("IA3IRQ2",     0x0B62, 12,     0x000      );
  fRegs[kIA3IRQ3] =   SimpleReg_t("IA3IRQ3",     0x0B63, 12,     0x000      );
  fRegs[kIA3IRQ4] =   SimpleReg_t("IA3IRQ4",     0x0B64, 12,     0x000      );
  fRegs[kIA3IRQ5] =   SimpleReg_t("IA3IRQ5",     0x0B65, 12,     0x000      );
  fRegs[kIA3IRQ6] =   SimpleReg_t("IA3IRQ6",     0x0B66, 12,     0x000      );
  fRegs[kIA3IRQ7] =   SimpleReg_t("IA3IRQ7",     0x0B67, 12,     0x000      );
  fRegs[kIA3IRQ8] =   SimpleReg_t("IA3IRQ8",     0x0B68, 12,     0x000      );
  fRegs[kIA3IRQ9] =   SimpleReg_t("IA3IRQ9",     0x0B69, 12,     0x000      );
  fRegs[kIA3IRQA] =   SimpleReg_t("IA3IRQA",     0x0B6A, 12,     0x000      );
  fRegs[kIA3IRQB] =   SimpleReg_t("IA3IRQB",     0x0B6B, 12,     0x000      );
  fRegs[kIA3IRQC] =   SimpleReg_t("IA3IRQC",     0x0B6C, 12,     0x000      );
  fRegs[kIRQSW3]  =   SimpleReg_t("IRQSW3",      0x0B6D, 13,     0x1FFF     );
  fRegs[kIRQHW3]  =   SimpleReg_t("IRQHW3",      0x0B6E, 13,     0x0000     );
  fRegs[kIRQHL3]  =   SimpleReg_t("IRQHL3",      0x0B6F, 13,     0x0000     );
  fRegs[kCTGDINI] =   SimpleReg_t("CTGDINI",     0x0B80, 32,     0x00000000 );  // Global Counter/Timer
  fRegs[kCTGCTRL] =   SimpleReg_t("CTGCTRL",     0x0B81, 12,     0xE3F      );
  fRegs[kC08CPU0] =   SimpleReg_t("C08CPU0",     0x0C00, 32,     0x00000000 );  // CPU constants
  fRegs[kC09CPU0] =   SimpleReg_t("C09CPU0",     0x0C01, 32,     0x00000000 );
  fRegs[kC10CPU0] =   SimpleReg_t("C10CPU0",     0x0C02, 32,     0x00000000 );
  fRegs[kC11CPU0] =   SimpleReg_t("C11CPU0",     0x0C03, 32,     0x00000000 );
  fRegs[kC12CPUA] =   SimpleReg_t("C12CPUA",     0x0C04, 32,     0x00000000 );
  fRegs[kC13CPUA] =   SimpleReg_t("C13CPUA",     0x0C05, 32,     0x00000000 );
  fRegs[kC14CPUA] =   SimpleReg_t("C14CPUA",     0x0C06, 32,     0x00000000 );
  fRegs[kC15CPUA] =   SimpleReg_t("C15CPUA",     0x0C07, 32,     0x00000000 );
  fRegs[kC08CPU1] =   SimpleReg_t("C08CPU1",     0x0C08, 32,     0x00000000 );
  fRegs[kC09CPU1] =   SimpleReg_t("C09CPU1",     0x0C09, 32,     0x00000000 );
  fRegs[kC10CPU1] =   SimpleReg_t("C10CPU1",     0x0C0A, 32,     0x00000000 );
  fRegs[kC11CPU1] =   SimpleReg_t("C11CPU1",     0x0C0B, 32,     0x00000000 );
  fRegs[kC08CPU2] =   SimpleReg_t("C08CPU2",     0x0C10, 32,     0x00000000 );
  fRegs[kC09CPU2] =   SimpleReg_t("C09CPU2",     0x0C11, 32,     0x00000000 );
  fRegs[kC10CPU2] =   SimpleReg_t("C10CPU2",     0x0C12, 32,     0x00000000 );
  fRegs[kC11CPU2] =   SimpleReg_t("C11CPU2",     0x0C13, 32,     0x00000000 );
  fRegs[kC08CPU3] =   SimpleReg_t("C08CPU3",     0x0C18, 32,     0x00000000 );
  fRegs[kC09CPU3] =   SimpleReg_t("C09CPU3",     0x0C19, 32,     0x00000000 );
  fRegs[kC10CPU3] =   SimpleReg_t("C10CPU3",     0x0C1A, 32,     0x00000000 );
  fRegs[kC11CPU3] =   SimpleReg_t("C11CPU3",     0x0C1B, 32,     0x00000000 );
  fRegs[kNMOD]    =   SimpleReg_t("NMOD",        0x0D40, 6,      0x08       );  // NI interface
  fRegs[kNDLY]    =   SimpleReg_t("NDLY",        0x0D41, 30,     0x24924924 );
  fRegs[kNED]     =   SimpleReg_t("NED",         0x0D42, 16,     0xA240     );
  fRegs[kNTRO]    =   SimpleReg_t("NTRO",        0x0D43, 18,     0x3FFFC    );
  fRegs[kNRRO]    =   SimpleReg_t("NRRO",        0x0D44, 18,     0x3FFFC    );
  fRegs[kNES]     =   SimpleReg_t("NES",         0x0D45, 32,     0x00000000 );
  fRegs[kNTP]     =   SimpleReg_t("NTP",         0x0D46, 32,     0x0000FFFF );
  fRegs[kNBND]    =   SimpleReg_t("NBND",        0x0D47, 16,     0x6020     );
  fRegs[kNP0]     =   SimpleReg_t("NP0",         0x0D48, 11,     0x44C      );
  fRegs[kNP1]     =   SimpleReg_t("NP1",         0x0D49, 11,     0x44C      );
  fRegs[kNP2]     =   SimpleReg_t("NP2",         0x0D4A, 11,     0x44C      );
  fRegs[kNP3]     =   SimpleReg_t("NP3",         0x0D4B, 11,     0x44C      );
  fRegs[kNCUT]    =   SimpleReg_t("NCUT",        0x0D4C, 32,     0xFFFFFFFF );
  fRegs[kTPPT0]   =   SimpleReg_t("TPPT0",       0x3000, 7,      0x01       );  // Filter and Preprocessor
  fRegs[kTPFS]    =   SimpleReg_t("TPFS",        0x3001, 7,      0x05       );
  fRegs[kTPFE]    =   SimpleReg_t("TPFE",        0x3002, 7,      0x14       );
  fRegs[kTPPGR]   =   SimpleReg_t("TPPGR",       0x3003, 7,      0x15       );
  fRegs[kTPPAE]   =   SimpleReg_t("TPPAE",       0x3004, 7,      0x1E       );
  fRegs[kTPQS0]   =   SimpleReg_t("TPQS0",       0x3005, 7,      0x00       );
  fRegs[kTPQE0]   =   SimpleReg_t("TPQE0",       0x3006, 7,      0x0A       );
  fRegs[kTPQS1]   =   SimpleReg_t("TPQS1",       0x3007, 7,      0x0B       );
  fRegs[kTPQE1]   =   SimpleReg_t("TPQE1",       0x3008, 7,      0x14       );
  fRegs[kEBD]     =   SimpleReg_t("EBD",         0x3009, 3,      0x0        );
  fRegs[kEBAQA]   =   SimpleReg_t("EBAQA",       0x300A, 7,      0x00       );
  fRegs[kEBSIA]   =   SimpleReg_t("EBSIA",       0x300B, 7,      0x20       );
  fRegs[kEBSF]    =   SimpleReg_t("EBSF",        0x300C, 1,      0x1        );
  fRegs[kEBSIM]   =   SimpleReg_t("EBSIM",       0x300D, 1,      0x1        );
  fRegs[kEBPP]    =   SimpleReg_t("EBPP",        0x300E, 1,      0x1        );
  fRegs[kEBPC]    =   SimpleReg_t("EBPC",        0x300F, 1,      0x1        );
  fRegs[kEBIS]    =   SimpleReg_t("EBIS",        0x3014, 10,     0x005      );
  fRegs[kEBIT]    =   SimpleReg_t("EBIT",        0x3015, 12,     0x028      );
  fRegs[kEBIL]    =   SimpleReg_t("EBIL",        0x3016, 8,      0xF0       );
  fRegs[kEBIN]    =   SimpleReg_t("EBIN",        0x3017, 1,      0x1        );
  fRegs[kFLBY]    =   SimpleReg_t("FLBY",        0x3018, 1,      0x0        );
  fRegs[kFPBY]    =   SimpleReg_t("FPBY",        0x3019, 1,      0x0        );
  fRegs[kFGBY]    =   SimpleReg_t("FGBY",        0x301A, 1,      0x0        );
  fRegs[kFTBY]    =   SimpleReg_t("FTBY",        0x301B, 1,      0x0        );
  fRegs[kFCBY]    =   SimpleReg_t("FCBY",        0x301C, 1,      0x0        );
  fRegs[kFPTC]    =   SimpleReg_t("FPTC",        0x3020, 2,      0x3        );
  fRegs[kFPNP]    =   SimpleReg_t("FPNP",        0x3021, 9,      0x078      );
  fRegs[kFPCL]    =   SimpleReg_t("FPCL",        0x3022, 1,      0x1        );
  fRegs[kFGTA]    =   SimpleReg_t("FGTA",        0x3028, 12,     0x014      );
  fRegs[kFGTB]    =   SimpleReg_t("FGTB",        0x3029, 12,     0x80C      );
  fRegs[kFGCL]    =   SimpleReg_t("FGCL",        0x302A, 1,      0x1        );
  fRegs[kFTAL]    =   SimpleReg_t("FTAL",        0x3030, 10,     0x0F6      );
  fRegs[kFTLL]    =   SimpleReg_t("FTLL",        0x3031, 9,      0x11D      );
  fRegs[kFTLS]    =   SimpleReg_t("FTLS",        0x3032, 9,      0x0D3      );
  fRegs[kFCW1]    =   SimpleReg_t("FCW1",        0x3038, 8,      0x1E       );
  fRegs[kFCW2]    =   SimpleReg_t("FCW2",        0x3039, 8,      0xD4       );
  fRegs[kFCW3]    =   SimpleReg_t("FCW3",        0x303A, 8,      0xE6       );
  fRegs[kFCW4]    =   SimpleReg_t("FCW4",        0x303B, 8,      0x4A       );
  fRegs[kFCW5]    =   SimpleReg_t("FCW5",        0x303C, 8,      0xEF       );
  fRegs[kTPFP]    =   SimpleReg_t("TPFP",        0x3040, 9,      0x037      );
  fRegs[kTPHT]    =   SimpleReg_t("TPHT",        0x3041, 14,     0x00A0     );
  fRegs[kTPVT]    =   SimpleReg_t("TPVT",        0x3042, 6,      0x00       );
  fRegs[kTPVBY]   =   SimpleReg_t("TPVBY",       0x3043, 1,      0x0        );
  fRegs[kTPCT]    =   SimpleReg_t("TPCT",        0x3044, 5,      0x08       );
  fRegs[kTPCL]    =   SimpleReg_t("TPCL",        0x3045, 5,      0x01       );
  fRegs[kTPCBY]   =   SimpleReg_t("TPCBY",       0x3046, 1,      0x1        );
  fRegs[kTPD]     =   SimpleReg_t("TPD",         0x3047, 4,      0xF        );
  fRegs[kTPCI0]   =   SimpleReg_t("TPCI0",       0x3048, 5,      0x00       );
  fRegs[kTPCI1]   =   SimpleReg_t("TPCI1",       0x3049, 5,      0x00       );
  fRegs[kTPCI2]   =   SimpleReg_t("TPCI2",       0x304A, 5,      0x00       );
  fRegs[kTPCI3]   =   SimpleReg_t("TPCI3",       0x304B, 5,      0x00       );
  fRegs[kADCMSK]  =   SimpleReg_t("ADCMSK",      0x3050, 21,     0x1FFFFF   );
  fRegs[kADCINB]  =   SimpleReg_t("ADCINB",      0x3051, 2,      0x2        );
  fRegs[kADCDAC]  =   SimpleReg_t("ADCDAC",      0x3052, 5,      0x10       );
  fRegs[kADCPAR]  =   SimpleReg_t("ADCPAR",      0x3053, 18,     0x195EF    );
  fRegs[kADCTST]  =   SimpleReg_t("ADCTST",      0x3054, 2,      0x0        );
  fRegs[kSADCAZ]  =   SimpleReg_t("SADCAZ",      0x3055, 1,      0x1        );
  fRegs[kFGF0]    =   SimpleReg_t("FGF0",        0x3080, 9,      0x000      );
  fRegs[kFGF1]    =   SimpleReg_t("FGF1",        0x3081, 9,      0x000      );
  fRegs[kFGF2]    =   SimpleReg_t("FGF2",        0x3082, 9,      0x000      );
  fRegs[kFGF3]    =   SimpleReg_t("FGF3",        0x3083, 9,      0x000      );
  fRegs[kFGF4]    =   SimpleReg_t("FGF4",        0x3084, 9,      0x000      );
  fRegs[kFGF5]    =   SimpleReg_t("FGF5",        0x3085, 9,      0x000      );
  fRegs[kFGF6]    =   SimpleReg_t("FGF6",        0x3086, 9,      0x000      );
  fRegs[kFGF7]    =   SimpleReg_t("FGF7",        0x3087, 9,      0x000      );
  fRegs[kFGF8]    =   SimpleReg_t("FGF8",        0x3088, 9,      0x000      );
  fRegs[kFGF9]    =   SimpleReg_t("FGF9",        0x3089, 9,      0x000      );
  fRegs[kFGF10]   =   SimpleReg_t("FGF10",       0x308A, 9,      0x000      );
  fRegs[kFGF11]   =   SimpleReg_t("FGF11",       0x308B, 9,      0x000      );
  fRegs[kFGF12]   =   SimpleReg_t("FGF12",       0x308C, 9,      0x000      );
  fRegs[kFGF13]   =   SimpleReg_t("FGF13",       0x308D, 9,      0x000      );
  fRegs[kFGF14]   =   SimpleReg_t("FGF14",       0x308E, 9,      0x000      );
  fRegs[kFGF15]   =   SimpleReg_t("FGF15",       0x308F, 9,      0x000      );
  fRegs[kFGF16]   =   SimpleReg_t("FGF16",       0x3090, 9,      0x000      );
  fRegs[kFGF17]   =   SimpleReg_t("FGF17",       0x3091, 9,      0x000      );
  fRegs[kFGF18]   =   SimpleReg_t("FGF18",       0x3092, 9,      0x000      );
  fRegs[kFGF19]   =   SimpleReg_t("FGF19",       0x3093, 9,      0x000      );
  fRegs[kFGF20]   =   SimpleReg_t("FGF20",       0x3094, 9,      0x000      );
  fRegs[kFGA0]    =   SimpleReg_t("FGA0",        0x30A0, 6,      0x00       );
  fRegs[kFGA1]    =   SimpleReg_t("FGA1",        0x30A1, 6,      0x00       );
  fRegs[kFGA2]    =   SimpleReg_t("FGA2",        0x30A2, 6,      0x00       );
  fRegs[kFGA3]    =   SimpleReg_t("FGA3",        0x30A3, 6,      0x00       );
  fRegs[kFGA4]    =   SimpleReg_t("FGA4",        0x30A4, 6,      0x00       );
  fRegs[kFGA5]    =   SimpleReg_t("FGA5",        0x30A5, 6,      0x00       );
  fRegs[kFGA6]    =   SimpleReg_t("FGA6",        0x30A6, 6,      0x00       );
  fRegs[kFGA7]    =   SimpleReg_t("FGA7",        0x30A7, 6,      0x00       );
  fRegs[kFGA8]    =   SimpleReg_t("FGA8",        0x30A8, 6,      0x00       );
  fRegs[kFGA9]    =   SimpleReg_t("FGA9",        0x30A9, 6,      0x00       );
  fRegs[kFGA10]   =   SimpleReg_t("FGA10",       0x30AA, 6,      0x00       );
  fRegs[kFGA11]   =   SimpleReg_t("FGA11",       0x30AB, 6,      0x00       );
  fRegs[kFGA12]   =   SimpleReg_t("FGA12",       0x30AC, 6,      0x00       );
  fRegs[kFGA13]   =   SimpleReg_t("FGA13",       0x30AD, 6,      0x00       );
  fRegs[kFGA14]   =   SimpleReg_t("FGA14",       0x30AE, 6,      0x00       );
  fRegs[kFGA15]   =   SimpleReg_t("FGA15",       0x30AF, 6,      0x00       );
  fRegs[kFGA16]   =   SimpleReg_t("FGA16",       0x30B0, 6,      0x00       );
  fRegs[kFGA17]   =   SimpleReg_t("FGA17",       0x30B1, 6,      0x00       );
  fRegs[kFGA18]   =   SimpleReg_t("FGA18",       0x30B2, 6,      0x00       );
  fRegs[kFGA19]   =   SimpleReg_t("FGA19",       0x30B3, 6,      0x00       );
  fRegs[kFGA20]   =   SimpleReg_t("FGA20",       0x30B4, 6,      0x00       );
  fRegs[kFLL00]   =   SimpleReg_t("FLL00",       0x3100, 6,      0x00       );  // non-linearity table, 64 x 6 bits
  fRegs[kFLL01]   =   SimpleReg_t("FLL01",       0x3101, 6,      0x00       );
  fRegs[kFLL02]   =   SimpleReg_t("FLL02",       0x3102, 6,      0x00       );
  fRegs[kFLL03]   =   SimpleReg_t("FLL03",       0x3103, 6,      0x00       );
  fRegs[kFLL04]   =   SimpleReg_t("FLL04",       0x3104, 6,      0x00       );
  fRegs[kFLL05]   =   SimpleReg_t("FLL05",       0x3105, 6,      0x00       );
  fRegs[kFLL06]   =   SimpleReg_t("FLL06",       0x3106, 6,      0x00       );
  fRegs[kFLL07]   =   SimpleReg_t("FLL07",       0x3107, 6,      0x00       );
  fRegs[kFLL08]   =   SimpleReg_t("FLL08",       0x3108, 6,      0x00       );
  fRegs[kFLL09]   =   SimpleReg_t("FLL09",       0x3109, 6,      0x00       );
  fRegs[kFLL0A]   =   SimpleReg_t("FLL0A",       0x310A, 6,      0x00       );
  fRegs[kFLL0B]   =   SimpleReg_t("FLL0B",       0x310B, 6,      0x00       );
  fRegs[kFLL0C]   =   SimpleReg_t("FLL0C",       0x310C, 6,      0x00       );
  fRegs[kFLL0D]   =   SimpleReg_t("FLL0D",       0x310D, 6,      0x00       );
  fRegs[kFLL0E]   =   SimpleReg_t("FLL0E",       0x310E, 6,      0x00       );
  fRegs[kFLL0F]   =   SimpleReg_t("FLL0F",       0x310F, 6,      0x00       );
  fRegs[kFLL10]   =   SimpleReg_t("FLL10",       0x3110, 6,      0x00       );
  fRegs[kFLL11]   =   SimpleReg_t("FLL11",       0x3111, 6,      0x00       );
  fRegs[kFLL12]   =   SimpleReg_t("FLL12",       0x3112, 6,      0x00       );
  fRegs[kFLL13]   =   SimpleReg_t("FLL13",       0x3113, 6,      0x00       );
  fRegs[kFLL14]   =   SimpleReg_t("FLL14",       0x3114, 6,      0x00       );
  fRegs[kFLL15]   =   SimpleReg_t("FLL15",       0x3115, 6,      0x00       );
  fRegs[kFLL16]   =   SimpleReg_t("FLL16",       0x3116, 6,      0x00       );
  fRegs[kFLL17]   =   SimpleReg_t("FLL17",       0x3117, 6,      0x00       );
  fRegs[kFLL18]   =   SimpleReg_t("FLL18",       0x3118, 6,      0x00       );
  fRegs[kFLL19]   =   SimpleReg_t("FLL19",       0x3119, 6,      0x00       );
  fRegs[kFLL1A]   =   SimpleReg_t("FLL1A",       0x311A, 6,      0x00       );
  fRegs[kFLL1B]   =   SimpleReg_t("FLL1B",       0x311B, 6,      0x00       );
  fRegs[kFLL1C]   =   SimpleReg_t("FLL1C",       0x311C, 6,      0x00       );
  fRegs[kFLL1D]   =   SimpleReg_t("FLL1D",       0x311D, 6,      0x00       );
  fRegs[kFLL1E]   =   SimpleReg_t("FLL1E",       0x311E, 6,      0x00       );
  fRegs[kFLL1F]   =   SimpleReg_t("FLL1F",       0x311F, 6,      0x00       );
  fRegs[kFLL20]   =   SimpleReg_t("FLL20",       0x3120, 6,      0x00       );
  fRegs[kFLL21]   =   SimpleReg_t("FLL21",       0x3121, 6,      0x00       );
  fRegs[kFLL22]   =   SimpleReg_t("FLL22",       0x3122, 6,      0x00       );
  fRegs[kFLL23]   =   SimpleReg_t("FLL23",       0x3123, 6,      0x00       );
  fRegs[kFLL24]   =   SimpleReg_t("FLL24",       0x3124, 6,      0x00       );
  fRegs[kFLL25]   =   SimpleReg_t("FLL25",       0x3125, 6,      0x00       );
  fRegs[kFLL26]   =   SimpleReg_t("FLL26",       0x3126, 6,      0x00       );
  fRegs[kFLL27]   =   SimpleReg_t("FLL27",       0x3127, 6,      0x00       );
  fRegs[kFLL28]   =   SimpleReg_t("FLL28",       0x3128, 6,      0x00       );
  fRegs[kFLL29]   =   SimpleReg_t("FLL29",       0x3129, 6,      0x00       );
  fRegs[kFLL2A]   =   SimpleReg_t("FLL2A",       0x312A, 6,      0x00       );
  fRegs[kFLL2B]   =   SimpleReg_t("FLL2B",       0x312B, 6,      0x00       );
  fRegs[kFLL2C]   =   SimpleReg_t("FLL2C",       0x312C, 6,      0x00       );
  fRegs[kFLL2D]   =   SimpleReg_t("FLL2D",       0x312D, 6,      0x00       );
  fRegs[kFLL2E]   =   SimpleReg_t("FLL2E",       0x312E, 6,      0x00       );
  fRegs[kFLL2F]   =   SimpleReg_t("FLL2F",       0x312F, 6,      0x00       );
  fRegs[kFLL30]   =   SimpleReg_t("FLL30",       0x3130, 6,      0x00       );
  fRegs[kFLL31]   =   SimpleReg_t("FLL31",       0x3131, 6,      0x00       );
  fRegs[kFLL32]   =   SimpleReg_t("FLL32",       0x3132, 6,      0x00       );
  fRegs[kFLL33]   =   SimpleReg_t("FLL33",       0x3133, 6,      0x00       );
  fRegs[kFLL34]   =   SimpleReg_t("FLL34",       0x3134, 6,      0x00       );
  fRegs[kFLL35]   =   SimpleReg_t("FLL35",       0x3135, 6,      0x00       );
  fRegs[kFLL36]   =   SimpleReg_t("FLL36",       0x3136, 6,      0x00       );
  fRegs[kFLL37]   =   SimpleReg_t("FLL37",       0x3137, 6,      0x00       );
  fRegs[kFLL38]   =   SimpleReg_t("FLL38",       0x3138, 6,      0x00       );
  fRegs[kFLL39]   =   SimpleReg_t("FLL39",       0x3139, 6,      0x00       );
  fRegs[kFLL3A]   =   SimpleReg_t("FLL3A",       0x313A, 6,      0x00       );
  fRegs[kFLL3B]   =   SimpleReg_t("FLL3B",       0x313B, 6,      0x00       );
  fRegs[kFLL3C]   =   SimpleReg_t("FLL3C",       0x313C, 6,      0x00       );
  fRegs[kFLL3D]   =   SimpleReg_t("FLL3D",       0x313D, 6,      0x00       );
  fRegs[kFLL3E]   =   SimpleReg_t("FLL3E",       0x313E, 6,      0x00       );
  fRegs[kFLL3F]   =   SimpleReg_t("FLL3F",       0x313F, 6,      0x00       );
  fRegs[kPASADEL] =   SimpleReg_t("PASADEL",     0x3158, 8,      0xFF       );  // end of non-lin table
  fRegs[kPASAPHA] =   SimpleReg_t("PASAPHA",     0x3159, 6,      0x3F       );
  fRegs[kPASAPRA] =   SimpleReg_t("PASAPRA",     0x315A, 6,      0x0F       );
  fRegs[kPASADAC] =   SimpleReg_t("PASADAC",     0x315B, 8,      0x80       );
  fRegs[kPASACHM] =   SimpleReg_t("PASACHM",     0x315C, 19,     0x7FFFF    );
  fRegs[kPASASTL] =   SimpleReg_t("PASASTL",     0x315D, 8,      0xFF       );
  fRegs[kPASAPR1] =   SimpleReg_t("PASAPR1",     0x315E, 1,      0x0        );
  fRegs[kPASAPR0] =   SimpleReg_t("PASAPR0",     0x315F, 1,      0x0        );
  fRegs[kSADCTRG] =   SimpleReg_t("SADCTRG",     0x3161, 1,      0x0        );
  fRegs[kSADCRUN] =   SimpleReg_t("SADCRUN",     0x3162, 1,      0x0        );
  fRegs[kSADCPWR] =   SimpleReg_t("SADCPWR",     0x3163, 3,      0x7        );
  fRegs[kL0TSIM]  =   SimpleReg_t("L0TSIM",      0x3165, 14,     0x0050     );
  fRegs[kSADCEC]  =   SimpleReg_t("SADCEC",      0x3166, 7,      0x00       );
  fRegs[kSADCMC]  =   SimpleReg_t("SADCMC",      0x3170, 8,      0xC0       );
  fRegs[kSADCOC]  =   SimpleReg_t("SADCOC",      0x3171, 8,      0x19       );
  fRegs[kSADCGTB] =   SimpleReg_t("SADCGTB",     0x3172, 32,     0x37737700 );
  fRegs[kSEBDEN]  =   SimpleReg_t("SEBDEN",      0x3178, 3,      0x0        );
  fRegs[kSEBDOU]  =   SimpleReg_t("SEBDOU",      0x3179, 3,      0x0        );
  fRegs[kTPL00]   =   SimpleReg_t("TPL00",       0x3180, 5,      0x00       );  // pos table, 128 x 5 bits
  fRegs[kTPL01]   =   SimpleReg_t("TPL01",       0x3181, 5,      0x00       );
  fRegs[kTPL02]   =   SimpleReg_t("TPL02",       0x3182, 5,      0x00       );
  fRegs[kTPL03]   =   SimpleReg_t("TPL03",       0x3183, 5,      0x00       );
  fRegs[kTPL04]   =   SimpleReg_t("TPL04",       0x3184, 5,      0x00       );
  fRegs[kTPL05]   =   SimpleReg_t("TPL05",       0x3185, 5,      0x00       );
  fRegs[kTPL06]   =   SimpleReg_t("TPL06",       0x3186, 5,      0x00       );
  fRegs[kTPL07]   =   SimpleReg_t("TPL07",       0x3187, 5,      0x00       );
  fRegs[kTPL08]   =   SimpleReg_t("TPL08",       0x3188, 5,      0x00       );
  fRegs[kTPL09]   =   SimpleReg_t("TPL09",       0x3189, 5,      0x00       );
  fRegs[kTPL0A]   =   SimpleReg_t("TPL0A",       0x318A, 5,      0x00       );
  fRegs[kTPL0B]   =   SimpleReg_t("TPL0B",       0x318B, 5,      0x00       );
  fRegs[kTPL0C]   =   SimpleReg_t("TPL0C",       0x318C, 5,      0x00       );
  fRegs[kTPL0D]   =   SimpleReg_t("TPL0D",       0x318D, 5,      0x00       );
  fRegs[kTPL0E]   =   SimpleReg_t("TPL0E",       0x318E, 5,      0x00       );
  fRegs[kTPL0F]   =   SimpleReg_t("TPL0F",       0x318F, 5,      0x00       );
  fRegs[kTPL10]   =   SimpleReg_t("TPL10",       0x3190, 5,      0x00       );
  fRegs[kTPL11]   =   SimpleReg_t("TPL11",       0x3191, 5,      0x00       );
  fRegs[kTPL12]   =   SimpleReg_t("TPL12",       0x3192, 5,      0x00       );
  fRegs[kTPL13]   =   SimpleReg_t("TPL13",       0x3193, 5,      0x00       );
  fRegs[kTPL14]   =   SimpleReg_t("TPL14",       0x3194, 5,      0x00       );
  fRegs[kTPL15]   =   SimpleReg_t("TPL15",       0x3195, 5,      0x00       );
  fRegs[kTPL16]   =   SimpleReg_t("TPL16",       0x3196, 5,      0x00       );
  fRegs[kTPL17]   =   SimpleReg_t("TPL17",       0x3197, 5,      0x00       );
  fRegs[kTPL18]   =   SimpleReg_t("TPL18",       0x3198, 5,      0x00       );
  fRegs[kTPL19]   =   SimpleReg_t("TPL19",       0x3199, 5,      0x00       );
  fRegs[kTPL1A]   =   SimpleReg_t("TPL1A",       0x319A, 5,      0x00       );
  fRegs[kTPL1B]   =   SimpleReg_t("TPL1B",       0x319B, 5,      0x00       );
  fRegs[kTPL1C]   =   SimpleReg_t("TPL1C",       0x319C, 5,      0x00       );
  fRegs[kTPL1D]   =   SimpleReg_t("TPL1D",       0x319D, 5,      0x00       );
  fRegs[kTPL1E]   =   SimpleReg_t("TPL1E",       0x319E, 5,      0x00       );
  fRegs[kTPL1F]   =   SimpleReg_t("TPL1F",       0x319F, 5,      0x00       );
  fRegs[kTPL20]   =   SimpleReg_t("TPL20",       0x31A0, 5,      0x00       );
  fRegs[kTPL21]   =   SimpleReg_t("TPL21",       0x31A1, 5,      0x00       );
  fRegs[kTPL22]   =   SimpleReg_t("TPL22",       0x31A2, 5,      0x00       );
  fRegs[kTPL23]   =   SimpleReg_t("TPL23",       0x31A3, 5,      0x00       );
  fRegs[kTPL24]   =   SimpleReg_t("TPL24",       0x31A4, 5,      0x00       );
  fRegs[kTPL25]   =   SimpleReg_t("TPL25",       0x31A5, 5,      0x00       );
  fRegs[kTPL26]   =   SimpleReg_t("TPL26",       0x31A6, 5,      0x00       );
  fRegs[kTPL27]   =   SimpleReg_t("TPL27",       0x31A7, 5,      0x00       );
  fRegs[kTPL28]   =   SimpleReg_t("TPL28",       0x31A8, 5,      0x00       );
  fRegs[kTPL29]   =   SimpleReg_t("TPL29",       0x31A9, 5,      0x00       );
  fRegs[kTPL2A]   =   SimpleReg_t("TPL2A",       0x31AA, 5,      0x00       );
  fRegs[kTPL2B]   =   SimpleReg_t("TPL2B",       0x31AB, 5,      0x00       );
  fRegs[kTPL2C]   =   SimpleReg_t("TPL2C",       0x31AC, 5,      0x00       );
  fRegs[kTPL2D]   =   SimpleReg_t("TPL2D",       0x31AD, 5,      0x00       );
  fRegs[kTPL2E]   =   SimpleReg_t("TPL2E",       0x31AE, 5,      0x00       );
  fRegs[kTPL2F]   =   SimpleReg_t("TPL2F",       0x31AF, 5,      0x00       );
  fRegs[kTPL30]   =   SimpleReg_t("TPL30",       0x31B0, 5,      0x00       );
  fRegs[kTPL31]   =   SimpleReg_t("TPL31",       0x31B1, 5,      0x00       );
  fRegs[kTPL32]   =   SimpleReg_t("TPL32",       0x31B2, 5,      0x00       );
  fRegs[kTPL33]   =   SimpleReg_t("TPL33",       0x31B3, 5,      0x00       );
  fRegs[kTPL34]   =   SimpleReg_t("TPL34",       0x31B4, 5,      0x00       );
  fRegs[kTPL35]   =   SimpleReg_t("TPL35",       0x31B5, 5,      0x00       );
  fRegs[kTPL36]   =   SimpleReg_t("TPL36",       0x31B6, 5,      0x00       );
  fRegs[kTPL37]   =   SimpleReg_t("TPL37",       0x31B7, 5,      0x00       );
  fRegs[kTPL38]   =   SimpleReg_t("TPL38",       0x31B8, 5,      0x00       );
  fRegs[kTPL39]   =   SimpleReg_t("TPL39",       0x31B9, 5,      0x00       );
  fRegs[kTPL3A]   =   SimpleReg_t("TPL3A",       0x31BA, 5,      0x00       );
  fRegs[kTPL3B]   =   SimpleReg_t("TPL3B",       0x31BB, 5,      0x00       );
  fRegs[kTPL3C]   =   SimpleReg_t("TPL3C",       0x31BC, 5,      0x00       );
  fRegs[kTPL3D]   =   SimpleReg_t("TPL3D",       0x31BD, 5,      0x00       );
  fRegs[kTPL3E]   =   SimpleReg_t("TPL3E",       0x31BE, 5,      0x00       );
  fRegs[kTPL3F]   =   SimpleReg_t("TPL3F",       0x31BF, 5,      0x00       );
  fRegs[kTPL40]   =   SimpleReg_t("TPL40",       0x31C0, 5,      0x00       );
  fRegs[kTPL41]   =   SimpleReg_t("TPL41",       0x31C1, 5,      0x00       );
  fRegs[kTPL42]   =   SimpleReg_t("TPL42",       0x31C2, 5,      0x00       );
  fRegs[kTPL43]   =   SimpleReg_t("TPL43",       0x31C3, 5,      0x00       );
  fRegs[kTPL44]   =   SimpleReg_t("TPL44",       0x31C4, 5,      0x00       );
  fRegs[kTPL45]   =   SimpleReg_t("TPL45",       0x31C5, 5,      0x00       );
  fRegs[kTPL46]   =   SimpleReg_t("TPL46",       0x31C6, 5,      0x00       );
  fRegs[kTPL47]   =   SimpleReg_t("TPL47",       0x31C7, 5,      0x00       );
  fRegs[kTPL48]   =   SimpleReg_t("TPL48",       0x31C8, 5,      0x00       );
  fRegs[kTPL49]   =   SimpleReg_t("TPL49",       0x31C9, 5,      0x00       );
  fRegs[kTPL4A]   =   SimpleReg_t("TPL4A",       0x31CA, 5,      0x00       );
  fRegs[kTPL4B]   =   SimpleReg_t("TPL4B",       0x31CB, 5,      0x00       );
  fRegs[kTPL4C]   =   SimpleReg_t("TPL4C",       0x31CC, 5,      0x00       );
  fRegs[kTPL4D]   =   SimpleReg_t("TPL4D",       0x31CD, 5,      0x00       );
  fRegs[kTPL4E]   =   SimpleReg_t("TPL4E",       0x31CE, 5,      0x00       );
  fRegs[kTPL4F]   =   SimpleReg_t("TPL4F",       0x31CF, 5,      0x00       );
  fRegs[kTPL50]   =   SimpleReg_t("TPL50",       0x31D0, 5,      0x00       );
  fRegs[kTPL51]   =   SimpleReg_t("TPL51",       0x31D1, 5,      0x00       );
  fRegs[kTPL52]   =   SimpleReg_t("TPL52",       0x31D2, 5,      0x00       );
  fRegs[kTPL53]   =   SimpleReg_t("TPL53",       0x31D3, 5,      0x00       );
  fRegs[kTPL54]   =   SimpleReg_t("TPL54",       0x31D4, 5,      0x00       );
  fRegs[kTPL55]   =   SimpleReg_t("TPL55",       0x31D5, 5,      0x00       );
  fRegs[kTPL56]   =   SimpleReg_t("TPL56",       0x31D6, 5,      0x00       );
  fRegs[kTPL57]   =   SimpleReg_t("TPL57",       0x31D7, 5,      0x00       );
  fRegs[kTPL58]   =   SimpleReg_t("TPL58",       0x31D8, 5,      0x00       );
  fRegs[kTPL59]   =   SimpleReg_t("TPL59",       0x31D9, 5,      0x00       );
  fRegs[kTPL5A]   =   SimpleReg_t("TPL5A",       0x31DA, 5,      0x00       );
  fRegs[kTPL5B]   =   SimpleReg_t("TPL5B",       0x31DB, 5,      0x00       );
  fRegs[kTPL5C]   =   SimpleReg_t("TPL5C",       0x31DC, 5,      0x00       );
  fRegs[kTPL5D]   =   SimpleReg_t("TPL5D",       0x31DD, 5,      0x00       );
  fRegs[kTPL5E]   =   SimpleReg_t("TPL5E",       0x31DE, 5,      0x00       );
  fRegs[kTPL5F]   =   SimpleReg_t("TPL5F",       0x31DF, 5,      0x00       );
  fRegs[kTPL60]   =   SimpleReg_t("TPL60",       0x31E0, 5,      0x00       );
  fRegs[kTPL61]   =   SimpleReg_t("TPL61",       0x31E1, 5,      0x00       );
  fRegs[kTPL62]   =   SimpleReg_t("TPL62",       0x31E2, 5,      0x00       );
  fRegs[kTPL63]   =   SimpleReg_t("TPL63",       0x31E3, 5,      0x00       );
  fRegs[kTPL64]   =   SimpleReg_t("TPL64",       0x31E4, 5,      0x00       );
  fRegs[kTPL65]   =   SimpleReg_t("TPL65",       0x31E5, 5,      0x00       );
  fRegs[kTPL66]   =   SimpleReg_t("TPL66",       0x31E6, 5,      0x00       );
  fRegs[kTPL67]   =   SimpleReg_t("TPL67",       0x31E7, 5,      0x00       );
  fRegs[kTPL68]   =   SimpleReg_t("TPL68",       0x31E8, 5,      0x00       );
  fRegs[kTPL69]   =   SimpleReg_t("TPL69",       0x31E9, 5,      0x00       );
  fRegs[kTPL6A]   =   SimpleReg_t("TPL6A",       0x31EA, 5,      0x00       );
  fRegs[kTPL6B]   =   SimpleReg_t("TPL6B",       0x31EB, 5,      0x00       );
  fRegs[kTPL6C]   =   SimpleReg_t("TPL6C",       0x31EC, 5,      0x00       );
  fRegs[kTPL6D]   =   SimpleReg_t("TPL6D",       0x31ED, 5,      0x00       );
  fRegs[kTPL6E]   =   SimpleReg_t("TPL6E",       0x31EE, 5,      0x00       );
  fRegs[kTPL6F]   =   SimpleReg_t("TPL6F",       0x31EF, 5,      0x00       );
  fRegs[kTPL70]   =   SimpleReg_t("TPL70",       0x31F0, 5,      0x00       );
  fRegs[kTPL71]   =   SimpleReg_t("TPL71",       0x31F1, 5,      0x00       );
  fRegs[kTPL72]   =   SimpleReg_t("TPL72",       0x31F2, 5,      0x00       );
  fRegs[kTPL73]   =   SimpleReg_t("TPL73",       0x31F3, 5,      0x00       );
  fRegs[kTPL74]   =   SimpleReg_t("TPL74",       0x31F4, 5,      0x00       );
  fRegs[kTPL75]   =   SimpleReg_t("TPL75",       0x31F5, 5,      0x00       );
  fRegs[kTPL76]   =   SimpleReg_t("TPL76",       0x31F6, 5,      0x00       );
  fRegs[kTPL77]   =   SimpleReg_t("TPL77",       0x31F7, 5,      0x00       );
  fRegs[kTPL78]   =   SimpleReg_t("TPL78",       0x31F8, 5,      0x00       );
  fRegs[kTPL79]   =   SimpleReg_t("TPL79",       0x31F9, 5,      0x00       );
  fRegs[kTPL7A]   =   SimpleReg_t("TPL7A",       0x31FA, 5,      0x00       );
  fRegs[kTPL7B]   =   SimpleReg_t("TPL7B",       0x31FB, 5,      0x00       );
  fRegs[kTPL7C]   =   SimpleReg_t("TPL7C",       0x31FC, 5,      0x00       );
  fRegs[kTPL7D]   =   SimpleReg_t("TPL7D",       0x31FD, 5,      0x00       );
  fRegs[kTPL7E]   =   SimpleReg_t("TPL7E",       0x31FE, 5,      0x00       );
  fRegs[kTPL7F]   =   SimpleReg_t("TPL7F",       0x31FF, 5,      0x00       );
  fRegs[kMEMRW]   =   SimpleReg_t("MEMRW",       0xD000, 7,      0x79       );  // end of pos table
  fRegs[kMEMCOR]  =   SimpleReg_t("MEMCOR",      0xD001, 9,      0x000      );
  fRegs[kDMDELA]  =   SimpleReg_t("DMDELA",      0xD002, 4,      0x8        );
  fRegs[kDMDELS]  =   SimpleReg_t("DMDELS",      0xD003, 4,      0x8        );



  for(Int_t iAddr = 0; iAddr < fgkDmemWords; iAddr++) {
     
     if(iAddr == fgkDmemAddrDeflCorr - fgkDmemStartAddress) {
	fDmem[iAddr] = new UInt_t[fgkDmemSizeSmIndividual];
	fDmemDepth[iAddr] = fgkDmemSizeSmIndividual;
     }

     else if(iAddr == fgkDmemAddrNdrift - fgkDmemStartAddress) {
	fDmem[iAddr] = new UInt_t[fgkDmemSizeSmRocIndividual];
	fDmemDepth[iAddr] = fgkDmemSizeSmRocIndividual;
     }

     else if(iAddr >= fgkDmemAddrDeflCutStart-fgkDmemStartAddress && iAddr <= fgkDmemAddrDeflCutEnd-fgkDmemStartAddress) {
	fDmem[iAddr]  = new UInt_t[fgkDmemSizeSmIndividual];   
	fDmemDepth[iAddr] = fgkDmemSizeSmIndividual;
     }

     else if(iAddr >= fgkDmemAddrTrackletStart-fgkDmemStartAddress && iAddr <= fgkDmemAddrTrackletEnd-fgkDmemStartAddress) {
	fDmem[iAddr]  = new UInt_t[fgkDmemSizeTotalIndividual];
	fDmemDepth[iAddr] = fgkDmemSizeTotalIndividual;
     }

     else if(iAddr >= fgkDmemAddrLUTStart-fgkDmemStartAddress && iAddr <= fgkDmemAddrLUTEnd-fgkDmemStartAddress) {
	fDmem[iAddr]  = new UInt_t;   // same value for all MCMs
 	fDmemDepth[iAddr] = fgkDmemSizeUniform;
    }

     else if(iAddr == fgkDmemAddrLUTcor0-fgkDmemStartAddress) {
	fDmem[iAddr]  = new UInt_t[fgkDmemSizeSmIndividual];   
	fDmemDepth[iAddr]  = fgkDmemSizeSmIndividual;
     }

     else if(iAddr == fgkDmemAddrLUTcor1-fgkDmemStartAddress) {
	fDmem[iAddr]  = new UInt_t[fgkDmemSizeSmIndividual];   
	fDmemDepth[iAddr]  = fgkDmemSizeSmIndividual;
     }
	
     else if(iAddr == fgkDmemAddrLUTnbins-fgkDmemStartAddress) {
	fDmem[iAddr]  = new UInt_t;   // same value for all MCMs
	fDmemDepth[iAddr] = fgkDmemSizeUniform;
     }

     else if(iAddr == fgkDmemAddrLUTLength-fgkDmemStartAddress) {
	fDmem[iAddr]  = new UInt_t;   // same value for all MCMs
	fDmemDepth[iAddr] = fgkDmemSizeUniform;
     }

     else {
	fDmem[iAddr] = NULL;
	fDmemDepth[iAddr] = fgkDmemSizeEmpty;
     }
     
  }

  InitRegs();
  ResetDmem();
}


AliTRDtrapConfig* AliTRDtrapConfig::Instance()
{
  // return a pointer to an instance of this class

  if (!fgInstance) {
    fgInstance = new AliTRDtrapConfig();
    fgInstance->LoadConfig();
  }

  return fgInstance;
}


AliTRDtrapConfig::~AliTRDtrapConfig()
{
  for(Int_t iAddr = 0; iAddr < fgkDmemWords; iAddr++) {
     if(iAddr == fgkDmemAddrDeflCorr - fgkDmemStartAddress)
	delete [] fDmem[iAddr];

     else if(iAddr == fgkDmemAddrNdrift - fgkDmemStartAddress)
	delete [] fDmem[iAddr];

     else if(iAddr >= fgkDmemAddrDeflCutStart-fgkDmemStartAddress && iAddr <= fgkDmemAddrDeflCutEnd-fgkDmemStartAddress)
	delete [] fDmem[iAddr];

     else if(iAddr >= fgkDmemAddrTrackletStart-fgkDmemStartAddress && iAddr <= fgkDmemAddrTrackletEnd-fgkDmemStartAddress)
	delete [] fDmem[iAddr];

     else if(iAddr >= fgkDmemAddrLUTStart-fgkDmemStartAddress && iAddr <= fgkDmemAddrLUTEnd-fgkDmemStartAddress)
	delete fDmem[iAddr];

     else if(iAddr == fgkDmemAddrLUTcor0-fgkDmemStartAddress)
	delete [] fDmem[iAddr];

     else if(iAddr == fgkDmemAddrLUTcor1-fgkDmemStartAddress)
	delete [] fDmem[iAddr];
	
     else if(iAddr == fgkDmemAddrLUTnbins-fgkDmemStartAddress)
	delete fDmem[iAddr];

     else if(iAddr == fgkDmemAddrLUTLength-fgkDmemStartAddress)
	delete fDmem[iAddr];
  }
}


void AliTRDtrapConfig::InitRegs()
{
   // Reset the content of all TRAP registers to the reset values (see TRAP User Manual)

   for (Int_t iReg = 0; iReg < kLastReg; iReg++) {

     fRegisterValue[iReg].individualValue = 0x0;

     fRegisterValue[iReg].globalValue = GetRegResetValue((TrapReg_t) iReg);
     fRegisterValue[iReg].state = RegValue_t::kGlobal;
   }
}


void AliTRDtrapConfig::ResetRegs()
{
   // Reset the content of all TRAP registers to the reset values (see TRAP User Manual)

   for (Int_t iReg = 0; iReg < kLastReg; iReg++) {
      if(fRegisterValue[iReg].state == RegValue_t::kIndividual) {
	if (fRegisterValue[iReg].individualValue) {
	  delete [] fRegisterValue[iReg].individualValue;
	  fRegisterValue[iReg].individualValue = 0x0;
	}
      }

      fRegisterValue[iReg].globalValue = GetRegResetValue((TrapReg_t) iReg);
      fRegisterValue[iReg].state = RegValue_t::kGlobal;
      //    printf("%-8s: 0x%08x\n", GetRegName((TrapReg_t) iReg), fRegisterValue[iReg].globalValue);
   }
}


void AliTRDtrapConfig::ResetDmem()
{
     for(Int_t iAddr = 0; iAddr < fgkDmemWords; iAddr++) {
	if(fDmemDepth[iAddr] == 0)
	   continue;
	for(Int_t j=0; j < fDmemDepth[iAddr]; j++) {
	   fDmem[iAddr][j]=0;
	}
     }
}


Int_t AliTRDtrapConfig::GetTrapReg(TrapReg_t reg, Int_t det, Int_t rob, Int_t mcm) const
{
  // get the value of an individual TRAP register 
  // if it is individual for TRAPs a valid TRAP has to be specified

  if ((reg < 0) || (reg >= kLastReg)) {
    AliError("Non-existing register requested");
    return -1;
  }
  else {
    if (fRegisterValue[reg].state == RegValue_t::kGlobal) {
      return fRegisterValue[reg].globalValue;
    }
    else if (fRegisterValue[reg].state == RegValue_t::kIndividual) {
       if((det >= 0 && det < AliTRDgeometry::Ndet()) && 
          (rob >= 0 && rob < AliTRDfeeParam::GetNrobC1()) && 
          (mcm >= 0 && mcm < fgkMaxMcm)) {
         return fRegisterValue[reg].individualValue[det*AliTRDfeeParam::GetNrobC1()*fgkMaxMcm + rob*fgkMaxMcm + mcm];
       }
       else {
         AliError("Invalid MCM specified or register is individual");
         return -1;
       }
    }
    else {  // should never be reached
      AliError("MCM register status neither kGlobal nor kIndividual");
      return -1;
    }
  }
  return -1;
}


Bool_t AliTRDtrapConfig::SetTrapReg(TrapReg_t reg, Int_t value)
{
  // set a global value for the given TRAP register,
  // i.e. the same value for all TRAPs

   if (fRegisterValue[reg].state == RegValue_t::kGlobal) {
      fRegisterValue[reg].globalValue = value;
      return kTRUE;
   }
   else {
      AliError("Register has individual values");
   }
   return kFALSE;
}


Bool_t AliTRDtrapConfig::SetTrapReg(TrapReg_t reg, Int_t value, Int_t det)
{
  // set a global value for the given TRAP register,
  // i.e. the same value for all TRAPs

   if (fRegisterValue[reg].state == RegValue_t::kGlobal) {
      fRegisterValue[reg].globalValue = value;
      return kTRUE;
   }
   else if (fRegisterValue[reg].state == RegValue_t::kIndividual) {
      // if the register is in idividual mode but a broadcast is requested, the selected register is 
      // set to value for all MCMs on the chamber

      if( (det>=0 && det<AliTRDgeometry::Ndet())) {
	 for(Int_t rob=0; rob<AliTRDfeeParam::GetNrobC1(); rob++) {
	    for(Int_t mcm=0; mcm<fgkMaxMcm; mcm++)
	       fRegisterValue[reg].individualValue[det*AliTRDfeeParam::GetNrobC1()*fgkMaxMcm + rob*fgkMaxMcm + mcm] = value;
	 }
      }
      else {
	 AliError(Form("Invalid detector number: %i\n", det));
	 return kFALSE;
      }
   }
   else {  // should never be reached
      AliError("MCM register status neither kGlobal nor kIndividual");
      return kFALSE;
   }
   
   return kFALSE;
}


Bool_t AliTRDtrapConfig::SetTrapReg(TrapReg_t reg, Int_t value, Int_t det, Int_t rob, Int_t mcm)
{
  // set the value for the given TRAP register of an individual MCM 

  //std::cout << "-- reg: 0x" << std::hex << fRegs[reg].addr << 
  //std::dec << ", data " << value << ", det " << det << ", rob " << rob << ", mcm " << mcm << std::endl;

   if( (det >= 0 && det < AliTRDgeometry::Ndet()) && 
       (rob >= 0 && rob < AliTRDfeeParam::GetNrobC1()) && 
       (mcm >= 0 && mcm < fgkMaxMcm) ) {
     if (fRegisterValue[reg].state == RegValue_t::kGlobal) {
	Int_t defaultValue = fRegisterValue[reg].globalValue;
	
	fRegisterValue[reg].state = RegValue_t::kIndividual;
	fRegisterValue[reg].individualValue = new Int_t[AliTRDgeometry::Ndet()*AliTRDfeeParam::GetNrobC1()*fgkMaxMcm];

	for(Int_t i = 0; i < AliTRDgeometry::Ndet()*AliTRDfeeParam::GetNrobC1()*fgkMaxMcm; i++)
	   fRegisterValue[reg].individualValue[i] = defaultValue; // set the requested register of all MCMs to the value previously stored

	fRegisterValue[reg].individualValue[det*AliTRDfeeParam::GetNrobC1()*fgkMaxMcm + rob*fgkMaxMcm + mcm] = value;
     }
     else if (fRegisterValue[reg].state == RegValue_t::kIndividual) {
	fRegisterValue[reg].individualValue[det*AliTRDfeeParam::GetNrobC1()*fgkMaxMcm + rob*fgkMaxMcm + mcm] = value;
     }
     else {  // should never be reached
	AliError("MCM register status neither kGlobal nor kIndividual");
	return kFALSE;
     }
  }
   else {
      AliError(Form("Invalid value for det, ROB or MCM selected: %i, %i, %i", det, rob, mcm));
     return kFALSE;
   }

  return kTRUE;
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

  return -1;
}


Bool_t AliTRDtrapConfig::Poke(Int_t addr, UInt_t value, Int_t det, Int_t rob, Int_t mcm)
{
  // writing to given address

  if ( (addr >= fgkDmemStartAddress) && 
       (addr < (fgkDmemStartAddress + fgkDmemWords)) ) {
    AliDebug(2, Form("DMEM 0x%08x : %i", addr, value));
    SetDmem(addr, value, det, rob, mcm);
    return kTRUE;
  }
  else {
    TrapReg_t mcmReg = GetRegByAddress(addr);
    if ( mcmReg >= 0 && mcmReg < kLastReg) {
      AliDebug(2, Form("Register: %s : %i\n", GetRegName(mcmReg), value));
      SetTrapReg(mcmReg, (UInt_t) value, det, rob, mcm);
      return kTRUE;
    }
  }
  
  return kFALSE;
}


Bool_t AliTRDtrapConfig::SetDmem(Int_t addr, UInt_t value)
{
  // Set the content of the given DMEM address 

   addr = addr - fgkDmemStartAddress;

   if(addr < 0 || addr >=  fgkDmemWords) {
      AliDebug(5, Form("No DMEM address: 0x%08x", addr+fgkDmemStartAddress));
      return kFALSE;
   }

   switch(fDmemDepth[addr]) {
   case fgkDmemSizeEmpty:
      AliDebug(5, Form("DMEM address %i not active", addr));
      return kFALSE;
      break;
   case fgkDmemSizeUniform:
      if(fDmem[addr][0]!=0)
	 AliDebug(5, Form("Warning: Setting new value to DMEM 0x%08x", addr+fgkDmemStartAddress));

      fDmem[addr][0] = value;
      break;
   case fgkDmemSizeSmIndividual:
      for(Int_t i=0; i<fgkDmemSizeSmIndividual; i++) {
	 if(fDmem[addr][i]!=0)
	    AliDebug(5, Form("Warning: Setting new value to DMEM 0x%08x", addr+fgkDmemStartAddress));

	 fDmem[addr][i]=value;
      }
      break;
   case fgkDmemSizeTotalIndividual:
      for(Int_t i=0; i<fgkDmemSizeTotalIndividual; i++) {
	 if(fDmem[addr][i]!=0)
	    AliDebug(5, Form("Warning: Setting new value to DMEM 0x%08x", addr+fgkDmemStartAddress));

	 fDmem[addr][i]=value;
      }
      break;
   case fgkDmemSizeSmRocIndividual:
      for(Int_t i=0; i<fgkDmemSizeSmRocIndividual; i++) {
	 if(fDmem[addr][i]!=0)
	    AliDebug(5, Form("Warning: Setting new value to DMEM 0x%08x", addr+fgkDmemStartAddress));

	 fDmem[addr][i]=value;
      }
      break;
   default:
      AliError(Form("Invalid selection type"));
      break;
   }

   return kTRUE;
}


Bool_t AliTRDtrapConfig::SetDmem(Int_t addr, UInt_t value, Int_t det, Int_t rob, Int_t mcm)
{
  // Set the content of the given DMEM address 

   addr = addr - fgkDmemStartAddress;
   Int_t roc = det%30;
   Int_t loc;
   
   if(addr < 0 || addr >=  fgkDmemWords) {
      AliError(Form("No DMEM address: 0x%08x", addr+fgkDmemStartAddress));
      return kFALSE;
   }

   Int_t detFactor=8*16;
   Int_t robFactor=16;

   switch(fDmemDepth[addr]) {
   case fgkDmemSizeEmpty:
      AliError(Form("DMEM address 0x%08x not active", addr+fgkDmemStartAddress));
      return kFALSE;
      break;
   case fgkDmemSizeUniform:
      if(fDmem[addr][0]!=0)
	 AliDebug(5, Form("Warning: Setting new value to DMEM 0x%08x", addr+fgkDmemStartAddress));

      fDmem[addr][0] = value;
      break;
   case fgkDmemSizeSmIndividual:
      loc = detFactor*roc + robFactor*rob + mcm;
      if(loc < fgkDmemSizeSmIndividual) {
	 if(fDmem[addr][loc]!=0)
	    AliDebug(5, Form("Warning: Setting new value to DMEM 0x%08x", addr+fgkDmemStartAddress));

	 fDmem[addr][loc] = value;
      }
      else {
	 AliError(Form("DMEM sub-address %i out of scope", loc));
	 return kFALSE;
      }
      break;
   case fgkDmemSizeTotalIndividual:
      loc = detFactor*det + robFactor*rob + mcm;
      if(loc < fgkDmemSizeTotalIndividual) {
	 if(fDmem[addr][loc]!=0)
	    AliDebug(5, Form("Warning: Setting new value to DMEM 0x%08x", addr+fgkDmemStartAddress));
	 
	 fDmem[addr][loc]=value;
      }
      else {
	 AliError(Form("DMEM sub-address %i out of scope", loc));
	 return kFALSE;
      }
      break;
   case fgkDmemSizeSmRocIndividual:
      if(det < fgkDmemSizeSmRocIndividual) {
	 if(fDmem[addr][det]!=0)
	    AliDebug(5, Form("Warning: Setting new value to DMEM 0x%08x", addr+fgkDmemStartAddress));

	 fDmem[addr][det]=value;
      }
      else {
	 AliError(Form("DMEM sub-address %i out of scope", det));
	 return kFALSE;
      }
      
      break;
   default:
      AliError(Form("Invalid selection type"));
      return kFALSE;
      break;
   }
   
   return kTRUE;
}


UInt_t AliTRDtrapConfig::GetDmemUnsigned(Int_t addr) const
{
   addr = addr - fgkDmemStartAddress;
   if(addr >=  fgkDmemWords) {
      AliError(Form("No DMEM address: 0x%08x", addr+fgkDmemStartAddress));
      return 0;
   }

   if(fDmemDepth[addr] == fgkDmemSizeUniform)
      return fDmem[addr][0];
   else {
      AliError(Form("No global DMEM value at 0x%08x", addr+fgkDmemStartAddress));
      return 0;
   }
   return 0;
}

UInt_t AliTRDtrapConfig::GetDmemUnsigned(Int_t addr, Int_t det, Int_t rob, Int_t mcm) const
{
   addr = addr - fgkDmemStartAddress;
   Int_t roc = det%30;
   Int_t loc;
   
   if(addr < 0 || addr >=  fgkDmemWords) {
      AliError(Form("No DMEM address: 0x%08x", addr+fgkDmemStartAddress));
      return 0;
   }

   Int_t detFactor=8*16;
   Int_t robFactor=16;

   switch(fDmemDepth[addr]) {
   case fgkDmemSizeEmpty:
      AliError(Form("DMEM address 0x%08x not active", addr+fgkDmemStartAddress));
      return 0;
      break;
   case fgkDmemSizeUniform:
      return fDmem[addr][0];
      break;
   case fgkDmemSizeSmIndividual:
      loc = detFactor*roc + robFactor*rob + mcm;
      if(loc < fgkDmemSizeSmIndividual) {
	 return fDmem[addr][loc];
      }
      else {
	 AliError(Form("DMEM sub-address %i out of scope", loc));
	 return 0;
      }
      break;
   case fgkDmemSizeTotalIndividual:
      loc = detFactor*det + robFactor*rob + mcm;
      if(loc < fgkDmemSizeTotalIndividual) {
	 return fDmem[addr][loc];
      }
      else {
	 AliError(Form("DMEM sub-address %i out of scope", loc));
	 return 0;
      }
      break;
   case fgkDmemSizeSmRocIndividual:
      if(det < fgkDmemSizeSmRocIndividual) {
	 return fDmem[addr][det];
      }
      else {
	 AliError(Form("DMEM sub-address %i out of scope", det));
	 return 0;
      }
      break;
   default:
      AliError(Form("Invalid selection type"));
      return 0;
      break;
   }
   
   return 0;
}


Bool_t AliTRDtrapConfig::LoadConfig()
{
  // load a set of TRAP register values (configuration)
  // here a default set is implemented for testing
  // for a detailed description of the registers see the TRAP manual

  // no. of timebins
  SetTrapReg(kC13CPUA, 24); 

  // pedestal filter
  SetTrapReg(kFPNP, 4*10);
  SetTrapReg(kFPTC, 0);
  SetTrapReg(kFPBY, 0); // bypassed!
  
  // gain filter
  for (Int_t adc = 0; adc < 20; adc++) {
    SetTrapReg(TrapReg_t(kFGA0+adc), 40);
    SetTrapReg(TrapReg_t(kFGF0+adc), 15);
  }
  SetTrapReg(kFGTA, 20);
  SetTrapReg(kFGTB, 2060);
  SetTrapReg(kFGBY, 0);  // bypassed!

  // tail cancellation
  SetTrapReg(kFTAL, 267);
  SetTrapReg(kFTLL, 356);
  SetTrapReg(kFTLS, 387);
  SetTrapReg(kFTBY, 0);

  // tracklet calculation
  SetTrapReg(kTPQS0, 5);
  SetTrapReg(kTPQE0, 10);
  SetTrapReg(kTPQS1, 11);
  SetTrapReg(kTPQE1, 20);
  SetTrapReg(kTPFS, 5);
  SetTrapReg(kTPFE, 20);
  SetTrapReg(kTPVBY, 0);
  SetTrapReg(kTPVT, 10);
  SetTrapReg(kTPHT, 150);
  SetTrapReg(kTPFP, 40);
  SetTrapReg(kTPCL, 1);
  SetTrapReg(kTPCT, 10);

  // ndrift (+ 5 binary digits)
  SetDmem(0xc025, 20 << 5);
  // deflection + tilt correction
  SetDmem(0xc022, 0); 
  // deflection range table
  for (Int_t iTrklCh = 0; iTrklCh < 18; iTrklCh++) {
    SetDmem(0xc030 + 2 * iTrklCh, -64); // min. deflection
    SetDmem(0xc031 + 2 * iTrklCh,  63); // max. deflection
  }
  
  // hit position LUT
  const UShort_t lutPos[128] = {
    0,  1,  1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8,  9,  9, 10, 10, 11, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15,
    16, 16, 16, 17, 17, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24, 24, 25, 25, 25, 26, 26, 26, 26,
    27, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 26,
    26, 26, 26, 25, 25, 25, 24, 24, 23, 23, 22, 22, 21, 21, 20, 20, 19, 18, 18, 17, 17, 16, 15, 14, 13, 12, 11, 10,  9,  8,  7,  7};
  for (Int_t iCOG = 0; iCOG < 128; iCOG++)
    SetTrapReg((TrapReg_t) (kTPL00 + iCOG), lutPos[iCOG]);

  // event buffer
  SetTrapReg(kEBSF, 1);  // 0: store filtered; 1: store unfiltered
  // zs applied to data stored in event buffer (sel. by EBSF)
  SetTrapReg(kEBIS, 15 << 2); // single indicator threshold (plus two digits)
  SetTrapReg(kEBIT, 30 << 2); // sum indicator threshold (plus two digits)
  SetTrapReg(kEBIL, 0xf0);   // lookup table
  SetTrapReg(kEBIN, 0);      // neighbour sensitivity

  // raw data
  SetTrapReg(kNES, (0x0000 << 16) | 0x1000);

  return kTRUE;
}


Bool_t  AliTRDtrapConfig::LoadConfig(Int_t det, TString filename)
{
   // load a TRAP configuration from a file
   // The file format is the format created by the standalone 
   // command coder: scc / show_cfdat 
   // which are two tools to inspect/export configurations from wingDB

  ResetRegs(); // does not really make sense here???

  std::ifstream infile;
  infile.open(filename.Data(), std::ifstream::in);
  if (!infile.is_open()) {
    AliError("Can not open MCM configuration file");
    return kFALSE;
  }

  Int_t cmd, extali, addr, data;
  Int_t no;
  char tmp;
  
  while(infile.good()) {
    cmd=-1;
    extali=-1;
    addr=-1;
    data=-1;
    infile >> std::skipws >> no >> tmp >> cmd >> addr >> data >> extali;
    //      std::cout << "no: " << no << ", cmd " << cmd << ", extali " << extali << ", addr " << addr << ", data " << data <<  endl;
    
    if(cmd!=-1 && extali!=-1 && addr != -1 && data!= -1) {
      AddValues(det, cmd, extali, addr, data);
    }
    else if(!infile.eof() && !infile.good()) {
      infile.clear();
      infile.ignore(256, '\n');
    }
    
    if(!infile.eof())
      infile.clear();
  }
  
  infile.close();
  
  return kTRUE;
}


Bool_t AliTRDtrapConfig::ReadPackedConfig(Int_t hc, UInt_t *data, Int_t size) 
{
  // Read the packed configuration from the passed memory block
  //
  // To be used to retrieve the TRAP configuration from the 
  // configuration as sent in the raw data. 

  AliDebug(1, "Reading packed configuration");

  Int_t det = hc/2;

  Int_t idx = 0;
  Int_t err = 0;
  Int_t step, bwidth, nwords, exitFlag, bitcnt;
  
  UShort_t caddr;
  UInt_t dat, msk, header, dataHi;
  
  while (idx < size && *data != 0x00000000) {
    
    Int_t rob = (*data >> 28) & 0x7;
    Int_t mcm = (*data >> 24) & 0xf;

    AliDebug(1, Form("Config of det. %3i MCM %i:%02i (0x%08x)", det, rob, mcm, *data));
    data++;
    
    while (idx < size && *data != 0x00000000) {
      
      header = *data;
      data++;
      idx++;
      
      AliDebug(5, Form("read: 0x%08x", header));
      
      if (header & 0x01) // single data
	{
	  dat   = (header >>  2) & 0xFFFF;       // 16 bit data 
	  caddr = (header >> 18) & 0x3FFF;    // 14 bit address 
	  
	  if (caddr != 0x1FFF)  // temp!!! because the end marker was wrong
	    {
	      if (header & 0x02) // check if > 16 bits
		{
		  dataHi = *data;
		  AliDebug(5, Form("read: 0x%08x", dataHi));
		  data++;
		  idx++;
		  err += ((dataHi ^ (dat | 1)) & 0xFFFF) != 0;
		  dat = (dataHi & 0xFFFF0000) | dat;
		}
	      AliDebug(5, Form("addr=0x%04x (%s) data=0x%08x\n", caddr, GetRegName(GetRegByAddress(caddr)), dat));
	      if ( ! Poke(caddr, dat, det, rob, mcm) )
		AliDebug(5, Form("(single-write): non-existing address 0x%04x containing 0x%08x\n", caddr, header));
	      if (idx > size)
		{
		  AliDebug(5, Form("(single-write): no more data, missing end marker\n"));
		  return -err;
		}
	    }
	  else
	    {
	      AliDebug(5, Form("(single-write): address 0x%04x => old endmarker?\n", caddr));
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
	  
	  if (exitFlag) 
	    break;
	  
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
			header = *data;
			AliDebug(5, Form("read 0x%08x", header));
			data++;
			idx++;
			err += (header & 1);
			header = header >> 1;
			bitcnt = 31 - bwidth;
		      }
		    AliDebug(5, Form("addr=0x%04x (%s) data=0x%08x\n", caddr, GetRegName(GetRegByAddress(caddr)), header & msk));
		    if ( ! Poke(caddr, header & msk, det, rob, mcm) )
		      AliDebug(5, Form("(single-write): non-existing address 0x%04x containing 0x%08x\n", caddr, header));
		    
		    caddr += step;
		    header = header >> bwidth;
		    if (idx >= size)
		      {
			AliDebug(5, Form("(block-write): no end marker! %d words read\n", idx));
			return -err;
		      }
		  }
		break;
	      } // end case 5-15                                         
	    case 31:
	      {
		while (nwords > 0)
		  {
		    header = *data;
		    AliDebug(5, Form("read 0x%08x", header));
		    data++;
		    idx++;
		    nwords--;
		    err += (header & 1);
		    
		    AliDebug(5, Form("addr=0x%04x (%s) data=0x%08x", caddr, GetRegName(GetRegByAddress(caddr)), header >> 1));
		    if ( ! Poke(caddr, header >> 1, det, rob, mcm) )
		      AliDebug(5, Form("(single-write): non-existing address 0x%04x containing 0x%08x\n", caddr, header));
		    
		    caddr += step;
		    if (idx >= size)
		      {
			AliDebug(5, Form("no end marker! %d words read", idx));
			return -err;
		      }
		  }
		break;
	      }
	    default: return err;
	    } // end switch 
	} // end block case
    }
  } // end while
  AliDebug(5, Form("no end marker! %d words read", idx));
  return -err; // only if the max length of the block reached!                       
}


Bool_t AliTRDtrapConfig::PrintTrapReg(TrapReg_t reg, Int_t det, Int_t rob, Int_t mcm) const
{
  // print the value stored in the given register
  // if it is individual a valid MCM has to be specified

  if (fRegisterValue[reg].state == RegValue_t::kGlobal) {
    printf("%s (%i bits) at 0x%08x is 0x%08x and resets to: 0x%08x (currently global mode)\n", 
           GetRegName((TrapReg_t) reg),
           GetRegNBits((TrapReg_t) reg),
           GetRegAddress((TrapReg_t) reg),
           fRegisterValue[reg].globalValue,
           GetRegResetValue((TrapReg_t) reg));
  }
  else if (fRegisterValue[reg].state == RegValue_t::kIndividual) {
    if((det >= 0 && det < AliTRDgeometry::Ndet()) && 
       (rob >= 0 && rob < AliTRDfeeParam::GetNrobC1()) && 
       (mcm >= 0 && mcm < fgkMaxMcm)) {
      printf("%s (%i bits) at 0x%08x is 0x%08x and resets to: 0x%08x (currently individual mode)\n", 
             GetRegName((TrapReg_t) reg),
             GetRegNBits((TrapReg_t) reg),
             GetRegAddress((TrapReg_t) reg),
             fRegisterValue[reg].individualValue[det*AliTRDfeeParam::GetNrobC1()*fgkMaxMcm + rob*fgkMaxMcm + mcm],
             GetRegResetValue((TrapReg_t) reg));
    }
    else {
      AliError("Register value is MCM-specific: Invalid detector, ROB or MCM requested");
      return kFALSE;
    }
  }
  else {  // should never be reached
    AliError("MCM register status neither kGlobal nor kIndividual");
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


Bool_t AliTRDtrapConfig::AddValues(UInt_t det, UInt_t cmd, UInt_t extali, Int_t addr, UInt_t data)
{
   // transfer the informations provided by LoadConfig to the internal class variables

  if(cmd != fgkScsnCmdWrite) {
    AliError(Form("Invalid command received: %i", cmd));
    return kFALSE;
  }

  TrapReg_t mcmReg = GetRegByAddress(addr);
  Int_t rocType = AliTRDgeometry::GetStack(det) == 2 ? 0 : 1;

  static const int mcmListSize=40;  // 40 is more or less arbitrary
  Int_t mcmList[mcmListSize];

  // configuration registers
  if(mcmReg >= 0 && mcmReg < kLastReg) {
    
    for(Int_t linkPair=0; linkPair<fgkMaxLinkPairs; linkPair++) {
      if(AliTRDfeeParam::ExtAliToAli(extali, linkPair, rocType, mcmList, mcmListSize)!=0) {
	Int_t i=0;
        while(mcmList[i] != -1 && i<mcmListSize) {
          if(mcmList[i]==127)
            SetTrapReg( (TrapReg_t) mcmReg, data, det);
          else
            SetTrapReg( (TrapReg_t) mcmReg, data, det, (mcmList[i]>>7), (mcmList[i]&0x7F));
          i++;
        }
      }
    }
    return kTRUE;
  }
  // DMEM
  else if ( (addr >= fgkDmemStartAddress) && 
	    (addr < (fgkDmemStartAddress + fgkDmemWords))) {
    for(Int_t linkPair=0; linkPair<fgkMaxLinkPairs; linkPair++) {
      if(AliTRDfeeParam::ExtAliToAli(extali, linkPair, rocType, mcmList, mcmListSize)!=0) {
        Int_t i=0;
        while(mcmList[i] != -1 && i < mcmListSize) {
          if(mcmList[i] == 127)
	     SetDmem(addr, data, det, 0, 127);
          else
	     SetDmem(addr, data, det, mcmList[i] >> 7, mcmList[i] & 0x7f);
          i++;
        }
      }
    }
    return kTRUE;
  }
  else 
    return kFALSE;
}


AliTRDtrapConfig::TrapReg_t AliTRDtrapConfig::GetRegByAddress(Int_t address) const
{
  // get register by its address
  // used for reading of configuration data as sent to real FEE

  TrapReg_t mcmReg = kLastReg;
  Int_t reg  = 0;
  do {
    if(fRegs[reg].fAddr == address)
      mcmReg = (TrapReg_t) reg;
    reg++;
  }  while (mcmReg == kLastReg && reg < kLastReg);

  return mcmReg;
}


void AliTRDtrapConfig::PrintMemDatx(ostream &os, Int_t addr) const
{
   PrintMemDatx(os, addr, 0, 0, 127);
}

void AliTRDtrapConfig::PrintMemDatx(ostream &os, Int_t addr, Int_t det, Int_t rob, Int_t mcm) const
{
   if(addr < fgkDmemStartAddress || addr >= fgkDmemStartAddress+fgkDmemWords) {
      AliError(Form("Invalid DMEM address 0x%08x!", addr));
      return;
   }
   PrintDatx(os, addr, GetDmemUnsigned(addr, det, rob, mcm), rob, mcm);
}


void AliTRDtrapConfig::PrintMemDatx(ostream &os, TrapReg_t reg) const
{
   PrintMemDatx(os, reg, 0, 0, 127);
}


void AliTRDtrapConfig::PrintMemDatx(ostream &os, TrapReg_t reg, Int_t det, Int_t rob, Int_t mcm) const
{
   if(reg>= kLastReg) {
      AliError(Form("Invalid register %i!", reg));
      return;
   }
   PrintDatx(os, GetRegAddress(reg), GetTrapReg(reg, det, rob, mcm), rob, mcm);
}


void AliTRDtrapConfig::PrintDatx(ostream &os, UInt_t addr, UInt_t data, Int_t rob, Int_t mcm) const
{
   os << std::setw(5) << 10 
      << std::setw(8) << addr
      << std::setw(12) << data;
   if(mcm==127)
      os << std::setw(8) << 127;
   else
      os << std::setw(8) << AliTRDfeeParam::AliToExtAli(rob, mcm);
   
   os << std::endl;
}
