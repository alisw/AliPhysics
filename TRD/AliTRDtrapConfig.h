#ifndef ALITRDTRAPCONFIG_H
#define ALITRDTRAPCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// 
// Class holding the configuration of the tracklet processor 
// in the TRD FEE
//


#include <TObject.h>
#include <TString.h>
#include <fstream>

class AliTRDtrapConfig : public TObject
{
 public:
  ~AliTRDtrapConfig(); 

  static AliTRDtrapConfig* Instance();

  enum TrapReg_t { kSML0, 
		  kSML1, 
		  kSML2,
		  kSMMODE,
		  kNITM0,
		  kNITM1,
		  kNITM2,
		  kNIP4D,
		  kCPU0CLK,
		  kCPU1CLK,
		  kCPU2CLK,
		  kCPU3CLK,
		  kNICLK,
		  kFILCLK,
		  kPRECLK,
		  kADCEN,
		  kNIODE,
		  kNIOCE,
		  kNIIDE,
		  kNIICE,
		  kARBTIM,      
		  kIA0IRQ0,     
		  kIA0IRQ1,     
		  kIA0IRQ2,     
		  kIA0IRQ3,     
		  kIA0IRQ4,     
		  kIA0IRQ5,     
		  kIA0IRQ6,     
		  kIA0IRQ7,     
		  kIA0IRQ8,     
		  kIA0IRQ9,     
		  kIA0IRQA,     
		  kIA0IRQB,     
		  kIA0IRQC,     
		  kIRQSW0,      
		  kIRQHW0,      
		  kIRQHL0,      
		  kIA1IRQ0,     
		  kIA1IRQ1,     
		  kIA1IRQ2,     
		  kIA1IRQ3,     
		  kIA1IRQ4,     
		  kIA1IRQ5,     
		  kIA1IRQ6,     
		  kIA1IRQ7,     
		  kIA1IRQ8,     
		  kIA1IRQ9,     
		  kIA1IRQA,     
		  kIA1IRQB,     
		  kIA1IRQC,     
		  kIRQSW1,      
		  kIRQHW1,      
		  kIRQHL1,      
		  kIA2IRQ0,     
		  kIA2IRQ1,     
		  kIA2IRQ2,     
		  kIA2IRQ3,     
		  kIA2IRQ4,     
		  kIA2IRQ5,     
		  kIA2IRQ6,     
		  kIA2IRQ7,     
		  kIA2IRQ8,     
		  kIA2IRQ9,     
		  kIA2IRQA,     
		  kIA2IRQB,     
		  kIA2IRQC,     
		  kIRQSW2,      
		  kIRQHW2,      
		  kIRQHL2,      
		  kIA3IRQ0,     
		  kIA3IRQ1,     
		  kIA3IRQ2,     
		  kIA3IRQ3,     
		  kIA3IRQ4,     
		  kIA3IRQ5,     
		  kIA3IRQ6,     
		  kIA3IRQ7,     
		  kIA3IRQ8,     
		  kIA3IRQ9,     
		  kIA3IRQA,     
		  kIA3IRQB,     
		  kIA3IRQC,     
		  kIRQSW3,      
		  kIRQHW3,      
		  kIRQHL3,      
		  kCTGDINI,     
		  kCTGCTRL,     
		  kC08CPU0,     
		  kC09CPU0,     
		  kC10CPU0,     
		  kC11CPU0,     
		  kC12CPUA,     
		  kC13CPUA,     
		  kC14CPUA,     
		  kC15CPUA,     
		  kC08CPU1,     
		  kC09CPU1,     
		  kC10CPU1,     
		  kC11CPU1,     
		  kC08CPU2,     
		  kC09CPU2,     
		  kC10CPU2,     
		  kC11CPU2,     
		  kC08CPU3,     
		  kC09CPU3,     
		  kC10CPU3,     
		  kC11CPU3,     
		  kNMOD,        
		  kNDLY,        
		  kNED,         
		  kNTRO,        
		  kNRRO,        
		  kNES,         
		  kNTP,         
		  kNBND,        
		  kNP0,         
		  kNP1,         
		  kNP2,         
		  kNP3,         
		  kNCUT,        
		  kTPPT0,       
		  kTPFS,        
		  kTPFE,        
		  kTPPGR,       
		  kTPPAE,       
		  kTPQS0,       
		  kTPQE0,       
		  kTPQS1,       
		  kTPQE1,       
		  kEBD,         
		  kEBAQA,       
		  kEBSIA,       
		  kEBSF,        
		  kEBSIM,       
		  kEBPP,        
		  kEBPC,        
		  kEBIS,        
		  kEBIT,        
		  kEBIL,        
		  kEBIN,        
		  kFLBY,        
		  kFPBY,        
		  kFGBY,        
		  kFTBY,        
		  kFCBY,        
		  kFPTC,        
		  kFPNP,        
		  kFPCL,        
		  kFGTA,        
		  kFGTB,        
		  kFGCL,        
		  kFTAL,        
		  kFTLL,        
		  kFTLS,        
		  kFCW1,        
		  kFCW2,        
		  kFCW3,        
		  kFCW4,        
		  kFCW5,        
		  kTPFP,        
		  kTPHT,        
		  kTPVT,        
		  kTPVBY,       
		  kTPCT,        
		  kTPCL,        
		  kTPCBY,       
		  kTPD,         
		  kTPCI0,       
		  kTPCI1,       
		  kTPCI2,       
		  kTPCI3,       
		  kADCMSK,      
		  kADCINB,      
		  kADCDAC,      
		  kADCPAR,      
		  kADCTST,      
		  kSADCAZ,      
		  kFGF0,        
		  kFGF1,        
		  kFGF2,        
		  kFGF3,        
		  kFGF4,        
		  kFGF5,        
		  kFGF6,        
		  kFGF7,        
		  kFGF8,        
		  kFGF9,        
		  kFGF10,       
		  kFGF11,       
		  kFGF12,       
		  kFGF13,       
		  kFGF14,       
		  kFGF15,       
		  kFGF16,       
		  kFGF17,       
		  kFGF18,       
		  kFGF19,       
		  kFGF20,       
		  kFGA0,        
		  kFGA1,        
		  kFGA2,        
		  kFGA3,        
		  kFGA4,        
		  kFGA5,        
		  kFGA6,        
		  kFGA7,        
		  kFGA8,        
		  kFGA9,        
		  kFGA10,       
		  kFGA11,       
		  kFGA12,       
		  kFGA13,       
		  kFGA14,       
		  kFGA15,       
		  kFGA16,       
		  kFGA17,       
		  kFGA18,       
		  kFGA19,       
		  kFGA20,       
		  kFLL00,       
		  kFLL01,       
		  kFLL02,       
		  kFLL03,       
		  kFLL04,       
		  kFLL05,       
		  kFLL06,       
		  kFLL07,       
		  kFLL08,       
		  kFLL09,       
		  kFLL0A,       
		  kFLL0B,       
		  kFLL0C,       
		  kFLL0D,       
		  kFLL0E,       
		  kFLL0F,       
		  kFLL10,       
		  kFLL11,       
		  kFLL12,       
		  kFLL13,       
		  kFLL14,       
		  kFLL15,       
		  kFLL16,       
		  kFLL17,       
		  kFLL18,       
		  kFLL19,       
		  kFLL1A,       
		  kFLL1B,       
		  kFLL1C,       
		  kFLL1D,       
		  kFLL1E,       
		  kFLL1F,       
		  kFLL20,       
		  kFLL21,       
		  kFLL22,       
		  kFLL23,       
		  kFLL24,       
		  kFLL25,       
		  kFLL26,       
		  kFLL27,       
		  kFLL28,       
		  kFLL29,       
		  kFLL2A,       
		  kFLL2B,       
		  kFLL2C,       
		  kFLL2D,       
		  kFLL2E,       
		  kFLL2F,       
		  kFLL30,       
		  kFLL31,       
		  kFLL32,       
		  kFLL33,       
		  kFLL34,       
		  kFLL35,       
		  kFLL36,       
		  kFLL37,       
		  kFLL38,       
		  kFLL39,       
		  kFLL3A,       
		  kFLL3B,       
		  kFLL3C,       
		  kFLL3D,       
		  kFLL3E,       
		  kFLL3F,       
		  kPASADEL,     
		  kPASAPHA,     
		  kPASAPRA,     
		  kPASADAC,     
		  kPASACHM,     
		  kPASASTL,     
		  kPASAPR1,     
		  kPASAPR0,     
		  kSADCTRG,     
		  kSADCRUN,     
		  kSADCPWR,     
		  kL0TSIM,      
		  kSADCEC,      
		  kSADCMC,      
		  kSADCOC,      
		  kSADCGTB,     
		  kSEBDEN,      
		  kSEBDOU,      
		  kTPL00,       
		  kTPL01,       
		  kTPL02,       
		  kTPL03,       
		  kTPL04,       
		  kTPL05,       
		  kTPL06,       
		  kTPL07,       
		  kTPL08,       
		  kTPL09,       
		  kTPL0A,       
		  kTPL0B,       
		  kTPL0C,       
		  kTPL0D,       
		  kTPL0E,       
		  kTPL0F,       
		  kTPL10,       
		  kTPL11,       
		  kTPL12,       
		  kTPL13,       
		  kTPL14,       
		  kTPL15,       
		  kTPL16,       
		  kTPL17,       
		  kTPL18,       
		  kTPL19,       
		  kTPL1A,       
		  kTPL1B,       
		  kTPL1C,       
		  kTPL1D,       
		  kTPL1E,       
		  kTPL1F,       
		  kTPL20,       
		  kTPL21,       
		  kTPL22,       
		  kTPL23,       
		  kTPL24,       
		  kTPL25,       
		  kTPL26,       
		  kTPL27,       
		  kTPL28,       
		  kTPL29,       
		  kTPL2A,       
		  kTPL2B,       
		  kTPL2C,       
		  kTPL2D,       
		  kTPL2E,       
		  kTPL2F,       
		  kTPL30,       
		  kTPL31,       
		  kTPL32,       
		  kTPL33,       
		  kTPL34,       
		  kTPL35,       
		  kTPL36,       
		  kTPL37,       
		  kTPL38,       
		  kTPL39,       
		  kTPL3A,       
		  kTPL3B,       
		  kTPL3C,       
		  kTPL3D,       
		  kTPL3E,       
		  kTPL3F,       
		  kTPL40,       
		  kTPL41,       
		  kTPL42,       
		  kTPL43,       
		  kTPL44,       
		  kTPL45,       
		  kTPL46,       
		  kTPL47,       
		  kTPL48,       
		  kTPL49,       
		  kTPL4A,       
		  kTPL4B,       
		  kTPL4C,       
		  kTPL4D,       
		  kTPL4E,       
		  kTPL4F,       
		  kTPL50,       
		  kTPL51,       
		  kTPL52,       
		  kTPL53,       
		  kTPL54,       
		  kTPL55,       
		  kTPL56,       
		  kTPL57,       
		  kTPL58,       
		  kTPL59,       
		  kTPL5A,       
		  kTPL5B,       
		  kTPL5C,       
		  kTPL5D,       
		  kTPL5E,       
		  kTPL5F,       
		  kTPL60,       
		  kTPL61,       
		  kTPL62,       
		  kTPL63,       
		  kTPL64,       
		  kTPL65,       
		  kTPL66,       
		  kTPL67,       
		  kTPL68,       
		  kTPL69,       
		  kTPL6A,       
		  kTPL6B,       
		  kTPL6C,       
		  kTPL6D,       
		  kTPL6E,       
		  kTPL6F,       
		  kTPL70,       
		  kTPL71,       
		  kTPL72,       
		  kTPL73,       
		  kTPL74,       
		  kTPL75,       
		  kTPL76,       
		  kTPL77,       
		  kTPL78,       
		  kTPL79,       
		  kTPL7A,       
		  kTPL7B,       
		  kTPL7C,       
		  kTPL7D,       
		  kTPL7E,       
		  kTPL7F,       
		  kMEMRW,       
		  kMEMCOR,      
		  kDMDELA,      
		  kDMDELS,      
		  kLastReg };   // enum of all TRAP registers, to be used for access to them

  const char* GetRegName(TrapReg_t reg)       const { return reg >= 0 && reg < kLastReg ? fRegs[reg].fName.Data() : ""; }
  UShort_t    GetRegAddress(TrapReg_t reg)    const { return reg >= 0 && reg < kLastReg ? fRegs[reg].fAddr : 0; }
  UShort_t    GetRegNBits(TrapReg_t reg)      const { return reg >= 0 && reg < kLastReg ? fRegs[reg].fNbits : 0; }
  UInt_t      GetRegResetValue(TrapReg_t reg) const { return reg >= 0 && reg < kLastReg ? fRegs[reg].fResetValue : 0; }

  TrapReg_t          GetRegByAddress(Int_t address) const;

  Int_t  GetTrapReg(TrapReg_t reg, Int_t det = -1, Int_t rob = -1, Int_t mcm = -1) const;
  Bool_t PrintTrapReg(TrapReg_t reg, Int_t det = -1, Int_t rob = -1, Int_t mcm = -1) const;
  Bool_t PrintTrapAddr(Int_t addr, Int_t det = -1, Int_t rob = -1, Int_t mcm = -1) const;

  Bool_t SetTrapReg(TrapReg_t reg, Int_t value);
  Bool_t SetTrapReg(TrapReg_t reg, Int_t value, Int_t det);
  Bool_t SetTrapReg(TrapReg_t reg, Int_t value, Int_t det, Int_t rob, Int_t mcm);

  UInt_t Peek(Int_t addr, Int_t det, Int_t rob, Int_t mcm) const;
  Bool_t Poke(Int_t addr, UInt_t value, Int_t det, Int_t rob, Int_t mcm);

  void InitRegs();
  void ResetRegs();
  void ResetDmem();

  // DMEM 
  Bool_t SetDmem(Int_t addr, UInt_t value);
  //  Bool_t SetDmem(Int_t addr, UInt_t value, Int_t det);
  Bool_t SetDmem(Int_t addr, UInt_t value, Int_t det, Int_t rob, Int_t mcm);
  Bool_t SetDmem(Int_t addr, Int_t value) { return SetDmem(addr, (UInt_t) value); }
  //  Bool_t SetDmem(Int_t addr, Int_t value, Int_t det) { return SetDmem(addr, (UInt_t) value, det); }
  Bool_t SetDmem(Int_t addr, Int_t value, Int_t det, Int_t rob, Int_t mcm) { return SetDmem(addr, (UInt_t) value, det, rob, mcm); }

  UInt_t GetDmemUnsigned(Int_t addr) const;
  UInt_t GetDmemUnsigned(Int_t addr, Int_t det, Int_t rob, Int_t mcm) const;

  void PrintMemDatx(ostream &os, Int_t addr) const;
  void PrintMemDatx(ostream &os, Int_t addr, Int_t det, Int_t rob, Int_t mcm) const;
  void PrintMemDatx(ostream &os, TrapReg_t reg) const;
  void PrintMemDatx(ostream &os, TrapReg_t reg, Int_t det, Int_t rob, Int_t mcm) const;
  void PrintDatx(ostream &os, UInt_t addr, UInt_t data, Int_t rob, Int_t mcm) const;
 
  // PID
  void SetPIDscale(Double_t sq0, Double_t sq1) {fScaleQ0=sq0; fScaleQ1=sq1;}  // The two scaling factors are not available in the TRAP but needed by the DCS board
  void GetPIDscale(Double_t scale[2]) {scale[0] = fScaleQ0; scale[1]=fScaleQ1;}  // The two scaling factors are not available in the TRAP but needed by the DCS board

  // configuration handling
  Bool_t LoadConfig();

  Bool_t ReadPackedConfig(Int_t det, UInt_t *data, Int_t size);

  Bool_t AddValues(UInt_t det, UInt_t cmd, UInt_t extali, Int_t addr, UInt_t data);

  // DMEM addresses
  static const Int_t fgkDmemAddrLUTcor0       = 0xC02A;
  static const Int_t fgkDmemAddrLUTcor1       = 0xC028;
  static const Int_t fgkDmemAddrLUTnbins      = 0xC029;
					      
  static const Int_t fgkDmemAddrLUTStart      = 0xC100; // LUT start address
  static const Int_t fgkDmemAddrLUTEnd        = 0xC3FF; // maximum possible end address for the LUT table
  static const Int_t fgkDmemAddrLUTLength     = 0xC02B; // address where real size of the LUT table is stored
					      
  static const Int_t fgkDmemAddrTrackletStart = 0xC0E0; // Storage area for tracklets, start address
  static const Int_t fgkDmemAddrTrackletEnd   = 0xC0E3; // Storage area for tracklets, end address

  static const Int_t fgkDmemAddrDeflCorr      = 0xc022; // DMEM address of deflection correction
  static const Int_t fgkDmemAddrNdrift        = 0xc025; // DMEM address of Ndrift
  static const Int_t fgkDmemAddrDeflCutStart  = 0xc030; // DMEM start address of deflection cut 
  static const Int_t fgkDmemAddrDeflCutEnd    = 0xc055; // DMEM end address of deflection cut

  // DMEM memory in simulation;
  static const Int_t fgkDmemStartAddress;           // start address in TRAP GIO
  static const Int_t fgkDmemWords = 0xc400;         // number of words in DMEM

 protected:
  static AliTRDtrapConfig *fgInstance;  // pointer to instance (singleton)

  struct SimpleReg_t {
    TString   fName;       // Name of the register
    UShort_t  fAddr;       // Address in GIO of TRAP
    UShort_t  fNbits;      // Number of bits, from 1 to 32
    UInt_t    fResetValue; // reset value
    SimpleReg_t(const char *nnn = 0, UShort_t a = 0, UShort_t n = 0, UInt_t r = 0) : 
      fName(nnn), fAddr(a), fNbits(n), fResetValue(r) {}
  };

  struct RegValue_t {
    enum {
      kInvalid = 0,
      kGlobal,
      kIndividual
    } state; // mode of storage (global or per MCM)
    union {
      Int_t globalValue;
      Int_t *individualValue;
    };
  };

  // configuration registers
  SimpleReg_t fRegs[kLastReg];          // array of TRAP registers
  RegValue_t fRegisterValue[kLastReg];  // array of TRAP register values in use

  static const UInt_t fgkScsnCmdWrite=10;  // Command number for the write command 
  static const Int_t fgkMaxLinkPairs=4;    // number of linkpairs used during configuration
  static const Int_t fgkMaxMcm;            // max. no. of MCMs to be treated
  static const Int_t fgkMcmlistSize=256;     // list of MCMs to which a value has to be written

  // DMEM
  UInt_t* fDmem[fgkDmemWords]; // DMEM storage
  Int_t fDmemDepth[fgkDmemWords]; // memory space indicator for fDmem

  static const Int_t fgkDmemSizeEmpty=0;
  static const   Int_t fgkDmemSizeUniform = 1;
  static const   Int_t fgkDmemSizeSmIndividual = 30*8*16;   // storage for each MCM within one supermodule
  static const   Int_t fgkDmemSizeTotalIndividual = 540*6*8*16;  // one individual value for each and every MCM in the TRD
  static const   Int_t fgkDmemSizeSmRocIndividual = 540; // one individual value for each chamber in TRD

  // Online PID
  Double_t fScaleQ0;  // scaling factor for the x-axis of the PID table
  Double_t fScaleQ1;  // scaling factor for the y-axis of the PID table

  AliTRDtrapConfig(); // private constructor due to singleton implementation

/* not yet used
  struct BlockDescr_t {
    UShort_t addr;
    UShort_t nregs;
    UShort_t nbits;
    UShort_t step;
  };

  struct CmdReg_t {
    char *name;
    UShort_t addr;
  };

  enum DbankProp_t { kDBankEmpty = 0,
		    kDBankHeader,
		    kDBankData,
		    kDBankNoB,
		    kDBankCRC,
		    kDBankEHeader,
		    kDBankSCSNData };
*/


 private:
  AliTRDtrapConfig& operator=(const AliTRDtrapConfig &rhs); // not implemented
  AliTRDtrapConfig(const AliTRDtrapConfig& cfg); // not implemented

  ClassDef(AliTRDtrapConfig, 2);
};

#endif

