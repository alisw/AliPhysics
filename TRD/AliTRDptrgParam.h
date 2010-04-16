#ifndef ALITRDPTRGPARAM_H
#define ALITRDPTRGPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// --------------------------------------------------------
// 
// Singleton class to hold the parameters steering the PTRG 
// 
// --------------------------------------------------------
#include "TObject.h"

  enum AliTRDptrgFEBType_t{ kUndefined = (Int_t)0, 
                            kTZERO = (Int_t)1, 
                            kVZERO = (Int_t)2 };
  enum AliTRDptrgOperatingMode_t{ kHits = (Int_t)0, kDigits = (Int_t)1 };
  enum AliTRDptrgFEBPosition_t{ kB = (Int_t)0, 
                                kA = (Int_t)1, 
                                kC = (Int_t)2,
                                kUnknown = (Int_t)3 };


class AliTRDptrgParam : public TObject {
 public:
  struct AliTRDptrgPTmasks {
    Bool_t fCBA[2]; // contribute CB-A look up results to pretrigger decision?
    Bool_t fCBC[2]; // contribute CB-C look up results to pretrigger decision?
    Bool_t fLUTs[2]; // CB-B look up results contribution to pretrigger decision
    Bool_t fTLMU[8]; // TLMU output signal contribution to pretrigger decisions

    AliTRDptrgPTmasks() {
      fCBA[0] = kFALSE;
      fCBA[1] = kFALSE;
      fCBC[0] = kFALSE;
      fCBC[1] = kFALSE;
      fLUTs[0] = kFALSE;
      fLUTs[1] = kFALSE;
      for (Int_t i = 0; i < 8; i++) {
        fTLMU[i] = kFALSE;
      } 
    }
  };
  virtual ~AliTRDptrgParam();

  static AliTRDptrgParam *Instance(); // Singleton
  static void Terminate();  // delete Singleton
  
  void LoadStandardConfiguration(); // initialize with standard values
  Bool_t LoadConfigurationFromFile(TString filename); // load file 
  
  Int_t GenerateLUTs(); // generates all LUTs

  // --- GETTER FUNCTIONS -----------------------------------------------------
  // -- TLMU --
  const UInt_t* GetTLMUInputMask() const { return this->fTLMUInputMask; };
  // get TLMU input mask
   
  UInt_t** GetTLMUcmatrices() const { return this->fTLMUcmatrices; };
  // get TLMU coincidence matrices
  
  UInt_t** GetTLMUmultiplicity() const { return this->fTLMUmultiplicity; };
  // get TLMU multiplicity slices

  Int_t** GetTLMUoutput() const { return this->fTLMUoutput; };
  // get TLMU output mux configuration

  // -- T0 --
  UInt_t* GetFEBT0Thresholds(AliTRDptrgFEBPosition_t FEBposition) const;
  // get T0 FEB Thresholds

  Int_t* GetFEBT0LUT(AliTRDptrgFEBPosition_t FEBposition, Int_t iLUT); 
  // get T0 FEB LUTs
 
  // -- V0 --
  UInt_t* GetFEBV0Thresholds(AliTRDptrgFEBPosition_t FEBposition, Int_t iCard) const;
  // get V0 FEB Thresholds

  Int_t* GetFEBV0LUT(AliTRDptrgFEBPosition_t FEBposition, Int_t iCard, 
                     Int_t iLUT);
  // get V0 FEB LUTs

  Int_t* GetCBLUT(UInt_t CB, Int_t LUTid);
  // returns the corresponding LUT (control boxes only)

  const AliTRDptrgPTmasks* GetPTmasks() const { return &fPTmasks; };
  // returns the list containing the information which CB-B inputs are masked
  // out or forwarded as pre trigger output to the CTP


  Int_t CheckVariables() const; // returns -1 if a variable is already deleted

 protected:
  UInt_t GetMultiplicity(UInt_t BitVector) const; 
  // returns the multiplicity ('1's) 
  
  UInt_t GetMultiplicity(Int_t BitVector) const;  
  // returns the multiplicity ('1's)

  // helper functions for configuration file reading
  // -----------------------------------------------
  Bool_t ParseTLMU(TString identifier, TString value);
  // parses the TLMU configuration parameters

  Bool_t ParseCBB(TString identifier, TString value);
  // parses the CBB configuration parameters
  
  Bool_t ParseCBAC(TString identifier, TString value);
  // parses the CB-A and CB-C configuration parameters

  Bool_t ParseFEB(TString identifier, TString value);
  // parses the FEB configuration parameters

  Bool_t ParseMultiplicityCondition(TString condition, UInt_t* threshold,
                                    UInt_t* mask);
  // parses a multiplicity condition "M(#mask#)>#threshold#"

  UInt_t BinaryTStringToInt(TString number) const;
  // converts TString containing a binary number to a unsigned integer

  void SplitUpValues(TString value, TObjArray& arr);
  // splits a value string which contains multiple values seperated by ' ' 
  // and '\t'

  TString CleanTString(TString string);
  // removes ' ' and '\t' in a TString

  void PrepareLine(TString line, TString& identifier, TString& value);
  // divides identifier and value (seperator is the first ' ' or '\t'  

  
  // (helper) functions for conversion of logical equations into LUTs 
  // ----------------------------------------------------------------
  Int_t LookUp(TString* const identifier) const; // translates an identifier used in a
  // logical equation into an address bit of the corresponding LUT
  
  void MergeResults(TArrayI*& partResult1, TArrayI*& partResult2, 
                    TArrayI*& results, TArrayI*& signalsInvolved1, 
                    TArrayI*& signalsInvolved2, TArrayI*& signalsInvolved,
                    Bool_t useOR);
  // merges the results of to logical equation parts
  
  void ConvertLogicalEqToBitVectors(TString eq, TArrayI*& results,
                                    TArrayI*& signalsInvolved);
  // converts logical equations to bit vectors
  // neglected input signals are for now assumed to be 0!
  
  void CheckSignalsInvolved(TArrayI*& results, TArrayI*& signalsInvolved,
			    Int_t inputWidth);
  // adds all signal combinations needed to behave correctly in every state of
  // neglected signals
  
  Int_t* GenerateLUTbasedOnEq(TString eq, Int_t inputWidth, Int_t initValue);
  // generates a lut based on a logical functions (uses the functions above)

  static AliTRDptrgParam *fgInstance; // instance pointer

  // TLMU configuration --------------------------------------------------------
  UInt_t fTLMUInputMask[18]; // masks TOF-to-TRD bits
  UInt_t fTLMUInputStretch; // designates how long TLMU input is stretched   
  UInt_t** fTLMUcmatrices; // [matrix][section] unsigned int values
  // Bits 0..17 identify supermodules, bits equal 1 are checked for coincidence

  UInt_t** fTLMUmultiplicity; // [slice][0 = lower bound, 1 = upper bound]
  // use a lower bound above 576 to disable
  
  Int_t** fTLMUoutput; // [output][0 = cmatrix, 1 = multslice] 
  // output bit assignment, -1 disables

  // T0 ------------------------------------------------------------------------
  // [position][channel] 12 channels at A and C side
  UInt_t** fFEBT0Thresholds; // threshold for analog value discrimination
  
  // [position][LUT][0 = threshold, 1 = bitmask] 2 LUTs at A and C side  
  UInt_t*** fFEBT0Multiplicities; // multiplicity threshold for T0
  Int_t*** fFEBT0LUTs; // look up tables [position][LUT][entry]
    
  // V0 ------------------------------------------------------------------------
  // [position][feb][channel] 4x8 channels per side (A and C)
  UInt_t*** fFEBV0Thresholds; // threshold for analog value discrimation

  // [position][feb][LUT][0 = threshold, 1 = bitmask] 2 LUTs per FEB 
  // (4 per Side) at each side ( A and C)
  UInt_t**** fFEBV0Multiplicities; // multiplicity threshold for V0   
  Int_t**** fFEBV0LUTs; // look up tables [position][feb][LUT][entry] 

  // CB-{A/B/C}
  // 0 = B, 1 = A, 2 = C
  Int_t*** fCBLUTs; // control box look up tables

  // CB-A ----------------------------------------------------------------------
  TString fCBALUTequX; // logical equation used for generation of LUT X of CB-A
  TString fCBALUTequY; // logical equation used for generation of LUT Y of CB-A

  // CB-C ----------------------------------------------------------------------
  TString fCBCLUTequX; // logical equation used for generation of LUT X of CB-C
  TString fCBCLUTequY; // logical equation used for generation of LUT Y of CB-C

  // CBB -----------------------------------------------------------------------
  TString fCBBLUTequX; // logical equation used for generation of LUT X of CB-B 
  TString fCBBLUTequY; // logical equation used for generation of LUT Y of CB-B

  // CTP -----------------------------------------------------------------------
  // PT mask
  AliTRDptrgPTmasks fPTmasks; 
  // masks usage of internal signals for the pretrigger wake up signal
                              

  // CBB-LUT to TriggerInput assignment
  
  // class state ---------------------------------------------------------------

 private:
  AliTRDptrgParam();			     // instance only via Instance()
  AliTRDptrgParam(const AliTRDptrgParam &rhs); // not implemented
  AliTRDptrgParam& operator=(const AliTRDptrgParam &rhs); // not implemented

  ClassDef(AliTRDptrgParam, 1);
};
#endif
