#ifndef ALITRDMCMSIMCONFIGHANDLER_H
#define ALITRDMCMSIMCONFIGHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////
//                                                            //
//  Multi Chip Module Simulation Configuration Handler Class  //
//                                                            //
////////////////////////////////////////////////////////////////


#include <TObject.h>
#include "AliTRDltuParam.h"

class AliTRDtrapConfig;

class AliTRDtrapConfigHandler : public TObject {
 public:
                    AliTRDtrapConfigHandler();
  virtual          ~AliTRDtrapConfigHandler();

  void ResetMCMs();                                           // Reset all trap registers and DMEM of the MCMs
  Int_t LoadConfig();                                         // load a default configuration suitable for simulation
  Int_t LoadConfig(TString filename);                         // load a TRAP configuration from a file

  void ProcessLTUparam(Int_t dest, Int_t addr, UInt_t data);  // Process the LTU parameters
  void PrintGeoTest();                                        // Prints some information about the geometry. Only for debugging

  // UInt_t peek(Int_t rob, Int_t mcm, Int_t addr);   // not implemented yet
  // Int_t poke(Int_t rob, Int_t mcm, Int_t addr, UInt_t value);   // not implemented yet

 private:
  void  ConfigureDyCorr(Int_t det);                             // deflection length correction due to Lorentz angle and tilted pad correction
  void  ConfigureDRange(Int_t det);                             // deflection range LUT,  range calculated according to B-field (in T) and pt_min (in GeV/c)
  void  ConfigureNTimebins(Int_t det);                          // timebins in the drift region
  void  ConfigurePIDcorr(Int_t det);                            // Calculate the mcm individual correction factors for the PID

  Double_t Square(Double_t val) { return val*val; };          // returns the square of a given number

  AliTRDtrapConfigHandler(const AliTRDtrapConfigHandler &h);             // not implemented
  AliTRDtrapConfigHandler &operator=(const AliTRDtrapConfigHandler &h);  // not implemented

  static const UInt_t fgkScsnCmdWrite=10;  // SCSN command for the write command
  static const UInt_t fgkScsnCmdRestr=18;  // SCSN command to restrict commands to specified chambers
  static const UInt_t fgkScsnLTUparam=27;  // extended SCSN command for the LTU configuration

  static const Int_t fgkMCMperROBCol = 4;  // MCMs per ROB column
  static const Int_t fgkPadsPerMCM = 18;   // readout pads per MCM
  static const Int_t fgkMCMperROBRow = 4;  // MCMs per ROB row

  AliTRDltuParam     ltuParam;             // ltuParam class for the actual calculation of the parameters

  UInt_t fRestrictiveMask;                 // mask to restrict subsequent commands to specified chambers

  ClassDef(AliTRDtrapConfigHandler,0)
};

#endif
