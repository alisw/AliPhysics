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
#include "AliTRDCalOnlineGainTable.h"
#include "AliTRDtrapConfig.h"

class AliTRDtrapConfigHandler : public TObject {
 public:
                    AliTRDtrapConfigHandler(AliTRDtrapConfig *cfg);
  virtual          ~AliTRDtrapConfigHandler();

  void Init();                                                // Set DMEM allocation modes
  void ResetMCMs();                                           // Reset all trap registers and DMEM of the MCMs
  Int_t LoadConfig();                                         // load a default configuration suitable for simulation
  Int_t LoadConfig(TString filename);                         // load a TRAP configuration from a file
  Int_t SetGaintable(AliTRDCalOnlineGainTable const &gtbl);   // Set a gain table to correct Q0 and Q1 for PID

  void ProcessLTUparam(Int_t dest, Int_t addr, UInt_t data);  // Process the LTU parameters
  void PrintGeoTest();                                        // Prints some information about the geometry. Only for debugging

  // UInt_t peek(Int_t rob, Int_t mcm, Int_t addr);   // not implemented yet
  // Int_t poke(Int_t rob, Int_t mcm, Int_t addr, UInt_t value);   // not implemented yet

 private:
  Bool_t AddValues(UInt_t det, UInt_t cmd, UInt_t extali, Int_t addr, UInt_t data);

  void  ConfigureDyCorr(Int_t det);                             // deflection length correction due to Lorentz angle and tilted pad correction
  void  ConfigureDRange(Int_t det);                             // deflection range LUT,  range calculated according to B-field (in T) and pt_min (in GeV/c)
  void  ConfigureNTimebins(Int_t det);                          // timebins in the drift region
  void  ConfigurePIDcorr(Int_t det);                            // Calculate the mcm individual correction factors for the PID

  Double_t Square(Double_t val) const { return val*val; };          // returns the square of a given number

  AliTRDtrapConfigHandler(const AliTRDtrapConfigHandler &h);             // not implemented
  AliTRDtrapConfigHandler &operator=(const AliTRDtrapConfigHandler &h);  // not implemented

  static const UInt_t fgkScsnCmdReset=6;     // SCSN command for reset
  static const UInt_t fgkScsnCmdPause=8;     // SCSN command to pause
  static const UInt_t fgkScsnCmdRead=9;      // SCSN command to read
  static const UInt_t fgkScsnCmdWrite=10;    // SCSN command to write
  static const UInt_t fgkScsnCmdPtrg=12;     // SCSN command for pretrigger
  static const UInt_t fgkScsnCmdRobPower=16; // SCSN command to switch ROB power
  static const UInt_t fgkScsnCmdRobReset=17; // SCSN command for ROB reset

  static const UInt_t fgkScsnCmdRestr=18;    // SCSN command to restrict commands to specified chambers
  static const UInt_t fgkScsnCmdTtcRx=19;    // SCSN command to configure TTCrx
  static const UInt_t fgkScsnCmdHwPtrg=20;   // SCSN command to issue pretrigger pulse
  static const UInt_t fgkScsnCmdSetHC=22;    // SCSN command to set HC ID
  static const UInt_t fgkScsnCmdMcmTemp=24;  // SCSN command for MCM temperature sensors
  static const UInt_t fgkScsnCmdPM=25;       // SCSN command for patchmaker
  static const UInt_t fgkScsnCmdOri=26;      // SCSN command for ORI configuration
  static const UInt_t fgkScsnLTUparam=27;    // extended SCSN command for the LTU configuration

  static const Int_t fgkMCMperROBCol = 4;  // MCMs per ROB column
  static const Int_t fgkPadsPerMCM = 18;   // readout pads per MCM
  static const Int_t fgkMCMperROBRow = 4;  // MCMs per ROB row

  static const Int_t fgkMaxLinkPairs=4;    // number of linkpairs used during configuration
  static const Int_t fgkMcmlistSize=256;     // list of MCMs to which a value has to be written

  AliTRDltuParam     ltuParam;             // ltuParam class for the actual calculation of the parameters

  UInt_t fRestrictiveMask;                 // mask to restrict subsequent commands to specified chambers

  AliTRDtrapConfig *fTrapConfig;           // pointer to TRAP config in use
  AliTRDCalOnlineGainTable fGtbl;          // gain table

  ClassDef(AliTRDtrapConfigHandler,0)
};

#endif
