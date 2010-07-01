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

class AliTRDgeometry;
class AliTRDtrapConfig;
class AliTRDpadPlane;


class AliTRDtrapConfigHandler : public TObject {
 public:
                    AliTRDtrapConfigHandler();
  virtual          ~AliTRDtrapConfigHandler();

  void ResetMCMs();                                           // Reset all trap registers and DMEM of the MCMs
  Int_t LoadConfig(TString filename, Int_t det);              // load a TRAP configuration from a file

  void ProcessLTUparam(Int_t dest, Int_t addr, UInt_t data);  // Process the LTU parameters
  void PrintGeoTest();                                        // Prints some information about the geometry. Only for debugging

  // UInt_t peek(Int_t rob, Int_t mcm, Int_t addr);   // not implemented yet
  // Int_t poke(Int_t rob, Int_t mcm, Int_t addr, UInt_t value);   // not implemented yet

 private:
  void  ConfigureDyCorr();                                    // deflection length correction due to Lorentz angle and tilted pad correction
  void  ConfigureDRange();                                    // deflection range LUT,  range calculated according to B-field (in T) and pt_min (in GeV/c)
  void  ConfigureNTimebins();                                 // timebins in the drift region
  void  ConfigurePIDcorr();                                   // Calculate the mcm individual correction factors for the PID

  Int_t GetPadPosNonRot(Int_t rob, Int_t mcm, Int_t channel, Double_t trackCoor[3]);    // calcutate the gobal coordinates for an mcm channel in the supermodule at position -0.5
  void GetLocalPadPos(AliTRDpadPlane *plane, Int_t rob, Int_t mcm, Int_t channel, Double_t result[2]); // calculate the local coordinates for an mcm channel

  Double_t Square(Double_t val);  // returns the square of a given number

  AliTRDtrapConfigHandler(const AliTRDtrapConfigHandler &h);             // not implemented
  AliTRDtrapConfigHandler &operator=(const AliTRDtrapConfigHandler &h);  // not implemented



  static const UInt_t fgkScsnCmdWrite=10;  // SCSN command for the write command
  static const UInt_t fgkScsnLTUparam=27;  // extended SCSN command for the LTU configuration

  static const Int_t fgkDyMaxCut = 63;     // Maximum value of the deflection cut
  static const Int_t fgkDyMinCut = -64;    // Minimum value of the deflection cut

  static const Int_t fgkMCMperROBCol = 4;  // MCMs per ROB column
  static const Int_t fgkPadsPerMCM = 18;   // readout pads per MCM
  static const Int_t fgkMCMperROBRow = 4;  // MCMs per ROB row


  AliTRDgeometry *fGeo;                    // Pointer to the AliTRDgeometry class

  Int_t fDet;   // detector number (0 - 539)

  Double_t fBField;                        // value of the L3 magnet field
  Double_t fOmegaTau;                      // ometa tau
  Double_t fPtMin;                         // lower p_t threshold for the tracks which should pass the deflection cut
  Int_t fNTimebins;                        // Number of time bins in the drift region (only relevant for GTU)
  UInt_t fScaleQ0;                         // scale parameter to map the accumulated charge in the first time window to a memory address
  UInt_t fScaleQ1;                         // scale parameter to map the accumulated charge in the second time window to a memory address
  Bool_t fPidTracklengthCorr;              // Factor to correct the accumulated charge for track length effects
  Bool_t fTiltCorr;                        // tilting correction


  ClassDef(AliTRDtrapConfigHandler,0)
};


#endif

