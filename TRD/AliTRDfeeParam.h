#ifndef ALITRDFEEPARAM_H
#define ALITRDFEEPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD front end electronics parameters class                            //
//  Contains all FEE (MCM, TRAP, PASA) related                            //
//  parameters, constants, and mapping.                                   //
//                                                                        //
//  Author:                                                               //
//    Ken Oyama (oyama@physi.uni-heidelberg.de)                           //
//                                                                        //
//  many things now configured by AliTRDtrapConfig reflecting             //
//  the real memory structure of the TRAP (Jochen)                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class TRootIoCtor;

class AliTRDCommonParam;
class AliTRDpadPlane;
class AliTRDgeometry;

//_____________________________________________________________________________
class AliTRDfeeParam : public TObject
{

 public:

  AliTRDfeeParam(TRootIoCtor *);
  AliTRDfeeParam(const AliTRDfeeParam &p);
  virtual           ~AliTRDfeeParam();
  AliTRDfeeParam    &operator=(const AliTRDfeeParam &p);
  virtual void       Copy(TObject &p) const;

  static AliTRDfeeParam *Instance();  // Singleton
  static void            Terminate();

  // Translation from MCM to Pad and vice versa
  virtual Int_t    GetPadRowFromMCM(Int_t irob, Int_t imcm) const;
  virtual Int_t    GetPadColFromADC(Int_t irob, Int_t imcm, Int_t iadc) const;
  virtual Int_t    GetExtendedPadColFromADC(Int_t irob, Int_t imcm, Int_t iadc) const;
  virtual Int_t    GetMCMfromPad(Int_t irow, Int_t icol) const;
  virtual Int_t    GetMCMfromSharedPad(Int_t irow, Int_t icol) const;
  virtual Int_t    GetROBfromPad(Int_t irow, Int_t icol) const;
  virtual Int_t    GetROBfromSharedPad(Int_t irow, Int_t icol) const;
  virtual Int_t    GetRobSide(Int_t irob) const;
  virtual Int_t    GetColSide(Int_t icol) const;

  // SCSN-related
  static  UInt_t   AliToExtAli(Int_t rob, Int_t aliid);  // Converts the MCM-ROB combination to the extended MCM ALICE ID (used to address MCMs on the SCSN Bus)
  static  Int_t    ExtAliToAli( UInt_t dest, UShort_t linkpair, UShort_t rocType, Int_t *list, Int_t listSize);  // translates an extended MCM ALICE ID to a list of MCMs
  static  Short_t  ChipmaskToMCMlist( UInt_t cmA, UInt_t cmB, UShort_t linkpair, Int_t *mcmList, Int_t listSize );
  static  Short_t  GetRobAB( UShort_t robsel, UShort_t linkpair );  // Returns the chamber side (A=0, B=0) of a ROB

  // geometry
  static  Float_t  GetSamplingFrequency() { return (Float_t)fgkLHCfrequency / 4000000.0; }
  static  Int_t    GetNmcmRob()           { return fgkNmcmRob;      }
  static  Int_t    GetNmcmRobInRow()      { return fgkNmcmRobInRow; }
  static  Int_t    GetNmcmRobInCol()      { return fgkNmcmRobInCol; }
  static  Int_t    GetNrobC0()            { return fgkNrobC0;       }
  static  Int_t    GetNrobC1()            { return fgkNrobC1;       }
  static  Int_t    GetNadcMcm()           { return fgkNadcMcm;      }
  static  Int_t    GetNcol()              { return fgkNcol;         }
  static  Int_t    GetNcolMcm()           { return fgkNcolMcm;      }
  static  Int_t    GetNrowC0()            { return fgkNrowC0;       }
  static  Int_t    GetNrowC1()            { return fgkNrowC1;       }

  // tracklet simulation
	  Bool_t   GetTracklet()         const { return fgTracklet; } 
  static  void     SetTracklet(Bool_t trackletSim = kTRUE) { fgTracklet = trackletSim; }

  // Concerning raw data format
          Int_t    GetRAWversion() const                    { return fRAWversion;        }
          void     SetRAWversion( Int_t rawver );

 protected:

  static AliTRDfeeParam *fgInstance;         // Singleton instance
  static Bool_t          fgTerminated;       // Defines if this class has already been terminated

  AliTRDCommonParam     *fCP;                // TRD common parameters class

  // Remark: ISO C++ allows initialization of static const values only for integer.

  // Basic Geometrical numbers
  static const Int_t    fgkLHCfrequency      = 40079000 ; // [Hz] LHC clock (should be moved to STEER?)
  static const Int_t    fgkNmcmRob           = 16;        // Number of MCMs per ROB         (old fgkMCMmax)
  static const Int_t    fgkNmcmRobInRow      = 4;         // Number of MCMs per ROB in row dir. (old fgkMCMrow)
  static const Int_t    fgkNmcmRobInCol      = 4;         // Number of MCMs per ROB in col dir. (old fgkMCMrow)
  static const Int_t    fgkNrobC0            = 6;         // Number of ROBs per C0 chamber  (old fgkROBmaxC0)
  static const Int_t    fgkNrobC1            = 8;         // Number of ROBs per C1 chamber  (old fgkROBmaxC1)
  static const Int_t    fgkNadcMcm           = 21;        // Number of ADC channels per MCM (old fgkADCmax)
  static const Int_t    fgkNcol              = 144;       // Number of pads per padplane row(old fgkColmax)
  static const Int_t    fgkNcolMcm           = 18;        // Number of pads per MCM         (old fgkPadmax)
  static const Int_t    fgkNrowC0            = 12;        // Number of Rows per C0 chamber  (old fgkRowmaxC0)
  static const Int_t    fgkNrowC1            = 16;        // Number of Rows per C1 chamber  (old fgkRowmaxC1)

 // Tracklet  processing on/off 
  static       Bool_t   fgTracklet; // tracklet processing

  // For raw production
               Int_t    fRAWversion;                      // Raw data production version
  static const Int_t    fgkMaxRAWversion      = 3;        // Maximum raw version number supported

 private:

  AliTRDfeeParam();

  ClassDef(AliTRDfeeParam,4)                              // The TRD front end electronics parameter

};
#endif

