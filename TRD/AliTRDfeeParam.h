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

  static  UInt_t   AliToExtAli(Int_t rob, Int_t aliid);  // Converts the MCM-ROB combination to the extended MCM ALICE ID (used to address MCMs on the SCSN Bus)
  static  Int_t    ExtAliToAli( UInt_t dest, UShort_t linkpair, UShort_t rocType, Int_t *list, Int_t listSize);  // translates an extended MCM ALICE ID to a list of MCMs
  static  Short_t  ChipmaskToMCMlist( UInt_t cmA, UInt_t cmB, UShort_t linkpair, Int_t *mcmList, Int_t listSize );
  static  Short_t  GetRobAB( UShort_t robsel, UShort_t linkpair );  // Returns the chamber side (A=0, B=0) of a ROB

  static  Float_t  GetSamplingFrequency() { return (Float_t)fgkLHCfrequency / 4000000.0; }
  static  Int_t    GetNmcmRob()           { return fgkNmcmRob;      }
  static  Int_t    GetNmcmRobInRow()      { return fgkNmcmRobInRow; }
  static  Int_t    GetNmcmRobInCol()      { return fgkNmcmRobInCol; }
  static  Int_t    GetNrobC0()            { return fgkNrobC0;       }
  static  Int_t    GetNrobC1()            { return fgkNrobC1;       }
  static  Int_t    GetNadcMcm()           { return fgkNadcMcm;      }
  // static  Int_t    GetNtimebin()       { return fgkNtimebin;     }
  static  Int_t    GetNcol()              { return fgkNcol;         }
  static  Int_t    GetNcolMcm()           { return fgkNcolMcm;      }
  static  Int_t    GetNrowC0()            { return fgkNrowC0;       }
  static  Int_t    GetNrowC1()            { return fgkNrowC1;       }

  // static  Int_t    GetADCpedestal()    { return fgkADCpedestal;  }
  // static  Int_t    GetADCnoise()       { return fgkADCnoise;     }
  static  Int_t    GetADCDAC()            { return fgkADCDAC;       }

  static  Bool_t   IsPFon()               { return fgkPFon;         }
  static  Bool_t   IsGFon()               { return fgkGFon;         }
  static  Bool_t   IsTFon()               { return fgkTFon;         }

  static  Int_t    GetPFtimeConstant()    { return fgkPFtimeConstant;   }
  static  Int_t    GetPFeffectPedestal()  { return fgkPFeffectPedestal; }

  //new
  static  Int_t    GetQacc0Start()        {  return fgkPREPqAcc0Start; }
  static  Int_t    GetQacc0End()          {  return fgkPREPqAcc0End; }
  static  Int_t    GetQacc1Start()        {  return fgkPREPqAcc1Start; }
  static  Int_t    GetQacc1End()          {  return fgkPREPqAcc1End; }
          Float_t  GetMinClusterCharge() const {  return fgkMinClusterCharge; }
  static  Int_t    GetLinearFitStart()    {  return fgkPREPLinearFitStart; }
  static  Int_t    GetLinearFitEnd()      {  return fgkPREPLinearFitEnd;  }

  //        Float_t  GetClusThr()           { return fClusThr; };
  //        Float_t  GetPadThr() const { return fPadThr; };
  //        Int_t    GetTailCancelation() const { return fTCOn; };
  //        Int_t    GetNexponential() const { return fTCnexp; };
  //virtual void     GetFilterParam(Float_t &r1, Float_t &r2, Float_t &c1, Float_t &c2, Float_t &ped) const;
  //        Int_t    GetFilterType() const { return fFilterType; };

  static  Int_t    GetTFtype()            { return fgkTFtype;       }
  //static  Int_t    GetTFnExp()            { return fgkTFnExp;       }
          Int_t    GetTFnExp()      const { return fTFnExp;         }
          Float_t  GetTFr1()        const { return fTFr1;           }
          Float_t  GetTFr2()        const { return fTFr2;           }
          Float_t  GetTFc1()        const { return fTFc1;           }
          Float_t  GetTFc2()        const { return fTFc2;           }

 // for tracklets
	  Bool_t   GetTracklet()         const { return fgTracklet; } 
  static  void     SetTracklet(Bool_t trackletSim = kTRUE) { fgTracklet = trackletSim; }
          Int_t    GetMaxNrOfTracklets() const { return fgkMaxNrOfTracklets; } 
	  Bool_t   GetMCTrackletOutput() const { return fgkMCTrackletOutput; }

  static  Float_t  GetTFattPar()          { return ((Float_t) fgkTFattPar1) / ((Float_t) fgkTFattPar2); }
          Float_t  GetTFf0()        const { return 1.0 + fgkTFon*(-1.0+GetTFattPar()); }   // 1 if TC off

          void     SetEBsglIndThr(Int_t val);  
          Int_t    GetEBsglIndThr() const { return fEBsglIndThr;    }

          void     SetEBsumIndThr(Int_t val);
          Int_t    GetEBsumIndThr() const { return fEBsumIndThr;    }

          void     SetEBindLUT(Int_t val);
          Int_t    GetEBindLUT()    const { return fEBindLUT;       }

          void     SetEBignoreNeighbour(Int_t val);
          Int_t    GetEBignoreNeighbour() const             { return fEBignoreNeighbour; }

  // Concerning raw data format
          Int_t    GetRAWversion() const                    { return fRAWversion;        }
          void     SetRAWversion( Int_t rawver );
          Bool_t   GetRAWstoreRaw() const                   { return fRAWstoreRaw;       }
          void     SetRAWstoreRaw( Bool_t storeraw )        { fRAWstoreRaw = storeraw;   }

          void     SetXenon();
          void     SetArgon();

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
  // static const Int_t    fgkNtimebin       = 24;        // Number of Time bins should come from calibDB
  static const Int_t    fgkNcol              = 144;       // Number of pads per padplane row(old fgkColmax)
  static const Int_t    fgkNcolMcm           = 18;        // Number of pads per MCM         (old fgkPadmax)
  static const Int_t    fgkNrowC0            = 12;        // Number of Rows per C0 chamber  (old fgkRowmaxC0)
  static const Int_t    fgkNrowC1            = 16;        // Number of Rows per C1 chamber  (old fgkRowmaxC1)

  // ADC intrinsic parameters
  static const Int_t    fgkADCDAC            = 0;         // 5 bit ADC gain parameter

  // TRAP filter global setup
  static const Bool_t   fgkPFon              = kTRUE;     // Pedestal Filter enable/disable flag.
  static const Bool_t   fgkGFon              = kFALSE;    // Gain correction Filter enable/disable flag
  static const Bool_t   fgkTFon              = kTRUE;     // Tail cancelation Filter enable/disable flag (old name fTCOn)

  // PF setup
  static const Int_t    fgkPFtimeConstant    =  0;        // 0 for fastest, 3 for slowest (no effect, probably)
  static const Int_t    fgkPFeffectPedestal  = 10;        // [in ADC units] the desired baseline (Additive)

  // GF setup
  static const Int_t    fgkGFnoise           =  0;        // Noise level increased by gain filter x 100 [in ADC] (to be measured)

  // TF setup
  static const Int_t    fgkTFtype            = 1;         // TC type (0=analog, 1=digital, 2=MI, 3=close to electronics) (old name fFilterType)

  // OLD TF setup (calculated from above)  (valid only for fgkTFsimType = 0 or 1)
  //static const Int_t    fgkTFnExp          = 1;           // Number of exponential for simType 0 and 1
               Int_t    fTFnExp;                            // Number of exponential for simType 0 and 1

 // Tracklet  processing on/off 
  static       Bool_t   fgTracklet; // tracklet processing

  static const Int_t    fgkMaxNrOfTracklets = 4;          // Max. nr of tracklet words for one mcm

  // additional tracklet folder structure output, 
  // containing all necessary Monte Carlo information; maybe this should go somewhere else;
  static const Bool_t   fgkMCTrackletOutput = kTRUE;      // Default should be kTRUE

  // following need Instance because initialized in constructor
               Float_t  fTFr1;                            // Time constant [us] long (old name fR1)
               Float_t  fTFr2;                            // Time constant [us] short(old name fR2)
               Float_t  fTFc1;                            // Weight long  (old name fC1)
               Float_t  fTFc2;                            // Weight short (old name fC2)

  // here is for TRAP simulation (not yet used)
  static const Int_t    fgkTFdecayWeightL     = 270;      // 0 to 1024 corresponds to 0 to 0.5
  static const Int_t    fgkTFdecayParL        = 348;      // 0 to 511 corresponds to 0.75 to 1
  static const Int_t    fgkTFdecayParS        = 449;      // 0 to 511 correponds to 0.25 to 0.5
  static const Int_t    fgkTFattPar1          = 45;       // attenuationParameter = fgkTFattenuationParameter1/fgkTFattenuationParameter2
  static const Int_t    fgkTFattPar2          = 14;       //                      = -alphaL/ln(lambdaL)-(1-alphaL)/ln(lambdaS)

  // ZS parameters
               Int_t    fEBsglIndThr;                     // EBIS in ADC units
               Int_t    fEBsumIndThr;                     // EBIT in ADC units
               Int_t    fEBindLUT;                        // EBIL lookup table
               Int_t    fEBignoreNeighbour;               // EBIN 0:include neighbor

  // Charge accumulators
  static const Int_t    fgkPREPqAcc0Start     =  5;       // Preprocessor Charge Accumulator 0 Start
  static const Int_t    fgkPREPqAcc0End       = 10;       // Preprocessor Charge Accumulator 0 End
  static const Int_t    fgkPREPqAcc1Start     = 11;       // Preprocessor Charge Accumulator 1 Start
  static const Int_t    fgkPREPqAcc1End       = 20;       // Preprocessor Charge Accumulator 1 End
  static const Int_t    fgkMinClusterCharge   = 20;       // Hit detection [in ADC units]

  //new
  static const Int_t    fgkPREPLinearFitStart = 5;        // Time constants for linear fit
  static const Int_t    fgkPREPLinearFitEnd   = 20;       // Time constants for linear fit

  // OLD TRAP processing parameters calculated from above
  //static const Float_t  fClusThr;                       // Cluster threshold
  //static const Float_t  fPadThr;                        // Pad threshold

  // For raw production
               Int_t    fRAWversion;                      // Raw data production version
  static const Int_t    fgkMaxRAWversion      = 3;        // Maximum raw version number supported
               Bool_t   fRAWstoreRaw;                     // Store unfiltered data for raw data stream

 private:

  AliTRDfeeParam();

  ClassDef(AliTRDfeeParam,3)                              // The TRD front end electronics parameter

};
#endif

