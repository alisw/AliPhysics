// -*- mode: c++ -*-
#ifndef ALICALOCONSTANTS_H
#define ALICALOCONSTANTS_H

/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2010     *
 *                                                                        *
 * Primary Author: Per Thomas Hille  <perthomas.hille@yale.edu>           *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to   perthomas.hille@yale.edu                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


//
// Constants used by the HLT
// ALICE Offlinw
// and by EMCAL standalone debug tools
//
//
namespace CaloConstants
{
  const int MAXHOSTS=20; // related to the emcal debug online display
  const int TIMEBINS     = 256;       // number of sampling bins of the raw RO signal (we typically use 15-50; max is 1k+) 
  const double TIMEBINWITH = 100E-9 ;   // each sample is 100 ns
  const double TIMEBINMAX  =  TIMEBINS*TIMEBINWITH; 
  //  const double TAU = 2.35;
  //  const int  ORDER = 2;
  
  const int OVERFLOWCUT = 950;
 
  const double HGLGFACTOR = 16;
  
  // const double ECENTRALHIT = 0.85; //Percentage of total enegry contain in a single tower for a central hit 

  namespace ALTROConstants
  {
    const int ALTROMAXSAMPLES = 1008;    // The maximum number of samples of the ALTRO
    const int ALTROMAXPRESAMPLES = 15;   // Maximum number of presamles from the ALTRO chip     
    const int NALTROS        =   4;      // Number of ALTROs per frontend card
    const int NALTROCHANNELS =  16;      // Number of readout channels per ALTRO chip
    const int MINHARDWAREADDRESS = -2;   // Smallest possible HW address ( in offline )
    const int MAXHARDWAREADDRESS = 4096; // Max harware address,  ( its to high ) 
    const int MAXBINVALUE = 1023;        // Max possible ALTRO ADC value ( 10 bit )
    const int NGAINS         =   2;      // Number of gains ( high + low )
    const int HIGHGAIN    =   1;         // Mnemonic for High Gain
    const int LOWGAIN     =   0;         // Mnemonic for Low Gain
    const int HG = HIGHGAIN;             // Abbrevation for HIGHGAIN
    const int LG = LOWGAIN;              // Abbrevation for LOWGAIN
  }

  //FEE constants common to PHOS EMCAL
  const int CSPSPERFEE       =   32;    // Charge Sensitive Preamplifiers (CSPs) per FEE
  const int NBRANCHES        =    2;    // Branches per RCU   
  const int MAXHWADDRESSES   = 4096;    // Highest possible harware address
  
  namespace EMCALConstants
  {
    // const int NZROWSMOD      =  48;   // Number of rows per module
    // const int NXCOLUMNSMOD   =  24;   // Number of columns per module 
    const double ECENTRALHIT = 0.845678; //Percentage of total enegry contain in a single tower for a central hit  
    
    const int NZROWSMOD      =  24;   // Number of rows per module
    const int NXCOLUMNSMOD   =  48;   // Number of columns per module 
    
    const int NROWSMOD     = NZROWSMOD;   // Number of rows per module
    const int NCOLUMNSMOD  = NXCOLUMNSMOD;   // Number of columns per module 
    
    //   const int NZROWSMOD      =  24;   // Number of rows per module
    //   const int NXCOLUMNSMOD   =  48;   // Number of columns per module 
    
    const int NRCUSPERSECTOR = 4;     // Number of RCUs per sector
    const int NMODULES    =    10;    // Number of modules of the EMCAL detector
    const int NRCUSPERMODULE =  2 ;   // Number of RCUs per Module
    const int NFEECS         =  9;    // Number of Frontend cards per branch*/
    const int NZROWSRCU     =   48;   // Number of Rows per RCU
    const int NXCOLUMNSRCU  =   16;   // Number of columns per RCU
    const int ORDER  = 2; // Order of shaping stages of the signal conditioning unit
    const double TAU = 2.35;  // approximate shaping time
    
  }

  namespace PHOSConstants
  {
    const int NZROWSMOD      =  56;   // Number of rows per module       
    const int NXCOLUMNSMOD   =  64;   // Number of columns per module
    const int NMODULES    =     5;    // Number of modules of the PHOS detector
    const int NRCUSPERMODULE =  4 ;   // Number of RCUs per Module
    const int NFEECS         =  14;   // Number of Frontend cards per branch
  }


  //namespace FitAlgorithm
  // {
  //  enum fitAlgorithm { kStandard = 0, kCrude = 1, kPeakFinder = 2, kNeuralNet = 3, kFastFit= 4,
  //			kLogFit = 5, kLMS = 6,  kLMSOffline = 7, kFakeAltro = 9, kNONE = 8}; // possible return values
  // }
  
  namespace FitAlgorithm
  {
    enum fitAlgorithm { kStandard = 0, kCrude = 1, kPeakFinder = 2, kNeuralNet = 3, kFastFit= 4, kFakeAltro = 9, kNONE = 8}; // possible return values
  }

  
  
  namespace ReturnCodes
  {
    enum kReturnCode {kFitPar=1, kDummy=-1, kCrude=-9, kNoFit=-99, kInvalid=-9999};  // possible return values
  }

  namespace PeakFinderConstants
  {
    const int  MAXSTART = 3;      //  Start max 3 samples into the digitized array
    const int  SAMPLERANGE = 15;  //  Use maximum 15 samples for the Peak-Finder
  }
}

//For easier notation
namespace ALTRO =   CaloConstants::ALTROConstants;       // For easier notation
namespace Algo  =   CaloConstants::FitAlgorithm;         // For easier notation 
namespace Ret   =   CaloConstants::ReturnCodes;          // For easier notation
namespace PF    =   CaloConstants::PeakFinderConstants;  // For easier notation
namespace CALO  =   CaloConstants;                       // For easier notation
namespace EMCAL =   CaloConstants::EMCALConstants;       // For easier notation
namespace PHOS  =   CaloConstants::PHOSConstants;        // For easier notation


#endif

