/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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
/* $Id$ */
/** @file    AliFMDPreprocessor.cxx
    @author  Hans Hjersing Dalsgaard <canute@nbi.dk>
    @date    Mon Mar 27 12:39:09 2006
    @brief   Shuttle "preprocessor" for the FMD
*/
//___________________________________________________________________
//
// The class processes data points from DCS (via Amanada), and DAQ DA
// files (via FXS) to make calibration data for the FMD. 
//
// Data points: 
//    *  Nothing yet. 
//
// DAQ FXS file:
//    * pedestals - a (ASCII) Comma Separated Values files with the
//                  fields 
//                       rcu	 DDL number 
//                       board   FEC board number 
//                       chip    ALTRO chip number on FEC
//                       channel ALTRO channel number
//                       strip   VA1 strip number
//                       sample  Sample number
//                       ped     Mean of ADC spectra
//                       noise   Spread of ADC spectra
//                       mu      Mean of Gaussian fit to ADC spectra
//                       sigma   Variance of Gaussian fit to ADC spectra
//                       chi2    Chi^2 per degrees of freedom of fit
//    * Gains     - a (ASCII) Comma Separated Values files with the
//                  fields 
//                       rcu	 DDL number 
//                       board   FEC board number 
//                       chip    ALTRO chip number on FEC
//                       channel ALTRO channel number
//                       strip   VA1 strip number
//                       gain    Slope of gain
//                       error   Error on gain
//                       chi2    Chi^2 per degrees of freedom of fit
//                  
// Latest changes by Christian Holm Christensen
//

// #include <iostream>

#include <fstream>
#include "AliFMDPreprocessor.h"
#include "AliFMDCalibPedestal.h"
#include "AliFMDCalibGain.h"
#include "AliFMDParameters.h"
#include "AliCDBMetaData.h"
#include "AliCDBManager.h"
// #include "AliDCSValue.h"
#include "AliLog.h"
#include <TTimeStamp.h>
// #include <TFile.h>
#include <TObjString.h>
#include <TString.h>
#include <TNamed.h>


ClassImp(AliFMDPreprocessor)
#if 0 // Do not remove - here to make Emacs happy
;
#endif 

//____________________________________________________
UInt_t AliFMDPreprocessor::Process(TMap* /* dcsAliasMap */)
{
  // Main member function. 
  // Parameters: 
  //    dcsAliassMap   Map of DCS data point aliases.
  // Return 
  //    ? 

  // Do we need this ?
  // if(!dcsAliasMap) return 1;
  // 
  // Invoking the cdb manager and the FMD parameters class
  // AliCDBManager* cdb   = AliCDBManager::Instance();
  // cdb->SetDefaultStorage("local://$ALICE_ROOT");
  // cdb->SetRun(0);
  AliFMDParameters* pars = AliFMDParameters::Instance();
  pars->Init(false, AliFMDParameters::kAltroMap);
  
  //Creating calibration objects
  TList*               pedFiles     = GetFileSources(kDAQ,"pedestal");
  TList*               gainFiles    = GetFileSources(kDAQ, "gain");
  AliFMDCalibPedestal* calibPed     = GetPedestalCalibration(pedFiles);
  AliFMDCalibGain*     calibGain    = GetGainCalibration(gainFiles);
  
  
  //Storing Calibration objects  
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Hans H. Dalsgaard");
  metaData.SetComment("Preprocessor stores pedestals and gains for the FMD.");
  
  Bool_t resultPed = kFALSE, resultGain = kFALSE;
  if(calibPed)  resultPed  = Store("Calib","Pedestal", calibPed, &metaData);
  if(calibGain) resultGain = Store("Calib","PulseGain", calibGain, &metaData);
  if (calibPed)  delete calibPed;
  if (calibGain) delete calibGain;
  
  return (resultPed && resultGain ? 0 : 1);
}

//____________________________________________________________________
AliFMDCalibPedestal* 
AliFMDPreprocessor::GetPedestalCalibration(TList* pedFiles)
{
  // Read DAQ DA produced CSV files of pedestals, and return a
  // calibration object. 
  // Parameters:
  //   pedFiles     List of pedestal files 
  // Return 
  //   A pointer to a newly allocated AliFMDCalibPedestal object, or
  //   null in case of errors. 
  if(!pedFiles) return 0;

  AliFMDCalibPedestal* calibPed = new AliFMDCalibPedestal();
  AliFMDParameters*    pars     = AliFMDParameters::Instance();
  TIter                iter(pedFiles);
  TObjString*          fileSource;
  
  while((fileSource = dynamic_cast<TObjString*>(iter.Next()))) {
    const Char_t* filename = GetFile(kDAQ, "pedestal", fileSource->GetName());
    std::ifstream in(filename);
    if(!in) {
      AliError(Form("File %s not found!", filename));
      continue;
    }

    // Get header (how long is it ?)
    TString header;
    header.ReadLine(in);
    header.ToLower();
    if(!header.Contains("pedestal")) {
      AliError("File header is not from pedestal!");
      continue;
    }
    Log("File contains data from pedestals");
    
    // Read columns line
    int lineno = 2;
    header.ReadLine(in);
    
    // Loop until EOF
    while(!in.eof()) {
      if(in.bad()) { 
	AliError(Form("Bad read at line %d in %s", lineno, filename));
	break;
      }
      UInt_t ddl=2, board, chip, channel, strip, tb;
      Float_t ped, noise, mu, sigma, chi2ndf;
      Char_t c[10];
	  
      in // >> ddl      >> c[0] 
	 >> board    >> c[1]
	 >> chip     >> c[2]
	 >> channel  >> c[3]
	 >> strip    >> c[4]
	 >> tb       >> c[5]
	 >> ped      >> c[6]
	 >> noise    >> c[7]
	 >> mu       >> c[8]
	 >> sigma    >> c[9]
	 >> chi2ndf;
      lineno++;
      // Ignore trailing garbage 
      if (strip > 127) continue;
      
      //Setting the pedestals via the hardware address
      UShort_t det, sec, str;
      Char_t ring;
	  
      pars->Hardware2Detector(ddl,board,chip,channel,det,ring,sec,str);
      strip += str;
      calibPed->Set(det,ring,sec,strip,ped,noise);
    }
  }
  return calibPed;
}	

//____________________________________________________________________
AliFMDCalibGain* 
AliFMDPreprocessor::GetGainCalibration(TList* gainFiles)
{
  // Read DAQ DA produced CSV files of pedestals, and return a
  // calibration object. 
  // Parameters:
  //   pedFiles     List of pedestal files 
  // Return 
  //   A pointer to a newly allocated AliFMDCalibPedestal object, or
  //   null in case of errors. 
  if(!gainFiles) return 0;
  
  AliFMDCalibGain*  calibGain  = new AliFMDCalibGain();
  AliFMDParameters* pars       = AliFMDParameters::Instance();
  TIter             iter(gainFiles);
  TObjString*       fileSource;
  while((fileSource = dynamic_cast<TObjString *>(iter.Next()))) {
    const Char_t* filename = GetFile(kDAQ, "gain", fileSource->GetName());
    std::ifstream in(filename);
    if(!in) {
      AliError(Form("File %s not found!", filename));
      continue;
    }

    //Get header (how long is it ?)
    TString header;
    header.ReadLine(in);
    header.ToLower();
    if(!header.Contains("gain")) {
      AliError("File header is not from gain!");
      continue;
    }
    Log("File contains data from pulse gain");

    // Read column headers
    header.ReadLine(in);

    int lineno  = 2;
    // Read until EOF 
    while(!in.eof()) {
      if(in.bad()) { 
	AliError(Form("Bad read at line %d in %s", lineno, filename));
	break;
      }
      UInt_t ddl=2, board, chip, channel, strip;
      Float_t gain,error,  chi2ndf;
      Char_t c[7];
	      
      in // >> ddl      >> c[0] 
	 >> board    >> c[1]
	 >> chip     >> c[2]
	 >> channel  >> c[3]
	 >> strip    >> c[4]
	 >> gain     >> c[5]
	 >> error    >> c[6]
	 >> chi2ndf;
      lineno++;
      // Ignore trailing garbage
      if(strip > 127) continue;
      
      //Setting the pedestals via the hardware address
      UShort_t det, sec, str;
      Char_t ring;
      pars->Hardware2Detector(ddl,board,chip,channel,det,ring,sec,str);

      strip += str;
      calibGain->Set(det,ring,sec,strip,gain);
    }
  }
  return calibGain;
}

//____________________________________________________________________
//
// EOF
//
