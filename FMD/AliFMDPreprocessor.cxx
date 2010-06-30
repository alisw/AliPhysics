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
// See also 
//
//   http://aliceinfo.cern.ch/Offline/Activities/Shuttle.html
//
// Latest changes by Christian Holm Christensen
//

 #include <iostream>

#include <fstream>
#include "AliFMDPreprocessor.h"
#include "AliFMDCalibPedestal.h"
#include "AliFMDCalibGain.h"
#include "AliFMDCalibStripRange.h"
#include "AliFMDCalibSampleRate.h"
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
AliFMDPreprocessor::AliFMDPreprocessor(AliShuttleInterface* shuttle)
  : AliPreprocessor("FMD", shuttle)
{
  AddRunType("PHYSICS");
  AddRunType("STANDALONE");
  AddRunType("PEDESTAL");
  AddRunType("GAIN");
}


//____________________________________________________
Bool_t AliFMDPreprocessor::GetAndCheckFileSources(TList*&     list,
						  Int_t       system, 
						  const char* id) 
{
  // Convinience function 
  // Parameters: 
  //   list     On return, list of files. 
  //   system   Alice system (DAQ, DCS, ...)
  //   id       File id
  // Return:
  //   kTRUE on success. 
  list = GetFileSources(system, id);
  if (!list) { 
    TString sys;
    switch (system) { 
    case kDAQ: sys = "DAQ";     break;
    case kDCS: sys = "DCS";     break;
    default:   sys = "unknown"; break;
    }
    Log(Form("Failed to get file sources for %s/%s", sys.Data(), system));
    return kFALSE;
  }
  return kTRUE;
}

//____________________________________________________
AliCDBEntry* 
AliFMDPreprocessor::GetFromCDB(const char* second, const char* third)
{
  return GetFromOCDB(second, third);
}


//____________________________________________________
UInt_t AliFMDPreprocessor::Process(TMap* /* dcsAliasMap */)
{
  // Main member function. 
  // Parameters: 
  //    dcsAliassMap   Map of DCS data point aliases.
  // Return 
  //    0 on success, >0 otherwise 
  Bool_t resultPed   = kTRUE;
  Bool_t resultGain  = kTRUE;
  Bool_t resultRange = kTRUE;
  Bool_t resultRate  = kTRUE;
  Bool_t resultZero  = kTRUE;
  Bool_t infoCalib   = kTRUE;
  Bool_t resultDead  = kTRUE;
  // Do we need this ?
  // if(!dcsAliasMap) return 1;
  // 
  // Invoking the cdb manager and the FMD parameters class
  // AliCDBManager* cdb   = AliCDBManager::Instance();
  // cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  // cdb->SetRun(0);
  
  // Get the run type 
  TString runType(GetRunType()); 
  
  AliFMDParameters* pars = AliFMDParameters::Instance();
  if(runType.Contains("PEDESTAL", TString::kIgnoreCase))
    pars->Init(this, false, AliFMDParameters::kAltroMap|AliFMDParameters::kPulseGain);
  else if(runType.Contains("GAIN", TString::kIgnoreCase))
    pars->Init(this, false, AliFMDParameters::kAltroMap|AliFMDParameters::kPedestal);
  else
    pars->Init(this, false, AliFMDParameters::kAltroMap);
  
  // This is if the SOR contains Fee parameters, and we run a DA to
  // extract these parameters.   The same code could work if we get
  // the information from DCS via the FXS 
  TList* files = 0;
  AliFMDCalibSampleRate*      calibRate  = 0;
  AliFMDCalibStripRange*      calibRange = 0;
  AliFMDCalibZeroSuppression* calibZero  = 0;

  if (GetAndCheckFileSources(files, kDAQ,pars->GetConditionsShuttleID()))
    infoCalib = GetInfoCalibration(files, calibRate, calibRange, calibZero);
 
  resultRate  = (!calibRate  ? kFALSE : kTRUE);
  resultRange = (!calibRange ? kFALSE : kTRUE);
  resultZero  = (!calibZero  ? kFALSE : kTRUE);
  

  
  //Creating calibration objects
  AliFMDCalibPedestal* calibPed  = 0;
  AliFMDCalibGain*     calibGain = 0;
  AliFMDCalibDeadMap*  calibDead = 0;
  if (runType.Contains("PEDESTAL", TString::kIgnoreCase)) { 
    if (GetAndCheckFileSources(files, kDAQ, pars->GetPedestalShuttleID())) {
      if(files->GetSize())
	calibPed = GetPedestalCalibration(files);
    }
    resultPed = (calibPed ? kTRUE : kFALSE);
  }
  if (runType.Contains("GAIN", TString::kIgnoreCase)) {
    if (GetAndCheckFileSources(files, kDAQ, pars->GetGainShuttleID())) {
      if(files->GetSize())
	calibGain = GetGainCalibration(files);
    }
    resultGain = (calibGain ? kTRUE : kFALSE);
  }
  if(runType.Contains("PEDESTAL", TString::kIgnoreCase) || runType.Contains("GAIN", TString::kIgnoreCase))
    calibDead = GetDeadChannelMap(calibPed,calibGain);
  
  //Storing Calibration objects  
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Hans H. Dalsgaard");
  metaData.SetComment("Preprocessor stores pedestals and gains for the FMD.");
  
  if(calibPed)  { 
    resultPed  = Store("Calib","Pedestal", calibPed, &metaData, 0, kTRUE);
    delete calibPed;
  }
  if(calibGain) { 
    resultGain = Store("Calib","PulseGain", calibGain, &metaData, 0, kTRUE);
    delete calibGain;
  }
  if(calibRange) { 
    resultRange = Store("Calib","StripRange", calibRange, &metaData, 0, kTRUE);
    delete calibRange;
  }
  if(calibRate) { 
    resultRate = Store("Calib","SampleRate", calibRate, &metaData, 0, kTRUE);
    delete calibRate;
  }
  if(calibZero) { 
    resultZero = Store("Calib","ZeroSuppression", calibZero,&metaData,0,kTRUE);
    delete calibZero;
  }
  if(calibDead) { 
    resultDead = Store("Calib","Dead", calibDead,&metaData,0,kTRUE);
    delete calibDead;
  }

  Bool_t success = (resultPed && resultGain  && resultRange && 
		    resultRate  && resultZero && resultDead && infoCalib);
  
  Log(Form("FMD preprocessor was %s", (success ? "successful" : "failed")));
  return (success ? 0 : 1);
}

//____________________________________________________________________
Bool_t
AliFMDPreprocessor::GetInfoCalibration(TList* files, 
				       AliFMDCalibSampleRate*&      s,
				       AliFMDCalibStripRange*&      r, 
				       AliFMDCalibZeroSuppression*& z)
{
  // Get info calibrations. 
  // Parameters:
  //     files List of files. 
  //     s     On return, newly allocated object 
  //     r     On return, newly allocated object 
  //     z     On return, newly allocated object 
  // Return: 
  //     kTRUE on success
  if (!files) return kFALSE; // Should really be false
  if (files->GetEntries() <= 0) return kFALSE;
  
  s = new AliFMDCalibSampleRate();
  r = new AliFMDCalibStripRange();
  z = new AliFMDCalibZeroSuppression();
  
  AliFMDParameters*    pars     = AliFMDParameters::Instance();
  TIter                iter(files);
  TObjString*          fileSource;

  while((fileSource = dynamic_cast<TObjString*>(iter.Next()))) {
    const Char_t* filename = GetFile(kDAQ, pars->GetConditionsShuttleID(), fileSource->GetName());
    std::ifstream in(filename);
    if(!in) {
      Log(Form("File %s not found!", filename));
      continue;
    }
    s->ReadFromFile(in);
    r->ReadFromFile(in);
  }
  return kTRUE;
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
    const Char_t* filename = GetFile(kDAQ, pars->GetPedestalShuttleID(), fileSource->GetName());
    std::ifstream in(filename);
    if(!in) {
      Log(Form("File %s not found!", filename));
      continue;
    }
    // Loop until EOF
    int lineno = 0;
    char cc;
    while((cc = in.peek())!=EOF) {
      if(in.bad()) { 
	Log(Form("Bad read at line %d in %s", lineno, filename));
	break;
      }
      if (cc == '#') { 
	TString line;
	line.ReadLine(in);
	lineno++;
	if (lineno == 1) {
	  line.ToLower();
	  if(!line.Contains(pars->GetPedestalShuttleID())) {
	    Log(Form("File header is not from pedestal!: %s", line.Data()));
	    break;
	  }
	  Log("File contains data from pedestals");
	}
	continue;
      }
      UShort_t det, sec, strip;
      Char_t ring;
      Float_t ped, noise, mu, sigma, chi2ndf;
      Char_t c[8];
	  
      in >> det      >> c[0] 
	 >> ring     >> c[1]
	 >> sec      >> c[2]
	 >> strip    >> c[3]
	 >> ped      >> c[4]
	 >> noise    >> c[5]
	 >> mu       >> c[6]
	 >> sigma    >> c[7]
	 >> chi2ndf;
      lineno++;
      
      // Ignore trailing garbage 
      // if (strip > 127) continue;
      
      //Setting DDL to comply with the FMD in DAQ
      // UInt_t FmdDDLBase = 3072; 
      // ddl = ddl - FmdDDLBase;
      //Setting the pedestals via the hardware address
      
      
      // pars->Hardware2Detector(ddl,board,chip,channel,det,ring,sec,str);
      // strip += str;
     
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
    const Char_t* filename = GetFile(kDAQ, pars->GetGainShuttleID(), fileSource->GetName());
    std::ifstream in(filename);
    if(!in) {
      Log(Form("File %s not found!", filename));
      continue;
    }
    // Loop until EOF                                                                                                                                   
    int lineno = 0;
    char cc;
    while((cc = in.peek())!=EOF) {
      if(in.bad()) {
        Log(Form("Bad read at line %d in %s", lineno, filename));
        break;
      }
      if (cc == '#') {
        TString line;
        line.ReadLine(in);
        lineno++;
        if (lineno == 1) {
          line.ToLower();
          if(!line.Contains(pars->GetGainShuttleID())) {
            Log(Form("File header is not from gains!: %s", line.Data()));
            break;
          }
          Log("File contains data from gains");
        }
	continue;
      }
      UShort_t det, sec, strip;
      Char_t ring;
      
      Float_t gain,error,  chi2ndf;
      Char_t c[6];
      
      in >> det      >> c[0] 
	 >> ring     >> c[1]
	 >> sec      >> c[2]
	 >> strip    >> c[3]
	 >> gain     >> c[4]
	 >> error    >> c[5]
	 >> chi2ndf;
      lineno++;
      // Ignore trailing garbage
      //if(strip > 127) continue;
      
      //Setting DDL to comply with the FMD in DAQ
      // UInt_t FmdDDLBase = 3072; 
      // ddl = ddl - FmdDDLBase;
      //Setting the pedestals via the hardware address
      //   pars->Hardware2Detector(ddl,board,chip,channel,det,ring,sec,str);

      // strip += str;
      calibGain->Set(det,ring,sec,strip,gain);
    }
  }
  return calibGain;
}
//____________________________________________________________________
AliFMDCalibDeadMap*    
AliFMDPreprocessor::GetDeadChannelMap(AliFMDCalibPedestal* pedcalib,
				      AliFMDCalibGain*     gaincalib) {
  //creating dead channel map. '0' means 51200 entries
  AliFMDCalibDeadMap* deadmap = new AliFMDCalibDeadMap(0);
  //deadmap->Reset(kTRUE);
  Float_t noise = 0;
  Float_t gain  = 0;
  
  AliFMDParameters* pars = AliFMDParameters::Instance();
  //Looping over the channels.
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      
      for(UShort_t sec =0; sec < nsec;  sec++) {
	
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  
	  Bool_t isDead = kFALSE;
	  if(pedcalib)
	    noise = pedcalib->Width(det, ring, sec, strip);
	  else 
	    noise = pars->GetPedestalWidth(det, ring, sec, strip);
       
	  if(gaincalib)
	    gain  = gaincalib->Value(det, ring, sec, strip);
	  else 
	    gain  = pars->GetPulseGain(det, ring, sec, strip);
	  
	  //marking these channels dead.
	  if (gain < 0.5 || gain > 5 || noise > 10 || noise == 0) isDead = kTRUE;
	  
	  deadmap->operator()(det, ring, sec, strip) = isDead;
	}
      }
    }
  }

  return deadmap;
}
//____________________________________________________________________
//
// EOF
//
