#ifndef ALIMUONCHECK_H
#define ALIMUONCHECK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup evaluation
/// \class AliMUONCheck
/// \brief Class for data quality control
/// 
//  Author Frederic Yermia, INFN Torino

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliESDEvent;
class TTree;
class TH1F ;

class AliMUONCheck : public TObject
{
public:
  AliMUONCheck(const char* galiceFile, const char* esdFile,
               Int_t firstEvent=0, Int_t lastEvent=-1, const char* outDir="");
  AliMUONCheck(const char* galiceFile, const char* galiceFileSim, const char* esdFile,
               Int_t firstEvent=0, Int_t lastEvent=-1, const char* outDir="");
  virtual ~AliMUONCheck();
 
  void CheckESD(Bool_t pdc06TriggerResponse= false);
  void CheckKine();
  void CheckTrackRef();
  void CheckOccupancy(Bool_t perDetEle =kFALSE) const;  
  
  void SetEventsToCheck(Int_t firstEvent, Int_t lastEvent);
  void SetOutFileName(const TString& outFileName) { fOutFileName = outFileName; } 

private:
  /// Not implemented
  AliMUONCheck(const AliMUONCheck& rhs);
  /// Not implemented
  AliMUONCheck& operator=(const AliMUONCheck& rhs);
  
private:
  static const TString& GetDefaultOutFileName(); 

  TString fFileName;   //!<! File (galice.root) to read from fro reconstructed data
  TString fFileNameSim; //!<! File (galiceSim.root) for simulated data
  TString fesdFileName; //!<! File (AliESDs.root) to read from
 
  const char* fkOutDir;  //!<! output data directory
  TString fOutFileName;  //!<! output file name 
  
  Int_t   fFirstEvent;  //!<! First event to consider
  Int_t   fLastEvent;   //!<! Last event to consider

  ClassDef(AliMUONCheck,0) // Dumper of MUON related data
}; 

#endif
