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

class AliMUONData;
class AliRunLoader;
class AliLoader;
class AliESD;
class TTree;
class TH1F ;

class AliMUONCheck : public TObject
{
public:
  AliMUONCheck(const char* galiceFile, const char* esdFile,
               Int_t firstEvent=0, Int_t lastEvent=-1, const char* outDir="");
  virtual ~AliMUONCheck();
 
  /// Return true if contains valid data
  Bool_t IsValid() const { return (fData!=0); }
  
  void CheckESD(Bool_t pdc06TriggerResponse= false);
  void CheckKine();
  void CheckTrackRef();
  void CheckOccupancy(Bool_t perDetEle =kFALSE) const;  
  void CheckRecTracks() const;
  
  void SetEventsToCheck(Int_t firstEvent, Int_t lastEvent);

private:
  /// Not implemented
  AliMUONCheck(const AliMUONCheck& rhs);
  /// Not implemented
  AliMUONCheck& operator=(const AliMUONCheck& rhs);
  
private:
  TString fFileName;   //!< File (galice.root) to read from
  TString fesdFileName; //!< File (AliESDs.root) to read from
 
  const char* foutDir;  //!< output data directory
  
  Int_t   fFirstEvent;  //!< First event to consider
  Int_t   fLastEvent;   //!< Last event to consider

  AliRunLoader* fRunLoader; //!< AliRunLoader pointer
  AliLoader*    fLoader; //!< MUON loader pointer

  AliMUONData*  fData;  //!< AliMUONData pointer (to access containers)

  TTree   * fTree ;     //!< pointer to the analyzed TTree or TChain
  AliESD  * fESD ;      //!< Declaration of leave types

  ClassDef(AliMUONCheck,0) // Dumper of MUON related data
}; 

#endif
