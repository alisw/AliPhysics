#ifndef ALIMUONMCHVIEWAPPLICATION_H
#define ALIMUONMCHVIEWAPPLICATION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONMchViewApplication
/// \brief Main class for the mchview program
///
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TRint
#   include <TRint.h>
#endif
#ifndef ROOT_RQ_OBJECT
#   include <RQ_OBJECT.h>
#endif

class TGMainFrame;

class AliMUONMchViewApplication : public TRint
{
  RQ_OBJECT("AliMUONMchViewApplication")

public:
  AliMUONMchViewApplication(const char* name, int* argc, char** argv, 
                            Float_t wfraction, Float_t hfraction);
  virtual ~AliMUONMchViewApplication();

  void HandleMenu(Int_t i);

  /// Return the version number of the mchview application
  static const char* Version() { return "0.9"; }
  
private:
  /// Not implemented
  AliMUONMchViewApplication(const AliMUONMchViewApplication& rhs);
  /// Not implemented
  AliMUONMchViewApplication& operator=(const AliMUONMchViewApplication& rhs);
  
  void CreateMenuBar(UInt_t w);
  void Save();
  void Save(const char* filename);
  void Open();
  void Open(const char* filename);
  
private:
  TGMainFrame* fMainFrame; ///< pointer to our mainframe

  static const Int_t fgkFILESAVEAS; ///< File/Save As... menu
  static const Int_t fgkFILEOPEN; ///< File/Open... menu
  static const Int_t fgkFILEEXIT; ///< File/Exit menu

  ClassDef(AliMUONMchViewApplication,1) // mchview application
};

#endif
