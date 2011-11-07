#ifndef ALIITSPREPROCESSORSDD_H
#define ALIITSPREPROCESSORSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////
//                                                //
//  Class for SDD data preprocessing              //
//  Origin: E. Crescio - crescio@to.infn.it       //
//          F. Prino   - prino@to.infn.t          //
//                                                //
////////////////////////////////////////////////////


#include "AliPreprocessor.h"
#include "AliITSDDLModuleMapSDD.h"
#include "AliITSDriftSpeedSDD.h"


class AliITSPreprocessorSDD : public AliPreprocessor { 
 

 public:
 
  AliITSPreprocessorSDD( AliShuttleInterface* shuttle);
  virtual ~AliITSPreprocessorSDD(){;}

  enum {kNumberOfSDD = 260};    // number of SDD modules
  enum {kNumberOfSDDLay3 = 84};    // number of SDD modules on layer 3
  enum {kNumberOfDDL = 24};     // number of DDLs in SDD
  enum {kModulesPerDDL = 12};   // number of modules in each DDL 
  enum {kNumberOfChannels = 512}; // number of channels per module

 protected:      


  
  virtual UInt_t Process(TMap* dcsAliasMap);
  UInt_t ProcessPulser(AliITSDDLModuleMapSDD* ddlmap);
  UInt_t ProcessInjector(AliITSDDLModuleMapSDD* ddlmap);
  Bool_t ProcessDCSDataPoints(TMap* dcsAliasMap);
  Double_t* RescaleDriftSpeedModule(TList* theList, Int_t iBadMod, Int_t iBadSide, AliITSDriftSpeedSDD* refSpeed) const;

  static const TString fgkNameHistoPedestals; //name of ped. histo
  static const TString fgkNameHistoNoise;  //name of noise histo
  ClassDef(AliITSPreprocessorSDD,4)  // Alice ITS-SDD preprocessor.

 };



#endif

    
