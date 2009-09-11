#ifndef ALIITSDDLMODULEMAPSDD_H
#define ALIITSDDLMODULEMAPSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id:  $ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to store SDD DDL mapping in the OCDB                    //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include<TObject.h>
#include<AliITSgeomTGeo.h>


class AliITSDDLModuleMapSDD : public TObject {

 public:

  AliITSDDLModuleMapSDD();
  AliITSDDLModuleMapSDD(Char_t *ddlmapfile);
  virtual ~AliITSDDLModuleMapSDD(){};

  void SetDefaultMap();
  void SetDec07part1Map();
  void SetDec07part2Map();
  void SetFeb08Map();
  void SetJun08Map();
  void SetJun09Map();
  void SetDDLMapElement(Int_t iDDL, Int_t iChan, Int_t iMod){fDDLModuleMap[iDDL][iChan]=iMod;}
  void SetDDLMap(AliITSDDLModuleMapSDD* ddlmap);
  void ReadDDLModuleMap(Char_t *ddlmapfile);

  Int_t GetModuleNumber(UInt_t iDDL, UInt_t iChan) const {return fDDLModuleMap[iDDL][iChan];}
  void FindInDDLMap(Int_t modIndex, Int_t &iDDL, Int_t &iCarlos) const;
  void FindInDDLMap(Int_t lay, Int_t lad, Int_t det, Int_t &iDDL, Int_t &iCarlos) const {
    FindInDDLMap(AliITSgeomTGeo::GetModuleIndex(lay,lad,det),iDDL,iCarlos);
  }
  static Int_t GetNDDLs(){return kDDLsNumber;}
  static Int_t GetNModPerDDL(){return kModulesPerDDL;}


  void PrintDDLMap() const ;

 protected:
  
  enum {kDDLsNumber = 24};      // number of DDLs in SDD
  enum {kModulesPerDDL = 12};   // number of modules in each DDL 

  Int_t fDDLModuleMap[kDDLsNumber][kModulesPerDDL]; //  mapping DDL/module -> module number

  ClassDef(AliITSDDLModuleMapSDD,1);
};
#endif
