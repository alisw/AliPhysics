#ifndef ALIPMDEMPDISCRIMINATOR_H
#define ALIPMDEMPDISCRIMINATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Date   : March 03 2004                             //
//  Does photon and Hadron discrimination              //
//                                                     //
//-----------------------------------------------------//

#include "AliPMDDiscriminator.h"
class TClonesArray;
class TFile;
class TObjArray;
class TTree;
class TNtuple;

class AliPMDrecdata;
class AliPMDclupid;

class AliPMDEmpDiscriminator : public AliPMDDiscriminator
{

 public:

  AliPMDEmpDiscriminator();
  virtual ~AliPMDEmpDiscriminator();

  void Discrimination(TObjArray *pmdcontin, TObjArray *pmdcontout);


  ClassDef(AliPMDEmpDiscriminator,2) // To run PMD discrimination
};
#endif

