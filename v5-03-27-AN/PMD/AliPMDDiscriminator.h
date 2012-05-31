#ifndef ALIPMDDISCRIMINATOR_H
#define ALIPMDDISCRIMINATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Date   : March 03 2004                             //
//  Does photon and Hadron discrimination              //
//                                                     //
//-----------------------------------------------------//

class TClonesArray;
class TFile;
class TObjArray;
class TTree;
class TNtuple;

class AliPMDrecpoint1;
class AliPMDclupid;

class AliPMDDiscriminator : public TObject
{

 public:

  AliPMDDiscriminator(){};
  virtual ~AliPMDDiscriminator(){};

  virtual void Discrimination(TObjArray *pmdcontin, TObjArray *pmdcontout) = 0;

  ClassDef(AliPMDDiscriminator,3) // Base class for PMD discrimination
};
#endif

