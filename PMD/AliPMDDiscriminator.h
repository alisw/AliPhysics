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

class AliPMDcluster;
class AliPMDclupid;

class AliPMDDiscriminator
{

 public:

  AliPMDDiscriminator();
  virtual ~AliPMDDiscriminator();

  void Discrimination(TObjArray *pmdcontin, TObjArray *pmdcontout);
  void EmpDiscrimination(TObjArray *pmdcontin, TObjArray *pmdcontout);
  void NNDiscrimination();

  void SetDiscrimination(Int_t idiscrim);

 protected:

  Int_t   fDiscrim;       // To switch on different discrimination method
  Int_t   fDebug;         // Debugging switch (0:NO, 1:YES)

  ClassDef(AliPMDDiscriminator,1) // To run PMD discrimination
};
#endif

