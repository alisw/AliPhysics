#ifndef ALIJETKINEREADERHEADER_H
#define ALIJETKINEREADERHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet Kinematics Reader Header
// Header for the Kine reader in the jet analysis
// Author: Andreas Morsch (andreas.morsch@cern.ch)

#include "AliJetReaderHeader.h"
 
class AliJetKineReaderHeader : public AliJetReaderHeader
{

 public:
  AliJetKineReaderHeader();
  virtual ~AliJetKineReaderHeader();
  
  // Setters
  void SetFastSimTPC(Bool_t flag = kTRUE) {fFastSimTPC = flag;} // if TPC fast simulation
  void SetFastSimEMCAL(Bool_t flag = kTRUE) {fFastSimEMCAL = flag;} // if EMCAL fast simulation
  void SetChargedOnly(Bool_t flag = kTRUE) {fChargedOnly = flag;} // for charged particles only, no smearing and no acceptance cuts

  // Getter
  Bool_t  FastSimTPC() const  {return fFastSimTPC;}
  Bool_t  FastSimEMCAL() const  {return fFastSimEMCAL;}
  Bool_t  ChargedOnly() const  {return fChargedOnly;}

	  
 protected:
  //parameters set by user
  Bool_t fFastSimTPC;   // TPC fast simulation flag
  Bool_t fFastSimEMCAL; // EMCAL fast simulation flag
  Bool_t fChargedOnly;  // Charged particle only flag

  ClassDef(AliJetKineReaderHeader,2);
};
 
#endif
