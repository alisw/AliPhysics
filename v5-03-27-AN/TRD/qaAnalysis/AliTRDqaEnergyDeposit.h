#ifndef ALITRDQAENERGYDEPOSIT_H
#define ALITRDQAENERGYDEPOSIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id: AliTRDqaEnergyDeposit.h  $ */

//
// This class is a part of a package of high level QA monitoring for TRD.
// In this class the dEdX is analyzed as a function of the particle type,
// momentum and chamber.
//
// S. Radomski
// radomski@physi.uni-heidelberg.de
// March 2008
//

#include "AliAnalysisTask.h"  

class TTree; 
class AliESDEvent; 
class TH1D; 
class TH2D;
class AliExternalTrackParam;

class AliTRDqaEnergyDeposit : public AliAnalysisTask {

public:
  AliTRDqaEnergyDeposit();
  AliTRDqaEnergyDeposit(const char *name);
  AliTRDqaEnergyDeposit(const AliTRDqaEnergyDeposit & trd);
  AliTRDqaEnergyDeposit &operator=(const AliTRDqaEnergyDeposit & /*g*/) { return *this; };
  virtual ~AliTRDqaEnergyDeposit() {}
   
  virtual void Exec(Option_t * opt = "");
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "");

protected:
 
  TTree        *fChain;             //!pointer to the analyzed TTree or TChain
  AliESDEvent  *fESD;               //! Declaration of leave types

  TObjArray *fOutputContainer;      //! output data container

  // histogrms

  //TH2D *fSignalPt[2];       // pt-dedx distribution for pos and neg
  TH2D *fSignalPtSum[2];      // pt-dedx distribution for pos and neg
  TH2D *fSignalPtType[2*5];   // weight with the PID probability 
  TH1D *fProb[2*5];           // probabilities
  //TH2D *fSignalPtChamber[2*540];

  TH2D *fSignalPtPure[2*5];   // dedx for a clean sample 

  void FillElectrons() const;
  void FillPions() const {}
  void FillKaons() const {}
  void FillProtons() const {}

  ClassDef(AliTRDqaEnergyDeposit, 0); // a TRD analysis task 
};
#endif // ALITRDQAENERGYDEPOSIT_H
