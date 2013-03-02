#ifndef ALITPCQADATAMAKERSIM_H
#define ALITPCQADATAMAKERSIM_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: $ */

/*
  Based on AliPHOSQADataMaker
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  P. Christiansen, Lund, January 2008
*/


// --- ROOT system ---
#include <TH1.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include <AliLog.h>
#include <AliQADataMakerSim.h>
#include <AliRawReader.h>

#include <AliTPCdataQA.h>

class AliTPCQADataMakerSim: public AliQADataMakerSim {

public:
  
  enum HHitType_t   {kNhits=0, kElectrons, kRadius, kPrimPerCm, kElectronsPerCm} ; 
  enum HDigitType_t {kDigitsADC=0} ; 
  AliTPCQADataMakerSim() ;          // ctor
  AliTPCQADataMakerSim(const AliTPCQADataMakerSim& qadm) ;   
  AliTPCQADataMakerSim& operator = (const AliTPCQADataMakerSim& qadm) ;
  virtual ~AliTPCQADataMakerSim() { ; } // dtor
  
private:
  virtual void   StartOfDetectorCycle() {}; // empty 
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray** list) ;

  // Digits QA
  virtual void   InitDigits();
  virtual void   MakeDigits(TTree *digitTree);
  virtual void   MakeDigits() {AliWarning("Method not implemented\n");}

  // Hits QA
  virtual void   InitHits();
  virtual void   MakeHits(TTree *hitTree);
  virtual void   MakeHits() {AliWarning("Method not implemented\n");}

  // SDigits QA (empty)
  virtual void   InitSDigits() {}
  virtual void   MakeSDigits(TTree* ) {AliWarning("Method not implemented\n");}
  virtual void   MakeSDigits() {AliWarning("Method not implemented\n");}

  ClassDef(AliTPCQADataMakerSim,1)  // TPC Sim Quality Assurance Data Maker 
};

#endif // ALITPCQADATAMAKERSIM_H
