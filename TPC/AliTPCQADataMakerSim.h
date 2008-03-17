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
  AliTPCQADataMakerSim() ;          // ctor
  AliTPCQADataMakerSim(const AliTPCQADataMakerSim& qadm) ;   
  AliTPCQADataMakerSim& operator = (const AliTPCQADataMakerSim& qadm) ;
  virtual ~AliTPCQADataMakerSim() { ; } // dtor
  
private:
  virtual void   StartOfDetectorCycle() {}; // empty 
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX, TObjArray* list) ;

  // Digits QA
  virtual void   InitDigits();
  virtual void   MakeDigits(TTree *digitTree);

  // Hits QA
  virtual void   InitHits();
  virtual void   MakeHits(TTree *hitTree);

  // SDigits QA (empty)
  virtual void   InitSDigits() {;}
  virtual void   MakeSDigits(TTree* ) {  AliInfo("AliTPCQADataMakerSim::MakeSDigits() method not supported"); }

  TH1F* fHistDigitsADC;    //! Digit ADC distribution

  TH1F* fHistHitsNhits;    //! Number of primary interactions
  TH1F* fHistHitsElectrons;//! Electrons per interaction
  TH1F* fHistHitsRadius;   //! Radius of interactions
  TH1F* fHistHitsPrimPerCm;//! Primaries per cm
  TH1F* fHistHitsElectronsPerCm;//! Electrons per cm

  ClassDef(AliTPCQADataMakerSim,1)  // TPC Sim Quality Assurance Data Maker 
};

#endif // ALITPCQADATAMAKERSIM_H
