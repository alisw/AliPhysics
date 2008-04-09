#ifndef AliITSQASSDDataMakerSim_H
#define AliITSQASSDDataMakerSim_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  contained in a DB
//
//
//  W. Ferrarese + P. Cerello Feb 2008

#include "AliQA.h"
class AliITSQADataMakerSim;
class TObjArray;
class TClonesArray;

class AliITSQASSDDataMakerSim : public TObject {

public:
  AliITSQASSDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim); //ctor
  AliITSQASSDDataMakerSim(const AliITSQASSDDataMakerSim& qadm);
  AliITSQASSDDataMakerSim& operator = (const AliITSQASSDDataMakerSim& qac);

  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list);
  virtual ~AliITSQASSDDataMakerSim() {;}   // dtor
  virtual void InitDigits();
  virtual void InitSDigits();
  virtual void InitHits();
  virtual void MakeDigits(TClonesArray * /*digits*/){;}
  virtual void MakeSDigits(TClonesArray * /*sdigits*/){;}
  virtual void MakeHits (TClonesArray * /*hits*/){;}
  virtual void MakeDigits(TTree * digits);
  virtual void MakeSDigits(TTree * sdigits);
  virtual void MakeHits(TTree * hits);
  const Int_t Digits() { return fSSDhDigits; }
  const Int_t SDigits() { return fSSDhSDigits; }
  const Int_t Hits() { return fSSDhHits; }

private:

  AliITSQADataMakerSim *fAliITSQADataMakerSim;//pointer to the main ctor
  Int_t   fSSDhDigits;                        //number of booked SSD Digits histograms;
  Int_t   fSSDhSDigits;                       //number of booked SSD SDigits histograms;
  Int_t   fSSDhHits;                          //number of booked SSD Hits histograms;
  Int_t   fDigitsOffset;                      // number of histo booked when SSD start
  Int_t   fSDigitsOffset;                     // number of histo booked when SSD start
  Int_t   fHitsOffset;                        // number of histo booked when SSD start
  ClassDef(AliITSQASSDDataMakerSim,1)      // description 

};

#endif


