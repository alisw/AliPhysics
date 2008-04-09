#ifndef AliITSQASDDDataMakerSim_H
#define AliITSQASDDDataMakerSim_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//
//  Checks the quality assurance. 
//  By comparing with reference data
//  contained in a DB
//
//
//  W. Ferrarese + P. Cerello Feb 2008

/* $Id$ */

#include "AliQA.h"
class AliITSQADataMakerSim;
class AliRunLoader;
class AliRun;
class TObjArray;
class TClonesArray;

class AliITSQASDDDataMakerSim : public TObject {

public:
  AliITSQASDDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim); //ctor
  AliITSQASDDDataMakerSim(const AliITSQASDDDataMakerSim& qadm);
  AliITSQASDDDataMakerSim& operator = (const AliITSQASDDDataMakerSim& qac);

  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list);
  virtual ~AliITSQASDDDataMakerSim() {;}   // dtor
  virtual void InitDigits();
  virtual void InitSDigits();
  virtual void InitHits();
  virtual void MakeDigits(TClonesArray * /*digits*/){;}
  virtual void MakeSDigits(TClonesArray * /*sdigits*/){;}
  virtual void MakeHits(TClonesArray * /*hits*/){;}
  virtual void MakeDigits(TTree * digits);
  virtual void MakeSDigits(TTree * sdigits);
  virtual void MakeHits(TTree * hits);
  const Int_t Digits() { return fSDDhDigits; }
  const Int_t SDigits() { return fSDDhSDigits; }
  const Int_t Hits() { return fSDDhHits; }

private:

  AliITSQADataMakerSim *fAliITSQADataMakerSim;//pointer to the main ctor
  Int_t   fSDDhDigits;                        //number of booked SDD Digits histograms;
  Int_t   fSDDhSDigits;                       //number of booked SDD SDigits histograms;
  Int_t   fSDDhHits;                          //number of booked SDD Hits histograms;
  Int_t   fDigitsOffset;                      // number of histo booked when SDD start
  Int_t   fSDigitsOffset;                     // number of histo booked when SDD start
  Int_t   fHitsOffset;                        // number of histo booked when SDD start
  ClassDef(AliITSQASDDDataMakerSim,1)      // description 

};

#endif


