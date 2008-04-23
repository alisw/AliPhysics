#ifndef AliITSQASPDDataMakerSim_H
#define AliITSQASPDDataMakerSim_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//
//  Checks the quality assurance. 
//  By comparing with reference data
//  contained in a DB
//
//
//  W. Ferrarese + P. Cerello Feb 2008
//


/* $Id$ */

#include "AliQA.h"
class AliITSQADataMakerSim;
class TObjArray;
class TClonesArray;

class AliITSQASPDDataMakerSim : public TObject {

public:
  AliITSQASPDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim); //ctor
  AliITSQASPDDataMakerSim(const AliITSQASPDDataMakerSim& qadm);
  AliITSQASPDDataMakerSim& operator = (const AliITSQASPDDataMakerSim& qac);

  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list);
  virtual ~AliITSQASPDDataMakerSim() {;}   // dtor
  virtual void InitDigits();
  virtual void InitSDigits();
  virtual void InitHits();
  virtual void MakeDigits(TClonesArray * /*digits*/){;}
  virtual void MakeSDigits(TClonesArray * /*sdigits*/){;}
  virtual void MakeHits(TClonesArray * /*hits*/){;}
  virtual void MakeDigits(TTree * digits);
  virtual void MakeSDigits(TTree * sdigits);
  virtual void MakeHits(TTree * hits);
  const Int_t Digits() { return fSPDhDigits; }
  const Int_t SDigits() { return fSPDhSDigits; }
  const Int_t Hits() { return fSPDhHits; }

private:

  AliITSQADataMakerSim *fAliITSQADataMakerSim;//pointer to the main ctor
  Int_t   fSPDhDigits;                        //number of booked SPD Digits histograms;
  Int_t   fSPDhSDigits;                       //number of booked SPD SDigits histograms;
  Int_t   fSPDhHits;                          //number of booked SPD Hits histograms;
  Int_t   fDigitsOffset;                      // number of histo booked when SPD start
  Int_t   fSDigitsOffset;                     // number of histo booked when SPD start
  Int_t   fHitsOffset;                        // number of histo booked when SPD start
  ClassDef(AliITSQASPDDataMakerSim,1)      // description 

};

#endif

