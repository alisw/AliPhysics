#ifndef ALIITSQASSDCHECKER_H
#define ALIITSQASSDCHECKER_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  INFN Torino
//  P. Cerello - apr 2008
//


// --- ROOT system ---
class TFile ; 
class TH2F ;  
class TCanvas;

// --- AliRoot header files ---
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliQACheckerBase.h"


class AliITSLoader ; 

class AliITSQASSDChecker: public TObject {

public:
  AliITSQASSDChecker():fSubDetOffset(0),fStepBitSSD(NULL),fLowSSDValue(NULL),fHighSSDValue(NULL), fImage(NULL) {;}          // ctor
  AliITSQASSDChecker& operator = (const AliITSQASSDChecker& qac) ; //operator =
  virtual ~AliITSQASSDChecker() {if(fStepBitSSD) delete[] fStepBitSSD ;if(fLowSSDValue)delete[]fLowSSDValue;if(fHighSSDValue) delete[]fHighSSDValue; } // dtor
  virtual Double_t Check(AliQAv1::ALITASK_t /*index*/, const TObjArray * list, const AliDetectorRecoParam * recoParam);

  void SetStepBit(const Double_t *steprange);
  Double_t *GetStepBit() const {return fStepBitSSD;};

  void CheckRaws(TH1* histo);
  void CheckRecPoints(TH1 * /* histo */) const {return;}
  void SetTaskOffset(Int_t TaskOffset);
  void SetSSDLimits(const Float_t *lowvalue, const Float_t * highvalue);

  virtual Bool_t   MakeSSDImage( TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode) ;
  //virtual  Bool_t   MakeSSDImage( TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode) 
  //{ AliInfo(Form("Use default MakeImage method for the  %s for  task %s mode %s \n",list[0]->ClassName(), AliQAv1::GetTaskName(task).Data(), AliQAv1::GetModeName(mode))); return kFALSE;} 
  Bool_t MakeSSDRawsImage(TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode );//{AliInfo("The method for raw image has been called\n");}



private:
  
  AliITSQASSDChecker(const AliITSQASSDChecker& qac):TObject(),fSubDetOffset(qac.fSubDetOffset),fStepBitSSD(qac.fStepBitSSD),fLowSSDValue(qac.fLowSSDValue),fHighSSDValue(qac.fHighSSDValue), fImage(qac.fImage) {;} // cpy ctor   
  Int_t fSubDetOffset;            // checking operation starting point
  Double_t *fStepBitSSD;          // step size 
  Float_t *fLowSSDValue;          // low limit
  Float_t *fHighSSDValue;         // high limit

  TCanvas **    fImage          ; //[AliRecoParam::kNSpecies] 
  
  ClassDef(AliITSQASSDChecker,3)  // description 

};

#endif // AliITSQASSDChecker_H
