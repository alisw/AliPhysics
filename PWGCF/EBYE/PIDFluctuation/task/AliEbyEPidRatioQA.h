#ifndef ALIEBYEPIDRATIOQA_H
#define ALIEBYEPIDRATIOQA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    //
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//=========================================================================//


#include "THnSparse.h"
#include "AliEbyEPidRatioBase.h"

class AliEbyEPidRatioQA : public AliEbyEPidRatioBase {

 public:

  AliEbyEPidRatioQA();
  virtual ~AliEbyEPidRatioQA();
  virtual void Process();
  THnSparseF* GetHnQA()  {return fHnQA;}

 private:

  AliEbyEPidRatioQA(const AliEbyEPidRatioQA&); // not implemented
  AliEbyEPidRatioQA& operator=(const AliEbyEPidRatioQA&); // not implemented
  virtual void CreateHistograms();

  THnSparseF           *fHnQA;                  //! THnSparseF : tracks for QA

  ClassDef(AliEbyEPidRatioQA, 1);
};

#endif
