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
//        Copied from NetParticle Classes
//        Origin: Authors: Jochen Thaeder <jochen@thaeder.de>
//                         Michael Weber <m.weber@cern.ch>
//=========================================================================//


#include "THnSparse.h"
#include "AliEbyEPidRatioBase.h"

class AliEbyEPidRatioQA : public AliEbyEPidRatioBase {

 public:

  AliEbyEPidRatioQA();
  virtual ~AliEbyEPidRatioQA();
  virtual void Process();
  THnSparseF* GetHnQAPid()  {return fHnQAa;}
  THnSparseF* GetHnQADca()  {return fHnQAb;}

 private:

  AliEbyEPidRatioQA(const AliEbyEPidRatioQA&); // not implemented
  AliEbyEPidRatioQA& operator=(const AliEbyEPidRatioQA&); // not implemented
  virtual void CreateHistograms();

  THnSparseF           *fHnQAa;        //! THnSparseF : tracks for QA
  THnSparseF           *fHnQAb;        //! THnSparseF : tracks for QA

  ClassDef(AliEbyEPidRatioQA, 1);
};

#endif
