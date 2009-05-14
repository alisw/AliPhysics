//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTQADATAMAKERREC_H
#define ALIHLTQADATAMAKERREC_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTQADataMakerRec.h
    @author Matthias Richter
    @date   2009-05-14
    @brief  Container for the HLT offline QA
*/

#include "AliQADataMakerRec.h"

class AliHLTQADataMakerRec: public AliQADataMakerRec {

public:
  
  AliHLTQADataMakerRec();
  virtual ~AliHLTQADataMakerRec();
  
private:
  /** copy constructor prohibited */
  AliHLTQADataMakerRec(const AliHLTQADataMakerRec&);   
  /** assignment operator prohibited */
  AliHLTQADataMakerRec& operator = (const AliHLTQADataMakerRec&);

  virtual void   StartOfDetectorCycle();
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray** list);

  ClassDef(AliHLTQADataMakerRec,0)  // HLT Quality Assurance Data Maker for reconstruction
};

#endif // ALIHLTQADATAMAKERREC_H
