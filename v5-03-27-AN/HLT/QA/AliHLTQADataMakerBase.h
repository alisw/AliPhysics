//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTQADATAMAKERBASE_H
#define ALIHLTQADATAMAKERBASE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTQADataMakerBase.h
    @author Matthias Richter
    @date   2010-03-10
    @brief  Base class for HLT detector QA data makers
*/

#include "AliQADataMakerRec.h"

class AliHLTQADataMakerBase: public AliQADataMakerRec {
 public:
  AliHLTQADataMakerBase();
  virtual ~AliHLTQADataMakerBase();

 protected:
  TObjArray** GetESDsQAList() {return fESDsQAList;}
  TObjArray** GetRawsQAList() {return fRawsQAList;}
  TObjArray** GetRecPointsQAList() {return fRecPointsQAList;}
  TObjArray** GetDigitsQAList() {return fDigitsQAList;}

  /** specific Exec handler which which can handle both the Esd and
   * HLTEsd in order to call MakeESDs with two parameters
   */
  virtual void Exec(AliQAv1::TASKINDEX_t task, TObject * data);

  /** dummy function, required by the QA framework, however
   * HLT QA is done in the other MakeESDs function
   */
  virtual void MakeESDs(AliESDEvent * esd);

  /** specific function with the two ESDs as argument
   */
  virtual void MakeESDs(AliESDEvent * esd, AliESDEvent* hltesd);

  friend class AliHLTQADataMakerRec;  
 private:
  /** copy constructor prohibited */
  AliHLTQADataMakerBase(const AliHLTQADataMakerBase&);   
  /** assignment operator prohibited */
  AliHLTQADataMakerBase& operator = (const AliHLTQADataMakerBase&);

  ClassDef(AliHLTQADataMakerBase,0)  // Base class for HLT QA Data Makers
};

#endif // ALIHLTQADATAMAKERBASE_H
