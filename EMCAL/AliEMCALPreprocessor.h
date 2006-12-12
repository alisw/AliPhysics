#ifndef ALIEMCALPREPROCESSOR_H
#define ALIEMCALPREPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* History of cvs commits:
 *
 * $Log$
 * Revision 1.1  2006/12/07 16:32:16  gustavo
 * First shuttle code, online calibration histograms producer, EMCAL preprocessor
 *
 *
 */
///////////////////////////////////////////////////////////////////////////////
// Class AliEMCALPreprocessor
///////////////////////////////////////////////////////////////////////////////


#include "AliPreprocessor.h"

class AliEMCALPreprocessor : public AliPreprocessor {
public:

  AliEMCALPreprocessor();
  AliEMCALPreprocessor(AliShuttleInterface* shuttle);

protected:

  virtual UInt_t Process(TMap* valueSet);

  ClassDef(AliEMCALPreprocessor,0);

};

#endif
