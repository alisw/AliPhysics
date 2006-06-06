#ifndef ALI_DEFAULT_PREPROCESSOR_H
#define ALI_DEFAULT_PREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Default preprocessor that writes the retrieved DCS values to CDB
//

#include <AliPreprocessor.h>

class AliDefaultPreprocessor: public AliPreprocessor {

public:
	AliDefaultPreprocessor(const char* detector, AliShuttleInterface* shuttle);
	virtual ~AliDefaultPreprocessor();

	virtual UInt_t Process(TMap* dcsAliasMap);

private:
  ClassDef(AliDefaultPreprocessor, 0);
};

#endif
