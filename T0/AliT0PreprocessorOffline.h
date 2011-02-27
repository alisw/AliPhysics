#ifndef ALI_T0_PREPROCESSOROFFLINE_H
#define ALI_T0_PREPRECESSOROFFLINE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id: AliT0Preprocessor.h 39249 2010-02-28 19:09:52Z alla $ */


// T0 preprocessor. 
#include "TNamed.h"

class AliT0PreprocessorOffline: public TNamed 
{
  public:
  AliT0PreprocessorOffline();  
  virtual ~AliT0PreprocessorOffline();
  void	CalibOffsetChannels(TString FileName, Int_t ustartRun, Int_t uendRun, TString ocdbStorage);
 
  private:
	AliT0PreprocessorOffline(const AliT0PreprocessorOffline & proc); // copy constructor	
	AliT0PreprocessorOffline& operator=(const AliT0PreprocessorOffline&); //operator
 
	ClassDef(AliT0PreprocessorOffline, 1)
};


#endif
