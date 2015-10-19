#ifndef ALI_T0_PREPROCESSOROFFLINE_H
#define ALI_T0_PREPROCESSOROFFLINE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id: AliT0Preprocessor.h 39249 2010-02-28 19:09:52Z alla $ */


// T0 preprocessor. 
#include "TNamed.h"
class AliCDBStorage;

class AliT0PreprocessorOffline: public TNamed 
{
  public:
  AliT0PreprocessorOffline();  
  virtual ~AliT0PreprocessorOffline();
  void	CalibOffsetChannels(TString FileName, Int_t ustartRun, Int_t uendRun, AliCDBStorage* ocdbStorage);
  void	CalibT0sPosition(TString FileName, Int_t ustartRun, Int_t uendRun, AliCDBStorage* ocdbStorage);
  void  Process(TString FileName, Int_t ustartRun, Int_t uendRun, AliCDBStorage* ocdbStorage);
  void setDArun(Int_t runnumber) {fNewDArun = runnumber; };
  Int_t GetStatus() const;

  private:
  AliT0PreprocessorOffline(const AliT0PreprocessorOffline & proc); // copy constructor	
  AliT0PreprocessorOffline& operator=(const AliT0PreprocessorOffline&); //operator
  Int_t startRun;                         // start Run - used to make fast selection in THnSparse
  Int_t endRun;                           // end   Run - used to make fast selection in THnSparse
  Int_t startTime;                        // startTime - used to make fast selection in THnSparse
  Int_t endTime;                          // endTime   - used to make fast selection in THnSparse
  AliCDBStorage*  ocdbStorage;            // OCDB storage
  Int_t fNewDArun;                         // run number with new DA
  Int_t fStatusDelay;                     //status time delay calibration  
  Int_t fStatusAdjust;                   // status time adjust calibration

  ClassDef(AliT0PreprocessorOffline, 3)
};


#endif
