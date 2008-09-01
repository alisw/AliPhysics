#ifndef ALITRDcalib_H
#define ALITRDcalib_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDcalib.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"

class AliESDEvent;
class AliESDfriend;
class TList;
class AliTRDCalibraFillHisto;


class AliTRDcalib : public AliAnalysisTask{
public:

  AliTRDcalib(const Char_t *name = "TRD Track Info");
  ~AliTRDcalib(){};
  
  void   ConnectInputData(Option_t *);
  void   CreateOutputObjects();
  Int_t  GetDebugLevel() const {return fdebugLevel;} 
  Int_t  GetLow() { return flow; }
  Int_t  GetHigh() { return fhigh; }
  Bool_t GetFillZero() { return ffillZero; }
  void   Exec(Option_t *);
  void   SetDebugLevel(Int_t level)    { fdebugLevel = level; };
  void   SetLow(Int_t low)             { flow = low; };
  void   SetHigh(Int_t high)           { fhigh = high; };
  void   SetFillZero(Bool_t fillZero)  { ffillZero = fillZero; };
  void   SetSpecificStorage(const char *specificstorage) { fspecificstorage = specificstorage; };
  void   Terminate(Option_t *);

private:

  AliTRDcalib(const AliTRDcalib&);
  AliTRDcalib& operator=(const AliTRDcalib&);

private:

  AliESDEvent      *fESD;                         // ESD event
  AliESDfriend     *fESDfriend;                   // ESD friends

  TList            *fListHist;                    //! list of histograms
  AliTRDCalibraFillHisto *fTRDCalibraFillHisto;   //calibration analyse object

  Int_t             flow;                         // low number of clusters
  Int_t             fhigh;                        // high number of clusters
  Bool_t            ffillZero;                    // fill with zeros
  Int_t             fdebugLevel;                  // debug level
  const char *      fspecificstorage;             // number of timebins

   
  ClassDef(AliTRDcalib, 1)          // entry to TRD analysis
};
#endif
