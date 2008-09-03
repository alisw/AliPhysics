#ifndef ALITRDTRACKINGEFFICIENCYCOMBINED_H
#define ALITRDTRACKINGEFFICIENCYCOMBINED_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackingEfficiencyCombined.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"

class TObjArray;
class TTreeSRedirector;

class AliTRDtrackingEfficiencyCombined : public AliAnalysisTask{
public:
  AliTRDtrackingEfficiencyCombined(const char *name = "combined tracking efficiency");
  virtual ~AliTRDtrackingEfficiencyCombined();
  
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *);
  virtual void Terminate(Option_t *);
  
  void SetDebugLevel(Int_t debugLevel) { fDebugLevel = debugLevel;}
  int GetDebugLevel() const { return fDebugLevel; }
  
private:
  AliTRDtrackingEfficiencyCombined(const AliTRDtrackingEfficiencyCombined &);
  AliTRDtrackingEfficiencyCombined& operator=(const AliTRDtrackingEfficiencyCombined &);
  
  TObjArray *fObjectContainer;	  	       //! Container for output histograms
  TObjArray *fTrackInfos;			       //! Input Container
  Int_t fDebugLevel;                             // Debug level
  TTreeSRedirector *fDebugStream;	 	       // Debug streamer
  
  ClassDef(AliTRDtrackingEfficiencyCombined, 1); // Combined tracking efficiency
};
		
#endif
