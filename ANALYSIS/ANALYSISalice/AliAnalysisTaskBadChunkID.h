/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Modified version of AliAnalysisTaskCheckCascade.h
// Used bits of code from AliAnalysisTaskCheckPerformanceStrange
//
// --- David Dobrigkeit Chinellato
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliAnalysisTaskBadChunkID_H
#define AliAnalysisTaskBadChunkID_H

class TList;
class TH1F;
class AliESDEvent;

#include "TString.h"
//#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskBadChunkID : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskBadChunkID();
  AliAnalysisTaskBadChunkID(const char *name);
  virtual ~AliAnalysisTaskBadChunkID();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
private:
  // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
  // your data member object is created on the worker nodes and streaming is not needed.
  // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
  
  TList  *fList;  //! List of Control plots
  TTree  *fTree;  //! Output Tree
  
  TH1F      *fHistNEvents;    //! Keeping track of N(events)
  
  Int_t fRunNumber;       //! Run Number
  TString fFileName;      //! Chunk Number
  Int_t fNGlobalTracks;   //! Number of kITSrefit tracks
  Int_t fNTracks;         //! Number of tracks
  
  //Objects Controlling Task Behaviour
  
  AliAnalysisTaskBadChunkID(const AliAnalysisTaskBadChunkID&);            // not implemented
  AliAnalysisTaskBadChunkID& operator=(const AliAnalysisTaskBadChunkID&); // not implemented
  
  ClassDef(AliAnalysisTaskBadChunkID, 11);
};

#endif
