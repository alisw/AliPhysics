#ifndef AliMeanVertexCalibTask_cxx
#define AliMeanVertexCalibTask_cxx

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Class AliMeanVertexCalibTask
// AliAnalysisTask to extract from ESD the information on primary vertex
// reconstruction in order to compute the MeanVertex object
//
// Author:D.Caffarri, davide.caffarri@pd.infn.it  
//        A.Dainese,  andrea.dainese@pd.infn.it
//*************************************************************************

class TList;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliMeanVertexCalibTask : public AliAnalysisTaskSE
{
  
 public: 
  AliMeanVertexCalibTask(const char *name = "AliMeanVertexCalibTask");
  virtual ~AliMeanVertexCalibTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void           SetOnlyITSTPCTracks() {fOnlyITSTPCTracks=kFALSE;}
  void           SetOnlyITSSATracks() {fOnlyITSSATracks=kTRUE;}

  
    
  
 private:    
  AliESDEvent *fESD;            // ESD object
  TList       *fOutput;         //! list send on output slot 0

  Bool_t       fOnlyITSTPCTracks; // only ITS-TPC tracks to redo ITSTPC vertex
  Bool_t       fOnlyITSSATracks;  // only ITS-SA tracks to redo ITSTPC vertex
  
  
  AliMeanVertexCalibTask(const AliMeanVertexCalibTask&);
  AliMeanVertexCalibTask& operator=(const AliMeanVertexCalibTask&);
 
  AliESDVertex* ReconstructPrimaryVertex(Bool_t constr=kFALSE, Int_t mode=0) const;
  
  ClassDef(AliMeanVertexCalibTask, 1);
  
};

#endif


