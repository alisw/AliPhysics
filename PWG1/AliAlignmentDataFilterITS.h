#ifndef ALIALIGNMENTDATAFILTERITS_H
#define ALIALIGNMENTDATAFILTERITS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAlignmentDataFilterITS
// AliAnalysisTaskSE to extract from ESD tracks the AliTrackPointArrays
// with ITS points for selected tracks. This are the input data for alignment
// Author: A.Dainese, andrea.dainese@ln.infn.it
//*************************************************************************

class TTree;
class TNtuple;
class TList;
class TH1F;
class TH2F;

#include "AliAnalysisTask.h"

class AliAlignmentDataFilterITS : public AliAnalysisTask
{
 public:

  AliAlignmentDataFilterITS();
  AliAlignmentDataFilterITS(const char *name);
  virtual ~AliAlignmentDataFilterITS();


  // Implementation of interface methods
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
 private:

  void FilterCosmic(const AliESDEvent *esd);
  void FilterCollision(const AliESDEvent *esd);

  AliAlignmentDataFilterITS(const AliAlignmentDataFilterITS &source);
  AliAlignmentDataFilterITS& operator=(const AliAlignmentDataFilterITS& source); 

  AliESDEvent  *fESD;        // ESD object
  AliESDfriend *fESDfriend;  // ESD friend object
  TList   *fListOfHistos;    //! list of histos: output slot 1
  TTree   *fspTree;          //! output tree with space points: output slot 0
  TH1F    *fHistNpoints;     //! output histogram
  TH1F    *fHistPt;          //! output histogram
  TH2F    *fHistLayer0;      //! output histogram
  TH2F    *fHistLayer1;      //! output histogram
  TH2F    *fHistLayer2;      //! output histogram
  TH2F    *fHistLayer3;      //! output histogram
  TH2F    *fHistLayer4;      //! output histogram
  TH2F    *fHistLayer5;      //! output histogram
  TNtuple *fntExtra;         //! output QA ntuple  
  TNtuple *fntCosmicMatching;//! output QA ntuple  

  ClassDef(AliAlignmentDataFilterITS,0); // AliAnalysisTaskSE to extract ITS points for alignment
};

#endif

