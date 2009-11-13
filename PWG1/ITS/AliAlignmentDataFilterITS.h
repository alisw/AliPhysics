#ifndef ALIALIGNMENTDATAFILTERITS_H
#define ALIALIGNMENTDATAFILTERITS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAlignmentDataFilterITS
// AliAnalysisTask to extract from ESD tracks the AliTrackPointArrays
// with ITS points for selected tracks. This are the input data for alignment
// Author: A.Dainese, andrea.dainese@pd.infn.it
//*************************************************************************

class TTree;
class TNtuple;
class TList;
class TH1F;
class TH2F;
class TObjString;

class AliTrackPointArray;

#include <TString.h>
#include "AliITSReconstructor.h"
#include "AliITSRecoParam.h"
#include "AliAnalysisTask.h"

class AliAlignmentDataFilterITS : public AliAnalysisTask
{
 public:

  AliAlignmentDataFilterITS(const char *name="filterITS");
  virtual ~AliAlignmentDataFilterITS();


  // Implementation of interface methods
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
  void SetOnlySPDFO(Bool_t set=kTRUE) {fOnlySPDFO=set;}
  void SetGeometryFileName(TString name="geometry.root") {fGeometryFileName=name;}
  void SetITSRecoParam(AliITSRecoParam *rp) {fITSRecoParam=rp;}
  static Int_t WriteTrackPointsInIdealGeom(Char_t *fin="AliTrackPoints.root", 
					   Char_t *fout="AliTrackPoints_IdGeom.root",
					   Char_t *fmis="Run0_999999999_v3_s0.root",
					   Char_t *fgeo="geometry.root",
					   Bool_t prn=0);

 private:

  void FilterCosmic(const AliESDEvent *esd);
  void FilterCollision(const AliESDEvent *esd);
  const AliITSRecoParam *GetRecoParam() const;

  AliAlignmentDataFilterITS(const AliAlignmentDataFilterITS &);
  AliAlignmentDataFilterITS& operator=(const AliAlignmentDataFilterITS&);


  Bool_t fOnlySPDFO;         // only SPDtriggered events
  TString fGeometryFileName; // where to find the geometry.root
  AliITSRecoParam *fITSRecoParam;  // keeps the settings for the filter
  AliESDEvent  *fESD;        // ESD object
  AliESDfriend *fESDfriend;  // ESD friend object
  TList   *fListOfHistos;    //! list of histos: output slot 1
  TTree   *fspTree;          //! output tree with space points: output slot 0
  TH1F    *fHistNevents;     //! output histogram
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

  ClassDef(AliAlignmentDataFilterITS,2); // AliAnalysisTask to extract ITS points for alignment
};

#endif

