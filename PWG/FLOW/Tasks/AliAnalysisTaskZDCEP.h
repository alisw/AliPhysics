/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/***********************************
 * ZDC Event Plane                 *
 *                                 *
 * author: Jacopo Margutti         *
 * email:  jacopo.margutti@cern.ch *
 ***********************************/

#ifndef  AliAnalysisTaskZDCEP_H
#define  AliAnalysisTaskZDCEP_H

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"

class AliAnalysisUtils;
class AliAnalysisTaskSE;
class AliMultSelection;
class TProfile;
class TProfile2D;
class TProfile3D;
class TH3D;
class AliFlowEvent;

class  AliAnalysisTaskZDCEP : public  AliAnalysisTaskSE
{
public:
  //  two  class  constructors
  AliAnalysisTaskZDCEP();
  AliAnalysisTaskZDCEP(const  char *name);
  //  class  destructor
  virtual ~AliAnalysisTaskZDCEP();
  //  called  once  at  beginning  or  runtime
  virtual void UserCreateOutputObjects();
  //  called  for  each  event
  virtual void UserExec(Option_t* option);
  // get centrality bin
  virtual Int_t GetCenBin(Double_t Centrality);
  //  called  at  end  of  analysis
  virtual void Terminate(Option_t* option);

  enum DataSet {
    k2015o_pass1_pass1pidfix,
    k2015o_muon_calo_pass1,
  };

  void SetZDCCalibList(TList* const wlist) {this->fZDCCalibList = wlist;}
  TList* GetZDCCalibList() const {return this->fZDCCalibList;}
  void SetTowerEqList(TList* const wlist) {this->fTowerEqList = wlist;}
  TList* GetTowerEqList() const {return this->fTowerEqList;}
  void GetZDCQVectors(Double_t QAX, Double_t QAY, Double_t QCX, Double_t QCY);
  void SetDataSet(DataSet set) {this->fDataSet = set;};
  DataSet GetDataSet() const {return this->fDataSet;}

private:

  TList* fOutputList;             //! list containing ZDC q-vectors
  TList* fHistList;               //! list for calibration histograms
  TList* fQAList;                 //! output list for QA histograms (slot 2)
  Double_t fZDCGainAlpha;         //
  TList *fZDCCalibList;           // list for ZDC Q-vector re-centering
  TList *fTowerEqList;            // list for ZDC gain equalization
  TProfile* fZDCQHist[4];         //!
  TProfile3D* fZDCVtxHist[4];     //!
  TProfile2D* fZDCEcomTotHist[4]; //!
  TH3D *fZDCVtxFitHist[4];        //!
  TH1D *fZDCVtxFitCenProjHist[4][3]; //!
  TProfile3D *fZDCVtxCenHistMagPol[10][8]; //!
  TProfile3D* fZDCVtxCenHist[10][4]; //!
  TH1D* fCRCZDCQVecDummyEZDCBins[10]; //!
  TH2D* fZDCQvec2Ddis[10][2];      //!
  TProfile3D* fZDCCenVtxZ;         //!

  // QA histograms
  TList *fQAListMagPol;            //! QA list per magnet polarity
  TProfile*   fQVecCen[4][2];      //!
  TProfile3D* fQVecVtx[4][2];      //!
  TProfile*   fQVecCorCen[4][2];   //!
  // TProfile* fQVecDeltaC[4][2];     //!
  // TProfile* fQVecCorDeltaC[4][2];  //!
  TH2D* fQvecC2Ddis[2]; //!
  TH2D* fQvecA2Ddis[2]; //!
  TH1D* fEventCounter;             //!
  TH1D* fCentralityHisto;        //!

  TH3D *fZDCQVecVtxCenEZDC3D[10][10][4]; //!
  TH1D *fTowerGainEq[2][5];              //!

  AliFlowVector* fZDCFlowVect[2]; //! ZDC q-vectors
  Int_t fCachedRunNum;            //
  DataSet fDataSet;               //
  Int_t fnRun;                    //
  const static Int_t fnRunMax = 200; //
  TList *fQVecListRun[fnRunMax];        //! run-by-run list
  TProfile2D* fQVecRbRCen[fnRunMax];    //!
  TProfile3D* fQVecRbRVtxZ[fnRunMax];   //!
  const static Int_t fCRCnTow = 5;
  TProfile *fZNCTower[fnRunMax][fCRCnTow]; //! ZNC tower spectra
  TProfile *fZNATower[fnRunMax][fCRCnTow]; //! ZNA tower spectra
  TH2D *fZNCTowerSpec[fCRCnTow];  //! ZNC tower spectra
  TH2D *fZNATowerSpec[fCRCnTow];  //! ZNA tower spectra
  TArrayI fRunList;               // run list
  TArrayD fAvVtxPosX;             // average vx position vs run number
  TArrayD fAvVtxPosY;             // average vy position vs run number
  TArrayD fAvVtxPosZ;             // average vz position vs run number
  Bool_t fbIsMagnetPolarityNegative;     //
  AliFlowEvent* fFlowEvent;       // flowevent

  AliAnalysisUtils* fAnalysisUtils; //!
  AliMultSelection* fMultSelection; //!


  AliAnalysisTaskZDCEP(const AliAnalysisTaskZDCEP&);
  AliAnalysisTaskZDCEP& operator=(const AliAnalysisTaskZDCEP&);

  ClassDef(AliAnalysisTaskZDCEP,8);
};

#endif
