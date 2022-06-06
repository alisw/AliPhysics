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
// Modified version of AliAnalysisTaskMCPredictions.h from David Dobrigkeit Chinellato
//
// --- Francesca Ercolessi
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliAnalysisTaskMCPredictionsStrgVsMultVsZDC_H
#define AliAnalysisTaskMCPredictionsStrgVsMultVsZDC_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliPPVsMultUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMCPredictionsStrgVsMultVsZDC : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskMCPredictionsStrgVsMultVsZDC();
  AliAnalysisTaskMCPredictionsStrgVsMultVsZDC(const char *name, Float_t lCenterOfMassEnergy = 13000, Bool_t kDoPythia = kTRUE, Bool_t kDoEPOS = kFALSE);
  virtual ~AliAnalysisTaskMCPredictionsStrgVsMultVsZDC();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  Bool_t IsEPOSLHC() const;

  void SetSelectINELgtZERO ( Bool_t lOpt ) { fkSelectINELgtZERO = lOpt; }
 
  //---------------------------------------------------------------------------------------  
private:
  // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
  // your data member object is created on the worker nodes and streaming is not needed.
  // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
  TList  *fListHist;  //! List of Cascade histograms
  
  //Histograms 
  TH1D *fHistEventCounter; //! histogram for event counting
  TH1D *fHistV0MMult; //! 
  TH1D *fHistV0AMult; //! 
  TH1D *fHistV0CMult; //! 
  TH1D *fHistMult05; //! 
  TH1D *fHistMult08; //! 
  TH1D *fHistMult08to15; //! 
  TH1D *fHistSPDClusters; //! 
  TH1D *fHistNMPI; //! 
  TH1D *fHistQ2; //! 
  TH1D *fHistb; //! 
  TH1D *fHistLeadingE; //!
  TH1D *fHistEffEnergy; //! 
  TH2D *f2DHistINELgt0SPDV0M; //!
  TH2D *f2DHistLeadingESPDV0M; //!
  TH2D *f2DHistEffEnergySPDV0M; //!
  TH2D *f2DHistNchSPDV0M; //!
  TH2D *f2DHistNMPISPDV0M; //!
  TH2D *f2DHistQ2SPDV0M; //!
  TH2D *f2DHistbSPDV0M; //!
  TH1D *fHistPt[21]; //! 
  TH2D *f2DHistPartSPDV0M[21]; //!
  TH2D *f2DHistAvPtSPDV0M[21]; //!
  
  //Bool
  Bool_t fkSelectINELgtZERO;  
  Bool_t fkDoPythia; 
  Bool_t fkDoEPOS;  

  //Variables
  Float_t fCenterOfMassEnergy; 
  Int_t fkNSpecies; //! 
  
  AliAnalysisTaskMCPredictionsStrgVsMultVsZDC(const AliAnalysisTaskMCPredictionsStrgVsMultVsZDC&);            // not implemented
  AliAnalysisTaskMCPredictionsStrgVsMultVsZDC& operator=(const AliAnalysisTaskMCPredictionsStrgVsMultVsZDC&); // not implemented
  
  ClassDef(AliAnalysisTaskMCPredictionsStrgVsMultVsZDC, 1);
};

#endif

