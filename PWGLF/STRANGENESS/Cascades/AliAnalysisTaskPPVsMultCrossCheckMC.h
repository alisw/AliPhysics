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

#ifndef AliAnalysisTaskPPVsMultCrossCheckMC_H
#define AliAnalysisTaskPPVsMultCrossCheckMC_H

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

class AliAnalysisTaskPPVsMultCrossCheckMC : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskPPVsMultCrossCheckMC();
    AliAnalysisTaskPPVsMultCrossCheckMC(const char *name);
    virtual ~AliAnalysisTaskPPVsMultCrossCheckMC();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
    
    void SetPureMonteCarlo( Bool_t lSet = kTRUE ) { lPureMonteCarlo = lSet; } 
    void SetCheckVtxZMC   ( Bool_t lSet = kTRUE ) { fCheckVtxZMC = lSet;    }
    void SetAlternateMCSelection ( Bool_t lSet = kTRUE ) { fAlternateMCSelection = lSet; }
    void SetSkipPS ( Bool_t lSet = kTRUE ) { fSkipPS = lSet; } 
    void SetUseRecoVtxZ ( Bool_t lSet = kTRUE ) { fUseRecoVtxZ = lSet; }
    //Task Configuration: trigger selection 
    void SetSelectedTriggerClass(AliVEvent::EOfflineTriggerTypes trigType) { fTrigType = trigType;}
    void SetSelectedTriggerClass(TString trigName) { fkSelectTriggerByName = kTRUE; fTrigName = trigName;}

    void SetUseMultSelection ( Bool_t lUseMultSelection = kTRUE) {
        fkMultSelection = lUseMultSelection;
    }
   
 
    Double_t GetV0MAmplitude ( AliESDEvent *lInputEvent ) const;
    
//---------------------------------------------------------------------------------------

private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;  //! List of Cascade histograms

    AliPIDResponse *fPIDResponse;     // PID response object
    AliESDtrackCuts *fESDtrackCuts;   // ESD track cuts used for primary track definition
    AliPPVsMultUtils *fPPVsMultUtils; //
    AliAnalysisUtils *fUtils; //

    //Histograms (Desired objects in this cross-checking task) 
    TH1F *fHistEventCounter; //! histogram for event counting
    
    //Will attempt to read ESD or not... 
    Bool_t lPureMonteCarlo;
    Bool_t fCheckVtxZMC;
    Bool_t fAlternateMCSelection;
    Bool_t fSkipPS;
    Bool_t fUseRecoVtxZ;
    Bool_t fkMultSelection;
    Bool_t fkSelectTriggerByName;  // to select trigger by name (if it's not availble in AliVEvent)
  
    AliVEvent::EOfflineTriggerTypes fTrigType; // trigger type
    TString   fTrigName; // trigger name (if it's not available in AliVEvent) 
 
    //Basic Histograms for counting events as a function of V0M percentiles...
    TH1F *fHistV0M_DataSelection; //!
    TH1F *fHistV0M_MCSelection;   //!
    TH1F *fHistV0MAmplitude_DataSelection; //!
    TH1F *fHistV0MAmplitude_MCSelection;   //!
    TH1F *fHistV0MTrue_DataSelection; //!
    TH1F *fHistV0MTrue_MCSelection;   //!
    //// same for Ntrackl08   
    TH1F *fHistTracklets08Cent_DataSelection; //!
    TH1F *fHistTracklets08Cent_MCSelection;   //!
    TH1F *fHistTracklets08_DataSelection; //!
    TH1F *fHistTracklets08_MCSelection;   //!
    TH1F *fHistTracklets08True_DataSelection; //!
    TH1F *fHistTracklets08True_MCSelection;   //!
    //// same for Ntrackl0815   
    TH1F *fHistTracklets0815Cent_DataSelection; //!
    TH1F *fHistTracklets0815Cent_MCSelection;   //!
    TH1F *fHistTracklets0815_DataSelection; //!
    TH1F *fHistTracklets0815_MCSelection;   //!
    TH1F *fHistTracklets0815True_DataSelection; //!
    TH1F *fHistTracklets0815True_MCSelection;   //!
 
    TH2F *fHistV0MVsMidRapidityTrue_DataSelection; //!
    TH2F* fHistV0MAmplitudeVsMidRapidityTrue_DataSelection; //!
    TH2F *fHistV0MTrueVsMidRapidityTrue_DataSelection; //!
    ////
    TH2F *fHistTracklets08CentVsMidRapidityTrue_DataSelection; //!
    TH2F* fHistTracklets08VsMidRapidityTrue_DataSelection; //!
    TH2F *fHistTracklets08TrueVsMidRapidityTrue_DataSelection; //!
    ////
    TH2F *fHistTracklets0815CentVsMidRapidityTrue_DataSelection; //!
    TH2F* fHistTracklets0815VsMidRapidityTrue_DataSelection; //!
    TH2F *fHistTracklets0815TrueVsMidRapidityTrue_DataSelection; //!
       
    TH2F *fHistV0MVsMidRapidityTrue_MCSelection; //!
    TH2F* fHistV0MAmplitudeVsMidRapidityTrue_MCSelection; //!
    TH2F *fHistV0MTrueVsMidRapidityTrue_MCSelection; //!
    ////
    TH2F *fHistTracklets08CentVsMidRapidityTrue_MCSelection; //!
    TH2F* fHistTracklets08VsMidRapidityTrue_MCSelection; //!
    TH2F *fHistTracklets08TrueVsMidRapidityTrue_MCSelection; //!
    ////
    TH2F *fHistTracklets0815CentVsMidRapidityTrue_MCSelection; //!
    TH2F* fHistTracklets0815VsMidRapidityTrue_MCSelection; //!
    TH2F *fHistTracklets0815TrueVsMidRapidityTrue_MCSelection; //!
   
 
    //Desired Spectra: pi/K/p/K0/Lambda/Xi/Omega/Phi/K*
    //Data Selection Spectra
    //1-dimensional with transverse momentum (integrated in multiplicity)
    TH1F *fHistPt_Generated[9];     //! for keeping track of base spectra
    TH1F *fHistPt_DataSelection[9]; //!
    TH1F *fHistPt_MCSelection[9];   //!
    
    //2-dimensional with unchecked V0M percentile...
    TH2F *fHistPtVsV0M_Generated[9]; //! 9 spectra
    TH2F *fHistPtVsV0M_DataSelection[9]; //! 9 spectra
    TH2F *fHistPtVsV0M_MCSelection[9];   //! 9 spectra
    
    //2-dimensional with V0M Amplitudes
    TH2F *fHistPtVsV0MAmplitude_Generated[9]; //! 9 spectra
    TH2F *fHistPtVsV0MAmplitude_DataSelection[9]; //! 9 spectra
    TH2F *fHistPtVsV0MAmplitude_MCSelection[9];   //! 9 spectra
    
    //2-dimensional with true generated counts in V0M acceptance
    TH2F *fHistPtVsV0MTrue_Generated[9]; //! 9 spectra
    TH2F *fHistPtVsV0MTrue_DataSelection[9]; //! 9 spectra
    TH2F *fHistPtVsV0MTrue_MCSelection[9];   //! 9 spectra
    
    ////
    //2-dimensional with unchecked Tracklets08 percentile...
    TH2F *fHistPtVsTracklets08Cent_Generated[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets08Cent_DataSelection[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets08Cent_MCSelection[9];   //! 9 spectra

    //2-dimensional with Tracklets08 Amplitudes
    TH2F *fHistPtVsTracklets08_Generated[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets08_DataSelection[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets08_MCSelection[9];   //! 9 spectra

    //2-dimensional with true generated counts in Tracklets08 acceptance
    TH2F *fHistPtVsTracklets08True_Generated[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets08True_DataSelection[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets08True_MCSelection[9];   //! 9 spectra
    /// 
     //2-dimensional with unchecked Tracklets08 percentile...
    TH2F *fHistPtVsTracklets0815Cent_Generated[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets0815Cent_DataSelection[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets0815Cent_MCSelection[9];   //! 9 spectra

    //2-dimensional with Tracklets0815 Amplitudes
    TH2F *fHistPtVsTracklets0815_Generated[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets0815_DataSelection[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets0815_MCSelection[9];   //! 9 spectra

    //2-dimensional with true generated counts in Tracklets0815 acceptance
    TH2F *fHistPtVsTracklets0815True_Generated[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets0815True_DataSelection[9]; //! 9 spectra
    TH2F *fHistPtVsTracklets0815True_MCSelection[9];   //! 9 spectra

    
    AliAnalysisTaskPPVsMultCrossCheckMC(const AliAnalysisTaskPPVsMultCrossCheckMC&);            // not implemented
    AliAnalysisTaskPPVsMultCrossCheckMC& operator=(const AliAnalysisTaskPPVsMultCrossCheckMC&); // not implemented

    ClassDef(AliAnalysisTaskPPVsMultCrossCheckMC, 1);
};

#endif
