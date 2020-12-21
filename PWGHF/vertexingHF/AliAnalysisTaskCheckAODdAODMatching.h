#ifndef ALIANALYSISTASKCHECKAODDAODMATCHING_H
#define ALIANALYSISTASKCHECKAODDAODMATCHING_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. */

/////////////////////////////////////////////////////////////////////////////////////////
// \class AliAnalysisTaskCheckAODdAODMatching                                          //
// \brief analysis task for tagging the AOD files with AOD-dAOD mismatches             //
// \author: F. Grosa, fabrizio.grosa@cern.ch                                           //
// \author: F. Prino, francesco.prino@cern.ch                                          //
// \author: C. Terrevoli, cristina.terrevoli@cern.ch                                   //
/////////////////////////////////////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TTree.h>
#include <TList.h>
#include <TString.h>

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"

using namespace std;

class AliAnalysisTaskCheckAODdAODMatching : public AliAnalysisTaskSE {

public:

    enum MismatchStatus
    {
        kMismatchNevents = BIT(0), // mismatches tagged with number of events
        kMismatchTProcessID = BIT(1), // mismatches tagged with TProcessID
        kMismatchCand = BIT(2), // mismatches tagged with checks on tracks and candidates
    };

    AliAnalysisTaskCheckAODdAODMatching();
    AliAnalysisTaskCheckAODdAODMatching(const char *name);
    virtual ~AliAnalysisTaskCheckAODdAODMatching();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);

private:

    void MapAODtracks();

    TList* fOutput;                                                                 //!<! output list for histograms
    TH1F* fHistFiles;                                                               //!<! histo with event info
    TTree* fTreeMismatch;                                                           //!<! TTree with names of AOD files with file names and mismatch status

    AliAODEvent* fAOD;                                                              /// AOD event
    int fAODMap[1000000];                                                           /// AOD map
    TString fPrevInputFileName;                                                     /// AOD input file name of previous event
    unsigned char fStatus;                                                          /// bit map with mismatch status
    int fRunNumber;                                                                 /// run number
    int fNevents;                                                                   /// number of events in AOD file

    ClassDef(AliAnalysisTaskCheckAODdAODMatching, 1);
};

#endif
