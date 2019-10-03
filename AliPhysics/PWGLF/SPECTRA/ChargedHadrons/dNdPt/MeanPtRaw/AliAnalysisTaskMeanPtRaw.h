#ifndef AliAnalysisTaskMeanPtRaw_cxx
#define AliAnalysisTaskMeanPtRaw_cxx

// Generation of raw <pT> vs mult 
// by Philipp Luettig, 03.02.2015
// last modified: 08.05.2015

class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class THnSparse;
class TList;
class TH1;

#include <vector>


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMeanPtRaw : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskMeanPtRaw(const char *name = "AliAnalysisTaskMeanPtRaw");
    virtual ~AliAnalysisTaskMeanPtRaw() {}

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void AddAliESDtrackCut(AliESDtrackCuts *esdTrackCuts, const char *name = "esdTrackCut");
    void SetMaxVertexZ(Float_t vZ) {fMaxVertexZ = vZ;}

  private:
    AliESDEvent                 *fESD;                  //ESD object
    TList                       *fOutputList;           // List where all the output files are stored

    // Histogram
    THnSparse					*fPtVsMultRaw; // pt vs multiplicity raw
    TH1I						*fTrackCutName; // identifier of track cuts

    std::vector<AliESDtrackCuts*>    fESDTrackCuts; //[fNTrackCuts] Esd track cuts
    std::vector<TString>			fTrackCutNames;
    Int_t fNTrackCuts; // number of tracks cuts sets

    //variables
    Float_t                     fMaxVertexZ;            // Maximum value for Vertex Z position
    Int_t						fEventCount;
    
    
    AliAnalysisTaskMeanPtRaw(const AliAnalysisTaskMeanPtRaw&); // not implemented
    AliAnalysisTaskMeanPtRaw& operator=(const AliAnalysisTaskMeanPtRaw&); // not implemented

    ClassDef(AliAnalysisTaskMeanPtRaw, 1); // example of analysis
};

#endif

