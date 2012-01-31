//-------------------------------------------------------------------------
//     Task for the Analysis Framework 
// Creates a PWG2 specific AOD, from the ESD, adding information required
// for PWG2 analysis. 
//  - Puts the per-track information into the AliPWG2AODTrack container, 
//    together with the link to the original AliAODTrack
//
//     Author: Adam Kisiel, OSU, Adam.Kisiel@cern.ch
//-------------------------------------------------------------------------
#ifndef ALIANALYSISTASKPWG2ESDFILTER_H
#define ALIANALYSISTASKPWG2ESDFILTER_H
 
#include <AliAnalysisTask.h>

class AliESDEvent;
class TChain;
class AliAODEvent;
class AliAnalysisFilter;
class TTree;
class TClonesArray;
class TList;

class AliAnalysisTaskPWG2ESDfilter : public AliAnalysisTask
{
 public:
    AliAnalysisTaskPWG2ESDfilter();
    AliAnalysisTaskPWG2ESDfilter(const char* name);
    AliAnalysisTaskPWG2ESDfilter(const AliAnalysisTaskPWG2ESDfilter &task); 
    virtual ~AliAnalysisTaskPWG2ESDfilter() {;}

    AliAnalysisTaskPWG2ESDfilter& operator=(const AliAnalysisTaskPWG2ESDfilter &task);
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option = "");
    virtual void CreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void Exec(Option_t *option);
    virtual void Terminate(Option_t *option);
    // Setters
    virtual void SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
    virtual void SetKinkFilter (AliAnalysisFilter*  KinkF) {fKinkFilter  =  KinkF;}
    virtual void SetV0Filter   (AliAnalysisFilter*    V0F) {fV0Filter    =    V0F;}
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    
 private:
    Int_t              fDebug;         //  Debug flag
    TTree*             fTree;          //! chained files
    AliESDEvent*       fESD;           //! ESD
    AliAODEvent*       fAOD;           //! AOD event 
    TTree*             fTreeA;         //! AOD tree
    AliAnalysisFilter* fTrackFilter;   //  Track Filter
    AliAnalysisFilter* fKinkFilter;    //  Kink  Filter
    AliAnalysisFilter* fV0Filter;      //  V0    Filter    
    TClonesArray*      fPWG2AODTracks; //! container for PWG2 specific information

    ClassDef(AliAnalysisTaskPWG2ESDfilter, 1); // Analysis task for standard ESD filtering
};
 
#endif
