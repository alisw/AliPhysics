//-------------------------------------------------------------------------
//     Task for the Analysis Framework 
// Updates the an already created AOD with the PWG2 information taken from 
// the ESD.
//  - Puts the per-track information into the AliPWG2AODTrack container, 
//    together with the link to the original AliAODTrack
//
//     Author: Adam Kisiel, OSU, Adam.Kisiel@cern.ch
//-------------------------------------------------------------------------
#ifndef ALIANALYSISTASKPWG2AODUPDATE_H
#define ALIANALYSISTASKPWG2AODUPDATE_H
 
#include <AliAnalysisTaskSE.h>

class AliESDEvent;
class AliAODEvent;
class TClonesArray;

class AliAnalysisTaskPWG2AODUpdate : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskPWG2AODUpdate();
    AliAnalysisTaskPWG2AODUpdate(const char* name);
    AliAnalysisTaskPWG2AODUpdate(const AliAnalysisTaskPWG2AODUpdate &task); 
    virtual ~AliAnalysisTaskPWG2AODUpdate() { ; }

    AliAnalysisTaskPWG2AODUpdate& operator=(const AliAnalysisTaskPWG2AODUpdate &task);
    // Implementation of interface methods
    virtual void LocalInit() {Init();}

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t* /*option*/);
    
 private:
    AliESDEvent*       fESD;           //! ESD
    AliAODEvent*       fAOD;           //! AOD event 
    TClonesArray*      fPWG2AODTracks; //! container for PWG2 specific information

    ClassDef(AliAnalysisTaskPWG2AODUpdate, 1); // Analysis task for standard ESD filtering
};
 
#endif
