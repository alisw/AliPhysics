#ifndef AliCaloTriggerMimicHelper_H
#define AliCaloTriggerMimicHelper_H

#include "AliAnalysisTaskSE.h"
#include "TChain.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSTriggerUtils.h"
#include <iostream>
#include <vector>
#include <map>

class AliPHOSTriggerUtils;
using namespace std;

class AliCaloTriggerMimicHelper : public AliAnalysisTaskSE {

  public:
    AliCaloTriggerMimicHelper(const char *name="AliCaloTriggerMimicHelper", Int_t clusterType = 0, Int_t isMC=0);
    enum phosTriggerType{kPHOSAny,kPHOSL0,kPHOSL1low,kPHOSL1med,kPHOSL1high} ;
    //Uncopyable & operator=(const Uncopyable&);

    virtual ~AliCaloTriggerMimicHelper();                            //virtual destructor
    void UserCreateOutputObjects();

    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    void SetLightOutput( Bool_t flag )                      { fDoLightOutput = flag                         ; }
    void SetPHOSTrigger(phosTriggerType t=kPHOSL0)          { fPHOSTrigger=t                                ; }
    phosTriggerType GetPHOSTrigger()                        { return fPHOSTrigger                           ; }
    void SetRunNumber(int run)                              { fRunNumber=run; fForceRun=kTRUE               ; }             //Use given run number, don't read from data
    void SetEventChosenByTrigger( Bool_t flag )             { fEventChosenByTrigger = flag                  ; }
    void SetEventChosenByTriggerTrigUtils( Bool_t flag )    { fEventChosenByTriggerTrigUtils = flag         ; }
    Bool_t GetEventChosenByTrigger()                        { return fEventChosenByTrigger                  ; }
    Bool_t GetEventChosenByTriggerTrigUtils()               { return fEventChosenByTriggerTrigUtils         ; }
    AliPHOSGeometry*  GetGeomPHOS()                         { return fGeomPHOS                              ; }
    AliPHOSTriggerUtils* GetAliPHOSTriggerUtils()           { return fPHOSTrigUtils                         ; }
    AliVEvent* GetCurrentEvent()                            { return fInputEvent                            ; }
    void SetDebugOutput( Int_t flag )                       { fDoDebugOutput = flag                         ; }
    void SetTriggerHelperRunMode( Int_t flag )              { fTriggerHelperRunMode = flag                  ; }
    Int_t GetTriggerHelperRunMode()                         { return fTriggerHelperRunMode                  ; }
    TList* GetTriggerMimicHelperHistograms()                { return fOutputList                            ; }
    Int_t IsClusterIDTriggered(Int_t ClusterID)             { return fMapClusterToTriggered[ClusterID]      ; }
    Int_t IsClusterIDBadMapTrigger(Int_t ClusterID)         { return fMapClusterToTriggerMap[ClusterID]     ; }

  private:
    AliCaloTriggerMimicHelper (const AliCaloTriggerMimicHelper&);             // not implemented
    AliCaloTriggerMimicHelper & operator=(const AliCaloTriggerMimicHelper&);  // not implemented

    // private methods
    void SetClusterType(Int_t iClusterType)     { fClusterType = iClusterType                 ; }
    void SetTriggerDataOrMC(AliVCluster * clu, Bool_t isMCPhoton);


    // basic variables/objects
    TString                 fNameOfClassObject;                         // name of this class object
    Int_t                   fClusterType;                               // EMCal(1), PHOS(2) or not running (0)
    Int_t                   fRunNumber;                                 // current run number
    Int_t                   nModules;
    Int_t                   fNMaxPHOSModules;
    Int_t                   nMaxCellsPHOS;


    phosTriggerType         fPHOSTrigger;                               // Kind of PHOS trigger: L0,L1
    AliPHOSTriggerUtils *   fPHOSTrigUtils ;                            //! utils to analyze PHOS trigger
    AliPHOSGeometry*        fGeomPHOS;                                  //! pointer to PHOS geometry

    Bool_t                  fDoLightOutput;                             // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    Bool_t                  fForceRun ;                                 // use fixed run number, dont read from data
    Bool_t                  fIsMC ;                                     // Is this is MC
    Int_t                   fTriggerHelperRunMode;                      //0 is standard
    Bool_t                  fEventChosenByTrigger;                      //!
    Bool_t                  fEventChosenByTriggerTrigUtils;             //!
    Int_t                   fCurrentClusterTriggered;                   //
    Int_t                   fCurrentClusterTriggeredTrigUtils;          //
    Int_t                   fCurrentClusterTriggerBadMapResult;         //
    Int_t                   fDoDebugOutput;                             //

    map<Int_t,Int_t>   fMapClusterToTriggered;                     //! connects a given cluster ID with trigger bad map
    map<Int_t,Int_t>   fMapClusterToTriggerMap;                    //! connects a given cluster ID with trigger bad map

    TList*                  fOutputList;                                //!
    TH1I*                   fHist_Event_Accepted;                       //!
    TH1I*                   fHist_Triggered_wEventFlag;                 //!
    TH1I*                   fHist_Cluster_Accepted;                     //!
    TH1I*                   fHist_nModues;                              //!
    TH1I*                   fHist_cellID_All;                           //!
    TH1I*                   fHist_cellID_isAccepted;                    //!
    TH1I*                   fHist_relID0_All;                           //!
    TH1I*                   fHist_relID0_cellIDwasAccepted;             //!
    TH1I*                   fHist_relID0_isAccepted;                    //!

    ClassDef(AliCaloTriggerMimicHelper, 5);
};

#endif
