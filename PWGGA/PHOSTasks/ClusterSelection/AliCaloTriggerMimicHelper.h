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
    void SetMinEnergyToTrigger(Double_t flag)               { minEnergyToTrigger = flag                     ; }
    void SetMinCellsToTrigger(Int_t flag)                   { minCellsToTrigger = flag                      ; }
    void SetEfficiencyChoiceOption_TriggerHelper(Int_t flag) { fEfficiencyChoiceOption_TriggerHelper = flag ; }
    Int_t GetTriggerHelperRunMode()                         { return fTriggerHelperRunMode                  ; }
    TList* GetTriggerMimicHelperHistograms()                { return fOutputList                            ; }
    Bool_t HasEventFlagPassed()                             { return fEventFlagPassed                       ; }
    Int_t IsClusterIDTriggered(Int_t ClusterID)             { return fMapClusterIDToHaveTriggered[ClusterID]      ; }//return 0 for not triggered clusters
    Int_t IsClusterIDBadMapTrigger(Int_t ClusterID)         { return fMapClusterIDToIsInTriggerMap[ClusterID]     ; }//return 0 for bad clusters
    Int_t IsTriggeredClusterIDInBadDDL(Int_t ClusterID)     { return fMapTriggeredClusterInBadDDL[ClusterID]; } //2 return Bad DDLs, 1 include maybe bad DDLs as well

  private:
    AliCaloTriggerMimicHelper (const AliCaloTriggerMimicHelper&);             // not implemented
    AliCaloTriggerMimicHelper & operator=(const AliCaloTriggerMimicHelper&);  // not implemented

    // private methods
    void SetClusterType(Int_t iClusterType)     { fClusterType = iClusterType                 ; }
    void SetTriggerDataOrMC(AliVCluster * clu, Bool_t isMCPhoton);
    Int_t WhichDDL(Int_t module=0, Int_t cellx=0);
    Int_t WhichTRU(Int_t cellx=0, Int_t cellz=0);
    Int_t WhichTRUChannel(Int_t cellx, Int_t cellz, Int_t &chX, Int_t &chZ);
    Int_t IsDDLBad(Int_t iDDL=0, Int_t iRun=0);                         //returns 0 for good DDLs, 1 for maybe bad DDLs and 2 for bad DDLs


    // basic variables/objects
    TString                 fNameOfClassObject;                         // name of this class object
    Int_t                   fClusterType;                               // EMCal(1), PHOS(2) or not running (0)
    Int_t                   fRunNumber;                                 // current run number
    Int_t                   nModules;
    Int_t                   fNMaxPHOSModules;
    Int_t                   nMaxCellsPHOS;
    Int_t                   maxRows;
    Int_t                   maxColumns;
    Int_t                   maxCellsModule;
    Int_t                   startDDLNumber;
    Int_t                   endDDLNumber;
    Int_t                   maxNumberOfDDLs;
    Int_t                   startTRU_Number;
    Int_t                   endTRU_Number;
    Int_t                   maxNumberOfTRUs;



    phosTriggerType         fPHOSTrigger;                               // Kind of PHOS trigger: L0,L1
    AliPHOSTriggerUtils *   fPHOSTrigUtils ;                            //! utils to analyze PHOS trigger
    AliPHOSGeometry*        fGeomPHOS;                                  //! pointer to PHOS geometry

    Bool_t                  fDoLightOutput;                             // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    Bool_t                  fForceRun ;                                 // use fixed run number, dont read from data
    Bool_t                  fIsMC ;                                     // Is this is MC
    Int_t                   fTriggerHelperRunMode;                      //0 uses L0 events, 1 uses events which have no L0, 3 uses all events
    Bool_t                  fEventFlagPassed;                           //Only true if the event selection has passed
    Bool_t                  fEventChosenByTrigger;                      //!
    Bool_t                  fEventChosenByTriggerTrigUtils;             //!
    Int_t                   fCurrentClusterTriggered;                   //
    Int_t                   fCurrentClusterTriggeredTrigUtils;          //
    Int_t                   fCurrentClusterTriggerBadMapResult;         //
    Int_t                   fCurrentTriggeredClusterInBadDDL;           //
    Int_t                   fDoDebugOutput;                             //
    Double_t                minEnergyToTrigger;                         //
    Int_t                   minCellsToTrigger;                         //
    Int_t                   fEfficiencyChoiceOption_TriggerHelper;      //Sets, which Turnoncurves will be used

    map<Int_t,Int_t>        fMapClusterIDToHaveTriggered;               //! connects a given cluster ID with the Information if Cluster has triggered
    map<Int_t,Int_t>        fMapClusterIDToIsInTriggerMap;              //! connects a given cluster ID with the information, if it is in trigger good map

    map<Int_t,Int_t>        fMapTriggeredClusterInBadDDL;               //! connects a given cluster ID with trigger bad map

    TList*                  fOutputList;                                //!
    Bool_t                  fdo_fHist_Event_Accepted;                    //Turn On or Off if Histograms are created and used
    TH1I*                   fHist_Event_Accepted;                       //!
    Bool_t                  fdo_fHist_Triggered_wEventFlag;              //Turn On or Off if Histograms are created and used
    TH1I*                   fHist_Triggered_wEventFlag;                 //!
    Bool_t                  fdo_fHist_Cluster_Accepted;                  //Turn On or Off if Histograms are created and used
    TH1I*                   fHist_Cluster_Accepted;                     //!
    Bool_t                  fdo_fHist_cellID;                            //Turn On or Off if Histograms are created and used
    TH1I*                   fHist_cellID_All;                           //!
    TH1I*                   fHist_cellID_isAccepted;                    //!
    Bool_t                  fdo_fHist_relID;                             //Turn On or Off if Histograms are created and used
    TH1I*                   fHist_relID0_All;                           //!
    TH1I*                   fHist_relID0_cellIDwasAccepted;             //!
    TH1I*                   fHist_relID0_isAccepted;                    //!
    Bool_t                  fdo_fHist_GammaClusE;                        //Turn On or Off if Histograms are created and used
    TH1D*                   fHist_GammaClusE_Trig;                      //!
    TH1D*                   fHist_GammaClusE_notTrig;                   //!
    Bool_t                  fdo_TriggeredClusters_ColumnVsRow_overThresh;//Turn On or Off if Histograms are created and used
    TH2I**                  fHist_TriggeredClusters_ColumnVsRow_overThresh;//!
    Bool_t                  fdo_TriggeredClusters_ColumnVsRow_underThresh;//Turn On or Off if Histograms are created and used
    TH2I**                  fHist_TriggeredClusters_ColumnVsRow_underThresh;//!
    Bool_t                  fdo_Any_4x4_Distance;                       //Turn On or Off if Histograms are created and used
    Bool_t                  fdo_4x4_Distance_All;                       //Turn On or Off if Histograms are created and used
    TH2I*                   fHist_4x4_Distance_All;                     //!
    Bool_t                  fdo_Tr4x4_Distance_Triggered;               //Turn On or Off if Histograms are created and used
    TH2I*                   fHist_Tr4x4_Distance_Triggered;             //!
    Bool_t                  fdo_Tr4x4_Distance_notTriggered;            //Turn On or Off if Histograms are created and used
    TH2I*                   fHist_Tr4x4_Distance_notTriggered;          //!
    Bool_t                  fdo_Any_TRU;                                //Turn On or Off if Histograms are created and used
    Bool_t                  fdo_ClusEVsTiming_TRU;                      //Turn On or Off if Histograms are created and used
    TH2D***                 fHist_ClusEVsTiming_TRU;                    //!
    Bool_t                  fdo_ClusEVsTiming_TRU_Trig;                 //Turn On or Off if Histograms are created and used
    TH2D***                 fHist_ClusEVsTiming_TRU_Trig;               //!
    Bool_t                  fdo_ClusEVsTiming_TRU_notTrig;              //Turn On or Off if Histograms are created and used
    TH2D***                 fHist_ClusEVsTiming_TRU_notTrig;            //!
    Bool_t                  fdo_TRU_Numbers;                            //Turn On or Off if Histograms are created and used
    TH1I**                  fHist_TRU_Numbers;                          //!
    Bool_t                  fdo_TRU_Channels;                           //Turn On or Off if Histograms are created and used
    TH1I**                  fHist_TRU_Channels;                         //!
    Bool_t                  fdo_TRU_ChannelsXZ;                         //Turn On or Off if Histograms are created and used
    TH2I**                  fHist_TRU_ChannelsXZ;                       //!

    Double_t                fEnergyThreshold_ColumnVsRow;

    ClassDef(AliCaloTriggerMimicHelper, 14);
};

#endif
