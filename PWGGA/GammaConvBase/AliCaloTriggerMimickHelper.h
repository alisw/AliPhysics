#ifndef AliCaloTriggerMimickHelper_H
#define AliCaloTriggerMimickHelper_H

#include "AliAnalysisTaskSE.h"
#include "TChain.h"
#include "AliPHOSGeometry.h"
#include "../PHOSTasks/ClusterSelection/AliPHOSTriggerUtils.h"
#include <iostream>

class AliPHOSTriggerUtils;
using namespace std;

class AliCaloTriggerMimickHelper : public AliAnalysisTaskSE {

  public:
    AliCaloTriggerMimickHelper(const char *name="AliCaloTriggerMimickHelper", Int_t clusterType = 0, Bool_t isMC=kFALSE);
    AliCaloTriggerMimickHelper(const char *name="AliCaloTriggerMimickHelper", Int_t clusterType = 0, Int_t isMC=0);
    enum phosTriggerType{kPHOSAny,kPHOSL0,kPHOSL1low,kPHOSL1med,kPHOSL1high} ;
    //Uncopyable & operator=(const Uncopyable&);

    virtual ~AliCaloTriggerMimickHelper();                            //virtual destructor
    void UserCreateOutputObjects();

    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    void SetLightOutput( Bool_t flag )              { fDoLightOutput = flag                         ; }
    void SetPHOSTrigger(phosTriggerType t=kPHOSL0)  { fPHOSTrigger=t                                ; }
    phosTriggerType GetPHOSTrigger()                { return fPHOSTrigger                           ; }
    void SetRunNumber(int run)                      { fRunNumber=run; fForceRun=kTRUE               ; }             //Use given run number, don't read from data
    void SetEventChosenByTrigger( Bool_t flag )     { fEventChosenByTrigger = flag                  ; }
    Bool_t GetEventChosenByTrigger()                { return fEventChosenByTrigger                  ; }
    AliPHOSGeometry*  GetGeomPHOS()                 { return fGeomPHOS                              ; }

  private:
    AliCaloTriggerMimickHelper (const AliCaloTriggerMimickHelper&);             // not implemented
    AliCaloTriggerMimickHelper & operator=(const AliCaloTriggerMimickHelper&);  // not implemented

    // private methods
    void SetClusterType(Int_t iClusterType)     { fClusterType = iClusterType                 ; }
    void SetTriggerDataOrMC(AliVCluster * clu, Bool_t isMCPhoton);


    // basic variables/objects
    Int_t                   fClusterType;                               //! EMCal(1), PHOS(2) or not running (0)
    Int_t                   fRunNumber;                                 // current run number
    Int_t                   nModules;
    Int_t                   fNMaxPHOSModules;
    Int_t                   nMaxCellsPHOS;


    phosTriggerType         fPHOSTrigger;                               // Kind of PHOS trigger: L0,L1
    AliPHOSTriggerUtils *   fPHOSTrigUtils ;                            //! utils to analyze PHOS trigger
    AliPHOSGeometry*        fGeomPHOS;                                  // pointer to PHOS geometry

    Bool_t                  fDoLightOutput;                             // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    Bool_t                  fForceRun ;                                 // use fixed run number, dont read from data
    Bool_t                  fIsMC ;                                     //! Is this is MC
    Bool_t                  fEventChosenByTrigger;                      //!

    ClassDef(AliCaloTriggerMimickHelper, 1);
};

#endif
