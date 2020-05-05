#ifndef AliTriggerMimickHelper_H
#define AliTriggerMimickHelper_H

#include "AliAnalysisTaskSE.h"
#include "AliPHOSTriggerUtils.h"

using namespace std;

class AliTriggerMimickHelper : public AliAnalysisTaskSE {

  public:
    AliTriggerMimickHelper(const char *name="CaloTrackMatcher_0_0", Int_t clusterType = 0, Bool_t isMC);
    enum phosTriggerType{kPHOSAny,kPHOSL0,kPHOSL1low,kPHOSL1med,kPHOSL1high} ;
    //Uncopyable & operator=(const Uncopyable&);

    virtual ~AliTriggerMimickHelper();                            //virtual destructor
    void UserCreateOutputObjects();

    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    void SetLightOutput( Bool_t flag )              { fDoLightOutput = flag                         ; }
    void SetPHOSTrigger(phosTriggerType t=kPHOSL0)  { fPHOSTrigger=t                                ; }
    phosTriggerType GetPHOSTrigger()                { return fPHOSTrigger                           ; }
    void SetRunNumber(int run)                      { fRunNumber=run; fForceRun=kTRUE               ; }             //Use given run number, don't read from data
    void SetEventChosenByTrigger( Bool_t flag )     { fEventChosenByTrigger = flag                  ; }
    Bool_t GetEventChosenByTrigger()                { return fEventChosenByTrigger                  ; }

  private:
    AliTriggerMimickHelper (const AliTriggerMimickHelper&);             // not implemented
    AliTriggerMimickHelper & operator=(const AliTriggerMimickHelper&);  // not implemented

    // private methods
    void SetClusterType(Int_t iClusterType)     { fClusterType = iClusterType                 ; }
    void SetTriggerDataOrMC(AliVCluster * clu, Bool_t isMCPhoton);


    // basic variables/objects
    Int_t                   fClusterType;                               //! EMCal(1), PHOS(2) or not running (0)
    Int_t                   fRunNumber;                                 // current run number

    phosTriggerType         fPHOSTrigger;                               // Kind of PHOS trigger: L0,L1
    AliPHOSTriggerUtils *   fPHOSTrigUtils ;                            //! utils to analyze PHOS trigger

    Bool_t                  fDoLightOutput;                             // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    Bool_t                  fForceRun ;                                 // use fixed run number, dont read from data
    Bool_t                  fIsMC ;                                     //! Is this is MC
    Bool_t                  fEventChosenByTrigger;                      //!

    ClassDef(AliTriggerMimickHelper, 1);
};

#endif
