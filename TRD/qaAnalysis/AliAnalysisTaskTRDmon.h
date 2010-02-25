//
//
// Analysis task used for TRD monitoring
// 
// Authors:    
//         Anton Andronic <A.Andronic@gsi.de>
//         Ionut Arsene   <I.C.Arsene@gsi.de>


#ifndef ALIANALYSISTASKTRDMON_H
#define ALIANALYSISTASKTRDMON_H

#include "AliAnalysisTask.h"

class AliESDEvent;
class AliESDfriend;
class AliMCEvent;
class TList;
class TH1F;
class TH2F;
class TProfile;

class AliAnalysisTaskTRDmon : public AliAnalysisTask {
 public:
    AliAnalysisTaskTRDmon(const char * name = "Data analysis");
    AliAnalysisTaskTRDmon(const AliAnalysisTaskTRDmon &task);
    ~AliAnalysisTaskTRDmon(){};

    AliAnalysisTaskTRDmon& operator=(const AliAnalysisTaskTRDmon &task);    

    virtual void ConnectInputData(Option_t *);
    virtual void CreateOutputObjects();
    virtual void Exec(Option_t *);
    virtual void Terminate(Option_t *);
    //    Bool_t IsRunningTerminate() const { return TestBit(kIsRunningTerminate); };
    //    void SetRunTerminate(Bool_t runTerminate) {SetBit(kIsRunningTerminate, runTerminate); };
    void SetTriggerName(const Char_t* triggerName) {fEventTriggerName=triggerName;}
    void SetIsCollisionEvent(Bool_t flag=kTRUE) {fIsCollisionEvent=flag;}
    const TString& GetTriggerName() const {return fEventTriggerName;}
    Bool_t GetIsCollisionEvent() const {return fIsCollisionEvent;}

 private:
    //    enum{kIsRunningTerminate = BIT(14)};   // not less
    AliESDEvent *fESD;               //!
    AliESDfriend *fESDfriend;        //!
    TString fEventTriggerName;       //  trigger class for the events to be analyzed (all events analyzed if empty string)
    Bool_t fIsCollisionEvent;        //  flag to determine if analyzed events are collisions or not
    TObjArray *fOutStorage;          //!
    TH1F *fHzvert1;                  //!
    TH1F *fHzvert2;                  //!
    TH1F *fHntracks;                 //!
    TH1F *fHntracks2;                //!
    TH1F *fHntracks3;                //!
    TH1F *fHdca;                     //!
    TH1F *fHdcaz;                    //!
    TH1F *fHpt;                      //!
    TH1F *fHpt2;                     //!
    TH1F *fHpt3;                     //!
    TH1F *fHpt3n;                    //!
    TH1F *fHpt4;                     //!
    TH1F *fHpt4n;                    //!
    TH1F *fHtheta;                   //!
    TH1F *fHphi;                     //!
    TH1F *fHtpccl;                   //!
    TH1F *fHtpccl2;                  //!
    TH2F *fHdedxp;                   //! dEdx-p TPC
    TH2F *fHetaphi;                  //! eta-phi tracks TPC
    TH2F *fHetancl;                  //! eta-Nclusters TPC
    TH2F *fHphincl;                  //! phi-Nclusters TPC
    TH1F *fHtrdtr;                   //!
    TH1F *fHtrdtr2;                  //!
    TProfile *fHph;                  //! <PH>
    TH2F *fHph2d;                    //! PH 2d
    TH1F *fHncltrkl;                 //!
    TH1F *fHntrkl;                   //!
    TH1F *fHsm;                      //!
    TH2F *fHetantr;                  //! Ntracklets-eta
    TH2F *fHphintr;                  //! Ntracklets-phi
    TH1F *fHcltime;                  //!
    TH1F *fHcldiff;                  //!
    TH2F *fHxyA;                     //! x-y side A (TPC)
    TH2F *fHxyC;                     //! x-y side C (TPC)
    TH1F *fHnFriendTracks;           //!
    TH1F *fHnCalibObjects;           //!

    ClassDef(AliAnalysisTaskTRDmon, 1)
};
#endif
	
