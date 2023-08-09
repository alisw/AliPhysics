#ifndef AliAnalysisTaskMultspec_evbased_H
#define AliAnalysisTaskMultspec_evbased_H

//class AliESDpid;
//class AliESDEvent;

#include "TString.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "THistManager.h"
#include "TH1D.h"
#include "TH2D.h"
#include "AliMCEvent.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskMultspec.h"

class AliAnalysisTaskMultspec_evbased : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskMultspec_evbased();
    AliAnalysisTaskMultspec_evbased(const char *name, int particle, TString lExtraOptions = "", bool ismc = kTRUE, bool removePythiaGen = kTRUE);
    virtual ~AliAnalysisTaskMultspec_evbased();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
    
    void SetIsMC(bool mc) {fisMC=mc;}
    void SetParticle(int part) {if(part<kallparts_Mult) fpart = part; else{printf("Invalid particle. Aborting..."); return;}}
    void SetRemovePythiaGen(bool pythiagen) {fRemovePythiaGen=pythiagen;}

  private:
    
    //outputs
    THistManager *fHistos;                               //!<! Output histos
    
    //fillers
    V0filler_Mult *ffillV0 = nullptr;                               //!<! Transient V0filler
    Cascfiller_Mult *ffillCasc = nullptr;                           //!<! Transient Cascfiller

    //objects retreived from input handler
    AliPIDResponse *fPIDResponse;                              //!
    
    //AliEventCuts object
    AliEventCuts fEventCuts;                                   //

    //pile-up rejection flag
    bool fRejectPileupEvts;                                    //

    //decide if to store cascade tree
    int fpart;                                                 // which particle are we analyzing?

    //flag for MC handeling
    bool fisMC;                                                //

    //flag for remove Pythia events
    bool fRemovePythiaGen;  

    //functions to allow flushing part of code out of UserExec
    bool ApplyCut(int parttype);
    void DataPosting();
    
    AliAnalysisTaskMultspec_evbased(const AliAnalysisTaskMultspec_evbased&);            // not implemented
    AliAnalysisTaskMultspec_evbased& operator=(const AliAnalysisTaskMultspec_evbased&); // not implemented

    ClassDef(AliAnalysisTaskMultspec_evbased, 1);


};

#endif
