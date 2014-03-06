#ifndef ALIANALYSISTASKESDNUCLEXFILTER_H
#define ALIANALYSISTASKESDNUCLEXFILTER_H


// Create a new AOD starting from the general AOD. This Task can be used also strating 
//from ESD changing the input handler. (Method to be testeted on the grid)
// filtering of the ESD. 
//
// Authors: S. Bufalino (stefania.bufalino@cern.ch)
//          R. Lea      (ramona.lea@cern.ch)
// Based on AliAnalysisTaskESDMuonFilter.h  

#ifndef ALIANALYSISTASKSE_H
#  include "AliAnalysisTaskSE.h"
#  include "AliESDtrack.h"
#  include "AliAODTrack.h"
#  include "AliAODPid.h"
#  include "AliAODNuclExReplicator.h"
//#include "AliAOD3LH.h"

#endif
#ifndef ALIANALYSISCUTS_H
#  include "AliAnalysisCuts.h"
#endif

class AliAnalysisFilter;
class AliESDtrack;
class AliPIDResponse;

class AliAnalysisTaskESDNuclExFilter : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskESDNuclExFilter(Bool_t onlyMuon=kTRUE, Bool_t keepAllEvents=kTRUE, Int_t mcMode=0, Int_t nsigmaTrk1=3 ,Int_t nsigmaTrk2 = 3, Int_t partType1 = 2,Int_t partType2 = 7);
    AliAnalysisTaskESDNuclExFilter(const char* name, Bool_t onlyMuon=kTRUE, Bool_t keepAllEvents=kTRUE, Int_t mcMode=0, Int_t nsigmaTrk1=3 ,Int_t nsigmaTrk2 = 3, Int_t partType1 = 2,Int_t partType2 = 7);
    virtual ~AliAnalysisTaskESDNuclExFilter() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

    virtual void ConvertESDtoAOD();
  

    // Setters
    virtual void SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
    void SetWriteMuonAOD(Bool_t enableMuonAOD){fEnableMuonAOD = enableMuonAOD;}
    void SetWriteDimuonAOD(Bool_t enableDimuonAOD){fEnableDimuonAOD = enableDimuonAOD;}

    void PrintTask(Option_t *option="", Int_t indent=0) const;
 

 private:
    AliAnalysisTaskESDNuclExFilter(const AliAnalysisTaskESDNuclExFilter&);
    AliAnalysisTaskESDNuclExFilter& operator=(const AliAnalysisTaskESDNuclExFilter&);
    void AddFilteredAOD(const char* aodfilename, const char* title, Bool_t toMerge);

    AliAnalysisFilter* fTrackFilter; //  Track Filter
    Bool_t fEnableMuonAOD; // flag for enabling Muon AOD production
    Bool_t fEnableDimuonAOD; // flag for enabling Dimuon AOD production
    Bool_t fOnlyMuon; // flag for disabling branches irrelevant for (most) muon analyses
    Bool_t fKeepAllEvents; // keep even events where there's no muons (to get e.g. unbiased vertex distribution)
    Int_t  fMCMode; // whether and how we're filtering MC data
    
    Int_t fnSigmaTrk1;
    Int_t fnSigmaTrk2;
    Int_t fpartType1;
    Int_t fpartType2;

    AliAODNuclExReplicator* murep;

    AliPIDResponse  *fPIDResponse;                  //! PID response object

  ClassDef(AliAnalysisTaskESDNuclExFilter, 5); // Analysis task for standard ESD filtering
};
 
/* class AliAnalysisNonMuonTrackCuts : public AliAnalysisCuts */
/* { */
/* public: */
/*   AliAnalysisNonMuonTrackCuts(); */
/*   virtual ~AliAnalysisNonMuonTrackCuts() {} */
/*   /\* virtual Bool_t IsSelected(TObject* obj); *\/ */
/*   /\* virtual Bool_t IsSelected(TList*   /\\* list *\\/ ) { return kTRUE; } *\/ */

/*   ClassDef(AliAnalysisNonMuonTrackCuts,1); // Select muon spectrometer tracks */
/* }; */

/* class AliAnalysisNonPrimaryVertices : public AliAnalysisCuts */
/* { */
/* public: */
/*   AliAnalysisNonPrimaryVertices(); */
/*   virtual ~AliAnalysisNonPrimaryVertices() {} */
/*   /\* virtual Bool_t IsSelected(TObject* obj); *\/ */
/*   /\* virtual Bool_t IsSelected(TList*   /\\* list *\\/ ) { return kTRUE; } *\/ */
  
/*   ClassDef(AliAnalysisNonPrimaryVertices,1); // Select primary vertices */
/* }; */

#endif
