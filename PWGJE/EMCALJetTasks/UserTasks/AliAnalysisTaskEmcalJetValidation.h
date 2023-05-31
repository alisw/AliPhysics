/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskEmcalJetValidation_H
#define AliAnalysisTaskEmcalJetValidation_H

#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskSE.h"
//#include "AliPIDResponse.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"
#include "AliFiducialCut.h"
#include "AliEMCALRecoUtils.h"
#include <map>
#include <string>
#include <vector>

class TList;
class AliESDEvent;
class AliEmcalJet;
class AliAODEvent;
class TH1F;
class TH1D;
class TH2F;
class TH2D;
class AliEmcalList;
class AliEmcalJetFinder;
class AliMCEvent;
class AliFJWrapper;
class AliMultSelection;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliVParticle;
class AliLog;
class AliOADBContainer;
class AliESDtrackCuts;

class AliAnalysisTaskEmcalJetValidation : public AliAnalysisTaskSE
//class AliAnalysisTaskEmcalJetValidation : public AliAnalysisTaskEmcalJet
{
public:
                            AliAnalysisTaskEmcalJetValidation();
                            AliAnalysisTaskEmcalJetValidation(const char *name);
    virtual                 ~AliAnalysisTaskEmcalJetValidation();

    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);

    static AliAnalysisTaskEmcalJetValidation* AddTask(TString suffix = "", UInt_t trigger= AliVEvent::kINT7, TString jsonconfigfile="", Bool_t readMC=kFALSE);

    void                    ExecOnceLocal();
    void                    SetJetR(Double_t jr){  fJetR = jr; };                                //sets jet radius
    void                    SetAlgo(AliJetContainer::EJetAlgo_t algo){ fJetAlgo = algo; }        // sets  antikt/kt
    void                    SetGhostArea(Double_t ga){ fGhostArea = ga; }                        // sets ghost area
    void                    SetRecoScheme(AliJetContainer::ERecoScheme_t rs){ fRecoScheme =rs;}  // recombination scheme
    void                    InitFromJson(TString filename);                                      //initialization from json file



    std::string             GetJsonString(const char* jsonFileName, const char* section, const char* key);
    int                     GetJsonBool(const char* jsonFileName, const char* section, const char* key);
    int                     GetJsonInteger(const char* jsonFileName, const char* section, const char* key);
    float                   GetJsonFloat(const char* jsonFileName, const char* section, const char* key);
    float*                  GetJsonArray(const char* jsonFileName, const char* section, const char* key, int& size);
    float** GetJsonMatrix(const char* jsonFileName, const char* section, const char* key, int& size1, int& size2);

private:
    AliESDEvent*            fESD;           //! input event
    TList*                  fOutputList;    //! output list
    TH1F*                   fHistJetPt;     //! jet Pt
    TH1F*	                  fHistNEvents;   //! histogram for total number of events
    TH1F*                   fHistNEventVtx; //! event vertex distribution
    TH1F*                   fHistJetPhi;    //! jet Phi
    TH1F*                   fHistJetEta;    //! jet Eta
    TH1F*                   fHistTrackPt;   //! track Pt
    TH1F*                   fHistTrackPhi;  //! track Phi
    TH1F*                   fHistTrackEta;  //! track Eta
    TH1F*                   fHistNTracksAll; //! total number of esd tracks
    TH1F*                   fHistMinCrossedRowsTPC; //! number of esd tracks satisfying the trackcut MinNCrossedRowsTPC(70)
    TH1F*                   fHistMaxChi2PerClusterTPC; //! for MaxChi2PerClusterTPC(4.0) track cut
    TH1F*                   fHistRatioCrossedRowsOverFindableCLustersTPC; //! for SetMinRatioCrossedRowsOverFindableClustersTPC(0.8) track cut
    TH1F*                   fHistMaxChi2PerClusterITS; //! for MaxChi2PerClusterITS(36.0) track cut

    AliFJWrapper*           fFastJetWrapper;  //! utility

    AliESDtrackCuts*        fTrackCuts;  //! track cuts

    Bool_t fInitializedLocal;       //!  flag which marks the first access to  ExecOnceLocal()
    Double_t fMinPt;                //  minimum track/jet pt
    Double_t fJetEtaRange;          //  fiducial cut on jets
    Double_t fJetR;                 //  fiducial cut on jets

    AliJetContainer::EJetAlgo_t fJetAlgo;   //  antikt/kt
    Double_t fGhostArea;                    //  ghost area
    AliJetContainer::ERecoScheme_t fRecoScheme;  // recombination scheme


    AliAnalysisTaskEmcalJetValidation(const AliAnalysisTaskEmcalJetValidation&); // not implemented
    AliAnalysisTaskEmcalJetValidation& operator=(const AliAnalysisTaskEmcalJetValidation&); // not implemented

    ClassDef(AliAnalysisTaskEmcalJetValidation, 18);
};

#endif
