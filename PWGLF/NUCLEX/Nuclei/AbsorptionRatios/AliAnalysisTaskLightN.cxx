/*
 * AliAnalysisTaskFemtoPlotration.cxx
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger
 */

#include "AliAnalysisTaskLightN.h"
//#include "AliLog.h"
ClassImp(AliAnalysisTaskLightN)

AliAnalysisTaskLightN::AliAnalysisTaskLightN()
:AliAnalysisTaskSE()
,fTrackBufferSize(0)
,fMVPileUp(false)
,fEvtCutQA(false)
,fIsMC(false)
,fname()
,fAnalysisParticle()
,fQA()
,fEvtCutsParticle()
,fEvtHistListParticle(0)
,fTrackCutsProton()
,fAntiTrackCutsProton()
,fTrackCutsDeuteron()
,fAntiTrackCutsDeuteron()
,fTrackCutHistListProton(0)
,fTrackCutHistMCListProton(0)
,fAntiTrackCutHistListProton(0)
,fAntiTrackCutHistMCListProton(0)
,fTrackCutHistListDeuteron(0)
,fTrackCutHistMCListDeuteron(0)
,fAntiTrackCutHistListDeuteron(0)
,fAntiTrackCutHistMCListDeuteron(0)
{}

AliAnalysisTaskLightN::AliAnalysisTaskLightN(const char *name,bool isMC)
:AliAnalysisTaskSE(name)
,fTrackBufferSize(0)
,fMVPileUp(false)
,fEvtCutQA(false)
,fIsMC(isMC)
,fname(name)
,fAnalysisParticle()
,fQA()
,fEvtCutsParticle()
,fEvtHistListParticle(0)
,fTrackCutsProton()
,fAntiTrackCutsProton()
,fTrackCutsDeuteron()
,fAntiTrackCutsDeuteron()
,fTrackCutHistListProton(0)
,fTrackCutHistMCListProton(0)
,fAntiTrackCutHistListProton(0)
,fAntiTrackCutHistMCListProton(0)
,fTrackCutHistListDeuteron(0)
,fTrackCutHistMCListDeuteron(0)
,fAntiTrackCutHistListDeuteron(0)
,fAntiTrackCutHistMCListDeuteron(0)
{
    DefineOutput(1, TList::Class());  //Output for the Event Class and Pair Cleaner
    DefineOutput(2, TList::Class());  //Output for the Event Cuts
    DefineOutput(3, TList::Class());  //Output for the Track Cuts (Proton)
    DefineOutput(4, TList::Class());  //Output for the Antitrack Cuts (Proton)
    DefineOutput(5, TList::Class());  //Output for the Track Cuts(Deuteron)
    DefineOutput(6, TList::Class());  //Output for the AntiTrack Cuts (Deuteron)
    if (fIsMC) {
        DefineOutput(7, TList::Class());  //Output for the Track Cut MC Info (Proton)
        DefineOutput(8, TList::Class()); //Output for the AntiTrack Cut MC Info (Proton)
        DefineOutput(9, TList::Class()); //Output for the Track Cut MC Info (Deuteron)
        DefineOutput(10, TList::Class()); //Output for the AntiTrack Cut MC Info (Deuteron)
    }
}

AliAnalysisTaskLightN::~AliAnalysisTaskLightN() {
    if (fAnalysisParticle) {
        delete fAnalysisParticle;
    }
}

void AliAnalysisTaskLightN::UserCreateOutputObjects() {
    
    
    //Initialization for the analyis of all particles
    fAnalysisParticle=new AliLightNAnalysis();
    fAnalysisParticle->SetTrackBufferSize(fTrackBufferSize);
    fAnalysisParticle->SetMVPileUp(fMVPileUp);
    fAnalysisParticle->SetEvtCutQA(fEvtCutQA);
    //Set the event cuts
    if (fEvtCutsParticle) {
        fAnalysisParticle->SetEventCutsParticle(fEvtCutsParticle);
    } else {
        AliFatal("Event cuts missing!");
    }
    
    //Set the track cuts for all particles to be analyzed
    if (fTrackCutsProton) {
        fAnalysisParticle->SetTrackCutsProton(fTrackCutsProton);
    } else {
        AliFatal("TrackCuts missing!");
    }
    if (fAntiTrackCutsProton) {
        fAnalysisParticle->SetAntiTrackCutsProton(fAntiTrackCutsProton);
    } else {
        AliFatal("Antitrack cuts missing");
    }
    if (fTrackCutsDeuteron) {
        fAnalysisParticle->SetTrackCutsDeuteron(fTrackCutsDeuteron);
    } else {
        AliFatal("TrackCuts missing!");
    }
    if (fAntiTrackCutsDeuteron) {
        fAnalysisParticle->SetAntiTrackCutsDeuteron(fAntiTrackCutsDeuteron);
    } else {
        AliFatal("Antitrack cuts missing");
    }
    
    //Initialize the cuts and the histograms
    fAnalysisParticle->Init();
    
    //Outputs for the event selection
    if (fAnalysisParticle->GetQAList()) {
        fQA=fAnalysisParticle->GetQAList();
    } else {
        AliFatal("QA Histograms not available");
    }
    if (fAnalysisParticle->GetEventCutHists()) {
        fEvtHistListParticle=fAnalysisParticle->GetEventCutHists();
    } else {
        AliFatal("Event Cut Histograms not available");
    }
    
    //Outputs for (anti)proton
    if (fAnalysisParticle->GetTrackCutHistsProton()) {
        fTrackCutHistListProton=fAnalysisParticle->GetTrackCutHistsProton();
    } else {
        AliFatal("Event Cut Histograms not available");
    }
    if (fAnalysisParticle->GetAntitrackCutHistsProton()) {
        fAntiTrackCutHistListProton=fAnalysisParticle->GetAntitrackCutHistsProton();
    } else {
        AliFatal("Event Cut Histograms not available");
    }
    PostData(1,fQA);
    PostData(2,fEvtHistListParticle);
    PostData(3,fTrackCutHistListProton);
    PostData(4,fAntiTrackCutHistListProton);
    
    if (fIsMC) {
        if (fTrackCutsProton->GetMCQAHists()) {
            fTrackCutHistMCListProton=fTrackCutsProton->GetMCQAHists();
        } else {
            AliFatal("No Track Cut MC Histograms!");
        }
        if (fAntiTrackCutsProton->GetMCQAHists()) {
            fAntiTrackCutHistMCListProton=fAntiTrackCutsProton->GetMCQAHists();
        } else {
            AliFatal("No Antitrack Cut MC Histograms!");
        }
        PostData(7,fTrackCutHistMCListProton);
        PostData(8,fAntiTrackCutHistMCListProton);
    }
    
    
    //Outputs for (anti)deuterons
    if (fAnalysisParticle->GetTrackCutHistsDeuteron()) {
        fTrackCutHistListDeuteron=fAnalysisParticle->GetTrackCutHistsDeuteron();
    } else {
        AliFatal("Event Cut Histograms not available");
    }
    if (fAnalysisParticle->GetAntitrackCutHistsDeuteron()) {
        fAntiTrackCutHistListDeuteron=fAnalysisParticle->GetAntitrackCutHistsDeuteron();
    } else {
        AliFatal("Event Cut Histograms not available");
    }
    PostData(5,fTrackCutHistListDeuteron);
    PostData(6,fAntiTrackCutHistListDeuteron);
    
    if (fIsMC) {
        if (fTrackCutsDeuteron->GetMCQAHists()) {
            fTrackCutHistMCListDeuteron=fTrackCutsDeuteron->GetMCQAHists();
        } else {
            AliFatal("No Track Cut MC Histograms!");
        }
        if (fAntiTrackCutsDeuteron->GetMCQAHists()) {
            fAntiTrackCutHistMCListDeuteron=fAntiTrackCutsDeuteron->GetMCQAHists();
        } else {
            AliFatal("No Antitrack Cut MC Histograms!");
        }
        PostData(9,fTrackCutHistMCListDeuteron);
        PostData(10,fAntiTrackCutHistMCListDeuteron);
    }
    
}

void AliAnalysisTaskLightN::UserExec(Option_t *) {
    AliAODEvent *Event=static_cast<AliAODEvent*>(fInputEvent);
    if (!Event) {
        AliFatal("No Input Event");
    } else {
        
        //Execute the Analysis ((anti)p and (anti)d)
        fAnalysisParticle->Make(Event);
        if (fAnalysisParticle->GetQAList()) {
            fQA=fAnalysisParticle->GetQAList();
        } else {
            AliFatal("QA Histograms not available");
        }
        if (fAnalysisParticle->GetEventCutHists()) {
            fEvtHistListParticle=fAnalysisParticle->GetEventCutHists();
        } else {
            AliFatal("Event Cut Histograms not available");
        }
        if (fAnalysisParticle->GetTrackCutHistsProton()) {
            fTrackCutHistListProton=fAnalysisParticle->GetTrackCutHistsProton();
        } else {
            AliFatal("Track Cut Histograms not available");
        }
        if (fAnalysisParticle->GetAntitrackCutHistsProton()) {
            fAntiTrackCutHistListProton=fAnalysisParticle->GetAntitrackCutHistsProton();
        } else {
            AliFatal("AntiTrack Cut Histograms not available");
        }
        PostData(1,fQA);
        PostData(2,fEvtHistListParticle);
        PostData(3,fTrackCutHistListProton);
        PostData(4,fAntiTrackCutHistListProton);
        if (fIsMC) {
            if (fTrackCutsProton->GetMCQAHists()) {
                fTrackCutHistMCListProton=fTrackCutsProton->GetMCQAHists();
            } else {
                AliFatal("No Track Cut MC Histograms!");
            }
            if (fAntiTrackCutsProton->GetMCQAHists()) {
                fAntiTrackCutHistMCListProton=fAntiTrackCutsProton->GetMCQAHists();
            } else {
                AliFatal("No Antitrack Cut MC Histograms!");
            }
            PostData(7,fTrackCutHistMCListProton);
            PostData(8,fAntiTrackCutHistMCListProton);
        }
        
        if (fAnalysisParticle->GetTrackCutHistsDeuteron()) {
            fTrackCutHistListDeuteron=fAnalysisParticle->GetTrackCutHistsDeuteron();
        } else {
            AliFatal("Track Cut Histograms not available");
        }
        if (fAnalysisParticle->GetAntitrackCutHistsDeuteron()) {
            fAntiTrackCutHistListDeuteron=fAnalysisParticle->GetAntitrackCutHistsDeuteron();
        } else {
            AliFatal("AntiTrack Cut Histograms not available");
        }
        
        PostData(5,fTrackCutHistListDeuteron);
        PostData(6,fAntiTrackCutHistListDeuteron);
        if (fIsMC) {
            if (fTrackCutsDeuteron->GetMCQAHists()) {
                fTrackCutHistMCListDeuteron=fTrackCutsDeuteron->GetMCQAHists();
            } else {
                AliFatal("No Track Cut MC Histograms!");
            }
            if (fAntiTrackCutsDeuteron->GetMCQAHists()) {
                fAntiTrackCutHistMCListDeuteron=fAntiTrackCutsDeuteron->GetMCQAHists();
            } else {
                AliFatal("No Antitrack Cut MC Histograms!");
            }
            PostData(9,fTrackCutHistMCListDeuteron);
            PostData(10,fAntiTrackCutHistMCListDeuteron);
        }
    }
}
