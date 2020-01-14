#include "TList.h"
#include "TChain.h"
#include "TH1F.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskVdmStability.h"


AliAnalysisTaskVdmStability::AliAnalysisTaskVdmStability():
    AliAnalysisTaskSE(),
    fAOD{0},
    fOutputList{0},
    fHist{0}
{
    // ROOT IO constructor, don't allocate memory here!
}
AliAnalysisTaskVdmStability::AliAnalysisTaskVdmStability(const char* taskname):
    AliAnalysisTaskSE(taskname),
    fAOD{0},
    fOutputList{0},
    fHist{0}
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
AliAnalysisTaskVdmStability::~AliAnalysisTaskVdmStability()
{
    
}
//defíne output
void AliAnalysisTaskVdmStability::UserCreateOutputObjects()
{
    // create a new TList that OWNS its objects
    fOutputList = new TList();
    fOutputList->SetOwner(true);

    // create our histo and add it to the list
    fHist = new TH1F("fHist", "fHist", 100, 0, 100);
    fOutputList->Add(fHist);

    // add the list to our output file
    PostData(1,fOutputList);
}

//Event loop
void AliAnalysisTaskVdmStability::UserExec(Option_t*)
{
    // get an event from the analysis manager
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());

    // check if there actually is an event
    if(!fAOD)
        ::Fatal("AliAnalysisTaskMyTask::UserExec", "No AOD event found, check the event handler.");

    // let's loop over the tracks and fill our histogram

    // first we get the number of tracks
    int iTracks{fAOD->GetNumberOfTracks()};

    // and then loop over them
    for(int i{0}; i < iTracks; i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track) continue;

        // here we do some track selection
       // if(!track->TestFilterbit(128) continue;

           // fill our histogram
           fHist->Fill(track->Pt());
           }
           // and save the data gathered in this iteration
           PostData(1, fOutputList);
}
           
 void AliAnalysisTaskVdmStability::Terminate(Option_t*)
{
               
}
