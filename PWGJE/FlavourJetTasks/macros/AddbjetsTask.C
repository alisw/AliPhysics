
class AliAnalysisDataContainer;

#include "TString.h"
#include "AliAnalysisTaskbjets.h"
#include "AliAnalysisManager.h"
//#include "AliAnalysisTaskParticleInJet.h"

AliAnalysisTaskbjets* AddbjetsTask(TString name = "name")

/*TString ntracksMC = "mcparticles",
 TString nclusters = "",
 //const char *ntracksMC = "tracksMC";
 //const char *ntracks = "Tracks";
 TString ntracks = "tracks",
 TString njets = "Jet_AKTChargedR040_tracks_pT0150_E_scheme",
 TString njetsMC = "Jet_AKTChargedR040_mcparticles_pT0150_E_scheme",
 Double_t jetradius =0.4,
 TString type = "TPC",
 TString nrho = "",
 Bool_t isMC = kTRUE,
 TString taskname = "AliAnalysisTaskEmcalJetBJetTaggingIP",
 TString nrhoMC = "RhoMC",*/



{
 TString ntracksMC = "mcparticles";
 TString nclusters = "";
 //const char *ntracksMC = "tracksMC";
 //const char *ntracks = "Tracks";
 TString ntracks = "tracks";
 TString njets = "Jet_AKTChargedR040_tracks_pT0150_E_scheme";
 TString njetsMC = "Jet_AKTChargedR040_mcparticles_pT0150_E_scheme";
 Double_t jetradius =0.4;
 TString type = "TPC";
 TString nrho = "";
 Bool_t isMC = kTRUE;
 TString taskname = "AliAnalysisTaskEmcalJetBJetTaggingIP";
 TString nrhoMC = "RhoMC";
 
 
 
 //Printf("Check done %i",__LINE__);
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    
    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    
    
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskbjets* task = new AliAnalysisTaskbjets(name.Data());   
    if(!task) return 0x0;
    
    
    // Setup input containers
    //==============================================================================
    Printf("%s :: Setting up input containers.",taskname.Data());
    AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks.Data());
    AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters.Data());
    TString strType(type.Data());
    AliJetContainer *jetCont = task->AddJetContainer(njets.Data(),strType,jetradius);
    
    Printf("Check done %i",__LINE__);
    
    if(jetCont) {
        jetCont->SetRhoName(nrho.Data());
        jetCont->ConnectParticleContainer(trackCont);
        jetCont->ConnectClusterContainer(clusterCont);
        //DefineCutsTaskpp(jetCont, jetradius);
    }
    Printf("Check done %i",__LINE__);

    if(isMC){
        AliParticleContainer *trackContMC   = task->AddParticleContainer(ntracksMC.Data());
        AliJetContainer *jetContMC = task->AddJetContainer(njetsMC.Data(),strType,jetradius);
        
        if(jetContMC) {
            jetContMC->SetRhoName(nrhoMC.Data());
            jetContMC->ConnectParticleContainer(trackContMC);
            jetContMC->SetIsParticleLevel(kTRUE);
            //jetContMC->SetMaxTrackPt(1000);
            //DefineCutsTaskpp(jetContMC, jetradius);
        }
    }
    
   //=====================================================================================
    Printf("Check done %i",__LINE__);
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    Printf("Check done %i",__LINE__);
    return task;
    
} 
