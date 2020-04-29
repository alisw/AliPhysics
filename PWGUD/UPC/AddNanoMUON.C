AliAnalysisTaskNanoMUON* AddNanoMUON(const char* suffix = "")
{
TString name = "NanoMUON";
name += suffix;     
// get the manager via the static access member. since it's static, you don't need
// to create an instance of the class here to call the function
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
// get is MC information
Bool_t fIsMC;
if( mgr->GetMCtruthEventHandler() ) fIsMC = kTRUE;

// by default, a file is open for writing. here, we get the filename
TString fileName = AliAnalysisManager::GetCommonFileName();
fileName += ":NanoMUON";      // create a subfolder in the file
// now we create an instance of your task
AliAnalysisTaskNanoMUON* task = new AliAnalysisTaskNanoMUON(name.Data());   
if(!task) return 0x0;

// set task options
task->SetMC(fIsMC);

// add your task to the manager
mgr->AddTask(task);
// your task needs input: here we connect the manager to your task
mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
// same for the output
TString fRecTreeName = "fRecTree";
fRecTreeName += suffix;
TString fOutputListName = "fOutputList";
fOutputListName += suffix;
TString fGenTreeName = "fGenTree";
fGenTreeName += suffix;
TString fTrgTreeName = "fTrgTree";
fTrgTreeName += suffix;

mgr->ConnectOutput(task,1,mgr->CreateContainer(fRecTreeName.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer,fileName.Data()));
mgr->ConnectOutput(task,2,mgr->CreateContainer(fOutputListName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
mgr->ConnectOutput(task,3,mgr->CreateContainer(fGenTreeName.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer,fileName.Data()));
mgr->ConnectOutput(task,4,mgr->CreateContainer(fTrgTreeName.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer,fileName.Data()));

// in the end, this macro returns a pointer to your task. this will be convenient later on
// when you will run your analysis in an analysis train on grid
return task;
}
