#ifndef __CINT__
#endif
void runProcess()
{

    TStopwatch timer;
    timer.Start();
    
    Bool_t boolProcess = kTRUE;
    
    gROOT->LoadMacro("PWG2resonancesUtils.C");
    gROOT->LoadMacro("ProcessConfig.C");

    // number of events
    Long64_t numOfEvents = 10000;
    // number of events to skip
    Long64_t numOfEventsSkip = 330000*0;

    TString macro;
    macro = "AliRsnTrain.C";

//     boolProcess = RunLocaly(macro,numOfEvents,numOfEventsSkip);
    boolProcess = RunOnProof(macro,numOfEvents,numOfEventsSkip);
//     boolProcess = RunOnAliEn(macro,numOfEvents,numOfEventsSkip);

    timer.Stop();
    timer.Print();

    Info("","Process ended %s.",((boolProcess)?"OK":"with an ERROR"));
    
}

