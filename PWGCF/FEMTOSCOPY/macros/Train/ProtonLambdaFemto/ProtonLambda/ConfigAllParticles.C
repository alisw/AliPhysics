#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoV0PairCut.h"
#include "AliFemtoV0TrackPairCut.h"
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoCorrFctnNonIdDR.h"
#include "AliFemtoAvgSepCorrFctn.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliESDtrack.h"
#endif

const int nMult = 10;
int runMult[nMult] = {1, 1, 1, 1, 1, 1, 0, 0, 0, 0};
int multBins[nMult+1] = {0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900};

AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis);

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(bool mcAnalysis=false, bool sepCuts=false, int year=2015)
{
    // create analysis managers
    AliFemtoManager* Manager=new AliFemtoManager();
    AliFemtoModelManager *modelMgr = new AliFemtoModelManager();
    
    // add event reader
    if(year==2015)
    {
        AliFemtoEventReaderAODMultSelection* Reader2015 = GetReader2015(mcAnalysis);
        Manager->SetEventReader(Reader2015);
    }
    else if(year==2011)
    {
        AliFemtoEventReaderAODChain* Reader2011 = GetReader2011(mcAnalysis);
        Manager->SetEventReader(Reader2011);
    }
    
    // declare necessary objects
    AliFemtoVertexMultAnalysis    *femtoAnalysis[nMult];
    AliFemtoModelCorrFctn         *modelCF[nMult];
    
    // setup analysis
    int anIter = 0;
    for (int imult=0; imult<nMult; imult++)
    {
        if (!runMult[imult]) continue;
        
        // create new analysis
        AliFemtoVertexMultAnalysis *analysis = new AliFemtoVertexMultAnalysis(8, -8.0, 8.0, 4, multBins[imult], multBins[imult+1]);
        analysis->SetNumEventsToMix(10);
        analysis->SetMinSizePartCollection(1);
        analysis->SetVerboseMode(kTRUE);
        femtoAnalysis[imult] = analysis;
        
        // set event cut
        AliFemtoBasicEventCut *eventCut = new AliFemtoBasicEventCut();
        eventCut->SetVertZPos(-8,8);
        eventCut->SetEPVZERO(-TMath::Pi()/2.,TMath::Pi()/2.);
        femtoAnalysis[imult]->SetEventCut(eventCut);
        
        
        // set other analysis options and cuts
        femtoAnalysis[imult]->SetFirstParticleCut(new AliFemtoMCTrackCut());
        femtoAnalysis[imult]->SetSecondParticleCut(new AliFemtoMCTrackCut());
        femtoAnalysis[imult]->SetEnablePairMonitors(true);
        femtoAnalysis[imult]->SetPairCut(new AliFemtoDummyPairCut());
        
        // add correlation function on model data
        modelCF[imult] = new AliFemtoModelCorrFctn(Form("CF_qinv_Model_M%i",imult),100,0.0,2.0);
        modelCF[imult]->SetFillkT(true);
        modelCF[imult]->ConnectToManager(modelMgr);
        femtoAnalysis[imult]->AddCorrFctn(modelCF[imult]);
        femtoAnalysis[imult]->SetEnablePairMonitors(true);
        
        Manager->AddAnalysis(femtoAnalysis[imult]);
        
        
    }
    return Manager;
}

AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis)
{
    AliFemtoEventReaderAODMultSelection* Reader = new AliFemtoEventReaderAODMultSelection();
    Reader->SetFilterBit(7);
    Reader->SetReadV0(1);
    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
    Reader->SetEPVZERO(kTRUE);
    Reader->SetCentralityFlattening(kTRUE);
    if(mcAnalysis) Reader->SetReadMC(kTRUE);
    
    return Reader;
}

AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis)
{
    AliFemtoEventReaderAODChain* Reader = new AliFemtoEventReaderAODChain();
    Reader->SetFilterBit(7);
    Reader->SetReadV0(1);
//    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
    Reader->SetEPVZERO(kTRUE);
    Reader->SetCentralityFlattening(kTRUE);
    if(mcAnalysis) Reader->SetReadMC(kTRUE);
    
    return Reader;
}

