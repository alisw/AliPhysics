///////////////////////////////////////////////////////////////////////////////
//                         AddTaskJetFlowToyMC                               //
// Author: Redmer A. Bertens, Utrecht University, 2013/4, rbertens@cern.ch   //
///////////////////////////////////////////////////////////////////////////////

/* AddTask macro for jet flow toy mc task
 * task uses an afterburner to tune vn in the pico track 
 * selection which can be used by a jet finder 
 * 
 * task can either generate new tracks or use an existing event and 
 * add vn using the afterburner technique
 *
 * the namespace AliAnalysisTaskJetFlowMC contains a number of useful
 * functions which may be called when doing a local analysis (on an analysis 
 * train one would like to avoid calling macro's from within this macro, but 
 * rather load the macros separately using the namespace functions as examples)
*/

class AliAnalysisDataContainer;
class AliAnalysisTaskJetFlowMC;
class AliGenerator;

AliAnalysisTaskJetFlowMC* AddTaskJetFlowMC(
  const char *outputTracks      = "JetFlowToyMC",
  const char *inputTracks       = "PicoTracks",
  const char *name              ="AliAnalysisTaskJetFlowMC",
  Bool_t doQA                   = kFALSE,
  Bool_t doDecay                = kFALSE,       // be sure to load pythia libs
  Bool_t doEmbedding            = kFALSE,       // not to be used on train
  Double_t ptHardPythiaMin      = 0.,           // pt hard bin lower bound
  Double_t ptHardPythiaMax      = 10.           // pt hard bin upper bound

  )
{  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)                             return 0x0;
  if (!mgr->GetInputEventHandler())     return 0x0;
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":";
  fileName += name;
        // create the task
  AliAnalysisTaskJetFlowMC *task = new AliAnalysisTaskJetFlowMC(name, doQA);
  task->SetTracksOutName(outputTracks);
  task->SetTracksInName(inputTracks);
        // connect input and output
  mgr->AddTask(task);
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  if(doQA) {
      // this task only produces output when the qa flag is set to true
      mgr->ConnectOutput (task, 1, mgr->CreateContainer(Form("%s_container", fileName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  }
  // to decay tracks using as a default pythia
  if(doDecay) task->SetDecayer(TaskJetFlowMC::GetDecayer(kTRUE));
  // to embed pythia jets to the events
  if(doEmbedding) AliJetEmbeddingFromGenTask* eTask = TaskJetFlowMC::EmbedGeneratedJets(TaskJetFlowMC::GetPythiaGenerator(2760., ptHardPythiaMin, ptHardPythiaMax), outputTracks);
  return task;
}

namespace TaskJetFlowMC {

    TF1* GetSpectrum() {
        // full spectrum used for ALICE SIMULATION PLOTS ALI-SIMUL-75145 ALI-SIMUL-75171
        // combination of boltzmann spectrum and hard jet spectrum
        TF1* fspectrum = new TF1("fspectrum", "[0]*(TMath::Power([1], 2)*x*TMath::Exp(-[1]*x))+(x>1)*[2]*(1.17676e-01*TMath::Sqrt(0.1396*0.1396+x*x)*TMath::Power(1.+1./[3]/8.21795e-01*TMath::Sqrt(0.1396*0.1396+x*x),-1.*[3]))*(1/(1 + TMath::Exp(-(x - [4])/[5])))", .2, 200.);
        fspectrum->SetParameters(2434401791.20528, 2.98507, 10069622.25117, 5.50000, 2.80000, 0.20000);
        return fspectrum;   
    }

    TF1* GetThermalSpectrum() {
        // pure boltzmann part of thermal spectrum
        TF1* boltzmann = new TF1("boltzmann", "[0]*(TMath::Power([1], 2)*x*TMath::Exp(-[1]*x))");
        boltzmann->SetParameters(2434401791.20528, 2.98507);
        return boltzmann;
    }

    TVirtualDecayer* GetDecayer(Bool_t local = kTRUE) {
        // setup a decayer
        if(local) gSystem->Load("$ALICE_ROOT/lib/tgt_linuxx8664gcc/libpythia6");
        TPythia6Decayer* decayer = new TPythia6Decayer();
        decayer->SetForceDecay(TPythia6Decayer::kHardonicD);
        return decayer;
    }

    AliGenerator* GetPythiaGenerator(
            Float_t e_cms = 2760.,
            Double_t ptHardMin = 0., Double_t ptHardMax = 11.,
            Int_t tune = 2, Int_t cr = 1) {
        // setup a pythia6 generator
        gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCGenPythia.C");
        return AddMCGenPythia(e_cms, ptHardMin, ptHardMax, tune, cr);
    }

    AliJetEmbeddingFromGenTask* EmbedGeneratedJets(AliGenerator* gen, const char* outputTracks) {
        // generate pythia evnets on the fly and embed them to the event
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbeddingFromGen.C");
        return AddTaskJetEmbeddingFromGen(
                gen,                            // generator
                outputTracks,                   // tracks name
                "JetEmbeddingFromGenTask",      // task name
                0.150,                          // min pt
                200,                            // max pt
                -0.9,                           // min eta
                0.9,                            // max eta
                0,                              // min phi
                TMath::Pi() * 2,                // max phi
                kFALSE,                         // copy tracks
                kTRUE);                         // do qa plots
    }

}
