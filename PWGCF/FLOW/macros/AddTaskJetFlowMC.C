///////////////////////////////////////////////////////////////////////////////
//                         AddTaskJetFlowToyMC                               //
//   Author: Redmer A. Bertens, Utrecht University, 2013, rbertens@cern.ch   //
///////////////////////////////////////////////////////////////////////////////

/* AddTask macro for jet flow toy mc task
 * task uses an afterburner to tune vn in the pico track 
 * selection which can be used by a jet finder 
 * note that this task does not generate MC particles, it changes
 * the azimuthal distribution of already available tracks
*/

class AliAnalysisDataContainer;
class AliAnalysisTaskJetFlowMC;

AliAnalysisTaskJetFlowMC* AddTaskJetFlowMC(
  const char *outputTracks      = "JetFlowToyMC",
  const char *inputTracks       = "PicoTracks",
  const char *name              ="AliAnalysisTaskJetFlowMC",
  Bool_t doQA                   = kFALSE,
  Bool_t doDecay                = kFALSE
  )
{  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)                             return 0x0;
  if (!mgr->GetInputEventHandler())     return 0x0;
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":";
  fileName += name;
        // create the task
  AliAnalysisTaskJetFlowMC *eTask = new AliAnalysisTaskJetFlowMC(name, doQA);
  eTask->SetTracksOutName(outputTracks);
  eTask->SetTracksInName(inputTracks);
        // connect input and output
  mgr->AddTask(eTask);
  mgr->ConnectInput  (eTask, 0, mgr->GetCommonInputContainer());
  if(doQA) {
      // this task only produces output when the qa flag is set to true
      mgr->ConnectOutput (eTask, 1, mgr->CreateContainer(Form("%s_container", fileName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  }
  if(doDecay) task->SetDecayer(TaskJetFlowMC::GetDecayer(kTRUE));
  return eTask;
}

namespace TaskJetFlowMC {

    TF1* GetSpectrum() {
        // thermal spectrum used for ALICE SIMULATION PLOTS ALI-SIMUL-75145 ALI-SIMUL-75171
        TF1* spectrum = new TF1("fspectrum", "[0]*(TMath::Power([1], 2)*x*TMath::Exp(-[1]*x))+(x>1)*[2]*(1.17676e-01*TMath::Sqrt(0.1396*0.1396+x*x)*TMath::Power(1.+1./[3]/8.21795e-01*TMath::Sqrt(0.1396*0.1396+x*x),-1.*[3]))*(1/(1 + TMath::Exp(-(x - [4])/[5])))", .2, 200.);
        fspectrum->SetParameters(2434401791.20528 ,2.98507 ,10069622.25117 ,5.50000 ,2.80000 ,0.20000 );
        return fspectrum;   
    }

    TVirtualDecayer* GetDecayer(Bool_t local = kTRUE) {
        // setup a decayer
        if(local) gSystem->Load("$ALICE_ROOT/lib/tgt_linuxx8664gcc/libpythia6");
        TPythia6Decayer* decayer = new TPythia6Decayer();
        decayer->SetForceDecay(TPythia6Decayer::kHardonicD);
        return decayer;
    }

}
