// $Id$

#include "AliEmcalClusTrackMatcherTask.h"
//#include <TClonesArray.h>
//#include <TParticle.h>
//#include "AliAODJet.h"
//#include "AliAnalysisManager.h"
//#include "AliESDtrack.h"
//#include "AliESDtrackCuts.h"
//#include "AliBambooFJWrapper.h"
//#include "AliESDCaloCluster.h"

ClassImp(AliEmcalClusTrackMatcherTask)


//________________________________________________________________________
AliEmcalClusTrackMatcherTask::AliEmcalClusTrackMatcherTask(const char *name) : 
  AliAnalysisTaskSE("AliEmcalClusTrackMatcherTask"),
  fTracksName("Tracks"),
  fCaloName("CaloClusters")
{
  // Standard constructor.
  if (!name)
    return;
  SetName(name);
//  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,CaloClusters,Tracks";
}

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::~AliEmcalClusTrackMatcherTask()
{
  // Destructor
}

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::UserCreateOutputObjects()
{
  // Create user objects.

}

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.

#if 0
  // add jets to event if not yet there
  if (!(InputEvent()->FindListObject(fJetsName)))
    InputEvent()->AddObject(fJets);

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  TClonesArray *tracks = 0;
  TClonesArray *clus   = 0;
  TList *l = InputEvent()->GetList();
  if ((fType==0)||(fType==1)) {
    if (fPrimTracksName == "Tracks")
      am->LoadBranch("Tracks");
    tracks = dynamic_cast<TClonesArray*>(l->FindObject(fPrimTracksName));
    if (!tracks) {
      AliError(Form("Pointer to tracks %s == 0", fPrimTracksName.Data() ));
      return;
    }
  }
  if ((fType==0)||(fType==2)) {
    if (fCaloName == "CaloClusters")
      am->LoadBranch("CaloClusters");
    clus = dynamic_cast<TClonesArray*>(l->FindObject(fCaloName));
    if (!clus) {
      AliError(Form("Pointer to clus %s == 0", fCaloName.Data() ));
      return;
    }
    if (!(InputEvent()->FindListObject(Form("%s_neutrals",fJetsName.Data()))))
      InputEvent()->AddObject(fNeutrals);
  }
#endif
}

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.

}
