// $Id$
//
// Track/cluster matcher.
//
//

#include <TClonesArray.h>
#include <TString.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliPicoTrack.h"
#include "AliLog.h"

#include "AliEmcalClusTrackMatcherTask.h"

ClassImp(AliEmcalClusTrackMatcherTask)

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::AliEmcalClusTrackMatcherTask(const char *name) : 
  AliAnalysisTaskSE("AliEmcalClusTrackMatcherTask"),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fDoClusTrack(1),
  fDoTrackClus(0)
{
  // Standard constructor.

  if (!name)
    return;
  SetName(name);
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,CaloClusters,Tracks";
}

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::~AliEmcalClusTrackMatcherTask()
{
  // Destructor.
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

  TList *l = InputEvent()->GetList();
  if (!l) 
    return;

  TClonesArray *tracks = dynamic_cast<TClonesArray*>(l->FindObject(fTracksName));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
    return;
  }

  TClonesArray *clus = dynamic_cast<TClonesArray*>(l->FindObject(fCaloName));
  if (!clus) {
    AliError(Form("Pointer to clus %s == 0", fCaloName.Data() ));
    return;
  }

  const Int_t Ntrks = tracks->GetEntries();
  const Int_t Ncls  = clus->GetEntries();

  if (fDoClusTrack) {
    for(Int_t i=0; i < Ncls; ++i) {
      AliVCluster *c = dynamic_cast<AliVCluster*>(clus->At(i));
      if (!c)
        continue;
      c->SetEmcCpvDistance(-1);
      c->SetTrackDistance(999,999);
      Double_t dEtaMin  = 1e9;
      Double_t dPhiMin  = 1e9;
      Double_t dRMin    = 1e9;
      Int_t    imin     = -1;
      for(Int_t t = 0; t<Ntrks; ++t) {
        AliVTrack *track = dynamic_cast<AliVTrack*>(tracks->At(t));
      if (!track)
        continue;
      if (!track->IsEMCAL())
          continue;
        Double_t etadiff=999;
        Double_t phidiff=999;
        AliPicoTrack::GetEtaPhiDiff(track,c,phidiff,etadiff);
        Double_t dR = TMath::Sqrt(etadiff*etadiff+phidiff*phidiff);
	if(dR<dRMin) {
          dEtaMin = etadiff;
          dPhiMin = phidiff;
	  dRMin=dR;
          imin = t;
        }
      }
      c->SetEmcCpvDistance(imin);
      c->SetTrackDistance(dPhiMin, dEtaMin);
    }
  }

  if (fDoTrackClus) {
    for(Int_t t = 0; t<Ntrks; ++t) {
      AliVTrack *track = dynamic_cast<AliVTrack*>(tracks->At(t));
      if (!track)
        continue;
      if (!track->IsEMCAL())
        continue;
      Double_t dEtaMin  = 1e9;
      Double_t dPhiMin  = 1e9;
      Double_t dRMin    = 1e9;
      Int_t    imin     = -1;
      for(Int_t i=0; i < Ncls; ++i) {
        AliVCluster *c = dynamic_cast<AliVCluster*>(clus->At(i));
        if (!c)
          continue;
        Double_t etadiff=999;
        Double_t phidiff=999;
        AliPicoTrack::GetEtaPhiDiff(track,c,phidiff,etadiff);
        Double_t dR = TMath::Sqrt(etadiff*etadiff+phidiff*phidiff);
	if(dR<dRMin){
          dEtaMin = etadiff;
          dPhiMin = phidiff;
	  dRMin   = dR;
          imin = i;
        }
      }
      track->SetEMCALcluster(imin);
    }
  }
}

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
