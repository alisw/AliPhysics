///
/// \file AliFemtoBasicTrackCut.cxx
///

#include "AliFemtoBasicTrackCut.h"
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoBasicTrackCut);
  /// \endcond
#endif

AliFemtoBasicTrackCut::AliFemtoBasicTrackCut():
  fCharge(0),
  fNTracksPassed(0),
  fNTracksFailed(0)
{
  // Default constructor
  fCharge = 1;  // takes both charges 0

  fNSigmaPion[0] = -100.0;
  fNSigmaPion[1] = 100.0;

  fNSigmaKaon[0] = -100.0;
  fNSigmaKaon[1] = 100.0;

  fNSigmaProton[0] = -100.0;
  fNSigmaProton[1] = 100.0;

  fNHits[0] = 10;
  fNHits[1] = 180;

  fPt[0] = 0.0;
  fPt[1] = 100.0;  // 100

  fRapidity[0] = -2;
  fRapidity[1] = 2;  //-2 2

  fDCA[0] = -1.0;
  fDCA[1] = 20.0;

}
//------------------------------
bool AliFemtoBasicTrackCut::Pass(const AliFemtoTrack* track)
{
  // test the particle and return true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria

  bool goodPID = true;

  /* ----- NOT DOING PID CUTS !!!! ------
  bool goodPID = ((track->NSigmaPion()   > fNSigmaPion[0]) &&
                  (track->NSigmaPion()   < fNSigmaPion[1]) &&
                  (track->NSigmaKaon()   > fNSigmaKaon[0]) &&
                  (track->NSigmaKaon()   < fNSigmaKaon[1]) &&
                  (track->NSigmaProton() > fNSigmaProton[0]) &&
                  (track->NSigmaProton() < fNSigmaProton[1]));
  ----- NOT DOING PID CUTS !!!! ------ */

  // if user requests "charge=0" then that means ignore charge
  if (fCharge != 0) {
    goodPID = (goodPID && (track->Charge() == fCharge));
  }

  if (!goodPID) {
    fNTracksFailed++;
    return goodPID;
  }

  const float tEnergy = ::sqrt(track->P().Mag2() + fMass * fMass),
                mom_z = track->P().z(),
            tRapidity = 0.5*::log((tEnergy + mom_z) / (tEnergy - mom_z)),
                   pt = track->Pt();

  // const int hits = track->NHits();
  // const float dca = track->DCAxy();

  bool goodTrack = (fPt[0] <= pt && pt < fPt[1])
                && (fRapidity[0] <= tRapidity && tRapidity <= fRapidity[1]);
                // && (fDCA[0] <= dca && dca <= fDCA[1])
                // && (fNHits[0] <= hits && hits < fNHits[1]);

  //       (track->PidProbPion()>0.5) && //moje
  //       (track->PidProbMuon()<0.47) && //moje
  //       (track->Label()>0); //moje

  goodTrack ? fNTracksPassed++ : fNTracksFailed++;
  return goodTrack;
}
//------------------------------
AliFemtoString AliFemtoBasicTrackCut::Report()
{
  // construct report
  TString report;

  report += TString::Format("Particle mass:\t%E\n", Mass());
  report += TString::Format("Particle charge:\t%d\n", fCharge);
  report += TString::Format("Particle Nsigma from pion:\t%E - %E\n", fNSigmaPion[0], fNSigmaPion[1]);
  report += TString::Format("Particle Nsigma from kaon:\t%E - %E\n", fNSigmaKaon[0], fNSigmaKaon[1]);
  report += TString::Format("Particle Nsigma from proton:\t%E - %E\n", fNSigmaProton[0], fNSigmaProton[1]);
  report += TString::Format("Particle #hits:\t%d - %d\n", fNHits[0], fNHits[1]);
  report += TString::Format("Particle pT:\t%E - %E\n", fPt[0], fPt[1]);
  report += TString::Format("Particle rapidity:\t%E - %E\n", fRapidity[0], fRapidity[1]);
  report += TString::Format("Particle DCA:\t%E - %E\n", fDCA[0], fDCA[1]);
  report += TString::Format("Number of tracks which passed:\t%ld  Number which failed:\t%ld\n", fNTracksPassed, fNTracksFailed);

  return AliFemtoString((const char *)report);
}

TList *AliFemtoBasicTrackCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *settings_list = new TList();

  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.mass=%f", Mass())
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.charge=%i", fCharge)
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.nsigmapion.minimum=%f", fNSigmaPion[0])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.nsigmapion.maximum=%f", fNSigmaPion[1])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.nsigmakaon.minimum=%f", fNSigmaKaon[0])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.nsigmakaon.maximum=%f", fNSigmaKaon[1])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.nsigmaproton.minimum=%f", fNSigmaProton[0])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.nsigmaproton.maximum=%f", fNSigmaProton[1])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.nhits.minimum=%i", fNHits[0])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.nhits.maximum=%i", fNHits[1])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.pt.minimum=%f", fPt[0])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.pt.maximum=%f", fPt[1])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.rapidity.minimum=%f", fRapidity[0])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.rapidity.maximum=%f", fRapidity[1])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.dca.minimum=%f", fDCA[0])
  ));
  settings_list->Add(new TObjString(
    TString::Format("AliFemtoBasicTrackCut.dca.maximum=%f", fDCA[1])
  ));

  return settings_list;
}
