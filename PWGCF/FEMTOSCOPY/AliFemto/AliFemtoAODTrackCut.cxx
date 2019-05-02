///
/// \file AliFemtoAODTrackCut.cxx
///

/**
 * $Log$
 *
 * Revision 1.3  2007/05/22 09:01:42  akisiel
 * Add the possibiloity to save cut settings in the ROOT file
 *
 * Revision 1.2  2007/05/21 10:38:25  akisiel
 * More coding rule conformance
 *
 * Revision 1.1  2007/05/16 10:25:06  akisiel
 * Making the directory structure of AliFemtoUser flat. All files go into one common directory
 *
 * Revision 1.4  2007/05/03 09:46:10  akisiel
 * Fixing Effective C++ warnings
 *
 * Revision 1.3  2007/04/27 07:25:59  akisiel
 * Make revisions needed for compilation from the main AliRoot tree
 *
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.4  2007-04-03 16:00:08  mchojnacki
 * Changes to iprove memory managing
 *
 * Revision 1.3  2007/03/13 15:30:03  mchojnacki
 * adding reader for simulated data
 *
 * Revision 1.2  2007/03/08 14:58:03  mchojnacki
 * adding some alice stuff
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 */

#include "AliFemtoAODTrackCut.h"
#include <cstdio>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoAODTrackCut);
  /// \endcond
#endif



// electron
// 0.13 - 1.8
// 0       7.594129e-02    8.256141e-03
// 1       -5.535827e-01   8.170825e-02
// 2       1.728591e+00    3.104210e-01
// 3       -2.827893e+00   5.827802e-01
// 4       2.503553e+00    5.736207e-01
// 5       -1.125965e+00   2.821170e-01
// 6       2.009036e-01    5.438876e-02

// pion
// 0.13 - 2.0
// 0       1.063457e+00    8.872043e-03
// 1       -4.222208e-01   2.534402e-02
// 2       1.042004e-01    1.503945e-02

// kaon
// 0.18 - 2.0
// 0       -7.289406e-02   1.686074e-03
// 1       4.415666e-01    1.143939e-02
// 2       -2.996790e-01   1.840964e-02
// 3       6.704652e-02    7.783990e-03

// proton
// 0.26 - 2.0
// 0       -3.730200e-02   2.347311e-03
// 1       1.163684e-01    1.319316e-02
// 2       8.354116e-02    1.997948e-02
// 3       -4.608098e-02   8.336400e-03


AliFemtoAODTrackCut::AliFemtoAODTrackCut():
  fCharge(0),
  fLabel(false),
  fMaxchiNdof(1000.0f),
  fMaxSigmaToVertex(1000.0f),
  fNTracksPassed(0),
  fNTracksFailed(0),
  fMostProbable(0)
{
  /// Default constructor
  fPt[0] = 0.0;
  fPt[1] = 100.0;  //100

  fRapidity[0] = -2.0;
  fRapidity[1] = 2.0;

  fPidProbElectron[0] = -1.0;
  fPidProbElectron[1] = 2.0;

  fPidProbPion[0] = -1.0;
  fPidProbPion[1] = 2.0;

  fPidProbKaon[0] = -1.0;
  fPidProbKaon[1] = 2.0;

  fPidProbProton[0] = -1.0;
  fPidProbProton[1] = 2.0;

  fPidProbMuon[0] = -1.0;
  fPidProbMuon[1] = 2.0;
}
//------------------------------
AliFemtoAODTrackCut::~AliFemtoAODTrackCut()
{
  /// noop

}
//------------------------------
bool AliFemtoAODTrackCut::Pass(const AliFemtoTrack *track)
{
  /// test the particle and return
  /// true if it meets all the criteria
  /// false if it doesn't meet at least one of the criteria

  float tMost[5];

  if (((track->ITSchi2() + track->TPCchi2()) / (track->ITSncls() + track->TPCncls())) > fMaxchiNdof) {
    return false;
  }

  if (fMaxSigmaToVertex < track->SigmaToVertex()) {
    return false;
  }

  if (fLabel && track->Label() < 0) {
    fNTracksFailed++;
    return false;
  }
  if (fCharge != 0 && track->Charge() != fCharge) {
    fNTracksFailed++;
    return false;
  }

  float tEnergy = ::sqrt(track->P().Mag2() + fMass * fMass);
  float tRapidity = 0.5 *::log((tEnergy + track->P().z()) / (tEnergy - track->P().z()));
  float tPt = ::sqrt((track->P().x()) * (track->P().x()) + (track->P().y()) * (track->P().y()));

  if ((tRapidity < fRapidity[0]) || (tRapidity > fRapidity[1])) {
    fNTracksFailed++;
    return false;
  }
  if ((tPt < fPt[0]) || (tPt > fPt[1])) {
    fNTracksFailed++;
    return false;
  }


  if ((track->PidProbElectron() < fPidProbElectron[0]) || (track->PidProbElectron() > fPidProbElectron[1])) {
    fNTracksFailed++;
    return false;
  }
  if ((track->PidProbPion() < fPidProbPion[0]) || (track->PidProbPion() > fPidProbPion[1])) {
    fNTracksFailed++;
    //cout<<"No Go Through the cut"<<endl;
    //cout<<fPidProbPion[0]<<" < pi ="<<track->PidProbPion()<<" <"<<fPidProbPion[1]<<endl;
    return false;
  }
  if ((track->PidProbKaon() < fPidProbKaon[0]) || (track->PidProbKaon() > fPidProbKaon[1])) {
    fNTracksFailed++;
    //cout<<"No Go Through the cut"<<endl;
    //cout<<fPidProbKaon[0]<<" < k ="<<track->PidProbKaon()<<" <"<<fPidProbKaon[1]<<endl;
    return false;
  }
  if ((track->PidProbProton() < fPidProbProton[0]) || (track->PidProbProton() > fPidProbProton[1])) {
    fNTracksFailed++;
    //cout<<"No Go Through the cut"<<endl;
    //cout<<fPidProbProton[0]<<" < p  ="<<track->PidProbProton()<<" <"<<fPidProbProton[1]<<endl;
    return false;
  }
  if ((track->PidProbMuon() < fPidProbMuon[0]) || (track->PidProbMuon() > fPidProbMuon[1])) {
    fNTracksFailed++;
    //cout<<"No Go Through the cut"<<endl;
    //cout<<fPidProbMuon[0]<<" <  mi="<<track->PidProbMuon()<<" <"<<fPidProbMuon[1]<<endl;
    return false;
  }

  if (fMostProbable) {
    tMost[0] = track->PidProbElectron() * PidFractionElectron(track->P().Mag());
    tMost[1] = 0.0;
    tMost[2] = track->PidProbPion() * PidFractionPion(track->P().Mag());
    tMost[3] = track->PidProbKaon() * PidFractionKaon(track->P().Mag());
    tMost[4] = track->PidProbProton() * PidFractionProton(track->P().Mag());
    int imost = 0;
    float ipidmax = 0.0;
    for (int ip = 0; ip < 5; ip++)
      if (tMost[ip] > ipidmax) {
        ipidmax = tMost[ip];
        imost = ip;
      };
    if (imost != fMostProbable) return false;
  }

  // cout<<"Go Through the cut"<<endl;
  // cout<<fLabel<<" Label="<<track->Label()<<endl;
  // cout<<fCharge<<" Charge="<<track->Charge()<<endl;
  // cout<<fPt[0]<<" < Pt ="<<Pt<<" <"<<fPt[1]<<endl;
  //cout<<fRapidity[0]<<" < Rapidity ="<<tRapidity<<" <"<<fRapidity[1]<<endl;
  //cout<<fPidProbElectron[0]<<" <  e="<<track->PidProbElectron()<<" <"<<fPidProbElectron[1]<<endl;
  //cout<<fPidProbPion[0]<<" <  pi="<<track->PidProbPion()<<" <"<<fPidProbPion[1]<<endl;
  //cout<<fPidProbKaon[0]<<" <  k="<<track->PidProbKaon()<<" <"<<fPidProbKaon[1]<<endl;
  //cout<<fPidProbProton[0]<<" <  p="<<track->PidProbProton()<<" <"<<fPidProbProton[1]<<endl;
  //cout<<fPidProbMuon[0]<<" <  mi="<<track->PidProbMuon()<<" <"<<fPidProbMuon[1]<<endl;
  fNTracksPassed++;
  return true;
}
//------------------------------
AliFemtoString AliFemtoAODTrackCut::Report()
{
  /// Prepare report from the execution

  AliFemtoString report;
  report += Form("Particle mass:\t%E\n", this->Mass());
  report += Form("Particle charge:\t%d\n", fCharge);
  report += Form("Particle pT:\t%E - %E\n", fPt[0], fPt[1]);
  report += Form("Particle rapidity:\t%E - %E\n", fRapidity[0], fRapidity[1]);
  report += Form("Number of tracks which passed:\t%ld  Number which failed:\t%ld\n", fNTracksPassed, fNTracksFailed);
  return report;
}

TList *AliFemtoAODTrackCut::ListSettings()
{
  /// return a list of settings in a writable form

  TList *tListSetttings = new TList();
  tListSetttings->AddVector(
    new TObjString(Form("AliFemtoAODTrackCut.mass=%f", this->Mass())),
    new TObjString(Form("AliFemtoAODTrackCut.charge=%i", fCharge)),
    new TObjString(Form("AliFemtoAODTrackCut.pidprobpion.minimum=%f", fPidProbPion[0])),
    new TObjString(Form("AliFemtoAODTrackCut.pidprobpion.maximum=%f", fPidProbPion[1])),
    new TObjString(Form("AliFemtoAODTrackCut.pidprobkaon.minimum=%f", fPidProbKaon[0])),
    new TObjString(Form("AliFemtoAODTrackCut.pidprobkaon.maximum=%f", fPidProbKaon[1])),
    new TObjString(Form("AliFemtoAODTrackCut.pidprobproton.minimum=%f", fPidProbProton[0])),
    new TObjString(Form("AliFemtoAODTrackCut.pidprobproton.maximum=%f", fPidProbProton[1])),
    new TObjString(Form("AliFemtoAODTrackCut.pidprobelectron.minimum=%f", fPidProbElectron[0])),
    new TObjString(Form("AliFemtoAODTrackCut.pidprobelectron.maximum=%f", fPidProbElectron[1])),
    new TObjString(Form("AliFemtoAODTrackCut.pidprobMuon.minimum=%f", fPidProbMuon[0])),
    new TObjString(Form("AliFemtoAODTrackCut.pidprobMuon.maximum=%f", fPidProbMuon[1])),
    new TObjString(Form("AliFemtoAODTrackCut.pt.minimum=%f", fPt[0])),
    new TObjString(Form("AliFemtoAODTrackCut.pt.maximum=%f", fPt[1])),
    new TObjString(Form("AliFemtoAODTrackCut.rapidity.minimum=%f", fRapidity[0])),
    new TObjString(Form("AliFemtoAODTrackCut.rapidity.maximum=%f", fRapidity[1])),
    new TObjString(Form("AliFemtoAODTrackCut.maxchindof=%f", fMaxchiNdof)),
    new TObjString(Form("AliFemtoAODTrackCut.maxsigmatovertex=%f", fMaxSigmaToVertex)),
    nullptr);

  if (fMostProbable) {
    const char *type = (fMostProbable == 2) ? "Pion"
                     : (fMostProbable == 3) ? "Kaon"
                     : (fMostProbable == 4) ? "Proton"
                                            : "Unknown";

    tListSetttings->Add(new TObjString(Form("AliFemtoAODTrackCut.mostprobable=%s", type)));
  }

  return tListSetttings;
}

// electron
// 0.13 - 1.8
// 0       7.594129e-02    8.256141e-03
// 1       -5.535827e-01   8.170825e-02
// 2       1.728591e+00    3.104210e-01
// 3       -2.827893e+00   5.827802e-01
// 4       2.503553e+00    5.736207e-01
// 5       -1.125965e+00   2.821170e-01
// 6       2.009036e-01    5.438876e-02
float AliFemtoAODTrackCut::PidFractionElectron(float mom) const
{
  /// Provide a parameterized fraction of electrons dependent on momentum

  if (mom < 0.13 || 1.8 < mom) return 0.0;
  return (7.594129e-02
          - 5.535827e-01 * mom
          + 1.728591e+00 * mom * mom
          - 2.827893e+00 * mom * mom * mom
          + 2.503553e+00 * mom * mom * mom * mom
          - 1.125965e+00 * mom * mom * mom * mom * mom
          + 2.009036e-01 * mom * mom * mom * mom * mom * mom);
}

// pion
// 0.13 - 2.0
// 0       1.063457e+00    8.872043e-03
// 1       -4.222208e-01   2.534402e-02
// 2       1.042004e-01    1.503945e-02
float AliFemtoAODTrackCut::PidFractionPion(float mom) const
{
  /// Provide a parameterized fraction of pions dependent on momentum

  if (mom < 0.13 || 2.0 < mom) return 0.0;
  return (1.063457e+00
          - 4.222208e-01 * mom
          + 1.042004e-01 * mom * mom);
}

// kaon
// 0.18 - 2.0
// 0       -7.289406e-02   1.686074e-03
// 1       4.415666e-01    1.143939e-02
// 2       -2.996790e-01   1.840964e-02
// 3       6.704652e-02    7.783990e-03
float AliFemtoAODTrackCut::PidFractionKaon(float mom) const
{
  /// Provide a parameterized fraction of kaons dependent on momentum

  if (mom < 0.18 || 2.0 < mom) return 0.0;
  return (-7.289406e-02
          + 4.415666e-01 * mom
          - 2.996790e-01 * mom * mom
          + 6.704652e-02 * mom * mom * mom);
}

// proton
// 0.26 - 2.0
// 0       -3.730200e-02   2.347311e-03
// 1       1.163684e-01    1.319316e-02
// 2       8.354116e-02    1.997948e-02
// 3       -4.608098e-02   8.336400e-03
float AliFemtoAODTrackCut::PidFractionProton(float mom) const
{
  /// Provide a parameterized fraction of protons dependent on momentum

  if (mom < 0.26 || 2.0 < mom) return 0.0;
  return (-3.730200e-02
          + 1.163684e-01 * mom
          + 8.354116e-02 * mom * mom
          - 4.608098e-02 * mom * mom * mom);
}
