// main17.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates
// (a) how to use UserHooks to regularize onium cross section for pT -> 0,
// (b) how decays could be handled externally.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// A derived class to do J/psi decays.

class JpsiDecay : public DecayHandler {

public:

  // Constructor.
  JpsiDecay(ParticleData* pdtPtrIn, Rndm* rndmPtrIn) {times = 0;
    pdtPtr = pdtPtrIn; rndmPtr = rndmPtrIn;}

  // Routine for doing the decay.
  bool decay(vector<int>& idProd, vector<double>& mProd,
    vector<Vec4>& pProd, int iDec, const Event& event);

private:

  // Count number of times JpsiDecay is called.
  int times;

  // Pointer to the particle data table.
  ParticleData* pdtPtr;

  // Pointer to the random number generator.
  Rndm* rndmPtr;

};

//--------------------------------------------------------------------------

// The actual J/psi decay routine.
// Not intended for realism, just to illustrate the principles.

bool JpsiDecay::decay(vector<int>& idProd, vector<double>& mProd,
  vector<Vec4>& pProd, int iDec, const Event& event) {

  // Always do decay J/psi -> mu+ mu-; store the muons.
  idProd.push_back(-13);
  idProd.push_back(13);

  // Muon mass(es), here from Pythia tables, also stored.
  double mMuon = pdtPtr->m0(13);
  mProd.push_back(mMuon);
  mProd.push_back(mMuon);

  // Calculate muon energy and momentum in J/psi rest frame.
  double eMuon = 0.5 * mProd[0];
  double pAbsMuon = sqrt(eMuon * eMuon - mMuon * mMuon);

  // Assume decay angles isotropic in rest frame.
  double cosTheta = 2. * rndmPtr->flat() - 1.;
  double sinTheta = sqrt(max(0., 1. - cosTheta * cosTheta));
  double phi = 2. * M_PI * rndmPtr->flat();
  double pxMuon = pAbsMuon * sinTheta * cos(phi);
  double pyMuon = pAbsMuon * sinTheta * sin(phi);
  double pzMuon = pAbsMuon * cosTheta;

  // Define mu+ and mu- four-vectors in the J/psi rest frame.
  Vec4 pMuPlus(   pxMuon,  pyMuon,  pzMuon, eMuon);
  Vec4 pMuMinus( -pxMuon, -pyMuon, -pzMuon, eMuon);

  // Boost them by velocity vector of the J/psi mother and store.
  pMuPlus.bst(pProd[0]);
  pMuMinus.bst(pProd[0]);
  pProd.push_back(pMuPlus);
  pProd.push_back(pMuMinus);

  // Print message the first few times, to show that it works.
  if (times++ < 10) {
    int iMother = event[iDec].mother1();
    int idMother = event[iMother].id();
    cout << "\n J/psi decay performed, J/psi in line " << iDec
         << ", mother id = " << idMother << "\n";
  }

  // Done
  return true;

}

//==========================================================================

int main() {

  // Number of events to generate and to list. Max number of errors.
  int nEvent = 2000;
  int nList  = 2;
  int nAbort = 5;

  // Pythia generator.
  Pythia pythia;

  // Initialization for charmonium (singlet+octet) production at the LHC.
  pythia.readString("Charmonium:all = on");
  pythia.readString("Beams:eCM = 7000.");

  // Normally cutoff at pTHat = 1, but push it lower combined with dampening.
  pythia.readString("PhaseSpace:pTHatMin = 0.5");
  pythia.readString("PhaseSpace:pTHatMinDiverge = 0.5");

  // Set up to do a user veto and send it in.
  // First argument: multiplies the pT0 of multiparton interactions
  // to define the pT dampeing scale.
  // Second argument: how many powers of alpha_strong to
  // reweight with new (larger) argument.
  // Third argument: choice of process scale two different ways;
  // probably does not make much difference.
  // See "User Hooks" in manual for detail on SuppressSmallPT.
  UserHooks* oniumUserHook = new SuppressSmallPT( 1., 3, false);
  pythia.setUserHooksPtr( oniumUserHook);

  // A class to do J/psi decays externally.
  DecayHandler* handleDecays = new JpsiDecay(&pythia.particleData,
    &pythia.rndm);

  // The list of particles the class can handle.
  vector<int> handledParticles;
  handledParticles.push_back(443);

  // Hand pointer and list to Pythia.
  pythia.setDecayPtr( handleDecays, handledParticles);

  // Switch off automatic event listing in favour of manual.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Initialization.
  pythia.init();

  // Book histograms.
  Hist pThard("pTHat of hard subprocess", 100, 0., 50.);
  Hist pTJPsi("pT of J/Psi", 100, 0., 50.);

  // Begin event loop.
  int iList = 0;
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Histogram pThard spectrum of process.
    double pTHat = pythia.info.pTHat();
    pThard.fill( pTHat );

    // Look for event with externally handled decays.
    bool externalDecay = false;
    for (int i = 0; i < pythia.event.size(); ++i) {
      int status = pythia.event[i].statusAbs();
      if (status == 93 || status == 94) {externalDecay = true; break;}
    }

    // List first few events with external decay.
    if (externalDecay && ++iList <= nList) {
      pythia.process.list();
      pythia.event.list();
    }

    // Histogram pT spectrum of J/Psi.
   for (int i = 0; i < pythia.event.size(); ++i)
   if (pythia.event[i].id() == 443) pTJPsi.fill( pythia.event[i].pT() );

  // End of event loop.
  }

  // Final statistics. Print histograms.
  pythia.stat();
  cout << pThard << pTJPsi;

  // Done.
  delete handleDecays;
  delete oniumUserHook;
  return 0;
}
