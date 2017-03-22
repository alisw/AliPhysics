// ProcessLevel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the ProcessLevel class.

#include "Pythia8/ProcessLevel.h"

namespace Pythia8 {

//==========================================================================

// The ProcessLevel class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Allow a few failures in final construction of events.
const int ProcessLevel::MAXLOOP = 5;

//--------------------------------------------------------------------------

// Destructor.

ProcessLevel::~ProcessLevel() {

  // Run through list of first hard processes and delete them.
  for (int i = 0; i < int(containerPtrs.size()); ++i)
    delete containerPtrs[i];

  // Run through list of second hard processes and delete them.
  for (int i = 0; i < int(container2Ptrs.size()); ++i)
    delete container2Ptrs[i];

}

//--------------------------------------------------------------------------

// Main routine to initialize generation process.

bool ProcessLevel::init( Info* infoPtrIn, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
  BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  Couplings* couplingsPtrIn, SigmaTotal* sigmaTotPtrIn, bool doLHA,
  SLHAinterface* slhaInterfacePtrIn, UserHooks* userHooksPtrIn,
  vector<SigmaProcess*>& sigmaPtrs, vector<PhaseSpace*>& phaseSpacePtrs,
  ostream& os) {

  // Store input pointers for future use.
  infoPtr          = infoPtrIn;
  particleDataPtr  = particleDataPtrIn;
  rndmPtr          = rndmPtrIn;
  beamAPtr         = beamAPtrIn;
  beamBPtr         = beamBPtrIn;
  couplingsPtr     = couplingsPtrIn;
  sigmaTotPtr      = sigmaTotPtrIn;
  userHooksPtr     = userHooksPtrIn;
  slhaInterfacePtr = slhaInterfacePtrIn;

  // Send on some input pointers.
  resonanceDecays.init( infoPtr, particleDataPtr, rndmPtr);

  // Set up SigmaTotal. Store sigma_nondiffractive for future use.
  sigmaTotPtr->init( infoPtr, settings, particleDataPtr);
  int    idA = infoPtr->idA();
  int    idB = infoPtr->idB();
  double eCM = infoPtr->eCM();
  sigmaTotPtr->calc( idA, idB, eCM);
  sigmaND = sigmaTotPtr->sigmaND();

  // Options to allow second hard interaction and resonance decays.
  doSecondHard  = settings.flag("SecondHard:generate");
  doSameCuts    = settings.flag("PhaseSpace:sameForSecond");
  doResDecays   = settings.flag("ProcessLevel:resonanceDecays");
  startColTag   = settings.mode("Event:startColTag");

  // Second interaction not to be combined with biased phase space.
  if (doSecondHard && userHooksPtr != 0
  && userHooksPtr->canBiasSelection()) {
    infoPtr->errorMsg("Error in ProcessLevel::init: "
      "cannot combine second interaction with biased phase space");
    return false;
  }

  // Mass and pT cuts for two hard processes.
  mHatMin1      = settings.parm("PhaseSpace:mHatMin");
  mHatMax1      = settings.parm("PhaseSpace:mHatMax");
  if (mHatMax1 < mHatMin1) mHatMax1 = eCM;
  pTHatMin1     = settings.parm("PhaseSpace:pTHatMin");
  pTHatMax1     = settings.parm("PhaseSpace:pTHatMax");
  if (pTHatMax1 < pTHatMin1) pTHatMax1 = eCM;
  mHatMin2      = settings.parm("PhaseSpace:mHatMinSecond");
  mHatMax2      = settings.parm("PhaseSpace:mHatMaxSecond");
  if (mHatMax2 < mHatMin2) mHatMax2 = eCM;
  pTHatMin2     = settings.parm("PhaseSpace:pTHatMinSecond");
  pTHatMax2     = settings.parm("PhaseSpace:pTHatMaxSecond");
  if (pTHatMax2 < pTHatMin2) pTHatMax2 = eCM;

  // Check whether mass and pT ranges agree or overlap.
  cutsAgree     = doSameCuts;
  if (mHatMin2 == mHatMin1 && mHatMax2 == mHatMax1 && pTHatMin2 == pTHatMin1
      && pTHatMax2 == pTHatMax1) cutsAgree = true;
  cutsOverlap   = cutsAgree;
  if (!cutsAgree) {
    bool mHatOverlap = (max( mHatMin1, mHatMin2)
                      < min( mHatMax1, mHatMax2));
    bool pTHatOverlap = (max( pTHatMin1, pTHatMin2)
                       < min( pTHatMax1, pTHatMax2));
    if (mHatOverlap && pTHatOverlap) cutsOverlap = true;
  }

  // Set up containers for all the internal hard processes.
  SetupContainers setupContainers;
  setupContainers.init(containerPtrs, infoPtr, settings, particleDataPtr,
                       couplingsPtr);

  // Append containers for external hard processes, if any.
  if (sigmaPtrs.size() > 0) {
    for (int iSig = 0; iSig < int(sigmaPtrs.size()); ++iSig)
      containerPtrs.push_back( new ProcessContainer(sigmaPtrs[iSig],
        true, phaseSpacePtrs[iSig]) );
  }

  // Append single container for Les Houches processes, if any.
  if (doLHA) {
    SigmaProcess* sigmaPtr = new SigmaLHAProcess();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );

    // Store location of this container, and send in LHA pointer.
    iLHACont = containerPtrs.size() - 1;
    containerPtrs[iLHACont]->setLHAPtr(lhaUpPtr);
  }

  // If no processes found then refuse to do anything.
  if ( int(containerPtrs.size()) == 0) {
    infoPtr->errorMsg("Error in ProcessLevel::init: "
      "no process switched on");
    return false;
  }

  // Check whether pT-based weighting in 2 -> 2 is requested.
  if (settings.flag("PhaseSpace:bias2Selection")) {
    bool bias2Sel = false;
    if (sigmaPtrs.size() == 0 && !doLHA && !doSecondHard) {
      bias2Sel = true;
      for (int i = 0; i < int(containerPtrs.size()); ++i) {
        if (containerPtrs[i]->nFinal() != 2) bias2Sel = false;
        int code = containerPtrs[i]->code();
        if (code > 100 && code < 110) bias2Sel = false;
      }
    }
    if (!bias2Sel) {
      infoPtr->errorMsg("Error in ProcessLevel::init: "
        "requested event weighting not possible");
      return false;
    }
  }

  // Check that SUSY couplings were indeed initialized where necessary.
  bool hasSUSY = false;
  for (int i = 0; i < int(containerPtrs.size()); ++i)
    if (containerPtrs[i]->isSUSY()) hasSUSY = true;

  // If SUSY processes requested but no SUSY couplings present
  if(hasSUSY && !couplingsPtr->isSUSY) {
    infoPtr->errorMsg("Error in ProcessLevel::init: "
      "SUSY process switched on but no SUSY couplings found");
    return false;
  }

  // Fill SLHA blocks SMINPUTS and MASS from PYTHIA SM parameter values.
  slhaInterfacePtr->pythia2slha(particleDataPtr);

  // Initialize each process.
  int numberOn = 0;
  for (int i = 0; i < int(containerPtrs.size()); ++i)
    if (containerPtrs[i]->init(true, infoPtr, settings, particleDataPtr,
      rndmPtr, beamAPtr, beamBPtr, couplingsPtr, sigmaTotPtr,
      &resonanceDecays, slhaInterfacePtr, userHooksPtr)) ++numberOn;

  // Sum maxima for Monte Carlo choice.
  sigmaMaxSum = 0.;
  for (int i = 0; i < int(containerPtrs.size()); ++i)
    sigmaMaxSum += containerPtrs[i]->sigmaMax();

  // Option to pick a second hard interaction: repeat as above.
  int number2On = 0;
  if (doSecondHard) {
    setupContainers.init2(container2Ptrs, settings);
    if ( int(container2Ptrs.size()) == 0) {
      infoPtr->errorMsg("Error in ProcessLevel::init: "
        "no second hard process switched on");
      return false;
    }
    for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
      if (container2Ptrs[i2]->init(false, infoPtr, settings, particleDataPtr,
        rndmPtr, beamAPtr, beamBPtr, couplingsPtr, sigmaTotPtr,
        &resonanceDecays, slhaInterfacePtr, userHooksPtr)) ++number2On;
    sigma2MaxSum = 0.;
    for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
      sigma2MaxSum += container2Ptrs[i2]->sigmaMax();
  }

  // Printout during initialization is optional.
  if (settings.flag("Init:showProcesses")) {

    // Construct string with incoming beams and for cm energy.
    string collision = "We collide " + particleDataPtr->name(idA)
      + " with " + particleDataPtr->name(idB) + " at a CM energy of ";
    string pad( 51 - collision.length(), ' ');

    // Print initialization information: header.
    os << "\n *-------  PYTHIA Process Initialization  ---------"
       << "-----------------*\n"
       << " |                                                   "
       << "               |\n"
       << " | " << collision << scientific << setprecision(3)
       << setw(9) << eCM << " GeV" << pad << " |\n"
       << " |                                                   "
       << "               |\n"
       << " |---------------------------------------------------"
       << "---------------|\n"
       << " |                                                   "
       << " |             |\n"
       << " | Subprocess                                    Code"
       << " |   Estimated |\n"
       << " |                                                   "
       << " |    max (mb) |\n"
       << " |                                                   "
       << " |             |\n"
       << " |---------------------------------------------------"
       << "---------------|\n"
       << " |                                                   "
       << " |             |\n";

    // Loop over existing processes: print individual process info.
    map<int, double> sigmaMaxM;
    map<int, string> nameM;
    for (int i = 0; i < int(containerPtrs.size()); ++i) {
      int code = containerPtrs[i]->code();
      nameM[code] = containerPtrs[i]->name();
      sigmaMaxM[code] = containerPtrs[i]->sigmaMax() > sigmaMaxM[code] ?
        containerPtrs[i]->sigmaMax() : sigmaMaxM[code];
    }
    for (map<int, string>::iterator i = nameM.begin(); i != nameM.end(); ++i)
      os << " | " << left << setw(45) << i->second
         << right << setw(5) << i->first << " | "
         << scientific << setprecision(3) << setw(11)
         << sigmaMaxM[i->first] << " |\n";

    // Loop over second hard processes, if any, and repeat as above.
    if (doSecondHard) {
      os << " |                                                   "
         << " |             |\n"
         << " |---------------------------------------------------"
         <<"---------------|\n"
         << " |                                                   "
         <<" |             |\n";
      sigmaMaxM.clear();
      nameM.clear();
      for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2) {
        int code = container2Ptrs[i2]->code();
        nameM[code] = container2Ptrs[i2]->name();
        sigmaMaxM[code] = container2Ptrs[i2]->sigmaMax() > sigmaMaxM[code] ?
          container2Ptrs[i2]->sigmaMax() : sigmaMaxM[code];
      }
      for (map<int, string>::iterator i2 = nameM.begin(); i2 != nameM.end();
           ++i2)
        os << " | " << left << setw(45) << i2->second
           << right << setw(5) << i2->first << " | "
           << scientific << setprecision(3) << setw(11)
           << sigmaMaxM[i2->first] << " |\n";
    }

    // Listing finished.
    os << " |                                                     "
       << "             |\n"
       << " *-------  End PYTHIA Process Initialization ----------"
       <<"-------------*" << endl;
  }

  // If sum of maxima vanishes then refuse to do anything.
  if ( numberOn == 0  || sigmaMaxSum <= 0.) {
    infoPtr->errorMsg("Error in ProcessLevel::init: "
      "all processes have vanishing cross sections");
    return false;
  }
  if ( doSecondHard && (number2On == 0  || sigma2MaxSum <= 0.) ) {
    infoPtr->errorMsg("Error in ProcessLevel::init: "
      "all second hard processes have vanishing cross sections");
    return false;
  }

  // If two hard processes then check whether some (but not all) agree.
  allHardSame  = true;
  noneHardSame = true;
  if (doSecondHard) {
    bool foundMatch = false;

    // Check for each first process if matched in second.
    for (int i = 0; i < int(containerPtrs.size()); ++i) {
      foundMatch = false;
      if (cutsOverlap)
      for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
        if (container2Ptrs[i2]->code() == containerPtrs[i]->code())
          foundMatch = true;
      containerPtrs[i]->isSame( foundMatch );
      if (!foundMatch)  allHardSame = false;
      if ( foundMatch) noneHardSame = false;
    }

    // Check for each second process if matched in first.
    for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2) {
      foundMatch = false;
      if (cutsOverlap)
      for (int i = 0; i < int(containerPtrs.size()); ++i)
        if (containerPtrs[i]->code() == container2Ptrs[i2]->code())
          foundMatch = true;
      container2Ptrs[i2]->isSame( foundMatch );
      if (!foundMatch)  allHardSame = false;
      if ( foundMatch) noneHardSame = false;
    }
  }

  // Concluding classification, including cuts.
  if (!cutsAgree) allHardSame = false;
  someHardSame = !allHardSame && !noneHardSame;

  // Reset counters for average impact-parameter enhancement.
  nImpact       = 0;
  sumImpactFac  = 0.;
  sum2ImpactFac = 0.;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Main routine to generate the hard process.

bool ProcessLevel::next( Event& process) {

  // Generate the next event with two or one hard interactions.
  bool physical = (doSecondHard) ? nextTwo( process) : nextOne( process);

  // Check that colour assignments make sense.
  if (physical) physical = checkColours( process);

  // Done.
  return physical;
}

//--------------------------------------------------------------------------

// Generate (= read in) LHA input of resonance decay only.

bool ProcessLevel::nextLHAdec( Event& process) {

  // Read resonance decays from LHA interface.
  infoPtr->setEndOfFile(false);
  if (!lhaUpPtr->setEvent()) {
    infoPtr->setEndOfFile(true);
    return false;
  }

  // Store LHA output in standard event record format.
  containerLHAdec.constructDecays( process);

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Accumulate and update statistics (after possible user veto).

void ProcessLevel::accumulate() {

  // Increase number of accepted events.
  containerPtrs[iContainer]->accumulate();

  // Provide current generated cross section estimate.
  long   nTrySum    = 0;
  long   nSelSum    = 0;
  long   nAccSum    = 0;
  double sigmaSum   = 0.;
  double delta2Sum  = 0.;
  double sigSelSum  = 0.;
  double weightSum  = 0.;
  int    codeNow;
  long   nTryNow, nSelNow, nAccNow;
  double sigmaNow, deltaNow, sigSelNow, weightNow;
  map<int, bool> duplicate;
  for (int i = 0; i < int(containerPtrs.size()); ++i)
  if (containerPtrs[i]->sigmaMax() != 0.) {
    codeNow         = containerPtrs[i]->code();
    nTryNow         = containerPtrs[i]->nTried();
    nSelNow         = containerPtrs[i]->nSelected();
    nAccNow         = containerPtrs[i]->nAccepted();
    sigmaNow        = containerPtrs[i]->sigmaMC();
    deltaNow        = containerPtrs[i]->deltaMC();
    sigSelNow       = containerPtrs[i]->sigmaSelMC();
    weightNow       = containerPtrs[i]->weightSum();
    nTrySum        += nTryNow;
    nSelSum        += nSelNow;
    nAccSum        += nAccNow;
    sigmaSum       += sigmaNow;
    delta2Sum      += pow2(deltaNow);
    sigSelSum      += sigSelNow;
    weightSum      += weightNow;
    if (!doSecondHard) {
      if (!duplicate[codeNow])
        infoPtr->setSigma( codeNow, containerPtrs[i]->name(),
          nTryNow, nSelNow, nAccNow, sigmaNow, deltaNow, weightNow);
      else
        infoPtr->addSigma( codeNow, nTryNow, nSelNow, nAccNow, sigmaNow,
          deltaNow);
      duplicate[codeNow] = true;
    }
  }

  // Normally only one hard interaction. Then store info and done.
  if (!doSecondHard) {
    double deltaSum = sqrtpos(delta2Sum);
    infoPtr->setSigma( 0, "sum", nTrySum, nSelSum, nAccSum, sigmaSum, deltaSum,
      weightSum);
    return;
  }

  // Increase counter for a second hard interaction.
  container2Ptrs[i2Container]->accumulate();

  // Update statistics on average impact factor.
  ++nImpact;
  sumImpactFac     += infoPtr->enhanceMPI();
  sum2ImpactFac    += pow2(infoPtr->enhanceMPI());

  // Cross section estimate for second hard process.
  double sigma2Sum  = 0.;
  double sig2SelSum = 0.;
  for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
  if (container2Ptrs[i2]->sigmaMax() != 0.) {
    nTrySum        += container2Ptrs[i2]->nTried();
    sigma2Sum      += container2Ptrs[i2]->sigmaMC();
    sig2SelSum     += container2Ptrs[i2]->sigmaSelMC();
  }

  // Average impact-parameter factor and error.
  double invN       = 1. / max(1, nImpact);
  double impactFac  = max( 1., sumImpactFac * invN);
  double impactErr2 = ( sum2ImpactFac * invN / pow2(impactFac) - 1.) * invN;

  // Cross section estimate for combination of first and second process.
  // Combine two possible ways and take average.
  double sigmaComb  = 0.5 * (sigmaSum * sig2SelSum + sigSelSum * sigma2Sum);
  sigmaComb        *= impactFac / sigmaND;
  if (allHardSame) sigmaComb *= 0.5;
  double deltaComb  = sqrtpos(2. / nAccSum + impactErr2) * sigmaComb;

  // Store info and done.
  infoPtr->setSigma( 0, "sum", nTrySum, nSelSum, nAccSum, sigmaComb, deltaComb,
    weightSum);

}

//--------------------------------------------------------------------------

// Print statistics on cross sections and number of events.

void ProcessLevel::statistics(bool reset, ostream& os) {

  // Special processing if two hard interactions selected.
  if (doSecondHard) {
    statistics2(reset, os);
    return;
  }

  // Header.
  os << "\n *-------  PYTHIA Event and Cross Section Statistics  ------"
     << "-------------------------------------------------------*\n"
     << " |                                                            "
     << "                                                     |\n"
     << " | Subprocess                                    Code |       "
     << "     Number of events       |      sigma +- delta    |\n"
     << " |                                                    |       "
     << "Tried   Selected   Accepted |     (estimated) (mb)   |\n"
     << " |                                                    |       "
     << "                            |                        |\n"
     << " |------------------------------------------------------------"
     << "-----------------------------------------------------|\n"
     << " |                                                    |       "
     << "                            |                        |\n";

  // Reset sum counters.
  long   nTrySum   = 0;
  long   nSelSum   = 0;
  long   nAccSum   = 0;
  double sigmaSum  = 0.;
  double delta2Sum = 0.;

  // Reset process maps.
  map<int, string> nameM;
  map<int, long> nTryM, nSelM, nAccM;
  map<int, double> sigmaM, delta2M;
  vector<ProcessContainer*> lheContainerPtrs;

  // Loop over existing processes.
  for (int i = 0; i < int(containerPtrs.size()); ++i)
  if (containerPtrs[i]->sigmaMax() != 0.) {

    // Read info for process. Sum counters.
    nTrySum       += containerPtrs[i]->nTried();
    nSelSum       += containerPtrs[i]->nSelected();
    nAccSum       += containerPtrs[i]->nAccepted();
    sigmaSum      += containerPtrs[i]->sigmaMC();
    delta2Sum     += pow2(containerPtrs[i]->deltaMC());

    // Skip Les Houches containers.
    if (containerPtrs[i]->code() == 9999) {
      lheContainerPtrs.push_back(containerPtrs[i]);
      continue;
    }

    // Internal process info.
    int code = containerPtrs[i]->code();
    nameM[code]   = containerPtrs[i]->name();
    nTryM[code]  += containerPtrs[i]->nTried();
    nSelM[code]  += containerPtrs[i]->nSelected();
    nAccM[code]  += containerPtrs[i]->nAccepted();
    sigmaM[code] += containerPtrs[i]->sigmaMC();
    delta2M[code]+= pow2(containerPtrs[i]->deltaMC());
  }

  // Print internal process info.
  for (map<int, string>::iterator i = nameM.begin(); i != nameM.end(); ++i) {
    int code = i->first;
    os << " | " << left << setw(45) << i->second
       << right << setw(5) << code << " | "
       << setw(11) << nTryM[code] << " " << setw(10) << nSelM[code] << " "
       << setw(10) << nAccM[code] << " | " << scientific << setprecision(3)
       << setw(11) << sigmaM[code]
       << setw(11) << sqrtpos(delta2M[code]) << " |\n";
  }

  // Print Les Houches process info.
  for (int i = 0; i < int(lheContainerPtrs.size()); ++i) {
    ProcessContainer *ptr = lheContainerPtrs[i];
    os << " | " << left << setw(45) << ptr->name()
       << right << setw(5) << ptr->code() << " | "
       << setw(11) << ptr->nTried() << " " << setw(10) << ptr->nSelected()
       << " " << setw(10) << ptr->nAccepted() << " | " << scientific
       << setprecision(3) << setw(11) << ptr->sigmaMC() << setw(11)
       << ptr->deltaMC() << " |\n";

    // Print subdivision by user code for Les Houches process.
    for (int j = 0; j < ptr->codeLHASize(); ++j)
      os << " |    ... whereof user classification code " << setw(10)
         << ptr->subCodeLHA(j) << " | " << setw(11) << ptr->nTriedLHA(j)
         << " " << setw(10) << ptr->nSelectedLHA(j) << " " << setw(10)
         << ptr->nAcceptedLHA(j) << " |                        | \n";
  }

  // Print summed process info.
  os << " |                                                    |       "
     << "                            |                        |\n"
     << " | " << left << setw(50) << "sum" << right << " | " << setw(11)
     << nTrySum << " " << setw(10) << nSelSum << " " << setw(10)
     << nAccSum << " | " << scientific << setprecision(3) << setw(11)
     << sigmaSum << setw(11) << sqrtpos(delta2Sum) << " |\n";

  // Listing finished.
  os << " |                                                            "
     << "                                                     |\n"
     << " *-------  End PYTHIA Event and Cross Section Statistics -----"
     << "-----------------------------------------------------*" << endl;

  // Optionally reset statistics contants.
  if (reset) resetStatistics();

}

//--------------------------------------------------------------------------

// Reset statistics on cross sections and number of events.

void ProcessLevel::resetStatistics() {

  for (int i = 0; i < int(containerPtrs.size()); ++i)
    containerPtrs[i]->reset();
  if (doSecondHard)
  for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
    container2Ptrs[i2]->reset();

}

//--------------------------------------------------------------------------

// Generate the next event with one interaction.

bool ProcessLevel::nextOne( Event& process) {

  // Update CM energy for phase space selection.
  double eCM = infoPtr->eCM();
  for (int i = 0; i < int(containerPtrs.size()); ++i)
    containerPtrs[i]->newECM(eCM);

  // Outer loop in case of rare failures.
  bool physical = true;
  for (int loop = 0; loop < MAXLOOP; ++loop) {
    if (!physical) process.clear();
    physical = true;

    // Loop over tries until trial event succeeds.
    for ( ; ; ) {

      // Pick one of the subprocesses.
      double sigmaMaxNow = sigmaMaxSum * rndmPtr->flat();
      int iMax = containerPtrs.size() - 1;
      iContainer = -1;
      do sigmaMaxNow -= containerPtrs[++iContainer]->sigmaMax();
      while (sigmaMaxNow > 0. && iContainer < iMax);

      // Do a trial event of this subprocess; accept or not.
      if (containerPtrs[iContainer]->trialProcess()) break;

      // Check for end-of-file condition for Les Houches events.
      if (infoPtr->atEndOfFile()) return false;
    }

    // Update sum of maxima if current maximum violated.
    if (containerPtrs[iContainer]->newSigmaMax()) {
      sigmaMaxSum = 0.;
      for (int i = 0; i < int(containerPtrs.size()); ++i)
        sigmaMaxSum += containerPtrs[i]->sigmaMax();
    }

    // Construct kinematics of acceptable process.
    containerPtrs[iContainer]->constructState();
    if ( !containerPtrs[iContainer]->constructProcess( process) )
      physical = false;

    // Do all resonance decays.
    if ( physical && doResDecays
      && !containerPtrs[iContainer]->decayResonances( process) )
      physical = false;

    // Retry process for unphysical states.
    for (int i =1; i < process.size(); ++i)
      if (process[i].e() < 0.) {
        infoPtr->errorMsg("Error in ProcessLevel::nextOne: "
          "Constructed particle with negative energy.");
        physical = false;
      }

    // Add any junctions to the process event record list.
    if (physical) findJunctions( process);

    // Outer loop should normally work first time around.
    if (physical) break;
  }

  // Done.
  return physical;
}

//--------------------------------------------------------------------------

// Generate the next event with two hard interactions.

bool ProcessLevel::nextTwo( Event& process) {

  // Update CM energy for phase space selection.
  double eCM = infoPtr->eCM();
  for (int i = 0; i < int(containerPtrs.size()); ++i)
    containerPtrs[i]->newECM(eCM);
  for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
    container2Ptrs[i2]->newECM(eCM);

  // Outer loop in case of rare failures.
  bool physical = true;
  for (int loop = 0; loop < MAXLOOP; ++loop) {
    if (!physical) process.clear();
    physical = true;

    // Loop over both hard processes to find consistent common kinematics.
    for ( ; ; ) {

      // Loop internally over tries for hardest process until succeeds.
      for ( ; ; ) {

        // Pick one of the subprocesses.
        double sigmaMaxNow = sigmaMaxSum * rndmPtr->flat();
        int iMax = containerPtrs.size() - 1;
        iContainer = -1;
        do sigmaMaxNow -= containerPtrs[++iContainer]->sigmaMax();
        while (sigmaMaxNow > 0. && iContainer < iMax);

        // Do a trial event of this subprocess; accept or not.
        if (containerPtrs[iContainer]->trialProcess()) break;

        // Check for end-of-file condition for Les Houches events.
        if (infoPtr->atEndOfFile()) return false;
      }

      // Update sum of maxima if current maximum violated.
      if (containerPtrs[iContainer]->newSigmaMax()) {
        sigmaMaxSum = 0.;
        for (int i = 0; i < int(containerPtrs.size()); ++i)
          sigmaMaxSum += containerPtrs[i]->sigmaMax();
      }

      // Loop internally over tries for second hardest process until succeeds.
      for ( ; ; ) {

        // Pick one of the subprocesses.
        double sigma2MaxNow = sigma2MaxSum * rndmPtr->flat();
        int i2Max = container2Ptrs.size() - 1;
        i2Container = -1;
        do sigma2MaxNow -= container2Ptrs[++i2Container]->sigmaMax();
        while (sigma2MaxNow > 0. && i2Container < i2Max);

        // Do a trial event of this subprocess; accept or not.
        if (container2Ptrs[i2Container]->trialProcess()) break;
      }

      // Update sum of maxima if current maximum violated.
      if (container2Ptrs[i2Container]->newSigmaMax()) {
        sigma2MaxSum = 0.;
        for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
          sigma2MaxSum += container2Ptrs[i2]->sigmaMax();
      }

      // Pick incoming flavours (etc), needed for PDF reweighting.
      containerPtrs[iContainer]->constructState();
      container2Ptrs[i2Container]->constructState();

      // Check whether common set of x values is kinematically possible.
      double xA1      = containerPtrs[iContainer]->x1();
      double xB1      = containerPtrs[iContainer]->x2();
      double xA2      = container2Ptrs[i2Container]->x1();
      double xB2      = container2Ptrs[i2Container]->x2();
      if (xA1 + xA2 >= 1. || xB1 + xB2 >= 1.) continue;

      // Reset beam contents. Naive parton densities for second interaction.
      // (Subsequent procedure could be symmetrized, but would be overkill.)
      beamAPtr->clear();
      beamBPtr->clear();
      int    idA1     = containerPtrs[iContainer]->id1();
      int    idB1     = containerPtrs[iContainer]->id2();
      int    idA2     = container2Ptrs[i2Container]->id1();
      int    idB2     = container2Ptrs[i2Container]->id2();
      double Q2Fac1   = containerPtrs[iContainer]->Q2Fac();
      double Q2Fac2   = container2Ptrs[i2Container]->Q2Fac();
      double pdfA2Raw = beamAPtr->xf( idA2, xA2,Q2Fac2);
      double pdfB2Raw = beamBPtr->xf( idB2, xB2,Q2Fac2);

      // Remove partons in first interaction from beams.
      beamAPtr->append( 3, idA1, xA1);
      beamAPtr->xfISR( 0, idA1, xA1, Q2Fac1);
      beamAPtr->pickValSeaComp();
      beamBPtr->append( 4, idB1, xB1);
      beamBPtr->xfISR( 0, idB1, xB1, Q2Fac1);
      beamBPtr->pickValSeaComp();

      // Reevaluate pdf's for second interaction and weight by reduction.
      double pdfA2Mod = beamAPtr->xfMPI( idA2, xA2,Q2Fac2);
      double pdfB2Mod = beamBPtr->xfMPI( idB2, xB2,Q2Fac2);
      double wtPdfMod = (pdfA2Mod * pdfB2Mod) / (pdfA2Raw * pdfB2Raw);
      if (wtPdfMod < rndmPtr->flat()) continue;

      // Reduce by a factor of 2 for identical processes when others not,
      // and when in same phase space region.
      bool toLoop = false;
      if ( someHardSame && containerPtrs[iContainer]->isSame()
        && container2Ptrs[i2Container]->isSame()) {
        if (cutsAgree) {
          if (rndmPtr->flat() > 0.5) toLoop = true;
        } else {
        double mHat1 = containerPtrs[iContainer]->mHat();
        double pTHat1 = containerPtrs[iContainer]->pTHat();
        double mHat2 = container2Ptrs[i2Container]->mHat();
        double pTHat2 = container2Ptrs[i2Container]->pTHat();
        if (mHat1 > mHatMin2 && mHat1 < mHatMax2
           && pTHat1 > pTHatMin2 && pTHat1 < pTHatMax2
           && mHat2 > mHatMin1 && mHat2 < mHatMax1
           && pTHat2 > pTHatMin1 && pTHat2 < pTHatMax1
           && rndmPtr->flat() > 0.5) toLoop = true;
        }
      }
      if (toLoop) continue;

      // If come this far then acceptable event.
      break;
    }

    // Construct kinematics of acceptable processes.
    Event process2;
    process2.init( "(second hard)", particleDataPtr, startColTag);
    process2.initColTag();
    if ( !containerPtrs[iContainer]->constructProcess( process) )
      physical = false;
    if (physical && !container2Ptrs[i2Container]->constructProcess( process2,
      false) ) physical = false;

    // Do all resonance decays.
    if ( physical && doResDecays
      &&  !containerPtrs[iContainer]->decayResonances( process) )
      physical = false;
    if ( physical && doResDecays
      &&  !container2Ptrs[i2Container]->decayResonances( process2) )
      physical = false;

    // Append second hard interaction to normal process object.
    if (physical) combineProcessRecords( process, process2);

    // Add any junctions to the process event record list.
    if (physical) findJunctions( process);

    // Outer loop should normally work first time around.
    if (physical) break;
  }

  // Done.
  return physical;
}

//--------------------------------------------------------------------------

// Append second hard interaction to normal process object.
// Complication: all resonance decay chains must be put at the end.

void ProcessLevel::combineProcessRecords( Event& process, Event& process2) {

  // Find first event record size, excluding resonances.
  int nSize = process.size();
  int nHard = 5;
  while (nHard < nSize && process[nHard].mother1() == 3) ++nHard;

  // Save resonance products temporarily elsewhere.
  vector<Particle> resProd;
  if (nSize > nHard) {
    for (int i = nHard; i < nSize; ++i) resProd.push_back( process[i] );
    process.popBack(nSize - nHard);
  }

  // Find second event record size, excluding resonances.
  int nSize2 = process2.size();
  int nHard2 = 5;
  while (nHard2 < nSize2 && process2[nHard2].mother1() == 3) ++nHard2;

  // Find amount of necessary position and colour offset for second process.
  int addPos  = nHard  - 3;
  int addCol  = process.lastColTag() - startColTag;

  // Loop over all particles (except beams) from second process.
  for (int i = 3; i < nSize2; ++i) {

    // Offset mother and daughter pointers and colour tags of particle.
    process2[i].offsetHistory( 2, addPos, 2, addPos);
    process2[i].offsetCol( addCol);

    // Append hard-process particles from process2 to process.
    if (i < nHard2) process.append( process2[i] );
  }

  // Reinsert resonance decay chains of first hard process.
  int addPos2 = nHard2 - 3;
  if (nHard < nSize) {

    // Offset daughter pointers of unmoved mothers.
    for (int i = 5; i < nHard; ++i)
      process[i].offsetHistory( 0, 0, nHard - 1, addPos2);

    // Modify history of resonance products when restoring.
    for (int i = 0; i < int(resProd.size()); ++i) {
      resProd[i].offsetHistory( nHard - 1, addPos2, nHard - 1, addPos2);
      process.append( resProd[i] );
    }
  }

  // Insert resonance decay chains of second hard process.
  if (nHard2 < nSize2) {
    int nHard3  = nHard + nHard2 - 3;
    int addPos3 = nSize - nHard;

    // Offset daughter pointers of second-process mothers.
    for (int i = nHard + 2; i < nHard3; ++i)
      process[i].offsetHistory( 0, 0, nHard3 - 1, addPos3);

    // Modify history of second-process resonance products and insert.
    for (int i = nHard2; i < nSize2; ++i) {
      process2[i].offsetHistory( nHard3 - 1, addPos3, nHard3 - 1, addPos3);
      process.append( process2[i] );
    }
  }

  // Store PDF scale for second interaction.
  process.scaleSecond( process2.scale() );

}

//--------------------------------------------------------------------------

// Add any junctions to the process event record list.
// Also check that do not doublebook if called repeatedly.

void ProcessLevel::findJunctions( Event& junEvent) {

  // Check all hard vertices for BNV
  for (int i = 1; i<junEvent.size(); i++) {

    // Ignore colorless particles and stages before hard-scattering
    // final state.
    if (abs(junEvent[i].status()) <= 21 || junEvent[i].colType() == 0)
      continue;
    vector<int> motherList   = junEvent[i].motherList();
    int iMot1 = motherList[0];
    vector<int> sisterList = junEvent[iMot1].daughterList();

    // Check baryon number of vertex.
    int barSum = 0;
    map<int,int> colVertex, acolVertex;

    // Loop over mothers (enter with crossed colors and negative sign).
    for (unsigned int indx = 0; indx < motherList.size(); indx++) {
      int iMot = motherList[indx];
      if ( abs(junEvent[iMot].colType()) == 1 )
        barSum -= junEvent[iMot].colType();
      else if ( abs(junEvent[iMot].colType()) == 3)
        barSum -= 2*junEvent[iMot].colType()/3;
      int col  = junEvent[iMot].acol();
      int acol  = junEvent[iMot].col();

      // If unmatched (so far), add end. Else erase matching parton.
      if (col > 0) {
        if (acolVertex.find(col) == acolVertex.end() ) colVertex[col] = iMot;
        else acolVertex.erase(col);
      } else if (col < 0) {
        if (colVertex.find(-col) == colVertex.end() ) acolVertex[-col] = iMot;
        else colVertex.erase(-col);
      }
      if (acol > 0) {
        if (colVertex.find(acol) == colVertex.end()) acolVertex[acol] = iMot;
        else colVertex.erase(acol);
      } else if (acol < 0) {
        if (acolVertex.find(-acol) == acolVertex.end())
             colVertex[-acol] = iMot;
        else acolVertex.erase(-acol);
      }
    }

    // Loop over sisters.
    for (unsigned int indx = 0; indx < sisterList.size(); indx++) {
      int iDau = sisterList[indx];
      if ( abs(junEvent[iDau].colType()) == 1 )
        barSum += junEvent[iDau].colType();
      else if ( abs(junEvent[iDau].colType()) == 3)
        barSum += 2*junEvent[iDau].colType()/3;
      int col  = junEvent[iDau].col();
      int acol  = junEvent[iDau].acol();

      // If unmatched (so far), add end. Else erase matching parton.
      if (col > 0) {
        if (acolVertex.find(col) == acolVertex.end() ) colVertex[col] = iDau;
        else acolVertex.erase(col);
      } else if (col < 0) {
        if (colVertex.find(-col) == colVertex.end() ) acolVertex[-col] = iDau;
        else colVertex.erase(-col);
      }
      if (acol > 0) {
        if (colVertex.find(acol) == colVertex.end()) acolVertex[acol] = iDau;
        else colVertex.erase(acol);
      } else if (acol < 0) {
        if (acolVertex.find(-acol) == acolVertex.end())
             colVertex[-acol] = iDau;
        else acolVertex.erase(-acol);
      }

    }

    // Skip if baryon number conserved in this vertex.
    if (barSum == 0) continue;

    // Check and skip any junctions that have already been added.
    for (int iJun = 0; iJun < junEvent.sizeJunction(); ++iJun) {
      // Remove the tags corresponding to each of the 3 existing junction legs.
      for (int j = 0; j < 3; ++j) {
        int colNow = junEvent.colJunction(iJun, j);
        if (junEvent.kindJunction(iJun) % 2 == 1) colVertex.erase(colNow);
        else acolVertex.erase(colNow);
      }
    }

    // Skip if no junction colors remain.
    if (colVertex.size() == 0 && acolVertex.size() == 0) continue;

    // If baryon number violated, is B = +1 or -1 (larger values not handled).
    int kindJun = 0;
    if (colVertex.size() == 3 && acolVertex.size() == 0) kindJun = 1;
    else if (colVertex.size() == 0 && acolVertex.size() == 3) kindJun = 2;
    else {
      infoPtr->errorMsg("Error in ProcessLevel::findJunctions: "
                        "N(unmatched (anti)colour tags) != 3");
      return;
    }

    // From now on, use colJun as shorthand for colVertex or acolVertex.
    map<int,int> colJun = (kindJun == 1) ? colVertex : acolVertex;

    // Order so incoming tags appear first in colVec, outgoing tags last.
    vector<int> colVec;
    for (map<int,int>::iterator it = colJun.begin();
         it != colJun.end(); it++) {
      int col  = it->first;
      int iCol = it->second;
      for (unsigned int indx = 0; indx < motherList.size(); indx++) {
        if (iCol == motherList[indx]) {
          kindJun += 2;
          colVec.insert(colVec.begin(),col);
        }
      }
      if (colVec.size() == 0 || colVec[0] != col) colVec.push_back(col);
    }

    // Add junction with these tags.
    junEvent.appendJunction( kindJun, colVec[0], colVec[1], colVec[2]);

  }

}
//--------------------------------------------------------------------------

// Check that colours match up.

bool ProcessLevel::checkColours( Event& process) {

  // Variables and arrays for common usage.
  bool physical = true;
  bool match;
  int colType, col, acol, iPos, iNow, iNowA;
  vector<int> colTags, colPos, acolPos;

  // Check that each particle has the kind of colours expected of it.
  for (int i = 0; i < process.size(); ++i) {
    colType = process[i].colType();
    col     = process[i].col();
    acol    = process[i].acol();
    if      (colType ==  0 && (col != 0 || acol != 0)) physical = false;
    else if (colType ==  1 && (col <= 0 || acol != 0)) physical = false;
    else if (colType == -1 && (col != 0 || acol <= 0)) physical = false;
    else if (colType ==  2 && (col <= 0 || acol <= 0)) physical = false;
    // Colour-sextet assignments (colType == 3 for sextet, -3 for antisextet).
    // Sextet (two colours) represented by (colour, negative anticolour) tags.
    // Antisextet (two anticolours) by (negative colour, anticolour) tags.
    else if (colType ==  3 && (col <= 0 || acol >= 0)) physical = false;
    else if (colType == -3 && (col >= 0 || acol <= 0)) physical = false;
    // All other cases
    else if (colType == -2 || colType < -3 || colType > 3) physical = false;

    // Add to the list of colour tags.
    if (col > 0) {
      match = false;
      for (int ic = 0; ic < int(colTags.size()) ; ++ic)
        if (col == colTags[ic]) match = true;
      if (!match) colTags.push_back(col);
    } else if (acol > 0) {
      match = false;
      for (int ic = 0; ic < int(colTags.size()) ; ++ic)
        if (acol == colTags[ic]) match = true;
      if (!match) colTags.push_back(acol);
    }
    // Colour sextets : map negative colour -> anticolour and vice versa
    if (col < 0) {
      match = false;
      for (int ic = 0; ic < int(colTags.size()) ; ++ic)
        if (-col == colTags[ic]) match = true;
      if (!match) colTags.push_back(-col);
    } else if (acol < 0) {
      match = false;
      for (int ic = 0; ic < int(colTags.size()) ; ++ic)
        if (-acol == colTags[ic]) match = true;
      if (!match) colTags.push_back(-acol);
    }
  }

  // Warn and give up if particles did not have the expected colours.
  if (!physical) {
    infoPtr->errorMsg("Error in ProcessLevel::checkColours: "
      "incorrect colour assignment");
    return false;
  }

  // Remove (anti)colours coming from an (anti)junction.
  for (int iJun = 0; iJun < process.sizeJunction(); ++iJun) {
    for (int j = 0; j < 3; ++j) {
      int colJun = process.colJunction(iJun, j);
      for (int ic = 0; ic < int(colTags.size()) ; ++ic)
        if (colJun == colTags[ic]) {
          colTags[ic] = colTags[colTags.size() - 1];
          colTags.pop_back();
          break;
        }
    }
  }

  // Loop through all colour tags and find their positions (by sign).
  for (int ic = 0; ic < int(colTags.size()); ++ic) {
    col = colTags[ic];
    colPos.resize(0);
    acolPos.resize(0);
    for (int i = 0; i < process.size(); ++i) {
      if (process[i].col() == col || process[i].acol() == -col)
        colPos.push_back(i);
      if (process[i].acol() == col || process[i].col() == -col)
        acolPos.push_back(i);
    }

    // Trace colours back through decays; remove daughters.
    while (colPos.size() > 1) {
      iPos = colPos.size() - 1;
      iNow = colPos[iPos];
      if ( process[iNow].mother1() == colPos[iPos - 1]
        && process[iNow].mother2() == 0) colPos.pop_back();
      else break;
    }
    while (acolPos.size() > 1) {
      iPos = acolPos.size() - 1;
      iNow = acolPos[iPos];
      if ( process[iNow].mother1() == acolPos[iPos - 1]
        && process[iNow].mother2() == 0) acolPos.pop_back();
      else break;
    }

    // Now colour should exist in only 2 copies.
    if (colPos.size() + acolPos.size() != 2) physical = false;

    // If both colours or both anticolours then one mother of the other.
    else if (colPos.size() == 2) {
      iNow = colPos[1];
      if ( process[iNow].mother1() != colPos[0]
        && process[iNow].mother2() != colPos[0] ) physical = false;
    }
    else if (acolPos.size() == 2) {
      iNowA = acolPos[1];
      if ( process[iNowA].mother1() != acolPos[0]
        && process[iNowA].mother2() != acolPos[0] ) physical = false;
    }

    // If one of each then should have same mother(s), or point to beams.
    else {
      iNow  = colPos[0];
      iNowA = acolPos[0];
      if ( process[iNow].status() == -21 &&  process[iNowA].status() == -21 );
      else if ( (process[iNow].mother1() != process[iNowA].mother1())
             || (process[iNow].mother2() != process[iNowA].mother2()) )
        physical = false;
    }

  }

  // Error message if problem found. Done.
  if (!physical) infoPtr->errorMsg("Error in ProcessLevel::checkColours: "
                   "unphysical colour flow");
  return physical;

}

//--------------------------------------------------------------------------

// Print statistics when two hard processes allowed.

void ProcessLevel::statistics2(bool reset, ostream& os) {

  // Average impact-parameter factor and error.
  double invN          = 1. / max(1, nImpact);
  double impactFac     = max( 1., sumImpactFac * invN);
  double impactErr2    = ( sum2ImpactFac * invN / pow2(impactFac) - 1.) * invN;

  // Derive scaling factor to be applied to first set of processes.
  double sigma2SelSum  = 0.;
  int    n2SelSum      = 0;
  for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2) {
    sigma2SelSum      += container2Ptrs[i2]->sigmaSelMC();
    n2SelSum          += container2Ptrs[i2]->nSelected();
  }
  double factor1       = impactFac * sigma2SelSum / sigmaND;
  double rel1Err       = sqrt(1. / max(1, n2SelSum) + impactErr2);
  if (allHardSame) factor1 *= 0.5;

  // Derive scaling factor to be applied to second set of processes.
  double sigma1SelSum  = 0.;
  int    n1SelSum      = 0;
  for (int i = 0; i < int(containerPtrs.size()); ++i) {
    sigma1SelSum      += containerPtrs[i]->sigmaSelMC();
    n1SelSum          += containerPtrs[i]->nSelected();
  }
  double factor2       = impactFac * sigma1SelSum / sigmaND;
  if (allHardSame) factor2 *= 0.5;
  double rel2Err       = sqrt(1. / max(1, n1SelSum) + impactErr2);

  // Header.
  os << "\n *-------  PYTHIA Event and Cross Section Statistics  ------"
     << "--------------------------------------------------*\n"
     << " |                                                            "
     << "                                                |\n"
     << " | Subprocess                               Code |            "
     << "Number of events       |      sigma +- delta    |\n"
     << " |                                               |       Tried"
     << "   Selected   Accepted |     (estimated) (mb)   |\n"
     << " |                                               |            "
     << "                       |                        |\n"
     << " |------------------------------------------------------------"
     << "------------------------------------------------|\n"
     << " |                                               |            "
     << "                       |                        |\n"
     << " | First hard process:                           |            "
     << "                       |                        |\n"
     << " |                                               |            "
     << "                       |                        |\n";

  // Reset sum counters.
  long   nTrySum   = 0;
  long   nSelSum   = 0;
  long   nAccSum   = 0;
  double sigmaSum  = 0.;
  double delta2Sum = 0.;

  // Reset process maps.
  map<int, string> nameM;
  map<int, long> nTryM, nSelM, nAccM;
  map<int, double> sigmaM, delta2M;

  // Loop over existing first processes.
  for (int i = 0; i < int(containerPtrs.size()); ++i)
  if (containerPtrs[i]->sigmaMax() != 0.) {

    // Read info for process. Sum counters.
    int code = containerPtrs[i]->code();
    nTrySum       += containerPtrs[i]->nTried();
    nSelSum       += containerPtrs[i]->nSelected();
    nAccSum       += containerPtrs[i]->nAccepted();
    sigmaSum      += containerPtrs[i]->sigmaMC() * factor1;
    delta2Sum     += pow2(containerPtrs[i]->deltaMC() * factor1);
    nameM[code]    = containerPtrs[i]->name();
    nTryM[code]   += containerPtrs[i]->nTried();
    nSelM[code]   += containerPtrs[i]->nSelected();
    nAccM[code]   += containerPtrs[i]->nAccepted();
    sigmaM[code]  += containerPtrs[i]->sigmaMC() * factor1;
    delta2M[code] += pow2(containerPtrs[i]->deltaMC() * factor1);
    delta2M[code] += pow2(containerPtrs[i]->sigmaMC() * factor1 * rel1Err);
  }

  // Print first process info.
  for (map<int, string>::iterator i = nameM.begin(); i != nameM.end(); ++i) {
    int code = i->first;
    os << " | " << left << setw(40) << i->second
       << right << setw(5) << code << " | "
       << setw(11) << nTryM[code] << " " << setw(10) << nSelM[code] << " "
       << setw(10) << nAccM[code] << " | " << scientific << setprecision(3)
       << setw(11) << sigmaM[code]
       << setw(11) << sqrtpos(delta2M[code]) << " |\n";
  }

  // Print summed info for first processes.
  delta2Sum       += pow2( sigmaSum * rel1Err );
  os << " |                                               |            "
     << "                       |                        |\n"
     << " | " << left << setw(45) << "sum" << right << " | " << setw(11)
     << nTrySum << " " << setw(10) << nSelSum << " " << setw(10)
     << nAccSum << " | " << scientific << setprecision(3) << setw(11)
     << sigmaSum << setw(11) << sqrtpos(delta2Sum) << " |\n";

  // Separation lines to second hard processes.
  os << " |                                               |            "
     << "                       |                        |\n"
     << " |------------------------------------------------------------"
     << "------------------------------------------------|\n"
     << " |                                               |            "
     << "                       |                        |\n"
     << " | Second hard process:                          |            "
     << "                       |                        |\n"
     << " |                                               |            "
     << "                       |                        |\n";

  // Reset sum counters.
  nTrySum   = 0;
  nSelSum   = 0;
  nAccSum   = 0;
  sigmaSum  = 0.;
  delta2Sum = 0.;

  // Reset process maps.
  nameM.clear();
  nTryM.clear();
  nSelM.clear();
  nAccM.clear();
  sigmaM.clear();
  delta2M.clear();

  // Loop over existing second processes.
  for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
  if (container2Ptrs[i2]->sigmaMax() != 0.) {

    // Read info for process. Sum counters.
    int code = container2Ptrs[i2]->code();
    nTrySum       += container2Ptrs[i2]->nTried();
    nSelSum       += container2Ptrs[i2]->nSelected();
    nAccSum       += container2Ptrs[i2]->nAccepted();
    sigmaSum      += container2Ptrs[i2]->sigmaMC() * factor2;
    delta2Sum     += pow2(container2Ptrs[i2]->deltaMC() * factor2);
    nameM[code]    = container2Ptrs[i2]->name();
    nTryM[code]   += container2Ptrs[i2]->nTried();
    nSelM[code]   += container2Ptrs[i2]->nSelected();
    nAccM[code]   += container2Ptrs[i2]->nAccepted();
    sigmaM[code]  += container2Ptrs[i2]->sigmaMC() * factor2;
    delta2M[code] += pow2(container2Ptrs[i2]->deltaMC() * factor2);
    delta2M[code] += pow2(container2Ptrs[i2]->sigmaMC() * factor2 * rel2Err);
  }

  // Print second process info.
  for (map<int, string>::iterator i2 = nameM.begin(); i2 != nameM.end();
    ++i2) {
    int code = i2->first;
    os << " | " << left << setw(40) << i2->second
       << right << setw(5) << code << " | "
       << setw(11) << nTryM[code] << " " << setw(10) << nSelM[code] << " "
       << setw(10) << nAccM[code] << " | " << scientific << setprecision(3)
       << setw(11) << sigmaM[code]
       << setw(11) << sqrtpos(delta2M[code]) << " |\n";
  }

  // Print summed info for second processes.
  delta2Sum       += pow2( sigmaSum * rel2Err );
  os << " |                                               |            "
     << "                       |                        |\n"
     << " | " << left << setw(45) << "sum" << right << " | " << setw(11)
     << nTrySum << " " << setw(10) << nSelSum << " " << setw(10)
     << nAccSum << " | " << scientific << setprecision(3) << setw(11)
     << sigmaSum << setw(11) << sqrtpos(delta2Sum) << " |\n";

  // Print information on how the two processes were combined.
  os << " |                                               |            "
     << "                       |                        |\n"
     << " |------------------------------------------------------------"
     << "------------------------------------------------|\n"
     << " |                                                            "
     << "                                                |\n"
     << " | Uncombined cross sections for the two event sets were "
     << setw(10) << sigma1SelSum << " and " << sigma2SelSum << " mb, "
     << "respectively, combined  |\n"
     << " | using a sigma(nonDiffractive) of " << setw(10) << sigmaND
     << " mb and an impact-parameter enhancement factor of "
     << setw(10) << impactFac << ".   |\n";

  // Listing finished.
  os << " |                                                            "
     << "                                                |\n"
     << " *-------  End PYTHIA Event and Cross Section Statistics -----"
     << "------------------------------------------------*" << endl;

  // Optionally reset statistics contants.
  if (reset) resetStatistics();

}

//==========================================================================

} // end namespace Pythia8
