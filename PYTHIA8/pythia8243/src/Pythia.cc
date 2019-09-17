// Pythia.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Pythia class.

#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"

// Access time information.
#include <ctime>

// Allow string and character manipulation.
#include <cctype>

namespace Pythia8 {

//==========================================================================

// The Pythia class.

//--------------------------------------------------------------------------

// The current Pythia (sub)version number, to agree with XML version.
const double Pythia::VERSIONNUMBERHEAD = PYTHIA_VERSION;
const double Pythia::VERSIONNUMBERCODE = 8.243;

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of tries to produce parton level from given input.
const int Pythia::NTRY          = 10;

// Negative integer to denote that no subrun has been set.
const int Pythia::SUBRUNDEFAULT = -999;

//--------------------------------------------------------------------------

// Constructor.

Pythia::Pythia(string xmlDir, bool printBanner) :
  couplingsPtr(), mergingPtr(),
  mergingHooksPtr(), heavyIonsPtr(), hiHooksPtr(), doProcessLevel(),
  doPartonLevel(), doHadronLevel(), doSoftQCDall(), doSoftQCDinel(),
  doCentralDiff(), doDiffraction(), doSoftQCD(), doVMDsideA(), doVMDsideB(),
  doHardDiff(), doResDec(), doFSRinRes(), decayRHadrons(), abortIfVeto(),
  checkEvent(), checkHistory(), nErrList(), epTolErr(), epTolWarn(), mTolErr(),
  mTolWarn(), beamHasGamma(), beamAisResGamma(), beamBisResGamma(),
  beamAhasResGamma(), beamBhasResGamma(), gammaMode(), isUnresolvedA(),
  isUnresolvedB(), showSaV(), showMaD(), doReconnect(), forceHadronLevelCR(),
  idA(), idB(), frameType(), boostType(), nCount(), nShowLHA(), nShowInfo(),
  nShowProc(), nShowEvt(), reconnectMode(), mA(), mB(), pxA(), pxB(), pyA(),
  pyB(), pzA(), pzB(), eA(), eB(), pzAcm(), pzBcm(), eCM(), betaZ(), gammaZ(),
  nErrEvent(), pdfAPtr(), pdfBPtr(), pdfHardAPtr(), pdfHardBPtr(),
  pdfPomAPtr(), pdfPomBPtr(), pdfGamAPtr(), pdfGamBPtr(), pdfHardGamAPtr(),
  pdfHardGamBPtr(), pdfUnresAPtr(), pdfUnresBPtr(), pdfUnresGamAPtr(),
  pdfUnresGamBPtr(), pdfGamFluxAPtr(), pdfGamFluxBPtr(), pdfVMDAPtr(),
  pdfVMDBPtr(), useNewPdfA(), useNewPdfB(), useNewPdfHard(), useNewPdfPomA(),
  useNewPdfPomB(), useNewPdfGamA(), useNewPdfGamB(), useNewPdfHardGamA(),
  useNewPdfHardGamB(), useNewPdfUnresA(), useNewPdfUnresB(),
  useNewPdfUnresGamA(), useNewPdfUnresGamB(), useNewPdfVMDA(), useNewPdfVMDB(),
  hasUserHooksVector(), doLHA(), useNewLHA(), lhaUpPtr(), decayHandlePtr(),
  userHooksPtr(), hasUserHooks(), doVetoProcess(), doVetoPartons(),
  retryPartonLevel(), beamShapePtr(), useNewBeamShape(), doMomentumSpread(),
  doVertexSpread(), doVarEcm(), eMinPert(), eWidthPert(), timesDecPtr(),
  timesPtr(), spacePtr(), useNewTimesDec(), useNewTimes(), useNewSpace(),
  partonVertexPtr(), useNewPartonVertex(), hasMerging(), hasOwnMerging(),
  hasMergingHooks(), hasOwnMergingHooks(), doMerging(), hasHeavyIons(),
  hasOwnHeavyIons(), doHeavyIons() {

  // Initialise / reset pointers and global variables.
  initPtrs();

  // Find path to data files, i.e. xmldoc directory location.
  // Environment variable takes precedence, then constructor input,
  // and finally the pre-processor constant XMLDIR.
  xmlPath = "";
  const char* PYTHIA8DATA = "PYTHIA8DATA";
  char* envPath = getenv(PYTHIA8DATA);
  if (envPath != 0 && *envPath != '\0') {
    int i = 0;
    while (*(envPath+i) != '\0') xmlPath += *(envPath+(i++));
  } else {
    if (xmlDir[ xmlDir.length() - 1 ] != '/') xmlDir += "/";
    xmlPath = xmlDir;
    ifstream xmlFile((xmlPath + "Index.xml").c_str());
    //if (!xmlFile.good()) xmlPath = XMLDIR;
    xmlFile.close();
  }
  if (xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";

  // Read in files with all flags, modes, parms and words.
  settings.initPtr( &info);
  string initFile = xmlPath + "Index.xml";
  isConstructed = settings.init( initFile);
  if (!isConstructed) {
    info.errorMsg("Abort from Pythia::Pythia: settings unavailable");
    return;
  }

  // Save XML path in settings.
  settings.addWord( "xmlPath", xmlPath);

  // Check that XML and header version numbers match code version number.
  if (!checkVersion()) return;

  // Read in files with all particle data.
  particleData.initPtr( &info, &settings, &rndm, couplingsPtr);
  string dataFile = xmlPath + "ParticleData.xml";
  isConstructed = particleData.init( dataFile);
  if (!isConstructed) {
    info.errorMsg("Abort from Pythia::Pythia: particle data unavailable");
    return;
  }

  // Write the Pythia banner to output.
  if (printBanner) banner();

  // Not initialized until at the end of the init() call.
  isInit = false;
  info.addCounter(0);

  HeavyIons::addSpecialSettings(settings);

}

//--------------------------------------------------------------------------

// Constructor from pre-initialised ParticleData and Settings objects.

Pythia::Pythia(Settings& settingsIn, ParticleData& particleDataIn,
  bool printBanner) :
  couplingsPtr(), mergingPtr(), mergingHooksPtr(),
  heavyIonsPtr(), hiHooksPtr(), doProcessLevel(), doPartonLevel(),
  doHadronLevel(), doSoftQCDall(), doSoftQCDinel(), doCentralDiff(),
  doDiffraction(), doSoftQCD(), doVMDsideA(), doVMDsideB(), doHardDiff(),
  doResDec(), doFSRinRes(), decayRHadrons(), abortIfVeto(), checkEvent(),
  checkHistory(), nErrList(), epTolErr(), epTolWarn(), mTolErr(), mTolWarn(),
  beamHasGamma(), beamAisResGamma(), beamBisResGamma(), beamAhasResGamma(),
  beamBhasResGamma(), gammaMode(), isUnresolvedA(), isUnresolvedB(),
  showSaV(), showMaD(), doReconnect(), forceHadronLevelCR(), idA(), idB(),
  frameType(), boostType(), nCount(), nShowLHA(), nShowInfo(), nShowProc(),
  nShowEvt(), reconnectMode(), mA(), mB(), pxA(), pxB(), pyA(), pyB(),
  pzA(), pzB(), eA(), eB(), pzAcm(), pzBcm(), eCM(), betaZ(), gammaZ(),
  nErrEvent(), pdfAPtr(), pdfBPtr(), pdfHardAPtr(), pdfHardBPtr(),
  pdfPomAPtr(), pdfPomBPtr(), pdfGamAPtr(), pdfGamBPtr(), pdfHardGamAPtr(),
  pdfHardGamBPtr(), pdfUnresAPtr(), pdfUnresBPtr(), pdfUnresGamAPtr(),
  pdfUnresGamBPtr(), pdfGamFluxAPtr(), pdfGamFluxBPtr(), pdfVMDAPtr(),
  pdfVMDBPtr(), useNewPdfA(), useNewPdfB(), useNewPdfHard(), useNewPdfPomA(),
  useNewPdfPomB(), useNewPdfGamA(), useNewPdfGamB(), useNewPdfHardGamA(),
  useNewPdfHardGamB(), useNewPdfUnresA(), useNewPdfUnresB(),
  useNewPdfUnresGamA(), useNewPdfUnresGamB(), useNewPdfVMDA(), useNewPdfVMDB(),
  hasUserHooksVector(), doLHA(), useNewLHA(), lhaUpPtr(), decayHandlePtr(),
  userHooksPtr(), hasUserHooks(), doVetoProcess(), doVetoPartons(),
  retryPartonLevel(), beamShapePtr(), useNewBeamShape(), doMomentumSpread(),
  doVertexSpread(), doVarEcm(), eMinPert(), eWidthPert(), timesDecPtr(),
  timesPtr(), spacePtr(), useNewTimesDec(), useNewTimes(), useNewSpace(),
  partonVertexPtr(), useNewPartonVertex(), hasMerging(), hasOwnMerging(),
  hasMergingHooks(), hasOwnMergingHooks(), doMerging(),hasHeavyIons(),
  hasOwnHeavyIons(), doHeavyIons() {

  // Initialise / reset pointers and global variables.
  initPtrs();

  // Copy XML path from existing Settings database.
  xmlPath = settingsIn.word("xmlPath");

  // Copy settings database and redirect pointers.
  settings = settingsIn;
  settings.initPtr( &info);
  isConstructed = settings.getIsInit();
  if (!isConstructed) {
    info.errorMsg("Abort from Pythia::Pythia: settings unavailable");
    return;
  }

  // Check XML and header version numbers match code version number.
  if (!checkVersion()) return;

  // Copy particleData database and redirect pointers.
  particleData = particleDataIn;
  particleData.initPtr( &info, &settings, &rndm, couplingsPtr);
  isConstructed = particleData.getIsInit();
  if (!isConstructed) {
    info.errorMsg("Abort from Pythia::Pythia: particle data unavailable");
    return;
  }

  // Write the Pythia banner to output.
  if (printBanner) banner();

  // Not initialized until at the end of the init() call.
  isInit = false;
  info.addCounter(0);

}

//--------------------------------------------------------------------------

// Constructor from string streams.

Pythia::Pythia( istream& settingsStrings, istream& particleDataStrings,
  bool printBanner) :
  couplingsPtr(), mergingPtr(), mergingHooksPtr(), heavyIonsPtr(),
  hiHooksPtr(), doProcessLevel(), doPartonLevel(), doHadronLevel(),
  doSoftQCDall(), doSoftQCDinel(), doCentralDiff(), doDiffraction(),
  doSoftQCD(), doVMDsideA(), doVMDsideB(), doHardDiff(), doResDec(),
  doFSRinRes(), decayRHadrons(), abortIfVeto(), checkEvent(), checkHistory(),
  nErrList(), epTolErr(), epTolWarn(), mTolErr(), mTolWarn(), beamHasGamma(),
  beamAisResGamma(), beamBisResGamma(), beamAhasResGamma(), beamBhasResGamma(),
  gammaMode(), isUnresolvedA(), isUnresolvedB(), showSaV(), showMaD(),
  doReconnect(), forceHadronLevelCR(), idA(), idB(), frameType(), boostType(),
  nCount(), nShowLHA(), nShowInfo(), nShowProc(), nShowEvt(), reconnectMode(),
  mA(), mB(), pxA(), pxB(), pyA(), pyB(), pzA(), pzB(), eA(), eB(), pzAcm(),
  pzBcm(), eCM(), betaZ(), gammaZ(), nErrEvent(), pdfAPtr(), pdfBPtr(),
  pdfHardAPtr(), pdfHardBPtr(), pdfPomAPtr(), pdfPomBPtr(), pdfGamAPtr(),
  pdfGamBPtr(), pdfHardGamAPtr(), pdfHardGamBPtr(), pdfUnresAPtr(),
  pdfUnresBPtr(), pdfUnresGamAPtr(), pdfUnresGamBPtr(), pdfGamFluxAPtr(),
  pdfGamFluxBPtr(), pdfVMDAPtr(), pdfVMDBPtr(), useNewPdfA(), useNewPdfB(),
  useNewPdfHard(), useNewPdfPomA(), useNewPdfPomB(), useNewPdfGamA(),
  useNewPdfGamB(), useNewPdfHardGamA(), useNewPdfHardGamB(), useNewPdfUnresA(),
  useNewPdfUnresB(), useNewPdfUnresGamA(), useNewPdfUnresGamB(),
  useNewPdfVMDA(), useNewPdfVMDB(), hasUserHooksVector(), doLHA(),
  useNewLHA(), lhaUpPtr(), decayHandlePtr(), userHooksPtr(), hasUserHooks(),
  doVetoProcess(), doVetoPartons(), retryPartonLevel(), beamShapePtr(),
  useNewBeamShape(), doMomentumSpread(), doVertexSpread(), doVarEcm(),
  eMinPert(), eWidthPert(), timesDecPtr(), timesPtr(), spacePtr(),
  useNewTimesDec(), useNewTimes(), useNewSpace(), partonVertexPtr(),
  useNewPartonVertex(), hasMerging(), hasOwnMerging(), hasMergingHooks(),
  hasOwnMergingHooks(), doMerging(), hasHeavyIons(), hasOwnHeavyIons(),
  doHeavyIons() {

  // Initialise / reset pointers and global variables.
  initPtrs();

  // Copy settings database
  settings.init( settingsStrings );

  // Reset pointers to pertain to this PYTHIA object.
  settings.initPtr( &info);
  isConstructed = settings.getIsInit();
  if (!isConstructed) {
    info.errorMsg("Abort from Pythia::Pythia: settings unavailable");
    return;
  }

  // Check XML and header version numbers match code version number.
  if (!checkVersion()) return;

  // Read in files with all particle data.
  particleData.initPtr( &info, &settings, &rndm, couplingsPtr);
  isConstructed = particleData.init( particleDataStrings );
  if (!isConstructed) {
    info.errorMsg("Abort from Pythia::Pythia: particle data unavailable");
    return;
  }

  // Write the Pythia banner to output.
  if (printBanner) banner();

  // Not initialized until at the end of the init() call.
  isInit = false;
  info.addCounter(0);

}

//--------------------------------------------------------------------------

// Destructor.

Pythia::~Pythia() {

  // Delete the PDF's created with new.
  if (useNewPdfHard && pdfHardAPtr != pdfAPtr) delete pdfHardAPtr;
  if (useNewPdfHard && pdfHardBPtr != pdfBPtr) delete pdfHardBPtr;
  if (useNewPdfA) delete pdfAPtr;
  if (useNewPdfB) delete pdfBPtr;
  if (useNewPdfPomA) delete pdfPomAPtr;
  if (useNewPdfPomB) delete pdfPomBPtr;
  if (useNewPdfGamA) delete pdfGamAPtr;
  if (useNewPdfGamB) delete pdfGamBPtr;
  if (useNewPdfUnresA) delete pdfUnresAPtr;
  if (useNewPdfUnresB) delete pdfUnresBPtr;
  if (useNewPdfUnresGamA) delete pdfUnresGamAPtr;
  if (useNewPdfUnresGamB) delete pdfUnresGamBPtr;
  if (useNewPdfVMDA) delete pdfVMDAPtr;
  if (useNewPdfVMDB) delete pdfVMDBPtr;

  // Delete the Les Houches object created with new.
  if (useNewLHA) delete lhaUpPtr;

  // Delete vector of UserHooks (but not the UserHooks themselves).
  if (hasUserHooksVector) delete userHooksPtr;

  // Delete the Merging object created with new.
  if (hasOwnMerging) delete mergingPtr;

  // Delete the MergingHooks object created with new.
  if (hasOwnMergingHooks) delete mergingHooksPtr;

  // Delete the HeavyIons object created with new.
  if (hasOwnHeavyIons) delete heavyIonsPtr;

  // Delete the BeamShape object created with new.
  if (useNewBeamShape) delete beamShapePtr;

  // Delete the timelike and spacelike showers created with new.
  if (useNewTimesDec) delete timesDecPtr;
  if (useNewTimes && !useNewTimesDec) delete timesPtr;
  if (useNewSpace) delete spacePtr;

  // Delete the parton vertex object created with new.
  if (useNewPartonVertex) delete partonVertexPtr;

}

//--------------------------------------------------------------------------

// Initialise new Pythia object (common code called by constructors).

void Pythia::initPtrs() {

  // Initial values for pointers to PDF's.
  useNewPdfA         = false;
  useNewPdfB         = false;
  useNewPdfHard      = false;
  useNewPdfPomA      = false;
  useNewPdfPomB      = false;
  useNewPdfGamA      = false;
  useNewPdfGamB      = false;
  useNewPdfHardGamA  = false;
  useNewPdfHardGamB  = false;
  useNewPdfUnresA    = false;
  useNewPdfUnresB    = false;
  useNewPdfUnresGamA = false;
  useNewPdfUnresGamB = false;
  useNewPdfVMDA      = false;
  useNewPdfVMDB      = false;
  pdfAPtr            = 0;
  pdfBPtr            = 0;
  pdfHardAPtr        = 0;
  pdfHardBPtr        = 0;
  pdfPomAPtr         = 0;
  pdfPomBPtr         = 0;
  pdfGamAPtr         = 0;
  pdfGamBPtr         = 0;
  pdfHardGamAPtr     = 0;
  pdfHardGamBPtr     = 0;
  pdfUnresAPtr       = 0;
  pdfUnresBPtr       = 0;
  pdfUnresGamAPtr    = 0;
  pdfUnresGamBPtr    = 0;
  pdfVMDAPtr         = 0;
  pdfVMDBPtr         = 0;

  // Initial pointers to externally provided photon fluxes.
  pdfGamFluxAPtr    = 0;
  pdfGamFluxBPtr    = 0;

  // Initial values for pointers to Les Houches Event objects.
  doLHA              = false;
  useNewLHA          = false;
  lhaUpPtr           = 0;

  //Initial value for couplings pointer
  couplingsPtr       = &couplings;

  // Initial value for pointer to external decay handler.
  decayHandlePtr     = 0;

  // Initial value for pointer to user hooks.
  userHooksPtr       = 0;
  hasUserHooksVector = false;

  // Initial value for pointer to merging hooks.
  doMerging          = false;
  hasMerging         = false;
  hasOwnMerging      = false;
  mergingPtr         = 0;
  hasMergingHooks    = false;
  hasOwnMergingHooks = false;
  mergingHooksPtr    = 0;

  // Initial value for pointer to HeavyIons objects.
  doHeavyIons        = false;
  hasHeavyIons       = false;
  hasOwnHeavyIons    = false;
  heavyIonsPtr       = 0;
  hiHooksPtr         = 0;

  // Initial value for pointer to beam shape.
  useNewBeamShape    = false;
  beamShapePtr       = 0;

  // Initial values for pointers to timelike and spacelike showers.
  useNewTimesDec     = false;
  useNewTimes        = false;
  useNewSpace        = false;
  timesDecPtr        = 0;
  timesPtr           = 0;
  spacePtr           = 0;

  // Initial values for pointer parton-vertex setting.
  useNewPartonVertex = false;
  partonVertexPtr    = 0;

}

//--------------------------------------------------------------------------

// Check for consistency of version numbers (called by constructors).

bool Pythia::checkVersion() {

  // Check that XML version number matches code version number.
  double versionNumberXML = parm("Pythia:versionNumber");
  isConstructed = (abs(versionNumberXML - VERSIONNUMBERCODE) < 0.0005);
  if (!isConstructed) {
    ostringstream errCode;
    errCode << fixed << setprecision(3) << ": in code " << VERSIONNUMBERCODE
            << " but in XML " << versionNumberXML;
    info.errorMsg("Abort from Pythia::Pythia: unmatched version numbers",
      errCode.str());
    return false;
  }

  // Check that header version number matches code version number.
  isConstructed = (abs(VERSIONNUMBERHEAD - VERSIONNUMBERCODE) < 0.0005);
  if (!isConstructed) {
    ostringstream errCode;
    errCode << fixed << setprecision(3) << ": in code " << VERSIONNUMBERCODE
            << " but in header " << VERSIONNUMBERHEAD;
    info.errorMsg("Abort from Pythia::Pythia: unmatched version numbers",
      errCode.str());
    return false;
  }

  // All is well that ends well.
  return true;

}

//--------------------------------------------------------------------------

// Read in one update for a setting or particle data from a single line.

bool Pythia::readString(string line, bool warn) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // If empty line then done.
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return true;

  // If Settings input stretches over several lines then continue with it.
  if (settings.unfinishedInput()) return settings.readString(line, warn);

  // If first character is not a letter/digit, then taken to be a comment.
  int firstChar = line.find_first_not_of(" \n\t\v\b\r\f\a");
  if (!isalnum(line[firstChar])) return true;

  // Send on particle data to the ParticleData database.
  if (isdigit(line[firstChar])) {
    bool passed = particleData.readString(line, warn);
    if (passed) particleDataBuffer << line << endl;
    return passed;
  }

  // Everything else sent on to Settings.
  return settings.readString(line, warn);

}

//--------------------------------------------------------------------------

// Read in updates for settings or particle data from user-defined file.

bool Pythia::readFile(string fileName, bool warn, int subrun) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Open file for reading.
  const char* cstring = fileName.c_str();
  ifstream is(cstring);
  if (!is.good()) {
    info.errorMsg("Error in Pythia::readFile: did not find file", fileName);
    return false;
  }

  // Hand over real work to next method.
  return readFile( is, warn, subrun);

}

//--------------------------------------------------------------------------

// Read in updates for settings or particle data
// from user-defined stream (or file).

bool Pythia::readFile(istream& is, bool warn, int subrun) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Read in one line at a time.
  string line;
  bool isCommented = false;
  bool accepted = true;
  int subrunNow = SUBRUNDEFAULT;
  while ( getline(is, line) ) {

    // Check whether entering, leaving or inside commented-commands section.
    int commentLine = readCommented( line);
    if      (commentLine == +1)  isCommented = true;
    else if (commentLine == -1)  isCommented = false;
    else if (isCommented) ;

    else {
      // Check whether entered new subrun.
      int subrunLine = readSubrun( line, warn);
      if (subrunLine >= 0) subrunNow = subrunLine;

      // Process the line if in correct subrun.
      if ( (subrunNow == subrun || subrunNow == SUBRUNDEFAULT)
         && !readString( line, warn) ) accepted = false;
    }

  // Reached end of input file.
  };
  return accepted;

}

//--------------------------------------------------------------------------

// Routine to pass in pointers to PDF's. Usage optional.

bool Pythia::setPDFPtr( PDF* pdfAPtrIn, PDF* pdfBPtrIn, PDF* pdfHardAPtrIn,
  PDF* pdfHardBPtrIn, PDF* pdfPomAPtrIn, PDF* pdfPomBPtrIn,
  PDF* pdfGamAPtrIn, PDF* pdfGamBPtrIn, PDF* pdfHardGamAPtrIn,
  PDF* pdfHardGamBPtrIn, PDF* pdfUnresAPtrIn, PDF* pdfUnresBPtrIn,
  PDF* pdfUnresGamAPtrIn, PDF* pdfUnresGamBPtrIn, PDF* pdfVMDAPtrIn,
  PDF* pdfVMDBPtrIn) {

  // Delete any PDF's created in a previous init call.
  if (useNewPdfHard && pdfHardAPtr != pdfAPtr) delete pdfHardAPtr;
  if (useNewPdfHard && pdfHardBPtr != pdfBPtr) delete pdfHardBPtr;
  if (useNewPdfA) delete pdfAPtr;
  if (useNewPdfB) delete pdfBPtr;
  if (useNewPdfPomA) delete pdfPomAPtr;
  if (useNewPdfPomB) delete pdfPomBPtr;
  if (useNewPdfGamA) delete pdfGamAPtr;
  if (useNewPdfGamB) delete pdfGamBPtr;
  if (useNewPdfUnresA) delete pdfUnresAPtr;
  if (useNewPdfUnresB) delete pdfUnresBPtr;
  if (useNewPdfUnresGamA) delete pdfUnresGamAPtr;
  if (useNewPdfUnresGamB) delete pdfUnresGamBPtr;
  if (useNewPdfHardGamA && pdfHardGamAPtr != pdfGamAPtr) delete pdfHardGamAPtr;
  if (useNewPdfHardGamB && pdfHardGamBPtr != pdfGamBPtr) delete pdfHardGamBPtr;
  if (useNewPdfVMDA) delete pdfVMDAPtr;
  if (useNewPdfVMDB) delete pdfVMDBPtr;

  // Reset pointers to be empty.
  useNewPdfA         = false;
  useNewPdfB         = false;
  useNewPdfHard      = false;
  useNewPdfPomA      = false;
  useNewPdfPomB      = false;
  useNewPdfGamA      = false;
  useNewPdfGamB      = false;
  useNewPdfHardGamA  = false;
  useNewPdfHardGamB  = false;
  useNewPdfUnresA    = false;
  useNewPdfUnresB    = false;
  useNewPdfUnresGamA = false;
  useNewPdfUnresGamB = false;
  useNewPdfVMDA      = false;
  useNewPdfVMDB      = false;
  pdfAPtr            = 0;
  pdfBPtr            = 0;
  pdfHardAPtr        = 0;
  pdfHardBPtr        = 0;
  pdfPomAPtr         = 0;
  pdfPomBPtr         = 0;
  pdfGamAPtr         = 0;
  pdfGamBPtr         = 0;
  pdfHardGamAPtr     = 0;
  pdfHardGamBPtr     = 0;
  pdfUnresAPtr       = 0;
  pdfUnresBPtr       = 0;
  pdfUnresGamAPtr    = 0;
  pdfUnresGamBPtr    = 0;
  pdfVMDAPtr         = 0;
  pdfVMDBPtr         = 0;

  // Switch off external PDF's by zero as input.
  if (pdfAPtrIn == 0 && pdfBPtrIn == 0) return true;

  // The two PDF objects cannot be one and the same.
  if (pdfAPtrIn == pdfBPtrIn) return false;

  // Save pointers.
  pdfAPtr       = pdfAPtrIn;
  pdfBPtr       = pdfBPtrIn;

  // By default same pointers for hard-process PDF's.
  pdfHardAPtr   = pdfAPtrIn;
  pdfHardBPtr   = pdfBPtrIn;

  // Optionally allow separate pointers for hard process.
  if (pdfHardAPtrIn != 0 && pdfHardBPtrIn != 0) {
    if (pdfHardAPtrIn == pdfHardBPtrIn) return false;
    pdfHardAPtr = pdfHardAPtrIn;
    pdfHardBPtr = pdfHardBPtrIn;
  }

  // Optionally allow pointers for Pomerons in the proton.
  if (pdfPomAPtrIn != 0 && pdfPomBPtrIn != 0) {
    if (pdfPomAPtrIn == pdfPomBPtrIn) return false;
    pdfPomAPtr  = pdfPomAPtrIn;
    pdfPomBPtr  = pdfPomBPtrIn;
  }

  // Optionally allow pointers for Gammas in the leptons.
  if (pdfGamAPtrIn != 0 && pdfGamBPtrIn != 0) {
    if (pdfGamAPtrIn == pdfGamBPtrIn) return false;
    pdfGamAPtr  = pdfGamAPtrIn;
    pdfGamBPtr  = pdfGamBPtrIn;
  }

  // Optionally allow pointers for Hard PDFs for photons in the leptons.
  if (pdfHardGamAPtrIn != 0 && pdfHardGamBPtrIn != 0) {
    if (pdfHardGamAPtrIn == pdfHardGamBPtrIn) return false;
    pdfHardGamAPtr  = pdfHardGamAPtrIn;
    pdfHardGamBPtr  = pdfHardGamBPtrIn;
  }

  // Optionally allow pointers for unresolved PDFs.
  if (pdfUnresAPtrIn != 0 && pdfUnresBPtrIn != 0) {
    if (pdfUnresAPtrIn == pdfUnresBPtrIn) return false;
    pdfUnresAPtr = pdfUnresAPtrIn;
    pdfUnresBPtr = pdfUnresBPtrIn;
  }

  // Optionally allow pointers for unresolved PDFs for photons from leptons.
  if (pdfUnresGamAPtrIn != 0 && pdfUnresGamBPtrIn != 0) {
    if (pdfUnresGamAPtrIn == pdfUnresGamBPtrIn) return false;
    pdfUnresGamAPtr = pdfUnresGamAPtrIn;
    pdfUnresGamBPtr = pdfUnresGamBPtrIn;
  }

  // Optionally allow pointers for VMD in the gamma.
  if (pdfVMDAPtrIn != 0 && pdfVMDBPtrIn != 0) {
    if (pdfVMDAPtrIn == pdfVMDBPtrIn) return false;
    pdfVMDAPtr  = pdfVMDAPtrIn;
    pdfVMDBPtr  = pdfVMDBPtrIn;
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Routine to pass in pointers to PDF's. Usage optional.

bool Pythia::setPDFAPtr( PDF* pdfAPtrIn ) {

  // Delete any PDF's created in a previous init call.
  if (useNewPdfHard && pdfHardAPtr != pdfAPtr) delete pdfHardAPtr;
  if (useNewPdfHard && pdfHardBPtr != pdfBPtr) delete pdfHardBPtr;
  if (useNewPdfA) delete pdfAPtr;
  if (useNewPdfB) delete pdfBPtr;
  if (useNewPdfPomA) delete pdfPomAPtr;
  if (useNewPdfPomB) delete pdfPomBPtr;
  if (useNewPdfGamA) delete pdfGamAPtr;
  if (useNewPdfGamB) delete pdfGamBPtr;
  if (useNewPdfUnresA) delete pdfUnresAPtr;
  if (useNewPdfUnresB) delete pdfUnresBPtr;
  if (useNewPdfUnresGamA) delete pdfUnresGamAPtr;
  if (useNewPdfUnresGamB) delete pdfUnresGamBPtr;
  if (useNewPdfHardGamA && pdfHardGamAPtr != pdfGamAPtr) delete pdfHardGamAPtr;
  if (useNewPdfHardGamB && pdfHardGamBPtr != pdfGamBPtr) delete pdfHardGamBPtr;
  if (useNewPdfVMDA) delete pdfVMDAPtr;
  if (useNewPdfVMDB) delete pdfVMDBPtr;

  // Reset pointers to be empty.
  useNewPdfA         = false;
  useNewPdfB         = false;
  useNewPdfHard      = false;
  useNewPdfPomA      = false;
  useNewPdfPomB      = false;
  useNewPdfGamA      = false;
  useNewPdfGamB      = false;
  useNewPdfHardGamA  = false;
  useNewPdfHardGamB  = false;
  useNewPdfUnresA    = false;
  useNewPdfUnresB    = false;
  useNewPdfUnresGamA = false;
  useNewPdfUnresGamB = false;
  useNewPdfVMDA      = false;
  useNewPdfVMDB      = false;
  pdfAPtr            = 0;
  pdfBPtr            = 0;
  pdfHardAPtr        = 0;
  pdfHardBPtr        = 0;
  pdfPomAPtr         = 0;
  pdfPomBPtr         = 0;
  pdfGamAPtr         = 0;
  pdfGamBPtr         = 0;
  pdfHardGamAPtr     = 0;
  pdfHardGamBPtr     = 0;
  pdfUnresAPtr       = 0;
  pdfUnresBPtr       = 0;
  pdfUnresGamAPtr    = 0;
  pdfUnresGamBPtr    = 0;
  pdfVMDAPtr         = 0;
  pdfVMDBPtr         = 0;

  // Switch off external PDF's by zero as input.
  if (pdfAPtrIn == 0) return true;

  // Save pointers.
  pdfAPtr       = pdfAPtrIn;
  // By default same pointers for hard-process PDF's.
  pdfHardAPtr   = pdfAPtrIn;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Routine to pass in pointers to PDF's. Usage optional.

bool Pythia::setPDFBPtr( PDF* pdfBPtrIn ) {

  // Delete any PDF's created in a previous init call.
  if (useNewPdfHard && pdfHardAPtr != pdfAPtr) delete pdfHardAPtr;
  if (useNewPdfHard && pdfHardBPtr != pdfBPtr) delete pdfHardBPtr;
  if (useNewPdfA) delete pdfAPtr;
  if (useNewPdfB) delete pdfBPtr;
  if (useNewPdfPomA) delete pdfPomAPtr;
  if (useNewPdfPomB) delete pdfPomBPtr;
  if (useNewPdfGamA) delete pdfGamAPtr;
  if (useNewPdfGamB) delete pdfGamBPtr;
  if (useNewPdfUnresA) delete pdfUnresAPtr;
  if (useNewPdfUnresB) delete pdfUnresBPtr;
  if (useNewPdfUnresGamA) delete pdfUnresGamAPtr;
  if (useNewPdfUnresGamB) delete pdfUnresGamBPtr;
  if (useNewPdfHardGamA && pdfHardGamAPtr != pdfGamAPtr) delete pdfHardGamAPtr;
  if (useNewPdfHardGamB && pdfHardGamBPtr != pdfGamBPtr) delete pdfHardGamBPtr;
  if (useNewPdfVMDA) delete pdfVMDAPtr;
  if (useNewPdfVMDB) delete pdfVMDBPtr;

  // Reset pointers to be empty.
  useNewPdfA         = false;
  useNewPdfB         = false;
  useNewPdfHard      = false;
  useNewPdfPomA      = false;
  useNewPdfPomB      = false;
  useNewPdfGamA      = false;
  useNewPdfGamB      = false;
  useNewPdfHardGamA  = false;
  useNewPdfHardGamB  = false;
  useNewPdfUnresA    = false;
  useNewPdfUnresB    = false;
  useNewPdfUnresGamA = false;
  useNewPdfUnresGamB = false;
  useNewPdfVMDA      = false;
  useNewPdfVMDB      = false;
  pdfAPtr            = 0;
  pdfBPtr            = 0;
  pdfHardAPtr        = 0;
  pdfHardBPtr        = 0;
  pdfPomAPtr         = 0;
  pdfPomBPtr         = 0;
  pdfGamAPtr         = 0;
  pdfGamBPtr         = 0;
  pdfHardGamAPtr     = 0;
  pdfHardGamBPtr     = 0;
  pdfUnresAPtr       = 0;
  pdfUnresBPtr       = 0;
  pdfUnresGamAPtr    = 0;
  pdfUnresGamBPtr    = 0;
  pdfVMDAPtr         = 0;
  pdfVMDBPtr         = 0;

  // Switch off external PDF's by zero as input.
  if (pdfBPtrIn == 0) return true;

  // Save pointers.
  pdfBPtr       = pdfBPtrIn;
  // By default same pointers for hard-process PDF's.
  pdfHardBPtr   = pdfBPtrIn;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Routine to initialize with the variable values of the Beams kind.

bool Pythia::init() {

  // Check that constructor worked.
  isInit = false;
  if (!isConstructed) {
    info.errorMsg("Abort from Pythia::init: constructor "
      "initialization failed");
    return false;
  }

  // Early catching of heavy ion mode.
  doHeavyIons = HeavyIons::isHeavyIon(settings) ||
                settings.mode("HeavyIon:mode") == 2;
  if ( doHeavyIons ) {
    if ( !heavyIonsPtr ) {
      heavyIonsPtr = new Angantyr(*this);
      hasOwnHeavyIons = true;
    }
    if ( hiHooksPtr ) heavyIonsPtr->setHIUserHooksPtr(hiHooksPtr);
    if ( !heavyIonsPtr->init() ) doHeavyIons = false;
  }

  // Early readout, if return false or changed when no beams.
  doProcessLevel = settings.flag("ProcessLevel:all");

  // Check that changes in Settings and ParticleData have worked.
  if (settings.unfinishedInput()) {
    info.errorMsg("Abort from Pythia::init: opening { not matched by "
      "closing }");
    return false;
  }
  if (settings.readingFailed()) {
    info.errorMsg("Abort from Pythia::init: some user settings "
      "did not make sense");
    return false;
  }
  if (particleData.readingFailed()) {
    info.errorMsg("Abort from Pythia::init: some user particle data "
      "did not make sense");
    return false;
  }

  // Initialize the random number generator.
  if ( settings.flag("Random:setSeed") )
    rndm.init( settings.mode("Random:seed") );

  // Find which frame type to use.
  info.addCounter(1);
  frameType = mode("Beams:frameType");

  // Set up values related to CKKW-L merging.
  bool doUserMerging     = settings.flag("Merging:doUserMerging");
  bool doMGMerging       = settings.flag("Merging:doMGMerging");
  bool doKTMerging       = settings.flag("Merging:doKTMerging");
  bool doPTLundMerging   = settings.flag("Merging:doPTLundMerging");
  bool doCutBasedMerging = settings.flag("Merging:doCutBasedMerging");
  // Set up values related to unitarised CKKW merging
  bool doUMEPSTree       = settings.flag("Merging:doUMEPSTree");
  bool doUMEPSSubt       = settings.flag("Merging:doUMEPSSubt");
  // Set up values related to NL3 NLO merging
  bool doNL3Tree         = settings.flag("Merging:doNL3Tree");
  bool doNL3Loop         = settings.flag("Merging:doNL3Loop");
  bool doNL3Subt         = settings.flag("Merging:doNL3Subt");
  // Set up values related to unitarised NLO merging
  bool doUNLOPSTree      = settings.flag("Merging:doUNLOPSTree");
  bool doUNLOPSLoop      = settings.flag("Merging:doUNLOPSLoop");
  bool doUNLOPSSubt      = settings.flag("Merging:doUNLOPSSubt");
  bool doUNLOPSSubtNLO   = settings.flag("Merging:doUNLOPSSubtNLO");
  bool doXSectionEst     = settings.flag("Merging:doXSectionEstimate");
  doMerging = doUserMerging || doMGMerging || doKTMerging
    || doPTLundMerging || doCutBasedMerging || doUMEPSTree || doUMEPSSubt
    || doNL3Tree || doNL3Loop || doNL3Subt || doUNLOPSTree
    || doUNLOPSLoop || doUNLOPSSubt || doUNLOPSSubtNLO || doXSectionEst;
  doMerging = doMerging || settings.flag("Merging:doMerging");

  // Set up MergingHooks object.
  bool inputMergingHooks = (mergingHooksPtr != 0);
  if (doMerging && !inputMergingHooks) {
    if (hasOwnMergingHooks && mergingHooksPtr) delete mergingHooksPtr;
    mergingHooksPtr = new MergingHooks();
    hasOwnMergingHooks = true;
  }

  hasMergingHooks  = (mergingHooksPtr != 0);
  // Merging hooks required for merging. If no merging hooks pointer is
  // available, exit.
  if (doMerging && !hasMergingHooks) {
    info.errorMsg("Abort from Pythia::init: "
      "no merging hooks object has been provided");
    return false;
  } else if (doMerging) {
    mergingHooksPtr->setLHEInputFile("");
  }

  // Set up Merging object.
  bool inputMerging = (mergingPtr != 0);
  if (doMerging && !inputMerging) {
    if (hasOwnMerging && mergingPtr) delete mergingPtr;
    mergingPtr = new Merging();
    hasOwnMerging = true;
  }
  hasMerging  = (mergingPtr != 0);

  // Initialization with internal processes: read in and set values.
  if (frameType < 4 ) {
    doLHA     = false;
    boostType = frameType;
    idA       = mode("Beams:idA");
    idB       = mode("Beams:idB");
    eCM       = parm("Beams:eCM");
    eA        = parm("Beams:eA");
    eB        = parm("Beams:eB");
    pxA       = parm("Beams:pxA");
    pyA       = parm("Beams:pyA");
    pzA       = parm("Beams:pzA");
    pxB       = parm("Beams:pxB");
    pyB       = parm("Beams:pyB");
    pzB       = parm("Beams:pzB");

   // Initialization with a Les Houches Event File or an LHAup object.
  } else {
    doLHA     = true;
    boostType = 2;
    string lhef        = word("Beams:LHEF");
    string lhefHeader  = word("Beams:LHEFheader");
    bool   readHeaders = flag("Beams:readLHEFheaders");
    bool   setScales   = flag("Beams:setProductionScalesFromLHEF");
    bool   skipInit    = flag("Beams:newLHEFsameInit");
    int    nSkipAtInit = mode("Beams:nSkipLHEFatInit");

    // For file input: renew file stream or (re)new Les Houches object.
    if (frameType == 4) {
      const char* cstring1 = lhef.c_str();
      bool useExternal = (lhaUpPtr && !useNewLHA && lhaUpPtr->useExternal());
      if (!useExternal && useNewLHA && skipInit)
        lhaUpPtr->newEventFile(cstring1);
      else if (!useExternal) {
        if (useNewLHA) delete lhaUpPtr;
        // Header is optional, so use NULL pointer to indicate no value.
        const char* cstring2 = (lhefHeader == "void")
          ? NULL : lhefHeader.c_str();
        lhaUpPtr   = new LHAupLHEF(&info, cstring1, cstring2, readHeaders,
                         setScales);
        useNewLHA  = true;
      }

      // Check that file was properly opened.
      if (!lhaUpPtr->fileFound()) {
        info.errorMsg("Abort from Pythia::init: "
          "Les Houches Event File not found");
        return false;
      }

    // For object input: at least check that not null pointer.
    } else {
      if (lhaUpPtr == 0) {
        info.errorMsg("Abort from Pythia::init: "
          "LHAup object not found");
        return false;
      }

      // LHAup object generic abort using fileFound() routine.
      if (!lhaUpPtr->fileFound()) {
        info.errorMsg("Abort from Pythia::init: "
          "LHAup initialisation error");
        return false;
      }
    }

    // Send in pointer to info. Store or replace LHA pointer in other classes.
    lhaUpPtr->setPtr( &info);
    processLevel.setLHAPtr( lhaUpPtr);

    // If second time around, only with new file, then simplify.
    // Optionally skip ahead a number of events at beginning of file.
    if (skipInit) {
      isInit = true;
      info.addCounter(2);
      if (nSkipAtInit > 0) lhaUpPtr->skipEvent(nSkipAtInit);
      return true;
    }

    // Set up values related to merging hooks.
    if (frameType == 4 || frameType == 5) {

      // Store the name of the input LHEF for merging.
      if (doMerging) {
        string lhefIn = (frameType == 4) ? lhef : "";
        mergingHooksPtr->setLHEInputFile( lhefIn);
      }

    }

    // Set LHAinit information (in some external program).
    if ( !lhaUpPtr->setInit()) {
      info.errorMsg("Abort from Pythia::init: "
        "Les Houches initialization failed");
      return false;
    }

    // Extract beams from values set in an LHAinit object.
    idA = lhaUpPtr->idBeamA();
    idB = lhaUpPtr->idBeamB();
    int idRenameBeams = settings.mode("LesHouches:idRenameBeams");
    if (abs(idA) == idRenameBeams) idA = 16;
    if (abs(idB) == idRenameBeams) idB = -16;
    if (idA == 0 || idB == 0) doProcessLevel = false;
    eA  = lhaUpPtr->eBeamA();
    eB  = lhaUpPtr->eBeamB();

    // Optionally skip ahead a number of events at beginning of file.
    if (nSkipAtInit > 0) lhaUpPtr->skipEvent(nSkipAtInit);
  }

  // Set up values related to user hooks.
  hasUserHooks     = (userHooksPtr != 0);
  doVetoProcess    = false;
  doVetoPartons    = false;
  retryPartonLevel = false;
  if (hasUserHooks) {
    userHooksPtr->initPtr( &info, &settings, &particleData, &rndm, &beamA,
      &beamB, &beamPomA, &beamPomB, couplingsPtr, &partonSystems, &sigmaTot);
    if (!userHooksPtr->initAfterBeams()) {
      info.errorMsg("Abort from Pythia::init: could not initialise UserHooks");
      return false;
    }
    doVetoProcess    = userHooksPtr->canVetoProcessLevel();
    doVetoPartons    = userHooksPtr->canVetoPartonLevel();
    retryPartonLevel = userHooksPtr->retryPartonLevel();
  }

  // Back to common initialization. Reset error counters.
  nErrEvent = 0;
  info.setTooLowPTmin(false);
  info.sigmaReset();

  // Initialize data members extracted from database.
  doPartonLevel    = settings.flag("PartonLevel:all");
  doHadronLevel    = settings.flag("HadronLevel:all");
  doCentralDiff    = settings.flag("SoftQCD:centralDiffractive");
  doSoftQCDall     = settings.flag("SoftQCD:all");
  doSoftQCDinel    = settings.flag("SoftQCD:inelastic");
  doDiffraction    = settings.flag("SoftQCD:singleDiffractive")
                  || settings.flag("SoftQCD:doubleDiffractive")
                  || doSoftQCDall || doSoftQCDinel || doCentralDiff;
  doSoftQCD        = doDiffraction ||
                     settings.flag("SoftQCD:nonDiffractive") ||
                     settings.flag("SoftQCD:elastic");
  doHardDiff       = settings.flag("Diffraction:doHard");
  doResDec         = settings.flag("ProcessLevel:resonanceDecays");
  doFSRinRes       = doPartonLevel && settings.flag("PartonLevel:FSR")
                  && settings.flag("PartonLevel:FSRinResonances");
  decayRHadrons    = settings.flag("RHadrons:allowDecay");
  doMomentumSpread = settings.flag("Beams:allowMomentumSpread");
  doVertexSpread   = settings.flag("Beams:allowVertexSpread");
  doVarEcm         = settings.flag("Beams:allowVariableEnergy");
  eMinPert         = settings.parm("Beams:eMinPert");
  eWidthPert       = settings.parm("Beams:eWidthPert");
  abortIfVeto      = settings.flag("Check:abortIfVeto");
  checkEvent       = settings.flag("Check:event");
  checkHistory     = settings.flag("Check:history");
  nErrList         = settings.mode("Check:nErrList");
  epTolErr         = settings.parm("Check:epTolErr");
  epTolWarn        = settings.parm("Check:epTolWarn");
  mTolErr          = settings.parm("Check:mTolErr");
  mTolWarn         = settings.parm("Check:mTolWarn");

  // Warn/abort for unallowed process and beam combinations.
  bool doHardProc  = !settings.onlySoftQCD() || doLHA;
  if (doSoftQCD && doHardProc) {
    info.errorMsg("Warning from Pythia::init: "
      "should not combine softQCD processes with hard ones");
  }
  if (doVarEcm) doMomentumSpread = false;
  if (doVarEcm && doHardProc) {
    info.errorMsg("Abort from Pythia::init: "
      "variable energy only works for softQCD processes");
    return false;
  }

  // Find out if beam is or has a resolved photon beam.
  beamHasGamma     = settings.flag("PDF:lepton2gamma");
  gammaMode        = settings.mode("Photon:ProcessType");
  bool beamAneedResGamma = (gammaMode == 1) || (gammaMode == 2)
    || (gammaMode == 0);
  bool beamBneedResGamma = (gammaMode == 1) || (gammaMode == 3)
    || (gammaMode == 0);
  beamAisResGamma  = beamAneedResGamma && idA == 22;
  beamBisResGamma  = beamBneedResGamma && idB == 22;
  bool isChargedLeptonA = (abs(idA) == 11 || abs(idA) == 13 || abs(idA) == 15);
  bool isChargedLeptonB = (abs(idB) == 11 || abs(idB) == 13 || abs(idB) == 15);
  beamAhasResGamma = beamAneedResGamma && beamHasGamma && isChargedLeptonA;
  beamBhasResGamma = beamBneedResGamma && beamHasGamma && isChargedLeptonB;
  doVMDsideA = doSoftQCD && (beamAisResGamma || beamAhasResGamma);
  doVMDsideB = doSoftQCD && (beamBisResGamma || beamBhasResGamma);

  // Turn off central diffraction for VMD processes.
  if (doVMDsideA || doVMDsideB) {
    if (doCentralDiff) {
      info.errorMsg("Warning in Pythia::init: "
        "Central diffractive events not implemented for gamma + p/gamma");
      return false;
    }
    if (doSoftQCDall) {
      info.errorMsg("Warning in Pythia::init: "
        "Central diffractive events not implemented for gamma + p/gamma");
      settings.flag("SoftQCD:all", false);
      settings.flag("SoftQCD:elastic", true);
      settings.flag("SoftQCD:nonDiffractive", true);
      settings.flag("SoftQCD:singleDiffractive", true);
      settings.flag("SoftQCD:doubleDiffractive", true);
    }
    if (doSoftQCDinel) {
      info.errorMsg("Warning in Pythia::init: "
        "Central diffractive events not implemented for gamma + p/gamma");
      settings.flag("SoftQCD:inelastic", false);
      settings.flag("SoftQCD:nonDiffractive", true);
      settings.flag("SoftQCD:singleDiffractive", true);
      settings.flag("SoftQCD:doubleDiffractive", true);
    }
  }

  // Initialise merging hooks.
  if ( doMerging && (hasMergingHooks || hasOwnMergingHooks) ) {
    mergingHooksPtr->initPtr( &settings, &info, &particleData, &partonSystems);
    mergingHooksPtr->init();
  }

  // Check that combinations of settings are allowed; change if not.
  checkSettings();

  // Initialize the SM couplings (needed to initialize resonances).
  couplingsPtr->init( settings, &rndm );

  // Initialize SLHA interface (including SLHA/BSM couplings).
  bool useSLHAcouplings = false;
  slhaInterface = SLHAinterface();
  slhaInterface.setPtr( &info);
  slhaInterface.init( settings, &rndm, couplingsPtr, &particleData,
    useSLHAcouplings, particleDataBuffer );
  if (useSLHAcouplings) couplingsPtr = slhaInterface.couplingsPtr;

  // Reset couplingsPtr to the correct memory address.
  particleData.initPtr( &info, &settings, &rndm, couplingsPtr);
  if (hasUserHooks) userHooksPtr->initPtr( &info, &settings, &particleData,
    &rndm, &beamA, &beamB, &beamPomA, &beamPomB, couplingsPtr,
    &partonSystems, &sigmaTot);

  // Set headers to distinguish the two event listing kinds.
  int startColTag = settings.mode("Event:startColTag");
  process.init("(hard process)", &particleData, startColTag);
  event.init("(complete event)", &particleData, startColTag);

  // Final setup stage of particle data, notably resonance widths.
  particleData.initWidths( resonancePtrs);

  // Set up R-hadrons particle data, where relevant.
  rHadrons.init( &info, settings, &particleData, &rndm);

  // Set up and initialize setting of parton production vertices.
  if (settings.flag("PartonVertex:setVertex")) {
    if (partonVertexPtr == 0) {
      partonVertexPtr    = new PartonVertex();
      useNewPartonVertex = true;
    }
    partonVertexPtr->initPtr( &info, &settings, &rndm);
    partonVertexPtr->init();
  }

  // Set up objects for timelike and spacelike showers.
  if (timesDecPtr == 0 || timesPtr == 0) {
    TimeShower* timesNow = new SimpleTimeShower();
    if (timesDecPtr == 0) {
      timesDecPtr    = timesNow;
      useNewTimesDec = true;
    }
    if (timesPtr == 0) {
      timesPtr    = timesNow;
      useNewTimes = true;
    }
  }
  if (spacePtr == 0) {
    spacePtr    = new SimpleSpaceShower();
    useNewSpace = true;
  }

  // Initialize pointers in showers.
  timesPtr->initPtr( &info, &settings, &particleData, &rndm, couplingsPtr,
    &partonSystems, userHooksPtr, mergingHooksPtr, partonVertexPtr);
  timesDecPtr->initPtr( &info, &settings, &particleData, &rndm, couplingsPtr,
    &partonSystems, userHooksPtr, mergingHooksPtr, partonVertexPtr);
  spacePtr->initPtr( &info, &settings, &particleData, &rndm, couplingsPtr,
    &partonSystems, userHooksPtr, mergingHooksPtr, partonVertexPtr);

  // Set up values related to beam shape.
  if (beamShapePtr == 0) {
    beamShapePtr    = new BeamShape();
    useNewBeamShape = true;
  }
  beamShapePtr->init( settings, &rndm);

  // Check that beams and beam combination can be handled.
  if (!checkBeams()) {
    info.errorMsg("Abort from Pythia::init: "
      "checkBeams initialization failed");
    return false;
  }

  // Do not set up beam kinematics when no process level.
  if (!doProcessLevel) boostType = 1;
  else {

    // Set up beam kinematics.
    if (!initKinematics()) {
      info.errorMsg("Abort from Pythia::init: "
        "kinematics initialization failed");
      return false;
    }

    // Set up pointers to PDFs.
    if (!initPDFs()) {
      info.errorMsg("Abort from Pythia::init: PDF initialization failed");
      return false;
    }

    // Set up the two beams and the common remnant system.
    StringFlav* flavSelPtr = hadronLevel.getStringFlavPtr();
    beamA.init( idA, pzAcm, eA, mA, &info, settings, &particleData, &rndm,
      pdfAPtr, pdfHardAPtr, isUnresolvedA, flavSelPtr);
    beamB.init( idB, pzBcm, eB, mB, &info, settings, &particleData, &rndm,
      pdfBPtr, pdfHardBPtr, isUnresolvedB, flavSelPtr);

    // Init also unresolved PDF pointers for photon beams when needed.
    if ( ( beamA.isGamma() || beamAhasResGamma )
        && ( gammaMode == 0 || gammaMode == 3 || gammaMode == 4 ) )
      beamA.initUnres( pdfUnresAPtr);
    if ( ( beamB.isGamma() || beamBhasResGamma )
        && ( gammaMode == 0 || gammaMode == 2 || gammaMode == 4 ) )
      beamB.initUnres( pdfUnresBPtr);

    // Optionally set up new alternative beams for these Pomerons.
    if ( doDiffraction || doHardDiff ) {
      beamPomA.init( 990,  0.5 * eCM, 0.5 * eCM, 0., &info, settings,
        &particleData, &rndm, pdfPomAPtr, pdfPomAPtr, false, flavSelPtr);
      beamPomB.init( 990, -0.5 * eCM, 0.5 * eCM, 0., &info, settings,
        &particleData, &rndm, pdfPomBPtr, pdfPomBPtr, false, flavSelPtr);
    }

    // Initialise VMD beams from gammas (in leptons). Use pion PDF for VMDs.
    if (doVMDsideA)
      beamVMDA.init( 111,  0.5 * eCM, 0.5 * eCM, 0., &info, settings,
        &particleData, &rndm, pdfVMDAPtr, pdfVMDAPtr, false, flavSelPtr);
    if (doVMDsideB)
      beamVMDB.init( 111,  0.5 * eCM, 0.5 * eCM, 0., &info, settings,
        &particleData, &rndm, pdfVMDBPtr, pdfVMDBPtr, false, flavSelPtr);

    // Optionally set up photon beams from lepton beams if resolved photons.
    if (beamHasGamma && gammaMode < 4) {
      if ( !(beamA.isHadron()) )
        beamGamA.init( 22,  0.5 * eCM, 0.5 * eCM, 0., &info, settings,
          &particleData, &rndm, pdfGamAPtr, pdfHardGamAPtr, !beamAisResGamma,
          flavSelPtr);
      if ( !(beamB.isHadron()) )
        beamGamB.init( 22,  0.5 * eCM, 0.5 * eCM, 0., &info, settings,
          &particleData, &rndm, pdfGamBPtr, pdfHardGamBPtr, !beamBisResGamma,
          flavSelPtr);

      // Initialize also unresolved PDFs for relevant processes.
      if ( gammaMode == 0 || gammaMode == 3 )
        beamGamA.initUnres( pdfUnresGamAPtr);
      if ( gammaMode == 0 || gammaMode == 2 )
        beamGamB.initUnres( pdfUnresGamBPtr);
    }

  }

  // Send info/pointers to process level for initialization.
  if ( doProcessLevel && !processLevel.init( &info, settings, &particleData,
    &rndm, &beamA, &beamB, &beamGamA, &beamGamB, &beamVMDA, &beamVMDB,
    couplingsPtr, &sigmaTot, doLHA, &slhaInterface, userHooksPtr,
    sigmaPtrs, phaseSpacePtrs) ) {
    info.errorMsg("Abort from Pythia::init: "
      "processLevel initialization failed");
    return false;
  }

  // Initialize timelike showers already here, since needed in decays.
  // The pointers to the beams are needed by some external plugin showers.
  timesDecPtr->init( &beamA, &beamB);

  // Alternatively only initialize resonance decays.
  if ( !doProcessLevel) processLevel.initDecays( &info, settings,
    &particleData, &rndm, lhaUpPtr);

  // Send info/pointers to parton level for initialization.
  if ( doPartonLevel && doProcessLevel && !partonLevel.init( &info, settings,
    &particleData, &rndm, &beamA, &beamB, &beamPomA, &beamPomB, &beamGamA,
    &beamGamB, &beamVMDA, &beamVMDB, couplingsPtr, &partonSystems, &sigmaTot,
    timesDecPtr, timesPtr, spacePtr, &rHadrons, userHooksPtr, mergingHooksPtr,
    partonVertexPtr, false) ) {
    info.errorMsg("Abort from Pythia::init: "
      "partonLevel initialization failed" );
    return false;
  }

  // Make pointer to shower available for merging machinery.
  if ( doMerging && (hasMergingHooks || hasOwnMergingHooks) )
    mergingHooksPtr->setShowerPointer(&partonLevel);

  // Alternatively only initialize final-state showers in resonance decays.
  if ( !doProcessLevel || !doPartonLevel) partonLevel.init( &info, settings,
    &particleData, &rndm, 0, 0, 0, 0, 0, 0, 0, 0, couplingsPtr, &partonSystems,
    0, timesDecPtr, 0, 0, &rHadrons, 0, 0, partonVertexPtr, false);

  // Send info/pointers to parton level for trial shower initialization.
  if ( doMerging && !trialPartonLevel.init( &info, settings, &particleData,
      &rndm, &beamA, &beamB, &beamPomA, &beamPomB, &beamGamA, &beamGamB,
      &beamVMDA, &beamVMDB, couplingsPtr, &partonSystems, &sigmaTot,
      timesDecPtr, timesPtr, spacePtr, &rHadrons, userHooksPtr,
      mergingHooksPtr, partonVertexPtr, true) ) {
    info.errorMsg("Abort from Pythia::init: "
      "trialPartonLevel initialization failed");
    return false;
  }

  // Initialise the merging wrapper class.
  if (doMerging ) {
    mergingPtr->initPtr( &settings, &info, &particleData, &rndm,
      &beamA, &beamB, mergingHooksPtr, &trialPartonLevel, couplingsPtr );
    mergingPtr->init();
  }

  // Send info/pointers to hadron level for initialization.
  // Note: forceHadronLevel() can come, so we must always initialize.
  if ( !hadronLevel.init( &info, settings, &particleData, &rndm,
    couplingsPtr, timesDecPtr, &rHadrons, decayHandlePtr,
    handledParticles, userHooksPtr) ) {
    info.errorMsg("Abort from Pythia::init: "
      "hadronLevel initialization failed");
    return false;
  }

  // Optionally check particle data table for inconsistencies.
  if ( settings.flag("Check:particleData") )
    particleData.checkTable( settings.mode("Check:levelParticleData") );

  // Optionally show settings and particle data, changed or all.
  bool showCS  = settings.flag("Init:showChangedSettings");
  bool showAS  = settings.flag("Init:showAllSettings");
  bool showCPD = settings.flag("Init:showChangedParticleData");
  bool showCRD = settings.flag("Init:showChangedResonanceData");
  bool showAPD = settings.flag("Init:showAllParticleData");
  int  show1PD = settings.mode("Init:showOneParticleData");
  bool showPro = settings.flag("Init:showProcesses");
  if (showCS)      settings.listChanged();
  if (showAS)      settings.listAll();
  if (show1PD > 0) particleData.list(show1PD);
  if (showCPD)     particleData.listChanged(showCRD);
  if (showAPD)     particleData.listAll();

  // Listing options for the next() routine.
  nCount       = settings.mode("Next:numberCount");
  nShowLHA     = settings.mode("Next:numberShowLHA");
  nShowInfo    = settings.mode("Next:numberShowInfo");
  nShowProc    = settings.mode("Next:numberShowProcess");
  nShowEvt     = settings.mode("Next:numberShowEvent");
  showSaV      = settings.flag("Next:showScaleAndVertex");
  showMaD      = settings.flag("Next:showMothersAndDaughters");

  // Init colour reconnection and junction splitting.
  colourReconnection.init( &info, settings, &rndm, &particleData,
    &beamA, &beamB, &partonSystems);
  junctionSplitting.init(&info, settings, &rndm, &particleData);

  // Flags for colour reconnection.
  doReconnect        = settings.flag("ColourReconnection:reconnect");
  reconnectMode      = settings.mode("ColourReconnection:mode");
  forceHadronLevelCR = settings.flag("ColourReconnection:forceHadronLevelCR");

  // Succeeded.
  isInit = true;
  info.addCounter(2);
  if (useNewLHA && showPro) lhaUpPtr->listInit();
  return true;

}

//--------------------------------------------------------------------------

// Check that combinations of settings are allowed; change if not.

void Pythia::checkSettings() {

  // Double rescattering not allowed if ISR or FSR.
  if ((settings.flag("PartonLevel:ISR") || settings.flag("PartonLevel:FSR"))
    && settings.flag("MultipartonInteractions:allowDoubleRescatter")) {
    info.errorMsg("Warning in Pythia::checkSettings: "
        "double rescattering switched off since showering is on");
    settings.flag("MultipartonInteractions:allowDoubleRescatter", false);
  }

  // Collisions with direct photon(s).
  if ( (idA == 22 && !beamAisResGamma) || (idB == 22 && !beamBisResGamma) ) {

    // Turn MPIs off.
    if ( settings.flag("PartonLevel:MPI") ) {
      info.errorMsg("Warning in Pythia::checkSettings: "
        "MPIs turned off for collision with unresolved photon");
      settings.flag("PartonLevel:MPI", false);
    }

    // Check that no soft QCD processes initialized with direct photons.
    if ( settings.flag("SoftQCD:nonDiffractive") ) {
      info.errorMsg("Warning in Pythia::checkSettings: "
        "Soft QCD processes turned off for collision with unresolved photon");
      settings.flag("SoftQCD:nonDiffractive", false);
    }

  }

  // Lepton-lepton/hadron collisions with direct photon(s).
  if ( ( (abs(idA) > 10 && abs(idA) < 17)
    && !beamAhasResGamma && beamHasGamma)
    || ( (abs(idB) > 10 && abs(idB) < 17)
    && !beamBhasResGamma && beamHasGamma) ) {

    // Turn MPIs off.
    if ( settings.flag("PartonLevel:MPI") ) {
      info.errorMsg("Warning in Pythia::checkSettings: MPIs turned off for "
                    "collision with unresolved photon");
      settings.flag("PartonLevel:MPI", false);
    }

    // Check that no soft QCD processes initialized with direct photons.
    if ( settings.flag("SoftQCD:nonDiffractive") ) {
      info.errorMsg("Warning in Pythia::checkSettings: "
        "Soft QCD processes turned off for collision with unresolved photon");
      settings.flag("SoftQCD:nonDiffractive", false);
    }

  }

}

//--------------------------------------------------------------------------

// Check that beams and beam combination can be handled. Set up unresolved.

bool Pythia::checkBeams() {

  // Absolute flavours. If not to do process level then no check needed.
  int idAabs = abs(idA);
  int idBabs = abs(idB);
  if (!doProcessLevel) return true;

  // Neutrino beams always unresolved, charged lepton ones conditionally.
  bool isLeptonA    = (idAabs > 10 && idAabs < 17);
  bool isLeptonB    = (idBabs > 10 && idBabs < 17);
  bool isUnresLep   = !settings.flag("PDF:lepton");
  bool isGammaA     = idAabs == 22;
  bool isGammaB     = idBabs == 22;
  isUnresolvedA     = (isLeptonA && isUnresLep);
  isUnresolvedB     = (isLeptonB && isUnresLep);

  // Also photons may be unresolved.
  if ( idAabs == 22 && !beamAisResGamma ) isUnresolvedA = true;
  if ( idBabs == 22 && !beamBisResGamma ) isUnresolvedB = true;

  // If photons from beam particle, beam not unresolved.
  if ( beamAhasResGamma ) isUnresolvedA = false;
  if ( beamBhasResGamma ) isUnresolvedB = false;

  // Equate Dark Matter "beams" with incoming neutrinos.
  if (idAabs > 50 && idAabs < 61) isLeptonA = isUnresolvedA = true;
  if (idBabs > 50 && idBabs < 61) isLeptonB = isUnresolvedB = true;

  // Lepton-lepton collisions.
  if (isLeptonA && isLeptonB ) {

    // Photon-photon collision from lepton beams.
    if (beamHasGamma) {

      // Non-diffractive events only for resolved photon-photon.
      if ( (!beamAhasResGamma || !beamBhasResGamma)
        && settings.flag("SoftQCD:nonDiffractive") ) {
        info.errorMsg("Error in Pythia::init: Soft QCD only with resolved"
          " photons with lepton beams.");
        return false;

      }
      // Otherwise photon-photon within lepton beams OK.
      else return true;
    }

    // Lepton-lepton collisions OK (including neutrinos) if both (un)resolved
    else if (isUnresolvedA == isUnresolvedB) return true;
  }

  // MBR model only implemented for pp/ppbar/pbarp collisions.
  int PomFlux     = settings.mode("SigmaDiffractive:PomFlux");
  if (PomFlux == 5) {
    bool ispp       = (idAabs == 2212 && idBabs == 2212);
    bool ispbarpbar = (idA == -2212 && idB == -2212);
    if (ispp && !ispbarpbar) return true;
    info.errorMsg("Error in Pythia::init: cannot handle this beam combination"
      " with PomFlux == 5");
    return false;
  }

  // Hadron-hadron collisions OK, with Pomeron counted as hadron.
  bool isHadronA = (idAabs == 2212) || (idAabs == 2112) || (idA == 111)
                || (idAabs == 211)  || (idA == 990);
  bool isHadronB = (idBabs == 2212) || (idBabs == 2112) || (idB == 111)
                || (idBabs == 211)  || (idB == 990);
  int modeUnresolvedHadron = settings.mode("BeamRemnants:unresolvedHadron");
  if (isHadronA && modeUnresolvedHadron%2 == 1) isUnresolvedA = true;
  if (isHadronB && modeUnresolvedHadron > 1)    isUnresolvedB = true;
  if (isHadronA && isHadronB) {
    // lepton2gamma flag with hadron beams may cause problems.
    if (beamHasGamma) {
      info.errorMsg("Error in Pythia::init: lepton2gamma should be off for"
        " hadron+hadron collision");
      return false;
    } else {
      return true;
    }
  }

  // Photon-photon collisions.
  if ( (idAabs == 22) && (idBabs == 22) ) {

    // No non-diffractive events for unresolved photon-photon.
    if ( ( !beamAisResGamma || !beamBisResGamma )
         && settings.flag("SoftQCD:nonDiffractive") ) {
      info.errorMsg("Error in Pythia::init: Soft QCD only with resolved"
        " photons.");
    }

    // lepton2gamma flag with photon beams may cause problems.
    if (beamHasGamma) {
      info.errorMsg("Error in Pythia::init: lepton2gamma should be off for"
        " hadron+hadron collision");
      return false;
    }

    // Otherwise OK.
    else return true;
  }

  // Gamma+hadron mode OK.
  if ( (isGammaA && isHadronB) || (isGammaB && isHadronA) )
    return true;

  // Lepton-hadron collisions OK for DIS processes or LHEF input,
  // although still primitive. Also e+p with real photons.
  if ( (isLeptonA && isHadronB) || (isHadronA && isLeptonB) ) {
    bool doDIS = settings.flag("WeakBosonExchange:all")
              || settings.flag("WeakBosonExchange:ff2ff(t:gmZ)")
              || settings.flag("WeakBosonExchange:ff2ff(t:W)")
              || !settings.flag("Check:beams")
              || (frameType == 4);
    if (doDIS || beamHasGamma ) return true;
  }

  // Allow to explicitly omit beam check for LHEF input.
  if ( settings.mode("Beams:frameType") == 4
    && !settings.flag("Check:beams")) return true;

  // If no case above then failed.
  info.errorMsg("Error in Pythia::init: cannot handle this beam combination");
  return false;

}

//--------------------------------------------------------------------------

// Calculate kinematics at initialization. Store beam four-momenta.

bool Pythia::initKinematics() {

  // Find masses. Initial guess that we are in CM frame.
  mA       = particleData.m0(idA);
  mB       = particleData.m0(idB);
  betaZ    = 0.;
  gammaZ   = 1.;

  // Collinear beams not in CM frame: find CM energy.
  if (boostType == 2) {
    eA     = max(eA, mA);
    eB     = max(eB, mB);
    pzA    = sqrt(eA*eA - mA*mA);
    pzB    = -sqrt(eB*eB - mB*mB);
    pAinit = Vec4( 0., 0., pzA, eA);
    pBinit = Vec4( 0., 0., pzB, eB);
    eCM    = sqrt( pow2(eA + eB) - pow2(pzA + pzB) );

    // Find boost to rest frame.
    betaZ  = (pzA + pzB) / (eA + eB);
    gammaZ = (eA + eB) / eCM;
    if (abs(betaZ) < 1e-10) boostType = 1;
  }

  // Completely general beam directions: find CM energy.
  else if (boostType == 3) {
    eA     = sqrt( pxA*pxA + pyA*pyA + pzA*pzA + mA*mA);
    eB     = sqrt( pxB*pxB + pyB*pyB + pzB*pzB + mB*mB);
    pAinit = Vec4( pxA, pyA, pzA, eA);
    pBinit = Vec4( pxB, pyB, pzB, eB);
    eCM = (pAinit + pBinit).mCalc();

    // Find boost+rotation needed to move from/to CM frame.
    MfromCM.reset();
    MfromCM.fromCMframe( pAinit, pBinit);
    MtoCM = MfromCM;
    MtoCM.invert();
  }

  // Fail if CM energy below beam masses.
  if (eCM < mA + mB) {
    info.errorMsg("Error in Pythia::initKinematics: too low energy");
    return false;
  }

  // Set up CM-frame kinematics with beams along +-z axis.
  pzAcm    = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
           * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
  pzBcm    = -pzAcm;
  eA       = sqrt(mA*mA + pzAcm*pzAcm);
  eB       = sqrt(mB*mB + pzBcm*pzBcm);

  // If in CM frame then store beam four-vectors (else already done above).
  if (boostType != 2 && boostType != 3) {
    pAinit = Vec4( 0., 0., pzAcm, eA);
    pBinit = Vec4( 0., 0., pzBcm, eB);
  }

  // Store main info for access in process generation.
  info.setBeamA( idA, pzAcm, eA, mA);
  info.setBeamB( idB, pzBcm, eB, mB);
  info.setECM( eCM);

  // Must allow for generic boost+rotation when beam momentum spread.
  if (doMomentumSpread) boostType = 3;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Set up pointers to PDFs.

bool Pythia::initPDFs() {

  // Delete any PDF's created in a previous init call.
  if (useNewPdfHard) {
    if (pdfHardAPtr != pdfAPtr) {
      delete pdfHardAPtr;
      pdfHardAPtr = 0;
    }
    if (pdfHardBPtr != pdfBPtr) {
      delete pdfHardBPtr;
      pdfHardBPtr = 0;
    }
    useNewPdfHard = false;
  }
  if (useNewPdfA) {
    delete pdfAPtr;
    useNewPdfA    = false;
    pdfAPtr       = 0;
  }
  if (useNewPdfB) {
    delete pdfBPtr;
    useNewPdfB    = false;
    pdfBPtr       = 0;
  }
  if (useNewPdfPomA) {
    delete pdfPomAPtr;
    useNewPdfPomA = false;
    pdfPomAPtr    = 0;
  }
  if (useNewPdfPomB) {
    delete pdfPomBPtr;
    useNewPdfPomB = false;
    pdfPomBPtr    = 0;
  }
  if (useNewPdfGamA) {
    delete pdfGamAPtr;
    useNewPdfGamA = false;
    pdfGamAPtr    = 0;
  }
  if (useNewPdfGamB) {
    delete pdfGamBPtr;
    useNewPdfGamB = false;
    pdfGamBPtr    = 0;
  }
  if (useNewPdfHardGamA) {
    delete pdfHardGamAPtr;
    useNewPdfHardGamA = false;
    pdfHardGamAPtr    = 0;
  }
  if (useNewPdfHardGamB) {
    delete pdfHardGamBPtr;
    useNewPdfHardGamB = false;
    pdfHardGamBPtr    = 0;
  }
  if (useNewPdfUnresA) {
    delete pdfUnresAPtr;
    useNewPdfUnresA = false;
    pdfUnresAPtr    = 0;
  }
  if (useNewPdfUnresB) {
    delete pdfUnresBPtr;
    useNewPdfUnresB = false;
    pdfUnresBPtr    = 0;
  }
  if (useNewPdfUnresGamA) {
    delete pdfUnresGamAPtr;
    useNewPdfUnresGamA = false;
    pdfUnresGamAPtr    = 0;
  }
  if (useNewPdfUnresGamB) {
    delete pdfUnresGamBPtr;
    useNewPdfUnresGamB = false;
    pdfUnresGamBPtr    = 0;
  }
  if (useNewPdfVMDA) {
    delete pdfVMDAPtr;
    useNewPdfVMDA = false;
    pdfVMDAPtr    = 0;
  }
  if (useNewPdfVMDB) {
    delete pdfVMDBPtr;
    useNewPdfVMDB = false;
    pdfVMDBPtr    = 0;
  }

  // Optionally set up photon PDF's for lepton -> gamma collisions. Done before
  // the main PDFs so that the gamma pointer can be used for the main PDF
  // (lepton). Both set also in case that only one of the photons is resolved.
  bool setupGammaBeams = (settings.flag("PDF:lepton2gamma")
    && (gammaMode < 4) );
  if (setupGammaBeams) {
    if ( (abs(idA) == 11 || abs(idA) == 13 || abs(idA) == 15)
        && pdfGamAPtr == 0 ) {
      pdfGamAPtr = getPDFPtr(22, 1, "A");
      if (!pdfGamAPtr->isSetup()) return false;
      useNewPdfGamA = true;

      // Set also unresolved photon beam when also unresolved photons.
      if (gammaMode != 1) {
        pdfUnresGamAPtr = getPDFPtr(22, 1, "A", false);
        if (!pdfUnresGamAPtr->isSetup()) return false;
        useNewPdfUnresGamA = true;
      }

      // Set up optional hard photon PDF pointers.
      if (settings.flag("PDF:useHard")) {
        pdfHardGamAPtr = getPDFPtr(22, 2);
        if (!pdfHardGamAPtr->isSetup()) return false;
        useNewPdfHardGamA = true;
      } else pdfHardGamAPtr = pdfGamAPtr;
    }
    if ( (abs(idB) == 11 || abs(idB) == 13 || abs(idB) == 15)
        && pdfGamBPtr == 0 ) {
      pdfGamBPtr = getPDFPtr(22, 1, "B");
      if (!pdfGamBPtr->isSetup()) return false;
      useNewPdfGamB = true;

      // Set also unresolved photon beam when also unresolved photons.
      if (gammaMode != 1) {
        pdfUnresGamBPtr = getPDFPtr(22, 1, "B", false);
        if (!pdfUnresGamBPtr->isSetup()) return false;
        useNewPdfUnresGamB = true;
      }

      // Set up optional hard photon PDF pointers.
      if (settings.flag("PDF:useHard")) {
        pdfHardGamBPtr = getPDFPtr(22, 2, "B");
        if (!pdfHardGamBPtr->isSetup()) return false;
        useNewPdfHardGamB = true;
      } else pdfHardGamBPtr = pdfGamBPtr;
    }
  }

  // Set up the PDF's, if not already done.
  if (pdfAPtr == 0) {
    pdfAPtr     = getPDFPtr(idA);
    if (pdfAPtr == 0 || !pdfAPtr->isSetup()) {
      info.errorMsg("Error in Pythia::init: "
        "could not set up PDF for beam A");
      return false;
    }
    pdfHardAPtr = pdfAPtr;
    useNewPdfA  = true;
  }
  if (pdfBPtr == 0) {
    pdfBPtr     = getPDFPtr(idB, 1, "B");
    if (pdfBPtr == 0 || !pdfBPtr->isSetup()) {
      info.errorMsg("Error in Pythia::init: "
        "could not set up PDF for beam B");
      return false;
    }
    pdfHardBPtr = pdfBPtr;
    useNewPdfB  = true;
  }

  // Optionally set up separate PDF's for hard process.
  if (settings.flag("PDF:useHard") && useNewPdfA && useNewPdfB) {
    pdfHardAPtr = getPDFPtr(idA, 2);
    if (!pdfHardAPtr->isSetup()) return false;
    pdfHardBPtr = getPDFPtr(idB, 2, "B");
    if (!pdfHardBPtr->isSetup()) return false;
    useNewPdfHard = true;
  }

  // Optionally use nuclear modifications for hard process PDFs.
  if (settings.flag("PDF:useHardNPDFA")) {
    int idANucleus = settings.mode("PDF:nPDFBeamA");
    pdfHardAPtr = getPDFPtr(idANucleus, 2, "A");
    if (!pdfHardAPtr->isSetup()) {
      info.errorMsg("Error in Pythia::init: "
        "could not set up nuclear PDF for beam A");
      return false;
    }
    useNewPdfHard = true;
  }
  if (settings.flag("PDF:useHardNPDFB")) {
    int idBNucleus = settings.mode("PDF:nPDFBeamB");
    pdfHardBPtr = getPDFPtr(idBNucleus, 2, "B");
    if (!pdfHardBPtr->isSetup()) {
      info.errorMsg("Error in Pythia::init: "
        "could not set up nuclear PDF for beam B");
      return false;
    }
    useNewPdfHard = true;
  }

  // Optionally set up additional unresolved PDFs for photon beams.
  if ( (idA == 22 || idB == 22) && gammaMode != 1 ) {
    if ( idA == 22 && pdfUnresAPtr == 0 ) {
      pdfUnresAPtr = getPDFPtr(idA, 1, "A", false);
      if (!pdfUnresAPtr->isSetup()) return false;
      useNewPdfUnresA = true;
    }
    if ( idB == 22 && pdfUnresBPtr == 0 ) {
      pdfUnresBPtr = getPDFPtr(idB, 1, "B", false);
      if (!pdfUnresBPtr->isSetup()) return false;
      useNewPdfUnresB = true;
    }
  }

  // Optionally set up additional unresolved PDFs for photon beam from lepton.
  if ( (abs(idA) == 11 || abs(idA) == 13 || abs(idA) == 15)
        && beamHasGamma && gammaMode != 1 ) {
    if ( pdfUnresAPtr == 0 ) {
      pdfUnresAPtr = getPDFPtr(idA, 1, "A", false);
      if (!pdfUnresAPtr->isSetup()) return false;
      useNewPdfUnresA = true;
    }
  }
  if ( (abs(idB) == 11 || abs(idB) == 13 || abs(idB) == 15)
        && beamHasGamma && gammaMode != 1 ) {
    if ( pdfUnresBPtr == 0 ) {
      pdfUnresBPtr = getPDFPtr(idB, 1, "B", false);
      if (!pdfUnresBPtr->isSetup()) return false;
      useNewPdfUnresB = true;
    }
  }

  // Optionally set up Pomeron PDF's for diffractive physics.
  if ( doDiffraction || doHardDiff) {
    if (pdfPomAPtr == 0) {
      pdfPomAPtr    = getPDFPtr(990);
      useNewPdfPomA = true;
    }
    if (pdfPomBPtr == 0) {
      pdfPomBPtr    = getPDFPtr(990);
      useNewPdfPomB = true;
    }
  }

  // Optionally set up VMD PDF's for photon physics.
  if ( doSoftQCD && (doVMDsideA || doVMDsideB)) {
    if (pdfVMDAPtr == 0) {
      pdfVMDAPtr    = getPDFPtr(111);
      useNewPdfVMDA = true;
    }
    if (pdfVMDBPtr == 0) {
      pdfVMDBPtr    = getPDFPtr(111);
      useNewPdfVMDB = true;
    }
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Main routine to generate the next event, using internal machinery.

bool Pythia::next() {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Check if we the generation is taken over by the HeavyIons object.
  // Allows HeavyIons::next to call next for this Pythia object
  // without going into a loop.
  if ( doHeavyIons ) {
    doHeavyIons = false;
    bool ok = heavyIonsPtr->next();
    doHeavyIons = true;
    return ok;
  }

  // Regularly print how many events have been generated.
  int nPrevious = info.getCounter(3);
  if (nCount > 0 && nPrevious > 0 && nPrevious%nCount == 0)
    cout << "\n Pythia::next(): " << nPrevious
         << " events have been generated " << endl;

  // Set/reset info counters specific to each event.
  info.addCounter(3);
  for (int i = 10; i < 13; ++i) info.setCounter(i);

  // Simpler option when no hard process, i.e. mainly hadron level.
  if (!doProcessLevel) {

    // Optionally fetch in resonance decays from LHA interface.
    if (doLHA && !processLevel.nextLHAdec( event)) {
      if (info.atEndOfFile()) info.errorMsg("Abort from Pythia::next:"
        " reached end of Les Houches Events File");
      return false;
    }

    // Reset info and partonSystems arrays (while event record contains data).
    info.clear();
    partonSystems.clear();

    // Set correct energy for system.
    Vec4 pSum = 0.;
    for (int i = 1; i < event.size(); ++i)
      if (event[i].isFinal()) pSum += event[i].p();
    event[0].p( pSum );
    event[0].m( pSum.mCalc() );

    // Generate hadronization and decays.
    bool status = (doHadronLevel) ? forceHadronLevel() : true;
    if (status) info.addCounter(4);
    if (status && nPrevious < nShowEvt) event.list(showSaV, showMaD);
    return status;
  }

  // Reset arrays.
  info.clear();
  process.clear();
  event.clear();
  partonSystems.clear();
  beamA.clear();
  beamB.clear();
  beamPomA.clear();
  beamPomB.clear();
  beamGamA.clear();
  beamGamB.clear();
  beamVMDA.clear();
  beamVMDB.clear();

  // Pick current beam valence flavours (for pi0, K0S, K0L, Pomeron).
  beamA.newValenceContent();
  beamB.newValenceContent();
  if ( doDiffraction || doHardDiff) {
    beamPomA.newValenceContent();
    beamPomB.newValenceContent();
  }
  if (doVMDsideA) beamVMDA.newValenceContent();
  if (doVMDsideB) beamVMDB.newValenceContent();

  // Can only generate event if initialization worked.
  if (!isInit) {
    info.errorMsg("Abort from Pythia::next: "
      "not properly initialized so cannot generate events");
    return false;
  }

  // Pick beam momentum spread and beam vertex.
  if (doMomentumSpread || doVertexSpread) beamShapePtr->pick();

  // Recalculate kinematics when beam momentum spread.
  if (doMomentumSpread || doVarEcm) nextKinematics();

  // Outer loop over hard processes; only relevant for user-set vetoes.
  for ( ; ; ) {

    info.addCounter(10);
    bool hasVetoed = false;
    bool hasVetoedDiff = false;

    // Provide the hard process that starts it off. Only one try.
    info.clear();
    process.clear();
    partonSystems.clear();

    // Reset the event information. Necessary if the previous event was read
    // from LHEF, while the current event is not read from LHEF.
    info.setLHEF3EventInfo();

    if ( !processLevel.next( process) ) {
      if (doLHA && info.atEndOfFile()) info.errorMsg("Abort from "
        "Pythia::next: reached end of Les Houches Events File");
      else info.errorMsg("Abort from Pythia::next: "
        "processLevel failed; giving up");
      return false;
    }

    info.addCounter(11);

    // Update tried and selected events immediately after next event was
    // generated. Note: This does not accumulate cross section.
    processLevel.accumulate(false);

    // Possibility for a user veto of the process-level event.
    if (doVetoProcess) {
      hasVetoed = userHooksPtr->doVetoProcessLevel( process);
      if (hasVetoed) {
        if (abortIfVeto) return false;
        continue;
      }
    }

    // Possibility to perform matrix element merging for this event.
    if (doMerging) {
      int veto = mergingPtr->mergeProcess( process );
      // Apply possible merging scale cut.
      if (veto == -1) {
        hasVetoed = true;
        if (abortIfVeto) return false;
        continue;
      // Exit because of vanishing no-emission probability.
      } else if (veto == 0) {
        event = process;
        break;
      }

      // Redo resonance decays after the merging, in case the resonance
      // structure has been changed because of reclusterings.
      if (veto == 2 && doResDec) processLevel.nextDecays( process);
    }

    // Possibility to stop the generation at this stage.
    if (!doPartonLevel) {
      boostAndVertex( true, true);
      processLevel.accumulate();
      event.scale( process.scale() );
      event.scaleSecond( process.scaleSecond() );
      info.addCounter(4);
      if (doLHA && nPrevious < nShowLHA) lhaUpPtr->listEvent();
      if (nPrevious < nShowInfo) info.list();
      if (nPrevious < nShowProc) process.list(showSaV, showMaD);
      return true;
    }

    // Save spare copy of process record in case of problems.
    Event processSave = process;
    int sizeMPI       = info.sizeMPIarrays();
    info.addCounter(12);
    for (int i = 14; i < 19; ++i) info.setCounter(i);

    // Allow up to ten tries for parton- and hadron-level processing.
    bool physical   = true;
    for (int iTry = 0; iTry < NTRY; ++iTry) {

      info.addCounter(14);
      physical  = true;
      hasVetoed = false;

      // Restore original process record if problems.
      if (iTry > 0) process = processSave;
      if (iTry > 0) info.resizeMPIarrays( sizeMPI);

      // Reset event record and (extracted partons from) beam remnants.
      event.clear();
      beamA.clear();
      beamB.clear();
      beamPomA.clear();
      beamPomB.clear();
      beamGamA.clear();
      beamGamB.clear();
      beamVMDA.clear();
      beamVMDB.clear();
      partonSystems.clear();

      // Parton-level evolution: ISR, FSR, MPI.
      if ( !partonLevel.next( process, event) ) {

        // Abort event generation if parton level is set to abort.
        if (info.getAbortPartonLevel()) return false;

        // Skip to next hard process for failure owing to deliberate veto,
        // or alternatively retry for the same hard process.
        hasVetoed = partonLevel.hasVetoed();
        if (hasVetoed) {
          if (retryPartonLevel) {
            --iTry;
            continue;
          }
          if (abortIfVeto) return false;
          break;
        }

        // If hard diffractive event has been discarded retry partonLevel.
        hasVetoedDiff = partonLevel.hasVetoedDiff();
        if (hasVetoedDiff) {
          info.errorMsg("Warning in Pythia::next: "
            "discarding hard diffractive event from partonLevel; try again");
          break;
        }

        // Else make a new try for other failures.
        info.errorMsg("Error in Pythia::next: "
          "partonLevel failed; try again");
        physical = false;
        continue;
      }
      info.addCounter(15);

      // Possibility for a user veto of the parton-level event.
      if (doVetoPartons) {
        hasVetoed = userHooksPtr->doVetoPartonLevel( event);
        if (hasVetoed) {
          if (abortIfVeto) return false;
          break;
        }
      }

      // Boost to lab frame (before decays, for vertices).
      boostAndVertex( true, true);

      // Possibility to stop the generation at this stage.
      if (!doHadronLevel) {
        processLevel.accumulate();
        partonLevel.accumulate();
        event.scale( process.scale() );
        event.scaleSecond( process.scaleSecond() );
        // Optionally check final event for problems.
        if (checkEvent && !check()) {
          info.errorMsg("Abort from Pythia::next: "
            "check of event revealed problems");
          return false;
        }
        info.addCounter(4);
        if (doLHA && nPrevious < nShowLHA) lhaUpPtr->listEvent();
        if (nPrevious < nShowInfo) info.list();
        if (nPrevious < nShowProc) process.list(showSaV, showMaD);
        if (nPrevious < nShowEvt)  event.list(showSaV, showMaD);
        return true;
      }

      // Hadron-level: hadronization, decays.
      info.addCounter(16);
      if ( !hadronLevel.next( event) ) {
        info.errorMsg("Error in Pythia::next: "
          "hadronLevel failed; try again");
        physical = false;
        continue;
      }

      // If R-hadrons have been formed, then (optionally) let them decay.
      if (decayRHadrons && rHadrons.exist() && !doRHadronDecays()) {
        info.errorMsg("Error in Pythia::next: "
          "decayRHadrons failed; try again");
        physical = false;
        continue;
      }
      info.addCounter(17);

      // Optionally check final event for problems.
      if (checkEvent && !check()) {
        info.errorMsg("Error in Pythia::next: "
          "check of event revealed problems");
        physical = false;
        continue;
      }

      // Stop parton- and hadron-level looping if you got this far.
      info.addCounter(18);
      break;
    }

    // If event vetoed then to make a new try.
    if (hasVetoed || hasVetoedDiff)  {
      if (abortIfVeto) return false;
      continue;
    }

    // If event failed any other way (after ten tries) then give up.
    if (!physical) {
      info.errorMsg("Abort from Pythia::next: "
        "parton+hadronLevel failed; giving up");
      return false;
    }

    // Process- and parton-level statistics. Event scale.
    processLevel.accumulate();
    partonLevel.accumulate();
    event.scale( process.scale() );
    event.scaleSecond( process.scaleSecond() );

    // End of outer loop over hard processes. Done with normal option.
    info.addCounter(13);
    break;
  }

  // List events.
  if (doLHA && nPrevious < nShowLHA) lhaUpPtr->listEvent();
  if (nPrevious < nShowInfo) info.list();
  if (nPrevious < nShowProc) process.list(showSaV,showMaD);
  if (nPrevious < nShowEvt)  event.list(showSaV, showMaD);

  // Done.
  info.addCounter(4);
  return true;

}

//--------------------------------------------------------------------------

// Variant of the main event-generation routine, for variable CM energies.

bool Pythia::next(double eCMin) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Check that generation has been initialized for variable energy.
  if (!doVarEcm) {
    info.errorMsg("Abort from Pythia::next: "
      "generation not initialized for variable energies");
    return false;
  }

  // Check that the frameType matches the input provided.
  if (frameType != 1) {
    info.errorMsg("Abort from Pythia::next: "
      "input parameters do not match frame type");
    return false;
  }

  // Save input value.
  eCM = eCMin;

  // Call regular next method for event generation.
  return next();

}

//--------------------------------------------------------------------------

// Variant of the main event-generation routine, for variable beam energies.

bool Pythia::next(double eAin, double eBin) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Check that generation has been initialized for variable energy.
  if (!doVarEcm) {
    info.errorMsg("Abort from Pythia::next: "
      "generation not initialized for variable energies");
    return false;
  }

  // Check that the frameType matches the input provided.
  if (frameType != 2) {
    info.errorMsg("Abort from Pythia::next: "
      "input parameters do not match frame type");
    return false;
  }

  // Save input values.
  eA = eAin;
  eB = eBin;

  // Call regular next method for event generation.
  return next();

}

//--------------------------------------------------------------------------

// Variant of the main event-generation routine, for variable beam momenta.

bool Pythia::next(double pxAin, double pyAin, double pzAin,
                  double pxBin, double pyBin, double pzBin) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Check that generation has been initialized for variable energy.
  if (!doVarEcm) {
    info.errorMsg("Abort from Pythia::next: "
      "generation not initialized for variable energies");
    return false;
  }

  // Check that the frameType matches the input provided.
  if (frameType != 3) {
    info.errorMsg("Abort from Pythia::next: "
      "input parameters do not match frame type");
    return false;
  }

  // Save input value.
  pxA = pxAin;
  pyA = pyAin;
  pzA = pzAin;
  pxB = pxBin;
  pyB = pyBin;
  pzB = pzBin;

  // Call regular next method for event generation.
  return next();

}

//--------------------------------------------------------------------------

// Generate only the hadronization/decay stage, using internal machinery.
// The "event" instance should already contain a parton-level configuration.

bool Pythia::forceHadronLevel(bool findJunctions) {

  // Can only generate event if initialization worked.
  if (!isInit) {
    info.errorMsg("Abort from Pythia::forceHadronLevel: "
      "not properly initialized so cannot generate events");
    return false;
  }

  // Check whether any junctions in system. (Normally done in ProcessLevel.)
  // Avoid it if there are no final-state coloured partons.
  if (findJunctions) {
    event.clearJunctions();
    for (int i = 0; i < event.size(); ++i)
    if (event[i].isFinal()
    && (event[i].col() != 0 || event[i].acol() != 0)) {
      processLevel.findJunctions( event);
      break;
    }
  }

  // Allow for CR before the hadronization.
  if (forceHadronLevelCR) {

    // Setup parton system for SK-I and SK-II colour reconnection.
    // Require all final state particles to have the Ws as mothers.
    if (reconnectMode == 3 || reconnectMode == 4) {
      partonSystems.clear();
      partonSystems.addSys();
      partonSystems.addSys();
      partonSystems.setInRes(0, 3);
      partonSystems.setInRes(1, 4);
      for (int i = 5; i < event.size(); ++i) {
        if (event[i].mother1() - 3 < 0 || event[i].mother1() - 3 > 1) {
          info.errorMsg("Error in Pythia::forceHadronLevel: "
            " Event is not setup correctly for SK-I or SK-II CR");
          return false;
        }
        partonSystems.addOut(event[i].mother1() - 3,i);
      }
    }

    // save spare copy of event in case of failure.
    Event spareEvent = event;
    bool colCorrect = false;

    // Allow up to ten tries for CR.
    for (int iTry = 0; iTry < NTRY; ++ iTry) {
      colourReconnection.next(event, 0);
      if (junctionSplitting.checkColours(event)) {
        colCorrect = true;
        break;
      }
      else event = spareEvent;
    }

    if (!colCorrect) {
      info.errorMsg("Error in Pythia::forceHadronLevel: "
        "Colour reconnection failed.");
      return false;
    }
  }

  // Save spare copy of event in case of failure.
  Event spareEvent = event;

  // Allow up to ten tries for hadron-level processing.
  bool physical = true;
  for (int iTry = 0; iTry < NTRY; ++ iTry) {
    physical = true;

    // Check whether any resonances need to be handled at process level.
    if (doResDec) {
      process = event;
      processLevel.nextDecays( process);

      // Allow for showers if decays happened at process level.
      if (process.size() > event.size()) {
        if (doFSRinRes) {
          partonLevel.setupShowerSys( process, event);
          partonLevel.resonanceShowers( process, event, false);
        } else event = process;
      }
    }

    // Hadron-level: hadronization, decays.
    if (hadronLevel.next( event)) break;

    // If failure then warn, restore original configuration and try again.
    info.errorMsg("Error in Pythia::forceHadronLevel: "
      "hadronLevel failed; try again");
    physical = false;
    event    = spareEvent;
  }

  // Done for simpler option.
  if (!physical)  {
    info.errorMsg("Abort from Pythia::forceHadronLevel: "
      "hadronLevel failed; giving up");
    return false;
  }

  // Optionally check final event for problems.
  if (checkEvent && !check()) {
    info.errorMsg("Abort from Pythia::forceHadronLevel: "
      "check of event revealed problems");
    return false;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Recalculate kinematics for each event when beam momentum has a spread.

void Pythia::nextKinematics() {

  // Momentum spread: read out momentum shift to give current beam momenta.
  if (doMomentumSpread) {
    pAnow = pAinit + beamShapePtr->deltaPA();
    pAnow.e( sqrt(pAnow.pAbs2() + mA*mA) );
    pBnow = pBinit + beamShapePtr->deltaPB();
    pBnow.e( sqrt(pBnow.pAbs2() + mB*mB) );
    eCM   = (pAnow + pBnow).mCalc();

  // For variable energy in rest frame only need new eCM value, already set.
  } else if (frameType == 1) {

  // Variable energy but collinear beams: give current beam momenta.
  } else if (frameType == 2) {
    pAnow = Vec4( 0., 0.,  sqrtpos( eA*eA - mA*mA), eA);
    pBnow = Vec4( 0., 0., -sqrtpos( eB*eB - mB*mB), eB);
    eCM   = (pAnow + pBnow).mCalc();

  // Variable three-momenta stored and energy calculated.
  } else if (frameType == 3) {
    pAnow = Vec4( pxA, pyA, pzA, sqrt(pxA*pxA + pyA*pyA + pzA*pzA + mA*mA) );
    pBnow = Vec4( pxB, pyB, pzB, sqrt(pxB*pxB + pyB*pyB + pzB*pzB + mB*mB) );
    eCM   = (pAnow + pBnow).mCalc();

  // Other possibilites not supported.
  } else {
    info.errorMsg("Error from Pythia::nextKinematics: unsupported frameType");
    return;
  }

  // Construct CM frame kinematics.
  pzAcm = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
        * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
  pzBcm = -pzAcm;
  eA    = sqrt(mA*mA + pzAcm*pzAcm);
  eB    = sqrt(mB*mB + pzBcm*pzBcm);

  // Set relevant info for other classes to use.
  info.setBeamA( idA, pzAcm, eA, mA);
  info.setBeamB( idB, pzBcm, eB, mB);
  info.setECM( eCM);
  beamA.newPzE( pzAcm, eA);
  beamB.newPzE( pzBcm, eB);

  // Set boost/rotation matrices from/to CM frame.
  if (frameType != 1) {
    MfromCM.reset();
    MfromCM.fromCMframe( pAnow, pBnow);
    MtoCM = MfromCM;
    MtoCM.invert();
  }

}

//--------------------------------------------------------------------------

// Boost from CM frame to lab frame, or inverse. Set production vertex.

void Pythia::boostAndVertex( bool toLab, bool setVertex) {

  // Boost process from CM frame to lab frame.
  if (toLab) {
    if      (boostType == 2) process.bst(0., 0., betaZ, gammaZ);
    else if (boostType == 3) process.rotbst(MfromCM);

    // Boost nonempty event from CM frame to lab frame.
    if (event.size() > 0) {
      if      (boostType == 2) event.bst(0., 0., betaZ, gammaZ);
      else if (boostType == 3) event.rotbst(MfromCM);
    }

  // Boost process from lab frame to CM frame.
  } else {
    if      (boostType == 2) process.bst(0., 0., -betaZ, gammaZ);
    else if (boostType == 3) process.rotbst(MtoCM);

    // Boost nonempty event from lab frame to CM frame.
    if (event.size() > 0) {
      if      (boostType == 2) event.bst(0., 0., -betaZ, gammaZ);
      else if (boostType == 3) event.rotbst(MtoCM);
    }
  }

  // Set production vertex; assumes particles are in lab frame and at origin.
  if (setVertex && doVertexSpread) {
    Vec4 vertex = beamShapePtr->vertex();
    for (int i = 0; i < process.size(); ++i) process[i].vProdAdd( vertex);
    for (int i = 0; i < event.size(); ++i) event[i].vProdAdd( vertex);
  }

}

//--------------------------------------------------------------------------

// Perform R-hadron decays, either as part of normal evolution or forced.

bool Pythia::doRHadronDecays( ) {

  // Check if R-hadrons exist to be processed.
  if ( !rHadrons.exist() ) return true;

  // Do the R-hadron decay itself.
  if ( !rHadrons.decay( event) ) return false;

  // Perform showers in resonance decay chains.
  if ( !partonLevel.resonanceShowers( process, event, false) ) return false;

  // Subsequent hadronization and decays.
  if ( !hadronLevel.next( event) ) return false;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Print statistics on event generation.

void Pythia::stat() {

  if ( doHeavyIons ) {
    heavyIonsPtr->stat();
    return;
  }

  // Read out settings for what to include.
  bool showPrL = settings.flag("Stat:showProcessLevel");
  bool showPaL = settings.flag("Stat:showPartonLevel");
  bool showErr = settings.flag("Stat:showErrors");
  bool reset   = settings.flag("Stat:reset");

  // Statistics on cross section and number of events.
  if (doProcessLevel) {
    if (showPrL) processLevel.statistics(false);
    if (reset)   processLevel.resetStatistics();
  }

  // Statistics from other classes, currently multiparton interactions.
  if (showPaL) partonLevel.statistics(false);
  if (reset)   partonLevel.resetStatistics();

  // Merging statistics.
  if (doMerging) mergingPtr->statistics();

  // Summary of which and how many warnings/errors encountered.
  if (showErr) info.errorStatistics();
  if (reset)   info.errorReset();

}

//--------------------------------------------------------------------------

// Write the Pythia banner, with symbol and version information.

void Pythia::banner() {

  // Read in version number and last date of change.
  double versionNumber = settings.parm("Pythia:versionNumber");
  int versionDate = settings.mode("Pythia:versionDate");
  string month[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

  // Get date and time.
  time_t t = time(0);
  char dateNow[12];
  strftime(dateNow,12,"%d %b %Y",localtime(&t));
  char timeNow[9];
  strftime(timeNow,9,"%H:%M:%S",localtime(&t));

  cout << "\n"
       << " *-------------------------------------------"
       << "-----------------------------------------* \n"
       << " |                                           "
       << "                                         | \n"
       << " |  *----------------------------------------"
       << "--------------------------------------*  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   PPP   Y   Y  TTTTT  H   H  III    A  "
       << "    Welcome to the Lund Monte Carlo!  |  | \n"
       << " |  |   P  P   Y Y     T    H   H   I    A A "
       << "    This is PYTHIA version " << fixed << setprecision(3)
       << setw(5) << versionNumber << "      |  | \n"
       << " |  |   PPP     Y      T    HHHHH   I   AAAAA"
       << "    Last date of change: " << setw(2) << versionDate%100
       << " " << month[ (versionDate/100)%100 - 1 ]
       << " " << setw(4) << versionDate/10000 <<  "  |  | \n"
       << " |  |   P       Y      T    H   H   I   A   A"
       << "                                      |  | \n"
       << " |  |   P       Y      T    H   H  III  A   A"
       << "    Now is " << dateNow << " at " << timeNow << "    |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   Christian Bierlich;  Department of As"
       << "tronomy and Theoretical Physics,      |  | \n"
       << " |  |      Lund University, Solvegatan 14A, S"
       << "E-223 62 Lund, Sweden;                |  | \n"
       << " |  |      e-mail: christian.bierlich@thep.lu"
       << ".se                                   |  | \n"
       << " |  |   Nishita Desai;  Department of Theoret"
       << "ical Physics, Tata Institute,         |  | \n"
       << " |  |      Homi Bhabha Road, Mumbai 400005, I"
       << "ndia;                                 |  | \n"
       << " |  |      e-mail: desai@theory.tifr.res.in  "
       << "                                      |  | \n"
       << " |  |   Ilkka Helenius;  Department of Physic"
       << "s, University of Jyvaskyla,           |  | \n"
       << " |  |      P.O. Box 35, FI-40014 University o"
       << "f Jyvaskyla, Finland;                 |  | \n"
       << " |  |      e-mail: ilkka.m.helenius@jyu.fi   "
       << "                                      |  | \n"
       << " |  |   Philip Ilten;  School of Physics and "
       << "Astronomy,                            |  | \n"
       << " |  |      University of Birmingham, Birmingh"
       << "am, B152 2TT, UK;                     |  | \n"
       << " |  |      e-mail: philten@cern.ch           "
       << "                                      |  | \n"
       << " |  |   Leif Lonnblad;  Department of Astrono"
       << "my and Theoretical Physics,           |  | \n"
       << " |  |      Lund University, Solvegatan 14A, S"
       << "E-223 62 Lund, Sweden;                |  | \n"
       << " |  |      e-mail: leif.lonnblad@thep.lu.se  "
       << "                                      |  | \n"
       << " |  |   Stephen Mrenna;  Computing Division, "
       << "Simulations Group,                    |  | \n"
       << " |  |      Fermi National Accelerator Laborat"
       << "ory, MS 234, Batavia, IL 60510, USA;  |  | \n"
       << " |  |      e-mail: mrenna@fnal.gov           "
       << "                                      |  | \n"
       << " |  |   Stefan Prestel;  Department of Astron"
       << "omy and Theoretical Physics,          |  | \n"
       << " |  |      Lund University, Solvegatan 14A, S"
       << "E-223 62 Lund, Sweden;                |  | \n"
       << " |  |      e-mail: stefan.prestel@thep.lu.se "
       << "                                      |  | \n"
       << " |  |   Christine O. Rasmussen;  Department o"
       << "f Astronomy and Theoretical Physics,  |  | \n"
       << " |  |      Lund University, Solvegatan 14A, S"
       << "E-223 62 Lund, Sweden;                |  | \n"
       << " |  |      e-mail: christine.rasmussen@thep.l"
       << "u.se                                  |  | \n"
       << " |  |   Torbjorn Sjostrand;  Department of As"
       << "tronomy and Theoretical Physics,      |  | \n"
       << " |  |      Lund University, Solvegatan 14A, S"
       << "E-223 62 Lund, Sweden;                |  | \n"
       << " |  |      e-mail: torbjorn@thep.lu.se       "
       << "                                      |  | \n"
       << " |  |   Peter Skands;  School of Physics,    "
       << "                                      |  | \n"
       << " |  |      Monash University, PO Box 27, 3800"
       << " Melbourne, Australia;                |  | \n"
       << " |  |      e-mail: peter.skands@monash.edu   "
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   The main program reference is 'An Int"
       << "roduction to PYTHIA 8.2',             |  | \n"
       << " |  |   T. Sjostrand et al, Comput. Phys. Com"
       << "mun. 191 (2015) 159                   |  | \n"
       << " |  |   [arXiv:1410.3012 [hep-ph]]           "
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   The main physics reference is the 'PY"
       << "THIA 6.4 Physics and Manual',         |  | \n"
       << " |  |   T. Sjostrand, S. Mrenna and P. Skands"
       << ", JHEP05 (2006) 026 [hep-ph/0603175]  |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   An archive of program versions and do"
       << "cumentation is found on the web:      |  | \n"
       << " |  |   http://www.thep.lu.se/Pythia         "
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   This program is released under the GN"
       << "U General Public Licence version 2.   |  | \n"
       << " |  |   Please respect the MCnet Guidelines f"
       << "or Event Generator Authors and Users. |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   Disclaimer: this program comes withou"
       << "t any guarantees.                     |  | \n"
       << " |  |   Beware of errors and use common sense"
       << " when interpreting results.           |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   Copyright (C) 2019 Torbjorn Sjostrand"
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  *----------------------------------------"
       << "--------------------------------------*  | \n"
       << " |                                           "
       << "                                         | \n"
       << " *-------------------------------------------"
       << "-----------------------------------------* \n" << endl;

}

//--------------------------------------------------------------------------

// Check for lines in file that mark the beginning of new subrun.

int Pythia::readSubrun(string line, bool warn) {

  // If empty line then done.
  int subrunLine = SUBRUNDEFAULT;
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos)
    return subrunLine;

  // If first character is not a letter, then done.
  string lineNow = line;
  int firstChar = lineNow.find_first_not_of(" \n\t\v\b\r\f\a");
  if (!isalpha(lineNow[firstChar])) return subrunLine;

  // Replace an equal sign by a blank to make parsing simpler.
  while (lineNow.find("=") != string::npos) {
    int firstEqual = lineNow.find_first_of("=");
    lineNow.replace(firstEqual, 1, " ");
  }

  // Get first word of a line.
  istringstream splitLine(lineNow);
  string name;
  splitLine >> name;

  // Replace two colons by one (:: -> :) to allow for such mistakes.
  while (name.find("::") != string::npos) {
    int firstColonColon = name.find_first_of("::");
    name.replace(firstColonColon, 2, ":");
  }

  // Convert to lowercase. If no match then done.
  if (toLower(name) != "main:subrun") return subrunLine;

  // Else find new subrun number and return it.
  splitLine >> subrunLine;
  if (!splitLine) {
    if (warn) cout << "\n PYTHIA Warning: Main:subrun number not"
        << " recognized; skip:\n   " << line << endl;
    subrunLine = SUBRUNDEFAULT;
  }
  return subrunLine;

}

//--------------------------------------------------------------------------

// Check for lines in file that mark the beginning or end of commented section.
// Return +1 for beginning, -1 for end, 0 else.

int Pythia::readCommented(string line) {

  // If less than two nontrivial characters on line then done.
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return 0;
  int firstChar = line.find_first_not_of(" \n\t\v\b\r\f\a");
  if (int(line.size()) < firstChar + 2) return 0;

  // If first two nontrivial characters are /* or */ then done.
  if (line.substr(firstChar, 2) == "/*") return +1;
  if (line.substr(firstChar, 2) == "*/") return -1;

  // Else done.
  return 0;

}

//--------------------------------------------------------------------------

// Check that the final event makes sense: no unknown id codes;
// charge and energy-momentum conserved.

bool Pythia::check() {

  // Reset.
  bool physical     = true;
  bool listVertices = false;
  bool listHistory  = false;
  bool listSystems  = false;
  bool listBeams    = false;
  iErrId.resize(0);
  iErrCol.resize(0);
  iErrEpm.resize(0);
  iErrNan.resize(0);
  iErrNanVtx.resize(0);
  Vec4 pSum;
  double chargeSum  = 0.;

  // Incoming beams counted with negative momentum and charge.
  if (doProcessLevel) {
    // Incoming particles will be at position "1" and "2" in most cases.
    // However, need to be careful how to find incoming particles after
    // QED radiation in DIS-type collisions. Thus, first find both incoming
    // particles.
    int iA = 1;
    int iB = 2;
    if (!beamHasGamma) {
      if (beamA.isLepton() && beamB.isHadron())
        { iA = beamA[0].iPos(); iB = 2; }
      if (beamB.isLepton() && beamA.isHadron())
        { iB = beamB[0].iPos(); iA = 1; }
      int iPos = 0;
      while ( beamA.isHadron() && iPos < beamB.size()
           && beamA.id() == beamB[iPos].id() )
        { iA = beamA[iPos].iPos(); iPos++;}
      iPos = 0;
      while ( beamB.isHadron() && iPos < beamB.size()
           && beamB.id() == beamB[iPos].id() )
        { iB = beamB[iPos].iPos(); iPos++; }
    }
    // Count incoming momentum and charge.
    pSum      = - (event[iA].p() + event[iB].p());
    chargeSum = - (event[1].charge() + event[2].charge());

  // If no ProcessLevel then sum final state of process record.
  } else if (process.size() > 0) {
    pSum = - process[0].p();
    for (int i = 0; i < process.size(); ++i)
      if (process[i].isFinal()) chargeSum -= process[i].charge();

  // If process not filled, then use outgoing primary in event.
  } else {
    pSum = - event[0].p();
    for (int i = 1; i < event.size(); ++i)
      if (event[i].statusAbs() < 10 || event[i].statusAbs() == 23)
        chargeSum -= event[i].charge();
  }
  double eLab = abs(pSum.e());

  // Loop over particles in the event.
  for (int i = 0; i < event.size(); ++i) {

    // Look for any unrecognized particle codes.
    int id = event[i].id();
    if (id == 0 || !particleData.isParticle(id)) {
      ostringstream errCode;
      errCode << ", i = " << i << ", id = " << id;
      info.errorMsg("Error in Pythia::check: "
        "unknown particle code", errCode.str());
      physical = false;
      iErrId.push_back(i);

    // Check that colour assignments are the expected ones.
    } else {
      int colType = event[i].colType();
      int col     = event[i].col();
      int acol    = event[i].acol();
      if ( event[i].statusAbs() / 10 == 8 ) acol = col = 0;
      if ( (colType ==  0 && (col  > 0 || acol  > 0))
        || (colType ==  1 && (col <= 0 || acol  > 0))
        || (colType == -1 && (col  > 0 || acol <= 0))
        || (colType ==  2 && (col <= 0 || acol <= 0)) ) {
        ostringstream errCode;
        errCode << ", i = " << i << ", id = " << id << " cols = " << col
                << " " << acol;
        info.errorMsg("Error in Pythia::check: "
          "incorrect colours", errCode.str());
        physical = false;
        iErrCol.push_back(i);
      }
    }

    // Some intermediate shower partons excepted from (E, p, m) consistency.
    bool checkMass = event[i].statusAbs() != 49 && event[i].statusAbs() != 59;

    // Look for particles with mismatched or not-a-number energy/momentum/mass.
    if (abs(event[i].px()) >= 0. && abs(event[i].py()) >= 0.
      && abs(event[i].pz()) >= 0.  && abs(event[i].e()) >= 0.
      && abs(event[i].m()) >= 0.) {
      double errMass = abs(event[i].mCalc() - event[i].m())
        / max( 1.0, event[i].e());
      if (checkMass && errMass > mTolErr) {
        info.errorMsg("Error in Pythia::check: "
          "unmatched particle energy/momentum/mass");
        physical = false;
        iErrEpm.push_back(i);
      } else if (checkMass && errMass > mTolWarn) {
        info.errorMsg("Warning in Pythia::check: "
          "not quite matched particle energy/momentum/mass");
      }
    } else {
      info.errorMsg("Error in Pythia::check: "
        "not-a-number energy/momentum/mass");
      physical = false;
      iErrNan.push_back(i);
    }

    // Look for particles with not-a-number vertex/lifetime.
    if (abs(event[i].xProd()) >= 0. && abs(event[i].yProd()) >= 0.
      && abs(event[i].zProd()) >= 0.  && abs(event[i].tProd()) >= 0.
      && abs(event[i].tau()) >= 0.) ;
    else {
      info.errorMsg("Error in Pythia::check: "
        "not-a-number vertex/lifetime");
      physical     = false;
      listVertices = true;
      iErrNanVtx.push_back(i);
    }

    // Add final-state four-momentum and charge.
    if (event[i].isFinal()) {
      pSum      += event[i].p();
      chargeSum += event[i].charge();
    }

  // End of particle loop.
  }

  // Check energy-momentum/charge conservation.
  double epDev = abs(pSum.e()) + abs(pSum.px()) + abs(pSum.py())
    + abs(pSum.pz());
  if (epDev > epTolErr * eLab) {
    info.errorMsg("Error in Pythia::check: energy-momentum not conserved");
    physical = false;
  } else if (epDev > epTolWarn * eLab) {
    info.errorMsg("Warning in Pythia::check: "
      "energy-momentum not quite conserved");
  }
  if (abs(chargeSum) > 0.1) {
    info.errorMsg("Error in Pythia::check: charge not conserved");
    physical = false;
  }

  // Check that beams and event records agree on incoming partons.
  // Only meaningful for resolved beams.
  if (info.isResolved() && !info.hasUnresolvedBeams())
  for (int iSys = 0; iSys < beamA.sizeInit(); ++iSys) {
    int eventANw  = partonSystems.getInA(iSys);
    int eventBNw  = partonSystems.getInB(iSys);
    // For photon beams from leptons make sure to use correct beams.
    int beamANw   = ( beamA.getGammaMode() == 0 || !beamHasGamma
                 || (beamA.getGammaMode() == 2 && beamB.getGammaMode() == 2)) ?
                 beamA[iSys].iPos() : beamGamA[iSys].iPos();
    int beamBNw   = ( beamB.getGammaMode() == 0 || !beamHasGamma
                 || (beamB.getGammaMode() == 2 && beamA.getGammaMode() == 2)) ?
                 beamB[iSys].iPos() : beamGamB[iSys].iPos();
    if (eventANw != beamANw || eventBNw != beamBNw) {
      info.errorMsg("Error in Pythia::check: "
        "event and beams records disagree");
      physical    = false;
      listSystems = true;
      listBeams   = true;
    }
  }

  // Check that mother and daughter information match for each particle.
  vector<int> noMot;
  vector<int> noDau;
  vector< pair<int,int> > noMotDau;
  if (checkHistory) {

    // Loop through the event and check that there are beam particles.
    bool hasBeams = false;
    for (int i = 0; i < event.size(); ++i) {
      int status = event[i].status();
      if (abs(status) == 12) hasBeams = true;

      // Check that mother and daughter lists not empty where not expected to.
      vector<int> mList = event[i].motherList();
      vector<int> dList = event[i].daughterList();
      if (mList.size() == 0 && abs(status) != 11 && abs(status) != 12)
        noMot.push_back(i);
      if (dList.size() == 0 && status < 0 && status != -11)
        noDau.push_back(i);

      // Check that the particle appears in the daughters list of each mother.
      for (int j = 0; j < int(mList.size()); ++j) {
        if ( event[mList[j]].daughter1() <= i
          && event[mList[j]].daughter2() >= i ) continue;
        vector<int> dmList = event[mList[j]].daughterList();
        bool foundMatch = false;
        for (int k = 0; k < int(dmList.size()); ++k)
        if (dmList[k] == i) {
          foundMatch = true;
          break;
        }
        if (!hasBeams && mList.size() == 1 && mList[0] == 0) foundMatch = true;
        if (!foundMatch) {
          bool oldPair = false;
          for (int k = 0; k < int(noMotDau.size()); ++k)
          if (noMotDau[k].first == mList[j] && noMotDau[k].second == i) {
            oldPair = true;
            break;
          }
          if (!oldPair) noMotDau.push_back( make_pair( mList[j], i) );
        }
      }

      // Check that the particle appears in the mothers list of each daughter.
      for (int j = 0; j < int(dList.size()); ++j) {
        if ( event[dList[j]].statusAbs() > 80
          && event[dList[j]].statusAbs() < 90
          && event[dList[j]].mother1() <= i
          && event[dList[j]].mother2() >= i) continue;
        vector<int> mdList = event[dList[j]].motherList();
        bool foundMatch = false;
        for (int k = 0; k < int(mdList.size()); ++k)
        if (mdList[k] == i) {
          foundMatch = true;
          break;
        }
        if (!foundMatch) {
          bool oldPair = false;
          for (int k = 0; k < int(noMotDau.size()); ++k)
          if (noMotDau[k].first == i && noMotDau[k].second == dList[j]) {
            oldPair = true;
            break;
          }
          if (!oldPair) noMotDau.push_back( make_pair( i, dList[j]) );
        }
      }
    }

    // Warn if any errors were found.
    if (noMot.size() > 0 || noDau.size() > 0 || noMotDau.size() > 0) {
      info.errorMsg("Error in Pythia::check: "
        "mismatch in daughter and mother lists");
      physical    = false;
      listHistory = true;
    }
  }

  // Done for sensible events.
  if (physical) return true;

  // Print (the first few) flawed events: local info.
  if (nErrEvent < nErrList) {
    cout << "\n PYTHIA erroneous event info: \n";
    if (iErrId.size() > 0) {
      cout << " unknown particle codes in lines ";
      for (int i = 0; i < int(iErrId.size()); ++i)
        cout << iErrId[i] << " ";
      cout << "\n";
    }
    if (iErrCol.size() > 0) {
      cout << " incorrect colour assignments in lines ";
      for (int i = 0; i < int(iErrCol.size()); ++i)
        cout << iErrCol[i] << " ";
      cout << "\n";
    }
    if (iErrEpm.size() > 0) {
      cout << " mismatch between energy/momentum/mass in lines ";
      for (int i = 0; i < int(iErrEpm.size()); ++i)
        cout << iErrEpm[i] << " ";
      cout << "\n";
    }
    if (iErrNan.size() > 0) {
      cout << " not-a-number energy/momentum/mass in lines ";
      for (int i = 0; i < int(iErrNan.size()); ++i)
        cout << iErrNan[i] << " ";
      cout << "\n";
    }
    if (iErrNanVtx.size() > 0) {
      cout << " not-a-number vertex/lifetime in lines ";
      for (int i = 0; i < int(iErrNanVtx.size()); ++i)
        cout << iErrNanVtx[i] << " ";
      cout << "\n";
    }
    if (epDev > epTolErr * eLab) cout << scientific << setprecision(3)
      << " total energy-momentum non-conservation = " << epDev << "\n";
    if (abs(chargeSum) > 0.1) cout << fixed << setprecision(2)
      << " total charge non-conservation = " << chargeSum << "\n";
    if (noMot.size() > 0) {
      cout << " missing mothers for particles ";
      for (int i = 0; i < int(noMot.size()); ++i) cout << noMot[i] << " ";
      cout << "\n";
    }
    if (noDau.size() > 0) {
      cout << " missing daughters for particles ";
      for (int i = 0; i < int(noDau.size()); ++i) cout << noDau[i] << " ";
      cout << "\n";
    }
    if (noMotDau.size() > 0) {
      cout << " inconsistent history for (mother,daughter) pairs ";
      for (int i = 0; i < int(noMotDau.size()); ++i)
        cout << "(" << noMotDau[i].first << "," << noMotDau[i].second << ") ";
      cout << "\n";
    }

    // Print (the first few) flawed events: standard listings.
    info.list();
    event.list(listVertices, listHistory);
    if (listSystems) partonSystems.list();
    if (listBeams) beamA.list();
    if (listBeams) beamB.list();
  }

  // Update error counter. Done also for flawed event.
  ++nErrEvent;
  return false;

}

//--------------------------------------------------------------------------

// Routine to set up a PDF pointer.

PDF* Pythia::getPDFPtr(int idIn, int sequence, string beam, bool resolved) {

  // Temporary pointer to be returned.
  PDF* tempPDFPtr = 0;

  // One option is to treat a Pomeron like a pi0.
  if (idIn == 990 && settings.word("PDF:PomSet") == "2") idIn = 111;

  // Proton beam, normal or hard choice. Also used for neutron.
  if (abs(idIn) == 2212 || abs(idIn) == 2112) {
    string pWord = settings.word("PDF:p"
      + string(sequence == 1 ? "" : "Hard") + "Set" + beam);
    if (pWord == "void" && sequence != 1 && beam == "B")
      pWord = settings.word("PDF:pHardSet");
    if (pWord == "void") pWord = settings.word("PDF:pSet");
    istringstream pStream(pWord);
    int pSet = 0;
    pStream >> pSet;

    // Use internal LHAgrid1 implementation for LHAPDF6 files.
    if (pSet == 0 && pWord.length() > 9
      && toLower(pWord).substr(0,9) == "lhagrid1:")
      tempPDFPtr = new LHAGrid1(idIn, pWord, xmlPath, &info);

    // Use sets from LHAPDF.
    else if (pSet == 0) tempPDFPtr = new LHAPDF(idIn, pWord, &info);

    // Use internal sets.
    else if (pSet == 1) tempPDFPtr = new GRV94L(idIn);
    else if (pSet == 2) tempPDFPtr = new CTEQ5L(idIn);
    else if (pSet <= 6)
      tempPDFPtr = new MSTWpdf(idIn, pSet - 2, xmlPath, &info);
    else if (pSet <= 12)
      tempPDFPtr = new CTEQ6pdf(idIn, pSet - 6, 1., xmlPath, &info);
    else if (pSet <= 16)
      tempPDFPtr = new NNPDF(idIn, pSet - 12, xmlPath, &info);
    else if (pSet <= 22)
      tempPDFPtr = new LHAGrid1(idIn, pWord, xmlPath, &info);
    else tempPDFPtr = 0;
  }

  // Pion beam (or, in one option, Pomeron beam).
  else if (abs(idIn) == 211 || idIn == 111) {
    string piWord = settings.word("PDF:piSet" + beam);
    if (piWord == "void" && beam == "B") piWord = settings.word("PDF:piSet");
    istringstream piStream(piWord);
    int piSet = 0;
    piStream >> piSet;

    // If VMD process then scale PDF accordingly:
    // f_a^VMD = alphaEM * (1/f_rho^2 + 1/f_omega^2 + 1/f_phi^2 + 1/f_J/psi)
    //         * f_a^pi0.
    // COR: New value here includes J/psi
    double rescale = (doVMDsideA || doVMDsideB) ? 0.0046549 : 1.;

    // Use internal LHAgrid1 implementation for LHAPDF6 files.
    if (piSet == 0 && piWord.length() > 9
      && toLower(piWord).substr(0,9) == "lhagrid1:")
      tempPDFPtr = new LHAGrid1(idIn, piWord, xmlPath, &info);

    // Use sets from LHAPDF.
    else if (piSet == 0) tempPDFPtr = new LHAPDF(idIn, piWord, &info);

    // Use internal set.
    else if (piSet == 1) tempPDFPtr = new GRVpiL(idIn, rescale);
    else tempPDFPtr = 0;
  }

  // Pomeron beam, if not treated like a pi0 beam.
  else if (idIn == 990) {
    string pomWord = settings.word("PDF:PomSet");
    double rescale = settings.parm("PDF:PomRescale");
    istringstream pomStream(pomWord);
    int pomSet = 0;
    pomStream >> pomSet;

    // Use internal LHAgrid1 implementation for LHAPDF6 files.
    if (pomSet == 0 && pomWord.length() > 9
      && toLower(pomWord).substr(0,9) == "lhagrid1:")
      tempPDFPtr = new LHAGrid1(idIn, pomWord, xmlPath, &info);

    // Use sets from LHAPDF.
    else if (pomSet == 0) tempPDFPtr = new LHAPDF(idIn, pomWord, &info);

    // A generic Q2-independent parametrization.
    else if (pomSet == 1) {
      double gluonA      = settings.parm("PDF:PomGluonA");
      double gluonB      = settings.parm("PDF:PomGluonB");
      double quarkA      = settings.parm("PDF:PomQuarkA");
      double quarkB      = settings.parm("PDF:PomQuarkB");
      double quarkFrac   = settings.parm("PDF:PomQuarkFrac");
      double strangeSupp = settings.parm("PDF:PomStrangeSupp");
      tempPDFPtr = new PomFix( 990, gluonA, gluonB, quarkA, quarkB,
        quarkFrac, strangeSupp);
    }

    // The H1 Q2-dependent parametrizations. Initialization requires files.
    else if (pomSet == 3 || pomSet == 4)
      tempPDFPtr = new PomH1FitAB( 990, pomSet - 2, rescale, xmlPath, &info);
    else if (pomSet == 5)
      tempPDFPtr = new PomH1Jets( 990, 1, rescale, xmlPath, &info);
    else if (pomSet == 6)
      tempPDFPtr = new PomH1FitAB( 990, 3, rescale, xmlPath, &info);
    // The parametrizations of Alvero, Collins, Terron and Whitmore.
    else if (pomSet > 6 && pomSet < 11)  {
      tempPDFPtr = new CTEQ6pdf( 990, pomSet + 4, rescale, xmlPath, &info);
      info.errorMsg("Warning: Pomeron flux parameters forced for ACTW PDFs");
      settings.mode("SigmaDiffractive:PomFlux", 4);
      double pomFluxEps = (pomSet == 10) ? 0.19 : 0.14;
      settings.parm("SigmaDiffractive:PomFluxEpsilon", pomFluxEps);
      settings.parm("SigmaDiffractive:PomFluxAlphaPrime", 0.25);
    }
    else if (pomSet == 11 )
      tempPDFPtr = new PomHISASD(990, getPDFPtr(2212), settings, &info);
    else if (pomSet >= 12 && pomSet <= 15)
      tempPDFPtr = new LHAGrid1(idIn, "1" + pomWord, xmlPath, &info);
    else tempPDFPtr = 0;
  }

  // Set up nuclear PDFs.
  else if (idIn > 100000000) {

    // Which nPDF set to use.
    int nPDFSet = (beam == "B") ? settings.mode("PDF:nPDFSetB")
                                : settings.mode("PDF:nPDFSetA");

    // Temporary pointer for storing proton PDF pointer.
    PDF* tempProtonPDFPtr = (beam == "B") ? pdfHardBPtr : pdfHardAPtr;
    if (nPDFSet == 0) tempPDFPtr = new Isospin(idIn, tempProtonPDFPtr);
    else if (nPDFSet == 1 || nPDFSet == 2)  tempPDFPtr = new EPS09(idIn,
      nPDFSet, 1, xmlPath, tempProtonPDFPtr, &info);
    else if (nPDFSet == 3) tempPDFPtr = new EPPS16(idIn, 1, xmlPath,
      tempProtonPDFPtr, &info);
    else tempPDFPtr = 0;
  }

  // Photon beam, either point-like (unresolved) or resolved.
  else if (abs(idIn) == 22) {

    // For unresolved beam use the point-like PDF.
    if (!resolved) {
      tempPDFPtr = new GammaPoint(idIn);
    } else {
      int gammaSet = settings.mode("PDF:GammaSet");

      // Point-like beam if unresolved photons.
      bool beamAisPoint = ( !beamAisResGamma && !beamAhasResGamma );
      bool beamBisPoint = ( !beamBisResGamma && !beamBhasResGamma );
      bool beamIsPoint = ( beamAisPoint && !(beam == "B") )
                      || ( beamBisPoint && (beam == "B") );

      // Use different PDFs for hard process.
      if ( sequence == 2) {

        // Find the name or number of the hard PDF set.
        string gmWord = settings.word("PDF:GammaHardSet");
        int gmSet     = 0;
        if (gmWord == "void") gmSet = settings.mode("PDF:GammaSet");
        else {
          istringstream gmStream(gmWord);
          gmStream >> gmSet;
        }

        // Use sets from LHAPDF. Only available for hard processes.
        if (gmSet == 0 && !beamIsPoint) {
          tempPDFPtr = new LHAPDF(idIn, gmWord, &info);
          return tempPDFPtr;
        }

        // Or set up an internal set (though currently only one).
        gammaSet = gmSet;
      }

      // Set up the PDF.
      if      (beamIsPoint)   tempPDFPtr = new GammaPoint(idIn);
      else if (gammaSet == 1) tempPDFPtr = new CJKL(idIn, &rndm);
      else                    tempPDFPtr = 0;
    }
  }

  // Lepton beam: neutrino, resolved charged lepton or unresolved ditto.
  // Also photon inside lepton PDFs.
  else if (abs(idIn) > 10 && abs(idIn) < 17) {
    if (abs(idIn)%2 == 0) tempPDFPtr = new NeutrinoPoint(idIn);

    // Set up resolved photon inside lepton for beam A.
    if ( beamAhasResGamma && !(beam == "B") && resolved ) {

      // Find the pre-set photon PDF, hard or normal.
      PDF* tempGammaPDFPtr = 0;
      if ( sequence == 2) tempGammaPDFPtr = pdfHardGamAPtr;
      else                tempGammaPDFPtr = pdfGamAPtr;

      // Get the mass of lepton and maximum virtuality of the photon.
      double m2beam     = pow2(particleData.m0(idIn));
      double Q2maxGamma = settings.parm("Photon:Q2max");

      // Initialize the gamma-inside-lepton PDFs with internal photon flux.
      if (settings.mode("PDF:lepton2gammaSet") == 1) {
        tempPDFPtr = new Lepton2gamma(idIn, m2beam, Q2maxGamma,
          tempGammaPDFPtr, &info, &rndm);

      // Initialize the gamma-inside-lepton PDFs with external photon flux.
      // Requires that the pointer to the flux set.
      } else if ( settings.mode("PDF:lepton2gammaSet") == 2 ) {
        PDF* tempGammaFluxPtr = pdfGamFluxAPtr;
        if ( tempGammaFluxPtr != 0) tempPDFPtr = new EPAexternal(idIn, m2beam,
          tempGammaFluxPtr, tempGammaPDFPtr, &settings, &info, &rndm );
        else {
          tempPDFPtr = 0;
          info.errorMsg("Error in Pythia::getPDFPtr: "
            "No external photon flux provided with PDF:lepton2gammaSet == 2");
        }
      } else tempPDFPtr = 0;

    // Set up resolved photon inside lepton for beam B.
    } else if ( beamBhasResGamma && (beam == "B") && resolved ) {

      // Find the pre-set photon PDF, hard or normal.
      PDF* tempGammaPDFPtr = 0;
      if ( sequence == 2) tempGammaPDFPtr = pdfHardGamBPtr;
      else                tempGammaPDFPtr = pdfGamBPtr;

      // Get the mass of lepton and maximum virtuality of the photon.
      double m2beam     = pow2(particleData.m0(idIn));
      double Q2maxGamma = settings.parm("Photon:Q2max");

      // Initialize the gamma-inside-lepton PDFs with internal photon flux.
      if (settings.mode("PDF:lepton2gammaSet") == 1) {
        tempPDFPtr = new Lepton2gamma(idIn, m2beam, Q2maxGamma,
          tempGammaPDFPtr, &info, &rndm);

      // Initialize the gamma-inside-lepton PDFs with external photon flux.
      } else if ( settings.mode("PDF:lepton2gammaSet") == 2 ) {
        PDF* tempGammaFluxPtr = pdfGamFluxBPtr;
        if ( tempGammaFluxPtr != 0) tempPDFPtr = new EPAexternal(idIn, m2beam,
          tempGammaFluxPtr, tempGammaPDFPtr, &settings, &info, &rndm );
        else {
          tempPDFPtr = 0;
          info.errorMsg("Error in Pythia::getPDFPtr: "
            "No external photon flux provided with PDF:lepton2gammaSet == 2");
        }
      } else tempPDFPtr = 0;

    // Usual lepton PDFs.
    } else if (settings.flag("PDF:lepton")) {
      double m2beam = pow2(particleData.m0(idIn));
      double Q2maxGamma = settings.parm("Photon:Q2max");
      if (settings.mode("PDF:lepton2gammaSet") == 1 ) {
        tempPDFPtr = new Lepton(idIn, Q2maxGamma, &info, &rndm);

      // External photon flux for direct-photon processes.
      } else if (settings.mode("PDF:lepton2gammaSet") == 2 ) {
        PDF* tempGammaPDFPtr  = 0;
        PDF* tempGammaFluxPtr = (beam == "B") ?
          pdfGamFluxBPtr : pdfGamFluxAPtr;
        if ( tempGammaFluxPtr != 0) tempPDFPtr = new EPAexternal(idIn, m2beam,
          tempGammaFluxPtr, tempGammaPDFPtr, &settings, &info, &rndm );
        else {
          tempPDFPtr = 0;
          info.errorMsg("Error in Pythia::getPDFPtr: "
            "No external photon flux provided with PDF:lepton2gammaSet == 2");
        }
      } else tempPDFPtr = 0;
    }
    else tempPDFPtr = new LeptonPoint(idIn);

  // Dark matter beam set up as pointlike lepton.
  } else if (abs(idIn) > 50 && abs(idIn) < 60) {
    tempPDFPtr = new LeptonPoint(idIn);
  }

  // Optionally allow extrapolation beyond x and Q2 limits.
  if (tempPDFPtr)
    tempPDFPtr->setExtrapolate( settings.flag("PDF:extrapolate") );

  // Done.
  return tempPDFPtr;
}

//==========================================================================

} // end namespace Pythia8
