// SLHAinterface.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// Main authors of this file: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "Pythia8/SLHAinterface.h"

namespace Pythia8 {

//==========================================================================

// The SLHAinterface class.

//--------------------------------------------------------------------------

// Initialize and switch to SUSY couplings if reading SLHA spectrum.

void SLHAinterface::init( Settings& settings, Rndm* rndmPtr,
  Couplings* couplingsPtrIn, ParticleData* particleDataPtr,
  bool& useSLHAcouplings) {

  // Initialize SLHA couplingsPtr to PYTHIA one by default
  couplingsPtr     = couplingsPtrIn;
  useSLHAcouplings = false;

  // Check if SUSY couplings need to be read in
  if( !initSLHA(settings, particleDataPtr))
    infoPtr->errorMsg("Error in SLHAinterface::init: "
      "Could not read SLHA file");

  // SLHA sets isSUSY flag to tell us if there was an SLHA SUSY spectrum
  if (couplingsPtr->isSUSY) {
    // Initialize the derived SUSY couplings class (SM first, then SUSY)
    coupSUSY.init( settings, rndmPtr);
    coupSUSY.initSUSY(&slha, infoPtr, particleDataPtr, &settings);
    // Switch couplingsPtr to point to the derived class
    // and tell PYTHIA to use it
    couplingsPtr = (Couplings *) &coupSUSY;
    useSLHAcouplings = true;
  }

}

//--------------------------------------------------------------------------

// Initialize SUSY Les Houches Accord data.

bool SLHAinterface::initSLHA(Settings& settings,
  ParticleData* particleDataPtr) {

  // Error and warning prefixes for this method
  string errPref  = "Error in SLHAinterface::initSLHA: ";
  string warnPref = "Warning in SLHAinterface::initSLHA: ";
  string infoPref = "Info from SLHAinterface::initSLHA: ";

  // Initial and settings values.
  int    ifailLHE    = 1;
  int    ifailSpc    = 1;
  int    readFrom    = settings.mode("SLHA:readFrom");
  string lhefFile    = settings.word("Beams:LHEF");
  string lhefHeader  = settings.word("Beams:LHEFheader");
  string slhaFile    = settings.word("SLHA:file");
  int    verboseSLHA = settings.mode("SLHA:verbose");
  bool   slhaUseDec  = settings.flag("SLHA:useDecayTable");
  bool   noSLHAFile  = ( slhaFile == "none" || slhaFile == "void"
                      || slhaFile == ""     || slhaFile == " " );

  // Set internal data members
  meMode      = settings.mode("SLHA:meMode");

  // No SUSY by default
  couplingsPtr->isSUSY = false;

  // Option with no SLHA read-in at all.
  if (readFrom == 0) return true;

  // First check LHEF header (if reading from LHEF)
  if (readFrom == 1) {
    // Check if there is a <slha> tag in the LHEF header
    // Note: if the <slha> tag is NOT inside the <header>, it will be ignored.
    string slhaInHeader( infoPtr->header("slha") );
    if (slhaInHeader == "" && noSLHAFile) return true;
    // If there is an <slha> tag, read file.
    if (lhefHeader != "void")
      ifailLHE = slha.readFile(lhefHeader, verboseSLHA, slhaUseDec );
    else if (lhefFile != "void")
      ifailLHE = slha.readFile(lhefFile, verboseSLHA, slhaUseDec );
    else if (noSLHAFile) {
      istringstream slhaInHeaderStream(slhaInHeader);
      ifailLHE = slha.readFile(slhaInHeaderStream, verboseSLHA, slhaUseDec );
    }
  }

  // If LHEF read successful, everything needed should already be ready
  if (ifailLHE == 0) {
    ifailSpc = 0;
    couplingsPtr->isSUSY = true;
    // If no LHEF file or no SLHA info in header, read from SLHA:file
  } else {
    lhefFile = "void";
    if ( noSLHAFile ) return true;
    ifailSpc = slha.readFile(slhaFile,verboseSLHA, slhaUseDec);
  }

  // In case of problems, print error and fail init.
  if (ifailSpc != 0) {
    infoPtr->errorMsg(errPref + "problem reading SLHA file", slhaFile);
    return false;
  } else {
    couplingsPtr->isSUSY = true;
  }

  // Check spectrum for consistency. Switch off SUSY if necessary.
  ifailSpc = slha.checkSpectrum();

  // ifail >= 1 : no MODSEL found -> don't switch on SUSY
  if (ifailSpc == 1) {
    // no SUSY, but MASS ok
    couplingsPtr->isSUSY = false;
    infoPtr->errorMsg(infoPref +
      "No MODSEL found, keeping internal SUSY switched off");
  } else if (ifailSpc >= 2) {
    // no SUSY, but problems
    infoPtr->errorMsg(warnPref + "Problem with SLHA MASS or QNUMBERS.");
    couplingsPtr->isSUSY = false;
  }
  // ifail = 0 : MODSEL found, spectrum OK
  else if (ifailSpc == 0) {
    // Print spectrum. Done.
    slha.printSpectrum(0);
  }
  else if (ifailSpc < 0) {
    infoPtr->errorMsg(warnPref + "Problem with SLHA spectrum.",
      "\n Only using masses and switching off SUSY.");
    settings.flag("SUSY:all", false);
    couplingsPtr->isSUSY = false;
    slha.printSpectrum(ifailSpc);
  }

  // SLHA1 : SLHA2 compatibility
  // Check whether scalar particle masses are ordered
  bool isOrderedQ = true;
  bool isOrderedL = true;
  int idSdown[6]={1000001,1000003,1000005,2000001,2000003,2000005};
  int idSup[6]={1000002,1000004,1000006,2000002,2000004,2000006};
  int idSlep[6]={1000011,1000013,1000015,2000011,2000013,2000015};
  for (int j=0;j<=4;j++) {
    if (slha.mass(idSlep[j+1]) < slha.mass(idSlep[j]))
      isOrderedL  = false;
    if (slha.mass(idSup[j+1]) < slha.mass(idSup[j]))
      isOrderedQ  = false;
    if (slha.mass(idSdown[j+1]) < slha.mass(idSdown[j]))
      isOrderedQ  = false;
  }

  // If ordered, change sparticle labels to mass-ordered enumeration
  for (int i=1;i<=6;i++) {
    ostringstream indx;
    indx << i;
    string uName = "~u_"+indx.str();
    string dName = "~d_"+indx.str();
    string lName = "~e_"+indx.str();
    if (isOrderedQ) {
      particleDataPtr->names(idSup[i],uName,uName+"bar");
      particleDataPtr->names(idSdown[i],dName,dName+"bar");
    }
    if (isOrderedL) particleDataPtr->names(idSlep[i],lName+"-",lName+"+");
  }

  // NMSSM spectrum (modify existing Higgs names and add particles)
  if ( (ifailSpc == 1 || ifailSpc == 0) &&  slha.modsel(3) >= 1 ) {
    // Modify Higgs names
    particleDataPtr->name(25,"H_1");
    particleDataPtr->name(35,"H_2");
    particleDataPtr->name(45,"H_3");
    particleDataPtr->name(36,"A_1");
    particleDataPtr->name(46,"A_2");
    particleDataPtr->name(1000045,"~chi_50");
  }

  // SLHA2 spectrum with flavour mixing (modify squark and/or slepton names)
  if ( (ifailSpc == 1 || ifailSpc == 0) &&  slha.modsel(6) >= 1 ) {
    // Squark flavour violation
    if ( (slha.modsel(6) == 1 || slha.modsel(6) >= 3)
         && slha.usqmix.exists() && slha.dsqmix.exists() ) {
      // Modify squark names
      particleDataPtr->names(1000001,"~d_1","~d_1bar");
      particleDataPtr->names(1000002,"~u_1","~u_1bar");
      particleDataPtr->names(1000003,"~d_2","~d_2bar");
      particleDataPtr->names(1000004,"~u_2","~u_2bar");
      particleDataPtr->names(1000005,"~d_3","~d_3bar");
      particleDataPtr->names(1000006,"~u_3","~u_3bar");
      particleDataPtr->names(2000001,"~d_4","~d_4bar");
      particleDataPtr->names(2000002,"~u_4","~u_4bar");
      particleDataPtr->names(2000003,"~d_5","~d_5bar");
      particleDataPtr->names(2000004,"~u_5","~u_5bar");
      particleDataPtr->names(2000005,"~d_6","~d_6bar");
      particleDataPtr->names(2000006,"~u_6","~u_6bar");
    }
    // Slepton flavour violation
    if ( (slha.modsel(6) == 2 || slha.modsel(6) >= 3)
         && slha.selmix.exists()) {
      // Modify slepton names
      particleDataPtr->names(1000011,"~e_1-","~e_1+");
      particleDataPtr->names(1000013,"~e_2-","~e_2+");
      particleDataPtr->names(1000015,"~e_3-","~e_3+");
      particleDataPtr->names(2000011,"~e_4-","~e_4+");
      particleDataPtr->names(2000013,"~e_5-","~e_5+");
      particleDataPtr->names(2000015,"~e_6-","~e_6+");
    }
    // Neutrino flavour violation
    if ( (slha.modsel(6) == 2 || slha.modsel(6) >= 3)
         && slha.upmns.exists()) {
      // Modify neutrino names (note that SM processes may not use UPMNS)
      particleDataPtr->names(12,"nu_1","nu_1bar");
      particleDataPtr->names(14,"nu_2","nu_2bar");
      particleDataPtr->names(16,"nu_3","nu_3bar");
    }
    // Sneutrino flavour violation
    if ( (slha.modsel(6) == 2 || slha.modsel(6) >= 3)
         && slha.snumix.exists()) {
      // Modify sneutrino names
      particleDataPtr->names(1000012,"~nu_1","~nu_1bar");
      particleDataPtr->names(1000014,"~nu_2","~nu_2bar");
      particleDataPtr->names(1000016,"~nu_3","~nu_3bar");
    }
    // Optionally allow for separate scalar and pseudoscalar sneutrinos
    if ( slha.snsmix.exists() && slha.snamix.exists() ) {
      // Scalar sneutrinos
      particleDataPtr->names(1000012,"~nu_S1","~nu_S1bar");
      particleDataPtr->names(1000014,"~nu_S2","~nu_S2bar");
      particleDataPtr->names(1000016,"~nu_S3","~nu_S3bar");
      // Add the pseudoscalar sneutrinos
      particleDataPtr->addParticle(1000017, "~nu_A1", "~nu_A1bar",1, 0., 0);
      particleDataPtr->addParticle(1000018, "~nu_A2", "~nu_A2bar",1, 0., 0);
      particleDataPtr->addParticle(1000019, "~nu_A3", "~nu_A3bar",1, 0., 0);
    }
  }

  // SLHA2 spectrum with RPV
  if ( (ifailSpc == 1 || ifailSpc == 0) &&  slha.modsel(4) >= 1 ) {
    if ( slha.rvnmix.exists() ) {
      // Neutralinos -> neutrinos
      // Maintain R-conserving names since mass-ordering unlikely to change.
      particleDataPtr->names(12,"nu_1","nu_1bar");
      particleDataPtr->names(14,"nu_2","nu_2bar");
      particleDataPtr->names(16,"nu_3","nu_3bar");
      particleDataPtr->name(1000022,"~chi_10");
      particleDataPtr->name(1000023,"~chi_20");
      particleDataPtr->name(1000025,"~chi_30");
      particleDataPtr->name(1000035,"~chi_40");
    }
    if ( slha.rvumix.exists() && slha.rvvmix.exists() ) {
      // Charginos -> charged leptons (note sign convention)
      // Maintain R-conserving names since mass-ordering unlikely to change.
      particleDataPtr->names(11,"e-","e+");
      particleDataPtr->names(13,"mu-","mu+");
      particleDataPtr->names(15,"tau-","tau+");
      particleDataPtr->names(1000024,"~chi_1+","~chi_1-");
      particleDataPtr->names(1000037,"~chi_2+","~chi_2-");
    }
    if ( slha.rvhmix.exists() ) {
      // Sneutrinos -> higgses (general mass-ordered names)
      particleDataPtr->name(25,"H_10");
      particleDataPtr->name(35,"H_20");
      particleDataPtr->names(1000012,"H_30","H_30");
      particleDataPtr->names(1000014,"H_40","H_40");
      particleDataPtr->names(1000016,"H_50","H_50");
    }
    if ( slha.rvamix.exists() ) {
      // Sneutrinos -> higgses (general mass-ordered names)
      particleDataPtr->name(36,"A_10");
      particleDataPtr->names(1000017,"A_20","A_20");
      particleDataPtr->names(1000018,"A_30","A_30");
      particleDataPtr->names(1000019,"A_40","A_40");
    }
    if ( slha.rvlmix.exists() ) {
      // sleptons -> charged higgses (note sign convention)
      particleDataPtr->names(37,"H_1+","H_1-");
      particleDataPtr->names(1000011,"H_2-","H_2+");
      particleDataPtr->names(1000013,"H_3-","H_3+");
      particleDataPtr->names(1000015,"H_4-","H_4+");
      particleDataPtr->names(2000011,"H_5-","H_5+");
      particleDataPtr->names(2000013,"H_6-","H_6+");
      particleDataPtr->names(2000015,"H_7-","H_7+");
    }
  }

  // SLHA2 spectrum with CPV
  if ( (ifailSpc == 1 || ifailSpc == 0) &&  slha.modsel(5) >= 1 ) {
    // no scalar/pseudoscalar distinction
    particleDataPtr->name(25,"H_10");
    particleDataPtr->name(35,"H_20");
    particleDataPtr->name(36,"H_30");
  }


  // Import qnumbers
  vector<int> isQnumbers;
  bool foundLowCode = false;
  if ( (ifailSpc == 1 || ifailSpc == 0) && slha.qnumbers.size() > 0) {
    for (int iQnum=0; iQnum < int(slha.qnumbers.size()); iQnum++) {
      // Always use positive id codes
      int id = abs(slha.qnumbers[iQnum](0));
      ostringstream idCode;
      idCode << id;
      if (particleDataPtr->isParticle(id)) {
        infoPtr->errorMsg(warnPref + "ignoring QNUMBERS", "for id = "
          + idCode.str() + " (already exists)", true);
      } else {
        int qEM3    = slha.qnumbers[iQnum](1);
        int nSpins  = slha.qnumbers[iQnum](2);
        int colRep  = slha.qnumbers[iQnum](3);
        int hasAnti = slha.qnumbers[iQnum](4);
        // Translate colRep to PYTHIA colType
        int colType = 0;
        if (colRep == 3) colType = 1;
        else if (colRep == -3) colType = -1;
        else if (colRep == 8) colType = 2;
        else if (colRep == 6) colType = 3;
        else if (colRep == -6) colType = -3;
        // Default name: PDG code
        string name, antiName;
        ostringstream idStream;
        idStream<<id;
        name     = idStream.str();
        antiName = "-"+name;
        if (iQnum < int(slha.qnumbersName.size())) {
          name = slha.qnumbersName[iQnum];
          antiName = slha.qnumbersAntiName[iQnum];
          if (antiName == "") {
            if (name.find("+") != string::npos) {
              antiName = name;
              antiName.replace(antiName.find("+"),1,"-");
            } else if (name.find("-") != string::npos) {
              antiName = name;
              antiName.replace(antiName.find("-"),1,"+");
            } else {
              antiName = name+"bar";
            }
          }
        }
        if ( hasAnti == 0) {
          antiName = "";
          particleDataPtr->addParticle(id, name, nSpins, qEM3, colType);
        } else {
          particleDataPtr->addParticle(id, name, antiName, nSpins, qEM3,
            colType);
        }
        // Store list of particle codes added by QNUMBERS
        isQnumbers.push_back(id);
        if (id < 1000000) foundLowCode = true;
      }
    }
  }
  // Inform user that BSM particles should ideally be assigned id codes > 1M
  if (foundLowCode)
    infoPtr->errorMsg(warnPref
      + "using QNUMBERS for id codes < 1000000 may clash with SM.");

  // Import mass spectrum.
  bool   keepSM            = settings.flag("SLHA:keepSM");
  double minMassSM         = settings.parm("SLHA:minMassSM");
  bool   allowUserOverride = settings.flag("SLHA:allowUserOverride");
  vector<int> idModified;
  if (ifailSpc == 1 || ifailSpc == 0) {

    // Start at beginning of mass array
    int    id = slha.mass.first();
    // Loop through to update particle data.
    for (int i = 1; i <= slha.mass.size() ; i++) {
      // Step to next mass entry
      if (i>1) id = slha.mass.next();
      ostringstream idCode;
      idCode << id;
      double mass = abs(slha.mass(id));

      // Check if this ID was added by qnumbers
      bool isInternal = true;
      for (unsigned int iq = 0; iq<isQnumbers.size(); ++iq)
        if (id == isQnumbers[iq]) isInternal = false;

      // Ignore masses for known SM particles or particles with
      // default masses < minMassSM; overwrite masses for rest.
      if (keepSM && (id < 25 || (id > 80 && id < 1000000)) && isInternal)
        infoPtr->errorMsg(warnPref + "ignoring MASS entry", "for id = "
          + idCode.str()
          + " (SLHA:keepSM. Use id > 1000000 for new particles)", true);
      else if (id < 1000000 && particleDataPtr->m0(id) < minMassSM
        && isInternal) {
        infoPtr->errorMsg(warnPref + "ignoring MASS entry", "for id = "
          + idCode.str() + " (m0 < SLHA:minMassSM)", true);
      }

      // Also ignore SLHA mass values if user has already set
      // a different value and is allowed to override them.
      else if (allowUserOverride && particleDataPtr->hasChanged(id)) {
        ostringstream mValue;
        mValue << particleDataPtr->m0(id);
        infoPtr->errorMsg(warnPref + "keeping user mass",
          "for id = " + idCode.str() + ", m0 = " + mValue.str(), true);
        idModified.push_back(id);
      }
      else {
        particleDataPtr->m0(id,mass);
        idModified.push_back(id);
      }
    };

  }

  // Update decay data.
  for (int iTable=0; iTable < int(slha.decays.size()); iTable++) {

    // Pointer to this SLHA table
    LHdecayTable* slhaTable=&(slha.decays[iTable]);

    // Extract ID and create pointer to corresponding particle data object
    int idRes     = slhaTable->getId();
    ostringstream idCode;
    idCode << idRes;
    ParticleDataEntry* particlePtr
      = particleDataPtr->particleDataEntryPtr(idRes);

    // Check if this ID was added by qnumbers
    bool isInternal = true;
    for (unsigned int iq = 0; iq<isQnumbers.size(); ++iq)
      if (idRes == isQnumbers[iq]) isInternal = false;

    // Ignore decay channels for known SM particles or particles with
    // default masses < minMassSM; overwrite masses for rest.
    if (keepSM && (idRes < 25 || (idRes > 80 && idRes < 1000000))
        && isInternal) {
      infoPtr->errorMsg(warnPref + "ignoring DECAY table", "for id = "
        + idCode.str()
        + " (SLHA:keepSM. Use id > 1000000 for new particles)", true);
      continue;
    }
    else if (idRes < 1000000 && particleDataPtr->m0(idRes) < minMassSM
             && !particleDataPtr->hasChanged(idRes) && isInternal) {
      infoPtr->errorMsg(warnPref + "ignoring DECAY table", "for id = "
        + idCode.str() + " (m0 < SLHA:minMassSM)", true);
      continue;
    }

    // Verbose output. Let user know we are using these tables.
    if (verboseSLHA <= 1)
      infoPtr->errorMsg(infoPref + "importing SLHA decay table(s)","");
    else
      infoPtr->errorMsg(infoPref + "importing SLHA decay table","for id = "
        +idCode.str(),true);

    // Extract and store total width (absolute value, neg -> switch off)
    double widRes         = abs(slhaTable->getWidth());
    double pythiaMinWidth = settings.parm("ResonanceWidths:minWidth");
    if (widRes > 0. && widRes < pythiaMinWidth) {
      infoPtr->errorMsg(warnPref + "forcing width = 0 ","for id = "
        + idCode.str() + " (width < ResonanceWidths:minWidth)" , true);
      widRes = 0.0;
    }
    particlePtr->setMWidth(widRes);

    // Set lifetime in mm for displaced vertex calculations
    // (convert GeV^-1 to mm)
    if (widRes > 0.) {
      double decayLength = 1.97e-13/widRes;
      particlePtr->setTau0(decayLength);

      // Reset decay table of the particle. Allow decays, treat as resonance.
      if (slhaTable->size() > 0) {
        particlePtr->clearChannels();
        particleDataPtr->mayDecay(idRes,true);
        particleDataPtr->isResonance(idRes,true);
      } else {
        infoPtr->errorMsg(warnPref + "empty DECAY table ","for id = "
          + idCode.str() + " (total width provided but no branching"
          + " fractions)", true);
      }
    }
    // Reset to stable if width <= 0.0
    else {
      particlePtr->clearChannels();
      particleDataPtr->mayDecay(idRes,false);
      particleDataPtr->isResonance(idRes,false);
    }

    // Set initial minimum mass.
    double brWTsum   = 0.;
    double massWTsum = 0.;

    // Loop over SLHA channels, import into Pythia, treating channels
    // with negative branching fractions as having the equivalent positive
    // branching fraction, but being switched off for this run
    for (int iChannel=0 ; iChannel<slhaTable->size(); iChannel++) {
      LHdecayChannel slhaChannel = slhaTable->getChannel(iChannel);
      double brat      = slhaChannel.getBrat();
      vector<int> idDa = slhaChannel.getIdDa();
      if (idDa.size() >= 9) {
        infoPtr->errorMsg(errPref + "max number of DECAY products is 8");
      } else if (idDa.size() <= 1) {
        infoPtr->errorMsg(errPref + "min number of DECAY products is 2");
      }
      else {
        int onMode = 1;
        if (brat < 0.0) onMode = 0;
        int meModeNow = meMode;

        // Check phase space, including margin ~ sqrt(sum(widths^2))
        double massSum(0.);
        double widSqSum = pow2(widRes);
        int nDa = idDa.size();
        for (int jDa=0; jDa<nDa; ++jDa) {
          massSum  += particleDataPtr->m0( idDa[jDa] );
          widSqSum +=  pow2(particleDataPtr->mWidth( idDa[jDa] ));
        }
        double deltaM = particleDataPtr->m0(idRes) - massSum;
        // Negative mass difference: intrinsically off shell
        if (onMode == 1 && brat > 0.0 && deltaM < 0.) {
          // String containing decay name
          ostringstream errCode;
          errCode << idRes <<" ->";
          for (int jDa=0; jDa<nDa; ++jDa) errCode<<" "<<idDa[jDa];
          // Could mass fluctuations at all give the needed deltaM ?
          if (abs(deltaM) > 100. * sqrt(widSqSum)) {
            infoPtr->errorMsg(warnPref + "switched off DECAY mode",
              ": " + errCode.str()+" (too far off shell)",true);
            onMode = 0;
          }
          // If ~ OK within fluctuations
          else {
            // Ignore user-selected meMode
            if (meModeNow != 100) {
              infoPtr->errorMsg(warnPref + "adding off shell DECAY mode",
                ": "+errCode.str()+" (forced meMode = 100)",true);
              meModeNow = 100;
            } else {
              infoPtr->errorMsg(warnPref + "adding off shell DECAY mode",
                errCode.str(), true);
            }
          }
        }
        // Branching-ratio-weighted average mass in decay.
        brWTsum   += abs(brat);
        massWTsum += abs(brat) * massSum;

        // Add channel
        int id0 = idDa[0];
        int id1 = idDa[1];
        int id2 = (idDa.size() >= 3) ? idDa[2] : 0;
        int id3 = (idDa.size() >= 4) ? idDa[3] : 0;
        int id4 = (idDa.size() >= 5) ? idDa[4] : 0;
        int id5 = (idDa.size() >= 6) ? idDa[5] : 0;
        int id6 = (idDa.size() >= 7) ? idDa[6] : 0;
        int id7 = (idDa.size() >= 8) ? idDa[7] : 0;
        particlePtr->addChannel(onMode,abs(brat),meModeNow,
                                id0,id1,id2,id3,id4,id5,id6,id7);

      }
    }

    // Set minimal mass, but always below nominal one.
    if (slhaTable->size() > 0) {
      double massAvg = massWTsum / brWTsum;
      double massMin = min( massAvg, particlePtr->m0()) ;
      particlePtr->setMMin(massMin);
    }

    // Add to list of particles that have been modified
    idModified.push_back(idRes);

  }

  // Sanity check of all decay tables with modified MASS or DECAY info
  for (int iMod = 0; iMod < int(idModified.size()); ++iMod) {
    int id = idModified[iMod];
    ostringstream idCode;
    idCode << id;
    ParticleDataEntry* particlePtr
      = particleDataPtr->particleDataEntryPtr(id);
    double m0  = particlePtr->m0();
    double wid = particlePtr->mWidth();
    // Always set massless particles stable
    if (m0 <= 0.0 && (wid > 0.0 || particlePtr->mayDecay())) {
      infoPtr->errorMsg(warnPref + "massless particle forced stable"," id = "
        + idCode.str(), true);
      particlePtr->clearChannels();
      particlePtr->setMWidth(0.0);
      particlePtr->setMayDecay(false);
      particleDataPtr->isResonance(id,false);
      continue;
    }
    // Declare zero-width particles to be stable (for now)
    if (wid == 0.0 && particlePtr->mayDecay()) {
      particlePtr->setMayDecay(false);
      continue;
    }
    // Check at least one on-shell channel is available
    double mSumMin = 10. * m0;
    int nChannels = particlePtr->sizeChannels();
    if (nChannels >= 1) {
      for (int iChannel=0; iChannel<nChannels; ++iChannel) {
        DecayChannel channel = particlePtr->channel(iChannel);
        if (channel.onMode() <= 0) continue;
        int nProd = channel.multiplicity();
        double mSum = 0.;
        for (int iDa = 0; iDa < nProd; ++iDa) {
          int idDa   = channel.product(iDa);
          mSum += particleDataPtr->m0(idDa);
        }
        mSumMin = min(mSumMin, mSum);
      }
      // Require at least one on-shell channel
      if (mSumMin > m0) {
        infoPtr->errorMsg(warnPref + "particle forced stable"," id = "
          + idCode.str() + " (no on-shell decay channels)", true);
        particlePtr->setMWidth(0.0);
        particlePtr->setMayDecay(false);
        continue;
      }
      else {
      // mMin: lower cutoff on Breit-Wigner: default is mMin = m0 - 5*Gamma
      // (User is allowed to specify a lower value if desired.)
      // Increase minimum if needed to ensure at least one channel on shell
        double mMin = min(particlePtr->mMin(), max(0.0,m0 - 5.*wid));
        mMin = max(mSumMin,mMin);
        particlePtr->setMMin(mMin);
      }
    }
  }

  return true;

}

//--------------------------------------------------------------------------

// Initialize SLHA blocks SMINPUTS and MASS from PYTHIA SM parameter values.
// E.g., to make sure that there are no important unfilled entries

void SLHAinterface::pythia2slha(ParticleData* particleDataPtr) {

  // Initialize block SMINPUTS.
  string blockName = "sminputs";
  double mZ = particleDataPtr->m0(23);
  slha.set(blockName,1,1.0/couplingsPtr->alphaEM(pow2(mZ)));
  slha.set(blockName,2,couplingsPtr->GF());
  slha.set(blockName,3,couplingsPtr->alphaS(pow2(mZ)));
  slha.set(blockName,4,mZ);
  // b mass (should be running mass, here pole for time being)
  slha.set(blockName,5,particleDataPtr->m0(5));
  slha.set(blockName,6,particleDataPtr->m0(6));
  slha.set(blockName,7,particleDataPtr->m0(15));
  slha.set(blockName,8,particleDataPtr->m0(16));
  slha.set(blockName,11,particleDataPtr->m0(11));
  slha.set(blockName,12,particleDataPtr->m0(12));
  slha.set(blockName,13,particleDataPtr->m0(13));
  slha.set(blockName,14,particleDataPtr->m0(14));
  // Force 3 lightest quarks massless
  slha.set(blockName,21,double(0.0));
  slha.set(blockName,22,double(0.0));
  slha.set(blockName,23,double(0.0));
  // c mass (should be running mass, here pole for time being)
  slha.set(blockName,24,particleDataPtr->m0(4));

  // Initialize block MASS.
  blockName = "mass";
  int id = 1;
  int count = 0;
  while (particleDataPtr->nextId(id) > id) {
    slha.set(blockName,id,particleDataPtr->m0(id));
    id = particleDataPtr->nextId(id);
    ++count;
    if (count > 10000) {
      infoPtr->errorMsg("Error in SLHAinterface::pythia2slha(): "
        "encountered infinite loop when saving mass block");
      break;
    }
  }

}

//==========================================================================

} // end namespace Pythia8



