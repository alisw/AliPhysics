//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed 
//      for the BaBar collaboration.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Module: EvtCheckDecays
//
// Description:  Holds code to conduct various checks on the 
//      EvtDecayTable::decaytable()
//
// Modification history:
//      Abi Soffer             Nov 29, 2007, created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGen/EvtCheckDecays.hh"

#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtId.hh"
#include <vector>
#include <string.h>
#include <iostream>
using namespace std;



void EvtCheckDecays::checkConj(bool compareArguments) {
  std::ostream & str = report(INFO, "EvtGen: checkConj") 
    << "Checking CP-conjugation of all decays:" << endl;

  const std::vector<EvtParticleDecayList> & decaytable
    = EvtDecayTable::decaytable();

  vector<unsigned int> skipEntries;

  // Loop over the entris (=particles) in the decay table:
  for(size_t e = 0; e < decaytable.size(); ++e){
    bool diffFound = false;

    // Skip entries that were already found to be conjugates of other
    // particles:
    bool skipThis = false;
    for (size_t s = 0; s < skipEntries.size(); ++s) {
      if (skipEntries[s] == e) {
	skipThis = true;
	break;
      }
    }
    if (skipThis) {
      continue;
    }

    // Start working with this particle's decays list:
    const EvtParticleDecayList & decList = decaytable[e];

    if (decList.getNMode() == 0) { // an undecaying entry in the table. ignore:
      continue;
    }

    // Get the decaying particle's ID (somewhat non-intuitively, this
    // is available from the individual decays, but not the decay list):
    const EvtId parentId = decList.getDecayModel(0)->getParentId();

    if (selfConj(parentId)) { // ignore self-conjugate particles:
      continue;
    }

    // Find the charge conjugate particle by looping over decaytable again:
    const EvtParticleDecayList * decList2 = 0;
    bool conjFound = false;

    // Need to create parentId2 out of the loop in which it gets filled witn
    // meaningful information. Since EVtId has no argumentless constructor,
    // initialize it to parentId. This won't cause logic problems, since
    // the variable conjFound is used to determine if parentId2 is meaningful:
    EvtId parentId2 = parentId;

    // the loop starts at e + 1, so as not to double count:
    for(size_t e2 = e + 1; e2 < EvtPDL::entries(); ++e2){ 
      decList2 = &(decaytable[e2]);
      
      if (decList2->getNMode() == 0) { // an undecaying entry. ignore:
	continue;
      }

      parentId2 = decList2->getDecayModel(0)->getParentId();
      if (parentId2.isConjugate(parentId)) { 
	conjFound = true;
	skipEntries.push_back(e2);
	break; // found the conjugate. decList2 is the conj of decList
      }
    }

    // Check if conjugate was found:
    if (false == conjFound) {
      report(INFO, "EvtGen: checkConj") 
	<< "No conjugate particle decays found for " << EvtPDL::name(parentId) 
	<< endl;
    }
    else { // conjugate found. 
      const std::string name  = EvtPDL::name(parentId);
      const std::string name2 = EvtPDL::name(parentId2);

      // Analyze the two decays and compare them:
      // First, compare the number of decay modes:
      if (decList.getNMode() != decList2->getNMode()) {
	report(INFO, "EvtGen: checkConj") 
	  << "*** Unequal numbers of decays: " 
	  << decList.getNMode() << " for " << name
	  << " and " << decList2->getNMode() << " for " << name2
	  << "." << endl;

	// Find the particle with the more decays:
	const EvtParticleDecayList * decListLonger = &decList;
	string nameLonger = name;
	string nameShorter = name2;
	if (decList.getNMode() < decList2->getNMode()) {
	  decListLonger = decList2;
	  nameLonger = name2;
	  nameShorter = name;
	}

	const int ndiff = abs(decList.getNMode() - decList2->getNMode());
	str << " The last " << ndiff  << " " << nameLonger
	    << " decays (the ones missing from " << nameShorter
	    << ", assuming" << endl
	    << " there are no other problems) are " << endl;

	for (int m = decListLonger->getNMode() - ndiff; 
	     m < decListLonger->getNMode(); 
	     ++m) {
	  const EvtDecayBase * dec  = decListLonger->getDecayModel(m);
	  str << "  ";
	  dec->printInfo();
	}
	
	diffFound = true;
      }      
      
      // Compare each decay mode. Expect them to be entered in the same order:
      for (int m=0; m < decList.getNMode() && m < decList2->getNMode(); ++m) {
	const EvtDecayBase * dec  = decList.getDecayModel(m);
	const EvtDecayBase * dec2 = decList2->getDecayModel(m);

	// Compare BRs:
	if (dec->getBranchingFraction() != dec2->getBranchingFraction()) {
	  report(INFO, "EvtGen: checkConj") 
	    << "*** Unequal BRs: " << endl << "  " 
	    << dec->getBranchingFraction() << " for decay ";
	  dec->printInfo();
	  str << "O and " << endl << "  "
	      << dec2->getBranchingFraction() << " for decay ";
	  dec2->printInfo();
	  str << endl;
	  diffFound = true;
	}
	
	// Compare number of daughters:
	if (dec->getNDaug() != dec2->getNDaug()) {
	  report(INFO, "EvtGen: checkConj") 
	    << "*** Unequal #'s of daughters: " << endl << "  " 
	    << dec->getNDaug() << " for decay ";
	  dec->printInfo();
	  str << " and " << endl << "  "
	      << dec2->getNDaug() << " for decay ";
	  dec2->printInfo();
	  str << endl;
	  diffFound = true;
	}


	// Compare daughter ID's. Here we allow for the daughter order 
	// not to match, since it's a common practice to write DECAY.DEC
	// without much regard to carefully conjugating the daughter order:

	vector<int> usedDtrIndices; // to avoid double-considering a dtr2
	for (int d = 0; d < dec->getNDaug(); ++d) { // loop on dtrs of dec
	  const EvtId dtr = dec->getDaug(d);
	  const bool sConj = selfConj(dtr);
	  
	  bool dtrMatchFound = false;
	  
	  for (int d2 = 0; d2 < dec2->getNDaug(); ++d2) { // on dtrs of dec2
	    // Skip if this index has been already found for this decay:
	    bool skipDtr = false;
	    for (size_t s = 0; s < usedDtrIndices.size(); ++s) {
	      if (usedDtrIndices[s] == d2) {
		skipDtr = true;
		break;
	      }
	      if (skipDtr) {
		continue; // to next value of d2
	      }
	    } // end loop on s

	    const EvtId dtr2 = dec2->getDaug(d2);
	    // See if there is a matching daughter: Either dtr is
	    // self-conjugate and we find an identical particle among
	    // dec2's daughters, or it is no self-conjugate and we
	    // find its conjugate:
	    if (dtr.isConjugate(dtr2)                && !sConj ||
		dec->getDaug(d) == dec2->getDaug(d2) &&  sConj) {
	      dtrMatchFound = true; 
	      usedDtrIndices.push_back(d2); // don't consider this d2 again
	      break; // out of d2 loop
	    }
	  } // end d2 loop 

	  if (false == dtrMatchFound) {
	    // Couldn't find a match for dtr among the daughters of dec2, so:
	    report(INFO, "EvtGen: checkConj") 
	      << "*** Daughter #" << d << " in decay"
	      << endl << "  "; 
	    dec->printInfo();
	    str << " has no conjugate in decay " << endl << "  ";
	    dec2->printInfo();
	    str << endl;
	    diffFound = true;
	  }
	}     

	// Compare models:
	if (dec->getModelName() != dec2->getModelName()) {
	  report(INFO, "EvtGen: checkConj") 
	    << "*** Unequal model names in decays: " << endl << "  ";
	  dec->printInfo();
	  str << " and " << endl << "  ";
	  dec2->printInfo();
	  str << endl;
	  diffFound = true;
	}

	// Compare numbers of arguments:
	if (dec->getNArg() != dec2->getNArg()) {
	  report(INFO, "EvtGen: checkConj") 
	    << "*** Unequal numbers of arguments: " << endl << "  "
	    << dec->getNArg() << " in decay ";	  
	  dec->printInfo();
	  str << " and " << endl << "  "
	      << dec2->getNArg() << " in decay ";
	  dec2->printInfo();
	    str << endl;
	    diffFound = true;
	}

	// Argument value comparison may not always be desired, since
	// the argument values aren't always stored by EvtDecayBase, causing
	// many printouts:
	if (compareArguments) {	
	  // Compare argument values. First, check if we can tell
	  // the numbers of arguments correctly:
	  if ( dec->getNArg() != dec->getNStoredArg() ||
	       dec2->getNArg() != dec2->getNStoredArg()) {
	    report(INFO, "EvtGen: checkConj") 
	      << "=== Not sure about number of arguments in decay" 
	      << endl << "  ";
	    dec->printInfo();
	    str << " and " << endl << "  ";
	    dec2->printInfo();
	    str << "  Argument value comparison may not be complete." << endl;
	  }
	  
	  // Now do the actual argument values comparison:
	  for (int a = 0; a < dec->getNArg() && a < dec2->getNArg()
		 && a < dec->getNStoredArg() && a < dec2->getNStoredArg(); 
	       ++a) {
	    if (dec->getStoredArg(a) != dec2->getStoredArg(a)) {
	      report(INFO, "EvtGen: checkConj") 
		<< "*** Unequal arguments #: " << a << ":" 
		<< endl << "  " 
		<< dec->getStoredArg(a) << " for decay ";
	      dec->printInfo();
	      str << " and " << endl << "  "
		  << dec2->getStoredArg(a) << " for decay ";
	      dec2->printInfo();
	      str << endl;
	      diffFound = true;
	    } // end if arguments unequal
	  } // end loop on arguments
	} // end if to compare arguments
      } // end loop on decays
    } // end if conjugate found

    if (diffFound) { // mark diff between particle decays
      report(INFO, "EvtGen: checkConj") 
	<< "-------------------------------------------------------" << endl;
    }

  } // end loop on particle entries
}      
      


bool EvtCheckDecays::selfConj(const EvtId & id) {
  string name = EvtPDL::name(id);
  if (name == "vpho" ||
      name == "gamma" ||
      name == "g" ||
      name == "specflav" ||
      name == "phasespa" ||
      name == "pi0" ||
      name == "pi(2S)0" ||
      name == "eta" ||
      name == "eta(2S)" ||
      name == "eta(1405)" ||
      name == "eta(1475)" ||
      name == "phi(1680)" ||
      name == "eta'" ||
      name == "rho0" ||
      name == "rho(2S)0" ||
      name == "rho(3S)0" ||
      name == "omega" ||
      name == "phi" ||
      name == "a_00" ||
      name == "f_0" ||
      name == "f'_0" ||
      name == "b_10" ||
      name == "h_1" ||
      name == "h'_1" ||
      name == "a_10" ||
      name == "f_1" ||
      name == "f'_1" ||
      name == "a_10" ||
      name == "a_20" ||
      name == "f_2" ||
      name == "f_0(1500)" ||
      name == "f'_2" ||
      name == "K_S0" ||
      name == "K_L0" ||
      name == "eta_c" ||
      name == "eta_c(2S)" ||
      name == "J/psi" ||
      name == "psi(2S)" ||
      name == "psi(3770)" ||
      name == "psi(4040)" ||
      name == "psi(4160)" ||
      name == "psi(4415)" ||
      name == "h_c" ||
      name == "chi_c0" ||
      name == "chi_c1" ||
      name == "chi_c2" ||
      name == "eta_b" ||
      name == "eta_b(2S)" ||
      name == "eta_b(3S)" ||
      name == "Upsilon" ||
      name == "Upsilon(2S)" ||
      name == "Upsilon(3S)" ||
      name == "Upsilon(4S)" ||
      name == "Upsilon(5S)" ||
      name == "h_b" ||
      name == "h_b(2P)" ||
      name == "h_b(3P)" ||
      name == "chi_b0" ||
      name == "chi_b1" ||
      name == "chi_b2" ||
      name == "chi_b0(2P)" ||
      name == "chi_b1(2P)" ||
      name == "chi_b2(2P)" ||
      name == "chi_b0(3P)" ||
      name == "chi_b1(3P)" ||
      name == "chi_b2(3P)" ||
      name == "eta_b2(1D)" ||
      name == "eta_b2(2D)" ||
      name == "Upsilon_1(1D)" ||
      name == "Upsilon_2(1D)" ||
      name == "Upsilon_3(1D)" ||
      name == "Upsilon_1(2D)" ||
      name == "Upsilon_2(2D)" ||
      name == "Upsilon_3(2D)" ||
      name == "Xu0" ||
      name == "sigma_0") {
    return true;
  }

  return false;
}
