// Info.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Info class.

#include "Info.h"

namespace Pythia8 {

//==========================================================================

// Info class.
// This class contains a mixed bag of information on the event generation 
// activity, especially on the current subprocess properties.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of times the same error message will be repeated at most.
const int Info::TIMESTOPRINT = 1; 

// LHA convention with cross section in pb may require conversion from mb.
const double Info::CONVERTMB2PB = 1e9;

//--------------------------------------------------------------------------

// List (almost) all information currently set.

void Info::list(ostream& os) const {

  // Header and beam info.
  os << "\n --------  PYTHIA Info Listing  ------------------------"
     << "---------------- \n \n" 
     << scientific << setprecision(3) 
     << " Beam A: id = " << setw(6) << idASave << ", pz = " << setw(10) 
     << pzASave << ", e = " << setw(10) << eASave << ", m = " << setw(10) 
     << mASave << ".\n"
     << " Beam B: id = " << setw(6) << idBSave << ", pz = " << setw(10) 
     << pzBSave << ", e = " << setw(10) << eBSave << ", m = " << setw(10) 
     << mBSave << ".\n\n";

  // Done if no subprocess has been defined.
  if (codeSave == 0 && nFinalSave == 0) {
    os << " No process has been set; something must have gone wrong! \n"
       << "\n --------  End PYTHIA Info Listing  --------------------"
       << "----------------" << endl; 
    return;
  }

  // Colliding parton info.
  if (isRes) {
    os << " In 1: id = " << setw(4) << id1pdfSave[0] << ", x = "
       << setw(10) << x1pdfSave[0] << ", pdf = " << setw(10) << pdf1Save[0] 
       << " at Q2 = " << setw(10) << Q2FacSave[0] << ".\n"  
       << " In 2: id = " << setw(4) << id2pdfSave[0] << ", x = " 
       << setw(10) << x2pdfSave[0] << ", pdf = " << setw(10) << pdf2Save[0] 
       << " at same Q2.\n";  
    bool matchIdX = true;
    if (id1pdfSave[0] != id1Save[0] || id2pdfSave[0] != id2Save[0]) 
      matchIdX = false;
    if (abs(x1pdfSave[0] - x1Save[0]) > 1e-4 * x1Save[0]) matchIdX = false;
    if (abs(x2pdfSave[0] - x2Save[0]) > 1e-4 * x2Save[0]) matchIdX = false;
    if (!matchIdX) os << " Warning: above flavour/x info does not match"
       << " incoming partons in event!\n";
    os << "\n";
  }

  // Process name and code.
  os << ((isRes && !hasSubSave[0]) ? " Subprocess " : " Process ") << nameSave 
     << " with code " << codeSave << " is 2 -> " << nFinalSave << ".\n";

  // Subprocess name and code for minimum bias processes.
  if (hasSubSave[0])
    os << " Subprocess " << nameSubSave[0] << " with code " << codeSubSave[0] 
       << " is 2 -> " << nFinalSubSave[0] << ".\n";

  // Process-type-specific kinematics information.
  if ( isRes && nFinalSave == 1) 
    os << " It has sHat = " << setw(10) << sH[0] << ".\n";  
  else if ( isRes && nFinalSave == 2)  
    os << " It has sHat = " << setw(10) << sH[0] << ",    tHat = " 
       << setw(10) << tH[0] << ",    uHat = " << setw(10) << uH[0] << ",\n"
       << "       pTHat = " << setw(10) << pTH[0] << ",   m3Hat = " 
       << setw(10) << m3H[0] << ",   m4Hat = " << setw(10) << m4H[0] << ",\n"
       << "    thetaHat = " << setw(10) << thetaH[0] << ",  phiHat = " 
       << setw(10) << phiH[0] << ".\n";  
  else if ( nFinalSave == 2)  
    os << " It has s = " << setw(10) << sH[0] << ",    t = " << setw(10) 
       << tH[0] << ",    u = " << setw(10) << uH[0] << ",\n"
       << "       pT = " << setw(10) << pTH[0] << ",   m3 = " << setw(10) 
       << m3H[0] << ",   m4 = " << setw(10) << m4H[0] << ",\n" 
       << "    theta = " << setw(10) << thetaH[0] << ",  phi = " << setw(10) 
       << phiH[0] << ".\n";
  else if ( isRes && nFinalSave == 3)  
    os << " It has sHat = " << setw(10) << sH[0] << ", <pTHat> = " 
       << setw(10) << pTH[0] << ".\n";  
  else if ( nFinalSave == 3)  
    os << " It has s = " << setw(10) << sH[0] << ",  t_A = " << setw(10) 
       << tH[0] << ",  t_B = " << setw(10) << uH[0] << ",\n" 
       << "     <pT> = " << setw(10) << pTH[0] << ".\n";    

  // Couplings.
  if (isRes) os << "     alphaEM = " << setw(10) << alphaEMSave[0] 
    << ",  alphaS = " << setw(10) << alphaSSave[0] << "    at Q2 = " 
    << setw(10) << Q2RenSave[0] << ".\n"; 

  // Diffractive subsystems.
  for (int iDS = 1; iDS < 4; ++iDS) if (id1Save[iDS] != 0) {
    if (iDS == 1) os << "\n Diffractive system on side A: \n";
    if (iDS == 2) os << "\n Diffractive system on side B: \n";
    if (iDS == 3) os << "\n Central diffractive system: \n";
    os << " In 1: id = " << setw(4) << id1pdfSave[iDS] << ", x = "
       << setw(10) << x1pdfSave[iDS] << ", pdf = " << setw(10) 
       << pdf1Save[iDS] << " at Q2 = " << setw(10) << Q2FacSave[iDS] 
       << ".\n" << " In 2: id = " << setw(4) << id2pdfSave[iDS] 
       << ", x = " << setw(10) << x2pdfSave[iDS] << ", pdf = " 
       << setw(10) << pdf2Save[iDS] << " at same Q2.\n";  
    os << " Subprocess " << nameSubSave[iDS] << " with code " 
       << codeSubSave[iDS] << " is 2 -> " << nFinalSubSave[iDS] << ".\n";
    if (nFinalSubSave[iDS] == 1) 
      os << " It has sHat = " << setw(10) << sH[iDS] << ".\n";  
    else if (nFinalSubSave[iDS] == 2)  
      os << " It has sHat = " << setw(10) << sH[iDS] << ",    tHat = " 
         << setw(10) << tH[iDS] << ",    uHat = " << setw(10) << uH[iDS] 
         << ",\n" << "       pTHat = " << setw(10) << pTH[iDS] 
         << ",   m3Hat = " << setw(10) << m3H[iDS] << ",   m4Hat = " 
         << setw(10) << m4H[iDS] << ",\n" << "    thetaHat = " << setw(10) 
         << thetaH[iDS] << ",  phiHat = "  << setw(10) << phiH[iDS] << ".\n";
      os << "     alphaEM = " << setw(10) << alphaEMSave[iDS] 
      << ",  alphaS = " << setw(10) << alphaSSave[iDS] << "    at Q2 = " 
      << setw(10) << Q2RenSave[iDS] << ".\n";  
  }

  // Impact parameter.
  if (bIsSet) os << "\n Impact parameter b = " << setw(10) << bMPISave 
    << " gives enhancement factor = " << setw(10) << enhanceMPISave 
    << ".\n";

  // Multiparton interactions and shower evolution.
  if (evolIsSet) os << " Max pT scale for MPI = " << setw(10) << pTmaxMPISave
    << ", ISR = " << setw(10) << pTmaxISRSave << ", FSR = " << setw(10) 
    << pTmaxISRSave << ".\n Number of MPI = " << setw(5) << nMPISave 
    << ", ISR = " << setw(5) << nISRSave << ", FSRproc = " << setw(5) 
    << nFSRinProcSave << ", FSRreson = " << setw(5) << nFSRinResSave 
    << ".\n"; 
       
  // Listing finished.
  os << "\n --------  End PYTHIA Info Listing  --------------------"
     << "----------------" << endl; 

}

//--------------------------------------------------------------------------

// Event weight and accumulated weight.

double Info::weight() const { return (abs(lhaStrategySave) == 4) 
  ? CONVERTMB2PB * weightSave : weightSave;
} 
  
double Info::weightSum() const {return (abs(lhaStrategySave) == 4) 
  ? CONVERTMB2PB * wtAccSum : wtAccSum;
}

//--------------------------------------------------------------------------
  
// Print a message the first few times. Insert in database.
 
  void Info::errorMsg(string messageIn, string extraIn, bool showAlways, 
    ostream& os) {
   
  // Recover number of times message occured. Also inserts new string.
  int times = messages[messageIn];
  ++messages[messageIn];

  // Print message the first few times.
  if (times < TIMESTOPRINT || showAlways) os << " PYTHIA " 
    << messageIn << " " << extraIn << endl;

}

//--------------------------------------------------------------------------

// Provide total number of errors/aborts/warnings experienced to date.

int Info::errorTotalNumber() {

  int nTot = 0;
  for ( map<string, int>::iterator messageEntry = messages.begin();
    messageEntry != messages.end(); ++messageEntry)
    nTot += messageEntry->second;
  return nTot;

}

//--------------------------------------------------------------------------

// Print statistics on errors/aborts/warnings.

void Info::errorStatistics(ostream& os) {

  // Header.
  os << "\n *-------  PYTHIA Error and Warning Messages Statistics  "
     << "----------------------------------------------------------* \n"
     << " |                                                       "
     << "                                                          | \n"
     << " |  times   message                                      "
     << "                                                          | \n" 
     << " |                                                       "
     << "                                                          | \n";

  // Loop over all messages
  map<string, int>::iterator messageEntry = messages.begin();
  if (messageEntry == messages.end()) 
    os << " |      0   no errors or warnings to report              "
       << "                                                          | \n";
  while (messageEntry != messages.end()) {
    // Message printout.
    string temp = messageEntry->first;
    int len = temp.length();
    temp.insert( len, max(0, 102 - len), ' ');
    os << " | " << setw(6) << messageEntry->second << "   " 
       << temp << " | \n";
    ++messageEntry;
  } 

  // Done. 
  os << " |                                                       "
     << "                                                          | \n"
     << " *-------  End PYTHIA Error and Warning Messages Statistics"
     << "  ------------------------------------------------------* " 
     << endl;

}

//==========================================================================

} // end namespace Pythia8
