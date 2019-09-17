// Info.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Info class.

#include "Pythia8/Info.h"
#include <limits>

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

void Info::list() const {

  // Header and beam info.
  cout << "\n --------  PYTHIA Info Listing  ------------------------"
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
    cout << " No process has been set; something must have gone wrong! \n"
         << "\n --------  End PYTHIA Info Listing  --------------------"
         << "----------------" << endl;
    return;
  }

  // Colliding parton info.
  if (isRes) {
    cout << " In 1: id = " << setw(4) << id1pdfSave[0] << ", x = "
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
    if (!matchIdX) cout << " Warning: above flavour/x info does not match"
         << " incoming partons in event!\n";
    cout << "\n";
  }

  // Process name and code.
  cout << ((isRes && !hasSubSave[0]) ? " Subprocess " : " Process ")
       << nameSave << " with code " << codeSave << " is 2 -> "
       << nFinalSave << ".\n";

  // Subprocess name and code for nondiffractive processes.
  if (hasSubSave[0])
    cout << " Subprocess " << nameSubSave[0] << " with code " << codeSubSave[0]
         << " is 2 -> " << nFinalSubSave[0] << ".\n";

  // Process-type-specific kinematics information.
  if ( isRes && nFinalSave == 1)
    cout << " It has sHat = " << setw(10) << sH[0] << ".\n";
  else if ( isRes && nFinalSave == 2)
    cout << " It has sHat = " << setw(10) << sH[0] << ",    tHat = "
         << setw(10) << tH[0] << ",    uHat = " << setw(10) << uH[0] << ",\n"
         << "       pTHat = " << setw(10) << pTH[0] << ",   m3Hat = "
         << setw(10) << m3H[0] << ",   m4Hat = " << setw(10) << m4H[0] << ",\n"
         << "    thetaHat = " << setw(10) << thetaH[0] << ",  phiHat = "
         << setw(10) << phiH[0] << ".\n";
    else if ( nFinalSave == 2)
    cout << " It has s = " << setw(10) << sH[0] << ",    t = " << setw(10)
         << tH[0] << ",    u = " << setw(10) << uH[0] << ",\n"
         << "       pT = " << setw(10) << pTH[0] << ",   m3 = " << setw(10)
         << m3H[0] << ",   m4 = " << setw(10) << m4H[0] << ",\n"
         << "    theta = " << setw(10) << thetaH[0] << ",  phi = " << setw(10)
         << phiH[0] << ".\n";
  else if ( isRes && nFinalSave == 3)
    cout << " It has sHat = " << setw(10) << sH[0] << ", <pTHat> = "
         << setw(10) << pTH[0] << ".\n";
  else if ( nFinalSave == 3)
    cout << " It has s = " << setw(10) << sH[0] << ",  t_A = " << setw(10)
         << tH[0] << ",  t_B = " << setw(10) << uH[0] << ",\n"
         << "     <pT> = " << setw(10) << pTH[0] << ".\n";

  // Couplings.
  if (isRes) cout << "     alphaEM = " << setw(10) << alphaEMSave[0]
    << ",  alphaS = " << setw(10) << alphaSSave[0] << "    at Q2 = "
    << setw(10) << Q2RenSave[0] << ".\n";

  // Diffractive subsystems.
  for (int iDS = 1; iDS < 4; ++iDS) if (id1Save[iDS] != 0) {
    if (iDS == 1) cout << "\n Diffractive system on side A: \n";
    if (iDS == 2) cout << "\n Diffractive system on side B: \n";
    if (iDS == 3) cout << "\n Central diffractive system: \n";
    cout << " In 1: id = " << setw(4) << id1pdfSave[iDS] << ", x = "
         << setw(10) << x1pdfSave[iDS] << ", pdf = " << setw(10)
         << pdf1Save[iDS] << " at Q2 = " << setw(10) << Q2FacSave[iDS]
         << ".\n" << " In 2: id = " << setw(4) << id2pdfSave[iDS]
         << ", x = " << setw(10) << x2pdfSave[iDS] << ", pdf = "
         << setw(10) << pdf2Save[iDS] << " at same Q2.\n";
    cout << " Subprocess " << nameSubSave[iDS] << " with code "
         << codeSubSave[iDS] << " is 2 -> " << nFinalSubSave[iDS] << ".\n";
    if (nFinalSubSave[iDS] == 1) {
      cout << " It has sHat = " << setw(10) << sH[iDS] << ".\n";
    } else if (nFinalSubSave[iDS] == 2) {
      cout << " It has sHat = " << setw(10) << sH[iDS] << ",    tHat = "
           << setw(10) << tH[iDS] << ",    uHat = " << setw(10) << uH[iDS]
           << ",\n" << "       pTHat = " << setw(10) << pTH[iDS]
           << ",   m3Hat = " << setw(10) << m3H[iDS] << ",   m4Hat = "
           << setw(10) << m4H[iDS] << ",\n" << "    thetaHat = " << setw(10)
           << thetaH[iDS] << ",  phiHat = "  << setw(10) << phiH[iDS] << ".\n";
    }
    cout << "     alphaEM = " << setw(10) << alphaEMSave[iDS]
         << ",  alphaS = " << setw(10) << alphaSSave[iDS] << "    at Q2 = "
         << setw(10) << Q2RenSave[iDS] << ".\n";
  }

  // Impact parameter.
  if (bIsSet) cout << "\n Impact parameter b = " << setw(10) << bMPISave
    << " gives enhancement factor = " << setw(10) << enhanceMPISave
    << ".\n";

  // Multiparton interactions and shower evolution.
  if (evolIsSet) cout << " Max pT scale for MPI = " << setw(10) << pTmaxMPISave
    << ", ISR = " << setw(10) << pTmaxISRSave << ", FSR = " << setw(10)
    << pTmaxISRSave << ".\n Number of MPI = " << setw(5) << nMPISave
    << ", ISR = " << setw(5) << nISRSave << ", FSRproc = " << setw(5)
    << nFSRinProcSave << ", FSRreson = " << setw(5) << nFSRinResSave
    << ".\n";

  // Listing finished.
  cout << "\n --------  End PYTHIA Info Listing  --------------------"
       << "----------------" << endl;

}

//--------------------------------------------------------------------------

// Event weights and accumulated weight.

double Info::weight(int iWeight) const {
  double wt = (iWeight <= 0 || iWeight >= int(weightSave.size()))
    ? weightSave[0] : weightSave[iWeight];
  return (abs(lhaStrategySave) == 4) ? CONVERTMB2PB * wt : wt;
}

double Info::weightSum() const {
  return (abs(lhaStrategySave) == 4) ? CONVERTMB2PB * wtAccSum : wtAccSum;
}
//--------------------------------------------------------------------------

// Uncertainty variations initialization

void Info::initUncertainties(vector<string>* variationListIn, bool isISR) {
  size_t vNames = weightLabelSave.size();
  externalVariations.clear();
  externalVarNames.clear();
  externalGroupNames.clear();
  externalMap.clear();
  initialNameSave.clear();
  externalVariations.push_back("Baseline");
  initialNameSave.push_back("Baseline");
  for(vector<string>::const_iterator v=variationListIn->begin();
      v != variationListIn->end(); ++v) {
    string line = *v;
    // Remove initial blank spaces
    while (line.find(" ") == 0) line.erase( 0, 1);
    size_t pos=0;
    // Search for pdf:family keyword for SpaceShower
    if( isISR && ((pos = line.find("isr:pdf:family")) != string::npos) ) {
      size_t posEnd = line.find(" ",pos);
      if( posEnd == string::npos ) posEnd = line.size();
      for(size_t i=0; i < vNames; ++i) {
        string local = weightLabelSave[i];
        if( local.find("isr:pdf:member") != string::npos ) {
          size_t iEqual = local.find("=")+1;
          string nMember = local.substr(iEqual,local.size());
          nMember.append(" ");
          string tmpLine = line;
          tmpLine.replace(pos,posEnd-pos,local);
          size_t iBlank = line.find_first_of(" ");
          tmpLine.replace(iBlank,1,nMember);
          externalVariations.push_back(tmpLine);
          initialNameSave.push_back(line.substr(0,line.find_first_of(" ")));
        }
      }
    } else {
      externalVariations.push_back(line);
      initialNameSave.push_back(line.substr(0,line.find_first_of(" ")));
    }
  }
  externalVariationsSize = externalVariations.size();
  size_t nNames = externalVariationsSize;
  externalVarNames.resize(nNames);
  externalVarNames[0].push_back("Baseline");
  externalGroupNames.resize(nNames);
  externalGroupNames[0]="Baseline";
  for(size_t iWeight = 0; iWeight < nNames; ++iWeight) {
    string uVarString = toLower(externalVariations[iWeight]);
    size_t firstBlank  = uVarString.find_first_of(" ");
    size_t endLine = uVarString.size();
    if( firstBlank > endLine) continue;
    externalGroupNames[iWeight] = uVarString.substr(0,firstBlank);
    uVarString  = uVarString.substr(firstBlank+1,endLine);
    size_t pos = 0;
    while ((pos = uVarString.find(" ")) != string::npos) {
      string token = uVarString.substr(0, pos);
      externalVarNames[iWeight].push_back(token);
      uVarString.erase(0, pos + 1);
    }
    if (uVarString == "" || uVarString == " ") continue;
    externalVarNames[iWeight].push_back(uVarString);
  }
  externalMap.resize(nNames);
  for(size_t iWeight = 0; iWeight < vNames; ++iWeight) {
    for(size_t iV = 0; iV < nNames; ++iV) {
      for(size_t iW = 0; iW < externalVarNames[iV].size(); ++iW) {
        if( externalVarNames[iV][iW] == weightLabelSave[iWeight] ) {
          externalMap[iV].push_back(iWeight);
        } else if ( isISR && externalVarNames[iV][iW].find("isr:pdf:family")
        != string::npos && weightLabelSave[iWeight].find("isr:pdf:member")
        != string::npos ) {
          externalMap[iV].push_back(iWeight);
        }
      }
    }
  }
}

//--------------------------------------------------------------------------

// List of all hard processes switched on.

vector<int> Info::codesHard() {
  vector<int> codesNow;
  for (map<int, long>::iterator nTryEntry = nTryM.begin();
    nTryEntry != nTryM.end(); ++nTryEntry)
      codesNow.push_back( nTryEntry->first );
  return codesNow;
}

//--------------------------------------------------------------------------

// Print a message the first few times. Insert in database.

  void Info::errorMsg(string messageIn, string extraIn, bool showAlways) {

  // Recover number of times message occured. Also inserts new string.
  int times = messages[messageIn];
  ++messages[messageIn];

  // Print message the first few times.
  if (times < TIMESTOPRINT || showAlways) cout << " PYTHIA "
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

void Info::errorStatistics() {

  // Header.
  cout << "\n *-------  PYTHIA Error and Warning Messages Statistics  "
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
    cout << " |      0   no errors or warnings to report              "
         << "                                                          | \n";
  while (messageEntry != messages.end()) {
    // Message printout.
    string temp = messageEntry->first;
    int len = temp.length();
    temp.insert( len, max(0, 102 - len), ' ');
    cout << " | " << setw(6) << messageEntry->second << "   "
         << temp << " | \n";
    ++messageEntry;
  }

  // Done.
  cout << " |                                                       "
       << "                                                          | \n"
       << " *-------  End PYTHIA Error and Warning Messages Statistics"
       << "  ------------------------------------------------------* "
       << endl;

}

//--------------------------------------------------------------------------

// Return a list of all header key names

vector < string > Info::headerKeys() {
  vector < string > keys;
  for (map < string, string >::iterator it = headers.begin();
      it != headers.end(); it++)
    keys.push_back(it->first);
  return keys;
}

//--------------------------------------------------------------------------

// Reset the LHEF3 objects read from the init and header blocks.

void Info::setLHEF3InitInfo() {
  initrwgt     = 0;
  generators   = 0;
  weightgroups = 0;
  init_weights = 0;
  headerBlock  = "";
}

//--------------------------------------------------------------------------

// Set the LHEF3 objects read from the init and header blocks.

void Info::setLHEF3InitInfo( int LHEFversionIn, LHAinitrwgt *initrwgtIn,
  vector<LHAgenerator> *generatorsIn,
  map<string,LHAweightgroup> *weightgroupsIn,
  map<string,LHAweight> *init_weightsIn, string headerBlockIn ) {
  LHEFversionSave = LHEFversionIn;
  initrwgt        = initrwgtIn;
  generators      = generatorsIn;
  weightgroups    = weightgroupsIn;
  init_weights    = init_weightsIn;
  headerBlock     = headerBlockIn;
}

//--------------------------------------------------------------------------

// Reset the LHEF3 objects read from the event block.

void Info::setLHEF3EventInfo() {
  eventAttributes    = 0;
  weights_detailed   = 0;
  weights_compressed = 0;
  scales             = 0;
  weights            = 0;
  rwgt               = 0;
  weights_detailed_vector.resize(0);
  eventComments      = "";
  eventWeightLHEF    = 1.0;
}

//--------------------------------------------------------------------------

// Set the LHEF3 objects read from the event block.

void Info::setLHEF3EventInfo( map<string, string> *eventAttributesIn,
   map<string,double> *weights_detailedIn,
   vector<double> *weights_compressedIn,
   LHAscales *scalesIn, LHAweights *weightsIn,
   LHArwgt *rwgtIn, vector<double> weights_detailed_vecIn,
   string eventCommentsIn, double eventWeightLHEFIn ) {
   eventAttributes    = eventAttributesIn;
   weights_detailed   = weights_detailedIn;
   weights_compressed = weights_compressedIn;
   scales             = scalesIn;
   weights            = weightsIn;
   rwgt               = rwgtIn;
   weights_detailed_vector = weights_detailed_vecIn;
   eventComments      = eventCommentsIn;
   eventWeightLHEF    = eventWeightLHEFIn;
}

//--------------------------------------------------------------------------

// Retrieve events tag information.

string Info::getEventAttribute(string key, bool doRemoveWhitespace) {
  if (!eventAttributes) return "";
  if ( eventAttributes->find(key) != eventAttributes->end() ) {
    string res = (*eventAttributes)[key];
    if (doRemoveWhitespace)
      res.erase (remove (res.begin(), res.end(), ' '), res.end());
    return res;
  }
  return "";
}

//--------------------------------------------------------------------------

// Retrieve LHEF version.

int Info::LHEFversion() { return LHEFversionSave;}

//--------------------------------------------------------------------------

// Retrieve initrwgt tag information.

unsigned int Info::getInitrwgtSize() {
  if (!initrwgt) return 0;
  return initrwgt->weights.size();
}

//--------------------------------------------------------------------------

// Retrieve generator tag information.

unsigned int Info::getGeneratorSize() {
  if (!generators) return 0;
  return generators->size();
}

string Info::getGeneratorValue(unsigned int n) {
  if (!generators || generators->size() < n+1) return "";
  return (*generators)[n].contents;
}

string Info::getGeneratorAttribute( unsigned int n, string key,
  bool doRemoveWhitespace) {
  if (!generators || generators->size() < n+1) return "";
  string res("");
  if ( key == "name") {
    res = (*generators)[n].name;
  } else if ( key == "version") {
    res = (*generators)[n].version;
  } else if ( (*generators)[n].attributes.find(key)
           != (*generators)[n].attributes.end() ) {
    res = (*generators)[n].attributes[key];
  }
  if (doRemoveWhitespace && res != "")
    res.erase (remove (res.begin(), res.end(), ' '), res.end());
  return res;
}

//--------------------------------------------------------------------------

// Retrieve rwgt tag information.

unsigned int Info::getWeightsDetailedSize() {
  if (!weights_detailed) return 0;
  return weights_detailed->size();
}

double Info::getWeightsDetailedValue(string n) {
  if (weights_detailed->empty()
    || weights_detailed->find(n) == weights_detailed->end())
    return std::numeric_limits<double>::quiet_NaN();
  return (*weights_detailed)[n];
}

string Info::getWeightsDetailedAttribute(string n, string key,
  bool doRemoveWhitespace) {
  if (!rwgt || rwgt->wgts.find(n) == rwgt->wgts.end())
    return "";
  string res("");
  if ( key == "id") {
    res = rwgt->wgts[n].id;
  } else if ( rwgt->wgts[n].attributes.find(key)
           != rwgt->wgts[n].attributes.end() ) {
    res = rwgt->wgts[n].attributes[key];
  }
  if (doRemoveWhitespace && res != "")
    res.erase (remove (res.begin(), res.end(), ' '), res.end());
  return res;
}

//--------------------------------------------------------------------------

// Retrieve weights tag information.

unsigned int Info::getWeightsCompressedSize() {
  if (!weights_compressed) return 0;
  return weights_compressed->size();
}

double Info::getWeightsCompressedValue(unsigned int n) {
  if (weights_compressed->empty() || weights_compressed->size() < n+1)
    return std::numeric_limits<double>::quiet_NaN();
  return (*weights_compressed)[n];
}

string Info::getWeightsCompressedAttribute(string key,
  bool doRemoveWhitespace) {
  if (!weights || weights->attributes.find(key) == weights->attributes.end())
    return "";
  string res("");
  if ( weights->attributes.find(key)
           != weights->attributes.end() ) {
    res = weights->attributes[key];
  }
  if (doRemoveWhitespace && res != "")
    res.erase (remove (res.begin(), res.end(), ' '), res.end());
  return res;
}

//--------------------------------------------------------------------------

// Retrieve scales tag information.

string Info::getScalesValue(bool doRemoveWhitespace) {
  if (!scales) return "";
  string res = scales->contents;
  if (doRemoveWhitespace && res != "")
    res.erase (remove (res.begin(), res.end(), ' '), res.end());
  return res;
}

double Info::getScalesAttribute(string key) {
  if (!scales) return std::numeric_limits<double>::quiet_NaN();
  double res = std::numeric_limits<double>::quiet_NaN();
  if ( key == "muf") {
    res = scales->muf;
  } else if ( key == "mur") {
    res = scales->mur;
  } else if ( key == "mups") {
    res = scales->mups;
  } else if ( key == "SCALUP") {
    res = scales->SCALUP;
  } else if ( scales->attributes.find(key)
           != scales->attributes.end() ) {
    res = scales->attributes[key];
  }
  return res;
}

//==========================================================================

} // end namespace Pythia8
