// LHEF3.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file is written by Stefan Prestel.
// It contains the main class for LHEF 3.0 functionalities.
// Function definitions.

#include "Pythia8/LHEF3.h"

namespace Pythia8 {

//==========================================================================

// The XMLTag struct is used to represent all information within an XML tag.
// It contains the attributes as a map, any sub-tags as a vector of pointers
// to other XMLTag objects, and any other information as a single string.

//--------------------------------------------------------------------------

// Constants.
const XMLTag::pos_t XMLTag::end = string::npos;

//==========================================================================

// The LHAweights struct.

//--------------------------------------------------------------------------

// Construct from XML tag.

LHAweights::LHAweights(const XMLTag & tag) {
  for ( map<string,string>::const_iterator it = tag.attr.begin();
    it != tag.attr.end(); ++it ) {
    string v = it->second.c_str();
    attributes[it->first] = v;
  }

  contents = tag.contents;

  istringstream iss(tag.contents);
  double w;
  while ( iss >> w ) weights.push_back(w);
}

//--------------------------------------------------------------------------

// Print out.

void LHAweights::list(ostream & file) const {
  file << "<weights";
  for ( map<string,string>::const_iterator it = attributes.begin();
        it != attributes.end(); ++it )
    file << " " << it->first << "=\"" << it->second << "\"";
  file << ">";
  for ( int j = 0, M = weights.size(); j < M; ++j ) file << " " << weights[j];
  file << "</weights>" << endl;
}

//==========================================================================

// The LHAscales struct: Collect different scales relevant for an event.

//--------------------------------------------------------------------------

// Construct from an XML-tag.

LHAscales::LHAscales(const XMLTag & tag, double defscale)
  : muf(defscale), mur(defscale), mups(defscale), SCALUP(defscale) {
  for ( map<string,string>::const_iterator it = tag.attr.begin();
        it != tag.attr.end(); ++it ) {
    double v = atof(it->second.c_str());
    if ( it->first == "muf" ) muf = v;
    else if ( it->first == "mur" ) mur = v;
    else if ( it->first == "mups" ) mups = v;
    else attributes.insert(make_pair(it->first, v));
  }
  contents = tag.contents;
}

//--------------------------------------------------------------------------

// Print out the corresponding XML-tag.

void LHAscales::list(ostream & file) const {
  file << "<scales";
  file << " muf=\"" << muf << "\"";
  file << " mur=\"" << mur << "\"";
  file << " mups=\"" << mups << "\"";
  for ( map<string,double>::const_iterator it = attributes.begin();
        it != attributes.end(); ++it )
    file << " " << it->first << "=\"" << it->second << "\"";
  file << ">" << contents;
  file << "</scales>" << endl;
}

//==========================================================================

// The LHAgenerator struct: Collect generator information for an event file.

//--------------------------------------------------------------------------

// Construct from an XML-tag

LHAgenerator::LHAgenerator(const XMLTag & tag, string defname)
  : name(defname), version(defname), contents(defname) {
  for ( map<string,string>::const_iterator it = tag.attr.begin();
        it != tag.attr.end(); ++it ) {
    if ( it->first == "name" ) name = it->second;
    else if ( it->first == "version" ) version = it->second;
    else attributes.insert(make_pair(it->first, it->second));
  }
  contents = tag.contents;
}

//--------------------------------------------------------------------------

// Print out the corresponding XML-tag.

void LHAgenerator::list(ostream & file) const {
  file << "<generator";
  if ( name    != "" ) file << " name=\""    << name    << "\"";
  if ( version != "" ) file << " version=\"" << version << "\"";
  for ( map<string,string>::const_iterator it = attributes.begin();
        it != attributes.end(); ++it )
    file << " " << it->first << "=\"" << it->second << "\"";
  file << " >";
  file << contents;
  file << "</generator>" << endl;
}

//==========================================================================

// The LHAwgt struct: Collect the wgt information.

//--------------------------------------------------------------------------

// Construct from an XML-tag

LHAwgt::LHAwgt(const XMLTag & tag, double defwgt)
  : id(""), contents(defwgt) {
  for ( map<string,string>::const_iterator it = tag.attr.begin();
        it != tag.attr.end(); ++it ) {
    if ( it->first == "id" ) id = it->second;
    else attributes.insert(make_pair(it->first, it->second));
  }
  contents = atof(tag.contents.c_str());
}

//--------------------------------------------------------------------------

// Print out the corresponding XML-tag.

void LHAwgt::list(ostream & file) const {
  file << "<wgt";
  if ( id    != "" ) file << " id=\""    << id << "\"";
  for ( map<string,string>::const_iterator it = attributes.begin();
        it != attributes.end(); ++it )
    file << " " << it->first << "=\"" << it->second << "\"";
  file << " >";
  file << contents;
  file << "</wgt>" << endl;
}

//==========================================================================

// The LHAweight struct: Collect the weight information.

//--------------------------------------------------------------------------

// Construct from an XML-tag.

LHAweight::LHAweight(const XMLTag & tag, string defname)
  : id(defname), contents(defname) {
  for ( map<string,string>::const_iterator it = tag.attr.begin();
        it != tag.attr.end(); ++it ) {
    if ( it->first == "id" ) id = it->second;
    else attributes.insert(make_pair(it->first, it->second));
  }
  contents = tag.contents;
}

//--------------------------------------------------------------------------

// Print out the corresponding XML-tag.

void LHAweight::list(ostream & file) const {
  file << "<weight";
  if ( id  != "" ) file << " id=\""    << id << "\"";
  for ( map<string,string>::const_iterator it = attributes.begin();
        it != attributes.end(); ++it )
    file << " " << it->first << "=\"" << it->second << "\"";
  file << " >";
  file << contents;
  file << "</weight>" << endl;
}

//==========================================================================

// The LHAweightgroup struct: The LHAweightgroup assigns a group-name to a set
// of LHAweight objects.

//--------------------------------------------------------------------------

// Construct a group of LHAweight objects from an XML tag and
// insert them in the given vector.

LHAweightgroup::LHAweightgroup(const XMLTag & tag) {

  for ( map<string,string>::const_iterator it = tag.attr.begin();
        it != tag.attr.end(); ++it ) {
    if ( it->first == "name" ) name = it->second;
    else attributes.insert(make_pair(it->first,it->second));
  }
  if ( name=="" ) {
    string key("type");
    if( attributes.find(key) != attributes.end() ) {
      name = attributes[key];
    }
  }

  contents = tag.contents;

  // Now add the weight's step by step.
  string s;
  vector<XMLTag*> tags = XMLTag::findXMLTags(tag.contents, &s);
  for ( int i = 0, N = tags.size(); i < N; ++i ) {
    const XMLTag & tagnow = *tags[i];
    LHAweight wt(tagnow);
    weights.insert(make_pair(wt.id, wt));
    weightsKeys.push_back(wt.id);
  }
  for ( int i = 0, N = tag.tags.size(); i < N; ++i ) {
    const XMLTag & tagnow = *tag.tags[i];
    const LHAweight & wt(tagnow);
    weights.insert(make_pair(wt.id, wt));
    weightsKeys.push_back(wt.id);
  }

  for ( int i = 0, N = tags.size(); i < N; ++i ) if (tags[i]) delete tags[i];

}

//--------------------------------------------------------------------------

// Print out the corresponding XML-tag.

void LHAweightgroup::list(ostream & file) const {
  file << "<weightgroup";
  if ( name != "" ) file << " name=\"" << name << "\"";
  for ( map<string,string>::const_iterator it = attributes.begin();
        it != attributes.end(); ++it )
    file << " " << it->first << "=\"" << it->second << "\"";
  file << " >\n";
  for ( map<string,LHAweight>::const_iterator it = weights.begin();
        it != weights.end(); ++it ) it->second.list(file);
  file << "</weightgroup>" << endl;
}

//==========================================================================

// The LHArwgt struct: Assigns a group-name to a set of LHAwgt objects.

//--------------------------------------------------------------------------

// Construct a group of LHAwgt objects from an XML tag and
// insert them in the given vector.

LHArwgt::LHArwgt(const XMLTag & tag) {

  for ( map<string,string>::const_iterator it = tag.attr.begin();
        it != tag.attr.end(); ++it ) {
    string v = it->second.c_str();
    attributes[it->first] = v;
  }
  contents = tag.contents;

  // Now add the wgt's step by step.
  string s;
  vector<XMLTag*> tags = XMLTag::findXMLTags(tag.contents, &s);
  for ( int i = 0, N = tags.size(); i < N; ++i ) {
    const XMLTag & tagnow = *tags[i];
    LHAwgt wt(tagnow);
    wgts.insert(make_pair(wt.id, wt));
    wgtsKeys.push_back(wt.id);
  }
  for ( int i = 0, N = tag.tags.size(); i < N; ++i ) {
    const XMLTag & tagnow = *tag.tags[i];
    LHAwgt wt(tagnow);
    wgts.insert(make_pair(wt.id, wt));
    wgtsKeys.push_back(wt.id);
  }

  for ( int i = 0, N = tags.size(); i < N; ++i ) if (tags[i]) delete tags[i];

}

//--------------------------------------------------------------------------

// Print out the corresponding XML-tag.

void LHArwgt::list(ostream & file) const {
  file << "<rwgt";
  for ( map<string,string>::const_iterator it = attributes.begin();
        it != attributes.end(); ++it )
    file << " " << it->first << "=\"" << it->second << "\"";
  file << " >\n";
  for ( map<string,LHAwgt>::const_iterator it = wgts.begin();
        it != wgts.end(); ++it ) it->second.list(file);
  file << "</rwgt>" << endl;
}

//==========================================================================

// The LHAinitrwgt assigns a group-name to a set of LHAweightgroup objects.

//--------------------------------------------------------------------------

// Construct a group of LHAweightgroup objects from an XML tag and
// insert them in the given vector.

LHAinitrwgt::LHAinitrwgt(const XMLTag & tag) {
  for ( map<string,string>::const_iterator it = tag.attr.begin();
        it != tag.attr.end(); ++it ) {
    string v = it->second.c_str();
    attributes[it->first] = v;
  }
  contents = tag.contents;

  // Now add the wgt's step by step.
  string s;
  vector<XMLTag*> tags = XMLTag::findXMLTags(tag.contents, &s);
  for ( int i = 0, N = tags.size(); i < N; ++i ) {
    const XMLTag & tagnow = *tags[i];
    if ( tagnow.name == "weightgroup" ) {
      LHAweightgroup wgroup(tagnow);
      string wgname = wgroup.name;
      // if still no name, use integer as a key
      if (wgname=="") {
        stringstream iss;
        iss << i;
        wgname=iss.str();
      }
      weightgroups.insert(make_pair(wgname, wgroup));
      weightgroupsKeys.push_back(wgname);
      string ss;
      vector<XMLTag*> tags2 = XMLTag::findXMLTags(tagnow.contents, &ss);
      for ( int k = 0, M = tags2.size(); k < M; ++k ) {
        const XMLTag & tagnow2 = *tags2[k];
        if ( tagnow2.name == "weight" ) {
          LHAweight wt(tagnow2);
          string wtname = wt.id;
          weights.insert(make_pair(wtname, wt));
          weightsKeys.push_back(wtname);
        }
      }
      for ( int j = 0, M = tags2.size(); j < M; ++j )
        if (tags2[j]) delete tags2[j];
    } else if ( tagnow.name == "weight" ) {
      LHAweight wt(tagnow);
      string wtname = wt.id;
      weights.insert(make_pair(wtname, wt));
      weightsKeys.push_back(wtname);
    }
  }

  // Now add the wgt's step by step.
  for ( int i = 0, N = tag.tags.size(); i < N; ++i ) {
    const XMLTag & tagnow = *tag.tags[i];
    if ( tagnow.name == "weightgroup" ) {
      LHAweightgroup wgroup(tagnow);
      string wgname = wgroup.name;
      weightgroups.insert(make_pair(wgname, wgroup));
      weightgroupsKeys.push_back(wgname);
      string ss;
      vector<XMLTag*> tags2 = XMLTag::findXMLTags(tagnow.contents, &ss);
      for ( int k = 0, M = tags2.size(); k < M; ++k ) {
        const XMLTag & tagnow2 = *tags2[k];
        if ( tagnow2.name == "weight" ) {
          LHAweight wt(tagnow2);
          string wtname = wt.id;
          weights.insert(make_pair(wtname, wt));
          weightsKeys.push_back(wtname);
        }
      }
      for ( int k = 0, M = tagnow.tags.size(); k < M; ++k ) {
        const XMLTag & tagnow2 = *tagnow.tags[k];
        if ( tagnow2.name == "weight" ) {
          LHAweight wt(tagnow2);
          string wtname = wt.id;
          weights.insert(make_pair(wtname, wt));
          weightsKeys.push_back(wtname);
        }
      }
      for ( int j = 0, M = tags2.size(); j < M; ++j )
        if (tags2[j]) delete tags2[j];
    } else if ( tagnow.name == "weight" ) {
      LHAweight wt(tagnow);
      string wtname = wt.id;
      weights.insert(make_pair(wtname, wt));
      weightsKeys.push_back(wtname);
    }
  }

  for ( int i = 0, N = tags.size(); i < N; ++i ) if (tags[i]) delete tags[i];

}

//--------------------------------------------------------------------------

// Print out the corresponding XML-tag.

void LHAinitrwgt::list(ostream & file) const {
  file << "<initrwgt";
  for ( map<string,string>::const_iterator it = attributes.begin();
        it != attributes.end(); ++it )
    file << " " << it->first << "=\"" << it->second << "\"";
  file << " >\n";
  for ( map<string,LHAweightgroup>::const_iterator it = weightgroups.begin();
        it != weightgroups.end(); ++it ) it->second.list(file);
  for ( map<string,LHAweight>::const_iterator it = weights.begin();
        it != weights.end(); ++it ) it->second.list(file);
  file << "</initrwgt>" << endl;
}

//==========================================================================

// The HEPRUP class is a simple container for the Les Houches file init block.

void HEPRUP::clear() {
  IDBMUP = make_pair(0,0);
  EBMUP = make_pair(0,0);
  PDFGUP = make_pair(0,0);
  PDFSUP = make_pair(0,0);
  IDWTUP = -1;
  NPRUP = 0;
  XSECUP.resize(0);
  XERRUP.resize(0);
  XMAXUP.resize(0);
  LPRUP.resize(0);
  initrwgt.clear();
  generators.resize(0);
  weightgroups.clear();
  weights.clear();

}

//==========================================================================

// The HEPEUP class is a simple container corresponding to the Les Houches
// accord (<A HREF="http://arxiv.org/abs/hep-ph/0109068">hep-ph/0109068</A>)
// common block with the same name. The members are named in the same
// way as in the common block. However, fortran arrays are represented
// by vectors, except for the arrays of length two which are
// represented by pair objects.

//--------------------------------------------------------------------------

// Copy information from the given HEPEUP.

HEPEUP & HEPEUP::setEvent(const HEPEUP & x) {

  NUP                = x.NUP;
  IDPRUP             = x.IDPRUP;
  XWGTUP             = x.XWGTUP;
  XPDWUP             = x.XPDWUP;
  SCALUP             = x.SCALUP;
  AQEDUP             = x.AQEDUP;
  AQCDUP             = x.AQCDUP;
  IDUP               = x.IDUP;
  ISTUP              = x.ISTUP;
  MOTHUP             = x.MOTHUP;
  ICOLUP             = x.ICOLUP;
  PUP                = x.PUP;
  VTIMUP             = x.VTIMUP;
  SPINUP             = x.SPINUP;
  heprup             = x.heprup;
  scalesSave         = x.scalesSave;
  weightsSave        = x.weightsSave;
  weights_detailed   = x.weights_detailed;
  weights_compressed = x.weights_compressed;
  rwgtSave           = x.rwgtSave;
  attributes         = x.attributes;
  return *this;

}

//--------------------------------------------------------------------------

// Reset the HEPEUP object.

void HEPEUP::reset() {
  NUP = 0;
  weights_detailed.clear();
  weights_compressed.clear();
  weightsSave.clear();
  rwgtSave.clear();
  scalesSave.clear();
  attributes.clear();
}

//--------------------------------------------------------------------------

// Assuming the NUP variable, corresponding to the number of
// particles in the current event, is correctly set, resize the
// relevant vectors accordingly.

void HEPEUP::resize() {
  IDUP.resize(NUP);
  ISTUP.resize(NUP);
  MOTHUP.resize(NUP);
  ICOLUP.resize(NUP);
  PUP.resize(NUP, vector<double>(5));
  VTIMUP.resize(NUP);
  SPINUP.resize(NUP);
}

//==========================================================================

// The Reader class is initialized with a stream from which to read a
// version 1/2 Les Houches Accord event file. In the constructor of
// the Reader object the optional header information is read and then
// the mandatory init is read. After this the whole header block
// including the enclosing lines with tags are available in the public
// headerBlock member variable. Also the information from the init
// block is available in the heprup member variable and any additional
// comment lines are available in initComments. After each successful
// call to the readEvent() function the standard Les Houches Accord
// information about the event is available in the hepeup member
// variable and any additional comments in the eventComments
// variable. A typical reading sequence would look as follows:

//--------------------------------------------------------------------------

// Used internally in the constructors to read header and init blocks.
bool Reader::init() {

  bool readingHeader = false;
  bool readingInit = false;

  // Make sure we are reading a LHEF file:
  getLine();

  if ( currentLine.find("<LesHouchesEvents" ) == string::npos )
    return false;
  version = 0;
  if ( currentLine.find("version=\"1" ) != string::npos )
    version = 1;
  else if ( currentLine.find("version=\"2" ) != string::npos )
    version = 2;
  else if ( currentLine.find("version=\"3" ) != string::npos )
    version = 3;
  else
    return false;

  // Clear all members.
  outsideBlock="";
  headerBlock="";
  headerComments="";
  heprup.clear();
  initComments="";
  hepeup.clear();
  eventComments="";

  // Loop over all lines until we hit the </init> tag.
  while ( getLine() && currentLine.find("</init>") == string::npos ) {
    if ( currentLine.find("<header") != string::npos
      && currentLine.find("#") == string::npos) {
      // We have hit the header block, so we should dump this and
      // all following lines to headerBlock until we hit the end of
      // it.
      readingHeader = true;
      headerBlock = currentLine + "\n";
    }
    else if ( ( currentLine.find("<init>") != string::npos
      || currentLine.find("<init ") != string::npos )
      && currentLine.find("#") == string::npos) {
      // We have hit the init block, so we should expect to find the
      // standard information in the following.
      readingInit = true;

      // The first line tells us how many lines to read next.
      getLine();
      istringstream iss(currentLine);
      if ( !( iss >> heprup.IDBMUP.first >> heprup.IDBMUP.second
                  >> heprup.EBMUP.first >> heprup.EBMUP.second
                  >> heprup.PDFGUP.first >> heprup.PDFGUP.second
                  >> heprup.PDFSUP.first >> heprup.PDFSUP.second
                  >> heprup.IDWTUP >> heprup.NPRUP ) ) {
        heprup.NPRUP = -42;
        return false;
      }
      heprup.resize();

      for ( int i = 0; i < heprup.NPRUP; ++i ) {
        getLine();
        istringstream isss(currentLine);
        if ( !( isss >> heprup.XSECUP[i] >> heprup.XERRUP[i]
                    >> heprup.XMAXUP[i] >> heprup.LPRUP[i] ) ) {
          heprup.NPRUP = -42;
          return false;
        }
      }
    }
    else if ( currentLine.find("</header>") != string::npos
      && currentLine.find("#") == string::npos) {
      // The end of the header block. Dump this line as well to the
      // headerBlock and we're done.
      readingHeader = false;
      headerBlock += currentLine + "\n";
    }
    else if ( readingHeader ) {
      // We are in the process of reading the header block. Dump the
      // line to headerBlock.
      headerBlock += currentLine + "\n";
      headerComments += currentLine + "\n";
    }
    else if ( readingInit ) {
      // Here we found a comment line. Dump it to initComments.
      initComments += currentLine + "\n";
    }
    else {
      // We found some other stuff outside the standard tags.
      outsideBlock += currentLine + "\n";
    }
  }

  if ( file == NULL ) heprup.NPRUP = -42;

  // Scan the header block for XML tags
  string leftovers;
  vector<XMLTag*> tags1 = XMLTag::findXMLTags(headerComments, &leftovers);
  if ( leftovers.find_first_not_of(" \t\n") == string::npos )
    leftovers="";

  for ( int i = 0, N = tags1.size(); i < N; ++i ) {
    const XMLTag & tag = *tags1[i];

    if ( tag.name == "initrwgt" ) {
      LHAinitrwgt irwgt(tag);
      heprup.initrwgt = irwgt;
      for ( int j = 0, M = tag.tags.size(); j < M; ++j ) {
        XMLTag & ctag = *tag.tags[j];
        if ( ctag.name == "weightgroup" ) {
          LHAweightgroup wgroup(ctag);
          string wgname = wgroup.name;
          heprup.weightgroups.insert(make_pair(wgname, wgroup));

          string ss;
          vector<XMLTag*> tags2 = XMLTag::findXMLTags(ctag.contents, &ss);
          for ( int k = 0, O = tags2.size(); k < O; ++k ) {
            const XMLTag & tagnow2 = *tags2[k];
            if ( tagnow2.name == "weight" ) {
              LHAweight wt(tagnow2);
              string wtname = wt.id;
              heprup.weights.insert(make_pair(wtname, wt));
            }
          }
          for ( int k = 0, O = ctag.tags.size(); k < O; ++k ) {
            const XMLTag & tagnow2 = *ctag.tags[k];
            if ( tagnow2.name == "weight" ) {
              LHAweight wt(tagnow2);
              string wtname = wt.id;
              heprup.weights.insert(make_pair(wtname, wt));
            }
          }
        } else if ( ctag.name == "weight" ) {
          string tname = ctag.attr["id"];
          heprup.weights.insert(make_pair(tname, LHAweight(ctag)));
        }
      }
    }
  }

  heprup.generators.clear();
  // Scan the init block for XML tags
  leftovers="";
  vector<XMLTag*> tags2 = XMLTag::findXMLTags(initComments, &leftovers);
  if ( leftovers.find_first_not_of(" \t\n") == string::npos )
    leftovers="";

  for ( int i = 0, N = tags2.size(); i < N; ++i ) {
    const XMLTag & tag = *tags2[i];
    if ( tag.name == "generator" ) {
      heprup.generators.push_back(LHAgenerator(tag));
    }
  }

  for ( int i = 0, N = tags1.size(); i < N; ++i )
    if (tags1[i]) delete tags1[i];
  for ( int i = 0, N = tags2.size(); i < N; ++i )
    if (tags2[i]) delete tags2[i];

  // Done
  return true;

}

//--------------------------------------------------------------------------

// Read an event from the file and store it in the hepeup
// object. Optional comment lines are stored in the eventComments
// member variable. return true if the read was successful.

bool Reader::readEvent(HEPEUP * peup) {

  HEPEUP & eup = (peup? *peup: hepeup);
  eup.clear();
  eup.heprup = &heprup;
  weights_detailed_vec.clear();

  // Check if the initialization was successful. Otherwise we will
  // not read any events.
  if ( heprup.NPRUP < 0 ) return false;
  eventComments = "";
  outsideBlock = "";
  eup.NUP = 0;

  // Keep reading lines until we hit the next event or the end of
  // the event block. Save any inbetween lines. Exit if we didn't
  // find an event.
  while ( getLine() && currentLine.find("<event") == string::npos )
    outsideBlock += currentLine + "\n";

  // Get event attributes.
  if (currentLine != "") {
    string eventLine(currentLine);
    eventLine += "</event>";
    vector<XMLTag*> evtags = XMLTag::findXMLTags(eventLine);
    XMLTag & evtag = *evtags[0];
    for ( map<string,string>::const_iterator it = evtag.attr.begin();
          it != evtag.attr.end(); ++it ) {
      eup.attributes.insert(make_pair(it->first,it->second));
    }
    for ( int i = 0, N = evtags.size(); i < N; ++i )
      if (evtags[i]) delete evtags[i];
  }

  if ( !getLine()  ) return false;

  // We found an event. The first line determines how many
  // subsequent particle lines we have.
  istringstream iss(currentLine);
  if ( !( iss >> eup.NUP >> eup.IDPRUP >> eup.XWGTUP
              >> eup.SCALUP >> eup.AQEDUP >> eup.AQCDUP ) )
    return false;
  eup.resize();

  // Read all particle lines.
  for ( int i = 0; i < eup.NUP; ++i ) {
    if ( !getLine() ) return false;
    istringstream isss(currentLine);
    if ( !( isss >> eup.IDUP[i] >> eup.ISTUP[i]
                >> eup.MOTHUP[i].first >> eup.MOTHUP[i].second
                >> eup.ICOLUP[i].first >> eup.ICOLUP[i].second
                >> eup.PUP[i][0] >> eup.PUP[i][1] >> eup.PUP[i][2]
                >> eup.PUP[i][3] >> eup.PUP[i][4]
                >> eup.VTIMUP[i] >> eup.SPINUP[i] ) )
      return false;
  }

  // Now read any additional comments.
  while ( getLine() && currentLine.find("</event>") == string::npos )
    eventComments += currentLine + "\n";

  if ( file == NULL ) return false;

  eup.scalesSave = LHAscales(eup.SCALUP);

  // Scan the init block for XML tags
  string leftovers;
  vector<XMLTag*> tags = XMLTag::findXMLTags(eventComments, &leftovers);
  if ( leftovers.find_first_not_of(" \t\n") == string::npos )
    leftovers="";

  eventComments = "";
  istringstream f(leftovers);
  string l;
  while (getline(f, l)) {
     size_t p = l.find_first_not_of(" \t");
     l.erase(0, p);
     p = l.find_last_not_of(" \t");
     if (string::npos != p) l.erase(p+1);
     if (l.find_last_not_of("\n") != string::npos)
       eventComments += l + "\n";
  }

  for ( int i = 0, N = tags.size(); i < N; ++i ) {
    XMLTag & tag = *tags[i];

    if ( tag.name == "weights" ) {
      LHAweights wts(tag);
      eup.weightsSave = wts;

      for ( int k = 0, M = int(wts.weights.size()); k < M; ++k ) {
        eup.weights_compressed.push_back(wts.weights[k]);
      }

    }
    else if ( tag.name == "scales" ) {
      eup.scalesSave = LHAscales(tag, eup.SCALUP);
    }
    else if ( tag.name == "rwgt" ) {
      LHArwgt rwgt0(tag);
      eup.rwgtSave = rwgt0;
      string s;
      vector<XMLTag*> tags2 = XMLTag::findXMLTags(rwgt0.contents, &s);
      for ( int k = 0, M = tags2.size(); k < M; ++k ) {
        const XMLTag & tagnow = *tags2[k];
        if ( tagnow.name == "wgt" ) {
          LHAwgt wt(tagnow);
          eup.weights_detailed.insert(make_pair(wt.id, wt.contents));
          weights_detailed_vec.push_back(wt.contents);
        }
      }
      for ( int k = 0, M = tag.tags.size(); k < M; ++k ) {
        const XMLTag & tagnow = *tag.tags[k];
        if ( tagnow.name == "wgt" ) {
          LHAwgt wt(tagnow);
          eup.weights_detailed.insert(make_pair(wt.id, wt.contents));
          weights_detailed_vec.push_back(wt.contents);
        }
      }
    }
  }

  for ( int i = 0, N = tags.size(); i < N; ++i ) if (tags[i]) delete tags[i];

  return true;

}

//==========================================================================

// The Writer class is initialized with a stream to which to write a
// version 3.0 Les Houches Accord event file.

//--------------------------------------------------------------------------

// Write out an optional header block followed by the standard init
// block information together with any comment lines.

void Writer::init() {

  // Write out the standard XML tag for the event file.
  if ( version == 1 )
    file << "<LesHouchesEvents version=\"1.0\">" << endl;
  else
    file << "<LesHouchesEvents version=\"3.0\">" << endl;

  file << setprecision(8);

  // Print headercomments and header init information.
  file << "<header>" << endl;
  file << hashline(headerStream.str(),true) << std::flush;
  if ( version != 1 ) heprup.initrwgt.list(file);
  file << "</header>" << endl;

  file << "<init>"<< endl
       << " " << setw(8) << heprup.IDBMUP.first
       << " " << setw(8) << heprup.IDBMUP.second
       << " " << setw(14) << heprup.EBMUP.first
       << " " << setw(14) << heprup.EBMUP.second
       << " " << setw(4) << heprup.PDFGUP.first
       << " " << setw(4) << heprup.PDFGUP.second
       << " " << setw(4) << heprup.PDFSUP.first
       << " " << setw(4) << heprup.PDFSUP.second
       << " " << setw(4) << heprup.IDWTUP
       << " " << setw(4) << heprup.NPRUP << endl;
  heprup.resize();
  for ( int i = 0; i < heprup.NPRUP; ++i )
    file << " " << setw(14) << heprup.XSECUP[i]
         << " " << setw(14) << heprup.XERRUP[i]
         << " " << setw(14) << heprup.XMAXUP[i]
         << " " << setw(6) << heprup.LPRUP[i] << endl;

  if ( version == 1 ) {
    file << hashline(initStream.str(),true) << std::flush
         << "</init>" << endl;
    initStream.str("");
    return;
  }

  for ( int i = 0, N = heprup.generators.size(); i < N; ++i ) {
    heprup.generators[i].list(file);
  }

  file << hashline(initStream.str(),true) << std::flush
       << "</init>" << endl;
  initStream.str("");
}

//--------------------------------------------------------------------------

// Write out the event stored in hepeup, followed by optional
// comment lines.

bool Writer::writeEvent(HEPEUP * peup, int pDigits) {

  HEPEUP & eup = (peup? *peup: hepeup);

  file << "<event";
  for ( map<string,string>::const_iterator it = eup.attributes.begin();
        it != eup.attributes.end(); ++it )
    file << " " << it->first << "=\"" << it->second << "\"";
  file << ">" << std::flush << endl;
  file << " " << setw(4) << eup.NUP
       << " " << setw(6) << eup.IDPRUP
       << " " << setw(14) << eup.XWGTUP
       << " " << setw(14) << eup.SCALUP
       << " " << setw(14) << eup.AQEDUP
       << " " << setw(14) << eup.AQCDUP << endl;
  eup.resize();

  for ( int i = 0; i < eup.NUP; ++i )
    file << " " << setw(8) << eup.IDUP[i]
         << " " << setw(2) << eup.ISTUP[i]
         << " " << setw(4) << eup.MOTHUP[i].first
         << " " << setw(4) << eup.MOTHUP[i].second
         << " " << setw(4) << eup.ICOLUP[i].first
         << " " << setw(4) << eup.ICOLUP[i].second
         << " " << setw(pDigits) << eup.PUP[i][0]
         << " " << setw(pDigits) << eup.PUP[i][1]
         << " " << setw(pDigits) << eup.PUP[i][2]
         << " " << setw(pDigits) << eup.PUP[i][3]
         << " " << setw(pDigits) << eup.PUP[i][4]
         << " " << setw(1) << eup.VTIMUP[i]
         << " " << setw(1) << eup.SPINUP[i] << endl;

  // Write event comments.
  file << hashline(eventStream.str()) << std::flush;
  eventStream.str("");

  if ( version != 1 ) {
    eup.rwgtSave.list(file);
    eup.weightsSave.list(file);
    eup.scalesSave.list(file);
  }

  file << "</event>" << endl;

  if ( !file ) return false;

  return true;

}

//--------------------------------------------------------------------------

// Write out an event as a string.

string Writer::getEventString(HEPEUP * peup) {

  HEPEUP & eup = (peup? *peup: hepeup);

  stringstream helper;

  helper << "<event";
  for ( map<string,string>::const_iterator it = eup.attributes.begin();
        it != eup.attributes.end(); ++it )
    helper << " " << it->first << "=\"" << it->second << "\"";
  helper << ">" << std::flush << endl;
  helper << " " << setw(4) << eup.NUP
       << " " << setw(6) << eup.IDPRUP
       << " " << setw(14) << eup.XWGTUP
       << " " << setw(14) << eup.SCALUP
       << " " << setw(14) << eup.AQEDUP
       << " " << setw(14) << eup.AQCDUP << endl;
  eup.resize();

  for ( int i = 0; i < eup.NUP; ++i ) {
    helper << " " << setw(8) << eup.IDUP[i]
         << " " << setw(2) << eup.ISTUP[i]
         << " " << setw(4) << eup.MOTHUP[i].first
         << " " << setw(4) << eup.MOTHUP[i].second
         << " " << setw(6) << eup.ICOLUP[i].first
         << " " << setw(6) << eup.ICOLUP[i].second
         << fixed
         << setprecision(15)
         << " " << setw(22) << eup.PUP[i][0]
         << " " << setw(22) << eup.PUP[i][1]
         << " " << setw(22) << eup.PUP[i][2]
         << " " << setw(22) << eup.PUP[i][3]
         << " " << setw(22) << eup.PUP[i][4]
         << " " << setw(6) << eup.VTIMUP[i]
         << " " << setw(6) << eup.SPINUP[i] << endl;
  }

  // Write event comments.
  helper << hashline(eventStream.str()) << std::flush;
  eventStream.str("");

  if ( version != 1 ) {
    eup.rwgtSave.list(helper);
    eup.weightsSave.list(helper);
    eup.scalesSave.list(helper);
  }

  helper << "</event>" << endl;

  string helperString = helper.str();
  //event = helperString.c_str();

  return helperString;

}

//--------------------------------------------------------------------------

// Make sure that each line in the string s starts with a
// #-character and that the string ends with a new-line.

string Writer::hashline(string s, bool comment) {
  string ret;
  istringstream is(s);
  string ss;
  while ( getline(is, ss) ) {
    if ( comment )
        ss = "# " + ss;
    ret += ss + '\n';
  }
  return ret;
}

//==========================================================================

}
