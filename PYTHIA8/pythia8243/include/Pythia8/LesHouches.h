// LesHouches.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Les Houches Accord user process information.
// LHAProcess: stores a single process; used by the other classes.
// LHAParticle: stores a single particle; used by the other classes.
// LHAup: base class for initialization and event information.
// LHAupLHEF: derived class for reading from an Les Houches Event File.
// Code for interfacing with Fortran commonblocks is found in LHAFortran.h.

#ifndef Pythia8_LesHouches_H
#define Pythia8_LesHouches_H

#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/Streams.h"

namespace Pythia8 {

//==========================================================================

// A class for the processes stored in LHAup.

class LHAProcess {

public:

  // Constructors.
  LHAProcess() : idProc(0), xSecProc(0.), xErrProc(0.), xMaxProc(0.) { }
  LHAProcess(int idProcIn, double xSecIn, double xErrIn, double xMaxIn) :
    idProc(idProcIn), xSecProc(xSecIn), xErrProc(xErrIn),
    xMaxProc(xMaxIn) { }

  // Process properties.
  int    idProc;
  double xSecProc, xErrProc, xMaxProc;

} ;

//==========================================================================

// A class for the particles stored in LHAup.

class LHAParticle {

public:

  // Constructors.
  LHAParticle() : idPart(0), statusPart(0), mother1Part(0),
    mother2Part(0), col1Part(0), col2Part(0), pxPart(0.), pyPart(0.),
    pzPart(0.), ePart(0.), mPart(0.), tauPart(0.), spinPart(9.),
    scalePart(-1.) { }
  LHAParticle(int idIn, int statusIn, int mother1In, int mother2In,
    int col1In, int col2In, double pxIn, double pyIn, double pzIn,
    double eIn, double mIn, double tauIn, double spinIn,
    double scaleIn) :
    idPart(idIn), statusPart(statusIn), mother1Part(mother1In),
    mother2Part(mother2In), col1Part(col1In), col2Part(col2In),
    pxPart(pxIn), pyPart(pyIn), pzPart(pzIn), ePart(eIn), mPart(mIn),
    tauPart(tauIn), spinPart(spinIn), scalePart(scaleIn) { }

  // Particle properties.
  int    idPart, statusPart, mother1Part, mother2Part, col1Part, col2Part;
  double pxPart, pyPart, pzPart, ePart, mPart, tauPart, spinPart,
         scalePart;

} ;

//==========================================================================

// LHAup is base class for initialization and event information
// from an external parton-level generator.

class LHAup {

public:

  // Destructor.
  virtual ~LHAup() {}

  // Set info pointer.
  void setPtr(Info* infoPtrIn) {infoPtr = infoPtrIn;}

  // Method to be used for LHAupLHEF derived class.
  virtual void newEventFile(const char*) {}
  virtual bool fileFound() {return true;}
  virtual bool useExternal() {return false;}

  // A pure virtual method setInit, wherein all initialization information
  // is supposed to be set in the derived class. Can do this by reading a
  // file or some other way, as desired. Returns false if it did not work.
  virtual bool setInit() = 0;

  // Give back info on beams.
  int    idBeamA()       const {return idBeamASave;}
  int    idBeamB()       const {return idBeamBSave;}
  double eBeamA()        const {return eBeamASave;}
  double eBeamB()        const {return eBeamBSave;}
  int    pdfGroupBeamA() const {return pdfGroupBeamASave;}
  int    pdfGroupBeamB() const {return pdfGroupBeamBSave;}
  int    pdfSetBeamA()   const {return pdfSetBeamASave;}
  int    pdfSetBeamB()   const {return pdfSetBeamBSave;}

  // Give back weight strategy.
  int    strategy()      const {return strategySave;}

  // Give back info on processes.
  int    sizeProc()      const {return processes.size();}
  int    idProcess(int proc) const {return processes[proc].idProc;}
  double xSec(int proc)  const {return processes[proc].xSecProc;}
  double xErr(int proc)  const {return processes[proc].xErrProc;}
  double xMax(int proc)  const {return processes[proc].xMaxProc;}
  double xSecSum()       const {return xSecSumSave;}
  double xErrSum()       const {return xErrSumSave;}

  // Print the initialization info; useful to check that setting it worked.
  void   listInit();

  // A pure virtual method setEvent, wherein information on the next event
  // is supposed to be set in the derived class.
  // Strategies +-1 and +-2: idProcIn is the process type, selected by PYTHIA.
  // Strategies +-3 and +-4: idProcIn is dummy; process choice is made locally.
  // The method can find the next event by a runtime interface to another
  // program, or by reading a file, as desired.
  // The method should return false if it did not work.
  virtual bool setEvent(int idProcIn = 0) = 0;

  // Give back process number, weight, scale, alpha_em, alpha_s.
  int    idProcess()       const {return idProc;}
  double weight()          const {return weightProc;}
  double scale()           const {return scaleProc;}
  double alphaQED()        const {return alphaQEDProc;}
  double alphaQCD()        const {return alphaQCDProc;}

  // Give back info on separate particle.
  int    sizePart()        const {return particles.size();}
  int    id(int part)      const {return particles[part].idPart;}
  int    status(int part)  const {return particles[part].statusPart;}
  int    mother1(int part) const {return particles[part].mother1Part;}
  int    mother2(int part) const {return particles[part].mother2Part;}
  int    col1(int part)    const {return particles[part].col1Part;}
  int    col2(int part)    const {return particles[part].col2Part;}
  double px(int part)      const {return particles[part].pxPart;}
  double py(int part)      const {return particles[part].pyPart;}
  double pz(int part)      const {return particles[part].pzPart;}
  double e(int part)       const {return particles[part].ePart;}
  double m(int part)       const {return particles[part].mPart;}
  double tau(int part)     const {return particles[part].tauPart;}
  double spin(int part)    const {return particles[part].spinPart;}
  double scale(int part)   const {return particles[part].scalePart;}

  // Give back info on flavour and x values of hard-process initiators.
  int    id1()             const {return id1Save;}
  int    id2()             const {return id2Save;}
  double x1()              const {return x1Save;}
  double x2()              const {return x2Save;}

  // Optional: give back info on parton density values of event.
  bool   pdfIsSet()        const {return pdfIsSetSave;}
  int    id1pdf()          const {return id1pdfSave;}
  int    id2pdf()          const {return id2pdfSave;}
  double x1pdf()           const {return x1pdfSave;}
  double x2pdf()           const {return x2pdfSave;}
  double scalePDF()        const {return scalePDFSave;}
  double pdf1()            const {return pdf1Save;}
  double pdf2()            const {return pdf2Save;}

  // Optional: give back info on parton shower scales.
  bool   scaleShowersIsSet() const {return scaleShowersIsSetSave;}
  double scaleShowers(int i) const {return scaleShowersSave[i];}

  // Print the info; useful to check that reading an event worked.
  void   listEvent();

  // Skip ahead a number of events, which are not considered further.
  // Mainly intended for debug when using the LHAupLHEF class.
  virtual bool skipEvent(int nSkip) {
    for (int iSkip = 0; iSkip < nSkip; ++iSkip) if (!setEvent()) return false;
    return true;}

  // Four routines to write a Les Houches Event file in steps.
  virtual bool openLHEF(string fileNameIn);
  virtual bool closeLHEF(bool updateInit = false);
  bool   initLHEF();
  bool   eventLHEF(bool verbose = true);

  // Get access to the Les Houches Event file name.
  string getFileName()     const {return fileName;}

protected:

  // Constructor. Sets default to be that events come with unit weight.
  LHAup(int strategyIn = 3) : infoPtr(), nupSave(), idprupSave(), xwgtupSave(),
    scalupSave(), aqedupSave(), aqcdupSave(), xSecSumSave(), xErrSumSave(),
    getPDFSave(), getScale(), getScaleShowers(), id1InSave(), id2InSave(),
    id1pdfInSave(), id2pdfInSave(), x1InSave(), x2InSave(), x1pdfInSave(),
    x2pdfInSave(), scalePDFInSave(), pdf1InSave(), pdf2InSave(),
    scaleShowersInSave(), fileName("void"), dateNow(), timeNow(),
    strategySave(strategyIn), idBeamASave(), idBeamBSave(), eBeamASave(),
    eBeamBSave(), pdfGroupBeamASave(), pdfGroupBeamBSave(), pdfSetBeamASave(),
    pdfSetBeamBSave(), idProc(), weightProc(), scaleProc(), alphaQEDProc(),
    alphaQCDProc(), pdfIsSetSave(), scaleShowersIsSetSave(false), id1Save(),
    id2Save(), id1pdfSave(), id2pdfSave(), x1Save(), x2Save(), x1pdfSave(),
    x2pdfSave(), scalePDFSave(), pdf1Save(), pdf2Save(), scaleShowersSave() {
    processes.reserve(10); particles.reserve(20);
    setBeamA( 0, 0., 0, 0); setBeamB( 0, 0., 0, 0); }

  // Allow conversion from mb to pb.
  static const double CONVERTMB2PB;

  // Pointer to various information on the generation.
  Info* infoPtr;

  // Input beam info.
  void setBeamA(int idIn, double eIn, int pdfGroupIn = 0, int pdfSetIn = 0)
    { idBeamASave = idIn; eBeamASave = eIn; pdfGroupBeamASave = pdfGroupIn;
    pdfSetBeamASave = pdfSetIn;}
  void setBeamB(int idIn, double eIn, int pdfGroupIn = 0, int pdfSetIn = 0)
    { idBeamBSave = idIn; eBeamBSave = eIn; pdfGroupBeamBSave = pdfGroupIn;
    pdfSetBeamBSave = pdfSetIn;}

  // Input process weight strategy.
  void setStrategy(int strategyIn) {strategySave = strategyIn;}

  // Input process info.
  void addProcess(int idProcIn, double xSecIn = 1., double xErrIn = 0.,
    double xMaxIn = 1.) { processes.push_back( LHAProcess( idProcIn,
    xSecIn, xErrIn, xMaxIn)); }

  // Possibility to update some cross section info at end of run.
  void setXSec(int iP, double xSecIn) {processes[iP].xSecProc = xSecIn;}
  void setXErr(int iP, double xErrIn) {processes[iP].xErrProc = xErrIn;}
  void setXMax(int iP, double xMaxIn) {processes[iP].xMaxProc = xMaxIn;}

  // Input info on the selected process.
  void setProcess(int idProcIn = 0, double weightIn = 1., double
    scaleIn = 0., double alphaQEDIn = 0.0073, double alphaQCDIn = 0.12) {
    idProc = idProcIn; weightProc = weightIn; scaleProc = scaleIn;
    alphaQEDProc = alphaQEDIn; alphaQCDProc = alphaQCDIn;
    // Clear particle list. Add empty zeroth particle for correct indices.
    particles.clear(); addParticle(0); pdfIsSetSave = false;
    scaleShowersIsSetSave = false;}

  // Input particle info, one particle at the time.
  void addParticle(LHAParticle particleIn) {
    particles.push_back(particleIn);}
  void addParticle(int idIn, int statusIn = 0, int mother1In = 0,
    int mother2In = 0, int col1In = 0, int col2In = 0, double pxIn = 0.,
    double pyIn = 0., double pzIn = 0., double eIn = 0., double mIn = 0.,
    double tauIn = 0., double spinIn = 9., double scaleIn = -1.) {
    particles.push_back( LHAParticle( idIn, statusIn, mother1In, mother2In,
    col1In, col2In, pxIn, pyIn, pzIn, eIn, mIn, tauIn, spinIn,
    scaleIn) ); }

  // Input info on flavour and x values of hard-process initiators.
  void setIdX(int id1In, int id2In, double x1In, double x2In)
    { id1Save = id1In; id2Save = id2In; x1Save = x1In; x2Save = x2In;}

  // Optionally input info on parton density values of event.
  void setPdf(int id1pdfIn, int id2pdfIn, double x1pdfIn, double x2pdfIn,
    double scalePDFIn, double pdf1In, double pdf2In, bool pdfIsSetIn)
    { id1pdfSave = id1pdfIn; id2pdfSave = id2pdfIn; x1pdfSave = x1pdfIn;
    x2pdfSave = x2pdfIn; scalePDFSave = scalePDFIn; pdf1Save = pdf1In;
    pdf2Save = pdf2In; pdfIsSetSave = pdfIsSetIn;}

  // Optionally input info on parton shower starting scale; two for DPS.
  void setScaleShowers( double scaleIn1, double scaleIn2 = 0.)
    { scaleShowersSave[0] = scaleIn1; scaleShowersSave[1] = scaleIn2;
    scaleShowersIsSetSave = true;}

  // Three routines for LHEF files, but put here for flexibility.
  bool setInitLHEF(istream& is, bool readHeaders = false);
  bool setNewEventLHEF(istream& is);
  bool setOldEventLHEF();

  // Helper routines to open and close a file handling GZIPSUPPORT:
  //   ifstream ifs;
  //   istream *is = openFile("myFile.txt", ifs);
  //   -- Process file using is --
  //   closeFile(is, ifs);
  istream* openFile(const char *fn, ifstream &ifs);
  void     closeFile(istream *&is, ifstream &ifs);

  // LHAup is a friend class to infoPtr, but derived classes
  // are not. This wrapper function can be used by derived classes
  // to set headers in the Info class.
  void setInfoHeader(const string &key, const string &val) {
    infoPtr->setHeader(key, val); }

  // Event properties from LHEF files, for repeated use.
  int    nupSave, idprupSave;
  double xwgtupSave, scalupSave, aqedupSave, aqcdupSave, xSecSumSave,
         xErrSumSave;
  vector<LHAParticle> particlesSave;
  bool   getPDFSave, getScale, getScaleShowers;
  int    id1InSave, id2InSave, id1pdfInSave, id2pdfInSave;
  double x1InSave, x2InSave, x1pdfInSave, x2pdfInSave, scalePDFInSave,
         pdf1InSave, pdf2InSave, scaleShowersInSave[2];

  // File to which to write Les Houches Event File information.
  string fileName;
  fstream osLHEF;
  char dateNow[12];
  char timeNow[9];

private:

  // Event weighting and mixing strategy.
  int strategySave;

  // Beam particle properties.
  int    idBeamASave, idBeamBSave;
  double eBeamASave, eBeamBSave;
  int    pdfGroupBeamASave, pdfGroupBeamBSave, pdfSetBeamASave,
         pdfSetBeamBSave;

  // The process list, stored as a vector of processes.
  vector<LHAProcess> processes;

  // Store info on the selected process.
  int    idProc;
  double weightProc, scaleProc, alphaQEDProc, alphaQCDProc;

  // The particle list, stored as a vector of particles.
  vector<LHAParticle> particles;

  // Info on initiators and optionally on parton density values of event.
  bool   pdfIsSetSave, scaleShowersIsSetSave;
  int    id1Save, id2Save, id1pdfSave, id2pdfSave;
  double x1Save, x2Save, x1pdfSave, x2pdfSave, scalePDFSave, pdf1Save,
         pdf2Save, scaleShowersSave[2];

};

//==========================================================================

// A derived class with information read from a Les Houches Event File.

class LHAupLHEF : public LHAup {

public:

  // Constructor.
  LHAupLHEF(Pythia8::Info* infoPtrIn, istream* isIn, istream* isHeadIn,
    bool readHeadersIn = false, bool setScalesFromLHEFIn = false ) :
    infoPtr(infoPtrIn), filename(""), headerfile(""),
    is(isIn), is_gz(NULL), isHead(isHeadIn), isHead_gz(NULL),
    readHeaders(readHeadersIn), reader(is),
    setScalesFromLHEF(setScalesFromLHEFIn), hasExtFileStream(true),
    hasExtHeaderStream(true) {}

  LHAupLHEF(Pythia8::Info* infoPtrIn, const char* filenameIn,
    const char* headerIn = NULL, bool readHeadersIn = false,
    bool setScalesFromLHEFIn = false ) :
    infoPtr(infoPtrIn), filename(filenameIn), headerfile(headerIn),
    is(NULL), is_gz(NULL), isHead(NULL), isHead_gz(NULL),
    readHeaders(readHeadersIn), reader(filenameIn),
    setScalesFromLHEF(setScalesFromLHEFIn), hasExtFileStream(false),
    hasExtHeaderStream(false) {
    is = (openFile(filenameIn, ifs));
    isHead = (headerfile == NULL) ? is : openFile(headerfile, ifsHead);
    is_gz = new igzstream(filename);
    isHead_gz = (headerfile == NULL) ? is_gz : new igzstream(headerfile);
  }

  // Destructor.
  ~LHAupLHEF() {
     // Close files
     closeAllFiles();
  }

  // Helper routine to correctly close files.
  void closeAllFiles() {

    if (!hasExtHeaderStream && isHead_gz != is_gz) isHead_gz->close();
    if (isHead_gz != is_gz) delete isHead_gz;
    if (is_gz) is_gz->close();
    if (is_gz) delete is_gz;

    // Close header file if separate, and close main file.
    if (!hasExtHeaderStream && isHead != is) closeFile(isHead, ifsHead);
    if (!hasExtFileStream) closeFile(is, ifs);
  }

  // Want to use new file with events, but without reinitialization.
  void newEventFile(const char* filenameIn) {
    // Close files and then open new file.
    closeAllFiles();
    is    = (openFile(filenameIn, ifs));
    is_gz = new igzstream(filenameIn);
    // Re-initialise Les Houches file reader.
    reader.setup(filenameIn);
    // Set isHead to is to keep expected behaviour in
    // fileFound() and closeAllFiles().
    isHead    = is;
    isHead_gz = is_gz;
  }

  // Confirm that file was found and opened as expected.
  bool fileFound() {return (useExternal() || (isHead->good() && is->good()));}
  bool useExternal() {return (hasExtHeaderStream && hasExtFileStream);}

  // Routine for doing the job of reading and setting initialization info.
  bool setInit() {
    return setInitLHEF(*isHead, readHeaders);
  }

  // Routine for doing the job of reading and setting initialization info.
  bool setInitLHEF( istream & isIn, bool readHead);

  // Routine for doing the job of reading and setting info on next event.
  bool setEvent(int = 0) {
    if (!setNewEventLHEF()) return false;
    return setOldEventLHEF();
  }

  // Skip ahead a number of events, which are not considered further.
  bool skipEvent(int nSkip) {
    for (int iSkip = 0; iSkip < nSkip; ++iSkip)
      if (!setNewEventLHEF()) return false;
     return true;
  }

  // Routine for doing the job of reading and setting info on next event.
  bool setNewEventLHEF();

  // Update cross-section information at the end of the run.
  bool updateSigma() {return true;}

protected:

  // Used internally to read a single line from the stream.
  bool getLine(string & line, bool header = true) {
#ifdef GZIPSUPPORT
    if      ( isHead_gz &&  header && !getline(*isHead_gz, line)) return false;
    else if ( is_gz     && !header && !getline(*is_gz, line))     return false;
    if      (header && !getline(*isHead, line)) return false;
    else if (!header && !getline(*is, line))    return false;
#else
    if      (header && !getline(*isHead, line)) return false;
    else if (!header && !getline(*is, line))    return false;
#endif
    // Replace single by double quotes
    replace(line.begin(),line.end(),'\'','\"');
    return true;
  }

private:

  Info*  infoPtr;
  const char* filename;
  const char* headerfile;

  // File from which to read (or a stringstream).
  // Optionally also a file from which to read the LHEF header.
  istream  *is;
  igzstream  *is_gz;
  ifstream  ifs;
  istream  *isHead;
  igzstream  *isHead_gz;
  ifstream  ifsHead;

  // Flag to read headers or not
  bool readHeaders;

  Reader reader;

  // Flag to set particle production scales or not.
  bool setScalesFromLHEF, hasExtFileStream, hasExtHeaderStream;

};

//==========================================================================

// A derived class with information read from PYTHIA 8 itself, for output.

class LHAupFromPYTHIA8 : public LHAup {

public:

  // Constructor.
  LHAupFromPYTHIA8(Event* processPtrIn, Info* infoPtrIn) {
    processPtr = processPtrIn; infoPtr = infoPtrIn;}

  // Destructor.
  ~LHAupFromPYTHIA8() {}

  // Routine for doing the job of reading and setting initialization info.
  bool setInit();

  // Routine for doing the job of reading and setting info on next event.
  bool setEvent(int = 0);

  // Update cross-section information at the end of the run.
  bool updateSigma();

private:

  // Pointers to process event record and further information.
  Event* processPtr;
  Info*  infoPtr;

};

//==========================================================================

// A derived class with LHEF 3.0 information read from PYTHIA 8 itself, for
// output.

class LHEF3FromPythia8 : public LHAup {

public:

  // Constructor.
  LHEF3FromPythia8(Event* eventPtrIn, Settings* settingsPtrIn,
    Info* infoPtrIn, ParticleData* particleDataPtrIn, int pDigitsIn = 15,
    bool writeToFileIn = true) :
    eventPtr(eventPtrIn),settingsPtr(settingsPtrIn), infoPtr(infoPtrIn),
    particleDataPtr(particleDataPtrIn), writer(osLHEF), pDigits(pDigitsIn),
    writeToFile(writeToFileIn) {}

  // Routine for reading, setting and printing the initialisation info.
  bool setInit();

  // Routine for reading, setting and printing the next event.
  void setEventPtr(Event* evPtr) { eventPtr = evPtr; }
  bool setEvent(int = 0);
  string getEventString() { return writer.getEventString(&hepeup); }

  // Function to open the output file.
  bool openLHEF(string fileNameIn);

  // Function to close (and possibly update) the output file.
  bool closeLHEF(bool updateInit = false);

private:

  // Pointer to event that should be printed.
  Event* eventPtr;

  // Pointer to settings and info objects.
  Settings* settingsPtr;
  Info* infoPtr;
  ParticleData* particleDataPtr;

  // LHEF3 writer
  Writer writer;

  // Number of digits to set width of double write out
  int  pDigits;
  bool writeToFile;

  // Some internal init and event block objects for convenience.
  HEPRUP heprup;
  HEPEUP hepeup;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_LesHouches_H
