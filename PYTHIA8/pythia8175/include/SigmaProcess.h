// SigmaProcess.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for hard-process differential cross sections.
// SigmaProcess: base class for cross sections.
// Sigma0Process: base class for unresolved processes, derived from above.
// Sigma1Process: base class for 2 -> 1 processes, derived from above.
// Sigma2Process: base class for 2 -> 2 processes, derived from above.
// Sigma3Process: base class for 2 -> 3 processes, derived from above.
// SigmaLHAProcess: wrapper class for Les Houches Accord external input.
// Actual physics processes are found in separate files:
// SigmaQCD for QCD processes;
// SigmaEW for electroweak processes (including photon production);
// SigmaOnia for charmonium and bottomonium processes;
// SigmaHiggs for Higgs processes;
// SigmaSUSY for supersymmetric production;
// SigmaLeftRightSym for  processes in left-right-symmetric scenarios;
// SigmaLeptoquark for leptoquark production.
// SigmaExtraDim for processes in scenarios for extra dimensions;
// SigmaGeneric for generic scalar/fermion/vector pair production.

#ifndef Pythia8_SigmaProcess_H
#define Pythia8_SigmaProcess_H

#include "Basics.h"
#include "BeamParticle.h"
#include "Event.h"
#include "Info.h"
#include "LesHouches.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "PythiaComplex.h"
#include "PythiaStdlib.h"
#include "ResonanceWidths.h"
#include "Settings.h"
#include "SigmaTotal.h"
#include "StandardModel.h"
#include "SusyLesHouches.h"

namespace Pythia8 {

//==========================================================================

// InBeam is a simple helper class for partons and their flux in a beam.

class InBeam {

public:

  // Constructor.
  InBeam( int idIn = 0) : id(idIn), pdf(0.) {}

  // Values.
  int    id;
  double pdf;

};

//==========================================================================

// InPair is a simple helper class for colliding parton pairs and their flux.

class InPair {

public:

  // Constructor.
  InPair( int idAIn = 0, int idBIn = 0) : idA(idAIn), idB(idBIn),
    pdfA(0.), pdfB(0.), pdfSigma(0.) {}

  // Values.
  int    idA, idB;
  double pdfA, pdfB, pdfSigma;

};
 
//==========================================================================

// SigmaProcess is the base class for cross section calculations.

class SigmaProcess {

public:

  // Destructor.
  virtual ~SigmaProcess() {}

  // Perform simple initialization and store pointers.
  void init(Info* infoPtrIn, Settings* settingsPtrIn,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn, Couplings* couplings, 
	    SigmaTotal* sigmaTotPtrIn = 0, SusyLesHouches* slhaPtr = 0); 

  // Store or replace Les Houches pointer.
  void setLHAPtr( LHAup* lhaUpPtrIn) {lhaUpPtr = lhaUpPtrIn;}  

  // Initialize process. Only used for some processes.
  virtual void initProc() {}

  // Set up allowed flux of incoming partons. Default is no flux.
  virtual bool initFlux();

  // Input and complement kinematics for resolved 2 -> 1 process.
  // Usage: set1Kin( x1in, x2in, sHin).
  virtual void set1Kin( double , double , double ) {} 

  // Input and complement kinematics for resolved 2 -> 2 process.
  // Usage: set2Kin( x1in, x2in, sHin, tHin, m3in, m4in, runBW3in, runBW4in).
  virtual void set2Kin( double , double , double , double , double , 
    double, double, double ) {} 

  // Ditto, but for Multiparton Interactions applications, so different input.
  // Usage: set2KinMPI( x1in, x2in, sHin, tHin, uHin, 
  //                   alpSin, alpEMin, needMasses, m3in, m4in)
  virtual void set2KinMPI( double , double , double , double , 
    double , double , double , bool , double , double ) {}

  // Input and complement kinematics for resolved 2 -> 3 process.
  // Usage: set3Kin( x1in, x2in, sHin, p3prel, p4prel, p5prel, 
  //                 m3in, m4in, m5in, runBW3in, runBW4in, runBW5in); 
  virtual void set3Kin( double , double , double , Vec4 , Vec4 , Vec4 , 
    double , double , double , double , double , double ) {}

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin() {}

  // Evaluate sigma for unresolved, sigmaHat(sHat) for 2 -> 1 processes, 
  // d(sigmaHat)/d(tHat) for (resolved) 2 -> 2 processes, and |M|^2 for 
  // 2 -> 3 processes. Answer in "native" units, either mb or GeV^-2. 
  virtual double sigmaHat() {return 0.;}

  // Wrapper to sigmaHat, to (a) store current incoming flavours and 
  // (b) convert from GeV^-2 to mb where required.
  // For 2 -> 1/2 also (c) convert from from |M|^2 to d(sigmaHat)/d(tHat).
  virtual double sigmaHatWrap(int id1in = 0, int id2in = 0) {
    id1 = id1in; id2 = id2in; 
    return ( convert2mb() ? CONVERT2MB * sigmaHat() : sigmaHat() ); }

  // Convolute above with parton flux and K factor. Sum over open channels. 
  virtual double sigmaPDF();

  // Select incoming parton channel and extract parton densities (resolved).
  void pickInState(int id1in = 0, int id2in = 0);

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol() {}

  // Perform kinematics for a Multiparton Interaction, in its rest frame.
  virtual bool final2KinMPI( int = 0, int = 0, Vec4 = 0., Vec4 = 0.,
    double = 0., double = 0.) {return true;}

  // Evaluate weight for simultaneous flavours (only gamma*/Z0 gamma*/Z0).
  // Usage: weightDecayFlav( process).
  virtual double weightDecayFlav( Event&) {return 1.;} 

  // Evaluate weight for decay angular configuration.
  // Usage: weightDecay( process, iResBeg, iResEnd), where 
  // iResBeg <= i < iResEnd is range of sister partons to test decays of.
  virtual double weightDecay( Event&, int, int) {return 1.;}

  // Set scale, when that is missing for an external LHA process.
  virtual void setScale() {} 

  // Process name and code, and the number of final-state particles.
  virtual string name()            const {return "unnamed process";}
  virtual int    code()            const {return 0;}
  virtual int    nFinal()          const {return 2;}

  // Need to know which incoming partons to set up interaction for.
  virtual string inFlux()          const {return "unknown";}

  // Need to know whether to convert cross section answer from GeV^-2 to mb.
  virtual bool   convert2mb()      const {return true;}

  // For 2 -> 2 process optional conversion from |M|^2 to d(sigmaHat)/d(tHat).
  virtual bool   convertM2()       const {return false;}

  // Special treatment needed for Les Houches processes.
  virtual bool   isLHA()           const {return false;}

  // Special treatment needed for elastic and diffractive processes.
  virtual bool   isMinBias()       const {return false;}
  virtual bool   isResolved()      const {return true;}
  virtual bool   isDiffA()         const {return false;}
  virtual bool   isDiffB()         const {return false;}
  virtual bool   isDiffC()         const {return false;}

  // Special treatment needed for SUSY processes.
  virtual bool   isSUSY()          const {return false;}  

  // Special treatment needed if negative cross sections allowed.
  virtual bool   allowNegativeSigma() const {return false;}

  // Flavours in 2 -> 2/3 processes where masses needed from beginning. 
  // (For a light quark masses will be used in the final kinematics,
  // but not at the matrix-element level. For a gluon no masses at all.) 
  virtual int    id3Mass()         const {return 0;}
  virtual int    id4Mass()         const {return 0;}
  virtual int    id5Mass()         const {return 0;}

  // Special treatment needed if process contains an s-channel resonance.
  virtual int    resonanceA()      const {return 0;}
  virtual int    resonanceB()      const {return 0;}

  // 2 -> 2 and 2 -> 3 processes only through s-channel exchange.
  virtual bool   isSChannel()      const {return false;}

  // NOAM: Insert an intermediate resonance in 2 -> 1 -> 2 (or 3) listings.
  virtual int    idSChannel()      const {return 0;}

  // QCD 2 -> 3 processes need special phase space selection machinery.
  virtual bool   isQCD3body()      const {return false;}

  // Special treatment in 2 -> 3 with two massive propagators.
  virtual int    idTchan1()        const {return 0;}
  virtual int    idTchan2()        const {return 0;}
  virtual double tChanFracPow1()   const {return 0.3;}
  virtual double tChanFracPow2()   const {return 0.3;}
  virtual bool   useMirrorWeight() const {return false;}

  // Special process-specific gamma*/Z0 choice if >=0 (e.g. f fbar -> H0 Z0).
  virtual int    gmZmode()         const {return -1;}

  // Tell whether tHat and uHat are swapped (= same as swap 3 and 4).
  bool swappedTU()          const {return swapTU;}
  
  // Give back particle properties: flavours, colours, masses, or all.
  int    id(int i)          const {return idSave[i];}
  int    col(int i)         const {return colSave[i];} 
  int    acol(int i)        const {return acolSave[i];}
  double m(int i)           const {return mSave[i];}
  Particle getParton(int i) const {return parton[i];}

  // Give back couplings and parton densities. Not all known for minbias.
  double Q2Ren()            const {return Q2RenSave;}
  double alphaEMRen()       const {return alpEM;}
  double alphaSRen()        const {return alpS;}
  double Q2Fac()            const {return Q2FacSave;}
  double pdf1()             const {return pdf1Save;}
  double pdf2()             const {return pdf2Save;}

  // Give back angles; relevant only for multipe-interactions processes.
  double thetaMPI()         const {return atan2( sinTheta, cosTheta);}
  double phiMPI()           const {return phi;}
  double sHBetaMPI()        const {return sHBeta;}
  double pT2MPI()           const {return pT2Mass;}
  double pTMPIFin()         const {return pTFin;}

  // Save and load kinematics for trial interactions
  void saveKin() {
    for (int i = 0; i < 6; i++) { partonT[i] = parton[i]; 
      mSaveT[i] = mSave[i]; }
    pTFinT = pTFin; phiT = phi; cosThetaT = cosTheta; sinThetaT = sinTheta; }
  void loadKin() {
    for (int i = 0; i < 6; i++) { parton[i] = partonT[i]; 
    mSave[i] = mSaveT[i]; }
    pTFin = pTFinT; cosTheta = cosThetaT; sinTheta = sinThetaT; phi = phiT;
  }
  void swapKin() {
    for (int i = 0; i < 6; i++) { swap(parton[i], partonT[i]);
                                  swap(mSave[i], mSaveT[i]); }
    swap(pTFin, pTFinT); swap(cosTheta, cosThetaT);
    swap(sinTheta, sinThetaT); swap(phi, phiT); }

protected:

  // Constructor.
  SigmaProcess() : infoPtr(0), settingsPtr(0), particleDataPtr(0),
    rndmPtr(0), beamAPtr(0), beamBPtr(0), couplingsPtr(0), sigmaTotPtr(0),
    slhaPtr(0), lhaUpPtr(0) {for (int i = 0; i < 6; ++i) mSave[i] = 0.;}

  // Constants: could only be changed in the code itself.
  static const double CONVERT2MB, MASSMARGIN, COMPRELERR;
  static const int    NCOMPSTEP;

  // Pointer to various information on the generation.
  Info*           infoPtr;
 
  // Pointer to the settings database.
  Settings*       settingsPtr;

  // Pointer to the particle data table.
  ParticleData*   particleDataPtr;

  // Pointer to the random number generator.
  Rndm*           rndmPtr;

  // Pointers to incoming beams.
  BeamParticle*   beamAPtr;
  BeamParticle*   beamBPtr;

  // Pointer to Standard Model couplings, including alphaS and alphaEM.
  Couplings*      couplingsPtr;
  
  // Pointer to the total/elastic/diffractive cross section object.
  SigmaTotal*     sigmaTotPtr;

  // Pointer to the SLHA object.
  SusyLesHouches* slhaPtr;

  // Pointer to LHAup for generating external events.
  LHAup*          lhaUpPtr;

  // Initialization data, normally only set once.
  int    nQuarkIn, renormScale1, renormScale2, renormScale3, renormScale3VV, 
         factorScale1, factorScale2, factorScale3, factorScale3VV;
  double Kfactor, mcME, mbME, mmuME, mtauME, renormMultFac, renormFixScale, 
         factorMultFac, factorFixScale;

  // CP violation parameters for Higgs sector, normally only set once.
  int    higgsH1parity, higgsH2parity, higgsA3parity;
  double higgsH1eta, higgsH2eta, higgsA3eta;  

  // Information on incoming beams.
  int    idA, idB;
  double mA, mB; 
  bool   isLeptonA, isLeptonB, hasLeptonBeams;

  // Partons in beams, with PDF's.
  vector<InBeam> inBeamA;
  vector<InBeam> inBeamB;
  void addBeamA(int idIn) {inBeamA.push_back(InBeam(idIn));}
  void addBeamB(int idIn) {inBeamB.push_back(InBeam(idIn));}
  int sizeBeamA() const {return inBeamA.size();}
  int sizeBeamB() const {return inBeamB.size();}
 
  // Allowed colliding parton pairs, with pdf's.
  vector<InPair> inPair;
  void addPair(int idAIn, int idBIn) {
    inPair.push_back(InPair(idAIn, idBIn));}
  int sizePair() const {return inPair.size();}

  // Store common subprocess kinematics quantities.
  double mH, sH, sH2;

  // Store Q2 renormalization and factorization scales, and related values.
  double Q2RenSave, alpEM, alpS, Q2FacSave, x1Save, x2Save, pdf1Save, 
         pdf2Save, sigmaSumSave;

  // Store flavour, colour, anticolour, mass, angles and the whole particle.
  int      id1, id2, id3, id4, id5;
  int      idSave[6], colSave[6], acolSave[6];
  double   mSave[6], cosTheta, sinTheta, phi, sHMass, sHBeta, pT2Mass, pTFin;
  Particle parton[6];

  // Minimal set of saved kinematics for trial interactions when
  // using the x-dependent matter profile of multiparton interactions.
  Particle partonT[6];
  double   mSaveT[6], pTFinT, cosThetaT, sinThetaT, phiT;

  // Calculate and store all modified masses and four-vectors 
  // intended for matrix elements. Return false if failed.
  virtual bool setupForME() {return true;}
  bool     setupForMEin();
  double   mME[5];
  Vec4     pME[5];

  // Store whether tHat and uHat are swapped (= same as swap 3 and 4).
  bool swapTU;

  // Set flavour, colour and anticolour.
  void setId( int id1in = 0, int id2in = 0, int id3in = 0, int id4in = 0,
    int id5in = 0) {idSave[1] = id1in; idSave[2] = id2in; idSave[3] = id3in; 
    idSave[4] = id4in; idSave[5] = id5in;}
  void setColAcol( int col1 = 0, int acol1 = 0, 
    int col2 = 0, int acol2 = 0, int col3 = 0, int acol3 = 0, 
    int col4 = 0, int acol4 = 0, int col5 = 0, int acol5 = 0) {
    colSave[1] = col1; acolSave[1] = acol1; colSave[2] = col2; 
    acolSave[2] = acol2; colSave[3] = col3; acolSave[3] = acol3; 
    colSave[4] = col4; acolSave[4] = acol4; colSave[5] = col5; 
    acolSave[5] = acol5; }
  void swapColAcol() { swap(colSave[1], acolSave[1]); 
    swap(colSave[2], acolSave[2]); swap(colSave[3], acolSave[3]); 
    swap(colSave[4], acolSave[4]); swap(colSave[5], acolSave[5]);}
  void swapCol1234() { swap(colSave[1], colSave[2]); 
    swap(colSave[3], colSave[4]); swap(acolSave[1], acolSave[2]); 
    swap(acolSave[3], acolSave[4]);}
  void swapCol12() { swap(colSave[1], colSave[2]); 
    swap(acolSave[1], acolSave[2]);}
  void swapCol34() { swap(colSave[3], colSave[4]); 
    swap(acolSave[3], acolSave[4]);}

  // Common code for top and Higgs secondary decay angular weights.
  double weightTopDecay( Event& process, int iResBeg, int iResEnd);
  double weightHiggsDecay( Event& process, int iResBeg, int iResEnd);

};
 
//==========================================================================

// Sigma0Process is the base class for unresolved and minimum-bias processes. 
// It is derived from SigmaProcess.

class Sigma0Process : public SigmaProcess {

public:

  // Destructor.
  virtual ~Sigma0Process() {}

  // Number of final-state particles.
  virtual int    nFinal() const {return 2;};

  // No partonic flux to be set up.
  virtual bool   initFlux() {return true;}

  // Evaluate sigma for unresolved processes. 
  virtual double sigmaHat() {return 0.;}

  // Since no PDF's there is no difference from above. 
  virtual double sigmaPDF() {return sigmaHat();}

  // Answer for these processes already in mb, so do not convert.
  virtual bool convert2mb() const {return false;}

protected:

  // Constructor.
  Sigma0Process() {}

};
 
//==========================================================================

// Sigma1Process is the base class for 2 -> 1 processes.
// It is derived from SigmaProcess.

class Sigma1Process : public SigmaProcess {

public:

  // Destructor.
  virtual ~Sigma1Process() {}

  // Number of final-state particles.
  virtual int    nFinal() const {return 1;};

  // Input and complement kinematics for resolved 2 -> 1 process.
  virtual void   set1Kin( double x1in, double x2in, double sHin) {
    store1Kin( x1in, x2in, sHin); sigmaKin();} 

  // Evaluate sigmaHat(sHat) for resolved 2 -> 1 processes. 
  virtual double sigmaHat() {return 0.;}

  // Wrapper to sigmaHat, to (a) store current incoming flavours, 
  // (b) convert from GeV^-2 to mb where required, and
  // (c) convert from |M|^2 to d(sigmaHat)/d(tHat) where required.
  virtual double sigmaHatWrap(int id1in = 0, int id2in = 0); 

protected:

  // Constructor.
  Sigma1Process() {}

  // Store kinematics and set scales for resolved 2 -> 1 process.
  virtual void   store1Kin( double x1in, double x2in, double sHin);

  // Calculate modified masses and four-vectors for matrix elements.
  virtual bool   setupForME();

};
 
//==========================================================================

// Sigma2Process is the base class for 2 -> 2 processes.
// It is derived from SigmaProcess.

class Sigma2Process : public SigmaProcess {

public:

  // Destructor.
  virtual ~Sigma2Process() {}

  // Number of final-state particles.
  virtual int    nFinal() const {return 2;};

  // Input and complement kinematics for resolved 2 -> 2 process.
  virtual void   set2Kin( double x1in, double x2in, double sHin, 
    double tHin, double m3in, double m4in, double runBW3in, 
    double runBW4in) { store2Kin( x1in, x2in, sHin, tHin, m3in, m4in, 
    runBW3in, runBW4in); sigmaKin();}

  // Ditto, but for Multiparton Interactions applications, so different input.
  virtual void   set2KinMPI( double x1in, double x2in, double sHin, 
    double tHin, double uHin, double alpSin, double alpEMin, 
    bool needMasses, double m3in, double m4in) {
    store2KinMPI( x1in, x2in, sHin, tHin, uHin, alpSin, alpEMin, 
    needMasses, m3in, m4in); sigmaKin();}

  // Evaluate d(sigmaHat)/d(tHat) for resolved 2 -> 2 processes. 
  virtual double sigmaHat() {return 0.;}

  // Wrapper to sigmaHat, to (a) store current incoming flavours, 
  // (b) convert from GeV^-2 to mb where required, and
  // (c) convert from |M|^2 to d(sigmaHat)/d(tHat) where required.
  virtual double sigmaHatWrap(int id1in = 0, int id2in = 0) {
    id1 = id1in; id2 = id2in; double sigmaTmp = sigmaHat(); 
    if (convertM2())  sigmaTmp /= 16. * M_PI * sH2; 
    if (convert2mb()) sigmaTmp *= CONVERT2MB; return sigmaTmp;}

  // Perform kinematics for a Multiparton Interaction, in its rest frame.
  virtual bool   final2KinMPI( int i1Res = 0, int i2Res = 0, Vec4 p1Res = 0., 
    Vec4 p2Res = 0., double m1Res = 0., double m2Res = 0.);

protected:

  // Constructor.
  Sigma2Process() : tH(0.), uH(0.), tH2(0.), uH2(0.), m3(0.), s3(0.), 
    m4(0.), s4(0.), pT2(0.), runBW3(0.), runBW4(0.) {}

  // Store kinematics and set scales for resolved 2 -> 2 process.
  virtual void   store2Kin( double x1in, double x2in, double sHin, 
    double tHin, double m3in, double m4in, double runBW3in, 
    double runBW4in);
  virtual void   store2KinMPI( double x1in, double x2in, double sHin, 
    double tHin, double uHin, double alpSin, double alpEMin, 
    bool needMasses, double m3in, double m4in);

  // Calculate modified masses and four-vectors for matrix elements.
  virtual bool   setupForME();

  // Store subprocess kinematics quantities.
  double tH, uH, tH2, uH2, m3, s3, m4, s4, pT2, runBW3, runBW4;

};
 
//==========================================================================

// Sigma3Process is the base class for 2 -> 3 processes.
// It is derived from SigmaProcess.

class Sigma3Process : public SigmaProcess {

public:

  // Destructor.
  virtual ~Sigma3Process() {}

  // Number of final-state particles.
  virtual int    nFinal() const {return 3;};

  // Input and complement kinematics for resolved 2 -> 3 process.
  virtual void   set3Kin( double x1in, double x2in, double sHin, 
    Vec4 p3cmIn, Vec4 p4cmIn, Vec4 p5cmIn, double m3in, double m4in, 
    double m5in, double runBW3in, double runBW4in, double runBW5in) { 
    store3Kin( x1in, x2in, sHin, p3cmIn, p4cmIn, p5cmIn, m3in, m4in, m5in,
    runBW3in, runBW4in, runBW5in); sigmaKin();}

  // Evaluate d(sigmaHat)/d(tHat) for resolved 2 -> 3 processes. 
  virtual double sigmaHat() {return 0.;}

protected:

  // Constructor.
  Sigma3Process() {}

  // Store kinematics and set scales for resolved 2 -> 3 process.
  virtual void   store3Kin( double x1in, double x2in, double sHin, 
    Vec4 p3cmIn, Vec4 p4cmIn, Vec4 p5cmIn, double m3in, double m4in, 
    double m5in, double runBW3in, double runBW4in, double runBW5in);

  // Calculate modified masses and four-vectors for matrix elements.
  virtual bool   setupForME();

  // Store subprocess kinematics quantities.
  double m3, s3, m4, s4, m5, s5, runBW3, runBW4, runBW5;
  Vec4   p3cm, p4cm, p5cm;

};
 
//==========================================================================

// SigmaLHAProcess is a wrapper class for Les Houches Accord external input.
// It is derived from SigmaProcess.

class SigmaLHAProcess : public SigmaProcess {

public:

  // Constructor.
  SigmaLHAProcess() {}

  // Destructor.
  virtual ~SigmaLHAProcess() {}

  // No partonic flux to be set up.
  virtual bool   initFlux() {return true;}

  // Dummy function: action is put in PhaseSpaceLHA.
  virtual double sigmaPDF() {return 1.;}

  // Evaluate weight for decay angular configuration, where relevant.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd);

  // Set scale, when that is missing for an external LHA process.
  virtual void   setScale(); 

  // Info on the subprocess.
  virtual string name()     const {return "Les Houches User Process(es)";}
  virtual int    code()     const {return 9999;}

  // Number of final-state particles depends on current process choice.
  virtual int    nFinal()   const;
 
  // Answer for these processes not in GeV^-2, so do not do this conversion.
  virtual bool   convert2mb() const {return false;}

  // Ensure special treatment of Les Houches processes.
  virtual bool   isLHA()    const {return true;}

  // Special treatment needed if negative cross sections allowed.
  virtual bool   allowNegativeSigma() const {
    return (lhaUpPtr->strategy() < 0);}

private:

};
 
//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaProcess_H
 
