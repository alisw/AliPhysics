// ColourReconnection.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the Colour reconnection handling.
// Reconnect the colours between the partons before hadronization.
// It Contains the following classes:
// ColourDipole, ColourParticle, ColourJunction, ColourReconnection.

#ifndef Pythia8_ColourReconnection_H
#define Pythia8_ColourReconnection_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/StringFragmentation.h"
#include "Pythia8/PartonDistributions.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StringLength.h"

namespace Pythia8 {

//==========================================================================

// Contain a single colour chain. It always start from a quark and goes to
// an anti quark or from an anti-junction to at junction
// (or possible combinations).

class ColourDipole {

public:

  // Constructor.
  ColourDipole( int colIn = 0, int iColIn = 0, int iAcolIn = 0,
    int colReconnectionIn = 0, bool isJunIn = false, bool isAntiJunIn = false,
    bool isActiveIn = true, bool isRealIn = false) : col(colIn), iCol(iColIn),
    iAcol(iAcolIn), colReconnection(colReconnectionIn), isJun(isJunIn),
    isAntiJun(isAntiJunIn),isActive(isActiveIn), isReal(isRealIn)
    {leftDip = 0; rightDip = 0; iColLeg = 0; iAcolLeg = 0; printed = false;}

  double mDip(Event & event) {
    if (isJun || isAntiJun) return 1E9;
    else return m(event[iCol].p(),event[iAcol].p());
  }

  // Members.
  int    col, iCol, iAcol, iColLeg, iAcolLeg, colReconnection;
  bool   isJun, isAntiJun, isActive, isReal, printed, inChain;
  ColourDipole *leftDip, *rightDip;
  vector<ColourDipole *> colDips, acolDips;
  double p1p2;

  // Printing function, mainly intended for debugging.
  void print();

};

//==========================================================================

// Junction class. In addition to the normal junction class, also contains a
// list of dipoles connected to it.

class ColourJunction : public Junction {

public:

  ColourJunction(const Junction& ju) : Junction(ju) {
      for(int i = 0;i < 3;++i) {
        dips[i] = 0; dipsOrig[i] = 0;}
  }
  ColourJunction(const ColourJunction& ju) : Junction(Junction(ju)) {
    for(int i = 0;i < 3;++i) {
      dips[i] = ju.dips[i]; dipsOrig[i] = ju.dipsOrig[i];}
  }
  ColourJunction& operator=( const ColourJunction& ju) {
    this->Junction::operator=(ju);
    for(int i = 0;i < 3;++i) {
      dips[i] = ju.dips[i]; dipsOrig[i] = ju.dipsOrig[i];
    }
    return (*this);
  }

  ColourDipole * getColDip(int i) {return dips[i];}
  void setColDip(int i, ColourDipole * dip) {dips[i] = dip;}
  ColourDipole * dips[3];
  ColourDipole * dipsOrig[3];
  void print();

};

//==========================================================================

// TrialReconnection class.

//--------------------------------------------------------------------------

class TrialReconnection {

public:

  TrialReconnection(ColourDipole* dip1In = 0, ColourDipole* dip2In = 0,
    ColourDipole* dip3In = 0, ColourDipole* dip4In = 0, int modeIn = 0,
    double lambdaDiffIn = 0) {
    dips.push_back(dip1In); dips.push_back(dip2In);
    dips.push_back(dip3In); dips.push_back(dip4In);
    mode = modeIn; lambdaDiff = lambdaDiffIn;
  }

  void list() {
    cout << "mode: " << mode << " " << "lambdaDiff: " << lambdaDiff << endl;
    for (int i = 0;i < int(dips.size()) && dips[i] != 0;++i) {
      cout << "   "; dips[i]->print(); }
  }

  vector<ColourDipole*> dips;
  int mode;
  double lambdaDiff;

};

//==========================================================================

// ColourParticle class.

//--------------------------------------------------------------------------

class ColourParticle : public Particle {

public:

 ColourParticle(const Particle& ju) : Particle(ju) {}

  vector<vector<ColourDipole *> > dips;
  vector<bool> colEndIncluded, acolEndIncluded;
  vector<ColourDipole *> activeDips;
  bool isJun;
  int junKind;

  // Printing functions, intended for debugging.
  void  list();
  void listActiveDips();
  void print();

};

//==========================================================================

// The ColourReconnection class handles the colour reconnection.

//--------------------------------------------------------------------------

class ColourReconnection {

public:

  // Constructor
  ColourReconnection() {}

  // Initialization.
  bool init( Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn,
    ParticleData* particleDataPtrIn, BeamParticle* beamAPtrIn,
    BeamParticle* beamBPtrIn, PartonSystems* partonSystemsPtrIn);

  // New beams possible for handling of hard diffraction.
  void reassignBeamPtrs( BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn)
    {beamAPtr = beamAPtrIn; beamBPtr = beamBPtrIn;}

  // Do colour reconnection for current event.
  bool next( Event & event, int oldSize);

private:

  // Constants: could only be changed in the code itself.
  static const double MINIMUMGAIN, MINIMUMGAINJUN, HBAR;
  static const int MAXRECONNECTIONS;

  // Variables needed.
  bool   allowJunctions, sameNeighbourCol, singleReconOnly, lowerLambdaOnly;
  int    nSys, nReconCols, swap1, swap2, reconnectMode, flipMode,
         timeDilationMode;
  double eCM, sCM, pT0, pT20Rec, pT0Ref, ecmRef, ecmPow, reconnectRange,
         m0, m0sqr, m2Lambda, fracGluon, dLambdaCut, timeDilationPar,
         timeDilationParGeV, tfrag, blowR, blowT, rHadron, kI;

  // List of current dipoles.
  vector<ColourDipole*> dipoles, usedDipoles;
  vector<ColourJunction> junctions;
  vector<ColourParticle> particles;
  vector<TrialReconnection> junTrials, dipTrials;
  vector<vector<int> > iColJun;
  map<int,double> formationTimes;

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to particle data table.
  ParticleData*  particleDataPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointers to the two incoming beams.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;

  // Pointer to information on subcollision parton locations.
  PartonSystems* partonSystemsPtr;

  // This is only to access the function call junctionRestFrame.
  StringFragmentation stringFragmentation;

  // This class is used to calculate the string length.
  StringLength stringLength;

  // Do colour reconnection for the event using the new model.
  bool nextNew( Event & event, int oldSize);

  // Simple test swap between two dipoles.
  void swapDipoles(ColourDipole* dip1, ColourDipole* dip2, bool back = false);

  // Setup the dipoles.
  void setupDipoles( Event& event, int iFirst = 0);

  // Form pseuparticle of a given dipole (or junction system).
  void makePseudoParticle( ColourDipole* dip, int status,
    bool setupDone = false);

  // Find the indices in the particle list of the junction and also their
  // respectively leg numbers.
  bool getJunctionIndices(ColourDipole* dip, int &iJun, int &i0, int &i1,
    int &i2, int &junLeg0, int &junLeg1, int &junLeg2);

  // Form all possible pseudoparticles.
  void makeAllPseudoParticles(Event & event, int iFirst = 0);

  // Update all colours in the event.
  void updateEvent( Event& event, int iFirst = 0);

  double calculateStringLength( ColourDipole* dip,
    vector<ColourDipole*> & dips);

  // Calculate the string length for two event indices.
  double calculateStringLength( int i, int j);

  // Calculate the length of a single junction
  // given the 3 entries in the particle list.
  double calculateJunctionLength(int i, int j, int k);

  // Calculate the length of a double junction,
  // given the 4 entries in the particle record.
  // First two are expected to be the quarks and second two the anti quarks.
  double calculateDoubleJunctionLength( int i, int j, int k, int l);

  // Find all the particles connected in the junction.
  // If a single junction, the size of iParticles should be 3.
  // For multiple junction structures, the size will increase.
  bool findJunctionParticles( int iJun, vector<int>& iParticles,
    vector<bool> &usedJuns, int &nJuns, vector<ColourDipole*> &dips);

  // Do a single trial reconnection.
  void singleReconnection( ColourDipole* dip1, ColourDipole* dip2);

  // Do a single trial reconnection to form a junction.
  void singleJunction(ColourDipole* dip1, ColourDipole* dip2);

  // Do a single trial reconnection to form a junction.
  void singleJunction(ColourDipole* dip1, ColourDipole* dip2,
    ColourDipole* dip3);

  // Print the chain containing the dipole.
  void listChain(ColourDipole* dip);

  // Print all the chains.
  void listAllChains();

  // Print dipoles, intended for debuggning purposes.
  void listDipoles( bool onlyActive = false, bool onlyReal = false);

  // Print particles, intended for debugging purposes.
  void listParticles();

  // Print junctions, intended for debugging purposes.
  void listJunctions();

  // Check that the current dipole setup is consistent. Debug purpose only.
  void checkDipoles();

  // Check that the current dipole setup is consistent. Debug purpose only.
  void checkRealDipoles(Event& event, int iFirst);

  // Calculate the invariant mass of a dipole.
  double mDip(ColourDipole* dip);

  // Find the neighbour to anti colour side. Return false if the dipole is
  // connected to a junction or the new particle has a junction inside of it.
  bool findAntiNeighbour(ColourDipole*& dip);

  // Find the neighbour to colour side. Return false if the dipole is
  // connected to a junction or the new particle has a junction inside of it.
  bool findColNeighbour(ColourDipole*& dip);

  // Check that trials do not contain junctions / unusable pseudoparticles.
  bool checkJunctionTrials();

  // Store all dipoles connected to the ones used in the junction connection.
  void storeUsedDips(TrialReconnection& trial);

  // Change colour structure to describe the reconnection in juncTrial.
  void doJunctionTrial(Event& event, TrialReconnection& juncTrial);

  // Change colour structure to describe the reconnection in juncTrial.
  void doDipoleTrial(TrialReconnection& trial);

  // Change colour structure if it is three dipoles forming a junction system.
  void doTripleJunctionTrial(Event& event, TrialReconnection& juncTrial);

  // Calculate the difference between the old and new lambda.
  double getLambdaDiff(ColourDipole* dip1,
    ColourDipole* dip2, ColourDipole* dip3, ColourDipole* dip4, int mode);

  // Calculate the difference between the old and new lambda (dipole swap).
  double getLambdaDiff(ColourDipole* dip1, ColourDipole* dip2);

  // Update the list of dipole trial swaps to account for latest swap.
  void updateDipoleTrials();

  // Update the list of dipole trial swaps to account for latest swap.
  void updateJunctionTrials();

  // Check whether up to four dipoles are 'causally' connected.
  bool checkTimeDilation(ColourDipole* dip1 = 0, ColourDipole* dip2 = 0,
    ColourDipole* dip3 = 0, ColourDipole* dip4 = 0);

  // Check whether two four momenta are casually connected.
  bool checkTimeDilation(Vec4 p1, Vec4 p2, double t1, double t2);

  // Find the momentum of the dipole.
  Vec4 getDipoleMomentum(ColourDipole* dip);

  // Find all particles connected to a junction system (particle list).
  void addJunctionIndices(int iSinglePar, vector<int> &iPar,
    vector<int> &usedJuncs);

  // Find all the formation times.
  void setupFormationTimes( Event & event);

  // Get the mass of all partons connected to a junction system (event list).
  double getJunctionMass(Event & event, int col);

  // Find all particles connected to a junction system (event list).
  void addJunctionIndices(Event & event, int iSinglePar,
    vector<int> &iPar, vector<int> &usedJuncs);

  // The old MPI-based scheme.
  bool reconnectMPIs( Event& event, int oldSize);

  // Vectors and methods needed for the new gluon-move model.

  // Array of (indices of) all final coloured particles.
  vector<int> iReduceCol, iExpandCol;

  // Array of all lambda distances between coloured partons.
  int nColMove;
  vector<double> lambdaijMove;

  // Function to return lambda value from array.
  double lambda12Move( int i, int j) {
    int iAC = iReduceCol[i]; int jAC = iReduceCol[j];
    return lambdaijMove[nColMove * min( iAC, jAC) + max( iAC, jAC)];
  }

  // Function to return lambda(i,j) + lambda(i,k) - lambda(j,k).
  double lambda123Move( int i, int j, int k) {
    int iAC = iReduceCol[i]; int jAC = iReduceCol[j]; int kAC = iReduceCol[k];
    return lambdaijMove[nColMove * min( iAC, jAC) + max( iAC, jAC)]
         + lambdaijMove[nColMove * min( iAC, kAC) + max( iAC, kAC)]
         - lambdaijMove[nColMove * min( jAC, kAC) + max( jAC, kAC)];
  }

  // The new gluon-move scheme.
  bool reconnectMove( Event& event, int oldSize);

  // The common part for both Type I and II reconnections in e+e..
  bool reconnectTypeCommon( Event& event, int oldSize);

  // The e+e- type I CR model.
  map<double,pair<int,int> > reconnectTypeI( Event& event,
    vector<vector<ColourDipole> > &dips, Vec4 decays[2]);
  //  bool reconnectTypeI( Event& event, int oldSize);

  // The e+e- type II CR model.
  map<double,pair<int,int> > reconnectTypeII( Event& event,
    vector<vector<ColourDipole> > &dips, Vec4 decays[2]);
  //    bool reconnectTypeII( Event& event, int oldSize);

  // calculate the determinant of 3 * 3 matrix.
  double determinant3(vector<vector< double> >& vec);
};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ColourReconnection_H
