// History.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file is written by Stefan Prestel.
// It contains the main class for matrix element merging.
// Header file for the Clustering and History classes.

#ifndef Pythia8_History_H
#define Pythia8_History_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonLevel.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SimpleWeakShowerMEs.h"

namespace Pythia8 {

//==========================================================================

// Declaration of Clustering class.
// This class holds information about one radiator, recoiler, emitted system.
// This class is a container class for History class use.

class Clustering {

public:

   // The emitted parton location.
  int emitted;
  // The emittor parton
  int emittor;
  // The recoiler parton.
  int recoiler;
  // The colour connected recoiler (Can be different for ISR)
  int partner;
  // The scale associated with this clustering.
  double pTscale;
  // The flavour of the radiator prior to the emission.
  int flavRadBef;
  // Spin of the radiator (-1 is left handed, +1 is right handed).
  int spinRad;
  // Spin of the emitted  (-1 is left handed, +1 is right handed).
  int spinEmt;
  // Spin of the recoiler (-1 is left handed, +1 is right handed).
  int spinRec;
  // Spin of the radiator before the splitting.
  int spinRadBef;
  // The radiator before the splitting.
  int radBef;
  // The recoiler before the splitting.
  int recBef;

  bool hasProbSet;
  double prob;

  // Default constructor
  Clustering() : emitted(0), emittor(0), recoiler(0), partner(0), pTscale(),
    flavRadBef(0), spinRad(9), spinEmt(9), spinRec(9), spinRadBef(9),
    radBef(0), recBef(0), hasProbSet(false), prob(-1.) {}

  // Default destructor
  ~Clustering(){}

  // Copy constructor
  Clustering( const Clustering& inSystem ) :
    emitted(inSystem.emitted), emittor(inSystem.emittor),
    recoiler(inSystem.recoiler), partner(inSystem.partner),
    pTscale(inSystem.pTscale), flavRadBef(inSystem.flavRadBef),
    spinRad(inSystem.spinRad), spinEmt(inSystem.spinEmt),
    spinRec(inSystem.spinRec), spinRadBef(inSystem.spinRad),
    radBef(inSystem.radBef), recBef(inSystem.recBef),
    hasProbSet(inSystem.hasProbSet), prob(inSystem.prob) {}

  // Constructor with input
  Clustering( int emtIn, int radIn, int recIn, int partnerIn,
    double pTscaleIn, int flavRadBefIn = 0, int spinRadIn = 9,
    int spinEmtIn = 9, int spinRecIn = 9, int spinRadBefIn = 9,
    int radBefIn = 0, int recBefIn = 0, bool hasProbIn = false,
    double probIn = -1.)
    : emitted(emtIn), emittor(radIn), recoiler(recIn), partner(partnerIn),
      pTscale(pTscaleIn), flavRadBef(flavRadBefIn), spinRad(spinRadIn),
      spinEmt(spinEmtIn), spinRec(spinRecIn), spinRadBef(spinRadBefIn),
      radBef(radBefIn), recBef(recBefIn), hasProbSet(hasProbIn),
      prob(probIn) {}

  // Function to return pythia pT scale of clustering
  double pT() const { return pTscale; }

  // print for debug
  void list() const;

};

//==========================================================================

// Declaration of History class
//
// A History object represents an event in a given step in the CKKW-L
// clustering procedure. It defines a tree-like recursive structure,
// where the root node represents the state with n jets as given by
// the matrix element generator, and is characterized by the member
// variable mother being null. The leaves on the tree corresponds to a
// fully clustered paths where the original n-jets has been clustered
// down to the Born-level state. Also states which cannot be clustered
// down to the Born-level are possible - these will be called
// incomplete. The leaves are characterized by the vector of children
// being empty.

class History {

public:

  // The only constructor. Default arguments are used when creating
  // the initial history node. The \a depth is the maximum number of
  // clusterings requested. \a scalein is the scale at which the \a
  // statein was clustered (should be set to the merging scale for the
  // initial history node. \a beamAIn and beamBIn are needed to
  // calcutate PDF ratios, \a particleDataIn to have access to the
  // correct masses of particles. If \a isOrdered is true, the previous
  // clusterings has been ordered. \a is the PDF ratio for this
  // clustering (=1 for FSR clusterings). \a probin is the accumulated
  // probabilities for the previous clusterings, and \ mothin is the
  // previous history node (null for the initial node).
  History( int depthIn,
           double scalein,
           Event statein,
           Clustering c,
           MergingHooks* mergingHooksPtrIn,
           BeamParticle beamAIn,
           BeamParticle beamBIn,
           ParticleData* particleDataPtrIn,
           Info* infoPtrIn,
           PartonLevel* showersIn,
           CoupSM* coupSMPtrIn,
           bool isOrdered,
           bool isStronglyOrdered,
           bool isAllowed,
           bool isNextInInput,
           double probin,
           History * mothin);

  // The destructor deletes each child.
  ~History() {
    for ( int i = 0, N = children.size(); i < N; ++i ) delete children[i];
  }

  void clear() {
    map<double,History *>().swap(paths);
    map<double,History *>().swap(goodBranches);
    map<double,History *>().swap(badBranches);
    for ( int i = 0, N = children.size(); i < N; ++i ) delete children[i];
    vector<History *>().swap(children);
  }

  void clearPaths() {
    bool allClear=true;
    for ( int i = 0, N = children.size(); i < N; ++i )
      if (children[i]->state.size()!= 0) allClear = false;
    if (allClear) state.free();
    return;
  }

  // Function to project paths onto desired paths.
  bool projectOntoDesiredHistories();

  // For CKKW-L, NL3 and UMEPS:
  // In the initial history node, select one of the paths according to
  // the probabilities. This function should be called for the initial
  // history node.
  // IN  trialShower*    : Previously initialised trialShower object,
  //                       to perform trial showering and as
  //                       repository of pointers to initialise alphaS
  //     PartonSystems* : PartonSystems object needed to initialise
  //                      shower objects
  // OUT double         : (Sukadov) , (alpha_S ratios) , (PDF ratios)
  double weightTREE(PartonLevel* trial, AlphaStrong * asFSR,
    AlphaStrong * asISR, AlphaEM * aemFSR, AlphaEM * aemISR, double RN);


  // For default NL3:
  // Return weight of virtual correction and subtractive for NL3 merging
  double weightLOOP(PartonLevel* trial, double RN);
  // Return O(\alpha_s)-term of CKKWL-weight for NL3 merging
  double weightFIRST(PartonLevel* trial, AlphaStrong* asFSR,
    AlphaStrong * asISR, AlphaEM * aemFSR, AlphaEM * aemISR, double RN,
    Rndm* rndmPtr);


  // For UMEPS:
  double weight_UMEPS_TREE(PartonLevel* trial, AlphaStrong * asFSR,
    AlphaStrong * asISR, AlphaEM * aemFSR, AlphaEM * aemISR, double RN);
  double weight_UMEPS_SUBT(PartonLevel* trial, AlphaStrong * asFSR,
    AlphaStrong * asISR, AlphaEM * aemFSR, AlphaEM * aemISR, double RN);


  // For unitary NL3:
  double weight_UNLOPS_TREE(PartonLevel* trial, AlphaStrong * asFSR,
    AlphaStrong * asISR, AlphaEM * aemFSR, AlphaEM * aemISR, double RN,
    int depthIn = -1);
  double weight_UNLOPS_SUBT(PartonLevel* trial, AlphaStrong * asFSR,
    AlphaStrong * asISR, AlphaEM * aemFSR, AlphaEM * aemISR, double RN,
    int depthIn = -1);
  double weight_UNLOPS_LOOP(PartonLevel* trial, AlphaStrong * asFSR,
     AlphaStrong * asISR, AlphaEM * aemFSR, AlphaEM * aemISR, double RN,
     int depthIn = -1);
  double weight_UNLOPS_SUBTNLO(PartonLevel* trial, AlphaStrong * asFSR,
    AlphaStrong * asISR, AlphaEM * aemFSR, AlphaEM * aemISR, double RN,
    int depthIn = -1);
  double weight_UNLOPS_CORRECTION( int order, PartonLevel* trial,
    AlphaStrong* asFSR, AlphaStrong * asISR, AlphaEM * aemFSR,
    AlphaEM * aemISR, double RN, Rndm* rndmPtr );

  // Function to check if any allowed histories were found
  bool foundAllowedHistories() {
    return (children.size() > 0 && foundAllowedPath); }
  // Function to check if any ordered histories were found
  bool foundOrderedHistories() {
    return (children.size() > 0 && foundOrderedPath); }
  // Function to check if any ordered histories were found
  bool foundCompleteHistories() {
    return (children.size() > 0 && foundCompletePath); }

  // Function to set the state with complete scales for evolution
  void getStartingConditions( const double RN, Event& outState );
  // Function to get the state with complete scales for evolution
  bool getClusteredEvent( const double RN, int nSteps, Event& outState);
  // Function to get the first reclustered state above the merging scale.
  bool getFirstClusteredEventAboveTMS( const double RN, int nDesired,
    Event& process, int & nPerformed, bool updateProcess = true );
  // Function to return the depth of the history (i.e. the number of
  // reclustered splittings)
  int nClusterings();
  int nOrdered(double maxscale);
  vector<double> scales();

  // Function to get the lowest multiplicity reclustered event
  Event lowestMultProc( const double RN) {
    // Return lowest multiplicity state
    return (select(RN)->state);
  }

  // Calculate and return pdf ratio
  double getPDFratio( int side, bool forSudakov, bool useHardPDF,
                      int flavNum, double xNum, double muNum,
                      int flavDen, double xDen, double muDen);

  // Envelope function that calls the recursive getWeakProb.
  double getWeakProb();

  // Recursive function that returns the weak probability for the given path.
  // Mode refers to which ME correction to use, 1 = sChannel, 2 = gluon
  // channel, 3 = double quark t-channel, 4 is double quark u-channel.
  double getWeakProb(vector<int>& mode, vector<Vec4>& mom,
     vector<int> fermionLines);

  // return the weak probability of a single reclustering.
  // Mode refers to which ME correction to use, 1 = sChannel, 2 = gluon
  // channel, 3 = double quark t-channel, 4 is double quark u-channel.
  double getSingleWeakProb(vector<int> &mode, vector<Vec4> &mom,
    vector<int> fermionLines);

  // Find map between indecies in the current state and the state after
  // the splitting.
  // NOT IMPLEMENTED FOR MULTIPLE W/Z/GAMMA (NEED TO HAVE A WAY TO
  // IDENTIFY THEM).
  void findStateTransfer(map<int,int> &transfer);

  // Function to print the history that would be chosen from the random number
  // RN. Mainly for debugging.
  void printHistory( const double RN );
  // Function to print the states in a history, starting from the hard process.
  // Mainly for debugging.
  void printStates();

  // Make Pythia class friend
  friend class Pythia;
  // Make Merging class friend
  friend class Merging;

private:

  // Number of trial emission to use for calculating the average number of
  // emissions
  static const int NTRIAL;

  // Function to set all scales in the sequence of states. This is a
  // wrapper routine for setScales and setEventScales methods
  void setScalesInHistory();

  // Function to find the index (in the mother histories) of the
  // child history, thus providing a way access the path from both
  // initial history (mother == 0) and final history (all children == 0)
  // IN vector<int> : The index of each child in the children vector
  //                  of the current history node will be saved in
  //                  this vector
  // NO OUTPUT
  void findPath(vector<int>& out);

  // Functions to set the  parton production scales and enforce
  // ordering on the scales of the respective clusterings stored in
  // the History node:
  // Method will start from lowest multiplicity state and move to
  // higher states, setting the production scales the shower would
  // have used.
  // When arriving at the highest multiplicity, the method will switch
  // and go back in direction of lower states to check and enforce
  // ordering for unordered histories.
  // IN vector<int> : Vector of positions of the chosen child
  //                  in the mother history to allow to move
  //                  in direction initial->final along path
  //    bool        : True: Move in direction low->high
  //                       multiplicity and set production scales
  //                  False: Move in direction high->low
  //                       multiplicity and check and enforce
  //                       ordering
  // NO OUTPUT
  void setScales( vector<int> index, bool forward);

  // Function to find a particle in all higher multiplicity events
  // along the history path and set its production scale to the input
  // scale
  // IN  int iPart       : Parton in refEvent to be checked / rescaled
  //     Event& refEvent : Reference event for iPart
  //     double scale    : Scale to be set as production scale for
  //                       unchanged copies of iPart in subsequent steps
  void scaleCopies(int iPart, const Event& refEvent, double rho);

  // Function to set the OVERALL EVENT SCALES [=state.scale()] to
  // the scale of the last clustering
  // NO INPUT
  // NO OUTPUT
  void setEventScales();

  // Function to print information on the reconstructed scales in one path.
  // For debug only.
  void printScales() { if ( mother ) mother->printScales();
    cout << " size " << state.size() << " scale " << scale << " clusterIn "
      << clusterIn.pT() << " state.scale() " << state.scale() << endl; }

  // Function to project paths onto desired paths.
  bool trimHistories();
  // Function to tag history for removal.
  void remove(){ doInclude = false; }
  // Function to return flag of allowed histories to choose from.
  bool keep(){ return doInclude; }
  // Function implementing checks on a paths, for deciding if the path should
  // be considered valid or not.
  bool keepHistory();
  // Function to check if a path is ordered in evolution pT.
  bool isOrderedPath( double maxscale );

  bool followsInputPath();

  // Function to check if all reconstucted states in a path pass the merging
  // scale cut.
  bool allIntermediateAboveRhoMS( double rhoms, bool good = true );
  // Function to check if any ordered paths were found (and kept).
  bool foundAnyOrderedPaths();

  // Functions to return the z value of the last ISR splitting
  // NO INPUT
  // OUTPUT double : z value of last ISR splitting in history
  double zISR();

  // Functions to return the z value of the last FSR splitting
  // NO INPUT
  // OUTPUT double : z value of last FSR splitting in history
  double zFSR();

  // Functions to return the pT scale of the last ISR splitting
  // NO INPUT
  // OUTPUT double : pT scale of last ISR splitting in history
  double pTISR();

  // Functions to return the pT scale of the last FSR splitting
  // NO INPUT
  // OUTPUT double : pT scale of last FSR splitting in history
  double pTFSR();

  // Functions to return the event with nSteps additional partons
  // INPUT  int   : Number of splittings in the event,
  //                as counted from core 2->2 process
  // OUTPUT Event : event with nSteps additional partons
  Event clusteredState( int nSteps);

  // Function to choose a path from all paths in the tree
  // according to their splitting probabilities
  // IN double    : Random number
  // OUT History* : Leaf of history path chosen
  History * select(double rnd);

  // For a full path, find the weight calculated from the ratio of
  // couplings, the no-emission probabilities, and possible PDF
  // ratios. This function should only be called for the last history
  // node of a full path.
  // IN  TimeShower : Already initialised shower object to be used as
  //                  trial shower
  //     double     : alpha_s value used in ME calculation
  //     double     : Maximal mass scale of the problem (e.g. E_CM)
  //     AlphaStrong: Initialised shower alpha_s object for FSR alpha_s
  //                  ratio calculation
  //     AlphaStrong: Initialised shower alpha_s object for ISR alpha_s
  //                  ratio calculation (can be different from previous)
  double weightTree(PartonLevel* trial, double as0, double aem0,
    double maxscale, double pdfScale, AlphaStrong * asFSR, AlphaStrong * asISR,
    AlphaEM * aemFSR, AlphaEM * aemISR, double& asWeight, double& aemWeight,
    double& pdfWeight);

  // Function to return the \alpha_s-ratio part of the CKKWL weight.
  double weightTreeALPHAS( double as0, AlphaStrong * asFSR,
    AlphaStrong * asISR, int njetMax = -1 );
  // Function to return the \alpha_em-ratio part of the CKKWL weight.
  double weightTreeALPHAEM( double aem0, AlphaEM * aemFSR,
    AlphaEM * aemISR, int njetMax = -1 );
  // Function to return the PDF-ratio part of the CKKWL weight.
  double weightTreePDFs( double maxscale, double pdfScale, int njetMax = -1 );
  // Function to return the no-emission probability part of the CKKWL weight.
  double weightTreeEmissions( PartonLevel* trial, int type, int njetMin,
    int njetMax, double maxscale );

  // Function to generate the O(\alpha_s)-term of the CKKWL-weight
  double weightFirst(PartonLevel* trial, double as0, double muR,
    double maxscale, AlphaStrong * asFSR, AlphaStrong * asISR, Rndm* rndmPtr );

  // Function to generate the O(\alpha_s)-term of the \alpha_s-ratios
  // appearing in the CKKWL-weight.
  double weightFirstALPHAS( double as0, double muR, AlphaStrong * asFSR,
    AlphaStrong * asISR);
  // Function to generate the O(\alpha_em)-term of the \alpha_em-ratios
  // appearing in the CKKWL-weight.
  double weightFirstALPHAEM( double aem0, double muR, AlphaEM * aemFSR,
    AlphaEM * aemISR);
  // Function to generate the O(\alpha_s)-term of the PDF-ratios
  // appearing in the CKKWL-weight.
  double weightFirstPDFs( double as0, double maxscale, double pdfScale,
    Rndm* rndmPtr );
  // Function to generate the O(\alpha_s)-term of the no-emission
  // probabilities appearing in the CKKWL-weight.
  double weightFirstEmissions(PartonLevel* trial, double as0, double maxscale,
    AlphaStrong * asFSR, AlphaStrong * asISR, bool fixpdf, bool fixas );

  // Function to return the default factorisation scale of the hard process.
  double hardFacScale(const Event& event);
  // Function to return the default renormalisation scale of the hard process.
  double hardRenScale(const Event& event);

  // Perform a trial shower using the \a pythia object between
  // maxscale down to this scale and return the corresponding Sudakov
  // form factor.
  // IN  trialShower : Shower object used as trial shower
  //     double     : Maximum scale for trial shower branching
  // OUT  0.0       : trial shower emission outside allowed pT range
  //      1.0       : trial shower successful (any emission was below
  //                  the minimal scale )
  double doTrialShower(PartonLevel* trial, int type, double maxscale,
    double minscale = 0.);

  // Function to bookkeep the indices of weights generated in countEmissions
  bool updateind(vector<int> & ind, int i, int N);

  // Function to count number of emissions between two scales for NLO merging
  vector<double> countEmissions(PartonLevel* trial, double maxscale,
    double minscale, int showerType, double as0, AlphaStrong * asFSR,
    AlphaStrong * asISR, int N, bool fixpdf, bool fixas);

  // Function to integrate PDF ratios between two scales over x and t,
  // where the PDFs are always evaluated at the lower t-integration limit
  double monteCarloPDFratios(int flav, double x, double maxScale,
           double minScale, double pdfScale, double asME, Rndm* rndmPtr);

  // Default: Check if a ordered (and complete) path has been found in
  // the initial node, in which case we will no longer be interested in
  // any unordered paths.
  bool onlyOrderedPaths();

  // Check if a strongly ordered (and complete) path has been found in the
  // initial node, in which case we will no longer be interested in
  // any unordered paths.
  bool onlyStronglyOrderedPaths();

  // Check if an allowed (according to user-criterion) path has been found in
  // the initial node, in which case we will no longer be interested in
  // any forbidden paths.
  bool onlyAllowedPaths();

  // When a full path has been found, register it with the initial
  // history node.
  // IN  History : History to be registered as path
  //     bool    : Specifying if clusterings so far were ordered
  //     bool    : Specifying if path is complete down to 2->2 process
  // OUT true if History object forms a plausible path (eg prob>0 ...)
  bool registerPath(History & l, bool isOrdered, bool isStronglyOrdered,
         bool isAllowed, bool isComplete);

  // For the history-defining state (and if necessary interfering
  // states), find all possible clusterings.
  // NO INPUT
  // OUT vector of all (rad,rec,emt) systems
  vector<Clustering> getAllQCDClusterings();

  // For one given state, find all possible clusterings.
  // IN  Event : state to be investigated
  // OUT vector of all (rad,rec,emt) systems in the state
  vector<Clustering> getQCDClusterings( const Event& event);

  // Function to construct (rad,rec,emt) triples from the event
  // IN  int   : Position of Emitted in event record for which
  //             dipoles should be constructed
  //     int   : Colour topogy to be tested
  //             1= g -> qqbar, causing 2 -> 2 dipole splitting
  //             2= q(bar) -> q(bar) g && g -> gg,
  //              causing a 2 -> 3 dipole splitting
  //     Event : event record to be checked for ptential partners
  // OUT vector of all allowed radiator+recoiler+emitted triples
  vector<Clustering> findQCDTriple (int emtTagIn, int colTopIn,
                       const Event& event, vector<int> posFinalPartn,
                       vector <int> posInitPartn );

  vector<Clustering> getAllEWClusterings();
  vector<Clustering> getEWClusterings( const Event& event);
  vector<Clustering> findEWTripleW( int emtTagIn, const Event& event,
                       vector<int> posFinalPartn, vector<int> posInitPartn );
  vector<Clustering> findEWTripleZ( int emtTagIn, const Event& event,
                       vector<int> posFinalPartn, vector<int> posInitPartn );

  vector<Clustering> getAllSQCDClusterings();
  vector<Clustering> getSQCDClusterings( const Event& event);
  vector<Clustering> findSQCDTriple (int emtTagIn, int colTopIn,
                       const Event& event, vector<int> posFinalPartn,
                       vector <int> posInitPartn );

  // Function to attach (spin-dependent duplicates of) a clustering.
  void attachClusterings (vector<Clustering>& clus, int iEmt, int iRad,
    int iRec, int iPartner, double pT, const Event& event);

  // Calculate and return the probability of a clustering.
  // IN  Clustering : rad,rec,emt - System for which the splitting
  //                  probability should be calcuated
  // OUT splitting probability
  double getProb(const Clustering & SystemIn);

  // Set up the beams (fill the beam particles with the correct
  // current incoming particles) to allow calculation of splitting
  // probability.
  // For interleaved evolution, set assignments dividing PDFs into
  // sea and valence content. This assignment is, until a history path
  // is chosen and a first trial shower performed, not fully correct
  // (since content is chosen form too high x and too low scale). The
  // assignment used for reweighting will be corrected after trial
  // showering
  void setupBeams();

  // Calculate the PDF ratio used in the argument of the no-emission
  // probability.
  double pdfForSudakov();

  // Calculate the hard process matrix element to include in the selection
  // probability.
  double hardProcessME( const Event& event);

  // Perform the clustering of the current state and return the
  // clustered state.
  // IN Clustering : rad,rec,emt triple to be clustered to two partons
  // OUT clustered state
  Event cluster( Clustering & inSystem);

  // Function to get the flavour of the radiator before the splitting
  // for clustering
  // IN  int   : Position of the radiator after the splitting, in the event
  //     int   : Position of the emitted after the splitting, in the event
  //     Event : Reference event
  // OUT int   : Flavour of the radiator before the splitting
  int getRadBeforeFlav(const int radAfter, const int emtAfter,
        const Event& event);

  // Function to get the spin of the radiator before the splitting
  // IN int  : Spin of the radiator after the splitting
  //    int  : Spin of the emitted after the splitting
  // OUT int : Spin of the radiator before the splitting
  int getRadBeforeSpin(const int radAfter, const int emtAfter,
        const int spinRadAfter, const int spinEmtAfter,
        const Event& event);

  // Function to get the colour of the radiator before the splitting
  // for clustering
  // IN  int   : Position of the radiator after the splitting, in the event
  //     int   : Position of the emitted after the splitting, in the event
  //     Event : Reference event
  // OUT int   : Colour of the radiator before the splitting
  int getRadBeforeCol(const int radAfter, const int emtAfter,
        const Event& event);

  // Function to get the anticolour of the radiator before the splitting
  // for clustering
  // IN  int   : Position of the radiator after the splitting, in the event
  //     int   : Position of the emitted after the splitting, in the event
  //     Event : Reference event
  // OUT int   : Anticolour of the radiator before the splitting
  int getRadBeforeAcol(const int radAfter, const int emtAfter,
        const Event& event);

  // Function to get the parton connected to in by a colour line
  // IN  int   : Position of parton for which partner should be found
  //     Event : Reference event
  // OUT int   : If a colour line connects the "in" parton with another
  //             parton, return the Position of the partner, else return 0
  int getColPartner(const int in, const Event& event);

  // Function to get the parton connected to in by an anticolour line
  // IN  int   : Position of parton for which partner should be found
  //     Event : Reference event
  // OUT int   : If an anticolour line connects the "in" parton with another
  //             parton, return the Position of the partner, else return 0
  int getAcolPartner(const int in, const Event& event);

  // Function to get the list of partons connected to the particle
  // formed by reclusterinf emt and rad by colour and anticolour lines
  // IN  int          : Position of radiator in the clustering
  // IN  int          : Position of emitted in the clustering
  //     Event        : Reference event
  // OUT vector<int>  : List of positions of all partons that are connected
  //                    to the parton that will be formed
  //                    by clustering emt and rad.
  vector<int> getReclusteredPartners(const int rad, const int emt,
    const Event& event);

  // Function to extract a chain of colour-connected partons in
  // the event
  // IN     int          : Type of parton from which to start extracting a
  //                       parton chain. If the starting point is a quark
  //                       i.e. flavType = 1, a chain of partons that are
  //                       consecutively connected by colour-lines will be
  //                       extracted. If the starting point is an antiquark
  //                       i.e. flavType =-1, a chain of partons that are
  //                       consecutively connected by anticolour-lines
  //                       will be extracted.
  // IN      int         : Position of the parton from which a
  //                       colour-connected chain should be derived
  // IN      Event       : Refernence event
  // IN/OUT  vector<int> : Partons that should be excluded from the search.
  // OUT     vector<int> : Positions of partons along the chain
  // OUT     bool        : Found singlet / did not find singlet
  bool getColSinglet(const int flavType, const int iParton,
    const Event& event, vector<int>& exclude,
    vector<int>& colSinglet);

  // Function to check that a set of partons forms a colour singlet
  // IN  Event       : Reference event
  // IN  vector<int> : Positions of the partons in the set
  // OUT bool        : Is a colour singlet / is not
  bool isColSinglet( const Event& event, vector<int> system);
  // Function to check that a set of partons forms a flavour singlet
  // IN  Event       : Reference event
  // IN  vector<int> : Positions of the partons in the set
  // IN  int         : Flavour of all the quarks in the set, if
  //                   all quarks in a set should have a fixed flavour
  // OUT bool        : Is a flavour singlet / is not
  bool isFlavSinglet( const Event& event,
    vector<int> system, int flav=0);

  // Function to properly colour-connect the radiator to the rest of
  // the event, as needed during clustering
  // IN  Particle& : Particle to be connected
  //     Particle  : Recoiler forming a dipole with Radiator
  //     Event     : event to which Radiator shall be appended
  // OUT true               : Radiator could be connected to the event
  //     false              : Radiator could not be connected to the
  //                          event or the resulting event was
  //                          non-valid
  bool connectRadiator( Particle& radiator, const int radType,
                        const Particle& recoiler, const int recType,
                        const Event& event );

  // Function to find a colour (anticolour) index in the input event
  // IN  int col       : Colour tag to be investigated
  //     int iExclude1 : Identifier of first particle to be excluded
  //                     from search
  //     int iExclude2 : Identifier of second particle to be excluded
  //                     from  search
  //     Event event   : event to be searched for colour tag
  //     int type      : Tag to define if col should be counted as
  //                      colour (type = 1) [->find anti-colour index
  //                                         contracted with col]
  //                      anticolour (type = 2) [->find colour index
  //                                         contracted with col]
  // OUT int           : Position of particle in event record
  //                     contraced with col [0 if col is free tag]
  int FindCol(int col, int iExclude1, int iExclude2,
              const Event& event, int type, bool isHardIn);

  // Function to in the input event find a particle with quantum
  // numbers matching those of the input particle
  // IN  Particle : Particle to be searched for
  //     Event    : Event to be searched in
  // OUT int      : > 0 : Position of matching particle in event
  //                < 0 : No match in event
  int FindParticle( const Particle& particle, const Event& event,
    bool checkStatus = true );

  // Function to check if rad,emt,rec triple is allowed for clustering
  // IN int rad,emt,rec,partner : Positions (in event record) of the three
  //                      particles considered for clustering, and the
  //                      correct colour-connected recoiler (=partner)
  //    Event event     : Reference event
  bool allowedClustering( int rad, int emt, int rec, int partner,
    const Event& event );

  // Function to check if rad,emt,rec triple is results in
  // colour singlet radBefore+recBefore
  // IN int rad,emt,rec : Positions (in event record) of the three
  //                      particles considered for clustering
  //    Event event     : Reference event
  bool isSinglett( int rad, int emt, int rec, const Event& event );

  // Function to check if event is sensibly constructed: Meaning
  // that all colour indices are contracted and that the charge in
  // initial and final states matches
  // IN  event : event to be checked
  // OUT TRUE  : event is properly construced
  //     FALSE : event not valid
  bool validEvent( const Event& event );

  // Function to check whether two clusterings are identical, used
  // for finding the history path in the mother -> children direction
  bool equalClustering( Clustering clus1 , Clustering clus2 );

  // Chose dummy scale for event construction. By default, choose
  //     sHat     for 2->Boson(->2)+ n partons processes and
  //     M_Boson  for 2->Boson(->)             processes
  double choseHardScale( const Event& event ) const;

  // If the state has an incoming hadron return the flavour of the
  // parton entering the hard interaction. Otherwise return 0
  int getCurrentFlav(const int) const;

   // If the state has an incoming hadron return the x-value for the
   // parton entering the hard interaction. Otherwise return 0.
  double getCurrentX(const int) const;

  double getCurrentZ(const int rad, const int rec, const int emt,
    int idRadBef = 0) const;

  // Function to compute "pythia pT separation" from Particle input
  double pTLund(const Event& event, int radAfterBranch, int emtAfterBranch,
    int recAfterBranch, int showerType, int idRadBef = 0);

  // Function to return the position of the initial line before (or after)
  // a single (!) splitting.
  int posChangedIncoming(const Event& event, bool before);

  // Function to give back the ratio of PDFs and PDF * splitting kernels
  // needed to convert a splitting at scale pdfScale, chosen with running
  // PDFs, to a splitting chosen with PDFs at a fixed scale mu. As needed to
  // properly count emissions.
  double pdfFactor( const Event& event, const int type, double pdfScale,
    double mu );

  // Function giving the product of splitting kernels and PDFs so that the
  // resulting flavour is given by flav. This is used as a helper routine
  // to dgauss
  double integrand(int flav, double x, double scaleInt, double z);

  // Function providing a list of possible new flavours after a w emssion
  // from the input flavour.
  vector<int> posFlavCKM(int flav);

  // Check if the new flavour structure is possible.
  // If clusType is 1 final clustering is assumed, otherwise initial clustering
  // is assumed.
  bool checkFlavour(vector<int>& flavCounts, int flavRad, int flavRadBef,
    int clusType);

  // Check if the weak recoil structure is allowed.
  bool checkWeakRecoils(map<int,int> &allowedRecoils, bool isFirst = false);

  // Find the recoiler for an ISR scattered weak particle.

  int findISRRecoiler();

  // Reverse the boost carried out by the ISR.
  // The three first momenta are given by the ME,
  // the last two are filled in by this function.
  void reverseBoostISR(Vec4& pMother, Vec4& pSister, Vec4& pPartner,
    Vec4& pDaughter, Vec4& pRecoiler, int sign, double eCM, double& phi);

  // Check if an event reclustered into a 2 -> 2 dijet.
  // (Only enabled if W reclustering is used).
  bool isQCD2to2(const Event& event);

  // Check if an event reclustered into a 2 -> 1 Drell-Yan.
  // (Only enabled if W reclustering is used).
  bool isEW2to1(const Event& event);

  // Set selected child indices.
  void setSelectedChild();

  // Setup the weak dipole showers.
  void setupSimpleWeakShower(int nSteps);

  // Update weak dipoles after an emission.
  void transferSimpleWeakShower(vector<int> &mode, vector<Vec4> &mom,
    vector<int> fermionLines, vector<pair<int,int> > &dipoles, int nSteps);

  // Update the weak modes after an emission.
  vector<int> updateWeakModes(vector<int>& weakModes,
    map<int,int>& stateTransfer);

  // Update the fermion lines for the 2 -> 2 process. This is needed for
  // the weak probabilities.
  vector<int> updateWeakFermionLines(vector<int> fermionLines,
    map<int,int>& stateTransfer);

  // Update the list of weak dipoles. This is needed to setup the PS correctly.
  vector<pair<int,int> > updateWeakDipoles(vector<pair<int,int> > &dipoles,
    map<int,int>& stateTransfer);

  // Setup the hard process information needed for calculating weak
  // probabilities and setting up the shower.
  void setupWeakHard(vector<int>& mode, vector<int>& fermionLines,
    vector<Vec4>& mom);

  // Functions to retrieve scale information from external showers.
  double getShowerPluginScale(const Event& event, int rad, int emt, int rec,
    string key, double scalePythia);

  //----------------------------------------------------------------------//
  // Class members.
  //----------------------------------------------------------------------//

  // The state of the event correponding to this step in the
  // reconstruction.
  Event state;

  // The previous step from which this step has been clustered. If
  // null, this is the initial step with the n-jet state generated by
  // the matrix element.
  History * mother;

  // The different steps corresponding to possible clusterings of this
  // state.
  vector<History *> children;

  // After a path is selected, store the child index.
  int selectedChild;

  // The different paths which have been reconstructed indexed with
  // the (incremental) corresponding probability. This map is empty
  // unless this is the initial step (mother == 0).
  map<double,History *> paths;

  // The sum of the probabilities of the full paths. This is zero
  // unless this is the initial step (mother == 0).
  double sumpath;

  // The different allowed paths after projection, indexed with
  // the (incremental) corresponding probability. This map is empty
  // unless this is the initial step (mother == 0).
  map<double,History *> goodBranches, badBranches;
  // The sum of the probabilities of allowed paths after projection. This is
  // zero unless this is the initial step (mother == 0).
  double sumGoodBranches, sumBadBranches;

  // This is set true if an ordered (and complete) path has been found
  // and inserted in paths.
  bool foundOrderedPath;

  // This is set true if a strongly ordered (and complete) path has been found
  // and inserted in paths.
  bool foundStronglyOrderedPath;

  // This is set true if an allowed (according to a user criterion) path has
  // been found and inserted in paths.
  bool foundAllowedPath;

  // This is set true if a complete (with the required number of
  // clusterings) path has been found and inserted in paths.
  bool foundCompletePath;

  // The scale of this step, corresponding to clustering which
  // constructed the corresponding state (or the merging scale in case
  // mother == 0).
  double scale;

  // Flag indicating if a clustering in the construction of all histories is
  // the next clustering demanded by inout clusterings in LesHouches 2.0
  // accord.
  bool nextInInput;

  // The probability associated with this step and the previous steps.
  double prob;

  // The partons and scale of the last clustering.
  Clustering clusterIn;
  int iReclusteredOld, iReclusteredNew;

  // Flag to include the path amongst allowed paths.
  bool doInclude;

  // Pointer to MergingHooks object to get all the settings.
  MergingHooks* mergingHooksPtr;

   // The default constructor is private.
  History() : mother(), selectedChild(), sumpath(), sumGoodBranches(),
    sumBadBranches(), foundOrderedPath(), foundStronglyOrderedPath(),
    foundAllowedPath(), foundCompletePath(), scale(), nextInInput(), prob(),
    iReclusteredOld(), iReclusteredNew(), doInclude(), mergingHooksPtr(),
    particleDataPtr(), infoPtr(), showers(), coupSMPtr(), sumScalarPT() {}

  // The copy-constructor is private.
  History(const History &) {}

  // The assignment operator is private.
  History & operator=(const History &) {
    return *this;
  }

  // BeamParticle to get access to PDFs
  BeamParticle beamA;
  BeamParticle beamB;
  // ParticleData needed to initialise the shower AND to get the
  // correct masses of partons needed in calculation of probability
  ParticleData* particleDataPtr;

  // Info object to have access to all information read from LHE file
  Info* infoPtr;

  // Class for calculation weak shower ME.
  SimpleWeakShowerMEs simpleWeakShowerMEs;

  // Pointer to showers, to simplify external clusterings.
  PartonLevel* showers;

  // Pointer to standard model couplings.
  CoupSM* coupSMPtr;

  // Minimal scalar sum of pT used in Herwig to choose history
  double sumScalarPT;

  double probMaxSave;
  double probMax() {
    if (mother) return mother->probMax();
    return probMaxSave;
  }
  void updateProbMax(double probIn, bool isComplete = false) {
    if (mother) mother->updateProbMax(probIn, isComplete);
    if (!isComplete && !foundCompletePath) return;
    if (abs(probIn) > probMaxSave) probMaxSave = probIn;
  }

  int depth, minDepthSave;
  int minDepth() {
    if ( mother ) return mother->minDepth();
    return minDepthSave;
  }
  void updateMinDepth(int depthIn) {
    if ( mother ) return mother->updateMinDepth(depthIn);
    minDepthSave = (minDepthSave>0) ? min(minDepthSave,depthIn) : depthIn;
  }

  int nMaxOrd;
  int nMaxOrdered() {
    if ( mother ) return mother->nMaxOrdered();
    return nMaxOrd;
  }
  void updateNmaxOrdered(int nord) {
    if ( mother ) mother->updateNmaxOrdered(nord);
    nMaxOrd = max(nMaxOrd,nord);
  }
  int maxDepth() {
    if ( mother ) return mother->maxDepth();
    return depth;
  }
  int npaths() {
    if ( mother ) return mother->npaths();
    return paths.size();
  }

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_History_H
