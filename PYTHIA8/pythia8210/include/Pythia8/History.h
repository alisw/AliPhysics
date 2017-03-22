// History.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
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
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/PartonLevel.h"

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

  // Default constructor
  Clustering(){
    emitted   = 0;
    emittor   = 0;
    recoiler  = 0;
    partner   = 0;
    pTscale   = 0.0;
  }

  // Default destructor
  ~Clustering(){}

  // Copy constructor
  Clustering( const Clustering& inSystem ){
    emitted  = inSystem.emitted;
    emittor  = inSystem.emittor;
    recoiler = inSystem.recoiler;
    partner  = inSystem.partner;
    pTscale  = inSystem.pTscale;
  }

  // Constructor with input
  Clustering( int emtIn, int radIn, int recIn, int partnerIn,
    double pTscaleIn ){
    emitted  = emtIn;
    emittor  = radIn;
    recoiler = recIn;
    partner  = partnerIn;
    pTscale  = pTscaleIn;
  }

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
  History( int depth,
           double scalein,
           Event statein,
           Clustering c,
           MergingHooks* mergingHooksPtrIn,
           BeamParticle beamAIn,
           BeamParticle beamBIn,
           ParticleData* particleDataPtrIn,
           Info* infoPtrIn,
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
                    AlphaStrong * asISR, double RN);

  // For default NL3:
  // Return weight of virtual correction and subtractive for NL3 merging
  double weightLOOP(PartonLevel* trial, double RN);
  // Return O(\alpha_s)-term of CKKWL-weight for NL3 merging
  double weightFIRST(PartonLevel* trial, AlphaStrong* asFSR,
                  AlphaStrong* asISR, double RN, Rndm* rndmPtr );

  // For UMEPS:
  double weight_UMEPS_TREE(PartonLevel* trial, AlphaStrong * asFSR,
                    AlphaStrong * asISR, double RN);
  double weight_UMEPS_SUBT(PartonLevel* trial, AlphaStrong * asFSR,
                    AlphaStrong * asISR, double RN);

  // For unitary NL3:
  double weight_UNLOPS_TREE(PartonLevel* trial, AlphaStrong * asFSR,
                    AlphaStrong * asISR, double RN);
  double weight_UNLOPS_SUBT(PartonLevel* trial, AlphaStrong * asFSR,
                    AlphaStrong * asISR, double RN);
  double weight_UNLOPS_LOOP(PartonLevel* trial, double RN);
  double weight_UNLOPS_SUBTNLO(PartonLevel* trial, double RN);
  double weight_UNLOPS_CORRECTION( int order, PartonLevel* trial,
                  AlphaStrong* asFSR, AlphaStrong* asISR,
                  double RN, Rndm* rndmPtr );

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

  // Function to get the lowest multiplicity reclustered event
  Event lowestMultProc( const double RN) {
    // Return lowest multiplicity state
    return (select(RN)->state);
  }

  // Calculate and return pdf ratio
  double getPDFratio( int side, bool forSudakov, bool useHardPDF,
                      int flavNum, double xNum, double muNum,
                      int flavDen, double xDen, double muDen);

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
  double weightTree(PartonLevel* trial, double as0, double maxscale,
    double pdfScale, AlphaStrong * asFSR, AlphaStrong * asISR,
    double& asWeight, double& pdfWeight);

  // Function to return the \alpha_s-ratio part of the CKKWL weight.
  double weightTreeALPHAS( double as0, AlphaStrong * asFSR,
    AlphaStrong * asISR );
  // Function to return the PDF-ratio part of the CKKWL weight.
  double weightTreePDFs( double maxscale, double pdfScale );
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
  vector<Clustering> findQCDTriple (int EmtTagIn, int colTopIn,
                       const Event& event, vector<int> PosFinalPartn,
                       vector <int> PosInitPartn );

  vector<Clustering> getAllEWClusterings();
  vector<Clustering> getEWClusterings( const Event& event);
  vector<Clustering> findEWTriple( int EmtTagIn, const Event& event,
                       vector<int> PosFinalPartn );

  vector<Clustering> getAllSQCDClusterings();
  vector<Clustering> getSQCDClusterings( const Event& event);
  vector<Clustering> findSQCDTriple (int EmtTagIn, int colTopIn,
                       const Event& event, vector<int> PosFinalPartn,
                       vector <int> PosInitPartn );

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
  Event cluster(const Clustering & inSystem);

  // Function to get the flavour of the radiator before the splitting
  // for clustering
  // IN  int   : Position of the radiator after the splitting, in the event
  //     int   : Position of the emitted after the splitting, in the event
  //     Event : Reference event
  // OUT int   : Flavour of the radiator before the splitting
  int getRadBeforeFlav(const int RadAfter, const int EmtAfter,
        const Event& event);

  // Function to get the colour of the radiator before the splitting
  // for clustering
  // IN  int   : Position of the radiator after the splitting, in the event
  //     int   : Position of the emitted after the splitting, in the event
  //     Event : Reference event
  // OUT int   : Colour of the radiator before the splitting
  int getRadBeforeCol(const int RadAfter, const int EmtAfter,
        const Event& event);

  // Function to get the anticolour of the radiator before the splitting
  // for clustering
  // IN  int   : Position of the radiator after the splitting, in the event
  //     int   : Position of the emitted after the splitting, in the event
  //     Event : Reference event
  // OUT int   : Anticolour of the radiator before the splitting
  int getRadBeforeAcol(const int RadAfter, const int EmtAfter,
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
  bool connectRadiator( Particle& Radiator, const int RadType,
                        const Particle& Recoiler, const int RecType,
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

  double getCurrentZ(const int rad, const int rec, const int emt) const;

  // Function to compute "pythia pT separation" from Particle input
  double pTLund(const Particle& RadAfterBranch,const Particle& EmtAfterBranch,
                const Particle& RecAfterBranch, int ShowerType);

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
  History() {}

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

  // Minimal scalar sum of pT used in Herwig to choose history
  double sumScalarPT;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_History_H
