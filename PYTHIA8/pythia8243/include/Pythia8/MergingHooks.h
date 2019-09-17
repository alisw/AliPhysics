// MergingHooks.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file is written by Stefan Prestel.
// Header file to allow user access to program at different stages.
// HardProcess: Container class for the hard process to be merged. Holds the
//              bookkeeping of particles not be be reclustered
// MergingHooks: Steering class for matrix element merging. Some functions can
//               be redefined in a derived class to have access to the merging

#ifndef Pythia8_MergingHooks_H
#define Pythia8_MergingHooks_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"


namespace Pythia8 {

class PartonLevel;

//==========================================================================

// Declaration of hard process class
// This class holds information on the desired hard 2->2 process
// for the merging.
// This class is a container class for History class use.

class HardProcess {

public:

   // Flavour of the first incoming particle
  int hardIncoming1;
  // Flavour of the second incoming particle
  int hardIncoming2;
  // Flavours of the outgoing particles
  vector<int> hardOutgoing1;
  vector<int> hardOutgoing2;
  // Flavour of intermediate bosons in the hard 2->2
  vector<int> hardIntermediate;

  // Current reference event
  Event state;
  // Potential positions of outgoing particles in reference event
  vector<pair<int,int> > PosOutgoing1;
  vector<pair<int,int> > PosOutgoing2;
  // Potential positions of intermediate bosons in reference event
  vector<pair<int,int> > PosIntermediate;

  // Information on merging scale read from LHE file
  double tms;

  // Default constructor
  HardProcess() : hardIncoming1(), hardIncoming2(), tms(){}
  // Default destructor
  virtual ~HardProcess(){}

  // Copy constructor
  HardProcess( const HardProcess& hardProcessIn )
    : state(hardProcessIn.state),
      tms(hardProcessIn.tms) {
    hardIncoming1 = hardProcessIn.hardIncoming1;
    hardIncoming2 = hardProcessIn.hardIncoming2;
    for(int i =0; i < int(hardProcessIn.hardOutgoing1.size());++i)
      hardOutgoing1.push_back( hardProcessIn.hardOutgoing1[i]);
    for(int i =0; i < int(hardProcessIn.hardOutgoing2.size());++i)
      hardOutgoing2.push_back( hardProcessIn.hardOutgoing2[i]);
    for(int i =0; i < int(hardProcessIn.hardIntermediate.size());++i)
      hardIntermediate.push_back( hardProcessIn.hardIntermediate[i]);
    for(int i =0; i < int(hardProcessIn.PosOutgoing1.size());++i)
      PosOutgoing1.push_back( hardProcessIn.PosOutgoing1[i]);
    for(int i =0; i < int(hardProcessIn.PosOutgoing2.size());++i)
      PosOutgoing2.push_back( hardProcessIn.PosOutgoing2[i]);
    for(int i =0; i < int(hardProcessIn.PosIntermediate.size());++i)
      PosIntermediate.push_back( hardProcessIn.PosIntermediate[i]);
  }

  // Constructor with path to LHE file
  HardProcess( string LHEfile, ParticleData* particleData) : hardIncoming1(),
    hardIncoming2(), tms() {
    state = Event();
    state.init("(hard process)", particleData);
    translateLHEFString(LHEfile);
  }

  // Constructor with core process input
  virtual void initOnProcess( string process, ParticleData* particleData);

  // Constructor with path to LHE file input
  void initOnLHEF( string LHEfile, ParticleData* particleData);

  // Function to access the LHE file and read relevant information
  void translateLHEFString( string LHEpath);

  // Function to translate the process string (in MG/ME notation)
  virtual void translateProcessString( string process);

  // Function to clear hard process information
  void clear();

  // Function to check whether the sets of candidates Pos1, Pos2, together
  // with the proposed candidate iPos give an allowed hard process state
  virtual bool allowCandidates(int iPos, vector<pair<int,int> > Pos1,
    vector<pair<int,int> > Pos2, const Event& event);
  // Function to identify the hard subprocess in the current event
  virtual void storeCandidates( const Event& event, string process);
  // Function to check if the particle event[iPos] matches any of
  // the stored outgoing particles of the hard subprocess
  virtual bool matchesAnyOutgoing(int iPos, const Event& event);
  // Function to check if instead of the particle event[iCandidate], another
  // particle could serve as part of the hard process. Assumes that iCandidate
  // is already stored as part of the hard process.
  virtual bool findOtherCandidates(int iPos, const Event& event,
    bool doReplace);
  // Function to exchange a stored hard process candidate with another choice.
  virtual bool exchangeCandidates( vector<int> candidates1,
    vector<int> candidates2, map<int,int> further1, map<int,int> further2);

  // Function to get the number of coloured final state partons in the
  // hard process
  int nQuarksOut();
  // Function to get the number of uncoloured final state particles in the
  // hard process
  int nLeptonOut();
  // Function to get the number of electroweak final state bosons in the
  // hard process
  int nBosonsOut();

  // Function to get the number of coloured initial state partons in the
  // hard process
  int nQuarksIn();
  // Function to get the number of uncoloured initial state particles in the
  // hard process
  int nLeptonIn();
  // Function to report if a resonace decay was found in the 2->2 sub-process
  // of the  current state
  int hasResInCurrent();
  // Function to report the number of resonace decays in the 2->2 sub-process
  // of the  current state
  int nResInCurrent();
  // Function to report if a resonace decay was found in the 2->2 hard process
  bool hasResInProc();
  // Function to print the hard process (for debug)
  void list() const;
  // Function to print the hard process candidates in the
  // Matrix element state (for debug)
  void listCandidates() const;

};

//==========================================================================

// MergingHooks is base class for user input to the merging procedure.

class MergingHooks {

public:

  // Constructor.
  MergingHooks() : useShowerPluginSave(), showers(),
    doUserMergingSave(false),
    doMGMergingSave(false),
    doKTMergingSave(false),
    doPTLundMergingSave(false),
    doCutBasedMergingSave(false), includeMassiveSave(),
    enforceStrongOrderingSave(),
    orderInRapiditySave(), pickByFullPSave(), pickByPoPT2Save(),
    includeRedundantSave(), pickBySumPTSave(), allowColourShufflingSave(),
    resetHardQRenSave(), resetHardQFacSave(), unorderedScalePrescipSave(),
    unorderedASscalePrescipSave(), unorderedPDFscalePrescipSave(),
    incompleteScalePrescipSave(), ktTypeSave(), nReclusterSave(),
    nQuarksMergeSave(), nRequestedSave(), scaleSeparationFactorSave(),
    nonJoinedNormSave(), fsrInRecNormSave(), herwigAcollFSRSave(),
    herwigAcollISRSave(), pT0ISRSave(), pTcutSave(),
    doNL3TreeSave(false),
    doNL3LoopSave(false),
    doNL3SubtSave(false),
    doUNLOPSTreeSave(false),
    doUNLOPSLoopSave(false),
    doUNLOPSSubtSave(false),
    doUNLOPSSubtNLOSave(false),
    doUMEPSTreeSave(false),
    doUMEPSSubtSave(false),
    doEstimateXSection(false), applyVeto(),
    doRemoveDecayProducts(false), nInProcessNow(-1), muMISave(),
    kFactor0jSave(), kFactor1jSave(),
    kFactor2jSave(), tmsValueSave(), tmsValueNow(), DparameterSave(),
    nJetMaxSave(), nJetMaxNLOSave(), nJetMinWTASave(-1),
    doOrderHistoriesSave(true),
    doCutOnRecStateSave(false),
    doWeakClusteringSave(false),
    doSQCDClusteringSave(false), muFSave(), muRSave(), muFinMESave(),
    muRinMESave(),
    doIgnoreEmissionsSave(true),
    doIgnoreStepSave(true), pTsave(), weightCKKWL1Save(), weightCKKWL2Save(),
    nMinMPISave(), weightCKKWLSave(), weightFIRSTSave(), nJetMaxLocal(),
    nJetMaxNLOLocal(),
    hasJetMaxLocal(false),
    includeWGTinXSECSave(false), nHardNowSave(), nJetNowSave(),
    tmsHardNowSave(), tmsNowSave() {
      inputEvent = Event(); resonances.resize(0); infoPtr = 0; settingsPtr = 0;
      particleDataPtr = 0; partonSystemsPtr = 0; useOwnHardProcess = false;
      hardProcess = 0; stopScaleSave= 0.0; }

  // Make History class friend to allow access to advanced switches
  friend class History;
  // Make Pythia class friend
  friend class Pythia;
  // Make PartonLevel class friend
  friend class PartonLevel;
  // Make SpaceShower class friend
  friend class SpaceShower;
  // Make TimeShower class friend
  friend class TimeShower;
  // Make Merging class friend
  friend class Merging;

  //----------------------------------------------------------------------//
  // Functions that allow user interference
  //----------------------------------------------------------------------//

  // Destructor.
  virtual ~MergingHooks();
  // Function encoding the functional definition of the merging scale
  virtual double tmsDefinition( const Event& event){ return event[0].e();}
  // Function to dampen weights calculated from histories with lowest
  // multiplicity reclustered events that do not pass the ME cuts
  virtual double dampenIfFailCuts( const Event& inEvent ) {
    // Dummy statement to avoid compiler warnings
    if(false) cout << inEvent[0].e();
    return 1.;
  }
  // Hooks to disallow states in the construction of all histories, e.g.
  // because jets are below the merging scale or fail the matrix element cuts
  // Function to allow interference in the construction of histories
  virtual bool canCutOnRecState() { return doCutOnRecStateSave; }
  // Function to check reclustered state while generating all possible
  // histories
  // Function implementing check of reclustered events while constructing
  // all possible histories
  virtual bool doCutOnRecState( const Event& event ) {
    // Dummy statement to avoid compiler warnings.
    if(false) cout << event[0].e();
    // Count number of final state partons.
    int nPartons = 0;
    for( int i=0; i < int(event.size()); ++i)
      if(  event[i].isFinal()
        && (event[i].isGluon() || event[i].isQuark()) )
        nPartons++;
    // For gg -> h, allow only histories with gluons in initial state
    if( hasEffectiveG2EW() && nPartons < 2){
      if(event[3].id() != 21 && event[4].id() != 21)
        return true;
    }
    return false;
  }
  // Function to allow not counting a trial emission.
  virtual bool canVetoTrialEmission() { return false;}
  // Function to check if trial emission should be rejected.
  virtual bool doVetoTrialEmission( const Event&, const Event& ) {
    return false; }

  // Function to calculate the hard process matrix element.
  virtual double hardProcessME( const Event& inEvent ) {
    // Dummy statement to avoid compiler warnings.
    if ( false ) cout << inEvent[0].e();
    return 1.; }

  // Functions for internal use inside Pythia source code
  // Initialize.
  virtual void init();
  // Functions for internal use inside Pythia source code
  // Initialize.
  void initPtr( Settings* settingsPtrIn, Info* infoPtrIn,
    ParticleData* particleDataPtrIn, PartonSystems* partonSystemsPtrIn)
    { settingsPtr = settingsPtrIn; infoPtr = infoPtrIn;
      particleDataPtr = particleDataPtrIn;
      partonSystemsPtr = partonSystemsPtrIn;}

  //----------------------------------------------------------------------//
  // Simple output functions
  //----------------------------------------------------------------------//

  // Function returning the value of the merging scale.
  double tms() {
    if(doCutBasedMergingSave) return 0.;
    //else return tmsValueSave;
    else return tmsValueNow;
  }
  double tmsCut() {
    if(doCutBasedMergingSave) return 0.;
    else return tmsValueSave;
  }
  void tms( double tmsIn ) { tmsValueNow = tmsIn; }

  // Function returning the value of the Delta R_{ij} cut for
  // cut based merging scale definition.
  double dRijMS() {
    return ((tmsListSave.size() == 3) ? tmsListSave[0] : 0.);
  }
  // Function returning the value of the pT_{i} cut for
  // cut based merging scale definition.
  double pTiMS() {
    return ((tmsListSave.size() == 3) ? tmsListSave[1] : 0.);
  }
  // Function returning the value of the pT_{i} cut for
  // cut based merging scale definition.
  double QijMS() {
    return ((tmsListSave.size() == 3) ? tmsListSave[2] : 0.);
  }
  // Function returning the value of the maximal number of merged jets.
  int nMaxJets() { return (hasJetMaxLocal) ? nJetMaxLocal : nJetMaxSave;}
  // Function returning the value of the maximal number of merged jets,
  // for which NLO corrections are available.
  int nMaxJetsNLO()
    { return (hasJetMaxLocal) ? nJetMaxNLOLocal : nJetMaxNLOSave;}
  int nMinJetWTA() { return nJetMinWTASave;}
  // Function to return hard process string.
  string getProcessString() { return processNow;}
  // Function to return the number of outgoing partons in the core process
  int nHardOutPartons(){ return hardProcess->nQuarksOut();}
  // Function to return the number of outgoing leptons in the core process
  int nHardOutLeptons(){ return hardProcess->nLeptonOut();}
  // Function to return the number of outgoing electroweak bosons in the core
  // process.
  int nHardOutBosons(){ return hardProcess->nBosonsOut();}
  // Function to return the number of incoming partons (hadrons) in the core
  // process.
  int nHardInPartons(){ return hardProcess->nQuarksIn();}
  // Function to return the number of incoming leptons in the core process.
  int nHardInLeptons(){ return hardProcess->nLeptonIn();}
  // Function to report the number of resonace decays in the 2->2 sub-process
  // of the  current state.
  int nResInCurrent(){ return hardProcess->nResInCurrent();}
  // Function to determine if user defined merging should be applied.
  bool doUserMerging(){ return doUserMergingSave;}
  // Function to determine if automated MG/ME merging should be applied.
  bool doMGMerging() { return doMGMergingSave;}
  // Function to determine if KT merging should be applied.
  bool doKTMerging() { return doKTMergingSave;}
  // Function to determine if PTLund merging should be applied.
  bool doPTLundMerging() { return doPTLundMergingSave;}
  // Function to determine if cut based merging should be applied.
  bool doCutBasedMerging() { return doCutBasedMergingSave;}
  bool doCKKWLMerging() { return (doUserMergingSave || doMGMergingSave
    || doKTMergingSave || doPTLundMergingSave || doCutBasedMergingSave); }
  // Functions to determine if and which part of  UMEPS merging
  // should be applied
  bool doUMEPSTree() { return doUMEPSTreeSave;}
  bool doUMEPSSubt() { return doUMEPSSubtSave;}
  bool doUMEPSMerging() { return (doUMEPSTreeSave || doUMEPSSubtSave);}
  // Functions to determine if and which part of  NL3 merging
  // should be applied
  bool doNL3Tree() { return doNL3TreeSave;}
  bool doNL3Loop() { return doNL3LoopSave;}
  bool doNL3Subt() { return doNL3SubtSave;}
  bool doNL3Merging() { return (doNL3TreeSave || doNL3LoopSave
                             || doNL3SubtSave); }
  // Functions to determine if and which part of UNLOPS merging
  // should be applied
  bool doUNLOPSTree() { return doUNLOPSTreeSave;}
  bool doUNLOPSLoop() { return doUNLOPSLoopSave;}
  bool doUNLOPSSubt() { return doUNLOPSSubtSave;}
  bool doUNLOPSSubtNLO() { return doUNLOPSSubtNLOSave;}
  bool doUNLOPSMerging() { return (doUNLOPSTreeSave || doUNLOPSLoopSave
                             || doUNLOPSSubtSave || doUNLOPSSubtNLOSave); }
  // Return the number clustering steps that have actually been done.
  int nRecluster() { return nReclusterSave;}

  // Return number of requested additional jets on top of the Born process.
  int nRequested() { return nRequestedSave;}

  //----------------------------------------------------------------------//
  // Output functions to analyse/prepare event for merging
  //----------------------------------------------------------------------//

  // Function to check if event contains an emission not present in the hard
  // process.
  bool isFirstEmission(const Event& event);

  // Function to allow effective gg -> EW boson couplings.
  bool hasEffectiveG2EW() {
    if ( getProcessString().compare("pp>h") == 0 ) return true;
    return false; }

  // Function to allow effective gg -> EW boson couplings.
  bool allowEffectiveVertex( vector<int> in, vector<int> out) {
    if ( getProcessString().compare("ta+ta->jj") == 0
      || getProcessString().compare("ta-ta+>jj") == 0 ) {
      int nInFermions(0), nOutFermions(0), nOutBosons(0);
      for (int i=0; i < int(in.size()); ++i)
        if (abs(in[i])<20) nInFermions++;
      for (int i=0; i < int(out.size()); ++i) {
        if (abs(out[i])<20) nOutFermions++;
        if (abs(out[i])>20) nOutBosons++;
      }
      return (nInFermions%2==0 && nOutFermions%2==0);
    }
    return false;
  }

  // Return event stripped from decay products.
  Event bareEvent( const Event& inputEventIn, bool storeInputEvent );
  // Write event with decay products attached to argument.
  bool reattachResonanceDecays( Event& process );

  // Check if particle at position iPos is part of the hard sub-system
  bool isInHard( int iPos, const Event& event);

  // Function to return the number of clustering steps for the current event
  virtual int getNumberOfClusteringSteps(const Event& event,
    bool resetNjetMax = false);

  //----------------------------------------------------------------------//
  // Functions to steer contruction of histories
  //----------------------------------------------------------------------//

  // Function to force preferred picking of ordered histories. By default,
  // unordered histories will only be considered if no ordered histories
  // were found.
  void orderHistories( bool doOrderHistoriesIn) {
    doOrderHistoriesSave = doOrderHistoriesIn; }
  // Function to force cut on reconstructed states internally, as needed
  // for gg -> Higgs to ensure that e.g. uu~ -> Higgs is not constructed.
  void allowCutOnRecState( bool doCutOnRecStateIn) {
    doCutOnRecStateSave = doCutOnRecStateIn; }

  // Function to allow final state clusterings of weak bosons
  void doWeakClustering( bool doWeakClusteringIn ) {
    doWeakClusteringSave = doWeakClusteringIn; }

  //----------------------------------------------------------------------//
  // Functions used as default merging scales
  //----------------------------------------------------------------------//

  // Function to check if the input particle is a light jet, i.e. should be
  // checked against the merging scale defintion.
  bool checkAgainstCut( const Particle& particle);
  // Function to return the value of the merging scale function in the
  // current event.
  virtual double tmsNow( const Event& event );
  // Find the minimal Lund pT between coloured partons in the event
  double rhoms( const Event& event, bool withColour);
  // Function to calculate the minimal kT in the event
  double kTms(const Event & event);
  // Find the if the event passes the Delta R_{ij}, pT_{i} and Q_{ij} cuts on
  // the matrix element, as needed for cut-based merging scale definition
  double cutbasedms( const Event& event );

  //----------------------------------------------------------------------//
  // Functions to steer shower evolution (public to allow for PS plugin)
  //----------------------------------------------------------------------//

  // Flag to indicate trial shower usage.
  void doIgnoreEmissions( bool doIgnoreIn ) {
    doIgnoreEmissionsSave = doIgnoreIn;
  }
  // Function to allow not counting a trial emission.
  virtual bool canVetoEmission() { return !doIgnoreEmissionsSave; }
  // Function to check if emission should be rejected.
  virtual bool doVetoEmission( const Event& );

  //----------------------------------------------------------------------//
  // Functions used as clusterings / probabilities
  //----------------------------------------------------------------------//

  bool useShowerPluginSave;
  virtual bool useShowerPlugin() { return useShowerPluginSave; }

  //----------------------------------------------------------------------//
  // Functions to retrieve if merging weight should countin the internal
  // cross section and the event weight.
  //----------------------------------------------------------------------//

  bool includeWGTinXSEC() { return includeWGTinXSECSave;}

  //----------------------------------------------------------------------//
  // Functions to retrieve event veto information
  //----------------------------------------------------------------------//

  int nHardNow()      { return nHardNowSave; }
  double tmsHardNow() { return tmsHardNowSave; }
  int nJetsNow()      { return nJetNowSave; }
  double tmsNow()     { return tmsNowSave;}

  void setHardProcessPtr(HardProcess* hardProcIn) { hardProcess = hardProcIn; }

  //----------------------------------------------------------------------//
  // The members, switches etc.
  //----------------------------------------------------------------------//

  // Helper class doing all the core process book-keeping
  bool useOwnHardProcess;
  HardProcess* hardProcess;

  // Pointer to various information on the generation.
  Info*          infoPtr;

  Settings* settingsPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

  // Pointer to the particle systems.
  PartonSystems* partonSystemsPtr;

  PartonLevel* showers;
  void setShowerPointer(PartonLevel* psIn ) {showers = psIn;}

  // AlphaS objects for alphaS reweighting use
  AlphaStrong AlphaS_FSRSave;
  AlphaStrong AlphaS_ISRSave;
  AlphaEM AlphaEM_FSRSave;
  AlphaEM AlphaEM_ISRSave;

  // Saved path to LHE file for more automated merging
  string lheInputFile;

  // Flags for merging procedure definition.
  bool   doUserMergingSave, doMGMergingSave, doKTMergingSave,
         doPTLundMergingSave, doCutBasedMergingSave,
         includeMassiveSave, enforceStrongOrderingSave, orderInRapiditySave,
         pickByFullPSave, pickByPoPT2Save, includeRedundantSave,
         pickBySumPTSave, allowColourShufflingSave, resetHardQRenSave,
         resetHardQFacSave;
  int    unorderedScalePrescipSave, unorderedASscalePrescipSave,
         unorderedPDFscalePrescipSave, incompleteScalePrescipSave,
         ktTypeSave, nReclusterSave, nQuarksMergeSave, nRequestedSave;

  double scaleSeparationFactorSave, nonJoinedNormSave,
         fsrInRecNormSave, herwigAcollFSRSave, herwigAcollISRSave,
         pT0ISRSave, pTcutSave;
  bool   doNL3TreeSave, doNL3LoopSave, doNL3SubtSave;
  bool   doUNLOPSTreeSave, doUNLOPSLoopSave, doUNLOPSSubtSave,
         doUNLOPSSubtNLOSave;
  bool   doUMEPSTreeSave, doUMEPSSubtSave;

  // Flag to only do phase space cut, rejecting events below the tms cut.
  bool   doEstimateXSection;

  bool applyVeto;

  // Save input event in case decay products need to be detached.
  Event inputEvent;
  vector< pair<int,int> > resonances;
  bool doRemoveDecayProducts;
  int nInProcessNow;

  // Starting scale for attaching MPI.
  double muMISave;

  // Precalculated K-factors.
  double kFactor0jSave;
  double kFactor1jSave;
  double kFactor2jSave;

  // Saved members.
  double tmsValueSave, tmsValueNow, DparameterSave;
  int nJetMaxSave;
  int nJetMaxNLOSave;
  int nJetMinWTASave;

  string processSave, processNow;

  // List of cut values to used to define a merging scale. Ordering:
  // 0: DeltaR_{jet_i,jet_j,min}
  // 1: p_{T,jet_i,min}
  // 2: Q_{jet_i,jet_j,min}
  vector<double> tmsListSave;

  // INTERNAL Hooks to allow construction of all histories,
  // including un-ordered ones
  bool doOrderHistoriesSave;

  // INTERNAL Hooks to disallow states in the construction of all histories,
  // e.g. because jets are below the merging scale, of to avoid the
  // construction of uu~ -> Higgs histories.
  bool doCutOnRecStateSave;

  // INTERNAL Hooks to allow clustering W bosons.
  bool doWeakClusteringSave, doSQCDClusteringSave;

  // Store / get first scale in PDF's that Pythia should have used
  double muFSave;
  double muRSave;

  // Store / get factorisation scale used in matrix element calculation.
  double muFinMESave;
  double muRinMESave;

  // Flag to indicate trial shower usage.
  bool doIgnoreEmissionsSave;
  // Flag to indicate if events should be vetoed.
  bool doIgnoreStepSave;
  // Stored weights in case veot needs to be revoked
  double pTsave;
  double weightCKKWL1Save, weightCKKWL2Save;
  int nMinMPISave;
  // Save CKKW-L weight / O(\alpha_s) weight.
  double weightCKKWLSave, weightFIRSTSave;

  // Local copies of nJetMax inputs, if recalculation is necessary.
  int nJetMaxLocal;
  int nJetMaxNLOLocal;
  bool hasJetMaxLocal;

  // Event veto and hard process information, if veto should not applied be
  // directly, but is up to the user.
  bool includeWGTinXSECSave;
  int nHardNowSave, nJetNowSave;
  double tmsHardNowSave, tmsNowSave;

  //----------------------------------------------------------------------//
  // Generic setup functions
  //----------------------------------------------------------------------//

  // Function storing candidates for the hard process in the current event
  // Needed in order not to cluster members of the core process
  void storeHardProcessCandidates(const Event& event){
    hardProcess->storeCandidates(event,getProcessString());
  }

  // Function to set the path to the LHE file, so that more automated merging
  // can be used.
  // Remove "_1.lhe" suffix from LHE file name.
  // This is done so that the HarsProcess class can access both the +0 and +1
  // LHE files to find both the merging scale and the core process string
  // Store.
  void setLHEInputFile( string lheFile) {
    lheInputFile = lheFile.substr(0,lheFile.size()-6); }

  //----------------------------------------------------------------------//
  // Functions for output of class members.
  //----------------------------------------------------------------------//

  // Return AlphaStrong objects
  AlphaStrong* AlphaS_FSR() { return &AlphaS_FSRSave;}
  AlphaStrong* AlphaS_ISR() { return &AlphaS_ISRSave;}
  AlphaEM* AlphaEM_FSR() { return &AlphaEM_FSRSave;}
  AlphaEM* AlphaEM_ISR() { return &AlphaEM_ISRSave;}

  // Functions to return advanced merging switches
  // Include masses in definition of evolution pT and splitting kernels
  bool includeMassive() { return includeMassiveSave;}
  // Prefer strongly ordered histories
  bool enforceStrongOrdering() { return enforceStrongOrderingSave;}
  // Prefer histories ordered in rapidity and evolution pT
  bool orderInRapidity() { return orderInRapiditySave;}
  // Pick history probabilistically by full (correct) splitting probabilities
  bool pickByFull() { return pickByFullPSave;}
  // Pick history probabilistically, with easier form of probabilities
  bool pickByPoPT2() { return pickByPoPT2Save;}
  // Include redundant terms (e.g. PDF ratios) in the splitting probabilities
  bool includeRedundant(){ return includeRedundantSave;}
  // Pick by winner-takes-it-all, with minimum sum of scalar evolution pT
  bool pickBySumPT(){ return pickBySumPTSave;}

  // Prescription for combined scale of unordered emissions
  // 0 : use larger scale
  // 1 : use smaller scale
  int unorderedScalePrescip() { return unorderedScalePrescipSave;}
  // Prescription for combined scale used in alpha_s for unordered emissions
  // 0 : use combined emission scale in alpha_s weight for both (!) splittings
  // 1 : use original reclustered scales of each emission in alpha_s weight
  int unorderedASscalePrescip() { return unorderedASscalePrescipSave;}
  // Prescription for combined scale used in PDF ratios for unordered
  // emissions
  // 0 : use combined emission scale in PDFs for both (!) splittings
  // 1 : use original reclustered scales of each emission in PDF ratiost
  int unorderedPDFscalePrescip() { return unorderedPDFscalePrescipSave;}
  // Prescription for starting scale of incomplete histories
  // 0: use factorization scale
  // 1: use sHat
  // 2: use s
  int incompleteScalePrescip() { return incompleteScalePrescipSave;}

  // Allow swapping one colour index while reclustering
  bool allowColourShuffling() { return allowColourShufflingSave;}

  // Allow use of dynamical renormalisation scale of the core 2-> 2 process.
  bool resetHardQRen() { return resetHardQRenSave; }
  // Allow use of dynamical factorisation scale of the core 2-> 2 process.
  bool resetHardQFac() { return resetHardQFacSave; }

  // Factor by which two scales should differ to be classified strongly
  // ordered.
  double scaleSeparationFactor() { return scaleSeparationFactorSave;}
  // Absolute normalization of splitting probability for non-joined
  // evolution.
  double nonJoinedNorm() { return nonJoinedNormSave;}
  // Absolute normalization of splitting probability for final state
  // splittings with initial state recoiler
  double fsrInRecNorm() { return fsrInRecNormSave;}
  // Factor multiplying scalar evolution pT for FSR splitting, when picking
  // history by minimum scalar pT (see Jonathan Tully's thesis)
  double herwigAcollFSR() { return herwigAcollFSRSave;}
  // Factor multiplying scalar evolution pT for ISR splitting, when picking
  // history by minimum scalar pT (see Jonathan Tully's thesis)
  double herwigAcollISR() { return herwigAcollISRSave;}
  // ISR regularisation scale
  double pT0ISR() { return pT0ISRSave;}
  // Shower cut-off scale
  double pTcut() { return pTcutSave;}

  // MI starting scale.
  void muMI( double mu) { muMISave = mu; }
  double muMI() { return muMISave;}

  // Full k-Factor for current event
  double kFactor(int njet = 0) {
    return (njet == 0) ? kFactor0jSave
          :(njet == 1) ? kFactor1jSave
          : kFactor2jSave;
  }
  // O(\alhpa_s)-term of the k-Factor for current event
  double k1Factor( int njet = 0) {
    return (kFactor(njet) - 1)/infoPtr->alphaS();
  }

  // Function to return if construction of histories is biased towards ordered
  // histories.
  bool orderHistories() { return doOrderHistoriesSave;}

  // INTERNAL Hooks to disallow states in the construction of all histories,
  // e.g. because jets are below the merging scale, of to avoid the
  // construction of uu~ -> Higgs histories.
  bool allowCutOnRecState() { return doCutOnRecStateSave;}

  // INTERNAL Hooks to allow clustering W bosons.
  bool doWeakClustering() { return doWeakClusteringSave;}
  // INTERNAL Hooks to allow clustering clustering of gluons to squarks.
  bool doSQCDClustering() { return doSQCDClusteringSave;}

  // Store / get first scale in PDF's that Pythia should have used
  double muF() { return (muFSave > 0.) ? muFSave : infoPtr->QFac();}
  double muR() { return (muRSave > 0.) ? muRSave : infoPtr->QRen();}
  // Store / get factorisation scale used in matrix element calculation.
  double muFinME() {
    // Start with checking the event attribute called "muf".
    string mus = infoPtr->getEventAttribute("muf2",true);
    double mu  = (mus.empty()) ? 0. : atof((char*)mus.c_str());
    mu = sqrt(mu);
    // Check the scales tag of the event.
    if (infoPtr->scales) mu = infoPtr->getScalesAttribute("muf");
    return (mu > 0.) ? mu : (muFinMESave > 0.) ? muFinMESave : infoPtr->QFac();
  }
  double muRinME() {
    // Start with checking the event attribute called "mur2".
    string mus = infoPtr->getEventAttribute("mur2",true);
    double mu  = (mus.empty()) ? 0. : atof((char*)mus.c_str());
    mu = sqrt(mu);
    // Check the scales tag of the event.
    if (infoPtr->scales) mu = infoPtr->getScalesAttribute("mur");
    return (mu > 0.) ? mu : (muRinMESave > 0.) ? muRinMESave : infoPtr->QRen();
  }


  //----------------------------------------------------------------------//
  // Functions to steer merging code
  //----------------------------------------------------------------------//

  // Flag to indicate if events should be vetoed.
  void doIgnoreStep( bool doIgnoreIn ) { doIgnoreStepSave = doIgnoreIn; }
  // Function to allow event veto.
  virtual bool canVetoStep() { return !doIgnoreStepSave; }

  // Stored weights in case veto needs to be revoked
  void storeWeights( double weight ){ weightCKKWL1Save = weightCKKWL2Save
     = weight; }

  // Function to check event veto.
  virtual bool doVetoStep( const Event& process, const Event& event,
    bool doResonance = false );

  // Set starting scales
  virtual bool setShowerStartingScales( bool isTrial, bool doMergeFirstEmm,
    double& pTscaleIn, const Event& event,
    double& pTmaxFSRIn, bool& limitPTmaxFSRin,
    double& pTmaxISRIn, bool& limitPTmaxISRin,
    double& pTmaxMPIIn, bool& limitPTmaxMPIin );

  // Set shower stopping scale. Necessary to e.g. avoid accumulation of
  // incorrect (low-pT) shower weights through trial showering.
  double stopScaleSave;
  void setShowerStoppingScale( double scale = 0.) { stopScaleSave = scale;}
  double getShowerStoppingScale() { return stopScaleSave;}

  void nMinMPI( int nMinMPIIn ) { nMinMPISave = nMinMPIIn; }
  int nMinMPI() { return nMinMPISave;}

  //----------------------------------------------------------------------//
  // Functions for internal merging scale definions
  //----------------------------------------------------------------------//

  // Function to calculate the kT separation between two particles
  double kTdurham(const Particle& RadAfterBranch,
    const Particle& EmtAfterBranch, int Type, double D );
  // Function to compute "pythia pT separation" from Particle input
  double rhoPythia(const Event& event, int rad, int emt, int rec,
    int ShowerType);

  // Function to find a colour (anticolour) index in the input event,
  // used to find colour-connected recoilers
  int findColour(int col, int iExclude1, int iExclude2,
    const Event& event, int type, bool isHardIn);
  // Function to compute Delta R separation from 4-vector input
  double deltaRij(Vec4 jet1, Vec4 jet2);

  //----------------------------------------------------------------------//
  // Functions for weight management
  //----------------------------------------------------------------------//

  // Function to get the CKKW-L weight for the current event
  double getWeightNLO() { return (weightCKKWLSave - weightFIRSTSave);}
  // Return CKKW-L weight.
  double getWeightCKKWL() { return weightCKKWLSave; }
  // Return O(\alpha_s) weight.
  double getWeightFIRST() { return weightFIRSTSave; }
  // Set CKKW-L weight.
  void setWeightCKKWL( double weightIn){
    weightCKKWLSave = weightIn;
    if ( !includeWGTinXSEC() ) infoPtr->setWeightCKKWL(weightIn); }
  // Set O(\alpha_s) weight.
  void setWeightFIRST( double weightIn){
    weightFIRSTSave = weightIn;
    infoPtr->setWeightFIRST(weightIn); }


  //----------------------------------------------------------------------//
  // Functions and members to store the event veto information
  //----------------------------------------------------------------------//

  // Set CKKWL veto information.
  void setEventVetoInfo(int nJetNowIn, double tmsNowIn) {
    nJetNowSave = nJetNowIn; tmsNowSave = tmsNowIn; }

  // Set the hard process information.
  void setHardProcessInfo(int nHardNowIn, double tmsHardNowIn) {
    nHardNowSave = nHardNowIn; tmsHardNowSave = tmsHardNowIn; }

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_MergingHooks_H
