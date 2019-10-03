// 
// This class inspects the event 
//
#ifndef ALIFMDEVENTINSPECTOR_H
#define ALIFMDEVENTINSPECTOR_H
/**
 * @file   AliFMDEventInspector.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:02:48 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_aod
 */
#include <TNamed.h>
#include <TAxis.h>
#include <TList.h>
#include "AliDisplacedVertexSelection.h"
class AliESDEvent;
class AliOADBPhysicsSelection;
class TH2D;
class TH1D;
class TH1I;
class TH1F;
class TH2F;
class TH2I;
class TAxis;
class TVector3;
// class TList;

/** 
 * This class inspects the event 
 *
 * @par Input:
 *   - AliESDFMD object possibly corrected for sharing
 *
 * @par Output:
 *   - A histogram of v_z of events with triggers. 
 *   - A histogram of v_z of events with vertex and triggers 
 *   - A histogram of trigger counters 
 * 
 * Note, that these are added to the master output list 
 *
 * @par Corrections used: 
 *   - None
 *
 * @ingroup pwglf_forward_algo 
 */
class AliFMDEventInspector : public TNamed
{
public:
  /** 
   * Return codes 
   */
  enum ECodes {
    /** all ok - bin 1 */
    kOk = 0,
    /** No ESD event - bin 2 */
    kNoEvent = 0x1, 
    /** No triggers found - bin 3 */
    kNoTriggers = 0x2, 
    /** No SPD data - bin 4 */ 
    kNoSPD = 0x4, 
    /** No FMD data - bin 5 */
    kNoFMD = 0x8, 
    /** No vertex found - bin 6 */
    kNoVertex = 0x10, 
    /** Vertex out of range - bin 7 */
    kBadVertex = 0x20
  };
  /** 
   * Trigger bins 
   */
  enum ETrgBins { 
    kInel, 
    kInelGt0, 
    kNSD,
    kV0AND, 
    kEmpty, 
    kA, 
    kB, 
    kC, 
    kE,
    kPileUp,
    kMCNSD,
    kSatellite,
    kSpdOutlier,
    kOffline
  };
  /**
   * Ways to determine the primary vertex 
   */
  enum EVtxType { 
    kNormal, 
    kpA2012, 
    kpA2013, 
    kPWGUD, 
    kDisplaced
  };
  /** 
   * Which kinds of pile-up to use 
   */
  enum EPileupType { 
    kSPD          = 0x1, 
    kTracks       = 0x2,
    kOutOfBunch   = 0x4,
    kSPDBins      = 0x8
  };
  /** 
   * Centrality methods 
   */
  enum ECentMethod { 
    kV0Multiplicity, 
    kV0Amplitude, 
    kV0Charge, 
    kFMDRough, 
    kNTracks, 
    kLTracks, 
    kCL0, 
    kCL1, 
    kCND, 
    kNParticles,
    kNeutrons,
    kV0vsFMD, 
    kV0vsNTracks, 
    kZEMvsZDC,
    kDefaultCent 
  };

  /** 
   * Collision systems
   */
  enum ECollisionSystem { 
    kUnknown,
    kPP, 
    kPbPb
  };
  enum EVtxStatus { 
    kVtxOK = 1, 
    kNoVtx, 
    kNoSPDVtx, 
    kFewContrib, 
    kUncertain,
    kNotVtxZ
  }; 
  enum ETrgStatus {
    kNoTrgWords=1,
    kPP2760Fast, 
    kMUON,
    kTriggered,
    kMinBias,           // CINT1
    kMinBiasNoSPD,      // CINT5
    kV0AndTrg,          // CINT7
    kMinBiasAD,         // CINT10
    kHighMult,
    kCentral, 
    kSemiCentral, 
    kDiffractive,
    kADOr,              // 
    kADAnd,             // 
    kUser,
    kOther
  };
  /** 
   * Folder name 
   */
  static const char* fgkFolderName;
  /** 
   * Constructor 
   */
  AliFMDEventInspector();
  /** 
   * Constructor 
   * 
   * @param name Name of object
   */
  AliFMDEventInspector(const char* name);
  /** 
   * Destructor 
   */
  virtual ~AliFMDEventInspector();

  /** 
   * Initialize the object 
   * 
   * @param vtxAxis Vertex axis in use 
   */
  virtual void SetupForData(const TAxis& vtxAxis);
  /** 
   * Process the event 
   * 
   * @param event     Input event 
   * @param triggers  On return, the triggers fired 
   * @param lowFlux   On return, true if the event is considered a low-flux 
   *                  event (according to the setting of fLowFluxCut) 
   * @param ivz       On return, the found vertex bin (1-based).  A zero
   *                  means outside of the defined vertex range
   * @param ip        On return, coordinates of interaction point
   * @param cent      On return, the centrality (in percent) or < 0 
   *                  if not found
   * @param nClusters On return, number of SPD clusters in @f$ |\eta|<1@f$ 
   * 
   * @return 0 (or kOk) on success, otherwise a bit mask of error codes 
   */
  UInt_t Process(const AliESDEvent* event, 
		 UInt_t&            triggers,
		 Bool_t&            lowFlux,
		 UShort_t&          ivz, 
		 TVector3&          ip,
		 Double_t&          cent,
		 UShort_t&          nClusters);
  /** 
   * Define the output histograms.  These are put in a sub list of the
   * passed list.   The histograms are merged before the parent task calls 
   * AliAnalysisTaskSE::Terminate 
   * 
   * @param dir Directory to add to 
   */
  void CreateOutputObjects(TList* dir);
  /** 
   * Set the number of SPD tracklets for which we consider the event a
   * low-flux event or not .
   * 
   * @param c Cut (default 1000)
   */
  void SetLowFluxCut(Int_t c) { fLowFluxCut = c; }
  /** 
   * Set the maximum error on @f$ v_z@f$
   * 
   * @param c Maximum error (in centimeters)
   */
  void SetMaxVzErr(Double_t c=0.1) { fMaxVzErr = c; }
  /** 
   * Set the vertex method to use 
   * 
   * @param t Method 
   */
  void SetVertexMethod(EVtxType t) { fVtxMethod = t; }
  /** 
   * Use the first physics vtx code.   
   * 
   * @param use Use it or not 
   * @deprecated Use SetVertexMethod 
   */
  void SetUseFirstPhysicsVtx(Bool_t use) {
    SetVertexMethod(use ? kPWGUD : kNormal); }
 /** 
   * Use the first physics vtx code.   
   * 
   * @param use Use it or not 
   * @deprecated Use SetVertexMethod 
   */
  void SetpA2012Vtx(Bool_t use) {
    SetVertexMethod(use ? kpA2012 : kNormal); }
 /** 
   * Use the 2012 pA vtx code.   
   * 
   * @param use Use it or not 
   */

  void SetUseV0AndForNSD(Bool_t use=true) {fUseV0AND = use; }
  /** 
   * Set the minimum number of contributors for a 2nd pile-up vertex 
   * 
   * @param nContrib Least number of contributors 
   */
  void SetMinPileupContributors(UShort_t nContrib=3) 
  { 
    fMinPileupContrib = nContrib;
  }
  /** 
   * Set minimum distance from primary vertex to 2nd pile-up vertex 
   * 
   * @param cm Distance (in centimeters)
   */
  void SetMinPileupDistance(Double_t cm=0.8)
  {
    fMinPileupDistance = cm;
  }
  /** 
   * Enable selection of displaced vertices. 
   * 
   * @param use whether to use
   * @deprecated Use SetVertexMethod 
   */
  void SetUseDisplacedVertices(Bool_t use=true) { 
    SetVertexMethod(use ? kDisplaced : kNormal); }  
  Bool_t IsUseDisplacedVertices() const { return AllowDisplaced(); }
  /** 
   * Set the lower centrality cut - if negative, do not use 
   *
   * @deprecated We should accept all events in the AOD pass
   * 
   * @param mincent Lower cut on centrality
   */
  void SetMinCentrality(Double_t mincent=-1.0);
  /** 
   * Set the upper centrality cut - if negative, do not use 
   *
   * @deprecated We should accept all events in the AOD pass
   * 
   * @param maxcent Upper cut on centrality
   */
  void SetMaxCentrality(Double_t maxcent=-1.0);
  /** 
   * Set the centrality method to use.  Possible values are 
   *
   * - VOM      - VZERO multiplicity 
   * - V0A      - VZERO amplitude
   * - V0C      - VZERO charge
   * - FMD      - FMD scaled energy loss
   * - TRK      - Number of tracks
   * - TKL      - Number of tracks
   * - CL0      - 
   * - CL1      - 
   * - CND      - 
   * - NPA      - Neutral particles 
   * - ZNA      - ZDC neutron amplitude 
   * - V0MvsFMD - VZERO versus FMD 
   * - TKLvsVOM - Tracks versus VZERO 
   * - ZEMvsZDC - ZDC 
   * 
   * @param m 
   */
  void SetCentralityMethod(const TString& m) { fCentMethod = m; }
  /** 
   * Set the centrality method
   * 
   * @param m Method identifier 
   */
  void SetCentralityMethod(ECentMethod m);
  /** 
   * Pile-up flags.  A bit mask of EPileupType 
   * - 0x1 Check for number of vertices from SPD 
   * - 0x2 Check for number of vertices from tracks 
   * - 0x4 Check out-of-bunch pile-up 
   * 
   * @param flags A bit mask of EPileupType 
   */
  void SetPileupFlags(UShort_t flags=0x5) { fPileupFlags = flags; }
  /** 
   * Set the debug level.  The higher the value the more output 
   * 
   * @param dbg Debug level 
   */
  void SetDebug(Int_t dbg=1) { fDebug = dbg; }
  /** 
   * Set whether this is MC or not.  Needed by energy loss fitter task
   * that never instantices AliFMDMCEventInspector.  In particular, we
   * need this to make sure we ignore the FAST partition flag in MC
   * for 2.76TeV pp.
   * 
   * @param isMC If true, assume MC input 
   */
  void SetMC(Bool_t isMC=true) { fMC = isMC; }
  Bool_t IsMC() const { return fMC; }
  /** 
   * Fetch our histograms from the passed list 
   * 
   * @param d             Input
   * @param hEventsTr     On return, pointer to histogram, or null
   * @param hEventsTrVtx  On return, pointer to histogram, or null
   * @param hEventsAcc    On return, pointer to histogram, or null
   * @param hTriggers     On return, pointer to histogram, or null
   * 
   * @return true on success, false otherwise 
   */
  Bool_t FetchHistograms(const TList* d, 
			 TH1I*& hEventsTr, 
			 TH1I*& hEventsTrVtx, 
			 TH1I*& hEventsAcc,
			 TH1I*& hTriggers) const;
  /** 
   * Read the collision system, collision energy, and L3 field setting
   * from the ESD
   * 
   * @param esd ESD to get information from 
   * 
   * @return true on success, false 
   */
  Bool_t ReadRunDetails(const AliESDEvent* esd);
  /** 
   * Get the collision system (one of the value in ECollisionSystem)
   * 
   * @return Collision system 
   */
  UShort_t GetCollisionSystem() const { return fCollisionSystem; }
  /** 
   * Get the center of mass energy (per nucleon pair) in GeV 
   * 
   * @return center of mass energy (per nucleon pair) in GeV 
   */
  UShort_t GetEnergy() const { return fEnergy; }
  /** 
   * Get the magnetic field setting of the L3 magnet in kilo Gauss. 
   * 
   * @return The magnetic field setting 
   */
  Short_t  GetField() const { return fField; }
  /** 
   * Get the current run number 
   * 
   * @return The current run number 
   */
  ULong_t GetRunNumber() const { return fRunNumber; }
  /** 
   * Get the production year.
   * 
   * - For real data, this is the year of the data taking 
   * - For MC this is the year the production is anchored to.
   * 
   * @return A two-digit year (post millennium), or -1
   */
  Short_t GetProductionYear() const { return fProdYear; }
  /** 
   * Get the production period. 
   * 
   * - For real data, this is the period identifier of the data taking 
   * - For MC data, this is the period identifier the production was
   *   anchored to
   * 
   * @return Period identifier or null
   */
  Char_t  GetProductionPeriod() const { return fProdLetter; } 
  /** 
   * Get the AliROOT revision used for this production
   * 
   * @return SVN revision number or -1
   */
  Short_t GetProductionRevision() const { return fProdSVN; } 
  /** 
   * Check if the production was an MC production anchored in some
   * real data.
   * 
   * @return true if this (MC) production was anchored 
   */
  Bool_t  IsProductionMC() const { return fProdMC; }
  
  /** 
   * Print information
   * 
   * @param option Not used 
   */
  void Print(Option_t* option="") const;
  // getter for fmincentrality
  // Double_t GetMinCentrality() const { return fMinCent;}
  // gettter for fmaxcentrality 
  // Double_t GetMaxCentrality() const { return fMaxCent;}
  /** 
   * Store information about running conditions in output list 
   * 
   * 3 TNamed objects are defined.   The names are fixed, but the
   * title is a string representation of the information, and the
   * unique ID contains the identifier 
   *
   * - sys   Contains the collision system string and identifier. 
   * - sNN   Contains the center-of-mass energy per nucleon (GeV)
   * - field Contains the L3 magnetic field (kG)
   * - run   Contains the run number
   * 
   */
  virtual void StoreInformation();
  /** 
   * Store - if possible - production information in a sub-list 
   * 
   */
  virtual void StoreProduction();
  /** 
   * Return a string representing the return code 
   * 
   * @param mask Code 
   * 
   * @return String representation 
   */
  static const char* CodeString(UInt_t mask);
protected:
  /** 
   * Copy constructor - not implemented
   * 
   * @param o Object to copy from 
   */
  AliFMDEventInspector(const AliFMDEventInspector& o)
    : TNamed(o),
      fHEventsTr(0), 
      fHEventsTrVtx(0),
      fHEventsAccepted(0),
      fHEventsAcceptedXY(0),
      fHTriggers(0),
      fHType(0),
      fHWords(0),
      fHCent(0),
      fHCentVsQual(0),
      fHStatus(0),
      fHVtxStatus(0),
      fHTrgStatus(0), 
      fHPileup(0),
      fLowFluxCut(1000),
      fMaxVzErr(0.2),
      fList(0),
      fEnergy(0),
      fField(999), 
      fCollisionSystem(kUnknown),
      fDebug(0),
      fCentAxis(0),
      fVtxAxis(10,-10,10),
      fVtxMethod(kNormal),
      // fUseFirstPhysicsVertex(false),
      fUseV0AND(false),
      fPileupFlags(kSPD|kTracks|kOutOfBunch|kSPDBins),
      fMinPileupContrib(3), 
      fMinPileupDistance(0.8),
      // fUseDisplacedVertices(false),
      fDisplacedVertex(),
      fCollWords(),
      fBgWords(),
      fCentMethod("default"),
      fMinCent(-1.0),
      fMaxCent(-1.0),
      // fUsepA2012Vertex(false),
      fRunNumber(0),
      fMC(false),
      fProdYear(-1),
      fProdLetter('?'),
      fProdPass(-1),
      fProdSVN(-1),
      fProdMC(false)
  {}
  /** 
   * Assignement operator - not implemented
   * 
   * @return Reference to this object
   */
  AliFMDEventInspector& operator=(const AliFMDEventInspector&){return *this;}
  Bool_t AllowDisplaced() const { return fVtxMethod == kDisplaced; }
  /** 
   * Cache the configure trigger classes from the physis selection.  
   * 
   * @param cache   where to cache the trigger class. 
   * @param classes List of configured classes. 
   * @param o       Object from OADB with config
   */
  void CacheConfiguredTriggerClasses(TList& cache, 
				     const TList* classes,
				     AliOADBPhysicsSelection* o);
  /** 
   * Read the trigger information from the ESD event 
   * 
   * @param esd        ESD event 
   * @param triggers   On return, contains the trigger bits 
   * @param nClusters  On return, number of SPD clusters in @f$ |\eta|<1@f$ 
   * 
   * @return @c true on success, @c false otherwise 
   */
  Bool_t ReadTriggers(const AliESDEvent& esd, UInt_t& triggers, 
		      UShort_t& nClusters);
  /** 
   * Possible extra check for p-Pb from 2012/13 - require V0AND.
   * 
   * @param esd Event structure 
   * 
   * @return true on success
   */
  Bool_t CheckpAExtraV0(const AliESDEvent& esd) const;
  /** 
   * Check, for the @f$\sqrt{s}=2.76GeV@f$ pp run wether this event
   * was in the fast partition, and if so, filter it out.
   * 
   * @param fastonly Event was in fast-only partition 
   * 
   * @return true if event was in the fast-only partition, for the run
   * period.
   */
  virtual Bool_t CheckFastPartition(bool fastonly) const;
  /** 
   * Check if we have an INEL&gt;0 trigger 
   *  
   * @param esd        Data 
   * @param nClusters  On return, number of clusters
   * @param triggers   Trigger mask to be filled
   * 
   * @return true if we have an INEL&gt;0 trigger 
   */
  virtual Bool_t CheckINELGT0(const AliESDEvent& esd, UShort_t& nClusters, 
			      UInt_t& triggers) const;
  /** 
   * Check if this is a pile-up event
   * 
   * @param esd       Data
   * @param triggers Trigger mask to be filled
   * 
   * @return true if this is a pile-up event
   */
  virtual Bool_t CheckPileup(const AliESDEvent& esd, UInt_t& triggers) const;
  /** 
   * Check for multi-vertex pile-up 
   * 
   * @param esd ESD event 
   * @param checkOtherBC Also check other BC's 
   * 
   * @return true if multiple vertices found 
   */
  virtual Bool_t CheckMultiVertex(const AliESDEvent& esd, 
				  Bool_t checkOtherBC=false) const;
  /** 
   * Check if we have a cosmic trigger.  These should be filtered out. 
   * 
   * @param trigStri Trigger string 
   * 
   * @return true if we have a cosmic trigger
   */
  virtual Bool_t CheckCosmics(const TString& trigStri) const;
  /** 
   * Check if the trigger string corresponds to an empty event 
   * 
   * @param trigStr  Trigger string 
   * @param triggers Trigger mask to be filled
   * 
   * @return true if the trigger string corresponds to an empty event 
   */
  virtual Bool_t CheckEmpty(const TString& trigStr, UInt_t& triggers) const;
  /** 
   * Check the trigger words to see if we have a B, A, C, or E event. 
   * 
   * @param esd       Data
   * @param triggers  Trigger mask to be filled
   * 
   * @return always true
   */
  virtual Bool_t CheckWords(const AliESDEvent& esd, UInt_t& triggers) const;
  /** 
   * Read the vertex information from the ESD event 
   * 
   * @param esd  ESD event 
   * @param ip   On return, the coordinates of the IP
   * 
   * @return @c true on success, @c false otherwise 
   */
  Bool_t ReadVertex(const AliESDEvent& esd, TVector3& ip);
  /** 
   * Check the vertex using the method used in PWG-UD.  That is 
   *
   * - Check we have a vertex and status is OK
   * - Check if we have an SPD vertex and that it's status is OK 
   * - Check if the vertex is from the Z-vertexer, and if it is, 
   *   - Check that the dispersion and resolution is OK 
   * 
   * @param esd Data 
   * @param ip  On return, the coordinates of the IP
   * 
   * @return status
   */
  virtual EVtxStatus CheckPWGUDVertex(const AliESDEvent& esd, 
				      TVector3& ip) const;
  /** 
   * Check the vertex. That is
   *
   * - Check if we have an SPD vertex and that it's status is OK 
   * - Check that we have enough contributors 
   * - Check that the reslution is OK 
   * 
   * @param esd Data 
   * @param ip  On return, the coordinates of the IP
   * 
   * @return status
   */
  virtual EVtxStatus CheckpA2012Vertex(const AliESDEvent& esd, 
				       TVector3& ip) const;
  /** 
   * Check the vertex. That is
   *
   * - Check if we have an normal vertex and that it's status is OK 
   * - Check that we have enough contributors 
   * - Check if we have an SPD vertex
   * - Check that we have enough contributors 
   * - If from Z-vertexer, check the resolution
   * - Check that the two found vertices are within 0.5cm of each other
   * 
   * @param esd Data 
   * @param ip  On return, the coordinates of the IP
   * 
   * @return status
   */
  virtual EVtxStatus CheckpA2013Vertex(const AliESDEvent& esd, 
				       TVector3& ip) const;
  /** 
   * Check the vertex for pA 2012 settings. That is
   *
   * 
   * @param esd Data 
   * @param ip  On return, the coordinates of the IP
   * 
   * @return true if the vertex was found and met the requirements
   */
  virtual EVtxStatus CheckVertex(const AliESDEvent& esd, TVector3& ip) const;
  /** 
   * Read centrality from event 
   * 
   * @param esd  Event 
   * @param cent On return, the centrality or negative if not found
   * @param qual On return, centrality quality flag
   * 
   * @return False on error, true otherwise 
   */
  virtual Bool_t ReadCentrality(const AliESDEvent& esd, Double_t& cent,
				UShort_t& qual) const;

  TH1I*    fHEventsTr;            //! Histogram of events w/trigger
  TH1I*    fHEventsTrVtx;         //! Events w/trigger and vertex 
  TH1I*    fHEventsAccepted;      //! Events w/trigger and vertex in range 
  TH2D*    fHEventsAcceptedXY;    //! XY vtx with trigger and Z vertex in range 
  TH1I*    fHTriggers;            //! Triggers
  TH1I*    fHType;                //! Type (low/high flux) of event
  TH1I*    fHWords;               //! Trigger words 
  TH1F*    fHCent;                //! Centrality 
  TH2F*    fHCentVsQual;          //! Centrality vs quality 
  TH1I*    fHStatus;              //! Event processing status 
  TH1I*    fHVtxStatus;           //! Vertex processing status 
  TH1I*    fHTrgStatus;           //! Trigger processing status 
  TH1I*    fHPileup;              //! Pile-up status 
  Int_t    fLowFluxCut;           //  Low flux cut
  Double_t fMaxVzErr;             //  Maximum error on v_z
  TList*   fList;                 //! Histogram container 
  UShort_t fEnergy;               // CMS energy (per nucleon pair) [GeV]
  Short_t  fField;                // L3 magnetic field [kG]
  UShort_t fCollisionSystem;      //  Collision system
  Int_t    fDebug;                //  Debug level 
  TAxis*   fCentAxis;             // Centrality axis used in histograms
  TAxis    fVtxAxis;              // IP_z Axis 
  EVtxType fVtxMethod;            // Vertex method to use 
  // Bool_t   fUseFirstPhysicsVertex;//Use the vtx code from p+p first physics
  Bool_t   fUseV0AND;             // Use the vtx code from p+p first physics
  UShort_t fPileupFlags;          // Which pileups to use 
  UShort_t fMinPileupContrib;     // Min contributors to 2nd pile-up IP
  Double_t fMinPileupDistance;    // Min distance of 2nd pile-up IP
  // Bool_t   fUseDisplacedVertices; // Analyze displaced vertices?
  AliDisplacedVertexSelection fDisplacedVertex; //Displaced vertex selector
  TList    fCollWords;            //! Configured collision words 
  TList    fBgWords;              //! Configured background words 
  TString  fCentMethod;           // Centrality method
  Double_t fMinCent;              // min centrality
  Double_t fMaxCent;              // max centrailty
  // Bool_t   fUsepA2012Vertex;      // flag to use pA2012 Veretx selection
  ULong_t  fRunNumber;            // Current run number 
  Bool_t   fMC;                   // Is this MC input
  Short_t  fProdYear;             // Production year 
  Char_t   fProdLetter;           // Production letter 
  Short_t  fProdPass;             // Pass number 
  Int_t    fProdSVN;              // AliROOT revision used in production
  Bool_t   fProdMC;               // True if anchor production

  ClassDef(AliFMDEventInspector,16); // Inspect the event 
};

#endif
// Local Variables:
//   mode: C++
// End:



