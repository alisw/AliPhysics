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
 * @ingroup pwg2_forward_aod
 */
#include <TNamed.h>
#include <TAxis.h>
class AliESDEvent;
class TH2D;
class TH1D;
class TH1I;
class TH1F;
class TH2F;
class TAxis;
class TList;

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
 * @ingroup pwg2_forward_algo 
 * @ingroup pwg2_forward_aod
 */
class AliFMDEventInspector : public TNamed
{
public:
  /** 
   * Return codes 
   */
  enum ECodes {
    /** all ok */
    kOk = 0,
    /** No ESD event */
    kNoEvent = 0x1, 
    /** No triggers found */
    kNoTriggers = 0x2, 
    /** No SPD data */ 
    kNoSPD = 0x4, 
    /** No FMD data */
    kNoFMD = 0x8, 
    /** No vertex found */
    kNoVertex = 0x10, 
    /** Vertex out of range */
    kBadVertex = 0x20
  };
  /** 
   * Trigger bins 
   */
  enum ETrgBins { 
    kInel, 
    kInelGt0, 
    kNSD, 
    kEmpty, 
    kA, 
    kB, 
    kC, 
    kE,
    kPileUp,
    kMCNSD,
    kOffline
  };
  /** 
   * Collision systems
   */
  enum ECollisionSystem { 
    kUnknown,
    kPP, 
    kPbPb
  };
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
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDEventInspector(const AliFMDEventInspector& o);
  /** 
   * Destructor 
   */
  virtual ~AliFMDEventInspector();
  /** 
   * Assignement operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object
   */
  AliFMDEventInspector& operator=(const AliFMDEventInspector&);

  /** 
   * Initialize the object 
   * 
   * @param vtxAxis Vertex axis in use 
   */
  virtual void Init(const TAxis& vtxAxis);
  /** 
   * Process the event 
   * 
   * @param event     Input event 
   * @param triggers  On return, the triggers fired 
   * @param lowFlux   On return, true if the event is considered a low-flux 
   *                  event (according to the setting of fLowFluxCut) 
   * @param ivz       On return, the found vertex bin (1-based).  A zero
   *                  means outside of the defined vertex range
   * @param vz        On return, the z position of the interaction
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
		 Double_t&          vz,
		 Double_t&          cent,
		 UShort_t&          nClusters);
  /** 
   * Define the output histograms.  These are put in a sub list of the
   * passed list.   The histograms are merged before the parent task calls 
   * AliAnalysisTaskSE::Terminate 
   * 
   * @param dir Directory to add to 
   */
  void DefineOutput(TList* dir);
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
   * Use the first physics vtx code.   
   * 
   * @param use Use it or not 
   */
  void SetUseFirstPhysicsVtx(Bool_t use) {fUseFirstPhysicsVertex = use; }
  /** 
   * Set the debug level.  The higher the value the more output 
   * 
   * @param dbg Debug level 
   */
  void SetDebug(Int_t dbg=1) { fDebug = dbg; }
  /** 
   * Fetch our histograms from the passed list 
   * 
   * @param d             Input
   * @param hEventsTr     On return, pointer to histogram, or null
   * @param hEventsTrVtx  On return, pointer to histogram, or null
   * @param hTriggers     On return, pointer to histogram, or null
   * 
   * @return true on success, false otherwise 
   */
  Bool_t FetchHistograms(const TList* d, 
			 TH1I*& hEventsTr, 
			 TH1I*& hEventsTrVtx, 
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
   * Print information
   * 
   * @param option Not used 
   */
  void Print(Option_t* option="") const;
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
   * @param runNo Run number - read off from ESD event
   */
  virtual void StoreInformation(Int_t runNo);
protected:
  /** 
   * Read the trigger information from the ESD event 
   * 
   * @param esd        ESD event 
   * @param triggers   On return, contains the trigger bits 
   * @param nClusters  On return, number of SPD clusters in @f$ |\eta|<1@f$ 
   * 
   * @return @c true on success, @c false otherwise 
   */
  Bool_t ReadTriggers(const AliESDEvent* esd, UInt_t& triggers, 
		      UShort_t& nClusters);
  /** 
   * Read the vertex information from the ESD event 
   * 
   * @param esd  ESD event 
   * @param vz   On return, the vertex Z position 
   * 
   * @return @c true on success, @c false otherwise 
   */
  Bool_t ReadVertex(const AliESDEvent* esd, Double_t& vz, Double_t& vx, Double_t& vy);
  /** 
   * Read centrality from event 
   * 
   * @param esd  Event 
   * @param cent On return, the centrality or negative if not found
   * @param qual On return, centrality quality flag
   * 
   * @return False on error, true otherwise 
   */
  virtual Bool_t ReadCentrality(const AliESDEvent* esd, Double_t& cent,
				UShort_t& qual) const;

  TH1I*    fHEventsTr;    //! Histogram of events w/trigger
  TH1I*    fHEventsTrVtx; //! Events w/trigger and vertex 
  TH1I*    fHEventsAccepted; //! Events w/trigger and vertex in range 
  TH2D*    fHEventsAcceptedXY; //! XY vtx with trigger and Z vertex in range 
  TH1I*    fHTriggers;    //! Triggers
  TH1I*    fHType;        //! Type (low/high flux) of event
  TH1I*    fHWords;       //! Trigger words 
  TH1F*    fHCent;        //! Centrality 
  TH2F*    fHCentVsQual;  //! Centrality vs quality 
  Int_t    fLowFluxCut;   //  Low flux cut
  Double_t fMaxVzErr;     //  Maximum error on v_z
  TList*   fList;         //! Histogram container 
  UShort_t fEnergy;       // CMS energy (per nucleon pair) [GeV]
  Short_t  fField;        // L3 magnetic field [kG]
  UShort_t fCollisionSystem; //  Collision system
  Int_t    fDebug;        //  Debug level 
  TAxis*   fCentAxis;     // Centrality axis used in histograms
  TAxis    fVtxAxis;
  Bool_t   fUseFirstPhysicsVertex; //Use the vtx code from p+p first physics
  ClassDef(AliFMDEventInspector,3); // Inspect the event 
};

#endif
// Local Variables:
//   mode: C++
// End:



