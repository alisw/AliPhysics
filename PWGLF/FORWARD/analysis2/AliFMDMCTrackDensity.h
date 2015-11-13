#ifndef ALIFMDMCTRACKDENSITY_MC
#define ALIFMDMCTRACKDENSITY_MC
#include "AliForwardUtil.h"
#include "AliBaseMCTrackDensity.h"
class TH1D;
class AliESDFMD;

/**
 * A class to calculate the particle density from track references.
 * This code is used both in AliForwardMCCorrectionsTask and
 * AliFMDMCDensity calculator. 
 * 
 * @par Input: 
 *    - AliESDFMD object  - from reconstruction
 *    - Kinematics
 *    - Track-References
 *
 * @par Output: 
 *    - AliESDFMD object  - content is # of track references/strip
 *
 * @par Corrections used: 
 *    - None
 *
 * @par Histograms: 
 *    - Incident angle vs number of track references
 *    - Incident angle vs number of strips/cluster
 *
 * @ingroup pwglf_forward_algo
 * @ingroup pwglf_forward_mc
 * @ingroup pwglf_forward_aod
 */
class AliFMDMCTrackDensity : public AliBaseMCTrackDensity 
{
public:
  /** 
   * Default constructor.  Do not use - for ROOT I/O system use only 
   */
  AliFMDMCTrackDensity();
  /** 
   * Normal constructor 
   * 
   * @param name Not used
   */
  AliFMDMCTrackDensity(const char* name);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDMCTrackDensity(const AliFMDMCTrackDensity& o);
  /** 
   * Assignment operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  AliFMDMCTrackDensity& operator=(const AliFMDMCTrackDensity& o);
  /** 
   * Destructor. 
   */
  virtual ~AliFMDMCTrackDensity() {}

  /** 
   * Set maximum number of strips per 'cluster' 
   * 
   * @param n  Maximum number of strips per 'cluster' 
   */
  void SetMaxConsequtiveStrips(UShort_t n) { fMaxConsequtiveStrips = n; }
   /** 
   * Set minimum dE for strip to be part of 'cluster' 
   * 
   * @param v  Set minimum dE for strip to be part of 'cluster' 
   */


 void SetLowCutvalue(Double_t v) { fLowCutvalue = v; }
 /** 
   * Loops over all the particles in the passed event.  If @a primary
   * is not null, then that histogram is filled with the primary
   * particle information - irrespective of whether the particle
   * actually hits the FMD or not.  For each track (primary or
   * secondary, unless only primary information is requested - see
   * SetUseOnlyPrimary) loop over all track references to that
   * particle and check if they come from the FMD.  In that case,
   * figure out which strip(s) to assign the track to, and fill the @a
   * hits structure.
   * 
   * @param esd      FMD ESD structure 
   * @param event    MC event 
   * @param vz       IP z-coordinate
   * @param output   Output of FMD hits
   * @param primary  Primary information, if available. 
   * 
   * @return true 
   */
  Bool_t Calculate(const AliESDFMD&    esd, 
		   const AliMCEvent&   event, 
		   const TVector3&     ip, 
		   AliESDFMD&          output,
		   TH2D*               primary);
  /** 
   * Define ouputs 
   * 
   * @param list List to add outputs to
   */
  void CreateOutputObjects(TList* list);
  
  void Print(Option_t* option="") const;
protected:
  /** 
   * Must be defined to return the track-reference ID for this detector
   * 
   * @return Detector id set on track references
   */
  Int_t GetDetectorId() const;
  /** 
   * Process a track reference 
   * 
   * @param particle Particle 
   * @param mother   Ultimate mother (if not primary)
   * @param ref      Reference 
   * 
   * @return 0 if no output should be generated for this reference, or
   * pointer to track-reference to produce output for.
   */
  AliTrackReference* ProcessRef(AliMCParticle* particle, 
				const AliMCParticle* mother, 
				AliTrackReference* ref);
  /** 
   * Called at before loop over track references
   * 
   */
  void BeginTrackRefs();
  /** 
   * Called at before loop over track references
   * 
   * @param nRefs Number of references 
   */
  void EndTrackRefs(Int_t nRefs);
  /** 
   * Store a particle hit in Base<i>dr</i>[<i>s,t</i>] in @a output
   * 
   * 
   * @param particle  Particle to store
   * @param mother    Ultimate mother of particle 
   * @param ref       Longest track reference
   *
   * @return weight
   */  
  Double_t StoreParticle(AliMCParticle*       particle, 
			 const AliMCParticle* mother,
			 AliTrackReference*   ref) const;
  /** 
   * Structure holding the state of the `tracker' 
   * 
   */
  mutable struct State 
  {
    Double_t angle;            // Angle 
    UShort_t oldDetector;      // Last detector
    Char_t   oldRing;          // Last ring
    UShort_t oldSector;        // Last sector
    UShort_t oldStrip;         // Last strip 
    UShort_t startStrip;       // First strip 
    UShort_t nRefs;            // Number of references
    UShort_t nStrips;          // Number of strips 
    UShort_t count;            // Count of hit strips 
    AliTrackReference* longest; //! Longest track through 
    /** 
     * Clear this state
     * 
     * @param alsoCount If true, also clear count 
     */
    void Clear(Bool_t alsoCount=false);
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this object 
     */
    State& operator=(const State& o);
  } fState; //! State 
  
    
  UShort_t   fMaxConsequtiveStrips; // Max 'cluster' size
  Double_t   fLowCutvalue; 	    // cut on minimal value of the dE to be part of cluster 
  TH1D*      fNr;                   // Number of track-refs per cluster
  TH1D*      fNt;                   // Size of cluster in strips 
  TH1D*      fNc;                   // Number of clusters per track
  TH2D*      fNcr;                  // Number of clusters per track
  AliESDFMD* fOutput;               //! Output ESD object

  ClassDef(AliFMDMCTrackDensity,6); // Calculate track-ref density
};

#endif
// Local Variables:
//  mode: C++ 
// End:
