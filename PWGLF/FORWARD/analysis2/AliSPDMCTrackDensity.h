#ifndef ALISPDMCTRACKDENSITY_MC
#define ALISPDMCTRACKDENSITY_MC
#include <AliBaseMCTrackDensity.h>

/**
 * A class to calculate the particle density from track references.
 * This code is used both in AliForwardMCCorrectionsTask and
 * AliSPDMCDensity calculator. 
 * 
 * @par Input: 
 *    - AliMultiplicity object  - from reconstruction
 *    - Kinematics
 *    - Track-References
 *
 * @par Output: 
 *    - AliESDSPD object  - content is # of track references/strip
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
class AliSPDMCTrackDensity : public AliBaseMCTrackDensity
{
public:
  /** 
   * Default constructor.  Do not use - for ROOT I/O system use only 
   */
  AliSPDMCTrackDensity();
  /** 
   * Normal constructor 
   * 
   * @param name Not used
   */
  AliSPDMCTrackDensity(const char* name);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliSPDMCTrackDensity(const AliSPDMCTrackDensity& o);
  /** 
   * Assignment operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  AliSPDMCTrackDensity& operator=(const AliSPDMCTrackDensity& o);
  /** 
   * Destructor. 
   */
  virtual ~AliSPDMCTrackDensity() {}

  /** 
   * Loops over all the particles in the passed event.  If @a primary
   * is not null, then that histogram is filled with the primary
   * particle information - irrespective of whether the particle
   * actually hits the SPD or not.  For each track (primary or
   * secondary, unless only primary information is requested - see
   * SetUseOnlyPrimary) loop over all track references to that
   * particle and check if they come from the SPD.  In that case,
   * figure out which @f$(\eta,\varphi)@f$-bin to assign the track to,
   * and fill the @a output histogram
   * 
   * @param event    MC event 
   * @param vz       IP z--coordinate
   * @param output   Output of SPD hits
   * @param primary  Primary information, if available. 
   * 
   * @return true 
   */
  Bool_t Calculate(const AliMCEvent&   event, 
		   Double_t            vz,
		   TH2D&               output,
		   TH2D*               primary);
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
  Bool_t CheckTrackRef(AliTrackReference* /*ref*/) const;
  /** 
   * Store a particle hit in Base<i>dr</i>[<i>s,t</i>] in @a output
   * 
   * 
   * @param particle  Particle to store
   * @param mother    Ultimate mother of particle 
   * @param longest   Longest track reference
   * @param vz        Z coordinate of IP
   * @param nC        Total number of track-references in this sector  
   * @param nT 	      Number of distint strips hit in this sector
   * @param output    Output structure 
   */  
  Double_t StoreParticle(AliMCParticle*       particle, 
			 const AliMCParticle* mother,
			 AliTrackReference*   ref) const;
  Double_t fMinR;             // Min radius 
  Double_t fMaxR;             // Max radius 
  Double_t fMinZ;             // Min z
  Double_t fMaxZ;             // Max z
  AliTrackReference* fStored; //! Last stored
  TH2D*              fOutput; //! Output 

  ClassDef(AliSPDMCTrackDensity,3); // Calculate track-ref density
};

#endif
// Local Variables:
//  mode: C++ 
// End:
