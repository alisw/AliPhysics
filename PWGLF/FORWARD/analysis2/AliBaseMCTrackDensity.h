#ifndef ALIBaseMCTRACKDENSITY_MC
#define ALIBaseMCTRACKDENSITY_MC
#include <TNamed.h>
#include "AliForwardFlowWeights.h"
class TList;
class TH1D;
class TH2D;
class TF1;
class TGraph;
class AliMCEvent;
class AliMCParticle;
class AliTrackReference;
class AliStack;

/**
 * A class to calculate the particle density from track references.
 * This code is used both in AliForwardMCCorrectionsTask and
 * AliBaseMCDensity calculator. 
 * 
 * @par Input: 
 *    - Kinematics
 *    - Track-References
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
class AliBaseMCTrackDensity : public TNamed
{
public:
  /** 
   * Default constructor.  Do not use - for ROOT I/O system use only 
   */
  AliBaseMCTrackDensity();
  /** 
   * Normal constructor 
   * 
   * @param name Not used
   */
  AliBaseMCTrackDensity(const char* name);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliBaseMCTrackDensity(const AliBaseMCTrackDensity& o);
  /** 
   * Assignment operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  AliBaseMCTrackDensity& operator=(const AliBaseMCTrackDensity& o);
  /** 
   * Destructor. 
   */
  virtual ~AliBaseMCTrackDensity() {}

  /** 
   * Set whether to only consider primaries 
   * 
   * @param use If true, consider only primaries
   */
  void SetUseOnlyPrimary(Bool_t use) { fUseOnlyPrimary = use; }
  /** 
   * Whether to `after-burn' the particles with flow 
   * 
   * @param use 
   */  
  void SetUseFlowWeights(Bool_t use) { fUseFlowWeights = use; }
  /** 
   * Set whether to print debug messages.  Please note this will
   * produce a lot of output. 
   * 
   * @param debug Whether to enable debug messages or not 
   */
  void SetDebug(Bool_t debug=true) { fDebug = debug; }
  /** 
   * Define ouputs 
   * 
   * @param list List to add outputs to
   */
  virtual void DefineOutput(TList* list);
  /** 
   * Print information to standard out
   * 
   * @param option Not used
   */
  virtual void Print(Option_t* option="") const;
protected:
  /** 
   * Loops over all the particles in the passed event.  If @a primary
   * is not null, then that histogram is filled with the primary
   * particle information - irrespective of whether the particle
   * actually hits the Base or not.  For each track (primary or
   * secondary, unless only primary information is requested - see
   * SetUseOnlyPrimary) loop over all track references to that
   * particle and check if they come from the Base.  In that case,
   * figure out which strip(s) to assign the track to, and fill the @a
   * hits structure.
   * 
   * @param event    MC event 
   * @param vz       IP z-coordinate
   * @param primary  Primary information, if available. 
   * 
   * @return true 
   */
  Bool_t ProcessTracks(const AliMCEvent&   event, 
		       Double_t            vz,
		       TH2D*               primary);
  /** 
   * Process a single track 
   * 
   * @param particle   Particle 
   * @param mother     Ultimate mother
   * 
   * @return true on succsess 
   */
  Bool_t ProcessTrack(AliMCParticle* particle, 
		      const AliMCParticle* mother);
  /** 
   * Must be defined to return the track-reference ID for this detector
   * 
   * @return Detector id set on track references
   */
  virtual Int_t GetDetectorId() const = 0;
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
  virtual AliTrackReference* ProcessRef(AliMCParticle* particle, 
					const AliMCParticle* mother, 
					AliTrackReference* ref) = 0;
  /** 
   * Called at before loop over track references
   * 
   */
  virtual void BeginTrackRefs() {}
  /** 
   * Check a track reference 
   * 
   * @return true if the track reference should be used
   */
  virtual Bool_t CheckTrackRef(AliTrackReference*) const { return true; }
  /** 
   * Called at before loop over track references
   * 
   */
  virtual void EndTrackRefs(Int_t /*nRefs*/) {}
  /** 
   * Store a particle hit in Base<i>dr</i>[<i>s,t</i>] in @a output
   * 
   * 
   * @param particle  Particle to store
   * @param mother    Ultimate mother of particle 
   * @param ref      Longest track reference
   *
   * @return Weight factor
   */  
  virtual Double_t StoreParticle(AliMCParticle*       particle, 
				 const AliMCParticle* mother,
				 AliTrackReference*   ref) const;
  /** 
   * Get the parameters of this event 
   * 
   * @param event Event
   * 
   * @return true if found, false otherwise 
   */
  Bool_t GetCollisionParameters(const AliMCEvent& event);
  /** 
   * Get incident angle of this track reference
   * 
   * @param ref Track reference
   * 
   * @return incident angle (in radians)
   */
  Double_t GetTrackRefTheta(const AliTrackReference* ref) const;
  /** 
   * Get ultimate mother of a track 
   * 
   * @param iTr   Track number 
   * @param event Event
   * 
   * @return Pointer to mother or null 
   */
  const AliMCParticle* GetMother(Int_t iTr, const AliMCEvent& event) const;
  /** 
   * Calculate flow weight 
   *
   * @param eta  Pseudo rapidity 
   * @param pt   Transverse momemtum 
   * @param phi  Azimuthal angle 
   * @param id   Particle PDG code
   *
   * @return Flow weight for the particle
   */
  Double_t CalculateWeight(Double_t eta, Double_t pt, 
			   Double_t phi, Int_t id) const;

  Bool_t   fUseOnlyPrimary;       // Only use primaries 
  Bool_t   fUseFlowWeights;       // Wether to use flow weights 
  TH2D*    fBinFlow;              // eta,phi bin flow 
  TH2D*    fEtaBinFlow;           // dEta vs eta of strip
  TH2D*    fPhiBinFlow;           // dPhi vs phi of strip
  TH1D*    fNRefs;                // Number of track-references per track
  AliForwardFlowWeights fWeights; // Flow weights
  Double_t fVz;                   // IP z-coordinate of this event
  Double_t fB;                    // Impact parameter of this event
  Double_t fPhiR;                 // Reaction plane  of this event
  Bool_t   fDebug;                // Debug flag

  ClassDef(AliBaseMCTrackDensity,2); // Calculate track-ref density
};

#endif
// Local Variables:
//  mode: C++ 
// End:
