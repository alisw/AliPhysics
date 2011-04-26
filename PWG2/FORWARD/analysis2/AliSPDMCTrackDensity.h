#ifndef ALISPDMCTRACKDENSITY_MC
#define ALISPDMCTRACKDENSITY_MC
#include <TNamed.h>
class TList;
class AliMCEvent;
class AliMultiplicity;
class AliMCParticle;
class AliTrackReference;
class TH3D;
class TH2D;
class TH1D;

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
 * @ingroup pwg2_forward_algo
 * @ingroup pwg2_forward_mc
 * @ingroup pwg2_forward_aod
 */
class AliSPDMCTrackDensity : public TNamed
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
   * Set whether to only consider primaries 
   * 
   * @param use If true, consider only primaries
   */
  void SetUseOnlyPrimary(Bool_t use) { fUseOnlyPrimary = use; }
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
  /** 
   * Define ouputs 
   * 
   * @param list List to add outputs to
   */
  void DefineOutput(TList* list);
  
  void Print(Option_t* option="") const;
protected:
  /** 
   * Store a particle hit in FMD<i>dr</i>[<i>s,t</i>] in @a output
   * 
   * 
   * @param particle  Particle to store
   * @param mother    Ultimate mother 
   * @param output    Output structure 
   */  
  void StoreParticle(AliMCParticle*       particle, 
		     const AliMCParticle* mother,
		     Int_t                refNo,
		     Double_t             vz, 
		     TH2D&                output) const;
  /** 
   * Get ultimate mother of a track 
   * 
   * @param iTr   Track number 
   * @param stack Stack 
   * 
   * @return Pointer to mother or null 
   */
  const AliMCParticle* GetMother(Int_t iTr, const AliMCEvent& event) const;
  /** 
   * Get incident angle of this track reference
   * 
   * @param ref Track reference
   * @param vz  Z coordinate of the IP
   * 
   * @return incident angle (in radians)
   */
  Double_t GetTrackRefTheta(const AliTrackReference* ref,
			    Double_t vz) const;
  Bool_t   fUseOnlyPrimary;       // Only use primaries 
  Double_t fMinR;             // Min radius 
  Double_t fMaxR;             // Max radius 
  Double_t fMinZ;             // Min z
  Double_t fMaxZ;             // Max z
  TH2D*    fRZ;               // Location in (r,z)
  TH3D*    fXYZ;              // Location in (x,y,z)
  TH1D*    fNRefs;            // Refs per track 
  TH2D*    fBinFlow;          // eta,phi bin flow 

  ClassDef(AliSPDMCTrackDensity,1); // Calculate track-ref density
};

#endif
// Local Variables:
//  mode: C++ 
// End:
