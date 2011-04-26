#ifndef ALIFMDMCTRACKDENSITY_MC
#define ALIFMDMCTRACKDENSITY_MC
#include "AliForwardUtil.h"
#include <TNamed.h>
class TList;
class TH1D;
class TH2D;
class AliMCEvent;
class AliESDFMD;
class AliMCParticle;
class AliTrackReference;
class AliStack;

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
 * @ingroup pwg2_forward_algo
 * @ingroup pwg2_forward_mc
 * @ingroup pwg2_forward_aod
 */
class AliFMDMCTrackDensity : public TNamed
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
   * Set whether to only consider primaries 
   * 
   * @param use If true, consider only primaries
   */
  void SetUseOnlyPrimary(Bool_t use) { fUseOnlyPrimary = use; }
  
  /** 
   * Set whether to print debug messages.  Please note this will
   * produce a lot of output. 
   * 
   * @param debug Whether to enable debug messages or not 
   */
  void SetDebug(Bool_t debug=true) { fDebug = debug; }
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
		   Double_t            vz,
		   AliESDFMD&          output,
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
   * @param mother    Ultimate mother of particle 
   * @param longest   Longest track reference
   * @param vz        Z coordinate of IP
   * @param nC        Total number of track-references in this sector  
   * @param nT 	      Number of distint strips hit in this sector
   * @param output    Output structure 
   */  
  void StoreParticle(AliMCParticle* particle, 
		     const AliMCParticle* mother,
		     Int_t          longest,
		     Double_t       vz,
		     UShort_t       nC, 
		     UShort_t       nT,
		     AliESDFMD&     output) const;
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
  /** 
   * Get ultimate mother of a track 
   * 
   * @param iTr   Track number 
   * @param stack Stack 
   * 
   * @return Pointer to mother or null 
   */
  const AliMCParticle* GetMother(Int_t iTr, const AliMCEvent& event) const;
  Bool_t   fUseOnlyPrimary;       // Only use primaries 
  UShort_t fMaxConsequtiveStrips; // Max 'cluster' size
  TH1D*    fNr;                   // Number of track-refs per cluster
  TH1D*    fNt;                   // Size of cluster in strips 
  TH2D*    fBinFlow;              // eta,phi bin flow 
  TH2D*    fEtaBinFlow;           // dEta vs eta of strip
  TH2D*    fPhiBinFlow;           // dPhi vs phi of strip
  Bool_t   fDebug;                // Debug flag

  ClassDef(AliFMDMCTrackDensity,1); // Calculate track-ref density
};

#endif
// Local Variables:
//  mode: C++ 
// End:
