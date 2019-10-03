#ifndef ALIFMDMCTRACKELOSS_MC
#define ALIFMDMCTRACKELOSS_MC
#include "AliForwardUtil.h"
#include "AliBaseMCTrackDensity.h"
#include "AliFMDFloatMap.h"
class TH1D;
class TH2;
class TH2D;
class AliESDFMD;
class TClonesArray;
class TTree;

/**
 * A class to calculate the particle eloss from track references.
 * This code is used both in AliForwardMCCorrectionsTask and
 * AliFMDMCEloss calculator. 
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
class AliFMDMCTrackELoss : public AliBaseMCTrackDensity 
{
public:
  /** 
   * Default constructor.  Do not use - for ROOT I/O system use only 
   */
  AliFMDMCTrackELoss();
  /** 
   * Normal constructor 
   * 
   * @param name Not used
   */
  AliFMDMCTrackELoss(const char* name);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDMCTrackELoss(const AliFMDMCTrackELoss& o){;}
  /** 
   * Assignment operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  AliFMDMCTrackELoss& operator=(const AliFMDMCTrackELoss& o){return *this;}
  /** 
   * Destructor. 
   */
  virtual ~AliFMDMCTrackELoss() {}

  /** 
   * @{ 
   * @name Options 
   */
  /** 
   * Whether to make an output nTree 
   * 
   * @param use If true, make an nTree of hits 
   */
  void SetUseTree(Bool_t use=true) { fUseTree = use; }
  /** 
   * Set maximum number of strips per 'cluster' 
   * 
   * @param n  Maximum number of strips per 'cluster' 
   */
  void SetMaxConsequtiveStrips(UShort_t n) { fMaxConsequtiveStrips = n; }
  /* @} */

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
   * @param ip       IP coordinates
   * @param cent     Event centrality 
   * 
   * @return true 
   */
  Bool_t Calculate(const AliESDFMD&    esd, 
		   const AliMCEvent&   event, 
		   const TVector3&     ip, 
		   Double_t            cent);
  /** 
   * Define ouputs 
   * 
   * @param list List to add outputs to
   */
  void CreateOutputObjects(TList* list);

  /** 
   * @{ 
   * @name Access to cache 
   */
  /** 
   * Clear caches 
   */
  void Clear(Option_t* opt="");
  /**
   * Get reference to cache of eloss per primary 
   *
   * @return Reference to cache of eloss per primary 
   */
  AliFMDFloatMap& GetPrimaries() { return fPrimaries; }
  /**
   * Get constant reference to cache of eloss per primary 
   *
   * @return Constant reference to cache of eloss per primary 
   */
  const AliFMDFloatMap& GetPrimaries() const { return fPrimaries; }
  /**
   * Get reference to cache of eloss per secondary
   *
   * @return Reference to cache of eloss per secondary
   */
  AliFMDFloatMap& GetSecondaries() { return fSecondaries; }
  /**
   * Get constant reference to cache of eloss per secondary
   *
   * @return Constant reference to cache of eloss per secondary
   */
  const AliFMDFloatMap& GetSecondaries() const { return fSecondaries; }
  /**
   * Get reference to cache of eloss for all 
   *
   * @return Reference to cache of eloss for all 
   */
  AliFMDFloatMap& GetAll() { return fAll; }
  /**
   * Get constant reference to cache of eloss for all 
   *
   * @return Constant reference to cache of eloss for all 
   */
  const AliFMDFloatMap& GetAll() const { return fAll; }
  /**
   * Get reference to cache of psuedo-rapidity (@f$ \eta@f$)
   *
   * @return Reference to cache of psuedo-rapidity (@f$ \eta@f$)
   */
  AliFMDFloatMap& GetEta() { return fEta; }
  /**
   * Get constant reference to cache of psuedo-rapidity (@f$ \eta@f$)
   *
   * @return Constant reference to cache of psuedo-rapidity (@f$ \eta@f$)
   */
  const AliFMDFloatMap& GetEta() const { return fEta; }

  /** 
   * Get pointer to NTutple if defined
   * 
   * @return Pointer to nTree or null if not enabled 
   */
  TTree* GetTree() const { return fTree; }
  TClonesArray* GetHits() const { return fHits; }

  TH2* GetBetaGammadEdx() const { return fBetaGammadEdx; }
  TH2* GetBetaGammaEta() const { return fBetaGammaEta; } 
  TH2* GetDEdxEta() const { return fDEdxEta; }
  /* @} */

  /** 
   * Print this task 
   * 
   * @param option Not used 
   */
  void Print(Option_t* option="") const;


  /** 
   * Structure to hold hit information 
   */
  struct Hit : public TObject { 
    Double_t fGamma;
    Double_t fBeta;
    Double_t fEta;
    Double_t fDe;
    Double_t fDx;
    UInt_t   fDetId;
    Int_t    fPdg;
    Bool_t   fPrimary;
    mutable UShort_t fDetector; //! 
    mutable Char_t   fRing;     //! 
    mutable UShort_t fSector;   //!
    mutable UShort_t fStrip;    //! 
    
    Hit();
    void Decode() const; 
    UShort_t D() const { Decode(); return fDetector; }
    Char_t   R() const { Decode(); return fRing;     }
    UShort_t S() const { Decode(); return fSector;   }
    UShort_t T() const { Decode(); return fStrip;    }
    Double_t DeDx() const { return (fDx > 0 ? fDe/fDx : 0); }
    Double_t Eta() const { return fEta; }
    Double_t BetaGamma() const { return fBeta*fGamma; }
    UInt_t AbsPdg() const;
    Bool_t IsPrimary() const { return fPrimary; }
    Bool_t IsPion() const { return AbsPdg() == 211; }
    Bool_t IsKaon() const { return AbsPdg() == 321; }
    Bool_t IsProton() const { return AbsPdg() == 2212; }
    Bool_t IsElectron() const { return AbsPdg() == 11 || AbsPdg() == 13; }

    ClassDef(Hit,1);
  };
  /** 
   * Structure to hold event information 
   */
  struct Event  { 
    Double_t fIpZ;
    Double_t fCent;
  };
    
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
  //__________________________________________________________________
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
    Double_t de;               // Summed energy loss 
    Double_t dx;               // Summed dx 
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
  Event fEvent;
    

  UShort_t   fMaxConsequtiveStrips; // Max 'cluster' size
  Bool_t     fUseTree;             // Whether to make an nTree 
  TClonesArray* fHits; 
  TTree*     fTree;  
  TH1D*      fNr;                   // Number of track-refs per cluster
  TH1D*      fNt;                   // Size of cluster in strips 
  TH1D*      fNc;                   // Number of clusters per track
  TH2D*      fNcr;                  // Number of clusters per track
  TH2*       fBetaGammadEdx;
  TH2*       fBetaGammaEta; 
  TH2*       fDEdxEta;
  mutable AliFMDFloatMap fPrimaries;        // Cache of eloss per primary 
  mutable AliFMDFloatMap fSecondaries;      // Cache of eloss per secondary
  mutable AliFMDFloatMap fAll;              // Cache of eloss for all 
  AliFMDFloatMap fEta;              // Cache of pseudo-rapidity


  ClassDef(AliFMDMCTrackELoss,1); // Calculate track-ref density
};

#endif
// Local Variables:
//  mode: C++ 
// End:
