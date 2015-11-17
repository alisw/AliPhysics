#ifndef ALIFMDESDFIXER_H
#define ALIFMDESDFIXER_H 
#include <TObject.h>
#include <TBits.h>
class TH1;
class TList;
class AliESDFMD;
class TVector3;

/**
 * Class to fix up an ESD object for various small issues. 
 *
 * @par Input:
 *    - AliESDFMD object - from reconstruction 
 *    - Ip z coordinate - from reconstruction 
 *    - Reco noise factor - Assumed noise factor used in reconstruction
 *
 * @par Output:
 *    - The same AliESDFMD object but modified 
 *
 * @par Corrections used
 *    - AliFMDCorrNoiseGain if the noise correction is enabled 
 * 
 * @par Histograms
 *    - Change in @f$\Delta@f$ due to noise correction
 *    - Number of channels declared dead 
 *    - Change in @f$\eta@f$  
 *
 * @ingroup pwglf_forward_algo 
 * @ingroup pwglf_forward_aod
 */
class AliFMDESDFixer : public TObject 
{
public:
  /**
   * Default CTOR - for ROOT I/O - do not use 
   */
  AliFMDESDFixer();

  /** 
   * User CTOR 
   *
   * @param name Dummy argument 
   */
  AliFMDESDFixer(const char* name);
  /** 
   * Get name of object 
   * 
   * @return Always the same 
   */
  const char* GetName() const { return "fmdESDFixer"; }
  /**
   * Create output objects 
   */
  void CreateOutputObjects(TList* l);
  /** 
   * Fix the ESD object 
   * 
   * @param esd ESD object 
   * @param ip  IP coordinates
   */
  void Fix(AliESDFMD& esd, const TVector3& ip);

  /** 
   * @{ 
   * @name Noise suppression fix 
   */
  /** 
   * Set the factor assumed to be used in the reconstruction.  If this
   * is set way high (>=4) then this corrector will effectively be
   * disabled.
   * 
   * @param f 
   */
  void SetRecoNoiseFactor(Int_t f) { fRecoFactor = f; }
  void SetMaxNoiseCorrection(Double_t x) { fMaxNoiseCorr = x; }
  /** 
   * Get the factor assumed in reconstruction pass
   * 
   * @return factor assumed in reconstruction pass
   */
  Int_t GetRecoNoiseFactor() const { return fRecoFactor; }
  /** 
   * Check if we're using the noise correction. 
   * 
   * @return true if fRecoFactor < 4
   */
  Bool_t IsUseNoiseCorrection() const { return fRecoFactor < 4; }
  /** 
   * Find the target noise factor 
   * 
   * @param esd   ESD object. 
   * @param check If true, also check for correction object 
   * 
   * @return Needed noise factor, or 0 or less if no correction is needed
   */
  Int_t FindTargetNoiseFactor(const AliESDFMD& esd, Bool_t check=true) const;
  /* @} */

  /** 
   * @{ 
   * @name Recalculations 
   */
  /** 
   * In case of a displaced vertices recalculate @f$\eta@f$ and angle
   * correction
   * 
   * @param use recalculate or not
   * 
   */
  void SetRecalculateEta(Bool_t use) { fRecalculateEta = use; }
  /* @} */

  /**
   * @{ 
   * @name Special treatmeant of invalid signals 
   */
  /** 
   * Set whether to consider invalid multiplicities as null (or empty)
   * signal. 
   * 
   * @param flag If true, count invalids as empty
   */
  void SetInvalidIsEmpty(Bool_t flag) { fInvalidIsEmpty = flag; }
  /** @} */

  /** 
   * @{ 
   * @name Dead strip handling 
   */
  /** 
   * Add a dead strip
   * 
   * @param d  Detector
   * @param r  Ring 
   * @param s  Sector 
   * @param t  Strip
   */
  void AddDead(UShort_t d, Char_t r, UShort_t s, UShort_t t);
  /** 
   * Add a dead region in a detector ring
   * 
   * @param d   Detector
   * @param r   Ring
   * @param s1  First sector (inclusive)
   * @param s2  Last sector (inclusive)
   * @param t1  First strip (inclusive)
   * @param t2  Last strip (inclusive)
   */
  void AddDeadRegion(UShort_t d, Char_t r, UShort_t s1, UShort_t s2, 
		     UShort_t t1, UShort_t t2);
  /** 
   * Add dead strips from a script.  The script is supposed to accept
   * a pointer to this object (AliFMDSharingFilter) and then call
   * AddDead or AddDeadRegion as needed.
   * 
   * @code 
   * void deadstrips(AliFMDSharingFilter* filter)
   * { 
   *   filter->AddDead(...);
   *   // ... and so on 
   * }
   * @endcode 
   *
   * @param script The script to read dead strips from. 
   */
  void AddDead(const Char_t* script);
  /* @} */

  void Print(Option_t* option="") const;
protected:
  AliFMDESDFixer(const AliFMDESDFixer&);
  AliFMDESDFixer& operator=(const AliFMDESDFixer&);
  
  /** 
   * @{ 
   * @name Worker functions 
   */
  /** 
   * Check if a strip is marked as dead 
   * 
   * @param d   Detector 
   * @param r   Ring
   * @param s   Sector 
   * @param t   Strip 
   * 
   * @return true if dead 
   */
  virtual Bool_t IsDead(UShort_t d, Char_t r, UShort_t s, UShort_t t) const;
  /** 
   * Possibly raise a strip from the dead or kill it 
   * 
   * @param d   Detector 
   * @param r   Ring
   * @param s   Sector 
   * @param t   Strip 
   * @param m   Multiplicity
   * 
   * @return true if this was killed
   */
  Bool_t CheckDead(UShort_t d, Char_t r, UShort_t s, UShort_t t, Double_t& m);
  /** 
   * Re-calculate eta and correct the multiplicity accordingly 
   * 
   * @param d          Detector 
   * @param r          Ring
   * @param s          Sector 
   * @param t          Strip 
   * @param ip         Ip coordinates
   * @param mult       In/out multiplicity
   * @param eta        In/out eta 
   * @param cosTheta   On return, the cosine of theta or null
   */
  void RecalculateEta(UShort_t d, Char_t r, UShort_t s, UShort_t t, 
		      const TVector3& ip, Double_t& mult, Double_t& eta, 
		      Double_t& cosTheta);
  /** 
   * Correct for noise suppression
   * 
   * @param f            Factor to apply
   * @param c            Correction 
   * @param cosTheta     Cosine to theta 
   * @param mult         In/Out multiplity 
   *
   * @return true if signal is good, otherwise false 
   */
  Bool_t NoiseCorrect(Int_t f, Double_t c, Double_t cosTheta, Double_t& mult);

  Int_t    fRecoFactor;      // Noise factor used in Reco
  Double_t fMaxNoiseCorr;    // If noise corr above this, flag as dead 
  Bool_t   fRecalculateEta;  // Whether to recalc eta and angle cor (disp vtx)
  TBits    fXtraDead;        // List of extra dead channels
  Bool_t   fHasXtraDead;     // Whether we have xtra dead channels
  Bool_t   fInvalidIsEmpty;  // Consider kInvalidMult as zero 
  TH1*     fNoiseChange;     // Diagnostics
  TH1*     fEtaChange;       // Diagnostics
  TH1*     fDeadChange;      // Diagnostics
  
  ClassDef(AliFMDESDFixer,1);  // Fix FMD ESD object for issues 
};

#endif
// 
// EOF
// 
  
