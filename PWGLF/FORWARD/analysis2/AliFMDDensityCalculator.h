// This class calculates the inclusive charged particle density
// in each for the 5 FMD rings. 
//
#ifndef ALIFMDDENSITYCALCULATOR_H
#define ALIFMDDENSITYCALCULATOR_H
/**
 * @file   AliFMDDensityCalculator.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:02:09 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_aod
 */
#include <TNamed.h>
#include <TList.h>
#include <TArrayI.h>
#include <TVector3.h>
#include "AliForwardUtil.h"
#include "AliFMDMultCuts.h"
#include "AliPoissonCalculator.h"
class AliESDFMD;
class TH2D;
class TH1D;
class TProfile;
class AliFMDCorrELossFit;

/** 
 * This class calculates the inclusive charged particle density
 * in each for the 5 FMD rings. 
 *
 * @par Input:
 *   - AliESDFMD object possibly corrected for sharing
 *
 * @par Output:
 *   - 5 RingHistos objects - each with a number of vertex dependent 
 *     2D histograms of the inclusive charge particle density 
 * 
 * @par Corrections used: 
 *   - AliFMDAnaCalibEnergyDistribution 
 *   - AliFMDDoubleHitCorrection 
 *   - AliFMDDeadCorrection 
 *
 * @ingroup pwglf_forward_algo
 * @ingroup pwglf_forward_aod
 */
class AliFMDDensityCalculator : public TNamed
{
public:
  /**
   * How to correct for the missing phi coverage at the corners of the
   * sensors 
   * 
   */
  enum { 
    /** No correction */
    kPhiNoCorrect,
    /** Correct the calculated number charged particles */
    kPhiCorrectNch,
    /** Correct the energy loss */
    kPhiCorrectELoss
  };
  /** 
   * Folder name 
   */
  static const char* fgkFolderName;
  /** 
   * Constructor 
   */
  AliFMDDensityCalculator();
  /** 
   * Constructor 
   * 
   * @param name Name of object
   */
  AliFMDDensityCalculator(const char* name);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDDensityCalculator(const AliFMDDensityCalculator& o);
  /** 
   * Destructor 
   */
  virtual ~AliFMDDensityCalculator();
  /** 
   * Assignement operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object
   */
  AliFMDDensityCalculator& operator=(const AliFMDDensityCalculator& o);
  /** 
   * Initialize this sub-algorithm
   * 
   * @param etaAxis Not used 
   */
  virtual void SetupForData(const TAxis& etaAxis);
  /** 
   * Do the calculations 
   * 
   * @param fmd      AliESDFMD object (possibly) corrected for sharing
   * @param hists    Histogram cache
   * @param lowFlux  Low flux flag. 
   * @param cent     Centrality 
   * @param ip       Coordinates of interaction point
   * 
   * @return true on successs 
   */
  virtual Bool_t Calculate(const AliESDFMD&        fmd, 
			   AliForwardUtil::Histos& hists, 
			   Bool_t   		   lowFlux, 
			   Double_t  		   cent=-1, 
			   const TVector3&         ip=TVector3(1024,1024,0));
  /** 
   * Scale the histograms to the total number of events 
   * 
   * @param dir     where to put the output
   * @param output  Output list
   * @param nEvents Number of events 
   */
  virtual void Terminate(const TList* dir, TList* output, Int_t nEvents);
  /** 
   * Output diagnostic histograms to directory 
   * 
   * @param dir List to write in
   */  
  virtual void CreateOutputObjects(TList* dir);
  /** 
   * Set the debug level.  The higher the value the more output 
   * 
   * @param dbg Debug level 
   */
  void SetDebug(Int_t dbg=1) { fDebug = dbg; }
  /** 
   * Maximum particle weight to use 
   * 
   * @param m 
   */
  void SetMaxParticles(UShort_t m) { fMaxParticles = m; }  
  /** 
   * Set whether to use poisson statistics to estimate the 
   * number of particles that has hit within a region.  If this is true, 
   * then the average charge particle density is given by 
   * @f[
   *  \lambda = -\log\left(\frac{N_e}{N_t}\right)
   * @f]
   * where $N_e$ is the number of strips within the region that has no
   * hits over threshold, and $N_t$ is the total number of strips in the 
   * region/ 
   * 
   * @param u Whether to use poisson statistics to estimate the 
   * number of particles that has hit within a region.
   */
  void SetUsePoisson(Bool_t u) { fUsePoisson = u; }
  /** 
   * In case of a displaced vertices recalculate eta and angle correction
   * 
   * @param use recalculate or not
   * 
   */
  void SetRecalculateEta(Bool_t use) { fRecalculateEta = use; }
  /** 
   * In case of a displaced vertices recalculate eta and angle correction
   * 
   * @param use recalculate or not
   * 
   */
  void SetRecalculatePhi(Bool_t use) { fRecalculatePhi = use; }
  /** 
   * Set whether to use the phi acceptance correction. 
   * 
   * How the phi acceptance is used depends on the value passed.  
   * - 0:  No phi acceptance 
   * - 1:  Phi acceptance correction done to estimate of particles 
   * - 2:  Phi acceptance correction done to energy deposited 
   *
   * @param u If >0, use the phi acceptance (default is false)
   */
   
  void SetUsePhiAcceptance(UShort_t u=kPhiCorrectNch) { fUsePhiAcceptance = u; }
  /** 
   * Set the luming factors used in the Poisson method
   * 
   * @param eta Must be 1 or larger 
   * @param phi Must be 1 or larger 
   */
  void SetLumping(Int_t eta, Int_t phi) { 
    fEtaLumping = (eta < 1 ? 1 : eta); 
    fPhiLumping = (phi < 1 ? 1 : phi); 
  }
  /** 
   * Get the multiplicity cut.  If the user has set fMultCut (via
   * SetMultCut) then that value is used.  If not, then the lower
   * value of the fit range for the enery loss fits is returned.
   * 
   * @param d      Detector 
   * @param r      Ring 
   * @param eta    Psuedo-rapidity
   * @param errors Factor in errors
   *
   * @return Lower cut on multiplicity
   */
  Double_t GetMultCut(UShort_t d, Char_t r, Double_t eta, 
		      Bool_t errors=true) const;
  /** 
   * Get the multiplicity cut.  If the user has set fMultCut (via
   * SetMultCut) then that value is used.  If not, then the lower
   * value of the fit range for the enery loss fits is returned.
   * 
   * @param d      Detector 
   * @param r      Ring 
   * @param ieta   Psuedo-rapidity bin
   * @param errors Factor in errors
   *
   * @return Lower cut on multiplicity
   */
  Double_t GetMultCut(UShort_t d, Char_t r, Int_t ieta, 
		      Bool_t errors=true) const;
  /** 
   * Print information 
   * 
   * @param option Print options 
   *   - max  Print max weights 
   */
  void Print(Option_t* option="") const;
  /** 
   * Get the cuts used 
   * 
   * @return Reference to cuts object
   */
  AliFMDMultCuts& GetCuts() { return fCuts; }
  /** 
   * Set the cuts to use 
   * 
   * @param c Cuts to use 
   */
  void SetCuts(const AliFMDMultCuts& c) { fCuts = c; }
protected:
  /** 
   * Find the max weight to use for FMD<i>dr</i> in eta bin @a iEta
   * 
   * @param cor   Correction
   * @param d     Detector 
   * @param r     Ring 
   * @param iEta  Eta bin 
   *
   * @return The maximum weight 
   */
  Int_t FindMaxWeight(const AliFMDCorrELossFit* cor,
		      UShort_t d, Char_t r, Int_t iEta) const;

  /** 
   * Find the max weights and cache them 
   * 
   * @param axis Default @f$\eta@f$ axis from parent task 
   */  
  void CacheMaxWeights(const TAxis& axis);
  /** 
   * Find the (cached) maximum weight for FMD<i>dr</i> in 
   * @f$\eta@f$ bin @a iEta
   * 
   * @param d     Detector
   * @param r     Ring
   * @param iEta  Eta bin
   * 
   * @return max weight or <= 0 in case of problems 
   */
  Int_t GetMaxWeight(UShort_t d, Char_t r, Int_t iEta) const;
  /** 
   * Find the (cached) maximum weight for FMD<i>dr</i> iat
   * @f$\eta@f$ 
   * 
   * @param d     Detector
   * @param r     Ring
   * @param eta   Eta bin
   * 
   * @return max weight or <= 0 in case of problems 
   */
  Int_t GetMaxWeight(UShort_t d, Char_t r, Float_t eta) const;

  /** 
   * Get the number of particles corresponding to the signal mult
   * 
   * @param mult     Signal
   * @param d        Detector
   * @param r        Ring 
   * @param eta      Pseudo-rapidity 
   * @param lowFlux  Low-flux flag 
   * 
   * @return The number of particles 
   */
  virtual Float_t NParticles(Float_t  mult, 
			     UShort_t d, 
			     Char_t   r, 
			     Float_t  eta, 
			     Bool_t   lowFlux) const;
  /** 
   * Get the inverse correction factor.  This consist of
   * 
   * - acceptance correction (corners of sensors) 
   * - double hit correction (for low-flux events) 
   * - dead strip correction 
   * 
   * @param d        Detector
   * @param r        Ring 
   * @param t        Strip 
   * @param eta      Pseudo-rapidity 
   * @param lowFlux  Low-flux flag 
   * 
   * @return the correction factor 
   */
  virtual Float_t Correction(UShort_t d, Char_t r, UShort_t t, 
			     Float_t eta, Bool_t lowFlux) const;
  /** 
   * Get the acceptance correction for strip @a t in an ring of type @a r
   * 
   * @param r  Ring type ('I' or 'O')
   * @param t  Strip number 
   * 
   * @return Inverse acceptance correction 
   */
  virtual Float_t AcceptanceCorrection(Char_t r, UShort_t t) const;
  /** 
   * Generate the acceptance corrections 
   * 
   * @param r Ring to generate for 
   * 
   * @return Newly allocated histogram of acceptance corrections
   */
  virtual TH1D*   GenerateAcceptanceCorrection(Char_t r) const;
  /** 
   * Internal data structure to keep track of the histograms
   */
  struct RingHistos : public AliForwardUtil::RingHistos
  { 
    /** 
     * Default CTOR
     */
    RingHistos();
    /** 
     * Constructor
     * 
     * @param d detector
     * @param r ring 
     */
    RingHistos(UShort_t d, Char_t r);
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    RingHistos(const RingHistos& o);
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this 
     */
    RingHistos& operator=(const RingHistos& o);
    /** 
     * Destructor 
     */
    ~RingHistos();
    /** 
     * Initialize the object 
     * 
     * @param eAxis 
     */
    void SetupForData(const TAxis& eAxis);
    /** 
     * Make output 
     * 
     * @param dir Where to put it 
     */
    void CreateOutputObjects(TList* dir);
    /** 
     * Scale the histograms to the total number of events 
     * 
     * @param dir     Where the output is 
     * @param nEvents Number of events 
     */
    void Terminate(TList* dir, Int_t nEvents);
    TList*    fList;
    TH2D*     fEvsN;           // Correlation of Eloss vs uncorrected Nch
    TH2D*     fEvsM;           // Correlation of Eloss vs corrected Nch
    TProfile* fEtaVsN;         // Average uncorrected Nch vs eta
    TProfile* fEtaVsM;         // Average corrected Nch vs eta
    TProfile* fCorr;           // Average correction vs eta
    TH2D*     fDensity;        // Distribution inclusive Nch
    TH2D*     fELossVsPoisson; // Correlation of energy loss vs Poisson N_ch
    TH1D*     fDiffELossPoisson;// Relative difference to Poisson
    AliPoissonCalculator fPoisson; // Calculate density using Poisson method
    TH1D*     fELoss;          // Energy loss as seen by this 
    TH1D*     fELossUsed;      // Energy loss in strips with signal 
    Double_t  fMultCut;        // If set, use this
    TH1D*     fTotal;          // Total number of strips per eta
    TH1D*     fGood;           // Number of good strips per eta
    TH2D*     fPhiAcc;         // Phi acceptance vs IpZ
    TH1D*     fPhiBefore;      // Phi before re-calce 
    TH1D*     fPhiAfter;       // Phi after re-calc
    ClassDef(RingHistos,9);
  };
  /** 
   * Get the ring histogram container 
   * 
   * @param d Detector
   * @param r Ring 
   * 
   * @return Ring histogram container 
   */
  RingHistos* GetRingHistos(UShort_t d, Char_t r) const;
  TList    fRingHistos;    //  List of histogram containers
  TH1D*    fSumOfWeights;  //  Histogram
  TH1D*    fWeightedSum;   //  Histogram
  TH1D*    fCorrections;   //  Histogram
  UShort_t fMaxParticles;  //  Maximum particle weight to use 
  Bool_t   fUsePoisson;    //  If true, then use poisson statistics 
  UShort_t fUsePhiAcceptance; // Whether to correct for corners 
  TH1D*    fAccI;          //  Acceptance correction for inner rings
  TH1D*    fAccO;          //  Acceptance correction for outer rings
  TArrayI  fFMD1iMax;      //  Array of max weights 
  TArrayI  fFMD2iMax;      //  Array of max weights 
  TArrayI  fFMD2oMax;      //  Array of max weights 
  TArrayI  fFMD3iMax;      //  Array of max weights 
  TArrayI  fFMD3oMax;      //  Array of max weights 
  TH2D*    fMaxWeights;    //  Histogram of max weights
  TH2D*    fLowCuts;       //  Histogram of low cuts
  Int_t    fEtaLumping;    //  How to lump eta bins for Poisson 
  Int_t    fPhiLumping;    //  How to lump phi bins for Poisson 
  Int_t    fDebug;         //  Debug level 
  AliFMDMultCuts fCuts;    // Cuts
  Bool_t fRecalculateEta;  // Whether to recalc eta and angle correction (disp vtx)
  Bool_t fRecalculatePhi;  // Whether to correct for (X,Y) offset

  ClassDef(AliFMDDensityCalculator,9); // Calculate Nch density 
};

#endif
// Local Variables:
//   mode: C++
// End:

