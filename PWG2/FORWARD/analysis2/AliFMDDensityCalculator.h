#ifndef ALIROOT_PWG2_FORWARD_ANALYSIS2_ALIFMDDENSITYCALCULATOR_H
#define ALIROOT_PWG2_FORWARD_ANALYSIS2_ALIFMDDENSITYCALCULATOR_H
#include <TNamed.h>
#include <TList.h>
#include "AliForwardUtil.h"
class AliESDFMD;
class TH2D;

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
 * @ingroup pwg2_forward_analysis 
 */
class AliFMDDensityCalculator : public TNamed
{
public:
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
  AliFMDDensityCalculator& operator=(const AliFMDDensityCalculator&);
  /** 
   * Do the calculations 
   * 
   * @param fmd      AliESDFMD object (possibly) corrected for sharing
   * @param hists    Histogram cache
   * @param vtxBin   Vertex bin 
   * @param lowFlux  Low flux flag. 
   * 
   * @return true on successs 
   */
  virtual Bool_t Calculate(const AliESDFMD& fmd, 
			   AliForwardUtil::Histos& hists, 
			   Int_t vtxBin, Bool_t lowFlux);
  /** 
   * Scale the histograms to the total number of events 
   * 
   * @param nEvents Number of events 
   */
  void ScaleHistograms(Int_t nEvents);
  /** 
   * Output diagnostic histograms to directory 
   * 
   * @param dir List to write in
   */  
  void Output(TList* dir);
protected:
  /** 
   * Get the number of particles corresponding to the signal mult
   * 
   * @param mult     Signal
   * @param d        Detector
   * @param r        Ring 
   * @param s        Sector 
   * @param t        Strip (not used)
   * @param v        Vertex bin 
   * @param eta      Pseudo-rapidity 
   * @param lowFlux  Low-flux flag 
   * 
   * @return The number of particles 
   */
  virtual Float_t NParticles(Float_t mult, 
			     UShort_t d, Char_t r, UShort_t s, UShort_t t, 
			     Int_t v, Float_t eta, Bool_t lowFlux) const;
  /** 
   * Get the inverse correction factor.  This consist of
   * 
   * - acceptance correction (corners of sensors) 
   * - double hit correction (for low-flux events) 
   * - dead strip correction 
   * 
   * @param d        Detector
   * @param r        Ring 
   * @param s        Sector 
   * @param t        Strip (not used)
   * @param v        Vertex bin 
   * @param eta      Pseudo-rapidity 
   * @param lowFlux  Low-flux flag 
   * 
   * @return 
   */
  virtual Float_t Correction(UShort_t d, Char_t r, UShort_t s, UShort_t t, 
			     Int_t v, Float_t eta, Bool_t lowFlux) const;
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
   * Internal data structure to keep track of the histograms
   */
  struct RingHistos : public TObject 
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
    void Output(TList* dir);
    UShort_t  fDet;          // Detector
    Char_t    fRing;         // Ring
    TH2D*     fEvsN;         // Correlation of Eloss vs uncorrected Nch
    TH2D*     fEvsM;         // Correlation of Eloss vs corrected Nch
    TH2D*     fDensity;      // Distribution inclusive Nch
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

  TList    fRingHistos;    // List of histogram containers
  Double_t fMultCut;       // Low cut on scaled energy loss

  ClassDef(AliFMDDensityCalculator,1); // Calculate Nch density 
};

#endif
// Local Variables:
//   mode: C++
// End:

