#ifndef ALIFMDCORRNOISEGAIN_H
#define ALIFMDCORRNOISEGAIN_H
#include <AliFMDFloatMap.h>

/**
 * Get the noise calibration.  That is, the ratio 
 *
 * @f[ 
 *   \frac{\sigma_{i}}{g_{i}k}
 * @f] 
 *
 * where @f$ k@f$ is a constant determined by the electronics of
 * units DAC/MIP, and @f$ \sigma_i, g_i@f$ are the noise and gain of
 * the @f$ i @f$ strip respectively.
 *
 * This correction is needed because some of the reconstructed data
 * (what which have an AliESDFMD class version less than or equal to
 * 3) used the wrong zero-suppression factor.  The zero suppression
 * factor used by the on-line electronics was 4, but due to a coding
 * error in the AliFMDRawReader a zero suppression factor of 1 was
 * assumed during the reconstruction.  This shifts the zero of the
 * energy loss distribution artificially towards the left (lover
 * valued signals).
 *
 * So let's assume the real zero-suppression factor is @f$ f@f$ while
 * the zero suppression factor @f$ f'@f$ assumed in the reconstruction
 * was (wrongly) lower.  The number of ADC counts @f$ c_i'@f$ used in
 * the reconstruction can be calculated from the reconstructed signal
 * @f$ m_i'@f$ by
 *
 * @f[
 *    c_i' = m_i \times g_i \times k / \cos\theta_i
 * @f] 
 *
 * where @f$\theta_i@f$ is the incident angle of the @f$ i@f$ strip. 
 * 
 * This number of counts used the wrong noise factor @f$ f'@f$ so to
 * correct to the on-line value, we need to do
 *
 * @f[ 
 *   c_i = c_i' - \lfloor f'\times n_i\rfloor + \lfloor f\times n_i\rfloor
 * @f] 
 * 
 * which gives the correct number of ADC counts over the pedestal. To
 * convert back to the scaled energy loss signal we then need to
 * calculate (noting that @f$ f,f'@f$ are integers)
 *
 * @f{eqnarray}{
 *    m_i &=& \frac{c_i \times \cos\theta_i}{g_i \times k}\\ 
 *    &=& \left(c_i' - \lfloor f'\times n_i\rfloor + 
 *           \lfloor f\times n_i\rfloor\right)\frac{\cos\theta}{g_i \times k}\\
 *    &=& \left(\frac{m_i'\times g_i\times k}{\cos\theta} -
 *            \lfloor f'\times n_i\rfloor + \lfloor f\times n_i\rfloor\right)
 *         \frac{\cos\theta}{g_i \times k}\\
 *    &=& m_i' + \frac{1}{g_i \times k}
 *         \left(\lfloor f\times n_i\rfloor-
 *             \lfloor f'\times n_i\rfloor\right)\cos\theta\\
 *    &=& m_i' + \frac{\lfloor n_i\rfloor}{g_i \times k}
 *        \left(f-f'\right)\cos\theta
 * @f}
 * 
 */
class AliFMDCorrNoiseGain : public TObject 
{
public:
  /**
   * Default constructor 
   */
  AliFMDCorrNoiseGain() : fValues(0) { fValues.Reset(-1); }
  /** 
   * Constructor from a float map 
   * 
   * @param map Construct from this map 
   */
  AliFMDCorrNoiseGain(const AliFMDFloatMap& map) : fValues(map) {}
  /** 
   * Get the noise value for a particular strip 
   * 
   * @param d  Detector
   * @param r  Ring 
   * @param s  Sector 
   * @param t  Strip 
   * 
   * @return Noise value for strip 
   */
  Float_t Get(UShort_t d, Char_t r, UShort_t s, UShort_t t) const 
  { 
    return fValues(d, r, s, t);
  }
  /** 
   * Set the value for a strip. 
   * 
   * @param d Detector 
   * @param r Ring 
   * @param s Sector
   * @param t Strip
   * @param x Value 
   */
  void Set(UShort_t d, Char_t r, UShort_t s, UShort_t t, Float_t x) 
  { 
    fValues(d, r, s, t) = x;
  }
  /** 
   * Get a reference to the noise map 
   * 
   * @return Noise map 
   */
  const AliFMDFloatMap& Values() { return fValues; }
protected: 
  AliFMDFloatMap fValues; // The noise-gain map 
  ClassDef(AliFMDCorrNoiseGain,1); // Clone of AliFMDCalibPedestal
};

#endif
// Local Variables:
//   mode: C++
// End:

 
