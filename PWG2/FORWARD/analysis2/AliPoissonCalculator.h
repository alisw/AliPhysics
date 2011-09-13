#ifndef ALIPOISSONCALCULATOR_H
#define ALIPOISSONCALCULATOR_H
#include <TNamed.h>
#include <TList.h>
class TH2D;
class TH1D;
class TBrowser;
//class TList;

/** 
 * A class to calculate the multiplicity in @f$(\eta,\varphi)@f$ bins
 * using Poisson statistics. 
 *
 * The input is assumed to be binned in @f$(\eta,\varphi)@f$ as
 * described by the 2D histogram passwd to the Reset member function.  
 *
 * The data is grouped in to regions as defined by the parameters
 * fEtaLumping and fPhiLumping.  The total number of cells and number
 * of empty cells is then calculate in each region.  The mean
 * multiplicity over the region is then determined as 
 *
 * @f[
 * \langle m\rangle = -\log\left(\frac{e}{t}\right)
 * @f]
 * where @f$ e@f$ is the number of empty cells and @f$t@f$ is the
 * total number of cells in the region.  A correction for counting
 * statistics, is then applied 
 * @f{eqnarray*}{
 *    c &=& \frac{1}{1 - \exp{-\langle m\rangle}}\\ &=& 
 *          \frac{1}{1 - \frac{e}{t}}
 * @f}
 * and the final number in each cell is then 
 * @f$h_i c \langle m\rangle@f$ 
 * where @f$h_i@f$ is the number of hits in the cell @f$i@f$ 
 * 
 */
class AliPoissonCalculator : public TNamed
{
public:
  /** 
   * Constructor 
   */
  AliPoissonCalculator(); 
  /** 
   * Constructor 
   * 
   */
  AliPoissonCalculator(const char*/*, UShort_t d, Char_t r*/);
  /** 
   * Copy constructor
   * 
   * @param o Object to copy from
   */
  AliPoissonCalculator(const AliPoissonCalculator& o);

  /** 
   * Destructor 
   */
  virtual ~AliPoissonCalculator();
  /** 
   * Assignment operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object
   */  
  AliPoissonCalculator& operator=(const AliPoissonCalculator& o);
  /** 
   * Set the number of eta bins to group into a region
   * 
   * @param n Number of eta bins per region
   */  
  void SetEtaLumping(UShort_t n) { fEtaLumping = n; } //*MENU*
  /** 
   * Set the actual object
   * 
   * @param v vtxbin
   */  
  void SetObject(UShort_t d, Char_t r, UShort_t v, Double_t cent); //*MENU*
  /** 
   * Set the number of phi bins to group into a region
   * 
   * @param n Number of phi bins per region
   */  
  void SetPhiLumping(UShort_t n) { fPhiLumping = n; } //*MENU*
  /** 
   * Intialize this object 
   * 
   * @param etaLumping If larger than 0, set the eta lumping to this
   * @param phiLumping If larger than 0, set the phi lumping to this
   */
  void Init(UShort_t d=-1, Char_t r='I',Int_t etaLumping=-1, Int_t phiLumping=-1);
  /** 
   * Output stuff to the passed list
   * 
   * @param d List to add output histograms to 
   */
  void Output(TList* d);
  /** 
   * Reset the cache histogram 
   * 
   * @param base Basic histogram 
   */
  void Reset(const TH2D* base);
  /** 
   * Fill in an observation 
   * 
   * @param eta     Eta value 
   * @param phi     Phi value
   * @param hit     True if hit 
   * @param weight  Weight if this 
   */
  void Fill(UShort_t strip, UShort_t sec, Bool_t hit, Double_t weight=1);
  /** 
   * Take the average occupancy over many events or not
   * 
   * @param use use or not
   */
  void SetUseAverageOverEvents(Bool_t use)  {fRunningAverage = use;}
  /** 
   * Calculate result and store in @a output
   * 
   * @return The result histogram (fBase overwritten)
   */
  TH2D* Result();
  /** 
   * @return Always true 
   */
  Bool_t IsFolder() const { return kTRUE; }
  /** 
   * Print information 
   * 
   * @param option Not used
   */
  void Print(const Option_t* option="") const;
  /** 
   * Browse this object
   * 
   * @param b Object to browse 
   */
  void Browse(TBrowser* b);

  TH2D* GetEmptyVsTotal() const { return fEmptyVsTotal; }
  TH1D* GetMean() const { return fMean; }
  TH1D* GetOccupancy() const { return fOcc; }
  TH2D* GetCorrection() const { return fCorr; }
protected:
  /** 
   * Clean up allocated space 
   * 
   */
  void CleanUp();
  /** 
   * Calculate the mean 
   *
   * This is based on the fact that for a Poisson
   * @f[ 
   *   P(n;\lambda) = \frac{-\lambda^n e^{-\lambda}}{n!}
   * @f]
   * we have the probability for 0 observation 
   * @f[
   *   P(0;\lambda) = e^{-\lambda} = \frac{N_{empty}}{N_{total}}
   * @f] 
   * and so we get that the mean is the defined region is 
   * @f[
   *   \lambda = -\log\left(\frac{N_{empty}}{N_{total}}\right)
   * @f]
   * 
   * Note the boundary conditions
   * - @f$N_{total}=0 \rightarrow\lambda=0@f$
   * - @f$N_{empty}<\epsilon\rightarrow N_{empty} = \epsilon@f$ 
   * 
   * @param empty Number of empty bins 
   * @param total Total number of bins 
   * 
   * @return The mean in the defined region 
   */
  Double_t CalculateMean(Double_t empty, Double_t total) const;
  /** 
   * The mean @f$\lambda@f$ calculated above is not the full story.
   * In addition it needs to be corrected using the expression 
   * @f[
   *   \frac{1}{1-e^{\lambda}} =
   *     \frac{1}{1-\frac{N_{empty}}{N_{total}}}
   * @f]
   * 
   * Note the boundary conditions
   * - @f$N_{total}=0 \rightarrow\lambda=0@f$
   * - @f$|N_{total}-N_{empty}|<\epsilon\rightarrow N_{empty} =
   * N_{total}-\epsilon@f$  
   *
   * @param empty Number of empty bins 
   * @param total Total number of bins 
   * 
   * @return The correction to the mean. 
   */
  Double_t CalculateCorrection(Double_t empty, Double_t total) const;
  UShort_t fEtaLumping;   // Grouping of eta bins 
  UShort_t fPhiLumping;   // Grouping of phi bins 
  TH2D*    fTotal;        // Total number of strips in a region
  TH2D*    fEmpty;        // Total number of strips in a region
  TH2D*    fBasic;        // Total number basic hits in a region
  TH2D*    fEmptyVsTotal; // Empty versus total cells 
  TH1D*    fMean;         // Mean calculated by poisson method 
  TH1D*    fOcc;          // Histogram of occupancies 
  TH2D*    fCorr;         // Correction as a function of mean 
  TList    fTotalList;    // List of total hists
  TList    fEmptyList;    // List of empty hists
  Bool_t   fRunningAverage; // Whether to take average per event or not
  ClassDef(AliPoissonCalculator,1) // Calculate N_ch using Poisson
};

#endif
// Local Variables:
//   mode: C++
// End:

