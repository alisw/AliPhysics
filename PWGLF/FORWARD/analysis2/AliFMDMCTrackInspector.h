#ifndef ALIFMDMCTRACKINSPECTOR_H
#define ALIFMDMCTRACKINSPECTOR_H
#include "AliFMDEnergyFitter.h"
#include "AliFMDMCTrackELoss.h"
#include <TArrayF.h>
class AliMCEvent;
class AliESDEvent;

/**
 * Class to fit the simulated energy loss in the FMD
 * 
 * @ingroup pwglf_forward_mc
 * @ingroup pwglf_forward_eloss
 */
class AliFMDMCTrackInspector : public AliFMDEnergyFitter
{
public:
  /** 
   * Constructor - do not use
   */
  AliFMDMCTrackInspector();
  /** 
   * Constructor 
   * 
   * @param title    Title of object - not significant 
   */
  AliFMDMCTrackInspector(const char* title);
  /**
   * Destructor
   */
  virtual ~AliFMDMCTrackInspector();

  /** 
   * Define the output histograms.  These are put in a sub list of the
   * passed list.   The histograms are merged before the parent task calls 
   * AliAnalysisTaskSE::Terminate 
   * 
   * @param dir Directory to add to 
   */
  virtual void CreateOutputObjects(TList* dir);
  /** 
   * Set-up before processing an event 
   * 
   * @param mcInput MC input 
   * 
   * @return true
   */
  virtual Bool_t PreEvent(const AliMCEvent& mcInput);
  /** 
   * Process a single event 
   * 
   * @param esdInput ESD input 
   * @param mcInput  MC input
   * @param cent     Event centrality 
   * 
   * @return true
   */
  virtual Bool_t Event(const AliESDEvent&  esdInput,
		       const AliMCEvent&   mcInput,
		       Double_t            cent=-1);
  /** 
   * Post-process accumulated signals 
   * 
   * @return true 
   */
  virtual Bool_t PostEvent();

  /** 
   * Get reference to the tracker of energy loss 
   * 
   * @return Reference to tracker of energy loss 
   */  
  AliFMDMCTrackELoss& GetTracker() { return fTracker; }
  /** 
   * Get constant reference to the tracker of energy loss 
   * 
   * @return Constant reference to tracker of energy loss 
   */  
  const AliFMDMCTrackELoss& GetTracker() const { return fTracker; }
protected:
  /**
   * copy constructor - not implemented
   */
  AliFMDMCTrackInspector(const AliFMDMCTrackInspector&);
  /**
   * Assignment operator - not implemented
   * 
   */
  AliFMDMCTrackInspector& operator=(const AliFMDMCTrackInspector&);

  /**
   * Container of ring histograms 
   * 
   */
public:
  struct RingHistos : public AliFMDEnergyFitter::RingHistos
  {
    /** 
     * Default CTOR - do not use 
     */
    RingHistos();
    /** 
     * User CTOR 
     * 
     * @param d Detector number 
     * @param r Ring identifier 
     */
    RingHistos(UShort_t d, Char_t r);
    /** 
     * DTOR
     */
    ~RingHistos() {}
    /** 
     * Copy constructor - not defined
     * 
     * @param o Object to copy from 
     */
    RingHistos(const RingHistos& o){;}
    /** 
     * Assignment operator  - not defined
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this 
     */
    RingHistos& operator=(const RingHistos& o){return *this;}
    /** 
     * Create a bin array of increasing bins. This overload uses the
     * service AliFMDEncodedEdx::Spec::FillBinArray. 
     * 
     * @param nBins Number of bins - ignored
     * @param low   Low cut - ignored
     * @param high  High cut - ignored
     * 
     * @return Array of bin boundaries 
     */
    TArrayD MakeIncreasingAxis(Int_t nBins, 
			       Double_t low,
			       Double_t high) const;
    /** 
     * Initialise object 
     * 
     * @param eAxis      Eta axis
     * @param cAxis      Centrality axis 
     * @param maxDE      Max energy loss to consider 
     * @param nDEbins    Number of bins 
     * @param useIncrBin Whether to use an increasing bin size 
     */
    virtual void SetupForData(const TAxis& eAxis, 
			      const TAxis& cAxis,
			      Double_t     maxDE=10, 
			      Int_t        nDEbins=300, 
			      Bool_t       useIncrBin=true);
    /** 
     * Fill in observation 
     * 
     * @param flag 0 - fill all, 1 - primary, 2 - secondary
     * @param eta  Eta of particle observations
     * @param mult Scaled energy loss 
     */
    virtual void FillMC(UShort_t flag, Double_t eta, Double_t mult);
    /** 
     * Do scaling of histogram before fitting.  This can be
     * overwritten to do some smoothing or the like. By default, this
     * simply scales to the bin width.
     * 
     * @param dist Histogram to scale. 
     */     
    virtual void  Scale(TH1* dist) const;
    /** 
     * Fit the final distributions - called via Terminate 
     * 
     * @param dir           Containing directory
     * @param lowCut        Lower cut on @f$\Delta/\Delta_{mip}@f$ 
     * @param nParticles    Max. number of particle peaks to fit
     * @param minEntries    Least number of entries required before fitting
     * @param minusBins     Number of bins below the 1st peak we start fitting
     * @param relErrorCut   Largest relative error on paramters
     * @param chi2nuCut     Largest value of the @f$\chi^2/\nu@f$ 
     * @param minWeight     Least weight to consider 
     * @param regCut        When to regalurize 
     * @param residuals     How to do residuals - if at all 
     * 
     * @return List of histograms of parameters 
     */
    TObjArray* Fit(TList*           dir, 
		   Double_t         lowCut, 
		   UShort_t         nParticles,
		   UShort_t         minEntries,
		   UShort_t         minusBins, 
		   Double_t         relErrorCut, 
		   Double_t         chi2nuCut,
		   Double_t         minWeight,
		   Double_t         regCut,
		   EResidualMethod  residuals) const;
    
    TH2* fPrimary;   // @f$\Delta@f$ vs @f$\eta@f$ for primaries
    TH2* fSecondary; // @f$\Delta@f$ vs @f$\eta@f$ for second.
  public:
    TH2* fBetaGammadEdx;
    TH2* fBetaGammaEta;
    TH2* fDedxEta;
    ClassDef(RingHistos,1); // Cache of histograms per ring 
  };
protected:
  /** 
   * Create a container of histograms for a single ring
   * 
   * @param d Detector 
   * @param r Ring 
   * 
   * @return Newly allocated container 
   */
  AliFMDEnergyFitter::RingHistos* CreateRingHistos(UShort_t d, Char_t r) const;
  /** Our 'tracker' */
  AliFMDMCTrackELoss fTracker;
  /** Cache of current MC IP */
  TArrayF          fIp;
  /** Cache of number of tracks */
  Int_t            fNTrack;
  /** Cache of numbr of primaries */
  Int_t            fNPrimary;

  ClassDef(AliFMDMCTrackInspector,1);
};

#endif
// Local Variables:
//  mode: C++
// End:
