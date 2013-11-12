#ifndef ALIFMDMCHITENERGYFITTER_H
#define ALIFMDMCHITENERGYFITTER_H
#include "AliFMDEnergyFitter.h"
#include "AliFMDFloatMap.h"
#include "TArrayF.h"
class AliMCAuxHandler;
class AliMCEvent;
class AliESDEvent;
class TNtuple;
class TClonesArray;

/**
 * Class to fit the simulated energy loss in the FMD
 * 
 * @ingroup pwglf_forward_mc
 * @ingroup pwglf_forward_eloss
 */
class AliFMDMCHitEnergyFitter : public AliFMDEnergyFitter
{
public:
  /** 
   * Constructor - do not use
   */
  AliFMDMCHitEnergyFitter();
  /** 
   * Constructor 
   * 
   * @param title    Title of object - not significant 
   * @param useTuple If true, also make an NTuple of hits on output 3
   */
  AliFMDMCHitEnergyFitter(const char* title, Bool_t useTuple=false);
  /**
   * Destructor
   */
  virtual ~AliFMDMCHitEnergyFitter();

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
   * @param handler  Handler of additional input
   * 
   * @return true
   */
  virtual Bool_t Event(const AliESDEvent&  esdInput,
		       const AliMCEvent&   mcInput,
		       AliMCAuxHandler&    handler);
  /** 
   * Accumulate signals from MC hits
   * 
   * @param mcInput  MC input event
   * @param hits     Clones array of hits 
   * 
   * @return 
   */
  virtual Bool_t AccumulateHits(const AliMCEvent&   mcInput,
				const TClonesArray& hits);
  /** 
   * Post-process accumulated signals 
   * 
   * @param esdInput ESD event 
   * 
   * @return true 
   */
  virtual Bool_t PostEvent(const AliESDEvent& esdInput);
  
  TNtuple* GetTuple() { return fTuple; }
protected:
  /**
   * copy constructor - not implemented
   */
  AliFMDMCHitEnergyFitter(const AliFMDMCHitEnergyFitter&);
  /**
   * Assignment operator - not implemented
   * 
   */
  AliFMDMCHitEnergyFitter& operator=(const AliFMDMCHitEnergyFitter&);

  /**
   * Container of ring histograms 
   * 
   */
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
    RingHistos(const RingHistos& o);
    /** 
     * Assignment operator  - not defined
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this 
     */
    RingHistos& operator=(const RingHistos& o);
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
    TH2* fKind;      // Particle kind
    ClassDef(RingHistos,1); // Cache of histograms per ring 
  };
  /** 
   * Create a container of histograms for a single ring
   * 
   * @param d Detector 
   * @param r Ring 
   * 
   * @return Newly allocated container 
   */
  AliFMDEnergyFitter::RingHistos* CreateRingHistos(UShort_t d, Char_t r) const;
  /** Cache of per-strip energy loss of primaries */
  AliFMDFloatMap   fSumPrimary;
  /** Cache of per-strip energy loss of secondaries */
  AliFMDFloatMap   fSumSecondary;
  /** Cache of current MC IP */
  TArrayF          fIp;
  /** Cache of number of tracks */
  Int_t            fNTrack;
  /** Cache of numbr of primaries */
  Int_t            fNPrimary;
  /** Output nTuple of per-hit information */
  TNtuple*         fTuple;

  ClassDef(AliFMDMCHitEnergyFitter,1);
};

#endif
// Local Variables:
//  mode: C++
// End:
