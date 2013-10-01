//
// Task to analyse the AOD for for dN/deta in the forward regions 
//
#ifndef ALIMCTRUTHDNDETATASK_H
#define ALIMCTRUTHDNDETATASK_H
/**
 * @file   AliMCTruthdNdetaTask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:04:54 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_dndeta
 */
#include "AliBasedNdetaTask.h"
class TList;
class TH2D;
class TH1D;

/**
 * Tasks to determine @f$ dN/d\eta@f$ in the forward regions
 *
 * @ingroup pwglf_forward_tasks_dndeta
 * @ingroup pwglf_forward_dndeta
 */
class AliMCTruthdNdetaTask : public AliBasedNdetaTask
{
public:
  /** 
   * Constructor 
   * 
   */
  AliMCTruthdNdetaTask();
  /** 
   * Constructor
   * 
   * @param name    Name of task 
   */
  AliMCTruthdNdetaTask(const char* name);
  /**
   * Destructor
   * 
   */
  virtual ~AliMCTruthdNdetaTask() {}
  /** 
   * Called at end of event processing.
   *
   * This is called once in the master 
   * 
   * @param option Not used 
   */
  virtual Bool_t Finalize();
protected:
  /** 
   * Copy constructor 
   *
   * @param o Object to copy from 
   */
  AliMCTruthdNdetaTask(const AliMCTruthdNdetaTask& o);
  /** 
   * Assigmement operator
   * 
   * @return Reference to this
   */
  AliMCTruthdNdetaTask& operator=(const AliMCTruthdNdetaTask&);

  /** 
   * Retrieve the histogram 
   * 
   * @param aod AOD event 
   * @param mc  Whether to get the MC histogram or not
   * 
   * @return Retrieved histogram or null
   */
  TH2D* GetHistogram(const AliAODEvent& aod, Bool_t mc);
  /** 
   * Get the marker style 
   * 
   * @return Marker style 
   */
  virtual Int_t GetMarker() const { return GetMarkerStyle(kStar); }
  /** 
   * Get the colour to use for markers (only pp - in PbPb we use a rainbow)
   * 
   * @return Marker colour 
   */
  virtual Int_t GetColor() const { return kGray+2; }
  /** 
   * Make a new centrality bin
   * 
   * @param name   Histogram names
   * @param l      Lower cut
   * @param h      Upper cut
   * 
   * @return Newly allocated object (of our type)
   */
  AliBasedNdetaTask::CentralityBin* 
  MakeCentralityBin(const char* name, Short_t l, Short_t h) const;

  /**
   * Class that holds data for a single centrality bin 
   * 
   */
  class CentralityBin : public AliBasedNdetaTask::CentralityBin 
  {
  public:
    /** 
     * Constructor 
     */
    CentralityBin() : AliBasedNdetaTask::CentralityBin(), fSumTruth(0) {}
    /** 
     * Constructor 
     * 
     * @param name Name used for histograms (e.g., Forward)
     * @param low  Lower centrality cut in percent 
     * @param high Upper centrality cut in percent 
     */
    CentralityBin(const char* name, Short_t low, Short_t high)
      : AliBasedNdetaTask::CentralityBin(name, low, high), 
	fSumTruth(0)
    {}
    /** 
     * Copy constructor 
     * 
     * @param other Object to copy from 
     */
    CentralityBin(const CentralityBin& other) 
      : AliBasedNdetaTask::CentralityBin(other), 
	fSumTruth(other.fSumTruth) 
    {}
    /** 
     * Destructor 
     */
    virtual ~CentralityBin() {}
    /** 
     * Assignement operator 
     * 
     * 
     * @return 
     */
    CentralityBin& operator=(const CentralityBin&) { return *this; }
    /** 
     * Process an event
     * 
     * @param forward     Forward data (for trigger, vertex, & centrality)
     * @param triggerMask Trigger mask 
     * @param isZero      True if this is a zero bin event 
     * @param vzMin       Minimum IP z coordinate
     * @param vzMax       Maximum IP z coordinate
     * @param data        Data histogram 
     * @param mc          MC histogram
     */
    virtual Bool_t ProcessEvent(const AliAODForwardMult* forward, 
				Int_t                    triggerMask,
				Bool_t                   isZero,
				Double_t                 vzMin, 
				Double_t                 vzMax, 
				const TH2D*              data, 
				const TH2D*              mc);
    /** 
     * End of processing 
     * 
     * @param sums        List of sums
     * @param results     Output list of results
     * @param scheme      Normalisation scheme options
     * @param shapeCorr   Shape correction or nil
     * @param trigEff     Trigger efficiency 
     * @param trigEff0    0-bin trigger efficiency 
     * @param symmetrice  Whether to symmetrice the results
     * @param rebin       Whether to rebin the results
     * @param rootProj    If true, use TH2::ProjectionX
     * @param corrEmpty   Whether to correct for empty bins
     * @param cutEdges    Whether to cut edges when rebinning
     * @param triggerMask Trigger mask 
     * @param color       Marker colour 
     * @param marker      Marker style 
     * @param mclist      List of MC results 
     * @param truthlist   List of MC truth results 
     */
    virtual void End(TList*      sums, 
		     TList*      results,
		     UShort_t    scheme,
		     const TH2F* shapeCorr, 
		     Double_t    trigEff,
		     Double_t    trigEff0,
		     Bool_t      symmetrice,
		     Int_t       rebin, 
		     Bool_t      rootProj,
		     Bool_t      corrEmpty, 
		     Bool_t      cutEdges, 
		     Int_t       triggerMask,
		     Int_t       marker,
		     Int_t       color,
		     TList*      mclist,
		     TList*      truthlist);
  protected: 
    TH2D*           fSumTruth;    //  Sum of primary histograms
    ClassDef(CentralityBin,2); // A centrality bin     
  };
  Bool_t fHasData; // whether we actually have data or not 
  ClassDef(AliMCTruthdNdetaTask,2); // Determine multiplicity in forward region
};

#endif
//
// Local Variables:
//  mode: C++
// End:
//
