//
// Task to analyse the AOD for for dN/deta in the forward regions 
//
#ifndef ALIFORWARDDNDETATASK_H
#define ALIFORWARDDNDETATASK_H
/**
 * @file   AliForwarddNdetaTask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:04:54 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_dndeta
 */
#include "AliBasedNdetaTask.h"
class TList;
class TH2D;
class TH1D;

/**
 * Tasks to determine @f$ dN/d\eta@f$ in the forward regions
 *
 * @ingroup pwg2_forward_tasks_dndeta
 * @ingroup pwg2_forward_dndeta
 */
class AliForwarddNdetaTask : public AliBasedNdetaTask
{
public:
  /** 
   * Constructor 
   * 
   */
  AliForwarddNdetaTask();
  /** 
   * Constructor
   * 
   * @param name    Name of task 
   */
  AliForwarddNdetaTask(const char* name);
  /**
   * Destructor
   * 
   */
  virtual ~AliForwarddNdetaTask() {}
  /** 
   * Called at each event 
   *
   * This is called once in the master 
   * 
   * @param option Not used 
   */
  virtual void UserExec(Option_t* option);
protected:
  /** 
   * Copy constructor 
   */
  AliForwarddNdetaTask(const AliForwarddNdetaTask& o);
  /** 
   * Assigmement operator
   * 
   * @return Reference to this
   */
  AliForwarddNdetaTask& operator=(const AliForwarddNdetaTask&) { return *this; }

  /** 
   * Retrieve the histogram 
   * 
   * @param aod AOD event 
   * @param mc  Whether to get the MC histogram or not
   * 
   * @return Retrieved histogram or null
   */
  TH2D* GetHistogram(const AliAODEvent* aod, Bool_t mc);
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

  struct CentralityBin : public AliBasedNdetaTask::CentralityBin 
  {
    /** 
     * Constructor 
     */
    CentralityBin() : AliBasedNdetaTask::CentralityBin(), fSumPrimary(0) {}
    /** 
     * Constructor 
     * 
     * @param name Name used for histograms (e.g., Forward)
     * @param low  Lower centrality cut in percent 
     * @param high Upper centrality cut in percent 
     */
    CentralityBin(const char* name, Short_t low, Short_t high)
      : AliBasedNdetaTask::CentralityBin(name, low, high), 
	fSumPrimary(0)
    {}
    /** 
     * Copy constructor 
     * 
     * @param other Object to copy from 
     */
    CentralityBin(const CentralityBin& other) 
      : AliBasedNdetaTask::CentralityBin(other), 
	fSumPrimary(other.fSumPrimary) 
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
     * @param vzMin       Minimum IP z coordinate
     * @param vzMax       Maximum IP z coordinate
     * @param primary     MC truth histogram
     */
    virtual void ProcessPrimary(const AliAODForwardMult* forward, 
				Int_t triggerMask,
				Double_t vzMin, Double_t vzMax, 
				const TH2D* primary);
    /** 
     * End of processing 
     * 
     * @param sums        List of sums
     * @param results     Output list of results
     * @param scheme      Normalisation scheme options
     * @param shapeCorr   Shape correction or nil
     * @param trigEff     Trigger efficiency 
     * @param symmetrice  Whether to symmetrice the results
     * @param rebin       Whether to rebin the results
     * @param rootProj    If true, use TH2::ProjectionX
     * @param corrEmpty   Whether to correct for empty bins
     * @param cutEdges    Whether to cut edges when rebinning
     * @param triggerMask Trigger mask 
     * @param color       Marker colour 
     * @param marker      Marker style 
     */
    virtual void End(TList*      sums, 
		     TList*      results,
		     UShort_t    scheme,
		     const TH1*  shapeCorr, 
		     Double_t    trigEff,
		     Bool_t      symmetrice,
		     Int_t       rebin, 
		     Bool_t      rootProj,
		     Bool_t      corrEmpty, 
		     Bool_t      cutEdges, 
		     Int_t       triggerMask,
		     Int_t       color, 
		     Int_t       marker);
  protected: 
    TH2D*           fSumPrimary;    //  Sum of primary histograms
    ClassDef(CentralityBin,1); // A centrality bin     
  };

  ClassDef(AliForwarddNdetaTask,1); // Determine multiplicity in forward region
};

#endif
//
// Local Variables:
//  mode: C++
// End:
//
