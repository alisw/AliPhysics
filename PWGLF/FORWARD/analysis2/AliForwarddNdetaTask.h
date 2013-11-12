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
 * @ingroup pwglf_forward_dndeta
 */
#include "AliBasedNdetaTask.h"
class TList;
class TH2D;
class TH1D;

/**
 * Tasks to determine @f$ dN/d\eta@f$ in the forward regions
 *
 * @image html alice-int-2012-040-performance_spdfmdvzero.png "dN/deta in PbPb"
 * 
 * @ingroup pwglf_forward_tasks_dndeta
 * @ingroup pwglf_forward_dndeta
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
protected:
  /** 
   * Copy constructor 
   *
   * @param o object to copy from 
   */
  AliForwarddNdetaTask(const AliForwarddNdetaTask& o);
  /** 
   * Assigmement operator
   * 
   * @return Reference to this
   */
  AliForwarddNdetaTask& operator=(const AliForwarddNdetaTask&);

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
   * Get the colour to use for markers (only pp - in PbPb we use a rainbow)
   * 
   * @return Marker colour 
   */
  virtual Int_t GetColor() const { return kRed+2; }
  /** 
   * Massage data histograms for certain vertices in the satellite analysis 
   * 
   * @param vtx 
   * @param data 
   * @param mcData 
   */
  virtual void CheckEventData(Double_t vtx, 
			      TH2*     data, 
			      TH2*     mcData);
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
   * A structure holding the per-centrality bin information 
   * 
   */
  class CentralityBin : public AliBasedNdetaTask::CentralityBin 
  {
  public:
    /** 
     * Constructor 
     */
    CentralityBin() : AliBasedNdetaTask::CentralityBin() {}
    /** 
     * Constructor 
     * 
     * @param name Name used for histograms (e.g., Forward)
     * @param low  Lower centrality cut in percent 
     * @param high Upper centrality cut in percent 
     */
    CentralityBin(const char* name, Short_t low, Short_t high)
      : AliBasedNdetaTask::CentralityBin(name, low, high)
    {}
    /** 
     * Copy constructor 
     * 
     * @param other Object to copy from 
     */
    CentralityBin(const CentralityBin& other);
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
    CentralityBin& operator=(const CentralityBin&);
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
    ClassDef(CentralityBin,3); // A centrality bin     
  };

  ClassDef(AliForwarddNdetaTask,3); // Determine multiplicity in forward region
};

#endif
//
// Local Variables:
//  mode: C++
// End:
//
