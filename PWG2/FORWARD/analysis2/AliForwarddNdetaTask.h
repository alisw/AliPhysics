//
// Task to analyse the AOD for for dN/deta in the forward regions 
//
#ifndef ALIFORWARDDNDETATASK_H
#define ALIFORWARDDNDETATASK_H
#include "AliBasedNdetaTask.h"
class TList;
class TH2D;
class TH1D;

/**
 * Task to determine the 
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
   * @param maxVtx  Set @f$v_z@f$ range
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
     * @param shapeCorr   Shape correction or nil
     * @param trigEff     Trigger efficiency 
     * @param symmetrice  Whether to symmetrice the results
     * @param rebin       Whether to rebin the results
     * @param corrEmpty   Whether to correct for empty bins
     * @param cutEdges    Whether to cut edges when rebinning
     * @param vzMin       Minimum IP z coordinate
     * @param vzMax 	  Maximum IP z coordinate
     * @param triggerMask Trigger mask 
     */
    virtual void End(TList*      sums, 
		     TList*      results,
		     const TH1*  shapeCorr, 
		     Double_t    trigEff,
		     Bool_t      symmetrice,
		     Int_t       rebin, 
		     Bool_t      corrEmpty, 
		     Bool_t      cutEdges, 
		     Double_t    vzMin, 
		     Double_t    vzMax, 
		     Int_t       triggerMask);
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
