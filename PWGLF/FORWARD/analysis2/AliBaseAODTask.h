#ifndef ALIBASEAODTASK_H
#define ALIBASEAODTASK_H
#include <AliAnalysisTaskSE.h>
#include <TAxis.h>
class AliAODEvent;
class AliAODForwardMult;
class AliAODCentralMult;
class AliAODMultEventClass;
class TList;

/**
 * Base class for reading in AOD stuff 
 * 
 */
class AliBaseAODTask : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor (for I/O - do not use)
   */
  AliBaseAODTask(); 
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   * @param title Class name used in configuration script 
   */
  AliBaseAODTask(const char* name,
		 const char* title);
  /** 
   * Destructor
   */
  virtual ~AliBaseAODTask() {} 
  /** 
   * Configure this task via a macro 
   * 
   * @param macro Macro to configure va 
   * 
   * @return true on success, false otherwise
   */
  virtual Bool_t Configure(const char* macro="-default-");
  /** 
   * @{ 
   * @name Set parameters 
   */
  /** 
   * Set the vertex range to use 
   * 
   * @param min Minimum (in centermeter)
   * @param max Maximum (in centermeter)
   */  
  void SetIpZRange(Double_t min, Double_t max) { fMinIpZ=min; fMaxIpZ=max; }
  void SetIpZMin(Double_t min) { fMinIpZ = min; }
  void SetIpZMax(Double_t max) { fMaxIpZ = max; }
  /** 
   * Set the trigger maskl 
   * 
   * @param mask Trigger mask
   */
  void SetTriggerMask(UInt_t mask);
  /** 
   * Set the trigger mask 
   * 
   * @param mask trigger mask 
   */
  void SetTriggerMask(const char* mask);
  /** 
   * Set mask of events to filter out 
   * 
   * @param mask The fitler mask 
   */
  void SetFilterMask(UInt_t mask);
  void SetFilterMask(const char* mask);
  /** 
   * Set the centrality bins to use. 
   * 
   * @code 
   *   UShort_t bins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
   *   task->SetCentralityBins(11, bins);
   * @endcode 
   * 
   * @param n     Number of bins (elements in @a bins minus 1)
   * @param bins  Bin limits 
   */
  void SetCentralityAxis(UShort_t n, Short_t* bins);
  /** 
   * Set the centrality bins to use. 
   * 
   * @code 
   *   UShort_t bins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
   *   task->SetCentralityBins(11, bins);
   * @endcode 
   * 
   * @param n     Number of bins (elements in @a bins minus 1)
   * @param bins  Bin limits 
   */
  void SetCentralityAxis(UShort_t n, Double_t* bins);
  /** 
   * Define a single centrality bin from @a low to @a high 
   * 
   * @param low  Lower bound 
   * @param high Upper bound
   */
  void SetCentralityAxis(Short_t low, Short_t high);
  /** 
   * Set the centrality axis to use based on a string.  The bin edges
   * are separated by colons.
   * 
   * @param bins String of bin edges
   */
  void SetCentralityAxis(const char* bins);
  /* @} */
  /** 
   * @{ 
   * @name Interface member functions 
   */
  /** 
   * Add this task to the manager and connect the outputs.  If @a
   * sumFile is null or the empty string, then the sum container is
   * stored in the default output file of the manager.  If @a resFile
   * is null or the empty string, then it is set to @a resFile if
   * defined, otherwise to the default output file of the manager.
   * 
   * @param sumFile Output file for sums
   * @param resFile Output file for sums
   * 
   * @return true on success 
   */
  virtual Bool_t Connect(const char* sumFile=0, const char* resFile=0);
  /** 
   * Book output objects. Derived class should define this to book
   * output objects on the processing output list @c fList before the
   * actual event processing.  This is called on the master and on
   * each slave.
   * 
   * If this member function returns false, the execution is stopped
   * with a fatal signal.
   *
   * @return true on success. 
   */
  virtual Bool_t Book() = 0;
  /** 
   * Called after reading in the first event. Here we can setup stuff
   * depending on the conditions we're running under.
   * 
   * @return true on success.  If this returns false, then we turn the
   * task into a zombie and we do no more processing.
   */
  virtual Bool_t PreData() { return true; }
  /** 
   * Called before processing a single event - should not do anything
   * but clear data, etc.
   * 
   * @return true on success
   */
  virtual Bool_t PreEvent() { return true; } 
  /** 
   * Process a single event
   * 
   * @param aod Input event 
   * 
   * @return true on success 
   */
  virtual Bool_t Event(AliAODEvent& aod) = 0;
  /** 
   * Called after processing a single event - should not do anything
   * but clear data, etc.
   * 
   * @return true on success
   */
  virtual Bool_t PostEvent() { return true; } 
  /** 
   * Do the final analysis on the merged output. 
   * 
   * @return true on success
   */
  virtual Bool_t Finalize() = 0;
  /* @} */
     
  /** 
   * @{ 
   * @name Utilities 
   */
  /** 
   * Print information 
   * 
   * @param option Not used
   */
  void Print(Option_t* option="") const;
  /* @} */  
protected:
  /** 
   * Copyt constructor - not defined
   *
   * @param o Object to copy from 
   */
  AliBaseAODTask(const AliBaseAODTask& o); 
  /** 
   * Assignment operatoer - not defined
   *
   * @param o Object to assign from 
   *
   * @return reference to this object
   */
  AliBaseAODTask& operator=(const AliBaseAODTask& o); 
  /** @{ 
   *  @name Task interface 
   */
  /** 
   * Initialise on master - does nothing
   * 
   */
  virtual void   Init() {}
  /** 
   * Create output objects.  
   *
   * This is called once per slave process 
   */
  virtual void UserCreateOutputObjects();
  /** 
   * Process a single event 
   * 
   * @param option Not used
   */
  virtual void UserExec(Option_t* option);
  /** 
   * Called at end of event processing.
   *
   * This is called once in the master 
   * 
   * @param option Not used 
   */
  virtual void Terminate(Option_t* option);
  /* @} */
  /** 
   * Get the forward object from the AOD 
   * 
   * @param aod AOD event 
   * @param mc   If true, for MC 
   * @param verb If truem be verbose
   * 
   * @return Forward object or null
   */
  AliAODForwardMult* GetForward(const AliAODEvent& aod, 
				Bool_t mc=false, 
				Bool_t verb=true);
  /** 
   * Get the Mult event class object from the AOD 
   * 
   * @param aod AOD event 
   * @param verb If truem be verbose
   * 
   * @return Forward object or null
   */
  AliAODMultEventClass* GetMultClass(const AliAODEvent& aod, 
				      Bool_t verb=true);
  /** 
   * Get the central object from the AOD 
   * 
   * @param aod  AOD event 
   * @param mc   If true, for MC 
   * @param verb If truem be verbose
   * 
   * @return Central object or null
   */
  AliAODCentralMult* GetCentral(const AliAODEvent& aod, 
				Bool_t mc=false, 
				Bool_t verb=true);
  /** 
   * Get the histogram of primary particles 
   * 
   * @param aod AOD event 
   * 
   * @return Pointer to primary particles, or null
   */
  TH2D* GetPrimary(const AliAODEvent& aod);
  /** 
   * Store information about the job on output 
   * 
   * @param forward Forward object
   */
  virtual void StoreInformation(AliAODForwardMult& forward);
  /** 
   * Check if the event corresponds to the selected trigger(s),
   * vertex, and centrality.  Derived classes can overload this to
   * enable event processing - even if the event is not within cuts.
   * 
   * @param forward Forward object
   * 
   * @return true if the event is within the cuts. 
   */
  virtual Bool_t CheckEvent(const AliAODForwardMult& forward);
  /** 
   * Check if we have centrality bins defined
   * 
   * @return true if we have one or more centrality bins 
   */
  Bool_t HasCentrality() const 
  { 
    return (fCentAxis.GetNbins() >= 1 && 
	    fCentAxis.GetXbins() && 
	    fCentAxis.GetXbins()->GetArray()); 
  }
  /** 
   * Get the centrality.  The trigger mask of the forward object is
   * not modified
   * 
   * @param event    Our event 
   * @param forward  Our FMD event 
   * @param qual     On return, the quality flag
   *  
   * @return The centrality percentage 
   */
  virtual Double_t GetCentrality(AliAODEvent&       event,
				 AliAODForwardMult* forward,
				 Int_t&             qual);
  /**
   * Get the centrality.  If the quality is bad, set the corresponding
   * bit on the forward object.
   *
   * @param event    Our event 
   * @param forward  Our FMD event 
   */
  virtual Double_t GetCentrality(AliAODEvent&       event,
				 AliAODForwardMult* forward);
  /** 
   * Get the Z coordinate of the interaction point 
   * 
   * @param event    Our event
   * @param forward  Our FMD event
   * 
   * @return The z coordinate of the interaction point 
   */
  virtual Double_t GetIpZ(AliAODEvent& event,
			  AliAODForwardMult* forward);
  /** 
   * Get the IPs (x,y) coordinates 
   * 
   * @param aod Input event 
   * @param x   On return, the X coordinate 
   * @param y   On return, the Y coordinate 
   * 
   * @return true on success, false otherwise 
   */
  virtual Bool_t GetIpXY(AliAODEvent& aod, Double_t& x, Double_t& y);
  /** 
   * Get the name of the default configuration script to use.
   * Sub-classes can override this to give another default
   * configuration script.  Note, it should problably only return the
   * base name (not full path) of the script.
   * 
   * @return Name of the configuration script to use. 
   */
  virtual const char* DefaultConfig() const { return "dNdetaConfig.C"; }

  UInt_t   fTriggerMask;   // Trigger mask 
  UInt_t   fFilterMask;    // Events to filter out 
  Double_t fMinIpZ;        // Least z--coordiante of interaction point
  Double_t fMaxIpZ;        // Largest z--coordiante of interaction point
  TAxis    fCentAxis;      // Centrality axis 
  TH1I*    fTriggers;      // Histogram of triggers
  TH1I*    fEventStatus;   // Histogram of event selection 
  TH1D*    fVertex;        // Vertex distribution of all events 
  TH1D*    fCent;          // Centrality distribution of all events
  TH1D*    fAccVertex;     // Vertex distribution of accepted events 
  TH2D*    fAccVertexXY;   // Vertex (x,y) distribution of accepted events 
  TH1D*    fAccCent;       // Centrality distribution of accepted events
  Bool_t   fFirstEvent;    // Information stored or not 
  Bool_t   fCloneList;     // Wether to clone sum list for results
  TList*   fSums;          // Output list of sums
  TList*   fResults;       // Output list of results

  ClassDef(AliBaseAODTask,2)
};
#endif
//
// Local Variables:
//  mode: C++
// End:
//
