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
   * Set the @f$ \mathrm{IP}_z@f$ axis 
   * 
   * @param n   Number of bins
   * @param min Least value 
   * @param max Largest value 
   */
  void SetIPzAxis(Int_t n, Double_t min, Double_t max)
  {
    SetAxis(fIPzAxis, n, min, max);
  }
  /** 
   * Set the @f$ \mathrm{IP}_z@f$ axis 
   * 
   * @param n   Number of bins
   * @param max Largest absolute value 
   */
  void SetIPzAxis(Int_t n, Double_t max)
  {
    SetAxis(fIPzAxis, n, max);
  }
  /** 
   * Set the interaction point Z axis 
   * 
   * @param spec bin specification  
   */
  void SetIPzAxis(const TString& spec)
  {
    SetAxis(fIPzAxis, spec);
  }
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
  /** 
   * Set the very least centrality to consider.  This is for cases
   * where the centrality calibration of simulated data doesn't really
   * match the real-data one, but we want to keep the original
   * centrality bins.  E.g., for LHC15k1b1, the calibration of the
   * DPMJet centrality does not give the same mean number of tracklets
   * for the very most central bin.  We therefor need to rule out the
   * event with a very large number of tracklets from the sample by
   * setting this parameter to for example 0.1.  If this parameter is
   * set to a negative value (default) it is not considered.
   * 
   * @param x Absolute lowest centrality to consider 
   */
  void SetAbsMinCent(Double_t x=-1) { fAbsMinCent = x; }
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
  /** 
   * Set mask of events to filter out 
   *
   * @param mask 
   */
  void SetFilterMask(const char* mask);
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
   * @{ 
   * @name Service functions to set axis limits 
   */
  /** 
   * Fix axis attributes
   * 
   * @param axis  Axis to fix 
   * @param title Possible title for axis 
   */
  static void FixAxis(TAxis& axis, const char* title=0);
  /** 
   * Set an axis based on bin borders 
   * 
   * @param axis    Axis to set 
   * @param n       Number of bins 
   * @param borders Bin borders (n+1 entries)
   */
  static void SetAxis(TAxis& axis, Int_t n, Double_t* borders);
  /** 
   * Set an axis based on test string of specs.  The token separator
   * is given in @a sep.
   * 
   * @param axis Axis to set 
   * @param spec Specification
   * @param sep  Token separate 
   */
  static void SetAxis(TAxis& axis, const TString& spec, const char* sep=":,");
  /** 
   * Set axis with least and largest values
   * 
   * @param axis Axis to set
   * @param n    Number of bins 
   * @param l    Least value 
   * @param h    Largest value 
   */
  static void SetAxis(TAxis& axis, Int_t n, Double_t l, Double_t h);
  /** 
   * Set a symmetric axis 
   * 
   * @param axis Axis to set
   * @param n    Number of bins 
   * @param m    Maximum absolute value 
   */
  static void SetAxis(TAxis& axis, Int_t n, Double_t m);
  /** 
   * Print axis 
   * 
   * @param axis Axis to print 
   * @param nSig Number of significant digits 
   * @param alt  Alternative 
   */
  static void PrintAxis(const TAxis& axis, Int_t nSig=2, const char* alt=0);
  /* @} */
  /** 
   * @{ 
   * @name Get event information 
   */
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
   * Store the analysis trains name on the output pointed to by slot
   * number @a no.
   * 
   * @param no Output slot 
   */
  virtual Bool_t StoreTrainName(Int_t no);
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
  /* @} */
  /** 
   * Get the name of the default configuration script to use.
   * Sub-classes can override this to give another default
   * configuration script.  Note, it should problably only return the
   * base name (not full path) of the script.
   * 
   * @return Name of the configuration script to use. 
   */
  virtual const char* DefaultConfig() const { return "dNdetaConfig.C"; }

  /** Trigger mask */
  UInt_t   fTriggerMask;   
  /** Events to filter out */
  UInt_t   fFilterMask;    
  /** Centrality axis */
  TAxis    fCentAxis;      
  /** The absolute minimum centrality to consider  - for MC with poor match*/
  Double_t   fAbsMinCent;
  /** Collision point axis */
  TAxis    fIPzAxis;       
  /** Histogram of triggers */
  TH1I*    fTriggers;      
  /** Histogram of event selection */
  TH1I*    fEventStatus;   
  /** Vertex distribution of all events */
  TH1D*    fVertex;        
  /** Centrality distribution of all events */
  TH1D*    fCent;          
  /** Vertex distribution of accepted events */
  TH1D*    fAccVertex;     
  /** Vertex (x,y) distribution of accepted events */
  TH2D*    fAccVertexXY;   
  /** Centrality distribution of accepted events */
  TH1D*    fAccCent;       
  /** Information stored or not */
  Bool_t   fFirstEvent;    
  /** Wether to clone sum list for results */
  Bool_t   fCloneList;     
  /** Output list of sums */
  TList*   fSums;          
  /** Output list of results */
  TList*   fResults;       

  ClassDef(AliBaseAODTask,3)
};
#endif
//
// Local Variables:
//  mode: C++
// End:
//
