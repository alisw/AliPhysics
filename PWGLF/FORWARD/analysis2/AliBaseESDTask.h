//
// Base class for FMD ESD input - sub classes can use the services to
// easily setup a task to process the FMD ESD data 
//
#ifndef ALIFMDESDTASK_H
#define ALIFMDESDTASK_H

#include <AliAnalysisTaskSE.h>
class AliFMDEventInspector;
class AliCorrectionManagerBase;
class AliAODHandler;
class AliESDEvent;
class TAxis;
class TList;

/**
 * Base class for tasks that analyse the FMD ESD.  This wraps a
 * single-event analysis task, and provides a modified interface to
 * implement for the sub-classes:
 *
 * @code
 * class MyFMDESDTask : public AliBaseESDTask 
 * { 
 * public: 
 *   MyFMDESDTask() : AliBaseESDTask() {}
 *   MyFMDESDTask(const char* name) : AliBaseESDTask(name)
 *   { 
 *   }
 *   AliFMDEventInspector& GetEventInspector() { return fInspector; }
 *   const AliFMDEventInspector& GetEventInspector() const{ return fInspector; }
 *   Bool_t Book() 
 *   {
 *     fNeededCorrections = AliForwardCorrectionManager::kELossFits;
 *     fSharingFilter.CreateUserObject(fList);
 *     return true;
 *   }
 *   Bool_t PreData(const TAxis& vertex, const TAxis& eta)
 *   {
 *     fSharingFilter.SetupForData(eta);
 *     return true;
 *   }
 *   Bool_t PreEvent() { fESDFMD.Clear(); return true; } 
 *   Bool_t PostEvent() { return true; } 
 *   Bool_t Event(AliESDEvent& esd)
 *   {
 *     Bool_t   lowFlux   = kFALSE;
 *     UInt_t   triggers  = 0;
 *     UShort_t ivz       = 0;
 *     TVector3 ip;
 *     Double_t cent      = -1;
 *     UShort_t nClusters = 0;
 *     UInt_t   found     = fEventInspector.Process(esd, triggers, lowFlux, 
 *                                                  ivz, ip, cent, nClusters);
 *     if (found & AliFMDEventInspector::kNoEvent)    return false;
 *     if (found & AliFMDEventInspector::kNoTriggers) return false;
 *     if (found & AliFMDEventInspector::kNoSPD)      return;
 *     if (found & AliFMDEventInspector::kNoFMD)      return;
 *     if (found & AliFMDEventInspector::kNoVertex)   return;
 *     if (triggers & AliAODForwardMult::kPileUp)     return;
 *     if (found & AliFMDEventInspector::kBadVertex)  return;
 * 
 *     Bool_t ret = fSharingFilter.Filter(esd, lowFlux, fESDFMD, ip.Z());
 *     return ret;
 *   }
 *   Bool_t Finalize()
 *   {
 *     GetSharingFilter().Terminate(fList,fResults,Int_t(nTr));
 *     return true;
 *   }
 * protected:
 *   AliFMDEventInsepctor fInspector;
 *   AliFMDSharingFilter  fSharingFilter;
 *   AliESDFMD            fESDFMD;
 * };
 * @endcode 
 * 
 */
class AliBaseESDTask : public AliAnalysisTaskSE
{
public:
  /** 
   * Default (I/O) constructor - do not use directly
   */
  AliBaseESDTask();
  /** 
   * User constructor 
   * 
   * @param name  Name of the task 
   * @param title Class name used in configuration script 
   * @param manager Correction manager 
   */
  AliBaseESDTask(const char* name, const char* title,
		 AliCorrectionManagerBase* manager);
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
  virtual Bool_t Connect(const char* sumFile=0, 
			 const char* resFile=0)
  {
    return Connect(sumFile, resFile, false);
  }
  /** 
   * Add this task to the manager and connect the outputs.  If @a
   * sumFile is null or the empty string, then the sum container is
   * stored in the default output file of the manager.  If @a resFile
   * is null or the empty string, then it is set to @a resFile if
   * defined, otherwise to the default output file of the manager.
   * 
   * @param sumFile Output file for sums
   * @param resFile Output file for sums
   * @param old     Use old names
   * 
   * @return true on success 
   */
  virtual Bool_t Connect(const char* sumFile, 
			 const char* resFile,
			 Bool_t      old);
  /** 
   * Called when initializing the train 
   */
  virtual Bool_t Setup() { return true; }
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
   * Called on first event _before_ reading corrections.  Here, the
   * user class can do additional checking to see if the some (more or
   * less) corrections are needed.
   * 
   * @param esd Event 
   */
  virtual void PreCorrections(const AliESDEvent* esd);
  /** 
   * Called after reading in the first event. Here we can setup stuff
   * depending on the conditions we're running under.
   * 
   * @return true on success.  If this returns false, then we turn the
   * task into a zombie and we do no more processing.
   */
  virtual Bool_t PreData(const TAxis& vertex, const TAxis& eta);
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
   * @param esd Input event 
   * 
   * @return true on success 
   */
  virtual Bool_t Event(AliESDEvent& esd) = 0;
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
  virtual Bool_t Finalize() { return true; }
  /** 
   * @{ 
   * @name Utility methods 
   */
  /** 
   * Print information 
   * 
   * @param option Not used
   */
  virtual void Print(Option_t* option="") const;
  /** 
   * Set the debug level 
   * 
   * @param dbg 
   */
  virtual void SetDebug(Int_t dbg);
  /** 
   * Overload super class method for setting debug level to call our
   * SetDebug member function.
   * 
   * @param dbg Debug level (0: no output, 1: essentials, 3: a whole lot)
   */
  virtual void SetDebugLevel(Int_t dbg) 
  { 
    AliAnalysisTaskSE::SetDebugLevel(dbg); 
    SetDebug(dbg);
  }
  void SetIPzMethod(const char* str);
  /* @} */
  // --- Configuration etc -------------------------------------------
  /** @{ 
   * @name Access sub-components 
   */
  /** 
   * Configure this task via a macro 
   * 
   * @param macro Macro to configure va 
   * 
   * @return true on success, false otherwise
   */
  virtual Bool_t Configure(const char* macro="-default-");
  /** 
   * Get a reference to the event inspector. User must override this
   * to return proper object
   * 
   * @return Reference to the event inspector 
   */
  virtual AliFMDEventInspector& GetEventInspector() = 0;
  /** 
   * Get a reference to the event inspector. User must override this
   * to return proper object
   * 
   * @return Reference to the event inspector 
   */
  virtual const AliFMDEventInspector& GetEventInspector() const = 0;
  /* @} */
protected:
  /** 
   * Copy constructor - left undefined
   * 
   * @param o Object to copy from 
   */
  AliBaseESDTask(const AliBaseESDTask& o);
  /** 
   * Assignment operator - left undefined 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object.
   */
  AliBaseESDTask& operator=(const AliBaseESDTask& o);
  // --- Customize ---------------------------------------------------
  /** 
   * Evaluate wether this is for MC. User class can override this 
   * 
   * @return true if we're to initialize corrections for MC input
   */
  virtual Bool_t IsMC() const { return false; }
  /** 
   * Set the default eta axis to use in case we didn't get one from
   * the read-in corretions.  Override this if the sub class should go
   * on even without a valid eta axis from the corrections (e.g. QA
   * task)
   * 
   * @return null
   */
  virtual TAxis* DefaultEtaAxis() const;
  /** 
   * Set the default eta axis to use in case we didn't get one from
   * the read-in corretions.  Override this if the sub class should go
   * on even without a valid eta axis from the corrections (e.g. QA
   * task)
   * 
   * @return null
   */
  virtual TAxis* DefaultVertexAxis() const;
  /** 
   * Get the correction mananger.  Derived class should overload this
   * to return the proper object.
   * 
   * @return Pointer to correction manager
   */
  virtual AliCorrectionManagerBase* GetManager() const { return fCorrManager; }
  /** 
   * Get the correction mananger.  Derived class should overload this
   * to return the proper object.
   * 
   * @return Pointer to correction manager
   */
  virtual AliCorrectionManagerBase* GetManager() { return fCorrManager; }

  // --- Task methods ------------------------------------------------
  /** 
   * @{ 
   * @name Task interface methods 
   */
  /** 
   * Initialize the task 
   */
  void LocalInit();
  /** 
   * Create output objects 
   */
  void UserCreateOutputObjects();
  /** 
   * Process each event 
   *
   * @param option Not used
   */  
  void UserExec(Option_t* option);
  /** 
   * End of job
   * 
   * @param option Not used 
   */
  void Terminate(Option_t* option);
  /** 
   * @} 
   */
  // --- Services for derived classes --------------------------------
  /**
   * Create output branches - called from UserCreateOutputObjects
   */
  virtual void CreateBranches(AliAODHandler*) {}
  /**
   * Mark this event as one to store in the AOD 
   * 
   */
  virtual void MarkEventForStore() const;
  /** 
   * Check if all needed corrections are there and accounted for.  If not,
   * do a Fatal exit 
   * 
   * @param what Which corrections is needed
   * 
   * @return true if all present, false otherwise
   */  
  virtual Bool_t CheckCorrections(UInt_t what) const;
  /** 
   * Read corrections
   * 
   * 
   * @param pe  On return, the eta axis
   * @param pv  On return ,the vertex axis 
   * @param mc  True assume MC input
   * @param sat True if we need for satellite interactions too 
   * 
   * @return true ons succcss
   */
  virtual Bool_t ReadCorrections(const TAxis*& pe, 
				 const TAxis*& pv,
				 Bool_t mc=false,
				 Bool_t sat=false);
  /**
   * Get the ESD event. IF this is the first event, initialise
   *
   * @return Pointer to ESD event structore 
   */
  virtual AliESDEvent* GetESDEvent();
  /** 
   * Store the analysis trains name on the output pointed to by slot
   * number @a no.
   * 
   * @param no Output slot 
   */
  virtual Bool_t StoreTrainName(Int_t no);
  /** 
   * Get default configuration script name
   * 
   * @return Script name
   */
  virtual const char* DefaultConfig() const 
  {
    return "ForwardAODConfig.C";
  }

  // --- Members -----------------------------------------------------
  Bool_t fFirstEvent;        // Wheter we're waiting for the first event
  TList* fList;              // Output list 
  TList* fResults;           // Results list 
  UInt_t fNeededCorrections; // Set this to bit-mask of corrections we need
  UInt_t fExtraCorrections;  // Set this to bit-mask of corrections we'd like
  Bool_t fCloneList;         // Result list is a clone of sum list
private:
  /**
   * A pointer to the corrections manager.  This is here to make the
   * corrections manager persistent - that is, when we write the
   * analysis train to a file (as done in PROOF) we should also write
   * down the corrections mananger.   This pointer ensures that. 
   * 
   */
  AliCorrectionManagerBase* fCorrManager; // Pointer to corrections manager

  ClassDef(AliBaseESDTask,1);
};
inline void AliBaseESDTask::PreCorrections(const AliESDEvent*) {}
#endif
// Local Variables:
//   mode: C++
// End:
