#ifndef ALIFMDMCHITENERGYFITTERTASK_H
#define ALIFMDMCHITENERGYFITTERTASK_H
#include <AliBaseESDTask.h>
#include <AliFMDMCHitEnergyFitter.h>
#include <AliFMDEventInspector.h>
class AliMCAuxHandler;

/** 
 * This task is designed to read in MC truth information about the
 * energy loss in each strip, and then fit the distributions from
 * secondaries and primaries separately.
 * 
 * Then (this is not implemented yet) it tries to deconvolve the
 * contributions from secondaries and primaries separately from the
 * total sum energy loss distribution, and in that way estimate the
 * secondary to primary charge particle ratios per @f$\eta@f$ bin.  If
 * the same procedure is applied to real data, then we can have an
 * estimate - from data - of the ratio of secondaries to primaries and
 * thus perhaps get a clearer picture of the secondary particle
 * contamination.
 *
 * The the function fitted to the (scaled) energy loss
 * (@f$\Delta/\Delta_{mip}@f$ distributions is 
 *
 * @f[
 *    F(\Delta;\Delta_p,\xi,\sigma,\mathbf{a}) = 
 *    \sum_{i=1}^{N}a_i f(\Delta;\Delta_i,\xi_i,\sigma_i)
 * @f]
 * where @f$ a@f$ is of length @f$ N@f$ and 
 * @f[
 *   f(\Delta;\Delta,\xi,\sigma) = 
 *     \int_{-\infty}^{+\infty}d\Delta' L(\Delta,\Delta',\xi)
 *     \frac{1}{\sqrt{2\pi\sigma^2}}
 *     e^{-\frac{(\Delta'-\Delta_mp^2}{2\sigma^2}}
 * @f]
 * and 
 * @f[ 
 *     \Delta_i = i(\Delta_1+\xi\log(i))\\
 *     \xi_i    = i\xi_1\\
 *     \sigma_i = \sqrt{i}\sigma_1\\
 *     a_1      = 1\\
 *     \Delta_p = \Delta_1\\
 *     \xi      = \xi_1\\
 *     \sigma   = \sigma_1
 * @f]
 *
 * See also AliLandauGausFitter 
 */
class AliFMDMCHitEnergyFitterTask : public AliBaseESDTask
{
public:
  /** 
   * Default CTOR - do not use
   */
  AliFMDMCHitEnergyFitterTask()
    : AliBaseESDTask(),
      fEventInspector(), 
      fEnergyFitter(), 
      fHitHandler(0)
  {}
  /** 
   * CTOR
   * 
   * @param name     Name - not used
   * @param useTuple Whether to use store an NTuple
   */
  AliFMDMCHitEnergyFitterTask(const char* name, 
			      Bool_t      useTuple=false);
  /** 
   * DTOR
   */
  ~AliFMDMCHitEnergyFitterTask() {}

  /** 
   * Called when setting up the train on the client - i.e., called
   * once before the job hits the worker nodes
   * 
   * @return true
   */  
  Bool_t Setup();
  /** 
   * Called at start-up on the clients - i.e., called once per worker
   * node before the event processing.
   * 
   * @return true
   */
  Bool_t Book();
  /** 
   * Called after the first event was seen - i.e., called once per
   * worker node just after the first event was seen, but before event
   * processing.
   * 
   * @param ipz Interaction point Z--coordinate axis to use 
   * @param eta @f$\eta=-\log[\tan^{-1}(\theta/2)]@f$ axis to use
   * 
   * @return true
   */
  Bool_t PreData(const TAxis& ipz, const TAxis& eta);
  /** 
   * Process a single event 
   * 
   * @param esd ESD input event 
   * 
   * @return true on success
   */
  Bool_t Event(AliESDEvent& esd);
  /** 
   * Finalize the task.  This is called once after all event
   * processing and after all outputs have been merged.  This is
   * called on the master/client once.
   * 
   * @return true on success
   */
  Bool_t Finalize() { fEnergyFitter.Fit(fResults); return true; }
  /** 
   * Print information to standard out
   * 
   * @param option Passed to sub objects as-is
   */
  void   Print(Option_t* option="") const;
  /** 
   * Get the event inspector 
   * 
   * @return Reference to the event inspector 
   */
  AliFMDEventInspector& GetEventInspector() { return fEventInspector; }
  /** 
   * Get the event inspector 
   * 
   * @return Constant reference to the event inspector 
   */
  const AliFMDEventInspector& GetEventInspector() const 
  { 
    return fEventInspector; 
  }
  /** 
   * Get the energy fitter 
   * 
   * @return Reference to the energy fitter 
   */
  AliFMDMCHitEnergyFitter& GetEnergyFitter() { return fEnergyFitter; }
  /** 
   * Get the energy fitter 
   * 
   * @return Constant reference to the energy fitter 
   */
  const AliFMDMCHitEnergyFitter& GetEnergyFitter() const 
  { 
    return fEnergyFitter; 
  }
  void SetDebug(Int_t dbg) { 
    AliBaseESDTask::SetDebug(dbg); GetEnergyFitter().SetDebug(dbg); }
protected:
  /** 
   * Copy CTOR - not Implemented
   * 
   * @param o Object to copy from 
   */  
  AliFMDMCHitEnergyFitterTask(const AliFMDMCHitEnergyFitterTask& o);
  /** 
   * Assignment operator - not implemented
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  AliFMDMCHitEnergyFitterTask& operator=(const AliFMDMCHitEnergyFitterTask&o);
  /** 
   * Get the default @f$\eta=-\log[\tan^{-1}(\theta/2)]@f$ axis to use 
   * 
   * @return Pointer to static object 
   */  
  TAxis* DefaultEtaAxis() const;
  /** 
   * Get the default interaction point Z--coordinate axis to use 
   * 
   * @return Pointer to static object 
   */  
  TAxis* DefaultVertexAxis() const;
  
  AliFMDEventInspector    fEventInspector; // The event inspector
  AliFMDMCHitEnergyFitter fEnergyFitter;   // The energy loss fitter 
  AliMCAuxHandler*        fHitHandler;     // Handler for reading hits in

  ClassDef(AliFMDMCHitEnergyFitterTask,1); // Task to fit Delta from MC hits
};

#endif
// Local Variables:
//  mode: C++
// End:
