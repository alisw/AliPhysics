#ifndef ALICOPYHEADERTASK_H
#define ALICOPYHEADERTASK_H
/**
 * @file   AliCopyHeaderTask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 10:56:38 2011
 * 
 * @brief  A task to copy the ESD header to AOD
 * 
 * @ingroup pwglf_forward_tasks 
 */
#include <AliAnalysisTaskSE.h>
class AliESDVertex;

/**
 * Task to copy header from ESD to AOD 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 * @ingroup pwglf_forward_aod
 */
class AliCopyHeaderTask : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor
   * 
   * @param name Name
   */
  AliCopyHeaderTask(const char* name="header") 
    : AliAnalysisTaskSE(name),
      fCalculateRefMult(true)
  {
    fBranchNames = "ESD:AliESDHeader.,AliESDRun.";
  }
  /** 
   * Copy constructor 
   * 
   * @param other Object to copy from 
   */
  AliCopyHeaderTask(const AliCopyHeaderTask& other) 
    : AliAnalysisTaskSE(other),
      fCalculateRefMult(other.fCalculateRefMult)
  {}
  /** 
   * Destructor
   */
  virtual ~AliCopyHeaderTask() {}
  /** 
   * Assignmen operator 
   * 
   * @param other Object to assing from 
   * 
   * @return Reference to this object 
   */
  AliCopyHeaderTask& operator=(const AliCopyHeaderTask& other) 
  {
    if (this == &other) return *this;
    AliAnalysisTaskSE::operator=(other);
    fCalculateRefMult = other.fCalculateRefMult;
    return *this;
  }
  /** 
   * @{ 
   * @name Implementation of interface methods
   */
  /** 
   * Create output objects.  Called at the beginning of each
   * (Grid/Proof) sub-job.
   * 
   */
  virtual void   UserCreateOutputObjects() {}
  /** 
   * Called at the beginning of the master job 
   * 
   */  
  virtual void   Init() {}
  /** 
   * Called at start of each (Grid/Proof) sub-job
   */
  virtual void   LocalInit() {Init();}
  /** 
   * Called on each event. 
   * 
   * @param option Not used 
   */
  virtual void   UserExec(Option_t *option);
  /** 
   * Called at the end of processing on the merged output. 
   * 
   * @param option Not used 
   */
  virtual void   Terminate(Option_t *option);
  /* @} */
  /** 
   * Set whether to calculate the reference multiplicity.  This is
   * optional because the calculation could be processor/memory heavy.
   * 
   * @param calc If true, calculate the reference multiplicity (number
   * of tracks and tracklets) in @f$|\eta|<0.5@f$ and
   * @f$|\eta|<0.8@f$.
   */
  virtual void SetCalculateRefMult(Bool_t calc=true) { fCalculateRefMult=calc;}
  /** 
   * Connect this task to the train. 
   * 
   * @return true on success. 
   */
  virtual Bool_t Connect();
protected:
  /** 
   * Copy an ESD primary vertex to the AOD 
   * 
   * @param aod   Where to copy to
   * @param vtx   Vertex (if any) to copy to
   * @param type  Type of vertex 
   */
  void CopyVertex(AliAODEvent& aod, const AliESDVertex* vtx, Int_t type); 

  Bool_t fCalculateRefMult; // Whether to calculate reference multiplicity
  ClassDef(AliCopyHeaderTask,2); // Task to copy header from ESD to AOD
};

#endif
/* 
 * Local Variables:
 *  mode: C++ 
 * End:
 */
