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
    : AliAnalysisTaskSE(name)
  {}
  /** 
   * Copy constructor 
   * 
   * @param other Object to copy from 
   */
  AliCopyHeaderTask(const AliCopyHeaderTask& other) 
    : AliAnalysisTaskSE(other)
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
    AliAnalysisTaskSE::operator=(other);
    return *this;
  }
  /** 
   * @{ 
   * @name Implementation of interface methods
   */
  virtual void   UserCreateOutputObjects() {}
  virtual void   Init() {}
  virtual void   LocalInit() {Init();}
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *option);
  /* @} */

  ClassDef(AliCopyHeaderTask,1); // Task to copy header from ESD to AOD
};

#endif
/* 
 * Local Variables:
 *  mode: C++ 
 * End:
 */
