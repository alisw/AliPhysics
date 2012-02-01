#ifndef ALICOPYHEADERTASK_H
#define ALICOPYHEADERTASK_H
/**
 * @file   AliCopyHeaderTask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 10:56:38 2011
 * 
 * @brief  A task to copy the ESD header to AOD
 * 
 * @ingroup pwg2_forward_tasks 
 */
#include <AliAnalysisTaskSE.h>

/**
 * Task to copy header from ESD to AOD 
 * 
 * @ingroup pwg2_forward_scripts_tasks
 * @ingroup pwg2_forward_aod
 */
class AliCopyHeaderTask : public AliAnalysisTaskSE
{
public:
  AliCopyHeaderTask(const char* name="header") 
    : AliAnalysisTaskSE(name)
  {}
  AliCopyHeaderTask(const AliCopyHeaderTask& other) 
    : AliAnalysisTaskSE(other)
  {}
  virtual ~AliCopyHeaderTask() {}
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

  ClassDef(AliCopyHeaderTask,1);
};

#endif
/* 
 * Local Variables:
 *  mode: C++ 
 * End:
 */
