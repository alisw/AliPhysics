#ifndef ALIFMDANALYSISTASKSE_H
#define ALIFMDANALYSISTASKSE_H

#include "AliAnalysisTaskSE.h"
#include "AliFMDAnalysisTaskSharing.h"
#include "AliFMDAnalysisTaskDensity.h"
#include "AliFMDAnalysisTaskBackgroundCorrection.h"
#include "AliFMDAnalysisTaskDndeta.h"
#include "AliFMDAnalysisTaskBFCorrelation.h"

class AliFMDAnaParameters;

/** @defgroup FMD_ana Analysis tasks 
    @brief Various classes to do analysis tasks 
*/

/**
 * @class AliFMDAnalysisTaskSE
 * @brief Collected analysis task 
 * @ingroup FMD_ana
 *
 * Collector of various analysis tasks.  It will do the full chain of 
 * analysis tasks:
 *
 * - AliFMDAnalysisTaskSharing
 * - AliFMDAnalysisTaskDensity
 * - AliFMDAnalysisTaskBackgroundCorrection
 * - AliFMDAnalysisTaskDndeta
 * - AliFMDAnalysisTaskBFCorrelation
 */

class AliFMDAnalysisTaskSE : public AliAnalysisTaskSE
{
public:				// 
  AliFMDAnalysisTaskSE();
  AliFMDAnalysisTaskSE(const char* name);
  virtual ~AliFMDAnalysisTaskSE() {;}
  AliFMDAnalysisTaskSE(const AliFMDAnalysisTaskSE& o) 
    : AliAnalysisTaskSE(),
      fListOfHistos(o.fListOfHistos),
      fSharing(o.fSharing),
      fDensity(o.fDensity),
      fBackground(o.fBackground),
      fDndeta(o.fDndeta),   
      fBFCorrelation(o.fBFCorrelation), 
      fParams(o.fParams)
  {}
  AliFMDAnalysisTaskSE& operator=(const AliFMDAnalysisTaskSE&) { return *this; }
  
  virtual void Init();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* /*option*/);
  void         Terminate(Option_t */*option*/);
  void         SetParams(AliFMDAnaParameters* p) { fParams = p; }
  void         Print(Option_t* option="") const ;
  
private:
  
  TList*                                 fListOfHistos;  // Output list
  AliFMDAnalysisTaskSharing              fSharing;       // Sharing task
  AliFMDAnalysisTaskDensity              fDensity;       // Density task
  AliFMDAnalysisTaskBackgroundCorrection fBackground;    // Background task
  AliFMDAnalysisTaskDndeta               fDndeta;        // dN/deta task
  AliFMDAnalysisTaskBFCorrelation        fBFCorrelation; // BF correlation task
  AliFMDAnaParameters*                   fParams;        // Analysis parameters
  
  ClassDef(AliFMDAnalysisTaskSE, 1);

};
#endif
// Local Variables:
//  mode: C++ 
// End:
