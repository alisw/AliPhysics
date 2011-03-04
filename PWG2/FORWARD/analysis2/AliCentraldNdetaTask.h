//
// Task to analyse the AOD for for dN/deta in the central regions 
//
#ifndef ALICENTRALDNDETATASK_H
#define ALICENTRALDNDETATASK_H
#include "AliBasedNdetaTask.h"
class TList;
class TH2D;
class TH1D;

/**
 * Task to determine the 
 */
class AliCentraldNdetaTask : public AliBasedNdetaTask
{
public:
  /** 
   * Constructor 
   * 
   */
  AliCentraldNdetaTask() : AliBasedNdetaTask() {}
  /** 
   * Constructor
   * 
   * @param name    Name of task 
   * @param maxVtx  Set @f$v_z@f$ range
   */
 AliCentraldNdetaTask(const char*) 
   : AliBasedNdetaTask("Central") 
  { 
    fSymmetrice = false; 
    fCorrEmpty  = false;
  }
  /**
   * Destructor
   * 
   */
  virtual ~AliCentraldNdetaTask() {}
protected:
  /** 
   * Get the histogram 
   * 
   * @param aod 
   * @param mc 
   * 
   * @return 
   */
  TH2D* GetHistogram(const AliAODEvent* aod, Bool_t mc=false);

  ClassDef(AliCentraldNdetaTask,1); // Determine multiplicity in central area
};

#endif
// Local Variables:
//   mode: C++
// End:
 
