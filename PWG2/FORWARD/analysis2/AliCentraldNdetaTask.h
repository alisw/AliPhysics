//
// Task to analyse the AOD for for dN/deta in the central regions 
//
#ifndef ALICENTRALDNDETATASK_H
#define ALICENTRALDNDETATASK_H
/**
 * @file   AliCentraldNdetaTask.h
 * @author Hans Hjersing Dalsgaard
 * @date   Wed Mar 23 13:59:26 2011
 * 
 * @brief  
 * 
 * @ingroup pwg2_forward_dndeta
 * 
 */
#include "AliBasedNdetaTask.h"
class TList;
class TH2D;
class TH1D;

/**
 * Tasks to determine @f$ dN/d\eta@f$ in the forward regions
 *
 * @ingroup pwg2_forward_tasks_dndeta
 * @ingroup pwg2_forward_dndeta
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
   * @param name Name of task - ignored
   */
  AliCentraldNdetaTask(const char* name);
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
  /** 
   * Get the colour to use for markers
   * 
   * @return Marker colour 
   */
  virtual Int_t GetColor() const { return kRed+1; }
  /** 
   * Get the marker style 
   * 
   * @return Marker style 
   */
  virtual Int_t GetMarker() const { return 21; }

  ClassDef(AliCentraldNdetaTask,1); // Determine multiplicity in central area
};

#endif
// Local Variables:
//   mode: C++
// End:
 
