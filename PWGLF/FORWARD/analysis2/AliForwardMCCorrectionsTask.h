// 
// Calculate the corrections in the forward regions
// 
#ifndef ALIFORWARDMCCORRECTIONS_H
#define ALIFORWARDMCCORRECTIONS_H
/**
 * @file   AliForwardMCCorrectionsTask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:05:51 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_aod
 */
#include "AliBaseMCCorrectionsTask.h"
#include <AliESDFMD.h>
#include "AliFMDMCTrackDensity.h"
#include "AliForwardUtil.h"
#include <TH1I.h>
class AliESDEvent;
class AliFMDCorrSecondaryMap;
class TList;
class TVector3;

/** 
 * Calculate the simulation-based corrections in the forward regions
 * 
 * @image html alice-int-2012-040-secondary_origin.png "Fraction of secondaries"
 *
 * @par Inputs: 
 *   - AliESDEvent (for steering only)
 *   - AliTrackReference
 *   - Kinematics 
 *   - Geometry 
 *
 * @par Outputs: 
 *   - AliFMDCorrSecondaryMap
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 *   - None
 * 
 * @ingroup pwglf_forward_tasks
 * @ingroup pwglf_forward_mc
 * @ingroup pwglf_forward_aod
 * 
 */
class AliForwardMCCorrectionsTask : public AliBaseMCCorrectionsTask
{
public:
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   */
  AliForwardMCCorrectionsTask(const char* name);
  /** 
   * Constructor
   */
  AliForwardMCCorrectionsTask();
  /** 
   * @{ 
   * @name Interface methods 
   */
  /** 
   * Called before the event processing 
   * 
   * @return true on success
   */
  virtual Bool_t PreEvent();
  /** 
   * @} 
   */
  /** 
   * Print this object 
   * 
   * @param option   Not used
   */
  void Print(Option_t* option="") const;
  /** 
   * Get a reference to the track density calculator 
   * 
   * @return Reference to the track density calculator 
   */
  AliBaseMCTrackDensity& GetTrackDensity() { return fTrackDensity; }
  /** 
   * Get a reference to the track density calculator 
   * 
   * @return Reference to the track density calculator 
   */
  const AliBaseMCTrackDensity& GetTrackDensity() const { return fTrackDensity; }
protected: 
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardMCCorrectionsTask(const AliForwardMCCorrectionsTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardMCCorrectionsTask& operator=(const AliForwardMCCorrectionsTask& o);
  /**
   * A vertex bin.   These are only used internally and are never
   * streamed.
   * 
   */
  struct VtxBin : public AliBaseMCCorrectionsTask::VtxBin
  {
    /** 
     * Constructor 
     */
    VtxBin();
    /** 
     * Constructor
     *  
     * @param low       Lower @f$v_z@f$ bound
     * @param high      Upper @f$v_z@f$ bound
     * @param etaAxis   @f$\eta@f$ axis to use 
     */
    VtxBin(Double_t low, Double_t high, const TAxis& etaAxis);
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    VtxBin(const VtxBin& o){;}
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this object
     */
    VtxBin& operator=(const VtxBin& o){return *this;}
    /** 
     * Declare output in passed list 
     * 
     * @param list List to put output in 
     */
    TList* CreateOutputObjects(TList* list);
    /** 
     * Calculate the background correction
     * 
     * @param hits      Summed hits (track-refs)
     * @param primary   Summed primaries 
     * 
     * @return Background correction
     */
    TH2D* MakeBg(const TH2D* hits, const TH2D* primary) const;
    /** 
     * End of job process 
     * 
     * @param o   List to add output to 
     * @param i   Input 
     * @param iVz Vertex bin
     * @param map Corrections map 
     */
    void Terminate(const TList*            i, 
		TList*                  o,
		UShort_t                iVz, 
		AliFMDCorrSecondaryMap* map);

    AliForwardUtil::Histos fHists;    // Cache of per-ring histograms
    // ClassDef(VtxBin,2); // Vertex bin 
  };
  /** 
   * Create a vertex bin 
   * 
   * @param low     Low cut on @f$IP_{z}@f$ 
   * @param high    High cut on @f$IP_{z}@f$ 
   * 
   * @return Newly created vertex bin
   */
  AliBaseMCCorrectionsTask::VtxBin* CreateVtxBin(Double_t low, Double_t high);
  /** 
   * Process an ESD event
   * 
   * @param esd   ESD event 
   * @param mc    MC event
   * @param bin   Vertex bin 
   * @param ip    @f$IP_{z}@f$ 
   * 
   * @return true on success
   */
  Bool_t ProcessESD(const AliESDEvent& esd, const AliMCEvent& mc, 
		    AliBaseMCCorrectionsTask::VtxBin& bin,
		    const TVector3& ip);
  /** 
   * Create corrections objects and store them in passed list
   * 
   * @param results Output list 
   */
  virtual void CreateCorrections(TList* results);
  /** 
   * Do the final processing of a vertex bin 
   * 
   * @param bin       Vertex bin
   * @param iVz       Vertex bin number 
   * 
   * @return true on successd
   */
  virtual Bool_t FinalizeVtxBin(AliBaseMCCorrectionsTask::VtxBin* bin, 
				UShort_t     iVz);


  AliFMDMCTrackDensity    fTrackDensity; // Get the track density 
  AliESDFMD               fESDFMD;       // Cache object
  AliFMDCorrSecondaryMap* fSecCorr;
  ClassDef(AliForwardMCCorrectionsTask,4) // Forward corrections class
};

#endif
// Local Variables:
//  mode: C++
// End:

