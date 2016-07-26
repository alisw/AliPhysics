// 
// Calculate the corrections in the base regions
// 
#ifndef ALIBASEMCCORRECTIONS_H
#define ALIBASEMCCORRECTIONS_H
/**
 * @file   AliBaseMCCorrectionsTask.h
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Mar 23 14:05:51 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_base_aod
 */
#include <AliBaseESDTask.h>
#include <AliESDFMD.h>
#include "AliFMDMCEventInspector.h"
#include <TH1I.h>
class AliBaseMCTrackDensity;
class AliCorrectionManagerBase;
class AliESDEvent;
class TH2D;
class TH1D;
class TList;
class TVector3;

/** 
 * Calculate the corrections in the base regions
 * 
 * @par Inputs: 
 *   - AliESDEvent 
 *
 * @par Outputs: 
 *   - AliAODBaseMult 
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 * 
 * @ingroup pwglf_base_tasks
 * @ingroup pwglf_base_mc
 * @ingroup pwglf_base_aod
 * 
 */
class AliBaseMCCorrectionsTask : public AliBaseESDTask
{
public:
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   * @param m    Manager 
   */
  AliBaseMCCorrectionsTask(const char* name,
			   AliCorrectionManagerBase* m);
  /** 
   * Constructor
   */
  AliBaseMCCorrectionsTask();
  /** 
   * @{ 
   * @name Interface methods 
   */
  /** 
   * Create output objects 
   * 
   * @return true on success
   */
  virtual Bool_t Book();
  /** 
   * Process each event 
   *
   * @param esd ESD event
   * 
   * @return true on success
   */  
  virtual Bool_t Event(AliESDEvent& esd);
  /** 
   * End of job
   * 
   * @return true on success
   */
  virtual Bool_t Finalize();
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
   * Set the vertex axis to use
   * 
   * @param nBins Number of bins
   * @param vzMin Least @f$z@f$ coordinate of interation point
   * @param vzMax Largest @f$z@f$ coordinate of interation point
   */
  void SetVertexAxis(Int_t nBins, Double_t vzMin, Double_t vzMax=-1000000);
  /** 
   * Set the vertex axis to use
   * 
   * @param axis Axis
   */
  void SetVertexAxis(const TAxis& axis);
  /** 
   * Set the eta axis to use
   * 
   * @param nBins Number of bins
   * @param etaMin Least @f$\eta@f$ 
   * @param etaMax Largest @f$\eta@f$ 
   */
  void SetEtaAxis(Int_t nBins, Double_t etaMin, Double_t etaMax=-1000000);
  /** 
   * Set the eta axis to use
   * 
   * @param axis Axis
   */
  void SetEtaAxis(const TAxis& axis);
  /** 
   * Set-up for satellite collisions 
   * 
   * @param sat If true, set-up for satellites 
   */
  void SetSatellite(Bool_t sat);
  /** 
   * Get a reference to the track density calculator 
   * 
   * @return Reference to the track density calculator 
   */
  virtual AliBaseMCTrackDensity& GetTrackDensity() = 0;
  /** 
   * Get a reference to the track density calculator 
   * 
   * @return Reference to the track density calculator 
   */
  virtual const AliBaseMCTrackDensity& GetTrackDensity() const = 0;
  /** 
   * Get a reference to the event inspector
   * 
   * @return Reference to the event inspector 
   */
  AliFMDEventInspector& GetEventInspector() { return fInspector; }
  /** 
   * Get a reference to the event inspector
   * 
   * @return Reference to the event inspector 
   */
  const AliFMDEventInspector& GetEventInspector() const { return fInspector;}
  /**
   * setter for the fUseESDVertexCoordinate flag 
   */
  void SetUseESDVertex(Bool_t use){ fUseESDVertex = use;}
  /**
   * setter for the fCalculateafterESDeventcuts flag 
   */
  void SetAfterEventSel(Bool_t use){ fAfterEventSel = use; }
protected: 
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliBaseMCCorrectionsTask(const AliBaseMCCorrectionsTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliBaseMCCorrectionsTask& operator=(const AliBaseMCCorrectionsTask& o);
  /**
   * A vertex bin.  These are only used internally and are not
   * streamed. We can therefore leave them out of the streamer (i.e.,
   * not define a dictionary for them).
   * 
   */
  struct VtxBin : public TNamed
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
     * @param nPhi      Number of @f$\varphi@f$ bins 
     */
    VtxBin(Double_t low, Double_t high, const TAxis& etaAxis, UShort_t nPhi);
    virtual ~VtxBin() {}
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
     * Get bin name
     * 
     * @param low       Lower @f$v_z@f$ bound
     * @param high      Upper @f$v_z@f$ bound
     * 
     * @return Bin name 
     */
    static const char* BinName(Double_t low, Double_t high);
    /** 
     * Declare output in passed list 
     * 
     * @param list List to put output in 
     */
    virtual TList* CreateOutputObjects(TList* list);
    TH2D*                  fPrimary;  // Cache or primary 
    TH1D*                  fCounts;   // Event count 

    // ClassDef(VtxBin,1); // Vertex bin 
  };
  virtual TAxis* DefaultVertexAxis() const { return const_cast<TAxis*>(&fVtxAxis); }
  virtual TAxis* DefaultEtaAxis() const { return const_cast<TAxis*>(&fEtaAxis); }
  /** 
   * Create a vertex bin 
   * 
   * @param low      Low cut on @f$IP_{z}@f$
   * @param high     High cut on @f$IP_{z}@f$ 
   * 
   * @return Newly create vertex bin 
   */
  virtual VtxBin* CreateVtxBin(Double_t low, Double_t high) = 0;
  /** 
   * Process an ESD event 
   * 
   * @param esd ESD event 
   * @param mc  MC event 
   * @param bin Vertex bin
   * @param ip  @f$IP_{z}@f$ 
   * 
   * @return true on success 
   */
  virtual Bool_t ProcessESD(const AliESDEvent& esd, 
			    const AliMCEvent& mc, 
			    VtxBin& bin,
			    const TVector3& ip) = 0;
  /** 
   * Create corrections objects and store them in passed list
   * 
   * @param results Output list 
   */
  virtual void CreateCorrections(TList* results) = 0;
  /** 
   * Do the final processing of a vertex bin 
   * 
   * @param bin       Vertex bin
   * @param iVz       Vertex bin number 
   * 
   * @return true on successd
   */
  virtual Bool_t FinalizeVtxBin(VtxBin*      bin, 
				UShort_t     iVz) = 0;
  /** 
   * Define our vertex bins 
   * 
   * @param list List to read or add binst from/to
   */
  void DefineBins(TList* list);

  AliFMDMCEventInspector fInspector; // Event inspector 

  TObjArray* fVtxBins;      // Vertex bins 
  TH1I*      fHEvents;      // All Events
  TH1I*      fHEventsTr;    // Histogram of events w/trigger
  TH1I*      fHEventsTrVtx; // Events w/trigger and vertex 
  TAxis      fVtxAxis;      // Vertex axis 
  TAxis      fEtaAxis;      // Eta axis 
  Bool_t     fUseESDVertex; // if true use Z vertex from ESD in calculations
  Bool_t     fAfterEventSel; //if true corr. be calc. after events selection

  ClassDef(AliBaseMCCorrectionsTask,1) // Base corrections class
};

#endif
// Local Variables:
//  mode: C++
// End:
