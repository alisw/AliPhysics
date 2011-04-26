// 
// Calculate the corrections in the central regions
// 
#ifndef ALICENTRALMCCORRECTIONS_H
#define ALICENTRALMCCORRECTIONS_H
/**
 * @file   AliCentralMCCorrectionsTask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:05:51 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_central_aod
 */
#include <AliAnalysisTaskSE.h>
#include "AliFMDMCEventInspector.h"
#include "AliSPDMCTrackDensity.h"
#include <TH1I.h>
class AliCentralCorrSecondaryMap;
class AliESDEvent;
class TH2D;
class TH1D;
class TList;


/** 
 * Calculate the corrections in the central regions
 * 
 * @par Inputs: 
 *   - AliESDEvent 
 *
 * @par Outputs: 
 *   - AliAODCentralMult 
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 * 
 * @ingroup pwg2_central_tasks
 * @ingroup pwg2_central_mc
 * @ingroup pwg2_central_aod
 * 
 */
class AliCentralMCCorrectionsTask : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   */
  AliCentralMCCorrectionsTask(const char* name);
  /** 
   * Constructor
   */
  AliCentralMCCorrectionsTask();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliCentralMCCorrectionsTask(const AliCentralMCCorrectionsTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliCentralMCCorrectionsTask& operator=(const AliCentralMCCorrectionsTask& o);
  /** 
   * @{ 
   * @name Interface methods 
   */
  /** 
   * Initialize the task 
   * 
   */
  virtual void Init();
  /** 
   * Create output objects 
   * 
   */
  virtual void UserCreateOutputObjects();
  /** 
   * Process each event 
   *
   * @param option Not used
   */  
  virtual void UserExec(Option_t* option);
  /** 
   * End of job
   * 
   * @param option Not used 
   */
  virtual void Terminate(Option_t* option);
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
   * Set the number of phi bins to use 
   * 
   * @param nBins 
   */
  void SetNPhiBins(UShort_t nBins) { fNPhiBins = nBins; }
  /** 
   * Get a reference to the track density calculator 
   * 
   * @return Reference to the track density calculator 
   */
  AliSPDMCTrackDensity& GetTrackDensity() { return fTrackDensity; }
  /** 
   * Get a reference to the track density calculator 
   * 
   * @return Reference to the track density calculator 
   */
  const AliSPDMCTrackDensity& GetTrackDensity() const { return fTrackDensity; }
  /** 
   * Get a reference to the event inspector
   * 
   * @return Reference to the event inspector 
   */
  AliFMDMCEventInspector& GetEventInspector() { return fInspector; }
  /** 
   * Get a reference to the event inspector
   * 
   * @return Reference to the event inspector 
   */
  const AliFMDMCEventInspector& GetEventInspector() const { return fInspector;}
protected: 
  /**
   * A vertex bin 
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
     */
    VtxBin(Double_t low, Double_t high, const TAxis& etaAxis, UShort_t nPhi);
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    VtxBin(const VtxBin& o);
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this object
     */
    VtxBin& operator=(const VtxBin& o);
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
    void DefineOutput(TList* list);
    /** 
     * End of job process 
     * 
     * @param o List to add output to 
     */
    void Finish(const TList* i, 
		TList* o,
		UShort_t iVz, 
		AliCentralCorrSecondaryMap* map);

    TH2D* fHits;    // Cache of per-ring histograms
    TH2D* fPrimary;  // Cache or primary 
    TH1D* fCounts;   // Event count 

    ClassDef(VtxBin,1); // Vertex bin 
  };
  /** 
   * Define our vertex bins 
   * 
   * @param list List to read or add binst from/to
   */
  void DefineBins(TList* list);

  AliFMDMCEventInspector fInspector; // Event inspector 
  AliSPDMCTrackDensity   fTrackDensity; // Get the track density 

  TObjArray* fVtxBins;      // Vertex bins 
  Bool_t     fFirstEvent;   // First event flag 
  TH1I*      fHEvents;      // All Events
  TH1I*      fHEventsTr;    // Histogram of events w/trigger
  TH1I*      fHEventsTrVtx; // Events w/trigger and vertex 
  TAxis      fVtxAxis;      // Vertex axis 
  TAxis      fEtaAxis;      // Eta axis 
  TList*     fList;         // Output list 
  UShort_t   fNPhiBins;     // Nunber of phi bins
  ClassDef(AliCentralMCCorrectionsTask,1) // Central corrections class
};

#endif
// Local Variables:
//  mode: C++
// End:

