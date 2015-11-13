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
 * @ingroup pwglf_central_aod
 */
#include "AliBaseMCCorrectionsTask.h"
#include "AliSPDMCTrackDensity.h"
class AliCentralCorrSecondaryMap;
class AliCentralCorrAcceptance;


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
 * @ingroup pwglf_central_tasks
 * @ingroup pwglf_central_mc
 * @ingroup pwglf_central_aod
 * 
 */
class AliCentralMCCorrectionsTask : public AliBaseMCCorrectionsTask
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
   * @{ 
   * @name Interface methods 
   */
  /** 
   * Print this object 
   * 
   * @param option   Not used
   */
  void Print(Option_t* option="") const;

  /** 
   * Set the number of phi bins to use 
   * 
   * @param nBins 
   */
  void SetNPhiBins(UShort_t nBins) { fNPhiBins = nBins; }
  /** 
   * Whether to make effective corrections
   * 
   * @param e if true, make effective correction
   */
  void SetEffectiveCorrection(Bool_t e) { fEffectiveCorr = e; }
  /** 
   * Set the maximum @f$|\eta|@f$ to accept. 
   * 
   * @param maxEta maximum @f$|\eta|@f$
   */
  void SetEtaCut(Double_t maxEta=1.9) { fEtaCut = maxEta; }
  /** 
   * If a particular phi bin has less then this fraction of the
   * largest signal in the corresponding @f$\eta@f$ slice, then it is
   * not used.
   * 
   * @param least Lower bound on fraction of largest signal in this
   * @f$\eta@f$ slice
   */
  void SetAcceptanceCut(Double_t least=0.8) { fCorrCut = least; }
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
protected: 
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
   * A vertex bin.  These are only used internally and are never
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
     * @param nPhi      Number of phi bins 
     */
    VtxBin(Double_t low, Double_t high, const TAxis& etaAxis, UShort_t nPhi);
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
     * End of job process 
     * 
     * @param o         List to add output to 
     * @param i         Input list
     * @param iVz       Vertex bin 
     * @param effective Make an effective correction 
     * @param etaCut    Maximum @f$|\eta|@f$ to use 
     * @param accCut    Cut on acceptance 
     * @param acorr     Acceptance correction 
     * @param map       Corrections map 
     */
    void Terminate(const TList* i, 
		   TList* o,
		   UShort_t iVz, 
		   Bool_t effective,
		   Double_t etaCut, 
		   Double_t accCut,
		   AliCentralCorrSecondaryMap* map,
		   AliCentralCorrAcceptance* acorr);
    
    TH2D* fHits;     // Cache of MC-truth hits
    TH2D* fClusters; // Cache of reconstructed hits

    // ClassDef(VtxBin,3); // Vertex bin 
  };
  /** 
   * Define our vertex bins 
   * 
   * @param list List to read or add binst from/to
   */
  void DefineBins(TList* list);
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
  virtual Bool_t FinalizeVtxBin(AliBaseMCCorrectionsTask::VtxBin*      bin, 
				UShort_t     iVz);


  AliSPDMCTrackDensity        fTrackDensity; // Get the track density 
  AliCentralCorrSecondaryMap* fSecCorr;
  AliCentralCorrAcceptance*   fAccCorr;

  UShort_t   fNPhiBins;      // Nunber of phi bins
  Bool_t     fEffectiveCorr; // Whether to make effective corrections
  Double_t   fEtaCut;        // Maximum Eta
  Double_t   fCorrCut;       // Correction cut
  ClassDef(AliCentralMCCorrectionsTask,3) // Central corrections class
};

#endif
// Local Variables:
//  mode: C++
// End:

