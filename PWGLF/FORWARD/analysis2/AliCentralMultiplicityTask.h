// 
// Base class for classes that calculate the multiplicity in the
// SPD clusters event-by-event
// 
#ifndef ALICENTRALMULTIPLICITYTASK_H
#define ALICENTRALMULTIPLICITYTASK_H
/**
 * @file   AliCentralMultiplicityTask.h
 * @author Hans Hjersing Dalsgaard
 * @date   Wed Mar 23 14:00:03 2011
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_aod
 * 
 */
#include <AliAnalysisTaskSE.h>
#include "AliFMDEventInspector.h"
#include "AliAODCentralMult.h"
class AliCentralCorrectionManager;
class AliESDEvent;
class AliMultiplicity;
class TH2D;
class TList;
class TTree;
class TObjArray;

/** 
 * Class that calculates the multiplicity in the
 * central region event-by-event
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
 * @ingroup pwglf_forward_tasks
 * @ingroup pwglf_forward_aod
 * 
 */
class AliCentralMultiplicityTask : public AliAnalysisTaskSE
{
public:
  /** 
   * @{ 
   * @name Interface methods 
   */
   /** 
   * Constructor 
   * 
   * @param name Name of task 
   */
  AliCentralMultiplicityTask(const char* name); 
  /** 
   * Constructor 
   *
   * Reserved for ROOT's I/O system - do not use
   */
  AliCentralMultiplicityTask();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliCentralMultiplicityTask(const AliCentralMultiplicityTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliCentralMultiplicityTask& operator=(const AliCentralMultiplicityTask& o);
  /** 
   * Configure this task via a macro 
   * 
   * @param macro Macro to configure va 
   * 
   * @return true on success, false otherwise
   */
  virtual Bool_t Configure(const char* macro="CentralAODConfig.C");
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
   * Print information 
   * 
   * @param option Not used
   */
  virtual void Print(Option_t* option="") const;
  /** 
   * Set whether to use the secondary corrections 
   * 
   * @param use Whether to use secondary corrections 
   */
  virtual void SetUseSecondary(Bool_t use) { fUseSecondary = use; }
  /** 
   * Set whether to use the acceptance corrections 
   * 
   * @param use Whether to use acceptance corrections 
   */
  virtual void SetUseAcceptance(Bool_t use) { fUseAcceptance = use; }
  /** 
   * Set whether to make diagnostics or not
   * 
   * @param use If true, store some extra diagnostic histograms
   */
  virtual void SetMakeDiagnostics(Bool_t use=true) { fStore = use; }
  /** 
   * Get the event inspector
   * 
   * @return Reference to used event inspector
   */
  AliFMDEventInspector& GetInspector() { return fInspector; }
  /** 
   * Get the event inspector
   * 
   * @return Reference to used event inspector
   */
  const AliFMDEventInspector& GetInspector() const { return fInspector; }
  /** 
   * Get the event inspector
   * 
   * @return Reference to used event inspector
   */
  AliFMDEventInspector& GetEventInspector() { return fInspector; }
  /** 
   * Get the event inspector
   * 
   * @return Reference to used event inspector
   */
  const AliFMDEventInspector& GetEventInspector() const { return fInspector; }

protected:
  /** 
   * Get the ESD event and initialise manager on first event if not
   * done already
   * 
   * @return Pointer to valid ESD event object 
   */
  virtual AliESDEvent* GetESDEvent();
  /** 
   * Mark this event for storage in AOD output
   * 
   */
  virtual void MarkEventForStore() const;
  /** 
   * Process the ESD SPD information 
   * 
   * @param hist    Histogram to fill
   * @param spdmult SPD multiplicity object
   */
  virtual void ProcessESD(TH2D& hist, 
			  const AliMultiplicity* spdmult) const;
  /** 
   * Find our eta limits
   * 
   */
  virtual void FindEtaLimits();
  /**
   * A vertex bin. 
   *
   * Used to store and apply corrections and fiducial cuts
   */
  struct VtxBin : public TObject
  {
    VtxBin(Int_t iVz=0, Double_t minIpZ=0, Double_t maxIpZ=0);
    VtxBin(const VtxBin& o);
    VtxBin& operator=(const VtxBin& o);
    
    const char* GetName() const;
    void SetupForData(TList* l, TH2* coverage, Bool_t store=true);
    void Correct(TH2D&  aodHist,
		 Bool_t useSecondary,
		 Bool_t useAcceptance,
		 Bool_t sum=true) const;
    void Print(Option_t* option="") const;

    Int_t        fId;     // Vertex bin number 
    Double_t     fMinIpZ; // Least value of ipZ 
    Double_t     fMaxIpZ; // Largest value of ipZ 
    Int_t        fEtaMin; // Smallest eta bin to use 
    Int_t        fEtaMax; // Largest eta bin to use 
    TH2*         fSec;    // Our secondary correction
    TH1*         fAcc;    // Our acceptance correction 
    mutable TH2* fHits;   // Diagnostics sum 

    ClassDef(VtxBin,1);
  };
    
protected: 
  /** 
   * Make a simple @f$\frac{dN_{ch}}{d\eta}@f$ estimate. 
   * 
   * @param input   Sum list
   * @param output  Output list 
   * @param nTr     On return, the number of events w/triggers
   * @param nTrVtx  On return, the number of events w/triggers+vertex
   * @param nAcc    On return, the number of accepted events 
   * 
   * @return true on success, false otherwise 
   */
  Bool_t MakeSimpledNdeta(const TList* input, 
			  TList*       output,
			  Double_t&    nTr, 
			  Double_t&    nTrVtx, 
			  Double_t&    nAcc);
  AliFMDEventInspector   fInspector;        // Inspect events 
  TList*                 fList;             // Output list
  AliAODCentralMult      fAODCentral;       // Output object
  Bool_t                 fUseSecondary;     // Whether to secondary map
  Bool_t                 fUseAcceptance;    // Whether to use acceptance corr.
  Bool_t                 fFirstEventSeen;   // Have we seen first event     
  Int_t                  fIvz;              // Event's vertex bin 
  TH2D*                  fNClusterTracklet; // # of clusters vs tracklets 
  TH2D*                  fClusterPerTracklet; // Clusters per tracklet. 
  TH1D*                  fNCluster;         //! Number of clusters 
  TH1D*                  fNTracklet;        //! number of tracklets 
  TObjArray*             fVtxList;          //! Array of vertex bins
  Bool_t                 fStore;            // Store diagnostics
  TH2D*                  fHData;            // Sum of signals 
private:
  AliCentralCorrectionManager* fCorrManager; 
  ClassDef(AliCentralMultiplicityTask,4)    // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

