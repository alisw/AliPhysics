#ifndef ALIROOT_PWG2_FORWARD_ALIFORWARDCORRECTIONS_H
#define ALIROOT_PWG2_FORWARD_ALIFORWARDCORRECTIONS_H
#include <AliAnalysisTaskSE.h>
#include "AliForwardUtil.h"
#include "AliFMDSharingFilter.h"
#include "AliFMDDensityCalculator.h"
#include "AliFMDCorrections.h"
#include "AliFMDHistCollector.h"
#include "AliAODForwardMult.h"
#include <AliESDFMD.h>
#include <TH1I.h>
class AliFMDAnaParameters;
class AliESDEvent;
class TH2D;
class TH3D;
class TList;
class TTree;


/** 
 * Calculate the corrections in the forward regions
 * 
 * @par Inputs: 
 *   - AliESDEvent 
 *
 * @par Outputs: 
 *   - AliAODForwardMult 
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 * 
 * @ingroup pwg2_forward_analysis 
 * 
 */
class AliForwardMCCorrectionsTask : public AliAnalysisTaskSE
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
  void         Print(Option_t* option="") const;

  void SetVertexAxis(Int_t nBins, Double_t vzMin, Double_t vzMax=-1000000);
  void SetVertexAxis(const TAxis& axis);
  void SetEtaAxis(Int_t nBins, Double_t etaMin, Double_t etaMax=-1000000);
  void SetEtaAxis(const TAxis& axis);
protected: 
  TH2D*  GetVertexProj(Int_t v, TH3D* src) const;
  TH3D* Make3D(const char* name, const char* title, Int_t nPhi) const;
  TH1D* Make1D(const char* name, const char* title) const;
  void  FillPrimary(Bool_t gotInel, Bool_t gotVtx, 
		    Double_t vz, Double_t eta, Double_t phi);
  void FillStrip(UShort_t d, Char_t r, 
		 Double_t vz, Double_t eta, Double_t phi,
		 Bool_t first);
  TH1I*  fHEvents;           // All Events
  TH1I*  fHEventsTr;         // Histogram of events w/trigger
  TH1I*  fHEventsTrVtx;      // Events w/trigger and vertex 
  TH1I*  fHEventsVtx;        // Events w/vertex 
  TH1I*  fHTriggers;         // Triggers
  TH3D*  fPrimaryInnerAll;   // Distribution of primaries - all events
  TH3D*  fPrimaryOuterAll;   // Distribution of primaries - all events
  TH3D*  fPrimaryInnerTrVtx; // Distribution of primaries - trg+vtx events
  TH3D*  fPrimaryOuterTrVtx; // Distribution of primaries - trg+vtx events
  TH3D*  fHitsFMD1i;         // Distribution of FMD1i hits 
  TH3D*  fHitsFMD2i;         // Distribution of FMD2i hits 
  TH3D*  fHitsFMD2o;         // Distribution of FMD2o hits 
  TH3D*  fHitsFMD3i;         // Distribution of FMD3i hits 
  TH3D*  fHitsFMD3o;         // Distribution of FMD3o hits 
  TH1D*  fStripsFMD1i;       // Distribution of FMD1i # strips hit
  TH1D*  fStripsFMD2i;       // Distribution of FMD2i # strips hit 
  TH1D*  fStripsFMD2o;       // Distribution of FMD2o # strips hit 
  TH1D*  fStripsFMD3i;       // Distribution of FMD3i # strips hit 
  TH1D*  fStripsFMD3o;       // Distribution of FMD3o # strips hit 
  TAxis  fVtxAxis;           // Vertex axis 
  TAxis  fEtaAxis;           // Eta axis 

  TList* fList; // Output list 

  ClassDef(AliForwardMCCorrectionsTask,1) // Forward corrections class
};

#endif
// Local Variables:
//  mode: C++
// End:

