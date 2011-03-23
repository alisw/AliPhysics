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
 * @ingroup pwg2_forward_aod
 */
#include <AliAnalysisTaskSE.h>
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
 * @ingroup pwg2_forward_tasks
 * @ingroup pwg2_forward_mc
 * @ingroup pwg2_forward_aod
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
  /** 
   * Print this object 
   * 
   * @param option   Not used
   */
  void         Print(Option_t* option="") const;

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
protected: 
  /** 
   * Get vertex project
   * 
   * @param v   Vertex bin 
   * @param src Source 3D histogram 
   * 
   * @return 2D projection of the V'th bin
   */
  TH2D*  GetVertexProj(Int_t v, TH3D* src) const;
  /** 
   * Make a 3D histogram
   * 
   * @param name   Name 
   * @param title  Title 
   * @param nPhi   Number of phi bins
   * 
   * @return Histogram
   */
  TH3D* Make3D(const char* name, const char* title, Int_t nPhi) const;
  /** 
   * Make 1D histogram
   * 
   * @param name   Name 
   * @param title  Title
   * 
   * @return Histogram
   */
  TH1D* Make1D(const char* name, const char* title) const;
  /** 
   * Fill in primary information
   * 
   * @param gotInel   Got INEL trigger from ESD
   * @param gotVtx    Got vertex Z from ESD 
   * @param vz        @f$z@f$ coordinate of interation point
   * @param eta       Pseudo rapidity 
   * @param phi       Azimuthal angle
   */
  void  FillPrimary(Bool_t gotInel, Bool_t gotVtx, 
		    Double_t vz, Double_t eta, Double_t phi);
  /** 
   * Fill in per-strip information
   * 
   * @param d         Detector
   * @param r         Ring
   * @param vz        @f$z@f$ coordinate of interation point
   * @param eta       Pseudo rapidity 
   * @param phi       Azimuthal angle
   * @param first     First fill in this event
   */
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

