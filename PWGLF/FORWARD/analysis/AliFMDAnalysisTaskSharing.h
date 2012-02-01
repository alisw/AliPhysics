#ifndef ALIFMDANALYSISTASKSHARING_H
#define ALIFMDANALYSISTASKSHARING_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               
 **/
 
#include "AliAnalysisTask.h"

#include "AliESDFMD.h"
//#include "TTree.h"
#include "AliESDEvent.h"
#include "AliFMDFloatMap.h"
class TChain;
class AliAODEvent;
class AliESDVertex;

/**
 * Do the sharing correction. 
 * 
 * This is the task to do the FMD sharing or hit merging.
 * It reads the input ESDFMD data and posts an ESDFMD object to
 * the tasks that must be performed after this task ie.
 * Density, BackgroundCorrection and Dndeta.
 *
 * Inputs: An AliESDFMD object
 *
 * Output:
 *   An AliESDFMD object, but with hits merged according to algorithm.
 *
 * Used correction objects: 
 *   Energy distribution fits (MPV and width of 1st landau)
 *   Hard low cut on 'mult' of 0.3
 *
 * Simplications: 
 *   Remove all diagnostics histograms except a few.  The histograms 
 *   needed for subsequent use are the ones for the sharing correction - 
 *   but only from MC data. 
 * 
 *   Remove calculation vertex efficiency.  This is best taken care of
 *   elsewhere. 
 *   
 * Open issues: 
 *   The ESD signal is un-angle-corrected and the after merging, 
 *   Re-angle-corrected.  I think this is wrong and will cause a 
 *   problem for the low-eta (low-theta) bins where the correction is
 *   largets.   Essentially, a particle that traverses two strips at 
 *   low theta will have a relatively large path through and will, 
 *   all things equal, deposite more energy.  The sharing filter may 
 *   then not pick this hit as steming from the same particle, but 
 *   but rather from 2 particles.   This is especially true of the MPV
 *   and width of the 1st Landau is determined from a full
 *   detector/ring spectra. 
 *
 * @ingroup FMD_ana
 * 
 * 
 */
class AliFMDAnalysisTaskSharing : public AliAnalysisTask
{
public:
  /** 
   * Constructor
   */
  AliFMDAnalysisTaskSharing();
  /** 
   * Constrictor
   * 
   * @param name Name of task 
   * @param SE   Whether we're run from an SE task or not
   */
  AliFMDAnalysisTaskSharing(const char* name, Bool_t SE = kTRUE);
  /** 
   * Destructor
   * 
   */
  virtual ~AliFMDAnalysisTaskSharing() {;}
  /** 
   * Copy constructor
   * 
   * @param o Object to copy from 
   */
  AliFMDAnalysisTaskSharing(const AliFMDAnalysisTaskSharing& o) 
    : AliAnalysisTask(),
      fDebug(o.fDebug),
      fESD(o.fESD),
      // fOutputESD(),
      foutputESDFMD(o.foutputESDFMD),
      fSharedThis(o.fSharedThis),
      fSharedPrev(o.fSharedPrev),
      fDiagList(),
      fStandalone(o.fStandalone),
      fEsdVertex(o.fEsdVertex),
      fStatus(o.fStatus),
      fLastTrackByStrip(o.fLastTrackByStrip),
      fLastOrbit(o.fLastOrbit) {}
  /** 
   * Assignment operator
   * 
   * @return Reference to this object
   */
  AliFMDAnalysisTaskSharing& 
  operator=(const AliFMDAnalysisTaskSharing&) { return *this; }
    
  /**
   * @{ 
   * @name Implementation of interface methods
   */
  virtual void ConnectInputData(Option_t *option = "");
  virtual void CreateOutputObjects();
  virtual void Init() {}
  virtual void LocalInit() {Init();}
  virtual void Exec(Option_t */*option*/);
  virtual void Terminate(Option_t* /* option*/);
  virtual void SetDebugLevel(Int_t level) {fDebug = level;}
  /** 
   * @} 
   */
  /** 
   * Get the multiplicity of a strip
   * 
   * @param mult    Previous(?) multiplicty
   * @param eta     Pseudo rapidity of strip
   * @param Eprev   Previous energy deposition
   * @param Enext   Next energy deposition
   * @param det     Detector
   * @param ring    Ring
   * @param sec     Sector
   * @param strip   Strip
   * 
   * @return 
   */
  Float_t GetMultiplicityOfStrip(Float_t  mult, 
				 Float_t  eta, 
				 Float_t  Eprev, 
				 Float_t  Enext, 
				 UShort_t det, 
				 Char_t   ring, 
				 UShort_t sec, 
				 UShort_t strip);
  // void GetVertex(Double_t* vertexXYZ) ;
  /** 
   * Set the Output data
   * 
   * @param fmd Output data
   */      
  void SetFMDData(AliESDFMD* fmd) {foutputESDFMD = fmd;}
  /** 
   * Set the output list 
   * 
   * @param outlist 
   */
  void SetOutputList(TList* outlist) {fDiagList = outlist;}
  /** 
   * Set the vertex
   * 
   * @param vertex 
   */
  void SetVertex(AliESDVertex* vertex) {fEsdVertex = vertex;}
  /** 
   * Set the input data
   *
   * @param esd Input
   */
  void SetInputESD(AliESDEvent* esd) {fESD = esd;}
  /** 
   * Get status flag
   * 
   * @return @c true on success
   */
  Bool_t GetEventStatus() const {return fStatus;}
  /** 
   * Get the vertex efficiency from data.  This is calculated as  
   * 
   * @f[
   *   e_{vtx} = \frac{# events with vertex}{# of events with trigger}
   * @f]
   * 
   * @return 
   */
  Float_t GetVtxEfficiencyFromData() ;
  /** 
   * Get the vertex efficiency from the data for NSD triggers
   * 
   * @return 
   */
  Float_t GetNSDVtxEfficiencyFromData() ;
    
 private:
  /** 
   * Calculate eta from theta
   * 
   * @param eta Input eta
   * 
   * @return Theta corresponding to eta
   */
  Float_t Eta2Theta(Float_t eta) const ;
  /** 
   * Get the psuedo-rapidity of a strip
   * 
   * @param det     Detector
   * @param ring    Ring
   * @param sector  Sector
   * @param strip   Strip
   * @param zvtx    Vertex position along beam-axis
   * 
   * @return Eta
   */
  Double_t EtaFromStrip(UShort_t det, 
			Char_t   ring, 
			UShort_t sector, 
			UShort_t strip, 
			Double_t zvtx);
  /** 
   * Process a primary particle (MC only)
   * 
   */
  void ProcessPrimary();
    
  Int_t          fDebug;            //  Debug flag
  AliESDEvent*   fESD;              //! ESD
  AliESDFMD*     foutputESDFMD;     // the output ESDFMD object
  Bool_t         fSharedThis;       // was this strip shared?
  Bool_t         fSharedPrev;       // was the previous strip shared?
  TList*         fDiagList;         // list of diag histos
  Bool_t         fStandalone;       // do we run standalone or in SE task
  AliESDVertex*  fEsdVertex;        // vtx info from the ESD
  Bool_t         fStatus;           // event status
  AliFMDFloatMap fLastTrackByStrip; // the last track to hit this strip
  UInt_t         fLastOrbit;
  
  ClassDef(AliFMDAnalysisTaskSharing, 0); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//  mode: C++
// End Variables:
