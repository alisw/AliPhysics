// 
// This class inspects the event 
//
#ifndef ALIFMDMCEVENTINSPECTOR_H
#define ALIFMDMCEVENTINSPECTOR_H
/**
 * @file   AliFMDMCEventInspector.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:03:40 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_aod
 */
#include "AliFMDEventInspector.h"
class AliMCEvent;
class TH2F;
class AliStack;

/** 
 * This class inspects the event 
 *
 * @par Input:
 *   - AliESDFMD object possibly corrected for sharing
 *
 * @par Output:
 *   - A histogram of v_z of events with triggers. 
 *   - A histogram of v_z of events with vertex and triggers 
 *   - A histogram of trigger counters 
 * 
 * Note, that these are added to the master output list 
 *
 * @par Corrections used: 
 *   - None
 *
 * @ingroup pwg2_forward_algo 
 * @ingroup pwg2_forward_mc
 * @ingroup pwg2_forward_aod
 */
class AliFMDMCEventInspector : public AliFMDEventInspector
{
public:
  /** 
   * Constructor 
   */
  AliFMDMCEventInspector();
  /** 
   * Constructor 
   * 
   * @param name Name of object
   */
  AliFMDMCEventInspector(const char* name);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDMCEventInspector(const AliFMDMCEventInspector& o);
  /** 
   * Destructor 
   */
  virtual ~AliFMDMCEventInspector();
  /** 
   * Assignement operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object
   */
  AliFMDMCEventInspector& operator=(const AliFMDMCEventInspector&);

  /** 
   * Initialize the object 
   * 
   * @param vtxAxis Vertex axis in use 
   */
  void Init(const TAxis& vtxAxis);
  /** 
   * Process MC truth event.  Note, returned values are the MC truth
   * values
   * 
   * @param event     Input event 
   * @param triggers  On return, the triggers fired 
   * @param ivz       On return, the found vertex bin (1-based).  A zero
   *                  means outside of the defined vertex range
   * @param vz        On return, the z position of the interaction
   * @param b         On return, impact parameter [fm] (if available)
   * @param npart     On return, number of participants (if available)
   * @param nbin      On return, number of binary collisions (if available)
   * @param phiR      On return, reaction plane angle (if available)
   * 
   * @return 0 (or kOk) on success, otherwise a bit mask of error codes 
   */
  UInt_t ProcessMC(AliMCEvent*       event, 
		   UInt_t&           triggers,
		   UShort_t&         ivz, 
		   Double_t&         vz,
		   Double_t&         b,
		   Int_t&            npart, 
		   Int_t&            nbin,
		   Double_t&         phiR);
  /** 
   * Compare the result of analysing the ESD for 
   * the inclusive charged particle density to analysing 
   * MC truth 
   * 
   * @param vz       Found @f$ v_z@f$
   * @param trueVz   True  @f$ v_z@f$
   * @param cent     Centrality
   * @param b        Impact parameter (if available)
   * @param npart    Number of participants (if available)
   * @param nbin     Number of binary collisions (if available)
   * 
   * @return true
   */
  virtual Bool_t CompareResults(Double_t vz,    Double_t trueVz, 
				Double_t cent,  Double_t b,
				Int_t    npart, Int_t    nbin);
  /** 
   * Store information about running conditions in output list 
   * 
   * The 3 TNamed objects from AliFMDEventInspector::StoreInformation
   * are defined.  In addition, a fourth TNamed object is defined.
   * The presence of this indicate MC data.
   *
   * - mc    Nothing special, and unique id is 1
   * 
   * @param runNo Run number 
   */
  virtual void StoreInformation(Int_t runNo);
  /** 
   * Read the production details 
   * 
   * @param event MC event
   */
  virtual void ReadProductionDetails(AliMCEvent* event);
protected:
  /** 
   * Read centrality from event 
   * 
   * @param esd  Event 
   * @param cent On return, the centrality or negative if not found
   * @param qual Quality flag 
   * 
   * @return False on error, true otherwise 
   */
  virtual Bool_t ReadCentrality(const AliESDEvent* esd, Double_t& cent,
				UShort_t& qual) const;
  /** 
   * Check if the event is single diffractive 
   * 
   * @param stack  Stack of MC particles 
   * @param xiMin  Lower cut off
   * @param xiMax  Upper cut off 
   * 
   * @return 
   */
  Bool_t IsSingleDiffractive(AliStack* stack,
			     Double_t xiMin=0, 
			     Double_t xiMax=1./81) const;
  TH1F* fHVertex;  // Histogram of vertex 
  TH1F* fHPhiR;    // Histogram of event plane 
  TH1F* fHB;       // Histogram of impact parameter 
  TH2F* fHBvsPart; // Impact parameter vs # participants 
  TH2F* fHBvsBin;  // Impact parameter vs # participants 
  TH2F* fHBvsCent; // Impact parameter vs centrality
  TH2F* fHVzComp;  // True vs reconstructed vz
  TH2F* fHCentVsPart; // Centrality versus # participants 
  TH2F* fHCentVsBin;  // Centrality versus # binary collisions 
  TString fProduction; // Production information 
  ClassDef(AliFMDMCEventInspector,3); // Inspect the event 
};

#endif
// Local Variables:
//   mode: C++
// End:



