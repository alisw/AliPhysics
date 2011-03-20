// 
// This class inspects the event 
//
#ifndef ALIFMDMCEVENTINSPECTOR_H
#define ALIFMDMCEVENTINSPECTOR_H
#include "AliFMDEventInspector.h"
class AliMCEvent;
class TH2F;

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
   * @param lowFlux   On return, true if the event is considered a low-flux 
   *                  event (according to the setting of fLowFluxCut) 
   * @param ivz       On return, the found vertex bin (1-based).  A zero
   *                  means outside of the defined vertex range
   * @param vz        On return, the z position of the interaction
   * @param cent      On return, the centrality (in percent) or < 0 
   *                  if not found
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
   * @param esd 
   * @param mc 
   * 
   * @return 
   */
  virtual Bool_t CompareResults(Double_t vz,    Double_t trueVz, 
				Double_t cent,  Double_t b,
				Int_t    npart, Int_t    nbin);
protected:
  /** 
   * Read centrality from event 
   * 
   * @param esd  Event 
   * @param cent On return, the centrality or negative if not found
   * 
   * @return False on error, true otherwise 
   */
  virtual Bool_t ReadCentrality(const AliESDEvent* esd, Double_t& cent,
				UShort_t& qual) const;

  TH1F* fHVertex;  // Histogram of vertex 
  TH1F* fHPhiR;    // Histogram of event plane 
  TH1F* fHB;       // Histogram of impact parameter 
  TH2F* fHBvsPart; // Impact parameter vs # participants 
  TH2F* fHBvsBin;  // Impact parameter vs # participants 
  TH2F* fHBvsCent; // Impact parameter vs centrality
  TH2F* fHVzComp;  // True vs reconstructed vz
  TH2F* fHCentVsPart; // Centrality versus # participants 
  TH2F* fHCentVsBin;  // Centrality versus # binary collisions 
  ClassDef(AliFMDMCEventInspector,2); // Inspect the event 
};

#endif
// Local Variables:
//   mode: C++
// End:



