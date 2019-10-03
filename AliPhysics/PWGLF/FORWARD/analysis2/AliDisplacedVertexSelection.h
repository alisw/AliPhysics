/**
 * @file   AliDisplacedVertexSelection.h
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 00:53:03 2013
 * 
 * @brief  Selection of events from satellite interactions 
 * 
 * @ingroup pwglf_forward_aod
 * 
 */
#ifndef ALIDISPLACEDVERTEXSELECTION_H
#define ALIDISPLACEDVERTEXSELECTION_H
#include <TObject.h>
class AliESDEvent;
class AliMCEvent;
class TH1D;
class TH2D;
class TList;

/** 
 * Selection of events from satellite interactions 
 *
 * @ingroup pwglf_forward_aod
 */
class AliDisplacedVertexSelection : public TObject 
{
public:
  /** 
   * Constructor 
   */
  AliDisplacedVertexSelection();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliDisplacedVertexSelection(const AliDisplacedVertexSelection& o);
  /** 
   * Assingment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliDisplacedVertexSelection& operator=(const AliDisplacedVertexSelection& o);
  /** 
   * Define the output 
   * 
   * @param l     List to add output to
   * @param name  Name of the list 
   * @param mc    True if we're looking at MC data
   */
  void SetupForData(TList* l, const char* name=0, Bool_t mc=false);
  /** 
   * Print information 
   * 
   * @param option  Not used 
   */
  void Print(Option_t* option="") const;
  /** 
   * Process an ESD event to get the information 
   * 
   * @param esd ESD event 
   * 
   * @return true on success
   */
  Bool_t Process(const AliESDEvent* esd);
  /** 
   * Process an MC event to find true satellite vertex
   * 
   * @param mcevent MC event
   * 
   * @return true if found or not MC input, false in case of problems
   */
  Bool_t ProcessMC(const AliMCEvent* mcevent);
  /**
   * Check if this event is marked as a satellite interaction 
   *
   * @return true if the found vertex isn't invalid
   */
  Bool_t IsSatellite() const { return fVertexZ != kInvalidVtxZ; }
  /** 
   * Get the interaction point Z-coordinate from ZDC timing. 
   * 
   * 
   * @return Interaction point Z-coordinate
   */
  Double_t GetVertexZ() const { return fVertexZ; }
  /** 
   * Return the centrality percentile 
   * 
   * 
   * @return Centrality percentile (ZDC vs ZEM)
   */
  Double_t GetCentralityPercentile() const { return fCent; }
protected:
  Bool_t CheckOutlier(Int_t ivtx, const AliESDEvent* esd) const;
  Float_t GetZemCorr(Int_t k, Bool_t minusminus) const;
  
  enum { 
    kMaxK        = 10,
    kInvalidVtxZ = 9999
  };
  TH1D*     fDeltaTdc;
  TH1D*     fSumTdc;
  TH1D*     fZdcEnergy;
  TH1D*     fZemEnergy;
  TH2D*     fCorrelationZemZdc;
  TH2D*     fCorrelationSumDelta;  
  Double_t  fVertexZ;  // Interaction point Z-coordinate
  Double_t  fCent;     // Centrality percentile
  TH1D*     fHVertexZ; // Histogram of vertices 
  TH1D*     fHCent;    // Histogram of centrality 
  Bool_t    fMC; // MC flag
  
  ClassDef(AliDisplacedVertexSelection,4); // Satelitte collisions 
};

#endif
// Local Variables: 
//  mode: C++
// End:
