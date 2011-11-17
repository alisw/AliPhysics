// Code to analyse dN/deta from the forward analysis
// This can plot the results 
// Also works for MC data 
#ifndef ALIFMDDNDETA_H
#define ALIFMDDNDETA_H

#include "TObject.h"
// #include "TList.h"
#include "TString.h"
class TList;
class TH1F;
class TH3F;
class TProfile3D;
class TProfile2D;
class TH3D;
/**
 * This class creates dN/deta for the FMD from the analysis objects.
 * The contents of this class should probably go into a task at some point
 */
class AliFMDDndeta : public TObject
{
 public:
  /** 
   * 
   * 
   */
  AliFMDDndeta();
  /** 
   * 
   * 
   * @param o 
   */
  AliFMDDndeta(const AliFMDDndeta& o) 
    : TObject(),
      fList(0),
      fMultList(),
      fNbinsToCut(o.fNbinsToCut),
      fVtxCut1(o.fVtxCut1),
      fVtxCut2(o.fVtxCut2),
      fIsInit(o.fIsInit),
      fIsGenerated(),
      fPrimEvents(o.fPrimEvents),
      fEvents(o.fEvents),
      fPrimdNdeta(o.fPrimdNdeta),
      fDrawAll(kFALSE)
  {
    for (int i = 0; i < 5; i++) fAnalysisNames[i] = "";
  }
  /** 
   * 
   * 
   * 
   * @return 
   */  
  AliFMDDndeta& operator=(const AliFMDDndeta& /*o*/) 
  {
    // Assignment operator 
    return (*this);
  }
  /**
   * Enumeration of analysis types
   * 
   */
  enum Analysis {kHits, kHitsTrVtx, kMult, kMultTrVtx, kMultNSD};

  /** 
   * Initialise 
   * 
   * @param filename 
   */  
  void Init(const Char_t* filename);
  /** 
   * Initialise 
   * 
   * @param list 
   */
  void Init(TList* list);
  /** 
   * Generate the multiplicity for analysis type @a what
   * 
   * @param what 
   */
  void GenerateMult(Analysis what);
  /** 
   * Draw the result.
   * 
   * @param what 
   * @param rebin 
   * @param realdata 
   * @param filename 
   */
  void DrawDndeta(Analysis what, 
		  Int_t rebin = 1, 
		  Bool_t realdata = kFALSE, 
		  TString filename = "none");
  /** 
   * Set the number of bins to cut
   * 
   * @param nbins 
   */
  void SetNbinsToCut(Int_t nbins) {fNbinsToCut = nbins;}
  /** 
   * Set the vertex cut
   * 
   * @param vtxcut 
   */
  void SetVtxCut1(Int_t vtxcut) {fVtxCut1 = vtxcut;}
  /** 
   * Set the vertex cut
   * 
   * @param vtxcut 
   */
  void SetVtxCut2(Int_t vtxcut) {fVtxCut2 = vtxcut;}
  /** 
   * Whether to draw all 
   * 
   * @param drawall 
   */
  void SetDrawAll(Bool_t drawall) {fDrawAll = drawall;}
  /** 
   * Create sharing efficiency from file 
   * 
   * @param filename 
   * @param store 
   */
  void CreateSharingEfficiency(const Char_t* filename, Bool_t store = kFALSE);
  /** 
   * Get the list of multiplicities for a given type of analysis. 
   * 
   * @param what 
   * 
   * @return 
   */
  TList* GetMultList(Analysis what) const {return fMultList[what];}
 private:
  void GenerateHits(Analysis what);
  void SetNames(Analysis what);
  const char* GetAnalysisName(Analysis what, UShort_t det, 
			      Char_t ring, Int_t vtxbin);
  const char* GetPrimName(Analysis what, UShort_t det, 
			  Char_t ring, Int_t vtxbin);
  void   RebinHistogram(TH1F* hist, Int_t rebin);
  TList*  fList;                         // A list of input histograms
  TList*  fMultList[5];                  // A list of mult histograms 
  Int_t   fNbinsToCut;                   // The number of bins to cut
  Int_t   fVtxCut1;                      // Vtx low
  Int_t   fVtxCut2;                      // Vtx high
  Bool_t  fIsInit;                       // Are we init ? 
  Bool_t  fIsGenerated[5];               // Have we generated ?
  TString fPrimEvents;                   // Number of prim events
  TString fEvents;                       // Number of events
  TString fPrimdNdeta;                   // the primary dNdeta from MC
  TString fAnalysisNames[5];             // Names of analysis
  // TProfile3D*   fDataObject;                  // New data object
  Bool_t fDrawAll;                        //Draw relevant or all
  //TH3D*   fDataObject;                  // New data object
  
  ClassDef(AliFMDDndeta,2);
};


#endif
// Local Variables:
//  mode: C++
// End:
