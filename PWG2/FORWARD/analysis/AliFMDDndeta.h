#ifndef ALIFMDDNDETA_H
#define ALIFMDDNDETA_H

#include "TObject.h"
#include "TList.h"
#include "TString.h"
class TH1F;

//
// This class creates dN/deta for the FMD from the analysis objects.
// The contents of this class should probably go into a task at some point
//

class AliFMDDndeta : public TObject
{

 public:
  AliFMDDndeta();
 AliFMDDndeta(const AliFMDDndeta& o) : 
  TObject(),
    fList(0),
    fMultList(),
    fNbinsToCut(o.fNbinsToCut),
    fVtxCut1(o.fVtxCut1),
    fVtxCut2(o.fVtxCut2),
    fIsInit(o.fIsInit),
    fIsGenerated(),
    fPrimEvents(o.fPrimEvents),
    fEvents(o.fEvents),
    fPrimdNdeta(fPrimdNdeta)
    {}
  
  AliFMDDndeta& operator=(const AliFMDDndeta& /*o*/) 
    {
      // Assignment operator 
      
      return (*this);
    }
  
  enum Analysis {kHits, kHitsTrVtx, kMult, kMultTrVtx};
  
  void Init(const Char_t* filename);
  void Init(TList* list);
  void GenerateMult(Analysis what);
  void DrawDndeta(Analysis what, Int_t rebin = 1, Bool_t realdata = kFALSE);
  void SetNbinsToCut(Int_t nbins) {fNbinsToCut = nbins;}
  void SetVtxCut1(Int_t vtxcut) {fVtxCut1 = vtxcut;}
  void SetVtxCut2(Int_t vtxcut) {fVtxCut2 = vtxcut;}
  void CreateSharingEfficiency(const Char_t* filename, Bool_t store = kFALSE);
  TList* GetMultList() {return &fMultList;}
 private:
  void GenerateHits(Analysis what);
  void SetNames(Analysis what);
  const char* GetAnalysisName(Analysis what, UShort_t det, Char_t ring, Int_t vtxbin);
  const char* GetPrimName(Analysis what, UShort_t det, Char_t ring, Int_t vtxbin);
  void   RebinHistogram(TH1F* hist, Int_t rebin);
  TList* fList;                          // A list of input histograms
  TList  fMultList;                      // A list of mult histograms 
  Int_t  fNbinsToCut;                    // The number of bins to cut
  Int_t  fVtxCut1;                       // Vtx low
  Int_t  fVtxCut2;                       // Vtx high
  Bool_t fIsInit;                        // Are we init ? 
  Bool_t fIsGenerated[4];                // Have we generated ?
  TString fPrimEvents;                   // Number of prim events
  TString fEvents;                       // Number of events
  TString fPrimdNdeta;                   // the primary dNdeta from MC
  
  ClassDef(AliFMDDndeta,2);
};


#endif
// Local Variables:
//  mode: C++
// End Variables;
