#ifndef ALIDIELECTRONEVTVSTRKHIST_H
#define ALIDIELECTRONEVTVSTRKHIST_H

/**
 * @Author: Pascal Dillenseger <pascaldillenseger>
 * @Date:   2017-08-09, 17:39:28
 * @Email:  pdillens@cern.ch
 * @Filename: AliDielectronEvtVsTrkHist.h
 * @Last modified by:   pascaldillenseger
 * @Last modified time: 2017-09-12, 14:22:00
 */

#include <TArray.h>
#include <THashList.h>
#include <THnSparse.h>
#include <TNamed.h>
#include <TObject.h>

#include "AliVEvent.h"
#include "AliVTrack.h"

#include "AliDielectronHistos.h"
#include "AliDielectronVarManager.h"


//______________________________________________
class AliDielectronEvtVsTrkHist : public TNamed {

public:
  // Sparse enumerator
  enum eEvtVsTrkSparses
  {
    kSparseMatchEffITSTPC = 0,
    kSparseNMaxValues
  };

  // Histogram enumerator
  enum eEvtVsTrk3DHistos
  {
    k3DHistoNMaxValues = 0
  };

// Init
  AliDielectronEvtVsTrkHist(const char *name = "EvtVsTrkHist", const char *title = "EvtVsTrkHist");
  virtual ~AliDielectronEvtVsTrkHist();

  void Init();

  Float_t GetVarValueEvent(const AliVEvent *ev, Int_t var);
  Float_t GetTrackValue(AliVTrack *trk, Int_t var);

  const char* GetSparseUniqueNameInfo( Int_t iName);
  const char* Get3DHistoUniqueNameInfo( Int_t iName);

  void FillHistograms(const AliVEvent *ev);

  void SetHistogramList(AliDielectronHistos *dieHistos) {fHistoList = (THashList*) dieHistos->GetHistogramList()->FindObject("EvtVsTrk"); AliDielectronEvtVsTrkHist::Init();} // All histograms from the EvtVsTrk list are used
  Bool_t SetHistogramObj (TObject *obj);
  Bool_t SetSparseObj(TObject *obj);

  void SetSparseVars(Int_t nSparse, Int_t nAxis, Int_t type);
  void SetEventplaneAngles(const AliVEvent *ev);

  void SetPIDResponse(AliPIDResponse *pidResponse) {fPIDResponse = pidResponse;}

  Bool_t IsSelectedITS(AliVTrack *trk);
  Bool_t IsSelectedTPC(AliVTrack *trk);
  Bool_t IsSelectedKinematics(AliVTrack *trk);

  void CalculateMatchingEfficiency();

private:

  Int_t  fEventNumber;   // Event number in ESD file
  Bool_t fIsAODEvent;
  Bool_t fSameEvent;     // Check if it is still the same event in the loop to not fill the event vars for every track

  Int_t  fNObjs;       // Number of obj to be filled
  Int_t  fNHistos;       // Number of histos to be filled
  Int_t  fNSparses;       // Number of sparses to be filled

  Int_t fNtracksITS;    // to be removed only for crosscheck
  Int_t fNtracksTPC;    // to be removed only for crosscheck

  THashList *fHistoList;
  // TH3F *f3DHisoObjs[AliDielectronEvtVsTrkHist::k3DHistoNMaxValues]; // Activate if a 3D histogram is defined
  THnSparseD *fSparseObjs[AliDielectronEvtVsTrkHist::kSparseNMaxValues];

  Int_t fSparseVars[AliDielectronEvtVsTrkHist::kSparseNMaxValues][20];
  Bool_t fSparseIsTrackVar[AliDielectronEvtVsTrkHist::kSparseNMaxValues][20];
  Float_t fEventPlaneAngle[2]; // 0 = TPC, 1 = V0C

  AliPIDResponse *fPIDResponse; // Used for the pid cut

  THnSparseD *fMatchEffITS;  // Sparse object to normalize the matching eff
  THnSparseD *fMatchEffTPC;  // Sparse object to normalize the matching eff

  ClassDef(AliDielectronEvtVsTrkHist,1);
};

#endif
