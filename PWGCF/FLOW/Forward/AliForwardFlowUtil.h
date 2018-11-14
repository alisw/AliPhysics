#ifndef AliForwardFlowUtil_cxx
#define AliForwardFlowUtil_cxx
/**
 * @file AliForwardFlowUtil.h
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief
 *
 * @ingroup pwgcf_forward_flow
 */
#include "TObject.h"
#include "TString.h"
#include "TH3.h"
#include "TH2.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliTrackReference.h"
#include "AliMCParticle.h"
#include "AliAODVertex.h"

class AliForwardFlowUtil : public TObject {
  typedef std::vector< Double_t > edgeContainer;

 public:
   AliForwardFlowUtil();

  Bool_t ExtraEventCutFMD(TH2D& forwarddNdedp, double cent, Bool_t mc);
  void FillFromTrackrefs(TH2D*& cen, TH2D*& fwd) const;
  void FillFromPrimaries(TH2D*& cen, TH2D*& fwd) const;
  void FillFromTracklets(TH2D*& cen) const;
  void FillFromTracks(TH2D*& cen, Int_t tracktype) const;
  Double_t GetZ();
  Double_t GetCentrality(TString centrality_estimator);
  // Check if a given particle itself hit the FMD. If so, return the
    // (first) track reference of such a hit
  AliTrackReference* IsHitTPC(AliMCParticle* p);
  AliTrackReference* IsHitFMD(AliMCParticle* p);
  AliVEvent* fevent; //!
  AliAODEvent* fAODevent; //!
  AliMCEvent* fMCevent; //!

private:
  ClassDef(AliForwardFlowUtil, 2);
};
#endif
