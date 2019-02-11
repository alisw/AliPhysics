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
#include "AliForwardSettings.h"

class AliForwardFlowUtil : public TObject {
  typedef std::vector< Double_t > edgeContainer;

 public:
   AliForwardFlowUtil();


  Bool_t ExtraEventCutFMD(TH2D& forwarddNdedp, double cent, Bool_t mc,TH2D* hOutliers);
  void FillData(TH2D*& refDist, TH2D*& centralDist, TH2D*& forwardDist);

  void FillFromTrackrefs(TH2D*& cen, TH2D*& fwd) const;
  void FillFromTrackrefs(TH2D*& fwd) const;
  void FillFromPrimaries(TH2D*& cen, TH2D*& fwd) const;
  void FillFromPrimariesAOD(TH2D*& cen, TH2D*& fwd) const;
  void FillFromPrimaries(TH2D*& cen) const;
  void FillFromPrimariesAOD(TH2D*& cen) const;
  void FillFromTracklets(TH2D*& cen) const;
  void FillFromTracks(TH2D*& cen, UInt_t tracktype) const;

  void FillFromTrackrefs(TH3D*& cen, TH3D*& fwd, Double_t zvertex) const;
  void FillFromPrimaries(TH3D*& cen, TH3D*& fwd, Double_t zvertex) const;
  void FillFromTracklets(TH3D*& cen, Double_t zvertex) const;
  void FillFromTracks(TH3D*& cen, Int_t tracktype, Double_t zvertex) const;
  Double_t GetZ();
  Double_t GetCentrality(TString centrality_estimator);
  // Check if a given particle itself hit the FMD. If so, return the
    // (first) track reference of such a hit
  AliTrackReference* IsHitTPC(AliMCParticle* p);
  AliTrackReference* IsHitFMD(AliMCParticle* p);

  void MakeFakeHoles(TH2D& forwarddNdedp);
  AliVEvent* fevent; //!
  AliAODEvent* fAODevent; //!
  AliMCEvent* fMCevent; //!
  Bool_t mc; //!
  TH1F* dNdeta;
  AliForwardSettings fSettings;
  Double_t maxpt;//!
  Double_t minpt;//!
  Bool_t dodNdeta;//!
private:
  ClassDef(AliForwardFlowUtil, 2);
};
#endif
