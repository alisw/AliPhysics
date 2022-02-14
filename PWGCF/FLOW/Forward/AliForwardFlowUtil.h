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
#include "AliFMDStripIndex.h"
#include "AliFMDEncodedEdx.h"
#include "AliFMDMCTrackDensity.h"

class AliForwardFlowUtil : public TObject {
  typedef std::vector< Double_t > edgeContainer;

 public:
  AliForwardFlowUtil();
  // ~AliForwardFlowUtil();                                       // destructor
  // AliForwardFlowUtil(const AliForwardFlowUtil &L);             // copy constructor
  // AliForwardFlowUtil & operator=(const AliForwardFlowUtil &L); // assignment
  
  Int_t GetNUARunNumber(Int_t runnumber);
  Bool_t IsGoodRun(Int_t runnumber);
  Bool_t XeXe_Run(Int_t runnumber);
  Bool_t PbPb_lowIR_Run(Int_t runnumber);
  Bool_t PbPb_highIR_Run(Int_t runnumber);
  Bool_t pPb_Run(Int_t runnumber);
  Bool_t ExtraEventCutFMD(TH2D& forwarddNdedp, double cent, Bool_t mc,TH2D* hOutliers);
  void FillData(TH2D*& refDist, TH2D*& centralDist, TH2D*& forwardDist);
  void FillDataCentral(TH2D*& centralDist);

  // ESD
  void FillFromTrackrefsITS(TH2D*& fwd) ;
  void FillFromTrackrefsFMD(TH2D*& fwd) ;
  void FillFromTrackrefsFMDperTR(TH2D*& fwd) ;
  void FillFromPrimariesFMD(TH2D*& fwd) const;
  void FillFromPrimariesFMDperTR(TH2D*& fwd) ;
  void FillFromPrimariesTPC(TH2D*& cen) const;
  void FillFromPrimariesSPD(TH2D*& cen) const;
  void FillFromPrimariesITS(TH2D*& cen) const;

  // AOD
  void FillFromForwardClusters(TH2D*& fwd);
  void FillFromCentralClusters(TH2D*& cen) const;
  void FillFromTracklets(TH2D*& cen) const;
  void FillFromTracks(TH2D*& cen, UInt_t tracktype) const;
  void FillFromPrimariesAODITS(TH2D*& cen) const;
  void FillFromPrimariesAODSPD(TH2D*& cen) const;
  void FillFromPrimariesAODFMD(TH2D*& fwd) const;
  void FillFromPrimariesAODTPC(TH2D*& cen) const;

  // unused
  void FillFromPrimaries(TH3D*& cen, TH3D*& fwd, Double_t zvertex) const;
  void FillFromTracklets(TH3D*& cen, Double_t zvertex) const;
  void FillFromTracks(TH3D*& cen, Int_t tracktype, Double_t zvertex) const;
  void FillFromPrimaries(TH2D*& cen, TH2D*& fwd) const;
  void FillFromPrimaries(TH2D*& cen) const;
  void FillFromPrimariesAOD(TH2D*& cen, TH2D*& fwd) const;
  void FillFromPrimariesAOD(TH2D*& cen) const;
  AliMCParticle* GetMother(AliMCParticle* p);
  Bool_t IsRedefinedPhysicalPrimary(AliMCParticle* p);

  Bool_t ProcessTrackITS(AliMCParticle* particle,TH2D*& cen);

  Double_t GetZ();
  Double_t GetCentrality(TString centrality_estimator);
  // Check if a given particle itself hit the FMD. If so, return the
    // (first) track reference of such a hit
  AliTrackReference* IsHitTPC(AliMCParticle* p);
  AliTrackReference* IsHitFMD(AliMCParticle* p);
  Double_t GetTrackRefEta(AliMCParticle* p);
  Double_t GetTrackRefEta(AliTrackReference* ref);
  Double_t GetTrackRefPhi(AliMCParticle* p);
  Double_t GetTrackRefPhi(AliTrackReference* ref);

  void MakeFakeHoles(TH2D& forwarddNdedp);
  Bool_t FMDAcceptanceExistMC(Double_t eta,Double_t phi,Double_t vertex);

  AliVEvent* fevent; //!
  AliAODEvent* fAODevent; //!
  AliMCEvent* fMCevent; //!
  TH1F* dNdeta; //!
  AliForwardSettings fSettings;
  Double_t minpt;//!
  Double_t maxpt;//!
  Bool_t dodNdeta;//!
  AliFMDMCTrackDensity* fTrackDensity; //!

  mutable struct State
  {
    Double_t angle;            // Angle
    UShort_t oldDetector;      // Last detector
    Char_t   oldRing;          // Last ring
    UShort_t oldSector;        // Last sector
    UShort_t oldStrip;         // Last strip
    UShort_t startStrip;       // First strip
    UShort_t nRefs;            // Number of references
    UShort_t nStrips;          // Number of strips
    UShort_t count;            // Count of hit strips
    AliTrackReference* longest; //! Longest track through
    /**
     * Clear this state
     *
     * @param alsoCount If true, also clear count
    */
    void Clear(Bool_t alsoCount=false);
    /**
     * Assignment operator
     *
     * @param o Object to assign from
     *
     * @return Reference to this object
     */
    State& operator=(const State& o);
  } fState; //! State
  Double_t fLowCutvalue;
  Bool_t            fTrackGammaToPi0;
  AliTrackReference*  ProcessRef(AliMCParticle* particle,AliTrackReference* ref,TH2D*& fwd);

  void BeginTrackRefs();
  void EndTrackRefs();

  void StoreParticle(AliMCParticle* particle,AliTrackReference* ref,TH2D*& fwd) ;


  Bool_t ProcessTrack(AliMCParticle* particle, TH2D*& fwd);

  Double_t GetTrackRefTheta(const AliTrackReference* ref) const;

  AliTrackReference* fStored; //! Last stored



// For a new value newValue, compute the new count, new mean, the new M2.
// mean accumulates the mean of the entire dataset
// M2 aggregates the squared distance from the mean
// count aggregates the number of samples seen so far



private:
  ClassDef(AliForwardFlowUtil, 1);
};
#endif
