#ifndef AliAnalysisTaskEtaPhigg_cxx
#define AliAnalysisTaskEtaPhigg_cxx

// Class for analysis of photon ggHBT correlations
// Authors: D.Peresunko

class THashList;
class TH3F;
class AliPHOSGeometry;
class AliCaloPhoton;
class AliAODTrack;
class AliAODEvent;
class AliEPFlattener;
class AliV0ReaderV1;
class AliConvEventCuts;
class AliConversionPhotonCuts;
class AliEMCALGeometry;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEtaPhigg : public AliAnalysisTaskSE
{
 public:
  enum CutList { kDefault, kDisp, kCPV, kBoth, kDistance1, kDistance2, kDistance3 };

  AliAnalysisTaskEtaPhigg(const char* name = "AliAnalysisTaskEtaPhigg");
  virtual ~AliAnalysisTaskEtaPhigg() {}

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

 protected:
  AliAnalysisTaskEtaPhigg(const AliAnalysisTaskEtaPhigg&);            // not implemented
  AliAnalysisTaskEtaPhigg& operator=(const AliAnalysisTaskEtaPhigg&); // not implemented

  void FillHistogram(const char* key, Double_t x) const;                         // Fill 1D histogram witn name key
  void FillHistogram(const char* key, Double_t x, Double_t y) const;             // Fill 2D histogram witn name key
  void FillHistogram(const char* key, Double_t x, Double_t y, Double_t z) const; // Fill 3D histogram witn name key
  void FillHistogram(const char* key, Double_t x, Double_t y, Double_t z,
                     Double_t w) const; // Fill 3D histogram witn name key

  Int_t ConvertRunNumber(Int_t run);
  Bool_t PairCut(const AliCaloPhoton* ph1, const AliCaloPhoton* ph2, Int_t cut) const;
  Bool_t IsGoodChannel(Int_t ix, Int_t iz); // for CPV in Mod3
  Int_t FindTrackMatching(int det, TVector3& locPHOS, int mod, float& dxMin, float& dzMin);
  int TestCPV(int mod, double e, double xPHOS, double zPHOS, double dxPHOS, double dzPHOS, int itr);
  void SetCutNames();

 protected:
  static constexpr int kCuts = 190;    // PID cuts
  static constexpr int kCentBins = 7;  // centrality
  static constexpr int kVtxBins = 20;  // z-vertex
  static constexpr int kPRBins = 6;    // Reaction plane
  static constexpr int kKtbins = 4;    // Kt bins for OSL param
  static constexpr double kgMass = 0.; // photon mass
  char fcut[kCuts][30];
  double fKtBins[kKtbins + 1] = { 0.2, 0.4, 0.7, 1., 2. };

  THashList* fOutputContainer;                      //! final histogram container
  AliAODEvent* fEvent;                              //!
  TList* fPHOSEvents[kVtxBins][kCentBins][kPRBins]; //! Containers for events with PHOS photons [zvtx][fCenBin][irp]
  TClonesArray* fPHOSEvent;                         //! PHOS photons in current event
  TClonesArray* fCPVEvent;                          //! CPV event

  // Reaction plain for v2
  AliEPFlattener* fV0AFlat; //!
  AliEPFlattener* fV0CFlat; //!
  Double_t fRP;             //!
  int fMF;

  Int_t fRunNumber;    // Current run number
  Float_t fCentrality; //! Centrality of the currecnt event
  Int_t fCenBin;       //! Current centrality bin

  AliPHOSGeometry* fPHOSGeo; //! PHOS geometry
  Int_t fEventCounter;       // number of analyzed events

  TH2C* fBadMap; //!

  TH2F* fhReQinv[kCentBins][kCuts];    //!
  TH2F* fhMiQinv[kCentBins][kCuts];    //!
  TH2F* fhReQinvCut[kCentBins][kCuts]; //!
  TH2F* fhMiQinvCut[kCentBins][kCuts]; //!
  TH2F* fhReq[kCentBins][kCuts];       //!
  TH2F* fhMiq[kCentBins][kCuts];       //!
  TH2F* fhReqCut[kCentBins][kCuts];    //!
  TH2F* fhMiqCut[kCentBins][kCuts];    //!

  TH3F* fhReOSL[kCentBins][kKtbins];    //!
  TH3F* fhMiOSL[kCentBins][kKtbins];    //!
  TH3F* fhReOSLCTS[kCentBins][kKtbins]; //!
  TH3F* fhMiOSLCTS[kCentBins][kKtbins]; //!

  ClassDef(AliAnalysisTaskEtaPhigg, 2); // PHOS analysis task
};

#endif
