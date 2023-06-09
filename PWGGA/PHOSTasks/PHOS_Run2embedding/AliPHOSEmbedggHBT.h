#ifndef AliPHOSEmbedggHBT_cxx
#define AliPHOSEmbedggHBT_cxx

// Class for analysis of photon dEta-dPhi correlations in embedded data
// Authors: D.Peresunko

class THashList;
class AliPHOSGeometry;
class AliCaloPhoton;
class AliAODTrack;
class AliAODEvent;
class AliEPFlattener;
class AliV0ReaderV1;
class AliConvEventCuts;
class AliConversionPhotonCuts;
class AliEMCALGeometry;

#include "AliAnalysisTaskEtaPhigg.h"

class AliPHOSEmbedggHBT : public AliAnalysisTaskEtaPhigg
{
 public:
  AliPHOSEmbedggHBT() = default;
  AliPHOSEmbedggHBT(const char* name);
  virtual ~AliPHOSEmbedggHBT();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*) {}

 protected:
  AliPHOSEmbedggHBT(const AliPHOSEmbedggHBT&);            // not implemented
  AliPHOSEmbedggHBT& operator=(const AliPHOSEmbedggHBT&); // not implemented

  Bool_t SecondaryPi0Cut(const AliCaloPhoton* ph1, const AliCaloPhoton* ph2) const;
  Int_t JetRejection(Int_t module) const;
  Bool_t IsGoodChannel(Int_t ix, Int_t iz); // for CPV in Mod3
  void ReclusterizeCPV();
  int CommonParent(const AliCaloPhoton* p1, const AliCaloPhoton* p2) const;

 protected:
  TList* fMCEvents = nullptr;           //! PHOS photons in current event
  TClonesArray* fSignalEvent = nullptr; //! PHOS photons in current event
  TList* fSignalEvents = nullptr;       //!

  TH2F* fhReQinvMCprim = nullptr;       //!
  TH2F* fhReQinvMC = nullptr;           //!
  TH2F* fhMiQinvMCprim = nullptr;       //!
  TH2F* fhMiQinvMC = nullptr;           //!
  TH2F* fhReqMCprim = nullptr;          //!
  TH2F* fhReqMC = nullptr;              //!
  TH2F* fhMiqMCprim = nullptr;          //!
  TH2F* fhMiqMC = nullptr;              //!
  TH2F* fhReQinvConv[kCentBins][kCuts]; //!
  TH2F* fhReqConv[kCentBins][kCuts];    //!
  TH2F* fhReQinvSignal[kCuts];          //!
  TH2F* fhReQinvSignalConv[kCuts];      //!
  TH2F* fhMiQinvSignal[kCuts];          //!
  TH2F* fhReqSignal[kCuts];             //!
  TH2F* fhReqSignalConv[kCuts];         //!
  TH2F* fhMiqSignal[kCuts];             //!

  ClassDef(AliPHOSEmbedggHBT, 1);
};

#endif
