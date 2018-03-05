#ifndef ALIFEMTOCutValues_H
#define ALIFEMTOCutValues_H
//
//author O.Arnold

#include <iostream>
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TBits.h"
#include "TObject.h"
#include "TVector2.h"
#include "AliESDtrack.h"
#include "TVector3.h"
#include "TDatabasePDG.h"

using namespace std;


class AliFemtoCutValues
{
 private:

  //Event:
  Float_t fEvtZvtxLow; //lower z-vertex selection
  Float_t fEvtZvtxUp; //upper z-vertex selection

  //V0s
  Float_t fV0PIDThresholdPtTPCLow; //At which pt-value the PID for V0s starts
  Float_t fV0PIDThresholdPtTPCUp; //At which pt-value the PID for V0s stops
  Float_t fV0SelWidth;
  Float_t fV0EtaRange;
  Float_t fAntiCutK0sLow;
  Float_t fAntiCutK0sUp;
  Float_t fV0Nsigma;
  Float_t fV0Pointing;
  Float_t fV0RxyLow;
  Float_t fV0RxyUp;
  Float_t fV0Decayvtx;
  Float_t fV0DCAtrackPV;
  Float_t fV0DCAtracksV0decay;
  Int_t fV0PtBinsCPA; //Number of pt bins to do cosine pointing angle differentially
  Int_t fV0PtBinsInvMass; //Number of pt bins to do invariant mass differentially
  Int_t fV0TPCCluster;

  //Protons
  Float_t fProtonPIDThresholdPtTPCLow; //At which pt-value the PID for Protons starts
  Float_t fProtonPIDThresholdPtTOFUp; //At which pt-value the PID for Protons stops
  Float_t fProtonPIDTPCTOFSwitch; //Value at which additionally also the TOF is used for PID
  Float_t fProtonDCAxyCut;
  Float_t fProtonDCAzCut;
  Float_t fProtonEtaRange;
  Float_t fProtonNsigma;
  Int_t fProtonPtBinsDCA; //Number of pt bins to do DCA plots differentially
  Int_t fProtonPtBinsPurity; //Number of pt bins to do Purity plots differentially
  Int_t fProtonTPCCluster;
  Int_t fTrackFilterBit;

  //functions
  Int_t fFindPtBin(Double_t ptVal,Double_t ptLow,Double_t ptUp,Int_t nBins);


  TDatabasePDG fPDGDatabase; //!

 public:

  enum systematics
  {
    kDefault,
    kProtonVariationLowerPtThresholdUp,
    kProtonVariationLowerPtThresholdDown,
    kProtonVariationEtaRangeUp,
    kProtonVariationEtaRangeDown,
    kProtonVariationNsigmaUp,
    kProtonVariationNsigmaDown,
    kProtonTPCClusterUp,
    kProtonVariationFilterBitGlobal,
    kV0VariationLowerPtThresholdDown,
    kV0VariationLowerPtThresholdUp,
    kV0VariationCosinePointingUp,
    kV0VariationNsigmaDown,
    kV0VariationTPCClusterUp,
    kV0VariationEtaRangeUp,
    kV0VariationEtaRangeDown,
    kV0VariationDCAatV0DecayDown,
    kV0VariationDCADaughtersToPVUp
  };

  AliFemtoCutValues();
  AliFemtoCutValues(systematics cutType);
  virtual ~AliFemtoCutValues();
  AliFemtoCutValues(const AliFemtoCutValues &obj);
  AliFemtoCutValues &operator=(const AliFemtoCutValues &obj);

  //Int_t fFindPtBin(Double_t ptVal,Double_t ptLow,Double_t ptUp,Int_t nBins);

  //Event
  Double_t GetEvtzVertexLow() {return fEvtZvtxLow;}
  Double_t GetEvtzVertexUp() {return fEvtZvtxUp;}

  //V0s
  Int_t FindV0PtBinInvMass(Double_t V0pt){return fFindPtBin(V0pt,fV0PIDThresholdPtTPCLow,fV0PIDThresholdPtTPCUp,fV0PtBinsInvMass);}
  Int_t FindV0PtBinCPA(Double_t V0pt){return fFindPtBin(V0pt,fV0PIDThresholdPtTPCLow,fV0PIDThresholdPtTPCUp,fV0PtBinsCPA);}

  //Setter
  void SetCutVariation(systematics cutType);

  //Getter
  Int_t GetV0PtBinsCPA(){return fV0PtBinsCPA;}
  Int_t GetV0PtBinsInvMass(){return fV0PtBinsInvMass;}
  Int_t GetV0TPCCluster(){return fV0TPCCluster;}
  Int_t GetFilterBit(){return fTrackFilterBit;}

  Double_t GetV0PIDWidth(){return fV0SelWidth;}
  Double_t GetV0PIDthrPtLow(){return fV0PIDThresholdPtTPCLow;}
  Double_t GetV0PIDthrPtUp(){return fV0PIDThresholdPtTPCUp;}
  Double_t GetHadronMasses(Int_t Pdgcode) {return (fPDGDatabase.GetParticle(Pdgcode))->Mass();}
  Double_t GetV0EtaRange(){return fV0EtaRange;}
  Double_t GetK0sAntiCutLow(){return fAntiCutK0sLow;}
  Double_t GetK0sAntiCutUp(){return fAntiCutK0sUp;}
  Double_t GetV0nsigma(){return fV0Nsigma;}
  Double_t GetV0pointing(){return fV0Pointing;}
  Double_t GetV0rxylow(){return fV0RxyLow;}
  Double_t GetV0rxyup(){return fV0RxyUp;}
  Double_t GetV0decayvtxcut(){return fV0Decayvtx;}
  Double_t GetV0DCAtrPV(){return fV0DCAtrackPV;}
  Double_t GetV0DCAtracksV0decay(){return fV0DCAtracksV0decay;}

  //Protons
  Int_t FindProtonPtBinDCA(Double_t Protonpt){return fFindPtBin(Protonpt,fProtonPIDThresholdPtTPCLow,fProtonPIDThresholdPtTOFUp,fProtonPtBinsDCA);}
  Int_t FindProtonPtBinPurity(Double_t Protonpt){return fFindPtBin(Protonpt,fProtonPIDTPCTOFSwitch,fProtonPIDThresholdPtTOFUp,fProtonPtBinsPurity);}

  Int_t GetProtonTPCCluster(){return fProtonTPCCluster;}

  Double_t GetProtonPIDthrPtLow(){return fProtonPIDThresholdPtTPCLow;}
  Double_t GetProtonPIDthrPtUp(){return fProtonPIDThresholdPtTOFUp;}
  Double_t GetProtonPIDTOFTPCSwitch(){return fProtonPIDTPCTOFSwitch;}
  Double_t GetProtonDCAxyCut(){return fProtonDCAxyCut;}
  Double_t GetProtonDCAzCut(){return fProtonDCAzCut;}
  Double_t GetProtonEtaRange(){return fProtonEtaRange;}
  Double_t GetProtonnsigma(){return fProtonNsigma;}

  Int_t GetPtBinsPurity(){return fProtonPtBinsPurity;}
  Int_t GetPtBinsDCA(){return fProtonPtBinsDCA;}

  ClassDef(AliFemtoCutValues, 2)
};


#endif
