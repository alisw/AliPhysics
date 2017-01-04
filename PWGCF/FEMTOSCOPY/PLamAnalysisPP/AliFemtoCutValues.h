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
  Double_t fEvtZvtxLow; //lower z-vertex selection
  Double_t fEvtZvtxUp; //upper z-vertex selection

  //V0s
  Double_t fV0PIDThresholdPtTPCLow; //At which pt-value the PID for V0s starts
  Double_t fV0PIDThresholdPtTPCUp; //At which pt-value the PID for V0s stops
  Double_t fV0SelWidth;
  Double_t fV0EtaRange;
  Double_t fAntiCutK0sLow;
  Double_t fAntiCutK0sUp;
  Double_t fV0Nsigma;
  Double_t fV0Pointing;
  Double_t fV0RxyLow;
  Double_t fV0RxyUp;
  Double_t fV0Decayvtx;
  Double_t fV0DCAtrackPV;
  Double_t fV0DCAtracksV0decay;
  Int_t fV0PtBinsCPA; //Number of pt bins to do cosine pointing angle differentially
  Int_t fV0PtBinsInvMass; //Number of pt bins to do invariant mass differentially
  Int_t fV0TPCCluster;
  
  //Protons
  Double_t fProtonPIDThresholdPtTPCLow; //At which pt-value the PID for Protons starts
  Double_t fProtonPIDThresholdPtTOFUp; //At which pt-value the PID for Protons stops
  Double_t fProtonPIDTPCTOFSwitch; //Value at which additionally also the TOF is used for PID
  Double_t fProtonDCAxyCut;
  Double_t fProtonDCAzCut;
  Double_t fProtonEtaRange;
  Double_t fProtonNsigma;
  Int_t fProtonPtBinsDCA; //Number of pt bins to do DCA plots differentially
  Int_t fProtonPtBinsPurity; //Number of pt bins to do Purity plots differentially
  Int_t fProtonTPCCluster;

  //functions
  Int_t fFindPtBin(Double_t ptVal,Double_t ptLow,Double_t ptUp,Int_t nBins);
  
  
  TDatabasePDG *fPDGDatabase; //!
  
 public:
  
  AliFemtoCutValues();
  virtual ~AliFemtoCutValues();
  AliFemtoCutValues(const AliFemtoCutValues &obj);
  AliFemtoCutValues &operator=(const AliFemtoCutValues &obj);
  
  //Int_t fFindPtBin(Double_t ptVal,Double_t ptLow,Double_t ptUp,Int_t nBins);

  //Event
  Double_t GetEvtzVertexLow() {return fEvtZvtxLow;};
  Double_t GetEvtzVertexUp() {return fEvtZvtxUp;};

  //V0s
  Int_t FindV0PtBinInvMass(Double_t V0pt){return fFindPtBin(V0pt,fV0PIDThresholdPtTPCLow,fV0PIDThresholdPtTPCUp,fV0PtBinsInvMass);};
  Int_t FindV0PtBinCPA(Double_t V0pt){return fFindPtBin(V0pt,fV0PIDThresholdPtTPCLow,fV0PIDThresholdPtTPCUp,fV0PtBinsCPA);};

  //Getter
  Int_t GetV0PtBinsCPA(){return fV0PtBinsCPA;}
  Int_t GetV0PtBinsInvMass(){return fV0PtBinsInvMass;}
  Int_t GetV0TPCCluster(){return fV0TPCCluster;};

  Double_t GetV0PID_width(){return fV0SelWidth;}
  Double_t GetV0PIDthrPtLow(){return fV0PIDThresholdPtTPCLow;};
  Double_t GetV0PIDthrPtUp(){return fV0PIDThresholdPtTPCUp;};
  Double_t GetHadronMasses(Int_t Pdgcode) {return fPDGDatabase->GetParticle(Pdgcode)->Mass();};
  Double_t GetV0EtaRange(){return fV0EtaRange;};
  Double_t GetK0sAntiCutLow(){return fAntiCutK0sLow;};
  Double_t GetK0sAntiCutUp(){return fAntiCutK0sUp;};
  Double_t GetV0nsigma(){return fV0Nsigma;};
  Double_t GetV0pointing(){return fV0Pointing;};
  Double_t GetV0rxylow(){return fV0RxyLow;};
  Double_t GetV0rxyup(){return fV0RxyUp;};
  Double_t GetV0decayvtxcut(){return fV0Decayvtx;};
  Double_t GetV0DCAtrPV(){return fV0DCAtrackPV;};
  Double_t GetV0DCAtracksV0decay(){return fV0DCAtracksV0decay;};
  
  //Protons
  Int_t FindProtonPtBinDCA(Double_t Protonpt){return fFindPtBin(Protonpt,fProtonPIDThresholdPtTPCLow,fProtonPIDThresholdPtTOFUp,fProtonPtBinsDCA);};
  Int_t FindProtonPtBinPurity(Double_t Protonpt){return fFindPtBin(Protonpt,fProtonPIDTPCTOFSwitch,fProtonPIDThresholdPtTOFUp,fProtonPtBinsPurity);};

  Int_t GetProtonTPCCluster(){return fProtonTPCCluster;};
  
  Double_t GetProtonPIDthrPtLow(){return fProtonPIDThresholdPtTPCLow;}
  Double_t GetProtonPIDthrPtUp(){return fProtonPIDThresholdPtTOFUp;}
  Double_t GetProtonPIDTOFTPCSwitch(){return fProtonPIDTPCTOFSwitch;}
  Double_t GetProtonDCAxyCut(){return fProtonDCAxyCut;}
  Double_t GetProtonDCAzCut(){return fProtonDCAzCut;}
  Double_t GetProtonEtaRange(){return fProtonEtaRange;}
  Double_t GetProtonnsigma(){return fProtonNsigma;};
  
  Int_t GetPtBinsPurity(){return fProtonPtBinsPurity;};
  Int_t GetPtBinsDCA(){return fProtonPtBinsDCA;};
  
  //ClassDef(AliFemtoCutValues, 1);
};


#endif
