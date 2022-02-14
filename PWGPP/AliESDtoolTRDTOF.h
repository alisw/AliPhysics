#ifndef ALIESDTOOLTRDTOF_H
#define ALIESDTOOLTRDTOF_H

class AliPIDResponse;
class TTreeSRedirector;
class TTreeStream;
class TTree;
class TGraph;
class TH1F;
class AliExternalTrackParam;
class AliESDEvent;
class AliESDfriend;
class AliESDtrack;
class TTree;
//class TVectorD;
#include "TNamed.h"
#include "TVectorD.h"
#include "AliTRDgeometry.h"

class AliESDtoolTRDTOF : public TNamed {
public:
  AliESDtoolTRDTOF();
  void     CacheTRDGeom();
  void     MakeActiveMapFromV0(TTree * treeV0, Int_t nPoints=10000, Double_t threshold=0.3);
  Double_t ChacheTrackInfo();
  //
  Double_t GetDet(Int_t iLayer){ return fVecDet[iLayer%7];}
  Double_t GetSec(Int_t iLayer){ return fVecSec[iLayer%7];}
  Double_t GetZ(Int_t iLayer){ return fVecZ[iLayer%7];}
  Double_t GetdSec(Int_t iLayer){ return fVecdSec[iLayer%7];}
  Double_t GetdEdge(Int_t iLayer){ return fVecdEdge[iLayer%7];}
  Double_t GetStatus(Int_t iLayer){ return fVecStatus[iLayer%7];}
  Double_t isActive(Int_t iLayer){ return fVecActive[iLayer%7];}
  Double_t isNotActive(Int_t iLayer) {return (fVecDeadR[iLayer%7]+fVecDeadZ[iLayer%7]+fVecDeadDet[iLayer%7])>0;}
  Double_t isDeadR(Int_t iLayer){ return fVecDeadR[iLayer%7];}
  Double_t isDeadZ(Int_t iLayer){ return fVecDeadZ[iLayer%7];}
  Double_t isDeadDet(Int_t iLayer){ return fVecDeadDet[iLayer%7];}
  Double_t TRDFound(){return nFound;}
  Double_t TRDFindable(){return nFindable;}
  public:
  AliESDtrack *fESDtrack;
  TTree * fTree;
  const Int_t kNDetectors=640; // number of TRD+TOF detectors
  const Int_t occuCut=8;       // occupancy cut
  const Double_t kMarginR=5;
  const Double_t kMarginZ=3;
  //
  TVectorD fActiveMap;
  TH1F* fHistoDetector;
  TVectorD fVecSec, fVecZ, fVecDet, fVecdSec, fVecdEdge, fVecActive, fVecStatus,fVecDeadZ,fVecDeadR, fVecDeadDet0,fVecDeadDet;
  Double_t nFindable;
  Double_t nFound;
  AliTRDgeometry geom;
  TMatrixD fZBoundary0;
  TMatrixD fZBoundary1;
  //
  TTreeSRedirector * fStreamer;                  /// streamer
  static AliESDtoolTRDTOF* fgInstance;                /// instance of the tool -needed in order to use static functions (for TTreeFormula)

  private:
  AliESDtoolTRDTOF(AliESDtoolTRDTOF&);
  AliESDtoolTRDTOF &operator=(const AliESDtoolTRDTOF&);

};

#endif