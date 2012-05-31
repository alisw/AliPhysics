#ifndef AliMFTClusterQA_H
#define AliMFTClusterQA_H

#include "TObject.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliLoader.h"
#include "AliMFT.h"
#include "TClonesArray.h"
#include "AliMFTCluster.h"
#include "AliMFTSegmentation.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "AliLog.h"
#include "TString.h"

//====================================================================================================================================================
//
// Class for the analysis of the MFT clusters (a.k.a. rec points). Few QA histograms are created
//
// Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

class AliMFTClusterQA : public TObject {
  
public:

  AliMFTClusterQA();
  virtual ~AliMFTClusterQA() {;}

  void Init(Char_t *readDir, Char_t *outDir, Int_t nEventsToAnalyze);
  Bool_t LoadNextEvent();
  void BookHistos();
  void Terminate();

private:

  AliMFTClusterQA(const AliMFTClusterQA& obj);
  AliMFTClusterQA& operator=(const AliMFTClusterQA& other);

protected:

  static const Int_t fNMaxPlanes = AliMFTConstants::fNMaxPlanes;

  TH1D *fHistNClustersPerEvent[fNMaxPlanes], *fHistNPixelsPerCluster[fNMaxPlanes];
  TH1D *fHistClusterSizeX[fNMaxPlanes], *fHistClusterSizeY[fNMaxPlanes], *fHistClusterRadialPosition[fNMaxPlanes];
  TH2D *fClusterScatterPlotXY[fNMaxPlanes];

  AliLoader *fMFTLoader;
  AliRunLoader *fRunLoader;
  AliMFT *fMFT;

  Int_t fNPlanes, fNEvents, fEv;

  TFile *fFileOut;

  TString fReadDir, fOutDir;

  ClassDef(AliMFTClusterQA, 1); 

};

//====================================================================================================================================================
 
#endif


