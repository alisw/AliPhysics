#ifndef ALIHLTFASTJETMONITOR_H
#define ALIHLTFASTJETMONITOR_H

#include "TH1F.h"
#include "TObjArray.h"
#include "TString.h"
#include "AliHLTScalars.h"

class AliHLTFastJetMonitor : public TObject
{

 public:
  // constructor
  AliHLTFastJetMonitor();

  // destructor
  virtual ~AliHLTFastJetMonitor();

  // make histos
  Int_t MakeHisto(AliHLTScalars *scalar);

  // retrieve histograms
  TObjArray* GetHistograms();

 private:
  TObjArray *hList;
  TH1F      *hTracksPt;
  TH1F      *hClusterEn;
  TH1F      *hClusterEta;
  TH1F      *hClusterPhi;
  TH1F      *hJetsPt;

  AliHLTFastJetMonitor(const AliHLTFastJetMonitor &);
  AliHLTFastJetMonitor & operator = (const AliHLTFastJetMonitor &);

  ClassDef(AliHLTFastJetMonitor, 0);

};

#endif


