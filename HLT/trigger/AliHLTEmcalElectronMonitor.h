#ifndef ALIHLTEMCALELECTRONMONITOR_H
#define ALIHLTEMCALELECTRONMONITOR_H

#include "TH1F.h"
#include "TObjArray.h"
#include "TString.h"
#include "AliHLTScalars.h"

class AliHLTEmcalElectronMonitor : public TObject
{

 public:
  // constructor
  AliHLTEmcalElectronMonitor();

  // destructor
  virtual ~AliHLTEmcalElectronMonitor();

  // make histos
  Int_t MakeHisto(AliHLTScalars *scalar);

  // retrieve histograms
  TObjArray* GetHistograms();

 private:
  TObjArray *hList;
  TH1F      *hTracksPt;
  TH1F      *hClusterEn;
  TH1F      *hdEta;
  TH1F      *hdPhi;
  TH1F      *hdR;
  TH1F      *hEoverP;
  
  AliHLTEmcalElectronMonitor(const AliHLTEmcalElectronMonitor &);
  AliHLTEmcalElectronMonitor & operator = (const AliHLTEmcalElectronMonitor &);

  ClassDef(AliHLTEmcalElectronMonitor, 0);

};

#endif


