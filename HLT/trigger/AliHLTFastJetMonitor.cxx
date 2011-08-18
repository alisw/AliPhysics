#include "AliHLTFastJetMonitor.h"
#include "AliHLTScalars.h"
#include "TString.h"
#include "TMath.h"


ClassImp(AliHLTFastJetMonitor);

AliHLTFastJetMonitor::AliHLTFastJetMonitor():
  hList(NULL),
  hTracksPt(NULL),
  hClusterEn(NULL),
  hClusterEta(NULL),
  hClusterPhi(NULL),
  hJetsPt(NULL)
{

  // book histograms
  hList = new TObjArray;

  hTracksPt   = new TH1F("hTracksPt",   "Tracks pT (GeV/c)", 300, 0, 300);
  hList->Add(hTracksPt);

  hClusterEn  = new TH1F("hClusterEn",  "Cluster Energy (GeV)", 300, 0, 300);
  hList->Add(hClusterEn);

  hClusterEta = new TH1F("hClusterEta", "Cluster position in #eta", 300, -2*TMath::Pi(),2*TMath::Pi());
  hList->Add(hClusterEta);

  hClusterPhi = new TH1F("hClusterPhi", "Cluster position in #phi", 300, -TMath::Pi(),TMath::Pi());
  hList->Add(hClusterPhi);
}
//___________________________________________________________________________________________________________________________________________________

AliHLTFastJetMonitor::~AliHLTFastJetMonitor()
{

  // default destructor

}
//___________________________________________________________________________________________________________________________________________________

TObjArray* AliHLTFastJetMonitor::GetHistograms()
{

  // pointer to histogram objects
  
  return hList;

}
//___________________________________________________________________________________________________________________________________________________

Int_t AliHLTFastJetMonitor::MakeHisto(AliHLTScalars *scalar)
{

  // make the histograms
  
  hTracksPt->Fill( scalar->GetScalar("TracksPt").Value() );
  
  hClusterEn->Fill( scalar->GetScalar("ClusterEn").Value() );
  
  hClusterEta->Fill( scalar->GetScalar("ClusterEta").Value() );
  
  hClusterPhi->Fill( scalar->GetScalar("ClusterPhi").Value() );
  
  hJetsPt->Fill( scalar->GetScalar("JetsPt").Value() );

  return 0;

}
  
