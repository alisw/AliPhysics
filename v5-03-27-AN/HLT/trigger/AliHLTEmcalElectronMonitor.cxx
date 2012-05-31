#include "AliHLTEmcalElectronMonitor.h"
#include "AliHLTScalars.h"
#include "TString.h"
#include "TMath.h"


ClassImp(AliHLTEmcalElectronMonitor);

AliHLTEmcalElectronMonitor::AliHLTEmcalElectronMonitor():
  hList(NULL),
  hTracksPt(NULL),
  hClusterEn(NULL),
  hdEta(NULL),
  hdPhi(NULL),
  hdR(NULL),
  hEoverP(NULL)
{

  // book histograms
  hList = new TObjArray;

  hTracksPt   = new TH1F("hTracksPt","Tracks pT (GeV/c)", 500, 0, 100);
  hList->Add(hTracksPt);

  hClusterEn  = new TH1F("hClusterEn","Cluster Energy (GeV)", 500, 0, 100);
  hList->Add(hClusterEn);

  hdEta = new TH1F("hdEta", "#Delta #eta (Cluster-Track)", 200,-.1,.1);
  hList->Add(hdEta);

  hdPhi = new TH1F("hdPhi", "#Delta #phi (Cluster-Track)", 200,-.1,.1);
  hList->Add(hdPhi);
  
  hdR = new TH1F("hdR","#Delta R (Track-Cluster);#Delta R(#sqrt{#Delta #eta ^{2} +#Delta #Phi ^{2}})",200,0.,.1);
  hList->Add(hdR);
  
  hEoverP=new TH1F("hEoverP","E/P for matched tracks;E/P",200,0.,10.);
  hList->Add(hEoverP);
  
}
//___________________________________________________________________________________________________________________________________________________

AliHLTEmcalElectronMonitor::~AliHLTEmcalElectronMonitor()
{

  // default destructor

}
//___________________________________________________________________________________________________________________________________________________

TObjArray* AliHLTEmcalElectronMonitor::GetHistograms()
{

  // pointer to histogram objects
  
  return hList;

}
//___________________________________________________________________________________________________________________________________________________

Int_t AliHLTEmcalElectronMonitor::MakeHisto(AliHLTScalars *scalar)
{

  // make the histograms
  
  hTracksPt->Fill( scalar->GetScalar("TracksPt").Value() );
  
  hClusterEn->Fill( scalar->GetScalar("ClusterEn").Value() );
  
  hdEta->Fill( scalar->GetScalar("dEta").Value() );
  
  hdPhi->Fill( scalar->GetScalar("dPhi").Value() );
  
  hdR->Fill( scalar->GetScalar("dR").Value() );
  
  hEoverP->Fill( scalar->GetScalar("EoverP").Value() );
  

  return 0;

}
  
