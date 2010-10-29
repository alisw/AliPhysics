#include "AliAnalysisMultPbTrackHistoManager.h"
#include "AliLog.h"
#include "TH1.h"
#include "TH3D.h"
#include "TH1I.h"
#include "TROOT.h"
#include <iostream>
using namespace std;
ClassImp(AliAnalysisMultPbTrackHistoManager)

const char * AliAnalysisMultPbTrackHistoManager::kStatStepNames[]     = { "All Events", "After centrality selection", "With Vertex" };
const char * AliAnalysisMultPbTrackHistoManager::kHistoPtEtaVzNames[] = { "hGenPtEtaVz", "hRecPtEtaVz", "hRecPtEtaVzPrim", 
									  "hRecPtEtaVzSecWeak", "hRecPtEtaVzSecMaterial", "hRecPtEtaVzFake"};
const char * AliAnalysisMultPbTrackHistoManager::kHistoDCANames[]     = { "hGenDCA", "hRecDCA", "hRecDCAPrim", "hRecDCASecWeak","hRecDCASecMaterial", "hRecDCAFake"};



AliAnalysisMultPbTrackHistoManager::AliAnalysisMultPbTrackHistoManager() : AliHistoListWrapper(), fHNameSuffix(""){ 
  // standard ctor

}

AliAnalysisMultPbTrackHistoManager::AliAnalysisMultPbTrackHistoManager(const char * name, const char * title): AliHistoListWrapper(name,title), fHNameSuffix("")  {
  //named ctor
};

AliAnalysisMultPbTrackHistoManager::AliAnalysisMultPbTrackHistoManager(const AliAnalysisMultPbTrackHistoManager& obj) : AliHistoListWrapper (obj) {
  // copy ctor
  AliError("Not Implemented");
};

AliAnalysisMultPbTrackHistoManager::~AliAnalysisMultPbTrackHistoManager() {
  // dtor

}

TH3D * AliAnalysisMultPbTrackHistoManager::GetHistoPtEtaVz(Histo_t id) {
  // Returns a 3D histo of Pt/eta/vtx. It it does not exist, books it.

  TH3D * h = (TH3D*) GetHisto(kHistoPtEtaVzNames[id]);
  if (!h) {
    h = BookHistoPtEtaVz(kHistoPtEtaVzNames[id], Form("Pt Eta Vz distribution (%s)",kHistoPtEtaVzNames[id]));
  }

  return h;

}

TH1D * AliAnalysisMultPbTrackHistoManager::GetHistoDCA(Histo_t id) {
  // Returns a 3D histo of Pt/eta/vtx. It it does not exist, books it.

  TH1D * h = (TH1D*) GetHisto(kHistoDCANames[id]);
  if (!h) {
    h = BookHistoDCA(kHistoDCANames[id], Form("Pt Eta Vz distribution (%s)",kHistoDCANames[id]));
  }

  return h;

}

TH1D * AliAnalysisMultPbTrackHistoManager::GetHistoPt (Histo_t id, 
						       Float_t minEta, Float_t maxEta, 
						       Float_t minVz, Float_t maxVz, 
						       Bool_t scaleWidth) {
  // Returns a projection of the 3D histo pt/eta/vz.
  // WARNING: since that is a histo, the requested range will be discretized to the binning.
  // Always avoids under (over) flows
  // If scaleWidth = kTRUE, the projection is scaled for the bin width (default)

  TH3D * h3D = GetHistoPtEtaVz(id);

  // Get range in terms of bin numners.  If the float range is
  // less than -11111 take the range from the first to the last bin (i.e. no
  // under/over-flows)
  Int_t min1 = minEta  < -11111 ? 1 : h3D ->GetYaxis()->FindBin(minEta);
  Int_t min2  = minVz  < -11111 ? 1 : h3D ->GetZaxis()->FindBin(minVz) ;

  Int_t max1 = maxEta  < -11111 ? h3D->GetNbinsY() : h3D ->GetYaxis()->FindBin(maxEta-0.00001);
  Int_t max2  = maxVz  < -11111 ? h3D->GetNbinsZ() : h3D ->GetZaxis()->FindBin(maxVz -0.00001);


  TString hname = h3D->GetName();
  hname = hname +  "_pt_" + long (min1)  + "_" + long(max1) + "_" + long (min2)  + "_" + long(max2);

  
  if (gROOT->FindObjectAny(hname.Data())){
    AliError(Form("An object called %s already exists",hname.Data()));
  }

  TH1D * h = h3D->ProjectionX(hname.Data(), min1, max1, min2, max2, "E");
  if(scaleWidth) h ->Scale(1.,"width");

  return h;

}

TH1D * AliAnalysisMultPbTrackHistoManager::GetHistoVz (Histo_t id, 
						       Float_t minPt, Float_t maxPt,
						       Float_t minEta, Float_t maxEta,
						       Bool_t scaleWidth) { 
  // Returns a projection of the 3D histo pt/eta/vz.
  // WARNING: since that is a histo, the requested range will be discretized to the binning.
  // Always avoids under (over) flows
  // If scaleWidth = kTRUE, the projection is scaled for the bin width (default)

  TH3D * h3D = GetHistoPtEtaVz(id);

  // Get range in terms of bin numners.  If the float range is
  // less than -11111 take the range from the first to the last bin (i.e. no
  // under/over-flows)
  Int_t min1  = minPt  < -11111 ? 1 : h3D ->GetXaxis()->FindBin(minPt) ;
  Int_t min2  = minEta < -11111 ? 1 : h3D ->GetYaxis()->FindBin(minEta);

  Int_t max1  = maxPt  < -11111 ? h3D->GetNbinsX() : h3D ->GetXaxis()->FindBin(maxPt -0.00001);
  Int_t max2  = maxEta < -11111 ? h3D->GetNbinsY() : h3D ->GetYaxis()->FindBin(maxEta-0.00001);


  TString hname = h3D->GetName();
  hname = hname +  "_Vz_" + long (min1)  + "_" + long(max1) + "_" + long (min2)  + "_" + long(max2);

  if (gROOT->FindObjectAny(hname.Data())){
    AliError(Form("An object called %s already exists",hname.Data()));
  }

  TH1D * h = h3D->ProjectionZ(hname.Data(), min1, max1, min2, max2, "E");
  if(scaleWidth) h ->Scale(1.,"width");
  return h;


}

TH1D * AliAnalysisMultPbTrackHistoManager::GetHistoEta (Histo_t id, 
							Float_t minPt, Float_t maxPt, 
							Float_t minVz, Float_t maxVz,
							Bool_t scaleWidth) {
  // Returns a projection of the 3D histo pt/eta/vz.
  // WARNING: since that is a histo, the requested range will be discretized to the binning.
  // Always avoids under (over) flows
  // If scaleWidth = kTRUE, the projection is scaled for the bin width (default)

  TH3D * h3D = GetHistoPtEtaVz(id);

  // Get range in terms of bin numners.  If the float range is
  // less than -11111 take the range from the first to the last bin (i.e. no
  // under/over-flows)
  Int_t min1 = minPt < -11111 ? 1 : h3D ->GetXaxis()->FindBin(minPt) ;
  Int_t min2 = minVz < -11111 ? 1 : h3D ->GetYaxis()->FindBin(minVz);

  Int_t max1 = maxPt < -11111 ? h3D->GetNbinsX() : h3D ->GetXaxis()->FindBin(maxPt -0.00001);
  Int_t max2 = maxVz < -11111 ? h3D->GetNbinsY() : h3D ->GetYaxis()->FindBin(maxVz-0.00001);

  TString hname = h3D->GetName();
  hname = hname +  "_Eta_" + long (min1)  + "_" + long(max1) + "_" + long (min2)  + "_" + long(max2);

  if (gROOT->FindObjectAny(hname.Data())){
    AliError(Form("An object called %s already exists",hname.Data()));
  }

  TH1D * h = h3D->ProjectionY(hname.Data(), min1, max1, min2, max2, "E");
  if(scaleWidth) h ->Scale(1.,"width");
  return h;
}


TH1I * AliAnalysisMultPbTrackHistoManager::GetHistoStats() {
  // Returns histogram with event statistiscs (processed events at each step)

  TH1I * h =  (TH1I*) GetHisto("hStats");
  if (!h) h = BookHistoStats();
  return h;
  
} 



TH1 * AliAnalysisMultPbTrackHistoManager::GetHisto(const char * name) {
  //Search list for histo
  // TODO: keep track of histo id rather than searching by name?
  return (TH1*) fList->FindObject(TString(name)+fHNameSuffix);

}

void AliAnalysisMultPbTrackHistoManager::ScaleHistos(Double_t nev, Option_t * option) {
  // Scales all histos in the list for nev
  // option can be used to pass further options to TH1::Scale
  TH1 * h = 0;
  TIter iter = fList->MakeIterator();
  while ((h = (TH1*) iter.Next())) {
    if (!h->InheritsFrom("TH1")) {
      AliFatal (Form("%s does not inherits from TH1, cannot scale",h->GetName()));
    }
    cout << "Scaling " << h->GetName() << " " << nev << endl;
    
    h->Scale(1./nev,option);
  }

}

TH3D * AliAnalysisMultPbTrackHistoManager::BookHistoPtEtaVz(const char * name, const char * title) {
  // Books a 3D histo of Pt/eta/vtx
  // TODO: make the binning settable, variable binning?

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TString hname = name;
  hname+=fHNameSuffix;

  AliInfo(Form("Booking %s",hname.Data()));
  

  TH3D * h = new TH3D (hname,title, 
		       50,0,10, // Pt
		       20,-1, 1,   // Eta
		       10,-10,10 // Vz
		       );

  h->SetYTitle("#eta");
  h->SetXTitle("p_{T}");
  h->SetZTitle("V_{z} (cm)");
  h->Sumw2();
  
  fList->Add(h);

  TH1::AddDirectory(oldStatus);
  return h;
}

TH1D * AliAnalysisMultPbTrackHistoManager::BookHistoDCA(const char * name, const char * title) {
  // Books a DCA histo 

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TString hname = name;
  hname+=fHNameSuffix;

  AliInfo(Form("Booking %s",hname.Data()));
  

  TH1D * h = new TH1D (hname,title, 100,0,50);

  h->SetXTitle("#Delta DCA");
  h->Sumw2();
  
  fList->Add(h);

  TH1::AddDirectory(oldStatus);
  return h;
}

TH1I * AliAnalysisMultPbTrackHistoManager::BookHistoStats() {
  // Books histogram with event statistiscs (processed events at each step)

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  AliInfo(Form("Booking Stat histo"));

  TH1I * h = new TH1I (TString("hStats")+fHNameSuffix, "Number of processed events", kNStatBins, -0.5, kNStatBins-0.5);
  for(Int_t istatbin = 0; istatbin < kNStatBins; istatbin++){
    h->GetXaxis()->SetBinLabel(istatbin+1,kStatStepNames[istatbin]);
  }
  TH1::AddDirectory(oldStatus);
  fList->Add(h);
  return h;
}




 
