#include "AliAnalysisMultPbTrackHistoManager.h"
#include "AliLog.h"
#include "TH1.h"
#include "TH3D.h"
#include "TH1I.h"
#include "TROOT.h"
#include "TMCProcess.h"
#include "AliMCParticle.h"
#include "TMap.h"
#include <iostream>
#include "TH2D.h"
#include "TDatabasePDG.h"
#include "TParameter.h"

using namespace std;

ClassImp(AliAnalysisMultPbTrackHistoManager)

const char * AliAnalysisMultPbTrackHistoManager::kStatStepNames[]     = { "All Events", "After centrality selection",  "After physics Selection", "With Vertex (quality cuts)", "After ZDC cut", "After vertex range cut (10 cm)" };
const char * AliAnalysisMultPbTrackHistoManager::kHistoPtEtaVzNames[] = { "hGenPtEtaVz", "hRecPtEtaVz", "hRecPtEtaVzPrim", 
									  "hRecPtEtaVzSecWeak", "hRecPtEtaVzSecMaterial", "hRecPtEtaVzFake"};
const char * AliAnalysisMultPbTrackHistoManager::kHistoDCANames[]     = { "hGenDCA", "hRecDCA", "hRecDCAPrim", "hRecDCASecWeak","hRecDCASecMaterial", "hRecDCAFake"};
const char * AliAnalysisMultPbTrackHistoManager::kHistoPrefix[]     = { "hGen", "hRec", "hRecPrim", "hRecSecWeak","hRecSecMaterial", "hRecFake", "hRecHighestMeanPt", "hRecLowestMeanPt"};
const char * AliAnalysisMultPbTrackHistoManager::kSpeciesName[]     = { "pi+", "K+", "p", "l+",  "pi-", "K-", "barp", "l-", "Others"};


AliAnalysisMultPbTrackHistoManager::AliAnalysisMultPbTrackHistoManager() : AliHistoListWrapper(), fHNameSuffix(""), fParticleSpecies(0){ 
  // standard ctor

}

AliAnalysisMultPbTrackHistoManager::AliAnalysisMultPbTrackHistoManager(const char * name, const char * title): AliHistoListWrapper(name,title), fHNameSuffix(""), fParticleSpecies(0)  {
  //named ctor
};

AliAnalysisMultPbTrackHistoManager::AliAnalysisMultPbTrackHistoManager(const AliAnalysisMultPbTrackHistoManager& obj) : AliHistoListWrapper (obj) {
  // copy ctor
  AliError("Not Implemented");
};

AliAnalysisMultPbTrackHistoManager::~AliAnalysisMultPbTrackHistoManager() {
  // dtor

}

TH2D * AliAnalysisMultPbTrackHistoManager::GetHistoElectronCutQA() {
  // Get a p vs dE/dx plot, to check the histos we are rejecting
  TString name = "histoElectronCutQA";

  TH2D * h = (TH2D*) GetHisto(name);
  
  if (!h) {
    AliInfo(Form("Booking %s", (name+fHNameSuffix).Data()));    
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    h = new TH2D (name+fHNameSuffix, name+fHNameSuffix, 100, 0, 3, 50, 0, 200);
    fList->Add(h);
    TH1::AddDirectory(oldStatus);

  }

  return h;
}


TH3D * AliAnalysisMultPbTrackHistoManager::GetHistoPtEtaVz(Histo_t id, Int_t particle) {
  // Returns a 3D histo of Pt/eta/vtx. It it does not exist, books it.

  TString name = kHistoPtEtaVzNames[id];

  if(particle >= 0) name += kSpeciesName[particle];

  TH3D * h = (TH3D*) GetHisto(name.Data());
  if (!h) {
    h = BookHistoPtEtaVz(name.Data(), Form("Pt Eta Vz distribution (%s)",kHistoPtEtaVzNames[id]));
  }

  return h;

}

TH2D * AliAnalysisMultPbTrackHistoManager::GetHistoDCA(Histo_t id) {
  // Returns a 3D histo of Pt/eta/vtx. It it does not exist, books it.

  TH2D * h = (TH2D*) GetHisto(kHistoDCANames[id]);
  if (!h) {
    h = BookHistoDCA(kHistoDCANames[id], Form("DCA vs pt distribution (%s)",kHistoDCANames[id]));
  }

  return h;

}

TH2D * AliAnalysisMultPbTrackHistoManager::GetHistoV0vsNtracks(Histo_t id) {
  // Returns a histo of V0 vs ntracks
  

  TString name = TString(kHistoPrefix[id])+"_V0vsNtracks";

  TH2D * h =  (TH2D*) GetHisto(name);
  if (!h) {
    name+=fHNameSuffix;
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    AliInfo(Form("Booking histo %s",name.Data()));

    h = new TH2D (name.Data(), Form("V0 vs Ntracks (%s)",kHistoPrefix[id]), 300,0,3000, 300, 0, 20000);			 
    h->SetXTitle("N_{Tracks}");
    h->SetYTitle("V0 Amplitude");
    TH1::AddDirectory(oldStatus);
    fList->Add(h);


  }
  return h;
  

}


TH1D * AliAnalysisMultPbTrackHistoManager::GetHistoMult(Histo_t id) {
  // Returns a 3D histo of Pt/eta/vtx. It it does not exist, books it.

  TString name = kHistoPrefix[id];
  name += "Mult";
  TH1D * h = (TH1D*) GetHisto(name.Data());
  if (!h) {
    h = BookHistoMult(name.Data(), Form("Multiplicity distribution (%s)",kHistoPrefix[id]));
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

  Int_t min1 = minEta  < -11111 ? -1 : h3D ->GetYaxis()->FindBin(minEta);
  Int_t min2  = minVz  < -11111 ? -1 : h3D ->GetZaxis()->FindBin(minVz) ;

  Int_t max1 = maxEta  < -11111 ? h3D->GetNbinsY() : h3D ->GetYaxis()->FindBin(maxEta-0.00001);
  Int_t max2  = maxVz  < -11111 ? h3D->GetNbinsZ() : h3D ->GetZaxis()->FindBin(maxVz -0.00001);
  // Int_t max1 = maxEta  < -11111 ? -1 : h3D ->GetYaxis()->FindBin(maxEta-0.00001);
  // Int_t max2  = maxVz  < -11111 ? -1 : h3D ->GetZaxis()->FindBin(maxVz -0.00001);


  TString hname = h3D->GetName();
  hname = hname +  "_pt_" + long (min1)  + "_" + long(max1) + "_" + long (min2)  + "_" + long(max2);

  
  if (gROOT->FindObjectAny(hname.Data())){
    AliError(Form("An object called %s already exists,adding suffix",hname.Data()));
    hname += "_2";
  }

  TH1D * h = h3D->ProjectionX(hname.Data(), min1, max1, min2, max2, "E");
  if(scaleWidth) h ->Scale(1.,"width");

  return h;

}


TH2D * AliAnalysisMultPbTrackHistoManager::GetHistoPtVz (Histo_t id, 
							 Float_t minEta, Float_t maxEta, 
							 Bool_t scaleWidth) {
  // Returns a projection of the 3D histo pt/eta/vz.
  // WARNING: since that is a histo, the requested range will be discretized to the binning.
  // Always avoids under (over) flows
  // If scaleWidth = kTRUE, the projection is scaled for the bin width (default)

  // FIXME: what do I do here for the scaling?

  TH3D * h3D = GetHistoPtEtaVz(id);

  // Get range in terms of bin numners.  If the float range is
  // less than -11111 take the range from the first to the last bin (i.e. no
  // under/over-flows)
  Int_t min1 = minEta  < -11111 ? -1 : h3D ->GetYaxis()->FindBin(minEta);
  Int_t max1 = maxEta  < -11111 ? -1 : h3D ->GetYaxis()->FindBin(maxEta-0.00001);


  TString hname = h3D->GetName();
  hname = hname +  "_ptvz_" + long (min1)  + "_" + long(max1);
  
  if (gROOT->FindObjectAny(hname.Data())){
    AliError(Form("An object called %s already exists,adding suffix",hname.Data()));
    hname += "_2";
  }

  h3D->GetYaxis()->SetRange(min1,max1);
  
  TH2D * h =  (TH2D*) h3D->Project3D("zxe");
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
    AliError(Form("An object called %s already exists, adding suffix",hname.Data()));
    hname+="_2";
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
    AliError(Form("An object called %s already exists, adding suffix",hname.Data()));
    hname+="_2";
  }

  TH1D * h = h3D->ProjectionY(hname.Data(), min1, max1, min2, max2, "E");
  if(scaleWidth) h ->Scale(1.,"width");
  return h;
}

TH1D * AliAnalysisMultPbTrackHistoManager::GetHistoProcess(Histo_t id) {

  // Returns histogram with particle specties

  TString name = TString(kHistoPrefix[id])+"_Process";

  TH1D * h =  (TH1D*) GetHisto(name);
  if (!h) {
    name+=fHNameSuffix;
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    AliInfo(Form("Booking histo %s",name.Data()));

    h = new TH1D (name.Data(), Form("Particle production process (%s)",kHistoPrefix[id]), kPNoProcess+1, -0.5, kPNoProcess+1-0.5);			 
    Int_t nbin = kPNoProcess+1;
    for(Int_t ibin = 0; ibin < nbin; ibin++){
      h->GetXaxis()->SetBinLabel(ibin+1,TMCProcessName[ibin]);      
    }
    TH1::AddDirectory(oldStatus);
    fList->Add(h);


  }
  return h;
  

}


TH1D * AliAnalysisMultPbTrackHistoManager::GetHistoSpecies(Histo_t id) {

  // Returns histogram with particle specties

  TString name = TString(kHistoPrefix[id])+"_Species";

  TH1D * h =  (TH1D*) GetHisto(name);
  if (!h) {
    name+=fHNameSuffix;
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    AliInfo(Form("Booking histo %s",name.Data()));

    h = new TH1D (name.Data(), Form("Particle species (%s)",kHistoPrefix[id]), kNPart+1, -0.5, kNPart+1-0.5);			 
    Int_t nbin = kNPart+1;
    for(Int_t ibin = 0; ibin < nbin; ibin++){
      h->GetXaxis()->SetBinLabel(ibin+1,kSpeciesName[ibin]);      
    }
    TH1::AddDirectory(oldStatus);
    fList->Add(h);


  }
  return h;
  

}


TH1D * AliAnalysisMultPbTrackHistoManager::GetHistoMeanPt(Histo_t id) {
  // mean pt computed event by event
  TString name = TString(kHistoPrefix[id])+"_MeanPt";

  TH1D * h =  (TH1D*) GetHisto(name);
  if (!h) {
    name+=fHNameSuffix;
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    AliInfo(Form("Booking histo %s",name.Data()));

    h = new TH1D (name.Data(), Form("Pt Event by Event(%s)",kHistoPrefix[id]), 100, 0, 1);			 
    h->Sumw2();
    TH1::AddDirectory(oldStatus);
    fList->Add(h);


  }
  return h;

}

TH1D * AliAnalysisMultPbTrackHistoManager::GetHistoPtEvent(Histo_t id) {

  // Returns of dNdpt, used to compute mean pt event by event

  TString name = TString(kHistoPrefix[id])+"_PtEvent";

  TH1D * h =  (TH1D*) GetHisto(name);
  if (!h) {
    name+=fHNameSuffix;
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    AliInfo(Form("Booking histo %s",name.Data()));
    const Int_t nptbins = 68;
    const Double_t binsPt[] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,40.0,45.0,50.0};

    h = new TH1D (name.Data(), Form("Pt Event by Event(%s)",kHistoPrefix[id]), nptbins,binsPt);			 
    h->SetBit(kKeepMaxMean,1);
    if(id == kHistoRecLowestMeanPt ) h->SetBit(kKeepMinMean,1);
    h->Sumw2();
    TH1::AddDirectory(oldStatus);
    fList->Add(h);


  }
  return h;
  

}


TH1D * AliAnalysisMultPbTrackHistoManager::GetHistoVzEvent(Histo_t id) {

  // Returns histogram with Vz of the event

  TString name = TString(kHistoPrefix[id])+"_VzEvent";

  TH1D * h =  (TH1D*) GetHisto(name);
  if (!h) {
    name+=fHNameSuffix;
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    AliInfo(Form("Booking histo %s",name.Data()));

    h = new TH1D (name.Data(), Form("Vz of the event (%s)",kHistoPrefix[id]), 10, -10, 10);			 
    h->Sumw2();
    h->SetXTitle("V_{z}^{event} (cm)");
    TH1::AddDirectory(oldStatus);
    fList->Add(h);


  }
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
    if (h->InheritsFrom("TH1I")) {
      AliInfo (Form("Not scaling integer histo %s",h->GetName()));
      continue;
    }
    AliInfo(Form("Scaling %s, nev %2.2f", h->GetName(), nev));
    
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
  
  // Binning from Jacek task
  // const Int_t nptbins = 49;
  // const Double_t binsPt[] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0};//,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,40.0,45.0,50.0};
  const Int_t nptbins = 68;
  const Double_t binsPt[] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,40.0,45.0,50.0};
  
  const Int_t netabins=20;
  Double_t binsEta[netabins+1];
  Float_t minEta = -1;
  Float_t maxEta =  1;
  Float_t etaStep = (maxEta-minEta)/netabins;
  for(Int_t ibin = 0; ibin < netabins; ibin++){    
    binsEta[ibin]   = minEta + ibin*etaStep;
    binsEta[ibin+1] = minEta + ibin*etaStep + etaStep;
  }

  const Int_t nvzbins=10;
  Double_t binsVz[nvzbins+1];
  Float_t minVz = -10;
  Float_t maxVz =  10;
  Float_t vzStep = (maxVz-minVz)/nvzbins;
  for(Int_t ibin = 0; ibin < nvzbins; ibin++){    
    binsVz[ibin]   = minVz + ibin*vzStep;
    binsVz[ibin+1] = minVz + ibin*vzStep + vzStep;
  }
 

  TH3D * h = new TH3D (hname,title, 
		       nptbins,  binsPt,
		       netabins, binsEta,
		       nvzbins,  binsVz
		       );

  h->SetXTitle("p_{T} (GeV)");
  h->SetYTitle("#eta");
  h->SetZTitle("V_{z}^{tracks} (cm)");
  h->Sumw2();
  
  TH1::AddDirectory(oldStatus);
  fList->Add(h);
  return h;
}

TH2D * AliAnalysisMultPbTrackHistoManager::BookHistoDCA(const char * name, const char * title) {
  // Books a DCA histo 

  const Int_t nptbins = 68;
  const Double_t binsPt[] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,40.0,45.0,50.0};
  const Int_t ndcabins=100;
  Double_t binsDca[ndcabins+1];
  Float_t minDca = -3;
  Float_t maxDca =  3;
  //  const Double_t binsDCA[] = {-3,-2.8,-2.6,-2.4,-2.2,-2.0,-1.9,};
  Float_t dcaStep = (maxDca-minDca)/ndcabins;
  for(Int_t ibin = 0; ibin < ndcabins; ibin++){    
    binsDca[ibin]   = minDca + ibin*dcaStep;
    binsDca[ibin+1] = minDca + ibin*dcaStep + dcaStep;
  }
 



  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TString hname = name;
  hname+=fHNameSuffix;

  AliInfo(Form("Booking %s",hname.Data()));
  

#if defined WEIGHTED_DCA
  TH1D * h = new TH1D (hname,title, 200,0,200);
  h->SetXTitle("#Delta DCA");
#elif defined TRANSVERSE_DCA
  TH2D * h = new TH2D (hname,title, ndcabins,binsDca,nptbins,binsPt);
  h->SetYTitle("p_{T} (GeV)");
  h->SetXTitle("d_{0} r#phi (cm)");
#endif
  h->Sumw2();
  
  fList->Add(h);

  TH1::AddDirectory(oldStatus);
  return h;
}
TH1D * AliAnalysisMultPbTrackHistoManager::BookHistoMult(const char * name, const char * title) {
  // Books a multiplicity histo 

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TString hname = name;
  hname+=fHNameSuffix;

  AliInfo(Form("Booking %s",hname.Data()));
  

  TH1D * h = new TH1D (hname,title, 600, 0,6000);

  h->SetXTitle("N tracks");
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

Int_t AliAnalysisMultPbTrackHistoManager::GetLocalParticleID(AliMCParticle * part) {
  // returns the local code (Part_t)
  
  Int_t pdgcode = part->PdgCode();
  switch(pdgcode) {

  case 211:
    return kPartPiPlus;
    break;
  case -211:
    return kPartPiMinus;
    break;
  case 2212:
    return kPartP;
    break;
  case -2212:
    return kPartPBar;
    break;
  case 321:
    return kPartKPlus;
    break;
  case -321:
    return kPartKMinus;
    break;
  case -11:
    return kPartLMinus;
    break;
  case 11:
    return kPartLPlus;
    break;
  case -13:
    return kPartLMinus;
    break;
  case 13:
    return kPartLPlus;
    break;
  default:
    return kPartOther;
  }
}


void AliAnalysisMultPbTrackHistoManager::FillSpeciesMap(Int_t pdgCode) {

  // Fills map of species vs ID
  if(!fParticleSpecies) fParticleSpecies = new TMap;
  
  TObjString * key = new TObjString(Form("%s",TDatabasePDG::Instance()->GetParticle(pdgCode)->GetName()));
  TParameter<Long64_t>* count = dynamic_cast<TParameter<Long64_t>*> (fParticleSpecies->GetValue(key));
  if (!count)
  {
    count = new TParameter<Long64_t>(key->String().Data(), 0);
    fParticleSpecies->Add(key, count);
  }
  count->SetVal(count->GetVal() + 1);

}

 
Long64_t AliAnalysisMultPbTrackHistoManager::Merge(TCollection* list)
{
  // Merge a list of AliHistoListWrapper objects with this.
  // Returns the number of merged objects (including this).

  // We have to make sure that all the list contain the same histos in
  // the same order. We thus also have to sort the list (sorting is
  // done by name in TList).

  AliInfo("Merging");
  if (!list)

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  TList collections;

  Int_t count = 0;

  while ((obj = iter->Next())) {
    cout << "1" << endl;
    //    AliHistoListWrapper* entry = dynamic_cast<AliHistoListWrapper*> (obj);
    AliAnalysisMultPbTrackHistoManager* entry = dynamic_cast<AliAnalysisMultPbTrackHistoManager*> (obj);
    cout << "2 " << endl;
    if (entry == 0) 
      continue;
    cout << "3 " << entry->GetName() << endl;

    // Merge map fParticleSpecies        
    if(!entry->fParticleSpecies){
      AliInfo("Cannot get fParticleSpecies");
      continue;
    }
    TIterator* iter2 = entry->fParticleSpecies->MakeIterator();
    TObjString* obj2 = 0;
    cout << "4" << endl;
    while ((obj2 = dynamic_cast<TObjString*> (iter2->Next())))
    {
      cout << "5" << endl;
      TParameter<Long64_t>* param2 = static_cast<TParameter<Long64_t>*> (entry->fParticleSpecies->GetValue(obj2));
      
      TParameter<Long64_t>* param1 = dynamic_cast<TParameter<Long64_t>*> (fParticleSpecies->GetValue(obj2));
      cout << " - (other)" << param2->GetName() << " " << param2->GetVal();
      if (param1)
      {
	cout << " (this) " << param1->GetName() << " " << param1->GetVal() << endl;
        param1->SetVal(param1->GetVal() + param2->GetVal());
      }
      else
      {
	cout << "" << endl;	
        param1 = dynamic_cast<TParameter<Long64_t>*> (param2->Clone());
        fParticleSpecies->Add(new TObjString(obj2->String()), param1);
      }
    }    
    delete iter2;

  }

  // merge histograms
  // We need a separate loop as this one could be restarted if needed.
  iter->Reset();
  while ((obj = iter->Next())) {
    AliAnalysisMultPbTrackHistoManager* entry = dynamic_cast<AliAnalysisMultPbTrackHistoManager*> (obj);
    if (entry == 0) 
      continue;

    Bool_t foundDiffinThisIterStep = kFALSE;

    //    Printf("%d - %s",count, obj->GetName());

    TList * hlist = entry->GetList();

    // Check if all histos in this fList are also in the one from entry and viceversa
    // Use getters to automatically book non defined histos    

    Bool_t areListsDifferent=kTRUE;
    Int_t iloop = 0;
    Int_t max_loops = hlist->GetSize() + fList->GetSize(); // In the worst case all of the histos will be different...
    while(areListsDifferent) {
      if(iloop>max_loops) AliFatal("Infinite Loop?");
      iloop++;
      // sort
      hlist->Sort();
      fList->Sort();
      // loop over the largest 
      TObject * hist =0;
      TIterator * iterlist = 0;
      TList * thislist  = 0; // the list over which I'm iterating
      TList * otherlist = 0; // the other

      if (hlist->GetSize() >= fList->GetSize()) { 
	thislist  = hlist;
	otherlist = fList;
      }
      else{
	thislist  = fList;
	otherlist = hlist;	
      }
      iterlist = thislist->MakeIterator();

      while ((hist= iterlist->Next())){ 
	if(!otherlist->FindObject(hist->GetName())){
	  AliInfo(Form("Adding object %s",hist->GetName()));	  
	  TH1 * hclone =  (TH1*) hist->Clone();
	  if (!hclone->InheritsFrom("TH1")) AliFatal(Form("Found a %s. This class only supports objects inheriting from TH1",hclone->ClassName()));
	  hclone->Reset();
	  otherlist->Add(hclone);
	  foundDiffinThisIterStep=kTRUE;
	}
      }

      // re-sort before checking
      hlist->Sort();
      fList->Sort();

      // check if everything is fine    
      areListsDifferent=kFALSE;
      if (hlist->GetSize() == fList->GetSize()) {	
	Int_t nhist =  fList->GetSize();
	for(Int_t ihist = 0; ihist < nhist; ihist++){
	  if(strcmp(fList->At(ihist)->GetName(),hlist->At(ihist)->GetName())) areListsDifferent = kTRUE;
	}
      } else {
	areListsDifferent=kTRUE;
      }
    }

    // last check: if something is not ok die loudly 
    if (hlist->GetSize() != fList->GetSize()) {
      AliFatal("Mismatching size!");
    }
    Int_t nhist =  fList->GetSize();
    for(Int_t ihist = 0; ihist < nhist; ihist++){
      if(strcmp(fList->At(ihist)->GetName(),hlist->At(ihist)->GetName())){
	AliFatal(Form("Mismatching histos: %s -> %s", fList->At(ihist)->GetName(),hlist->At(ihist)->GetName()));
      } else {
	// Specific merging strategies, according to the bits...
	if (fList->At(ihist)->TestBits(kKeepMinMean|kKeepMaxMean)){
	  TH1 * h1 = (TH1*) fList->At(ihist);
	  TH1 * h2 = (TH1*) hlist->At(ihist);
	  if (h1->GetEntries()>0 && h2->GetEntries()>0) {
	    if (h1->TestBit(kKeepMinMean)) {
	      AliInfo(Form("Keeping only the histo with the lowest mean [%s][%f][%f]",h1->GetName(), h1->GetMean(), h2->GetMean()));
	      if(h1->GetMean() > h2->GetMean()) h1->Reset();
	      else h2->Reset();
	    }
	    if (h1->TestBit(kKeepMaxMean)) {
	      AliInfo(Form("Keeping only the histo with the highest mean [%s][%f][%f]",h1->GetName(), h1->GetMean(), h2->GetMean()));
	      if(h2->GetMean() > h1->GetMean()) h1->Reset();
	      else h2->Reset();
	    }
	  }
	}
      }
    }
    
    if (foundDiffinThisIterStep){
      iter->Reset(); // We found a difference: previous lists could
		     // also be affected... We start from scratch
      collections.Clear();
      count = 0;
    }
    else {
      
      collections.Add(hlist);
      
      count++;
    }
  }

  fList->Merge(&collections);
  
  delete iter;

  AliInfo("Merged");

  return count+1;
}
