//____________________________________________________________________
//
// $Id$
//
// Script that contains a class to get the media where a certain track
// was created.   This is used for background studies. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <TH2D.h>
#include <AliFMDHit.h>
#include <AliFMDInput.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <TArrayI.h>
#include <AliRun.h>
#include <TNtuple.h>
#include <AliDetector.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMedium.h>
#include <TParticle.h>
#include <THStack.h>
/** @class Media 
    @brief Media description 
    @ingroup FMD_script
 */
struct Media : public TNamed
{
  TArrayI         fMeds;
  TNtuple*        fCount;
  TH1D*           fHist;
  Media(const char* name, const char* title) 
    : TNamed(name, title), fMeds(0) 
  {
    fHist  = new TH1D(Form("h%s",name), title, 120, -6, 6);
    // fHist->SetFillStyle(3001);
    fCount = new TNtuple(GetName(), GetTitle(), 
			 "E:DET:RNG:SEC:STR:X:Y:Z:EDEP:PDG");
  }
  Media(const char* name, TArrayI* a) 
    : TNamed(name, "Media information"), fMeds(a->fN, a->fArray)
  {
    fHist  = new TH1D(Form("h%s",name), GetTitle(), 120, -6, 6);
    fHist->SetFillStyle(3001);
    fHist->SetFillColor(fMeds.fN > 0 ? fMeds[0] : 0);
    fCount = new TNtuple(GetName(), GetTitle(), 
			 "E:DET:RNG:SEC:STR:X:Y:Z:EDEP:PDG");
    std::cout << "Media list " << name << ":" << std::endl;
    for (Int_t i = 0; i < fMeds.fN; i++) {
      std::cout << " " << fMeds[i] << std::flush;
      if (fMeds[i] == 0 && i > 0) break;
    }
    std::cout << std::endl;
  }
  Bool_t HasMedia(Int_t id) 
  {
    if (fMeds.fN <= 0) return kFALSE;
    for (Int_t i = 0; i < fMeds.fN; i++) {
      if (fMeds[i] == id) return kTRUE;
      if (fMeds[i] == 0 && i > 0) break;
    }
    return kFALSE;
  }
  void Fill(Int_t ev, AliFMDHit* hit) 
  {
    Float_t x[10];
    x[0] = ev;
    x[1] = hit->Detector();
    x[2] = Int_t(hit->Ring());
    x[3] = hit->Sector();
    x[4] = hit->Strip();
    x[5] = hit->X();
    x[6] = hit->Y();
    x[7] = hit->Z();
    x[8] = hit->Edep();
    x[9] = hit->Pdg();
    // return;
    fCount->Fill(x);    
    // r = TMath::Sqrt(hit->X() * hit->X() + hit->Y() * hit->Y());
    // if(r!=0) Float_t wgt=1./(2. * 3.1415*r);
    Double_t r = TMath::Sqrt(hit->X()*hit->X()+hit->Y()*hit->Y());
    Double_t t = TMath::ATan2(r, hit->Z());
    Double_t e = -TMath::Log(TMath::Tan(t/2));
    fHist->Fill(e);
  }
};


//____________________________________________________________________
/** @class GetMedia
    @brief Get media where a particle is produced
    @code 
    aliroot
    Root> gAlice->Init("$ALICE_ROOT/FMD/Config.C");
    Root> .x $ALICE_ROOT/FMD/scriptsWriteMedArrays.C
    Root> gAlice->Run();
    Root> .q
    aliroot
    Root> .L Compile.C
    Root> Compile("GetMedia.C")
    Root> GetMedia c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class GetMedia : public AliFMDInput
{
private:
  TString    fModList;
  TObjArray  fMedia;
  Media*     fOther;
  Media*     fAll;
  Media*     fPrim;
  Int_t      fEv;
  THStack*   fStack;
  TFile*     fOutput;
public:
  //__________________________________________________________________
  GetMedia(const char* modlist="FMD:ITS:BODY:ABSO:T0:PIPE", 
	   const char* output="media.root") 
    :  fModList(modlist)
  { 
    AddLoad(kGeometry);
    AddLoad(kHits);
    AddLoad(kKinematics);

    fOutput = TFile::Open(output, "RECREATE");
    fStack  = new THStack("summed", "Summed hits");
    fOther  = AddMedia("other", "Unknown media");
    fAll    = AddMedia("All", "All media");
    fPrim   = AddMedia("Primary", "Primary particles");
    fEv     = 0;

    
    fAll->fHist->SetFillColor(0);
    fAll->fHist->SetFillStyle(3004);
    fPrim->fHist->SetFillColor(1);
    fPrim->fHist->SetFillStyle(3001);
    fStack->Add(fPrim->fHist);
  }
  //__________________________________________________________________
  Media* AddMedia(const char* name, const char* title, TArrayI* a=0)
  {
    Media* media = 0;
    if (!a) media = new Media(name, title);
    else    media = new Media(name, a);
    if (a)  {
      fMedia.Add(media);
      fStack->Add(media->fHist);
    }
    return media;
  }
  //__________________________________________________________________
  Media* FindMedia(Int_t med) 
  {
    for (Int_t i = 0; i < fMedia.GetEntries(); i++) {
      Media* media = static_cast<Media*>(fMedia.At(i));
      if (!media) continue;
      if (media->HasMedia(med)) return media;
    }
    return 0;
  }
  //__________________________________________________________________
  Bool_t Init()
  {
    Bool_t ret = AliFMDInputHits::Init();
    if (!gGeoManager) {
      Error("Init", "No geometry defined - make sure you have that");
      return kFALSE;
    }
    if (!fRun) {
      Error("Init", "gAlice not defined");
      return kFALSE;
    }
    Int_t now = 0;
    Int_t colon;
    TFile* file = TFile::Open("medid.root", "READ");
    if (!file) {
      Error("Init", "couldn't open medid.root");
      return kFALSE;
    }
    fOutput->cd();
    while ((colon = fModList.Index(":", now)) >= 0) {
      TString d(fModList(now, colon-now));
      now = colon+1;
      TArrayI* a = 0;
      file->GetObject(d.Data(), a);
      if (!a) {
	Warning("Init", "No medium array for %s", d.Data());
	continue;
      }
      AddMedia(d.Data(),0, a);
    }
    if (now < fModList.Length()) {
      colon = fModList.Length();
      TString d(fModList(now, colon-now));
      TArrayI* a = 0;
      file->GetObject(d.Data(), a);
      if (!a) 
	Warning("Init", "No medium array for %s", d.Data());
      else 
	AddMedia(d.Data(),0, a);
    }
    file->Close();
    if (fMedia.GetEntries() <= 0) return kFALSE;
    return ret;
  }  
  //__________________________________________________________________
  /** Begining of event
      @param ev Event number
      @return @c false on error */
  Bool_t Begin(Int_t ev) 
  {
    fEv = ev;
    return AliFMDInputHits::Begin(ev);
  }
  //__________________________________________________________________
  Bool_t ProcessHit(AliFMDHit* hit, TParticle* track) 
  {
    if (!hit || !track) {
      std::cout << "No hit or track (hit: " << hit << ", track: " 
		<< track << ")" << std::endl;
      return kFALSE;
    }
    fAll->Fill(fEv, hit);
    if (track->IsPrimary()) { 
      fPrim->Fill(fEv, hit);
      return kTRUE;
    }
      
    // Get production vertex 
    Double_t vx = track->Vx();
    Double_t vy = track->Vy();
    Double_t vz = track->Vz();
    // Get node 
    TGeoNode* prodNode = gGeoManager->FindNode(vx,vy,vz);
    if (!prodNode) { 
      Warning("ProcessHit", "didn't find a node for production vertex");
      return kTRUE;
    }
    //Info("ProcessHit", "Production vertex in node %s", prodNode->GetName());

    // Get volume 
    TGeoVolume* prodVol = prodNode->GetVolume();
    if (!prodVol) { 
      Warning("ProcessHit", "didn't find a volume for production vertex");
      return kTRUE;
    }
    //Info("ProcessHit", "Production vertex in volume %s",prodVol->GetName());
    
    // Med medium 
    TGeoMedium* prodMed = prodVol->GetMedium();
    if (!prodMed) { 
      Warning("ProcessHit", "didn't find a medium for production vertex");
      return kTRUE;
    }
    // Info("ProcessHit", "Production vertex in medium %s %d", 
    //      prodMed->GetName(), prodMed->GetId());
    
    
    Media* media = FindMedia(prodMed->GetId());
    if (!media) { 
      Warning("ProcessHit", "Media not found %s (%d)", 
	      prodMed->GetName(), prodMed->GetId());
      media = fOther; 
    }
    // Info("ProcessHit", "Adding to %s", media->GetName());
    media->Fill(fEv, hit);
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t Finish()
  {
    fOutput->cd();
    fStack->Write();
    fStack->Draw();
    fOutput->Write();
    fOutput->Close();
    fOutput = 0;
    fMedia.Delete();
    return kTRUE;
  }
};

//____________________________________________________________________
//
// EOF
//
