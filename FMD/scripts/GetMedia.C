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

/** @class Media 
    @brief Media description 
    @ingroup FMD_script
 */
struct Media : public TNamed
{
  TArrayI*        fMeds;
  TNtuple*        fCount;
  Media(const char* name, const char* title) 
    : TNamed(name, title), fMeds(0) 
  {
    fCount = new TNtuple(GetName(), "E/I:DET:RNG:SEC:STR:X/F:Y:Z:EDEP:PDG/I");
  }
  Media(const char* name) 
    : TNamed(name, "Media information"), fMeds(0), fCount(0)
  {
    AliDetector* det = gAlice->GetDetector(name);
    if (!det) {
      Warning("Media", "Detector %s not found in gAlice", name);
      return;
    }
    fMeds = det->GetIdtmed();
    if (!fMeds) {
      Warning("Media", "No mediums for detector %s found", name);
      return;
    }
    fCount = new TNtuple(GetName(), "E/I:DET:RNG:SEC:STR:X/F:Y:Z:EDEP:PDG/I");
  }
};


//____________________________________________________________________
/** @class GetMedia
    @brief Get media where a particle is produced
    @code 
    Root> .L Compile.C
    Root> Compile("GetMedia.C")
    Root> GetMedia c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class GetMedia : public AliFMDInputHits
{
private:
  TString    fModList;
  TObjArray  fMedia;
  Media*     fOther;
  Media*     fAll;
  Int_t      fEv;
  TFile*     fOutput;
public:
  //__________________________________________________________________
  GetMedia(const char* modlist="FMD:ITS:BODY:ABSO:START:PIPE", 
	   const char* output="media.root") 
  { 
    fOutput = TFile::Open(output, "RECREATE");
    fOther  = new Media("other", "Unknown media"),
    fAll    = new Media("All", "All media")
    fEv     = 0;
    AddLoad(kKinematics);
  }
  //__________________________________________________________________
  Media* FindMedia(Int_t med) 
  {
    TIter next(&fMedia);
    Media* media = 0;
    while ((media == static_cast<Media*>(next()))) {
      if (!media->fMeds) continue;
      TArrayI& ids = *(media->fMeds);
      for (Int_t i = 0; i < ids.fN; i++) 
	if (med == ids[i]) return media;
    }
    return 0;
  }
  //__________________________________________________________________
  Bool_t Init()
  {
    if (!gGeoManager) {
      Error("GetMedia", "No geometry defined - make sure you have that");
      return kFALSE;
    }
    if (!gAlice) {
      Error("GetMedia", "gAlice not defined");
      return kFALSE;
    }
    Int_t now = 0;
    Int_t colon;
    while ((colon = fModList.Index(":", now)) >= 0) {
      fMedia.Add(new Media(fModList(now, colon-now)));
      now = colon+1;
    }
    if (now < fModList.Length()) 
      fMedia.Add(new Media(now, fModList.Length()-now));
    if (fMedia.GetEntries() <= 0) return kFALSE;
    return AliFMDInputHits::Init();
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
    // Get production vertex 
    Double_t vx = track->Vx();
    Double_t vy = track->Vy();
    Double_t vz = track->Vz();
    
    fAll->fCounts->Fill(fEv, hit->Detector(),Int_t(hit->Ring()),
			hit->Sector(), hit->Strip(), 
			hit->X(), hit->Y(), hit->Z(), hit->Edep(), 
			hit->Pdg());
    // Get node 
    TGeoNode* prodNode = gGeoManager->FindNode(vx,vy,vz);
    if (!prodNode) return kTRUE;
    
    // Get volume 
    TGeoVolume* prodVol = prodNode->GetVolume();
    if (!prodVol) return kTRUE;
    
    // Med medium 
    TGeoMedium* prodMed = prodVol->GetMedium();
    if (!prodMed) return kTRUE;
    
    Media* media = FindMedia(prodMed->GetUniqueID());
    if (media) media = fOther; 
    
    // r = TMath::Sqrt(hit->X() * hit->X() + hit->Y() * hit->Y());
    // if(r!=0) Float_t wgt=1./(2. * 3.1415*r);
    media->fCounts->Fill(fEv, hit->Detector(),Int_t(hit->Ring()),
			 hit->Sector(), hit->Strip(), 
			 hit->X(), hit->Y(), hit->Z(), hit->Edep(), 
			 hit->Pdg());
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t Finish()
  {
    fOutput->Write();
    fOutput->Close();
    return kTRUE;
  }
};

//____________________________________________________________________
//
// EOF
//
