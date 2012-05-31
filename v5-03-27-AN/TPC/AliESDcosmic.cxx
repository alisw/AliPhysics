#include "TMath.h"
#include "TArrayI.h"
#include "TDatabasePDG.h"
#include "TTreeStream.h"
#include "TVectorD.h"


#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliTracker.h"
#include "AliESDcosmic.h"

ClassImp(AliESDcosmic)

AliESDcosmic::AliESDcosmic():
    TObject(),
    fESD(0),
    fTracks(0),
    fTracksAcorde(0),
    fPair(0),
    fDebugStreamer(0)
{
  //
  //
  //
}


AliESDcosmic::~AliESDcosmic(){
  //
  //
  //
  delete fPair;
  delete fTracksAcorde;
  
}

void AliESDcosmic::ProcessEvent(AliESDEvent* event){
  //
  //
  //
  fESD=event;
  TString  string =  fESD->GetFiredTriggerClasses();
  if (!(string.Contains("ASL"))) return;
  //
  const Float_t kCutMinDir=0.95;
  const Float_t kMaxD=1000.;
  Int_t ntracks=event->GetNumberOfTracks(); 
  if (ntracks<=0) return;
  if (!fTracks) fTracks= new TClonesArray("AliESDtrack",0);
  fPair = new TArrayI(ntracks);

  for (Int_t i=0;i<ntracks;++i) {
    AliESDtrack *track0 = event->GetTrack(i); 
    AliESDtrack *trackPair=0;
    (*fPair)[i]=-1;
    // track0 - choosen upper part
    if (!track0) continue;
    if (!track0->GetOuterParam()) continue; 
    Double_t dir0[3];
    track0->GetDirection(dir0); 
    Float_t minDist=kMaxD;
    //
    for (Int_t j=i;j<ntracks;++j) {
      AliESDtrack *track1 = event->GetTrack(j);   
      //track 1 lower part
      if (!track1) continue;
      if (!track1->GetOuterParam()) continue;
      Double_t dir1[3];
      track1->GetDirection(dir1);
      Float_t dir = (dir0[0]*dir1[0] + dir0[1]*dir1[1] + dir0[2]*dir1[2]);
      if (TMath::Abs(dir)<kCutMinDir) continue;               // direction vector product 
      //
      // calculate distance
      Float_t dy = (track0->GetY()+track1->GetY());
      Float_t sy2 = track0->GetSigmaY2()+track1->GetSigmaY2();
      Float_t dphi = (track0->GetAlpha()-track1->GetAlpha()-TMath::Pi());
      Float_t sphi2  = track0->GetSigmaSnp2()+track1->GetSigmaSnp2();
      Float_t dtheta = (track0->GetTgl()-track1->GetTgl());
      Float_t stheta2 = track0->GetSigmaTgl2()+track1->GetSigmaTgl2();
      Float_t normDist = TMath::Sqrt(dy*dy/sy2+dphi*dphi/sphi2+dtheta*dtheta/stheta2);
      //
      if (normDist>minDist) continue;    
      minDist = normDist;
      trackPair=track1;
      (*fPair)[i]=j;
    }
    //
  }
  PropagateToAcorde();
}


void  AliESDcosmic::PropagateToAcorde(){
  //
  //
  //
  const Double_t kRL3=510;   // radius of L3 magnet
  const Double_t kxAcorde=850.;
  //
  if (!fTracksAcorde) fTracksAcorde= new TClonesArray("AliExternalTrackParam",0);
  Int_t ntracks=fESD->GetNumberOfTracks(); 
  Int_t counter=0;
  //
  for (Int_t i=0; i<ntracks;i++){
    Int_t index = i;
    AliESDtrack * upperTrack = fESD->GetTrack(i);
    if (upperTrack->GetOuterParam()==0) continue;
    Double_t gxyz[3];
    upperTrack->GetOuterParam()->GetXYZ(gxyz);
    if ((*fPair)[i]>0){
      AliESDtrack * track2 = fESD->GetTrack((*fPair)[i]);
      if (track2->GetOuterParam()){
	Double_t gxyz2[3];
	track2->GetOuterParam()->GetXYZ(gxyz2);
	if (gxyz2[1]>gxyz[1]) {
	  upperTrack=track2;
	  index=(*fPair)[i];
	}
      }
    }
    //
    AliExternalTrackParam *upper = (AliExternalTrackParam *)(upperTrack->GetOuterParam()->Clone());
    Bool_t isOK = upper->PropagateTo(kRL3,fESD->GetMagneticField());
    upper->GetXYZ(gxyz);
    if (gxyz[1]<0) continue;
    for (Int_t iter=0; iter<20;iter++){
      upper->GetXYZ(gxyz);
      Double_t galpha = TMath::ATan2(gxyz[1],gxyz[0]);
      Double_t alpha=galpha;
      galpha*=180/TMath::Pi();
      if (iter>1){
	if (galpha<45.)   alpha = TMath::Pi()/8;
	if (galpha>135.) alpha =  TMath::Pi()*(1-1/8.);
	if (galpha>45.&&galpha<135.) alpha = TMath::Pi()/2.;
      }
      if (isOK) upper->Rotate(alpha);
      if (isOK) isOK = upper->PropagateTo(kxAcorde,0);
    }
    if (isOK) {
      new ((*fTracksAcorde)[index]) AliExternalTrackParam(*upper);
      counter++;
    }    
  }
}

void AliESDcosmic::DumpToTree(){
  //
  //
  //
  TTreeSRedirector * cstream = fDebugStreamer;
  if (!cstream) return;
  if (!fESD) return;
  if (!fTracksAcorde) return;
  Int_t ntracks0 =fESD->GetNumberOfTracks(); 
  Int_t ntracks  = fTracksAcorde->GetEntries();
  Float_t mag = fESD->GetMagneticField();
  Int_t   run = fESD->GetRunNumber();
  Int_t   event = fESD->GetEventNumberInFile();
  AliESDHeader* header = fESD->GetHeader();
  
  (*cstream)<<"eventInfo"<<
    "run="<<run<<
    "event="<<event<<
    "header="<<header<<
    "mag="<<mag<<
    "ntracks0="<<ntracks0<<
    "ntracks="<<ntracks<<    
    "\n";
    
  for (Int_t i=0;i<ntracks;i++){
    if (!fTracksAcorde->At(i)) continue;
    TVectorD gxyz(3);
    TVectorD gpxyz(3);
    AliExternalTrackParam * param = (AliExternalTrackParam *)fTracksAcorde->At(i);
    param->GetXYZ(gxyz.GetMatrixArray());
    param->GetPxPyPz(gpxyz.GetMatrixArray());
    (*cstream) << "esdCosmic" <<  
      "run="<<run<<
      "event="<<event<<
      "header="<<header<<
      "mag="<<mag<<
      "ntracks0="<<ntracks0<<
      "ntracks="<<ntracks<<      
      "tr.="<<param<<
      "p.="<<&gxyz<<
      "m.="<<&gpxyz<<
      "\n";
  }
    
}

