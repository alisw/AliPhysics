//////////////////////////////////////////////////////////////////////////////
//
// A much longer description of event info
//
//////////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>
#include <TMath.h>

#include <AliLog.h>
#include <AliVParticle.h>
#include "AliMEStrackInfo.h"
#include "AliMESeventInfo.h"

ClassImp(AliMESeventInfo)
ClassImp(AliMESeventInfo::AliMESevShape)

//______________________________________________________________
AliMESeventInfo::AliMESeventInfo()
  : TObject()
  ,fQuality(0)
  ,fVertexZ(0.)
  ,fEvShape()
{
  //
  // Constructor
  //
// 	memset(fMultiplicity, 0, kNmult*sizeof(Int_t));
	memset(fMultiplicity, 0, kNmult*sizeof(Double_t));
}

//______________________________________________________________
void AliMESeventInfo::Clear(Option_t *)
{
  // Reset info
  fQuality = 0;
  fVertexZ = 0.;
//   memset(fMultiplicity, 0, kNmult*sizeof(Int_t));
  memset(fMultiplicity, 0, kNmult*sizeof(Double_t));
  fEvShape.fSphericity = 0.;
  fEvShape.fThrust[0] = 0.; fEvShape.fThrust[1] = 0.;
  fEvShape.fRecoil=0.;
  fEvShape.fDir[0] = 0.; fEvShape.fDir[1] = 0.;
  memset(fEvShape.fFW, 0, FW_MAX_ORDER*sizeof(Double_t));
  fEvShape.fPxyLead[0] = 0.; fEvShape.fPxyLead[1] = 0.;
}

//______________________________________________________________
// fill ev shape object
// directivity
void AliMESeventInfo::MakeDirectivity(TObjArray* tracks){

  Double_t rv[2] = {-2., -2.};

  if(!tracks->GetEntries()){
    AliDebug(2, "Failed event shape estimation. No tracks in event.");
    return;
  }
  if((rv[0] = Directivity(tracks, kTRUE)) < 0. ){
    AliDebug(2, "Failed D+ estimation");
    // return kFALSE;
  }
  if((rv[1] = Directivity(tracks, kFALSE)) < 0. ){
    AliDebug(2, "Failed D- estimation");
    // return kFALSE;
  }
  memcpy(fEvShape.fDir, rv, 2*sizeof(Double_t));

  return;
}

// thrust
Bool_t AliMESeventInfo::MakeThrust(TObjArray* tracks){

    if(!tracks->GetEntries()){
      AliInfo("Failed event shape estimation. No tracks in event.");
      return kFALSE;
    }
    Double_t rv[2] = {0.};
    if(!Thrust(tracks, rv)){
        AliInfo("Failed Thrust estimation");
        //return kFALSE;
    }
    memcpy(fEvShape.fThrust, rv, 2*sizeof(Double_t));

    return kTRUE;
}

// sphericity
void AliMESeventInfo::MakeSphericity(TObjArray* tracks){

    Double_t rv = -1.;

    if(!tracks->GetEntries()){
      AliDebug(2, "Failed event shape estimation. No tracks in event.");
      return;
    }
    if((rv = Sphericity(tracks)) < 0. ){
        AliDebug(2, "Failed Sphericity estimation");
        //return kFALSE;
    }
    fEvShape.fSphericity = rv;

    return;
}

// recoil
Bool_t AliMESeventInfo::MakeRecoil(TObjArray* tracks){

    if(!tracks->GetEntries()){
      AliInfo("Failed event shape estimation. No tracks in event.");
      return kFALSE;
    }
    Double_t rv = 0.;
    if((rv = Recoil(tracks)) < 0. ){
    AliInfo("Failed Recoil estimation");
    //return kFALSE;
    }
    fEvShape.fRecoil = rv;

    return kTRUE;
}

// fox wolfram moments
Bool_t AliMESeventInfo::MakeFoxWolframMoments(TObjArray* tracks){

    if(!tracks->GetEntries()){
      AliInfo("Failed event shape estimation. No tracks in event.");
      return kFALSE;
    }
    Double_t fw[FW_MAX_ORDER]={0.};
    if(!TransverseFoxWolframMoments(tracks, fw)){
    AliInfo("Failed Fox-Wolfram moments estimation");
    //return kFALSE;
    }
    memcpy(fEvShape.fFW, fw, FW_MAX_ORDER*sizeof(Double_t));

    return kTRUE;
}

// get leading particle
Bool_t AliMESeventInfo::FindLeadingParticle(TObjArray* tracks){

    if(!tracks->GetEntries()){
      AliDebug(2, "Failed event shape estimation. No tracks in event.");
      return kFALSE;
    }
    Double_t rv[2] = {0.};
    if(!LeadingParticleDirection(tracks, rv)){
      AliDebug(2, "Failed LeadingParticleDirection estimation");
      return kFALSE;
    }
    memcpy(fEvShape.fPxyLead, rv, 2*sizeof(Double_t));

    return kTRUE;
}

//______________________________________________________________
Double_t AliMESeventInfo::Directivity(TObjArray* tracks, Bool_t etaSign)
{
  // compute directivity
  Int_t ntracks(0);
  if(!(ntracks=tracks->GetEntries())) return -2.;

  Double_t dirx(0.), diry(0.), dir(0.);
  AliMEStrackInfo *track(NULL);
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    // printf("part %i:\t", iTracks);
    if(!(track = dynamic_cast<AliMEStrackInfo*> (tracks->At(iTracks)))) continue;
    // printf("primary: %i, eta: %f\n", track->HasOrigin(AliMEStrackInfo::kPrimary), track->Eta());
    if(! track->HasOrigin(AliMEStrackInfo::kPrimary) ) continue;
    if(track->Eta()*(etaSign?1.:-1.) < 0.) continue;
    if(TMath::Abs(track->Eta()) > 0.8) continue;
    dirx+=track->Px();
    diry+=track->Py();
    dir +=track->Pt();
    // printf("track selected\n\n");
  }
  return dir>kAlmost0?(TMath::Sqrt(dirx*dirx+diry*diry)/dir):-2.;
}

//______________________________________________________________
Bool_t AliMESeventInfo::LeadingParticleDirection(TObjArray* tracks, Double_t pxy[2])
{
  // compute leading particle direction

  memset(pxy, 0, 2*sizeof(Double_t));
  Int_t ntracks(0);
  if(!(ntracks=tracks->GetEntries())) return kFALSE;

  Double_t ptmax(0.);
  Double_t etamax(0.);
  Double_t phimax(0.);
  // AliVParticle *track(NULL);
  AliMEStrackInfo *track(NULL);
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    // if(!(track = dynamic_cast<AliVParticle*> (tracks->At(iTracks)))) continue;
    if(!(track = dynamic_cast<AliMEStrackInfo*> (tracks->At(iTracks)))) continue;
    if(TMath::Abs(track->Eta()) >= 0.8) continue;
    if(! track->HasOrigin(AliMEStrackInfo::kPrimary) ) continue;
    if(track->Pt()<ptmax) continue;
    pxy[0] = track->Px();
    pxy[1] = track->Py();
    ptmax  = track->Pt();
    etamax = track->Eta();
    phimax = track->Phi();
  }

  return kTRUE;
}


//______________________________________________________________
Bool_t AliMESeventInfo::Thrust(TObjArray* tracks, Double_t t[2])
{
  // compute thrust value[0] and direction[1]

  memset(t, 0, 2*sizeof(Double_t));
  Int_t ntracks(0);
  if(!(ntracks=tracks->GetEntries())) return kFALSE;

  Double_t deltaPhi=0.05*TMath::Pi()/180., //the resolution of the thrust orientation
           absPtSum(0.);

  AliVParticle *track(NULL);
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    if(!(track = dynamic_cast<AliVParticle*>(tracks->At(iTracks)))) continue;
    absPtSum += track->Pt();
  }

  for (Double_t nphi(0.); nphi<=TMath::TwoPi(); nphi+=deltaPhi){//nphi loop
    Double_t thru(0.);
    for(Int_t j=0;j<ntracks;j++){
      if(!(track = dynamic_cast<AliVParticle*> (tracks->At(j)))) continue;
      thru+=track->Pt()*TMath::Abs(TMath::Cos(track->Phi()-nphi));
    }
    if(thru>t[0]){
      t[0] = thru;
      t[1] = nphi;
    }
  }//nphi loop

  if(absPtSum<kAlmost0) return kFALSE;
  t[0] /= absPtSum;
  return kTRUE;
}

/*
//______________________________________________________________
Bool_t AliMESeventInfo::Thrust(TObjArray* tracks, Double_t t[2])
{
    return kFALSE;
}
*/
//______________________________________________________________
Double_t AliMESeventInfo::Sphericity(TObjArray* tracks)
{
    // compute sphericity

  Int_t ntracks(0);
  if(!(ntracks=tracks->GetEntries())) return -1.;

  Double_t a(0.), b(0.), c(0.), d(0.);
  AliMEStrackInfo *track(NULL);
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    if(!(track = dynamic_cast<AliMEStrackInfo*> (tracks->At(iTracks)))) continue;
    if(! track->HasOrigin(AliMEStrackInfo::kPrimary) ) continue;
    if(TMath::Abs(track->Eta()) > 0.8) continue;
    a+=track->Px()*track->Px();
    b+=track->Px()*track->Py();
    d+=track->Py()*track->Py();
  }
  c = b;
  Double_t ad=a+d,
           delta=TMath::Sqrt((ad)*(ad)-4.*(a*d-b*c)),
           lambda1=(ad-delta)/2.;
           //lambda2=(ad+delta)/2.;
  return ad>kAlmost0?(2.*lambda1/ad):-1;
}
/*
//______________________________________________________________
Double_t AliMESeventInfo::Sphericity(TObjArray* tracks)
{
    return 0;
}
*/

//______________________________________________________________
Double_t AliMESeventInfo::Recoil(TObjArray* tracks)
{
  // compute recoil

  Int_t ntracks(0);
  if(!(ntracks=tracks->GetEntries())) return -1.;

  Double_t dirx(0.), diry(0.), dir(0.);
  AliVParticle *track(NULL);
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    if(!(track = dynamic_cast<AliVParticle*> (tracks->At(iTracks)))) continue;
    dirx+=track->Px();
    diry+=track->Py();
    dir +=track->Pt();
  }
  return dir>kAlmost0?(TMath::Sqrt(dirx*dirx+diry*diry)/dir):-1.;
}
/*
//______________________________________________________________
Double_t AliMESeventInfo::Recoil(TObjArray* tracks)
{
    return 0;
}
*/

//______________________________________________________________
Bool_t AliMESeventInfo::TransverseFoxWolframMoments(TObjArray* tracks, Double_t fw[FW_MAX_ORDER])
{
  memset(fw, 0, FW_MAX_ORDER*sizeof(Double_t));
  Int_t ntracks(0);
  if(!(ntracks=tracks->GetEntries())) return kFALSE;
  Double_t legendre[FW_MAX_ORDER+1] = {0.}, sumcl[FW_MAX_ORDER+1] = {0.};
  legendre[0]=1.;  // always unity

  AliVParticle *ti(NULL), *tj(NULL);
  for (Int_t i(0); i < ntracks; i++) {
    if(!(ti = dynamic_cast<AliVParticle*> (tracks->At(i)))) continue;
    Double_t pti(ti->Pt()),
             phii(ti->Phi());
    for (Int_t j(0); j < ntracks; j++) {
      if(!(tj = dynamic_cast<AliVParticle*> (tracks->At(j)))) continue;
      Double_t ptj(tj->Pt()),
               phij(tj->Phi()),
               cosPhiij(TMath::Cos(phii-phij)),
               ptij(pti*ptj*(i==j?1.:2.));

      legendre[1]=cosPhiij;

      for(Int_t order=0;order<FW_MAX_ORDER+1;order++){
        if(order>1){
          Double_t alfa=cosPhiij*legendre[order-1];
          legendre[order]=alfa+(alfa-legendre[order-2])*(1-1./order);
        }
        sumcl[order]+=ptij*legendre[order];
      }
    }
  }
  if(sumcl[0]<kAlmost0) return kFALSE;
  for(Int_t order=0;order<FW_MAX_ORDER;order++) fw[order]=sumcl[order+1]/sumcl[0];
  return kTRUE;
}
/*
//______________________________________________________________
Bool_t AliMESeventInfo::TransverseFoxWolframMoments(TObjArray* tracks, Double_t fw[FW_MAX_ORDER])
{
return kFALSE;
}
*/
//______________________________________________________________
void AliMESeventInfo::Print(Option_t */*o*/) const
{
  //
  // Dump event info to terminal
  //
  printf("Event Info  ::\n"
         "PileUp       : [%c]\n"
         "Multiplicity : [%g %g %g]\n"
         "Vertex       : [%c] Type[%s]\n"
         "Trigger      : MB[%c] HM[%c]\n",
    IsPileUp()?'y':'n',
    fMultiplicity[kGlob08], fMultiplicity[kComb], fMultiplicity[kNoClSPD],
    HasVertex()?'y':'n', HasVertexGlobal()?"Global":"ITS",
    HasTriggerMB()?'y':'n', HasTriggerHM()?'y':'n');
  printf("Event Shape  : D+[%f] D-[%f] T[%f %f] S[%f] R[%f] Leading(px, py)[%f %f]\n",
    fEvShape.fDir[0], fEvShape.fDir[1], fEvShape.fThrust[0], fEvShape.fThrust[1],  fEvShape.fSphericity,  fEvShape.fRecoil, fEvShape.fPxyLead[0], fEvShape.fPxyLead[1]);
  printf("             : FW[");
  for(Int_t ifw(0); ifw<FW_MAX_ORDER; ifw++) printf("%f ", fEvShape.fFW[ifw]); printf("]\n");
}

//______________________________________________________________
void AliMESeventInfo::SetEvShape(const AliMESevShape& ev)
{
  //
  // Set event shape from source
  //
  fEvShape = AliMESevShape(ev);
}

//______________________________________________________________
void AliMESeventInfo::SetEvShape(Double_t dir[2], Double_t sfr, Double_t tr[2], Double_t rec, Double_t fw[FW_MAX_ORDER], Double_t pxy[2])
{
  //
  // Set event shape from scratch
  //
  fEvShape = AliMESevShape(dir, sfr, tr, rec, fw, pxy);
}

//______________________________________________________________
AliMESeventInfo::AliMESevShape::AliMESevShape()
  : TNamed()
  ,fSphericity(0.)
  ,fRecoil(0.)
{
  //
  // Constructor
  //
  fDir[0] = 0.; fDir[1] = 0.;
  fThrust[0] = 0.; fThrust[1] = 0.;
  memset(fFW, 0, FW_MAX_ORDER*sizeof(Double_t));
  fPxyLead[0] = 0.; fPxyLead[1] = 0.;
}

//______________________________________________________________
AliMESeventInfo::AliMESevShape::AliMESevShape(const char* name, const char* title)
  : TNamed(name, title)
  ,fSphericity(0.)
  ,fRecoil(0.)
{
  //
  // Constructor
  //
  fDir[0] = 0.; fDir[1] = 0.;
  fThrust[0] = 0.; fThrust[1] = 0.;
  memset(fFW, 0, FW_MAX_ORDER*sizeof(Double_t));
  fPxyLead[0] = 0.; fPxyLead[1] = 0.;
}

//______________________________________________________________
AliMESeventInfo::AliMESevShape::AliMESevShape(Double_t dir[2], Double_t sphr, Double_t tr[2], Double_t rec, Double_t fw[FW_MAX_ORDER], Double_t pxy[2])
  : TNamed()
  ,fSphericity(sphr)
  ,fRecoil(rec)
{
  //
  // Constructor
  //
  memcpy(fDir, dir, 2*sizeof(Double_t));
  memcpy(fThrust, tr, 2*sizeof(Double_t));
  memcpy(fFW, fw, FW_MAX_ORDER*sizeof(Double_t));
  memcpy(fPxyLead, pxy, 2*sizeof(Double_t));
}

//______________________________________________________________
AliMESeventInfo::AliMESevShape::AliMESevShape(const AliMESevShape &evs)
  : TNamed((const TNamed&)evs)
  ,fSphericity(evs.fSphericity)
  ,fRecoil(evs.fRecoil)
{
  //
  // Constructor
  //
  memcpy(fDir, evs.fDir, 2*sizeof(Double_t));
  memcpy(fThrust, evs.fThrust, 2*sizeof(Double_t));
  memcpy(fFW, evs.fFW, FW_MAX_ORDER*sizeof(Double_t));
  memcpy(fPxyLead, evs.fPxyLead, 2*sizeof(Double_t));
}
