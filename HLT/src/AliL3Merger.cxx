
//Author:        Uli Frankenfeld
//Last Modified: 06.12.2000

#include "AliL3Logging.h"
#include <math.h>
#include <iostream.h>
#include "AliL3Merger.h"
#include "AliL3Track.h"
#include "AliL3TrackSegmentData.h"
#include "AliL3Transform.h"
#include "AliL3TrackArray.h"

#ifdef use_root //use root ntuple for slow merge
#include <TNtuple.h>
#include <TTree.h>
#include <TFile.h>
#endif
//_____________________________________________________________
//
// The L3 merger base class
//

ClassImp(AliL3Merger)

AliL3Merger::AliL3Merger(){
  //Default constructor
  fTransformer= 0;
  SetArray(0);
}


AliL3Merger::AliL3Merger(Int_t ntrackarrays){
  //Constructor.
  SetArray(ntrackarrays);
  fCurrentTracks=0;
}

AliL3Merger::~AliL3Merger(){
  //Destructor
  DeleteArray();
}

void AliL3Merger::DeleteArray(){
  for(Int_t i=0; i<fNIn;i++)
    delete fInTrack[i];
  delete[] (Byte_t*) fInTrack;
  delete fOutTrack;
}

void AliL3Merger::SetArray(Int_t nin){
  fNIn = nin;
  fInTrack = (AliL3TrackArray **) new Byte_t[fNIn*sizeof(AliL3TrackArray *)];
  for(Int_t i=0; i<fNIn;i++){
    fInTrack[i] = new AliL3TrackArray();
  }
  fOutTrack= new AliL3TrackArray();
}

void AliL3Merger::Reset(){
  for(Int_t i=0; i<fNIn;i++){
    fInTrack[i]->Reset();
  }
  fOutTrack->Reset();
}

void AliL3Merger::FillTracks(Int_t ntracks, AliL3TrackSegmentData* tr){
  //Read tracks from shared memory (or memory)
  AliL3TrackArray *destination = GetInTracks(fCurrentTracks);
  if(Is2Global())
    destination->FillTracks(ntracks, tr, fSlice, fTransformer);
  else
    destination->FillTracks(ntracks, tr);
}

void AliL3Merger::AddAllTracks(){
  for(Int_t i=0; i<GetNIn();i++){
    AliL3TrackArray *in = GetInTracks(i);
    AliL3TrackArray *out = GetOutTracks();
    out->AddTracks(in);
  }
}

void AliL3Merger::SortGlobalTracks(AliL3Track **tracks, Int_t ntrack){
  AliL3Track **tmp = new AliL3Track*[ntrack]; 
  for(Int_t i=0;i<ntrack;i++) tmp[i] = tracks[i];
  Int_t *t = new Int_t[ntrack];
  for(Int_t i=0;i<ntrack;i++) t[i]=-1;

  for(Int_t j=0;j<ntrack;j++){
    Double_t minr=300;
    Int_t    mini=0;
    for(Int_t i=0;i<ntrack;i++){
     if(!tracks[i]) continue;
       Double_t rr=pow(tracks[i]->GetFirstPointX(),2)+pow(tracks[i]->GetFirstPointY(),2);
       Double_t r=sqrt(rr);
       if(r<minr){
         minr=r;
         mini=i;
       }
    }
    t[j]=mini;
    tracks[mini]=0;
  }
  for(Int_t i=0;i<ntrack;i++) tracks[i] = tmp[t[i]];
  delete[] t;
  delete[] tmp;
}


void AliL3Merger::SortTracks(AliL3Track **tracks, Int_t ntrack){
  AliL3Track **tmp = new  AliL3Track*[ntrack];
  for(Int_t i=0;i<ntrack;i++) tmp[i] = tracks[i];
  Int_t *t = new Int_t[ntrack];
  for(Int_t i=0;i<ntrack;i++) t[i]=-1;

  for(Int_t j=0;j<ntrack;j++){
    Double_t minx=300; 
    Int_t    mini=0;
    for(Int_t i=0;i<ntrack;i++){
     if(!tracks[i]) continue;
       if(tracks[i]->GetFirstPointX()<minx){
         minx=tracks[i]->GetFirstPointX();
         mini=i;
       }     
    }
    t[j]=mini;  
    tracks[mini]=0;
  }
  for(Int_t i=0;i<ntrack;i++) tracks[i] = tmp[t[i]];
  delete[] t;
  delete[] tmp;
}

void AliL3Merger::AddTrack(AliL3TrackArray *mergedtrack,AliL3Track *track){
  AliL3Track *t[1];
  t[0] = track;
  MultiMerge(mergedtrack,t,1);
}

AliL3Track * AliL3Merger::MergeTracks(AliL3TrackArray *mergedtrack,AliL3Track *t0,AliL3Track *t1){
  AliL3Track *t[2];
  t[0] = t0; 
  t[1] = t1;
  SortTracks(t,2);
  return MultiMerge(mergedtrack,t,2);
}

AliL3Track * AliL3Merger::MultiMerge(AliL3TrackArray *mergedtracks,AliL3Track **tracks, Int_t ntrack){
// merge the tracks!!

//check npoints
  Int_t nps = 0;
  for(Int_t i=0;i<ntrack;i++){
    nps+=tracks[i]->GetNHits();
  }
  if(nps>174){
    LOG(AliL3Log::kWarning,"AliL3Merger::MultiMerge","Adding Points")
    <<AliL3Log::kDec<<"Too many Points: "<<nps<<ENDLOG;
    return 0;
  }

  //create new track
  AliL3Track *newtrack = mergedtracks->NextTrack();
  //copy points
  UInt_t nn[174];
  nps = 0;

//  for(Int_t i=0;i<ntrack;i++){
  for(Int_t i=ntrack-1;i>=0;i--){
    memcpy(&nn[nps],tracks[i]->GetHitNumbers(),tracks[i]->GetNHits()*sizeof(UInt_t));
    nps+=tracks[i]->GetNHits();
  }
  AliL3Track *tpf=tracks[0];
  AliL3Track *tpl=tracks[ntrack-1];

  newtrack->SetNHits(nps);
  newtrack->SetHits(nps,nn);
  newtrack->SetFirstPoint(tpf->GetFirstPointX(),tpf->GetFirstPointY(),tpf->GetFirstPointZ());
  newtrack->SetLastPoint(tpl->GetLastPointX(),tpl->GetLastPointY(),tpl->GetLastPointZ());
  newtrack->SetPt(tpf->GetPt());
  newtrack->SetPsi(tpf->GetPsi());
  newtrack->SetTgl(tpf->GetTgl());
  newtrack->SetCharge(tpf->GetCharge());
  return newtrack;
}

void* AliL3Merger::GetNtuple(char *varlist){
#ifdef use_root
  TNtuple* nt = new TNtuple("ntuple","ntuple",varlist);
  return (void*) nt;
#else
  return 0;
#endif
}

void* AliL3Merger::GetNtuple(){
#ifdef use_root
  TNtuple* nt = new TNtuple("ntuple","ntuple",
       "dx:dy:dz:dk:dpsi:dtgl:dq:disx:disy:disz:dis:n0:n1:diff:drx:dry:drz");
  return (void*) nt;
#else
  return 0;
#endif
}

Bool_t AliL3Merger::WriteNtuple(char *filename, void* nt){
#ifdef use_root
  TNtuple *ntuple=(TNtuple *) nt;
  TFile *f = new TFile(filename,"RECREATE");
  ntuple->Write();
  f->Close();
  delete ntuple; 
  return kTRUE; 
#else
  return kFALSE;
#endif
}

void AliL3Merger::FillNtuple(void *nt,AliL3Track *innertrack,AliL3Track *outertrack){
  Float_t data[17];
  if(outertrack->IsPoint()&&innertrack->IsPoint()){
    data[0] =Float_t(innertrack->GetPointX()-outertrack->GetPointX());
    data[1] =Float_t(innertrack->GetPointY()-outertrack->GetPointY());
    data[2] =Float_t(innertrack->GetPointZ()-outertrack->GetPointZ());
    data[3] =Float_t(innertrack->GetKappa()-outertrack->GetKappa());
    Double_t psi= innertrack->GetPointPsi() - outertrack->GetPointPsi();
    if(psi>PI) psi-=2*PI;
    if(psi<-PI)psi+=2*PI;
    data[4] =Float_t(psi);
    data[5] =Float_t(innertrack->GetTgl()-outertrack->GetTgl());
    data[6] =Float_t(innertrack->GetCharge()-outertrack->GetCharge());
    data[7] =Float_t(innertrack->GetLastPointX()-outertrack->GetFirstPointX());
    data[8] =Float_t(innertrack->GetLastPointY()-outertrack->GetFirstPointY());
    data[9] =Float_t(innertrack->GetLastPointZ()-outertrack->GetFirstPointZ());
    data[10] =sqrt(pow(data[7],2)+pow(data[8],2)+pow(data[9],2));
    data[11]= outertrack->GetNHits();
    data[12]= innertrack->GetNHits();
    data[13] = Float_t(TrackDiff(innertrack,outertrack));
#ifdef use_root
    TNtuple *ntuple = (TNtuple *) nt;
    ntuple->Fill(data);
#endif
  }
}

void AliL3Merger::FillNtuple(void *nt,Float_t *data){
#ifdef use_root
    TNtuple *ntuple = (TNtuple *) nt;
    ntuple->Fill(data);
#endif
}

Double_t AliL3Merger::GetAngle(Double_t a1,Double_t a2){
  Double_t da = a1 - a2 +4*PI;
  da = fmod(da,2*PI);
  if(da>PI) da = 2*PI -da;
  return da;
}

void AliL3Merger::SetParameter(Double_t maxy, Double_t maxz, Double_t maxkappa, Double_t maxpsi, Double_t maxtgl){
  fMaxY = maxy;
  fMaxZ = maxz;
  fMaxKappa = maxkappa;
  fMaxPsi = maxpsi;
  fMaxTgl = maxtgl;
}

Bool_t AliL3Merger::IsTrack(AliL3Track *innertrack,AliL3Track *outertrack){

  if(innertrack->GetCharge()!=outertrack->GetCharge()) return kFALSE;
  if( (!innertrack->IsPoint()) || (!outertrack->IsPoint()) )  return kFALSE; 
  if(innertrack->GetNHits()+outertrack->GetNHits()>174) return kFALSE;

  if(fabs(innertrack->GetPointY()-outertrack->GetPointY()) >fMaxY) return kFALSE;
  if(fabs(innertrack->GetPointZ()-outertrack->GetPointZ()) >fMaxZ) return kFALSE;
  if(fabs(innertrack->GetKappa()-outertrack->GetKappa())   >fMaxKappa) return kFALSE;
  if(GetAngle(innertrack->GetPointPsi(),outertrack->GetPointPsi()) >fMaxPsi) return kFALSE;
  if(fabs(innertrack->GetTgl()-outertrack->GetTgl()) >fMaxTgl) return kFALSE;
  //if no rejection up to this point: merge!!
  return kTRUE;
}

Bool_t AliL3Merger::IsRTrack(AliL3Track *innertrack,AliL3Track *outertrack){
  return IsTrack(innertrack,outertrack);
}

Double_t AliL3Merger::TrackDiff(AliL3Track *innertrack,AliL3Track *outertrack){
  Double_t diff =-1;
  Double_t x[4],y[4],z[4],dy[4],dz[4];
  AliL3Track *tracks[2]; 

  tracks[0] = innertrack;
  tracks[1] = outertrack;
  SortGlobalTracks(tracks,2);
  innertrack = tracks[0]; 
  outertrack = tracks[1];

  x[0] = innertrack->GetFirstPointX();
  x[1] = innertrack->GetLastPointX();
  x[2] = outertrack->GetFirstPointX();
  x[3] = outertrack->GetLastPointX();

  y[0] = innertrack->GetFirstPointY();
  y[1] = innertrack->GetLastPointY();
  y[2] = outertrack->GetFirstPointY();
  y[3] = outertrack->GetLastPointY();

  z[0] = innertrack->GetFirstPointZ();
  z[1] = innertrack->GetLastPointZ();
  z[2] = outertrack->GetFirstPointZ();
  z[3] = outertrack->GetLastPointZ();


  outertrack->CalculatePoint(x[0]);
  if(!outertrack->IsPoint()) return diff;
  dy[0] = fabs(y[0] - outertrack->GetPointY());
  dz[0] = fabs(z[0] - outertrack->GetPointZ());

  outertrack->CalculatePoint(x[1]);
  if(!outertrack->IsPoint()) return diff;
  dy[1] = fabs(y[1] - outertrack->GetPointY());
  dz[1] = fabs(z[1] - outertrack->GetPointZ());

  innertrack->CalculatePoint(x[2]);
  if(!innertrack->IsPoint()) return diff;
  dy[2] = fabs(y[2] - innertrack->GetPointY());
  dz[2] = fabs(z[2] - innertrack->GetPointZ());

  innertrack->CalculatePoint(x[3]);
  if(!innertrack->IsPoint()) return diff;
  dy[3] = fabs(y[3] - innertrack->GetPointY());
  dz[3] = fabs(z[3] - innertrack->GetPointZ());

  diff=0;
  for(Int_t i=0;i<4;i++)
    diff+=sqrt(dy[i]*dy[i]+dz[i]*dz[i]);
  return diff; 
}



void AliL3Merger::PrintDiff(AliL3Track *innertrack,AliL3Track *outertrack){
  if(!innertrack->IsPoint()||!outertrack->IsPoint()){
    cerr<<"AliL3Merger::PrintDiff: No Points"<<endl;
    cerr<<"---------------------------"<<endl;
    return;
  } 
  Double_t dx = innertrack->GetPointX()-outertrack->GetPointX();
  Double_t dy = innertrack->GetPointY()-outertrack->GetPointY();
  Double_t dz = innertrack->GetPointZ()-outertrack->GetPointZ();
  Double_t dk = innertrack->GetKappa()-outertrack->GetKappa();
    Double_t dpsi= innertrack->GetPointPsi() - outertrack->GetPointPsi();
    if(dpsi>PI) dpsi-=2*PI;
    if(dpsi<-PI)dpsi+=2*PI;
//  Double_t dpsi = GetAngle(innertrack->GetPointPsi(),outertrack->GetPointPsi());
  Double_t dtgl= innertrack->GetTgl()-outertrack->GetTgl();
  Double_t dq =innertrack->GetCharge()-outertrack->GetCharge();

  fprintf(stderr,"dx: %4f dy: %4f dz: %4f dk: %4f dpsi: %4f dtgl: %4f dq: %4f\n",
      dx,dy,dz,dk,dpsi,dtgl,dq);
cerr<<"---------------------------"<<endl;
//  cerr<<endl; 

}

void AliL3Merger::Print(){
  // print some infos
  for(Int_t i=0; i<fNIn; i++){
    AliL3TrackArray *ttt= GetInTracks(i);
    for(Int_t j =0;j<ttt->GetNTracks();j++){
      AliL3Track *track=ttt->GetCheckedTrack(j);
      if(!track) continue;
      track->CalculateHelix();
//      Double_t angle = atan2(track->GetLastPointY(),track->GetLastPointX());
//      if(angle<0) angle+=PI;
      if(track->CalculatePoint(135))
//      if(!track->CalculateEdgePoint(angle)) cerr<<"**************"<<endl;     
//      if(track->CalculatePoint(track->GetLastPointX()))
//      if(track->CalculatePoint(0))
      {
//      PrintTrack(track);
//      track->CalculateReferencePoint(PI/180.);
      track->CalculateReferencePoint(0.001);
      fprintf(stderr,"npt: %3d dx: %8.5f dy: %8.5f dz: %8.5f\n",
        track->GetNHits(),(float)track->GetPointX()-track->GetPointX(),
                          (float)track->GetPointY()-track->GetPointY(),
                          (float)track->GetPointZ()-track->GetPointZ());

      cerr<<"---------------------------"<<endl;
      }
    }  
  }
}

void AliL3Merger::PrintTrack(AliL3Track *track){
  fprintf(stderr,"npt: %3d pt: %.2f psi: %.2f tgl: %5.2f q: %2d\n",
    track->GetNHits(),track->GetPt(),track->GetPsi(),
    track->GetTgl(),track->GetCharge());
  fprintf(stderr,
    "x1: %6.2f y1: %6.2f z1: %6.2f xl: %6.2f yl: %6.2f zl: %6.2f\n",
    track->GetFirstPointX(),track->GetFirstPointY(),track->GetFirstPointZ(),
    track->GetLastPointX(),track->GetLastPointY(),track->GetLastPointZ());
  if(track->IsPoint()){
    fprintf(stderr,
      "R: %.2f Xc: %.2f Yc: %.2f Xp: %.2f Yp: %.2f Zp: %.2f Psip: %.2f\n",
      track->GetRadius(),track->GetCenterX(),track->GetCenterY(),
      track->GetPointX(),track->GetPointY(),track->GetPointZ(),
      track->GetPointPsi());
  }
}




