// $Id$
#include "AliHLTD0toKpi.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "AliESDtrack.h"
#include "TVector3.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "TObjArray.h"
#include "AliVertexerTracks.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliExternalTrackParam.h"
#include "AliKFVertex.h"

ClassImp(AliHLTD0toKpi)

AliHLTD0toKpi::AliHLTD0toKpi() 
{
}

Double_t AliHLTD0toKpi::InvMass(AliExternalTrackParam* d1, AliExternalTrackParam* d2)
{
  Double_t mpi=TDatabasePDG::Instance()->GetParticle(211)->Mass();
  Double_t mK=TDatabasePDG::Instance()->GetParticle(321)->Mass();

  Double_t energy[2]; 
  energy[1] = TMath::Sqrt(mK*mK+d1->GetP()*d1->GetP());
  energy[0] = TMath::Sqrt(mpi*mpi+d2->GetP()*d2->GetP());

  Double_t p1[3],p2[3];
  d1->GetPxPyPz(p1);
  d2->GetPxPyPz(p2);
  
  Double_t momTot2 = (p1[0]+p2[0])*(p1[0]+p2[0])+
                     (p1[1]+p2[1])*(p1[1]+p2[1])+
                     (p1[2]+p2[2])*(p1[2]+p2[2]);

  return TMath::Sqrt((energy[0]+energy[1])*(energy[0]+energy[1])-momTot2);

}
void AliHLTD0toKpi::cosThetaStar(AliExternalTrackParam* d1, AliExternalTrackParam* d2,Double_t &D0,Double_t &D0bar)
{
  Double_t mD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t mpi=TDatabasePDG::Instance()->GetParticle(211)->Mass();
  Double_t mK=TDatabasePDG::Instance()->GetParticle(321)->Mass();

  Double_t pStar = TMath::Sqrt(TMath::Power(mD0*mD0-mK*mK-mpi*mpi,2.)-4.*mK*mK*mpi*mpi)/(2.*mD0);
 
  Double_t px = d1->Px()+d2->Px();
  Double_t py = d1->Py()+d2->Py();
  Double_t pz = d1->Pz()+d2->Pz();
  Double_t p = TMath::Sqrt(px*px+py*py+pz*pz);
  Double_t energy = TMath::Sqrt(p*p+mD0*mD0);

  Double_t beta = p/energy;
  Double_t gamma = energy/mD0;
  
  Double_t qL;
  TVector3 mom(d1->Px(),d1->Py(),d1->Pz());
  TVector3 momD(px,py,pz);
  qL = mom.Dot(momD)/momD.Mag();

  D0 = (qL/gamma-beta*TMath::Sqrt(pStar*pStar+mK*mK))/pStar;
  
  TVector3 mom2(d2->Px(),d2->Py(),d2->Pz());
  TVector3 momD2(px,py,pz);
  qL = mom2.Dot(momD2)/momD2.Mag();

  D0bar = (qL/gamma-beta*TMath::Sqrt(pStar*pStar+mK*mK))/pStar;

}
Double_t AliHLTD0toKpi::pointingAngle(AliExternalTrackParam* n, AliExternalTrackParam* p, Double_t *pv, Double_t *sv)
{

  TVector3 mom(n->Px()+p->Px(),n->Py()+p->Py(),n->Pz()+p->Pz());
  TVector3 flight(sv[0]-pv[0],sv[1]-pv[1],sv[2]-pv[2]);
  
  double pta = mom.Angle(flight);

  return TMath::Cos(pta); 
}

AliAODVertex* AliHLTD0toKpi::ReconstructSecondaryVertex(TObjArray *trkArray, Double_t b, AliESDVertex *v, bool useKF)
{
  
  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;
  
  if(!useKF){
    AliVertexerTracks *vertexer = new AliVertexerTracks(b);
    vertexer->SetVtxStart(v);
    //if(isESD){vertexESD = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);}
    UShort_t *id = new UShort_t[2];
    AliHLTGlobalBarrelTrack *t1 = (AliHLTGlobalBarrelTrack*) trkArray->At(0);
    AliHLTGlobalBarrelTrack *t2 = (AliHLTGlobalBarrelTrack*) trkArray->At(1);
    id[0]=(UShort_t) t1->GetID();
    id[1]=(UShort_t) t2->GetID();
    vertexESD = (AliESDVertex*)vertexer->VertexForSelectedTracks(trkArray,id);
    delete id;
    delete vertexer; vertexer=NULL;
    
    if(!vertexESD) return vertexAOD;
    
    if(vertexESD->GetNContributors()!=trkArray->GetEntriesFast()) { 
      //AliDebug(2,"vertexing failed"); 
      delete vertexESD; vertexESD=NULL;
      return vertexAOD;
    }
  }
  else{
    AliKFParticle::SetField(b);
    
    AliKFVertex vertexKF;
    
    Int_t nTrks = trkArray->GetEntriesFast();
    for(Int_t i=0; i<nTrks; i++) {
      AliESDtrack *esdTrack = (AliESDtrack*)trkArray->At(i);
      AliKFParticle daughterKF(*esdTrack,211);
      vertexKF.AddDaughter(daughterKF);
    }
    vertexESD = new AliESDVertex(vertexKF.Parameters(),
				 vertexKF.CovarianceMatrix(),
				 vertexKF.GetChi2(),
				 vertexKF.GetNContributors());
  }
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  //dispersion = vertexESD->GetDispersion();
  delete vertexESD; vertexESD=NULL;

  Int_t nprongs= trkArray->GetEntriesFast();
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);

  return vertexAOD;

}
Double_t AliHLTD0toKpi::Pt(AliExternalTrackParam* d1, AliExternalTrackParam* d2)
{
  Double_t p1[3],p2[3];
  d1->GetPxPyPz(p1);
  d2->GetPxPyPz(p2);
  
  Double_t pt2 = (p1[0]+p2[0])*(p1[0]+p2[0]) + (p1[1]+p2[1])*(p1[1]+p2[1]);

  return TMath::Sqrt(pt2);
}
