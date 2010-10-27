// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Ovrebekk <st05886@alf.uib.no>,                *                  
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************


/** @file  AliAnalysisTaskD0Trigger.cxx  
    @author Gaute Ovrebekk
    @date 
    @brief An analysis task for the D0 Trigger.    
*/

class AliAnalysisTask;
class AliAnalysisManager;

#include "TH1F.h"
#include "TCanvas.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "TMath.h"
#include "AliExternalTrackParam.h"
#include "AliAODVertex.h"
#include "TDatabasePDG.h"
#include "AliESDtrack.h"
#include "TVector3.h"
#include "AliVertexerTracks.h"
#include "AliKFVertex.h"
#include "TDatabasePDG.h"
#include "TVector3.h"

#include "AliAnalysisTaskD0Trigger.h"

ClassImp(AliAnalysisTaskD0Trigger)

AliAnalysisTaskD0Trigger::AliAnalysisTaskD0Trigger()
:
AliAnalysisTaskSE()
  , fOutputList(0)
  , fPtMin(0.0)
  , fdca(0.0)
  , finvMass(0.0)     
  , fcosThetaStar(0.0)
  , fd0(0.0) 
  , fd0d0(0.0)
  , fcosPoint(0.0)
  , fD0PDG(TDatabasePDG::Instance()->GetParticle(421)->Mass())
  , fD0massHLT(NULL)
  , fD0ptHLT(NULL)
  , fD0massOFF(NULL)
  , fD0ptOFF(NULL)
  , fPos()
  , fNeg()
  , ftwoTrackArray(NULL)
  , fTotalD0HLT(0)
  , fTotalD0OFF(0)
  , fField(0)
  , fNevents(0)
  , fuseKF(false)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container

  //DefineOutput(1, TList::Class());
}

AliAnalysisTaskD0Trigger::AliAnalysisTaskD0Trigger(const char *name,float cuts[7])
:
AliAnalysisTaskSE(name)
  , fOutputList(0)
  , fPtMin(cuts[0])
  , fdca(cuts[1])
  , finvMass(cuts[2])     
  , fcosThetaStar(cuts[3])
  , fd0(cuts[4]) 
  , fd0d0(cuts[5])
  , fcosPoint(cuts[6])
  , fD0PDG(TDatabasePDG::Instance()->GetParticle(421)->Mass())
  , fD0massHLT(NULL)
  , fD0ptHLT(NULL)
  , fD0massOFF(NULL)
  , fD0ptOFF(NULL)
  , fPos()
  , fNeg()
  , ftwoTrackArray(NULL)
  , fTotalD0HLT(0)
  , fTotalD0OFF(0)
  , fField(0)
  , fNevents(0)
  , fuseKF(false)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container

  DefineOutput(1, TList::Class());
}
void AliAnalysisTaskD0Trigger::UserCreateOutputObjects(){
  // Create histograms

  OpenFile(1);

  fOutputList = new TList();
  fOutputList->SetName(GetName());

  fD0massHLT = new TH1F("hMassHLT","D^{0} mass plot from HLT reconstruction",100,1.7,2);
  fD0ptHLT = new TH1F("hPtHLT","D^{0} Pt plot HLT reconstruction",20,0,20);

  fD0massOFF = new TH1F("hMassOFF","D^{0} mass plot Offline reconstruction",100,1.7,2);
  fD0ptOFF = new TH1F("hPtOFF","D^{0} Pt plot Offline reconstruction",20,0,20);

  fOutputList->Add(fD0massHLT);
  fOutputList->Add(fD0ptHLT);
  fOutputList->Add(fD0massOFF);
  fOutputList->Add(fD0ptOFF);
  
}

void AliAnalysisTaskD0Trigger::NotifyRun(){
  // see header file of AliAnalysisTask for documentation
}

void AliAnalysisTaskD0Trigger::UserExec(Option_t *){
  // see header file of AliAnalysisTask for documentation

  AliESDEvent *esdOFF = dynamic_cast<AliESDEvent*>(InputEvent());
  
  if(!esdOFF){
    Printf("ERROR: fESD not available");
    return;
  }
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(fInputHandler);
  AliESDEvent *esdHLT = NULL;   
  if(esdH) esdHLT = esdH->GetHLTEvent();
    
  if(!esdHLT){
    Printf("ERROR: HLTesd not available");
    return;
  }
  
  fNevents++;
  
  ftwoTrackArray = new TObjArray(2);
  Int_t nD0=0;

  fField = esdHLT->GetMagneticField();
  const AliESDVertex* pvHLT = esdHLT->GetPrimaryVertexTracks();
  AliESDVertex *pVertexHLT =  new AliESDVertex(*pvHLT);
  if(pVertexHLT->GetNContributors()<2){
    //Printf("ERROR: Contributiors to HLT Primary vertex is to low or not set");
    return;
  }

  for(Int_t it=0;it<esdHLT->GetNumberOfTracks();it++){
    SingleTrackSelect(esdHLT->GetTrack(it),pVertexHLT);
  }
  
  RecD0(nD0,pVertexHLT,true);
  fTotalD0HLT+=nD0;
  
  fPos.clear();
  fNeg.clear();

  nD0=0;

  fField = esdOFF->GetMagneticField();
  const AliESDVertex* pvOFF = esdOFF->GetPrimaryVertexTracks();
  AliESDVertex *pVertexOFF =  new AliESDVertex(*pvOFF);
  if(pVertexOFF->GetNContributors()<2){
    //Printf("ERROR: Contributiors to OFFLINE Primary vertex is to low or not set");
    return;
  }
  for(Int_t it=0;it<esdOFF->GetNumberOfTracks();it++){
    SingleTrackSelect(esdOFF->GetTrack(it),pVertexOFF);
  }
  
  RecD0(nD0,pVertexOFF,false);  
  fTotalD0OFF+=nD0;
  
  fPos.clear();
  fNeg.clear();

  // Post output data.
  PostData(1, fOutputList);
  ftwoTrackArray->Clear();
  delete pVertexHLT;
  delete pVertexOFF;
}

void AliAnalysisTaskD0Trigger::Terminate(Option_t *){
  Printf("Event Number: %d",fNevents);
  Printf("Total Number of D0 foundfor HLT: %d",fTotalD0HLT);
  Printf("Total Number of D0 found for OFFLINE: %d",fTotalD0OFF);
}

void AliAnalysisTaskD0Trigger::SingleTrackSelect(AliExternalTrackParam* t, AliESDVertex *pV){
  // Offline har || på disse kuttene på de to henfallsproduktene 
  Double_t pv[3];
  pV->GetXYZ(pv);
  
  if(t->Pt()<fPtMin){return;}
  if(TMath::Abs(t->GetD(pv[0],pv[1],fField)) > fd0){return;}
  
  if(t->Charge()>0){
    fPos.push_back(t);
  }
  else{
    fNeg.push_back(t);
  }
}

void AliAnalysisTaskD0Trigger::RecD0(Int_t& nD0, AliESDVertex *pV,bool isHLT){
  // Reconstructing D0
  Double_t starD0,starD0bar,xdummy,ydummy; 
  Double_t d0[2];
  Double_t svpos[3];
  Double_t pvpos[3];
  ftwoTrackArray->Clear();

  if(!pV){
    Printf("No Primary Vertex");
    return;
  }
  pV->GetXYZ(pvpos);
    
  for(UInt_t ip=0;ip<fPos.size();ip++){
    AliExternalTrackParam *tP=fPos[ip];
    for(UInt_t in=0;in<fNeg.size();in++){
      AliExternalTrackParam *tN=fNeg[in];
          
      tP->PropagateToDCA(pV,fField,kVeryBig);  //do I need this??????
      tN->PropagateToDCA(pV,fField,kVeryBig);
      
      Double_t dcatPtN = tP->GetDCA(tN,fField,xdummy,ydummy);
      if(dcatPtN>fdca) { continue; }
      
      ftwoTrackArray->AddAt(tP,0);
      ftwoTrackArray->AddAt(tN,1);
      AliAODVertex *vertexp1n1 = ReconstructSecondaryVertex(ftwoTrackArray,fField,pV,fuseKF);
      if(!vertexp1n1) { 
        ftwoTrackArray->Clear();
        continue; 
      }
      
      vertexp1n1->GetXYZ(svpos);
      
      tP->PropagateToDCA(vertexp1n1,fField,kVeryBig); 
      tN->PropagateToDCA(vertexp1n1,fField,kVeryBig);
      
      if((TMath::Abs(InvMass(tN,tP)-fD0PDG)) > finvMass && TMath::Abs((InvMass(tP,tN))-fD0PDG) > finvMass){continue;}
      CosThetaStar(tN,tP,starD0,starD0bar);
      if(TMath::Abs(starD0) > fcosThetaStar && TMath::Abs(starD0bar) > fcosThetaStar){continue;}
      d0[0] = tP->GetD(pvpos[0],pvpos[1],fField);
      d0[1] = tN->GetD(pvpos[0],pvpos[1],fField);
      if((d0[0]*d0[1]) > fd0d0){continue;}
      if(PointingAngle(tN,tP,pvpos,svpos) < fcosPoint){continue;}
      
      if(isHLT){
	fD0massHLT->Fill(InvMass(tN,tP));
	fD0massHLT->Fill(InvMass(tP,tN));
	fD0ptHLT->Fill(Pt(tP,tN));
      }
      else{
	fD0massOFF->Fill(InvMass(tN,tP));
	fD0massOFF->Fill(InvMass(tP,tN));
	fD0ptOFF->Fill(Pt(tP,tN));
      }
      nD0++;
      delete vertexp1n1;
    }
  }
}

Double_t AliAnalysisTaskD0Trigger::InvMass(AliExternalTrackParam* d1, AliExternalTrackParam* d2)
{
  // Calculating the invariant mass
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

void AliAnalysisTaskD0Trigger::CosThetaStar(AliExternalTrackParam* d1, AliExternalTrackParam* d2,Double_t &D0,Double_t &D0bar)
{
  //Calculating the decay angle
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

Double_t AliAnalysisTaskD0Trigger::PointingAngle(AliExternalTrackParam* n, AliExternalTrackParam* p, Double_t *pv, Double_t *sv)
{
  // Calcutating the pointing angle
  TVector3 mom(n->Px()+p->Px(),n->Py()+p->Py(),n->Pz()+p->Pz());
  TVector3 flight(sv[0]-pv[0],sv[1]-pv[1],sv[2]-pv[2]);
  
  double pta = mom.Angle(flight);

  return TMath::Cos(pta); 
}

AliAODVertex* AliAnalysisTaskD0Trigger::ReconstructSecondaryVertex(TObjArray *trkArray, Double_t b, const AliESDVertex *v, bool useKF)
{
  // Finding the vertex of the two tracks in trkArray
  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;
  
  if(!useKF){
    AliVertexerTracks *vertexer = new AliVertexerTracks(b);
    AliESDVertex* vertex =  const_cast<AliESDVertex*>(v);
    vertexer->SetVtxStart(vertex);
    //if(isESD){vertexESD = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);}
    UShort_t *id = new UShort_t[2];
    AliExternalTrackParam *t1 = (AliExternalTrackParam*) trkArray->At(0);
    AliExternalTrackParam *t2 = (AliExternalTrackParam*) trkArray->At(1);
    id[0]=(UShort_t) t1->GetID();
    id[1]=(UShort_t) t2->GetID();
    vertexESD = (AliESDVertex*)vertexer->VertexForSelectedTracks(trkArray,id);
    delete [] id;
    delete vertexer; vertexer=NULL;
    
    if(!vertexESD) return vertexAOD;
    
    if(vertexESD->GetNContributors()!=trkArray->GetEntriesFast()) { 
      //AliDebug(2,"vertexing failed"); 
      delete vertexESD; vertexESD=NULL;
      delete vertex;
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

Double_t AliAnalysisTaskD0Trigger::Pt(AliExternalTrackParam* d1, AliExternalTrackParam* d2)
{
  //Calculating pT of the two tracks
  Double_t p1[3],p2[3];
  d1->GetPxPyPz(p1);
  d2->GetPxPyPz(p2);
  
  Double_t pt2 = (p1[0]+p2[0])*(p1[0]+p2[0]) + (p1[1]+p2[1])*(p1[1]+p2[1]);

  return TMath::Sqrt(pt2);
}
