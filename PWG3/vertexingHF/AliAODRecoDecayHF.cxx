/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// Base class for AOD reconstructed heavy-flavour decay
//
// Author: A.Dainese, andrea.dainese@lnl.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <TVector3.h>
#include <TRandom.h>
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODEvent.h"
#include "AliVertexerTracks.h"
#include "AliExternalTrackParam.h"
#include "AliKFVertex.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"

ClassImp(AliAODRecoDecayHF)

//--------------------------------------------------------------------------
AliAODRecoDecayHF::AliAODRecoDecayHF() :
  AliAODRecoDecay(),
  fOwnPrimaryVtx(0x0),
  fEventPrimaryVtx(),
  fListOfCuts(),
  fd0err(0x0), 
  fProngID(0x0),
  fSelectionMap(0)
{
  //
  // Default Constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF::AliAODRecoDecayHF(AliAODVertex *vtx2,Int_t nprongs,Short_t charge,
				     Double_t *px,Double_t *py,Double_t *pz,
				     Double_t *d0,Double_t *d0err) :
  AliAODRecoDecay(vtx2,nprongs,charge,px,py,pz,d0),
  fOwnPrimaryVtx(0x0),
  fEventPrimaryVtx(),
  fListOfCuts(),
  fd0err(0x0),
  fProngID(0x0),
  fSelectionMap(0)
{
  //
  // Constructor with AliAODVertex for decay vertex
  //
  fd0err = new Double_t[GetNProngs()];
  for(Int_t i=0; i<GetNProngs(); i++) fd0err[i] = d0err[i];
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF::AliAODRecoDecayHF(AliAODVertex *vtx2,Int_t nprongs,Short_t charge,
				     Double_t *d0,Double_t *d0err) :
  AliAODRecoDecay(vtx2,nprongs,charge,d0),
  fOwnPrimaryVtx(0x0),
  fEventPrimaryVtx(),
  fListOfCuts(),
  fd0err(0x0),
  fProngID(0x0),
  fSelectionMap(0)
{
  //
  // Constructor with AliAODVertex for decay vertex and without prongs momenta
  //
  fd0err = new Double_t[GetNProngs()];
  for(Int_t i=0; i<GetNProngs(); i++) fd0err[i] = d0err[i];
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF::AliAODRecoDecayHF(Double_t vtx1[3],Double_t vtx2[3],
				     Int_t nprongs,Short_t charge,
				     Double_t *px,Double_t *py,Double_t *pz,
				     Double_t *d0) :
  AliAODRecoDecay(0x0,nprongs,charge,px,py,pz,d0),
  fOwnPrimaryVtx(0x0),
  fEventPrimaryVtx(),
  fListOfCuts(),
  fd0err(0x0),
  fProngID(0x0), 
  fSelectionMap(0)
{
  //
  // Constructor that can used for a "MC" object
  //

  fOwnPrimaryVtx = new AliAODVertex(vtx1);

  AliAODVertex *vtx = new AliAODVertex(vtx2);
  SetOwnSecondaryVtx(vtx);

}
//--------------------------------------------------------------------------
AliAODRecoDecayHF::AliAODRecoDecayHF(const AliAODRecoDecayHF &source) :
  AliAODRecoDecay(source),
  fOwnPrimaryVtx(0x0),
  fEventPrimaryVtx(source.fEventPrimaryVtx),
  fListOfCuts(source.fListOfCuts),
  fd0err(0x0),
  fProngID(0x0),
  fSelectionMap(source.fSelectionMap)
{
  //
  // Copy constructor
  //
  if(source.GetOwnPrimaryVtx()) fOwnPrimaryVtx = new AliAODVertex(*(source.GetOwnPrimaryVtx()));

  if(source.GetNProngs()>0) {
    fd0err = new Double_t[GetNProngs()];
    memcpy(fd0err,source.fd0err,GetNProngs()*sizeof(Double_t));
    if(source.fProngID) {
      fProngID = new UShort_t[GetNProngs()];
      memcpy(fProngID,source.fProngID,GetNProngs()*sizeof(UShort_t));
    }
  }
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF &AliAODRecoDecayHF::operator=(const AliAODRecoDecayHF &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliAODRecoDecay::operator=(source);

  fEventPrimaryVtx = source.fEventPrimaryVtx;
  fListOfCuts = source.fListOfCuts;
  fSelectionMap = source.fSelectionMap;

  if(source.GetOwnPrimaryVtx()) fOwnPrimaryVtx = new AliAODVertex(*(source.GetOwnPrimaryVtx()));

  if(source.GetNProngs()>0) {
    if(source.fd0err) {
      fd0err = new Double_t[GetNProngs()];
      memcpy(fd0err,source.fd0err,GetNProngs()*sizeof(Double_t));
    }
    if(source.fProngID) {
      fProngID = new UShort_t[GetNProngs()];
      memcpy(fProngID,source.fProngID,GetNProngs()*sizeof(UShort_t));
    }
  }
  return *this;
}
//--------------------------------------------------------------------------
AliAODRecoDecayHF::~AliAODRecoDecayHF() {
  //  
  // Default Destructor
  //
  if(fOwnPrimaryVtx) delete fOwnPrimaryVtx;
  if(fd0err) delete [] fd0err;
  if(fProngID) delete [] fProngID;
}
//---------------------------------------------------------------------------
AliKFParticle *AliAODRecoDecayHF::ApplyVertexingKF(Int_t *iprongs,Int_t nprongs,Int_t *pdgs,Bool_t topoCostraint, Double_t bzkG, Double_t *mass) const {
  //
  // Applies the KF vertexer 
  // Int_t iprongs[nprongs] = indices of the prongs to be used from the vertexer
  // Int_t pdgs[nprongs] = pdgs assigned to the prongs, needed to define the AliKFParticle
  // Bool_t topoCostraint = if kTRUE, the topological constraint is applied
  // Double_t bzkG = magnetic field
  // Double_t mass[2] = {mass, sigma} for the mass constraint (if mass[0]>0 the constraint is applied).
  //

  AliKFParticle::SetField(bzkG);
  AliKFParticle *vertexKF=0;
  
  AliKFVertex copyKF;
  Int_t nt=0,ntcheck=0;

  Double_t pos[3]={0.,0.,0.};
  if(!fOwnPrimaryVtx) {
    printf("AliAODRecoDecayHF::ApplyVertexingKF(): cannot apply because primary vertex is not found\n");
    return vertexKF;
  }
  fOwnPrimaryVtx->GetXYZ(pos);
  Int_t contr=fOwnPrimaryVtx->GetNContributors();
  Double_t covmatrix[6]={0.,0.,0.,0.,0.,0.};
  fOwnPrimaryVtx->GetCovarianceMatrix(covmatrix);
  Double_t chi2=fOwnPrimaryVtx->GetChi2();
  AliESDVertex primaryVtx2(pos,covmatrix,chi2,contr,"Vertex");
 

  if(topoCostraint){
   copyKF=AliKFVertex(primaryVtx2);
   nt=primaryVtx2.GetNContributors();
   ntcheck=nt;
  }

  vertexKF = new AliKFParticle();
  for(Int_t i= 0;i<nprongs;i++){
    Int_t ipr=iprongs[i];
    AliAODTrack *aodTrack = (AliAODTrack*)GetDaughter(ipr);
    if(!aodTrack) {
      printf("AliAODRecoDecayHF::ApplyVertexingKF(): no daughters available\n");
      delete vertexKF; vertexKF=NULL;
      return vertexKF;
    }
    AliKFParticle daughterKF(*aodTrack,pdgs[i]);
    vertexKF->AddDaughter(daughterKF);
    
    if(topoCostraint && nt>0){
      //Int_t index=(Int_t)GetProngID(ipr);
      if(!aodTrack->GetUsedForPrimVtxFit()) continue;
      copyKF -= daughterKF;
      ntcheck--;
    }
  }
  
  if(topoCostraint){
    if(ntcheck>0) {
      copyKF += (*vertexKF);
      vertexKF->SetProductionVertex(copyKF);
    }
 }
  
  if(mass[0]>0.){
    vertexKF->SetMassConstraint(mass[0],mass[1]);
  }
  
  return vertexKF;
}
//---------------------------------------------------------------------------
AliAODVertex* AliAODRecoDecayHF::RemoveDaughtersFromPrimaryVtx(AliAODEvent *aod) {
  //
  // This method returns a primary vertex without the daughter tracks of the 
  // candidate and it recalculates the impact parameters and errors.
  // 
  // The output vertex is created with "new". The user has to 
  // set it to the candidate with SetOwnPrimaryVtx(), unset it at the end 
  // of processing with UnsetOwnPrimaryVtx() and delete it.
  // If a NULL pointer is returned, the removal failed (too few tracks left).
  //
  // For the moment, the primary vertex is recalculated from scratch without
  // the daughter tracks.
  //

  AliAODVertex *vtxAOD = aod->GetPrimaryVertex();
  if(!vtxAOD) return 0;
  TString title=vtxAOD->GetTitle();
  if(!title.Contains("VertexerTracks")) return 0;



  AliVertexerTracks *vertexer = new AliVertexerTracks(aod->GetMagneticField());

  Int_t ndg = GetNDaughters();

  vertexer->SetITSMode();
  vertexer->SetMinClusters(3);
  vertexer->SetConstraintOff();

  if(title.Contains("WithConstraint")) {
    Float_t diamondcovxy[3];
    aod->GetDiamondCovXY(diamondcovxy);
    Double_t pos[3]={aod->GetDiamondX(),aod->GetDiamondY(),0.};
    Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
    AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
    vertexer->SetVtxStart(diamond);
    delete diamond; diamond=NULL;
  }

  Int_t skipped[10];
  Int_t nTrksToSkip=0,id;
  AliAODTrack *t = 0;
  for(Int_t i=0; i<ndg; i++) {
    t = (AliAODTrack*)GetDaughter(i);
    id = (Int_t)t->GetID();
    if(id<0) continue;
    skipped[nTrksToSkip++] = id;
  }
  vertexer->SetSkipTracks(nTrksToSkip,skipped);
  AliESDVertex *vtxESDNew = vertexer->FindPrimaryVertex(aod);

  delete vertexer; vertexer=NULL;

  if(!vtxESDNew) return 0;
  if(vtxESDNew->GetNContributors()<=0) { 
    delete vtxESDNew; vtxESDNew=NULL;
    return 0;
  }

  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vtxESDNew->GetXYZ(pos); // position
  vtxESDNew->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vtxESDNew->GetChi2toNDF();
  delete vtxESDNew; vtxESDNew=NULL;

  AliAODVertex *vtxAODNew = new AliAODVertex(pos,cov,chi2perNDF);

  RecalculateImpPars(vtxAODNew,aod);

  return vtxAODNew;
}
//-----------------------------------------------------------------------------------
void AliAODRecoDecayHF::RecalculateImpPars(AliAODVertex *vtxAODNew,AliAODEvent* aod) {
  //
  // now recalculate the daughters impact parameters
  //
  Double_t dz[2],covdz[3];
  for(Int_t i=0; i<GetNDaughters(); i++) {
    AliAODTrack *t = (AliAODTrack*)GetDaughter(i);
    AliExternalTrackParam etp; etp.CopyFromVTrack(t);
    if(etp.PropagateToDCA(vtxAODNew,aod->GetMagneticField(),3.,dz,covdz)) {
      fd0[i]    = dz[0];
      fd0err[i] = TMath::Sqrt(covdz[0]);
    }
  }

  return;
}
//-----------------------------------------------------------------------------------
void AliAODRecoDecayHF::Misalign(TString misal) {
  //
  // Method to smear the impact parameter of the duaghter tracks
  // and the sec. vtx position accordinlgy 
  // Useful to study effect of misalignment.
  // The starting point are parameterisations of the impact parameter resolution
  // from MC and data 
  // Errors on d0 and vtx are not recalculated (to be done)
  //
  if(misal=="null")return;
  Double_t pard0rphiMC[3]={36.7,36.,1.25};// d0(pt)=[0]+[1]/(pt^[2]); in micron, conversion to cm is done later
  Double_t pard0rphimisal[3]={0,0,0};
  Double_t pard0zMC[3]={85.,130.,0.7};// d0(pt)=[0]+[1]/(pt^[2]); in micron, conversion to cm is done later
  Double_t pard0zmisal[3]={0,0,0};
  if(misal=="data") {
    //use this to reproduce data d0(pt) trend for pions
    pard0rphimisal[0]=37.;
    pard0rphimisal[1]=37.5;
    pard0rphimisal[2]=1.25;
    pard0zmisal[0]=96.;
    pard0zmisal[1]=131.;
    pard0zmisal[2]=0.7;
  }
  else if(misal=="resB") {
    // temporary values: asymptotic value larger by a factor 1.2 w.r.t. MC
    pard0rphimisal[0]=44.4;
    pard0rphimisal[1]=37.5;
    pard0rphimisal[2]=1.25;
    pard0zmisal[0]=115.2;
    pard0zmisal[1]=131.;
    pard0zmisal[2]=0.7;
  }
  else if(misal=="resC") {
    // temporary values: slightly larger asymptotic value, larger values at low pt
    pard0rphimisal[0]=40.;
    pard0rphimisal[1]=40.;
    pard0rphimisal[2]=1.3;
    pard0zmisal[0]=125.;
    pard0zmisal[1]=131.;
    pard0zmisal[2]=0.85;
  }
  else printf("AliAODRecoDecayHF::Misalign():  wrong misalign type specified \n");
 

  AliAODVertex *evVtx=0x0,*secVtx=0x0;
  Double_t evVtxPos[3]={-9999.,-9999.,-9999.},secVtxPos[3]={-9999.,9999.,9999.};
  if(fOwnPrimaryVtx)fOwnPrimaryVtx->GetXYZ(evVtxPos);
  else {
    evVtx=(AliAODVertex*)(fEventPrimaryVtx.GetObject());
    evVtx->GetXYZ(evVtxPos);
  }
  secVtx=(AliAODVertex*)GetSecondaryVtx();
  secVtx->GetXYZ(secVtxPos);
  
  TVector3 v2v1(secVtxPos[0]-evVtxPos[0],secVtxPos[1]-evVtxPos[1],0.);

  Double_t sigmarphinull,sigmarphimisal,sigmarphiadd;
  Double_t sigmaznull,sigmazmisal,sigmazadd;
  Double_t deltad0rphi[10],deltad0z[10];
  
  // loop on the two prongs
  for(Int_t i=0; i<fNProngs; i++) { 
    sigmarphinull = pard0rphiMC[0]+pard0rphiMC[1]/TMath::Power(PtProng(i),pard0rphiMC[2]);
    sigmarphimisal = pard0rphimisal[0]+pard0rphimisal[1]/TMath::Power(PtProng(i),pard0rphimisal[2]);
    if(sigmarphimisal>sigmarphinull) {
      sigmarphiadd = TMath::Sqrt(sigmarphimisal*sigmarphimisal-
				 sigmarphinull*sigmarphinull);
      deltad0rphi[i] = gRandom->Gaus(0.,sigmarphiadd);
    } else {
      deltad0rphi[i] = 0.;
    }

    sigmaznull =  pard0zMC[0]+pard0zMC[1]/TMath::Power(PtProng(i),pard0zMC[2]);
    sigmazmisal = pard0zmisal[0]+pard0zmisal[1]/TMath::Power(PtProng(i),pard0zmisal[2]);
    if(sigmazmisal>sigmaznull) {
      sigmazadd = TMath::Sqrt(sigmazmisal*sigmazmisal-
			      sigmaznull*sigmaznull);
      deltad0z[i] = gRandom->Gaus(0.,sigmazadd);
    } else {
      deltad0z[i] = 0.;
    }

    TVector3 pxy(fPx[i],fPy[i],0.);
    TVector3 pxycrossv2v1=pxy.Cross(v2v1);
    if( pxycrossv2v1.Z()*fd0[i] > 0 ) {
      secVtxPos[0]+=1.e-4*deltad0rphi[i]*(-fPy[i])/PtProng(i);// e-4: conversion to cm
      secVtxPos[1]+=1.e-4*deltad0rphi[i]*(+fPx[i])/PtProng(i);    
    } else {
      secVtxPos[0]+=1.e-4*deltad0rphi[i]*(+fPy[i])/PtProng(i);
      secVtxPos[1]+=1.e-4*deltad0rphi[i]*(-fPx[i])/PtProng(i);    
    }
    
    // change d0rphi
    fd0[i] += 1.e-4*deltad0rphi[i]; // e-4: conversion to cm
    // change secondary vertex z
    secVtxPos[2]+=0.5e-4*deltad0z[i];
  }
  secVtx->SetX(secVtxPos[0]);
  secVtx->SetY(secVtxPos[1]);
  secVtx->SetZ(secVtxPos[2]);

  return;
}
