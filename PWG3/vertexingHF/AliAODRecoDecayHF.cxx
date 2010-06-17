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

/////////////////////////////////////////////////////////////
//
// Base class for AOD reconstructed heavy-flavour decay
//
// Author: A.Dainese, andrea.dainese@lnl.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <TVector3.h>
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
  fProngID(0x0) 
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
  fProngID(0x0) 
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
  fProngID(0x0) 
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
  fProngID(0x0) 
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
  fProngID(0x0)
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

  if(source.GetOwnPrimaryVtx()) fOwnPrimaryVtx = new AliAODVertex(*(source.GetOwnPrimaryVtx()));

  if(source.GetNProngs()>0) {
    fd0err = new Double_t[GetNProngs()];
    memcpy(fd0err,source.fd0err,GetNProngs()*sizeof(Double_t));
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
  vertexer->SetMinClusters(4);
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

  // now recalculate the daughters impact parameters
  AliExternalTrackParam *etp = 0;
  Double_t dz[2],covdz[3];
  for(Int_t i=0; i<ndg; i++) {
    t = (AliAODTrack*)GetDaughter(i);
    etp = new AliExternalTrackParam(t);
    if(etp->PropagateToDCA(vtxAODNew,aod->GetMagneticField(),3.,dz,covdz)) {
      fd0[i]    = dz[0];
      fd0err[i] = TMath::Sqrt(covdz[0]);
    }
    delete etp; etp=NULL;
  }

  return vtxAODNew;
}
