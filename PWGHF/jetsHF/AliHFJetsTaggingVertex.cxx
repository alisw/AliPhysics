/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// ******************************************
// Manager class for HF jet analysis   
// Author: andrea.rossi@cern.ch, elena.bruna@cern.ch,min.jung.kweon@cern.ch,linus.feldkamp@cern.ch
// *******************************************

#include "TClonesArray.h"
#include "Riostream.h"
#include "AliAODJet.h"
#include "AliEmcalJet.h"
#include "AliVParticle.h"
#include "AliKFVertex.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliRDHFJetsCuts.h"
#include "AliRDHFJetsCutsVertex.h"
#include "AliAODVertex.h"
#include "AliAODEvent.h"
#include "AliVertexerTracks.h"
#include "AliVTrack.h"
#include "AliHFJetsTagging.h"
#include "AliHFJetsTaggingVertex.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliPicoTrack.h"


ClassImp(AliHFJetsTaggingVertex)

using std::cout;
using std::endl;

AliHFJetsTaggingVertex::AliHFJetsTaggingVertex():
AliHFJetsTagging(),
  fCutsHFjets(0x0){

  fTrackArray=new TObjArray();

}

AliHFJetsTaggingVertex::AliHFJetsTaggingVertex(const char* name):
  AliHFJetsTagging(name),
  fCutsHFjets(0x0){

  fTrackArray=new TObjArray();

}

AliHFJetsTaggingVertex::~AliHFJetsTaggingVertex(){
  if(fTrackArray) {
    delete fTrackArray;
    fTrackArray=0;
  }

  if(fCutsHFjets) {
    delete fCutsHFjets;
    //fCutsHFjets=0;
  }

}

AliHFJetsTaggingVertex &AliHFJetsTaggingVertex::operator=(const AliHFJetsTaggingVertex &c)
{
  // assigment operator

  if (this != &c)
    ((AliHFJetsTaggingVertex &) c).Copy(*this);

  return *this;
}

//___________________________________________________________________________
Int_t AliHFJetsTaggingVertex::FindVertices(const AliEmcalJet *jet, AliAODTrack **fAODTrackInfoP, AliAODTrack **fAODTrackInfoN, TClonesArray *fTrackArrayIn, AliAODEvent* aod, AliESDVertex* v1, Double_t magzkG , TClonesArray *arrVertices, Double_t *arrDispersion){

  AliInfo(MSGINFO("+++ Executing FindVertices +++"));
  
  Int_t nprong=fCutsHFjets->GetNprongs();

  Int_t nvert=0;
  if (!fCutsHFjets->IsJetSelected(jet)) {
    AliDebug(AliLog::kDebug,Form(MSGDEBUG("--> Jet not selected in FindVertices, pt=%f, eta=%f"), jet->Pt(),jet->Eta()));
    return -1;
  }

  Int_t iVerticesHF=0;

  arrVertices->Clear();
  TClonesArray &verticesHF = *arrVertices;
  Double_t dispersion=0.;
  for(Int_t ii=0;ii<5000;ii++)arrDispersion[ii]=-999;

  // TRefArray* reftracks=(TRefArray*)jet->GetRefTracks();
  fTrackArray->Clear();
  AliESDtrackCuts *esdtrcuts=fCutsHFjets->GetTrackCuts();
  // Int_t ntrks=reftracks->GetEntriesFast();
  Int_t ntrks=jet->GetNumberOfTracks();
  cout<<"ntrks="<<ntrks<<endl;
  if(ntrks<2){
    return 0;
  }
  Int_t up  = 0;
  Int_t up2 = 0;
  if (nprong == 2) {
    up  = ntrks-1;
    up2 = ntrks;
  }
  else if (nprong == 3) {
    up  = ntrks-2;
    up2 = ntrks-1;
  }
  // Int_t reftracks= jet->GetNumberOfTracks();
  //make array of ESD tracks, then needed for fTrackArray
  TObjArray* tarresd=new TObjArray();
  tarresd->SetOwner(kTRUE);

  for(Int_t j=0; j < ntrks; j++) {

    AliVTrack* trkjet=((AliPicoTrack*)jet->TrackAt(j,fTrackArrayIn))->GetTrack();
	      
    AliAODTrack *tmpTr;
    if (trkjet->GetID()>-1) tmpTr=fAODTrackInfoP[trkjet->GetID()];
    else tmpTr=fAODTrackInfoN[TMath::Abs(trkjet->GetID())];
    // const AliAODTrack *nTrack=fGTI[v0->GetNegID()];
    AliESDtrack *esdt = new AliESDtrack(tmpTr);

    Double_t point[16]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,999.,-999.,-999.,-999.,-999.,-999.};
    for(Int_t jj=0;jj<10;jj++){
      if(tmpTr->TestFilterBit(TMath::Power(2,jj))){
    	point[jj]=1;      
      }
    }
    Printf(MSGINFO("\n \n %d **** filterbit =  %d"), j, tmpTr->GetFilterMap());
    for (Int_t ii = 0; ii<10; ii++)
      Printf("%d       bit = %.1f", ii, point[ii]);
    
    tarresd->Add(esdt);
  }

  for (Int_t it=0; it < up; it++) {
    AliVTrack   *trkjet = ((AliPicoTrack *) jet->TrackAt(it,fTrackArrayIn))->GetTrack();
    AliAODTrack *tmpTr;
    if (trkjet->GetID()>-1) tmpTr=fAODTrackInfoP[trkjet->GetID()];
    else tmpTr=fAODTrackInfoN[TMath::Abs(trkjet->GetID())];
    // AliAODTrack* tmpTr = (AliAODTrack*)(jet->GetRefTracks()->At(it));
    //    Int_t id=(Int_t)TMath::Abs(tmpTr->GetID());

    if(!fCutsHFjets->IsDaughterSelected(tmpTr,v1,esdtrcuts))continue;

    //AliESDtrack *esdt1 = new AliESDtrack(tmpTr);
    AliESDtrack *esdt1 =(AliESDtrack*)tarresd->At(it);

    if(nprong<2){
      //cannot find vertices w/ less than two tracks

    }
    else
      for(Int_t it2=it+1;it2<up2;it2++){
	AliVTrack* trkjet2=((AliPicoTrack*)jet->TrackAt(it2,fTrackArrayIn))->GetTrack();
	AliAODTrack* tmpTr2;

	if (trkjet2->GetID()>-1) tmpTr2=fAODTrackInfoP[trkjet2->GetID()];
	else tmpTr2=fAODTrackInfoN[TMath::Abs(trkjet2->GetID())];
	// = (AliAODTrack*)(jet->GetRefTracks()->At(it2));
        if(!fCutsHFjets->IsDaughterSelected(tmpTr2,v1,esdtrcuts)) continue;

        AliESDtrack *esdt2 =(AliESDtrack*)tarresd->At(it2);

        //fill TClonesArray of tracks
        fTrackArray->AddAt(esdt1,0);
        fTrackArray->AddAt(esdt2,1);

        if(nprong==2) {

          AliAODVertex* vert = ReconstructSecondaryVertex(fTrackArray,v1,magzkG,dispersion);

          //printf("vertex done \n");
          if(vert){
            vert->AddDaughter(tmpTr);
            vert->AddDaughter(tmpTr2);
            if(!fCutsHFjets->IsVertexSelected(vert,aod,magzkG,dispersion))continue;
            nvert++;
         
            new(verticesHF[iVerticesHF])AliAODVertex(*vert);
            arrDispersion[iVerticesHF]=dispersion;
            iVerticesHF++;
	    
            fTrackArray->Clear();
          }
        }
        if (nprong >= 3) {
          for(Int_t it3=it2+1;it3<ntrks;it3++){
            AliVTrack* trkjet3=((AliPicoTrack*)jet->TrackAt(it3,fTrackArrayIn))->GetTrack();
            AliAODTrack* tmpTr3;
            if (trkjet3->GetID()>-1) tmpTr3=fAODTrackInfoP[trkjet3->GetID()];
            else tmpTr3=fAODTrackInfoN[TMath::Abs(trkjet3->GetID())];
            if(!fCutsHFjets->IsDaughterSelected(tmpTr3,v1,esdtrcuts)) continue;
            AliESDtrack *esdt3 =(AliESDtrack*)tarresd->At(it3);
            
            fTrackArray->AddAt(esdt3,2);
            AliAODVertex* vert = ReconstructSecondaryVertex(fTrackArray,v1,magzkG,dispersion);
            if(vert){
              vert->AddDaughter(tmpTr);
              vert->AddDaughter(tmpTr2);
              vert->AddDaughter(tmpTr3);
              if(!fCutsHFjets->IsVertexSelected(vert,aod,magzkG,dispersion))continue;
              nvert++;
              //AliAODVertex* vert3prong= new(verticesHF[iVerticesHF])AliAODVertex(*vert);
              new(verticesHF[iVerticesHF])AliAODVertex(*vert);
              arrDispersion[iVerticesHF]=dispersion;
              // cout<<"================================ vertex 3 prong -> "<<vert->GetX()<<"  "<<vert->GetY()<<"  "<<vert->GetZ()<<"  "<<endl;
              iVerticesHF++;
              
            }	    

          }

        }

	if (nprong > 3) {

          //only 2- and 3- track vertices at the moment
    
        }

      }

  }
  fTrackArray->Clear();

  tarresd->Clear();
  delete tarresd;
  return nvert;
}

//-----------------------------------------------------------------------------
AliAODVertex* AliHFJetsTaggingVertex::ReconstructSecondaryVertex(TObjArray *trkArray, AliESDVertex* v1, Double_t magzkG ,Double_t &dispersion) const
{
  //cout<<"in AliHFJetsTaggingVertex::ReconstructSecondaryVertex"<<endl;

  // Secondary vertex reconstruction with AliVertexerTracks or AliKFParticle
  //AliCodeTimerAuto("",0);

  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;

  AliVertexerTracks* vertexerTracks=new AliVertexerTracks(magzkG);
  Int_t nTrks = trkArray->GetEntriesFast();

  Int_t secVtxWithKF=fCutsHFjets->GetSecVtxWithKF();
  if(!secVtxWithKF) { // AliVertexerTracks

    vertexerTracks->SetVtxStart(v1);//--> need primary vertex!!!

    vertexESD = (AliESDVertex*)vertexerTracks->VertexForSelectedESDTracks(trkArray);


    if(!vertexESD) {
      // printf("CANNOT BUILD THE VTX \n");
      return vertexAOD;
    }
    //    vertexESD->Print();
    
    if(vertexESD->GetNContributors()!=trkArray->GetEntriesFast()) {

      //AliDebug(2,"vertexing failed"); 
      delete vertexESD; vertexESD=NULL;
      return vertexAOD;
    }

    Double_t vertRadius2=vertexESD->GetX()*vertexESD->GetX()+vertexESD->GetY()*vertexESD->GetY();
    if(vertRadius2>8.){
      // printf(" vertex outside beam pipe (%f), reject candidate to avoid propagation through material",vertRadius2);
      delete vertexESD; vertexESD=NULL;
      return vertexAOD;
    }

  } else { // Kalman Filter vertexer (AliKFParticle)

    AliKFParticle::SetField(magzkG); //-->  need magnetic field!!!!

    AliKFVertex vertexKF;

    //Int_t nTrks = trkArray->GetEntriesFast();
    nTrks = trkArray->GetEntriesFast();

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
  dispersion = vertexESD->GetDispersion();
  delete vertexESD; vertexESD=NULL;

  Int_t nprongs= trkArray->GetEntriesFast();
  //cout<<"in AliHFJetsTaggingVertex::ReconstructSecondaryVertex --> nprongs = "<<nprongs<<endl;
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);
  //cout<<"in AliHFJetsTaggingVertex::ReconstructSecondaryVertex --> vertexAOD = "<<vertexAOD<<endl;
  //  vertexAOD->Print();
  delete vertexerTracks;
  return vertexAOD;

}

Double_t AliHFJetsTaggingVertex::GetVertexInvariantMass(AliAODVertex *vtx,Double_t massParticle){
  Double_t pxyz[3];
  Double_t pxyzSum[4]={0.,0.,0.,0.};

  for(Int_t jp=0;jp<vtx->GetNDaughters();jp++){
    AliAODTrack *tr=(AliAODTrack*)vtx->GetDaughter(jp);

    tr->GetPxPyPz(pxyz);
    pxyzSum[1]+=pxyz[0];
    pxyzSum[2]+=pxyz[1];
    pxyzSum[3]+=pxyz[2];
    pxyzSum[0]+=TMath::Sqrt(massParticle*massParticle+pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[2]*pxyz[2]);//pion mass assumed
  }

  return TMath::Sqrt(pxyzSum[0]*pxyzSum[0]-pxyzSum[1]*pxyzSum[1]-pxyzSum[2]*pxyzSum[2]-pxyzSum[3]*pxyzSum[3]);
}

void AliHFJetsTaggingVertex::GetVtxPxy(AliAODVertex *vtx,Double_t *pxyzSum){
  Double_t pxyz[3];
  pxyzSum[0]=0;
  pxyzSum[1]=0;
  pxyzSum[2]=0;
  for(Int_t jp=0;jp<vtx->GetNDaughters();jp++){
    AliAODTrack *tr=(AliAODTrack*)vtx->GetDaughter(jp);
 
    tr->GetPxPyPz(pxyz);
    pxyzSum[0]+=pxyz[0];
    pxyzSum[1]+=pxyz[1];
    pxyzSum[2]+=pxyz[2];
  }

  return;
}

Int_t AliHFJetsTaggingVertex::GetNTracksFromCommonVertex(AliAODVertex *vtx,const TClonesArray *mcPart,Int_t &mcVert,Double_t &xVtxMC,Double_t &yVtxMC,Int_t &nfromBandD,Int_t &nfromD,Int_t &nfromPromptD){
  Int_t label;

  Double_t *x=new Double_t[vtx->GetNDaughters()];
  Double_t *y=new Double_t[vtx->GetNDaughters()];
  Double_t *z=new Double_t[vtx->GetNDaughters()];
  Int_t *ntrks=new Int_t[vtx->GetNDaughters()];
  Int_t *fromB=new Int_t[vtx->GetNDaughters()];
  Int_t *fromD=new Int_t[vtx->GetNDaughters()];
  Int_t *fromPromptD=new Int_t[vtx->GetNDaughters()];
  Int_t *fromDfromB=new Int_t[vtx->GetNDaughters()];
  Int_t *fromSameB=new Int_t[vtx->GetNDaughters()];
  Int_t *fromDfromSameB=new Int_t[vtx->GetNDaughters()];
  Int_t *labmoth=new Int_t[vtx->GetNDaughters()];
  Int_t *labmothD=new Int_t[vtx->GetNDaughters()];
  Int_t *labgrmoth=new Int_t[vtx->GetNDaughters()];
  Int_t *labgrmothD=new Int_t[vtx->GetNDaughters()];
  Int_t *labgrgrmothD=new Int_t[vtx->GetNDaughters()];
  Int_t *labgrgrmothB=new Int_t[vtx->GetNDaughters()];
  Int_t *labgrgrmothB2=new Int_t[vtx->GetNDaughters()];
  mcVert=0;
  Int_t maxTrks=0;
  Bool_t vtxfound=kFALSE;
  Bool_t vtxPromptDfound=kFALSE;
  Bool_t vtxBDfound=kFALSE;
  Int_t mcBDVert=-1;
  Int_t mcPromptDVert=-1;
  for(Int_t jp=0;jp<vtx->GetNDaughters();jp++){
    fromB[jp]=0;
    fromD[jp]=0;
    fromPromptD[jp]=0;
    fromDfromB[jp]=0;
    fromSameB[jp]=0;
    fromDfromSameB[jp]=0;
    labmoth[jp]=0;
    labgrmoth[jp]=0;
    labmothD[jp]=0;
    labgrmothD[jp]=0;
    labgrgrmothD[jp]=0;
    labgrgrmothB[jp]=0;
    labgrgrmothB2[jp]=0;

  }

  for(Int_t jp=0;jp<vtx->GetNDaughters();jp++){
    ntrks[jp]=0;

    vtxfound=kFALSE;
    vtxBDfound=kFALSE;
    AliAODTrack *tr=(AliAODTrack*)vtx->GetDaughter(jp);
    label=TMath::Abs(tr->GetLabel());
    if(label<1){

      // Printf("Vertex with proton from beam ??? \n");
      continue;
    }

    AliAODMCParticle *mcp=(AliAODMCParticle*)mcPart->At(label);
    if(!mcp)continue;

    labmoth[jp]=mcp->GetMother();

    AliAODMCParticle* moth=(AliAODMCParticle*)mcPart->At(labmoth[jp]);
    if(!moth)continue;

    if(TMath::Abs(moth->GetPdgCode()==5) || (TMath::Abs(moth->GetPdgCode())>500 && TMath::Abs(moth->GetPdgCode())<600) || (TMath::Abs(moth->GetPdgCode())>5000 && TMath::Abs(moth->GetPdgCode())<6000)){
      fromB[jp]=1;
      // cout<<"B found -> moth="<<moth->GetPdgCode()<<"  labmoth="<<labmoth[jp]<<endl;
    }

    if((TMath::Abs(moth->GetPdgCode())>400 && TMath::Abs(moth->GetPdgCode())<500)|| (TMath::Abs(moth->GetPdgCode())>4000 && TMath::Abs(moth->GetPdgCode())<5000)){
      labmothD[jp]=labmoth[jp];
      fromD[jp]=1;
      // cout<<"mother D found -> mothD="<<moth->GetPdgCode()<<"  labmothD="<<labmoth[jp]<<endl;
    }

    label=moth->GetMother();
    AliAODMCParticle* grmothD=(AliAODMCParticle*)mcPart->At(label);
    if(grmothD){
      if((TMath::Abs(grmothD->GetPdgCode())>400 && TMath::Abs(grmothD->GetPdgCode())<500) || (TMath::Abs(grmothD->GetPdgCode())>4000 && TMath::Abs(grmothD->GetPdgCode())<4000)){
        fromD[jp]=1;
        labgrmothD[jp]=label;
        // cout<<jp<<"  grandmother D found -> moth="<<moth->GetPdgCode()<<"  grmoth="<<grmothD->GetPdgCode()<<"  labgrmothD="<<labgrmothD[jp]<<endl;
      }

      Int_t label2=grmothD->GetMother();
      AliAODMCParticle* grgrmothD=(AliAODMCParticle*)mcPart->At(label2);
      if(grgrmothD){
        if((TMath::Abs(grgrmothD->GetPdgCode())>400 && TMath::Abs(grgrmothD->GetPdgCode())<500) || (TMath::Abs(grgrmothD->GetPdgCode())>4000 && TMath::Abs(grgrmothD->GetPdgCode())<4000)){
          fromD[jp]=1;
          labgrgrmothD[jp]=label2;
          // cout<<jp<<"  -->grandgrandmother D found -> moth="<<moth->GetPdgCode()<<"  grgrmoth="<<grgrmothD->GetPdgCode()<<"  labgrgrmothD="<<labgrgrmothD[jp]<<endl;

        }
      }
    }
    if((TMath::Abs(moth->GetPdgCode())>400 && TMath::Abs(moth->GetPdgCode())<500) || (TMath::Abs(moth->GetPdgCode())>4000 && TMath::Abs(moth->GetPdgCode())<5000)){
      labgrmoth[jp]=moth->GetMother();
      AliAODMCParticle* grmoth=(AliAODMCParticle*)mcPart->At(labgrmoth[jp]);      if(grmoth){
        if(TMath::Abs(grmoth->GetPdgCode()==5) ||(TMath::Abs(grmoth->GetPdgCode())>500 && TMath::Abs(grmoth->GetPdgCode())<600) || (TMath::Abs(grmoth->GetPdgCode())>5000 && TMath::Abs(grmoth->GetPdgCode())<6000)){
          fromDfromB[jp]=1;
          fromD[jp]=0;
          // cout<<"D found -> moth="<<moth->GetPdgCode()<<"  grmoth="<<grmoth->GetPdgCode()<<"  labgrmoth="<<labgrmoth[jp]<<endl;
        }

        labgrgrmothB[jp]=grmoth->GetMother();
        AliAODMCParticle* grgrmoth=(AliAODMCParticle*)mcPart->At(labgrgrmothB[jp]);
        if(grgrmoth){
          if((TMath::Abs(grgrmoth->GetPdgCode()==5) || (TMath::Abs(grgrmoth->GetPdgCode())>500 && TMath::Abs(grgrmoth->GetPdgCode())<600)) || ((TMath::Abs(grgrmoth->GetPdgCode())>5000 && TMath::Abs(grgrmoth->GetPdgCode())<6000))){
            fromDfromB[jp]=1;
            fromD[jp]=0;
            // cout<<"** D found -> moth="<<moth->GetPdgCode()<<"  gratgrmoth="<<grgrmoth->GetPdgCode()<<"  labgreatgrmoth="<<labgrgrmothB[jp]<<endl;
          }
          labgrgrmothB2[jp]=grgrmoth->GetMother();
          AliAODMCParticle* grgrmoth2=(AliAODMCParticle*)mcPart->At(labgrgrmothB2[jp]);
          if(grgrmoth2){
            if((TMath::Abs(grgrmoth2->GetPdgCode()==5) || (TMath::Abs(grgrmoth2->GetPdgCode())>500 && TMath::Abs(grgrmoth2->GetPdgCode())<600)) || ((TMath::Abs(grgrmoth2->GetPdgCode())>5000 && TMath::Abs(grgrmoth2->GetPdgCode())<6000))){
              fromDfromB[jp]=1;
              fromD[jp]=0;
            }

          }
        }
      }
    }
    for(Int_t jv=0;jv<mcVert;jv++){
      if(TMath::Abs(x[jv]-mcp->Xv())+TMath::Abs(y[jv]-mcp->Yv())+TMath::Abs(z[jv]-mcp->Zv())<0.0001){
        ntrks[jv]++;
        vtxfound=kTRUE;
        break;
      }

    }
    for(Int_t jv=0;jv<jp;jv++){
      if ((fromB[jv]==fromB[jp] && fromB[jp]==1 && labmoth[jv]==labmoth[jp])||(fromDfromB[jv]==1 && fromB[jp]==1 && labgrmoth[jv]==labmoth[jp])){
        fromSameB[mcBDVert]++;
        vtxBDfound=kTRUE;
        break;
      }

      if ((fromB[jv]==1 && fromDfromB[jp]==1 && labmoth[jv]==labgrmoth[jp]) || (fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jp]==1 && labgrmoth[jv]==labgrmoth[jp]) || (fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jp]==1 && labgrmoth[jv]==labgrgrmothB[jp]) || (fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jp]==1 && labgrgrmothB[jv]==labgrmoth[jp]) || (fromB[jv]==1 && fromDfromB[jp]==1 && labmoth[jv]==labgrgrmothB[jp]) || (fromDfromB[jv]==1  && fromB[jp]==1 && labgrgrmothB[jv]==labmoth[jp])   ||  (fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jp]==1 && labgrmoth[jv]==labgrgrmothB2[jp]) || (fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jp]==1 && labgrgrmothB2[jv]==labgrmoth[jp]) || (fromB[jv]==1 && fromDfromB[jp]==1 && labmoth[jv]==labgrgrmothB2[jp]) || (fromDfromB[jv]==1  && fromB[jp]==1 && labgrgrmothB2[jv]==labmoth[jp])){

        fromDfromSameB[mcBDVert]++;
        vtxBDfound=kTRUE;

        break;
      }
    }

    if(!vtxBDfound){
      if(fromB[jp]==1 || fromDfromB[jp]==1){
        mcBDVert++;
        if(fromB[jp]==1)fromSameB[mcBDVert]++;
        if(fromDfromB[jp]==1)fromDfromSameB[mcBDVert]++;
      }
    }

    for(Int_t jv=0;jv<jp;jv++){
      if(fromB[jv]==fromB[jp] && fromB[jp]==0 && fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jv]==0 && fromD[jv]==fromD[jp] && fromD[jp]==1){
        if(labgrmothD[jv]==labmothD[jp] || labmothD[jv]==labgrmothD[jp] || labmothD[jv]==labmothD[jp] || labgrgrmothD[jv]==labmothD[jp] || labgrgrmothD[jp]==labmothD[jv] ||  labgrgrmothD[jv]==labgrmothD[jp] || labgrgrmothD[jp]==labgrmothD[jv]){
          fromPromptD[mcPromptDVert]++;
          vtxPromptDfound=kTRUE;
          break;
        }
      }
    }


    if(!vtxfound){
      x[mcVert]=mcp->Xv();
      y[mcVert]=mcp->Yv();
      z[mcVert]=mcp->Zv();

      mcVert++;
    }
    if(!vtxPromptDfound){
      if(fromD[jp]==1 &&fromB[jp]==0 && fromDfromB[jp]==0){

        mcPromptDVert++;
        fromPromptD[mcPromptDVert]++;


      }
    }


  }
  for(Int_t jv=0;jv<mcBDVert+1;jv++){
    AliDebug(AliLog::kDebug,Form("vertexBD %d, ntracksB=%d, ntracksD=%d",jv,fromSameB[jv],fromDfromSameB[jv]));
  }

  for(Int_t jv=0;jv<mcPromptDVert+1;jv++){
    AliDebug(AliLog::kDebug,Form("vertexPromptD %d, ntracks=%d",jv,fromPromptD[jv]));
  }

  for(Int_t jv=0;jv<mcVert;jv++){
    if(ntrks[jv]>maxTrks){
      maxTrks=ntrks[jv];
      xVtxMC=x[jv];
      yVtxMC=y[jv];
    }
  }

  Int_t max=0;
  for(Int_t jv=0;jv<mcBDVert+1;jv++){
    Int_t sum=fromSameB[jv]+fromDfromSameB[jv];
    if (sum>max){
      max=sum;
      nfromBandD=max;
      nfromD=fromDfromSameB[jv];
    }
  }

  Int_t maxD=0;
  for(Int_t jv=0;jv<mcPromptDVert+1;jv++){

    if (fromPromptD[jv]>maxD){
      maxD=fromPromptD[jv];
    }
  }
  nfromPromptD=maxD;
  delete x;
  delete y;
  delete z;
  delete ntrks;
  delete fromB;
  delete fromDfromB;
  delete fromSameB;
  delete fromDfromSameB;
  delete labmoth;
  delete labgrmoth;
  return maxTrks;
}
