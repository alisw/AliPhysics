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
// Modified by: ycorrale@cern.ch
// *******************************************

//--Root--
#include <TClonesArray.h>

//--AliRoot--
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliEmcalJet.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliKFVertex.h"
#include "AliPicoTrack.h"
#include "AliVertexerTracks.h"

//--AliHFJetsClass--
#include "AliHFJetsTaggingVertex.h"

//_____________________________________________________________________________________

ClassImp(AliHFJetsTaggingVertex)

//_____________________________________________________________________________________
AliHFJetsTaggingVertex::AliHFJetsTaggingVertex() : AliHFJetsTagging(), fCutsHFjets(NULL) {

  fTrackArray = new TObjArray();
}

//_____________________________________________________________________________________
AliHFJetsTaggingVertex::AliHFJetsTaggingVertex(const char *name) : AliHFJetsTagging(name),
fCutsHFjets(NULL) {
  
  fTrackArray = new TObjArray();
}

//_____________________________________________________________________________________
AliHFJetsTaggingVertex::~AliHFJetsTaggingVertex() {
  
  if (fTrackArray) {
    delete fTrackArray; fTrackArray = NULL;
  }
  
  if (fCutsHFjets) {
    delete fCutsHFjets; fCutsHFjets = NULL;
  }
}

//_____________________________________________________________________________________
AliHFJetsTaggingVertex &AliHFJetsTaggingVertex::operator=(const AliHFJetsTaggingVertex &c)
{
  // assigment operator
  if ( this != &c )
    ( (AliHFJetsTaggingVertex &) c).Copy( *this );

  return *this;
}

//_____________________________________________________________________________________
Int_t AliHFJetsTaggingVertex::FindVertices(const AliEmcalJet *jet,
                                           map_AliAODTrk     *fAODTrackInfo,
                                           TClonesArray      *fTrackArrayIn,
                                           AliAODEvent       *aod,
                                           AliESDVertex      *v1,
                                           Double_t           magzkG,
                                           TClonesArray      *arrVtxHF,
                                           vctr_pair_dbl_int &arrVtxDisp)
{
  AliInfo(MSGINFO("+++ Executing FindVertices +++"));
  
  if (!fCutsHFjets->IsJetSelected(jet)) {
   
    AliDebug(AliLog::kDebug, Form(MSGDEBUG("--> Jet not selected in FindVertices, pt=%f, eta=%f"),
                                  jet->Pt(),
                                  jet->Eta()));
    return -1;
  }
  
  Int_t nSecndVxtHF = 0;

  arrVtxHF->Clear();
  arrVtxDisp.clear();
  
  Double_t vtxRes = 0.;
  AliESDtrackCuts *esdtrcuts = fCutsHFjets->GetTrackCuts();
  
  Int_t nTrksInJet = jet->GetNumberOfTracks();
  AliInfoF(MSGINFO("nTrksInJet = %d \n"), nTrksInJet);
  
  if (nTrksInJet < 2) {
    AliWarning(MSGWARNING("Cannot find vertices w/ only one track"));
    return 0;
  } 
 

  Int_t nProngTrack = fCutsHFjets->GetNprongs();
  if ( nProngTrack < 2 || nProngTrack > 3 ) {
    AliWarning(MSGWARNING("Cannot find vertices w/ less(more) than two(three) tracks"));
    return 0;
  }
  //make array of ESD tracks, then needed for fTrackArray
  typedef vector< pair<Int_t, AliESDtrack *> > vctr_pair_int_esdTrk;
  vctr_pair_int_esdTrk arrESDtrkInfo;
  arrESDtrkInfo.reserve(nTrksInJet);
  
  for (Int_t j = 0; j < nTrksInJet; ++j) {
    
    AliAODTrack *jTrk   = ((AliAODTrack*)jet->TrackAt(j, fTrackArrayIn));
    //utilize dynamic cast and then check pointer
    Int_t jTrkID = jTrk->GetID();
    
    AliInfoF(MSGINFO("Track index  %d"), jTrkID);
    if (jTrkID < 0) {

      AliInfoF(MSGINFO("Track with index < 0 %d"), jTrkID);
      continue;
    }
   
    
    // AliAODTrack *tmpAODtrk = (* fAODTrackInfo)[jTrkID].first;
    // printf("tmpAODtrk =%x \n",tmpAODtrk);
 
    // if (!tmpAODtrk) {
  
    //   AliWarning(MSGWARNING("AliPicoTrack %d not found on AOD event. Please check the physics selection for Emcal"));
    //   continue;
    // }
    
    if (!fCutsHFjets->IsDaughterSelected(jTrk, v1, esdtrcuts)) continue;
    
    AliESDtrack *tmpESDtrk = new AliESDtrack(jTrk);
    
    arrESDtrkInfo.push_back(make_pair(jTrkID, tmpESDtrk));
  }
  
  Int_t nGoodTrks = (Int_t)arrESDtrkInfo.size();
  AliInfoF(MSGINFO("Number of good tracks = %d"), nGoodTrks);
  
  if (nGoodTrks < 2) {
    AliWarning(MSGWARNING("Cannot find vertices w/ only one good track"));
    return 0;
  }

  Int_t up = nGoodTrks - ((nProngTrack == 2) ? 1 : 2);
  
  for (Int_t it1 = 0; it1 < up; ++it1) {
    
    Int_t        jTrkID_1 = (arrESDtrkInfo.at(it1)).first;
    AliESDtrack *esdTrk_1 = (arrESDtrkInfo.at(it1)).second;
    
    fTrackArray->Clear();
    fTrackArray->AddAt(esdTrk_1, 0);

    for (Int_t it2 = it1+1; it2 < up+1; ++it2) {
      
      Int_t        jTrkID_2 = (arrESDtrkInfo.at(it2)).first;
      AliESDtrack *esdTrk_2 = (arrESDtrkInfo.at(it2)).second;

      fTrackArray->AddAt(esdTrk_2, 1);
      
      if (nProngTrack == 2) {
        
        AliAODVertex *vert = ReconstructSecondaryVertex(fTrackArray, v1, magzkG, vtxRes);
        
        if (vert) {
          
          vert->AddDaughter((* fAODTrackInfo)[jTrkID_1].first);
          vert->AddDaughter((* fAODTrackInfo)[jTrkID_2].first);
          
          if ( !fCutsHFjets->IsVertexSelected(vert, aod, magzkG, vtxRes) ) continue;
          
          
          Int_t nTrkBelongToV0 = (* fAODTrackInfo)[jTrkID_1].second +
                                 (* fAODTrackInfo)[jTrkID_2].second;
          
          new ((* arrVtxHF)[nSecndVxtHF]) AliAODVertex(* vert);
          arrVtxDisp.push_back(make_pair(vtxRes, nTrkBelongToV0));
          nSecndVxtHF++;
          
        } // end if (vert)
        
      } // end if ( nProngTrack == 2 )
      else {
        
        for (Int_t it3 = it2+1; it3 < nGoodTrks; ++it3) {
          
          Int_t        jTrkID_3 = (arrESDtrkInfo.at(it3)).first;
          AliESDtrack *esdTrk_3 = (arrESDtrkInfo.at(it3)).second;
          
          fTrackArray->AddAt(esdTrk_3, 2);
          
          AliAODVertex *vert = ReconstructSecondaryVertex(fTrackArray, v1, magzkG, vtxRes);
          
          if (vert) {
            
            vert->AddDaughter((* fAODTrackInfo)[jTrkID_1].first);
            vert->AddDaughter((* fAODTrackInfo)[jTrkID_2].first);
            vert->AddDaughter((* fAODTrackInfo)[jTrkID_3].first);

            if (!fCutsHFjets->IsVertexSelected(vert, aod, magzkG, vtxRes)) continue;

            Int_t nTrkBelongToV0 = (* fAODTrackInfo)[jTrkID_1].second +
                                   (* fAODTrackInfo)[jTrkID_2].second +
                                   (* fAODTrackInfo)[jTrkID_3].second;
            
            new ((* arrVtxHF)[nSecndVxtHF]) AliAODVertex(* vert);
            arrVtxDisp.push_back(make_pair(vtxRes, nTrkBelongToV0));
            nSecndVxtHF++;
            
          } // end if (vert)
          
        } // end for it3
        
      } // end else
      
    } // end for it2
    
  } // end for it1
  
  fTrackArray->Clear();

  for (vctr_pair_int_esdTrk::iterator it = arrESDtrkInfo.begin(); it != arrESDtrkInfo.end(); ++it) {
    delete (* it).second;
  }
  
  return nSecndVxtHF;
}

//_____________________________________________________________________________________
AliAODVertex *AliHFJetsTaggingVertex::ReconstructSecondaryVertex(TObjArray    *trkArray,
                                                                AliESDVertex *v1,
                                                                Double_t      magzkG,
                                                                Double_t     &vtxRes) const
{
  //cout<<"in AliHFJetsTaggingVertex::ReconstructSecondaryVertex"<<endl;

  // Secondary vertex reconstruction with AliVertexerTracks or AliKFParticle
  //AliCodeTimerAuto("",0);

  AliESDVertex      *vertexESD = NULL;
  AliVertexerTracks *vertexerTracks = new AliVertexerTracks(magzkG);
  
  Int_t nProngTrks   = trkArray->GetEntriesFast();
  Int_t secVtxWithKF = fCutsHFjets->GetSecVtxWithKF();
  
  if (!secVtxWithKF) { // AliVertexerTracks

    vertexerTracks->SetVtxStart(v1); //--> need primary vertex!!!

    vertexESD = (AliESDVertex *)vertexerTracks->VertexForSelectedESDTracks(trkArray);


    if (!vertexESD) {
      // printf("CANNOT BUILD THE VTX \n");
      return NULL;
    }
    //    vertexESD->Print();
    
    if(vertexESD->GetNContributors() != nProngTrks) {

      //AliDebug(2,"vertexing failed"); 
      delete vertexESD; vertexESD = NULL;
      return NULL;
    }
    
    Double_t vertRadius2 = vertexESD->GetX()*vertexESD->GetX() + vertexESD->GetY()*vertexESD->GetY();
    if (vertRadius2 > 8.) {
      // printf(" vertex outside beam pipe (%f), reject candidate to avoid propagation through material",vertRadius2);
      delete vertexESD; vertexESD = NULL;
      return NULL;
    }
  }
  else { // Kalman Filter vertexer (AliKFParticle)

    AliKFParticle::SetField(magzkG); //-->  need magnetic field!!!!

    AliKFVertex vertexKF;

    //Int_t nTrks = trkArray->GetEntriesFast();
    for (Int_t i = 0; i < nProngTrks; ++i) {
      
      AliESDtrack *esdTrack = (AliESDtrack *)trkArray->At(i);

      AliKFParticle daughterKF(*esdTrack, 211);
      vertexKF.AddDaughter(daughterKF);
    }
    vertexESD = new AliESDVertex(vertexKF.Parameters(),
                                 vertexKF.CovarianceMatrix(),
                                 vertexKF.GetChi2(),
                                 vertexKF.GetNContributors());
  }
  
  // convert to AliAODVertex
  Double_t pos[3], cov[6], chi2xNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2xNDF = vertexESD->GetChi2toNDF();
  vtxRes   = vertexESD->GetDispersion();

  delete vertexESD; vertexESD = NULL;
  delete vertexerTracks;
  
  return (new AliAODVertex(pos, cov, chi2xNDF, NULL, -1, AliAODVertex::kUndef, nProngTrks));
}

//_____________________________________________________________________________________
void AliHFJetsTaggingVertex::GetVtxPxy(AliAODVertex *vtx, Double_t *pxyzSum) {
  
  Double_t pxyz[3];
  pxyzSum[0]=0;
  pxyzSum[1]=0;
  pxyzSum[2]=0;
  
  for (Int_t jp = 0; jp < vtx->GetNDaughters(); ++jp) {
    
    AliAODTrack *tr = (AliAODTrack *)vtx->GetDaughter(jp);
    
    tr->GetPxPyPz(pxyz);
    pxyzSum[0] += pxyz[0];
    pxyzSum[1] += pxyz[1];
    pxyzSum[2] += pxyz[2];
  }
  
  return;
}

//_____________________________________________________________________________________
Double_t AliHFJetsTaggingVertex::GetVertexInvariantMass(AliAODVertex *vtx, Double_t mass)
{
  Double_t pxyz[3];
  Double_t pxyzSum[4] = {0., 0., 0., 0.};
  
  for (Int_t jp = 0; jp < vtx->GetNDaughters(); ++jp) {
    
    AliAODTrack *tr = (AliAODTrack *)vtx->GetDaughter(jp);

    tr->GetPxPyPz(pxyz);
    pxyzSum[1] += pxyz[0];
    pxyzSum[2] += pxyz[1];
    pxyzSum[3] += pxyz[2];
    pxyzSum[0] += TMath::Sqrt(mass*mass+pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]+pxyz[2]*pxyz[2]); //pion mass assumed
  }

  return TMath::Sqrt(pxyzSum[0]*pxyzSum[0]-pxyzSum[1]*pxyzSum[1]-pxyzSum[2]*pxyzSum[2]-pxyzSum[3]*pxyzSum[3]);
}


//_____________________________________________________________________________________
Int_t AliHFJetsTaggingVertex::GetNTracksFromCommonVertex(AliAODVertex       *vtx,
                                                        const TClonesArray *mcPartArray,
                                                        Int_t              &mcVert,
                                                        Double_t           &xVtxMC,
                                                        Double_t           &yVtxMC,
                                                        Int_t              &nfromBandD,
                                                        Int_t              &nfromD,
                                                        Int_t              &nfromPromptD)
{
  Int_t label;
  Int_t nProngTrk = vtx->GetNDaughters();
  
  Double_t *x        = new Double_t[nProngTrk];
  Double_t *y        = new Double_t[nProngTrk];
  Double_t *z        = new Double_t[nProngTrk];
  
  Int_t *ntrks          = new Int_t[nProngTrk];
  Int_t *fromB          = new Int_t[nProngTrk];
  Int_t *fromD          = new Int_t[nProngTrk];
  Int_t *fromPromptD    = new Int_t[nProngTrk];
  Int_t *fromDfromB     = new Int_t[nProngTrk];
  Int_t *fromSameB      = new Int_t[nProngTrk];
  Int_t *fromDfromSameB = new Int_t[nProngTrk];
  Int_t *labmoth        = new Int_t[nProngTrk];
  Int_t *labmothD       = new Int_t[nProngTrk];
  Int_t *labgrmoth      = new Int_t[nProngTrk];
  Int_t *labgrmothD     = new Int_t[nProngTrk];
  Int_t *labgrgrmothD   = new Int_t[nProngTrk];
  Int_t *labgrgrmothB   = new Int_t[nProngTrk];
  Int_t *labgrgrmothB2  = new Int_t[nProngTrk];
  
  mcVert = 0;
  Int_t maxTrks = 0;
  
  Bool_t vtxfound        = kFALSE;
  Bool_t vtxPromptDfound = kFALSE;
  Bool_t vtxBDfound      = kFALSE;
  Int_t mcBDVert=-1;
  Int_t mcPromptDVert=-1;
  
  for ( Int_t jp = 0; jp < nProngTrk; jp++ ) {
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
  
  for ( Int_t jp=0; jp < nProngTrk; jp++ ) {
    ntrks[jp] = 0;
    
    vtxfound   = kFALSE;
    vtxBDfound = kFALSE;
    
    AliAODTrack *tr = (AliAODTrack *)vtx->GetDaughter(jp);
    label = TMath::Abs(tr->GetLabel());
    if ( label < 1 ) {
      // Printf("Vertex with proton from beam ??? \n");
      continue;
    }
    
    AliAODMCParticle *mcPart = (AliAODMCParticle *)mcPartArray->At(label);
    if ( !mcPart ) continue;
    
/*    labmoth[jp] = mcPart->GetMother();
    
    AliAODMCParticle *moth = ( labmoth[jp] > -1 ) ? (AliAODMCParticle*)mcPartArray->At(labmoth[jp]) : 0;
    if ( !moth ) continue;
    Int_t mothPDG = TMath::Abs( moth->GetPdgCode() );

    label = moth->GetMother();
    AliAODMCParticle *grMothD = ( label > -1) ? (AliAODMCParticle *)mcPartArray->At(label) : 0;
    Int_t grMothPDG = TMath::Abs( grMothD->GetPdgCode() );
    
    //Check if is a B meson
    if ( ( mothPDG == 5 ) || ( mothPDG > 500 && mothPDG < 600) || ( mothPDG > 5000 && mothPDG < 6000 ) ) {
      
      fromB[jp] = 1;
      // cout<<"B found -> moth="<<moth->GetPdgCode()<<"  labmoth="<<labmoth[jp]<<endl;
    
    }
    
    if ( ( mothPDG > 400 && mothPDG < 500) || ( mothPDG > 4000 && mothPDG < 5000 ) ) {
      
      labmothD[jp] = labmoth[jp];
      fromD[jp]    = 1;

      //Check if is B decay product
      if ( grMothD ) {
        
        if( ( grMothPDG == 5 ) || ( grMothPDG > 500 && grMothPDG < 600) || ( grMothPDG > 5000 && grMothPDG < 6000 ) ) {

          fromDfromB[jp] = 1;
          fromD[jp]      = 0;
        
        }
        

      }
      
    }
 */
    for (Int_t jv = 0; jv < mcVert; jv++ ) {
      
      if ( TMath::Abs( x[jv] - mcPart->Xv() ) + TMath::Abs( y[jv] - mcPart->Yv() ) + TMath::Abs(z[jv]-mcPart->Zv() ) < 0.0001 ) {
        ntrks[jv]++;
        vtxfound = kTRUE;
        break;
      }
      
    }
    
/*    for (Int_t jv = 0; jv < jp; jv++ ) {
      if ( (fromB[jv] == fromB[jp] && fromB[jp] == 1 && labmoth[jv]==labmoth[jp] ) || (fromDfromB[jv] == 1 && fromB[jp] == 1 && labgrmoth[jv]==labmoth[jp] ) ) {
        fromSameB[mcBDVert]++;
        vtxBDfound = kTRUE;
        break;
      }

      if ((fromB[jv]==1 && fromDfromB[jp]==1 && labmoth[jv]==labgrmoth[jp]) || (fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jp]==1 && labgrmoth[jv]==labgrmoth[jp]) || (fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jp]==1 && labgrmoth[jv]==labgrgrmothB[jp]) || (fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jp]==1 && labgrgrmothB[jv]==labgrmoth[jp]) || (fromB[jv]==1 && fromDfromB[jp]==1 && labmoth[jv]==labgrgrmothB[jp]) || (fromDfromB[jv]==1  && fromB[jp]==1 && labgrgrmothB[jv]==labmoth[jp])   ||  (fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jp]==1 && labgrmoth[jv]==labgrgrmothB2[jp]) || (fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jp]==1 && labgrgrmothB2[jv]==labgrmoth[jp]) || (fromB[jv]==1 && fromDfromB[jp]==1 && labmoth[jv]==labgrgrmothB2[jp]) || (fromDfromB[jv]==1  && fromB[jp]==1 && labgrgrmothB2[jv]==labmoth[jp])){

        fromDfromSameB[mcBDVert]++;
        vtxBDfound=kTRUE;

        break;
      }
    }

    if ( !vtxBDfound ) {
      if ( fromB[jp] == 1 || fromDfromB[jp] == 1 ) {
        mcBDVert++;
        if(fromB[jp]==1)fromSameB[mcBDVert]++;
        if(fromDfromB[jp]==1)fromDfromSameB[mcBDVert]++;
      }
    }

    for (Int_t jv = 0; jv < jp; jv++){
      if(fromB[jv]==fromB[jp] && fromB[jp]==0 && fromDfromB[jv]==fromDfromB[jp] && fromDfromB[jv]==0 && fromD[jv]==fromD[jp] && fromD[jp]==1){
        if(labgrmothD[jv]==labmothD[jp] || labmothD[jv]==labgrmothD[jp] || labmothD[jv]==labmothD[jp] || labgrgrmothD[jv]==labmothD[jp] || labgrgrmothD[jp]==labmothD[jv] ||  labgrgrmothD[jv]==labgrmothD[jp] || labgrgrmothD[jp]==labgrmothD[jv]){
          fromPromptD[mcPromptDVert]++;
          vtxPromptDfound=kTRUE;
          break;
        }
      }
    }

*/
    if (!vtxfound){
      x[mcVert] = mcPart->Xv();
      y[mcVert] = mcPart->Yv();
      z[mcVert] = mcPart->Zv();

      mcVert++;
    }
    
    if ( !vtxPromptDfound ) {
      if (fromD[jp] == 1 && fromB[jp] == 0 && fromDfromB[jp] ==0 ) {

        mcPromptDVert++;
        fromPromptD[mcPromptDVert]++;

      }
    }
  }
  
  for (Int_t jv = 0; jv < mcBDVert+1; jv++ ) {
    AliDebug(AliLog::kDebug,Form("vertexBD %d, ntracksB=%d, ntracksD=%d", jv, fromSameB[jv], fromDfromSameB[jv]));
  }
  
  for ( Int_t jv = 0; jv < mcPromptDVert + 1; jv++) {
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
