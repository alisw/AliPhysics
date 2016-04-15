/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliRDHFJetsCuts.cxx 60376 2013-01-18 14:11:18Z fprino $ */

/////////////////////////////////////////////////////////////
//
// Base class for cuts on AOD reconstructed heavy-flavour Jets
//
// Authors: Andrea Rossi (andrea.rossi@cern.ch), 
// Elena Bruna (elena.bruna@to.infn.it),  
// Sarah LaPointe (s.lapointe@cern.ch)
/////////////////////////////////////////////////////////////
#include <Riostream.h>

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliAODRecoDecayHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODRecoDecay.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliVertexerTracks.h"
#include "AliRDHFJetsCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "TRandom.h"
#include "AliRDHFJetsCutsVertex.h"


using std::cout;
using std::endl;

ClassImp(AliRDHFJetsCutsVertex)

//--------------------------------------------------------------------------
AliRDHFJetsCutsVertex::AliRDHFJetsCutsVertex(const Char_t* name, const Char_t* title) : 
  AliRDHFJetsCuts(name,title),
  fNprongs(0),
  fSecVtxWithKF(0),
  fMinPtHardestTrack(0),
  fImpPar(0),
  fDistPrimSec(0),
  fCosp(0),
  fInvMassCut(0),
  fSigvert(1000),
  fChi2(1000),
  fIsElec(0)
{
  //
  // Default Constructor
  //
  fCosp=-1.; 
}
//--------------------------------------------------------------------------
AliRDHFJetsCutsVertex::AliRDHFJetsCutsVertex(const AliRDHFJetsCutsVertex &source) :
  AliRDHFJetsCuts(source),
  fNprongs(source.fNprongs),
  fSecVtxWithKF(source.fSecVtxWithKF),
  fMinPtHardestTrack(source.fMinPtHardestTrack),
  fImpPar(source.fImpPar),
  fDistPrimSec(source.fDistPrimSec),
  fCosp(source.fCosp),
  fInvMassCut(source.fInvMassCut),
  fSigvert(source.fSigvert),
  fChi2(source.fChi2),
  fIsElec(source.fIsElec)
{
  //
  // Copy constructor
 
 
}
//--------------------------------------------------------------------------
AliRDHFJetsCutsVertex &AliRDHFJetsCutsVertex::operator=(const AliRDHFJetsCutsVertex &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFJetsCuts::operator=(source);

  fNprongs=source.fNprongs;
  fMinPtHardestTrack=source.fMinPtHardestTrack;
  fSecVtxWithKF=source.fSecVtxWithKF;
  fImpPar=source.fImpPar;
  fDistPrimSec=source.fDistPrimSec;
  fCosp=source.fCosp;
  fInvMassCut=source.fInvMassCut;
  fSigvert=source.fSigvert;
  fChi2=source.fChi2;
  fIsElec=source.fIsElec;

  return *this;
}
//--------------------------------------------------------------------------
AliRDHFJetsCutsVertex::~AliRDHFJetsCutsVertex() {
  //  
  // Default Destructor
  //


}
//--------------------------------------------------------------------------
Int_t AliRDHFJetsCutsVertex::IsVertexSelected(AliAODVertex* vert, AliAODEvent* aod, Double_t magzkG ,Double_t dispersion,Double_t massParticle){

  
  Double_t chi2 = vert->GetChi2perNDF();

  //Int_t ncontr=vert->GetNContributors();
  AliAODVertex *vtx1 = (AliAODVertex *)aod->GetPrimaryVertex();
  Double_t declen = vert->DistanceToVertex(vtx1);

  if ( declen < fDistPrimSec ) return 0; //minimum distance primary-secondary required
  if ( dispersion > fSigvert ) return 0; //maximum dispersion of tracks around sec vertex required
  if ( chi2 > fChi2 ) return 0; //maximum value of vertex chi2 required

  Int_t d0cut=0;
  Int_t sumCh = 0;

  Double_t maxpt = 0;
  Double_t pxyz[3];
  Double_t pxyzSum[4] = {0., 0., 0., 0.};
  Double_t pos[3];
  vtx1->GetPosition(pos);
  
  Double_t cosp=CosPointingAngle(vert,pos);
  if ( cosp < fCosp ) return 0; //cut on cos theta point required on sec vertex
  
  for ( Int_t jp = 0; jp < vert->GetNDaughters(); jp++ ) {
    AliAODTrack *tr = (AliAODTrack *)vert->GetDaughter(jp);
    
    tr->GetPxPyPz(pxyz);
    pxyzSum[1] += pxyz[0];
    pxyzSum[2] += pxyz[1];
    pxyzSum[3] += pxyz[2];
    pxyzSum[0] += TMath::Sqrt( massParticle*massParticle + pxyz[0]*pxyz[0] + pxyz[1]*pxyz[1] + pxyz[2]*pxyz[2] );//pion mass assumed

    Int_t charge = tr->Charge();
    sumCh += charge;

    Double_t pt = tr->Pt();
    if ( pt > maxpt) maxpt = pt;
    
    Double_t d0z0[2], covd0z0[3];
    tr->PropagateToDCA(vtx1, magzkG, kVeryBig, d0z0, covd0z0);
    if ( TMath::Abs(d0z0[0]) > fImpPar ) d0cut++;//d0z0[0]=d0xy, d0z0[1]=d0z,
    
  }

  // sll 06.15 make sure the charge of the vertex does not sum -3 or 3
  //  if ( TMath::Abs(sumCh) > 2 ) return 0;
  
  Double_t invmass = TMath::Sqrt( pxyzSum[0]*pxyzSum[0] - pxyzSum[1]*pxyzSum[1] - pxyzSum[2]*pxyzSum[2] - pxyzSum[3]*pxyzSum[3]);
 
  if ( d0cut < vert->GetNDaughters() - 1 ) return 0;   //cut on d0xy required on *all tracks* in the vertex
 
  if ( invmass < fInvMassCut ) return 0; //minimum value on inv mass of the vertex required

  if ( maxpt < fMinPtHardestTrack ) return 0; //minimum pt of the hardest track required
    
  if (fIsElec) {
    Bool_t fFlagElec = kFALSE;
    IsElecInVert(aod, vert, fFlagElec);
    
    if (!fFlagElec){
      printf("--> Jet vertex not selected --- no electron");
      return 0;
    }
  }
  return 1;
}

//---------------------------------------------------------------------------
void AliRDHFJetsCutsVertex::IsElecInVert(AliAODEvent *aod, AliAODVertex *jvert, Bool_t &fFlagElec){
  //
  // Convert to ESDtrack, get emcal info of tracks in vertex
  //
  Bool_t flagElec = kFALSE;

  for(Int_t ktracks = 0; ktracks<jvert->GetNDaughters(); ktracks++){

    AliAODTrack *aodtrack=(AliAODTrack*)jvert->GetDaughter(ktracks);
    AliESDtrack *esdTrack = new AliESDtrack(aodtrack);

    Double_t fClsE = -999, p = -999, fEovp=-999,  fdEdx=-999;
    // Track extrapolation to EMCAL
    Int_t fClsId = esdTrack->GetEMCALcluster();
    p = esdTrack->P();

    if(fClsId >0) {
      AliVCluster *cluster = aod->GetCaloCluster(fClsId);
      if ( TMath::Abs( cluster->GetTrackDx() ) < 0.05 && TMath::Abs( cluster->GetTrackDz() ) < 0.05 ) {
        
        fClsE = cluster->E();
        fEovp = fClsE/p;
        fdEdx=esdTrack->GetTPCsignal();
        // for the time being select dedx
        
        if ((fEovp > 0.8 && fEovp < 1.1) && (fdEdx > 75 && fdEdx < 105)  && !flagElec) flagElec=kTRUE;
      }
    }
  }
  fFlagElec = flagElec;

}

//--------------------------------------------------------------------------
Double_t AliRDHFJetsCutsVertex::CosPointingAngle(AliAODVertex* vtx1, Double_t point[3])
{
  //
  // Cosine of pointing angle in space assuming it is produced at "point"
  //
  Double_t px=0.,py=0.,pz=0.; 
  for(Int_t i=0;i<vtx1->GetNDaughters();i++) {
    AliAODTrack *tr=(AliAODTrack*)vtx1->GetDaughter(i);
    px+=tr->Px(); 
    py+=tr->Py(); 
    pz+=tr->Pz(); 
    
  }

  TVector3 mom(px,py,pz);
  TVector3 fline(vtx1->GetX()-point[0],
		 vtx1->GetY()-point[1],
		 vtx1->GetZ()-point[2]);

  Double_t ptot2 = mom.Mag2()*fline.Mag2();
  if(ptot2 <= 0) {
    return 0.0;
  } else {
    Double_t cos = mom.Dot(fline)/TMath::Sqrt(ptot2);
    if(cos >  1.0) cos =  1.0;
    if(cos < -1.0) cos = -1.0;
    return cos;
  }
}
