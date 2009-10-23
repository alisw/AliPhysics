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
/**************************************************************************
 *                                                                        *
 * Secondary vertexing construction Class                                 *
 *  Construct secondary vertex from Beauty hadron with electron and       *
 *  hadrons, then apply selection criteria                                *
 *                                                                        *
 * Authors:                                                               *
 *   MinJung Kweon <minjung@physi.uni-heidelberg.de>                      *
 *                                                                        *
 **************************************************************************/

#include <iostream>
#include <TH2F.h>
#include <TH3F.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TParticle.h>

#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include "AliHFEsecVtx.h"
#include <AliKFParticle.h>
#include <AliKFVertex.h>
#include <AliLog.h>
#include <AliPID.h>
#include <AliStack.h>

ClassImp(AliHFEsecVtx)

//_______________________________________________________________________________________________
AliHFEsecVtx::AliHFEsecVtx():
         fESD1(0x0)
        ,fStack(0x0)
        ,fNparents(0)
        ,fHistTagged()
        ,fPairTagged(0)
        ,fdistTwoSecVtx(-1)
        ,fcosPhi(-1)
        ,fsignedLxy(-1)
        ,finvmass(-1)
        ,finvmassSigma(-1)
        ,fBTagged(0)
        ,fBElec(0)
{ 
        // Default constructor

        Init();

}

//_______________________________________________________________________________________________
AliHFEsecVtx::AliHFEsecVtx(const AliHFEsecVtx &p):
         TObject(p)
        ,fESD1(0x0)
        ,fStack(0x0)
        ,fNparents(p.fNparents)
        ,fHistTagged()
        ,fPairTagged(p.fPairTagged)
        ,fdistTwoSecVtx(p.fdistTwoSecVtx)
        ,fcosPhi(p.fcosPhi)
        ,fsignedLxy(p.fsignedLxy)
        ,finvmass(p.finvmass)
        ,finvmassSigma(p.finvmassSigma)
        ,fBTagged(p.fBTagged)
        ,fBElec(p.fBElec)
{ 
        // Copy constructor
}

//_______________________________________________________________________________________________
AliHFEsecVtx&
AliHFEsecVtx::operator=(const AliHFEsecVtx &)
{
  // Assignment operator

  AliInfo("Not yet implemented.");
  return *this;
}

//_______________________________________________________________________________________________
AliHFEsecVtx::~AliHFEsecVtx()
{
        // Destructor

        //cout << "Analysis Done." << endl;
}

//__________________________________________
void AliHFEsecVtx::Init()
{

        // set pdg code and index

        fNparents = 7;

        fParentSelect[0][0] =  411;
        fParentSelect[0][1] =  421;
        fParentSelect[0][2] =  431;
        fParentSelect[0][3] = 4122;
        fParentSelect[0][4] = 4132;
        fParentSelect[0][5] = 4232;
        fParentSelect[0][6] = 4332;

        fParentSelect[1][0] =  511;
        fParentSelect[1][1] =  521;
        fParentSelect[1][2] =  531;
        fParentSelect[1][3] = 5122;
        fParentSelect[1][4] = 5132;
        fParentSelect[1][5] = 5232;
        fParentSelect[1][6] = 5332;

/*
        fid[0][0] = 0;
        fid[0][1] = 1;
        fid[0][2] = 2;

        fid[1][0] = 0;
        fid[1][1] = 1;
        fid[1][2] = 3;

        fid[2][0] = 0;
        fid[2][1] = 2;
        fid[2][2] = 3;

        fid[3][0] = 1;
        fid[3][1] = 2;
        fid[3][2] = 3;

        fia[0][0][0] = 0;
        fia[0][0][1] = 1;
        fia[0][1][0] = 0;
        fia[0][1][1] = 2;
        fia[0][2][0] = 1;
        fia[0][2][1] = 2;

        fia[1][0][0] = 0;
        fia[1][0][1] = 1;
        fia[1][1][0] = 0;
        fia[1][1][1] = 3;
        fia[1][2][0] = 1;
        fia[1][2][1] = 3;

        fia[2][0][0] = 0;
        fia[2][0][1] = 2;
        fia[2][1][0] = 0;
        fia[2][1][1] = 3;
        fia[2][2][0] = 2;
        fia[2][2][1] = 3;

        fia[3][0][0] = 1;
        fia[3][0][1] = 2;
        fia[3][1][0] = 1;
        fia[3][1][1] = 3;
        fia[3][2][0] = 2;
        fia[3][2][1] = 3;
*/

} 

//__________________________________________
void AliHFEsecVtx::ResetTagVar()
{
        // reset tag variables

        fdistTwoSecVtx = -1;
        fcosPhi = -1;
        fsignedLxy = -1;
        finvmass = -1;
        finvmassSigma = -1;
        fBTagged = kFALSE;
        fBElec = kFALSE;
}

//__________________________________________
void AliHFEsecVtx::InitAnaPair()
{
        // initialize pair tagging variables

        fPairTagged = 0;
        for (Int_t i=0; i<20; i++){
                fpairedTrackID[i] = -1;
                fpairedChi2[i] = -1;
                fpairedInvMass[i] = -1;
                fpairedSignedLxy[i] = -1;
        }

        fid[0][0] = 0;
        fid[0][1] = 1;
        fid[0][2] = 2;

        fid[1][0] = 0;
        fid[1][1] = 1;
        fid[1][2] = 3;

        fid[2][0] = 0;
        fid[2][1] = 2;
        fid[2][2] = 3;

        fid[3][0] = 1;
        fid[3][1] = 2;
        fid[3][2] = 3;

        fia[0][0][0] = 0;
        fia[0][0][1] = 1;
        fia[0][1][0] = 0;
        fia[0][1][1] = 2;
        fia[0][2][0] = 1;
        fia[0][2][1] = 2;

        fia[1][0][0] = 0;
        fia[1][0][1] = 1;
        fia[1][1][0] = 0;
        fia[1][1][1] = 3;
        fia[1][2][0] = 1;
        fia[1][2][1] = 3;

        fia[2][0][0] = 0;
        fia[2][0][1] = 2;
        fia[2][1][0] = 0;
        fia[2][1][1] = 3;
        fia[2][2][0] = 2;
        fia[2][2][1] = 3;

        fia[3][0][0] = 1;
        fia[3][0][1] = 2;
        fia[3][1][0] = 1;
        fia[3][1][1] = 3;
        fia[3][2][0] = 2;
        fia[3][2][1] = 3;

}

//_______________________________________________________________________________________________
void AliHFEsecVtx::CreateHistograms(TString hnopt)
{ 
        // create histograms

        fkSourceLabel[kAll]="all";
        fkSourceLabel[kDirectCharm]="directCharm";
        fkSourceLabel[kDirectBeauty]="directBeauty";
        fkSourceLabel[kBeautyCharm]="beauty2charm";
        fkSourceLabel[kGamma]="gamma";
        fkSourceLabel[kPi0]="pi0";
        fkSourceLabel[kElse]="others";
        fkSourceLabel[kBeautyGamma]="beauty22gamma";
        fkSourceLabel[kBeautyPi0]="beauty22pi0";
        fkSourceLabel[kBeautyElse]="beauty22others";


        TString hname;
        for (Int_t isource = 0; isource < 10; isource++ ){

           hname=hnopt+"InvMass_"+fkSourceLabel[isource];
           fHistPair[isource].fInvMass = new TH2F(hname,hname+";invMass;invMassSigma",120,-2,10,100,0,10);
           hname=hnopt+"InvMassCut1_"+fkSourceLabel[isource];
           fHistPair[isource].fInvMassCut1 = new TH2F(hname,hname+";invMass;invMassSigma",120,-2,10,100,0,10);
           hname=hnopt+"InvMassCut2_"+fkSourceLabel[isource];
           fHistPair[isource].fInvMassCut2 = new TH2F(hname,hname+";invMass;invMassSigma",120,-2,10,100,0,10);
           hname=hnopt+"KFChi2_"+fkSourceLabel[isource];
           fHistPair[isource].fKFChi2 = new TH1F(hname,hname,200,0,20);
           hname=hnopt+"OpenAngle_"+fkSourceLabel[isource];
           fHistPair[isource].fOpenAngle = new TH1F(hname,hname,100,0,3.14);
           hname=hnopt+"CosOpenAngle_"+fkSourceLabel[isource];
           fHistPair[isource].fCosOpenAngle = new TH1F(hname,hname,100,-1.1,1.1);
           hname=hnopt+"SignedLxy_"+fkSourceLabel[isource];
           fHistPair[isource].fSignedLxy = new TH2F(hname,hname,1000,-5,5,120,-2,10);
           hname=hnopt+"KFIP_"+fkSourceLabel[isource];
           fHistPair[isource].fKFIP = new TH1F(hname,hname,1000,-5,5);
           hname=hnopt+"IPMax_"+fkSourceLabel[isource];
           fHistPair[isource].fIPMax= new TH1F(hname,hname,500,0,5);

        }

        hname=hnopt+"pt_beautyelec";
        fHistTagged.fPtBeautyElec= new TH1F(hname,hname,250,0,50);
        hname=hnopt+"pt_taggedelec";
        fHistTagged.fPtTaggedElec= new TH1F(hname,hname,250,0,50);
        hname=hnopt+"pt_righttaggedelec";
        fHistTagged.fPtRightTaggedElec = new TH1F(hname,hname,250,0,50);
        hname=hnopt+"pt_wrongtaggedelec";
        fHistTagged.fPtWrongTaggedElec = new TH1F(hname,hname,250,0,50);

        hname=hnopt+"InvmassBeautyElecSecVtx";
        fHistTagged.fInvmassBeautyElecSecVtx= new TH1F(hname,hname,120,-2,10);
        hname=hnopt+"InvmassNotBeautyElecSecVtx";
        fHistTagged.fInvmassNotBeautyElecSecVtx= new TH1F(hname,hname,120,-2,10);

        hname=hnopt+"SignedLxyBeautyElecSecVtx";
        fHistTagged.fSignedLxyBeautyElecSecVtx= new TH1F(hname,hname,1000,-5,5);
        hname=hnopt+"SignedLxyNotBeautyElecSecVtx";
        fHistTagged.fSignedLxyNotBeautyElecSecVtx= new TH1F(hname,hname,1000,-5,5);

        hname=hnopt+"DistTwoVtxBeautyElecSecVtx";
        fHistTagged.fDistTwoVtxBeautyElecSecVtx= new TH1F(hname,hname,1000,-5,5);
        hname=hnopt+"DistTwoVtxNotBeautyElecSecVtx";
        fHistTagged.fDistTwoVtxNotBeautyElecSecVtx= new TH1F(hname,hname,1000,-5,5);

        hname=hnopt+"CosPhiBeautyElecSecVtx";
        fHistTagged.fCosPhiBeautyElecSecVtx= new TH1F(hname,hname,100,-1.1,1.1);
        hname=hnopt+"CosPhiNotBeautyElecSecVtx";
        fHistTagged.fCosPhiNotBeautyElecSecVtx= new TH1F(hname,hname,100,-1.1,1.1);

        hname=hnopt+"Chi2BeautyElecSecVtx";
        fHistTagged.fChi2BeautyElecSecVtx= new TH1F(hname,hname,200,0,20);
        hname=hnopt+"Chi2NotBeautyElecSecVtx";
        fHistTagged.fChi2NotBeautyElecSecVtx= new TH1F(hname,hname,200,0,20);

        hname=hnopt+"InvmassBeautyElec2trkVtx";
        fHistTagged.fInvmassBeautyElec2trkVtx= new TH1F(hname,hname,120,-2,10);
        hname=hnopt+"InvmassNotBeautyElec2trkVtx";
        fHistTagged.fInvmassNotBeautyElec2trkVtx= new TH1F(hname,hname,120,-2,10);

        hname=hnopt+"SignedLxyBeautyElec2trkVtx";
        fHistTagged.fSignedLxyBeautyElec2trkVtx= new TH1F(hname,hname,1000,-5,5);
        hname=hnopt+"SignedLxyNotBeautyElec2trkVtx";
        fHistTagged.fSignedLxyNotBeautyElec2trkVtx= new TH1F(hname,hname,1000,-5,5);

        hname=hnopt+"Chi2BeautyElec2trkVtx";
        fHistTagged.fChi2BeautyElec2trkVtx= new TH1F(hname,hname,200,0,20);
        hname=hnopt+"Chi2NotBeautyElec2trkVtx";
        fHistTagged.fChi2NotBeautyElec2trkVtx= new TH1F(hname,hname,200,0,20);
}

//_______________________________________________________________________________________________
void AliHFEsecVtx::AnaPair(AliESDtrack* track1, AliESDtrack* track2, Int_t index2)
{
        // calculate e-h pair characteristics and tag pair 

        Int_t sourcePart = PairCode(track1,track2);

        // get KF particle input pid
        Int_t pdg1 = GetMCPID(track1);
        Int_t pdg2 = GetMCPID(track2);
        

        // create KF particle of pair
        AliKFParticle::SetField(fESD1->GetMagneticField());
        AliKFParticle kfTrack1(*track1, pdg1);
        AliKFParticle kfTrack2(*track2, pdg2);

        AliKFParticle kfSecondary(kfTrack1,kfTrack2);

        // copy primary vertex from ESD
        AliKFVertex primVtxCopy(*(fESD1->GetPrimaryVertex()));
        if( primVtxCopy.GetNDF() <1 ) return;

        //primary vertex point
        Double_t pvx = primVtxCopy.GetX();
        Double_t pvy = primVtxCopy.GetY();
        //Double_t pvz = primVtxCopy.GetZ();

        //secondary vertex point from kf particle
        Double_t kfx = kfSecondary.GetX();
        Double_t kfy = kfSecondary.GetY();
        //Double_t kfz = kfSecondary.GetZ();
        
        //momentum at the decay point from kf particle
        Double_t kfpx = kfSecondary.GetPx();
        Double_t kfpy = kfSecondary.GetPy();
        //Double_t kfpz = kfSecondary.GetPz();
        

        Double_t dx = kfx-pvx;
        Double_t dy = kfy-pvy;



        // discriminating variables ----------------------------------------------------------

        // invariant mass of the KF particle
        Double_t invmass = -1;
        Double_t invmassSigma = -1;
        kfSecondary.GetMass(invmass,invmassSigma);

        // chi2 of the KF particle
        Double_t kfchi2 = TMath::Sqrt(TMath::Abs(kfSecondary.GetChi2()/kfSecondary.GetNDF()));

        // opening angle between two particles in XY plane
        Double_t phi = kfTrack1.GetAngleXY(kfTrack2);
        Double_t cosphi = TMath::Cos(phi);

        // projection of kf vertex vector to the kf momentum direction 
        Double_t costheta = ( dx*kfpx + dy*kfpy)/TMath::Sqrt(dx*dx+dy*dy)*TMath::Sqrt(kfpx*kfpx + kfpy*kfpy);
        Double_t signedLxy = TMath::Sqrt(dx*dx+dy*dy)*costheta;

        // DCA from primary to e-h KF particle (impact parameter of KF particle)
        Double_t vtx[2]={pvx, pvy}; 
        Double_t kfip = kfSecondary.GetDistanceFromVertexXY(vtx);


	      Float_t dcaR=-1; 
        Float_t dcaR1=-1, dcaR2=-1;
        Float_t dcaZ1=-1, dcaZ2=-1;
        track1->GetImpactParameters(dcaR1,dcaZ1);
        track2->GetImpactParameters(dcaR2,dcaZ2);

	      if (TMath::Abs(dcaR1) >= TMath::Abs(dcaR2)) dcaR=dcaR1;
	      else dcaR=dcaR2;

        // fill histograms 
        fHistPair[sourcePart].fInvMass->Fill(invmass,invmassSigma);
        fHistPair[sourcePart].fKFChi2->Fill(kfchi2);
        fHistPair[sourcePart].fOpenAngle->Fill(phi);
        fHistPair[sourcePart].fCosOpenAngle->Fill(cosphi);
        fHistPair[sourcePart].fSignedLxy->Fill(signedLxy,invmass);
        fHistPair[sourcePart].fKFIP->Fill(kfip);
        fHistPair[sourcePart].fIPMax->Fill(TMath::Abs(dcaR));

        // pair cuts 
        if( kfchi2 >2. ) return;

        if ( signedLxy > 0.05 && cosphi > 0.5) fHistPair[sourcePart].fInvMassCut1->Fill(invmass,invmassSigma);

        // pair tagging if it passed the above cuts
        if(signedLxy > 0. && cosphi > 0.) fHistPair[sourcePart].fInvMassCut2->Fill(invmass,invmassSigma);


        // pair tagging condition
        if ( signedLxy > 0.0 && cosphi > 0) { // testing loose cut
        //if ( signedLxy > 0.06 && cosphi > 0) {
          fpairedTrackID[fPairTagged] = index2;
          fpairedChi2[fPairTagged] = kfchi2;
          fpairedInvMass[fPairTagged] = invmass;
          fpairedSignedLxy[fPairTagged] = signedLxy;
          fPairTagged++;
        }

}

//_______________________________________________________________________________________________
void AliHFEsecVtx::RunSECVTX(AliESDtrack *track)
{
        // run secondary vertexing algorithm and do tagging

        ResetTagVar();

        Int_t imclabel = TMath::Abs(track->GetLabel());
        if(imclabel<0) return;
        TParticle* mcpart = fStack->Particle(imclabel);
        Int_t esource = GetElectronSource(imclabel);
        if (esource == kDirectBeauty || esource == kBeautyCharm || esource == kBeautyGamma || esource == kBeautyPi0 || esource == kBeautyElse){
                fHistTagged.fPtBeautyElec->Fill(mcpart->Pt());
		fBElec = kTRUE;
        }


        if (fPairTagged >= 4) {
          FindSECVTXCandid4Tracks(track);         
        }
        else if (fPairTagged == 3) {
          FindSECVTXCandid3Tracks(track);         
        }
        else if (fPairTagged == 2) {
          FindSECVTXCandid2Tracks(track);         
        }
        else if (fPairTagged == 1) {
          ApplyPairTagCut();      
        }


        if (fBTagged) {
                fHistTagged.fPtTaggedElec->Fill(mcpart->Pt());
                if (esource == kDirectBeauty || esource == kBeautyCharm || esource == kBeautyGamma || esource == kBeautyPi0 || esource == kBeautyElse){
                        fHistTagged.fPtRightTaggedElec->Fill(mcpart->Pt());
                }
                else fHistTagged.fPtWrongTaggedElec->Fill(mcpart->Pt());
        }

}

//_______________________________________________________________________________________________
void AliHFEsecVtx::ApplyPairTagCut()
{
        // apply tagging cut for e-h pair

        if (fBElec){
                fHistTagged.fInvmassBeautyElec2trkVtx->Fill(fpairedInvMass[0]);
                fHistTagged.fSignedLxyBeautyElec2trkVtx->Fill(fpairedSignedLxy[0]);
                fHistTagged.fChi2BeautyElec2trkVtx->Fill(fpairedChi2[0]);
        }
        else {
                fHistTagged.fInvmassNotBeautyElec2trkVtx->Fill(fpairedInvMass[0]);
                fHistTagged.fSignedLxyNotBeautyElec2trkVtx->Fill(fpairedSignedLxy[0]);
                fHistTagged.fChi2NotBeautyElec2trkVtx->Fill(fpairedChi2[0]);
        }
                
        if (fpairedChi2[0] > 2.0) return;
        if (fpairedInvMass[0] < 1.5) return;
        if (fpairedSignedLxy[0] < 0.05) return;

        fBTagged = kTRUE;
}

//_______________________________________________________________________________________________
Bool_t AliHFEsecVtx::ApplyPairTagCut(Int_t id) 
{
        // apply tagging cut for e-h pair of given indexed hadron

        if (fpairedChi2[id] > 2.0) return kFALSE;
        if (fpairedInvMass[id] < 1.5) return kFALSE;
        if (fpairedSignedLxy[id] < 0.05) return kFALSE;

        fBTagged = kTRUE;
        return kTRUE;
                
}

//_______________________________________________________________________________________________
Bool_t AliHFEsecVtx::ApplyTagCut(Double_t chi2) 
{
        // apply tagging cut for secondary vertex

        if (chi2 > 2.0) return kFALSE;
        if (finvmass < 1.5) return kFALSE;
        if (fsignedLxy < 0.05) return kFALSE;
        if (fcosPhi < 0.90) return kFALSE;
        if (fdistTwoSecVtx > 0.1) return kFALSE;
        
        fBTagged = kTRUE;
        return kTRUE;
}

//_______________________________________________________________________________________________
void AliHFEsecVtx::ApplyVtxTagCut(Double_t chi2, Int_t id1, Int_t id2)
{
        // apply tagging cut for e-h pair of given indexed hadron

        if(!ApplyTagCut(chi2)){
                if(!ApplyPairTagCut(id1)) ApplyPairTagCut(id2);
        }
        
}

//_______________________________________________________________________________________________
void AliHFEsecVtx::FindSECVTXCandid4Tracks(AliESDtrack *track) 
{
        // find secondary vertex for >= 4 e-h pairs

        // sort pair in increasing order (kFALSE-increasing order)
        Int_t index[20];
        Int_t indexA[4];
        Int_t indexB[3];
        Double_t sevchi2[4];
        AliESDtrack *htrack[4];

        if(fPairTagged>20) return; // protection

        TMath::Sort(fPairTagged,fpairedChi2,index,kFALSE);

        // select 4 partner tracks retruning smallest pair chi2 

        for (Int_t i=0; i<4; i++){
                htrack[i] = fESD1->GetTrack(fpairedTrackID[index[i]]);
        }
        
        // calculated chi2 with 1 electron and 3 partner tracks 
        for (Int_t i=0; i<4; i++){
                sevchi2[i] = GetSecVtxChi2(track, htrack[fid[i][0]], htrack[fid[i][1]], htrack[fid[i][2]]);
        }

        // select 3 partner tracks retruning smallest pair chi2
        // [think] if two smallest chi2 are similar, have to think about better handling of selection
        TMath::Sort(4,sevchi2,indexA,kFALSE);

        // calculated chi2 with 1 electron and 2 partner tracks 
        for (Int_t i=0; i<3; i++){
                sevchi2[i] = GetSecVtxChi2(track, htrack[fia[indexA[0]][i][0]], htrack[fia[indexA[0]][i][1]]);
        }

        // select 2 partner tracks retruning smallest pair chi2
        TMath::Sort(3,sevchi2,indexB,kFALSE);

        // calculate secondary vertex quality variables with 1 electron and 2 hadrons
        CalcSECVTXProperty(track,htrack[fia[indexA[0]][indexB[0]][0]],htrack[fia[indexA[0]][indexB[0]][1]]);

        if (fBElec){
                fHistTagged.fChi2BeautyElecSecVtx->Fill(sevchi2[indexB[0]]);
        }
        else {
                fHistTagged.fChi2NotBeautyElecSecVtx->Fill(sevchi2[indexB[0]]);
        }

        ApplyVtxTagCut(sevchi2[indexB[0]],index[fia[indexA[0]][indexB[0]][0]],index[fia[indexA[0]][indexB[0]][1]]);

}

//_______________________________________________________________________________________________
void AliHFEsecVtx::FindSECVTXCandid3Tracks(AliESDtrack *track) 
{
        // find secondary vertex for 3 e-h pairs

        // sort pair in increasing order (kFALSE-increasing order)
        Int_t indexA[1] = { 0 };
        Int_t indexB[3];
        Double_t sevchi2[3];
        AliESDtrack *htrack[3];

        // select 4 partner tracks retruning smallest pair chi2 

        for (Int_t i=0; i<3; i++){
                htrack[i] = fESD1->GetTrack(fpairedTrackID[i]);
        }
        
        // calculated chi2 with 1 electron and 2 partner tracks 
        for (Int_t i=0; i<3; i++){
                sevchi2[i] = GetSecVtxChi2(track, htrack[fia[indexA[0]][i][0]], htrack[fia[indexA[0]][i][1]]);
        }

        // select 2 partner tracks retruning smallest pair chi2
        TMath::Sort(3,sevchi2,indexB,kFALSE);

        // calculate secondary vertex quality variables with 1 electron and 2 hadrons
        CalcSECVTXProperty(track,htrack[fia[indexA[0]][indexB[0]][0]],htrack[fia[indexA[0]][indexB[0]][1]]);

        if (fBElec){
                fHistTagged.fChi2BeautyElecSecVtx->Fill(sevchi2[indexB[0]]);
        }
        else {
                fHistTagged.fChi2NotBeautyElecSecVtx->Fill(sevchi2[indexB[0]]);
        }

        ApplyVtxTagCut(sevchi2[indexB[0]],fia[indexA[0]][indexB[0]][0],fia[indexA[0]][indexB[0]][1]);
}

//_______________________________________________________________________________________________
void AliHFEsecVtx::FindSECVTXCandid2Tracks(AliESDtrack *track) 
{
        // find secondary vertex for 2 e-h pairs

        Double_t sevchi2[1];
        AliESDtrack *htrack[2];

        for (Int_t i=0; i<2; i++){
                htrack[i] = fESD1->GetTrack(fpairedTrackID[i]);
        }
        
        sevchi2[0] = GetSecVtxChi2(track, htrack[0], htrack[1]);

        // calculate secondary vertex quality variables with 1 electron and 2 hadrons
        CalcSECVTXProperty(track,htrack[0],htrack[1]);

        if (fBElec){
                fHistTagged.fChi2BeautyElecSecVtx->Fill(sevchi2[0]);
        }
        else {
                fHistTagged.fChi2NotBeautyElecSecVtx->Fill(sevchi2[0]);
        }

        ApplyVtxTagCut(sevchi2[0],0,1);
}

//_______________________________________________________________________________________________
void AliHFEsecVtx::CalcSECVTXProperty(AliESDtrack* track1, AliESDtrack* track2, AliESDtrack* track3)
{
        // calculate secondary vertex properties

        Int_t pdg1 = GetMCPID(track1);
        Int_t pdg2 = GetMCPID(track2);
        Int_t pdg3 = GetMCPID(track3);

        AliKFParticle::SetField(fESD1->GetMagneticField());
        AliKFParticle kfTrack1(*track1, pdg1);
        AliKFParticle kfTrack2(*track2, pdg2);
        AliKFParticle kfTrack3(*track3, pdg3);

        AliKFParticle kfSecondary12(kfTrack1,kfTrack2);
        AliKFParticle kfSecondary13(kfTrack1,kfTrack3);
        AliKFParticle kfSecondary(kfTrack1,kfTrack2,kfTrack3);

        // copy primary vertex from ESD
        AliKFVertex primVtxCopy(*(fESD1->GetPrimaryVertex()));        
        //printf("primary ndf= %d\n",primVtxCopy.GetNDF());
        if( primVtxCopy.GetNDF() <1 ) return;

        Double_t kdx12 = kfSecondary12.GetX()-primVtxCopy.GetX();
        Double_t kdy12 = kfSecondary12.GetY()-primVtxCopy.GetY();
        //Double_t kdz12 = kfSecondary12.GetZ()-primVtxCopy.GetZ();

        Double_t kdx13 = kfSecondary13.GetX()-primVtxCopy.GetX();
        Double_t kdy13 = kfSecondary13.GetY()-primVtxCopy.GetY();
        //Double_t kdz13 = kfSecondary13.GetZ()-primVtxCopy.GetZ();

        Double_t kdx = kfSecondary.GetX()-primVtxCopy.GetX();
        Double_t kdy = kfSecondary.GetY()-primVtxCopy.GetY();
        //Double_t kdz = kfSecondary.GetZ()-primVtxCopy.GetZ();

        // calculate distance and angle between two secvtxes 
        fdistTwoSecVtx = TMath::Sqrt((kdx12-kdx13)*(kdx12-kdx13) + (kdy12-kdy13)*(kdy12-kdy13));
        fcosPhi = ( kdx12*kdx13 + kdy12*kdy13 ) / ( TMath::Sqrt(kdx12*kdx12+kdy12*kdy12)*TMath::Sqrt(kdx13*kdx13+kdy13*kdy13) );
        //Double_t lengthdiff = TMath::Abs(TMath::Sqrt(kdx12*kdx12+kdy12*kdy12) - TMath::Sqrt(kdx13*kdx13+kdy13*kdy13));

        // calculate angle between secondary vertex vector and secondary particle momentum vector in transverse plane
        Double_t cosTheta = ( kdx*kfSecondary.GetPx() + kdy*kfSecondary.GetPy()) / TMath::Sqrt(kdx*kdx+kdy*kdy)*TMath::Sqrt(kfSecondary.GetPx()*kfSecondary.GetPx()+kfSecondary.GetPy()*kfSecondary.GetPy());
        // calculate signed Lxy
        fsignedLxy = TMath::Sqrt(kdx*kdx+kdy*kdy)*cosTheta;

        // calculate invariant mass of the kf particle
        kfSecondary.GetMass(finvmass,finvmassSigma);

	if (fBElec){
		fHistTagged.fInvmassBeautyElecSecVtx->Fill(finvmass);
		fHistTagged.fSignedLxyBeautyElecSecVtx->Fill(fsignedLxy);
		fHistTagged.fDistTwoVtxBeautyElecSecVtx->Fill(fdistTwoSecVtx);
		fHistTagged.fCosPhiBeautyElecSecVtx->Fill(fcosPhi);
	}
	else {
		fHistTagged.fInvmassNotBeautyElecSecVtx->Fill(finvmass);
		fHistTagged.fSignedLxyNotBeautyElecSecVtx->Fill(fsignedLxy);
		fHistTagged.fDistTwoVtxNotBeautyElecSecVtx->Fill(fdistTwoSecVtx);
		fHistTagged.fCosPhiNotBeautyElecSecVtx->Fill(fcosPhi);
	}

}

//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::GetMCPID(AliESDtrack *track) 
{
        // return mc pid

        Int_t label = TMath::Abs(track->GetLabel());
        TParticle* mcpart = fStack->Particle(label);
        if ( !mcpart ) return 0;
        Int_t pdgCode = mcpart->GetPdgCode();

        return pdgCode;
}

//_______________________________________________________________________________________________
Double_t AliHFEsecVtx::GetSecVtxChi2(AliESDtrack* track1, AliESDtrack* track2, AliESDtrack* track3, AliESDtrack* track4)
{
        // return 4 track secondary vertex chi2

        Int_t pdg1 = GetMCPID(track1);
        Int_t pdg2 = GetMCPID(track2);
        Int_t pdg3 = GetMCPID(track3);
        Int_t pdg4 = GetMCPID(track4);

        AliKFParticle::SetField(fESD1->GetMagneticField());
        AliKFParticle kfTrack1(*track1, pdg1);
        AliKFParticle kfTrack2(*track2, pdg2);
        AliKFParticle kfTrack3(*track3, pdg3);
        AliKFParticle kfTrack4(*track4, pdg4);
        AliKFParticle kfSecondary(kfTrack1,kfTrack2,kfTrack3,kfTrack4);

        return TMath::Sqrt(TMath::Abs(kfSecondary.GetChi2()/kfSecondary.GetNDF()));

}

//_______________________________________________________________________________________________
Double_t AliHFEsecVtx::GetSecVtxChi2(AliESDtrack* track1, AliESDtrack* track2, AliESDtrack* track3)
{
        // return 3 track secondary vertex chi2

        Int_t pdg1 = GetMCPID(track1);
        Int_t pdg2 = GetMCPID(track2);
        Int_t pdg3 = GetMCPID(track3);

        AliKFParticle::SetField(fESD1->GetMagneticField());
        AliKFParticle kfTrack1(*track1, pdg1);
        AliKFParticle kfTrack2(*track2, pdg2);
        AliKFParticle kfTrack3(*track3, pdg3);
        AliKFParticle kfSecondary(kfTrack1,kfTrack2,kfTrack3);

        return TMath::Sqrt(TMath::Abs(kfSecondary.GetChi2()/kfSecondary.GetNDF()));

}

//_______________________________________________________________________________________________
Double_t AliHFEsecVtx::GetSecVtxChi2(AliESDtrack* track1, AliESDtrack* track2)
{
        // return 2 track secondary vertex chi2

        Int_t pdg1 = GetMCPID(track1);
        Int_t pdg2 = GetMCPID(track2);

        AliKFParticle::SetField(fESD1->GetMagneticField());
        AliKFParticle kfTrack1(*track1, pdg1);
        AliKFParticle kfTrack2(*track2, pdg2);
        AliKFParticle kfSecondary(kfTrack1,kfTrack2);

        return TMath::Sqrt(TMath::Abs(kfSecondary.GetChi2()/kfSecondary.GetNDF()));

}

//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::PairOrigin(AliESDtrack* track1, AliESDtrack* track2)
{

        //
        // return pdg code of the origin(source) of the pair 
        // 
        //
        // ---*---*---*-----ancester A----- track1
        //                        |____*______ 
        //                             |______ track2
        // => if they originated from same ancester, 
        //    then return "the absolute value of pdg code of ancester A"
        //
        // ---*---*---B hadron-----ancester A----- track1
        //                               |____*______ 
        //                                    |______ track2
        // => if they originated from same ancester, and this ancester originally comes from B hadrons
        //    then return -1*"the absolute value of pdg code of ancester A"
        //
        // caution : it can also return parton pdg code if it originated from same string or gluon spliting. 
        //           

        if (track1->GetLabel()<0 || track2->GetLabel()<0) return 0;
        TParticle* part1 = fStack->Particle(TMath::Abs(track1->GetLabel()));
        TParticle* part2 = fStack->Particle(TMath::Abs(track2->GetLabel()));
        if (!(part1)) return 0;
        if (!(part2)) return 0;

        TParticle* part1Crtgen = part1; // copy track into current generation particle
        TParticle* part2Crtgen = part2; // copy track into current generation particle


        Int_t sourcePDG = 0;

        // if the two tracks' mother's label is same, get the mother info
        // in case of charm, check if it originated from beauty
        for (Int_t i=0; i<100; i++){ // iterate 100
                Int_t iLabel = part1Crtgen->GetFirstMother(); //label of mother of current generation for 1st partilce
                if (iLabel < 0) return 0;

                for (Int_t j=0; j<100; j++){ // iterate 100
                        Int_t jLabel = part2Crtgen->GetFirstMother(); //label of mother of current generation for 2nd partilce
                        if (jLabel < 0) return 0; // if jLabel == -1

                        if (iLabel == jLabel){ // check if two tracks are originated from same mother
                                TParticle* thismother = fStack->Particle(jLabel); // if yes, get "thismother" info 
                                sourcePDG = abs(thismother->GetPdgCode()); // get the pdg code of "this mother"

                                // check ancester to see if it is originally from beauty 
                                for (Int_t k=0; k<10; k++){ // check up to 10 ancesters
                                        Int_t ancesterLabel = thismother->GetFirstMother();
                                        if (ancesterLabel < 0) return sourcePDG; // if "thismoter" doesn't have mother anymore, return thismother's pdg  

                                        TParticle* thisancester = fStack->Particle(ancesterLabel);
                                        Int_t ancesterPDG = abs(thisancester->GetPdgCode());

                                        for (Int_t l=0; l<fNparents; l++){
                                                if (abs(ancesterPDG)==fParentSelect[1][l]){
                                                        sourcePDG = -1*sourcePDG; // multiply -1 for charm from bottom
                                                        return sourcePDG;
                                                }
                                        }
                                        thismother = thisancester;
                                }

                        }
                        part2Crtgen = fStack->Particle(jLabel); // if their mother is different, go up to previous generation of 2nd particle
                }
                part1Crtgen = fStack->Particle(iLabel); // if their mother is different, go up to previous generation of 2nd particle

                // if you don't find additionional(2nd particle) track originated from a given beauty hadron, break to save time
                Int_t motherPDGtmp = abs(part1Crtgen->GetPdgCode());
                for (Int_t l=0; l<fNparents; l++){
                        if (abs(motherPDGtmp)==fParentSelect[1][l]) return sourcePDG;
                }

        }

        return sourcePDG;

}

//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::PairCode(AliESDtrack* track1,AliESDtrack* track2)
{

        //           
        // return pair code which is predefinded as:
        //  kDirectCharm, kDirectBeauty, kBeautyCharm, kGamma, kPi0, kElse, kBeautyGamma, kBeautyPi0, kBeautyElse
        //           

        Int_t pairOriginsPDG = PairOrigin(track1,track2);

        Int_t sourcePart = kElse;

        if (pairOriginsPDG < 0) {
                sourcePart = kBeautyElse;
        }
        for (Int_t i=0; i<fNparents; i++){
                if (abs(pairOriginsPDG)==fParentSelect[0][i]) {
                        if (pairOriginsPDG>0) sourcePart = kDirectCharm;
                        if (pairOriginsPDG<0) {
                                sourcePart = kBeautyCharm;
                        }
                }
                if (abs(pairOriginsPDG)==fParentSelect[1][i]) {
                        if (pairOriginsPDG>0) {
                                sourcePart = kDirectBeauty;
                        }
                        if (pairOriginsPDG<0)  return kElse;
                }
        }
        if (pairOriginsPDG == 22) sourcePart = kGamma;
        if (pairOriginsPDG == -22) {
                sourcePart = kBeautyGamma;
        }
        if (pairOriginsPDG == 111) sourcePart = kPi0;
        if (pairOriginsPDG == -111) {
                sourcePart = kBeautyPi0;
        }

        return sourcePart;

}

//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::GetElectronSource(Int_t iTrack) 
{

    // return decay electron's origin 

    if (iTrack < 0) {
      AliDebug(1, "Stack label is negative, return\n");
      return -1;
    }

    TParticle* mcpart = fStack->Particle(iTrack);

    if ( abs(mcpart->GetPdgCode()) != 11 ) return -1; // check if it is electron !!!!!!!!!!!!!!!!!

    Int_t iLabel = mcpart->GetFirstMother();
    if (iLabel<0){
      AliDebug(1, "Stack label is negative, return\n");
      return -1;
    }

    Int_t origin = -1;
    Bool_t isFinalOpenCharm = kFALSE;

    TParticle *partMother = fStack->Particle(iLabel);
    Int_t maPdgcode = partMother->GetPdgCode();

    // if the mother is charmed hadron  
    if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) {

         for (Int_t i=0; i<fNparents; i++){
            if (abs(maPdgcode)==fParentSelect[0][i]){
              isFinalOpenCharm = kTRUE;
            }
         }
         if (!isFinalOpenCharm) return -1;


          // iterate until you find B hadron as a mother or become top ancester 
          for (Int_t i=1; i<100; i++){ // check back to the 100 generation older

             Int_t jLabel = partMother->GetFirstMother();
             if (jLabel == -1){
               origin = kDirectCharm;
               return origin;
             }
             if (jLabel < 0){ // safety protection
               AliDebug(1, "Stack label is negative, return\n");
               return -1;
             }

             // if there is an ancester
             TParticle* grandMa = fStack->Particle(jLabel);
             Int_t grandMaPDG = grandMa->GetPdgCode();

             for (Int_t j=0; j<fNparents; j++){
                if (abs(grandMaPDG)==fParentSelect[1][j]){

                  origin = kBeautyCharm;
                  return origin;
                }
             }

             partMother = grandMa;
          } // end of iteration 
    } // end of if
    else if ( int(abs(maPdgcode)/100.) == kBeauty || int(abs(maPdgcode)/1000.) == kBeauty ) {
         for (Int_t i=0; i<fNparents; i++){
            if (abs(maPdgcode)==fParentSelect[1][i]){
              origin = kDirectBeauty;
              return origin;
            }
         }
    } // end of if


    //============ gamma ================
    else if ( abs(maPdgcode) == 22 ) {
         origin = kGamma;

          // iterate until you find B hadron as a mother or become top ancester 
          for (Int_t i=1; i<100; i++){ // check back to the 100 generation older

             Int_t jLabel = partMother->GetFirstMother();
             if (jLabel == -1){
               origin = kGamma;
               return origin;
             }
             if (jLabel < 0){ // safety protection
               AliDebug(1, "Stack label is negative, return\n");
               return -1;
             }

             // if there is an ancester
             TParticle* grandMa = fStack->Particle(jLabel);
             Int_t grandMaPDG = grandMa->GetPdgCode();

             for (Int_t j=0; j<fNparents; j++){
                if (abs(grandMaPDG)==fParentSelect[1][j]){
                  origin = kBeautyGamma;
                  return origin;
                }
             }

             partMother = grandMa;
          } // end of iteration 

         return origin;
    } // end of if

    //============ pi0 ================
    else if ( abs(maPdgcode) == 111 ) {
         origin = kPi0;

          // iterate until you find B hadron as a mother or become top ancester 
          for (Int_t i=1; i<100; i++){ // check back to the 100 generation older

             Int_t jLabel = partMother->GetFirstMother();
             if (jLabel == -1){
               origin = kPi0;
               return origin;
             }
             if (jLabel < 0){ // safety protection
               AliDebug(1, "Stack label is negative, return\n");
               return -1;
             }

             // if there is an ancester
             TParticle* grandMa = fStack->Particle(jLabel);
             Int_t grandMaPDG = grandMa->GetPdgCode();

             for (Int_t j=0; j<fNparents; j++){
                if (abs(grandMaPDG)==fParentSelect[1][j]){
                  origin = kBeautyPi0;
                  return origin;
                }
             }

             partMother = grandMa;
          } // end of iteration 

         return origin;
    } // end of if


    else {
        origin = kElse;

          // iterate until you find B hadron as a mother or become top ancester 
          for (Int_t i=1; i<100; i++){ // check back to the 100 generation older

             Int_t jLabel = partMother->GetFirstMother();
             if (jLabel == -1){
               origin = kElse;
               return origin;
             }
             if (jLabel < 0){ // safety protection
               AliDebug(1, "Stack label is negative, return\n");
               return -1;
             }

             // if there is an ancester
             TParticle* grandMa = fStack->Particle(jLabel);
             Int_t grandMaPDG = grandMa->GetPdgCode();

             for (Int_t j=0; j<fNparents; j++){
                if (abs(grandMaPDG)==fParentSelect[1][j]){
                  origin = kBeautyElse;
                  return origin;
                }
             }

             partMother = grandMa;
          } // end of iteration 

    }

    return origin;

}

/*
//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::GetElectronSource(Int_t iTrackLabel) 
{

        //           
        // decay electron's origin 
        //           
             
        if (iTrackLabel < 0) {
                AliDebug(1, "Stack label is negative, return\n");
                return -1; 
        }           

        TParticle* mcpart = fStack->Particle(iTrackLabel);
        Int_t iLabel = mcpart->GetFirstMother();
        if (iLabel<0){
                AliDebug(1, "Stack label is negative, return\n");
                return -1;
        }

        Int_t origin = -1;
        Bool_t isFinalOpenCharm = kFALSE;

        TParticle *partMother = fStack->Particle(iLabel);
        Int_t maPdgcode = partMother->GetPdgCode(); // get mother's pdg code

        //beauty --------------------------
        if ( int(abs(maPdgcode)/100.) == kBeauty || int(abs(maPdgcode)/1000.) == kBeauty ) {
                for (Int_t i=0; i<fNparents; i++){
                        if (abs(maPdgcode)==fParentSelect[1][i]){
                                origin = kDirectBeauty;
                                return origin;
                        }
                        else return -1; // this track is originated beauties not in the final B hadron list => excited beauties
                }
        } // end of if

        //charm --------------------------
        else if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) {
        
                for (Int_t i=0; i<fNparents; i++){
                        if (abs(maPdgcode)==fParentSelect[0][i])
                                isFinalOpenCharm = kTRUE;
                }
                if (!isFinalOpenCharm) return -1; // this track is originated charms not in the final D hadron list => excited charms   
                                                  // to prevent any possible double counting  

                for (Int_t i=0; i<100; i++){ // iterate 100 until you find B hadron as a mother or become top ancester

                        Int_t jLabel = partMother->GetFirstMother();
                        if (jLabel == -1){
                                origin = kDirectCharm;
                                return origin;
                        }
                        if (jLabel < 0){ // safety protection even though not really necessary here
                                AliDebug(1, "Stack label is negative, return\n");
                                return -1;
                        }

                        // if there is an ancester, check if it in the final B hadron list 
                        TParticle* grandMa = fStack->Particle(jLabel);
                        Int_t grandMaPDG = grandMa->GetPdgCode();

                        for (Int_t j=0; j<fNparents; j++){
                                if (abs(grandMaPDG)==fParentSelect[1][j]){
                                origin = kBeautyCharm;
                                return origin;
                                }
                        } 

                        partMother = grandMa;
                } // end of iteration 
        } // end of if

        //gamma --------------------------
        else if ( abs(maPdgcode) == 22 ) {
                origin = kGamma;

                // iterate until you find B hadron as a mother or become top ancester 
                for (Int_t i=0; i<100; i++){ // check back to the 100 generation older

                        Int_t jLabel = partMother->GetFirstMother();
                        if (jLabel == -1){
                                origin = kGamma;
                                return origin;
                        }
                        if (jLabel < 0){ // safety protection
                                AliDebug(1, "Stack label is negative, return\n");
                                return -1;
                        }

                        // if there is an ancester
                        TParticle* grandMa = fStack->Particle(jLabel);
                        Int_t grandMaPDG = grandMa->GetPdgCode();

                        for (Int_t j=0; j<fNparents; j++){
                                if (abs(grandMaPDG)==fParentSelect[1][j]){
                                        origin = kBeautyGamma;
                                        return origin;
                                }
                        }

                        partMother = grandMa;
                } // end of iteration 

                return origin;
        } // end of if

        //pi0 --------------------------
        else if ( abs(maPdgcode) == 111 ) {
                origin = kPi0;

                // iterate until you find B hadron as a mother or become top ancester 
                for (Int_t i=0; i<100; i++){ // check back to the 100 generation older

                        Int_t jLabel = partMother->GetFirstMother();
                        if (jLabel == -1){
                                origin = kPi0;
                                return origin;
                        }
                        if (jLabel < 0){ // safety protection
                                AliDebug(1, "Stack label is negative, return\n");
                                return -1;
                        }

                        // if there is an ancester
                        TParticle* grandMa = fStack->Particle(jLabel);
                        Int_t grandMaPDG = grandMa->GetPdgCode();

                        for (Int_t j=0; j<fNparents; j++){
                                if (abs(grandMaPDG)==fParentSelect[1][j]){
                                        origin = kBeautyPi0;
                                        return origin;
                                }
                        }

                        partMother = grandMa;
                } // end of iteration 

                return origin;
        } // end of if


        //else --------------------------
        else {
                origin = kElse;

                // iterate until you find B hadron as a mother or become top ancester 
                for (Int_t i=0; i<100; i++){ // check back to the 100 generation older

                        Int_t jLabel = partMother->GetFirstMother();
                        if (jLabel == -1){
                                origin = kElse;
                                return origin;
                        }
                        if (jLabel < 0){ // safety protection
                                AliDebug(1, "Stack label is negative, return\n");
                                return -1;
                        }

                        // if there is an ancester
                        TParticle* grandMa = fStack->Particle(jLabel);
                        Int_t grandMaPDG = grandMa->GetPdgCode();

                        for (Int_t j=0; j<fNparents; j++){
                                if (abs(grandMaPDG)==fParentSelect[1][j]){
                                        origin = kBeautyElse;
                                        return origin;
                                }
                        }

                        partMother = grandMa;
                } // end of iteration 

        }

        return origin;

}
*/

//_______________________________________________________________________________________________
Bool_t AliHFEsecVtx::SingleTrackCut(AliESDtrack* track) const
{
        // test cuts 

        //if (track->Pt() < 1.0) return kFALSE;
        //if (TMath::Abs(track->Eta()) > 0.9) return kFALSE;
        //if (!(track->GetStatus() & AliESDtrack::kITSrefit)) return kFALSE;
        //if (!(track->GetStatus() & AliESDtrack::kTPCrefit)) return kFALSE;
        if (!(TESTBIT(track->GetITSClusterMap(),0))) return kFALSE; // ask hit on the first pixel layer
        //if (!(TESTBIT(track->GetITSClusterMap(),0) | TESTBIT(track->GetITSClusterMap(),1))) return kFALSE;


/*
        Float_t dcaR=-1;
        Float_t dcaZ=-1;
        track->GetImpactParameters(dcaR,dcaZ);
        if (TMath::Abs(TMath::Sqrt(dcaR*dcaR + dcaZ*dcaZ)) < 0.005) return kFALSE;
        if (TMath::Abs(TMath::Sqrt(dcaR*dcaR + dcaZ*dcaZ)) > 0.3) return kFALSE;
*/
        return kTRUE;
}
