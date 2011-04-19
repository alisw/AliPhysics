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
//
// QA class of primary vertex study for Heavy Flavor electrons
// this has functionality to reject electrons from primary vertex
// and check primary vertex characteristics
//  
// Authors:
//   MinJung Kweon <minjung@physi.uni-heidelberg.de>
//


#include <TH2F.h>
#include <TParticle.h>

#include "AliESDEvent.h"
#include <AliESDtrack.h>
#include <AliMCEvent.h>

#include <AliLog.h>
#include <AliKFParticle.h>
#include <AliKFVertex.h>

#include "AliHFEpriVtx.h"


ClassImp(AliHFEpriVtx)

//_______________________________________________________________________________________________
AliHFEpriVtx::AliHFEpriVtx():
        fESD1(0x0)
        ,fMCEvent(0x0)
        ,fNtrackswoPid(0)
        ,fHNtrackswoPid(0x0)
        ,fNESDprimVtxContributor(0x0)
        ,fNESDprimVtxIndices(0x0)
        ,fDiffDCAvsPt(0x0)
        ,fDiffDCAvsNt(0x0)
        ,fNsectrk2prim(0)
        ,fPVxRe(-999.)
        ,fPVyRe(-999.)
        ,fPVzRe(-999.)
{ 
        //
        // Default constructor
        //

        Init();

}

//_______________________________________________________________________________________________
AliHFEpriVtx::AliHFEpriVtx(const AliHFEpriVtx &p):
         TObject(p)
        ,fESD1(0x0)
        ,fMCEvent(0x0)
        ,fNtrackswoPid(p.fNtrackswoPid)
        ,fHNtrackswoPid(0x0)
        ,fNESDprimVtxContributor(0x0)
        ,fNESDprimVtxIndices(0x0)
        ,fDiffDCAvsPt(0x0)
        ,fDiffDCAvsNt(0x0)
        ,fNsectrk2prim(p.fNsectrk2prim)
        ,fPVxRe(p.fPVxRe)
        ,fPVyRe(p.fPVyRe)
        ,fPVzRe(p.fPVzRe)
{
        //
        // Copy constructor
        //
}

//_______________________________________________________________________________________________
AliHFEpriVtx&
AliHFEpriVtx::operator=(const AliHFEpriVtx &)
{
        //
        // Assignment operator
        //

        AliInfo("Not yet implemented.");
        return *this;
}

//_______________________________________________________________________________________________
AliHFEpriVtx::~AliHFEpriVtx()
{
        //
        // Destructor
        //

        AliInfo("Analysis Done.");
}

//__________________________________________
void AliHFEpriVtx::Init()
{
        //
        // initialize counters
        //

        fNtrackswoPid = 0;
        for (int i=0; i<10; i++){       
                fPrimVtx[i].fNtrackCount = 0 ;
                fPrimVtx[i].fNprimVtxContributorCount = 0 ;
        }
}

//_______________________________________________________________________________________________
void AliHFEpriVtx::CreateHistograms(TString hnopt)
{ 
        //
        // create histograms
        //

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

           hname=hnopt+"ntracks_"+fkSourceLabel[isource];
           fPrimVtx[isource].fNtracks = new TH1F(hname,hname,50,0,50);
           hname=hnopt+"nPrimVtxContributor_"+fkSourceLabel[isource];
           fPrimVtx[isource].fNprimVtxContributor = new TH1F(hname,hname,100,0,100);
           hname=hnopt+"PtElec_"+fkSourceLabel[isource];
           fPrimVtx[isource].fPtElec = new TH1F(hname,hname,250,0,50);
           hname=hnopt+"PtElecContributor_"+fkSourceLabel[isource];
           fPrimVtx[isource].fPtElecContributor = new TH1F(hname,hname,250,0,50);

        }

        hname=hnopt+"ntrackswopid";
        fHNtrackswoPid = new TH1F(hname,hname,50,0,50);
        hname=hnopt+"nESDprimVtxContributor";
        fNESDprimVtxContributor = new TH1I(hname,hname,100,0,100);
        hname=hnopt+"nESDprimVtxIndices";
        fNESDprimVtxIndices= new TH1I(hname,hname,100,0,100);
        hname=hnopt+"diffDCAvsPt";
        fDiffDCAvsPt = new TH2F(hname,hname,250,0,50,500,0,1);
        hname=hnopt+"diffDCAvsNt";
        fDiffDCAvsNt = new TH2F(hname,hname,100,0,100,500,0,1);

}

//_______________________________________________________________________________________________
void AliHFEpriVtx::CountNtracks(Int_t sourcePart, Int_t recpid, Double_t recprob)
{
        //
        // count number of tracks passed certain cuts
        //

        fNtrackswoPid++;

        if (!recpid && recprob>0.5)     
        fPrimVtx[kAll].fNtrackCount++;
        if(sourcePart<0) return;
        fPrimVtx[sourcePart].fNtrackCount++;

}

//_______________________________________________________________________________________________
void AliHFEpriVtx::FillNtracks()
{
        //
        // count number of tracks passed certain cuts
        //

        fHNtrackswoPid->Fill(fNtrackswoPid);
        for (int i=0; i<10; i++){
          fPrimVtx[i].fNtracks->Fill(fPrimVtx[i].fNtrackCount);
        }

}

//_______________________________________________________________________________________________
Int_t AliHFEpriVtx::GetMCPID(AliESDtrack const *track)
{
        //
        // get MC pid
        //

	AliMCParticle *mctrack = NULL;
	if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(track->GetLabel()))))) return 0; 
	TParticle *mcpart = mctrack->Particle();

        if ( !mcpart ) return 0;
        Int_t pdgCode = mcpart->GetPdgCode();

        return pdgCode;
}

//_______________________________________________________________________________________________
void AliHFEpriVtx::GetNPriVxtContributor() 
{
        //
        // count number of primary vertex contributor
        //

        const AliESDVertex *primvtx = fESD1->GetPrimaryVertex();
        fNESDprimVtxContributor->Fill(primvtx->GetNContributors());
        fNESDprimVtxIndices->Fill(primvtx->GetNIndices());
}

//_______________________________________________________________________________________________
void AliHFEpriVtx::CountPriVxtElecContributor(AliESDtrack *ESDelectron, Int_t sourcePart, Int_t recpid, Double_t recprob) 
{
        //
        // count number of electrons contributing to the primary vertex
        //


        if (recpid || recprob<0.5) return;

        // get track id of our selected electron
        Int_t elecTrkID = ESDelectron->GetID();

	AliMCParticle *mctrack = NULL;
	if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(ESDelectron->GetLabel()))))) return; 
	TParticle *mcpart = mctrack->Particle();


        if(!mcpart){
          AliDebug(1, "no mc particle, return\n");
          return;
        }

        AliKFParticle::SetField(fESD1->GetMagneticField());
        AliKFParticle kfElectron(*ESDelectron,11);
 
        // prepare kfprimary vertex
        AliKFVertex kfESDprimary;

        // Reconstructed Primary Vertex (with ESD tracks)
        const AliESDVertex *primvtx = fESD1->GetPrimaryVertex();
//         Int_t nt = primvtx->GetNContributors();
        Int_t n=primvtx->GetNIndices();

        if (n>0 && primvtx->GetStatus()){

                kfESDprimary = AliKFVertex(*primvtx);
                Double_t dcaBFrmElec = kfElectron.GetDistanceFromVertexXY(kfESDprimary);

                UShort_t *priIndex = primvtx->GetIndices();

                fPrimVtx[kAll].fPtElec->Fill(mcpart->Pt());
                if(sourcePart>=0) fPrimVtx[sourcePart].fPtElec->Fill(mcpart->Pt());

                for (Int_t i=0;i<n;i++){

                        Int_t idx = Int_t(priIndex[i]);
                        if (idx == elecTrkID){
                                fPrimVtx[kAll].fNprimVtxContributorCount++;
                                fPrimVtx[kAll].fPtElecContributor->Fill(mcpart->Pt());
                                if(sourcePart<0) continue;
                                fPrimVtx[sourcePart].fNprimVtxContributorCount++;
                                fPrimVtx[sourcePart].fPtElecContributor->Fill(mcpart->Pt());

                                kfESDprimary -= kfElectron;
                                Double_t dcaAFrmElec = kfElectron.GetDistanceFromVertexXY(kfESDprimary);
                                fDiffDCAvsPt->Fill(mcpart->Pt(),dcaBFrmElec-dcaAFrmElec);
                                fDiffDCAvsNt->Fill(n,dcaBFrmElec-dcaAFrmElec);
                        }
                } 
        }  

}

//_______________________________________________________________________________________________
void AliHFEpriVtx::FillNprimVtxContributor() const
{
        //
        // Fill histogram with number of electrons contributing to the primary vertex
        //

        for (int i=0; i<10; i++){
          fPrimVtx[i].fNprimVtxContributor->Fill(fPrimVtx[i].fNprimVtxContributorCount);
        }

}

//_______________________________________________________________________________________________
Double_t AliHFEpriVtx::GetDistanceFromRecalVertexXY(const AliESDtrack * const ESDelectron) 
{
        //
        // return recalculated DCA after removing input track from the primary vertex
        //

        // get track id of our selected electron
        Int_t elecTrkID = ESDelectron->GetID();

        AliKFParticle::SetField(fESD1->GetMagneticField());
        AliKFParticle kfElectron(*ESDelectron,11);
 
        // prepare kfprimary vertex
        AliKFVertex kfESDprimary;

        // Reconstructed Primary Vertex (with ESD tracks)
        const AliESDVertex *primvtx = fESD1->GetPrimaryVertex();
        Int_t n=primvtx->GetNIndices();

        if (n>0 && primvtx->GetStatus()){

                kfESDprimary = AliKFVertex(*primvtx);
                UShort_t *priIndex = primvtx->GetIndices();

                for (Int_t i=0;i<n;i++){

                        Int_t idx = Int_t(priIndex[i]);
                        if (idx == elecTrkID){
                                kfESDprimary -= kfElectron;
                                Double_t dcaAFrmElec = kfElectron.GetDistanceFromVertexXY(kfESDprimary);

                                return dcaAFrmElec;
                        }
                } 
        }  
	      return -1;

}

void AliHFEpriVtx::RecalcPrimvtx(Int_t nkftrk, const Int_t * const trkid, const AliKFParticle * const kftrk)
{
        //
        // recalculate primary vertex after removing the input track
        //

        const AliESDVertex *primvtx = fESD1->GetPrimaryVertex();

        AliKFVertex kfESDprimary;
        Int_t n = primvtx->GetNIndices();
        fNsectrk2prim = 0;
        fPVxRe = -999.;
        fPVyRe = -999.;
        fPVyRe = -999.;

        if (n>0 && primvtx->GetStatus()){
          kfESDprimary = AliKFVertex(*primvtx);
          UShort_t *priIndex = primvtx->GetIndices();
          for (Int_t j=0; j<nkftrk; j++){
            for (Int_t i=0;i<n;i++){
              Int_t idx = Int_t(priIndex[i]);
              if (idx == trkid[j]){
                kfESDprimary -= kftrk[j];
                fNsectrk2prim++;
              }
            }
          }
        }

        fPVxRe = kfESDprimary.GetX();
        fPVyRe = kfESDprimary.GetY();
        fPVzRe = kfESDprimary.GetZ();

}


//_______________________________________________________________________________________________
void AliHFEpriVtx::RecalcPrimvtx(const AliESDtrack * const ESDelectron)
{
        //
        // recalculate primary vertex after removing the input track
        //

        // get track id of our selected electron
        Int_t elecTrkID = ESDelectron->GetID();

        AliKFParticle::SetField(fESD1->GetMagneticField());
        AliKFParticle kfElectron(*ESDelectron,11);

        const AliESDVertex *primvtx = fESD1->GetPrimaryVertex();

        AliKFVertex kfESDprimary;
        Int_t n = primvtx->GetNIndices();
        fPVxRe = -999.;
        fPVyRe = -999.;
        fPVyRe = -999.;
        
        if (n>0 && primvtx->GetStatus()){
          kfESDprimary = AliKFVertex(*primvtx);
          UShort_t *priIndex = primvtx->GetIndices();
          for (Int_t i=0;i<n;i++){
            Int_t idx = Int_t(priIndex[i]);
            if (idx == elecTrkID){
              kfESDprimary -= kfElectron;
            }
          }
        }     
        
        fPVxRe = kfESDprimary.GetX();
        fPVyRe = kfESDprimary.GetY();
        fPVzRe = kfESDprimary.GetZ();

}
