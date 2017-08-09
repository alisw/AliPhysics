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


////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Class for the Selection of Non-Heavy-Flavour-Electrons trought 	  //
//	the invariant mass method. The selection can be done from two	  //
//  diferent algorithms, which can be choosed calling the function    //
//  "SetAlgorithm(TString Algorithm)".								  //
//                                                                    //
//  		Author: Elienos Pereira de Oliveira Filho				  //
//					(University of São Paulo)               	      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TH1F.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVEvent.h"
#include "AliESDtrackCuts.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliPIDResponse.h"
#include "AliPID.h"
#include "AliExternalTrackParam.h"
#include "AliSelectNonHFE.h"
#include "AliKFParticle.h"
#include "AliLog.h"
#include "stdio.h"
#include "iostream"
#include "fstream"

ClassImp(AliSelectNonHFE)
//________________________________________________________________________
AliSelectNonHFE::AliSelectNonHFE(const char *name, const Char_t *title)
: TNamed(name, title)
,fTrackCuts(0)
,fAlgorithm("DCA")
,fAngleCut(999)
,fdcaCut(999)
,fTPCnSigmaMin(-3)
,fTPCnSigmaMax(3)
,fMassCut(0.5)
,fChi2OverNDFCut(999)
,fPtMin(0.3)
,fIsLS(kFALSE)
,fIsULS(kFALSE)
,fIsAOD(kFALSE)
,fHasPtCut(kFALSE)
,fNLS(0)
,fNULS(0)
,fTpcNcls(50)
,fweight1(1)
,fweight2(1)
,fLSPartner(0)
,fULSPartner(0)
,fHistMass(0)
,fHistMassBack(0)
,fHistDCA(0)
,fHistDCABack(0)
,fHistAngle(0)
,fHistAngleBack(0)
,fPIDResponse(0)
,fNClusITS(2)
,fUseGlobalTracks(kFALSE)
,fRequirePointOnITS(kFALSE)
,fRequireITSAndTPCRefit(kTRUE)
,fUseDCAPartnerCut(kFALSE)
,fDCAcutxyPartner(1.00)
,fDCAcutzPartner(2.00)
,fUseEtaCutForPart(kFALSE)
,fEtaCutMin(-0.8)
,fEtaCutMax(0.8)
,fRequireTPCNclusForPID(kFALSE)
,fTpcNclsPID(60)
,fUseTender(kFALSE)

{
    //
    // Constructor
    //
    
    fTrackCuts = new AliESDtrackCuts();
    //Configure Default Track Cuts
    fTrackCuts->SetAcceptKinkDaughters(kFALSE);
    fTrackCuts->SetRequireTPCRefit(kTRUE);
    fTrackCuts->SetEtaRange(-0.9,0.9);
    fTrackCuts->SetRequireSigmaToVertex(kTRUE);
    fTrackCuts->SetMaxChi2PerClusterTPC(4.0);
    fTrackCuts->SetMinNClustersTPC(50);
    fTrackCuts->SetPtRange(0.3,1e10);
    
    
}

//________________________________________________________________________
AliSelectNonHFE::AliSelectNonHFE()
: TNamed()
,fTrackCuts(0)
,fAlgorithm("DCA")
,fAngleCut(999)
,fdcaCut(999)
,fTPCnSigmaMin(-3)
,fTPCnSigmaMax(3)
,fMassCut(0.5)
,fChi2OverNDFCut(999)
,fPtMin(0.3)
,fIsLS(kFALSE)
,fIsULS(kFALSE)
,fIsAOD(kFALSE)
,fHasPtCut(kFALSE)
,fNLS(0)
,fNULS(0)
,fTpcNcls(50)
,fweight1(1)
,fweight2(1)
,fLSPartner(0)
,fULSPartner(0)
,fHistMass(0)
,fHistMassBack(0)
,fHistDCA(0)
,fHistDCABack(0)
,fHistAngle(0)
,fHistAngleBack(0)
,fPIDResponse(0)
,fNClusITS(2)
,fUseGlobalTracks(kFALSE)
,fRequirePointOnITS(kFALSE)
,fRequireITSAndTPCRefit(kTRUE)
,fUseDCAPartnerCut(kFALSE)
,fDCAcutxyPartner(1.00)
,fDCAcutzPartner(2.00)
,fUseEtaCutForPart(kFALSE)
,fEtaCutMin(-0.8)
,fEtaCutMax(0.8)
,fRequireTPCNclusForPID(kFALSE)
,fTpcNclsPID(60)
,fUseTender(kFALSE)


{
    //
    // Constructor
    //
    
    fTrackCuts = new AliESDtrackCuts();
    //Configure Default Track Cuts
    fTrackCuts->SetAcceptKinkDaughters(kFALSE);
    fTrackCuts->SetRequireTPCRefit(kTRUE);
    fTrackCuts->SetEtaRange(-0.9,0.9);
    fTrackCuts->SetRequireSigmaToVertex(kTRUE);
    fTrackCuts->SetMaxChi2PerClusterTPC(4.0);
    fTrackCuts->SetMinNClustersTPC(50);
    fTrackCuts->SetPtRange(0.3,1e10);
    
    
}

//_________________________________________
AliSelectNonHFE::~AliSelectNonHFE()
{
    //
    // Destructor
    //
    
    if(fTrackCuts) delete fTrackCuts;
    if(fLSPartner) delete [] fLSPartner;
    if(fULSPartner) delete [] fULSPartner;
    
    
}

//__________________________________________
void AliSelectNonHFE::FindNonHFE(Int_t iTrack1, AliVParticle *Vtrack1, AliVEvent *fVevent, TClonesArray  *fTracks_tender, Bool_t fUseTender)
{
    
	
	AliVTrack *track1 = dynamic_cast<AliVTrack*>(Vtrack1);
    AliESDtrack *etrack1 = dynamic_cast<AliESDtrack*>(Vtrack1);
	AliAODTrack *atrack1 = dynamic_cast<AliAODTrack*>(Vtrack1);
    
    AliExternalTrackParam extTrackParam1;
    extTrackParam1.CopyFromVTrack(track1);

    //
    // Find non HFE electrons
    //
    
    //Magnetic Field
    Double_t bfield = fVevent->GetMagneticField();
    
    //Second Track loop
    
    fIsULS = kFALSE; 	//Non-HFE Unlike signal Flag
    fIsLS = kFALSE; 	//Non-HFE like signal Flag
    fNULS = 0; 		//Non-HFE Unlike signal Flag
    fNLS = 0; 		//Non-HFE like signal Flag
    
    if(fLSPartner) delete [] fLSPartner;
    if(fULSPartner) delete [] fULSPartner;
    fLSPartner = new int [100]; 	//store the partners index
    fULSPartner = new int [100];	//store the partners index
	
	
	//=============================
	Int_t NTracks=0;
	
		
	if(fUseTender){
		NTracks = fTracks_tender->GetEntries();
	}
	else{
		NTracks=fVevent->GetNumberOfTracks();
	}
    //=============================
	
	
    for(Int_t iTrack2 = 0; iTrack2 < NTracks; iTrack2++)
    {
		
		//It will work for EMCal framework if the "fTracks_tender" and flag fUseTender is passed to SelectNonHFE!
        if(iTrack1==iTrack2) continue;
		
		AliVParticle* Vtrack2 = 0x0;
		if(!fUseTender) Vtrack2  = fVevent->GetTrack(iTrack2);
					
		if(fUseTender){
			Vtrack2 = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTrack2));
		}
		
		if (!Vtrack2)
		{
			printf("ERROR: Could not receive track %d\n", iTrack2);
			continue;
		}
		 
		
	
        
        AliVTrack *track2 = dynamic_cast<AliVTrack*>(Vtrack2);
        AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
        AliESDtrack *etrack2 = dynamic_cast<AliESDtrack*>(Vtrack2);
		
			
		
        AliExternalTrackParam extTrackParam2;
        extTrackParam2.CopyFromVTrack(track2);

        
        //Second track cuts
        if(fIsAOD)
        {
            //AOD Filter Bit
            if (fUseGlobalTracks)
            {
                if(!atrack2->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //Same as trigger
            }
            else
            {
             if(!atrack2->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue; //Old one
            }
            
            //ITS and TPC refit
            if (fRequireITSAndTPCRefit)
            {
                if((!(atrack2->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack2->GetStatus()&AliESDtrack::kTPCrefit))))
                    continue;
            }
            
            //Eta cut
            if (fUseEtaCutForPart)
                if(atrack2->Eta() < fEtaCutMin || atrack2->Eta() > fEtaCutMax)
                    continue;
            
            //NClusters on TPC
            if(atrack2->GetTPCNcls() < fTpcNcls) continue;
            
            //TPC NClusters for PID
            if (fRequireTPCNclusForPID)
                if(atrack2->GetTPCsignalN() < fTpcNclsPID) continue;
            
            
            //Number of Clusters on ITS
            if (fRequirePointOnITS)
                if (atrack2->GetITSNcls() < fNClusITS) continue;   //Add minimum number of clusters on the ITS
            
            if (fUseDCAPartnerCut)
            {
            //Calculate DCA
                Double_t d0z0[2], cov[3];
                const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
                if(atrack2->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov)){
                    if(TMath::Abs(d0z0[0]) > fDCAcutxyPartner || TMath::Abs(d0z0[1]) > fDCAcutzPartner ) continue;
                }
            }
            
        }
        else
        {
            if(!fTrackCuts->AcceptTrack(etrack2)) continue;
        }
        
        //Second track pid
        Double_t tpcNsigma2 = fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kElectron);
        if(tpcNsigma2<fTPCnSigmaMin || tpcNsigma2>fTPCnSigmaMax) continue;
        
        //Pt Cut
        if((track2->Pt() < fPtMin) && (fHasPtCut)) continue;
		
		
		
		
		
		
        
        if(fAlgorithm=="DCA")
        {
            //Variables
            Double_t p1[3];
            Double_t p2[3];
            Double_t xt1; //radial position track 1 at the DCA point
            Double_t xt2; //radial position track 2 at the DCA point
            Double_t dca12; //DCA 1-2
            Bool_t hasdcaT1;
            Bool_t hasdcaT2;
            
            if(fIsAOD)
            {
                //DCA track1-track2
                dca12 = extTrackParam2.GetDCA(&extTrackParam1,bfield,xt2,xt1);
                
                //Momento of the track extrapolated to DCA track-track
                //Track1
                hasdcaT1 = extTrackParam1.GetPxPyPzAt(xt1,bfield,p1);
                //Track2
                hasdcaT2 = extTrackParam2.GetPxPyPzAt(xt2,bfield,p2);
            }
            else
            {
                //DCA track1-track2
                dca12 = etrack2->GetDCA(etrack1,bfield,xt2,xt1);
                
                //Momento of the track extrapolated to DCA track-track
                //Track1
                hasdcaT1 = etrack1->GetPxPyPzAt(xt1,bfield,p1);
                //Track2
                hasdcaT2 = etrack2->GetPxPyPzAt(xt2,bfield,p2);
            }
            
            if(!hasdcaT1 || !hasdcaT2) AliWarning("It could be a problem in the extrapolation");
            
            //added by Cris to not take same track  for tender case (not necessary anymore, since now for emcal framework the TClones array with track and flag of tender can be passed)
			/*
            if(p1[0]==p2[0] && p1[1]==p2[1] && p1[2]==p2[2]){
				printf("Track %d was rejected when combined with main track %d", iTrack2, iTrack1);
				printf("Checking ID: id1 =%f, id2=%f \n", id1, id2);
                continue;
                
            }
			 */
			
            
            
            //track1-track2 Invariant Mass
            Double_t eMass = 0.000510998910; //Electron mass in GeV
            Double_t pP1 = sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]); //Track 1 momentum
            Double_t pP2 = sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]); //Track 1 momentum
            
            //Double_t E1 = sqrt(eMass*eMass+pP1*pP1);
            //Double_t E2 = sqrt(eMass*eMass+pP2*pP2);
            //Double_t imass = sqrt(2*eMass*eMass+2*(E1*E2-(p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2])));
            //Double_t angle = TMath::ACos((p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2])/(pP1*pP2));
            
            TLorentzVector v1(p1[0],p1[1],p1[2],sqrt(eMass*eMass+pP1*pP1));
            TLorentzVector v2(p2[0],p2[1],p2[2],sqrt(eMass*eMass+pP2*pP2));
            Double_t imass = (v1+v2).M(); //Invariant Mass
            Double_t angle = v1.Angle(v2.Vect()); //Opening Angle (Total Angle)
            
            Float_t fCharge1 = track1->Charge();
            Float_t fCharge2 = track2->Charge();
            
            if(imass<fMassCut && angle<fAngleCut && dca12<fdcaCut)
            {
                if(fCharge1*fCharge2<0)
                {
                    fIsULS=kTRUE;
                    fULSPartner[fNULS] = iTrack2;
                    fNULS++;
                }
                if(fCharge1*fCharge2>0)
                {
                    fIsLS=kTRUE;
                    fLSPartner[fNLS] = iTrack2;
                    fNLS++;
                }
            }
            
            //Fill some histograms
            if(fCharge1*fCharge2<0 && fHistMass) fHistMass->Fill(imass, fweight1);
            if(fCharge1*fCharge2>0 && fHistMassBack) fHistMassBack->Fill(imass, fweight2);
            
            if(fCharge1*fCharge2<0 && fHistDCA) fHistDCA->Fill(dca12);
            if(fCharge1*fCharge2>0 && fHistDCABack) fHistDCABack->Fill(dca12);
            
            if(fCharge1*fCharge2<0 && fHistAngle) fHistAngle->Fill(angle);
            if(fCharge1*fCharge2>0 && fHistAngleBack) fHistAngleBack->Fill(angle);
            
        }
        else if(fAlgorithm=="KF")
        {
			
			
			
			
            Int_t fPDGtrack1 = 11;
            Int_t fPDGtrack2 = 11;
            
            Float_t fCharge1 = track1->Charge();
            Float_t fCharge2 = track2->Charge();
            
            if(fCharge1>0) fPDGtrack1 = -11;
            if(fCharge2>0) fPDGtrack2 = -11;
            
            AliKFParticle::SetField(bfield);
            AliKFParticle fKFtrack1(*track1, fPDGtrack1);
            AliKFParticle fKFtrack2(*track2, fPDGtrack2);
            AliKFParticle fRecoGamma(fKFtrack1, fKFtrack2);
            
            //Reconstruction Cuts
            if(fRecoGamma.GetNDF()<1) continue;
            Double_t chi2OverNDF = fRecoGamma.GetChi2()/fRecoGamma.GetNDF();
            if(TMath::Sqrt(TMath::Abs(chi2OverNDF))>fChi2OverNDFCut) continue;
            
            //Invariant Mass
            Double_t imass; 
            Double_t width;
            fRecoGamma.GetMass(imass,width);
            
            //Opening Angle (Total Angle)
            Double_t angle = fKFtrack1.GetAngle(fKFtrack2);
            
            //Fill some histograms	
            if(fCharge1*fCharge2<0 && fHistAngle) fHistAngle->Fill(angle);
            if(fCharge1*fCharge2>0 && fHistAngleBack) fHistAngleBack->Fill(angle);
            
            if(angle>fAngleCut) continue;
            
            if(imass<fMassCut)
            {
                if(fCharge1*fCharge2<0)
                {
                    fIsULS=kTRUE;
                    fULSPartner[fNULS] = iTrack2;
                    fNULS++;
                }
                if(fCharge1*fCharge2>0)
                {
                    fIsLS=kTRUE;
                    fLSPartner[fNLS] = iTrack2;
                    fNLS++;
                }
            }
            
            //Fill some histograms
            if(fCharge1*fCharge2<0 && fHistMass) fHistMass->Fill(imass, fweight1);
            if(fCharge1*fCharge2>0 && fHistMassBack) fHistMassBack->Fill(imass, fweight2);
            
        }
        else
        {
            AliError( Form("Error: %s is not a valid algorithm option.",(const char*)fAlgorithm));
            return;
        }
        
    }
    
    
    
    
    return;
}
