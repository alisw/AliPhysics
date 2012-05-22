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

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"

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
  ,fAlgorithm("MA")
  ,fAngleCut(999)
  ,fdcaCut(999)
  ,fdEdxMin(62)
  ,fdEdxMax(100)
  ,fMassCut(0.5)
  ,fChi2OverNDFCut(999)
  ,fIsLS(kFALSE)
  ,fIsULS(kFALSE)
  ,fNLS(0)
  ,fNULS(0)
  ,fLSPartner(0)
  ,fULSPartner(0)
  ,fHistMass(0)
  ,fHistMassBack(0)
  ,fHistDCA(0)
  ,fHistDCABack(0)
  ,fHistAngle(0)
  ,fHistAngleBack(0)
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
  ,fAlgorithm("MA")
  ,fAngleCut(999)
  ,fdcaCut(999)
  ,fdEdxMin(62)
  ,fdEdxMax(100)
  ,fMassCut(0.5)
  ,fChi2OverNDFCut(999)
  ,fIsLS(kFALSE)
  ,fIsULS(kFALSE)
  ,fNLS(0)
  ,fNULS(0)
  ,fLSPartner(0)
  ,fULSPartner(0)
  ,fHistMass(0)
  ,fHistMassBack(0)
  ,fHistDCA(0)
  ,fHistDCABack(0)
  ,fHistAngle(0)
  ,fHistAngleBack(0)
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
void AliSelectNonHFE::FindNonHFE(Int_t iTrack1, AliESDtrack *track1, AliESDEvent *fESD)
{	
  //
  // Find non HFE electrons
  //


  //Magnetic Field
  Double_t bfield = fESD->GetMagneticField();
  
  //Second Track loop
  
  fIsULS = kFALSE; 	//Non-HFE Unlike signal Flag
  fIsLS = kFALSE; 	//Non-HFE like signal Flag
  fNULS = 0; 		//Non-HFE Unlike signal Flag
  fNLS = 0; 		//Non-HFE like signal Flag
  
  if(fLSPartner) delete [] fLSPartner;
  if(fULSPartner) delete [] fULSPartner;
  fLSPartner = new int [100]; 	//store the partners index
  fULSPartner = new int [100];	//store the partners index
  
  for(Int_t iTrack2 = 0; iTrack2 < fESD->GetNumberOfTracks(); iTrack2++) 
    {
		if(iTrack1==iTrack2) continue;
		
      AliESDtrack* track2 = fESD->GetTrack(iTrack2);
      if (!track2) 
	{
	  printf("ERROR: Could not receive track %d\n", iTrack2);
	  continue;
	}
      
      //Second track cuts
      Double_t dEdx2 = track2->GetTPCsignal();
      if(dEdx2<fdEdxMin || dEdx2>fdEdxMax) continue;
      if(!fTrackCuts->AcceptTrack(track2)) continue;
      
      if(fAlgorithm=="MA")
	{
	  //Variables
	  Double_t p1[3];
	  Double_t p2[3];
	  Double_t xt1; //radial position track 1 at the DCA point
	  Double_t xt2; //radial position track 2 at the DCA point
	  //DCA track1-track2
	  Double_t dca12 = track2->GetDCA(track1,bfield,xt2,xt1);
	  
	  //Momento of the track extrapolated to DCA track-track	
	  //Track1
	  Bool_t hasdcaT1 = track1->GetPxPyPzAt(xt1,bfield,p1);
	  //Track2
	  Bool_t hasdcaT2 = track2->GetPxPyPzAt(xt2,bfield,p2);
	  
	  if(!hasdcaT1 || !hasdcaT2) AliWarning("It could be a problem in the extrapolation");
	  
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
	  if(fCharge1*fCharge2<0 && fHistMass) fHistMass->Fill(imass);
	  if(fCharge1*fCharge2>0 && fHistMassBack) fHistMassBack->Fill(imass);
	  
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
	  if(fCharge1*fCharge2<0 && fHistMass) fHistMass->Fill(imass);
	  if(fCharge1*fCharge2>0 && fHistMassBack) fHistMassBack->Fill(imass);
	  
	}
      else
	{
	  AliError( Form("Error: %s is not a valid algorithm option.",(const char*)fAlgorithm));
	  return;
	}
      
    }
  
  return;
}
