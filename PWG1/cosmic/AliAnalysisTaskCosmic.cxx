#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TList.h"
#include "TTree.h"
#include "TParticle.h"
#include "TParticlePDG.h"


#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliAnalysisTaskCosmic.h"
#include "AliESDInputHandler.h"

ClassImp(AliAnalysisTaskCosmic)

//________________________________________________________________________
AliAnalysisTaskCosmic::AliAnalysisTaskCosmic(const char *name) 
    : AliAnalysisTaskSE(name),
      fHists(0),
      fhDZvsZ(0),
      fhDZvsPhi(0),
      fhCh1Ch2(0),
      fhPh1Ph2(0),
      fhCl1Cl2G(0),
      fhCl1Cl2B(0)
{
  //
  // Constructor
  //
    for (Int_t i = 0; i < 6; i++) {
	fhPt[i]      = 0;
	fhTheta[i]   = 0;
	fhPhi[i]     = 0;
	fhDPhi[i]    = 0;
	fhDTheta[i]  = 0;
	fhDZ[i]      = 0;
	fhDX[i]      = 0;
	fhDY[i]      = 0;
	fhDPt[i]     = 0;
	fhD1ovPt[i]  = 0;
	fpDPt[i]     = 0;
	fhDPtovPt[i] = 0;
	fpDPtS[i]    = 0;
    }
    
  DefineOutput(1,  TList::Class());
}


//________________________________________________________________________
void AliAnalysisTaskCosmic::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
    TString ext[6] = {"PC", "NC", "PZ", "NZ", "_Good", "_Bad"};
    
    char name[12];
    
    fHists = new TList();

    for (Int_t i = 0; i < 6; i++) {
	// Pt
	sprintf(name, "fhPt%2s", ext[i].Data());
	fhPt[i]   = new TH1F(name,   " pT distribution",        800, 0., 200.);
	fhPt[i]->SetXTitle("p_{T} [Gev]");
	// Phi
	sprintf(name, "fhPhi%2s", ext[i].Data());
	fhPhi[i]   = new TH1F(name,   "Phi distribution",        62,  0., 2. * TMath::Pi());
	fhPhi[i]->SetXTitle("#phi [rad]");
	// Theta
	sprintf(name, "fhTheta%2s", ext[i].Data());
	fhTheta[i]   = new TH1F(name,   "Theta distribution",    62,  0., TMath::Pi());    
	fhTheta[i]->SetXTitle("#theta [rad]");
	// Delta Phi
	sprintf(name, "fhDPhi%2s", ext[i].Data());
	fhDPhi[i]   = new TH1F(name,   "DeltaPhi distribution",        320, -0.4, 0.4);
	fhDPhi[i]->SetXTitle("#Delta#phi [rad]");
	// Delta Theta
	sprintf(name, "fhDTheta%2s", ext[i].Data());
	fhDTheta[i]   = new TH1F(name,   "DeltaTheta distribution",    320, -0.4, 0.4);    
	fhDTheta[i]->SetXTitle("#Delta#theta [rad]");
	// Delta Z
	sprintf(name, "fhDZ%2s", ext[i].Data());
	fhDZ[i]   = new TH1F(name,   "DeltaZ distribution",        200, -10., 10.);
	fhDZ[i]->SetXTitle("#DeltaZ [cm]");
	// Delta X
	sprintf(name, "fhDX%2s", ext[i].Data());
	fhDX[i]   = new TH1F(name,   "DeltaX distribution",        200, -10., 10.);
	fhDX[i]->SetXTitle("#DeltaX [cm]");
	// Delta Y
	sprintf(name, "fhDY%2s", ext[i].Data());
	fhDY[i]   = new TH1F(name,   "DeltaY distribution",        200, -10, 10.);
	fhDY[i]->SetXTitle("#DeltaY [cm]");
	// Delta Pt
	sprintf(name, "fhDPt%2s", ext[i].Data());
	fhDPt[i]   = new TH1F(name,   "DeltaPt distribution",        200, -20., 20.);
	fhDPt[i]->SetXTitle("#Delta p_{T}  [GeV]");

	// Delta 1/Pt
	sprintf(name, "fhD1ovPt%2s", ext[i].Data());
	fhD1ovPt[i]   = new TH1F(name,   "Delta 1/Pt distribution",        200, -1., 1.);
	fhD1ovPt[i]->SetXTitle("#Delta 1/Pt");

	// Delta Pt over Pt
	sprintf(name, "fhDPtovPt%2s", ext[i].Data());
	fhDPtovPt[i]   = new TH1F(name,   "DeltaPt/Pt distribution",        200, -2., 2.);
	fhDPtovPt[i]->SetXTitle("#DeltaPt/Pt");



	// Delta Pt/ Pt vs Pt
	sprintf(name, "fpDPt%2s", ext[i].Data());
	fpDPt[i] = new TProfile(name, "#Delta Pt / Pt", 20, 0., 20., -1, 1., "S");
	fpDPt[i]->SetXTitle("p_{T} [GeV]");
	fpDPt[i]->SetYTitle("#Delta 1/p_{T} [GeV^{-1}]");
	// Delta Pt error
	sprintf(name, "fpDPtS%2s", ext[i].Data());
	fpDPtS[i] = new TProfile(name, "#Delta Pt / <sigma>", 20, 0., 20., 0., 10.);
	fpDPtS[i]->SetXTitle("p_{T}");
	fpDPtS[i]->SetYTitle("#Delta p_{T} / <#sigma_{p_{T}}>");


	fHists->Add(fhPt[i]);
	fHists->Add(fhPhi[i]);
	fHists->Add(fhTheta[i]);
	fHists->Add(fhDPhi[i]);
	fHists->Add(fhDPt[i]);
	fHists->Add(fhD1ovPt[i]);
	fHists->Add(fhDPtovPt[i]);
	fHists->Add(fhDTheta[i]);
	fHists->Add(fhDZ[i]);
	fHists->Add(fhDX[i]);
	fHists->Add(fhDY[i]);
	fHists->Add(fpDPt[i]);
	fHists->Add(fpDPtS[i]);
    }

    fhDZvsZ      = new TH2F("fhDZvsZ",    "dz vs z", 500, -250., 250., 100, -10., 10.);
    fhDZvsZ->SetXTitle("z_{in} * sign(z_{in}) * sign(z_{out}) [cm]");
    fhDZvsZ->SetYTitle("#Delta z [cm]");
    
    
    fhDZvsPhi    = new TH2F("fhDZvsPhi", "dz vs phi", 36, 0., TMath::Pi(),  50, -2., 2.);

    fhCh1Ch2   = new TH2F("fCh1Ch2",   "ch1 vs ch2", 8, -2., 2., 8, -2., 2.);
    fhPh1Ph2   = new TH2F("fPh1Ph2",   "ph1 vs ph2", 128, 0., 2. * TMath::Pi(), 128, 0., 2. * TMath::Pi());
    fhCl1Cl2G  = new TH2F("fCl1Cl2G",   "#cl vs #cl", 200, 0., 200., 200, 0., 200);
    fhCl1Cl2B  = new TH2F("fCl1Cl2B",   "#cl vs #cl", 200, 0., 200., 200, 0., 200);    

    fHists->Add(fhDZvsZ);
    fHists->Add(fhDZvsPhi);
    fHists->Add(fhCh1Ch2);
    fHists->Add(fhPh1Ph2);
    fHists->Add(fhCl1Cl2G);
    fHists->Add(fhCl1Cl2B);
    fHists->SetOwner();

}

//________________________________________________________________________
void AliAnalysisTaskCosmic::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (!fInputEvent) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  AliESDEvent* esdE = (AliESDEvent*)fInputEvent;
  if (esdE->GetNumberOfTracks() != 2) return;
  
  for (Int_t iTracks = 0; iTracks < esdE->GetNumberOfTracks(); iTracks++) {
      AliESDtrack* track = esdE->GetTrack(iTracks);


      
//      const AliExternalTrackParam * track = trackG->GetTPCInnerParam();
//      if (!track) continue;
      // Pt cut
      Float_t pt     = track->Pt();
      Float_t charge = track->Charge();
      Float_t phi    = track->Phi();
      if (phi > 0. && phi  < TMath::Pi()) charge*=(-1.);



      // TPC track
      UInt_t status = track->GetStatus();
      if (status&AliESDtrack::kTPCrefit ==0) continue;

      Int_t nClustersTPC = track->GetTPCclusters(0);
      if (nClustersTPC < 50) continue;
      
      // 

      Float_t z      = track->GetZ();
      Float_t x      = track->Xv();
      Float_t y      = track->Yv();
      Float_t theta  = track->Theta();


      const AliExternalTrackParam * trackOut = track->GetOuterParam();
      Float_t zOut = 0.;
      if (trackOut)zOut = trackOut->Zv();
      

      if (charge > 0) {
	  fhPt[kPosC]   ->Fill(pt);
      } else {
	  fhPt[kNegC]   ->Fill(pt);
      }

      //      if ((TMath::Abs(esdE->GetCurrentL3()) > 1.e-3) && (pt < 1. || pt > 50.)) continue;

      
      //
      if (charge > 0) {
	  fhPhi[kPosC]  ->Fill(phi);
	  fhTheta[kPosC]->Fill(theta);
      } else {
	  fhPhi[kNegC]  ->Fill(phi);
	  fhTheta[kNegC]->Fill(theta);
      }

      
      if (z > 0) {
	  fhPt[kPosZ]   ->Fill(pt);
	  fhPhi[kPosZ]  ->Fill(phi);
	  fhTheta[kPosZ]->Fill(theta);
      } else {
	  fhPt[kNegZ]   ->Fill(pt);
	  fhPhi[kNegZ]  ->Fill(phi);
	  fhTheta[kNegZ]->Fill(theta);
      }



      // Tracks coming from above
      if (phi > 0. &&  phi < TMath::Pi()) {
	  
//	  printf("Track#1 %5d %5d %13.3f %13.3f \n", Int_t(Entry()), iTracks, phi, zOut);
	  
	  Float_t dphiMin = 999.;
	  Float_t rMin    = 999.;
	  Int_t   jMin    = -1;
	  // Search for a matching track (in dphi-dtheta)
	  for (Int_t jTracks = 0; jTracks < esdE->GetNumberOfTracks(); jTracks++) {
	      if (jTracks == iTracks) continue;

	      AliESDtrack* track2 = esdE->GetTrack(jTracks);

	      UInt_t status = track2->GetStatus();
	      if (status&AliESDtrack::kTPCrefit ==0) continue;

	      Int_t nClustersTPC2 = track2->GetTPCclusters(0);
	      if (nClustersTPC2 < 50) continue;
	      //	      if ((TMath::Abs(esdE->GetCurrentL3()) > 1.e-3) && (track2->Pt() < 1. || track2->Pt() > 50.)) continue;
//	      const AliExternalTrackParam * track2 = trackG2->GetTPCInnerParam();
//	      if (!track2) continue;
	      
	      Float_t phi2   = track2->Phi() - TMath::Pi();
	      Float_t theta2 = TMath::Pi()   - track2->Theta();
	      Float_t dphi   = phi2   - phi;
	      Float_t dtheta = theta2 - theta;
	      
	      if (dphi >  TMath::Pi()) dphi -= 2. * TMath::Pi();
	      if (dphi < -TMath::Pi()) dphi += 2. * TMath::Pi();

	      Float_t dR = TMath::Sqrt(dphi * dphi + dtheta * dtheta);

	      if (dR < rMin) {
		  rMin    = dR;
		  dphiMin = dphi;
		  jMin = jTracks;
	      }

	  } // tracks 2
	  if (jMin != -1) {
	      // we found a matching track candidate ...
	      AliESDtrack* track2 = esdE->GetTrack(jMin);
	      const AliExternalTrackParam * trackOut2 = track2->GetOuterParam();
	      Float_t zOut2 = 0.;
	      if (trackOut2) zOut2 = trackOut2->Zv();


	      Float_t theta2 = - track2->Theta() + TMath::Pi();
	      Float_t z2 = track2->GetZ();
	      Float_t x2 = track2->Xv();
	      Float_t y2 = track2->Yv();
	      Float_t charge2 = track2->Charge();
	      Float_t dz  = z2 - z;
	      Float_t dx  = x2 - x;
	      Float_t dy  = y2 - y;
	      Float_t pt1 = track->Pt();
	      Float_t pt2 = track2->Pt();
	      Float_t dpt = pt2 - pt1;
	      Float_t d1pt = 1./pt2 - 1./pt1;
	      
	      Float_t ptm = 0.5 * (pt1 + pt2);

	      Int_t nClustersTPC2 = track2->GetTPCclusters(0);

	      if (charge > 0.) {
		  fhDPhi[kPosC]  ->Fill(dphiMin);
		  fhDTheta[kPosC]->Fill(theta2 - theta);
	      } else {
		  fhDPhi[kNegC]  ->Fill(dphiMin);
		  fhDTheta[kNegC]->Fill(theta2 - theta);
	      }

	      if (z > 0.) {
		  fhDPhi[kPosZ]  ->Fill(dphiMin);
		  fhDTheta[kPosZ]->Fill(theta2 - theta);
	      } else {
		  fhDPhi[kNegZ]  ->Fill(dphiMin);
		  fhDTheta[kNegZ]->Fill(theta2 - theta);
	      }

	      if (TMath::Abs(dpt)/ptm > 0.5) {
		  fhDPhi[kBad]  ->Fill(dphiMin);
		  fhDTheta[kBad]->Fill(theta2 - theta);
	      } else {
		  fhDPhi[kGood]  ->Fill(dphiMin);
		  fhDTheta[kGood]->Fill(theta2 - theta);
	      }
	      // Good matches ...	      
//	      if (TMath::Abs(rMin < 0.04) && (charge == charge2) && TMath::Abs(dz) < 0.5 &&  TMath::Abs(dy) && TMath::Abs(dx) < 0.5 ) 
//	      if (TMath::Abs(rMin < 0.04) && (charge == charge2) && (zOut * zOut2) < 0.) 
	      if (TMath::Abs(rMin < 0.04) && (charge == charge2)) 
	      {

		  if (TMath::Abs(dpt)/ptm  < .1) {
		      fhCl1Cl2G->Fill(nClustersTPC, nClustersTPC2);
		  }
		  
		  if (TMath::Abs(dpt)/ptm  > .5) fhCl1Cl2B->Fill(nClustersTPC, nClustersTPC2);
		  fhPh1Ph2->Fill(track->Phi(), track2->Phi());


		  Double_t sigmaPt1 = TMath::Sqrt(track->GetSigma1Pt2());
		  Double_t sigmaPt2 = TMath::Sqrt(track->GetSigma1Pt2());
		  Double_t sigmaPt = 0.5 * TMath::Sqrt(sigmaPt1 * sigmaPt1 + sigmaPt2 * sigmaPt2);
		  if (TMath::Abs(dpt)/ptm > 0.2 && pt > 10. &&  charge > 0.)
//		      printf("Track#2 %5d %5d %13.3f %13.3f %13.3f %13.3f\n", Entry(), jMin, track2->Phi(), dz, dpt, zOut2);	      		  
		  if (zOut != 0. && zOut2 != 0.) fhDZvsZ->Fill(TMath::Abs(zOut) * zOut/zOut2, dz);
		  if (zOut * zOut2 > 0. && zOut < 0) fhDZvsPhi->Fill(phi, dz);
		  
		  fhCh1Ch2->Fill(charge, track2->Charge());
		  
		  
		  if (charge > 0.) {
		      fhDPt[kPosC]->Fill(dpt);
		      fhD1ovPt[kPosC]->Fill(d1pt);
		      fpDPt[kPosC]->Fill(ptm, d1pt);
		      if (ptm > 5. && ptm < 15.) fhDPtovPt[kPosC]->Fill(dpt/ptm);
		      fhDZ[kPosC]->Fill(dz);
		      fhDX[kPosC]->Fill(dx);
		      fhDY[kPosC]->Fill(dy);
		      fpDPtS[kPosC]->Fill(ptm, 2. * TMath::Abs(1./pt1 - 1./pt2)/sigmaPt);
		  } else {
		      fhDPt[kNegC]->Fill(dpt);
		      fhD1ovPt[kNegC]->Fill(d1pt);
		      fpDPt[kNegC]->Fill(ptm, d1pt);
		      if (ptm > 5. && ptm < 15.) fhDPtovPt[kNegC]->Fill(dpt/ptm);
		      fhDZ[kNegC]->Fill(dz);
		      fhDX[kNegC]->Fill(dx);
		      fhDY[kNegC]->Fill(dy);
		      fpDPtS[kNegC]->Fill(ptm, 2. * TMath::Abs(1./pt1 - 1./pt2)/sigmaPt);
		  }

		  if (z > 0.) {
		      fhDPt[kPosZ]->Fill(dpt);
		      fhD1ovPt[kPosZ]->Fill(d1pt);
		      fpDPt[kPosZ]->Fill(ptm, d1pt);
		      fhDZ[kPosZ]->Fill(dz);
		      fhDX[kPosZ]->Fill(dx);
		      fhDY[kPosZ]->Fill(dy);
		      fpDPtS[kPosZ]->Fill(ptm, 2. * TMath::Abs(1./pt1 - 1./pt2)/sigmaPt);
		  } else {
		      fhDPt[kNegZ]->Fill(dpt);
		      fhD1ovPt[kNegZ]->Fill(d1pt);
		      fpDPt[kNegZ]->Fill(ptm, d1pt);
		      fhDZ[kNegZ]->Fill(dz);
		      fhDX[kNegZ]->Fill(dx);
		      fhDY[kNegZ]->Fill(dy);
		      fpDPtS[kNegZ]->Fill(ptm, 2. * TMath::Abs(1./pt1 - 1./pt2)/sigmaPt);
		  } 


		  if (TMath::Abs(dpt)/ptm > 0.5) {
		      fhDPt[kBad]->Fill(dpt);
		      fpDPt[kBad]->Fill(ptm, TMath::Abs(dpt)/ptm);
		      fhDZ[kBad]->Fill(dz);
		      fhDX[kBad]->Fill(dx);
		      fhDY[kBad]->Fill(dy);
		      fpDPtS[kBad]->Fill(ptm, 2. * TMath::Abs(1./pt1 - 1./pt2)/sigmaPt);
		  } else {
		      fhDPt[kGood]->Fill(dpt);
		      fpDPt[kGood]->Fill(ptm, TMath::Abs(dpt)/ptm);
		      fhDZ[kGood]->Fill(dz);
		      fhDX[kGood]->Fill(dx);
		      fhDY[kGood]->Fill(dy);
		      fpDPtS[kGood]->Fill(ptm, 2. * TMath::Abs(1./pt1 - 1./pt2)/sigmaPt);
		  } 
		  
	      } // good matches
	  } // found possible match
      } // upper
  } // tracks 1
  PostData(1, fHists);
}      

//________________________________________________________________________
void AliAnalysisTaskCosmic::Terminate(Option_t *) 
{
}  
