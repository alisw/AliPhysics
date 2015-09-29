/**************************************************************************
 * Author: Michael Weber                                                  *
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

//-------------------------------------------------------------------------
//                          Class AliBalanceEbyE
//   This is the class for the Balance Function analysis on an EbyE basis
//
//    Origin: m.weber@cern.ch
//    Based on AliBalancePsi
//-------------------------------------------------------------------------

#include "TH2D.h"
#include "TH3D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TParticle.h"


#include "AliAnalysisTaskTriggeredBF.h"
#include "AliBalanceEbyE.h"
#include "AliLog.h"
#include "AliVParticle.h"


ClassImp(AliBalanceEbyE)

//____________________________________________________________________//
AliBalanceEbyE::AliBalanceEbyE() :
  TObject(), 
  fShuffle(kFALSE),
  fAnalysisLevel("ESD"),
  fAnalyzedEvents(0) ,
  fHistPN(0),
  fHistNP(0),
  fHistPP(0),
  fHistNN(0),
  fHistBF(0),
  fHistBFSum(0),
  fHistHBTbefore(0),
  fHistHBTafter(0),
  fHistConversionbefore(0),
  fHistConversionafter(0),
  fHistPsiMinusPhi(0),
  fHistResonancesBefore(0),
  fHistResonancesRho(0),
  fHistResonancesK0(0),
  fHistResonancesLambda(0),
  fHistQbefore(0),
  fHistQafter(0),
  fPsiInterval(15.),
  fDeltaEtaMax(2.0),
  fResonancesCut(kFALSE),
  fHBTCut(kFALSE),
  fHBTCutValue(0.02),
  fConversionCut(kFALSE),
  fInvMassCutConversion(0.04),
  fQCut(kFALSE),
  fDeltaPtMin(0.0)
{
  // Default constructor
}

//____________________________________________________________________//
AliBalanceEbyE::AliBalanceEbyE(const AliBalanceEbyE& balance):
  TObject(balance), fShuffle(balance.fShuffle), 
  fAnalysisLevel(balance.fAnalysisLevel),
  fAnalyzedEvents(balance.fAnalyzedEvents), 
  fHistPN(balance.fHistPN),
  fHistNP(balance.fHistNP),
  fHistPP(balance.fHistPP),
  fHistNN(balance.fHistNN),
  fHistBF(balance.fHistBF),
  fHistBFSum(balance.fHistBFSum),
  fHistHBTbefore(balance.fHistHBTbefore),
  fHistHBTafter(balance.fHistHBTafter),
  fHistConversionbefore(balance.fHistConversionbefore),
  fHistConversionafter(balance.fHistConversionafter),
  fHistPsiMinusPhi(balance.fHistPsiMinusPhi),
  fHistResonancesBefore(balance.fHistResonancesBefore),
  fHistResonancesRho(balance.fHistResonancesRho),
  fHistResonancesK0(balance.fHistResonancesK0),
  fHistResonancesLambda(balance.fHistResonancesLambda),
  fHistQbefore(balance.fHistQbefore),
  fHistQafter(balance.fHistQafter),
  fPsiInterval(balance.fPsiInterval),
  fDeltaEtaMax(balance.fDeltaEtaMax),
  fResonancesCut(balance.fResonancesCut),
  fHBTCut(balance.fHBTCut),
  fHBTCutValue(balance.fHBTCutValue),
  fConversionCut(balance.fConversionCut),
  fInvMassCutConversion(balance.fInvMassCutConversion),
  fQCut(balance.fQCut),
  fDeltaPtMin(balance.fDeltaPtMin)
{
  //copy constructor
}

//____________________________________________________________________//
AliBalanceEbyE::~AliBalanceEbyE() {
  // Destructor

  delete fHistPN;
  delete fHistNP;
  delete fHistPP;
  delete fHistNN;
  delete fHistBF;
  delete fHistHBTbefore;
  delete fHistHBTafter;
  delete fHistConversionbefore;
  delete fHistConversionafter;
  delete fHistPsiMinusPhi;
  delete fHistResonancesBefore;
  delete fHistResonancesRho;
  delete fHistResonancesK0;
  delete fHistResonancesLambda;
  delete fHistQbefore;
  delete fHistQafter;
    
}

//____________________________________________________________________//
void AliBalanceEbyE::InitHistograms() {
  // single particle histograms

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);


  // =========================================================
  // Create the Output objects 
  // =========================================================

  // QA histograms
  fHistHBTbefore        = new TH2D("fHistHBTbefore","before HBT cut",200,0,2,200,0,2.*TMath::Pi());
  fHistHBTafter         = new TH2D("fHistHBTafter","after HBT cut",200,0,2,200,0,2.*TMath::Pi());
  fHistConversionbefore = new TH3D("fHistConversionbefore","before Conversion cut;#Delta#eta;#Delta#phi;M_{inv}^{2}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistConversionafter  = new TH3D("fHistConversionafter","after Conversion cut;#Delta#eta;#Delta#phi;M_{inv}^{2}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistPsiMinusPhi      = new TH2D("fHistPsiMinusPhi","",4,-0.5,3.5,100,0,2.*TMath::Pi());
  fHistResonancesBefore = new TH3D("fHistResonancesBefore","before resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistResonancesRho    = new TH3D("fHistResonancesRho","after #rho resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistResonancesK0     = new TH3D("fHistResonancesK0","after #rho, K0 resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistResonancesLambda = new TH3D("fHistResonancesLambda","after #rho, K0, Lambda resonance cut;#Delta#eta;#Delta#phi;M_{inv}",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistQbefore          = new TH3D("fHistQbefore","before momentum difference cut;#Delta#eta;#Delta#phi;|#Delta p_{T}| (GeV/c)",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);
  fHistQafter           = new TH3D("fHistQafter","after momentum difference cut;#Delta#eta;#Delta#phi;|#Delta p_{T}| (GeV/c)",50,-2.0,2.0,50,-TMath::Pi()/2.,3.*TMath::Pi()/2.,300,0,1.5);


  // EbyE BF histograms
  fHistPN = new TH2F("fHistPN","fHistPN",40,-1.6,1.6,72,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fHistNP = new TH2F("fHistNP","fHistNP",40,-1.6,1.6,72,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fHistPP = new TH2F("fHistPP","fHistPP",40,-1.6,1.6,72,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fHistNN = new TH2F("fHistNN","fHistNN",40,-1.6,1.6,72,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fHistBF = new TH2F("fHistBF","fHistBF",40,-1.6,1.6,72,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fHistBFSum = new TH3F("fHistBFSum","fHistBFSum",20,0,100,40,-1.6,1.6,72,-0.5*TMath::Pi(),1.5*TMath::Pi());
  //fHistBFSum = new TH2D("fHistBFSum","fHistBFSum",40,-1.6,1.6,72,-0.5*TMath::Pi(),1.5*TMath::Pi());


  TH1::AddDirectory(oldStatus);

}

//____________________________________________________________________//
void AliBalanceEbyE::CalculateBalance(Double_t gReactionPlane,
				     TObjArray *particles, 
				     TObjArray *particlesMixed,
				     Float_t bSign,
				     Double_t kMultorCent,
				     Double_t vertexZ) {

  // Calculates the balance function
  fAnalyzedEvents++;

  Double_t trackVariablesSingle[3];
  Double_t trackVariablesPair[6];

  if (!particles){
    AliWarning("particles TObjArray is NULL pointer --> return");
    return;
  }
  
  // define end of particle loops
  Int_t iMax = particles->GetEntriesFast();
  Int_t jMax = iMax;
  if (particlesMixed)
    jMax = particlesMixed->GetEntriesFast();

  // Eta() is extremely time consuming, therefore cache it for the inner loop here:
  TObjArray* particlesSecond = (particlesMixed) ? particlesMixed : particles;

  TArrayF secondEta(jMax);
  TArrayF secondPhi(jMax);
  TArrayF secondPt(jMax);
  TArrayS secondCharge(jMax);
  TArrayD secondCorrection(jMax);

  for (Int_t i=0; i<jMax; i++){
    secondEta[i] = ((AliVParticle*) particlesSecond->At(i))->Eta();
    secondPhi[i] = ((AliVParticle*) particlesSecond->At(i))->Phi();
    secondPt[i]  = ((AliVParticle*) particlesSecond->At(i))->Pt();
    secondCharge[i]  = (Short_t)((AliVParticle*) particlesSecond->At(i))->Charge();
    secondCorrection[i]  = (Double_t)((AliBFBasicParticle*) particlesSecond->At(i))->Correction();   
  }
  
  //TLorenzVector implementation for resonances
  TLorentzVector vectorMother, vectorDaughter[2];
  TParticle pPion, pProton, pRho0, pK0s, pLambda;
  pPion.SetPdgCode(211); //pion
  pRho0.SetPdgCode(113); //rho0
  pK0s.SetPdgCode(310); //K0s
  pProton.SetPdgCode(2212); //proton
  pLambda.SetPdgCode(3122); //Lambda
  Double_t gWidthForRho0 = 0.01;
  Double_t gWidthForK0s = 0.01;
  Double_t gWidthForLambda = 0.006;
  Double_t nSigmaRejection = 3.0;

  // count triggers per event
  Int_t triggersP = 0;
  Int_t triggersN = 0;

  // reset ebye histograms
  fHistPN->Reset();
  fHistNP->Reset();
  fHistPP->Reset();
  fHistNN->Reset();
  fHistBF->Reset();

  // 1st particle loop
  for (Int_t i = 0; i < iMax; i++) {
    AliBFBasicParticle* firstParticle = (AliBFBasicParticle*) particles->At(i); 
    
    // some optimization
    Float_t firstEta = firstParticle->Eta();
    Float_t firstPhi = firstParticle->Phi();
    Float_t firstPt  = firstParticle->Pt();
    Float_t firstCorrection  = firstParticle->Correction();

    // Event plane (determine psi bin)
    Double_t gPsiMinusPhi    =   0.;
    Double_t gPsiMinusPhiBin = -10.;
    gPsiMinusPhi   = TMath::Abs(firstPhi - gReactionPlane);
    //in-plane
    if((gPsiMinusPhi <= 7.5*TMath::DegToRad())||
       ((172.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 187.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 0.0;
    //intermediate
    else if(((37.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 52.5*TMath::DegToRad()))||
	    ((127.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 142.5*TMath::DegToRad()))||
	    ((217.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 232.5*TMath::DegToRad()))||
	    ((307.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 322.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 1.0;
    //out of plane
    else if(((82.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 97.5*TMath::DegToRad()))||
	    ((262.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 277.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 2.0;
    //everything else
    else 
      gPsiMinusPhiBin = 3.0;
    
    fHistPsiMinusPhi->Fill(gPsiMinusPhiBin,gPsiMinusPhi);

    Short_t  charge1 = (Short_t) firstParticle->Charge();

    if(charge1 > 0)
      triggersP ++;
    else if(charge1 < 0)
      triggersN ++;
    
    trackVariablesSingle[0] = kMultorCent;
    trackVariablesSingle[1] = firstPt;
    trackVariablesSingle[2] = vertexZ;
    
    // 2nd particle loop
    for(Int_t j = 0; j < jMax; j++) {   

      if(!particlesMixed && j == i) continue; // no auto correlations (only for non mixing)

      // pT,Assoc < pT,Trig (not used for EbyE)
      // if(firstPt < secondPt[j]) continue;

      Short_t charge2 = secondCharge[j];
      
      trackVariablesPair[0]    =  trackVariablesSingle[0];
      trackVariablesPair[1]    =  firstEta - secondEta[j];  // delta eta
      trackVariablesPair[2]    =  firstPhi - secondPhi[j];  // delta phi
      //if (trackVariablesPair[2] > 180.)   // delta phi between -180 and 180 
      //trackVariablesPair[2] -= 360.;
      //if (trackVariablesPair[2] <  - 180.) 
      //trackVariablesPair[2] += 360.;
      if (trackVariablesPair[2] > TMath::Pi()) // delta phi between -pi and pi 
	trackVariablesPair[2] -= 2.*TMath::Pi();
      if (trackVariablesPair[2] <  - TMath::Pi()) 
	trackVariablesPair[2] += 2.*TMath::Pi();
      if (trackVariablesPair[2] <  - TMath::Pi()/2.) 
      trackVariablesPair[2] += 2.*TMath::Pi();
      
      trackVariablesPair[3]    =  firstPt;      // pt trigger
      trackVariablesPair[4]    =  secondPt[j];  // pt
      trackVariablesPair[5]    =  vertexZ;      // z of the primary vertex
      
      //Exclude resonances for the calculation of pairs by looking 
      //at the invariant mass and not considering the pairs that 
      //fall within 3sigma from the mass peak of: rho0, K0s, Lambda
      if(fResonancesCut) {
	if (charge1 * charge2 < 0) {

	  //rho0
	  vectorDaughter[0].SetPtEtaPhiM(firstPt,firstEta,firstPhi,pPion.GetMass());
	  vectorDaughter[1].SetPtEtaPhiM(secondPt[j],secondEta[j],secondPhi[j],pPion.GetMass());
	  vectorMother = vectorDaughter[0] + vectorDaughter[1];
	  fHistResonancesBefore->Fill(trackVariablesPair[1],trackVariablesPair[2],vectorMother.M());
	  if(TMath::Abs(vectorMother.M() - pRho0.GetMass()) <= nSigmaRejection*gWidthForRho0)
	    continue;
	  fHistResonancesRho->Fill(trackVariablesPair[1],trackVariablesPair[2],vectorMother.M());
	  
	  //K0s
	  if(TMath::Abs(vectorMother.M() - pK0s.GetMass()) <= nSigmaRejection*gWidthForK0s)
	    continue;
	  fHistResonancesK0->Fill(trackVariablesPair[1],trackVariablesPair[2],vectorMother.M());
	  
	  
	  //Lambda
	  vectorDaughter[0].SetPtEtaPhiM(firstPt,firstEta,firstPhi,pPion.GetMass());
	  vectorDaughter[1].SetPtEtaPhiM(secondPt[j],secondEta[j],secondPhi[j],pProton.GetMass());
	  vectorMother = vectorDaughter[0] + vectorDaughter[1];
	  if(TMath::Abs(vectorMother.M() - pLambda.GetMass()) <= nSigmaRejection*gWidthForLambda)
	    continue;
	  
	  vectorDaughter[0].SetPtEtaPhiM(firstPt,firstEta,firstPhi,pProton.GetMass());
	  vectorDaughter[1].SetPtEtaPhiM(secondPt[j],secondEta[j],secondPhi[j],pPion.GetMass());
	  vectorMother = vectorDaughter[0] + vectorDaughter[1];
	  if(TMath::Abs(vectorMother.M() - pLambda.GetMass()) <= nSigmaRejection*gWidthForLambda)
	    continue;
	  fHistResonancesLambda->Fill(trackVariablesPair[1],trackVariablesPair[2],vectorMother.M());
	
	}//unlike-sign only
      }//resonance cut

      // HBT like cut
      //if(fHBTCut){ // VERSION 3 (all pairs)
      if(fHBTCut && charge1 * charge2 > 0){  // VERSION 2 (only for LS)
	//if( dphi < 3 || deta < 0.01 ){   // VERSION 1
	//  continue;
	
	Double_t deta = firstEta - secondEta[j];
	Double_t dphi = firstPhi - secondPhi[j];
	// VERSION 2 (Taken from DPhiCorrelations)
	// the variables & cuthave been developed by the HBT group 
	// see e.g. https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700
	fHistHBTbefore->Fill(deta,dphi);
	
	// optimization
	if (TMath::Abs(deta) < fHBTCutValue * 2.5 * 3) //fHBTCutValue = 0.02 [default for dphicorrelations]
	  {
	    // phi in rad

	    Float_t phi1rad = firstPhi;
	    Float_t phi2rad = secondPhi[j];
	    
	    // check first boundaries to see if is worth to loop and find the minimum
	    Float_t dphistar1 = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, 0.8, bSign);
	    Float_t dphistar2 = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, 2.5, bSign);
	    
	    const Float_t kLimit = fHBTCutValue * 3;
	    
	    Float_t dphistarminabs = 1e5;
	    //Float_t dphistarmin = 1e5;
	    
	    if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0 ) {
	      for (Double_t rad=0.8; rad<2.51; rad+=0.01) {
		Float_t dphistar = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, rad, bSign);
		Float_t dphistarabs = TMath::Abs(dphistar);
		
		if (dphistarabs < dphistarminabs) {
		  //dphistarmin = dphistar;
		  dphistarminabs = dphistarabs;
		}
	      }
	      
	      if (dphistarminabs < fHBTCutValue && TMath::Abs(deta) < fHBTCutValue) {
		continue;
	      }
	    }
	  }
	fHistHBTafter->Fill(deta,dphi);
      }//HBT cut
	
      // conversions
      if(fConversionCut) {
	if (charge1 * charge2 < 0) {
	  Double_t deta = firstEta - secondEta[j];
	  Double_t dphi = firstPhi - secondPhi[j];
	  
	  Float_t m0 = 0.510e-3;
	  Float_t tantheta1 = 1e10;
	  
	  // phi in rad
	  Float_t phi1rad = firstPhi;
	  Float_t phi2rad = secondPhi[j];
	  
	  if (firstEta < -1e-10 || firstEta > 1e-10)
	    tantheta1 = 2 * TMath::Exp(-firstEta) / ( 1 - TMath::Exp(-2*firstEta));
	  
	  Float_t tantheta2 = 1e10;
	  if (secondEta[j] < -1e-10 || secondEta[j] > 1e-10)
	    tantheta2 = 2 * TMath::Exp(-secondEta[j]) / ( 1 - TMath::Exp(-2*secondEta[j]));
	  
	  Float_t e1squ = m0 * m0 + firstPt * firstPt * (1.0 + 1.0 / tantheta1 / tantheta1);
	  Float_t e2squ = m0 * m0 + secondPt[j] * secondPt[j] * (1.0 + 1.0 / tantheta2 / tantheta2);
	  
	  Float_t masssqu = 2 * m0 * m0 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( firstPt * secondPt[j] * ( TMath::Cos(phi1rad - phi2rad) + 1.0 / tantheta1 / tantheta2 ) ) );

	  fHistConversionbefore->Fill(deta,dphi,masssqu);
	  
	  if (masssqu < fInvMassCutConversion*fInvMassCutConversion){
	    continue;
	  }
	  fHistConversionafter->Fill(deta,dphi,masssqu);
	}
      }//conversion cut

      // momentum difference cut - suppress femtoscopic effects
      if(fQCut){ 

	Double_t ptDifference = TMath::Abs( firstPt - secondPt[j]);

	fHistQbefore->Fill(trackVariablesPair[1],trackVariablesPair[2],ptDifference);
	if(ptDifference < fDeltaPtMin) continue;
	fHistQafter->Fill(trackVariablesPair[1],trackVariablesPair[2],ptDifference);

      }

      if( charge1 > 0 && charge2 < 0)  fHistPN->Fill(trackVariablesPair[1],trackVariablesPair[2]); 
      else if( charge1 < 0 && charge2 > 0)  fHistNP->Fill(trackVariablesPair[1],trackVariablesPair[2]); 
      else if( charge1 > 0 && charge2 > 0)  fHistPP->Fill(trackVariablesPair[1],trackVariablesPair[2]); 
      else if( charge1 < 0 && charge2 < 0)  fHistNN->Fill(trackVariablesPair[1],trackVariablesPair[2]); 
      else {
	continue;
      }

    }//end of 2nd particle loop
  }//end of 1st particle loop

  // build balance function event-by-event
  fHistPN->Scale(1./(Double_t)triggersP);
  fHistNP->Scale(1./(Double_t)triggersN);
  fHistPP->Scale(1./(Double_t)triggersP);
  fHistNN->Scale(1./(Double_t)triggersN);
  fHistBF->Add(fHistPN,fHistNP,1.,1.);
  fHistBF->Add(fHistPP,-1.);
  fHistBF->Add(fHistNN,-1.);
  fHistBF->Scale(0.5);

  // add to sum (2D)
  //fHistBFSum->Add(fHistBF);

  // add to sum (3D)
  for(Int_t iCent = 0; iCent < fHistBFSum->GetNbinsX(); iCent++){
    if(trackVariablesSingle[0] >= fHistBFSum->GetXaxis()->GetBinLowEdge(iCent+1) && trackVariablesSingle[0] < fHistBFSum->GetXaxis()->GetBinUpEdge(iCent+1)){
      for(Int_t iEta = 0; iEta < fHistBFSum->GetNbinsY(); iEta++){
	for(Int_t iPhi = 0; iPhi < fHistBFSum->GetNbinsZ(); iPhi++){
	  
	  Float_t oldContent   = fHistBFSum->GetBinContent(iCent+1,iEta+1,iPhi+1);
	  Float_t eventContent = fHistBF->GetBinContent(iEta+1,iPhi+1);
	  Float_t newContent   = oldContent + eventContent;

	  Float_t oldError     = fHistBFSum->GetBinError(iCent+1,iEta+1,iPhi+1);
	  Float_t eventError   = fHistBF->GetBinError(iEta+1,iPhi+1);
	  Float_t newError     = TMath::Sqrt(oldError*oldError + eventError*eventError);

	  fHistBFSum->SetBinContent(iCent+1,iEta+1,iPhi+1,newContent);
	  fHistBFSum->SetBinError(iCent+1,iEta+1,iPhi+1,newError);
	  
	}
      }
    }
  }
}  



//____________________________________________________________________//
Float_t AliBalanceEbyE::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign) { 
  //
  // calculates dphistar
  //
  Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  
  static const Double_t kPi = TMath::Pi();
  
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;
  
  return dphistar;
}
