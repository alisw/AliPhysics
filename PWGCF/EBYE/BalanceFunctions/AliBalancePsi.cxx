/**************************************************************************
 * Author: Panos Christakoglou.                                           *
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

/* $Id: AliBalancePsi.cxx 54125 2012-01-24 21:07:41Z miweber $ */

//-----------------------------------------------------------------
//           Balance Function class
//   This is the class to deal with the Balance Function wrt Psi analysis
//   Origin: Panos Christakoglou, Nikhef, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------


//ROOT
#include <Riostream.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TString.h>

#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliTHn.h"
#include "AliAnalysisTaskTriggeredBF.h"

#include "AliBalancePsi.h"

ClassImp(AliBalancePsi)

//____________________________________________________________________//
AliBalancePsi::AliBalancePsi() :
  TObject(), 
  fShuffle(kFALSE),
  fAnalysisLevel("ESD"),
  fAnalyzedEvents(0) ,
  fCentralityId(0) ,
  fCentStart(0.),
  fCentStop(0.),
  fHistP(0),
  fHistN(0),
  fHistPN(0),
  fHistNP(0),
  fHistPP(0),
  fHistNN(0),
  fHistHBTbefore(0),
  fHistHBTafter(0),
  fHistConversionbefore(0),
  fHistConversionafter(0),
  fHistPsiMinusPhi(0),
  fPsiInterval(15.),
  fDeltaEtaMax(2.0),
  fHBTCut(kFALSE),
  fConversionCut(kFALSE),
  fEventClass("EventPlane"){
  // Default constructor
}

//____________________________________________________________________//
AliBalancePsi::AliBalancePsi(const AliBalancePsi& balance):
  TObject(balance), fShuffle(balance.fShuffle), 
  fAnalysisLevel(balance.fAnalysisLevel),
  fAnalyzedEvents(balance.fAnalyzedEvents), 
  fCentralityId(balance.fCentralityId),
  fCentStart(balance.fCentStart),
  fCentStop(balance.fCentStop),
  fHistP(balance.fHistP),
  fHistN(balance.fHistN),
  fHistPN(balance.fHistPN),
  fHistNP(balance.fHistNP),
  fHistPP(balance.fHistPP),
  fHistNN(balance.fHistNN),
  fHistHBTbefore(balance.fHistHBTbefore),
  fHistHBTafter(balance.fHistHBTafter),
  fHistConversionbefore(balance.fHistConversionbefore),
  fHistConversionafter(balance.fHistConversionafter),
  fHistPsiMinusPhi(balance.fHistPsiMinusPhi),
  fPsiInterval(balance.fPsiInterval),
  fDeltaEtaMax(balance.fDeltaEtaMax),
  fHBTCut(balance.fHBTCut),
  fConversionCut(balance.fConversionCut),
  fEventClass("EventPlane"){
  //copy constructor
}

//____________________________________________________________________//
AliBalancePsi::~AliBalancePsi() {
  // Destructor
  delete fHistP;
  delete fHistN;
  delete fHistPN;
  delete fHistNP;
  delete fHistPP;
  delete fHistNN;

  delete fHistHBTbefore;
  delete fHistHBTafter;
  delete fHistConversionbefore;
  delete fHistConversionafter;
  delete fHistPsiMinusPhi;
    
}

//____________________________________________________________________//
void AliBalancePsi::InitHistograms() {
  // single particle histograms

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  Int_t anaSteps   = 1;       // analysis steps
  Int_t iBinSingle[kTrackVariablesSingle];        // binning for track variables
  Double_t* dBinsSingle[kTrackVariablesSingle];   // bins for track variables  
  TString axisTitleSingle[kTrackVariablesSingle]; // axis titles for track variables
  
  // two particle histograms
  Int_t iBinPair[kTrackVariablesPair];         // binning for track variables
  Double_t* dBinsPair[kTrackVariablesPair];    // bins for track variables  
  TString axisTitlePair[kTrackVariablesPair];  // axis titles for track variables
  /**********************************************************
   
  ======> Modification: Change Event Classification Scheme
    
  ---> fEventClass == "EventPlane"
   
   Default operation with Event Plane 
   
  ---> fEventClass == "Multiplicity"
   
   Work with kTPCITStracklet multiplicity (from GetReferenceMultiplicity)
   
  ---> fEventClass == "Centrality" 
   
   Work with Centrality Bins

  ***********************************************************/
   
  //--- Multiplicity Bins ------------------------------------
    const Int_t kMultBins = 8;
    //A first rough attempt at four bins
    Double_t kMultBinLimits[kMultBins+1]={0,10,20,30,40,50,60,70,80};
  //----------------------------------------------------------
    
  //--- Centrality Bins --------------------------------------
    const Int_t kNCentralityBins = 9;
    Double_t centralityBins[kNCentralityBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};
  //----------------------------------------------------------
    
  //--- Event Plane Bins -------------------------------------
    //Psi_2: -0.5->0.5 (in plane), 0.5->1.5 (intermediate), 1.5->2.5 (out of plane), 2.5->3.5 (rest)
    const Int_t kNPsi2Bins = 4;
    Double_t psi2Bins[kNPsi2Bins+1] = {-0.5,0.5,1.5,2.5,3.5};
  //----------------------------------------------------------
    
  //Depending on fEventClass Variable, do one thing or the other...
    if(fEventClass == "Multiplicity"){
        iBinSingle[0]       = kMultBins;
        dBinsSingle[0]      = kMultBinLimits;
        axisTitleSingle[0]  = "kTPCITStracklet multiplicity";
        iBinPair[0]       = kMultBins;
        dBinsPair[0]      = kMultBinLimits;
        axisTitlePair[0]  = "kTPCITStracklet multiplicity";
    }
    if(fEventClass == "Centrality"){
        iBinSingle[0]       = kNCentralityBins;
        dBinsSingle[0]      = centralityBins;
        axisTitleSingle[0]  = "Centrality percentile [%]";
        iBinPair[0]       = kNCentralityBins;
        dBinsPair[0]      = centralityBins;
        axisTitlePair[0]  = "Centrality percentile [%]";
    }
    if(fEventClass == "EventPlane"){
        iBinSingle[0]       = kNPsi2Bins;
        dBinsSingle[0]      = psi2Bins;
        axisTitleSingle[0]  = "#varphi - #Psi_{2} (a.u.)";
        iBinPair[0]       = kNPsi2Bins;
        dBinsPair[0]      = psi2Bins;
        axisTitlePair[0]  = "#varphi - #Psi_{2} (a.u.)";
    }
  
   // delta eta
  const Int_t kNDeltaEtaBins = 80;
  Double_t deltaEtaBins[kNDeltaEtaBins+1];
  for(Int_t i = 0; i < kNDeltaEtaBins+1; i++)
    deltaEtaBins[i] = - fDeltaEtaMax + i * 2 * fDeltaEtaMax / (Double_t)kNDeltaEtaBins;
  iBinPair[1]       = kNDeltaEtaBins;
  dBinsPair[1]      = deltaEtaBins;
  axisTitlePair[1]  = "#Delta#eta"; 

   // delta phi
  const Int_t kNDeltaPhiBins = 72;
  Double_t deltaPhiBins[kNDeltaPhiBins+1];
  for(Int_t i = 0; i < kNDeltaPhiBins+1; i++){
    //deltaPhiBins[i] = -180.0 + i * 5.;
    deltaPhiBins[i] = -TMath::Pi()/2. + i * 5.*TMath::Pi()/180.;
  } 
  iBinPair[2]       = kNDeltaPhiBins;
  dBinsPair[2]      = deltaPhiBins;
  axisTitlePair[2]  = "#Delta#varphi (rad)"; 

  // pt(trigger-associated)
  const Int_t kNPtBins = 16;
  Double_t ptBins[kNPtBins+1] = {0.2,0.6,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.,12.,15.,20.};
  //for(Int_t i = 0; i < kNPtBins+1; i++){
  //ptBins[i] = 0.2 + i * 0.5;
  //} 
  iBinSingle[1]       = kNPtBins;
  dBinsSingle[1]      = ptBins;
  axisTitleSingle[1]  = "p_{T,trig.} (GeV/c)"; 

  iBinPair[3]       = kNPtBins;
  dBinsPair[3]      = ptBins;
  axisTitlePair[3]  = "p_{T,trig.} (GeV/c)"; 

  iBinPair[4]       = kNPtBins;
  dBinsPair[4]      = ptBins;
  axisTitlePair[4]  = "p_{T,assoc.} (GeV/c)";   

  TString histName;
  //+ triggered particles
  histName = "fHistP"; 
  if(fShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistP = new AliTHn(histName.Data(),histName.Data(),anaSteps,kTrackVariablesSingle,iBinSingle);
  for (Int_t j=0; j<kTrackVariablesSingle; j++) {
    fHistP->SetBinLimits(j, dBinsSingle[j]);
    fHistP->SetVarTitle(j, axisTitleSingle[j]);
  }

  //- triggered particles
  histName = "fHistN"; 
  if(fShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistN = new AliTHn(histName.Data(),histName.Data(),anaSteps,kTrackVariablesSingle,iBinSingle);
  for (Int_t j=0; j<kTrackVariablesSingle; j++) {
    fHistN->SetBinLimits(j, dBinsSingle[j]);
    fHistN->SetVarTitle(j, axisTitleSingle[j]);
  }
  
  //+- pairs
  histName = "fHistPN";
  if(fShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistPN = new AliTHn(histName.Data(),histName.Data(),anaSteps, kTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fHistPN->SetBinLimits(j, dBinsPair[j]);
    fHistPN->SetVarTitle(j, axisTitlePair[j]);
  }

  //-+ pairs
  histName = "fHistNP";
  if(fShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistNP = new AliTHn(histName.Data(),histName.Data(),anaSteps, kTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fHistNP->SetBinLimits(j, dBinsPair[j]);
    fHistNP->SetVarTitle(j, axisTitlePair[j]);
  }

  //++ pairs
  histName = "fHistPP";
  if(fShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistPP = new AliTHn(histName.Data(),histName.Data(),anaSteps, kTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fHistPP->SetBinLimits(j, dBinsPair[j]);
    fHistPP->SetVarTitle(j, axisTitlePair[j]);
  }

  //-- pairs
  histName = "fHistNN";
  if(fShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistNN = new AliTHn(histName.Data(),histName.Data(),anaSteps, kTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fHistNN->SetBinLimits(j, dBinsPair[j]);
    fHistNN->SetVarTitle(j, axisTitlePair[j]);
  }
  AliInfo("Finished setting up the AliTHn");

  // QA histograms
  fHistHBTbefore        = new TH2D("fHistHBTbefore","before HBT cut",200,0,2,200,0,2.*TMath::Pi());
  fHistHBTafter         = new TH2D("fHistHBTafter","after HBT cut",200,0,2,200,0,2.*TMath::Pi());
  fHistConversionbefore = new TH2D("fHistConversionbefore","before Conversion cut",200,0,2,200,0,2.*TMath::Pi());
  fHistConversionafter  = new TH2D("fHistConversionafter","after Conversion cut",200,0,2,200,0,2.*TMath::Pi());
  fHistPsiMinusPhi     = new TH2D("fHistPsiMinusPhi","",4,-0.5,3.5,100,0,2.*TMath::Pi());

  TH1::AddDirectory(oldStatus);

}

//____________________________________________________________________//
void AliBalancePsi::CalculateBalance(Double_t gReactionPlane,
				     TObjArray *particles, 
				     TObjArray *particlesMixed,
                     Float_t bSign,
                     Double_t kMultorCent) {
  // Calculates the balance function
  fAnalyzedEvents++;
    
  // Initialize histograms if not done yet
  if(!fHistPN){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("This works only in local mode --> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    InitHistograms();
  }

  Double_t trackVariablesSingle[kTrackVariablesSingle];
  Double_t trackVariablesPair[kTrackVariablesPair];

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
    secondCorrection[i]  = (Double_t)((AliBFBasicParticle*) particlesSecond->At(i))->Correction();   //==========================correction
  }
  
  // 1st particle loop
  for (Int_t i = 0; i < iMax; i++) {
    //AliVParticle* firstParticle = (AliVParticle*) particles->At(i);
    AliBFBasicParticle* firstParticle = (AliBFBasicParticle*) particles->At(i); //==========================correction
    
    // some optimization
    Float_t firstEta = firstParticle->Eta();
    Float_t firstPhi = firstParticle->Phi();
    Float_t firstPt  = firstParticle->Pt();
    Float_t firstCorrection  = firstParticle->Correction();//==========================correction

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
    
    trackVariablesSingle[0]    =  gPsiMinusPhiBin;
    trackVariablesSingle[1]    =  firstPt;
      if(fEventClass=="Multiplicity" || fEventClass == "Centrality" ) trackVariablesSingle[0] = kMultorCent;
    
    //fill single particle histograms
    if(charge1 > 0)      fHistP->Fill(trackVariablesSingle,0,firstCorrection); //==========================correction
    else if(charge1 < 0) fHistN->Fill(trackVariablesSingle,0,firstCorrection);  //==========================correction
    
    // 2nd particle loop
    for(Int_t j = 0; j < jMax; j++) {   

      if(!particlesMixed && j == i) continue; // no auto correlations (only for non mixing)

      // pT,Assoc < pT,Trig
      if(firstPt < secondPt[j]) continue;

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
      //	trackVariablesPair[5]    =  fCentrality;  // centrality

      // HBT like cut
      if(fHBTCut){ // VERSION 3 (all pairs)
        //if(fHBTCut && charge1 * charge2 > 0){  // VERSION 2 (only for LS)
	//if( dphi < 3 || deta < 0.01 ){   // VERSION 1
	//  continue;
	
	Double_t deta = firstEta - secondEta[j];
	Double_t dphi = firstPhi - secondPhi[j];
	// VERSION 2 (Taken from DPhiCorrelations)
	// the variables & cuthave been developed by the HBT group 
	// see e.g. https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700
	fHistHBTbefore->Fill(deta,dphi);
	
	// optimization
	if (TMath::Abs(deta) < 0.02 * 2.5 * 3) //twoTrackEfficiencyCutValue = 0.02 [default for dphicorrelations]
	  {
	    // phi in rad
	    //Float_t phi1rad = firstPhi*TMath::DegToRad();
	    //Float_t phi2rad = secondPhi[j]*TMath::DegToRad();
	    Float_t phi1rad = firstPhi;
	    Float_t phi2rad = secondPhi[j];
	    
	    // check first boundaries to see if is worth to loop and find the minimum
	    Float_t dphistar1 = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, 0.8, bSign);
	    Float_t dphistar2 = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, 2.5, bSign);
	    
	    const Float_t kLimit = 0.02 * 3;
	    
	    Float_t dphistarminabs = 1e5;
	    Float_t dphistarmin = 1e5;
	    
	    if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0 ) {
	      for (Double_t rad=0.8; rad<2.51; rad+=0.01) {
		Float_t dphistar = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, rad, bSign);
		Float_t dphistarabs = TMath::Abs(dphistar);
		
		if (dphistarabs < dphistarminabs) {
		  dphistarmin = dphistar;
		  dphistarminabs = dphistarabs;
		}
	      }
	      
	      if (dphistarminabs < 0.02 && TMath::Abs(deta) < 0.02) {
		//AliInfo(Form("HBT: Removed track pair %d %d with [[%f %f]] %f %f %f | %f %f %d %f %f %d %f", i, j, deta, dphi, dphistarminabs, dphistar1, dphistar2, phi1rad, pt1, charge1, phi2rad, pt2, charge2, bSign));
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
	  fHistConversionbefore->Fill(deta,dphi);
	  
	  Float_t m0 = 0.510e-3;
	  Float_t tantheta1 = 1e10;
	  
	  // phi in rad
	  //Float_t phi1rad = firstPhi*TMath::DegToRad();
	  //Float_t phi2rad = secondPhi[j]*TMath::DegToRad();
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
	  
	  if (masssqu < 0.04*0.04){
	    //AliInfo(Form("Conversion: Removed track pair %d %d with [[%f %f] %f %f] %d %d <- %f %f  %f %f   %f %f ", i, j, deta, dphi, masssqu, charge1, charge2,eta1,eta2,phi1,phi2,pt1,pt2));
	    continue;
	  }
	  fHistConversionafter->Fill(deta,dphi);
	}
      }//conversion cut
      
      if( charge1 > 0 && charge2 < 0)  fHistPN->Fill(trackVariablesPair,0,firstCorrection*secondCorrection[j]); //==========================correction
      else if( charge1 < 0 && charge2 > 0)  fHistNP->Fill(trackVariablesPair,0,firstCorrection*secondCorrection[j]);//==========================correction 
      else if( charge1 > 0 && charge2 > 0)  fHistPP->Fill(trackVariablesPair,0,firstCorrection*secondCorrection[j]);//==========================correction 
      else if( charge1 < 0 && charge2 < 0)  fHistNN->Fill(trackVariablesPair,0,firstCorrection*secondCorrection[j]);//==========================correction 
      else {
	//AliWarning(Form("Wrong charge combination: charge1 = %d and charge2 = %d",charge,charge2));
	continue;
      }
    }//end of 2nd particle loop
  }//end of 1st particle loop
}  

//____________________________________________________________________//
TH1D *AliBalancePsi::GetBalanceFunctionHistogram(Int_t iVariableSingle,
						 Int_t iVariablePair,
						 Double_t psiMin, 
						 Double_t psiMax,
						 Double_t ptTriggerMin,
						 Double_t ptTriggerMax,
						 Double_t ptAssociatedMin,
						 Double_t ptAssociatedMax) {
  //Returns the BF histogram, extracted from the 6 AliTHn objects 
  //(private members) of the AliBalancePsi class.
  //iVariableSingle: 0(phi-Psi), 1(pt-trigger)
  //iVariablePair: 0(phi-Psi) 1(Delta eta), 2(Delta phi), 3(pt-trigger), 4(pt-associated

  // security checks
  if(psiMin > psiMax-0.00001)                   psiMin          = psiMax-0.00001;
  if(ptTriggerMin > ptTriggerMax-0.00001)       ptTriggerMin    = ptTriggerMax-0.00001;
  if(ptAssociatedMin > ptAssociatedMax-0.00001) ptAssociatedMin = ptAssociatedMax-0.00001;

  // Psi_2
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistP->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistN->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.)) {
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
  }

  //Printf("P:%lf - N:%lf - PN:%lf - NP:%lf - PP:%lf - NN:%lf",fHistP->GetEntries(0),fHistN->GetEntries(0),fHistPN->GetEntries(0),fHistNP->GetEntries(0),fHistPP->GetEntries(0),fHistNN->GetEntries(0));

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH1D* hTemp1 = (TH1D*)fHistPN->Project(0,iVariablePair);
  TH1D* hTemp2 = (TH1D*)fHistNP->Project(0,iVariablePair);
  TH1D* hTemp3 = (TH1D*)fHistPP->Project(0,iVariablePair);
  TH1D* hTemp4 = (TH1D*)fHistNN->Project(0,iVariablePair);
  TH1D* hTemp5 = (TH1D*)fHistP->Project(0,iVariableSingle);
  TH1D* hTemp6 = (TH1D*)fHistN->Project(0,iVariableSingle);

  TH1D *gHistBalanceFunctionHistogram = 0x0;
  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)&&(hTemp5)&&(hTemp6)) {
    gHistBalanceFunctionHistogram = (TH1D*)hTemp1->Clone();
    gHistBalanceFunctionHistogram->Reset();
    
    switch(iVariablePair) {
    case 1:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta#eta");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta#eta)");
      break;
    case 2:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta#varphi (rad)");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta#varphi)");
      break;
    default:
      break;
    }

    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1->Add(hTemp3,-1.);
    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp2->Add(hTemp4,-1.);
    hTemp2->Scale(1./hTemp6->GetEntries());
    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5);

    //normalize to bin width
    gHistBalanceFunctionHistogram->Scale(1./((Double_t)gHistBalanceFunctionHistogram->GetXaxis()->GetBinWidth(1)*(Double_t)gHistBalanceFunctionHistogram->GetYaxis()->GetBinWidth(1)));
  }

  return gHistBalanceFunctionHistogram;
}

//____________________________________________________________________//
TH1D *AliBalancePsi::GetBalanceFunctionHistogram2pMethod(Int_t iVariableSingle,
							 Int_t iVariablePair,
							 Double_t psiMin, 
							 Double_t psiMax,
							 Double_t ptTriggerMin,
							 Double_t ptTriggerMax,
							 Double_t ptAssociatedMin,
							 Double_t ptAssociatedMax,
							 AliBalancePsi *bfMix) {
  //Returns the BF histogram, extracted from the 6 AliTHn objects 
  //after dividing each correlation function by the Event Mixing one 
  //(private members) of the AliBalancePsi class.
  //iVariableSingle: 0(phi-Psi), 1(pt-trigger)
  //iVariablePair: 0(phi-Psi) 1(Delta eta), 2(Delta phi), 3(pt-trigger), 4(pt-associated

  // security checks
  if(psiMin > psiMax-0.00001)                   psiMin          = psiMax-0.00001;
  if(ptTriggerMin > ptTriggerMax-0.00001)       ptTriggerMin    = ptTriggerMax-0.00001;
  if(ptAssociatedMin > ptAssociatedMax-0.00001) ptAssociatedMin = ptAssociatedMax-0.00001;
  
  // Psi_2
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistP->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistN->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.)) {
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
  }

  // ============================================================================================
  // the same for event mixing
  AliTHn *fHistPMix = bfMix->GetHistNp();
  AliTHn *fHistNMix = bfMix->GetHistNn();
  AliTHn *fHistPNMix = bfMix->GetHistNpn();
  AliTHn *fHistNPMix = bfMix->GetHistNnp();
  AliTHn *fHistPPMix = bfMix->GetHistNpp();
  AliTHn *fHistNNMix = bfMix->GetHistNnn();

  // Psi_2
  fHistPMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPNMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNPMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPPMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNNMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistPMix->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNMix->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPNMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNPMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPPMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNNMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.)) {
    fHistPNMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNPMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistPPMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNNMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
  }
  // ============================================================================================

  //Printf("P:%lf - N:%lf - PN:%lf - NP:%lf - PP:%lf - NN:%lf",fHistP->GetEntries(0),fHistN->GetEntries(0),fHistPN->GetEntries(0),fHistNP->GetEntries(0),fHistPP->GetEntries(0),fHistNN->GetEntries(0));

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH1D* hTemp1 = (TH1D*)fHistPN->Project(0,iVariablePair);
  TH1D* hTemp2 = (TH1D*)fHistNP->Project(0,iVariablePair);
  TH1D* hTemp3 = (TH1D*)fHistPP->Project(0,iVariablePair);
  TH1D* hTemp4 = (TH1D*)fHistNN->Project(0,iVariablePair);
  TH1D* hTemp5 = (TH1D*)fHistP->Project(0,iVariableSingle);
  TH1D* hTemp6 = (TH1D*)fHistN->Project(0,iVariableSingle);

  // ============================================================================================
  // the same for event mixing
  TH2D* hTemp1Mix = (TH2D*)fHistPNMix->Project(0,iVariablePair);
  TH2D* hTemp2Mix = (TH2D*)fHistNPMix->Project(0,iVariablePair);
  TH2D* hTemp3Mix = (TH2D*)fHistPPMix->Project(0,iVariablePair);
  TH2D* hTemp4Mix = (TH2D*)fHistNNMix->Project(0,iVariablePair);
  TH1D* hTemp5Mix = (TH1D*)fHistPMix->Project(0,iVariableSingle);
  TH1D* hTemp6Mix = (TH1D*)fHistNMix->Project(0,iVariableSingle);
  // ============================================================================================

  TH1D *gHistBalanceFunctionHistogram = 0x0;
  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)&&(hTemp5)&&(hTemp6)) {
    gHistBalanceFunctionHistogram = (TH1D*)hTemp1->Clone();
    gHistBalanceFunctionHistogram->Reset();
    
    switch(iVariablePair) {
    case 1:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta#eta");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta#eta)");
      break;
    case 2:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta#varphi (rad)");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta#varphi)");
      break;
    default:
      break;
    }

    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1Mix->Sumw2();
    hTemp2Mix->Sumw2();
    hTemp3Mix->Sumw2();
    hTemp4Mix->Sumw2();

    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp3->Scale(1./hTemp5->GetEntries());
    hTemp2->Scale(1./hTemp6->GetEntries());
    hTemp4->Scale(1./hTemp6->GetEntries());
    hTemp1Mix->Scale(1./hTemp5Mix->GetEntries());
    hTemp3Mix->Scale(1./hTemp5Mix->GetEntries());
    hTemp2Mix->Scale(1./hTemp6Mix->GetEntries());
    hTemp4Mix->Scale(1./hTemp6Mix->GetEntries());

    hTemp1->Divide(hTemp1Mix);
    hTemp2->Divide(hTemp2Mix);
    hTemp3->Divide(hTemp3Mix);
    hTemp4->Divide(hTemp4Mix);

    hTemp1->Add(hTemp3,-1.);
    hTemp2->Add(hTemp4,-1.);

    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5);

    //normalize to bin width
    gHistBalanceFunctionHistogram->Scale(1./((Double_t)gHistBalanceFunctionHistogram->GetXaxis()->GetBinWidth(1)*(Double_t)gHistBalanceFunctionHistogram->GetYaxis()->GetBinWidth(1)));
  }

  return gHistBalanceFunctionHistogram;
}

//____________________________________________________________________//
TH2D *AliBalancePsi::GetBalanceFunctionDeltaEtaDeltaPhi(Double_t psiMin, 
							Double_t psiMax,
							Double_t ptTriggerMin,
							Double_t ptTriggerMax,
							Double_t ptAssociatedMin,
							Double_t ptAssociatedMax) {
  //Returns the BF histogram in Delta eta vs Delta phi, 
  //extracted from the 6 AliTHn objects 
  //(private members) of the AliBalancePsi class.
  //iVariableSingle: 0(phi-Psi), 1(pt-trigger)
  //iVariablePair: 0(phi-Psi) 1(Delta eta), 2(Delta phi), 3(pt-trigger), 4(pt-associated

  // security checks
  if(psiMin > psiMax-0.00001)                   psiMin          = psiMax-0.00001;
  if(ptTriggerMin > ptTriggerMax-0.00001)       ptTriggerMin    = ptTriggerMax-0.00001;
  if(ptAssociatedMin > ptAssociatedMax-0.00001) ptAssociatedMin = ptAssociatedMax-0.00001;

  TString histName = "gHistBalanceFunctionHistogram2D";

  // Psi_2
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistP->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistN->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.)) {
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
  }

  //AliInfo(Form("P:%lf - N:%lf - PN:%lf - NP:%lf - PP:%lf - NN:%lf",fHistP->GetEntries(0),fHistN->GetEntries(0),fHistPN->GetEntries(0),fHistNP->GetEntries(0),fHistPP->GetEntries(0),fHistNN->GetEntries(0)));

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH2D* hTemp1 = (TH2D*)fHistPN->Project(0,1,2);
  TH2D* hTemp2 = (TH2D*)fHistNP->Project(0,1,2);
  TH2D* hTemp3 = (TH2D*)fHistPP->Project(0,1,2);
  TH2D* hTemp4 = (TH2D*)fHistNN->Project(0,1,2);
  TH1D* hTemp5 = (TH1D*)fHistP->Project(0,1);
  TH1D* hTemp6 = (TH1D*)fHistN->Project(0,1);

  TH2D *gHistBalanceFunctionHistogram = 0x0;
  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)&&(hTemp5)&&(hTemp6)) {
    gHistBalanceFunctionHistogram = (TH2D*)hTemp1->Clone();
    gHistBalanceFunctionHistogram->Reset();
    gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta#eta");   
    gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("#Delta#varphi (rad)");
    gHistBalanceFunctionHistogram->GetZaxis()->SetTitle("B(#Delta#eta,#Delta#varphi)");   
    
    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1->Add(hTemp3,-1.);
    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp2->Add(hTemp4,-1.);
    hTemp2->Scale(1./hTemp6->GetEntries());
    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5);

    //normalize to bin width
    gHistBalanceFunctionHistogram->Scale(1./((Double_t)gHistBalanceFunctionHistogram->GetXaxis()->GetBinWidth(1)*(Double_t)gHistBalanceFunctionHistogram->GetYaxis()->GetBinWidth(1)));
  }

  return gHistBalanceFunctionHistogram;
}

//____________________________________________________________________//
TH2D *AliBalancePsi::GetBalanceFunctionDeltaEtaDeltaPhi2pMethod(Double_t psiMin, 
								Double_t psiMax,
								Double_t ptTriggerMin,
								Double_t ptTriggerMax,
								Double_t ptAssociatedMin,
								Double_t ptAssociatedMax,
								AliBalancePsi *bfMix) {
  //Returns the BF histogram in Delta eta vs Delta phi,
  //after dividing each correlation function by the Event Mixing one 
  //extracted from the 6 AliTHn objects 
  //(private members) of the AliBalancePsi class.
  //iVariableSingle: 0(phi-Psi), 1(pt-trigger)
  //iVariablePair: 0(phi-Psi) 1(Delta eta), 2(Delta phi), 3(pt-trigger), 4(pt-associated
  TString histName = "gHistBalanceFunctionHistogram2D";

  if(!bfMix){
    AliError("balance function object for event mixing not available");
    return NULL;
  }

  // security checks
  if(psiMin > psiMax-0.00001)                   psiMin          = psiMax-0.00001;
  if(ptTriggerMin > ptTriggerMax-0.00001)       ptTriggerMin    = ptTriggerMax-0.00001;
  if(ptAssociatedMin > ptAssociatedMax-0.00001) ptAssociatedMin = ptAssociatedMax-0.00001;

  // Psi_2
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistP->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistN->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.)) {
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
  }


  // ============================================================================================
  // the same for event mixing
  AliTHn *fHistPMix = bfMix->GetHistNp();
  AliTHn *fHistNMix = bfMix->GetHistNn();
  AliTHn *fHistPNMix = bfMix->GetHistNpn();
  AliTHn *fHistNPMix = bfMix->GetHistNnp();
  AliTHn *fHistPPMix = bfMix->GetHistNpp();
  AliTHn *fHistNNMix = bfMix->GetHistNnn();

  // Psi_2
  fHistPMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPNMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNPMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPPMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNNMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistPMix->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNMix->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPNMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNPMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPPMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNNMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.)) {
    fHistPNMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNPMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistPPMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNNMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
  }
  // ============================================================================================


  //AliInfo(Form("P:%lf - N:%lf - PN:%lf - NP:%lf - PP:%lf - NN:%lf",fHistP->GetEntries(0),fHistN->GetEntries(0),fHistPN->GetEntries(0),fHistNP->GetEntries(0),fHistPP->GetEntries(0),fHistNN->GetEntries(0)));

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH2D* hTemp1 = (TH2D*)fHistPN->Project(0,1,2);
  TH2D* hTemp2 = (TH2D*)fHistNP->Project(0,1,2);
  TH2D* hTemp3 = (TH2D*)fHistPP->Project(0,1,2);
  TH2D* hTemp4 = (TH2D*)fHistNN->Project(0,1,2);
  TH1D* hTemp5 = (TH1D*)fHistP->Project(0,1);
  TH1D* hTemp6 = (TH1D*)fHistN->Project(0,1);

  // ============================================================================================
  // the same for event mixing
  TH2D* hTemp1Mix = (TH2D*)fHistPNMix->Project(0,1,2);
  TH2D* hTemp2Mix = (TH2D*)fHistNPMix->Project(0,1,2);
  TH2D* hTemp3Mix = (TH2D*)fHistPPMix->Project(0,1,2);
  TH2D* hTemp4Mix = (TH2D*)fHistNNMix->Project(0,1,2);
  TH1D* hTemp5Mix = (TH1D*)fHistPMix->Project(0,1);
  TH1D* hTemp6Mix = (TH1D*)fHistNMix->Project(0,1);
  // ============================================================================================

  TH2D *gHistBalanceFunctionHistogram = 0x0;
  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)&&(hTemp5)&&(hTemp6)) {
    gHistBalanceFunctionHistogram = (TH2D*)hTemp1->Clone();
    gHistBalanceFunctionHistogram->Reset();
    gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta#eta");   
    gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("#Delta#varphi (rad)");
    gHistBalanceFunctionHistogram->GetZaxis()->SetTitle("B(#Delta#eta,#Delta#varphi)");   
    
    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1Mix->Sumw2();
    hTemp2Mix->Sumw2();
    hTemp3Mix->Sumw2();
    hTemp4Mix->Sumw2();

    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp3->Scale(1./hTemp5->GetEntries());
    hTemp2->Scale(1./hTemp6->GetEntries());
    hTemp4->Scale(1./hTemp6->GetEntries());
    hTemp1Mix->Scale(1./hTemp5Mix->GetEntries());
    hTemp3Mix->Scale(1./hTemp5Mix->GetEntries());
    hTemp2Mix->Scale(1./hTemp6Mix->GetEntries());
    hTemp4Mix->Scale(1./hTemp6Mix->GetEntries());

    hTemp1->Divide(hTemp1Mix);
    hTemp2->Divide(hTemp2Mix);
    hTemp3->Divide(hTemp3Mix);
    hTemp4->Divide(hTemp4Mix);

    hTemp1->Add(hTemp3,-1.);
    hTemp2->Add(hTemp4,-1.);

    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5);
  
    //normalize to bin width
    gHistBalanceFunctionHistogram->Scale(1./((Double_t)gHistBalanceFunctionHistogram->GetXaxis()->GetBinWidth(1)*(Double_t)gHistBalanceFunctionHistogram->GetYaxis()->GetBinWidth(1)));
  }

  return gHistBalanceFunctionHistogram;
}

//____________________________________________________________________//
TH1D *AliBalancePsi::GetBalanceFunction1DFrom2D2pMethod(Bool_t bPhi,
                                                        Double_t psiMin,
                                                        Double_t psiMax,
                                                        Double_t ptTriggerMin,
                                                        Double_t ptTriggerMax,
                                                        Double_t ptAssociatedMin,
                                                        Double_t ptAssociatedMax,
                                                        AliBalancePsi *bfMix) {
  //Returns the BF histogram in Delta eta OR Delta phi,
  //after dividing each correlation function by the Event Mixing one
  // (But the division is done here in 2D, this was basically done to check the Integral)
  //extracted from the 6 AliTHn objects
  //(private members) of the AliBalancePsi class.
  //iVariableSingle: 0(phi-Psi), 1(pt-trigger)
  //iVariablePair: 0(phi-Psi) 1(Delta eta), 2(Delta phi), 3(pt-trigger), 4(pt-associated
  TString histName = "gHistBalanceFunctionHistogram1D";

  if(!bfMix){
    AliError("balance function object for event mixing not available");
    return NULL;
  }

  // security checks
  if(psiMin > psiMax-0.00001)                   psiMin          = psiMax-0.00001;
  if(ptTriggerMin > ptTriggerMax-0.00001)       ptTriggerMin    = ptTriggerMax-0.00001;
  if(ptAssociatedMin > ptAssociatedMax-0.00001) ptAssociatedMin = ptAssociatedMax-0.00001;

  // Psi_2
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistP->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistN->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.)) {
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
  }


  // ============================================================================================
  // the same for event mixing
  AliTHn *fHistPMix = bfMix->GetHistNp();
  AliTHn *fHistNMix = bfMix->GetHistNn();
  AliTHn *fHistPNMix = bfMix->GetHistNpn();
  AliTHn *fHistNPMix = bfMix->GetHistNnp();
  AliTHn *fHistPPMix = bfMix->GetHistNpp();
  AliTHn *fHistNNMix = bfMix->GetHistNnn();

  // Psi_2
  fHistPMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);
  fHistNMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);
  fHistPNMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);
  fHistNPMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);
  fHistPPMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);
  fHistNNMix->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001);

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistPMix->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNMix->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPNMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNPMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPPMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNNMix->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.)) {
    fHistPNMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNPMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistPPMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNNMix->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
  }
  // ============================================================================================


  //AliInfo(Form("P:%lf - N:%lf - PN:%lf - NP:%lf - PP:%lf - NN:%lf",fHistP->GetEntries(0),fHistN->GetEntries(0),fHistPN->GetEntries(0),fHistNP->GetEntries(0),fHistPP->GetEntries(0),fHistNN->GetEntries(0)));

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH2D* hTemp1 = (TH2D*)fHistPN->Project(0,1,2);
  TH2D* hTemp2 = (TH2D*)fHistNP->Project(0,1,2);
  TH2D* hTemp3 = (TH2D*)fHistPP->Project(0,1,2);
  TH2D* hTemp4 = (TH2D*)fHistNN->Project(0,1,2);
  TH1D* hTemp5 = (TH1D*)fHistP->Project(0,1);
  TH1D* hTemp6 = (TH1D*)fHistN->Project(0,1);

  // ============================================================================================
  // the same for event mixing
  TH2D* hTemp1Mix = (TH2D*)fHistPNMix->Project(0,1,2);
  TH2D* hTemp2Mix = (TH2D*)fHistNPMix->Project(0,1,2);
  TH2D* hTemp3Mix = (TH2D*)fHistPPMix->Project(0,1,2);
  TH2D* hTemp4Mix = (TH2D*)fHistNNMix->Project(0,1,2);
  TH1D* hTemp5Mix = (TH1D*)fHistPMix->Project(0,1);
  TH1D* hTemp6Mix = (TH1D*)fHistNMix->Project(0,1);
  // ============================================================================================

  TH1D *gHistBalanceFunctionHistogram = 0x0;
  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)&&(hTemp5)&&(hTemp6)) {

    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1Mix->Sumw2();
    hTemp2Mix->Sumw2();
    hTemp3Mix->Sumw2();
    hTemp4Mix->Sumw2();

    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp3->Scale(1./hTemp5->GetEntries());
    hTemp2->Scale(1./hTemp6->GetEntries());
    hTemp4->Scale(1./hTemp6->GetEntries());

    hTemp1Mix->Scale(1./hTemp5Mix->GetEntries());
    hTemp3Mix->Scale(1./hTemp5Mix->GetEntries());
    hTemp2Mix->Scale(1./hTemp6Mix->GetEntries());
    hTemp4Mix->Scale(1./hTemp6Mix->GetEntries());
  
    hTemp1->Divide(hTemp1Mix);
    hTemp2->Divide(hTemp2Mix);
    hTemp3->Divide(hTemp3Mix);
    hTemp4->Divide(hTemp4Mix);

    // now only project on one axis
    TH1D *h1DTemp1 = NULL;
    TH1D *h1DTemp2 = NULL;
    TH1D *h1DTemp3 = NULL;
    TH1D *h1DTemp4 = NULL;
    if(!bPhi){
      h1DTemp1 = (TH1D*)hTemp1->ProjectionX(Form("%s_projX",hTemp1->GetName()));
      h1DTemp2 = (TH1D*)hTemp2->ProjectionX(Form("%s_projX",hTemp2->GetName()));
      h1DTemp3 = (TH1D*)hTemp3->ProjectionX(Form("%s_projX",hTemp3->GetName()));
      h1DTemp4 = (TH1D*)hTemp4->ProjectionX(Form("%s_projX",hTemp4->GetName()));
    }
    else{
      h1DTemp1 = (TH1D*)hTemp1->ProjectionY(Form("%s_projX",hTemp1->GetName()));
      h1DTemp2 = (TH1D*)hTemp2->ProjectionY(Form("%s_projX",hTemp2->GetName()));
      h1DTemp3 = (TH1D*)hTemp3->ProjectionY(Form("%s_projX",hTemp3->GetName()));
      h1DTemp4 = (TH1D*)hTemp4->ProjectionY(Form("%s_projX",hTemp4->GetName()));
    }

    gHistBalanceFunctionHistogram = (TH1D*)h1DTemp1->Clone(Form("%s_clone",h1DTemp1->GetName()));
    gHistBalanceFunctionHistogram->Reset();
    if(!bPhi){
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta#eta");  
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta#eta)");  
    }
    else{
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta#varphi (rad)");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta#varphi)");  
    }

    h1DTemp1->Add(h1DTemp3,-1.);
    h1DTemp2->Add(h1DTemp4,-1.);

    gHistBalanceFunctionHistogram->Add(h1DTemp1,h1DTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5);
  }

  return gHistBalanceFunctionHistogram;
}


//____________________________________________________________________//
TH2D *AliBalancePsi::GetCorrelationFunctionPN(Double_t psiMin, 
					      Double_t psiMax,
					      Double_t ptTriggerMin,
					      Double_t ptTriggerMax,
					      Double_t ptAssociatedMin,
					      Double_t ptAssociatedMax) {
  //Returns the 2D correlation function for (+-) pairs

  // security checks
  if(psiMin > psiMax-0.00001)                   psiMin          = psiMax-0.00001;
  if(ptTriggerMin > ptTriggerMax-0.00001)       ptTriggerMin    = ptTriggerMax-0.00001;
  if(ptAssociatedMin > ptAssociatedMax-0.00001) ptAssociatedMin = ptAssociatedMax-0.00001;

  // Psi_2: axis 0
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistP->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.))
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);

  //fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(-0.5,2.5); 
  //fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(-0.5,2.5); 

  //TH2D *gHistTest = dynamic_cast<TH2D *>(fHistP->Project(0,0,1));
  //TCanvas *c1 = new TCanvas("c1","");
  //c1->cd();
  //if(!gHistTest){
  //AliError("Projection of fHistP = NULL");
  //return gHistTest;
  //}
  //else{
  //gHistTest->DrawCopy("colz");
  //}

  //0:step, 1: Delta eta, 2: Delta phi
  TH2D *gHist = dynamic_cast<TH2D *>(fHistPN->Project(0,1,2));
  if(!gHist){
    AliError("Projection of fHistPN = NULL");
    return gHist;
  }

  //AliInfo(Form("Entries (test): %lf",(Double_t)(gHistTest->GetEntries())));
  //AliInfo(Form("Entries (1D): %lf",(Double_t)(fHistP->Project(0,1)->GetEntries())));
  //AliInfo(Form("Entries (2D): %lf",(Double_t)(fHistPN->Project(0,1,2)->GetEntries())));
  
  //TCanvas *c2 = new TCanvas("c2","");
  //c2->cd();
  //fHistPN->Project(0,1,2)->DrawCopy("colz");

  if((Double_t)(fHistP->Project(0,1)->GetEntries())!=0)
    gHist->Scale(1./(Double_t)(fHistP->Project(0,1)->GetEntries()));

  //normalize to bin width
  gHist->Scale(1./((Double_t)gHist->GetXaxis()->GetBinWidth(1)*(Double_t)gHist->GetYaxis()->GetBinWidth(1)));
    
  return gHist;
}

//____________________________________________________________________//
TH2D *AliBalancePsi::GetCorrelationFunctionNP(Double_t psiMin, 
					      Double_t psiMax,
					      Double_t ptTriggerMin,
					      Double_t ptTriggerMax,
					      Double_t ptAssociatedMin,
					      Double_t ptAssociatedMax) {
  //Returns the 2D correlation function for (+-) pairs

  // security checks
  if(psiMin > psiMax-0.00001)                   psiMin          = psiMax-0.00001;
  if(ptTriggerMin > ptTriggerMax-0.00001)       ptTriggerMin    = ptTriggerMax-0.00001;
  if(ptAssociatedMin > ptAssociatedMax-0.00001) ptAssociatedMin = ptAssociatedMax-0.00001;

  // Psi_2: axis 0
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
    
  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistN->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.))
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);

  //0:step, 1: Delta eta, 2: Delta phi
  TH2D *gHist = dynamic_cast<TH2D *>(fHistNP->Project(0,1,2));
  if(!gHist){
    AliError("Projection of fHistPN = NULL");
    return gHist;
  }

  //Printf("Entries (1D): %lf",(Double_t)(fHistN->Project(0,2)->GetEntries()));
  //Printf("Entries (2D): %lf",(Double_t)(fHistNP->Project(0,2,3)->GetEntries()));
  if((Double_t)(fHistN->Project(0,1)->GetEntries())!=0)
    gHist->Scale(1./(Double_t)(fHistN->Project(0,1)->GetEntries()));

  //normalize to bin width
  gHist->Scale(1./((Double_t)gHist->GetXaxis()->GetBinWidth(1)*(Double_t)gHist->GetYaxis()->GetBinWidth(1)));
    
  return gHist;
}

//____________________________________________________________________//
TH2D *AliBalancePsi::GetCorrelationFunctionPP(Double_t psiMin, 
					      Double_t psiMax,
					      Double_t ptTriggerMin,
					      Double_t ptTriggerMax,
					      Double_t ptAssociatedMin,
					      Double_t ptAssociatedMax) {
  //Returns the 2D correlation function for (+-) pairs

  // security checks
  if(psiMin > psiMax-0.00001)                   psiMin          = psiMax-0.00001;
  if(ptTriggerMin > ptTriggerMax-0.00001)       ptTriggerMin    = ptTriggerMax-0.00001;
  if(ptAssociatedMin > ptAssociatedMax-0.00001) ptAssociatedMin = ptAssociatedMax-0.00001;

  // Psi_2: axis 0
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistP->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.))
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
      
  //0:step, 1: Delta eta, 2: Delta phi
  TH2D *gHist = dynamic_cast<TH2D *>(fHistPP->Project(0,1,2));
  if(!gHist){
    AliError("Projection of fHistPN = NULL");
    return gHist;
  }

  //Printf("Entries (1D): %lf",(Double_t)(fHistP->Project(0,2)->GetEntries()));
  //Printf("Entries (2D): %lf",(Double_t)(fHistPP->Project(0,2,3)->GetEntries()));
  if((Double_t)(fHistP->Project(0,1)->GetEntries())!=0)
    gHist->Scale(1./(Double_t)(fHistP->Project(0,1)->GetEntries()));

  //normalize to bin width
  gHist->Scale(1./((Double_t)gHist->GetXaxis()->GetBinWidth(1)*(Double_t)gHist->GetYaxis()->GetBinWidth(1)));
  
  return gHist;
}

//____________________________________________________________________//
TH2D *AliBalancePsi::GetCorrelationFunctionNN(Double_t psiMin, 
					      Double_t psiMax,
					      Double_t ptTriggerMin,
					      Double_t ptTriggerMax,
					      Double_t ptAssociatedMin,
					      Double_t ptAssociatedMax) {
  //Returns the 2D correlation function for (+-) pairs

  // security checks
  if(psiMin > psiMax-0.00001){
    AliError("psiMax <= psiMin");
    return NULL;
  }
  if(ptTriggerMin > ptTriggerMax-0.00001){
    AliError("ptTriggerMax <= ptTriggerMin");
    return NULL;
  }
  if(ptAssociatedMin > ptAssociatedMax-0.00001){
    AliError("ptAssociatedMax <= ptAssociatedMin");
    return NULL;
  }

  // Psi_2: axis 0
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistN->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.))
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    
  //0:step, 1: Delta eta, 2: Delta phi
  TH2D *gHist = dynamic_cast<TH2D *>(fHistNN->Project(0,1,2));
  if(!gHist){
    AliError("Projection of fHistPN = NULL");
    return gHist;
  }

  //Printf("Entries (1D): %lf",(Double_t)(fHistN->Project(0,2)->GetEntries()));
  //Printf("Entries (2D): %lf",(Double_t)(fHistNN->Project(0,2,3)->GetEntries()));
  if((Double_t)(fHistN->Project(0,1)->GetEntries())!=0)
    gHist->Scale(1./(Double_t)(fHistN->Project(0,1)->GetEntries()));

  //normalize to bin width
  gHist->Scale(1./((Double_t)gHist->GetXaxis()->GetBinWidth(1)*(Double_t)gHist->GetYaxis()->GetBinWidth(1)));
    
  return gHist;
}

//____________________________________________________________________//
TH2D *AliBalancePsi::GetCorrelationFunctionChargeIndependent(Double_t psiMin, 
							     Double_t psiMax,
							     Double_t ptTriggerMin,
							     Double_t ptTriggerMax,
							     Double_t ptAssociatedMin,
							     Double_t ptAssociatedMax) {
  //Returns the 2D correlation function for the sum of all charge combination pairs

  // security checks
  if(psiMin > psiMax-0.00001)                   psiMin          = psiMax-0.00001;
  if(ptTriggerMin > ptTriggerMax-0.00001)       ptTriggerMin    = ptTriggerMax-0.00001;
  if(ptAssociatedMin > ptAssociatedMax-0.00001) ptAssociatedMin = ptAssociatedMax-0.00001;

  // Psi_2: axis 0
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax-0.00001); 

  // pt trigger
  if((ptTriggerMin != -1.)&&(ptTriggerMax != -1.)) {
    fHistN->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistP->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRangeUser(ptTriggerMin,ptTriggerMax-0.00001);
  }

  // pt associated
  if((ptAssociatedMin != -1.)&&(ptAssociatedMax != -1.))
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRangeUser(ptAssociatedMin,ptAssociatedMax-0.00001);
    
  //0:step, 1: Delta eta, 2: Delta phi
  TH2D *gHistNN = dynamic_cast<TH2D *>(fHistNN->Project(0,1,2));
  if(!gHistNN){
    AliError("Projection of fHistNN = NULL");
    return gHistNN;
  }
  TH2D *gHistPP = dynamic_cast<TH2D *>(fHistPP->Project(0,1,2));
  if(!gHistPP){
    AliError("Projection of fHistPP = NULL");
    return gHistPP;
  }
  TH2D *gHistNP = dynamic_cast<TH2D *>(fHistNP->Project(0,1,2));
  if(!gHistNP){
    AliError("Projection of fHistNP = NULL");
    return gHistNP;
  }
  TH2D *gHistPN = dynamic_cast<TH2D *>(fHistPN->Project(0,1,2));
  if(!gHistPN){
    AliError("Projection of fHistPN = NULL");
    return gHistPN;
  }

  // sum all 2 particle histograms
  gHistNN->Add(gHistPP);
  gHistNN->Add(gHistNP);
  gHistNN->Add(gHistPN);

  // divide by sum of + and - triggers
  if((Double_t)(fHistN->Project(0,1)->GetEntries())!=0 && (Double_t)(fHistP->Project(0,1)->GetEntries())!=0)
    gHistNN->Scale(1./(Double_t)(fHistN->Project(0,1)->GetEntries() + fHistN->Project(0,1)->GetEntries()));

  //normalize to bin width
  gHistNN->Scale(1./((Double_t)gHistNN->GetXaxis()->GetBinWidth(1)*(Double_t)gHistNN->GetYaxis()->GetBinWidth(1)));
  
  return gHistNN;
}

//____________________________________________________________________//

Bool_t AliBalancePsi::GetMomentsAnalytical(TH1D* gHist,
					   Double_t &mean, Double_t &meanError,
					   Double_t &sigma, Double_t &sigmaError,
					   Double_t &skewness, Double_t &skewnessError,
					   Double_t &kurtosis, Double_t &kurtosisError) {
  //
  // helper method to calculate the moments and errors of a TH1D anlytically
  //
  
  Bool_t success = kFALSE;
  mean          = -1.;
  meanError     = -1.;
  sigma         = -1.;
  sigmaError    = -1.;
  skewness      = -1.;
  skewnessError = -1.;
  kurtosis      = -1.;
  kurtosisError = -1.;

  if(gHist){

    // ----------------------------------------------------------------------
    // basic parameters of histogram

    Int_t fNumberOfBins = gHist->GetNbinsX();
    //    Int_t fBinWidth     = gHist->GetBinWidth(1); // assume equal binning

    
    // ----------------------------------------------------------------------
    // first calculate the mean

    Double_t fWeightedAverage   = 0.;
    Double_t fNormalization = 0.;

    for(Int_t i = 1; i <= fNumberOfBins; i++) {

      fWeightedAverage   += gHist->GetBinContent(i) * gHist->GetBinCenter(i);
      fNormalization     += gHist->GetBinContent(i);

    }  
    
    mean = fWeightedAverage / fNormalization;


    // ----------------------------------------------------------------------
    // then calculate the higher moments

    Double_t fDelta  = 0.;
    Double_t fDelta2 = 0.;
    Double_t fDelta3 = 0.;
    Double_t fDelta4 = 0.;

    for(Int_t i = 1; i <= fNumberOfBins; i++) {
      fDelta  += gHist->GetBinContent(i) * gHist->GetBinCenter(i) - mean;
      fDelta2 += TMath::Power((gHist->GetBinContent(i) * gHist->GetBinCenter(i) - mean),2);
      fDelta3 += TMath::Power((gHist->GetBinContent(i) * gHist->GetBinCenter(i) - mean),3);
      fDelta4 += TMath::Power((gHist->GetBinContent(i) * gHist->GetBinCenter(i) - mean),4);
    }
    
    sigma    = fDelta2 / fNormalization;
    skewness = fDelta3 / fNormalization / TMath::Power(sigma,3);
    kurtosis = fDelta4 / fNormalization / TMath::Power(sigma,4) - 3;

    success = kTRUE;    
  }


  return success;
}


//____________________________________________________________________//
Float_t AliBalancePsi::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign) { 
  //
  // calculates dphistar
  //
  Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  
  static const Double_t kPi = TMath::Pi();
  
  // circularity
//   if (dphistar > 2 * kPi)
//     dphistar -= 2 * kPi;
//   if (dphistar < -2 * kPi)
//     dphistar += 2 * kPi;
  
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;
  
  return dphistar;
}

