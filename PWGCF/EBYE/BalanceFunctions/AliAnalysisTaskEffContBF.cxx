#include "TChain.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h" 
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"

#include "AliAnalysisTaskEffContBF.h"

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency of the Balance Function 
// for single particles and pairs
// 
// ---------------------------------------------------------------------

ClassImp(AliAnalysisTaskEffContBF)

//________________________________________________________________________
AliAnalysisTaskEffContBF::AliAnalysisTaskEffContBF(const char *name) 
  : AliAnalysisTaskSE(name), 
    fAOD(0),
    fAODtrackCutBit(128),
    fArrayMC(0), 
    fQAList(0), 
    fOutputList(0), 
    fHistEventStats(0), 
    fHistCentrality(0),
    fHistNMult(0), 
    fHistVz(0), 
    fHistContaminationSecondaries(0),
    fHistContaminationPrimaries(0),
    fHistGeneratedEtaPtPhiPlus(0), 
    fHistSurvivedEtaPtPhiPlus(0),
    fHistGeneratedEtaPtPhiMinus(0),
    fHistSurvivedEtaPtPhiMinus(0),
    fHistGeneratedEtaPtPlusControl(0),
    fHistSurvivedEtaPtPlusControl(0),
    fHistGeneratedEtaPtMinusControl(0),
    fHistSurvivedEtaPtMinusControl(0),
    fHistGeneratedEtaPtPlusPlus(0),
    fHistSurvivedEtaPtPlusPlus(0),
    fHistGeneratedEtaPtMinusMinus(0),
    fHistSurvivedEtaPtMinusMinus(0),
    fHistGeneratedEtaPtPlusMinus(0),
    fHistSurvivedEtaPtPlusMinus(0),
    fHistGeneratedPhiEtaPlusPlus(0),
    fHistSurvivedPhiEtaPlusPlus(0),
    fHistGeneratedPhiEtaMinusMinus(0),
    fHistSurvivedPhiEtaMinusMinus(0),
    fHistGeneratedPhiEtaPlusMinus(0),
    fHistSurvivedPhiEtaPlusMinus(0),
    fCentralityEstimator("V0M"), 
    fCentralityPercentileMin(0.0), 
    fCentralityPercentileMax(5.0), 
    fVxMax(3.0), 
    fVyMax(3.0),
    fVzMax(10.), 
    fMinNumberOfTPCClusters(80),
    fMaxChi2PerTPCCluster(4.0),
    fMaxDCAxy(3.0),
    fMaxDCAz(3.0),
    fMinPt(0.0),
    fMaxPt(20.0),
    fMinEta(-0.8), 
    fMaxEta(0.8),
    fEtaRangeMin(0.0), 
    fEtaRangeMax(1.6), 
    fPtRangeMin(0.0), 
    fPtRangeMax(20.0), 
    fPhiRangeMin(0.0),
    fPhiRangeMax(360.),
    fdPhiRangeMax(180.), 
    fEtaBin(100), //=100 (BF) 16
    fdEtaBin(64), //=64 (BF)  16
    fPtBin(100), //=100 (BF)  36
    fPhiBin(100), //=100 (BF)
    fdPhiBin(90), //=90 (BF)
    fHistSurvived4EtaPtPhiPlus(0),
    fHistSurvived8EtaPtPhiPlus(0)
  
{   
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEffContBF::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  fQAList = new TList();
  fQAList->SetName("QAList");
  fQAList->SetOwner();
  
  fOutputList = new TList();
  fOutputList->SetName("OutputList");
  fOutputList->SetOwner();
  
  //Event stats.
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fQAList->Add(fHistEventStats);

  //====================================================//
  Int_t ptBin = 36;
  Int_t etaBin = 16;
  Int_t phiBin = 100;

  Double_t nArrayPt[37]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0};
  Double_t nArrayEta[17]={-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}; 
  Double_t nArrayPhi[101]={0, 3.6, 7.2, 10.8, 14.4, 18., 21.6, 25.2, 28.8, 32.4, 36., 39.6, 43.2, 46.8, 50.4, 54., 57.6, 61.2, 64.8, 68.4, 72., 75.6, 79.2, 82.8, 86.4, 90., 93.6, 97.2, 100.8, 104.4, 108., 111.6, 115.2, 118.8, 122.4, 126., 129.6, 133.2, 136.8, 140.4, 144., 147.6, 151.2, 154.8, 158.4, 162., 165.6, 169.2, 172.8, 176.4, 180., 183.6, 187.2, 190.8, 194.4, 198., 201.6, 205.2, 208.8, 212.4, 216., 219.6, 223.2, 226.8, 230.4, 234, 237.6, 241.2, 244.8, 248.4, 252., 255.6, 259.2, 262.8, 266.4, 270., 273.6, 277.2, 280.8, 284.4, 288., 291.6, 295.2, 298.8, 302.4, 306., 309.6, 313.2, 316.8, 320.4, 324., 327.6, 331.2, 334.8, 338.4, 342., 345.6, 349.2, 352.8, 356.4, 360.};

  Int_t detaBin = 16;
  Int_t dphiBin = 90;
  Double_t nArrayDPhi[91]={0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,176,178,180};
  Double_t nArrayDEta[17]={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6}; 
  //====================================================//

  //AOD analysis
  fHistCentrality = new TH1F("fHistCentrality",";Centrality bin;Events",
			     20,0.5,20.5);
  fQAList->Add(fHistCentrality);
  
  //multiplicity (good MC tracks)
  TString histName;
  histName = "fHistNMult";
  fHistNMult = new TH1F(histName.Data(), 
			";N_{mult.}",
			200,0,20000);
  fQAList->Add(fHistNMult);

  //Vz addition+++++++++++++++++++++++++++++
  fHistVz = new TH1F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Entries",100,-20.,20.);
  fQAList->Add(fHistVz);

  //Contamination for Secondaries 
  fHistContaminationSecondaries = new TH3D("fHistContaminationSecondaries","Secondaries;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistContaminationSecondaries);

  //Contamination for Primaries
  fHistContaminationPrimaries = new TH3D("fHistContaminationPrimaries","Primaries;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistContaminationPrimaries);
  
  //eta vs pt for MC positives
  fHistGeneratedEtaPtPhiPlus = new TH3D("fHistGeneratedEtaPtPhiPlus",
					"Generated positive primaries;#eta;p_{T} (GeV/c);#phi",
					etaBin,nArrayEta, ptBin, nArrayPt,phiBin, nArrayPhi);
  // fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax,fPhiBin,nArrayPhi);
  fOutputList->Add(fHistGeneratedEtaPtPhiPlus);
  fHistSurvivedEtaPtPhiPlus = new TH3D("fHistSurvivedEtaPtPhiPlus",
				       "Survived positive primaries;#eta;p_{T} (GeV/c);#phi",
				       etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistSurvivedEtaPtPhiPlus);
  
  //eta vs pt for MC negatives
  fHistGeneratedEtaPtPhiMinus = new TH3D("fHistGeneratedEtaPtPhiMinus",
					 "Generated positive primaries;#eta;p_{T} (GeV/c);#phi",
					 etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistGeneratedEtaPtPhiMinus);
  fHistSurvivedEtaPtPhiMinus = new TH3D("fHistSurvivedEtaPtPhiMinus",
					"Survived positive primaries;#eta;p_{T} (GeV/c);#phi",
					etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistSurvivedEtaPtPhiMinus);
 
  //eta vs pt for MC positives (control)
  fHistGeneratedEtaPtPlusControl = new TH2F("fHistGeneratedEtaPtPlusControl",
					    "Generated positive primaries;#eta;p_{T} (GeV/c)",
					    etaBin,nArrayEta,ptBin,nArrayPt);
  fOutputList->Add(fHistGeneratedEtaPtPlusControl);
  fHistSurvivedEtaPtPlusControl = new TH2F("fHistSurvivedEtaPtPlusControl",
					   "Survived positive primaries;#eta;p_{T} (GeV/c)",
					   etaBin,nArrayEta,ptBin,nArrayPt);
  fOutputList->Add(fHistSurvivedEtaPtPlusControl);
  
  //eta vs pt for MC negatives (control)
  fHistGeneratedEtaPtMinusControl = new TH2F("fHistGeneratedEtaPtMinusControl",
					     "Generated positive primaries;#eta;p_{T} (GeV/c)",
					     etaBin,nArrayEta,ptBin,nArrayPt);
  fOutputList->Add(fHistGeneratedEtaPtMinusControl);
  fHistSurvivedEtaPtMinusControl = new TH2F("fHistSurvivedEtaPtMinusControl",
					    "Survived positive primaries;#eta;p_{T} (GeV/c)",
					    etaBin,nArrayEta,ptBin,nArrayPt);
  fOutputList->Add(fHistSurvivedEtaPtMinusControl);
  
  //eta vs pt for MC ++
  fHistGeneratedEtaPtPlusPlus = new TH2F("fHistGeneratedEtaPtPlusPlus",
					 "Generated ++ primaries;#Delta#eta;p_{T} (GeV/c)",
					 detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistGeneratedEtaPtPlusPlus);
  fHistSurvivedEtaPtPlusPlus = new TH2F("fHistSurvivedEtaPtPlusPlus",
					"Survived ++ primaries;#Delta#eta;p_{T} (GeV/c)",
					detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistSurvivedEtaPtPlusPlus);
  
  //eta vs pt for MC --
  fHistGeneratedEtaPtMinusMinus = new TH2F("fHistGeneratedEtaPtMinusMinus",
					   "Generated -- primaries;#Delta#eta;p_{T} (GeV/c)",
					   detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistGeneratedEtaPtMinusMinus);
  fHistSurvivedEtaPtMinusMinus = new TH2F("fHistSurvivedEtaPtMinusMinus",
					  "Survived -- primaries;#Delta#eta;p_{T} (GeV/c)",
					  detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistSurvivedEtaPtMinusMinus);
 
  //eta vs pt for MC +-
  fHistGeneratedEtaPtPlusMinus = new TH2F("fHistGeneratedEtaPtPlusMinus",
					  "Generated +- primaries;#Delta#eta;p_{T} (GeV/c)",
					  detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistGeneratedEtaPtPlusMinus);
  fHistSurvivedEtaPtPlusMinus = new TH2F("fHistSurvivedEtaPtPlusMinus",
					 "Survived +- primaries;#Delta#eta;p_{T} (GeV/c)",
					 detaBin,nArrayDEta,ptBin,nArrayPt);
  fOutputList->Add(fHistSurvivedEtaPtPlusMinus);
 
  //=============================//
  //phi vs eta for MC ++
  fHistGeneratedPhiEtaPlusPlus = new TH2F("fHistGeneratedPhiEtaPlusPlus",
					  "Generated ++ primaries;#Delta#phi",
					  dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistGeneratedPhiEtaPlusPlus);
  fHistSurvivedPhiEtaPlusPlus = new TH2F("fHistSurvivedPhiEtaPlusPlus",
					 "Survived ++ primaries;#Delta#phi;#Delta#eta",
					 dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistSurvivedPhiEtaPlusPlus);
  
  //phi vs eta for MC --
  fHistGeneratedPhiEtaMinusMinus = new TH2F("fHistGeneratedPhiEtaMinusMinus",
					    "Generated -- primaries;#Delta#phi;#Delta#eta",
					    dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistGeneratedPhiEtaMinusMinus);
  fHistSurvivedPhiEtaMinusMinus = new TH2F("fHistSurvivedPhiEtaMinusMinus",
					   "Survived -- primaries;#Delta#phi;#Delta#eta",
					   dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistSurvivedPhiEtaMinusMinus);
  
  //phi vs eta for MC +-
  fHistGeneratedPhiEtaPlusMinus = new TH2F("fHistGeneratedPhiEtaPlusMinus",
					   "Generated +- primaries;#Delta#phi;#Delta#eta",
					   dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistGeneratedPhiEtaPlusMinus);
  fHistSurvivedPhiEtaPlusMinus = new TH2F("fHistSurvivedPhiEtaPlusMinus",
					  "Survived +- primaries;#Delta#phi;#Delta#eta",
					  dphiBin,nArrayDPhi,detaBin,nArrayDEta);
  fOutputList->Add(fHistSurvivedPhiEtaPlusMinus);
  
  fHistSurvived4EtaPtPhiPlus = new TH3F("fHistSurvived4EtaPtPhiPlus",
                                        "Survived4 + primaries;#eta;p_{T} (GeV/c);#phi",
                                        etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistSurvived4EtaPtPhiPlus);
  fHistSurvived8EtaPtPhiPlus = new TH3F("fHistSurvived8EtaPtPhiPlus",
					"Survived8 + primaries;#eta;p_{T} (GeV/c);#phi",
					etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistSurvived8EtaPtPhiPlus);
  
  //fQAList->Print();
  //fOutputList->Print(); 
  PostData(1, fQAList);
  PostData(2, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskEffContBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;
  }

  fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));   
  if (!fArrayMC)  
    AliFatal("No array of MC particles found !!!"); // MW  no AliFatal use return values   
  
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
  
  // arrays for 2 particle histograms
  Int_t nMCLabelCounter         = 0;
  const Int_t maxMCLabelCounter = 20000;
  
  Double_t eta[maxMCLabelCounter];
  Double_t pt[maxMCLabelCounter];
  Double_t phi[maxMCLabelCounter];
  Int_t level[maxMCLabelCounter];
  Int_t charge[maxMCLabelCounter];
  
  //AliInfo(Form("%d %d",mcEvent->GetNumberOfTracks(),fAOD->GetNumberOfTracks()));

  fHistEventStats->Fill(1); //all events
  
  //Centrality stuff
  AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
  Double_t nCentrality = 0;
  AliCentrality *centrality = header->GetCentralityP();
  
  if(centrality->IsEventInCentralityClass(fCentralityPercentileMin,
					  fCentralityPercentileMax,
					  fCentralityEstimator.Data())) {
    //if(centrality){
    nCentrality =centrality->GetCentralityPercentile(fCentralityEstimator.Data());
    //Printf("Centrality: %lf",centrality->GetCentralityPercentile(fCentralityEstimator.Data()));
    
    fHistEventStats->Fill(2); //triggered + centrality
    fHistCentrality->Fill(nCentrality+1);
    
    //Printf("Centrality selection: %lf - %lf",fCentralityPercentileMin,fCentralityPercentileMax);
    
    if(fAnalysisMode.CompareTo("TPC") == 0 ) {
      const AliAODVertex *vertex = fAOD->GetPrimaryVertex(); 
      if(vertex) {
	if(vertex->GetNContributors() > 0) {
	  Double32_t fCov[6];    
	  vertex->GetCovarianceMatrix(fCov);   
	  if(fCov[5] != 0) {
	    fHistEventStats->Fill(3); //events with a proper vertex
	    if(TMath::Abs(vertex->GetX()) < fVxMax) {    // antes Xv
	      //Printf("X Vertex: %lf", vertex->GetX());
	      //Printf("Y Vertex: %lf", vertex->GetY());
	      if(TMath::Abs(vertex->GetY()) < fVyMax) {  // antes Yv
		if(TMath::Abs(vertex->GetZ()) < fVzMax) {  // antes Zv
		  //Printf("Z Vertex: %lf", vertex->GetZ());
 
		  fHistEventStats->Fill(4); //analyzed events
		  fHistVz->Fill(vertex->GetZ()); 
		  
		  //++++++++++++++++++CONTAMINATION++++++++++++++++++//
		  Int_t nGoodAODTracks = fAOD->GetNumberOfTracks();
		  Int_t nMCParticles = mcEvent->GetNumberOfTracks();
		  TArrayI labelMCArray(nMCParticles);
		  
		  for(Int_t jTracks = 0; jTracks < nGoodAODTracks; jTracks++) {
		    AliAODTrack* track = fAOD->GetTrack(jTracks);
		    if(!track) continue;
		    
		    if (!track->TestFilterBit(fAODtrackCutBit))
		      continue;
		    
		    //acceptance
		    if(TMath::Abs(track->Eta()) > fMaxEta) 
		      continue;
		    if((track->Pt() > fMaxPt)||(track->Pt() <  fMinPt)) 
		      continue;
		    if((track->Phi() > fPhiRangeMax)||(track->Phi() < fPhiRangeMin)) 
		      continue;
		    
		    Double_t phiRad = track->Phi(); 
		    Double_t phiDeg = phiRad*TMath::RadToDeg();
		    
		    Int_t label = TMath::Abs(track->GetLabel());
		    if(label > nMCParticles) continue;
		    AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) mcEvent->GetTrack(label); 
		    //fHistContaminationPrimaries->Fill(track->Eta(),track->Pt(),phiDeg);
		    //if (!(AODmcTrack->IsPhysicalPrimary())) {
		    //fHistContaminationSecondaries->Fill(track->Eta(),track->Pt(),phiDeg);
		    //}
		    if (AODmcTrack->IsPhysicalPrimary()) {
		      fHistContaminationPrimaries->Fill(track->Eta(),track->Pt(),phiDeg);
		    }
		    else{
		      fHistContaminationSecondaries->Fill(track->Eta(),track->Pt(),phiDeg);
		    }
		  }
		  //++++++++++++++++++CONTAMINATION++++++++++++++++++//
		  
		  //++++++++++++++++++EFFICIENCY+++++++++++++++++++++//
		  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
		    AliAODMCParticle *mcTrack = (AliAODMCParticle*) mcEvent->GetTrack(iTracks); 
		    if (!mcTrack) {
		      Printf("ERROR: Could not receive track %d (mc loop)", iTracks);
		      continue;
		    }
		    
		    //exclude particles generated out of the acceptance
		    Double_t vz = mcTrack->Zv();
		    if (TMath::Abs(vz) > 50.) continue;
		    //acceptance
		    if(TMath::Abs(mcTrack->Eta()) > fMaxEta) 
		      continue;
		    if((mcTrack->Pt() > fMaxPt)||(mcTrack->Pt() < fMinPt)) 
		      continue;
		    if((mcTrack->Phi() > fPhiRangeMax)||(mcTrack->Phi() < fPhiRangeMin)) 
		      continue;
		    
		    if(!mcTrack->IsPhysicalPrimary()) continue;   
		    
		    Short_t gMCCharge = mcTrack->Charge();
		    //Double_t phiRad = particle->Phi(); //antes
		    Double_t phiRad = mcTrack->Phi(); 
		    Double_t phiDeg = phiRad*TMath::RadToDeg();
		    
		    if(gMCCharge > 0)
		      fHistGeneratedEtaPtPhiPlus->Fill(mcTrack->Eta(),
						       mcTrack->Pt(),
						       phiDeg);
		    else if(gMCCharge < 0)
		      fHistGeneratedEtaPtPhiMinus->Fill(mcTrack->Eta(),
							mcTrack->Pt(),
							phiDeg);
		    
		    Bool_t labelTPC = kTRUE;
		    if(labelTPC) {
		      labelMCArray.AddAt(iTracks,nMCLabelCounter);
		      if(nMCLabelCounter >= maxMCLabelCounter){
			AliWarning(Form("MC Label Counter > Limit (%d) --> stop loop here",maxMCLabelCounter));
			break;
		      }
		      //fill the arrays for 2 particle analysis
		      eta[nMCLabelCounter]    = mcTrack->Eta();
		      pt[nMCLabelCounter]     = mcTrack->Pt();
		      phi[nMCLabelCounter]    = mcTrack->Phi()*TMath::RadToDeg();
		      charge[nMCLabelCounter] = gMCCharge;
		      
		      level[nMCLabelCounter]  = 1;
		      nMCLabelCounter += 1;
		    }  
		    //}//primaries //antes
		  }//loop over MC particles
		  
		  fHistNMult->Fill(nMCLabelCounter);
		  
		  //AOD track loop
		  Int_t nGoodTracks = fAOD->GetNumberOfTracks();   
		  TArrayI labelArray(nGoodTracks);
		  Int_t labelCounter = 0;
		  
		  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {		  
		    AliAODTrack *trackAOD = static_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));    
		    if(!trackAOD) continue;
		    
		    if (!trackAOD->TestFilterBit(fAODtrackCutBit)) 
		      continue; 
		    
		    Int_t label = TMath::Abs(trackAOD->GetLabel()); 
		    if(IsLabelUsed(labelArray,label)) continue;
		    labelArray.AddAt(label,labelCounter);
		    labelCounter += 1;
		    
		    Bool_t iFound = kFALSE;
		    Int_t mcGoods = nMCLabelCounter;
		    for (Int_t k = 0; k < mcGoods; k++) {
		      Int_t mcLabel = labelMCArray.At(k);
		      iFound = kFALSE;
		      
		      if (mcLabel != TMath::Abs(label)) continue;
		      if(mcLabel != label) continue;		    
		      if(label > trackAOD->GetLabel()) continue; 
		      
		      //acceptance
		      if(TMath::Abs(trackAOD->Eta()) > fMaxEta) 
			continue;
		      if((trackAOD->Pt() > fMaxPt)||(trackAOD->Pt() <  fMinPt)) 
			continue;
		      if((trackAOD->Phi() > fPhiRangeMax)||(trackAOD->Phi() < fPhiRangeMin)) 
			continue;
		      
		      level[k]  = 2;
		      
		      Short_t gCharge = trackAOD->Charge();
		      Double_t phiRad = trackAOD->Phi();
		      Double_t phiDeg = phiRad*TMath::RadToDeg();
		      
		      if(TMath::Abs(trackAOD->Eta()) < fMaxEta && trackAOD->Pt() > fMinPt&&trackAOD->Pt() < fMaxPt){ 
			
			level[k]  = 3;
			
			if(gCharge > 0)
			  fHistSurvivedEtaPtPhiPlus->Fill(trackAOD->Eta(),
							  trackAOD->Pt(),
							  phiDeg);
			else if(gCharge < 0)
			  fHistSurvivedEtaPtPhiMinus->Fill(trackAOD->Eta(),
							   trackAOD->Pt(),
							   phiDeg);
		      }//tracks		   
		    }//end of mcGoods
		  }//AOD track loop
		  
		  labelMCArray.Reset();
		  labelArray.Reset();	       
		  
		}//Vz cut
	      }//Vy cut
	    }//Vx cut
	  }//Vz resolution
	}//number of contributors
      }//valid vertex
    }//TPC analysis mode
  }//centrality  // NO ESTABA
  
  
  // Here comes the 2 particle analysis
  // loop over all good MC particles
  for (Int_t i = 0; i < nMCLabelCounter ; i++) {
    
    // control 1D histograms (charge might be different?)
    if(charge[i] > 0){
      if(level[i] > 0) fHistGeneratedEtaPtPlusControl->Fill(eta[i],pt[i]);
      if(level[i] > 2) fHistSurvivedEtaPtPlusControl->Fill(eta[i],pt[i]);
    }
    else if(charge[i] < 0){
      if(level[i] > 0) fHistGeneratedEtaPtMinusControl->Fill(eta[i],pt[i]);
      if(level[i] > 2) fHistSurvivedEtaPtMinusControl->Fill(eta[i],pt[i]);
    }
    
    
    for (Int_t j = i+1; j < nMCLabelCounter ; j++) {
      
      if(charge[i] > 0 && charge[j] > 0 ){
	if(level[i] > 0 && level[j] > 0) {  
	  fHistGeneratedEtaPtPlusPlus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	    fHistGeneratedPhiEtaPlusPlus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	}
	if(level[i] > 2 && level[j] > 2) {
	  fHistSurvivedEtaPtPlusPlus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	    fHistSurvivedPhiEtaPlusPlus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	}
      }
      
      else if(charge[i] < 0 && charge[j] < 0 ){
	if(level[i] > 0 && level[j] > 0) {
	  fHistGeneratedEtaPtMinusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);	    
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	    fHistGeneratedPhiEtaMinusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));  	    
	}
	if(level[i] > 2 && level[j] > 2) {
	  fHistSurvivedEtaPtMinusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	    fHistSurvivedPhiEtaMinusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	}
      }
      
      else if((charge[i] > 0 && charge[j] < 0)||(charge[i] < 0 && charge[j] > 0)){
	if(level[i] > 0 && level[j] > 0) {
	  fHistGeneratedEtaPtPlusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	    fHistGeneratedPhiEtaPlusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));	
	}
	if(level[i] > 2 && level[j] > 2) {
	  fHistSurvivedEtaPtPlusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	    fHistSurvivedPhiEtaPlusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	}	
      }
    }
  }
  
  
  // } //antes
  
  //Checking Entries for generated and survived
  /*TH1F* GeneratedEtaPt = (TH1F*)fHistGeneratedEtaPtPlusMinus->ProjectionX("GeneratedEtaPt",5,15);
    Int_t xGeneratedPt = fHistGeneratedEtaPtPlusMinus->GetNbinsX();
    for (Int_t h=1; h < xGeneratedPt+1; h++){
    Double_t binEntriesGenerated = GeneratedEtaPt->GetBinContent(h);
    Printf("binEntriesGenerated: %lf - xGeneratedPt: %d",binEntriesGenerated,h);  
    }
    TH1F* GeneratedPhiEta = (TH1F*)fHistGeneratedPhiEtaPlusMinus->ProjectionY("GeneratedPhiEta",5,15);
    Int_t yGeneratedPhi = fHistGeneratedPhiEtaPlusMinus->GetNbinsY();
    for (Int_t h=1; h < yGeneratedPhi+1; h++){
    Double_t binEntriesGenerated = GeneratedPhiEta->GetBinContent(h);
    Printf("binEntriesGenerated: %lf - yGeneratedPhi: %d",binEntriesGenerated,h);  
    }*/
  
  /*TH1F* SurvivedEtaPt = (TH1F*)fHistSurvivedEtaPtPlusMinus->ProjectionX("SurvivedEtaPt",5,15);
    Int_t xSurvivedPt = fHistSurvivedEtaPtPlusMinus->GetNbinsX();
    for (Int_t h=1; h < xSurvivedPt+1; h++){
    Double_t binEntriesSurvived = SurvivedEtaPt->GetBinContent(h);
    Printf("binEntriesSurvived: %lf - xSurvivedPt: %d",binEntriesSurvived,h);
    }
    TH1F* SurvivedPhiEta = (TH1F*)fHistSurvivedPhiEtaPlusMinus->ProjectionY("SurvivedPhiEta",5,15);
    Int_t ySurvivedPhi = fHistSurvivedPhiEtaPlusMinus->GetNbinsY();
    for (Int_t h=1; h < ySurvivedPhi+1; h++){
    Double_t binEntriesSurvived = SurvivedPhiEta->GetBinContent(h);
    Printf("binEntriesSurvived: %lf - ySurvivedPhi: %d",binEntriesSurvived,h);
    }*/  
}

//________________________________________________________________________
void AliAnalysisTaskEffContBF::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//____________________________________________________________________//
Bool_t AliAnalysisTaskEffContBF::IsLabelUsed(TArrayI labelArray, Int_t label) {
  //Checks if the label is used already
  Bool_t status = kFALSE;
  for(Int_t i = 0; i < labelArray.GetSize(); i++) {
    if(labelArray.At(i) == label)
      status = kTRUE;
  }

  return status;
}
