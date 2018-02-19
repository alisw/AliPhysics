#include "TChain.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliCentrality.h"
#include "AliGenEventHeader.h"

#include "AliLog.h"
#include "AliAnalysisTaskPIDMCEffCont.h"

// ---------------------------------------------------------------------
ClassImp(AliAnalysisTaskPIDMCEffCont)

//________________________________________________________________________
AliAnalysisTaskPIDMCEffCont::AliAnalysisTaskPIDMCEffCont(const char *name)
: AliAnalysisTaskSE(name),
fAOD(0),
fArrayMC(0),
fQAList(0),
fOutputList(0),
fHistEventStats(0),
fHistCentrality(0),
fHistNMult(0),
fHistVz(0),
fHistPt(0),
fHistPhi(0),
fHistEta(0),
fHistNSigmaTPCvsPtbeforePID(0),
fHistNSigmaTPCvsPtafterPID(0),
fHistContaminationSecondariesPlus(0),
fHistContaminationSecondariesMinus(0), 
fHistContaminationPrimariesPlus(0),
fHistContaminationPrimariesMinus(0), 
fHistGeneratedEtaPtPhiPlus(0),
fHistSurvivedEtaPtPhiPlus(0),
fHistGeneratedEtaPtPhiMinus(0),
fHistSurvivedEtaPtPhiMinus(0),
fUseCentrality(kFALSE),
fCentralityEstimator("V0M"),
fCentralityPercentileMin(0.0),
fCentralityPercentileMax(5.0),
fInjectedSignals(kTRUE),
fElectronRejection(kFALSE),
fElectronOnlyRejection(kFALSE),
fElectronRejectionNSigma(-1.),
fElectronRejectionMinPt(0.),
fElectronRejectionMaxPt(1000.),
fVxMax(3.0),
fVyMax(3.0),
fVzMax(10.),
fAODTrackCutBit(128),
fMinNumberOfTPCClusters(80),
fMaxChi2PerTPCCluster(4.0),
fMaxDCAxy(3.0),
fMaxDCAz(3.0),
fPtMin(0.0),
fPtMax(20.0),
fEtaMin(-0.8),
fEtaMax(0.8),
fUsePIDNsigma(kFALSE),
fUsePIDPDG(kFALSE),
fDCAextended(kFALSE),
fPIDResponse(0),
fPIDNSigma(3),
fParticleOfInterest(kPion),
fHistDCAXYptprimminus(0),
fHistDCAXYptprimplus(0),
fHistDCAXYptsecminusweak(0),
fHistDCAXYptsecplusweak(0),
fHistDCAXYptsecminusmat(0),
fHistDCAXYptsecplusmat(0),
fHistDCAXYptchargedminus(0),
fHistDCAXYptchargedplus(0),
fHistDCAXYptprimminus_ext(0),
fHistDCAXYptprimplus_ext(0),
fHistDCAXYptsecminusweak_ext(0),
fHistDCAXYptsecplusweak_ext(0),
fHistDCAXYptsecminusmat_ext(0),
fHistDCAXYptsecplusmat_ext(0),
fHistDCAXYptchargedminus_ext(0),
fHistDCAXYptchargedplus_ext(0)
{
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}

void AliAnalysisTaskPIDMCEffCont::UserCreateOutputObjects() {
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
    Int_t ptBin = 34;
    Int_t etaBin = 16;
    Int_t phiBin = 100;
    
    Double_t nArrayPt[35]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0};
    Double_t nArrayEta[17]={-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    
    Double_t nArrayPhi[phiBin+1];
    for(Int_t iBin = 0; iBin <= phiBin; iBin++)
        nArrayPhi[iBin] = iBin*TMath::TwoPi()/phiBin;
    
    Int_t detaBin = 16;
    Int_t dphiBin = 100;
    Double_t nArrayDPhi[dphiBin+1];
    for(Int_t iBin = 0; iBin <= dphiBin; iBin++)
        nArrayDPhi[iBin] = iBin*TMath::TwoPi()/dphiBin;
    Double_t nArrayDEta[17]={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
    //====================================================//
    
    //AOD analysis
    fHistCentrality = new TH1F("fHistCentrality",";Centrality bin;Events",
                               1001,-0.5,100.5);
    fQAList->Add(fHistCentrality);
    
    //multiplicity (good MC tracks)
    TString histName;
    histName = "fHistNMult";
    fHistNMult = new TH1F(histName.Data(),
                          ";N_{mult.}",
                          200,0,20000);
    fQAList->Add(fHistNMult);
    
    //Vz addition+++++++++++++++++++++++++++++
    //fHistVz = new TH1F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Entries",100,-20.,20.);
    fHistVz = new TH2F("fHistVz","Primary vertex distribution - z coordinate;Centrality (%);V_{z} (cm);Entries",10,0,100,100,-20.,20.);
    fQAList->Add(fHistVz);
 
    fHistPt = new TH1F("fHistPt","Pt distribution;p_{T} (GeV/c);Counts",100,0,10);
    fQAList->Add(fHistPt);

    fHistPhi = new TH1F("fHistPhi","Phi distribution;Phi;Number Of Entries",100,0,2*(TMath::Pi()));
    fQAList->Add(fHistPhi);
  
    fHistEta = new TH1F("fHistEta","Eta distribution;Eta;Number Of Entries",100,-1,1);
    fQAList->Add(fHistEta);
   
    fHistContaminationSecondariesPlus = new TH3D("fHistContaminationSecondariesPlus","Secondaries;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
    fOutputList->Add(fHistContaminationSecondariesPlus);
    
    fHistContaminationSecondariesMinus = new TH3D("fHistContaminationSecondariesMinus","Secondaries;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
    fOutputList->Add(fHistContaminationSecondariesMinus);
    
    //Contamination for Primaries
    fHistContaminationPrimariesPlus = new TH3D("fHistContaminationPrimariesPlus","Primaries;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
    fOutputList->Add(fHistContaminationPrimariesPlus);
    
    fHistContaminationPrimariesMinus = new TH3D("fHistContaminationPrimariesMinus","Primaries;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
    fOutputList->Add(fHistContaminationPrimariesMinus);
    
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

    fHistDCAXYptprimminus = new TH3F("fHistDCAxyptprimminus","DCA_{xy} vs pt for primaries (negative);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-0.5,0.5);
    fHistDCAXYptprimplus = new TH3F("fHistDCAxyptprimplus","DCA_{xy} vs pt for primaries (positive);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-0.5,0.5);
    fHistDCAXYptsecminusweak = new TH3F("fHistDCAxyptsecminusweak","DCA_{xy} vs pt for seondaries from weak decays (negative);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-0.5,0.5);
    fHistDCAXYptsecplusweak = new TH3F("fHistDCAxyptsecplusweak","DCA_{xy} vs pt for secondaries from weak decays (positive);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-0.5,0.5);
    fHistDCAXYptsecminusmat = new TH3F("fHistDCAxyptsecminusmat","DCA_{xy} vs pt for secondaries from material (negative);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-0.5,0.5);
    fHistDCAXYptsecplusmat = new TH3F("fHistDCAxyptsecplusmat","DCA_{xy} vs pt for secondaries from material (positive);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-0.5,0.5);

    fHistDCAXYptchargedminus = new TH3F("fHistDCAxychargedminus","DCA_{xy} vs pt for charged particles (negative);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);
    fHistDCAXYptchargedplus = new TH3F("fHistDCAxychargedplus","DCA_{xy} vs pt for charged particles (positive);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);
     
    fHistDCAXYptprimminus_ext = new TH3F("fHistDCAxyptprimminusext","DCA_{xy} vs pt for primaries (negative);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);
    fHistDCAXYptprimplus_ext = new TH3F("fHistDCAxyptprimplusext","DCA_{xy} vs pt for primaries (positive);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);
    fHistDCAXYptsecminusweak_ext = new TH3F("fHistDCAxyptsecminusweakext","DCA_{xy} vs pt for seondaries from weak decays (negative);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);
    fHistDCAXYptsecplusweak_ext = new TH3F("fHistDCAxyptsecplusweakext","DCA_{xy} vs pt for secondaries from weak decays (positive);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);
    fHistDCAXYptsecminusmat_ext = new TH3F("fHistDCAxyptsecminusmatext","DCA_{xy} vs pt for secondaries from material (negative);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);
    fHistDCAXYptsecplusmat_ext = new TH3F("fHistDCAxyptsecplusmatext","DCA_{xy} vs pt for secondaries from material (positive);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);

    fHistDCAXYptchargedminus_ext = new TH3F("fHistDCAxychargedminusext","DCA_{xy} vs pt for charged particles (negative);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);
    fHistDCAXYptchargedplus_ext = new TH3F("fHistDCAxychargedplusext","DCA_{xy} vs pt for charged particles (positive);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);

    fOutputList->Add(fHistDCAXYptprimminus);
    fOutputList->Add(fHistDCAXYptprimplus);
    fOutputList->Add(fHistDCAXYptsecminusweak);
    fOutputList->Add(fHistDCAXYptsecplusweak);
    fOutputList->Add(fHistDCAXYptsecminusmat);
    fOutputList->Add(fHistDCAXYptsecplusmat);
    fOutputList->Add(fHistDCAXYptchargedminus);
    fOutputList->Add(fHistDCAXYptchargedplus);  
    fOutputList->Add(fHistDCAXYptprimminus_ext);
    fOutputList->Add(fHistDCAXYptprimplus_ext);
    fOutputList->Add(fHistDCAXYptsecminusweak_ext);
    fOutputList->Add(fHistDCAXYptsecplusweak_ext);
    fOutputList->Add(fHistDCAXYptsecminusmat_ext);
    fOutputList->Add(fHistDCAXYptsecplusmat_ext);
    fOutputList->Add(fHistDCAXYptchargedminus_ext);
    fOutputList->Add(fHistDCAXYptchargedplus_ext);    
     
    //fQAList->Print();
    //fOutputList->Print(); 
 
    PostData(1, fQAList);
    PostData(2, fOutputList);
}



void AliAnalysisTaskPIDMCEffCont::UserExec(Option_t*)
{
    // event processing
   
    
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
        AliError("ERROR: Could not retrieve MC event");
        return;
    }
    
    
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
    Double_t nCentrality = 0;
    if(fUseCentrality) {
        
        AliMultSelection *MultSelection = 0x0;
        MultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
        if( !MultSelection) {
            //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
            AliWarning("AliMultSelection object not found!");
        }
        else{
            nCentrality = MultSelection->GetMultiplicityPercentile(fCentralityEstimator.Data());
        }
    }
           
    if(fUsePIDNsigma) {
        fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
        if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");
    }    


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
                            fHistVz->Fill(nCentrality,vertex->GetZ());
     	
                    if((nCentrality >= fCentralityPercentileMin) &&
     			    (nCentrality < fCentralityPercentileMax)){

               		     fHistEventStats->Fill(2); //triggered + centrality
                	     fHistCentrality->Fill(nCentrality);		    
	                       
                            //++++++++++++++++++CONTAMINATION++++++++++++++++++//
                            Int_t nGoodAODTracks = fAOD->GetNumberOfTracks();
                            Int_t nMCParticles = mcEvent->GetNumberOfTracks();
                            TArrayI labelMCArray(nMCParticles);
			                                
                            for(Int_t jTracks = 0; jTracks < nGoodAODTracks; jTracks++) {
                                AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(jTracks));
                                if(!track) AliFatal("Not a standard AOD");
                                if(!track) continue;
  
                                if (!track->TestFilterBit(fAODTrackCutBit))
                                    continue;
                                
    
                               //acceptance
                                if(TMath::Abs(track->Eta()) > fEtaMax)
                                    continue;
                                if((track->Pt() > fPtMax)||(track->Pt() <  fPtMin))
                                    continue;
                                
                                Double_t phiRad = track->Phi();
                                
                                Int_t label = TMath::Abs(track->GetLabel());
                                if(label > nMCParticles) continue;
                                AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) mcEvent->GetTrack(label);
                                Short_t gAODmcCharge = AODmcTrack->Charge();
                                
                                // =============================================================================================
                                // Skip tracks that come from injected signals
                                
				if (fInjectedSignals)
                                {
                                    
                                TString genname_cont;
                                Bool_t yesno=mcEvent->GetCocktailGenerator(label,genname_cont);
                                //if(yesno) Printf("Reconstructed: cocktail header name is %s", genname_cont.Data());
                                if(!genname_cont.Contains("Hijing_0"))
                                continue;
                                //Printf("Reconstructed: cocktail header name is %s", genname_cont.Data());
                                }
                                
                                Float_t ptreco  = track->Pt();
                                Float_t eta = track->Eta();
                                Float_t phi = track->Phi();
                                Int_t IDrec = AODmcTrack->PdgCode();
				                                
				Double_t  dca[2] = {0.0,0.0};
                                Double_t  cov[3] = {0.0,0.0,0.0};
		
				AliAODTrack copy(*track);

				if (fAODTrackCutBit==768){
                                if (track->TestFilterBit(256)){
                                if (!copy.PropagateToDCA(vertex,fAOD->GetMagneticField(),300.,dca,cov))
                                continue;
                                }
                                }

                                else {
                                if (!copy.PropagateToDCA(vertex,fAOD->GetMagneticField(),300.,dca,cov))
                                continue;
                                }				

                                if(fUsePIDNsigma) {

                                Float_t probMis = fPIDResponse->GetTOFMismatchProbability(track);
                                if (probMis < 0.01) { //if u want to re duce mismatch using also TPC

                                    Double_t nSigmaPionTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
                                    Double_t nSigmaKaonTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
                                    Double_t nSigmaProtonTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));
                                    
                                    Double_t nSigmaPionTOF   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion));
                                    Double_t nSigmaKaonTOF   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon));
                                    Double_t nSigmaProtonTOF = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton));
                                    
                                    
                                    if (fParticleOfInterest == kPion){
                                        
                                        if( ptreco > 0.2 && ptreco < 0.6 ){
                                            if( nSigmaPionTPC > fPIDNSigma || nSigmaKaonTPC < fPIDNSigma || nSigmaProtonTPC < fPIDNSigma )
                                                continue;
                                        }
                                        else if(ptreco > 0.6){
                                            if( nSigmaPionTPC > fPIDNSigma )
                                                continue;
                                            
                                            if( nSigmaPionTOF > fPIDNSigma || nSigmaKaonTOF < fPIDNSigma || nSigmaProtonTOF < fPIDNSigma )
                                                continue;
                                        }
                                    } //end of pion case
                                    
                                    
                                    else if(fParticleOfInterest == kKaon){
                                        
                                        if( ptreco > 0.2 && ptreco < 0.4 ){
                                            
                                            if( nSigmaPionTPC < fPIDNSigma || nSigmaKaonTPC > fPIDNSigma || nSigmaProtonTPC < fPIDNSigma )
                                                continue;
                                        }
                                        else if(ptreco >= 0.4  && ptreco <= 2.5){
                                            
                                            if( nSigmaKaonTPC > fPIDNSigma )
                                                continue;
                                            
                                            if( nSigmaPionTOF < fPIDNSigma || nSigmaKaonTOF > fPIDNSigma || nSigmaProtonTOF < fPIDNSigma )
                                                continue;
                                        }
                                    } //end of the kaon case
                                    
                                    
                                    else if (fParticleOfInterest == kProton){
                         
                                        if( ptreco > 0.2 && ptreco < 0.6 ){
                                            if( nSigmaPionTPC < fPIDNSigma || nSigmaKaonTPC < fPIDNSigma || nSigmaProtonTPC > fPIDNSigma )
                                                continue;
                                        }
                                        
                                        else if(ptreco > 0.6  && ptreco < 4.0 ){
                                            if( nSigmaProtonTPC > fPIDNSigma )
                                                
                                                continue;
                                            
                                            if( nSigmaPionTOF < fPIDNSigma || nSigmaKaonTOF < fPIDNSigma || nSigmaProtonTOF > fPIDNSigma )
                                                
                                                continue;
                                        }
                                    }//end of the proton case
                                   
                                }//end of probability case
                                
                                }//end of nsigma PID
                              
                                
                                if(fUsePIDPDG) {
                            
					if (fParticleOfInterest == kProton){
                                        if(TMath::Abs(IDrec)!=2212)
					continue;	
					}
                                
					else if (fParticleOfInterest == kKaon){
                                        
					if(TMath::Abs(IDrec)!=321)
                                        continue;
				    	}
                                    
					else if (fParticleOfInterest == kPion){
                                        
					if(TMath::Abs(IDrec)!=211)
                                        continue;
					}
                                }//end of PDG PID
                                
                                fHistPhi->Fill(phi);
                                fHistEta->Fill(eta);
                                fHistPt->Fill(ptreco);                                

                                if(gAODmcCharge < 0){
                                    
					if (fAODTrackCutBit==768){
                                        	if (fDCAextended){
                                            		if (track->TestFilterBit(512))
                                                		fHistDCAXYptchargedminus_ext->Fill(track->Pt(),track->Eta(),track->DCA());
                                            		else if (track->TestFilterBit(256))
                                                		fHistDCAXYptchargedminus_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                        	}
                                        	else if (!fDCAextended) {
                                                	if (track->TestFilterBit(512))
                                                		fHistDCAXYptchargedminus->Fill(track->Pt(),track->Eta(),track->DCA());
                                                	else if (track->TestFilterBit(256))
                                                       		fHistDCAXYptchargedminus->Fill(track->Pt(),track->Eta(),dca[0]);
                                       		}
                                    	}
                                 	else {
                                        	if (fDCAextended)
                                            		fHistDCAXYptchargedminus_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                        	else if (!fDCAextended)
                                            		fHistDCAXYptchargedminus->Fill(track->Pt(),track->Eta(),dca[0]);
                                    	}
                                }
                                
				else if (gAODmcCharge > 0){
					if (fAODTrackCutBit==768){
                                        	if (fDCAextended){
                                            		if (track->TestFilterBit(512))
                                                		fHistDCAXYptchargedplus_ext->Fill(track->Pt(),track->Eta(),track->DCA());
                                            		else if (track->TestFilterBit(256))
                                                		fHistDCAXYptchargedplus_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                        	}
                                        	else if (!fDCAextended){
                                        		if (track->TestFilterBit(512))
                                                		fHistDCAXYptchargedplus->Fill(track->Pt(),track->Eta(),track->DCA());
                                            		else if (track->TestFilterBit(256))
                                                		fHistDCAXYptchargedplus->Fill(track->Pt(),track->Eta(),dca[0]);
                                        	}
                                    	}	
                                    	else {
                                    		if (fDCAextended)
                                            		fHistDCAXYptchargedplus_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                        	else if (!fDCAextended)
                                            		fHistDCAXYptchargedplus->Fill(track->Pt(),track->Eta(),dca[0]);
                                    	}
                                }
                                
                                if (AODmcTrack->IsPhysicalPrimary()) {
					if(gAODmcCharge > 0){
                                        	fHistContaminationPrimariesPlus->Fill(track->Eta(),track->Pt(),phiRad);
                                            		if (fAODTrackCutBit==768){
                                                		if (fDCAextended){
                                                    			if (track->TestFilterBit(512))
                                                        			fHistDCAXYptprimplus_ext->Fill(track->Pt(),track->Eta(),track->DCA());
                                                    			else if (track->TestFilterBit(256))
                                                        			fHistDCAXYptprimplus_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                		}		
                                                		else if (!fDCAextended){
                                                    			if (track->TestFilterBit(512))
                                                        			fHistDCAXYptprimplus->Fill(track->Pt(),track->Eta(),track->DCA());
                                                    			else if (track->TestFilterBit(256))
                                                        			fHistDCAXYptprimplus->Fill(track->Pt(),track->Eta(),dca[0]);
                                                		}
                                            		}	
                                            		else {
                                                		if (fDCAextended)
                                                    			fHistDCAXYptprimplus_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                		else if (!fDCAextended)
                                              	     			fHistDCAXYptprimplus->Fill(track->Pt(),track->Eta(),dca[0]);
                                            		}
                                    	}
                                    
                                     	else if(gAODmcCharge < 0){
				    		fHistContaminationPrimariesMinus->Fill(track->Eta(),track->Pt(),phiRad);
                                        		if (fAODTrackCutBit==768){
                                                		if (fDCAextended){
                                                    			if (track->TestFilterBit(512))
                                                        			fHistDCAXYptprimminus_ext->Fill(track->Pt(),track->Eta(),track->DCA());
                                                    			else if (track->TestFilterBit(256))
                                                        			fHistDCAXYptprimminus_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                		}
                                                		else if (!fDCAextended){
                                                    			if (track->TestFilterBit(512))
                                                        			fHistDCAXYptprimminus->Fill(track->Pt(),track->Eta(),track->DCA());
                                                    			else if (track->TestFilterBit(256))
                                                        			fHistDCAXYptprimminus->Fill(track->Pt(),track->Eta(),dca[0]);
                                                		}
                                            		}
                                            		else {
                                                		if (fDCAextended)
                                                    			fHistDCAXYptprimminus_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                		else if (!fDCAextended)
                                                    			fHistDCAXYptprimminus->Fill(track->Pt(),track->Eta(),dca[0]);
                                            		}
                                       }
                                }
                                
                                else {
                                
					if(gAODmcCharge > 0){
                                        	fHistContaminationSecondariesPlus->Fill(track->Eta(),track->Pt(),phiRad);
                                        		if(AODmcTrack->IsSecondaryFromWeakDecay()){
	                                        		if (fAODTrackCutBit==768){
                                                			if (fDCAextended){
                                                    				if (track->TestFilterBit(512))
                                                        				fHistDCAXYptsecplusweak_ext->Fill(track->Pt(),track->Eta(),track->DCA());
                                                    				else if (track->TestFilterBit(256))
                                                        				fHistDCAXYptsecplusweak_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                			}
                                                			else if (!fDCAextended){
                                                    				if (track->TestFilterBit(512))
                                                        				fHistDCAXYptsecplusweak->Fill(track->Pt(),track->Eta(),track->DCA());
                                                    				else if (track->TestFilterBit(256))
                                                        				fHistDCAXYptsecplusweak->Fill(track->Pt(),track->Eta(),dca[0]);
                                                			}
                                            			}
                                            			else {
                                                			if (fDCAextended)
                                                    				fHistDCAXYptsecplusweak_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                			else if (!fDCAextended)
                                                    				fHistDCAXYptsecplusweak->Fill(track->Pt(),track->Eta(),dca[0]);
                                           			}
                                        		}

                                       		 	else if(AODmcTrack->IsSecondaryFromMaterial()){
                                            			if (fAODTrackCutBit==768){
                                                			if (fDCAextended){
                                                    				if (track->TestFilterBit(512))
                                                        				fHistDCAXYptsecplusmat_ext->Fill(track->Pt(),track->Eta(),track->DCA());
                                                    				else if (track->TestFilterBit(256))
                                                        				fHistDCAXYptsecplusmat_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                			}
                                                			else if (!fDCAextended){
                                                    				if (track->TestFilterBit(512))
                                                        				fHistDCAXYptsecplusmat->Fill(track->Pt(),track->Eta(),track->DCA());
                                                    				else if (track->TestFilterBit(256))
                                                        				fHistDCAXYptsecplusmat->Fill(track->Pt(),track->Eta(),dca[0]);
                                               				}			
                                            			}
                                            			else {
                                                			if (fDCAextended)
                                                    				fHistDCAXYptsecplusmat_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                			else if (!fDCAextended)
                                                    				fHistDCAXYptsecplusmat->Fill(track->Pt(),track->Eta(),dca[0]);
                                            			}
                                        		}
                                   	}
                                    
                                   	else if(gAODmcCharge < 0){
                                        	fHistContaminationSecondariesMinus->Fill(track->Eta(),track->Pt(),phiRad);
                                           		if(AODmcTrack->IsSecondaryFromWeakDecay()){
                                                		if (fAODTrackCutBit==768){
                                                    			if (fDCAextended){
                                                        			if (track->TestFilterBit(512))
                                                            				fHistDCAXYptsecminusweak_ext->Fill(track->Pt(),track->Eta(),track->DCA());
                                                        			else if (track->TestFilterBit(256))
                                                      			      		fHistDCAXYptsecminusweak_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                    			}		
                                                    			else if (!fDCAextended){
                                                        			if (track->TestFilterBit(512))
                                                            				fHistDCAXYptsecminusweak->Fill(track->Pt(),track->Eta(),track->DCA());
                                                        			else if (track->TestFilterBit(256))
                                                            				fHistDCAXYptsecminusweak->Fill(track->Pt(),track->Eta(),dca[0]);
                                                    			}			
                                                		}
                                                		else {
                                                    			if (fDCAextended)
                                                        			fHistDCAXYptsecminusweak_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                    			else if (!fDCAextended)
                                                        			fHistDCAXYptsecminusweak->Fill(track->Pt(),track->Eta(),dca[0]);
                                                		}
                                            		}	
                                       
                                            		else if(AODmcTrack->IsSecondaryFromMaterial()){
                                                		if (fAODTrackCutBit==768){
                                                    			if (fDCAextended){
                                                        			if (track->TestFilterBit(512))
                                                            				fHistDCAXYptsecminusmat_ext->Fill(track->Pt(),track->Eta(),track->DCA());
                                                        			else if (track->TestFilterBit(256))
                                                            				fHistDCAXYptsecminusmat_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                    			}		
                                                   	 		else if (!fDCAextended){
                                                        			if (track->TestFilterBit(512))
                                                            				fHistDCAXYptsecminusmat->Fill(track->Pt(),track->Eta(),track->DCA());
                                                        			else if (track->TestFilterBit(256))
                                                            				fHistDCAXYptsecminusmat->Fill(track->Pt(),track->Eta(),dca[0]);
                                                    			}
                                                		}
                                                		else {
                                                    			if (fDCAextended)
                                                        			fHistDCAXYptsecminusmat_ext->Fill(track->Pt(),track->Eta(),dca[0]);
                                                    			else if (!fDCAextended)
                                                        			fHistDCAXYptsecminusmat->Fill(track->Pt(),track->Eta(),dca[0]);
                                                		}
                                            		}
                                   	}
                                }
                            }//loop over tracks
                            

                            //++++++++++++++++++EFFICIENCY+++++++++++++++++++++//
                            for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
                                AliAODMCParticle *mcTrack = (AliAODMCParticle*) mcEvent->GetTrack(iTracks);
                                if (!mcTrack) {
                                    AliError(Form("ERROR: Could not receive track %d (mc loop)", iTracks));
                                    continue;
                                }
                                
                                if (fInjectedSignals){
                                TString genname;
                                Bool_t yesno=mcEvent->GetCocktailGenerator(iTracks,genname);
                                //if(!yesno) Printf("Generated: no cocktail header list was found for this event");
                                //if(yesno) Printf("Generated: cocktail header name is %s", genname.Data());
                                if(!genname.Contains("Hijing_0"))
                                continue;
                                //Printf("Generated: cocktail header name is %s", genname.Data());
                                }
			

                                //exclude particles generated out of the acceptance
                                Double_t vz = mcTrack->Zv();
                                if (TMath::Abs(vz) > 50.) continue;
                                //acceptance
                                if(TMath::Abs(mcTrack->Eta()) > fEtaMax)
                                    continue;
                                if((mcTrack->Pt() > fPtMax)||(mcTrack->Pt() < fPtMin))
                                    continue;
                                
                                if(!mcTrack->IsPhysicalPrimary()) continue;
                                
                                Short_t gMCCharge = mcTrack->Charge();
                                Double_t phiRad = mcTrack->Phi();
                        
                                Int_t IDgen = mcTrack->PdgCode();
				
                                if(fUsePIDNsigma || fUsePIDPDG) {

                                    if (fParticleOfInterest == kProton){
					
                                        if(TMath::Abs(IDgen)!=2212)
                                        continue;
          
                                    }//end of proton case
                                    
                                    else if (fParticleOfInterest == kKaon){
                                        
					if(TMath::Abs(IDgen)!=321)
                                        continue;
                                        
                                    }//end of kaon case
                                    
                                    else if (fParticleOfInterest == kPion){
                                       
                                        if(TMath::Abs(IDgen)!=211)
                                        continue;
                               	    }//end of pion case
                                
				}//end PID
                                
                                    
                                if(gMCCharge > 0)
                                
                                      fHistGeneratedEtaPtPhiPlus->Fill(mcTrack->Eta(),
                                                                     mcTrack->Pt(),
                                                                     phiRad);
                                else if(gMCCharge < 0)
                                      fHistGeneratedEtaPtPhiMinus->Fill(mcTrack->Eta(),
                                                                      mcTrack->Pt(),
                                                                      phiRad);
                                
                                
                                
                                    
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
                                    phi[nMCLabelCounter]    = mcTrack->Phi();
                                    charge[nMCLabelCounter] = gMCCharge;
                                    
                                    level[nMCLabelCounter]  = 1;
                                    nMCLabelCounter += 1;
                                }  
                            }//loop over MC particles
                            
                            fHistNMult->Fill(nMCLabelCounter);
                            
                            //AOD track loop
                            Int_t nGoodTracks = fAOD->GetNumberOfTracks();   
                            TArrayI labelArray(nGoodTracks);
                            Int_t labelCounter = 0;
                            
                            for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
                                AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));    
                                if(!trackAOD) continue;
                                
                                //track cuts
                                if (!trackAOD->TestFilterBit(fAODTrackCutBit)) 
                                    continue;
                                
                                Int_t label = TMath::Abs(trackAOD->GetLabel()); 
                                if(IsLabelUsed(labelArray,label)) continue;
                                labelArray.AddAt(label,labelCounter);
                                labelCounter += 1;
                                
				AliAODMCParticle *AODmcTracksurv = (AliAODMCParticle*) mcEvent->GetTrack(label);

                                Int_t mcGoods = nMCLabelCounter;
                                for (Int_t k = 0; k < mcGoods; k++) {
                                    Int_t mcLabel = labelMCArray.At(k);
                                    
                                    if (mcLabel != TMath::Abs(label)) continue;
                                    if(mcLabel != label) continue;		    
               
                                    //acceptance
                                    if(TMath::Abs(trackAOD->Eta()) > fEtaMax) 
                                        continue;
                                    if((trackAOD->Pt() > fPtMax)||(trackAOD->Pt() <  fPtMin)) 
                                        continue;

                                    Short_t gCharge = trackAOD->Charge();
                                    Double_t phiRad = trackAOD->Phi();
                		    Float_t ptsurv  = trackAOD->Pt();
                      		    Int_t IDsurv = AODmcTracksurv->PdgCode();

                                    if(fUsePIDNsigma) {
                                 
					Float_t probMis = fPIDResponse->GetTOFMismatchProbability(trackAOD);
                                	if (probMis < 0.01) { 
       
                                        Double_t nSigmaPionTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,AliPID::kPion));
                                        Double_t nSigmaKaonTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,AliPID::kKaon));
                                        Double_t nSigmaProtonTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,AliPID::kProton));
                                        
                                        Double_t nSigmaPionTOF   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackAOD,AliPID::kPion));
                                        Double_t nSigmaKaonTOF   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackAOD,AliPID::kKaon));
                                        Double_t nSigmaProtonTOF = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackAOD,AliPID::kProton));
                                        
                                        
                                        if (fParticleOfInterest == kPion){
                                            
                                            if( ptsurv > 0.2 && ptsurv < 0.6 ){
                                                if( nSigmaPionTPC > fPIDNSigma || nSigmaKaonTPC < fPIDNSigma || nSigmaProtonTPC < fPIDNSigma )
                                                    continue;
                                            }
                                            else if(ptsurv > 0.6){
                                                if( nSigmaPionTPC > fPIDNSigma )
                                                    continue;
                                                
                                                if( nSigmaPionTOF > fPIDNSigma || nSigmaKaonTOF < fPIDNSigma || nSigmaProtonTOF < fPIDNSigma )
                                                    continue;
                                            }
                                        } //end of pion case
                                        
                                        
                                        else if(fParticleOfInterest == kKaon){
                                            
                                            if( ptsurv > 0.2 && ptsurv < 0.4 ){
                                                
                                                if( nSigmaPionTPC < fPIDNSigma || nSigmaKaonTPC > fPIDNSigma || nSigmaProtonTPC < fPIDNSigma )
                                                    continue;
                                            }
                                            else if(ptsurv >= 0.4  && ptsurv <= 2.5){
                                                
                                                if( nSigmaKaonTPC > fPIDNSigma )
                                                    continue;
                                                
                                                if( nSigmaPionTOF < fPIDNSigma || nSigmaKaonTOF > fPIDNSigma || nSigmaProtonTOF < fPIDNSigma )
                                                    continue;
                                            }
                                        } //end of the kaon case
                                        
                                        
                                        else if (fParticleOfInterest == kProton){
                                 
                                            if( ptsurv > 0.2 && ptsurv < 0.6 ){
                                                if( nSigmaPionTPC < fPIDNSigma || nSigmaKaonTPC < fPIDNSigma || nSigmaProtonTPC > fPIDNSigma )
                                                    continue;
                                            }
                                            
                                            else if(ptsurv > 0.6  && ptsurv < 4.0 ){
                                                if( nSigmaProtonTPC > fPIDNSigma )
                                                    
                                                    continue;
                                                
                                                if( nSigmaPionTOF < fPIDNSigma || nSigmaKaonTOF < fPIDNSigma || nSigmaProtonTOF > fPIDNSigma )
                                                    
                                                    continue;
                                            }
                                        }//end of the proton case
                                        
					}//end of probability check

					}//end of nsigma PID
                                    
					if(fUsePIDPDG) {

                                        if (fParticleOfInterest == kProton){
                                        if(TMath::Abs(IDsurv)!=2212)
                                        continue;
                                        }

                                        else if (fParticleOfInterest == kKaon){

                                        if(TMath::Abs(IDsurv)!=321)
                                        continue;
                                        }

                                        else if (fParticleOfInterest == kPion){

                                        if(TMath::Abs(IDsurv)!=211)
                                        continue;
                                        }
 	                                }//end of PDG PID

					      
                                        if(TMath::Abs(trackAOD->Eta()) < fEtaMax && trackAOD->Pt() > fPtMin&&trackAOD->Pt() < fPtMax){ 
                                        level[k]  = 2;
                                        
                                        if(gCharge > 0)
                                            fHistSurvivedEtaPtPhiPlus->Fill(trackAOD->Eta(),
                                                                            trackAOD->Pt(),
                                                                            phiRad);
                                        else if(gCharge < 0)
                                            fHistSurvivedEtaPtPhiMinus->Fill(trackAOD->Eta(),
                                                                             trackAOD->Pt(),
                                                                             phiRad);
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
    
}
}

void AliAnalysisTaskPIDMCEffCont::Terminate(Option_t *) {
    // Draw result to the screen
    // Called once at the end of the query
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available\n");
        return;
    }
}

Bool_t AliAnalysisTaskPIDMCEffCont::IsLabelUsed(TArrayI labelArray, Int_t label) {
    //Checks if the label is used already
    Bool_t status = kFALSE;
    for(Int_t i = 0; i < labelArray.GetSize(); i++) {
        if(labelArray.At(i) == label)
            status = kTRUE;
    }
    
    return status;
}

    
    
    


