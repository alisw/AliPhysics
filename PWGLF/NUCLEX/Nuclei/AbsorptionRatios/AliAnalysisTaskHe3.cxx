#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TProfile.h"
#include "TNtuple.h"

#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliInputEventHandler.h>
#include <AliAODTrack.h>
#include "AliPIDResponse.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliAODHeader.h"

#include "AliAnalysisTaskHe3.h"

class AliAnalysisTaskHe3;
using namespace std;

ClassImp(AliAnalysisTaskHe3)

//________________________________________________________________________
AliAnalysisTaskHe3::AliAnalysisTaskHe3()
: AliAnalysisTaskSE()
,fVEvent(0)
,fAODEvent(0)
,fOutputList(0)
,fOutputEvent(0)
,fOutputProtons(0)
,fOutputAProtons(0)
,fOutputDeuterons(0)
,fOutputADeuterons(0)
,fHistTrackCuts(0)
,fEventStat(0)
,fVertexZ(0)
,fVtxContrib(0)
,fSPDVtxResol(0)
,fVtxDisplacement(0)
,fMultV0(0)
,fCentV0M(0)
,fCentV0MZoomed(0)
,fNch(0)
,fNchHeader(0)
,fNtupleHe3(0)
,fNtupleAHe3(0)
,fHistsProton()
,fHistsAProton()
,fHistsDeuteron()
,fHistsADeuteron()
,fPIDResponse(0)
,fFilterBit(256)
,fLowPCut(0.0)
,fHighPCut(1e30)
,fEtaCut(1.0)
,fMinClIts(0)
,fMaxDCAxyCut(10.)
,fMaxDCAxyAna(10.)
,fMaxDCAz(10.)
,fMaxTPCnSigma(10.)
,fMaxTOFnSigma(10.)
,fUseTOFPidCut(kFALSE)
,fMomTOFProt(10.)
,fMomTOFDeut(10.)
,kAnalyseAllParticles(kTRUE)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//________________________________________________________________________
AliAnalysisTaskHe3::AliAnalysisTaskHe3(const char* name)
: AliAnalysisTaskSE(name)
,fVEvent(0)
,fAODEvent(0)
,fOutputList(0)
,fOutputEvent(0)
,fOutputProtons(0)
,fOutputAProtons(0)
,fOutputDeuterons(0)
,fOutputADeuterons(0)
,fHistTrackCuts(0)
,fEventStat(0)
,fVertexZ(0)
,fVtxContrib(0)
,fSPDVtxResol(0)
,fVtxDisplacement(0)
,fMultV0(0)
,fCentV0M(0)
,fCentV0MZoomed(0)
,fNch(0)
,fNchHeader(0)
,fNtupleHe3(0)
,fNtupleAHe3(0)
,fHistsProton()
,fHistsAProton()
,fHistsDeuteron()
,fHistsADeuteron()
,fPIDResponse(0)
,fFilterBit(256)
,fLowPCut(0.0)
,fHighPCut(1e30)
,fEtaCut(1.0)
,fMinClIts(0)
,fMaxDCAxyCut(10.)
,fMaxDCAxyAna(10.)
,fMaxDCAz(10.)
,fMaxTPCnSigma(10.)
,fMaxTOFnSigma(10.)
,fUseTOFPidCut(kFALSE)
,fMomTOFProt(10.)
,fMomTOFDeut(10.)
,kAnalyseAllParticles(kTRUE)
{
    // constructor
    
    fHistsProton.clear();
    fHistsAProton.clear();
    fHistsDeuteron.clear();
    fHistsADeuteron.clear();
	fHistsHe3.clear();
	fHistsAHe3.clear();
    
    DefineInput(0, TChain::Class());    // a 'chain' of events created by the analysis manager automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
    // one can add more output objects by calling DefineOutput(2, classname::Class())
    // make sure to connect them properly in AddTask and to call PostData() for all of them!
    
}

//________________________________________________________________________
AliAnalysisTaskHe3::~AliAnalysisTaskHe3()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;
    }
}

//________________________________________________________________________
void AliAnalysisTaskHe3::UserCreateOutputObjects()
{
    // Called once at the start of the analysis
    
    // retrieve PID object from the input handler
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    // create output lists
    fOutputList = new TList();
    fOutputList->SetName("Output");
    fOutputList->SetOwner(kTRUE); // the list is owner of all objects it contains and will delete them
    
    fOutputEvent = new TList();
    fOutputEvent->SetName("Event");
    fOutputEvent->SetOwner(kTRUE);
    
    fOutputProtons = new TList();
    fOutputProtons->SetName("Protons");
    fOutputProtons->SetOwner(kTRUE);
    
    fOutputAProtons = new TList();
    fOutputAProtons->SetName("AProtons");
    fOutputAProtons->SetOwner(kTRUE);
    
    fOutputDeuterons = new TList();
    fOutputDeuterons->SetName("Deuterons");
    fOutputDeuterons->SetOwner(kTRUE);
    
    fOutputADeuterons = new TList();
    fOutputADeuterons->SetName("ADeuterons");
    fOutputADeuterons->SetOwner(kTRUE);

	fOutputHe3 = new TList();
	fOutputHe3->SetName("He3");
	fOutputHe3->SetOwner(kTRUE);

	fOutputAHe3 = new TList();
	fOutputAHe3->SetName("AHe3");
	fOutputAHe3->SetOwner(kTRUE);
    kAnalyseAllParticles=kTRUE;
    fOutputList->Add(fOutputEvent);
	if (kAnalyseAllParticles) {
			fOutputList->Add(fOutputProtons);
			fOutputList->Add(fOutputAProtons);
			fOutputList->Add(fOutputDeuterons);
			fOutputList->Add(fOutputADeuterons);
	}
	fOutputList->Add(fOutputHe3);
	fOutputList->Add(fOutputAHe3);
    
    // create event histograms
    fEventStat = new TH1I("fEventStat", "Event statistics", kNbinsEvent, 0, kNbinsEvent);
    fEventStat->GetXaxis()->SetBinLabel(1,"After Phys. Sel. and trigger");
    fEventStat->GetXaxis()->SetBinLabel(2,"kINT7 selected (cross-check)");
    fEventStat->GetXaxis()->SetBinLabel(3,"DAQ imcomplete (cross-check)");
    fEventStat->GetXaxis()->SetBinLabel(4,"V0 timing (cross-check)");
    fEventStat->GetXaxis()->SetBinLabel(5,"Cluster-tracklet cut (cross-check)");
    fEventStat->GetXaxis()->SetBinLabel(6,"Vertex Z position");
    fEventStat->GetXaxis()->SetBinLabel(7,"N of vtx contrib");
    fEventStat->GetXaxis()->SetBinLabel(8,"SPD pile-up in mult. bins");
    fEventStat->GetXaxis()->SetBinLabel(9,"Vertex displacement");
    fEventStat->GetXaxis()->SetBinLabel(10,"Vertex Z resolution");
    fOutputEvent->Add(fEventStat);
    
    fVertexZ = new TH1F("fVertexZ", "Vertex Z distribution", 300, -15., 15.);
    fVertexZ->GetXaxis()->SetTitle("vertex Z, cm");
    fVertexZ->Sumw2();
    fOutputEvent->Add(fVertexZ);
    
    fVtxContrib = new TH1F("fVtxContrib", "N vertex contributors", 200, 0., 200.);
    fVtxContrib->GetXaxis()->SetTitle("N vertex contrib.");
    fOutputEvent->Add(fVtxContrib);
    
    fSPDVtxResol = new TH1F("fSPDVtxResol", "SPD vertex resolution", 100, 0., 1.);
    fSPDVtxResol->GetXaxis()->SetTitle("SPD vertex resolution, cm");
    fOutputEvent->Add(fSPDVtxResol);
    
    fVtxDisplacement = new TH2F("fVtxDisplacement", "SPD vertex displacement", 300, -15., 15.,300, -15., 15.);
    fVtxDisplacement->GetXaxis()->SetTitle("SPD vertex position, cm");
    fVtxDisplacement->GetYaxis()->SetTitle("Track vertex position, cm");
    fOutputEvent->Add(fVtxDisplacement);
    
    fMultV0 = new TH1F("fMultV0", "Mult V0 distrobution", 2000, 0., 2000.);
    fMultV0->GetXaxis()->SetTitle("V0 mult.");
    fOutputEvent->Add(fMultV0);
    
    fCentV0M = new TH1F("fCentV0M", "Centrality V0M", 351, -50., 301.);
    fCentV0M->GetXaxis()->SetTitle("V0M percentile");
    fOutputEvent->Add(fCentV0M);

	fCentV0MZoomed = new TH1F("fCentV0MZoomed", "Centrality V0M", 100, 0.0, 1.0);
	fCentV0MZoomed->GetXaxis()->SetTitle("V0M percentile");
	fOutputEvent->Add(fCentV0MZoomed);

	fNch = new TH1F("fNch", "Nch tracks per event", 1001, -0.5, 1000.5);
	fNch->GetXaxis()->SetTitle("Nch tracks per event");
	fOutputEvent->Add(fNch);

	fNchHeader = new TH1F("fNchHeader", "Nch tracks per event from header", 1001, -0.5, 1000.5);
	fNch->GetXaxis()->SetTitle("Nch tracks per event");
	fOutputEvent->Add(fNchHeader);


	fNtupleHe3 = new TNtuple("fNtupleHe3", "fNtupleHe3", "p:pt:TPCSignal:TPCnSigmaHe3:Charge:TOF_beta:DCAxy:DCAz:TOFm2:TPCNClusters:ITSNClusters:TPCClusters4dEdx:Eta:ITSnSigmaHe3:Chi2TPC:Chi2ITS:TPCCrossedRows:pTPC");
	fNtupleAHe3 = new TNtuple("fNtupleAHe3", "fNtupleAHe3", "p:pt:TPCSignal:TPCnSigmaHe3:Charge:TOF_beta:DCAxy:DCAz:TOFm2:TPCNClusters:ITSNClusters:TPCClusters4dEdx:Eta:ITSnSigmaHe3:Chi2TPC:Chi2ITS:TPCCrossedRows:pTPC");


    // track cuts config
    fHistTrackCuts = new TProfile("fHistTrackCuts", "TrackCuts config", 10, 0, 10);
    fHistTrackCuts->GetXaxis()->SetBinLabel(1,"Filter Bit");
    fHistTrackCuts->GetXaxis()->SetBinLabel(2,"low #it{p} cut");
    fHistTrackCuts->GetXaxis()->SetBinLabel(3,"#eta cut");
    fHistTrackCuts->GetXaxis()->SetBinLabel(4,"Min. N ITS cl");
    fHistTrackCuts->GetXaxis()->SetBinLabel(5,"Max. DCAxy cut");
    fHistTrackCuts->GetXaxis()->SetBinLabel(6,"Max. DCAxy ana");
    fHistTrackCuts->GetXaxis()->SetBinLabel(7,"Max. DCAz");
    fHistTrackCuts->GetXaxis()->SetBinLabel(8,"TPCnSigma cut");
    fHistTrackCuts->GetXaxis()->SetBinLabel(9,"TOFnSigma cut");
    fHistTrackCuts->GetXaxis()->SetBinLabel(10,"TOFmomentum prot");
    fHistTrackCuts->GetXaxis()->SetBinLabel(11,"TOFmomentum deut");
    fOutputList->Add(fHistTrackCuts);
    
    fHistTrackCuts->Fill(0.0,fFilterBit);
    fHistTrackCuts->Fill(1,fLowPCut);
    fHistTrackCuts->Fill(2,fEtaCut);
    fHistTrackCuts->Fill(3,fMinClIts);
    fHistTrackCuts->Fill(4,fMaxDCAxyCut);
    fHistTrackCuts->Fill(5,fMaxDCAxyAna);
    fHistTrackCuts->Fill(6,fMaxDCAz);
    fHistTrackCuts->Fill(7,fMaxTPCnSigma);
    fHistTrackCuts->Fill(8,fMaxTOFnSigma);
    fHistTrackCuts->Fill(9,fMomTOFProt);
    fHistTrackCuts->Fill(10,fMomTOFDeut);
    
    // (anti)proton histograms
	// CreateHistosTrack() creates a list of histograms which are added to the vector to be accessed later. Since the addition of more histograms changes the numberings, these have to be handled carefully.
    if (kAnalyseAllParticles) {
			CreateHistosTrack(fHistsProton);  
			for (Int_t i=0;i<(int)(fHistsProton.size()); i++){
				fOutputProtons->Add(fHistsProton.at(i));
			}
			
			CreateHistosTrack(fHistsAProton);
			for (Int_t i=0;i<(int)(fHistsProton.size()); i++){
				fOutputAProtons->Add(fHistsAProton.at(i));
			}
			
			// (anti)deuteron histograms
			CreateHistosTrack(fHistsDeuteron);
			for (Int_t i=0;i<(int)(fHistsDeuteron.size()); i++){
				fOutputDeuterons->Add(fHistsDeuteron.at(i));
			}
			
			CreateHistosTrack(fHistsADeuteron);
			for (Int_t i=0;i<(int)(fHistsADeuteron.size()); i++){
				fOutputADeuterons->Add(fHistsADeuteron.at(i));
			}
	}

	CreateHistosTrack(fHistsHe3);
	for (Int_t i=0;i<(int)(fHistsHe3.size()); i++){
		fOutputHe3->Add(fHistsHe3.at(i));
	}
	fOutputHe3->Add(fNtupleHe3);
    
	CreateHistosTrack(fHistsAHe3);
	for (Int_t i=0;i<(int)(fHistsAHe3.size()); i++){
		fOutputAHe3->Add(fHistsAHe3.at(i));
	}
	fOutputAHe3->Add(fNtupleAHe3);

    PostData(1, fOutputList);
    // postdata will notify the analysis manager of changes / updates to the
    // fOutputList object. the manager will in the end take care of writing output to file
    // so it needs to know what's in the output
	
}

//________________________________________________________________________
void AliAnalysisTaskHe3::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    
    fVEvent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVEvent) {
        //printf("fVEvent not available\n");
        return;
    }
    
    fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAODEvent) {
        printf("fAODEvent not available\n");
        return;
    }
    
    fEventStat->Fill(kSelectedEvents); // all events after trigger and physics selection
    // only for cross-check: is kINT7 selected
    UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isINT7selected = fSelectMask and AliVEvent::kINT7;
    if (!isINT7selected){
        PostData(1, fOutputList);
        return;
    }
    fEventStat->Fill(kINT7selected);
    
    // only for cross-check: Incomplete events from DAQ
    if (fVEvent->IsIncompleteDAQ()) { // probably not implemented for AOD?
        PostData(1, fOutputList);
        return;
    }
    fEventStat->Fill(kDAQincomplete);
    
    //only for cross-check: VZERO timing desicion
    Int_t fV0ADecision;
    Int_t fV0CDecision;
    AliVVZERO* vzeroData = fVEvent->GetVZEROData();
    fV0ADecision = vzeroData->GetV0ADecision();
    fV0CDecision = vzeroData->GetV0CDecision();
    if (!fV0ADecision || !fV0CDecision) {
        PostData(1, fOutputList);
        return;
    }
    fEventStat->Fill(kV0timing);
    
    //only for cross-check: SPD clusters vs tracklets cut
    AliAnalysisUtils *utils = new AliAnalysisUtils();
    //utils->SetASPDCvsTCut(100.); // 65 by default
    //utils->SetBSPDCvsTCut(5.); // 4 by default
    Bool_t ClustersVsTrackletBG = utils->IsSPDClusterVsTrackletBG(fVEvent);
    if (ClustersVsTrackletBG) {
        PostData(1, fOutputList);
        return;
    }
    fEventStat->Fill(kClusterTrackletCut);
    
    // vertex filter
    const AliVVertex *vertex = fVEvent->GetPrimaryVertex();
    Bool_t lVertexAcceptable = (TMath::Abs(vertex->GetZ()) <= 10.);
    Bool_t lVertexNcontrib   = (vertex->GetNContributors() > 0);
    
    if (!lVertexAcceptable){
        PostData(1, fOutputList);
        return;
    }
    fEventStat->Fill(kVertexZ);
    
    if (!lVertexNcontrib){
        PostData(1, fOutputList);
        return;
    }
    fEventStat->Fill(kVtxNcontrib);
    
    // SPD pile-up in mult. bins
    if (fVEvent->IsPileupFromSPDInMultBins()){
        PostData(1, fOutputList);
        return;
    }
    fEventStat->Fill(kSPDPileUp);
    
    // spd vertex and its resolution
    const AliVVertex* vtxSPD = fVEvent->GetPrimaryVertexSPD();
    Double_t cov[6]={0};
    vtxSPD->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
//    if (TMath::Abs(vtxSPD->GetZ() - vertex->GetZ()) > 0.3){
//        PostData(1, fOutputList);
//        return;
//    }
    fEventStat->Fill(kVtxDisplace);
    
//    if (vtxSPD->IsFromVertexerZ() && zRes > 0.3 ) {
//        PostData(1, fOutputList);
//        return;
//    }
    fEventStat->Fill(kVtxRes);

	AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
	
    
     
    // V0 information
    Int_t MultV0A = 0;
    Int_t MultV0C = 0;
    for(Int_t i=0; i<32; ++i) {
        MultV0A += vzeroData->GetMultiplicityV0A(i);
        MultV0C += vzeroData->GetMultiplicityV0C(i);
    }
    fMultV0->Fill(MultV0A + MultV0C);
    
    // centrality distribution
    Float_t lPercentile = 300;
    AliMultSelection *MultSelection = 0x0;
    MultSelection = (AliMultSelection * ) fVEvent->FindListObject("MultSelection");
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    }else{
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
    }

	int Nch = 0;
    
    // track loop
    for (Int_t iTrack = 0; iTrack < fAODEvent->GetNumberOfTracks(); iTrack++){
        
        AliAODTrack* track = (AliAODTrack*)fAODEvent->GetTrack(iTrack);
        if (!track) { Printf("ERROR: Could not receive track %d", iTrack); continue; }
        // To do: track label for AOD data?
        Int_t label = track->GetLabel();
        Double_t trackP   = track->P();
        Double_t trackEta = track->Eta();
        
        // ===================== track cuts =====================
        // filter bit, pt and eta
		// todo: implement method so that if event doesnt have anything that passes track cuts, it doesnt fill Nch. This should solve the Nch 0 issue. Or increment Nch before track cuts
        if (!(track->TestFilterBit(fFilterBit)))          continue;
        if (trackP < fLowPCut || trackP > fHighPCut)      continue;
        if (TMath::Abs(trackEta) > fEtaCut)               continue;
        if (track->GetITSNcls() < fMinClIts)              continue;


        
        // DCA information
        double DCAxy = -99.;
        double DCAz  = -99.;
        double dcaVals[2] = {-99., -99.};
        double covar[3]={0.,0.,0.};
        AliAODTrack copy(*track);
        const AliVVertex *AODeventVtx = copy.GetAODEvent()->GetPrimaryVertex();
        Double_t AODeventMagneticF = copy.GetAODEvent()->GetMagneticField();
        if (copy.PropagateToDCA(AODeventVtx, AODeventMagneticF,10, dcaVals, covar)){
            DCAxy = dcaVals[0];
            DCAz  = dcaVals[1];
        } else {
            DCAxy = -99.; // track->DCA();
            DCAz  = -99.; // track->ZAtDCA();
        }
        
        if (TMath::Abs(DCAz ) > fMaxDCAz)        continue;
        if (TMath::Abs(DCAxy) > fMaxDCAxyCut)    continue;
		// ===================== Multiplicity counter ==========
		// After track cuts and secondary (DCA) cut

		if (TMath::Abs(track->Charge()) > 0.) { Nch++;} 
		//cout<< "Nch now : "<< Nch<<endl;



        // ===================== PID selection =====================
        // TPC
        Float_t nSigmaTPCprot = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        Float_t nSigmaTPCdeut = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
		Float_t nSigmaTPCHe3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);
        
        // ITS (not used for selection)
        Float_t nSigmaITSprot = fPIDResponse->NumberOfSigmasITS(track, AliPID::kProton);
        Float_t nSigmaITSdeut = fPIDResponse->NumberOfSigmasITS(track, AliPID::kDeuteron);
        
        // TOF
        Bool_t isTOF = (track->GetStatus() & AliVTrack::kTOFout)
                    && (track->GetStatus() & AliVTrack::kTIME);
        
        Float_t nSigmaTOFpion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
        Float_t nSigmaTOFprot = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
        Float_t nSigmaTOFdeut = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);
        
        // TO DO: check TOF matching in data and MC
        // ESD data: GetTOFsignalDz(), GetTOFsignalDx()
        
        Bool_t TOFpid_prot = isTOF;
        Bool_t TOFpid_deut = isTOF;
        
        if (fUseTOFPidCut){
            
            TOFpid_prot = (isTOF && TMath::Abs(nSigmaTOFprot) < fMaxTOFnSigma);
            TOFpid_deut = (isTOF && TMath::Abs(nSigmaTOFdeut) < fMaxTOFnSigma);
        
        }
        
        // protons
        Bool_t isProtonSelected   = (TMath::Abs(nSigmaTPCprot) < fMaxTPCnSigma) &&
                                    (!(trackP > fMomTOFProt && !TOFpid_prot));
        
        // deuterons
        Bool_t isDeuteronSelected = (TMath::Abs(nSigmaTPCdeut) < fMaxTPCnSigma) &&
                                    (!(trackP > fMomTOFDeut && !TOFpid_deut));
        
		Bool_t isHelium3Selected = (nSigmaTPCHe3 > -10.0) and (nSigmaTPCHe3 < 7); //Helium 3 is triggered in TPC only with a relaxed condition of 5, EDIT: cut reduced to 3 sigma due to large includion of background in the 5 sigma cut. See 16qt set from 260919. Edit2: considered relaxing condition to -10<nSigma<5 in order to better fit the TPCnSigma plots, but it is not necesary is the TPC Signal is cut at 100.
		//Edit 3: Condition relaxed to 10 in order to fit larger background in low statistics bins. 

        // ===================== fill track histos =====================
        if (isProtonSelected && kAnalyseAllParticles){
            if (track->Charge() > 0){
                FillHistosTrack(fHistsProton,track);
                
            } else if (track->Charge() < 0){
                FillHistosTrack(fHistsAProton,track);
                
            }
        }
        
        if (isDeuteronSelected && kAnalyseAllParticles){
            if (track->Charge() > 0){
                FillHistosTrack(fHistsDeuteron,track);
                
            } else if (track->Charge() < 0){
                FillHistosTrack(fHistsADeuteron,track);
                
            }
        }

		if (isHelium3Selected){
				Float_t vars[18];
			   	vars[0] = track->P();
			    	vars[1] = track->Pt();
			   	vars[2] = track->GetTPCsignal();
				vars[3] = nSigmaTPCHe3;
				vars[4] = track->Charge();
				vars[5] = GetTOFBeta(track);
				vars[6] = DCAxy;
				vars[7] = DCAz;
				vars[8] = GetMass2TOF(GetTOFBeta(track), track);
				vars[9] = track->GetTPCNcls();
				vars[10] = track->GetITSNcls();
				vars[11] = track->GetTPCsignalN();
				vars[12] = track->Eta();
				vars[13] = fPIDResponse->NumberOfSigmasITS(track, AliPID::kHe3);
				vars[14] = track->Chi2perNDF();
				vars[15] = track->GetITSchi2();
				vars[16] = track->GetTPCClusterInfo(2, 1);
				vars[17] = track->GetTPCmomentum();
				
			if (track->Charge() > 0){
				FillHistosTrack(fHistsHe3, track);
				if (track->GetTPCsignal() >100) {fNtupleHe3->Fill(vars);}
			} else if (track->Charge() < 0) {
				FillHistosTrack(fHistsAHe3,track);
				if (track->GetTPCsignal() > 100) {fNtupleAHe3->Fill(vars);} //The TPC Signal requirement is to exclude the large amount of well seperated background which would make the Ntuple too large to compile & analyse. The histograms do not have this requirement and the seperation of the two signals can be seen in them.
			}
		}
        
    } // track loop
	fNch->Fill(Nch);//Filled here to check how many events were completely excluded by the filterbit cut. Just cut off the zero bin for the correct number. Normalization can come from Vertex z distr.

	if (Nch>0) { //filled here to exclude filterbit bias with normalization
			// fill event histograms
			fVertexZ->Fill(vertex->GetZ());
			fVtxContrib->Fill(vertex->GetNContributors());
			fSPDVtxResol->Fill(zRes);
			fVtxDisplacement->Fill(vtxSPD->GetZ(),vertex->GetZ());
			fCentV0M->Fill(lPercentile);
			fCentV0MZoomed->Fill(lPercentile);		
			fNchHeader->Fill(header->GetRefMultiplicity());
	}
  

    PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHe3::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}

//________________________________________________________________________
void AliAnalysisTaskHe3::CreateHistosTrack(vector<TH1*> &histos)
{
	
	/* This function creates the hitograms which are pushed into the vector and can then be accessed (i.e. filled) */
		/* ####### histos vector histogram order #######
		 *
		 * fHistP					0
		 * fHistPt					1
		 * fHistEtaPhi				2
		 * fHistTPCCrossedRows		3
		 * fHistTPCClusters			4
		 * fHistTPCCRoverFind		5
		 * fHistTPCFracShared		6
		 * fHistTPCSignalN			7
		 * fHistTPCSignalNfrac		8
		 * fHistITSnCls				9
		 * fHistChi2				10
		 * fHistDCAxy				11
		 * fHistDCAz				12
		 * fHistDCAxyDCAz			13
		 * fHistTPCSignal			14
		 * fHistTPCnSigmaProt		15
		 * fHistTPCnSigmaDeut		16
		 * fHistTPCnSigmaHe3		17
		 * fHistTOFSignal			18
		 * fHistTOFnSigmaProt		19
		 * fHistTOFnSigmaDeut		20
		 * fHistTOFnSigmaPion		21
		 * fHistTOFnSigmaHe3		22
		 * fHistTOFmass2_DCAxy_p	23
		 * fHistITSSignal			24
		 * fHistITSnSigmaProt		25
		 * fHistITSnSigmaDeut		26
		 * fHistITSnSigmaHe3		27
		 * fHistTPCSignalChargeCorrected	28
		 * fHistTPCnSigmaHe3ChargeCorrected	29
		 * fHistTOFSignalChargeCorrected	30
		 * fHistTOFnSigmaHe3ChargeCorrected	31
		 * fHistTOFmass2_DCAxy_pChargeCorrected	32
		 * fHistITSSignalChargeCorrected		33
		 * fHistITSnSigmaHe3ChargeCorrected		34
		 *
		 */
		 
    TH1F* fHistP = new TH1F("fHistP","p distribution;#it{p}/Z, GeV/#it{c};counts",100, 0., 10.);
    fHistP->Sumw2();
    histos.push_back(fHistP);
    
    TH1F* fHistPt = new TH1F("fHistPt","pt distribution;#it{p}_{T}, GeV/#it{c};counts",100, 0., 10.);
    fHistPt->Sumw2();
    histos.push_back(fHistPt);

    TH2F* fHistEtaPhi = new TH2F("fHistEtaPhi","eta-phi distribution;#eta;#varphi",100, -1.0, 1.0, 360, 0, 2*TMath::Pi());
    histos.push_back(fHistEtaPhi);
    
    // tracking
    TH1F* fHistTPCCrossedRows = new TH1F("fHistTPCCrossedRows","TPC crossed rows;N TPC crossed rows;counts",200,0,200);
    histos.push_back(fHistTPCCrossedRows);
    
    TH1F* fHistTPCClusters = new TH1F("fHistTPCClusters","TPC clusters;N TPC clusters;counts",200,0,200);
    histos.push_back(fHistTPCClusters);
    
    TH1F* fHistTPCCRoverFind = new TH1F("fHistTPCCRoverFind","TPC crossed rows / findable;TPC crossed rows / findable;counts",150,0,1.5);
    histos.push_back(fHistTPCCRoverFind);
    
    TH1F* fHistTPCFracShared = new TH1F("fHistTPCFracShared","TPC shared clusters;TPC shared clusters;counts",250,-1.0,1.5);
    histos.push_back(fHistTPCFracShared);
    
    TH1F* fHistTPCSignalN = new TH1F("fHistTPCSignalN","Number of PID clusters TPC;N of TPC PID clusters;counts",200,0.,200.);
    histos.push_back(fHistTPCSignalN);
    
    TH1F* fHistTPCSignalNfrac = new TH1F("fHistTPCSignalNfrac","Fraction of PID clusters TPC;Fraction of TPC PID clusters;counts",120,0.,1.2);
    histos.push_back(fHistTPCSignalNfrac);
    
    TH1F* fHistITSnCls = new TH1F("fHistITSnCls","Number of ITS clusters;N of ITS clusters;counts",10,0.,10.);
    histos.push_back(fHistITSnCls);
    
    TH1F* fHistChi2 = new TH1F("fHistChi2","track chi2;#chi^{2};counts",100,0,10);
    histos.push_back(fHistChi2);
    
    TH1F* fHistDCAxy = new TH1F("fHistDCAxy","DCA xy;DCA_{xy};counts",400,-2.0,2.0);
    histos.push_back(fHistDCAxy);
    
    TH1F* fHistDCAz = new TH1F("fHistDCAz","DCA z;DCA_{z};counts",400,-2.0,2.0);
    histos.push_back(fHistDCAz);
    
    TH2F* fHistDCAxyDCAz = new TH2F("fHistDCAxyDCAz","DCAxy vs DCAz;DCA_{xy};DCA_{z}",400,-2.0,2.0, 400, -2.0, 2.0);
    histos.push_back(fHistDCAxyDCAz);
    
    // PID histos
    TH2F* fHistTPCSignal = new TH2F("fHistTPCSignal","TPC dE/dx;#it{p}/Z, GeV/#it{c};TPC dE/dx",500,0.,10.0,1500,0.,1500.);
    histos.push_back(fHistTPCSignal);
    
    TH2F* fHistTPCnSigmaProt = new TH2F("fHistTPCnSigmaProt","TPC nSigma prot;#it{p}/Z, GeV/#it{c};TPCn#sigma_{prot}",100,0.,10.0,100,-5.,5.);
    histos.push_back(fHistTPCnSigmaProt);
    
    TH2F* fHistTPCnSigmaDeut = new TH2F("fHistTPCnSigmaDeut","TPC nSigma deut;#it{p}/Z, GeV/#it{c};TPCn#sigma_{deut}",100,0.,10.0,100,-5.,5.);
    histos.push_back(fHistTPCnSigmaDeut);
    
    TH2F* fHistTPCnSigmaHe3 = new TH2F("fHistTPCnSigmaHe3","TPC nSigma He3;#it{p}/Z, GeV/#it{c};TPCn#sigma_{He3}",100,0.,10.0,100,-5.,5.);
	histos.push_back(fHistTPCnSigmaHe3);

    TH2F* fHistTOFSignal = new TH2F("fHistTOFSignal","TOF signal;#it{p}/Z, GeV/#it{c};TOF signal",500,0.,10.0,200,0.,2.);
    histos.push_back(fHistTOFSignal);
    
    TH2F* fHistTOFnSigmaProt = new TH2F("fHistTOFnSigmaProt","TOF nSigma prot;#it{p}/Z, GeV/#it{c};TOFn#sigma_{prot}",100,0.,10.0,100,-5.,5.);
    histos.push_back(fHistTOFnSigmaProt);
    
    TH2F* fHistTOFnSigmaDeut = new TH2F("fHistTOFnSigmaDeut","TOF nSigma deut;#it{p}/Z, GeV/#it{c};TOFn#sigma_{deut}",100,0.,10.0,100,-5.,5.);
    histos.push_back(fHistTOFnSigmaDeut);
    
    TH2F* fHistTOFnSigmaPion = new TH2F("fHistTOFnSigmaPion","TOF nSigma pion;#it{p}/Z, GeV/#it{c};TOFn#sigma_{pion}",100,0.,10.0,100,-5.,5.);
    histos.push_back(fHistTOFnSigmaPion);
    
    TH2F* fHistTOFnSigmaHe3 = new TH2F("fHistTOFnSigmaHe3","TOF nSigma He3;#it{p}/Z, GeV/#it{c};TOFn#sigma_{He3}",100,0.,10.0,100,-5.,5.);
	histos.push_back(fHistTOFnSigmaHe3);

    TH3F* fHistTOFmass2_DCAxy_p = new TH3F("fHistTOFmass2_DCAxy_p","TOF mass2 vs DCA_{xy} vs #it{p}/Z;#it{p}/Z, GeV/#it{c};TOF mass^{2} / Z^{2}, (GeV/#it{c}^{2})^{2};DCA_{xy}, cm;",50,0.,5.,1000,0.0,10.0,200,-1.0,1.0);
    histos.push_back(fHistTOFmass2_DCAxy_p);
    
    TH2F* fHistITSSignal = new TH2F("fHistITSSignal","ITS signal;#it{p}/Z, GeV/#it{c};ITS signal",500,0.,10.0,750,0.,750.);
    histos.push_back(fHistITSSignal);
    
    TH2F* fHistITSnSigmaProt = new TH2F("fHistITSnSigmaProt","ITS nSigma prot;#it{p}/Z, GeV/#it{c};ITSn#sigma_{prot}",100,0.,10.0,200,-10.,10.);
    histos.push_back(fHistITSnSigmaProt);
    
    TH2F* fHistITSnSigmaDeut = new TH2F("fHistITSnSigmaDeut","ITS nSigma deut;#it{p}/Z, GeV/#it{c};ITSn#sigma_{deut}",100,0.,10.0,200,-10.,10.);
    histos.push_back(fHistITSnSigmaDeut);
    
    TH2F* fHistITSnSigmaHe3 = new TH2F("fHistITSnSigmaHe3","ITS nSigma He3;#it{p}/Z, GeV/#it{c};ITSn#sigma_{He3}",100,0.,10.0,200,-10.,10.);
	histos.push_back(fHistITSnSigmaHe3);

    TH3F* fHistTOFmass2_DCAxy_pt = new TH3F("fHistTOFmass2_DCAxy_pt","TOF mass2 vs DCA_{xy} vs #it{p_T}/Z;#it{p}/Z, GeV/#it{c};TOF mass^{2} / Z^{2}, (GeV/#it{c}^{2})^{2};DCA_{xy}, cm;",50,0.,5.,1000,0.0,10.0,200,-1.0,1.0);
    histos.push_back(fHistTOFmass2_DCAxy_pt);
	
    TH2F* fHistDCAxyVsPt = new TH2F("fHistDCAxyVsPt","DCA xy;DCA_{xy};#it{p_T}",400,-2.0,2.0, 200, 0, 4);
    histos.push_back(fHistDCAxyVsPt);

	/* The charge corrected histograms are not working  correctly, need to fix them

   	TH2F* fHistTPCSignalChargeCorrected = new TH2F("fHistTPCSignalChargeCorrected","TPC dE/dx;#it{p}, GeV/#it{c};TPC dE/dx",500,0.,10.0,1500,0.,1500.);
    histos.push_back(fHistTPCSignalChargeCorrected);

    TH2F* fHistTPCnSigmaHe3ChargeCorrected = new TH2F("fHistTPCnSigmaHe3ChargeCorrected","TPC nSigma He3;#it{p}, GeV/#it{c};TPCn#sigma_{He3}",100,0.,10.0,100,-5.,5.);
	histos.push_back(fHistTPCnSigmaHe3ChargeCorrected);

    TH2F* fHistTOFSignalChargeCorrected = new TH2F("fHistTOFSignalChargeCorrected","TOF signal;#it{p}, GeV/#it{c};TOF signal",500,0.,10.0,200,0.,2.);
    histos.push_back(fHistTOFSignalChargeCorrected);

    TH2F* fHistTOFnSigmaHe3ChargeCorrected = new TH2F("fHistTOFnSigmaHe3ChargeCorrected","TOF nSigma He3;#it{p}/Z, GeV/#it{c};TOFn#sigma_{He3}",100,0.,10.0,100,-5.,5.);
	histos.push_back(fHistTOFnSigmaHe3ChargeCorrected);

    TH3F* fHistTOFmass2_DCAxy_pChargeCorrected = new TH3F("fHistTOFmass2_DCAxy_pChargeCorrected","TOF mass2 vs DCA_{xy} vs #it{p}/Z;#it{p}, GeV/#it{c};TOF mass^{2} / Z^{2}, (GeV/#it{c}^{2})^{2};DCA_{xy}, cm;",50,0.,5.,1000,0.0,10.0,200,-1.0,1.0);
    histos.push_back(fHistTOFmass2_DCAxy_pChargeCorrected);

    TH2F* fHistITSSignalChargeCorrected = new TH2F("fHistITSSignalChargeCorrected","ITS signal;#it{p}, GeV/#it{c};ITS signal",500,0.,10.0,750,0.,750.);
    histos.push_back(fHistITSSignalChargeCorrected);

    TH2F* fHistITSnSigmaHe3ChargeCorrected = new TH2F("fHistITSnSigmaHe3ChargeCorrected","ITS nSigma He3;#it{p}, GeV/#it{c};ITSn#sigma_{He3}",100,0.,10.0,200,-10.,10.);
	histos.push_back(fHistITSnSigmaHe3ChargeCorrected);
	*/
}

//________________________________________________________________________
void AliAnalysisTaskHe3::FillHistosTrack(vector<TH1*> &histos, AliAODTrack *track)
{
    Double_t trackP   = track->P();
	Double_t trackPCC = track->P()*track->Charge();
    Double_t trackPt  = track->Pt();
    Double_t trackEta = track->Eta();
    Double_t trackPhi = track->Phi();
    
    // TOF information
    Float_t fbetaTOF  = GetTOFBeta(track);
    Float_t fmass2TOF = GetMass2TOF(fbetaTOF,track);
    Float_t fmass2TOFCC = GetMass2TOF(fbetaTOF,track) * track->Charge() * track->Charge();
    
    // PID information
    Float_t nSigmaTPCprot = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    Float_t nSigmaITSprot = fPIDResponse->NumberOfSigmasITS(track, AliPID::kProton);
    Float_t nSigmaTOFprot = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
    
    Float_t nSigmaTPCdeut = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    Float_t nSigmaITSdeut = fPIDResponse->NumberOfSigmasITS(track, AliPID::kDeuteron);
    Float_t nSigmaTOFdeut = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);

    Float_t nSigmaTPCHe3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);
    Float_t nSigmaITSHe3 = fPIDResponse->NumberOfSigmasITS(track, AliPID::kHe3);
    Float_t nSigmaTOFHe3 = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kHe3);
    
    Float_t nSigmaTOFpion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);

    // DCA information
    double DCAxy = -99.;
    double DCAz  = -99.;
    double dcaVals[2] = {-99., -99.};
    double covar[3]={0.,0.,0.};
    AliAODTrack copy(*track);
    const AliVVertex *AODeventVtx = copy.GetAODEvent()->GetPrimaryVertex();
    Double_t AODeventMagneticF = copy.GetAODEvent()->GetMagneticField();
    if (copy.PropagateToDCA(AODeventVtx, AODeventMagneticF,10, dcaVals, covar)){
        DCAxy = dcaVals[0];
        DCAz  = dcaVals[1];
    } else {
        DCAxy = -99.; // track->DCA();
        DCAz  = -99.; // track->ZAtDCA();
    }
    
    // if DCAxy > cut for final analysis -> fill only 3d histogram
    if (TMath::Abs(DCAxy) > fMaxDCAxyAna){
        
        ((TH3F*)(histos.at(23)))->Fill(trackP,fmass2TOF,DCAxy);
        return;
        
    } else {
        
        histos.at(0)->Fill(trackP);
        histos.at(1)->Fill(trackPt);
        histos.at(2)->Fill(trackEta,trackPhi);
        histos.at(3)->Fill(track->GetTPCClusterInfo(2, 1));
        histos.at(4)->Fill(track->GetTPCNcls());
        if (!track->GetTPCNclsF()) {
            histos.at(5)->Fill(0);
        } else {
            histos.at(5)->Fill(track->GetTPCClusterInfo(2, 1)/
                               double(track->GetTPCNclsF()));
        }
        if (!track->GetTPCNcls()){
            histos.at(6)->Fill(-1);
        } else {
            histos.at(6)->Fill(track->GetTPCnclsS()/double(track->GetTPCNcls()));
        }
        histos.at(7)->Fill(track->GetTPCsignalN());
        if (!track->GetTPCNcls()){
            histos.at(8)->Fill(-1);
        } else {
            histos.at(8)->Fill(track->GetTPCsignalN()/double(track->GetTPCNcls())); 
        }
        histos.at(9)->Fill(track->GetITSNcls());
        histos.at(10)->Fill(track->Chi2perNDF());
        histos.at(11)->Fill(DCAxy);
        histos.at(12)->Fill(DCAz);
        histos.at(13)->Fill(DCAxy,DCAz);
        
        // fill PID histos
        histos.at(14)->Fill(trackP,track->GetTPCsignal());
        histos.at(15)->Fill(trackP,nSigmaTPCprot);
        histos.at(16)->Fill(trackP,nSigmaTPCdeut);
		histos.at(17)->Fill(trackP,nSigmaTPCHe3);
        histos.at(18)->Fill(trackP,fbetaTOF);
        histos.at(19)->Fill(trackP,nSigmaTOFprot);
        histos.at(20)->Fill(trackP,nSigmaTOFdeut);
        histos.at(21)->Fill(trackP,nSigmaTOFpion);
		histos.at(22)->Fill(trackP,nSigmaTOFHe3);
        ((TH3F*)(histos.at(23)))->Fill(trackP,fmass2TOF,DCAxy);
        histos.at(24)->Fill(trackP,track->GetITSsignal());
		//cout<<nSigmaITSprot<<endl; //FAST & CENT_woSDD have lower ITS requirement 
        histos.at(25)->Fill(trackP,nSigmaITSprot);
        histos.at(26)->Fill(trackP,nSigmaITSdeut);
		histos.at(27)->Fill(trackP,nSigmaITSHe3);
	    ((TH3F*)(histos.at(28)))->Fill(trackPt,fmass2TOF,DCAxy);
	    histos.at(29)->Fill(DCAxy, trackPt);
		/*
		histos.at(28)->Fill(trackPCC, track->GetTPCsignal());
		histos.at(29)->Fill(trackPCC, nSigmaTPCHe3);
		histos.at(30)->Fill(trackPCC, fbetaTOF);
		histos.at(31)->Fill(trackPCC, nSigmaTOFHe3);
		cout << endl<< "### HERE ###" << endl;
		cout<<histos.at(32)->GetName() <<endl;
		((TH3F*)(histos.at(32)))->Fill(trackPCC, fmass2TOFCC, DCAxy);
		histos.at(33)->Fill(trackPCC, track->GetITSsignal());
		histos.at(34)->Fill(trackPCC, nSigmaITSHe3);
		*/
        
    }
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskHe3::GetTOFBeta(AliAODTrack *track)
{
    float beta = -999;
    double integratedTimes[9] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    
    track->GetIntegratedTimes(integratedTimes);
    
    const float c = 2.99792457999999984e-02;
    float p = track->P();
    float l = integratedTimes[0]*c;
    
    float trackT0 = fPIDResponse->GetTOFResponse().GetStartTime(p);
    
    float timeTOF = track->GetTOFsignal() - trackT0;
    if(timeTOF > 0){
        beta  = l/timeTOF/c;
    }
    return beta;
}

//________________________________________________________________________
Float_t AliAnalysisTaskHe3::GetMass2TOF(Float_t beta, AliAODTrack *track)
{
    Float_t p = track->P();
    Float_t mass2sq = -999;
    if(!(beta==0)){
        mass2sq = ((1/(beta*beta))-1)*(p*p);
    }
    return mass2sq;
}
