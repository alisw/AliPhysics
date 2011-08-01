//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis
//  - MC output
//  implementation file
//
//*-- Author: Marcelo G. Munhoz (USP)
//_________________________________________________________________________

#include "AliAnalysisEmEtReconstructed.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "TH2F.h"
#include "TParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "TList.h"
#include "AliESDCaloCluster.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "AliEMCALTrack.h"
#include "AliESDtrackCuts.h"
#include "AliEMCALGeometry.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "TGeoManager.h"

using namespace std;

ClassImp(AliAnalysisEmEtReconstructed);


// ctor
AliAnalysisEmEtReconstructed::AliAnalysisEmEtReconstructed():AliAnalysisEtReconstructed()
,fAllRectotETDep(0)
,fElectronMatchtotETDep(0)
,fNeutralRectotET(0)
,fTotEMRectotET(0)

,fMuonMatchtotETDep(0), fPionMatchtotETDep(0), fKaonMatchtotETDep(0), fProtonMatchtotETDep(0)
,fTotChargedMatchtotETDep(0)

,fTotalRectotETDep(0)

,fESD(0)
,fGeoUt(0)

//,fHistAllRecEtaEDepETDep(0) 
//,fHistAllRecEtaETDep(0) 
,fHistAllRecETDep(0) 
,fHistAllRec(0) 
,fHistAllRectotETDep(0) 

/*
,fHistElectronMatchEtaEDepETDep(0) 
,fHistElectronMatchEtaPtETDep(0) 
,fHistElectronMatchEtaETDep(0) 
,fHistElectronMatchEtaPt(0) 

,fHistElectronRec_ResEDep_ETDep(0) 
,fHistElectronRec_ResPt_ETDep(0) 
,fHistElectronRec_ResEDep(0) 
,fHistElectronRec_ResPt(0) 
*/
,fHistElectronRecETDep(0) 
,fHistElectronRec(0) 
,fHistElectronMatchtotETDep(0) 
,fHistElectronRecdEdxP(0)

/*
,fHistNeutralRec_EtaE_ET(0)  
//,fHistNeutralRec_EtaPt_ET(0)  
,fHistNeutralRec_EtaET(0)  
,fHistNeutralRec_EtaE(0)  
//,fHistNeutralRec_EtaPt(0)  
*/
,fHistNeutralRectotET(0)  

,fHistTotEMRectotET(0)

/*
,fHistMuonMatchEtaEDepETDep(0) 
,fHistMuonMatchEtaPtETDep(0) 
,fHistMuonMatchEtaETDep(0) 
,fHistMuonMatchEtaPt(0) 

,fHistMuonRecResEDepETDep(0) 
,fHistMuonRecResPtETDep(0) 
,fHistMuonRecResEDep(0) 
,fHistMuonRecResPt(0) 
*/
,fHistMuonRecETDep(0) 
,fHistMuonRec(0) 
,fHistMuonMatchtotETDep(0) 
,fHistMuonRecdEdxP(0)

/*
,fHistPionMatchEtaEDepETDep(0) 
,fHistPionMatchEtaPtETDep(0) 
,fHistPionMatchEtaETDep(0) 
,fHistPionMatchEtaPt(0) 

,fHistPionRecResEDepETDep(0) 
,fHistPionRecResPtETDep(0) 
,fHistPionRecResEDep(0) 
,fHistPionRecResPt(0) 
*/
,fHistPionRecETDep(0) 
,fHistPionRec(0) 
,fHistPionMatchtotETDep(0) 
,fHistPionRecdEdxP(0)

/*
,fHistKaonMatchEtaEDepETDep(0) 
,fHistKaonMatchEtaPtETDep(0) 
,fHistKaonMatchEtaETDep(0) 
,fHistKaonMatchEtaPt(0) 

,fHistKaonRecResEDepETDep(0) 
,fHistKaonRecResPtETDep(0) 
,fHistKaonRecResEDep(0) 
,fHistKaonRecResPt(0) 
*/
,fHistKaonRecETDep(0) 
,fHistKaonRec(0) 
,fHistKaonMatchtotETDep(0) 
,fHistKaonRecdEdxP(0)

/*
,fHistProtonMatchEtaEDepETDep(0) 
,fHistProtonMatchEtaPtETDep(0) 
,fHistProtonMatchEtaETDep(0) 
,fHistProtonMatchEtaPt(0) 

,fHistProtonRecResEDepETDep(0) 
,fHistProtonRecResPtETDep(0) 
,fHistProtonRecResEDep(0) 
,fHistProtonRecResPt(0) 
*/
,fHistProtonRecETDep(0) 
,fHistProtonRec(0) 
,fHistProtonMatchtotETDep(0) 
,fHistProtonRecdEdxP(0)

,fHistTotChargedMatchtotETDep(0)

,fHistTotalRectotETDep(0)

,fHistDeltaRZ(0)
{//constructor
	fHistogramNameSuffix = TString("EmcalRec");
	
	fResCut = 0.02;
	//fResCut = fEmcalTrackDistanceCut;
	
	TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));
	//TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));
	TGeoManager::Import("geometry.root");
	//fGeoUt = new AliEMCALGeometry("EMCAL_FIRSTYEAR","EMCAL");
}

// dtor
AliAnalysisEmEtReconstructed::~AliAnalysisEmEtReconstructed() 
{//Destructor
  //Marcelo, are you sure you clean up all memory?
}

Int_t AliAnalysisEmEtReconstructed::AnalyseEvent(AliVEvent* ev)
{ // analyse MC and real event info
	if(!ev){
		Printf("ERROR: Event does not exist");   
		return 0;
	}
	
	fESD = dynamic_cast<AliESDEvent*>(ev);
	
	fGeoUt = new AliEMCALGeometry("EMCAL_FIRSTYEAR","EMCAL");
	//fGeoUt = new AliEMCALGeometry("EMCAL_COMPLETE1","EMCAL");
	fGeoUt->SetMisalMatrix(fESD->GetEMCALMatrix(0),0);
	
	ResetEventValues();
	
	// get all emcal clusters
	TRefArray* caloClusters = new TRefArray();
	fESD->GetEMCALClusters( caloClusters );
	
	Int_t nCluster = caloClusters->GetEntries();
	
	Float_t pos[3] = {0};
	TVector3 caloPos(0,0,0);
	TVector3 trackPos(0,0,0);
	Double_t res=0, maxPid=-99;
	Double_t xCluster[4]={0}, xCharged[7]={0};
		
	// loop the clusters
	for (int iCluster = 0; iCluster < nCluster; iCluster++ ) 
	{		
		// Retrieve calo cluster information
		AliESDCaloCluster* caloCluster = ( AliESDCaloCluster* )caloClusters->At( iCluster );
		Float_t caloE = caloCluster->E();
		caloCluster->GetPosition(pos);		
		caloPos.SetXYZ(pos[0],pos[1],pos[2]);
		
		// look for track that matches calo cluster  
		AliESDtrack *track = FindMatch(caloCluster, res);
		
		// Retrieve track PID
		if (track)
			maxPid = GetTrackPID(track);
		else
			maxPid = -99;
		
		// calculate ET
		Double_t etDep = CalculateTransverseEnergy(caloCluster);
		
		// All clusters
		//fHistAllRecEtaEDepETDep->Fill(caloE,caloPos.Eta(),etDep);
		//fHistAllRecEtaETDep->Fill(etDep,caloPos.Eta());
		
		xCluster[0] = caloE;
		xCluster[1] = caloPos.Eta();
		xCluster[2] = TMath::RadToDeg()*caloPos.Phi();
		xCluster[3] = caloCluster->GetNCells();
		fAllRectotETDep += etDep;		
		
		fHistAllRecETDep->Fill(xCluster,etDep);
		fHistAllRec->Fill(xCluster);

		if (track)
		{
			xCharged[0] = track->Eta();
			xCharged[1] = track->Pt();
		}
		else
		{
			xCharged[0] = -99;
			xCharged[1] = -99;
		}
		xCharged[2] = caloE;
		xCharged[3] = caloPos.Eta();
		xCharged[4] = TMath::RadToDeg()*caloPos.Phi();
		xCharged[5] = caloCluster->GetNCells();
		xCharged[6] = res;
		
		Bool_t isCharged = kFALSE;
		
		if (maxPid == AliPID::kProton)
		{
			/*
			fHistProtonRecResEDepETDep->Fill(caloE,res,etDep);
			fHistProtonRecResPtETDep->Fill(track->Pt(),res,etDep);							
			fHistProtonRecResEDep->Fill(caloE,res,etDep);
			fHistProtonRecResPt->Fill(track->Pt(),res,etDep);							
			*/
			
			fHistProtonRecETDep->Fill(xCharged,etDep);
			fHistProtonRec->Fill(xCharged);

			fHistProtonRecdEdxP->Fill(track->P(),track->GetTPCsignal());
			
			if ((res>0.) && (res<fResCut))
			{
				/*
				fHistProtonMatchEtaEDepETDep->Fill(caloE,track->Eta(),etDep);
				fHistProtonMatchEtaPtETDep->Fill(track->Pt(),track->Eta(),etDep);							
				fHistProtonMatchEtaETDep->Fill(etDep,track->Eta());
				fHistProtonMatchEtaPt->Fill(track->Pt(),track->Eta());							
				 */
				fProtonMatchtotETDep += etDep;

				isCharged = kTRUE;
			}
		}
		else if (maxPid == AliPID::kPion)
		{
			/*
			fHistPionRecResEDepETDep->Fill(caloE,res,etDep);
			fHistPionRecResPtETDep->Fill(track->Pt(),res,etDep);							
			fHistPionRecResEDep->Fill(caloE,res);
			fHistPionRecResPt->Fill(track->Pt(),res);							
			*/
			
			fHistPionRecETDep->Fill(xCharged,etDep);
			fHistPionRec->Fill(xCharged);

			fHistPionRecdEdxP->Fill(track->P(),track->GetTPCsignal());
			
			if ((res>0.) && (res<fResCut))
			{
				/*
				fHistPionMatchEtaEDepETDep->Fill(caloE,track->Eta(),etDep);
				fHistPionMatchEtaPtETDep->Fill(track->Pt(),track->Eta(),etDep);							
				fHistPionMatchEtaETDep->Fill(etDep,track->Eta());
				fHistPionMatchEtaPt->Fill(track->Pt(),track->Eta());							
				*/ 

				fPionMatchtotETDep += etDep;
				isCharged = kTRUE;
			}
		}
		else if (maxPid == AliPID::kKaon)
		{
			/*
			fHistKaonRecResEDepETDep->Fill(caloE,Res,etDep);
			fHistKaonRecResPtETDep->Fill(track->Pt(),res,etDep);							
			fHistKaonRecResEDep->Fill(caloE,res);
			fHistKaonRecResPt->Fill(track->Pt(),res);							
			*/
			
			fHistKaonRecETDep->Fill(xCharged,etDep);
			fHistKaonRec->Fill(xCharged);

			fHistKaonRecdEdxP->Fill(track->P(),track->GetTPCsignal());
			
			if ((res>0.) && (res<fResCut))
			{
				/*
				fHistKaonMatchEtaEDepETDep->Fill(caloE,track->Eta(),etDep);
				fHistKaonMatchEtaPtETDep->Fill(track->Pt(),track->Eta(),etDep);							
				fHistKaonMatchEtaETDep->Fill(etDep,track->Eta());
				fHistKaonMatchEtaPt->Fill(track->Pt(),track->Eta());							
				*/

				fKaonMatchtotETDep += etDep;
				isCharged = kTRUE;
			}
		}
		else if (maxPid == AliPID::kMuon)
		{
			/*
			fHistMuonRecResEDepETDep->Fill(caloE,res,etDep);
			fHistMuonRecResPtETDep->Fill(track->Pt(),res,etDep);	
			fHistMuonRecResEDep->Fill(caloE,res);
			fHistMuonRecResPt->Fill(track->Pt(),res);	
			*/
			
			fHistMuonRecETDep->Fill(xCharged,etDep);
			fHistMuonRec->Fill(xCharged);

			fHistMuonRecdEdxP->Fill(track->P(),track->GetTPCsignal());
			
			if ((res>0.) && (res<fResCut))
			{
				/*
				fHistMuonMatchEtaEDepETDep->Fill(caloE,track->Eta(),etDep);
				fHistMuonMatchEtaPtETDep->Fill(track->Pt(),track->Eta(),etDep);							
				fHistMuonMatchEtaETDep->Fill(etDep,track->Eta());
				fHistMuonMatchEtaPt->Fill(track->Pt(),track->Eta());							
				*/

				fMuonMatchtotETDep += etDep;						
				isCharged = kTRUE;
			}
		}
		else if (maxPid == AliPID::kElectron)
		{
			/*
			fHistElectronRec_ResEDep_ETDep->Fill(caloE,res,etDep);
			fHistElectronRec_ResPt_ETDep->Fill(track->Pt(),res,etDep);							
			fHistElectronRec_ResEDep->Fill(caloE,res);
			fHistElectronRec_ResPt->Fill(track->Pt(),res);							
			*/
			
			fHistElectronRecETDep->Fill(xCharged,etDep);
			fHistElectronRec->Fill(xCharged);

			fHistElectronRecdEdxP->Fill(track->P(),track->GetTPCsignal());
			
			if ((res>0.) && (res<fResCut))
			{
				/*
				fHistElectronMatchEtaEDepETDep->Fill(caloE,track->Eta(),etDep);
				fHistElectronMatchEtaPtETDep->Fill(track->Pt(),track->Eta(),etDep);							
				fHistElectronMatchEtaETDep->Fill(etDep,track->Eta());
				fHistElectronMatchEtaPt->Fill(track->Pt(),track->Eta());							
				*/

				fElectronMatchtotETDep += etDep;
				isCharged = kTRUE;
			}
		}
		
		if (!isCharged)
		{
			/*
			 fHistNeutralRec_EtaE_ET->Fill(caloE,caloPos.Eta(),etDep);
			//fHistNeutralRec_EtaPt_ET->Fill(caloPos.Pt(),caloPos.Eta(),etDep);							
			fHistNeutralRec_EtaET->Fill(etDep,caloPos.Eta());
			fHistNeutralRec_EtaE->Fill(caloE,caloPos.Eta());
			//fHistNeutralRec_EtaPt->Fill(caloPos.Pt(),caloPos.Eta());	
			 */
			fNeutralRectotET += etDep;			
		}
		
	} // end of loop over clusters	
	
	fTotEMRectotET = fElectronMatchtotETDep + fNeutralRectotET;
	fTotChargedMatchtotETDep = fMuonMatchtotETDep + fPionMatchtotETDep + fKaonMatchtotETDep + fProtonMatchtotETDep;
	fTotalRectotETDep = fTotEMRectotET + fTotChargedMatchtotETDep;
	
	fHistAllRectotETDep->Fill(fAllRectotETDep);
	
	fHistElectronMatchtotETDep->Fill(fElectronMatchtotETDep); 
	fHistNeutralRectotET->Fill(fNeutralRectotET);
	
	fHistTotEMRectotET->Fill(fTotEMRectotET);
	
	fHistMuonMatchtotETDep->Fill(fMuonMatchtotETDep); 
	fHistPionMatchtotETDep->Fill(fPionMatchtotETDep); 
	fHistKaonMatchtotETDep->Fill(fKaonMatchtotETDep); 
	fHistProtonMatchtotETDep->Fill(fProtonMatchtotETDep); 
	fHistTotChargedMatchtotETDep->Fill(fTotChargedMatchtotETDep);
	
	fHistTotalRectotETDep->Fill(fTotalRectotETDep);
	
	delete fGeoUt;
	delete caloClusters;
	
	return 0;    
}

void AliAnalysisEmEtReconstructed::Init()
{ // init
    AliAnalysisEt::Init();
}


void AliAnalysisEmEtReconstructed::ResetEventValues()
{ // reset event values
	AliAnalysisEt::ResetEventValues();
	
	// collision geometry defaults for p+p:
	fAllRectotETDep = 0;
	
	fElectronMatchtotETDep = 0;
	fNeutralRectotET = 0;
	
	fTotEMRectotET = 0;
	
	fMuonMatchtotETDep = 0; fPionMatchtotETDep = 0; fKaonMatchtotETDep = 0; fProtonMatchtotETDep = 0;
	fTotChargedMatchtotETDep = 0;
	
	fTotalRectotETDep = 0;
}


void AliAnalysisEmEtReconstructed::CreateHistograms()
{ // histogram related additions
	//AliAnalysisEt::CreateHistograms();
	
	//fHistAllRecEtaEDepETDep = CreateEtaEHisto2D("fHistAllRecEtaEDepETDep_","MC E_{T}, all particles","E_{T}(GeV)");
	//fHistAllRecEtaETDep = CreateEtaEtHisto2D("fHistAllRecEtaETDep_","MC all particles","#");
	
	fHistAllRecETDep = CreateClusterHistoSparse("fHistAllRecETDep_","E_{T}, all particles");
	fHistAllRec = CreateClusterHistoSparse("fHistAllRec_","counts, all particles");
	TString histname = "fHistAllRectotETDep_" + fHistogramNameSuffix;
	fHistAllRectotETDep = new TH1F(histname.Data(),"total ET, all particles",fgNumOfEBins, fgEAxis);

	/*
	fHistElectronMatchEtaEDepETDep = CreateEtaEHisto2D("fHistElectronMatchEtaEDepETDep_","MC E_{T}, primary Electrons, tracking matched","E_{T} dep (GeV)");
	fHistElectronMatchEtaPtETDep = CreateEtaPtHisto2D("fHistElectronMatchEtaPtETDep_","MC E_{T}, primary Electrons","E_{T} dep(GeV)");
	fHistElectronMatchEtaETDep = CreateEtaEtHisto2D("fHistElectronMatchEtaETDep_","MC primary Electrons","#");
	fHistElectronMatchEtaPt = CreateEtaPtHisto2D("fHistElectronMatchEtaPt_","MC E_{T}, primary Electrons","#");
	 */
	fHistElectronRecETDep = CreateChargedPartHistoSparse("fHistElectronRecETDep_","E_{T}, electrons");
	fHistElectronRec = CreateChargedPartHistoSparse("fHistElectronRec_","counts, electrons");
	histname = "fHistElectronMatchtotETDep_" + fHistogramNameSuffix;
	fHistElectronMatchtotETDep = new TH1F(histname.Data(),"total ET, MC primary Electrons",fgNumOfEBins, fgEAxis);
	/*
	fHistElectronRec_ResEDep_ETDep = CreateResEHisto2D("fHistElectronRec_ResEDep_ETDep_","MC E_{T}, primary Electrons","E_{T} dep (GeV)");
	fHistElectronRec_ResPt_ETDep = CreateResPtHisto2D("fHistElectronRec_ResPt_ETDep_","MC E_{T}, primary Electrons","E_{T} dep (GeV)");
	fHistElectronRec_ResEDep = CreateResEHisto2D("fHistElectronRec_ResEDep_","MC primary Electrons","#");
	fHistElectronRec_ResPt  = CreateResPtHisto2D("fHistElectronRec_ResPt_","MC primary Electrons","#");	
	 */
	histname = "fHistElectronRecdEdxP_" + fHistogramNameSuffix;
	fHistElectronRecdEdxP = new TH2F(histname,"TPC dEdx vs P",100,0.,10.,100,0.,200.);

	/*
	fHistNeutralRec_EtaE_ET = CreateEtaEHisto2D("fHistNeutralRec_EtaE_ET_","MC E_{T}, primary Neutrals","E_{T}(GeV)"); 
	//fHistNeutralRec_EtaPt_ET = CreateEtaPtHisto2D("fHistNeutralRec_EtaPt_ET_","MC E_{T}, primary Neutrals","E_{T}(GeV)"); 
	fHistNeutralRec_EtaET = CreateEtaEtHisto2D("fHistNeutralRec_EtaET_","MC primary Neutrals","#"); 
	fHistNeutralRec_EtaE = CreateEtaEHisto2D("fHistNeutralRec_EtaE_","MC primary Neutrals","#"); 
	//fHistNeutralRec_EtaPt = CreateEtaPtHisto2D("fHistNeutralRec_EtaPt_","MC primary Neutrals","#"); 
	*/
	histname = "fHistNeutralRectotET_" + fHistogramNameSuffix;
	fHistNeutralRectotET = new TH1F(histname.Data(),"total ET, neutral particles",fgNumOfEBins, fgEAxis);
	
	histname = "fHistTotEMRectotET_" + fHistogramNameSuffix;
	fHistTotEMRectotET = new TH1F(histname.Data(),"total electromagnetic ET",fgNumOfEBins, fgEAxis);

	/*
	fHistMuonMatchEtaEDepETDep = CreateEtaEHisto2D("fHistMuonMatchEtaEDepETDep_","MC E_{T}, primary Muons, tracking matched","E_{T} dep (GeV)");
	fHistMuonMatchEtaPtETDep = CreateEtaPtHisto2D("fHistMuonMatchEtaPtETDep_","MC E_{T}, primary Muons","E_{T} dep(GeV)");
	fHistMuonMatchEtaETDep = CreateEtaEtHisto2D("fHistMuonMatchEtaETDep_","MC primary Muons","#");
	fHistMuonMatchEtaPt = CreateEtaPtHisto2D("fHistMuonMatchEtaPt_","MC E_{T}, primary Muons","#");
	*/
	fHistMuonRecETDep = CreateChargedPartHistoSparse("fHistMuonRecETDep_","E_{T}, muons");
	fHistMuonRec = CreateChargedPartHistoSparse("fHistMuonRec_","counts, muons");
	histname = "fHistMuonMatchtotETDep_" + fHistogramNameSuffix;
	fHistMuonMatchtotETDep = new TH1F(histname.Data(),"total ET, Muons",fgNumOfEBins, fgEAxis);
	/*
	fHistMuonRecResEDepETDep = CreateResEHisto2D("fHistMuonRecResEDepETDep_","MC E_{T}, primary Muons","E_{T} dep (GeV)");
	fHistMuonRecResPtETDep = CreateResPtHisto2D("fHistMuonRecResPtETDep_","MC E_{T}, primary Muons","E_{T} dep (GeV)");
	fHistMuonRecResEDep = CreateResEHisto2D("fHistMuonRecResEDep_","MC primary Muons","#");
	fHistMuonRecResPt  = CreateResPtHisto2D("fHistMuonRecResPt_","MC primary Muons","#");
	*/
	histname = "fHistMuonRecdEdxP_" + fHistogramNameSuffix;
	fHistMuonRecdEdxP = new TH2F(histname,"TPC dEdx vs P",100,0.,10.,100,0.,200.);
	
	/*
	fHistPionMatchEtaEDepETDep = CreateEtaEHisto2D("fHistPionMatchEtaEDepETDep_","MC E_{T}, primary Pions, tracking matched","E_{T} dep (GeV)");
	fHistPionMatchEtaPtETDep = CreateEtaPtHisto2D("fHistPionMatchEtaPtETDep_","MC E_{T}, primary Pions","E_{T} dep(GeV)");
	fHistPionMatchEtaETDep = CreateEtaEtHisto2D("fHistPionMatchEtaETDep_","MC primary Pions","#");
	fHistPionMatchEtaPt = CreateEtaPtHisto2D("fHistPionMatchEtaPt_","MC E_{T}, primary Pions","#");
	*/
	fHistPionRecETDep = CreateChargedPartHistoSparse("fHistPionRecETDep_","E_{T}, pions");
	fHistPionRec = CreateChargedPartHistoSparse("fHistPionRec_","counts, pions");
	histname = "fHistPionMatchtotETDep_" + fHistogramNameSuffix;
	fHistPionMatchtotETDep = new TH1F(histname.Data(),"total ET, Pions",fgNumOfEBins, fgEAxis);
	/*
	fHistPionRecResEDepETDep = CreateResEHisto2D("fHistPionRecResEDepETDep_","MC E_{T}, primary Pions","E_{T} dep (GeV)");
	fHistPionRecResPtETDep = CreateResPtHisto2D("fHistPionRecResPtETDep_","MC E_{T}, primary Pions","E_{T} dep (GeV)");
	fHistPionRecResEDep = CreateResEHisto2D("fHistPionRecResEDep_","MC primary Pions","#");
	fHistPionRecResPt  = CreateResPtHisto2D("fHistPionRecResPt_","MC primary Pions","#");
	*/
	histname = "fHistPionRecdEdxP_" + fHistogramNameSuffix;
	fHistPionRecdEdxP = new TH2F(histname,"TPC dEdx vs P",100,0.,10.,100,0.,200.);
	
	/*
	fHistKaonMatchEtaEDepETDep = CreateEtaEHisto2D("fHistKaonMatchEtaEDepETDep_","MC E_{T}, primary Kaons, tracking matched","E_{T} dep (GeV)");
	fHistKaonMatchEtaPtETDep = CreateEtaPtHisto2D("fHistKaonMatchEtaPtETDep_","MC E_{T}, primary Kaons","E_{T} dep(GeV)");
	fHistKaonMatchEtaETDep = CreateEtaEtHisto2D("fHistKaonMatchEtaETDep_","MC primary Kaons","#");
	fHistKaonMatchEtaPt = CreateEtaPtHisto2D("fHistKaonMatchEtaPt_","MC primary Kaons","#");
	*/
	fHistKaonRecETDep = CreateChargedPartHistoSparse("fHistKaonRecETDep_","E_{T}, kaons");
	fHistKaonRec = CreateChargedPartHistoSparse("fHistKaonRec_","counts, kaons");
	histname = "fHistKaonMatchtotETDep_" + fHistogramNameSuffix;
	fHistKaonMatchtotETDep = new TH1F(histname.Data(),"total ET, Kaons",fgNumOfEBins, fgEAxis);
	/*
	fHistKaonRecResEDepETDep = CreateResEHisto2D("fHistKaonRecResEDepETDep_","MC E_{T}, primary Kaons","E_{T} dep (GeV)");
	fHistKaonRecResPtETDep = CreateResPtHisto2D("fHistKaonRecResPtETDep_","MC E_{T}, primary Kaons","E_{T} dep (GeV)");
	fHistKaonRecResEDep = CreateResEHisto2D("fHistKaonRecResEDep_","MC primary Kaons","#");
	fHistKaonRecResPt  = CreateResPtHisto2D("fHistKaonRecResPt_","MC primary Kaons","#");	
	*/
	histname = "fHistKaonRecdEdxP_" + fHistogramNameSuffix;
	fHistKaonRecdEdxP = new TH2F(histname,"TPC dEdx vs P",100,0.,10.,100,0.,200.);
	
	/*
	fHistProtonMatchEtaEDepETDep = CreateEtaEHisto2D("fHistProtonMatchEtaEDepETDep_","MC E_{T}, primary Protons, tracking matched","E_{T} dep (GeV)");
	fHistProtonMatchEtaPtETDep = CreateEtaPtHisto2D("fHistProtonMatchEtaPtETDep_","MC E_{T}, primary Protons","E_{T} dep(GeV)");
	fHistProtonMatchEtaETDep = CreateEtaEtHisto2D("fHistProtonMatchEtaETDep_","MC primary Protons","#");
	fHistProtonMatchEtaPt = CreateEtaPtHisto2D("fHistProtonMatchEtaPt_","MC primary Protons","#");
	*/
	fHistProtonRecETDep = CreateChargedPartHistoSparse("fHistProtonRecETDep_","E_{T}, protons");
	fHistProtonRec = CreateChargedPartHistoSparse("fHistProtonRec_","counts, protons");
	histname = "fHistProtonMatchtotETDep_" + fHistogramNameSuffix;
	fHistProtonMatchtotETDep = new TH1F(histname.Data(),"total ET, Protons",fgNumOfEBins, fgEAxis);
	/*
	fHistProtonRecResEDepETDep = CreateResEHisto2D("fHistProtonRecResEDepETDep_","MC E_{T}, primary Protons","E_{T} dep (GeV)");
	fHistProtonRecResPtETDep = CreateResPtHisto2D("fHistProtonRecResPtETDep_","MC E_{T}, primary Protons","E_{T} dep (GeV)");
	fHistProtonRecResEDep = CreateResEHisto2D("fHistProtonRecResEDep_","MC primary Protons","#");
	fHistProtonRecResPt  = CreateResPtHisto2D("fHistProtonRecResPt_","MC primary Protons","#");
	*/
	histname = "fHistProtonRecdEdxP_" + fHistogramNameSuffix;
	fHistProtonRecdEdxP = new TH2F(histname,"TPC dEdx vs P",100,0.,10.,100,0.,200.);
	
	histname = "fHistTotChargedMatchtotETDep_" + fHistogramNameSuffix;
	fHistTotChargedMatchtotETDep = new TH1F(histname.Data(),"total ET, charged particles",fgNumOfEBins, fgEAxis);
	
	histname = "fHistTotalRectotETDep_" + fHistogramNameSuffix;
	fHistTotalRectotETDep = new TH1F(histname.Data(),"total ET, all particles",fgNumOfEBins, fgEAxis);		
	
	histname = "fHistDeltaRZ_" + fHistogramNameSuffix;
	fHistDeltaRZ = new TH2F(histname,"#Delta#phi vs #Delta#eta (track projection - cluster position)",200,-0.1,0.1,200,-0.1,0.1);
}

void AliAnalysisEmEtReconstructed::FillOutputList(TList *list)
{//Function for filling the output list
	//AliAnalysisEt::FillOutputList(list);
	
	//list->Add(fHistAllRecEtaEDepETDep); 
	//list->Add(fHistAllRecEtaETDep); 
	list->Add(fHistAllRecETDep); 
	list->Add(fHistAllRec); 
	list->Add(fHistAllRectotETDep); 

	/*
	list->Add(fHistElectronMatchEtaEDepETDep); 
	list->Add(fHistElectronMatchEtaPtETDep); 
	list->Add(fHistElectronMatchEtaETDep); 
	list->Add(fHistElectronMatchEtaPt); 
	
	list->Add(fHistElectronRec_ResEDep_ETDep); 
	list->Add(fHistElectronRec_ResPt_ETDep); 
	list->Add(fHistElectronRec_ResEDep); 
	list->Add(fHistElectronRec_ResPt); 
	*/
	
	list->Add(fHistElectronRecETDep); 
	list->Add(fHistElectronRec); 
	list->Add(fHistElectronMatchtotETDep); 
	list->Add(fHistElectronRecdEdxP);
	
	/*
	list->Add(fHistNeutralRec_EtaE_ET);  
	//list->Add(fHistNeutralRec_EtaPt_ET);  
	list->Add(fHistNeutralRec_EtaET);  
	list->Add(fHistNeutralRec_EtaE);  
	//list->Add(fHistNeutralRec_EtaPt);  
	*/
	list->Add(fHistNeutralRectotET);  
	
	list->Add(fHistTotEMRectotET); 
	
	/*
	list->Add(fHistMuonMatchEtaEDepETDep); 
	list->Add(fHistMuonMatchEtaPtETDep); 
	list->Add(fHistMuonMatchEtaETDep); 
	list->Add(fHistMuonMatchEtaPt); 
	
	list->Add(fHistMuonRecResEDepETDep); 
	list->Add(fHistMuonRecResPtETDep); 
	list->Add(fHistMuonRecResEDep); 
	list->Add(fHistMuonRecResPt); 
	 */
	list->Add(fHistMuonRecETDep); 
	list->Add(fHistMuonRec); 
	list->Add(fHistMuonRecdEdxP);
	list->Add(fHistMuonMatchtotETDep); 

	/*
	list->Add(fHistPionMatchEtaEDepETDep); 
	list->Add(fHistPionMatchEtaPtETDep); 
	list->Add(fHistPionMatchEtaETDep); 
	list->Add(fHistPionMatchEtaPt); 
	
	list->Add(fHistPionRecResEDepETDep); 
	list->Add(fHistPionRecResPtETDep); 
	list->Add(fHistPionRecResEDep); 
	list->Add(fHistPionRecResPt); 
	*/
	list->Add(fHistPionRecETDep); 
	list->Add(fHistPionRec); 
	list->Add(fHistPionMatchtotETDep); 
	list->Add(fHistPionRecdEdxP);
	
	/*
	list->Add(fHistKaonMatchEtaEDepETDep); 
	list->Add(fHistKaonMatchEtaPtETDep); 
	list->Add(fHistKaonMatchEtaETDep); 
	list->Add(fHistKaonMatchEtaPt); 
	
	list->Add(fHistKaonRecResEDepETDep); 
	list->Add(fHistKaonRecResPtETDep); 
	list->Add(fHistKaonRecResEDep); 
	list->Add(fHistKaonRecResPt); 
	*/
	
	list->Add(fHistKaonRecETDep); 
	list->Add(fHistKaonRec); 
	list->Add(fHistKaonMatchtotETDep); 
	list->Add(fHistKaonRecdEdxP);
	
	/*
	list->Add(fHistProtonMatchEtaEDepETDep); 
	list->Add(fHistProtonMatchEtaPtETDep); 
	list->Add(fHistProtonMatchEtaETDep); 
	list->Add(fHistProtonMatchEtaPt); 
	
	list->Add(fHistProtonRecResEDepETDep); 
	list->Add(fHistProtonRecResPtETDep); 
	list->Add(fHistProtonRecResEDep); 
	list->Add(fHistProtonRecResPt); 
	*/

	list->Add(fHistProtonRecETDep); 
	list->Add(fHistProtonRec); 
	list->Add(fHistProtonMatchtotETDep); 
	list->Add(fHistProtonRecdEdxP);
	
	list->Add(fHistTotChargedMatchtotETDep); 
	list->Add(fHistTotalRectotETDep); 
	
	list->Add(fHistDeltaRZ);
}

//________________________________________________________________________
// project to a EMCal radius
Bool_t AliAnalysisEmEtReconstructed::GetTrackProjection(AliExternalTrackParam *trackParam, TVector3 &trackPos)
{//Gets the projection of the track
    Bool_t proj = kFALSE;
    Double_t emcalR = fGeoUt->GetEMCGeometry()->GetIPDistance();
	
    if (trackParam) //it is constructed from TParticle
    {
        Double_t trkPos[3] = {0};
		
        //Assume the track is a pion with mass 0.139GeV/c^2
        //Extrapolation step is 1cm
        if(!AliTrackerBase::PropagateTrackToBxByBz(trackParam, emcalR, 0.139, 1, kTRUE, 0.8) ) return proj;
		
        trackParam->GetXYZ(trkPos);
		
        trackPos.SetXYZ(trkPos[0],trkPos[1],trkPos[2]);
		
        proj = kTRUE;               
    }
	
    return proj;
}

//________________________________________________________________________
// project to a cluster position
Bool_t AliAnalysisEmEtReconstructed::GetTrackProjection(AliEMCALTrack* emcTrack, TVector3 &trackPos, TVector3 clusPos)
{//project to a cluster position
	Bool_t proj = kFALSE;
	
	if (emcTrack)
	{	
		Double_t trkPos[3] = {0};
		
		emcTrack->PropagateToGlobal(clusPos.X(),clusPos.Y(),clusPos.Z(),0.,0.);
		emcTrack->GetXYZ(trkPos);
		
		trackPos.SetXYZ(trkPos[0],trkPos[1],trkPos[2]);
		
		proj = kTRUE;
	}
	
	return proj;
}

//________________________________________________________________________	
AliESDtrack* AliAnalysisEmEtReconstructed::FindMatch(const AliESDCaloCluster *caloCluster, Double_t& resMin)
{//find a matched track
	Double_t res=0;
	resMin=999;
	
	TVector3 caloPos(0,0,0);
	Float_t pos[3] = {0};
	caloCluster->GetPosition(pos);		
	caloPos.SetXYZ(pos[0],pos[1],pos[2]);
	
	// loop over tracks
	TVector3 trackPos(0,0,0);
	TVector3 trackMatchPos(0,0,0);
	AliEMCALTrack *emcTrack = 0;	
	AliESDtrack *trackMatch = 0;	
	
	TObjArray* list = fEsdtrackCutsITSTPC->GetAcceptedTracks(fESD);;
	Int_t nGoodTracks = list->GetEntries();
	
	for (Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++) 
	{
		AliESDtrack *track = dynamic_cast<AliESDtrack*> (list->At(iTrack));
        if (!track)
        {
            AliError(Form("ERROR: Could not get track %d", iTrack));
            continue;
        }
		
		emcTrack = new AliEMCALTrack(*track);
		
		if (GetTrackProjection(emcTrack,trackPos,caloPos))
		{
			res = sqrt(pow(trackPos.Phi()-caloPos.Phi(),2)+pow(trackPos.Eta()-caloPos.Eta(),2));
			
			if (res < resMin)
			{
				resMin = res;
				trackMatch = track;
				trackMatchPos.SetXYZ(trackPos.X(),trackPos.Y(),trackPos.Z());
			}
		}
		
		delete emcTrack;
	}		
	
	fHistDeltaRZ->Fill(trackMatchPos.Phi()-caloPos.Phi(),trackMatchPos.Eta()-caloPos.Eta());

	return trackMatch;
}

//________________________________________________________________________	
Double_t AliAnalysisEmEtReconstructed::GetTrackPID(const AliESDtrack *track) const
{//Get the default track ID
	const Double_t *pidWeights = track->PID();
	Int_t maxpid = -1;
	Double_t maxpidweight = 0;
	
	if (pidWeights)
	{
		for (Int_t p =0; p < AliPID::kSPECIES; p++)
		{
			if (pidWeights[p] > maxpidweight)
			{
				maxpidweight = pidWeights[p];
				maxpid = p;
			}
		}
	}
	
	return maxpid;
}
