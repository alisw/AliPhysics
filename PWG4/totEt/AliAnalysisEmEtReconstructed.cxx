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

							    ,fHistAllRecETDep(0) 
							    ,fHistAllRec(0) 
							    ,fHistAllRectotETDep(0) 

							    ,fHistElectronRecETDep(0) 
							    ,fHistElectronRec(0) 
							    ,fHistElectronMatchtotETDep(0) 
							    ,fHistElectronRecdEdxP(0)

							    ,fHistNeutralRectotET(0)  

							    ,fHistTotEMRectotET(0)

							    ,fHistMuonRecETDep(0) 
							    ,fHistMuonRec(0) 
							    ,fHistMuonMatchtotETDep(0) 
							    ,fHistMuonRecdEdxP(0)

							    ,fHistPionRecETDep(0) 
							    ,fHistPionRec(0) 
							    ,fHistPionMatchtotETDep(0) 
							    ,fHistPionRecdEdxP(0)

							    ,fHistKaonRecETDep(0) 
							    ,fHistKaonRec(0) 
							    ,fHistKaonMatchtotETDep(0) 
							    ,fHistKaonRecdEdxP(0)

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
}

// dtor
AliAnalysisEmEtReconstructed::~AliAnalysisEmEtReconstructed() 
{//Destructor 
  delete fGeoUt;

  delete fHistAllRecETDep;
  delete fHistAllRec;
  delete fHistAllRectotETDep;
	
  delete fHistElectronRecETDep;
  delete fHistElectronRec;
  delete fHistElectronMatchtotETDep; 
	
  delete fHistElectronRecdEdxP;

  delete fHistNeutralRectotET;  

  delete fHistTotEMRectotET;

  delete fHistMuonRecETDep;
  delete fHistMuonRec;
  delete fHistMuonMatchtotETDep; 

  delete fHistMuonRecdEdxP;
	
  delete fHistPionRecETDep;
  delete fHistPionRec;
  delete fHistPionMatchtotETDep; 

  delete fHistPionRecdEdxP;

  delete fHistKaonRecETDep;
  delete fHistKaonRec;
  delete fHistKaonMatchtotETDep; 

  delete fHistKaonRecdEdxP;
	
  delete fHistProtonRecETDep;
  delete fHistProtonRec;
  delete fHistProtonMatchtotETDep; 

  delete fHistProtonRecdEdxP;
	
  delete fHistTotChargedMatchtotETDep;
	
  delete fHistTotalRectotETDep;
	
  //few checks
  delete fHistDeltaRZ;
	
}

Int_t AliAnalysisEmEtReconstructed::AnalyseEvent(AliVEvent* ev)
{ // analyse MC and real event info
  if(!ev){
    AliError("ERROR: Event does not exist");   
    return 0;
  }
	
  fESD = dynamic_cast<AliESDEvent*>(ev);

  if(!fGeoUt){
    fGeoUt = AliEMCALGeometry::GetInstance("EMCAL_FIRSTYEARV1");//new AliEMCALGeometry("EMCAL_FIRSTYEAR","EMCAL");
    AliInfo("Creating new AliEMCALGeometry");
  }
  //fGeoUt = new AliEMCALGeometry("EMCAL_COMPLETE1","EMCAL");
  if(!fGeoUt){
    AliInfo("No fGeoUt!");
  }
  else{
    if(!fESD->GetEMCALMatrix(0)){
      AliInfo("No matrix!");
    }
    else{
      fGeoUt->SetMisalMatrix(fESD->GetEMCALMatrix(0),0);
    }
  }
	
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
		
	
      if(fMakeSparse){
	fHistAllRecETDep->Fill(xCluster,etDep);
	fHistAllRec->Fill(xCluster);
      }

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
			
	  if(fMakeSparse){
	    fHistProtonRecETDep->Fill(xCharged,etDep);
	    fHistProtonRec->Fill(xCharged);
	  }

	  fHistProtonRecdEdxP->Fill(track->P(),track->GetTPCsignal());
			
	  if ((res>0.) && (res<fResCut))
	    {
	      fProtonMatchtotETDep += etDep;

	      isCharged = kTRUE;
	    }
	}
      else if (maxPid == AliPID::kPion)
	{
			
	  if(fMakeSparse){
	    fHistPionRecETDep->Fill(xCharged,etDep);
	    fHistPionRec->Fill(xCharged);
	  }

	  fHistPionRecdEdxP->Fill(track->P(),track->GetTPCsignal());
			
	  if ((res>0.) && (res<fResCut))
	    {
	      fPionMatchtotETDep += etDep;
	      isCharged = kTRUE;
	    }
	}
      else if (maxPid == AliPID::kKaon)
	{
			
	  if(fMakeSparse){
	    fHistKaonRecETDep->Fill(xCharged,etDep);
	    fHistKaonRec->Fill(xCharged);
	  }

	  fHistKaonRecdEdxP->Fill(track->P(),track->GetTPCsignal());
			
	  if ((res>0.) && (res<fResCut))
	    {

	      fKaonMatchtotETDep += etDep;
	      isCharged = kTRUE;
	    }
	}
      else if (maxPid == AliPID::kMuon)
	{
			
	  if(fMakeSparse){
	    fHistMuonRecETDep->Fill(xCharged,etDep);
	    fHistMuonRec->Fill(xCharged);
	  }
	  fHistMuonRecdEdxP->Fill(track->P(),track->GetTPCsignal());
			
	  if ((res>0.) && (res<fResCut))
	    {

	      fMuonMatchtotETDep += etDep;						
	      isCharged = kTRUE;
	    }
	}
      else if (maxPid == AliPID::kElectron)
	{
			
	  if(fMakeSparse){
	    fHistElectronRecETDep->Fill(xCharged,etDep);
	    fHistElectronRec->Fill(xCharged);
	  }

	  fHistElectronRecdEdxP->Fill(track->P(),track->GetTPCsignal());
			
	  if ((res>0.) && (res<fResCut))
	    {
	      fElectronMatchtotETDep += etDep;
	      isCharged = kTRUE;
	    }
	}
		
      if (!isCharged)
	{
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
	
  if(fMakeSparse){
    fHistAllRecETDep = CreateClusterHistoSparse("fHistAllRecETDep_","E_{T}, all particles");
    fHistAllRec = CreateClusterHistoSparse("fHistAllRec_","counts, all particles");
  }
  TString histname = "fHistAllRectotETDep_" + fHistogramNameSuffix;
  fHistAllRectotETDep = new TH1F(histname.Data(),"total ET, all particles",fgNumOfEBins, fgEAxis);


  if(fMakeSparse){
    fHistElectronRecETDep = CreateChargedPartHistoSparse("fHistElectronRecETDep_","E_{T}, electrons");
    fHistElectronRec = CreateChargedPartHistoSparse("fHistElectronRec_","counts, electrons");
  }
  histname = "fHistElectronMatchtotETDep_" + fHistogramNameSuffix;
  fHistElectronMatchtotETDep = new TH1F(histname.Data(),"total ET, MC primary Electrons",fgNumOfEBins, fgEAxis);

  histname = "fHistElectronRecdEdxP_" + fHistogramNameSuffix;
  fHistElectronRecdEdxP = new TH2F(histname,"TPC dEdx vs P",100,0.,10.,100,0.,200.);

  histname = "fHistNeutralRectotET_" + fHistogramNameSuffix;
  fHistNeutralRectotET = new TH1F(histname.Data(),"total ET, neutral particles",fgNumOfEBins, fgEAxis);
	
  histname = "fHistTotEMRectotET_" + fHistogramNameSuffix;
  fHistTotEMRectotET = new TH1F(histname.Data(),"total electromagnetic ET",fgNumOfEBins, fgEAxis);

  if(fMakeSparse){
    fHistMuonRecETDep = CreateChargedPartHistoSparse("fHistMuonRecETDep_","E_{T}, muons");
    fHistMuonRec = CreateChargedPartHistoSparse("fHistMuonRec_","counts, muons");
  }
  histname = "fHistMuonMatchtotETDep_" + fHistogramNameSuffix;
  fHistMuonMatchtotETDep = new TH1F(histname.Data(),"total ET, Muons",fgNumOfEBins, fgEAxis);

  histname = "fHistMuonRecdEdxP_" + fHistogramNameSuffix;
  fHistMuonRecdEdxP = new TH2F(histname,"TPC dEdx vs P",100,0.,10.,100,0.,200.);
	
  if(fMakeSparse){
    fHistPionRecETDep = CreateChargedPartHistoSparse("fHistPionRecETDep_","E_{T}, pions");
    fHistPionRec = CreateChargedPartHistoSparse("fHistPionRec_","counts, pions");
  }
  histname = "fHistPionMatchtotETDep_" + fHistogramNameSuffix;
  fHistPionMatchtotETDep = new TH1F(histname.Data(),"total ET, Pions",fgNumOfEBins, fgEAxis);
  histname = "fHistPionRecdEdxP_" + fHistogramNameSuffix;
  fHistPionRecdEdxP = new TH2F(histname,"TPC dEdx vs P",100,0.,10.,100,0.,200.);
	
  if(fMakeSparse){
    fHistKaonRecETDep = CreateChargedPartHistoSparse("fHistKaonRecETDep_","E_{T}, kaons");
    fHistKaonRec = CreateChargedPartHistoSparse("fHistKaonRec_","counts, kaons");
  }
  histname = "fHistKaonMatchtotETDep_" + fHistogramNameSuffix;
  fHistKaonMatchtotETDep = new TH1F(histname.Data(),"total ET, Kaons",fgNumOfEBins, fgEAxis);

  histname = "fHistKaonRecdEdxP_" + fHistogramNameSuffix;
  fHistKaonRecdEdxP = new TH2F(histname,"TPC dEdx vs P",100,0.,10.,100,0.,200.);
	
  if(fMakeSparse){
    fHistProtonRecETDep = CreateChargedPartHistoSparse("fHistProtonRecETDep_","E_{T}, protons");
    fHistProtonRec = CreateChargedPartHistoSparse("fHistProtonRec_","counts, protons");
  }
  histname = "fHistProtonMatchtotETDep_" + fHistogramNameSuffix;
  fHistProtonMatchtotETDep = new TH1F(histname.Data(),"total ET, Protons",fgNumOfEBins, fgEAxis);

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
	
	
  if(fMakeSparse){
    list->Add(fHistAllRecETDep); 
    list->Add(fHistAllRec); 
  }
  list->Add(fHistAllRectotETDep); 

  if(fMakeSparse){
    list->Add(fHistElectronRecETDep); 
    list->Add(fHistElectronRec); 
  }
  list->Add(fHistElectronMatchtotETDep); 
  list->Add(fHistElectronRecdEdxP);
	
	
  list->Add(fHistTotEMRectotET); 

  list->Add(fHistMuonRec); 
  list->Add(fHistMuonRecdEdxP);
  list->Add(fHistMuonMatchtotETDep); 

  if(fMakeSparse){
    list->Add(fHistPionRecETDep); 
    list->Add(fHistPionRec); 
  }
  list->Add(fHistPionMatchtotETDep); 
  list->Add(fHistPionRecdEdxP);
	
  if(fMakeSparse){
    list->Add(fHistKaonRecETDep); 
    list->Add(fHistKaonRec); 
  }
  list->Add(fHistKaonMatchtotETDep); 
  list->Add(fHistKaonRecdEdxP);

  if(fMakeSparse){
    list->Add(fHistProtonRecETDep); 
    list->Add(fHistProtonRec); 
  }
  list->Add(fHistProtonMatchtotETDep); 
  list->Add(fHistProtonRecdEdxP);
	
  list->Add(fHistTotChargedMatchtotETDep); 
  list->Add(fHistTotalRectotETDep); 
	
  list->Add(fHistDeltaRZ);
}

//________________________________________________________________________
//project to a EMCal radius
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
//project to a cluster position
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
