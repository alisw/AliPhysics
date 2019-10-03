//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis
//  - reconstruction output
//  implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEtReconstructed.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include "AliEMCALTrack.h"
#include "AliESDCaloCluster.h"
#include "TVector3.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliVParticle.h"
#include "TDatabasePDG.h"
#include "TList.h"
#include "AliESDpid.h"
#include <iostream>
#include "TH3F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH1I.h"
#include "TFile.h"
#include "AliAnalysisHadEtCorrections.h"
#include "AliAnalysisEtSelector.h"
#include "AliLog.h"
#include "AliCentrality.h"
#include "AliPHOSGeoUtils.h"
#include "AliPHOSGeometry.h"
#include "AliAnalysisEtRecEffCorrection.h"
#include "AliESDpid.h"


using namespace std;

ClassImp(AliAnalysisEtReconstructed);


AliAnalysisEtReconstructed::AliAnalysisEtReconstructed() :
        AliAnalysisEt()
	,fQATree(0)
	,fMakeQATree(0)
	,fClusterMultiplicity(0)
	,fTrackMultiplicity(0)
	,fEventID(0)
        ,fCorrections(0)
        ,fPidCut(0)
	,nChargedHadronsMeasured(0)
	,nChargedHadronsTotal(0)
        ,fHistChargedPionEnergyDeposit(0)
        ,fHistProtonEnergyDeposit(0)
        ,fHistAntiProtonEnergyDeposit(0)
        ,fHistChargedKaonEnergyDeposit(0)
        ,fHistMuonEnergyDeposit(0)
        ,fHistRemovedEnergy(0)
        ,fGeomCorrection(1.0)
        ,fEMinCorrection(1.0/0.687)
	,fRecEffCorrection(1.0)
	,fClusterPositionAccepted(0)
	,fClusterPositionAll(0)
	,fClusterPositionAcceptedEnergy(0)
	,fClusterPositionAllEnergy(0)
	,fClusterEnergy(0)
	,fClusterEnergyCent(0)
	,fClusterEnergyModifiedTrackMatchesCent(0)
	,fClusterEnergyCentMatched(0)
	,fClusterEnergyCentNotMatched(0)
	,fClusterEt(0)
	,fHistChargedEnergyRemoved(0)
	,fHistNeutralEnergyRemoved(0)
	,fHistGammaEnergyAdded(0)
	,fHistMatchedTracksEvspTvsCent(0)
	,fHistMatchedTracksEvspTvsCentEffCorr(0)
	,fHistMatchedTracksEvspTvsCentEffTMCorr(0)
	,fHistPeripheralMatchedTracksEvspTvsCentEffTMCorr(0)
	,fHistMatchedTracksEvspTvsCentEffTMCorr500MeV(0)
	,fHistFoundHadronsvsCent(0)
	,fHistNotFoundHadronsvsCent(0)
	,fHistFoundHadronsEtvsCent(0)
	,fHistNotFoundHadronsEtvsCent(0)
	,fHistFoundHadronsvsCent500MeV(0)
	,fHistNotFoundHadronsvsCent500MeV(0)
	,fHistFoundHadronsEtvsCent500MeV(0)
	,fHistNotFoundHadronsEtvsCent500MeV(0)
	,fHistNominalRawEt(0)
	,fHistNominalNonLinHighEt(0)
	,fHistNominalNonLinLowEt(0)
	,fHistNominalEffHighEt(0)
	,fHistNominalEffLowEt(0)
	,fHistTotRawEtEffCorr(0)
	,fHistTotRawEt(0)
	,fHistTotRawEtEffCorr500MeV(0)
	,fHistTotAllRawEt(0)
	,fHistTotAllRawEtEffCorr(0)
	,fHistNUsedClusters(0)
	,fHistNClustersPhosVsEmcal(0)
	,fHistClusterSizeVsCent(0)
	,fHistMatchedClusterSizeVsCent(0)
	,fHistTotAllRawEtVsTotalPt(0)
	,fHistTotAllRawEtVsTotalPtVsCent(0)
	,fHistTotMatchedRawEtVsTotalPtVsCent(0)
	,fHistPIDProtonsTrackMatchedDepositedVsNch(0)
	,fHistPIDAntiProtonsTrackMatchedDepositedVsNch(0)
	,fHistPIDProtonsTrackMatchedDepositedVsNcl(0)
	,fHistPIDAntiProtonsTrackMatchedDepositedVsNcl(0)
	,fHistPiKPTrackMatchedDepositedVsNch(0)
						      //,
	,fHistPIDProtonsTrackMatchedDepositedVsNchNoEff(0)
	,fHistPIDAntiProtonsTrackMatchedDepositedVsNchNoEff(0)
	,fHistPIDProtonsTrackMatchedDepositedVsNclNoEff(0)
	,fHistPIDAntiProtonsTrackMatchedDepositedVsNclNoEff(0)
	,fHistPiKPTrackMatchedDepositedVsNchNoEff(0)
	,fHistCentVsNchVsNclReco(0)
	,fHistRawSignalReco(0)
	,fHistEffCorrSignalReco(0)
	,fHistRecoRCorrVsPtVsCent(0)
{

}

AliAnalysisEtReconstructed::~AliAnalysisEtReconstructed()
{//destructor
  if(fMakeQATree){
    //fQATree->Clear();
    delete fQATree;
  }
    delete fCorrections;
    delete fHistChargedPionEnergyDeposit; /** Energy deposited in calorimeter by charged pions */
    delete fHistProtonEnergyDeposit; /** Energy deposited in calorimeter by protons */
    delete fHistAntiProtonEnergyDeposit; /** Energy deposited in calorimeter by anti-protons */
    delete fHistChargedKaonEnergyDeposit; /** Energy deposited in calorimeter by charged kaons */
    delete fHistMuonEnergyDeposit; /** Energy deposited in calorimeter by muons */

    delete fHistRemovedEnergy; // removed energy
    delete fClusterPositionAccepted;
    delete fClusterPositionAll;
    delete fClusterPositionAcceptedEnergy;
    delete fClusterPositionAllEnergy;
    delete fClusterEnergy;
    delete fClusterEnergyCent;
    delete fClusterEnergyModifiedTrackMatchesCent;
    delete fClusterEnergyCentMatched;
    delete fClusterEnergyCentNotMatched;
    delete fClusterEt;
    delete fHistChargedEnergyRemoved;
    delete fHistNeutralEnergyRemoved;
    delete fHistGammaEnergyAdded;
    delete fHistMatchedTracksEvspTvsCent;
    delete fHistMatchedTracksEvspTvsCentEffCorr;
    delete fHistMatchedTracksEvspTvsCentEffTMCorr;
    delete fHistPeripheralMatchedTracksEvspTvsCentEffTMCorr;
    delete fHistMatchedTracksEvspTvsCentEffTMCorr500MeV;
    delete fHistFoundHadronsvsCent;
    delete fHistNotFoundHadronsvsCent;
    delete fHistFoundHadronsEtvsCent;
    delete fHistNotFoundHadronsEtvsCent;
    delete fHistFoundHadronsvsCent500MeV;
    delete fHistNotFoundHadronsvsCent500MeV;
    delete fHistFoundHadronsEtvsCent500MeV;
    delete fHistNotFoundHadronsEtvsCent500MeV;
    delete fHistNominalRawEt;
    delete fHistNominalNonLinHighEt;
    delete fHistNominalNonLinLowEt;
    delete fHistNominalEffHighEt;
    delete fHistNominalEffLowEt;
    delete fHistTotRawEtEffCorr;
    delete fHistTotRawEt;
    delete fHistTotAllRawEt;
    delete fHistTotAllRawEtEffCorr;
    delete fHistTotRawEtEffCorr500MeV;
    delete fHistNUsedClusters;
    delete fHistNClustersPhosVsEmcal;
    delete fHistClusterSizeVsCent;
    delete fHistMatchedClusterSizeVsCent;
    delete fHistTotAllRawEtVsTotalPt;
    delete fHistTotAllRawEtVsTotalPtVsCent;
    delete fHistTotMatchedRawEtVsTotalPtVsCent;
    delete fHistPIDProtonsTrackMatchedDepositedVsNch;
    delete fHistPIDAntiProtonsTrackMatchedDepositedVsNch;
    delete fHistPIDProtonsTrackMatchedDepositedVsNcl;
    delete fHistPIDAntiProtonsTrackMatchedDepositedVsNcl;
    delete fHistPiKPTrackMatchedDepositedVsNch;
    delete fHistPIDProtonsTrackMatchedDepositedVsNchNoEff;
    delete fHistPIDAntiProtonsTrackMatchedDepositedVsNchNoEff;
    delete fHistPIDProtonsTrackMatchedDepositedVsNclNoEff;
    delete fHistPIDAntiProtonsTrackMatchedDepositedVsNclNoEff;
    delete fHistPiKPTrackMatchedDepositedVsNchNoEff;
    delete fHistCentVsNchVsNclReco;
    delete fHistRawSignalReco;
    delete fHistEffCorrSignalReco;
    delete fHistRecoRCorrVsPtVsCent;
}

Int_t AliAnalysisEtReconstructed::AnalyseEvent(AliVEvent* ev)
{

    //AliAnalysisEt::AnalyseEvent(ev);
    // analyse ESD event
    ResetEventValues();
    if (!ev) {
        AliFatal("ERROR: Event does not exist");
        return 0;
    }

    AliESDEvent *event = dynamic_cast<AliESDEvent*>(ev);
    if (!event) {
        AliFatal("ERROR: ESD Event does not exist");
        return 0;
    }
    if(!fSelector){
        AliFatal("ERROR: fSelector does not exist");
        return 0;
    }
    fSelector->SetEvent(event);
    
    Int_t cent = -1;
    fCentrality = event->GetCentrality();
    if (fCentrality && cent)
    {
        cent = fCentrality->GetCentralityClass5("V0M");
        fCentClass = fCentrality->GetCentralityClass5("V0M");
    }


  //for PID
  //AliESDpid *pID = new AliESDpid();
  //pID->MakePID(event);
  Float_t etPIDProtons = 0.0;
  Float_t etPIDAntiProtons = 0.0;
  Float_t etPiKPMatched = 0.0;
  Float_t etPIDProtonsNoEff = 0.0;
  Float_t etPIDAntiProtonsNoEff = 0.0;
  Float_t etPiKPMatchedNoEff = 0.0;
  Float_t multiplicity = fEsdtrackCutsTPC->GetReferenceMultiplicity(event,kTRUE);
  fTrackMultiplicity = multiplicity;
  fEventID = event->GetPeriodNumber();//This is not the event id

    Float_t totalMatchedPt = 0.0;
    Float_t totalPt = 0.0;
    TObjArray* list  = fEsdtrackCutsTPC->GetAcceptedTracks(event);
    Int_t nGoodTracks = list->GetEntries();
    for (Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++){
	AliESDtrack *track = dynamic_cast<AliESDtrack*> (list->At(iTrack));
	if (!track)
	  {
	    Printf("ERROR: Could not get track %d", iTrack);
	    continue;
	  }
	else{
	  totalPt +=track->Pt();
	  //pID->MakeITSPID(track);


	}
    }

    //TRefArray *caloClusters = fSelector->GetClusters();//just gets the correct set of clusters - does not apply any cuts
    //Float_t fClusterMult = caloClusters->GetEntries();

    Float_t nominalRawEt = 0;
    Float_t totEt500MeV = 0;
    Float_t nonlinHighRawEt = 0;
    Float_t nonlinLowRawEt = 0;
    Float_t effHighRawEt = 0;
    Float_t effLowRawEt = 0;
    Float_t uncorrEt = 0;
    Float_t rawSignal = 0;
    Float_t effCorrSignal = 0;

    nChargedHadronsMeasured = 0.0;
    nChargedHadronsTotal = 0.0;
    Float_t nChargedHadronsEtMeasured = 0.0;
    Float_t nChargedHadronsEtTotal = 0.0;
    Float_t nChargedHadronsMeasured500MeV = 0.0;
    Float_t nChargedHadronsTotal500MeV = 0.0;
    Float_t nChargedHadronsEtMeasured500MeV = 0.0;
    Float_t nChargedHadronsEtTotal500MeV = 0.0;
    Float_t fTotAllRawEt = 0.0;
    Float_t fTotRawEt = 0.0;
    Float_t fTotRawEtEffCorr = 0.0;
    Float_t fTotAllRawEtEffCorr = 0.0;
    Int_t nPhosClusters = 0;
    Int_t nEmcalClusters = 0;

    Int_t nUsedClusters = 0;

    TRefArray *caloClusters = fSelector->GetClusters();
    Int_t nCluster = caloClusters->GetEntries();
    fClusterMultiplicity = nCluster;
    //if we are making the QA tree and the cluster multiplicity (for PHOS) is less than expected, fill the QA tree so we know what event it was
    if(fMakeQATree && fClusterMultiplicity < 5e-3*fTrackMultiplicity-1.5) fQATree->Fill();

    for (int iCluster = 0; iCluster < nCluster; iCluster++ )
    {
        AliESDCaloCluster* cluster = ( AliESDCaloCluster* )caloClusters->At( iCluster );
        if (!cluster)
        {
            AliError(Form("ERROR: Could not get cluster %d", iCluster));
            continue;
        }
// 	if(!fSelector->CutGeometricalAcceptance(*cluster)){
// 	  Float_t pos[3];
// 	  cluster->GetPosition(pos);
// 	  TVector3 cp(pos);
// 	  cout<<"rejected cluster phi "<<cp.Phi()*180.0/TMath::Pi()<<" eta "<<cp.Eta()<<endl;
// 	}
        int x = 0;
	fCutFlow->Fill(x++);//fills 0
	if(cluster->IsEMCAL()) nEmcalClusters++;
	else nPhosClusters++;
	if(!fSelector->IsDetectorCluster(*cluster)) continue;
	fCutFlow->Fill(x++);//fills 1
	if(!fSelector->PassMinEnergyCut(*cluster)) continue;
	fCutFlow->Fill(x++);//fills 2
        if (!fSelector->PassDistanceToBadChannelCut(*cluster)) continue;
	fCutFlow->Fill(x++);//fills 3
        if (!fSelector->CutGeometricalAcceptance(*cluster)) continue;
	//fCutFlow->Fill(x++);
        Float_t pos[3];

        cluster->GetPosition(pos);
        TVector3 cp(pos);
	fClusterPositionAll->Fill(cp.Phi(), cp.PseudoRapidity());
	Float_t fReconstructedE = cluster->E();
	Float_t lostEnergy = 0.0;
	Float_t lostTrackPt = 0.0;
	fClusterPositionAllEnergy->Fill(cp.Phi(), cp.PseudoRapidity(),GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE);

	//if(TMath::Abs(cp.Eta())> fCuts->fCuts->GetGeometryEmcalEtaAccCut() || cp.Phi() >  fCuts->GetGeometryEmcalPhiAccMaxCut()*TMath::Pi()/180. ||  cp.Phi() >  fCuts->GetGeometryEmcalPhiAccMinCut()*TMath::Pi()/180.) continue;//Do not accept if cluster is not in the acceptance
	fTotAllRawEt += TMath::Sin(cp.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
	fTotAllRawEtEffCorr +=GetCorrectionModification(*cluster,0,0,cent)* CorrectForReconstructionEfficiency(*cluster,cent);

	fClusterEnergyCent->Fill(GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE,cent);
        Bool_t matched = kTRUE;//default to no track matched
	Bool_t countasmatched = kFALSE;
	Bool_t correctedcluster = kFALSE;

	Int_t trackMatchedIndex = cluster->GetTrackMatchedIndex();//find the index of the matched track
	matched = !(fSelector->PassTrackMatchingCut(*cluster));//PassTrackMatchingCut is false if there is a matched track
	if(matched){//if the track match is good (, is the track good?
	  if(trackMatchedIndex < 0) matched=kFALSE;//If the index is bad, don't count it
	  if(matched){
	    AliESDtrack *track = event->GetTrack(trackMatchedIndex);
	    //if this is a good track, accept track will return true.  The track matched is good, so not track matched is false
	    matched = fEsdtrackCutsTPC->AcceptTrack(track);//If the track is bad, don't count it
	    if(matched){//if it is still matched see if the track p was less than the energy
	      Float_t rcorr = TMath::Sin(cp.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
	      fHistRecoRCorrVsPtVsCent->Fill(rcorr,track->Pt(), fCentClass);
	      if(fSelector->PassMinEnergyCut( (fReconstructedE - fsub* track->P())*TMath::Sin(cp.Theta()) )  ){//if more energy was deposited than the momentum of the track  and more than one particle led to the cluster
		// 	      if(fReconstructedE - fsub* track->P() > 0.0){
		//cout<<"match corrected"<<endl;
		fReconstructedE = fReconstructedE - fsub* track->P();
		matched = kFALSE;
		countasmatched = kTRUE;
		correctedcluster = kTRUE;
		lostEnergy = fsub* track->P();
		lostTrackPt = track->Pt();
		fClusterEnergyModifiedTrackMatchesCent->Fill(GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE,cent);
	      }
// 	      else{
// 		    cerr<<"match passed ";
// 		    cerr<<"E "<<fReconstructedE<<" fsubmeanhad "<<fsub<<" p "<< track->P();
// 		    if(correctedcluster) cout<<" corrected";
// 		    cerr<<endl;
// 	      }
	    }
	  }
	}


        if (matched)
        {
	  
            if (cluster->GetNTracksMatched() > 0 && trackMatchedIndex>=0)
            {
                AliVTrack *track = event->GetTrack(trackMatchedIndex);
                if (!track) {
                    AliError("Error: track does not exist");
                }
                else {
		  totalMatchedPt +=track->Pt();
		  fClusterEnergyCentMatched->Fill(GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE,cent);
		  fHistMatchedClusterSizeVsCent->Fill(cluster->GetNCells(),cent);

		  float eff = fTmCorrections->TrackMatchingEfficiency(track->Pt(),cent);
		  if(TMath::Abs(eff)<1e-5) eff = 1.0;
		  //cout<<"pt "<<track->Pt()<<" eff "<<eff<<" total "<<nChargedHadronsTotal<<endl;
		  nChargedHadronsMeasured++;
		  nChargedHadronsTotal += 1/eff;
		  Double_t effCorrEt = GetCorrectionModification(*cluster,0,0,cent) * CorrectForReconstructionEfficiency(*cluster,fReconstructedE,cent);
		  nChargedHadronsEtMeasured+= TMath::Sin(cp.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
		  //One efficiency is the gamma efficiency and the other is the track matching efficiency.
		  nChargedHadronsEtTotal+= 1/eff *TMath::Sin(cp.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
		  //cout<<"nFound "<<1<<" nFoundTotal "<<1/eff<<" etMeas "<<TMath::Sin(cp.Theta())*fReconstructedE<<" ET total "<< 1/eff *TMath::Sin(cp.Theta())*fReconstructedE<<endl;

		  Float_t nSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion); 
		  Float_t nSigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton); 
		  bool isProton = (nSigmaPion>3.0 && nSigmaProton<3.0 && track->Pt()<0.9);
		  //cout<<"NSigmaProton "<<nSigmaProton<<endl;
		  etPiKPMatched += effCorrEt;
		  etPiKPMatchedNoEff  +=TMath::Sin(cp.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
		  if(isProton){
		    if(track->Charge()>0){
		      etPIDProtons += effCorrEt;
		      etPIDProtonsNoEff +=TMath::Sin(cp.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
		    }
		    else{
		      etPIDAntiProtonsNoEff +=TMath::Sin(cp.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
		      etPIDAntiProtons += effCorrEt;
		    }
		  }
		  if(TMath::Sin(cp.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE>0.5){
		    nChargedHadronsMeasured500MeV++;
		    nChargedHadronsTotal500MeV += 1/eff;
		    nChargedHadronsEtMeasured500MeV+= TMath::Sin(cp.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
		    nChargedHadronsEtTotal500MeV+= 1/eff *TMath::Sin(cp.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
		  }
		  cluster->GetPosition(pos);	  
		  TVector3 p2(pos);
		  uncorrEt += TMath::Sin(p2.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
		  if(correctedcluster || fReconstructedE <fsubmeanhade* track->P() ){//if more energy was deposited than the momentum of the track  and more than one particle led to the cluster and the corrected energy is greater than zero
		    fHistMatchedTracksEvspTvsCent->Fill(track->P(),TMath::Sin(cp.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE,cent);
		    fHistMatchedTracksEvspTvsCentEffCorr->Fill(track->P(),effCorrEt,cent);
		    //Weighed by the number of tracks we didn't find
		    fHistMatchedTracksEvspTvsCentEffTMCorr->Fill(track->P(), effCorrEt,cent, (1/eff-1) );
		    if(cent<16 && cent>11){//centralities 60-80% where false track matches are low
		      for(int cbtest = 0; cbtest<20; cbtest++){//then we calculate the deposit matched to hadrons with different centrality bins' efficiencies
			float efftest = fTmCorrections->TrackMatchingEfficiency(track->Pt(),cbtest);
			if(TMath::Abs(efftest)<1e-5) efftest = 1.0;
			Double_t effCorrEttest = GetCorrectionModification(*cluster,0,0,cent)*CorrectForReconstructionEfficiency(*cluster,fReconstructedE,cbtest);
			fHistPeripheralMatchedTracksEvspTvsCentEffTMCorr->Fill(track->P(), effCorrEttest,cbtest, (1/efftest-1) );
		      }
		    }
		    if(uncorrEt>=0.5) fHistMatchedTracksEvspTvsCentEffTMCorr500MeV->Fill(track->P(), effCorrEt,cent, (1/eff-1) );
		  }
		  const Double_t *pidWeights = track->PID();
		  
		  Double_t maxpidweight = 0;
		  Int_t maxpid = 0;
		  
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
		      if (fCuts->GetHistMakeTreeDeposit() && fDepositTree)
                        {
			  fEnergyDeposited =GetCorrectionModification(*cluster,0,0,cent)* fReconstructedE;
			  fMomentumTPC = track->P();
			  fCharge = track->Charge();
			  fParticlePid = maxpid;
			  fPidProb = maxpidweight;
			  AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
			  if (!esdTrack) {
			    AliError("Error: track does not exist");
			  }
			  else {
			    if (esdTrack) fTrackPassedCut = fEsdtrackCutsTPC->AcceptTrack(esdTrack);
			    fDepositTree->Fill();
			  }
                        }
		      
		      if (maxpidweight > fPidCut)
                        {
			  //Float_t dist = TMath::Sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
			  
			  //Float_t theta = TMath::ATan(pos[2]/dist)+TMath::Pi()/2;
			  
                            //Float_t et = fReconstructedE * TMath::Sin(theta);
                            if (maxpid == AliPID::kProton)
                            {

                                if (track->Charge() == 1)
                                {
                                    fHistProtonEnergyDeposit->Fill(GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE, track->E());
                                }
                                else if (track->Charge() == -1)
                                {
                                    fHistAntiProtonEnergyDeposit->Fill(GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE, track->E());
                                }
                            }
                            else if (maxpid == AliPID::kPion)
                            {
                                fHistChargedPionEnergyDeposit->Fill(GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE, track->E());
                            }
                            else if (maxpid == AliPID::kKaon)
                            {
                                fHistChargedKaonEnergyDeposit->Fill(GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE, track->E());
                            }
                            else if (maxpid == AliPID::kMuon)
                            {
                                fHistMuonEnergyDeposit->Fill(GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE, track->E());
                            }
                        }
                    }
                }
            }
            //continue;
        } // distance
        else{//these are clusters which were not track matched
	  fCutFlow->Fill(x++);
	  //std::cout << x++ << std::endl;
	  
	  //if (fReconstructedE >  fSingleCellEnergyCut && cluster->GetNCells() == fCuts->GetCommonSingleCell()) continue;
	  //if (fReconstructedE < fClusterEnergyCut) continue;
	  cluster->GetPosition(pos);
	  
	  TVector3 p2(pos);
	  if(countasmatched){//These are tracks where we partially subtracted the energy but we subtracted some energy
	    float eff = fTmCorrections->TrackMatchingEfficiency(lostTrackPt,cent);
	    if(TMath::Abs(eff)<1e-5) eff = 1.0;
	    //cout<<"pt "<<track->Pt()<<" eff "<<eff<<" total "<<nChargedHadronsTotal<<endl;
	    nChargedHadronsMeasured++;
	    nChargedHadronsTotal += 1/eff;
	    //Double_t effCorrEt = CorrectForReconstructionEfficiency(*cluster,lostEnergy,cent);
	    nChargedHadronsEtMeasured+= TMath::Sin(cp.Theta())*lostEnergy;
	    //One efficiency is the gamma efficiency and the other is the track matching efficiency.
	    nChargedHadronsEtTotal+= 1/eff *TMath::Sin(cp.Theta())*lostEnergy;	      
	  }
	  fClusterPositionAccepted->Fill(p2.Phi(), p2.PseudoRapidity());
	  fClusterPositionAcceptedEnergy->Fill(p2.Phi(), p2.PseudoRapidity(),GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE);
	  fClusterEnergy->Fill(GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE);
	  fClusterEnergyCentNotMatched->Fill(GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE,cent);
	  fHistClusterSizeVsCent->Fill(cluster->GetNCells(),cent);
	  fClusterEt->Fill(TMath::Sin(p2.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE,cent);
	  uncorrEt += TMath::Sin(p2.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
	  float myuncorrEt = TMath::Sin(p2.Theta())*GetCorrectionModification(*cluster,0,0,cent)*fReconstructedE;
	  fTotRawEt += myuncorrEt;
	  nUsedClusters++;
	  
	  Double_t effCorrEt = CorrectForReconstructionEfficiency(*cluster,fReconstructedE,cent)*GetCorrectionModification(*cluster,0,0,cent);
	  rawSignal += myuncorrEt;
	  effCorrSignal +=effCorrEt;
	  //cout<<"cluster energy "<<fReconstructedE<<" eff corr Et "<<effCorrEt<<endl;
	  fTotRawEtEffCorr += effCorrEt;
	  fTotNeutralEt += effCorrEt;
	  nominalRawEt += effCorrEt;
	  if(myuncorrEt>=0.5){
	    totEt500MeV += effCorrEt;
	    //cout<<"test "<<myuncorrEt<<"> 0.5"<<endl;
	  }
	  else{
	    //cout<<"test "<<myuncorrEt<<"< 0.5"<<endl;
	  }
	  nonlinHighRawEt += effCorrEt*GetCorrectionModification(*cluster,1,0,cent);
	  nonlinLowRawEt += effCorrEt*GetCorrectionModification(*cluster,-1,0,cent);
	  effHighRawEt += effCorrEt*GetCorrectionModification(*cluster,0,1,cent);
	  effLowRawEt += effCorrEt*GetCorrectionModification(*cluster,0,-1,cent);
	  fNeutralMultiplicity++;
        }
        fMultiplicity++;
    }
    
    fHistNUsedClusters->Fill(nUsedClusters,cent);
    // cout<<" cent "<<cent<<" total "<<nCluster<<" used "<<nUsedClusters<<endl;

    fHistRawSignalReco->Fill(rawSignal);
    fHistEffCorrSignalReco->Fill(effCorrSignal);

    fHistNClustersPhosVsEmcal->Fill(nPhosClusters,nEmcalClusters,cent);
    fChargedEnergyRemoved = GetChargedContribution(fNeutralMultiplicity);
    fNeutralEnergyRemoved = GetNeutralContribution(fNeutralMultiplicity);
    fHistChargedEnergyRemoved->Fill(fChargedEnergyRemoved, fNeutralMultiplicity);
    fHistNeutralEnergyRemoved->Fill(fNeutralEnergyRemoved, fNeutralMultiplicity);
    
    fGammaEnergyAdded = GetGammaContribution(fNeutralMultiplicity);
    fHistGammaEnergyAdded->Fill(fGammaEnergyAdded, fNeutralMultiplicity);

    //Double_t removedEnergy = GetChargedContribution(cent) + GetNeutralContribution(cent) + GetGammaContribution(cent) + GetSecondaryContribution(cent);//fNeutralMultiplicity
    Double_t removedEnergy = GetHadronContribution(cent) + GetNeutronContribution(cent) + GetKaonContribution(cent) + GetSecondaryContribution(cent);//fNeutralMultiplicity
    //cout<<" centbin "<<cent<<" removed energy ch "<< GetHadronContribution(cent)<<" n " << GetNeutronContribution(cent) <<" kaon "<< GetKaonContribution(cent)<<" secondary "<< GetSecondaryContribution(cent);//fNeutralMultiplicity
    //cout<<" test min et "<<fTmCorrections->GetMinEtCorrection(cent);
    //cout<<" test neutral "<<fTmCorrections->NeutralContr(cent);
    //cout<<" test neutron "<<fTmCorrections->NeutralContr(cent);
    fHistRemovedEnergy->Fill(removedEnergy);
    
    fTotNeutralEtAcc = fTotNeutralEt;
    //fHistTotRawEtEffCorr->Fill(fTotNeutralEt,cent);
    fHistTotRawEtEffCorr->Fill(fTotRawEtEffCorr,cent);
    fHistTotRawEt->Fill(fTotRawEt,cent);
    fHistTotAllRawEt->Fill(fTotAllRawEt,cent);
    fHistTotAllRawEtVsTotalPt->Fill(fTotAllRawEt,totalPt);
    fHistTotAllRawEtVsTotalPtVsCent->Fill(fTotAllRawEt,totalPt,cent);
    fHistTotMatchedRawEtVsTotalPtVsCent->Fill(fTotAllRawEt,totalMatchedPt,cent);
    fHistTotAllRawEtEffCorr->Fill(fTotAllRawEtEffCorr,cent);
    //cout<<"fTotAllRawEtEffCorr "<<fTotAllRawEtEffCorr<<" fTotAllRawEt "<<fTotAllRawEt<<" fTotRawEtEffCorr "<<fTotRawEtEffCorr<<"("<<fTotNeutralEt<<")"<<" fTotRawEt "<<fTotRawEt<<endl;
    //cout<<"uncorr "<<uncorrEt<<" raw "<<nominalRawEt<<" tot raw "<<fTotNeutralEt;
    //cout<<" raw "<<fTotNeutralEt<<" removed "<<removedEnergy<<" etmin "<<GetMinEtCorrection(cent)<<" final ";
    if(GetMinEtCorrection(cent)>0) fTotNeutralEt =  (fTotNeutralEt - removedEnergy)/GetMinEtCorrection(cent);
    //cout<<fTotNeutralEt<<endl;
    //cout<<" tot corr "<<fTotNeutralEt<<endl;
    fTotEt = fTotChargedEt + fTotNeutralEt;
// Fill the histograms...0
    FillHistograms();
    //std::cout << "fTotNeutralEt: " << fTotNeutralEt << ", Contribution from non-removed charged: " << GetChargedContribution(fNeutralMultiplicity) << ", neutral: " << GetNeutralContribution(fNeutralMultiplicity) << ", gammas: " << GetGammaContribution(fNeutralMultiplicity) << ", multiplicity: " << fNeutralMultiplicity<< std::endl;
    //cout<<"cent "<<cent<<" cluster mult "<<fClusterMult<<" fTotNeutralEt "<<fTotNeutralEt<<" nominalRawEt "<<nominalRawEt<<endl;
    fHistNominalRawEt->Fill(nominalRawEt,cent);
    fHistTotRawEtEffCorr500MeV->Fill(totEt500MeV,cent);
    fHistNominalNonLinHighEt->Fill(nonlinHighRawEt,cent);
    fHistNominalNonLinLowEt->Fill(nonlinLowRawEt,cent);
    fHistNominalEffHighEt->Fill(effHighRawEt,cent);
    fHistNominalEffLowEt->Fill(effLowRawEt,cent);
    fHistFoundHadronsvsCent->Fill(nChargedHadronsMeasured,cent);
    fHistNotFoundHadronsvsCent->Fill(nChargedHadronsTotal-nChargedHadronsMeasured,cent);
    fHistFoundHadronsEtvsCent->Fill(nChargedHadronsEtMeasured,cent);
    fHistNotFoundHadronsEtvsCent->Fill(nChargedHadronsEtTotal-nChargedHadronsEtMeasured,cent);
    //cout<<"found "<<nChargedHadronsMeasured<<" total "<<nChargedHadronsTotal<<" not found "<<nChargedHadronsTotal-nChargedHadronsMeasured<<" found "<< nChargedHadronsMeasured500MeV<<" not found "<< nChargedHadronsTotal500MeV-nChargedHadronsMeasured500MeV <<" total "<<nChargedHadronsTotal500MeV<<endl;
    fHistFoundHadronsvsCent500MeV->Fill(nChargedHadronsMeasured500MeV,cent);
    fHistNotFoundHadronsvsCent500MeV->Fill(nChargedHadronsTotal500MeV-nChargedHadronsMeasured500MeV,cent);
    fHistFoundHadronsEtvsCent500MeV->Fill(nChargedHadronsEtMeasured500MeV,cent);
    fHistNotFoundHadronsEtvsCent500MeV->Fill(nChargedHadronsEtTotal500MeV-nChargedHadronsEtMeasured500MeV,cent);
//     cout<<"Number of hadrons measured:  "<<nChargedHadronsMeasured<<" Estimated total number of hadrons "<<nChargedHadronsTotal<<" ET in track matched hadrons "<<
//       nChargedHadronsEtMeasured;
//     if(nChargedHadronsMeasured>0)cout<<" ("<<nChargedHadronsEtMeasured/nChargedHadronsMeasured<<") ";
//     cout<<" ET in all hadrons ";
//     cout<<nChargedHadronsEtTotal;
//     if(nChargedHadronsTotal>0) cout<<" ("<<nChargedHadronsEtTotal/nChargedHadronsTotal<<") ";
//     cout<<endl;
    fHistPIDProtonsTrackMatchedDepositedVsNch->Fill(etPIDProtons,multiplicity);
    fHistPIDAntiProtonsTrackMatchedDepositedVsNch->Fill(etPIDAntiProtons,multiplicity);
    fHistPIDProtonsTrackMatchedDepositedVsNcl->Fill(etPIDProtons,nCluster);
    fHistPIDAntiProtonsTrackMatchedDepositedVsNcl->Fill(etPIDAntiProtons,nCluster);
    fHistPIDProtonsTrackMatchedDepositedVsNchNoEff->Fill(etPIDProtonsNoEff,multiplicity);
    fHistPIDAntiProtonsTrackMatchedDepositedVsNchNoEff->Fill(etPIDAntiProtonsNoEff,multiplicity);
    fHistPIDProtonsTrackMatchedDepositedVsNclNoEff->Fill(etPIDProtonsNoEff,nCluster);
    fHistPIDAntiProtonsTrackMatchedDepositedVsNclNoEff->Fill(etPIDAntiProtonsNoEff,nCluster);
    fHistCentVsNchVsNclReco->Fill(cent,multiplicity,nCluster);
    fHistPiKPTrackMatchedDepositedVsNch->Fill(etPiKPMatched,multiplicity);
    fHistPiKPTrackMatchedDepositedVsNchNoEff->Fill(etPiKPMatchedNoEff,multiplicity);
    //delete pID;
    return 0;
}

bool AliAnalysisEtReconstructed::CheckGoodVertex(AliVParticle* track)
{ // check vertex

    Float_t bxy = 999.;
    Float_t bz = 999.;
    if (!track) {
        AliError("ERROR: no track");
        return kFALSE;
    }
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
    if (!esdTrack) {
        AliError("ERROR: no track");
        return kFALSE;
    }
    esdTrack->GetImpactParametersTPC(bxy,bz);


    bool status = (TMath::Abs(track->Xv()) < fCuts->GetReconstructedVertexXCut()) &&
                  (TMath::Abs(track->Yv()) < fCuts->GetReconstructedVertexYCut()) &&
                  (TMath::Abs(track->Zv()) < fCuts->GetReconstructedVertexZCut()) &&
                  (TMath::Abs(bxy) < fCuts->GetReconstructedIPxyCut()) &&
                  (TMath::Abs(bz) < fCuts->GetReconstructedIPzCut());

    return status;
}

void AliAnalysisEtReconstructed::Init()
{ // Init
    AliAnalysisEt::Init();
    fPidCut = fCuts->GetReconstructedPidCut();
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));
    if (!fCorrections) {
        cout<<"Warning!  You have not set corrections.  Your code will crash.  You have to set the corrections."<<endl;
    }
}

bool AliAnalysisEtReconstructed::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{ // propagate track to detector radius
    if (!track) {
        cout<<"Warning: track empty"<<endl;
        return kFALSE;
    }
    AliESDtrack *esdTrack= dynamic_cast<AliESDtrack*>(track);
    if (!esdTrack) {
        AliError("ERROR: no ESD track");
        return kFALSE;
    }
    // Printf("Propagating track: eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());

    Bool_t prop = esdTrack->PropagateTo(fDetectorRadius, magField);

    // if (prop) Printf("Track propagated, eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());
    return prop && fSelector->CutGeometricalAcceptance(*esdTrack);
}

void AliAnalysisEtReconstructed::FillOutputList(TList* list)
{ // add some extra histograms to the ones from base class
    AliAnalysisEt::FillOutputList(list);

    if(fMakeQATree){
        list->Add(fQATree);
    }
    list->Add(fHistChargedPionEnergyDeposit);
    list->Add(fHistProtonEnergyDeposit);
    list->Add(fHistAntiProtonEnergyDeposit);
    list->Add(fHistChargedKaonEnergyDeposit);
    list->Add(fHistMuonEnergyDeposit);

    list->Add(fHistRemovedEnergy);
    list->Add(fClusterPositionAccepted);
    list->Add(fClusterPositionAll);
    list->Add(fClusterPositionAcceptedEnergy);
    list->Add(fClusterPositionAllEnergy);
    list->Add(fClusterEnergy);
    list->Add(fClusterEnergyCent);
    list->Add(fClusterEnergyModifiedTrackMatchesCent);
    list->Add(fClusterEnergyCentMatched);
    list->Add(fClusterEnergyCentNotMatched);
    list->Add(fClusterEt);
    
    list->Add(fHistChargedEnergyRemoved);
    list->Add(fHistNeutralEnergyRemoved);
    list->Add(fHistGammaEnergyAdded);
    list->Add(fHistMatchedTracksEvspTvsCent);
    list->Add(fHistMatchedTracksEvspTvsCentEffCorr);
    list->Add(fHistMatchedTracksEvspTvsCentEffTMCorr);
    list->Add(fHistPeripheralMatchedTracksEvspTvsCentEffTMCorr);
    list->Add(fHistMatchedTracksEvspTvsCentEffTMCorr500MeV);
    list->Add(fHistFoundHadronsvsCent);
    list->Add(fHistNotFoundHadronsvsCent);
    list->Add(fHistFoundHadronsEtvsCent);
    list->Add(fHistNotFoundHadronsEtvsCent);
    list->Add(fHistFoundHadronsvsCent500MeV);
    list->Add(fHistNotFoundHadronsvsCent500MeV);
    list->Add(fHistFoundHadronsEtvsCent500MeV);
    list->Add(fHistNotFoundHadronsEtvsCent500MeV);
    list->Add(fHistNominalRawEt);
    list->Add(fHistNominalNonLinHighEt);
    list->Add(fHistNominalNonLinLowEt);
    list->Add(fHistNominalEffHighEt);
    list->Add(fHistNominalEffLowEt);
    list->Add(fHistTotRawEtEffCorr);
    list->Add(fHistTotRawEtEffCorr500MeV);
    list->Add(fHistTotAllRawEtEffCorr);
    list->Add(fHistTotRawEt);
    list->Add(fHistTotAllRawEt);
    list->Add(fHistNUsedClusters);
    list->Add(fHistNClustersPhosVsEmcal);
    list->Add(fHistClusterSizeVsCent);
    list->Add(fHistMatchedClusterSizeVsCent);
    list->Add(fHistTotAllRawEtVsTotalPt);
    list->Add(fHistTotAllRawEtVsTotalPtVsCent);
    list->Add(fHistTotMatchedRawEtVsTotalPtVsCent);
    list->Add(fHistPIDProtonsTrackMatchedDepositedVsNch);
    list->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNch);
    list->Add(fHistPIDProtonsTrackMatchedDepositedVsNcl);
    list->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNcl);
    list->Add(fHistPiKPTrackMatchedDepositedVsNch);
    list->Add(fHistPIDProtonsTrackMatchedDepositedVsNchNoEff);
    list->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNchNoEff);
    list->Add(fHistPIDProtonsTrackMatchedDepositedVsNclNoEff);
    list->Add(fHistPIDAntiProtonsTrackMatchedDepositedVsNclNoEff);
    list->Add(fHistPiKPTrackMatchedDepositedVsNchNoEff);
    list->Add(fHistCentVsNchVsNclReco);
    list->Add(fHistRawSignalReco);
    list->Add(fHistEffCorrSignalReco);
    list->Add(fHistRecoRCorrVsPtVsCent);
}

void AliAnalysisEtReconstructed::CreateHistograms()
{ // add some extra histograms to the ones from base class
    AliAnalysisEt::CreateHistograms();
    if(fMakeQATree){
        fQATree = new TTree("fQATree", "fQATree");
        fQATree->Branch("fClusterMultiplicity", &fClusterMultiplicity, "fClusterMultiplicity/I");
        fQATree->Branch("fTrackMultiplicity", &fTrackMultiplicity, "fTrackMultiplicity/I");
        fQATree->Branch("fEventID",&fEventID,"fEventID/I");
        fQATree->Branch("fCentClass",&fCentClass,"fCentClass/I");
    }

    Float_t scale = 1;//scale up histograms if EMCal 2011 so we have the right range
    if(fDataSet==2011   && !fHistogramNameSuffix.Contains("P")){
      scale = 2.5;
    }
    Int_t nbinsEt = 1000*scale;
    Double_t minEt = 0;
    Double_t maxEt = 10*scale;

    // possibly change histogram limits
//     if (fCuts) {
//         nbinsEt = fCuts->GetHistNbinsParticleEt();
//         minEt = fCuts->GetHistMinParticleEt();
//         maxEt = fCuts->GetHistMaxParticleEt();
//     }

    TString histname;
    histname = "fHistChargedPionEnergyDeposit" + fHistogramNameSuffix;
    fHistChargedPionEnergyDeposit = new TH2F(histname.Data(), "Energy deposited by #pi^{+/-}", nbinsEt, minEt, maxEt, nbinsEt, minEt, maxEt);
    fHistChargedPionEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    fHistChargedPionEnergyDeposit->SetYTitle("Energy of track");

    histname = "fHistProtonEnergyDeposit" + fHistogramNameSuffix;
    fHistProtonEnergyDeposit = new TH2F(histname.Data(), "Energy deposited by protons", nbinsEt, minEt, maxEt, nbinsEt, minEt, maxEt);
    fHistProtonEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    fHistProtonEnergyDeposit->SetYTitle("Energy of track");

    histname = "fHistAntiProtonEnergyDeposit" + fHistogramNameSuffix;
    fHistAntiProtonEnergyDeposit = new TH2F(histname.Data(), "Energy deposited by anti-protons", nbinsEt, minEt, maxEt, nbinsEt, minEt, maxEt);
    fHistAntiProtonEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    fHistAntiProtonEnergyDeposit->SetYTitle("Energy of track");

    histname = "fHistChargedKaonEnergyDeposit" + fHistogramNameSuffix;
    fHistChargedKaonEnergyDeposit = new TH2F(histname.Data(), "Energy deposited by K^{+/-}", nbinsEt, minEt, maxEt, nbinsEt, minEt, maxEt);
    fHistChargedKaonEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    fHistChargedKaonEnergyDeposit->SetYTitle("Energy of track");

    histname = "fHistMuonEnergyDeposit" + fHistogramNameSuffix;
    fHistMuonEnergyDeposit = new TH2F(histname.Data(), "Energy deposited by #mu^{+/-}", nbinsEt, minEt, maxEt, nbinsEt, minEt, maxEt);
    fHistMuonEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    fHistMuonEnergyDeposit->SetYTitle("Energy of track");

    histname = "fHistRemovedEnergy" + fHistogramNameSuffix;
    fHistRemovedEnergy = new TH1F(histname.Data(), histname.Data(), 1000, 0, 20);
    //fHistMuonEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    //fHistMuonEnergyDeposit->SetYTitle("Energy of track");

    histname = "fClusterPositionAccepted" + fHistogramNameSuffix;
    fClusterPositionAccepted = new TH2D(histname.Data(), "Position of accepted neutral clusters",300, -TMath::Pi(),TMath::Pi(), 100, -0.7 , 0.7);
    fClusterPositionAccepted->SetXTitle("#phi");
    fClusterPositionAccepted->SetYTitle("#eta");

    histname = "fClusterPositionAll" + fHistogramNameSuffix;
    fClusterPositionAll = new TH2D(histname.Data(), "Position of accepted neutral clusters",300, -TMath::Pi(),TMath::Pi(), 100, -0.7 , 0.7);
    fClusterPositionAll->SetXTitle("#phi");
    fClusterPositionAll->SetYTitle("#eta");

    histname = "fClusterPositionAcceptedEnergy" + fHistogramNameSuffix;
    fClusterPositionAcceptedEnergy = new TH2D(histname.Data(), "Position of accepted neutral clusters",300, -TMath::Pi(),TMath::Pi(), 100, -0.7 , 0.7);
    fClusterPositionAcceptedEnergy->SetXTitle("#phi");
    fClusterPositionAcceptedEnergy->SetYTitle("#eta");

    histname = "fClusterPositionAllEnergy" + fHistogramNameSuffix;
    fClusterPositionAllEnergy = new TH2D(histname.Data(), "Position of accepted neutral clusters",300, -TMath::Pi(),TMath::Pi(), 100, -0.7 , 0.7);
    fClusterPositionAllEnergy->SetXTitle("#phi");
    fClusterPositionAllEnergy->SetYTitle("#eta");

    histname = "fClusterEnergy" + fHistogramNameSuffix;
    fClusterEnergy = new TH1F(histname.Data(), histname.Data(), 100, 0, 5);
    fClusterEnergy->SetYTitle("Number of clusters");
    fClusterEnergy->SetXTitle("Energy of cluster");

    histname = "fClusterEnergyCent" + fHistogramNameSuffix;
    fClusterEnergyCent = new TH2F(histname.Data(), histname.Data(), 100, 0, 5,20,-0.5,19.5);
    fClusterEnergyCent->SetXTitle("Energy of cluster");
    fClusterEnergyCent->SetYTitle("Centrality Bin");
    fClusterEnergyCent->SetZTitle("Number of clusters");

    histname = "fClusterEnergyModifiedTrackMatchesCent" + fHistogramNameSuffix;
    fClusterEnergyModifiedTrackMatchesCent = new TH2F(histname.Data(), histname.Data(), 100, 0, 5,20,-0.5,19.5);
    fClusterEnergyModifiedTrackMatchesCent->SetXTitle("Energy of cluster");
    fClusterEnergyModifiedTrackMatchesCent->SetYTitle("Centrality Bin");
    fClusterEnergyModifiedTrackMatchesCent->SetZTitle("Number of clusters");

    histname = "fClusterEnergyCentMatched" + fHistogramNameSuffix;
    fClusterEnergyCentMatched = new TH2F(histname.Data(), histname.Data(), 100, 0, 5,20,-0.5,19.5);
    fClusterEnergyCentMatched->SetXTitle("Energy of cluster");
    fClusterEnergyCentMatched->SetYTitle("Centrality Bin");
    fClusterEnergyCentMatched->SetZTitle("Number of Clusters");

    histname = "fClusterEnergyCentNotMatched" + fHistogramNameSuffix;
    fClusterEnergyCentNotMatched = new TH2F(histname.Data(), histname.Data(), 100, 0, 5,20,-0.5,19.5);
    fClusterEnergyCentNotMatched->SetXTitle("Energy of cluster");
    fClusterEnergyCentNotMatched->SetYTitle("Centrality Bin");
    fClusterEnergyCentNotMatched->SetZTitle("Number of clusters");

    histname = "fClusterEt" + fHistogramNameSuffix;
    fClusterEt = new TH2F(histname.Data(), histname.Data(), 100, 0, 5,20,-0.5,19.5);
    fClusterEt->SetXTitle("Number of clusters");
    fClusterEt->SetYTitle("E_{T} of cluster");

    histname = "fHistChargedEnergyRemoved" + fHistogramNameSuffix;
    fHistChargedEnergyRemoved = new TH2D(histname.Data(), histname.Data(), 1000, .0, 30, 100, -0.5 , 99.5);

    histname = "fHistNeutralEnergyRemoved" + fHistogramNameSuffix;
    fHistNeutralEnergyRemoved = new TH2D(histname.Data(), histname.Data(), 1000, .0, 30, 100, -0.5 , 99.5);

    histname = "fHistGammaEnergyAdded" + fHistogramNameSuffix;
    fHistGammaEnergyAdded = new TH2D(histname.Data(), histname.Data(), 1000, .0, 30, 100, -0.5 , 99.5);

    fHistMatchedTracksEvspTvsCent = new TH3F("fHistMatchedTracksEvspTvsCent", "fHistMatchedTracksEvspTvsCent",100, 0, 3,100,0,3,20,-0.5,19.5);
    fHistMatchedTracksEvspTvsCentEffCorr = new TH3F("fHistMatchedTracksEvspTvsCentEffCorr", "fHistMatchedTracksEvspTvsCentEffCorr",100, 0, 3,100,0,3,20,-0.5,19.5);
    fHistMatchedTracksEvspTvsCentEffTMCorr = new TH3F("fHistMatchedTracksEvspTvsCentEffTMCorr", "fHistMatchedTracksEvspTvsCentEffTMCorr",100, 0, 3,100,0,3,20,-0.5,19.5);
    fHistPeripheralMatchedTracksEvspTvsCentEffTMCorr = new TH3F("fHistPeripheralMatchedTracksEvspTvsCentEffTMCorr", "fHistPeripheralMatchedTracksEvspTvsCentEffTMCorr",100, 0, 3,100,0,3,20,-0.5,19.5);
    fHistMatchedTracksEvspTvsCentEffTMCorr500MeV = new TH3F("fHistMatchedTracksEvspTvsCentEffTMCorr500MeV", "fHistMatchedTracksEvspTvsCentEffTMCorr500MeV",100, 0, 3,100,0,3,20,-0.5,19.5);

    float max = 200*scale;
    if(fHistogramNameSuffix.Contains("P")){max = 100;}
    fHistFoundHadronsvsCent = new TH2F("fHistFoundHadronsvsCent","fHistFoundHadronsvsCent",100,0,max,20,-0.5,19.5);
       fHistNotFoundHadronsvsCent = new TH2F("fHistNotFoundHadronsvsCent","fHistNotFoundHadronsvsCent",100,0,max,20,-0.5,19.5);
    fHistFoundHadronsEtvsCent = new TH2F("fHistFoundHadronsEtvsCent","fHistFoundHadronsEtvsCent",100,0,max,20,-0.5,19.5);
       fHistNotFoundHadronsEtvsCent = new TH2F("fHistNotFoundHadronsEtvsCent","fHistNotFoundHadronsEtvsCent",100,0,max,20,-0.5,19.5);
    fHistFoundHadronsvsCent500MeV = new TH2F("fHistFoundHadronsvsCent500MeV","fHistFoundHadronsvsCent500MeV",100,0,max,20,-0.5,19.5);
    fHistNotFoundHadronsvsCent500MeV = new TH2F("fHistNotFoundHadronsvsCent500MeV","fHistNotFoundHadronsvsCent500MeV",100,0,max,20,-0.5,19.5);
    fHistFoundHadronsEtvsCent500MeV = new TH2F("fHistFoundHadronsEtvsCent500MeV","fHistFoundHadronsEtvsCent500MeV",100,0,max,20,-0.5,19.5);
    fHistNotFoundHadronsEtvsCent500MeV = new TH2F("fHistNotFoundHadronsEtvsCent500MeV","fHistNotFoundHadronsEtvsCent500MeV",100,0,max,20,-0.5,19.5);

    Int_t nbinsForDistributions = 384;//2^7*3 - allows for easy rebinning
    fHistTotRawEtEffCorr = new TH2F("fHistTotRawEtEffCorr","fHistTotRawEtEffCorr",nbinsForDistributions,0,250*scale,20,-0.5,19.5);
    fHistTotRawEt = new TH2F("fHistTotRawEt","fHistTotRawEt",nbinsForDistributions,0,250*scale,20,-0.5,19.5);
    fHistTotRawEtEffCorr500MeV = new TH2F("fHistTotRawEtEffCorr500MeV","fHistTotRawEtEffCorr500MeV",nbinsForDistributions,0,250*scale,20,-0.5,19.5);
    fHistTotAllRawEt = new TH2F("fHistTotAllRawEt","fHistTotAllRawEt",nbinsForDistributions,0,250*scale,20,-0.5,19.5);
    fHistTotAllRawEtEffCorr = new TH2F("fHistTotAllRawEtEffCorr","fHistTotAllRawEtEffCorr",nbinsForDistributions,0,250*scale,20,-0.5,19.5);
    fHistNUsedClusters = new TH2F("fHistNUsedClusters","fHistNUsedClusters",500,0,500,20,-0.5,19);
    fHistNClustersPhosVsEmcal = new TH3F("fHistNClustersPhosVsEmcal","fHistNClustersPhosVsEmcal",50,0,50,250*scale,0,250*scale,20,-0.5,19);
    fHistClusterSizeVsCent = new TH2F("fHistClusterSizeVsCent","fHistClusterSizeVsCent",10,0.5,10.5,20,-0.5,19.5);
    fHistMatchedClusterSizeVsCent = new TH2F("fHistMatchedClusterSizeVsCent","fHistMatchedClusterSizeVsCent",10,0.5,10.5,20,-0.5,19.5);
    fHistTotAllRawEtVsTotalPt = new TH2F("fHistTotAllRawEtVsTotalPt","fHistTotAllRawEtVsTotalPt",125,0,250*scale,200,0,2000);
    fHistTotAllRawEtVsTotalPtVsCent = new TH3F("fHistTotAllRawEtVsTotalPtVsCent","fHistTotAllRawEtVsTotalPtVsCent",125,0,250*scale,200,0,2000,20,-0.5,19.5);
    fHistTotMatchedRawEtVsTotalPtVsCent = new TH3F("fHistTotMatchedRawEtVsTotalPtVsCent","fHistTotMatchedRawEtVsTotalPtVsCent",nbinsForDistributions,0,250*scale,100,0,200,20,-0.5,19.5);
    
    maxEt = 500*scale;
    histname = "fHistNominalRawEt" + fHistogramNameSuffix;
    fHistNominalRawEt = new TH2D(histname.Data(), histname.Data(),nbinsEt,minEt,maxEt,20,-0.5,19.5);
    histname = "fHistNominalNonLinHighEt" + fHistogramNameSuffix;
    fHistNominalNonLinHighEt = new TH2D(histname.Data(), histname.Data(),nbinsEt,minEt,maxEt,20,-0.5,19.5);
    histname = "fHistNominalNonLinLowEt" + fHistogramNameSuffix;
    fHistNominalNonLinLowEt = new TH2D(histname.Data(), histname.Data(),nbinsEt,minEt,maxEt,20,-0.5,19.5);
    histname = "fHistNominalEffHighEt" + fHistogramNameSuffix;
    fHistNominalEffHighEt = new TH2D(histname.Data(), histname.Data(),nbinsEt,minEt,maxEt,20,-0.5,19.5);
    histname = "fHistNominalEffLowEt" + fHistogramNameSuffix;
    fHistNominalEffLowEt = new TH2D(histname.Data(), histname.Data(),nbinsEt,minEt,maxEt,20,-0.5,19.5);

      Float_t maxEtRange = 25*scale;
      Float_t maxEtRangeHigh = 125*scale;
      Float_t minEtRange = 0;
      Int_t nbinsMult = 100;
      Float_t maxMult = 3000;
      Float_t minMult = 0;
      Int_t nbinsCl = 250*scale;
      Float_t maxCl = 500*scale;
      Float_t minCl = 0;
    fHistPIDProtonsTrackMatchedDepositedVsNch = new TH2F("fHistPIDProtonsTrackMatchedDepositedVsNch","PID'd protons deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsMult,minMult,maxMult);
    fHistPIDAntiProtonsTrackMatchedDepositedVsNch = new TH2F("fHistPIDAntiProtonsTrackMatchedDepositedVsNch","PID'd #bar{p} E_{T} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsMult,minMult,maxMult);
    fHistPIDProtonsTrackMatchedDepositedVsNcl = new TH2F("fHistPIDProtonsTrackMatchedDepositedVsNcl","PID'd protons deposited in calorimeter vs cluster multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsCl,minCl,maxCl);
    fHistPIDAntiProtonsTrackMatchedDepositedVsNcl = new TH2F("fHistPIDAntiProtonsTrackMatchedDepositedVsNcl","PID'd #bar{p} E_{T} deposited in calorimeter vs cluster multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsCl,minCl,maxCl);
    fHistPiKPTrackMatchedDepositedVsNch = new TH2F("fHistPiKPTrackMatchedDepositedVsNch","PiKP track matched",nbinsEt,minEtRange,maxEtRangeHigh,nbinsMult,minMult,maxMult);

    fHistPIDProtonsTrackMatchedDepositedVsNchNoEff = new TH2F("fHistPIDProtonsTrackMatchedDepositedVsNchNoEff","PID'd protons deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsMult,minMult,maxMult);
    fHistPIDAntiProtonsTrackMatchedDepositedVsNchNoEff = new TH2F("fHistPIDAntiProtonsTrackMatchedDepositedVsNchNoEff","PID'd #bar{p} E_{T} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsMult,minMult,maxMult);
    fHistPIDProtonsTrackMatchedDepositedVsNclNoEff = new TH2F("fHistPIDProtonsTrackMatchedDepositedVsNclNoEff","PID'd protons deposited in calorimeter vs cluster multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsCl,minCl,maxCl);
    fHistPIDAntiProtonsTrackMatchedDepositedVsNclNoEff = new TH2F("fHistPIDAntiProtonsTrackMatchedDepositedVsNclNoEff","PID'd #bar{p} E_{T} deposited in calorimeter vs cluster multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsCl,minCl,maxCl);
    fHistPiKPTrackMatchedDepositedVsNchNoEff = new TH2F("fHistPiKPTrackMatchedDepositedVsNchNoEff","PiKP track matched",nbinsEt,minEtRange,maxEtRangeHigh,nbinsMult,minMult,maxMult);


    fHistCentVsNchVsNclReco = new TH3F("fHistCentVsNchVsNclReco","Cent bin vs Nch Vs NCl",20,-0.5,19.5,nbinsMult,minMult,maxMult,nbinsCl,minCl,maxCl);

   fHistRawSignalReco = new TH1F("fHistRawSignalReco","fHistRawSignalReco",20,-0.5,19.5);
   fHistEffCorrSignalReco = new TH1F("fHistEffCorrSignalReco","fHistEffCorrSignalReco",20,-0.5,19.5);
   fHistRecoRCorrVsPtVsCent = new TH3F("fHistRecoRCorrVsPtVsCent","fHistRecoRCorrVsPtVsCent",72,0,2,50,0,10,20,-0.5,19.5);

}
Double_t AliAnalysisEtReconstructed::ApplyModifiedCorrections(const AliESDCaloCluster& cluster,Int_t nonLinCorr, Int_t effCorr, Int_t cent)
{
  Float_t pos[3];
  cluster.GetPosition(pos);
  TVector3 cp(pos);
  Double_t corrEnergy = fReCorrections->CorrectedEnergy(cluster.E(),cent);
  
  Double_t factorNonLin = GetCorrectionModification(cluster, nonLinCorr,effCorr,cent);

    cout<<"Warning:  This function should not get called!"<<endl;
  //std::cout << "Original energy: " << cluster.E() << ", corrected energy: " << corrEnergy << std::endl;
  return TMath::Sin(cp.Theta())*corrEnergy*factorNonLin;
}

Double_t AliAnalysisEtReconstructed::GetCorrectionModification(const AliESDCaloCluster& cluster,Int_t nonLinCorr, Int_t effCorr, Int_t cent){//nonLinCorr 0 = nominal 1 = high -1 = low, effCorr  0 = nominal 1 = high -1 = low
  if(nonLinCorr==0){
    cout<<"Warning:  This function should not get called!"<<endl;//this statement is basically here to avoid a compilation warning
  }
  if(effCorr==0){
    cout<<"Warning:  This function should not get called!"<<endl;//this statement is basically here to avoid a compilation warning
  }
  return cluster.E()*cent;
}
