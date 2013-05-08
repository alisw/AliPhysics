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


using namespace std;

ClassImp(AliAnalysisEtReconstructed);


AliAnalysisEtReconstructed::AliAnalysisEtReconstructed() :
        AliAnalysisEt()
        ,fCorrections(0)
        ,fPidCut(0)
        ,fHistChargedPionEnergyDeposit(0)
        ,fHistProtonEnergyDeposit(0)
        ,fHistAntiProtonEnergyDeposit(0)
        ,fHistChargedKaonEnergyDeposit(0)
        ,fHistMuonEnergyDeposit(0)
        ,fHistRemovedEnergy(0)
        ,fGeomCorrection(1.0)
        ,fEMinCorrection(1.0/0.687)
	,fRecEffCorrection(1.0)
	,fClusterPosition(0)
	,fClusterEnergy(0)
	,fClusterEt(0)
	,fHistChargedEnergyRemoved(0)
	,fHistNeutralEnergyRemoved(0)
	,fHistGammaEnergyAdded(0)
	,fHistMatchedTracksEvspTvsMult(0)
	,fHistMatchedTracksEvspTvsMultEffCorr(0)
	,fHistFoundHadronsvsCent(0)
	,fHistNotFoundHadronsvsCent(0)
	,fHistFoundHadronsEtvsCent(0)
	,fHistNotFoundHadronsEtvsCent(0)
	,fHistNominalRawEt(0)
	,fHistNominalNonLinHighEt(0)
	,fHistNominalNonLinLowEt(0)
	,fHistNominalEffHighEt(0)
	,fHistNominalEffLowEt(0)
{

}

AliAnalysisEtReconstructed::~AliAnalysisEtReconstructed()
{//destructor
    delete fCorrections;
    delete fHistChargedPionEnergyDeposit; /** Energy deposited in calorimeter by charged pions */
    delete fHistProtonEnergyDeposit; /** Energy deposited in calorimeter by protons */
    delete fHistAntiProtonEnergyDeposit; /** Energy deposited in calorimeter by anti-protons */
    delete fHistChargedKaonEnergyDeposit; /** Energy deposited in calorimeter by charged kaons */
    delete fHistMuonEnergyDeposit; /** Energy deposited in calorimeter by muons */

    delete fHistRemovedEnergy; // removed energy
    delete fClusterPosition;
    delete fClusterEnergy;
    delete fClusterEt;
    delete fHistChargedEnergyRemoved;
    delete fHistNeutralEnergyRemoved;
    delete fHistGammaEnergyAdded;
    delete fHistMatchedTracksEvspTvsMult;
    delete fHistMatchedTracksEvspTvsMultEffCorr;
    delete fHistFoundHadronsvsCent;
    delete fHistNotFoundHadronsvsCent;
    delete fHistFoundHadronsEtvsCent;
    delete fHistNotFoundHadronsEtvsCent;
    delete fHistNominalRawEt;
    delete fHistNominalNonLinHighEt;
    delete fHistNominalNonLinLowEt;
    delete fHistNominalEffHighEt;
    delete fHistNominalEffLowEt;
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

    TRefArray *caloClusters = fSelector->GetClusters();
    Float_t fClusterMult = caloClusters->GetEntries();

    Float_t nominalRawEt = 0;
    Float_t nonlinHighRawEt = 0;
    Float_t nonlinLowRawEt = 0;
    Float_t effHighRawEt = 0;
    Float_t effLowRawEt = 0;

    Float_t nChargedHadronsMeasured = 0.0;
    Float_t nChargedHadronsTotal = 0.0;
    Float_t nChargedHadronsEtMeasured = 0.0;
    Float_t nChargedHadronsEtTotal = 0.0;


    for (Int_t iCluster = 0; iCluster < event->GetNumberOfCaloClusters(); iCluster++)
    {
        AliESDCaloCluster* cluster = event->GetCaloCluster(iCluster);
        if (!cluster)
        {
            AliError(Form("ERROR: Could not get cluster %d", iCluster));
            continue;
        }
        int x = 0;
	fCutFlow->Fill(x++);
	if(!fSelector->IsDetectorCluster(*cluster)) continue;
	fCutFlow->Fill(x++);
	if(!fSelector->PassMinEnergyCut(*cluster)) continue;
	fCutFlow->Fill(x++);
        if (!fSelector->PassDistanceToBadChannelCut(*cluster)) continue;
	fCutFlow->Fill(x++);

        Float_t pos[3];

        cluster->GetPosition(pos);
        TVector3 cp(pos);

        Bool_t matched = kTRUE;//default to no track matched
	Int_t trackMatchedIndex = cluster->GetTrackMatchedIndex();//find the index of the matched track
	matched = fSelector->PassTrackMatchingCut(*cluster);
	if(!matched){
	  if(trackMatchedIndex < 0) matched=kTRUE;
	  AliESDtrack *track = event->GetTrack(trackMatchedIndex);
	  //if this is a good track, accept track will return true.  The track matched is good, so not track matched is false
	  matched = !(fEsdtrackCutsTPC->AcceptTrack(track));
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
		  nChargedHadronsMeasured++;
		  nChargedHadronsTotal += 1/fTmCorrections->TrackMatchingEfficiency(track->Pt(),fClusterMult);
		  Double_t effCorrEt = CorrectForReconstructionEfficiency(*cluster);
		  nChargedHadronsEtMeasured+= TMath::Sin(cp.Theta())*effCorrEt;
		  nChargedHadronsEtTotal+= 1/fTmCorrections->TrackMatchingEfficiency(track->Pt(),fClusterMult) *effCorrEt;
		  fHistMatchedTracksEvspTvsMult->Fill(track->P(),TMath::Sin(cp.Theta())*cluster->E(),fClusterMult);
		  fHistMatchedTracksEvspTvsMultEffCorr->Fill(track->P(),CorrectForReconstructionEfficiency(*cluster),fClusterMult);
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
                            fEnergyDeposited = cluster->E();
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

                            //Float_t et = cluster->E() * TMath::Sin(theta);
                            if (maxpid == AliPID::kProton)
                            {

                                if (track->Charge() == 1)
                                {
                                    fHistProtonEnergyDeposit->Fill(cluster->E(), track->E());
                                }
                                else if (track->Charge() == -1)
                                {
                                    fHistAntiProtonEnergyDeposit->Fill(cluster->E(), track->E());
                                }
                            }
                            else if (maxpid == AliPID::kPion)
                            {
                                fHistChargedPionEnergyDeposit->Fill(cluster->E(), track->E());
                            }
                            else if (maxpid == AliPID::kKaon)
                            {
                                fHistChargedKaonEnergyDeposit->Fill(cluster->E(), track->E());
                            }
                            else if (maxpid == AliPID::kMuon)
                            {
                                fHistMuonEnergyDeposit->Fill(cluster->E(), track->E());
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
	  
	  //if (cluster->E() >  fSingleCellEnergyCut && cluster->GetNCells() == fCuts->GetCommonSingleCell()) continue;
	  //if (cluster->E() < fClusterEnergyCut) continue;
	  cluster->GetPosition(pos);
	  
	    TVector3 p2(pos);
	    
	    fClusterPosition->Fill(p2.Phi(), p2.PseudoRapidity());
	    fClusterEnergy->Fill(cluster->E());
	    fClusterEt->Fill(TMath::Sin(p2.Theta())*cluster->E());

	    Double_t effCorrEt = CorrectForReconstructionEfficiency(*cluster);
	    fTotNeutralEt += effCorrEt;
	    nominalRawEt += effCorrEt;
	    nonlinHighRawEt += effCorrEt*GetCorrectionModification(*cluster,1,0);
	    nonlinLowRawEt += effCorrEt*GetCorrectionModification(*cluster,-1,0);
	    effHighRawEt += effCorrEt*GetCorrectionModification(*cluster,0,1);
	    effLowRawEt += effCorrEt*GetCorrectionModification(*cluster,0,-1);
            fNeutralMultiplicity++;
        }
        fMultiplicity++;
    }
    
    fChargedEnergyRemoved = GetChargedContribution(fNeutralMultiplicity);
    fNeutralEnergyRemoved = GetNeutralContribution(fNeutralMultiplicity);
    fHistChargedEnergyRemoved->Fill(fChargedEnergyRemoved, fNeutralMultiplicity);
    fHistNeutralEnergyRemoved->Fill(fNeutralEnergyRemoved, fNeutralMultiplicity);
    
    fGammaEnergyAdded = GetGammaContribution(fNeutralMultiplicity);
    fHistGammaEnergyAdded->Fill(fGammaEnergyAdded, fNeutralMultiplicity);

    Double_t removedEnergy = GetChargedContribution(fNeutralMultiplicity) + GetNeutralContribution(fNeutralMultiplicity) + GetGammaContribution(fNeutralMultiplicity) + GetSecondaryContribution(fNeutralMultiplicity);
    fHistRemovedEnergy->Fill(removedEnergy);
    
    fTotNeutralEt = fGeomCorrection * fEMinCorrection * (fTotNeutralEt - removedEnergy);
    fTotEt = fTotChargedEt + fTotNeutralEt;
// Fill the histograms...0
    FillHistograms();
    //std::cout << "fTotNeutralEt: " << fTotNeutralEt << ", Contribution from non-removed charged: " << GetChargedContribution(fNeutralMultiplicity) << ", neutral: " << GetNeutralContribution(fNeutralMultiplicity) << ", gammas: " << GetGammaContribution(fNeutralMultiplicity) << ", multiplicity: " << fNeutralMultiplicity<< std::endl;
    //cout<<"cent "<<cent<<" cluster mult "<<fClusterMult<<" fTotNeutralEt "<<fTotNeutralEt<<" nominalRawEt "<<nominalRawEt<<endl;
    fHistNominalRawEt->Fill(nominalRawEt,cent);
    fHistNominalNonLinHighEt->Fill(nonlinHighRawEt,cent);
    fHistNominalNonLinLowEt->Fill(nonlinLowRawEt,cent);
    fHistNominalEffHighEt->Fill(effHighRawEt,cent);
    fHistNominalEffLowEt->Fill(effLowRawEt,cent);
    fHistFoundHadronsvsCent->Fill(nChargedHadronsMeasured,cent);
    fHistNotFoundHadronsvsCent->Fill(nChargedHadronsTotal-nChargedHadronsMeasured,cent);
    fHistFoundHadronsEtvsCent->Fill(nChargedHadronsEtMeasured,cent);
    fHistNotFoundHadronsEtvsCent->Fill(nChargedHadronsEtTotal-nChargedHadronsEtMeasured,cent);
//     cout<<"Number of hadrons measured:  "<<nChargedHadronsMeasured<<" Estimated total number of hadrons "<<nChargedHadronsTotal<<" ET in track matched hadrons "<<
//       nChargedHadronsEtMeasured;
//     if(nChargedHadronsMeasured>0)cout<<" ("<<nChargedHadronsEtMeasured/nChargedHadronsMeasured<<") ";
//     cout<<" ET in all hadrons ";
//     cout<<nChargedHadronsEtTotal;
//     if(nChargedHadronsTotal>0) cout<<" ("<<nChargedHadronsEtTotal/nChargedHadronsTotal<<") ";
//     cout<<endl;
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

    list->Add(fHistChargedPionEnergyDeposit);
    list->Add(fHistProtonEnergyDeposit);
    list->Add(fHistAntiProtonEnergyDeposit);
    list->Add(fHistChargedKaonEnergyDeposit);
    list->Add(fHistMuonEnergyDeposit);

    list->Add(fHistRemovedEnergy);
    list->Add(fClusterPosition);
    list->Add(fClusterEnergy);
    list->Add(fClusterEt);
    
    list->Add(fHistChargedEnergyRemoved);
    list->Add(fHistNeutralEnergyRemoved);
    list->Add(fHistGammaEnergyAdded);
    list->Add(fHistMatchedTracksEvspTvsMult);
    list->Add(fHistMatchedTracksEvspTvsMultEffCorr);
    list->Add(fHistFoundHadronsvsCent);
    list->Add(fHistNotFoundHadronsvsCent);
    list->Add(fHistFoundHadronsEtvsCent);
    list->Add(fHistNotFoundHadronsEtvsCent);
    list->Add(fHistNominalRawEt);
    list->Add(fHistNominalNonLinHighEt);
    list->Add(fHistNominalNonLinLowEt);
    list->Add(fHistNominalEffHighEt);
    list->Add(fHistNominalEffLowEt);
}

void AliAnalysisEtReconstructed::CreateHistograms()
{ // add some extra histograms to the ones from base class
    AliAnalysisEt::CreateHistograms();

    Int_t nbinsEt = 1000;
    Double_t minEt = 0;
    Double_t maxEt = 10;

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

    histname = "fClusterPosition" + fHistogramNameSuffix;
    fClusterPosition = new TH2D(histname.Data(), "Position of accepted neutral clusters",1000, -2.0, -.5, 1000, -.13 , 0.13);
    fClusterPosition->SetXTitle("Energy deposited in calorimeter");
    fClusterPosition->SetYTitle("Energy of track");

    histname = "fClusterEnergy" + fHistogramNameSuffix;
    fClusterEnergy = new TH1F(histname.Data(), histname.Data(), 100, 0, 5);
    fClusterEnergy->SetXTitle("Number of clusters");
    fClusterEnergy->SetYTitle("Energy of cluster");

    histname = "fClusterEt" + fHistogramNameSuffix;
    fClusterEt = new TH1F(histname.Data(), histname.Data(), 100, 0, 5);
    fClusterEt->SetXTitle("Number of clusters");
    fClusterEt->SetYTitle("E_{T} of cluster");

    histname = "fHistChargedEnergyRemoved" + fHistogramNameSuffix;
    fHistChargedEnergyRemoved = new TH2D(histname.Data(), histname.Data(), 1000, .0, 30, 100, -0.5 , 99.5);

    histname = "fHistNeutralEnergyRemoved" + fHistogramNameSuffix;
    fHistNeutralEnergyRemoved = new TH2D(histname.Data(), histname.Data(), 1000, .0, 30, 100, -0.5 , 99.5);

    histname = "fHistGammaEnergyAdded" + fHistogramNameSuffix;
    fHistGammaEnergyAdded = new TH2D(histname.Data(), histname.Data(), 1000, .0, 30, 100, -0.5 , 99.5);

    fHistMatchedTracksEvspTvsMult = new TH3F("fHistMatchedTracksEvspTvsMult", "fHistMatchedTracksEvspTvsMult",100, 0, 3,100,0,3,10,0,100);
    fHistMatchedTracksEvspTvsMultEffCorr = new TH3F("fHistMatchedTracksEvspTvsMultEffCorr", "fHistMatchedTracksEvspTvsMultEffCorr",100, 0, 3,100,0,3,10,0,100);
    fHistFoundHadronsvsCent = new TH2F("fHistFoundHadronsvsCent","fHistFoundHadronsvsCent",100,0,100,20,0,20);
    fHistNotFoundHadronsvsCent = new TH2F("fHistNotFoundHadronsvsCent","fHistNotFoundHadronsvsCent",100,0,100,20,0,20);
    fHistFoundHadronsEtvsCent = new TH2F("fHistFoundHadronsEtvsCent","fHistFoundHadronsEtvsCent",100,0,200,20,0,20);
    fHistNotFoundHadronsEtvsCent = new TH2F("fHistNotFoundHadronsEtvsCent","fHistNotFoundHadronsEtvsCent",100,0,200,20,0,20);
    
    maxEt = 100;
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

}
Double_t AliAnalysisEtReconstructed::ApplyModifiedCorrections(const AliESDCaloCluster& cluster,Int_t nonLinCorr, Int_t effCorr)
{
  Float_t pos[3];
  cluster.GetPosition(pos);
  TVector3 cp(pos);
  Double_t corrEnergy = fReCorrections->CorrectedEnergy(cluster.E());
  
  Double_t factorNonLin = GetCorrectionModification(cluster, nonLinCorr,effCorr);

  //std::cout << "Original energy: " << cluster.E() << ", corrected energy: " << corrEnergy << std::endl;
  return TMath::Sin(cp.Theta())*corrEnergy*factorNonLin;
}

Double_t AliAnalysisEtReconstructed::GetCorrectionModification(const AliESDCaloCluster& cluster,Int_t nonLinCorr, Int_t effCorr){//nonLinCorr 0 = nominal 1 = high -1 = low, effCorr  0 = nominal 1 = high -1 = low
  if(nonLinCorr==0){
    cout<<"Warning:  This function should not get called!"<<endl;//this statement is basically here to avoid a compilation warning
  }
  if(effCorr==0){
    cout<<"Warning:  This function should not get called!"<<endl;//this statement is basically here to avoid a compilation warning
  }
  return cluster.E();
}
