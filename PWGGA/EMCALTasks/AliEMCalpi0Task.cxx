// /**************************************************************************
/// * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
/// *             /home/amorreal/Alice/Work/taskEMCal/AliEMCalpi0Task.cxx
/// *
/// * Author: The ALICE Off-line Project.                                    *
/// * Contributors are mentioned in the code where appropriate.              *
/// *                                                                        *
/// * Permission to use, copy, modify and distribute this software and its   *
/// * documentation strictly for non-commercial purposes is hereby granted   *
/// * without fee, provided that the above copyright notice appears in all   *
/// * copies and that both the copyright notice and this permission notice   *
/// * appear in the supporting documentation. The authors make no claims     *
/// * about the suitability of this software for any purpose. It is          *
/// * provided"as is" without express or implied warranty.                  *
/// **************************************************************************/
/// \class AliEMCalpi0Task
/// AliAnalysisTaskSE to analyse data from AODs (pi0 and Eta in the EMCal).
/// The Task reads as input AOD (or ESDs) then it produces a .root file:
/// general Histos and Tree
///
/// \author Astrid Morreale
#include "AliEMCalpi0Task.h"
#include "Alipi0EventStatStruct.h"
#include "Alipi0ClusterStatStruct.h"
#include "Alipi0DiClusterStatStruct.h"
#include "Alipi0mixedDiClusterStatStruct.h"
#include <Riostream.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TList.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TRandom3.h>


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliEMCALGeometry.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliCentrality.h"
#include "AliEMCALRecoUtils.h"
#include "AliExternalTrackParam.h"


/// ROOT includes
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <TH2F.h>
#include <TArrayI.h>
#include <TArrayF.h>
#include <TObjArray.h>

/// STEER includes
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliExternalTrackParam.h"
#include "AliESDfriendTrack.h"
#include "AliTrackerBase.h"


// EMCAL includes
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliTrackerBase.h"
#include "AliEMCALCalibTimeDepCorrection.h"
#include "AliEMCALPIDUtils.h"

#include <fstream>
#include <cassert>

#include <TChain.h>
#include <TError.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TLorentzVector.h>

#include <Riostream.h>

///\cond CLASSIMP
ClassImp(AliEMCalpi0Task)
///\endcond

//________________________________________________________________________
AliEMCalpi0Task::AliEMCalpi0Task(const char *name) :
AliAnalysisTaskSE(name),
fEvent(0),
fEventStatStruct(0x0),
fClusterStatArray( 0x0 ),
fdiClusterStatArray( 0x0 ),
fmixedDiClusterStatArray(0x0),
fClusterStatCount( 0 ),
fdiClusterStatCount( 0 ),
fmixedDiClusterStatCount( 0 )

{

    for (Int_t i=0; i < kNtype; i++)
    {
        fPool[i]             = 0;

    }
    /// Define input and output slots here
    /// Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliEMCalpi0Task::~AliEMCalpi0Task()
{

    /// destructor
    ///only if they are initialized to 0, or to something valid
    delete fEventStatStruct;
    delete fClusterStatArray;
    delete fdiClusterStatArray;
    delete fmixedDiClusterStatArray;
    /// delete the 'running' pool array, which is either zero,
    /// or corresponds with the last created pool during event processing.
    for (Int_t i=0; i < kNtype; i++)
    { delete fPool[i]; }

}

//________________________________________________________________________
void AliEMCalpi0Task::UserCreateOutputObjects( void )
{
    /// get AOD output handler
    AliAODHandler* handler = dynamic_cast<AliAODHandler*>( AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler() );
    if( handler )
    {



        fEventStatStruct = new EventStatStruct();
        handler->AddBranch("EventStatStruct", &fEventStatStruct);

        ///single clusters
        fClusterStatArray = new TClonesArray( "ClusterStatStruct", 0 );
        fClusterStatArray->SetName( "clusterStatArray" );
        handler->AddBranch("TClonesArray", &fClusterStatArray);
        fClusterStatCount = 0;

        ///di-clusters
        fdiClusterStatArray = new TClonesArray( "diClusterStatStruct", 0 );
        fdiClusterStatArray->SetName( "diClusterStatArray" );
        handler->AddBranch("TClonesArray", &fdiClusterStatArray);
        fdiClusterStatCount = 0;

        ///di-clusters mixed
        fmixedDiClusterStatArray = new TClonesArray( "mixedDiClusterStatStruct", 0 );
        fmixedDiClusterStatArray->SetName( "mixedDiClusterStatArray" );
        handler->AddBranch("TClonesArray", &fmixedDiClusterStatArray);
        fmixedDiClusterStatCount = 0;


    }
}



///________________________________________________________________________
void AliEMCalpi0Task::UserExec( Option_t* )
{

    /// increase event number
    AliInfo( Form( "event: %i", fEvent++ ) );
    /// ESD Filter analysis task executed for each event

    AliESDEvent* esd    = dynamic_cast<AliESDEvent*>(InputEvent());
    AliAODEvent* aod    = dynamic_cast< AliAODEvent*>(InputEvent());
    AliVEvent  * event  = InputEvent();


    fClusterStatArray->Clear();
    fClusterStatCount = 0;

    fdiClusterStatArray->Clear();
    fdiClusterStatCount = 0;

    fmixedDiClusterStatArray->Clear();
    fmixedDiClusterStatCount = 0;

    fEventStatStruct->Clear();


    if(esd){
        if( esd->IsPileupFromSPD(3,0.8) ) return;
    }


    ///Remove events with exotic clusters.
    ///Not used for MB, done at tender level if necessary.
    ///Can also be done at analysis level since Exo variable is saved.
    TRefArray *caloClusArr=new TRefArray();
    event->GetEMCALClusters(caloClusArr);
    const Int_t kNumber =caloClusArr->GetEntries();

    for( Int_t iclus = 0; iclus < kNumber; iclus++ ){

        AliVCluster*c=(AliVCluster*) caloClusArr->At(iclus);
        if(!c){
            return;
        }
        if(!c->IsEMCAL()){
            return;
        }

    }

    const UInt_t Mask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

    /// event characterization
    if(esd){
        fEventStatStruct->runNumber     = esd->GetRunNumber();
        TString triggerClasses = esd->GetFiredTriggerClasses();
    }

    if(aod){
        fEventStatStruct->runNumber     = aod->GetRunNumber();
        /// Get triggered classes, bordel....
        TString triggerClasses = aod->GetFiredTriggerClasses();
    }

    /// ULong64_t triggerMask  = event->GetTriggerMask();
    /// verification triggered classes that fired.

    fEventStatStruct->isMB           = Mask & AliVEvent::kMB;
    fEventStatStruct->isAnyINT       = Mask & AliVEvent::kAnyINT;
    fEventStatStruct->isCentral      = Mask & AliVEvent::kCentral;
    fEventStatStruct->isSemiCentral  = Mask & AliVEvent::kSemiCentral;
    fEventStatStruct->isEga          = Mask & AliVEvent::kEMCEGA;
    fEventStatStruct->kAllMB = (
        fEventStatStruct->isMB ||
        fEventStatStruct->isCentral ||
        fEventStatStruct->isSemiCentral );




    AliCentrality *esdCent;
    AliCentrality *aodCent;
    /// Centrality
    if(esd){
        esdCent     = esd->GetCentrality();
        fEventStatStruct->CentralityVZERO  = esdCent->GetCentralityPercentile("V0M");
        fEventStatStruct->CentralitySPD    = esdCent->GetCentralityPercentile("CL1");
    }

    if(aod){
        aodCent     = aod->GetCentrality();
        fEventStatStruct-> CentralityVZERO  = aodCent->GetCentralityPercentile("V0M");
        fEventStatStruct-> CentralitySPD    = aodCent->GetCentralityPercentile("CL1");

        Int_t nTracks=aod->GetNumberOfTracks();
        fEventStatStruct->multiplicity =nTracks;

        AliInfo( Form( "===NTRACKS===:%i",fEventStatStruct->multiplicity ) );
    }


    ///Pass the geometry transformation matrix from ESDs to geometry
    AliEMCALGeometry  *geom        =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
    AliVCaloCells     *Cells       =  event->GetEMCALCells();
    AliEMCALRecoUtils RecoUtils;

    Int_t    absID1    = -1;        Int_t    absID2 = -1;
    Int_t    ieta1     = -1;        Int_t    ieta2  = -1;
    Int_t    iphi1     = -1;        Int_t    iphi2  = -1;
    Int_t    iSM1      = -1;        Int_t    iSM2   = -1;
    Bool_t  shared1;      Bool_t  shared2;

    ///get reconstructed vertex position
    Double_t vertex_position[3];
    if(esd)
    {
        esd->GetVertex()->GetXYZ(vertex_position);
        fEventStatStruct->  vX     = vertex_position[0];
        fEventStatStruct->  vY     = vertex_position[1];
        fEventStatStruct->  vZ     = vertex_position[2];
    }

    if(aod){

        AliAODVertex* Vtx = aod->GetPrimaryVertex();
        fEventStatStruct->  nV     = Vtx->GetNContributors();
        fEventStatStruct->  vX     = Vtx->GetX();
        fEventStatStruct->  vY     = Vtx->GetY();
        fEventStatStruct->  vZ     = Vtx->GetZ();
        AliInfo( Form( "nV: %i, vZ:  %.3f", Vtx->GetNContributors(), Vtx->GetZ()) );

    }


    /// Cut on vertex at analysis level.
    /// if(vZ>-15.0 && vZ<15.0){//vertex cut
    ///array temporaire pour passe plus tard dans le boucles
    TRefArray caloClustersArr = TRefArray();
    event->GetEMCALClusters( &caloClustersArr );

    const Int_t kNumberOfEMCALClusters =caloClustersArr.GetEntries();

    /// pool to store clusters to be mixed
    /// it is not deleted directly. Instead, it is assigned to fPool, the previously assigned pool
    /// being deleted (if any), at that time.
    TObjArray *newPool = new TObjArray(kNumberOfEMCALClusters);
    newPool->SetOwner();

    TVector3 pos;
    pos -= vertex_position;
    ///Double_t r1 = pos.Mag();

    ///Fill the pool with all clusters
    Int_t nGoodClusters = 0;
    /*      const Float_t etaCut = 0.65;
    const Float_t phiCutMin = 1.3962634;
    const Float_t phiCutMax = 3.14159265359;
    */

    ///boucle sur tous les clusters

    for( Int_t iclus = 0; iclus < kNumberOfEMCALClusters; iclus++ )
    {

        AliVCluster*c1=(AliVCluster*) caloClustersArr.At(iclus);
        if (!c1) continue;
        if (!c1->IsEMCAL()) continue;



        Int_t bc;
        if(esd) bc = esd->GetBunchCrossNumber();
        if(aod) bc = aod->GetBunchCrossNumber();

        //AliVCaloCells     *Cells       =  event->GetEMCALCells();
        //AliEMCALRecoUtils RecoUtils;
        //for(Int_t icell=0;icell<Cells->GetNumberOfCells();icell++){
        //===============================
        if(c1->GetNCells()<1) continue;
        if(c1->E()<2) continue;
        //if(c1->GetM02()>1.0) continue;
        if(c1->GetDistanceToBadChannel()<1) continue;
        //===============================
        TLorentzVector pii;
        c1->GetMomentum(pii, vertex_position);
        TLorentzVector *pimx = new TLorentzVector;
        c1->GetMomentum(*pimx, vertex_position);

        // add TLorentz vector in the pool array
        // it will be deleted when the array is deleted
        newPool->Add(pimx);
        ++nGoodClusters;

        //create empty clusterStat object
        ClusterStatStruct clusterStat;

        Int_t id= -1;;
        //id = TMath::Abs(Cells->GetCellNumber(icell));
        Float_t Emax    = GetMaxCellEnergy( c1, id);
        //AliEMCALGeometry  *geom        =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
        Float_t tcell=Cells->GetCellTime(id);
        Float_t Ecross =  RecoUtils.GetECross(id, tcell,Cells, bc);
        Float_t Exo    = 1.0 - Ecross/Emax;

        clusterStat.EMax    = Emax;
        clusterStat.Tcell  = tcell;
        clusterStat.ECross = Ecross;
        clusterStat.ExoC   = Exo;


        AliInfo( Form( "TIME CELL %f",Exo) );

        //characteristiques cluster
        clusterStat.NCells                 = c1->GetNCells();
        clusterStat.distanceBadChannel     = c1->GetDistanceToBadChannel();
        clusterStat.Ecluster               = c1->E();
        clusterStat.NCellscluster          = c1->GetNCells();
        clusterStat.M20cluster             = c1->GetM20();
        clusterStat.M02cluster             = c1->GetM02();
        clusterStat.NCluscluster           = kNumberOfEMCALClusters;
        clusterStat.isEMCALcluster         = c1->IsEMCAL();
        clusterStat.dispersioncluster      = c1->GetDispersion();
        clusterStat.chi2cluster            = c1->Chi2();
        clusterStat.distBadChannelcluster  = c1->GetDistanceToBadChannel();
        clusterStat.phicluster             = pii.Phi();
        clusterStat.etacluster             = pii.Eta();
        clusterStat.ptcluster              = pii.Pt();
        RecoUtils.GetMaxEnergyCell(geom, Cells, c1, absID1, iSM1, ieta1, iphi1, shared1);

        Int_t trackindex=c1->GetNTracksMatched();
        if(trackindex>0)    clusterStat.TrackFlag1 =1;

        AliInfo( Form( "-----------matched Track: %i",  trackindex) );



        AliInfo( Form( "pT : %f", pii.Pt()) );
        for (Int_t iclus2 = iclus+1;  iclus2< kNumberOfEMCALClusters; iclus2++)
        {

            AliVCluster *c2 = (AliVCluster*) caloClustersArr.At(iclus2);
            if (!c2) continue;
            if (!c2->IsEMCAL()) continue;

            //======================
            if (c2->GetNCells()<1) continue;
            if (c2->E()<2) continue;
            //if (c2->GetM02()>1.0) continue;
            if (c2->GetDistanceToBadChannel()<1) continue;
            //======================
            Int_t trackindex2=c2->GetNTracksMatched();
            if(trackindex2>0)    clusterStat.TrackFlag2 =1;

            // Float_t en2 = c2->E();
            ///create empty diclusterStat object
            diClusterStatStruct diclusterStat;

            ///Set flags for studies
            if(c1->GetNCells()>2 && c2->GetNCells()>2)
            {
                diclusterStat.NCellFlag =kTRUE;
            }
            if(c1->E()>2 && c2->E()>2)
            {
                diclusterStat.EcFlag = kTRUE;
            }
            if(c1->GetM02()<0.5 && c2->GetM02()<0.5)
            {
                diclusterStat.M02Flag   = kTRUE;
            }
            if(c1->GetDistanceToBadChannel()>2 &&c2->GetDistanceToBadChannel()>2)
            {
                diclusterStat.D2BadChFlag = kTRUE;
            }


            TLorentzVector pjj;
            c2->GetMomentum(pjj, vertex_position);
            TLorentzVector pion;
            pion   = pii + pjj;
            ///remplissage des pions
            diclusterStat.piE      = pion.E();
            diclusterStat.piphi    = pion.Phi();
            diclusterStat.pieta    = pion.Eta();
            diclusterStat.ptpi     = pion.Pt();
            diclusterStat.pipx     = pion.Px();
            diclusterStat.pipy     = pion.Py();
            diclusterStat.pipz     = pion.Pz();
            diclusterStat.asympi   = TMath::Abs(pii.E()-pjj.E())/(pii.E()+pjj.E());
            diclusterStat.masspi   = pion.M();
            RecoUtils.GetMaxEnergyCell(geom, Cells, c2, absID2, iSM2,  ieta2,  iphi2, shared2);
            new((*fdiClusterStatArray)[fdiClusterStatCount]) diClusterStatStruct( diclusterStat );
            fdiClusterStatCount++;
        }
        /// copy clusterStat into TClonesArray, at last position
        new((*fClusterStatArray)[fClusterStatCount]) ClusterStatStruct( clusterStat );
        fClusterStatCount++;

    }


    ///fill mixed event, can be also done at analysis level for more segmentation options
    int poolIndex = 0;
    if( fEventStatStruct->CentralityVZERO>0  && fEventStatStruct->CentralityVZERO<=10 ) poolIndex =0;
    if( fEventStatStruct->CentralityVZERO>10 && fEventStatStruct->CentralityVZERO<=20 ) poolIndex =1;
    if( fEventStatStruct->CentralityVZERO>20 && fEventStatStruct->CentralityVZERO<=30 ) poolIndex =2;
    if( fEventStatStruct->CentralityVZERO>20 && fEventStatStruct->CentralityVZERO<=50 ) poolIndex =3;
    if( fEventStatStruct->CentralityVZERO>30 && fEventStatStruct->CentralityVZERO<=50 ) poolIndex =4;
    if( fEventStatStruct->CentralityVZERO>50 && fEventStatStruct->CentralityVZERO<=90 ) poolIndex =5;
    if( fEventStatStruct->CentralityVZERO>0  && fEventStatStruct->CentralityVZERO<=20 ) poolIndex =6;
    if( fEventStatStruct->CentralityVZERO>20 && fEventStatStruct->CentralityVZERO<=40 ) poolIndex =7;
    if( fEventStatStruct->CentralityVZERO>40 && fEventStatStruct->CentralityVZERO<=60 ) poolIndex =8;
    if( fEventStatStruct->CentralityVZERO>60 && fEventStatStruct->CentralityVZERO<=90 ) poolIndex =9;




    if( !fPool[poolIndex] )
    {
        fPool[poolIndex] = newPool;

    }
    else {

        Int_t nGoodClusters2 = fPool[poolIndex]->GetEntriesFast();
        //GetEntries();

        for (Int_t i=0; i<nGoodClusters; ++i)
        {

            TLorentzVector * pi = static_cast<TLorentzVector*>(newPool->At(i));
            for (Int_t j=0; j<nGoodClusters2; ++j)
            {

                TLorentzVector *pj= static_cast<TLorentzVector*>(fPool[poolIndex]->At(j));
                FillMixed(*pi,*pj);

            } //random cluster
        } //current cluster

        /// delete previous pool and assign to new pool.
        delete fPool[poolIndex];
        fPool[poolIndex] = newPool;

    }

    //}//15 cm vertex cut

    /// save AOD
    AliAODHandler* handler = dynamic_cast<AliAODHandler*>( AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler() );
    if( handler )
    {

        //AliAODEvent* aod = handler->GetAOD();
        //AliAODHeader* header = aod->GetHeader();
        //header->SetRunNumber(aod->GetRunNumber());
        AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
    } else AliInfo( "Error: invalid output event handler" );

}//end process



//________________________________________________________________________
Double_t AliEMCalpi0Task ::GetMaxCellEnergy(const AliVCluster *cluster, Int_t &id) const
{
    /// Get maximum energy of attached cell.
    AliESDEvent* esd  = dynamic_cast<AliESDEvent*>(InputEvent());
    AliAODEvent* aod =dynamic_cast< AliAODEvent*>(InputEvent());

    // AliEMCALGeometry  *fGeom        =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");

    id = -1;

    AliVCaloCells *cells = 0;
    if(esd){
        cells = esd->GetEMCALCells();
    }

    if(!esd){
        cells = aod->GetEMCALCells();
    }

    if (!cells)
        return 0;

    Double_t maxe = 0;
    Int_t ncells = cluster->GetNCells();
    for (Int_t i=0; i<ncells; i++) {
        Double_t e = cells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
        if (e>maxe) {
            maxe = e;
            id   = cluster->GetCellAbsId(i);
        }
    }
    return maxe;
}
//________________________________________________________________
void AliEMCalpi0Task::FillMixed( const TLorentzVector& p1, const TLorentzVector& p2)

{
    ///verification triggered classes that fired.
    //not used but may be useful for documentation sake
    const UInt_t eventSelectionMask( ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() );
    AliESDEvent* esd     =dynamic_cast<AliESDEvent*>(InputEvent());
    AliAODEvent* aod     =dynamic_cast< AliAODEvent*>(InputEvent());
    //AliVEvent  * event   = InputEvent();

    ///create empty mixedDiClusterStat object
    mixedDiClusterStatStruct  mixedDiClusterStat;

    // ULong64_t triggerMask = event->GetTriggerMask();

    mixedDiClusterStat.isMBmx          = (eventSelectionMask&AliVEvent::kMB);
    mixedDiClusterStat.isAnyINTmx      = (eventSelectionMask&AliVEvent::kAnyINT);
    mixedDiClusterStat.isCentralmx     = (eventSelectionMask&AliVEvent::kCentral);
    mixedDiClusterStat.isSemiCentralmx = (eventSelectionMask&AliVEvent::kSemiCentral);
    mixedDiClusterStat.isEgamx         = (eventSelectionMask&AliVEvent::kEMCEGA);

    mixedDiClusterStat.kAllMBmx= (
        mixedDiClusterStat.isMBmx ||
        mixedDiClusterStat.isAnyINTmx||
        mixedDiClusterStat.isCentralmx||
        mixedDiClusterStat.isSemiCentralmx );


    /// Centrality
    AliCentrality *esdCent;
    if(esd){
        esdCent     = esd->GetCentrality();
    }
    if(!esd){
        esdCent     = aod->GetCentrality();
    }


    mixedDiClusterStat.centMixedV0  = esdCent->GetCentralityPercentile("V0M");
    mixedDiClusterStat.centMixedSPD = esdCent->GetCentralityPercentile("CL1");

    TLorentzVector Mxpion;
    Mxpion = p1 + p2;

    mixedDiClusterStat.Mxmass      = Mxpion.M();
    mixedDiClusterStat.Mxpt        = Mxpion.Pt();

    new((*fmixedDiClusterStatArray)[fmixedDiClusterStatCount]) mixedDiClusterStatStruct( mixedDiClusterStat );
    fmixedDiClusterStatCount++;
}

//________________________________________________________________________
void AliEMCalpi0Task::Terminate(const Option_t*)
{}
