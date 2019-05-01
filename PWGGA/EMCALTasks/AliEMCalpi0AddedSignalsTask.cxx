/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* AliEMCalpi0AddedSignalsTask.cxx
*
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided"as is" without express or implied warranty.                  *
**************************************************************************/

//-----------------------------------------------------------------------------
/// \class AliEMCalpi0AddedSignalsTask
/// AliAnalysisTaskSE to analyse data from ESDs (pi0 and Eta in the EMCal).
/// The Task reads as input ESDs then it produces a .root file:
/// general Histos and Ntuple:Tree
///
/// \author Astrid Morreale
// Generator names:
//LHC14a1a: "hijing_0":"pi0_1":"eta_2":"pi0EMC_3":"pi0PHS_4":"etaEMC_5"
//LHC14a1 b and c: "hijing_0":"BOX":"PARAM":"Pythia":
//-----------------------------------------------------------------------------
#include "AliEMCalpi0AddedSignalsTask.h"
#include "Alipi0EventStatStruct.h"
#include "AlietaClusterMCStatStruct.h"
#include "Alipi0ClusterMCStatStruct.h"
#include "AliHijingClusterMCStatStruct.h"
#include "Alipi0ClusterStatStruct.h"
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
#include "AliMultiplicity.h"
#include "AliCentrality.h"
#include "AliEMCALRecoUtils.h"
#include "AliExternalTrackParam.h"

//Generator includes
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include <AliGenPythiaEventHeader.h>
#include <AliGenCocktailEventHeader.h>


// ROOT includes
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <TH2F.h>
#include <TArrayI.h>
#include <TArrayF.h>
#include <TObjArray.h>

// STEER includes
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
#include <set>

///\cond CLASSIMP
ClassImp(AliEMCalpi0AddedSignalsTask)
///\endcond

//________________________________________________________________________
AliEMCalpi0AddedSignalsTask::AliEMCalpi0AddedSignalsTask(const char *name) :
    AliAnalysisTaskSE(name),
    fEvent(0),
    fEventStatStruct(0x0),
    fClusterStatArray( 0x0 ),
    fClusterMCStatArray( 0x0 ),
    fClusterpiMCStatArray( 0x0 ),
    fClusterHijingMCStatArray( 0x0 ),
    fClusterStatCount( 0 ),
    fClusterMCStatCount( 0 ),
    fClusterpiMCStatCount( 0 ),
    fClusterHijingMCStatCount( 0 )
{
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliEMCalpi0AddedSignalsTask::~AliEMCalpi0AddedSignalsTask()
{

    /// destructor
    //only if they are initialized to 0, or to something valid
    delete fEventStatStruct;
    delete fClusterMCStatArray;
    delete fClusterpiMCStatArray;
    delete fClusterHijingMCStatArray;
    delete fClusterStatArray;

}

//________________________________________________________________________
void AliEMCalpi0AddedSignalsTask::UserCreateOutputObjects( void )
{
    // get AOD output handler
    AliAODHandler* handler = dynamic_cast<AliAODHandler*>( AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler() );
    if( handler )
    {

        //Add branch that access triggers, no array is needed since its one entry per event
        //  delete fEventStatStruct;
        fEventStatStruct = new EventStatStruct();
        handler->AddBranch("EventStatStruct", &fEventStatStruct);

        //MC information
        fClusterMCStatArray = new TClonesArray( "ClusterMCStatStruct", 0 );
        fClusterMCStatArray->SetName( "clusterMCStatArray" );
        handler->AddBranch("TClonesArray", &fClusterMCStatArray);
        fClusterMCStatCount = 0;

        //MC eta information
        fClusterpiMCStatArray = new TClonesArray( "ClusterpiMCStatStruct", 0 );
        fClusterpiMCStatArray->SetName( "clusterpiMCStatArray" );
        handler->AddBranch("TClonesArray", &fClusterpiMCStatArray);
        fClusterpiMCStatCount = 0;

        //MC hijing information
        fClusterHijingMCStatArray = new TClonesArray( "ClusterHijingMCStatStruct", 0 );
        fClusterHijingMCStatArray->SetName( "clusterHijingMCStatArray" );
        handler->AddBranch("TClonesArray", &fClusterHijingMCStatArray);
        fClusterHijingMCStatCount = 0;

        //single clusters
        fClusterStatArray = new TClonesArray( "ClusterStatStruct", 0 );
        fClusterStatArray->SetName( "clusterStatArray" );
        handler->AddBranch("TClonesArray", &fClusterStatArray);
        fClusterStatCount = 0;

    }
}

//________________________________________________________________________
void AliEMCalpi0AddedSignalsTask::UserExec( Option_t* )
{

    // increase event number
    AliInfo( Form( "event: %i", fEvent++ ) );
    // ESD Filter analysis task executed for each event

    AliESDEvent* esd        = dynamic_cast<AliESDEvent*>(InputEvent());
    AliAODEvent* aod        = dynamic_cast< AliAODEvent*>(InputEvent());
    AliVEvent  * event      = InputEvent();
    const Float_t etaCut    = 0.65;
    const Float_t phiCutMin = 1.3962634;
    const Float_t phiCutMax = 3.14159265359;
    // get the MC event
    AliMCEvent* mcEvent( MCEvent() );
    if( !mcEvent )
    {
        AliError( "Invalid MC event. Skipping." );
        return;
    }


    // get event header
    AliGenEventHeader *eventHeader = mcEvent->GenEventHeader();
    AliGenCocktailEventHeader *cocktailEventHeader = dynamic_cast<AliGenCocktailEventHeader *>(eventHeader);
    if( !cocktailEventHeader )
    {
        AliError( "Invalid MC event header. Skipping" );
        return ;
    }

    ///LHC14a1a: "hijing_0":"pi0_1":"eta_2":"pi0EMC_3":"pi0PHS_4":"etaEMC_5"
    // ref generator name
    //
    const TString refGeneratorName( "PARAM" );
    const TString refGeneratorName2( "Hijing" );

    //==========
    //-----------
    const TString HijingGeneratorName( "hijing_0");
    const TString pionGeneratorName( "pi0EMC_3");
    const TString etaGeneratorName( "etaEMC_5");
    //---------
    //==========

    fClusterMCStatArray->Clear();
    fClusterMCStatCount = 0;

    fClusterpiMCStatArray->Clear();
    fClusterpiMCStatCount = 0;

    fClusterHijingMCStatArray->Clear();
    fClusterHijingMCStatCount = 0;

    fClusterStatArray->Clear();
    fClusterStatCount = 0;



    fEventStatStruct->Clear();

    // first and last particle index
    Int_t firstParticle = 0;
    // Int_t lastParticle = 0;

    Int_t firstHijingParticle = 0;
    Int_t lastHijingParticle = 0;

    Int_t firstpiParticle = 0;
    Int_t lastpiParticle = 0;

    Int_t firstetaParticle = 0;
    Int_t lastetaParticle = 0;


    // loop over generators
    TList *genHeaders = cocktailEventHeader->GetHeaders();
    const Int_t nGenerators( genHeaders->GetEntries() );
    for(Int_t iGen = 0; iGen < nGenerators; ++iGen )
    {

        AliGenEventHeader * genHeader = (AliGenEventHeader*)genHeaders->At( iGen );
        const TString generatorName = genHeader->GetName();

        AliInfo( Form( "Generator index %i, name :\"%s\"", iGen, genHeader->GetName()));

        //0==Hijing, pi0==3, eta ==5
        //else Hijing, PARAM
        //if( generatorName == (HijingGeneratorName||refGeneratorName2) && iGen==0)
        if(( (generatorName == HijingGeneratorName)||(generatorName ==refGeneratorName2) ) && iGen==0)
        {

            firstHijingParticle=firstParticle;
            lastHijingParticle = firstHijingParticle + genHeader->NProduced();

        }
        else if (( (generatorName == pionGeneratorName)||(generatorName == refGeneratorName) )&& iGen==3)
        {
            firstpiParticle=firstParticle;
            lastpiParticle = firstpiParticle + genHeader->NProduced();


        }

        else if(( (generatorName == etaGeneratorName)||(generatorName ==refGeneratorName) )&& iGen==5)
        {
            firstetaParticle=firstParticle;
            lastetaParticle = firstetaParticle + genHeader->NProduced();


        }

        firstParticle += genHeader->NProduced();


    }

    if( firstHijingParticle!=lastHijingParticle )
    {

        //        AliInfo( Form( " generator named \"%s\", First particle: %i, last: %i", HijingGeneratorName.Data(), firstHijingParticle, lastHijingParticle ) );
        AliInfo( Form( " generator named \"%s\", First particle: %i, last: %i", refGeneratorName2.Data(), firstHijingParticle, lastHijingParticle ) );

    }

    else {
        //AliError( Form( "No generator named \"%s\" found in event. Skipping", HijingGeneratorName.Data() ) );
        AliError( Form( "No generator named \"%s\" found in event. Skipping", refGeneratorName2.Data() ) );
        return;
    }

    if( firstpiParticle!=lastpiParticle )
    {

        // AliInfo( Form( " generator named \"%s\", First particle: %i, last: %i", pionGeneratorName.Data(), firstpiParticle, lastpiParticle ) );
        AliInfo( Form( " generator named \"%s\", First particle: %i, last: %i", refGeneratorName.Data(), firstpiParticle, lastpiParticle ) );

    }
    else {
        // AliError( Form( "No generator named \"%s\" found in event. Skipping", pionGeneratorName.Data() ) );
        AliError( Form( "No generator named \"%s\" found in event. Skipping", refGeneratorName.Data() ) );
        return;

    }


    if( firstetaParticle!=lastetaParticle )
    {

        //  AliInfo( Form( " generator named \"%s\", First particle: %i, last: %i", etaGeneratorName.Data(), firstetaParticle, lastetaParticle ) );
        AliInfo( Form( " generator named \"%s\", First particle: %i, last: %i", refGeneratorName.Data(), firstetaParticle, lastetaParticle ) );

    }


    else {
        //AliError( Form( "No generator named \"%s\" found in event. Skipping", etaGeneratorName.Data() ) );
        AliError( Form( "No generator named \"%s\" found in event. Skipping", refGeneratorName.Data() ) );
        return;

    }

    // get particle stack
    AliStack* stack = mcEvent->Stack();
    if( !stack )
    {
        AliError( "Invalid particle stack. Skipping." );
        return;
    }

    //create empty clusterStat object
    ClusterMCStatStruct       clusterMCStat;
    ClusterpiMCStatStruct     clusterpiMCStat;
    ClusterHijingMCStatStruct clusterHijingMCStat;


    Int_t totatetas=0;
    //==========================================================================
    //Start with Etas only
    //===========================================================================
    // print relevant particles
    for( Int_t iPart = firstetaParticle; iPart < lastetaParticle; iPart++ )
    {
        TParticle* particle( stack->Particle(iPart) );

        if(
            particle->GetPdgCode()==221
            &&fabs(particle->Eta()) <=etaCut
            &&particle->GetNDaughters()==2
            &&(particle->Phi() >phiCutMin|| particle->Phi()< phiCutMax)
            )
        {

            TParticle* photon1( stack->Particle(particle->GetFirstDaughter()));
            TParticle* photon2( stack->Particle(particle->GetLastDaughter()));
            totatetas++;
            Int_t piddaughter1     =photon1->GetPdgCode();
            Int_t piddaughter2     =photon2->GetPdgCode();

            Int_t pid              =particle->GetPdgCode();
            // Int_t pidmother        =particle->GetMother(0);

            clusterMCStat.daughter1=piddaughter1;
            clusterMCStat.daughter2=piddaughter2;
            clusterMCStat.EnergyMC =particle->Energy();
            clusterMCStat.pTMC     =particle->Pt();
            clusterMCStat.pidMC    =pid;
            clusterMCStat.parentID =iPart;
            clusterMCStat.etaMC    =particle->Eta();
            clusterMCStat.phiMC    =particle->Phi();
            clusterMCStat.thetaMC  =particle->Theta();
            clusterMCStat.vXmc     =particle->Vx();
            clusterMCStat.vYmc     =particle->Vy();
            clusterMCStat.vZmc     =particle->Vz();

            new((*fClusterMCStatArray)[fClusterMCStatCount]) ClusterMCStatStruct( clusterMCStat );
            fClusterMCStatCount++;
        }

    }

    //==========================================================================
    //pi0s only
    //===========================================================================
    // print relevant particles
    for( Int_t iPart = firstpiParticle; iPart < lastpiParticle; iPart++ )
    {
        TParticle* particle( stack->Particle(iPart) );

        if(
            particle->GetPdgCode()==111
            &&fabs(particle->Eta()) <=etaCut
            &&particle->GetNDaughters()==2
            &&(particle->Phi() >phiCutMin|| particle->Phi()< phiCutMax))
        {

            TParticle* photon1( stack->Particle(particle->GetFirstDaughter()));
            TParticle* photon2( stack->Particle(particle->GetLastDaughter()));
            totatetas++;
            Int_t piddaughter1     =photon1->GetPdgCode();
            Int_t piddaughter2     =photon2->GetPdgCode();

            Int_t pid              =particle->GetPdgCode();
            // Int_t pidmother        =particle->GetMother(0);
            clusterpiMCStat.daughter1=piddaughter1;
            clusterpiMCStat.daughter2=piddaughter2;

            clusterpiMCStat.EnergyMC =particle->Energy();
            clusterpiMCStat.pTMC     =particle->Pt();
            clusterpiMCStat.pidMC    =pid;
            clusterpiMCStat.parentID =iPart;
            clusterpiMCStat.etaMC    =particle->Eta();
            clusterpiMCStat.phiMC    =particle->Phi();
            clusterpiMCStat.thetaMC  =particle->Theta();
            clusterpiMCStat.vXmc     =particle->Vx();
            clusterpiMCStat.vYmc     =particle->Vy();
            clusterpiMCStat.vZmc     =particle->Vz();

            new((*fClusterpiMCStatArray)[fClusterpiMCStatCount]) ClusterpiMCStatStruct( clusterpiMCStat );
            fClusterpiMCStatCount++;
        }

    }

   //__________________________________________
   //Do hijing only
   //__________________________________________
   // print relevant particles
    for( Int_t iPart = firstHijingParticle; iPart < lastHijingParticle; iPart++ )
    {
        TParticle* particle( stack->Particle(iPart) );

        if(
            fabs(particle->Eta()) <=etaCut
            &&(particle->Phi() >phiCutMin|| particle->Phi()< phiCutMax))
        {

            Int_t piddaughter1=-1;
            Int_t piddaughter2=-1;

            if(particle->GetNDaughters()>=2)
            {
                TParticle* photon1( stack->Particle(particle->GetFirstDaughter()));
                TParticle* photon2( stack->Particle(particle->GetLastDaughter()));
                piddaughter1             =photon1->GetPdgCode();
                piddaughter2             =photon2->GetPdgCode();
            }

            Int_t pid                      =particle->GetPdgCode();
            // Int_t pidmother                =particle->GetMother(0);
            clusterHijingMCStat.daughter1=piddaughter1;
            clusterHijingMCStat.daughter2=piddaughter2;

            clusterHijingMCStat.EnergyMC =particle->Energy();
            clusterHijingMCStat.pTMC     =particle->Pt();
            clusterHijingMCStat.pidMC    =pid;
            clusterHijingMCStat.parentID =iPart;
            clusterHijingMCStat.etaMC    =particle->Eta();
            clusterHijingMCStat.phiMC    =particle->Phi();
            clusterHijingMCStat.thetaMC  =particle->Theta();
            clusterHijingMCStat.vXmc     =particle->Vx();
            clusterHijingMCStat.vYmc     =particle->Vy();
            clusterHijingMCStat.vZmc     =particle->Vz();

            new((*fClusterHijingMCStatArray)[fClusterHijingMCStatCount]) ClusterHijingMCStatStruct( clusterHijingMCStat );
            fClusterHijingMCStatCount++;
        }

    }

    //____________________________________________________________
    //Reconstruction

    // === Physics Selection Task ===
    // bitwise operation is used against.
    // extraction de la masque de selection
    if(esd)
    {
        if( esd->IsPileupFromSPD(3,0.8) ) return;
    }


    //Remove events with exotic clusters
    TRefArray *caloClusArr=new TRefArray();
    event->GetEMCALClusters(caloClusArr);

    // const Int_t kNumber =caloClusArr->GetEntries();


    const UInt_t Mask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

    // event characterization
    if(esd)
    {
        fEventStatStruct->runNumber     = esd->GetRunNumber();
        TString triggerClasses = esd->GetFiredTriggerClasses();
    }

    if(aod)
    {
        fEventStatStruct->runNumber     = aod->GetRunNumber();

        // Get triggered classes, bordel....
        TString triggerClasses = aod->GetFiredTriggerClasses();
    }

    // ULong64_t triggerMask  = event->GetTriggerMask();
    ////verification triggered classes that fired.

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

    // Centrality
    if(esd)
    {
        esdCent     = esd->GetCentrality();
        //fleche pour pointers, otherwise point.
        fEventStatStruct->CentralityVZERO  = esdCent->GetCentralityPercentile("V0M");
        fEventStatStruct->CentralitySPD    = esdCent->GetCentralityPercentile("CL1");
        // number of particles
        // const Int_t nTracks( stack->GetNtrack() );
        fEventStatStruct->multiplicity =esd->GetNumberOfTracks();
        AliInfo( Form( "nTracks: %i",  fEventStatStruct->multiplicity) );

    }

    if(aod)
    {
        aodCent     = aod->GetCentrality();
        //fleche pour pointers, otherwise point.
        fEventStatStruct->CentralityVZERO  = aodCent->GetCentralityPercentile("V0M");
        fEventStatStruct->CentralitySPD    = aodCent->GetCentralityPercentile("CL1");
    }

    //Pass the geometry transformation matrix from ESDs to geometry
    AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
    AliVCaloCells     *Cells       =  event->GetEMCALCells();
    AliEMCALRecoUtils RecoUtils;

    //get reconstructed vertex position
    Double_t vertex_position[3];
    if(esd)esd->GetVertex()->GetXYZ(vertex_position);
    fEventStatStruct->  vX     = vertex_position[0];
    fEventStatStruct->  vY     = vertex_position[1];
    fEventStatStruct->  vZ     = vertex_position[2];


    if(aod)
    {
        fEventStatStruct->   vX     =0.0;
        fEventStatStruct->   vY     =0.0;
        fEventStatStruct->   vZ     =0.0;
    }

    if( fEventStatStruct->vZ>-15.&& fEventStatStruct->vZ<15. )
    {
        //variable size array
        std::set<int> particleLabels;

        //vertex cut
        //array temporaire pour passe plus tard dans le boucles
        TRefArray caloClustersArr = TRefArray();
        event->GetEMCALClusters( &caloClustersArr );

        const Int_t kNumberOfEMCALClusters =caloClustersArr.GetEntries();
        AliInfo( Form( "Number of clusters: %i", kNumberOfEMCALClusters) );


        TVector3 pos;
        pos -= vertex_position;

        //Double_t r1 = pos.Mag();

        Int_t nGoodClusters = 0;
        Int_t hijingclusters=0;
        Int_t piclusters=0;
        Int_t etaclusters=0;

        //boucle sur tous les clusters
        for( Int_t iclus = 0; iclus < kNumberOfEMCALClusters; iclus++ )
        {

            //first cluster
            // cast and check
            AliVCluster*c1=(AliVCluster*) caloClustersArr.At(iclus);
            if (!c1) continue;

            if (!c1->IsEMCAL()) continue;
            //================================
            if(c1->GetNCells()<2)continue;
            if(c1->GetDistanceToBadChannel()<2)continue;
            if(c1->E()<1)continue;
            if(c1->GetM02()>0.8)continue;
            if(c1->GetNTracksMatched()>0) continue;

            //create empty clusterStat object
            ClusterStatStruct clusterStat;

            /////////////////////////////////
            // AliVCaloCells *cells2 = 0;
            //cells2 = esd->GetEMCALCells();
            //if (!cells2) continue;
            //Int_t ncells2 = c1->GetNCells();
            ////////////////////////////////


            Int_t bc = 0;
            if(esd) bc = esd->GetBunchCrossNumber();
            if(aod) bc = aod->GetBunchCrossNumber();

            Int_t id= -1;;
            //id = TMath::Abs(Cells->GetCellNumber(icell));
            Float_t Emax    = GetMaxCellEnergy( c1, id);
            //AliEMCALGeometry  *geom        =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
            Float_t tcell=Cells->GetCellTime(id);
            Float_t Ecross =  RecoUtils.GetECross(id, tcell,Cells, bc);
            Float_t Exo    = 1.0 - Ecross/Emax;
            AliInfo( Form( "TIME CELL %f",Exo) );

            clusterStat.EMax    = Emax;
            clusterStat.Tcell  = tcell;
            clusterStat.ECross = Ecross;
            clusterStat.ExoC   = Exo;

            Int_t    recoPDG;
            Int_t    recoMotherPDG;

            // get contributors
            const Int_t nLabels( c1->GetNLabels());
            // look over labels
            Bool_t validHijing( kFALSE );
            Bool_t validpi( kFALSE );
            Bool_t valideta( kFALSE );

            //highest contributor
            const Int_t mainLabelIndex( c1->GetLabel() );

            if( mainLabelIndex < 0 ) continue;

            // get matching MC particle
            AliMCParticle *mainMCParticle( static_cast<AliMCParticle*>( mcEvent->GetTrack( mainLabelIndex ) ) );
            const Int_t iPart( mainMCParticle->Label() );

            // get matching particle in stack
            TParticle* mainParticle( stack->Particle(iPart) );

            // get energy ratio and store
            clusterStat.fraction= mainParticle->Energy()/c1->E();
            AliInfo( Form( "Cluster Energy %f", clusterStat.fraction ) );

            // reset number of pion, eta and hijing contributors to this cluster
            Int_t NpiContributeurs  =0;
            Int_t NetaContributeurs =0;
            Int_t NHiContributeurs  =0;




            // loop over all contributors
            // const Int_t nLabels( c1->GetNLabels());
            for( Int_t iLabel = 0; iLabel < nLabels; ++iLabel )
            {

                const Int_t labelIndex( c1->GetLabelAt( iLabel ) );
                if( labelIndex < 0 ) continue;

                // get matching MC particle
                AliMCParticle *mcParticle( static_cast<AliMCParticle*>( mcEvent->GetTrack( labelIndex ) ) );
                const Int_t iPart( mcParticle->Label() );



                //=================================
                //Hijing
                if( iPart >= firstHijingParticle && iPart < lastHijingParticle )
                {

                    //  AliInfo( Form( "First cluster:  %i has contributor from hijing generator", iclus ) );

                    // get corresponding particle:
                    TParticle* particle( stack->Particle(iPart) );
                    if(particle->GetMother(0)==-1)
                    {
                        clusterStat.RecoMomPDG=-99999;
                        clusterStat.RecoPDG   =particle->GetPdgCode();

                    }

                    else  {

                        TParticle* motherreco( stack->Particle(particle->GetMother(0)));
                        clusterStat.RecoMomPDG =motherreco->GetPdgCode();
                        clusterStat.RecoPDG    =particle->GetPdgCode();
                        /*  AliInfo( Form( "iPart 1: %i,  particle id: %i, parent id: %i,  parent pdg: %i",
                        iPart,
                        recoMotherPDG,
                        particle->GetMother(0),
                        recoMotherPDG
                        ));
                        */


                    } // parent test

                    // if current label is highest contributor, set flag
                    if( labelIndex == mainLabelIndex )validHijing = kTRUE;

                    // increment the number of hijing contributors to this cluster
                    NHiContributeurs++;
                    // if particle has not been already counted, add it to the set, and increment number of hijing clusters
                    if( particleLabels.find( labelIndex ) == particleLabels.end() )hijingclusters++;

                }

                //=================================
                //pi0
                if( iPart >= firstpiParticle && iPart < lastpiParticle )
                {
                    AliInfo( Form( "First cluster:  %i has contributor from pi0 generator", iclus ) );

                    // get corresponding particle:
                    TParticle* particle( stack->Particle(iPart) );
                    if(particle->GetMother(0)==-1)
                    {
                        clusterStat.RecoMomPDG   =-99999;
                        clusterStat.RecoPDG   =particle->GetPdgCode();

                    }

                    if(particle->GetMother(0)>-1)
                    {

                        TParticle* motherreco( stack->Particle(particle->GetMother(0)));
                        recoPDG       =particle->GetPdgCode();
                        recoMotherPDG =motherreco->GetPdgCode();
                        AliInfo( Form( "iPart 1: %i,  particle id: %i, parent id: %i,  parent pdg: %i",
                            iPart,
                            recoMotherPDG,
                            particle->GetMother(0),
                            recoMotherPDG
                            ));

                        clusterStat.RecoMomPDG  = recoMotherPDG;
                        clusterStat.RecoPDG     =  recoPDG;

                        //highest contributeur
                        if( labelIndex == mainLabelIndex )validpi = kTRUE;
                        NpiContributeurs++;
                        // if particle has not been already counted, add it to the set, and increment number of hijing clusters
                        if( particleLabels.find( labelIndex ) == particleLabels.end() ) piclusters++;

                    }

                }

                //=================================
                //eta
                if( iPart >= firstetaParticle && iPart < lastetaParticle )
                {
                    AliInfo( Form( "First cluster:  %i has contributor from eta generator", iclus ) );

                    // get corresponding particle:
                    TParticle* particle( stack->Particle(iPart) );
                    if(particle->GetMother(0)==-1)
                    {
                        clusterStat.RecoMomPDG   =-99999;
                        clusterStat.RecoPDG   =particle->GetPdgCode();

                    }
                    if(particle->GetMother(0)>-1)
                    {

                        TParticle* motherreco( stack->Particle(particle->GetMother(0)));
                        recoPDG       =particle->GetPdgCode();
                        recoMotherPDG =motherreco->GetPdgCode();
                        AliInfo( Form( "iPart 1: %i,  particle id: %i, parent id: %i,  parent pdg: %i",
                            iPart,
                            recoMotherPDG,
                            particle->GetMother(0),
                            recoMotherPDG
                            ));

                        clusterStat.RecoMomPDG  =recoMotherPDG;
                        clusterStat.RecoPDG     =recoPDG;

                        if( labelIndex == mainLabelIndex ) valideta = kTRUE;
                        NetaContributeurs++;
                        if( particleLabels.find( labelIndex ) == particleLabels.end() )  etaclusters++;



                    } // parent

                } // ipart test

                // add label to the particle set, to avoid double-counting
                particleLabels.insert( labelIndex );



            } // labels loop

            /////////////////////////////////////////////////////////////////////////////////////////////////

            clusterStat.validFlagHijing   = validHijing;
            clusterStat.validFlagpi       = validpi;
            clusterStat.validFlageta      = valideta;

            NpiContributeurs  =piclusters;
            NetaContributeurs =etaclusters;
            NHiContributeurs  =hijingclusters;

            clusterStat.Nlabels = c1->GetNLabels();
            clusterStat.Npi    = NpiContributeurs;
            clusterStat.Neta   = NetaContributeurs;
            clusterStat.Nhijing = NHiContributeurs;
            //================================

            TLorentzVector pii;
            c1->GetMomentum(pii, vertex_position);
            TLorentzVector *pimx = new TLorentzVector;
            c1->GetMomentum(*pimx, vertex_position);
            if(pii.Eta()>etaCut ||pii.Eta()<-etaCut )continue;
            if(pii.Phi() <phiCutMin  || pii.Phi()> phiCutMax)continue;


            ++nGoodClusters;

            //characteristiques cluster
            clusterStat.NCells                 = c1->GetNCells();
            clusterStat.distanceBadChannel     = c1->GetDistanceToBadChannel();
            clusterStat.Ecluster               = c1->E();
            clusterStat.M20cluster             = c1->GetM20();
            clusterStat.M02cluster             = c1->GetM02();
            clusterStat.isEMCALcluster         = c1->IsEMCAL();
            clusterStat.dispersioncluster      = c1->GetDispersion();
            clusterStat.chi2cluster            = c1->Chi2();
            clusterStat.distBadChannelcluster  = c1->GetDistanceToBadChannel();
            clusterStat.phicluster             = pii.Phi();
            clusterStat.etacluster             = pii.Eta();
            clusterStat.ptcluster              = pii.Pt();
            /////////// RecoUtils.GetMaxEnergyCell(geom, Cells, c1, absID1, iSM1, ieta1, iphi1, shared1);

            Int_t trackindex=c1->GetNTracksMatched();
            if(trackindex>0)clusterStat.TrackFlag =kTRUE;
            AliInfo( Form( "-----------matched Track: %i",  trackindex) );

            // copy clusterStat into TClonesArray, at last position
            new((*fClusterStatArray)[fClusterStatCount]) ClusterStatStruct( clusterStat );
            fClusterStatCount++;
        } //Premier cluster

        //
        fEventStatStruct->EMCalHijingClusters  =hijingclusters;
        fEventStatStruct->EMCalClusters=nGoodClusters;
        fEventStatStruct->NTotalClusters= kNumberOfEMCALClusters;

    }

    // save AOD
    AliAODHandler* handler = dynamic_cast<AliAODHandler*>( AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler() );
    if( handler )
    {

        // AliAODEvent* aod = handler->GetAOD();
        //AliAODHeader* header = aod->GetHeader();
        //header->SetRunNumber(aod->GetRunNumber());
        AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
    } else AliInfo( "Error: invalid output event handler" );

}//end process



//________________________________________________________________________
Double_t AliEMCalpi0AddedSignalsTask ::GetMaxCellEnergy(const AliVCluster *cluster, Int_t &id) const
{
    // Get maximum energy of attached cell.
    AliESDEvent* esd  = dynamic_cast<AliESDEvent*>(InputEvent());
    AliAODEvent* aod =dynamic_cast< AliAODEvent*>(InputEvent());

    // AliEMCALGeometry  *fGeom        =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");

    id = -1;

    AliVCaloCells *cells = 0;
    if(esd)
    {
        cells = esd->GetEMCALCells();
    }

    if(aod)
    {
        cells = aod->GetEMCALCells();
    }

    if (!cells)
        return 0;

    Double_t maxe = 0;
    Int_t ncells = cluster->GetNCells();
    for (Int_t i=0; i<ncells; i++)
    {
        Double_t e = cells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
        if (e>maxe)
        {
            maxe = e;
            id   = cluster->GetCellAbsId(i);
        }
    }
    return maxe;
}

//________________________________________________________________________
void AliEMCalpi0AddedSignalsTask::Terminate(const Option_t*)
{}

