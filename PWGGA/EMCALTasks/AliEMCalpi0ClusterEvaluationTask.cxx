/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*             /home/amorreal/Alice/Work/taskEMCal/AliEMCalpi0ClusterEvaluationTask.cxx
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
//
//-----------------------------------------------------------------------------
/// \class AliEMCalpi0ClusterEvaluationTask
/// AliAnalysisTaskSE to analyse data from ESDs (pi0 and Eta in the EMCal).
/// The Task reads as input ESDs then it produces a .root file:
/// general Histos and Ntuple:Tree
///
/// \author Astrid Morreale
//-----------------------------------------------------------------------------
#include "AliEMCalpi0ClusterEvaluationTask.h"

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
#include "AliEMCALCalibTimeDepCorrection.h" // Run dependent
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
ClassImp(AliEMCalpi0ClusterEvaluationTask)
///\endcond

//________________________________________________________________________
AliEMCalpi0ClusterEvaluationTask::AliEMCalpi0ClusterEvaluationTask(const char *name) :
AliAnalysisTaskSE(name),

fEvent(0),
ega0(0),      ega1(0),      ega2(0),      ega3(0),       ega4(0),
ega5(0),      ega6(0),      ega7(0),      ega8(0),       ega9(0),
mb0(0),       mb1(0),       mb2(0),       mb3(0),        mb4(0),
mb5(0),       mb6(0),       mb7(0),       mb8(0),        mb9(0),
allmb0(0),    allmb1(0),    allmb2(0),    allmb3(0),     allmb4(0),
allmb5(0),    allmb6(0),    allmb7(0),    allmb8(0),     allmb9(0),
cent0(0),     cent1(0),     cent2(0),     cent3(0),      cent4(0),
cent5(0),     cent6(0),     cent7(0),     cent8(0),      cent9(0),
semicent0(0), semicent1(0), semicent2(0), semicent3(0),  semicent4(0),
semicent5(0), semicent6(0), semicent7(0), semicent8(0),  semicent9(0),
all(0),       allmb(0),      mb(0),       central(0),   semicentral(0),
ega(0),

kAllMB(0),  isPileup(0),   isMB(0),        isAnyINT(0),       isCentral(0),  isSemiCentral(0), isEga(0),
isMBmx(0),  isAnyINTmx(0), isCentralmx(0), isSemiCentralmx(0), isEgamx(0),  kAllMBmx(0),
trigger(0), CentralityVZERO(0), CentralitySPD(0), runNumber(0), selectionMask(0),
vX(0),     vY(0),  vZ(0),

Ecluster(0),NCellscluster(0),M20cluster(0), M02cluster(0),NCluscluster(0),isEMCALcluster(0),
dispersioncluster(0),chi2cluster(0),distBadChannelcluster(0),phicluster(0),etacluster(0),
ptcluster(0), crossEnergy(0),
piE(0), piphi(0), pieta(0), ptpi(0), pipx(0), pipy(0), pipz(0), asympi(0), masspi(0),
fHistList(0)
{

    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    InitHistPointers();
}
//________________________________________________________________________
void AliEMCalpi0ClusterEvaluationTask::InitHistPointers() {

    for (Int_t i=0; i < kNtype; i++)
    {
        fMasspi0EGA[i]        = 0;
        fMassMixedEGA[i]      = 0;
        fEventsEGA[i]         = 0;

        fMasspi0MB[i]         = 0;
        fMassMixedMB[i]       = 0;
        fEventsMB[i]          = 0;

        fMasspi0AllMB[i]      = 0;
        fMassMixedAllMB[i]    = 0;
        fEventsAllMB[i]       = 0;


        fMasspi0Cent[i]      = 0;
        fMassMixedCent[i]    = 0;
        fEventsCent[i]       = 0;

        fMasspi0SemiCent[i]  = 0;
        fMassMixedSemiCent[i]= 0;
        fEventsSemiCent[i]   = 0;


        fpTMB[i]             = 0;
        fpTEGA[i]            = 0;
        fpTAllMB[i]          = 0;
        fpTkCent[i]          = 0;
        fpTkSemiCent[i]      = 0;


        fPool[i]             = 0;
    }

    fCentrality              = 0;
    fCentralityMB            = 0;
    fCentralityEGA           = 0;
    fCentralityCent          = 0;
    fCentralitySemiCent      = 0;
    fDispersion              = 0;
    fexo                     = 0;
    fTriggers                = 0;
    fshower                  = 0;

}

//________________________________________________________________________
AliEMCalpi0ClusterEvaluationTask::~AliEMCalpi0ClusterEvaluationTask()
{
    /// destructor
    //only if they are initialized to 0, or to something valid

    // Deleting the list (fHistList) also deletes all histograms stored in it
    // because of fHistList->SetOwner(). It is therefore not necessary to delete
    // the historgrams individually
    delete fHistList;

    // delete the 'running' pool array, which is either zero,
    // or corresponds with the last created pool during event processing.
    for (Int_t i=0; i < kNtype; i++)
    { delete fPool[i]; }

}

//________________________________________________________________________
void AliEMCalpi0ClusterEvaluationTask::UserCreateOutputObjects( void )
{
    fHistList = new TList();
    fHistList->SetOwner();

    for (Int_t i=0; i < kNtype; i++)
    {
        fMasspi0EGA[i]   = new TH2F(Form("MassEGA_%i",i),     Form("L1g m_{#gamma#gamma} centrality bin %i ;p_{T} (GeV/c)",i), 200,0,1, 150, 0, 30);
        fHistList->Add( fMasspi0EGA[i]  );
        fMassMixedEGA[i] = new TH2F(Form("MixedMassEGA_%i",i),Form("L1g mixed event m_{#gamma#gamma} centrality bin %i; p_{T} (GeV/c)",i), 200,0, 1, 150, 0, 30);
        fHistList->Add( fMassMixedEGA[i] );

        fEventsEGA[i]    = new TH1F(Form("EventsEGA_%i",i),   Form("L1g events in centrality bin %i",i),1, 0.5, 1.5);
        fHistList->Add( fEventsEGA[i]  );

        fMasspi0MB[i]    = new TH2F(Form("MassMB_%i",i),      Form("MB m_{#gamma#gamma} centrality bin %i ;p_{T} (GeV/c)",i),200,0,1, 150, 0, 30);
        fHistList->Add( fMasspi0MB[i]  );
        fMassMixedMB[i]  = new TH2F(Form("MixedMassMB_%i",i), Form("MB mixed event m_{#gamma#gamma} centrality bin %i;p_{T} (GeV/c)",i), 200,0, 1, 150, 0, 30);
        fHistList->Add(fMassMixedMB[i]);

        fEventsMB[i]     = new TH1F(Form("EventsMB_%i",i),    Form("MB events in centrality bin %i",i),1, 0.5, 1.5);
        fHistList->Add( fEventsMB[i]  );

        fMasspi0AllMB[i]  = new TH2F(Form("MassAllMB_%i",i),    Form("C+SC+AnyInt m_{#gamma#gamma} centrality bin %i ;p_{T} (GeV/c)",i),200,0,1, 150, 0, 30);
        fHistList->Add( fMasspi0AllMB[i]  );

        fMassMixedAllMB[i]= new TH2F(Form("MixedMassAllMB_%i",i),Form("C+SC+AnyInt mixed event m_{#gamma#gamma} centrality bin %i;p_{T} (GeV/c)",i), 200,0, 1, 150, 0, 30);
        fHistList->Add(fMassMixedAllMB[i]);

        fEventsAllMB[i]   = new TH1F(Form("EventsAllMB_%i",i),  Form("C+SC+ AnyInt Events in centrality bin %i",i),1, 0.5, 1.5);
        fHistList->Add( fEventsAllMB[i]  );

        fpTMB[i]   = new TH1F(Form("pTMB_%i",i),  Form("pT in centrality bin %i",i),150, 1, 30);
        fHistList->Add( fpTMB[i]  );

        fpTAllMB[i]   = new TH1F(Form("pTAllMB_%i",i),  Form("pT in centrality bin %i",i),150, 1, 30);
        fHistList->Add( fpTAllMB[i]  );

        fpTEGA[i]  = new TH1F(Form("pTEGA_%i",i),  Form("pT in centrality bin %i",i),150, 1, 30);
        fHistList->Add( fpTEGA[i]  );

        fpTkCent[i]  = new TH1F(Form("pTkCent_%i",i),  Form("pT in centrality bin %i",i),150, 1, 30);
        fHistList->Add( fpTkCent[i]  );

        fpTkSemiCent[i]  = new TH1F(Form("pTkSemiCent_%i",i),  Form("pT in centrality bin %i",i),150, 1, 30);
        fHistList->Add( fpTkSemiCent[i]  );

        //refined kCentral and SemiCentral

        fMasspi0SemiCent[i]  = new TH2F(Form("MassSemiCent_%i",i),    Form("SemiCent m_{#gamma#gamma} centrality bin %i ;p_{T} (GeV/c)",i),200,0,1, 150, 0, 30);
        fHistList->Add( fMasspi0SemiCent[i]  );

        fMassMixedSemiCent[i]= new TH2F(Form("MixedMassSemiCent_%i",i),Form("SemiCent mixed event m_{#gamma#gamma} centrality bin %i;p_{T} (GeV/c)",i), 200,0, 1, 150, 0, 30);
        fHistList->Add(fMassMixedSemiCent[i]);

        fEventsSemiCent[i]   = new TH1F(Form("EventsSemiCent_%i",i),  Form("SemiCent Events in centrality bin %i",i),1, 0.5, 1.5);
        fHistList->Add( fEventsSemiCent[i]  );

        fMasspi0Cent[i]  = new TH2F(Form("MassCent_%i",i),    Form("C+SC+AnyInt m_{#gamma#gamma} centrality bin %i ;p_{T} (GeV/c)",i),200,0,1, 150, 0, 30);
        fHistList->Add( fMasspi0Cent[i]  );

        fMassMixedCent[i]= new TH2F(Form("MixedMassCent_%i",i),Form("C+SC+AnyInt mixed event m_{#gamma#gamma} centrality bin %i;p_{T} (GeV/c)",i), 200,0, 1, 150, 0, 30);
        fHistList->Add(fMassMixedCent[i]);

        fEventsCent[i]   = new TH1F(Form("EventsCent_%i",i),  Form("C+SC+ AnyInt Events in centrality bin %i",i),1, 0.5, 1.5);
        fHistList->Add(fEventsCent[i]);


    }

    fCentrality   = new TH1F("centrality", "centrality",         100, 0, 100);
    fHistList->Add(fCentrality);

    fCentralityMB   = new TH1F("centralityMB", "centrality",     100, 0, 100);
    fHistList->Add(fCentralityMB);

    fCentralityEGA   = new TH1F("centralityEGA", "centrality",   100, 0, 100);
    fHistList->Add(fCentralityEGA);

    fCentralityCent   = new TH1F("centralityCent", "centrality", 100, 0, 100);
    fHistList->Add(fCentralityCent);

    fCentralitySemiCent   = new TH1F("centralitySemiCent", "centrality", 100, 0, 100);
    fHistList->Add(fCentralitySemiCent);

    fDispersion   = new TH1F("dispersion", "dispersion", 100, -1, 2);
    fHistList->Add(fDispersion);

    fexo   = new TH1F("exo", "exo", 100, -1, 2);
    fHistList->Add(fexo);

    fshower   = new TH1F("shower", "shower", 100, -1, 2);
    fHistList->Add(fshower);


    fTriggers     = new TH1F("triggers",   "triggers",   10,  0, 10);
    fTriggers->GetXaxis()->SetBinLabel( 1,"All");
    fTriggers->GetXaxis()->SetBinLabel( 2,"AllMB");
    fTriggers->GetXaxis()->SetBinLabel( 3,"MB");
    fTriggers->GetXaxis()->SetBinLabel( 4,"kCentral");
    fTriggers->GetXaxis()->SetBinLabel( 5,"kSemiCentral");
    fTriggers->GetXaxis()->SetBinLabel( 6,"EGA");
    fTriggers->GetXaxis()->SetBinLabel( 7,"MB/ALL");
    fTriggers->GetXaxis()->SetBinLabel( 8,"MB/EGA");
    fTriggers->GetXaxis()->SetBinLabel( 9,"MB/KCentral");
    fTriggers->GetXaxis()->SetBinLabel( 10,"MB/kSemiCentral");
    fHistList->Add(fTriggers);

    PostData(1, fHistList);

}

//________________________________________________________________________


//________________________________________________________________________
void AliEMCalpi0ClusterEvaluationTask::UserExec( Option_t* )
{

    // increase event number
    AliInfo( Form( "event: %i", fEvent++ ) );
    fTriggers->SetBinContent(7,   fEvent);
    // ESD Filter analysis task executed for each event

    AliESDEvent* esd    = dynamic_cast<AliESDEvent*>(InputEvent());
    AliAODEvent* aod    = dynamic_cast< AliAODEvent*>(InputEvent());
    AliVEvent  * event  = InputEvent();
    // === Physics Selection Task ===
    // bitwise operation is used against.
    // extraction de la masque de selection
    if(esd){
    isPileup       = esd->IsPileupFromSPD(3,0.8);
    if(isPileup) return;
           }

    //Remove events with exotic clusters
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

        Int_t id= -1;;
        Double_t Emax   = GetMaxCellEnergy( c, id);

        AliVCaloCells     *Cells       =  event->GetEMCALCells();
       // AliEMCALGeometry  *geom        =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
        AliEMCALRecoUtils RecoUtils;
        Int_t bc;
        if(esd)  bc = esd->GetBunchCrossNumber();
        if(aod) bc = aod->GetBunchCrossNumber();

        Double_t tcell=0;
        Double_t Ecross =  RecoUtils.GetECross(id, tcell,Cells, bc);
        Double_t Exo    = 1.0 - Ecross/Emax;
        fexo->Fill(Exo);
        if((Exo)>1){ return;}
        fshower->Fill(c->GetM02());

    }
    //888888888888888888888888888888888888888888888888888

    const UInt_t Mask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

    // event characterization
    if(esd){
        runNumber     = esd->GetRunNumber();
        TString triggerClasses = esd->GetFiredTriggerClasses();
    }

    if(aod){
        runNumber     = aod->GetRunNumber();
       // Get triggered classes, bordel....
    TString triggerClasses = aod->GetFiredTriggerClasses();
            }


   // ULong64_t triggerMask  = event->GetTriggerMask();

    ////verification triggered classes that fired.

    isMB           = Mask & AliVEvent::kMB;
    isAnyINT       = Mask & AliVEvent::kAnyINT;
    isCentral      = Mask & AliVEvent::kCentral;
    isSemiCentral  = Mask & AliVEvent::kSemiCentral;
    isEga          = Mask & AliVEvent::kEMCEGA;
    /*
    isMB             = (Mask== AliVEvent::kMB)? 1 : 0;
    isAnyINT         = (Mask== AliVEvent::kAnyINT)? 1 : 0;
    isCentral        = (Mask== AliVEvent::kCentral)? 1 : 0;
    isSemiCentral    = (Mask== AliVEvent::kSemiCentral)? 1 : 0;
    isEga            = (Mask== AliVEvent::kEMCEGA)? 1 : 0;
    */
    //if( isMB ||isCentral||isSemiCentral ){ kAllMB=kTRUE;}
    //else kAllMB = kFALSE;
    kAllMB = (isMB ||isCentral||isSemiCentral);

  AliCentrality *esdCent;
  AliCentrality *aodCent;
    // Centrality
if(esd){

    esdCent     = esd->GetCentrality();

    //fleche pour pointers, otherwise point.
    CentralityVZERO  = esdCent->GetCentralityPercentile("V0M");
    CentralitySPD    = esdCent->GetCentralityPercentile("CL1");

}

if(aod){
    aodCent     = aod->GetCentrality();
    //fleche pour pointers, otherwise point.
    CentralityVZERO  = aodCent->GetCentralityPercentile("V0M");
    CentralitySPD    = aodCent->GetCentralityPercentile("CL1");
}

    all++;
    if(kAllMB){allmb++;fCentrality->Fill(CentralityVZERO);}
    if(isMB){mb++;fCentralityMB->Fill(CentralityVZERO);}
    if(isEga){ega++;fCentralityEGA->Fill(CentralityVZERO);}
    if(isCentral){central++;fCentralityCent->Fill(CentralityVZERO);}
    if(isSemiCentral){semicentral++;fCentralitySemiCent->Fill(CentralityVZERO);}


    fTriggers->SetBinContent(1,  all);
    fTriggers->SetBinContent(2,  allmb);
    fTriggers->SetBinContent(3,  mb);
    fTriggers->SetBinContent(4,  central);
    fTriggers->SetBinContent(5,  semicentral);
    fTriggers->SetBinContent(6,  ega);
    //     fTriggers->SetBinContent(7,  mb/ega);
    //     fTriggers->SetBinContent(8,  mb/central);
    //     fTriggers->SetBinContent(9,  mb/semicentral);
    //     fTriggers->SetBinContent(10, ega/mb);




    // Fill per centrality class event count
    if(CentralityVZERO>0  && CentralityVZERO<=5 )
    {
        if(isEga){ega0++;fEventsEGA[0]->Fill(ega0);}
        if(isMB){mb0++;fEventsMB[0]->Fill(mb0);}
        if(kAllMB){allmb0++;fEventsAllMB[0]->Fill(allmb0);}
        if(isCentral){cent0++;fEventsCent[0]->Fill(cent0);}
        if(isSemiCentral){semicent0++;fEventsSemiCent[0]->Fill(semicent0);}
    }
    if(CentralityVZERO>5  && CentralityVZERO<=10){if(isEga){ega1++;fEventsEGA[1]->Fill(ega1);}if(isMB){mb1++;fEventsMB[1]->Fill(mb1);}if(kAllMB){allmb1++;fEventsAllMB[1]->Fill(allmb1);} if(isCentral){cent1++;fEventsCent[1]->Fill(cent1);} if(isSemiCentral){semicent1++;fEventsSemiCent[1]->Fill(semicent1);}}
    if(CentralityVZERO>10 && CentralityVZERO<=20){if(isEga){ega2++;fEventsEGA[2]->Fill(ega2);}if(isMB){mb2++;fEventsMB[2]->Fill(mb2);}if(kAllMB){allmb2++;fEventsAllMB[2]->Fill(allmb2);} if(isCentral){cent2++;fEventsCent[2]->Fill(cent2);} if(isSemiCentral){semicent2++;fEventsSemiCent[2]->Fill(semicent2);}}
    if(CentralityVZERO>20 && CentralityVZERO<=40){if(isEga){ega3++;fEventsEGA[3]->Fill(ega3);}if(isMB){mb3++;fEventsMB[3]->Fill(mb3);}if(kAllMB){allmb3++;fEventsAllMB[3]->Fill(allmb3);} if(isCentral){cent3++;fEventsCent[3]->Fill(cent3);} if(isSemiCentral){semicent3++;fEventsSemiCent[3]->Fill(semicent3);}}
    if(CentralityVZERO>40 && CentralityVZERO<=60){if(isEga){ega4++;fEventsEGA[4]->Fill(ega4);}if(isMB){mb4++;fEventsMB[4]->Fill(mb4);}if(kAllMB){allmb4++;fEventsAllMB[4]->Fill(allmb4);} if(isCentral){cent4++;fEventsCent[4]->Fill(cent4);} if(isSemiCentral){semicent4++;fEventsSemiCent[4]->Fill(semicent4);}}
    if(CentralityVZERO>60 && CentralityVZERO<=80){if(isEga){ega5++;fEventsEGA[5]->Fill(ega5);}if(isMB){mb5++;fEventsMB[5]->Fill(mb5);}if(kAllMB){allmb5++;fEventsAllMB[5]->Fill(allmb5);} if(isCentral){cent5++;fEventsCent[5]->Fill(cent5);} if(isSemiCentral){semicent5++;fEventsSemiCent[5]->Fill(semicent5);}}
    if(CentralityVZERO>0  && CentralityVZERO<=10){if(isEga){ega6++;fEventsEGA[6]->Fill(ega6);}if(isMB){mb6++;fEventsMB[6]->Fill(mb6);}if(kAllMB){allmb6++;fEventsAllMB[6]->Fill(allmb6);} if(isCentral){cent6++;fEventsCent[6]->Fill(cent6);} if(isSemiCentral){semicent6++;fEventsSemiCent[6]->Fill(semicent6);}}
    if(CentralityVZERO>0  && CentralityVZERO<=20){if(isEga){ega7++;fEventsEGA[7]->Fill(ega7);}if(isMB){mb7++;fEventsMB[7]->Fill(mb7);}if(kAllMB){allmb7++;fEventsAllMB[7]->Fill(allmb7);} if(isCentral){cent7++;fEventsCent[7]->Fill(cent7);} if(isSemiCentral){semicent7++;fEventsSemiCent[7]->Fill(semicent7);}}
    if(CentralityVZERO>40 && CentralityVZERO<=50){if(isEga){ega8++;fEventsEGA[8]->Fill(ega8);}if(isMB){mb8++;fEventsMB[8]->Fill(mb8);}if(kAllMB){allmb8++;fEventsAllMB[8]->Fill(allmb8);} if(isCentral){cent8++;fEventsCent[8]->Fill(cent8);} if(isSemiCentral){semicent8++;fEventsSemiCent[8]->Fill(semicent8);}}
    if(CentralityVZERO>50 && CentralityVZERO<=60){if(isEga){ega9++;fEventsEGA[9]->Fill(ega9);}if(isMB){mb9++;fEventsMB[9]->Fill(mb9);}if(kAllMB){allmb9++;fEventsAllMB[9]->Fill(allmb9);} if(isCentral){cent9++;fEventsCent[9]->Fill(cent9);} if(isSemiCentral){semicent9++;fEventsSemiCent[9]->Fill(semicent9);}}




    //Pass the geometry transformation matrix from ESDs to geometry
    AliEMCALGeometry  *geom        =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
    AliVCaloCells     *Cells       =  event->GetEMCALCells();
    AliEMCALRecoUtils RecoUtils;

    Int_t    absID1    = -1;        Int_t    absID2 = -1;
    Int_t    ieta1     = -1;        Int_t    ieta2  = -1;
    Int_t    iphi1     = -1;        Int_t    iphi2  = -1;
    Int_t    iSM1      = -1;        Int_t    iSM2   = -1;
    Bool_t  shared1;      Bool_t  shared2;

    //get reconstructed vertex position
    Double_t vertex_position[3];
    if(esd)esd->GetVertex()->GetXYZ(vertex_position);
    vX     = vertex_position[0];
    vY     = vertex_position[1];
    vZ     = vertex_position[2];


    if(aod){
    vX     =0.0;
    vY     =0.0;
    vZ     =0.0;
}


    // cout<<vZ<<endl;
    if(vZ>-15.||vZ<15.){//vertex cut

        //array temporaire pour passe plus tard dans le boucles
        TRefArray caloClustersArr = TRefArray();
        event->GetEMCALClusters( &caloClustersArr );

        const Int_t kNumberOfEMCALClusters =caloClustersArr.GetEntries();

        // pool to store clusters to be mixed
        // it is not deleted directly. Instead, it is assigned to fPool, the previously assigned pool
        // being deleted (if any), at that time.
        TObjArray *newPool = new TObjArray(kNumberOfEMCALClusters);
        newPool->SetOwner();

        TVector3 pos;
        pos -= vertex_position;
        //Double_t r1 = pos.Mag();

        //Fill the pool with all clusters
        Int_t nGoodClusters = 0;
        //boucle sur tous les clusters
        for( Int_t iclus = 0; iclus < kNumberOfEMCALClusters-1; iclus++ )
        {//first cluster

            AliVCluster*c1=(AliVCluster*) caloClustersArr.At(iclus);
            if (!c1) continue;
            if (!c1->IsEMCAL()) continue;
            if (c1->GetNCells()<2) continue;
           //if (c1->E()<2) continue;
           //if(c1->GetM02()>0.5) continue;
            if (c1->GetDistanceToBadChannel()<2) continue;

            TLorentzVector pii;
            c1->GetMomentum(pii, vertex_position);
            TLorentzVector *pimx = new TLorentzVector;
            c1->GetMomentum(*pimx, vertex_position);

            // add TLorentz vector in the pool array
            // it will be deleted when the array is deleted
            newPool->Add(pimx);
            ++nGoodClusters;

            //characteristiques cluster
            Ecluster               = c1->E();
            NCellscluster          = c1->GetNCells();
            M20cluster             = c1->GetM20();
            M02cluster             = c1->GetM02();
            NCluscluster           = kNumberOfEMCALClusters;
            isEMCALcluster         = c1->IsEMCAL();
            dispersioncluster      = c1->GetDispersion();
            chi2cluster            = c1->Chi2();
            distBadChannelcluster  = c1->GetDistanceToBadChannel();
            phicluster             = pii.Phi();
            etacluster             = pii.Eta();
            ptcluster              = pii.Pt();
            RecoUtils.GetMaxEnergyCell(geom, Cells, c1, absID1, iSM1, ieta1, iphi1, shared1);

            fDispersion->Fill(dispersioncluster);

            for (Int_t iclus2 = iclus+1;  iclus2< kNumberOfEMCALClusters-1; iclus2++)
            {//second cluster

                AliVCluster *c2 = (AliVCluster*) caloClustersArr.At(iclus2);
                if (!c2) continue;
                if (!c2->IsEMCAL()) continue;
                if (c2->GetNCells()<2) continue;
                //if (c2->E()<2) continue;//for evaluation of cuts
                //if(c2->GetM02()>0.5) continue;//for evaluation of cuts
                if (c2->GetDistanceToBadChannel()<2) continue;
               // Float_t en2 = c2->E();

                TLorentzVector pjj;
                c2->GetMomentum(pjj, vertex_position);
                TLorentzVector pion;
                pion   = pii + pjj;
                //remplissage des pions
                piE      = pion.E();
                piphi    = pion.Phi();
                pieta    = pion.Eta();
                ptpi     = pion.Pt();
                pipx     = pion.Px();
                pipy     = pion.Py();
                pipz     = pion.Pz();
                asympi     = TMath::Abs(pii.E()-pjj.E())/(pii.E()+pjj.E());
                masspi   = pion.M();

                if(CentralityVZERO>0  && CentralityVZERO<=5 ){if(isEga){fpTEGA[0]->Fill(ptpi);fMasspi0EGA[0]->Fill(masspi, ptpi);}if(isMB){fpTMB[0]->Fill(ptpi); fMasspi0MB[0]->Fill(masspi, ptpi);}if(kAllMB){fpTAllMB[0]->Fill(ptpi);fMasspi0AllMB[0]->Fill(masspi, ptpi);}if(isCentral){fpTkCent[0]->Fill(ptpi);fMasspi0Cent[0]->Fill(masspi, ptpi);}if(isSemiCentral){fpTkSemiCent[0]->Fill(ptpi);fMasspi0SemiCent[0]->Fill(masspi, ptpi);}}
                if(CentralityVZERO>5  && CentralityVZERO<=10){if(isEga){fpTEGA[1]->Fill(ptpi);fMasspi0EGA[1]->Fill(masspi, ptpi);}if(isMB){fpTMB[1]->Fill(ptpi); fMasspi0MB[1]->Fill(masspi, ptpi);}if(kAllMB){fpTAllMB[1]->Fill(ptpi);fMasspi0AllMB[1]->Fill(masspi, ptpi);}if(isCentral){fpTkCent[1]->Fill( ptpi);fMasspi0Cent[1]->Fill(masspi, ptpi);} if(isSemiCentral){fpTkSemiCent[1]->Fill( ptpi);fMasspi0SemiCent[1]->Fill(masspi, ptpi);}}
                if(CentralityVZERO>10 && CentralityVZERO<=20){if(isEga){fpTEGA[2]->Fill(ptpi);fMasspi0EGA[2]->Fill(masspi, ptpi);}if(isMB){fpTMB[2]->Fill(ptpi); fMasspi0MB[2]->Fill(masspi, ptpi);}if(kAllMB){fpTAllMB[2]->Fill(ptpi);fMasspi0AllMB[2]->Fill(masspi, ptpi);}if(isCentral){fpTkCent[2]->Fill( ptpi);fMasspi0Cent[2]->Fill(masspi, ptpi);} if(isSemiCentral){fpTkSemiCent[2]->Fill( ptpi);fMasspi0SemiCent[2]->Fill(masspi, ptpi);}}
                if(CentralityVZERO>20 && CentralityVZERO<=40){if(isEga){fpTEGA[3]->Fill(ptpi);fMasspi0EGA[3]->Fill(masspi, ptpi);}if(isMB){fpTMB[3]->Fill(ptpi); fMasspi0MB[3]->Fill(masspi, ptpi);}if(kAllMB){fpTAllMB[3]->Fill(ptpi);fMasspi0AllMB[3]->Fill(masspi, ptpi);}if(isCentral){fpTkCent[3]->Fill( ptpi);fMasspi0Cent[3]->Fill(masspi, ptpi);} if(isSemiCentral){fpTkSemiCent[3]->Fill( ptpi);fMasspi0SemiCent[3]->Fill(masspi, ptpi);}}
                if(CentralityVZERO>40 && CentralityVZERO<=60){if(isEga){fpTEGA[4]->Fill(ptpi);fMasspi0EGA[4]->Fill(masspi, ptpi);}if(isMB){fpTMB[4]->Fill(ptpi); fMasspi0MB[4]->Fill(masspi, ptpi);}if(kAllMB){fpTAllMB[4]->Fill(ptpi);fMasspi0AllMB[4]->Fill(masspi, ptpi);}if(isCentral){fpTkCent[4]->Fill( ptpi);fMasspi0Cent[4]->Fill(masspi, ptpi);} if(isSemiCentral){fpTkSemiCent[4]->Fill( ptpi);fMasspi0SemiCent[4]->Fill(masspi, ptpi);}}
                if(CentralityVZERO>60 && CentralityVZERO<=80){if(isEga){fpTEGA[5]->Fill(ptpi);fMasspi0EGA[5]->Fill(masspi, ptpi);}if(isMB){fpTMB[5]->Fill(ptpi); fMasspi0MB[5]->Fill(masspi, ptpi);}if(kAllMB){fpTAllMB[5]->Fill(ptpi);fMasspi0AllMB[5]->Fill(masspi, ptpi);}if(isCentral){fpTkCent[5]->Fill( ptpi);fMasspi0Cent[5]->Fill(masspi, ptpi);} if(isSemiCentral){fpTkSemiCent[5]->Fill( ptpi);fMasspi0SemiCent[5]->Fill(masspi, ptpi);}}
                if(CentralityVZERO>0  && CentralityVZERO<=10){if(isEga){fpTEGA[6]->Fill(ptpi);fMasspi0EGA[6]->Fill(masspi, ptpi);}if(isMB){fpTMB[6]->Fill(ptpi); fMasspi0MB[6]->Fill(masspi, ptpi);}if(kAllMB){fpTAllMB[6]->Fill(ptpi);fMasspi0AllMB[6]->Fill(masspi, ptpi);}if(isCentral){fpTkCent[6]->Fill( ptpi);fMasspi0Cent[6]->Fill(masspi, ptpi);} if(isSemiCentral){fpTkSemiCent[6]->Fill( ptpi);fMasspi0SemiCent[6]->Fill(masspi, ptpi);}}
                if(CentralityVZERO>0  && CentralityVZERO<=20){if(isEga){fpTEGA[7]->Fill(ptpi);fMasspi0EGA[7]->Fill(masspi, ptpi);}if(isMB){fpTMB[7]->Fill(ptpi); fMasspi0MB[7]->Fill(masspi, ptpi);}if(kAllMB){fpTAllMB[7]->Fill(ptpi);fMasspi0AllMB[7]->Fill(masspi, ptpi);}if(isCentral){fpTkCent[7]->Fill( ptpi);fMasspi0Cent[7]->Fill(masspi, ptpi);} if(isSemiCentral){fpTkSemiCent[7]->Fill( ptpi);fMasspi0SemiCent[7]->Fill(masspi, ptpi);}}
                if(CentralityVZERO>40 && CentralityVZERO<=50){if(isEga){fpTEGA[8]->Fill(ptpi);fMasspi0EGA[8]->Fill(masspi, ptpi);}if(isMB){fpTMB[8]->Fill(ptpi); fMasspi0MB[8]->Fill(masspi, ptpi);}if(kAllMB){fpTAllMB[8]->Fill(ptpi);fMasspi0AllMB[8]->Fill(masspi, ptpi);}if(isCentral){fpTkCent[8]->Fill( ptpi);fMasspi0Cent[8]->Fill(masspi, ptpi);} if(isSemiCentral){fpTkSemiCent[8]->Fill( ptpi);fMasspi0SemiCent[8]->Fill(masspi, ptpi);}}
                if(CentralityVZERO>50 && CentralityVZERO<=60){if(isEga){fpTEGA[9]->Fill(ptpi);fMasspi0EGA[9]->Fill(masspi, ptpi);}if(isMB){fpTMB[9]->Fill(ptpi); fMasspi0MB[9]->Fill(masspi, ptpi);}if(kAllMB){fpTAllMB[9]->Fill(ptpi);fMasspi0AllMB[9]->Fill(masspi, ptpi);}if(isCentral){fpTkCent[9]->Fill( ptpi);fMasspi0Cent[9]->Fill(masspi, ptpi);} if(isSemiCentral){fpTkSemiCent[9]->Fill( ptpi);fMasspi0SemiCent[9]->Fill(masspi, ptpi);}}

                RecoUtils.GetMaxEnergyCell(geom, Cells, c2, absID2, iSM2,  ieta2,  iphi2, shared2);

            } //Deuxieme cluster

        } //Premier cluster

        //     //-----------------------------------------
        //fill mixed event
        int poolIndex = 0;
        if(CentralityVZERO>0  && CentralityVZERO<=5 ) poolIndex =0;
        if(CentralityVZERO>5  && CentralityVZERO<=10) poolIndex =1;
        if(CentralityVZERO>10 && CentralityVZERO<=20) poolIndex =2;
        if(CentralityVZERO>20 && CentralityVZERO<=40) poolIndex =3;
        if(CentralityVZERO>40 && CentralityVZERO<=60) poolIndex =4;
        if(CentralityVZERO>60 && CentralityVZERO<=80) poolIndex =5;
        if(CentralityVZERO>0  && CentralityVZERO<=10) poolIndex =6;
        if(CentralityVZERO>0  && CentralityVZERO<=20) poolIndex =7;
        if(CentralityVZERO>40 && CentralityVZERO<=50) poolIndex =8;
        if(CentralityVZERO>50 && CentralityVZERO<=60) poolIndex =9;


        if( !fPool[poolIndex] )
        {
            fPool[poolIndex] = newPool;

        } else {

            Int_t nGoodClusters2 = fPool[poolIndex]->GetEntries();
            for (Int_t i=0; i<nGoodClusters; ++i)
            {


                TLorentzVector * pi = static_cast<TLorentzVector*>(newPool->At(i));
                for (Int_t j=0; j<nGoodClusters2; ++j)
                {

                    TLorentzVector *pj= static_cast<TLorentzVector*>(fPool[poolIndex]->At(j));
                    FillMixed(*pi,*pj);


                } //random cluster
            } //current cluster

            // delete previous pool and assign to new pool.
            delete fPool[poolIndex];
            fPool[poolIndex] = newPool;

        }
    }//15 cm vertex cut
    PostData(1, fHistList);

}//end process

//________________________________________________________________________
Double_t AliEMCalpi0ClusterEvaluationTask ::GetMaxCellEnergy(const AliVCluster *cluster, Int_t &id) const
{
    // Get maximum energy of attached cell.
    AliESDEvent* esd  = dynamic_cast<AliESDEvent*>(InputEvent());
    AliAODEvent* aod =dynamic_cast< AliAODEvent*>(InputEvent());




   // AliEMCALGeometry  *fGeom        =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");

    id = -1;

    AliVCaloCells *cells = 0;
    if(esd){
    cells = esd->GetEMCALCells();
     }

    if(!esd){
      if (aod)cells = aod->GetEMCALCells();
	   else cells = 0;
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
void AliEMCalpi0ClusterEvaluationTask::FillMixed( const TLorentzVector& p1, const TLorentzVector& p2)

{
    //verification triggered classes that fired.
    const UInt_t eventSelectionMask( ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() );
    AliESDEvent* esd     = dynamic_cast<AliESDEvent*>(InputEvent());
    AliAODEvent* aod     =dynamic_cast< AliAODEvent*>(InputEvent());
    //AliVEvent  * event   = InputEvent();


   // ULong64_t triggerMask = event->GetTriggerMask();

    isMBmx         = (eventSelectionMask&AliVEvent::kMB);
    isAnyINTmx      = (eventSelectionMask&AliVEvent::kAnyINT);
    isCentralmx     = (eventSelectionMask&AliVEvent::kCentral);
    isSemiCentralmx = (eventSelectionMask&AliVEvent::kSemiCentral);
    isEgamx         = (eventSelectionMask&AliVEvent::kEMCEGA);



    //     isMBmx          = (eventSelectionMask== AliVEvent::kMB)? 1 : 0;
    //     isAnyINTmx      = (eventSelectionMask== AliVEvent::kAnyINT)? 1 : 0;
    //     isCentralmx     = (eventSelectionMask== AliVEvent::kCentral)? 1 : 0;
    //     isSemiCentralmx = (eventSelectionMask== AliVEvent::kSemiCentral)? 1 : 0;
    //     isEgamx         = (eventSelectionMask== AliVEvent::kEMCEGA)? 1 : 0;
    kAllMBmx= (isMBmx || isAnyINTmx||isCentralmx||isSemiCentralmx );

    // Centrality
     AliCentrality *esdCent;
    if(esd){
     esdCent     = esd->GetCentrality();
     }
    if(!esd){
     if (aod){
		esdCent     = aod->GetCentrality();
	 } else {
		esdCent = NULL;
		return; 
	 } 	 
     }

    //fleche pour pointers, otherwise point.
    CentralityVZERO  = esdCent->GetCentralityPercentile("V0M");
    CentralitySPD    = esdCent->GetCentralityPercentile("CL1");

    TLorentzVector Mxpion;
    Mxpion = p1 + p2;
    Double_t Mxmass      = Mxpion.M();
    Double_t Mxpt        = Mxpion.Pt();

    if(Mxpt>=4){
        if(CentralityVZERO>0  && CentralityVZERO<=5 ){if(isEgamx){fMassMixedEGA[0]->Fill(Mxmass, Mxpt);}if(isMBmx){fMassMixedMB[0]->Fill(Mxmass, Mxpt);}if(kAllMB){fMassMixedAllMB[0]->Fill(Mxmass, Mxpt);} if(isCentralmx){fMassMixedCent[0]->Fill(Mxmass, Mxpt);}if(isSemiCentralmx){fMassMixedSemiCent[0]->Fill(Mxmass, Mxpt);}}
        if(CentralityVZERO>5  && CentralityVZERO<=10){if(isEgamx){fMassMixedEGA[1]->Fill(Mxmass, Mxpt);}if(isMBmx){fMassMixedMB[1]->Fill(Mxmass, Mxpt);}if(kAllMB){fMassMixedAllMB[1]->Fill(Mxmass, Mxpt);} if(isCentralmx){fMassMixedCent[1]->Fill(Mxmass, Mxpt);}if(isSemiCentralmx){fMassMixedSemiCent[1]->Fill(Mxmass, Mxpt);}}
        if(CentralityVZERO>10 && CentralityVZERO<=20){if(isEgamx){fMassMixedEGA[2]->Fill(Mxmass, Mxpt);}if(isMBmx){fMassMixedMB[2]->Fill(Mxmass, Mxpt);}if(kAllMB){fMassMixedAllMB[2]->Fill(Mxmass, Mxpt);} if(isCentralmx){fMassMixedCent[2]->Fill(Mxmass, Mxpt);}if(isSemiCentralmx){fMassMixedSemiCent[2]->Fill(Mxmass, Mxpt);}}
        if(CentralityVZERO>20 && CentralityVZERO<=40){if(isEgamx){fMassMixedEGA[3]->Fill(Mxmass, Mxpt);}if(isMBmx){fMassMixedMB[3]->Fill(Mxmass, Mxpt);}if(kAllMB){fMassMixedAllMB[3]->Fill(Mxmass, Mxpt);} if(isCentralmx){fMassMixedCent[3]->Fill(Mxmass, Mxpt);}if(isSemiCentralmx){fMassMixedSemiCent[3]->Fill(Mxmass, Mxpt);}}
        if(CentralityVZERO>40 && CentralityVZERO<=60){if(isEgamx){fMassMixedEGA[4]->Fill(Mxmass, Mxpt);}if(isMBmx){fMassMixedMB[4]->Fill(Mxmass, Mxpt);}if(kAllMB){fMassMixedAllMB[4]->Fill(Mxmass, Mxpt);} if(isCentralmx){fMassMixedCent[4]->Fill(Mxmass, Mxpt);}if(isSemiCentralmx){fMassMixedSemiCent[4]->Fill(Mxmass, Mxpt);}}
        if(CentralityVZERO>60 && CentralityVZERO<=80){if(isEgamx){fMassMixedEGA[5]->Fill(Mxmass, Mxpt);}if(isMBmx){fMassMixedMB[5]->Fill(Mxmass, Mxpt);}if(kAllMB){fMassMixedAllMB[5]->Fill(Mxmass, Mxpt);} if(isCentralmx){fMassMixedCent[5]->Fill(Mxmass, Mxpt);}if(isSemiCentralmx){fMassMixedSemiCent[5]->Fill(Mxmass, Mxpt);}}
        if(CentralityVZERO>0  && CentralityVZERO<=10){if(isEgamx){fMassMixedEGA[6]->Fill(Mxmass, Mxpt);}if(isMBmx){fMassMixedMB[6]->Fill(Mxmass, Mxpt);}if(kAllMB){fMassMixedAllMB[6]->Fill(Mxmass, Mxpt);} if(isCentralmx){fMassMixedCent[6]->Fill(Mxmass, Mxpt);}if(isSemiCentralmx){fMassMixedSemiCent[6]->Fill(Mxmass, Mxpt);}}
        if(CentralityVZERO>0  && CentralityVZERO<=20){if(isEgamx){fMassMixedEGA[7]->Fill(Mxmass, Mxpt);}if(isMBmx){fMassMixedMB[7]->Fill(Mxmass, Mxpt);}if(kAllMB){fMassMixedAllMB[7]->Fill(Mxmass, Mxpt);} if(isCentralmx){fMassMixedCent[7]->Fill(Mxmass, Mxpt);}if(isSemiCentralmx){fMassMixedSemiCent[7]->Fill(Mxmass, Mxpt);}}
        if(CentralityVZERO>40 && CentralityVZERO<=50){if(isEgamx){fMassMixedEGA[8]->Fill(Mxmass, Mxpt);}if(isMBmx){fMassMixedMB[8]->Fill(Mxmass, Mxpt);}if(kAllMB){fMassMixedAllMB[8]->Fill(Mxmass, Mxpt);} if(isCentralmx){fMassMixedCent[8]->Fill(Mxmass, Mxpt);}if(isSemiCentralmx){fMassMixedSemiCent[8]->Fill(Mxmass, Mxpt);}}
        if(CentralityVZERO>50 && CentralityVZERO<=60){if(isEgamx){fMassMixedEGA[9]->Fill(Mxmass, Mxpt);}if(isMBmx){fMassMixedMB[9]->Fill(Mxmass, Mxpt);}if(kAllMB){fMassMixedAllMB[9]->Fill(Mxmass, Mxpt);} if(isCentralmx){fMassMixedCent[9]->Fill(Mxmass, Mxpt);}if(isSemiCentralmx){fMassMixedSemiCent[9]->Fill(Mxmass, Mxpt);} }
    }
    //PostData(1, fHistList);
}

//________________________________________________________________________
void AliEMCalpi0ClusterEvaluationTask::Terminate(const Option_t*)
{
    fHistList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fHistList) {
        printf("ERROR: Output list not available\n");
        return;
    }

}
