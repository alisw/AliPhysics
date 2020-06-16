/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TEventList.h"
#include "TObject.h"
#include "TNamed.h"
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTRDEfficiency.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <typeinfo>

#include "TApplication.h"
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "AliESDTrdTrack.h"
#include "AliESDtrackCuts.h"

#include "AliAODConversionPhoton.h"
#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODTrdTrack.h"

#include "AliTRDgeometry.h"
#include "AliCDBManager.h"
#include <AliGeomManager.h>
#include <TGeoGlobalMagField.h>


class AliAnalysisTRDEfficiency;    // your analysis class
class AliV0ReaderV1;
 // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTRDEfficiency) // classimp: necessary for root

AliAnalysisTRDEfficiency::AliAnalysisTRDEfficiency() : AliAnalysisTaskSE(), 
    fESD(0), fOutputList(0), fHistPt(0), 
    v0reader(0),
    fhBhqu(0), 
    fhpt(0), 
    fhsag1(0),
    fhsag2(0),
    fhlabel(0), 
    
    fhgevent1(0),
    fhgevent2(0),
    fhgevent3(0),
    fhgevent4(0),
    fhgevent5(0),    
    fhgevent6(0),
    fhgevent7(0),
    fhgevent8(0),    
    fhgevent9(0),
    
    fhg(0),
    fhgdghtr(0), 
    fhgtest(0),
    flst(0x0),
    //fPIDResponse(0), 
    online(0)
    //esdfilter(0)
    
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTRDEfficiency::AliAnalysisTRDEfficiency(const char* name) : AliAnalysisTaskSE(name),
    fESD(0), fOutputList(0), fHistPt(0),
    v0reader(0),
    fhBhqu(0), 
    fhpt(0), 
    fhsag1(0), 
    fhsag2(0),
    fhlabel(0), 
    
    fhgevent1(0),
    fhgevent2(0),
    fhgevent3(0),
    fhgevent4(0),
    fhgevent5(0),    
    fhgevent6(0),
    fhgevent7(0),
    fhgevent8(0),    
    fhgevent9(0),
    
    fhg(0),
    fhgdghtr(0), 
    fhgtest(0),
    flst(0x0),
    //fPIDResponse(0), 
    online(0)
    //esdfilter(0)
    
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTRDEfficiency::~AliAnalysisTRDEfficiency()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}


//_______________________________________________________________________
// Bool_t AliAnalysisTRDEfficiency::UserNotify()
//{
//  return kTRUE;
//}

//_____________________________________________________________________________
void AliAnalysisTRDEfficiency::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)

    
    // example of a histogram
    fHistPt = new TH1F("fHistPt", "fHistPt",5000, 0, 5000);       // create your histogram
    
    v0reader = new AliV0ReaderV1();
    v0reader->SetUseOwnXYZCalculation(kTRUE);

    v0reader->SetCreateAODs(kFALSE);// AOD Output

    v0reader->SetUseAODConversionPhoton(kTRUE);

    v0reader->SetProduceV0FindingEfficiency(kFALSE);
    //v0reader->SetEventCuts("00000003");
    //v0reader->SetConversionCuts("06400008d00100001100000000");
    v0reader->SetEventCuts("00010113");
    v0reader->SetConversionCuts("00200009227300008250400000");

    //esdfilter = new AliAnalysisTaskESDfilter();
    //esdfilter->UserCreateOutputObjects();
    
    cout << "before the UserCreateOutputObjects" << endl;
    v0reader->UserCreateOutputObjects();
    cout << "before the notify" << endl;
    //v0reader->
    fhBhqu  = new TH1F("fhBhqu", "fhBhqu",  20000, 0, 20000);
    fhpt    = new TH1F("fhpt","fhpt",       1000, -50, 50);
    fhsag1  = new TH1F("fhsag1", "fhsag1",    200, -10, 10);
    fhsag2  = new TH1F("fhsag2", "fhsag2",    200, -10, 10);
    fhlabel = new TH1F("fhlabel", "fhlabel",10001, -1, 10000);

    Int_t dim1 = 18;
    Int_t bin1[18]   = {1000, 300, 200, 100, 2, 7,      // electron trd pt, pid, sagitta, rating, has trdtrack, # of tracklets
                        1000, 300, 200, 100, 2, 7,      // positron trd pt, pid, sagitta, rating, has trdtrack, # of tracklets
                        2000,2000, 100, 100, 2, 2};        // photon pt, R, eta, phi, in HQU event, in CINT7-T
    Double_t min1[18]={-100, 0, -10, 0, 0, 0,
                       -100, 0, -10, 0, 0, 0,
                       0   , 0, -2, -1, 0, 0};
    Double_t max1[18]={100, 300, 10, 1, 2, 7,
                       100, 300, 10, 1, 2, 7,
                       100, 200, 2,  7, 2, 2};
    fhg = new THnSparseD("fhg", "fhg", dim1, bin1, min1, max1);
                       
    flst = new TObjArray;
    
    
    fOutputList->Add(fHistPt);          // don't forget to add it to the list! the list will be written to file, so if you want
                                        // your histogram in the output file, add it to the list!
    fOutputList->Add(fhBhqu);
    fOutputList->Add(fhpt);
    fOutputList->Add(fhsag1);
    fOutputList->Add(fhsag2);
    fOutputList->Add(fhlabel);
    
    fOutputList->Add(fhg);
    
    //AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    //if (man) {
    //    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    //    if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
    //}
    
    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}


//_____________________________________________________________________________
void AliAnalysisTRDEfficiency::UserExec(Option_t *)
{
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());    // get an event (called fAOD) from the input file
                                                        // there's another event format (ESD) which works in a similar wya
                                                        // but is more cpu/memory unfriendly. for now, we'll stick with aod's
    
    if(!fESD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event
        // example part: i'll show how to loop over the tracks in an event 
        // and extract some information from them which we'll store in a histogram
    if(fESD) {
        if (!TGeoGlobalMagField::Instance()->GetField()) fESD->InitMagneticField();
    }
    TString clss = fESD->GetFiredTriggerClasses();
    if (!clss.Contains("CINT7-B-NOPF-CENT")) return;          // Only grab these min. bias events
    //if (!clss.Contains("CINT7HQU")) return;          // Only grab these min. bias events
    //if (!clss.Contains("CCALIB1-B-NOPF-CENT")) return;          // Only grab these min. bias events
    
    if ( !v0reader->GetReconstructedGammas() ) v0reader->Init();

    //cout << "before the process event" << endl;
    Bool_t processEvent = v0reader->ProcessEvent(fESD, NULL);       // processing event here
    cout << "fired classes " << clss << ", event file " << (fESD->GetHeader())->GetEventNumberInFile() << ", process event " << processEvent << endl;
    AliConvEventCuts *eventcuts = v0reader->GetEventCuts();
    //AliConversionPhotonCuts *photoncuts = v0reader_>GetPhotonCuts();
    cout << "event is selected " << eventcuts->EventIsSelected(fESD,NULL) << endl;
    cout << "event quality     " << eventcuts->GetEventQuality() << endl;
    
    if (!processEvent) return;

    cout << "event file " << (fESD->GetHeader())->GetEventNumberInFile() << ", fired classes " << clss << endl;

    TClonesArray *lst = v0reader->GetReconstructedGammas();
    
    Int_t year=2018;
    //AliCDBManager* man = AliCDBManager::Instance();
    //if (0) {
    //    man->SetDefaultStorage
    //    (Form("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/%d/OCDB/",year));
    //} else {
    //    man->SetDefaultStorage
    //    (Form("alien://folder=/alice/data/%d/OCDB/",year));
    //}
    //man->SetCacheFlag(kTRUE);
    //man->SetRun(fESD->GetRunNumber());
    online = new AliTRDonlineTrackMatching();
    
    //esdfilter->ConvertESDtoAOD();
    
    for (Int_t i = 0; i < lst->GetEntries(); i++){
        
        if (!lst->At(i)) continue;
        AliAODConversionPhoton *photon = (AliAODConversionPhoton*)lst->At(i);  // get photons
        Int_t v0index = photon->GetV0Index();
        Int_t pindex = photon->GetTrackLabelPositive();
        Int_t nindex = photon->GetTrackLabelNegative();
        
        // initializing photon variables
        AliESDtrack *ptrack = fESD->GetTrack(pindex);
        AliESDtrack *ntrack = fESD->GetTrack(nindex);  // according to v0reader info this is the correct one? sure
        AliESDv0 *V0 = fESD->GetV0(v0index);
        Double_t ppt, ppid, psag, prat, ptrd, ptkl = 0;
        Double_t npt, npid, nsag, nrat, ntrd, ntkl = 0;
        Double_t gpt  = photon->GetPhotonPt();
        Double_t gR   = photon->GetConversionRadius();
        Double_t geta = photon->GetPhotonEta();
        Double_t gphi = photon->GetPhotonPhi();
        Double_t ghqu = 0;
        Double_t gcnt = 0;
        
        for (Int_t j = 0; j < fESD->GetNumberOfTrdTracks(); j++){
            AliESDTrdTrack *trdtrack = fESD->GetTrdTrack(j);
            if (!trdtrack) continue;
            if (trdtrack->GetTrackMatch() == (AliVTrack*)ptrack){   // positive track match
                Double_t rating = GetRating(V0, ptrack, trdtrack);
                Double_t lstp[4] = {trdtrack->Pt(), ptrack->Pt(), photon->GetConversionRadius(), rating};
                ptrd = 1;
                ppt  = trdtrack->Pt();
                ppid = trdtrack->GetPID();
                psag = GetSagitta(trdtrack);
                prat = rating;
                ptkl = trdtrack->GetNTracklets();
                //fhgdghtr->Fill(lstp, 1);
            }
            if (trdtrack->GetTrackMatch() == (AliVTrack*)ntrack){   // negative track match
                Double_t rating = GetRating(V0, ntrack, trdtrack);
                Double_t lstn[4] = {trdtrack->Pt(), ntrack->Pt(), photon->GetConversionRadius(), rating};
                ntrd = 1;
                npt  = trdtrack->Pt();
                npid = trdtrack->GetPID();
                nsag = GetSagitta(trdtrack);
                nrat = rating;
                ntkl = trdtrack->GetNTracklets();
                //fhgdghtr->Fill(lstn, 1);
            }
        }
        
        
        if (clss.Contains("HQU")) ghqu = 1;     // basically only HQU positive cases go through
        if (clss.Contains("CINT7-T")) gcnt = 1; // basically only CINT7-T positive cases go through
        
        Double_t lstfhg[18] = {ppt, ppid, psag, prat, ptrd, ptkl,
                               npt, npid, nsag, nrat, ntrd, ntkl,
                               gpt, gR  , geta, gphi, ghqu, gcnt};
        
        fhg->Fill(lstfhg, 1);
        
    }
    
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}

//_____________________________________________________________________________
Bool_t AliAnalysisTRDEfficiency::GetTrackCuts(AliESDtrack *track){
    //
    static AliESDtrackCuts *EsdTrackCuts = 0x0;
    EsdTrackCuts = new AliESDtrackCuts();
    // TPC; clusterCut = 1, cutAcceptanceEdges = kTRUE, removeDistortedRegions = kFALSE
    EsdTrackCuts->AliESDtrackCuts::SetMinNCrossedRowsTPC(70);
    EsdTrackCuts->AliESDtrackCuts::SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    EsdTrackCuts->SetCutGeoNcrNcl(2., 130., 1.5, 0.0, 0.0);  // only dead zone and not clusters per length
    //EsdTrackCuts->AliESDtrackCuts::SetCutOutDistortedRegionsTPC(kTRUE);
    EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterTPC(4);
    EsdTrackCuts->AliESDtrackCuts::SetAcceptKinkDaughters(kFALSE);
    EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
    // ITS; selPrimaries = 1
    EsdTrackCuts->AliESDtrackCuts::SetRequireITSRefit(kTRUE);
    EsdTrackCuts->AliESDtrackCuts::SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                      AliESDtrackCuts::kAny);
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
    EsdTrackCuts->AliESDtrackCuts::SetMaxChi2TPCConstrainedGlobal(36);
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(2);
    EsdTrackCuts->AliESDtrackCuts::SetDCAToVertex2D(kFALSE);
    EsdTrackCuts->AliESDtrackCuts::SetRequireSigmaToVertex(kFALSE);
    EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterITS(36);
    
    return EsdTrackCuts->AcceptTrack(track);
}

Double_t AliAnalysisTRDEfficiency::GetSagitta(AliESDTrdTrack *trdtrack){
    Int_t b = trdtrack->GetB();
    Int_t c = trdtrack->GetC();
    Float_t sag = 0;
        
    // NOTE: yeeted this code from AliAnalysisTaskCheckTRDTriggerJPsi
    //in case of no gtu simulation -> return maximum 0.5
    if(b==0 && c==0) sag = 0.5;
    else{
        Int_t tmp = (((b & 0xfff) << 12) ^ 0x800000) - 0x800000;
        tmp += (c & 0xfff);
        Float_t invPtDev = tmp * 0.000001;
        sag = invPtDev;
    }
    return sag;
}

//_____________________________________________________________________________
Double_t AliAnalysisTRDEfficiency::GetRating(AliESDv0 *v0, AliESDtrack *track, AliESDTrdTrack *trdtrack){
    
    Double_t x, y, z, R;
    v0->GetXYZ(x, y, z);
    //cout << "~~~~ x " << x << "  y  " << y << "  z " << z << endl;
    //cout << "::: xv " << v0->Xv()  << " yv  " << v0->Yv()  << " zv " << v0->Zv()  << endl;
                    
    R = TMath::Sqrt(x*x + y*y);
    //cout << "************ yippy, here it is ************" << endl;
    online->ProcessEvent(fESD, kFALSE, trdtrack->GetLabel());
    online->EstimateTrackDistance(track, trdtrack, fESD->GetMagneticField(), &y, &z); // get ym and zm
    Double_t rating = online->RateTrackMatch(y, z, track->GetSignedPt(), -1*trdtrack->Pt());
    return rating;
}

//_____________________________________________________________________________
void AliAnalysisTRDEfficiency::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
