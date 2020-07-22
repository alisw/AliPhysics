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
    fESD(0), 
    fOutputList(0),
    fV0Reader(0),
    fV0ReaderName("V0ReaderV1"),

    fHistGamma(0),
    fHistPi0(0),
    fHistPi0bkg(0),
    fHistEvent(0),
    fHistEventTrigger(0),
    
    fPIDResponse(0), 
    //esdfilter(0)
    fOnline(0),
    fConvsel(0)
    
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTRDEfficiency::AliAnalysisTRDEfficiency(const char* name) : AliAnalysisTaskSE(name),
    fESD(0), 
    fOutputList(0),
    fV0Reader(0),
    fV0ReaderName("V0ReaderV1"),

    fHistGamma(0),
    fHistPi0(0),
    fHistPi0bkg(0),
    fHistEvent(0),
    fHistEventTrigger(0),
    
    fPIDResponse(0), 
    //esdfilter(0)
    fOnline(0),
    fConvsel(0)
    
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
    //fHistPt = new TH1F("fHistPt", "fHistPt",5000, 0, 5000);       // create your histogram
    //cout << "before everything" << endl;
    //TString fV0ReaderName("V0ReaderV1");
    //TString fV0ReaderName("V0ReaderV1_00000003_00200009227300008250400000");
    //fV0Reader = new AliV0ReaderV1();
    
    fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
    
    //cout << "######## Before the v0reader checks #########" << endl;
	if(!fV0Reader){ printf("Error: No V0 Reader");return;} // GetV0Reader
	
	//cout << "??????????  after the return  ???????????" << endl;
    //fV0Reader->SetUseOwnXYZCalculation(kTRUE);
    //fV0Reader->SetCreateAODs(kFALSE);// AOD Output
    //fV0Reader->SetUseAODConversionPhoton(kTRUE);
    //fV0Reader->SetProduceV0FindingEfficiency(kFALSE);
    //fV0Reader->SetEventCuts("00000003");
    //fV0Reader->SetConversionCuts("06400008d00100001100000000");
    //fV0Reader->SetEventCuts("00010113");
    //fV0Reader->SetConversionCuts("00200009227300008250400000");
    //cout << "~~~~~~~~~~~~~ After the v0reader ~~~~~~~~~~~~~~~" << endl;
    //esdfilter = new AliAnalysisTaskESDfilter();
    //esdfilter->UserCreateOutputObjects();
    //fConvsel = new AliConversionSelection("00010113", "00200009227300008250400000", "0152101500000000");
    
    TString eventcuts = (fV0Reader->GetEventCuts())->GetCutNumber();
    TString photoncuts = (fV0Reader->GetConversionCuts())->GetCutNumber();
    fConvsel = new AliConversionSelection(eventcuts.Data(), photoncuts.Data(), "0152101500000000");
    fV0Reader->UserCreateOutputObjects();

    fOnline = new AliTRDonlineTrackMatching();
    
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
    fHistGamma = new THnSparseD("fHistGamma", "fHistGamma", dim1, bin1, min1, max1);
              
    // pi0
    Int_t dim2 = 38;
    Int_t bin2[38]   ={1000, 300, 200, 100, 2, 7,      // electron trd pt, pid, sagitta, rating, has trdtrack, # of tracklets
                        1000, 300, 200, 100, 2, 7,      // positron trd pt, pid, sagitta, rating, has trdtrack, # of tracklets
                        2000,2000, 100, 100,
                        1000, 300, 200, 100, 2, 7,      // electron trd pt, pid, sagitta, rating, has trdtrack, # of tracklets
                        1000, 300, 200, 100, 2, 7,      // positron trd pt, pid, sagitta, rating, has trdtrack, # of tracklets
                        2000,2000, 100, 100,            // photon pt, R, eta, phi,
                        2, 2, 100, 200, 200, 400};      // HQU, CINT7-T, M, pT, R, z
    Double_t min2[38]={-100, 0, -10, 0, 0, 0,
                       -100, 0, -10, 0, 0, 0,
                       0   , 0, -2, -1,
                       -100, 0, -10, 0, 0, 0,
                       -100, 0, -10, 0, 0, 0,
                       0   , 0, -2, -1,
                       0, 0, 0, 0, 0, -200};
    Double_t max2[38]={100, 300, 10, 1, 2, 7,
                       100, 300, 10, 1, 2, 7,
                       100, 200, 2,  7,
                       100, 300, 10, 1, 2, 7,
                       100, 300, 10, 1, 2, 7,
                       100, 200, 2,  7,
                       2, 2, 1, 200, 200, 200};
        
    fHistPi0    = new THnSparseD("fHistPi0", "fHistPi0", dim2, bin2, min2, max2);
    fHistPi0bkg = new THnSparseD("fHistPi0bkg", "fHistPi0bkg", dim2, bin2, min2, max2);

    // event
    Int_t dim3 = 12;
    Int_t bin3[12]   = {10, 30, 20, 2, 2, 2,         // All events,  # of photons, # of pi0s, extra,...,extra
                        2, 2, 2, 20, 10, 0};         // 
    Double_t min3[12]={0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0};
    Double_t max3[12]={10, 30, 20, 2, 2, 2,
                       2, 2, 2, 20, 10, 0};
    fHistEvent         = new THnSparseD("fHistEvent", "fHistEvent", dim3, bin3, min3, max3);
    fHistEventTrigger  = new TH1F("fHistEventTrigger", "fHistEventTrigger", 100, 0 , 100);
    
    
    //fOutputList->Add(fHistPt);          // don't forget to add it to the list! the list will be written to file, so if you want
                                        // your histogram in the output file, add it to the list!

    fOutputList->Add(fHistGamma);
    fOutputList->Add(fHistPi0);
    fOutputList->Add(fHistPi0bkg);
    fOutputList->Add(fHistEvent);
    fOutputList->Add(fHistEventTrigger);
    
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if (man) {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
    }
    
    if(fV0Reader)
		if((AliConvEventCuts*)fV0Reader->GetEventCuts())
			if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
				fOutputList->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
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
    
    if ( !fV0Reader->GetReconstructedGammas() ) fV0Reader->Init();

    Bool_t processEvent = fV0Reader->ProcessEvent(fESD, NULL);       // processing event here
    //AliConvEventCuts *eventcuts = fV0Reader->GetEventCuts();
    
    if (!processEvent) return;
    TClonesArray *lst = fV0Reader->GetReconstructedGammas();
    
    fConvsel->ProcessEvent(lst, fESD, NULL);    
    Double_t photons = lst->GetEntries();
    Double_t pi0s = fConvsel->GetNumberOfPi0s();
    //esdfilter->ConvertESDtoAOD();
    Int_t pev = 0;
    Int_t nev = 0;
    //cout << "checking if this is running" << endl;
    //cout << "Event cuts  " << (fV0Reader->GetEventCuts())->GetCutNumber() << endl;
    //cout << "Photon cuts " << (fV0Reader->GetConversionCuts())->GetCutNumber() << endl;
    /************    Photons   ************/
    for (Int_t i = 0; i < lst->GetEntries(); i++){
        
        if (!lst->At(i)) continue;
        AliAODConversionPhoton *photon = (AliAODConversionPhoton*)lst->At(i);  // get photons
        Int_t v0index = photon->GetV0Index();
        Int_t pindex = photon->GetTrackLabelPositive();
        Int_t nindex = photon->GetTrackLabelNegative();
        
        // initializing photon variables
        AliESDtrack *ptrack = fESD->GetTrack(pindex);
        AliESDtrack *ntrack = fESD->GetTrack(nindex);  // according to fV0Reader info this is the correct one? sure
        AliESDv0 *V0 = fESD->GetV0(v0index);
        Double_t ppt = 0, ppid = 0, psag = 0, prat = 0, ptrd = 0, ptkl = 0;
        Double_t npt = 0, npid = 0, nsag = 0, nrat = 0, ntrd = 0, ntkl = 0;
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
                Int_t fLayerMask = trdtrack->GetLayerMask();
                ptrd = (fLayerMask >> 0) & 1;
                ppt  = trdtrack->Pt();
                ppid = trdtrack->GetPID();
                psag = GetSagitta(trdtrack);
                prat = rating;
                ptkl = trdtrack->GetNTracklets();
            }
            if (trdtrack->GetTrackMatch() == (AliVTrack*)ntrack){   // negative track match
                Double_t rating = GetRating(V0, ntrack, trdtrack);
                Int_t fLayerMask = trdtrack->GetLayerMask();
                ntrd = (fLayerMask >> 0) & 1;
                npt  = trdtrack->Pt();
                npid = trdtrack->GetPID();
                nsag = GetSagitta(trdtrack);
                nrat = rating;
                ntkl = trdtrack->GetNTracklets();
            }
            if (ppt || npt){   // TODO: check this
                Int_t tmp = GetEventCuts(trdtrack, clss);
                if (tmp > pev){
                    pev = tmp;
                }
                if (tmp > nev){
                    nev = tmp;
                }
            }
        }
        
        
        if (clss.Contains("HQU")) ghqu = 1;     // basically only HQU positive cases go through
        if (clss.Contains("CINT7-T")) gcnt = 1; // basically only CINT7-T positive cases go through
        
        Double_t lstfhg[18] = {ppt, ppid, psag, prat, ptrd, ptkl,
                               npt, npid, nsag, nrat, ntrd, ntkl,
                               gpt, gR  , geta, gphi, ghqu, gcnt};
        
        fHistGamma->Fill(lstfhg, 1);
        
    }
    
    Int_t top = 0;
    if (pev < nev) top = nev;
    else top = pev;
    for (Double_t i = 0; i < top+1; i++){
        Double_t lstevent2[12]={i, photons, pi0s, 0., 0., 0., 0., 0., 0., 0., 0., 0.};
        fHistEvent->Fill(lstevent2, 1);
    }
    
    /**********  pi0  **********/
    for (Int_t i = 0; i < fConvsel->GetNumberOfPi0s(); i++){
        AliAODConversionMother *pi0 = fConvsel->GetPi0(i);
        if (!pi0) continue;
        
        Double_t tmp0[16] = {0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0};
        checkPi0(lst, pi0, tmp0, 0);
                             
        Double_t tmp1[16] = {0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0};
        checkPi0(lst, pi0, tmp1, 1);
        
        Double_t cnt=clss.Contains("CINT7-T");
        Double_t hqu=clss.Contains("HQU");
        Double_t M = pi0->M();
        Double_t pt= pi0->Pt();
        Double_t R = pi0->GetProductionRadius();
        Double_t z = pi0->GetProductionZ();
        //lstpi0 = {};
        Double_t lstfhpi0[38] ={tmp0[0], tmp0[1], tmp0[2], tmp0[3], tmp0[4], tmp0[5],
                                tmp0[6], tmp0[7], tmp0[8], tmp0[9], tmp0[10],tmp0[11],
                                tmp0[12],tmp0[13],tmp0[14],tmp0[15],
                                
                                tmp1[0], tmp1[1], tmp1[2], tmp1[3], tmp1[4], tmp1[5],
                                tmp1[6], tmp1[7], tmp1[8], tmp1[9], tmp1[10],tmp1[11],
                                tmp1[12],tmp1[13],tmp1[14],tmp1[15],
                                
                                cnt, hqu, M, pt, R, z};
        fHistPi0->Fill(lstfhpi0, 1);
    }
    
    /********************  background pi0  ******************/ // not working
    for (Int_t i = 0; i < fConvsel->GetNumberOfBGs(); i++){
        AliAODConversionMother *bg = fConvsel->GetBG(i);
        if (!bg) continue;
        
        Double_t tmp0[16] = {0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0};
        checkPi0(lst, bg, tmp0, 0);
                             
        Double_t tmp1[16] = {0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0};
        checkPi0(lst, bg, tmp1, 1);
        
        Double_t cnt=clss.Contains("CINT7-T");
        Double_t hqu=clss.Contains("HQU");
        Double_t M = bg->M();
        Double_t pt= bg->Pt();
        Double_t R = bg->GetProductionRadius();
        Double_t z = bg->GetProductionZ();
        //lstpi0 = {};
        Double_t lstfhpi0bkg[38]={tmp0[0], tmp0[1], tmp0[2], tmp0[3], tmp0[4], tmp0[5],
                                  tmp0[6], tmp0[7], tmp0[8], tmp0[9], tmp0[10],tmp0[11],
                                  tmp0[12],tmp0[13],tmp0[14],tmp0[15],
                                
                                  tmp1[0], tmp1[1], tmp1[2], tmp1[3], tmp1[4], tmp1[5],
                                  tmp1[6], tmp1[7], tmp1[8], tmp1[9], tmp1[10],tmp1[11],
                                  tmp1[12],tmp1[13],tmp1[14],tmp1[15],
                            
                                  cnt, hqu, M, pt, R, z};
        fHistPi0bkg->Fill(lstfhpi0bkg, 1);
    }
    
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}

//_____________________________________________________________________________
Bool_t AliAnalysisTRDEfficiency::checkPi0(TClonesArray *lst, AliAODConversionMother *pi0, Double_t tmp[16], Int_t lbl){
    
    Int_t label = pi0->GetLabel(lbl);
    AliAODConversionPhoton *photon = (AliAODConversionPhoton*)lst->At(label);  // get photons
    if (!photon) return kFALSE;
    
    Int_t v0index = photon->GetV0Index();
    Int_t pindex = photon->GetTrackLabelPositive();
    Int_t nindex = photon->GetTrackLabelNegative();
    
    // let's try it (is this correct)
    AliESDtrack *ptrack = fESD->GetTrack(pindex);
    AliESDtrack *ntrack = fESD->GetTrack(nindex);  // according to fV0Reader info this is the correct one? sure
    AliESDv0 *V0 = fESD->GetV0(v0index);
    Double_t ppt = 0; 
    Double_t ppid = 0; 
    Double_t psag = 0;
    Double_t prat = 0; 
    Double_t ptrd = 0;
    Double_t ptkl = 0;
    Double_t npt = 0; 
    Double_t npid = 0;
    Double_t nsag = 0;
    Double_t nrat = 0;
    Double_t ntrd = 0;
    Double_t ntkl = 0;
    Double_t gpt  = photon->GetPhotonPt();
    Double_t gR   = photon->GetConversionRadius();
    Double_t geta = photon->GetPhotonEta();
    Double_t gphi = photon->GetPhotonPhi();        

    for (Int_t j = 0; j < fESD->GetNumberOfTrdTracks(); j++){
        AliESDTrdTrack *trdtrack = fESD->GetTrdTrack(j);
        if (!trdtrack) continue;
        if (trdtrack->GetTrackMatch() == (AliVTrack*)ptrack){   // positive track match
            Int_t fLayerMask = trdtrack->GetLayerMask();
            Double_t rating = GetRating(V0, ptrack, trdtrack);
            ptrd = (fLayerMask >> 0) & 1;       // does this give a boolean to track being in the 0th layer
            ppt  = trdtrack->Pt();
            ppid = trdtrack->GetPID();
            psag = GetSagitta(trdtrack);
            prat = rating;
            ptkl = trdtrack->GetNTracklets();
        }
        if (trdtrack->GetTrackMatch() == (AliVTrack*)ntrack){   // negative track match
            Int_t fLayerMask = trdtrack->GetLayerMask();
            Double_t rating = GetRating(V0, ntrack, trdtrack);
            ntrd = (fLayerMask >> 0) & 1;
            npt  = trdtrack->Pt();
            npid = trdtrack->GetPID();
            nsag = GetSagitta(trdtrack);
            nrat = rating;
            ntkl = trdtrack->GetNTracklets();
        }
    }
    tmp[0] = ppt;
    tmp[1] = ppid;
    tmp[2] = psag;
    tmp[3] = prat;
    tmp[4] = ptrd;
    tmp[5] = ptkl;
    
    tmp[6] = npt;
    tmp[7] = npid;
    tmp[8] = nsag;
    tmp[9] = nrat;
    tmp[10]= ntrd;
    tmp[11]= ntkl;
                         
    tmp[12]= gpt;
    tmp[13]= gR;
    tmp[14]= geta;
    tmp[15]= gphi;
                         
    return tmp;
}

Int_t AliAnalysisTRDEfficiency::GetEventCuts(AliESDTrdTrack *trdtrack, TString clss){
    
    if (trdtrack->Pt() < 2)         // pt cut
        return 1;  // only trd track
    
    if (TMath::Abs(GetSagitta(trdtrack))>0.2) // sag cut
        return 2; // pass pt
    
    if (trdtrack->GetPID() < 164)   // pid cut
        return 3; // pass sag

    if (trdtrack->GetNTracklets() < 5)// # of tracklet cuts
        return 4; // pass pid

    Int_t fLayerMask = trdtrack->GetLayerMask(); // layer 0 cuts
    if (( (fLayerMask >> 0) & 1) == 0 )
        return 5; // 
    
    if (!clss.Contains("CINT7-T"))
        return 6;
    
    if (!clss.Contains("HQU"))
        return 7;
    
    return 8;
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
    
    Double_t x = 0, y = 0, z = 0;
    v0->GetXYZ(x, y, z);    
    //R = TMath::Sqrt(x*x + y*y);
    fOnline->ProcessEvent(fESD, kFALSE, trdtrack->GetLabel());
    fOnline->EstimateTrackDistance(track, trdtrack, fESD->GetMagneticField(), &y, &z); // get ym and zm
    Double_t rating = fOnline->RateTrackMatch(y, z, track->GetSignedPt(), -1*trdtrack->Pt());
    return rating;
}

//_____________________________________________________________________________
void AliAnalysisTRDEfficiency::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
