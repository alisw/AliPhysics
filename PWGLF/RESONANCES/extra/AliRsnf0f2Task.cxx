/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

//==================================================================
// Simple class for f0(980) and f2(1270) analyses. 
// by Beomkyu KIM
//==================================================================

#include "TFile.h"
#include "AliRsnf0f2Task.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliCentrality.h"
Double_t        pi = TMath::Pi();
const char*     tsign[]={"PN","PP","NN"}; //P=Positive charge, N=Negative
Double_t        ptbins[]={0.5,1,2,3,5,7,9,11,15}; //GeV/c
const UInt_t    nptbins = sizeof(ptbins)/sizeof(ptbins[0])-1;
Double_t        centbins[]={ 0,10,20,40,60,80,100}; // Percentile
const UInt_t    ncentbins = sizeof(centbins)/sizeof(centbins[0])-1;

AliRsnf0f2RunTable::AliRsnf0f2RunTable() :
    fCollisionType(kUnknownCollType)
{;}

AliRsnf0f2RunTable::AliRsnf0f2RunTable(Int_t runnumber) 
{
    if (runnumber>=114737 && runnumber<=130850) fCollisionType = kPP; //LHC10bcde
    else if (runnumber>=144871 && runnumber<=146860) fCollisionType=kPP;//LHC11a
    else if (runnumber>=136851 && runnumber<=139517) fCollisionType=kAA;//LHC10h
    else if (runnumber>=167813 && runnumber<=170595) fCollisionType=kAA;//LHC11h
    else if (runnumber>=188356 && runnumber<=188503) fCollisionType=kPA;//LHC12g
    else if (runnumber>=189122 && runnumber<=192732) fCollisionType=kPA;//LHC12h
    else if (runnumber>=195344 && runnumber<=195483) fCollisionType=kPA;//LHC13b
    else if (runnumber>=195529 && runnumber<=195677) fCollisionType=kPA;//LHC13c
    else if (runnumber>=195724 && runnumber<=195872) fCollisionType=kPA;//LHC13d
    else if (runnumber>=195955 && runnumber<=195872) fCollisionType=kPA;//LHC13e
    else if (runnumber>=197669 && runnumber<=200000) fCollisionType=kPA;//LHC13g
    else fCollisionType=kUnknownCollType;
}

//___________________________________________________________________
AliRsnf0f2Task::AliRsnf0f2Task()
:AliAnalysisTaskSE("AliRsnf0f2Task")
    , fOption() 
    , fOutput(NULL) 
    , fTrigger(NULL)
    , fTrackCuts(NULL)
    , fEsd(NULL)
    , fAod(NULL)
    , IsFirstEvent(kTRUE)
    , fRunTable(NULL)
    , fCent(-999)
    , goodtrackindices()
    , fPIDResponse(NULL)
    , fEventNumbers(NULL)
    , fPhiEta(NULL)
    , fMass2D   ( sizeof(tsign)/sizeof(tsign[0]), std::vector<TH2D*> (ncentbins)) 
{
    DefineOutput (1, TList::Class());
}
//___________________________________________________________________
AliRsnf0f2Task::AliRsnf0f2Task
( 
      const char *name
    , const char *option
)
:AliAnalysisTaskSE(name)
    , fOption(option)
    , fOutput(NULL)
    , fTrigger(NULL)
    , fTrackCuts(NULL)
    , fEsd(NULL)
    , fAod(NULL)
    , IsFirstEvent(kTRUE)
    , fRunTable(NULL)
    , fCent(-999)
    , goodtrackindices()
    , fPIDResponse(NULL)
    , fEventNumbers(NULL)
    , fPhiEta(NULL)
    , fMass2D   ( sizeof(tsign)/sizeof(tsign[0]), std::vector<TH2D*> (ncentbins)) 
{
    DefineOutput (1, TList::Class());
}
//___________________________________________________________________
AliRsnf0f2Task::AliRsnf0f2Task
(
      const AliRsnf0f2Task& ap
)
    : fOption(ap.fOption)
    , fOutput(ap.fOutput)
    , fTrigger(ap.fTrigger)
    , fTrackCuts(ap.fTrackCuts)
    , fEsd(ap.fEsd)
    , fAod(ap.fAod)
    , IsFirstEvent(ap.IsFirstEvent)
    , fRunTable(ap.fRunTable)
    , fCent(ap.fCent)
    , goodtrackindices(ap.goodtrackindices)
    , fPIDResponse(ap.fPIDResponse)
    , fEventNumbers(ap.fEventNumbers)
    , fPhiEta(ap.fPhiEta)
    , fMass2D   ( ap.fMass2D) 
{
}
//___________________________________________________________________
AliRsnf0f2Task& AliRsnf0f2Task::operator = 
(
      const AliRsnf0f2Task& ap
)
{
    // assignment operator

    this->~AliRsnf0f2Task();
    new(this) AliRsnf0f2Task(ap);
    return *this;
}
//___________________________________________________________________
AliRsnf0f2Task::~AliRsnf0f2Task()
{
    delete fOutput;
    delete fTrigger;
    delete fTrackCuts;
    delete fPIDResponse;
    delete fEventNumbers;
    delete fRunTable; 
    delete fPhiEta;
    for (UInt_t i=0; i<fMass2D.size(); i++){ // likesign bins
        for (UInt_t j=0; j<fMass2D[i].size(); j++) { // cent bins
            delete fMass2D[i][j];
        }
    }
}

//___________________________________________________________________
void AliRsnf0f2Task::UserCreateOutputObjects()

{
    // Histograms container
    fOutput = new TList();
    fOutput->SetOwner(kTRUE);

    // Offline triggers -----------------------------------------------------
    fTrigger = new AliTriggerAnalysis; // offline trigger
    //-----------------------------------------------------------------------

    // TrackCuts for strangeness measure-------------------------------------
    fTrackCuts = new AliESDtrackCuts();
    {
        fTrackCuts -> GetStandardITSTPCTrackCuts2010(1,0);
        fTrackCuts -> SetEtaRange(-0.8,0.8);
        fTrackCuts -> SetPtRange(0.4, 1e10);
    }

   
    for (UInt_t i=0; i<fMass2D.size(); i++) {
        for (UInt_t j=0; j<fMass2D[i].size(); j++) {
            fMass2D[i][j] = new TH2D(Form("hMassVsPt%sCentBin%d",tsign[i],j)
                ,Form("hMassVsPt%sCentBin_%2.0f_%2.0f",tsign[i]
                ,centbins[j],centbins[j+1]),400,0.2,4.2,145,0.5,15);
        }
    }

    int NBINS=100;
    Float_t LogBins[NBINS+1];
    Float_t LimL=0.001;
    Float_t LimH=7000;
    Float_t logBW = (log(LimH)-log(LimL))/NBINS;
    for(int ij=0;ij<=NBINS;ij++) LogBins[ij]=LimL*exp(ij*logBW);
    
    
    const char *llname[]={"All","PS","PSPileUp","GoodZ","GoodZCut"};
    std::vector <TString> lname (llname, llname + sizeof(llname)/sizeof(llname[0]));
    for (UInt_t i=0; i<ncentbins; i++ ){
        lname.push_back(Form("Cent%2.0f_%2.0f",centbins[i],centbins[i+1]));
    }
    fEventNumbers = new TH1D("EventNumbers","",lname.size(),0,lname.size());
    for (UInt_t i=0; i< lname.size(); i++){
        fEventNumbers -> GetXaxis()-> SetBinLabel(i+1,lname[i]);
    }
    fOutput->Add(fEventNumbers);
    
    //-----------------------------------------------------------
    fPhiEta = new TH2D("hPhiEta","",180,0,2*pi,20,-1,1);
    fOutput -> Add(fPhiEta);

    

    PostData(1, fOutput);
}

//___________________________________________________________________
void AliRsnf0f2Task::UserExec(Option_t* )
{


    // Pointer to a event----------------------------------------------------
    AliVEvent *event = InputEvent();
    if (!event) {
         Printf("ERROR: Could not retrieve event"); 
         return; 
    }
    // ----------------------------------------------------------------------

    // connect to ESD tree --------------------------------------------------
    Int_t runnumber;
    if (event->IsA()==AliESDEvent::Class()){
        fEsd = dynamic_cast<AliESDEvent*>(event);
        if (!fEsd) return;
    }
    if (event->IsA()==AliAODEvent::Class()){
        fAod = dynamic_cast<AliAODEvent*>(event);
        if (!fAod) return;
    }

    if (fEsd) runnumber = fEsd->GetRunNumber();
    else runnumber = fAod->GetRunNumber();
    
    if( IsFirstEvent ) {
        fRunTable = new AliRsnf0f2RunTable(runnumber);
        IsFirstEvent = kFALSE;
    }

    // centrality
    fCent = -999;
    //if(fRunTable->IsHeavyIon() || fRunTable->IsPA()){
    if(fRunTable->IsAA()){
        AliCentrality *cent = event->GetCentrality();
        if( ! cent ) return;
        fCent = cent->GetCentralityPercentile("V0M");
    } else {
        fCent = 0; // Centrality 0-10% bin for pp or pA collisions. 
    }

    

    // Pointer to a MC event-------------------------------------------------
    //AliMCEvent *mcEvent = MCEvent();
    // ----------------------------------------------------------------------


    // Load InputHandler for each event---------------------------------------
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)
        AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
    // -----------------------------------------------------------------------

    fPIDResponse = (AliPIDResponse*) inputHandler->GetPIDResponse();
    if(!fPIDResponse){
        printf("AliRsnf0f2Task No PIDd\n");
    }


    Bool_t IsMinimumBias = kFALSE;
    fEventNumbers -> Fill("All",1);
    Bool_t IsMBOR = kFALSE;
    Bool_t IsMBAND= kFALSE;

    if(fRunTable->IsAA() || fRunTable->IsPA()){
        IsMinimumBias = (inputHandler -> IsEventSelected()) 
            & (AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
    } else {
        IsMinimumBias = inputHandler -> IsEventSelected() & AliVEvent::kMB ;
    }

    if (IsMinimumBias)  fEventNumbers -> Fill("PS",1);
    // Reject pile-up events--------------------------------------------------
    //if (!IsMC && event->IsPileupFromSPD(3.,0.8,3.,2.,5.)) return;
    if (fRunTable->IsPP() && event->IsPileupFromSPD(3.,0.8,3.,2.,5.)) return;
    if (IsMinimumBias) fEventNumbers -> Fill ("PSPileUp",1);
    // -----------------------------------------------------------------------


    AliVVertex* trackVtx = (fAod) ?
        (AliVVertex*)fAod->GetPrimaryVertex() :
        (AliVVertex*)fEsd->GetPrimaryVertexTracks();
    AliVVertex* spdVtx = (fAod) ?
        (AliVVertex*)fAod->GetPrimaryVertexSPD() :
        (AliVVertex*)fEsd->GetPrimaryVertexSPD();


    Bool_t IsGoodVertex = kFALSE;
    Bool_t IsGoodVertexCut = kFALSE;
    if (fabs(trackVtx->GetZ()-spdVtx->GetZ())<0.5) {
        IsGoodVertex = kTRUE; 
        if (IsMinimumBias) fEventNumbers -> Fill("GoodZ",1);
    }
    if (IsGoodVertex && fabs(trackVtx->GetZ())<10.) {
        IsGoodVertexCut = kTRUE;
        if (IsMinimumBias) {
            fEventNumbers -> Fill("GoodZCut",1);
            UInt_t centbin;
            if (FindCentBin(fCent, centbin)){
                fEventNumbers->Fill(Form("Cent%2.0f_%2.0f",centbins[centbin]
                        ,centbins[centbin+1]),1);
            } 
        }
    }

    if (IsMinimumBias && IsGoodVertexCut){
        if (this -> GoodTracksSelection())
            this -> FillTracks();
    }


    PostData(1, fOutput);   
}

Bool_t AliRsnf0f2Task::GoodTracksSelection(){

    UInt_t ntracks = 0;
    if (fEsd) ntracks = fEsd ->GetNumberOfTracks();
    else ntracks = fAod -> GetNumberOfTracks();
    goodtrackindices.clear();

    AliVTrack * track;

    for (UInt_t it = 0; it<ntracks; it++){
        if (fEsd){
            track = (AliESDtrack*) fEsd ->GetTrack(it);
            if (!track) continue;
            if (!fTrackCuts->AcceptTrack((AliESDtrack*) track)) continue;
            fPhiEta->Fill(track->Phi(),track->Eta());
        }
        else {
            track = (AliAODTrack*) fAod ->GetTrack(it);
            if (!track) continue;
            //hybrid track : AOD 086 -> Filter bit 272
            //hybrid track : AOD 160 -> Filter bit 768
            //hybrid track : AOD 145 -> Filter bit 768
            //hybrid track : AOD 115 -> Filter bit 768
            //cut parameters (numbers) will be modified to a secure way.
            if( ! ((AliAODTrack*) track)->TestFilterBit(768)) continue;
            if (fabs(track->Eta())>0.8) continue;
            if (track->Pt()<0.4) continue;
            fPhiEta->Fill(track->Phi(),track->Eta());
        }
        if (GetPID(fPIDResponse, track) != kPion) continue;
        goodtrackindices.push_back(it);
    }
    return goodtrackindices.size();
}

void AliRsnf0f2Task::FillTracks(){
    Double_t        pionmass = AliPID::ParticleMass(AliPID::kPion);
    Double_t        kaonmass = AliPID::ParticleMass(AliPID::kKaon);

    AliVTrack *track1, *track2;
    // charged track, pion, kaon
    TLorentzVector temp1,temp2;
    TLorentzVector vecsum;
    UInt_t ptbin ;
    UInt_t centbin;
    UInt_t ntracks = goodtrackindices.size();
    if (!FindCentBin(fCent,centbin)) return;
    for (UInt_t  it = 0; it < ntracks-1; it++) {
        //for (Int_t itrack = 0; itrack < ntracks; itrack++) {
        track1 = fEsd ? fEsd->GetTrack(goodtrackindices[it]) :
            fAod->GetTrack(goodtrackindices[it]) ;
        if (!track1) continue;
        temp1.SetXYZM(track1->Px(),track1->Py(), track1->Pz(), pionmass);
        for (UInt_t jt = it+1; jt < ntracks; jt++) {
            track2 = fEsd ? fEsd->GetTrack(goodtrackindices[jt]) :
                fAod->GetTrack(goodtrackindices[jt]) ;
            if (!track2) continue;
            temp2.SetXYZM(track2->Px(),track2->Py(), track2->Pz(), pionmass);
            vecsum = temp1+temp2; // two pion vector sum
            if (fabs(vecsum.Rapidity())>0.5) continue; //rapidity cut
            if (track1->Charge()*track2->Charge()==-1){
                fMass2D[kPN][centbin] -> Fill(vecsum.M(),vecsum.Pt());
            }
            if (track1->Charge()==1 && track2->Charge()==1){
                fMass2D[kPP][centbin] -> Fill(vecsum.M(),vecsum.Pt());
            }
            if (track1->Charge()==-1 && track2->Charge()==-1){
                fMass2D[kNN][centbin] -> Fill(vecsum.M(),vecsum.Pt());
            }
        }

    }
}
 

//___________________________________________________________________
void AliRsnf0f2Task::FinishTaskOutput()
{
    OpenFile(1);
    
    TDirectory *cwd = gDirectory;
    TDirectory *nwd;


    for (UInt_t i=0; i<fMass2D.size(); i++){ // sign bins 
        nwd = cwd->mkdir(tsign[i]);
        nwd -> cd();
        for (UInt_t j=0; j<fMass2D[i].size(); j++){ // cent bins
            fMass2D[i][j]->Write();
        }
        cwd -> cd();
    }
}
//___________________________________________________________________
void AliRsnf0f2Task::Terminate(Option_t*)
{

}
//___________________________________________________________________

Int_t AliRsnf0f2Task::GetPID(AliPIDResponse *pid, const AliVTrack *trk){
    if (!pid) return kUnknown; // no pid available

    Double_t sigmas[] ={-999,-999,-999,-999};

    Int_t ipid = kUnknown;
    Double_t lsigma = 3;
    sigmas[kPion] = pid -> NumberOfSigmasTPC(trk,AliPID::kPion);
    sigmas[kKaon] = pid -> NumberOfSigmasTPC(trk,AliPID::kKaon);
    sigmas[kProton] = pid -> NumberOfSigmasTPC(trk,AliPID::kProton);
    sigmas[kElectron] = pid -> NumberOfSigmasTPC(trk,AliPID::kElectron);
    for (int i=0; i<kUnknown; i++){
        if (fabs(sigmas[i]) < lsigma) {
            lsigma = fabs(sigmas[i]);
            ipid = i;
        }
    }

    // derive information, whether tof pid is available
    if (0){
        const Bool_t ka = !(trk->GetStatus() & AliESDtrack::kTOFmismatch);
        const Bool_t kb =  (trk->GetStatus() & AliESDtrack::kTOFpid);
        const Bool_t ktof = ka && kb;
    }

   if (lsigma>3 ) return kUnknown;
   else  return ipid; 

}
Bool_t AliRsnf0f2Task::FindPtBin(Double_t pt, UInt_t & bin){
    for (UInt_t i=0; i<nptbins; i++) 
        if (pt>=ptbins[i] && pt<ptbins[i+1]) {
            bin = i; 
            return kTRUE;
        }
    return kFALSE;
}
Bool_t AliRsnf0f2Task::FindCentBin(Double_t cent, UInt_t & bin){
    for (UInt_t i=0; i<ncentbins; i++) 
        if (cent>=centbins[i] && cent<centbins[i+1]) {
            bin = i; 
            return kTRUE;
        }
    return kFALSE;
}
