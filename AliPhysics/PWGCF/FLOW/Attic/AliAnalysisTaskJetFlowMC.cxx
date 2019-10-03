// Simple class to generate toy events which can be fed into the jet finder
//
// Author: Redmer Alexander Bertens, Utrecht University, 2013
// rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl

// root includes
#include "TChain.h"
#include "TList.h"
#include "TClonesArray.h"
#include "TArrayI.h"
#include "TRandom3.h"
#include "TParticle.h"
#include "TVirtualMCDecayer.h"
// aliroot includes
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskJetFlowMC.h"
#include "AliLog.h"
#include "AliPicoTrack.h"

class AliAnalysisTaskJetFlowMC;
using namespace std;

ClassImp(AliAnalysisTaskJetFlowMC)

//_____________________________________________________________________________
AliAnalysisTaskJetFlowMC::AliAnalysisTaskJetFlowMC() : AliAnalysisTaskSE("AliAnalysisTaskJetFlowMC"), fQA(kFALSE), fTracksOutName("JetFlowMC"),  fTracksInName("PicoTrack"), fTracksIn(0),  fTracksOut(0), fReuseTracks(kFALSE), fMult(2200), fCenBin(-1), fCentralityClasses(0), fFuncVn(0), fOutputList(0), fTrackSpectrum(0), fRandomizeEta(kTRUE), fJetSpectrumSF(0), fNoOfSFJets(0), fHistIntV2(0), fHistIntV3(0), fFlowFluctuations(-10), fMaxNumberOfIterations(100), fPsi2(-10), fPsi3(-10), fPrecisionPhi(1e-10), fDetectorType(kFixedEP), fHistSFJetSpectrum(0), fHistSFJetEtaPhi(0), fDecayer(0x0), fDecayerIterations(1), fDecayerResults(0x0) {
    // default constructor for root IO
    for(Int_t i(0); i < 10; i++) {
        fFuncDiffV2[i]                  = 0x0;
        fFuncDiffV3[i]                  = 0x0;
        fHistDiffV2[i]                  = 0x0;
        fHistDiffV3[i]                  = 0x0;
        fHistOriginalSpectrum[i]        = 0x0;
        fHistOriginalEtaPhi[i]          = 0x0;
        fHistToySpectrum[i]             = 0x0;
        fHistToyEtaPhi[i]               = 0x0;
        fHistOriginalDeltaPhi[i]        = 0x0;
        fHistToyDeltaPhi[i]             = 0x0;
        fHistToyVn[i]                   = 0x0;
    }
    for(Int_t i(0); i < 25; i++) fDecayerCache[i] = 0x0;
}
//_____________________________________________________________________________
AliAnalysisTaskJetFlowMC::AliAnalysisTaskJetFlowMC(const char *name, Bool_t qa, Int_t seed) : AliAnalysisTaskSE(name), fQA(qa), fTracksOutName("JetFlowMC"), fTracksInName("PicoTrack"), fTracksIn(0), fTracksOut(0), fReuseTracks(kFALSE), fMult(2200), fCenBin(-1), fCentralityClasses(0), fFuncVn(0), fOutputList(0), fTrackSpectrum(0), fRandomizeEta(kTRUE), fJetSpectrumSF(0), fNoOfSFJets(0), fHistIntV2(0), fHistIntV3(0), fFlowFluctuations(-10), fMaxNumberOfIterations(100), fPsi2(-10), fPsi3(-10), fPrecisionPhi(1e-10), fDetectorType(kFixedEP), fHistSFJetSpectrum(0), fHistSFJetEtaPhi(0), fDecayer(0x0), fDecayerIterations(1), fDecayerResults(0x0) {
    // constructor
    DefineInput(0, TChain::Class());
    if(fQA) DefineOutput(1, TList::Class());
    for(Int_t i(0); i < 10; i++) {
        fFuncDiffV2[i]                  = 0x0;
        fFuncDiffV3[i]                  = 0x0;
        fHistDiffV2[i]                  = 0x0;
        fHistDiffV3[i]                  = 0x0;
        fHistOriginalSpectrum[i]        = 0x0;
        fHistOriginalEtaPhi[i]          = 0x0;
        fHistToySpectrum[i]             = 0x0;
        fHistToyEtaPhi[i]               = 0x0;
        fHistOriginalDeltaPhi[i]        = 0x0;
        fHistToyDeltaPhi[i]             = 0x0;
        fHistToyVn[i]                   = 0x0;
    }
    // by default construction replace gRandom
    if(seed > -1) {
        if(gRandom) delete gRandom;
        gRandom = new TRandom3(seed);
    }
    for(Int_t i(0); i < 25; i++) fDecayerCache[i] = 0x0;
}
//_____________________________________________________________________________
AliAnalysisTaskJetFlowMC::~AliAnalysisTaskJetFlowMC()
{
    // desctructor, claim ownership of deleted objects by setting member pointers to NULL
    if(fOutputList) {
      delete fOutputList;
      fOutputList= NULL;
    }
    if(fDecayer) {
        delete fDecayer;
        fDecayer = NULL;
        for(Int_t i(0); i < 25; i++) {
            if(fDecayerCache[i]) {
                fDecayerCache[i]->Delete();
                delete fDecayerCache[i];
                fDecayerCache[i] = NULL;
            }
        }
        if(fDecayerResults) {
            fDecayerResults->Delete();
            delete fDecayerResults;
            fDecayerResults = NULL;
        }
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlowMC::UserCreateOutputObjects()
{
    // Create my user objects.
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);
    fTracksOut = new TClonesArray("AliPicoTrack");
    fTracksOut->SetName(fTracksOutName);
    if(fQA) {
        fOutputList = new TList();
        fOutputList->SetOwner(kTRUE);
    }
    if(!fCentralityClasses) { // classes must be defined at this point
        Int_t c[] = {-10, 110};
        fCentralityClasses = new TArrayI(sizeof(c)/sizeof(c[0]), c);
    }
    if(fHistIntV2 && !fHistIntV3) {      // define function
        fFuncVn = new TF1("kV2", "[0]*([1]+[2]*[3]*TMath::Cos([2]*(x-[4])))", 0, TMath::TwoPi());
        fFuncVn->SetParameter(0, 1.);        // normalization
        fFuncVn->SetParameter(3, 0.2);       // v2
        fFuncVn->FixParameter(1, 1.);        // constant
        fFuncVn->FixParameter(2, 2.);        // constant
    } else if (!fHistIntV2 && fHistIntV3) {
        fFuncVn = new TF1("kV3", "[0]*([1]+[2]*[3]*TMath::Cos([2]*(x-[4])))", 0, TMath::TwoPi());
        fFuncVn->SetParameter(0, 1.);        // normalization
        fFuncVn->SetParameter(3, 0.2);       // v3
        fFuncVn->FixParameter(1, 1.);        // constant
        fFuncVn->FixParameter(2, 3.);        // constant
    } else if (fHistIntV2 && fHistIntV3) {
        fFuncVn = new TF1("kCombined", "[0]*([1]+[2]*([3]*TMath::Cos([2]*(x-[4]))+[7]*TMath::Cos([5]*(x-[6]))))", 0, TMath::TwoPi());
        fFuncVn->SetParameter(0, 1.);       // normalization
        fFuncVn->SetParameter(3, 0.2);      // v2
        fFuncVn->FixParameter(1, 1.);       // constant
        fFuncVn->FixParameter(2, 2.);       // constant
        fFuncVn->FixParameter(5, 3.);       // constant
        fFuncVn->SetParameter(7, 0.2);      // v3
    }
    // add the generator objects that have been added to the task
    if(fTrackSpectrum && fQA)  fOutputList->Add(fTrackSpectrum);
    if(fJetSpectrumSF && fQA)  fOutputList->Add(fJetSpectrumSF);
    if(fHistIntV2 && fQA)      fOutputList->Add(fHistIntV2);
    if(fHistIntV3 && fQA)      fOutputList->Add(fHistIntV3);
    // create the QA histos
    if(fQA) {
        for(Int_t i(0); i < fCentralityClasses->GetSize()-1; i++) {
            fHistOriginalSpectrum[i] = BookTH1F("fHistOriginalSpectrum", "p_{t} [GeV/c]", 200, 0, 200, i);
            fHistOriginalEtaPhi[i] = BookTH2F("fHistOriginalEtaPhi", "#eta", "#varphi", 100, -1., 1., 100, 0, TMath::TwoPi(), i);
            fHistToySpectrum[i] = BookTH1F("fHistToySpectrum", "p_{t} [GeV/c]", 200, 0, 200, i);
            fHistToyEtaPhi[i] = BookTH2F("fHistToyEtaPhi", "#eta", "#varphi", 100, -1., 1., 100, 0, TMath::TwoPi(), i);
            fHistOriginalDeltaPhi[i] = BookTH1F("fHistOriginalDeltaPhi", "#varphi - #Psi", 100, 0., TMath::Pi(), i);
            fHistToyDeltaPhi[i] = BookTH1F("fHistToyDeltaPhi", "#varphi - #Psi", 100, 0., TMath::Pi(), i);
            fHistToyVn[i] = BookTH2F("fHistToyVn", "p_{t} [GeV/c]", "v_{n}", 1000, 0, 200, 200, 0, .2, i);
            // add to outputlist
            if(fFuncDiffV2[i]) fOutputList->Add(fFuncDiffV2[i]);
            if(fFuncDiffV3[i]) fOutputList->Add(fFuncDiffV3[i]);
            if(fHistDiffV2[i]) fOutputList->Add(fHistDiffV2[i]);
            if(fHistDiffV3[i]) fOutputList->Add(fHistDiffV3[i]);
        }
    }
    if(fJetSpectrumSF && fQA) {
        fHistSFJetSpectrum = BookTH1F("fHistSFJetSpectrum", "p_{t} SF jets [GeV/c]", 100, 0, 200);
        fHistSFJetEtaPhi = BookTH2F("fHistSFJetEtaPhi", "#eta", "#varphi", 100, -1., 1., 100, 0, TMath::TwoPi());
    }
    // reset the random seed
    if(gRandom) {
        delete gRandom;
        gRandom = new TRandom3(0);
    }
    if(!fQA) return;
    // if presente initialize the decayer
    if(fDecayer) {
        fDecayer->Init();
        // this may seem a bit amtibious but will reduce new / delete calls significantly
        // to avoid memory fragmentation
        for(Int_t i(0); i < 25; i++) fDecayerCache[i] = new TClonesArray("TParticle", 10);
        if(fDecayerIterations > 3) {
            printf(" -> Warning, fDecayerIterations set to %i. This number is too high, reverting to 3\n", fDecayerIterations);
            fDecayerIterations = 3;
        }
        fDecayerResults = new TClonesArray("AliPicoTrack", 100);
    }
    // post output
    fOutputList->Sort();
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
TH1F* AliAnalysisTaskJetFlowMC::BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Int_t c, Bool_t append)
{
    // book a TH1F and connect it to the output container
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);    
    if(!fOutputList) return 0x0;
    TString title(name);
    if(c!=-1) { // format centrality dependent histograms accordingly
        name = Form("%s_%i", name, c);
        title += Form("_%i-%i", fCentralityClasses->At(c), fCentralityClasses->At(1+c));
    }
    title += Form(";%s;[counts]", x);
    TH1F* histogram = new TH1F(name, title.Data(), bins, min, max);
    histogram->Sumw2();
    if(append) fOutputList->Add(histogram);
    return histogram;   
}
//_____________________________________________________________________________
TH2F* AliAnalysisTaskJetFlowMC::BookTH2F(const char* name, const char* x, const char*y, Int_t binsx, Double_t minx, Double_t maxx, Int_t binsy, Double_t miny, Double_t maxy, Int_t c, Bool_t append)
{
    // book a TH2F and connect it to the output container
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);    
    if(!fOutputList) return 0x0;
    TString title(name);
    if(c!=-1) { // format centrality dependent histograms accordingly
        name = Form("%s_%i", name, c);
        title += Form("_%i-%i", fCentralityClasses->At(c), fCentralityClasses->At(1+c));
    }
    title += Form(";%s;%s", x, y);
    TH2F* histogram = new TH2F(name, title.Data(), binsx, minx, maxx, binsy, miny, maxy);
    histogram->Sumw2();
    if(append) fOutputList->Add(histogram);
    return histogram;   
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlowMC::UserExec(Option_t *) 
{
    // user exec, called for each event.
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);    
    if(!AliAnalysisManager::GetAnalysisManager()) return;
    // retrieve tracks from input, only necessary when 'reuse' option is set, otherwise tracks will ge generated
    if (fReuseTracks && !fTracksIn) { 
        fTracksIn = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksInName));
        if (!fTracksIn) {
          AliError(Form("Could not retrieve tracks %s!", fTracksInName.Data())); 
          return;
        }
    }
    // get the centrality bin for qa plots
    Double_t cent(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"));
    fCenBin = -1;
    for(Int_t i(0); i < fCentralityClasses->GetSize()-1; i++) {
        if(cent >= fCentralityClasses->At(i) && cent <= fCentralityClasses->At(1+i)) {
            fCenBin = i;
            break; }
    }
    if(fCenBin < 0) return;
    // add tracks to event if not yet there
    fTracksOut->Delete();
    if (!(InputEvent()->FindListObject(fTracksOutName))) InputEvent()->AddObject(fTracksOut);
    fTracksOut->Clear();
    // get the event plane
    CalculateEventPlane();
    Int_t nacc(0);
    if(fReuseTracks) {
        const Int_t Ntracks = fTracksIn->GetEntriesFast();
        for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
            AliPicoTrack* track = static_cast<AliPicoTrack*>(fTracksIn->At(iTracks));
            if(!track) continue;
            Double_t phi(track->Phi()), pt((fTrackSpectrum) ? GetTrackPt() : track->Pt()), eta(fRandomizeEta ? GetTrackEta() : track->Eta());
            // fill qa histo's before applying any (possible) afterburner
            if(fQA) FillHistogramsOriginalData(pt, eta, phi);
            if(fHistDiffV2[fCenBin] || fFuncDiffV2[fCenBin])        V2AfterBurner(phi, eta, pt);
            else if(fHistDiffV3[fCenBin] || fFuncDiffV3[fCenBin])   V3AfterBurner(phi, eta, pt);
            else if(fHistIntV2 || fHistIntV3)                       SampleVnFromTF1(phi);
            if(fDecayer) nacc = InsertDecayDaughters(track, fTracksOut);
            else {
                /*AliPicoTrack *picotrack =*/ new ((*fTracksOut)[nacc]) AliPicoTrack(pt, eta, phi, track->Charge(), track->GetLabel(), 4, track->GetTrackEtaOnEMCal(), track->GetTrackPhiOnEMCal(), track->GetTrackPtOnEMCal(), 1); 
                nacc++;
            }
        }
    } else {
        Double_t pt(0), eta(0), phi(0);
        Short_t charge(0);
        for (Int_t iTracks = 0; iTracks < fMult; ++iTracks) {
            pt = GetTrackPt();
            eta = gRandom->Uniform(-.9, .9);
            phi = gRandom->Uniform(0., TMath::TwoPi());
            charge = (gRandom->Uniform() < 0.5) ? -1 : 1;
            // fill qa histo's before applying any (possible) afterburner
            if(fQA) FillHistogramsOriginalData(pt, eta, phi);
            if(fHistDiffV2[fCenBin] || fFuncDiffV2[fCenBin])        V2AfterBurner(phi, eta, pt);
            else if(fHistDiffV3[fCenBin] || fFuncDiffV3[fCenBin])   V3AfterBurner(phi, eta, pt);
            else if(fHistIntV2 || fHistIntV3)                       SampleVnFromTF1(phi);        
            if(fDecayer) nacc = InsertDecayDaughters(pt, phi, eta, -999, charge, fTracksOut);
            else {
                /*AliPicoTrack *picotrack =*/ new ((*fTracksOut)[nacc]) AliPicoTrack(pt, eta, phi, 1, 0, 4, eta, phi, pt, 1); 
                nacc++;
            }
        }
    }
    if(fJetSpectrumSF && fNoOfSFJets > 0) InjectSingleFragmentationJetSpectrum(nacc);
    if(fQA) PostData(1, fOutputList);
    if(fDebug > 0) PrintInfo();
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlowMC::V2AfterBurner(Double_t &phi, Double_t &eta, Double_t &pt) const
{
    // similar to AliFlowTrackSimple::AddV2, except for the flow fluctuations
    if(fDebug > 1) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);    
    phi = gRandom->Uniform(0, TMath::TwoPi());
    Double_t phi0(phi), v2(GetV2(pt)), f(0.), fp(0.), phiprev(0.);
    if(TMath::AreEqualAbs(v2, 0, 1e-5) && fQA) { 
        FillHistogramsToyData(pt, eta, phi, v2);
        return;
    }
    // introduce flow fluctuations (gaussian)
    if(fFlowFluctuations > -10) GetFlowFluctuation(v2);
    for (Int_t i=0; i!=fMaxNumberOfIterations; ++i) {
        phiprev=phi; //store last value for comparison
        f =  phi-phi0+v2*TMath::Sin(2.*(phi-fPsi2));
        fp = 1.0+2.0*v2*TMath::Cos(2.*(phi-fPsi2)); //first derivative
        phi -= f/fp;
        if (TMath::AreEqualAbs(phiprev,phi,fPrecisionPhi)) break; 
    }
    if(fQA) FillHistogramsToyData(pt, eta, phi, v2);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlowMC::V3AfterBurner(Double_t &phi, Double_t &eta, Double_t &pt) const
{
    // similar to AliFlowTrackSimple::AddV3, except for the flow fluctuations
    if(fDebug > 1) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);    
    phi = gRandom->Uniform(0, TMath::TwoPi());
    Double_t phi0(phi), v3(GetV3(pt)), f(0.), fp(0.), phiprev(0.);
    if(TMath::AreEqualAbs(v3, 0, 1e-5) && fQA) {
        FillHistogramsToyData(pt, eta, phi, v3);
        return;
    }
    // introduce flow fluctuations (gaussian)
    if(fFlowFluctuations > -10) GetFlowFluctuation(v3);
    for (Int_t i=0; i<fMaxNumberOfIterations; i++)  {
        phiprev=phi; //store last value for comparison
        f =  phi-phi0+2./3.*v3*TMath::Sin(3.*(phi-fPsi3));
        fp = 1.0+2.0*v3*TMath::Cos(3.*(phi-fPsi3)); //first derivative
        phi -= f/fp;
        if (TMath::AreEqualAbs(phiprev,phi,fPrecisionPhi)) break;
    }
    if(fQA) FillHistogramsToyData(pt, eta, phi, v3);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlowMC::InjectSingleFragmentationJetSpectrum(Int_t nacc)
{
    // inject single fragmentation jet spectrum to the tclones array, note that emcal params 
    // equal the barrel kinematics to pass the track and jet cuts later on
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);    
    for(Int_t i(nacc); i < (nacc + fNoOfSFJets); i++) {
        Double_t eta(gRandom->Uniform(-.5, .5)), phi(gRandom->Uniform(0, TMath::TwoPi())), pt(fJetSpectrumSF->GetRandom());
        /*AliPicoTrack *picotrack =*/ new ((*fTracksOut)[i]) AliPicoTrack(pt, eta, phi, +1, 0, 0, eta, phi, pt, 0);
        if(fQA) {
            fHistSFJetSpectrum->Fill(pt);
            fHistSFJetEtaPhi->Fill(eta, phi);
        }
        ++i;
    }
    nacc = 0;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlowMC::CalculateEventPlane() {
    // grab the event plane orientation from the AliVEvent header
    if(fDebug > 0) printf("__FILE__ = %s \n __LINE __ %i , __FUNC__ %s \n ", __FILE__, __LINE__, __func__);    
    Double_t a(0), b(0), e(0), f(0);
    switch (fDetectorType) {
        case kVZEROA : {
            fPsi2 = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 8, 2, e, f);
            fPsi3 = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 8, 3, a, b);
            break;
        }
        case kVZEROC : {
            fPsi2 = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 9, 2, e, f);
            fPsi3 = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 9, 3, a, b);
            break;
        }
        case kVZEROComb : {
            fPsi2 = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 10, 2, e, f) ;
            fPsi3 = InputEvent()->GetEventplane()->CalculateVZEROEventPlane(InputEvent(), 10, 3, a, b);
            break;
        }
        case kFixedEP : {
            // for fixed EP the EP is fixed to a constant value
            fPsi2 = 0.;
            fPsi3 = 1.; 
            break;
        }
        default : break;
    }
    // if requested, pass the event plane to the phi distribution
    if(fHistIntV2 && fHistIntV3) {
        fFuncVn->SetParameter(4, fPsi2);        fFuncVn->SetParameter(6, fPsi3);
        Double_t v2(fHistIntV2->GetBinContent(fHistIntV2->GetXaxis()->FindBin(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"))));
        Double_t v3(fHistIntV2->GetBinContent(fHistIntV3->GetXaxis()->FindBin(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"))));
        if(fFlowFluctuations > -10) {
            GetFlowFluctuation(v2);             GetFlowFluctuation(v3);
        }
        fFuncVn->SetParameter(3, v2);           fFuncVn->SetParameter(7, v3);
    } else if (fHistIntV2) {
        fFuncVn->SetParameter(4, fPsi2);
        Double_t v2(fHistIntV2->GetBinContent(fHistIntV2->GetXaxis()->FindBin(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"))));
        if(fFlowFluctuations > -10) GetFlowFluctuation(v2);
        fFuncVn->SetParameter(3, v2);
    } else if (fHistIntV3) {
        fFuncVn->SetParameter(4, fPsi3);
        Double_t v3(fHistIntV3->GetBinContent(fHistIntV2->GetXaxis()->FindBin(InputEvent()->GetCentrality()->GetCentralityPercentile("V0M"))));
        if(fFlowFluctuations > -10) GetFlowFluctuation(v3);
        fFuncVn->SetParameter(3, v3);
    }
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskJetFlowMC::InsertDecayDaughters(
        AliPicoTrack* mother,
        TClonesArray* daughters)
{                
    // wrapper function to return decay daughters
    return InsertDecayDaughters(
            mother->Pt(),
            mother->Phi(),
            mother->Eta(),
            mother->M(),
            mother->Charge(),
            daughters);
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskJetFlowMC::InsertDecayDaughters(
        Double_t pt,
        Double_t phi,
        Double_t eta,
        Double_t mass,
        Short_t charge,
        TClonesArray* daughters)
{                
    // use TVirtualMCDecayer to force decay, daughters will be 
    // temporarily be stored in the TClonesArray as TParticle's
    // but finally appended to the daughters TClonesArray array as AliPicoTracks
    Int_t cacheCounter(0);
    // bookkeeping variable, saves the status of the tclonesarray that is used
    // [0] = unused
    // [i] = correspond's to the i-th decay loop
    // [fDecayerIterations] corresponds to the array that 
    // will eventually be returned (so all final state daughters)
    static Int_t arrayStatus[25];              
    for(Int_t i(0); i < 25; i++) arrayStatus[i] = 0;
    // allow for 'cascaded' decay, i.e. do multiple iterations
    for(Int_t nestLevel(0); nestLevel < fDecayerIterations; nestLevel++) {
        // for the first loop there is no TClonesArray available so grab the
        // function arguments
        if(nestLevel == 0) {
            TLorentzVector pMother(     // 4-vector for mother
                pt*TMath::Cos(phi),     // px
                pt*TMath::Sin(phi),     // py
                pt*TMath::SinH(eta),    // pz
                TMath::Sqrt(mass*mass+pt*pt*TMath::CosH(eta)*TMath::CosH(eta)));        // total energy
            // make sure input array is empty
            fDecayerCache[cacheCounter]->Delete();
            // FIXME we need a fix for the PDG code ... 
            Int_t pdgCode(211); // default to charged pion 
            if(mass == 0.13957) pdgCode = (charge > 0 ) ? 211 : -211;
            // decay the mother and import the daughters
            fDecayer->Decay(pdgCode, &pMother);
            // fill first cache level with imported particles and set the decayer flag
            fDecayer->ImportParticles(fDecayerCache[cacheCounter]);
            arrayStatus[cacheCounter] = 1;        // means that arrayStatus[0] is ready for decay in the next loop
            nestLevel++;
            cacheCounter++;
        } else {
            // in subsequent loops check which arrays need to be decayed
            for(Int_t j(0); j < 25; j++) {
                // check if the array is supposed to be decayed at this nest level
                if(arrayStatus[j] == nestLevel) {
                    // grab each particle from the array 
                    for(Int_t k(0); k < fDecayerCache[j]->GetEntries(); k++) {
                        // 0 is the mother. only decay the mother if there are no daughters
                        if(k == 0 && fDecayerCache[j]->GetEntries() > 1) continue;
                        TParticle* daughter = static_cast<TParticle*>(fDecayerCache[j]->At(k));
                        if(!daughter) continue;
                        TLorentzVector pMother(     // 4-vector for mother
                            daughter->Px(),
                            daughter->Py(),
                            daughter->Pz(),
                            daughter->Energy());        // total energy
                        // make sure input array is empty
                        fDecayerCache[cacheCounter]->Delete();
                        // FIXME we need a fix for the PDG code ... 
                        // decay the mother and import the daughters
                        fDecayer->Decay(daughter->GetPdgCode(), &pMother);
                        fDecayer->ImportParticles(fDecayerCache[cacheCounter]);
                        // update the cache value
                        cacheCounter++;
                    }
                }
                // update the next level
                nestLevel++;
            }
        }
    }
    // what is left now is appending the freshly created tracks to the final track array
    Int_t offset(daughters->GetEntries());
    for(Int_t i(0); i < 25; i++) {
        // only take the final state arrays
        if(arrayStatus[i] == fDecayerIterations-1) {
            for(Int_t j(0); j < fDecayerCache[j]->GetEntries(); j++) {
                // if a mother has daughters, only push the daughters
                if(j == 0 && fDecayerCache[j]->GetEntries() > 1) continue;
                TParticle* daughter = static_cast<TParticle*>(fDecayerCache[j]->At(j));
                if(daughter) /*AliPicoTrack *picotrack =*/ new ((*daughters)[offset]) AliPicoTrack(daughter->Pt(), daughter->Eta(), daughter->Phi(), (daughter->GetPdgCode() > 0) ? 1 : -1, 0, 0, daughter->Eta(), daughter->Phi(), daughter->Pt(), 0);
                offset++;
            }
        }
    }
    return offset;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlowMC::Terminate(Option_t *)
{
    // terminate
}
//_____________________________________________________________________________
void AliAnalysisTaskJetFlowMC::PrintInfo() const
{
    // print info
    printf(" > No of retrieved tracks from %s \n \t %i \n", fTracksInName.Data(), fTracksIn->GetEntriesFast());
    printf(" > No of created tracks in %s \n \t %i \n", fTracksOutName.Data(), fTracksOut->GetEntriesFast());

}
//_____________________________________________________________________________
