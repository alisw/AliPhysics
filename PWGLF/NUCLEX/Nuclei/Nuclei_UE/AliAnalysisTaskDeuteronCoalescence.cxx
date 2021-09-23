#include "AliAnalysisTaskDeuteronCoalescence.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "TLorentzVector.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;
ClassImp(AliAnalysisTaskDeuteronCoalescence)

//____________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDeuteronCoalescence::AliAnalysisTaskDeuteronCoalescence():
AliAnalysisTaskSE(),
fESDeventCuts(),
fESDevent(nullptr),
fMCEvent(nullptr),
fMCstack(nullptr),
fMCEventHandler(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fAverage_Nch_Transv(7.2),
hProtonWeights(nullptr),
hNumberOfEvents(nullptr),
hTransverseMult(nullptr),
hRtDistribution(nullptr),
hProtonsINELgtZERO(nullptr),
hProtonsINELgtZERO_reshaped(nullptr),
hProtons_Toward(nullptr),
hProtons_Transv(nullptr),
hProtons_Away(nullptr),
hNeutronsINELgtZERO(nullptr),
hNeutronsINELgtZERO_reshaped(nullptr),
hNeutrons_Toward(nullptr),
hNeutrons_Transv(nullptr),
hNeutrons_Away(nullptr),
hRparticles(nullptr),
hRapidityProtons(nullptr),
hRapidityNeutrons(nullptr)
{}
//____________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDeuteronCoalescence::AliAnalysisTaskDeuteronCoalescence(const char *name):
AliAnalysisTaskSE(name),
fESDeventCuts(),
fESDevent(nullptr),
fMCEvent(nullptr),
fMCstack(nullptr),
fMCEventHandler(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fAverage_Nch_Transv(7.2),
hProtonWeights(nullptr),
hNumberOfEvents(nullptr),
hTransverseMult(nullptr),
hRtDistribution(nullptr),
hProtonsINELgtZERO(nullptr),
hProtonsINELgtZERO_reshaped(nullptr),
hProtons_Toward(nullptr),
hProtons_Transv(nullptr),
hProtons_Away(nullptr),
hNeutronsINELgtZERO(nullptr),
hNeutronsINELgtZERO_reshaped(nullptr),
hNeutrons_Toward(nullptr),
hNeutrons_Transv(nullptr),
hNeutrons_Away(nullptr),
hRparticles(nullptr),
hRapidityProtons(nullptr),
hRapidityNeutrons(nullptr)
{
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//____________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDeuteronCoalescence::~AliAnalysisTaskDeuteronCoalescence()  {
    
    fOutputList->Clear();
    delete fESDevent;
    delete fMCEvent;
    delete fMCstack;
    delete fMCEventHandler;
    delete fOutputList;
    delete fQAList;
    delete hProtonWeights;
    delete hNumberOfEvents;
    delete hTransverseMult;
    delete hRtDistribution;
    delete hProtonsINELgtZERO;
    delete hProtonsINELgtZERO_reshaped;
    delete hProtons_Toward;
    delete hProtons_Transv;
    delete hProtons_Away;
    delete hNeutronsINELgtZERO;
    delete hNeutronsINELgtZERO_reshaped;
    delete hNeutrons_Toward;
    delete hNeutrons_Transv;
    delete hNeutrons_Away;
    delete hRparticles;
    delete hRapidityProtons;
    delete hRapidityNeutrons;
}
//____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronCoalescence::UserCreateOutputObjects()  {
    
    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner();
    fQAList     -> SetOwner();
    
    //Event Selection
    fESDeventCuts.AddQAplotsToList(fQAList);
    fESDeventCuts.SetManualMode();
    fESDeventCuts.fRequireTrackVertex = kTRUE;
    fESDeventCuts.fMinVtz = -10.f;
    fESDeventCuts.fMaxVtz = +10.f;
    fESDeventCuts.fMaxDeltaSpdTrackAbsolute = 0.5f;
    fESDeventCuts.fMaxResolutionSPDvertex = 0.25f;
    fESDeventCuts.fTriggerMask = (AliVEvent::kINT7);
    fESDeventCuts.fRejectDAQincomplete = kTRUE;
    fESDeventCuts.fSPDpileupMinContributors = 3;
    fESDeventCuts.fSPDpileupMinZdist = 0.8;
    fESDeventCuts.fSPDpileupNsigmaZdist = 3.0;
    fESDeventCuts.fSPDpileupNsigmaDiamXY = 2.0;
    fESDeventCuts.fSPDpileupNsigmaDiamZ = 5.0;
    fESDeventCuts.fTrackletBGcut = kTRUE;
    
    //Event Counter
    hNumberOfEvents = new TH1D ("hNumberOfEvents","",20,0,20);
    hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    
    //Transverse Multiplicity
    hTransverseMult = new TH1D ("hTransverseMult","",200,0,200);
    hTransverseMult -> Sumw2();
    fOutputList -> Add(hTransverseMult);
   
    //R_{T} Distribution
    hRtDistribution = new TH1D ("hRtDistribution","",1000,0,10);
    hRtDistribution -> Sumw2();
    fOutputList -> Add(hRtDistribution);

    //p_{T} Intervals
    Double_t pt_proton[] = {0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,
        1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,10.0};
    Double_t pt_deuteron[]={0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8};
    const Int_t nPtProton   = sizeof(pt_proton)/sizeof(Double_t)-1;
    const Int_t nPtDeuteron = sizeof(pt_deuteron)/sizeof(Double_t)-1;
    
    //p_{T} Spectra in Events INEL>0
    hProtonsINELgtZERO            = new TH1D ("hProtonsINELgtZERO","",nPtProton,pt_proton);
    hProtonsINELgtZERO_reshaped   = new TH1D ("hProtonsINELgtZERO_reshaped","",nPtProton,pt_proton);
    hNeutronsINELgtZERO           = new TH1D ("hNeutronsINELgtZERO","",nPtProton,pt_proton);
    hNeutronsINELgtZERO_reshaped  = new TH1D ("hNeutronsINELgtZERO_reshaped","",nPtProton,pt_proton);
    hProtonsINELgtZERO            -> Sumw2();
    hNeutronsINELgtZERO           -> Sumw2();
    hProtonsINELgtZERO_reshaped   -> Sumw2();
    hNeutronsINELgtZERO_reshaped  -> Sumw2();
    fOutputList -> Add(hProtonsINELgtZERO);
    fOutputList -> Add(hNeutronsINELgtZERO);
    fOutputList -> Add(hProtonsINELgtZERO_reshaped);
    fOutputList -> Add(hNeutronsINELgtZERO_reshaped);

    
    //p_{T} Spectra Deuterons
    for (Int_t i=0 ; i<10 ; i++)  {
        
        hDeuteronsINELgtZERO[i] = new TH1D (Form("hDeuteronsINELgtZERO[%d]",i),"",nPtDeuteron,pt_deuteron);
        hDeuterons_Toward[i]    = new TH1D (Form("hDeuterons_Toward[%d]",i),"",400,0,4);
        hDeuterons_Transv[i]    = new TH1D (Form("hDeuterons_Transv[%d]",i),"",400,0,4);
        hDeuterons_Away[i]      = new TH1D (Form("hDeuterons_Away[%d]",i),"",400,0,4);
        hDeuteronsINELgtZERO[i] -> Sumw2();
        hDeuterons_Toward[i]    -> Sumw2();
        hDeuterons_Transv[i]    -> Sumw2();
        hDeuterons_Away[i]      -> Sumw2();
        fOutputList -> Add(hDeuteronsINELgtZERO[i]);
        fOutputList -> Add(hDeuterons_Toward[i]);
        fOutputList -> Add(hDeuterons_Transv[i]);
        fOutputList -> Add(hDeuterons_Away[i]);
    }

    //R_{T} Intervals
    Double_t Rt_Intervals[]={0.0,0.5,1.5,5.0};
    Int_t nRtIntervals=sizeof(Rt_Intervals)/sizeof(Double_t)-1;
    
    //Proton Spectra in the Azimuthal Regions
    hProtons_Toward   = new TH1D ("hProtons_Toward","",400,0,4);
    hProtons_Transv   = new TH1D ("hProtons_Transv","",400,0,4);
    hProtons_Away     = new TH1D ("hProtons_Away","",400,0,4);
    hProtons_Toward -> Sumw2();
    hProtons_Transv -> Sumw2();
    hProtons_Away   -> Sumw2();
    fOutputList -> Add(hProtons_Toward);
    fOutputList -> Add(hProtons_Transv);
    fOutputList -> Add(hProtons_Away);
    
    //Neutron Spectra in the Azimuthal Regions
    hNeutrons_Toward   = new TH1D ("hNeutrons_Toward","",400,0,4);
    hNeutrons_Transv   = new TH1D ("hNeutrons_Transv","",400,0,4);
    hNeutrons_Away     = new TH1D ("hNeutrons_Away","",400,0,4);
    hNeutrons_Toward -> Sumw2();
    hNeutrons_Transv -> Sumw2();
    hNeutrons_Away   -> Sumw2();
    //fOutputList -> Add(hNeutrons_Toward);
    //fOutputList -> Add(hNeutrons_Transv);
    //fOutputList -> Add(hNeutrons_Away);

    
    //QA Histograms & Debug
    hRparticles = new TH1D ("hRparticles","",1000,0,1e-5);
    hRparticles -> Sumw2();
    fOutputList -> Add(hRparticles);
    
    //Rapidity Distributions of Protons & Neutrons that form Deuteron in |y|<0.5
    hRapidityProtons  = new TH1D ("hRapidityProtons","",200,-2,2);
    hRapidityNeutrons = new TH1D ("hRapidityNeutrons","",200,-2,2);
    hRapidityProtons  -> Sumw2();
    hRapidityNeutrons -> Sumw2();
    fOutputList -> Add(hRapidityProtons);
    fOutputList -> Add(hRapidityNeutrons);


    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronCoalescence::UserExec(Option_t *)  {
    
    
    //Coalescence Momentum
    Double_t p0[]={0.230,0.235,0.240,0.245,0.250,0.255,0.260,0.265,0.270};
    const Int_t nTrials = sizeof(p0)/sizeof(Double_t);
    
    //Get Input Event (INEL>0 Selection)
    if (!GetEvent()) return;
   
    //Protons and Neutrons IDs
    vector<Int_t> proton_ID;
    vector<Int_t> neutron_ID;
    vector<Int_t> neutron_status;

    //Loop over Generated Particles
    for ( Int_t i=0; i<fMCEvent->GetNumberOfTracks(); i++ )  {

        //Get Primary Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(i);
        if (!particle) continue;
        if (!particle->IsPhysicalPrimary()) continue;
        if (IsInjectedParticle(particle))   continue;
        if (TMath::Abs(particle->Y())>1.0)  continue;
        
        //Variables
        Double_t pt = particle->Pt();
        Double_t protWeight = GetProtonWeight (pt);
        
        //Store Protons and Neutrons (no rapidity cut)
        if (particle->PdgCode()==-2212) { proton_ID.push_back(i); }
        if (particle->PdgCode()==-2112) { neutron_ID.push_back(i); neutron_status.push_back(0);}
        
        //Protons and Neutrons p_{T} Spectra from Pythia
        if (TMath::Abs(particle->Y())>0.5) continue;
        if (particle->PdgCode()==-2212) {hProtonsINELgtZERO  -> Fill(pt); hProtonsINELgtZERO_reshaped  -> Fill(pt,protWeight);}
        if (particle->PdgCode()==-2112) {hNeutronsINELgtZERO -> Fill(pt); hNeutronsINELgtZERO_reshaped -> Fill(pt,protWeight);}
    }
    
    
    //Save Deuteron Properties
    const Int_t nDeut_max=20;
    Double_t pt_deuteron     [nTrials][nDeut_max];
    Double_t phi_deuteron    [nTrials][nDeut_max];
    Double_t weight_deuteron [nTrials][nDeut_max];
    Int_t nDeuterons[nTrials];
    
    //Initialize Counters
    for (Int_t iTrial=0 ; iTrial<nTrials ; iTrial++) {nDeuterons[iTrial]=0;}
    
    //Loop over different trials
    for (Int_t iTrial=0 ; iTrial<nTrials ; iTrial++)  {
        
        //Reset Neutron Status
        for (Int_t in=0 ; in<(Int_t)neutron_ID.size() ; in++)  {neutron_status[in]=0;}
        
        for (Int_t ip=0 ; ip<(Int_t)proton_ID.size() ; ip++)  {
            
            AliMCParticle *proton = (AliMCParticle*) fMCEvent->GetTrack(proton_ID[ip]);
            TVector3 p_proton (proton->Px(),proton->Py(),proton->Pz());

            for (Int_t in=0 ; in<(Int_t)neutron_ID.size() ; in++)  {

                if (neutron_status[in]==1) continue;//Skip already used neutrons
                AliMCParticle *neutron = (AliMCParticle*) fMCEvent->GetTrack(neutron_ID[in]);
                TVector3 p_neutron (neutron->Px(),neutron->Py(),neutron->Pz());
                
                //Deuteron
                TVector3 p_deuteron = p_proton+p_neutron;
                Double_t deltaP = (p_proton-p_neutron).Mag();
                Double_t deutWeight = GetDeuteronWeight (p_proton.Pt(),p_neutron.Pt());

                //Rapidities
                Double_t yp = proton->Y();
                Double_t yn = neutron->Y();

                //Coalescence
                if (deltaP<p0[iTrial]) {
                    
                    neutron_status[in]=1;
                    Double_t y = GetRapidity (p_deuteron);
                    if (TMath::Abs(y)<0.5) {
                        
                        hDeuteronsINELgtZERO[iTrial] -> Fill(p_deuteron.Pt(),deutWeight);
                        
                        //Store Deuteron Properties
                        phi_deuteron[iTrial][nDeuterons[iTrial]] = TVector2::Phi_0_2pi(p_deuteron.Phi());
                        pt_deuteron[iTrial][nDeuterons[iTrial]]  = p_deuteron.Pt();
                        weight_deuteron[iTrial][nDeuterons[iTrial]] = deutWeight;
                        nDeuterons[iTrial]++;
                        
                        //Rapidity Distributions of Protons and Neutrons
                        hRapidityProtons  -> Fill (yp);
                        hRapidityNeutrons -> Fill (yn);
                    }
                    break;
                }
            }
        }
    }
    
    
    //Selection of Leading Particle
    Int_t lp = GetLeadingParticle();
    AliMCParticle *leading_particle = (AliMCParticle*) fMCEvent->GetTrack(lp);
    Double_t pt_leading  = leading_particle->Pt();
    Double_t phi_leading = TVector2::Phi_0_2pi(leading_particle->Phi());
    if (pt_leading<5.0) return;
    hNumberOfEvents -> Fill(17.5);

    //Transverse Multiplicity
    Int_t nParticles_Transv = GetTransverseMultiplicity(lp);
    hTransverseMult -> Fill(nParticles_Transv);
    
    //R_{T}
    Double_t Rt = static_cast<Double_t>(nParticles_Transv)/fAverage_Nch_Transv;
    hRtDistribution -> Fill (Rt);
    
    
    //Fill Proton Spectra in Azimuthal Region
    for (Int_t i=0 ; i<(Int_t)proton_ID.size() ; i++)  {
        
        //Get Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(proton_ID[i]);
        if (TMath::Abs(particle->Y())>0.5) continue;
        
        //Variables
        Double_t pt = particle->Pt();
        Double_t phi_particle = TVector2::Phi_0_2pi(particle->Phi());
        Double_t protWeight = GetProtonWeight (pt);

        //Fill p_{T} Spectra Protons
        if (IsParticleInTowardRegion(phi_particle,phi_leading))     hProtons_Toward->Fill(pt,protWeight);
        if (IsParticleInTransverseRegion(phi_particle,phi_leading)) hProtons_Transv->Fill(pt,protWeight);
        if (IsParticleInAwayRegion(phi_particle,phi_leading))       hProtons_Away->Fill(pt,protWeight);
    }
    
   
    //Fill Deuteron Spectra in Azimuthal Region
    for (Int_t iTrial=0 ; iTrial<nTrials ; iTrial++)  {
        for (Int_t i=0 ; i<nDeuterons[iTrial] ; i++)  {

            //Variables
            Double_t pt = pt_deuteron[iTrial][i];
            Double_t phi_particle = phi_deuteron[iTrial][i];
            Double_t deutWeight = weight_deuteron[iTrial][i];
                  
            //Fill p_{T} Spectra Protons
            if (IsParticleInTowardRegion(phi_particle,phi_leading))     hDeuterons_Toward[iTrial]->Fill(pt,deutWeight);
            if (IsParticleInTransverseRegion(phi_particle,phi_leading)) hDeuterons_Transv[iTrial]->Fill(pt,deutWeight);
            if (IsParticleInAwayRegion(phi_particle,phi_leading))       hDeuterons_Away[iTrial]  ->Fill(pt,deutWeight);
        }
    }
    
    
    //Post Output Data
    PostData(1, fOutputList);
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronCoalescence::GetEvent ()  {
    
    //Get ESD Event
    fESDevent = dynamic_cast <AliESDEvent*>(InputEvent());
    if (!fESDevent) return kFALSE;
    hNumberOfEvents -> Fill(0.5);

    //Get MC Event
    fMCEvent = MCEvent();
    if (!fMCEvent) return kFALSE;
    hNumberOfEvents -> Fill(1.5);
       
    //Get MC Event Handler
    fMCEventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!fMCEventHandler) return kFALSE;
    hNumberOfEvents -> Fill(2.5);

    //Standard Event Cuts
    if (!fESDeventCuts.AcceptEvent(fESDevent)) {
        PostData(2, fQAList);
        return kFALSE;
    }
    hNumberOfEvents -> Fill(3.5);
       
    //Reject Events with Incomplete DAQ
    if (fESDevent->IsIncompleteDAQ()) return kFALSE;
    hNumberOfEvents -> Fill(4.5);
       
    //V0 Timing Decision
    AliVVZERO *vzeroData = fESDevent->GetVZEROData();
    if (!(vzeroData->GetV0ADecision()) || !(vzeroData->GetV0CDecision())) return kFALSE;
    hNumberOfEvents -> Fill(5.5);
       
    //Pileup Rejection
    Int_t nClustersLayer0 = fESDevent->GetNumberOfITSClusters(0);
    Int_t nClustersLayer1 = fESDevent->GetNumberOfITSClusters(1);
    Int_t nTracklets      = fESDevent->GetMultiplicity()->GetNumberOfTracklets();
    if ((nClustersLayer0 + nClustersLayer1) > 65.0 + (Double_t)nTracklets*4.0) return kFALSE;
    hNumberOfEvents -> Fill(6.5);

    //Primary Vertex Tracks
    AliESDVertex *vertex_tracks = (AliESDVertex*) fESDevent->GetPrimaryVertexTracks();
    if (!vertex_tracks) return kFALSE;
    hNumberOfEvents -> Fill(7.5);
       
    //Vertex Contributors Tracks
    if ( vertex_tracks->GetNContributors() < 1 ) return kFALSE;
    hNumberOfEvents -> Fill(8.5);
       
    //Primary Vertex SPD
    AliESDVertex *vertex_SPD = (AliESDVertex*) fESDevent->GetPrimaryVertexSPD();
    if (!vertex_SPD) return kFALSE;
    hNumberOfEvents -> Fill(9.5);
       
    //Vertex Contributors SPD
    if (vertex_SPD->GetNContributors() < 1 ) return kFALSE;
    hNumberOfEvents -> Fill(10.5);
       
    //SPD Pile-up in Mult Bins
    if (fESDevent->IsPileupFromSPDInMultBins()) return kFALSE;
    hNumberOfEvents -> Fill(11.5);
       
    //Cut on Z-Vertex Resolution
    if (TMath::Abs(vertex_SPD->GetZ() - vertex_tracks->GetZ()) > 0.3) return kFALSE;
    hNumberOfEvents -> Fill(12.5);

    //Primary Vertex Selection
    if (vertex_tracks->GetZ() < -10.0) return kFALSE;
    if (vertex_tracks->GetZ() > +10.0) return kFALSE;
    hNumberOfEvents -> Fill(13.5);
              
    //Multiplicity
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
    if(!multiplicitySelection) return kFALSE;
    hNumberOfEvents -> Fill(14.5);
       
    //Selection of Multiplicity Range
    Double_t mult_percentile = multiplicitySelection->GetMultiplicityPercentile("V0M");
    if (mult_percentile <  0.0)   return kFALSE;
    if (mult_percentile >  100.0) return kFALSE;
    hNumberOfEvents -> Fill(15.5);
    
    //Selection of INEL>0 Events
    if (!IsINELgtZERO ()) return kFALSE;
    hNumberOfEvents -> Fill(16.5);
         
    return kTRUE;
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronCoalescence::IsINELgtZERO ()  {
    
    //Initialization
    Bool_t isEventINELgtZERO = (kFALSE);
       
    //Loop over Generated Particles
    for ( Int_t i=0; i<fMCEvent->GetNumberOfTracks(); i++ )  {

        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(i);
        if (!particle) continue;
        if (particle->Charge()==0)        continue;
        if (IsInjectedParticle(particle)) continue;
        if (particle->MCStatusCode()<=0)  continue;

        Bool_t isFSParticle = (kFALSE);
        if (TMath::Abs(particle->PdgCode())==211)  isFSParticle=kTRUE;
        if (TMath::Abs(particle->PdgCode())==321)  isFSParticle=kTRUE;
        if (TMath::Abs(particle->PdgCode())==2212) isFSParticle=kTRUE;
        if (TMath::Abs(particle->PdgCode())==11)   isFSParticle=kTRUE;
        if (isFSParticle && TMath::Abs(particle->Eta())<1.0) {isEventINELgtZERO = kTRUE; break;}
    }
    
    return isEventINELgtZERO;
}
//____________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDeuteronCoalescence::GetProtonWeight (Double_t pt)  {
    
    //Initialization
    Double_t w(1);
    
    Int_t ibin = hProtonWeights->FindBin(pt);
    w = hProtonWeights->GetBinContent(ibin);
    return w;
}
//____________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDeuteronCoalescence::GetDeuteronWeight (Double_t pt_prot, Double_t pt_neut)  {
    
    //Initialization
    Double_t w(1);
    
    Double_t wp = hProtonWeights->GetBinContent(hProtonWeights->FindBin(pt_prot));
    Double_t wn = hProtonWeights->GetBinContent(hProtonWeights->FindBin(pt_neut));
    w = (3.0/4.0)*wp*wn;

    return w;
}
//____________________________________________________________________________________________________________________________________________________
Int_t AliAnalysisTaskDeuteronCoalescence::GetLeadingParticle ()  {
    
    //Initialization
    Int_t ID_leading_particle(0);
    Double_t pt_max(0);
    
    //Loop over Generated Particles
    for ( Int_t i=0; i<fMCEvent->GetNumberOfTracks(); i++ )  {

        //Get Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(i);
        if (!particle) continue;
        if (particle->Charge()==0) continue;
        if (IsInjectedParticle(particle))    continue;
        if (TMath::Abs(particle->Eta())>0.8) continue;
        if (particle->MCStatusCode()<=0)  continue;
        if (particle->Pt()<0.15) continue;
        
        //Primary Particle
        Double_t x = particle->Xv();
        Double_t y = particle->Yv();
        Double_t z = particle->Zv();
        Double_t r = TMath::Sqrt(x*x+y*y+z*z);
        hRparticles -> Fill(r);
        //if (r>1e-10) continue;
        
        Bool_t isFSParticle = (kFALSE);
        if (TMath::Abs(particle->PdgCode())==211)  isFSParticle=kTRUE;
        if (TMath::Abs(particle->PdgCode())==321)  isFSParticle=kTRUE;
        if (TMath::Abs(particle->PdgCode())==2212) isFSParticle=kTRUE;
        if (TMath::Abs(particle->PdgCode())==11)   isFSParticle=kTRUE;
        if (!isFSParticle) continue;
        
        if (particle->Pt() > pt_max)  {
            
            pt_max = particle->Pt();
            ID_leading_particle = i;
        }
    }
    
    return ID_leading_particle;
}
//____________________________________________________________________________________________________________________________________________________
Int_t AliAnalysisTaskDeuteronCoalescence::GetTransverseMultiplicity (Int_t leading_particle_ID)  {
    
    //Initialization
    Int_t mult_Transverse(0);
    
    //Leading Particle
    AliMCParticle *leading_particle = (AliMCParticle*) fMCEvent->GetTrack(leading_particle_ID);
    
    //Loop over Generated Particles
    for ( Int_t i=0; i<fMCEvent->GetNumberOfTracks(); i++ )  {

        //Get Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(i);
        if (!particle) continue;
        if (particle->Charge()==0) continue;
        if (IsInjectedParticle(particle))    continue;
        if (TMath::Abs(particle->Eta())>0.8) continue;
        if (particle->MCStatusCode()<=0)     continue;
        if (particle->Pt()<0.15) continue;
               
        //Primary Particle
        Double_t x = particle->Xv();
        Double_t y = particle->Yv();
        Double_t z = particle->Zv();
        Double_t r = TMath::Sqrt(x*x+y*y+z*z);
        //if (r>1e-10) continue;
        
        Bool_t isFSParticle = (kFALSE);
        if (TMath::Abs(particle->PdgCode())==211)  isFSParticle=kTRUE;
        if (TMath::Abs(particle->PdgCode())==321)  isFSParticle=kTRUE;
        if (TMath::Abs(particle->PdgCode())==2212) isFSParticle=kTRUE;
        if (TMath::Abs(particle->PdgCode())==11)   isFSParticle=kTRUE;
        if (!isFSParticle) continue;
                   
        Double_t phi_particle = TVector2::Phi_0_2pi(particle->Phi());
        Double_t phi_leading  = TVector2::Phi_0_2pi(leading_particle->Phi());

        if (IsParticleInTransverseRegion (phi_particle,phi_leading)) mult_Transverse++;
    }
    
    return mult_Transverse;
    
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronCoalescence::IsParticleInTowardRegion (Double_t phi, Double_t phi_leading)  {
    
    //Initialization
    Bool_t isInTowardRegion = kFALSE;
              
    //DeltaPhi
    Double_t delta_phi = (180.0/TMath::Pi())*TVector2::Phi_0_2pi (phi-phi_leading);
       
    if (delta_phi>=0.0   && delta_phi<60.0)   isInTowardRegion=kTRUE;
    if (delta_phi>=300.0 && delta_phi<=360.0) isInTowardRegion=kTRUE;

    return isInTowardRegion;
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronCoalescence::IsParticleInTransverseRegion (Double_t phi, Double_t phi_leading)  {
    
    //Initialization
    Bool_t isInTransverseRegion = kFALSE;
              
    //DeltaPhi
    Double_t delta_phi = (180.0/TMath::Pi())*TVector2::Phi_0_2pi (phi-phi_leading);
       
    if (delta_phi>=60.0  && delta_phi<120.0) isInTransverseRegion=kTRUE;
    if (delta_phi>=240.0 && delta_phi<300.0) isInTransverseRegion=kTRUE;

    return isInTransverseRegion;
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronCoalescence::IsParticleInAwayRegion (Double_t phi, Double_t phi_leading)  {
    
    //Initialization
    Bool_t isInAwayRegion = kFALSE;
              
    //DeltaPhi
    Double_t delta_phi = (180.0/TMath::Pi())*TVector2::Phi_0_2pi (phi-phi_leading);
    if (delta_phi>=120.0 && delta_phi<240.0) isInAwayRegion=kTRUE;

    return isInAwayRegion;
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronCoalescence::IsInjectedParticle (AliMCParticle *particle)  {
    
    //Initialization
    Bool_t isInjected=kFALSE;
    
    if (TMath::Abs(particle->PdgCode())==1000010020) isInjected=kTRUE;//Deuteron
    if (TMath::Abs(particle->PdgCode())==1000020030) isInjected=kTRUE;//Helium-3
    if (TMath::Abs(particle->PdgCode())==1000010030) isInjected=kTRUE;//Triton
    if (TMath::Abs(particle->PdgCode())==1000020040) isInjected=kTRUE;//Helium-4
    if (TMath::Abs(particle->PdgCode())==1010010030) isInjected=kTRUE;//Hypertriton

    return isInjected;
}
//____________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDeuteronCoalescence::GetRapidity (TVector3 momentum)  {
    
    Double_t y(0);
    Double_t px = momentum.X();
    Double_t py = momentum.Y();
    Double_t pz = momentum.Z();
    Double_t m  = 1.87561294257;
    TLorentzVector P;
    P.SetXYZM (px,py,pz,m);
    y = P.Rapidity();
    
    return y;
}
//____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronCoalescence::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//____________________________________________________________________________________________________________________________________________________
