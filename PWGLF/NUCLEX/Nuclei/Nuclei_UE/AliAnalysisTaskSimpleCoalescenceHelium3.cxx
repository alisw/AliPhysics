#include "AliAnalysisTaskSimpleCoalescenceHelium3.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "TLorentzVector.h"
#include "AliMCParticle.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;
ClassImp(AliAnalysisTaskSimpleCoalescenceHelium3)

//_______________________________________________________________________________________________________________________________________
AliAnalysisTaskSimpleCoalescenceHelium3::AliAnalysisTaskSimpleCoalescenceHelium3():
AliAnalysisTaskSE(),
fAODevent(nullptr),
fMCEvent(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fProtonWeights(nullptr),
hNumberOfEvents(nullptr),
hProtons(nullptr),
hProtons_reshaped(nullptr),
hNeutrons(nullptr),
hNeutrons_reshaped(nullptr),
hHelium3{},
hDeltaP_pn(nullptr),
hDeltaP_pp(nullptr)
{}
//_______________________________________________________________________________________________________________________________________
AliAnalysisTaskSimpleCoalescenceHelium3::AliAnalysisTaskSimpleCoalescenceHelium3(const char *name):
AliAnalysisTaskSE(name),
fAODevent(nullptr),
fMCEvent(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fProtonWeights(nullptr),
hNumberOfEvents(nullptr),
hProtons(nullptr),
hProtons_reshaped(nullptr),
hNeutrons(nullptr),
hNeutrons_reshaped(nullptr),
hHelium3{},
hDeltaP_pn(nullptr),
hDeltaP_pp(nullptr)
{
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//_______________________________________________________________________________________________________________________________________
AliAnalysisTaskSimpleCoalescenceHelium3::~AliAnalysisTaskSimpleCoalescenceHelium3()  {
    
    fOutputList->Clear();
    delete fAODevent;
    delete fMCEvent;
    delete fOutputList;
    delete fQAList;
    delete fProtonWeights;
    delete hNumberOfEvents;
    delete hProtons;
    delete hProtons_reshaped;
    delete hNeutrons;
    delete hNeutrons_reshaped;
    delete hDeltaP_pn;
    delete hDeltaP_pp;

    for (Int_t i=0 ; i<100 ; i++) { delete hHelium3[i]; }

}
//_______________________________________________________________________________________________________________________________________
void AliAnalysisTaskSimpleCoalescenceHelium3::UserCreateOutputObjects()  {
    
    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner();
    fQAList     -> SetOwner();
    
    //Event Counter
    hNumberOfEvents = new TH1D ("hNumberOfEvents","",20,0,20);
    hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    
    //p_{T} Intervals
    Double_t pt_proton[] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,
        1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,10.0};
    Double_t pt_helium3[] = {1.0,1.5,2.0,2.5,3.0,4.0,5.0};
    const Int_t nPtProton  = sizeof(pt_proton)/sizeof(Double_t)-1;
    const Int_t nPtHelium3 = sizeof(pt_helium3)/sizeof(Double_t)-1;
    
    //p_{T} Spectra in Events INEL>0
    hProtons            = new TH1D ("hProtons","",nPtProton,pt_proton);
    hProtons_reshaped   = new TH1D ("hProtons_reshaped","",nPtProton,pt_proton);
    hNeutrons           = new TH1D ("hNeutrons","",nPtProton,pt_proton);
    hNeutrons_reshaped  = new TH1D ("hNeutrons_reshaped","",nPtProton,pt_proton);
    hProtons            -> Sumw2();
    hNeutrons           -> Sumw2();
    hProtons_reshaped   -> Sumw2();
    hNeutrons_reshaped  -> Sumw2();
    fOutputList -> Add(hProtons);
    fOutputList -> Add(hNeutrons);
    fOutputList -> Add(hProtons_reshaped);
    fOutputList -> Add(hNeutrons_reshaped);
    
    //p_{T} Spectra Helium3
    for (Int_t i=0 ; i<100 ; i++)  {
        
        hHelium3[i] = new TH1D (Form("hHelium3[%d]",i),"",nPtHelium3,pt_helium3);
        hHelium3[i] -> Sumw2();
        fOutputList -> Add(hHelium3[i]);
    }
   
    //DeltaP
    hDeltaP_pn = new TH1D ("hDeltaP_pn","",1000,0,1);
    hDeltaP_pp = new TH1D ("hDeltaP_pp","",1000,0,1);
    hDeltaP_pn -> Sumw2();
    hDeltaP_pp -> Sumw2();
    fOutputList -> Add(hDeltaP_pn);
    fOutputList -> Add(hDeltaP_pp);
    
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//_______________________________________________________________________________________________________________________________________
void AliAnalysisTaskSimpleCoalescenceHelium3::UserExec(Option_t *)  {
    
    //Coalescence Momentum
    /*
    Double_t p0[] = {
        0.200,0.205,0.210,0.215,0.220,0.225,0.230,0.235,0.240,0.245,0.250,0.255,0.260,0.265,0.270,0.275,0.280,0.285,0.290,0.295,
        0.300,0.305,0.310,0.315,0.320,0.325,0.330,0.335,0.340,0.345,0.350,0.355,0.360,0.365,0.370,0.375,0.380,0.385,0.390,0.395,0.400
    };*/
    
    //Coalescence Momentum: pp Pairs
    Double_t p0_pp[] = {
        0.22,0.23,0.24,0.25,0.26,0.27,0.28,
        0.22,0.23,0.24,0.25,0.26,0.27,0.28,
        0.22,0.23,0.24,0.25,0.26,0.27,0.28,
        0.22,0.23,0.24,0.25,0.26,0.27,0.28,
        0.22,0.23,0.24,0.25,0.26,0.27,0.28,
        0.22,0.23,0.24,0.25,0.26,0.27,0.28,
        0.22,0.23,0.24,0.25,0.26,0.27,0.28
    };
    
    //Coalescence Momentum: pn Pairs
    Double_t p0_pn[] = {
        0.22,0.22,0.22,0.22,0.22,0.22,0.22,
        0.23,0.23,0.23,0.23,0.23,0.23,0.23,
        0.24,0.24,0.24,0.24,0.24,0.24,0.24,
        0.25,0.25,0.25,0.25,0.25,0.25,0.25,
        0.26,0.26,0.26,0.26,0.26,0.26,0.26,
        0.27,0.27,0.27,0.27,0.27,0.27,0.27,
        0.28,0.28,0.28,0.28,0.28,0.28,0.28
    };
    
    const Int_t nTrials = sizeof(p0_pp)/sizeof(Double_t);
    
    //Particle Masses [GeV/c^2]
    Double_t mp = 0.93827208816;//Proton
    Double_t mn = 0.93956542052;//Neutron
    Double_t md = 1.87561294257;//Deuteron
    Double_t mh = 2.80839160743;//Helium3 (physics.nist.gov/cgi-bin/cuu/Value?mhc2mev)

    //Get Input Event (INEL>0 Selection)
    if (!GetEvent()) return;
   
    //Protons and Neutrons IDs
    vector<Int_t> proton_ID;
    vector<Int_t> neutron_ID;
    vector<Int_t> proton_status;
    vector<Int_t> neutron_status;

    //Loop over Generated Particles
    for ( Int_t i=0; i<fMCEvent->GetNumberOfTracks(); i++ )  {

        //Get Primary Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(i);
        if (!particle) continue;
        if (!particle->IsPhysicalPrimary()) continue;
        if (IsInjectedParticle(particle))   continue;
        if (TMath::Abs(particle->Y())>1.5)  continue;
        
        //Variables
        Double_t pt = particle->Pt();
        Double_t wp = GetProtonWeight (pt);
        
        //Store Protons and Neutrons
        if (particle->PdgCode()==-2212) { proton_ID.push_back(i);  proton_status.push_back(0);}
        if (particle->PdgCode()==-2112) { neutron_ID.push_back(i); neutron_status.push_back(0);}
        
        //Protons and Neutrons p_{T} Spectra
        if (TMath::Abs(particle->Y())>0.5) continue;
        if (particle->PdgCode()==-2212) {hProtons  -> Fill(pt); hProtons_reshaped  -> Fill(pt,wp);}
        if (particle->PdgCode()==-2112) {hNeutrons -> Fill(pt); hNeutrons_reshaped -> Fill(pt,wp);}
    }
    
    //Skip Events with not Enough Nucleons
    if ((Int_t)neutron_ID.size()<1) return;
    if ((Int_t)proton_ID.size()<2)  return;

    
    //Helium3 Counters
    Int_t nHelium3[nTrials];

    //Initialize Helium3 Counters
    for (Int_t iTrial=0 ; iTrial<nTrials ; iTrial++) {nHelium3[iTrial]=0;}
    
    
    //Coalescence
    for (Int_t iTrial=0 ; iTrial<nTrials ; iTrial++)  {
        
        //Reset Neutron Status
        for (Int_t in=0 ; in<(Int_t)neutron_ID.size() ; in++)  {neutron_status[in]=0;}
        
        for (Int_t ip1=0 ; ip1<(Int_t)proton_ID.size() ; ip1++)  {
            for (Int_t ip2=ip1+1 ; ip2<(Int_t)proton_ID.size() ; ip2++)  {
                
                //Proton 4-Momenta
                AliMCParticle *proton1 = (AliMCParticle*) fMCEvent->GetTrack(proton_ID[ip1]);
                AliMCParticle *proton2 = (AliMCParticle*) fMCEvent->GetTrack(proton_ID[ip2]);
                TLorentzVector p_proton1,p_proton2;
                p_proton1.SetXYZM(proton1->Px(),proton1->Py(),proton1->Pz(),mp);
                p_proton2.SetXYZM(proton2->Px(),proton2->Py(),proton2->Pz(),mp);

                for (Int_t in=0 ; in<(Int_t)neutron_ID.size() ; in++)  {
                    
                    //Neutron 4-Momentum
                    AliMCParticle *neutron = (AliMCParticle*) fMCEvent->GetTrack(neutron_ID[in]);
                    TLorentzVector p_neutron;
                    p_neutron.SetXYZM(neutron->Px(),neutron->Py(),neutron->Pz(),mn);
                    
                    //Helium3 4-Momentum
                    TLorentzVector p_helium3;
                    Double_t px_helium3 = proton1->Px()+proton2->Px()+neutron->Px();
                    Double_t py_helium3 = proton1->Py()+proton2->Py()+neutron->Py();
                    Double_t pz_helium3 = proton1->Pz()+proton2->Pz()+neutron->Pz();
                    p_helium3.SetXYZM(px_helium3,py_helium3,pz_helium3,mh);
                    Double_t beta_x = p_helium3.Px()/p_helium3.E();
                    Double_t beta_y = p_helium3.Py()/p_helium3.E();
                    Double_t beta_z = p_helium3.Pz()/p_helium3.E();
                    TVector3 beta (beta_x,beta_y,beta_z);
                    
                    //Lorentz Transformation
                    TLorentzVector p_proton1_prime = LorentzTransform (p_proton1,beta);
                    TLorentzVector p_proton2_prime = LorentzTransform (p_proton2,beta);
                    TLorentzVector p_neutron_prime = LorentzTransform (p_neutron,beta);
                    
                    //Variables
                    Double_t deltaP1 = (p_proton1_prime-p_proton2_prime).P();
                    Double_t deltaP2 = (p_proton1_prime-p_neutron_prime).P();
                    Double_t deltaP3 = (p_proton2_prime-p_neutron_prime).P();

                    //DeltaP
                    if (iTrial==0) {
                        
                        hDeltaP_pp->Fill(deltaP1);
                        hDeltaP_pn->Fill(deltaP2);
                    }

                    //Weight
                    Double_t wh = GetHelium3Weight (p_proton1.Pt(),p_proton2.Pt(),p_neutron.Pt());
                    
                    //Skip already used Neutrons
                    if (neutron_status[in]==1) continue;
                    
                    //Coalescence Condition
                    if (ThreeBodyCoalescence(deltaP1,deltaP2,deltaP3,p0_pp[iTrial],p0_pn[iTrial]))  {
                        
                        neutron_status[in]=1;
                        Double_t y = p_helium3.Rapidity();
                        if (TMath::Abs(y)<0.5) {
                            
                            hHelium3[iTrial] -> Fill(p_helium3.Pt(),wh);
                            nHelium3[iTrial]++;
                        }
                        break;
                    }
                }
            }
        }
        
    }
    
    //Post Output Data
    PostData(1, fOutputList);
}
//_______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskSimpleCoalescenceHelium3::GetEvent ()  {
    
    //Get AOD Event
    fAODevent = dynamic_cast <AliAODEvent*>(InputEvent());
    if (!fAODevent) return kFALSE;
    hNumberOfEvents -> Fill(0.5);

    //Get MC Event
    fMCEvent = MCEvent();
    if (!fMCEvent) return kFALSE;
    hNumberOfEvents -> Fill(1.5);

    //Selection of INEL>0 Events
    if (!IsINELgtZERO ()) return kFALSE;
    hNumberOfEvents -> Fill(2.5);

    return kTRUE;
}
//_______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskSimpleCoalescenceHelium3::IsINELgtZERO ()  {
    
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
//_______________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskSimpleCoalescenceHelium3::GetProtonWeight (Double_t pt)  {
    
    //Initialization
    Double_t wp(1);
    wp = fProtonWeights->Eval(pt);
    
    return wp;
}
//_______________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskSimpleCoalescenceHelium3::GetHelium3Weight (Double_t pt_prot1, Double_t pt_prot2, Double_t pt_neut)  {
    
    //Initialization
    Double_t wh(1);
    
    Double_t wp1 = fProtonWeights->Eval(pt_prot1);
    Double_t wp2 = fProtonWeights->Eval(pt_prot2);
    Double_t wn  = fProtonWeights->Eval(pt_neut);
    
    Double_t spin_factor = 1.0/4.0;
    wh = spin_factor*wp1*wp2*wn;

    return wh;
}
//_______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskSimpleCoalescenceHelium3::ThreeBodyCoalescence (Double_t deltaP1, Double_t deltaP2, Double_t deltaP3, Double_t p0_pp, Double_t p0_pn)  {
    
    //Initialization
    Bool_t isBoundStateFormed=(kFALSE);
    
    if ((deltaP1<p0_pp)&&(deltaP2<p0_pn)&&(deltaP3<p0_pn)) isBoundStateFormed=kTRUE;
    //if ((deltaP1<p0)&&(deltaP2<p0)&&(deltaP3<p0)) isBoundStateFormed=kTRUE;
   
    return isBoundStateFormed;
}
//_______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskSimpleCoalescenceHelium3::IsInjectedParticle (AliMCParticle *particle)  {
    
    //Initialization
    Bool_t isInjected=kFALSE;
    
    if (TMath::Abs(particle->PdgCode())==1000010020) isInjected=kTRUE;//Deuteron
    if (TMath::Abs(particle->PdgCode())==1000020030) isInjected=kTRUE;//Helium-3
    if (TMath::Abs(particle->PdgCode())==1000010030) isInjected=kTRUE;//Triton
    if (TMath::Abs(particle->PdgCode())==1000020040) isInjected=kTRUE;//Helium-4
    if (TMath::Abs(particle->PdgCode())==1010010030) isInjected=kTRUE;//Hypertriton

    return isInjected;
}
//_______________________________________________________________________________________________________________________________________
TLorentzVector AliAnalysisTaskSimpleCoalescenceHelium3::LorentzTransform (TLorentzVector R, TVector3 beta_vect)  {
    
    //Inizialization
    TLorentzVector R_prime (0.0, 0.0, 0.0, 0.0);
    
    //Beta Components
    Double_t Bx = beta_vect.X();
    Double_t By = beta_vect.Y();
    Double_t Bz = beta_vect.Z();
    
    //Beta & Gamma
    Double_t beta  = TMath::Sqrt(Bx*Bx + By*By + Bz*Bz);
    if (beta>=1.0) { return R_prime; }
    Double_t gamma = 1.0/TMath::Sqrt(1.0-(beta*beta));
    
    //Coordinates in the Lab System
    Double_t t = R.T();
    Double_t x = R.X();
    Double_t y = R.Y();
    Double_t z = R.Z();
    
    //Coordinates in the Deuteron Frame
    Double_t t_prime = gamma*t - gamma*Bx*x - gamma*By*y - gamma*Bz*z;
    Double_t x_prime = -gamma*Bx*t + (1.0+(gamma-1.0)*Bx*Bx/(beta*beta))*x + (gamma-1.0)*(Bx*By/(beta*beta))*y + (gamma-1.0)*(Bx*Bz/(beta*beta))*z;
    Double_t y_prime = -gamma*By*t + (gamma-1.0)*(Bx*By/(beta*beta))*x + (1.0+(gamma-1.0)*By*By/(beta*beta))*y + (gamma-1.0)*(By*Bz/(beta*beta))*z;
    Double_t z_prime = -gamma*Bz*t + (gamma-1.0)*(Bx*Bz/(beta*beta))*x + (gamma-1.0)*(By*Bz/(beta*beta))*y + (1.0+(gamma-1.0)*Bz*Bz/(beta*beta))*z;

    //Set Coordinates
    R_prime.SetXYZT(x_prime,y_prime,z_prime,t_prime);
    
    return R_prime;
}
//_______________________________________________________________________________________________________________________________________
void AliAnalysisTaskSimpleCoalescenceHelium3::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_______________________________________________________________________________________________________________________________________
