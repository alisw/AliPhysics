#include "AliAnalysisTaskSimpleCoalescenceDeuteronInJets.h"
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
#include "TH1I.h"
#include "TH2D.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;
ClassImp(AliAnalysisTaskSimpleCoalescenceDeuteronInJets)

//_______________________________________________________________________________________________________________________________________
AliAnalysisTaskSimpleCoalescenceDeuteronInJets::AliAnalysisTaskSimpleCoalescenceDeuteronInJets():
AliAnalysisTaskSE(),
fAODevent(nullptr),
fMCEvent(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fProtonWeights(nullptr),
hNumberOfEvents(nullptr),
hProtons{},
hNeutrons{},
hProtons_pythia{},
hNeutrons_pythia{},
hDeuterons{},
hDeuterons_pythia{},
hDeuterons_ptoverA{},
hDeuterons_ptoverA_pythia{},
hDeltaP_Lab{},
hDeltaP_CM{},
hPtProtonsCoal{},
hPtNeutronsCoal{},
hTheta_Lab{},
hTheta_CM{},
hNumberOfParticlesInJet(nullptr),
hNumberOfProtonsInJet(nullptr),
hNumberOfAntiProtonsInJet(nullptr),
hEventsVsMultiplicity(nullptr),
fMaximumPt(5.0),
fJetRadius(0.5)
{}
//_______________________________________________________________________________________________________________________________________
AliAnalysisTaskSimpleCoalescenceDeuteronInJets::AliAnalysisTaskSimpleCoalescenceDeuteronInJets(const char *name):
AliAnalysisTaskSE(name),
fAODevent(nullptr),
fMCEvent(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fProtonWeights(nullptr),
hNumberOfEvents(nullptr),
hProtons{},
hNeutrons{},
hProtons_pythia{},
hNeutrons_pythia{},
hDeuterons{},
hDeuterons_pythia{},
hDeuterons_ptoverA{},
hDeuterons_ptoverA_pythia{},
hDeltaP_Lab{},
hDeltaP_CM{},
hPtProtonsCoal{},
hPtNeutronsCoal{},
hTheta_Lab{},
hTheta_CM{},
hNumberOfParticlesInJet(nullptr),
hNumberOfProtonsInJet(nullptr),
hNumberOfAntiProtonsInJet(nullptr),
hEventsVsMultiplicity(nullptr),
fMaximumPt(5.0),
fJetRadius(0.5)
{
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//_______________________________________________________________________________________________________________________________________
AliAnalysisTaskSimpleCoalescenceDeuteronInJets::~AliAnalysisTaskSimpleCoalescenceDeuteronInJets()  {
    
    fOutputList->Clear();
    delete fAODevent;
    delete fMCEvent;
    delete fOutputList;
    delete fQAList;
    delete fProtonWeights;
    delete hNumberOfEvents;
    delete hNumberOfParticlesInJet;
    delete hNumberOfProtonsInJet;
    delete hNumberOfAntiProtonsInJet;
    delete hEventsVsMultiplicity;
    
    for (Int_t i=0 ; i<6 ; i++) {
        
        delete hProtons[i];
        delete hNeutrons[i];
        delete hProtons_pythia[i];
        delete hNeutrons_pythia[i];
        delete hDeuterons[i];
        delete hDeuterons_pythia[i];
        delete hDeuterons_ptoverA[i];
        delete hDeuterons_ptoverA_pythia[i];
       
        delete hDeltaP_Lab[i];
        delete hDeltaP_CM[i];
        delete hPtProtonsCoal[i];
        delete hPtNeutronsCoal[i];
        delete hTheta_Lab[i];
        delete hTheta_CM[i];
    }
}
//_______________________________________________________________________________________________________________________________________
void AliAnalysisTaskSimpleCoalescenceDeuteronInJets::UserCreateOutputObjects()  {
    
    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner();
    fQAList     -> SetOwner();
    
    //Event Counter
    hNumberOfEvents = new TH1D ("hNumberOfEvents","",20,0,20);
    hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    
        
    for (Int_t i=0 ; i<6 ; i++)  {
        
        //p_{T} Spectra of Nucleons: Re-weighted
        hProtons[i]  = new TH1D (Form("hProtons[%d]",i),"",1000,0,10);
        hNeutrons[i] = new TH1D (Form("hNeutrons[%d]",i),"",1000,0,10);
        hProtons[i]  -> Sumw2();
        hNeutrons[i] -> Sumw2();
        fOutputList  -> Add(hProtons[i]);
        fOutputList  -> Add(hNeutrons[i]);
        
        //p_{T} Spectra of Nucleons: from PYTHIA
        hProtons_pythia[i]  = new TH1D (Form("hProtons_pythia[%d]",i),"",1000,0,10);
        hNeutrons_pythia[i] = new TH1D (Form("hNeutrons_pythia[%d]",i),"",1000,0,10);
        hProtons_pythia[i]  -> Sumw2();
        hNeutrons_pythia[i] -> Sumw2();
        fOutputList         -> Add(hProtons_pythia[i]);
        fOutputList         -> Add(hNeutrons_pythia[i]);

        //p_{T} Spectra Deuterons
        hDeuterons[i]                = new TH1D (Form("hDeuterons[%d]",i),"",1000,0,10);
        hDeuterons_pythia[i]         = new TH1D (Form("hDeuterons_pythia[%d]",i),"",1000,0,10);
        hDeuterons_ptoverA[i]        = new TH1D (Form("hDeuterons_ptoverA[%d]",i),"",1000,0,10);
        hDeuterons_ptoverA_pythia[i] = new TH1D (Form("hDeuterons_ptoverA_pythia[%d]",i),"",1000,0,10);
        hDeuterons[i]                -> Sumw2();
        hDeuterons_pythia[i]         -> Sumw2();
        hDeuterons_ptoverA[i]        -> Sumw2();
        hDeuterons_ptoverA_pythia[i] -> Sumw2();
        fOutputList -> Add(hDeuterons[i]);
        fOutputList -> Add(hDeuterons_pythia[i]);
        fOutputList -> Add(hDeuterons_ptoverA[i]);
        fOutputList -> Add(hDeuterons_ptoverA_pythia[i]);
        
        //DeltaP Lab
        hDeltaP_Lab[i]  = new TH1D (Form("hDeltaP_Lab[%d]",i),"",1000,0,1);
        hDeltaP_Lab[i]  -> Sumw2();
        fOutputList -> Add(hDeltaP_Lab[i]);
        
        //DeltaP CM
        hDeltaP_CM[i]  = new TH1D (Form("hDeltaP_CM[%d]",i),"",1000,0,1);
        hDeltaP_CM[i]  -> Sumw2();
        fOutputList -> Add(hDeltaP_CM[i]);
        
        //p_{T} Spectrum Coalescing Protons
        hPtProtonsCoal[i]  = new TH1D (Form("hPtProtonsCoal[%d]",i),"",500,0,5);
        hPtProtonsCoal[i]  -> Sumw2();
        fOutputList -> Add(hPtProtonsCoal[i]);

        //p_{T} Spectrum Coalescing Neutrons
        hPtNeutronsCoal[i]  = new TH1D (Form("hPtNeutronsCoal[%d]",i),"",500,0,5);
        hPtNeutronsCoal[i]  -> Sumw2();
        fOutputList -> Add(hPtNeutronsCoal[i]);

        //Theta Lab
        hTheta_Lab[i]  = new TH1D (Form("hTheta_Lab[%d]",i),"",3600,0,360);
        hTheta_Lab[i]  -> Sumw2();
        fOutputList -> Add(hTheta_Lab[i]);

        //Theta CM
        hTheta_CM[i]  = new TH1D (Form("hTheta_CM[%d]",i),"",3600,0,360);
        hTheta_CM[i]  -> Sumw2();
        fOutputList -> Add(hTheta_CM[i]);
    }
    
    //General Histograms
    hNumberOfParticlesInJet = new TH1I ("hNumberOfParticlesInJet","",200,0,200);
    hNumberOfParticlesInJet -> Sumw2();
    fOutputList -> Add(hNumberOfParticlesInJet);
    
    hNumberOfProtonsInJet = new TH1I ("hNumberOfProtonsInJet","",200,0,200);
    hNumberOfProtonsInJet -> Sumw2();
    fOutputList -> Add(hNumberOfProtonsInJet);

    hNumberOfAntiProtonsInJet = new TH1I ("hNumberOfAntiProtonsInJet","",200,0,200);
    hNumberOfAntiProtonsInJet -> Sumw2();
    fOutputList -> Add(hNumberOfAntiProtonsInJet);

    
    //Histogram to select Multiplicity Bin
    Double_t multiplicity_interval[]={0,5,10,20,30,200};
    Int_t nMultBins  = sizeof(multiplicity_interval)/sizeof(Double_t)-1;
    hEventsVsMultiplicity = new TH1D ("hEventsVsMultiplicity","",nMultBins,multiplicity_interval);
    fOutputList -> Add(hEventsVsMultiplicity);

  
    //Post Data
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//_______________________________________________________________________________________________________________________________________
void AliAnalysisTaskSimpleCoalescenceDeuteronInJets::UserExec(Option_t *)  {
    
    //Coalescence Momentum [GeV/c]
    Double_t p0 = 0.2304;
    
    //Particle Masses [GeV/c^2]
    const Double_t mp = 0.93827208816;//Proton
    const Double_t mn = 0.93956542052;//Neutron
    const Double_t md = 1.87561294257;//Deuteron
    
    //Get Input Event (INEL>0 Selection)
    if (!GetEvent()) return;
    
    
    //************************************************************************************************
    //List of Particle Indices
    vector<Int_t> particle_ID;
    
    //Leading Particle ID
    Int_t leading_ID;

    //Flags
    Bool_t containsProton(kFALSE);
    Bool_t containsNeutron(kFALSE);
    
    //Maximum p_{T}
    Double_t pt_max(0);
    
    //Loop over Generated Particles
    for (Int_t i=0; i<fMCEvent->GetNumberOfTracks(); i++)  {
        
        //Get MC Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(i);
        if (!particle) continue;
        if (IsInjectedParticle(particle))        continue;
        if (particle->MCStatusCode()<=0)         continue;
        if (particle->IsSecondaryFromMaterial()) continue;
        if (TMath::Abs(particle->Eta())>0.8)     continue;
        
        //PDG Selection
        Int_t pdg = TMath::Abs(particle->PdgCode());
        if ((pdg!=11) && (pdg!=211) && (pdg!=321) && (pdg!=2212) && (pdg!=2112)) continue;
        
        //Nucleon Selection
        if (pdg==2212) containsProton=kTRUE;
        if (pdg==2112) containsNeutron=kTRUE;

        //Select Leading Particle
        if (particle->Charge()!=0 && particle->Pt()>pt_max) { leading_ID = i; pt_max = particle->Pt(); }
        
        //Store Array Element
        particle_ID.push_back(i);
    }
   
    //Selection of Events with Enough Nucleons
    if ((Int_t)particle_ID.size()<2) return;
    if (!containsProton)   return;
    if (pt_max<fMaximumPt) return;
    
    //4-Momentum of Leading Particle
    TLorentzVector P_leading (0,0,0,0);
    AliMCParticle *leading_particle = (AliMCParticle*) fMCEvent->GetTrack(leading_ID);
    P_leading.SetPxPyPzE(leading_particle->Px(),leading_particle->Py(),leading_particle->Pz(),leading_particle->E());
    Double_t pt_leading = P_leading.Pt();
    
    //Array of Particles inside Jet
    vector<Int_t> jet_particle_ID;
    jet_particle_ID.push_back(leading_ID);
    
    //Labels
    Int_t exit(0);
    Int_t nPartAssociated(0);
    
    //Jet Finder
    do {
        
        //Initialization
        Double_t distance_jet_min(1e+08);
        Double_t distance_bkg_min(1e+08);
        Int_t label_jet_particle(0);
        Int_t i_jet_particle(0);
        
        for (Int_t i=0 ; i<(Int_t)particle_ID.size() ; i++)  {
            
            //Skip Leading Particle & Elements already associated to the Jet
            if (particle_ID[i]==leading_ID) continue;
            if (particle_ID[i]==-1)         continue;
            
            //Get Particle 4-Momentum
            TLorentzVector P_particle(0,0,0,0);
            AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(particle_ID[i]);
            P_particle.SetPxPyPzE(particle->Px(),particle->Py(),particle->Pz(),particle->E());
                             
            //Variables
            Double_t one_over_pt2_part = 1.0/(P_particle.Pt()*P_particle.Pt());
            Double_t one_over_pt2_lead = 1.0/(P_leading.Pt()*P_leading.Pt());
            Double_t deltaY   = P_particle.Rapidity()-P_leading.Rapidity();
            Double_t deltaPhi = TVector2::Phi_0_2pi(P_particle.Phi()-P_leading.Phi());
            Double_t min      = Minimum (one_over_pt2_part,one_over_pt2_lead);
            Double_t Delta2   = deltaY*deltaY + deltaPhi*deltaPhi;
            
            //Distances
            Double_t distance_jet = min*Delta2/(fJetRadius*fJetRadius);
            Double_t distance_bkg = one_over_pt2_part;
            
            //Find Minimum Distance Jet
            if (distance_jet<distance_jet_min)  {
                distance_jet_min=distance_jet;
                label_jet_particle=particle_ID[i];
                i_jet_particle = i;
            }
            
            //Find Minimum Distance Bkg
            if (distance_bkg<distance_bkg_min)  {
                distance_bkg_min=distance_bkg;
            }
        }
        
        if (distance_jet_min<=distance_bkg_min)  {
            
            //Add Particle to Jet
            jet_particle_ID.push_back(label_jet_particle);

            //Update 4-Momentum of Leading Particle
            TLorentzVector P_i(0,0,0,0);
            AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(label_jet_particle);
            Double_t px_i = particle->Px();
            Double_t py_i = particle->Py();
            Double_t pz_i = particle->Pz();
            Double_t E_i  = particle->E();
            P_i.SetPxPyPzE(px_i,py_i,pz_i,E_i);
            P_leading = P_leading + P_i;
            
            //Remove Element
            particle_ID[i_jet_particle] = -1;
            nPartAssociated++;
        }
            
        if (nPartAssociated>=((Int_t)particle_ID.size()-1)) exit=1;
        if (distance_jet_min>distance_bkg_min) exit=2;
            
    } while (exit==0);
    
    //************************************************************************************************
  
    //Charged-Particle Counter
    Int_t nParticlesJet(0);
    Int_t nProtons(0);
    Int_t nAntiProtons(0);

    //Loop over Particles inside Jets
    for (Int_t i=0 ; i<(Int_t)jet_particle_ID.size() ; i++)  {
        
        //Get Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(jet_particle_ID[i]);
        if (particle->Charge()!=0) nParticlesJet++;
        if (particle->PdgCode()==2212) nProtons++;
        if (particle->PdgCode()==-2212) nAntiProtons++;
    }
    
    //Fill Multiplicity
    hNumberOfParticlesInJet -> Fill(nParticlesJet);
    hNumberOfProtonsInJet   -> Fill(nProtons);
    hNumberOfAntiProtonsInJet   -> Fill(nAntiProtons);

    //Selection of Multiplicity Bin
    Int_t iMult = hEventsVsMultiplicity->FindBin(nParticlesJet);
    hEventsVsMultiplicity->Fill(nParticlesJet);
    
    //Protons and Neutrons IDs
    vector<Int_t> proton_ID;
    vector<Int_t> neutron_ID;
    vector<Int_t> neutron_status;
    
    for (Int_t i=0 ; i<(Int_t)jet_particle_ID.size() ; i++)  {
        
        //Get Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(jet_particle_ID[i]);
        if (!particle->IsPhysicalPrimary()) continue;

        //Variables
        Double_t pt = particle->Pt();
        Double_t wp = GetProtonWeight (pt);
        Int_t pdg = TMath::Abs(particle->PdgCode());
        if (pt>5.0) continue;
        
        //Store Proton ID
        if (pdg==2212) {
            proton_ID.push_back(jet_particle_ID[i]);
            if (TMath::Abs(particle->Y())<0.5) {
                hProtons[0]            -> Fill(particle->Pt(),wp);
                hProtons_pythia[0]     -> Fill(particle->Pt());
                hProtons[iMult]        -> Fill(particle->Pt(),wp);
                hProtons_pythia[iMult] -> Fill(particle->Pt());
            }
        }
        
        //Store Neutron ID
        if (pdg==2112) {
            neutron_ID.push_back(jet_particle_ID[i]);
            neutron_status.push_back(0);
            if (TMath::Abs(particle->Y())<0.5) {
                hNeutrons[0]            -> Fill(particle->Pt(),wp);
                hNeutrons_pythia[0]     -> Fill(particle->Pt());
                hNeutrons[iMult]        -> Fill(particle->Pt(),wp);
                hNeutrons_pythia[iMult] -> Fill(particle->Pt());
            }
        }
    }
        
    //Skip Events with not Enough Nucleons
    if ((Int_t)proton_ID.size()<1)  return;
    if ((Int_t)neutron_ID.size()<1) return;


    //Deuteron Counter
    Int_t nDeuterons(0);

    //Coalescence
    for (Int_t ip=0 ; ip<(Int_t)proton_ID.size() ; ip++)  {
            
        //Proton 4-Momentum
        AliMCParticle *proton = (AliMCParticle*) fMCEvent->GetTrack(proton_ID[ip]);
        TLorentzVector p_proton;
        p_proton.SetXYZM(proton->Px(),proton->Py(),proton->Pz(),mp);
        Int_t pdg_proton = proton->PdgCode();
        
        for (Int_t in=0 ; in<(Int_t)neutron_ID.size() ; in++)  {

            //Neutron 4-Momentum
            AliMCParticle *neutron = (AliMCParticle*) fMCEvent->GetTrack(neutron_ID[in]);
            TLorentzVector p_neutron;
            p_neutron.SetXYZM(neutron->Px(),neutron->Py(),neutron->Pz(),mn);
            Int_t pdg_neutron = neutron->PdgCode();

            //Deuteron 4-Momentum
            TLorentzVector p_deuteron;
            p_deuteron.SetXYZM(proton->Px()+neutron->Px(),proton->Py()+neutron->Py(),proton->Pz()+neutron->Pz(),md);
            Double_t beta_x = p_deuteron.Px()/p_deuteron.E();
            Double_t beta_y = p_deuteron.Py()/p_deuteron.E();
            Double_t beta_z = p_deuteron.Pz()/p_deuteron.E();
            TVector3 beta (beta_x,beta_y,beta_z);
            
            //Beta
            Double_t beta_d = beta.Mag();
            if (beta_d>=1) continue;
            
            //Lorentz Transformation
            TLorentzVector p_proton_prime  = LorentzTransform (p_proton,beta);
            TLorentzVector p_neutron_prime = LorentzTransform (p_neutron,beta);
                
            //Variables
            Double_t deltaP     = (p_proton_prime-p_neutron_prime).P();
            Double_t deltaP_Lab = (p_proton-p_neutron).P();
            Double_t deutWeight = GetDeuteronWeight (p_proton.Pt(),p_neutron.Pt());
            Double_t angle_Lab  = (180.0/TMath::Pi())* p_proton.Angle(p_neutron.Vect());
            Double_t angle_CM   = (180.0/TMath::Pi())* p_proton_prime.Angle(p_neutron_prime.Vect());


            //Skip already used Neutrons
            if (neutron_status[in]==1) continue;
                
            //Skip Particle-Antiparticle Combinations
            if (pdg_proton*pdg_neutron<0) continue;
                
            //DeltaP
            hDeltaP_CM[0]      -> Fill(deltaP);
            hDeltaP_CM[iMult]  -> Fill(deltaP);
            hDeltaP_Lab[0]     -> Fill(deltaP_Lab);
            hDeltaP_Lab[iMult] -> Fill(deltaP_Lab);

            //Theta
            hTheta_Lab[0]      -> Fill(angle_Lab);
            hTheta_Lab[iMult]  -> Fill(angle_Lab);
            hTheta_CM[0]       -> Fill(angle_CM);
            hTheta_CM[iMult]   -> Fill(angle_CM);

            //Coalescence Condition
            if (TwoBodyCoalescence(deltaP,p0))  {
                    
                neutron_status[in]=1;
                Double_t y = p_deuteron.Rapidity();
                if (TMath::Abs(y)<0.5) {
                   
                    //Deuteron Spectra (with Weights)
                    hDeuterons[0]             -> Fill(p_deuteron.Pt(),deutWeight);
                    hDeuterons_ptoverA[0]     -> Fill(p_deuteron.Pt()/2.0,deutWeight);
                    hDeuterons[iMult]         -> Fill(p_deuteron.Pt(),deutWeight);
                    hDeuterons_ptoverA[iMult] -> Fill(p_deuteron.Pt()/2.0,deutWeight);
                    
                    //Deuteron Spectra (no Weights)
                    hDeuterons_pythia[0]             -> Fill(p_deuteron.Pt(),(3.0/4.0));
                    hDeuterons_ptoverA_pythia[0]     -> Fill(p_deuteron.Pt()/2.0,(3.0/4.0));
                    hDeuterons_pythia[iMult]         -> Fill(p_deuteron.Pt(),(3.0/4.0));
                    hDeuterons_ptoverA_pythia[iMult] -> Fill(p_deuteron.Pt()/2.0,(3.0/4.0));

                    //Proton and Neutron p_{T} Spectra
                    hPtProtonsCoal[0]      -> Fill(p_proton.Pt());
                    hPtProtonsCoal[iMult]  -> Fill(p_proton.Pt());
                    hPtNeutronsCoal[0]     -> Fill(p_neutron.Pt());
                    hPtNeutronsCoal[iMult] -> Fill(p_neutron.Pt());
                    nDeuterons++;
                }
                break;
            }
        }
    }
    
    
    //Post Output Data
    PostData(1, fOutputList);
}
//_______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskSimpleCoalescenceDeuteronInJets::GetEvent ()  {
    
    //Get AOD Event
    fAODevent = dynamic_cast <AliAODEvent*>(InputEvent());
    if (!fAODevent) return kFALSE;
    hNumberOfEvents -> Fill(0.5);

    //Get MC Event
    fMCEvent = MCEvent();
    if (!fMCEvent) return kFALSE;
    hNumberOfEvents -> Fill(1.5);

    //Selection of INEL>0 Events
    //if (!IsINELgtZERO ()) return kFALSE;
    hNumberOfEvents -> Fill(2.5);

    return kTRUE;
}
//_______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskSimpleCoalescenceDeuteronInJets::IsINELgtZERO ()  {
    
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
Double_t AliAnalysisTaskSimpleCoalescenceDeuteronInJets::GetProtonWeight (Double_t pt)  {
    
    //Initialization
    Double_t wp(1);
    wp = fProtonWeights->Eval(pt);
    
    return wp;
}
//_______________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskSimpleCoalescenceDeuteronInJets::GetDeuteronWeight (Double_t pt_prot, Double_t pt_neut)  {
    
    //Initialization
    Double_t wd(1);
    
    Double_t wp = fProtonWeights->Eval(pt_prot);
    Double_t wn = fProtonWeights->Eval(pt_neut);
    
    Double_t spin_factor = 3.0/4.0;
    wd = spin_factor*wp*wn;

    return wd;
}
//_______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskSimpleCoalescenceDeuteronInJets::TwoBodyCoalescence (Double_t deltaP, Double_t p0)  {
    
    //Initialization
    Bool_t isBoundStateFormed=(kFALSE);
    
    if (deltaP<p0) isBoundStateFormed=kTRUE;
   
    return isBoundStateFormed;
}
//_______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskSimpleCoalescenceDeuteronInJets::IsInjectedParticle (AliMCParticle *particle)  {
    
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
Double_t  AliAnalysisTaskSimpleCoalescenceDeuteronInJets::Minimum (Double_t x1, Double_t x2)  {
    
    Double_t x_min(x1);
    if (x1<x2) x_min = x1;
    if (x1>x2) x_min = x2;

    return x_min;
}
//_______________________________________________________________________________________________________________________________________
TLorentzVector AliAnalysisTaskSimpleCoalescenceDeuteronInJets::LorentzTransform (TLorentzVector R, TVector3 beta_vect)  {
    
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
void AliAnalysisTaskSimpleCoalescenceDeuteronInJets::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_______________________________________________________________________________________________________________________________________
