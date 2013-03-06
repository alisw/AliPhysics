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

/*****************************************************************
  AliFlowOnTheFlyEventGenerator

  Create events with realstic spectra and v_n values
  store as AliFlowEventSimple
  
  usage: see $ALICE_ROOT/PWGCF/FLOW/macros/GenerateOnTheFlyEvents.C

  origin:       Redmer Alexander Bertens, rbertens@cern.ch, rbertens@nikhef.nl, rbertens@uu.nl
                Carlos Eugenio Perez Lara
                Andrea Dubla
******************************************************************/

#include "Riostream.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TList.h"
#include "TTree.h"
#include "TParticle.h"
#include "TMath.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TString.h"
#include "TProfile.h"
#include "TParameter.h"
#include "TBrowser.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowEventSimple.h"
#include "TRandom.h"
#include "TVirtualMCDecayer.h"
#include "AliFlowOnTheFlyEventGenerator.h"

using std::endl;
using std::cout;

static Int_t OnTheFlyEventGeneratorCounter = 0;

ClassImp(AliFlowOnTheFlyEventGenerator)
//_____________________________________________________________________________
AliFlowOnTheFlyEventGenerator::AliFlowOnTheFlyEventGenerator() : fGenerators(0), fEmbedMe(0), fFlowEvent(0), fDecayer(0), fAddV2Mothers(0), fAddV3Mothers(0), fAddV2Daughters(0), fAddV3Daughters(0), fPsi2(0), fPsi3(0), fPrecisionPhi(.001), fMaxNumberOfIterations(100), fQA(0), fFF(0) {/* constructor used for IO by root */}
//_____________________________________________________________________________
AliFlowOnTheFlyEventGenerator::AliFlowOnTheFlyEventGenerator(Bool_t qa, Int_t ff, Int_t mult, TVirtualMCDecayer* decayer, Bool_t a, Bool_t b, Bool_t c, Bool_t d) : fGenerators(0), fEmbedMe(0), fFlowEvent(0), fDecayer(0), fAddV2Mothers(0), fAddV3Mothers(0), fAddV2Daughters(0), fAddV3Daughters(0), fPsi2(0), fPsi3(0), fPrecisionPhi(.001), fMaxNumberOfIterations(100), fQA(qa), fFF(ff) {
    // contructor
    if(!fFlowEvent) fFlowEvent = new AliFlowEventSimple(mult, (AliFlowEventSimple::ConstructionMethod)0, 0x0, 0, TMath::TwoPi(), -.9, .9);
    if(!fGenerators) {
        fGenerators = new TObjArray();
        fGenerators->SetOwner(kTRUE);
        InitGenerators();
    }
    fDecayer = decayer;         // decayer: user has ownership of decayer (see dtor)
    fAddV2Mothers = a;
    fAddV3Mothers = b;
    fAddV2Daughters = c; 
    fAddV3Daughters = d;
    if(fAddV3Mothers||fAddV3Daughters) printf(" > WARNING - v3 is not impelmented yet ! \n <");
}
//_____________________________________________________________________________
AliFlowOnTheFlyEventGenerator::~AliFlowOnTheFlyEventGenerator() 
{
    // destructor
    if(fGenerators)     delete fGenerators; 
    if(fFlowEvent)      delete fFlowEvent;
}
//_____________________________________________________________________________
AliFlowOnTheFlyEventGenerator::NaiveFlowAndSpectrumGenerator* AliFlowOnTheFlyEventGenerator::Find(Short_t pdg, Bool_t make)  
{ 
    // return instance of generator class for given fPDG value
    for(int i(0); i < fGenerators->GetEntriesFast(); i++) {
        if(((NaiveFlowAndSpectrumGenerator*)fGenerators->At(i))->GetPDGCode()==pdg) {
            return (NaiveFlowAndSpectrumGenerator*)(fGenerators->At(i));
        }
    }
    if(make) {
        AliFlowOnTheFlyEventGenerator::NaiveFlowAndSpectrumGenerator* _tmp = new AliFlowOnTheFlyEventGenerator::NaiveFlowAndSpectrumGenerator(pdg, fQA, fFF);
        fGenerators->Add(_tmp);
        return _tmp;
    }
    return 0x0;
}
//_____________________________________________________________________________
void AliFlowOnTheFlyEventGenerator::SetPtSpectrum(const char* func, Short_t species)
{
    // set pt spectrum for species, override default
    NaiveFlowAndSpectrumGenerator* gen = Find(species, kTRUE);
    TF1* spectrum = gen->GetPtSpectrum();
    TString name = "";
    if(spectrum) {
        name = spectrum->GetName();
        delete spectrum;
    }
    else name = Form("pt_%i", (int)species);
    gen->SetPtSpectrum(new TF1(name.Data(), func, 0., 20));
}
//_____________________________________________________________________________
void AliFlowOnTheFlyEventGenerator::SetPtDependentV2(const char* func, Short_t species)
{
    // set pt dependent v2 for species, override default
    NaiveFlowAndSpectrumGenerator* gen = Find(species, kTRUE);
    TF1* v2 = gen->GetDifferentialV2();
    TString name = "";
    if(v2) { // if v2 is already present, get its name and release its memory
        name = v2->GetName();
        delete v2;
    }
    else name = Form("v2_%i", (int)species);
    gen->SetDifferentialV2(new TF1(name.Data(),func, 0., 20));
}
//_____________________________________________________________________________
void AliFlowOnTheFlyEventGenerator::SetPtDependentV3(const char* func, Short_t species)
{
    // set pt dependent v2 for species, override default
    NaiveFlowAndSpectrumGenerator* gen = Find(species, kTRUE);
    TF1* v3 = gen->GetDifferentialV3();
    TString name = "";
    if(v3) { // if v3 is already present, get its name and release its memory
        name = v3->GetName();
        delete v3;
    }
    else name = Form("v3_%i", (int)species);
    gen->SetDifferentialV3(new TF1(name.Data(), func, 0., 20));
}
//_____________________________________________________________________________
AliFlowEventSimple* AliFlowOnTheFlyEventGenerator::GenerateOnTheFlyEvent(TClonesArray* event, Int_t nSpecies, Int_t species[], Int_t mult[], Int_t bg, Bool_t fluc)
{
    // generate a flow event according to settings provided in GENERATOR INTERFACE
    fPsi2 = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
    // generate a random number which will be used for flow fuctuations (if requested)
    Double_t fluc_poi(-1), fluc_rp(-1);
    if(fFF == 1 || fFF == 3 ) fluc_poi = gRandom->Uniform();
    if(fFF == 2) fluc_rp = gRandom->Uniform();
    if(fFF == 3) fluc_rp = fluc_poi;    // poi and rp fluctuations are fully correlated
    // calculate the total multiplicity for the event
    Int_t multPiKPr[] = {(int)(.7*bg/2.), (int)(.2*bg/2.), (int)(.1*bg/2.)};
    Int_t fPDGPiKPr[] = {211, 321, 2212};
    Int_t totalRPMultiplicity(0), totalPOIMultiplicity(0);
    for(int i(0); i < nSpecies; i++) totalPOIMultiplicity+=mult[i];
    for(int i(0); i < 3; i++)        totalRPMultiplicity+=multPiKPr[i];
    Int_t totalMultiplicity(totalRPMultiplicity+totalPOIMultiplicity);
    // generate the particles of interest. if requested, a vn boost is given to the primary particles
    for(int i(0); i < nSpecies; i++) {
        if(fluc) GenerateOnTheFlyTracks(gRandom->Uniform(mult[i]-TMath::Sqrt(mult[i]), mult[i]+TMath::Sqrt(mult[i])), species[i], event, fluc_poi);
        else GenerateOnTheFlyTracks(mult[i], species[i], event, fluc_poi);
    }
    // generate reference particles. if requested, a vn boost is given to the primary particles
    for(int i(0); i < 3; i++) {
        if(fluc) {
            GenerateOnTheFlyTracks(gRandom->Uniform(multPiKPr[i]-TMath::Sqrt(multPiKPr[i]), multPiKPr[i]+TMath::Sqrt(mult[i])), fPDGPiKPr[i], event, fluc_rp);
            GenerateOnTheFlyTracks(gRandom->Uniform(multPiKPr[i]-TMath::Sqrt(multPiKPr[i]), multPiKPr[i]+TMath::Sqrt(mult[i])), -1*fPDGPiKPr[i], event, fluc_rp);
        }
        else {
            GenerateOnTheFlyTracks(multPiKPr[i], fPDGPiKPr[i], event, fluc_rp);
            GenerateOnTheFlyTracks(multPiKPr[i], -1*fPDGPiKPr[i], event, fluc_rp);
        }
    }
    // if requested, decay the event
    if(fDecayer)        DecayOnTheFlyTracks(event); 
    // if an embedded event is given to the generator at this point, embed it
    if(fEmbedMe) {
        event->AbsorbObjects(fEmbedMe); // event has ownership of embedded event
        fEmbedMe = 0x0;                 // reset to aviod embeddding more than once
    }
    // if requested, the secondaries are given a vn boost
    if(fAddV2Daughters) AddV2(event);
    // convert the event to an aliflowsimple event
    // the tagging (rp) is done in this step, all tracks are tagged as rp's. 
    // daughters are made orphans (i.e. primaries with secondaries are rejected from the sample)
    // converting clears the event tclones array 
    return ConvertTClonesToFlowEvent(event, totalMultiplicity);
}
//_____________________________________________________________________________
AliFlowEventSimple* AliFlowOnTheFlyEventGenerator::ConvertTClonesToFlowEvent(TClonesArray* event, Int_t totalMultiplicity)
{
    // convert the content of a tclones array to an aliflowevent simple for output and tag particles as poi's and rp's
    // first step, see if there's already a flow event available in memory. if not create one
    if(!fFlowEvent) fFlowEvent = new AliFlowEventSimple(totalMultiplicity, (AliFlowEventSimple::ConstructionMethod)0, 0x0, 0, TMath::TwoPi(), -.9, .9);
    // if there is a flow event, clear the members without releasing the allocated memory
    else fFlowEvent->ClearFast();
    fFlowEvent->SetMCReactionPlaneAngle(fPsi2);
    // prepare a track
    TParticle* pParticle = 0x0;
    Int_t iSelParticlesRP(0);
    for (Int_t i(0); i < event->GetEntries(); i++) {    // begin the loop over the input tracks
        pParticle = (TParticle*)event->At(i);           // get the particle 
        if (!pParticle) continue;                       // skip if empty slot (no particle)
        if (pParticle->GetNDaughters()!=0) continue;    // see if the particle has daughters (if so, reject it)      
        AliFlowTrackSimple* pTrack = new AliFlowTrackSimple(pParticle);         // allocate space for the new flow track
        pTrack->SetWeight(pParticle->Pz());                                     // ugly hack: store pz here ...
        pTrack->SetID(pParticle->GetPdgCode());                                 // set pid code as id
        pTrack->SetForRPSelection(kTRUE);                                       // tag ALL particles as RP's, 
        // ugly hack 2: set charge to -1 for primaries and 1 for secondaries
        if(pParticle->GetFirstMother()==-1) pTrack->SetCharge(-1);
        else pTrack->SetCharge(1);
        iSelParticlesRP++;
        fFlowEvent->AddTrack(pTrack);
    }
    fFlowEvent->SetNumberOfRPs(iSelParticlesRP);
    // all trakcs have safely been copied so we can clear the event
    event->Clear();
    OnTheFlyEventGeneratorCounter = 0;
    return fFlowEvent;
}
//_____________________________________________________________________________
void AliFlowOnTheFlyEventGenerator::AddV2(TParticle* part, Double_t v2, Double_t fluc)
{
    // afterburner, adds v2, uses Newton-Raphson iteration
    Double_t phi = part->Phi();
    Double_t phi0=phi;
    Double_t f=0.;
    Double_t fp=0.;
    Double_t phiprev=0.;
    // introduce flow fluctuations (gaussian)
    if(fluc > -.5) {
        v2+=TMath::Sqrt(2*(v2*.25)*(v2*.25))*TMath::ErfInverse(2*fluc-1);
        if(fQA) Find(part->GetPdgCode(), kTRUE)->FillV2(part->Pt(), v2);
    }
    for (Int_t i=0; i!=fMaxNumberOfIterations; ++i) {
        phiprev=phi; //store last value for comparison
        f =  phi-phi0+v2*TMath::Sin(2.*(phi-fPsi2));
        fp = 1.0+2.0*v2*TMath::Cos(2.*(phi-fPsi2)); //first derivative
        phi -= f/fp;
        if (TMath::AreEqualAbs(phiprev,phi,fPrecisionPhi)) break;
    }
    part->SetMomentum( part->Pt()*TMath::Cos(phi), part->Pt()*TMath::Sin(phi), part->Pz(), part->Energy() );
}
//_____________________________________________________________________________
void AliFlowOnTheFlyEventGenerator::AddV2(TClonesArray* event)
{
    // afterburner, adds v2 for different species to all tracks in an event
    Double_t fluc(gRandom->Uniform());  // get a random number in case of fluctuations
    // FIXME at the moment no distincition between mothers and daughters
    TParticle *part;
    for(Int_t nTrack=0; nTrack!=event->GetEntriesFast(); ++nTrack) {
        part = (TParticle*) event->At(nTrack);
        // for each track, call overloaded addv2 function with the correct generator
        // create a generator in the case where the decayer has introduced a new particle species
	if(part)
	  AddV2(part, Find(part->GetPdgCode(), kTRUE)->GetV2(part->Pt()), fluc);
    }
}
//_____________________________________________________________________________
void AliFlowOnTheFlyEventGenerator::GenerateOnTheFlyTracks(Int_t mult, Int_t pid, TClonesArray* event, Double_t fluc) 
{
    // generate tracks for the event. the tracks that are generated here are the mothers
    Double_t ph, et, pt, mm, px, py, pz, ee;
    // access the descired generator through the static find method. if the particle species is unknown, we add it to the generator
    NaiveFlowAndSpectrumGenerator* generator = Find(pid, kTRUE);
    Int_t _tempCounter = OnTheFlyEventGeneratorCounter;
    for(int i=_tempCounter; i<=mult+_tempCounter; ++i) {
        OnTheFlyEventGeneratorCounter++;
        TParticle* part = (TParticle*)event->ConstructedAt(i);
        part->SetStatusCode(1);
        part->SetFirstMother(-1);
        part->SetLastMother(-1);
        part->SetFirstDaughter(-1);
        part->SetLastDaughter(-1);
        ph = gRandom->Uniform(0,TMath::TwoPi());
        et = gRandom->Uniform(-1.0,+1.0);
        pt = generator->GetPt();
        part->SetPdgCode( pid );
        mm =  part->GetMass();
        px = pt*TMath::Cos(ph);
        py = pt*TMath::Sin(ph);
        pz = pt*TMath::SinH(et);
        ee = TMath::Sqrt( mm*mm + pt*pt+pz*pz );
        part->SetMomentum( px, py, pz, ee );
        // if requested add v2 to the sample of mothers
        if(fAddV2Mothers) AddV2(part, generator->GetV2(pt), fluc);
    }
}
//_____________________________________________________________________________
void AliFlowOnTheFlyEventGenerator::DecayOnTheFlyTracks(TClonesArray *event) 
{
    // decay the tracks using a decayer that is set in the GenerateEventsOnTheFly.C macro
    if(!fDecayer) return;       // shouldn't happen ...
    fDecayer->Init();
    Int_t kf;
    TClonesArray*  arr = new TClonesArray("TParticle",10);
    Int_t secondaries=0;
    Int_t nStart=0;
    Int_t nTracks = event->GetEntriesFast();
    TParticle *part, *part1;
    TClonesArray &in = *event;
    for(Int_t loop=0; loop!=3; ++loop) {
        secondaries=0;
        for(Int_t nTrack=nStart; nTrack!=nTracks; ++nTrack) {
            arr->Clear();
            part = (TParticle*) event->At( nTrack );
            if(!part) continue;
            kf = TMath::Abs(part->GetPdgCode());
            if(kf==22 )  ForceGammaDecay(arr, part);     // if gamma, we decay it by hand
            else {      // else we let pythia do it
                TLorentzVector pmom( part->Px(), part->Py(), part->Pz(), part->Energy() );
                fDecayer->Decay(kf,&pmom);
                fDecayer->ImportParticles(arr);
            }
            Int_t nDaughters = arr->GetEntries();
            if(nDaughters<2) continue;
            for(Int_t nDaughter=1; nDaughter!=nDaughters; ++nDaughter) {
                part1 = (TParticle*) (arr->At(nDaughter));
                if(!part1) continue; // was part1, mistake?
	        new(in[nTracks+secondaries]) TParticle( part1->GetPdgCode(),
						        part1->GetStatusCode(),
						        part1->GetFirstMother(),
						        part1->GetSecondMother(),
						        part1->GetFirstDaughter(),
						        part1->GetLastDaughter(),
						        part1->Px(),
						        part1->Py(),
						        part1->Pz(),
						        part1->Energy(),
						        part1->Vx(),
						        part1->Vy(),
						        part1->Vz(),
						        part1->T());
	        secondaries++;
                if(nDaughter==1) part->SetFirstDaughter(nTracks+secondaries);
                else if ((nDaughters-1)==nDaughter) part->SetLastDaughter(nTracks+secondaries);
                //else part->SetDaughter(nDaughter,nTracks+secondaries);
                }
            }
        nStart = nTracks;
        nTracks = event->GetEntries();
    }
    delete arr;
    return;
}
//_____________________________________________________________________________
void AliFlowOnTheFlyEventGenerator::ForceGammaDecay(TClonesArray* arr, TParticle* part)
{
    // force gama decay
    const Float_t mass1(gRandom->Uniform(0.0000001, .1));   // FIXME randomly sample photon 'mass'
    const Float_t mass_d1(.000511);                         // daughter mass (electron)
    const Float_t mass_d2(.000511);
    Double_t p_dec1=sqrt((TMath::Abs(mass1*mass1-(mass_d1+mass_d2)*(mass_d1+mass_d2))*(mass1*mass1-(mass_d1-mass_d2)*(mass_d1-mass_d2))))/2/mass1;      // energy
    TLorentzVector *p_d1=new TLorentzVector();              // lorentz vectors for daughters
    TLorentzVector *p_d2=new TLorentzVector();
    Double_t yy_m=part->Y();                                // pseudo rapidity of mother
    Double_t pt_m=part->Pt();                               // pt of mother
    Double_t mt_m=sqrt(mass1*mass1+pt_m*pt_m);              // transverse mass of mother
    Double_t pz_m=mt_m*TMath::SinH(yy_m);                   // MAGIC
    Double_t phi_m=gRandom->Rndm()*2*TMath::Pi();
    Double_t px_m=pt_m*TMath::Cos(phi_m);
    Double_t py_m=pt_m*TMath::Sin(phi_m);
    Double_t costh_d=2*gRandom->Rndm()-1;
    Double_t phi_d=gRandom->Rndm()*2*TMath::Pi();
    Double_t pz_d=p_dec1*costh_d;
    Double_t pt_d=p_dec1*sqrt(1-costh_d*costh_d);
    Double_t px_d=pt_d*TMath::Cos(phi_d);
    Double_t py_d=pt_d*TMath::Sin(phi_d);
    p_d1->SetPxPyPzE(px_d,py_d,pz_d,sqrt(mass_d1*mass_d1+p_dec1*p_dec1));
    p_d2->SetPxPyPzE(-px_d,-py_d,-pz_d,sqrt(mass_d2*mass_d2+p_dec1*p_dec1));
    Double_t gamma_b=sqrt(mass1*mass1+pz_m*pz_m+pt_m*pt_m)/mass1;
    Double_t bx=px_m/gamma_b/mass1;    
    Double_t by=py_m/gamma_b/mass1;
    Double_t bz=pz_m/gamma_b/mass1;
    p_d1->Boost(bx,by,bz);
    p_d2->Boost(bx,by,bz);
    TParticle* daughter_1 = (TParticle*)arr->ConstructedAt(0);
    TParticle* daughter_2 = (TParticle*)arr->ConstructedAt(1);
    daughter_1->SetStatusCode(1);
    daughter_1->SetFirstMother(-1);
    daughter_1->SetLastMother(-1);
    daughter_1->SetFirstDaughter(-1);
    daughter_1->SetLastDaughter(-1);
    daughter_1->SetPdgCode(11);
    TLorentzVector& ref_p_d1 = *p_d1;
    daughter_1->SetMomentum(ref_p_d1);
    daughter_2->SetStatusCode(1);
    daughter_2->SetFirstMother(-1);
    daughter_2->SetLastMother(-1);
    daughter_2->SetFirstDaughter(-1);
    daughter_2->SetLastDaughter(-1);
    daughter_2->SetPdgCode(-11);
    TLorentzVector& pp = *p_d2;
    daughter_2->SetMomentum(pp);
}
//_____________________________________________________________________________
void AliFlowOnTheFlyEventGenerator::InitGenerators()
{
   // initializes generators for a list of common particles
   // new generators will be created 'on the fly' if they are encountered
   // if a generator is requested for a non-existent pdg value, things will get unstable ...
   Short_t PDGcode[] = {11, -11, // e-               
                        12, // v_e
                        13, -13, // mu-
                        14, // v_mu
                        15, -15, // t-
                        16, // v_t
                        21, // g
                        22, // gamma
                        111, // pi0
                        211, -211, // pi+
                        113, // rho0
                        213, -213, // rho+
                        221, // eta
                        331, // eta prime
                        130, // k0l
                        310, // k0s
                        311, // k0
                        313, // k*
                        321, -321, // k+
                        323, -323, // k*+
                        223, // omega
                        411, -411, // d+
                        413, -413, // d*+
                        421, // d0
                        423, // d*0
                        333, // phi
                        443, // jpsi 
                        2112, // neutron
                        2212, -2212, // proton
                        3122, // lambda0
                        3312, -3312, // xi- (cascade)
                        3314, -3314, // xi*-
                        3322, // xi0
                        3324, // xi*0
                        3334}; // Omeg
   for(int i(0); i < 47; i++) fGenerators->Add(new NaiveFlowAndSpectrumGenerator(PDGcode[i], fQA, fFF));
}
//_____________________________________________________________________________
void AliFlowOnTheFlyEventGenerator::PrintGenerators()
{
   // make sure all went ok
   if(!fGenerators) {
       printf(" > PANIC <\n");
       return;
   }
   printf(" > %i species available in generator\n", fGenerators->GetEntriesFast());
   for(int i(0); i < fGenerators->GetEntriesFast(); i++) 
       printf("   - PDG: %i \n", ((NaiveFlowAndSpectrumGenerator*)fGenerators->At(i))->GetPDGCode());
   printf(" more can be created 'on the fly', but beware of non-existent pdg values !");
}
//_____________________________________________________________________________
void AliFlowOnTheFlyEventGenerator::DoGeneratorQA(Bool_t v2, Bool_t v3) 
{
    // this function loops over all the species that are available in the fGenerators, which are
    // a) created by the user (in the event generator) as mother particles or
    // b) introduced by pythia (the decayer) as daughter particles
    // c) made by the InitGenerators() function
    // 'empty species' (i.e. species for which a generator exists but which were not actually sampled) are omitted
    // the output histograms are written to a file named "GeneratorfQA.root"
    if(!fQA) {
        printf(" > Request has been made for QA plots but QA histograms have been disabled !\n");
        return;
    }
    TFile *QAfile = new TFile("GeneratorfQA.root","RECREATE"); 
    for(int i(0); i < fGenerators->GetEntriesFast(); i++) {
        NaiveFlowAndSpectrumGenerator* gen = (NaiveFlowAndSpectrumGenerator*)fGenerators->At(i);
        if(!gen) continue;               // skip if we failed to get a generator
        TH1F* QApt = (TH1F*)gen->GetQAType(0); 
        TH2F* QAv2 = (TH2F*)gen->GetQAType(1);
        TH2F* QAv3 = (TH2F*)gen->GetQAType(2);
        if((!QApt)||(!QAv2)||(!QAv3)) {
            printf(" > couldn't read qa histogrmas for species %i <\n", gen->GetPDGCode());
            continue;
        }
        if(QApt->GetEntries()==0&&QAv2->GetEntries()==0&&QAv3->GetEntries()==0) continue; // skip if no tracks of this species have been sampled
        printf(" > saving QA histograms for sampled species %i <\n", gen->GetPDGCode());
        QAfile->mkdir(Form("fPDG_%i", gen->GetPDGCode()));        // create a directory in the output file
        QAfile->cd(Form("fPDG_%i", gen->GetPDGCode()));              // cd into this directory
        if(!QApt->GetEntries()==0) { // only plot the pt fSpectrum if this guy was generated
            if(!gen->GetPtSpectrum()->Integral(0,20)==0) {
                QApt->Scale(gen->GetPtSpectrum()->Integral(0,20)/QApt->Integral(),"width");
                QApt->Write();
            }
        }
        gen->GetPtSpectrum()->SetNpx(10000); // otherwise tf1 plot is very ugly
        gen->GetPtSpectrum()->Draw();
        gen->GetPtSpectrum()->Write();
        if(v2) {
            QAv2->Draw();
            gen->GetDifferentialV2()->SetNpx(10000);
            gen->GetDifferentialV2()->Draw();
            gen->GetDifferentialV2()->Write();
            QAv2->Write();
        }
        if(v3) {
            QAv3->Draw();
            gen->GetDifferentialV3()->Draw();
            gen->GetDifferentialV3()->SetNpx(10000);
            gen->GetDifferentialV3()->Write();
            QAv3->Write();
        }
    }
    QAfile->Close();
}
//_____________________________________________________________________________

