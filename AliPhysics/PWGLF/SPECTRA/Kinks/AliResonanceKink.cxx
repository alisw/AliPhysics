/**************************************************************************
 * Author: Paraskevi Ganoti, University of Athens (pganoti@phys.uoa.gr)   *
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
 
//----------------------------------------------------------------------------------------------------------------
//                        class AliResonanceKink
//        Example of an analysis task for reconstructing resonances having at least one kaon-kink in their decay 
//        products. It provides basic plots as well as plots helping to calculate the corrections.
//        Usage: To analyse a resonance having a kaon in its decay products, one has to modify the integer 
//        variables resonancePDG, fdaughter1pdg and fdaughter2pdg accordingly as well as daughter1pdgMass  and daughter2pdgMass.
//        Also, depending on the analysis mode (ESD or MC), fAnalysisType in the constructor must also be changed 
//-----------------------------------------------------------------------------------------------------------------

#include "TH2D.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TF1.h"
#include "TList.h"

#include "AliMCEvent.h"
#include "AliResonanceKink.h"
#include "AliESDkink.h"
#include "AliStack.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliExternalTrackParam.h"
#include "AliMCParticle.h"

ClassImp(AliResonanceKink)

//________________________________________________________________________
AliResonanceKink::AliResonanceKink() 
  : TObject(), fDebug(0), fListOfHistos(0), fOpeningAngle(0), fInvariantMass(0), fInvMassTrue(0), fPhiBothKinks(0), fetabins(0), floweta(0), fupeta(0), fRecPt(0), fRecEta(0), fRecEtaPt(0), fSimPt(0), fSimEta(0), fSimEtaPt(0), fSimPtKink(0), fSimEtaKink(0),  fSimEtaPtKink(0), f1(0), f2(0), fAnalysisType(), fvtxz(0), fNbins(0), fLowX(0), fHighX(0), fdaughter1pdg(0), fdaughter2pdg(0), fresonancePDGcode(0), fMaxNSigmaToVertex(0), fMinPtTrackCut(0), fMaxDCAxy(0), fMaxDCAzaxis(0), 
fMinTPCclusters(0),fMaxChi2PerTPCcluster(0), fMaxCov0(0), fMaxCov2(0), fMaxCov5(0) , fMaxCov9(0), fMaxCov14(0) //, fTPCrefitFlag(kFALSE)
, fInvmassPt(0), fInvmassPtTrue(0), fMCInvmassPt(0), fMCInvmassPtTrue(0), fminKinkRadius(0), fmaxKinkRadius(0), fminQt(0), fmaxQt(0), fptbins(0), flowpt(0), fupperpt(0), fmaxAbsEtaCut(0)
{
  // Constructor
}

//________________________________________________________________________
AliResonanceKink::AliResonanceKink(Int_t nbins, Float_t nlowx, Float_t nhighx, Int_t netabins, Float_t nloweta, Float_t nupeta, Int_t nptbins, Float_t nlowpt, Float_t nupperpt, Int_t daughter1, Int_t daughter2, Int_t resonancePDG) 
  : TObject(), fDebug(0), fListOfHistos(0), fOpeningAngle(0), fInvariantMass(0), fInvMassTrue(0), fPhiBothKinks(0), fetabins(netabins), floweta(nloweta), fupeta(nupeta), fRecPt(0), fRecEta(0), fRecEtaPt(0), fSimPt(0), fSimEta(0), fSimEtaPt(0), fSimPtKink(0), fSimEtaKink(0),  fSimEtaPtKink(0), f1(0), f2(0), fAnalysisType(), fvtxz(0), fNbins(nbins), fLowX(nlowx), fHighX(nhighx), fdaughter1pdg(daughter1), fdaughter2pdg(daughter2), fresonancePDGcode(resonancePDG), fMaxNSigmaToVertex(0), fMinPtTrackCut(0), 
fMaxDCAxy(0), fMaxDCAzaxis(0), fMinTPCclusters(0), fMaxChi2PerTPCcluster(0), fMaxCov0(0), fMaxCov2(0), fMaxCov5(0), fMaxCov9(0), fMaxCov14(0) //, fTPCrefitFlag(kFALSE)
, fInvmassPt(0), fInvmassPtTrue(0), fMCInvmassPt(0), fMCInvmassPtTrue(0), fminKinkRadius(0), fmaxKinkRadius(0), fminQt(0), fmaxQt(0), fptbins(nptbins), flowpt(nlowpt), fupperpt(nupperpt), fmaxAbsEtaCut(0)
{
   // Constructor
  
   fOpeningAngle=new TH1D("fOpeningAngle"," ", 100,-1.0,1.0);

   fInvariantMass=new TH1D("fInvariantMass"," ",fNbins,fLowX,fHighX);
   fInvMassTrue=new TH1D("fInvMassTrue"," ",fNbins,fLowX,fHighX);
   fPhiBothKinks=new TH1D("fPhiBothKinks"," ",fNbins,fLowX,fHighX);  // Applicable for phi(1020)

   fRecPt=new TH1D("fRecPt"," ", nptbins, nlowpt, nupperpt);
   fRecEta=new TH1D("fRecEta"," ",netabins, nloweta, nupeta);
   fRecEtaPt=new TH2D("fRecEtaPt"," ", nptbins, nlowpt, nupperpt, netabins, nloweta, nupeta); 
   fSimPt=new TH1D("fSimPt"," ", nptbins, nlowpt, nupperpt);
   fSimEta=new TH1D("fSimEta"," ", netabins, nloweta, nupeta); 
   fSimEtaPt=new TH2D("fSimEtaPt"," ", nptbins, nlowpt, nupperpt, netabins, nloweta, nupeta);
   fSimPtKink=new TH1D("fSimPtKink"," ", nptbins, nlowpt, nupperpt);
   fSimEtaKink=new TH1D("fSimEtaKink"," ", netabins, nloweta, nupeta);
   fSimEtaPtKink=new TH2D("fSimEtaPtKink"," ", nptbins, nlowpt, nupperpt, netabins, nloweta, nupeta);                
   
   f1=new TF1("f1","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",1.1,10.0);
   f1->SetParameter(0,0.493677);
   f1->SetParameter(1,0.9127037);
   f1->SetParameter(2,TMath::Pi());

   f2=new TF1("f2","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",0.1,10.0);
   f2->SetParameter(0,0.13957018);
   f2->SetParameter(1,0.2731374);
   f2->SetParameter(2,TMath::Pi());
   
   fvtxz=new TH1D("fvtxz"," ", 100,-20.0,20.0);
   
   fInvmassPt=new TH2D("fInvmassPt"," ",fNbins,fLowX,fHighX,100,0.0,10.0);
   fInvmassPtTrue=new TH2D("fInvmassPtTrue"," ",fNbins,fLowX,fHighX,100,0.0,10.0);  
   fMCInvmassPt=new TH2D("fMCInvmassPt"," ",fNbins,fLowX,fHighX,100,0.0,10.0);
   fMCInvmassPtTrue=new TH2D("fMCInvmassPtTrue"," ",fNbins,fLowX,fHighX,100,0.0,10.0); 
}

//________________________________________________________________________
AliResonanceKink:: ~AliResonanceKink()
{
 //  Destructor
 if(fOpeningAngle) delete fOpeningAngle;
 if(fInvariantMass) delete fInvariantMass;
 if(fInvMassTrue) delete fInvMassTrue;
 if(fPhiBothKinks) delete fPhiBothKinks;
 if(fRecPt) delete fRecPt;
 if(fRecEta) delete fRecEta;
 if(fRecEtaPt) delete fRecEtaPt;
 if(fSimPt) delete fSimPt;
 if(fSimEta) delete fSimEta;
 if(fSimEtaPt) delete fSimEtaPt;
 if(fSimPtKink) delete fSimPtKink;
 if(fSimEtaKink) delete fSimEtaKink;
 if(fSimEtaPtKink) delete fSimEtaPtKink;
 if(fvtxz) delete fvtxz;  
 if(fInvmassPt) delete fInvmassPt; 
 if(fInvmassPtTrue) delete fInvmassPtTrue; 
 if(fMCInvmassPt) delete fMCInvmassPt; 
 if(fMCInvmassPtTrue) delete fMCInvmassPtTrue;            
}
//________________________________________________________________________
TList* AliResonanceKink::GetHistogramList()
{
  // Adding histograms to the list
  fListOfHistos=new TList();
 
  fListOfHistos->Add(fOpeningAngle);
  fListOfHistos->Add(fInvariantMass);
  fListOfHistos->Add(fInvMassTrue);
  fListOfHistos->Add(fPhiBothKinks);
  fListOfHistos->Add(fRecPt);    
  fListOfHistos->Add(fRecEta);   
  fListOfHistos->Add(fRecEtaPt);    
  fListOfHistos->Add(fSimPt);    
  fListOfHistos->Add(fSimEta);   
  fListOfHistos->Add(fSimEtaPt);     
  fListOfHistos->Add(fSimPtKink);    
  fListOfHistos->Add(fSimEtaKink);   
  fListOfHistos->Add(fSimEtaPtKink);                                                           
  fListOfHistos->Add(fvtxz);
  fListOfHistos->Add(fInvmassPt);
  fListOfHistos->Add(fInvmassPtTrue);
  fListOfHistos->Add(fMCInvmassPt);
  fListOfHistos->Add(fMCInvmassPtTrue);     
     
  return fListOfHistos;
}

//________________________________________________________________________
void AliResonanceKink::InitOutputHistograms(Int_t nbins, Float_t nlowx, Float_t nhighx, Int_t netabins, Float_t nloweta, Float_t nupeta, Int_t nptbins, Float_t nlowpt, Float_t nupperpt)
{
  //  Initialisation of the output histograms
  fNbins=nbins; 
  fLowX=nlowx; 
  fHighX=nhighx;
    
  fOpeningAngle=new TH1D("fOpeningAngle"," ", 100,-1.0,1.0);

  fInvariantMass=new TH1D("fInvariantMass"," ",fNbins,fLowX,fHighX);
  fInvMassTrue=new TH1D("fInvMassTrue"," ",fNbins,fLowX,fHighX);
  fPhiBothKinks=new TH1D("fPhiBothKinks"," ",fNbins,fLowX,fHighX);  // Applicable for phi(1020)

  fRecPt=new TH1D("fRecPt"," ", nptbins, nlowpt, nupperpt);
  fRecEta=new TH1D("fRecEta"," ", netabins, nloweta, nupeta);
  fRecEtaPt=new TH2D("fRecEtaPt"," ", nptbins, nlowpt, nupperpt, netabins, nloweta, nupeta); 
  fSimPt=new TH1D("fSimPt"," ", nptbins, nlowpt, nupperpt);
  fSimEta=new TH1D("fSimEta"," ",netabins, nloweta, nupeta); 
  fSimEtaPt=new TH2D("fSimEtaPt"," ", nptbins, nlowpt, nupperpt, netabins, nloweta, nupeta);
  fSimPtKink=new TH1D("fSimPtKink"," ", nptbins, nlowpt, nupperpt);
  fSimEtaKink=new TH1D("fSimEtaKink"," ", netabins, nloweta, nupeta);
  fSimEtaPtKink=new TH2D("fSimEtaPtKink"," ", nptbins, nlowpt, nupperpt, netabins, nloweta, nupeta);                
   
  f1=new TF1("f1","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",1.1,10.0);
  f1->SetParameter(0,0.493677);
  f1->SetParameter(1,0.9127037);
  f1->SetParameter(2,TMath::Pi());

  f2=new TF1("f2","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",0.1,10.0);
  f2->SetParameter(0,0.13957018);
  f2->SetParameter(1,0.2731374);
  f2->SetParameter(2,TMath::Pi());
   
  fvtxz=new TH1D("fvtxz"," ", 100,-20.0,20.0);
  
  fInvmassPt=new TH2D("fInvmassPt"," ",fNbins,fLowX,fHighX, nptbins, nlowpt, nupperpt);
  fInvmassPtTrue=new TH2D("fInvmassPtTrue"," ",fNbins,fLowX,fHighX, nptbins, nlowpt, nupperpt);  
  fMCInvmassPt=new TH2D("fMCInvmassPt"," ",fNbins,fLowX,fHighX,nptbins, nlowpt, nupperpt);
  fMCInvmassPtTrue=new TH2D("fMCInvmassPtTrue"," ",fNbins,fLowX,fHighX, nptbins, nlowpt, nupperpt);
  
}

//________________________________________________________________________
void AliResonanceKink::Analyse(AliESDEvent* esd, AliMCEvent* mcEvent) 
{
  // Main loop
  // Called for each event
  Int_t resonancePDGcode, antiresonancePDGcode;
  Double_t daughter1pdgMass, daughter2pdgMass;
  
  if (fdaughter1pdg==kKPlus)  {
    resonancePDGcode=fresonancePDGcode;
    antiresonancePDGcode=-fresonancePDGcode;
    daughter1pdgMass=TDatabasePDG::Instance()->GetParticle(fdaughter1pdg)->Mass();
    daughter2pdgMass=TDatabasePDG::Instance()->GetParticle(fdaughter2pdg)->Mass(); 
  }
  
  if (fdaughter1pdg!=kKPlus)  {
    resonancePDGcode=-fresonancePDGcode;
    antiresonancePDGcode=fresonancePDGcode;
    daughter1pdgMass=TDatabasePDG::Instance()->GetParticle(fdaughter2pdg)->Mass();
    daughter2pdgMass=TDatabasePDG::Instance()->GetParticle(fdaughter1pdg)->Mass();   
  }  // to ensure that daughter1pdgMass has always the kaon mass
  
  if (fdaughter1pdg==fdaughter2pdg)  {
    resonancePDGcode=fresonancePDGcode;
    antiresonancePDGcode=fresonancePDGcode;
  }  
  if (!esd) {
    Printf("ERROR: esd not available");
    return;
  }
      if (!mcEvent) {
      return;
     }  
  
  AliStack* stack=mcEvent->Stack();

  if(fAnalysisType == "MC") {

  const AliESDVertex* vertex = GetEventVertex(esd);
  if(!vertex) return;
  Double_t vtx[3];
  vertex->GetXYZ(vtx);
  fvtxz->Fill(vtx[2]);

  for (Int_t iMc = 0; iMc < stack->GetNprimary(); ++iMc)
  {
    TParticle* particle = stack->Particle(iMc);

    if (!particle)
    {
      if (fDebug > 0) Printf("UNEXPECTED: particle with label %d not found in stack (mc loop)", iMc);
      continue;
    }

       if(TMath::Abs(particle->GetPdgCode())==fresonancePDGcode) {
       
       if(particle->GetNDaughters()>2) continue;
      
       Int_t firstD=particle->GetFirstDaughter();
       Int_t lastD=particle->GetLastDaughter();

       if((firstD<0)||(firstD>=stack->GetNtrack())) continue;
       if((lastD<0)||(lastD>=stack->GetNtrack())) continue;
       
       TParticle *daughterParticle1 = 0;
       TParticle *daughterParticle2 = 0;
       AliMCParticle   *mcDaughter1 = 0;
       AliMCParticle   *mcDaughter2 = 0;       
               
        if(fdaughter1pdg==kKPlus) {  
          daughterParticle1=stack->Particle(firstD);
          daughterParticle2=stack->Particle(lastD);
	  mcDaughter1= (AliMCParticle*) mcEvent->GetTrack(firstD);  
	  mcDaughter2= (AliMCParticle*) mcEvent->GetTrack(lastD);  	  
        }
        else 
	    if(fdaughter2pdg==kKPlus) {
              daughterParticle1=stack->Particle(lastD);
              daughterParticle2=stack->Particle(firstD); 
	      mcDaughter1= (AliMCParticle*) mcEvent->GetTrack(lastD);
	      mcDaughter2= (AliMCParticle*) mcEvent->GetTrack(firstD); 
             }    //to ensure that the first daughter is always the kaon
	     
       if(TMath::Abs(daughterParticle1->GetPdgCode())!=321) continue;
       
       TParticle      * daughters1Daughter=0;
       TParticle      * daughters2Daughter=0;       
       Int_t mcProcessDaughters1Daughter = -999;
       Int_t mcProcessDaughters2Daughter = -999;       
       AliMCParticle *mcDaughters1Daughter = 0;
       
       if(mcDaughter1->Charge()==0) continue;
       if(mcDaughter2->Charge()==0) continue;       //accept resonance decays in two charged daughters
       Int_t nDecayKaonDaughter=-99; 
       for(Int_t ia=0; ia<daughterParticle1->GetNDaughters(); ia++) {
       if(((daughterParticle1->GetFirstDaughter()+ia)>0)&&((daughterParticle1->GetFirstDaughter()+ia)<stack->GetNtrack())) {
	   daughters1Daughter=stack->Particle(daughterParticle1->GetFirstDaughter()+ia);
           mcProcessDaughters1Daughter=daughters1Daughter->GetUniqueID(); 
	   if(mcProcessDaughters1Daughter==4) {
	     nDecayKaonDaughter=daughterParticle1->GetFirstDaughter()+ia;
	     break;
	   }
        }
       }
       
       Int_t nProcessDaughter=-99; 
       for(Int_t ib=0; ib<daughterParticle2->GetNDaughters(); ib++) {
        if(((daughterParticle2->GetFirstDaughter()+ib)>0)&&((daughterParticle2->GetFirstDaughter()+ib)<stack->GetNtrack())) {
	   daughters2Daughter=stack->Particle(daughterParticle2->GetFirstDaughter()+ib);
           mcProcessDaughters2Daughter=daughters2Daughter->GetUniqueID(); 
	   if((mcProcessDaughters2Daughter==4)||(mcProcessDaughters2Daughter==13)) {
	     nProcessDaughter=mcProcessDaughters2Daughter;     
	     break;
	   }  
        }
       }  
       
       Int_t numberOfCharged=0;
       if((mcProcessDaughters1Daughter==4)&&(nDecayKaonDaughter>=0)) {
         for(Int_t ic=nDecayKaonDaughter; ic<=daughterParticle1->GetLastDaughter(); ic++) {
	   if ((ic>=0)&&(ic<stack->GetNtrack())) mcDaughters1Daughter= dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(ic));
	    else continue;
	    if(mcDaughters1Daughter && mcDaughters1Daughter->Charge()!=0) numberOfCharged=numberOfCharged+1;
         }
       }
       
       	 if(numberOfCharged>=2) continue; // leave out kaon decays to more than one charged daughter
         if ((particle->Pt()>fMinPtTrackCut)&&(TMath::Abs(particle->Eta())<fmaxAbsEtaCut)) {
	  
         fSimEta->Fill(particle->Eta());
	 fSimEtaPt->Fill(particle->Pt(), particle->Eta());
	 fSimPt->Fill(particle->Pt());
         fMCInvmassPtTrue->Fill(particle->GetMass(), particle->Pt());	   

         if((daughterParticle1->Pt()>fMinPtTrackCut)&&(TMath::Abs(daughterParticle1->Eta())<fmaxAbsEtaCut)&&(daughterParticle2->Pt()>fMinPtTrackCut)&&(TMath::Abs(daughterParticle2->Eta())<fmaxAbsEtaCut)) {
	   if((mcProcessDaughters1Daughter==4)&&(daughters1Daughter->R()>fminKinkRadius)&&(daughters1Daughter->R()<fmaxKinkRadius)&&( (nProcessDaughter<0)||((daughters2Daughter->R()>fminKinkRadius)&&(nProcessDaughter>0)))) { //below are the findable
	     fSimPtKink->Fill(particle->Pt());
	     fSimEtaKink->Fill(particle->Eta());
	     fSimEtaPtKink->Fill(particle->Pt(), particle->Eta());
	     fMCInvmassPt->Fill(particle->GetMass(), particle->Pt());
	   }
         }
       
       } // for the generated spectra 
             
          }   //for the particular resonance
  } //MC loop
  
  } // end fAnalysisType==MC 
  
  if(fAnalysisType == "ESD") {
  const AliESDVertex* vertex = GetEventVertex(esd);
  if(!vertex) return;
  Double_t vtx[3];
  vertex->GetXYZ(vtx);
  fvtxz->Fill(vtx[2]);
  Double_t ptrackpos[3], ptrackneg[3];
  
  TLorentzVector p4pos, anp4pos;
  TLorentzVector p4neg, anp4neg;
  TLorentzVector p4comb, anp4comb;
  
  for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* trackpos = esd->GetTrack(iTracks);
    if (!trackpos) {
      if (fDebug > 0) Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    if (trackpos->GetSign() < 0) continue;
    
    AliExternalTrackParam *tpcTrackpos = (AliExternalTrackParam *)trackpos->GetTPCInnerParam();
    if(!tpcTrackpos) continue;
    ptrackpos[0]=tpcTrackpos->Px();
    ptrackpos[1]=tpcTrackpos->Py();   
    ptrackpos[2]=tpcTrackpos->Pz();  
    
    Bool_t firstLevelAcceptPosTrack=IsAcceptedForKink(esd, vertex, trackpos);
    if(firstLevelAcceptPosTrack==kFALSE) continue;
    
    TVector3 posTrackMom(ptrackpos[0],ptrackpos[1],ptrackpos[2]);
  	
    TParticle * partpos = stack->Particle(TMath::Abs(trackpos->GetLabel()));
    if (!partpos) continue;
    Int_t pdgpos = partpos->GetPdgCode();
    Int_t mumlabelpos=partpos->GetFirstMother();
    mumlabelpos = TMath::Abs(mumlabelpos);
    TParticle * mumpos=stack->Particle(mumlabelpos);
    if (!mumpos) continue;
    Int_t mumpdgpos = mumpos->GetPdgCode();
    
    Int_t indexKinkPos=trackpos->GetKinkIndex(0);
    
    if(indexKinkPos>0) continue;
    
    Bool_t posKaonKinkFlag=0;
    
    if(indexKinkPos<0) {
      posKaonKinkFlag=IsKink(esd, indexKinkPos, posTrackMom);
    
      if(posKaonKinkFlag==1) anp4pos.SetVectM(posTrackMom,daughter1pdgMass);
      if(posKaonKinkFlag==0) continue;
    }
    
    if(indexKinkPos==0) {

    Bool_t secondLevelAcceptPosTrack=IsAcceptedForTrack(esd, vertex, trackpos);
    if(secondLevelAcceptPosTrack==kFALSE) continue;

      p4pos.SetVectM(posTrackMom, daughter2pdgMass);
    
    }
	
      for (Int_t j=0; j<esd->GetNumberOfTracks(); j++) {
        if(iTracks==j) continue;
        AliESDtrack* trackneg=esd->GetTrack(j);
        if (trackneg->GetSign() > 0) continue;
	
        AliExternalTrackParam *tpcTrackneg = (AliExternalTrackParam *)trackneg->GetTPCInnerParam();
        if(!tpcTrackneg) continue;
        ptrackneg[0]=tpcTrackneg->Px();
        ptrackneg[1]=tpcTrackneg->Py();   
        ptrackneg[2]=tpcTrackneg->Pz();  
    
        Bool_t firstLevelAcceptNegTrack=IsAcceptedForKink(esd, vertex, trackneg);
        if(firstLevelAcceptNegTrack==kFALSE) continue;	

        TVector3 negTrackMom(ptrackneg[0],ptrackneg[1],ptrackneg[2]);
	
        TParticle * partneg = stack->Particle(TMath::Abs(trackneg->GetLabel()));
        if (!partneg) continue;
        Int_t pdgneg = partneg->GetPdgCode();
	Int_t mumlabelneg=partneg->GetFirstMother();
        mumlabelneg = TMath::Abs(mumlabelneg);
        TParticle * mumneg=stack->Particle(mumlabelneg);
        if (!mumneg) continue;
        Int_t mumpdgneg = mumneg->GetPdgCode();
	
	Int_t indexKinkNeg=trackneg->GetKinkIndex(0);
	
	if(indexKinkNeg>0) continue;
	
	Bool_t negKaonKinkFlag=0;
	
	if(indexKinkNeg<0) {
	  negKaonKinkFlag=IsKink(esd, indexKinkNeg, negTrackMom);
	
	  if(negKaonKinkFlag==1) p4neg.SetVectM(negTrackMom,daughter1pdgMass);
	  if(negKaonKinkFlag==0) continue;
	}
	
	if(indexKinkNeg==0)  {
 
	   Bool_t secondLevelAcceptNegTrack=IsAcceptedForTrack(esd, vertex, trackneg);
           if(secondLevelAcceptNegTrack==kFALSE) continue;  
	  
	  anp4neg.SetVectM(negTrackMom, daughter2pdgMass);
	
        }
	
	Double_t openingAngle=(ptrackpos[0]*ptrackneg[0]+ptrackpos[1]*ptrackneg[1]+ptrackpos[2]*ptrackneg[2])/(posTrackMom.Mag()*negTrackMom.Mag());

        if((posKaonKinkFlag==1)&&(negKaonKinkFlag==1)) {
	  p4comb=anp4pos;
	  p4comb+=p4neg;
	  if((p4comb.Vect().Pt()>fMinPtTrackCut)&&(TMath::Abs(anp4pos.Vect().Eta())<fmaxAbsEtaCut)&&(TMath::Abs(p4neg.Vect().Eta())<fmaxAbsEtaCut)&&(p4comb.Vect().Eta()<fmaxAbsEtaCut)) {	  
	    if(openingAngle>0.6) fPhiBothKinks->Fill(p4comb.M());
	  }
	}
		
	if(negKaonKinkFlag==1) {
	  p4comb=p4pos;
          p4comb+=p4neg;
	  
	  if(p4comb.Vect().Pt()<=fMinPtTrackCut) continue;
	  
	  if((TMath::Abs(p4pos.Vect().Eta())<fmaxAbsEtaCut)&&(TMath::Abs(p4neg.Vect().Eta())<fmaxAbsEtaCut)&&(p4comb.Vect().Eta()<fmaxAbsEtaCut)) {	  
	    fInvariantMass->Fill(p4comb.M());
	    fInvmassPt->Fill(p4comb.M(), p4comb.Vect().Pt());
	    if ((mumpdgpos==(antiresonancePDGcode))&&(mumpdgneg==(antiresonancePDGcode))&&(mumlabelpos==mumlabelneg)
            &&(pdgpos==fdaughter2pdg)&&(pdgneg==(-fdaughter1pdg))&&(TMath::Abs(trackpos->GetLabel())>=0)&&(TMath::Abs(trackneg->GetLabel())>=0)&&(mumlabelpos>=0)&&(mumlabelneg>=0)) {
              fOpeningAngle->Fill(openingAngle);
              fInvMassTrue->Fill(p4comb.M());
	      fInvmassPtTrue->Fill(p4comb.M(), p4comb.Vect().Pt());

	      fRecPt->Fill(p4comb.Vect().Pt());
	      fRecEta->Fill(p4comb.Vect().Eta());
	      fRecEtaPt->Fill(p4comb.Vect().Perp(),p4comb.Vect().Eta());

	    }

           }
	  
	}
	
	if(posKaonKinkFlag==1) {
          anp4comb=anp4pos;
          anp4comb+=anp4neg;  
	  
	  if(anp4comb.Vect().Pt()<=fMinPtTrackCut) continue;	  
	  
          if((TMath::Abs(anp4neg.Vect().Eta())<fmaxAbsEtaCut)&&(TMath::Abs(anp4pos.Vect().Eta())<fmaxAbsEtaCut)&&(anp4comb.Vect().Eta()<fmaxAbsEtaCut)) {	  
	    fInvariantMass->Fill(anp4comb.M());
	    fInvmassPt->Fill(anp4comb.M(), anp4comb.Vect().Pt());
	    if ((mumpdgpos==resonancePDGcode)&&(mumpdgneg==resonancePDGcode)&&(mumlabelpos==mumlabelneg)
            &&(pdgpos==fdaughter1pdg)&&(pdgneg==(-fdaughter2pdg))&&(TMath::Abs(trackpos->GetLabel())>=0)&&(TMath::Abs(trackneg->GetLabel())>=0)&&(mumlabelpos>=0)  &&(mumlabelneg>=0)) {
              fOpeningAngle->Fill(openingAngle);
              fInvMassTrue->Fill(anp4comb.M());
	      fInvmassPtTrue->Fill(anp4comb.M(), anp4comb.Vect().Pt());
 	
	      fRecPt->Fill(anp4comb.Vect().Pt());
	      fRecEta->Fill(anp4comb.Vect().Eta());
	      fRecEtaPt->Fill(anp4comb.Vect().Pt(), anp4comb.Vect().Eta());
	   }

         }

	}
	 
      } //inner track loop 

  } //outer track loop 
  
  } // end fAnalysisType == ESD
    
}      

//____________________________________________________________________//
void AliResonanceKink::Analyse(AliESDEvent* esd) 
{
  // Main loop
  // Called for each event
  Double_t daughter1pdgMass, daughter2pdgMass;
  
  if (fdaughter1pdg==kKPlus)  {
    daughter1pdgMass=TDatabasePDG::Instance()->GetParticle(fdaughter1pdg)->Mass();
    daughter2pdgMass=TDatabasePDG::Instance()->GetParticle(fdaughter2pdg)->Mass(); 
  }
  
  if (fdaughter1pdg!=kKPlus)  {
    daughter1pdgMass=TDatabasePDG::Instance()->GetParticle(fdaughter2pdg)->Mass();
    daughter2pdgMass=TDatabasePDG::Instance()->GetParticle(fdaughter1pdg)->Mass();   
  }  // to ensure that daughter1pdgMass has always the kaon mass
  
  if (!esd) {
    Printf("ERROR: esd not available");
    return;
  }
 
 if(fAnalysisType == "DATA") {
  const AliESDVertex* vertex = GetEventVertex(esd);
  if(!vertex) return;
  Double_t vtx[3];
  vertex->GetXYZ(vtx);
  fvtxz->Fill(vtx[2]);
  Double_t ptrackpos[3], ptrackneg[3];
  
  TLorentzVector p4pos, anp4pos;
  TLorentzVector p4neg, anp4neg;
  TLorentzVector p4comb, anp4comb;
  
  for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* trackpos = esd->GetTrack(iTracks);
    if (!trackpos) {
      if (fDebug > 0) Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    if (trackpos->GetSign() < 0) continue;
    
    AliExternalTrackParam *tpcTrackpos = (AliExternalTrackParam *)trackpos->GetTPCInnerParam();
    if(!tpcTrackpos) continue;
    ptrackpos[0]=tpcTrackpos->Px();
    ptrackpos[1]=tpcTrackpos->Py();   
    ptrackpos[2]=tpcTrackpos->Pz();  
    
    Bool_t firstLevelAcceptPosTrack=IsAcceptedForKink(esd, vertex, trackpos);
    if(firstLevelAcceptPosTrack==kFALSE) continue;
    
    TVector3 posTrackMom(ptrackpos[0],ptrackpos[1],ptrackpos[2]);
    
    Int_t indexKinkPos=trackpos->GetKinkIndex(0);
    
    if(indexKinkPos>0) continue;
    
    Bool_t posKaonKinkFlag=0;
    
    if(indexKinkPos<0) {
      posKaonKinkFlag=IsKink(esd, indexKinkPos, posTrackMom);
    
      if(posKaonKinkFlag==1) anp4pos.SetVectM(posTrackMom,daughter1pdgMass);
      if(posKaonKinkFlag==0) continue;
    }
    
    if(indexKinkPos==0) {

    Bool_t secondLevelAcceptPosTrack=IsAcceptedForTrack(esd, vertex, trackpos);
    if(secondLevelAcceptPosTrack==kFALSE) continue;

      p4pos.SetVectM(posTrackMom, daughter2pdgMass);
    
    }
	
      for (Int_t j=0; j<esd->GetNumberOfTracks(); j++) {
        if(iTracks==j) continue;
        AliESDtrack* trackneg=esd->GetTrack(j);
        if (trackneg->GetSign() > 0) continue;
	
        AliExternalTrackParam *tpcTrackneg = (AliExternalTrackParam *)trackneg->GetTPCInnerParam();
        if(!tpcTrackneg) continue;
        ptrackneg[0]=tpcTrackneg->Px();
        ptrackneg[1]=tpcTrackneg->Py();   
        ptrackneg[2]=tpcTrackneg->Pz();  
    
        Bool_t firstLevelAcceptNegTrack=IsAcceptedForKink(esd, vertex, trackneg);
        if(firstLevelAcceptNegTrack==kFALSE) continue;	

        TVector3 negTrackMom(ptrackneg[0],ptrackneg[1],ptrackneg[2]);
	
	Int_t indexKinkNeg=trackneg->GetKinkIndex(0);
	
	if(indexKinkNeg>0) continue;
	
	Bool_t negKaonKinkFlag=0;
	
	if(indexKinkNeg<0) {
	  negKaonKinkFlag=IsKink(esd, indexKinkNeg, negTrackMom);
	
	  if(negKaonKinkFlag==1) p4neg.SetVectM(negTrackMom,daughter1pdgMass);
	  if(negKaonKinkFlag==0) continue;
	}
	
	if(indexKinkNeg==0)  {
 
	   Bool_t secondLevelAcceptNegTrack=IsAcceptedForTrack(esd, vertex, trackneg);
           if(secondLevelAcceptNegTrack==kFALSE) continue;  
	  
	  anp4neg.SetVectM(negTrackMom, daughter2pdgMass);
	
        }
	
	Double_t openingAngle=(ptrackpos[0]*ptrackneg[0]+ptrackpos[1]*ptrackneg[1]+ptrackpos[2]*ptrackneg[2])/(posTrackMom.Mag()*negTrackMom.Mag());

        if((posKaonKinkFlag==1)&&(negKaonKinkFlag==1)) {
	  p4comb=anp4pos;
	  p4comb+=p4neg;
	  if((p4comb.Vect().Pt()>fMinPtTrackCut)&&(TMath::Abs(anp4pos.Vect().Eta())<fmaxAbsEtaCut)&&(TMath::Abs(p4neg.Vect().Eta())<fmaxAbsEtaCut)&&(p4comb.Vect().Eta()<fmaxAbsEtaCut)) {	  
	    if(openingAngle>0.6) fPhiBothKinks->Fill(p4comb.M());
	  }
	}
		
	if(negKaonKinkFlag==1) {
	  p4comb=p4pos;
          p4comb+=p4neg;
	  
	  if(p4comb.Vect().Pt()<=fMinPtTrackCut) continue;
	  
	  if((TMath::Abs(p4pos.Vect().Eta())<fmaxAbsEtaCut)&&(TMath::Abs(p4neg.Vect().Eta())<fmaxAbsEtaCut)&&(p4comb.Vect().Eta()<fmaxAbsEtaCut)) {	  
	    fInvariantMass->Fill(p4comb.M());
	    fInvmassPt->Fill(p4comb.M(), p4comb.Vect().Pt());
	    fRecPt->Fill(p4comb.Vect().Pt());
	    fRecEta->Fill(p4comb.Vect().Eta());
	    fRecEtaPt->Fill(p4comb.Vect().Perp(),p4comb.Vect().Eta());
           }
	}
	
	if(posKaonKinkFlag==1) {
          anp4comb=anp4pos;
          anp4comb+=anp4neg;  
	  
	  if(anp4comb.Vect().Pt()<=fMinPtTrackCut) continue;	  
	  
          if((TMath::Abs(anp4neg.Vect().Eta())<fmaxAbsEtaCut)&&(TMath::Abs(anp4pos.Vect().Eta())<fmaxAbsEtaCut)&&(anp4comb.Vect().Eta()<fmaxAbsEtaCut)) {	  
	    fInvariantMass->Fill(anp4comb.M());
	    fInvmassPt->Fill(anp4comb.M(), anp4comb.Vect().Pt());
	    fRecPt->Fill(anp4comb.Vect().Pt());
	    fRecEta->Fill(anp4comb.Vect().Eta());
	    fRecEtaPt->Fill(anp4comb.Vect().Pt(), anp4comb.Vect().Eta()); 
          }

	}
	 
      } //inner track loop 

  } //outer track loop 
  
 } // end fAnalysisType == DATA

}
//____________________________________________________________________//
Float_t AliResonanceKink::GetSigmaToVertex(AliESDtrack* esdTrack) const {
  // Calculates the number of sigma to the vertex.
  
  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];

    esdTrack->GetImpactParametersTPC(b,bCov);
  
  if (bCov[0]<=0 || bCov[2]<=0) {
    bCov[0]=0; bCov[2]=0;
  }
  
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);
  
  if (bRes[0] == 0 || bRes[1] ==0) return -1;
  
  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));
  
  if (TMath::Exp(-d * d / 2) < 1e-10) return 1000;
  
  d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  
  return d;
}

//________________________________________________________________________
const AliESDVertex* AliResonanceKink::GetEventVertex(const AliESDEvent* esd) const
{
  // Get the vertex 
  
  const AliESDVertex* vertex = esd->GetPrimaryVertex();

  if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>2)) return vertex;
  else
  { 
     vertex = esd->GetPrimaryVertexSPD();
     if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>2)) return vertex;
     else
     return 0;
  }
}

//________________________________________________________________________

 Bool_t AliResonanceKink::IsAcceptedForKink(AliESDEvent *localesd,
            const AliESDVertex *localvertex, AliESDtrack* localtrack) {
   // Apply the selections for each kink

  Double_t gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
  Double_t dca3D = 0.0;
  
  AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)localtrack->GetTPCInnerParam();
  if(!tpcTrack) {
    gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0;
    dca[0] = -100.; dca[1] = -100.; dca3D = -100.;
    cov[0] = -100.; cov[1] = -100.; cov[2] = -100.;
  }
  else {
    gPt = tpcTrack->Pt();
    gPx = tpcTrack->Px();
    gPy = tpcTrack->Py();
    gPz = tpcTrack->Pz();
    tpcTrack->PropagateToDCA(localvertex,
    	       localesd->GetMagneticField(),100.,dca,cov);
  }
  
  if(GetSigmaToVertex(localtrack) > fMaxNSigmaToVertex) {
      if (fDebug > 1) Printf("IsAcceptedKink: Track rejected because it has a %lf sigmas to vertex TPC (max. requested: %lf)",   GetSigmaToVertex(localtrack),fMaxNSigmaToVertex);
      return kFALSE;
  }
  
  if(TMath::Abs(dca[0]) > fMaxDCAxy) {
      if (fDebug > 1) Printf("IsAcceptedKink: Track rejected because it has a value of dca(xy) (TPC) of %lf (max. requested: %lf)", TMath::Abs(dca[0]), fMaxDCAxy);
      return kFALSE;
  }
    
  if(TMath::Abs(dca[1]) > fMaxDCAzaxis) {
      if (fDebug > 1) Printf("IsAcceptedKink: Track rejected because it has a value of dca(z) of %lf (max. requested: %lf)", TMath::Abs(dca[1]), fMaxDCAzaxis);
      return kFALSE;
  }
  
  if(gPt <= fMinPtTrackCut) {
      if (fDebug > 1) Printf("IsAcceptedKink: Track rejected because it has a min value of pt of %lf (min. requested: %lf)", gPt, fMinPtTrackCut);
      return kFALSE;
  } 
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliResonanceKink::IsAcceptedForTrack(AliESDEvent *localesd,                                                                                                                                          const AliESDVertex *localvertex, AliESDtrack *localtrack) {
   // Apply the selections for each track

  Double_t gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
  Double_t dca3D = 0.0;
  
  AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)localtrack->GetTPCInnerParam();
  if(!tpcTrack) {
    gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0;
    dca[0] = -100.; dca[1] = -100.; dca3D = -100.;
    cov[0] = -100.; cov[1] = -100.; cov[2] = -100.;
  }
  else {
    gPt = tpcTrack->Pt();
    gPx = tpcTrack->Px();
    gPy = tpcTrack->Py();
    gPz = tpcTrack->Pz();
    tpcTrack->PropagateToDCA(localvertex,
    	       localesd->GetMagneticField(),100.,dca,cov);
  }
  
  Int_t fcls[200];
  Int_t nClustersTPC=localtrack->GetTPCclusters(fcls);
  Float_t chi2perTPCcluster=-1.0;
  if(nClustersTPC!=0) chi2perTPCcluster=(localtrack->GetTPCchi2())/Float_t(nClustersTPC);
  
  Double_t extCov[15];
  localtrack->GetExternalCovariance(extCov);
  
  if((localtrack->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because of no refited in TPC");
      return kFALSE;
  } 

  if(nClustersTPC < fMinTPCclusters) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of nclusters (TPC) of %d (min. requested: %d)", nClustersTPC, fMinTPCclusters);
      return kFALSE;
  } 
  
  if(chi2perTPCcluster > fMaxChi2PerTPCcluster) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of chi2perTPCcluster of %lf (max. requested: %lf)", chi2perTPCcluster, fMaxChi2PerTPCcluster);
      return kFALSE;
  } 

  if(extCov[0] > fMaxCov0) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of cov[0] of %lf (max. requested: %lf)", extCov[0], fMaxCov0);
      return kFALSE;
  }
  
  if(extCov[2] > fMaxCov2) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of cov[2] of %lf (max. requested: %lf)", extCov[2], fMaxCov2);
      return kFALSE;
  }
    
  if(extCov[5] > fMaxCov5) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of cov[5] of %lf (max. requested: %lf)", extCov[5], fMaxCov5);
      return kFALSE;
  }  
  
  if(extCov[9] > fMaxCov9) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of cov[9] of %lf (max. requested: %lf)", extCov[9], fMaxCov9);
      return kFALSE;
  }  
  
  if(extCov[14] > fMaxCov14) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of cov[14] of %lf (max. requested: %lf)", extCov[14], fMaxCov14);
      return kFALSE;
  } 
 
  return kTRUE;

}

//________________________________________________________________________
Bool_t AliResonanceKink::IsKink(AliESDEvent *localesd, Int_t kinkIndex, TVector3 trackMom) 
{
   // Test some kinematical criteria for each kink

	 AliESDkink *kink=localesd->GetKink(TMath::Abs(kinkIndex)-1);
	 const TVector3 motherMfromKink(kink->GetMotherP());
	 const TVector3 daughterMKink(kink->GetDaughterP());
	 Float_t qt=kink->GetQt();

	 Double_t maxDecAngKmu=f1->Eval(motherMfromKink.Mag(),0.,0.,0.);
	 Double_t maxDecAngpimu=f2->Eval(motherMfromKink.Mag(),0.,0.,0.);

         Float_t kinkAngle=TMath::RadToDeg()*kink->GetAngle(2);
	 
	 Float_t energyDaughterMu=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.105658*0.105658);
	 Float_t p1XM= motherMfromKink.Px();
         Float_t p1YM= motherMfromKink.Py();
         Float_t p1ZM= motherMfromKink.Pz();
         Float_t p2XM= daughterMKink.Px();
         Float_t p2YM= daughterMKink.Py();
         Float_t p2ZM= daughterMKink.Pz();
         Float_t p3Daughter=TMath::Sqrt(((p1XM-p2XM)*(p1XM-p2XM))+((p1YM-p2YM)*(p1YM-p2YM))+((p1ZM-p2ZM)*(p1ZM-p2ZM)));
         Double_t invariantMassKmu= TMath::Sqrt((energyDaughterMu+p3Daughter)*(energyDaughterMu+p3Daughter)-motherMfromKink.Mag()*motherMfromKink.Mag());

         if((kinkAngle>maxDecAngpimu)&&(qt>fminQt)&&(qt<fmaxQt)&&((kink->GetR()>fminKinkRadius)&&(kink->GetR()<fmaxKinkRadius))&&(TMath::Abs(trackMom.Eta())<fmaxAbsEtaCut)&&(invariantMassKmu<0.6)) {

           if(trackMom.Mag()<=1.1) {
		return kTRUE;
           }
	   else 
	   if (kinkAngle<maxDecAngKmu) {
		return kTRUE;
	   }
	 }
	 return kFALSE;
}
