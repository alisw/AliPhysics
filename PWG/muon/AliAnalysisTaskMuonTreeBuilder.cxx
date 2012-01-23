#ifndef ALIANALYSISTASKMUONTREEBUILDER_CXX
#define ALIANALYSISTASKMUONTREEBUILDER_CXX

/* $Id$ */

#include "AliAnalysisTaskMuonTreeBuilder.h"
#include "AliMCEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliESDHeader.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1I.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "TChain.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliESDInputHandler.h"
#include "TCanvas.h"
#include "AliPhysicsSelection.h"
#include "AliMultiplicity.h"

//	Analysis task for muon-dimuon analysis
//	Works for real and MC events
//	Works with the corresponding AddTask macro
//	Includes a flag for physics selection
//	
//	author: L. Bianchi - Universita' & INFN Torino

ClassImp(AliAnalysisTaskMuonTreeBuilder)

//__________________________________________________________________________
AliAnalysisTaskMuonTreeBuilder::AliAnalysisTaskMuonTreeBuilder() :
  fNevt(0),
  fBeamEnergy(5000.),
  fOutput(0x0),
  fOutputTree(0x0),
  fIsMC(kFALSE),
  fIsSelected(kFALSE),
  fNumMuonTracks(0),
  fNumSPDTracklets(0),
  fNumContributors(0),
  fNumDimuons(0)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskMuonTreeBuilder::AliAnalysisTaskMuonTreeBuilder(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fNevt(0),
  fBeamEnergy(5000.),
  fOutput(0x0),
  fOutputTree(0x0),
  fIsMC(kFALSE),
  fIsSelected(kFALSE),
  fNumMuonTracks(0),
  fNumSPDTracklets(0),
  fNumContributors(0),
  fNumDimuons(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskMuonTreeBuilder","Calling Constructor");

  DefineOutput(1,TTree::Class());
}

//___________________________________________________________________________
AliAnalysisTaskMuonTreeBuilder& AliAnalysisTaskMuonTreeBuilder::operator=(const AliAnalysisTaskMuonTreeBuilder& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
    fNevt = c.fNevt ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskMuonTreeBuilder::AliAnalysisTaskMuonTreeBuilder(const AliAnalysisTaskMuonTreeBuilder& c) :
  AliAnalysisTaskSE(c),
  fNevt(c.fNevt),
  fBeamEnergy(c.fBeamEnergy),
  fOutput(c.fOutput),
  fOutputTree(c.fOutputTree),
  fIsMC(kFALSE),
  fIsSelected(kFALSE),
  fNumMuonTracks(c.fNumMuonTracks),
  fNumSPDTracklets(c.fNumSPDTracklets),
  fNumContributors(c.fNumContributors),
  fNumDimuons(c.fNumDimuons)
{
  //
  // Copy Constructor										FIDUCIAL REGIONS?
  //
}

//___________________________________________________________________________
AliAnalysisTaskMuonTreeBuilder::~AliAnalysisTaskMuonTreeBuilder() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskMuonTreeBuilder","Calling Destructor");
}

//___________________________________________________________________________
void AliAnalysisTaskMuonTreeBuilder::UserCreateOutputObjects(){
	
//  
//  Creating User-Defined Output Objects
//  
  
 // TREE OUTPUT----------------------------------------------------------
 OpenFile(1);
 fOutputTree = new TTree("krec","Tree of reconstructed muons");

 fOutputTree->Branch("IsSelected",&fIsSelected,"IsSelected/O");
 fOutputTree->Branch("FiredTriggerClasses",fTrigClass,"FiredTriggerClasses/C");

 fOutputTree->Branch("NumMuonTracks",&fNumMuonTracks,"NumMuonTracks/I");
 fOutputTree->Branch("NumContributors",&fNumContributors,"NumContributors/I");
 fOutputTree->Branch("NumSPDTracklets",&fNumSPDTracklets,"NumSPDTraclets/I");
 fOutputTree->Branch("Vertex",fVertex,"Vertex[3]/D");
 fOutputTree->Branch("pT",fpT,"pT[10]/D");
 fOutputTree->Branch("E",fE,"E[10]/D");
 fOutputTree->Branch("px",fpx,"px[10]/D");
 fOutputTree->Branch("py",fpy,"py[10]/D");
 fOutputTree->Branch("pz",fpz,"pz[10]/D");
 fOutputTree->Branch("pxUncorr",fpxUncorr,"pxUncorr[10]/D");
 fOutputTree->Branch("pyUncorr",fpyUncorr,"pyUncorr[10]/D");
 fOutputTree->Branch("pzUncorr",fpzUncorr,"pzUncorr[10]/D");
 fOutputTree->Branch("y",fy,"y[10]/D");
 fOutputTree->Branch("eta",feta,"eta[10]/D");
 fOutputTree->Branch("phi",fphi,"phi[10]/D");
 fOutputTree->Branch("MatchTrig",fMatchTrig,"MatchTrig[10]/I");
 fOutputTree->Branch("TrackChi2",fTrackChi2,"TrackChi2[10]/D");
 fOutputTree->Branch("MatchTrigChi2",fMatchTrigChi2,"MatchTrigChi2[10]/D");
 fOutputTree->Branch("DCA",fDCA,"DCA[10]/D");
 fOutputTree->Branch("Charge",fCharge,"Charge[10]/S");
 fOutputTree->Branch("MuFamily",fMuFamily,"MuFamily[10]/I");
 fOutputTree->Branch("RAtAbsEnd",fRAtAbsEnd,"RAtAbsEnd[10]/D");
 
 fOutputTree->Branch("NumDimuons",&fNumDimuons,"NumDimuons/I");
 fOutputTree->Branch("DimuConstituent",fDimuonConstituent,"DimuonConstituent[45][2]/I");
 fOutputTree->Branch("pTdimuon",fpTdimuon,"pTdimuon[45]/D");
 fOutputTree->Branch("pxdimuon",fpxdimuon,"pxdimuon[45]/D");
 fOutputTree->Branch("pydimuon",fpydimuon,"pydimuon[45]/D");
 fOutputTree->Branch("pzdimuon",fpzdimuon,"pzdimuon[45]/D");
 fOutputTree->Branch("ydimuon",fydimuon,"ydimuon[45]/D");
 fOutputTree->Branch("Imassdimuon",fiMassdimuon,"iMassdimuon[45]/D");
 fOutputTree->Branch("costCS",fcostCS,"costCS[45]/D");
 fOutputTree->Branch("costHE",fcostHE,"costHE[45]/D");
 fOutputTree->Branch("phiCS",fphiCS,"phiCS[45]/D");
 fOutputTree->Branch("phiHE",fphiHE,"phiHE[45]/D");

 fOutputTree->Branch("PDG",fPDG,"PDG[10]/I");
 fOutputTree->Branch("PDGmother",fPDGmother,"PDGmother[10]/I");
 fOutputTree->Branch("PDGdimu",fPDGdimu,"PDGdimu[45]/I");

} 



//_________________________________________________
void AliAnalysisTaskMuonTreeBuilder::UserExec(Option_t *)
{
//  
//  User Exec
//

  fNevt++;


  fNumMuonTracks=0; 
  fNumSPDTracklets=666;
  fNumContributors=666;
  fNumDimuons=0;
  fIsSelected=kFALSE;
  fVertex[0]=666.; fVertex[1]=666.; fVertex[2]=666.;
  for(Int_t i=0; i<10;i++){
    fpT[i]=666.;
    fE[i]=666.;
    fpx[i]=666; 
    fpy[i]=666; 
    fpz[i]=666; 
    fpxUncorr[i]=666; 
    fpyUncorr[i]=666; 
    fpzUncorr[i]=666; 
    fy[i]=666.; 
    feta[i]=666.; 
    fphi[i]=666.;
    fMatchTrig[i]=666; 
    fTrackChi2[i]=666.; 
    fMatchTrigChi2[i]=666.;
    fDCA[i]=666.;
    fPDG[i]=666;
    fPDGmother[i]=666;
    fCharge[i]=666;
    fMuFamily[i]=666;
    fRAtAbsEnd[i]=666;
  }
  for(Int_t i=0; i<45;i++){  
    fDimuonConstituent[i][0]=666;  fDimuonConstituent[i][1]=666;
    fpTdimuon[i]=666.; 
    fpxdimuon[i]=666.; 
    fpydimuon[i]=666.; 
    fpzdimuon[i]=666.; 
    fydimuon[i]=666.; 
    fiMassdimuon[i]=666.;
    fcostCS[i]=666.; 
    fcostHE[i]=666.; 
    fphiCS[i]=666.; 
    fphiHE[i]=666.;
    fPDGdimu[i]=666;
  } 


////////
//// ESD
////////
  
  AliESDEvent *fESD = 0x0;
  AliMCEvent*  mcEvent  = 0x0;

  if(fIsMC){
    if (!fMCEvent) {
      Error("UserExec","NO MC EVENT FOUND!");
      return;
    }
  }

  fESD = dynamic_cast<AliESDEvent*>(InputEvent()); 
  if ( ! fESD ) {
    AliError("Cannot get input event");
    return; 
  }

  fIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  fNumMuonTracks = fESD->GetNumberOfMuonTracks() ;
   Int_t loopEnd = fNumMuonTracks;
  if(!fIsMC) {
    TString cla = fESD->GetFiredTriggerClasses();
    snprintf(fTrigClass,100,"%s",cla.Data());
  }

  if(fNumMuonTracks>0 && fIsMC){
    mcEvent = MCEvent();
  }
  fNumContributors = fESD->GetPrimaryVertexSPD()->GetNContributors();
  fNumSPDTracklets = fESD->GetMultiplicity()->GetNumberOfTracklets();
  fVertex[0]=fESD->GetPrimaryVertexSPD()->GetX();
  fVertex[1]=fESD->GetPrimaryVertexSPD()->GetY();
  fVertex[2]=fESD->GetPrimaryVertexSPD()->GetZ();
  printf("fVertex : %f - %f - %f\n",fVertex[0],fVertex[1],fVertex[2]);

  Int_t numdimu = 0;
  for (Int_t j = 0; j<loopEnd; j++) { 
    AliESDMuonTrack* mu1 = new AliESDMuonTrack(*(fESD->GetMuonTrack(j)));
    if (!mu1->ContainTrackerData()) {fNumMuonTracks=fNumMuonTracks-1; continue;}
    Double_t charge1 = mu1->Charge();
    fCharge[j] = mu1->Charge();
    fpT[j]   = mu1->Pt();
    fpx[j]  = mu1->Px();
    fpy[j]  = mu1->Py();
    fpz[j]  = mu1->Pz();
    fpxUncorr[j]  = mu1->PxUncorrected();
    fpyUncorr[j]  = mu1->PyUncorrected();
    fpzUncorr[j]  = mu1->PzUncorrected();
    fy[j]    = mu1->Y();
    feta[j]    = mu1->Eta();
    fphi[j]  = Phideg(mu1->Phi());
    Double_t emu1    = mu1->E();
    fE[j] = emu1;
    fDCA[j] = mu1->GetDCA();
    fTrackChi2[j] = mu1->GetChi2();
    fMatchTrig[j]   = mu1->GetMatchTrigger();
    fMatchTrigChi2[j]= mu1->GetChi2MatchTrigger();
    fRAtAbsEnd[j]=mu1->GetRAtAbsorberEnd();

    AliMCParticle* mcTrack = 0x0;
    if(fIsMC){
    if(mu1->GetLabel()==-1) continue;
    mcTrack = (AliMCParticle*)mcEvent->GetTrack(mu1->GetLabel());
    fPDG[j] = mcTrack->PdgCode();
    if (mcTrack->GetMother()==-1) continue;
    fPDGmother[j] = ((AliMCParticle*)mcEvent->GetTrack(mcTrack->GetMother()))->PdgCode();
    if (TMath::Abs(fPDG[j])==13) fMuFamily[j] = FindMuFamily(mcTrack,mcEvent);
    }
    for (Int_t jj = j+1; jj<loopEnd; jj++) {
      AliESDMuonTrack* mu2 = new AliESDMuonTrack(*(fESD->GetMuonTrack(jj)));
      if (!mu2->ContainTrackerData()) continue;

      Double_t pxmu2  = mu2->Px();
      Double_t pymu2  = mu2->Py();
      Double_t pzmu2  = mu2->Pz();
      Double_t emu2   = mu2->E();
      //Double_t charge2= mu2->Charge();

      fpTdimuon[numdimu] = TMath::Sqrt((fpx[j]+pxmu2)*(fpx[j]+pxmu2)+(fpy[j]+pymu2)*(fpy[j]+pymu2));
      fpxdimuon[numdimu] = fpx[j]+pxmu2;
      fpydimuon[numdimu] = fpy[j]+pymu2;
      fpzdimuon[numdimu] = fpz[j]+pzmu2;
      fydimuon[numdimu] = Rap((emu1+emu2),(fpz[j]+pzmu2));
      fiMassdimuon[numdimu] = Imass(emu1,fpx[j],fpy[j],fpz[j],emu2,pxmu2,pymu2,pzmu2);
      fcostCS[numdimu]=CostCS(fpx[j],fpy[j],fpz[j],emu1,charge1,pxmu2,pymu2,pzmu2,emu2);
      fcostHE[numdimu]=CostHE(fpx[j],fpy[j],fpz[j],emu1,charge1,pxmu2,pymu2,pzmu2,emu2);
      fphiCS[numdimu] = PhiCS(fpx[j],fpy[j],fpz[j],emu1,charge1,pxmu2,pymu2,pzmu2,emu2);
      fphiHE[numdimu] = PhiHE(fpx[j],fpy[j],fpz[j],emu1,charge1,pxmu2,pymu2,pzmu2,emu2);

      fDimuonConstituent[numdimu][0]=j;  fDimuonConstituent[numdimu][1]=jj;

      numdimu++;
      
      if(fIsMC){
      if(mu2->GetLabel()==-1) continue;
      if(TMath::Abs(fPDG[j])==13 && TMath::Abs(((AliMCParticle*)mcEvent->GetTrack(mu2->GetLabel()))->PdgCode())==13) fPDGdimu[numdimu-1]=FindDimuFamily(mcTrack,(AliMCParticle*)mcEvent->GetTrack(mu2->GetLabel()),mcEvent);
      else fPDGdimu[numdimu-1]=-3;
      }

      delete mu2;
    }
    fNumDimuons=numdimu;
    

    delete mu1;
    }        
  fOutputTree->Fill();
  PostData(1,fOutputTree);
}

//________________________________________________________________________
Int_t AliAnalysisTaskMuonTreeBuilder::FindDimuFamily(AliMCParticle* mcTrack1,AliMCParticle* mcTrack2, AliMCEvent* mcEvent) const
{
  // finds the family of the dimuon (works only if the 2 muons are real muons (not hadrons))
  Int_t familynumber;

  if(mcTrack1->GetMother()==mcTrack2->GetMother()) familynumber = TMath::Abs(((AliMCParticle*)mcEvent->GetTrack(mcTrack1->GetMother()))->PdgCode());
  else{
    Int_t familymu1 = FindMuFamily(mcTrack1,mcEvent);
    Int_t familymu2 = FindMuFamily(mcTrack2,mcEvent);
    if(familymu1==5 && familymu2==5) familynumber=5;							//bb dimuon
    else if(familymu1==4 && familymu2==4) familynumber=4;						//cc dimuon
    else if((familymu1==4 && familymu2==5)||(familymu2==4 && familymu1==5)) familynumber=45;		//bc dimuon
    else if (familymu1==-2 || familymu2==-2 || familymu1==-3 || familymu2==-3) familynumber=-2;		//hadron dimuon (at least 1 hadron involved)
    else familynumber=-1;
  }

  return familynumber;
}

//________________________________________________________________________
Int_t AliAnalysisTaskMuonTreeBuilder::FindMuFamily(AliMCParticle* mcTrack, AliMCEvent* mcEvent) const
{
  // finds the family of the muon
  Int_t imother = mcTrack->GetMother();
  if ( imother<0 ) return -1; // Drell-Yan Muon

  Int_t igrandma = imother;

  AliMCParticle* motherPart = (AliMCParticle*)mcEvent->GetTrack(imother);
  Int_t motherPdg = motherPart->PdgCode();

  // Track is an heavy flavor muon
  Int_t absPdg = TMath::Abs(motherPdg);
  if(absPdg/100==5 || absPdg/1000==5) return 5;
  if(absPdg/100==4 || absPdg/1000==4){
    Int_t newMother = -1;
    igrandma = imother;
    Int_t absGrandMotherPdg = TMath::Abs(motherPart->PdgCode());
    while ( absGrandMotherPdg > 10 ) {
      igrandma = ((AliMCParticle*)mcEvent->GetTrack(igrandma))->GetMother();
      if( igrandma < 0 ) break;
      absGrandMotherPdg = TMath::Abs(((AliMCParticle*)mcEvent->GetTrack(igrandma))->PdgCode());
    }

    if (absGrandMotherPdg==5) newMother = 5; // Charm from beauty
    else if (absGrandMotherPdg==4) newMother = 4;

    if(newMother<0) {
      //AliWarning("Mother not correctly found! Set to charm!\n");
      newMother = 4;
    }

    return newMother;
  }

  Int_t nPrimaries = mcEvent->Stack()->GetNprimary();
  // Track is a bkg. muon
  if (imother<nPrimaries) {
    return -2;			//is a primary
  }
  else {
    return -3;			//is a secondary
  }

}

//________________________________________________________________________
Double_t AliAnalysisTaskMuonTreeBuilder::Imass(Double_t e1, Double_t px1, Double_t py1, Double_t pz1,
				   Double_t e2, Double_t px2, Double_t py2, Double_t pz2) const
{
// invariant mass calculation
    Double_t imassrec = TMath::Sqrt((e1+e2)*(e1+e2)-((px1+px2)*(px1+px2)+
                                    (py1+py2)*(py1+py2)+(pz1+pz2)*(pz1+pz2)));
    return imassrec;
}

//________________________________________________________________________
Double_t AliAnalysisTaskMuonTreeBuilder::Rap(Double_t e, Double_t pz) const
{
// calculate rapidity
    Double_t rap;
    if(e>TMath::Abs(pz)){
	rap = 0.5*TMath::Log((e+pz)/(e-pz));
	return rap;
    }
    else{
	rap = 666.;
	return rap;
    }
}

//________________________________________________________________________
Double_t AliAnalysisTaskMuonTreeBuilder::Phideg(Double_t phi) const
{
// calculate Phi in range [-180,180] 
    Double_t phideg;
    
	phideg = phi-TMath::Pi();
	phideg = phideg*57.296;
	return phideg;
}

//________________________________________________________________________
Double_t AliAnalysisTaskMuonTreeBuilder::CostCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
// Cosine of the theta decay angle (mu+) in the Collins-Soper frame

  TLorentzVector pMu1CM, pMu2CM, pProjCM, pTargCM, pDimuCM; // In the CM. frame
  TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
  TVector3 beta,zaxisCS;
  Double_t mp=0.93827231;
  //
  // --- Fill the Lorentz vector for projectile and target in the CM frame
  //
  pProjCM.SetPxPyPzE(0.,0.,-fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); 
  pTargCM.SetPxPyPzE(0.,0.,fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); 
  //
  // --- Get the muons parameters in the CM frame 
  //
  pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
  pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  pDimuCM=pMu1CM+pMu2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  beta=(-1./pDimuCM.E())*pDimuCM.Vect();
  if(beta.Mag()>=1) return 666.;
  pMu1Dimu=pMu1CM;
  pMu2Dimu=pMu2CM;
  pProjDimu=pProjCM;
  pTargDimu=pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  
  //Debugging part -------------------------------------
  Double_t debugProj[4]={0.,0.,0.,0.};
  Double_t debugTarg[4]={0.,0.,0.,0.};
  Double_t debugMu1[4]={0.,0.,0.,0.};
  Double_t debugMu2[4]={0.,0.,0.,0.};
  pMu1Dimu.GetXYZT(debugMu1);
  pMu2Dimu.GetXYZT(debugMu2);
  pProjDimu.GetXYZT(debugProj);
  pTargDimu.GetXYZT(debugTarg);
  if (debugProj[0]!=debugProj[0] ||debugProj[1]!=debugProj[1] || debugProj[2]!=debugProj[2] ||debugProj[3]!=debugProj[3]) return 666; 
  if (debugTarg[0]!=debugTarg[0] ||debugTarg[1]!=debugTarg[1] || debugTarg[2]!=debugTarg[2] ||debugTarg[3]!=debugTarg[3]) return 666; 
  if (debugMu1[0]!=debugMu1[0] ||debugMu1[1]!=debugMu1[1] || debugMu1[2]!=debugMu1[2] ||debugMu1[3]!=debugMu1[3]) return 666; 
  if (debugMu2[0]!=debugMu2[0] ||debugMu2[1]!=debugMu2[1] || debugMu2[2]!=debugMu2[2] ||debugMu2[3]!=debugMu2[3]) return 666; 
  //----------------------------------------------------

  // --- Determine the z axis for the CS angle 
  zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  				     
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  Double_t cost;
  
  if(charge1>0) {cost = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());}
  else {cost = zaxisCS.Dot((pMu2Dimu.Vect()).Unit());}
  
  return cost;
}

//________________________________________________________________________
Double_t AliAnalysisTaskMuonTreeBuilder::CostHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
// Cosine of the theta decay angle (mu+) in the Helicity frame
  
  TLorentzVector pMu1CM, pMu2CM, pDimuCM; // In the CM frame 
  TLorentzVector pMu1Dimu, pMu2Dimu; // In the dimuon rest frame
  TVector3 beta,zaxisCS;
  //
  // --- Get the muons parameters in the CM frame
  //
  pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
  pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  pDimuCM=pMu1CM+pMu2CM;
  //
  // --- Translate the muon parameters in the dimuon rest frame
  //
  beta=(-1./pDimuCM.E())*pDimuCM.Vect();
  if(beta.Mag()>=1) return 666.;
  pMu1Dimu=pMu1CM;
  pMu2Dimu=pMu2CM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  
  //Debugging part -------------------------------------
  Double_t debugMu1[4]={0.,0.,0.,0.};
  Double_t debugMu2[4]={0.,0.,0.,0.};
  pMu1Dimu.GetXYZT(debugMu1);
  pMu2Dimu.GetXYZT(debugMu2);
  if (debugMu1[0]!=debugMu1[0] ||debugMu1[1]!=debugMu1[1] || debugMu1[2]!=debugMu1[2] ||debugMu1[3]!=debugMu1[3]) return 666; 
  if (debugMu2[0]!=debugMu2[0] ||debugMu2[1]!=debugMu2[1] || debugMu2[2]!=debugMu2[2] ||debugMu2[3]!=debugMu2[3]) return 666; 
  //----------------------------------------------------
 
  // --- Determine the z axis for the calculation of the polarization angle (i.e. the direction of the dimuon in the CM system)
  TVector3 zaxis;
  zaxis=(pDimuCM.Vect()).Unit();
  
  // --- Calculation of the polarization angle (angle between mu+ and the z axis defined above)
  Double_t cost;
  if(charge1>0) {cost = zaxis.Dot((pMu1Dimu.Vect()).Unit());} 
  else {cost = zaxis.Dot((pMu2Dimu.Vect()).Unit());} 
  return cost;
}

//________________________________________________________________________
Double_t AliAnalysisTaskMuonTreeBuilder::PhiCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
// Phi decay angle (mu+) in the Collins-Soper frame

   TLorentzVector pMu1CM, pMu2CM, pProjCM, pTargCM, pDimuCM; // In the CM frame
   TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
   TVector3 beta,yaxisCS, xaxisCS, zaxisCS;
   Double_t mp=0.93827231;
   
   // --- Fill the Lorentz vector for projectile and target in the CM frame
   pProjCM.SetPxPyPzE(0.,0.,-fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); 
   pTargCM.SetPxPyPzE(0.,0.,fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); 
   
   // --- Get the muons parameters in the CM frame 
   pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
   pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
   
   // --- Obtain the dimuon parameters in the CM frame
   pDimuCM=pMu1CM+pMu2CM;
   
   // --- Translate the dimuon parameters in the dimuon rest frame
   beta=(-1./pDimuCM.E())*pDimuCM.Vect();
   if(beta.Mag()>=1) return 666.;
   pMu1Dimu=pMu1CM;
   pMu2Dimu=pMu2CM;
   pProjDimu=pProjCM;
   pTargDimu=pTargCM;
   pMu1Dimu.Boost(beta);
   pMu2Dimu.Boost(beta);
   pProjDimu.Boost(beta);
   pTargDimu.Boost(beta);

   //Debugging part -------------------------------------
   Double_t debugProj[4]={0.,0.,0.,0.};
   Double_t debugTarg[4]={0.,0.,0.,0.};
   Double_t debugMu1[4]={0.,0.,0.,0.};
   Double_t debugMu2[4]={0.,0.,0.,0.};
   pMu1Dimu.GetXYZT(debugMu1);
   pMu2Dimu.GetXYZT(debugMu2);
   pProjDimu.GetXYZT(debugProj);
   pTargDimu.GetXYZT(debugTarg);
   if (debugProj[0]!=debugProj[0] ||debugProj[1]!=debugProj[1] || debugProj[2]!=debugProj[2] ||debugProj[3]!=debugProj[3]) return 666; 
   if (debugTarg[0]!=debugTarg[0] ||debugTarg[1]!=debugTarg[1] || debugTarg[2]!=debugTarg[2] ||debugTarg[3]!=debugTarg[3]) return 666; 
   if (debugMu1[0]!=debugMu1[0] ||debugMu1[1]!=debugMu1[1] || debugMu1[2]!=debugMu1[2] ||debugMu1[3]!=debugMu1[3]) return 666; 
   if (debugMu2[0]!=debugMu2[0] ||debugMu2[1]!=debugMu2[1] || debugMu2[2]!=debugMu2[2] ||debugMu2[3]!=debugMu2[3]) return 666; 
   //----------------------------------------------------

   // --- Determine the z axis for the CS angle 
   zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
   yaxisCS=(((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
   xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit();
 
   Double_t phi=0.;
   if(charge1>0) {
       phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
   } else {
       phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxisCS),((pMu2Dimu.Vect()).Dot(xaxisCS)));
   }
   if (phi>TMath::Pi()) phi=phi-TMath::Pi();
   
   return phi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskMuonTreeBuilder::PhiHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
// Phi decay angle (mu+) in the Helicity frame
  TLorentzVector pMu1Lab, pMu2Lab, pProjLab, pTargLab, pDimuLab; // In the lab. frame 
  TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
  TVector3 beta,xaxis, yaxis,zaxis;
  Double_t mp=0.93827231;

  // --- Get the muons parameters in the LAB frame
  pMu1Lab.SetPxPyPzE(px1,py1,pz1,e1);
  pMu2Lab.SetPxPyPzE(px2,py2,pz2,e2);
  
  // --- Obtain the dimuon parameters in the LAB frame
  pDimuLab=pMu1Lab+pMu2Lab;
  zaxis=(pDimuLab.Vect()).Unit();
  
  // --- Translate the muon parameters in the dimuon rest frame
  beta=(-1./pDimuLab.E())*pDimuLab.Vect();
  if(beta.Mag()>=1.) return 666.;

  pProjLab.SetPxPyPzE(0.,0.,-fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); // proiettile
  pTargLab.SetPxPyPzE(0.,0.,fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); // bersaglio

  pProjDimu=pProjLab;
  pTargDimu=pTargLab;

  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  
  yaxis=((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
  xaxis=(yaxis.Cross(zaxis)).Unit();
  
  pMu1Dimu=pMu1Lab;
  pMu2Dimu=pMu2Lab;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  
  //Debugging part -------------------------------------
  Double_t debugProj[4]={0.,0.,0.,0.};
  Double_t debugTarg[4]={0.,0.,0.,0.};
  Double_t debugMu1[4]={0.,0.,0.,0.};
  Double_t debugMu2[4]={0.,0.,0.,0.};
  pMu1Dimu.GetXYZT(debugMu1);
  pMu2Dimu.GetXYZT(debugMu2);
  pProjDimu.GetXYZT(debugProj);
  pTargDimu.GetXYZT(debugTarg);
  if (debugProj[0]!=debugProj[0] ||debugProj[1]!=debugProj[1] || debugProj[2]!=debugProj[2] ||debugProj[3]!=debugProj[3]) return 666; 
  if (debugTarg[0]!=debugTarg[0] ||debugTarg[1]!=debugTarg[1] || debugTarg[2]!=debugTarg[2] ||debugTarg[3]!=debugTarg[3]) return 666; 
  if (debugMu1[0]!=debugMu1[0] ||debugMu1[1]!=debugMu1[1] || debugMu1[2]!=debugMu1[2] ||debugMu1[3]!=debugMu1[3]) return 666; 
  if (debugMu2[0]!=debugMu2[0] ||debugMu2[1]!=debugMu2[1] || debugMu2[2]!=debugMu2[2] ||debugMu2[3]!=debugMu2[3]) return 666; 
  //----------------------------------------------------
  
  Double_t phi=0.;
   if(charge1 > 0) {
      phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
     } else { 
      phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxis),(pMu2Dimu.Vect()).Dot(xaxis));
   }  
   return phi;
}

//________________________________________________________________________
void AliAnalysisTaskMuonTreeBuilder::Terminate(Option_t *) 
{
// Terminate
}

#endif
