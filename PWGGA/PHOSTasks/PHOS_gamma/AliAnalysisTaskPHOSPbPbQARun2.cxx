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
/* $Id$ */

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskPHOSPbPbQARun2.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliCentrality.h"
#include "AliLog.h"
#include "TObjArray.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "AliMultSelection.h"
#include "AliAODMCParticle.h"

// Stripped-down version of Dmitri Peressounko' AliAnalysisTaskPi0Flow class
// used for the fast QA of PbPb data.
//...
// Author: Boris Polishchuk 
// Date   : 19.10.2011

ClassImp(AliAnalysisTaskPHOSPbPbQARun2)

//________________________________________________________________________
AliAnalysisTaskPHOSPbPbQARun2::AliAnalysisTaskPHOSPbPbQARun2() : AliAnalysisTaskSE(),
  fOutputContainer(0),fPHOSEvent(0),fCentrality(0),fCenBin(0),
  fPHOSGeo(0),fEventCounter(0), fMCArray(0)
{
  //Default constructor

  for(Int_t i=0; i<1;i++){
    for(Int_t j=0; j<5;j++)
      fPHOSEvents[i][j]=0 ;
  }

  // Initialize the PHOS geometry 
  fPHOSGeo = AliPHOSGeometry::GetInstance("Run2") ;

}

//________________________________________________________________________
AliAnalysisTaskPHOSPbPbQARun2::AliAnalysisTaskPHOSPbPbQARun2(const char *name) 
: AliAnalysisTaskSE(name),
  fOutputContainer(0),
  fPHOSEvent(0),
  fCentrality(0),fCenBin(0),
  fPHOSGeo(0),
  fEventCounter(0),
  fMCArray(0)
{
  // Constructor
  for(Int_t i=0; i<1; i++) {
    for(Int_t j=0; j<5; j++)
	fPHOSEvents[i][j]=0 ;
  }

  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());

  // Initialize the PHOS geometry
  fPHOSGeo = AliPHOSGeometry::GetInstance("Run2") ;

}

//________________________________________________________________________
void AliAnalysisTaskPHOSPbPbQARun2::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // AOD histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }

  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);

  fOutputContainer->Add(new TH1I("hSelEvent", "Event selection", 10, 0, 10));

  fOutputContainer->Add(new TH1I("hEvCenBins", "events per centrality bin", 5, 0, 5));

  fOutputContainer->Add(new TH2F("hCenPHOS","Centrality vs PHOSclusters", 100,0.,100.,200,0.,200.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSCells","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH2F("hCenTrack","Centrality vs tracks", 100,0.,100.,100,0.,15000.)) ;  


  fOutputContainer->Add(new TH1F("hZvertex","vertex z coordinate",400, -20.,+20.));
  fOutputContainer->Add(new TH1F("hZvertexSPD","SPD vertex z coordinate",400, -20.,+20.));
  fOutputContainer->Add(new TH1F("hDistanceVSPD","distance between vertices",400, -20.,+20.));
  fOutputContainer->Add(new TH1F("hNContributors","N of primary tracks from the primary vertex", 1e4, 0., 1e4));
 
  fOutputContainer->Add(new TH1I("hNPileupVtx", "Number of pileup vertices", 10, 0, 10));
  fOutputContainer->Add(new TH1F("hZPileupVtx", "Location of pileup vertices", 400, -20, 20));

  //pi0 spectrum
  Int_t nPtPhot = 300 ;
  Double_t ptPhotMax = 30 ;
  Int_t nM       = 500;
  Double_t mMin  = 0.0;
  Double_t mMax  = 1.0;
  char key[55] ;

  for(Int_t cent=0; cent<5; cent++){
    snprintf(key,55, "hClusterE_cen%d", cent);
    fOutputContainer->Add(new TH1F(key, "cluster energy", nPtPhot, 0., ptPhotMax));

    snprintf(key,55, "hClusterPt_cen%d", cent);
    fOutputContainer->Add(new TH1F(key, "cluster p_{T}", nPtPhot, 0., ptPhotMax));

    snprintf(key,55, "hNcellvsEcl_cen%d", cent);
    fOutputContainer->Add(new TH2F(key, "ncells vs cluster energy", nPtPhot, 0., ptPhotMax, 40, 0, 40));

    snprintf(key,55,"hPi0All_cen%d",cent) ;
    TString  mhTitle =  "Invariant mass, All clusters;M_{#gamma#gamma};p_{T}";
    fOutputContainer->Add(new TH2F(key, mhTitle , nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

    snprintf(key,55,"hMiPi0All_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key, mhTitle ,nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

    snprintf(key,55,"hPi0SingleAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,mhTitle ,nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

    snprintf(key,55,"hMiPi0SingleAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key, mhTitle, nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    
    snprintf(key,55,"hGammaMC_cen%d", cent); 
    fOutputContainer->Add(new  TH2F(key, "MC #gamma;p_{T};y", nPtPhot, 0, ptPhotMax, 240, -1.2, 1.2));
    
    for (Int_t sm = 1; sm < 5; sm ++) {

       snprintf(key,55, "hClusterEM%d_cen%d", sm, cent);
       fOutputContainer->Add(new TH1F(key, "cluster energy", nPtPhot, 0., ptPhotMax));

       snprintf(key,55, "hClusterPtM%d_cen%d", sm, cent);
       fOutputContainer->Add(new TH1F(key, "cluster p_{T}", nPtPhot, 0., ptPhotMax));

       snprintf(key,55, "hNcellvsEclM%d_cen%d", sm, cent);
       fOutputContainer->Add(new TH2F(key, "ncells vs cluster energy", nPtPhot, 0., ptPhotMax, 40, 0, 40));

       snprintf(key,55,"hPi0AllSM%d_cen%d",sm,cent) ;
       fOutputContainer->Add(new TH2F(key, mhTitle ,nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

       snprintf(key,55,"hMiPi0AllSM%d_cen%d",sm,cent) ;
       fOutputContainer->Add(new TH2F(key,mhTitle,nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

       snprintf(key,55,"hPi0SingleAllSM%d_cen%d",sm,cent) ;
       fOutputContainer->Add(new TH2F(key, mhTitle, nM, mMin, mMax, nPtPhot, 0., ptPhotMax));

       snprintf(key,55,"hMiPi0SingleAllSM%d_cen%d",sm,cent) ;
       fOutputContainer->Add(new TH2F(key, mhTitle, nM, mMin, mMax, nPtPhot, 0., ptPhotMax));
    }

  }

  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisTaskPHOSPbPbQARun2::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze AOD/AOD  
  //

  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());

  if (!event) {
    Printf("ERROR: Could not retrieve event");
    return;
  }

  FillHistogram("hSelEvent", 0.5);

  const AliAODVertex *aodVertex =  event->GetPrimaryVertex();
  const AliAODVertex *aodVertexSPD  = event->GetPrimaryVertexSPD();

  FillHistogram("hSelEvent", 1.5);

  if (!aodVertex) return;

  FillHistogram("hSelEvent", 2.5);

  const Int_t ncont = aodVertex->GetNContributors();
  FillHistogram("hNContributors", ncont);
  if (ncont < 1) return;

  FillHistogram("hSelEvent", 3.5);

  const Double_t vtx[3] = {aodVertex->GetX(), aodVertex->GetY(), aodVertex->GetZ()}; // vertex coordinated
  FillHistogram("hZvertex", vtx[2]);
  FillHistogram("hZvertexSPD", aodVertexSPD->GetZ());
  FillHistogram("hDistanceVSPD", vtx[2] - aodVertexSPD->GetZ());
  if (TMath::Abs(vtx[2]) > 10.) return;

  FillHistogram("hSelEvent", 4.5);

  if (event->IsPileupFromSPD()) {
    FillHistogram("hNPileupVtx", event->GetNumberOfPileupVerticesSPD());
    for (Int_t puVtx = 0;  puVtx < event->GetNumberOfPileupVerticesSPD(); puVtx++) {
      Double_t dZpileup = aodVertexSPD->GetZ() - event->GetPileupVertexSPD(puVtx)->GetZ();
      FillHistogram("hZPileupVtx", dZpileup);
    }
    return;
  }

  FillHistogram("hSelEvent", 5.5);

  char key[55] ;  

  Int_t zvtx=0 ;

  AliMultSelection *multSelection =static_cast<AliMultSelection*>(event->FindListObject("MultSelection"));
  if(multSelection) fCentrality = multSelection->GetMultiplicityPercentile("V0M");

  if(fCentrality < 20.)  fCenBin = 0;
  if(fCentrality >= 20.) fCenBin = 1;
  if(fCentrality >= 40.) fCenBin = 2;
  if(fCentrality >= 60.) fCenBin = 3;
  if(fCentrality >= 80.) fCenBin = 4;

  FillHistogram("hEvCenBins", fCenBin + 0.5);

  //printf("event nr %d, centrality %.3f [%d]\n", fEventCounter, fCentrality, fCenBin);

  if(!fPHOSEvents[zvtx][fCenBin]) 
    fPHOSEvents[zvtx][fCenBin]=new TList() ;

  TList * prevPHOS = fPHOSEvents[zvtx][fCenBin] ;

  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton", 200) ;

  Int_t multClust = event->GetNumberOfCaloClusters();
  AliAODCaloCells * cells = event->GetPHOSCells() ;

  FillHistogram("hCenPHOSCells",fCentrality,cells->GetNumberOfCells()) ;
  FillHistogram("hCenTrack",fCentrality,event->GetNumberOfTracks()) ;

  Int_t inPHOS = 0;

  for (Int_t i=0; i<multClust; i++) {

    AliAODCaloCluster *clu = event->GetCaloCluster(i);

    if (!clu->IsPHOS() || clu->GetCoreEnergy() < 0.3) continue;
    if (clu->GetNCells() < 3) continue;

    if (clu->GetType() != AliVCluster::kPHOSNeutral) continue; //select calorimeter clusters

    fMCArray = (TClonesArray*)event->FindListObject(AliAODMCParticle::StdBranchName());
    if (!fMCArray && TMath::Abs(clu->GetTOF()) > 12.5e-9) continue; // TOF cut for real data only!

    Float_t  position[3];
    clu->GetPosition(position);
    TVector3 global(position) ;
    Int_t relId[4] ;
    fPHOSGeo->GlobalPos2RelId(global,relId) ;
    Int_t mod  =  relId[0] ;
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;   

    if (mod < 1 || mod > 4) {
      Printf("Wrong module number %d", mod) ;
      continue ;
    }  

    if (clu->GetM02() < 0.1) continue;    

    FillHistogram(Form("hClusterE_cen%d", fCenBin), clu->E());
    FillHistogram(Form("hClusterEM%d_cen%d", mod, fCenBin), clu->E());
    FillHistogram(Form("hNcellvsEcl_cen%d", fCenBin), clu->E(), clu->GetNCells()); 
    FillHistogram(Form("hNcellvsEclM%d_cen%d", mod, fCenBin), clu->E(), clu->GetNCells()); 

    //ratio of core energy to total cluster energy
    Double_t ecore = clu->GetCoreEnergy();
    Double_t efull = clu->E();
    Double_t r     = ecore/efull;

    TLorentzVector pv1 ;
    clu->GetMomentum(pv1 ,vtx);
    pv1 = pv1*r;

    if (inPHOS>=fPHOSEvent->GetSize()){
      fPHOSEvent->Expand(inPHOS+50) ;
    }

    new((*fPHOSEvent)[inPHOS]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod) ;
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());

    ph->SetEMCx(float(cellX)) ;
    ph->SetEMCz(float(cellZ)) ;

    inPHOS++ ;
  }

  FillHistogram("hCenPHOS",fCentrality, inPHOS) ;

  //photons
  for (Int_t i=0; i<inPHOS-1; i++) {
    AliCaloPhoton * ph=(AliCaloPhoton*)fPHOSEvent->At(i) ;
    Int_t sm = ph->Module();
    FillHistogram(Form("hClusterPt_cen%d", fCenBin), ph->Pt());
    FillHistogram(Form("hClusterPtM%d_cen%d", sm, fCenBin), ph->Pt());
  }

  //pi0
  for (Int_t i1=0; i1<inPHOS-1; i1++) {

    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    Int_t sm1 = ph1->Module();

    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;

      Int_t sm2 = ph2->Module();      
      TLorentzVector p12  = *ph1  + *ph2;

      snprintf(key,55,"hPi0All_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt()) ; 

      snprintf(key,55,"hPi0SingleAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,ph1->Pt()) ; 
      FillHistogram(key,p12.M() ,ph2->Pt()) ; 

      if (sm1==sm2) {
	snprintf(key,55,"hPi0AllSM%d_cen%d",sm1,fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ; 

	snprintf(key,55,"hPi0SingleAllSM%d_cen%d",sm1,fCenBin) ;
        FillHistogram(key,p12.M() ,ph1->Pt()) ; 
        FillHistogram(key,p12.M() ,ph2->Pt()) ; 
      }

    } // end of loop i2
  } // end of loop i1

  //now mixed
  for (Int_t i1=0; i1 <inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    Int_t sm1 = ph1->Module();

    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;

      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
	AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;

	Int_t sm2 = ph2->Module();
	TLorentzVector p12  = *ph1  + *ph2;

	snprintf(key,55,"hMiPi0All_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;

	snprintf(key,55,"hMiPi0SingleAll_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,ph1->Pt()) ;
	FillHistogram(key,p12.M() ,ph2->Pt()) ;

	if (sm1==sm2) {
	  snprintf(key,55,"hMiPi0AllSM%d_cen%d",sm1,fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ; 

	  snprintf(key,55,"hMiPi0SingleAllSM%d_cen%d",sm1,fCenBin) ;
	  FillHistogram(key,p12.M() ,ph1->Pt()) ;
	  FillHistogram(key,p12.M() ,ph2->Pt()) ;
	}
      } // end of loop i2
    }
  } // end of loop i1


  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed
  if(fPHOSEvent->GetEntriesFast()>0){
    prevPHOS->AddFirst(fPHOSEvent) ;
    fPHOSEvent=0;
    if(prevPHOS->GetSize()>100){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last()) ;
      prevPHOS->RemoveLast() ;
      delete tmp ;
    }
  }

  if (fMCArray) {
    for (Int_t i = 0; i < fMCArray->GetEntriesFast(); i++) {
      auto particle =  (AliAODMCParticle*) fMCArray->At(i);
      Int_t pdg = particle->GetPdgCode();
      Double_t xv = particle->Xv(), yv = particle->Yv();
      Double_t rv = TMath::Hypot(xv, yv);
      if (particle->IsSecondaryFromMaterial()) continue;
      if (pdg == 22 && rv < 1.0 && particle->GetLabel() > -1) {
         Int_t iMother = particle->GetMother();
	 if (iMother > -1) {
           auto mparticle = (AliAODMCParticle*)fMCArray->At(iMother);
	   Int_t pdgm = mparticle->GetPdgCode();
	   if (pdgm == 111 && mparticle->GetLabel() > -1) {
             Int_t iMother2 = mparticle->GetMother();
	     if (iMother2 > -1) {
	       auto mparticle2 = (AliAODMCParticle*)fMCArray->At(iMother2);
	       Int_t pdg2 = mparticle2->GetPdgCode();
	       if (pdg2 == 310 || pdg2 == 130) continue;
	     }  
	   }
         }
	 FillHistogram(Form("hGammaMC_cen%d", fCenBin), particle->Pt(), particle->Y());
      }
    }
  }
  
  fEventCounter++;
}

//_____________________________________________________________________________
void AliAnalysisTaskPHOSPbPbQARun2::FillHistogram(const char * key,Double_t x)const
{
  //FillHistogram
  TH1I * tmpI = dynamic_cast<TH1I*>(fOutputContainer->FindObject(key)) ;
  if(tmpI){
    tmpI->Fill(x) ;
    return ;
  }
  TH1F * tmpF = dynamic_cast<TH1F*>(fOutputContainer->FindObject(key)) ;
  if(tmpF){
    tmpF->Fill(x) ;
    return ;
  }
  TH1D * tmpD = dynamic_cast<TH1D*>(fOutputContainer->FindObject(key)) ;
  if(tmpD){
    tmpD->Fill(x) ;
    return ;
  }
  AliInfo(Form("can not find histogram <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskPHOSPbPbQARun2::FillHistogram(const char * key,Double_t x,Double_t y)const{
  //FillHistogram
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1F")){
    ((TH1F*)tmp)->Fill(x,y) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y) ;
    return ;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName())) ;
}
