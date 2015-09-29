#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"

#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAnalysisTaskNeutralMesons.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "TGeoManager.h"

#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliStack.h"
#include "TParticle.h"

#include "AliOADBContainer.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

ClassImp(AliAnalysisTaskNeutralMesons)

//
//________________________________________________________________________
//
AliAnalysisTaskNeutralMesons::AliAnalysisTaskNeutralMesons()
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskNeutralMesons")
  , fESD(0)
  , fMC(0)
  , fGeom(0x0)
  , fRecoUtil(0x0)
  , fOutputList(0x0)
  , fRejectExoticCluster(kTRUE)
  , fRemoveBadChannel(kTRUE)
  , fhEclus(0x0)
  , fhEclusSM0(0x0)
  , fhEclusSM1(0x0)
  , fhEclusSM2(0x0)
  , fhEclusSM3(0x0)
  , fhEclusSM4(0x0)
  , fhEclusSM5(0x0)
  , fhEclusSM6(0x0)
  , fhEclusSM7(0x0)
  , fhEclusSM8(0x0)
  , fhEclusSM9(0x0)
  , fhEVsTime(0x0)
  , fhEVsNcells(0x0)
  , fhZvertex(0x0)
  , f(0x0)
  , f1(0x0)
  , f2(0x0)
  , f3(0x0)
  , fPi0(0x0)
  , fEta(0x0)
  , fhZvertexAll(0x0)
  , fclusterList(0)
  , fparticle(0)
  , fhpi0Decay(0)
  , fhpi0Inter(0)
  , fhRecDecay(0)
  , fhRecInter(0)
  , fPrim(0)
  , fPrimNoTRD(0)
  , fK0(0)
  , fPi0mixed(0x0)
  , fEtamixed(0x0)
  , fPool(0)
  , fhHitMap(0)
  , fNcl(0)
  , fNeventCls(0)
  , fTrpi0(0)
  , fK0Pt(0)
  , fhLowm(0)
  , fhHighm(0)
  , fPi0Conv(0x0)
  , fhLowmConv(0)
  , fhHighmConv(0)
  , fhPhotonEtaPhi(0)
  , fMCtype(1)
  , fMCpart(111)
  , fPrimOldGeo(0)
  , fGeoName(" ")
  , fLowEnergCut(0.3)
  , fNcellsInClusterFlag(kFALSE)
  , fBranchName(" ")
{
  
	// Constructor	
  fClusterZvtx=new TArrayD(1500);
  fNEvClusters=new TArrayD(1500);
  DefineInput (0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//
//________________________________________________________________________
//
AliAnalysisTaskNeutralMesons::AliAnalysisTaskNeutralMesons(const char *name)
  : AliAnalysisTaskSE(name)
  , fESD(0)
  , fMC(0)
  , fGeom(0x0)
  , fRecoUtil(0x0)
  , fRejectExoticCluster(kTRUE)
  , fRemoveBadChannel(kFALSE)
  , fhEclus(0x0)
  , fhEclusSM0(0x0)
  , fhEclusSM1(0x0)
  , fhEclusSM2(0x0)
  , fhEclusSM3(0x0)
  , fhEclusSM4(0x0)
  , fhEclusSM5(0x0)
  , fhEclusSM6(0x0)
  , fhEclusSM7(0x0)
  , fhEclusSM8(0x0)
  , fhEclusSM9(0x0)
  , fhEVsTime(0x0)
  , fhEVsNcells(0x0)
  , fhZvertex(0x0)
  , f(0x0)
  , f1(0x0)
  , f2(0x0)
  , f3(0x0)
  , fPi0(0x0)
  , fEta(0x0)
  , fhZvertexAll(0x0)
  , fclusterList(0)
  , fparticle(0)
  , fhpi0Decay(0)
  , fhpi0Inter(0)
  , fhRecDecay(0)
  , fhRecInter(0)
  , fPrim(0)
  , fPrimNoTRD(0)
  , fK0(0)
  , fPi0mixed(0x0)
  , fEtamixed(0x0)
  , fPool(0)
  , fhHitMap(0)
  , fNcl(0)
  , fNeventCls(0)
  , fTrpi0(0)
  , fK0Pt(0)
  , fhLowm(0)
  , fhHighm(0)
  , fPi0Conv(0x0)
  , fhLowmConv(0)
  , fhHighmConv(0)
  , fhPhotonEtaPhi(0)
  , fMCtype(1)
  , fMCpart(111)
  , fPrimOldGeo(0)
  , fGeoName(" ")
  , fLowEnergCut(0.3)
  , fNcellsInClusterFlag(kFALSE)
  , fBranchName(" ")
 {
  
	// Constructor	
  //DefineInput (0, TChain::Class());
  fClusterZvtx=new TArrayD(1500);
  fNEvClusters=new TArrayD(1500);
  DefineOutput(1, TList::Class());
}

//
//________________________________________________________________________
//
AliAnalysisTaskNeutralMesons::~AliAnalysisTaskNeutralMesons()
{
  //Destructor

  if(fRecoUtil) delete fRecoUtil;
  if(fOutputList) delete fOutputList;
  if(fGeom) delete fGeom;

}

//
//________________________________________________________________________
//
void AliAnalysisTaskNeutralMesons::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  //AliInfo("USER CREATE OUTPUT OBJECTS");
  
  fOutputList = new TList();
  //fOutputList->SetOwner(1);
  
  fhEclus = new TH1D("fhEclus", "cluster energy", 300, 0., 30.);
  fhEclusSM0 = new TH1D("fhEclusSM0", "cluster energy 0", 300, 0., 30.);
  fhEclusSM1 = new TH1D("fhEclusSM1", "cluster energy 1", 300, 0., 30.);
  fhEclusSM2 = new TH1D("fhEclusSM2", "cluster energy 2", 300, 0., 30.);
  fhEclusSM3 = new TH1D("fhEclusSM3", "cluster energy 3", 300, 0., 30.);
  fhEclusSM4 = new TH1D("fhEclusSM4", "cluster energy 4", 300, 0., 30.);
  fhEclusSM5 = new TH1D("fhEclusSM5", "cluster energy 5", 300, 0., 30.);
  fhEclusSM6 = new TH1D("fhEclusSM6", "cluster energy 6", 300, 0., 30.);
  fhEclusSM7 = new TH1D("fhEclusSM7", "cluster energy 7", 300, 0., 30.);
  fhEclusSM8 = new TH1D("fhEclusSM8", "cluster energy 8", 300, 0., 30.);
  fhEclusSM9 = new TH1D("fhEclusSM9", "cluster energy 9", 300, 0., 30.);

  fhEVsTime = new TH2D("fhEVsTime", ";cell Energy (GeV);", 400, 0., 20., 1200, -300., 300.);
  fhEVsNcells = new TH2D("fhEVsNcells", ";cell Energy (GeV);", 500, 0., 50., 200, -0.5, 199.5);
  fhZvertex = new TH1D("fhZvertex", "Vertex Z distribution", 140, -35., 35.);
  fPi0 = new TH2D("fPi0",";M_{#gamma#gamma} [GeV]; pt_{#gamma#gamma} [GeV/c]",1000,0.0,1.0,1500,0.0,30.0);
  fEta = new TH2D("fEta",";M_{#gamma#gamma} [GeV]; pt_{#gamma#gamma} [GeV/c]",1000,0.0,1.0,1500,0.0,30.0);
  fhZvertexAll = new TH1D("fhZvertexAll", "Vertex Z distribution", 140, -35., 35.);
  fhpi0Decay = new TH2D("fhpi0Decay","  ",1000,0.0,500.0,120,-30.0,30.0);
  fhpi0Inter = new TH2D("fhpi0Inter","  ",1000,0.0,500.0,120,-30.0,30.0);
  fhRecDecay = new TH2D("fhRecDecay",";M_{#gamma#gamma} [GeV]; pt_{#gamma#gamma} [GeV/c]",1000,0.0,1.0,1500,0.0,30.0);
  fhRecInter = new TH2D("fhRecInter",";M_{#gamma#gamma} [GeV]; pt_{#gamma#gamma} [GeV/c]",1000,0.0,1.0,1500,0.0,30.0);
  fPrim = new TH1D("fPrim", "  ", 2000,0.0,40.0);
  fPrimNoTRD = new TH1D("fPrimNoTRD", "  ", 2000,0.0,40.0);

  fK0 = new TH1D("fK0", "  ", 2000,0.0,40.0);
  fPi0mixed = new TH2D("fPi0mixed",";M_{#gamma#gamma} [GeV]; pt_{#gamma#gamma} [GeV/c]",1000,0.0,1.0,1500,0.0,30.0);
  fEtamixed = new TH2D("fEtamixed",";M_{#gamma#gamma} [GeV]; pt_{#gamma#gamma} [GeV/c]",1000,0.0,1.0,1500,0.0,30.0); 
  fhHitMap = new TH2D("fhHitMap", "Vertex Z distribution out of 10", 140, -0.7, 0.7, 1000, 80., 180.);
  fNcl = new TH1D("fNcl", "clusters in the event", 20, -0.5, 19.5);
  fK0Pt = new TH1D("fK0Pt", "  ", 2000,0.0,40.0);
  fhLowm = new TH2D("fhLowm",";M_{#gamma#gamma} [GeV]; pt_{#gamma#gamma} [GeV/c]",1000,0.0,1.0,1500,0.0,30.0);
  fhHighm = new TH2D("fhHighm",";M_{#gamma#gamma} [GeV]; pt_{#gamma#gamma} [GeV/c]",1000,0.0,1.0,1500,0.0,30.0);
  fPi0Conv = new TH2D("fPi0Conv",";M_{#gamma#gamma} [GeV]; pt_{#gamma#gamma} [GeV/c]",1000,0.0,1.0,1500,0.0,30.0);
  fhLowmConv = new TH2D("fhLowmConv",";M_{#gamma#gamma} [GeV]; pt_{#gamma#gamma} [GeV/c]",1000,0.0,1.0,1500,0.0,30.0);
  fhHighmConv = new TH2D("fhHighmConv",";M_{#gamma#gamma} [GeV]; pt_{#gamma#gamma} [GeV/c]",1000,0.0,1.0,1500,0.0,30.0);			      
  fhPhotonEtaPhi = new TH2D("fhPhotonEtaPhi", "  ", 140, -0.7, 0.7, 200, 80., 180.);
  fPrimOldGeo = new TH1D("fPrimOldGeo", "  ", 2000,0.0,40.0);
   
  fOutputList->Add(fhEclus);
  fOutputList->Add(fhEclusSM0);
  fOutputList->Add(fhEclusSM1);
  fOutputList->Add(fhEclusSM2);
  fOutputList->Add(fhEclusSM3);
  fOutputList->Add(fhEclusSM4);
  fOutputList->Add(fhEclusSM5);
  fOutputList->Add(fhEclusSM6);
  fOutputList->Add(fhEclusSM7);
  fOutputList->Add(fhEclusSM8);
  fOutputList->Add(fhEclusSM9);
  fOutputList->Add(fhEVsTime);
  fOutputList->Add(fhEVsNcells);
  fOutputList->Add(fhZvertex);
  fOutputList->Add(fPi0);
  fOutputList->Add(fEta);
  fOutputList->Add(fhZvertexAll);
  fOutputList->Add(fhpi0Decay);
  fOutputList->Add(fhpi0Inter);
  fOutputList->Add(fhRecDecay);
  fOutputList->Add(fhRecInter);
  fOutputList->Add(fPrim);
  fOutputList->Add(fPrimNoTRD);
  fOutputList->Add(fK0);
  fOutputList->Add(fPi0mixed);
  fOutputList->Add(fEtamixed);
  fOutputList->Add(fhHitMap);
  fOutputList->Add(fNcl);
  fOutputList->Add(fK0Pt);
  fOutputList->Add(fhLowm);
  fOutputList->Add(fhHighm);
  fOutputList->Add(fPi0Conv);
  fOutputList->Add(fhLowmConv);
  fOutputList->Add(fhHighmConv);
  fOutputList->Add(fhPhotonEtaPhi);
  fOutputList->Add(fPrimOldGeo);

  fGeom =  AliEMCALGeometry::GetInstance(fGeoName);

  fRecoUtil = new AliEMCALRecoUtils();

  f1=new TF1("f1", "(8+50.*TMath::Exp(-(x-0.45)/0.95))", 0., 30.);
  f2=new TF1("f2", "(-12-TMath::Exp(-(x-3.8)/0.98))", 0., 30.);

  f3=new TF1("f3", "(1.)/([0]*(1./(1.+[1]*exp(-x/[2]))*1./(1.+[3]*exp((x-[4])/[5]))))", 0.1, 50.); // sim Non Lin
  f3->SetParameters(9.81039e-01, 1.13508e-01, 1.00173e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01);
 
  PostData(1, fOutputList);

}
//________________________________________________________________________
void AliAnalysisTaskNeutralMesons::UserExec(Option_t *) 
{
  // Main loop, Called for each event

   fESD = InputEvent();
  if (!fESD) 
    {
      AliError("No event found");
      return;
    }

  Double_t vertex[3] ; fESD->GetPrimaryVertex()->GetXYZ(vertex) ;
 
  fMC = MCEvent();
  if (!fMC) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
  
  AliStack* stack=fMC->Stack();

  if (fMCtype==2) {
    Bool_t jetPt=ComparePtHardAndJetPt();
    Bool_t clusterPt=ComparePtHardAndClusterPt();
    
    if(jetPt==kFALSE) return;
    if(clusterPt==kFALSE) return;
   
  }	

 if (!EsdVertexOk()) return;

  for (Int_t iParticle = 0; iParticle < stack->GetNprimary(); iParticle++) {
    fparticle = stack->Particle(iParticle);
    if (!fparticle) continue;

 if (TMath::Abs(fparticle->GetPdgCode())==310) {
        if (TMath::Abs(fparticle->Y())<0.9) {
   //             if(((TMath::RadToDeg()*fparticle->Phi())>80.)&&((TMath::RadToDeg()*fparticle->Phi())<180.)) {
                        fK0Pt->Fill(fparticle->Pt());
     //           }
        }
    }

   if (TMath::Abs(fparticle->GetPdgCode())==311) {
        if (TMath::Abs(fparticle->Eta())<0.7) {
		if(((TMath::RadToDeg()*fparticle->Phi())>80.)&&((TMath::RadToDeg()*fparticle->Phi())<180.)) {
			fK0->Fill(fparticle->Pt());
   		}
	}
    }

    if (TMath::Abs(fparticle->GetPdgCode())==fMCpart) {

      if (TMath::Abs(fparticle->Eta())<0.7) {
	if(((TMath::RadToDeg()*fparticle->Phi())>80.)&&((TMath::RadToDeg()*fparticle->Phi())<180.)) {
	  fPrim->Fill(fparticle->Pt());
          if (((TMath::RadToDeg()*fparticle->Phi())>80.)&&((TMath::RadToDeg()*fparticle->Phi())<140.)) fPrimNoTRD->Fill(fparticle->Pt());
	  if (((TMath::RadToDeg()*fparticle->Phi())>80.)&&((TMath::RadToDeg()*fparticle->Phi())<120.)) fPrimOldGeo->Fill(fparticle->Pt());
	}
      }
    }
  }
  
  for (Int_t iParticle = 0; iParticle < stack->GetNtrack(); iParticle++) {
    
   fparticle = stack->Particle(iParticle);
   if (!fparticle) continue;

    // if (fparticle->GetPdgCode()==311) Printf("K0");
    // if (fparticle->GetPdgCode()==310) Printf("K0s");
    // if (fparticle->GetPdgCode()==130) Printf("K0L");

    if (fparticle->GetPdgCode()==111) {
	   if (TMath::Abs(fparticle->Eta())<0.7) {
	     if(((TMath::RadToDeg()*fparticle->Phi())>80.)&&((TMath::RadToDeg()*fparticle->Phi())<180.)) {
	       
	       if ((fparticle->GetUniqueID()==4)) fhpi0Decay->Fill(fparticle->R(), fparticle->Vz());
	       if ((fparticle->GetUniqueID()==13)) fhpi0Inter->Fill(fparticle->R(), fparticle->Vz());
		 
	     }
	   }
    }
	   
  }

  fclusterList = dynamic_cast<TClonesArray*> (AODEvent()->FindListObject(fBranchName)); 
  Int_t kNumberOfEMCALClusters   = fclusterList->GetEntries() ;

  
     for(Int_t i=0; i<kNumberOfEMCALClusters; i++)
       {
	 AliAODCaloCluster * cluster = (AliAODCaloCluster *) fclusterList->At(i);
	 MakeClusterCorrections(cluster);
    }
 
     TObjArray *newPool = new TObjArray(300);
     TArrayD *newVtxPool = new TArrayD(300);
     TArrayD *newClsPool = new TArrayD(300);
 
     newPool->SetOwner();
     Int_t nGoodClusters = 0;
     Int_t nselcl = 0;
 
     for (Int_t i=0; i<kNumberOfEMCALClusters; i++) {
       AliAODCaloCluster *ci = (AliAODCaloCluster *) fclusterList->At(i);
       if (!IsGoodCluster(ci))
  	continue;
       nselcl=nselcl+1;
     }

     for (Int_t i=0; i<kNumberOfEMCALClusters; i++) {
       AliAODCaloCluster *ci = (AliAODCaloCluster *) fclusterList->At(i);
       if (!IsGoodCluster(ci))
	 continue;
       
       TLorentzVector *pi = new TLorentzVector;
       ci->GetMomentum(*pi,vertex);
       newPool->Add(pi);
       newVtxPool->SetAt(vertex[2],i);
       newClsPool->SetAt(nselcl,i);
       
       ++nGoodClusters;
     }
     
     Int_t poolFlag=0;
     if ((fPool)&&((fPool->GetEntries()+nGoodClusters)>1499)) {
       delete fPool;
       fPool = newPool;
       poolFlag=1;
       for (Int_t ka=0; ka<1500; ka++) {
	 fClusterZvtx->SetAt(0.,ka);
	 fNEvClusters->SetAt(0.,ka);
       }
       for (Int_t ki=0; ki<nGoodClusters; ki++) {
	 fClusterZvtx->SetAt(Double_t(newVtxPool->At(ki)),ki);
	 fNEvClusters->SetAt(newClsPool->At(ki),ki);
       }
     }

     if (!fPool) {
       fPool = newPool;
       for (Int_t kp=0; kp<nGoodClusters; kp++) {
	 if (kp>299) continue;
	 fClusterZvtx->SetAt(Double_t(newVtxPool->At(kp)),kp);
	 fNEvClusters->SetAt(newClsPool->At(kp),kp); 
       }
     } else {

       Int_t nGoodClusters2 = fPool->GetEntries();
       for (Int_t i=0; i<nGoodClusters; i++) {
	 TLorentzVector *pi = static_cast<TLorentzVector*>(newPool->At(i));
	 for (Int_t j=0; j<nGoodClusters2; j++) {
	   TLorentzVector *pj = static_cast<TLorentzVector*>(fPool->At(j));
	   if (((TMath::Abs(vertex[2])-fClusterZvtx->At(j))<5.)&&(TMath::Abs(nGoodClusters-fNEvClusters->At(j))<5))
	     if (poolFlag==0) MakeInvMassMixed(*pi,*pj);
	 }
       }
       for (Int_t k=0; k<nGoodClusters; k++) {
 	
	 AliAODCaloCluster *ck = (AliAODCaloCluster *) fclusterList->At(k);
	 if (!IsGoodCluster(ck))
	   continue;
	 
 	TLorentzVector *pk = new TLorentzVector;
	ck->GetMomentum(*pk,vertex);
 	fPool->Add(pk);
 	fClusterZvtx->SetAt(newVtxPool->At(k),fPool->GetEntries()-1+k);
 	fNEvClusters->SetAt(newClsPool->At(k),fPool->GetEntries()-1+k);
       }
     }    
     delete newVtxPool;
     delete newClsPool;
 
     fNeventCls=0;
     fTrpi0=0;

     fhZvertexAll->Fill(vertex[2]);

     for(Int_t icl=0; icl<kNumberOfEMCALClusters; icl++)
       {
	 
	 AliAODCaloCluster * cluster = (AliAODCaloCluster *) fclusterList->At(icl);
	 if(!IsGoodCluster(cluster)) continue;

	 fNeventCls=fNeventCls+1;
	 fhEVsTime->Fill(cluster->E(), cluster->GetTOF()*TMath::Power(10,9));
	 fhEVsNcells->Fill(cluster->E(), cluster->GetNCells());
	 TLorentzVector pii;
	 cluster->GetMomentum(pii,vertex);
	 fhHitMap->Fill(pii.Eta(), TMath::RadToDeg()*pii.Phi());
	 
	 for(Int_t jcl=icl+1; jcl<kNumberOfEMCALClusters; jcl++)
	   {
	     AliAODCaloCluster * cluster2 = (AliAODCaloCluster *) fclusterList->At(jcl);
	     if(!IsGoodCluster(cluster2)) continue;
	  
	     TLorentzVector pjj;
	     cluster2->GetMomentum(pjj,vertex);
	     MakeInvMass(pii,pjj);
	     
	   } // Inner cluster loop
	 
       } // outer cluster loop
     
     PostData(1, fOutputList);
}

//________________________________________________________________________
//
void AliAnalysisTaskNeutralMesons::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));

  if (!fOutputList) {
    Printf("ERROR: fList not available");
    return;
  }  
}

//________________________________________________________________________
//
Bool_t AliAnalysisTaskNeutralMesons::EsdVertexOk() const
{
  const AliVVertex* vtx = fESD->GetPrimaryVertex();
  Int_t nContributors = vtx->GetNContributors();
   
  Double_t vertexA[3] ; fESD->GetPrimaryVertex()->GetXYZ(vertexA) ;
  if ( nContributors < 1) return kFALSE;
  if (TMath::Abs(vertexA[2]) > 10.) return kFALSE;
  fhZvertex->Fill(vertexA[2]);
  return kTRUE;
}

//________________________________________________________________________
//
Bool_t AliAnalysisTaskNeutralMesons::IsGoodCluster(AliVCluster *cluster)
{
  if(!cluster) return kFALSE;
  if (!cluster->IsEMCAL()) return kFALSE;

  Int_t absID = -1, nSupMod=-1, iphi=-1, ieta=-1;
  Bool_t lshare = kFALSE;
  AliVCaloCells *lcells = (AliVCaloCells*)fESD->GetEMCALCells();
  fRecoUtil->GetMaxEnergyCell(fGeom, lcells, cluster, absID, nSupMod, ieta, iphi, lshare);

  if (fNcellsInClusterFlag==kTRUE) if (cluster->GetNCells()<2) return kFALSE;

  if (cluster->E()<fLowEnergCut) return kFALSE;

  //  if((cluster->E()>3.)&&fRejectExoticCluster && fRecoUtil->IsExoticCluster(cluster, (AliVCaloCells*)fESD->GetEMCALCells())) return kFALSE;

 if ((iphi==0)||(iphi==23)||(ieta==0)||(ieta==47)) return kFALSE;

 //const Int_t nbadChannels=237;
 //const Int_t nbadChannels=98;  //11a
 const Int_t nbadChannels=457;  //11c
 //const Int_t nbadChannels=69;  //11b

    // 11d
 //Int_t badChannelArray[nbadChannels] = {103, 1263, 1275, 1860, 2117, 2298, 2776, 3664, 3665, 3666, 3667, 3668, 3669, 3670, 3671, 3672, 3673, 3674, 3675, 3676, 3677, 3678, 3679, 3712, 3713, 3714, 3715, 3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724, 3725, 3726, 3727, 3764, 6111, 6800, 6801, 6802, 6803, 6804, 6805, 6806, 6807, 6808, 6809, 6810, 6811, 6812, 6813, 6814, 6815, 7430, 7491, 8352, 8353, 8354, 8356, 8357, 8362, 8808, 8810, 8812, 8814, 9217, 9361, 9704, 9769, 9802, 9837, 9839, 9888, 9946, 10117, 11462, 74, 152, 759, 1059, 1175, 1204, 1288, 1376, 1382, 1386, 1519, 1967, 2026, 2047, 2112, 2114, 2115, 2116, 2118, 2119, 2120, 2123, 2124, 2125, 2350, 2506, 2540, 2793, 2891, 2985, 3135, 3503, 4377, 4817, 5600, 5601, 5602, 5603, 5612, 5613, 5614, 5615, 5648, 5649, 5650, 5651, 5660, 5661, 5662, 5663, 5836, 6104, 6481, 7371, 7375, 7425, 7572, 7874, 9269, 9302, 9389, 9696, 9697, 9698, 9699, 9700, 9701, 9702, 9703, 9705, 9706, 9707, 9708, 9709, 9710, 9711, 9750, 9758, 9792, 9793, 9794, 9795, 9798, 9800, 9801, 9803, 9804, 9815, 9824, 9825, 9828, 9829, 9830, 9831, 9832, 9833, 9834, 9835, 9836, 9838, 9872, 9874, 9875, 9878, 9882, 9883, 9889, 9890, 9891, 9892, 9893, 9894, 9896, 9897, 9898, 9899, 9900, 9901, 9902, 9903, 9927, 9936, 9937, 9938, 9939, 9940, 9941, 9942, 9943, 9945, 9947, 9948, 9949, 9950, 9951, 10086, 10112, 10113, 10114, 10115, 10116, 10118, 10119, 10120, 10121, 10122, 10123, 10124, 10125, 10138, 10718, 10723, 10771, 11042, 11091, 11363};

  //Int_t badChannelArray[nbadChannels] = {2776, 8353, 8357, 9361, 11462, 74, 103, 152, 759, 917, 1059, 1175, 1204, 1212, 1275, 1288, 1366, 1376, 1382, 1386, 1414, 1519, 1704, 1738, 1836, 1837, 1844, 1860, 1967, 2022, 2026, 2047, 2112, 2114, 2115, 2116, 2120, 2123, 2210, 2298, 2424, 2487, 2506, 2534, 2540, 2586, 2793, 2888, 2891, 2915, 3135, 3236, 3732, 3748, 3914, 3974, 4129, 4530, 4543, 4817, 5648, 5649, 5650, 5651, 5660, 6104, 6111, 6331, 6481, 6811, 7089, 7371, 7417, 7425, 7430, 7457, 7491, 7874, 8352, 8353, 8354, 8355, 8356, 8358, 8360, 8361, 8362, 8365, 9217, 9269, 9769, 9891, 9895, 9897, 9927, 10203, 10363, 10771};

  //11c
 Int_t badChannelArray[nbadChannels] = {9849,9835,9831,9823,9786,9784,9782,9774,9764,9718,9381,9338,9314,9301,9218,9201,8813,8811,8809,8807,8760,8718,866,8361,8355,8353,8306,8222,7873,7571,7440,7429,7374,7324,7320,7038,6893,6891,6889,6887,6885,6883,6881,6879,6845,6844,6843,6842,6841,6840,6839,6838,6837,6836,6835,6834,6833,6832,6831,6830,6813,6811,6809,6807,6805,6803,6801,6799,6797,6796,6795,6794,6793,6792,6791,6790,6789,6788,6787,6786,6785,6784,6783,6782,6764,6762,6760,6758,6756,6754,6752,6750,6749,6748,6747,6746,6745,6744,6743,6742,6741,6740,6739,6738,6737,6736,6735,6734,6701,6700,6699,6698,6697,6696,6695,6694,6693,6692,6691,6690,6689,6688,6687,6686,6677,6653,6652,6651,6650,6649,6648,6647,6646,6645,6644,6643,6642,6641,6640,6639,6638,6605,6604,6603,6602,6601,6600,6599,6598,6597,6596,6595,6594,6593,6592,6591,6590,6557,6556,6555,6554,6553,6552,6551,6550,6549,6548,6547,6546,6545,6544,6543,6542,6510,6509,6508,6507,6506,6505,6504,6503,6502,6501,6500,6499,6498,6497,6496,6495,6494,6461,6460,6459,6458,6457,6456,6455,6454,6453,6452,6451,6450,6449,6448,6447,6446,6413,6412,6411,6410,6409,6408,6407,6406,6405,6404,6403,6402,6401,6400,6399,6398,6374,6365,6364,6363,6362,6361,6360,6359,6358,6357,6356,6355,6354,6353,6352,6351,6350,6339,6339,6318,6316,6314,6312,6310,6308,6306,6304,6302,6137,6129,6086,6075,6060,6045,6043,6041,6039,6037,6035,6033,6031,6014,5997,5996,5995,5994,5993,5992,5991,5990,5989,5988,5987,5986,5985,5984,5983,5982,5949,5948,5947,5946,5945,5944,5943,5942,5941,5940,5939,5938,5937,5936,5935,5934,5918,5901,5900,5899,5898,5897,5896,5895,5894,5893,5892,5891,5890,5889,5888,5887,5886,5853,5852,5851,5850,5849,5848,5847,5846,5845,5844,5843,5842,5841,5840,5839,5838,5805,5804,5803,5802,5801,5800,5799,5798,5797,5796,5795,5794,5793,5792,5791,5790,5559,52,4106,4070,4025,3763,3725,3723,3721,3719,3717,3715,3713,3711,3677,3676,3675,3674,3673,3672,3671,3670,3669,3668,3667,3666,3665,3664,3663,3628,3626,3624,3622,3620,3618,3616,3614,3543,3516,323,319,3084,2934,2934,2777,2775,2505,2297,2209,2119,2115,2113,2094,2073,2072,2066,2064,2025,2020,2018,1996,1972,1968,1916,1911,1910,1859,1719,1544,151,1468,1413,1383,1381,1375,1365,1287,1275,1224,1212,11461,11267,11040,10727,10717,10672,10575,10121,10119,10117,10115,10113,10111,10085,1008,10077,10072,10070,10068,10066,10064,10062};
 
  //11b  
 //Int_t badChannelArray[nbadChannels] = {1275,1860,2117,2298, 2776, 6800, 6801, 6802, 6803, 6804, 6805, 6806, 6807,6808, 6809, 6810, 6811, 6812, 6813, 6814, 6815, 7491, 8352, 8353,8356, 8357, 8808, 8810, 8812, 8814, 103, 152, 1263, 1276, 1288, 1384, 1519, 1595,1967,2071,2112, 2114,2115,2116, 2120,2123,2506,2540,2778,2793,3135,3764, 4157,4174,5085,5767,6111,6425,6481,7089, 7425,7430,7874, 7878, 8811, 9217,9769,9815, 10723};

  const UShort_t* cellList=cluster->GetCellsAbsId();
  const Int_t nCells=cluster->GetNCells();
  
  for(Int_t iCell = 0; iCell<nCells; iCell++)
    {
    
    UShort_t cellIn=cellList[iCell];
    
    for(Int_t i=0; i<nbadChannels; i++) {
      //   if (cellIn==badChannelArray[i]) return kFALSE;
    }

  }// cell cluster loop

  return kTRUE;
}

////________________________________________________________________________
//
void AliAnalysisTaskNeutralMesons::MakeClusterCorrections(AliVCluster *cluster)
{
  Int_t absID = -1, nSupMod=-1, iphi=-1, ieta=-1;
  Bool_t lshare = kFALSE;
  AliVCaloCells *lcells = (AliVCaloCells*)fESD->GetEMCALCells();
  fRecoUtil->GetMaxEnergyCell(fGeom, lcells, cluster, absID, nSupMod, ieta, iphi, lshare);
  
    // Non linearity
  Double_t clE=cluster->E();
  cluster->SetE(clE*f3->Eval(clE));

}

////________________________________________________________________________
//
void AliAnalysisTaskNeutralMesons::MakeInvMass(const TLorentzVector& p1, const TLorentzVector& p2)
{ 
  // Fill histogram.
    
  TLorentzVector pion;
  pion = p1 + p2;

  Double_t mass = pion.M();
  Double_t pt   = pion.Pt();

  if ((mass<0.199)&&(mass>0.085)&&(pt>14.)) {
    if (fTrpi0==0) fNcl->Fill(fNeventCls);
    fTrpi0=1;
  }

  fPi0->Fill(mass,pt);

  Double_t eta=pion.Eta();
  Double_t phi=pion.Phi();
  if (mass<0.05) fhPhotonEtaPhi->Fill(eta, TMath::RadToDeg()*phi);
   
  if (pt<2.) fhLowm->Fill(mass,pt);
  if (pt>4.) fhHighm->Fill(mass,pt);

  Double_t asym = (p1.E()-p2.E())/(p1.E()+p2.E());
    
  if (asym<0.7) {
    fEta->Fill(mass,pt);
 
  }

}

//_____________________________________________________

void AliAnalysisTaskNeutralMesons::MakeInvMassMixed(const TLorentzVector& p1, const TLorentzVector& p2)
{ 
  // Fill histogram.
    
  TLorentzVector pion;
  pion = p1 + p2;

  Double_t mass = pion.M();
  Double_t pt   = pion.Pt();
  
  // investigate time

  fPi0mixed->Fill(mass,pt);

  Double_t asym = (p1.E()-p2.E())/(p1.E()+p2.E());
    
  if (asym<0.7) {
   // fEta->Fill(mass,pt);
    fEtamixed->Fill(mass,pt);
 
  }
}
//________________________________________________
Bool_t AliAnalysisTaskNeutralMesons::ComparePtHardAndJetPt()
{
  // Check the event, if the requested ptHard is much smaller than the jet pT, then there is a problem.
  // Only for PYTHIA.
  
  if(!strcmp(GetGenEventHeader()->ClassName(), "AliGenPythiaEventHeader"))
  {
    TParticle * jet =  0;
    AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) GetGenEventHeader();
    Int_t nTriggerJets =  pygeh->NTriggerJets();
    Float_t ptHard = pygeh->GetPtHard();
    
    Float_t tmpjet[]={0,0,0,0};
    for(Int_t ijet = 0; ijet< nTriggerJets; ijet++)
      {
	pygeh->TriggerJet(ijet, tmpjet);
	jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
	
	//Compare jet pT and pt Hard
	if(jet->Pt() > 4.* ptHard)
	  {
	    return kFALSE;
	  }
      }
    
    if(jet) delete jet;
  }
  
  return kTRUE ;
  
}

//____________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesons::ComparePtHardAndClusterPt()
{
  // Check the event, if the requested ptHard is smaller than the calorimeter cluster E, then there is a problem.
  // Only for PYTHIA.
  
  if(!strcmp(GetGenEventHeader()->ClassName(), "AliGenPythiaEventHeader"))
  {
    AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) GetGenEventHeader();
    Float_t ptHard = pygeh->GetPtHard();
    
    Int_t nclusters = fESD->GetNumberOfCaloClusters();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++)
      {
	AliVCluster * clus = fESD->GetCaloCluster(iclus) ;
	Float_t ecluster = clus->E();
	
	if(ecluster >1.5*ptHard)
	  {
	    return kFALSE;
	  }
      }
    
  }
  
  return kTRUE ;
  
}

//______________________________________________________________

AliGenEventHeader* AliAnalysisTaskNeutralMesons::GetGenEventHeader() const
{

  AliGenEventHeader * eventHeader = fMC->GenEventHeader();

  return eventHeader ;

}



