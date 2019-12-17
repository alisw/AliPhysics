//isputows
#include <numeric>
#include <functional>
#include <set>

#include <TChain.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultSelection.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAnalysisTaskStronglyIntensiveCorrTree.h"


ClassImp(AliAnalysisTaskStronglyIntensiveCorrTree);
ClassImp(TrackInfoCorr);

AliAnalysisTaskStronglyIntensiveCorrTree::AliAnalysisTaskStronglyIntensiveCorrTree() : AliAnalysisTaskSE()
  , fTree(NULL)
  , fTreeMC(NULL)
  , fIsMC(kFALSE)
  , fTrackFilter(0)
  , fTrigger(0)
  , fPtMin(0.2), fPtMax(5)
  , fPhiMin(0.), fPhiMax(TMath::TwoPi())
  , fMagneticField(0)
  , fRunNumber(0)
  , fEventCuts(0)
  , fHistEta(0)
  , fHistPt(0)
  , fHistPhi(0)
  , fHistDcaX(0)
  , fHistDcaY(0)
  , fHistDcaZ(0)
  , fHist2D_EtaDcaX(0)
  , fHist2D_EtaDcaY(0)
  , fHist2D_EtaDcaZ(0)
  , fHist2D_EtaPhi(0)
  , fHist2D_EtaPt(0)
  , fEventStatistics(0)
  , fHistEtaMC(0)
  , fHistPtMC(0)
  , fHistPhiMC(0)
  , fHistDcaXMC(0)
  , fHistDcaYMC(0)
  , fHistDcaZMC(0)
  , fHist2D_EtaDcaXMC(0)
  , fHist2D_EtaDcaYMC(0)
  , fHist2D_EtaDcaZMC(0)
  , fHist2D_EtaPhiMC(0)
  , fHist2D_EtaPtMC(0)
  , fHistEtaPrim(0)
  , fHistPtPrim(0)
  , fHistPhiPrim(0)
  , fHistDcaXPrim(0)
  , fHistDcaYPrim(0)
  , fHistDcaZPrim(0)
  , fHist2D_EtaDcaXPrim(0)
  , fHist2D_EtaDcaYPrim(0)
  , fHist2D_EtaDcaZPrim(0)
  , fHist2D_EtaPhiPrim(0)
  , fHist2D_EtaPtPrim(0)
  , fVertexZ(0){
  for (Int_t i=0; i<15; ++i) {
    if (i<4) fCent[i] = 0;
     fNf[i] = fNb[i] = fNf_MC[i]=fNb_MC[i] = fNf_MCPrim[i]=fNb_MCPrim[i]= 0;

  }
}
////////-----
AliAnalysisTaskStronglyIntensiveCorrTree::AliAnalysisTaskStronglyIntensiveCorrTree(const char *name) : AliAnalysisTaskSE(name)
  , fTree(NULL)
  , fTreeMC(NULL)
  , fIsMC(kFALSE)
  , fTrackFilter(0)
  , fTrigger(0)
  , fPtMin(0.2), fPtMax(5)
  , fPhiMin(0.), fPhiMax(TMath::TwoPi())
  , fMagneticField(0)
  , fRunNumber(0)
  , fEventCuts(0)
  , fHistEta(0)
  , fHistPt(0)
  , fHistPhi(0)
  , fHistDcaX(0)
  , fHistDcaY(0)
  , fHistDcaZ(0)
  , fHist2D_EtaDcaX(0)
  , fHist2D_EtaDcaY(0)
  , fHist2D_EtaDcaZ(0)
  , fHist2D_EtaPhi(0)
  , fHist2D_EtaPt(0)
  , fEventStatistics(0)
  , fHistEtaMC(0)
  , fHistPtMC(0)
  , fHistPhiMC(0)
  , fHistDcaXMC(0)
  , fHistDcaYMC(0)
  , fHistDcaZMC(0)
  , fHist2D_EtaDcaXMC(0)
  , fHist2D_EtaDcaYMC(0)
  , fHist2D_EtaDcaZMC(0)
  , fHist2D_EtaPhiMC(0)
  , fHist2D_EtaPtMC(0)
  , fHistEtaPrim(0)
  , fHistPtPrim(0)
  , fHistPhiPrim(0)
  , fHistDcaXPrim(0)
  , fHistDcaYPrim(0)
  , fHistDcaZPrim(0)
  , fHist2D_EtaDcaXPrim(0)
  , fHist2D_EtaDcaYPrim(0)
  , fHist2D_EtaDcaZPrim(0)
  , fHist2D_EtaPhiPrim(0)
  , fHist2D_EtaPtPrim(0)
  , fVertexZ(0){
  for (Int_t i=0; i<15; ++i) {
          if (i<4) fCent[i] = 0;
                    fNf[i] = fNb[i] = fNf_MC[i]=fNb_MC[i] = fNf_MCPrim[i]=fNb_MCPrim[i]= 0;

  }
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
  if(fIsMC){
          DefineOutput(3, TTree::Class());
          DefineOutput(4, TTree::Class());
  }
}

AliAnalysisTaskStronglyIntensiveCorrTree::~AliAnalysisTaskStronglyIntensiveCorrTree() {
  if (NULL != fTree) {
    delete fTree;
    fTree = NULL;
  }
  if (NULL != fTreePrim) {
    delete fTreePrim;
    fTreePrim = NULL;
  }
  if (NULL != fTreeMC) {
    delete fTreeMC;
    fTreeMC = NULL;
  }
   if(fOutputList) 	  delete fOutputList;
}

void AliAnalysisTaskStronglyIntensiveCorrTree::UserCreateOutputObjects() {
 

  TDirectory *owd = NULL;
  owd = gDirectory;
  
  OpenFile(1);
  
  
  fTree = new TTree;
  fTree->SetName("TreeLRC");

  fTree->Branch("RunNumber",  &fRunNumber);
  fTree->Branch("VertexZ",    &fVertexZ);
  fTree->Branch("Centrality", &fCent, "V0M/F:TRK:ZEMvsZDC:b");
  fTree->Branch("MagneticFiled", &fMagneticField);
  fTree->Branch("Nf", &fNf, "n[15]/S");
  fTree->Branch("Nb", &fNb, "n[15]/S");
  if(fIsMC){
          OpenFile(2);
          fTreePrim = new TTree;
          fTreePrim->SetName("TreeLRC_MCPrim");

          fTreePrim->Branch("RunNumber",  &fRunNumber);
          fTreePrim->Branch("VertexZ",    &fVertexZ);
          fTreePrim->Branch("Centrality", &fCent, "V0M/F:TRK:ZEMvsZDC:b");
          fTreePrim->Branch("MagneticFiled", &fMagneticField);
          fTreePrim->Branch("Nf", &fNf_MCPrim, "n[15]/S");
          fTreePrim->Branch("Nb", &fNb_MCPrim, "n[15]/S");
  }
  // Event statistics
  fEventStatistics = new TH1I("fEventStatistics","",10,0,10);
  
  fHistEta  = new TH1F ( "fHistEta" ,  "fHistEta " ,  40 , -1.0 ,1.0) ;
  fHistPt   = new TH1F ( "fHistPt" ,   "fHistPt " ,   220 , -0.5 ,10.5) ;
  fHistPhi  = new TH1F ( "fHistPhi" ,  "fHistPhi " ,  160 , 0 ,TMath::TwoPi()) ;
  fHistDcaX = new TH1F ( "fHistDcaX" , "fHistDcaX " , 400 , -10.,10.) ;
  fHistDcaY = new TH1F ( "fHistDcaY" , "fHistDcaY " , 400 , -10.,10.) ;
  fHistDcaZ = new TH1F ( "fHistDcaZ" , "fHistDcaZ " , 400 , -10.,10.) ;
  
  
  fHist2D_EtaDcaX   = new TH2F ( "fHist2D_EtaDcaX" , "fHist2D_EtaDcaX"     , 40 , -1.0 ,1.0, 400 , -10.,10.);
  fHist2D_EtaDcaY   = new TH2F ( "fHist2D_EtaDcaY" , "fHist2D_EtaDcaY"     , 40 , -1.0 ,1.0, 400 , -10.,10.);
  fHist2D_EtaDcaZ   = new TH2F ( "fHist2D_EtaDcaZ" , "fHist2D_EtaDcaZ"     , 40 , -1.0 ,1.0, 400 , -10.,10.) ;
  fHist2D_EtaPt     = new TH2F ( "fHist2D_EtaPt" ,   "fHist2D_EtaPt"       , 40 , -1.0 ,1.0, 220 , -.5,10.5) ;
  fHist2D_EtaPhi    = new TH2F ( "fHist2D_EtaPhi" ,  "fHist2D_EtaEtaPhi"   , 40 , -1.0 ,1.0, 160 ,  0 ,TMath::TwoPi()) ;
  
  if(fIsMC){
          fHistEtaMC  = new TH1F ( "fHistEtaMC" ,  "fHistEtaMC" ,  40 , -1.0 ,1.0) ;
          fHistPtMC   = new TH1F ( "fHistPtMC" ,   "fHistPtMC" ,   220 , -0.5 ,10.5) ;
          fHistPhiMC  = new TH1F ( "fHistPhiMC" ,  "fHistPhiMC" ,  160 , 0 ,TMath::TwoPi()) ;
 
          fHist2D_EtaPtMC     = new TH2F ( "fHist2D_EtaPtMC" ,   "fHist2D_EtaPtMC"       , 40 , -1.0 ,1.0, 220 , -.5,10.5) ;
          fHist2D_EtaPhiMC    = new TH2F ( "fHist2D_EtaPhiMC" ,  "fHist2D_EtaEtaPhiMC"   , 40 , -1.0 ,1.0, 160 ,  0 ,TMath::TwoPi()) ;
          
          fHistEtaPrim  = new TH1F ( "fHistEtaPrim" ,  "fHistEtaPrim" ,  40 , -1.0 ,1.0) ;
          fHistPtPrim   = new TH1F ( "fHistPtPrim" ,   "fHistPtPrim" ,   220 , -0.5 ,10.5) ;
          fHistPhiPrim  = new TH1F ( "fHistPhiPrim" ,  "fHistPhiPrim" ,  160 , 0 ,TMath::TwoPi()) ;
          fHistDcaXPrim = new TH1F ( "fHistDcaXPrim" , "fHistDcaXPrim" , 400 , -10.,10.) ;
          fHistDcaYPrim = new TH1F ( "fHistDcaYPrim" , "fHistDcaYPrim" , 400 , -10.,10.) ;
          fHistDcaZPrim = new TH1F ( "fHistDcaZPrim" , "fHistDcaZPrim" , 400 , -10.,10.) ;
  
  
          fHist2D_EtaDcaXPrim   = new TH2F ( "fHist2D_EtaDcaXPrim" , "fHist2D_EtaDcaXPrim" , 40 , -1.0 ,1.0, 400 , -10.,10.) ;
          fHist2D_EtaDcaYPrim   = new TH2F ( "fHist2D_EtaDcaYPrim" , "fHist2D_EtaDcaYPrim" , 40 , -1.0 ,1.0, 400 , -10.,10.) ;
          fHist2D_EtaDcaZPrim   = new TH2F ( "fHist2D_EtaDcaZPrim" , "fHist2D_EtaDcaZPrim" , 40 , -1.0 ,1.0, 400, -10.,10.) ;
          fHist2D_EtaPtPrim     = new TH2F ( "fHist2D_EtaPtPrim",   "fHist2D_EtaPtPrim"    , 40 , -1.0 ,1.0, 220 , -.5,10.5) ;
          fHist2D_EtaPhiPrim    = new TH2F ( "fHist2D_EtaPhiPrim",  "fHist2D_EtaEtaPhiPrim", 40 , -1.0 ,1.0, 160 ,  0 ,TMath::TwoPi()) ;
  }
  fOutputList=new TList();
  fOutputList->SetOwner();
  fOutputList->Add(fHistEta);
  fOutputList->Add(fHistPt);
  fOutputList->Add(fHistPhi);
  fOutputList->Add(fHistDcaX);
  fOutputList->Add(fHistDcaY);
  fOutputList->Add(fHistDcaZ);
  fOutputList->Add(fHist2D_EtaDcaX);
  fOutputList->Add(fHist2D_EtaDcaY);
  fOutputList->Add(fHist2D_EtaDcaZ);
  fOutputList->Add(fHist2D_EtaPt);
  fOutputList->Add(fHist2D_EtaPhi);
  fOutputList->Add(fEventStatistics);
  
  if(fIsMC){
          fOutputList->Add(fHistEtaMC);
          fOutputList->Add(fHistPtMC);
          fOutputList->Add(fHistPhiMC);
          fOutputList->Add(fHistDcaXMC);
          fOutputList->Add(fHistDcaYMC);
          fOutputList->Add(fHistDcaZMC);
          fOutputList->Add(fHist2D_EtaDcaXMC);
          fOutputList->Add(fHist2D_EtaDcaYMC);
          fOutputList->Add(fHist2D_EtaDcaZMC);
          fOutputList->Add(fHist2D_EtaPtMC);
          fOutputList->Add(fHist2D_EtaPhiMC);
  
          fOutputList->Add(fHistEtaPrim);
          fOutputList->Add(fHistPtPrim);
          fOutputList->Add(fHistPhiPrim);
          fOutputList->Add(fHistDcaXPrim);
          fOutputList->Add(fHistDcaYPrim);
          fOutputList->Add(fHistDcaZPrim);
          fOutputList->Add(fHist2D_EtaDcaXPrim);
          fOutputList->Add(fHist2D_EtaDcaYPrim);
          fOutputList->Add(fHist2D_EtaDcaZPrim);
          fOutputList->Add(fHist2D_EtaPtPrim);
          fOutputList->Add(fHist2D_EtaPhiPrim);
          
          OpenFile(3);
          fTreeMC = new TTree;
          fTreeMC->SetName("TreeLRC_MC");
          fTreeMC->Branch("RunNumber",  &fRunNumber);
          fTreeMC->Branch("VertexZ",    &fVertexZ);
          fTreeMC->Branch("Centrality", &fCent, "V0M/F:TRK:ZEMvsZDC:b");
          fTreeMC->Branch("MagneticFiled", &fMagneticField);
          fTreeMC->Branch("Nf", &fNf_MC, "n[15]/S");
          fTreeMC->Branch("Nb", &fNb_MC, "n[15]/S");
       
  }

  if (NULL != owd) owd->cd();

  PostData(1, fTree);
  PostData (2,  fOutputList);
  
 if (fIsMC){
          PostData(3,fTreePrim);
          PostData(4, fTreeMC);
 }
}

void AliAnalysisTaskStronglyIntensiveCorrTree::UserExec(Option_t* ) {
  AliAnalysisManager* pManager(AliAnalysisManager::GetAnalysisManager());
  if (NULL == pManager) return;

  AliInputEventHandler* pInputHandler(dynamic_cast<AliInputEventHandler*>(pManager->GetInputEventHandler()));
  if (NULL == pInputHandler) return;

  AliAODEvent* pAOD(dynamic_cast<AliAODEvent*>(InputEvent()));
  if (NULL == pAOD) return;
  fEventStatistics->Fill("before cuts",1);
  
  AliAODHeader *pAODHeader = dynamic_cast<AliAODHeader*>(pAOD->GetHeader());
  if (NULL == pAODHeader) return;
  fEventStatistics->Fill("after aod check",1);
  
  AliAODMCHeader *pAODMCHeader(dynamic_cast<AliAODMCHeader*>(pAOD->FindListObject(AliAODMCHeader::StdBranchName())));
  const Bool_t fIsMC(NULL != pAODMCHeader);
  

  TClonesArray *arrayMC(NULL);
  if (fIsMC) {
          arrayMC = dynamic_cast<TClonesArray*>(pAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
          if (NULL == arrayMC) return;
  }

  // physics selection
  //if (!pInputHandler->IsEventSelected()) return;
  
  UInt_t fSelectMask= fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask & fTrigger;
  
  //cout<<"Trigger "<<fTrigger<<endl;
  //cout<<"FilterBit "<<fTrackFilter<<endl;
  
  if(!isINT7selected) { return;}
  fEventStatistics->Fill("physics selection",1);
   
  fRunNumber = pAOD->GetRunNumber();
  AliMultSelection *multSelection =static_cast<AliMultSelection*>(pAOD->FindListObject("MultSelection"));
  if(multSelection) fCent[0] = multSelection->GetMultiplicityPercentile("V0M");
  if(multSelection) fCent[1] = multSelection->GetMultiplicityPercentile("TRK");
  if(multSelection) fCent[2] = multSelection->GetMultiplicityPercentile("ZEMvsZDC");
  
 
  fMagneticField=pAOD->GetMagneticField();
  //pAOD->GetPrimaryVertex()->Print();

  // vertex selection -- data
  const Int_t nVertex(pAOD->GetNumberOfVertices());
  if (0 == nVertex) return;
  
  const AliAODVertex* pVertex(pAOD->GetPrimaryVertex());
  if (NULL == pVertex) return;
  fEventStatistics->Fill("found primary vertex",1);
  
  if (!fEventCuts.AcceptEvent(fInputEvent)) return;
    fEventStatistics->Fill("AliEventCuts",1);
    
    
  const Int_t nTracksPrimary(pVertex->GetNContributors());
  const Int_t nDaughters(pVertex->GetNDaughters());
  if (nTracksPrimary < 1 && nDaughters < 1) return;

  fVertexZ = pVertex->GetZ();
   Double_t v[3] = { 0. };
  pVertex->GetXYZ(v);
  // data
  for (Int_t j(0); j<15; ++j)
    fNf[j] = fNb[j] = 0;
  if(abs(fVertexZ)<10.){fEventStatistics->Fill("Z vertex cut",1);}
 
  if(fVertexZ<10. && fVertexZ>-10 && fCent[0]>0. && fCent[0] <80.0 ){fEventStatistics->Fill("Centrality 0-80",1);}
 
  TObjArray* tracksData(GetAcceptedTracks(pAOD, arrayMC,0));
  if (NULL != tracksData) {
          for (Long64_t i(0), n(tracksData->GetEntries()); i<n; ++i) {
                    TrackInfoCorr *lp = dynamic_cast<TrackInfoCorr*>(tracksData->At(i));
                    if (NULL == lp) continue;
                    const Double_t eta(lp->Eta());
                    const Double_t pt(lp->Pt());
                    const Double_t dcaX(lp->XAtDCA()-v[0]);
                    const Double_t dcaY(lp->YAtDCA()-v[1]);
                    const Double_t dcaZ(lp->ZAtDCA());
                    const Double_t phi(lp->Phi());
                    if(fVertexZ<10. && fVertexZ>-10 && eta >-0.8 && eta <0.8 ){
                              fHistEta->Fill(eta);
        	                    fHistPt->Fill(pt);
        	                    fHistPhi->Fill(phi);
        	                    fHistDcaX->Fill(dcaX);
        	                    fHistDcaY->Fill(dcaY);
        	                    fHistDcaZ->Fill(dcaZ);
        
        	                    fHist2D_EtaDcaX->Fill(eta,dcaX);
        	                    fHist2D_EtaDcaY->Fill(eta,dcaY);
        	                    fHist2D_EtaDcaZ->Fill(eta,dcaZ);
        	                    fHist2D_EtaPt->Fill(eta,pt);
        	                    fHist2D_EtaPhi->Fill(eta,phi);
                    }
      
                    for (Int_t j=0; j<15; ++j) {
	          if (eta >=  0.1*j     && eta <  0.1*j+0.2)
	                    ++fNf[j];
	          if (eta >=  -0.1*j-0.2 && eta <  -0.1*j)
	                    ++fNb[j];
                    }
          }
   
          fTree->Fill();
          delete tracksData;
  }
  if(fIsMC) {
          for (Int_t j(0); j<15; ++j)
          fNf_MCPrim[j] = fNb_MCPrim[j] = 0;

          TObjArray* tracksMCPrim(GetAcceptedTracks(pAOD, arrayMC,1));
          if (NULL != tracksMCPrim) {
                    for (Long64_t i(0), n(tracksMCPrim->GetEntries()); i<n; ++i) {
                              TrackInfoCorr *lp = dynamic_cast<TrackInfoCorr*>(tracksMCPrim->At(i));
                              if (NULL == lp) continue;
                              const Double_t eta(lp->Eta());
                              const Double_t pt(lp->Pt());
                              const Double_t dcaX(lp->XAtDCA());
                              const Double_t dcaY(lp->YAtDCA());
                              const Double_t dcaZ(lp->ZAtDCA());
                              const Double_t phi(lp->Phi());
      
                              if(fVertexZ<10. && fVertexZ>-10 && eta >-0.8 && eta <0.8 && fIsMC){
        	                              fHistEtaPrim->Fill(eta);
        	                              fHistPtPrim->Fill(pt);
        	                              fHistPhiPrim->Fill(phi);
        	                              fHistDcaXPrim->Fill(dcaX);
        	                              fHistDcaYPrim->Fill(dcaY);
        	                              fHistDcaZPrim->Fill(dcaZ);
        
        	                              fHist2D_EtaDcaXPrim->Fill(eta,dcaX);
        	                              fHist2D_EtaDcaYPrim->Fill(eta,dcaY);
        	                              fHist2D_EtaDcaZPrim->Fill(eta,dcaZ);
        	                              fHist2D_EtaPtPrim->Fill(eta,pt);
        	                              fHist2D_EtaPhiPrim->Fill(eta,phi);
                              }
                              for (Int_t j=0; j<15; ++j) {
	                              if (eta >=  0.1*j     && eta <  0.1*j+0.2)
	                                        ++fNf_MCPrim[j];
	                              if (eta >=  -0.1*j-0.2 && eta <  -0.1*j)
	                                        ++fNb_MCPrim[j];
                              }
          }
    //printf("tracksData: (%5d) ", tracksMCPrim->GetEntries());
    
          fTreePrim->Fill();
          delete tracksMCPrim;
          }
  }
  // MC
  if (fIsMC && NULL != arrayMC) {
          for (Int_t j(0); j<15; ++j)
          fNf_MC[j] = fNb_MC[j] = 0;
    
          TObjArray* tracksMC(GetAcceptedTracks(arrayMC,1));
          if (NULL != tracksMC) {      
                    for (Long64_t i(0), n(tracksMC->GetEntries()); i<n; ++i) {
	                    TrackInfoCorr *lp = dynamic_cast<TrackInfoCorr*>(tracksMC->At(i));
	                    if (NULL == lp) continue;
	                    const Double_t eta(lp->Eta());
      	                    const Double_t pt(lp->Pt());
      	                    const Double_t phi(lp->Phi());
      	                    if(fVertexZ<10. && fVertexZ>-10 && eta >-0.8 && eta <0.8 ){
        	                              fHistEtaMC->Fill(eta);
        	                              fHistPtMC->Fill(pt);
        	                              fHistPhiMC->Fill(phi);
        	                              fHist2D_EtaPtMC->Fill(eta,pt);
        	                              fHist2D_EtaPhiMC->Fill(eta,phi);
      	                    }
	                    for (Int_t j=0; j<15; ++j) {
	                              if (eta >=  0.1*j     && eta <  0.1*j+0.2)
	                                        ++fNf_MC[j];
	                              if (eta >=  -0.1*j-0.2 && eta <  -0.1*j)
	                                         ++fNb_MC[j];
	                    }
                    }
     // printf("tracksMC  : (%5d) ", tracksMC->GetEntries());
                    //for (Int_t j=0; j<15; ++j) printf("%3d,%3d|", fNf_MC[j],fNb_MC[j]);
                   // printf("\n");
                    fTreeMC->Fill();
                    delete tracksMC;
          }
  }
  PostData(1, fTree);
 
  PostData(2,  fOutputList);
  
  if(fIsMC){
          PostData(3, fTreePrim);
          PostData(4, fTreeMC);
  }
}

void AliAnalysisTaskStronglyIntensiveCorrTree::Terminate(Option_t* ) {
  //
  fTree = dynamic_cast<TTree*>(GetOutputData(1));
  if (NULL == fTree) {
          AliFatal("NULL == fTree");
          return;
  }
  if(fIsMC){
          fTreePrim = dynamic_cast<TTree*>(GetOutputData(3));
          if (NULL == fTreePrim) {
                    AliFatal("NULL == fTree");
                    return;
          }
          fTreeMC = dynamic_cast<TTree*>(GetOutputData(4));
          if (NULL == fTreeMC) {
                    AliFatal("NULL == fTree");
                    return; // needed to avoid coverity warning
          }
  }
}

TObjArray* AliAnalysisTaskStronglyIntensiveCorrTree::GetAcceptedTracks(AliAODEvent*  pAOD,
							       TClonesArray* arrayMC,
							       Int_t fSelectPrimaryMCDataParticles ) {
  TObjArray* tracks= new TObjArray;
  tracks->SetOwner(kTRUE);

  const Long64_t N(pAOD->GetNumberOfTracks());
  AliDebug(5, Form("#tracks= %6lld", N));
  for (Long64_t i(0); i<N; ++i) {
    AliAODTrack* pAODTrack(dynamic_cast<AliAODTrack*>(pAOD->GetTrack(i)));
    if (NULL == pAODTrack) continue;
   
    if (!pAODTrack->TestFilterBit(fTrackFilter)) continue;
    //Printf("Filterbit %d", fTrackFilter);
   
    if (NULL != arrayMC) {
      //const Int_t label(pAODTrack->GetLabel());
      const Int_t label = TMath::Abs(pAODTrack->GetLabel()); 
      AliAODMCParticle* mcParticle((label >= 0) 
				   ? static_cast<AliAODMCParticle*>(arrayMC->At(label))
				   : NULL);
      if (label >= 0 && NULL == mcParticle)
	AliFatal("MC particle not found");

      switch (fSelectPrimaryMCDataParticles) {
      case -1:
	if (label < 0) continue;
	if (kTRUE  == mcParticle->IsPhysicalPrimary()) continue;
	break;
      case  0:
	
	break;
      case  1:
         // if (label < 0) continue;
	//if (NULL == mcParticle) continue;
	if (kFALSE == mcParticle->IsPhysicalPrimary()) continue;
	if(pAODTrack->GetLabel()<0)Printf("label=%d pAODTrackLabel=%d", label ,pAODTrack->GetLabel());
	break;
      default:
	AliFatal("fSelectPrimaryMCDataParticles != {-1,0,1}");
      }            
    }

    // select only charged tracks
    if (pAODTrack->Charge() == 0) continue;
    Double_t pos[3] = { 0. };
    pAODTrack->GetXYZ(pos);
    if (pAODTrack->Phi() < fPhiMin || pAODTrack->Phi() > fPhiMax) continue;
    if (pAODTrack->Pt()  < fPtMin  || pAODTrack->Pt()  > fPtMax)  continue;

    tracks->Add(new TrackInfoCorr(pAODTrack->Eta(), pAODTrack->Phi(),pAODTrack->Pt(), pAODTrack->XAtDCA(),pAODTrack->YAtDCA(), pAODTrack->ZAtDCA()));
  } // next track
  return tracks;
}

TObjArray* AliAnalysisTaskStronglyIntensiveCorrTree::GetAcceptedTracks(TClonesArray* tracksMC,
							       Int_t fSelectPrimaryMCParticles) {
  // for keeping track of MC labels
  std::set<Int_t> labelSet;

  TObjArray* tracks= new TObjArray;
  tracks->SetOwner(kTRUE);

  const Long64_t N(tracksMC->GetEntries());
  AliDebug(5, Form("#tracks= %6lld", N));
  for (Long64_t i(0); i<N; ++i) {
    AliAODMCParticle* pMCTrack(dynamic_cast<AliAODMCParticle*>(tracksMC->At(i)));
    if (NULL == pMCTrack) continue;    

    // no track filter selection for MC tracks    

    if (labelSet.find(pMCTrack->Label()) != labelSet.end()) {
     // Printf("Duplicate Label= %3d", pMCTrack->Label());
      continue;
    }
    labelSet.insert(pMCTrack->Label());    
    
  // Printf("isPrim = %d (%d)", pMCTrack->IsPhysicalPrimary(), fSelectPrimaryMCParticles);

    switch (fSelectPrimaryMCParticles) {
    case -1:
      if (kTRUE == pMCTrack->IsPhysicalPrimary()) continue;
      break;
    case  0:
      // NOP, take all MC tracks
      break;
    case  1:
      if (kFALSE == pMCTrack->IsPhysicalPrimary()) continue;
      break;
    default:
      AliFatal("fSelectPrimaryMCParticles != {-1,0,1}");
    }

    // select only charged tracks
    if (pMCTrack->Charge() == 0) continue;

    if (pMCTrack->Phi() < fPhiMin || pMCTrack->Phi() > fPhiMax) continue;
    if (pMCTrack->Pt()  < fPtMin  || pMCTrack->Pt()  > fPtMax)  continue;
   
    tracks->Add(new TrackInfoCorr(pMCTrack->Eta(), pMCTrack->Phi(),pMCTrack->Pt(),0, 0, 0));
  } // next track
  return tracks;
}
