// $Id$
//
//
//
//

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDUtils.h"
#include "AliESDInputHandler.h"
#include "AliCentrality.h"
#include "AliESDpid.h"
#include "AliKFParticle.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TParticle.h"


#include "AliESDtrackCuts.h"
#include "AliESDv0.h"
#include "AliV0vertexer.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliVCluster.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "TLorentzVector.h"


#include "AliAnalysisTaskEMCALPhoton.h"
#include "TFile.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEMCALPhoton)

//________________________________________________________________________
AliAnalysisTaskEMCALPhoton::AliAnalysisTaskEMCALPhoton() : 
  AliAnalysisTaskSE(), 
  fTrCuts(0),
  fPrTrCuts(0),
  fSelTracks(0),
  fSelPrimTracks(0),
  fPhotConvArray(0),
  fMyClusts(0),
  fMyAltClusts(0),
  fMyCells(0),
  fMyTracks(0),
  fMyMcParts(0),
  fHeader(0x0),
  fCaloClusters(0),
  fCaloClustersNew(0),
  fEMCalCells(0),
  fGeom(0x0),
  fTimeResTOF(0),
  fMipResponseTPC(0),
  fGeoName("EMCAL_COMPLETEV1"),
  fPeriod("LHC11d"),
  fIsTrain(0),
  fIsMC(0),
  fDebug(0),
  fIsGrid(0),
  fClusThresh(2.0),
  fClusterizer(0),
  fCaloClustersName("EmcalClusterv2"),
  fESD(0),
  fMCEvent(0),
  fStack(0x0),
  fOutputList(0),
  fTree(0),
  fNV0sBefAndAftRerun(0),
  fConversionVtxXY(0),
  fInvMassV0(0),
  fInvMassV0KF(0),
  fInvMassV0SS(0),
  fDedxPAll(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskEMCALPhoton::AliAnalysisTaskEMCALPhoton(const char *name) : 
  AliAnalysisTaskSE(name), 
  fTrCuts(0),
  fPrTrCuts(0),
  fSelTracks(0),
  fSelPrimTracks(0),
  fPhotConvArray(0),
  fMyClusts(0),
  fMyAltClusts(0),
  fMyCells(0),
  fMyTracks(0),
  fMyMcParts(0),
  fHeader(0),
  fCaloClusters(0),
  fCaloClustersNew(0),
  fEMCalCells(0),
  fGeom(0x0),
  fTimeResTOF(0),
  fMipResponseTPC(0),
  fGeoName("EMCAL_COMPLETEV1"),
  fPeriod("LHC11c"),
  fIsTrain(0),
  fIsMC(0),
  fDebug(0),
  fIsGrid(0),
  fClusThresh(2.),
  fClusterizer(0),
  fCaloClustersName("EmcalClusterv2"),
  fESD(0),
  fMCEvent(0),
  fStack(0x0),
  fOutputList(0),
  fTree(0),
  fNV0sBefAndAftRerun(0),
  fConversionVtxXY(0),
  fInvMassV0(0),
  fInvMassV0KF(0),
  fInvMassV0SS(0),
  fDedxPAll(0)
{
  // Constructor
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::UserCreateOutputObjects()
{
  // Create histograms, called once.
    
  fSelTracks = new TObjArray();
  
  fSelPrimTracks = new TObjArray();
  
  if (TClass::GetClass("AliPhotonHeaderObj"))
      TClass::GetClass("AliPhotonHeaderObj")->IgnoreTObjectStreamer();
  fHeader = new AliPhotonHeaderObj;

  if (TClass::GetClass("AliPhotonConvObj"))
      TClass::GetClass("AliPhotonConvObj")->IgnoreTObjectStreamer();
  fPhotConvArray = new TClonesArray("AliPhotonConvObj");
  
  if (TClass::GetClass("AliPhotonClusterObj"))
      TClass::GetClass("AliPhotonClusterObj")->IgnoreTObjectStreamer();
  fMyClusts = new TClonesArray("AliPhotonClusterObj");
  
  if (TClass::GetClass("AliPhotonCellObj"))
      TClass::GetClass("AliPhotonCellObj")->IgnoreTObjectStreamer();
  fMyCells = new TClonesArray("AliPhotonCellObj");
  
  if (TClass::GetClass("AliPhotonTrackObj"))
      TClass::GetClass("AliPhotonTrackObj")->IgnoreTObjectStreamer();
  fMyTracks = new TClonesArray("AliPhotonTrackObj");

  if (fClusterizer || fIsTrain){
    if(fClusterizer)
      fCaloClustersName = fClusterizer->GetNewClusterArrayName();
    else {
      if(fPeriod.Contains("10h") || fPeriod.Contains("11h"))
	fCaloClustersName = "EmcalClustersv1";
      else
	fCaloClustersName = "EmcalClustersv2";
    }
    if (TClass::GetClass("AliPhotonClusterObj"))
	TClass::GetClass("AliPhotonClusterObj")->IgnoreTObjectStreamer();
    fMyAltClusts = new TClonesArray("AliPhotonClusterObj");
  }
  cout<<fCaloClustersName<<endl;
  if(fIsMC){
    if (TClass::GetClass("AliPhotonMcPartObj"))
        TClass::GetClass("AliPhotonMcPartObj")->IgnoreTObjectStreamer();
    fMyMcParts = new TClonesArray("AliPhotonMcPartObj");
  }
 
  fCaloClusters = new TRefArray();
    
  fOutputList = new TList();
  fOutputList->SetOwner();// Container cleans up all histos (avoids leaks in merging) 

  if( !fTree){
    TFile *f = OpenFile(2); 
    TDirectory::TContext context(f);
    if (f) {
      f->SetCompressionLevel(2);
      fTree = new TTree("photon_ana_out", "StandaloneTree");
      fTree->SetDirectory(f);
      if (fIsTrain) {
        fTree->SetAutoFlush(-2*1024*1024);
        fTree->SetAutoSave(0);
      } else {
        fTree->SetAutoFlush(-32*1024*1024);
        fTree->SetAutoSave(0);
      }
    
    fTree->Branch("header", &fHeader, 16*1024, 99);    
    fTree->Branch("conversions", &fPhotConvArray, 8*16*1024, 99);
    fTree->Branch("clusters", &fMyClusts, 8*16*1024, 99);
    if(fClusterizer || fIsTrain)
      fTree->Branch(fCaloClustersName, &fMyAltClusts, 8*16*1024, 99);
    fTree->Branch("cells", &fMyCells, 8*16*1024, 99);
    fTree->Branch("IsoTracks", &fMyTracks, 8*16*1024, 99);
    if(fIsMC)
      fTree->Branch("mcparts", &fMyMcParts, 8*16*1024, 99);
    }
  }
  //if(fIsGrid)fOutputList->Add(fTree);
  fGeom = AliEMCALGeometry::GetInstance(fGeoName);
  
  
  fNV0sBefAndAftRerun = new TH2F("hNV0sBefAndAftRerun","check if the number of v0s change with rerun;old v0 n;new v0 n",50,0.5,50.5,50,0.5,50.5);
  fOutputList->Add(fNV0sBefAndAftRerun);
  
  fConversionVtxXY = new TH2F("hConversionVtxXY","x and y of conversion vertex candidates;x;y",1000,-100,100,1000,-100,100);
  fOutputList->Add(fConversionVtxXY);
  
  fInvMassV0 = new TH1F("hInvMassV0","v0->GetEffMass();v0->GetEffMass();dN/dM",400,0,4);
  fOutputList->Add(fInvMassV0);
  
  fInvMassV0KF = new TH1F("hInvMassV0KF","Inv. mass calculated from AliKFParticle made of V0 tracks;mass_{TrTr};dN/dM",400,0,4);
  fOutputList->Add(fInvMassV0KF);
      
  fInvMassV0SS = new TH1F("hInvMassV0SS","Inv. mass (same sign) calculated from AliKFParticle made of V0 tracks;mass_{TrTr};dN/dM",400,0,4);
  fOutputList->Add(fInvMassV0SS);
      
  fDedxPAll = new TH2F("hDedxPAll","dE/dx vs p (all selected tracks);p (GeV/c);dE/dx (a.u.)",400,0,40, 150,0,150);
  fOutputList->Add(fDedxPAll);
  
  
  PostData(1, fOutputList);
  PostData(2, fTree);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::UserExec(Option_t *) 
{
  // User exec, called once per event.

  Bool_t isSelected = 0;
  if(fPeriod.Contains("11")){
    if(fPeriod.Contains("11a"))
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC1);
    if(fPeriod.Contains("11c") ||fPeriod.Contains("11d") )
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC7);
    if(fPeriod.Contains("11h") )
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAny);//kEMCEGA);

  }
  if(fIsMC){
    isSelected = kTRUE;  
  }
  if(!isSelected)
    return;

  // Post output data.
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available, returning...\n");
    return;
  }
  
  AliESDVertex *pv = (AliESDVertex*)fESD->GetPrimaryVertex();
  if(!pv) 
    return;
  if(TMath::Abs(pv->GetZ())>15)
    return;

  TTree *tree = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetTree();
  TFile *inpfile = (TFile*)tree->GetCurrentFile();

  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track)
      continue;
    
    
    if (fTrCuts && fTrCuts->IsSelected(track)){
      fSelTracks->Add(track);
      fDedxPAll->Fill(track->P(), track->GetTPCsignal());
    }
    if (fPrTrCuts && fPrTrCuts->IsSelected(track))
      fSelPrimTracks->Add(track);
  } //track loop 

  fHeader->fInputFileName  = inpfile->GetName();
  fHeader->fTrClassMask    = fESD->GetHeader()->GetTriggerMask();
  fHeader->fTrCluster      = fESD->GetHeader()->GetTriggerCluster();
  AliCentrality *cent = InputEvent()->GetCentrality();
  fHeader->fV0Cent    = cent->GetCentralityPercentileUnchecked("V0M");
  fHeader->fCl1Cent   = cent->GetCentralityPercentileUnchecked("CL1");
  fHeader->fTrCent    = cent->GetCentralityPercentileUnchecked("TRK");
  if(!fIsTrain){
    for(Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
      if(fGeoName=="EMCAL_FIRSTYEARV1" && mod>3)
        break;
      fGeom->SetMisalMatrix(fESD->GetEMCALMatrix(mod), mod);
    }
  }
  fESD->GetEMCALClusters(fCaloClusters);
  fHeader->fNClus = fCaloClusters->GetEntries();
  fEMCalCells = fESD->GetEMCALCells();
  fHeader->fNCells = fEMCalCells->GetNumberOfCells();
  AliESDtrackCuts *fTrackCuts = new AliESDtrackCuts();
  fHeader->fTrackMult = fTrackCuts->GetReferenceMultiplicity(fESD);//kTrackletsITSTPC ,0.5); 

  fMCEvent = MCEvent();
  if(fMCEvent){
    fStack = (AliStack*)fMCEvent->Stack();
    fHeader->fNMcParts = fStack->GetNtrack();
  }

  
  FindConversions();
  FillMyCells();
  FillMyClusters();
  if(fCaloClustersNew)
    FillMyAltClusters();
  FillIsoTracks();
  if(fIsMC)
    GetMcParts();
  
  fTree->Fill();
  fSelTracks->Clear();
  fSelPrimTracks->Clear();
  fPhotConvArray->Clear();
  fMyClusts->Clear();
  if(fMyAltClusts)
    fMyAltClusts->Clear();
  fMyCells->Clear();
  fMyTracks->Clear();
  if(fMyMcParts)
    fMyMcParts->Clear();
  fCaloClusters->Clear();
  if(fCaloClustersNew)
    fCaloClustersNew->Clear();
  PostData(1, fOutputList);
  PostData(2, fTree);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::FindConversions() 
{
  // Find conversion.

  if(!fTrCuts)
    return;
  Int_t iconv = 0;
  Int_t nV0Orig = fESD->GetNumberOfV0s();
  Int_t nV0s = nV0Orig;
  Int_t ntracks = fESD->GetNumberOfTracks();
  if(!fPeriod.Contains("11h") && !fPeriod.Contains("10h")){
    fESD->ResetV0s();
    AliV0vertexer lV0vtxer;
    lV0vtxer.Tracks2V0vertices(fESD);
    nV0s = fESD->GetNumberOfV0s();
  }
  fNV0sBefAndAftRerun->Fill(nV0Orig, nV0s);
  for(Int_t iv0=0; iv0<nV0s; iv0++){
    AliESDv0 * v0 = fESD->GetV0(iv0);
    if(!v0)
      continue;
    Double_t vpos[3];
    fInvMassV0->Fill(v0->GetEffMass());
    v0->GetXYZ(vpos[0], vpos[1], vpos[2]);
    Int_t ipos = v0->GetPindex();
    Int_t ineg = v0->GetNindex();
    if(ipos<0 || ipos> ntracks)
      continue;
    if(ineg<0 || ineg> ntracks)
      continue;
    AliESDtrack *pos = static_cast<AliESDtrack*>(fESD->GetTrack(ipos));
    AliESDtrack *neg = static_cast<AliESDtrack*>(fESD->GetTrack(ineg));
    if(!pos || !neg)
      continue;
    if (!fTrCuts->IsSelected(pos) || !fTrCuts->IsSelected(neg))
      continue;
    /*if(pos->GetTPCsignal()<65 || neg->GetTPCsignal()<65)
      continue;*/
    const AliExternalTrackParam * paramPos = v0->GetParamP() ;
    const AliExternalTrackParam * paramNeg = v0->GetParamN() ;
    if(!paramPos || !paramNeg)
      continue;
    if(pos->GetSign() <0){//change tracks
      pos=neg ;
      neg=fESD->GetTrack(v0->GetPindex()) ;
      paramPos=paramNeg ;
      paramNeg=v0->GetParamP() ;
    }
    AliKFParticle negKF(*paramNeg,11);
    AliKFParticle posKF(*paramPos,-11);
    AliKFParticle photKF(negKF,posKF) ;
   
    if( neg->GetKinkIndex(0) > 0 || pos->GetKinkIndex(0) > 0) 
      continue ;
    
    if(pos->GetSign() == neg->GetSign()){ 
      fInvMassV0SS->Fill(photKF.GetMass());
      continue;
    }
    fInvMassV0KF->Fill(photKF.GetMass());
    TLorentzVector photLV(photKF.GetPx(), photKF.GetPy(),photKF.GetPz(), photKF.GetE()); 
    if(photLV.M()>0.05 || photLV.M()<0)
      continue;
    fConversionVtxXY->Fill(vpos[0], vpos[1]);//negKF.GetX(), negKF.GetY());
    Double_t convPhi = TMath::ATan2(photLV.Py(),photLV.Px());
    if(convPhi<0)
      convPhi+=TMath::TwoPi();
    TVector3 vecpos(vpos);
    Double_t v0Phi = 0;
    if(vecpos.Perp()>0)
      vecpos.Phi();
    if(v0Phi<0)
      v0Phi+=TMath::TwoPi();
    AliPhotonConvObj *myconv = static_cast<AliPhotonConvObj*>(fPhotConvArray->New(iconv++));
    myconv->fPt = photLV.Pt();
    myconv->fEta = photLV.Eta();
    myconv->fPhi = convPhi;
    myconv->fVR = vecpos.Perp();
    if(vecpos.Perp()>0)
      myconv->fVEta = vecpos.Eta();
    myconv->fVPhi = v0Phi;
    myconv->fMass = photLV.M();
    myconv->fMcLabel = -3; //WARNING: include the correct labeling
    //negative daughter
   myconv->fNegPt      =  negKF.GetPt();
   myconv->fNegEta     =  negKF.GetEta();
   Double_t trackPhi   =  negKF.GetPhi();
   if(trackPhi<0)
     trackPhi+=TMath::TwoPi();
   myconv->fNegPhi     =  trackPhi;
   myconv->fNegDedx    =  neg->GetTPCsignal();
   myconv->fNegMcLabel =  neg->GetLabel();
    //negative daughter
   myconv->fPosPt      =  posKF.GetPt();
   myconv->fPosEta     =  posKF.GetEta();
   trackPhi            =  posKF.GetPhi();
   if(trackPhi<0)
     trackPhi+=TMath::TwoPi();
   myconv->fPosPhi     =  trackPhi;
   myconv->fPosDedx    =  pos->GetTPCsignal();
   myconv->fPosMcLabel =  pos->GetLabel();
  }
  return;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::FillMyCells() 
{
  // Fill cells.

  if (!fEMCalCells)
    return;
  Int_t ncells = fEMCalCells->GetNumberOfCells();
  Int_t mcel = 0;
  for(Int_t icell = 0; icell<ncells; icell++){
    Int_t absID = TMath::Abs(fEMCalCells->GetCellNumber(icell));
    AliPhotonCellObj *mycell = static_cast<AliPhotonCellObj*>(fMyCells->New(mcel++));
    Float_t eta=-1, phi=-1;
    if(!fGeom){
      std::cout<<"+++fGeom not found! MyCells branch will not be filled for this event!+++\n";
      return;
    }
    if(!fGeom)
      return;
    /*if(!fIsMC)*/fGeom->EtaPhiFromIndex(absID,eta,phi);
    Float_t theta = 2*TMath::ATan(TMath::Exp(-eta));
    mycell->fAbsID = absID;
    mycell->fE = fEMCalCells->GetCellAmplitude(absID);
    mycell->fEt = fEMCalCells->GetCellAmplitude(absID)*TMath::Sin(theta);
    mycell->fEta = eta;
    mycell->fPhi = phi;
    mycell->fTime = fEMCalCells->GetCellTime(absID);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::FillMyClusters() 
{
  // Fill clusters.

  if (!fCaloClusters)
    return;
  Int_t nclus = fCaloClusters->GetEntries();
  Int_t mcl = 0;
  for(Int_t ic=0; ic < nclus; ic++){
    AliESDCaloCluster *clus = static_cast<AliESDCaloCluster*>(fCaloClusters->At(ic));
    if(!clus)
      continue;
    if(!clus->IsEMCAL())
      continue;
    if(clus->E() < fClusThresh)
      continue;
    Float_t pos[3];
    clus->GetPosition(pos);
    TVector3 cpos(pos);
    TString cellsAbsID;
    AliPhotonClusterObj *myclus = static_cast<AliPhotonClusterObj*>(fMyClusts->New(mcl++));
    myclus->fE       = clus->E();
    myclus->fEt      = clus->E()*TMath::Sin(cpos.Theta());
    myclus->fR       = cpos.Perp();
    myclus->fEta     = cpos.Eta();
    myclus->fPhi     = cpos.Phi();
    if(cpos.Phi()<0){
      myclus->fPhi+=TMath::TwoPi();
    }
    myclus->fN       = clus->GetNCells();
    Short_t  id = -1;
    myclus->fEmax    = GetMaxCellEnergy( clus, id); 
    myclus->fIdmax   = id;
    myclus->fTmax    = fEMCalCells->GetCellTime(id);
    myclus->fEcross  = GetCrossEnergy( clus, id);
    myclus->fDisp    = clus->GetDispersion();
    myclus->fM20     = clus->GetM20();
    myclus->fM02     = clus->GetM02();
    myclus->fTrDEta  = clus->GetTrackDz();
    myclus->fTrDPhi  = clus->GetTrackDx();
    myclus->fTrIso01 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.1, 0.);
    myclus->fTrIso02 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.2, 0.);
    myclus->fTrIso03 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.3, 0.);
    myclus->fTrIso04 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.4, 0.);
    myclus->fTrPhiBand01 = GetPhiBandEt( cpos.Eta(), cpos.Phi(), 0.1, 0.);
    myclus->fTrPhiBand02 = GetPhiBandEt( cpos.Eta(), cpos.Phi(), 0.2, 0.);
    myclus->fTrPhiBand03 = GetPhiBandEt( cpos.Eta(), cpos.Phi(), 0.3, 0.);
    myclus->fTrPhiBand04 = GetPhiBandEt( cpos.Eta(), cpos.Phi(), 0.4, 0.);
    for(Int_t icell=0;icell<myclus->fN;icell++){
      Int_t absID = clus->GetCellAbsId(icell);
      cellsAbsID.Append(Form("%d",absID));
      cellsAbsID.Append(";");
    }
    myclus->fCellsAbsId = cellsAbsID;
    myclus->fMcLabel = clus->GetLabel(); 
    Int_t matchIndex = clus->GetTrackMatchedIndex();
    if(matchIndex<0 || matchIndex>fESD->GetNumberOfTracks()){
      myclus->fTrEp = -1;
      continue;
    }
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(matchIndex));
    if(!track){
      myclus->fTrEp = -1;
      continue;
    }
    if(!fPrTrCuts){
      myclus->fTrEp = -1;
      continue;
    }
    if(!fPrTrCuts->IsSelected(track)){
      myclus->fTrEp = -1;
      continue;
    }
    myclus->fTrEp = clus->E()/track->P();
  }
  
}
//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::FillMyAltClusters() 
{
  // Fill clusters.

  if(!fCaloClustersNew)
    return;
  Int_t nclus = fCaloClustersNew->GetEntries();
  Int_t mcl = 0;
  for(Int_t ic=0; ic < nclus; ic++){
    AliESDCaloCluster *clus = static_cast<AliESDCaloCluster*>(fCaloClustersNew->At(ic));
    if(!clus)
      continue;
    if(!clus->IsEMCAL())
      continue;
    if(clus->E() < fClusThresh)
      continue;
    Float_t pos[3];
    clus->GetPosition(pos);
    TVector3 cpos(pos);
    TString cellsAbsID;
    AliPhotonClusterObj *myclus = static_cast<AliPhotonClusterObj*>(fMyAltClusts->New(mcl++));
    myclus->fE       = clus->E();
    myclus->fEt      = clus->E()*TMath::Sin(cpos.Theta());
    myclus->fR       = cpos.Perp();
    myclus->fEta     = cpos.Eta();
    myclus->fPhi     = cpos.Phi();
    if(cpos.Phi()<0){
      myclus->fPhi+=TMath::TwoPi();
    }
    myclus->fN       = clus->GetNCells();
    myclus->fDisp    = clus->GetDispersion();
    myclus->fM20     = clus->GetM20();
    myclus->fM02     = clus->GetM02();
    myclus->fTrDEta  = clus->GetTrackDz();
    myclus->fTrDPhi  = clus->GetTrackDx();
    myclus->fTrIso01 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.1, 0.);
    myclus->fTrIso02 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.2, 0.);
    myclus->fTrIso03 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.3, 0.);
    myclus->fTrIso04 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.4, 0.);
    for(Int_t icell=0;icell<myclus->fN;icell++){
      Int_t absID = clus->GetCellAbsId(icell);
      cellsAbsID.Append(Form("%d",absID));
      cellsAbsID.Append(";");
    }
    myclus->fCellsAbsId = cellsAbsID;
    myclus->fMcLabel = clus->GetLabel(); 
    Int_t matchIndex = clus->GetTrackMatchedIndex();
    if(matchIndex<0 || matchIndex>fESD->GetNumberOfTracks()){
      myclus->fTrEp = -1;
      continue;
    }
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(matchIndex));
    if(!track){
      myclus->fTrEp = -1;
      continue;
    }
    if(!fPrTrCuts){
      myclus->fTrEp = -1;
      continue;
    }
    if(!fPrTrCuts->IsSelected(track)){
      myclus->fTrEp = -1;
      continue;
    }
    myclus->fTrEp = clus->E()/track->P();
  }
  
}
//________________________________________________________________________
void  AliAnalysisTaskEMCALPhoton::FillIsoTracks()
{
  // Fill high pt tracks.

  if(!fSelPrimTracks)
    return;
  Int_t ntracks = fSelPrimTracks->GetEntries();
  Int_t imtr = 0;
  for(Int_t it=0;it<ntracks; it++){
    AliESDtrack *track = static_cast<AliESDtrack*>(fSelPrimTracks->At(it));
    if(!track)
      continue;
    if(track->Phi()<1.0 || track->Phi()>3.55)
      continue;
    if(TMath::Abs(track->Eta())>1.1)
      continue;
    /*if(track->Pt()<3)
      continue;*/
    AliPhotonTrackObj *mtr = static_cast<AliPhotonTrackObj*>(fMyTracks->New(imtr++));
    mtr->fPt = track->Pt();
    mtr->fEta = track->Eta();
    mtr->fPhi = track->Phi();
    mtr->fCharge = track->Charge();
    mtr->fDedx = track->GetTPCsignal();
    mtr->fMcLabel = track->GetLabel();
  }
}

//________________________________________________________________________
void  AliAnalysisTaskEMCALPhoton::GetMcParts()
{
  // Get MC particles.

  if (!fStack)
    return;

  const AliVVertex *evtVtx = fMCEvent->GetPrimaryVertex();
  if (!evtVtx)
    return;
  Int_t nTracks = fStack->GetNtrack();
  AliPhotonMcPartObj *mcp = 0;
  for (Int_t iTrack = 0; iTrack<nTracks; ++iTrack) {
    TParticle *mcP = static_cast<TParticle*>(fStack->Particle(iTrack));
    if (!mcP){
      mcp = static_cast<AliPhotonMcPartObj*>(fMyMcParts->New(iTrack));
      mcp->fLabel = iTrack;
      continue;
    }
    // primary particle
    Double_t dR = TMath::Sqrt((mcP->Vx()-evtVtx->GetX())*(mcP->Vx()-evtVtx->GetX()) + 
                              (mcP->Vy()-evtVtx->GetY())*(mcP->Vy()-evtVtx->GetY()) +
                              (mcP->Vz()-evtVtx->GetZ())*(mcP->Vz()-evtVtx->GetZ()));
    if(dR > 0.5){
      mcp = static_cast<AliPhotonMcPartObj*>(fMyMcParts->New(iTrack));
      mcp->fLabel = iTrack;
      continue;
    }
    
    // kinematic cuts
    Double_t pt = mcP->Pt() ;
    if (pt<0.5){
      mcp = static_cast<AliPhotonMcPartObj*>(fMyMcParts->New(iTrack));
      mcp->fLabel = iTrack;
      continue;
    }
    Double_t eta = mcP->Eta();
    if (TMath::Abs(eta)>0.7){
      mcp = static_cast<AliPhotonMcPartObj*>(fMyMcParts->New(iTrack));
      mcp->fLabel = iTrack;
      continue;
    }
    Double_t phi  = mcP->Phi();
    if (phi<1.0||phi>3.3){
      mcp->fLabel = iTrack;
      mcp = static_cast<AliPhotonMcPartObj*>(fMyMcParts->New(iTrack));
      continue;
    }
    // pion or eta meson or direct photon
    if(mcP->GetPdgCode() == 111) {
    } else if(mcP->GetPdgCode() == 221) {
    } else if(mcP->GetPdgCode() == 22 ) {
    } else {
      mcp = static_cast<AliPhotonMcPartObj*>(fMyMcParts->New(iTrack));
      mcp->fLabel = iTrack;
      continue;
    }
    
    bool checkIfAlreadySaved = false;
    for(Int_t imy=0;imy<fMyMcParts->GetEntries();imy++){
      AliPhotonMcPartObj *mymc = static_cast<AliPhotonMcPartObj*>(fMyMcParts->At(imy));
      if(!mymc)
	continue;
      if(imy == iTrack)
	checkIfAlreadySaved = true;
    }
    if(!checkIfAlreadySaved){
      FillMcPart(mcP,  iTrack);
    }
    for(Int_t id=mcP->GetFirstDaughter(); id <= mcP->GetLastDaughter(); id++){
      if(id<=mcP->GetMother(0))
	continue;
      if(id<0 || id>nTracks)
	continue;
      TParticle *mcD = static_cast<TParticle*>(fStack->Particle(id));
      if(!mcD)
	continue;
      FillMcPart(mcD, id);
      }
  }
}

//________________________________________________________________________
void  AliAnalysisTaskEMCALPhoton::FillMcPart(TParticle *mcP,  Int_t itrack)
{
  // Fill MC particles.

  if(!mcP)
    return;
  TVector3 vmcv(mcP->Vx(),mcP->Vy(), mcP->Vz());
  AliPhotonMcPartObj *mcp = static_cast<AliPhotonMcPartObj*>(fMyMcParts->New(itrack));
  mcp->fLabel = itrack ;
  mcp->fPdg = mcP->GetPdgCode() ;
  mcp->fPt = mcP->Pt() ;
  mcp->fEta = mcP->Eta() ;
  mcp->fPhi = mcP->Phi() ;
  mcp->fVR = vmcv.Perp();
  if(vmcv.Perp()>0){
    mcp->fVEta = vmcv.Eta() ;
    mcp->fVPhi = vmcv.Phi() ;
  }
  mcp->fMother = mcP->GetMother(0) ;
  mcp->fFirstD = mcP->GetFirstDaughter() ;
  mcp->fLastD = mcP->GetLastDaughter() ;
  mcp->fStatus = mcP->GetStatusCode();
  mcp->fIso = AliAnalysisTaskEMCALPhoton::GetMcIsolation(mcP, itrack, 0.4 , 0.2);
}
//________________________________________________________________________                                                                                                                                   
Double_t AliAnalysisTaskEMCALPhoton::GetMcIsolation(TParticle *mcP, Int_t itrack, Double_t radius, Double_t pt) const
{
  if (!fStack)
    return -1;
  if(!mcP)
    return -1;
  if(itrack<6 || itrack>8)
    return -1;
  if(mcP->GetPdgCode()!=22)
    return -1;
  Double_t sumpt=0;
  Float_t eta = mcP->Eta();
  Float_t phi = mcP->Phi();
  Float_t dR;
  Int_t nparts =  fStack->GetNtrack();
  for(Int_t ip = 0; ip<nparts; ip++){
    TParticle *mcisop = static_cast<TParticle*>(fStack->Particle(ip));
    if(!mcisop)
      continue;
    if(ip==itrack)
      continue;
    if(mcisop->GetStatusCode()!=1)
      continue;
    if(mcisop->Pt()<pt)
      continue;
    TVector3 vmcv(mcisop->Vx(),mcisop->Vy(), mcisop->Vz());  
    if(vmcv.Perp()>1)
      continue;
    dR = TMath::Sqrt((phi-mcisop->Phi())*(phi-mcisop->Phi())+(eta-mcisop->Eta())*(eta-mcisop->Eta()));
    if(dR>radius)
      continue;
    sumpt += mcisop->Pt();
  }
  return sumpt;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPhoton::GetTrackIsolation(Double_t cEta, Double_t cPhi, Double_t radius, Double_t pt) const
{
  // Compute isolation based on tracks.
  
  Double_t trkIsolation = 0;
  Double_t rad2 = radius*radius;
  Int_t ntrks = fSelPrimTracks->GetEntries();
  for(Int_t j = 0; j<ntrks; ++j) {
    AliVTrack *track = static_cast<AliVTrack*>(fSelPrimTracks->At(j));
    if (!track)
      continue;
    if (track->Pt()<pt)
      continue;
    
    Float_t eta = track->Eta();
    Float_t phi = track->Phi();
    Double_t phidiff = TVector2::Phi_mpi_pi(phi-cPhi);
    Double_t dist = (eta-cEta)*(eta-cEta)+phidiff*phidiff;
    if(dist>rad2)
      continue;
    trkIsolation += track->Pt();
  } 
  return trkIsolation;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPhoton::GetPhiBandEt(Double_t eta, Double_t phi, Double_t R, Double_t minpt) const
{
  // Get phi band.

  if(!fSelPrimTracks)
    return 0;
  Int_t ntracks = fSelPrimTracks->GetEntries();
  Double_t loweta = eta - R;
  Double_t upeta = eta + R;
  Double_t upphi = phi + R;
  Double_t lowphi = phi - R;
  Double_t et = 0;
  for(Int_t itr=0; itr<ntracks; itr++){
    AliESDtrack *track = static_cast<AliESDtrack*>(fSelPrimTracks->At(itr));
    if(!track)
      continue;
    if(track->Pt()<minpt)
      continue;
    if((track->Eta() < upeta) && (track->Eta() > loweta))
      continue;
    if((track->Phi() > upphi) || (track->Phi() < lowphi))
      continue;
    et+=track->Pt();
  }
  return et;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPhoton::GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax)
{
  // Calculate the energy of cross cells around the leading cell.

  AliVCaloCells *cells = 0;
  cells = fESD->GetEMCALCells();
  if (!cells)
    return 0;

  if (!fGeom)
    return 0;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Int_t iphis   = -1;
  Int_t ietas   = -1;

  Double_t crossEnergy = 0;

  fGeom->GetCellIndex(idmax,iSupMod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphis,ietas);

  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Int_t cellAbsId = cluster->GetCellAbsId(i);
    fGeom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
    Int_t aphidiff = TMath::Abs(iphi-iphis);
    if (aphidiff>1)
      continue;
    Int_t aetadiff = TMath::Abs(ieta-ietas);
    if (aetadiff>1)
      continue;
    if ( (aphidiff==1 && aetadiff==0) ||
	(aphidiff==0 && aetadiff==1) ) {
      crossEnergy += cells->GetCellAmplitude(cellAbsId);
    }
  }

  return crossEnergy;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPhoton ::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
{
  // Get maximum energy of attached cell.

  id = -1;

  AliVCaloCells *cells = 0;
  cells = fESD->GetEMCALCells();
  if (!cells)
    return 0;

  Double_t maxe = 0;
  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Double_t e = cells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
    if (e>maxe) {
      maxe = e;
      id   = cluster->GetCellAbsId(i);
    }
  }
  return maxe;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::Terminate(Option_t *) 
{
  // Called once at the end of the query
/*  if(fIsGrid)
    return;*/
  if (fTree) {
    printf("***tree %s being saved***\n",fTree->GetName());
    TFile *f = OpenFile(2);
    TDirectory::TContext context(f);
    if (f) 
      fTree->Write();
  }
}
