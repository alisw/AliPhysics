#include "TF1.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TMath.h"
#include "THashList.h"
#include "TChain.h"

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"

#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliCaloPhoton.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliPHOSGeometry.h"
#include "AliOADBContainer.h"
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskPHOSObjectCreator.h"

ClassImp(AliAnalysisTaskPHOSObjectCreator)

//________________________________________________________________________
AliAnalysisTaskPHOSObjectCreator::AliAnalysisTaskPHOSObjectCreator(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(0x0),
  fHistoMggvsEProbe(0x0),
  fHistoMggvsEPassingProbe(0x0),
  fEvent(0x0),
  fAODEvent(0x0),
  fESDEvent(0x0),
  fPHOSObjectArray(NULL),
  fPHOSGeo(0x0),
  fNonLinCorr(0x0),
  fUserNonLinCorr(0x0),
  fRunNumber(0),
  fUsePHOSTender(kTRUE),
  fIsMC(kFALSE),
  fBunchSpace(25.),
  fMCArrayESD(0x0),
  fMCArrayAOD(0x0),
  fIsM4Excluded(kTRUE),
  fIsSingleSim(kFALSE),
  fIsEmbedding(kFALSE)
{
  // Constructor
  for(Int_t i=0;i<3;i++){
    fVertex[i] = 0;
  }
 
  for(Int_t i=0;i<6;i++){
    fPHOSBadMap[i] = 0;
  }

  //Initialize non-linrarity correction
  fNonLinCorr = new TF1("nonlin","0.0241 + 1.0504*x + 0.000249*x*x",0.,100.);
  fUserNonLinCorr = new TF1("usernonlin","1.",0.,100.);
 
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, THashList::Class());

}
//________________________________________________________________________
AliAnalysisTaskPHOSObjectCreator::~AliAnalysisTaskPHOSObjectCreator()
{
  //destructor

  if(fPHOSObjectArray){
    delete fPHOSObjectArray;
    fPHOSObjectArray = 0x0;
  }

  if(fNonLinCorr){
    delete fNonLinCorr;
    fNonLinCorr = 0x0;
  }

  if(fUserNonLinCorr){
    delete fUserNonLinCorr;
    fUserNonLinCorr = 0x0;
  }

}
//________________________________________________________________________
void AliAnalysisTaskPHOSObjectCreator::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  const Int_t NpTgg = 101;
  Double_t pTgg[NpTgg]={};

  for(Int_t i=0;i<50;i++)     pTgg[i] = 0.1 * i;            //every 0.1 GeV/c, up to 5 GeV/c
  for(Int_t i=50;i<60;i++)    pTgg[i] = 0.5 * (i-50) + 5.0; //every 0.5 GeV/c, up to 10 GeV/c
  for(Int_t i=60;i<NpTgg;i++) pTgg[i] = 1.0 * (i-60) + 10.0;//every 1.0 GeV/c, up to 50 GeV/c

  fHistoMggvsEProbe = new TH2F("hMgg_STDCut_Probe","Probe #gamma for standard cluster cut",180,0,0.72,NpTgg-1,pTgg);
  fHistoMggvsEProbe->Sumw2();
  fHistoMggvsEProbe->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
  fHistoMggvsEProbe->SetYTitle("E_{#gamma} (GeV)");

  fHistoMggvsEPassingProbe = new TH2F("hMgg_STDCut_PassingProbe","Passing Probe #gamma for standard cluster cut",180,0,0.72,NpTgg-1,pTgg);
  fHistoMggvsEPassingProbe->Sumw2();
  fHistoMggvsEPassingProbe->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
  fHistoMggvsEPassingProbe->SetYTitle("E_{#gamma} (GeV)");

  fOutputContainer->Add(fHistoMggvsEProbe);
  fOutputContainer->Add(fHistoMggvsEPassingProbe);

  Init();//initialize PHOS object array

  PostData(1,fOutputContainer);
}
//________________________________________________________________________
void AliAnalysisTaskPHOSObjectCreator::UserExec(Option_t *option) 
{
  // Main loop
  // Called for each event

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fEvent){
    AliError("event is not available.");
    return;
  }

  //cout << "Object Creator usePHOSTender = " << fUsePHOSTender << " , fIsMC = " << fIsMC << endl;

  fESDEvent = dynamic_cast<AliESDEvent*>(fEvent);
  fAODEvent = dynamic_cast<AliAODEvent*>(fEvent);

  if(fIsMC){
    fMCArrayESD = 0x0;
    fMCArrayAOD = 0x0;
    if(fESDEvent){
      AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
      if(eventHandler){
        AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*> (eventHandler);
        if(mcEventHandler) fMCArrayESD = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack();
      }

      if(!fMCArrayESD) AliError("Could not get MC Stack!");

    }
    else if(fAODEvent){
      fMCArrayAOD = dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!fMCArrayAOD){
        AliError("Could not retrieve MC array!");
        return;
      }
    }

  }

  if(fRunNumber != fEvent->GetRunNumber()) { // Check run number
    fRunNumber = fEvent->GetRunNumber();
    SetGeometry();
  }

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  fVertex[0] = vVertex->GetX();
  fVertex[1] = vVertex->GetY();
  fVertex[2] = vVertex->GetZ();

  AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fEvent->GetPHOSCells());
  Int_t multClust = fEvent->GetNumberOfCaloClusters();

  fPHOSObjectArray->Clear();

  TLorentzVector p1,p1core;

  Double_t energy=0, tof=-999, M20=0, M02=0, R2=999, coreE=0, coreR2=999;
  Double_t TrackDx=0, TrackDz=0, TrackPt=0, r=999;
  Int_t TrackCharge = 0;
  Int_t trackindex=-1;

  Int_t relId[4] = {};
  Int_t module=-1, cellx=-1, cellz=-1;
  Int_t digMult=0;
  Float_t position[3] = {}; 
  AliVTrack *vtrack = 0x0;
  AliVCluster *cluster = 0x0;
  Double_t distance=0;
  Int_t inPHOS=0;

  for(Int_t iclu=0; iclu<multClust; iclu++){

    cluster = (AliVCluster*)fEvent->GetCaloCluster(iclu);

    if(cluster->GetType() != AliVCluster::kPHOSNeutral) continue;
    if(cluster->E() < 0.1) continue;//energy is set to 0 GeV in PHOS Tender, if its position is one th bad channel.//0.05 GeV is threshold of seed in a cluster by clustering algorithm.


    //printf("energy = %e , coreE = %e\n",cluster->E(),cluster->GetCoreEnergy());

    distance = cluster->GetDistanceToBadChannel();//in cm.

    //cout << "Ncell before = " << cluster->GetNCells() << endl;

    Int_t Ncell = cluster->GetNCells();
    for(Int_t i=0;i<cluster->GetNCells();i++){
      Int_t absId_tmp = cluster->GetCellAbsId(i);
      Double_t amp = cells->GetCellAmplitude(absId_tmp);

      if(amp < 2e-6) //less than 2 keV
        Ncell--;
    }
    //cluster->SetNCells(Ncell);

    cluster->GetPosition(position);
    TVector3 global1(position);
    fPHOSGeo->GlobalPos2RelId(global1,relId);

    module = relId[0];
    cellx  = relId[2];
    cellz  = relId[3];

    if(fIsM4Excluded && module==4) continue;

    if(module < 1 || 4 < module){
      AliError(Form("Wrong module number %d",module));
      return;
    }

    if(!fUsePHOSTender && !IsGoodChannel("PHOS",module,cellx,cellz)) continue;

    energy = cluster->E();
    digMult = cluster->GetNCells();
    tof = cluster->GetTOF();

    if(fUsePHOSTender){
      //When PHOSTender is applied, GetM20(), GetM02() returns M20,M02 of core of cluster respectively.
      M20 = cluster->GetM20();//M20 is short axis of elliptic shower shape.
      M02 = cluster->GetM02();//M02 is long  axis of elliptic shower shape.
      R2 = cluster->GetDispersion();//full dispersion
      coreE = cluster->GetCoreEnergy();
      coreR2 = cluster->Chi2();//core dispersion
      r = cluster->GetEmcCpvDistance();
    }
    else{//no tender
      //When PHOSTender is NOT applied, GetM20(), GetM02() returns M20,M02 of full cluster respectively.
      M20 = cluster->GetM20();//M20 is short axis of elliptic shower shape.
      M02 = cluster->GetM02();//M02 is long  axis of elliptic shower shape.
      R2 = TestLambda(energy,M20,M02);
      coreE = CoreEnergy(cluster,cells);
      EvalCoreLambdas(cluster,cells,M02,M20);
      coreR2 = TestCoreLambda(energy,M20,M02);

      vtrack = 0x0;
      if(fESDEvent){//for ESD
        trackindex = cluster->GetTrackMatchedIndex();
        if(trackindex > 0){
          vtrack = (AliVTrack*)(fEvent->GetTrack(trackindex));
        }//end of track matching
      }//end of track selection in ESD
      else if(fAODEvent){//for AOD
        if(cluster->GetNTracksMatched() > 0){
          vtrack = dynamic_cast<AliVTrack*>(cluster->GetTrackMatched(0));
        }//end of track matching
      }//end of AOD

      if(vtrack){
        TrackDx = cluster->GetTrackDx();
        TrackDz = cluster->GetTrackDz();
        TrackPt = vtrack->Pt();
        TrackCharge = vtrack->Charge();
        r = TestCPVRun2(TrackDx,TrackDz,TrackPt,TrackCharge);

      }//end of track matching
      else r = 999;//no matched track
    }

    cluster->GetMomentum(p1,fVertex);
    cluster->GetMomentum(p1core,fVertex);

    p1 *= fUserNonLinCorr->Eval(p1.E());
    p1core *= coreE/energy * fUserNonLinCorr->Eval(coreE); //use core energy in PbPb.

    if(p1.E() < 0.1) continue;//minimum energy cut after NL correciton. Note Ecore <= Efull

    new((*fPHOSObjectArray)[inPHOS]) AliCaloPhoton(p1.Px(),p1.Py(),p1.Pz(),p1.E());
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSObjectArray->At(inPHOS); 
    ph->SetCluster(cluster);
    ph->SetModule(module);
    ph->SetNCells(digMult); 
    ph->SetTime(tof);//unit of second
    ph->SetTOFBit(TMath::Abs(tof*1e+9) < fBunchSpace/2.);
    ph->SetDistToBadfp(distance/2.2);//in unit of cells with floating point. 2.2 cm is crystal size
    ph->SetEMCx((Double_t)position[0]);
    ph->SetEMCy((Double_t)position[1]);
    ph->SetEMCz((Double_t)position[2]);
    ph->SetWeight(1.);

    if(fIsMC){
      Bool_t sure = kTRUE;
      Int_t label = FindPrimary(ph,sure);
      ph->SetPrimary(label);
      //ph->SetPrimary(cluster->GetLabel());
    }

    ph->SetLambdas(M20,M02);
    ph->SetNsigmaCPV(r);
    ph->SetNsigmaFullDisp(TMath::Sqrt(R2));
    ph->SetNsigmaCoreDisp(TMath::Sqrt(coreR2));
    ph->SetMomV2(&p1core);//core energy

    //printf("energy = %e , coreE = %e\n",ph->Energy(),ph->GetMomV2()->Energy());

    inPHOS++;
  }

  EstimateSTDCutEfficiency(fPHOSObjectArray);

  const Int_t Nph = fPHOSObjectArray->GetEntries();

  for(Int_t iph=0;iph<Nph;iph++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSObjectArray->At(iph);
    AliVCluster *cluster = (AliVCluster*)ph->GetCluster();
    if(!PassSTDCut(cluster)){
      //AliInfo("cluster is rejected by the standard cluster cut.");
      fPHOSObjectArray->Remove(ph);
    }
  }
  fPHOSObjectArray->Compress();

  AliVEvent* input = InputEvent();
  TObject* outO = input->FindListObject("PHOSClusterArray");

  if(!outO){
    fPHOSObjectArray->SetName("PHOSClusterArray");
    input->AddObject(fPHOSObjectArray);
  }

}
//________________________________________________________________________
void AliAnalysisTaskPHOSObjectCreator::Terminate(Option_t *option) 
{
  //Called once at the end of the query
  //fUserNonLinCorr->Draw();

}
//________________________________________________________________________
void AliAnalysisTaskPHOSObjectCreator::Init() 
{
  //Called once at the end of the query
 
  //if(fPHOSObjectArray != NULL){
  if(fPHOSObjectArray){
    //delete fPHOSObjectArray;
    fPHOSObjectArray = NULL;
    //fPHOSObjectArray->Clear();
  }

  fPHOSObjectArray = new TClonesArray("AliCaloPhoton",200);
  fPHOSObjectArray->Delete();

}
//________________________________________________________________________
void AliAnalysisTaskPHOSObjectCreator::InitBadMap() 
{

  AliOADBContainer badmapContainer(Form("phosBadMap"));
  badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root","phosBadMap");
  TObjArray *maps = (TObjArray*)badmapContainer.GetObject(fRunNumber,"phosBadMap");
  if(!maps){
    AliError(Form("Can not read Bad map for run %d. \n You may choose to use your map with ForceUsingBadMap().",fRunNumber)) ;
  }
  else{
    AliInfo(Form("Setting PHOS bad map with name %s.",maps->GetName())) ;
    for(Int_t mod=0; mod<6;mod++){
      if(fPHOSBadMap[mod]) delete fPHOSBadMap[mod];
      TH2I * h = (TH2I*)maps->At(mod);
      if(h) fPHOSBadMap[mod]=new TH2I(*h);
    }
  }

}
//________________________________________________________________________
void AliAnalysisTaskPHOSObjectCreator::SetGeometry()
{
  // Initialize the PHOS geometry

  if(fUsePHOSTender){

    if(fRunNumber < 209122)//Run1
      fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;
    else//Run2
      fPHOSGeo = AliPHOSGeometry::GetInstance("Run2") ;
  }
  else{//PHOSTender is not applied.
    AliOADBContainer geomContainer("phosGeo");

    if(fRunNumber < 209122)//Run1
      fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;
    else//Run2
      fPHOSGeo = AliPHOSGeometry::GetInstance("Run2") ;

    if(fIsMC){ //use excatly the same geometry as in simulation, stored in esd
      if(fESDEvent){
        for(Int_t mod=0; mod<6; mod++) {
          const TGeoHMatrix * m = fESDEvent->GetPHOSMatrix(mod);
          if(m){
            fPHOSGeo->SetMisalMatrix(m,mod) ;
            printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo);
            m->Print();
          }
        }
      }
      else if(fAODEvent){ //To be fixed
        geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSMCGeometry.root","PHOSMCRotationMatrixes");
        TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
        for(Int_t mod=0; mod<6; mod++){
          if(!matrixes->At(mod)) continue;
          fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod);
          printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo);
          ((TGeoHMatrix*)matrixes->At(mod))->Print();
        }
      }
    }//end of MC
    else{ //Use best approaximation to real geometry
      geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
      TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
      for(Int_t mod=0; mod<6; mod++){
        if(!matrixes->At(mod)) continue;
        fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod);
        printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo);
        ((TGeoHMatrix*)matrixes->At(mod))->Print();
      }
    }//end of real data
  }

}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSObjectCreator::IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz)
{
  //Check if this channel belogs to the good ones

  if(strcmp(det,"PHOS")==0){
    if(mod>5 || mod<1){
      AliError(Form("No bad map for PHOS module %d ",mod)) ;
      return kTRUE ;
    }
    if(!fPHOSBadMap[mod]){
      AliError(Form("No Bad map for PHOS module %d",mod)) ;
      return kTRUE ;
    }
    if(fPHOSBadMap[mod]->GetBinContent(ix,iz)>0)
      return kFALSE ;
    else
      return kTRUE ;
  }
  else{
    AliError(Form("Can not find bad channels for detector %s ",det)) ;
  }
  return kTRUE ;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSObjectCreator::TestCPVRun2(Double_t dx, Double_t dz, Double_t pt, Int_t charge)
{
  //Parameterization of LHC15o period
  //_true if neutral_

  Double_t meanX=0.;
  Double_t meanZ=0.;

  Double_t sx = TMath::Min(5.2, 1.160 + 0.52 * TMath::Exp(-0.042 * pt*pt) + 5.1/TMath::Power(pt+0.62,3));
  Double_t sz = TMath::Min(3.3, 1.10  + 0.39 * TMath::Exp(-0.027 * pt*pt) + 0.70 /TMath::Power(pt+0.223,3));

  Double_t mf1 = fEvent->GetMagneticField(); //Positive for ++ and negative for --

  if(mf1<0.){ //field --
    meanZ = 0.077;
    if(charge>0)
      meanX =  TMath::Min(5.8, 0.2 + 0.7 * TMath::Exp(-0.019 * pt*pt) + 34./TMath::Power(pt+1.39,3));
    else
      meanX = -TMath::Min(5.8, 0.1 + 0.7 * TMath::Exp(-0.014 * pt*pt) + 30./TMath::Power(pt+1.36,3));
  }
  else{ //Field ++
    meanZ= 0.077;
    if(charge>0)
      meanX = -TMath::Min(5.8, 0.3 + 0.7 * TMath::Exp(-0.012 * pt*pt) + 35./TMath::Power(pt+1.43,3));
    else
      meanX =  TMath::Min(5.8, 0.2 + 0.6 * TMath::Exp(-0.014 * pt*pt) + 28./TMath::Power(pt+1.27,3));
  }

  Double_t rz=(dz-meanZ)/sz ;
  Double_t rx=(dx-meanX)/sx ;
  return TMath::Sqrt(rx*rx+rz*rz) ;

}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSObjectCreator::TestLambda(Double_t e,Double_t l1,Double_t l2)
{
  Double_t l2Mean = 1.53126+9.50835e+06/(1.+1.08728e+07*e+1.73420e+06*e*e) ;
  Double_t l1Mean = 1.12365+0.123770*TMath::Exp(-e*0.246551)+5.30000e-03*e ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*e+5.01265e+05*e*e)+9.000e-03*e;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*e+6.78604e-07*e*e)+9.000e-03*e;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*e) ;
  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma +
    0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
    0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;

  return R2;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSObjectCreator::TestCoreLambda(Double_t e, Double_t l1, Double_t l2)
{
  //Evaluates if lambdas correspond to photon cluster
  //Tuned using pp date
  //For core radius R=4.5

  Double_t   l1Mean  = 1.150200 + 0.097886/(1.+1.486645*e+0.000038*e*e) ;
  Double_t   l2Mean = 1.574706 + 0.997966*exp(-0.895075*e)-0.010666*e ;
  Double_t   l1Sigma = 0.100255 + 0.337177*exp(-0.517684*e)+0.001170*e ;
  Double_t   l2Sigma = 0.232580 + 0.573401*exp(-0.735903*e)-0.002325*e ;
  Double_t   c = -0.110983 -0.017353/(1.-1.836995*e+0.934517*e*e) ;

  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma +
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return R2 ;
}
//________________________________________________________________________
void AliAnalysisTaskPHOSObjectCreator::EvalCoreLambdas(AliVCluster *clu, AliVCaloCells *cells,Double_t &m02, Double_t &m20)
{
  //calculate dispecrsion of the cluster in the circle with radius distanceCut around the maximum

  const Double_t rCut=4.5;

  const Double32_t *elist = clu->GetCellsAmplitudeFraction();
  //Calculates the center of gravity in the local PHOS-module coordinates
  Float_t wtot = 0;
  const Int_t mulDigit=clu->GetNCells() ;
  Double_t xc[mulDigit] ;
  Double_t zc[mulDigit] ;
  Double_t wi[mulDigit] ;
  Double_t x = 0 ;
  Double_t z = 0 ;
  const Double_t logWeight=4.5 ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t relid[4] ;
    Float_t xi=0. ;
    Float_t zi=0. ;
    Int_t absId = clu->GetCellAbsId(iDigit) ;
    fPHOSGeo->AbsToRelNumbering(absId, relid) ;
    fPHOSGeo->RelPosInModule(relid, xi, zi);
    xc[iDigit]=xi ;
    zc[iDigit]=zi ;
    Double_t ei = elist[iDigit]*cells->GetCellAmplitude(absId) ;
    wi[iDigit]=0. ;
    if (clu->E()>0 && ei>0) {
      wi[iDigit] = TMath::Max( 0., logWeight + TMath::Log( ei / clu->E() ) ) ;
      Double_t w=wi[iDigit];
      x    += xc[iDigit] * w ;
      z    += zc[iDigit] * w ;
      wtot += w ;
    }
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }

  wtot = 0. ;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;
  Double_t xCut = 0.;
  Double_t zCut = 0.;

  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Double_t w=wi[iDigit];
    if(w>0.) {
      Double_t xi= xc[iDigit] ;
      Double_t zi= zc[iDigit] ;

      if((xi-x)*(xi-x)+(zi-z)*(zi-z) < rCut*rCut){
        xCut += w * xi ;
        zCut += w * zi ;
        dxx  += w * xi * xi ;
        dzz  += w * zi * zi ;
        dxz  += w * xi * zi ;
        wtot += w ;
      }
    }

  }
  if(wtot>0) {
    xCut/= wtot ;
    zCut/= wtot ;
    dxx /= wtot ;
    dzz /= wtot ;
    dxz /= wtot ;
    dxx -= xCut * xCut ;
    dzz -= zCut * zCut ;
    dxz -= xCut * zCut ;

    m02 =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    m20 =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
  }
  else {
    m20=m02=0.;
  }

}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSObjectCreator::CoreEnergy(AliVCluster *clu, AliVCaloCells *cells)
{
  //calculate energy of the cluster in the circle with radius distanceCut around the maximum

  //Can not use already calculated coordinates?
  //They have incidence correction...
  const Double_t distanceCut =3.5;//default value
  const Double_t logWeight=4.5;

  const Double32_t * elist = clu->GetCellsAmplitudeFraction() ;
// Calculates the center of gravity in the local PHOS-module coordinates
  Float_t wtot = 0;
  const Int_t mulDigit=clu->GetNCells() ;
  Double_t xc[mulDigit] ;
  Double_t zc[mulDigit] ;
  Double_t ei[mulDigit] ;
  Double_t x = 0 ;
  Double_t z = 0 ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iDigit), relid) ;
    fPHOSGeo->RelPosInModule(relid, xi, zi);
    xc[iDigit]=xi ;
    zc[iDigit]=zi ;
    ei[iDigit]=elist[iDigit]*cells->GetCellAmplitude(clu->GetCellsAbsId()[iDigit]);

    if( fDebug >= 3 )
      printf("%f ",ei[iDigit]);
    if (clu->E()>0 && ei[iDigit]>0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( ei[iDigit] / clu->E() ) ) ;
      x    += xc[iDigit] * w ;
      z    += zc[iDigit] * w ;
      wtot += w ;
    }
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }
  Double_t coreE=0. ;
  for(Int_t iDigit=0; iDigit < mulDigit; iDigit++) {
    Double_t distance = TMath::Sqrt((xc[iDigit]-x)*(xc[iDigit]-x)+(zc[iDigit]-z)*(zc[iDigit]-z)) ;
    if(distance < distanceCut)
      coreE += ei[iDigit] ;
  }

  //Apply non-linearity correction
  return fNonLinCorr->Eval(coreE);
}
//________________________________________________________________________
Double_t AliAnalysisTaskPHOSObjectCreator::CoreEnergy(AliVCluster *clu)
{
 //calculate energy of the cluster in the circle with radius distanceCut around the maximum

  //Can not use already calculated coordinates?
  //They have incidence correction...
  const Double_t distanceCut =3.5 ;
  const Double_t logWeight=4.5 ;

  Double32_t * elist = clu->GetCellsAmplitudeFraction() ;
// Calculates the center of gravity in the local PHOS-module coordinates
  Float_t wtot = 0;
  Double_t xc[100]={0} ;
  Double_t zc[100]={0} ;
  Double_t x = 0 ;
  Double_t z = 0 ;
  Int_t mulDigit=TMath::Min(100,clu->GetNCells()) ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iDigit), relid) ;
    fPHOSGeo->RelPosInModule(relid, xi, zi);
    xc[iDigit]=xi ;
    zc[iDigit]=zi ;
    if (clu->E()>0 && elist[iDigit]>0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( elist[iDigit] / clu->E() ) ) ;
      x    += xc[iDigit] * w ;
      z    += zc[iDigit] * w ;
      wtot += w ;
    }
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }
  Double_t coreE=0. ;
  for(Int_t iDigit=0; iDigit < mulDigit; iDigit++) {
    Double_t distance = TMath::Sqrt((xc[iDigit]-x)*(xc[iDigit]-x)+(zc[iDigit]-z)*(zc[iDigit]-z)) ;
    if(distance < distanceCut)
      coreE += elist[iDigit] ;
  }
  //Apply non-linearity correction
  return coreE ;

}
//________________________________________________________________________
void AliAnalysisTaskPHOSObjectCreator::EvalLambdas(AliVCluster * clu, Double_t &m02, Double_t &m20)
{
  //calculate dispecrsion of the cluster in the circle with radius distanceCut around the maximum

  const Double_t rCut=4.5 ;

  Double32_t * elist = clu->GetCellsAmplitudeFraction() ;
// Calculates the center of gravity in the local PHOS-module coordinates
  Float_t wtot = 0;
  Double_t xc[100]={0} ;
  Double_t zc[100]={0} ;
  Double_t x = 0 ;
  Double_t z = 0 ;
  Int_t mulDigit=TMath::Min(100,clu->GetNCells()) ;
  const Double_t logWeight=4.5 ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iDigit), relid) ;
    fPHOSGeo->RelPosInModule(relid, xi, zi);
    xc[iDigit]=xi ;
    zc[iDigit]=zi ;
    if (clu->E()>0 && elist[iDigit]>0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( elist[iDigit] / clu->E() ) ) ;
      x    += xc[iDigit] * w ;
      z    += zc[iDigit] * w ;
      wtot += w ;
    }
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }

  wtot = 0. ;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;
  Double_t xCut = 0. ;
  Double_t zCut = 0. ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    if (clu->E()>0 && elist[iDigit]>0.) {
        Double_t w = TMath::Max( 0., logWeight + TMath::Log( elist[iDigit] / clu->E() ) ) ;
        Double_t xi= xc[iDigit] ;
        Double_t zi= zc[iDigit] ;
  if((xi-x)*(xi-x)+(zi-z)*(zi-z) < rCut*rCut){
          xCut += w * xi ;
          zCut += w * zi ;
          dxx  += w * xi * xi ;
          dzz  += w * zi * zi ;
          dxz  += w * xi * zi ;
          wtot += w ;
  }
    }

  }

 if (wtot>0) {
    xCut/= wtot ;
    zCut/= wtot ;
    dxx /= wtot ;
    dzz /= wtot ;
    dxz /= wtot ;
    dxx -= xCut * xCut ;
    dzz -= zCut * zCut ;
    dxz -= xCut * zCut ;

    m02 =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    m20 =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
  }
  else {
    m20=m02=0.;
  }

}
//________________________________________________________________________
Int_t AliAnalysisTaskPHOSObjectCreator::FindTrackMatching(Int_t mod,TVector3 *locpos, Double_t &dx, Double_t &dz, Double_t &pt,Int_t &charge)
{
  //Find track with closest extrapolation to cluster
  AliESDEvent *esd = fESDEvent;
  AliAODEvent *aod = fAODEvent;

  if(!esd && !aod){
    AliError("Neither AOD nor ESD was found") ;
    return -1;
  }
  Double_t  magF =0.;
   if(esd)
     magF = esd->GetMagneticField();
   if(aod)
     magF = aod->GetMagneticField();

  Double_t magSign = 1.0;
  if(magF<0)magSign = -1.0;

  if (!TGeoGlobalMagField::Instance()->GetField()) {
    AliError("Margnetic filed was not initialized, use default") ;
    AliMagF* field = new AliMagF("Maps","Maps", magSign, magSign, AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }

  // *** Start the matching
  Int_t nt = 0;
  if(esd)
    nt = esd->GetNumberOfTracks();
  else
    nt = aod->GetNumberOfTracks();

  //Calculate actual distance to PHOS module
  TVector3 globaPos ;
  fPHOSGeo->Local2Global(mod, 0.,0., globaPos) ;
  const Double_t rPHOS = globaPos.Pt() ; //Distance to center of  PHOS module
  const Double_t kYmax = 72.+10. ; //Size of the module (with some reserve) in phi direction
  const Double_t kZmax = 64.+10. ; //Size of the module (with some reserve) in z direction
  const Double_t kAlpha0=330./180.*TMath::Pi() ; //First PHOS module angular direction
  const Double_t kAlpha= 20./180.*TMath::Pi() ; //PHOS module angular size
  Double_t minDistance = 1.e6;

  Double_t gposTrack[3] ;

  Double_t bz = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->SolenoidField();
  bz = TMath::Sign(0.5*kAlmost0Field,bz) + bz;

  Double_t b[3];
  Int_t itr=-1 ;
  AliESDtrack *esdTrack=0x0 ;
  AliAODTrack *aodTrack=0x0 ;
  Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
  for (Int_t i=0; i<nt; i++) {
      if(esd)
        esdTrack=esd->GetTrack(i);
      else
        aodTrack=(AliAODTrack*)aod->GetTrack(i);

      // Skip the tracks having "wrong" status (has to be checked/tuned)
      if(esd){
        ULong_t status = esdTrack->GetStatus();
        if((status & AliESDtrack::kTPCout) == 0) continue;
      }
      else{


      }


      //Continue extrapolation from TPC outer surface
      AliExternalTrackParam outerParam;
      if(esdTrack){
          outerParam = *(esdTrack->GetOuterParam());
      }
      if(aodTrack){
        aodTrack->GetPxPyPz(pxpypz);
        aodTrack->GetXYZ(xyz);
        aodTrack->GetCovarianceXYZPxPyPz(cv);
        outerParam.Set(xyz,pxpypz,cv,aodTrack->Charge());
      }

      Double_t z;
      if(!outerParam.GetZAt(rPHOS,bz,z))
        continue ;

      if (TMath::Abs(z) > kZmax)
        continue; // Some tracks miss the PHOS in Z

  //Direction to the current PHOS module
      Double_t phiMod=kAlpha0-kAlpha*mod ;
      if(!outerParam.RotateParamOnly(phiMod)) continue ; //RS use faster rotation if errors are not needed

      Double_t y;                       // Some tracks do not reach the PHOS
      if (!outerParam.GetYAt(rPHOS,bz,y)) continue; //    because of the bending

      if(TMath::Abs(y) < kYmax){
        outerParam.GetBxByBz(b) ;
        outerParam.PropagateToBxByBz(rPHOS,b);        // Propagate to the matching module
        //outerParam.CorrectForMaterial(...); // Correct for the TOF material, if needed
        outerParam.GetXYZ(gposTrack) ;
        TVector3 globalPositionTr(gposTrack) ;
        TVector3 localPositionTr ;
        fPHOSGeo->Global2Local(localPositionTr,globalPositionTr,mod) ;
        Double_t ddx = locpos->X()-localPositionTr.X();
        Double_t ddz = locpos->Z()-localPositionTr.Z();
        Double_t d2 = ddx*ddx + ddz*ddz;
        if(d2 < minDistance) {
    dx = ddx ;
      dz = ddz ;
    minDistance=d2 ;
    itr=i ;
          if(esdTrack){
            pt=esdTrack->Pt() ;
      charge=esdTrack->Charge() ;
          }
          else{
           pt=aodTrack->Pt() ;
           charge=aodTrack->Charge() ;
          }
        }
      }
    }//Scanned all tracks

   return itr ;

}
//________________________________________________________________________
void AliAnalysisTaskPHOSObjectCreator::DistanceToBadChannel(Int_t mod, TVector3 * locPos, Double_t &minDist)
{
 //Check if distance to bad channel was reduced
  Int_t range = minDist/2.2 +1 ; //Distance at which bad channels should be serached

  Int_t relid[4]={0,0,0,0} ;
  fPHOSGeo->RelPosToRelId(mod, locPos->X(), locPos->Z(), relid) ;
  Int_t xmin=TMath::Max(1,relid[2]-range) ;
  Int_t xmax=TMath::Min(64,relid[2]+range) ;
  Int_t zmin=TMath::Max(1,relid[3]-range) ;
  Int_t zmax=TMath::Min(56,relid[3]+range) ;

  Float_t x=0.,z=0.;
  for(Int_t ix=xmin;ix<=xmax;ix++){
    for(Int_t iz=zmin;iz<=zmax;iz++){
      if(fPHOSBadMap[mod] && fPHOSBadMap[mod]->GetBinContent(ix,iz)>0){ //Bad channel
        Int_t relidBC[4]={mod,0,ix,iz} ;
        fPHOSGeo->RelPosInModule(relidBC,x,z);
        Double_t dist = TMath::Sqrt((x-locPos->X())*(x-locPos->X()) + (z-locPos->Z())*(z-locPos->Z()));
        if(dist<minDist) minDist = dist;
      }
    }
  }

}
//________________________________________________________________________
Int_t AliAnalysisTaskPHOSObjectCreator::FindPrimary(AliCaloPhoton *ph,  Bool_t&sure)
{
  //Finds primary and estimates if it unique one?
  //First check can it be photon/electron
  AliVCluster *clu = (AliVCluster*)ph->GetCluster();
  const Double_t emFraction=0.9; //part of energy of cluster to be assigned to EM particle
  Int_t n = clu->GetNLabels();

  if(fESDEvent){
    return clu->GetLabel();
  }
  else if(fAODEvent){
    for(Int_t i=0;  i<n;  i++){
      Int_t label = clu->GetLabelAt(i);
      AliAODMCParticle *p =  (AliAODMCParticle*)fMCArrayAOD->At(label);
      Int_t pdg = p->PdgCode() ;
      if(pdg==22  ||  pdg==11 || pdg == -11){
        if(p->E()>emFraction*clu->E()){
          sure=kTRUE ;
          return label;
        }
      }
    }

    Double_t *Ekin = new Double_t[n];

    for(Int_t i=0;  i<n;  i++){
      Int_t label = clu->GetLabelAt(i);
      AliAODMCParticle* p = (AliAODMCParticle*)fMCArrayAOD->At(label);
      Ekin[i]=p->P() ;  // estimate of kinetic energy
      if(p->PdgCode()==-2212  ||  p->PdgCode()==-2112){
        Ekin[i]+=1.8  ;  //due to annihilation
      }
    }

    Int_t iMax=0;
    Double_t eMax=0.,eSubMax=0. ;
    for(Int_t i=0;  i<n;  i++){
      if(Ekin[i]>eMax){
        eSubMax=eMax;
        eMax=Ekin[i];
        iMax=i;
      }
    }
    if(eSubMax>0.8*eMax)//not obvious primary
      sure=kFALSE;
    else
      sure=kTRUE;

    delete[]  Ekin;

    return clu->GetLabelAt(iMax);

  }
  else{
    return clu->GetLabel();
  }
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPHOSObjectCreator::PassSTDCut(AliVCluster *cluster)
{
  if(cluster->GetM20() > 2.0) return kFALSE;
  if(cluster->E() > 1.0 && cluster->GetM02() < 0.1) return kFALSE;
  if(cluster->E() > 2.0 && cluster->GetM20() < 0.1) return kFALSE;
  return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskPHOSObjectCreator::EstimateSTDCutEfficiency(TClonesArray *array)
{
  const Int_t multClust = array->GetEntriesFast();

  TLorentzVector p12;
  Double_t m12=0;
  Double_t energy=0;
  Double_t weight = 1.;

  if(fIsSingleSim || fIsEmbedding){
    TF1 *f1 = new TF1("f1","[0]/TMath::TwoPi() * ([2]-1)*([2]-2)/([2]*[1]*([2]*[1] + 0.139*([2]-2) )) * TMath::Power(1+(TMath::Sqrt(x*x+0.139*0.139) - 0.139)/([2]*[1]),-[2])",0,100);//1/2pi x 1/Nev x 1/pT x d2N/dpTdy //TSallis pi0 in pp
    f1->SetNpx(1000);
    f1->SetParameters(2.70,0.132,6.64);
    AliAODMCParticle *p = (AliAODMCParticle*)fMCArrayAOD->At(0);//0 is always generated particle in single simulation.
    Double_t pT = p->Pt();
    weight = pT * f1->Eval(pT);

    delete f1;
    f1 = 0x0;
  }

  AliInfo(Form("weight is %e",weight));

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)array->At(i1);
    if(ph1->GetNsigmaCoreDisp() > 2.5) continue;
    if(ph1->Energy() < 1.0) continue;//to get high S/B

    for(Int_t i2=0;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)array->At(i2);
      AliVCluster *cluster = (AliVCluster*)ph2->GetCluster();

      if(i2==i1) continue;//reject same cluster combination

      p12 = *ph1 + *ph2;
      m12 = p12.M();
      energy = ph2->Energy();

      fHistoMggvsEProbe->Fill(m12,energy,weight);
      if(PassSTDCut(cluster)) fHistoMggvsEPassingProbe->Fill(m12,energy,weight);

    }//end of ph2

  }//end of ph1

}



