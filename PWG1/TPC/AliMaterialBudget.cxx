



//
//

// This class estiamtes the material budget of the inner detectors in ALICE based
// on the "upper/lower track matching"-method of the ALICE TPC.

//
//  
//

//




// ROOT includes
#include <TChain.h>
#include <TMath.h>
#include <TVectorD.h>
#include <TSystem.h>
#include <TFile.h>
#include "TGeoGlobalMagField.h"

// ALIROOT includes
#include <TTreeStream.h>
#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

#include <AliESD.h>
#include "AliMaterialBudget.h"
#include "AliGenInfoMaker.h"
#include "AliHelix.h"

//
#include "AliMCInfo.h"
#include "AliESDRecInfo.h"
#include "AliTPCParamSR.h"
#include "AliTracker.h"
#include "AliTPCseed.h"

// STL includes
#include <iostream>

using namespace std;

ClassImp(AliMaterialBudget)

//________________________________________________________________________
AliMaterialBudget::AliMaterialBudget() : 
  AliAnalysisTask(), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(0),
  fDebugOutputPath(),
  fListHist(0),
  fHistMult(0),
  fCutMaxD(5),        // maximal distance in rfi ditection
  fCutMaxDz(40),      // maximal distance in z ditection
  fCutTheta(0.03),    // maximal distan theta
  fCutMinDir(-0.99)   // direction vector products
{
  //
  // Default constructor (should not be used)
  //
}

AliMaterialBudget::AliMaterialBudget(const AliMaterialBudget& /*info*/) : 
  AliAnalysisTask(), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  //
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(),
  fDebugOutputPath(), 
  fListHist(0),
  fHistMult(0),
  fCutMaxD(5),        // maximal distance in rfi ditection
  fCutMaxDz(40),      // maximal distance in z ditection
  fCutTheta(0.03),    // maximal distan theta
  fCutMinDir(-0.99)   // direction vector products
{
  //
  // Default constructor 
  //
}



//________________________________________________________________________
AliMaterialBudget::AliMaterialBudget(const char *name) : 
  AliAnalysisTask(name, "AliMaterialBudget"), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(0),
  fDebugOutputPath(),
  fListHist(0),
  fHistMult(0),
  fCutMaxD(5),        // maximal distance in rfi ditection
  fCutMaxDz(40),      // maximal distance in z ditection
  fCutTheta(0.03),    // maximal distan theta
  fCutMinDir(-0.99)   // direction vector products
{
  //
  // Normal constructor
  //
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList
  DefineOutput(0, TList::Class());
  //
  //
}

AliMaterialBudget::~AliMaterialBudget(){
  //
  //
  //
  if (fDebugLevel>0)  printf("AliMaterialBudget::~AliMaterialBudget\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer=0;
}


//________________________________________________________________________
void AliMaterialBudget::ConnectInputData(Option_t *) 
{
  //
  // Connect the input data
  //
  if(fDebugLevel>3)
    cout << "AnalysisTaskTPCCluster::ConnectInputData()" << endl;

  TTree* tree=dynamic_cast<TTree*>(GetInputData(0));
  if (!tree) {
    //Printf("ERROR: Could not read chain from input slot 0");
  }
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdH) {
      //Printf("ERROR: Could not get ESDInputHandler");
    }
    else {
      esdH->SetActiveBranches("ESDfriend");
      fESD = esdH->GetEvent();
      //Printf("*** CONNECTED NEW EVENT ****");
    }  
  }
  AliMCEventHandler* mcinfo = (AliMCEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());  
  mcinfo->SetReadTR(kTRUE);
  
  fMCinfo = mcinfo->MCEvent();


}






//________________________________________________________________________
void AliMaterialBudget::CreateOutputObjects() 
{
  //
  // Connect the output objects
  //
  if(fDebugLevel>3)
    cout << "AnalysisTaskTPCCluster::CreateOutputObjects()" << endl;
  //
  fListHist = new TList();
  fListHist->SetOwner();
  //
  fHistMult = new TH1F("HistMult", "Number of Tracks per Event; number of tracks per event; number of tracks",501,-0.5,500.5);
  fListHist->Add(fHistMult);


}


//________________________________________________________________________
void AliMaterialBudget::Exec(Option_t *) {
  //
  // Execute analysis for current event 
  //

  if(fDebugLevel>3)
    cout << "AliMaterialBudget::Exec()" << endl;

  fHistMult->Fill(fESD->GetNumberOfTracks());
  //FindPairs(fESD); // nearly everything takes place in find pairs...

  // If MC has been connected   

  if (!fMCinfo){
    cout << "Not MC info\n" << endl;
  }else{
    ProcessMCInfo();
    //mcinfo->Print();
    //DumpInfo();
  }
  //
  PostData(0, fListHist);
}      




//________________________________________________________________________
void AliMaterialBudget::Terminate(Option_t *) {
    //
    // Terminate loop
    //
  if(fDebugLevel>3)
    printf("AliMaterialBudget: Terminate() \n");  
  //
  if (fDebugLevel>0) printf("AliMCtrackingTestTask::Terminate\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer = 0;
  return;
}



TTreeSRedirector *AliMaterialBudget::GetDebugStreamer(){
  //
  // Get Debug streamer
  // In case debug streamer not yet initialized and StreamLevel>0 create new one
  //
  if (fStreamLevel==0) return 0;
  if (fDebugStreamer) return fDebugStreamer;
  TString dsName;
  dsName=GetName();
  dsName+="Debug.root";
  dsName.ReplaceAll(" ","");
  fDebugStreamer = new TTreeSRedirector(dsName.Data());
  return fDebugStreamer;
}






AliExternalTrackParam * AliMaterialBudget::MakeTrack(const AliTrackReference* ref, TParticle*part)
{
  //
  // Make track out of the track ref
  // part - TParticle used to determine chargr
  // the covariance matrix - equal 0 - starting from ideal MC position
  if (!part->GetPDG()) return 0;
  Double_t xyz[3]={ref->X(),ref->Y(),ref->Z()};
  Double_t pxyz[3]={ref->Px(),ref->Py(),ref->Pz()};
  Int_t charge = TMath::Nint(part->GetPDG()->Charge()/3.);
  if (ref->X()*ref->Px()+ref->Y()*ref->Py() <0){
    pxyz[0]*=-1;
    pxyz[1]*=-1;
    pxyz[2]*=-1;
    charge*=-1;
  }
  Double_t cv[21];
  for (Int_t i=0; i<21;i++) cv[i]=0;
  AliExternalTrackParam * param = new AliExternalTrackParam(xyz,pxyz,cv,charge);
  return param;
}

Bool_t  AliMaterialBudget::PropagateToPoint(AliExternalTrackParam *param, Double_t *xyz, Double_t mass, Float_t step){
  // 
  // Propagate track to point xyz using 
  // AliTracker::PropagateTo functionality
  //
  //  param - track parameters
  //  xyz   - position to propagate
  //  mass  - particle mass
  //  step  - step to be used
  Double_t radius=TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  AliTracker::PropagateTrackToBxByBz(param, radius+step, mass, step, kTRUE,0.99);
  AliTracker::PropagateTrackToBxByBz(param, radius+0.5, mass, step*0.1, kTRUE,0.99);
  Double_t sxyz[3]={0,0,0};
  AliESDVertex vertex(xyz,sxyz);
  Bool_t isOK = param->PropagateToDCA(&vertex,AliTracker::GetBz(),10);
  return isOK;
}


void  AliMaterialBudget::ProcessMCInfo(){
  //
  //
  //
  //
  TParticle * particle= new TParticle;
  TClonesArray * trefs = new TClonesArray("AliTrackReference");
  const Double_t kPcut=0.05;
  const Double_t kMinDrITS = 2.;   // minimal distance between references
  const Double_t kMinDrTRD = 8.;   // minimal distance between references
  const Double_t kMinDrTOF = 10.;   // minimal distance between references
  //
  //
  // Process tracks
  //
  Int_t npart = fMCinfo->GetNumberOfTracks();
  if (npart==0) return;
  Double_t vertex[4]={0,0,0,0};
  fMCinfo->GetParticleAndTR(0, particle, trefs);
  if (particle){
    vertex[0]=particle->Vx();
    vertex[1]=particle->Vy();
    vertex[2]=particle->Vz();
    vertex[3]=particle->R();
  }
  //
  //
  AliTrackReference dummy,*pdummy= &dummy;
  AliExternalTrackParam epdummy,*pepdummy= &epdummy;
  Int_t nRefITS =0;
  Int_t nRefTPC =0;
  Int_t nRefTRD =0;
  Int_t nRefTOF =0;
  AliTrackReference * refITS0, *refITS1;
  AliTrackReference * refTPC0, *refTPC1;
  AliTrackReference * refTPCIn0, *refTPCIn1;
  AliTrackReference * refTRD0, *refTRD1;
  AliTrackReference * refTOF0, *refTOF1;
  AliTrackReference *refMinR;
  //
  for (Int_t ipart=0;ipart<npart;ipart++){
    //Int_t status = fMCinfo->GetParticleAndTR(ipart, particle, trefs);
    AliMCParticle * pp = (AliMCParticle*) fMCinfo->GetTrack(ipart);
    if (!pp) continue;
    if (particle->P()<kPcut) continue;
    Double_t mass = particle->GetMass();

    // RESET
    nRefITS =0;
    nRefTPC =0;
    nRefTRD =0;
    nRefTOF =0;
    refITS0=pdummy;    refITS1=pdummy;
    refTPC0=pdummy;    refTPC1=pdummy;
    refTPCIn0=pdummy;    refTPCIn1=pdummy;
    refTRD0=pdummy;    refTRD1=pdummy;
    refTOF0=pdummy;    refTOF1=pdummy;
    refMinR = pdummy;
    //
    Int_t nref = pp->GetNumberOfTrackReferences();
    if (nref==0) continue;
    for (Int_t iref = 0; iref < nref; iref++) { 
      AliTrackReference *ref = pp->GetTrackReference(iref);
      if (ref->DetectorId()==AliTrackReference::kDisappeared) continue;      
      //      if (ref.Px()*particle.Px()+ref.Py()*particle.Py()<0) break; // returning track
      //
      if (ref->DetectorId()==AliTrackReference::kITS){
	if (TMath::Abs(ref->R()-refITS1->R())>kMinDrITS) {
	  refITS1 = ref;
	  nRefITS++;
	}
	if (refITS0==pdummy) refITS0=ref;
      }
      if (ref->DetectorId()==AliTrackReference::kTPC){	
	nRefTPC++;
	refTPC1 = ref;
	if (refTPC0==pdummy) refTPC0=ref;
      }
      if (ref->DetectorId()==AliTrackReference::kTRD){
	if (TMath::Abs(ref->R()-refTRD1->R())>kMinDrTRD) {
	  refTRD1 = ref;
	  nRefTRD++;
	}
	if (refTRD0==pdummy) refTRD0=ref;
      }
      if (ref->DetectorId()==AliTrackReference::kTOF){
	if (TMath::Abs(ref->X()-refTOF1->X()) + TMath::Abs(ref->Y()-refTOF1->Y()) + TMath::Abs(ref->Z()-refTOF1->Z())>kMinDrTOF) {
	  refTOF1 = ref;
	  nRefTOF++;
	}
	if (refTOF0==pdummy) refTOF0=ref;
      }
      //
      // "find inner track ref"
      if (ref->DetectorId()==AliTrackReference::kTPC){	
	if (ref->Px()*ref->X()+ref->Py()*ref->Y()<0){
	  //  track in
	  if (refTPCIn0 == pdummy) refTPCIn0=ref;
	  if (refTPCIn0 != pdummy &&  refTPCIn0->R()>ref->R())
	    refTPCIn0=ref; 
	}
	if (ref->Px()*ref->X()+ref->Py()*ref->Y()>0){
	  //  track in
	  if (refTPCIn1 == pdummy) refTPCIn1=ref;
	  if (refTPCIn1 != pdummy &&  refTPCIn1->R()>ref->R())
	    refTPCIn1=ref; 
	}
      }


      if (refMinR==pdummy && ref->P()>0  ){
	refMinR=ref;
      }
      if (refMinR->R()>ref->R() && ref->P()>0 ) refMinR=ref;
      
    }
    //
    AliExternalTrackParam * trackMC = pepdummy;
    //track0->GetDZ(0,0,0,bz,dvertex0)
    Float_t dist[2]={0,0};
    AliMagF* field = (AliMagF*) TGeoGlobalMagField::Instance()->GetField();
    Double_t esdfield= fESD->GetMagneticField();
    Double_t xyz[3]={0,0,0};
    Double_t bxyz[3]={0,0,0};
    field->Field(xyz,bxyz);
    if (refMinR->P()>0) {
      trackMC = MakeTrack(refMinR,particle); 
      trackMC->GetDZ(0,0,0,bxyz[2],dist);
    }
    Double_t alphaTOF0 = TMath::ATan2(refTOF0->Y(),refTOF0->X());
    Double_t alphaTOF1 = TMath::ATan2(refTOF1->Y(),refTOF1->X());
    Int_t dsecTOF   = TMath::Nint(180*(alphaTOF0-alphaTOF1)/(TMath::Pi()*20.)-0.5);
    //
    // make the two different TPC tracks and propagate them to their DCA
    //
    Double_t dP = 0;
    Bool_t statusProp = kFALSE;
    Double_t dY = 0;
    Double_t dZ = 0;
    AliExternalTrackParam * track0  = pepdummy;
    AliExternalTrackParam * track1  = pepdummy;
    AliExternalTrackParam * otrack0 = pepdummy;
    AliExternalTrackParam * otrack1 = pepdummy;
    if (refTPCIn0!=pdummy && refTPCIn1!=pdummy) {
      track0 = MakeTrack(refTPCIn0,particle); 
      track1 = MakeTrack(refTPCIn1,particle);
      otrack0 = MakeTrack(refTPCIn0,particle); 
      otrack1 = MakeTrack(refTPCIn1,particle);
      dP = track0->P() - track1->P(); // momentum loss
      statusProp = AliMaterialBudget::PropagateCosmicToDCA(track0,track1,mass);
      if (statusProp) {
	dY = track0->GetY() - track1->GetY();
	dZ = track0->GetZ() - track1->GetZ();
      }
    }
    //
    TTreeSRedirector *pcstream = GetDebugStreamer();
    if (pcstream){      
      char name[100];
      for (Int_t id=0; id<3;id++){
	
	//	if (id==0) sprintf(name,"mcAll"); // all tracks: inconvenient to cut if on is only interest in tracks which reach the TPC
	if (id==0) continue; // require TPC
	if (id==1) sprintf(name,"mcITS");
	if (id==2) sprintf(name,"mcTPC");
	if (id==1&& nRefITS==0) continue;
	if (id==2&& nRefTPC==0) continue;

	(*pcstream)<<name<<
	  "ipart="<<ipart<<
	  "p.="<<particle<<
	  "mass="<<mass<<
	  "tbfield="<<bxyz[2]<<
	  "esdbfield="<<esdfield<<
	  // counter
	  "nref="<<nref<<
	  "nRefITS="<<nRefITS<<     
	  "nRefTPC="<<nRefTPC<<
	  "nRefTRD="<<nRefTRD<<
	  "nRefTOF="<<nRefTOF<<
	  //references
	  "refMinR.="<<refMinR<<
	  "refITS0.="<<refITS0<<
	  "refITS1.="<<refITS1<<
	  "refTPC0.="<<refTPC0<<
	  "refTPC1.="<<refTPC1<<
	  "refTPCIn0.="<<refTPCIn0<<
	  "refTPCIn1.="<<refTPCIn1<<
	  "refTRD0.="<<refTRD0<<
	  "refTRD1.="<<refTRD1<<
	  "refTOF0.="<<refTOF0<<
	  "refTOF1.="<<refTOF1<<
	  //trigger variables
	  "dsecTOF="<<dsecTOF<<    // delta TOF sectors
	  "aTOF0="<<alphaTOF0<<
	  "aTOF1="<<alphaTOF1<<
	  //
	  // track
	  "dr="<<dist[0]<<
	  "dz="<<dist[1]<<
	  //
	  // "two" TPC tracks 
	  "status="<<statusProp<<
	  "dP="<<dP<<
	  "otrack0.="<<otrack0<<
	  "otrack1.="<<otrack1<<
	  "track0.="<<track0<<
	  "track1.="<<track1<<
	  "dY="<<dY<<
	  "dZ="<<dZ<<
	  "\n";
      }
    }    
  }
}



void AliMaterialBudget::ProcessRefTracker(AliTrackReference* refIn,  AliTrackReference* refOut, TParticle*part,Int_t type){
  //
  // Test propagation from In to out
  //
  AliExternalTrackParam *param = 0;
  AliExternalTrackParam *paramMC = 0;
  Double_t xyzIn[3]={refIn->X(),refIn->Y(), refIn->Z()};
  Double_t mass = part->GetMass();
  Double_t step=1;
  //
  param=MakeTrack(refOut,part);
  paramMC=MakeTrack(refOut,part);
  if (!param) return;
  if (type<3) PropagateToPoint(param,xyzIn, mass, step);
  if (type==3) {
    AliTPCseed seed;
    seed.Set(param->GetX(),param->GetAlpha(),param->GetParameter(),param->GetCovariance());
    Float_t alpha= TMath::ATan2(refIn->Y(),refIn->X());
    if(seed.Rotate(alpha-seed.GetAlpha()) == kFALSE) return;
    seed.SetMass(mass);
    for (Float_t xlayer= seed.GetX(); xlayer>refIn->R(); xlayer-=step){
      seed.PropagateTo(xlayer);
    }
    seed.PropagateTo(refIn->R());
    param->Set(seed.GetX(),seed.GetAlpha(),seed.GetParameter(),seed.GetCovariance());
  }
  TTreeSRedirector *pcstream = GetDebugStreamer();
  TVectorD gpos(3);
  TVectorD gmom(3);
  param->GetXYZ(gpos.GetMatrixArray());
  param->GetPxPyPz(gmom.GetMatrixArray());
  if (pcstream){
    (*pcstream)<<"MC"<<
      "type="<<type<<
      "step="<<step<<
      "refIn.="<<refIn<<
      "refOut.="<<refOut<<
      "p.="<<part<<
      "par.="<<param<<   
      "parMC.="<<paramMC<<   
      "gpos.="<<&gpos<<
      "gmom.="<<&gmom<<
      "\n";
  }
}


void  AliMaterialBudget::FitTrackRefs(TParticle * part, TClonesArray * trefs){
  //
  //
  //
  //
  const Int_t kMinRefs=6;
  Int_t nrefs = trefs->GetEntries();
  if (nrefs<kMinRefs) return; // we should have enough references
  Int_t iref0 =-1;
  Int_t iref1 =-1;
  
  for (Int_t iref=0; iref<nrefs; iref++){
    AliTrackReference * ref = (AliTrackReference*)trefs->At(iref);
    if (!ref) continue;    
    Float_t dir = ref->X()*ref->Px()+ref->Y()*ref->Py();
    if (dir<0) break;
    if (ref->DetectorId()!=AliTrackReference::kTPC) continue;
    if (iref0<0) iref0 = iref;
    iref1 = iref;    
  }
  if (iref1-iref0<kMinRefs) return;
  Double_t covar[15];
  for (Int_t icov=0; icov<15; icov++) covar[icov]=0;
  covar[0]=1; 
  covar[2]=1; 
  covar[5]=1;
  covar[9]=1;
  covar[14]=1;

  AliTrackReference * refIn = (AliTrackReference*)trefs->At(iref0);
  AliTrackReference * refOut = (AliTrackReference*)trefs->At(iref1);
  AliExternalTrackParam *paramPropagate= MakeTrack(refIn,part);
  AliExternalTrackParam *paramUpdate   = MakeTrack(refIn,part);
  paramUpdate->AddCovariance(covar);
  Double_t mass = part->GetMass();
  Double_t charge = part->GetPDG()->Charge()/3.;
/*
  Float_t alphaIn= TMath::ATan2(refIn->Y(),refIn->X());
  Float_t radiusIn= refIn->R();
  Float_t alphaOut= TMath::ATan2(refOut->Y(),refOut->X());
  Float_t radiusOut= refOut->R();
*/
  Bool_t isOKP=kTRUE;
  Bool_t isOKU=kTRUE;
  AliMagF * field = (AliMagF*) TGeoGlobalMagField::Instance()->GetField();
  for (Int_t iref = iref0; iref<=iref1; iref++){
    AliTrackReference * ref = (AliTrackReference*)trefs->At(iref);
    Float_t alphaC= TMath::ATan2(ref->Y(),ref->X());
    Double_t pos[3] = {ref->X(), ref->Y(), ref->Z()};
    Double_t mag[3];
    field->Field(pos,mag);
    isOKP&=paramPropagate->Rotate(alphaC);
    isOKU&=paramUpdate->Rotate(alphaC);
    for (Float_t xref= paramPropagate->GetX(); xref<ref->R(); xref++){
      isOKP&=paramPropagate->PropagateTo(xref, mag[2]);
      isOKU&=paramUpdate->PropagateTo(xref, mag[2]);
    }
    isOKP&=paramPropagate->PropagateTo(ref->R(), mag[2]);
    isOKU&=paramUpdate->PropagateTo(ref->R(), mag[2]);
    Double_t clpos[2] = {0, ref->Z()};
    Double_t clcov[3] = { 0.005,0,0.005};
    isOKU&= paramUpdate->Update(clpos, clcov);  
  }
  TTreeSRedirector *pcstream = GetDebugStreamer();
  if (pcstream){
    TVectorD gposU(3);
    TVectorD gmomU(3);
    TVectorD gposP(3);
    TVectorD gmomP(3);
    paramUpdate->GetXYZ(gposU.GetMatrixArray());
    paramUpdate->GetPxPyPz(gmomU.GetMatrixArray());
    paramPropagate->GetXYZ(gposP.GetMatrixArray());
    paramPropagate->GetPxPyPz(gmomP.GetMatrixArray());

     (*pcstream)<<"MCupdate"<<
       "isOKU="<<isOKU<<
       "isOKP="<<isOKP<<
       "m="<<mass<<
       "q="<<charge<<
       "part.="<<part<<
       "refIn.="<<refIn<<
       "refOut.="<<refOut<<
       "pP.="<<paramPropagate<<
       "pU.="<<paramUpdate<<
       "gposU.="<<&gposU<<
       "gmomU.="<<&gmomU<<
       "gposP.="<<&gposP<<
       "gmomP.="<<&gmomP<<
       "\n";
   }
}



void AliMaterialBudget::FindPairs(AliESDEvent * event) {
  //
  // This function matches the cosmic tracks and calculates the energy loss in the material.
  // If accessible the "true" energy loss is determined with the MC track references.
  //
  //
  // Find cosmic pairs
  // 
  // Track0 is choosen in upper TPC part
  // Track1 is choosen in lower TPC part
  //
  if(fDebugLevel>3)
    cout << "AliMaterialBudget::FindPairs()" << endl;


  AliESDfriend *ESDfriend=static_cast<AliESDfriend*>(event->FindListObject("AliESDfriend"));
  Int_t ntracks=event->GetNumberOfTracks(); 
  TObjArray  tpcSeeds(ntracks);
  if (ntracks==0) return;
  Double_t vtxx[3]={0,0,0};
  Double_t svtxx[3]={0.000001,0.000001,100.};
  AliESDVertex vtx(vtxx,svtxx);
  //
  //track loop
  //
  for (Int_t i=0;i<ntracks;++i) {
   AliESDtrack *track = event->GetTrack(i);
     const AliExternalTrackParam * trackIn = track->GetInnerParam();
   const AliExternalTrackParam * trackOut = track->GetOuterParam();
   if (!trackIn) continue;
   if (!trackOut) continue;
   AliESDfriendTrack *friendTrack = ESDfriend->GetTrack(i);
   TObject *calibObject;
   AliTPCseed *seed = 0;
   for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
     if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
   }
   if (seed) tpcSeeds.AddAt(seed,i);
  }

  if (ntracks<2) return;
  //
  // Find pairs
  //
  for (Int_t i=0;i<ntracks;++i) {
    AliESDtrack *track0 = event->GetTrack(i);     
    // track0 - choosen upper part
    if (!track0) continue;
    if (!track0->GetOuterParam()) continue;
    if (track0->GetOuterParam()->GetAlpha()<0) continue;
    Double_t dir0[3];
    track0->GetDirection(dir0);    
    for (Int_t j=0;j<ntracks;++j) {
      if (i==j) continue;
      AliESDtrack *track1 = event->GetTrack(j);   
      //track 1 lower part
      if (!track1) continue;
      if (!track1->GetOuterParam()) continue;
      if (track1->GetOuterParam()->GetAlpha()>0) continue;
      //
      Double_t dir1[3];
      track1->GetDirection(dir1);
      
      AliTPCseed * seed0 = (AliTPCseed*) tpcSeeds.At(i);
      AliTPCseed * seed1 = (AliTPCseed*) tpcSeeds.At(j);
      if (! seed0) continue;
      if (! seed1) continue;
      //
      Float_t dir = (dir0[0]*dir1[0] + dir0[1]*dir1[1] + dir0[2]*dir1[2]);
      Float_t d0  = track0->GetLinearD(0,0);
      Float_t d1  = track1->GetLinearD(0,0);
      //
      // conservative cuts - convergence to be guarantied
      // applying before track propagation
      if (TMath::Abs(d0+d1)>fCutMaxD) continue;   // distance to the 0,0
      if (dir>fCutMinDir) continue;               // direction vector product
      Float_t bz = AliTracker::GetBz();
      Float_t dvertex0[2];   //distance to 0,0
      Float_t dvertex1[2];   //distance to 0,0 
      track0->GetDZ(0,0,0,bz,dvertex0);
      track1->GetDZ(0,0,0,bz,dvertex1);
      if (TMath::Abs(dvertex0[1])>250) continue;
      if (TMath::Abs(dvertex1[1])>250) continue;
      //
      //
      //
      Float_t dmax = TMath::Max(TMath::Abs(d0),TMath::Abs(d1));
      AliExternalTrackParam param0(*track0);
      AliExternalTrackParam param1(*track1);
      //
      // Propagate using Magnetic field and correct fo material budget
      //
      AliTracker::PropagateTrackToBxByBz(&param0,dmax+1,0.0005,3,kTRUE);
      AliTracker::PropagateTrackToBxByBz(&param1,dmax+1,0.0005,3,kTRUE);
      //
      // Propagate rest to the 0,0 DCA - z should be ignored
      //
      Bool_t b0 = param0.PropagateToDCA(&vtx,bz,1000);
      Bool_t b1 = param1.PropagateToDCA(&vtx,bz,1000);
      //      
      param0.GetDZ(0,0,0,bz,dvertex0);
      param1.GetDZ(0,0,0,bz,dvertex1);
      if (TMath::Abs(param0.GetZ()-param1.GetZ())>fCutMaxDz) continue;
      //
      Double_t xyz0[3];//,pxyz0[3];
      Double_t xyz1[3];//,pxyz1[3];
      param0.GetXYZ(xyz0);
      param1.GetXYZ(xyz1);
      Bool_t isPair = IsPair(&param0,&param1);
      //
      // HERE WE WILL PUT THE ACCESS TO THE MC TRACKS AND MATCH THESE !!!!
      //
      Int_t label0 = TMath::Abs(track0->GetLabel());
      AliMCParticle *mcParticle0 = (AliMCParticle*) fMCinfo->GetTrack(label0);
      TParticle *particle0 = mcParticle0->Particle();
      AliTrackReference *ref0 = GetFirstTPCTrackRef(mcParticle0); // get the first TPC track reference
      if (!ref0) continue;
      AliExternalTrackParam *paramMC0 = 0;
      paramMC0 = MakeTrack(ref0, particle0);
      //
      Int_t label1 = TMath::Abs(track1->GetLabel());
      AliMCParticle *mcParticle1 = (AliMCParticle*) fMCinfo->GetTrack(label1);
      TParticle *particle1 = mcParticle1->Particle();
      AliTrackReference *ref1 = GetFirstTPCTrackRef(mcParticle1); // get the first TPC track reference
      if (!ref1) continue;
      AliExternalTrackParam *paramMC1 = 0;
      paramMC1 = MakeTrack(ref1, particle1);
      //
      // ACCESS TOF INFORMATION
      Int_t nTrackRefTOF0 = 0;
      Int_t nTrackRefITS0 = 0;
      AliTrackReference * refLastTOF0 = 0;
      AliTrackReference * refFirstTOF0 = GetAllTOFinfo(mcParticle0, nTrackRefTOF0, nTrackRefITS0);
      Float_t alphaTOF0 = 0;
      if (refFirstTOF0) alphaTOF0 = refFirstTOF0->Alpha();
      //
      Int_t nTrackRefTOF1 = 0;
      Int_t nTrackRefITS1 = 0;
      AliTrackReference  * refLastTOF1 = 0;
      AliTrackReference  * refFirstTOF1 =GetAllTOFinfo(mcParticle1, nTrackRefTOF1, nTrackRefITS1);
      Float_t alphaTOF1 = 0;
      if (refFirstTOF1) alphaTOF1 = refFirstTOF1->Alpha();
      //cout << " STATUS "<<nTrackRefTOF0<<" "<<refFirstTOF0<<" "<<refLastTOF0<<" " <<refFirstTOF1<<" "<<refLastTOF1<<endl;
      //
      //
      //
      if (fStreamLevel>0){
	TTreeSRedirector * cstream =  GetDebugStreamer();
	AliExternalTrackParam *ip0 = (AliExternalTrackParam *)track0->GetInnerParam();
	AliExternalTrackParam *ip1 = (AliExternalTrackParam *)track1->GetInnerParam();
	AliExternalTrackParam *op0 = (AliExternalTrackParam *)track0->GetOuterParam();
	AliExternalTrackParam *op1 = (AliExternalTrackParam *)track1->GetOuterParam();
	//
	//
	//
	if (cstream) {
	  (*cstream) << "Track0" <<
	    "dir="<<dir<<               //  direction
	    "OK="<<isPair<<             //  will be accepted
	    "b0="<<b0<<                 //  propagate status
	    "b1="<<b1<<                 //  propagate status
	    //
	    "Particle.="<<particle0<<             // TParticle with generated momentum
	    "NTrackRefTOF0="<<nTrackRefTOF0<<     // Number of TOF track references upper
	    "NTrackRefTOF1="<<nTrackRefTOF1<<     // Number of TOF track references lower
	    "NTrackRefITS0="<<nTrackRefITS0<<     // Number of ITS track references upper
	    "NTrackRefITS1="<<nTrackRefITS1<<     // Number of ITS track references lower
      	    "Alpha0="<<alphaTOF0<<                   // alpha upper
	    "Alpha1="<<alphaTOF1<<                   // alpha lower
	    "RefFirstTOF0.="<<refFirstTOF0<<        // first tof reference upper
	    "RefLastTOF0.="<<refLastTOF0<<          // last tof reference upper
	    "RefFirstTOF1.="<<refFirstTOF1<<        // first tof reference lower
	    "RefLastTOF1.="<<refLastTOF1<<          // last tof reference lower
	    //
	    "Orig0.=" << track0 <<      //  original track  0
	    "Orig1.=" << track1 <<      //  original track  1
	    "Tr0.="<<&param0<<          //  track propagated to the DCA 0,0
	    "Tr1.="<<&param1<<          //  track propagated to the DCA 0,0	   
	    "Ip0.="<<ip0<<              //  inner param - upper
	    "Ip1.="<<ip1<<              //  inner param - lower
	    "Op0.="<<op0<<              //  outer param - upper
	    "Op1.="<<op1<<              //  outer param - lower
	    //
	    "paramTrackRef0.="<<paramMC0<< //  "ideal" MC true track parameters from track references - upper
	    "paramTrackRef1.="<<paramMC1<< //  "ideal" MC true track parameters from track references - lower
	    //
	    "v00="<<dvertex0[0]<<       //  distance using kalman
	    "v01="<<dvertex0[1]<<       // 
	    "v10="<<dvertex1[0]<<       //
	    "v11="<<dvertex1[1]<<       // 
	    "d0="<<d0<<                 //  linear distance to 0,0
	    "d1="<<d1<<                 //  linear distance to 0,0
	    //
	    //
	    //
	    "x00="<<xyz0[0]<<           // global position close to vertex
	    "x01="<<xyz0[1]<<
	    "x02="<<xyz0[2]<<
	    //
	    "x10="<<xyz1[0]<<           // global position close to vertex
	    "x11="<<xyz1[1]<<
	    "x12="<<xyz1[2]<<
	    //
	    "dir00="<<dir0[0]<<           // direction upper
	    "dir01="<<dir0[1]<<
	    "dir02="<<dir0[2]<<
	    //
	    "dir10="<<dir1[0]<<           // direction lower
	    "dir11="<<dir1[1]<<
	    "dir12="<<dir1[2]<<
	    //
	    //
	    "Seed0.=" << seed0 <<       //  original seed 0
	    "Seed1.=" << seed1 <<       //  original seed 1
	    //
	    "\n";
	}
      }
      delete paramMC0;
      delete paramMC1;
    }
  } 

  return;


}

Bool_t  AliMaterialBudget::IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1){
  //
  //
  /*
  // 0. Same direction - OPOSITE  - cutDir +cutT    
  TCut cutDir("cutDir","dir<-0.99")
  // 1. 
  TCut cutT("cutT","abs(Tr1.fP[3]+Tr0.fP[3])<0.03")
  //
  // 2. The same rphi 
  TCut cutD("cutD","abs(Tr0.fP[0]+Tr1.fP[0])<5")
  //
  //
  //
  TCut cutPt("cutPt","abs(Tr1.fP[4]+Tr0.fP[4])<1&&abs(Tr0.fP[4])+abs(Tr1.fP[4])<10");  
  // 1/Pt diff cut
  */
  const Double_t *p0 = tr0->GetParameter();
  const Double_t *p1 = tr1->GetParameter();
  if (TMath::Abs(p0[3]+p1[3])>fCutTheta) return kFALSE;
  if (TMath::Abs(p0[1]-p1[1])>fCutMaxDz) return kFALSE;
  if (TMath::Abs(p0[0]+p1[0])>fCutMaxD)  return kFALSE;
  
  Double_t d0[3], d1[3];
  tr0->GetDirection(d0);    
  tr1->GetDirection(d1);       
  if (d0[0]*d1[0] + d0[1]*d1[1] + d0[2]*d1[2] >fCutMinDir) return kFALSE;
  //
  return kTRUE;  
}
 

AliTrackReference * AliMaterialBudget::GetFirstTPCTrackRef(AliMCParticle *mcParticle) 
{
  // return first TPC track reference 
  if(!mcParticle) return 0;

  // find first track reference 
  // check direction to select proper reference point for looping tracks
  Int_t nTrackRef = mcParticle->GetNumberOfTrackReferences();
  AliTrackReference *ref = 0;
  AliTrackReference *refIn = 0;
  for (Int_t iref = 0; iref < nTrackRef; iref++) { 
    ref = mcParticle->GetTrackReference(iref);
    if(ref && (ref->DetectorId()==AliTrackReference::kTPC))
    {
      //Float_t dir = ref->X()*ref->Px()+ref->Y()*ref->Py();
      //if(dir < 0.) break;

      refIn = ref;
      break;
    }
  }
  
  return refIn;
}


AliTrackReference *  AliMaterialBudget::GetAllTOFinfo(AliMCParticle *mcParticle, Int_t &nTrackRef, Int_t &nTrackRefITS, Int_t retValue) {
  //
  //
  //

  if(!mcParticle) return 0;
  Int_t counter = 0;
  nTrackRef = 0;
  nTrackRefITS = 0;
  AliTrackReference *ref = 0;
  for (Int_t iref = 0; iref < mcParticle->GetNumberOfTrackReferences(); iref++) { 
    ref = mcParticle->GetTrackReference(iref);
    if(ref && (ref->DetectorId()==AliTrackReference::kTOF)) {
      counter = iref;
      nTrackRef++;
    }
    if(ref && (ref->DetectorId()==AliTrackReference::kITS)) nTrackRefITS++;    
  }
  if (nTrackRef ==0) return 0;
  if (retValue == 0) return mcParticle->GetTrackReference(counter - nTrackRef +1);
  return mcParticle->GetTrackReference(counter);
 
}


void AliMaterialBudget::FinishTaskOutput()
{
  //
  // According description in AliAnalisysTask this method is call
  // on the slaves before sending data
  //
  Terminate("slave");
  gSystem->Exec("pwd");
  RegisterDebugOutput();

}


void AliMaterialBudget::RegisterDebugOutput(){
  //
  //
  //
  //
  // store  - copy debug output to the destination position
  // currently ONLY for local copy
  TString dsName;
  dsName=GetName();
  dsName+="Debug.root";
  dsName.ReplaceAll(" ","");
  TString dsName2=fDebugOutputPath.Data();
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+=gSystem->HostName();
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+="/";
  dsName2+=gSystem->BaseName(gSystem->pwd());
  dsName2+="/";
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+=dsName;
  AliInfo(Form("copy %s\t%s\n",dsName.Data(),dsName2.Data()));
  printf("copy %s\t%s\n",dsName.Data(),dsName2.Data());
  TFile::Cp(dsName.Data(),dsName2.Data());
}


Bool_t AliMaterialBudget::PropagateCosmicToDCA(AliExternalTrackParam *param0, AliExternalTrackParam *param1, Double_t mass){
  //
  // param0 - upper part of cosmic track
  // param1 - lower part of cosmic track
  //
  // 0. Propagate both tracks to DCA to (0,0,0)
  // 1. After propagation to DCA rotate track param1 to corrdinate system of track1 <-> rotate param0 to coordinate system of param 1 ????
  // 2. Propagate track 1 to refernce x from track0
  //

  // step 0.

  Float_t d0  = param0->GetLinearD(0,0);
  Float_t d1  = param1->GetLinearD(0,0);
  Float_t dmax = TMath::Max(TMath::Abs(d0),TMath::Abs(d1));
  //
  // propagate in the beginning taking all material into account
  //
  AliTracker::PropagateTrackToBxByBz(param0,dmax+1.,mass,0.5,kTRUE,0.99,-1.);
  AliTracker::PropagateTrackToBxByBz(param1,dmax+1.,mass,0.5,kTRUE,0.99,1.);
  //
  Double_t vtxx[3]={0,0,0};
  Double_t svtxx[3]={0.000001,0.000001,100.};
  AliESDVertex vtx(vtxx,svtxx);
  //
  Bool_t b0 = param0->PropagateToDCA(&vtx,AliTracker::GetBz(),1000);
  Bool_t b1 = param1->PropagateToDCA(&vtx,AliTracker::GetBz(),1000);

  if (!(b0 && b1)) return 0;

  // step 1.
    
  Float_t dAlpha = param0->GetAlpha();
  param1->Rotate(dAlpha);

  // step 2.

  Float_t refX = param0->GetX();
  param1->PropagateTo(refX,AliTracker::GetBz());

  return kTRUE;
  
  
}


