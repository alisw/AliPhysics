//
// This class is the task to check the 
// Propagation and Update method used in the 
//               1. AliExternalTrackParam 
//               2. AliTracker
//
// Pure Monte-Carlo data used, not influence of detectors

// Input - TParticle + Array of track references - (Points alogn track trajectories)
// Output - Trees with track references - no histograms
//          MC tree -  test for material budget correction 
//                     see function ProcessRefTracker
//          MCupdate tree - test for correctness of propagation and update
//                     see function AliMCTrackingTestTask::FitTrackRefs
//
// Principle - Creates AliExternalTrackParam form 1 Track Refernece - 
//             Propagate it to other
// Magnetic field and the geometry has to be created before using it 

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
#include <TList.h>
#include <TTree.h>
// ALIROOT includes
#include <TTreeStream.h>
#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

#include <AliESD.h>
#include "AliTrackComparison.h"
#include "AliMCTrackingTestTask.h"
#include "AliGenInfoMaker.h"
#include "AliHelix.h"
#include "AliTrackPointArray.h"
#include "AliESDCaloCluster.h"

//
#include "AliMCInfo.h"
#include "AliComparisonObject.h"
#include "AliESDRecInfo.h"
#include "AliTPCParamSR.h"
#include "AliTracker.h"
#include "AliTPCseed.h"

#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGeomManager.h"

// STL includes
#include <iostream>

using namespace std;

ClassImp(AliMCTrackingTestTask)

#define USE_STREAMER 0
#define DEBUG 0

//________________________________________________________________________
AliMCTrackingTestTask::AliMCTrackingTestTask() : 
  AliAnalysisTask(), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fCurrentRun(-1),
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(0),
  fDebugOutputPath(),
  fOutList(NULL),
  fPitList(NULL),
  fCompList(NULL)
{
  //
  // Default constructor (should not be used)
  //
}

AliMCTrackingTestTask::AliMCTrackingTestTask(const AliMCTrackingTestTask& /*info*/) : 
  AliAnalysisTask(), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fCurrentRun(-1),
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(),
  fDebugOutputPath(),
  fOutList(NULL),
  fPitList(NULL),
  fCompList(NULL)
{
  //
  // Default constructor 
  //
}



//________________________________________________________________________
AliMCTrackingTestTask::AliMCTrackingTestTask(const char *name) : 
  AliAnalysisTask(name, "AliMCTrackingTestTask"), 
  fMCinfo(0),     //! MC event handler
  fESD(0),
  fCurrentRun(-1),
  fDebugStreamer(0),
  fStreamLevel(0),
  fDebugLevel(0),
  fDebugOutputPath(),
  fOutList(NULL),
  fPitList(NULL),
  fCompList(NULL)
{
  //
  // Normal constructor
  //
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());

  // create the list for comparison objects
  fCompList = new TList;
  //
  //
}

AliMCTrackingTestTask::~AliMCTrackingTestTask(){
  //
  //
  //
  if (fDebugLevel>0)  printf("AliMCTrackingTestTask::~AliMCTrackingTestTask\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer=0;
  if (fOutList)    delete fOutList;   fOutList   = 0;
  if (fCompList)   delete fCompList;  fCompList = 0; 
}


//_____________________________________________________________________________
Bool_t AliMCTrackingTestTask::AddComparisonObject(AliTrackComparison *cObj) 
{
  // add comparison object to the list
  if(cObj == 0) {
    Printf("ERROR: Could not add comparison object");
    return kFALSE;
  }

  // add object to the list
  fCompList->AddLast(cObj);
       
  return kTRUE;
}
//________________________________________________________________________
void AliMCTrackingTestTask::ConnectInputData(Option_t *) 
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
      fESD = esdH->GetEvent();
      //Printf("*** CONNECTED NEW EVENT ****");
    }  
  }
  AliMCEventHandler* mcinfo = (AliMCEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());  
  mcinfo->SetReadTR(kTRUE);
  
  fMCinfo = mcinfo->MCEvent();
}

//________________________________________________________________________
void AliMCTrackingTestTask::CreateOutputObjects() 
{
  //
  // Connect the output objects
  //
  if(fDebugLevel>3)
    cout << "AnalysisTaskTPCCluster::CreateOutputObjects()" << endl;
  
  if (!fOutList)
    fOutList = new TList;
  fOutList->SetOwner(kTRUE);
  fPitList = fOutList->MakeIterator();
  
  // create output list
  
  // add comparison objects to the output
  AliTrackComparison *cObj=0;
  Int_t count=0;
  TIterator *pitCompList = fCompList->MakeIterator();
  pitCompList->Reset();
  while(( cObj = (AliTrackComparison*)pitCompList->Next()) != NULL) {
    fOutList->Add(cObj);
    count++;
  }
  Printf("UserCreateOutputObjects(): Number of output comparison objects: %d \n", count);
  
  PostData(0, fOutList);  
}


//________________________________________________________________________
void AliMCTrackingTestTask::Exec(Option_t *) {
  //
  // Execute analysis for current event 
  //

#if DEBUG
  printf("New event!\n");
#endif

  if(fDebugLevel>3)
    cout << "AliMCTrackingTestTask::Exec()" << endl;
    


  // If MC has been connected   
  if (!fMCinfo){
    cout << "Not MC info\n" << endl;
  }else{
    ProcessMCInfo();
    // fMCinfo->Print();
    //DumpInfo();
  }

  PostData(0, fOutList);
}      


//________________________________________________________________________
void AliMCTrackingTestTask::Terminate(Option_t *) {
    //
    // Terminate loop
    //
  if(fDebugLevel>3)
    printf("AliMCTrackingTestTask: Terminate() \n");  
  //
  if (fDebugLevel>0) printf("AliMCtrackingTestTask::Terminate\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer = 0;

//   AliTrackComparison *cObj=0;
//   TIterator *pitCompList = fCompList->MakeIterator();
//   pitCompList->Reset();
//   while(( cObj = (AliTrackComparison*)pitCompList->Next()) != NULL) {
//     for(Int_t i=0; i<6; i++)
//       cObj->MakeDistortionMap(438,i);
//   }

  return;
}


//________________________________________________________________________
TTreeSRedirector *AliMCTrackingTestTask::GetDebugStreamer(){
  //
  // Get Debug streamer
  // In case debug streamer not yet initialized and StreamLevel>0 create new one
  //
  if (fStreamLevel==0) return 0;
  if (fDebugStreamer) return fDebugStreamer;
  TString dsName;
  dsName=GetName();
  dsName+="Debug.root";

  printf(" get debug streamer \n");

  dsName.ReplaceAll(" ","");
  fDebugStreamer = new TTreeSRedirector(dsName.Data());
  return fDebugStreamer;
}

//________________________________________________________________________
AliExternalTrackParam * AliMCTrackingTestTask::MakeTrack(const AliTrackReference* ref, TParticle*part)
{
  //
  // Make track out of the track ref
  // part - TParticle used to determine chargr
  // the covariance matrix - equal 0 - starting from ideal MC position
  Double_t xyz[3]={ref->X(),ref->Y(),ref->Z()};
  Double_t pxyz[3]={ref->Px(),ref->Py(),ref->Pz()};
  Double_t cv[21];
  for (Int_t i=0; i<21;i++) cv[i]=0;
  if (!part->GetPDG()) return 0;
  AliExternalTrackParam * param = new AliExternalTrackParam(xyz,pxyz,cv,TMath::Nint(part->GetPDG()->Charge()/3.));
  return param;
}

//________________________________________________________________________
Bool_t  AliMCTrackingTestTask::PropagateToPoint(AliExternalTrackParam *param, Double_t *xyz, Double_t mass, Float_t step){
  // 
  // Propagate track to point xyz using 
  // AliTracker::PropagateToBxByBz functionality
  //
  //  param - track parameters
  //  xyz   - position to propagate
  //  mass  - particle mass
  //  step  - step to be used
  Double_t radius=TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);

  AliTracker::PropagateTrackToBxByBz(param, radius+step, mass, step, kTRUE,0.99,-1);
  param->Rotate(alpha);
  Bool_t isOK = param->PropagateTo(radius,AliTracker::GetBz());

  return isOK;
}

//________________________________________________________________________
void  AliMCTrackingTestTask::ProcessMCInfo(){
  //
  //
  //
   //
#if DEBUG
  printf("ProcessMCInfo\n");
#endif

  const AliESDVertex *pVertex = fESD->GetPrimaryVertex();
  if(!pVertex) AliError("No primary vertex found!\n");
  Double_t vPos[3];
  pVertex->GetXYZ(vPos);

  TParticle * particle= new TParticle();
  TClonesArray * trefs = new TClonesArray("AliTrackReference");
  const Double_t kPcut=0.1;
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


  for (Int_t ipart=0;ipart<npart;ipart++){
    Int_t status = fMCinfo->GetParticleAndTR(ipart, particle, trefs);
    if (status<0 || !particle || !trefs) continue;
    Int_t nref = trefs->GetEntries();
    if (nref<5) continue;
    FitTrackRefs(particle,trefs);

    AliTrackReference * tpcIn=0;
    AliTrackReference * tpcOut=0;
    AliTrackReference * trdIn=0;
    AliTrackReference * trdOut=0;
    AliTrackReference * itsIn=0;
    AliTrackReference * itsOut=0;

    AliTrackReference * tofIn=0;
    AliTrackReference * tofOut=0;
    AliTrackReference * hmpidIn=0;
    AliTrackReference * hmpidOut=0;
    AliTrackReference * emcalIn=0;
    AliTrackReference * emcalOut=0;

    Double_t rmax=0;
    Double_t rmin=1000;
    for (Int_t iref=0;iref<nref;iref++){
      AliTrackReference * ref = (AliTrackReference*)trefs->At(iref);
      if (!ref) continue;
      
      Float_t dir = ref->X()*ref->Px()+ref->Y()*ref->Py();
      
      if (dir<0) break; // oposite direction - looping track - return back
      if (ref->P()<kPcut) continue;
      if (ref->R()<rmax) break;
      //if (ref->R()<rmin)  break; 
      //TPC
      if (ref->DetectorId()==AliTrackReference::kTPC){
	if (!tpcIn) {
	  tpcIn  = ref;
	}else{
	  if (ref->R()>tpcIn->R()) tpcOut = ref;
	}	
      }
      //ITS
      if (ref->DetectorId()==AliTrackReference::kITS){
	if (!itsIn) {
	  itsIn  = ref;
	}else{
	  if (ref->R()>itsIn->R()) itsOut = ref;
	}	
      }
      //TRD
      if (ref->DetectorId()==AliTrackReference::kTRD){
	if (!trdIn) {
	  trdIn  = ref;
	}else{
	  if (ref->R()>trdIn->R()) trdOut = ref;
	}	
      }      
      //TOF
      if (ref->DetectorId()==AliTrackReference::kTOF){
	if (!tofIn) {
	  tofIn  = ref;
	}else{
	  if (ref->R()>tofIn->R()) tofOut = ref;
	}	
      }      

      //HMPID
      if (ref->DetectorId()==AliTrackReference::kHMPID){
	if (!hmpidIn) {
	  hmpidIn  = ref;
	}else{
	  if (ref->R()>hmpidIn->R()) hmpidOut = ref;
	}	
      }      


      //EMCAL
      if (ref->DetectorId()==AliTrackReference::kEMCAL){
	if (!emcalIn) {
	  emcalIn  = ref;
	}else{
	  if (ref->R()>emcalIn->R()) emcalOut = ref;
	}	
      }      

      if (ref->R()<rmin) rmin=ref->R();
      if (ref->R()>rmax) rmax=ref->R();
    } // end ref loop

    // -----------------------------------------------
    // check for gammas
    
    // electron
    if ( TMath::Abs(particle->GetPdgCode()) == 11 ) {

      // findable
      if (IsFindable(ipart,70)) {

	// is from gamma conversion
	Int_t motherId = particle->GetFirstMother();
	if (motherId > 0) {
	  if (motherId < npart) {
	    TParticle* mother = (fMCinfo->Stack())->Particle(motherId);
	    if (mother && TMath::Abs(mother->GetPdgCode()) == 22) {
	      Int_t nDaughters = mother->GetNDaughters();
	      
	      for (Int_t idx=0; idx<nDaughters; ++idx) {
		Int_t daughterId = mother->GetDaughter(idx);
		if ( daughterId == ipart || daughterId >= npart )
		  continue;

		TParticle* daughter = (fMCinfo->Stack())->Particle(daughterId);
		if (daughter && TMath::Abs(daughter->GetPdgCode()) == 11) {
		  //Bool_t findable = IsFindable(daughterId,70);
#if USE_STREAMER
		  TTreeSRedirector *pcstream = GetDebugStreamer();
		  if (pcstream){
		    (*pcstream)<<"MCgamma"<<
		    "triggerP="<<particle<<      // trigger electron
		      "motherP="<<mother<<         // mother gamma
		      "daughterP="<<daughter<<     // daughter electron
		      "isFindable="<<findable<<     // 2nd is findable
		      "\n";
		  }
#endif		
		} 
	      }
	    }  
	  }
	}
      }
    }

    // -----------------------------------------------
    if (tpcIn && tpcOut) {
      ProcessRefTracker(tpcIn,tpcOut,particle,1);
      ProcessRefTracker(tpcIn,tpcOut,particle,3);
    }
    if (itsIn  && itsOut)  ProcessRefTracker(itsIn, itsOut, particle, 0);
    if (trdIn  && trdOut)  ProcessRefTracker(trdIn, trdOut, particle, 2);

    if (tpcOut && trdIn)   ProcessRefTracker(tpcOut,trdIn,  particle, 4);
    if (tpcOut && tofIn)   ProcessRefTracker(tpcOut,tofIn,  particle, 5);
    if (tpcOut && hmpidIn) ProcessRefTracker(tpcOut,hmpidIn,particle, 6);
    if (tpcOut && emcalIn) ProcessRefTracker(tpcOut,emcalIn,particle, 7);

    if (tpcOut && trdOut)  ProcessRefTracker(tpcOut,trdOut, particle, 8);
    if (trdIn  && tofIn)   ProcessRefTracker(trdIn, tofIn,  particle, 9);
    if (tofIn  && tofOut)  ProcessRefTracker(tofIn, tofOut, particle,10);


    // -----------------------------------------------
    //Test for new base tracking class
          
    AliTrackComparison *cObj=0;
    AliTrackPoint *point=new AliTrackPoint();
    AliESDCaloCluster *cluster=0;
    AliExternalTrackParam *tpcOuter=0;

    Double_t eMax=0;
    Int_t clsIndex=-1;

    for(Int_t iCluster=0; iCluster<fESD->GetNumberOfCaloClusters(); iCluster++)
      {
	cluster = fESD->GetCaloCluster(iCluster);
	if(!cluster->IsEMCAL()) continue;
	if(cluster->GetLabel()==ipart && cluster->E()>eMax)
	  {               
	    clsIndex=iCluster;
	    eMax=cluster->E();
	  }
      }
    
    if(clsIndex>-1)
      {
	if(cluster) cluster=0;
	cluster = fESD->GetCaloCluster(clsIndex);
	Float_t clusterPos[3];
	cluster->GetPosition(clusterPos);
	point->SetXYZ(clusterPos[0],clusterPos[1],clusterPos[2],0);
      }

    if(tpcOut)
      tpcOuter = MakeTrack(tpcOut,particle);
      

    Double_t mass = particle->GetMass();
    Int_t charge = TMath::Nint(particle->GetPDG()->Charge()/3.);
    fPitList->Reset();
    while(( cObj = (AliTrackComparison *)fPitList->Next()) != NULL) {
      TString objName(cObj->GetName());
      if(!objName.CompareTo("TPCOutToEMCalInElecCls"))
        { 
          if(TMath::Abs(particle->GetPdgCode())==11 && tpcOuter && point && cluster && emcalIn)
	    {
	      printf("\n\nTPCOutToEMCalInElecCls: ");
	      cout<<cObj->AddTracks(tpcOuter,point,mass,cluster->E(),vPos)<<endl;
	    }
        }
          
      if(!objName.CompareTo("TPCOutToEMCalInElec"))
        {
          if(TMath::Abs(particle->GetPdgCode())==11 && tpcOut && emcalIn)
	    {
	       printf("TPCOutToEMCalInElec: ");
	       cout<<cObj->AddTracks(tpcOut,emcalIn,mass,charge)<<endl;
	    }
        }

      if(!objName.CompareTo("TPCOutToEMCalInPion"))
        {
          if(TMath::Abs(particle->GetPdgCode())==211 && tpcOut && emcalIn)
            cObj->AddTracks(tpcOut,emcalIn,mass,charge);
        }

      if(!objName.CompareTo("TPCOutToTOFIn"))
        {
          if(tpcOut && tofIn)
            cObj->AddTracks(tpcOut,tofIn,mass,charge);
        }

      if(!objName.CompareTo("TPCOutToHMPIDIn"))
        {
          if(tpcOut && hmpidIn)
            cObj->AddTracks(tpcOut,hmpidIn,mass,charge);
        }
    }   
    //End of the test for new base tracking class
    // -----------------------------------------------
    delete point;

  }

  trefs->Clear("C");
  //delete particle;
  //delete tpcIn;

}

void AliMCTrackingTestTask::ProcessRefTracker(AliTrackReference* refIn,  AliTrackReference* refOut, TParticle*part,Int_t type){
  //
  // Test propagation from In to out
  //

#if DEBUG
  printf("ProcessRefTracker\n");
#endif

  AliExternalTrackParam *param = 0;
  AliExternalTrackParam *paramMC = 0;
  AliExternalTrackParam *paramDebug = 0;

  //  Double_t xyzIn[3]={refIn->X(),refIn->Y(), refIn->Z()};
  Double_t xyzOut[3]={refOut->X(),refOut->Y(), refOut->Z()};
  Double_t mass = part->GetMass();
  Double_t step=1;
  //
  param=MakeTrack(refIn,part);
  paramMC=MakeTrack(refOut,part);
  paramDebug=MakeTrack(refIn,part);

  if (!param) return;
  if (type!=3) PropagateToPoint(param,xyzOut, mass, step);
  //
#if 0
  /*
    if (type==3) {
    AliTPCseed seed;
    seed.Set(param->GetX(),param->GetAlpha(),param->GetParameter(),param->GetCovariance());
    Float_t alpha= TMath::ATan2(refIn->Y(),refIn->X());
    seed.Rotate(alpha-seed.GetAlpha());
    seed.SetMass(mass);
    for (Float_t xlayer= seed.GetX(); xlayer<refOut->R(); xlayer+=step){
      seed.PropagateTo(xlayer);
    }
    seed.PropagateTo(refOut->R());
    param->Set(seed.GetX(),seed.GetAlpha(),seed.GetParameter(),seed.GetCovariance());
  }
*/
#endif
#if USE_STREAMER
  TTreeSRedirector *pcstream = GetDebugStreamer();
  TVectorD gpos(3);
  TVectorD gmom(3);
  Bool_t isOK=kTRUE;
  isOK&=param->Rotate(paramMC->GetAlpha());
  isOK&=param->PropagateTo(paramMC->GetX(),AliTracker::GetBz());
  param->GetXYZ(gpos.GetMatrixArray());
  param->GetPxPyPz(gmom.GetMatrixArray());
  if (pcstream){
    (*pcstream)<<"MC"<<
      "isOK="<<isOK<<
      "type="<<type<<              // detector matching type
      "step="<<step<<              // propagation step length
      "refIn.="<<refIn<<           // starting track refernce
      "refOut.="<<refOut<<         // outer track reference
      "p.="<<part<<                // particle desription (TParticle)
      "par.="<<param<<             // AliExternalTrackParam create at starting point propagated to outer track ref radius
      "parMC.="<<paramMC<<         // AliExternalTrackParam created at the outer point  
      "gpos.="<<&gpos<<            // global position
      "gmom.="<<&gmom<<            // global momenta
      "\n";
  }
#endif
  delete param;
  delete paramMC;
  delete paramDebug;
}

 
Bool_t AliMCTrackingTestTask::IsFindable(Int_t label, Float_t minTrackLength ) {
  //
  // Find findable tracks
  //
  
  AliMCParticle *mcParticle = (AliMCParticle*) fMCinfo->GetTrack(label);
  if(!mcParticle) return kFALSE;
  
  Int_t counter; 
  Float_t tpcTrackLength = mcParticle->GetTPCTrackLength(AliTracker::GetBz(),0.05,counter,3.0); 
  //printf("tpcTrackLength %f \n", tpcTrackLength);
  
  return (tpcTrackLength>minTrackLength);    
}


void  AliMCTrackingTestTask::FitTrackRefs(TParticle * part, TClonesArray * trefs){
  //
  //
  //
  //

#if DEBUG
  printf("FitTrackRefs\n");
#endif

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
  //AliTrackReference * refOut = (AliTrackReference*)trefs->At(iref1);
  AliExternalTrackParam *paramPropagate= MakeTrack(refIn,part);
  AliExternalTrackParam *paramUpdate   = MakeTrack(refIn,part);
  paramUpdate->AddCovariance(covar);
  //Double_t mass = part->GetMass();
  //Double_t charge = part->GetPDG()->Charge()/3.;
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
#if USE_STREAMER
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
#endif
  delete paramPropagate;
  delete paramUpdate;
}




void AliMCTrackingTestTask::FinishTaskOutput()
{
  //
  // According description in AliAnalisysTask this method is call
  // on the slaves before sending data
  //
  Terminate("slave");
  gSystem->Exec("pwd");
  RegisterDebugOutput();

}


void AliMCTrackingTestTask::RegisterDebugOutput(){
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


/*
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chain = tool.MakeChain("mctracking.txt","MC",0,100);
  chain->Lookup();
  //
  //
  chain->SetAlias("pdg","(p.fPdgCode)");
  chain->SetAlias("dPRec","(refOut.P()-par.P())/refIn.P()");
  chain->SetAlias("dPMC","(refOut.P()-refIn.P())/refIn->P()");
  chain->SetAlias("dPtRec","(refOut.Pt()-par.Pt())/refIn.Pt()");
  chain->SetAlias("dPtMC","(refOut.Pt()-refIn.Pt())/refIn->Pt()");


  // ITS
  chain->Draw("-sqrt(dPRec):-sqrt(dPMC)","abs(pdg)!=11&&type==0&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  htemp->SetYTitle("#sqrt{#DeltaP_{rec}/P}");
  htemp->SetXTitle("#sqrt{#DeltaP_{mc}/P}");
  gPad->SaveAs("picLoss/dPcorr_ITS_step1.gif");
  gPad->SaveAs("picLoss/dPcorr_ITS_step1.eps");
  // TPC
  chain->Draw("-sqrt(dPRec):-sqrt(dPMC)","abs(pdg)!=11&&type==1&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  htemp->SetYTitle("#sqrt{#DeltaP_{rec}/P}");
  htemp->SetXTitle("#sqrt{#DeltaP_{mc}/P}");
  gPad->SaveAs("picLoss/dPcorr_TPC_step1.gif");
  gPad->SaveAs("picLoss/dPcorr_TPC_step1.eps");
  //
   // TPC
  chain->Draw("-sqrt(dPRec):-sqrt(dPMC)","abs(pdg)!=11&&type==3&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  htemp->SetYTitle("#sqrt{#DeltaP_{rec}/P}");
  htemp->SetXTitle("#sqrt{#DeltaP_{mc}/P}");
  gPad->SaveAs("picLoss/dPcorr_TPCseed_step1.gif");
  gPad->SaveAs("picLoss/dPcorr_TPCseed_step1.eps");


  // TRD
  chain->Draw("-sqrt(dPRec):-sqrt(dPMC)","abs(pdg)!=11&&type==2&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  htemp->SetYTitle("#sqrt{#DeltaP_{rec}/P}");
  htemp->SetXTitle("#sqrt{#DeltaP_{mc}/P}");
  gPad->SaveAs("picLoss/dPcorr_TRD_step1.gif");
  gPad->SaveAs("picLoss/dPcorr_TRD_step1.eps");

  //
  //
  //
  chain->Draw("(par.Pt()-refIn.Pt())/refIn.Pt()>>his(100,-0.02,0.02)","abs(pdg)!=11&&type==3&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  his->SetXTitle("(P_{trec}-P_{tmc})/P_{tmc}");
  gPad->SaveAs("picLoss/dPtcorr_TPCseed_step1_1D.eps");
  gPad->SaveAs("picLoss/dPtcorr_TPCseed_step1_1D.gif");

  chain->Draw("(par.P()-refIn.P())/refIn.P()>>his(100,-0.02,0.02)","abs(pdg)!=11&&type==3&&p.Pt()<0.5&&abs(p.R())<1&&abs(refOut.fZ)<220");
  his->SetXTitle("(P_{rec}-P_{mc})/P_{mc}");
  gPad->SaveAs("picLoss/dPcorr_TPCseed_step1_1D.eps");
  gPad->SaveAs("picLoss/dPcorr_TPCseed_step1_1D.gif");
*/
