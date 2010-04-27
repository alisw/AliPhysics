// $Id: AliHLTD0Trigger.cxx 
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Ovrebekk                                        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTD0Trigger.cxx
/// @author Gaute Ovrebekk
/// @date   2009-10-28
/// @brief  HLT trigger component for D0->Kpi

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTD0Trigger.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "TH1F.h"
#include "AliHLTD0toKpi.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliAODRecoDecay.h"
#include "AliHLTMCEvent.h"
#include "AliStack.h"
#include "TParticle.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTD0Trigger)

AliHLTD0Trigger::AliHLTD0Trigger()
  : AliHLTTrigger()
  , fPtMin(0.0)
  , fdca(0.0)
  , finvMass(0.0)     
  , fcosThetaStar(0.0)
  , fd0(0.0) 
  , fd0d0(0.0)
  , fcosPoint(0.0)
  , fplothisto(false)
  , fUseV0(false)
  , mD0PDG(TDatabasePDG::Instance()->GetParticle(421)->Mass())
  , fD0mass(NULL)
  , fD0pt(NULL)
  , fPos()
  , fNeg()
  , fd0calc(NULL)                  
  , ftwoTrackArray(NULL)
  , fTotalD0(0)
  , fTotalD0true(0)
  , fVertex(NULL)
  , fField(0)
  , fEvent(NULL)
  , fuseKF(false)
{
  
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

const char* AliHLTD0Trigger::fgkOCDBEntry="HLT/ConfigHLT/D0Trigger";

AliHLTD0Trigger::~AliHLTD0Trigger()
{
  //if(fd0calc){delete fd0calc;}  
  //if(fD0mass){delete fD0mass;}
  //if(ftwoTrackArray){delete ftwoTrackArray;}
  // see header file for class documentation
}
  
const char* AliHLTD0Trigger::GetTriggerName() const
{
  // see header file for class documentation
  return "D0Trigger";
}

AliHLTComponent* AliHLTD0Trigger::Spawn()
{
  // see header file for class documentation
  return new AliHLTD0Trigger;
}

int AliHLTD0Trigger::DoTrigger()
{
  // -- Iterator over Data Blocks --
  //const AliHLTComponentBlockData* iter = NULL;

  if (!IsDataEvent()) return 0;

  Int_t nD0=0;
  Int_t nD0true=0;
  TString description;

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeMCObject|kAliHLTDataOriginHLT); iter != NULL; iter = GetNextInputObject() ) {   
    fEvent = dynamic_cast<AliHLTMCEvent*>(const_cast<TObject*>( iter ) );
    if ( ! fEvent ) {
      HLTError( "No MC Event present!" );
      break;
    }
  }

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {   
    if(fUseV0){
      nD0=RecV0(iter);
    }
    else{
       AliESDEvent *event = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
       event->GetStdContent();
       fField = event->GetMagneticField();
       const AliESDVertex* pv = event->GetPrimaryVertexTracks();
       fVertex =  new AliESDVertex(*pv);
       if(fVertex->GetNContributors()<2){
	 HLTWarning("Contributors in ESD vertex to low or not been set");
	 continue;
       }
       for(Int_t it=0;it<event->GetNumberOfTracks();it++){
	 SingleTrackSelect(event->GetTrack(it));
       }
    
       RecD0(nD0,nD0true);       
    }
  }

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDVertex|kAliHLTDataOriginOut); 
	iter != NULL; iter = GetNextInputObject() ) { 
    fVertex = dynamic_cast<AliESDVertex*>(const_cast<TObject*>( iter ));
    if(!fVertex){
      HLTError("Vertex object is corrupted");
      //iResult = -EINVAL;    
    }
  }  
   
  for ( const AliHLTComponentBlockData* iter = GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginITS); 
	iter != NULL; iter = GetNextInputBlock() ) { 
    vector<AliHLTGlobalBarrelTrack> tracksVector;
    AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(iter->fPtr), iter->fSize, tracksVector);
    
    fField = GetBz();
    for(UInt_t i=0;i<tracksVector.size();i++){
      SingleTrackSelect(&tracksVector[i]);
    }
    RecD0(nD0,nD0true);
  }    

  fTotalD0+=nD0;
  fTotalD0true += nD0true;

  ftwoTrackArray->Clear();
  fPos.clear();
  fNeg.clear();
     
  HLTDebug("Number of D0 found: %d",nD0);
  HLTDebug("Number of True D0 found: %d",nD0true);
  HLTDebug("Total Number of D0 found: %d",fTotalD0);
  HLTDebug("Total Number of True D0 found: %d",fTotalD0true);

  if(fplothisto){
    PushBack( (TObject*) fD0mass, kAliHLTDataTypeHistogram,0);
    PushBack( (TObject*) fD0pt, kAliHLTDataTypeHistogram,0);
  }
  
  //if (iResult>=0) {
  if (1) {
   
    if (nD0>=1) {
      description.Form("Event contains %d D0(s)", nD0);
      SetDescription(description.Data());
      // Enable the central detectors for readout.
      GetReadoutList().Enable(
			      AliHLTReadoutList::kITSSPD |
			      AliHLTReadoutList::kITSSDD |
			      AliHLTReadoutList::kITSSSD |
			      AliHLTReadoutList::kTPC |
			      AliHLTReadoutList::kTRD |
			      AliHLTReadoutList::kTOF |
			      AliHLTReadoutList::kHMPID |
			      AliHLTReadoutList::kPHOS
			      );
      // Add the available HLT information for readout too.
      GetTriggerDomain().Add("CLUSTERS", "TPC ");
      TriggerEvent(true);
      return 0;
    }
    description.Form("No D0");
  } else {
    description.Form("No input blocks found");
  }
  SetDescription(description.Data());
  TriggerEvent(false);
  //return iResult;
  return 0;
}

int AliHLTD0Trigger::DoInit(int argc, const char** argv)
{
  
  fd0calc = new AliHLTD0toKpi();
  ftwoTrackArray = new TObjArray(2);

  fplothisto=false;
  // see header file for class documentation
  fD0mass = new TH1F("hMass","D^{0} mass plot",100,1.7,2);
  fD0pt = new TH1F("hPt","D^{0} Pt plot",20,0,20);
  // first configure the default
  int iResult=0;
  if (iResult>=0) iResult=ConfigureFromCDBTObjString(fgkOCDBEntry);

  // configure from the command line parameters if specified
  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);
  return iResult;
}

int AliHLTD0Trigger::DoDeinit()
{
  // see header file for class documentation
  if(fd0calc){delete fd0calc;}  
  if(fD0mass){delete fD0mass;}
  if(ftwoTrackArray){delete ftwoTrackArray;}
  if(fVertex){delete fVertex;}
  return 0;
}

int AliHLTD0Trigger::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{
  // see header file for class documentation

  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if (!entry || entry[0]==0) {
    entry=fgkOCDBEntry;
  }

  return ConfigureFromCDBTObjString(entry);
}

int AliHLTD0Trigger::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];
  // -minpt for decay
  if (argument.CompareTo("-pt")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fPtMin=argument.Atof();
    return 2;
  }    
  // minimum dca for decay tracks
  if (argument.CompareTo("-dca")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fdca=argument.Atof();
    return 2;
  }
  // inv. mass half width.
  if (argument.CompareTo("-invmass")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    finvMass=argument.Atof();
    return 2;
  }

// cos theta for decay angle
  if (argument.CompareTo("-costhetastar")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fcosThetaStar=argument.Atof();
    return 2;
  }

  // impact parameter for decay
  if (argument.CompareTo("-d0")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fd0=argument.Atof();
    return 2;
  }
  // product of impact parameter
  if (argument.CompareTo("-d0d0")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fd0d0=argument.Atof();
    return 2;
  }
  // product of impact parameter
  if (argument.CompareTo("-cospoint")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fcosPoint=argument.Atof();
    return 2;
  }
  if (argument.CompareTo("-plothistogram")==0) {
    fplothisto=true;
    return 1;
  }
  if (argument.CompareTo("-usev0")==0) {
    fUseV0=true;
    return 1;
  }
  if (argument.CompareTo("-useKF")==0) {
    fuseKF=true;
    return 1;
  }
  // unknown argument
  return -EINVAL;
}

void AliHLTD0Trigger::SingleTrackSelect(AliExternalTrackParam* t){
  // Offline har || på disse kuttene på de to henfallsproduktene 
  Double_t pv[3];
  fVertex->GetXYZ(pv);

  if(t->Pt()<fPtMin){return;}
  if(TMath::Abs(t->GetD(pv[0],pv[1],fField)) > fd0){return;}

  if(t->Charge()>0){
    fPos.push_back(t);
  }
  else{
    fNeg.push_back(t);
  }
}

void AliHLTD0Trigger::RecD0(Int_t& nD0, Int_t& nD0true){
 
  Double_t D0,D0bar,xdummy,ydummy; 
  Double_t d0[2];
  Double_t svpos[3];
  Double_t pvpos[3];
  
  if(!fVertex){
    HLTError("No Vertex is set");
    return;
  }
  fVertex->GetXYZ(pvpos);
    
  for(UInt_t ip=0;ip<fPos.size();ip++){
    AliExternalTrackParam *tP=fPos[ip];
    for(UInt_t in=0;in<fNeg.size();in++){
      AliExternalTrackParam *tN=fNeg[in];
          
      tP->PropagateToDCA(fVertex,fField,kVeryBig);  //do I need this??????
      tN->PropagateToDCA(fVertex,fField,kVeryBig);
      
      Double_t dcatPtN = tP->GetDCA(tN,fField,xdummy,ydummy);
      if(dcatPtN>fdca) { continue; }
      
      ftwoTrackArray->AddAt(tP,0);
      ftwoTrackArray->AddAt(tN,1);
      AliAODVertex *vertexp1n1 = fd0calc->ReconstructSecondaryVertex(ftwoTrackArray,fField,fVertex,fuseKF);
      if(!vertexp1n1) { 
	ftwoTrackArray->Clear();
	continue; 
      }
      
      vertexp1n1->GetXYZ(svpos);
      
      tP->PropagateToDCA(vertexp1n1,fField,kVeryBig); 
      tN->PropagateToDCA(vertexp1n1,fField,kVeryBig);
      
      /*
      Double_t px[2],py[2],pz[2];
      Double_t momentum[3];
      tP->GetPxPyPz(momentum);
      px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2];
      tN->GetPxPyPz(momentum);
      px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2];

      Short_t dummycharge=0;
      Double_t *dummyd0 = new Double_t[2];
      Int_t nprongs = 2;

      for(Int_t ipr=0;ipr<nprongs;ipr++) dummyd0[ipr]=0.;
      AliAODRecoDecay *rd = new AliAODRecoDecay(0x0,nprongs,dummycharge,px,py,pz,dummyd0);
      delete [] dummyd0; dummyd0=NULL;

      UInt_t pdg2[2],pdg2bar[2];
      Double_t mPDG,minv,minvbar;

      pdg2[0]=211; pdg2[1]=321; pdg2bar[0]=321; pdg2bar[1]=211;
      minv = rd->InvMass(nprongs,pdg2);
      minvbar = rd->InvMass(nprongs,pdg2bar);
      if(TMath::Abs(minv-mD0PDG)>finvMass && TMath::Abs(minv-mD0PDG)>finvMass) {continue; delete vertexp1n1; delete rd;}
      */
      if((TMath::Abs(fd0calc->InvMass(tN,tP)-mD0PDG)) > finvMass && TMath::Abs((fd0calc->InvMass(tP,tN))-mD0PDG) > finvMass){continue;}
      fd0calc->cosThetaStar(tN,tP,D0,D0bar);
      if(TMath::Abs(D0) > fcosThetaStar && TMath::Abs(D0bar) > fcosThetaStar){continue;}
      d0[0] = tP->GetD(pvpos[0],pvpos[1],fField);
      d0[1] = tN->GetD(pvpos[0],pvpos[1],fField);
      if((d0[0]*d0[1]) > fd0d0){continue;}
      if(fd0calc->pointingAngle(tN,tP,pvpos,svpos) < fcosPoint){continue;}
      
      if(fplothisto){
	//fD0mass->Fill(minv);
	//fD0mass->Fill(minvbar);
	fD0mass->Fill(fd0calc->InvMass(tN,tP));
	fD0mass->Fill(fd0calc->InvMass(tP,tN));
	fD0pt->Fill(fd0calc->Pt(tP,tN));
	/*
	if((fd0calc->InvMass(tN,tP) - mD0PDG) > finvMass){
	  fD0mass->Fill(fd0calc->InvMass(tN,tP));
	}
	else{
	  fD0mass->Fill(fd0calc->InvMass(tP,tN));
	  }
	*/
      }
      
      if(CheckTrackMC(tP,tN)){
	nD0true++;
      }

      nD0++;
      delete vertexp1n1;
    }
  }
}

Int_t AliHLTD0Trigger::RecV0(const TObject* iter){
  int nD0=0;
  Double_t d0[2];
  Double_t D0,D0bar; 
  Double_t svpos[3];
  Double_t pvpos[3];
  
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
  event->GetStdContent();
  Int_t nV0 = event->GetNumberOfV0s();
  Double_t field = event->GetMagneticField();
  const AliESDVertex* pv = event->GetPrimaryVertexTracks();
  pv->GetXYZ(pvpos);

  for (Int_t iv=0; iv<nV0; iv++) {
      
      AliESDtrack *tN=event->GetTrack( event->GetV0(iv)->GetNindex());
      AliESDtrack *tP=event->GetTrack( event->GetV0(iv)->GetPindex());      
      
      if(tN->Pt()<fPtMin && tP->Pt()<fPtMin){continue;} //||????
      
      d0[0] = tP->GetD(pvpos[0],pvpos[1],field);
      d0[1] = tN->GetD(pvpos[0],pvpos[1],field);
      
      if(d0[0]>fd0 && d0[0]>fd0){continue;} // make sure < or>
      
      event->GetV0(iv)->GetXYZ(svpos[0],svpos[1],svpos[2]);
      
      if(!tN->PropagateTo(svpos[0],field) && !tP->PropagateTo(svpos[0],field)){
      	HLTInfo("Tracks could not be propagated to secondary vertex");
      	continue;
      }
      
      Double_t tmp1, tmp2;
      if(tN->GetDCA(tP,field,tmp1,tmp2) > fdca){continue;}
      
      if((fd0calc->InvMass(tN,tP) - mD0PDG) > finvMass && (fd0calc->InvMass(tP,tN) - mD0PDG) > finvMass){continue;}
      fd0calc->cosThetaStar(tN,tP,D0,D0bar);
      if(D0 > fcosThetaStar && D0bar > fcosThetaStar){continue;}
      if((d0[0]*d0[1]) > fd0d0){continue;}
      if(fd0calc->pointingAngle(tN,tP,pvpos,svpos) < fcosPoint){continue;}
      
      nD0++;
      if((fd0calc->InvMass(tN,tP) - mD0PDG) > finvMass){
	fD0mass->Fill(fd0calc->InvMass(tN,tP));
      }
      else{
	fD0mass->Fill(fd0calc->InvMass(tP,tN));
      }
  }
  return nD0;
  
}
bool AliHLTD0Trigger::CheckTrackMC(AliExternalTrackParam* pt, AliExternalTrackParam* pn){

  if(!fEvent){return false;}

  int lP = pt->GetLabel();
  int lN = pn->GetLabel();

  if(lN>=0 && lP>=0){
    
    int imP = (fEvent->Particle(lP))->GetFirstMother();
    int imN = (fEvent->Particle(lN))->GetFirstMother();
    
    if(imP>=0 && imN>=0){
      TParticle * mP = fEvent->Particle(imP);
      TParticle * mN = fEvent->Particle(imN);
      if(fabs(mP->GetPdgCode())==421 && fabs(mN->GetPdgCode())==421 && imP == imN){
	return true;
      }
    }
  }
  return false;
}
