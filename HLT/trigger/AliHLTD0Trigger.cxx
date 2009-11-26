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
  , fPos()
  , fNeg()
  , fd0calc(NULL)                  
  , ftwoTrackArray(NULL)
  , fTotalD0(0)
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
  TString description;

    HLTDebug("Cuts: -pt:%f -dca:%f -invmass:%f -costhetastar:%f -d0:%f -d0d0:%f -cospoint:%f",fPtMin,fdca,finvMass,fcosThetaStar,fd0,fd0d0,fcosPoint);

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {   
    if(fUseV0){
      nD0=RecV0(iter);
    }
    else{
      nD0=RecESDTracks(iter);
    }
  }

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeTrack); iter != NULL; iter = GetNextInputObject() ) { 
    nD0=RecBarrelTracks(iter);
  }    

  fTotalD0+=nD0;
   
  ftwoTrackArray->Clear();
  fPos.clear();
  fNeg.clear();
  
  HLTDebug("Number of D0 found: %d",nD0);
  HLTDebug("Total Number of D0 found: %d",fTotalD0);

  if(fplothisto){PushBack( (TObject*) fD0mass, kAliHLTDataTypeHistogram,0);}
  
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

  // unknown argument
  return -EINVAL;
}

void AliHLTD0Trigger::SingleTrackSelect(AliESDtrack* t, Double_t b,Double_t* pv){
  // Offline har || på disse kuttene på de to henfallsproduktene 
    
  if(t->Pt()<fPtMin){return;}
  if(TMath::Abs(t->GetD(pv[0],pv[1],b)) > fd0){return;}

  if(t->Charge()>0){
    fPos.push_back(t);
  }
  else{
    fNeg.push_back(t);
  }
}

Int_t AliHLTD0Trigger::RecESDTracks(const TObject* iter){
  
  int nD0=0;
  Double_t D0,D0bar,xdummy,ydummy; 
  Double_t d0[2];
  Double_t svpos[3];
  Double_t pvpos[3];
  
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
  event->GetStdContent();
  Double_t field = event->GetMagneticField();
  const AliESDVertex* pv = event->GetPrimaryVertexTracks();
  pv->GetXYZ(pvpos);
  
  for(Int_t it=0;it<event->GetNumberOfTracks();it++){
    SingleTrackSelect(event->GetTrack(it),field,pvpos);
  }
  
  for(Int_t ip=0;ip<fPos.size();ip++){
    AliESDtrack *tP=fPos[ip];
    for(Int_t in=0;in<fNeg.size();in++){
      AliESDtrack *tN=fNeg[in];
      
      tP->PropagateToDCA(pv,field,kVeryBig);  //do I need this??????
      tN->PropagateToDCA(pv,field,kVeryBig);
      
      Double_t dcatPtN = tP->GetDCA(tN,field,xdummy,ydummy);
      if(dcatPtN>fdca) { continue; }
      
      ftwoTrackArray->AddAt(fPos[ip],0);
      ftwoTrackArray->AddAt(fNeg[in],1);
      AliAODVertex *vertexp1n1 = fd0calc->ReconstructSecondaryVertex(ftwoTrackArray,field,pv);
      if(!vertexp1n1) { 
	ftwoTrackArray->Clear();
	continue; 
      }
      
      vertexp1n1->GetXYZ(svpos);
      
      tP->PropagateToDCA(vertexp1n1,field,kVeryBig); 
      tN->PropagateToDCA(vertexp1n1,field,kVeryBig);
      
      if((TMath::Abs(fd0calc->InvMass(tN,tP)-mD0PDG)) > finvMass && TMath::Abs((fd0calc->InvMass(tP,tN))-mD0PDG) > finvMass){continue;}
      fd0calc->cosThetaStar(tN,tP,D0,D0bar);
      if(TMath::Abs(D0) > fcosThetaStar && TMath::Abs(D0bar) > fcosThetaStar){continue;}
      d0[0] = tP->GetD(pvpos[0],pvpos[1],field);
      d0[1] = tN->GetD(pvpos[0],pvpos[1],field);
      if((d0[0]*d0[1]) > fd0d0){continue;}
      if(fd0calc->pointingAngle(tN,tP,pvpos,svpos) < fcosPoint){continue;}
      
      if(fplothisto){
	if((fd0calc->InvMass(tN,tP) - mD0PDG) > finvMass){
	  fD0mass->Fill(fd0calc->InvMass(tN,tP));
	}
	else{
	  fD0mass->Fill(fd0calc->InvMass(tP,tN));
	}
      }
      nD0++;
      delete vertexp1n1;
    }
  }
  return nD0;
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

Int_t AliHLTD0Trigger::RecBarrelTracks(const TObject* iter){

  return 0;

}
