// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
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

/// @file   AliHLTTriggerCosmics.cxx
/// @author Kalliopi Kanaki
/// @date   2011-02-25
/// @brief  HLT trigger component for tagging cosmics tracks in the TPC

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTriggerCosmics.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTErrorGuard.h"
#include "AliTPCcalibTime.h"
#include "AliTracker.h"

#include "TObjArray.h"
#include "TObjString.h"
#include "TDatabasePDG.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerCosmics)

AliHLTTriggerCosmics::AliHLTTriggerCosmics()
  : AliHLTTrigger()
  , fName()
  , fTrackSelection()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

const char* AliHLTTriggerCosmics::fgkDefaultOCDBEntry="HLT/ConfigHLT/CosmicsTrigger";

AliHLTTriggerCosmics::~AliHLTTriggerCosmics(){
// see header file for class documentation
}

const char* AliHLTTriggerCosmics::GetTriggerName() const{
// see header file for class documentation

  if (!fName.IsNull())
    return fName.Data();
  else
    return "CosmicsTrigger";
}

AliHLTComponent* AliHLTTriggerCosmics::Spawn(){
  // see header file for class documentation
  return new AliHLTTriggerCosmics;
}

int AliHLTTriggerCosmics::DoTrigger(){
// see header file for class documentation
 
  if (!IsDataEvent()) {
    IgnoreEvent();  // dont generate any trigger decision.
  }

  int iResult=0;
  int numberOfCosmics=-1;

  const TObject *obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent *esd   = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  if(!esd){
     printf("Empty ESD\n"); 
     return 0;
  }
  esd->GetStdContent();
  
  Int_t startTime = esd->GetTimeStamp()-60*60*1;  //Start time one hour before first event, will make precise cuts later.
  Int_t   endTime = esd->GetTimeStamp()+60*60*23; //End time 23 hours after first event.
  fTrackSelection = new AliTPCcalibTime("calibTime","time dependent Vdrift calibration", startTime, endTime, 20*60);
  
  fTrackSelection->UpdateEventInfo(esd); // needed for getting the run number and time stamp information correct on the offline side
    
  //TArrayI clusterSideA(esd->GetNumberOfTracks());
  //TArrayI clusterSideC(esd->GetNumberOfTracks());
  Float_t bz = AliTracker::GetBz();
  Double_t vtxx[3]={0,0,0};
  Double_t svtxx[3]={0.000001,0.000001,100.};
  AliESDVertex vtx(vtxx,svtxx);

  for(Int_t i=0; i<esd->GetNumberOfTracks(); ++i){

    AliESDtrack *track0 = esd->GetTrack(i); // track 0 upper part   
    if(!track0) continue;    
    if(!track0->GetOuterParam()) continue;
    if(track0->GetOuterParam()->GetAlpha()<0) continue;
    Double_t d1[3];
    track0->GetDirection(d1);    

    for(Int_t j=0; j<esd->GetNumberOfTracks(); ++j){
      
      if(i==j) continue;
      AliESDtrack *track1 = esd->GetTrack(j); //track 1 lower part
      if(!track1) continue; 
      if(!track1->GetOuterParam()) continue;
      if( track0->GetTPCNcls() + track1->GetTPCNcls()< 80 /*kMinClusters*/) continue;
      
      //Int_t nAC = TMath::Max( TMath::Min(clusterSideA[i], clusterSideC[j]), TMath::Min(clusterSideC[i], clusterSideA[j]));
      //if(nAC<30/*kMinClustersCross*/) continue; 
      //Int_t nA0 = clusterSideA[i];
      //Int_t nC0 = clusterSideC[i];
      //Int_t nA1 = clusterSideA[j];
      //Int_t nC1 = clusterSideC[j];
      //      if (track1->GetOuterParam()->GetAlpha()>0) continue;
      //
      Double_t d2[3];
      track1->GetDirection(d2);
      
//       AliTPCseed * seed0 = (AliTPCseed*) tpcSeeds.At(i);
//       AliTPCseed * seed1 = (AliTPCseed*) tpcSeeds.At(j);
//       if (! seed0) continue; 
//       if (! seed1) continue;
      Float_t dir = (d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]);
      Float_t dist0  = track0->GetLinearD(0,0);
      Float_t dist1  = track1->GetLinearD(0,0);
      //
      // conservative cuts - convergence to be guarantied
      // applying before track propagation
      if (TMath::Abs(TMath::Abs(dist0)-TMath::Abs(dist1))>3 /*fCutMaxD*/) continue;   // distance to the 0,0
      if (TMath::Abs(dir)<TMath::Abs(-0.99/*fCutMinDir*/)) continue;               // direction vector product
                  
      Float_t dvertex0[2];   //distance to 0,0
      Float_t dvertex1[2];   //distance to 0,0 
      track0->GetDZ(0,0,0,bz,dvertex0);
      track1->GetDZ(0,0,0,bz,dvertex1);
      if (TMath::Abs(dvertex0[1])>250) continue;
      if (TMath::Abs(dvertex1[1])>250) continue;

      Float_t dmax = TMath::Max(TMath::Abs(dist0),TMath::Abs(dist1));
      AliExternalTrackParam param0(*track0);
      AliExternalTrackParam param1(*track1);
      //
      // Propagate using Magnetic field and correct for material budget
      //    
      
      AliTracker::PropagateTrackTo(&param0,dmax+1,TDatabasePDG::Instance()->GetParticle("e-")->Mass(),3,kTRUE);
      AliTracker::PropagateTrackTo(&param1,dmax+1,TDatabasePDG::Instance()->GetParticle("e-")->Mass(),3,kTRUE);
            
      // Propagate rest to the 0,0 DCA - z should be ignored

      param0.PropagateToDCA(&vtx,bz,1000);
      param1.PropagateToDCA(&vtx,bz,1000);
      param0.GetDZ(0,0,0,bz,dvertex0);
      param1.GetDZ(0,0,0,bz,dvertex1);
      
      Double_t xyz0[3];
      Double_t xyz1[3];
      param0.GetXYZ(xyz0);
      param1.GetXYZ(xyz1);
      Bool_t isPair  = fTrackSelection->IsPair(&param0,&param1);
      Bool_t isCross = fTrackSelection->IsCross(track0, track1);
      Bool_t isSame  = fTrackSelection->IsSame(track0, track1);

      if((isSame) || (isCross && isPair)){
	if( track0->GetTPCNcls() + track1->GetTPCNcls()> 80 ){
	   numberOfCosmics++;
	}
      }
    } // end 2nd order loop        
  } // end 1st order loop
  
  bool condition = false;
  TString description;
 
  if(numberOfCosmics>0){    
     description.Form("Event contains %d cosmics", numberOfCosmics);
     condition = true;
  } 
  else {
    if(IsDataEvent()) {
      description.Form("No input blocks found");
    } else {
      description.Form("No DataEvent found");
    }
  }
   
   SetDescription(description.Data());
 
  // add a specific trigger decision object with initialized name
  // the readout list however is fixed 
  AliHLTTriggerDecision decision(
				 condition,
				 GetTriggerName(),
				 GetReadoutList(),
				 GetDescription()
				 );
  TriggerEvent(&decision, kAliHLTDataTypeTObject|kAliHLTDataOriginOut);

  return iResult;
}

int AliHLTTriggerCosmics::DoInit(int argc, const char** argv){
// see header file for class documentation

  int iResult = 0;

  // check if the -triggername argument is used
  // the name of the trigger determines the following initialization
  vector<const char*> remainingArgs;
  for (int i=0; i<argc; i++) {
    if (strcmp(argv[i], "-triggername")==0) {
      if (++i<argc) fName=argv[i];
      else {
	HLTError("invalid parameter for argument '-triggername', string expected");
	return -EINVAL;
      }
      continue;
    }
    remainingArgs.push_back(argv[i]);
  }

  // get path from triggername, use default object otherwise
  TString cdbPath;
  if (!fName.IsNull()) {
    cdbPath="HLT/ConfigHLT/";
    cdbPath+=fName;
  } else {
    cdbPath=fgkDefaultOCDBEntry;
  }

  iResult = ConfigureFromCDBObject(cdbPath);

  // -- Configure from the command line parameters if specified
  if (iResult>=0 && argc>0) iResult = ConfigureFromArgumentString(remainingArgs.size(), &(remainingArgs[0]));

  return iResult;
}

int AliHLTTriggerCosmics::DoDeinit(){
// see header file for class documentation

  if(fTrackSelection) delete fTrackSelection; fTrackSelection = NULL;
  return 0;
}

int AliHLTTriggerCosmics::Reconfigure(const char* cdbEntry, const char* /*chainId*/){
// see header file for class documentation

  // configure from the specified antry or the default one
  TString cdbPath;
  if (!cdbEntry || cdbEntry[0]==0) {
    if (!fName.IsNull()) {
      cdbPath="HLT/ConfigHLT/";
      cdbPath+=fName;
    } else {
      cdbPath=fgkDefaultOCDBEntry;
    }
  } else {
    cdbPath=cdbEntry;
  }

  return ConfigureFromCDBObject(cdbPath);
}

int AliHLTTriggerCosmics::ReadPreprocessorValues(const char* /*modules*/){
// see header file for class documentation
  return 0;
}

Int_t AliHLTTriggerCosmics::ConfigureFromCDBObject(TString cdbPath){
// see header file for class documentation

  Int_t iResult = 0;
  TString arguments;

  // -- check for "-" and replace by "_._" in the path name
  cdbPath.ReplaceAll("-",1,"_._",3);

  TObject *pCDBObject = LoadAndExtractOCDBObject(cdbPath);
  if (pCDBObject) {
//     AliHLTESDTrackCuts *pCuts = dynamic_cast<AliHLTESDTrackCuts*>(pCDBObject);
//     if (pCuts) {
//       HLTInfo("Received AliHLTESDTrackCuts configuration object : \'%s\'", pCuts->GetTitle());
//       if (fHLTESDTrackCuts)
// 	delete fHLTESDTrackCuts;
//       fHLTESDTrackCuts = pCuts;
//     }
//     else {
//       TObjString *pString = dynamic_cast<TObjString*>(pCDBObject);
//       if(pString){
// 	HLTInfo("Received configuration object string: \'%s\'", pString->GetString().Data());
// 	arguments+=pString->GetString().Data();
//       } 
//       else{
// 	HLTError("Configuration object \"%s\" has wrong type, required AliHLTESDTrackCuts or TObjString", cdbPath.Data());
// 	iResult=-EINVAL;
//       }
//     }
  } 
  else {
    HLTError("Cannot fetch object \"%s\" from CDB", cdbPath.Data());
    iResult=-ENOENT;
  }
  
  if( iResult>=0 && !arguments.IsNull() ){
    const Char_t* array = arguments.Data();
    iResult = ConfigureFromArgumentString(1, &array);
  }

  return iResult;
}

//int AliHLTTriggerCosmics::ScanConfigurationArgument(int argc, const char** argv){
// see header file for class documentation

//   if (argc<=0) return 0;
//   int i=0;
//   TString argument=argv[i];
// 
//   if (!fHLTESDTrackCuts)
//     fHLTESDTrackCuts = new AliHLTESDTrackCuts("AliHLTESDTrackCuts","No track cuts");
// 
//   // -maxpt
//   if (argument.CompareTo("-maxpt")==0) {
//     if (++i>=argc) return -EPROTO;
//     argument=argv[i];
// 
//     Float_t minPt, maxPt;
//     fHLTESDTrackCuts->GetPtRange(minPt,maxPt);
//     maxPt = argument.Atof(); 
//     fHLTESDTrackCuts->SetPtRange(minPt,maxPt);
// 
//     TString title = fHLTESDTrackCuts->GetTitle();
//     if (!title.CompareTo("No track cuts")) title = "";
//     else title += " && ";
//     title += Form("p_t < %f", maxPt);
//     fHLTESDTrackCuts->SetTitle(title);
//     return 2;
//   }    
// 
//   // -minpt
//   if (argument.CompareTo("-minpt")==0) {
//     if (++i>=argc) return -EPROTO;
//     argument=argv[i];
// 
//     Float_t minPt, maxPt;
//     fHLTESDTrackCuts->GetPtRange(minPt,maxPt);
//     minPt = argument.Atof(); 
//     fHLTESDTrackCuts->SetPtRange(minPt,maxPt);
// 
//     TString title = fHLTESDTrackCuts->GetTitle();
//     if (!title.CompareTo("No track cuts")) title = "";
//     else title += " && ";
//     title += Form("p_t > %f", minPt);
//     fHLTESDTrackCuts->SetTitle(title);
//     return 2;
//   }    
// 
//   // -mintracks
//   if (argument.CompareTo("-mintracks")==0) {
//     if (++i>=argc) return -EPROTO;
//     argument=argv[i];
//     fMinTracks=argument.Atoi();
//     return 2;
//   }    
// 
//   // -min-ldca
//   // minimum longitudinal dca to vertex
//   if (argument.CompareTo("-min-ldca")==0) {
//     if (++i>=argc) return -EPROTO;
//     argument=argv[i];
// 
//     fHLTESDTrackCuts->SetMinDCAToVertexZ(argument.Atof());
//     TString title = fHLTESDTrackCuts->GetTitle();
//     if (!title.CompareTo("No track cuts")) title = "";
//     else title += " && ";
//     title += Form("DCAz > %f", argument.Atof());
//     fHLTESDTrackCuts->SetTitle(title);
//     return 2;
//   }
//   
//   // -max-ldca
//   // maximum longitudinal dca to vertex
//   if (argument.CompareTo("-max-ldca")==0) {
//     if (++i>=argc) return -EPROTO;
//     argument=argv[i];
// 
//     fHLTESDTrackCuts->SetMaxDCAToVertexZ(argument.Atof());
//     TString title = fHLTESDTrackCuts->GetTitle();
//     if (!title.CompareTo("No track cuts")) title = "";
//     else title += " && ";
//     title += Form("DCAz < %f", argument.Atof());
//     fHLTESDTrackCuts->SetTitle(title);
//     return 2;
//   }
// 
//   // -min-tdca
//   // minimum transverse dca to vertex
//   if (argument.CompareTo("-min-tdca")==0) {
//     if (++i>=argc) return -EPROTO;
//     argument=argv[i];
// 
//     fHLTESDTrackCuts->SetMinDCAToVertexXY(argument.Atof());
//     TString title = fHLTESDTrackCuts->GetTitle();
//     if (!title.CompareTo("No track cuts")) title = "";
//     else title += " && ";
//     title += Form("DCAr > %f", argument.Atof());
//     fHLTESDTrackCuts->SetTitle(title);
//     return 2;
//   }
//   
//   // -max-tdca
//   // maximum transverse dca to vertex
//   if (argument.CompareTo("-max-tdca")==0) {
//     if (++i>=argc) return -EPROTO;
//     argument=argv[i];
// 
//     fHLTESDTrackCuts->SetMaxDCAToVertexXY(argument.Atof());
//     TString title = fHLTESDTrackCuts->GetTitle();
//     if (!title.CompareTo("No track cuts")) title = "";
//     else title += " && ";
//     title += Form("DCAr < %f", argument.Atof());
//     fHLTESDTrackCuts->SetTitle(title);
//     return 2;
//   }
// 
//   // -- deprecated
// 
//   // -dca-reference
//   // reference point for the transverse and longitudinal dca cut
//   if (argument.CompareTo("-dca-reference")==0) {
//     if (++i>=argc) return -EPROTO;
//     HLTWarning("argument -dca-reference deprecated, ESDTrackCuts only allow for DCA to vertex");
//     return 2;
//   }
// 
//   // -solenoidBz
//   if (argument.CompareTo("-solenoidBz")==0) {
//     if (++i>=argc) return -EPROTO;
//     HLTWarning("argument -solenoidBz is deprecated, magnetic field set up globally (%f)", GetBz());
//     return 2;
//   }
// 
//   // unknown argument
//   return -EINVAL;
//}
