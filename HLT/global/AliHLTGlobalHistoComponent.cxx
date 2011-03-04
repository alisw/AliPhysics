// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/// @file   AliHLTGlobalHistoComponent.cxx
/// @author Matthias Richter
/// @date   2010-09-16
/// @brief  A histogramming component for global ESD properties based
///         on the AliHLTTTreeProcessor

#include "AliHLTGlobalHistoComponent.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "TTree.h"
#include "TString.h"
#include <cassert>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalHistoComponent)

AliHLTGlobalHistoComponent::AliHLTGlobalHistoComponent()
  : AliHLTTTreeProcessor()
  , fEvent(0)
  , fNofTracks(0)
  //, fNofV0s(0)
  //, fNofUPCpairs(0)
  , fNofContributors(0)
  , fVertexX(-99)
  , fVertexY(-99)
  , fVertexZ(-99)
  , fVertexStatus(kFALSE)
  , fMaxTrackCount(20000)
  , fTrackVariables()
  , fTrackVariablesInt()
  //, fV0Variables()
  //, fUPCVariables()
  //, fNEvents(0)
  //, fNGammas(0)
  //, fNKShorts(0)
  //, fNLambdas(0)
  //, fNPi0s(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTGlobalHistoComponent::~AliHLTGlobalHistoComponent(){
  // see header file for class documentation
  fTrackVariables.Reset();
  fTrackVariablesInt.Reset();
  //fV0Variables.Reset();
  //fUPCVariables.Reset();
}

void AliHLTGlobalHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list){
  // see header file for class documentation
  list.push_back(kAliHLTAllDataTypes);
}

AliHLTComponentDataType AliHLTGlobalHistoComponent::GetOutputDataType(){
  // see header file for class documentation
  return kAliHLTDataTypeHistogram|kAliHLTDataOriginOut;
}

TTree* AliHLTGlobalHistoComponent::CreateTree(int /*argc*/, const char** /*argv*/){
// create the tree and branches
 
  int iResult=0;
  TTree* pTree = new TTree("ESDproperties", "HLT ESD properties");
  if (!pTree) return NULL;

  const char* trackVariableNames = {
    // Note the black at the end of each name!
    "Track_pt "
    "Track_phi "
    "Track_eta "
    "Track_p "
    "Track_theta "
    "Track_Nclusters "
    "Track_status "
    "Track_charge "
    "Track_DCAr "
    "Track_DCAz "
    "Track_dEdx "
  };
  
  const char* trackIntVariableNames = {
    // Note the black at the end of each name!
    "Track_status "
  };
  
//   const char* V0VariableNames = {
//     // Note the black at the end of each name!
//     "V0_AP "
//     "V0_pt "  
//     "clust1 "
//     "clust2 "
//     "dev1 "
//     "dev2 "
//     "devPrim "
//     "length "
//     "sigmaLength "
//     "r "
//   };
//   
//   const char* UPCVariableNames = {
//     // Note the black at the end of each name!
//     "nClusters_1 "
//     "nClusters_2 "
//     "polarity_1 "
//     "polarity_2 "
//     "px_1 "
//     "py_1 "
//     "px_2 "
//     "py_2 "
//   };
    
  //int maxV0Count    = 100000;
  //int maxUPCCount   = 1;
     
  if ((iResult=fTrackVariables.Init(fMaxTrackCount, trackVariableNames))<0) {
    HLTError("failed to initialize internal structure for track properties (float)");
  }
  if ((iResult=fTrackVariablesInt.Init(fMaxTrackCount, trackIntVariableNames))<0) {
    HLTError("failed to initialize internal structure for track properties (int)");
  }
//   if ((iResult=fV0Variables.Init(maxV0Count, V0VariableNames))<0) {
//     HLTError("failed to initialize internal structure for V0 properties (float)");
//   }
//   if ((iResult=fUPCVariables.Init(maxUPCCount, UPCVariableNames))<0) {
//     HLTError("failed to initialize internal structure for UPC properties (float)");
//   }
  
  if (iResult>=0) {
    pTree->Branch("event",        &fEvent,           "event/I");
    pTree->Branch("trackcount",   &fNofTracks,       "trackcount/I");
    pTree->Branch("vertexX",      &fVertexX,         "vertexX/F");
    pTree->Branch("vertexY",      &fVertexY,         "vertexY/F");
    pTree->Branch("vertexZ",      &fVertexZ,         "vertexZ/F");
    //pTree->Branch("V0",           &fNofV0s,          "V0/I");
    //pTree->Branch("UPC",          &fNofUPCpairs,     "UPC/I");
    pTree->Branch("nContributors",&fNofContributors, "nContributors/I");
    pTree->Branch("vertexStatus", &fVertexStatus,    "vertexStatus/I");

    int i=0;
    // FIXME: this is a bit ugly since type 'f' and 'i' are specified
    // explicitely. Would be better to use a function like
    // AliHLTGlobalHistoVariables::GetType but could not get this working
    for (i=0; i<fTrackVariables.Variables(); i++) {
      TString specifier=fTrackVariables.GetKey(i);
      float* pArray=fTrackVariables.GetArray(specifier);
      specifier+="[trackcount]/f";
      pTree->Branch(fTrackVariables.GetKey(i), pArray, specifier.Data());
    }
    for (i=0; i<fTrackVariablesInt.Variables(); i++) {
      TString specifier=fTrackVariablesInt.GetKey(i);
      int* pArray=fTrackVariablesInt.GetArray(specifier);
      specifier+="[trackcount]/i";
      pTree->Branch(fTrackVariablesInt.GetKey(i), pArray, specifier.Data());
    }    
//     for (i=0; i<fV0Variables.Variables(); i++) {
//       TString specifier=fV0Variables.GetKey(i);
//       float* pArray=fV0Variables.GetArray(specifier);
//       specifier+="[V0]/f";
//       pTree->Branch(fV0Variables.GetKey(i), pArray, specifier.Data());
//     }
//     for (i=0; i<fUPCVariables.Variables(); i++) {
//       TString specifier=fUPCVariables.GetKey(i);
//       float* pArray=fUPCVariables.GetArray(specifier);
//       specifier+="[UPC]/f";
//       pTree->Branch(fUPCVariables.GetKey(i), pArray, specifier.Data());
//     }
  } else {
    delete pTree;
    pTree=NULL;
  }
  
  return pTree;
}

void AliHLTGlobalHistoComponent::FillHistogramDefinitions(){
  /// default histogram definitions
}

int AliHLTGlobalHistoComponent::FillTree(TTree* pTree, const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ){

  /// fill the tree from the ESD
  int iResult=0;
  if (!IsDataEvent()) return 0;

  ResetVariables();

  // fetch ESD from input stream
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  esd->GetStdContent();

  // fill track variables
  fNofTracks       = esd->GetNumberOfTracks();
  fVertexX         = esd->GetPrimaryVertexTracks()->GetX();
  fVertexY         = esd->GetPrimaryVertexTracks()->GetY();
  fVertexZ         = esd->GetPrimaryVertexTracks()->GetZ();
  //fNofV0s          = esd->GetNumberOfV0s();
  fNofContributors = esd->GetPrimaryVertexTracks()->GetNContributors();
  fVertexStatus    = esd->GetPrimaryVertexTracks()->GetStatus();
  //fNofUPCpairs     = 1;
  
  for (int i=0; i<fNofTracks; i++) {
    AliESDtrack *esdTrack = esd->GetTrack(i);
    if (!esdTrack) continue;
    
    Float_t DCAr, DCAz = -99;
    esdTrack->GetImpactParametersTPC(DCAr, DCAz);
    
    fTrackVariables.Fill("Track_pt"        , esdTrack->Pt()                      );
    fTrackVariables.Fill("Track_phi"       , esdTrack->Phi()*TMath::RadToDeg()   );
    fTrackVariables.Fill("Track_eta"       , esdTrack->Theta()                   );
    fTrackVariables.Fill("Track_p"         , esdTrack->P()                       );
    fTrackVariables.Fill("Track_theta"     , esdTrack->Theta()*TMath::RadToDeg() );
    fTrackVariables.Fill("Track_Nclusters" , esdTrack->GetTPCNcls()              );
    fTrackVariables.Fill("Track_status"    , esdTrack->GetStatus()               );
    fTrackVariables.Fill("Track_charge"    , esdTrack->Charge()                  );
    fTrackVariables.Fill("Track_DCAr"      , DCAr                        	 );
    fTrackVariables.Fill("Track_DCAz"      , DCAz                                );   
    fTrackVariables.Fill("Track_dEdx"      , esdTrack->GetTPCsignal()            );   
    fTrackVariablesInt.Fill("Track_status" , esdTrack->GetStatus()               );      
  }
  //HLTInfo("added parameters for %d tracks", fNofTracks);
  
//   if(fNofTracks==2){
//      AliESDtrack *esdTrack1 = esd->GetTrack(0);
//      if(!esdTrack1) return 0;
//      AliESDtrack *esdTrack2 = esd->GetTrack(1); 
//      if(!esdTrack2) return 0;
//      
//      if(esdTrack1->Charge()*esdTrack2->Charge()<0){
//      
//        fUPCVariables.Fill("nClusters_1", esdTrack1->GetTPCNcls() );
//        fUPCVariables.Fill("nClusters_2", esdTrack2->GetTPCNcls() );
//        fUPCVariables.Fill("polarity_1",  esdTrack1->Charge()     );
//        fUPCVariables.Fill("polarity_2",  esdTrack2->Charge()     );
//        fUPCVariables.Fill("px_1",        esdTrack1->Px()         );
//        fUPCVariables.Fill("py_1",        esdTrack1->Py()         );
//        fUPCVariables.Fill("px_2",        esdTrack2->Px()         );
//        fUPCVariables.Fill("py_2",        esdTrack2->Py()         );
//     } 
//   }
  
//  AliKFParticle::SetField( esd->GetMagneticField() );

  //const double kKsMass = 0.49767;
  //const double kLambdaMass = 1.11568;
  //const double kPi0Mass = 0.13498;

  //std::vector<AliKFParticle> vGammas;
  
  
//   for (int i=0; i<fNofV0s; i++) {
//     AliESDv0 *esdV0 = esd->GetV0(i);
//     if (!esdV0) continue;
//     
//     AliESDtrack *t1 = esd->GetTrack( esd->GetV0(i)->GetNindex());
//     AliESDtrack *t2 = esd->GetTrack( esd->GetV0(i)->GetPindex());      
// 
//     AliKFParticle kf1( *t1, 11 );
//     AliKFParticle kf2( *t2, 11 );
// 
//     AliKFVertex primVtx( *esd->GetPrimaryVertexTracks() );
//     double dev1 = kf1.GetDeviationFromVertex( primVtx );
//     double dev2 = kf2.GetDeviationFromVertex( primVtx );
//     
//     AliKFParticle v0( kf1, kf2 );
//     double devPrim = v0.GetDeviationFromVertex( primVtx );
//     primVtx+=v0;
//     v0.SetProductionVertex( primVtx );
// 
//     Double_t length, sigmaLength;
//     if( v0.GetDecayLength( length, sigmaLength ) ) continue;
// 
//     double dx = v0.GetX()-primVtx.GetX();
//     double dy = v0.GetY()-primVtx.GetY();
//     double r = sqrt(dx*dx + dy*dy);
//     
// 
//     // Armenteros-Podolanski plot
// 
//     double pt=0, ap=0;
//     //{
//       AliKFParticle kf01 = kf1, kf02 = kf2;
//       kf01.SetProductionVertex(v0);
//       kf02.SetProductionVertex(v0);
//       kf01.TransportToProductionVertex();
//       kf02.TransportToProductionVertex();      
//       double px1=kf01.GetPx(), py1=kf01.GetPy(), pz1=kf01.GetPz();
//       double px2=kf02.GetPx(), py2=kf02.GetPy(), pz2=kf02.GetPz();
//       double px = px1+px2, py = py1+py2, pz = pz1+pz2;
//       double p = sqrt(px*px+py*py+pz*pz);
//       double l1 = (px*px1 + py*py1 + pz*pz1)/p;
//       double l2 = (px*px2 + py*py2 + pz*pz2)/p;
//       pt = sqrt(px1*px1+py1*py1+pz1*pz1 - l1*l1);
//       ap = (l2-l1)/(l1+l2);
//     //}
//     
// //     if( 
// //        t1->GetTPCNcls()>=fAPCuts[0]
// //        && t2->GetTPCNcls()>=fAPCuts[0]
// //        && dev1>=fAPCuts[1]
// //        && dev2>=fAPCuts[1]
// //        && devPrim <= fAPCuts[2]
// //        && length >= fAPCuts[3]*sigmaLength
// //        && length >= fAPCuts[4]
// //        && r <= fAPCuts[5]
// //        )
// //        //{     
// //        //if( fAP ) fAP->Fill( ap, pt );
// //        //} 
// // 
// //     // Gamma finder
// // 
// //     bool isGamma = 0;
// //     
// //     if( 
// //        t1->GetTPCNcls()>=fGammaCuts[0]
// //        && t2->GetTPCNcls()>=fGammaCuts[0]
// //        && dev1>=fGammaCuts[1]
// //        && dev2>=fGammaCuts[1]
// //        && devPrim <= fGammaCuts[2]
// //        && length >= fGammaCuts[3]*sigmaLength
// //        && length >= fGammaCuts[4]
// //        && r <= fGammaCuts[5]
// //        ){
// //       double mass, error;	
// //       v0.GetMass(mass,error); 
// //       //if( fGamma ) fGamma->Fill( mass );
// // 
// //       if( TMath::Abs(mass)<=fGammaCuts[6]*error || TMath::Abs(mass)<=fGammaCuts[7] ){	
// //         AliKFParticle gamma = v0;
// //         gamma.SetMassConstraint(0);
// // //      if( fGammaXY
// // //          &&  t1->GetTPCNcls()>=60
// // //          && t2->GetTPCNcls()>=60
// // //          ) fGammaXY->Fill(gamma.GetX(), gamma.GetY());
// //         isGamma = 1;
// //         fNGammas++;
// //         vGammas.push_back( gamma );
// //       } 	   
// //     }
// //     
// //     if( isGamma ) continue;
// // 
// // 
// //     // KShort finder
// //     
// //     bool isKs = 0;
// //     
// //     if( 
// //        t1->GetTPCNcls()>=fKsCuts[0]
// //        && t2->GetTPCNcls()>=fKsCuts[0]
// //        && dev1>=fKsCuts[1]
// //        && dev2>=fKsCuts[1]
// //        && devPrim <= fKsCuts[2]
// //        && length >= fKsCuts[3]*sigmaLength
// //        && length >= fKsCuts[4]
// //        && r <= fKsCuts[5]
// //        ){     
// //     
// //       AliKFParticle piN( *t1, 211 );  
// //       AliKFParticle piP( *t2, 211 );  
// // 
// //       AliKFParticle Ks( piN, piP );
// //       Ks.SetProductionVertex( primVtx );
// // 
// //       double mass, error;
// //       Ks.GetMass( mass, error);
// //       //if( fKShort ) fKShort->Fill( mass );  
// //       if( TMath::Abs( mass - kKsMass )<=fKsCuts[6]*error || TMath::Abs( mass - kKsMass )<=fKsCuts[7] ){  
// //         isKs = 1;
// //         fNKShorts++;
// //       }
// //     }
// //     
// //     if( isKs ) continue;
// //     
// //     // Lambda finder 
// //     //printf("QQQQQQQQQQQQQQQQq :%f\n",fLambdaCuts[0]);
// //     if( 
// //        t1->GetTPCNcls()>=fLambdaCuts[0]
// //        && t2->GetTPCNcls()>=fLambdaCuts[0]
// //        && dev1>=fLambdaCuts[1]
// //        && dev2>=fLambdaCuts[1]
// //        && devPrim <= fLambdaCuts[2]
// //        && length >= fLambdaCuts[3]*sigmaLength
// //        && length >= fLambdaCuts[4]
// //        && r <= fLambdaCuts[5]
// //        && TMath::Abs( ap )>.4
// //        ){
// // 
// //       AliKFParticle kP, kpi;
// //       if( ap<0 ){ 
// //         kP = AliKFParticle( *t2, 2212 );
// //         kpi = AliKFParticle( *t1, 211 );
// //       } else {
// //         kP = AliKFParticle( *t1, 2212 );
// //         kpi = AliKFParticle( *t2, 211 );
// //       }
// // 
// //       AliKFParticle lambda = AliKFParticle( kP, kpi );
// //       lambda.SetProductionVertex( primVtx );  
// //       //double mass, error;
// //       lambda.GetMass( Lmass, Lerror);
// //       //if( fLambda ) fLambda->Fill( mass );
// //       if( TMath::Abs( Lmass - kLambdaMass )<=fLambdaCuts[6]*Lerror || TMath::Abs( Lmass - kLambdaMass )<=fLambdaCuts[7] ){
// //         fNLambdas++;
// //       }
// //    }
// 
//         
//     fV0Variables.Fill("V0_AP", ap);
//     fV0Variables.Fill("V0_pt", pt); 
//     fV0Variables.Fill("clust1", t1->GetTPCNcls()); 
//     fV0Variables.Fill("clust2", t2->GetTPCNcls()); 
//     fV0Variables.Fill("dev1", dev1); 
//     fV0Variables.Fill("dev2", dev2); 
//     fV0Variables.Fill("devPrim", devPrim); 
//     fV0Variables.Fill("length", length); 
//     fV0Variables.Fill("sigmaLength", sigmaLength); 
//     fV0Variables.Fill("r", r); 
//       
//  } // end of loop over V0s
  
  if (iResult<0) {
    // fill an empty event
    ResetVariables();
  }
  fEvent++;

  pTree->Fill();
  return iResult;
}

int AliHLTGlobalHistoComponent::ResetVariables(){
/// reset all filling variables
  fNofTracks=0;
  //fNofV0s=0;
  fTrackVariables.ResetCount();
  fTrackVariablesInt.ResetCount();
  //fV0Variables.ResetCount();
  //fUPCVariables.ResetCount();
  return 0;
}

AliHLTComponentDataType AliHLTGlobalHistoComponent::GetOriginDataType() const{
// get the origin of the output data
  return kAliHLTDataTypeHistogram|kAliHLTDataOriginHLT;
}

int AliHLTGlobalHistoComponent::ScanConfigurationArgument(int argc, const char** argv){
/// inherited from AliHLTComponent, scan argument

  if (argv==NULL || argc<1) return 0;

  int i=0;
  TString argument=argv[i];

  // -argument1 with one parameter
  if(argument.CompareTo("-max-track-count")==0){
    
    if (++i>=argc) return -EPROTO; // missing parameter
    argument = argv[i];
    fMaxTrackCount = argument.Atoi();
    
    HLTInfo("got %s with parameter %s", argument.Data(), argv[i]);
    return ++i; // two arguments scanned
  }

//   // -argument1 with one parameter
//   if (argument.CompareTo("-argument1")==0) {
//     if (++i>=argc) return -EPROTO; // missing parameter
//     HLTInfo("got %s with parameter %s", argument.Data(), argv[i]);
//     return ++i; // two arguments scanned
//   }
// 
//   // -argument2 without parameter
//   if (argument.CompareTo("-argument2")==0) {
//     HLTInfo("got %s", argument.Data());
//     return ++i; // one argument scanned
//   }
// 
  // no recognized argument, forward to base class
  return AliHLTTTreeProcessor::ScanConfigurationArgument(argc, argv);
}

int AliHLTGlobalHistoComponent::Reconfigure(const char* cdbEntry, const char* /*chainId*/){  
// see header file for class documentation

  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigHLT/";
    cdbPath+=GetComponentID();
  }

  return ConfigureFromCDBTObjString(cdbPath.Data());
}


