// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kenneth Aamodt <kenneth.aamodt@cern.ch>               *
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

/// @file   AliHLTTriggerGammaConversion.cxx
/// @author Kenneth Aamodt
/// @date   2009-11-01
/// @brief  HLT trigger component for gamma conversions.
///         

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTriggerGammaConversion.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "TObjArray.h"
#include "TObjString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerGammaConversion)

AliHLTTriggerGammaConversion::AliHLTTriggerGammaConversion() 
  : AliHLTTrigger()
  , fMaxInvMass(0.05)
  , fPtMax(0.0)
  , fPtMin(0.0)
  , fMaxDca(10.0)
  , fMaxR(200)
  , fNReconstructedGammas(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

const char* AliHLTTriggerGammaConversion::fgkOCDBEntry="HLT/ConfigHLT/GammaConversionTrigger";

AliHLTTriggerGammaConversion::~AliHLTTriggerGammaConversion()
{
  // see header file for class documentation
}

const char* AliHLTTriggerGammaConversion::GetTriggerName() const
{
  // see header file for class documentation
  return "GammaConversionTrigger";
}

AliHLTComponent* AliHLTTriggerGammaConversion::Spawn()
{
  // see header file for class documentation
  return new AliHLTTriggerGammaConversion;
}

int AliHLTTriggerGammaConversion::DoTrigger()
{
  // see header file for class documentation
  int iResult=0;
  fNReconstructedGammas=0;

  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {

    AliESDEvent *event = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
    event->GetStdContent();
    Int_t nV0 = event->GetNumberOfV0s();
    if(nV0<=0){
      continue;
    }
    AliKFParticle::SetField( event->GetMagneticField() );

    for (Int_t iv=0; iv<nV0; iv++) {
	
      AliESDtrack *t1=event->GetTrack( event->GetV0(iv)->GetNindex());
      AliESDtrack *t2=event->GetTrack( event->GetV0(iv)->GetPindex());      

      AliKFParticle kf1( *t1->GetInnerParam(), 11 );
      AliKFParticle kf2( *t2->GetInnerParam(), 11 );

      AliKFVertex primVtx( *event->GetPrimaryVertexTracks() );

      AliKFParticle v0( kf1, kf2 );
      primVtx+=v0;
      v0.SetProductionVertex( primVtx );

      if(kf1.GetDistanceFromParticle(kf2)>fMaxDca){
	continue;
      }

      double mass, error;       
      v0.GetMass(mass,error);      
      if( TMath::Abs(mass)>fMaxInvMass ){
	continue;
      }

      AliKFParticle gamma = v0;
      gamma.SetMassConstraint(0);
      
      double r= sqrt(gamma.GetX()*gamma.GetX()+gamma.GetY()*gamma.GetY());
      if(r>fMaxR){
	continue;
      }
      if(gamma.GetPt()<fPtMin){
	continue;
      }
      
      if (fPtMax>0.){
	if(gamma.GetPt()>fPtMax){
	  continue;
	}
      }
      
      fNReconstructedGammas++;
    }  
  } 
  TString description;
  TString  maxInvMass, ptcut, maxDca, maxR;
  maxInvMass.Form(" mass < %.03f GeV ,",fMaxInvMass);
  if (fPtMax>fPtMin) {
    ptcut.Form(" %.02f GeV/c <= pt <= %.02f GeV/c ,", fPtMin, fPtMax);
  } else {
    ptcut.Form(" pt >= %.02f GeV/c ,", fPtMin);
  }
  maxDca.Form(" dca <= %.04fcm ,",fMaxDca);
  maxR.Form(" r <= %.02fcm", fMaxR);
  
  if(fNReconstructedGammas>0){
    description.Form("Event contains %d gamma conversions,", fNReconstructedGammas);
    description += ptcut;
    description += maxDca;
    description += maxR;
    
    SetDescription(description.Data());
    
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
    TriggerEvent(true);
    return 0;
  }
  else{
    description.Form("No Gamma Conversions reconstructed that satisfy:");
    description += ptcut;
    description += maxDca;
    description += maxR;
  }    
  SetDescription(description.Data());
  TriggerEvent(false);
  return iResult;
}

int AliHLTTriggerGammaConversion::DoInit(int argc, const char** argv)
{
  // see header file for class documentation

  // first configure the default
  int iResult=0;
  iResult=ConfigureFromCDBTObjString(fgkOCDBEntry);

  // configure from the command line parameters if specified
  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);
  return iResult;
}

int AliHLTTriggerGammaConversion::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTriggerGammaConversion::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{
  // see header file for class documentation

  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if (!entry || entry[0]==0) {
    entry=fgkOCDBEntry;
  }

  return ConfigureFromCDBTObjString(entry);
}

int AliHLTTriggerGammaConversion::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -max-invmass
  if (argument.CompareTo("-max-invmass")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fMaxInvMass=argument.Atof();
    return 2;
  }    

  // -max-pt
  if (argument.CompareTo("-max-pt")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fPtMax=argument.Atof();
    return 2;
  }    

  // -min-pt
  if (argument.CompareTo("-min-pt")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fPtMin=argument.Atof();
    return 2;
  }    

  // -max-dca
  // maximum dca of v0
  if (argument.CompareTo("-max-dca")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fMaxDca=argument.Atof();
    return 2;
  }

  // -max-r
  // maximum radius of v0 in xy-plane
  if (argument.CompareTo("-max-r")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fMaxR=argument.Atof();
    return 2;
  }

  // unknown argument
  return -EINVAL;
}
