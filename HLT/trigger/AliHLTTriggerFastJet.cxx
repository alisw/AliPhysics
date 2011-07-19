//**************************************************************************                        
//* This file is property of and copyright by the ALICE HLT Project        *                        
//* ALICE Experiment at CERN, All rights reserved.                         *                        
//*                                                                        *                        
//* Primary Authors: leonidas.xaplanteris@gmail.com                        *                        
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

#include "AliHLTTriggerFastJet.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDCaloCluster.h"
#include "AliESDVZERO.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "TLorentzVector.h"
#include "TRefArray.h"
#include "TString.h"
#include "TMap.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/AreaDefinition.hh>

#include "AliFJWrapper.h"


ClassImp(AliHLTTriggerFastJet)

AliHLTTriggerFastJet::AliHLTTriggerFastJet(TString detector) :
AliHLTTrigger(),
  fEThreshold(0.0),
  fDetector(detector),
  fFastJetWrapper(NULL),
  fClustersRefs(NULL),
  fOCDBEntry(""),
  fOffset(10000),
  fInputDataType()
{
  fFastJetWrapper = new AliFJWrapper("FastJet","FastJet");

  fClustersRefs = new TRefArray();

}
//_____________________________________________________________

AliHLTTriggerFastJet::~AliHLTTriggerFastJet() {

  if (fFastJetWrapper)
    delete fFastJetWrapper;
  fFastJetWrapper = NULL;

  if (fClustersRefs)
    delete fClustersRefs;
  fClustersRefs = NULL;

}
//_____________________________________________________________

int AliHLTTriggerFastJet::DoInit(int argc, const char** argv) {
  
  int iResult = ConfigureFromCDBTObjString(fOCDBEntry);
  
  if ( iResult>=0 && argc>0 ) {
    iResult = ConfigureFromArgumentString(argc, argv);
    HLTImportant("Trigger threshold set from argument string :  %.02f GeV:", fEThreshold );
  } else if ( iResult>=0 ) {
    HLTImportant("Trigger threshold set from OCDB database entry : %.02f Gev:", fEThreshold );
  }
  return iResult;
}
//______________________________________________________________

int AliHLTTriggerFastJet::DoDeInit() {
  
  return 0;

}
//______________________________________________________________

Int_t AliHLTTriggerFastJet::DoTrigger() {

  Int_t iResult = 0;

  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  const TObject *obj = GetFirstInputObject( kAliHLTAllDataTypes, "AliESDEvent" );
  AliESDEvent   *esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  
  // check for MatchedTrack to avoid double counting
  

  if ( esd != NULL ) {
    esd->GetStdContent();
    
    Double_t vertex[3] = {0,0,0};
    esd->GetVertex()->GetXYZ(vertex);
    TLorentzVector gamma;
    
    // -- add VZERO
    AliESDVZERO *v0  = esd->GetVZEROData();
    Float_t MtotVOA = v0->GetMTotV0A();
    Float_t MtotVOC = v0->GetMTotV0A();

    // add the tracks to FastJet Wrapper
    for ( Int_t i = 0; i<esd->GetNumberOfTracks(); i++ ) {
      AliESDtrack *esdtrack = esd->GetTrack(i);
      fFastJetWrapper->AddInputVector(esdtrack->Px(),esdtrack->Py(),esdtrack->Pz(),esdtrack->P());
      AliESDtrack *tpctrack = AliESDtrackCuts::GetTPCOnlyTrack(esd,esdtrack->GetID());
      if ( tpctrack )
	fFastJetWrapper->AddInputVector(tpctrack->Px(),tpctrack->Py(),tpctrack->Pz(),tpctrack->P());
      delete tpctrack;
    }
    // add the Calo Clusters to FastJet Wrapper
    for ( Int_t j = 0; j<GetClustersFromEsd(esd,fClustersRefs); j++) {
      AliESDCaloCluster *cluster = static_cast<AliESDCaloCluster*>(fClustersRefs->At(j));
      if ( cluster->IsEMCAL() ) {
	cluster->GetMomentum(gamma,vertex);
	fFastJetWrapper->AddInputVector(gamma.Px(),gamma.Py(),gamma.Pz(),gamma.P(),fOffset);
	fOffset++;
	gamma.Clear();
      }
    }
    
    
  }
  
  fFastJetWrapper->SetStrategy(fastjet::Best);
  fFastJetWrapper->SetAlgorithm(fastjet::antikt_algorithm);
  fFastJetWrapper->SetRecombScheme(fastjet::BIpt_scheme);
  fFastJetWrapper->SetAreaType(fastjet::active_area);
  fFastJetWrapper->SetNRepeats(1);
  fFastJetWrapper->SetGhostArea(0.01);
  fFastJetWrapper->SetMaxRap(0.9);
  fFastJetWrapper->SetR(0.4);
  
  double median = 0.0;
  double sigma  = 0.0;

  fFastJetWrapper->Run();
  fFastJetWrapper->GetMedianAndSigma(median,sigma);

  if ( TriggerOnJet(fFastJetWrapper->GetSubtractedJetsPts(median)) )
       return iResult;
  
  TString description;
  description.Form(" No %s jets with energy >  %.02f GeV found! ", fDetector.Data(), fEThreshold);
  SetDescription(description.Data());
  TriggerEvent(kFALSE);
  
  return iResult;

}
//______________________________________________________________

template <class T>
Bool_t AliHLTTriggerFastJet::TriggerOnJet(T Jet) {

  for ( unsigned int ij = 0; ij<Jet.size(); ij++ ) {
    if ( Jet[ij] > fEThreshold ) {
      TString description;
      description.Form(" Event contains at least one %s jet with energy > %.02f GeV! ", fDetector.Data(), fEThreshold);
      SetDescription(description.Data());
      
      // Enable the detectors for readout
      GetReadoutList().Enable(AliHLTReadoutList::kEMCAL |
			      AliHLTReadoutList::kTPC);
      
      // Add the available HLT info for readout
      GetTriggerDomain().Add(kAliHLTAnyDataTypeID, fDetector.Data());
      
      // Set trigger desicion
      TriggerEvent(kTRUE);
      
      return(kTRUE);
    }
  }
  
  return kFALSE;

}
//______________________________________________________________

int AliHLTTriggerFastJet::Reconfigure(const char* cdbEntry, const char* /*chainId*/) {
  
  const char* entry = cdbEntry;
  if ( !entry || entry[0]==0 ) entry=fOCDBEntry;

  return ConfigureFromCDBTObjString(entry);

}

//______________________________________________________________

int AliHLTTriggerFastJet::ScanConfigurationArgument(int argc, const char** argv) {

  if ( argc<=0 ) return 0;
  int i = 0;
  TString argument = argv[i];

  if ( argument.CompareTo("-energy") == 0 ) {
    if (++i>=argc) return -EPROTO;
    argument = argv[i];
    fEThreshold=argument.Atof();
    return 2;
  }
  
  return -EINVAL;

}
//______________________________________________________________

void AliHLTTriggerFastJet::GetOutputDataSize(unsigned long &constBase, double &inputMultiplier) {
  
  constBase = sizeof(AliHLTTriggerDecision) + sizeof(AliHLTDomainEntry)*14;
  inputMultiplier = 1;

}
//______________________________________________________________

void AliHLTTriggerFastJet::GetOCDBObjectDescription(TMap* const targetMap) {
  
  if ( !targetMap ) return;
  targetMap->Add(new TObjString(fOCDBEntry),
		 new TObjString(Form("%s threshold trigger OCDB object", fDetector.Data()) ) 
		 );

}
//_______________________________________________________________
