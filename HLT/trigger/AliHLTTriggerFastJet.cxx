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
#include "AliHLTCaloClusterReader.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliESDVZERO.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "TLorentzVector.h"
#include "TRefArray.h"
#include "TString.h"
#include "TMap.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "AliHLTScalars.h"

#include "AliFJWrapper.h"


ClassImp(AliHLTTriggerFastJet)

AliHLTTriggerFastJet::AliHLTTriggerFastJet() :
AliHLTTrigger(),
  fEThreshold(0.0),
  fDetector("EMCAL"),
  fFastJetWrapper(NULL),
  fMakeStats(kFALSE),
  EsdTrackCuts(NULL),
  fOCDBEntry("HLT/ConfigHLT/EmcalJetTrigger")
{
  fFastJetWrapper = new AliFJWrapper("FastJet","FastJet");

  // check if the fMakeStats flag is set to use AliHLTScalars
  // to Add a new scalar:
  // scalars.Add(const char *name, const char *description, Double_t value) 
  
  //  if ( fMakeStats ) AliHLTScalars scalars;

}
//_____________________________________________________________

AliHLTTriggerFastJet::~AliHLTTriggerFastJet() {

  if (fFastJetWrapper)
    delete fFastJetWrapper;
  fFastJetWrapper = NULL;

  if (EsdTrackCuts)
    delete EsdTrackCuts;
  EsdTrackCuts = NULL;
  
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



AliHLTComponent* AliHLTTriggerFastJet::Spawn() {
  // see header file for class documentation
  return new AliHLTTriggerFastJet;
}


const char* AliHLTTriggerFastJet::GetTriggerName() const {
  // see header file for class documentation
  return "EmcalJetTrigger";
}

Int_t AliHLTTriggerFastJet::DoTrigger() {

  Int_t iResult = 0;
  AliHLTScalars scalars;
  
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  const TObject *obj = GetFirstInputObject( kAliHLTAllDataTypes, "AliESDEvent" );
  AliESDEvent   *esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  
  // check for MatchedTrack to avoid double counting
  


  if ( esd != NULL ) {
    esd->GetStdContent();
    
    Double_t vertex[3] = {0,0,0};
    vertex[0] = esd->GetVertex()->GetX();
    vertex[1] = esd->GetVertex()->GetY();
    vertex[2] = esd->GetVertex()->GetZ();
    cout << "Vertex " << vertex[0] << vertex[1] << vertex [2] << endl;
    //TLorentzVector gamma;
    
    // -- add VZERO
    // uncomment for usage (to avoid warning)    
    
    /*
    AliESDVZERO *v0  = esd->GetVZEROData();
    
    Float_t MtotVOA = v0->GetMTotV0A();
    Float_t MtotVOC = v0->GetMTotV0A();
    */
    Int_t nTracks = esd->GetNumberOfTracks();

    TObjArray *tracks = new TObjArray(nTracks);
    tracks->SetOwner(1);
    

    EsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    EsdTrackCuts->SetMinNClustersTPC(70);

    // add the tracks to FastJet Wrapper                                                          

    for ( Int_t i = 0; i < nTracks; i++ ) {

      const AliESDVertex *vtxTPC = esd->GetPrimaryVertexTPC();
      if (!vtxTPC) continue;
      cout << "Got primary vertex " << endl;

      AliESDtrack *esdtrack = esd->GetTrack(i);
      if (!esdtrack) continue;
      tracks->Add(esdtrack);
      cout << "Got esd track " << endl;

      AliESDtrack *tpctrack = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(esd),esdtrack->GetID());
      if (!tpctrack) continue;
      cout << "Got tpc track " << endl;

      AliExternalTrackParam exParam;
      Bool_t relate = tpctrack->RelateToVertexTPC(vtxTPC,esd->GetMagneticField(),kVeryBig,&exParam);
      if (!relate) {
	delete tpctrack;
	continue;
      }

      tpctrack->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());
      cout << "Set param to tpctrack " << endl;

      fFastJetWrapper->AddInputVector(tpctrack->Px(),tpctrack->Py(),tpctrack->Pz(),tpctrack->P());
      if ( fMakeStats ) scalars.Add("TracksPt","TPC tracks pT", tpctrack->Pt());
      cout << "Added track with P: " << tpctrack->P() << endl;
    }


//    // add the Calo Clusters to FastJet Wrapper
//    TRefArray *caloClustersArr = new TRefArray();
//    esd->GetEMCALClusters(caloClustersArr);
//    const Int_t nEMCALClusters = caloClustersArr->GetEntries();
//
//    for ( UInt_t ic = 0; ic < nEMCALClusters; ic++ )
//      {
//	AliESDCaloCluster *cluster = dynamic_cast<AliESDCaloCluster*>( caloClustersArr->At(ic) );
//        Double_t clusterEn = cluster->E();
//	Double_t subE = 0;
//
//        TArrayI *matched = cluster->GetTracksMatched();
//	for ( UInt_t jk = 0; jk < matched->GetSize(); jk++ )
//          {
//            if ( matched->At(jk) > -1 && matched->At(jk) < esd ->GetNumberOfTracks() )
//              {
//                AliESDtrack *track = dynamic_cast<AliESDtrack*>( tracks->At( (matched->At(jk)) ) );
//                if ( !track ) continue;
//                subE += track->P();
//              }
//          }
//        clusterEn -= subE;
//        if ( clusterEn < 0)
//          clusterEn = 0;
//        cluster->SetE( clusterEn );
//        if ( cluster->E() == 0 ) continue;
//
//        cluster->GetMomentum( gamma, vertex );
//        fFastJetWrapper->AddInputVector( gamma.Px(), gamma.Py(), gamma.Pz(), clusterEn );
//        if ( fMakeStats ) {
//           scalars.Add("ClusterEn","Clusters Energy",clusterEn);
//           scalars.Add("ClusterEta","Clusters Eta", gamma.Eta());
//           scalars.Add("ClusterPhi","Clusters Phi", gamma.Phi());
//        }
//        gamma.Clear();
//      }
//  }
    for ( Int_t j = 0; j< esd->GetNumberOfCaloClusters(); j++) {
      //AliESDCaloCluster *cluster = static_cast<AliESDCaloCluster*>(esd->GetCaloCluster(j));
      AliVCluster *cluster = dynamic_cast<AliVCluster*>(esd->GetCaloCluster(j));
      if ( cluster->IsEMCAL() ) {
        Int_t trackindex = cluster->GetTrackMatchedIndex();
   	if (trackindex<0) continue;
   	AliESDtrack *track = esd->GetTrack(trackindex);
   	if (!track) continue;
   	fFastJetWrapper->AddInputVector(track->Px(),track->Py(),track->Pz(),cluster->E());
   	cout << "Added cluster with E: " << cluster->E() << endl;
	if ( fMakeStats ) {
	  scalars.Add("ClusterEn" , "Clusters Energy", cluster->E());
	  scalars.Add("ClusterEta", "Clusters Eta"   , track->Eta());
	  scalars.Add("ClusterPhi", "Clusters Phi"   , track->Phi());
	}
      }
    }
  }
  
  fFastJetWrapper->SetupStrategyfromOpt("Best");
  fFastJetWrapper->SetupAlgorithmfromOpt("antikt");
  fFastJetWrapper->SetupSchemefromOpt("BIpt");
  fFastJetWrapper->SetupAreaTypefromOpt("active");
  fFastJetWrapper->SetNRepeats(1);
  fFastJetWrapper->SetGhostArea(0.01);
  fFastJetWrapper->SetMaxRap(0.9);
  fFastJetWrapper->SetR(0.4);
  
  double median = 0.0;
  double sigma  = 0.0;

  fFastJetWrapper->Run();
  fFastJetWrapper->GetMedianAndSigma(median,sigma);
  
  // checks if trigger criteria are fulfilled
  if ( TriggerOnJet(fFastJetWrapper->GetSubtractedJetsPts(median)) ) {
    if ( fMakeStats ) iResult = PushBack(&scalars, kAliHLTDataTypeEventStatistics|kAliHLTDataOriginHLT);
    return iResult;
  }
 
  // if we got to this point it meens the trigger criteria were not fulfilled
  TString description;
  description.Form(" No jets with energy >  %.02f GeV found! ", fDetector.Data(), fEThreshold);
  SetDescription(description.Data());
  TriggerEvent(kFALSE);
 
  if ( fMakeStats ) iResult = PushBack(&scalars, kAliHLTDataTypeEventStatistics|kAliHLTDataOriginHLT);
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
      
      // if fMakeStats true fill stats
      if ( fMakeStats )	scalars.Add("JetsPt","Jets pT", Jet[ij]);

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

  if ( argument.CompareTo("-makestats") == 0) {
    fMakeStats = kTRUE;
    return 2;
  }
  
  return -EINVAL;

}
//______________________________________________________________

void AliHLTTriggerFastJet::GetOutputDataTypes(AliHLTComponentDataTypeList &list) const {
  // return the output data types generated
  
  list.push_back(kAliHLTDataTypeTriggerDecision);
  list.push_back(kAliHLTDataTypeEventStatistics|kAliHLTDataOriginHLT);
}
//_______________________________________________________________

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
