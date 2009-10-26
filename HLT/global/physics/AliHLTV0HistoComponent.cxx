//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Ovrebekk <ovrebekk@ift.uib.no>                  *
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

/** @file   AliHLTV0HistoComponent.cxx
    @author Gaute Ovrebekk
    @brief  Component for ploting charge in clusters
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTV0HistoComponent.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include <TFile.h>
#include <TString.h>
#include "TObjString.h"
#include "TObjArray.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliHLTMessage.h"

//#include "AliHLTTPC.h"
//#include <stdlib.h>
//#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTV0HistoComponent)

AliHLTV0HistoComponent::AliHLTV0HistoComponent()
:
  fKShort(0),
  fLambda(0),
  fAP(0),
  fNKShorts(0),
  fNLambdas(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTV0HistoComponent::~AliHLTV0HistoComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTV0HistoComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "V0Histo";
}

void AliHLTV0HistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginOut );
}

AliHLTComponentDataType AliHLTV0HistoComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram;
}

void AliHLTV0HistoComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 80000;
  inputMultiplier = 1;
}

AliHLTComponent* AliHLTV0HistoComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTV0HistoComponent;
}

int AliHLTV0HistoComponent::DoInit( int argc, const char** argv )
{
  // init
  
  fKShort = new TH1F("hKShort","HLT KShort inv mass",80,0.4,.6); 
  fKShort->SetFillColor(kBlue);
  fKShort->SetStats(0);
  fLambda = new TH1F("hLambda","HLT Lambda inv mass",50,1.0,1.4); 
  fLambda->SetFillColor(kBlue);
  fLambda->SetStats(0);
  fAP = new TH2F("hAP","HLT Armenteros-Podolanski plot",60,-1.,1.,60,0.,0.3);
  fAP->SetMarkerStyle(8);
  fAP->SetMarkerSize(0.4);
  fAP->SetYTitle("p_{t}[GeV]");
  fAP->SetXTitle("(p^{+}_{L}-p^{-}_{L})/(p^{+}_{L}+p^{-}_{L})");
  fAP->SetStats(0);
  fAP->SetOption("COLZ");
  fNKShorts = 0;
  fNLambdas = 0;

  int iResult=0;
  TString configuration="";
  TString argument="";
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;
  }
  
  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  }  

  return iResult; 
}
  
int AliHLTV0HistoComponent::DoDeinit()
{
  // see header file for class documentation
  delete fKShort;
  delete fLambda;
  delete fAP;
  return 0;
}

int AliHLTV0HistoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {

    AliESDEvent *event = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
    event->GetStdContent();
    Int_t nV0 = event->GetNumberOfV0s();
    AliKFParticle::SetField( event->GetMagneticField() );
            
    for (Int_t iv=0; iv<nV0; iv++) {
      AliESDv0 *v0 = event->GetV0(iv);
	
      AliESDtrack *t1=event->GetTrack(v0->GetPindex());
      AliESDtrack *t2=event->GetTrack(v0->GetNindex());      

      AliKFParticle kf1( *t1->GetInnerParam(), 211 );
      AliKFParticle kf2( *t2->GetInnerParam(), 211 );

      AliKFVertex primVtx( *event->GetPrimaryVertexTracks() );
      double dev1 = kf1.GetDeviationFromVertex( primVtx );
      double dev2 = kf2.GetDeviationFromVertex( primVtx );
      
      AliKFParticle V0( kf1, kf2 );
      primVtx+=V0;
      V0.SetProductionVertex( primVtx );

      double dx = V0.GetX()-primVtx.GetX();
      double dy = V0.GetY()-primVtx.GetY();
      double r = sqrt(dx*dx + dy*dy);
      
      AliKFParticle kf01 = kf1, kf02 = kf2;
      kf01.SetProductionVertex(V0);
      kf02.SetProductionVertex(V0);
      kf01.TransportToProductionVertex();
      kf02.TransportToProductionVertex();
      
      double px1=kf01.GetPx(), py1=kf01.GetPy(), pz1=kf01.GetPz();
      double px2=kf02.GetPx(), py2=kf02.GetPy(), pz2=kf02.GetPz();
      double px = px1+px2, py = py1+py2, pz = pz1+pz2;
      double p = sqrt(px*px+py*py+pz*pz);
      double l1 = (px*px1 + py*py1 + pz*pz1)/p;
      double l2 = (px*px2 + py*py2 + pz*pz2)/p;
      double pt = sqrt(px1*px1+py1*py1+pz1*pz1 - l1*l1);
      double ap = (l1-l2)/(l1+l2);

      // AP plot
      //if( dev1>3.5 && dev2>3.5 ){
	if( fAP ) fAP->Fill( ap, pt );
	//}

      // kShort finder 
            
      if( fKShort ) fKShort->Fill( V0.GetMass() );
      fNKShorts++;

      // Lambda finder 
      if( r>=2. && fabs( ap )>= 0.4 ){                    
	fNLambdas++;
	AliKFParticle kP, kpi;
       if( ap<0 ){ 
	  kP = AliKFParticle( *t1->GetInnerParam(), 2212 );
	  kpi = AliKFParticle( *t2->GetInnerParam(), 211 );
	} else {
	  kP = AliKFParticle( *t2->GetInnerParam(), 2212 );
	  kpi = AliKFParticle( *t1->GetInnerParam(), 211 );
	}
	AliKFParticle Lambda = AliKFParticle( kP, kpi );
	Lambda.SetProductionVertex( primVtx );
	
	if( fLambda ) fLambda->Fill( Lambda.GetMass() );
      }
    }// V0's
  
    if( fKShort ){
      PushBack( (TObject*) fKShort, kAliHLTDataTypeHistogram,0);
    }
    if( fAP ){
      PushBack( (TObject*) fAP, kAliHLTDataTypeHistogram,0);
    }
    if( fLambda ){
      PushBack( (TObject*) fLambda, kAliHLTDataTypeHistogram, 0);
    }
  }  
  
  HLTInfo("V0Histo found %d K0*'s and %d Lambdas total", fNKShorts, fNLambdas );
  
  return 0;
}

int AliHLTV0HistoComponent::Configure(const char* arguments)
{
  
  int iResult=0;
  if (!arguments) return iResult;
  
  TString allArgs=arguments;
  TString argument;
  
  TObjArray* pTokens=allArgs.Tokenize(" ");
  
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
      
      if (argument.CompareTo("-plot-all")==0) {
	HLTInfo("Ploting charge of all clusters");
	//fPlotAll = kTRUE;
	continue;
      }
      
      else if (argument.CompareTo("-plot-trackclusters")==0) {
	HLTInfo("Ploting charge of clusters used on a track");
	//fPlotAll = kFALSE;
	continue;
      }
      else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  
  return iResult;
}

int AliHLTV0HistoComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path="HLT/ConfigTPC/KryptonHistoComponent";
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
	iResult=Configure(pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }

  return iResult;
}
