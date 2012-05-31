// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: S.Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de>    *
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
    @author Sergey Gorbunov
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
#include "TSystem.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliHLTMessage.h"
#include "TTimeStamp.h"

//#include "AliHLTTPC.h"
//#include <stdlib.h>
//#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTV0HistoComponent)

AliHLTV0HistoComponent::AliHLTV0HistoComponent() :
  fUID(0),
  fGamma(0),
  fKShort(0),
  fLambda(0),
  fPi0(0),
  fAP(0),
  fGammaXY(0),
  fNEvents(0),
  fNGammas(0),
  fNKShorts(0),
  fNLambdas(0),
  fNPi0s(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  for( int i=0; i<8; i++){
    fGammaCuts[i] = 0;
    fKsCuts[i] = 0;
    fLambdaCuts[i] = 0;
    fAPCuts[i] = 0;
  }
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
  return kAliHLTDataTypeHistogram  | kAliHLTDataOriginOut;
}

void AliHLTV0HistoComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 80000;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTV0HistoComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTV0HistoComponent;
}

int AliHLTV0HistoComponent::DoInit( int argc, const char** argv ) {
  // init

  fUID = 0;

  fGamma = new TH1F("hGamma","HLT:  #gamma inv mass",50,-.06,.2); 
  fGamma->SetFillColor(kGreen);
  fGamma->SetStats(0);
  fKShort = new TH1F("hKShort","HLT:  K_{s}^{0} inv mass",80,0.4,.6); 
  fKShort->SetFillColor(kGreen);
  fKShort->SetStats(0);
  fLambda = new TH1F("hLambda","HLT:  #Lambda^{0} inv mass",50,1.0,1.36); 
  fLambda->SetFillColor(kGreen);
  fLambda->SetStats(0);

  fPi0 = new TH1F("hPi0","HLT:  #Pi^{0} inv mass",50,0.0,0.38); 
  fPi0->SetFillColor(kGreen);
  fPi0->SetStats(0);
  
  fAP = new TH2F("hAP","HLT:  Armenteros-Podolanski plot",60,-1.,1.,60,-0.02,0.3);
  fAP->SetMarkerStyle(8);
  fAP->SetMarkerSize(0.4);
  fAP->SetYTitle("p_{t}[GeV]");
  fAP->SetXTitle("(p^{+}_{L}-p^{-}_{L})/(p^{+}_{L}+p^{-}_{L})");
  fAP->SetStats(0);
  fAP->SetOption("CONT4 Z");

  fGammaXY = new TH2F("hGammaXY","HLT:  #gamma conversions",100,-100.,100.,100,-100.,100.);
  fGammaXY->SetMarkerStyle(8);
  fGammaXY->SetMarkerSize(0.4);
  fGammaXY->SetYTitle("X[cm]");
  fGammaXY->SetXTitle("Y[cm]");
  fGammaXY->SetStats(0);
  fGammaXY->SetOption("COLZ");

  fNEvents =0;
  fNGammas = 0;
  fNKShorts = 0;
  fNLambdas = 0;
  fNPi0s = 0;

  // cuts: 
  // [0] == 0   --- N clusters on each daughter track
  // [1] == 2.5 --- (daughter-primVtx)/sigma >= cut
  // [2] == 3.5 --- (v0 - primVtx)/sigma <= cut
  // [3] == 3.0 --- (decay_length)/sigma >= cut
  // [4] == 0.0 --- (decay_length)[cm]   >= cut
  // [5] == 300.0 --- (v0 radius)[cm]    <= cut
  // [6] == 3.5  --- (v0 mass - true value)/sigma <= cut (for identification)
  // [7] == 0.05 --- (v0 mass - true value)       <= cut (for identification)

  fGammaCuts[0] = 0;
  fGammaCuts[1] = 2.5;
  fGammaCuts[2] = 3.5;
  fGammaCuts[3] = 3.0;
  fGammaCuts[4] = 0.0;
  fGammaCuts[5] = 300.0;
  fGammaCuts[6] = 3.5;
  fGammaCuts[7] = 0.05;

  fAPCuts[0] = 60;
  fAPCuts[1] = 2.5;
  fAPCuts[2] = 3.5;
  fAPCuts[3] = 3.0;
  fAPCuts[4] = 0.0;
  fAPCuts[5] = 50.0;
  fAPCuts[6] = 4.0;
  fAPCuts[7] = 0.05;

  fKsCuts[0] = 60;
  fKsCuts[1] = 2.5;
  fKsCuts[2] = 3.5;
  fKsCuts[3] = 3.0;
  fKsCuts[4] = 1.5;
  fKsCuts[5] = 50.0;
  fKsCuts[6] = 4.0;
  fKsCuts[7] = 0.015;

  fLambdaCuts[0] = 60;   // [0] == 60  --- N clusters on each daughter track
  fLambdaCuts[1] = 3.0;  // [1] == 3.0 --- (daughter-primVtx)/sigma >= cut
  fLambdaCuts[2] = 3.0;  // [2] == 3.0 --- (v0 - primVtx)/sigma <= cut
  fLambdaCuts[3] = 3.5;  // [3] == 3.5 --- (decay_length)/sigma >= cut
  fLambdaCuts[4] = 4.0;  // [4] == 0.0 --- (decay_length)[cm]   >= cut
  fLambdaCuts[5] = 50.0; // [5] == 300.0 --- (v0 radius)[cm]    <= cut
  fLambdaCuts[6] = 4.0;  // [6] == 3.5  --- (v0 mass - true value)/sigma <= cut (for identification)
  fLambdaCuts[7] = 0.03; // [7] == 0.05 --- (v0 mass - true value)       <= cut (for identification)

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
  delete fGamma;
  delete fKShort;
  delete fLambda;
  delete fAP;
  delete fGammaXY;
  fUID = 0;
  return 0;
}

int AliHLTV0HistoComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/)
{
  
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  if( fUID == 0 ){
    TTimeStamp t;
    fUID = ( gSystem->GetPid() + t.GetNanoSec())*10 + evtData.fEventID;
  }

  fNEvents++;

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {

    AliESDEvent *event = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
    event->GetStdContent();
    Int_t nV0 = event->GetNumberOfV0s();
    AliKFParticle::SetField( event->GetMagneticField() );

    const double kKsMass = 0.49767;
    const double kLambdaMass = 1.11568;
    const double kPi0Mass = 0.13498;

    std::vector<AliKFParticle> vGammas;

    for (Int_t iv=0; iv<nV0; iv++) {
	
      AliESDtrack *t1=event->GetTrack( event->GetV0(iv)->GetNindex());
      AliESDtrack *t2=event->GetTrack( event->GetV0(iv)->GetPindex());      

      AliKFParticle kf1( *t1, 11 );
      AliKFParticle kf2( *t2, 11 );

      AliKFVertex primVtx( *event->GetPrimaryVertexTracks() );
      double dev1 = kf1.GetDeviationFromVertex( primVtx );
      double dev2 = kf2.GetDeviationFromVertex( primVtx );
      
      AliKFParticle v0( kf1, kf2 );
      double devPrim = v0.GetDeviationFromVertex( primVtx );
      primVtx+=v0;
      v0.SetProductionVertex( primVtx );

      Double_t length, sigmaLength;
      if( v0.GetDecayLength( length, sigmaLength ) ) continue;

      double dx = v0.GetX()-primVtx.GetX();
      double dy = v0.GetY()-primVtx.GetY();
      double r = sqrt(dx*dx + dy*dy);

      // AP plot

      double pt=0, ap=0;
      {
	AliKFParticle kf01 = kf1, kf02 = kf2;
	kf01.SetProductionVertex(v0);
	kf02.SetProductionVertex(v0);
	kf01.TransportToProductionVertex();
	kf02.TransportToProductionVertex();      
	double px1=kf01.GetPx(), py1=kf01.GetPy(), pz1=kf01.GetPz();
	double px2=kf02.GetPx(), py2=kf02.GetPy(), pz2=kf02.GetPz();
	double px = px1+px2, py = py1+py2, pz = pz1+pz2;
	double p = sqrt(px*px+py*py+pz*pz);
	double l1 = (px*px1 + py*py1 + pz*pz1)/p;
	double l2 = (px*px2 + py*py2 + pz*pz2)/p;
	pt = sqrt(px1*px1+py1*py1+pz1*pz1 - l1*l1);
	ap = (l2-l1)/(l1+l2);
      }
      
      if( 
	 t1->GetTPCNcls()>=fAPCuts[0]
	 && t2->GetTPCNcls()>=fAPCuts[0]
	 && dev1>=fAPCuts[1]
	 && dev2>=fAPCuts[1]
	 && devPrim <= fAPCuts[2]
	 && length >= fAPCuts[3]*sigmaLength
	 && length >= fAPCuts[4]
	 && r <= fAPCuts[5]
	 ){	
	if( fAP ) fAP->Fill( ap, pt );
      } 

      // Gamma finder

      bool isGamma = 0;
      
      if( 
	 t1->GetTPCNcls()>=fGammaCuts[0]
	 && t2->GetTPCNcls()>=fGammaCuts[0]
	 && dev1>=fGammaCuts[1]
	 && dev2>=fGammaCuts[1]
	 && devPrim <= fGammaCuts[2]
	 && length >= fGammaCuts[3]*sigmaLength
	 && length >= fGammaCuts[4]
	 && r <= fGammaCuts[5]
	 ){
	double mass, error;       
	v0.GetMass(mass,error);	
	if( fGamma ) fGamma->Fill( mass );

	if( TMath::Abs(mass)<=fGammaCuts[6]*error || TMath::Abs(mass)<=fGammaCuts[7] ){	  
	  AliKFParticle gamma = v0;
	  gamma.SetMassConstraint(0);
	  if( fGammaXY
	      &&  t1->GetTPCNcls()>=60
	      && t2->GetTPCNcls()>=60
	      ) fGammaXY->Fill(gamma.GetX(), gamma.GetY());
	  isGamma = 1;
	  fNGammas++;
	  vGammas.push_back( gamma );
	}	     
      }
      
      if( isGamma ) continue;


      // KShort finder
      
      bool isKs = 0;
      
      if( 
	 t1->GetTPCNcls()>=fKsCuts[0]
	 && t2->GetTPCNcls()>=fKsCuts[0]
	 && dev1>=fKsCuts[1]
	 && dev2>=fKsCuts[1]
	 && devPrim <= fKsCuts[2]
	 && length >= fKsCuts[3]*sigmaLength
	 && length >= fKsCuts[4]
	 && r <= fKsCuts[5]
	 ){	
    
	AliKFParticle piN( *t1, 211 );	
	AliKFParticle piP( *t2, 211 );	
	
	AliKFParticle Ks( piN, piP );
	Ks.SetProductionVertex( primVtx );
	
	double mass, error;
	Ks.GetMass( mass, error);
	if( fKShort ) fKShort->Fill( mass );	
	if( TMath::Abs( mass - kKsMass )<=fKsCuts[6]*error || TMath::Abs( mass - kKsMass )<=fKsCuts[7] ){  
	  isKs = 1;
	  fNKShorts++;
	}
      }
      
      if( isKs ) continue;
      
      // Lambda finder 
     
      if( 
	 t1->GetTPCNcls()>=fLambdaCuts[0]
	 && t2->GetTPCNcls()>=fLambdaCuts[0]
	 && dev1>=fLambdaCuts[1]
	 && dev2>=fLambdaCuts[1]
	 && devPrim <= fLambdaCuts[2]
	 && length >= fLambdaCuts[3]*sigmaLength
	 && length >= fLambdaCuts[4]
	 && r <= fLambdaCuts[5]
	 && TMath::Abs( ap )>.4
	 ){
	
	AliKFParticle kP, kpi;
	if( ap<0 ){ 
	  kP = AliKFParticle( *t2, 2212 );
	  kpi = AliKFParticle( *t1, 211 );
	} else {
	  kP = AliKFParticle( *t1, 2212 );
	  kpi = AliKFParticle( *t2, 211 );
	}

	AliKFParticle lambda = AliKFParticle( kP, kpi );
	lambda.SetProductionVertex( primVtx );	
	double mass, error;
	lambda.GetMass( mass, error);
	if( fLambda ) fLambda->Fill( mass );
	if( TMath::Abs( mass - kLambdaMass )<=fLambdaCuts[6]*error || TMath::Abs( mass - kLambdaMass )<=fLambdaCuts[7] ){
	  fNLambdas++;
	}
      }

    }// V0's


    // Pi0 finder 

    for(UInt_t g1=0;g1<vGammas.size();g1++){
      for(UInt_t g2=g1+1;g2<vGammas.size();g2++){
	AliKFParticle pi0(vGammas.at(g1),vGammas.at(g2));
	double mass, error;
	pi0.GetMass(mass,error);
	fPi0->Fill(mass);
	if( TMath::Abs( mass - kPi0Mass )<=0.03 ){
	  fNPi0s++;
	}
      }
    }

  
    if( fGamma ) PushBack( (TObject*) fGamma, kAliHLTDataTypeHistogram,fUID);
    
    if( fKShort ) PushBack( (TObject*) fKShort, kAliHLTDataTypeHistogram,fUID);
    
    if( fLambda ) PushBack( (TObject*) fLambda, kAliHLTDataTypeHistogram, fUID);
 
    if( fPi0 ) PushBack( (TObject*) fPi0, kAliHLTDataTypeHistogram, fUID);
    
    if( fAP ) PushBack( (TObject*) fAP, kAliHLTDataTypeHistogram, fUID);    

    if( fGammaXY ) PushBack( (TObject*) fGammaXY, kAliHLTDataTypeHistogram, fUID);
  }  
  if( fNPi0s>0 ){
    HLTInfo("Found %d Gammas, %d KShorts, %d Lambdas and %d Pi0's in %d events", fNGammas, fNKShorts, fNLambdas, fNPi0s, fNEvents );    
  }
  else HLTInfo("Found %d Gammas, %d KShorts and %d Lambdas in %d events", fNGammas, fNKShorts, fNLambdas, fNEvents );
  
  return 0;
}

int AliHLTV0HistoComponent::Configure(const char* arguments)
{
  
  int iResult=0;
  if (!arguments) return iResult;
  
  TString allArgs=arguments;
  TString argument;
  
  TObjArray* pTokens=allArgs.Tokenize(" ");
  int bMissingParam=0;

  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;

      if (argument.CompareTo("-cutsGamma")==0) {
	TString spar = "";	
	for( int j=0; j<8; j++ ){
	  if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	  spar+=" ";
	  spar+=((TObjString*)pTokens->At(i))->GetString();
	  fGammaCuts[j]=((TObjString*)pTokens->At(i))->GetString().Atof();
	}
	if( !bMissingParam ){
	  HLTInfo("Gamma cuts are set to: %s", spar.Data());
	  continue;
	}
      } else if (argument.CompareTo("-cutsAP")==0) {
	TString spar = "";	
	for( int j=0; j<8; j++ ){
	  if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	  spar+=" ";
	  spar+=((TObjString*)pTokens->At(i))->GetString();
	  fAPCuts[j]=((TObjString*)pTokens->At(i))->GetString().Atof();
	}
	if( !bMissingParam ){
	  HLTInfo("AP cuts are set to: %s", spar.Data());
	  continue;
	}
      }
      else if (argument.CompareTo("-cutsKs")==0) {
	TString spar = "";	
	for( int j=0; j<8; j++ ){
	  if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	  spar+=" ";
	  spar+=((TObjString*)pTokens->At(i))->GetString();
	  fKsCuts[j]=((TObjString*)pTokens->At(i))->GetString().Atof();
	}
	if( !bMissingParam ){
	  HLTInfo("KShort cuts are set to: %s", spar.Data());
	  continue;
	}
      } else if (argument.CompareTo("-cutsLambda")==0) {
	TString spar = "";	
	for( int j=0; j<8; j++ ){
	  if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	  spar+=" ";
	  spar+=((TObjString*)pTokens->At(i))->GetString();
	  fLambdaCuts[j]=((TObjString*)pTokens->At(i))->GetString().Atof();
	}
	if( !bMissingParam ){
	  HLTInfo("Lampda cuts are set to: %s", spar.Data());
	  continue;
	}
      }else {
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
  const char* path="HLT/ConfigHLT/V0Histo";
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
