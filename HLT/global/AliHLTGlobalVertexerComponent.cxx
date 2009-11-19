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

/** @file   AliHLTGlobalVertexerComponent.cxx
    @author Sergey Gorbunov
    @brief  Component for reconstruct primary vertex and V0's
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTGlobalVertexerComponent.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include <TFile.h>
#include <TString.h>
#include "TObjString.h"
#include "TObjArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliHLTMessage.h"
#include "TMath.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalVertexerComponent)

AliHLTGlobalVertexerComponent::AliHLTGlobalVertexerComponent()
:
  fESD(0),
  fTrackInfos(0),
  fPrimaryVtx(),  
  fHistPrimVertexXY(0),
  fHistPrimVertexZX(0),
  fHistPrimVertexZY(0),
  fNEvents(0),
  fPlotHistograms(0),
  fFitTracksToVertex(1),
  fConstrainedTrackDeviation(4.),
  fV0DaughterPrimDeviation( 2.5 ),
  fV0PrimDeviation( 3.5 ),
  fV0Chi(3.5),
  fV0DecayLengthInSigmas(3.)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTGlobalVertexerComponent::~AliHLTGlobalVertexerComponent()
{
  // see header file for class documentation

  if( fTrackInfos ) delete[] fTrackInfos;
  if( fHistPrimVertexXY ) delete   fHistPrimVertexXY;
  if( fHistPrimVertexZX ) delete   fHistPrimVertexZX;
  if( fHistPrimVertexZY ) delete   fHistPrimVertexZY;
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTGlobalVertexerComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "GlobalVertexer";
}

void AliHLTGlobalVertexerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginOut );
}

AliHLTComponentDataType AliHLTGlobalVertexerComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTGlobalVertexerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)

{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginOut);
  tgtList.push_back(kAliHLTDataTypeHistogram);
  return tgtList.size();
}

void AliHLTGlobalVertexerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 80000;
  inputMultiplier = 2.;
}

AliHLTComponent* AliHLTGlobalVertexerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTGlobalVertexerComponent;
}

int AliHLTGlobalVertexerComponent::DoInit( int argc, const char** argv )
{
  // init
  
  fHistPrimVertexXY = new TH2F("primVertexXY","HLT:  Primary vertex distribution in XY",60,-2.,2.,60,-2.,2.);
  fHistPrimVertexXY->SetMarkerStyle(8);
  fHistPrimVertexXY->SetMarkerSize(0.4);
  fHistPrimVertexXY->SetYTitle("Y");
  fHistPrimVertexXY->SetXTitle("X");
  fHistPrimVertexXY->SetStats(0);
  //fHistPrimVertexXY->SetOption("COLZ");

  fHistPrimVertexZX = new TH2F("primVertexZX","HLT:  Primary vertex distribution in ZX",60,-15.,15.,60,-2.,2.);
  fHistPrimVertexZX->SetMarkerStyle(8);
  fHistPrimVertexZX->SetMarkerSize(0.4);
  fHistPrimVertexZX->SetYTitle("X");
  fHistPrimVertexZX->SetXTitle("Z");
  fHistPrimVertexZX->SetStats(0);
  //fHistPrimVertexZX->SetOption("COLZ");

  fHistPrimVertexZY = new TH2F("primVertexZY","HLT:  Primary vertex distribution in ZY",60,-10.,10.,60,-2.,2.);
  fHistPrimVertexZY->SetMarkerStyle(8);
  fHistPrimVertexZY->SetMarkerSize(0.4);
  fHistPrimVertexZY->SetYTitle("Y");
  fHistPrimVertexZY->SetXTitle("Z");
  fHistPrimVertexZY->SetStats(0);
  //fHistPrimVertexZY->SetOption("COLZ");

  fNEvents =0;

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
  
int AliHLTGlobalVertexerComponent::DoDeinit()
{
  // see header file for class documentation

  if( fTrackInfos ) delete[] fTrackInfos;
  if( fHistPrimVertexXY ) delete fHistPrimVertexXY;
  if( fHistPrimVertexZX ) delete fHistPrimVertexZX;
  if( fHistPrimVertexZY ) delete fHistPrimVertexZY;

  fTrackInfos = 0;
  fHistPrimVertexXY = 0;
  fHistPrimVertexZX = 0;
  fHistPrimVertexZY = 0;

  fPlotHistograms = 0;
  fFitTracksToVertex = 1;
  fConstrainedTrackDeviation = 4.;
  fV0DaughterPrimDeviation = 2.5 ;
  fV0PrimDeviation =3.5;
  fV0Chi = 3.5;
  fV0DecayLengthInSigmas = 3.;
  fNEvents = 0;
  return 0;
}

int AliHLTGlobalVertexerComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  int iResult = 0;

  fNEvents++;

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {

    AliESDEvent *event = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
    if( !event ) continue;
    event->GetStdContent();

    // primary vertex & V0's 
          
    SetESD( event );
    FindPrimaryVertex();
    FindV0s();
    const AliESDVertex *vPrim = event->GetPrimaryVertexTracks();

    // plot histograms

    if( fPlotHistograms ){
      if( vPrim && vPrim->GetNContributors()>=3 ){
	if( fHistPrimVertexXY ) fHistPrimVertexXY->Fill( vPrim->GetX(),vPrim->GetY()  );
	if( fHistPrimVertexZX ) fHistPrimVertexZX->Fill( vPrim->GetZ(),vPrim->GetX()  );
	if( fHistPrimVertexZY ) fHistPrimVertexZY->Fill( vPrim->GetZ(),vPrim->GetY()  );    
      }  
      if( fHistPrimVertexXY ) PushBack( fHistPrimVertexXY, kAliHLTDataTypeHistogram,0);
      if( fHistPrimVertexZX ) PushBack( fHistPrimVertexZX, kAliHLTDataTypeHistogram,0);
      if( fHistPrimVertexZY ) PushBack( fHistPrimVertexZY, kAliHLTDataTypeHistogram,0);
    }
    
    iResult = PushBack( event, kAliHLTDataTypeESDObject|kAliHLTDataOriginOut, 0);
    if( iResult<0 ) break;
  }
    
  return iResult;
}

int AliHLTGlobalVertexerComponent::Configure(const char* arguments)
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
      
      if (argument.CompareTo("-fitTracksToVertex")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("fitTracksToVertex is set set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fFitTracksToVertex=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      }
      else if (argument.CompareTo("-plotHistograms")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("plotHistograms is set set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fPlotHistograms=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      }
      else if (argument.CompareTo("-constrainedTrackDeviation")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("constrainedTrackDeviation is set set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fConstrainedTrackDeviation=((TObjString*)pTokens->At(i))->GetString().Atof();
	continue;
      }
      else if (argument.CompareTo("-v0DaughterPrimDeviation")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("v0DaughterPrimDeviation is set set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fV0DaughterPrimDeviation=((TObjString*)pTokens->At(i))->GetString().Atof();
	continue;
      }
      else if (argument.CompareTo("-v0PrimDeviation")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("v0PrimDeviation is set set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fV0PrimDeviation=((TObjString*)pTokens->At(i))->GetString().Atof();
	continue;
      }
      else if (argument.CompareTo("-v0Chi")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("v0Chi is set set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fV0Chi=((TObjString*)pTokens->At(i))->GetString().Atof();
	continue;
      }
      else if (argument.CompareTo("-v0DecayLengthInSigmas")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("v0DecayLengthInSigmas is set set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fV0DecayLengthInSigmas=((TObjString*)pTokens->At(i))->GetString().Atof();
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
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }  
  
  return iResult;
}

int AliHLTGlobalVertexerComponent::Reconfigure(const char* cdbEntry, const char* chainId)
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



void AliHLTGlobalVertexerComponent::SetESD( AliESDEvent *event )
{
  //* Fill fTrackInfo array

  if( fTrackInfos ) delete[] fTrackInfos;
  fTrackInfos = 0;
  fESD = event;

  AliKFParticle::SetField( fESD->GetMagneticField() );

  Int_t nESDTracks=event->GetNumberOfTracks(); 

  fTrackInfos = new AliESDTrackInfo[ nESDTracks ];

  for (Int_t iTr=0; iTr<nESDTracks; iTr++){ 
  
    AliESDTrackInfo &info = fTrackInfos[iTr];
    info.fOK = 0;
    info.fPrimUsedFlag = 0;
    
    //* track quality check

    AliESDtrack *pTrack = event->GetTrack(iTr);    
    if( !pTrack  ) continue;
    if (pTrack->GetKinkIndex(0)>0) continue;
    if ( !( pTrack->GetStatus()&AliESDtrack::kTPCin ) ) continue;
    
    //* Construct KFParticle for the track

    if(  pTrack->GetStatus()&AliESDtrack::kITSin ){
      info.fParticle = AliKFParticle( *pTrack, 211 );    
    } else {
      info.fParticle = AliKFParticle( *pTrack->GetInnerParam(), 211 );    
    }
    info.fOK = 1;
  }
}


void AliHLTGlobalVertexerComponent::FindPrimaryVertex(  )
{
  //* Find event primary vertex

  int nTracks = fESD->GetNumberOfTracks();

  const AliKFParticle **vSelected = new const AliKFParticle*[nTracks]; //* Selected particles for vertex fit
  Int_t *vIndex = new int [nTracks];                    //* Indices of selected particles
  Bool_t *vFlag = new bool [nTracks];                    //* Flags returned by the vertex finder

  fPrimaryVtx.Initialize();
  fPrimaryVtx.SetBeamConstraint(fESD->GetDiamondX(),fESD->GetDiamondY(),0,
				TMath::Sqrt(fESD->GetSigma2DiamondX()),TMath::Sqrt(fESD->GetSigma2DiamondY()),5.3);
  
  Int_t nSelected = 0;
  for( Int_t i = 0; i<nTracks; i++){ 
    if(!fTrackInfos[i].fOK ) continue;
    //if( fESD->GetTrack(i)->GetTPCNcls()<60  ) continue;
    const AliKFParticle &p = fTrackInfos[i].fParticle;
    Double_t chi = p.GetDeviationFromVertex( fPrimaryVtx );      
    if( chi > fConstrainedTrackDeviation ) continue;
    vSelected[nSelected] = &(fTrackInfos[i].fParticle);
    vIndex[nSelected] = i;
    nSelected++;  
  }
  fPrimaryVtx.ConstructPrimaryVertex( vSelected, nSelected, vFlag, fConstrainedTrackDeviation );
  for( Int_t i = 0; i<nSelected; i++){ 
    if( vFlag[i] ) fTrackInfos[vIndex[i]].fPrimUsedFlag = 1;
  }

  for( Int_t i = 0; i<nTracks; i++ ){
    AliESDTrackInfo &info = fTrackInfos[i];
    info.fPrimDeviation = info.fParticle.GetDeviationFromVertex( fPrimaryVtx );   
  }
  
  if( fPrimaryVtx.GetNContributors()>=3 ){
    AliESDVertex vESD( fPrimaryVtx.Parameters(), fPrimaryVtx.CovarianceMatrix(), fPrimaryVtx.GetChi2(), fPrimaryVtx.GetNContributors() );
    fESD->SetPrimaryVertexTracks( &vESD );

    // relate the tracks to vertex

    if( fFitTracksToVertex ){      
      for( Int_t i = 0; i<nTracks; i++ ){
	if( !fTrackInfos[i].fPrimUsedFlag ) continue;	  
	if( fTrackInfos[i].fPrimDeviation > fConstrainedTrackDeviation ) continue;
	fESD->GetTrack(i)->RelateToVertex( &vESD, fESD->GetMagneticField(),100. );
      }
    }

  } else {
    for( Int_t i = 0; i<nTracks; i++)
      fTrackInfos[i].fPrimUsedFlag = 0;
  }


  delete[] vSelected;
  delete[] vIndex;
  delete[] vFlag;
}


void AliHLTGlobalVertexerComponent::FindV0s(  )
{
  //* V0 finder

  int nTracks = fESD->GetNumberOfTracks();
  //AliKFVertex primVtx( *fESD->GetPrimaryVertexTracks() );
  AliKFVertex &primVtx = fPrimaryVtx;
  if( primVtx.GetNContributors()<3 ) return;

  bool *constrainedV0   = new bool[nTracks];
  for( Int_t iTr = 0; iTr<nTracks; iTr++ ){ 
    constrainedV0[iTr] = 0;
  }
  

  for( Int_t iTr = 0; iTr<nTracks; iTr++ ){ //* first daughter

    AliESDTrackInfo &info = fTrackInfos[iTr];
    if( !info.fOK ) continue;    
    if( info.fParticle.GetQ() >0 ) continue;    
    if( info.fPrimDeviation < fV0DaughterPrimDeviation ) continue;

    for( Int_t jTr = 0; jTr<nTracks; jTr++ ){  //* second daughter

      AliESDTrackInfo &jnfo = fTrackInfos[jTr];
      if( !jnfo.fOK ) continue;
      if( jnfo.fParticle.GetQ() < 0 ) continue;
      if( jnfo.fPrimDeviation < fV0DaughterPrimDeviation ) continue;
   
      //* construct V0 mother

      AliKFParticle v0( info.fParticle, jnfo.fParticle );     

      //* check V0 Chi^2
      
      if( v0.GetChi2()<0 || v0.GetChi2() > fV0Chi*fV0Chi*v0.GetNDF() ) continue;

      //* subtruct daughters from primary vertex 

      AliKFVertex primVtxCopy = primVtx;    
       
      if( info.fPrimUsedFlag ){	
	if( primVtxCopy.GetNContributors()<=2 ) continue;
	primVtxCopy -= info.fParticle;
      }
      if( jnfo.fPrimUsedFlag ){
	if( primVtxCopy.GetNContributors()<=2 ) continue;
	primVtxCopy -= jnfo.fParticle;
      }
      //* Check v0 Chi^2 deviation from primary vertex 

      if( v0.GetDeviationFromVertex( primVtxCopy ) > fV0PrimDeviation ) continue;

      //* Add V0 to primary vertex to improve the primary vertex resolution

      primVtxCopy += v0;      

      //* Set production vertex for V0

      v0.SetProductionVertex( primVtxCopy );

      //* Get V0 decay length with estimated error

      Double_t length, sigmaLength;
      if( v0.GetDecayLength( length, sigmaLength ) ) continue;

      //* Reject V0 if it decays too close[sigma] to the primary vertex

      if( length  < fV0DecayLengthInSigmas*sigmaLength ) continue;

      //* add ESD v0 
      
      AliESDv0 v0ESD( *fESD->GetTrack( iTr ), iTr, *fESD->GetTrack( jTr ), jTr );  
      fESD->AddV0( &v0ESD );

      // relate the tracks to vertex

      if( fFitTracksToVertex ){
	if( constrainedV0[iTr] || constrainedV0[jTr]
	    || info.fPrimDeviation < fConstrainedTrackDeviation || jnfo.fPrimDeviation < fConstrainedTrackDeviation ) continue;
	AliESDVertex vESD(v0.Parameters(), v0.CovarianceMatrix(), v0.GetChi2(), 2);
	fESD->GetTrack(iTr)->RelateToVertex( &vESD, fESD->GetMagneticField(),100. );
	fESD->GetTrack(jTr)->RelateToVertex( &vESD, fESD->GetMagneticField(),100. );
	constrainedV0[iTr] = 1;
	constrainedV0[jTr] = 1;	
      }
    }
  }
  delete[] constrainedV0;
}

