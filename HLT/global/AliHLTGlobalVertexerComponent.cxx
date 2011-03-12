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
#include "AliHLTDataTypes.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include <TFile.h>
#include <TString.h>
#include "TObjString.h"
#include "TObjArray.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliHLTMessage.h"
#include "TMath.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "TStopwatch.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalVertexerComponent)

AliHLTGlobalVertexerComponent::AliHLTGlobalVertexerComponent()
:
  fNTracks(0),
  fTrackInfos(0),
  fPrimaryVtx(),   
  fFitTracksToVertex(1),
  fConstrainedTrackDeviation(4.),
  fV0DaughterPrimDeviation( 2.5 ),
  fV0PrimDeviation( 3.5 ),
  fV0Chi(3.5),
  fV0DecayLengthInSigmas(3.),
  fV0TimeLimit(10.e-3), 
  fBenchmark("GlobalVertexer")
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
  list.push_back( kAliHLTDataTypeESDObject );
  list.push_back( kAliHLTDataTypeTrack|kAliHLTDataOriginITS );
  list.push_back( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC );
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
  tgtList.push_back( kAliHLTDataTypeGlobalVertexer|kAliHLTDataOriginOut);
  tgtList.push_back( kAliHLTDataTypeESDVertex|kAliHLTDataOriginOut);
  tgtList.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginOut);
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

  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
  fBenchmark.SetTimer(1,"convert");
  fBenchmark.SetTimer(2,"vprim");
  fBenchmark.SetTimer(3,"v0");
  fV0TimeLimit = 10.e-3;

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

  fTrackInfos = 0;
  fFitTracksToVertex = 1;
  fConstrainedTrackDeviation = 4.;
  fV0DaughterPrimDeviation = 2.5 ;
  fV0PrimDeviation =3.5;
  fV0Chi = 3.5;
  fV0DecayLengthInSigmas = 3.;
  fV0TimeLimit = 10.e-3;

  return 0;
}

int AliHLTGlobalVertexerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, 
					    AliHLTComponentTriggerData& /*trigData*/ )
{
  //cout<<"AliHLTGlobalVertexerComponent::DoEvent called"<<endl;
 
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  int iResult = 0;
  //cout<<"SG: GlobalVertexer:: DoEvent called"<<endl;
  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);
  
  vector< AliExternalTrackParam > tracks;
  vector< int > trackId;
  vector< pair<int,int> > v0s;

  AliESDEvent *event = 0; 

  for( const AliHLTComponentBlockData *i= GetFirstInputBlock(kAliHLTDataTypeESDObject); i!=NULL; i=GetNextInputBlock() ){
    fBenchmark.AddInput(i->fSize);
  }


  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {
    event = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
    if( !event ) continue;
    
    event->GetStdContent();
    Int_t nESDTracks=event->GetNumberOfTracks(); 

    for (Int_t iTr=0; iTr<nESDTracks; iTr++){   
      AliESDtrack *pTrack = event->GetTrack(iTr);    
      if( !pTrack  ) continue;
      if( pTrack->GetKinkIndex(0)>0) continue;
      if( !( pTrack->GetStatus()&AliESDtrack::kTPCin ) ) continue;
      tracks.push_back(*pTrack);
      trackId.push_back(iTr);
    }      
    break;
  }

  if( tracks.size()==0 ){
    
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginITS);
	 pBlock!=NULL; pBlock=GetNextInputBlock()) {

      fBenchmark.AddInput(pBlock->fSize);

      AliHLTTracksData* dataPtr = reinterpret_cast<AliHLTTracksData*>( pBlock->fPtr );
      int nTracks = dataPtr->fCount;

      AliHLTExternalTrackParam* currOutTrack = dataPtr->fTracklets;
      for( int itr=0; itr<nTracks; itr++ ){
	AliHLTGlobalBarrelTrack t(*currOutTrack);
	tracks.push_back( t );
	trackId.push_back( currOutTrack->fTrackID );
	unsigned int dSize = sizeof( AliHLTExternalTrackParam ) + currOutTrack->fNPoints * sizeof( unsigned int );
	currOutTrack = ( AliHLTExternalTrackParam* )( (( Byte_t * )currOutTrack) + dSize );
      }
    }
  }

  if( tracks.size()==0 ){
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
	 pBlock!=NULL; pBlock=GetNextInputBlock()) {

      fBenchmark.AddInput(pBlock->fSize);

      AliHLTTracksData* dataPtr = reinterpret_cast<AliHLTTracksData*>( pBlock->fPtr );
      int nTracks = dataPtr->fCount;      
      AliHLTExternalTrackParam* currOutTrack = dataPtr->fTracklets;
      for( int itr=0; itr<nTracks; itr++ ){
	AliHLTGlobalBarrelTrack t(*currOutTrack);
	tracks.push_back( t );
	trackId.push_back( currOutTrack->fTrackID );
	unsigned int dSize = sizeof( AliHLTExternalTrackParam ) + currOutTrack->fNPoints * sizeof( unsigned int );
	currOutTrack = ( AliHLTExternalTrackParam* )( (( Byte_t * )currOutTrack) + dSize );
      }
    }
  }

  
  // primary vertex & V0's 
          
  AliKFParticle::SetField( GetBz() );

  fBenchmark.Start(1);

  {  //* Fill fTrackInfo array

    if( fTrackInfos ) delete[] fTrackInfos;
    fTrackInfos = 0;
    fNTracks=tracks.size(); 
    fTrackInfos = new AliESDTrackInfo[ fNTracks ];
    for (Int_t iTr=0; iTr<fNTracks; iTr++){       
      AliESDTrackInfo &info = fTrackInfos[iTr];
      info.fOK = 1;
      info.fPrimUsedFlag = 0;
      info.fParticle = AliKFParticle( tracks[iTr], 211 );
      for( int i=0; i<8; i++ ) if( !finite(info.fParticle.GetParameter(i)) ) info.fOK = 0;
      for( int i=0; i<36; i++ ) if( !finite(info.fParticle.GetCovariance(i)) ) info.fOK = 0;
    }
  }

  fBenchmark.Stop(1);
  fBenchmark.Start(2);
  FindPrimaryVertex();
  fBenchmark.Stop(2);
  fBenchmark.Start(3);
  FindV0s( v0s );
  fBenchmark.Stop(3);

  int *buf = new int[sizeof(AliHLTGlobalVertexerData)/sizeof(int)+1 + fNTracks + 2*v0s.size()];
  AliHLTGlobalVertexerData *data = reinterpret_cast<AliHLTGlobalVertexerData*>(buf);

  if( data) {  // fill the output structure
        
    data->fFitTracksToVertex = fFitTracksToVertex;
    for( int i=0; i<3; i++ ) data->fPrimP[i] = fPrimaryVtx.Parameters()[i];
    for( int i=0; i<6; i++ ) data->fPrimC[i] = fPrimaryVtx.CovarianceMatrix()[i];
    data->fPrimChi2 = fPrimaryVtx.GetChi2();
    data->fPrimNContributors = fPrimaryVtx.GetNContributors();
    data->fNPrimTracks = 0;
    for( Int_t i = 0; i<fNTracks; i++ ){
      if( !fTrackInfos[i].fPrimUsedFlag ) continue;	  
      if( fTrackInfos[i].fPrimDeviation > fConstrainedTrackDeviation ) continue;
      data->fTrackIndices[ (data->fNPrimTracks)++ ] =  trackId[i];
    }
    int *listV0 = data->fTrackIndices + data->fNPrimTracks;
    data->fNV0s = v0s.size();
    for( int i=0; i<data->fNV0s; i++ ){
      listV0[2*i] = trackId[v0s[i].first];
      listV0[2*i+1] = trackId[v0s[i].second];
    }

    unsigned int blockSize = sizeof(AliHLTGlobalVertexerData) + (data->fNPrimTracks + 2*data->fNV0s)*sizeof(int);

    iResult = PushBack( reinterpret_cast<void*>(data), blockSize, kAliHLTDataTypeGlobalVertexer|kAliHLTDataOriginOut );  
    fBenchmark.AddOutput(blockSize);
  }  
  
  
  // output the vertex if found
  {
    if( iResult==0 && data && data->fPrimNContributors >=3 ){
      AliESDVertex vESD( data->fPrimP, data->fPrimC, data->fPrimChi2, data->fPrimNContributors );
      iResult = PushBack( (TObject*) &vESD, kAliHLTDataTypeESDVertex|kAliHLTDataOriginOut,0 );
      fBenchmark.AddOutput(GetLastObjectSize());
    }
  }

  // output the ESD event
  if( iResult==0 && event && data ){  
    FillESD( event, data ); 
    iResult = PushBack( event, kAliHLTDataTypeESDObject|kAliHLTDataOriginOut, 0);
    fBenchmark.AddOutput(GetLastObjectSize());
  }

  delete[] buf;

  fBenchmark.Stop(0);
  HLTInfo(fBenchmark.GetStatistics());
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
      else if (argument.CompareTo("-v0MinEventRate")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Minimum event rate for V0 finder is set set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	Double_t rate = ((TObjString*)pTokens->At(i))->GetString().Atof();
	fV0TimeLimit = (rate >0 ) ?1./rate :60; // 1 minute maximum time
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

  return 0; // no CDB path is set so far
  /*
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
*/
}


struct AliHLTGlobalVertexerDeviation
{
  int fI; // track index
  float fD; // deviation from vertex

  bool operator<(const AliHLTGlobalVertexerDeviation &a) const { return fD<a.fD; }
};


void AliHLTGlobalVertexerComponent::FindPrimaryVertex()
{
  //* Find event primary vertex

  fPrimaryVtx.Initialize();
  //fPrimaryVtx.SetBeamConstraint(fESD->GetDiamondX(),fESD->GetDiamondY(),0,
  //TMath::Sqrt(fESD->GetSigma2DiamondX()),TMath::Sqrt(fESD->GetSigma2DiamondY()),5.3);

  // select rough region (in sigmas) in which the vertex could be found, all tracks outside these limits are rejected
  // from the primary vertex finding
  fPrimaryVtx.SetBeamConstraint( 0, 0, 0, 3., 3., 5.3 );

  const AliKFParticle **vSelected = new const AliKFParticle*[fNTracks]; //* Selected particles for the vertex fit
  AliHLTGlobalVertexerDeviation *dev = new AliHLTGlobalVertexerDeviation[fNTracks];  

  Int_t nSelected = 0;
  
  for( Int_t i = 0; i<fNTracks; i++){ 
    if(!fTrackInfos[i].fOK ) continue;
    //if( fESD->GetTrack(i)->GetTPCNcls()<60  ) continue;
    const AliKFParticle &p = fTrackInfos[i].fParticle;
    Double_t chi = p.GetDeviationFromVertex( fPrimaryVtx );      
    if( chi > fConstrainedTrackDeviation ) continue;
    dev[nSelected].fI = i; 
    dev[nSelected].fD = chi; 
    vSelected[nSelected] = &(fTrackInfos[i].fParticle);
    nSelected++;  
  }    

  // fit

  while( nSelected>2 ){ 

    //* Primary vertex finder with rejection of outliers

    for( Int_t i = 0; i<nSelected; i++){ 
      vSelected[i] = &(fTrackInfos[dev[i].fI].fParticle);      
    }

    double xv = fPrimaryVtx.GetX();
    double yv = fPrimaryVtx.GetY();
    double zv = fPrimaryVtx.GetZ(); // values from previous iteration of calculations          
    fPrimaryVtx.Initialize();
    fPrimaryVtx.SetBeamConstraint( 0, 0, 0, 3., 3., 5.3 );
    fPrimaryVtx.SetVtxGuess( xv, yv, zv );    

    fPrimaryVtx.Construct( vSelected, nSelected, 0, -1, 1 ); // refilled for every iteration
    
    for( Int_t it=0; it<nSelected; it++ ){ 
      const AliKFParticle &p = fTrackInfos[dev[it].fI].fParticle;
      if( nSelected <= 20 ){
	AliKFVertex tmp =  fPrimaryVtx - p; // exclude the current track from the sample and recalculate the vertex
	dev[it].fD = p.GetDeviationFromVertex( tmp );
      } else {
	dev[it].fD = p.GetDeviationFromVertex( fPrimaryVtx );	
      }
    }
    sort(dev,dev+nSelected); // sort tracks with increasing chi2 (used for rejection)
       
    int nRemove = (int) ( 0.3*nSelected );  //remove 30% of the tracks (done for performance, only if there are more than 20 tracks)  
    if( nSelected - nRemove <=20 ) nRemove = 1;  // removal based on the chi2 of every track   
    int firstRemove = nSelected - nRemove;
    while( firstRemove<nSelected ){
      if( dev[firstRemove].fD >= fConstrainedTrackDeviation ) break;
      firstRemove++;
    }
    if( firstRemove>=nSelected ) break;
    nSelected = firstRemove;
  }

  for( Int_t i = 0; i<fNTracks; i++){ 
    fTrackInfos[i].fPrimUsedFlag = 0;
  }

  if( nSelected < 3 ){  // no vertex for fewer than 3 contributors  
    fPrimaryVtx.NDF() = -3;
    fPrimaryVtx.Chi2() = 0;
    nSelected = 0;
  }
  
  for( Int_t i = 0; i<nSelected; i++){ 
    AliESDTrackInfo &info = fTrackInfos[dev[i].fI];
    info.fPrimUsedFlag = 1;
    info.fPrimDeviation = dev[i].fD;
  }

  for( Int_t i = 0; i<fNTracks; i++ ){
    AliESDTrackInfo &info = fTrackInfos[i];
    if( info.fPrimUsedFlag ) continue;
    info.fPrimDeviation = info.fParticle.GetDeviationFromVertex( fPrimaryVtx );   
  }

  delete[] vSelected;
  delete[] dev;
}



void AliHLTGlobalVertexerComponent::FindV0s( vector<pair<int,int> > &v0s  )
{
  //* V0 finder

  AliKFVertex &primVtx = fPrimaryVtx;
  if( primVtx.GetNContributors()<3 ) return;

  TStopwatch timer;
  Int_t statN = 0;
  Bool_t run = 1;

  for( Int_t iTr = 0; iTr<fNTracks && run; iTr++ ){ //* first daughter

    AliESDTrackInfo &info = fTrackInfos[iTr];
    if( !info.fOK ) continue;    
    if( info.fParticle.GetQ() >0 ) continue;    
    if( info.fPrimDeviation < fV0DaughterPrimDeviation ) continue;

    for( Int_t jTr = 0; jTr<fNTracks; jTr++ ){  //* second daughter
      
      
      AliESDTrackInfo &jnfo = fTrackInfos[jTr];
      if( !jnfo.fOK ) continue;
      if( jnfo.fParticle.GetQ() < 0 ) continue;
      if( jnfo.fPrimDeviation < fV0DaughterPrimDeviation ) continue;

      // check the time once a while...

      if( (++statN)%100 ==0 ){ 
	if( timer.RealTime()>= fV0TimeLimit ){  run = 0; break; }
	timer.Continue();
      }

      //* check if the particles fit

      if( info.fParticle.GetDeviationFromParticle(jnfo.fParticle) > fV0Chi ) continue;

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
      
      //* keep v0 

      v0s.push_back(pair<int,int>(iTr,jTr));
    }
  }
}




void AliHLTGlobalVertexerComponent::FillESD( AliESDEvent *event, AliHLTGlobalVertexerData *data
)
{
  //* put output of a vertexer to the esd event

  Int_t nESDTracks = event->GetNumberOfTracks();

  const int *listPrim = data->fTrackIndices;
  const int *listV0 = data->fTrackIndices + data->fNPrimTracks;

  std::map<int,int> mapId;
  bool *constrainedToVtx   = new bool[nESDTracks];

  for( int i=0; i<nESDTracks; i++ ){
    constrainedToVtx[i] = 0;
    if( !event->GetTrack(i) ) continue;
    mapId[ event->GetTrack(i)->GetID() ] = i;
  }

  if( data->fPrimNContributors >=3 ){

    AliESDVertex vESD( data->fPrimP, data->fPrimC, data->fPrimChi2, data->fPrimNContributors );
    event->SetPrimaryVertexTPC( &vESD );
    event->SetPrimaryVertexTracks( &vESD );

    // relate tracks to the primary vertex

    if( data->fFitTracksToVertex ){
      for( Int_t i = 0; i<data->fNPrimTracks; i++ ){
	Int_t id = listPrim[ i ];
	map<int,int>::iterator it = mapId.find(id);
	if( it==mapId.end() ) continue;
	Int_t itr = it->second;
	event->GetTrack(itr)->RelateToVertex( &vESD, event->GetMagneticField(),100. );
	constrainedToVtx[ itr ] = 1;
      }
    }
  }

  //* add ESD v0s and relate tracks to v0s
 

  for( int i=0; i<data->fNV0s; i++ ){

    Int_t id1 = listV0[ 2*i ];
    Int_t id2 = listV0[ 2*i + 1];
    map<int,int>::iterator it = mapId.find(id1);
    if( it==mapId.end() ) continue;
    Int_t iTr = it->second;
    it = mapId.find(id2);
    if( it==mapId.end() ) continue;
    Int_t jTr = it->second;
   
    AliESDv0 v0( *event->GetTrack( iTr ), iTr, *event->GetTrack( jTr ), jTr );  
    event->AddV0( &v0 );

    // relate the tracks to the vertex  

    if( data->fFitTracksToVertex ){
      if( constrainedToVtx[iTr] || constrainedToVtx[jTr] ) continue;
      double pos[3];
      double sigma[3] = {.1,.1,.1};
      v0.XvYvZv(pos);
      AliESDVertex vESD(pos, sigma);
      event->GetTrack(iTr)->RelateToVertex( &vESD, event->GetMagneticField(),100. );
      event->GetTrack(jTr)->RelateToVertex( &vESD, event->GetMagneticField(),100. );
      constrainedToVtx[iTr] = 1;
      constrainedToVtx[jTr] = 1;    
    }
  }

  delete[] constrainedToVtx;
}
