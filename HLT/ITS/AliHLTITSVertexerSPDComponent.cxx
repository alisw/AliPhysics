// $Id: AliHLTITSVertexerSPDComponent.cxx 32659 2009-06-02 16:08:40Z sgorbuno $
// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de> *
//                  Ivan Kisel <kisel@kip.uni-heidelberg.de>                *
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************

///  @file   AliHLTITSVertexerSPDComponent.cxx
///  @author Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de>
///  @date   June 2009
///  @brief  An ITS tracker processing component for the HLT


/////////////////////////////////////////////////////
//                                                 //
// a ITS tracker processing component for the HLT  //
//                                                 //
/////////////////////////////////////////////////////

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTITSVertexerSPDComponent.h"
#include "AliHLTArray.h"
#include "AliExternalTrackParam.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "AliITStrackerHLT.h"
#include "AliHLTITSSpacePointData.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliHLTDataTypes.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliGeomManager.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliESDVertex.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp( AliHLTITSVertexerSPDComponent )
AliHLTITSVertexerSPDComponent::AliHLTITSVertexerSPDComponent()
    :
    fSolenoidBz( 0 ),
    fProduceHistos(1),
    fAutoCalibration(1000),
    fFullTime( 0 ),
    fRecoTime( 0 ),
    fNEvents( 0 ),    
    fHistoVertexXY(0),
    fHistoVertexX(0),
    fHistoVertexY(0),
  fHistoVertexZ(0)

{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  fDefRunVtx[0] = 0;
  fDefRunVtx[1] = 0;
  fDefRunVtx[2] = 0;
}

AliHLTITSVertexerSPDComponent::AliHLTITSVertexerSPDComponent( const AliHLTITSVertexerSPDComponent& )
    :
    AliHLTProcessor(),
    fSolenoidBz( 0 ),
    fProduceHistos(1),
    fAutoCalibration(1000),
    fFullTime( 0 ),
    fRecoTime( 0 ),
    fNEvents( 0 ),    
    fHistoVertexXY(0),
    fHistoVertexX(0),
    fHistoVertexY(0),
  fHistoVertexZ(0)
{
  // see header file for class documentation
  HLTFatal( "copy constructor untested" );
}

AliHLTITSVertexerSPDComponent& AliHLTITSVertexerSPDComponent::operator=( const AliHLTITSVertexerSPDComponent& )
{
  // see header file for class documentation
  HLTFatal( "assignment operator untested" );
  return *this;
}

AliHLTITSVertexerSPDComponent::~AliHLTITSVertexerSPDComponent()
{
  // see header file for class documentation
  
  delete fHistoVertexXY;
  delete fHistoVertexX;
  delete fHistoVertexY;
  delete fHistoVertexZ;
}

//
// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process
//

const char* AliHLTITSVertexerSPDComponent::GetComponentID()
{
  // see header file for class documentation
  return "ITSVertexerSPD";
}

void AliHLTITSVertexerSPDComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list )
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD );
}

AliHLTComponentDataType AliHLTITSVertexerSPDComponent::GetOutputDataType()
{
  // see header file for class documentation  
  return kAliHLTMultipleDataType;
}

int AliHLTITSVertexerSPDComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)

{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back( kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS);
  tgtList.push_back(kAliHLTDataTypeHistogram);
  return tgtList.size();
}

void AliHLTITSVertexerSPDComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // define guess for the output data size
  constBase = 10000;       // minimum size
  inputMultiplier = 0.5; // size relative to input
}

AliHLTComponent* AliHLTITSVertexerSPDComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTITSVertexerSPDComponent;
}

void AliHLTITSVertexerSPDComponent::SetDefaultConfiguration()
{
  // Set default configuration for the CA tracker component
  // Some parameters can be later overwritten from the OCDB

  fSolenoidBz = -5.00668;
  fProduceHistos=1;
  fAutoCalibration = 1000;
  fDefRunVtx[0] = 0;
  fDefRunVtx[1] = 0;
  fDefRunVtx[2] = 0;
  fFullTime = 0;
  fRecoTime = 0;
  fNEvents = 0;
}

int AliHLTITSVertexerSPDComponent::ReadConfigurationString(  const char* arguments )
{
  // Set configuration parameters for the CA tracker component from the string

  int iResult = 0;
  if ( !arguments ) return iResult;

  TString allArgs = arguments;
  TString argument;
  int bMissingParam = 0;

  TObjArray* pTokens = allArgs.Tokenize( " " );

  int nArgs =  pTokens ? pTokens->GetEntries() : 0;

  for ( int i = 0; i < nArgs; i++ ) {
    argument = ( ( TObjString* )pTokens->At( i ) )->GetString();
    if ( argument.IsNull() ) continue;

    if ( argument.CompareTo( "-solenoidBz" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fSolenoidBz = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      HLTInfo( "Magnetic Field set to: %f", fSolenoidBz );
      continue;
    }

    if ( argument.CompareTo( "-produceHistos" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fProduceHistos = ( ( TObjString* )pTokens->At( i ) )->GetString().Atoi();
      HLTInfo( "fProduceHistos set to: %d", fProduceHistos );
      continue;
    }

    if ( argument.CompareTo( "-runVertex" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fDefRunVtx[0] = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fDefRunVtx[1] = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fDefRunVtx[2] = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      HLTInfo( "Default run vertex is set to (%f,%f,%f)",fDefRunVtx[0],
	       fDefRunVtx[1],fDefRunVtx[2] );
      for( int i=0; i<3; i++){
	fRunVtx[i] = fDefRunVtx[i];
	fRunVtxNew[i] = fDefRunVtx[i];
      }
      fRunVtx[3] = 0.;
      fRunVtxNew[3] = 0.;
      continue;
    }

    if ( argument.CompareTo( "-beamDiamondCalibration" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fAutoCalibration = ( ( TObjString* )pTokens->At( i ) )->GetString().Atoi();
      HLTInfo( "N events for recalibration of the run vertex is set to: %d", fAutoCalibration );
      continue;
    }
    

    HLTError( "Unknown option \"%s\"", argument.Data() );
    iResult = -EINVAL;
  }
  delete pTokens;

  if ( bMissingParam ) {
    HLTError( "Specifier missed for parameter \"%s\"", argument.Data() );
    iResult = -EINVAL;
  }

  return iResult;
}


int AliHLTITSVertexerSPDComponent::ReadCDBEntry( const char* cdbEntry, const char* chainId )
{
  // see header file for class documentation

  const char* defaultNotify = "";

  if ( !cdbEntry ) {
    return 0;// need to add the HLT/ConfigITS/ITSTracker directory to cdb SG!!!
    cdbEntry = "HLT/ConfigITS/ITSTracker";
    defaultNotify = " (default)";
    chainId = 0;
  }

  HLTInfo( "configure from entry \"%s\"%s, chain id %s", cdbEntry, defaultNotify, ( chainId != NULL && chainId[0] != 0 ) ? chainId : "<none>" );
  AliCDBEntry *pEntry = AliCDBManager::Instance()->Get( cdbEntry );//,GetRunNo());

  if ( !pEntry ) {
    HLTError( "cannot fetch object \"%s\" from CDB", cdbEntry );
    return -EINVAL;
  }

  TObjString* pString = dynamic_cast<TObjString*>( pEntry->GetObject() );

  if ( !pString ) {
    HLTError( "configuration object \"%s\" has wrong type, required TObjString", cdbEntry );
    return -EINVAL;
  }

  HLTInfo( "received configuration object string: \"%s\"", pString->GetString().Data() );

  return  ReadConfigurationString( pString->GetString().Data() );
}


int AliHLTITSVertexerSPDComponent::Configure( const char* cdbEntry, const char* chainId, const char *commandLine )
{
  // Configure the component
  // There are few levels of configuration,
  // parameters which are set on one step can be overwritten on the next step

  //* read hard-coded values

  SetDefaultConfiguration();

  //* read the default CDB entry

  int iResult1 = ReadCDBEntry( NULL, chainId );

  //* read magnetic field

  int iResult2 = ReadCDBEntry( kAliHLTCDBSolenoidBz, chainId );

  //* read the actual CDB entry if required

  int iResult3 = ( cdbEntry ) ? ReadCDBEntry( cdbEntry, chainId ) : 0;

  //* read extra parameters from input (if they are)

  int iResult4 = 0;

  if ( commandLine && commandLine[0] != '\0' ) {
    HLTInfo( "received configuration string from HLT framework: \"%s\"", commandLine );
    iResult4 = ReadConfigurationString( commandLine );
  }

  // Initialise the tracker here

  return iResult1 ? iResult1 : ( iResult2 ? iResult2 : ( iResult3 ? iResult3 : iResult4 ) );
}



int AliHLTITSVertexerSPDComponent::DoInit( int argc, const char** argv )
{
  // Configure the ITS tracker component
  
  fProduceHistos = 1;
  fAutoCalibration = 1000;

  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry("");
  }
  AliGeomManager::ApplyAlignObjsFromCDB("ITS");
  
  TString arguments = "";
  for ( int i = 0; i < argc; i++ ) {
    if ( !arguments.IsNull() ) arguments += " ";
    arguments += argv[i];
  }

  int ret = Configure( NULL, NULL, arguments.Data() );

  if( fProduceHistos ){
    fHistoVertexXY = new TH2F("hITSvertexXY", "ITSSPD vertex in XY", 100,-2,2,100,-2,2);
    fHistoVertexX = new TH1F("hITSvertexX", "ITSSPD vertex X", 100,-2,2);
    fHistoVertexY = new TH1F("hITSvertexY", "ITSSPD vertex Y", 100,-2,2);
    fHistoVertexZ = new TH1F("hITSvertexZ", "ITSSPD vertex Z", 100,-15,15);
  }

  for( int i=0; i<3; i++){
    fRunVtx[i] = fDefRunVtx[i];
    fRunVtxNew[i] = fDefRunVtx[i];
  }
  fRunVtx[3] = 0.;
  fRunVtxNew[3] = 0.;

  return ret;
}


int AliHLTITSVertexerSPDComponent::DoDeinit()
{
  // see header file for class documentation

  delete fHistoVertexXY;
  delete fHistoVertexX;
  delete fHistoVertexY;
  delete fHistoVertexZ;

  fHistoVertexXY = 0;
  fHistoVertexX = 0;
  fHistoVertexY = 0;
  fHistoVertexZ = 0;
  fAutoCalibration = 1000;
  fProduceHistos = 1;
  return 0;
}



int AliHLTITSVertexerSPDComponent::Reconfigure( const char* cdbEntry, const char* chainId )
{
  // Reconfigure the component from OCDB .

  return Configure( cdbEntry, chainId, NULL );
}


int AliHLTITSVertexerSPDComponent::DoEvent
(
  const AliHLTComponentEventData& evtData,
  const AliHLTComponentBlockData* blocks,
  AliHLTComponentTriggerData& /*trigData*/,
  AliHLTUInt8_t* /*outputPtr*/,
  AliHLTUInt32_t& size,
  vector<AliHLTComponentBlockData>& /*outputBlocks*/ )
{
  //* process event

  AliHLTUInt32_t maxBufferSize = size;
  size = 0; // output size

  if (!IsDataEvent()) return 0;

  if ( evtData.fBlockCnt <= 0 ) {
    HLTWarning( "no blocks in event" );
    return 0;
  }


  TStopwatch timer;

  // Event reconstruction in ITS

  int iResult=0;

  int nBlocks = evtData.fBlockCnt;
  int nClustersTotal = 0;
 

  const int kNPhiBins = 20;
  const double kDPhi = TMath::TwoPi() / kNPhiBins;
  double vtxX = fRunVtx[0], vtxY = fRunVtx[1], vtxZ = fRunVtx[2];


  std::vector<AliHLTITSVertexerSPDComponent::AliHLTITSVZCluster> clusters[2][ kNPhiBins ];

  for (int ndx=0; ndx<nBlocks && iResult>=0; ndx++) {

    const AliHLTComponentBlockData* iter = blocks+ndx;
 
    // Read ITS SPD clusters

    if ( iter->fDataType == (kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD) ){

      AliHLTITSClusterData *inPtr=reinterpret_cast<AliHLTITSClusterData*>( iter->fPtr );
      int nClusters = inPtr->fSpacePointCnt;
      nClustersTotal+=nClusters;
      for( int icl=0; icl<nClusters; icl++ ){	
	AliHLTITSSpacePointData &d = inPtr->fSpacePoints[icl];
	if( d.fLayer>1 ) continue;// SPD only;
	if( d.fLayer<0 ) continue;// SPD only;
	Int_t lab[4] = { d.fTracks[0], d.fTracks[1], d.fTracks[2], d.fIndex };
	Int_t info[3] = { d.fNy, d.fNz, d.fLayer };
	Float_t hit[6] = { d.fY, d.fZ, d.fSigmaY2, d.fSigmaZ2, d.fQ, d.fSigmaYZ };
	if( d.fLayer==4 ) hit[5] = -hit[5];


	AliITSRecPoint p( lab, hit, info );

	if (!p.Misalign()){
	  HLTWarning("Can not misalign an SPD cluster");
	}
	Float_t xyz[3];      
	p.GetGlobalXYZ(xyz);      

	// get phi bin
	double phi = TMath::ATan2( xyz[1]-vtxY, xyz[0]-vtxX );
	if( phi<0 ) phi+=TMath::TwoPi();
	int iphi = (int ) phi/kDPhi;
	if( iphi<0 ) iphi = 0;
	if( iphi>=kNPhiBins ) iphi = kNPhiBins-1;
	AliHLTITSVZCluster c;
	c.fX = xyz[0];
	c.fY = xyz[1];
	c.fZ = xyz[2];
	clusters[d.fLayer][iphi].push_back( c );	
     }   
    }    
  }// end read input blocks
  

  // Reconstruct the vertex

  TStopwatch timerReco;
    

  const double kZBinMin = -10;
  const int kNZBins = 40;
  const double kZStep = 1./kNZBins;

  double histX[kNZBins];
  double histY[kNZBins];
  double histZ[kNZBins];
  double histW[kNZBins];
  int histN[kNZBins];
  
  // fit few times with decreasing of the radius cut

  double vtxDR2[2] = { 2.*2., .2*.2};
  double  maxW = 0;  
  int bestBin = -1;
  
  for( int iter = 0; iter<2; iter++ ){
    
    for( int i=0; i<kNZBins; i++ ){
      histX[i] = 0;
      histY[i] = 0;
      histZ[i] = 0;
      histW[i] = 0;
      histN[i] = 0;
    }
 
    for( int iPhi=0; iPhi<kNPhiBins; iPhi++ ){
      int nCluUp = clusters[1][iPhi].size();    
      int nCluDn = clusters[0][iPhi].size();    
      for( int icUp=0; icUp<nCluUp; icUp++ ){
	AliHLTITSVZCluster &cu = clusters[1][iPhi][icUp];
	double x0 = cu.fX - vtxX;
	double y0 = cu.fY - vtxY;
	double z0 = cu.fZ - vtxZ;
	double bestR2 = 1.e10;
	int bestDn=-1, bestBin =0;
	double bestV[3]={0,0,0};
	for( int icDn=0; icDn<nCluDn; icDn++ ){
	  AliHLTITSVZCluster &cd = clusters[0][iPhi][icDn];
	  double x1 = cd.fX - vtxX;
	  double y1 = cd.fY - vtxY;
	  double z1 = cd.fZ - vtxZ;
	  double dx = x1 - x0;
	  double dy = y1 - y0;
	  double l2 = 1./(dx*dx + dy*dy);
	  double a = x1*y0 - x0*y1;
	  double r2 = a*a*l2; 
	  if( r2>vtxDR2[iter] ) continue;
	  double xv = -dy*a*l2;
	  double yv =  dx*a*l2;
	  double zv = ( (x1*z0-x0*z1)*dx + (y1*z0-y0*z1)*dy )*l2;
	  int zbin = (int)((zv - kZBinMin)*kZStep);
	  if( zbin<0 || zbin>=kNZBins ) continue;
	  if( r2<bestR2 ){
	    bestR2 = r2;
	    bestDn = icDn;
	    bestV[0] = xv;
	    bestV[1] = yv;
	    bestV[2] = zv;
	    bestBin = zbin;
	  }
	}
	if( bestDn>=0 ){
	  double w = 1./(1.+bestR2);
	  histX[bestBin]+=bestV[0]*w;
	  histY[bestBin]+=bestV[1]*w;
	  histZ[bestBin]+=bestV[2]*w;
	  histW[bestBin]+=w;
	  histN[bestBin]+=1;
	}
      }
    }
    
    maxW = 0;  
    bestBin = -1;
    for( int i=0; i<kNZBins; i++ ){
      if( histW[i]>maxW ){
	bestBin = i;
	maxW = histW[i];
      }
    }
    if( bestBin<0 || histN[bestBin] <3 ){
      bestBin = -1;
      break;
    }
    vtxX +=histX[bestBin]/maxW;
    vtxY +=histY[bestBin]/maxW;
    vtxZ +=histZ[bestBin]/maxW;    
  }

  if( bestBin>=0 ){
    double nv = 1;
    double nNew = fRunVtxNew[3] + nv;
    double v[3] = {vtxX, vtxY, vtxZ};
    for( int i=0; i<3; i++){
      fRunVtxNew[i] = ( fRunVtxNew[i]*fRunVtxNew[3] + v[i]*nv )/nNew;
    }
    fRunVtxNew[3] = nNew;
  }  
  
  if( fAutoCalibration>0 && fRunVtxNew[3] >= fAutoCalibration ){
    for( int i=0; i<4; i++ ){
      fRunVtx[i] = fRunVtxNew[i];
      fRunVtxNew[i] = 0;
    }
    //cout<<"ITSVertexerSPD: set run vtx to "<<fRunVtx[0]<<" "<<fRunVtx[1]<<" "<<fRunVtx[2]<<endl;  
    if( fRunVtx[0]>3. ) fRunVtx[0] = 3;
    if( fRunVtx[0]<-3. ) fRunVtx[0] = -3;
    if( fRunVtx[1]>3. ) fRunVtx[1] = 3;
    if( fRunVtx[1]<-3. ) fRunVtx[1] = -3;
    if( fRunVtx[2]>30. ) fRunVtx[2] = 30;
    if( fRunVtx[2]<30. ) fRunVtx[2] = -30;
  }

  timerReco.Stop();
  
  // Fill output 

  if( bestBin>=0 ){
    double pos[3] = {vtxX, vtxY, vtxZ};
    double s = 400.E-4;
    double cov[6] = {s*s,0,s*s,0,0,s*s};
    AliESDVertex v(pos, cov, 0, histN[bestBin]);
    PushBack( (TObject*) &v, kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS,0 );

    //cout<<"ITSVertexerSPD: vtx found "<<vtxX<<" "<<vtxY<<" "<<vtxZ<<endl;

    if( fHistoVertexXY ) fHistoVertexXY->Fill( vtxX, vtxY );
    if( fHistoVertexX ) fHistoVertexX->Fill( vtxX );
    if( fHistoVertexY ) fHistoVertexY->Fill( vtxY );
    if( fHistoVertexZ ) fHistoVertexZ->Fill( vtxZ );
  }

  if( fHistoVertexXY ) PushBack( (TObject*) fHistoVertexXY, kAliHLTDataTypeHistogram,0);
  if( fHistoVertexX ) PushBack( (TObject*) fHistoVertexX, kAliHLTDataTypeHistogram,0);
  if( fHistoVertexY ) PushBack( (TObject*) fHistoVertexY, kAliHLTDataTypeHistogram,0);
  if( fHistoVertexZ ) PushBack( (TObject*) fHistoVertexZ, kAliHLTDataTypeHistogram,0);


  timer.Stop();
  fFullTime += timer.RealTime();
  fRecoTime += timerReco.RealTime();
  fNEvents++;

  // Set log level to "Warning" for on-line system monitoring
  int hz = ( int ) ( fFullTime > 1.e-10 ? fNEvents / fFullTime : 100000 );
  int hz1 = ( int ) ( fRecoTime > 1.e-10 ? fNEvents / fRecoTime : 100000 );
  HLTInfo( "ITS Z Vertexer: input %d clusters; time: full %d / reco %d Hz",
	   nClustersTotal, hz, hz1 );

  return iResult;
}
