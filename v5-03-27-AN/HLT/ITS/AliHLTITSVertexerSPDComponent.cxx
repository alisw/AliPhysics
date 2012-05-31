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
#include "AliESDVertex.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp( AliHLTITSVertexerSPDComponent )
AliHLTITSVertexerSPDComponent::AliHLTITSVertexerSPDComponent()
    :
    fZRange(40),
    fZBinSize(1),
    fFullTime( 0 ),
    fRecoTime( 0 ),
    fNEvents( 0 ),
    fSumW(0),
    fSumN(0),
    fZMin(0),
    fNZBins(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  for( int i=0; i<9; i++ ) fSum[i] = 0;

  fRunVtx[0] = 0;
  fRunVtx[1] = 0;
  fRunVtx[2] = 0;
}

AliHLTITSVertexerSPDComponent::AliHLTITSVertexerSPDComponent( const AliHLTITSVertexerSPDComponent& )
    :
    AliHLTProcessor(),
    fZRange(40),
    fZBinSize(1),
    fFullTime( 0 ),
    fRecoTime( 0 ),
    fNEvents( 0 ),
    fSumW(0),
    fSumN(0),
    fZMin(0),
    fNZBins(0)
{
  // see header file for class documentation

  for( int i=0; i<9; i++ ) fSum[i] = 0;

  fRunVtx[0] = 0;
  fRunVtx[1] = 0;
  fRunVtx[2] = 0;

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
  for( int i=0; i<9; i++ ){
    delete[] fSum[i];
    fSum[i] = 0;
  }
  delete[] fSumW;
  delete[] fSumN;
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

  fRunVtx[0] = 0;
  fRunVtx[1] = 0;
  fRunVtx[2] = 0;
  fZRange = 40;
  fZBinSize = 1;
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

    if ( argument.CompareTo( "-runVertex" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fRunVtx[0] = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fRunVtx[1] = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fRunVtx[2] = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      HLTInfo( "Default run vertex is set to (%f,%f,%f)",fRunVtx[0],
	       fRunVtx[1],fRunVtx[2] );
      continue;
    }

    if ( argument.CompareTo( "-zRange" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fZRange = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      HLTInfo( "Z range for the vertex search is set to +-%f cm", fZRange );
      continue;
    }

    if ( argument.CompareTo( "-zBinSize" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fZBinSize = ( ( TObjString* )pTokens->At( i ) )->GetString().Atof();
      HLTInfo( "Size of the Z bin for the vertex search is set to %f cm", fZBinSize );
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
    //cdbEntry = "HLT/ConfigITS/ITSTracker";
    //defaultNotify = " (default)";
    //chainId = 0;
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

  //* read the actual CDB entry if required

  int iResult2 = ( cdbEntry ) ? ReadCDBEntry( cdbEntry, chainId ) : 0;
  
  //* read extra parameters from input (if they are)

  int iResult3 = 0;

  if ( commandLine && commandLine[0] != '\0' ) {
    HLTInfo( "received configuration string from HLT framework: \"%s\"", commandLine );
    iResult3 = ReadConfigurationString( commandLine );
  }

  // Initialise the tracker here

  return iResult1 ? iResult1 : ( iResult2 ? iResult2 : iResult3 );
}



int AliHLTITSVertexerSPDComponent::DoInit( int argc, const char** argv )
{
  // Configure the ITS tracker component
  
  SetDefaultConfiguration();

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

  for( int i=0; i<9; i++ ) delete[] fSum[i];
  delete[] fSumW;
  delete[] fSumN;

  fZMin = -fZRange;
  fNZBins = ( (int) (2*fZRange/fZBinSize ))+1;

  for( int i=0; i<9; i++ ) fSum[i] = new double [fNZBins];

  fSumW = new double [fNZBins];
  fSumN = new int [fNZBins];

  return ret;
}


int AliHLTITSVertexerSPDComponent::DoDeinit()
{
  // see header file for class documentation

  for( int i=0; i<9; i++ ){    
    delete[] fSum[i];
    fSum[i] = 0;
  }
  delete[] fSumW;
  delete[] fSumN;
  fSumW = 0;
  fSumN = 0;
  
  SetDefaultConfiguration();

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

  //AliHLTUInt32_t maxBufferSize = size;
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
	int iphi = (int ) (phi/kDPhi);
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
    
  double zScale = 1./fZBinSize;

  // fit few times decreasing the radius cut

  double vtxDR2[5] = { 1.*1., .5*.5, .4*.4, .4*.4, .2*.2};
  double vtxDZ [5] = { 1., .5, .3, .3, .3};
  double  maxW = 0;
  int bestBin = -1;
  
  for( int iter = 0; iter<3; iter++ ){

    bool doZBins = (iter<2);
    int nSearchBins = doZBins ?fNZBins :1;
    {
      for( int i=0; i<nSearchBins; i++ ){
	fSum[0][i] = 0;
	fSum[1][i] = 0;
	fSum[2][i] = 0;
	fSum[3][i] = 0;
	fSum[4][i] = 0;
	fSum[5][i] = 0;
	fSum[6][i] = 0;
	fSum[7][i] = 0;
	fSum[8][i] = 0;
	fSumW[i] = 0;
	fSumN[i] = 0;
      }
    }

    for( int iPhi=0; iPhi<kNPhiBins; iPhi++ ){
      int nCluUp = clusters[1][iPhi].size();    
      int nCluDn = clusters[0][iPhi].size();    
      for( int icUp=0; icUp<nCluUp; icUp++ ){
	AliHLTITSVZCluster &cu = clusters[1][iPhi][icUp];
	double x0 = cu.fX - vtxX;
	double y0 = cu.fY - vtxY;
	double z0 = cu.fZ;// - vtxZ;
	double bestR2 = 1.e10;
	int bestDn=-1, bestBinDn =0;
	double bestP[3]={0,0,0};
	for( int icDn=0; icDn<nCluDn; icDn++ ){
	  AliHLTITSVZCluster &cd = clusters[0][iPhi][icDn];
	  double x1 = cd.fX - vtxX;
	  double y1 = cd.fY - vtxY;
	  double z1 = cd.fZ;// - vtxZ;
	  double dx = x1 - x0;
	  double dy = y1 - y0;
	  double l2 = 1./(dx*dx + dy*dy);
	  double a = x1*y0 - x0*y1;
	  double r2 = a*a*l2; 
	  if( r2 > vtxDR2[iter] ) continue;
	  if( r2 > bestR2 ) continue;
	  //double xv = -dy*a*l2;
	  //double yv =  dx*a*l2;
	  double zv = ( (x1*z0-x0*z1)*dx + (y1*z0-y0*z1)*dy )*l2;
	  int zbin;
	  if( doZBins ){
	    zbin = (int)((zv - fZMin)*zScale);
	    if( zbin<0 || zbin>=fNZBins ) continue;	  
	  } else {
	    zbin = 0;
	    if( TMath::Abs( zv - vtxZ ) > vtxDZ[ iter ] ) continue;
	  }
	  bestR2 = r2;
	  bestDn = icDn;
	  bestP[0] = x1;
	  bestP[1] = y1;
	  bestP[2] = z1;
	  bestBinDn = zbin;	  
	}
	if( bestDn < 0 ) continue;

	double dx = bestP[0] - x0;
	double dy = bestP[1] - y0;
	double dz = bestP[2] - z0;	  
	double w = 1./(1.+bestR2);
	  
	// Equations:
	//
	// A_i*x + B_i*y + C_i*z = D_i
	//
	// Sum[0] = sum A_i*A_i
	// Sum[1] = sum B_i*A_i
	// Sum[2] = sum B_i*B_i
	// Sum[3] = sum C_i*A_i
	// Sum[4] = sum C_i*B_i
	// Sum[5] = sum C_i*C_i	  
	// Sum[6] = sum A_i*D_i
	// Sum[7] = sum B_i*D_i
	// Sum[8] = sum C_i*D_i
	
	double n = w / TMath::Sqrt(dx*dx + dy*dy + dz*dz);	  
	dy*=n;
	dx*=n;
	dz*=n;
	double a, b, c, d;
	a = dy; b = -dx; c = 0; d = dy*x0 - dx*y0;
	{
	  fSum[0][bestBinDn]+= a*a;
	  fSum[1][bestBinDn]+= b*a;
	  fSum[2][bestBinDn]+= b*b;
	  //fSum[3][bestBinDn]+= c*a;
	  //fSum[4][bestBinDn]+= c*b;
	  //fSum[5][bestBinDn]+= c*c;
	  fSum[6][bestBinDn]+= a*d;
	  fSum[7][bestBinDn]+= b*d;
	  //fSum[8][bestBinDn]+= c*d;
	}
	a = dz; b = 0; c = -dx; d = dz*x0 - dx*z0;
	{
	  fSum[0][bestBinDn]+= a*a;
	  //fSum[1][bestBinDn]+= b*a;
	  //fSum[2][bestBinDn]+= b*b;
	  fSum[3][bestBinDn]+= c*a;
	  //fSum[4][bestBinDn]+= c*b;
	  fSum[5][bestBinDn]+= c*c;
	  fSum[6][bestBinDn]+= a*d;
	  //fSum[7][bestBinDn]+= b*d;
	  fSum[8][bestBinDn]+= c*d;
	}
	
	a = 0; b = dz; c = -dy; d = dz*y0 - dy*z0;
	{
	  //fSum[0][bestBinDn]+= a*a;
	  //fSum[1][bestBinDn]+= b*a;
	  fSum[2][bestBinDn]+= b*b;
	  //fSum[3][bestBinDn]+= c*a;
	  fSum[4][bestBinDn]+= c*b;
	  fSum[5][bestBinDn]+= c*c;
	  //fSum[6][bestBinDn]+= a*d;
	  fSum[7][bestBinDn]+= b*d;
	  fSum[8][bestBinDn]+= c*d;
	}

	a = dz; b = 0; c = -dx; d = dz*x0 - dx*z0;
	{
	  fSum[0][bestBinDn]+= a*a;
	  //fSum[1][bestBinDn]+= b*a;
	  //fSum[2][bestBinDn]+= b*b;
	  fSum[3][bestBinDn]+= c*a;
	  //fSum[4][bestBinDn]+= c*b;
	  fSum[5][bestBinDn]+= c*c;
	  fSum[6][bestBinDn]+= a*d;
	  //fSum[7][bestBinDn]+= b*d;
	  fSum[8][bestBinDn]+= c*d;
	}
	
	fSumW[bestBinDn]+=w;
	fSumN[bestBinDn]+=1;
      }    
    }
    
    maxW = 0;  
    bestBin = -1;
    for( int i=0; i<nSearchBins; i++ ){
      if( fSumN[i]>maxW ){
	bestBin = i;
	maxW = fSumN[i];
      }
    }
    if( bestBin<0 || fSumN[bestBin] < 3 ){
      bestBin = -1;
      break;
    }
    
    // calculate the vertex position in the best bin
    Double_t w[6] = { fSum[0][bestBin],
		      fSum[1][bestBin], fSum[2][bestBin],
		      fSum[3][bestBin], fSum[4][bestBin], fSum[5][bestBin] };
    Double_t wI[6];
    {
      wI[0] = w[2]*w[5] - w[4]*w[4];
      wI[1] = w[3]*w[4] - w[1]*w[5];
      wI[2] = w[0]*w[5] - w[3]*w[3];
      wI[3] = w[1]*w[4] - w[2]*w[3];
      wI[4] = w[1]*w[3] - w[0]*w[4];
      wI[5] = w[0]*w[2] - w[1]*w[1];	 
      
      Double_t s = ( w[0]*wI[0] + w[1]*wI[1] + w[3]*wI[3] );
      
      if( s<1.e-10 ){
	bestBin = -1;
	break;
      }
      s = 1./s;	  
      wI[0]*=s;
      wI[1]*=s;
      wI[2]*=s;
      wI[3]*=s;
      wI[4]*=s;
      wI[5]*=s;
    }    

    vtxX += wI[0]*fSum[6][bestBin] + wI[1]*fSum[7][bestBin] + wI[3]*fSum[8][bestBin];
    vtxY += wI[1]*fSum[6][bestBin] + wI[2]*fSum[7][bestBin] + wI[4]*fSum[8][bestBin];
    vtxZ = wI[3]*fSum[6][bestBin] + wI[4]*fSum[7][bestBin] + wI[5]*fSum[8][bestBin];
    //cout<<"SG: "<<iter<<": "<<vtxX<<" "<<vtxY<<" "<<vtxZ<<endl;
  }

  timerReco.Stop();
  
  // Fill output 

  if( bestBin>=0 ){
    double pos[3] = {vtxX, vtxY, vtxZ};
    double s = 400.E-4;
    double cov[6] = {s*s,0,s*s,0,0,s*s};
    AliESDVertex v(pos, cov, 0, fSumN[bestBin]);
    PushBack( (TObject*) &v, kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS,0 );

    //cout<<"ITSVertexerSPD: vtx found "<<vtxX<<" "<<vtxY<<" "<<vtxZ<<endl;
  }

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
