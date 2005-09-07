// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group
//-------------------------------------------------------------------------
//          Implementation of the AliHLTTPCHoughBaseTransformer class
//  that is the base class for AliHLTTPCHoughTransformer,
//  AliHLTTPCHoughTransformerVhdl, AliHLTTPCHoughTransformerGlobal,
//  AliHLTTPCHoughTransformerRow    
//-------------------------------------------------------------------------

#include "AliHLTTPCStandardIncludes.h"

#include "AliHLTTPCHoughBaseTransformer.h"

/** \class AliHLTTPCHoughBaseTransformer
<pre>
//_____________________________________________________________
// AliHLTTPCHoughBaseTransformer
//
// The base class for implementations of Hough Transform on ALICE TPC data.
//
// This is an abstract class, and is only meant to provide the interface
// to the different implementations.
//
</pre>
*/

ClassImp(AliHLTTPCHoughBaseTransformer)

AliHLTTPCHoughBaseTransformer::AliHLTTPCHoughBaseTransformer()
{
  //Default constructor
  fDigitRowData = 0;

  fSlice = 0;
  fPatch = 0;
  fLastPatch = -1;
  fLastTransformer = 0;
  fNEtaSegments =0;
  fEtaMin = 0;
  fEtaMax = 0;
  fLowerThreshold = 0;
  fUpperThreshold = 1023;
  fZVertex = 0.0;
}

AliHLTTPCHoughBaseTransformer::AliHLTTPCHoughBaseTransformer(Int_t slice,Int_t patch,Int_t netasegments,Float_t zvertex)
{
  //normal ctor
  fDigitRowData = 0;

  fSlice = 0;
  fPatch = 0;
  fLastPatch = -1;
  fNEtaSegments =0;
  fEtaMin = 0;
  fEtaMax = 0;
  fLowerThreshold = 3;
  fUpperThreshold = 1023;
  fZVertex = zvertex;

  Init(slice,patch,netasegments);
}

AliHLTTPCHoughBaseTransformer::~AliHLTTPCHoughBaseTransformer()
{
  //dtor
}

void AliHLTTPCHoughBaseTransformer::Init(Int_t slice,Int_t patch,Int_t netasegments,Int_t /*n_seqs*/)
{
  //Transformer init
  fSlice = slice;
  fPatch = patch;
  fLastPatch = -1;
  fNEtaSegments = netasegments;
  fEtaMin = 0;
  fEtaMax = fSlice < 18 ? 1. : -1.;
}
