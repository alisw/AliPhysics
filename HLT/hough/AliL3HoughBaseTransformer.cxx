// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3MemHandler.h"
#include "AliL3DigitData.h"
#include "AliL3Histogram.h"
#include "AliL3HoughBaseTransformer.h"

/** \class AliL3HoughBaseTransformer
<pre>
//_____________________________________________________________
// AliL3HoughBaseTransformer
//
// The base class for implementations of Hough Transform on ALICE TPC data.
//
// This is an abstract class, and is only meant to provide the interface
// to the different implementations.
//
</pre>
*/

ClassImp(AliL3HoughBaseTransformer)

AliL3HoughBaseTransformer::AliL3HoughBaseTransformer()
{
  //Default constructor
  fDigitRowData = 0;

  fSlice = 0;
  fPatch = 0;
  fNEtaSegments =0;
  fEtaMin = 0;
  fEtaMax = 0;
  fLowerThreshold = 0;
  fUpperThreshold = 1023;
  fZVertex = 0.0;
}

AliL3HoughBaseTransformer::AliL3HoughBaseTransformer(Int_t slice,Int_t patch,Int_t netasegments,Float_t zvertex)
{
  //normal ctor
  fDigitRowData = 0;

  fSlice = 0;
  fPatch = 0;
  fNEtaSegments =0;
  fEtaMin = 0;
  fEtaMax = 0;
  fLowerThreshold = 3;
  fUpperThreshold = 1023;
  fZVertex = zvertex;

  Init(slice,patch,netasegments);
}

AliL3HoughBaseTransformer::~AliL3HoughBaseTransformer()
{
  //dtor
}

void AliL3HoughBaseTransformer::Init(Int_t slice,Int_t patch,Int_t netasegments,Int_t /*n_seqs*/)
{
  //Transformer init
  fSlice = slice;
  fPatch = patch;
  fNEtaSegments = netasegments;
  fEtaMin = 0;
  fEtaMax = fSlice < 18 ? 1. : -1.;
}
