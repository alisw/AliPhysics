// Author: Svein Lindal <slindal@fys.uio.no>

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliHLTEveAny.h"
#include "AliHLTHOMERBlockDesc.h"
#include "TCanvas.h"

ClassImp(AliHLTEveAny);

AliHLTEveAny::AliHLTEveAny() : 
  AliHLTEveBase()
{
  // Constructor.
  SetDetector("Any detector");
}

AliHLTEveAny::~AliHLTEveAny()
{
  //Destructor, not implemented
}


void AliHLTEveAny::ProcessBlock(AliHLTHOMERBlockDesc * block) {
  //See header file for documentation

  if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) { 
    ProcessHistogram(block);
    
  } else {
    cout << "block of datatype " << block->GetDataType() << " with no parent"<<endl;
    
  }
}

void AliHLTEveAny::ProcessHistogram(AliHLTHOMERBlockDesc * block ) {
  if(!fCanvas) {
    fCanvas = CreateCanvas("Orphans", "Orphans");
    fCanvas->Divide(3, 2);
    SetMaxHistograms(6);
  }
  AddHistogramsToCanvas(block, fCanvas, fHistoCount);
}

void AliHLTEveAny::UpdateElements() {
  if(fCanvas) fCanvas->Update();
}

void AliHLTEveAny::ResetElements(){
  fHistoCount = 0;
}
