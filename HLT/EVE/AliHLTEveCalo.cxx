/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Svein Lindal <slindal@fys.uio.no   >                  *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTEveCalo.cxx
/// @author Svein Lindal <slindal@fys.uio.no>
/// @brief  Calorimeter base class for the HLT EVE display

#include "AliHLTEveCalo.h"
#include "AliHLTHOMERBlockDesc.h"
#include "TCanvas.h"
#include "AliHLTEveBase.h"
#include "TEveBoxSet.h"
#include "AliPHOSGeometry.h"
#include "TVector3.h"
#include "AliEveHOMERManager.h"
#include "TEveManager.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "TEveTrans.h"
#include "TString.h"

ClassImp(AliHLTEveCalo);

AliHLTEveCalo::AliHLTEveCalo(Int_t nm, TString name) : 
  AliHLTEveBase(), 
  fBoxSet(NULL),
  fElementList(NULL),
  fNModules(nm),
  fName(name)
{
  // Constructor.
}

AliHLTEveCalo::~AliHLTEveCalo()
{
  //Destructor
  if(fBoxSet)
    delete fBoxSet;
  fBoxSet = NULL;

  if(fElementList)
    delete fElementList;
  fElementList = NULL;
}


void AliHLTEveCalo::ProcessBlock(AliHLTHOMERBlockDesc * block) {
  //See header file for documentation

  if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) { 
    ProcessHistogram(block);
   
  } else {

    if( !fElementList ) {
      fElementList = CreateElementList();
      fEventManager->GetEveManager()->AddElement(fElementList);
    }
    
    if ( block->GetDataType().CompareTo("CALOCLUS") == 0 ){
      //cout <<"Skipping calo clusters"<<endl;
      ProcessClusters( block );
    }
    else if ( block->GetDataType().CompareTo("DIGITTYP") == 0 )
      ProcessDigits( block);
    
    else if ( block->GetDataType().CompareTo("CHANNELT") == 0 ) 
      ProcessClusters( block );
  }
}

void AliHLTEveCalo::ProcessHistogram(AliHLTHOMERBlockDesc * block ) {
  //See header file for documentation
  
  if(!fCanvas) {
    fCanvas = CreateCanvas(Form("%s QA", fName.Data()), Form("%s QA", fName.Data()));
    fCanvas->Divide(3, 2);
  }

  AddHistogramsToCanvas(block, fCanvas, fHistoCount);


}


void AliHLTEveCalo::ProcessDigits(AliHLTHOMERBlockDesc* block) {
  //See header file for documentation
  
  AliHLTCaloDigitDataStruct *ds = reinterpret_cast<AliHLTCaloDigitDataStruct*> (block->GetData());
  UInt_t nDigits = block->GetSize()/sizeof(AliHLTCaloDigitDataStruct);
    

  for(UInt_t i = 0; i < nDigits; i++, ds++) {

    Float_t x = (ds->fX - 32)* 2.2;
      Float_t z = (ds->fZ - 28) * 2.2;

      cout << "MODULE DIGITTYP  :" << ds->fModule;

    fBoxSet[4-ds->fModule].AddBox(x, 0, z, 2.2, ds->fEnergy*200, 2.2);
    fBoxSet[4-ds->fModule].DigitValue(static_cast<Int_t>(ds->fEnergy*10));
  }

}


void AliHLTEveCalo::ProcessClusters(AliHLTHOMERBlockDesc* block) {
  //See header file for documentation


  AliHLTCaloClusterHeaderStruct *dh = reinterpret_cast<AliHLTCaloClusterHeaderStruct*> (block->GetData());
  AliHLTCaloClusterReader * clusterReader = new AliHLTCaloClusterReader();
  clusterReader->SetMemory(dh);  

  AliHLTCaloClusterDataStruct * ds;


  
  while( (ds = clusterReader->NextCluster()) ){
    //    AddClusters(ds->fGlobalPos, ds->fModule, ds->fEnergy);
  }

  AliHLTCaloDigitDataStruct *dg = clusterReader->GetDigits();
  UInt_t nDigits = clusterReader->GetNDigits();;
  for(UInt_t i = 0; i < nDigits; i++, dg++) {
    AddDigits(dg->fX, dg->fZ, dg->fModule, dg->fEnergy);
  }
}

void AliHLTEveCalo::UpdateElements() {
  //See header file for documentation
  if(fCanvas) fCanvas->Update();

  if(fBoxSet) {
    for(int im = 0; im < fNModules; im++) {
      fBoxSet[im].ElementChanged();
    }
  }
}

void AliHLTEveCalo::ResetElements(){
  //See header file for documentation
  fHistoCount = 0;
  
  if ( fBoxSet ){
    for(int im = 0; im < fNModules; im++){
      cout<<"Resetting"<<endl;
      fBoxSet[im].Reset();   
    }
  }
}
