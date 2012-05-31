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

#include "TCollection.h"
#include "TObjArray.h"
#include "AliHLTEveCalo.h"
#include "AliHLTHOMERBlockDesc.h"
#include "TCanvas.h"
#include "AliHLTEveBase.h"
#include "TEveBoxSet.h"
#include "AliPHOSGeometry.h"
#include "TVector3.h"
#include "AliEveHLTEventManager.h"
#include "TEveManager.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "TEveTrans.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TRefArray.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"

ClassImp(AliHLTEveCalo);

AliHLTEveCalo::AliHLTEveCalo(Int_t nm, TString name) : 
  AliHLTEveBase(name), 
  fBoxSetDigits(NULL),
  fBoxSetClusters(NULL),
  fNModules(nm),
  fClustersArray(NULL),
  fName(name), 
  fPadTitles(NULL),
  fInvMassCanvas(NULL),
  fClusterReader(NULL)
{
  // Constructor.

  SetMaxHistograms(9);

  fPadTitles = new TString[GetMaxHistograms()];
 
  for(int i = 0; i < GetMaxHistograms(); i++) {
    fPadTitles[i] = "";
  }

  fClustersArray = new TRefArray();
  fClusterReader = new AliHLTCaloClusterReader();




}

AliHLTEveCalo::~AliHLTEveCalo()
{
  //Destructor
  if(fBoxSetDigits)
    delete fBoxSetDigits;
  fBoxSetDigits = NULL;
  
  if(fBoxSetClusters)
    delete fBoxSetClusters;
  fBoxSetClusters = NULL;


  if(fPadTitles)
    delete [] fPadTitles;
  fPadTitles = NULL;


  if(fClusterReader)
    delete fClusterReader;
  fClusterReader = NULL;

}


void AliHLTEveCalo::ProcessBlock(AliHLTHOMERBlockDesc * block) {
  //See header file for documentation

  if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) { 
    ProcessHistogram(block);
   
  } else {

    
    if ( block->GetDataType().CompareTo("CALOCLUS") == 0 ){
      //cout <<"Skipping calo clusters"<<endl;
      ProcessClusters( block );
    } else if ( block->GetDataType().CompareTo("DIGITTYP") == 0 ) {
      //ProcessDigits( block);
      //
    } else if ( block->GetDataType().CompareTo("CHANNELT") == 0 ) {
      //ProcessClusters( block );
    } else if (!block->GetDataType().CompareTo("ALIESDV0")) {
      ProcessEsdBlock(block);
    }
  }
}

void AliHLTEveCalo::ProcessHistogram(AliHLTHOMERBlockDesc * block ) {
  //See header file for documentation
  
  if(!fCanvas) {
    fCanvas = CreateCanvas(Form("%s QA", fName.Data()), Form("%s QA", fName.Data()));
    fCanvas->Divide(3, 3);
  }

  if(!fInvMassCanvas) {
    fInvMassCanvas = CreateCanvas(Form("%s IM", fName.Data()), Form("%s IM", fName.Data()));
    fInvMassCanvas->Divide(3, 2);
  }


  AddHistogramsToCanvas(block, fCanvas, fHistoCount);


}


// void AliHLTEveCalo::ProcessDigits(AliHLTHOMERBlockDesc* block) {
//   //See header file for documentation
  
//   AliHLTCaloDigitDataStruct *ds = reinterpret_cast<AliHLTCaloDigitDataStruct*> (block->GetData());
//   UInt_t nDigits = block->GetSize()/sizeof(AliHLTCaloDigitDataStruct);
    

//   for(UInt_t i = 0; i < nDigits; i++, ds++) {

//     Float_t x = (ds->fX - 32)* 2.2;
//       Float_t z = (ds->fZ - 28) * 2.2;


//     fBoxSetDigits[4-ds->fModule].AddBox(x, 0, z, 2.2, ds->fEnergy*200, 2.2);
//     fBoxSetDigits[4-ds->fModule].DigitValue(static_cast<Int_t>(ds->fEnergy*10));
//   }

// }

void AliHLTEveCalo::ProcessEsdBlock(AliHLTHOMERBlockDesc * block) {
  AliESDEvent * event = dynamic_cast<AliESDEvent*>(block->GetTObject());
  if (event) {
    ProcessEvent(event);
  } else {
    cout << "problem getting the event!"<<endl;
  }

}

void AliHLTEveCalo::ProcessEvent(AliESDEvent * event) {
  //see header file for documentation



  Int_t nClusters = GetClusters(event, fClustersArray);
  for(int ic = 0; ic < nClusters; ic++) {
    AliESDCaloCluster * cluster = dynamic_cast<AliESDCaloCluster*>(fClustersArray->At(ic));
    ProcessESDCluster(cluster);
  }
  
}


void AliHLTEveCalo::ProcessClusters(AliHLTHOMERBlockDesc* block) {
  //See header file for documentation

  AliHLTCaloClusterHeaderStruct *dh = reinterpret_cast<AliHLTCaloClusterHeaderStruct*> (block->GetData());
  fClusterReader->SetMemory(dh);  

  AliHLTCaloClusterDataStruct * ds;

  while( (ds = fClusterReader->NextCluster()) ){
     AddClusters(ds->fGlobalPos, ds->fModule, ds->fEnergy);
  }

  AliHLTCaloDigitDataStruct *dg = fClusterReader->GetDigits();
  UInt_t nDigits = fClusterReader->GetNDigits();;
  for(UInt_t i = 0; i < nDigits; i++, dg++) {
    AddDigits(dg->fX, dg->fZ, dg->fModule, dg->fEnergy);
  }
}

void AliHLTEveCalo::UpdateElements() {
  //See header file for documentation
  if(fCanvas) fCanvas->Update();
  if(fInvMassCanvas) fInvMassCanvas->Update();


  if(fBoxSetDigits) {
    for(int im = 0; im < fNModules; im++) {
      fBoxSetDigits[im].ElementChanged();
    }
  }

  if(fBoxSetClusters) {
    for(int im = 0; im < fNModules; im++) {
      fBoxSetClusters[im].ElementChanged();
    }
  }

}

void AliHLTEveCalo::ResetElements(){
  //See header file for documentation
  fHistoCount = 0;
  
  if ( fBoxSetDigits ){
    for(int im = 0; im < fNModules; im++){
      fBoxSetDigits[im].Reset();   
    }
  }

  if ( fBoxSetClusters ){
    for(int im = 0; im < fNModules; im++){
      fBoxSetClusters[im].Reset();   
    }
  }

}

Int_t AliHLTEveCalo::GetPadNumber(TString name) {


  for(int i = 0; i < GetMaxHistograms(); i++) {
    if (!fPadTitles[i].CompareTo(name)){
      return i+1;
    }
    else if (!fPadTitles[i].CompareTo("")) {
      //cout <<"in empty title"<<endl;
      fPadTitles[i] = name;
      return i+1;
    }
  }
  
  cout << "AliHLTEveCalo()->GetPadNUmber():returning one"<<endl;
  return 1;

}

void AliHLTEveCalo::DrawInvMassHistogram(TH1F * histo) {
  
  fInvMassCanvas->cd(++fHistoCount);

  histo->SetAxisRange(histo->GetXaxis()->GetBinLowEdge(histo->FindFirstBinAbove(0, 1) - 3), histo->GetXaxis()->GetBinUpEdge(histo->FindLastBinAbove(0, 1) + 3), "X");
  //histo->Fit("gaus", "", "", histo->GetXaxis()->GetBinLowEdge(histo->FindFirstBinAbove(0, 1)), histo->GetXaxis()->GetBinUpEdge(histo->FindLastBinAbove(0, 1)));

  histo->Draw();

}

void AliHLTEveCalo::AddHistogramsToCanvas(AliHLTHOMERBlockDesc * block, TCanvas * canvas, Int_t &/*cdCount*/) {
  //See header file for documentation

  if ( ! block->GetClassName().CompareTo("TObjArray")) {
    TIter next((TObjArray*)(block->GetTObject()));
    TObject *object;


   
    while (( object = (TObject*) next())) {

      TString name = static_cast<TH1*>(object)->GetName();
      if(name.Contains("InvMass")){
	DrawInvMassHistogram(static_cast<TH1F*>(object));
	
      } else {

	
	Int_t iPad = GetPadNumber(name);
	canvas->cd(iPad);
	
	
	//Check if histo is 2D histo
	TH2F* histo2 = dynamic_cast<TH2F*>(object);
	if(histo2){
	  
	  Int_t lb = histo2->FindLastBinAbove(0,1);
	  if(lb > -1) {
	    histo2->SetAxisRange(0, histo2->GetXaxis()->GetBinUpEdge(histo2->FindLastBinAbove(0, 1) + 3), "X");
	    histo2->SetAxisRange(0, histo2->GetYaxis()->GetBinUpEdge(histo2->FindLastBinAbove(0, 2) + 3), "Y");
	  }
	  
	  histo2->Draw("COLZ");
	}
	
	
	//Must be 1D histo
	else {
	  TH1F* histo = dynamic_cast<TH1F*>(object);
	  if (histo) {
	    
	    TString name2 = histo->GetName();
	    
	    if(name2.Contains("Energy")) {
	      histo->SetAxisRange(0, histo->GetXaxis()->GetBinUpEdge(histo->FindLastBinAbove(0, 1) + 3), "X");
	    }
	    
	    
	    histo->Draw();
	  } else {
	    cout <<"AliHLTEveCaloBase::AddHistogramsTocCanvas: Histogram neither TH1F nor TH2F"<<endl;
	  }
	}
      }
    }

  } else if ( ! block->GetClassName().CompareTo("TH1F")) {

    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
    if(histo) {
      
      Int_t iPad = GetPadNumber(histo->GetName());
      canvas->cd(iPad);
      histo->Draw();
    }
  
  
  } else if ( ! block->GetClassName().CompareTo("TH2F")) {
    
    TH2F *histo = reinterpret_cast<TH2F*>(block->GetTObject());
    if(histo) {
      
      Int_t iPad = GetPadNumber(histo->GetName());
      canvas->cd(iPad);
      histo->Draw();
    }


  }
  
  canvas->cd();
}

