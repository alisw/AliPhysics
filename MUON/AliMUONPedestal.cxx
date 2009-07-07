/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    * 
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliMUONPedestal.h"
#include "AliMUONErrorCounter.h"
#include "AliMUONVStore.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"

#include <TString.h>
#include <THashTable.h>
#include <THashList.h>
#include <TTimeStamp.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <Riostream.h>

#include <sstream>

//-----------------------------------------------------------------------------
/// \class AliMUONPedestal
///
/// Implementation of the pedestal computing
///
/// add
/// 
///
/// \author Alberto Baldisseri, JL Charvet (05/05/2009)
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONPedestal)
/// \endcond

//______________________________________________________________________________
AliMUONPedestal::AliMUONPedestal()
: TObject(),
fN(0),
fNEvents(0),
fRunNumber(0),
fNChannel(0),
fNManu(0),
fErrorBuspatchTable(new THashTable(100,2)),
fManuBuspatchTable(new THashTable(100,2)),
fDate(new TTimeStamp()),
fFilcout(0),
fPedestalStore(new AliMUON2DMap(kFALSE)),
fIndex(-1)
{
/// Default constructor
  sprintf(fHistoFileName," ");
  sprintf(fprefixDA," "); 
}
//  AliMUONPedestal& operator=(const AliMUONPedestal& other); Copy ctor

//______________________________________________________________________________
AliMUONPedestal::~AliMUONPedestal()
{
/// Destructor
  delete fPedestalStore;
}

//______________________________________________________________________________
void AliMUONPedestal::MakePed(Int_t busPatchId, Int_t manuId, Int_t channelId, Int_t charge)
{
  /// Compute pedestals values
  AliMUONVCalibParam* ped = 
    static_cast<AliMUONVCalibParam*>(fPedestalStore ->FindObject(busPatchId, manuId));

  if (!ped) {
    fNManu++;
    ped = new AliMUONCalibParamND(2, kNChannels,busPatchId, manuId, -1.); // put default wise -1, not connected channel
    fPedestalStore ->Add(ped);  
  }

  // Initialization for the first value
  if (ped->ValueAsDouble(channelId, 0) == -1) ped->SetValueAsDouble(channelId, 0, 0.);
  if (ped->ValueAsDouble(channelId, 1) == -1) ped->SetValueAsDouble(channelId, 1, 0.);

  Double_t pedMean  = ped->ValueAsDouble(channelId, 0) + (Double_t) charge;
  Double_t pedSigma = ped->ValueAsDouble(channelId, 1) + (Double_t) charge*charge;

  ped->SetValueAsDouble(channelId, 0, pedMean);
  ped->SetValueAsDouble(channelId, 1, pedSigma);

  char bpmanuname[256];
  AliMUONErrorCounter* manuCounter;
  sprintf(bpmanuname,"bp%dmanu%d",busPatchId,manuId);		
  if (!(manuCounter = (AliMUONErrorCounter*)GetManuBuspatchTable()->FindObject(bpmanuname)))
    {
      // New (buspatch,manu)
      manuCounter = new AliMUONErrorCounter(busPatchId,manuId);
      manuCounter->SetName(bpmanuname);
      GetManuBuspatchTable()->Add(manuCounter);
    }
  else
    {
      // Existing buspatch
      manuCounter->Increment();
    }	
}

//______________________________________________________________________________
void AliMUONPedestal::Finalize()
{
  Double_t pedMean;
  Double_t pedSigma;
  Int_t busPatchId;
  Int_t manuId;
  Int_t channelId;

  // print in logfile
  if (fErrorBuspatchTable->GetSize())
  {
    cout<<"\n* Buspatches with less statistics (due to parity errors)"<<endl;
    (*fFilcout)<<"\n* Buspatches with less statistics (due to parity errors)"<<endl;
    TIterator* iter = fErrorBuspatchTable->MakeIterator();
    AliMUONErrorCounter* parityerror;
    while((parityerror = (AliMUONErrorCounter*) iter->Next()))
    {
      cout<<"  bp "<<parityerror->BusPatch()<<": events used = "<<fNEvents-parityerror->Events()<<endl;
      (*fFilcout)<<"  bp "<<parityerror->BusPatch()<<": events used = "<<fNEvents-parityerror->Events()<<endl;
    }

  }

// iterator over pedestal
  TIter next(fPedestalStore ->CreateIterator());
  AliMUONVCalibParam* ped;

  while ( ( ped = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
  {
    busPatchId              = ped->ID0();
    manuId                  = ped->ID1();
    if(manuId==0)
      {
	cout << " !!! BIG WARNING: ManuId = " << manuId << " !!! in  BP = " << busPatchId << endl;
	(*fFilcout) << " !!! BIG WARNING: ManuId = " << manuId << " !!! in  BP = " << busPatchId << endl;
      }
    Int_t eventCounter;
    // Correct the number of events for buspatch with errors
    char bpname[256];
    AliMUONErrorCounter* errorCounter;
    sprintf(bpname,"bp%d",busPatchId);
    if ((errorCounter = (AliMUONErrorCounter*)fErrorBuspatchTable->FindObject(bpname)))
    {
      eventCounter = fNEvents - errorCounter->Events();
    }
    else
    {
      eventCounter = fNEvents;
    }

    Int_t occupancy;
    // value of (buspatch, manu) occupancy
    char bpmanuname[256];
    AliMUONErrorCounter* manuCounter;
    sprintf(bpmanuname,"bp%dmanu%d",busPatchId,manuId);
    manuCounter = (AliMUONErrorCounter*)fManuBuspatchTable->FindObject(bpmanuname);
    occupancy = manuCounter->Events()/64/eventCounter;
    if(occupancy>1)
    {
 	cout << " !!! BIG WARNING: ManuId = " << manuId << " !!! in  BP = " << busPatchId << " occupancy (>1) = " << occupancy << endl;
	(*fFilcout) << " !!! BIG WARNING: ManuId = " << manuId << " !!! in  BP = " << busPatchId << " occupancy (>1) = " << occupancy <<endl;
    }

    for (channelId = 0; channelId < ped->Size() ; ++channelId) 
    {
      pedMean  = ped->ValueAsDouble(channelId, 0);

      if (pedMean > 0) // connected channels
      {
	  ped->SetValueAsDouble(channelId, 0, pedMean/(Double_t)eventCounter);
	  pedMean  = ped->ValueAsDouble(channelId, 0);
	  pedSigma = ped->ValueAsDouble(channelId, 1);
	  ped->SetValueAsDouble(channelId, 1, TMath::Sqrt(TMath::Abs(pedSigma/(Double_t)eventCounter - pedMean*pedMean)));
	  if(manuId == 0 || occupancy>1)
	  {
	      ped->SetValueAsDouble(channelId, 0, kADCMax);
	      ped->SetValueAsDouble(channelId, 1, kADCMax);
	  }
      }
      else
      {
	  ped->SetValueAsDouble(channelId, 0, kADCMax);
	  ped->SetValueAsDouble(channelId, 1, kADCMax);
      }
    }
  }
}

//______________________________________________________________________________
void AliMUONPedestal::MakeASCIIoutput(ostream& out) const
{
/// put pedestal store in the output stream

  out<<"//===========================================================================" << endl;
  out<<"//                       Pedestal file calculated by MUONTRKda"<<endl;
  out<<"//===========================================================================" << endl;
  out<<"//       * Run           : " << fRunNumber << endl; 
  out<<"//       * Date          : " << fDate->AsString("l") <<endl;
  out<<"//       * Statictics    : " << fNEvents << endl;
  out<<"//       * # of MANUS    : " << fNManu << endl;
  out<<"//       * # of channels : " << fNChannel << endl;
  if (fErrorBuspatchTable->GetSize())
  {
    out<<"//"<<endl;
    out<<"//       * Buspatches with less statistics (due to parity errors)"<<endl;
    TIterator* iter = fErrorBuspatchTable->MakeIterator();
    AliMUONErrorCounter* parityerror;
    while((parityerror = (AliMUONErrorCounter*) iter->Next()))
    {
      out<<"//         bp "<<parityerror->BusPatch()<<" events used "<<fNEvents-parityerror->Events()<<endl;
    }
  }  
  out<<"//"<<endl;
  out<<"//---------------------------------------------------------------------------" << endl;
  out<<"//---------------------------------------------------------------------------" << endl;
  out<<"//      BP     MANU     CH.      MEAN    SIGMA"<<endl;
  out<<"//---------------------------------------------------------------------------" << endl;

  // iterator over pedestal
  TIter next(fPedestalStore ->CreateIterator());
  AliMUONVCalibParam* ped;
  
  while ( ( ped = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
  {
    Int_t busPatchId = ped->ID0();
    Int_t manuId = ped->ID1();
    
    for ( Int_t channelId = 0; channelId < ped->Size(); ++channelId ) 
    {
      Double_t pedMean  = ped->ValueAsDouble(channelId, 0);
      Double_t pedSigma = ped->ValueAsDouble(channelId, 1);
        
      out << "\t" << busPatchId << "\t" << manuId <<"\t"<< channelId << "\t" << pedMean <<"\t"<< pedSigma << endl;        
    }
  }
}

//______________________________________________________________________________
void AliMUONPedestal::MakeControlHistos()
{

  if (fIndex>=0) return; // Pedestal run (fIndex=-1)

  Double_t pedMean;
  Double_t pedSigma;
  Int_t busPatchId;
  Int_t manuId;
  Int_t channelId;

// histo
  TFile*  histoFile = 0;
  TTree* tree = 0;
  TH1F* pedMeanHisto = 0;
  TH1F* pedSigmaHisto = 0;
    
  sprintf(fHistoFileName,"%s.root",fprefixDA);
  histoFile = new TFile(fHistoFileName,"RECREATE","MUON Tracking pedestals");

  Char_t name[255];
  Char_t title[255];
  sprintf(name,"pedmean_allch");
  sprintf(title,"Pedestal mean all channels");
  Int_t nx = kADCMax+1;
  Int_t xmin = 0;
  Int_t xmax = kADCMax; 
  pedMeanHisto = new TH1F(name,title,nx,xmin,xmax);
  pedMeanHisto->SetDirectory(histoFile);

  sprintf(name,"pedsigma_allch");
  sprintf(title,"Pedestal sigma all channels");
  nx = 201;
  xmin = 0;
  xmax = 200; 
  pedSigmaHisto = new TH1F(name,title,nx,xmin,xmax);
  pedSigmaHisto->SetDirectory(histoFile);

  tree = new TTree("t","Pedestal tree");
  tree->Branch("bp",&busPatchId,"bp/I");
  tree->Branch("manu",&manuId,",manu/I");
  tree->Branch("channel",&channelId,",channel/I");
  tree->Branch("pedMean",&pedMean,",pedMean/D");
  tree->Branch("pedSigma",&pedSigma,",pedSigma/D");

  // iterator over pedestal
  TIter next(fPedestalStore ->CreateIterator());
  AliMUONVCalibParam* ped;
  
  while ( ( ped = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
  {
    busPatchId = ped->ID0();
    manuId = ped->ID1();
    
    for ( channelId = 0; channelId < ped->Size(); ++channelId ) 
    {
      pedMean  = ped->ValueAsDouble(channelId, 0);
      pedSigma = ped->ValueAsDouble(channelId, 1);
          
      pedMeanHisto->Fill(pedMean);
      pedSigmaHisto->Fill(pedSigma);
      tree->Fill();  
    }
  }
    
  histoFile->Write();  
  histoFile->Close(); 

}
