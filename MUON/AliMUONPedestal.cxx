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
fNEvents(0),
fRunNumber(0),
fNChannel(0),
fNManu(0),
fErrorBuspatchTable(new THashTable(100,2)),
fDate(new TTimeStamp()),
fFilcout(0),
fPedestalStore(new AliMUON2DMap(kFALSE)),
fIndex(-1)
{
/// Default constructor

//   sprintf(fOutFolder,".");
  sprintf(fHistoFileName,"");
  sprintf(fprefixDA,""); 
}
//  AliMUONPedestal& operator=(const AliMUONPedestal& other); Copy ctor

//______________________________________________________________________________
AliMUONPedestal::~AliMUONPedestal()
{
/// Destructor
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
}

//______________________________________________________________________________
TString AliMUONPedestal::WritePedHeader(void) 
{
///

  ostringstream stream;
  stream<<"//===========================================================================" << endl;
  stream<<"//                       Pedestal file calculated by MUONTRKda"<<endl;
  stream<<"//===========================================================================" << endl;
  stream<<"//       * Run           : " << fRunNumber << endl; 
  stream<<"//       * Date          : " << fDate->AsString("l") <<endl;
  stream<<"//       * Statictics    : " << fNEvents << endl;
  stream<<"//       * # of MANUS    : " << fNManu << endl;
  stream<<"//       * # of channels : " << fNChannel << endl;
  if (fErrorBuspatchTable->GetSize())
  {
    stream<<"//"<<endl;
    stream<<"//       * Buspatches with less statistics (due to parity errors)"<<endl;
    TIterator* iter = fErrorBuspatchTable->MakeIterator();
    AliMUONErrorCounter* parityerror;
    while((parityerror = (AliMUONErrorCounter*) iter->Next()))
    {
      stream<<"//         bp "<<parityerror->BusPatch()<<" events used "<<fNEvents-parityerror->Events()<<endl;
    }
  }  
  stream<<"//"<<endl;
  stream<<"//---------------------------------------------------------------------------" << endl;
  stream<<"//---------------------------------------------------------------------------" << endl;
  stream<<"//      BP     MANU     CH.      MEAN    SIGMA"<<endl;
  stream<<"//---------------------------------------------------------------------------" << endl;

  return TString(stream.str().c_str());
}

//______________________________________________________________________________
TString AliMUONPedestal::WritePedData(Int_t BP, Int_t Manu, Int_t ch, Double_t pedMean, Double_t pedSigma) 
{
///

  ostringstream stream("");
  stream << "\t" << BP << "\t" << Manu <<"\t"<< ch << "\t"
         << pedMean <<"\t"<< pedSigma << endl;
  return TString(stream.str().c_str());

}

//______________________________________________________________________________
void AliMUONPedestal::MakePedStore(TString shuttleFile_1 = "") 
{

  /// Store pedestals in ASCII files
  Double_t pedMean;
  Double_t pedSigma;
  ofstream fileout;
#ifdef ALI_AMORE
  ostringstream stringout; // String to be sent to AMORE_DB
#endif
  TString tempstring;  
  Int_t busPatchId;
  Int_t manuId;
  Int_t channelId;

// histo
  TFile*  histoFile = 0;
  TTree* tree = 0;
  TH1F* pedMeanHisto = 0;
  TH1F* pedSigmaHisto = 0;

  if (fIndex<0) // Pedestal run (fIndex=-1)
  {
    sprintf(fHistoFileName,"%s_%d.root",fprefixDA,fRunNumber);
    histoFile = new TFile(fHistoFileName,"RECREATE","MUON Tracking pedestals");

    Char_t name[255];
    Char_t title[255];
    sprintf(name,"pedmean_allch");
    sprintf(title,"Pedestal mean all channels");
    Int_t nx = 4096;
    Int_t xmin = 0;
    Int_t xmax = 4095; 
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
  }

  if (!shuttleFile_1.IsNull()) {
    fileout.open(shuttleFile_1.Data());
    tempstring = WritePedHeader();
    fileout << tempstring;
#ifdef ALI_AMORE
    stringout << tempstring;
#endif
  }
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

    for (channelId = 0; channelId < ped->Size() ; ++channelId) {
      pedMean  = ped->ValueAsDouble(channelId, 0);

      if (pedMean > 0) { // connected channels

        ped->SetValueAsDouble(channelId, 0, pedMean/(Double_t)eventCounter);

        pedMean  = ped->ValueAsDouble(channelId, 0);
        pedSigma = ped->ValueAsDouble(channelId, 1);

        ped->SetValueAsDouble(channelId, 1, TMath::Sqrt(TMath::Abs(pedSigma/(Double_t)eventCounter - pedMean*pedMean)));

        pedMean  = ped->ValueAsDouble(channelId, 0);
        pedSigma = ped->ValueAsDouble(channelId, 1);


        if (!shuttleFile_1.IsNull()) {
          tempstring = WritePedData(busPatchId,manuId,channelId,pedMean,pedSigma);
          fileout << tempstring;
#ifdef ALI_AMORE
          stringout << tempstring;
#endif
        }
        if(fIndex<0) // Pedestal Run
        {
          pedMeanHisto->Fill(pedMean);
          pedSigmaHisto->Fill(pedSigma);
          tree->Fill();
        }
      }
    }
  }

// file outputs
  if (!shuttleFile_1.IsNull())  fileout.close();

// Outputs to root file and AMORE DB
#ifdef ALI_AMORE
  //
  //Send objects to the AMORE DB
  //
  const char *role=gSystem->Getenv("AMORE_DA_NAME");
  if ( role ){
    amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
//  TObjString peddata(stringout.str().c_str());
    TObjString peddata(stringout.str().c_str());
    Int_t status =0;
    status = amoreDA.Send("Pedestals",&peddata);
    if ( status )
      cout << "Warning: Failed to write Pedestals in the AMORE database : " << status << endl;
  } 
  else {
    cout << "Warning: environment variable 'AMORE_DA_NAME' not set. Cannot write to the AMORE database" << endl;
  }
#endif
  if(fIndex<0) // Pedestal Run
  { 
    histoFile->Write();  
    histoFile->Close(); 
    delete fPedestalStore ;
  }

}
