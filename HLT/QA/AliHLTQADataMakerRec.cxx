// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Zhong-Bao Yin <Zhong-Bao.Yin@cern.ch>                 *
//*                  Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   AliHLTQADataMakerRec.cxx
    @author Zhongbao Yin, Matthias Richter
    @date   2009-05-14
    @brief  Container for the HLT offline QA
*/
#include "AliHLTQADataMakerRec.h"
#include "AliESDEvent.h"
#include <iostream>

#include "TH1F.h"
#include "TH2F.h"

#include "AliESDtrack.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTQADataMakerRec)

AliHLTQADataMakerRec::AliHLTQADataMakerRec()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTQADataMakerRec::~AliHLTQADataMakerRec()
{
  // see header file for class documentation
}

void AliHLTQADataMakerRec::Exec(AliQAv1::TASKINDEX_t task, TObject * data) 
{ 
  // special handling for esds
  if ( task == AliQAv1::kESDS ) {
    AliESDEvent * esd = NULL;
    AliESDEvent * hltesd = NULL;
    if (data->IsA() == AliESDEvent::Class()) {
      // silently skip this. Currently HLT QA is still called as
      // part of AliQAManager::RunOneEvent with the esd
      return;
    }
    if (data->InheritsFrom("TObjArray")) {
      TObjArray* array=dynamic_cast<TObjArray*>(data);
      if (array && array->GetEntriesFast()>0) {
	esd = dynamic_cast<AliESDEvent *>(array->At(0)) ;
      }
      if (array && array->GetEntriesFast()>1) {
	hltesd = dynamic_cast<AliESDEvent *>(array->At(1)) ;
      }
    } else {
      esd = static_cast<AliESDEvent *>(data) ; 
    }

    if (esd && strcmp(esd->ClassName(), "AliESDEvent") == 0) {
      if (hltesd) {
	MakeESDs(esd, hltesd);
      } else {
	AliError(Form("HLT ESD missing or wrong class type (%p), skipping HLT QA task kESDs", hltesd));
      }
    } else {
      AliError(Form("ESD missing or wrong class type (%p), skipping HLT QA task kESDSs", esd));
    }
  } else {
    // forward for all other types
    AliQADataMakerRec::Exec(task, data);
  }
}

void AliHLTQADataMakerRec::StartOfDetectorCycle()
{
  // see header file for class documentation
}

void AliHLTQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray** /*list*/)
{
  // see header file for class documentation
  
  if(GetESDsData(kPHLT)->GetEntries()){
    GetESDsData(kPRatio)->Sumw2();
    GetESDsData(kPRatio)->Add(GetESDsData(kPOffline));
    GetESDsData(kPRatio)->Divide(GetESDsData(kPHLT));
  }
  
  if(GetESDsData(kPHLTFired)->GetEntries()){ 
    GetESDsData(kPRatioFired)->Sumw2(); 
    GetESDsData(kPRatioFired)->Add(GetESDsData(kPOfflineFired)); 
    GetESDsData(kPRatioFired)->Divide(GetESDsData(kPHLTFired));
  } 
  

}

void AliHLTQADataMakerRec::MakeRaws(AliRawReader * rawReader)
{
  // see header file for class documentation
  if (!rawReader) return;
}

void AliHLTQADataMakerRec::InitESDs(){
  
  //create ESDs histograms in ESDs subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH2F * histESDMultiplicity =
    new TH2F("hESDMultiplicity", 
	     "Number of tracks from all events; Number of offline tracks; Number of HLT tracks", 
	     300, 0, 300, 300, 0, 300);
  histESDMultiplicity->Sumw2();
  Add2ESDsList(histESDMultiplicity, kMultiplicity, !expert, image);

  TH2F * histESDMultiplicityFired = 
    new TH2F("hESDMultiplicityFired",  
             "Number of tracks from HLT triggered events; Number of offline tracks; Number of HLT tracks",  
             300, 0, 300, 300, 0, 300);
  histESDMultiplicityFired->Sumw2();
  Add2ESDsList(histESDMultiplicityFired, kMultiplicityFired, !expert, image);

  TH2F * histESDNCls = 
    new TH2F("hESDNCls", "Mean number of TPC clusters from all events; Offline; HLT",
	     200, 0, 200, 200, 0, 200);
  histESDNCls->Sumw2();
  Add2ESDsList(histESDNCls, kNCls, !expert, image);
  
  TH2F * histESDNClsFired =  
    new TH2F("hESDNClsFired", "Mean number of TPC clusters from triggered events; Offline; HLT", 
             200, 0, 200, 200, 0, 200); 
  histESDNClsFired->Sumw2();
  Add2ESDsList(histESDNClsFired, kNClsFired, !expert, image);

  TH1F * histPHLT 
    = new TH1F("hPHLT", "P distribution for all events from HLT; P [GeV/c]", 
	       100, -0.5, 99.5);
  histPHLT->Sumw2();
  Add2ESDsList(histPHLT, kPHLT, !expert, image);

  TH1F * histPOffline 
    = new TH1F("hPOffline", 
	       "P distribution for all events from offline; P [GeV/c]",  
	       100, -0.5, 99.5);
  histPOffline->Sumw2();
  Add2ESDsList(histPOffline, kPOffline, !expert, image);
  
  TH1F * histPHLTFired 
    = new TH1F("hPHLTFired", 
	       "P distribution for fired events from HLT; P [GeV/c]",  
	       100, -0.5, 99.5);
  histPHLTFired->Sumw2();
  Add2ESDsList(histPHLTFired, kPHLTFired, !expert, image); 
 
  TH1F * histPOfflineFired 
    = new TH1F("hPOfflineFired",  
	       "P distribution for fired events from offline; P [GeV/c]",   
	       100, -0.5, 99.5);
  histPOfflineFired->Sumw2();
  Add2ESDsList(histPOfflineFired, kPOfflineFired, !expert, image);

  TH1F * histPRatio 
    = new TH1F("hPRatio", 
	       "Ratio of P distribution for all events; P [GeV/c]",
	       100, -0.5, 99.5);
  histPRatio->Sumw2();
  Add2ESDsList(histPRatio, kPRatio, !expert, image);
  
  TH1F * histPRatioFired  
    = new TH1F("hPRatioFired",  
               "Ratio of P distribution for fired events; P [GeV/c]",
	       100, -0.5, 99.5);
  histPRatioFired->Sumw2();
  Add2ESDsList(histPRatioFired, kPRatioFired, !expert, image);

  TH1F * histPtHLT  
    = new TH1F("hPtHLT", 
	       "P_{T} distribution for all events from HLT; P_{T} [GeV/c]",  
               200, -100.5, 99.5);
  histPtHLT->Sumw2();
  Add2ESDsList(histPtHLT, kPtHLT, !expert, image); 
 
  TH1F * histPtOffline  
    = new TH1F("hPtOffline",  
               "P_{T} distribution for all events from offline; P_{T} [GeV/c]",   
               200, -100.5, 99.5);
  histPtOffline->Sumw2();
  Add2ESDsList(histPtOffline, kPtOffline, !expert, image);

  TH1F * histPtHLTFired   
    = new TH1F("hPtHLTFired",  
               "P_{T} distribution for fired events from HLT; P_{T} [GeV/c]",   
               200, -100.5, 99.5);
  histPtHLTFired->Sumw2();
  Add2ESDsList(histPtHLTFired, kPtHLTFired, !expert, image);  
  
  TH1F * histPtOfflineFired   
    = new TH1F("hPtOfflineFired",   
               "P_{T} distribution for fired events from offline; P_{T} [GeV/c]", 
	       200, -100.5, 99.5);
  histPtOfflineFired->Sumw2();
  Add2ESDsList(histPtOfflineFired, kPtOfflineFired, !expert, image); 

  TH1F * histNClsPerTrkHLT
    = new TH1F("hNClsPerTrkHLT", "Clusters per HLT track; N cluster; Counts",
	       200, 0, 200);
  histNClsPerTrkHLT->Sumw2();
  Add2ESDsList(histNClsPerTrkHLT, kNClsPerTrkHLT, !expert, image);

  TH1F * histNClsPerTrkOffline 
    = new TH1F("hNClsPerTrkOffline", 
	       "Clusters per offline track; N cluster; Counts", 
               200, 0, 200); 
  histNClsPerTrkOffline->Sumw2();
  Add2ESDsList(histNClsPerTrkOffline, kNClsPerTrkOffline, !expert, image);

TH1F * histNClsPerTrkHLTFired 
  = new TH1F("hNClsPerTrkHLTFired", 
	     "Clusters per HLT track from HLT fired events; N cluster; Counts", 
	     200, 0, 200); 
 histNClsPerTrkHLTFired->Sumw2();
 Add2ESDsList(histNClsPerTrkHLTFired, kNClsPerTrkHLTFired, !expert, image); 
 
  TH1F * histNClsPerTrkOfflineFired  
    = new TH1F("hNClsPerTrkOfflineFired",  
               "Clusters per offline track from HLT fired events; N cluster; Counts",  
               200, 0, 200);  
  histNClsPerTrkOfflineFired->Sumw2();
  Add2ESDsList(histNClsPerTrkOfflineFired, kNClsPerTrkOfflineFired, !expert, image);

  TH1F * histPhiHLT = 
    new TH1F("hPhiHLT", "Phi distribution of HLT tracks; Phi; Counts", 
	     360, 0, 360);
  histPhiHLT->Sumw2();
  Add2ESDsList(histPhiHLT, kPhiHLT, !expert, image);

  TH1F * histPhiOffline =  
    new TH1F("hPhiOffline", 
	     "Phi distribution of offline tracks; Phi; Counts",  
             360, 0, 360); 
  histPhiOffline->Sumw2(); 
  Add2ESDsList(histPhiOffline, kPhiOffline, !expert, image);

  TH1F * histPhiHLTFired =  
    new TH1F("hPhiHLTFired", "Phi distribution of HLT tracks from HLT fired event ; Phi; Counts",  
             360, 0, 360); 
  histPhiHLTFired->Sumw2(); 
  Add2ESDsList(histPhiHLTFired, kPhiHLTFired, !expert, image); 
 
  TH1F * histPhiOfflineFired =   
    new TH1F("hPhiOfflineFired",  
             "Phi distribution of offline tracks from HLT fired events; Phi; Counts",   
             360, 0, 360);  
  histPhiOfflineFired->Sumw2();  
  Add2ESDsList(histPhiOfflineFired, kPhiOfflineFired, !expert, image); 

  TH1F * histEtaHLT =
    new TH1F("hEtaHLT", "Eta distribution of HLT tracks; Eta; Counts",
	     200, -1, 1);
  histEtaHLT->Sumw2();
  Add2ESDsList(histEtaHLT, kEtaHLT, !expert, image);

  TH1F * histEtaOffline = 
    new TH1F("hEtaOffline", "Eta distribution of offline tracks; Eta; Counts", 
             200, -1, 1); 
  histEtaHLT->Sumw2(); 
  Add2ESDsList(histEtaOffline, kEtaOffline, !expert, image);

  TH1F * histEtaHLTFired = 
    new TH1F("hEtaHLTFired", 
	     "Eta distribution of HLT tracks from HLT fired events; Eta; Counts", 
             200, -1, 1); 
  histEtaHLTFired->Sumw2(); 
  Add2ESDsList(histEtaHLTFired, kEtaHLTFired, !expert, image); 
 
  TH1F * histEtaOfflineFired =  
    new TH1F("hEtaOfflineFired", 
	     "Eta distribution of offline tracks from HLT fired events; Eta; Counts",  
             200, -1, 1);  
  histEtaHLTFired->Sumw2();  
  Add2ESDsList(histEtaOfflineFired, kEtaOfflineFired, !expert, image);

}

void AliHLTQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // see header file for class documentation
  
  // as an extension in the QA interface also the hlt esd can be sent
  // in order to preserve backward compatibility, a new function has been
  // introduced.
  //
  // NOTE: This function is not the place for HLT QA
  if (!esd) return;
}

void AliHLTQADataMakerRec::MakeESDs(AliESDEvent * esd, AliESDEvent* hltesd)
{
  // HLT QA on ESDs
  if (!esd || !hltesd) {
    AliError("invalid parameter: missing ESDs");
    return;
  }

  // make QA data from ESDs
 
  const Int_t nESDTracks = esd->GetNumberOfTracks();
  const Int_t nHLTesdTracks = hltesd->GetNumberOfTracks();
  GetESDsData(kMultiplicity)->Fill(nESDTracks, nHLTesdTracks);

  Int_t nClsHLT = 0;
  Int_t nClsOffline = 0;
  
  for(Int_t i = 0; i < nESDTracks; i++){
    AliESDtrack *  esdTrk = esd->GetTrack(i);
    GetESDsData(kPOffline)->Fill(esdTrk->P());
    GetESDsData(kPtOffline)->Fill(esdTrk->GetSignedPt());
    GetESDsData(kNClsPerTrkOffline)->Fill(esdTrk->GetTPCNcls());
    GetESDsData(kPhiOffline)->Fill(esdTrk->Phi()*TMath::RadToDeg());
    GetESDsData(kEtaOffline)->Fill(esdTrk->Eta());
    nClsOffline+=esdTrk->GetTPCNcls();
  }

  for(Int_t i = 0; i < nHLTesdTracks; i++){ 
    AliESDtrack *  hltEsdTrk = hltesd->GetTrack(i); 
    GetESDsData(kPHLT)->Fill(hltEsdTrk->P()); 
    GetESDsData(kPtHLT)->Fill(hltEsdTrk->GetSignedPt());
    GetESDsData(kNClsPerTrkHLT)->Fill(hltEsdTrk->GetTPCNcls());
    GetESDsData(kPhiHLT)->Fill(hltEsdTrk->Phi()*TMath::RadToDeg());
    GetESDsData(kEtaHLT)->Fill(hltEsdTrk->Eta());
    nClsHLT += hltEsdTrk->GetTPCNcls();
  } 
  
  if(nESDTracks)
    nClsOffline /= Float_t(nESDTracks);
  if(nHLTesdTracks)
    nClsHLT /= Float_t(nHLTesdTracks);

  GetESDsData(kNCls)->Fill(nClsOffline, nClsHLT);


  if(hltesd->IsHLTTriggerFired()){
    GetESDsData(kMultiplicityFired)->Fill(nESDTracks, nHLTesdTracks);
    GetESDsData(kNClsFired)->Fill(nClsOffline, nClsHLT);
    
    for(Int_t i = 0; i < nESDTracks; i++){ 
      AliESDtrack *  esdTrk = esd->GetTrack(i); 
      GetESDsData(kPOfflineFired)->Fill(esdTrk->P()); 
      GetESDsData(kPtOfflineFired)->Fill(esdTrk->GetSignedPt());
      GetESDsData(kNClsPerTrkOfflineFired)->Fill(esdTrk->GetTPCNcls());
      GetESDsData(kPhiOfflineFired)->Fill(esdTrk->Phi()*TMath::RadToDeg());
      GetESDsData(kEtaOfflineFired)->Fill(esdTrk->Eta());
    } 
   
    for(Int_t i = 0; i < nHLTesdTracks; i++){  
      AliESDtrack *  hltEsdTrk = hltesd->GetTrack(i);  
      GetESDsData(kPHLTFired)->Fill(hltEsdTrk->P());  
      GetESDsData(kPtHLTFired)->Fill(hltEsdTrk->GetSignedPt());
      GetESDsData(kNClsPerTrkHLTFired)->Fill(hltEsdTrk->GetTPCNcls());
      GetESDsData(kPhiHLTFired)->Fill(hltEsdTrk->Phi()*TMath::RadToDeg());
      GetESDsData(kEtaHLTFired)->Fill(hltEsdTrk->Eta());
    }  

    
  }
  
}
