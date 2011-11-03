//-*- Mode: C++ -*-
// $Id: AliHLTMultiplicityCorrelations.cxx  $
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <jochen@thaeder.de>                    *
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

/** @file   AliHLTMultiplicityCorrelations.cxx
    @author Jochen Thaeder <jochen@thaeder.de>
    @brief  Correlation plots for multiplicity studies
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"

#include "AliHLTMultiplicityCorrelations.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTMultiplicityCorrelations)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliHLTMultiplicityCorrelations::AliHLTMultiplicityCorrelations() :
  TNamed(),
  fHistList(NULL),
  fESDEvent(NULL),
  fESDZDC(NULL),
  fESDVZERO(NULL),
  fESDTrackCuts(NULL),
  fProcessTPC(kTRUE), fProcessSPD(kTRUE),
  fProcessVZERO(kTRUE), fProcessZDC(kTRUE),
  fEsdTracks(0), fEsdTracksA(0),
  fTpcTracks(0), fTpcTracksA(0),
  fTpcTracksRef(0),
  fVzeroMult(0.), fVzeroTriggerMult(0.),
  fSpdNClusters(0), fSpdNClustersInner(0), fSpdNClustersOuter(0),
  fVzeroBinning(1500), fVzeroBinningMin(0.), fVzeroBinningMax(30000.),
  fTpcBinning(800),fTpcBinningMin(0.),fTpcBinningMax(8000.),
  fZdcBinning(280),fZdcBinningMin(0.),fZdcBinningMax(140.),
  fZemBinning(100),fZemBinningMin(0.),fZemBinningMax(5.),
  fSpdBinning(750),fSpdBinningMin(0.),fSpdBinningMax(15000.) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
}
//##################################################################################
AliHLTMultiplicityCorrelations::AliHLTMultiplicityCorrelations(Char_t* name, Char_t* title) : 
  TNamed(name,title), 
  fHistList(NULL),
  fESDEvent(NULL),
  fESDZDC(NULL),
  fESDVZERO(NULL),
  fESDTrackCuts(NULL),
  fProcessTPC(kTRUE), fProcessSPD(kTRUE),
  fProcessVZERO(kTRUE), fProcessZDC(kTRUE),
  fEsdTracks(0), fEsdTracksA(0),
  fTpcTracks(0), fTpcTracksA(0),
  fTpcTracksRef(0),
  fVzeroMult(0.), fVzeroTriggerMult(0.),
  fSpdNClusters(0), fSpdNClustersInner(0), fSpdNClustersOuter(0),
  fVzeroBinning(1500), fVzeroBinningMin(0.), fVzeroBinningMax(30000.),
  fTpcBinning(800),fTpcBinningMin(0.),fTpcBinningMax(8000.),
  fZdcBinning(280),fZdcBinningMin(0.),fZdcBinningMax(140.),
  fZemBinning(100),fZemBinningMin(0.),fZemBinningMax(5.),
  fSpdBinning(750),fSpdBinningMin(0.),fSpdBinningMax(15000.) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
}

//##################################################################################
AliHLTMultiplicityCorrelations::~AliHLTMultiplicityCorrelations() {
  // see header file for class documentation

  if ( fHistList ) {
    fHistList->Clear();
    delete fHistList;
  }
  fHistList = NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                                   Initialize / Reset
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::Initialize() {
  // see header file for class documentation  

  Int_t iResult = 0;

  fHistList = new TList();
  fHistList->SetOwner(kTRUE);
  fHistList->SetName("MultiplicityCorrelations");
  iResult = SetupHistograms();
  
  if (fProcessZDC)   { HLTInfo ("Processing of ZDC enabled"); }
  if (fProcessTPC)   { HLTInfo ("Processing of TPC enabled"); }
  if (fProcessSPD)   { HLTInfo ("Processing of SPD enabled"); }
  if (fProcessVZERO) { HLTInfo ("Processing of VZERO enabled"); }

  return iResult;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::Initialize( const Char_t* listName) {
  // see header file for class documentation  

  Int_t iResult = 0;

  fHistList = new TList();
  fHistList->SetOwner(kTRUE);
  fHistList->SetName(Form("MultiplicityCorrelations_%s",listName));
  iResult = SetupHistograms();
  
  if (fProcessZDC)   { HLTInfo ("Processing of ZDC enabled"); }
  if (fProcessTPC)   { HLTInfo ("Processing of TPC enabled"); }
  if (fProcessSPD)   { HLTInfo ("Processing of SPD enabled"); }
  if (fProcessVZERO) { HLTInfo ("Processing of VZERO enabled"); }

  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Output - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::ProcessEvent( AliESDEvent *esd, AliESDVZERO* esdVZERO, 
						    Int_t nSpdClusters) {
  // see header file for class documentation  

  Int_t iResult = 0;

  if ( ! AddESDEvent(esd) ) {
    HLTWarning("No ESD event.");
    return -1;
  }
  
  if ( esdVZERO )
    fESDVZERO = esdVZERO;

  // -- TPC .. To be done before the others
  if (fESDEvent->GetNumberOfTracks() > 0 && fProcessTPC)
    iResult = ProcessTPC();
  
  // -- SPD
  fSpdNClusters = nSpdClusters;
  if (fProcessSPD)
    iResult = ProcessSPD();
    
  // -- VZERO
  if (fESDVZERO && fProcessVZERO)
    iResult = ProcessVZERO();

  // -- ZDC and Correlations
  if (fESDZDC && fProcessZDC)
    iResult = ProcessZDC();
 
  return iResult;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
///                                                                              ///
//////                             PRIVATE                                    //////
///                                                                              ///
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

/*
 * ---------------------------------------------------------------------------------
 *                          Setup / Initialize - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Bool_t AliHLTMultiplicityCorrelations::AddESDEvent( AliESDEvent* esd ) {
  // see header file for class documentation  

  fESDEvent = esd;
  
  // -- Check for ESD
  if ( !fESDEvent ) {
    HLTWarning("No ESD event present.");
    return kFALSE;
  }
  
  // -- Check for PrimaryVertex
  if ( !esd->GetPrimaryVertexTracks() ){
    HLTError("No Vertex present.");
    return kFALSE;
  }
  
  fESDZDC = esd->GetESDZDC();
  if ( !fESDZDC ) {
    HLTInfo("No ZDC information !");
  }
  
  fESDVZERO = esd->GetVZEROData();
  if ( !fESDVZERO ) {
    HLTInfo("No VZERO information !");
  }
  
  return kTRUE;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::SetupHistograms() {
  // see header file for class documentation  

  Int_t iResult = 0;
  
  if (fProcessVZERO) iResult = SetupVZERO();
  if (fProcessZDC)   iResult = SetupZDC();
  if (fProcessTPC)   iResult = SetupTPC();
  if (fProcessSPD)   iResult = SetupSPD();
  
  iResult = SetupCorrelations();

  return iResult;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::SetupVZERO() {
  // see header file for class documentation  

  // VzeroMult
  fHistList->Add(new TH1F("fVzeroMult",   "Multiplicity^{VZERO};Multiplicity^{VZERO};N_{Events}",   
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fVzeroMultAC", "Multiplicity^{VZERO} A vs C;Multiplicity^{VZERO}A ;Multiplicity^{VZERO} C", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax)); 

  // VzeroTriggerMult 
  fHistList->Add(new TH1F("fVzeroTriggerMult",   "Trigger Multiplicity^{VZERO};Trigger Multiplicity^{VZERO};N_{Events}",   
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax)); 
  fHistList->Add(new TH2F("fVzeroTriggerMultAC", "Trigger Multiplicity^{VZERO} A vs C;Trigger Multiplicity^{VZERO}A ;Trigger Multiplicity^{VZERO} C", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax)); 

  return 0;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::SetupZDC() {
  // see header file for class documentation  
  
  // E_{ZDC}
  fHistList->Add(new TH1F("fZdcEzdc",  "E_{ZDC};E_{ZDC} (TeV);N_{Events}",   
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax));
  fHistList->Add(new TH2F("fZdcEzdcAEzdcC", "E_{ZDC} A vs E_{ZDC} C;E_{ZDC} A (TeV); E_{ZDC} C (TeV)",   
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax, fZdcBinning,fZdcBinningMin,fZdcBinningMax));

  // E_{ZEM} vs E_{ZDC} 
  fHistList->Add(new TH2F("fZdcEzemEzdc", "E_{ZEM} vs E_{ZDC};E_{ZEM} (TeV); E_{ZDC} (TeV)",   
			  fZemBinning,fZemBinningMin,fZemBinningMax, fZdcBinning,fZdcBinningMin,fZdcBinningMax));

  return 0;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::SetupTPC() {
  // see header file for class documentation  

  // Multiplicty
  fHistList->Add(new TH1F("fTpcNch0", "TPC N_{ch} esdTracks; TPC N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH1F("fTpcNch1", "TPC N_{ch} accepted esdTracks; TPC N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH1F("fTpcNch2", "TPC N_{ch} tpcTracks; TPC N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH1F("fTpcNch3", "TPC N_{ch} accepted tpcTracks; TPC N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH1F("fTpcNch4", "TPC N_{ch} reference Tracks; TPC N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));

  return 0;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::SetupCorrelations() {
  // see header file for class documentation  

  // ----------------------------------------------------
  //   ZDC vs TPC
  // ----------------------------------------------------
  if (fProcessTPC && fProcessZDC) {

    // E_{ZDC} vs N_{ch}
    fHistList->Add(new TH2F("fCorrEzdcTpc",    "E_{ZDC} vs N_{ch}^{TPC};E_{ZDC} (TeV);N_{ch}^{TPC}",   
			    fZdcBinning,fZdcBinningMin,fZdcBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));

    fHistList->Add(new TH2F("fCorrEzdcRefTpc", "E_{ZDC} vs N_{ch}^{Ref TPC};E_{ZDC} (TeV);N_{ch}^{Ref TPC}",   
			    fZdcBinning,fZdcBinningMin,fZdcBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  }
  // ----------------------------------------------------
  //   ZDC vs VZERO
  // ----------------------------------------------------
  if (fProcessZDC && fProcessVZERO) {

    // E_{ZDC} vs Multiplicty VZERO
    fHistList->Add(new TH2F("fCorrEzdcVzero", 
			    "E_{ZDC} vs Multiplicity^{VZERO};E_{ZDC} (TeV);Multiplicity^{VZERO}",  
			    fZdcBinning,fZdcBinningMin,fZdcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
    fHistList->Add(new TH2F("fCorrEzdcTriggerVzero", 
			    "E_{ZDC} vs Trigger Multiplicity^{VZERO};E_{ZDC} (TeV);Trigger Multiplicity^{VZERO}",  
			    fZdcBinning,fZdcBinningMin,fZdcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  }
  // ----------------------------------------------------
  //   VZERO vs TPC
  // ----------------------------------------------------
  if ( fProcessTPC && fProcessVZERO ) {

    fHistList->Add(new TH2F("fCorrVzeroTpc", 
			    "Multiplicity^{VZERO} vs N_{ch}^{TPC};Multiplicity^{VZERO};N_{ch}^{TPC}", 
			    fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
    fHistList->Add(new TH2F("fCorrVzeroRefTpc", 
			    "Multiplicity^{VZERO} vs N_{ch}^{Ref TPC};Multiplicity^{VZERO};N_{ch}^{Ref TPC}", 
			    fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  
    fHistList->Add(new TH2F("fCorrTriggerVzeroTpc", 
			    "Trigger Multiplicity^{VZERO} vs N_{ch}^{TPC};Multiplicity^{VZERO};N_{ch}^{TPC}", 
			    fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
    fHistList->Add(new TH2F("fCorrTriggerVzeroRefTpc", 
			    "Trigger Multiplicity^{VZERO} vs N_{ch}^{Ref TPC};Multiplicity^{VZERO};N_{ch}^{Ref TPC}", 
			    fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  }
  // ----------------------------------------------------
  //   SPD vs TPC
  // ----------------------------------------------------
  if ( fProcessSPD && fProcessTPC ) {
    fHistList->Add(new TH2F("fCorrSpdOuterTpc", "N_{clusters}^{SPD}_{Outer} vs N_{ch}^{TPC};N_{clusters}^{SPD};N_{ch}^{TPC}",   
			    fSpdBinning,fSpdBinningMin,fSpdBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));  
    fHistList->Add(new TH2F("fCorrSpdOuterRefTpc", "N_{clusters}^{SPD}_{Outer} vs N_{ch}^{Ref TPC};N_{clusters}^{SPD};N_{ch}^{Ref TPC}",   
			    fSpdBinning,fSpdBinningMin,fSpdBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));  
    

  }
  // ----------------------------------------------------
  //   SPD vs VZERO
  // ----------------------------------------------------
  if ( fProcessSPD && fProcessVZERO ) {
    fHistList->Add(new TH2F("fCorrVzeroSpdOuter", 
			    "Multiplicity^{VZERO} vs N_{ch}^{SPD}_{Outer};Multiplicity^{VZERO};N^{SPD}", 
			    fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fSpdBinning,fSpdBinningMin,fSpdBinningMax));
    
    fHistList->Add(new TH2F("fCorrVzeroTriggerSpdOuter", 
			    "Trigger Multiplicity^{VZERO} vs N_{ch}^{SPD}_{Outer};Trigger Multiplicity^{VZERO};N^{SPD}", 
			    fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fSpdBinning,fSpdBinningMin,fSpdBinningMax));
    
  }
  // ----------------------------------------------------
  //   SPD vs ZDC
  // ----------------------------------------------------
  if ( fProcessSPD && fProcessZDC) {
    fHistList->Add(new TH2F("fCorrEzdcSpdOuter", "E_{ZDC} vs N_{ch}^{SPD}_{Outer};E_{ZDC} (TeV);N^{SPD}",   
			    fZdcBinning,fZdcBinningMin,fZdcBinningMax, fSpdBinning,fSpdBinningMin,fSpdBinningMax));
  }

  return 0;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::SetupSPD() {
  // see header file for class documentation  
  
  fHistList->Add(new TH1F("fSpdNClusters", "Multplicity_{SPD};Multplicity_{SPD};N_{Events}",   
			  fSpdBinning,fSpdBinningMin,fSpdBinningMax));

  fHistList->Add(new TH1F("fSpdNClustersInner", "Multplicity_{SPD} Layer 0;Multplicity_{SPD};N_{Events}",   
			  fSpdBinning,fSpdBinningMin,fSpdBinningMax));

  fHistList->Add(new TH1F("fSpdNClustersOuter", "Multplicity_{SPD} Layer 1;Multplicity_{SPD};N_{Events}",   
			  fSpdBinning,fSpdBinningMin,fSpdBinningMax));

  return 0;
}


/*
 * ---------------------------------------------------------------------------------
 *                               Process - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::ProcessTPC() {
  // see header file for class documentation  

  Int_t iResult = 0 ;

  fEsdTracks = 0;  
  fEsdTracksA = 0;  
  fTpcTracks = 0;  
  fTpcTracksA = 0;  
  fTpcTracksRef = 0;  
  
  // -----------------------------------------------------------------
  for (Int_t idx = 0; idx < fESDEvent->GetNumberOfTracks(); idx++) {
    AliESDtrack* esdTrack = fESDEvent->GetTrack(idx);
    if (!esdTrack) continue;
    ++fEsdTracks;
    
    if (!fESDTrackCuts->AcceptTrack(esdTrack)) continue;
    ++fEsdTracksA;
  }
  // -----------------------------------------------------------------
  for (Int_t idx = 0; idx < fESDEvent->GetNumberOfTracks(); idx++) {
    AliESDtrack* esdTrack = fESDEvent->GetTrack(idx);
    if (!esdTrack) continue;

    // -- TPC only
    const AliExternalTrackParam *tpcTrack = esdTrack->GetTPCInnerParam();
    if (!tpcTrack) continue;
    ++fTpcTracks;


    if (!fESDTrackCuts->AcceptTrack(esdTrack)) continue;
    ++fTpcTracksA;
  }
  // -----------------------------------------------------------------

  fTpcTracksRef = fESDTrackCuts->GetReferenceMultiplicity(fESDEvent,kTRUE);

  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch0")))->Fill(fEsdTracks);
  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch1")))->Fill(fEsdTracksA);
  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch2")))->Fill(fTpcTracks);
  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch3")))->Fill(fTpcTracksA);
  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch4")))->Fill(fTpcTracksRef);

  return iResult;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::ProcessVZERO() {
  // see header file for class documentation  

  Int_t iResult = 0 ;

  Float_t vzeroMultA = fESDVZERO->GetMTotV0A();
  Float_t vzeroMultC = fESDVZERO->GetMTotV0C();
  fVzeroMult = vzeroMultA + vzeroMultC;

  (static_cast<TH1F*>(fHistList->FindObject("fVzeroMult")))->Fill(fVzeroMult);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroMultAC")))->Fill(vzeroMultA,vzeroMultC);

  Float_t vzeroTriggerMultA = Float_t(fESDVZERO->GetTriggerChargeA());
  Float_t vzeroTriggerMultC = Float_t(fESDVZERO->GetTriggerChargeC());
  fVzeroTriggerMult  = vzeroTriggerMultA + vzeroTriggerMultC;

  (static_cast<TH1F*>(fHistList->FindObject("fVzeroTriggerMult")))->Fill(fVzeroTriggerMult);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroTriggerMultAC")))->Fill(vzeroTriggerMultA,vzeroTriggerMultC);

  // -- VZERO - TPC correlations
  if (fESDEvent->GetNumberOfTracks() > 0 && fProcessTPC) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroTpc")))->Fill(fVzeroMult, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroRefTpc")))->Fill(fVzeroMult, fTpcTracksRef);

    (static_cast<TH2F*>(fHistList->FindObject("fCorrTriggerVzeroTpc")))->Fill(fVzeroTriggerMult, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrTriggerVzeroRefTpc")))->Fill(fVzeroTriggerMult, fTpcTracksRef);
  }

  // -- VZERO - SPD correlations
  if (fESDEvent->GetNumberOfTracks() > 0 && fProcessSPD) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroSpdOuter")))->Fill(fVzeroMult, fSpdNClustersOuter);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroTriggerSpdOuter")))->Fill(fVzeroTriggerMult, fSpdNClustersOuter);
  }
  
  return iResult;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::ProcessZDC() {
  // see header file for class documentation  

  Int_t iResult = 0 ;

  Double_t zdcEznA = fESDZDC->GetZDCN2Energy() / 1000.;
  Double_t zdcEznC = fESDZDC->GetZDCN1Energy() / 1000.;
  
  Double_t zdcEzpA = fESDZDC->GetZDCP2Energy() / 1000.;
  Double_t zdcEzpC = fESDZDC->GetZDCP1Energy() / 1000.;
  
  Double_t zdcEA = (zdcEznA + zdcEzpA);
  Double_t zdcEC = (zdcEznC + zdcEzpC);
  Double_t zdcE = zdcEA + zdcEC;
  
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzdc")))->Fill(zdcE);
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzdcAEzdcC")))->Fill(zdcEA, zdcEC);

  Double_t zdcEzemA = fESDZDC->GetZDCEMEnergy(1) / 1000.;
  Double_t zdcEzemC = fESDZDC->GetZDCEMEnergy(0) / 1000.;
  Double_t zdcEzem  = zdcEzemA + zdcEzemC;
  
  (static_cast<TH2F*>(fHistList->FindObject("fZdcEzemEzdc")))->Fill(zdcEzem, zdcE);
 
  // -- ZDC - TPC correlations
  if (fESDEvent->GetNumberOfTracks() > 0 && fProcessTPC) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcTpc")))->Fill(zdcE, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcRefTpc")))->Fill(zdcE, fTpcTracksRef);
  }
  
  // -- ZDC - SPD correlations
  if (fProcessSPD) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcSpdOuter")))->Fill(zdcE, fSpdNClustersOuter);
  }

  // -- VZERO - ZDC correlations
  if (fESDVZERO && fProcessVZERO) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcVzero")))->Fill(zdcE, fVzeroMult);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcTriggerVzero")))->Fill(zdcE, fVzeroTriggerMult);
  }

  return iResult;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::ProcessSPD() {
  // see header file for class documentation
  
  (static_cast<TH2F*>(fHistList->FindObject("fSpdNClusters")))->Fill(fSpdNClusters);
  (static_cast<TH2F*>(fHistList->FindObject("fSpdNClustersInner")))->Fill(fSpdNClustersInner);
  (static_cast<TH2F*>(fHistList->FindObject("fSpdNClustersOuter")))->Fill(fSpdNClustersOuter);

  
  // -- SPD vs TPC correlations
  if (fProcessTPC) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrSpdOuterTpc")))->Fill(fSpdNClustersOuter,fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrSpdOuterRefTpc")))->Fill(fSpdNClustersOuter,fTpcTracksRef);
  }

  return 0;
}
