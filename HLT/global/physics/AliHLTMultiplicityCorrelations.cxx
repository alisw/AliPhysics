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
  fHistList(NULL),
  fESDEvent(NULL),
  fESDZDC(NULL),
  fESDVZERO(NULL),
  fESDTrackCuts(NULL),
  fEsdTracks(0), fEsdTracksA(0),
  fTpcTracks(0), fTpcTracksA(0),
  fVzeroMult(0.), fVzeroMultA(0.), fVzeroMultC(0.),
  fVzeroMultFlagged(0.), fVzeroMultFlaggedA(0.), fVzeroMultFlaggedC(0.),
  fVzeroBinning(350), fVzeroBinningMin(0.), fVzeroBinningMax(35000.),
  fTpcBinning(200),fTpcBinningMin(0.),fTpcBinningMax(8000.),
  fZdcBinning(280),fZdcBinningMin(0.),fZdcBinningMax(140.),
  fZemBinning(100),fZemBinningMin(0.),fZemBinningMax(5.),
  fZnpBinning(200),fZnpBinningMin(0.),fZnpBinningMax(100.) {
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

  iResult = SetupHistograms();

  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Output - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::ProcessEvent( AliESDEvent *esd ) {
  // see header file for class documentation  

  Int_t iResult = 0;

  if ( ! AddESDEvent(esd) ) {
    HLTWarning("No ESD event.");
    return -1;
  }

  // -- TPC .. To be done before the others
  if (fESDEvent->GetNumberOfTracks() > 0)
    iResult = ProcessTPC();

  // -- VZERO
  if (fESDVZERO)
    iResult = ProcessVZERO();

  // -- ZDC and Correlations
  if (fESDZDC)
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
  
  iResult = SetupVZERO();
  iResult = SetupZDC();
  iResult = SetupTPC();
  iResult = SetupCorrelations();

  return iResult;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::SetupVZERO() {
  // see header file for class documentation  

  // VzeroMult
  fHistList->Add(new TH1F("fVzeroMult",  "Multiplicity^{VZERO};Multiplicity^{VZERO};N_{Events}",   
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH1F("fVzeroMultA", "Multiplicity^{VZERO} A;Multiplicity^{VZERO};N_{Events}", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH1F("fVzeroMultC", "Multiplicity^{VZERO} C;Multiplicity^{VZERO};N_{Events}", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax)); 

  // Flagged VzeroMult
  fHistList->Add(new TH1F("fVzeroFlaggedMult",  "Multiplicity_{flagged}^{VZERO};Multiplicity^{VZERO};N_{Events}",   
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax)); 
  fHistList->Add(new TH1F("fVzeroFlaggedMultA", "Multiplicity_{flagged}^{VZERO} A;Multiplicity^{VZERO};N_{Events}", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax)); 
  fHistList->Add(new TH1F("fVzeroFlaggedMultC", "Multiplicity_{flagged}^{VZERO} C;Multiplicity^{VZERO};N_{Events}", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax)); 
  
  fHistList->Add(new TH1F("fVzeroTime",  "Time;Time;N_{Events}",   500, 0, 1000));
  fHistList->Add(new TH1F("fVzeroTimeA", "Time A;Time;N_{Events}", 500, 0, 1000));
  fHistList->Add(new TH1F("fVzeroTimeC", "Time B;Time;N_{Events}", 500, 0, 1000));

  fHistList->Add(new TH1F("fVzeroADC",  "ADC;ADC;N_{Events}",   500, 0, 150000));
  fHistList->Add(new TH1F("fVzeroADCA", "ADC A;ADC;N_{Events}", 500, 0, 100000));
  fHistList->Add(new TH1F("fVzeroADCC", "ADC B;ADC;N_{Events}", 500, 0, 100000));

  return 0;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::SetupZDC() {
  // see header file for class documentation  
  
  // E_{ZN}
  fHistList->Add(new TH1F("fZdcEzn",  "E_{ZN};E_{ZN} (TeV);N_{Events}",   
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax));
  fHistList->Add(new TH1F("fZdcEznA", "E_{ZN} A;E_{ZN} (TeV);N_{Events}", 
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax));
  fHistList->Add(new TH1F("fZdcEznC", "E_{ZN} C;E_{ZN} (TeV);N_{Events}", 
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax));
  
  // E_{ZP}
  fHistList->Add(new TH1F("fZdcEzp",  "E_{ZP};E_{ZP} (TeV);N_{Events}",   
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax));
  fHistList->Add(new TH1F("fZdcEzpA", "E_{ZP} A;E_{ZP} (TeV);N_{Events}", 
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax));
  fHistList->Add(new TH1F("fZdcEzpC", "E_{ZP} C;E_{ZP} (TeV);N_{Events}", 
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax));
  
  // E_{ZDC}
  fHistList->Add(new TH1F("fZdcEzdc",  "E_{ZDC};E_{ZDC} (TeV);N_{Events}",   
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax));
  fHistList->Add(new TH1F("fZdcEzdcA", "E_{ZDC} A;E_{ZDC} (TeV);N_{Events}", 
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax));
  fHistList->Add(new TH1F("fZdcEzdcC", "E_{ZDC} C;E_{ZDC} (TeV);N_{Events}", 
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax));

  // E_{ZEM}
  fHistList->Add(new TH1F("fZdcEzem",  "E_{ZEM};E_{ZEM} (TeV);N_{Events}",   
			  fZemBinning,fZemBinningMin,fZemBinningMax));
  fHistList->Add(new TH1F("fZdcEzemA", "E_{ZEM} A;E_{ZEM} (TeV);N_{Events}", 
			  fZemBinning,fZemBinningMin,fZemBinningMax));
  fHistList->Add(new TH1F("fZdcEzemC", "E_{ZEM} C;E_{ZEM} (TeV);N_{Events}", 
			  fZemBinning,fZemBinningMin,fZemBinningMax));
  
  // N_{Part}
  fHistList->Add(new TH1F("fZdcNpart",  "N_{Part} ;N_{Part};N_{Events}",  50,0,499));
  fHistList->Add(new TH1F("fZdcNpartA", "N_{Part} A;N_{Part};N_{Events}", 50,0,499));
  fHistList->Add(new TH1F("fZdcNpartC", "N_{Part} C;N_{Part};N_{Events}", 50,0,499));
  
  // b
  fHistList->Add(new TH1F("fZdcB",  "b;b {fm);N_{Events}",   31,0,30));
  fHistList->Add(new TH1F("fZdcBA", "b A;b (fm);N_{Events}", 31,0,30));
  fHistList->Add(new TH1F("fZdcBC", "b C;b (fm);N_{Events}", 31,0,30));

  // E_{ZEM} vs E_{ZDC} 
  fHistList->Add(new TH2F("fZdcEzemEzdc", "E_{ZEM} vs E_{ZDC};E_{ZEM} (TeV); E_{ZDC} (TeV)",   
			  fZemBinning,fZemBinningMin,fZemBinningMax, fZdcBinning,fZdcBinningMin,fZdcBinningMax));
  fHistList->Add(new TH2F("fZdcEzemEzdcA","E_{ZEM} vs E_{ZDC} A;E_{ZEM} (TeV); E_{ZDC} (TeV)", 
			  fZemBinning,fZemBinningMin,fZemBinningMax, fZdcBinning,fZdcBinningMin,fZdcBinningMax));
  fHistList->Add(new TH2F("fZdcEzemEzdcC","E_{ZEM} vs E_{ZDC} C;E_{ZEM} (TeV); E_{ZDC} (TeV)", 
			  fZemBinning,fZemBinningMin,fZemBinningMax, fZdcBinning,fZdcBinningMin,fZdcBinningMax));

  return 0;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::SetupTPC() {
  // see header file for class documentation  

  Int_t    n = 50;
  Double_t s = 0.1;
  Double_t e = 50.;
  
  // -- Create LogPtBinning
  // ------------------------
  Double_t logMin = TMath::Log10(s);
  Double_t logMax = TMath::Log10(e);
  Double_t binwidth = (logMax-logMin)/n;
  
  Double_t *logBinning = new Double_t[n+1];

  logBinning[0] = s;
  for (Int_t ii = 1; ii <= n; ii++)
    logBinning[ii] = s + TMath::Power(10, logMin + ii*binwidth);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  // dN_{ch} / dP_{T}
  fHistList->Add(new TH1F("fTpcPt0", "dN_{ch} / dP_{T} esdTracks; P_{T} (GeV/c);dN/dP_{T} (c/GeV)", n, logBinning));  
  fHistList->Add(new TH1F("fTpcPt1", "dN_{ch} / dP_{T} accepted esdTracks; P_{T} (GeV/c);dN/dP_{T} (c/GeV)", n, logBinning));  
  fHistList->Add(new TH1F("fTpcPt2", "dN_{ch} / dP_{T} tpcTracks; P_{T} (GeV/c);dN/dP_{T} (c/GeV)", n, logBinning));  
  fHistList->Add(new TH1F("fTpcPt3", "dN_{ch} / dP_{T} accepted tpcTracks; P_{T} (GeV/c);dN/dP_{T} (c/GeV)", n, logBinning));  

  // Multiplicty
  fHistList->Add(new TH1F("fTpcNch0", "N_{ch} esdTracks; N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH1F("fTpcNch1", "N_{ch} accepted esdTracks; N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH1F("fTpcNch2", "N_{ch} tpcTracks; N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH1F("fTpcNch3", "N_{ch} accepted tpcTracks; N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));

  return 0;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::SetupCorrelations() {
  // see header file for class documentation  

  // ----------------------------------------------------
  //   ZDC vs TPC
  // ----------------------------------------------------

  // E_{ZDC} vs N_{ch}
  fHistList->Add(new TH2F("fCorrEzdcNch", "E_{ZDC} vs N_{ch}^{TPC};E_{ZDC} (TeV);N_{ch}^{TPC}",   
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrEzdcANch","E_{ZDC} vs N_{ch}^{TPC} A;E_{ZDC} (TeV);N_{ch}^{TPC}", 
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrEzdcCNch","E_{ZDC} vs N_{ch}^{TPC} C;E_{ZDC} (TeV);N_{ch}^{TPC}", 
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  
  // E_{ZEM} vs N_{ch}
  fHistList->Add(new TH2F("fCorrEzemNch", "E_{ZEM} vs N_{ch}^{TPC};E_{ZEM} (TeV);N_{ch}^{TPC}",   
			  fZemBinning,fZemBinningMin,fZemBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrEzemANch","E_{ZEM} vs N_{ch}^{TPC} A;E_{ZEM} (TeV);N_{ch}^{TPC}", 
			  fZemBinning,fZemBinningMin,fZemBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrEzemCNch","E_{ZEM} vs N_{ch}^{TPC} C;E_{ZEM} (TeV);N_{ch}^{TPC}", 
			  fZemBinning,fZemBinningMin,fZemBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  
  // E_{ZP} vs N_{ch}
  fHistList->Add(new TH2F("fCorrEzpNch", "E_{ZP} vs N_{ch}^{TPC};E_{ZP} (TeV);N_{ch}^{TPC}",   
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrEzpANch","E_{ZP} vs N_{ch}^{TPC} A;E_{ZP} (TeV);N_{ch}^{TPC}", 
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrEzpCNch","E_{ZP} vs N_{ch}^{TPC} C;E_{ZP} (TeV);N_{ch}^{TPC}", 
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  
  // E_{ZN} vs N_{ch}
  fHistList->Add(new TH2F("fCorrEznNch", "E_{ZN} vs N_{ch}^{TPC};E_{ZN} (TeV);N_{ch}^{TPC}",   
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrEznANch","E_{ZN} vs N_{ch}^{TPC} A;E_{ZN} (TeV);N_{ch}^{TPC}", 
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrEznCNch","E_{ZN} vs N_{ch}^{TPC} C;E_{ZN} (TeV);N_{ch}^{TPC}", 
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  
  // N_{Part} vs N_{ch}
  fHistList->Add(new TH2F("fCorrZdcNpartNch", "N_{part} vs N_{ch}^{TPC};N_{Events};N_{ch}^{TPC}",   
			  50,0,499, fTpcBinning,fTpcBinningMin,fTpcBinningMax)); 
  fHistList->Add(new TH2F("fCorrZdcNpartANch","N_{part} vs N_{ch}^{TPC} A;N_{Events};N_{ch}^{TPC}", 
			  50,0,499, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrZdcNpartCNch","N_{part} vs N_{ch}^{TPC} C;N_{Events};N_{ch}^{TPC}", 
			  50,0,499, fTpcBinning,fTpcBinningMin,fTpcBinningMax));

  // b vs N_{ch}
  fHistList->Add(new TH2F("fCorrZdcbNch", "b_{ZDC} vs N_{ch}^{TPC};b (fm);N_{ch}^{TPC}",  
			  31,0,30, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrZdcbANch","b_{ZDC} vs N_{ch}^{TPC} A;b (fm);N_{ch}^{TPC}", 
			  31,0,30, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrZdcbCNch","b_{ZDC} vs N_{ch}^{TPC} C;b (fm);N_{ch}^{TPC}", 
			  31,0,30, fTpcBinning,fTpcBinningMin,fTpcBinningMax));

  // ----------------------------------------------------
  //   ZDC vs VZERO
  // ----------------------------------------------------

  // E_{ZDC} vs Multiplicty VZERO
  fHistList->Add(new TH2F("fCorrEzdcVzero", 
			  "E_{ZDC} vs Multiplicity^{VZERO};E_{ZDC} (TeV);Multiplicity^{VZERO}",  
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzdcVzeroA",
			  "E_{ZDC} vs Multiplicity^{VZERO} A;E_{ZDC} (TeV);Multiplicity^{VZERO}",
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzdcVzeroC",
			  "E_{ZDC} vs Multiplicity^{VZERO} C;E_{ZDC} (TeV);Multiplicity^{VZERO}",
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  
  // E_{ZEM} vs Multiplicty VZERO
  fHistList->Add(new TH2F("fCorrEzemVzero", 
			  "E_{ZEM} vs Multiplicity^{VZERO};E_{ZEM} (TeV);Multiplicity^{VZERO}",  
			  fZemBinning,fZemBinningMin,fZemBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzemVzeroA",
			  "E_{ZEM} vs Multiplicity^{VZERO} A;E_{ZEM} (TeV);Multiplicity^{VZERO}",
			  fZemBinning,fZemBinningMin,fZemBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzemVzeroC",
			  "E_{ZEM} vs Multiplicity^{VZERO} C;E_{ZEM} (TeV);Multiplicity^{VZERO}",
			  fZemBinning,fZemBinningMin,fZemBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  
  // E_{ZP} vs Multiplicty VZERO
  fHistList->Add(new TH2F("fCorrEzpVzero", 
			  "E_{ZP} vs Multiplicity^{VZERO};E_{ZP} (TeV);Multiplicity^{VZERO}",  
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzpVzeroA",
			  "E_{ZP} vs Multiplicity^{VZERO} A;E_{ZP} (TeV);Multiplicity^{VZERO}",
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzpVzeroC",
			  "E_{ZP} vs Multiplicity^{VZERO} C;E_{ZP} (TeV);Multiplicity^{VZERO}",
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  
  // E_{ZN} vs Multiplicty VZERO
  fHistList->Add(new TH2F("fCorrEznVzero", 
			  "E_{ZN} vs Multiplicity^{VZERO};E_{ZN} (TeV);Multiplicity^{VZERO}",   
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEznVzeroA",
			  "E_{ZN} vs Multiplicity^{VZERO} A;E_{ZN} (TeV);Multiplicity^{VZERO}",
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEznVzeroC",
			  "E_{ZN} vs Multiplicity^{VZERO} C;E_{ZN} (TeV);Multiplicity^{VZERO}",
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  
  // N_{Part} vs Multiplicty VZERO
  fHistList->Add(new TH2F("fCorrZdcNpartVzero", 
			  "N_{part} vs Multiplicity^{VZERO};N_{Events};Multiplicity^{VZERO}",   
			  50,0,499, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax)); 
  fHistList->Add(new TH2F("fCorrZdcNpartVzeroA",
			  "N_{part} vs Multiplicity^{VZERO} A;N_{Events};Multiplicity^{VZERO}", 
			  50,0,499, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrZdcNpartVzeroC",
			  "N_{part} vs Multiplicity^{VZERO} C;N_{Events};Multiplicity^{VZERO}",
			  50,0,499, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));

  // b vs Multiplicty VZERO
  fHistList->Add(new TH2F("fCorrZdcbVzero", "b_{ZDC} vs Multiplicity^{VZERO};b (fm);Multiplicity^{VZERO}",    
			  31,0,30, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrZdcbVzeroA","b_{ZDC} vs Multiplicity^{VZERO} A;b (fm);Multiplicity^{VZERO}",  
			  31,0,30, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrZdcbVzeroC","b_{ZDC} vs Multiplicity^{VZERO} C;b (fm);Multiplicity^{VZERO}",  
			  31,0,30, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));

  // -- -- -- -- -- -- -- 

  // E_{ZDC} vs Multiplicty VZERO flagged
  fHistList->Add(new TH2F("fCorrEzdcVzeroF", 
			  "E_{ZDC} vs Multiplicity_{flagged}^{VZERO};E_{ZDC} (TeV);Multiplicity^{VZERO}",   
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzdcVzeroFA",
			  "E_{ZDC} vs Multiplicity_{flagged}^{VZERO} A;E_{ZDC} (TeV);Multiplicity^{VZERO}", 
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzdcVzeroFC",
			  "E_{ZDC} vs Multiplicity_{flagged}^{VZERO} C;E_{ZDC} (TeV);Multiplicity^{VZERO}", 
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  
  // E_{ZEM} vs Multiplicty VZERO flagged
  fHistList->Add(new TH2F("fCorrEzemVzeroF", 
			  "E_{ZEM} vs Multiplicity_{flagged}^{VZERO};E_{ZEM} (TeV);Multiplicity^{VZERO}",   
			  fZemBinning,fZemBinningMin,fZemBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzemVzeroFA",
			  "E_{ZEM} vs Multiplicity_{flagged}^{VZERO} A;E_{ZEM} (TeV);Multiplicity^{VZERO}", 
			  fZemBinning,fZemBinningMin,fZemBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzemVzeroFC",
			  "E_{ZEM} vs Multiplicity_{flagged}^{VZERO} C;E_{ZEM} (TeV);Multiplicity^{VZERO}", 
			  fZemBinning,fZemBinningMin,fZemBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  
  // E_{ZP} vs Multiplicty VZERO flagged
  fHistList->Add(new TH2F("fCorrEzpVzeroF", 
			  "E_{ZP} vs Multiplicity_{flagged}^{VZERO};E_{ZP} (TeV);Multiplicity^{VZERO}",   
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzpVzeroFA",
			  "E_{ZP} vs Multiplicity_{flagged}^{VZERO} A;E_{ZP} (TeV);Multiplicity^{VZERO}", 
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEzpVzeroFC",
			  "E_{ZP} vs Multiplicity_{flagged}^{VZERO} C;E_{ZP} (TeV);Multiplicity^{VZERO}",
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  
  // E_{ZN} vs Multiplicty VZERO flagged
  fHistList->Add(new TH2F("fCorrEznVzeroF", 
			  "E_{ZN} vs Multiplicity_{flagged}^{VZERO};E_{ZN} (TeV);Multiplicity^{VZERO}",  
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEznVzeroFA",
			  "E_{ZN} vs Multiplicity_{flagged}^{VZERO} A;E_{ZN} (TeV);Multiplicity^{VZERO}", 
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrEznVzeroFC",
			  "E_{ZN} vs Multiplicity_{flagged}^{VZERO} C;E_{ZN} (TeV);Multiplicity^{VZERO}", 
			  fZnpBinning,fZnpBinningMin,fZnpBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  
  // N_{Part} vs Multiplicty VZERO flagged
  fHistList->Add(new TH2F("fCorrZdcNpartVzeroF", 
			  "N_{part} vs Multiplicity_{flagged}^{VZERO};N_{Events};Multiplicity^{VZERO}",   
			  50,0,499, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax)); 
  fHistList->Add(new TH2F("fCorrZdcNpartVzeroFA",
			  "N_{part} vs Multiplicity_{flagged}^{VZERO} A;N_{Events};Multiplicity^{VZERO}", 
			  50,0,499, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrZdcNpartVzeroFC",
			  "N_{part} vs Multiplicity_{flagged}^{VZERO} C;N_{Events};Multiplicity^{VZERO}", 
			  50,0,499, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));

  // b vs Multiplicty VZERO flagged
  fHistList->Add(new TH2F("fCorrZdcbVzeroF", 
			  "b_{ZDC} vs Multiplicity_{flagged}^{VZERO};b (fm);Multiplicity^{VZERO}",    
			  31,0,30, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrZdcbVzeroFA",
			  "b_{ZDC} vs Multiplicity_{flagged}^{VZERO} A;b (fm);Multiplicity^{VZERO}",  
			  31,0,30, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
  fHistList->Add(new TH2F("fCorrZdcbVzeroFC",
			  "b_{ZDC} vs Multiplicity_{flagged}^{VZERO} C;b (fm);Multiplicity^{VZERO}",  
			  31,0,30, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));

  // ----------------------------------------------------
  //   VZERO vs TPC
  // ----------------------------------------------------

  fHistList->Add(new TH2F("fCorrVzeroNch", 
			  "Multiplicity^{VZERO} vs N_{ch}^{TPC};Multiplicity^{VZERO};N_{ch}^{TPC}", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrVzeroANch", 
			  "Multiplicity^{VZERO} vs N_{ch}^{TPC} A;Multiplicity^{VZERO};N_{ch}^{TPC}", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrVzeroCNch", 
			  "Multiplicity^{VZERO} vs N_{ch}^{TPC} C;Multiplicity^{VZERO};N_{ch}^{TPC}", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));

  fHistList->Add(new TH2F("fCorrVzeroFNch", 
			  "Multiplicity^{VZERO} vs N_{ch}^{TPC};Multiplicity_{flagged}^{VZERO};N_{ch}^{TPC}", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrVzeroFANch", 
			  "Multiplicity^{VZERO} vs N_{ch}^{TPC} A;Multiplicity_{flagged}^{VZERO};N_{ch}^{TPC}", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH2F("fCorrVzeroFCNch", 
			  "Multiplicity^{VZERO} vs N_{ch}^{TPC} C;Multiplicity_{flagged}^{VZERO};N_{ch}^{TPC}", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));

  // ----------------------------------------------------
  //  
  // ----------------------------------------------------

  return 0;
}

/*
 * ---------------------------------------------------------------------------------
 *                               Process - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::ProcessTPC() {

  Int_t iResult = 0 ;

  TH1F* hPt0 = static_cast<TH1F*>(fHistList->FindObject("fTpcPt0")); // all
  TH1F* hPt1 = static_cast<TH1F*>(fHistList->FindObject("fTpcPt1")); // all accepted

  TH1F* hPt2 = static_cast<TH1F*>(fHistList->FindObject("fTpcPt2")); // GetTPCInnerParam
  TH1F* hPt3 = static_cast<TH1F*>(fHistList->FindObject("fTpcPt3")); // GetTPCInnerParam A

  fEsdTracks = 0;  
  fEsdTracksA = 0;  
  fTpcTracks = 0;  
  fTpcTracksA = 0;  
  
  // -----------------------------------------------------------------
  for (Int_t idx = 0; idx < fESDEvent->GetNumberOfTracks(); idx++) {

    // -- check  track
    AliESDtrack* esdTrack = fESDEvent->GetTrack(idx);
    if (!esdTrack) continue;

    ++fEsdTracks;
    hPt0->Fill(esdTrack->Pt());
    
    //______________________________________________
    if (!fESDTrackCuts->AcceptTrack(esdTrack)) continue;
    ++fEsdTracksA;

    hPt1->Fill(esdTrack->Pt());
  }
  // -----------------------------------------------------------------
  for (Int_t idx = 0; idx < fESDEvent->GetNumberOfTracks(); idx++) {

    // -- check  track
    AliESDtrack* esdTrack = fESDEvent->GetTrack(idx);
    if (!esdTrack) continue;

    // -- TPC only
    const AliExternalTrackParam *tpcTrack = esdTrack->GetTPCInnerParam();
    if (!tpcTrack) continue;
    ++fTpcTracks;
    hPt2->Fill(tpcTrack->Pt());
    
    //______________________________________________
    if (!fESDTrackCuts->AcceptTrack(esdTrack)) continue;
    ++fTpcTracksA;

    hPt3->Fill(tpcTrack->Pt());
  }
  // -----------------------------------------------------------------

  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch0")))->Fill(fEsdTracks);
  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch1")))->Fill(fEsdTracksA);
  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch2")))->Fill(fTpcTracks);
  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch3")))->Fill(fTpcTracksA);

  return iResult;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::ProcessVZERO() {

  Int_t iResult = 0 ;

  fVzeroMultA = fESDVZERO->GetMTotV0A();
  fVzeroMultC = fESDVZERO->GetMTotV0C();
  fVzeroMult  = fVzeroMultA + fVzeroMultC;

  (static_cast<TH1F*>(fHistList->FindObject("fVzeroMult")))->Fill(fVzeroMult);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroMultA")))->Fill(fVzeroMultA);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroMultC")))->Fill(fVzeroMultC);
  
  fVzeroMultFlaggedA = 0.;
  fVzeroMultFlaggedC = 0.;
  
  Double_t vzeroTimeA = 0.;
  Double_t vzeroTimeC = 0.;
  
  Double_t vzeroAdcA = 0.;
  Double_t vzeroAdcC = 0.;
  
  for (Int_t idx = 0; idx < 32; idx++) {
    
    if ( fESDVZERO->GetBBFlag(idx) )
      fVzeroMultFlaggedC += fESDVZERO->GetMultiplicityV0C(idx);
    if ( fESDVZERO->GetBBFlag(idx+32) )
      fVzeroMultFlaggedA += fESDVZERO->GetMultiplicityV0A(idx);
    
    vzeroTimeA += fESDVZERO->GetTimeV0A(idx);
    vzeroTimeC += fESDVZERO->GetTimeV0C(idx);
    
    vzeroAdcA += fESDVZERO->GetAdcV0A(idx);
    vzeroAdcC += fESDVZERO->GetAdcV0C(idx);
  }

  fVzeroMultFlagged = fVzeroMultFlaggedA + fVzeroMultFlaggedC;

  Double_t vzeroTime     = vzeroTimeA + vzeroTimeC;
  Double_t vzeroAdc      = vzeroAdcA + vzeroAdcC;
  
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroFlaggedMult")))->Fill(fVzeroMultFlagged);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroFlaggedMultA")))->Fill(fVzeroMultFlaggedA);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroFlaggedMultC")))->Fill(fVzeroMultFlaggedC);
  
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroTime")))->Fill(vzeroTime);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroTimeA")))->Fill(vzeroTimeA);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroTimeC")))->Fill(vzeroTimeC);
  
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroADC")))->Fill(vzeroAdc);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroADCA")))->Fill(vzeroAdcA);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroADCC")))->Fill(vzeroAdcC);

  // -- VZERO - TPC correlations
  if (fESDEvent->GetNumberOfTracks() > 0) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroNch")))->Fill(fVzeroMult, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroANch")))->Fill(fVzeroMultA, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroCNch")))->Fill(fVzeroMultC, fTpcTracksA);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroFNch")))->Fill(fVzeroMultFlagged, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroFANch")))->Fill(fVzeroMultFlaggedA, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroFCNch")))->Fill(fVzeroMultFlaggedC, fTpcTracksA);
  }

  return iResult;
}

//##################################################################################
Int_t AliHLTMultiplicityCorrelations::ProcessZDC() {

  Int_t iResult = 0 ;

  Double_t zdcEznA = fESDZDC->GetZDCN2Energy() / 1000.;
  Double_t zdcEznC = fESDZDC->GetZDCN1Energy() / 1000.;
  Double_t zdcEzn  = zdcEznA + zdcEznC;
  
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzn")))->Fill(zdcEzn);
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEznA")))->Fill(zdcEznA);
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEznC")))->Fill(zdcEznC);
  
  Double_t zdcEzpA = fESDZDC->GetZDCP2Energy() / 1000.;
  Double_t zdcEzpC = fESDZDC->GetZDCP1Energy() / 1000.;
  Double_t zdcEzp  = zdcEzpA + zdcEzpC;
  
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzp")))->Fill(zdcEzp);
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzpA")))->Fill(zdcEzpA);
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzpC")))->Fill(zdcEzpC);
  
  Double_t zdcEA = (zdcEznA + zdcEzpA);
  Double_t zdcEC = (zdcEznC + zdcEzpC);
  Double_t zdcE = zdcEA + zdcEC;
  
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzdc")))->Fill(zdcE);
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzdcA")))->Fill(zdcEA);
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzdcC")))->Fill(zdcEC);
  
  Double_t zdcEzemA = fESDZDC->GetZDCEMEnergy(1) / 1000.;
  Double_t zdcEzemC = fESDZDC->GetZDCEMEnergy(0) / 1000.;
  Double_t zdcEzem  = zdcEzemA + zdcEzemC;
  
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzem")))->Fill(zdcEzem);
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzemA")))->Fill(zdcEzemA);
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzemC")))->Fill(zdcEzemC);
 
  (static_cast<TH2F*>(fHistList->FindObject("fZdcEzemEzdc")))->Fill(zdcEzem, zdcE);
  (static_cast<TH2F*>(fHistList->FindObject("fZdcEzemEzdcA")))->Fill(zdcEzemA, zdcEA);
  (static_cast<TH2F*>(fHistList->FindObject("fZdcEzemEzdcC")))->Fill(zdcEzemC, zdcEC);
 
  (static_cast<TH1F*>(fHistList->FindObject("fZdcNpart")))->Fill(fESDZDC->GetZDCParticipants());
  (static_cast<TH1F*>(fHistList->FindObject("fZdcNpartA")))->Fill(fESDZDC->GetZDCPartSideA());
  (static_cast<TH1F*>(fHistList->FindObject("fZdcNpartC")))->Fill(fESDZDC->GetZDCPartSideC());
  
  (static_cast<TH1F*>(fHistList->FindObject("fZdcB")))->Fill(fESDZDC->GetImpactParameter());
  (static_cast<TH1F*>(fHistList->FindObject("fZdcBA")))->Fill(fESDZDC->GetImpactParamSideA());
  (static_cast<TH1F*>(fHistList->FindObject("fZdcBC")))->Fill(fESDZDC->GetImpactParamSideC());
 
  // -- ZDC - TPC correlations
  if (fESDEvent->GetNumberOfTracks() > 0) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcNch")))->Fill(zdcE, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcANch")))->Fill(zdcEA, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcCNch")))->Fill(zdcEC, fTpcTracksA);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzemNch")))->Fill(zdcEzem, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzemANch")))->Fill(zdcEzemA, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzemCNch")))->Fill(zdcEzemC, fTpcTracksA);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzpNch")))->Fill(zdcEzp, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzpANch")))->Fill(zdcEzpA, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzpCNch")))->Fill(zdcEzpC, fTpcTracksA);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEznNch")))->Fill(zdcEzn, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEznANch")))->Fill(zdcEznA, fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEznCNch")))->Fill(zdcEznC, fTpcTracksA);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcNpartNch")))->Fill(fESDZDC->GetZDCParticipants(), fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcNpartANch")))->Fill(fESDZDC->GetZDCPartSideA(),   fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcNpartCNch")))->Fill(fESDZDC->GetZDCPartSideC(),   fTpcTracksA);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcbNch")))->Fill(fESDZDC->GetImpactParameter(),   fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcbANch")))->Fill(fESDZDC->GetImpactParamSideA(), fTpcTracksA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcbCNch")))->Fill(fESDZDC->GetImpactParamSideC(), fTpcTracksA);
  }

  // -- VZERO - TPC correlations
  if (fESDVZERO) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcVzero")))->Fill(zdcE, fVzeroMult);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcVzeroA")))->Fill(zdcEA, fVzeroMultA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcVzeroC")))->Fill(zdcEC, fVzeroMultC);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzemVzero")))->Fill(zdcEzem, fVzeroMult);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzemVzeroA")))->Fill(zdcEzemA, fVzeroMultA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzemVzeroC")))->Fill(zdcEzemC, fVzeroMultC);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzpVzero")))->Fill(zdcEzp, fVzeroMult);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzpVzeroA")))->Fill(zdcEzpA, fVzeroMultA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzpVzeroC")))->Fill(zdcEzpC, fVzeroMultC);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEznVzero")))->Fill(zdcEzn, fVzeroMult);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEznVzeroA")))->Fill(zdcEznA, fVzeroMultA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEznVzeroC")))->Fill(zdcEznC, fVzeroMultC);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcNpartVzero")))->Fill(fESDZDC->GetZDCParticipants(), fVzeroMult);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcNpartVzeroA")))->Fill(fESDZDC->GetZDCPartSideA(),   fVzeroMultA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcNpartVzeroC")))->Fill(fESDZDC->GetZDCPartSideC(),   fVzeroMultC);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcbVzero")))->Fill(fESDZDC->GetImpactParameter(),   fVzeroMult);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcbVzeroA")))->Fill(fESDZDC->GetImpactParamSideA(), fVzeroMultA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcbVzeroC")))->Fill(fESDZDC->GetImpactParamSideC(), fVzeroMultC);

    // -- -- -- -- 

    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcVzeroF")))->Fill(zdcE, fVzeroMultFlagged);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcVzeroFA")))->Fill(zdcEA, fVzeroMultFlaggedA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcVzeroFC")))->Fill(zdcEC, fVzeroMultFlaggedC);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzemVzeroF")))->Fill(zdcEzem, fVzeroMultFlagged);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzemVzeroFA")))->Fill(zdcEzemA, fVzeroMultFlaggedA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzemVzeroFC")))->Fill(zdcEzemC, fVzeroMultFlaggedC);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzpVzeroF")))->Fill(zdcEzp, fVzeroMultFlagged);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzpVzeroFA")))->Fill(zdcEzpA, fVzeroMultFlaggedA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzpVzeroFC")))->Fill(zdcEzpC, fVzeroMultFlaggedC);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEznVzeroF")))->Fill(zdcEzn, fVzeroMultFlagged);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEznVzeroFA")))->Fill(zdcEznA, fVzeroMultFlaggedA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEznVzeroFC")))->Fill(zdcEznC, fVzeroMultFlaggedC);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcNpartVzero")))->Fill(fESDZDC->GetZDCParticipants(), fVzeroMultFlagged);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcNpartVzeroFA")))->Fill(fESDZDC->GetZDCPartSideA(),   fVzeroMultFlaggedA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcNpartVzeroFC")))->Fill(fESDZDC->GetZDCPartSideC(),   fVzeroMultFlaggedC);
    
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcbVzeroF")))->Fill(fESDZDC->GetImpactParameter(),   fVzeroMultFlagged);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcbVzeroFA")))->Fill(fESDZDC->GetImpactParamSideA(), fVzeroMultFlaggedA);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrZdcbVzeroFC")))->Fill(fESDZDC->GetImpactParamSideC(), fVzeroMultFlaggedC);
  }

  return iResult;
}
