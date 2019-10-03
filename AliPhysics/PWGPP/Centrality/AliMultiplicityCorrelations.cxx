//-*- Mode: C++ -*-

#if __GNUC__>= 3
   using namespace std;
#endif

#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"

#include "AliMultiplicityCorrelations.h"

// Task for HI Multiplicity correlation checks
// Author: Jochen Thaeder <jochen@thaeder.de>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliMultiplicityCorrelations)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliMultiplicityCorrelations::AliMultiplicityCorrelations() :
  TNamed(),
  fHistList(NULL),
  fIsMC(kFALSE),
  fESDEvent(NULL),
  fESDZDC(NULL), fESDVZERO(NULL), fESDMultiplicity(NULL),
  fESDTrackCuts(NULL), fESDTrackCuts2(NULL),
  fCleanSample(kFALSE), fCleanMinKeep(0.),fCleanMaxKeep(999999999.),
  fRunNo(-1), fCurrentRunNo(-1),
  fProcessTPC(kTRUE), fProcessSPD(kTRUE), fProcessVZERO(kTRUE), fProcessZDC(kTRUE),
  fEsdTracks(0), fEsdTracksA(0), fTpcTracks(0), fTpcTracksA(0),
  fVzeroMult(0.), fVzeroMultA(0.), fVzeroMultC(0.),
  fSpdNClusters(0), fSpdNClustersInner(0), fSpdNClustersOuter(0),
  fVzeroBinning(1500), fVzeroBinningMin(0.), fVzeroBinningMax(30000.),
  fTpcBinning(800),fTpcBinningMin(0.),fTpcBinningMax(8000.),
  fZdcBinning(1000),fZdcBinningMin(0.),fZdcBinningMax(6000.),
  fZemBinning(500),fZemBinningMin(0.),fZemBinningMax(2500.),
  fSpdBinning(750),fSpdBinningMin(0.),fSpdBinningMax(15000.) {
  // see header file for class documentation
  
}
//##################################################################################
AliMultiplicityCorrelations::AliMultiplicityCorrelations(Char_t* name, Char_t* title) : 
  TNamed(name,title), 
  fHistList(NULL),
  fIsMC(kFALSE),
  fESDEvent(NULL),
  fESDZDC(NULL), fESDVZERO(NULL), fESDMultiplicity(NULL),
  fESDTrackCuts(NULL), fESDTrackCuts2(NULL),       
  fCleanSample(kFALSE), fCleanMinKeep(0.),fCleanMaxKeep(999999999.),
  fRunNo(-1), fCurrentRunNo(-1),
  fProcessTPC(kTRUE), fProcessSPD(kTRUE),fProcessVZERO(kTRUE), fProcessZDC(kTRUE),
  fEsdTracks(0), fEsdTracksA(0), fTpcTracks(0), fTpcTracksA(0),
  fVzeroMult(0.), fVzeroMultA(0.), fVzeroMultC(0.),
  fSpdNClusters(0), fSpdNClustersInner(0), fSpdNClustersOuter(0),
  fVzeroBinning(1000), fVzeroBinningMin(0.), fVzeroBinningMax(20000.),
  fTpcBinning(800),fTpcBinningMin(0.),fTpcBinningMax(8000.),
  fZdcBinning(280),fZdcBinningMin(0.),fZdcBinningMax(140.),
  fZemBinning(100),fZemBinningMin(0.),fZemBinningMax(5.),
  fSpdBinning(750),fSpdBinningMin(0.),fSpdBinningMax(15000.) {
  // see header file for class documentation
  
}

//##################################################################################
AliMultiplicityCorrelations::~AliMultiplicityCorrelations() {
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
Int_t AliMultiplicityCorrelations::Initialize( const Char_t* listName) {
  // see header file for class documentation  

  Int_t iResult = 0;

  fHistList = new TList();
  fHistList->SetOwner(kTRUE);
  fHistList->SetName(Form("MultiplicityCorrelations_%s",listName));
  iResult = SetupHistograms();
  
  if (fProcessZDC)   { AliInfo("Processing of ZDC enabled"); }
  if (fProcessTPC)   { AliInfo("Processing of TPC enabled"); }
  if (fProcessSPD)   { AliInfo("Processing of SPD enabled"); }
  if (fProcessVZERO) { AliInfo("Processing of VZERO enabled"); }

  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Output - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliMultiplicityCorrelations::ProcessEvent(AliESDEvent *esd) {
  // see header file for class documentation  

  Int_t iResult = 0;

  if ( ! AddESDEvent(esd) ) {
    AliWarning("No ESD event.");
    return -1;
  }
  
  // -- TPC .. To be done before the others
  if (fESDEvent->GetNumberOfTracks() > 0 && fProcessTPC)
    iResult = ProcessTPC();
    
  // -- TPC / global track cuts
  if (iResult == -3)
    return iResult;
  
  // -- check if cleaning 
  if (iResult == -2)
    return iResult;
  
  // -- SPD
  if (fESDMultiplicity && fProcessSPD)
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
Bool_t AliMultiplicityCorrelations::AddESDEvent( AliESDEvent* esd ) {
  // see header file for class documentation  

  fESDEvent = esd;
  
  // -- Check for ESD
  if ( !fESDEvent ) {
    AliWarning("No ESD event present.");
    return kFALSE;
  }
    
  fESDZDC = esd->GetESDZDC();
  if ( !fESDZDC ) {
    AliInfo("No ZDC information !");
  }
  
  fESDVZERO = esd->GetVZEROData();
  if ( !fESDVZERO ) {
    AliInfo("No VZERO information !");
  }
  
  fESDMultiplicity = const_cast<AliMultiplicity*>(esd->GetMultiplicity());
  if ( !fESDMultiplicity ) {
    AliInfo("No Multiplicity information !");
  }

  if (fCurrentRunNo == esd->GetRunNumber())
    return kTRUE;
      
  fCurrentRunNo = esd->GetRunNumber();
  fRunNo = fCurrentRunNo;

  // CHANGE HERE FOR RUN RANGES
  if      ( fRunNo == 137162 ) fRunNo = 137161;
  else if ( fRunNo == 137365 ) fRunNo = 137366;
  else if ( fRunNo  > 137366 ) fRunNo = 137366;
  // CHANGE HERE FOR RUN RANGES

  return kTRUE;
}

//##################################################################################
Int_t AliMultiplicityCorrelations::SetupHistograms() {
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
Int_t AliMultiplicityCorrelations::SetupVZERO() {
  // see header file for class documentation  

  // VzeroMult
  fHistList->Add(new TH1F("fVzeroMult",  "Multiplicity^{VZERO};Multiplicity^{VZERO};N_{Events}",   
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));

  fHistList->Add(new TH1F("fVzeroMultUnCorr", "Multiplicity^{VZERO} uncorrected;Multiplicity^{VZERO};N_{Events}",   
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));

  fHistList->Add(new TH2F("fVzeroMultAC", "Multiplicity^{VZERO} A vs C;Multiplicity^{VZERO}A ;Multiplicity^{VZERO} C", 
			  fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax)); 

  return 0;
}

//##################################################################################
Int_t AliMultiplicityCorrelations::SetupZDC() {
  // see header file for class documentation  
    
  // E_{ZDC}
  fHistList->Add(new TH1F("fZdcEzdc",  "E_{ZDC};E_{ZDC} (TeV);N_{Events}",   
			  fZdcBinning,fZdcBinningMin,fZdcBinningMax));

  // E_{ZEM}
  fHistList->Add(new TH1F("fZdcEzem",  "E_{ZEM};E_{ZEM} (TeV);N_{Events}",   
			  fZemBinning,fZemBinningMin,fZemBinningMax));
  
  // E_{ZEM} vs E_{ZDC} 
  fHistList->Add(new TH2F("fZdcEzemEzdc", "E_{ZEM} vs E_{ZDC};E_{ZEM} (TeV); E_{ZDC} (TeV)",   
			  fZemBinning,fZemBinningMin,fZemBinningMax, fZdcBinning,fZdcBinningMin,fZdcBinningMax));

  return 0;
}

//##################################################################################
Int_t AliMultiplicityCorrelations::SetupTPC() {
  // see header file for class documentation  

  // Multiplicity
  fHistList->Add(new TH1F("fTpcNch0", "N_{ch} esdTracks; N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH1F("fTpcNch1", "N_{ch} accepted esdTracks; N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH1F("fTpcNch2", "N_{ch} tpcTracks; N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));
  fHistList->Add(new TH1F("fTpcNch3", "N_{ch} accepted tpcTracks; N_{ch};N_{Events}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax));

  fHistList->Add(new TH2F("fTpcCorrNch","N_{ch} accepted tpcTracks vs globalTracks; N_{tpc};N_{global}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));

  fHistList->Add(new TH1F("fTpcRatioNch","N_{ch} accepted tpcTracks/globalTracks; N_{ch};Ratio", 
			  201,0.,2.));

  fHistList->Add(new TH2F("fTpcCorrNchAll","N_{ch} accepted tpcTracks vs globalTracks - uncleaned ; N_{tpc};N_{global}", 
			  fTpcBinning,fTpcBinningMin,fTpcBinningMax, fTpcBinning,fTpcBinningMin,fTpcBinningMax));

  fHistList->Add(new TH1F("fTpcRatioNchAll","N_{ch} accepted tpcTracks/globalTracks - uncleaned; N_{ch};Ratio", 
			  201,0.,2.));

  return 0;
}

//##################################################################################
Int_t AliMultiplicityCorrelations::SetupCorrelations() {
  // see header file for class documentation  

  // ----------------------------------------------------
  //   TPC vs ZDC
  // ----------------------------------------------------
  if (fProcessTPC && fProcessZDC) {

    // N_{ch} vs E_{ZDC}
    fHistList->Add(new TH2F("fCorrEzdcNch", "N_{ch}^{TPC} vs E_{ZDC};N_{ch}^{TPC};E_{ZDC} (TeV)",   
			    fTpcBinning,fTpcBinningMin,fTpcBinningMax, fZdcBinning,fZdcBinningMin,fZdcBinningMax));


    // N_{ch} vs E_{ZEM}
    fHistList->Add(new TH2F("fCorrEzemNch", "N_{ch}^{TPC} vs E_{ZEM};N_{ch}^{TPC};E_{ZEM} (TeV)",   
			    fTpcBinning,fTpcBinningMin,fTpcBinningMax, fZdcBinning,fZdcBinningMin,fZdcBinningMax));
  }
  // ----------------------------------------------------
  //   ZDC vs VZERO
  // ----------------------------------------------------
  if (fProcessZDC && fProcessVZERO) {

    // E_{ZDC} vs Multiplicity VZERO
    fHistList->Add(new TH2F("fCorrEzdcVzero", 
			    "E_{ZDC} vs Multiplicity^{VZERO};E_{ZDC} (TeV);Multiplicity^{VZERO}",  
			    fZdcBinning,fZdcBinningMin,fZdcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));
    
    // E_{ZEM} vs Multiplicity VZERO
    fHistList->Add(new TH2F("fCorrEzemVzero", 
			    "E_{ZEM} vs Multiplicity^{VZERO};E_{ZEM} (TeV);Multiplicity^{VZERO}",  
			    fZemBinning,fZemBinningMin,fZemBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));

  }
  // ----------------------------------------------------
  //   TPC vs VZERO
  // ----------------------------------------------------
  if ( fProcessTPC && fProcessVZERO ) {

    fHistList->Add(new TH2F("fCorrVzeroNch", 
			    "N_{ch}^{TPC} vs Multiplicity^{VZERO};N_{ch}^{TPC};Multiplicity^{VZERO}", 
			    fTpcBinning,fTpcBinningMin,fTpcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));

    fHistList->Add(new TH2F("fCorrVzeroNchUnCorr", 
			    "N_{ch}^{TPC} vs Multiplicity^{VZERO} uncorrected;N_{ch}^{TPC};Multiplicity^{VZERO}", 
			    fTpcBinning,fTpcBinningMin,fTpcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));

    fHistList->Add(new TH2F("fCorrVzeroESDNch", 
			    "N_{ch}^{global} vs Multiplicity^{VZERO};N_{ch}^{TPC};Multiplicity^{VZERO}", 
			    fTpcBinning,fTpcBinningMin,fTpcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));

    fHistList->Add(new TH2F("fCorrVzeroESDNchUnCorr", 
			    "N_{ch}^{global} vs Multiplicity^{VZERO} uncorrected;N_{ch}^{TPC};Multiplicity^{VZERO}", 
			    fTpcBinning,fTpcBinningMin,fTpcBinningMax, fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax));

  }
  // ----------------------------------------------------
  //   TPC vs SPD
  // ----------------------------------------------------
  if ( fProcessSPD && fProcessTPC ) {

    fHistList->Add(new TH2F("fCorrSpdTpcNch", "N_{ch}^{TPC} vs N_{clusters}^{SPD};N_{ch}^{TPC};N_{clusters}^{SPD}",   
			    fTpcBinning,fTpcBinningMin,fTpcBinningMax, fSpdBinning,fSpdBinningMin,fSpdBinningMax));  

    fHistList->Add(new TH2F("fCorrSpdOuterTpcNch"," N_{ch}^{TPC} vs N_{clusters}^{SPD}_{Outer};N_{ch}^{TPC};N_{clusters}^{SPD}",   
			    fTpcBinning,fTpcBinningMin,fTpcBinningMax, fSpdBinning,fSpdBinningMin,fSpdBinningMax));  
  }
  // ----------------------------------------------------
  //   SPD vs VZERO
  // ----------------------------------------------------
  if ( fProcessSPD && fProcessVZERO ) {

    fHistList->Add(new TH2F("fCorrVzeroSpd", 
			    "Multiplicity^{VZERO} vs N_{ch}^{SPD};Multiplicity^{VZERO};N^{SPD}", 
			    fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fSpdBinning,fSpdBinningMin,fSpdBinningMax));
    
    fHistList->Add(new TH2F("fCorrVzeroSpdOuter", 
			    "Multiplicity^{VZERO} vs N_{ch}^{SPD}_{Outer};Multiplicity^{VZERO};N^{SPD}", 
			    fVzeroBinning,fVzeroBinningMin,fVzeroBinningMax, fSpdBinning,fSpdBinningMin,fSpdBinningMax));
  }
  // ----------------------------------------------------
  //   SPD vs ZDC
  // ----------------------------------------------------
  if ( fProcessSPD && fProcessZDC) {

    // E_{ZDC} vs Multiplicity SPD
    fHistList->Add(new TH2F("fCorrEzdcSpd", "E_{ZDC} vs N_{ch}^{SPD};E_{ZDC} (TeV);N^{SPD}",   
			    fZdcBinning,fZdcBinningMin,fZdcBinningMax, fSpdBinning,fSpdBinningMin,fSpdBinningMax));

    fHistList->Add(new TH2F("fCorrEzdcSpdOuter", "E_{ZDC} vs N_{ch}^{SPD};E_{ZDC} (TeV);N^{SPD}",   
			    fZdcBinning,fZdcBinningMin,fZdcBinningMax, fSpdBinning,fSpdBinningMin,fSpdBinningMax));
    
  }

  return 0;
}

//##################################################################################
Int_t AliMultiplicityCorrelations::SetupSPD() {
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
Int_t AliMultiplicityCorrelations::ProcessTPC() {
  // see header file for class documentation  

  Int_t iResult = 0 ;

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
    
    if (!fESDTrackCuts2->AcceptTrack(esdTrack)) continue;
    ++fEsdTracksA;
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

    if (!fESDTrackCuts->AcceptTrack(esdTrack)) continue;
    ++fTpcTracksA;
  }

  // -----------------------------------------------------------------

  if (fCleanSample) {
    if(fTpcTracksA < fCleanMinKeep || fTpcTracksA > fCleanMaxKeep)
      return -2;
  }

  Float_t nESD  = (Float_t) fEsdTracksA;
  Float_t nTPC  = (Float_t) fTpcTracksA;
  Float_t ratio = -1;

  if ( nESD != 0 )
    ratio = nTPC/nESD;

  (static_cast<TH2F*>(fHistList->FindObject("fTpcCorrNchAll")))->Fill(fTpcTracksA,fEsdTracksA);  
  (static_cast<TH1F*>(fHistList->FindObject("fTpcRatioNchAll")))->Fill(ratio);

  // Cleaning Cut tpcTracks and globalTracks
  if ( fEsdTracksA < (-39+0.5797*fTpcTracksA))
    return -3;

  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch0")))->Fill(fEsdTracks);
  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch1")))->Fill(fEsdTracksA);
  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch2")))->Fill(fTpcTracks);
  (static_cast<TH1F*>(fHistList->FindObject("fTpcNch3")))->Fill(fTpcTracksA);

  (static_cast<TH2F*>(fHistList->FindObject("fTpcCorrNch")))->Fill(fTpcTracksA,fEsdTracksA);
  (static_cast<TH1F*>(fHistList->FindObject("fTpcRatioNch")))->Fill(ratio);

  return iResult;
}

//##################################################################################
Int_t AliMultiplicityCorrelations::ProcessVZERO() {
  // see header file for class documentation  

  Int_t iResult = 0 ;
  Float_t vzeroCorr;;
  fVzeroMultA = fESDVZERO->GetMTotV0A();
  fVzeroMultC = fESDVZERO->GetMTotV0C();
  if (!fIsMC)
    fVzeroMult  = GetCorrVZERO(vzeroCorr);
  else
    fVzeroMult  = (fVzeroMultA+fVzeroMultC) * 0.85871;

  (static_cast<TH1F*>(fHistList->FindObject("fVzeroMult")))->Fill(fVzeroMult);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroMultUnCorr")))->Fill(fVzeroMultA+fVzeroMultC);
  (static_cast<TH1F*>(fHistList->FindObject("fVzeroMultAC")))->Fill(fVzeroMultA,fVzeroMultC);

  // -- VZERO - TPC correlations
  if (fESDEvent->GetNumberOfTracks() > 0 && fProcessTPC) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroNch")))->Fill(fTpcTracksA,fVzeroMult);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroNchUnCorr")))->Fill(fTpcTracksA,fVzeroMultA+fVzeroMultC);

    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroESDNch")))->Fill(fEsdTracksA,fVzeroMult);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroESDNchUnCorr")))->Fill(fEsdTracksA,fVzeroMultA+fVzeroMultC);
  }

  // -- VZERO - SPD correlations
  if (fESDEvent->GetNumberOfTracks() > 0 && fProcessSPD) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroSpd")))->Fill(fVzeroMult, fSpdNClusters);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrVzeroSpdOuter")))->Fill(fVzeroMult, fSpdNClustersOuter);
  }
  
  return iResult;
}

//##################################################################################
Int_t AliMultiplicityCorrelations::ProcessZDC() {
  // see header file for class documentation  

  Int_t iResult = 0 ;

  Double_t zdcEznA = fESDZDC->GetZDCN2Energy() / 8.;
  Double_t zdcEznC = fESDZDC->GetZDCN1Energy() / 8.;
  //  Double_t zdcEznA = fESDZDC->GetZDCN2Energy() / 1000.;
  //  Double_t zdcEznC = fESDZDC->GetZDCN1Energy() / 1000.;
  
  Double_t zdcEzpA = fESDZDC->GetZDCP2Energy() / 8.;
  Double_t zdcEzpC = fESDZDC->GetZDCP1Energy() / 8.;
  //  Double_t zdcEzpA = fESDZDC->GetZDCP2Energy() / 1000.;
  //  Double_t zdcEzpC = fESDZDC->GetZDCP1Energy() / 1000.;
  
  Double_t zdcEA = (zdcEznA + zdcEzpA);
  Double_t zdcEC = (zdcEznC + zdcEzpC);
  Double_t zdcE = zdcEA + zdcEC;
  
  Double_t zdcEzemA = fESDZDC->GetZDCEMEnergy(1) / 8.;
  Double_t zdcEzemC = fESDZDC->GetZDCEMEnergy(0) / 8.;
  //  Double_t zdcEzemA = fESDZDC->GetZDCEMEnergy(1) / 1000.;
  //  Double_t zdcEzemC = fESDZDC->GetZDCEMEnergy(0) / 1000.;
  Double_t zdcEzem  = zdcEzemA + zdcEzemC;

  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzdc")))->Fill(zdcE);
  
  (static_cast<TH1F*>(fHistList->FindObject("fZdcEzem")))->Fill(zdcEzem);
  (static_cast<TH2F*>(fHistList->FindObject("fZdcEzemEzdc")))->Fill(zdcEzem, zdcE);
  
  // -- ZDC - TPC correlations
  if (fESDEvent->GetNumberOfTracks() > 0 && fProcessTPC) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcNch")))->Fill(fTpcTracksA, zdcE);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzemNch")))->Fill(fTpcTracksA, zdcEzem);
  }
  
  // -- ZDC - SPD correlations
  if (fProcessSPD) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcSpd")))->Fill(zdcE, fSpdNClusters);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcSpdOuter")))->Fill(zdcE, fSpdNClustersOuter);
  }

  // -- VZERO - ZDC correlations
  if (fESDVZERO && fProcessVZERO) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzdcVzero")))->Fill(zdcE, fVzeroMult);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrEzemVzero")))->Fill(zdcEzem, fVzeroMult);
  }

  return iResult;
}

//##################################################################################
Int_t AliMultiplicityCorrelations::ProcessSPD() {
  // see header file for class documentation
  
   const AliESDVertex* vtxESD = fESDEvent->GetPrimaryVertexSPD();
   Float_t zvtx = vtxESD->GetZ();

  if (!fIsMC) {
    fSpdNClustersInner = GetCorrSPD2(fESDMultiplicity->GetNumberOfITSClusters(0),zvtx);
    fSpdNClustersOuter = GetCorrSPD2(fESDMultiplicity->GetNumberOfITSClusters(1),zvtx);
  }
  else {
    fSpdNClustersInner = fESDMultiplicity->GetNumberOfITSClusters(0);
    fSpdNClustersOuter = fESDMultiplicity->GetNumberOfITSClusters(1);
  }    

  fSpdNClusters      = fSpdNClustersOuter + fSpdNClustersInner;

  (static_cast<TH1F*>(fHistList->FindObject("fSpdNClusters")))->Fill(fSpdNClusters);
  (static_cast<TH1F*>(fHistList->FindObject("fSpdNClustersInner")))->Fill(fSpdNClustersInner);
  (static_cast<TH1F*>(fHistList->FindObject("fSpdNClustersOuter")))->Fill(fSpdNClustersOuter);
  
  // -- SPD vs TPC correlations
  if (fProcessTPC) {
    (static_cast<TH2F*>(fHistList->FindObject("fCorrSpdTpcNch")))->Fill(fTpcTracksA, fSpdNClusters);
    (static_cast<TH2F*>(fHistList->FindObject("fCorrSpdOuterTpcNch")))->Fill(fTpcTracksA, fSpdNClustersOuter);
  }

  return 0;
}

//##################################################################################
Float_t AliMultiplicityCorrelations::GetCorrVZERO(Float_t &v0CorrResc) {
  // correct V0 non-linearity, prepare a version rescaled to SPD2 corr

  Double_t *par0;
  Double_t *par1;
  Double_t *par2;
  
  Double_t par0R137161[64] = { 6.71e-02 , 6.86e-02 , 7.06e-02 , 6.32e-02 , 
			       5.91e-02 , 6.07e-02 , 5.78e-02 , 5.73e-02 , 5.91e-02 , 6.22e-02 , 
			       5.90e-02 , 6.11e-02 , 5.55e-02 , 5.29e-02 , 5.19e-02 , 5.56e-02 , 
			       6.25e-02 , 7.03e-02 , 5.64e-02 , 5.81e-02 , 4.57e-02 , 5.30e-02 , 
			       5.13e-02 , 6.43e-02 , 6.27e-02 , 6.48e-02 , 6.07e-02 , 1.01e-01 , 
			       6.68e-02 , 7.16e-02 , 6.36e-02 , 5.95e-02 , 2.52e-02 , 2.82e-02 , 
			       2.56e-02 , 2.86e-02 , 2.82e-02 , 2.10e-02 , 2.13e-02 , 2.32e-02 , 
			       2.75e-02 , 4.34e-02 , 3.78e-02 , 4.52e-02 , 4.11e-02 , 3.89e-02 , 
			       4.10e-02 , 3.73e-02 , 4.51e-02 , 5.07e-02 , 5.42e-02 , 4.74e-02 , 
			       4.33e-02 , 4.44e-02 , 4.64e-02 , 3.01e-02 , 6.38e-02 , 5.26e-02 , 
			       4.99e-02 , 5.26e-02 , 5.47e-02 , 3.84e-02 , 5.00e-02 , 5.20e-02 };
  Double_t par1R137161[64] = { -6.68e-05 , -7.78e-05 , -6.88e-05 , -5.92e-05 , 
			       -2.43e-05 , -3.54e-05 , -2.91e-05 , -1.99e-05 , -1.40e-05 , -4.01e-05 , 
			       -2.29e-05 , -3.68e-05 , -2.53e-05 , -2.44e-06 , -9.22e-06 , -1.51e-05 , 
			       -2.80e-05 , -2.34e-05 , -1.72e-05 , -1.81e-05 , -1.29e-05 , -2.65e-05 , 
			       -1.61e-05 , -2.86e-05 , -1.74e-05 , -4.23e-05 , -3.41e-05 , -1.05e-04 , 
			       -2.76e-05 , -4.71e-05 , -3.06e-05 , -2.32e-05 , -1.55e-06 , 2.15e-05 , 
			       1.40e-05 , 2.16e-05 , 1.21e-05 , 3.05e-06 , 1.67e-05 , -3.84e-06 , 
			       3.09e-06 , 1.50e-05 , 3.47e-06 , 4.87e-06 , -3.71e-07 , -1.75e-06 , 
			       -1.80e-06 , 9.99e-06 , -6.46e-06 , -4.91e-06 , 1.33e-05 , -2.52e-07 , 
			       -3.85e-06 , 4.94e-06 , -2.48e-07 , -1.20e-05 , 2.07e-06 , 6.12e-06 , 
			       -1.18e-06 , 4.54e-06 , -1.54e-05 , -1.25e-05 , 1.46e-06 , -6.67e-06 };
  Double_t par2R137161[64] = { 1.29e-08 , 1.51e-08 , 1.43e-08 , 1.11e-08 , 
			       5.04e-09 , 6.99e-09 , 5.58e-09 , 4.15e-09 , 4.00e-09 , 8.22e-09 , 
			       4.97e-09 , 7.66e-09 , 4.91e-09 , 1.10e-09 , 2.64e-09 , 3.64e-09 , 
			       5.76e-09 , 5.46e-09 , 3.38e-09 , 3.47e-09 , 2.43e-09 , 4.13e-09 , 
			       2.80e-09 , 5.80e-09 , 3.86e-09 , 7.46e-09 , 5.98e-09 , 2.58e-08 , 
			       5.50e-09 , 8.72e-09 , 5.23e-09 , 4.37e-09 , 2.33e-09 , -6.01e-10 , 
			       3.99e-11 , -2.02e-10 , 7.67e-10 , 2.03e-09 , 1.17e-10 , 2.56e-09 , 
			       1.16e-09 , -4.75e-10 , 1.28e-09 , 1.23e-09 , 1.62e-09 , 1.61e-09 , 
			       1.93e-09 , 2.97e-10 , 2.21e-09 , 2.16e-09 , 5.22e-10 , 1.03e-09 , 
			       1.56e-09 , 5.00e-10 , 1.01e-09 , 2.93e-09 , 1.05e-09 , 9.96e-11 , 
			       1.21e-09 , 7.45e-10 , 3.07e-09 , 2.31e-09 , 6.70e-10 , 1.89e-09 };
  
  Double_t par0R137366[64] = { 7.12e-02 , 7.34e-02 , 7.39e-02 , 6.54e-02 , 6.11e-02 , 6.31e-02 , 6.15e-02 , 
			       6.00e-02 , 6.10e-02 , 6.49e-02 , 6.17e-02 , 6.33e-02 , 6.00e-02 , 5.48e-02 , 
			       5.44e-02 , 5.81e-02 , 6.49e-02 , 7.07e-02 , 5.91e-02 , 6.18e-02 , 4.82e-02 , 
			       5.67e-02 , 5.36e-02 , 6.60e-02 , 6.37e-02 , 6.78e-02 , 6.31e-02 , 1.04e-01 , 
			       6.91e-02 , 7.32e-02 , 6.61e-02 , 6.16e-02 , 2.64e-02 , 2.81e-02 , 2.64e-02 , 
			       2.85e-02 , 2.87e-02 , 2.18e-02 , 2.19e-02 , 2.43e-02 , 2.81e-02 , 4.37e-02 , 
			       3.90e-02 , 4.66e-02 , 4.24e-02 , 4.09e-02 , 4.21e-02 , 3.88e-02 , 4.83e-02 , 
			       5.23e-02 , 5.44e-02 , 4.85e-02 , 4.42e-02 , 4.58e-02 , 4.74e-02 , 3.14e-02 , 
			       6.31e-02 , 5.30e-02 , 5.01e-02 , 5.33e-02 , 5.70e-02 , 3.95e-02 , 4.98e-02 , 5.31e-02 };
  Double_t par1R137366[64] = { -6.99e-05 , -6.99e-05 , -6.94e-05 , -6.55e-05 , -3.55e-05 , -4.50e-05 , 
			       -3.10e-05 , -2.81e-05 , -2.29e-05 , -3.89e-05 , -2.53e-05 , -4.25e-05 ,
			       -1.87e-05 , -2.01e-05 , -1.53e-05 , -2.14e-05 , -2.86e-05 , -4.70e-05 ,
			       -2.23e-05 , -3.30e-05 ,-9.74e-06 , -2.62e-05 , -1.76e-05 , -2.38e-05 , 
			       -2.40e-05 , -3.43e-05 , -2.75e-05 , -6.86e-05 ,-2.35e-05 , -4.45e-05 , 
			       -2.51e-05 , -2.20e-05 , -1.25e-16 , -2.04e-17 , -2.06e-17 , -3.74e-19 ,
			       -1.18e-18 , -2.02e-15 , -3.78e-06 , -1.26e-06 , -2.71e-06 , -6.23e-17 , 
			       -7.39e-08 , -1.76e-16 , -8.98e-06 , -4.10e-18 , -1.34e-05 , -1.06e-16 , 
			       -3.34e-06 , -1.04e-05 , -5.28e-06 , -7.34e-06 , -1.05e-05 , -7.68e-06 ,
			       -1.78e-05 , -1.19e-05 , -1.78e-05 , -1.34e-06 , -9.23e-06 , -3.34e-06 ,
			       -8.02e-06 , -1.39e-05 , -1.38e-05 , -1.40e-05 };
  Double_t par2R137366[64] = { 1.41e-08 , 1.47e-08 , 1.48e-08 , 1.24e-08 , 6.82e-09 , 8.73e-09 , 6.26e-09 , 
			       5.53e-09 , 5.40e-09 , 7.93e-09 , 5.49e-09 , 8.77e-09 , 4.21e-09 , 3.93e-09 , 
			       3.60e-09 , 4.67e-09 , 5.59e-09 , 8.81e-09 , 3.89e-09 , 6.19e-09 , 1.97e-09 , 
			       4.38e-09 , 3.26e-09 , 5.00e-09 , 4.58e-09 , 6.39e-09 , 5.03e-09 , 1.30e-08 , 
			       4.95e-09 , 8.26e-09 , 4.57e-09 , 4.10e-09 , 2.35e-09 , 2.30e-09 , 2.15e-09 , 
			       2.27e-09 , 2.17e-09 , 2.27e-09 , 2.97e-09 , 2.25e-09 , 1.69e-09 , 1.44e-09 , 
			       1.66e-09 , 1.75e-09 , 2.88e-09 , 1.82e-09 , 3.64e-09 , 1.80e-09 , 1.71e-09 , 
			       2.66e-09 , 3.01e-09 , 1.95e-09 , 2.64e-09 , 2.42e-09 , 3.68e-09 , 2.66e-09 , 
			       3.92e-09 , 1.18e-09 , 2.26e-09 , 1.57e-09 , 2.02e-09 , 2.71e-09 , 2.99e-09 , 3.04e-09 }; 
  
  

  if ( fRunNo == 137161 ) {
    par0=par0R137161;
    par1=par1R137161;
    par2=par2R137161;
  }  
  else {
    par0=par0R137366;
    par1=par1R137366;
    par2=par2R137366;
  }

  Float_t multCorr = 0;
  Float_t multCorr2 = 0;
  Float_t multChCorr[64];

  for(Int_t i = 0; i < 64; ++i) {
    Double_t b = (fESDVZERO->GetMultiplicity(i)*par1[i]-par0[i]);
    Double_t s = (b*b-4.*par2[i]*fESDVZERO->GetMultiplicity(i)*fESDVZERO->GetMultiplicity(i));
    Double_t n;
    if (s<0) {
      printf("FPE %d %.2f %.2f %.2e\n",i,fESDVZERO->GetMultiplicity(i),b,(b*b-4.*par2[i]*fESDVZERO->GetMultiplicity(i)*fESDVZERO->GetMultiplicity(i)));
      n = -b;
    }
    else {
      n = (-b + TMath::Sqrt(s));
    }
    multChCorr[i] = 2.*fESDVZERO->GetMultiplicity(i)/n*par0[i];
    multCorr += multChCorr[i];
    multCorr2 += (multChCorr[i]/par0[i]/64.);
  }
  v0CorrResc = multCorr2;
  return multCorr;
}


//##################################################################################
Float_t AliMultiplicityCorrelations::GetCorrSPD2(Float_t spd2raw,Float_t zv) const {
  // renormalize N spd2 clusters at given Zv to acceptance at Zv=0
  const double pars[] = {8.10030e-01,-2.80364e-03,-7.19504e-04};
  zv -= pars[0];
  float corr = 1 + zv*(pars[1] + zv*pars[2]);
  return corr>0 ? spd2raw/corr : -1;
}
