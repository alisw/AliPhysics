#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

#include "AliOADBContainer.h"
#include "AliOADBMuonTrackCutsParam.h"
#endif

// Needed libraries:
// gSystem->Load("libANALYSIS.so");gSystem->Load("libOADB.so");gSystem->Load("libANALYSISalice.so");gSystem->Load("libCORRFW.so");gSystem->Load("libPWGmuon.so");

//________________________________________________________
AliOADBMuonTrackCutsParam* CreateOADBObject ( TString periodName, Double_t meanDcaX, Double_t meanDcaY, Double_t meanDcaZ, Double_t meanPCorr23, Double_t meanPCorr310, Double_t sigmaPdca23, Double_t sigmaPdca310, Double_t nSigmaPdca, Double_t chi2NormCut, Double_t relPResolution, Double_t slopeResolution, Double_t sharpPtApt, Double_t sharpPtLpt, Double_t sharpPtHpt )
{
  AliOADBMuonTrackCutsParam* oadbObj = new AliOADBMuonTrackCutsParam ( periodName.Data() );
  oadbObj->SetMeanDCA ( meanDcaX, meanDcaY, meanDcaZ );
  oadbObj->SetMeanPCorr ( meanPCorr23, meanPCorr310 );
  oadbObj->SetSigmaPdca ( sigmaPdca23, sigmaPdca310 );
  oadbObj->SetNSigmaPdca ( nSigmaPdca );
  oadbObj->SetChi2NormCut ( chi2NormCut );
  oadbObj->SetRelPResolution ( relPResolution );
  oadbObj->SetSlopeResolution ( slopeResolution );
  oadbObj->SetSharpPtCut ( sharpPtApt, sharpPtLpt, sharpPtHpt );

  return oadbObj;
}

//________________________________________________________
Bool_t HasDefault ( AliOADBContainer* container )
{
  return ( container->GetDefaultList()->GetEntries() > 0 );
}

//________________________________________________________
void AddParams ( AliOADBContainer* container, AliOADBMuonTrackCutsParam* oadbObj, Int_t firstRun, Int_t lastRun, TString passName, AliOADBContainer* defaultContainer = 0x0 )
{
  TString oadbObjName = oadbObj->GetName();
  if ( ! oadbObjName.Contains(passName.Data()) ) oadbObjName += "_" + passName;
  oadbObj->SetName(oadbObjName.Data());
  container->AppendObject ( oadbObj, firstRun, lastRun, passName );
  if ( ! defaultContainer ) return;
  defaultContainer->AppendObject ( oadbObj->Clone(), firstRun, lastRun );
}

enum {kData, kDataDef, kMC, kMCDef, kNcontainers};
AliOADBContainer* containers[kNcontainers] = {0x0};

//________________________________________________________
void FillDataAndMC ( AliOADBMuonTrackCutsParam* oadbObj, Int_t firstRun, Int_t lastRun, TString passName, Bool_t isDefault, Bool_t resetDCAforMC = kTRUE )
{
  AliOADBContainer* defContainer = ( isDefault ) ? containers[kDataDef] : 0x0;
  AddParams ( containers[kData], oadbObj, firstRun, lastRun, passName, defContainer );
  AliOADBMuonTrackCutsParam* oadbObjMC = static_cast<AliOADBMuonTrackCutsParam*>(oadbObj->Clone());
  if ( resetDCAforMC ) oadbObjMC->SetMeanDCA ( 0., 0., 0. );
  defContainer = ( isDefault ) ? containers[kMCDef] : 0x0;
  AddParams ( containers[kMC], oadbObjMC, firstRun, lastRun, passName, defContainer );
}

//________________________________________________________
void buildMuonTrackCutsOADB ( )
{
  Double_t defMeanPCorr23 = 2.*1.5;
  Double_t defMeanPCorr310 = 2.*1.2;
  Double_t defSigmaPdca23 = 99.;
  Double_t defSigmaPdca310 = 54.;
  Double_t defNSigmaPdca = 6.;
  Double_t defChi2NormCut = 1.e6;
  Double_t lut4 = 4.2;
  Double_t lut2 = 1.7;
  Double_t lut1 = 1.0;
  Double_t lut05 = 0.5;

  AliOADBMuonTrackCutsParam* oadbObj = 0x0;

  TString baseContName = "MuonTrackCutsParam";
  TString contNameSuffix[kNcontainers] = {"data","data_def","MC","MC_def"};


  for ( Int_t icont=0; icont<kNcontainers; icont++ ) {
    containers[icont] = new AliOADBContainer(Form("%s_%s",baseContName.Data(),contNameSuffix[icont].Data()));
    oadbObj = CreateOADBObject ( "default",
                                0., 0., 0.,
                                defMeanPCorr23, defMeanPCorr310,
                                defSigmaPdca23, defSigmaPdca310, defNSigmaPdca,
                                defChi2NormCut,
                                5.e-4, 6.e-4,
                                0., lut1, lut4 );
    containers[icont]->AddDefaultObject(oadbObj);
  }



  oadbObj = CreateOADBObject ( "LHC12hi",
                              0.18943, -0.250591, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              80., defSigmaPdca310, defNSigmaPdca,
                              defChi2NormCut,
                              4.e-4, 5.e-4,
                              0., lut1, lut4 );
  FillDataAndMC ( oadbObj, 189576, 193341, "muon_calo_pass2", kTRUE );


  oadbObj = CreateOADBObject ( "LHC13def",
                              0.18943, -0.250591, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              80., defSigmaPdca310, defNSigmaPdca,
                              defChi2NormCut,
                              4.e-4, 5.e-4,
                              0., lut05, lut4 );
  FillDataAndMC ( oadbObj, 195682, 197388, "muon_pass2", kTRUE );


  oadbObj = CreateOADBObject ( "LHC13f",
                              0.18943, -0.250591, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              80., defSigmaPdca310, defNSigmaPdca,
                              defChi2NormCut,
                              4.e-4, 5.e-4,
                              0., lut05, lut4 );
  FillDataAndMC ( oadbObj, 196433, 197388, "muon_calo", kFALSE );


  oadbObj = CreateOADBObject ( "LHC11h_2",
                              -1.118002, -1.175119, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNSigmaPdca,
                              defChi2NormCut,
                              5.e-4, 6.e-4,
                              0., lut1, lut4 );
  FillDataAndMC ( oadbObj, 167706, 170593, "pass2", kFALSE );


  oadbObj = CreateOADBObject ( "LHC11c3defgh",
                              0.296514, -0.229262, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNSigmaPdca,
                              defChi2NormCut,
                              4.e-4, 5.e-4,
                              0., lut1, lut4 );
  FillDataAndMC ( oadbObj, 154726, 170593, "pass2_muon", kTRUE );


  oadbObj = CreateOADBObject ( "LHC11c2",
                              0.296514, -0.229262, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNSigmaPdca,
                              defChi2NormCut,
                              4.e-4, 5.e-4,
                              0., lut1, lut2 );
  FillDataAndMC ( oadbObj, 153059, 154495, "pass2_muon", kTRUE );


  oadbObj = CreateOADBObject ( "LHC11c1",
                               0.296514, -0.229262, 0.,
                               defMeanPCorr23, defMeanPCorr310,
                               defSigmaPdca23, defSigmaPdca310, defNSigmaPdca,
                               defChi2NormCut,
                               4.e-4, 5.e-4,
                               0., lut05, lut1 );
  FillDataAndMC ( oadbObj, 151661, 152935, "pass2_muon", kTRUE );


  oadbObj = CreateOADBObject ( "LHC11h",
                              -1.118002, -1.175119, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNSigmaPdca,
                              defChi2NormCut,
                              5.e-4, 6.e-4,
                              0., lut1, lut4 );
  FillDataAndMC ( oadbObj, 167706, 170593, "pass1_muon", kFALSE );


  oadbObj = CreateOADBObject ( "LHC11def",
                              -1.146336, -1.130535, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNSigmaPdca,
                              defChi2NormCut,
                              5.e-4, 6.e-4,
                              0., lut1, lut4 );
  FillDataAndMC ( oadbObj, 156620, 162717, "pass1", kFALSE );


  oadbObj = CreateOADBObject ( "LHC11a",
                              -1.146336, -1.130535, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNSigmaPdca,
                              defChi2NormCut,
                              5.e-4, 6.e-4,
                              0., lut05, lut1 );
  FillDataAndMC ( oadbObj, 146688, 146860, "pass1", kTRUE );


  oadbObj = CreateOADBObject ( "LHC10h",
                              -0.4599, -0.9172, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNSigmaPdca,
                              3.5,
                              5.e-4, 6.e-4,
                              0., lut05, lut1 );
  FillDataAndMC ( oadbObj, 137135, 139513, "pass1", kTRUE );


  oadbObj = CreateOADBObject ( "LHC10pp",
                              -0.4599, -0.9172, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNSigmaPdca,
                              defChi2NormCut,
                              5.e-4, 6.e-4,
                              0., lut05, lut1 );
  FillDataAndMC ( oadbObj, 114783, 136376, "pass1", kTRUE );


  TString oadbFilename = "$ALICE_ROOT/OADB/PWG/MUON/MuonTrackCuts.root";
  if ( ! gSystem->AccessPathName(gSystem->ExpandPathName(oadbFilename.Data())) ) {
    gSystem->Exec(Form("rm %s", oadbFilename.Data()));
  }

  for ( Int_t icont=0; icont<kNcontainers; icont++ ) {
    containers[icont]
    ->WriteToFile ( oadbFilename.Data() );
  }
}
