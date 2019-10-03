#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"

#include "AliOADBContainer.h"
#include "AliOADBMuonTrackCutsParam.h"
#endif

// Needed libraries:
// gSystem->Load("libANALYSIS");gSystem->Load("libOADB");gSystem->Load("libANALYSISalice");gSystem->Load("libCORRFW");gSystem->Load("libPWGmuon");

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

//________________________________________________________
AliOADBMuonTrackCutsParam* CloneOADBobject ( TString objName, const AliOADBMuonTrackCutsParam* inputOADBobj )
{
  return static_cast<AliOADBMuonTrackCutsParam*>(inputOADBobj->Clone(objName.Data()));
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
  // These parameters never changed since the beginning
  Double_t defMeanPCorr23 = 2.*1.5;
  Double_t defMeanPCorr310 = 2.*1.2;
  Double_t defChi2NormCut = 1.e6;

  // These parameters represent all of the possible LUT values
  Double_t lut4 = 4.2;
  Double_t lut2 = 1.7;
  Double_t lut1 = 1.0;
  Double_t lut05 = 0.5;
  Double_t aptCut = 0.;

  // This parameters are unchanged since RunII
  // The alignment procedudre has become quite under control
  // Once the alignment is performed, we always recover the same resolution values
  // and therefore the same pxDCA values
  Double_t defMeanDcaX = 0.;
  Double_t defMeanDcaY = 0.;
  Double_t defMeanDcaZ = 0.;
  Double_t defSigmaPdca23 = 80.;
  Double_t defSigmaPdca310 = 54.;
  Double_t defNsigmaPdca = 6.;
  Double_t defRelPresolution = 4.e-4;
  Double_t defSlopeResolution = 5.e-4;

  AliOADBMuonTrackCutsParam* oadbObj = 0x0;
  AliOADBMuonTrackCutsParam* cloneOadbObj = 0x0;

  TString baseContName = "MuonTrackCutsParam";
  TString contNameSuffix[kNcontainers] = {"data","data_def","MC","MC_def"};


  for ( Int_t icont=0; icont<kNcontainers; icont++ ) {
    containers[icont] = new AliOADBContainer(Form("%s_%s",baseContName.Data(),contNameSuffix[icont].Data()));
    oadbObj = CreateOADBObject ( "default",
                                defMeanDcaX, defMeanDcaY, defMeanDcaZ,
                                defMeanPCorr23, defMeanPCorr310,
                                defSigmaPdca23, defSigmaPdca310, defNsigmaPdca,
                                defChi2NormCut,
                                defRelPresolution, defSlopeResolution,
                                aptCut, lut05, lut4 ); // Change this line when the LUT changes
    // So that the default values are ok for new runs
    containers[icont]->AddDefaultObject(oadbObj);
  }

  // pPb 5 TeV + pPb 8 TeV 2016
  oadbObj = CreateOADBObject ( "LHC16t",
                              defMeanDcaX, defMeanDcaY, defMeanDcaZ,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNsigmaPdca,
                              defChi2NormCut,
                              defRelPresolution, defSlopeResolution,
                              aptCut, lut05, lut4 );
  FillDataAndMC ( oadbObj, 267161, 267166, "muon_calo_pass1", kTRUE );

  cloneOadbObj = CloneOADBobject("LHC16s",oadbObj);
  FillDataAndMC ( cloneOadbObj, 266405, 267131, "muon_calo_pass2", kTRUE );

  cloneOadbObj = CloneOADBobject("LHC16s",oadbObj);
  FillDataAndMC ( cloneOadbObj, 266405, 267131, "muon_calo_pass1", kFALSE );

  cloneOadbObj = CloneOADBobject("LHC16r",oadbObj);
  FillDataAndMC ( cloneOadbObj, 265589, 266318, "muon_calo_pass1", kTRUE );

  cloneOadbObj = CloneOADBobject("LHC16q",oadbObj);
  FillDataAndMC ( cloneOadbObj, 265304, 265525, "pass1_FAST", kTRUE );


//  // pp 13 TeV 2016
//  oadbObj = CreateOADBObject ( "LHC16fghijklmnop",
//                              defMeanDcaX, defMeanDcaY, defMeanDcaZ,
//                              defMeanPCorr23, defMeanPCorr310,
//                              defSigmaPdca23, defSigmaPdca310, defNsigmaPdca,
//                              defChi2NormCut,
//                              defRelPresolution, defSlopeResolution,
//                              aptCut, lut05, lut4 );
//  FillDataAndMC ( oadbObj, 253614, 264347, "muon_calo_pass1", kTRUE );
//
//  oadbObj = CreateOADBObject ( "LHC16de",
//                              defMeanDcaX, defMeanDcaY, defMeanDcaZ,
//                              defMeanPCorr23, defMeanPCorr310,
//                              defSigmaPdca23, defSigmaPdca310, defNsigmaPdca,
//                              defChi2NormCut,
//                              defRelPresolution, defSlopeResolution,
//                              aptCut, lut1, lut4 );
//  FillDataAndMC ( oadbObj, 252235, 253603, "muon_calo_pass1", kTRUE );


  // PbPb 5 TeV 2015
  oadbObj = CreateOADBObject ( "LHC15o",
                              defMeanDcaX, defMeanDcaY, defMeanDcaZ,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNsigmaPdca,
                              defChi2NormCut,
                              defRelPresolution, defSlopeResolution,
                              aptCut, lut1, lut4 );
  FillDataAndMC ( oadbObj, 244918, 246994, "muon_calo_pass1", kTRUE );


  // pp 5 TeV 2015
  oadbObj = CreateOADBObject ( "LHC15n",
                              defMeanDcaX, defMeanDcaY, defMeanDcaZ,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNsigmaPdca,
                              defChi2NormCut,
                              defRelPresolution, defSlopeResolution,
                              aptCut, lut05, lut4 );
  FillDataAndMC ( oadbObj, 244340, 244628, "muon_calo_pass2", kTRUE );


  // pp 13 TeV 2015
  oadbObj = CreateOADBObject ( "LHC15l",
                              defMeanDcaX, defMeanDcaY, defMeanDcaZ,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNsigmaPdca,
                              defChi2NormCut,
                              defRelPresolution, defSlopeResolution,
                              aptCut, lut1, lut4 );
  FillDataAndMC ( oadbObj, 239319, 241544, "muon_calo_pass2", kTRUE );

  cloneOadbObj = CloneOADBobject("LHC15ghij",oadbObj);
  FillDataAndMC ( cloneOadbObj, 228936, 238621, "muon_calo_pass2", kTRUE );


  // pPb 5.02 TeV 2013
  oadbObj = CreateOADBObject ( "LHC13def",
                              0.18943, -0.250591, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNsigmaPdca,
                              defChi2NormCut,
                              defRelPresolution, defSlopeResolution,
                              aptCut, lut05, lut4 );
  FillDataAndMC ( oadbObj, 195682, 197388, "muon_pass2", kTRUE );

  cloneOadbObj = CloneOADBobject("LHC13f",oadbObj);
  FillDataAndMC ( cloneOadbObj, 196433, 197388, "muon_calo", kFALSE );


  // pp 8 TeV 2012
  oadbObj = CreateOADBObject ( "LHC12hi",
                              0.18943, -0.250591, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              defSigmaPdca23, defSigmaPdca310, defNsigmaPdca,
                              defChi2NormCut,
                              defRelPresolution, defSlopeResolution,
                              aptCut, lut1, lut4 );
  FillDataAndMC ( oadbObj, 189576, 193341, "muon_calo_pass2", kTRUE );


  // PbPb 2.76 TeV 2011 + pp 7 TeV 2011
  oadbObj = CreateOADBObject ( "LHC11h_2",
                              -1.118002, -1.175119, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              99., 54., 6.,
                              defChi2NormCut,
                              5.e-4, 6.e-4,
                              aptCut, lut1, lut4 );
  FillDataAndMC ( oadbObj, 167706, 170593, "pass2", kFALSE );

  oadbObj = CreateOADBObject ( "LHC11c3defgh",
                              0.296514, -0.229262, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              99., 54., 6.,
                              defChi2NormCut,
                              defRelPresolution, defSlopeResolution,
                              aptCut, lut1, lut4 );
  FillDataAndMC ( oadbObj, 154726, 170593, "pass2_muon", kTRUE );

  cloneOadbObj = CloneOADBobject ("LHC11c2",oadbObj);
  cloneOadbObj->SetSharpPtCut(aptCut, lut1, lut2);
  FillDataAndMC ( cloneOadbObj, 153059, 154495, "pass2_muon", kTRUE );

  cloneOadbObj = CloneOADBobject ("LHC11c1",oadbObj);
  cloneOadbObj->SetSharpPtCut(aptCut, lut05, lut1);
  FillDataAndMC ( cloneOadbObj, 151661, 152935, "pass2_muon", kTRUE );

  oadbObj = CreateOADBObject ( "LHC11h",
                              -1.118002, -1.175119, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              99., 54., 6.,
                              defChi2NormCut,
                              5.e-4, 6.e-4,
                              aptCut, lut1, lut4 );
  FillDataAndMC ( oadbObj, 167706, 170593, "pass1_muon", kFALSE );

  oadbObj = CreateOADBObject ( "LHC11def",
                              -1.146336, -1.130535, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              99., 54., 6.,
                              defChi2NormCut,
                              5.e-4, 6.e-4,
                              aptCut, lut1, lut4 );
  FillDataAndMC ( oadbObj, 156620, 162717, "pass1", kFALSE );


  // pp 2.76 TeV 2011
  oadbObj = CreateOADBObject ( "LHC11a",
                              -1.146336, -1.130535, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              99., 54., 6.,
                              defChi2NormCut,
                              5.e-4, 6.e-4,
                              aptCut, lut05, lut1 );
  FillDataAndMC ( oadbObj, 146688, 146860, "pass1", kTRUE );


  // PbPb 2.76 TeV 2010
  oadbObj = CreateOADBObject ( "LHC10h",
                              -0.4599, -0.9172, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              99., 54., 6.,
                              3.5,
                              5.e-4, 6.e-4,
                              aptCut, lut05, lut1 );
  FillDataAndMC ( oadbObj, 137135, 139513, "pass1", kTRUE );


  // pp 7 TeV 2010
  oadbObj = CreateOADBObject ( "LHC10pp",
                              -0.4599, -0.9172, 0.,
                              defMeanPCorr23, defMeanPCorr310,
                              99., 54., 6.,
                              defChi2NormCut,
                              5.e-4, 6.e-4,
                              aptCut, lut05, lut1 );
  FillDataAndMC ( oadbObj, 114783, 136376, "pass1", kTRUE );



  TString oadbFilename = "$ALICE_PHYSICS/OADB/PWG/MUON/MuonTrackCuts.root";
  gSystem->ExpandPathName(oadbFilename);
  if ( ! gSystem->AccessPathName(oadbFilename.Data()) ) {
    gSystem->Exec(Form("rm %s", oadbFilename.Data()));
  }

  for ( Int_t icont=0; icont<kNcontainers; icont++ ) {
    containers[icont]
    ->WriteToFile ( oadbFilename.Data() );
  }

  printf("Produced OADB file: %s\n",oadbFilename.Data());
}

//_____________________________________________
void PrintOADB ( TString oadbFilename = "$ALICE_PHYSICS/OADB/PWG/MUON/MuonTrackCuts.root" )
{
  gSystem->ExpandPathName(oadbFilename);
  if ( gSystem->AccessPathName(oadbFilename.Data()) ) {
    printf("Error: cannot open %s\n",oadbFilename.Data());
    return;
  }
  TFile* file = TFile::Open(oadbFilename.Data());
  for ( Int_t icont=0; icont<file->GetListOfKeys()->GetEntries(); icont++ ) {
    TString currKey = file->GetListOfKeys()->At(icont)->GetName();
    AliOADBContainer* container = dynamic_cast<AliOADBContainer*>(file->Get(currKey.Data()));
    if ( ! container ) continue;
    printf("\n\nContainer: %s\n",currKey.Data());
    for ( Int_t idx=0; idx<container->GetNumberOfEntries(); idx++ ) {
      TObject* obj = container->GetObjectByIndex(idx);
      printf("\n%s  (%i - %i)\n",obj->GetName(),container->LowerLimit(idx),container->UpperLimit(idx));
      obj->Print();
    }
    printf("\nDefault:\n");
    container->GetDefaultList()->Print();
  }
}
