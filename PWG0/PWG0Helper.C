/* $Id$ */

// Helper macros can be found in this file
// A set of them can be used to connect to proof and execute selectors.

TVirtualProof* connectProof(const char* proofServer)
{
  TVirtualProof* proof = TProof::Open(proofServer);

  if (!proof)
  {
    printf("ERROR: PROOF connection not established.\n");
    return 0;
  }

  // enable the new packetizer
  //proof->AddInput(new TNamed("PROOF_Packetizer", "TPacketizerProgressive"));

  proof->ClearInput();

  return proof;
}

Bool_t prepareQuery(TString libraries, TString packages, Int_t useAliRoot)
{
  // if not proof load libraries
  if (!gProof)
  {
    TObjArray* librariesList = libraries.Tokenize(";");
    for (Int_t i=0; i<librariesList->GetEntries(); ++i)
    {
      TObjString* str = dynamic_cast<TObjString*> (librariesList->At(i));
      if (!str)
        continue;

      printf("Loading %s...", str->String().Data());
      Int_t result = CheckLoadLibrary(str->String());
      if (result < 0)
      {
        printf("failed\n");
        //return kFALSE;
      }
      else
        printf("succeeded\n");
    }
  }
  else
  {
    if (useAliRoot > 0)
      ProofEnableAliRoot(useAliRoot);

    TObjArray* packagesList = packages.Tokenize(";");
    for (Int_t i=0; i<packagesList->GetEntries(); ++i)
    {
      TObjString* str = dynamic_cast<TObjString*> (packagesList->At(i));
      if (!str)
        continue;

      /*if (!EnablePackageLocal(str->String()))
      {
        printf("Loading of package %s locally failed\n", str->String().Data());
        return kFALSE;
      }*/

      if (gProof->UploadPackage(Form("%s.par", str->String().Data())))
      {
        printf("Uploading of package %s failed\n", str->String().Data());
        return kFALSE;
      }

      if (gProof->EnablePackage(str->String()))
      {
        printf("Loading of package %s failed\n", str->String().Data());
        return kFALSE;
      }
    }
 }

 return kTRUE;
}

Int_t executeQuery(TChain* chain, TList* inputList, TString selectorName, const char* option = "")
{
  if (!gProof)
    chain->GetUserInfo()->AddAll(inputList);
  else
  {
    for (Int_t i=0; i<inputList->GetEntries(); ++i)
      gProof->AddInput(inputList->At(i));
  }

  TStopwatch timer;
  timer.Start();

  Long64_t result = -1;

  if (gProof)
    result = chain->MakeTDSet()->Process(selectorName, option);
  else
    result = chain->Process(selectorName, option);

  if (result < 0)
    printf("ERROR: Executing process failed with %d.\n", result);

  timer.Stop();
  timer.Print();

  return result;
}

void ProofEnableAliRoot(Int_t aliroot)
{
  // enables a locally deployed AliRoot in a PROOF cluster

  /* executes the following commands on each node:
     gSystem->Setenv("ALICE_ROOT", "/home/alicecaf/ALICE/aliroot-head")
     gSystem->AddIncludePath("/home/alicecaf/ALICE/aliroot-head/include");
     gSystem->SetDynamicPath(Form("%s:%s", gSystem->GetDynamicPath(), "/home/alicecaf/ALICE/aliroot-head/lib/tgt_linux"))
     gSystem->Load("libMinuit");
     gROOT->Macro("$ALICE_ROOT/macros/loadlibs.C");
  */

  const char* location = 0;
	const char* target = "tgt_linux";
  
  switch (aliroot)
  {
    case 1: location = "/home/alicecaf/ALICE/aliroot-v4-04-Release"; break;
    case 2: location = "/home/alicecaf/ALICE/aliroot-head"; break;
		case 11: location = "/data1/qfiete/aliroot-head"; target = "tgt_linuxx8664gcc"; break;
    default: return;
  }

  gProof->Exec(Form("gSystem->Setenv(\"ALICE_ROOT\", \"%s\")", location), kTRUE);
  gProof->AddIncludePath(Form("%s/include", location));
  gProof->AddDynamicPath(Form("%s/lib/%s", location, target));

  // load all libraries
  gProof->Exec("gSystem->Load(\"libMinuit\")");
  gProof->Exec("gROOT->Macro(\"$ALICE_ROOT/macros/loadlibs.C\")");
}

Bool_t EnablePackageLocal(const char* package)
{
  printf("Enabling package %s locally...\n", package);

  TString currentDir(gSystem->pwd());
  if (!gSystem->cd(package))
    return kFALSE;

  gROOT->ProcessLine(".x PROOF-INF/SETUP.C");
  gSystem->cd(currentDir);

  return kTRUE;
}

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;

  return gSystem->Load(library);
}

void redeployPackages(const char* proofServer, Bool_t localAliRoot = kTRUE)
{
  // deploys PWG0base and PWG0dep (the latter only when localAliRoot is true) that are expected in $ALICE_ROOT
  // when localAliRoot is false ESD.par is also deployed

  TProof::Reset(proofServer);
  TVirtualProof* proof = TProof::Open(proofServer);
  proof->ClearPackages();

  if (localAliRoot)
    ProofEnableAliRoot();
  else
  {
    proof->UploadPackage("$ALICE_ROOT/ESD.par");
    proof->EnablePackage("ESD");
  }

  proof->UploadPackage("$ALICE_ROOT/PWG0base.par");
  proof->EnablePackage("PWG0base");

  proof->UploadPackage("$ALICE_ROOT/PWG0dep.par");
  proof->EnablePackage("PWG0dep");
}
