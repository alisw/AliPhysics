#ifndef FASTSIM_H
#define FASTSIM_H
#include <TSelector.h>
#include <TQObject.h>
#ifndef __CINT__
# include "AliGenerator.h"
# include "AliRunLoader.h"
# include "AliStack.h"
# include "AliHeader.h"
# include "AliGenEventHeader.h"
# include "AliRun.h"
# include "AliCollisionGeometry.h"
# include "AliGenPythiaEventHeader.h"
# include "AliGenDPMjetEventHeader.h"
# include "AliGenGeVSimEventHeader.h"
# include "AliGenHerwigEventHeader.h"
# include "AliGenHijingEventHeader.h"
# include "AliGenCocktailEventHeader.h"
# include "AliPDG.h"
# include <TROOT.h>
# include <TString.h>
# include <TMath.h>
# include <TParticle.h>
# include <TH1.h>
# include <TTree.h>
# include <TClonesArray.h>
# include <TList.h>
# include <TProof.h>
# include <TChain.h>
# include <TParticlePDG.h>
# include <TStopwatch.h>
# include <TFile.h>
# include <TProofOutputFile.h>
# include <TCanvas.h>
# include <TTimer.h>
# include <TRandom.h>
# include <TUrl.h>
# include <TMacro.h>
# include <TSystemDirectory.h>
# include <TFileCollection.h>
# include <TPRegexp.h>
# include <fstream>
# include <TLeaf.h>
# include <algorithm>
# include "FastShortHeader.C"
# include "FastCentEstimators.C"
# include "FastMonitor.C"
#else
class AliGenerator;
class AliRunLoader;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class TH1;
class TTree;
class TClonesArray;
class TBrowser;
class TList;
class TFile;
class TProofOutputFile;
class TCanvas;
class TVirtualPad;
class TTimer;
class TUrl;
class TAxis;
class TParticle;
class TMacro;
class TLeaf;
class BCentEstimator;
class FastMonitor;
#endif


  
/** To get DPMJEt common block */
typedef struct {
   Double_t	rproj;
   Double_t	rtarg;
   Double_t	bimpac;
   Int_t	nwtsam;
   Int_t	nwasam;
   Int_t	nwbsam;
   Int_t	nwtacc;
   Int_t	nwtaac;
   Int_t	nwtbac;
   Int_t        ncp;
   Int_t        nct;
} DtglcpCommon;
DtglcpCommon* _dtglcp = 0;

//====================================================================
struct FastSimMonitor : public FastMonitor
{
  FastSimMonitor(TSelector* s=0)
    : FastMonitor(s, "FastMonitor")
  {
    Register("list/histograms/b",       "hist", 0,    false);
    Register("list/histograms/cent",    "hist", 0,    false);
    Register("list/histograms/timing",  "hist", 0x18, false);
    Register("list/histograms/dNdeta",  "",     0x18, false);
    Register("list/histograms/dNdy",    "",     0x18, false);
    Register("list/histograms/trigger", "e",    0x10, false);
    Register("list/estimators/rawV0MP", "e",    0x02, false);
    if (gProof)  gProof->AddFeedback("list");      
  }
};

//====================================================================
/** 
 * Run a event generator simulation 
 */
struct FastSim : public TSelector
{
  /** 
   * Constructor 
   * 
   * @param eg      Event generator 
   * @param runNo   Run number to simulate 
   * @param bMin    Lease impact parameter 
   * @param bMax    Largest impact parameter 
   * @param nEvents Number of events 
   * @param monitor Monitor frequency 
   */
  FastSim(const char* eg="",
	  ULong_t     runNo=0,
	  Double_t    bMin=0,
	  Double_t    bMax=20,
	  Long64_t    nEvents=0,
	  Int_t       monitor=0)
    : TSelector(),
      fEGName(eg),
      fRunNo(runNo),
      fBMin(bMin),
      fBMax(bMax),
      fGRP(0),
      fOverrides(0),
      fNEvents(nEvents),
      fIsTgtA(false),
      fIsProjA(false),
      fGenerator(0),
      fRunLoader(0),
      fStack(0),
      fHeader(0),
      fTree(0),
      fParticles(0),
      fList(0),
      fHEta(0),
      fHY(0),
      fHIpz(0),
      fHType(0),
      fHCent(0),
      fHB(0),
      fHPhiR(0),
      fHTime(0),
      fHTrig(0),
      fBEstimator(0),
      fCentEstimators(0),
      fProofFile(0),
      fAliceFile(0),
      fKineFile(0),
      fFile(0),
      fFileName(""),
      fVerbose(true),
      fMonitor(monitor)
  {}
  const char* FileName() const
  {
    static TString fn;
    if (fn.IsNull()) {
      if (!fFileName.IsNull())  fn = fFileName;
      else {
	TString egName = (fGenerator ? fGenerator->GetName() : "");
	if (egName.IsNull()) egName = fEGName;
	
	fn = Form("%s_%09d", egName.Data(), fRunNo);
	if (fGenerator) {
	  TString tgt, proj;
	  Int_t   tgtA, tgtZ, projA, projZ;
	  fGenerator->GetTarget(tgt, tgtA, tgtZ);
	  fGenerator->GetProjectile(proj, projA, projZ);
	  fn.Append(Form("_%s%s", tgt.Data(), proj.Data()));
	  fn.Append(Form("_%05d", Int_t(fGenerator->GetEnergyCMS())));
	}
	else {
	  Long_t en = gROOT->ProcessLine("grp->energy");
	  // Bool_t a1 = gROOT->ProcessLine("(GRPData*)%p->beam1.IsA()",fGRP);
	  // Bool_t a2 = gROOT->ProcessLine("(GRPData*)%p->beam2.IsA()",fGRP);
	  fn.Append(Form("_%s%s", (fIsTgtA ? "A" : "p"),
			 (fIsProjA ? "A" : "p")));
	  fn.Append(Form("_%05ld", (fGRP ? en : 0)));
	}

	if (fNEvents > 0) {
	  if (fNEvents >= 1000000)
	    fn.Append(Form("_%lldM", fNEvents/1000000));
	  else if (fNEvents >= 1000)
	    fn.Append(Form("_%lldk", fNEvents/1000));
	  else
	    fn.Append(Form("_%lld", fNEvents));
	}
	fn.Append(".root");
	fFileName = fn;
      }
    }
    return fn.Data();
  }
  /** 
   * Get the name of the selector 
   * 
   * @return Name of selector 
   */
  const char* GetName() const { return "FastSim"; }
  /** 
   * Get the title of the selector 
   * 
   * @return The title 
   */
  const char* GetTitle() const { return "ALICE Event Generator simulation"; }
  /** 
   * Get collision energy either form generator or GRP
   * 
   * @return Collision energy in GeV
   */
  virtual UShort_t GetSNN() const
  {
    Float_t fsNN = 0;
    if (fGenerator) fsNN = fGenerator->GetEnergyCMS();
    else            fsNN = gROOT->ProcessLine("grp->energy");
      
    UShort_t   sNN      = (TMath::Abs(fsNN -  2760) < 10 ?  2760 :
			   TMath::Abs(fsNN -  5023) < 10 ?  5023 :
			   TMath::Abs(fsNN -  5440) < 10 ?  5440 :
			   TMath::Abs(fsNN -  2360) < 10 ?  2360 :
			   TMath::Abs(fsNN -   900) < 10 ?   900 :
			   TMath::Abs(fsNN -  7000) < 10 ?  7000 :
			   TMath::Abs(fsNN -  8000) < 10 ?  8000 :
			   TMath::Abs(fsNN - 13000) < 10 ? 13000 :
			   0);
    return sNN;
  }
  /** 
   * Get the event generator title 
   * 
   * @return Event generator title 
   */
  virtual const char* GetEGTitle() const
  {
    
    if (fGenerator) return fGenerator->GetTitle();
    return "";
  }
  /** 
   * Create our outputs 
   * 
   * @return true on success 
   */
  Bool_t SetupOutput()
  {
    Printf("=== Output =============================\n"
	   " File name: %s\n"
	   "========================================", FileName());

    if (fVerbose) Info("SetupOutput", "First the file");
    Bool_t isProof = false;
    if (fInput && fInput->FindObject("PROOF_Ordinal"))
      isProof = true;
    if (isProof) {
      Info("SetupOutput", "Making Proof File");
      fProofFile = new TProofOutputFile(FileName(), "M");
      // TProofOutputFile::kMerge,
      // TProofOutputFile::kRemote);
      fFile = fProofFile->OpenFile("RECREATE");
    }
    else
      fFile = TFile::Open(FileName(), "RECREATE");

    UShort_t sNN = GetSNN();
    TString  tit = GetEGTitle();
    if (fVerbose) Info("SetupOutput", "Making our tree: %s", tit.Data());
    fTree      = new TTree("T", tit.Data());
    fParticles = new TClonesArray("TParticle");
    fTree->Branch("header", &fShortHead,
		  "run/i:event:ntgt:nproj:nbin:type:"
		  "ipx/D:ipy:ipz:b:c:phir:"
		  "nsnp/i:nsnt:nspp:nspt:mask");
    fTree->Branch("particles", &fParticles);
    fTree->SetAutoSave(500); // Save every on every 100 events 
    fTree->SetDirectory(fFile);
    fTree->SetAlias("primary", "(particles.fBits&(1<<14))");
    fTree->SetAlias("weak",    "(particles.fBits&(1<<15))");
    fTree->SetAlias("charged", "(particles.fBits&(1<<16))");
    fTree->SetAlias("pt",      "(sqrt(pow(particles.fPx,2)+"
		    /*       */"pow(particles.fPy,2)))");
    fTree->SetAlias("eta",     "(particles.Pt()<1e-100?"
		    "(particles.fPz>0?100:-100):"
		    "-log(tan(atan2(particles.Pt(),particles.fPz)/2)))");
    fTree->SetAlias("good",    "(primary&&charged&&abs(eta)<1000)");
    fTree->SetAlias("sd",      "(header.fType & 0x1)");
    fTree->SetAlias("dd",      "(header.fType & 0x2)");
    fTree->SetAlias("pion",    "(abs(particles.fPdgCode)==211)");
    fTree->SetAlias("kaon",    "(abs(particles.fPdgCode)==321)");
    fTree->SetAlias("proton",  "(abs(particles.fPdgCode)==2212)");
    fTree->SetAlias("electron","(abs(particles.fPdgCode)==11)");
    fTree->SetAlias("other",   "(!pion&&!kaon&&!proton&&!electron)");
    fTree->SetAlias("beta",    "(particles.P()/particle.Energy())");
    fTree->SetAlias("gamma",   "(1./sqrt(1-beta*beta))");
    fTree->SetAlias("npart",   "(header.ntgt+header.nproj)");
    fTree->SetAlias("v0a",     "(header.mask&0x1)");
    fTree->SetAlias("v0c",     "(header.mask&0x2)");
    fTree->SetAlias("ada",     "(header.mask&0x4)");
    fTree->SetAlias("adc",     "(header.mask&0x8)");
    fTree->SetAlias("eta1",    "(header.mask&0x10)");
    
    // Let's add the title of the generator to the output.  We make a
    // histogram because that can be merged.
    Info("SetupOutput", "Making generator histogram: %s",tit.Data());
    TH1* egTitle = new TH1I("eg", tit.Data(), 1, 0, 1);
    egTitle->SetDirectory(0);
    egTitle->Fill(0.5);
    egTitle->SetXTitle("N_{\\hbox{workers}}");
    egTitle->SetYTitle("Count");
    egTitle->SetFillColor(kYellow+4);
    egTitle->SetLineColor(kYellow+4);
    egTitle->SetMarkerColor(kYellow+4);
    egTitle->SetFillStyle(1001);
    egTitle->SetMarkerStyle(1);
    egTitle->SetLineStyle(1);
    egTitle->SetLineWidth(2);
    egTitle->SetMarkerSize(1);
    
    
    if (fVerbose) Info("SetupOutput", "Making histograms");
    Double_t maxEta = 10;
    Double_t dEta   = 10./200;
    fHEta = new TH1D("dNdeta", "Charged particle pseudo-rapidity density",
		     Int_t(2*maxEta/dEta+.5), -maxEta, +maxEta);
    fHEta->Sumw2();
    fHEta->SetXTitle("\\eta");
    fHEta->SetYTitle("\\hbox{d}N_{\\hbox{ch}}/\\hbox{d}\\eta");
    fHEta->SetMarkerColor(kRed+2);
    fHEta->SetMarkerStyle(20);
    fHEta->SetDirectory(0);

    fHY = new TH1D("dNdy", "Charged particle rapidity density",
		     Int_t(2*maxEta/dEta+.5), -maxEta, +maxEta);
    fHY->Sumw2();
    fHY->SetXTitle("\\mathit{y}");
    fHY->SetYTitle("\\hbox{d}N_{\\hbox{ch}}/\\hbox{d}y");
    fHY->SetMarkerColor(kRed+2);
    fHY->SetMarkerStyle(20);
    fHY->SetDirectory(0);
    
    fHIpz = new TH1D("ipZ", "Z-coordinate of interaction point",
		     10, -10, 10);
    fHIpz->SetMarkerColor(kGreen+2);
    fHIpz->SetFillColor(kGreen+2);
    fHIpz->SetFillStyle(3001);
    fHIpz->SetXTitle("IP_{#it{z}} [cm]");
    fHIpz->SetDirectory(0);

    fHType = new TH1D("type", "Diffractive", 3, .5, 3.5);
    fHType->SetMarkerColor(kOrange+2);
    fHType->SetFillColor(kOrange+2);
    fHType->SetFillStyle(3001);
    fHType->SetDirectory(0);
    fHType->GetXaxis()->SetBinLabel(1, "Non");
    fHType->GetXaxis()->SetBinLabel(2, "Single");
    fHType->GetXaxis()->SetBinLabel(3, "Double");

    fHCent = new TH1D("cent", "Centrality", 20, 0, 100);
    fHCent->SetMarkerColor(kPink+2);
    fHCent->SetFillColor(kPink+2);
    fHCent->SetFillStyle(3001);
    fHCent->SetDirectory(0);
    fHCent->SetYTitle("Events");
    fHCent->SetXTitle("Centrality [%]");
    
    fHB = new TH1D("b", "Impact parameter", 20, 0, 20);
    fHB->SetMarkerColor(kYellow+2);
    fHB->SetFillColor(kYellow+2);
    fHB->SetFillStyle(3001);
    fHB->SetYTitle("Events");
    fHB->SetXTitle("#it{b} [fm]");
    fHB->SetDirectory(0);

    fHPhiR = new TH1D("phiR", "Event plane angle", 360, 0, 360);
    fHPhiR->SetMarkerColor(kMagenta+2);
    fHPhiR->SetFillColor(kMagenta+2);
    fHPhiR->SetFillStyle(3001);
    fHPhiR->SetXTitle("#it{#Phi}_{R} [degrees]");
    fHPhiR->SetDirectory(0);

    fHTime = new TH1D("timing", "Timing of processing", 5,0.5,5.5);
    fHTime->SetMarkerColor(kBlue+2);
    fHTime->SetFillColor(kBlue+2);
    fHTime->SetFillStyle(3001);
    fHTime->SetYTitle("seconds / event");
    fHTime->GetXaxis()->SetBinLabel(1,"Reset");
    fHTime->GetXaxis()->SetBinLabel(2,"Generate");
    fHTime->GetXaxis()->SetBinLabel(3,"Header");
    fHTime->GetXaxis()->SetBinLabel(4,"Particles");
    fHTime->GetXaxis()->SetBinLabel(5,"Filling");
    fHTime->SetDirectory(0);

    fHTrig = new TH1D("trigger", "Trigger bits set", 6, -0.5, 5.5);
    fHTrig->SetMarkerColor(kMagenta+2);
    fHTrig->SetFillColor(kMagenta+2);
    fHTrig->SetFillStyle(3001);
    fHTrig->SetYTitle("Events");
    fHTrig->GetXaxis()->SetBinLabel(1, "All");
    fHTrig->GetXaxis()->SetBinLabel(2, "V0A");
    fHTrig->GetXaxis()->SetBinLabel(3, "V0C");
    fHTrig->GetXaxis()->SetBinLabel(4, "ADA");
    fHTrig->GetXaxis()->SetBinLabel(5, "ADC");
    fHTrig->GetXaxis()->SetBinLabel(6, "N_{ch}|_{|#eta|<1}#ge1");
    fHTrig->SetDirectory(0);
    
    fList = new TList;
    fList->SetName("list");
    fList->Add(egTitle);

    TList* histos = new TList;
    histos->SetName("histograms");
    histos->SetOwner(true);
    histos->Add(fHEta);
    histos->Add(fHY);
    histos->Add(fHIpz);
    histos->Add(fHType);
    histos->Add(fHCent);
    histos->Add(fHB);
    histos->Add(fHPhiR);
    histos->Add(fHTime);
    histos->Add(fHTrig);
    fList->Add(histos);
    
    TList* estimators = new TList;
    estimators->SetName("estimators");
    estimators->SetOwner(true);
    fList->Add(estimators);

    
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next()))) {
      Info("SetupOutput", "Setting up estimator %s", estimator->GetName());
      estimator->Setup(estimators, fTree,sNN,fIsTgtA,fIsProjA);
      estimator->SetVerbose(fVerbose);
      // estimator->Print("nah");
    }
    
    if (fVerbose) Info("SetupOutput", "Adding list ot outputs");
    fOutput->Add(fList);
    // fOutput->ls();
    
    return true;
  }
  /** 
   * Set our seed.  This will read a UInt_t worth of data from
   * /dev/urandom
   * 
   */
  virtual void SetupSeed()
  {
    std::ifstream in("/dev/urandom");
    UInt_t seed = 0;
    in.read(reinterpret_cast<char*>(&seed), sizeof(seed));
    in.close();
    Printf("=== Random =============================\n"
	   " Random number seed: %d\n"
	   "========================================", seed);
    gRandom->SetSeed(seed);
  }
  virtual Bool_t SetupGRP()
  {
    Printf(" === Setup ==============================");
    Printf("  Run #:                          %6d", fRunNo);
    Printf("  EG:     %30s", fEGName.Data());
    Printf("  B range:             %5.1ffm - %5.1ffm", fBMin, fBMax);
    Printf(" ========================================");
    Printf("Macro path: %s", gROOT->GetMacroPath());

    // --- Check if we shoud get the GRP line ------------------------
    if (!fGRP && fInput) {
      fGRP = fInput->FindObject("GRP");
      std::ofstream* pout = new std::ofstream("grp.dat");
      if (pout) {
	if (fVerbose) 
	  Info("SetupGRP", "Writing GRP line '%s' to \"grp.dat\"",
	       fGRP->GetTitle());
	std::ostream& out = *pout;
	out << fGRP->GetTitle() << std::endl;
	pout->close();
      }
    }
    if (fVerbose) 
      Info("SetupGRP", "Overrides: %p Input: %p", fOverrides, fInput);
    if (!fOverrides && fInput) {
      fOverrides = static_cast<TList*>(fInput->FindObject("overrides"));
      if (!fOverrides && fVerbose)
	Info("SetupGRP", "No GRP overrides found in input:");
    }
    
    // --- Load our settings -----------------------------------------
    if (fVerbose) Info("SetupGRP", "Loading scripts");
    // Check if we have the global "grp" already 
    if (gROOT->ProcessLine("grp") == 0) 
      gROOT->Macro(Form("GRP.C(%d)", fRunNo));
    if (fVerbose) Info("SetupGRP", "Perhaps override");
    OverrideGRP();
    gROOT->ProcessLine("grp->Print()");
    return true;
  }    
  /** 
   * Set-up the generator. 
   * 
   * @return true on success 
   */
  virtual Bool_t SetupGen()
  {
    if (fVerbose) Info("SetupGen", "Load base config");
    gROOT->Macro("BaseConfig.C");
    if (fVerbose) Info("SetupGen", "Load EG config");
    gROOT->Macro("EGConfig.C");

    gROOT->ProcessLine(Form("VirtualEGCfg::LoadGen(\"%s\")",fEGName.Data()));

    // --- Make our generator ----------------------------------------
    // Info("SetupGen", "Creating generator");
    TString egMk = Form("egCfg->MakeGenerator(\"%s\",%f,%f)",
			fEGName.Data(), fBMin, fBMax);
    Long64_t egPtr = gROOT->ProcessLine(egMk);
    if (egPtr == 0) {
      Error("SetupGen", "Failed to make generator");
      return false;
    }
    fGenerator = reinterpret_cast<AliGenerator*>(egPtr);
    TString tgt, proj;
    Int_t   tgtA=0, tgtZ=0, projA=0, projZ=0;
    fGenerator->GetTarget(tgt, tgtA, tgtZ);
    fGenerator->GetProjectile(proj, projA, projZ);
    fIsTgtA  = !(tgtA  == tgtZ  && tgtA == 1);
    fIsProjA = !(projA == projZ && projZ == 1);
    if (fVerbose) 
      Info("SetupGen", "tgt=%s (%3d,%2d) proj=%s (%3d,%2d) CMS=%fGeV",
	   tgt.Data(), tgtA, tgtZ, proj.Data(), projA, projZ,
	   fGenerator->GetEnergyCMS());
    Info("SetupGen", "Generator: %s", fGenerator->GetTitle());
    if (fFileName.IsNull()) FileName();

    return true;
  }
  /** 
   * Setup the generator etc. of the job 
   * 
   * @return true on success 
   */
  virtual Bool_t SetupRun()
  {
    // --- gAlice (bare ROOT) ----------------------------------------
    if (!gAlice)
      new AliRun("gAlice", "The ALICE Off-line framework");

    Long64_t nev = (fNEvents <= 0 ? 0xFFFFFFFF : fNEvents);
    Printf("=== Run ================================\n"
	   " Number of events: %lld\n"
	   "========================================", nev);
    TObject* ord      = (fInput ? fInput->FindObject("PROOF_Ordinal") : 0);
    UShort_t saveMode = 0;
    TString  post     = "";
    TString  dir      = "";
    if (ord) {
      TObject* save = fInput->FindObject("PROOF_SaveGALICE");
      if (save && fVerbose) {
	Info("SetupRun", "Got save option:");
	save->Print();
      }
      TString optSave(save ? save->GetTitle() : "split");
      optSave.ToLower();
      if       (optSave.EqualTo("none"))   saveMode = 0;
      else if  (optSave.EqualTo("merge"))  saveMode = 1;
      else if  (optSave.EqualTo("split"))  saveMode = 2;
      if (fProofFile && saveMode > 0) 
	dir  = fProofFile->GetDir(true);
      if (saveMode > 1)
	post = Form("_%s", ord->GetTitle());	
    }
    TString  galiceName(Form("%sgalice.root",dir.Data()));
    TString  kineName(Form("%sKinematics.root",dir.Data()));
    
    // --- Run-loader, stack, etc  -----------------------------------
    // Info("SetupRun", "Set-up run Loader");    
    fRunLoader = AliRunLoader::Open(galiceName, "FASTRUN", "RECREATE");
    fRunLoader->SetKineFileName(kineName);
    fRunLoader->SetCompressionLevel(2);
    fRunLoader->SetNumberOfEventsPerFile(nev);
    fRunLoader->LoadKinematics("RECREATE");
    fRunLoader->MakeTree("E");
    gAlice->SetRunLoader(fRunLoader);
    fRunLoader->MakeStack();
    fStack  = fRunLoader->Stack();
    fHeader = fRunLoader->GetHeader();

    // --- Initialize generator --------------------------------------
    // Info("SetupRun", "Initializing generator");
    fGenerator->Init();
    fGenerator->SetStack(fStack);

    if (saveMode < 1) {
      if (ord) 
	Info("SetupRun", "Not saving galice.root and Kinematics.root");
      return true;
    }
    
    TString aliceOut = Form("galice%s.root", post.Data());
    fAliceFile = new TProofOutputFile(aliceOut, "M");

    TString kineOut = Form("Kinematics%s.root", post.Data());
    fKineFile = new TProofOutputFile(kineOut, "M");
    
    return true;
  }
  /** 
   * Read the previously created grp.dat file 
   * 
   */
  Bool_t ReadGRPLine()
  {
    std::ifstream* pin = new std::ifstream("grp.dat");
    if (!pin) {
      Warning("ReadGRPLine", "Failed to open \"grp.dat\"");
      return false;
    }
    std::istream&  in  = *pin;
    TString line;
    TString env;
    do {
      line.ReadLine(in);
      if (line.IsNull()) continue;
      if (line.BeginsWith("#")) continue;
      env = line;
      break;
    } while (!in.eof());
    pin->close();

    if (env.IsNull()) {
      Warning("ReadGRPLine", "Got no line from \"grp.dat\"");
      return false;
    }
    
    fGRP = new TNamed("GRP",env.Data());
    if (fVerbose) Info("ReadGRPLine", "Read \"%s\"", env.Data());
    return true;
  }
  /** 
   * Possibly override settings from GRP. 
   * 
   */
  void OverrideGRP()
  {
    Long_t ret = gROOT->ProcessLine("grp");
    if (ret == 0) {
      Warning("OverrideGRP", "GRP not set yet, cannot override");
      return;
    }
    if (!fOverrides) {
      if (fVerbose) Info("OverrideGRP", "No overrides defined");
      return;
    }
    TIter next(fOverrides);
    TObject* o = 0;
    while ((o = next())) {
      if (fVerbose)
	Info("OverrideGRP", "Overriding GRP setting %s with %s",
	     o->GetName(), o->GetTitle());
      gROOT->ProcessLine(Form("grp->%s = %s;",
			      o->GetName(), o->GetTitle()));
    }
    Info("OverrideGRP", "After overriding:");
    gROOT->ProcessLine("grp->Print()");
  }
  /** 
   * Add an item to the list of things from GRP to override
   * 
   * @param field Field name of GRPData
   * @param value Field value of GRPData
   */
  void AddOverride(const TString& field, const TString& value)
  {
    if (!fOverrides) {
      fOverrides = new TList;
      fOverrides->SetName("overrides");
    }
    fOverrides->Add(new TNamed(field, value));
  }
  /** 
   * Set up job 
   * 
   */
  virtual void Init(TTree*)
  {
  }
  /** 
   * Set up job 
   * 
   */
  virtual void Begin(TTree*)
  {
    // Make a monitor
    // Info("Begin", "gProof=%p Nomonitor=%p",
    //      gProof, (gProof ? gProof->GetParameter("NOMONITOR") : 0));
    if (fVerbose) Info("Begin", "Called for FastSim");

    if (fMonitor > 0 && !gROOT->IsBatch()) {
      FastSimMonitor* m = new FastSimMonitor;
      m->Connect(fMonitor);
    }
    gROOT->Macro(Form("GRP.C(%d)", fRunNo));
    if (ReadGRPLine()) {
      if(gProof) {
	gProof->AddInput(fGRP);
	if (fOverrides) gProof->AddInput(fOverrides);
      }
    }
    if (fVerbose) Info("Begin", "Perhaps override");
    OverrideGRP();
    if (fVerbose) Info("Begin", "Defining centrality estimators");

    fBEstimator = new BCentEstimator;
    fCentEstimators = new TList;
    // fCentEstximators->Add(new V0CentEstimator(-1));           //V0C
    // fCentEstimators->Add(new V0CentEstimator( 0));           //V0M
    // fCentEstimators->Add(new V0CentEstimator(+1));           //V0A
    fCentEstimators->Add(fBEstimator);
    fCentEstimators->Add(new V0CentEstimator(-1,true));      //V0CP
    fCentEstimators->Add(new V0CentEstimator( 0,true));      //V0MP
    fCentEstimators->Add(new V0CentEstimator(+1,true));      //V0AP
    // fCentEstimators->Add(new ZNCentEstimator(-1,zN|zE));     //ZNCE
    // fCentEstimators->Add(new ZNCentEstimator(+1,zN|zE));     //ZNAE 
    fCentEstimators->Add(new ZNCentEstimator(-1,true,false)); //ZNCE
    fCentEstimators->Add(new ZNCentEstimator(+1,true,false)); //ZNAE 
    fCentEstimators->Add(new ZNCentEstimator(-1,true,false,true)); //ZNCEP
    fCentEstimators->Add(new ZNCentEstimator(+1,true,false,true)); //ZNAE 
    fCentEstimators->Add(new ZNCentEstimator(-1,true,true));  //ZNCS
    fCentEstimators->Add(new ZNCentEstimator(+1,true,true));  //ZNAS
    // fCentEstimators->Add(new ZNCentEstimator(-1,zE|zP));     //ZPCE
    // fCentEstimators->Add(new ZNCentEstimator(+1,zE|zP));     //ZPAE
    // fCentEstimators->Add(new ZNCentEstimator(-1,zS));        //ZPCS
    // fCentEstimators->Add(new ZNCentEstimator(+1,zS));        //ZPAS
    // fCentEstimators->Add(new RefMultEstimator(0.8));
    // fCentEstimators->Add(new RefMultEstimator(0.5));
#if 0
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next()))) {
      estimator->Print("nah");
    }
#endif 
    if (fVerbose) Info("Begin", "End of begin");
  }
  /** 
   * Set-up this sub-job 
   * 
   */
  void SlaveBegin(TTree*)
  {
     if (fVerbose) Info("SlavesBegin", "Called for FastSim");
    SetupSeed();
    SetupGRP();
    SetupGen();
    SetupOutput();
    SetupRun();
  }
  /* Reset internal caches etc. 
   * 
   * @param iEv Event number 
   * 
   * @return true on success
   */
  virtual Bool_t PreEvent(Long64_t iEv)
  {
    // --- Reset header ----------------------------------------------
    fShortHead.Reset(fRunNo, iEv);
    fParticles->Clear();
    // --- Reset header, etc.  ---------------------------------------
    fHeader->Reset(fRunNo, iEv);
    fRunLoader->SetEventNumber(iEv);
    fStack->Reset();
    fRunLoader->MakeTree("K");

    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->PreEvent();
    
    return true;
  }
    
  /** 
   * Process the event header 
   * 
   * @return true if the event should be diagnosed
   */
  virtual Bool_t ProcessHeader()
  {
    // --- Copy to short header --------------------------------------
    fShortHead.fRunNo   = fHeader->GetRun();
    fShortHead.fEventId = fHeader->GetEvent();
    TArrayF ip;
    fHeader->GenEventHeader()->PrimaryVertex(ip);
    fShortHead.fIpX     = ip[0];
    fShortHead.fIpY     = ip[1];
    fShortHead.fIpZ     = ip[2];

    // --- Check header type -----------------------------------------
    AliGenEventHeader* genHeader = fHeader->GenEventHeader();
    AliCollisionGeometry* geometry = 
      dynamic_cast<AliCollisionGeometry*>(genHeader);
    AliGenPythiaEventHeader* pythia    = 
      dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
    AliGenDPMjetEventHeader* dpm       = 
      dynamic_cast<AliGenDPMjetEventHeader*>(genHeader);
    AliGenGeVSimEventHeader* gev       = 
      dynamic_cast<AliGenGeVSimEventHeader*>(genHeader);
    AliGenHerwigEventHeader* herwig    = 
      dynamic_cast<AliGenHerwigEventHeader*>(genHeader);
    AliGenCocktailEventHeader* cocktail = 
      dynamic_cast<AliGenCocktailEventHeader*>(genHeader);
    AliGenHijingEventHeader* hijing = 
      dynamic_cast<AliGenHijingEventHeader*>(genHeader);
    if (cocktail) {
      TList* headers = cocktail->GetHeaders();
      if (!headers) Warning("", "No headers in cocktail!");
      TIter next(headers);
      AliGenEventHeader* header = 0;
      AliCollisionGeometry* geom = 0;
      while ((header = static_cast<AliGenEventHeader*>(next()))) {
	AliCollisionGeometry* g = dynamic_cast<AliCollisionGeometry*>(header);
	if (g) geom = g;
	hijing = dynamic_cast<AliGenHijingEventHeader*>(header);
      }
      if (geom) geometry = geom;
    }
    if (hijing)   fShortHead.fEG = FastShortHeader::kHijing;
    if (dpm)      fShortHead.fEG = FastShortHeader::kDPMJet;
    if (pythia)   fShortHead.fEG = FastShortHeader::kPythia;
    
    if (geometry) {
      fShortHead.fB          = geometry->ImpactParameter();
      fShortHead.fNtgt       = geometry->TargetParticipants();
      fShortHead.fNproj      = geometry->ProjectileParticipants();
      fShortHead.fNbin       = geometry->NN();
      fShortHead.fPhiR       = geometry->ReactionPlaneAngle();
      fShortHead.fNSpecNproj = geometry->ProjSpectatorsn(); 
      fShortHead.fNSpecPproj = geometry->ProjSpectatorsp();
      fShortHead.fNSpecNtgt  = geometry->TargSpectatorsn(); 
      fShortHead.fNSpecPtgt  = geometry->TargSpectatorsp();
    }
    // --- Determine diffraction flags -------------------------------
    Bool_t sd = false;
    Bool_t dd = false;
    if (pythia) {
      Int_t type = pythia->ProcessType();
      if (type < 100) { // pythia6
	switch (type) {
	case 92: case 93: sd = true; break;
	case 94:          dd = true; break;
	}
      }
      else {
	switch (type) { // Pythia8
	case 103: case 104: sd = true; break;
	case 105:           dd = true; break;
	}
      }
      fShortHead.fB     = pythia->GetImpactParameter();
      fShortHead.fNtgt  = 1;
      fShortHead.fNproj = 1;
      fShortHead.fNbin  = 1;
    }
    if (dpm) {
      // Try to figure out number of spectators in either side 
      UInt_t nSpecNproj = 0;
      UInt_t nSpecNtgt  = 0;
      UInt_t nSpecPproj = 0;
      UInt_t nSpecPtgt  = 0;
      Int_t  nPart      = fStack->GetNprimary();
      for (Int_t iPart = 0; iPart < nPart; iPart++) {
	TParticle*    particle  = fStack->Particle(iPart);
	Int_t         ks        = particle->GetStatusCode();
	Int_t         side      = 0;
	if (ks == 13) side = -1;
	if (ks == 14) side = +1;
	if (side == 0) continue;
	Int_t         kf        = particle->GetPdgCode();
	if (kf == kNeutron) {
	  if (side < 0) nSpecNproj++;
	  else          nSpecNtgt++;
	}
	else if (kf == kProton) { 
	  if (side < 0) nSpecPproj++;
	  else          nSpecPtgt++;
	}
      }
      fShortHead.fNSpecNtgt   = nSpecNtgt; 
      fShortHead.fNSpecNproj  = nSpecNproj;
      fShortHead.fNSpecPtgt   = nSpecPtgt;
      fShortHead.fNSpecPproj  = nSpecPproj;
      // fShortHead.Print();

      Int_t type = dpm->ProcessType();
#ifndef NO_DPMJET_TYPE
      switch (type) {
      case 5: case 6: sd = true;
      case 7:         dd = true;
      }      
#else
      static bool   first  = true;
      
      if (first) {
	Func_t add = gSystem->DynFindSymbol("*", "dtglcp_");
	if (!add) 
	  Warning("", "Didn't find dtglcp_");
	else 
	  _dtglcp = (DtglcpCommon*)add;
      }
      // The below - or rather a different implementation with some
      // errors - was proposed by Cvetan - I don't think it's right
      // though.  See also
      //
      //   https://cern.ch/twiki/pub/ALICE/PAPaperCentrality/normalization.pdf
      //   https://cern.ch/twiki/bin/view/ALICE/PAMCProductionStudies
      //
      Int_t nsd1=0, nsd2=0, ndd=0;
      Int_t npP = dpm->ProjectileParticipants();
      Int_t npT = dpm->TargetParticipants();
      // Get the numbeer of single and double diffractive participants
      dpm->GetNDiffractive(nsd1,nsd2,ndd);
      // Check if all partipants are single/double diffractive 
      if      ((ndd == 0) && ((npP == nsd1) || (npT == nsd2)))	sd = true;
      else if (ndd == (npP + npT))                       	dd = true;
      Int_t ncp = dpm->NN();
      Int_t nct = dpm->NNw();
      Int_t nwp = dpm->NwN();
      Int_t nwt = dpm->NwNw();
      Int_t nwtacc = _dtglcp->nwtacc;
      Int_t nwtsam = _dtglcp->nwtsam;
      if (first) {
	Printf("@ Npp sd1 Npt sd2  dd tpe Ncp Nct Nwp Nwt acc sam");
	first = false;
      }
      Printf("@ %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d",
	     npP, nsd1, npT, nsd2, ndd, type, ncp, nct, nwp, nwt,
	     nwtacc, nwtsam);
#endif 
    }
    if (gev) fShortHead.fPhiR = gev->GetEventPlane();
    if (herwig) {
      Int_t type = herwig->ProcessType();
      switch (type) {
      case 5: case 6: sd = true; break;
      }
      fShortHead.fNtgt  = 1;
      fShortHead.fNproj = 1;
      fShortHead.fNbin  = 1;
    }
    fShortHead.fType = (sd ? 0x1 : 0) | (dd ? 0x2 : 0);

    // --- Check centrality -----------------------------------------
    Double_t b  = fShortHead.fB;
    Double_t c  = fBEstimator->GetCentrality(b);
    Info("ProcessHeader", "b=%f isProjA=%d isTgtA=%d cms=%f",
	 b, fIsProjA, fIsTgtA, fGenerator ? fGenerator->GetEnergyCMS() : -1);
    if (c >= 0) fShortHead.fC = c;
    
    // --- Check if within vertex cut -------------------------------
    Bool_t selected = (fShortHead.fIpZ <= fHIpz->GetXaxis()->GetXmax() &&
		       fShortHead.fIpZ >= fHIpz->GetXaxis()->GetXmin());

    // --- Only update histograms if within IPz cut ------------------
    if (selected) {
      fHPhiR->Fill(fShortHead.fPhiR*TMath::RadToDeg());
      fHB->Fill(fShortHead.fB);
      fHIpz->Fill(fShortHead.fIpZ);
      if (dd) fHType->Fill(3);
      if (sd) fHType->Fill(2);
      if (!dd && !sd) fHType->Fill(1);
      fHCent->Fill(c);
      fHEta ->AddBinContent(0); // Count events
      fHY   ->AddBinContent(0); // Count events
      fHTime->AddBinContent(0); // Count events
      // fShortHead.Print();
    }
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->ProcessHeader(fShortHead);
    
    return selected;
  }
  /** 
   * Set trigger bits in header 
   * 
   * @param eta Eta of particle 
   */
  virtual void CheckTrigger(Double_t eta)
  {
    if (eta < +5.1 && eta > +2.8)
      fShortHead.fTrigMask |= FastShortHeader::kV0A;
    if (eta < -1.7 && eta > -3.7)
      fShortHead.fTrigMask |= FastShortHeader::kV0C;
    if (eta < +6.3 && eta > +4.8)
      fShortHead.fTrigMask |= FastShortHeader::kADA;
    if (eta < -4.9 && eta > -7.0)
      fShortHead.fTrigMask |= FastShortHeader::kADC;
    if (TMath::Abs(eta) < 1)
      fShortHead.fTrigMask |= FastShortHeader::kEta1;	
  }
  /** 
   * Process all particles 
   * 
   * @param selected True if particle information should be diagnosed
   * 
   * @return true on success
   */
  virtual Bool_t ProcessParticles(Bool_t selected)
  {
    Int_t nPart = fStack->GetNprimary();
    for (Int_t iPart = 0; iPart < nPart; iPart++) {
      TParticle*    particle  = fStack->Particle(iPart);
      TParticlePDG* pdg       = particle->GetPDG();
      Bool_t        primary   = fStack->IsPhysicalPrimary(iPart);
      Bool_t        weakDecay = fStack->IsSecondaryFromWeakDecay(iPart);
      Bool_t        charged   = (pdg &&  TMath::Abs(pdg->Charge()) > 0);
      if (primary)   particle->SetBit(BIT(14));
      if (weakDecay) particle->SetBit(BIT(15));
      if (charged)   particle->SetBit(BIT(16));

      new ((*fParticles)[iPart]) TParticle(*particle);

      TIter next(fCentEstimators);
      FastCentEstimator* estimator = 0;
      while ((estimator = static_cast<FastCentEstimator*>(next())))
	estimator->Process(particle);
      
      if (!selected || !charged || !primary) continue;

      Double_t y = particle->Y();
      if (y > fHY->GetXaxis()->GetXmin() && 
	  y < fHY->GetXaxis()->GetXmax())
	// Avoid filling under/overflow bins 
	fHY->Fill(y);
      
      Double_t pT    = particle->Pt();
      if (pT < 1e-10) continue; /// Along beam axis 
      Double_t pZ    = particle->Pz();
      Double_t theta = TMath::ATan2(pT, pZ);
      Double_t eta   = -TMath::Log(TMath::Tan(theta/2));
      CheckTrigger(eta);
      if (eta > fHEta->GetXaxis()->GetXmin() && 
	  eta < fHEta->GetXaxis()->GetXmax())
	// Avoid filling under/overflow bins 
	fHEta->Fill(eta);
    }
    return true;
  }
  /** 
   * Do final event processing (fill output)
   * 
   */
  virtual void PostEvent()
  {
    fHeader->SetNprimary(fStack->GetNprimary());
    fHeader->SetNtrack(fStack->GetNtrack());

    fHTrig->Fill(0);
    if (fShortHead.fTrigMask & FastShortHeader::kV0A)  fHTrig->Fill(1); 
    if (fShortHead.fTrigMask & FastShortHeader::kV0C)  fHTrig->Fill(2); 
    if (fShortHead.fTrigMask & FastShortHeader::kADA)  fHTrig->Fill(3); 
    if (fShortHead.fTrigMask & FastShortHeader::kADC)  fHTrig->Fill(4); 
    if (fShortHead.fTrigMask & FastShortHeader::kEta1) fHTrig->Fill(5); 
    
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->PostEvent();

    fTree->Fill();
    
    fStack->FinishEvent();
    fHeader->SetStack(fStack);
    
    fRunLoader->TreeE()->Fill();
    fRunLoader->WriteKinematics("OVERWRITE");
  }
  virtual void Generate() { fGenerator->Generate(); }
 
  /** 
   * Process one event 
   * 
   * @param iEv Event number 
   * 
   * @return true on success, false otherwize 
   */
  virtual Bool_t Process(Long64_t iEv)
  {
    // --- The stopwatch ---------------------------------------------
    TStopwatch timer;
    timer.Start();
    PreEvent(iEv);
    fHTime->Fill(1, timer.RealTime());
    
    // --- Generate event --------------------------------------------
    timer.Start();
    Generate();
    fHTime->Fill(2, timer.RealTime());

    // --- Process the header ----------------------------------------
    timer.Start();
    Bool_t selected = ProcessHeader();
    fHTime->Fill(3, timer.RealTime());
    
    // --- Loop over particles ---------------------------------------
    timer.Start();
    ProcessParticles(selected);
    fHTime->Fill(4, timer.RealTime());

    // --- Do final stuff --------------------------------------------    
    timer.Start();
    PostEvent();
    fHTime->Fill(5, timer.RealTime());

    return true;
  }
  virtual void FinishRun()
  {
    fGenerator->FinishRun();
    fRunLoader->WriteHeader("OVERWRITE");
    fGenerator->Write();
    fRunLoader->Write();

  }
  /** 
   * Finalize this sub-job  
   * 
   */
  void SlaveTerminate()
  {
    FinishRun();
    if (fFile) {
      if (fProofFile) {
	if (fVerbose) fProofFile->Print();
	fOutput->Add(fProofFile);
	fOutput->Add(new TH1F("filename", fFileName.Data(),1,0,1));
      }
      // Flush out tree 
      fFile->cd();
      fTree->Write(0, TObject::kOverwrite);
      fFile->Close();
      fFile->Delete();
      fFile = 0;
    }
    if (fAliceFile) {
      TFile* galice = GetGAlice();
      if (galice) {
	if (fVerbose) fAliceFile->Print();
	fAliceFile->AdoptFile(galice);
	fAliceFile->SetOutputFileName(fAliceFile->GetName());
	fOutput->Add(fAliceFile);
	galice->Write();
      }
    }
    if (fKineFile) {
      TFile* kine = GetKine();
      if (kine) {
	if (fVerbose) fKineFile->Print();
	fKineFile->AdoptFile(kine);
	fKineFile->SetOutputFileName(fKineFile->GetName());
	fOutput->Add(fKineFile);
	kine->Write();
      }
    }

    if (fVerbose) {
      Info("SlaveTerminate", "Content of output list");
      gROOT->IncreaseDirLevel();
      fOutput->ls();
      gROOT->DecreaseDirLevel();

      gSystem->Exec("echo \"Content of working directory\"");
      gSystem->Exec("ls -l1 | sed 's/^/  /'");
    }
  }
  /** 
   * Write a collection to disk, transforming sub-collections to
   * directories.
   * 
   * @param c 
   * @param dir 
   */
  void FlushList(TCollection* c, TDirectory* dir)
  {
    dir->cd();
    TIter next(c);
    TObject* o = 0;
    while ((o = next())) {
      if (o->IsA()->InheritsFrom(TCollection::Class())) {
	if (fVerbose) Info("FlushList", "Got collection: %s", c->GetName());
	TDirectory* cur = dir->mkdir(o->GetName());
	FlushList(static_cast<TCollection*>(o), cur);
	dir->cd();
	continue;
      }
      o->Write();
    }
    dir->cd();
  }
      
  /** 
   * Final processing of the data 
   * 
   */
  void Terminate()
  {
    if (gProof) gProof->ClearFeedback();

    if (!fList)
      fList = static_cast<TList*>(fOutput->FindObject("histograms"));
    if (!fList) {
      Error("Terminate", "No output list");
      return;
    }
    
    if (!fProofFile) {
      TObject* fn = fOutput->FindObject("filename");
      if (fn) fFileName  = fn->GetTitle();
      fProofFile =
	static_cast<TProofOutputFile*>(fOutput->FindObject(FileName()));
    }
    if (fProofFile) 
      fFile = fProofFile->OpenFile("UPDATE");
    if (!fFile)
      fFile = TFile::Open(FileName(),"UPDATE");

    TList* estimators = static_cast<TList*>(fList->FindObject("estimators"));
    TList* histos     = static_cast<TList*>(fList->FindObject("histograms"));
    if (!histos) {
      Warning("Terminate", "No histogram list found in output");
      fList->ls();
    }
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->Terminate(estimators);
	
    fHEta  = static_cast<TH1*>(histos->FindObject("dNdeta"));
    fHY    = static_cast<TH1*>(histos->FindObject("dNdy"));
    fHIpz  = static_cast<TH1*>(histos->FindObject("ipZ"));
    fHType = static_cast<TH1*>(histos->FindObject("type"));
    fHCent = static_cast<TH1*>(histos->FindObject("cent"));
    fHB    = static_cast<TH1*>(histos->FindObject("b"));
    fHPhiR = static_cast<TH1*>(histos->FindObject("phiR"));
    fHTime = static_cast<TH1*>(histos->FindObject("timing"));
    fHTrig = static_cast<TH1*>(histos->FindObject("trigger"));

    if (!(fHEta && fHY && fHIpz && fHType && fHB && fHPhiR && fHTime)) {
      Warning("Terminate", "Missing histograms (%p,%p,%p,%p,%p,%p,%p)",
	      fHEta, fHY, fHIpz, fHType, fHB, fHPhiR, fHTime);
      return;
    }

    Int_t nTotal = fHIpz->GetEntries();
    fHEta ->Scale(1./nTotal, "width");
    fHY   ->Scale(1./nTotal, "width");
    fHB   ->Scale(1./nTotal, "width");
    fHPhiR->Scale(1./nTotal, "width");
    fHTime->Scale(1./nTotal, "width");
    fHTrig->Scale(1./nTotal, "width");

    if (!fFile){
      Warning("Terminate", "No file to write to");
      return;
    }

    FlushList(fList, fFile); // ->Write();
    
    fTree = static_cast<TTree*>(fFile->Get("T"));
    if (!fTree)  Warning("Terminate", "No tree");

    if (fVerbose) fFile->ls();
    fFile->Close();

    MoveAliceFiles();
  }
  /** 
   * Retrieve the galice.root file from ROOT 
   * 
   * @return Pointer to file or null
   */
  TFile* GetGAlice()
  {
    if (!fRunLoader) return 0;
    TString galiceName = fRunLoader->GetFileName();
    TFile*  file       = gROOT->GetFile(galiceName);
    if (!file) {
      Warning("GetGAlice", "Didn't find galice file \"%s\"", galiceName.Data());
      gROOT->GetListOfFiles()->ls();
      return 0;
    }
    return file;
  }
  /** 
   * Retrieve the Kinematics.root file from ROOT 
   * 
   * @return Pointer to file or null
   */
  TFile* GetKine()
  {
    if (!fRunLoader) return 0;
    TString kineName   = "Kinematics.root";
    TString galiceName = fRunLoader->GetFileName();
    TString dir        = gSystem->DirName(galiceName);
    if (dir.EqualTo(".")) dir = "";
    if (!dir.IsNull() && dir[dir.Length()-1] != '/') dir.Append("/");
    kineName.Prepend(dir);

    TFile*  file = gROOT->GetFile(kineName);
    if (!file) {
      Warning("GetKine", "Didn't find kinematics file \"%s\"", kineName.Data());
      gROOT->GetListOfFiles()->ls();
      return 0;
    }
    return file;
  }
  /** 
   * Move retrieved ALICE files (galice.root and Kinematics.root) to
   * separate su-directories, and create a collection of the TE tree
   * stored in the galice.root files.
   * 
   */
  void MoveAliceFiles()
  {
    if (!fInput) return;

    TObject* save  = fInput->FindObject("PROOF_SaveGALICE");
    if (!save) return;
    
    TString  sMode = save->GetTitle();
    if (!sMode.EqualTo("split", TString::kIgnoreCase)) return;

    TList*            lst   = new TList;
    TSystemDirectory* dir   = new TSystemDirectory(".",
						   gSystem->WorkingDirectory());
    TList*            files = dir->GetListOfFiles();
    TSystemFile*      file  = 0;
    TIter             next(files);
    while ((file = static_cast<TSystemFile*>(next()))) {
      if (file->IsDirectory()) continue;
      TString fn(file->GetName());
      if (!fn.BeginsWith("galice") && !fn.BeginsWith("Kinematics"))
	continue;

      TPRegexp regex("(.*)_([^_]+)\\.root");
      TObjArray* matches = regex.MatchS(fn);
      if (matches->GetEntriesFast() < 3) {
	delete matches;
	continue;
      }
      TString ord = matches->At(2)->GetName();
      TString bse = matches->At(1)->GetName();

      if (gSystem->AccessPathName(ord,kFileExists))
	gSystem->MakeDirectory(ord);

      if (fVerbose) 
	Info("MoveAliceFiles", "Moving %s to %s/%s.root",
	     fn.Data(), ord.Data(), bse.Data());
      file->Move(Form("%s/%s.root", ord.Data(), bse.Data()));

      if (!bse.EqualTo("galice")) continue;
      TObjString* url = new TObjString(Form("file://%s/%s/%s.root?#TE",
					    file->GetTitle(),
					    ord.Data(),
					    bse.Data()));
      if (fVerbose) 
	Info("MoveAliceFiles", "Adding \"%s\" to file list",
	     url->GetName());
      lst->Add(url);
    }
    if (lst->GetEntries() <= 0) return;
    if (fVerbose) lst->ls();
    
    TFile* out   = TFile::Open("index.root","RECREATE");
    lst->Write("TE",TObject::kSingleKey);
    out->Write();
    out->Close();
    
  }
  /** 
   * Interface version used 
   * 
   * @return 1y
   */
  Int_t Version() const { return 1; }

  /**
   * @{ 
   * @name Parameters 
   */
  TString  fEGName;               // Name of event generator
  Int_t    fRunNo;                // Run to simulate 
  Double_t fBMin;                 // Least impact parameter 
  Double_t fBMax;                 // Largest impact parameter
  TObject* fGRP;                  //! GRP in one line
  TList*   fOverrides;            //! GRP setting to override
  Long64_t fNEvents;              //  Number of requested events
  Bool_t   fIsTgtA;               //! True if target beam is nuclei
  Bool_t   fIsProjA;              //! True if projectile beam is nuclei
  /* @} */
  /** 
   * @{ 
   * @name ALICE EG interface 
   */
  AliGenerator* fGenerator;       //! Event generator
  AliRunLoader* fRunLoader;       //! Loader of trees
  AliStack*     fStack;           //! Stack of particles
  AliHeader*    fHeader;          //! Header handler
  /* @} */
  /** 
   * @{ 
   * @name Custom output 
   */
  TTree*        fTree;            //! Custom tree 
  TClonesArray* fParticles;       //! List of particles
  /**
   * @{ 
   * @name Diagnostics 
   */
  TList* fList;                   //! List of outputs
  TH1*   fHEta;                   //! dN/deta
  TH1*   fHY;                     //! dN/dy
  TH1*   fHIpz;                   //! IPz histogram
  TH1*   fHType;                  //! Event type histogram
  TH1*   fHCent;                  //! Event type histogram
  TH1*   fHB;                     //! B histogram
  TH1*   fHPhiR;                  //! Reaction plane
  TH1*   fHTime;                  //! Timing
  TH1*   fHTrig;                  //! Trigger bits
  /* @} */
  /** 
   * @{ 
   * @name Centrality 
   */
  BCentEstimator* fBEstimator;     // Always present 
  TList* fCentEstimators;          // Centrality estimators
  /* @} */
  /**
   * @{ 
   * @name Output files 
   */
  TProofOutputFile* fProofFile;   //! Proof output file
  TProofOutputFile* fAliceFile;   //! 
  TProofOutputFile* fKineFile;    //! 
  TFile*            fFile;        //! Output file
  mutable TString   fFileName;    //! Output file name 
  /* @} */
  Bool_t fVerbose; // Verbosity 
  Int_t  fMonitor;
  
#ifndef __CINT__
  FastShortHeader fShortHead;
#endif 
  /** 
   * Run this selector as a normal process
   * 
   * @param nev        Number of events
   * @param run        Run number to anchor in
   * @param gen        Generator 
   * @param bMin       Least impact parameter [fm]
   * @param bMax       Largest impact parameter [fm]
   * @param monitor    Monitor frequency [s]
   * @param verbose    Be verbose 
   * @param overrides  GRP overrides 
   *
   * @return true on succes
   */
  static Bool_t  LocalRun(Long64_t       nev,
			  UInt_t         run,
			  const TString& gen,
			  Double_t       bMin,
			  Double_t       bMax,
			  Int_t          monitor,
			  Bool_t         verbose,
			  const TString& overrides="")
			  
  {
    FastSim* sim = new FastSim(gen,run,bMin,bMax,nev,monitor);
    SetOverrides(sim, overrides);
    sim->fVerbose = verbose;
    sim->Begin(0);
    sim->SlaveBegin(0);

    for (Long64_t i=0; i <nev; i++) {
      Printf("=== Event # %6lld/%6lld ==========================",
	     i+1, nev);
      sim->Process(i);
    }
    sim->SlaveTerminate();
    sim->Terminate();

    return true;
  }
  /**
   * Load needed libraries in a proof serssion 
   */
  static void ProofLoadLibs()
  {
    if (!gProof) return;

    // Remember to copy changes to RunFast.C
    TList clsLib;
    clsLib.Add(new TNamed("TVirtualMC",              "libVMC"));
    clsLib.Add(new TNamed("TLorentzVector",          "libPhysics"));
    clsLib.Add(new TNamed("TLinearFitter",           "libMinuit"));
    clsLib.Add(new TNamed("TTree",                   "libTree"));
    clsLib.Add(new TNamed("TProof",                  "libProof"));
    clsLib.Add(new TNamed("TGFrame",                 "libGui"));
    clsLib.Add(new TNamed("TSAXParser",              "libXMLParser"));
    clsLib.Add(new TNamed("AliVEvent",               "libSTEERBase"));
    clsLib.Add(new TNamed("AliESDEvent",             "libESD"));
    clsLib.Add(new TNamed("AliAODEvent",             "libAOD"));
    clsLib.Add(new TNamed("AliAnalysisManager",      "libANALYSIS"));
    clsLib.Add(new TNamed("AliCDBManager",           "libCDB"));
    clsLib.Add(new TNamed("AliRawVEvent",            "libRAWDatabase"));
    clsLib.Add(new TNamed("AliHit",                  "libSTEER"));
    clsLib.Add(new TNamed("AliGenMC",                "libEVGEN"));
    clsLib.Add(new TNamed("AliFastEvent",            "libFASTSIM"));

    TIter next(&clsLib);
    TObject* obj = 0;
    while ((obj = next())) {
      gProof->Exec(Form("gROOT->LoadClass(\"%s\",\"%s\");",
			obj->GetName(), obj->GetTitle()));
    }
  }
  /** 
   * Run this selector in PROOF(Lite)
   * 
   * @param url        Proof URL
   * @param nev        Number of events
   * @param run        Run number to anchor in
   * @param gen        Generator 
   * @param bMin       Least impact parameter [fm]
   * @param bMax       Largest impact parameter [fm]
   * @param monitor    Monitor frequency [s]
   * @param opt        Compilation options
   * @param verbose    Be verbose 
   * @param overrides  GRP overrides 
   * @param save       Where to save
   * 
   * @return true on succes
   */
  static Bool_t ProofRun(const TUrl&    url,
			 Long64_t       nev,
			 UInt_t         run,
			 const TString& gen,
			 Double_t       bMin,
			 Double_t       bMax,
			 Int_t          monitor=-1,
			 Bool_t         verbose=false,
			 const TString& overrides="",
			 const TString& save="none",
			 const char*    opt="")
  {
    TProof::Reset(url.GetUrl());
    TProof::Open(url.GetUrl());
    gProof->ClearCache();

    TString phy = gSystem->ExpandPathName("$(ALICE_PHYSICS)");
    TString ali = gSystem->ExpandPathName("$(ALICE_ROOT)");
    // TString fwd = gSystem->ExpandPathName("$ANA_SRC");
    TString fwd = phy + "/PWGLF/FORWARD/analysis2";

    gProof->AddIncludePath(Form("%s/include", ali.Data()));
    gProof->AddIncludePath(Form("%s/include", phy.Data()));
    ProofLoadLibs();
    gProof->Load(Form("%s/sim/GRP.C",fwd.Data()), true);
    gProof->Load(Form("%s/sim/BaseConfig.C",fwd.Data()), true);
    gProof->Load(Form("%s/sim/EGConfig.C",fwd.Data()), true);

    // gROOT->ProcessLine("gProof->SetLogLevel(5);");
    gProof->Load(Form("%s/sim/FastShortHeader.C", fwd.Data()));
    gProof->Load(Form("%s/sim/FastCentEstimators.C+%s",fwd.Data(),opt));
    gProof->Load(Form("%s/sim/FastMonitor.C+%s",fwd.Data(),opt));
    gProof->Load(Form("%s/sim/FastSim.C+%s", fwd.Data(), opt),true);
    gProof->SetParameter("PROOF_SaveGALICE", save);

    FastSim* sim = new FastSim(gen,run,bMin,bMax,nev,monitor);
    SetOverrides(sim, overrides);
    sim->fVerbose = verbose;
    gProof->Process(sim, nev, "");

    return true; // status >= 0;
  }
  /** 
   * Extract key value pair from string 
   * 
   * @param in  Input string 
   * @param key On return, the key 
   * @param val On return, the value
   * @param sep Separator between key an value 
   * 
   * @return false of separator isn't found in input 
   */
  static Bool_t Str2KeyVal(const TString& in,
			   TString&       key,
			   TString&       val,
			   const char     sep='=')
  {
    Int_t idx = in.Index(sep);
    if (idx == kNPOS) return false;

    key = in(0,idx);
    val = in(idx+1, in.Length()-idx-1);
    return true;
  }
  static void SetOverrides(FastSim* sim, const TString& override)
  {
    if (override.IsNull()) return;

    const char* valid[] = { "beamEnergy", // UInt_t [GeV]
			    "energy",     // UInt_t [GeV]
			    "period",     // String			     
			    "run",        // UInt_t
			    "beam1.a",    // UInt_t
			    "beam1.z",    // UInt_t
			    "beam2.a",    // UInt_t
			    "beam2.z",    // UInt_t
			    0 };
    TObjArray*  tokens = override.Tokenize(",");
    TObjString* token  = 0;
    TIter       next(tokens);
    while ((token = static_cast<TObjString*>(next()))) {
      TString& str = token->String();
      if (str.IsNull()) continue;

      TString  key, val;
      if (!Str2KeyVal(str,key,val, ':')) {
	Printf("Warning: FastSim::Run: incomplete override '%s'",str.Data());
	continue;
      }
      const char** pvalid = valid;
      while (*pvalid) {
	if (key.EqualTo(*pvalid, TString::kIgnoreCase)) {
	  break;
	}
	pvalid++;
      }
      if (!*pvalid) {
	Printf("Warning: FastSim::Run: Invalid override '%s'", key.Data());
	continue;
      }
      // Special case for a string 
      if (key.EqualTo("period",TString::kIgnoreCase))
	val = Form("\"%s\"", val.Data());
      sim->AddOverride(*pvalid, val);
    }
    // delete tokens;
  }
  /** 
   * Run a simulation. 
   * 
   * @a url is the execution URL of the form 
   * 
   * @verbatim 
   PROTOCOL://[HOST[:PORT]]/[?OPTIONS]
   @endverbatim 
   *
   * Where PROTOCOL is one of 
   * 
   * - local for local (single thread) execution 
   * - lite for Proof-Lite execution 
   * - proof for Proof exection 
   * 
   * HOST and PORT is only relevant for Proof. 
   *
   * Options is a list of & separated options 
   * 
   * - events=NEV  Set the number of events to process 
   * - run=RUNNO   Set the run number to anchor in 
   * - eg=NAME     Set the event generator 
   * - b=RANGE     Set the impact parameter range in fermi 
   * - monitor=SEC Set the update rate in seconds of monitor histograms 
   *
   * 
   * @param url Exection URL 
   * @param opt Optimization used when compiling 
   * 
   * @return true on success 
   */
  static Bool_t Run(const char*  url, const char* opt="")
  {
    Printf("Will run fast simulation with:\n\n\t%s\n\n",url);
    Long64_t     nev     = 10000;
    UInt_t       run     = 0;
    TString      eg      = "default";
    TString      override= "";
    TString      save    = "none";
    Double_t     bMin    = 0;
    Double_t     bMax    = 20;
    Int_t        monitor = -1;
    Bool_t       verbose = false;
    TUrl         u(url);
    TString      out;
    TObjArray*   opts    = TString(u.GetOptions()).Tokenize("&");
    TObjString*  token   = 0;
    TIter        nextToken(opts);
    while ((token = static_cast<TObjString*>(nextToken()))) {
      TString& str = token->String();
      if (str.IsNull()) continue;

      if (str.EqualTo("verbose")) { verbose = true; str = ""; }

      if (str.IsNull()) continue;
      
      TString  key, val;
      if (!Str2KeyVal(str,key,val)) {
	if (!out.IsNull()) out.Append("&");
	out.Append(str);
	continue;
      }

      if      (key.EqualTo("events"))   nev      = val.Atoll();
      else if (key.EqualTo("run"))      run      = val.Atoi();
      else if (key.EqualTo("eg"))       eg       = val;
      else if (key.EqualTo("override")) override = val;
      else if (key.EqualTo("save"))     save     = val;
      else if (key.EqualTo("monitor"))  monitor  = val.Atoi();
      else if (key.EqualTo("b")) {
	TString min, max;
	if (Str2KeyVal(val, min, max, '-')) {
	  bMin = min.Atof();
	  bMax = max.Atof();
	}
      }
      else {
	if (!out.IsNull()) out.Append("&");
	out.Append(str);
      }
    }
    opts->Delete();
    u.SetOptions(out);
    if (!u.IsValid()) {
      Printf("Error: FastSim::Run: URL %s is invalid", u.GetUrl());
      return false;
    }
    
    Bool_t isLocal = TString(u.GetProtocol()).EqualTo("local");

    Printf("Run %s for %lld events anchored at %d\n"
	   "  Impact paramter range:  %5.1f-%5.1f fm\n"
	   "  Monitor frequency:      %d sec\n"
	   "  Execution url:          %s",
	   eg.Data(), nev, run, bMin, bMax, monitor, u.GetUrl());


    TStopwatch timer;
    timer.Start();

    Bool_t ret = false;
    if (isLocal)
      ret = LocalRun(nev, run, eg, bMin, bMax, monitor, verbose, override);
    else 
      ret = ProofRun(u, nev, run, eg, bMin, bMax,
		     monitor, verbose, override, save, opt);
    timer.Print();

    return ret;
  }
		    
  ClassDef(FastSim,3); 
};


struct EPosSim : public FastSim
{
  EPosSim(UInt_t run=0, Int_t monitor=0)
    : FastSim("epos", run, 0, 20, 100000, monitor),
      fInTree(0),
      fInNTot(0),
      fInB(0),
      fInPDG(0),
      fInStatus(0),
      fInPx(0),
      fInPy(0),
      fInPz(0),
      fInE(0),
      fInM(0),
      fInNcollH(0),
      fInNpartP(0),
      fInNpartT(0),
      fInNcoll(0),
      fInNSpcPN(0),
      fInNSpcTN(0),
      fInNSpcPP(0),
      fInNSpcTP(0),
      fInPhiR(0)
  {
  }
  Bool_t SetupBranches()
  {
    fInNTot   = fInTree->GetLeaf("nPart");
    fInB      = fInTree->GetLeaf("ImpactParameter");
    fInPDG    = fInTree->GetLeaf("pdgid");
    fInStatus = fInTree->GetLeaf("status");
    fInPx     = fInTree->GetLeaf("px");
    fInPy     = fInTree->GetLeaf("py");
    fInPz     = fInTree->GetLeaf("pz");
    fInE      = fInTree->GetLeaf("E");
    fInM      = fInTree->GetLeaf("m");
    // These are probably EPOS-LHC specific 
    fInNcollH = fInTree->GetLeaf("Ncoll_hard");
    fInNpartP = fInTree->GetLeaf("Npart_proj");
    fInNpartT = fInTree->GetLeaf("Npart_targ");
    fInNcoll  = fInTree->GetLeaf("Ncoll");
    fInNSpcPN = fInTree->GetLeaf("Nspec_proj_neut");
    fInNSpcTN = fInTree->GetLeaf("Nspec_targ_neut");
    fInNSpcPP = fInTree->GetLeaf("Nspec_proj_prot");
    fInNSpcTP = fInTree->GetLeaf("Nspec_targ_prot");
    fInPhiR   = fInTree->GetLeaf("phiR");
    
    

    return (fInNTot && fInB && fInPDG && fInPx && fInPy && fInPz && fInE);
  }
  void Begin(TTree* tree)
  {
    FastSim::Begin(tree);
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next()))) {
      if (!estimator->IsA()->InheritsFrom(V0CentEstimator::Class()))
	continue;
      V0CentEstimator* v = static_cast<V0CentEstimator*>(estimator);
      v->Flip(false); // Flip detector acceptance of V0A/C
      // v->Print("nah");
    }
  }
  void Init(TTree* tree)
  {
    Info("Init", "Initializing with tree %p (%s)",
	 tree, (tree ? tree->ClassName() : ""));
    if (!tree) return;

    TFile* file = tree->GetCurrentFile();
    Info("Init", "Current file: (%p) %s", file,
	 (file ? file->GetName() : ""));
    
    fInTree = tree;
    if (!SetupBranches())
      Fatal("Init", "Failed to set-up branches");
    // if (!SetupEstimator())
    //   Fatal("Init", "Failed to set-up estimator");
  }    
  /** 
   * Called when the file changes 
   * 
   * @return true on success, false otherwise 
   */
  Bool_t Notify()
  {
    if (!fInTree) {
      Warning("Notify", "No tree set yet!");
      return false;
    }
    TFile* file = fInTree->GetCurrentFile();
    Info("Notify", "processing file: (%p) %s", file,
	 (file ? file->GetName() : ""));
    if (!file) return true;
    if (!SetupBranches()) {
      Warning("Notify", "Failed to set-up branches");
      return false;
    }
    return true;
  }
  const char* GetEGTitle() const
  {
    static TString tmp;
    if (!tmp.IsNull()) return tmp.Data();
    tmp = "EPOS-LHC ";
    Long_t ret =
      gROOT->ProcessLine("Form(\"%s(%d,%d)+%s(%d,%d) @ %5d b in[%4.1f,%4.1f]\","
			 "grp->beam1.Name(), grp->beam1.a, grp->beam1.z,"
			 "grp->beam2.Name(), grp->beam2.a, grp->beam2.z,"
			 "Int_t(grp->energy))");
    tmp.Append(reinterpret_cast<const char*>(ret));
    Printf("EG title set to %s", tmp.Data());
    return tmp.Data();
  }
  void SetupSeed() {}
  Bool_t SetupGen()
  {
    fIsTgtA  = gROOT->ProcessLine("grp->beam1.IsA()");
    fIsProjA = gROOT->ProcessLine("grp->beam2.IsA()");
    return true;
    
  }
  Bool_t SetupRun() { return true; }
  Bool_t PreEvent(Long64_t iEv)
  {
    Int_t read = fInTree->GetTree()->GetEntry(iEv);
    if (read <= 0) return false;
    
    // --- Reset header ----------------------------------------------
    fShortHead.Reset(fRunNo, iEv);
    fParticles->Clear();
    // Reset input

    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->PreEvent();
    
    return true;    
  }
  Bool_t ProcessHeader()
  {
    fShortHead.fIpX        = 0;
    fShortHead.fIpY        = 0;
    fShortHead.fIpZ        = 0;
    fShortHead.fB          = fInB->GetValue();
    fShortHead.fNtgt       = (fInNpartT ? fInNpartT->GetValue() : 0);
    fShortHead.fNproj      = (fInNpartP ? fInNpartP->GetValue() : 0);
    fShortHead.fNbin       = (fInNcoll  ? fInNcoll ->GetValue() : 0);
    fShortHead.fPhiR       = (fInPhiR   ? fInPhiR  ->GetValue() : 0);
    fShortHead.fNSpecNproj = (fInNSpcPN ? fInNSpcPN->GetValue() : 0);
    fShortHead.fNSpecNtgt  = (fInNSpcTN ? fInNSpcTN->GetValue() : 0);
    fShortHead.fNSpecPproj = (fInNSpcPP ? fInNSpcPP->GetValue() : 0);
    fShortHead.fNSpecPtgt  = (fInNSpcTP ? fInNSpcTP->GetValue() : 0);
    fShortHead.fEG         = FastShortHeader::kEPOS;
    
    Double_t c = fBEstimator->GetCentrality(fShortHead.fB);
    if (c >= 0) fShortHead.fC = c;

    // --- Check if within vertex cut -------------------------------
    Bool_t selected = (fShortHead.fIpZ <= fHIpz->GetXaxis()->GetXmax() &&
		       fShortHead.fIpZ >= fHIpz->GetXaxis()->GetXmin());

    // --- Only update histograms if within IPz cut ------------------
    if (selected) {
      fHPhiR->Fill(fShortHead.fPhiR*TMath::RadToDeg());
      fHB->Fill(fShortHead.fB);
      fHIpz->Fill(fShortHead.fIpZ);
      fHCent->Fill(c);
      // fShortHead.Print();
    }
    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->ProcessHeader(fShortHead);
    return selected;
  }
  virtual Bool_t ProcessParticles(Bool_t selected)
  {
    Int_t nTot = fInNTot->GetValue();
    for (Int_t iPart = 0; iPart < nTot; iPart++) {
      Int_t    status = fInStatus->GetValue(iPart);
      Int_t    pdg    = fInPDG->GetValue(iPart);
      Double_t px     = fInPx->GetValue(iPart);
      Double_t py     = fInPy->GetValue(iPart);
      Double_t pz     = fInPz->GetValue(iPart);
      // Double_t pz     = -fInPz->GetValue(iPart); // flip sign on pZ
      Double_t e      = fInE->GetValue(iPart);
      // Double_t m      = fInM->GetValue(iPart);

      TParticle* particle =
	new ((*fParticles)[iPart]) TParticle(pdg, status,-1,-1,-1,-1,
					     px, py, pz, e, 0, 0, 0, 0);
      TParticlePDG* pdgP      = particle->GetPDG();
      Bool_t        primary   = status == 1;
      Bool_t        weakDecay = false;
      Bool_t        charged   = (pdgP &&  TMath::Abs(pdgP->Charge()) > 0);
      if (primary)   particle->SetBit(BIT(14));
      if (weakDecay) particle->SetBit(BIT(15));
      if (charged)   particle->SetBit(BIT(16));

      TIter next(fCentEstimators);
      FastCentEstimator* estimator = 0;
      while ((estimator = static_cast<FastCentEstimator*>(next())))
	estimator->Process(particle);
      
      if (!selected || !charged || !primary) continue;
      fHY  ->Fill(particle->Y());
      Double_t pT    = particle->Pt();
      if (pT < 1e-10) continue; /// Along beam axis 
      Double_t pZ    = particle->Pz();
      Double_t theta = TMath::ATan2(pT, pZ);
      Double_t eta   = -TMath::Log(TMath::Tan(theta/2));
      CheckTrigger(eta);
      fHEta->Fill(eta);
    }
    return true;
  }
  void PostEvent()
  {
    fTree->Fill();

    TIter next(fCentEstimators);
    FastCentEstimator* estimator = 0;
    while ((estimator = static_cast<FastCentEstimator*>(next())))
      estimator->PostEvent();
  }    
  void Generate() {}
  void FinishRun() {}
  TTree* fInTree;   //!
  TLeaf* fInNTot;   //!
  TLeaf* fInB;      //!
  TLeaf* fInPDG;    //!
  TLeaf* fInStatus; //! 
  TLeaf* fInPx;     //!
  TLeaf* fInPy;     //! 
  TLeaf* fInPz;     //! 
  TLeaf* fInE;      //! 
  TLeaf* fInM;      //!
  TLeaf* fInNcollH; //!
  TLeaf* fInNpartP; //!
  TLeaf* fInNpartT; //!
  TLeaf* fInNcoll;  //!
  TLeaf* fInNSpcPN; //!
  TLeaf* fInNSpcTN; //!
  TLeaf* fInNSpcPP; //!
  TLeaf* fInNSpcTP; //!
  TLeaf* fInPhiR;   //!
  
  /** 
   * Run this selector as a normal process
   * 
   * @param nev        Number of events
   * @param run        Run number to anchor in
   * @param monitor    Monitor frequency [s]
   * @param verbose    Whether to be verbose 
   * 
   * @return true on succes
   */
  static Bool_t  LocalRun(Long64_t       nev,
			  UInt_t         run,
			  Int_t          monitor,
			  Bool_t         verbose)
			  
  {
    EPosSim* sim = new EPosSim(run, monitor);
    sim->fVerbose = verbose;
    sim->Begin(0);
    sim->SlaveBegin(0);

    for (Long64_t i=0; i <nev; i++) {
      Printf("=== Event # %6lld/%6lld ==========================",
	     i+1, nev);
      sim->Process(i);
    }
    sim->SlaveTerminate();
    sim->Terminate();

    return true;
  }
  /** 
   * Run this selector in PROOF(Lite)
   * 
   * @param url        Proof URL
   * @param opt        Compilation options
   * 
   * @return true on succes
   */
  static Bool_t SetupProof(const TUrl&    url,
			   const char*    opt="")
  {
    TProof::Reset(url.GetUrl());
    TProof::Open(url.GetUrl());
    gProof->ClearCache();

    TString phy = gSystem->ExpandPathName("$(ALICE_PHYSICS)");
    TString ali = gSystem->ExpandPathName("$(ALICE_ROOT)");
    // TString fwd = gSystem->ExpandPathName("$ANA_SRC");
    TString fwd = phy + "/PWGLF/FORWARD/analysis2";

    gProof->AddIncludePath(Form("%s/include", ali.Data()));
    gProof->AddIncludePath(Form("%s/include", phy.Data()));
    ProofLoadLibs();
    gProof->Load(Form("%s/sim/GRP.C",fwd.Data()), true);
    gProof->Load(Form("%s/sim/BaseConfig.C",fwd.Data()), true);
    gProof->Load(Form("%s/sim/EGConfig.C",fwd.Data()), true);

    // gROOT->ProcessLine("gProof->SetLogLevel(5);");
    gProof->Load(Form("%s/sim/FastMonitor.C+%s",fwd.Data(),opt));
    gProof->Load(Form("%s/sim/FastShortHeader.C", fwd.Data()));
    gProof->Load(Form("%s/sim/FastCentEstimators.C+%s",fwd.Data(),opt));
    gProof->Load(Form("%s/sim/FastSim.C+%s", fwd.Data(), opt),true);

    return true; // status >= 0;
  }
  /** 
   * Run a simulation. 
   * 
   * @a url is the execution URL of the form 
   * 
   * @verbatim 
   PROTOCOL://[HOST[:PORT]]/[?OPTIONS]
   @endverbatim 
   *
   * Where PROTOCOL is one of 
   * 
   * - local for local (single thread) execution 
   * - lite for Proof-Lite execution 
   * - proof for Proof exection 
   * 
   * HOST and PORT is only relevant for Proof. 
   *
   * Options is a list of & separated options 
   * 
   * - events=NEV  Set the number of events to process 
   * - run=RUNNO   Set the run number to anchor in 
   * - eg=NAME     Set the event generator 
   * - b=RANGE     Set the impact parameter range in fermi 
   * - monitor=SEC Set the update rate in seconds of monitor histograms 
   *
   * 
   * @param url Exection URL 
   * @param opt Optimization used when compiling 
   * 
   * @return true on success 
   */
  static Bool_t Run(const char*  url, const char* opt="")
  {
    Printf("Will run fast simulation with:\n\n\t%s\n\n",url);
    UInt_t       run     = 0;
    Long64_t     nev     = -1;
    Int_t        monitor = -1;
    Bool_t       verbose = false;
    TUrl         u(url);
    TString      out;
    TObjArray*   opts    = TString(u.GetOptions()).Tokenize("&");
    TObjString*  token   = 0;
    TIter        nextToken(opts);
    while ((token = static_cast<TObjString*>(nextToken()))) {
      TString& str = token->String();
      if (str.IsNull()) continue;

      if (str.EqualTo("verbose")) { verbose = true; str = ""; }

      if (str.IsNull()) continue;
      
      TString  key, val;
      if (!Str2KeyVal(str,key,val)) {
	if (!out.IsNull()) out.Append("&");
	out.Append(str);
	continue;
      }

      if      (key.EqualTo("run"))      run      = val.Atoi();
      else if (key.EqualTo("monitor"))  monitor  = val.Atoi();
      else if (key.EqualTo("events"))   nev      = val.Atoi();
      else {
	if (!out.IsNull()) out.Append("&");
	out.Append(str);
      }
    }
    opts->Delete();
    u.SetOptions(out);
    if (!u.IsValid()) {
      Printf("Error: FastSim::Run: URL %s is invalid", u.GetUrl());
      return false;
    }
    
    Printf("Run EPos anchored at %d\n"
	   "  Monitor frequency:      %d sec\n"
	   "  Execution url:          %s",
	   run, monitor, u.GetUrl());

    TString treeName = u.GetAnchor();
    if (treeName.IsNull()) treeName = "Particle";
    TFile*  file = TFile::Open(u.GetFile(), "READ");
    if (!file) {
      Printf("Error: FastAnalysis::Run: Failed to open %s",
	     u.GetFile());
      return false;
    }

    TChain*   chain = new TChain(treeName, treeName);
    TObject*  o = file->Get(treeName);
    if (!o) {
      Printf("Error: FastAnalysis::Run: Couldn't get %s from %s",
	     treeName.Data(), u.GetFile());
      file->Close();
      return false;
    }
    Int_t cret = 0;
    if (o->IsA()->InheritsFrom(TChain::Class()))
      cret = chain->Add(static_cast<TChain*>(o));
    else if (o->IsA()->InheritsFrom(TTree::Class()))
      cret = chain->AddFile(u.GetFile());
    else if (o->IsA()->InheritsFrom(TCollection::Class()))
      cret = chain->AddFileInfoList(static_cast<TCollection*>(o));
    else if (o->IsA()->InheritsFrom(TFileCollection::Class()))
      cret = chain->AddFileInfoList(static_cast<TFileCollection*>(o)
				    ->GetList());
    file->Close();
    if (cret <= 0 || chain->GetListOfFiles()->GetEntries() <= 0) {
      Printf("Error: FastAnalysis::Run: Failed to create chain");
      return false;
    }

    TString       proto    = u.GetProtocol();
    Bool_t        isProof  = (proto.EqualTo("proof") || proto.EqualTo("lite"));
    if (isProof) {
      if (!SetupProof(u,opt)) return false;
      chain->SetProof();
    }

    EPosSim* sim = new EPosSim(run, monitor);
    sim->fVerbose = verbose;
    if (nev < 0) nev = TChain::kBigNumber;
    
    TStopwatch timer;
    timer.Start();

    Long64_t ret = chain->Process(sim, "", nev, 0);
    timer.Print();

    return ret > 0;
  }
  Int_t Version() const { return 2; }
  ClassDef(EPosSim,1);
};

#endif
//
// EOF
// 
