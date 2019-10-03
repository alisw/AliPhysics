#include "TrainSetup.C"
#ifndef __CINT__
#include <AliESDInputHandlerRP.h>
// #include <AliMixInputEventHandler.h>
// #include <AliMixEventPool.h>
// #include <AliMixEventCutObj.h>
#include <TMacro.h>
#else
class AliAnalysisTaskSE;
#endif

/** 
 * Train definition to do Rubens dN/deta analysis in centrality bins. 
 * 
 * Run as 
 @verbatim 
 runTrain --class=TrackletTrain --name=NAME --url=URL 
 @endverbatim 
 *
 * See also 
 *
 * - https://twiki.cern.ch/twiki/bin/view/ALICE/FMDTrains
 * - http://hehi00.nbi.dk:8888/pwglfforward/train_setup_doc.html
 * 
 * This train uses a custom PAR file (RubensCode.par) of the classes
 *
 * - AliTrackletTaskMulti
 * - AliITSMultRecBg
 *
 * If these classes lived in a compiled AliPhysics library (Say
 * libPWGUD.so), the we wouldn't need that PAR file.
 */
struct MakeTrackletTrain : public TrainSetup
{
  /** 
   * Constructor.  This sets up the available options 
   * 
   * @param name Free form name 
   */
  MakeTrackletTrain(const char* name)
    : TrainSetup(name)
  {
    // Define all our options here
    fOptions.Add("trig",          "NAME",  "Trigger to use",           "V0AND");
    fOptions.Add("cent",          "METHOD","Centrality selector",       "V0M");
    fOptions.Add("cent-bins",     "Bins",  "Centrality bins",
		 "0-5-10-20-30-40-50-60-70-80");
    fOptions.Add("check-reco",             "Check reconstruction",      false);
    fOptions.Add("sig-n-dev",     "X",     "Cut on weighted distance",  25.);
    fOptions.Add("scale-dtheta",           "Scale dTheta",              true);
    fOptions.Add("cut-dtheta",             "Cut dTheta",                false);
    fOptions.Add("dphi",          "X",     "dPhi window",               0.06);
    fOptions.Add("dtheta",        "X",     "dTheta window",             0.025);
    fOptions.Add("phi-shift",     "X",     "Bending shift",             0.0045);
    fOptions.Add("outlier-phi",   "X",     "Phi outlier cut",           0.005);
    fOptions.Add("outlier-z-eta", "X",     "Z-Eta outlier cut",         0.05);
    fOptions.Add("phi-rot",       "X",     "Rotation BG angle",   TMath::Pi());
    fOptions.Add("inj-scale",     "X",     "Injection BG scale",        1.);
    fOptions.Add("remove-outliers",        "Whether to remove outliers",true);
    fOptions.Add("sig-dphi-s",    "X",     "Cut on dPhi-phiBent",       -1.);
    fOptions.Add("sig-n-std",     "X",     "Cut on weighted distance",  1.5);
    fOptions.Add("mc-v0-scale",   "X",     "Scaling of MC V0",          1.);
    fOptions.Add("eta-min",       "ETA",   "Least eta",                 -2);
    fOptions.Add("eta-max",       "ETA",   "Largest eta",               +2);
    fOptions.Add("deta",          "WIDTH", "Eta bin width",             0.25);
    fOptions.Add("ipz-min",       "CM",    "Least IPz",                 -15);
    fOptions.Add("ipz-max",       "CM",    "Largest IPz",               +15);
    fOptions.Add("fill-reco",              "Fill with new reco",        true);
    fOptions.Add("create-inj",             "Create injection BG",       true);
    fOptions.Add("create-rot",             "Create rotation BG",        false);
    fOptions.Add("reweight",      "OPTIONS","','-seperates list",       "");
    fOptions.SetDescription("Run tracklet analysis on ESD data files to "
			    "extract dNch/deta. The train should be run on "
			    "both real data and simulated data. When "
			    "processing simulated data, it is possible to "
			    "reweigh the particles according to predefined "
			    "histograms in pT.  These histograms should "
			    "be stored in separate files. The option to "
			    "select reweighting is 'reweight', which is "
			    "a comma separated list of options:\n"
			    "- pt, pid, str: Reweigh according to pT, PID or "
			    "weak decays\n"
			    "- up, down: Vary reweight up or down\n"
			    "- pi, K, p: For PID reweigh of particular specie");
    
  }
  /** 
   * Create the input handler.  This is overwritten from the base
   * class to allow using AliESDInputHandlerRP for rec. points., and
   * AliMixInputEventHandler if requested.
   * 
   * @param type Type of analysis 
   * @param needRec If we need rec-points
   * 
   * @return The input handler 
   */
  AliVEventHandler* CreateInputHandler(UShort_t type, Bool_t needRec=false)
  {
    
    Bool_t   fill_reco       = fOptions.AsBool("fill-reco");;
    Bool_t   create_inj      = fOptions.AsBool("create-inj");;
    Bool_t   create_rot      = fOptions.AsBool("create-rot");;
    needRec = fill_reco || create_inj || create_rot || needRec;

    Info("CreateInputHandler",
	 "fill_reco=%d create_inj=%d create_rot=%d needRec=%d",
	 fill_reco, create_inj, create_rot, needRec);
    return TrainSetup::CreateInputHandler(type, needRec);
  }
  /** 
   * Create the MC input handler.  Overwritten here to allow setting
   * the pre-read mode.
   * 
   * @param type Input type 
   * @param mc   True for MC 
   * 
   * @return The MC input handler 
   */
  AliVEventHandler* CreateMCHandler(UShort_t type, bool mc)
  {
    AliMCEventHandler* ret =
      static_cast<AliMCEventHandler*>(TrainSetup::CreateMCHandler(type,mc));
    if (ret) ret->SetPreReadMode(AliMCEventHandler::kLmPreRead);
    return ret;
  }
  /** 
   * Create output handler.  Overloaded here to set no output handler,
   * as this train does not define any AOD output.
   * 
   * @return Always null
   */
  virtual AliVEventHandler* CreateOutputHandler(UShort_t)
  {
    return 0;
  }
  /** 
   * Analyse the reweighting option string and set options on task
   * appropriately.  The string is a comma or space separated list of
   * what to reweigh and how to do it. 
   *
   * What to reweigh can be specfied using one or more of the strings 
   *
   * - pt  Reweight in pT 
   * - pid Reweight particle abundance of pi, K, proton 
   * - str Reweight particles from strange weak decays 
   *
   * How to reweigh can be specifed as 
   *
   * - + or up   Increase weights (for pt < 0.05 by +30%)
   * - - or down Decrease weights (for pt < 0.05 by -30%)
   * - If none of these are given, then the weights are used as is. 
   * 
   * If pid rewighting is done and one of up or down are given, then
   * one can specify which particle type to reweigh
   *
   * - pi or pion    Reweight (up or down) pions 
   * - K  or kaon    Reweight (up or down) kaons
   * - p  or proton  Reweight (up or down) protons 
   *
   * Note, if PID, with explicit selection of pions, and strangeness
   * reweighting are specified, then the up/down flag applies to both
   * PID reweighting and the strangeness reweighting
   * 
   * @param task The task to modify 
   */
  void SetupReweighting(AliAnalysisTaskSE* task)
  {
    TString sel = fOptions.AsString("reweight");
    sel.ToLower();
    if (sel.IsNull() || sel.BeginsWith("no"))
      return;

    TList       files;
    Int_t       what = 0;
    Int_t       opt  = 0;
    TObjArray*  tokens = sel.Tokenize(", ");
    TIter       next(tokens);
    TObjString* ostr;
    // First find what should be done 
    while ((ostr = static_cast<TObjString*>(next()))) {
      const TString& token = ostr->String();
      
      if      (token.EqualTo("pt"))   {
	what |= 0x1;
	files.Add(new TObjString("REWEIGHTpt.root"));
	Printf("Will reweigh in pT");
      }
      else if (token.EqualTo("pid")) {
	what |= 0x2;
	Printf("Will reweigh particle species");
      }
      else if (token.EqualTo("str"))  {
	what |= 0x4;
	Printf("Will reweight particles from strange weak decays");
      }
    }
    if (what == 0x0) return;
    
    // Now figure out how to do it 
    next.Reset();
    TString part;
    while ((ostr = static_cast<TObjString*>(next()))) {
      const TString& token = ostr->String();
      Int_t aOpt = TMath::Abs(opt);
      if      (token.EqualTo("up")   || token.EqualTo("+")) 
	opt = (aOpt==0 ? +1 : +aOpt); 
      else if (token.EqualTo("down") || token.EqualTo("-"))
	opt = (aOpt==0 ? -1 : -aOpt);
      else if (token.EqualTo("pi")   || token.EqualTo("pion")){
	  opt = 1; part = "pi";
      }
      else if (token.EqualTo("k")    || token.EqualTo("kaon")) {
	opt = 2; part = "ka";
      }
      else if (token.EqualTo("p")    || token.EqualTo("proton")) {	  
	opt = 3; part = "pr";
      }
    }
    if (opt != 0)
      Printf("Will reweigh %s (%c30%% for pT<0.05)",
	     opt < 0 ? "down" : "up", opt < 0 ? '-' : '+');
    if (what & 0x2) {
      if (!part.IsNull()) {
	Printf("Will reweight %s in particular", part.Data());
	part.Prepend("_");
	part.Append(opt < 0 ? "-" : "+");
      }
      files.Add(new TObjString(Form("REWEIGHTpid%s.root", part.Data())));
    }
    if (what & 0x4)
      files.Add(new TObjString(Form("REWEIGHTstr%s.root",
				    opt == -1 ? "-" :
				    opt == +1 ? "+" : "")));
    delete tokens;

    Printf("Setting reweighing flag=0x%x with option=%d", what, opt);
    SetOnTask(task, "ReweightStack", what);
    SetOnTask(task, "ReweightFlag",  opt);

    TIter nextF(&files);
    while ((ostr = static_cast<TObjString*>(nextF()))) {
      Printf("Loading reweighting file %s", ostr->GetName());
      fRailway->LoadAux(ostr->GetName());
    }
  }
  /** 
   * Create our task, and return it.  This uses the interpreter to
   * make the object.  
   * 
   * @param mgr Analysis manager 
   * 
   * @return Pointer to the task 
   */
  AliAnalysisTaskSE* CreateTask(AliAnalysisManager* mgr)
  {
    // Enable these lines to load the code from PWGUD directory.
    // These do not seem to be up-to-speed with the latest
    // developments, so for now, we use private scripts - sigh!
    //
    // Note, PWGLF also has these scripts, but it is not clear that
    // they are anymore current than the ones in PWGUD.
    // 
    // gROOT->SetMacroPath(Form("%s:$ALICE_PHYSICS/PWGUD/multVScentPbPb",
    // gROOT->GetMacroPath()));
    // gSystem->AddIncludePath("-I$ALICE_PHYSICS/PWGUD/multVScentPbPb");
    Info("CreateTasks", "Loading code");
    fRailway->LoadSource("FixPaths.C");
    // fRailway->LoadAux("AliITSMultRecBg.h");
    fRailway->LoadAux("AliTrackletTaskMulti.h");
    // fRailway->LoadSource("AliITSMultRecBg.cxx");
    fRailway->LoadSource("AliTrackletTaskMulti.cxx");

    // --- Create the task using interpreter -------------------------
    Bool_t             mc  = mgr->GetMCtruthEventHandler() != 0;
    const char*        cls = "AliTrackletTaskMulti";
    Long_t             ret = gROOT->ProcessLine(Form("new %s(\"%s\")",cls,cls));
    AliAnalysisTaskSE* task =reinterpret_cast<AliAnalysisTaskSE*>(ret);
    if (!task) return 0;

    // --- Add task to train -----------------------------------------
    mgr->AddTask(task);

    // --- Figure out the trigger options ----------------------------
    TString trg = fOptions.Get("trig");
    trg.ToUpper();
    UInt_t  sel = AliVEvent::kINT7;
    if      (trg.EqualTo("MB"))    sel = AliVEvent::kMB;
    else if (trg.EqualTo("V0AND")) sel = AliVEvent::kINT7;
    else if (trg.EqualTo("V0OR"))  sel = AliVEvent::kCINT5;
    else if (trg.EqualTo("ANY"))   sel = AliVEvent::kAny;

    // --- Fix up some options ---------------------------------------
    Double_t sig_dphi_s      = fOptions.AsDouble("sig-dphi-s",-1);
    Double_t sig_n_std       = fOptions.AsDouble("sig-n-std",1.5);
    Double_t dphi            = fOptions.AsDouble("dphi",0.06);
    if (sig_dphi_s<0) fOptions.Set("sig-dphi-s", TMath::Sqrt(sig_n_std)*dphi);

    // --- Set various options on task -------------------------------
    task->SelectCollisionCandidates(sel);
    SetOnTask(task, "UseMC", mc);
    SetOnTask(task, "TriggerSelection", sel);
    FromOption(task, "UseCentralityVar", 	"cent", 	"");
    FromOption(task, "CheckReconstructables",   "check-reco",	false);
    FromOption(task, "NStdDev",			"sig-n-dev",	25.);
    FromOption(task, "ScaleDThetaBySin2T",	"scale-dtheta",	false);
    FromOption(task, "CutOnDThetaX",	        "cut-dtheta",	false);
    FromOption(task, "PhiWindow",		"dphi",		0.06);
    FromOption(task, "ThetaWindow",		"dtheta",	0.025);
    FromOption(task, "PhiShift",		"phi-shift",	0.0045);
    FromOption(task, "PhiOverlapCut",		"outlier-phi",	0.005);
    FromOption(task, "ZetaOverlapCut",		"outlier-z-eta",0.05);
    FromOption(task, "PhiRot",			"phi-rot",	TMath::Pi());
    FromOption(task, "InjScale",		"inj-scale",	1.);
    FromOption(task, "RemoveOverlaps",		"remove-outliers",false);
    FromOption(task, "DPhiSCut",		"sig-dphi-s",	-1.);
    FromOption(task, "NStdCut",			"sig-n-std",	1.5);
    FromOption(task, "ScaleMCV0",		"mc-v0-scale",	1.);
    FromOption(task, "EtaMin",			"eta-min",	-2.);
    FromOption(task, "EtaMax",			"eta-max",	+2.);
    FromOption(task, "EtaBinWidth",		"deta",		+2.5);
    FromOption(task, "ZVertexMin",		"ipz-min",	-15.);
    FromOption(task, "ZVertexMax",		"ipz-max",	+15.);
    FromOption(task, "DoNormalReco",		"fill-reco",	false);
    FromOption(task, "DoInjection",		"create-inj",	false);
    FromOption(task, "DoRotation",		"create-rot",	false);
    SetupReweighting(task);

    // --- Set centrality bins ---------------------------------------
    TString centBins = fOptions.AsString("cent-bins");
    TObjArray* tokens = centBins.Tokenize("-");
    TArrayD array(tokens->GetEntries());
    for (Int_t i = 0; i < array.GetSize(); i++) {
      TObjString* ostr = static_cast<TObjString*>(tokens->At(i));
      TString&    str  = ostr->String();
      array[i]         = str.Atof();
    }
    gROOT->ProcessLine(Form("((%s)%p)->SetCentPercentiles((Double_t*)%p,%d)",
			    cls, task, array.GetArray(), array.GetSize()-1));

    // --- Connect I/O -----------------------------------------------
    AliAnalysisDataContainer *out =
      mgr->CreateContainer("clist", TList::Class(),
			   AliAnalysisManager::kOutputContainer,
			   (mc ? "trmc.root" : "trdt.root"));
    mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,out);

    return task;
  }
  /** 
   * Create our tasks.  
   * 
   * @param mgr 
   */
  void CreateTasks(AliAnalysisManager* mgr)
  {
    AliAnalysisTaskSE* task = CreateTask(mgr);
    if (!task) {
      Fatal("", "Failed to create the task");
      return;
    }

#if 0
    gSystem->RedirectOutput("settings.txt");
    Printf("=== Settings for %s train === ", mc ? "mc" : "dt");
    task->Dump();
    gSystem->RedirectOutput(0);
#endif

  }
  /** 
   * Get the train setup name 
   * 
   * @return The train setup name 
   */
  const char* ClassName() const { return "MakeTrackletTrain"; }
};
//
// EOF
//
