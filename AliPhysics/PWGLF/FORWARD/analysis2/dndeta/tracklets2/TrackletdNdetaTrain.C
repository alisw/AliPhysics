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
struct TrackletdNdetaTrain : public TrainSetup
{
  /** 
   * Constructor.  This sets up the available options 
   * 
   * @param name Free form name 
   */
  TrackletdNdetaTrain(const char* name)
    : TrainSetup(name)
  {
    // Define all our options here
    fOptions.Add("trig",          "NAME",  "Trigger to use",           "V0AND");
    fOptions.Add("cent",          "METHOD","Centrality selector",      "V0M");
    fOptions.Add("cent-bins",     "BINS",  "Centrality bins",
		 "0-5-10-20-30-40-50-60-70-80");
    fOptions.Add("max-delta",        "X",    "Cut on weighted distance",25.);
    fOptions.Add("tail-delta",       "X",    "Tail cut on distance", 5.);
    fOptions.Add("scale-dtheta",             "Scale dTheta" ,        true);
    fOptions.Add("cut-dtheta",               "Cut dTheta",           false);
    fOptions.Add("dphi-window",      "X",    "dPhi window",          0.06);
    fOptions.Add("dtheta-window",    "X",    "dTheta window",        0.025);
    fOptions.Add("dphi-shift",       "X",    "Bending shift",        0.0045);
    fOptions.Add("phi-overlap-cut",  "X",    "Phi overlap cut",      0.005);
    fOptions.Add("z-eta-overlap-cut","X",    "Z-Eta overlap cut",    0.05);
    fOptions.Add("phi-rotation",     "X",    "Rotation BG angle",TMath::Pi());
    fOptions.Add("inj-scale",        "X",    "Injection BG scale",   1.);
    fOptions.Add("shifted-dphi-cut", "X",    "Cut on dPhi-phiBent",  -1.);
    fOptions.Add("delta-cut",        "X",    "Cut on weighted distance",  1.5);
    fOptions.Add("eta-bins",         "BINS", "Eta bins",          "r16:-2:+2");
    fOptions.Add("ipz-bins",         "BINS", "IPz bins",          "u15");
    fOptions.Add("reconstruct",      "WHAT", "Reconstruction types",
		 "NORMAL,INJECTION");
    fOptions.Add("reweight",      "OPTIONS","','-seperates list",       "");
    fOptions.SetDescription("Run tracklet analysis on ESD data files to "
			    "extract dNch/deta. The train should be run on "
			    "both real data and simulated data.\n"
			    "BINS are specified as comma or colon separated "
			    "lists of bins.  If prefixed by 'r', then the "
			    "specification assumes the first number is the "
			    "number of bins, and the following one or two "
			    "numbers give the min/max.  If prefixed by 'u' "
			    "then unit bins are assumed, and the following "
			    "one or two numbers give the min/max.\n"
			    "When processing simulated data, it is possible to "
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
   * @param needRec if true, also get rec-points 
   * 
   * @return The input handler 
   */
  AliVEventHandler* CreateInputHandler(UShort_t type, Bool_t needRec=false)
  {
    TString  recMode = fOptions.AsString("reconstruct");
    recMode.ToLower();
    Bool_t needRP = (recMode.Contains("nor") ||
		     recMode.Contains("inj") ||
		     recMode.Contains("rot") || needRec);
    // Info("CreateInputHandler", "Will create for %d (%d)", type, needRP);
    return TrainSetup::CreateInputHandler(type, needRP);
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
    fRailway->LoadSource("AliTrackletdNdetaUtils.C");
    fRailway->LoadSource("AliTrackletdNdetaTask.C");

    // --- Create the task using interpreter -------------------------
    Bool_t             mc  = mgr->GetMCtruthEventHandler() != 0;
    const char*        nme = (mc ? "MidMCdNdeta" : "MiddNdeta");
    const char*        cls = (mc ?
			      "AliTrackletdNdetaMCTask" :
			      "AliTrackletdNdetaTask");
    Long_t             ret = gROOT->ProcessLine(Form("new %s(\"%s\")",cls,nme));
    AliAnalysisTaskSE* task =reinterpret_cast<AliAnalysisTaskSE*>(ret);
    if (!task) return 0;

    // --- Add task to train -----------------------------------------
    gROOT->ProcessLine(Form("((%s*)%p)->Connect()", cls, task));

    // --- Figure out the trigger options ----------------------------
    TString trg = fOptions.Get("trig");
    trg.ToUpper();
    UInt_t  sel = AliVEvent::kINT7;
    if      (trg.EqualTo("MB"))    sel = AliVEvent::kMB;
    else if (trg.EqualTo("V0AND")) sel = AliVEvent::kINT7;
    else if (trg.EqualTo("V0OR"))  sel = AliVEvent::kCINT5;
    else if (trg.EqualTo("ANY"))   sel = AliVEvent::kAny;

    
    // --- Set various options on task -------------------------------
    const char* defCent = "0-5-10-20-30-40-50-60-70-80-90";
    task->SelectCollisionCandidates(sel);
    FromOption(task, "ReconstructionMode",      "reconstruct",      "nor,inj");
    FromOption(task, "CentralityMethod", 	"cent", 	    "V0M");
    FromOption(task, "CentralityAxis",          "cent-bins",        defCent);
    FromOption(task, "EtaAxis",                 "eta-bins",         "r16:2");
    FromOption(task, "IPzAxis",                 "ipz-bins",         "u15");
    // FromOption(task, "CheckReconstructables",   "check-reco",    false);
    FromOption(task, "MaxDelta",		"max-delta",	    25.);
    FromOption(task, "TailDelta",		"tail-delta",	    5.);
    FromOption(task, "ScaleDTheta",	        "scale-dtheta",	    false);
    FromOption(task, "DPhiWindow",		"dphi-window",	    0.06);
    FromOption(task, "DThetaWindow",		"dtheta-window",    0.025);
    FromOption(task, "DPhiShift",		"dphi-shift",	    0.0045);
    FromOption(task, "PhiOverlapCut",		"phi-overlap-cut"  ,0.005);
    FromOption(task, "ZEtaOverlapCut",		"z-eta-overlap-cut",0.05);
    FromOption(task, "PhiRotation",		"phi-rotation",	TMath::Pi());
    // FromOption(task, "InjScale",		"inj-scale",	1.);
    FromOption(task, "ShiftedDPhiCut",		"shifted-dphi-cut",-1.);
    FromOption(task, "DeltaCut",		"delta-cut",	    1.5);
    if (fOptions.AsBool("cut-dtheta")) 
      FromOption(task, "ScaleDThetaCut",	"dtheta-window",    0.025);
    
    SetupReweighting(task);

    task->Print("");
    
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
  const char* ClassName() const { return "TrackletdNdetaTrain"; }
};
//
// EOF
//
