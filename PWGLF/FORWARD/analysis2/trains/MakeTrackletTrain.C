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
    fOptions.Add("eta-min",       "ETA",   "Least eta",                 -1.8);
    fOptions.Add("eta-max",       "ETA",   "Largest eta",               +1.8);
    fOptions.Add("ipz-min",       "CM",    "Least IPz",                 -10);
    fOptions.Add("ipz-max",       "CM",    "Largest IPz",               +10);
    fOptions.Add("mix-t-min",     "N",     "Min Trackets to mix",       1);
    fOptions.Add("mix-t-max",     "N",     "Max Trackets to mix",       20000);
    fOptions.Add("mix-t-bins",    "N",     "Tracklet mixing bins",      20000);
    fOptions.Add("mix-ipz-bins",  "N",     "IPz mixing bins",           20);
    fOptions.Add("cent",          "METHOD","Centrality selector",       "V0M");
    fOptions.Add("cent-bins",     "Bins",  "Centrality bins",
		 "0-5-10-20-30-40-50-60-70-80");
    fOptions.Add("sig-n-std",     "X",     "Cut on weighted distance",  1.5);
    fOptions.Add("sig-dphi-s",    "X",     "Cut on dPhi-phiBent",       -1);
    fOptions.Add("mc-v0-scale",   "X",     "Scaling of MC V0",          0.7520);
    fOptions.Add("fill-reco",              "Fill with new reco",        true);
    fOptions.Add("create-inj",             "Create injection BG",       true);
    fOptions.Add("create-rot",             "Create rotation BG",        false);
    fOptions.Add("create-mix",             "Create mixed BG",           false);
    fOptions.Add("phi-rot",       "X",     "Rotation BG angle",   TMath::Pi());
    fOptions.Add("inj-scale",     "X",     "Injection BG scale",        1.);
    fOptions.Add("scale-dtheta",           "Scale dTheta",              true);
    fOptions.Add("dphi",          "X",     "dPhi window",               0.06);
    fOptions.Add("dtheta",        "X",     "dTheta window",             0.025);
    fOptions.Add("phi-shift",     "X",     "Bending shift",             0.0045);
    fOptions.Add("remove-outliers",        "Whether to remove outliers",true);
    fOptions.Add("outlier-phi",   "X",     "Phi outlier cut",           0.005);
    fOptions.Add("outlier-z-eta", "X",     "Z-Eta outlier cut",         0.05);
    fOptions.Add("check-reco",             "Check reconstruction",      false);
    
  }
  /** 
   * Create the input handler.  This is overwritten from the base
   * class to allow using AliESDInputHandlerRP for rec. points., and
   * AliMixInputEventHandler if requested.
   * 
   * @param type Type of analysis 
   * 
   * @return The input handler 
   */
  AliVEventHandler* CreateInputHandler(UShort_t type)
  {
    
    Bool_t   fill_reco       = fOptions.AsBool("fill-reco");;
    Bool_t   create_inj      = fOptions.AsBool("create-inj");;
    Bool_t   create_rot      = fOptions.AsBool("create-rot");;
    Bool_t   create_mix      = fOptions.AsBool("create-mix");;
    bool     needRec = fill_reco || create_inj || create_rot ||  create_mix;

    AliVEventHandler* ret = 0;
    AliESDInputHandlerRP* esd = 0;
    if (type == Railway::kESD && needRec) 
      ret = esd = new AliESDInputHandlerRP();
    else
      ret = TrainSetup::CreateInputHandler(type);

    // If mixing requtest, set it up
    if (create_mix && esd) {
      Int_t    nMix     = 1;
      // Int_t    ntMin    = fOptions.AsInt("mix-t-min", 1);
      // Int_t    ntMax    = fOptions.AsInt("mix-t-max", 20000);
      // Int_t    ntBin    = fOptions.AsInt("mix-t-bins", 20000);
      Int_t    ipzBin   = fOptions.AsInt("mix-ipz-bins", 14);
      Double_t ipz_min  = fOptions.AsDouble("ipz-min",-7);;
      Double_t ipz_max  = fOptions.AsDouble("ipz-max",+7);;

      Info("", "Creating mix");
      // Execute this in a macro;
      TMacro m;
      m.AddLine("{");
      m.AddLine(Form("AliMixInputEventHandler* mixer = "
		     "new AliMixInputEventHandler(%d);", nMix));
      m.AddLine(Form("mixer->SetInputHandlerForMixing("
		     "(const AliInputEventHandler *const)%p);",esd));
      m.AddLine("AliMixEventPool* pool = new AliMixEventPool(\"pool\");");
      m.AddLine(Form("pool->AddCut(new AliMixEventCutObj("
		     "AliMixEventCutObj::kZVertex,%f,%f,%d));",
		     ipz_min,ipz_max, ipzBin));
      m.AddLine("mixer->SetEventPool(evPool)");
      m.AddLine(Form("((AliVEventHandler)%p)->SetMixingHandler(mixer);",
		     esd));
      m.AddLine("}");
      m.Exec();
    }

    Info("","Input handler: %p", ret);
    return ret;    
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
   * Create our task, and return it.  This uses the interpreter to
   * make the object.  
   * 
   * @param mgr Analysis manager 
   * 
   * @return Pointer to the task 
   */
  AliAnalysisTaskSE* CreateTask(AliAnalysisManager* mgr)
  {
    const char*         cls ="AliTrackletTaskMulti";
    Info("", "Creating task %p", gROOT->GetClass(cls));
    Long_t              ret =gROOT->ProcessLine(Form("new %s(\"%s\")",cls,cls));
    AliAnalysisTaskSE* task =reinterpret_cast<AliAnalysisTaskSE*>(ret);
    mgr->AddTask(task);

    AliAnalysisDataContainer *out =
      mgr->CreateContainer("clist", TList::Class(),
			   AliAnalysisManager::kOutputContainer,
			   "trbg.root");
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
    // fRailway->LoadLibrary("RubensCode",true,true);
    Double_t sig_dphi_s      = fOptions.AsDouble("sig-dphi-s",-1);
    Double_t sig_n_std       = fOptions.AsDouble("sig-n-std",1.5);
    Double_t dphi            = fOptions.AsDouble("dphi",0.06);

    if (sig_dphi_s<0)
      fOptions.Set("sig-dphi-s", TMath::Sqrt(sig_n_std)*dphi);

    // gROOT->SetMacroPath(Form("%s:$ALICE_PHYSICS/PWGUD/multVScentPbPb",
    //  gROOT->GetMacroPath()));
    //  gSystem->AddIncludePath("-I$ALICE_PHYSICS/PWGUD/multVScentPbPb");
    Info("CreateTasks", "Loading code");
    fRailway->LoadAux("AliITSMultRecBg.h");
    fRailway->LoadAux("AliTrackletTaskMulti.h");
    fRailway->LoadSource("AliITSMultRecBg.cxx");
    fRailway->LoadSource("AliTrackletTaskMulti.cxx");

    Bool_t             mc  = mgr->GetMCtruthEventHandler() != 0;
    AliAnalysisTaskSE* task = CreateTask(mgr);
    if (!task) {
      Fatal("", "Failed to create the task");
      return;
    }
    SetOnTask(task, "UseMC", mc);
    FromOption(task, "UseCentralityVar", 	"cent", 	"");
    FromOption(task, "EtaMin",			"eta-min",	-.5);
    FromOption(task, "EtaMax",			"eta-max",	+.5);
    FromOption(task, "ZVertexMin",		"ipz-min",	-7.);
    FromOption(task, "ZVertexMax",		"ipz-max",	+7.);
    FromOption(task, "NStdCut",			"sig-n-std",	1.5);
    FromOption(task, "DPhiSCut",		"sig-dphi-s",	-1.);
    FromOption(task, "ScaleMCV0",		"mc-v0-scale",	0.7520);
    FromOption(task, "DoNormalReco",		"fill-reco",	false);
    FromOption(task, "DoInjection",		"create-inj",	false);
    FromOption(task, "DoRotation",		"create-rot",	false);
    // FromOption(task, "DoMixing",		"create-mix",	false);
    FromOption(task, "PhiRot",			"phi-rot",	TMath::Pi());
    FromOption(task, "InjScale",		"inj-scale",	1.);
    FromOption(task, "ScaleDThetaBySin2T",	"scale-dtheta",	false);
    FromOption(task, "PhiWindow",		"dphi",		0.06);
    FromOption(task, "ThetaWindow",		"dtheta",	0.025);
    FromOption(task, "PhiShift",		"phi-shift",	0.0045);
    FromOption(task, "RemoveOverlaps",		"remove-outliers",false);
    FromOption(task, "PhiOverlapCut",		"outlier-phi",	0.005);
    FromOption(task, "ZetaOverlapCut",		"outlier-z-eta",0.05);
    FromOption(task, "SetCheckReconstructables","check-reco",	false);

    TString centBins = fOptions.AsString("cent-bins");
    TObjArray* tokens = centBins.Tokenize("-");
    TArrayD array(tokens->GetEntries());
    for (Int_t i = 0; i < array.GetSize(); i++) {
      TObjString* ostr = static_cast<TObjString*>(tokens->At(i));
      TString&    str  = ostr->String();
      array[i]         = str.Atof();
    }
    const char* cls = "AliTrackletTaskMulti";
    gROOT->ProcessLine(Form("((%s)%p)->SetCentPercentiles((Double_t*)%p,%d)",
			    cls, task, array.GetArray(), array.GetSize()-1));
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
