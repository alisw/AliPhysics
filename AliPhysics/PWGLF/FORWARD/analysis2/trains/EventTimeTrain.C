// EventTimeTrain.C
#ifndef __CINT__
#include <AliAnalysisManager.h>
#include <fstream>
#else
class AliAnalysisManager;
#endif
#include "TrainSetup.C"
#include "ParUtilities.C"

/** 
 * Train to record time of each event 
 * 
 * @ingroup pwglf_forward_eventtime
 */        
class EventTimeTrain : public TrainSetup
{
public:
  /** 
   * Constructor 
   * 
   * @param name The name 
   */
  EventTimeTrain(const char* name="eventTime") : TrainSetup(name)
  { 
    fOptions.Set("type", "ESD");
  }
  /** 
   * Create our tasks 
   * 
   */
  void CreateTasks(AliAnalysisManager*)
  {
    if (!ParUtilities::MakeScriptPAR(fRailway->Mode() == Railway::kLocal,
				     "EventTimeTask.C",
				     // Gui because of CDB 
				     // XMLParser because of CDB
				     // CDB because look-up of trigger config
				     "Gui,XMLParser,"
				     "STEERBase,CDB,ESD,AOD,ANALYSIS,OADB,"
				     "ANALYSISalice",
				     fRailway)) 
      Fatal("","Failed to make PAR");
    fRailway->LoadLibrary("EventTimeTask");
    gROOT->ProcessLine("EventTimeTask::Create()");
  }
  /** 
   * Do not create a physics selection
   */
  void CreatePhysicsSelection(Bool_t, AliAnalysisManager*) {}
  /** 
   * Do not create a centrality selection
   */
  void CreateCentralitySelection(Bool_t*) {}
  /** 
   * Do not create an output handler
   */
  AliVEventHandler* CreateOutputHandler(UShort_t) { return 0; }
  /** 
   * The train class name 
   * 
   * @return Train class name
   */
  const char* ClassName() const { return "EventTimeTrain"; }
  /** 
   * Overloaded to create new dNdeta.C and dndeta.sh in the output 
   * directory
   * 
   * @param asShellScript 
   */
  void SaveSetup(Bool_t asShellScript)
  {
    TrainSetup::SaveSetup(asShellScript);
    
    SaveSort();
  }
  void SaveSort()
  {
    std::ofstream o("Sort.C");
    o << "// Written by " << ClassName() << "\n"
      << "void Sort(const char* prefix=\"\",\n"
      << "          const char* fileName=\"time.root\",\n"
      << "          const char* outName=\"map.root\",\n"
      << "          const char* treeName=\"T\")\n"
      << "{\n"
      << "  gSystem->AddIncludePath(\"-DNO_TASK -I$ALICE_PHYSICS/include\");\n"
      << "  TString mac(\"EventTimeTask/EventTimeTask.C+g\");\n"
      << "  if (prefix && prefix[0] != '\\0') mac.Prepend(prefix);\n"
      << "  gROOT->LoadMacro(mac);\n"
      << "  EventTimeSorter s;\n"
      << "  if (!s.Run(fileName,outName,treeName)) return;\n"
      << "  s.Test(fileName,outName,treeName);\n"
      << "}\n"
      << std::endl;
    o.close();
  }
  void PostShellCode(std::ostream& f)
  {
    f << "  echo \"=== Sort results ...\"\n"
      << "  aliroot -l -b -q ${prefix}Sort.C\\(\\\"${prefix}\\\"\\)\n"
      << std::endl;
  }
};
// EOF
