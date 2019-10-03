#ifndef __CINT__
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisCuts.h"
//#include "AliFlowTrackSimple.h"      // added as hint for hidden library dependency to libPWGflowBase
//#include "AliFlowCandidateTrack.h"   // added as hint for hidden library dependency to libPWGflowTasks
//#include "AliCFContainer.h"          // added as hint for hidden library dependency to libCORRFW
//#include "AliAODRecoDecayHF2Prong.h" // added as hint for hidden library dependency to libPWGHFvertexingHF
#include "AliHFAssociatedTrackCuts.h"
#include "AliAnalysisTaskDxHFEParticleSelection.h"
#include "AliAnalysisTaskDxHFECorrelation.h"
#include "AliDxHFEParticleSelection.h"
#include "AliDxHFEParticleSelectionEl.h"
#include "AliDxHFEParticleSelectionD0.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliHFEcuts.h"
#include "AliLog.h"
#include "TObject.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "AliRDHFCutsD0toKpi.h"
using namespace std;
#endif

/// @file   AddTaskDxHFE.C
/// @author Matthias.Richter@ift.uib.no, Hege.Erdal@ift.uib.no
/// @date   2013-02-12
/// @brief  Add the ParticleSelection task to the manager
///
int AddTaskDxHFE()
{

  TString ofilename;
  Int_t system=0;
  TString taskOptions;
  Int_t Particle=AliAnalysisTaskDxHFEParticleSelection::kD0;


  /////////////////////////////////////////////////////////////////////////////    ofilename=
  
  gROOT->LoadMacro("AddTaskDxHFECorrelation.C");
  gROOT->LoadMacro("AddTaskDxHFECorrelationME.C");
  gROOT->LoadMacro("AddTaskDxHFEParticleSelection.C");
  //  gROOT->LoadMacro("AddSingleTrackEfficiencyTaskDxHFECorrelations.C");


  TString filename;
  TString configuration;
  if(gDirectory){
    const char* confObjectName="run_single_task_configuration";
    TObject* confObject=gDirectory->FindObject(confObjectName);
    if (confObject) {
      configuration=confObject->GetTitle();
      cout << endl <<"=======================================" << endl;
      cout << "configuration: " << configuration.Data() << endl;
      TObjArray* tokens=configuration.Tokenize(" ");
      if (tokens) {
	TIter next(tokens);
	TObject* token;
	while ((token=next())) {
	  TString argument=token->GetName();
	  if (argument.BeginsWith("file=")) {
	    filename=argument;
	  }
	  if (argument.BeginsWith("mc")) {
            filename+=" "+argument;
          }
	}	    
      }
      delete tokens;
    }
  }

  ////////////////////////////////////////////
  //          Setting up tasks              //
  ////////////////////////////////////////////


  ////////////////////////////////////////
  //      Correlation, mixed event      //
  ////////////////////////////////////////


  ///////////////////////////////////
  //        Default task           //
  //       4ITS - kFirst           //
  //   Regular correlation setup   //
  // run mode: reducedmode(default)//
  //                               //
  // fBit4   4ITSkF  maxTOF=1.0GeV //
  // system pPb      optTOF=4.0GeV //
  // Effmaps: D:no   el: yes       //
  // InvMass: 0.150  Trig:kINT7    //
  // TPC rejection above 1GeV:     // Not configurable yet
  // proton:on, nsig:-3<nsig<3     // Not configurable yet
  // pion:on, nsig: -4<nsig>4      // Not configurable yet
  //                               //
  // D0 selection: Default         //
  // D0 cutset 0                   // Not configurable yet
  ///////////////////////////////////

  taskOptions=filename+ " filterbit=4 maxPtCombinedPID=4. onlyTOFwhenpresent itsclusters=4 itsreq=kFirst extraname=4ITSkF useinvmasscut invmasscut=0.150 system=p-Pb triggermask=kINT7 TrackEffName=STE3D_DxHFECorrpPb_03May16PID_4T_1ppi_35Spi_4ITSkF.root "; // "; //fillD0scheme=both
  if(!AddTaskDxHFECorrelationME(taskOptions.Data()))
    return 0;
  

  ///////////////////////////////
  //      Electron purity      //
  ///////////////////////////////
  /*
  //Hadron
  taskOptions=filename+ " filterbit=4 maxPtCombinedPID=4. onlyTOFwhenpresent itsclusters=4 itsreq=kFirst extraname=4ITSkF useinvmasscut invmasscut=0.15 triggermask=kINT7 system=p-Pb mc elmcreco=afterfullpid particle=electron ElSelection=hadron ";
  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))
  return 0;
  
  //Standard
  taskOptions=filename+ " filterbit=4 maxPtCombinedPID=2.5 itsclusters=4 itsreq=kFirst extraname=4ITSkFnoelsel useinvmasscut invmasscut=0.15 triggermask=kINT7 system=p-Pb mc elmcreco=afterfullpid particle=electron ";
  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))
  return 0;
  */

  ///////////////////////////////////
  //       DxHFE eff default       //
  //    kINT7    PDG:electron      //
  //                               //
  ///////////////////////////////////

 /*if(!AddSingleTrackEfficiencyTaskDxHFECorrelations(kTRUE, "ElectronFbit0", AliPID::kElectron,11, AliVEvent::kINT7, kFALSE, AliCFSingleTrackEfficiencyTask::kSlow, AliSingleTrackEffCuts::kNoBayesianPID, "", ""))
     return 0;
 */

  /////////////////////////////////
  // Old examples below this line//
  /////////////////////////////////

  
  ///////////////////////////////////
  //       4ITS - kFirst           //
  //   Regular correlation setup   //
  // fBit0   4ITSkF  maxTOF=2.0GeV //
  // Effmaps: D:yes  el: yes       //
  // InvMass: 0.150  Trig:kAnyInt  //
  ///////////////////////////////////
  
  /* taskOptions=filename+ " filterbit=0 maxPtCombinedPID=2.0 itsclusters=4 itsreq=kFirst useinvmasscut invmasscut=0.150 triggermask=kAnyInt extraname=4ITSkF reducedmode system=p-Pb PromptD0EffName=D0Eff_From_c_wLimAcc_2D_pPb.root TrackEffName=STE3D_DxHFECorrpPb_14Feb16.root ";
  if(!AddTaskDxHFECorrelationME(taskOptions.Data()))
     return 0;
  */

  /* Outdated, from old effmap code
  ////////////////////////////////////////////////////////////////////////////////////
  //get nonHFE and HFE from physical primaries
  taskOptions=filename+ " filterbit=4 maxTOFpt=2 elsource=conv  extraname=HFE reducedmodewithmc";
  if(!AddTaskSingleTrackEfficiencyDxHFE(taskOptions.Data()))
     return 0;
  */
  /*
  //Get conversion electrons (need notusePhysPrim) with and without max radius of 2 cm  
  taskOptions=filename+ " filterbit=4 maxTOFpt=2 elsource=conv notusePhysPrim  extraname=ConvNoPhysPrimReq ";
  if(!AddTaskSingleTrackEfficiencyDxHFE(taskOptions.Data()))
     return 0;
*/
  // //================================================================
  // //  D0e (ME) with cuts on 120clusters on electron, trigger=D0, fill both, plus add D2H invmasstask

  // taskOptions=filename+" name=DxHFE trigger=D0 fillD0scheme=both runD0MassReference";

  // if(gDirectory) gDirectory->Clear();
  // if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  // if(!AddTaskDxHFECorrelationME(taskOptions.Data()))
  //   return 0;



  //==========================================//
  //                                          //
  //                                          //
  //           Without Inv. Mass              //
  //                                          //
  //                                          //
  //                                          //
  //==========================================//



  //================================================================
  //  D0e (ME) with cuts 
 
  ///////////////////////////
  //AddTaskDxHFECorrelation//
  ///////////////////////////
  /*  taskOptions=filename+" name=DxHFE extraname=ITS3kAnyCorr itsclusters=3 itsreq=kAny particle=electron p-Pb";
  
  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFECorrelation(taskOptions.Data()))
    return 0;
  
  taskOptions=filename+" name=DxHFE extraname=ITS4kFIRSTCorr itsclusters=4 itsreq=kFIRST particle=electron p-Pb";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFECorrelation(taskOptions.Data()))
    return 0;
  
  /////////////////////////////////
  //AddTaskDxHFEParticleSelection//
  /////////////////////////////////
  */
  /*
  taskOptions=filename+" name=DxHFE extraname=ITS3kAny itsclusters=3 itsreq=kAny particle=electron storelastcutstep system=p-Pb EMCALPID";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;
  */
  /*
  taskOptions=filename+" name=DxHFE extraname=ITS4kAny itsclusters=4 itsreq=kAny particle=electron storelastcutstep p-Pb";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;

  taskOptions=filename+" name=DxHFE extraname=ITS3kFIRST itsclusters=3 itsreq=kFIRST particle=electron storelastcutstep p-Pb";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;

  taskOptions=filename+" name=DxHFE extraname=ITS4kFIRST itsclusters=4 itsreq=kFIRST particle=electron storelastcutstep p-Pb";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;
*/

  //==========================================//
  //                                          //
  //                                          //
  //             With Inv. Mass               //
  //                                          //
  //                                          //
  //                                          //
  //==========================================//
  /*
  taskOptions=filename+" name=DxHFE trigger=D0 fillD0scheme=both extraname=DefaultIM particle=electron useinvmass";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))
    return 0;


  taskOptions=filename+" name=DxHFE trigger=D0 fillD0scheme=both extraname=ITS4kAnyIM itsclusters=4 itsreq=kAny particle=electron useinvmass";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;

  taskOptions=filename+" name=DxHFE trigger=D0 fillD0scheme=both extraname=ITS4kFIRSTIM itsclusters=4 itsreq=kFIRST particle=electron useinvmass";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
    return 0;

  taskOptions=filename+" name=DxHFE trigger=D0 fillD0scheme=both extraname=ITS4kAnyCAitsIM itsclusters=4 itsreq=kAny particle=electron elreco=afterhfeits useinvmass";

  if(gDirectory) gDirectory->Clear();
  if(gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", taskOptions.Data()));

  if(!AddTaskDxHFEParticleSelection(taskOptions.Data()))//ParticleSelection(taskOptions.Data()))
  return 0;*/
    return 1;
}

