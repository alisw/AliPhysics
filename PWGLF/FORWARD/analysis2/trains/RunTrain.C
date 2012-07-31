/**
 * @file   RunTrain.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue May 29 09:15:22 2012
 * 
 * @brief  Wrapper script to run trains 
 * 
 * 
 * @ingroup pwglf_forward_trains
 */
class TrainSetup;

/** Global pointer to train setup */
TrainSetup* trainObj = 0;

/** 
 * Function to run a train. 
 *
 * @verbatim 
Usage: RunTrain(CLASS,NAME,OPTIONS,RUNS,NEVENTS)

CLASS:       The train-setup class to use
NAME:        Name of the train
OPTIONS:     Comma separated list of one or more of:
             help                           Show this help
             par                            Use PAR files (PROOF and Grid)
             mc                             Assume simulation input
             debug                          Execute in debugger
             type=AOD|ESD                   Type of train
             mode=LOCAL|PROOF|GRID          Execution mode
             oper=TEST|TERMINATE|FULL|INIT  Operation mode
             date=YYYY-MM-DD HH:MM:SS       Set date string
             cluster=HOST                   PROOF cluster
             dataSet=NAME                   Data set (PROOF only)
             dataDir=DIRECTORY              Data directory
             pass=NUMBER                    ESD Pass (grid only)
             verb=NUMBER                    Verbosity
             root=TAG                       ROOT version (Grid)
             aliroot=TAG                    AliROOT version (Grid)
             alien=TAG                      AliEn API version (Grid)
             overwrite                      Allow overwrite
RUNS:        Comma separated list of run numbers, or file names of
             files that contain run numbers
NEVENTS:     Number of events to analyse
@endverbatim
 * 
 * Individual classes derived from TrainSetup may define additional
 * (or less) options.  Pass the option @b help to see the list of
 * available options
 *
 * Usually, the train defines the execution mode (ESD or AOD) and the
 * option @b type should not be given
 *
 * If a date is set (option @b date) either to a string of the form 
 * @c YYYY-MM-DD HH:MM:SS or to the string @c now, then the date will be 
 * appended to the job name. 
 * 
 * A sub-directory named according to the name passed will be made,
 * and relevant files will copied there.
 *
 * @param trainClass Class name of the train setup
 * @param trainName  Name of the train 
 * @param options    Options for the train 
 * @param runs       Runs to add to the train 
 * @param nEvents    Number of events to analyse 
 *
 * @ingroup pwglf_forward_trains_run
 */
void RunTrain(const char* trainClass, 
	      const char* trainName, 
	      const char* options="", 
	      const char* runs="", 
	      Int_t       nEvents=-1)
{
  const char* builder = 
    "$(ALICE_ROOT)/PWGLF/FORWARD/analysis2/trains/BuildTrain.C";
  gROOT->LoadMacro(builder);
  
  BuildTrain(trainClass);

  gROOT->ProcessLine(Form("trainObj = new %s(\"%s\")", 
			  trainClass, trainName));
  if (!trainObj) {
    Error("RunTrain", "Failed to make train %s of class %s", 
	  trainName, trainClass);
    gApplication->Terminate();
  }
  TrainSetup::Runner r(*trainObj);
  if (!r.Init(options)) return;
  if (r.IsHelpAsked()) { 
    std::ostream& o = std::cout;
    o << "Usage: RunTrain(CLASS,NAME,OPTIONS,RUNS,NEVENTS)\n\n"
      << "CLASS:       The train-setup class to use\n"
      << "NAME:        Name of the train\n"
      << "OPTIONS:     Comma separated list of one or more of:\n";
    r.PrintHelp(o, "           ");
    o << "RUNS:        Comma separated list of run numbers, or file names of\n"
      << "             files that contain run numbers\n"
      << "NEVENTS:     Number of events to analyse\n"
      << std::endl;
    return;
  }
  r.Run(runs, nEvents);
}
/*
 * EOF
 */
