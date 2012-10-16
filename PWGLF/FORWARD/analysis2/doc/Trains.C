#error This file is not for compilation 
/** 
 * @page train_setup_doc Using the TrainSetup facility 
 *
 * @tableofcontents
 *
 * @section train_setup_overview Overview 
 *
 * Users should define a class that derives from TrainSetup.  The class
 * should implement the member function TrainSetup::CreateTasks to add
 * needed tasks to the train
 * 
 * @code 
 * // MyTrain.C 
 * class MyTrain : public TrainSetup
 * {
 * public:
 *   MyTrain(const char* name="MyTrain")
 *     : TrainSetup(name), 
 *       fParameter(0)
 *   {
 *      // SetType(kAOD); // AOD input 
 *      // SetType(kESD); // ESD input 
 *      // Input unspecified - can be set later 
 *   }
 * protected:
 *   void CreateTasks(EMode mode, Bool_t par, AliAnalysisManager* mgr)
 *   {
 *     AliAnalysisManager::SetCommonFileName("my_analysis.root");
 *     LoadLibrary("MyAnalysis", par, true);
 *     Bool_t mc = mgr->GetMCtruthEventHandler() != 0;
 *     gROOT->Macro(Form("AddTaskMyAnalysis.C(%f)",fParameter));
 *   }
 *   const char* ClassName() const { return "MyTrain"; }
 *   void MakeOptions(Runner& r)
 *   { 
 *     TrainSetup::MakeOptions(r);
 *     r.Add(new Option("parameter", "My parameter", "VALUE"));
 *   }
 *   void SetOptions(Runner& r)
 *   {
 *     TrainSetup::SetOptions(r);
 *     Option* param = r.FindOption("parameter");
 *     if (param)  fParameter = param->AsDouble();
 *   }
 *   Double_t fParameter;
 * };
 * @endcode 
 *
 * @section train_setup_exec Execution of the train 
 *
 * A user defined TrainSetup class can then be run like 
 * 
 * @code 
 * Root> .x RunTrain.C("<class>", "<name>", "<options>", "<runs>", nEvents)
 * @endcode 
 * 
 * or using the program @b runTrain
 * 
 * @code
 * > runTrain --class=<class> --name=<name> [<options>] 
 * @endcode
 *
 * Here, 
 *
 * - &lt;class&gt; is the name of the user defined class deriving from
 *   TrainSetup.
 * 
 * - &lt;name&gt; is an arbitary name to give to the train. Note, an
 *   @e escaped @e name will be generated from this, which replaces
 *   all spaces and the like with '_' and (optionally) with the date
 *   and time appended. 
 *
 * - &lt;options&gt; is a list of options.  For RunTrain this is a
 *   comma separated list of options in the form
 *   &lt;option&gt;=&lt;value&gt; for value options and &lt;option&gt;
 *   for flags (booleans).  For @c runTrain, the options are of the
 *   traditional Unix long type: @c --&lt;option&gt;=&lt;value&gt; and
 *   @c --&lt;option&gt;.  The exact list of options for a given train
 *   can be listed by passing the option @b help.
 *
 * See also ::RunTrain and ::main
 *
 * In both cases, a new sub-directory called @e escaped @e name of the
 * train is created, and various files are copied there - depending on
 * the mode of execution.
 *
 * For local analysis, no aditional files are copied there, but the
 * output will be put there. 
 *
 * For PROOF analysis, the needed PAR files are copied there and
 * expanded.  The output of the job may end up in this directory if so
 * instructed.
 *
 * For Grid analysis, various JDL and steering scripts are copied to
 * this directory.
 *
 * In all cases, a file named @c rerun.C (and for @b runTrain:
 * rerun.sh) is generated in this sub-directory.  It contains the
 * setting used for the train and can easily be used to run merging
 * and terminate as needed.
 *
 * @section train_setup_proof_spec PROOF specifics
 *
 * Local and Grid jobs are in a sense very similar.  That is, the
 * individual Grid jobs are very much like Local jobs, in that they
 * always produce output files which albiet not after Terminate.
 *
 * PROOF jobs are very different.  In a PROOF analysis, each slave
 * only produces in memory output which is then sent via net
 * connections (sockets) to the master.  One therefore needs to be
 * very of output object ownership and the like.  
 *
 * Another major difference is that output files are generated within
 * the PROOF cluster, and are generally not accessible from the
 * outside.  For plain PROOF clusters in a local area network or
 * so-called <i>Lite</i> session, it is generally not a problem since
 * the files are accessible on the LAN or local machine for Lite
 * sessions.  However, for large scale analysis farms (AAFs), the
 * workers and masters are generally on a in-accessible sub-net, and
 * there's no direct access to the produced files.  Now, for normal
 * output files, like histogram files, etc. there are provisions for
 * this, which means the final merged output is sent back to the
 * client.  Special output, such as AODs, are however not merged nor
 * sent back to the user by default.  There are two ways to deal with this: 
 * 
 * <ol>
 *  <li> Register the output tree as a data set on the cluster.  This is useful if you need to process the results again on the cluster.</li>
 *  <li> Send the output to a (possibly custom) XRootd server.  This is useful if you need to process the output outside of the cluster</li> 
 * </ol>
 *
 * The first mode is specified by passing the option
 * <tt>dsname=</tt><i>&lt;name&gt;</i> in the cluster URI.  The created
 * dataset will normally be made in
 * <tt>/default/</tt><i>&lt;user&gt;</i><tt>/</tt><i>&lt;name&gt;</i>. If the
 * <tt>=</tt><i>&lt;name&gt;</i> part is left out, the <i>escaped name</i> of
 * the job will be used.  
 *
 * The second mode is triggered by passing the option
 * <tt>storage=<i>URI</i></tt> to the train setup.  The <i>URI</i>
 * should be of the form
 *
 * @code
 *   rootd://<host>[:<port>]/<path>
 * @endcode
 * 
 * where <i>&lt;host&gt;</i> is the name of a machine accessible by
 * the cluster, <i>&lt;port&gt;</i> is an optional port number (e.g.,
 * if different from 1093), and <i>&lt;path&gt;</i> is an absolute
 * path on <i>&lt;host&gt;</i>.
 * 
 * The XRootd process should be started (optionally by the user) on
 * <i>&lt;host&gt;</i> as
 *
 * @code
 *    xrootd -p <port> <path>
 * @endcode
 *
 * When running jobs on AAFs, one can use the Grid handler to set-up
 * aspects of the job.  However, sometimes it's desirable to leave the
 * Grid handler out.  To do that, pass the option <tt>plain</tt> in
 * the cluster URI.
 *
 * @section train_setup_input Specifying the input
 * @subsection train_setup_local Local data input
 * 
 * For both ESD and AOD input for local jobs, one must specify the
 * root of the sub-tree that holds the data.  That is, if - for
 * example - the data resides in a directory structure like
 *
 * <pre>
 *   /some/directory/&lt;run&gt;/&lt;seq&gt;/AliESDs.root 
 * </pre>
 *
 * then one should specify the input location like 
 *
 * @code 
 *   train->SetDataDir("/some/directory");
 * @endcode
 *
 * <tt>/some/directory</tt> is then search recursively for input files
 * that match the pattern given by the analysis type (ESD:
 * <tt>AliESDs.root</tt>, AOD: <tt>AliAOD.root</tt>).  The found files
 * are then chained together.  If MC input is specified, then the
 * companion files <tt>galice.root</tt>, <tt>Kinematics.root</tt>, and
 * <tt>TrackRefs.root</tt> must be found in the same directories as
 * the <tt>AliESDs.root</tt> files
 *
 * @subsection train_setup_proof PROOF input. 
 * 
 * The input data for a PROOF based analysis can be specified as per a
 * Local job if the cluster used is local, in which case the data must
 * be available to the slaves at the specified locations, or one can
 * specify a data-set name via
 * 
 * @code 
 *   train->SetDataSet("<data-set-name>");
 * @endcode 
 *
 * @b Note: For AAFs using the Grid Handler one <i>must</i> use data sets. 
 *
 * @subsection train_setup_grid_esd Grid ESD input. 
 *
 * Suppose the ESD files are stored on the Grid as 
 *
 * <pre>
 *   /alice/data/&lt;year&gt;/&lt;period&gt;/&lt;run&gt;/ESDs/pass&lt;no&gt;/&lt;year&gt;&lt;run&gt;&lt;chunk&gt;.&lt;part&gt;/AliESDs.root 
 * </pre> 
 *
 * where &lt;run&gt; is zero-padded by typically 3 '0's.  One should
 * specify the input location like
 *
 * @code 
 *   train->SetDataDir("/alice/data/<year>/<period>");
 *   train->SetDataPattern("ESDs/pass<no>/&ast;/");
 *   train->AddRun(<run>);
 * @endcode
 *
 * If a particular kind of pass is needed, say
 * <tt>pass&lt;no&gt;_MUON</tt>, one should do
 *
 * @code 
 *  train->SetDataPattern("ESDs/pass<no>_MUON/&ast;/");
 * @endcode
 * 
 * The AliEn analysis plug-in is then instructed to look for data files under 
 * 
 * <pre>
 *   /alice/data/&lt;year&gt;/&lt;period&gt;/&lt;run&gt;/ESDs/pass&lt;no&gt;/&nbsp;*&nbsp;/AliESDs.root 
 * </pre>
 *
 * for each added run. 
 *
 * For simulation output, the files are generally stored like 
 * 
 * <pre>
 *   /alice/sim/&lt;year&gt;/&lt;prod&gt;/&lt;run&gt;/&lt;seq&gt;/AliESDs.root 
 * </pre> 
 *
 * where &lt;run&gt; is generally @e not zero-padded. One should
 * specify the input location like
 *
 * @code 
 *   train->SetDataDir("/alice/data/<year>/<period>");
 *   train->SetDataPattern("*");
 *   train->AddRun(<run>);
 * @endcode
 *
 *
 * @subsection train_setup_grid_aod Grid AOD input 
 * 
 * Suppose your AOD files are placed in directories like 
 * 
 * <pre>
 *   /some/directory/&lt;run&gt;/&lt;seq&gt;/AliAOD.root
 * </pre>
 *
 * where &lt;run&gt; is zero-padded by typically 3 '0's.  One should
 * then specify the input as 
 *
 * @code 
 *   train->SetDataDir("/some/directory");
 *   train->SetDataPattern("*");
 *   train->AddRun(<run>);
 * @endcode
 * 
 * The AliEn analysis plug-in is then instructed to look for data files under 
 * 
 * <pre>
 *   /some/directory/&lt;run&gt;/&nbsp;*&nbsp;/AliAOD.root 
 * </pre>
 *
 * for each added run.  
 *
 * Suppose the AODs are in 
 *
 * <pre>
 *  /alice/data/&lt;year&gt;/&lt;period&gt;/&lt;run&gt;/ESDs/pass&lt;no&gt;/AOD&vers&gt;/&lt;seq&gt;/AliAOD.root 
 * </pre>
 * 
 *  @code 
 *   train->SetDataDir("/alice/data/<year>/<period>");
 *   train->SetDataPattern("ESDs/pass<no>/AOD<vers>/&ast;/");
 *   train->AddRun(<run>);
 * @endcode
 *
 * For simulation output, the files are generally stored like 
 * 
 * <pre>
 *   /alice/sim/&lt;year&gt;/&lt;prod&gt;/&lt;run&gt;/&lt;seq&gt;/AliAOD.root 
 * </pre> 
 *
 * where &lt;run&gt; is generally @e not zero-padded. One should
 * should specify the input location like 
 *
 * @code 
 *   train->SetDataDir("/alice/data/<year>/<period>");
 *   train->SetDataPattern("*");
 *   train->AddRun(<run>);
 * @endcode
 *
 * @section train_setup_other Other features  
 * @subsection train_setup_options Options interface 
 *
 * If the train does not depend on additional options or parameters,
 * the member functions TrainSetup::MakeOptions and
 * TrainSetup::SetOptions can be left un-overloaded in the derived
 * class.  However, options defined in this way can be set through the
 * command line of the program @b runTrain, and provides a great deal
 * of flexiblity.  The Option class provides means of translating the
 * passed string values to integers, doubles, booleans, and of course
 * strings.
 *
 * @subsection train_setup_aux Auxillary libraries, sources, and files
 *
 * Auxillary libraries should be loaded using 
 *
 * - TrainSetup::LoadLibrary(const char*,Bool_t,Bool_t)
 *
 * where first argument is the name of the library, the econd should
 * be true if the library should be loaded as a PAR file (PROOF and
 * Grid only), and the argument should be true if the library should
 * be loaded on the PROOF slaves/Grid workers too.
 *
 * If the train needs additional files, say a script for setting up
 * the tasks, or some data file, it can be passed on the the
 * PROOF/Grid workers using the member functions 
 *
 * - TrainSetup::AddExtraFile(const char*)
 * - TrainSetup::AddSource(const char*,bool)
 * 
 * @subsection train_setup_overload Overloading the behaviour 
 *
 * The base class TrainSetup tries to implement a sensible setup for a
 * given type of analysis, but some times a particular train needs a
 * bit of tweaking.  One can therefore overload the following functions 
 * 
 * - TrainSetup::CreateGridHandler()
 * - TrainSetup::CreateInputHandler(EType)
 * - TrainSetup::CreateMCHandler(EType,bool)
 * - TrainSetup::CreateOutputHandler(EType)
 * - TrainSetup::CreatePhysicsSelection(Bool_t,AliAnalysisManager*)
 * - TrainSetup::CreateCentralitySelection(Bool_t,AliAnalysisManager*)
 *
 * @section train_setup_scripts Tasks defined in scripts 
 *
 * A task can even be defined in a script, like for example a task like 
 * 
 * @include MyAnalysis.C 
 *
 * Our train set-up can then use the member function
 * TrainSetup::MakeScriptPAR to make a PAR file of the script and use
 * that to make a library loaded on the workers and then generate an
 * object of our task defined in the script.
 *
 * @include MyTrain.C
 *
 * This can allow for fast development and testing of analysis tasks
 * without having to wait for official tasks and builds of all of
 * AliROOT
 */
//
// EOF
//
