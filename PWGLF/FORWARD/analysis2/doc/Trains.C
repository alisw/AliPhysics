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
 *   {
 *      // fOptions.Set("type", "AOD"); // AOD input
 *      // fOptions.Set("type", "ESD"); // ESD input
 *      fOptions.Add("parameter", "VALUE", "Help on parameter", "value");
 *   }
 * protected:
 *   void CreateTasks(AliAnalysisManager* mgr)
 *   {
 *     AliAnalysisManager::SetCommonFileName("my_analysis.root");
 *     LoadLibrary("MyAnalysis", par, true);
 *     Bool_t mc = mgr->GetMCtruthEventHandler() != 0;
 *     Double_t param = fOptions.AsDouble("parameter");
 *     gROOT->Macro(Form("AddTaskMyAnalysis.C(%f)",param));
 *   }
 *   const char* ClassName() const { return "MyTrain"; }
 * };
 * @endcode 
 *
 * @section train_setup_params Paramters of the setup 
 *
 * Parameters of the user defined class deriving from TrainSetup is
 * best handled by adding options to the internal member @c fOptions
 * in the constructor e.g.,
 *
 * @code 
 *   fOptions.Add("<name>", "<dummy>", "<description>", "<default>");
 *   fOptions.Add("<name>", "<description>");
 * @endcode
 * 
 * The first form defined a parmater that has a value, while the
 * second form defines a flag (or toggle).  The values or flags can be
 * retrieved later by doing
 * 
 * @code 
 *    Double_t value = fOptions.AsDouble("<name>",<value if not set>);
 *    Int_t    value = fOptions.AsInt("<name>",<value if not set>);
 *    Long64_t value = fOptions.AsLong("<name>",<value if not set>);
 *    TString  value = fOptions.Get("<name>");
 *    Bool_t   value = fOptions.Has("<name>");
 * @endcode
 *
 * Parameters defined this way are directly accessible as options to
 * pass to either runTrain or RunTrain.C
 *
 * @section train_setup_exec Execution of the train 
 *
 * A user defined TrainSetup class can then be run like 
 * 
 * @code 
 * Root> .x RunTrain.C("<class>", "<name>", "<uri>", "<options>")
 * @endcode 
 * 
 * or using the program @b runTrain
 * 
 * @code
 * > runTrain --class=<class> --name=<name> --url=<uri> [<options>] 
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
 * - &lt;uri&gt; is the job execution URI which specified both the
 *   execution environment and the input data, as well as some options.  
 *   See more below. 
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
 * In all cases, a file named @c ReRun.C (and for @b runTrain:
 * rerun.sh) is generated in this sub-directory.  It contains the
 * setting used for the train and can easily be used to run merging
 * and terminate as needed.
 *
 * @section train_setup_url_spec Execution URI 
 *
 * This URI has the form 
 * @code 
 *  <protocol>://[[<user>@]<host>]/<input>[?<options>][#<treename>]
 * @endcode 
 * and specifies several things.
 * 
 * <dl>
 *   <dt><tt>&lt;protocol&gt;</tt></dt>
 *   <dd>One of 
 *     <dl>
 *       <dt><tt>local</tt></dt>
 *       <dd>Local analysis on local data executed sequentially on the
 *         local machine</dd>
 *       <dt><tt>lite</tt></dt>
 *       <dd>Proof-Lite analysis on local data executed in parallel on
 *         the local machine</dd>
 *       <dt><tt>proof</tt></dt> 
 *       <dd>Proof analysis on cluster data executed in parallel on a
 *         PROOF cluster</dd>
 *       <dt><tt>alien</tt></dt>
 *       <dd>Grid analysis on grid data executed on the Grid</dd>
 *     </dl>
 *   </dd>
 *   <dt><tt>[[&lt;user&gt;@]&lt;host&gt;]</tt></dt>
 *   <dd>Sets the master host for Proof analysis</dd>
 *   <dt><tt>&lt;input&gt;</tt></dt>
 *   <dd>Input data specification.  The exact form depends on the
 *     protocol used e.g., for local analysis it can be a single,
 *     while for other environments it could be a data set name, and
 *     so on.</dd>
 *   <dt><tt>&lt;options&gt;</tt></dt>
 *   <dd>Protocol specific options</dd>
 *   <dt><tt>&lt;treename&gt;</tt></dt>
 *   <dd>If specified, gives what data to analyse</dd>
 *  </dl>
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
 *  <li> Register the output tree as a data set on the cluster.  This
 *    is useful if you need to process the results again on the
 *    cluster.</li>
 *  <li> Send the output to a (possibly custom) XRootd server.  This
 *    is useful if you need to process the output outside of the
 *    cluster</li>
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
 * aspects of the job.  To enable the Grid handler, pass the option
 * <tt>plugin</tt> in the execution URI
 *
 * @section train_setup_input Specifying the input
 * @subsection train_setup_local Local and Lite data input
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
 *   local:///some/directory[?pattern=AliESDs.root][#esdTree]
 *   lite:///some/directory[?pattern=AliESDs.root][#esdTree]
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
 * The input data for a PROOF based analysis is specified as data set
 * names,
 * 
 * @code 
 *   proof://[<user>@]<host>/<data-set-name>[?options][#<treename>]
 * @endcode 
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
 *   alien:///alice/data/<year>/<period>?pattern=ESDs/pass<no>/&ast;/&run=<run>[#<treename>]
 * @endcode
 *
 * If a particular kind of pass is needed, say
 * <tt>pass&lt;no&gt;_MUON</tt>, one should do modify the
 * <tt>pattern</tt> option accordingly
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
 *   alien:///alice/data/<year>/<period>?pattern=&ast;/&mc&run=<run>[#<treename>]
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
 *   alien:///some/directory?pattern=&ast;/&run=<run>[#<treename>
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
 *   alien:///alice/data/<year>/<period>?pattern=ESDs/pass<no>/AOD<vers>/&ast;/&run=<run>[#<treename>]
 * @endcode
 *
 * @section train_setup_other Other features  
 *
 * @subsection train_setup_aux Auxillary libraries, sources, and files
 *
 * Auxillary libraries should be loaded using 
 *
 * - Helper::LoadLibrary(const char*)
 *
 * where the argument is the name of the library
 *
 * If the train needs additional files, say a script for setting up
 * the tasks, or some data file, it can be passed on the the
 * PROOF/Grid workers using the member functions 
 *
 * - Helper::LoadAux(const char*)
 * - Helper::LoadSource(const char*,bool)
 * 
 * @subsection train_setup_overload Overloading the behaviour 
 *
 * The base class TrainSetup tries to implement a sensible setup for a
 * given type of analysis, but some times a particular train needs a
 * bit of tweaking.  One can therefore overload the following functions 
 * 
 * - TrainSetup::CreateGridHandler()
 * - TrainSetup::CreateInputHandler(EType)
 * - TrainSetup::CreateMCHandler(UShort_t,bool)
 * - TrainSetup::CreateOutputHandler(UShort_t)
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
 * ParUtilities::MakeScriptPAR to make a PAR file of the script and use
 * that to make a library loaded on the workers and then generate an
 * object of our task defined in the script.
 *
 * @include MyTrain.C
 *
 * This can allow for fast development and testing of analysis tasks
 * without having to wait for official tasks and builds of all of
 * AliROOT
 *
 * @section train_setup_impl Implementation details
 *
 * @subsection train_setup_imp_helper Helpers 
 *
 * The specifics of the each possible execution environment and input
 * is handled by sub-classes of the base class Helper. Each of these
 * helpers define
 *
 * - URI options. 
 * - Steps to be done before the tasks are added to the train 
 * - How to load libraries, additional scripts and files 
 * - Steps to be done after the setup of tasks 
 * - How to execute the analysis 
 * 
 * Currently defined helpers are 
 *
 * - LocalHelper for local jobs 
 * - ProofHelper for jobs running on a PROOF farm 
 * - LiteHelper for jobs running in a PROOF-Lite session 
 * - AAFHelper Special kind of ProofHelper for jobs running on AAFs 
 * - AFFPluginHelper As AAFHelper, but uses the AliEn plugin 
 * - GridHelper for Grid jobs 
 */
//
// EOF
//
