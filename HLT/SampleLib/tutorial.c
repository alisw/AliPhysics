/************************************************************************
**
** ALICE HLT project
** Copyright (c) 2005
**
** This file is property of and copyright by the Experimental Nuclear 
** Physics Group, Dep. of Physics and Technology
** University of Bergen, Norway, 2004
** This file has been written by Matthias Richter,
** Matthias.Richter@ift.uib.no
**
** Permission to use, copy, modify and distribute this software and its  
** documentation strictly for non-commercial purposes is hereby granted  
** without fee, provided that the above copyright notice appears in all  
** copies and that both the copyright notice and this permission notice  
** appear in the supporting documentation. The authors make no claims    
** about the suitability of this software for any purpose. It is         
** provided "as is" without express or implied warranty.                 
**
*************************************************************************/

/** @file   tutorial.c
    @author Matthias Richter
    @date   
    @brief  HLT examples and tutorial. */

/** 
@defgroup alihlt_tutorial HLT examples and tutorial

The HLT analysis components can be run either in the AliRoot
framework (simulation and/or reconstruction) or the HLT online
framework.

We think of the HLT as a 'black box' with data input and output. In
addition there is access to calibration data from OCDB (or the local
HLT copy HCDB). All components can only work on the data they get as
input. As the different detector algorithms/components will run in 
separate processes and even on different machines, no data exchange
is possible via global data structures and variables.

HLT chains in the AliRoot framework are described by means of
AliHLTConfiguration.

-# @ref tut_hltsystem
   -# @ref tut_load_libraries
   -# @ref tut_dummy_chain
   -# @ref tut_tpc_sector
-# @ref tut_simulation
-# @ref tut_reconstruction
   -# @ref tut_module_agent
   -# @ref tut_reconstruction_sample
   -# @ref tut_reconstruction_custom
-# @ref tut_alirawreaderhlt
-# @ref tut_macros 

<br>
<hr width="75%">
<br>
@section tut_hltsystem Running Components in the HLT System

@subsection tut_load_libraries Library setup
Component libraries must be loader via the AliHLTComponentHandler
or AliHLTSystem::LoadComponentLibraries. You can run the following
macro from the AliRoot promt.
<pre>
{
  AliHLTSystem gHLT;
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");
}
</pre>

<br>
@subsection tut_dummy_chain Example: Running a dummy chain
The simplest chain consists of a publisher component, a processor
and a data sink. The AliHLTDummyComponent is a sample component
which just copies a fraction of the input data to the output.
You can run the following macro from the AliRoot promt.
<pre>
{
  AliHLTSystem gHLT;
  gHLT.LoadComponentLibraries("libAliHLTUtil.so libAliHLTSample.so");
  // The AliHLTFilePublisher (component Id \em 'FilePublisher' provides
  // the given file (see AliHLTFilePublisher for more options) to the
  // subsequent components in the chain.
  AliHLTConfiguration publisher("fp1", "FilePublisher", NULL, "-datatype 'DUMMYDAT' 'SMPL' -datafile some-data.dat");

  // The AliHLTDummyComponent (Id \em 'Dummy') just forwards a certain
  // fraction of the input to the output or just repeats the input data
  // if percentage > 100
  AliHLTConfiguration copy("cp", "Dummy", "fp1", "output_percentage 80");

  // The AliHLTFileWriter (Id 'FileWriter') is a data sink. It writes
  // all incoming data blocks to files. Several options available.
  AliHLTConfiguration sink1("sink1", "FileWriter", "cp", NULL);

  // here you specify the top most configuration of the chain. The
  // configuration depends on all the parents. The task lisy is build
  // according to that.
  gHLT.BuildTaskList("sink1");
  gHLT.Run();
}
</pre>
@note You have to specify a real file name instead of \em some-data.dat

<br>
@subsection tut_tpc_sector Example: One sector of the TPC
This example builds an analysis chain for TPC sector 0. It works on 
simulated data and assumes the ddl files to be present in the current
directory.
<pre>
{
  AliHLTSystem gHLT;
  // load the component library
  gHLT.LoadComponentLibraries("libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so");

  // data source components
  AliHLTConfiguration fp0("fp0", "FilePublisher", NULL, "-datatype 'DDL_RAW ' 'TPC ' -dataspec 0x00000000 "
                                                        "-datafile TPC_768.ddl");
  AliHLTConfiguration fp1("fp1", "FilePublisher", NULL, "-datatype 'DDL_RAW ' 'TPC ' -dataspec 0x00000101 "
                                                        "-datafile TPC_769.ddl");
  AliHLTConfiguration fp2("fp2", "FilePublisher", NULL, "-datatype 'DDL_RAW ' 'TPC ' -dataspec 0x00000202 "
                                                        "-datafile TPC_840.ddl");
  AliHLTConfiguration fp3("fp3", "FilePublisher", NULL, "-datatype 'DDL_RAW ' 'TPC ' -dataspec 0x00000303 "
                                                        "-datafile TPC_841.ddl");
  AliHLTConfiguration fp4("fp4", "FilePublisher", NULL, "-datatype 'DDL_RAW ' 'TPC ' -dataspec 0x00000404 "
                                                        "-datafile TPC_842.ddl");
  AliHLTConfiguration fp5("fp5", "FilePublisher", NULL, "-datatype 'DDL_RAW ' 'TPC ' -dataspec 0x00000505 "
                                                        "-datafile TPC_843.ddl");

  // cluster finders
  AliHLTConfiguration cf0("cf0", "TPCClusterFinderPacked", "fp0", "pp-run rawreadermode 4 timebins 446");
  AliHLTConfiguration cf1("cf1", "TPCClusterFinderPacked", "fp1", "pp-run rawreadermode 4 timebins 446");
  AliHLTConfiguration cf2("cf2", "TPCClusterFinderPacked", "fp2", "pp-run rawreadermode 4 timebins 446");
  AliHLTConfiguration cf3("cf3", "TPCClusterFinderPacked", "fp3", "pp-run rawreadermode 4 timebins 446");
  AliHLTConfiguration cf4("cf4", "TPCClusterFinderPacked", "fp4", "pp-run rawreadermode 4 timebins 446");
  AliHLTConfiguration cf5("cf5", "TPCClusterFinderPacked", "fp5", "pp-run rawreadermode 4 timebins 446");

  // tracker
  AliHLTConfiguration tracker("tracker", "TPCSliceTracker", "cf0 cf1 cf2 cf3 cf4 cf5", "-pp-run -bfield 0.5");

  // the data sink component
  AliHLTConfiguration writer("writer", "TPCEsdWriter", "tracker", "-datafile AliHLTTPCESDs.root");

  gHLT.BuildTaskList("writer");
  gHLT.Run();
}
</pre>

<br>
<hr width="75%">
<br>

@section tut_reconstruction AliRoot reconstruction
The integration into the AliRoot reconstruction works via the
@ref AliHLTReconstructor plugin. The intention is to run HLT analysis
chains in AliRoot in the same way as in the online framework, i.e.
the full components are run also from the offline framework rather
than just the algorithm hooked on by a special interface class.
By this one achieves the highest possible compatibility.

The AliRoot reconstruction consists mainly of three steps:
-# local event reconstruction: this is usually the place for digit/raw
data conversion to clusters/space points. All the events are processed
at once. If HLT analysis chain are executed in the local event
reconstruction, the chain must contain an output recorder as the ESD
is filled on an event by event basis and the corresponding method called
later.
-# event reconstruction: this is the reconstruction on an event by event
basis. Immediately after the reconstruction, the FillESD method is
called.
-# ESD fill: the reconstructed event is written to the ESD.

Regarding the HLT, all analysis is supposed to run on-line on the HLT farm.
Thus, only the processing of the HLTOUT data is necessary during the
default reconstruction. However, it is possible to run HLT chains embedded
into AliReconstruction mainly for the purpose of debugging and the
development cycle.

The actual chains to be run and the HLTOUT handlers to be applied depend
on the HLT library modules which are loaded to the system. There is a
default collection of libraries defined in AliHLTSystem::fgkHLTDefaultLibs.
The default libraries are loaded if nothing else is specified. The libraries
implement \em agents (AliHLTModuleAgent childs) describing the properties
of a module.

A specific library can be chosen like (provided you have a simulated
event in the current directory):
<pre>
{
  AliReconstruction rec;                 // the reconstruction instance
  rec.SetRunLocalReconstruction("HLT");  // run local rec only for HLT 
  rec.SetRunTracking("");                // switch off tracking
  rec.SetFillESD("HLT");                 // 
  rec.SetOption("HLT", "libAliHLTSample.so loglevel=0x7c");
  rec.Run();
}
</pre>

@subsection tut_module_agent The Module Agent
Each component library has to implement a module agent in order to be
hooked up to the AliRoot reconstruction or simulation. The agent defines
the features of the libraries and the configurations to be run during the
different steps of the reconstruction. The agent 
- can register all components of the library. This is an 
  alternative to the component registration via global objects (see 
  @ref alihltcomponent-handling).
- registers HLT configurations (see @ref AliHLTConfiguration)
- specifies the configurations to be run
- specifies additional component libraries required to run the
  configurations.
- provides a preprocessor (see AliHLTModulePreprocessor /
			   AliHLTPreprocessor)
- provides handlers and handler descriptions for HLTOUT data blocks.

Finally, one global object of the module agent has to be specified in
the source code. All registration and integration into the HLT system
is carried out automatically.

@see 
    @ref AliHLTModuleAgent for the interface description <br>
    @ref AliHLTAgentSample for a sample implementation

@subsection tut_reconstruction_sample The sample library
The libAliHLTSample library provides examples how to implement the
library agent (@ref AliHLTAgentSample), how to add configurations and
define HLT chains for reconstruction. 

The sample library is not part of the default libraries loaded by the
HLT steering during reconstruction. The example can be run by the
following macro macro above.

The agent defines the following chains:
-# a simple data copying consisting of a
   - @ref AliHLTFilePublisher  publishes some data generated before in /tmp
   - @ref AliHLTDummyComponent copies a fraction of the incoming data
   - @ref AliHLTFileWriter     writes the data input to a file
-# digit publishing from the TPCloader <br>
   This chain illustrates how data can be published from the AliRunLoader
   in order to be processed by another component (not in the sample chain).
   Finally, the @ref AliHLTSampleOfflineSinkComponent is component which is
   the backend and has again the AliRoot structures.
   - @ref AliHLTLoaderPublisherComponent
   - @ref AliHLTSampleOfflineSinkComponent

In the same way any other component library can be integrated into the
AliRoot reconstruction.

@subsection tut_reconstruction_custom Running a custom HLT chain
The default configurations from the library modules can be overridden by
custom configurations by means of options specified to AliReconstruction.
- <tt>config=\em macro</tt><br> a configuration macro. The macro is a normal
  ROOT macro defining HLT component configurations by means of
  AliHLTConfiguration
- <tt>chains=\em chains</tt><br> a comma separated list of chains to be run.
  A chain is defined by the topmost configuration.

\b Note: The file publisher needs a file to read, either you replace
\em some-data.dat with the path of an existing file or just create a
dummy file in the current working directory. Futhermore, there has to be at
least one simulated event since AliReconstruction relies on a couple of files
in the folder.
<pre>
{
  AliReconstruction rec;                 // the reconstruction instance
  rec.SetInput("./");                    // to be independent of galice.root
  rec.SetLoadAlignFromCDB(kFALSE);
  rec.SetFillTriggerESD(kFALSE);
  rec.SetRunQA(kFALSE);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT");  // run local rec only for HLT 
  rec.SetRunTracking("");                // switch off tracking
  rec.SetFillESD("HLT");                 // 
  rec.SetOption("HLT", "libAliHLTSample.so libAliHLTUtil.so "
                       "config=$ALICE_ROOT/HLT/exa/conf-sample.C "
                       "chains=sink");
  //rec.SetEventRange(0,0);
  rec.Run();
}
</pre>

@see
- conf-sample.C

<br>
<hr width="75%">
<br>

@section tut_simulation AliRoot simulation
In order to simulate the behavior of HLT analysis chains and to
include this functionality, HLT reconstruction can be embedded
into AliRoot simulation. As a matter of fact, HLT always reconstructs
data, <em><b>HLT simulation</b></em> means <em><b>HLT reconstruction
embedded into AliRoot</b></em>.

The HLT simulation is run at the last step of the AliSimulation, the
setup to be run depends on the available plugins as described in section
@ref tut_reconstruction. The options for the HLT simulation can be set
with the <tt>AliSimulation::SetRunHLT</tt> function.
<pre>
  AliSimulation sim;
  ...
  sim.SetRunHLT("libAliHLTSample.so loglevel=0x7c");
</pre>

Options:
- <tt>config=\em macro</tt><br> a configuration macro. The macro is a normal
  ROOT macro defining HLT component configurations by means of
  AliHLTConfiguration
- <tt>chains=\em chains</tt><br> a comma separated list of chains to be run.
  A chain is defined by the topmost configuration.
- <tt>rawfile=\em chains</tt><br> provide a raw reader to the HLT simulation.
  Some chains work solely on raw data. The data needs to be simulated before
  and a RawReader is created internally to provide the data to the source
  components.
- <tt>loglevel=\em 0x7c</tt><br> default loglevel is 0x79, only Warnings and
  higher are printed. 0x7c also makes the Info messages visible.

@see
  - sim-hlt-rawddl.C for example


<br>
<hr width="75%">
<br>

@section tut_alirawreaderhlt Replacing reconstruction input by data from the HLT
The HLTOUT data can contain data blocks which obey exactly the raw DDL 
format of a detector. E.g. selective readout components or loss-less
compression algorithms allow to provide a sub sample of the original data.
All data from the HLT is transferred via the 10 HLT DDL links, a redirection
mechanism is necessary to replace the original detector input by the data
from HLTOUT. The replacements works by means of the AliRawReaderHLT and
needs the following modules:
-# Implementation of an AliHLTOUTHandlerEquId child class<br>
   A handler of this type is necessary to determine the equipment Id of a
   data block from the data type and specification.
<pre>
  class AliHLTSampleRawDataHandler : public AliHLTOUTHandlerEquId {
  public:
    // ... constructors and destructor
    &nbsp;
    // overloaded AliHLTOUTHandlerEquId::ProcessData(AliHLTOUT*)
    int ProcessData(AliHLTOUT* pData);
  };
</pre>
  Alternatively the AliHLTOUTHandlerEquId can be used directly. It implements
  a default processing.
-# Adjust module agent
   The module agent needs to provide the handler and the description of the
   handler and has to implement the following functions: 
<pre>
   // see AliHLTAgentSample::GetHandlerDescription()
   int GetHandlerDescription(AliHLTComponentDataType,
                             AliHLTUInt32_t,
                             AliHLTOUTHandlerDesc&) const;
&nbsp;   
   // see AliHLTAgentSample::GetOutputHandler()
   AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType,
                                      AliHLTUInt32_t);
   &nbsp;   
   // see AliHLTAgentSample::DeleteOutputHandler()
   int DeleteOutputHandler(AliHLTOUTHandler*);
</pre>
-# Set the HLT input
   The AliReconstruction class handles the redirection transparently by
   use of the AliRawReaderHLT.
<pre>
   AliReconstruction rec;
   // ....
   rec.SetUseHLTData("ITSSDD");
</pre>
-# Run
   Run the reconstruction as normal

@see
   - AliHLTReconstructor
   - AliRawReaderHLT
   - rec-from-hltout.C

<br>
<hr width="75%">
<br>

@section tut_macros Example macros
 */


/* note pad

Making a new module/library

Automatic ROOT dictionary generation:
The automatic ROOT dictionary generation relies on the rule, that the main class
of a header file has the same name as the file (except the prefix).

Troubleshooting:
Error: link requested for unknown class <class name> <library>-LinkDef.h:<line no>
most likely there is no class <class name> defined in the header file <class name>.h*

 */
#error Not for compilation
//
// EOF
//
