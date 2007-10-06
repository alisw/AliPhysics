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

-# @ref tut_hltsystem
   -# @ref tut_load_libraries
   -# @ref tut_dummy_chain
   -# @ref tut_tpc_sector
-# @ref tut_reconstruction
   -# @ref tut_module_agent
   -# @ref tut_reconstruction_sample

<br>
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
  AliHLTConfiguration publisher("fp1", "FilePublisher", NULL, "-datafile some-data.dat");
  AliHLTConfiguration copy("cp", "Dummy", "fp1", "output_percentage 80");
  AliHLTConfiguration sink1("sink1", "FileWriter", "cp", NULL);
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
  gHLT.LoadComponentLibraries("libAliHLTUtil.so libAliHLTTPC.so");

  // data source components
  AliHLTConfiguration fp0("fp0", "FilePublisher", NULL, "-datatype DDL_RWPK 'TPC ' -dataspec 0x00000000"
			                                "-datafile TPC_768.ddl");
  AliHLTConfiguration fp1("fp1", "FilePublisher", NULL, "-datatype DDL_RWPK 'TPC ' -dataspec 0x00000101"
                                                        "-datafile TPC_769.ddl");
  AliHLTConfiguration fp2("fp2", "FilePublisher", NULL, "-datatype DDL_RWPK 'TPC ' -dataspec 0x00000202"
                                                        "-datafile TPC_840.ddl");
  AliHLTConfiguration fp3("fp3", "FilePublisher", NULL, "-datatype DDL_RWPK 'TPC ' -dataspec 0x00000303"
			                                "-datafile TPC_841.ddl");
  AliHLTConfiguration fp4("fp4", "FilePublisher", NULL, "-datatype DDL_RWPK 'TPC ' -dataspec 0x00000404"
			                                "-datafile TPC_842.ddl");
  AliHLTConfiguration fp5("fp5", "FilePublisher", NULL, "-datatype DDL_RWPK 'TPC ' -dataspec 0x00000505"
			                                "-datafile TPC_843.ddl");

  // cluster finders
  AliHLTConfiguration cf0("cf0", "TPCClusterFinderPacked", "fp0", "pp-run rawreadermode 4 timebins 446");
  AliHLTConfiguration cf1("cf1", "TPCClusterFinderPacked", "fp1", "pp-run rawreadermode 4 timebins 446");
  AliHLTConfiguration cf2("cf2", "TPCClusterFinderPacked", "fp2", "pp-run rawreadermode 4 timebins 446");
  AliHLTConfiguration cf3("cf3", "TPCClusterFinderPacked", "fp3", "pp-run rawreadermode 4 timebins 446");
  AliHLTConfiguration cf4("cf4", "TPCClusterFinderPacked", "fp4", "pp-run rawreadermode 4 timebins 446");
  AliHLTConfiguration cf5("cf5", "TPCClusterFinderPacked", "fp5", "pp-run rawreadermode 4 timebins 446");

  // tracker
  AliHLTConfiguration tracker("tracker", "TPCSliceTracker", "cf0 cf1 cf2 cf3 cf4 cf5", "pp-run bfield 0.5");

  // the data sink component
  AliHLTConfiguration writer("writer", "TPCEsdWriter", "tracker", "-datafile AliESDs.root");

  gHLT.BuildTaskList("writer");
  gHLT.Run();
}
</pre>
    
@section tut_reconstruction AliRoot reconstruction
The HLT analysis components can be run either in the AliRoot
reconstruction or the HLT online framework.
The integration into the AliRoot reconstruction works via the
@ref AliHLTReconstructor plugin. The intention is to run HLT analysis
chains in AliRoot in the same way as in the online framework, i.e.
the full components are run also from the offline framework rather
than just the algorithm hooked on by a special interface class.
By this one achieves the highest possible compatibility.

We think of the HLT as a 'black box' with data input and output. In
addition there is access to calibration data from OCDB (or the local
HLT copy HCDB). All components can only work on the data they get as
input. As the different detector algorithms/components will run in 
separated processes and even on different machines, no data exchange
is possible via global data structures and variables possible.

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

@subsection tut_module_agent The Module Agent
Each component library has to implement a module agent in order to be
hooked up to the AliRoot reconstruction system. The agent defines the
features of the libraries and the configurations to be run during the
different steps of the reconstruction. The agent 
- can register all components of the library. This is an 
  alternative to the component registration via global objects (see 
  @ref alihltcomponent-handling).
- registers HLT configurations (see @ref AliHLTConfiguration)
- specifies the configurations to be run
- specifies additional component libraries required to run the
  configurations.

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
following macro in AliRoot (provided you have a simulated event in the
current directory):
<pre>
{
  AliReconstruction rec;                 // the reconstruction instance
  rec.SetRunLocalReconstruction("HLT");  // run local rec only for HLT 
  rec.SetRunTracking("");                // switch off tracking
  rec.SetOption("HLT", "libAliHLTSample.so");
  rec.Run();
}
</pre>

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
