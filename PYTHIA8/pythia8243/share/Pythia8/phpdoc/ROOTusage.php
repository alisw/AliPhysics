<html>
<head>
<title>ROOT usage</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='ROOTusage.php'>
 
<h2>ROOT usage</h2> 
<ol id="toc">
  <li><a href="#section0">Standalone usage</a></li>
  <li><a href="#section1">Storing partial PYTHIA events in ROOT trees</a></li>
  <li><a href="#section2">PYTHIA as a plugin to ROOT</a></li>
</ol>

 
Many PYTHIA users wish to use <a href="http://root.cern.ch/">ROOT</a> 
to produce histograms, or even to run PYTHIA as a plugin to ROOT. 
This is possible. It is not a task supported by the PYTHIA team, 
however. All issues involving ROOT usage should be directed to the 
ROOT team, or to the local support team of your collaboration. 
Below some helpful hints have been collected. The text is based on 
contributions by Rene Brun, Andreas Morsch and Axel Naumann. 
Another example may be found in the 
<a href="http://vincia.hepforge.org">VINCIA</a> 
add-on program for parton showers, but this should also work for 
a PYTHIA standalone run. 
<br/><br/> 
Note that in all that follows, a Linux-type system with a Bash shell 
and GNU Make is assumed. In particular, for Mac OS X, the 
<code>LD_LIBRARY_PATH</code> should be replaced with 
<code>DYLD_LIBRARY_PATH</code> and the extension for shared libraries 
<code>.so</code> should be replaced with <code>.dylib</code>. 
 
<br/><br/><hr/> 
<a name="section0"></a> 
<h3>Standalone usage</h3> 
 
One can perform the generation and analysis of events in a completely 
standalone fashion, and only use ROOT to process the completed events. 
Two example programs are provided in the <code>examples</code> 
directory, with details provided below. 
 
The examples assume that ROOT is installed, that you have run 
<pre> 
    ./configure --with-root=root-installation-directory 
</pre> 
where you have to specify which is the ROOT installation directory, 
and subsequently run <code>make</code>. More fine-grained options are 
available with <code>configure</code>, if need be. 
 
<h4>Histogramming with ROOT</h4> 
 
An example of histogramming with ROOT is provided in 
<code>examples/main91.cc</code>. It may be compiled and run 
just like the other example programs. After PYTHIA has run, a ROOT 
histogram of the charged multiplicity in the events will be shown. 
This is now stored in the  <code>hist.root</code> file. If you can 
make this example work, the road should be open to do the same for 
all other histogramming needs. Specifically, you need to edit the 
<code>examples/Makefile</code> file to add the other programs to 
link as <code>main91.cc</code> currently does. 
 
<h4>Storing PYTHIA events in ROOT trees</h4> 
 
Instead of only generating histograms, it is possible to store entire 
PYTHIA events in ROOT trees. The <code>examples/main92</code> code 
provides an example of this and is comprised of the following files: 
<ul> 
  <li><code>main92.cc</code> is the main example program showing how 
  PYTHIA events can be stored in ROOT trees;</li> 
  <li><code>main92LinkDef.h</code> is used by Makefile to generate the 
  dictionary for all PYTHIA classes involved in the IO, as needed for 
  the example; and</li> 
  <li><code>main92.h</code> is a small include declaring the 
  <code>Pythia8</code> namespace as default.</li> 
</ul> 
 
<br/> 
The example may be compiled and run with as usual. Afterwards, the new 
<code>pytree.root</code> file will contain the PYTHIA events. Note 
that files can become quite large when many events are generated. To 
open these files within the ROOT interpreter the PYTHIA class 
dictionary must be loaded, <code>.L main92.so</code>. In compiled 
code, the PYTHIA class dictionary <code>main92.so</code> must be 
linked against, to either read or write PYTHIA events to a ROOT file. 
 
<h4>Error notice</h4> 
 
It appears that ROOTCINT cannot handle the <code>dlfcn.h</code> header 
in the current ROOT version. If you run into this problem with your 
ROOT installation, you could try to insert the following lines in 
your <code>PythiaStdlib.h</code> file: 
<pre> 
   // Stdlib header file for dynamic library loading. 
   #ifndef __CINT__ 
   #define dlsym __ 
   #include &lt;dlfcn.h&gt; 
   #undef dlsym 
   #endif 
</pre> 
 
<br/><hr/> 
<a name="section1"></a> 
<h3>Storing partial PYTHIA events in ROOT trees</h3> 
 
Instead of storing full PYTHIA events in ROOT trees, a common user case 
is to store only track information relevant to a particular analysis. 
The resulting ROOT trees will then be what is often referred to as 
"n-tuples". 
The advantage of this over the above method is a significant reduction 
of disk space used, as well as the possibility to construct trees 
resembling those familiar from the experiments' central MC production. 
 
The <code>examples/main93</code> example provides this - among other - 
functionality. As for the above example, it is split up in several 
files. 
<ul> 
  <li><code>main93.cc</code> is the main example program;</li> 
  <li><code>main93.cmnd</code> is a sample input command file;</li> 
  <li><code>main93LinkDef.h</code> is used by Makefile to generate the 
  dictionary for only the used PYTHIA classes involved in the IO, for 
  the example; and</li> 
  <li><code>main93.h</code> defines a "track" and an "event" class where 
  relevant event -and track information is defined.</li> 
</ul> 
 
<h4>Compiling the example</h4> 
 
The <code>main93</code> example is compatible with ROOT v.6 and above. 
One should have a working installation of ROOT, and then configure PYTHIA 
with: 
<pre> 
    ./configure --with-root=root-installation-directory 
</pre> 
One can then compile <code>main93</code> with the usual: 
<pre> 
    make main93 
</pre> 
provided that all ROOT paths are set correctly by eg. running: 
<pre> 
    source root-installation-directory/bin/thisroot.sh 
</pre> 
 
<h4>Running the example</h4> 
 
The <code>main93</code> example can be run with several command line options. 
Running: 
<pre> 
    ./main93 -h 
</pre> 
will display a help text showing these options. 
To produce events, the user needs to supply a command file with option 
<code>-c COMMAND-FILE.cmnd</code>. The example command file 
<code>main93.cmnd</code> is a good starting point. The crucial command to 
output ROOT trees is to set <code>Main:writeROOT = on</code>. 
 
The ROOT file will be named <code>pythia.root</code> per default. This can 
be changed by appending <code>-o ONAME</code> on the command line. 
 
<h4>Changing the event information</h4> 
 
The header file <code>main93.h</code> defines a simple event class and track 
class, which in turn defines the information stored to the tree. If a user 
wants to change this, either by adding more track information or imposing 
cuts corresponding to detector acceptance (thus reducing the file size), 
this can be done directly in this header file. Both the track class and the 
event class has <code>init</code> functions returning a boolean value, and 
by returning <code>false</code>, the track/event is rejected. 
The <code>main93</code> example must be recompiled after making any changes to 
the header file. 
 
<br/><br/><hr/> 
<a name="section2"></a> 
<h3>PYTHIA as a plugin to ROOT</h3> 
 
In more ROOT-centric applications, PYTHIA can be run as a ROOT plug-in. 
This requires a version of ROOT that has been 
<a href="http://root.cern.ch/drupal/content/installing-root-source"> 
installed from source</a>. The reason is that the interfaces depend on 
PYTHIA header files that are not distributed with ROOT. Installing ROOT 
is not more difficult than the PYTHIA installation, and some 
guidelines are provided below. 
 
<h4>Installation</h4> 
 
To be run as a plugin, PYTHIA must be compiled as a shared library. 
This is achieved by running the PYTHIA <code>configure</code> script 
with the <code>--enable-shared</code> option before <code>make</code> 
is run.<br/><br/> 
 
Define an environment variable for the path to your 
PYTHIA installation directory 
<pre> 
    export PYTHIA8=path_to_PYTHIA8_installation 
</pre> 
Before compiling ROOT, 
<a href="http://root.cern.ch/drupal/content/installing-root-source"> 
configure ROOT</a> by running the <code>configure</code> command 
including the following options 
<pre> 
    --enable-pythia8 
    --with-pythia8-incdir=$PYTHIA8/include/Pythia8 
    --with-pythia8-libdir=$PYTHIA8/lib 
</pre> 
In case ROOT has already been compiled before, it will only recompile 
the PYTHIA module and build the library <code>libEGPythia8</code>. 
 
<h4>Interfaces</h4> 
 
When running PYTHIA as a plugin, the exact interface structure becomes 
very relevant. ROOT provides two simple interfaces (wrappers) for 
PYTHIA 8. The code for these interfaces are located in 
<pre> 
    path_to_ROOT_source/montecarlo/pythia8 
</pre> 
<br/> 
The two interfaces are 
<ul> 
  <li><code><a href="http://root.cern.ch/root/html/TPythia8.html"> 
  TPythia8</a></code> is an implementation of the 
  <code><a href="http://root.cern.ch/root/html/TGenerator.html"> 
  TGenerator</a></code> interface for PYTHIA 8.<br/> 
  It allows you to use PYTHIA within a ROOT macro or as a plug-in 
  for a general-purpose particle generator based on this interface. The 
  main methods of the interface are 
  <ul> 
    <li><code>GenerateEvent()</code> which triggers the 
    generation of the next event, and </li> 
    <li><code>ImportParticles(TClonesArray* particles)</code> 
    which copies the native PYTHIA stack into a 
    <code><a href="http://root.cern.ch/root/html/TClonesArray.html"> 
    TClonesArray</a></code> of 
    <code><a href="http://root.cern.ch/root/html/TParticle.html"> 
    TParticles</a></code>.</li> 
  </ul> 
 
  In addition, some methods that are directly related to corresponding 
  PYTHIA methods are implemented 
  <ul> 
    <li><code>ReadString(const char* string)</code> &rarr; 
    <code>readString(...)</code></li> 
    <li><code>ReadConfigFile(const char* string)</code> &rarr; 
    <code>readFile(...)</code></li> 
    <li><code>Initialize(int idAin, int idBin, double ecms)</code> &rarr; 
    <code>init()</code> 
    <br/>Warning: this method will have to be updated for the 8.2 version! 
    </li> 
    <li><code>EventListing()</code> &rarr; 
    <code>event.list()</code></li> 
    <li><code>PrintStatistic()</code> &rarr; 
    <code>stat()</code> 
    <br/><b>Warning:</b> this method will have to be updated for 
    the 8.2 version! 
    </li> 
  </ul> 
 
  These methods provide already the basic PYTHIA functionality 
  interactively from the ROOT command line. However, this does not mean 
  that the usage of PYTHIA from within ROOT is restricted to these methods. 
  In compiled code, one can always obtain a pointer to the 
  <code>Pythia</code> instance e.g. 
  <pre> 
    TPythia8        *tp = new TPythia8(); 
    Pythia8::Pythia *p  = tp->Pythia8();</pre> 
  giving access to the full PYTHIA functionality. To access this 
  functionality in the CINT interpreter see the "Advanced usage" 
  section below.</li> 
 
  <li><code><a href="http://root.cern.ch/root/html/TPythia8Decayer.html"> 
  TPythia8Decayer</a></code> is an implementation of the 
  <code><a href="http://root.cern.ch/root/html/TVirtualMCDecayer.html"> 
  TVirtualMCDecayer</a></code> interface.<br/> 
  It allows you to use PYTHIA as a plug-in decayer for simulation 
  frameworks based on the Virtual Monte Carlo 
  (<a href="http://root.cern.ch/drupal/content/vmc">VMC</a>) interface 
  classes. The main methods of the interface are 
  <ul> 
    <li><code>TPythia8Decayer::Init()</code> for initialisation,</li> 
    <li><code>TPythia8Decayer::Decay(Int_t pdg, TLorentzVector* p)</code> 
    to decay a particle with PDG code <code>pdg</code> and 
    <a href="http://root.cern.ch/root/html/TLorentzVector.html"> 
    4-momentum</a> <code>p</code>, and </li> 
    <li><code>ImportParticles(TClonesArray* particles)</code> 
    to retrieve the decay products as 
    <code><a href="http://root.cern.ch/root/html/TParticle.html"> 
    TParticles</a></code> in the 
    <code><a href="http://root.cern.ch/root/html/TClonesArray.html"> 
    TClonesArray</a> particles</code>.</li> 
  </ul></li> 
</ul> 
 
<h4>An example</h4> 
 
A <a href="http://root.cern.ch/root/html/tutorials/pythia/pythia8.C.html"> 
basic example</a> for generating minimum-bias events with PYTHIA 8 inside 
a ROOT macro, and filling some histograms with the kinematics of the 
final-state particles is provided in either of the locations below 
<pre> 
    /path_to_ROOT_source/tutorials/pythia/pythia8.C 
    /path_to_ROOT_installation/share/doc/root/tutorials/pythia/pythia8.C 
</pre> 
<br/> 
Note that before executing this script 
<ul> 
  <li>the environment variables <code>PYTHIA8</code> and 
  <code>PYTHIA8DATA</code> must be setup correctly e.g. 
  <pre> 
    export PYTHIA8=/path_to_PYTHIA_installation 
    export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc 
  </pre></li> 
  <li>your LD_LIBRARY_PATH must contain the location of the 
  PYTHIA 8 shared library, e.g. 
  <pre> 
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path_to_PYTHIA8_installation/lib 
  </pre> 
  </li> 
</ul> 
The script can then be run with ROOT 
<pre> 
    root pythia8.C 
</pre> 
After execution, ROOT will display some histograms from the event 
generation. 
 
<h4>Advanced usage</h4> 
 
To access the full PYTHIA functionality from the CINT interpreter, 
a ROOT dictionary must be created. Currently that option has not been 
implemented as a standard option for PYTHIA 8.2, but it should be in 
the same spirit as what can be found in the 8.1 <code>rootexamples</code> 
directory. Also note that one dictionary is found in the 
<code>examples/main92LinkDef.h</code> file. 
 
This may then be loaded in ROOT giving full access to the full PYTHIA 8 
functionality, e.g. in an interactive session 
<pre> 
    gSystem->Load("path_to_PYTHIA8_installation/rootexamples/pythiaDict"); 
    Pythia8::Pythia *p = new Pythia8::Pythia(); 
    p->readString("SoftQCD:nonDiffractive = on"); 
</pre> 
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
