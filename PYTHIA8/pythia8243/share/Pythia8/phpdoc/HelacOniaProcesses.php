<html>
<head>
<title>HelacOnia Processes</title>
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

<form method='post' action='HelacOniaProcesses.php'>
 
<h2>HelacOnia Processes</h2> 
 
HelacOnia [<a href="Bibliography.php#refSha15" target="page">Sha15</a>] is an external package which provides 
automated calculations for heavy quarkonia production using NRQCD, 
similar in style to <?php $filepath = $_GET["filepath"];
echo "<a href='MadGraph5Processes.php?filepath=".$filepath."' target='page'>";?>MadGraph5</a> 
and the extension MadOnia, which is only available for MadGraph4. This 
can be useful when additional quarkonia processes other than the 
internal processes provided in <?php $filepath = $_GET["filepath"];
echo "<a href='OniaProcesses.php?filepath=".$filepath."' target='page'>";?>Onia</a> 
are needed, including matrix elements which are not spin-averaged, as 
well as the ability to produce <i>n</i>-leg matrix elements beyond 
the leading tree-level diagrams. The HelacOnia code can be downloaded 
from 
<br><a href="http://helac-phegas.web.cern.ch/helac-phegas/helac-onia.html" 
target="page">http://helac-phegas.web.cern.ch/helac-phegas/helac-onia.html</a>, 
</br>where only version 2 and above is compatible with PYTHIA. 
 
<p/> 
Within HelacOnia, events can automatically be passed to PYTHIA for 
additional processing, e.g. showering, MPI, and 
hadronization. However, in many cases it may be simpler to produce 
HelacOnia events directly in PYTHIA. The <code>LHAupHelaconia</code> 
class provided in <code>Pythia8Plugins/LHAHelaconia</code> is designed 
to provide such utility. Here we will describe how this can be used to 
wrap the HelacOnia generator as a PYTHIA Les Houches interface. 
 
<p/> 
Of course, HelacOnia can also output files of parton-level events 
according to the <?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>LHEF</a> standard, 
that can be read in and processed further by PYTHIA 8. This is the 
most commonly used approach, and requires no further description here. 
 
<a name="section0"></a> 
<h3>HelacOnia executable inside PYTHIA</h3> 
 
The <code>Pythia::setLHAupPtr(LHAup* lhaUpPtr)</code> method allows 
a Pythia generator to accept a pointer to an object derived from the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>LHAup</a></code> base class. 
Such an object will be initialized from within Pythia, and be called 
repeatedly to generate the next parton-level event, using the LHA 
specification as a standard to transfer the relevant information back 
to Pythia. Properly constructed, the operation of an <code>LHAup</code> 
object thus is almost completely hidden from the user, and generates 
events almost like an ordinary internal Pythia process. 
 
<p/> 
The <code>LHAupHelaconia</code> is precisely such a class, derived from 
<code>LHAup</code>, that contains the code needed to wrap a 
HelacOnia executable. Thereby the generation of HelacOnia 
processes from within Pythia becomes straightforward. An explicit 
example is provided in <code>main35.cc</code>. We describe some of the 
key elements used there and in the general case. 
 
<a name="anchor1"></a>
<p/><strong> LHAupHelaconia::LHAupHelaconia(Pythia* pythia, string dir = &quot;helaconiarun&quot;, string exe = &quot;ho_cluster&quot;) &nbsp;</strong> <br/>
creates an instance of the <code>LHAupHelaconia</code> class. 
<br/><code>argument</code><strong> pythia </strong>  :  pointer to the <code>Pythia</code> instance, 
such that some of its facilities can be used inside the interface. 
   
<br/><code>argument</code><strong> dir </strong> (<code>default = <strong>helaconiarun</strong></code>) :  the name of the run 
directory, into which HelacOnia puts its (intermediate) results. 
   
<br/><code>argument</code><strong> exe </strong> (<code>default = <strong>ho_cluster</strong></code>) :  the name of the HelacOnia 
executable that <code>LHAupHelaconia</code> is meant to wrap. In addition 
it may be necessary to prepend the full pathname of the executable: 
<code>"(something)/HELAC-Onia-2.0.1/cluster/bin/ho_cluster"</code>. 
   
   
 
<a name="anchor2"></a>
<p/><strong> bool LHAupHelaconia::readString(string line) &nbsp;</strong> <br/>
allows the user to send commands to HelacOnia. 
<br/><code>argument</code><strong> line </strong>  :  the command to be sent to HelacOnia. For 
example, the following will produce <i>J/psi</i> events events from 13 TeV 
proton proton collisions: <br/><code>readString("generate u u~ > 
cc~(3S11) g");</code> <br/> A special case is the generation of 
colour-octet states. In PYTHIA these are evolved to colour-singlet 
states through the emission of a soft gluon with the mass splitting 
set by <code>Onia:massSplit</code>. To ensure the colour-octet states 
in HelacOnia are produced with the correct masses needed for this 
splitting, the specific colour-octet state for the process must be 
set. For example: 
<br/><code>readString("generate u u~ > cc~(3S18) g");</code> 
<br/>requires that the colour-singlet state into which the 
colour-octet state should decay be set. This could be set via: 
<br/><code>readString("set state = 443");</code> 
<br/>for the case where a final state <i>J/psi</i> is 
requested. Note that this is not a command passed to HelacOnia, but 
rather a command which PYTHIA uses to set the heavy quark mass in 
HelacOnia and then translate the HelacOnia output to the correct 
colour-singlet state. 
   
   
 
<a name="anchor3"></a>
<p/><strong> void LHAupHelaconia::setEvents(int events) &nbsp;</strong> <br/>
the number of events to generate per HelacOnia run. Normally does not 
need to be set, but defaults to 10000. 
   
 
<a name="anchor4"></a>
<p/><strong> void LHAupHelaconia::setSeed(int seed, int runs = 30081) &nbsp;</strong> <br/>
the random seed (sequence), normally not needed to be set explicitly. 
If the random seed is negative (default of -1), then the HelacOnia 
seed is taken as the Pythia parameter <code>"Random:seed"</code>, which 
must be greater than 0. If the maximum number of allowed runs is exceeded 
(default of 30081) an error is thrown. The seed for a HelacOnia run is set as: 
<br/> (random seed - 1) * (maximum runs) + (number of runs) + 1. 
<br/>HelacOnia can only handle random seeds up to 30081 * 30081. So, with 
this strategy, one can generate Pythia jobs with seeds from 1 to 30081, 
with each job running HelacOnia less than 30081 times, and ensure a fully 
statistically independent sample. If more than 30081 jobs are needed, then 
the maximum allowed runs can be lowered accordingly, and if need be, 
setEvents can be used to increase the number of events generated per run. 
   
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
