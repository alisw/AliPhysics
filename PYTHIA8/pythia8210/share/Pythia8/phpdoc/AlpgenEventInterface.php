 
<html>
<head>
<title>ALPGEN Event Interface</title>
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

<form method='post' action='AlpgenEventInterface.php'>
 
<h2>ALPGEN Event Interface</h2> 
 
This manual page describes the ALPGEN [<a href="Bibliography.php" target="page">Man03</a>] event interface 
for PYTHIA8.   While future versions of 
ALPGEN will be able to write out events in LHEF format, previous 
versions always output events in an ALPGEN native format (a combination 
of a ".unw" and a "_unw.par" file). The ALPGEN component of this code 
contains a reader for this native format (for unweighted events), as 
well as parameter reading for both ALPGEN native and LHE file formats. 
The reader was designed to work together with an implementation of 
the ALPGEN-style parton-jet matching <code>JetMatchingAlpgen</code> 
described on the <?php $filepath = $_GET["filepath"];
echo "<a href='JetMatching.php?filepath=".$filepath."' target='page'>";?>Jet Matching</a> 
page. However, it will also work with a implementation of the 
Madgraph-style [<a href="Bibliography.php" target="page">Alw11</a>] parton-jet matching 
<code>JetMatchingMadgraph</code> also described on the 
<?php $filepath = $_GET["filepath"];
echo "<a href='JetMatching.php?filepath=".$filepath."' target='page'>";?>Jet Matching</a> page. 
A sensible choice of <code>JetMatching</code> parameters is needed 
when using ALPGEN files with Madgraph-style matching and vice versa. 
 
<p/> 
It should be noted that all the functionality described here is provided 
through external routines, and therefore the presence of these features is 
dependent on the main program being used. This structure allows for the 
easy extensibility of the merging scheme. The files of interest are located 
in the <code>include/Pythia8Plugins/</code> subdirectory: 
<ul> 
<li> 
<code>GeneratorInput.h</code> : provides three classes for the reading of 
ALPGEN event and parameter files. <code>LHAupAlpgen</code> is an 
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?><code>LHAup</code></a> derived 
class for reading in ALPGEN native format event files. 
<code>AlpgenPar</code> is a class for the parsing of ALPGEN parameter 
files, making the information available through a simple interface. 
<code>AlpgenHooks</code> is a 
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?><code>UserHooks</code></a> derived class that 
provides the <code>Alpgen:*</code> options, described below. Further 
technical details of these classes are given at the end of this manual 
page. 
</li> 
<li> 
<code>main32.cc, main32.cmnd</code> : a sample main program and card 
file showing the usage of previous file and an MLM <code>UserHooks</code> 
class. In combination, it reads in a sample ALPGEN (or Madgraph) event file 
while performing the MLM merging procedure as implemented in ALPGEN 
(or as in Madgraph). Some commented-out sets of options are provided 
in the card file, which can be activated to try different merging setups. 
</li> 
<li> 
<code>main32.unw, main32_unw.par</code> : an ALPGEN format event and 
parameter file containing 100 W + 3 jet events. It is not feasible 
to package large event files with the PYTHIA distribution, but this 
sample is enough to show the different components in action. 
</li> 
</ul> 
 
<h2>ALPGEN main options</h2> 
 
These following options are provided by the AlpgenHooks class, 
which must be loaded for this functionality to be present 
 
<p/> 
ALPGEN event files that have been written out in LHEF format should be 
read in through the normal LHEF machinery 
(see <?php $filepath = $_GET["filepath"];
echo "<a href='BeamParameters.php?filepath=".$filepath."' target='page'>";?>Beam Parameters</a>). Files in 
ALPGEN's native format, instead, may be processed using the 
<code>Alpgen:file</code> option below. When using this option, the 
ALPGEN parameter file is stored in the PYTHIA Info object under the key 
<code>AlpgenPar</code>, see the "Header information" section of the 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>Event Information</a> manual page for 
more details. Processes not implemented by the PYTHIA 6 interface 
supplied with ALPGEN are also not implemented here. 
 
<p/> 
When reading in ALPGEN native event files, some momenta are shifted by 
the file reader to ensure energy-momentum conservation. The magnitude of 
these shifts should be small (around the MeV level in the worst case) 
and warnings will be produced if they are above a set threshold. A large 
number of warnings may signify unexpected behaviour and should 
potentially be investigated. It is also known that certain event 
classes, for example an event with both light and heavy <i>b</i> 
quarks may give rise to these warnings. 
 
<p/> 
The ALPGEN file reader supports the reading of the event and parameter 
files in gzip format with file extensions ".unw.gz" and "_unw.par.gz" 
respectively. This requires the use of external libraries, however, and 
the <code>README</code> file in the main directory contains instructions 
on how to enable this. 
 
<p/> 
All other <code>Alpgen:*</code> options apply to both LHE and native 
file formats, and include options to guide the MLM merging procedure 
based on the parameters that are read in with the events file. 
 
<br/><br/><table><tr><td><strong>Alpgen:file  </td><td></td><td> <input type="text" name="1" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
This option is used to read in ALPGEN format event files. Using this option 
overrides any previously set beam options inside PYTHIA. The path to the 
files, not including any file extension, should be provided e.g. for input 
files <i>input_unw.par</i> and <i>input.unw</i>, the value 
<i>input</i> should be used. 
   
 
<br/><br/><strong>Alpgen:setLightMasses</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
When switched on, <i>c</i> and <i>b</i> quark masses provided 
by ALPGEN are set in the PYTHIA 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>particle database</a>. 
Since ALPGEN may set these two masses to vanish, the parton shower 
programs have been provided with some protection, but other parts of 
the code may not be as fortunate. You should therefore only switch on 
this option if you know what you are doing. 
   
 
<br/><br/><strong>Alpgen:setHeavyMasses</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
When switched on, <i>t</i>, <i>Z</i>, <i>W</i> and <i>H</i> 
masses provided by ALPGEN are set in the PYTHIA 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>particle database</a>. 
   
 
<br/><br/><strong>Alpgen:setMLM</strong>  <input type="radio" name="4" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="4" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
When switched on, the merging parameters (see below) are set according to 
the ALPGEN hard process cuts: 
<ul> 
<li> <code>JetMatching:eTjetMin = min(ptjmin + 5., 1.2 * ptjmin)</code>, </li> 
<li> <code>JetMatching:coneRadius = drjmin</code>, 
<li> <code>JetMatching:etaJetMax = etajmax</code>. 
</ul> 
where the <code>ptjmin</code>, <code>drjmin</code> and 
<code>etajmax</code> are the incoming ALPGEN parameters. Note that any 
existing values of these parameters are overwritten. 
   
 
<br/><br/><strong>Alpgen:setNjet</strong>  <input type="radio" name="5" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="5" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
When switched on, the <code>JetMatching:nJet</code> parameter (see below) 
is set to the incoming <code>njet</code> ALPGEN parameter. Note that any 
existing value of this parameter is overwritten. 
   
 
<h2>Class information</h2> 
 
Some more technical information about the different classes is given 
below. For clarity, some limited information on certain private methods 
is provided. 
 
<h3>LHAupAlpgen</h3> 
 
This class is derived from the 
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?><code>LHAup</code></a> base class, and 
uses the standard machinery to pass initialisation and event data to 
PYTHIA. These standard functions are not documented here. The complete 
parameter file is stored in the PYTHIA Info object, if given, under the 
key <code>AlpgenPar</code>. 
 
<a name="method1"></a>
<p/><strong>LHAupAlpgen::LHAupAlpgen(const char *baseFNin, Info *infoPtrIn = NULL) &nbsp;</strong> <br/>
The constructor for the class takes the base filename for the ALPGEN 
format files (without file extensions) and optionally a pointer to a 
PYTHIA Info class, used for warning/error message printing and for 
storing the ALPGEN parameter file. The event and 
parameter files are opened immediately, with the <code>AlpgenPar</code> 
class, described below, used to parse the parameter file. 
   
 
<a name="method2"></a>
<p/><strong>bool LHAupAlpgen::addResonances() &nbsp;</strong> <br/>
This is a private method used when an event is read in. The information 
read from the event file does not always contain a complete listing of 
all particles and four-momenta, and so various details must be 
reconstructed. Exactly which details are filled in can vary based on the 
ALPGEN process in question. 
   
 
<a name="method3"></a>
<p/><strong>bool LHAupAlpgen::rescaleMomenta() &nbsp;</strong> <br/>
This is another private method used when an event is read in. 
It shuffles and rescales momenta in an event to ensure energy-momentum 
conservation.  First, <i>pT</i> is made to balance by splitting any 
imbalance between all outgoing particles with their energies also 
scaled. Second, the <i>e/pZ</i> of the two incoming particles are 
scaled to balance the outgoing particles. Finally, any intermediate 
resonances are recalculated from their decay products. 
   
 
<h3>AlpgenPar</h3> 
 
This class parses an ALPGEN parameter file and makes the information 
available through a simple interface. The information is stored 
internally in key/value (string/double) format. All lines prior to: 
<pre>  ************** run parameters </pre> 
are ignored, and in the general case, a line e.g. 
<pre>  10   3.00000000000000        ! njets</pre> 
would be stored with key "njets" and value "3.0". The following lines 
are special cases where the line may be split or the key translated: 
<pre> 
  3 ! hard process code 
  0.000   4.700 174.300  80.419  91.188 120.000 ! mc,mb,mt,mw,mz,mh 
  912.905 0.0914176   ! Crosssection +- error (pb) 
  100 29787.4  ! unwtd events, lum (pb-1) Njob= 2 
</pre> 
In the first line, the key "hard process code" is translated to 
"hpc". In the second, the mass values are split and each given an entry 
in the internal store. In the third, the cross section and cross section 
error are stored under the keys "xsecup" and "xerrup" respectively. 
Finally, the number of events and luminosity are stored under the keys 
"nevent" and "lum" respectively. In the event that a duplicate key is 
present, with differing values, the stored value is overwritten and a 
warning given. 
 
<a name="method4"></a>
<p/><strong>AlpgenPar::AlpgenPar(Info *infoPtrIn = NULL) &nbsp;</strong> <br/>
The constructor does nothing except for store the PYTHIA Info 
pointer, if given. This is used for warning/error message printing. 
   
 
<a name="method5"></a>
<p/><strong>bool AlpgenPar::parse(const string paramStr) &nbsp;</strong> <br/>
This method parses an ALPGEN parameter file. The parameter file is 
passed as a single string, mainly intended to be read out from the 
PYTHIA Info object using the header information methods. 
   
 
<a name="method6"></a>
<p/><strong>bool AlpgenPar::haveParam(const string &amp;paramIn) &nbsp;</strong> <br/>
Method to check if a parameter with key <code>paramIn</code> is present. 
Returns true if present, else false. 
   
 
<a name="method7"></a>
<p/><strong>double AlpgenPar::getParam(const string &amp;paramIn) &nbsp;</strong> <br/>
   
<strong>int AlpgenPar::getParamAsInt(const string &amp;paramIn) &nbsp;</strong> <br/>
Return the parameter with key <code>paramIn</code> as a double or 
integer. The presence of a parameter should have already been checked 
using the <code>haveParam()</code> function above. If the parameter is 
not present, 0 is returned. 
   
 
<a name="method8"></a>
<p/><strong>void AlpgenPar::void printParams() &nbsp;</strong> <br/>
Method to print a list of stored parameters. 
   
 
<h3>AlpgenHooks</h3> 
 
This <?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?><code>UserHooks</code></a> derived class 
provides all the <code>Alpgen:*</code> options. It is provided as a 
UserHooks class such that the code works regardless of whether ALPGEN 
native or LHE file formats are used. It is declared with virtual 
inheritance so that it may be combine with other UserHooks classes, see 
the "Combining UserHooks" section below. 
 
<a name="method9"></a>
<p/><strong>AlpgenHooks(Pythia &amp;pythia) &nbsp;</strong> <br/>
The constructor takes a PYTHIA object as input, so that the beam 
parameter settings can be overridden if the <code>Alpgen:file</code> 
option is given. If this is the case, an <code>LHAupAlpgen</code> 
instance is automatically created and passed to PYTHIA. 
   
 
<a name="method10"></a>
<p/><strong>bool initAfterBeams() &nbsp;</strong> <br/>
This is the only UserHooks method that is overridden. It is called 
directly after PYTHIA has initialised the beams, and therefore the 
header information should be present in the PYTHIA Info object. The 
<code>AlpgenPar</code> class is used to parse ALPGEN parameters, if 
present, which are then used to set further PYTHIA settings. 
   
 
<input type="hidden" name="saved" value="1"/>

<?php
echo "<input type='hidden' name='filepath' value='".$_GET["filepath"]."'/>"?>

<table width="100%"><tr><td align="right"><input type="submit" value="Save Settings" /></td></tr></table>
</form>

<?php

if($_POST["saved"] == 1)
{
$filepath = $_POST["filepath"];
$handle = fopen($filepath, 'a');

if($_POST["1"] != "void")
{
$data = "Alpgen:file = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "Alpgen:setLightMasses = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "Alpgen:setHeavyMasses = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "on")
{
$data = "Alpgen:setMLM = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "on")
{
$data = "Alpgen:setNjet = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
