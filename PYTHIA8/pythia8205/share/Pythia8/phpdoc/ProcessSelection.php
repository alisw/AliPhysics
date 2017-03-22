<html>
<head>
<title>Process Selection</title>
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

<form method='post' action='ProcessSelection.php'>
 
<h2>Process Selection</h2> 
 
There is no way PYTHIA could contain all processes of interest, 
neither in terms of potential physics topics nor in terms of 
high-multiplicity final states. What exists is a reasonably 
complete setup of all <i>2 &rarr; 1</i> and <i>2 &rarr; 2</i> 
processes within the Standard Model, plus some examples of 
processes beyond that, again for low multiplicities. Combined with 
the PYTHIA parton showers, this should be enough to get a flying 
start in the study of many physics scenarios. 
Other processes could be fed in via the 
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches Accord</a> 
or be implemented as a 
<?php $filepath = $_GET["filepath"];
echo "<a href='SemiInternalProcesses.php?filepath=".$filepath."' target='page'>";?>Semi-Internal Process</a>. 
In the latter case the existing processes would act as obvious 
templates. 
 
<p/> 
By default all processes are switched off. You should switch on 
those you want to simulate. This may be done at two (occasionally 
three) levels, either for each individual process or for a group of 
processes. That is, a process is going to be generated either if its 
own flag or its group flag is on. There is no built-in construction 
to switch on a group and then switch off a few of its members. 
 
<p/> 
Each process is assigned an integer code. This code is not used in 
the internal administration of events (so having the same code for 
two completely different processes would not be a problem), but only 
intended to allow a simpler user separation of different processes. 
Also the process name is available, as a string. 
 
<p/> 
To ease navigation, the list of processes has been split into several 
separate pages, by main topic. The classification is hopefully 
intuitive, but by no means unambiguous. For instance, essentially 
all processes involve QCD, so the "QCD processes" are the ones that 
only involve QCD. (And also that is not completely true, once one 
includes all that may happen in multiparton interactions.) On these 
separate pages also appear the settings that are completely local 
to that particular process class, but not the ones that have a 
broader usage. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='QCDProcesses.php?filepath=".$filepath."' target='page'>";?>QCD Processes</a></h3> 
 
QCD processes fall in two main categories: soft and hard. The soft ones 
contain elastic, diffractive and "minimum-bias" events, together 
covering the total cross section. Hard processes are the normal 
<i>2 &rarr; 2</i> ones, including charm and bottom production. 
<br/>Reserved code range: 101 - 199. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='ElectroweakProcesses.php?filepath=".$filepath."' target='page'>";?>Electroweak Processes</a></h3> 
 
Prompt-photon, <i>gamma^*/Z^0</i> and <i>W^+-</i> production, 
plus a few processes with <i>t</i>-channel boson exchange. 
<br/>Reserved code range: 201 - 299. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='OniaProcesses.php?filepath=".$filepath."' target='page'>";?>Onia Processes</a></h3> 
 
Colour singlet and octet production of charmonium and bottomonium. 
<br/>Reserved code range: 401 - 499 for charmonium and 
501 - 599 for bottomonium. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='TopProcesses.php?filepath=".$filepath."' target='page'>";?>Top Processes</a></h3> 
 
Top production, singly or doubly. 
<br/>Reserved code range: 601 - 699. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='FourthGenerationProcesses.php?filepath=".$filepath."' target='page'>";?>Fourth-Generation 
Processes</a></h3> 
 
Production of hypothetical fourth-generation fermions. 
<br/>Reserved code range: 801 - 899. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='HiggsProcesses.php?filepath=".$filepath."' target='page'>";?>Higgs Processes</a></h3> 
 
Higgs production, within or beyond the Standard Model. 
See section on Left-Right-Symmetry processes for doubly charged Higgs bosons. 
<br/>Reserved code range: 901 - 999 for a Standard Model Higgs 
and 1001 - 1199 for MSSM Higgs bosons. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='SUSYProcesses.php?filepath=".$filepath."' target='page'>";?>SUSY Processes</a></h3> 
 
Production of supersymmetric particles, currently barely begun. 
<br/>Reserved code range: 1001 - 2999. (Whereof 1001 - 1199 
for Higgs bosons; see above.) 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='NewGaugeBosonProcesses.php?filepath=".$filepath."' target='page'>";?>New-Gauge-Boson 
Processes</a></h3> 
 
Production of new gauge bosons such as <i>Z'</i> and <i>W'</i>. 
<br/>Reserved code range: 3001 - 3099. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='LeftRightSymmetryProcesses.php?filepath=".$filepath."' target='page'>";?>Left-Right-Symmetry 
Processes</a></h3> 
 
Production of righthanded <i>Z_R</i> and <i>W_R</i> bosons and of 
doubly charged Higgs bosons. 
<br/>Reserved code range: 3101 - 3199. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='LeptoquarkProcesses.php?filepath=".$filepath."' target='page'>";?>Leptoquark Processes</a></h3> 
 
Production of a simple scalar leptoquark state. 
<br/>Reserved code range: 3201 - 3299. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='CompositenessProcesses.php?filepath=".$filepath."' target='page'>";?>Compositeness Processes</a></h3> 
 
Production of excited fermion states and contact-interaction modification 
to interactions between fermions (excluding technicolor; see below). 
<br/>Reserved code range: 4001 - 4099. 
 
<h3>Technicolor Processes</h3> 
 
Production of technicolor particles and modifications of QCD processes 
by technicolor interactions. Does not exist yet. 
<br/>Reserved code range: 4101 - 4199. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='HiddenValleyProcesses.php?filepath=".$filepath."' target='page'>";?>Hidden Valley Processes</a></h3> 
A scenario for the pair production of new particles with couplings 
under a new gauge group, with invisible gauge bosons. Radiation of 
these gauge bosons is included in the standard final-state parton 
shower. 
<br/>Reserved code range: 4901 - 4999. 
 
<h3><?php $filepath = $_GET["filepath"];
echo "<a href='ExtraDimensionalProcesses.php?filepath=".$filepath."' target='page'>";?>Extra-Dimensional 
Processes</a></h3> 
 
A vast area, here represented by the production of a Randall-Sundrum 
excited graviton state and a Kaluza-Klein gluon, a Kaluza-Klein tower 
of <i>gamma/Z^0</i> excitations in one TeV^- sized extra dimension, 
several Large Extra Dimension processes, and a few related Unparticle 
processes. 
<br/>Reserved code range: 5001 - 5099. 
 
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
 
