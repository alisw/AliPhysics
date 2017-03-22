<html>
<head>
<title>Program Classes</title>
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

<form method='post' action='ProgramClasses.php'>

<h2>Program Classes</h2>

The complete PYTHIA 8 package contains 
<?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFiles.php?filepath=".$filepath."' target='page'>";?>a multitude of classes</a>.
There is no reason to describe all of them, with all of their
methods, since most should not be touched by a normal user.
Nevertheless some of the crucial ones are described in detail,
as a help not only to advanced users but also to developers. 
We here provide a quick reference (still incomplete) which 
classes you can find described where on these pages. Normally
you have to scroll down to find the details, since the top of
the page contains information of more general interest.

<p/>
<table cellspacing="5">

<tr> 
<td><b>Class</b></td> 
<td><b>Reference</b></td> 
<td><b>Comment</b></td> 
</tr>

<tr> 
<td><code>BeamShape</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='BeamShape.php?filepath=".$filepath."' target='page'>";?>Beam Shape</a></td> 
<td>vertex and beam momentum spread</td> 
</tr>

<tr> 
<td><code>CellJet</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='EventAnalysis.php?filepath=".$filepath."' target='page'>";?>Event Analysis</a></td> 
<td>jet cone clustering analysis, intended for hadron collider 
topologies</td> 
</tr>

<tr> 
<td><code>ClusterJet</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='EventAnalysis.php?filepath=".$filepath."' target='page'>";?>Event Analysis</a></td> 
<td>jet clustering analysis, intended for <i>e^+e^-</i> 
collider topologies</td> 
</tr>

<tr> 
<td><code>CoupSM</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='StandardModelParameters.php?filepath=".$filepath."' target='page'>";?>Standard Model Parameters</a></td> 
<td>couplings and mixing matrices</td> 
</tr>

<tr> 
<td><code>CoupSUSY</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='SUSYLesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>SUSY Les Houches Accord</a></td> 
<td>compute and store SUSY couplings, using input from SusyLesHouches</td> 
</tr>

<tr> 
<td><code>DecayChannel</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>Particle Data Scheme</a></td> 
<td>the properties of a single decay channel of  particle species</td> 
</tr>

<tr> 
<td><code>Event</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='EventRecord.php?filepath=".$filepath."' target='page'>";?>Event Record</a></td> 
<td>the complete event record</td> 
</tr>

<tr> 
<td><code>HepMC::I_Pythia8</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='HepMCInterface.php?filepath=".$filepath."' target='page'>";?>HepMC Interface</a></td> 
<td>convert a PYTHIA event record to the HepMC format</td> 
</tr>

<tr> 
<td><code>Hist</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='Histograms.php?filepath=".$filepath."' target='page'>";?>Histograms</a></td> 
<td>a primitive built-in histogramming package</td> 
</tr>

<tr> 
<td><code>Info</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>Event Information</a></td> 
<td>various one-of-a-kind information on the current event</td> 
</tr>

<tr> 
<td><code>LHAPDFInterface</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='PartonDistributions.php?filepath=".$filepath."' target='page'>";?>Parton Distributions</a></td> 
<td>use PDF sets from the LHAPDF library</td> 
</tr>

<tr> 
<td><code>LHAup</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches Accord</a></td> 
<td>base class for Les Houches external hard-event input</td> 
</tr>

<tr> 
<td><code>Particle</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>Particle Properties</a></td> 
<td>the properties of a particle in the event record</td> 
</tr>

<tr> 
<td><code>ParticleData</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>Particle Data Scheme</a></td> 
<td>the database of particle species properties</td> 
</tr>

<tr> 
<td><code>ParticleDataEntry</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>Particle Data Scheme</a></td> 
<td>the properties of a particle species</td> 
</tr>

<tr> 
<td><code>PDF</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='PartonDistributions.php?filepath=".$filepath."' target='page'>";?>Parton Distributions</a></td> 
<td>base class for all PDFs</td> 
</tr>

<tr> 
<td><code>Pythia</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?>Program Flow</a></td> 
<td>the top-level class, that drives the generation process</td> 
</tr>

<tr> 
<td><code>ResonanceWidths</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='SemiInternalResonances.php?filepath=".$filepath."' target='page'>";?>Semi-Internal Resonances</a></td> 
<td>base class for all resonance-width treatment</td> 
</tr>

<tr> 
<td><code>Rndm</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='RandomNumbers.php?filepath=".$filepath."' target='page'>";?>RandomNumbers</a></td> 
<td>a primitive built-in histogramming package</td> 
</tr>

<tr> 
<td><code>RotBstMatrix</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='FourVectors.php?filepath=".$filepath."' target='page'>";?>Four-Vectors</a></td> 
<td>rotation and boosts of four-vectors</td> 
</tr>

<tr> 
<td><code>Settings</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='SettingsScheme.php?filepath=".$filepath."' target='page'>";?>Settings Scheme</a></td> 
<td>the database that regulates the behaviour of the program</td> 
</tr>

<tr> 
<td><code>SigmaProcess</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='SemiInternalProcesses.php?filepath=".$filepath."' target='page'>";?>Semi-Internal Processes</a></td> 
<td>base class for all hard-process cross sections</td> 
</tr>

<tr> 
<td><code>SlowJet</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='EventAnalysis.php?filepath=".$filepath."' target='page'>";?>Event Analysis</a></td> 
<td>simple jet clustering using the <i>kT</i>, anti-<i>kT</i> or
Cambridge/Aachen algorithms.</td> 
</tr>

<tr> 
<td><code>SpaceShower</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='ImplementNewShowers.php?filepath=".$filepath."' target='page'>";?>Implement New Showers</a></td> 
<td>spacelike, i..e. initial-state showers</td> 
</tr>

<tr> 
<td><code>Sphericity</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='EventAnalysis.php?filepath=".$filepath."' target='page'>";?>Event Analysis</a></td> 
<td>sphericity analysis of events</td> 
</tr>

<tr> 
<td><code>SusyLesHouches</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='SUSYLesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>SUSY Les Houches Accord</a></td> 
<td>reading in of SLHA information</td> 
</tr>

<tr> 
<td><code>Thrust</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='EventAnalysis.php?filepath=".$filepath."' target='page'>";?>Event Analysis</a></td> 
<td>thrust analysis of events</td> 
</tr>

<tr> 
<td><code>TimeShower</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='ImplementNewShowers.php?filepath=".$filepath."' target='page'>";?>Implement New Showers</a></td> 
<td>timelike, i..e. final-state showers</td> 
</tr>

<tr> 
<td><code>UserHooks</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>User Hooks</a></td> 
<td>get access to some steps of the generation process</td> 
</tr>

<tr> 
<td><code>Vec4</code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='FourVectors.php?filepath=".$filepath."' target='page'>";?>Four-Vectors</a></td> 
<td>four-vectors</td> 
</tr>

<tr> 
<td><code></code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='.php?filepath=".$filepath."' target='page'>";?></a></td> 
<td></td> 
</tr>

<tr> 
<td><code></code></td> 
<td><?php $filepath = $_GET["filepath"];
echo "<a href='.php?filepath=".$filepath."' target='page'>";?></a></td> 
<td></td> 
</tr>

</table>

<p/>
You can also access an extensive list of the public 
<?php $filepath = $_GET["filepath"];
echo "<a href='ProgramMethods.php?filepath=".$filepath."' target='page'>";?>methods</a> in these classes, 
which compensates for the brevity of the current page. 

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
