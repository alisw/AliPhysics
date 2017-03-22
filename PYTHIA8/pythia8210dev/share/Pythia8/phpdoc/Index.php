<html>
<head>
<title>Index</title>
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

<form method='post' action='Index.php'>
 
<img src='pythia99.gif' alt='Pythia logo' hspace=10/> 
 
<h2>PYTHIA 8 Index</h2> 
 
<h4>Program Overview</h4> 
 
<a href='Frontpage.php?filepath=".$filepath."' target='page'>Frontpage</a><br/> 
<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>Program Flow</a><br/> 
<a href='SettingsScheme.php?filepath=".$filepath."' target='page'>Settings Scheme</a><br/> 
<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>Particle Data Scheme</a><br/> 
<a href='ProgramFiles.php?filepath=".$filepath."' target='page'>Program Files</a><br/> 
<a href='ProgramClasses.php?filepath=".$filepath."' target='page'>Program Classes</a><br/> 
<a href='ProgramMethods.php?filepath=".$filepath."' target='page'>Program Methods</a><br/> 
<a href='SampleMainPrograms.php?filepath=".$filepath."' target='page'>Sample Main Programs</a><br/> 
 
<h4>Setup Run Tasks</h4> 
 
<?php
$filepath = "files/".$_GET["filename"];
$filename = $_GET["filename"];
echo "<a href='SaveSettings.php?returning=1&filename=".$filename."' target='page'><b>Save Settings</b></a><br/> 
<a href='MainProgramSettings.php?filepath=".$filepath."' target='page'>Main-Program Settings</a><br/> 
<a href='BeamParameters.php?filepath=".$filepath."' target='page'>Beam Parameters</a><br/> 
<a href='RandomNumberSeed.php?filepath=".$filepath."' target='page'>Random-Number Seed</a><br/> 
<a href='PDFSelection.php?filepath=".$filepath."' target='page'>PDF Selection</a><br/> 
<a href='MasterSwitches.php?filepath=".$filepath."' target='page'>Master Switches</a><br/> 
<a href='ProcessSelection.php?filepath=".$filepath."' target='page'>Process Selection</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='QCDProcesses.php?filepath=".$filepath."' target='page'>QCD</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='ElectroweakProcesses.php?filepath=".$filepath."' target='page'>Electroweak</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='OniaProcesses.php?filepath=".$filepath."' target='page'>Onia</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='TopProcesses.php?filepath=".$filepath."' target='page'>Top</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='FourthGenerationProcesses.php?filepath=".$filepath."' target='page'>Fourth Generation</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='HiggsProcesses.php?filepath=".$filepath."' target='page'>Higgs</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='SUSYProcesses.php?filepath=".$filepath."' target='page'>SUSY</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='NewGaugeBosonProcesses.php?filepath=".$filepath."' target='page'>New Gauge Bosons</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='LeftRightSymmetryProcesses.php?filepath=".$filepath."' target='page'>Left-Right Symmetry</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='LeptoquarkProcesses.php?filepath=".$filepath."' target='page'>Leptoquark</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='CompositenessProcesses.php?filepath=".$filepath."' target='page'>Compositeness</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='HiddenValleyProcesses.php?filepath=".$filepath."' target='page'>Hidden Valleys</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='ExtraDimensionalProcesses.php?filepath=".$filepath."' target='page'>Extra Dimensions</a><br/> 
<a href='ASecondHardProcess.php?filepath=".$filepath."' target='page'>A Second Hard Process</a><br/> 
<a href='PhaseSpaceCuts.php?filepath=".$filepath."' target='page'>Phase Space Cuts</a><br/> 
<a href='CouplingsAndScales.php?filepath=".$filepath."' target='page'>Couplings and Scales</a><br/> 
<a href='StandardModelParameters.php?filepath=".$filepath."' target='page'>Standard-Model Parameters</a><br/> 
<a href='TotalCrossSections.php?filepath=".$filepath."' target='page'>Total Cross Sections</a><br/> 
<a href='ResonanceDecays.php?filepath=".$filepath."' target='page'>Resonance Decays</a><br/> 
<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>Timelike Showers</a><br/> 
<a href='SpacelikeShowers.php?filepath=".$filepath."' target='page'>Spacelike Showers</a><br/> 
<a href='WeakShowers.php?filepath=".$filepath."' target='page'>Weak Showers</a><br/> 
<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>Multiparton Interactions</a><br/> 
<a href='BeamRemnants.php?filepath=".$filepath."' target='page'>Beam Remnants</a><br/> 
<a href='ColourReconnection.php?filepath=".$filepath."' target='page'>Colour Reconnection</a><br/> 
<a href='Diffraction.php?filepath=".$filepath."' target='page'>Diffraction</a><br/> 
<a href='Fragmentation.php?filepath=".$filepath."' target='page'>Fragmentation</a><br/> 
<a href='FlavourSelection.php?filepath=".$filepath."' target='page'>Flavour Selection</a><br/> 
<a href='ParticleDecays.php?filepath=".$filepath."' target='page'>Particle Decays</a><br/> 
<a href='RHadrons.php?filepath=".$filepath."' target='page'>R-hadrons</a><br/> 
<a href='BoseEinsteinEffects.php?filepath=".$filepath."' target='page'>Bose-Einstein Effects</a><br/> 
<a href='HadronScattering.php?filepath=".$filepath."' target='page'>Hadron Scattering</a><br/> 
<a href='ParticleData.php?filepath=".$filepath."' target='page'>Particle Data</a><br/> 
<a href='ErrorChecks.php?filepath=".$filepath."' target='page'>Error Checks</a><br/> 
<a href='Tunes.php?filepath=".$filepath."' target='page'>Tunes</a><br/> 

";?>
 
<h4>Study Output</h4> 
 
<?php
$filepath = "files/".$_GET["filename"];
$filename = $_GET["filename"];
echo "<a href='FourVectors.php?filepath=".$filepath."' target='page'>Four-Vectors</a><br/> 
<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>Particle Properties</a><br/> 
<a href='EventRecord.php?filepath=".$filepath."' target='page'>Event Record</a><br/> 
<a href='EventInformation.php?filepath=".$filepath."' target='page'>Event Information</a><br/> 
<a href='EventStatistics.php?filepath=".$filepath."' target='page'>Event Statistics</a><br/> 
<a href='EventAnalysis.php?filepath=".$filepath."' target='page'>Event Analysis</a><br/> 
<a href='Histograms.php?filepath=".$filepath."' target='page'>Histograms</a><br/> 
<a href='AdvancedUsage.php?filepath=".$filepath."' target='page'>Advanced Usage</a><br/> 

";?>
 
<h4>Link to Other Programs</h4> 
 
<?php
$filepath = "files/".$_GET["filename"];
$filename = $_GET["filename"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>Les Houches Accord</a><br/> 
<a href='SUSYLesHouchesAccord.php?filepath=".$filepath."' target='page'>SUSY Les Houches Accord</a><br/> 
<a href='HepMCInterface.php?filepath=".$filepath."' target='page'>HepMC Interface</a><br/> 
<a href='ProMCFiles.php?filepath=".$filepath."' target='page'>ProMC Files</a><br/> 
<a href='SemiInternalProcesses.php?filepath=".$filepath."' target='page'>Semi-Internal Processes</a><br/> 
<a href='SemiInternalResonances.php?filepath=".$filepath."' target='page'>Semi-Internal Resonances</a><br/> 
<a href='MadGraph5Processes.php?filepath=".$filepath."' target='page'>MadGraph 5 Processes</a><br/> 
<a href='AlpgenEventInterface.php?filepath=".$filepath."' target='page'>Alpgen Event Interface</a><br/> 
<a href='MatchingAndMerging.php?filepath=".$filepath."' target='page'>Matching and Merging</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='POWHEGMerging.php?filepath=".$filepath."' target='page'>POWHEG Merging</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='aMCatNLOMatching.php?filepath=".$filepath."' target='page'>aMC@NLO Matching</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>CKKW-L Merging</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='JetMatching.php?filepath=".$filepath."' target='page'>Jet Matching</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='UMEPSMerging.php?filepath=".$filepath."' target='page'>UMEPS Merging</a><br/> 
&nbsp;&nbsp;--&nbsp;&nbsp; 
<a href='NLOMerging.php?filepath=".$filepath."' target='page'>NLO Merging</a><br/> 
<a href='UserHooks.php?filepath=".$filepath."' target='page'>User Hooks</a><br/> 
<a href='HadronLevelStandalone.php?filepath=".$filepath."' target='page'>Hadron-Level Standalone</a><br/> 
<a href='ExternalDecays.php?filepath=".$filepath."' target='page'>External Decays</a><br/> 
<a href='BeamShape.php?filepath=".$filepath."' target='page'>Beam Shape</a><br/> 
<a href='PartonDistributions.php?filepath=".$filepath."' target='page'>Parton Distributions</a><br/> 
<a href='JetFinders.php?filepath=".$filepath."' target='page'>Jet Finders</a><br/> 
<a href='RandomNumbers.php?filepath=".$filepath."' target='page'>Random Numbers</a><br/> 
<a href='ImplementNewShowers.php?filepath=".$filepath."' target='page'>Implement New Showers</a><br/> 
<a href='RIVETusage.php?filepath=".$filepath."' target='page'>RIVET usage</a><br/> 
<a href='ROOTusage.php?filepath=".$filepath."' target='page'>ROOT usage</a><br/> 

";?>
 
<h4>Reference Materiel</h4> 
 
<a href='UpdateHistory.php?filepath=".$filepath."' target='page'>Update History</a><br/> 
<a href='Bibliography.php?filepath=".$filepath."' target='page'>Bibliography</a><br/> 
<a href='Glossary.php?filepath=".$filepath."' target='page'>Glossary</a><br/> 
<a href='Version.php?filepath=".$filepath."' target='page'>Version</a> 
 
<h4>Separate documents</h4> 
 
<a href="../pdfdoc/pythia8200.pdf" target="page">Introduction</a><br/> 
<a href="../pdfdoc/worksheet8200.pdf" target="page">Worksheet</a> 
 
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
