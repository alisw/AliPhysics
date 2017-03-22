<html>
<head>
<title>Master Switches</title>
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

<form method='post' action='MasterSwitches.php'>

<h2>Master Switches</h2>

Sometimes it may be convenient to omit certain aspects of the event 
generation chain. This cannot be motivated in a full-blown production
run, but can often be convenient for own understanding and for
debug purposes. The flags on this page allow just that.

<p/>
The event generation is subdivided into three levels: the process
level, the parton level and the hadron level, and flags are grouped
accordingly. 

<h3>Process Level</h3>

The <code>ProcessLevel</code> class administrates the initial step of 
the event generation, wherein the basic process is selected. Currently 
this is done either using some of the internal processes, or with 
Les Houches Accord input.

<p/>
There could not be a complete event without an initial process, so
it would not be a normal action to switch off this step. Furthermore,
without a process set, it is also not possible to carry out the tasks
on the parton level. It is still possible, however, to hadronize 
a parton-level configuration provided by some external program.

<br/><br/><strong>ProcessLevel:all</strong>  <input type="radio" name="1" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="1" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If off, do not attempt to carry out any generation at all on the 
process level. For the parton level only final-state radiation
is possible, using the <code>Pythia::forceTimeShower(...)</code> method.
Do allow parton configurations stored in the event record to hadronize 
and hadrons to decay, however, as set by the <code>HadronLevel</code> 
switches. Further details are found 
<?php $filepath = $_GET["filepath"];
echo "<a href='HadronLevelStandalone.php?filepath=".$filepath."' target='page'>";?>here</a>.
   

<p/>
For <code>ProcessLevel:all = on</code> one part of the event generation 
on this level may be switched off individually: 

<br/><br/><strong>ProcessLevel:resonanceDecays</strong>  <input type="radio" name="2" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="2" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch to allow resonance decays; on/off = true/false.
Normal hadrons and leptons do not count as resonances, so this is 
aimed specifically towards <i>Z^0, W^+-, t, h^0</i> and similar
objects beyond the Standard Model. Do not use this option if you
may produce coloured resonances and intend to allow hadronization,
since currently the program would not know how to handle this.
  

<p/>
It is possible to stop the generation immediately after the basic 
process has been selected, see <code>PartonLevel:all</code> below.

<h3>PartonLevel</h3>

The <code>PartonLevel</code> class administrates the middle step of the 
event generation, i.e. the evolution from an input (hard) process from
<code>ProcessLevel</code>, containing a few partons only, to a complete 
parton-level configuration to be handed on to <code>HadronLevel</code>. 
This step involves the application of initial- and final-state radiation, 
multiparton interactions and the structure of beam remnants.

<br/><br/><strong>PartonLevel:all</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If off then stop the generation after the hard process has been 
generated, but before the parton-level and hadron-level steps. 
The <code>process</code> record is filled, but the <code>event</code> 
one is then not.
  

<p/>
For <code>PartonLevel:all = on</code> some parts of the event generation 
on this level may be switched off individually: 

<br/><br/><strong>PartonLevel:MPI</strong>  <input type="radio" name="4" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="4" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for multiparton interactions; on/off = true/false.
Further options are found <?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?>here</a>.
  

<br/><br/><strong>PartonLevel:ISR</strong>  <input type="radio" name="5" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="5" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for initial-state radiation; on/off = true/false.
Further options are found <?php $filepath = $_GET["filepath"];
echo "<a href='SpacelikeShowers.php?filepath=".$filepath."' target='page'>";?>here</a>.
  

<br/><br/><strong>PartonLevel:FSR</strong>  <input type="radio" name="6" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="6" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for final-state radiation; on/off = true/false.
Further options are found <?php $filepath = $_GET["filepath"];
echo "<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>";?>here</a>.
If you leave this switch on, the following two switches allow 
more detailed control to switch off only parts of the showers. 
  

<br/><br/><strong>PartonLevel:FSRinProcess</strong>  <input type="radio" name="7" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="7" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Switch for final-state radiation in association with the hard process 
itself; on/off = true/false. In addition <code>PartonLevel:FSR</code>
must be on for these emissions to occur. 
  

<br/><br/><strong>PartonLevel:FSRinResonances</strong>  <input type="radio" name="8" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="8" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for final-state radiation in any resonance decays 
subsequent to the hard process itself; on/off = true/false. In addition 
<code>PartonLevel:FSR</code> must be on for these emissions to occur.
  

<p/>
Switching off all the above MPI/ISR/FSR switches is <b>not</b> equivalent 
to setting <code>PartonLevel:all = off</code>. In the former case a 
minimal skeleton of parton-level operations are carried out, such as 
tying together the scattered partons with the beam remnants into colour 
singlets, and storing this information in the <code>event</code> record. 
It is therefore possible to go on and hadronize the event, if desired. 
In the latter case <b>no</b> operations at all are carried out on the 
parton level, and therefore it is also not possible to go on to the 
hadron level.

<br/><br/><strong>PartonLevel:Remnants</strong>  <input type="radio" name="9" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="9" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for addition of beam remnants; on/off = true/false.  
Only intended for very special applications, and cannot be used to 
generate complete events. Specifically, unlike the other switches above,
the program will complain and possibly crash unlike you also set 
<code>HadronLevel:all = off</code> and <code>Check:event = off</code>.
  

<p/>
It is possible to stop the generation immediately after the parton level 
has been set up, see <code>HadronLevel:all</code> below.

<h3>HadronLevel</h3>

The <code>HadronLevel</code> class administrates the final step of the 
event generation, wherein the partonic configuration from 
<code>PartonLevel</code> is hadronized, including string fragmentation 
and secondary decays.

<p/>
Most of the code in this class itself deals with subdividing the partonic
content of the event into separate colour singlets, that can be
treated individually by the string fragmentation machinery. When a
junction and an antijunction are directly connected, it also breaks 
the string between the two, so that the topology can be reduced back 
to two separate one-junction systems, while still preserving the
expected particle flow in the junction-junction string region(s).

<br/><br/><strong>HadronLevel:all</strong>  <input type="radio" name="10" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="10" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If off then stop the generation after the hard process and 
parton-level activity has been generated, but before the 
hadron-level steps.
  

<p/>
For <code>HadronLevel:all = on</code> some parts of the event generation 
on this level may be switched off individually: 

<br/><br/><strong>HadronLevel:Hadronize</strong>  <input type="radio" name="11" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="11" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for hadronization; on/off = true/false.
Further options are found <?php $filepath = $_GET["filepath"];
echo "<a href='Fragmentation.php?filepath=".$filepath."' target='page'>";?>here</a>.
  

<br/><br/><strong>HadronLevel:Decay</strong>  <input type="radio" name="12" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="12" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for decays; on/off = true/false.
Further options are found <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDecays.php?filepath=".$filepath."' target='page'>";?>here</a>.
  

<br/><br/><strong>HadronLevel:BoseEinstein</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Master switch for the simulation of Bose-Einstein effects; 
on/off = true/false. Further options are found 
<?php $filepath = $_GET["filepath"];
echo "<a href='BoseEinsteinEffects.php?filepath=".$filepath."' target='page'>";?>here</a>.
  

<h3>Printing</h3>

<br/><br/><strong>Print:quiet</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Can be set on to avoid the printing during program execution, to the 
largest extent possible. This flag acts by setting the relevant values 
for <code>Init:showProcesses</code>, 
<code>Init:showMultipartonInteractions</code>,  
<code>Init:showChangedSettings</code>,  
<code>Init:showAllSettings</code>,  
<code>Init:showChangedParticleData</code>,  
<code>Init:showChangedResonanceData</code>,  
<code>Init:showAllParticleData</code>,  
<code>Init:showOneParticleData</code>,  
<code>Next:numberCount</code>,  
<code>Next:numberShowLHA</code>,  
<code>Next:numberShowInfo</code>,  
<code>Next:numberShowProcess</code>, and  
<code>Next:numberShowEvent</code>. 
The change is to off or 0 for <code>Print:quiet = off</code>,
and restores to the respective default value for <code>= on</code>. 
Those changes take effect immediately, so individual settings can be 
changed afterwards.   
  


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

if($_POST["1"] != "on")
{
$data = "ProcessLevel:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "on")
{
$data = "ProcessLevel:resonanceDecays = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "PartonLevel:all = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "on")
{
$data = "PartonLevel:MPI = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "on")
{
$data = "PartonLevel:ISR = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "on")
{
$data = "PartonLevel:FSR = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "on")
{
$data = "PartonLevel:FSRinProcess = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "on")
{
$data = "PartonLevel:FSRinResonances = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "on")
{
$data = "PartonLevel:Remnants = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "on")
{
$data = "HadronLevel:all = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "on")
{
$data = "HadronLevel:Hadronize = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "on")
{
$data = "HadronLevel:Decay = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "HadronLevel:BoseEinstein = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "Print:quiet = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
