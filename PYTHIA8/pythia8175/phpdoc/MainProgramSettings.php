<html>
<head>
<title>Main-Program and Related Settings</title>
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

<form method='post' action='MainProgramSettings.php'>

<h2>Main-Program and Related Settings</h2>

<h3>Introduction</h3>

The main program is up to the user to write. However, 
<?php $filepath = $_GET["filepath"];
echo "<a href='SampleMainPrograms.php?filepath=".$filepath."' target='page'>";?>sample main programs</a> 
are provided. In one such class of programs, key settings of the run 
are read in from a "cards file". For experimental collaborations
this is actually the most common way to run a program like PYTHIA.
The commands in such a file may be of two types<br/>
(a) instructions directly to <code>Pythia</code>, like which 
processes to generate, and<br/>
(b) instructions to the main program for what it should do, 
like how many events to generate, and how many events should 
be listed.<br/>
In principle these two kinds could be kept completely separate. 
However, to make life simpler, a number of useful main-program 
settings are defined on this page, so that they are recognized by 
the <code>Settings</code> machinery. They can thus be put among 
the other cards without distinction. It is up to you to decide which 
ones, if any, you actually want to use when you write your main program.

<p/>
To further reduce the necessary amount of main-program code, some of 
the tasks that you can steer from your main program can also be done 
internally. This in particular relates to some information printing.
To give an example, the <code>Main:numberToList</code> mode can be 
used by you in your main program to decide to list a few of the
generated events, whereas <code>Next:numberListEvent</code> is used 
internally in a <code>pythia.next()</code> call to do such a listing
automatically. Ultimately, in both cases, a 
<code>pythia.event.list()</code> call is the one that generates
the listing, explicitly in the main program in the former case,
implicitly called from <code>pythia.next()</code> in the latter.  

<p/>
The settings names on this page thus fall into four main groups
<ul>
<li><code>Init:...</code> denote actions that automatically may be
taken during the <code>pythia.init()</code> call.</li>
<li><code>Next:...</code> denote actions that automatically may be
taken during the <code>pythia.next()</code> call.</li>
<li><code>Stat:...</code> denote actions that automatically may be
taken during the <code>pythia.stat()</code> call.</li>
<li><code>Main:...</code> denote actions that you yourself 
have the freedom to make use of in your main program.</li>
</ul>
The use of several of the <code>Main:...</code> options is deprecated
in favour of the other possibilities.

<p/>
The <code>Main:...</code> options works like this. Once you have used 
the <code>pythia.readFile(fileName)</code> method to read in the cards 
file, where the values have been set, you can interrogate the 
<code>Settings</code> database to make the values available in your 
main program. A slight complication is that you need to use a different  
<code>Settings</code> method for each of the four possible return types 
that you want to extract. To save some typing the same method names are 
found directly in the <code>Pythia</code> class, and just send on to the
<code>Settings</code> ones to do the job, e.g.
<pre>
  bool   showCS = pythia.flag("Main:showChangedSettings");
  int    nEvent = pythia.mode("Main:numberOfEvents");
  double spare1 = pythia.parm("Main:spareParm1");
  string file   = pythia.word("Main:allSettingsFile"); 
</pre>

<h3>Main-program settings</h3>

The settings in this section <i>must</i> be under the control of the
user, i.e. there are no internal equivalents.

<br/><br/><table><tr><td><strong>Main:numberOfEvents  </td><td></td><td> <input type="text" name="1" value="1000" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1000</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to be generated.
  

<br/><br/><table><tr><td><strong>Main:timesAllowErrors  </td><td></td><td> <input type="text" name="2" value="10" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10</strong></code>)</td></tr></table>
Allow this many times that <code>pythia.next()</code> returns false, 
i.e. that an event is flawed, before aborting the run.
  

<h3>Initialization settings</h3>

<br/><br/><strong>Init:showProcesses</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print a list of all processes that will be simulated, with 
their estimated cross section maxima, as used for the 
subsequent Monte Carlo selection. Also print corresponding 
Les Houches initialization data, where relevant. 
  

<br/><br/><strong>Init:showMultipartonInteractions</strong>  <input type="radio" name="4" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="4" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print initialization information for the multiparton interactions 
machinery.
  

<br/><br/><strong>Init:showChangedSettings</strong>  <input type="radio" name="5" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="5" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print a list of the changed flag/mode/parameter/word settings.
  

<br/><br/><strong>Init:showAllSettings</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print a list of all flag/mode/parameter/word settings.
Warning: this will be a long list.
  

<br/><br/><strong>Init:showChangedParticleData</strong>  <input type="radio" name="7" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="7" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print a list of particle and decay data for those particles 
that were changed (one way or another).
  

<br/><br/><strong>Init:showChangedResonanceData</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
In the previous listing also include the resonances that are 
initialized at the beginning of a run and thus get new particle
data, even if these may well agree with the default ones. 
Warning: this will be a rather long list.
  

<br/><br/><strong>Init:showAllParticleData</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print a list of all particle and decay data.
Warning: this will be a long list.
  

<br/><br/><table><tr><td><strong>Init:showOneParticleData  </td><td></td><td> <input type="text" name="10" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Print particle and decay data for the particle with this particular 
identity code. Default means that no particle is printed.
  

<br/><br/><strong>Main:showChangedSettings</strong>  <input type="radio" name="11" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="11" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Deprecated. Print a list of the changed flag/mode/parameter/word settings.
  

<br/><br/><strong>Main:showAllSettings</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Deprecated. Print a list of all flag/mode/parameter/word settings.
Warning: this will be a long list.
  

<br/><br/><strong>Main:showChangedParticleData</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Deprecated. Print a list of particle and decay data for those particles 
that were changed (one way or another).
  

<br/><br/><strong>Main:showChangedResonanceData</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Deprecated. In the previous listing also include the resonances that are 
initialized at the beginning of a run and thus get new particle
data, even if these may well agree with the default ones. 
Warning: this will be a rather long list.
  

<br/><br/><strong>Main:showAllParticleData</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Deprecated. Print a list of all particle and decay data.
Warning: this will be a long list.
  

<br/><br/><table><tr><td><strong>Main:showOneParticleData  </td><td></td><td> <input type="text" name="16" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Deprecated. Print particle and decay data for the particle with this 
particular identity code. Default means that no particle is printed.
  

<br/><br/><strong>Main:writeChangedSettings</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Write a file with the changed flag/mode/parameter/word settings, in
a format appropriate to be read in at the beginning of a new  
run, using the <code>pythia.readFile(fileName)</code> method. 
  

<br/><br/><table><tr><td><strong>Main:changedSettingsFile  </td><td></td><td> <input type="text" name="18" value="currentSettings.cmnd" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>currentSettings.cmnd</strong></code>)</td></tr></table>
The name of the file to which the changed flag/mode/parameter/word
settings are written if <code>Main:writeChangedSettings</code>
is on. 
  

<br/><br/><strong>Main:writeAllSettings</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Write a file with all flag/mode/parameter/word settings, in
a format appropriate to be read in at the beginning of a new  
run, using the <code>pythia.readFile(fileName)</code> method. 
  

<br/><br/><table><tr><td><strong>Main:allSettingsFile  </td><td></td><td> <input type="text" name="20" value="allSettings.cmnd" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>allSettings.cmnd</strong></code>)</td></tr></table>
The name of the file to which a flag/mode/parameter/word 
settings are written if <code>Main:writeAllSettings</code>
is on. 
  

<h3>Event-generation settings</h3>

<br/><br/><table><tr><td><strong>Next:numberCount  </td><td></td><td> <input type="text" name="21" value="1000" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1000</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Print a line telling how many events have been generated so far,
once every <code>numberCount</code> events. If set zero then no
lines are ever printed. 

<br/><br/><table><tr><td><strong>Next:numberShowLHA  </td><td></td><td> <input type="text" name="22" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to list the Les Houches input information for,
where relevant.
  

<br/><br/><table><tr><td><strong>Next:numberShowInfo  </td><td></td><td> <input type="text" name="23" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to list the <code>Info</code> information for,
where relevant.
  

<br/><br/><table><tr><td><strong>Next:numberShowProcess  </td><td></td><td> <input type="text" name="24" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to list the <code>process</code> record for,
where relevant.
  

<br/><br/><table><tr><td><strong>Next:numberShowEvent  </td><td></td><td> <input type="text" name="25" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to list the <code>event</code> record for,
where relevant.
  

<br/><br/><strong>Next:showScaleAndVertex</strong>  <input type="radio" name="26" value="on"><strong>On</strong>
<input type="radio" name="26" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
In addition to the normal information in the listing of the 
<code>process</code> and <code>event</code> records, a second line
per particle provides information on the production scale,
particle polarization and production vertex. 
  

<br/><br/><strong>Next:showMothersAndDaughters</strong>  <input type="radio" name="27" value="on"><strong>On</strong>
<input type="radio" name="27" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
In addition to the normal information in the listing of the 
<code>process</code> and <code>event</code> records, further lines
list all the mothers and daughters of each particle. 
  

<br/><br/><table><tr><td><strong>Main:numberToList  </td><td></td><td> <input type="text" name="28" value="2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Deprecated. The number of events to list.
  

<br/><br/><table><tr><td><strong>Main:timesToShow  </td><td></td><td> <input type="text" name="29" value="50" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>50</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Deprecated. Print the number of events generated so far, this many times, 
i.e. once every <code>numberOfEvents/numberToShow</code> events.
  

<h3>Statistics</h3>

<br/><br/><strong>Stat:showProcessLevel</strong>  <input type="radio" name="30" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="30" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print the available statistics on number of generated events and
cross sections, where relevant.
  

<br/><br/><strong>Stat:showPartonLevel</strong>  <input type="radio" name="31" value="on"><strong>On</strong>
<input type="radio" name="31" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print the available statistics on number and types of multiparton
interactions, where relevant.
  

<br/><br/><strong>Stat:showErrors</strong>  <input type="radio" name="32" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="32" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print the available statistics on number and types of 
aborts, errors and warnings. 
  

<br/><br/><strong>Stat:reset</strong>  <input type="radio" name="33" value="on"><strong>On</strong>
<input type="radio" name="33" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Reset the statistics of the above three kinds. The default is that 
all stored statistics information is unaffected by the 
<code>pythia.stat()</code> call. Counters are automatically reset 
in each new <code>pythia.init()</code> call, however, so the only time 
the reset option makes a difference is if <code>stat()</code> 
is called several times in a (sub)run.  
  

<br/><br/><strong>Main:showAllStatistics</strong>  <input type="radio" name="34" value="on"><strong>On</strong>
<input type="radio" name="34" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print all available statistics or only the minimal set at the end 
of the run.
  

<h3>Subruns</h3>

You can use <?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?>subruns</a> to carry out
several tasks in the same run. In that case you will need repeated
instances of the first setting below in your command file, and could
additionally use the second and third as well.

<br/><br/><table><tr><td><strong>Main:subrun  </td><td></td><td> <input type="text" name="35" value="-999" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-999</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of the current subrun, a non-negative integer, put as
first line in a section of lines to be read for this particular subrun.
  

<br/><br/><strong>Main:LHEFskipInit</strong>  <input type="radio" name="36" value="on"><strong>On</strong>
<input type="radio" name="36" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If you read several Les Houches Event Files that you want to see 
considered as one single combined event sample you can set this flag
<code>on</code> after the first subrun to skip (most of) the  
(re-)initialization step.
  

<br/><br/><table><tr><td><strong>Main:numberOfSubruns  </td><td></td><td> <input type="text" name="37" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
The number of subruns you intend to use in your current run.  
Unlike the two settings above, <code>Pythia</code> itself will not
interpret this number, but you could e.g. have a loop in your main
program to loop over subruns from 0 through 
<code>numberOfSubruns - 1</code>. 
  

<h3>Spares</h3>

For currently unforeseen purposes, a few dummy settings are made 
available here. The user can set the desired value in a "cards file"
and then use that value in the main program as desired.

<br/><br/><strong>Main:spareFlag1</strong>  <input type="radio" name="38" value="on"><strong>On</strong>
<input type="radio" name="38" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
  

<br/><br/><strong>Main:spareFlag2</strong>  <input type="radio" name="39" value="on"><strong>On</strong>
<input type="radio" name="39" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
  

<br/><br/><strong>Main:spareFlag3</strong>  <input type="radio" name="40" value="on"><strong>On</strong>
<input type="radio" name="40" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
  

<br/><br/><table><tr><td><strong>Main:spareMode1  </td><td></td><td> <input type="text" name="41" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareMode2  </td><td></td><td> <input type="text" name="42" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareMode3  </td><td></td><td> <input type="text" name="43" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareParm1 </td><td></td><td> <input type="text" name="44" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareParm2 </td><td></td><td> <input type="text" name="45" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareParm3 </td><td></td><td> <input type="text" name="46" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareWord1  </td><td></td><td> <input type="text" name="47" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareWord2  </td><td></td><td> <input type="text" name="48" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareWord3  </td><td></td><td> <input type="text" name="49" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
  

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

if($_POST["1"] != "1000")
{
$data = "Main:numberOfEvents = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "10")
{
$data = "Main:timesAllowErrors = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "Init:showProcesses = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "on")
{
$data = "Init:showMultipartonInteractions = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "on")
{
$data = "Init:showChangedSettings = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "Init:showAllSettings = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "on")
{
$data = "Init:showChangedParticleData = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "Init:showChangedResonanceData = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "Init:showAllParticleData = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "0")
{
$data = "Init:showOneParticleData = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "on")
{
$data = "Main:showChangedSettings = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "Main:showAllSettings = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "Main:showChangedParticleData = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "Main:showChangedResonanceData = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "Main:showAllParticleData = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "0")
{
$data = "Main:showOneParticleData = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "Main:writeChangedSettings = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "currentSettings.cmnd")
{
$data = "Main:changedSettingsFile = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "Main:writeAllSettings = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "allSettings.cmnd")
{
$data = "Main:allSettingsFile = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "1000")
{
$data = "Next:numberCount = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "1")
{
$data = "Next:numberShowLHA = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "1")
{
$data = "Next:numberShowInfo = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "1")
{
$data = "Next:numberShowProcess = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "1")
{
$data = "Next:numberShowEvent = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "off")
{
$data = "Next:showScaleAndVertex = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "off")
{
$data = "Next:showMothersAndDaughters = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "2")
{
$data = "Main:numberToList = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "50")
{
$data = "Main:timesToShow = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "on")
{
$data = "Stat:showProcessLevel = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "off")
{
$data = "Stat:showPartonLevel = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "on")
{
$data = "Stat:showErrors = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "off")
{
$data = "Stat:reset = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "off")
{
$data = "Main:showAllStatistics = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "-999")
{
$data = "Main:subrun = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "off")
{
$data = "Main:LHEFskipInit = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "0")
{
$data = "Main:numberOfSubruns = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "off")
{
$data = "Main:spareFlag1 = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "off")
{
$data = "Main:spareFlag2 = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
if($_POST["40"] != "off")
{
$data = "Main:spareFlag3 = ".$_POST["40"]."\n";
fwrite($handle,$data);
}
if($_POST["41"] != "0")
{
$data = "Main:spareMode1 = ".$_POST["41"]."\n";
fwrite($handle,$data);
}
if($_POST["42"] != "0")
{
$data = "Main:spareMode2 = ".$_POST["42"]."\n";
fwrite($handle,$data);
}
if($_POST["43"] != "0")
{
$data = "Main:spareMode3 = ".$_POST["43"]."\n";
fwrite($handle,$data);
}
if($_POST["44"] != "0.")
{
$data = "Main:spareParm1 = ".$_POST["44"]."\n";
fwrite($handle,$data);
}
if($_POST["45"] != "0.")
{
$data = "Main:spareParm2 = ".$_POST["45"]."\n";
fwrite($handle,$data);
}
if($_POST["46"] != "0.")
{
$data = "Main:spareParm3 = ".$_POST["46"]."\n";
fwrite($handle,$data);
}
if($_POST["47"] != "void")
{
$data = "Main:spareWord1 = ".$_POST["47"]."\n";
fwrite($handle,$data);
}
if($_POST["48"] != "void")
{
$data = "Main:spareWord2 = ".$_POST["48"]."\n";
fwrite($handle,$data);
}
if($_POST["49"] != "void")
{
$data = "Main:spareWord3 = ".$_POST["49"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
