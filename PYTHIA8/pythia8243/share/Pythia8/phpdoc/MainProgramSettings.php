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
<ol id="toc">
  <li><a href="#section0">Introduction</a></li>
  <li><a href="#section1">Initialization settings</a></li>
  <li><a href="#section2">Event-generation settings</a></li>
  <li><a href="#section3">Statistics</a></li>
  <li><a href="#section4">Main-program settings</a></li>
  <li><a href="#section5">Subruns</a></li>
  <li><a href="#section6">Spares</a></li>
</ol>

 
<a name="section0"></a> 
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
However, to make life simpler, a few useful main-program 
settings are defined on this page, so that they are recognized by 
the <code>Settings</code> machinery. They can thus be put among 
the other cards without distinction. It is up to you to decide which 
ones, if any, you actually want to use when you write your main program. 
 
<p/> 
To further reduce the necessary amount of main-program code, some of 
the tasks that you can steer from your main program can also be done 
internally. This in particular relates to some information printing. 
To give an example, <code>pythia.event.list()</code> can be inserted 
to print an event, i.e. all the particles belonging to it. Given the 
length of these listings one would list only a few events at the 
beginning of the run, to get some feeling for the character of events. 
This could be achieved e.g. with a main-program statement<br/> 
<code>   if (iEvent &lt; 3) pythia.event.list()</code><br/> 
to list the first three events in a loop over <code>iEvent</code>, 
after <code>pythia.next()</code> has been used to generate the next 
event. Alternatively a <code>Next:numberShowEvent = 3</code> 
setting, e.g. in a command file, would achieve the same, by an 
internal call at the end of <code>pythia.next()</code>. 
 
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
 
<a name="section1"></a> 
<h3>Initialization settings</h3> 
 
<br/><br/><strong>Init:showProcesses</strong>  <input type="radio" name="1" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="1" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print a list of all processes that will be simulated, with 
their estimated cross section maxima, as used for the 
subsequent Monte Carlo selection. Also print corresponding 
Les Houches initialization data, where relevant. 
   
 
<br/><br/><strong>Init:showMultipartonInteractions</strong>  <input type="radio" name="2" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="2" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print initialization information for the multiparton interactions 
machinery. 
   
 
<br/><br/><strong>Init:showChangedSettings</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print a list of the changed flag/mode/parameter/word settings. 
   
 
<br/><br/><strong>Init:showAllSettings</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print a list of all flag/mode/parameter/word settings. 
Warning: this will be a long list. 
   
 
<br/><br/><strong>Init:showChangedParticleData</strong>  <input type="radio" name="5" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="5" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print a list of particle and decay data for those particles 
that were changed (one way or another). 
   
 
<br/><br/><strong>Init:showChangedResonanceData</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
In the previous listing also include the resonances that are 
initialized at the beginning of a run and thus get new particle 
data, even if these may well agree with the default ones. 
Warning: this will be a rather long list. 
   
 
<br/><br/><strong>Init:showAllParticleData</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print a list of all particle and decay data. 
Warning: this will be a long list. 
   
 
<br/><br/><table><tr><td><strong>Init:showOneParticleData  </td><td></td><td> <input type="text" name="8" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Print particle and decay data for the particle with this particular 
identity code. Default means that no particle is printed. 
   
 
<a name="section2"></a> 
<h3>Event-generation settings</h3> 
 
<br/><br/><table><tr><td><strong>Next:numberCount  </td><td></td><td> <input type="text" name="9" value="1000" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1000</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Print a line telling how many events have been generated so far, once 
every <code>numberCount</code> events. If set zero then no lines are 
ever printed. 
<br/>In <code>include/Pythia8Plugins/ProgressLog.h</code> an alternative 
method is implemented that intermittently prints out run progress 
information, reports on CPU usage and estimates when the run will end. 
It is used in the <code>main111.cc</code> example. 
   
 
<br/><br/><table><tr><td><strong>Next:numberShowLHA  </td><td></td><td> <input type="text" name="10" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to list the Les Houches input information for, 
where relevant. 
   
 
<br/><br/><table><tr><td><strong>Next:numberShowInfo  </td><td></td><td> <input type="text" name="11" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to list the <code>Info</code> information for, 
where relevant. 
   
 
<br/><br/><table><tr><td><strong>Next:numberShowProcess  </td><td></td><td> <input type="text" name="12" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to list the <code>process</code> record for, 
where relevant. 
   
 
<br/><br/><table><tr><td><strong>Next:numberShowEvent  </td><td></td><td> <input type="text" name="13" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to list the <code>event</code> record for, 
where relevant. 
   
 
<br/><br/><strong>Next:showScaleAndVertex</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
In addition to the normal information in the listing of the 
<code>process</code> and <code>event</code> records, a second line 
per particle provides information on the production scale, 
particle polarization and production vertex. 
   
 
<br/><br/><strong>Next:showMothersAndDaughters</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
In addition to the normal information in the listing of the 
<code>process</code> and <code>event</code> records, further lines 
list all the mothers and daughters of each particle. 
   
 
<a name="section3"></a> 
<h3>Statistics</h3> 
 
<br/><br/><strong>Stat:showProcessLevel</strong>  <input type="radio" name="16" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="16" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print the available statistics on number of generated events and 
cross sections, where relevant. 
   
 
<br/><br/><strong>Stat:showPartonLevel</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print the available statistics on number and types of multiparton 
interactions, where relevant. 
   
 
<br/><br/><strong>Stat:showErrors</strong>  <input type="radio" name="18" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="18" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print the available statistics on number and types of 
aborts, errors and warnings. 
   
 
<br/><br/><strong>Stat:reset</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Reset the statistics of the above three kinds. The default is that 
all stored statistics information is unaffected by the 
<code>pythia.stat()</code> call. Counters are automatically reset 
in each new <code>pythia.init()</code> call, however, so the only time 
the reset option makes a difference is if <code>stat()</code> 
is called several times in a (sub)run. 
   
 
<a name="section4"></a> 
<h3>Main-program settings</h3> 
 
The settings in this section <i>must</i> be under the control of the 
user, i.e. there are no internal equivalents. The first one is especially 
important and would be a standard feature of any separate command file. 
 
<br/><br/><table><tr><td><strong>Main:numberOfEvents  </td><td></td><td> <input type="text" name="20" value="1000" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1000</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to be generated. 
   
 
<br/><br/><table><tr><td><strong>Main:numberOfTriedEvents  </td><td></td><td> <input type="text" name="21" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to be tried during generation. Any number smaller than one 
means that the setting will be ignored. 
   
 
<br/><br/><table><tr><td><strong>Main:numberOfSelectedEvents  </td><td></td><td> <input type="text" name="22" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to be selected during generation. Any number smaller than 
one means that the setting will be ignored. 
   
 
<br/><br/><table><tr><td><strong>Main:numberOfAcceptedEvents  </td><td></td><td> <input type="text" name="23" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to be accepted during generation. Any number smaller than 
one means that the setting will be ignored. 
   
 
<br/><br/><table><tr><td><strong>Main:timesAllowErrors  </td><td></td><td> <input type="text" name="24" value="10" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10</strong></code>)</td></tr></table>
Allow this many times that <code>pythia.next()</code> returns false, 
i.e. that an event is flawed, before aborting the run. 
   
 
<p/> 
The <code>Main:...</code> options works like this. Once you have used 
the <code>pythia.readFile(fileName)</code> method to read in the cards 
file, where the values have been set, you can interrogate the 
<code>Settings</code> database to make the values available in your 
main program. A slight complication is that you need to use a different 
<code>Settings</code> method for each of the four possible return types 
that you want to extract. To save some typing the same method names are 
found directly in the <code>Pythia</code> class, and these just send on 
to the <code>Settings</code> ones to do the job, e.g.<br/> 
<code>   int nEvent = pythia.mode("Main:numberOfEvents"); </code> 
 
<p/> 
The area of subruns is covered separately below. A few spares are also 
defined after that, for unforeseen applications. 
 
<a name="section5"></a> 
<h3>Subruns</h3> 
 
You can use <?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?>subruns</a> to carry out 
several tasks in the same run. In that case you will need repeated 
instances of the first setting below in your command file, and could 
additionally use the second and third as well. 
 
<br/><br/><table><tr><td><strong>Main:subrun  </td><td></td><td> <input type="text" name="25" value="-999" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-999</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of the current subrun, a non-negative integer, put as 
first line in a section of lines to be read for this particular subrun. 
   
 
<br/><br/><strong>Main:LHEFskipInit</strong>  <input type="radio" name="26" value="on"><strong>On</strong>
<input type="radio" name="26" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If you read several Les Houches Event Files that you want to see 
considered as one single combined event sample you can set this flag 
<code>on</code> after the first subrun to skip (most of) the 
(re-)initialization step. 
   
 
<br/><br/><table><tr><td><strong>Main:numberOfSubruns  </td><td></td><td> <input type="text" name="27" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
The number of subruns you intend to use in your current run. 
Unlike the two settings above, <code>Pythia</code> itself will not 
interpret this number, but you could e.g. have a loop in your main 
program to loop over subruns from 0 through 
<code>numberOfSubruns - 1</code>. 
   
 
<a name="section6"></a> 
<h3>Spares</h3> 
 
For currently unforeseen purposes, a few dummy settings are made 
available here. The user can set the desired value in a "cards file" 
and then use that value in the main program as desired. 
 
<br/><br/><strong>Main:spareFlag1</strong>  <input type="radio" name="28" value="on"><strong>On</strong>
<input type="radio" name="28" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
   
 
<br/><br/><strong>Main:spareFlag2</strong>  <input type="radio" name="29" value="on"><strong>On</strong>
<input type="radio" name="29" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
   
 
<br/><br/><strong>Main:spareFlag3</strong>  <input type="radio" name="30" value="on"><strong>On</strong>
<input type="radio" name="30" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
   
 
<br/><br/><table><tr><td><strong>Main:spareMode1  </td><td></td><td> <input type="text" name="31" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
   
 
<br/><br/><table><tr><td><strong>Main:spareMode2  </td><td></td><td> <input type="text" name="32" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
   
 
<br/><br/><table><tr><td><strong>Main:spareMode3  </td><td></td><td> <input type="text" name="33" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
   
 
<br/><br/><table><tr><td><strong>Main:spareParm1 </td><td></td><td> <input type="text" name="34" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
   
 
<br/><br/><table><tr><td><strong>Main:spareParm2 </td><td></td><td> <input type="text" name="35" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
   
 
<br/><br/><table><tr><td><strong>Main:spareParm3 </td><td></td><td> <input type="text" name="36" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
   
 
<br/><br/><table><tr><td><strong>Main:spareWord1  </td><td></td><td> <input type="text" name="37" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
   
 
<br/><br/><table><tr><td><strong>Main:spareWord2  </td><td></td><td> <input type="text" name="38" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
   
 
<br/><br/><table><tr><td><strong>Main:spareWord3  </td><td></td><td> <input type="text" name="39" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
   
 
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
$data = "Init:showProcesses = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "on")
{
$data = "Init:showMultipartonInteractions = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "Init:showChangedSettings = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "Init:showAllSettings = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "on")
{
$data = "Init:showChangedParticleData = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "Init:showChangedResonanceData = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "Init:showAllParticleData = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "0")
{
$data = "Init:showOneParticleData = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "1000")
{
$data = "Next:numberCount = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "1")
{
$data = "Next:numberShowLHA = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "1")
{
$data = "Next:numberShowInfo = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "1")
{
$data = "Next:numberShowProcess = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "1")
{
$data = "Next:numberShowEvent = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "Next:showScaleAndVertex = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "Next:showMothersAndDaughters = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "on")
{
$data = "Stat:showProcessLevel = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "Stat:showPartonLevel = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "on")
{
$data = "Stat:showErrors = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "Stat:reset = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "1000")
{
$data = "Main:numberOfEvents = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "0")
{
$data = "Main:numberOfTriedEvents = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "0")
{
$data = "Main:numberOfSelectedEvents = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "0")
{
$data = "Main:numberOfAcceptedEvents = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "10")
{
$data = "Main:timesAllowErrors = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "-999")
{
$data = "Main:subrun = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "off")
{
$data = "Main:LHEFskipInit = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "0")
{
$data = "Main:numberOfSubruns = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "off")
{
$data = "Main:spareFlag1 = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "off")
{
$data = "Main:spareFlag2 = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "off")
{
$data = "Main:spareFlag3 = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "0")
{
$data = "Main:spareMode1 = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "0")
{
$data = "Main:spareMode2 = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "0")
{
$data = "Main:spareMode3 = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "0.")
{
$data = "Main:spareParm1 = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "0.")
{
$data = "Main:spareParm2 = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "0.")
{
$data = "Main:spareParm3 = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "void")
{
$data = "Main:spareWord1 = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "void")
{
$data = "Main:spareWord2 = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "void")
{
$data = "Main:spareWord3 = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
