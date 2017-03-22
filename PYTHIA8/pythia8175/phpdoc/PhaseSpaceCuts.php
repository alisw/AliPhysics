<html>
<head>
<title>Phase Space Cuts</title>
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

<form method='post' action='PhaseSpaceCuts.php'>

<h2>Phase Space Cuts</h2>

<code>PhaseSpace</code> is base class for all hard-process phase-space 
generators, either generic <i>2 -> 1</i> or <i>2 -> 2</i> ones, 
or specialized ones like for elastic and diffractive scattering.

<p/>
In it, it is possible to constrain the kinematics of most processes.
(Exceptions are "soft physics", i.e. minimum bias, elastic and 
diffractive processes. The Coulomb singularity for elastic scatterings,
if simulated, is <?php $filepath = $_GET["filepath"];
echo "<a href='TotalCrossSections.php?filepath=".$filepath."' target='page'>";?>handled separately</a>.) 
These constraints apply in the rest frame of the hard subprocess, and 
topologies normally would be changed e.g. by subsequent showering 
activity. The cross section of a process is adjusted to only 
correspond to the allowed phase space.

<p/>
The more particles in the final state, the more cuts could be applied.
Here we have tried to remain with the useful minimum, however. More
generic possibilities could be handled by the 
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>user hooks</a> facility. 

<h3>Cuts in all processes</h3>

<br/><br/><table><tr><td><strong>PhaseSpace:mHatMin </td><td></td><td> <input type="text" name="1" value="4." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The minimum invariant mass.
  

<br/><br/><table><tr><td><strong>PhaseSpace:mHatMax </td><td></td><td> <input type="text" name="2" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
The maximum invariant mass.
A value below <code>mHatMin</code> means there is no upper limit.
  

<h3>Cuts in <i>2 -> 1</i> processes</h3>

When a resonance <code>id</code> is produced, the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>mMin(id)</a></code> and 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>mMax(id)</a></code> 
methods restrict the allowed mass range
of this resonance. Therefore the allowed range is chosen to be the 
overlap of this range and the <code>mHatMin</code> to 
<code>mHatMax</code> range above. Most resonances by default have no 
upper mass limit, so effects mainly concern the lower limit. 
Should there be no overlap between the two ranges then the process 
will be switched off.

<h3>Cuts in <i>2 -> 2</i> processes</h3>

<br/><br/><table><tr><td><strong>PhaseSpace:pTHatMin </td><td></td><td> <input type="text" name="3" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The minimum invariant <i>pT</i>.
  

<br/><br/><table><tr><td><strong>PhaseSpace:pTHatMax </td><td></td><td> <input type="text" name="4" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
The maximum invariant <i>pT</i>.
A value below <code>pTHatMin</code> means there is no upper limit.
  

<br/><br/><table><tr><td><strong>PhaseSpace:pTHatMinDiverge </td><td></td><td> <input type="text" name="5" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.5</code>)</td></tr></table>
Extra <i>pT</i> cut to avoid the divergences of some processes 
in the limit <i>pT -> 0</i>. Specifically, if either or both
produced particles have a mass below <code>pTHatMinDiverge</code> 
then <i>pT</i> is limited from below by the larger of 
<code>pTHatMin</code> and <code>pTHatMinDiverge</code>.
  

<br/><br/><strong>PhaseSpace:useBreitWigners</strong>  <input type="radio" name="6" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="6" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allows masses to be selected according to Breit-Wigner shapes in 
<i>2 -> 2</i> processes, whenever particles have been declared 
with a nonvanishing width above the threshold below. In those cases 
also the limits below will be used for the mass selection. For 
<i>2 -> 1</i> processes the Breit-Wigner shape is part of the 
cross section itself, and therefore always included.
  

<br/><br/><table><tr><td><strong>PhaseSpace:minWidthBreitWigners </td><td></td><td> <input type="text" name="7" value="0.01" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01</strong></code>; <code>minimum = 1e-6</code>)</td></tr></table>
The minimum width a resonance must have for the mass to be dynamically
selected according to a Breit-Wigner shape, within the limits set below.
Only applies when <code>useBreitWigners</code> is on; else the nominal
mass value is always used.
  

<p/>
For a particle with a Breit-Wigner shape selected, according to the 
rules above and to the rules of the particle species itself, the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>mMin(id)</a></code> and 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>mMax(id)</a></code> 
methods restrict the allowed mass range of the particle, just like for 
the <i>2 -> 1 </i> processes.   

<h3>Cuts in <i>2 -> 3</i> processes</h3>

There are two main classes of <i>2 -> 3</i> processes. One is the 
processes such as <i>WW/ZZ</i>-fusion Higgs production, i.e.
<i>q q -> q q H</i>, where there are no special singularities 
associated with two partons in the final state being collinear,
or even for <i>pT -> 0</i>. For this class, no further cuts 
have been introduced than those already available for <i>2 -> 2</i> 
processes. Specifically, for now all three are restricted exactly the 
same way by <code>pTHatMin</code> and <code>pTHatMax</code>. As above, 
Breit-Wigner mass ranges can be restricted.

<p/>
The other <i>2 -> 3</i> event class is QCD processes, such as 
<i>g g -> g g g</i>. Here the soft and collinear singularities 
play a major role, and the phase space generation and cuts have 
been adapted to this. For this class, an alternative set of cuts 
is used, as outlined in the following. First of all the three
outgoing partons are ordered in falling <i>pT</i>, i.e. 
<i>pT_3 > pT_4 > pT_5</i> (where the labeling 3, 4, 5 of the outgoing 
partons is random, i.e. unrelated to the order specified in the
process name). The allowed ranges of <i>pT_3</i> and <i>pT_5</i>
can be specified, but obviously <i>pT_3max >= pT_5max</i> and
<i>pT_3min >= pT_5min</i>. The <i>pT_4</i> is not constrained 
explicitly, but is constructed from the vector sum of <i>pT_3</i>
and <i>pT_5</i>, subject to the constraint that it has to lie
between the two in magnitude. While the <i>pT</i> cuts take care
of singularities collinear with the incoming beams, it is also 
necessary to handle final-state singularities, when two outgoing
partons become collinear. This is done by requiring a minimal 
separation in <i>R</i>, where 
<i>R^2 = (Delta eta)^2 + (Delta phi)^2</i>. 
Finally, a note about efficiency. The QCD <i>2 -> 3</i> phase space 
is not set up to explicitly include <i>mHat</i> as one of the basic
variables. Such a cut is only done after a phase space point is already 
selected, which means that a narrow mass choice will slow down the 
program appreciably. Also narrow <i>pT_3</i> and <i>pT_5</i> bins
are likely to give inefficient generation, if it gives rise to
significant indirect restrictions on <i>pT_4</i>. 

<br/><br/><table><tr><td><strong>PhaseSpace:pTHat3Min </td><td></td><td> <input type="text" name="8" value="10." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The minimum invariant <i>pT</i> of the highest-<i>pT</i> parton in
QCD <i>2 -> 3</i> processes.
  

<br/><br/><table><tr><td><strong>PhaseSpace:pTHat3Max </td><td></td><td> <input type="text" name="9" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
The maximum invariant <i>pT</i> of the highest-<i>pT</i> parton in
QCD <i>2 -> 3</i> processes
A value below <code>pTHat3Min</code> means there is no upper limit.
  

<br/><br/><table><tr><td><strong>PhaseSpace:pTHat5Min </td><td></td><td> <input type="text" name="10" value="10." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The minimum invariant <i>pT</i> of the lowest-<i>pT</i> parton in
QCD <i>2 -> 3</i> processes.
  

<br/><br/><table><tr><td><strong>PhaseSpace:pTHat5Max </td><td></td><td> <input type="text" name="11" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
The maximum invariant <i>pT</i> of the lowest-<i>pT</i> parton in
QCD <i>2 -> 3</i> processes
A value below <code>pTHat5Min</code> means there is no upper limit.
  

<br/><br/><table><tr><td><strong>PhaseSpace:RsepMin </td><td></td><td> <input type="text" name="12" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The minimum separation <i>R</i> in <i>(eta, phi)</i> space between
any two outgoing partons in QCD <i>2 -> 3</i> processes.
  


<h3>Cuts for a second hard process</h3>

If you use the machinery that allows the generation of a specified 
<?php $filepath = $_GET["filepath"];
echo "<a href='ASecondHardProcess.php?filepath=".$filepath."' target='page'>";?>second hard process</a> then,
by default, the same phase space cuts will be used for it as listed
above. Optionally, however, you may use a second set of cuts, as 
described here. In this context "first" and "second" is merely a 
technical distinction; you are welcome e.g. to pick <i>pT</i> ranges 
such that the second interaction always has a larger <i>pT</i> than 
the first.

<br/><br/><strong>PhaseSpace:sameForSecond</strong>  <input type="radio" name="13" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="13" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
By default use the same cuts for a second hard process as for the 
first. If <code>off</code> then instead use the mass and <i>pT</i>
cuts below, where relevant. (The other cuts above still remain the same.)  
  

<br/><br/><table><tr><td><strong>PhaseSpace:mHatMinSecond </td><td></td><td> <input type="text" name="14" value="4." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The minimum invariant mass for a second interaction, if separate.
  

<br/><br/><table><tr><td><strong>PhaseSpace:mHatMaxSecond </td><td></td><td> <input type="text" name="15" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
The maximum invariant mass for a second interaction, if separate.
A value below <code>mHatMin</code> means there is no upper limit.
  

<br/><br/><table><tr><td><strong>PhaseSpace:pTHatMinSecond </td><td></td><td> <input type="text" name="16" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The minimum invariant <i>pT</i> for a second interaction, if separate.
  

<br/><br/><table><tr><td><strong>PhaseSpace:pTHatMaxSecond </td><td></td><td> <input type="text" name="17" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
The maximum invariant <i>pT</i> for a second interaction, if separate.
A value below <code>pTHatMin</code> means there is no upper limit.
  

<h3>Generation strategy and documentation</h3>

During the initialization stage a simplified function is found,
that is intended to be above the true cross-section behaviour
over the whole of phase space. It is chosen to be easily integrable 
and invertible. That way a trial phase space point can be selected 
according this simple function, and then be accepted by the ratio of
true to the simple function. For a good efficiency the ratio should be
close to unity,  yet never above it. This constrains the absolute 
normalization of the simple function. The initial search may fail to 
find the phase space point where the true-to-simple ratio is maximal,
however. This then can lead to subsequent maximum violations, where the 
ratio is above unity. Two alternative strategies are implemented to 
handle such situations, see below.

<br/><br/><strong>PhaseSpace:showSearch</strong>  <input type="radio" name="18" value="on"><strong>On</strong>
<input type="radio" name="18" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Possibility to print information on the search for phase-space 
coefficients that (in a multichannel approach) provides an analytical 
upper envelope of the differential cross section, and the 
corresponding upper estimate of the cross section. Of interest 
for crosschecks by expert users only. 
  

<br/><br/><strong>PhaseSpace:showViolation</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Possibility to print information whenever the assumed maximum 
differential cross section of a process is violated, i.e. when 
the initial maximization procedure did not find the true maximum.
Also, should negative cross sections occur, print whenever a more
negative value is encountered.
  

<br/><br/><strong>PhaseSpace:increaseMaximum</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Strategy for handling cases where a larger cross section is
obtained during the event generation than was assumed at initialization,
i.e. when a violation occurs.
<br/><b>off:</b>each event comes with a weight, which normally is unity
(as a consequence of the acceptance/rejection step), and is found in
<code><?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>Info::weight()</a></code>. 
For events which exceed the maximum instead the true-to-simple ratio 
is stored as event weight, which then is above unity. If the user so 
wishes this weight can then be carried along when event properties are
histogrammed. Since normally such violations should be rare and not 
too much above unity one could expect most users to ignore such issues
be default. Should maximum violations turn out to be frequent (visible
in the <code><?php $filepath = $_GET["filepath"];
echo "<a href='EventStatistics.php?filepath=".$filepath."' target='page'>";?>Pythia::statistics()</a></code>
output) the option exists to use the information.
<br/><b>on:</b>the maximum is increased whenever it is exceeded. Thus
events generated after this point will be "correctly" distributed,
while ones generated previously obviously then have had too high a 
relative weight. If violations occur early on and/or are small this
strategy should do a good job of correcting to the desired phase-space
distribution. This strategy may be more convenient for the normal user,
who would not wish to worry about event weights. It does have the
disadvantage that the raised maximum introduces an extra amount of
"history memory" to the generation sequence, so that it becomes less
easy to save-and-restore the <?php $filepath = $_GET["filepath"];
echo "<a href='RandomNumbers.php?filepath=".$filepath."' target='page'>";?>random-number 
state</a> for debugging purposes.  
  

<h3>Reweighting of <i>2 -> 2</i> processes</h3>

Events normally come with unit weight, i.e. are distributed across
the allowed phase space region according to the appropriate differential
cross sections. Sometimes it may be convenient to have an uneven 
distribution of events. The classical example here is that many cross
sections drop off with transverse momentum <i>pT</i>, such that few
events are generated at large <i>pT</i> scales. If one wants to 
plot the <i>pT</i> cross section, and all that comes with it, the 
statistical error will then degrade with increasing <i>pT</i> 
where fewer events end up. 

<p/>
One solution is to split the full <i>pT</i> range into several 
separate subranges, where the events of each subsample obtains a 
different overall normalization. Specifically, if you generate a
comparable number of events in each <i>pT</i> bin, such that 
larger <i>pT</i> bins are oversampled, these bins come with a
correspondingly reduced overall weight, that needs to be taken into
account when the bins are combined. The other is to have a continuously
increasing oversampling of events at larger <i>pT</i> scales, which
is compensated by a continuously decreasing weight for the event.

<p/>
Both of these solutions are supported. Specifically, for 
<i>2 -> 2</i> processes, the <i>pTHat</i> scale offers a 
convenient classification of the event. (Of course, two events 
starting out from the same <i>pTHat</i> scale will experience
different parton shower evolutions, etc., and may therefore look 
quite different at the end.) The two cuts 
<code>PhaseSpace:pTHatMin</code> and <code>PhaseSpace:pTHatMax</code>
therefore offers a way to slice a <i>pT</i> range into subranges,
see e.g. <code>main08.cc</code>. Alternatively the 
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>User Hooks</a> machinery offers the 
possibility for you to define your own reweighting of phase space
sampling, with a corresponding event weight, with 
<code>UserHooks::canBiasSelection</code> and related methods.

<p/>
As a simplified option, we here offer the possibility to bias the 
<i>2 -> 2</i> sampling by a power of <i>pTHat</i>, then with     
events having a weight the inverse of this. This fast track will only
work under a number of strict conditions, implemented to reduce the 
risk of abuse. (Whereas a <code>UserHooks</code> setup can be more 
flexible.) Specifically it will work if only high-<i>pT</i>
<i>2 -> 2</i> processes already implemented in PYTHIA are requested,
notably the <code>HardQCD</code> ones. That is, you cannot mix with 
<i>2 -> 1</i> or <i>2 -> 3</i> processes, nor with external 
processes (notably Les Houches input) or <code>SoftQCD</code> ones, 
and  you cannot use the option to define a 
<?php $filepath = $_GET["filepath"];
echo "<a href='ASecondHardProcess.php?filepath=".$filepath."' target='page'>";?>second hard process</a> in
the same event. Furthermore you have to be careful about the choice 
of <code>PhaseSpace:pTHatMin</code>, since a <i>pTHat = 0</i> 
event would come with an infinite weight. 

<br/><br/><strong>PhaseSpace:bias2Selection</strong>  <input type="radio" name="21" value="on"><strong>On</strong>
<input type="radio" name="21" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Possibility to switch on a biased phase space sampling, 
with compensatingly weighted events, for <i>2 -> 2</i> processes. 
Can only be used under the specific conditions explained in 
the paragraph above; under other conditions the initialization 
will abort. 
  

<br/><br/><table><tr><td><strong>PhaseSpace:bias2SelectionPow </td><td></td><td> <input type="text" name="22" value="4." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
If the above flag is on, then a <i>2 -> 2</i> process at a scale 
<i>pTHat</i> will be oversampled in phase space by an amount
<i>(pTHat/pTRef)^pow</i>, where you set the power <i>pow</i>
here. Events are assigned a compensating 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>weight</a> the inverse of this,
i.e. <code>Info::weight()</code> will return <i>(pTRef/pTHat)^pow</i>.
This weight should then be used in the histogramming of event properties.
The final overall normalization also involves the 
<code>Info::weightSum()</code> value.  
  

<br/><br/><table><tr><td><strong>PhaseSpace:bias2SelectionRef </td><td></td><td> <input type="text" name="23" value="10." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10.</strong></code>; <code>minimum = 1.</code>)</td></tr></table>
The reference scale <i>pTRef</i> introduced above, such that events
with this <i>pTHat</i> obtain unit weight in the reweighting procedure.
The value of this parameter has no impact on the final result of the
reweighting procedure, but is only there for convenience, i.e. to
give "reasonably-sized" weights.  
  

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

if($_POST["1"] != "4.")
{
$data = "PhaseSpace:mHatMin = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "-1.")
{
$data = "PhaseSpace:mHatMax = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "0.")
{
$data = "PhaseSpace:pTHatMin = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "-1.")
{
$data = "PhaseSpace:pTHatMax = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1.")
{
$data = "PhaseSpace:pTHatMinDiverge = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "on")
{
$data = "PhaseSpace:useBreitWigners = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "0.01")
{
$data = "PhaseSpace:minWidthBreitWigners = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "10.")
{
$data = "PhaseSpace:pTHat3Min = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "-1.")
{
$data = "PhaseSpace:pTHat3Max = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "10.")
{
$data = "PhaseSpace:pTHat5Min = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "-1.")
{
$data = "PhaseSpace:pTHat5Max = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "1.")
{
$data = "PhaseSpace:RsepMin = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "on")
{
$data = "PhaseSpace:sameForSecond = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "4.")
{
$data = "PhaseSpace:mHatMinSecond = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "-1.")
{
$data = "PhaseSpace:mHatMaxSecond = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "0.")
{
$data = "PhaseSpace:pTHatMinSecond = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "-1.")
{
$data = "PhaseSpace:pTHatMaxSecond = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "off")
{
$data = "PhaseSpace:showSearch = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "PhaseSpace:showViolation = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "PhaseSpace:increaseMaximum = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "off")
{
$data = "PhaseSpace:bias2Selection = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "4.")
{
$data = "PhaseSpace:bias2SelectionPow = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "10.")
{
$data = "PhaseSpace:bias2SelectionRef = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
