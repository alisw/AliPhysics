<html>
<head>
<title>Rope Hadronization</title>
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

<form method='post' action='RopeHadronization.php'>
 
<h2>Rope Hadronization</h2> 
<ol id="toc">
  <li><a href="#section0">Main settings</a></li>
  <li><a href="#section1">String shoving</a></li>
  <li><a href="#section2">Flavour Ropes</a></li>
</ol>

 
In collisions of protons, there are often tens of multiparton interactions, 
all producing Lund strings occupying the same area in transverse space of 
<i>~1</i> fm<i>^2</i>. The Rope Hadronization framework describes the 
interactions between such overlapping strings, by (a) allowing nearby strings 
to shove each other with an interaction potential derived from the colour 
superconductor analogy [<a href="Bibliography.php#refBie16b" target="page">Bie16b</a>], [<a href="Bibliography.php#refBie17" target="page">Bie17</a>] and (b) at 
hadronization time, colour charges at string endpoints and in gluon "kinks" 
can act together coherently to form a "rope", which is hadronized with a 
larger, effective string tension [<a href="Bibliography.php#refBie14" target="page">Bie14</a>]. 
The latter has noticeable effects on the flavour composition of the hadronic 
final state [<a href="Bibliography.php#refBie15" target="page">Bie15</a>], and this effect is denoted "flavour ropes" below. 
 
<p/> 
Since both models deal with string overlaps in transverse space, it is 
necessary to provide such information, as it is not present in the Pythia MPI 
model. The information is provided through the 
<?php $filepath = $_GET["filepath"];
echo "<a href='VertexInformation.php?filepath=".$filepath."' target='page'>";?>Parton Vertex</a> methods. The string 
shoving mechanism is exemplified in the <code>main101</code> example, 
and the flavour ropes in <code>main102</code>. 
A simpler version of flavour composition ropes exist [<a href="Bibliography.php#refBie16c" target="page">Bie16c</a>], 
which do not require vertex information. This can be enabled by a switch. 
 
<a name="section0"></a> 
<h3>Main settings</h3> 
 
The main settings are common for both the string shoving and the flavour rope 
models. 
 
<br/><br/><strong>Ropewalk:RopeHadronization</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Master switch for all aspects of rope hadronization. 
The Rope Hadronization framework is intended to work seamlessly with the rest 
of Pythia 8. It is, however, still a new model, and no Pythia tunes with ropes 
enabled exists yet. Therefore Rope Hadronization must be explicitly switched 
on, and it is up to the user to provide a sensible tune. 
   
 
<br/><br/><strong>Ropewalk:doShoving</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Enable the string shoving mechanism. In addition to this, the above 
<code>Ropewalk:RopeHadronization</code> flag must also be switched on. 
   
 
<br/><br/><strong>Ropewalk:doFlavour</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Enable the flavour ropes mechanism. In addition to this, the above 
<code>Ropewalk:RopeHadronization</code> flag must also be switched on. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:r0 </td><td></td><td> <input type="text" name="4" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
The transverse radius of a string, in units of fm. This can be viewed as an 
overall strength parameter of the Rope Hadronization framework, as all effects 
scale with increasing string overlap. Notice that the value for the string 
radius must be seen compared to the parameters which determines string 
placement in the transverse plane, as described in 
<?php $filepath = $_GET["filepath"];
echo "<a href='VetexInformation.php?filepath=".$filepath."' target='page'>";?>Vertex Information</a>. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:m0 </td><td></td><td> <input type="text" name="5" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.01</code>; <code>maximum = 5.0</code>)</td></tr></table>
Imposed lower mass cutoff to allow for calculation of rapidities of dipoles 
composed of massless gluons. 
   
 
<a name="section1"></a> 
<h3>String shoving</h3> 
 
The string shoving mechanism allows strings to push each other, before 
hadronization, as described in [<a href="Bibliography.php#refBie16b" target="page">Bie16b</a>]. 
 
<p/> 
String shoving divides the event up in many small rapidity slices (in 
the lab frame), and all string pieces in all slices are allowed to push 
each other with a force: 
 
<br/><i> 
f(d_\perp) = \frac{g_A \kappa d_\perp}{R^2} 
\exp\left(-\frac{d^2_\perp }{4R^2}\right), 
</i><br/> 
 
where <i>d_\perp</i> is the distance in transverse space between two string 
pieces, calculated dynamically using 
<?php $filepath = $_GET["filepath"];
echo "<a href='VetexInformation.php?filepath=".$filepath."' target='page'>";?>Vertex Information</a>. Model parameters are 
<i>g_A</i>, the amplitude of the shoving force, <i>R</i>, the string 
radius, and <i>g_E</i>, a parameter dividing the equilibrium string radius 
to account for the effect of strings starting out with a vanishing string 
radius. 
 
<p/> 
The model should be used with some caution. Simply switching it on, one will 
not retain full description of single particle observables in minimum bias pp 
collisions, as the excitation gluons will increase multiplicity. Besides normal 
tuning, one can use the parameter <code>FragmentationSystems:mJoin</code> to 
join the excitation gluons together, in order to recover single particle 
observables. 
 
<br/><br/><table><tr><td><strong>Ropewalk:rCutOff </td><td></td><td> <input type="text" name="6" value="6.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>6.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 100.</code>)</td></tr></table>
This parameter gives the maximum cut-off radius, at which strings stops 
interacting. The purpose of the parameter is to decrease computation time by 
not calculating arbitrarily small pushes. In pp collisions at LHC energies, 
no significant variation in the results is observed by increasing this value 
above the default. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:gAmplitude </td><td></td><td> <input type="text" name="7" value="5.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 100.</code>)</td></tr></table>
The amplitude of the shoving force. Note that many traditional 
Min Bias/UE observables, such as multiplicity and <i>p_\perp</i>, as 
well as transverse quantities, are sensitive to this parameter. As 
such, a change of this, in principle warrants a full retuning of the 
MPI framework. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:gExponent </td><td></td><td> <input type="text" name="8" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 100.</code>)</td></tr></table>
This value multiplies the string radius in the shoving 
force, allowing for a variation between string radius in the flavour rope 
treatment and the shoving treatment, if one wishes to run both simultaneously. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:deltay </td><td></td><td> <input type="text" name="9" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.01</code>; <code>maximum = 10.</code>)</td></tr></table>
This value gives the width of the rapidity slices in which the event is 
split before shoving. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:tShove </td><td></td><td> <input type="text" name="10" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 100.</code>)</td></tr></table>
The total shoving time in units of fm/c. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:deltat </td><td></td><td> <input type="text" name="11" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>; <code>minimum = 0.01</code>; <code>maximum = 100.0</code>)</td></tr></table>
The size of the steps taken in time during shoving. Since the whole 
event needs to be retraced after every time step, this should not be 
too small. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:tInit </td><td></td><td> <input type="text" name="12" value="1.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.5</strong></code>; <code>minimum = 0.</code>; <code>maximum = 100.</code>)</td></tr></table>
The strings are allowed to propagate for some time, given in fm/c 
by this parameter, before shoving takes place. This accounts for the 
fact that the strings are created with a vanishing transverse size, 
and only shove each other when their transverse size is large enough 
for interaction. Furthermore, the physical value of this parameter is 
largely connected to the values set for 
<?php $filepath = $_GET["filepath"];
echo "<a href='VetexInformation.php?filepath=".$filepath."' target='page'>";?>Vertex Information</a>. 
   
 
<br/><br/><strong>Ropewalk:shoveGluonLoops</strong>  <input type="radio" name="13" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="13" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow for shoving of strings which form a gluon loop. 
This is mainly a technical setting, and should be kept switched on, 
unless the user has a specific intention of switching it off. 
   
 
<br/><br/><strong>Ropewalk:shoveJunctionStrings</strong>  <input type="radio" name="14" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="14" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow for shoving of strings that includes a junction topology from 
eg. beam remnants. This is mainly a technical setting, and should be 
kept switched on, unless the user has a specific intention of 
switching it off. 
   
 
<br/><br/><strong>Ropewalk:shoveMiniStrings</strong>  <input type="radio" name="15" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="15" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow for shoving of ministrings. This is mainly a technical setting, and 
should be kept switched on, unless the user has a specific intention of 
switching it off. 
   
 
<br/><br/><strong>Ropewalk:limitMom</strong>  <input type="radio" name="16" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="16" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
It is possible to switch off shoving for dipoles with a <i>p_\perp</i> 
above a given value. This is intended as a cut-off to disallow string segments 
moving so fast that they would anyway escape shoving from soft strings to 
have gluonic excitations added to them. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:pTcut </td><td></td><td> <input type="text" name="17" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1000.</code>)</td></tr></table>
The value of <i>p_\perp</i> at which shoving is turned off, if the flag 
<code>Ropewalk:limitMom</code> is on. 
   
 
<a name="section2"></a> 
<h3>Flavour Ropes</h3> 
 
The Flavour Ropes mechanism allows strings situated close in impact parameter 
space to interact coherently, forming a rope, which hadronizes with a larger, 
effective string tension. The model is described in ref. [<a href="Bibliography.php#refBie14" target="page">Bie14</a>], 
building on an older idea by Biro et al. [<a href="Bibliography.php#refBir84" target="page">Bir84</a>]. 
 
<p/> 
In the flavour rope formalism, a rope is described as an SU(3) multiplet, 
characterized uniquely by two quantum numbers <i>p</i> and <i>q</i>. 
The quantum numbers are calculated, following ref. [<a href="Bibliography.php#refBir84" target="page">Bir84</a>], by 
a random walk procedure in colour space, taking <i>m, n</i> steps, 
where <i>m</i> and <i>n</i> signify the number of overlapping 
strings which are respectively parallel or anti-parallel to the 
hadronizing string. 
 
<p/> 
When the rope quantum numbers have been determined, the effective string 
tension is calculated per individual breaking, using a lattice QCD 
determination of the string tension [<a href="Bibliography.php#refBal00" target="page">Bal00</a>]. The effective 
string tension is then used to rescale the hadronization parameters 
described in the section on <?php $filepath = $_GET["filepath"];
echo "<a href='Fragmentation.php?filepath=".$filepath."' target='page'>";?>String 
Fragmentation</a>. One point to note regarding the rescaling is the 
fragmentation parameter <code>StringFlav:probQQtoQ</code>, describing 
baryon relative to meson production. Baryon production is, as suggested 
by eg. the popcorn hadronization model [<a href="Bibliography.php#refEde97" target="page">Ede97</a>],  complicated 
than meson production. The current modelling of this in the flavour ropes 
framework is limited, but intended to be extended in the future. 
 
<br/><br/><table><tr><td><strong>Ropewalk:beta </td><td></td><td> <input type="text" name="18" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.0</code>)</td></tr></table>
In the current implementation of the rope model, the theoretical ignorance 
about baryon production has been parameterized, assuming that the parameter 
<code>StringFlav:probQQtoQ</code> will factorize into two parts, 
one which will scale with effective string tension, one which will not. 
This parameter controls how large a fraction of the parameter will scale 
with string tension. 
   
 
<br/><br/><strong>Ropewalk:alwaysHighest</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Setting this flag on will skip the random walk procedure for flavour ropes, 
and assume that one always ends up in the highest possible SU(3) multiplet. 
This would be adequate for situations where all lower multiplets are assumed 
handled by colour reconnection and junction formation. 
   
 
<br/><br/><strong>Ropewalk:doBuffon</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Setting this flag on, enables a simpler treatment of flavour ropes. This is not 
reliant on vertex information, but string-string overlaps are decided randomly 
&aacute; la Buffon's needle [<a href="Bibliography.php#refBie16c" target="page">Bie16c</a>]: All strings are thrown randomly 
into a circular area in transverse space to estimate overlaps. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:stringProtonRatio </td><td></td><td> <input type="text" name="21" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.0</code>)</td></tr></table>
Only used if <code>Ropewalk:buffonRope</code> is enabled. The ratio of the 
string transverse area to a proton transverse area. Determines the amount of 
overlap in collisions. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:rapiditySpan </td><td></td><td> <input type="text" name="22" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.0</code>)</td></tr></table>
Only used if <code>Ropewalk:buffonRope</code> is enabled. Determines how far in 
rapidity from a string break overlaps are counted. 
   
 
<br/><br/><strong>Ropewalk:setFixedKappa</strong>  <input type="radio" name="23" value="on"><strong>On</strong>
<input type="radio" name="23" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Setting this flag gives the user the possibility to ignore the generator 
space-time information altogether, using only a provided string tension. 
This could be useful for (toy) studies of hadronization in very dense 
environments, such as central heavy ion collisions. 
   
 
<br/><br/><table><tr><td><strong>Ropewalk:presetKappa </td><td></td><td> <input type="text" name="24" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 100.0</code>)</td></tr></table>
The effective string tension is normally calculated dynamically using overlaps 
of strings, based on <?php $filepath = $_GET["filepath"];
echo "<a href='VertexInformation.php?filepath=".$filepath."' target='page'>";?>Parton Vertex</a> 
information. By setting <code>Ropewalk:setFixedKappa</code>, this information 
is ignored, and a preset value provided in the <code>presetKappa</code> 
variable is used. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:kappa </td><td></td><td> <input type="text" name="25" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 10.</code>)</td></tr></table>
A base value of the string tension can be added, and modified along with other 
parameters, to allow for studies of exotic quark production in the Rope model. 
   
 
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

if($_POST["1"] != "off")
{
$data = "Ropewalk:RopeHadronization = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "Ropewalk:doShoving = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "Ropewalk:doFlavour = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0.5")
{
$data = "Ropewalk:r0 = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "0.2")
{
$data = "Ropewalk:m0 = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "6.0")
{
$data = "Ropewalk:rCutOff = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "5.0")
{
$data = "Ropewalk:gAmplitude = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "1.0")
{
$data = "Ropewalk:gExponent = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "0.2")
{
$data = "Ropewalk:deltay = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "1.0")
{
$data = "Ropewalk:tShove = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "0.1")
{
$data = "Ropewalk:deltat = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "1.5")
{
$data = "Ropewalk:tInit = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "on")
{
$data = "Ropewalk:shoveGluonLoops = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "on")
{
$data = "Ropewalk:shoveJunctionStrings = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "on")
{
$data = "Ropewalk:shoveMiniStrings = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "on")
{
$data = "Ropewalk:limitMom = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "2.0")
{
$data = "Ropewalk:pTcut = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "0.2")
{
$data = "Ropewalk:beta = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "Ropewalk:alwaysHighest = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "Ropewalk:doBuffon = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "0.2")
{
$data = "Ropewalk:stringProtonRatio = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "0.5")
{
$data = "Ropewalk:rapiditySpan = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "off")
{
$data = "Ropewalk:setFixedKappa = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "0.")
{
$data = "Ropewalk:presetKappa = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "0.2")
{
$data = "StringFlav:kappa = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
