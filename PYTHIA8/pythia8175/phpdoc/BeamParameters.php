<html>
<head>
<title>Beam Parameters</title>
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

<form method='post' action='BeamParameters.php'>

<h2>Beam Parameters</h2>

The settings on this page relate to the beam identities and energies, 
to a beam momentum spread and to a beam interaction spot. 
As always, momenta and energies are to be given in units of GeV,
and of space and time in mm. 

<h3>Incoming beams</h3>

There are two ways to set the identities and energies of the two incoming 
beam particles. One is to use the <code>init()</code> method with no 
arguments. Then the settings variables below will be read and used. The 
alternative is to call <code><?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?>init(...)</a></code> 
with arguments that provide this information. Then you need not use the 
variables below (although it would still be possible). Note that, if nothing 
is done, you will default to LHC at 14 TeV.

<p/>
Currently the beam particles must be either a hadron pair or a lepton
pair. In the former category <i>p p</i> and <i>pbar p</i> 
combinations dominate, but it is also possible to combine with 
<i>pi^+</i>, <i>pi^-</i> and <i>pi^0</i>. In the latter 
<i>e^+ e^-</i> and <i>mu^+ mu^-</i> would be the most useful 
combinations, but also others should work if combined with an 
appropriate hard process.

<br/><br/><table><tr><td><strong>Beams:idA  </td><td></td><td> <input type="text" name="1" value="2212" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2212</strong></code>)</td></tr></table>
The PDG <code>id</code> code for the first incoming particle.
Allowed codes include 
<br/><i>2212 = p</i>, <i>-2212 = pbar</i>, 
<br/><i>2112 = n</i>, <i>-2112 = nbar</i>, 
<br/><i>211 = pi^+</i>, <i>-211 = pi^-</i>, <i>111 = pi^0</i>,  
<br/><i>990 = Pomeron</i> (used in diffractive machinery;
here mainly for debug purposes),
<br/><i>11 = e^-</i>, <i>-11 = e^+</i>, 
<br/><i>13 = mu^-</i>, <i>-13 = mu^+</i>, 
<br/>and a few more leptons/neutrinos in a few combinations.
  

<br/><br/><table><tr><td><strong>Beams:idB  </td><td></td><td> <input type="text" name="2" value="2212" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2212</strong></code>)</td></tr></table>
The PDG <code>id</code> code for the second incoming particle.
  
 
<br/><br/><table><tr><td><strong>Beams:frameType  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 5</code>)</td></tr></table>
Choice of frame for the two colliding particles. For options
1 - 3 the beam identities are specified above, while they are
obtained by the Les Houches information for options 4 and 5.
<br/>
<input type="radio" name="3" value="1" checked="checked"><strong>1 </strong>: the beams are colliding in their CM frame,  and therefore only the CM energy needs to be provided, see  <code>Beams:eCM</code> below. <br/>
<input type="radio" name="3" value="2"><strong>2 </strong>: the beams are back-to-back, but with different energies, see <code>Beams:eA</code> and <code>Beams:eB</code> below. This option could also be used for fixed-target configurations. <br/>
<input type="radio" name="3" value="3"><strong>3 </strong>: the beams are not back-to-back, and therefore the three-momentum of each incoming particle needs to be specified, see <code>Beams:pxA</code> through <code>Beams:pzB</code> below. <br/>
<input type="radio" name="3" value="4"><strong>4 </strong>: the beam and event information is stored in a  <aloc href="LesHouchesAccord">Les Houches Event File</aloc>,  see <code>Beams:LHEF</code> below. <br/>
<input type="radio" name="3" value="5"><strong>5 </strong>: the beam and event information is obtained by a pointer to an <code><aloc href="LesHouchesAccord">LHAup</aloc></code>  class instance.   <br/>

<br/><br/><table><tr><td><strong>Beams:eCM </td><td></td><td> <input type="text" name="4" value="14000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>14000.</strong></code>; <code>minimum = 10.</code>)</td></tr></table>
Collision CM energy, to be set if <code>Beams:frameType</code> = 1. 
  

<br/><br/><table><tr><td><strong>Beams:eA </td><td></td><td> <input type="text" name="5" value="7000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>7000.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The energy of the first incoming particle, moving in the 
<i>+z </i>direction, to be set if <code>Beams:frameType</code> = 2. 
If the particle energy is smaller than its mass
it is assumed to be at rest. 
  

<br/><br/><table><tr><td><strong>Beams:eB </td><td></td><td> <input type="text" name="6" value="7000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>7000.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The energy of the second incoming particle, moving in the 
<i>-z</i> direction, to be set if <code>Beams:frameType</code> = 2. 
If the particle energy is smaller than its mass
it is assumed to be at rest.
  

<br/><br/><table><tr><td><strong>Beams:pxA </td><td></td><td> <input type="text" name="7" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>p_x</i> component of the first incoming particle, 
to be set if <code>Beams:frameType</code> = 3. 
  

<br/><br/><table><tr><td><strong>Beams:pyA </td><td></td><td> <input type="text" name="8" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>p_y</i> component of the first incoming particle, 
to be set if <code>Beams:frameType</code> = 3. 
  

<br/><br/><table><tr><td><strong>Beams:pzA </td><td></td><td> <input type="text" name="9" value="7000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>7000.</strong></code>)</td></tr></table>
The <i>p_z</i> component of the first incoming particle, 
to be set if <code>Beams:frameType</code> = 3. 
  

<br/><br/><table><tr><td><strong>Beams:pxB </td><td></td><td> <input type="text" name="10" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>p_x</i> component of the second incoming particle, 
to be set if <code>Beams:frameType</code> = 3. 
  

<br/><br/><table><tr><td><strong>Beams:pyB </td><td></td><td> <input type="text" name="11" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>p_y</i> component of the second incoming particle, 
to be set if <code>Beams:frameType</code> = 3. 
  

<br/><br/><table><tr><td><strong>Beams:pzB </td><td></td><td> <input type="text" name="12" value="-7000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-7000.</strong></code>)</td></tr></table>
The <i>p_z</i> component of the second incoming particle, 
to be set if <code>Beams:frameType</code> = 3. 
  

<br/><br/><table><tr><td><strong>Beams:LHEF  </td><td></td><td> <input type="text" name="13" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
The name of a Les Houches Event File,
to be set if <code>Beams:frameType</code> = 4. 
  

<br/><br/><table><tr><td><strong>Beams:LHEFheader  </td><td></td><td> <input type="text" name="14" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
As some information in a Les Houches Event File init block is only known
at the end of generation, some programs choose to output this as a
separate file. If <code>Beams:LHEFheader</code> is given, information up
till the end of the init block is read from this file, with
the events themselves read as usual from the file given in
<code>Beams:LHEF</code>.
  

<br/><br/><strong>Beams:newLHEFsameInit</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow to begin reading events from a new LHEF or or a new 
<code>LHAup</code> instance without a completely new initialization. 
Only useful when <code>Beams:frameType</code> = 4 or 5.
  

<br/><br/><strong>Beams:readLHEFheaders</strong>  <input type="radio" name="16" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="16" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Read in LHEF header blocks and store them in the
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>Info</a> class. See also
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>LHAupLHEF</a> for more information.
  

<p/><code>mode&nbsp; </code><strong> Beams:nSkipLHEFatInit &nbsp;</strong> 
 (<code>default = <strong>0</strong></code>)<br/>
Skip the first <i>nSkip</i> events of the input stream 
(cf. the <code>LHAup::skipEvent(nSkip)</code> method).
Only used when <code>Beams:frameType</code> = 4 or 5.
  

<h3>Beam momentum spread</h3>

This framework currently is intended for a modest beam spread, such as
experienced at hadron colliders. Thus it can be safely assumed that the 
physics does not change over the CM energy range probed, so that the 
parameters of the physics initialization at the nominal energy can be
used as is. Currently it can <b>not</b> be used for the more extensive
energy spread expected at linear <i>e^+ e^-</i> colliders. Also,
any attempt to combine it with external Les Houches input of 
parton-level events is at own risk. 
     
<p/>
On this page you can set the momentum spread according to a simple 
Gaussian distribution. If you instead want a more sophisticated 
parametrization, you can write and link your own 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='BeamShape.php?filepath=".$filepath."' target='page'>";?>BeamShape</a></code> class.

<br/><br/><strong>Beams:allowMomentumSpread</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow the beam momenta to be smeared around their initialization
nominal values. 
  

<br/><br/><table><tr><td><strong>Beams:sigmaPxA </td><td></td><td> <input type="text" name="18" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width of a Gaussian distribution of the <i>p_x</i> spread of the
first incoming particle.
  

<br/><br/><table><tr><td><strong>Beams:sigmaPyA </td><td></td><td> <input type="text" name="19" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width of a Gaussian distribution of the <i>p_y</i> spread of the
first incoming particle.
  

<br/><br/><table><tr><td><strong>Beams:sigmaPzA </td><td></td><td> <input type="text" name="20" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width of a Gaussian distribution of the <i>p_z</i> spread of the
first incoming particle.
  

<br/><br/><table><tr><td><strong>Beams:maxDevA </td><td></td><td> <input type="text" name="21" value="5." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The triply Gaussian distribution <i>(p_x, p_y, p_z)</i> is restricted to 
a maximal total deviation from the nominal values <i>(p_x0, p_y0, p_z0)</i>
for the first incoming particle, like
<br/><i>
(p_x - p_x0)^2/sigma_px^2 + (p_y - p_y0)^2/sigma_py^2 +
(p_z - p_z0)^2/sigma_pz^2 < maxDev^2
</i><br/>
(Note the absence of a factor 2 in the denominator, unlike the Gaussians 
used to pick <i>(p_x, p_y, p_z)</i>.) 
  

<br/><br/><table><tr><td><strong>Beams:sigmaPxB </td><td></td><td> <input type="text" name="22" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width of a Gaussian distribution of the <i>p_x</i> spread of the
second incoming particle.
  

<br/><br/><table><tr><td><strong>Beams:sigmaPyB </td><td></td><td> <input type="text" name="23" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width of a Gaussian distribution of the <i>p_y</i> spread of the
second incoming particle.
  

<br/><br/><table><tr><td><strong>Beams:sigmaPzB </td><td></td><td> <input type="text" name="24" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width of a Gaussian distribution of the <i>p_z</i> spread of the
second incoming particle.
  

<br/><br/><table><tr><td><strong>Beams:maxDevB </td><td></td><td> <input type="text" name="25" value="5." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The triply Gaussian distribution <i>(p_x, p_y, p_z)</i> is restricted to 
a maximal total deviation from the nominal values <i>(p_x0, p_y0, p_z0)</i>, 
for the second incoming particle, like
<br/><i>
(p_x - p_x0)^2/sigma_px^2 + (p_y - p_y0)^2/sigma_py^2 +
(p_z - p_z0)^2/sigma_pz^2 < maxDev^2
</i><br/>
(Note the absence of a factor 2 in the denominator, unlike the Gaussians 
used to pick <i>(p_x, p_y, p_z)</i>.) 
  

<h3>Beam interaction vertex</h3>

On this page you can set the spread of the interaction vertex according to 
a simple Gaussian distribution. If you instead want a more sophisticated 
parametrization, you can write and link your own 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='BeamShape.php?filepath=".$filepath."' target='page'>";?>BeamShape</a></code> class.

<br/><br/><strong>Beams:allowVertexSpread</strong>  <input type="radio" name="26" value="on"><strong>On</strong>
<input type="radio" name="26" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow the interaction vertex of the two colliding beams to be smeared.
If off, then the vertex is set to be the origin.
  

<br/><br/><table><tr><td><strong>Beams:sigmaVertexX </td><td></td><td> <input type="text" name="27" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width of a Gaussian distribution of the <i>x</i> location of the
interaction vertex. 
  

<br/><br/><table><tr><td><strong>Beams:sigmaVertexY </td><td></td><td> <input type="text" name="28" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width of a Gaussian distribution of the <i>y</i> location of the
interaction vertex. 
  

<br/><br/><table><tr><td><strong>Beams:sigmaVertexZ </td><td></td><td> <input type="text" name="29" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width of a Gaussian distribution of the <i>z</i> location of the
interaction vertex. 
  

<br/><br/><table><tr><td><strong>Beams:maxDevVertex </td><td></td><td> <input type="text" name="30" value="5." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The triply Gaussian distribution of interaction vertex position 
<i>(x, y, z)</i> is restricted to a maximal total deviation from the 
origin, like
<br/><i>
x^2/sigma_x^2 + y^2/sigma_y^2 + z^2/sigma_z^2 < maxDevVertex^2
</i><br/>
(Note the absence of a factor 2 in the denominator, unlike the Gaussians 
used to pick <i>(x, y, z)</i>.) 
  

<br/><br/><table><tr><td><strong>Beams:sigmaTime </td><td></td><td> <input type="text" name="31" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width of a Gaussian distribution of the collision time (in units of
mm/c). Note that, if the above space parametrization is viewed as the
effect of two incoming beams along the <i>+-z</i> axis, with each beam
having a Gaussian spread, then the spread of the time would also become 
a Gaussian with the same width as the <i>z</i> one (times the 
velocity of the beams, which we expect is close to unity). For flexibility
we have not enforced any such relation, however. 
  

<br/><br/><table><tr><td><strong>Beams:maxDevTime </td><td></td><td> <input type="text" name="32" value="5." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The collision time is restricted to be in the range 
<i>|t| &lt; sigma_t * maxDevTime</i>. 
  

<p/>
The distributions above are all centered at the origin. It is also 
possible to shift the above distributions to be centered around another
nominal position. You must have <code>Beams:allowVertexSpread = on</code>
to use this possibility.

<br/><br/><table><tr><td><strong>Beams:offsetVertexX </td><td></td><td> <input type="text" name="33" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>x</i> location of the interaction vertex is centered at this value. 
  

<br/><br/><table><tr><td><strong>Beams:offsetVertexY </td><td></td><td> <input type="text" name="34" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>y</i> location of the interaction vertex is centered at this value. 
  

<br/><br/><table><tr><td><strong>Beams:offsetVertexZ </td><td></td><td> <input type="text" name="35" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>z</i> location of the interaction vertex is centered at this value. 
  

<br/><br/><table><tr><td><strong>Beams:offsetTime </td><td></td><td> <input type="text" name="36" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The time <i>t</i> of the interaction vertex is centered at this value. 
  

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

if($_POST["1"] != "2212")
{
$data = "Beams:idA = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "2212")
{
$data = "Beams:idB = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "1")
{
$data = "Beams:frameType = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "14000.")
{
$data = "Beams:eCM = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "7000.")
{
$data = "Beams:eA = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "7000.")
{
$data = "Beams:eB = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "0.")
{
$data = "Beams:pxA = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "0.")
{
$data = "Beams:pyA = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "7000.")
{
$data = "Beams:pzA = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "0.")
{
$data = "Beams:pxB = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "0.")
{
$data = "Beams:pyB = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "-7000.")
{
$data = "Beams:pzB = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "void")
{
$data = "Beams:LHEF = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "void")
{
$data = "Beams:LHEFheader = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "Beams:newLHEFsameInit = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "on")
{
$data = "Beams:readLHEFheaders = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "Beams:allowMomentumSpread = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "0.")
{
$data = "Beams:sigmaPxA = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "0.")
{
$data = "Beams:sigmaPyA = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "0.")
{
$data = "Beams:sigmaPzA = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "5.")
{
$data = "Beams:maxDevA = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "0.")
{
$data = "Beams:sigmaPxB = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "0.")
{
$data = "Beams:sigmaPyB = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "0.")
{
$data = "Beams:sigmaPzB = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "5.")
{
$data = "Beams:maxDevB = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "off")
{
$data = "Beams:allowVertexSpread = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "0.")
{
$data = "Beams:sigmaVertexX = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "0.")
{
$data = "Beams:sigmaVertexY = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "0.")
{
$data = "Beams:sigmaVertexZ = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "5.")
{
$data = "Beams:maxDevVertex = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "0.")
{
$data = "Beams:sigmaTime = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "5.")
{
$data = "Beams:maxDevTime = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "0.")
{
$data = "Beams:offsetVertexX = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "0.")
{
$data = "Beams:offsetVertexY = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "0.")
{
$data = "Beams:offsetVertexZ = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "0.")
{
$data = "Beams:offsetTime = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
