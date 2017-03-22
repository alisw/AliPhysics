<html>
<head>
<title>Left-Right-Symmetry Processes</title>
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

<form method='post' action='LeftRightSymmetryProcesses.php'>
 
<h2>Left-Right-Symmetry Processes</h2> 
 
At current energies, the world is left-handed, i.e. the Standard Model 
contains an <i>SU(2)_L</i> group. Left-right symmetry at some larger 
scale implies the need for an <i>SU(2)_R</i> group. Thus the particle 
content is expanded by right-handed <i>Z_R^0</i> and <i>W_R^+-</i> 
and right-handed neutrinos. The Higgs fields have to be in a triplet 
representation, leading to doubly-charged Higgs particles, one set for 
each of the two <i>SU(2)</i> groups. Also the number of neutral and 
singly-charged Higgs states is increased relative to the Standard Model, 
but a search for the lowest-lying states of this kind is no different 
from e.g. the freedom already accorded by the MSSM Higgs scenarios. 
 
<p/> 
PYTHIA implements the scenario of [<a href="Bibliography.php" target="page">Hui97</a>]. 
 
<p/> 
The <i>W_R^+-</i> has been implemented as a simple copy of the 
ordinary <i>W^+-</i>, with the exception that it couples to 
right-handed neutrinos instead of the ordinary left-handed ones. 
Thus the standard CKM matrix is used in the quark sector, and the 
same vector and axial coupling strengths, leaving only the mass as 
free parameter. The <i>Z_R^0</i> implementation (without interference 
with the photon or the ordinary <i>Z^0</i>) allows decays both to 
left- and right-handed neutrinos, as well as other fermions, according 
to one specific model ansatz. Obviously both the <i>W_R^+-</i> 
and the <i>Z_R^0</i> descriptions are  likely to be simplifications, 
but provide a starting point. 
 
<p/> 
For the doubly-charged Higgs bosons, the main decay modes implemented are 
<i>H_L^++ &rarr; W_L^+ W_L^+, l_i^+ l_j^+ </i> (<i>i, j</i> generation 
indices) and <i>H_R^++ &rarr; W_R^+ W_R^+, l_i^+ l_j^+</i>. 
 
<p/> 
The right-handed neutrinos can be allowed to decay further. Assuming them 
to have a mass below that of <i>W_R^+-</i>, they decay to three-body 
states via a virtual <i>W_R^+-</i>, <i>nu_Rl &rarr; l+- f fbar'</i>, 
where both lepton charges are allowed owing to the Majorana character 
of the neutrinos. If there is a significant mass splitting, also 
sequential decays <i>nu_Rl &rarr; l+- l'-+  nu'_Rl</i> are allowed. 
Currently the decays are isotropic in phase space. If the neutrino 
masses are close to or above the <i>W_R^</i> ones, this description 
has to be substituted by a sequential decay via a real <i>W_R^</i> 
(not implemented, but actually simpler to do than the one here). 
 
 
<h3>Production processes</h3> 
 
A few different production processes have been implemented, which normally 
would not overlap and therefore could be run together. 
 
<br/><br/><strong>LeftRightSymmmetry:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of implemented processes within a 
left-right-symmetric scenario. 
   
 
<br/><br/><strong>LeftRightSymmmetry:ffbar2ZR</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar &rarr; Z_R^0</i>. 
Code 3101. 
   
 
<br/><br/><strong>LeftRightSymmmetry:ffbar2WR</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i><f fbar' &rarr; W_R^+</i>. 
Code 3102. 
   
 
<br/><br/><strong>LeftRightSymmmetry:ll2HL</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>l_i l_j &rarr; H_L^--</i>. 
Code 3121. 
   
 
<br/><br/><strong>LeftRightSymmmetry:lgm2HLe</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>l_i gamma &rarr; H_L^-- e^+</i>. 
Code 3122. 
   
 
<br/><br/><strong>LeftRightSymmmetry:lgm2HLmu</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>l_i gamma &rarr; H_L^-- mu^+</i>. 
Code 3123. 
   
 
<br/><br/><strong>LeftRightSymmmetry:lgm2HLtau</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>l_i gamma &rarr; H_L^-- tau^+</i>. 
Code 3124. 
   
 
<br/><br/><strong>LeftRightSymmmetry:ff2HLff</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f_1 f_2 &rarr; H_L^-- f_3 f_4</i> via <i>WW</i> fusion. 
Code 3125. 
   
 
<br/><br/><strong>LeftRightSymmmetry:ffbar2HLHL</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar &rarr;  H_L^++ H_L^--</i>. 
Code 3126. 
   
 
<br/><br/><strong>LeftRightSymmmetry:ll2HR</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>l_i l_j &rarr; H_R^--</i>. 
Code 3141. 
   
 
<br/><br/><strong>LeftRightSymmmetry:lgm2HRe</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>l_i gamma &rarr; H_R^-- e^+</i>. 
Code 3142. 
   
 
<br/><br/><strong>LeftRightSymmmetry:lgm2HRmu</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>l_i gamma &rarr; H_R^-- mu^+</i>. 
Code 3143. 
   
 
<br/><br/><strong>LeftRightSymmmetry:lgm2HRtau</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>l_i gamma &rarr; H_R^-- tau^+</i>. 
Code 3144. 
   
 
<br/><br/><strong>LeftRightSymmmetry:ff2HRff</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f_1 f_2 &rarr; H_R^-- f_3 f_4</i> via <i>WW</i> fusion. 
Code 3145. 
   
 
<br/><br/><strong>LeftRightSymmmetry:ffbar2HRHR</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar &rarr;  H_R^++ H_R^--</i>. 
Code 3146. 
   
 
<h3>Parameters</h3> 
 
The basic couplings of the model are 
 
<br/><br/><table><tr><td><strong>LeftRightSymmmetry:gL </td><td></td><td> <input type="text" name="16" value="0.64" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.64</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
lefthanded coupling <i>g_L = e / sin(theta)</i>. 
   
 
<br/><br/><table><tr><td><strong>LeftRightSymmmetry:gR </td><td></td><td> <input type="text" name="17" value="0.64" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.64</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
righthanded coupling <i>g_R</i>, assumed the same as <i>g_L</i>. 
   
 
<br/><br/><table><tr><td><strong>LeftRightSymmmetry:vL </td><td></td><td> <input type="text" name="18" value="5." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5.</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
vacuum expectation value <i>v_L</i> (in GeV) for the left-triplet. 
   
 
<p/> 
The corresponding vacuum expectation value <i>v_R</i> is assumed 
given by <i>v_R = sqrt(2) M_WR / g_R</i> and is not stored explicitly. 
 
<p/> 
The Yukawa couplings of a lepton pair to a <i>H^--</i>, assumed the 
same for <i>H_L^--</i> and <i>H_R^--</i>, is described by a symmetric 
3-by-3 matrix. The default matrix is dominated by the diagonal elements 
and especially by the <i>tau tau</i> one. 
 
<br/><br/><table><tr><td><strong>LeftRightSymmmetry:coupHee </td><td></td><td> <input type="text" name="19" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Yukawa coupling for <i>H^-- &rarr; e- e-</i>. 
   
 
<br/><br/><table><tr><td><strong>LeftRightSymmmetry:coupHmue </td><td></td><td> <input type="text" name="20" value="0.01" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Yukawa coupling for <i>H^-- &rarr; mu- e-</i>. 
   
 
<br/><br/><table><tr><td><strong>LeftRightSymmmetry:coupHmumu </td><td></td><td> <input type="text" name="21" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Yukawa coupling for <i>H^-- &rarr; mu- mu-</i>. 
   
 
<br/><br/><table><tr><td><strong>LeftRightSymmmetry:coupHtaue </td><td></td><td> <input type="text" name="22" value="0.01" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Yukawa coupling for <i>H^-- &rarr; tau- e-</i>. 
   
 
<br/><br/><table><tr><td><strong>LeftRightSymmmetry:coupHtaumu </td><td></td><td> <input type="text" name="23" value="0.01" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Yukawa coupling for <i>H^-- &rarr; tau- mu-</i>. 
   
 
<br/><br/><table><tr><td><strong>LeftRightSymmmetry:coupHtautau </td><td></td><td> <input type="text" name="24" value="0.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.3</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Yukawa coupling for <i>H^-- &rarr; tau- tau-</i>. 
   
 
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
$data = "LeftRightSymmmetry:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "LeftRightSymmmetry:ffbar2ZR = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "LeftRightSymmmetry:ffbar2WR = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "LeftRightSymmmetry:ll2HL = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "LeftRightSymmmetry:lgm2HLe = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "LeftRightSymmmetry:lgm2HLmu = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "LeftRightSymmmetry:lgm2HLtau = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "LeftRightSymmmetry:ff2HLff = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "LeftRightSymmmetry:ffbar2HLHL = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "LeftRightSymmmetry:ll2HR = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "LeftRightSymmmetry:lgm2HRe = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "LeftRightSymmmetry:lgm2HRmu = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "LeftRightSymmmetry:lgm2HRtau = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "LeftRightSymmmetry:ff2HRff = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "LeftRightSymmmetry:ffbar2HRHR = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "0.64")
{
$data = "LeftRightSymmmetry:gL = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "0.64")
{
$data = "LeftRightSymmmetry:gR = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "5.")
{
$data = "LeftRightSymmmetry:vL = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "0.1")
{
$data = "LeftRightSymmmetry:coupHee = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "0.01")
{
$data = "LeftRightSymmmetry:coupHmue = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "0.1")
{
$data = "LeftRightSymmmetry:coupHmumu = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "0.01")
{
$data = "LeftRightSymmmetry:coupHtaue = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "0.01")
{
$data = "LeftRightSymmmetry:coupHtaumu = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "0.3")
{
$data = "LeftRightSymmmetry:coupHtautau = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
 
