<html>
<head>
<title>Dark Matter Processes</title>
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

<form method='post' action='DarkMatterProcesses.php'>
 
<h2>Dark Matter Processes</h2> 
<ol id="toc">
  <li><a href="#section0">Scalar Mediator</a></li>
  <li><a href="#section1">Vector Mediator <ei>Z'</ei></a></li>
  <li><a href="#section2">Drell-Yan production of charged co-annihilation partners </a></li>
</ol>

 
This page contains the production of Dirac fermion Dark Matter via new 
<i>s</i>-channel mediators. A summary of the physics scenarios and 
parameters involved can be found in [<a href="Bibliography.php#refDes18" target="page">Des18</a>]. Examples how these 
processes can be run are found in <code>main75.cc</code> and 
<code>main76.cc</code>. 
 
<p/> 
The particles in the scenarios considered here are a mix of established 
PDG ones (51 - 55) and extensions thereof ( 56 - 59), as follows: 
<table> 
<tr><td>51</td><td>Scalar DM (currently unused);</td></tr> 
<tr><td>52</td><td>Fermionic DM (<i>chi_1</i>);</td></tr> 
<tr><td>53</td><td>Vector DM (currently unused);</td></tr> 
<tr><td>54</td><td>Scalar (or pseudoscalar) mediator (<i>S</i>);</td></tr> 
<tr><td>55</td><td>Vector (or axial vector) mediator (<i>Z'</i>); </td></tr> 
<tr><td>56</td><td>Charged scalar partner (<i>l^~</i>);</td></tr> 
<tr><td>57</td><td>Singly charged partner (<i>chi^+</i>);</td></tr> 
<tr><td>58</td><td>Neutral partner (<i>chi_2</i>);</td></tr> 
<tr><td>59</td><td>Doubly charged partner (<i>chi^++</i>).</td></tr> 
</table> 
 
<a name="section0"></a> 
<h3>Scalar Mediator</h3> 
 
<br/><br/><strong>DM:gg2S2XX</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g &rarr; S &rarr; X Xbar</i>. 
Code 6011. S is assumed to be on-shell. 
   
 
<br/><br/><strong>DM:gg2S2XXj</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g &rarr;S &rarr; X Xbar j</i>. 
Code 6012. S is assumed to be on-shell. (Not yet available.) 
   
 
<p/> 
Fermion couplings to scalar S are assumed to be proportional to mass 
of the fermion and couplings are the factor multiplying SM Higgs 
coupling (i.e. sin(mixing) in case of portal models). 
 
<br/><br/><table><tr><td><strong>Sdm:vf </td><td></td><td> <input type="text" name="3" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>)</td></tr></table>
Scalar coupling of SM fermions. 
   
<br/><br/><table><tr><td><strong>Sdm:af </td><td></td><td> <input type="text" name="4" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
Pseudo-scalar coupling of SM fermions. 
   
<br/><br/><table><tr><td><strong>Sdm:vX </td><td></td><td> <input type="text" name="5" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>)</td></tr></table>
Scalar coupling of DM fermion. 
   
<br/><br/><table><tr><td><strong>Sdm:aX </td><td></td><td> <input type="text" name="6" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
Pseudo-scalar coupling of DM fermion. 
   
 
<a name="section1"></a> 
<h3>Vector Mediator <i>Z'</i></h3> 
 
The Vector mediator model assumes a simplified U(1) model with 
couplings to fermionic Dark Matter.  Both vector and axial couplings 
are possible.  Interference with <i>gamma/Z</i> is currently not 
implemented, therefore this should mainly be used when <i>Z' &rarr; X 
Xbar</i>.  However, a quick check of dijet or dilepton cross sections 
can be made by setting the mode <code>Zp:decayMode</code>. 
 
<br/><br/><table><tr><td><strong>Zp:decayMode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
Switch to change decay of the <ei>Z'</ei>. 
<br/>
<input type="radio" name="7" value="0" checked="checked"><strong>0 </strong>: <ei>X Xbar</ei>  <br/>
<input type="radio" name="7" value="1"><strong>1 </strong>: <ei>q qbar</ei> (dijets)  <br/>
<input type="radio" name="7" value="2"><strong>2 </strong>: <ei>l lbar</ei> (charged dileptons)  <br/>
<input type="radio" name="7" value="3"><strong>3 </strong>: <ei>nu nubar + X Xbar</ei> (invisible)  <br/>
 
<br/><br/><strong>DM:ffbar2Zp2XX</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar &rarr;Z'^0 &rarr; X Xbar</i>. 
Code 6001. <i>Z'</i> is assumed to be on-shell. 
   
 
<br/><br/><strong>DM:ffbar2ZpH</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar &rarr;Z'^0 H </i>. 
Code 6004. <i>Z'</i> is assumed to be on-shell. The coupling of the 
<i>Z'</i> to the SM Higgs is given by the parameter <code>Zp:coupH</code>. 
Interference with gamma/Z currently not implemented therefore this is 
only suitable when <i>Z' &rarr; X Xbar</i>.  This can be ensured 
using 
<br/><code>55:onMode = off</code> 
<br/><code>55:onIfAny = 52</code> 
   
 
<br/><br/><strong>DM:ffbar2Zp2XXj</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar &rarr;Z'^0 j &rarr; X Xbar j</i>. 
Code 6002. <i>Z'</i> is assumed to be on-shell. (Not yet available.) 
   
 
<br/><br/><strong>DM:qg2Zp2XXj</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g &rarr;Z'^0 j &rarr; X Xbar j</i>. 
Code 6003. <i>Z'</i> is assumed to be on-shell. (Not yet available.) 
   
 
<p/> 
The vector and axial couplings of fermionic DM to the <i>Z'^0</i> 
can be set freely. The couplings of quarks and leptons can either 
be chosen freely for a new <i>U(1)</i> or be given by kinetic 
mixing with the SM <i>U(1)_Y</i>. The SM fermion couplings are 
assumed to be universal, i.e. generation-independent. The choice of 
fixed axial and vector couplings implies a resonance width that 
increases linearly with the <i>Z'</i> mass. Also some overall 
coupling strengths can be chosen freely. 
 
<br/><br/><table><tr><td><strong>Zp:vX </td><td></td><td> <input type="text" name="12" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
Vector coupling of DM fermion. 
   
<br/><br/><table><tr><td><strong>Zp:aX </td><td></td><td> <input type="text" name="13" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
Axial coupling of DM fermion. 
   
 
<br/><br/><strong>Zp:kineticMixing</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Flag for the two main options to set the couplings of the <i>Z'</i> 
to SM quarks and leptons. In the default off option the overall 
coupling strength of the new <i>U(1)</i> gauge group is given by 
<code>Zp:gZp</code> and the separate fermion couplings by 
<code>Zp:vu</code> through <code>Zp:av</code>. In the alternative, 
with kinetic mixing on, the coupling to the DM is still given by 
<code>Zp:gZp</code>, but the mixing parameter <code>Zp:epsilon</code> 
now specifies how the separate fermion couplings are related to their 
<i>U(1)_Y</i> values. 
   
 
<br/><br/><table><tr><td><strong>Zp:gZp </td><td></td><td> <input type="text" name="15" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>)</td></tr></table>
Gauge coupling of a new <i>U(1)</i>. This parameter also sets the 
coupling of the DM to the <i>Z'</i>, whether kinetic mixing is on 
or not. 
   
 
<br/><br/><table><tr><td><strong>Zp:vu </td><td></td><td> <input type="text" name="16" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
Vector coupling of up-type quarks. 
   
<br/><br/><table><tr><td><strong>Zp:au </td><td></td><td> <input type="text" name="17" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
Axial coupling of up-type quarks. 
   
 
<br/><br/><table><tr><td><strong>Zp:vd </td><td></td><td> <input type="text" name="18" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
Vector coupling of down-type quarks. 
   
<br/><br/><table><tr><td><strong>Zp:ad </td><td></td><td> <input type="text" name="19" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
Axial coupling of down-type quarks. 
   
 
<br/><br/><table><tr><td><strong>Zp:vl </td><td></td><td> <input type="text" name="20" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
Vector coupling of charged leptons. 
   
<br/><br/><table><tr><td><strong>Zp:al </td><td></td><td> <input type="text" name="21" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
Axial coupling of charged leptons. 
   
 
<br/><br/><table><tr><td><strong>Zp:vv </td><td></td><td> <input type="text" name="22" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
Vector coupling of neutrinos. 
   
<br/><br/><table><tr><td><strong>Zp:av </td><td></td><td> <input type="text" name="23" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
Axial coupling of neutrinos. 
   
 
<br/><br/><table><tr><td><strong>Zp:epsilon </td><td></td><td> <input type="text" name="24" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>)</td></tr></table>
Kinetic mixing parameter between the dark <i>U(1)</i> and the SM 
hypercharge <i>U(1)</i>. In the implemented case the <i>Z'^0</i> 
mass is larger than the SM <i>Z^0</i> one. The fermionic current 
for <i>Z'^0</i> is described in  [<a href="Bibliography.php#refCli184" target="page">Cli184</a>]. 
   
 
<br/><br/><table><tr><td><strong>Zp:coupH </td><td></td><td> <input type="text" name="25" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>)</td></tr></table>
Coupling to the SM higgs when kinetic mixing is off. When on the 
coupling is instead set by <code>Zp:epsilon</code>. 
   
 
<a name="section2"></a> 
<h3>Drell-Yan production of charged co-annihilation partners </h3> 
 
We implement two models of co-annihilating Dark Matter where the 
co-annihilation partner carries EW charge and therefore can be 
produced via Drell-Yan production.  The underlying model and 
production process can be selected by choosing the parameter 
<code>DM:DYtype</code>. 
 
<p/> The first model consists of co-annihilation with a scalar with 
leptonic quantum numbers and which couples to a right-handed SM lepton 
and a Dirac fermion Dark Matter via Yukawa couplings.  It is possible 
to choose lepton flavour violating couplings. 
 
<p/> The next model is a generalisation of the mixed gaugino sector 
of Supersymmetry parametrised by one SM singlet and one SU(2) N-plet 
which mix to form Dark matter.  N = 2, 3 and 5 are supported by the 
code.  The main motivation for this choice is to provide a fully 
flexible implementation to calculate production of long-lived 
particles at the LHC.  The resultant spectrum consists of one neutral 
partner, one singly charged partner, and one doubly charged partner 
(in the case of the 5-plet).  The only free parameters are masses of 
the singlet and N-plet and the mixing suppression scale. This 
determines both production and decay of the particles and can cover a 
range of signatures including displaced leptons and vertices, 
long-lived, kinked or disappearing tracks. 
 
<br/><br/><strong>DM:qqbar2DY</strong>  <input type="radio" name="26" value="on"><strong>On</strong>
<input type="radio" name="26" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar &rarr; X_i Xbar_i </i>. 
Code 6020. 
   
 
<br/><br/><table><tr><td><strong>DM:DYtype  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 4</code>)</td></tr></table>
Select the co-annihilation partner to produce via Drell-Yan production. 
<br/>
<input type="radio" name="27" value="1" checked="checked"><strong>1 </strong>:  Scalar lepton production.  Coupling to electrons,  muons and taus can be fixed by <code>DM:yuk1</code>, <code>DM:yuk2</code>,  and <code>DM:yuk3</code> respectively. <br/>
<input type="radio" name="27" value="2"><strong>2 </strong>:  Production of charged partners (i.e. "charginos").<br/>
<input type="radio" name="27" value="3"><strong>3 </strong>:  Production of doubly charged partners.<br/>
<input type="radio" name="27" value="4"><strong>4 </strong>:  Production of neutral and singly-charged partners.<br/>
<!-- Further options to be added 
<input type="radio" name="27" value="5"><strong>5 </strong>:  Production of singly and doubly-charged partners <br/>
--> 
 
<p/> 
The Yukawa couplings can be set using 
<br/><br/><table><tr><td><strong>DM:yuk1 </td><td></td><td> <input type="text" name="28" value="" size="20"/> </td></tr></table>
Electron-DM Yukawa. 
   
<br/><br/><table><tr><td><strong>DM:yuk2 </td><td></td><td> <input type="text" name="29" value="" size="20"/> </td></tr></table>
Muon-DM Yukawa. 
   
<br/><br/><table><tr><td><strong>DM:yuk3 </td><td></td><td> <input type="text" name="30" value="" size="20"/> </td></tr></table>
Tau-DM Yukawa. 
   
 
<p/> 
The parameters for the singlet-N-plet model can be set via the following: 
 
<br/><br/><table><tr><td><strong>DM:Nplet  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>; <code>maximum = 3</code>)</td></tr></table>
<br/>
<option value = "1"> Doublet </option> 
<option value = "2"> Triplet </option> 
<option value = "3"> Quintuplet </option> 
 
<br/><br/><table><tr><td><strong>DM:M1 </td><td></td><td> <input type="text" name="32" value="" size="20"/> </td></tr></table>
Mass of the DM singlet state. 
   
 
<br/><br/><table><tr><td><strong>DM:M2 </td><td></td><td> <input type="text" name="33" value="" size="20"/> </td></tr></table>
Mass of the DM N-plet state. 
   
 
<br/><br/><table><tr><td><strong>DM:Lambda </td><td></td><td> <input type="text" name="34" value="" size="20"/> </td></tr></table>
The suppression scale of the mixing. The Wilson co-efficient is absorbed 
into the suppression scale as there is no independent measurement to 
disentangle it from the scale. 
   
 
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
$data = "DM:gg2S2XX = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "DM:gg2S2XXj = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "0.1")
{
$data = "Sdm:vf = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0")
{
$data = "Sdm:af = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1.0")
{
$data = "Sdm:vX = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0")
{
$data = "Sdm:aX = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "0")
{
$data = "Zp:decayMode = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "DM:ffbar2Zp2XX = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "DM:ffbar2ZpH = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "DM:ffbar2Zp2XXj = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "DM:qg2Zp2XXj = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "1.")
{
$data = "Zp:vX = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0.")
{
$data = "Zp:aX = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "Zp:kineticMixing = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "0.1")
{
$data = "Zp:gZp = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "1.")
{
$data = "Zp:vu = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "0.")
{
$data = "Zp:au = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "1.")
{
$data = "Zp:vd = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "0.")
{
$data = "Zp:ad = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "1.")
{
$data = "Zp:vl = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "0.")
{
$data = "Zp:al = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "1.")
{
$data = "Zp:vv = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "0.")
{
$data = "Zp:av = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "0.1")
{
$data = "Zp:epsilon = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "0.1")
{
$data = "Zp:coupH = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "off")
{
$data = "DM:qqbar2DY = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "1")
{
$data = "DM:DYtype = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "")
{
$data = "DM:yuk1 = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "")
{
$data = "DM:yuk2 = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "")
{
$data = "DM:yuk3 = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "2")
{
$data = "DM:Nplet = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "")
{
$data = "DM:M1 = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "")
{
$data = "DM:M2 = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "")
{
$data = "DM:Lambda = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
