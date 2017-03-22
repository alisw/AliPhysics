<html>
<head>
<title>Total Cross Sections</title>
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

<form method='post' action='TotalCrossSections.php'>
 
<h2>Total Cross Sections</h2> 
 
The <code>SigmaTotal</code> class returns the total, elastic, diffractive 
and nondiffractive cross sections in hadronic collisions, and also the 
slopes of the <i>d(sigma)/dt</i> distributions. Most of the parametrizations 
used are from [<a href="Bibliography.php" target="page">Sch94, Sch97</a>] which borrows some of the total cross 
sections from [<a href="Bibliography.php" target="page">Don92</a>]. If you use the MBR (Minimum Bias Rockefeller) 
model [<a href="Bibliography.php" target="page">Cie12</a>], <code>Diffraction:PomFlux = 5</code>, this model 
contains its own parametrizations of all cross sections in <i>p p</i> 
and <i>pbar p</i> collisions. 
 
<p/> 
There are strong indications that the currently implemented diffractive 
cross section parametrizations, which should be in reasonable agreement 
with data at lower energies, overestimate the diffractive rate at larger 
values. If you wish to explore this (or other) aspect, it is possible to 
override the cross section values in two different ways. The first offers 
(almost) complete freedom, but needs to be defined separately for each 
CM energy, while the second introduces a simpler parametrized damping. 
The two cannot be combined. Furthermore the Coulomb term for elastic 
scattering, which by default is off, can be switched on. 
 
<p/> 
The allowed combinations of incoming particles are <i>p + p</i>, 
<i>pbar + p</i>, <i>pi+ + p</i>, <i>pi- + p</i>, 
<i>pi0/rho0 + p</i>, <i>phi + p</i>, <i>J/psi + p</i>, 
<i>rho + rho</i>, <i>rho + phi</i>, <i>rho + J/psi</i>, 
<i>phi + phi</i>, <i>phi + J/psi</i>, <i>J/psi + J/psi</i>. 
The strong emphasis on vector mesons is related to the description 
of <i>gamma + p</i> and <i>gamma + gamma</i> interactions in a 
Vector Dominance Model framework (which will not be available for some 
time to come, so this is a bit of overkill). Nevertheless, the sections 
below, with allowed variations, are mainly intended to make sense for 
<i>p + p</i>. 
 
<h3>Central diffraction</h3> 
 
Central diffraction (CD), a.k.a. double Pomeron exchange (DPE), was not 
part of the framework in [<a href="Bibliography.php" target="page">Sch94</a>]. It has now been added for 
multiparticle states, i.e. excluding the resonance region below 1 GeV 
mass, as well as other exclusive states, but only for <i>p p</i> or 
<i>pbar p</i>. It uses the same proton-Pomeron vertex as in single 
diffraction, twice, to describe <i>x_Pomeron</i> and <i>t</i> spectra. 
This fixes the energy dependence, which has been integrated and 
parametrized. The absolute normalization has been left open, however. 
Furthermore, since CD has not been included in previous tunes to data, 
a special flag is available to reproduce the old behaviour (with due 
complications when one does not want to do this). 
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaAXB2TeV </td><td></td><td> <input type="text" name="1" value="1.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.5</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The CD cross section for <i>p p</i> and <i>pbar p</i> collisions, 
normalized to its value at 2 TeV CM energy, expressed in mb. The energy 
dependence is then parametrized, and behaves roughly like 
<i>ln^1.5(s)</i>. Is used for the options 
<code>Diffraction:PomFlux = 1 - 4</code>, while the MBR model 
(<code>= 5</code>) has its own parametrization. 
   
 
<br/><br/><strong>SigmaTotal:zeroAXB</strong>  <input type="radio" name="2" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="2" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
several existing <?php $filepath = $_GET["filepath"];
echo "<a href='Tunes.php?filepath=".$filepath."' target='page'>";?>tunes</a> do not include CD. 
An inclusion of a nonvanishing CD cross section directly affects 
the nondiffractive phenomenology, even if not dramatically, and so 
this flag is used to forcibly set the CD cross section to vanish 
in such tunes. You can switch CD back on <i>after</i> the selection of 
a tune, if you so wish, by resetting <code>SigmaTotal:zeroAXB = off</code>. 
This option has no effect for the MBR model 
(<code>Diffraction:PomFlux = 5</code>), where the CD cross section 
has been included from the onset. 
   
 
<h3>Set cross sections</h3> 
 
<br/><br/><strong>SigmaTotal:setOwn</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow a user to set own cross sections by hand; on/off = true/false. 
   
 
<p/> 
When <code>SigmaTotal:setOwn = on</code>, the user is expected to set 
values for the corresponding cross sections: 
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaTot </td><td></td><td> <input type="text" name="4" value="80." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>80.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Total cross section in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaEl </td><td></td><td> <input type="text" name="5" value="20." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>20.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Elastic cross section in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaXB </td><td></td><td> <input type="text" name="6" value="8." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>8.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Single Diffractive cross section <i>A + B &rarr; X + B</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaAX </td><td></td><td> <input type="text" name="7" value="8." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>8.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Single Diffractive cross section <i>A + B &rarr; A + X</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaXX </td><td></td><td> <input type="text" name="8" value="4." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Double Diffractive cross section <i>A + B &rarr; X_1 + X_2</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaAXB </td><td></td><td> <input type="text" name="9" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Central Diffractive cross section <i>A + B &rarr; A + X + B</i> in mb. 
   
 
<p/> 
Note that the total cross section subtracted by the elastic and various 
diffractive ones gives the inelastic nondiffractive cross section, 
which therefore is not set separately. If this cross section evaluates 
to be negative the internal parametrizations are used instead of the 
ones here. However, since the nondiffractive inelastic cross section 
is what makes up the minimum-bias event class, and plays a major role 
in the description of multiparton interactions, it is important that a 
consistent set is used. 
 
<h3>Dampen diffractive cross sections</h3> 
 
As already noted, unitarization effects may dampen the rise of diffractive 
cross sections relative to the default parametrizations. The settings 
here allows one way to introduce a dampening, which is used in some 
of the existing <?php $filepath = $_GET["filepath"];
echo "<a href='Tunes.php?filepath=".$filepath."' target='page'>";?>tunes</a>. 
 
<br/><br/><strong>SigmaDiffractive:dampen</strong>  <input type="radio" name="10" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="10" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow a user to dampen diffractive cross sections; on/off = true/false. 
   
 
<p/> 
When <code>SigmaDiffractive:dampen = on</code>, the three diffractive 
cross sections are damped so that they never can exceed the respective 
values below. Specifically, if the standard parametrization gives 
the cross section <i>sigma_old(s)</i> and a fixed <i>sigma_max</i> 
is set, the actual cross section becomes <i>sigma_new(s) 
= sigma_old(s) * sigma_max / (sigma_old(s) + sigma_max)</i>. 
This reduces to <i>sigma_old(s)</i> at low energies and to 
<i>sigma_max</i> at high ones. Note that the asymptotic value 
is approached quite slowly, however. 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:maxXB </td><td></td><td> <input type="text" name="11" value="65." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>65.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The above <i>sigma_max</i> for <i>A + B &rarr; X + B</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:maxAX </td><td></td><td> <input type="text" name="12" value="65." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>65.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The above <i>sigma_max</i> for <i>A + B &rarr; A + X</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:maxXX </td><td></td><td> <input type="text" name="13" value="65." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>65.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The above <i>sigma_max</i> for <i>A + B &rarr; X_1 + X_2</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:maxAXB </td><td></td><td> <input type="text" name="14" value="3." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The above <i>sigma_max</i> for <i>A + B &rarr; A + X + B</i> in mb. 
   
 
<p/> 
As above, a reduced diffractive cross section automatically translates 
into an increased nondiffractive one, such that the total (and elastic) 
cross section remains fixed. 
 
 
<h3>Set elastic cross section</h3> 
 
<p/> 
In the above option the <i>t</i> slopes are based on the internal 
parametrizations. In addition there is no Coulomb-term contribution 
to the elastic (or total) cross section, which of course becomes 
infinite if this contribution is included. If you have switched on 
<code>SigmaTotal:setOwn</code> you can further switch on a machinery 
to include the Coulomb term, including interference with the conventional 
strong-interaction Pomeron one [<a href="Bibliography.php" target="page">Ber87</a>]. Then the elastic cross 
section is no longer taken from <code>SigmaTotal:sigmaEl</code> but 
derived from the parameters below and <code>SigmaTotal:sigmaTot</code>, 
using the optical theorem. The machinery is only intended to be used for 
<i>p p</i> and <i>pbar p</i> collisions. The description of 
diffractive events, and especially their slopes, remains unchanged. 
 
<br/><br/><strong>SigmaElastic:setOwn</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>no</strong></code>)<br/>
Allow a user to set parameters for the normalization and shape of the 
elastic cross section the by hand; yes/no = true/false. 
   
 
<br/><br/><table><tr><td><strong>SigmaElastic:bSlope </td><td></td><td> <input type="text" name="16" value="18." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>18.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
the slope <i>b</i> of the strong-interaction term <i>exp(bt)</i>, 
in units of GeV^-2. 
   
 
<br/><br/><table><tr><td><strong>SigmaElastic:rho </td><td></td><td> <input type="text" name="17" value="0.13" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.13</strong></code>; <code>minimum = -1.</code>; <code>maximum = 1.</code>)</td></tr></table>
the ratio of the real to the imaginary parts of the nuclear scattering 
amplitude. 
   
 
<br/><br/><table><tr><td><strong>SigmaElastic:lambda </td><td></td><td> <input type="text" name="18" value="0.71" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.71</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 2.</code>)</td></tr></table>
the main parameter of the electric form factor 
<i>G(t) = lambda^2 / (lambda + |t|)^2</i>, in units of GeV^2. 
   
 
<br/><br/><table><tr><td><strong>SigmaElastic:tAbsMin </td><td></td><td> <input type="text" name="19" value="5e-5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5e-5</strong></code>; <code>minimum = 1e-10</code>)</td></tr></table>
since the Coulomb contribution is infinite a lower limit on 
<i>|t|</i> must be set to regularize the divergence, 
in units of GeV^2. 
   
 
<br/><br/><table><tr><td><strong>SigmaElastic:phaseConst </td><td></td><td> <input type="text" name="20" value="0.577" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.577</strong></code>)</td></tr></table>
The Coulomb term is taken to contain a phase factor 
<i>exp(+- i alpha phi(t))</i>, with + for <i>p p</i> and - for 
<i>pbar p</i>, where <i>phi(t) = - phaseConst - ln(-B t/2)</i>. 
This constant is model dependent [<a href="Bibliography.php" target="page">Cah82</a>]. 
   
 
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

if($_POST["1"] != "1.5")
{
$data = "SigmaTotal:sigmaAXB2TeV = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "on")
{
$data = "SigmaTotal:zeroAXB = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "SigmaTotal:setOwn = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "80.")
{
$data = "SigmaTotal:sigmaTot = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "20.")
{
$data = "SigmaTotal:sigmaEl = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "8.")
{
$data = "SigmaTotal:sigmaXB = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "8.")
{
$data = "SigmaTotal:sigmaAX = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "4.")
{
$data = "SigmaTotal:sigmaXX = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "1.")
{
$data = "SigmaTotal:sigmaAXB = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "on")
{
$data = "SigmaDiffractive:dampen = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "65.")
{
$data = "SigmaDiffractive:maxXB = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "65.")
{
$data = "SigmaDiffractive:maxAX = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "65.")
{
$data = "SigmaDiffractive:maxXX = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "3.")
{
$data = "SigmaDiffractive:maxAXB = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "no")
{
$data = "SigmaElastic:setOwn = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "18.")
{
$data = "SigmaElastic:bSlope = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "0.13")
{
$data = "SigmaElastic:rho = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "0.71")
{
$data = "SigmaElastic:lambda = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "5e-5")
{
$data = "SigmaElastic:tAbsMin = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "0.577")
{
$data = "SigmaElastic:phaseConst = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
