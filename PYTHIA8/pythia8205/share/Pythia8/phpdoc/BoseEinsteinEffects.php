<html>
<head>
<title>Bose-Einstein Effects</title>
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

<form method='post' action='BoseEinsteinEffects.php'>
 
<h2>Bose-Einstein Effects</h2> 
 
The <code>BoseEinstein</code> class performs shifts of momenta 
of identical particles to provide a crude estimate of 
Bose-Einstein effects. The algorithm is the BE_32 one described in 
[<a href="Bibliography.php" target="page">Lon95</a>], with a Gaussian parametrization of the enhancement. 
We emphasize that this approach is not based on any first-principles 
quantum mechanical description of interference phenomena; such 
approaches anyway have many problems to contend with. Instead a cruder 
but more robust approach is adopted, wherein BE effects are introduced 
after the event has already been generated, with the exception of the 
decays of long-lived particles. The trick is that momenta of identical 
particles are shifted relative to each other so as to provide an 
enhancement of pairs closely separated, which is compensated by a 
depletion of pairs in an intermediate region of separation. 
 
<p/> 
More precisely, the intended target form of the BE correlations in 
BE_32 is 
<br/><i> 
f_2(Q) = (1 + lambda * exp(-Q^2 R^2)) 
       * (1 + alpha * lambda * exp(-Q^2 R^2/9) * (1 - exp(-Q^2 R^2/4))) 
</i><br/> 
where <i>Q^2 = (p_1 + p_2)^2 - (m_1 + m_2)^2</i>. 
Here the strength <i>lambda</i> and effective radius <i>R</i> 
are the two main parameters. The first factor of the 
equation is implemented by pulling pairs of identical hadrons closer 
to each other. This is done in such a way that three-momentum is 
conserved, but at the price of a small but non-negligible negative 
shift in the energy of the event. The second factor compensates this 
by pushing particles apart. The negative <i>alpha</i> parameter is 
determined iteratively, separately for each event, so as to restore 
energy conservation. The effective radius parameter is here <i>R/3</i>, 
i.e. effects extend further out in <i>Q</i>. Without the dampening 
<i>(1 - exp(-Q^2 R^2/4))</i> in the second factor the value at the 
origin would become <i>f_2(0) = (1 + lambda) * (1 + alpha * lambda)</i>, 
with it the desired value <i>f_2(0) = (1 + lambda)</i> is restored. 
The end result can be viewed as a poor man's rendering of a rapidly 
dampened oscillatory behaviour in <i>Q</i>. 
 
<p/> 
Further details can be found in [<a href="Bibliography.php" target="page">Lon95</a>]. For instance, the 
target is implemented under the assumption that the initial distribution 
in <i>Q</i> can be well approximated by pure phase space at small 
values, and implicitly generates higher-order effects by the way 
the algorithm is implemented. The algorithm is applied after the decay 
of short-lived resonances such as the <i>rho</i>, but before the decay 
of longer-lived particles. 
 
<p/> 
This algorithm is known to do a reasonable job of describing BE 
phenomena at LEP. It has not been tested against data for hadron 
colliders, to the best of our knowledge, so one should exercise some 
judgment before using it. Therefore by default the master switch 
<?php $filepath = $_GET["filepath"];
echo "<a href='MasterSwitches.php?filepath=".$filepath."' target='page'>";?>HadronLevel:BoseEinstein</a> is off. 
Furthermore, the implementation found here is not (yet) as 
sophisticated as the one used at LEP2, in that no provision is made 
for particles from separate colour singlet systems, such as 
<i>W</i>'s and <i>Z</i>'s, interfering only at a reduced rate. 
 
<p/> 
<b>Warning:</b> The algorithm will create a new copy of each particle 
with shifted momentum by BE effects, with status code 99, while the 
original particle with the original momentum at the same time will be 
marked as decayed. This means that if you e.g. search for all 
<i>pi+-</i> in an event you will often obtain the same particle twice. 
One way to protect yourself from unwanted doublecounting is to 
use only particles with a positive status code, i.e. ones for which 
<code>event[i].isFinal()</code> is <code>true</code>. 
 
 
<h3>Main parameters</h3> 
 
Assuming you have set <code>HadronLevel:BoseEinstein = on</code>, 
you can regulate the detailed behaviour with the following settings. 
 
<br/><br/><strong>BoseEinstein:Pion</strong>  <input type="radio" name="1" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="1" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Include effects or not for identical <i>pi^+</i>, <i>pi^-</i> 
and <i>pi^0</i>. 
   
 
<br/><br/><strong>BoseEinstein:Kaon</strong>  <input type="radio" name="2" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="2" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Include effects or not for identical <i>K^+</i>, <i>K^-</i>, 
<i>K_S^0</i> and <i>K_L^0</i>. 
   
 
<br/><br/><strong>BoseEinstein:Eta</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Include effects or not for identical <i>eta</i> and <i>eta'</i>. 
   
 
<br/><br/><table><tr><td><strong>BoseEinstein:lambda </td><td></td><td> <input type="text" name="4" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 2.</code>)</td></tr></table>
The strength parameter for Bose-Einstein effects. On physical grounds 
it should not be above unity, but imperfections in the formalism 
used may require that nevertheless. 
   
 
<br/><br/><table><tr><td><strong>BoseEinstein:QRef </td><td></td><td> <input type="text" name="5" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.05</code>; <code>maximum = 1.</code>)</td></tr></table>
The size parameter of the region in <i>Q</i> space over which 
Bose-Einstein effects are significant.  Can be thought of as 
the inverse of an effective distance in normal space, 
<i>R = hbar / QRef</i>, with <i>R</i> as used in the above equation. 
That is, <i>f_2(Q) = (1 + lambda * exp(-(Q/QRef)^2)) * (...)</i>. 
   
 
<br/><br/><table><tr><td><strong>BoseEinstein:widthSep </td><td></td><td> <input type="text" name="6" value="0.02" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.02</strong></code>; <code>minimum = 0.001</code>; <code>maximum = 1.</code>)</td></tr></table>
Particle species with a width above this value (in GeV) are assumed 
to be so short-lived that they decay before Bose-Einstein effects 
are considered, while otherwise they do not. In the former case the 
decay products thus can obtain shifted momenta, in the latter not. 
The default has been picked such that both <i>rho</i> and 
<i>K^*</i> decay products would be modified. 
   
 
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
$data = "BoseEinstein:Pion = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "on")
{
$data = "BoseEinstein:Kaon = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "BoseEinstein:Eta = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "1.")
{
$data = "BoseEinstein:lambda = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "0.2")
{
$data = "BoseEinstein:QRef = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.02")
{
$data = "BoseEinstein:widthSep = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
 
