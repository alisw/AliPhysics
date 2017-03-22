<html>
<head>
<title>Hadron Scattering</title>
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

<form method='post' action='HadronScattering.php'>
 
<h2>Hadron Scattering</h2> 
 
A simple hadron scattering model. It is intended to take into account 
that the overlap of multiple strings at low transverse dimensions 
is likely to lead to some collective effects, not unlike those 
observed in heavy-ion collisions, even if not quite as pronounced. 
Specifically, it is assumed that the hadrons produced can scatter 
against each other on the way out, before the fragmenting system 
has had time to expand enough that the hadrons get free. Thereby 
heavier particles are shifted to higher transverse momenta, at the 
expense of the lighter ones. 
 
<br/><b>Warning:</b> This is still at an experimental level, 
and should not be used unless you know what you are doing. 
 
<br/><br/><strong>HadronScatter:scatter</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Master flag for hadron scattering. 
   
 
<br/><br/><strong>HadronScatter:afterDecay</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Perform hadron scattering before or after first round of decays, 
involving very short-lived particles like the <i>rho</i>. 
The default is to perform scattering directly after the 
string fragmentation, before any decays. 
   
 
<br/><br/><strong>HadronScatter:allowDecayProd</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow two hadrons with same parent hadron to scatter. 
   
 
<br/><br/><strong>HadronScatter:scatterRepeat</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow hadrons which have already scattered to scatter again. 
Even if switched on, the same pair can not scatter off each 
other twice. 
   
 
<h3>Hadron selection</h3> 
 
<br/><br/><table><tr><td><strong>HadronScatter:hadronSelect  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 0</code>)</td></tr></table>
Probability that a hadron is soft enough to scatter. 
(A high-<ei>pT</ei> hadron presumably being part of a jet, 
and thus produced away from the high-particle-density region 
at small transverse dimensions.) 
<br/>
<input type="radio" name="5" value="0" checked="checked"><strong>0 </strong>:   <ei>P = N exp(-pT^2 / 2 / sigma^2) /    ( (1 - k) exp(-pT^2 / 2 / sigma^2) + k pT0^p / (pT0^2 + pT^2)^(p/2), </ei>  with <ei>sigma = 2 StringPT:sigma</ei> and <ei>pT0</ei> the same as that  used in <ei>MultipartonInteractions</ei>.  <br/>
 
<br/><br/><table><tr><td><strong>HadronScatter:N </td><td></td><td> <input type="text" name="6" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.01</code>; <code>maximum = 1.0</code>)</td></tr></table>
<i>N</i> parameter as above. 
   
<br/><br/><table><tr><td><strong>HadronScatter:k </td><td></td><td> <input type="text" name="7" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.01</code>; <code>maximum = 1.0</code>)</td></tr></table>
<i>k</i> parameter as above. 
   
<br/><br/><table><tr><td><strong>HadronScatter:p </td><td></td><td> <input type="text" name="8" value="6" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>6</strong></code>; <code>minimum = 2</code>; <code>maximum = 30</code>)</td></tr></table>
<i>p</i> parameter as above. 
   
 
<h3>Scattering probability</h3> 
 
<br/><br/><table><tr><td><strong>HadronScatter:scatterProb  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Probability for a pair of hadrons to scatter. 
<br/>
<input type="radio" name="9" value="0" checked="checked"><strong>0 </strong>: All hadrons scatter with probability  <ei>j max(0, 1 - dR^2 / rMax^2)</ei>. Angular distribution  is picked flat in <ei>cos(theta).</ei><br/>
<input type="radio" name="9" value="1"><strong>1 </strong>: As option 0, above, but only <ei>pi-pi</ei>,  <ei>pi-K</ei> and <ei>pi-p</ei> scatterings are considered.  <br/>
<input type="radio" name="9" value="2"><strong>2 </strong>: Only <ei>pi-pi</ei>, <ei>pi-K</ei> and  <ei>pi-p</ei> scatterings are considered, with probability  given by <ei>(1 - exp(-j sigEl)) max(0, 1 - dR^2 / rMax^2)</ei>.  The elastic cross sections and angular distributions are taken  from the partial-wave distributions.  <br/>
 
<br/><br/><table><tr><td><strong>HadronScatter:j </td><td></td><td> <input type="text" name="10" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 10.0</code>)</td></tr></table>
<i>j</i> parameter as above. 
   
<br/><br/><table><tr><td><strong>HadronScatter:rMax </td><td></td><td> <input type="text" name="11" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 2.0</code>)</td></tr></table>
<i>rMax</i> parameter as above. 
   
 
<br/><br/><strong>HadronScatter:tile</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Use tiling in <i>(eta, phi)</i> to reduce number of pairwise tests. 
   
 
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
$data = "HadronScatter:scatter = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "HadronScatter:afterDecay = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "HadronScatter:allowDecayProd = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "HadronScatter:scatterRepeat = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "0")
{
$data = "HadronScatter:hadronSelect = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "1.0")
{
$data = "HadronScatter:N = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "1.0")
{
$data = "HadronScatter:k = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "6")
{
$data = "HadronScatter:p = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "0")
{
$data = "HadronScatter:scatterProb = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "0.5")
{
$data = "HadronScatter:j = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "0.5")
{
$data = "HadronScatter:rMax = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "HadronScatter:tile = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
