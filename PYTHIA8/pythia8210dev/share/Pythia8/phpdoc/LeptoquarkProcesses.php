<html>
<head>
<title>Leptoquark Processes</title>
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

<form method='post' action='LeptoquarkProcesses.php'>
 
<h2>Leptoquark Processes</h2> 
 
Leptoquarks arise in many scenarios, and can have widely different 
characteristics, with respect to spin, isospin am d flavour. 
The current implementation in no sense attempts to exhaust these 
possibilities, but only to encode one of the simplest possibilities, 
with a single scalar leptoquark, denoted <i>LQ</i> and assigned PDG 
code 42. The leptoquark is assumed to carry specific quark 
and lepton quantum numbers, by default <i>u</i> quark plus electron. 
These flavour numbers are conserved, i.e. a process such as 
<i>u e^- &rarr; LQ &rarr; d nu_e</i> is not allowed. 
 
<p/> 
Although only one leptoquark is implemented, its flavours may be 
changed arbitrarily to study the different possibilities. The 
flavours of the leptoquark are defined by the quark and lepton 
flavours in the decay mode list. Therefore, to change from the 
current <i>u e^-</i> to <i>c mu^+</i>, say, you only need 
a line 
<br/><code>pythia.readString("42:0:products = 4 -13");</code> 
<br/>in your main program, or the equivalent in a command file. 
The former must always be a quark, while the latter could be a lepton 
or an antilepton; a charge-conjugate partner is automatically defined 
by the program. At initialization, the charge is recalculated as a 
function of the flavours defined; also the leptoquark name is redefined 
to be of the type <code>LQ_q,l</code>, where actual quark and lepton 
flavours are displayed. 
 
<p/> 
The leptoquark is likely to be fairly long-lived, in which case it 
could have time to fragment into a mesonic- or baryonic-type state, which 
would decay later on. Currently this possibility is not handled; therefore 
the leptoquark is always assumed to decay before fragmentation. 
For that reason the leptoquark can also not be put stable. 
 
<h3>Production processes</h3> 
 
Four production processes have been implemented, which normally would 
not overlap and therefore could be run together. 
 
<br/><br/><strong>LeptoQuark:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of lowest-order <i>LQ</i> production 
processes, i.e. the four ones below. 
   
 
<br/><br/><strong>LeptoQuark:ql2LQ</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q l &rarr; LQ</i>. 
Code 3201. 
   
 
<br/><br/><strong>LeptoQuark:qg2LQl</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q g &rarr; LQ l</i>. 
Code 3202. 
   
 
<br/><br/><strong>LeptoQuark:gg2LQLQbar</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g &rarr; LQ LQbar</i>. 
Code 3203. 
   
 
<br/><br/><strong>LeptoQuark:qqbar2LQLQbar</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar &rarr; LQ LQbar</i>. 
Code 3204. 
   
 
<h3>Parameters</h3> 
 
In the above scenario the main free parameters are the leptoquark flavour 
content, set as already described, and the <i>LQ</i> mass, set as usual. 
In addition there is one further parameter. 
 
<br/><br/><table><tr><td><strong>LeptoQuark:kCoup </td><td></td><td> <input type="text" name="6" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
multiplicative factor in the <i>LQ &rarr; q l</i> squared Yukawa coupling, 
and thereby in the <i>LQ</i> width and the <i>q l &rarr; LQ</i> and 
other cross sections. Specifically, <i>lambda^2/(4 pi) = k alpha_em</i>, 
i.e. it corresponds to the $k$ factor of [<a href="Bibliography.php" target="page">Hew88</a>]. 
   
 
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
$data = "LeptoQuark:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "LeptoQuark:ql2LQ = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "LeptoQuark:qg2LQl = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "LeptoQuark:gg2LQLQbar = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "LeptoQuark:qqbar2LQLQbar = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "1.0")
{
$data = "LeptoQuark:kCoup = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
