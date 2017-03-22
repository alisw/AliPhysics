<html>
<head>
<title>Top Processes</title>
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

<form method='post' action='TopProcesses.php'>
 
<h2>Top Processes</h2> 
 
Different ways to produce top quarks, singly or in pairs. 
 
<br/><br/><strong>Top:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of top production. 
   
 
<br/><br/><strong>Top:gg2ttbar</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g &rarr; t tbar</i>. 
Code 601. 
   
 
<br/><br/><strong>Top:qqbar2ttbar</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar &rarr; t tbar</i> by gluon exchange. 
Code 602. 
   
 
<br/><br/><strong>Top:qq2tq(t:W)</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q q' &rarr; t q''</i> by <i>t</i>-channel exchange 
of a <i>W^+-</i> boson. 
Code 603. 
   
 
<br/><br/><strong>Top:ffbar2ttbar(s:gmZ)</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar &rarr; t tbar</i> by <i>s</i>-channel exchange 
of a <i>gamma^*/Z^0</i> boson. 
Code 604. 
   
 
<br/><br/><strong>Top:ffbar2tqbar(s:W)</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar' &rarr; t q''</i> by <i>s</i>-channel exchange 
of a <i>W^+-</i> boson. 
Code 605. 
   
 
<br/><br/><strong>Top:gmgm2ttbar</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>gamma gamma &rarr; t tbar</i>. 
Code 606. 
   
 
<p/> 
By default top always decays to a <i>W</i> and a down-type quark. 
It is possible to switch on the <i>t &rarr; H+ b</i> decay mode. 
Note that its partial width is calculated using the <i>tan(beta)</i> 
value stored in <code>HiggsHchg:tanBeta</code>, so that it can be used 
without having to read in a SUSY parameter file. For the <i>H+</i> to 
decay also <code>Higgs:useBSM = on</code> is necessary. 
 
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
$data = "Top:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "Top:gg2ttbar = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "Top:qqbar2ttbar = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "Top:qq2tq(t:W) = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "Top:ffbar2ttbar(s:gmZ) = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "Top:ffbar2tqbar(s:W) = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "Top:gmgm2ttbar = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
