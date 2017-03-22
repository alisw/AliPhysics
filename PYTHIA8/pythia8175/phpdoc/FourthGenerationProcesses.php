<html>
<head>
<title>Fourth-Generation Processes</title>
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

<form method='post' action='FourthGenerationProcesses.php'>

<h2>Fourth-Generation Processes</h2>

A fourth generation can be accommodated within the Standard Model, 
without the introduction of any new concepts. Many experimental 
constraints exist, but it has not been fully excluded. Therefore 
we offer a simple implementation, along the lines of the top. 
It could also be useful as a template for studies of other 
new particles with similar characteristics. 

<p/>
The fourth generation are given names as in the third, but with a prime,
i.e. <i>b'</i> with PDG code 7, <i>t'</i> with code 8, 
<i>tau'</i> with code 17, and <i>nu'_tau</i> with code 18.
Most important for you is to assign a mass hierarchy, to decide which
fermions can decay into which. The current implementation assumes that 
mass splittings are big enough that fourth-generation fermions can
decay to third-generation ones by the emission of an on-shell <i>W</i>.
To this end, the standard three-generation CKM mixing matrix has been 
extended to include a fourth generation, see below. Since no mixing has
been implemented in the neutrino sector it would be assumed that the
lighter of <i>tau'</i> and <i>nu'_tau</i> is stable. No decay modes 
have been implemented that go beyond the Standard Model, so 
modifications would be needed if e.g. also SUSY is included in the game.

<h3>Production processes</h3>

<h4>1) <i>b'</i> processes</h4>

Different ways to produce <i>b'</i> quarks, singly or in pairs.
For a <i>b' t'</i> pair see section 3 below.

<br/><br/><strong>FourthBottom:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of <i>b'</i> production.
Also includes the process <i>f fbar' -> t' b'bar</i> in section 3 below.
  

<br/><br/><strong>FourthBottom:gg2bPrimebPrimebar</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> b' b'bar</i>. 
Code 801.
  

<br/><br/><strong>FourthBottom:qqbar2bPrimebPrimebar</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar -> b' b'bar</i> by gluon exchange. 
Code 802.
  

<br/><br/><strong>FourthBottom:qq2bPrimeq(t:W)</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q q' -> b' q''</i> by <i>t</i>-channel exchange 
of a <i>W^+-</i> boson. 
Code 803.
  

<br/><br/><strong>FourthBottom:ffbar2bPrimebPrimebar(s:gmZ)</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> b' b'bar</i> by <i>s</i>-channel exchange 
of a <i>gamma^*/Z^0</i> boson. 
Code 804.
  

<br/><br/><strong>FourthBottom:ffbar2bPrimeqbar(s:W)</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar' -> b' qbar''</i> by <i>s</i>-channel exchange 
of a <i>W^+-</i> boson. Here <i>q''</i> is either <i>u</i> or
<i>c</i>.
Code 805.
  

<br/><br/><strong>FourthBottom:ffbar2bPrimetbar(s:W)</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar' -> b' tbar</i> by <i>s</i>-channel exchange 
of a <i>W^+-</i> boson. 
Code 806.
  

<h4>2) <i>t'</i> processes</h4>

Different ways to produce <i>t'</i> quarks, singly or in pairs.
For a <i>b' t'</i> pair see section 3 below.

<br/><br/><strong>FourthTop:all</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of <i>t'</i> production.
Also includes the process <i>f fbar' -> t' b'bar</i> in section 3 below.
  

<br/><br/><strong>FourthTop:gg2tPrimetPrimebar</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> t' t'bar</i>. 
Code 821.
  

<br/><br/><strong>FourthTop:qqbar2tPrimetPrimebar</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar -> t' t'bar</i> by gluon exchange. 
Code 822.
  

<br/><br/><strong>FourthTop:qq2tPrimeq(t:W)</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q q' -> t' q''</i> by <i>t</i>-channel exchange 
of a <i>W^+-</i> boson. 
Code 823.
  

<br/><br/><strong>FourthTop:ffbar2tPrimetPrimebar(s:gmZ)</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> t' t'bar</i> by <i>s</i>-channel exchange 
of a <i>gamma^*/Z^0</i> boson. 
Code 824.
  

<br/><br/><strong>FourthTop:ffbar2tPrimeqbar(s:W)</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar' -> t' qbar''</i> by <i>s</i>-channel exchange 
of a <i>W^+-</i> boson.
Code 825.
  

<h4>3) Pair-processes with different flavours</h4>
 
Different ways to produce two different fourth-generation fermions.

<br/><br/><strong>FourthPair:ffbar2tPrimebPrimebar(s:W)</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar' -> t' b'bar</i> by <i>s</i>-channel exchange 
of a <i>W^+-</i> boson.
Code 841.
  

<br/><br/><strong>FourthPair:ffbar2tauPrimenuPrimebar(s:W)</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar' -> tau' nu'_taubar</i> by <i>s</i>-channel 
exchange of a <i>W^+-</i> boson.
Code 842.
  

<p/>
Missing in this list is scatterings <i>q q' -> t' b'</i> by 
<i>t</i>-channel exchange of a <i>W^+-</i> boson, since currently
the matrix element for such processes have not been implemented for
two massive particles in the final state. Since this process would
involve two CKM-suppressed vertices it ought to be small. 

<h3>Parameters</h3>

The Cabibbo-Kobayashi-Maskawa matrix is extended by seven further values.
So as not to mess up the Standard Model, the normal 3 * 3 matrix is
kept unitary, and so the new off-diagonal elements lead to a slight
breaking of this. For exploratory studies this should be good enough;
more detailed 4 * 4 tunes to data would only make sense the day there
are evidence for the existence of a fourth generation. 

<br/><br/><table><tr><td><strong>FourthGeneration:VubPrime </td><td></td><td> <input type="text" name="16" value="0.001" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.001</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
The <i>V_ub'</i> matrix element in the 4 * 4 CKM matrix.
  

<br/><br/><table><tr><td><strong>FourthGeneration:VcbPrime </td><td></td><td> <input type="text" name="17" value="0.01" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
The <i>V_cb'</i> matrix element in the 4 * 4 CKM matrix.
  

<br/><br/><table><tr><td><strong>FourthGeneration:VtbPrime </td><td></td><td> <input type="text" name="18" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
The <i>V_tb'</i> matrix element in the 4 * 4 CKM matrix.
  

<br/><br/><table><tr><td><strong>FourthGeneration:VtPrimed </td><td></td><td> <input type="text" name="19" value="0.001" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.001</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
The <i>V_t'd</i> matrix element in the 4 * 4 CKM matrix.
  

<br/><br/><table><tr><td><strong>FourthGeneration:VtPrimes </td><td></td><td> <input type="text" name="20" value="0.01" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
The <i>V_t's</i> matrix element in the 4 * 4 CKM matrix.
  

<br/><br/><table><tr><td><strong>FourthGeneration:VtPrimeb </td><td></td><td> <input type="text" name="21" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
The <i>V_t'b</i> matrix element in the 4 * 4 CKM matrix.
  

<br/><br/><table><tr><td><strong>FourthGeneration:VtPrimebPrime </td><td></td><td> <input type="text" name="22" value="0.99" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.99</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
The <i>V_t'b'</i> matrix element in the 4 * 4 CKM matrix.
   

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
$data = "FourthBottom:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "FourthBottom:gg2bPrimebPrimebar = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "FourthBottom:qqbar2bPrimebPrimebar = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "FourthBottom:qq2bPrimeq(t:W) = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "FourthBottom:ffbar2bPrimebPrimebar(s:gmZ) = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "FourthBottom:ffbar2bPrimeqbar(s:W) = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "FourthBottom:ffbar2bPrimetbar(s:W) = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "FourthTop:all = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "FourthTop:gg2tPrimetPrimebar = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "FourthTop:qqbar2tPrimetPrimebar = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "FourthTop:qq2tPrimeq(t:W) = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "FourthTop:ffbar2tPrimetPrimebar(s:gmZ) = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "FourthTop:ffbar2tPrimeqbar(s:W) = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "FourthPair:ffbar2tPrimebPrimebar(s:W) = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "FourthPair:ffbar2tauPrimenuPrimebar(s:W) = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "0.001")
{
$data = "FourthGeneration:VubPrime = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "0.01")
{
$data = "FourthGeneration:VcbPrime = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "0.1")
{
$data = "FourthGeneration:VtbPrime = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "0.001")
{
$data = "FourthGeneration:VtPrimed = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "0.01")
{
$data = "FourthGeneration:VtPrimes = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "0.1")
{
$data = "FourthGeneration:VtPrimeb = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "0.99")
{
$data = "FourthGeneration:VtPrimebPrime = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->

