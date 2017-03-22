<html>
<head>
<title>Onia Processes</title>
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

<form method='post' action='OniaProcesses.php'>

<h2>Onia Processes</h2>

Production of J/psi or Upsilon, directly and via chi states and the 
colour-octet mechanism. 
In each process the square-bracketed expression specifies the state 
in spectroscopic notation, <i>(2S+1) L J</i>, followed by 
<i>(1)</i> for colour-singlet states and <i>(8)</i> for 
colour-octet ditto. 

<p/>
The original Fortran code for these processes has been contributed 
by Stefan Wolf [unpublished]. For the C++ version only the unpolarized
expressions are retained, since the theoretical predictions of the 
colour-octet model anyway do not agree with the experimental 
observations. Furthermore, the polarization effects are modest,
so isotropic decay is not a bad starting point. Such an event sample
can afterwards be reweighted at will by the user, to test various
assumptions.

<p/>
The description of  
<?php $filepath = $_GET["filepath"];
echo "<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>";?>final-state radiation</a>
is in this case based on some further model assumptions.

<p/>
Most of the processes below are divergent in the limit <i>pT -> 0</i>, 
and therefore a <i>pTmin</i> scale should be set. Comparisons with 
data indicate that this divergence can be tamed the same way as for 
the normal QCD <i>2 -> 2</i> cross sections [<a href="Bibliography.php" target="page">Bar06,Kra08</a>], 
which makes sense, since they are all dominated by the same kind of 
<i>t</i>-channel gluon exchange. It is therefore possible to use the 
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>SuppressSmallPT</a> user hook to impose a 
reweighting that cancels the low-<i>pT</i> divergence.

<p/>
An eikonalized description of these processes is included in the 
multiparton-interactions framework. Here the low-<i>pT</i> dampening
is automatic, and additionally the framework is more consistent
(e.g. with respect to energy-momentum constraints and the 
impact-parameter description) for events where the onium production 
is not the hardest subprocess, as would often be the case in the 
low-<i>pT</i> limit.

<h3>Charmonium</h3>

<br/><br/><strong>Charmonium:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of charmonium production.
  

<br/><br/><strong>Charmonium:gg2QQbar[3S1(1)]g</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> ccbar[3S1(1)] g</i>.
Code 401.
  

<br/><br/><strong>Charmonium:gg2QQbar[3P0(1)]g</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> ccbar[3P0(1)] g</i>.
Code 402.
  

<br/><br/><strong>Charmonium:gg2QQbar[3P1(1)]g</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> ccbar[3P1(1)] g</i>.
Code 403.
  

<br/><br/><strong>Charmonium:gg2QQbar[3P2(1)]g</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> ccbar[3P2(1)] g</i>.
Code 404.
  

<br/><br/><strong>Charmonium:qg2QQbar[3P0(1)]q</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> ccbar[3P0(1)] q</i>.
Code 405.
  

<br/><br/><strong>Charmonium:qg2QQbar[3P1(1)]q</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> ccbar[3P1(1)] q</i>.
Code 406.
  

<br/><br/><strong>Charmonium:qg2QQbar[3P2(1)]q</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> ccbar[3P2(1)] q</i>.
Code 407.
  

<br/><br/><strong>Charmonium:qqbar2QQbar[3P0(1)]g</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> ccbar[3P0(1)] g</i>.
Code 408.
  

<br/><br/><strong>Charmonium:qqbar2QQbar[3P1(1)]g</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> ccbar[3P1(1)] g</i>.
Code 409.
  

<br/><br/><strong>Charmonium:qqbar2QQbar[3P2(1)]g</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> ccbar[3P2(1)] g</i>.
Code 410.
  

<br/><br/><strong>Charmonium:gg2QQbar[3S1(8)]g</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> ccbar[3S1(8)] g</i>.
Code 411.
  

<br/><br/><strong>Charmonium:gg2QQbar[1S0(8)]g</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> ccbar[3S1(8)] g</i>.
Code 412.
  

<br/><br/><strong>Charmonium:gg2QQbar[3PJ(8)]g</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> ccbar[3S1(8)] g</i>.
Code 413.
  

<br/><br/><strong>Charmonium:qg2QQbar[3S1(8)]q</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> ccbar[3S1(8)] q</i>.
Code 414.
  

<br/><br/><strong>Charmonium:qg2QQbar[1S0(8)]q</strong>  <input type="radio" name="16" value="on"><strong>On</strong>
<input type="radio" name="16" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> ccbar[3S1(8)] q</i>.
Code 415.
  

<br/><br/><strong>Charmonium:qg2QQbar[3PJ(8)]q</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> ccbar[3S1(8)] q</i>.
Code 416.
  

<br/><br/><strong>Charmonium:qqbar2QQbar[3S1(8)]g</strong>  <input type="radio" name="18" value="on"><strong>On</strong>
<input type="radio" name="18" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> ccbar[3S1(8)] g</i>.
Code 417.
  

<br/><br/><strong>Charmonium:qqbar2QQbar[1S0(8)]g</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> ccbar[3S1(8)] g</i>.
Code 418.
  

<br/><br/><strong>Charmonium:qqbar2QQbar[3PJ(8)]g</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> ccbar[3S1(8)] g</i>.
Code 419.
  

<h3>Bottomonium</h3>

<br/><br/><strong>Bottomonium:all</strong>  <input type="radio" name="21" value="on"><strong>On</strong>
<input type="radio" name="21" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of charmonium production.
  

<br/><br/><strong>Bottomonium:gg2QQbar[3S1(1)]g</strong>  <input type="radio" name="22" value="on"><strong>On</strong>
<input type="radio" name="22" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> bbbar[3S1(1)] g</i>.
Code 501.
  

<br/><br/><strong>Bottomonium:gg2QQbar[3P0(1)]g</strong>  <input type="radio" name="23" value="on"><strong>On</strong>
<input type="radio" name="23" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> bbbar[3P0(1)] g</i>.
Code 502.
  

<br/><br/><strong>Bottomonium:gg2QQbar[3P1(1)]g</strong>  <input type="radio" name="24" value="on"><strong>On</strong>
<input type="radio" name="24" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> bbbar[3P1(1)] g</i>.
Code 503.
  

<br/><br/><strong>Bottomonium:gg2QQbar[3P2(1)]g</strong>  <input type="radio" name="25" value="on"><strong>On</strong>
<input type="radio" name="25" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> bbbar[3P2(1)] g</i>.
Code 504.
  

<br/><br/><strong>Bottomonium:qg2QQbar[3P0(1)]q</strong>  <input type="radio" name="26" value="on"><strong>On</strong>
<input type="radio" name="26" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> bbbar[3P0(1)] q</i>.
Code 505.
  

<br/><br/><strong>Bottomonium:qg2QQbar[3P1(1)]q</strong>  <input type="radio" name="27" value="on"><strong>On</strong>
<input type="radio" name="27" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> bbbar[3P1(1)] q</i>.
Code 506.
  

<br/><br/><strong>Bottomonium:qg2QQbar[3P2(1)]q</strong>  <input type="radio" name="28" value="on"><strong>On</strong>
<input type="radio" name="28" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> bbbar[3P2(1)] q</i>.
Code 507.
  

<br/><br/><strong>Bottomonium:qqbar2QQbar[3P0(1)]g</strong>  <input type="radio" name="29" value="on"><strong>On</strong>
<input type="radio" name="29" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> bbbar[3P0(1)] g</i>.
Code 508.
  

<br/><br/><strong>Bottomonium:qqbar2QQbar[3P1(1)]g</strong>  <input type="radio" name="30" value="on"><strong>On</strong>
<input type="radio" name="30" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> bbbar[3P1(1)] g</i>.
Code 509.
  

<br/><br/><strong>Bottomonium:qqbar2QQbar[3P2(1)]g</strong>  <input type="radio" name="31" value="on"><strong>On</strong>
<input type="radio" name="31" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> bbbar[3P2(1)] g</i>.
Code 510.
  

<br/><br/><strong>Bottomonium:gg2QQbar[3S1(8)]g</strong>  <input type="radio" name="32" value="on"><strong>On</strong>
<input type="radio" name="32" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> bbbar[3S1(8)] g</i>.
Code 511.
  

<br/><br/><strong>Bottomonium:gg2QQbar[1S0(8)]g</strong>  <input type="radio" name="33" value="on"><strong>On</strong>
<input type="radio" name="33" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> bbbar[3S1(8)] g</i>.
Code 512.
  

<br/><br/><strong>Bottomonium:gg2QQbar[3PJ(8)]g</strong>  <input type="radio" name="34" value="on"><strong>On</strong>
<input type="radio" name="34" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>g g -> bbbar[3S1(8)] g</i>.
Code 513.
  

<br/><br/><strong>Bottomonium:qg2QQbar[3S1(8)]q</strong>  <input type="radio" name="35" value="on"><strong>On</strong>
<input type="radio" name="35" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> bbbar[3S1(8)] q</i>.
Code 514.
  

<br/><br/><strong>Bottomonium:qg2QQbar[1S0(8)]q</strong>  <input type="radio" name="36" value="on"><strong>On</strong>
<input type="radio" name="36" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> bbbar[3S1(8)] q</i>.
Code 515.
  

<br/><br/><strong>Bottomonium:qg2QQbar[3PJ(8)]q</strong>  <input type="radio" name="37" value="on"><strong>On</strong>
<input type="radio" name="37" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q g -> bbbar[3S1(8)] q</i>.
Code 516.
  

<br/><br/><strong>Bottomonium:qqbar2QQbar[3S1(8)]g</strong>  <input type="radio" name="38" value="on"><strong>On</strong>
<input type="radio" name="38" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> bbbar[3S1(8)] g</i>.
Code 517.
  

<br/><br/><strong>Bottomonium:qqbar2QQbar[1S0(8)]g</strong>  <input type="radio" name="39" value="on"><strong>On</strong>
<input type="radio" name="39" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> bbbar[3S1(8)] g</i>.
Code 518.
  

<br/><br/><strong>Bottomonium:qqbar2QQbar[3PJ(8)]g</strong>  <input type="radio" name="40" value="on"><strong>On</strong>
<input type="radio" name="40" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i>q qbar -> bbbar[3S1(8)] g</i>.
Code 519.
  

<h3>Onium matrix elements</h3>

The implementation of charmonium and bottomonium production, including
the colour-octet production mechanism, requires information on NRQCD
matrix elements for the various wavefunctions involved. Default values
for these are encoded in the following ten variables. They
are taken from [<a href="Bibliography.php" target="page">Nas00</a>]; see also [<a href="Bibliography.php" target="page">Bar06</a>]. 

<br/><br/><table><tr><td><strong>Charmonium:OJpsi3S11 </td><td></td><td> <input type="text" name="41" value="1.16" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.16</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
<i>&lt;O(J/psi)[3S1(1)]&gt;</i>.
  

<br/><br/><table><tr><td><strong>Charmonium:OJpsi3S18 </td><td></td><td> <input type="text" name="42" value="0.0119" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0119</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
<i>&lt;O(J/psi)[3S1(8)]&gt;</i>.
  

<br/><br/><table><tr><td><strong>Charmonium:OJpsi1S08 </td><td></td><td> <input type="text" name="43" value="0.01" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
<i>&lt;O(J/psi)[1S0(8)]&gt;</i>.
  

<br/><br/><table><tr><td><strong>Charmonium:OJpsi3P08 </td><td></td><td> <input type="text" name="44" value="0.01" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
<i>&lt;O(J/psi)[3P0(8)]&gt;/m_c^2</i>.
  

<br/><br/><table><tr><td><strong>Charmonium:Ochic03P01 </td><td></td><td> <input type="text" name="45" value="0.05" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.05</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
<i>&lt;O(chi_c0)[3P0(8)]&gt;/m_c^2</i>.
  

<br/><br/><table><tr><td><strong>Bottomonium:OUpsilon3S11 </td><td></td><td> <input type="text" name="46" value="9.28" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>9.28</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
<i>&lt;O(Upsilon)[3S1(1)]&gt;</i>.
  

<br/><br/><table><tr><td><strong>Bottomonium:OUpsilon3S18 </td><td></td><td> <input type="text" name="47" value="0.15" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.15</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
<i>&lt;O(Upsilon)[3S1(8)]&gt;</i>.
  

<br/><br/><table><tr><td><strong>Bottomonium:OUpsilon1S08 </td><td></td><td> <input type="text" name="48" value="0.02" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.02</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
<i>&lt;O(Upsilon)[1S0(8)]&gt;</i>.
  

<br/><br/><table><tr><td><strong>Bottomonium:OUpsilon3P08 </td><td></td><td> <input type="text" name="49" value="0.02" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.02</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
<i>&lt;O(Upsilon)[3P0(8)]&gt;/m_b^2</i>.
  

<br/><br/><table><tr><td><strong>Bottomonium:Ochib03P01 </td><td></td><td> <input type="text" name="50" value="0.085" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.085</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
<i>&lt;O(chi_b0)[3P0(8)]&gt;/m_b^2</i>.
  


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
$data = "Charmonium:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "Charmonium:gg2QQbar[3S1(1)]g = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "Charmonium:gg2QQbar[3P0(1)]g = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "Charmonium:gg2QQbar[3P1(1)]g = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "Charmonium:gg2QQbar[3P2(1)]g = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "Charmonium:qg2QQbar[3P0(1)]q = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "Charmonium:qg2QQbar[3P1(1)]q = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "Charmonium:qg2QQbar[3P2(1)]q = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "Charmonium:qqbar2QQbar[3P0(1)]g = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "Charmonium:qqbar2QQbar[3P1(1)]g = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "Charmonium:qqbar2QQbar[3P2(1)]g = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "Charmonium:gg2QQbar[3S1(8)]g = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "Charmonium:gg2QQbar[1S0(8)]g = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "Charmonium:gg2QQbar[3PJ(8)]g = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "Charmonium:qg2QQbar[3S1(8)]q = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off")
{
$data = "Charmonium:qg2QQbar[1S0(8)]q = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "Charmonium:qg2QQbar[3PJ(8)]q = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "off")
{
$data = "Charmonium:qqbar2QQbar[3S1(8)]g = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "Charmonium:qqbar2QQbar[1S0(8)]g = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "Charmonium:qqbar2QQbar[3PJ(8)]g = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "off")
{
$data = "Bottomonium:all = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "off")
{
$data = "Bottomonium:gg2QQbar[3S1(1)]g = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "off")
{
$data = "Bottomonium:gg2QQbar[3P0(1)]g = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "off")
{
$data = "Bottomonium:gg2QQbar[3P1(1)]g = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "off")
{
$data = "Bottomonium:gg2QQbar[3P2(1)]g = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "off")
{
$data = "Bottomonium:qg2QQbar[3P0(1)]q = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "off")
{
$data = "Bottomonium:qg2QQbar[3P1(1)]q = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "off")
{
$data = "Bottomonium:qg2QQbar[3P2(1)]q = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "off")
{
$data = "Bottomonium:qqbar2QQbar[3P0(1)]g = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "off")
{
$data = "Bottomonium:qqbar2QQbar[3P1(1)]g = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "off")
{
$data = "Bottomonium:qqbar2QQbar[3P2(1)]g = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "off")
{
$data = "Bottomonium:gg2QQbar[3S1(8)]g = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "off")
{
$data = "Bottomonium:gg2QQbar[1S0(8)]g = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "off")
{
$data = "Bottomonium:gg2QQbar[3PJ(8)]g = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "off")
{
$data = "Bottomonium:qg2QQbar[3S1(8)]q = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "off")
{
$data = "Bottomonium:qg2QQbar[1S0(8)]q = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "off")
{
$data = "Bottomonium:qg2QQbar[3PJ(8)]q = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "off")
{
$data = "Bottomonium:qqbar2QQbar[3S1(8)]g = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "off")
{
$data = "Bottomonium:qqbar2QQbar[1S0(8)]g = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
if($_POST["40"] != "off")
{
$data = "Bottomonium:qqbar2QQbar[3PJ(8)]g = ".$_POST["40"]."\n";
fwrite($handle,$data);
}
if($_POST["41"] != "1.16")
{
$data = "Charmonium:OJpsi3S11 = ".$_POST["41"]."\n";
fwrite($handle,$data);
}
if($_POST["42"] != "0.0119")
{
$data = "Charmonium:OJpsi3S18 = ".$_POST["42"]."\n";
fwrite($handle,$data);
}
if($_POST["43"] != "0.01")
{
$data = "Charmonium:OJpsi1S08 = ".$_POST["43"]."\n";
fwrite($handle,$data);
}
if($_POST["44"] != "0.01")
{
$data = "Charmonium:OJpsi3P08 = ".$_POST["44"]."\n";
fwrite($handle,$data);
}
if($_POST["45"] != "0.05")
{
$data = "Charmonium:Ochic03P01 = ".$_POST["45"]."\n";
fwrite($handle,$data);
}
if($_POST["46"] != "9.28")
{
$data = "Bottomonium:OUpsilon3S11 = ".$_POST["46"]."\n";
fwrite($handle,$data);
}
if($_POST["47"] != "0.15")
{
$data = "Bottomonium:OUpsilon3S18 = ".$_POST["47"]."\n";
fwrite($handle,$data);
}
if($_POST["48"] != "0.02")
{
$data = "Bottomonium:OUpsilon1S08 = ".$_POST["48"]."\n";
fwrite($handle,$data);
}
if($_POST["49"] != "0.02")
{
$data = "Bottomonium:OUpsilon3P08 = ".$_POST["49"]."\n";
fwrite($handle,$data);
}
if($_POST["50"] != "0.085")
{
$data = "Bottomonium:Ochib03P01 = ".$_POST["50"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->

