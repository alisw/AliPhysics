<html>
<head>
<title>Electroweak Processes</title>
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

<form method='post' action='ElectroweakProcesses.php'>

<h2>Electroweak Processes</h2>

This page contains processes involving Prompt-photon, <i>gamma^*/Z^0</i> 
and <i>W^+-</i> production, plus a few with <i>t</i>-channel boson 
exchange. 

<h3>Prompt photon processes</h3>

This group collects the processes where one or two photons are
produced by the hard process. Additional sources of photons 
include parton showers and hadron decays. A <i>pT</i> cut
is required to stay away from the unphysical low-<i>pT</i> region.
An eikonalized description, intended to be valid at all <i>pT</i>,
is included as part of the multiparton-interactions framework.

<br/><br/><strong>PromptPhoton:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of all prompt photon processes, 
as listed separately in the following.
  

<br/><br/><strong>PromptPhoton:qg2qgamma</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> q gamma</i>.
Code 201.
  

<br/><br/><strong>PromptPhoton:qqbar2ggamma</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> g gamma</i>.
Code 202.
  

<br/><br/><strong>PromptPhoton:gg2ggamma</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> g gamma</i>.
<br/><b>Note:</b> This is a box graph. The full quark-mass dependence 
in the loop leads to very complicated expressions. The current 
implementation is based on assuming five massless quarks (see below), 
and thus is questionable at small (<i>pT &lt; m_b</i>) or large 
(<i>pT > m_t</i>) transverse momenta.
Code 203.
  

<br/><br/><strong>PromptPhoton:ffbar2gammagamma</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> gamma gamma</i>.
Code 204.
  

<br/><br/><strong>PromptPhoton:gg2gammagamma</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> gamma gamma</i>.
<br/><b>Note:</b> This is a box graph. The full quark-mass dependence 
in the loop leads to very complicated expressions. The current 
implementation is based on assuming five massless quarks (see below), 
and thus is questionable at small (<i>pT &lt; m_b</i>) or large 
(<i>pT > m_t</i>) transverse momenta.
Code 205.
  

<br/><br/><table><tr><td><strong>PromptPhoton:nQuarkLoop  </td><td></td><td> <input type="text" name="7" value="5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 3</code>; <code>maximum = 6</code>)</td></tr></table>
Number of quark flavours included in the box graphs responsible for 
<i>g g -> g gamma</i> and <i>g g-> gamma gamma</i> processes.
Owing to the complexity if the massive expressions, quarks are treated 
as massless. The default value should be applicable in the range of 
transverse momenta above the <i>b</i> mass but below the <i>t</i> one.
  

<h3>Weak boson processes</h3>

Under this heading we group processes involving the production
of a single electroweak gauge boson, i.e. a <i>gamma^*/Z^0</i>
or a <i>W^+-</i>, or a pair of them, or one of them in 
combination with a parton. Since the three sets are partly 
conflicting, each is associated with its own group flag.
In addition, <i>t</i>-channel exchange of such a boson 
between two fermions form a separate group.

<p/>
There is one flag that can be used to influence the <i>gamma^*/Z^0</i>
structure in all the processes below where it is produced, unless 
otherwise stated. 
<br/><br/><table><tr><td><strong>WeakZ0:gmZmode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Choice of full <ei>gamma^*/Z^0</ei> structure or not in relevant 
processes.
<br/>
<input type="radio" name="8" value="0" checked="checked"><strong>0 </strong>: full <ei>gamma^*/Z^0</ei> structure, with interference included.<br/>
<input type="radio" name="8" value="1"><strong>1 </strong>: only pure <ei>gamma^*</ei> contribution.<br/>
<input type="radio" name="8" value="2"><strong>2 </strong>: only pure <ei>Z^0</ei> contribution.<br/>
<br/><b>Note</b>: irrespective of the option used, the particle produced 
will always be assigned code 23 for <ei>Z^0</ei>, and open decay channels
is purely dictated by what is set for the <ei>Z^0</ei>. 

<h4>Boson exchange</h4>

The two processes in this subgroup is included as part of the 
multiparton-interactions framework.

<br/><br/><strong>WeakBosonExchange:all</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of <i>gamma^*/Z^0</i>
or <i>W^+-</i> exchange between two fermions.
  

<br/><br/><strong>WeakBosonExchange:ff2ff(t:gmZ)</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f f' -> f f'</i> via </i>gamma^*/Z^0</i>
<i>t</i>-channel exchange, with full interference
between the <i>gamma^*</i> and <i>Z^0</i>.
Code 211.
  

<br/><br/><strong>WeakBosonExchange:ff2ff(t:W)</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f_1 f_2 -> f_3 f_4</i> via </i>W^+-</i>
<i>t</i>-channel exchange.
Code 212.
  

<h4>Single boson</h4>

<br/><br/><strong>WeakSingleBoson:all</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of a single <i>gamma^*/Z^0</i>
or <i>W^+-</i> production.
  

<br/><br/><strong>WeakSingleBoson:ffbar2gmZ</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> gamma^*/Z^0</i>, with full interference
between the <i>gamma^*</i> and <i>Z^0</i>.
Code 221.
  

<br/><br/><strong>WeakSingleBoson:ffbar2W</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar' -> W^+-</i>.
Code 222.
  

<br/><br/><strong>WeakSingleBoson:ffbar2ffbar(s:gm)</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> gamma^* -> f' fbar'</i>. Subset of 
process 221, but written as a <i>2 -> 2</i> process, so that 
<i>pT</i> can be used as ordering variable, e.g. in multiparton 
interactions. Hardcoded for the final state being either of the 
five quark flavours or three lepton ones. Not included in the 
<code>WeakSingleBoson:all</code> set, but included in the 
multiparton-interactions framework. 
Code 223.
  

<h4>Boson pair</h4>

Note that, in the decay of the two vector bosons produced by an 
<i>f fbar -> V V</i> process, the full four-fermion correlations 
from the leading-order <i>f fbar -> V V -> 4f</i> matrix elements 
are included [<a href="Bibliography.php" target="page">Gun86</a>] (with some extensions by the authors).
The matrix elements are provided in the double-resonant approach, i.e. 
excludes graph like <i>f fbar -> V -> f fbar -> f fbar V -> 4f</i>.

<br/><br/><strong>WeakDoubleBoson:all</strong>  <input type="radio" name="16" value="on"><strong>On</strong>
<input type="radio" name="16" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of pair production of <i>gamma^*/Z^0</i>
and <i>W^+-</i>.
  
 
<br/><br/><strong>WeakDoubleBoson:ffbar2gmZgmZ</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar' -> gamma^*/Z^0 gamma^*/Z^0</i>. 
Code 231.
  

<br/><br/><strong>WeakDoubleBoson:ffbar2ZW</strong>  <input type="radio" name="18" value="on"><strong>On</strong>
<input type="radio" name="18" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar' -> Z^0 W^+-</i>. Note that here the 
<i>gamma^*</i> contribution is not (currently) included.
Code 232.
  
 
<br/><br/><strong>WeakDoubleBoson:ffbar2WW</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> W^+ W^-</i>.
Code 233.
  

<h4>Boson and parton</h4>

<br/><br/><strong>WeakBosonAndParton:all</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of production of a single electroweak 
gauge boson, i.e. a <i>gamma^*/Z^0</i> or a <i>W^+-</i>, in 
association with a parton, i.e. a quark, gluon, photon or lepton.
These processes give first-order corrections to the ones in the
<code>WeakSingleBoson</code> class, and both sets cannot be used
simultaneously without unphysical double-counting. The current class
should only be used to study the high-<i>pT</i> tail of the 
gauge-boson production processes (for LHC applications at least
<i>pT</i> > 20 GeV), while the ones in <code>WeakSingleBoson</code> 
should be used for inclusive production.
  
 
<br/><br/><strong>WeakBosonAndParton:qqbar2gmZg</strong>  <input type="radio" name="21" value="on"><strong>On</strong>
<input type="radio" name="21" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> gamma^*/Z^0 g</i>.
Code 241.
  
 
<br/><br/><strong>WeakBosonAndParton:qg2gmZq</strong>  <input type="radio" name="22" value="on"><strong>On</strong>
<input type="radio" name="22" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> gamma^*/Z^0 q </i>.
Code 242.
  
 
<br/><br/><strong>WeakBosonAndParton:ffbar2gmZgm</strong>  <input type="radio" name="23" value="on"><strong>On</strong>
<input type="radio" name="23" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> gamma^*/Z^0 gamma</i>.
Code 243.
  
 
<br/><br/><strong>WeakBosonAndParton:fgm2gmZf</strong>  <input type="radio" name="24" value="on"><strong>On</strong>
<input type="radio" name="24" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f gamma ->  gamma^*/Z^0 f</i>.
Code 244.
  
 
<br/><br/><strong>WeakBosonAndParton:qqbar2Wg</strong>  <input type="radio" name="25" value="on"><strong>On</strong>
<input type="radio" name="25" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> W^+- g</i>.
Code 251.
  
 
<br/><br/><strong>WeakBosonAndParton:qg2Wq</strong>  <input type="radio" name="26" value="on"><strong>On</strong>
<input type="radio" name="26" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> W^+- q</i>.
Code 252.
  
 
<br/><br/><strong>WeakBosonAndParton:ffbar2Wgm</strong>  <input type="radio" name="27" value="on"><strong>On</strong>
<input type="radio" name="27" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> W^+- gamma</i>.
Code 253.
  
 
<br/><br/><strong>WeakBosonAndParton:fgm2Wf</strong>  <input type="radio" name="28" value="on"><strong>On</strong>
<input type="radio" name="28" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f gamma -> W^+- f</i>.
Code 254.
  

<h3> Photon Collision Processes</h3>

A few electroweak two-photon production processes are available.
To use them, photon PDFs have to be defined for the incoming
beam particles. For proton beams an appropriate set would be
MRST QED 2004 [<a href="Bibliography.php" target="page">Mar05</a>], available in the LHAPDF library.

<br/><br/><strong>PhotonCollision:all</strong>  <input type="radio" name="29" value="on"><strong>On</strong>
<input type="radio" name="29" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of six processes presented below.
  

<br/><br/><strong>PhotonCollision:gmgm2qqbar</strong>  <input type="radio" name="30" value="on"><strong>On</strong>
<input type="radio" name="30" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>gamma gamma -> q qbar</i>, where <i>q</i> 
is a light quark (<i>u, d, s</i>) .
Code 261.
  

<br/><br/><strong>PhotonCollision:gmgm2ccbar</strong>  <input type="radio" name="31" value="on"><strong>On</strong>
<input type="radio" name="31" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>gamma gamma -> c cbar</i>.
Code 262.
  

<br/><br/><strong>PhotonCollision:gmgm2bbbar</strong>  <input type="radio" name="32" value="on"><strong>On</strong>
<input type="radio" name="32" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>gamma gamma -> b bbar</i>.
Code 263.
  

<br/><br/><strong>PhotonCollision:gmgm2ee</strong>  <input type="radio" name="33" value="on"><strong>On</strong>
<input type="radio" name="33" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>gamma gamma -> e+ e-</i>.
Code 264.
  

<br/><br/><strong>PhotonCollision:gmgm2mumu</strong>  <input type="radio" name="34" value="on"><strong>On</strong>
<input type="radio" name="34" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>gamma gamma -> mu+ mu-</i>.
Code 265.
  

<br/><br/><strong>PhotonCollision:gmgm2tautau</strong>  <input type="radio" name="35" value="on"><strong>On</strong>
<input type="radio" name="35" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>gamma gamma -> tau+ tau-</i>.
Code 266.
  

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
$data = "PromptPhoton:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "PromptPhoton:qg2qgamma = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "PromptPhoton:qqbar2ggamma = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "PromptPhoton:gg2ggamma = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "PromptPhoton:ffbar2gammagamma = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "PromptPhoton:gg2gammagamma = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "5")
{
$data = "PromptPhoton:nQuarkLoop = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "0")
{
$data = "WeakZ0:gmZmode = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "WeakBosonExchange:all = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "WeakBosonExchange:ff2ff(t:gmZ) = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "WeakBosonExchange:ff2ff(t:W) = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "WeakSingleBoson:all = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "WeakSingleBoson:ffbar2gmZ = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "WeakSingleBoson:ffbar2W = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "WeakSingleBoson:ffbar2ffbar(s:gm) = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off")
{
$data = "WeakDoubleBoson:all = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "WeakDoubleBoson:ffbar2gmZgmZ = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "off")
{
$data = "WeakDoubleBoson:ffbar2ZW = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "WeakDoubleBoson:ffbar2WW = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "WeakBosonAndParton:all = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "off")
{
$data = "WeakBosonAndParton:qqbar2gmZg = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "off")
{
$data = "WeakBosonAndParton:qg2gmZq = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "off")
{
$data = "WeakBosonAndParton:ffbar2gmZgm = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "off")
{
$data = "WeakBosonAndParton:fgm2gmZf = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "off")
{
$data = "WeakBosonAndParton:qqbar2Wg = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "off")
{
$data = "WeakBosonAndParton:qg2Wq = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "off")
{
$data = "WeakBosonAndParton:ffbar2Wgm = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "off")
{
$data = "WeakBosonAndParton:fgm2Wf = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "off")
{
$data = "PhotonCollision:all = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "off")
{
$data = "PhotonCollision:gmgm2qqbar = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "off")
{
$data = "PhotonCollision:gmgm2ccbar = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "off")
{
$data = "PhotonCollision:gmgm2bbbar = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "off")
{
$data = "PhotonCollision:gmgm2ee = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "off")
{
$data = "PhotonCollision:gmgm2mumu = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "off")
{
$data = "PhotonCollision:gmgm2tautau = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->

