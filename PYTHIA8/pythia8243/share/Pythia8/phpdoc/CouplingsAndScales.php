<html>
<head>
<title>Couplings and Scales</title>
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

<form method='post' action='CouplingsAndScales.php'>
 
<h2>Couplings and Scales</h2> 
<ol id="toc">
  <li><a href="#section0">Couplings and K factor</a></li>
  <li><a href="#section1">Renormalization scales</a></li>
  <li><a href="#section2">Factorization scales</a></li>
</ol>

 
Here is collected some possibilities to modify the scale choices 
of couplings and parton densities for all internally implemented 
hard processes. This is based on them all being derived from the 
<code>SigmaProcess</code> base class. The matrix-element coding is 
also used by the multiparton-interactions machinery, but there with a 
separate choice of <i>alpha_strong(M_Z^2)</i> value and running, 
and separate PDF scale choices. Also, in <i>2 &rarr; 2</i> and 
<i>2 &rarr; 3</i> processes where resonances are produced, their 
couplings and thereby their Breit-Wigner shapes are always evaluated 
with the resonance mass as scale, irrespective of the choices below. 
 
<p/> 
We stress that couplings and scales are set separately from the 
values on this page for 
<?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?>multiparton interactions</a>, 
<?php $filepath = $_GET["filepath"];
echo "<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>";?>timelike showers</a>, and 
<?php $filepath = $_GET["filepath"];
echo "<a href='SpacelikeShowers.php?filepath=".$filepath."' target='page'>";?>spacelike showers</a>. 
This allows a bigger flexibility, but also requires a bit more work 
e.g. if you insist on using the same <i>alpha_s</i> everywhere. 
 
<a name="section0"></a> 
<h3>Couplings and K factor</h3> 
 
The size of QCD cross sections is mainly determined by 
<br/><br/><table><tr><td><strong>SigmaProcess:alphaSvalue </td><td></td><td> <input type="text" name="1" value="0.13" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.13</strong></code>; <code>minimum = 0.06</code>; <code>maximum = 0.25</code>)</td></tr></table>
The <i>alpha_strong</i> value at scale <i>M_Z^2</i>. 
   
 
<p/> 
The actual value is then regulated by the running to the <i>Q^2</i> 
renormalization scale, at which <i>alpha_strong</i> is evaluated 
<br/><br/><table><tr><td><strong>SigmaProcess:alphaSorder  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Order at which <ei>alpha_strong</ei> runs, 
<br/>
<input type="radio" name="2" value="0"><strong>0 </strong>: zeroth order, i.e. <ei>alpha_strong</ei> is kept  fixed.<br/>
<input type="radio" name="2" value="1" checked="checked"><strong>1 </strong>: first order, which is the normal value.<br/>
<input type="radio" name="2" value="2"><strong>2 </strong>: second order. Since other parts of the code do  not go to second order there is no strong reason to use this option,  but there is also nothing wrong with it.<br/>
 
<p/> 
QED interactions are regulated by the <i>alpha_electromagnetic</i> 
value at the <i>Q^2</i> renormalization scale of an interaction. 
<br/><br/><table><tr><td><strong>SigmaProcess:alphaEMorder  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = -1</code>; <code>maximum = 1</code>)</td></tr></table>
The running of <ei>alpha_em</ei> used in hard processes. 
<br/>
<input type="radio" name="3" value="1" checked="checked"><strong>1 </strong>: first-order running, constrained to agree with  <code>StandardModel:alphaEMmZ</code> at the <ei>Z^0</ei> mass.  <br/>
<input type="radio" name="3" value="0"><strong>0 </strong>: zeroth order, i.e. <ei>alpha_em</ei> is kept  fixed at its value at vanishing momentum transfer.<br/>
<input type="radio" name="3" value="-1"><strong>-1 </strong>: zeroth order, i.e. <ei>alpha_em</ei> is kept  fixed, but at <code>StandardModel:alphaEMmZ</code>, i.e. its value  at the <ei>Z^0</ei> mass.  <br/>
 
<p/> 
In addition there is the possibility of a global rescaling of 
cross sections (which could not easily be accommodated by a 
changed <i>alpha_strong</i>, since <i>alpha_strong</i> runs) 
<br/><br/><table><tr><td><strong>SigmaProcess:Kfactor </td><td></td><td> <input type="text" name="4" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.5</code>; <code>maximum = 4.0</code>)</td></tr></table>
Multiply almost all cross sections by this common fix factor. Excluded 
are only unresolved processes, where cross sections are better 
<?php $filepath = $_GET["filepath"];
echo "<a href='TotalCrossSections.php?filepath=".$filepath."' target='page'>";?>set directly</a>, and 
multiparton interactions, which have a separate <i>K</i> factor 
<?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?>of their own</a>. 
This degree of freedom is primarily intended for hadron colliders, and 
should not normally be used for <i>e^+e^-</i> annihilation processes. 
   
 
<a name="section1"></a> 
<h3>Renormalization scales</h3> 
 
The <i>Q^2</i> renormalization scale can be chosen among a few different 
alternatives, separately for <i>2 &rarr; 1</i>, <i>2 &rarr; 2</i> and two 
different kinds of <i>2 &rarr; 3</i> processes. In addition a common 
multiplicative factor may be imposed. 
 
<br/><br/><table><tr><td><strong>SigmaProcess:renormScale1  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
The <ei>Q^2</ei> renormalization scale for <ei>2 &rarr; 1</ei> processes. 
The same options also apply for those <ei>2 &rarr; 2</ei> and 
<ei>2 &rarr; 3</ei> processes that have been specially marked as 
proceeding only through an <ei>s</ei>-channel resonance, by the 
<code>isSChannel()</code> virtual method of <code>SigmaProcess</code>. 
<br/>
<input type="radio" name="5" value="1" checked="checked"><strong>1 </strong>: the squared invariant mass, i.e. <ei>sHat</ei>.  <br/>
<input type="radio" name="5" value="2"><strong>2 </strong>: fix scale set in <code>SigmaProcess:renormFixScale</code>  below.  <br/>
 
<br/><br/><table><tr><td><strong>SigmaProcess:renormScale2  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>; <code>maximum = 6</code>)</td></tr></table>
The <ei>Q^2</ei> renormalization scale for <ei>2 &rarr; 2</ei> processes. 
<br/>
<input type="radio" name="6" value="1"><strong>1 </strong>: the smaller of the squared transverse masses of the two  outgoing particles, i.e. <ei>min(mT_3^2, mT_4^2) =  pT^2 + min(m_3^2, m_4^2)</ei>.  <br/>
<input type="radio" name="6" value="2" checked="checked"><strong>2 </strong>: the geometric mean of the squared transverse masses of  the two outgoing particles, i.e. <ei>mT_3 * mT_4 =  sqrt((pT^2 + m_3^2) * (pT^2 + m_4^2))</ei>.  <br/>
<input type="radio" name="6" value="3"><strong>3 </strong>: the arithmetic mean of the squared transverse masses of  the two outgoing particles, i.e. <ei>(mT_3^2 + mT_4^2) / 2 =  pT^2 + 0.5 * (m_3^2 + m_4^2)</ei>. Useful for comparisons  with PYTHIA 6, where this is the default.  <br/>
<input type="radio" name="6" value="4"><strong>4 </strong>: squared invariant mass of the system,  i.e. <ei>sHat</ei>. Useful for processes dominated by  <ei>s</ei>-channel exchange.  <br/>
<input type="radio" name="6" value="5"><strong>5 </strong>: fix scale set in <code>SigmaProcess:renormFixScale</code>  below.  <br/>
<input type="radio" name="6" value="6"><strong>6 </strong>: Use squared invariant momentum transfer <ei>-tHat</ei>.  This is a common choice for lepton-hadron scattering processes. In that  case <ei>-tHat=Q^2</ei>.  <br/>
 
<br/><br/><table><tr><td><strong>SigmaProcess:renormScale3  </td><td>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 1</code>; <code>maximum = 6</code>)</td></tr></table>
The <ei>Q^2</ei> renormalization scale for "normal" <ei>2 &rarr; 3</ei> 
processes, i.e excepting the vector-boson-fusion processes below. 
Here it is assumed that particle masses in the final state either match 
or are heavier than that of any <ei>t</ei>-channel propagator particle. 
(Currently only <ei>g g / q qbar &rarr; H^0 Q Qbar</ei> processes are 
implemented, where the "match" criterion holds.) 
<br/>
<input type="radio" name="7" value="1"><strong>1 </strong>: the smaller of the squared transverse masses of the three  outgoing particles, i.e. min(mT_3^2, mT_4^2, mT_5^2).  <br/>
<input type="radio" name="7" value="2"><strong>2 </strong>: the geometric mean of the two smallest squared transverse  masses of the three outgoing particles, i.e.  <ei>sqrt( mT_3^2 * mT_4^2 * mT_5^2 / max(mT_3^2, mT_4^2, mT_5^2) )</ei>.  <br/>
<input type="radio" name="7" value="3" checked="checked"><strong>3 </strong>: the geometric mean of the squared transverse masses of the  three outgoing particles, i.e. <ei>(mT_3^2 * mT_4^2 * mT_5^2)^(1/3)</ei>.  <br/>
<input type="radio" name="7" value="4"><strong>4 </strong>: the arithmetic mean of the squared transverse masses of  the three outgoing particles, i.e. <ei>(mT_3^2 + mT_4^2 + mT_5^2)/3</ei>.  <br/>
<input type="radio" name="7" value="5"><strong>5 </strong>: squared invariant mass of the system,  i.e. <ei>sHat</ei>.  <br/>
<input type="radio" name="7" value="6"><strong>6 </strong>: fix scale set in <code>SigmaProcess:renormFixScale</code>  below.  <br/>
 
<br/><br/><table><tr><td><strong>SigmaProcess:renormScale3VV  </td><td>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 1</code>; <code>maximum = 6</code>)</td></tr></table>
The <ei>Q^2</ei> renormalization scale for <ei>2 &rarr; 3</ei> 
vector-boson-fusion processes, i.e. <ei>f_1 f_2 &rarr; H^0 f_3 f_4</ei> 
with <ei>Z^0</ei> or <ei>W^+-</ei>  <ei>t</ei>-channel propagators. 
Here the transverse masses of the outgoing fermions do not reflect the 
virtualities of the exchanged bosons. A better estimate is obtained 
by replacing the final-state fermion masses by the vector-boson ones 
in the definition of transverse masses. We denote these combinations 
<ei>mT_Vi^2 = m_V^2 + pT_i^2</ei>. 
<br/>
<input type="radio" name="8" value="1"><strong>1 </strong>: the squared mass <ei>m_V^2</ei> of the exchanged  vector boson.  <br/>
<input type="radio" name="8" value="2"><strong>2 </strong>: the geometric mean of the two propagator virtuality  estimates, i.e. <ei>sqrt(mT_V3^2 * mT_V4^2)</ei>.  <br/>
<input type="radio" name="8" value="3" checked="checked"><strong>3 </strong>: the geometric mean of the three relevant squared  transverse masses, i.e. <ei>(mT_V3^2 * mT_V4^2 * mT_H^2)^(1/3)</ei>.  <br/>
<input type="radio" name="8" value="4"><strong>4 </strong>: the arithmetic mean of the three relevant squared  transverse masses, i.e. <ei>(mT_V3^2 + mT_V4^2 + mT_H^2)/3</ei>.  <br/>
<input type="radio" name="8" value="5"><strong>5 </strong>: squared invariant mass of the system,  i.e. <ei>sHat</ei>.  <br/>
<input type="radio" name="8" value="6"><strong>6 </strong>: fix scale set in <code>SigmaProcess:renormFixScale</code>  below.  <br/>
 
<br/><br/><table><tr><td><strong>SigmaProcess:renormMultFac </td><td></td><td> <input type="text" name="9" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.</code>)</td></tr></table>
The <i>Q^2</i> renormalization scale for <i>2 &rarr; 1</i>, 
<i>2 &rarr; 2</i> and <i>2 &rarr; 3</i> processes is multiplied by 
this factor relative to the scale described above (except for the options 
with a fix scale). Should be use sparingly for <i>2 &rarr; 1</i> processes. 
   
 
<br/><br/><table><tr><td><strong>SigmaProcess:renormFixScale </td><td></td><td> <input type="text" name="10" value="10000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10000.</strong></code>; <code>minimum = 1.</code>)</td></tr></table>
A fix <i>Q^2</i> value used as renormalization scale for 
<i>2 &rarr; 1</i>, <i>2 &rarr; 2</i> and <i>2 &rarr; 3</i> processes 
in some of the options above. 
   
 
<a name="section2"></a> 
<h3>Factorization scales</h3> 
 
Corresponding options exist for the <i>Q^2</i> factorization scale 
used as argument in PDF's. Again there is a choice of form for 
<i>2 &rarr; 1</i>, <i>2 &rarr; 2</i> and <i>2 &rarr; 3</i> processes 
separately. For simplicity we have let the numbering of options agree, 
for each event class separately, between normalization and factorization 
scales, and the description has therefore been slightly shortened. The 
default values are <b>not</b> necessarily the same, however. 
 
<br/><br/><table><tr><td><strong>SigmaProcess:factorScale1  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
The <ei>Q^2</ei> factorization scale for <ei>2 &rarr; 1</ei> processes. 
The same options also apply for those <ei>2 &rarr; 2</ei> and 
<ei>2 &rarr; 3</ei> processes that have been specially marked as 
proceeding only through an <ei>s</ei>-channel resonance. 
<br/>
<input type="radio" name="11" value="1" checked="checked"><strong>1 </strong>: the squared invariant mass, i.e. <ei>sHat</ei>.  <br/>
<input type="radio" name="11" value="2"><strong>2 </strong>: fix scale set in <code>SigmaProcess:factorFixScale</code>  below.  <br/>
 
<br/><br/><table><tr><td><strong>SigmaProcess:factorScale2  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 6</code>)</td></tr></table>
The <ei>Q^2</ei> factorization scale for <ei>2 &rarr; 2</ei> processes. 
<br/>
<input type="radio" name="12" value="1" checked="checked"><strong>1 </strong>: the smaller of the squared transverse masses of the two  outgoing particles.  <br/>
<input type="radio" name="12" value="2"><strong>2 </strong>: the geometric mean of the squared transverse masses of  the two outgoing particles.  <br/>
<input type="radio" name="12" value="3"><strong>3 </strong>: the arithmetic mean of the squared transverse masses of  the two outgoing particles. Useful for comparisons with PYTHIA 6, where  this is the default.  <br/>
<input type="radio" name="12" value="4"><strong>4 </strong>: squared invariant mass of the system,  i.e. <ei>sHat</ei>. Useful for processes dominated by  <ei>s</ei>-channel exchange.  <br/>
<input type="radio" name="12" value="5"><strong>5 </strong>: fix scale set in <code>SigmaProcess:factorFixScale</code>  below.  <br/>
<input type="radio" name="12" value="6"><strong>6 </strong>: Use squared invariant momentum transfer <ei>-tHat</ei>.  This is a common choice for lepton-hadron scattering processes. In that  case <ei>-tHat=Q^2</ei>.  <br/>
 
<br/><br/><table><tr><td><strong>SigmaProcess:factorScale3  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>; <code>maximum = 6</code>)</td></tr></table>
The <ei>Q^2</ei> factorization scale for "normal" <ei>2 &rarr; 3</ei> 
processes, i.e excepting the vector-boson-fusion processes below. 
<br/>
<input type="radio" name="13" value="1"><strong>1 </strong>: the smaller of the squared transverse masses of the three  outgoing particles.  <br/>
<input type="radio" name="13" value="2" checked="checked"><strong>2 </strong>: the geometric mean of the two smallest squared transverse  masses of the three outgoing particles.  <br/>
<input type="radio" name="13" value="3"><strong>3 </strong>: the geometric mean of the squared transverse masses of the  three outgoing particles.  <br/>
<input type="radio" name="13" value="4"><strong>4 </strong>: the arithmetic mean of the squared transverse masses of  the three outgoing particles.  <br/>
<input type="radio" name="13" value="5"><strong>5 </strong>: squared invariant mass of the system,  i.e. <ei>sHat</ei>.  <br/>
<input type="radio" name="13" value="6"><strong>6 </strong>: fix scale set in <code>SigmaProcess:factorFixScale</code>  below.  <br/>
 
<br/><br/><table><tr><td><strong>SigmaProcess:factorScale3VV  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>; <code>maximum = 6</code>)</td></tr></table>
The <ei>Q^2</ei> factorization scale for <ei>2 &rarr; 3</ei> 
vector-boson-fusion processes, i.e. <ei>f_1 f_2 &rarr; H^0 f_3 f_4</ei> 
with <ei>Z^0</ei> or <ei>W^+-</ei>  <ei>t</ei>-channel propagators. 
Here we again introduce the combinations <ei>mT_Vi^2 = m_V^2 + pT_i^2</ei> 
as replacements for the normal squared transverse masses of the two 
outgoing quarks. 
<br/>
<input type="radio" name="14" value="1"><strong>1 </strong>: the squared mass <ei>m_V^2</ei> of the exchanged  vector boson.  <br/>
<input type="radio" name="14" value="2" checked="checked"><strong>2 </strong>: the geometric mean of the two propagator virtuality  estimates.  <br/>
<input type="radio" name="14" value="3"><strong>3 </strong>: the geometric mean of the three relevant squared  transverse masses.  <br/>
<input type="radio" name="14" value="4"><strong>4 </strong>: the arithmetic mean of the three relevant squared  transverse masses.  <br/>
<input type="radio" name="14" value="5"><strong>5 </strong>: squared invariant mass of the system,  i.e. <ei>sHat</ei>.  <br/>
<input type="radio" name="14" value="6"><strong>6 </strong>: fix scale set in <code>SigmaProcess:factorFixScale</code>  below.  <br/>
 
<br/><br/><table><tr><td><strong>SigmaProcess:factorMultFac </td><td></td><td> <input type="text" name="15" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.</code>)</td></tr></table>
The <i>Q^2</i> factorization scale for <i>2 &rarr; 1</i>, 
<i>2 &rarr; 2</i> and <i>2 &rarr; 3</i> processes is multiplied by 
this factor relative to the scale described above (except for the options 
with a fix scale). Should be use sparingly for <i>2 &rarr; 1</i> processes. 
   
 
<br/><br/><table><tr><td><strong>SigmaProcess:factorFixScale </td><td></td><td> <input type="text" name="16" value="10000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10000.</strong></code>; <code>minimum = 1.</code>)</td></tr></table>
A fix <i>Q^2</i> value used as factorization scale for <i>2 &rarr; 1</i>, 
<i>2 &rarr; 2</i> and <i>2 &rarr; 3</i> processes in some of the options 
above. 
   
 
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

if($_POST["1"] != "0.13")
{
$data = "SigmaProcess:alphaSvalue = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "1")
{
$data = "SigmaProcess:alphaSorder = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "1")
{
$data = "SigmaProcess:alphaEMorder = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "1.0")
{
$data = "SigmaProcess:Kfactor = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1")
{
$data = "SigmaProcess:renormScale1 = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "2")
{
$data = "SigmaProcess:renormScale2 = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "3")
{
$data = "SigmaProcess:renormScale3 = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "3")
{
$data = "SigmaProcess:renormScale3VV = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "1.")
{
$data = "SigmaProcess:renormMultFac = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "10000.")
{
$data = "SigmaProcess:renormFixScale = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "1")
{
$data = "SigmaProcess:factorScale1 = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "1")
{
$data = "SigmaProcess:factorScale2 = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "2")
{
$data = "SigmaProcess:factorScale3 = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "2")
{
$data = "SigmaProcess:factorScale3VV = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "1.")
{
$data = "SigmaProcess:factorMultFac = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "10000.")
{
$data = "SigmaProcess:factorFixScale = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
