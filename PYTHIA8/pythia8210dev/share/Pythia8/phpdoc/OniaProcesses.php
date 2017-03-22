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
 
Production of any <i>3S1</i>, <i>3PJ</i>, and <i>3DJ</i> charmonium 
and bottomonium states via the colour-singlet and colour-octet 
mechanisms. This includes by default, but is not limited to, production of 
the <i>3S1</i> <i>J/psi</i> and <i>Upsilon</i> and their 
radially excited states, as well as the <i>3PJ</i> <i>chi</i> 
states and the <i>3D1</i> <i>psi(3770)</i>. In each process the 
heavy quark content, either <i>ccbar</i> or <i>bbbar</i>, is 
followed by a round-bracketed expression which specifies the physical 
state in spectroscopic notation, <i>(2S+1) L J</i>. Proceding this 
is a square-bracketed expression, also in spectroscopic notation, 
which specifies the Fock state through which the process occurs, 
where <i>(1)</i> indicates a colour-singlet state and <i>(8)</i> a 
colour-octet state. 
 
<p> The unphysical colour-octet states follow the <code>id</code> 
scheme of <i>99 n_q n_s n_r n_L n_J</i> where <i>n_q</i> is the 
quark flavour of the state and <i>n_s</i> is the colour-octet state 
type. Here <i>0</i> is <i>3S1</i>, <i>1</i> is <i>1S0</i>, 
and <i>2</i> is <i>3PJ</i>. All remaining numbers follow the 
standard PDG numbering scheme. If a physical state is requested 
without a corresponding colour-octet state, a colour-octet state is 
automatically added to the <code>ParticleData</code> 
when a colour-octet process is selected. The colour-octet state is 
created with a mass given by the mass of the physical state plus the 
singlet-octet mass splitting parameter <code>Onia:massSplit</code>, 
which is by default set at 200 MeV, and decays exclusively 
to a gluon and the physical state. If the user wishes to manually 
set the mass splitting for each colour-octet state individually 
then <code>Onia:forceMassSplit</code> can be set to <i>off</i>. 
By default the widths of the octet states are set to vanish. 
This is not realistic, given their presumably rather rapid decay, 
but a nonvanishing width is not likely to have any measurable 
consequences that go beyond what comes from viewing the singlet-octet 
mass splitting as an effective parameter. 
 
<p/> 
The original Fortran code for these processes has been contributed 
by Stefan Wolf [unpublished]. For the C++ version only the unpolarized 
expressions are retained, since the theoretical predictions of the 
colour-octet model anyway do not agree with the experimental 
observations. Furthermore, the polarization effects are modest, 
so isotropic decay is not a bad starting point. Such an event sample 
can afterwards be reweighted at will by the user, to test various 
assumptions. The expressions for the colour-singlet production of 
the <i>3S1</i> and <i>3PJ</i> states can be found 
in [<a href="Bibliography.php" target="page">Bai83</a>] and [<a href="Bibliography.php" target="page">Gas87</a>]. Colour-octet expressions can 
be found in [<a href="Bibliography.php" target="page">Cho96</a>] for the <i>1S0</i>, <i>3S1</i>, 
and <i>3PJ</i> states, and the matrix elements for the <i>3DJ</i> 
states are taken from [<a href="Bibliography.php" target="page">Yua98</a>]. 
 
<p/> 
The implementation of charmonium and bottomonium production, including 
the colour-octet production mechanism, requires information on 
long-distance NRQCD matrix elements for the various wavefunctions 
involved. Default values for these are encoded in the <i>O</i> 
parameters and are taken from [<a href="Bibliography.php" target="page">Nas00</a>]; see 
also [<a href="Bibliography.php" target="page">Bar07</a>]. The <i>3DJ</i> long-distance matrix elements 
are extracted from [<a href="Bibliography.php" target="page">Yua98</a>]. 
 
<p/> 
Note that states that differ only by the radial excitation number 
<i>n_r</i> share the same short-dinstence matrix elements. The 
program has therefore been written such that further radial excitations 
can be easily added by editing this file, without requiring a recompilation 
of the code. All related arrays must be expanded in exactly the same way, 
however, i.e. the code of the colour singlet state, the long-distance 
matrix elements and the individual process on/off switches. 
 
<p/> 
The description of 
<?php $filepath = $_GET["filepath"];
echo "<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>";?>final-state radiation</a> 
is in this case based on some further model assumptions. 
 
<p/> 
Most of the processes below are divergent in the limit <i>pT &rarr; 0</i>, 
and therefore a <i>pTmin</i> scale should be set. Comparisons with 
data indicate that this divergence can be tamed the same way as for 
the normal QCD <i>2 &rarr; 2</i> cross sections [<a href="Bibliography.php" target="page">Bar07,Kra08</a>], 
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
 
<br/><br/><strong>Onia:forceMassSplit</strong>  <input type="radio" name="1" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="1" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Force the mass splitting between the colour-singlet states and their 
corresponding colour-octet state to be <code>Onia:massSplit</code>. 
   
 
<br/><br/><table><tr><td><strong>Onia:massSplit </td><td></td><td> <input type="text" name="2" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
Mass splitting in GeV between the physical colour-singlet 
states and their corresponding colour-octet state. 
   
 
<br/><br/><strong>Onia:all</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of onia production. 
   
<br/><br/><strong>Onia:all(3S1)</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of <i>3S1</i> onia production, 
e.g. <i>J/psi</i> and <i>Upsilon</i>. 
   
<br/><br/><strong>Onia:all(3PJ)</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of <i>3PJ</i> onia production, 
e.g. <i>chi_c</i> and <i>chi_b</i>. 
   
<br/><br/><strong>Onia:all(3DJ)</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of <i>3DJ</i> onia production, 
e.g. <i>psi(3770)</i>. 
   
<br/><br/><strong>Charmonium:all</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of charmonium production, 
e.g. <i>J/psi</i> and <i>chi_c</i>. 
   
<br/><br/><strong>Bottomonium:all</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of bottomonium production, 
e.g. <i>Upsilon</i> and <i>chi_b</i>. 
   
 
<h3>Charmonium 3S1 States</h3> 
 
<b>Warning</b>: changed <code>fvec</code>, <code>mvec</code> or 
<code>pvec</code> values must be provided as a comma-separated list 
with the right number of elements, without any blanks inside the list. 
 
<br/><br/><table><tr><td><strong>Charmonium:states(3S1)  </td><td></td><td> <input type="text" name="9" value="443,100443" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>443,100443</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The <i>3S1</i> charmonium states that can be produced from the following 
processes. Note that all vectors within this section, 
either of flags or parameters, must be the same length as this 
vector. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:O(3S1)[3S1(1)] </td><td></td><td> <input type="text" name="10" value="1.16,0.76" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.16,0.76</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The colour-singlet long-distance matrix 
elements <i>&lt;O[3S1(1)]&gt;</i> for the <i>3S1</i> charmonium states. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:O(3S1)[3S1(8)] </td><td></td><td> <input type="text" name="11" value="0.0119,0.0050" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0119,0.0050</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The colour-octet long-distance matrix 
elements <i>&lt;O[3S1(8)]&gt;</i> for the <i>3S1</i> charmonium states. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:O(3S1)[1S0(8)] </td><td></td><td> <input type="text" name="12" value="0.01,0.004" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01,0.004</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The colour-octet long-distance matrix 
elements <i>&lt;O[1S0(8)]&gt;</i> for the <i>3S1</i> 
charmonium states. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:O(3S1)[3P0(8)] </td><td></td><td> <input type="text" name="13" value="0.01,0.004" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01,0.004</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The colour-octet long-distance matrix 
elements <i>&lt;O[3P0(8)]&gt;/m_Q^2</i> for the <i>3S1</i> charmonium 
states. The remaining <i>&lt;O[3PJ(8)]&gt;/m_Q^2</i> 
are calculated from these long-distance matrix elements. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:gg2ccbar(3S1)[3S1(1)]g  </td><td></td><td> <input type="text" name="14" value="off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off</strong></code>)</td></tr></table>
Colour-singlet production of <i>3S1</i> charmonium states via 
<i>g g &rarr; ccbar[3S1(1)] g</i>. 
Code 401. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:gg2ccbar(3S1)[3S1(8)]g  </td><td></td><td> <input type="text" name="15" value="off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> charmonium states via 
<i>g g &rarr; ccbar[3S1(8)] g</i>. 
Code 402. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qg2ccbar(3S1)[3S1(8)]q  </td><td></td><td> <input type="text" name="16" value="off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> charmonium states via 
<i>q g &rarr; ccbar[3S1(8)] q</i>. 
Code 403. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qqbar2ccbar(3S1)[3S1(8)]g  </td><td></td><td> <input type="text" name="17" value="off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> charmonium states via 
<i>q qbar &rarr; ccbar[3S1(8)] g</i>. 
Code 404. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:gg2ccbar(3S1)[1S0(8)]g  </td><td></td><td> <input type="text" name="18" value="off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> charmonium states via 
<i>g g &rarr; ccbar[1S0(8)] g</i>. 
Code 405. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qg2ccbar(3S1)[1S0(8)]q  </td><td></td><td> <input type="text" name="19" value="off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> charmonium states via 
<i>q g &rarr; ccbar[1S0(8)] q</i>. 
Code 406. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qqbar2ccbar(3S1)[1S0(8)]g  </td><td></td><td> <input type="text" name="20" value="off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> charmonium states via 
<i>q qbar &rarr; ccbar[1S0(8)] g</i>. 
Code 407. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:gg2ccbar(3S1)[3PJ(8)]g  </td><td></td><td> <input type="text" name="21" value="off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> charmonium states via 
<i>g g &rarr; ccbar[3PJ(8)] g</i>. 
Code 408. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qg2ccbar(3S1)[3PJ(8)]q  </td><td></td><td> <input type="text" name="22" value="off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> charmonium states via 
<i>q g &rarr; ccbar[3PJ(8)] q</i>. 
Code 409. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qqbar2ccbar(3S1)[3PJ(8)]g  </td><td></td><td> <input type="text" name="23" value="off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> charmonium states via 
<i>q qbar &rarr; ccbar[3SJ(8)] g</i>. 
Code 410. 
   
 
<h3>Charmonium 3PJ States</h3> 
 
<b>Warning</b>: changed <code>fvec</code>, <code>mvec</code> or 
<code>pvec</code> values must be provided as a comma-separated list 
with the right number of elements, without any blanks inside the list. 
 
<br/><br/><table><tr><td><strong>Charmonium:states(3PJ)  </td><td></td><td> <input type="text" name="24" value="10441,20443,445" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10441,20443,445</strong></code>)</td></tr></table>
The <i>3PJ</i> charmonium states that can be produced from the following 
processes. Note that all vectors within this section, 
either of flags or parameters, must be the same length as this 
vector. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:O(3PJ)[3P0(1)] </td><td></td><td> <input type="text" name="25" value="0.05,0.05,0.05" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.05,0.05,0.05</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The color-singlet long-distance matrix elements 
<i>&lt;O[3P0(1)]&gt;/m_Q^2</i> for the <i>3PJ</i> charmonium 
states. The remaining <i>&lt;O[3PJ(1)]&gt;/m_Q^2</i> 
are calculated from these long-distance matrix elements. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:O(3PJ)[3S1(8)] </td><td></td><td> <input type="text" name="26" value="0.0031,0.0031,0.0031" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0031,0.0031,0.0031</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The color-singlet long-distance matrix elements 
<i>&lt;O[3S1(8)]&gt;</i> for the <i>3PJ</i> charmonium states. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:gg2ccbar(3PJ)[3PJ(1)]g  </td><td></td><td> <input type="text" name="27" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-singlet production of <i>3PJ</i> charmonium states via 
<i>g g &rarr; ccbar[3PJ(1)] g</i>. 
Code 411. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qg2ccbar(3PJ)[3PJ(1)]q  </td><td></td><td> <input type="text" name="28" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-singlet production of <i>3PJ</i> charmonium states via 
<i>q g &rarr; ccbar[3PJ(1)] q</i>. 
Code 412. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qqbar2ccbar(3PJ)[3PJ(1)]g  </td><td></td><td> <input type="text" name="29" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-singlet production of <i>3PJ</i> charmonium states via 
<i>q qbar &rarr; ccbar[3PJ(1)] g</i>. 
Code 413. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:gg2ccbar(3PJ)[3S1(8)]g  </td><td></td><td> <input type="text" name="30" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3PJ</i> charmonium states via 
<i>g g &rarr; ccbar[3S1(8)] g</i>. 
Code 414. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qg2ccbar(3PJ)[3S1(8)]q  </td><td></td><td> <input type="text" name="31" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3PJ</i> charmonium states via 
<i>q g &rarr; ccbar[3S1(8)] q</i>. 
Code 415. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qqbar2ccbar(3PJ)[3S1(8)]g  </td><td></td><td> <input type="text" name="32" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3PJ</i> charmonium states via 
<i>q qbar &rarr; ccbar[3S1(8)] g</i>. 
Code 416. 
   
 
<h3>Charmonium 3DJ States</h3> 
 
<b>Warning</b>: changed <code>fvec</code>, <code>mvec</code> or 
<code>pvec</code> values must be provided as a comma-separated list 
with the right number of elements, without any blanks inside the list. 
 
<br/><br/><table><tr><td><strong>Charmonium:states(3DJ)  </td><td></td><td> <input type="text" name="33" value="30443" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>30443</strong></code>)</td></tr></table>
The <i>3DJ</i> charmonium states that can be produced from the following 
processes. Note that all vectors within this section, 
either of flags or parameters, must be the same length as this 
vector. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:O(3DJ)[3D1(1)] </td><td></td><td> <input type="text" name="34" value="0.161" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.161</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The color-singlet long-distance matrix elements 
<i>&lt;O[3D1(1)]&gt;</i> for the <i>3PJ</i> charmonium 
states. For a <i>3DJ</i> charmonium state where <i>J</i> is 
not <i>1</i> the long distance matrix 
element <i>&lt;O[3DJ(1)]&gt;</i> is calculated 
by <i>(2J+1)&lt;O[3D1(1)]/3&gt;</i> using leading order spin symmetry 
relations. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:O(3DJ)[3P0(8)] </td><td></td><td> <input type="text" name="35" value="0.01" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The colour-octet long-distance matrix 
elements <i>&lt;O[3P0(8)]&gt;/m_Q^2</i> for the 3DJ charmonium 
states. The remaining <i>&lt;O[3PJ(8)]&gt;/m_Q^2</i> 
are calculated from these long-distance matrix elements. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:gg2ccbar(3DJ)[3DJ(1)]g  </td><td></td><td> <input type="text" name="36" value="off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)</td></tr></table>
Colour-singlet production of <i>3PJ</i> charmonium states via 
<i>g g &rarr; ccbar[3DJ(1)] g</i>. 
Code 417. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:gg2ccbar(3DJ)[3PJ(8)]g  </td><td></td><td> <input type="text" name="37" value="off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)</td></tr></table>
Colour-octet production of <i>3DJ</i> charmonium states via 
<i>g g &rarr; ccbar[3PJ(8)] g</i>. 
Code 418. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qg2ccbar(3DJ)[3PJ(8)]q  </td><td></td><td> <input type="text" name="38" value="off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)</td></tr></table>
Colour-octet production of <i>3DJ</i> charmonium states via 
<i>q g &rarr; ccbar[3PJ(8)] q</i>. 
Code 419. 
   
 
<br/><br/><table><tr><td><strong>Charmonium:qqbar2ccbar(3DJ)[3PJ(8)]g  </td><td></td><td> <input type="text" name="39" value="off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)</td></tr></table>
Colour-octet production of <i>3DJ</i> charmonium states via 
<i>q qbar &rarr; ccbar[3PJ(8)] g</i>. 
Code 420. 
   
 
<h3>Bottomonium 3S1 States</h3> 
 
<b>Warning</b>: changed <code>fvec</code>, <code>mvec</code> or 
<code>pvec</code> values must be provided as a comma-separated list 
with the right number of elements, without any blanks inside the list. 
 
<br/><br/><table><tr><td><strong>Bottomonium:states(3S1)  </td><td></td><td> <input type="text" name="40" value="553,100553,200553" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>553,100553,200553</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The <i>3S1</i> bottomonium states that can be produced from the following 
processes. Note that all vectors within this section, 
either of flags or parameters, must be the same length as this 
vector. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:O(3S1)[3S1(1)] </td><td></td><td> <input type="text" name="41" value="9.28,4.63,3.54" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>9.28,4.63,3.54</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The colour-singlet long-distance matrix 
elements <i>&lt;O[3S1(1)]&gt;</i> for the <i>3S1</i> bottomonium states. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:O(3S1)[3S1(8)] </td><td></td><td> <input type="text" name="42" value="0.15,0.045,0.075" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.15,0.045,0.075</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The colour-octet long-distance matrix 
elements <i>&lt;O[3S1(8)]&gt;</i> for the <i>3S1</i> bottomonium states. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:O(3S1)[1S0(8)] </td><td></td><td> <input type="text" name="43" value="0.02,0.06,0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.02,0.06,0.1</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The colour-octet long-distance matrix 
elements <i>&lt;O[1S0(8)]&gt;</i> for the <i>3S1</i> 
bottomonium states. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:O(3S1)[3P0(8)] </td><td></td><td> <input type="text" name="44" value="0.02,0.06,0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.02,0.06,0.1</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The colour-octet long-distance matrix 
elements <i>&lt;O[3P0(8)]&gt;/m_Q^2</i> for the <i>3S1</i> bottomonium 
states. The remaining <i>&lt;O[3PJ(8)]&gt;/m_Q^2</i> 
are calculated from these long-distance matrix elements. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:gg2bbbar(3S1)[3S1(1)]g  </td><td></td><td> <input type="text" name="45" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-singlet production of <i>3S1</i> bottomonium states via 
<i>g g &rarr; bbbar[3S1(1)] g</i>. 
Code 501. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:gg2bbbar(3S1)[3S1(8)]g  </td><td></td><td> <input type="text" name="46" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> bottomonium states via 
<i>g g &rarr; bbbar[3S1(8)] g</i>. 
Code 502. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qg2bbbar(3S1)[3S1(8)]q  </td><td></td><td> <input type="text" name="47" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> bottomonium states via 
<i>q g &rarr; bbbar[3S1(8)] q</i>. 
Code 503. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qqbar2bbbar(3S1)[3S1(8)]g  </td><td></td><td> <input type="text" name="48" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> bottomonium states via 
<i>q qbar &rarr; bbbar[3S1(8)] g</i>. 
Code 504. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:gg2bbbar(3S1)[1S0(8)]g  </td><td></td><td> <input type="text" name="49" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> bottomonium states via 
<i>g g &rarr; bbbar[1S0(8)] g</i>. 
Code 505. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qg2bbbar(3S1)[1S0(8)]q  </td><td></td><td> <input type="text" name="50" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> bottomonium states via 
<i>q g &rarr; bbbar[1S0(8)] q</i>. 
Code 506. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qqbar2bbbar(3S1)[1S0(8)]g  </td><td></td><td> <input type="text" name="51" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> bottomonium states via 
<i>q qbar &rarr; bbbar[1S0(8)] g</i>. 
Code 507. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:gg2bbbar(3S1)[3PJ(8)]g  </td><td></td><td> <input type="text" name="52" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> bottomonium states via 
<i>g g &rarr; bbbar[3PJ(8)] g</i>. 
Code 508. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qg2bbbar(3S1)[3PJ(8)]q  </td><td></td><td> <input type="text" name="53" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> bottomonium states via 
<i>q g &rarr; bbbar[3PJ(8)] q</i>. 
Code 509. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qqbar2bbbar(3S1)[3PJ(8)]g  </td><td></td><td> <input type="text" name="54" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3S1</i> bottomonium states via 
<i>q qbar &rarr; bbbar[3SJ(8)] g</i>. 
Code 510. 
   
 
<h3>Bottomonium 3PJ States</h3> 
 
<b>Warning</b>: changed <code>fvec</code>, <code>mvec</code> or 
<code>pvec</code> values must be provided as a comma-separated list 
with the right number of elements, without any blanks inside the list. 
 
<br/><br/><table><tr><td><strong>Bottomonium:states(3PJ)  </td><td></td><td> <input type="text" name="55" value="10551,20553,555" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10551,20553,555</strong></code>)</td></tr></table>
The <i>3PJ</i> bottomonium states that can be produced from the following 
processes. Note that all vectors within this section, 
either of flags or parameters, must be the same length as this 
vector. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:O(3PJ)[3P0(1)] </td><td></td><td> <input type="text" name="56" value="0.085,0.085,0.085" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.085,0.085,0.085</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The color-singlet long-distance matrix elements 
<i>&lt;O[3P0(1)]&gt;/m_Q^2</i> for the <i>3PJ</i> bottomonium 
states. The remaining <i>&lt;O[3PJ(1)]&gt;/m_Q^2</i> 
are calculated from these long-distance matrix elements. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:O(3PJ)[3S1(8)] </td><td></td><td> <input type="text" name="57" value="0.04,0.04,0.04" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.04,0.04,0.04</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The color-singlet long-distance matrix elements 
<i>&lt;O[3S1(8)]&gt;</i> for the <i>3PJ</i> bottomonium states. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:gg2bbbar(3PJ)[3PJ(1)]g  </td><td></td><td> <input type="text" name="58" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-singlet production of <i>3PJ</i> bottomonium states via 
<i>g g &rarr; bbbar[3PJ(1)] g</i>. 
Code 511. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qg2bbbar(3PJ)[3PJ(1)]q  </td><td></td><td> <input type="text" name="59" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-singlet production of <i>3PJ</i> bottomonium states via 
<i>q g &rarr; bbbar[3PJ(1)] q</i>. 
Code 512. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qqbar2bbbar(3PJ)[3PJ(1)]g  </td><td></td><td> <input type="text" name="60" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-singlet production of <i>3PJ</i> bottomonium states via 
<i>q qbar &rarr; bbbar[3PJ(1)] g</i>. 
Code 513. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:gg2bbbar(3PJ)[3S1(8)]g  </td><td></td><td> <input type="text" name="61" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3PJ</i> bottomonium states via 
<i>g g &rarr; bbbar[3S1(8)] g</i>. 
Code 514. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qg2bbbar(3PJ)[3S1(8)]q  </td><td></td><td> <input type="text" name="62" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3PJ</i> bottomonium states via 
<i>q g &rarr; bbbar[3S1(8)] q</i>. 
Code 515. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qqbar2bbbar(3PJ)[3S1(8)]g  </td><td></td><td> <input type="text" name="63" value="off,off,off" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>off,off,off</strong></code>)</td></tr></table>
Colour-octet production of <i>3PJ</i> bottomonium states via 
<i>q qbar &rarr; bbbar[3S1(8)] g</i>. 
Code 516. 
   
 
<h3>Bottomonium 3DJ States</h3> 
 
<b>Warning</b>: changed <code>fvec</code>, <code>mvec</code> or 
<code>pvec</code> values must be provided as a comma-separated list 
with the right number of elements, without any blanks inside the list. 
 
<br/><br/><table><tr><td><strong>Bottomonium:states(3DJ)  </td><td></td><td> <input type="text" name="64" value="" size="20"/> </td></tr></table>
The <i>3DJ</i> bottomonium states that can be produced from the following 
processes. Currently, no <i>3DJ</i> states are included in the 
default <code>ParticleData</code> and so none are included here. Note 
that all vectors within this section, either of flags or parameters, 
must be the same length as this vector. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:O(3DJ)[3D1(1)] </td><td></td><td> <input type="text" name="65" value="" size="20"/>  &nbsp;&nbsp;(; <code>minimum = 0.0</code>)</td></tr></table>
The color-singlet long-distance matrix elements 
<i>&lt;O[3D1(1)]&gt;</i> for the <i>3PJ</i> bottomonium 
states. For a <i>3DJ</i> bottomonium state where <i>J</i> is 
not <i>1</i> the long distance matrix 
element <i>&lt;O[3DJ(1)]&gt;</i> is calculated 
by <i>(2J+1)&lt;O[3D1(1)]/3&gt;</i> using leading order spin symmetry 
relations. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:O(3DJ)[3P0(8)] </td><td></td><td> <input type="text" name="66" value="" size="20"/>  &nbsp;&nbsp;(; <code>minimum = 0.0</code>)</td></tr></table>
The colour-octet long-distance matrix 
elements <i>&lt;O[3P0(8)]&gt;/m_Q^2</i> for the 3DJ bottomonium 
states. The remaining <i>&lt;O[3PJ(8)]&gt;/m_Q^2</i> 
are calculated from these long-distance matrix elements. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:gg2bbbar(3DJ)[3DJ(1)]g  </td><td></td><td> <input type="text" name="67" value="" size="20"/> </td></tr></table>
Colour-singlet production of <i>3PJ</i> bottomonium states via 
<i>g g &rarr; bbbar[3DJ(1)] g</i>. 
Code 517. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:gg2bbbar(3DJ)[3PJ(8)]g  </td><td></td><td> <input type="text" name="68" value="" size="20"/> </td></tr></table>
Colour-octet production of <i>3DJ</i> bottomonium states via 
<i>g g &rarr; bbbar[3PJ(8)] g</i>. 
Code 518. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qg2bbbar(3DJ)[3PJ(8)]q  </td><td></td><td> <input type="text" name="69" value="" size="20"/> </td></tr></table>
Colour-octet production of <i>3DJ</i> bottomonium states via 
<i>q g &rarr; bbbar[3PJ(8)] q</i>. 
Code 519. 
   
 
<br/><br/><table><tr><td><strong>Bottomonium:qqbar2bbbar(3DJ)[3PJ(8)]g  </td><td></td><td> <input type="text" name="70" value="" size="20"/> </td></tr></table>
Colour-octet production of <i>3DJ</i> bottomonium states via 
<i>q qbar &rarr; bbbar[3PJ(8)] g</i>. 
Code 520. 
   
 
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
$data = "Onia:forceMassSplit = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0.2")
{
$data = "Onia:massSplit = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "Onia:all = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "Onia:all(3S1) = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "Onia:all(3PJ) = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "Onia:all(3DJ) = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "Charmonium:all = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "Bottomonium:all = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "443,100443")
{
$data = "Charmonium:states(3S1) = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "1.16,0.76")
{
$data = "Charmonium:O(3S1)[3S1(1)] = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "0.0119,0.0050")
{
$data = "Charmonium:O(3S1)[3S1(8)] = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "0.01,0.004")
{
$data = "Charmonium:O(3S1)[1S0(8)] = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0.01,0.004")
{
$data = "Charmonium:O(3S1)[3P0(8)] = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off,off")
{
$data = "Charmonium:gg2ccbar(3S1)[3S1(1)]g = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off,off")
{
$data = "Charmonium:gg2ccbar(3S1)[3S1(8)]g = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off,off")
{
$data = "Charmonium:qg2ccbar(3S1)[3S1(8)]q = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off,off")
{
$data = "Charmonium:qqbar2ccbar(3S1)[3S1(8)]g = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "off,off")
{
$data = "Charmonium:gg2ccbar(3S1)[1S0(8)]g = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off,off")
{
$data = "Charmonium:qg2ccbar(3S1)[1S0(8)]q = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off,off")
{
$data = "Charmonium:qqbar2ccbar(3S1)[1S0(8)]g = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "off,off")
{
$data = "Charmonium:gg2ccbar(3S1)[3PJ(8)]g = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "off,off")
{
$data = "Charmonium:qg2ccbar(3S1)[3PJ(8)]q = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "off,off")
{
$data = "Charmonium:qqbar2ccbar(3S1)[3PJ(8)]g = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "10441,20443,445")
{
$data = "Charmonium:states(3PJ) = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "0.05,0.05,0.05")
{
$data = "Charmonium:O(3PJ)[3P0(1)] = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "0.0031,0.0031,0.0031")
{
$data = "Charmonium:O(3PJ)[3S1(8)] = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "off,off,off")
{
$data = "Charmonium:gg2ccbar(3PJ)[3PJ(1)]g = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "off,off,off")
{
$data = "Charmonium:qg2ccbar(3PJ)[3PJ(1)]q = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "off,off,off")
{
$data = "Charmonium:qqbar2ccbar(3PJ)[3PJ(1)]g = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "off,off,off")
{
$data = "Charmonium:gg2ccbar(3PJ)[3S1(8)]g = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "off,off,off")
{
$data = "Charmonium:qg2ccbar(3PJ)[3S1(8)]q = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "off,off,off")
{
$data = "Charmonium:qqbar2ccbar(3PJ)[3S1(8)]g = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "30443")
{
$data = "Charmonium:states(3DJ) = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "0.161")
{
$data = "Charmonium:O(3DJ)[3D1(1)] = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "0.01")
{
$data = "Charmonium:O(3DJ)[3P0(8)] = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "off")
{
$data = "Charmonium:gg2ccbar(3DJ)[3DJ(1)]g = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "off")
{
$data = "Charmonium:gg2ccbar(3DJ)[3PJ(8)]g = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "off")
{
$data = "Charmonium:qg2ccbar(3DJ)[3PJ(8)]q = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "off")
{
$data = "Charmonium:qqbar2ccbar(3DJ)[3PJ(8)]g = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
if($_POST["40"] != "553,100553,200553")
{
$data = "Bottomonium:states(3S1) = ".$_POST["40"]."\n";
fwrite($handle,$data);
}
if($_POST["41"] != "9.28,4.63,3.54")
{
$data = "Bottomonium:O(3S1)[3S1(1)] = ".$_POST["41"]."\n";
fwrite($handle,$data);
}
if($_POST["42"] != "0.15,0.045,0.075")
{
$data = "Bottomonium:O(3S1)[3S1(8)] = ".$_POST["42"]."\n";
fwrite($handle,$data);
}
if($_POST["43"] != "0.02,0.06,0.1")
{
$data = "Bottomonium:O(3S1)[1S0(8)] = ".$_POST["43"]."\n";
fwrite($handle,$data);
}
if($_POST["44"] != "0.02,0.06,0.1")
{
$data = "Bottomonium:O(3S1)[3P0(8)] = ".$_POST["44"]."\n";
fwrite($handle,$data);
}
if($_POST["45"] != "off,off,off")
{
$data = "Bottomonium:gg2bbbar(3S1)[3S1(1)]g = ".$_POST["45"]."\n";
fwrite($handle,$data);
}
if($_POST["46"] != "off,off,off")
{
$data = "Bottomonium:gg2bbbar(3S1)[3S1(8)]g = ".$_POST["46"]."\n";
fwrite($handle,$data);
}
if($_POST["47"] != "off,off,off")
{
$data = "Bottomonium:qg2bbbar(3S1)[3S1(8)]q = ".$_POST["47"]."\n";
fwrite($handle,$data);
}
if($_POST["48"] != "off,off,off")
{
$data = "Bottomonium:qqbar2bbbar(3S1)[3S1(8)]g = ".$_POST["48"]."\n";
fwrite($handle,$data);
}
if($_POST["49"] != "off,off,off")
{
$data = "Bottomonium:gg2bbbar(3S1)[1S0(8)]g = ".$_POST["49"]."\n";
fwrite($handle,$data);
}
if($_POST["50"] != "off,off,off")
{
$data = "Bottomonium:qg2bbbar(3S1)[1S0(8)]q = ".$_POST["50"]."\n";
fwrite($handle,$data);
}
if($_POST["51"] != "off,off,off")
{
$data = "Bottomonium:qqbar2bbbar(3S1)[1S0(8)]g = ".$_POST["51"]."\n";
fwrite($handle,$data);
}
if($_POST["52"] != "off,off,off")
{
$data = "Bottomonium:gg2bbbar(3S1)[3PJ(8)]g = ".$_POST["52"]."\n";
fwrite($handle,$data);
}
if($_POST["53"] != "off,off,off")
{
$data = "Bottomonium:qg2bbbar(3S1)[3PJ(8)]q = ".$_POST["53"]."\n";
fwrite($handle,$data);
}
if($_POST["54"] != "off,off,off")
{
$data = "Bottomonium:qqbar2bbbar(3S1)[3PJ(8)]g = ".$_POST["54"]."\n";
fwrite($handle,$data);
}
if($_POST["55"] != "10551,20553,555")
{
$data = "Bottomonium:states(3PJ) = ".$_POST["55"]."\n";
fwrite($handle,$data);
}
if($_POST["56"] != "0.085,0.085,0.085")
{
$data = "Bottomonium:O(3PJ)[3P0(1)] = ".$_POST["56"]."\n";
fwrite($handle,$data);
}
if($_POST["57"] != "0.04,0.04,0.04")
{
$data = "Bottomonium:O(3PJ)[3S1(8)] = ".$_POST["57"]."\n";
fwrite($handle,$data);
}
if($_POST["58"] != "off,off,off")
{
$data = "Bottomonium:gg2bbbar(3PJ)[3PJ(1)]g = ".$_POST["58"]."\n";
fwrite($handle,$data);
}
if($_POST["59"] != "off,off,off")
{
$data = "Bottomonium:qg2bbbar(3PJ)[3PJ(1)]q = ".$_POST["59"]."\n";
fwrite($handle,$data);
}
if($_POST["60"] != "off,off,off")
{
$data = "Bottomonium:qqbar2bbbar(3PJ)[3PJ(1)]g = ".$_POST["60"]."\n";
fwrite($handle,$data);
}
if($_POST["61"] != "off,off,off")
{
$data = "Bottomonium:gg2bbbar(3PJ)[3S1(8)]g = ".$_POST["61"]."\n";
fwrite($handle,$data);
}
if($_POST["62"] != "off,off,off")
{
$data = "Bottomonium:qg2bbbar(3PJ)[3S1(8)]q = ".$_POST["62"]."\n";
fwrite($handle,$data);
}
if($_POST["63"] != "off,off,off")
{
$data = "Bottomonium:qqbar2bbbar(3PJ)[3S1(8)]g = ".$_POST["63"]."\n";
fwrite($handle,$data);
}
if($_POST["64"] != "")
{
$data = "Bottomonium:states(3DJ) = ".$_POST["64"]."\n";
fwrite($handle,$data);
}
if($_POST["65"] != "")
{
$data = "Bottomonium:O(3DJ)[3D1(1)] = ".$_POST["65"]."\n";
fwrite($handle,$data);
}
if($_POST["66"] != "")
{
$data = "Bottomonium:O(3DJ)[3P0(8)] = ".$_POST["66"]."\n";
fwrite($handle,$data);
}
if($_POST["67"] != "")
{
$data = "Bottomonium:gg2bbbar(3DJ)[3DJ(1)]g = ".$_POST["67"]."\n";
fwrite($handle,$data);
}
if($_POST["68"] != "")
{
$data = "Bottomonium:gg2bbbar(3DJ)[3PJ(8)]g = ".$_POST["68"]."\n";
fwrite($handle,$data);
}
if($_POST["69"] != "")
{
$data = "Bottomonium:qg2bbbar(3DJ)[3PJ(8)]q = ".$_POST["69"]."\n";
fwrite($handle,$data);
}
if($_POST["70"] != "")
{
$data = "Bottomonium:qqbar2bbbar(3DJ)[3PJ(8)]g = ".$_POST["70"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
