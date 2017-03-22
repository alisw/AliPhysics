<html>
<head>
<title>Extra-Dimensional Processes</title>
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

<form method='post' action='ExtraDimensionalProcesses.php'>

<h2>Extra-Dimensional Processes</h2>

Scenarios with extra dimensions (ED) allow a multitude of processes.
Currently three different categories of processes are implemented. 
The first involves the production of excited Kaluza Klein states 
within so-called Randall-Sundrum (RS) scenarios, the second is 
related to resonance production in TeV-1 sized extra dimensions 
and the third relates to phenomena from large extra dimensions (LED). 
Due to the close relation between the LED model and a so-called 
unparticle model, similar unparticle processes are also kept in this 
section.

<h3>Randall-Sundrum Resonances, production processes</h3>

The graviton (G*) and gluon (KKgluon*) resonance states are assigned 
PDG code 5100039 and 5100021 respectively. The G* processes are 
described in [<a href="Bibliography.php" target="page">Bij01</a>] and the KKgluon* process in [<a href="Bibliography.php" target="page">Ask11</a>]. 
Decays into fermion and boson pairs are handled with the correct 
angular distributions, while subsequent decays are handled 
isotropically.

<p/>
There are two lowest-order processes that together normally 
should be sufficient for a simulation of <i>G^*</i> production. 

<br/><br/><strong>ExtraDimensionsG*:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of lowest-order <i>G^*</i> production
processes, i.e. the two ones below.
  

<br/><br/><strong>ExtraDimensionsG*:gg2G*</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> G^*</i>. 
Code 5001.
  

<br/><br/><strong>ExtraDimensionsG*:ffbar2G*</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> G^*</i>. 
Code 5002.
  

<p/>
In addition there are three first-order processes included. These 
are of less interest, but can be used for dedicated studies of the 
high-<i>pT</i> tail of <i>G^*</i> production. As usual, it would 
be double counting to include the lowest-order and first-order 
processes simultaneously. Therefore the latter ones are not included 
with the <code>ExtraDimensionsG*:all = on</code> option. In this set 
of processes all decay angles are assumed isotropic.

<br/><br/><strong>ExtraDimensionsG*:gg2G*g</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> G^* g</i>. 
Code 5003.
  

<br/><br/><strong>ExtraDimensionsG*:qg2G*q</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q g -> G^* q</i>. 
Code 5004.
  

<br/><br/><strong>ExtraDimensionsG*:qqbar2G*g</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar -> G^* g</i>. 
Code 5005.
  

<p/>
There is also one process for the production of a gluon resonance.

<br/><br/><strong>ExtraDimensionsG*:qqbar2KKgluon*</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar -> g^*/KKgluon^*</i>. 
Code 5006.
  

<h3>Randall-Sundrum Resonances, parameters</h3>

In the above scenario the main free parameters are the masses, which 
are set as usual. In addition there are the following coupling parameters. 
The coupling <i>kappaMG</i> follows the conventions in [<a href="Bibliography.php" target="page">Bij01</a>], 
where as the flavour dependent couplings follow the conventions used in 
[<a href="Bibliography.php" target="page">Dav01</a>].

<br/><br/><strong>ExtraDimensionsG*:SMinBulk</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Parameter to choose between the two scenarios: 
<i>off</i>, SM on the TeV brane (common <i>kappaMG</i> coupling); 
<i>on</i>, SM in the ED bulk (flavour dependent couplings). 
This parameter is only relevant for the lowest-order graviton 
(<i>G*</i>) processes, where as the first-order processes 
corresponds to the <i>off</i> scenario.
  

<br/><br/><strong>ExtraDimensionsG*:VLVL</strong>  <input type="radio" name="9" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="9" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Parameter to specify Z/W coupling scenario: 
<i>off</i>, usual Z/W boson couplings; 
<i>on</i>, coupling only to longitudinal Z/W bosons. 
In both cases the <i>GZZ</i> and <i>GWW</i> values are used 
and this parameter is only relevant when <i>SMinBulk = on</i>. 
The formulas for longitudinal bosons should be appropriate up to 
<i>O(m_V/E_V)</i> corrections.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:kappaMG </td><td></td><td> <input type="text" name="10" value="0.054" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.054</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
dimensionless coupling, which enters quadratically in all partial
widths of the <i>G^*</i>. Is 
<i>kappa m_G* = sqrt(2) x_1 k / Mbar_Pl</i>,
where <i>x_1 = 3.83</i> is the first zero of the <i>J_1</i> Bessel 
function and <i>Mbar_Pl</i> is the modified Planck mass.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:Gll </td><td></td><td> <input type="text" name="11" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Coupling between graviton and leptons.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:Gqq </td><td></td><td> <input type="text" name="12" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Coupling between graviton and light quarks.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:Gbb </td><td></td><td> <input type="text" name="13" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Coupling between graviton and bottom quark.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:Gtt </td><td></td><td> <input type="text" name="14" value="0.001" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.001</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Coupling between graviton and top quark.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:Ggg </td><td></td><td> <input type="text" name="15" value="0.000013" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.000013</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Coupling between graviton and gluon.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:Ggmgm </td><td></td><td> <input type="text" name="16" value="0.000013" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.000013</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Coupling between graviton and gamma.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:GZZ </td><td></td><td> <input type="text" name="17" value="0.001" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.001</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Coupling between graviton and Z boson.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:GWW </td><td></td><td> <input type="text" name="18" value="0.001" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.001</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Coupling between graviton and W boson.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:Ghh </td><td></td><td> <input type="text" name="19" value="0.001" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.001</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Coupling between graviton and Higgs bosons.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:KKgqR </td><td></td><td> <input type="text" name="20" value="-0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-0.2</strong></code>)</td></tr></table>
Coupling between KK-gluon and a right-handed light quark.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:KKgqL </td><td></td><td> <input type="text" name="21" value="-0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-0.2</strong></code>)</td></tr></table>
Coupling between KK-gluon and a left-handed light quark.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:KKgbR </td><td></td><td> <input type="text" name="22" value="-0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-0.2</strong></code>)</td></tr></table>
Coupling between KK-gluon and a right-handed bottom quark.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:KKgbL </td><td></td><td> <input type="text" name="23" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>)</td></tr></table>
Coupling between KK-gluon and a left-handed bottom quark.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:KKgtR </td><td></td><td> <input type="text" name="24" value="4.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4.0</strong></code>)</td></tr></table>
Coupling between KK-gluon and a right-handed top quark.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:KKgtL </td><td></td><td> <input type="text" name="25" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>)</td></tr></table>
Coupling between KK-gluon and a left-handed top quark.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:KKintMode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Choice of full <ei>g^*/KK-gluon^*</ei> structure or not in relevant 
processes.
<br/>
<input type="radio" name="26" value="0" checked="checked"><strong>0 </strong>: full <ei>g^*/KK-gluon^*</ei> structure, with  interference included.<br/>
<input type="radio" name="26" value="1"><strong>1 </strong>: only pure <ei>gluon_{SM}</ei> contribution.<br/>
<input type="radio" name="26" value="2"><strong>2 </strong>: only pure <ei>gluon_{KK}</ei> contribution.<br/>

<h3>TeV^-1 Sized Extra Dimension, production processes</h3>

This section contains a processes involving the production
of electroweak KK gauge bosons, i.e. <i>gamma_{KK}/Z_{KK}</i>,
in one TeV^-1 sized extra dimension. The process is described 
in [<a href="Bibliography.php" target="page">Bel10</a>] and allows for individual final states to be 
specified.

<br/><br/><strong>ExtraDimensionsTEV:ffbar2ddbar</strong>  <input type="radio" name="27" value="on"><strong>On</strong>
<input type="radio" name="27" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> d dbar </i>,
Code 5061.
  

<br/><br/><strong>ExtraDimensionsTEV:ffbar2uubar</strong>  <input type="radio" name="28" value="on"><strong>On</strong>
<input type="radio" name="28" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> u ubar </i>,
Code 5062.
  

<br/><br/><strong>ExtraDimensionsTEV:ffbar2ssbar</strong>  <input type="radio" name="29" value="on"><strong>On</strong>
<input type="radio" name="29" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> s sbar </i>,
Code 5063.
  

<br/><br/><strong>ExtraDimensionsTEV:ffbar2ccbar</strong>  <input type="radio" name="30" value="on"><strong>On</strong>
<input type="radio" name="30" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> c cbar </i>,
Code 5064.
  

<br/><br/><strong>ExtraDimensionsTEV:ffbar2bbbar</strong>  <input type="radio" name="31" value="on"><strong>On</strong>
<input type="radio" name="31" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> b bbar </i>,
Code 5065.
  

<br/><br/><strong>ExtraDimensionsTEV:ffbar2ttbar</strong>  <input type="radio" name="32" value="on"><strong>On</strong>
<input type="radio" name="32" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> t tbar </i>,
Code 5066.
  

<br/><br/><strong>ExtraDimensionsTEV:ffbar2e+e-</strong>  <input type="radio" name="33" value="on"><strong>On</strong>
<input type="radio" name="33" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> e+ e- </i>,
Code 5071.
  

<br/><br/><strong>ExtraDimensionsTEV:ffbar2nuenuebar</strong>  <input type="radio" name="34" value="on"><strong>On</strong>
<input type="radio" name="34" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> nue nuebar </i>,
Code 5072.
  

<br/><br/><strong>ExtraDimensionsTEV:ffbar2mu+mu-</strong>  <input type="radio" name="35" value="on"><strong>On</strong>
<input type="radio" name="35" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> mu+ mu- </i>,
Code 5073.
  

<br/><br/><strong>ExtraDimensionsTEV:ffbar2numunumubar</strong>  <input type="radio" name="36" value="on"><strong>On</strong>
<input type="radio" name="36" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> numu numubar </i>,
Code 5074.
  

<br/><br/><strong>ExtraDimensionsTEV:ffbar2tau+tau-</strong>  <input type="radio" name="37" value="on"><strong>On</strong>
<input type="radio" name="37" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> tau+ tau- </i>,
Code 5075.
  

<br/><br/><strong>ExtraDimensionsTEV:ffbar2nutaunutaubar</strong>  <input type="radio" name="38" value="on"><strong>On</strong>
<input type="radio" name="38" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (gamma_{KK}/Z_{KK}) -> nutau nutaubar </i>,
Code 5076.
  

<h3>TeV^-1 Sized Extra Dimension, parameters</h3>

Irrespective of the parameter options used, the particle produced, 
<i>gamma_{KK}/Z_{KK}</i>, will always be assigned code 5000023.

<br/><br/><table><tr><td><strong>ExtraDimensionsTEV:gmZmode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
Choice of full <ei>gamma_{KK}/Z_{KK}</ei> structure or not in relevant 
processes.
<br/>
<input type="radio" name="39" value="0"><strong>0 </strong>: full <ei>gamma_{SM}/Z_{SM}</ei> structure, with  interference included.<br/>
<input type="radio" name="39" value="1"><strong>1 </strong>: only pure <ei>gamma_{SM}</ei> contribution.<br/>
<input type="radio" name="39" value="2"><strong>2 </strong>: only pure <ei>Z_{SM}</ei> contribution.<br/>
<input type="radio" name="39" value="3" checked="checked"><strong>3 </strong>: full <ei>gamma_{KK}/Z_{KK}</ei> structure, with  interference included.<br/>
<input type="radio" name="39" value="4"><strong>4 </strong>: only pure <ei>gamma_{KK}</ei> contribution, with  SM interference included.<br/>
<input type="radio" name="39" value="5"><strong>5 </strong>: only pure <ei>Z_{KK}</ei> contribution, with SM  interference included.<br/>

<br/><br/><table><tr><td><strong>ExtraDimensionsTEV:nMax </td><td></td><td> <input type="text" name="40" value="10" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10</strong></code>; <code>minimum = 1</code>; <code>maximum = 100</code>)</td></tr></table>
The number of included KK excitations.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsTEV:mStar </td><td></td><td> <input type="text" name="41" value="4000.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4000.0</strong></code>; <code>minimum = 1000.0</code>)</td></tr></table>
The KK mass <i>m^*</i>, given by the inverse of the single extra 
dimension radius.
  

<h3>Large Extra Dimensions, production processes</h3>

The LED graviton, where the KK-modes normally are summed and do not 
give rise to phenomena individually, is assigned PDG code 5000039. 
The graviton emission and virtual graviton exchange processes use 
the same implementation as the corresponding unparticle processes, 
which are all described in [<a href="Bibliography.php" target="page">Ask10</a>]. It is also possible to 
generate monojet events from scalar graviton emission as described 
in [<a href="Bibliography.php" target="page">Azu05</a>], by turning on the option <i>GravScalar</i>.

<p/>
<i>Note:</i> As discussed in [<a href="Bibliography.php" target="page">Ask09</a>], for the graviton or 
unparticle emission processes the underlying Breit-Wigner mass 
distribution should be matched to the graviton mass spectrum in order 
to achieve an optimal MC efficiency. 

<p/>
The following lowest order graviton emission processes are available.

<br/><br/><strong>ExtraDimensionsLED:monojet</strong>  <input type="radio" name="42" value="on"><strong>On</strong>
<input type="radio" name="42" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of lowest-order <i>G jet</i> emission
processes, i.e. the three ones below.
  

<br/><br/><strong>ExtraDimensionsLED:gg2Gg</strong>  <input type="radio" name="43" value="on"><strong>On</strong>
<input type="radio" name="43" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> G g</i>. 
Code 5021.
  

<br/><br/><strong>ExtraDimensionsLED:qg2Gq</strong>  <input type="radio" name="44" value="on"><strong>On</strong>
<input type="radio" name="44" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q g -> G q</i>. 
Code 5022.
  

<br/><br/><strong>ExtraDimensionsLED:qqbar2Gg</strong>  <input type="radio" name="45" value="on"><strong>On</strong>
<input type="radio" name="45" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar -> G g</i>. 
Code 5023.
  

<br/><br/><strong>ExtraDimensionsLED:ffbar2GZ</strong>  <input type="radio" name="46" value="on"><strong>On</strong>
<input type="radio" name="46" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> G Z</i>. 
Code 5024.
  

<br/><br/><strong>ExtraDimensionsLED:ffbar2Ggamma</strong>  <input type="radio" name="47" value="on"><strong>On</strong>
<input type="radio" name="47" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> G gamma</i>. This process corresponds 
to the photon limit of the <i>G Z</i> process, as described in 
[<a href="Bibliography.php" target="page">Ask09</a>].
Code 5025.
  

<p/>
The following LED processes with virtual graviton exchange are 
available.

<br/><br/><strong>ExtraDimensionsLED:ffbar2gammagamma</strong>  <input type="radio" name="48" value="on"><strong>On</strong>
<input type="radio" name="48" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (LED G*) -> gamma gamma</i>. If the 
graviton contribution is zero, the results corresponds to the 
SM contribution, i.e. equivalent to 
<code>PromptPhoton:ffbar2gammagamma</code>.
Code 5026.
  

<br/><br/><strong>ExtraDimensionsLED:gg2gammagamma</strong>  <input type="radio" name="49" value="on"><strong>On</strong>
<input type="radio" name="49" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> (LED G*) -> gamma gamma</i>. 
Code 5027.
  

<br/><br/><strong>ExtraDimensionsLED:ffbar2llbar</strong>  <input type="radio" name="50" value="on"><strong>On</strong>
<input type="radio" name="50" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (LED G*) -> l l </i>, where 
<i>l</i> is a charged lepton. If the graviton contribution 
is zero, the results corresponds to the SM contribution, i.e. 
similar to <code>WeakSingleBoson:ffbar2gmZ</code>. Does not 
include t-channel amplitude relevant for e^+e^- to e^+e^- 
and no K-factor is used.
Code 5028.
  

<br/><br/><strong>ExtraDimensionsLED:gg2llbar</strong>  <input type="radio" name="51" value="on"><strong>On</strong>
<input type="radio" name="51" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> (LED G*) -> l l</i>. 
Code 5029.
  

<p/> 
Dijet production including graviton exchange is also available, using 
the same effective theory approach as the LED G exchange processes 
above or including more detailed amplitudes in accordance with 
[<a href="Bibliography.php" target="page">Fra11</a>]. In case of the latter, the value of <i>LambdaT</i> 
is used as the value of the cut-off scale <i>Lambda</i>. For this 
reason the dijet processes only relates to the LED model and no 
unparticle versions are available. The processes are grouped together 
like their <i>HardQCD</i> equivalents and should therefore converge 
to the same results in the limit of an insignificant graviton 
contribution.

<p/>
<i>Warning:</i> These LED dijets processes are still being validated. 

<br/><br/><strong>ExtraDimensionsLED:dijets</strong>  <input type="radio" name="52" value="on"><strong>On</strong>
<input type="radio" name="52" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of lowest-order <i>jet jet</i> 
production processes with graviton exchange, i.e. the six ones 
below.
  

<br/><br/><strong>ExtraDimensionsLED:gg2DJgg</strong>  <input type="radio" name="53" value="on"><strong>On</strong>
<input type="radio" name="53" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> (LED G*) -> g g</i>. 
Code 5030.
  

<br/><br/><strong>ExtraDimensionsLED:gg2DJqqbar</strong>  <input type="radio" name="54" value="on"><strong>On</strong>
<input type="radio" name="54" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> (LED G*) -> q qbar</i>. Number of 
outgoing flavours specified by <i>nQuarkNew</i> parameter 
below.
Code 5031.
  

<br/><br/><strong>ExtraDimensionsLED:qg2DJqg</strong>  <input type="radio" name="55" value="on"><strong>On</strong>
<input type="radio" name="55" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q g -> (LED G*) -> q g</i> and 
<i>qbar g -> (LED G*) -> qbar g</i>. 
Code 5032.
  

<br/><br/><strong>ExtraDimensionsLED:qq2DJqq</strong>  <input type="radio" name="56" value="on"><strong>On</strong>
<input type="radio" name="56" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q q(bar)' -> (LED G*) -> q q(bar)'</i>. 
Including <i>q</i> and <i>qbar</i> of same or different 
flavours, but the outgoing flavours equals the incoming ones. 
Code 5033.
  

<br/><br/><strong>ExtraDimensionsLED:qqbar2DJgg</strong>  <input type="radio" name="57" value="on"><strong>On</strong>
<input type="radio" name="57" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar -> (LED G*) -> g g</i>. 
Code 5034.
  

<br/><br/><strong>ExtraDimensionsLED:qqbar2DJqqbarNew</strong>  <input type="radio" name="58" value="on"><strong>On</strong>
<input type="radio" name="58" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar -> (LED G*) -> q' qbar'</i>. Number of 
outgoing flavours specified by <i>nQuarkNew</i> parameter below.
Code 5035.
  

<h3>Large Extra Dimensions, parameters</h3>

<br/><br/><strong>ExtraDimensionsLED:GravScalar</strong>  <input type="radio" name="59" value="on"><strong>On</strong>
<input type="radio" name="59" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow the monojet processes to produce scalar graviton emission
instead of the default tensor one. The scalar option is according 
to the processes described in [<a href="Bibliography.php" target="page">Azu05</a>] and includes two 
coupling constants below.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:n  </td><td></td><td> <input type="text" name="60" value="2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>)</td></tr></table>
Number of extra dimensions.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:MD </td><td></td><td> <input type="text" name="61" value="2000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2000.</strong></code>; <code>minimum = 100.0</code>)</td></tr></table>
Fundamental scale of gravity in <i>D = 4 + n</i> dimensions.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:LambdaT </td><td></td><td> <input type="text" name="62" value="2000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2000.</strong></code>; <code>minimum = 100.0</code>)</td></tr></table>
Ultraviolet cutoff parameter for the virtual graviton exchange processes.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:NegInt  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
Allows to change sign of the interference terms in the graviton exchange 
processes, common in connection to using the <ei>Hewett</ei> convention 
<ref>Hew99</ref>. 
<br/>
<input type="radio" name="63" value="0" checked="checked"><strong>0 </strong>: 1<br/>
<input type="radio" name="63" value="1"><strong>1 </strong>: -1<br/>

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:CutOffMode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
Options for when the hard scale of the process (e.g. <ei>sHat</ei>) 
approaches or exceed the scale of validity of the low energy effective 
theory (e.g. <ei>M_D</ei>). <ei>Note:</ei> Option 1 only concerns the 
graviton emission processes and the form factor is currently not available 
for the scalar graviton processes.
<br/>
<input type="radio" name="64" value="0" checked="checked"><strong>0 </strong>: Do nothing, i.e. all values of <ei>sHat</ei> contribute. <br/>
<input type="radio" name="64" value="1"><strong>1 </strong>: Truncate contributing <ei>sHat</ei> region  (<ref>Ask09</ref>).<br/>
<input type="radio" name="64" value="2"><strong>2 </strong>: Form factor, using <ei>mu = renormScale2</ei> .<br/>
<input type="radio" name="64" value="3"><strong>3 </strong>: Form factor, using <ei>mu = E_jet</ei>.<br/>

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:t </td><td></td><td> <input type="text" name="65" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.001</code>)</td></tr></table>
Form factor parameter.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:g </td><td></td><td> <input type="text" name="66" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Coupling related to scalar graviton emission.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:c </td><td></td><td> <input type="text" name="67" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Coupling related to scalar graviton emission.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:nQuarkNew  </td><td></td><td> <input type="text" name="68" value="3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
Number of allowed outgoing new quark flavours in the above
<i>q qbar -> (LED G*) -> q' qbar'</i> and <i>g g -> (LED G*) -> q' qbar'</i> 
processes. Similar to <i>HardQCD:nQuarkNew</i> for the QCD processes.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:opMode  </td><td></td><td> <input type="text" name="69" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
Options to specify <i>S</i> function for LED dijet amplitudes.
<br/><code>option </code><strong> 0</strong> : Use detailed amplitude, as described in [<a href="Bibliography.php" target="page">Fra11</a>].  
<br/><code>option </code><strong> 1</strong> : Use conventional <i>LambdaT</i> parametrization, like the other LED processes.  
  

<h3>Unparticles, production processes</h3>

As mentioned above, the similar unparticle and graviton processes 
share the same implementations. The unparticle processes, however, 
only use the dedicated unparticle parameters below. The unparticle 
is also assigned the PDG code 5000039 and is therefore called 
<i>Graviton</i> in the event record. The graviton and unparticle 
emission as well as virtual graviton and unparticle exchange processes 
are described in [<a href="Bibliography.php" target="page">Ask10</a>].

<p/>
<i>Note:</i> As discussed in [<a href="Bibliography.php" target="page">Ask09</a>], for the graviton or 
unparticle emission processes the underlying Breit-Wigner mass 
distribution should be matched to the graviton mass spectrum in order 
to achieve an optimal MC efficiency. 

<p/>
The following unparticle emission processes are available.

<br/><br/><strong>ExtraDimensionsUnpart:monojet</strong>  <input type="radio" name="70" value="on"><strong>On</strong>
<input type="radio" name="70" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of lowest-order <i>U jet</i> emission
processes, i.e. the three ones below.
  

<br/><br/><strong>ExtraDimensionsUnpart:gg2Ug</strong>  <input type="radio" name="71" value="on"><strong>On</strong>
<input type="radio" name="71" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> U g</i>. 
Code 5045.
  

<br/><br/><strong>ExtraDimensionsUnpart:qg2Uq</strong>  <input type="radio" name="72" value="on"><strong>On</strong>
<input type="radio" name="72" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q g -> U q</i>. 
Code 5046.
  

<br/><br/><strong>ExtraDimensionsUnpart:qqbar2Ug</strong>  <input type="radio" name="73" value="on"><strong>On</strong>
<input type="radio" name="73" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar -> U g</i>. 
Code 5047.
  

<br/><br/><strong>ExtraDimensionsUnpart:ffbar2UZ</strong>  <input type="radio" name="74" value="on"><strong>On</strong>
<input type="radio" name="74" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> U Z</i>.  
Code 5041.
  

<br/><br/><strong>ExtraDimensionsUnpart:ffbar2Ugamma</strong>  <input type="radio" name="75" value="on"><strong>On</strong>
<input type="radio" name="75" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> U gamma</i>. This process corresponds 
to the photon limit of the <i>U Z</i> process, as described in 
[<a href="Bibliography.php" target="page">Ask09</a>].
Code 5042.
  

<p/>
The following processes with virtual unparticle exchange are available.

<br/><br/><strong>ExtraDimensionsUnpart:ffbar2gammagamma</strong>  <input type="radio" name="76" value="on"><strong>On</strong>
<input type="radio" name="76" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (U*) -> gamma gamma</i>. If the unparticle  
contribution is zero in the spin-2 case, the results corresponds to 
the SM contribution, i.e. equivalent to 
<code>PromptPhoton:ffbar2gammagamma</code>.
Code 5043.
  

<br/><br/><strong>ExtraDimensionsUnpart:gg2gammagamma</strong>  <input type="radio" name="77" value="on"><strong>On</strong>
<input type="radio" name="77" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> (U*) -> gamma gamma</i>. 
Code 5044.
  

<br/><br/><strong>ExtraDimensionsUnpart:ffbar2llbar</strong>  <input type="radio" name="78" value="on"><strong>On</strong>
<input type="radio" name="78" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (U*) -> l lbar </i>, where 
<i>l</i> is a charged lepton. If the unparticle contribution 
is zero, the results corresponds to the SM contribution, i.e. 
similar to <code>WeakSingleBoson:ffbar2gmZ</code>. Does not 
include t-channel amplitude relevant for e^+e^- to e^+e^- 
and no K-factor is used.
Code 5048.
  

<br/><br/><strong>ExtraDimensionsUnpart:gg2llbar</strong>  <input type="radio" name="79" value="on"><strong>On</strong>
<input type="radio" name="79" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> (U*) -> l lbar</i>. 
Code 5049.
  

<h3>Unparticles, parameters</h3>

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:spinU  </td><td></td><td> <input type="text" name="80" value="2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Unparticle spin.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:dU </td><td></td><td> <input type="text" name="81" value="1.4" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.4</strong></code>; <code>minimum = 1.0</code>)</td></tr></table>
Scale dimension parameter.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:LambdaU </td><td></td><td> <input type="text" name="82" value="2000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2000.</strong></code>; <code>minimum = 100.0</code>)</td></tr></table>
Unparticle renormalization scale.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:lambda </td><td></td><td> <input type="text" name="83" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Unparticle coupling to the SM fields.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:ratio </td><td></td><td> <input type="text" name="84" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 1.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
Ratio, <i>lambda'/lambda</i>, between the two possible coupling constants 
of the spin-2 ME. <b>Warning:</b> A <i>ratio</i> value different from one 
give rise to an IR divergence which makes the event generation very slow, so 
this values is fixed to <i>ratio = 1</i> for the moment.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:CutOffMode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
Options for when the hard scale of the process (e.g. <ei>sHat</ei>) 
approaches or exceed the scale of validity of the low energy effective 
theory (<ei>Lambda_U</ei>). This mode only concerns the unparticle 
emission processes. 
<br/>
<input type="radio" name="85" value="0" checked="checked"><strong>0 </strong>: Do nothing, i.e. all values of <ei>sHat</ei>  contribute.<br/>
<input type="radio" name="85" value="1"><strong>1 </strong>: Truncate contributing <ei>sHat</ei> region  (<ref>Ask09</ref>).<br/>

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:gXX  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Chiral unparticle couplings, <ei>gXX = gLL = gRR</ei>. Only relevant 
for lepton production from spin-1 unparticle exchange. 
<br/>
<input type="radio" name="86" value="0" checked="checked"><strong>0 </strong>: 1<br/>
<input type="radio" name="86" value="1"><strong>1 </strong>: -1<br/>
<input type="radio" name="86" value="2"><strong>2 </strong>: 0<br/>

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:gXY  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Chiral unparticle couplings, <ei>gXY = gLR = gRL</ei>. Only relevant 
for lepton production from spin-1 unparticle exchange. 
<br/>
<input type="radio" name="87" value="0" checked="checked"><strong>0 </strong>: 1<br/>
<input type="radio" name="87" value="1"><strong>1 </strong>: -1<br/>
<input type="radio" name="87" value="2"><strong>2 </strong>: 0<br/>

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
$data = "ExtraDimensionsG*:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "ExtraDimensionsG*:gg2G* = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "ExtraDimensionsG*:ffbar2G* = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "ExtraDimensionsG*:gg2G*g = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "ExtraDimensionsG*:qg2G*q = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "ExtraDimensionsG*:qqbar2G*g = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "ExtraDimensionsG*:qqbar2KKgluon* = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "ExtraDimensionsG*:SMinBulk = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "on")
{
$data = "ExtraDimensionsG*:VLVL = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "0.054")
{
$data = "ExtraDimensionsG*:kappaMG = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "0.0")
{
$data = "ExtraDimensionsG*:Gll = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "0.0")
{
$data = "ExtraDimensionsG*:Gqq = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0.0")
{
$data = "ExtraDimensionsG*:Gbb = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "0.001")
{
$data = "ExtraDimensionsG*:Gtt = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "0.000013")
{
$data = "ExtraDimensionsG*:Ggg = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "0.000013")
{
$data = "ExtraDimensionsG*:Ggmgm = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "0.001")
{
$data = "ExtraDimensionsG*:GZZ = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "0.001")
{
$data = "ExtraDimensionsG*:GWW = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "0.001")
{
$data = "ExtraDimensionsG*:Ghh = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "-0.2")
{
$data = "ExtraDimensionsG*:KKgqR = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "-0.2")
{
$data = "ExtraDimensionsG*:KKgqL = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "-0.2")
{
$data = "ExtraDimensionsG*:KKgbR = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "1.0")
{
$data = "ExtraDimensionsG*:KKgbL = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "4.0")
{
$data = "ExtraDimensionsG*:KKgtR = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "1.0")
{
$data = "ExtraDimensionsG*:KKgtL = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "0")
{
$data = "ExtraDimensionsG*:KKintMode = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2ddbar = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2uubar = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2ssbar = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2ccbar = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2bbbar = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2ttbar = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2e+e- = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2nuenuebar = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2mu+mu- = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2numunumubar = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2tau+tau- = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "off")
{
$data = "ExtraDimensionsTEV:ffbar2nutaunutaubar = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "3")
{
$data = "ExtraDimensionsTEV:gmZmode = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
if($_POST["40"] != "10")
{
$data = "ExtraDimensionsTEV:nMax = ".$_POST["40"]."\n";
fwrite($handle,$data);
}
if($_POST["41"] != "4000.0")
{
$data = "ExtraDimensionsTEV:mStar = ".$_POST["41"]."\n";
fwrite($handle,$data);
}
if($_POST["42"] != "off")
{
$data = "ExtraDimensionsLED:monojet = ".$_POST["42"]."\n";
fwrite($handle,$data);
}
if($_POST["43"] != "off")
{
$data = "ExtraDimensionsLED:gg2Gg = ".$_POST["43"]."\n";
fwrite($handle,$data);
}
if($_POST["44"] != "off")
{
$data = "ExtraDimensionsLED:qg2Gq = ".$_POST["44"]."\n";
fwrite($handle,$data);
}
if($_POST["45"] != "off")
{
$data = "ExtraDimensionsLED:qqbar2Gg = ".$_POST["45"]."\n";
fwrite($handle,$data);
}
if($_POST["46"] != "off")
{
$data = "ExtraDimensionsLED:ffbar2GZ = ".$_POST["46"]."\n";
fwrite($handle,$data);
}
if($_POST["47"] != "off")
{
$data = "ExtraDimensionsLED:ffbar2Ggamma = ".$_POST["47"]."\n";
fwrite($handle,$data);
}
if($_POST["48"] != "off")
{
$data = "ExtraDimensionsLED:ffbar2gammagamma = ".$_POST["48"]."\n";
fwrite($handle,$data);
}
if($_POST["49"] != "off")
{
$data = "ExtraDimensionsLED:gg2gammagamma = ".$_POST["49"]."\n";
fwrite($handle,$data);
}
if($_POST["50"] != "off")
{
$data = "ExtraDimensionsLED:ffbar2llbar = ".$_POST["50"]."\n";
fwrite($handle,$data);
}
if($_POST["51"] != "off")
{
$data = "ExtraDimensionsLED:gg2llbar = ".$_POST["51"]."\n";
fwrite($handle,$data);
}
if($_POST["52"] != "off")
{
$data = "ExtraDimensionsLED:dijets = ".$_POST["52"]."\n";
fwrite($handle,$data);
}
if($_POST["53"] != "off")
{
$data = "ExtraDimensionsLED:gg2DJgg = ".$_POST["53"]."\n";
fwrite($handle,$data);
}
if($_POST["54"] != "off")
{
$data = "ExtraDimensionsLED:gg2DJqqbar = ".$_POST["54"]."\n";
fwrite($handle,$data);
}
if($_POST["55"] != "off")
{
$data = "ExtraDimensionsLED:qg2DJqg = ".$_POST["55"]."\n";
fwrite($handle,$data);
}
if($_POST["56"] != "off")
{
$data = "ExtraDimensionsLED:qq2DJqq = ".$_POST["56"]."\n";
fwrite($handle,$data);
}
if($_POST["57"] != "off")
{
$data = "ExtraDimensionsLED:qqbar2DJgg = ".$_POST["57"]."\n";
fwrite($handle,$data);
}
if($_POST["58"] != "off")
{
$data = "ExtraDimensionsLED:qqbar2DJqqbarNew = ".$_POST["58"]."\n";
fwrite($handle,$data);
}
if($_POST["59"] != "off")
{
$data = "ExtraDimensionsLED:GravScalar = ".$_POST["59"]."\n";
fwrite($handle,$data);
}
if($_POST["60"] != "2")
{
$data = "ExtraDimensionsLED:n = ".$_POST["60"]."\n";
fwrite($handle,$data);
}
if($_POST["61"] != "2000.")
{
$data = "ExtraDimensionsLED:MD = ".$_POST["61"]."\n";
fwrite($handle,$data);
}
if($_POST["62"] != "2000.")
{
$data = "ExtraDimensionsLED:LambdaT = ".$_POST["62"]."\n";
fwrite($handle,$data);
}
if($_POST["63"] != "0")
{
$data = "ExtraDimensionsLED:NegInt = ".$_POST["63"]."\n";
fwrite($handle,$data);
}
if($_POST["64"] != "0")
{
$data = "ExtraDimensionsLED:CutOffMode = ".$_POST["64"]."\n";
fwrite($handle,$data);
}
if($_POST["65"] != "1.")
{
$data = "ExtraDimensionsLED:t = ".$_POST["65"]."\n";
fwrite($handle,$data);
}
if($_POST["66"] != "1.0")
{
$data = "ExtraDimensionsLED:g = ".$_POST["66"]."\n";
fwrite($handle,$data);
}
if($_POST["67"] != "1.0")
{
$data = "ExtraDimensionsLED:c = ".$_POST["67"]."\n";
fwrite($handle,$data);
}
if($_POST["68"] != "3")
{
$data = "ExtraDimensionsLED:nQuarkNew = ".$_POST["68"]."\n";
fwrite($handle,$data);
}
if($_POST["69"] != "0")
{
$data = "ExtraDimensionsLED:opMode = ".$_POST["69"]."\n";
fwrite($handle,$data);
}
if($_POST["70"] != "off")
{
$data = "ExtraDimensionsUnpart:monojet = ".$_POST["70"]."\n";
fwrite($handle,$data);
}
if($_POST["71"] != "off")
{
$data = "ExtraDimensionsUnpart:gg2Ug = ".$_POST["71"]."\n";
fwrite($handle,$data);
}
if($_POST["72"] != "off")
{
$data = "ExtraDimensionsUnpart:qg2Uq = ".$_POST["72"]."\n";
fwrite($handle,$data);
}
if($_POST["73"] != "off")
{
$data = "ExtraDimensionsUnpart:qqbar2Ug = ".$_POST["73"]."\n";
fwrite($handle,$data);
}
if($_POST["74"] != "off")
{
$data = "ExtraDimensionsUnpart:ffbar2UZ = ".$_POST["74"]."\n";
fwrite($handle,$data);
}
if($_POST["75"] != "off")
{
$data = "ExtraDimensionsUnpart:ffbar2Ugamma = ".$_POST["75"]."\n";
fwrite($handle,$data);
}
if($_POST["76"] != "off")
{
$data = "ExtraDimensionsUnpart:ffbar2gammagamma = ".$_POST["76"]."\n";
fwrite($handle,$data);
}
if($_POST["77"] != "off")
{
$data = "ExtraDimensionsUnpart:gg2gammagamma = ".$_POST["77"]."\n";
fwrite($handle,$data);
}
if($_POST["78"] != "off")
{
$data = "ExtraDimensionsUnpart:ffbar2llbar = ".$_POST["78"]."\n";
fwrite($handle,$data);
}
if($_POST["79"] != "off")
{
$data = "ExtraDimensionsUnpart:gg2llbar = ".$_POST["79"]."\n";
fwrite($handle,$data);
}
if($_POST["80"] != "2")
{
$data = "ExtraDimensionsUnpart:spinU = ".$_POST["80"]."\n";
fwrite($handle,$data);
}
if($_POST["81"] != "1.4")
{
$data = "ExtraDimensionsUnpart:dU = ".$_POST["81"]."\n";
fwrite($handle,$data);
}
if($_POST["82"] != "2000.")
{
$data = "ExtraDimensionsUnpart:LambdaU = ".$_POST["82"]."\n";
fwrite($handle,$data);
}
if($_POST["83"] != "1.0")
{
$data = "ExtraDimensionsUnpart:lambda = ".$_POST["83"]."\n";
fwrite($handle,$data);
}
if($_POST["84"] != "1.0")
{
$data = "ExtraDimensionsUnpart:ratio = ".$_POST["84"]."\n";
fwrite($handle,$data);
}
if($_POST["85"] != "0")
{
$data = "ExtraDimensionsUnpart:CutOffMode = ".$_POST["85"]."\n";
fwrite($handle,$data);
}
if($_POST["86"] != "0")
{
$data = "ExtraDimensionsUnpart:gXX = ".$_POST["86"]."\n";
fwrite($handle,$data);
}
if($_POST["87"] != "0")
{
$data = "ExtraDimensionsUnpart:gXY = ".$_POST["87"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->

