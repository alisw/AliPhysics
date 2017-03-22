<html>
<head>
<title>New-Gauge-Boson Processes</title>
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

<form method='post' action='NewGaugeBosonProcesses.php'>
 
<h2>New-Gauge-Boson Processes</h2> 
 
This page contains the production of new <i>Z'^0</i> and 
<i>W'^+-</i> gauge bosons, e.g. within the context of a new 
<i>U(1)</i> or <i>SU(2)</i> gauge group, and also a 
(rather speculative) horizontal gauge boson <i>R^0</i>. 
Left-right-symmetry scenarios also contain new gauge bosons, 
but are described 
<?php $filepath = $_GET["filepath"];
echo "<a href='LeftRightSymmetryProcesses.php?filepath=".$filepath."' target='page'>";?>separately</a>. 
 
<h3><i>Z'^0</i></h3> 
 
This group only contains one subprocess, with the full 
<i>gamma^*/Z^0/Z'^0</i> interference structure for couplings 
to fermion pairs. It is possible to pick only a subset, e.g, only 
the pure <i>Z'^0</i> piece. No higher-order processes are 
available explicitly, but the ISR showers contain automatic 
matching to the <i>Z'^0</i> + 1 jet matrix elements, as for 
the corresponding &nbsp <i>gamma^*/Z^0</i> process. 
 
<br/><br/><strong>NewGaugeBoson:ffbar2gmZZprime</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar &rarr;Z'^0</i>. 
Code 3001. 
   
 
<br/><br/><table><tr><td><strong>Zprime:gmZmode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 6</code>)</td></tr></table>
Choice of full <ei>gamma^*/Z^0/Z'^0</ei> structure or not in 
the above process. Note that, with the <ei>Z'^0</ei> part switched 
off, this process is reduced to what already exists among 
<aloc href="ElectroweakProcesses">electroweak processes</aloc>, 
so those options are here only for crosschecks. 
<br/>
<input type="radio" name="2" value="0" checked="checked"><strong>0 </strong>: full <ei>gamma^*/Z^0/Z'^0</ei> structure,  with interference included.<br/>
<input type="radio" name="2" value="1"><strong>1 </strong>: only pure <ei>gamma^*</ei> contribution.<br/>
<input type="radio" name="2" value="2"><strong>2 </strong>: only pure <ei>Z^0</ei> contribution.<br/>
<input type="radio" name="2" value="3"><strong>3 </strong>: only pure <ei>Z'^0</ei> contribution.<br/>
<input type="radio" name="2" value="4"><strong>4 </strong>: only the <ei>gamma^*/Z^0</ei> contribution,  including interference.<br/>
<input type="radio" name="2" value="5"><strong>5 </strong>: only the <ei>gamma^*/Z'^0</ei> contribution,  including interference.<br/>
<input type="radio" name="2" value="6"><strong>6 </strong>: only the <ei>Z^0/Z'^0</ei> contribution,  including interference.<br/>
<br/><b>Note</b>: irrespective of the option used, the particle produced 
will always be assigned code 32 for <ei>Z'^0</ei>, and open decay channels 
is purely dictated by what is set for the <ei>Z'^0</ei>. 
 
<p/> 
The couplings of the <i>Z'^0</i> to quarks and leptons can 
either be assumed universal, i.e. generation-independent, or not. 
In the former case eight numbers parametrize the vector and axial 
couplings of down-type quarks, up-type quarks, leptons and neutrinos, 
respectively. Depending on your assumed neutrino nature you may 
want to restrict your freedom in that sector, but no limitations 
are enforced by the program. The default corresponds to the same 
couplings as that of the Standard Model <i>Z^0</i>, with axial 
couplings <i>a_f = +-1</i> and vector couplings 
<i>v_f = a_f - 4 e_f sin^2(theta_W)</i>, with 
<i>sin^2(theta_W) = 0.23</i>. Without universality 
the same eight numbers have to be set separately also for the 
second and the third generation. The choice of fixed axial and 
vector couplings implies a resonance width that increases linearly 
with the <i>Z'^0</i> mass. 
 
<p/> 
By a suitable choice of the parameters, it is possible to simulate 
just about any imaginable <i>Z'^0</i> scenario, with full 
interference effects in cross sections and decay angular 
distributions and generation-dependent couplings; the default values 
should mainly be viewed as placeholders. The conversion 
from the coupling conventions in a set of different <i>Z'^0</i> 
models in the literature to those used in PYTHIA is described in 
[<a href="Bibliography.php" target="page">Cio08</a>]. 
 
<br/><br/><strong>Zprime:universality</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If on then you need only set the first-generation couplings 
below, and these are automatically also used for the second and 
third generation. If off, then couplings can be chosen separately 
for each generation. 
   
 
<p/> 
Here are the couplings always valid for the first generation, 
and normally also for the second and third by trivial analogy: 
 
<br/><br/><table><tr><td><strong>Zprime:vd </td><td></td><td> <input type="text" name="4" value="-0.693" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-0.693</strong></code>)</td></tr></table>
vector coupling of <i>d</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:ad </td><td></td><td> <input type="text" name="5" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
axial coupling of <i>d</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vu </td><td></td><td> <input type="text" name="6" value="0.387" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.387</strong></code>)</td></tr></table>
vector coupling of <i>u</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:au </td><td></td><td> <input type="text" name="7" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
axial coupling of <i>u</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:ve </td><td></td><td> <input type="text" name="8" value="-0.08" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-0.08</strong></code>)</td></tr></table>
vector coupling of <i>e</i> leptons. 
   
 
<br/><br/><table><tr><td><strong>Zprime:ae </td><td></td><td> <input type="text" name="9" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
axial coupling of <i>e</i> leptons. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vnue </td><td></td><td> <input type="text" name="10" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
vector coupling of <i>nu_e</i> neutrinos. 
   
 
<br/><br/><table><tr><td><strong>Zprime:anue </td><td></td><td> <input type="text" name="11" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
axial coupling of <i>nu_e</i> neutrinos. 
   
 
<p/> 
Here are the further couplings that are specific for 
a scenario with <code>Zprime:universality</code> switched off: 
 
<br/><br/><table><tr><td><strong>Zprime:vs </td><td></td><td> <input type="text" name="12" value="-0.693" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-0.693</strong></code>)</td></tr></table>
vector coupling of <i>s</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:as </td><td></td><td> <input type="text" name="13" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
axial coupling of <i>s</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vc </td><td></td><td> <input type="text" name="14" value="0.387" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.387</strong></code>)</td></tr></table>
vector coupling of <i>c</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:ac </td><td></td><td> <input type="text" name="15" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
axial coupling of <i>c</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vmu </td><td></td><td> <input type="text" name="16" value="-0.08" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-0.08</strong></code>)</td></tr></table>
vector coupling of <i>mu</i> leptons. 
   
 
<br/><br/><table><tr><td><strong>Zprime:amu </td><td></td><td> <input type="text" name="17" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
axial coupling of <i>mu</i> leptons. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vnumu </td><td></td><td> <input type="text" name="18" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
vector coupling of <i>nu_mu</i> neutrinos. 
   
 
<br/><br/><table><tr><td><strong>Zprime:anumu </td><td></td><td> <input type="text" name="19" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
axial coupling of <i>nu_mu</i> neutrinos. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vb </td><td></td><td> <input type="text" name="20" value="-0.693" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-0.693</strong></code>)</td></tr></table>
vector coupling of <i>b</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:ab </td><td></td><td> <input type="text" name="21" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
axial coupling of <i>b</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vt </td><td></td><td> <input type="text" name="22" value="0.387" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.387</strong></code>)</td></tr></table>
vector coupling of <i>t</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:at </td><td></td><td> <input type="text" name="23" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
axial coupling of <i>t</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vtau </td><td></td><td> <input type="text" name="24" value="-0.08" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-0.08</strong></code>)</td></tr></table>
vector coupling of <i>tau</i> leptons. 
   
 
<br/><br/><table><tr><td><strong>Zprime:atau </td><td></td><td> <input type="text" name="25" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
axial coupling of <i>tau</i> leptons. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vnutau </td><td></td><td> <input type="text" name="26" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
vector coupling of <i>nu_tau</i> neutrinos. 
   
 
<br/><br/><table><tr><td><strong>Zprime:anutau </td><td></td><td> <input type="text" name="27" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
axial coupling of <i>nu_tau</i> neutrinos. 
   
 
<p/> 
The coupling to the decay channel <i>Z'^0 &rarr; W^+ W^-</i> is 
more model-dependent. By default it is therefore off, but can be 
switched on as follows. 
 
<br/><br/><table><tr><td><strong>Zprime:coup2WW </td><td></td><td> <input type="text" name="28" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
the coupling <i>Z'^0 &rarr; W^+ W^-</i> is taken to be this number 
times <i>m_W^2 / m_Z'^2</i> times the <i>Z^0 &rarr; W^+ W^-</i> 
coupling. Thus a unit value corresponds to the 
<i>Z^0 &rarr; W^+ W^-</i> coupling, scaled down by a factor 
<i>m_W^2 / m_Z'^2</i>, and gives a <i>Z'^0</i> partial 
width into this channel that again increases linearly. If you 
cancel this behaviour, by letting <code>Zprime:coup2WW</code> be 
proportional to <i>m_Z'^2 / m_W^2</i>, you instead obtain a 
partial width that goes like the fifth power of the <i>Z'^0</i> 
mass. These two extremes correspond to the "extended gauge model" 
and the "reference model", respectively, of [<a href="Bibliography.php" target="page">Alt89</a>]. 
Note that this channel only includes the pure <i>Z'</i> part, 
while <i>f fbar &rarr; gamma^*/Z^*0 &rarr; W^+ W^-</i> is available 
as a separate electroweak process. 
   
 
Furthermore, we have left some amount of 
freedom in the choice of decay angular correlations in this 
channel, but obviously alternative shapes could be imagined. 
 
<br/><br/><table><tr><td><strong>Zprime:anglesWW </td><td></td><td> <input type="text" name="29" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
in the decay chain <i>Z'^0 &rarr; W^+ W^- &rarr;f_1 fbar_2 f_3 fbar_4</i> 
the decay angular distributions is taken to be a mixture of two 
possible shapes. This parameter gives the fraction that is distributed 
as in Higgs <i>h^0 &rarr; W^+ W^-</i> (longitudinal bosons), 
with the remainder (by default all) is taken to be the same as for 
<i>Z^0 &rarr; W^+ W^-</i> (a mixture of transverse and longitudinal 
bosons). 
   
 
<p/> 
A massive <i>Z'^0</i> is also likely to decay into Higgs bosons 
and potentially into other now unknown particles. Such possibilities 
clearly are quite model-dependent, and have not been included 
for now. 
 
<p/> 
Finally, to allow the exploration of more BSM physics scenarios, 
we include the possibility of the <i>Z'^0</i> (and hence the 
<i>gamma</i> and <i>Z^0</i>) coupling to a fourth generation of fermions. 
This provides redundancy with and extensions beyond those processes 
implemented as 
<?php $filepath = $_GET["filepath"];
echo "<a href='Fourth-Generation Processes.php?filepath=".$filepath."' target='page'>";?>fourth-generation processes</a>. 
By default, the decay channels for the fourth-generation and not included. 
They are enabled using: 
<br/><br/><strong>Zprime:coup2gen4</strong>  <input type="radio" name="30" value="on"><strong>On</strong>
<input type="radio" name="30" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<i> Z'^0 </i> couples to 4th generation fermions. 
   
 
<p/> 
Here are the further couplings that are specific for 
a scenario with <code>Zprime:universality</code> switched off: 
 
<br/><br/><table><tr><td><strong>Zprime:vbPrime </td><td></td><td> <input type="text" name="31" value="-0.693" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-0.693</strong></code>)</td></tr></table>
vector coupling of <i>b'</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:abPrime </td><td></td><td> <input type="text" name="32" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
axial coupling of <i>b'</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vtPrime </td><td></td><td> <input type="text" name="33" value="0.387" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.387</strong></code>)</td></tr></table>
vector coupling of <i>t'</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:atPrime </td><td></td><td> <input type="text" name="34" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
axial coupling of <i>t'</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vtauPrime </td><td></td><td> <input type="text" name="35" value="-0.08" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-0.08</strong></code>)</td></tr></table>
vector coupling of <i>tau'</i> leptons. 
   
 
<br/><br/><table><tr><td><strong>Zprime:atauPrime </td><td></td><td> <input type="text" name="36" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
axial coupling of <i>tau'</i> leptons. 
   
 
<br/><br/><table><tr><td><strong>Zprime:vnutauPrime </td><td></td><td> <input type="text" name="37" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
vector coupling of <i>nu_tau'</i> neutrinos. 
   
 
<br/><br/><table><tr><td><strong>Zprime:anutauPrime </td><td></td><td> <input type="text" name="38" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
axial coupling of <i>nu_tau'</i> neutrinos. 
   
 
 
<h3><i>W'^+-</i></h3> 
 
The <i>W'^+-</i> implementation is less ambitious than the 
<i>Z'^0</i>. Specifically, while indirect detection of a 
<i>Z'^0</i> through its interference contribution is 
a possible discovery channel in lepton colliders, there is no 
equally compelling case for <i>W^+-/W'^+-</i> interference 
effects being of importance for discovery, and such interference 
has therefore not been implemented for now. Related to this, a 
<i>Z'^0</i> could appear on its own in a new <i>U(1)</i> group, 
while <i>W'^+-</i> would have to sit in a <i>SU(2)</i> group 
and thus have a <i>Z'^0</i> partner that is likely to be found 
first. Only one process is implemented but, like for the 
<i>W^+-</i>, the ISR showers contain automatic matching to the 
<i>W'^+-</i> + 1 jet matrix elements. 
 
<br/><br/><strong>NewGaugeBoson:ffbar2Wprime</strong>  <input type="radio" name="39" value="on"><strong>On</strong>
<input type="radio" name="39" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar' &rarr; W'^+-</i>. 
Code 3021. 
   
 
<p/> 
The couplings of the <i>W'^+-</i> are here assumed universal, 
i.e. the same for all generations. One may set vector and axial 
couplings freely, separately for the <i>q qbar'</i> and the 
<i>l nu_l</i> decay channels. The defaults correspond to the 
<i>V - A</i> structure and normalization of the Standard Model 
<i>W^+-</i>, but can be changed to simulate a wide selection 
of models. One limitation is that, for simplicity, the same 
Cabibbo--Kobayashi--Maskawa quark mixing matrix is assumed as for 
the standard <i>W^+-</i>. Depending on your assumed neutrino 
nature you may want to restrict your freedom in the lepton sector, 
but no limitations are enforced by the program. 
 
<br/><br/><table><tr><td><strong>Wprime:vq </td><td></td><td> <input type="text" name="40" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
vector coupling of quarks. 
   
 
<br/><br/><table><tr><td><strong>Wprime:aq </td><td></td><td> <input type="text" name="41" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
axial coupling of quarks. 
   
 
<br/><br/><table><tr><td><strong>Wprime:vl </td><td></td><td> <input type="text" name="42" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
vector coupling of leptons. 
   
 
<br/><br/><table><tr><td><strong>Wprime:al </td><td></td><td> <input type="text" name="43" value="-1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.</strong></code>)</td></tr></table>
axial coupling of leptons. 
   
 
<p/> 
The coupling to the decay channel <i>W'^+- &rarr; W^+- Z^0</i> is 
more model-dependent, like for <i>Z'^0 &rarr; W^+ W^-</i> described 
above. By default it is therefore off, but can be 
switched on as follows. Furthermore, we have left some amount of 
freedom in the choice of decay angular correlations in this 
channel, but obviously alternative shapes could be imagined. 
 
<br/><br/><table><tr><td><strong>Wprime:coup2WZ </td><td></td><td> <input type="text" name="44" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
the coupling <i>W'^0 &rarr; W^+- Z^0</i> is taken to be this number 
times <i>m_W^2 / m_W'^2</i> times the <i>W^+- &rarr; W^+- Z^0</i> 
coupling. Thus a unit value corresponds to the 
<i>W^+- &rarr; W^+- Z^0</i> coupling, scaled down by a factor 
<i>m_W^2 / m_W'^2</i>, and gives a <i>W'^+-</i> partial 
width into this channel that increases linearly with the 
<i>W'^+-</i> mass. If you cancel this behaviour, by letting 
<code>Wprime:coup2WZ</code> be proportional to <i>m_W'^2 / m_W^2</i>, 
you instead obtain a partial width that goes like the fifth power 
of the <i>W'^+-</i> mass. These two extremes correspond to the 
"extended gauge model" and the "reference model", respectively, 
of [<a href="Bibliography.php" target="page">Alt89</a>]. 
   
 
<br/><br/><table><tr><td><strong>Wprime:anglesWZ </td><td></td><td> <input type="text" name="45" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
in the decay chain <i>W'^+- &rarr; W^+- Z^0 &rarr;f_1 fbar_2 f_3 fbar_4</i> 
the decay angular distributions is taken to be a mixture of two 
possible shapes. This parameter gives the fraction that is distributed 
as in Higgs <i>H^+- &rarr; W^+- Z^0</i> (longitudinal bosons), 
with the remainder (by default all) is taken to be the same as for 
<i>W^+- &rarr; W^+- Z^0</i> (a mixture of transverse and longitudinal 
bosons). 
   
 
<p/> 
A massive <i>W'^+-</i> is also likely to decay into Higgs bosons 
and potentially into other now unknown particles. Such possibilities 
clearly are quite model-dependent, and have not been included 
for now. 
 
<h3><i>R^0</i></h3> 
 
The <i>R^0</i> boson (particle code 41) represents one possible 
scenario for a horizontal gauge boson, i.e. a gauge boson 
that couples between the generations, inducing processes like 
<i>s dbar &rarr; R^0 &rarr; mu^- e^+</i>. Experimental limits on 
flavour-changing neutral currents forces such a boson to be fairly 
heavy. In spite of being neutral the antiparticle is distinct from 
the particle: one carries a net positive generation number and 
the other a negative one. This particular model has no new 
parameters beyond the <i>R^0</i> mass. Decays are assumed isotropic. 
For further details see [<a href="Bibliography.php" target="page">Ben85</a>]. 
 
<br/><br/><strong>NewGaugeBoson:ffbar2R0</strong>  <input type="radio" name="46" value="on"><strong>On</strong>
<input type="radio" name="46" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f_1 fbar_2 &rarr; R^0 &rarr; f_3 fbar_4</i>, where 
<i>f_1</i> and <i>fbar_2</i> are separated by <i>+-</i> one 
generation and similarly for <i>f_3</i> and <i>fbar_4</i>. 
Thus possible final states are e.g. <i>d sbar</i>, <i>u cbar</i> 
<i>s bbar</i>, <i>c tbar</i>, <i>e- mu+</i> and 
<i>mu- tau+</i>. 
Code 3041. 
   
 
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
$data = "NewGaugeBoson:ffbar2gmZZprime = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0")
{
$data = "Zprime:gmZmode = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "Zprime:universality = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "-0.693")
{
$data = "Zprime:vd = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "-1.")
{
$data = "Zprime:ad = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.387")
{
$data = "Zprime:vu = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "1.")
{
$data = "Zprime:au = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "-0.08")
{
$data = "Zprime:ve = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "-1.")
{
$data = "Zprime:ae = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "1.")
{
$data = "Zprime:vnue = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "1.")
{
$data = "Zprime:anue = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "-0.693")
{
$data = "Zprime:vs = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "-1.")
{
$data = "Zprime:as = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "0.387")
{
$data = "Zprime:vc = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "1.")
{
$data = "Zprime:ac = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "-0.08")
{
$data = "Zprime:vmu = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "-1.")
{
$data = "Zprime:amu = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "1.")
{
$data = "Zprime:vnumu = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "1.")
{
$data = "Zprime:anumu = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "-0.693")
{
$data = "Zprime:vb = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "-1.")
{
$data = "Zprime:ab = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "0.387")
{
$data = "Zprime:vt = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "1.")
{
$data = "Zprime:at = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "-0.08")
{
$data = "Zprime:vtau = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "-1.")
{
$data = "Zprime:atau = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "1.")
{
$data = "Zprime:vnutau = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "1.")
{
$data = "Zprime:anutau = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "0.")
{
$data = "Zprime:coup2WW = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "0.")
{
$data = "Zprime:anglesWW = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "off")
{
$data = "Zprime:coup2gen4 = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "-0.693")
{
$data = "Zprime:vbPrime = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "-1.")
{
$data = "Zprime:abPrime = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "0.387")
{
$data = "Zprime:vtPrime = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "1.")
{
$data = "Zprime:atPrime = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "-0.08")
{
$data = "Zprime:vtauPrime = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "-1.")
{
$data = "Zprime:atauPrime = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "1.")
{
$data = "Zprime:vnutauPrime = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "1.")
{
$data = "Zprime:anutauPrime = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "off")
{
$data = "NewGaugeBoson:ffbar2Wprime = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
if($_POST["40"] != "1.")
{
$data = "Wprime:vq = ".$_POST["40"]."\n";
fwrite($handle,$data);
}
if($_POST["41"] != "-1.")
{
$data = "Wprime:aq = ".$_POST["41"]."\n";
fwrite($handle,$data);
}
if($_POST["42"] != "1.")
{
$data = "Wprime:vl = ".$_POST["42"]."\n";
fwrite($handle,$data);
}
if($_POST["43"] != "-1.")
{
$data = "Wprime:al = ".$_POST["43"]."\n";
fwrite($handle,$data);
}
if($_POST["44"] != "0.")
{
$data = "Wprime:coup2WZ = ".$_POST["44"]."\n";
fwrite($handle,$data);
}
if($_POST["45"] != "0.")
{
$data = "Wprime:anglesWZ = ".$_POST["45"]."\n";
fwrite($handle,$data);
}
if($_POST["46"] != "off")
{
$data = "NewGaugeBoson:ffbar2R0 = ".$_POST["46"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
