<html>
<head>
<title>Higgs Processes</title>
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

<form method='post' action='HiggsProcesses.php'>

<h2>Higgs Processes</h2>

This page documents Higgs production within and beyond the Standard Model
(SM and BSM for short). This includes several different processes and, 
for the BSM scenarios, a large set of parameters that would only be fixed 
within a more specific framework such as MSSM. Three choices can be made 
irrespective of the particular model:

<br/><br/><strong>Higgs:cubicWidth</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
The partial width of a Higgs particle to a pair of gauge bosons,
<i>W^+ W^-</i> or <i>Z^0 Z^0</i>, depends cubically on the
Higgs mass. When selecting the Higgs according to a Breit-Wigner,
so that the actual mass <i>mHat</i> does not agree with the
nominal <i>m_Higgs</i> one, an ambiguity arises which of the 
two to use [<a href="Bibliography.php" target="page">Sey95</a>]. The default is to use a linear 
dependence on <i>mHat</i>, i.e. a width proportional to 
<i>m_Higgs^2 * mHat</i>, while <code>on</code> gives a 
<i>mHat^3</i> dependence. This does not affect the widths to 
fermions, which only depend linearly on <i>mHat</i>.
This flag is used both for SM and BSM Higgs bosons.
  

<br/><br/><strong>Higgs:runningLoopMass</strong>  <input type="radio" name="2" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="2" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
The partial width of a Higgs particle to a pair of gluons or photons,
or a <i>gamma Z^0</i> pair, proceeds in part through quark loops,
mainly <i>b</i> and <i>t</i>. There is some ambiguity what kind
of masses to use. Default is running MSbar ones, but alternatively 
fixed pole masses are allowed (as was standard in PYTHIA 6), which 
typically gives a noticeably higher cross section for these channels.
(For a decay to a pair of fermions, such as top, the running mass is
used for couplings and the fixed one for phase space.)
  

<br/><br/><strong>Higgs:clipWings</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
The Breit-Wigner shape of a Higgs is nontrivial, owing to the rapid 
width variation with the mass of a Higgs. This implies that a Higgs 
of low nominal mass may still acquire a non-negligible high-end tail.
The validity of the calculation may be questioned in these wings. 
With this option on, the <code>Higgs:wingsFac</code> value is used to 
cut away the wings.
  

<br/><br/><table><tr><td><strong>Higgs:wingsFac </td><td></td><td> <input type="text" name="4" value="50." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>50.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
With <code>Higgs:clipWings</code> on, all Higgs masses which deviate 
from the nominal one by more than <code>Higgs:wingsFac</code>
times the nominal width are forbidden. This is achieved by setting
the <code>mMin</code> and <code>mMax</code> values of the Higgs states
at initialization (but never so as to allow a wider range than already
set by the user, alternatively by the default values).   
  

<h3>Standard-Model Higgs, basic processes</h3>

This section provides the standard set of processes that can be
run together to provide a reasonably complete overview of possible
production channels for a single SM Higgs. 
The main parameter is the choice of Higgs mass, which can be set in the
normal <code>ParticleData</code> database; thereafter the properties 
within the SM are essentially fixed. 

<br/><br/><strong>HiggsSM:all</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of Higgs production within the Standard Model.
  

<br/><br/><strong>HiggsSM:ffbar2H</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0</i>, where <i>f</i> sums over available
flavours except top. Related to the mass-dependent Higgs point coupling 
to fermions, so at hadron colliders the bottom contribution will
dominate.
Code 901.
  

<br/><br/><strong>HiggsSM:gg2H</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0</i> via loop contributions primarily from
top.
Code 902.
  

<br/><br/><strong>HiggsSM:gmgm2H</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>gamma gamma -> H^0</i> via loop contributions primarily 
from top and <i>W</i>.
Code 903.
  

<br/><br/><strong>HiggsSM:ffbar2HZ</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0 Z^0</i> via <i>s</i>-channel <i>Z^0</i>
exchange.
Code 904.
  

<br/><br/><strong>HiggsSM:ffbar2HW</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0 W^+-</i> via <i>s</i>-channel <i>W^+-</i>
exchange.
Code 905.
  

<br/><br/><strong>HiggsSM:ff2Hff(t:ZZ)</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f f' -> H^0 f f'</i> via <i>Z^0 Z^0</i> fusion.
Code 906.
  

<br/><br/><strong>HiggsSM:ff2Hff(t:WW)</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f_1 f_2 -> H^0 f_3 f_4</i> via <i>W^+ W^-</i> fusion.
Code 907.
  

<br/><br/><strong>HiggsSM:gg2Httbar</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0 t tbar</i> via <i>t tbar</i> fusion
(or, alternatively put, Higgs radiation off a top line).
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 908.
  

<br/><br/><strong>HiggsSM:qqbar2Httbar</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> H^0 t tbar</i> via <i>t tbar</i> fusion
(or, alternatively put, Higgs radiation off a top line).
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 909.
  

<h3>Standard-Model Higgs, further processes</h3>

A number of further production processes has been implemented, that 
are specializations of some of the above ones to the high-<i>pT</i> 
region. The sets therefore could not be used simultaneously
without unphysical double-counting, as further explained below. 
They are not switched on by the <code>HiggsSM:all</code> flag, but 
have to be switched on for each separate process after due consideration.

<p/>
The first three processes in this section are related to the Higgs
point coupling to fermions, and so primarily are of interest for 
<i>b</i> quarks. It is here useful to begin by reminding that 
a process like <i>b bbar -> H^0</i> implies that a <i>b/bbar</i> 
is taken from each incoming hadron, leaving behind its respective
antiparticle. The initial-state showers will then add one 
<i>g -> b bbar</i> branching on either side, so that effectively
the process becomes <i>g g -> H0 b bbar</i>. This would be the
same basic process as the <i>g g -> H^0 t tbar</i> one used for top.
The difference is that (a) no PDF's are defined for top and 
(b) the shower approach would not be good enough to provide sensible
kinematics for the <i>H^0 t tbar</i> subsystem. By contrast, owing 
to the <i>b</i> being much lighter than the Higgs, multiple 
gluon emissions must be resummed for <i>b</i>, as is done by PDF's 
and showers, in order to obtain a sensible description of the total 
production rate,  when the <i>b</i> quarks predominantly are produced 
at small <i>pT</i> values.

<br/><br/><strong>HiggsSM:qg2Hq</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> H^0 q</i>. This process gives first-order 
corrections to the <i>f fbar -> H^0</i> one above, and should only be 
used to study  the high-<i>pT</i> tail, while <i>f fbar -> H^0</i> 
should be used for inclusive production. Only the dominant <i>c</i> 
and <i>b</i> contributions are included, and generated separately 
for technical reasons. Note that another first-order process would be 
<i>q qbar -> H^0 g</i>, which is not explicitly implemented here,
but is obtained from showering off the lowest-order process. It does not 
contain any <i>b</i> at large <i>pT</i>, however, so is less 
interesting for many applications. 
Code 911.

  
<br/><br/><strong>HiggsSM:gg2Hbbbar</strong>  <input type="radio" name="16" value="on"><strong>On</strong>
<input type="radio" name="16" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0 b bbar</i>. This process is yet one order 
higher of the <i>b bbar -> H^0</i> and <i>b g -> H^0 b</i> chain,
where now two quarks should be required above some large <i>pT</i>
threshold.
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 912.
  

<br/><br/><strong>HiggsSM:qqbar2Hbbbar</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> H^0 b bbar</i> via an <i>s</i>-channel
gluon, so closely related to the previous one, but typically less 
important owing to the smaller rate of (anti)quarks relative to 
gluons.
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 913.
  

<p/>
The second set of processes are predominantly first-order corrections 
to the <i>g g -> H^0</i> process, again dominated by the top loop.
We here only provide the kinematical expressions obtained in the 
limit that the top quark goes to infinity, but scaled to the 
finite-top-mass coupling in <i>g g -> H^0</i>. (Complete loop
expressions are available e.g. in PYTHIA 6.4 but are very lengthy.) 
This provides a reasonably accurate description for "intermediate" 
<i>pT</i> values, but fails when the <i>pT</i> scale approaches
the top mass. 
 
<br/><br/><strong>HiggsSM:gg2Hg(l:t)</strong>  <input type="radio" name="18" value="on"><strong>On</strong>
<input type="radio" name="18" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0 g</i> via loop contributions primarily 
from top.
Code 914.
  
 
<br/><br/><strong>HiggsSM:qg2Hq(l:t)</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> H^0 q</i> via loop contributions primarily 
from top. Not to be confused with the <code>HiggsSM:qg2Hq</code>
process above, with its direct fermion-to-Higgs coupling.
Code 915.
  
 
<br/><br/><strong>HiggsSM:qqbar2Hg(l:t)</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> H^0 g</i> via an <i>s</i>-channel gluon
and loop contributions primarily from top. Is strictly speaking a 
"new" process, not directly derived from <i>g g -> H^0</i>, and
could therefore be included in the standard mix without double-counting, 
but is numerically negligible.
Code 916.
  

<h3>Beyond-the-Standard-Model Higgs, introduction</h3>

Further Higgs multiplets arise in a number of scenarios. We here 
concentrate on the MSSM scenario with two Higgs doublets, but with 
flexibility enough that also other two-Higgs-doublet scenarios could 
be represented by a suitable choice of parameters. Conventionally the 
Higgs states are labeled <i>h^0, H^0, A^0</i> and <i>H^+-</i>.
If the scalar and pseudocalar states mix the resulting states are 
labeled <i>H_1^0, H_2^0, H_3^0</i>. In process names and parameter 
explanations both notations will be used, but for settings labels 
we have adapted the shorthand hybrid notation <code>H1</code> for
<i>h^0(H_1^0)</i>, <code>H2</code> for <i>H^0(H_2^0)</i> and
<code>A3</code> for <i>A^0(H_3^0)</i>. (Recall that the 
<code>Settings</code> database does not distinguish upper- and lowercase 
characters, so that the user has one thing less to worry about, but here 
it causes problems with <i>h^0</i> vs. <i>H^0</i>.) We leave the issue 
of mass ordering between <i>H^0</i> and <i>A^0</i> open, and thereby 
also that of <i>H_2^0</i> and <i>H_3^0</i>.

<br/><br/><strong>Higgs:useBSM</strong>  <input type="radio" name="21" value="on"><strong>On</strong>
<input type="radio" name="21" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Master switch to initialize and use the two-Higgs-doublet states. 
If off, only the above SM Higgs processes can be used, with couplings 
as predicted in the SM. If on, only the below BSM Higgs processes can 
be used, with couplings that can be set freely, also found further down 
on this page. 
  

<h3>Beyond-the-Standard-Model Higgs, basic processes</h3>

This section provides the standard set of processes that can be
run together to provide a reasonably complete overview of possible
production channels for a single neutral Higgs state in a two-doublet
scenarios such as MSSM. The list of processes for neutral states closely 
mimics the one found for the SM Higgs. Some of the processes 
vanish for a pure pseudoscalar <i>A^0</i>, but are kept for flexibility 
in cases of mixing with the scalar <i>h^0</i> and <i>H^0</i> states, 
or for use in the context of non-MSSM models. This should work well to 
represent e.g. that a small admixture of the "wrong" parity would allow 
a process such as <i>q qbar -> A^0 Z^0</i>, which otherwise is forbidden. 
However, note that the loop integrals e.g. for <i>g g -> h^0/H^0/A^0</i>
are hardcoded to be for scalars for the former two particles and for a
pseudoscalar for the latter one, so absolute rates would not be 
correctly represented in the case of large scalar/pseudoscalar mixing.  

<br/><br/><strong>HiggsBSM:all</strong>  <input type="radio" name="22" value="on"><strong>On</strong>
<input type="radio" name="22" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of Higgs production beyond the Standard Model,
as listed below.
  

<h4>1) <i>h^0(H_1^0)</i> processes</h4>

<br/><br/><strong>HiggsBSM:allH1</strong>  <input type="radio" name="23" value="on"><strong>On</strong>
<input type="radio" name="23" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of <i>h^0(H_1^0)</i> production processes.
  

<br/><br/><strong>HiggsBSM:ffbar2H1</strong>  <input type="radio" name="24" value="on"><strong>On</strong>
<input type="radio" name="24" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> h^0(H_1^0)</i>, where <i>f</i> sums over available
flavours except top.
Code 1001.
  

<br/><br/><strong>HiggsBSM:gg2H1</strong>  <input type="radio" name="25" value="on"><strong>On</strong>
<input type="radio" name="25" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> h^0(H_1^0)</i> via loop contributions primarily from
top.
Code 1002.
  

<br/><br/><strong>HiggsBSM:gmgm2H1</strong>  <input type="radio" name="26" value="on"><strong>On</strong>
<input type="radio" name="26" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>gamma gamma -> h^0(H_1^0)</i> via loop contributions 
primarily from top and <i>W</i>.
Code 1003.
  

<br/><br/><strong>HiggsBSM:ffbar2H1Z</strong>  <input type="radio" name="27" value="on"><strong>On</strong>
<input type="radio" name="27" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> h^0(H_1^0) Z^0</i> via <i>s</i>-channel 
<i>Z^0</i> exchange.
Code 1004.
  

<br/><br/><strong>HiggsBSM:ffbar2H1W</strong>  <input type="radio" name="28" value="on"><strong>On</strong>
<input type="radio" name="28" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> h^0(H_1^0) W^+-</i> via <i>s</i>-channel 
<i>W^+-</i> exchange.
Code 1005.
  

<br/><br/><strong>HiggsBSM:ff2H1ff(t:ZZ)</strong>  <input type="radio" name="29" value="on"><strong>On</strong>
<input type="radio" name="29" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f f' -> h^0(H_1^0) f f'</i> via <i>Z^0 Z^0</i> fusion.
Code 1006.
  

<br/><br/><strong>HiggsBSM:ff2H1ff(t:WW)</strong>  <input type="radio" name="30" value="on"><strong>On</strong>
<input type="radio" name="30" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f_1 f_2 -> h^0(H_1^0) f_3 f_4</i> via <i>W^+ W^-</i> 
fusion.
Code 1007.
  

<br/><br/><strong>HiggsBSM:gg2H1ttbar</strong>  <input type="radio" name="31" value="on"><strong>On</strong>
<input type="radio" name="31" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> h^0(H_1^0) t tbar</i> via <i>t tbar</i> fusion
(or, alternatively put, Higgs radiation off a top line).
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1008.
  

<br/><br/><strong>HiggsBSM:qqbar2H1ttbar</strong>  <input type="radio" name="32" value="on"><strong>On</strong>
<input type="radio" name="32" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> h^0(H_1^0) t tbar</i> via <i>t tbar</i> fusion
(or, alternatively put, Higgs radiation off a top line).
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1009.


<h4>2) <i>H^0(H_2^0)</i> processes</h4>

<br/><br/><strong>HiggsBSM:allH2</strong>  <input type="radio" name="33" value="on"><strong>On</strong>
<input type="radio" name="33" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of <i>H^0(H_2^0)</i> production processes.
  

<br/><br/><strong>HiggsBSM:ffbar2H2</strong>  <input type="radio" name="34" value="on"><strong>On</strong>
<input type="radio" name="34" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0(H_2^0)</i>, where <i>f</i> sums over available
flavours except top.
Code 1021.
  

<br/><br/><strong>HiggsBSM:gg2H2</strong>  <input type="radio" name="35" value="on"><strong>On</strong>
<input type="radio" name="35" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0(H_2^0)</i> via loop contributions primarily from
top.
Code 1022.
  

<br/><br/><strong>HiggsBSM:gmgm2H2</strong>  <input type="radio" name="36" value="on"><strong>On</strong>
<input type="radio" name="36" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>gamma gamma -> H^0(H_2^0)</i> via loop contributions primarily
from top and <i>W</i>.
Code 1023.
  

<br/><br/><strong>HiggsBSM:ffbar2H2Z</strong>  <input type="radio" name="37" value="on"><strong>On</strong>
<input type="radio" name="37" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0(H_2^0) Z^0</i> via <i>s</i>-channel 
<i>Z^0</i> exchange.
Code 1024.
  

<br/><br/><strong>HiggsBSM:ffbar2H2W</strong>  <input type="radio" name="38" value="on"><strong>On</strong>
<input type="radio" name="38" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0(H_2^0) W^+-</i> via <i>s</i>-channel 
<i>W^+-</i> exchange.
Code 1025.
  

<br/><br/><strong>HiggsBSM:ff2H2ff(t:ZZ)</strong>  <input type="radio" name="39" value="on"><strong>On</strong>
<input type="radio" name="39" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f f' -> H^0(H_2^0) f f'</i> via <i>Z^0 Z^0</i> fusion.
Code 1026.
  

<br/><br/><strong>HiggsBSM:ff2H2ff(t:WW)</strong>  <input type="radio" name="40" value="on"><strong>On</strong>
<input type="radio" name="40" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f_1 f_2 -> H^0(H_2^0) f_3 f_4</i> via <i>W^+ W^-</i> fusion.
Code 1027.
  

<br/><br/><strong>HiggsBSM:gg2H2ttbar</strong>  <input type="radio" name="41" value="on"><strong>On</strong>
<input type="radio" name="41" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0(H_2^0) t tbar</i> via <i>t tbar</i> fusion
(or, alternatively put, Higgs radiation off a top line).
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1028.
  

<br/><br/><strong>HiggsBSM:qqbar2H2ttbar</strong>  <input type="radio" name="42" value="on"><strong>On</strong>
<input type="radio" name="42" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> H^0(H_2^0) t tbar</i> via <i>t tbar</i> fusion
(or, alternatively put, Higgs radiation off a top line).
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1029.

<h4>3) <i>A^0(H_3^0)</i> processes</h4>

<br/><br/><strong>HiggsBSM:allA3</strong>  <input type="radio" name="43" value="on"><strong>On</strong>
<input type="radio" name="43" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of <i>A^0(H_3^0)</i> production processes.
  

<br/><br/><strong>HiggsBSM:ffbar2A3</strong>  <input type="radio" name="44" value="on"><strong>On</strong>
<input type="radio" name="44" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> A^0(H_3^0)</i>, where <i>f</i> sums over available
flavours except top.
Code 1041.
  

<br/><br/><strong>HiggsBSM:gg2A3</strong>  <input type="radio" name="45" value="on"><strong>On</strong>
<input type="radio" name="45" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> A^0(A_3^0)</i> via loop contributions primarily from
top.
Code 1042.
  

<br/><br/><strong>HiggsBSM:gmgm2A3</strong>  <input type="radio" name="46" value="on"><strong>On</strong>
<input type="radio" name="46" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>gamma gamma -> A^0(A_3^0)</i> via loop contributions primarily
from top and <i>W</i>.
Code 1043.
  

<br/><br/><strong>HiggsBSM:ffbar2A3Z</strong>  <input type="radio" name="47" value="on"><strong>On</strong>
<input type="radio" name="47" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> A^0(A_3^0) Z^0</i> via <i>s</i>-channel 
<i>Z^0</i> exchange.
Code 1044.
  

<br/><br/><strong>HiggsBSM:ffbar2A3W</strong>  <input type="radio" name="48" value="on"><strong>On</strong>
<input type="radio" name="48" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> A^0(A_3^0) W^+-</i> via <i>s</i>-channel 
<i>W^+-</i> exchange.
Code 1045.
  

<br/><br/><strong>HiggsBSM:ff2A3ff(t:ZZ)</strong>  <input type="radio" name="49" value="on"><strong>On</strong>
<input type="radio" name="49" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f f' -> A^0(A_3^0) f f'</i> via <i>Z^0 Z^0</i> fusion.
Code 1046.
  

<br/><br/><strong>HiggsBSM:ff2A3ff(t:WW)</strong>  <input type="radio" name="50" value="on"><strong>On</strong>
<input type="radio" name="50" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f_1 f_2 -> A^0(A_3^0) f_3 f_4</i> via <i>W^+ W^-</i> fusion.
Code 1047.
  

<br/><br/><strong>HiggsBSM:gg2A3ttbar</strong>  <input type="radio" name="51" value="on"><strong>On</strong>
<input type="radio" name="51" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> A^0(A_3^0) t tbar</i> via <i>t tbar</i> fusion
(or, alternatively put, Higgs radiation off a top line).
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1048.
  

<br/><br/><strong>HiggsBSM:qqbar2A3ttbar</strong>  <input type="radio" name="52" value="on"><strong>On</strong>
<input type="radio" name="52" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> A^0(A_3^0) t tbar</i> via <i>t tbar</i> fusion
(or, alternatively put, Higgs radiation off a top line).
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1049.

<h4>4) <i>H+-</i> processes</h4>

<br/><br/><strong>HiggsBSM:allH+-</strong>  <input type="radio" name="53" value="on"><strong>On</strong>
<input type="radio" name="53" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of <i>H^+-</i> production processes.
  

<br/><br/><strong>HiggsBSM:ffbar2H+-</strong>  <input type="radio" name="54" value="on"><strong>On</strong>
<input type="radio" name="54" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar' -> H^+-</i>, where <i>f, fbar'</i> sums over 
available incoming flavours. Since couplings are assumed 
generation-diagonal, in practice this means <i>c sbar -> H^+</i>
and <i>s cbar -> H^-</i>.
Code 1061.
  

<br/><br/><strong>HiggsBSM:bg2H+-t</strong>  <input type="radio" name="55" value="on"><strong>On</strong>
<input type="radio" name="55" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>b g -> H^+ tbar</i>. At hadron colliders this is the 
dominant process for single-charged-Higgs production.
Code 1062.
  

<h4>5) Higgs-pair processes</h4>

<br/><br/><strong>HiggsBSM:allHpair</strong>  <input type="radio" name="56" value="on"><strong>On</strong>
<input type="radio" name="56" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of Higgs pair-production processes.
  

<br/><br/><strong>HiggsBSM:ffbar2A3H1</strong>  <input type="radio" name="57" value="on"><strong>On</strong>
<input type="radio" name="57" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> A^0(H_3) h^0(H_1)</i>.
Code 1081.
  

<br/><br/><strong>HiggsBSM:ffbar2A3H2</strong>  <input type="radio" name="58" value="on"><strong>On</strong>
<input type="radio" name="58" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> A^0(H_3) H^0(H_2)</i>.
Code 1082.
  

<br/><br/><strong>HiggsBSM:ffbar2H+-H1</strong>  <input type="radio" name="59" value="on"><strong>On</strong>
<input type="radio" name="59" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^+- h^0(H_1)</i>.
Code 1083.
  

<br/><br/><strong>HiggsBSM:ffbar2H+-H2</strong>  <input type="radio" name="60" value="on"><strong>On</strong>
<input type="radio" name="60" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^+- H^0(H_2)</i>.
Code 1084.
  

<br/><br/><strong>HiggsBSM:ffbar2H+H-</strong>  <input type="radio" name="61" value="on"><strong>On</strong>
<input type="radio" name="61" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H+ H-</i>.
Code 1085.
  

<h3>Beyond-the-Standard-Model Higgs, further processes</h3>

This section mimics the above section on "Standard-Model Higgs, 
further processes", i.e. it contains higher-order corrections
to the processes already listed. The two sets therefore could not 
be used simultaneously without unphysical double-counting.
They are not controlled by any group flag, but have to be switched 
on for each separate process after due consideration. We refer to
the standard-model description for a set of further comments on
the processes.

<h4>1) <i>h^0(H_1^0)</i> processes</h4> 

<br/><br/><strong>HiggsBSM:qg2H1q</strong>  <input type="radio" name="62" value="on"><strong>On</strong>
<input type="radio" name="62" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> h^0 q</i>. This process gives first-order 
corrections to the <i>f fbar -> h^0</i> one above, and should only be 
used to study  the high-<i>pT</i> tail, while <i>f fbar -> h^0</i> 
should be used for inclusive production. Only the dominant <i>c</i> 
and <i>b</i> contributions are included, and generated separately 
for technical reasons. Note that another first-order process would be 
<i>q qbar -> h^0 g</i>, which is not explicitly implemented here,
but is obtained from showering off the lowest-order process. It does not 
contain any <i>b</i> at large <i>pT</i>, however, so is less 
interesting for many applications. 
Code 1011.
  

<br/><br/><strong>HiggsBSM:gg2H1bbbar</strong>  <input type="radio" name="63" value="on"><strong>On</strong>
<input type="radio" name="63" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> h^0 b bbar</i>. This process is yet one order
higher of the <i>b bbar -> h^0</i> and <i>b g -> h^0 b</i> chain,
where now two quarks should be required above some large <i>pT</i>
threshold.
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1012.
  

<br/><br/><strong>HiggsBSM:qqbar2H1bbbar</strong>  <input type="radio" name="64" value="on"><strong>On</strong>
<input type="radio" name="64" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> h^0 b bbar</i> via an <i>s</i>-channel
gluon, so closely related to the previous one, but typically less
important owing to the smaller rate of (anti)quarks relative to
gluons.
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1013.
  
 
<br/><br/><strong>HiggsBSM:gg2H1g(l:t)</strong>  <input type="radio" name="65" value="on"><strong>On</strong>
<input type="radio" name="65" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> h^0 g</i> via loop contributions primarily 
from top.
Code 1014.
  
 
<br/><br/><strong>HiggsBSM:qg2H1q(l:t)</strong>  <input type="radio" name="66" value="on"><strong>On</strong>
<input type="radio" name="66" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> h^0 q</i> via loop contributions primarily 
from top. Not to be confused with the <code>HiggsBSM:qg2H1q</code>
process above, with its direct fermion-to-Higgs coupling.
Code 1015.
  
 
<br/><br/><strong>HiggsBSM:qqbar2H1g(l:t)</strong>  <input type="radio" name="67" value="on"><strong>On</strong>
<input type="radio" name="67" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> h^0 g</i> via an <i>s</i>-channel gluon
and loop contributions primarily from top. Is strictly speaking a 
"new" process, not directly derived from <i>g g -> h^0</i>, and
could therefore be included in the standard mix without double-counting, 
but is numerically negligible.
Code 1016.
  

<h4>2) <i>H^0(H_2^0)</i> processes</h4>

<br/><br/><strong>HiggsBSM:qg2H2q</strong>  <input type="radio" name="68" value="on"><strong>On</strong>
<input type="radio" name="68" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> H^0 q</i>. This process gives first-order 
corrections to the <i>f fbar -> H^0</i> one above, and should only be 
used to study  the high-<i>pT</i> tail, while <i>f fbar -> H^0</i> 
should be used for inclusive production. Only the dominant <i>c</i> 
and <i>b</i> contributions are included, and generated separately 
for technical reasons. Note that another first-order process would be 
<i>q qbar -> H^0 g</i>, which is not explicitly implemented here,
but is obtained from showering off the lowest-order process. It does not 
contain any <i>b</i> at large <i>pT</i>, however, so is less 
interesting for many applications. 
Code 1031.
  

<br/><br/><strong>HiggsBSM:gg2H2bbbar</strong>  <input type="radio" name="69" value="on"><strong>On</strong>
<input type="radio" name="69" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0 b bbar</i>. This process is yet one order
higher of the <i>b bbar -> H^0</i> and <i>b g -> H^0 b</i> chain,
where now two quarks should be required above some large <i>pT</i>
threshold.
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1032.
  

<br/><br/><strong>HiggsBSM:qqbar2H2bbbar</strong>  <input type="radio" name="70" value="on"><strong>On</strong>
<input type="radio" name="70" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> H^0 b bbar</i> via an <i>s</i>-channel
gluon, so closely related to the previous one, but typically less
important owing to the smaller rate of (anti)quarks relative to
gluons.
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1033.
  
 
<br/><br/><strong>HiggsBSM:gg2H2g(l:t)</strong>  <input type="radio" name="71" value="on"><strong>On</strong>
<input type="radio" name="71" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0 g</i> via loop contributions primarily 
from top.
Code 1034.
  
 
<br/><br/><strong>HiggsBSM:qg2H2q(l:t)</strong>  <input type="radio" name="72" value="on"><strong>On</strong>
<input type="radio" name="72" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> H^0 q</i> via loop contributions primarily 
from top. Not to be confused with the <code>HiggsBSM:qg2H1q</code>
process above, with its direct fermion-to-Higgs coupling.
Code 1035.
  
 
<br/><br/><strong>HiggsBSM:qqbar2H2g(l:t)</strong>  <input type="radio" name="73" value="on"><strong>On</strong>
<input type="radio" name="73" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> H^0 g</i> via an <i>s</i>-channel gluon
and loop contributions primarily from top. Is strictly speaking a 
"new" process, not directly derived from <i>g g -> H^0</i>, and
could therefore be included in the standard mix without double-counting, 
but is numerically negligible.
Code 1036.
  

<h4>3) <i>A^0(H_3^0)</i> processes</h4>

<br/><br/><strong>HiggsBSM:qg2A3q</strong>  <input type="radio" name="74" value="on"><strong>On</strong>
<input type="radio" name="74" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> A^0 q</i>. This process gives first-order 
corrections to the <i>f fbar -> A^0</i> one above, and should only be 
used to study  the high-<i>pT</i> tail, while <i>f fbar -> A^0</i> 
should be used for inclusive production. Only the dominant <i>c</i> 
and <i>b</i> contributions are included, and generated separately 
for technical reasons. Note that another first-order process would be 
<i>q qbar -> A^0 g</i>, which is not explicitly implemented here,
but is obtained from showering off the lowest-order process. It does not 
contain any <i>b</i> at large <i>pT</i>, however, so is less 
interesting for many applications. 
Code 1051.
  

<br/><br/><strong>HiggsBSM:gg2A3bbbar</strong>  <input type="radio" name="75" value="on"><strong>On</strong>
<input type="radio" name="75" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> A^0 b bbar</i>. This process is yet one order
higher of the <i>b bbar -> A^0</i> and <i>b g -> A^0 b</i> chain,
where now two quarks should be required above some large <i>pT</i>
threshold.
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1052.
  

<br/><br/><strong>HiggsBSM:qqbar2A3bbbar</strong>  <input type="radio" name="76" value="on"><strong>On</strong>
<input type="radio" name="76" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> A^0 b bbar</i> via an <i>s</i>-channel
gluon, so closely related to the previous one, but typically less
important owing to the smaller rate of (anti)quarks relative to
gluons.
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 1053.
  
 
<br/><br/><strong>HiggsBSM:gg2A3g(l:t)</strong>  <input type="radio" name="77" value="on"><strong>On</strong>
<input type="radio" name="77" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> A^0 g</i> via loop contributions primarily 
from top.
Code 1054.
  
 
<br/><br/><strong>HiggsBSM:qg2A3q(l:t)</strong>  <input type="radio" name="78" value="on"><strong>On</strong>
<input type="radio" name="78" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> A^0 q</i> via loop contributions primarily 
from top. Not to be confused with the <code>HiggsBSM:qg2H1q</code>
process above, with its direct fermion-to-Higgs coupling.
Code 1055.
  
 
<br/><br/><strong>HiggsBSM:qqbar2A3g(l:t)</strong>  <input type="radio" name="79" value="on"><strong>On</strong>
<input type="radio" name="79" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> A^0 g</i> via an <i>s</i>-channel gluon
and loop contributions primarily from top. Is strictly speaking a 
"new" process, not directly derived from <i>g g -> A^0</i>, and
could therefore be included in the standard mix without double-counting, 
but is numerically negligible.
Code 1056.
  

<h3>Parameters for Beyond-the-Standard-Model Higgs production and decay</h3>

This section offers a big flexibility to set couplings of the various
Higgs states to fermions and gauge bosons, and also to each other.
The intention is that, for scenarios like MSSM, you should use standard 
input from the <?php $filepath = $_GET["filepath"];
echo "<a href='SUSYLesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>SUSY Les Houches 
Accord</a>, rather than having to set it all yourself. In other cases,
however, the freedom is there for you to use. Kindly note that some
of the internal calculations of partial widths from the parameters provided
do not include mixing between the scalar and pseudoscalar states. 

<p/>
Masses would be set in the <code>ParticleData</code> database,
while couplings are set below. When possible, the couplings of the Higgs 
states are normalized to the corresponding coupling within the SM. 
When not, their values within the MSSM are indicated, from which
it should be straightforward to understand what to use instead.
The exception is some couplings that vanish also in the MSSM, where the
normalization has been defined in close analogy with nonvanishing ones. 
Some parameter names are asymmetric but crossing can always be used,
i.e. the coupling for <i>A^0 -> H^0 Z^0</i> obviously is also valid
for <i>H^0 -> A^0 Z^0</i> and <i>Z^0 -> H^0 A^0</i>. 
Note that couplings usually appear quadratically in matrix elements.

<br/><br/><table><tr><td><strong>HiggsH1:coup2d </td><td></td><td> <input type="text" name="80" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>h^0(H_1^0)</i> coupling to down-type quarks.
  

<br/><br/><table><tr><td><strong>HiggsH1:coup2u </td><td></td><td> <input type="text" name="81" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>h^0(H_1^0)</i> coupling to up-type quarks.
  

<br/><br/><table><tr><td><strong>HiggsH1:coup2l </td><td></td><td> <input type="text" name="82" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>h^0(H_1^0)</i> coupling to (charged) leptons.
  

<br/><br/><table><tr><td><strong>HiggsH1:coup2Z </td><td></td><td> <input type="text" name="83" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>h^0(H_1^0)</i> coupling to <i>Z^0</i>.
  

<br/><br/><table><tr><td><strong>HiggsH1:coup2W </td><td></td><td> <input type="text" name="84" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>h^0(H_1^0)</i> coupling to <i>W^+-</i>.
  

<br/><br/><table><tr><td><strong>HiggsH1:coup2Hchg </td><td></td><td> <input type="text" name="85" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>h^0(H_1^0)</i> coupling to <i>H^+-</i> (in loops).
Is <i>sin(beta - alpha) + cos(2 beta) sin(beta + alpha) /
(2 cos^2theta_W)</i> in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsH2:coup2d </td><td></td><td> <input type="text" name="86" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>H^0(H_2^0)</i> coupling to down-type quarks.
  

<br/><br/><table><tr><td><strong>HiggsH2:coup2u </td><td></td><td> <input type="text" name="87" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>H^0(H_2^0)</i> coupling to up-type quarks.
  

<br/><br/><table><tr><td><strong>HiggsH2:coup2l </td><td></td><td> <input type="text" name="88" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>H^0(H_2^0)</i> coupling to (charged) leptons.
  

<br/><br/><table><tr><td><strong>HiggsH2:coup2Z </td><td></td><td> <input type="text" name="89" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>H^0(H_2^0)</i> coupling to <i>Z^0</i>.
  

<br/><br/><table><tr><td><strong>HiggsH2:coup2W </td><td></td><td> <input type="text" name="90" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>H^0(H_2^0)</i> coupling to <i>W^+-</i>.
  

<br/><br/><table><tr><td><strong>HiggsH2:coup2Hchg </td><td></td><td> <input type="text" name="91" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>H^0(H_2^0)</i> coupling to <i>H^+-</i> (in loops).
Is <i>cos(beta - alpha) + cos(2 beta) cos(beta + alpha) /
(2 cos^2theta_W)</i> in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsH2:coup2H1H1 </td><td></td><td> <input type="text" name="92" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>H^0(H_2^0)</i> coupling to a <i>h^0(H_1^0)</i> pair.
Is <i>cos(2 alpha) cos(beta + alpha) - 2 sin(2 alpha) 
sin(beta + alpha)</i> in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsH2:coup2A3A3 </td><td></td><td> <input type="text" name="93" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>H^0(H_2^0)</i> coupling to an <i>A^0(H_3^0)</i> pair.
Is <i>cos(2 beta) cos(beta + alpha)</i> in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsH2:coup2H1Z </td><td></td><td> <input type="text" name="94" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>H^0(H_2^0)</i> coupling to a <i>h^0(H_1^0) Z^0</i> pair.
Vanishes in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsH2:coup2A3H1 </td><td></td><td> <input type="text" name="95" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>H^0(H_2^0)</i> coupling to an <i>A^0(H_3^0) h^0(H_1^0)</i> pair.
Vanishes in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsH2:coup2HchgW </td><td></td><td> <input type="text" name="96" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>H^0(H_2^0)</i> coupling to a <i>H^+- W-+</i> pair.
Is <i>sin(beta - alpha)</i> in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsA3:coup2d </td><td></td><td> <input type="text" name="97" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>A^0(H_3^0)</i> coupling to down-type quarks.
  

<br/><br/><table><tr><td><strong>HiggsA3:coup2u </td><td></td><td> <input type="text" name="98" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>A^0(H_3^0)</i> coupling to up-type quarks.
  

<br/><br/><table><tr><td><strong>HiggsA3:coup2l </td><td></td><td> <input type="text" name="99" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>A^0(H_3^0)</i> coupling to (charged) leptons.
  

<br/><br/><table><tr><td><strong>HiggsA3:coup2H1Z </td><td></td><td> <input type="text" name="100" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>A^0(H_3^0)</i> coupling to a <i>h^0(H_1^0) Z^0</i> pair.
Is <i>cos(beta - alpha)</i> in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsA3:coup2H2Z </td><td></td><td> <input type="text" name="101" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>A^0(H_3^0)</i> coupling to a <i>H^0(H_2^0) Z^0</i> pair.
Is <i>sin(beta - alpha)</i> in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsA3:coup2Z </td><td></td><td> <input type="text" name="102" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>A^0(H_3^0)</i> coupling to <i>Z^0</i>.
Vanishes in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsA3:coup2W </td><td></td><td> <input type="text" name="103" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>A^0(H_3^0)</i> coupling to <i>W^+-</i>.
Vanishes in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsA3:coup2H1H1 </td><td></td><td> <input type="text" name="104" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>A^0(H_3^0)</i> coupling to a <i>h^0(H_1^0)</i> pair.
Vanishes in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsA3:coup2Hchg </td><td></td><td> <input type="text" name="105" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>A^0(H_3^0)</i> coupling to <i>H^+-</i>.
Vanishes in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsA3:coup2HchgW </td><td></td><td> <input type="text" name="106" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>A^0(H_3^0)</i> coupling to a <i>H^+- W-+</i> pair.
Is 1 in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsHchg:tanBeta </td><td></td><td> <input type="text" name="107" value="5." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5.</strong></code>)</td></tr></table>
The <i>tan(beta)</i> value, which leads to an enhancement of the 
<i>H^+-</i> coupling to down-type fermions and suppression to
up-type ones. The same angle also appears in many other places,
but this particular parameter is only used for the charged-Higgs case. 
  

<br/><br/><table><tr><td><strong>HiggsHchg:coup2H1W </td><td></td><td> <input type="text" name="108" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
The <i>H^+-</i> coupling to a <i>h^0(H_1^0) W^+-</i> pair.
Is <i>cos(beta - alpha)</i> in the MSSM.
  

<br/><br/><table><tr><td><strong>HiggsHchg:coup2H2W </td><td></td><td> <input type="text" name="109" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>H^+-</i> coupling to a <i>H^0(H_2^0) W^+-</i> pair.
Is <i>sin(beta - alpha)</i> in the MSSM.
  

<p/>
Another set of parameters are not used in the production stage but
exclusively for the description of angular distributions in decays.

<br/><br/><table><tr><td><strong>HiggsH1:parity  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
possibility to modify angular decay correlations in the decay of a 
<ei>h^0(H_1)</ei> decay <ei>Z^0 Z^0</ei> or <ei>W^+ W^-</ei> to four 
fermions. Currently it does not affect the partial width of the 
channels, which is only based on the above parameters.
<br/>
<input type="radio" name="110" value="0"><strong>0 </strong>: isotropic decays.<br/>
<input type="radio" name="110" value="1" checked="checked"><strong>1 </strong>: assuming the <ei>h^0(H_1)</ei> is a pure scalar  (CP-even), as in the MSSM.<br/>
<input type="radio" name="110" value="2"><strong>2 </strong>: assuming the <ei>h^0(H_1)</ei> is a pure pseudoscalar (CP-odd).<br/>
<input type="radio" name="110" value="3"><strong>3 </strong>: assuming the <ei>h^0(H_1)</ei> is a mixture of the two,  including the CP-violating interference term. The parameter <ei>eta</ei>, see below, sets the strength of the CP-odd admixture, with the interference term being proportional to <ei>eta</ei> and the CP-odd one to <ei>eta^2</ei>.<br/>

<br/><br/><table><tr><td><strong>HiggsH1:etaParity </td><td></td><td> <input type="text" name="111" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>eta</i> value of CP-violation in the 
<code>HiggsSM:parity = 3</code> option. 
  

<br/><br/><table><tr><td><strong>HiggsH2:parity  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
possibility to modify angular decay correlations in the decay of a 
<ei>H^0(H_2)</ei> decay <ei>Z^0 Z^0</ei> or <ei>W^+ W^-</ei> to four 
fermions. Currently it does not affect the partial width of the 
channels, which is only based on the above parameters.
<br/>
<input type="radio" name="112" value="0"><strong>0 </strong>: isotropic decays.<br/>
<input type="radio" name="112" value="1" checked="checked"><strong>1 </strong>: assuming the <ei>H^0(H_2)</ei> is a pure scalar  (CP-even), as in the MSSM.<br/>
<input type="radio" name="112" value="2"><strong>2 </strong>: assuming the <ei>H^0(H_2)</ei> is a pure pseudoscalar (CP-odd).<br/>
<input type="radio" name="112" value="3"><strong>3 </strong>: assuming the <ei>H^0(H_2)</ei> is a mixture of the two,  including the CP-violating interference term. The parameter <ei>eta</ei>, see below, sets the strength of the CP-odd admixture, with the interference term being proportional to <ei>eta</ei> and the CP-odd one to <ei>eta^2</ei>.<br/>

<br/><br/><table><tr><td><strong>HiggsH2:etaParity </td><td></td><td> <input type="text" name="113" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>eta</i> value of CP-violation in the 
<code>HiggsSM:parity = 3</code> option. 
  

<br/><br/><table><tr><td><strong>HiggsA3:parity  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
possibility to modify angular decay correlations in the decay of a 
<ei>A^0(H_3)</ei> decay <ei>Z^0 Z^0</ei> or <ei>W^+ W^-</ei> to four 
fermions. Currently it does not affect the partial width of the 
channels, which is only based on the above parameters.
<br/>
<input type="radio" name="114" value="0"><strong>0 </strong>: isotropic decays.<br/>
<input type="radio" name="114" value="1"><strong>1 </strong>: assuming the <ei>A^0(H_3)</ei> is a pure scalar  (CP-even).<br/>
<input type="radio" name="114" value="2" checked="checked"><strong>2 </strong>: assuming the <ei>A^0(H_3)</ei> is a pure pseudoscalar (CP-odd), as in the MSSM.<br/>
<input type="radio" name="114" value="3"><strong>3 </strong>: assuming the <ei>A^0(H_3)</ei> is a mixture of the two,  including the CP-violating interference term. The parameter <ei>eta</ei>, see below, sets the strength of the CP-odd admixture, with the interference term being proportional to <ei>eta</ei> and the CP-odd one to <ei>eta^2</ei>.<br/>

<br/><br/><table><tr><td><strong>HiggsA3:etaParity </td><td></td><td> <input type="text" name="115" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>eta</i> value of CP-violation in the 
<code>HiggsSM:parity = 3</code> option. 
  

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
$data = "Higgs:cubicWidth = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "on")
{
$data = "Higgs:runningLoopMass = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "Higgs:clipWings = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "50.")
{
$data = "Higgs:wingsFac = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "HiggsSM:all = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "HiggsSM:ffbar2H = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "HiggsSM:gg2H = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "HiggsSM:gmgm2H = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "HiggsSM:ffbar2HZ = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "HiggsSM:ffbar2HW = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "HiggsSM:ff2Hff(t:ZZ) = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "HiggsSM:ff2Hff(t:WW) = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "HiggsSM:gg2Httbar = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "HiggsSM:qqbar2Httbar = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "HiggsSM:qg2Hq = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off")
{
$data = "HiggsSM:gg2Hbbbar = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "HiggsSM:qqbar2Hbbbar = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "off")
{
$data = "HiggsSM:gg2Hg(l:t) = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "HiggsSM:qg2Hq(l:t) = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "HiggsSM:qqbar2Hg(l:t) = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "off")
{
$data = "Higgs:useBSM = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "off")
{
$data = "HiggsBSM:all = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "off")
{
$data = "HiggsBSM:allH1 = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "off")
{
$data = "HiggsBSM:ffbar2H1 = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "off")
{
$data = "HiggsBSM:gg2H1 = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "off")
{
$data = "HiggsBSM:gmgm2H1 = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "off")
{
$data = "HiggsBSM:ffbar2H1Z = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "off")
{
$data = "HiggsBSM:ffbar2H1W = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "off")
{
$data = "HiggsBSM:ff2H1ff(t:ZZ) = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "off")
{
$data = "HiggsBSM:ff2H1ff(t:WW) = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "off")
{
$data = "HiggsBSM:gg2H1ttbar = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "off")
{
$data = "HiggsBSM:qqbar2H1ttbar = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "off")
{
$data = "HiggsBSM:allH2 = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "off")
{
$data = "HiggsBSM:ffbar2H2 = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "off")
{
$data = "HiggsBSM:gg2H2 = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "off")
{
$data = "HiggsBSM:gmgm2H2 = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "off")
{
$data = "HiggsBSM:ffbar2H2Z = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "off")
{
$data = "HiggsBSM:ffbar2H2W = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "off")
{
$data = "HiggsBSM:ff2H2ff(t:ZZ) = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
if($_POST["40"] != "off")
{
$data = "HiggsBSM:ff2H2ff(t:WW) = ".$_POST["40"]."\n";
fwrite($handle,$data);
}
if($_POST["41"] != "off")
{
$data = "HiggsBSM:gg2H2ttbar = ".$_POST["41"]."\n";
fwrite($handle,$data);
}
if($_POST["42"] != "off")
{
$data = "HiggsBSM:qqbar2H2ttbar = ".$_POST["42"]."\n";
fwrite($handle,$data);
}
if($_POST["43"] != "off")
{
$data = "HiggsBSM:allA3 = ".$_POST["43"]."\n";
fwrite($handle,$data);
}
if($_POST["44"] != "off")
{
$data = "HiggsBSM:ffbar2A3 = ".$_POST["44"]."\n";
fwrite($handle,$data);
}
if($_POST["45"] != "off")
{
$data = "HiggsBSM:gg2A3 = ".$_POST["45"]."\n";
fwrite($handle,$data);
}
if($_POST["46"] != "off")
{
$data = "HiggsBSM:gmgm2A3 = ".$_POST["46"]."\n";
fwrite($handle,$data);
}
if($_POST["47"] != "off")
{
$data = "HiggsBSM:ffbar2A3Z = ".$_POST["47"]."\n";
fwrite($handle,$data);
}
if($_POST["48"] != "off")
{
$data = "HiggsBSM:ffbar2A3W = ".$_POST["48"]."\n";
fwrite($handle,$data);
}
if($_POST["49"] != "off")
{
$data = "HiggsBSM:ff2A3ff(t:ZZ) = ".$_POST["49"]."\n";
fwrite($handle,$data);
}
if($_POST["50"] != "off")
{
$data = "HiggsBSM:ff2A3ff(t:WW) = ".$_POST["50"]."\n";
fwrite($handle,$data);
}
if($_POST["51"] != "off")
{
$data = "HiggsBSM:gg2A3ttbar = ".$_POST["51"]."\n";
fwrite($handle,$data);
}
if($_POST["52"] != "off")
{
$data = "HiggsBSM:qqbar2A3ttbar = ".$_POST["52"]."\n";
fwrite($handle,$data);
}
if($_POST["53"] != "off")
{
$data = "HiggsBSM:allH+- = ".$_POST["53"]."\n";
fwrite($handle,$data);
}
if($_POST["54"] != "off")
{
$data = "HiggsBSM:ffbar2H+- = ".$_POST["54"]."\n";
fwrite($handle,$data);
}
if($_POST["55"] != "off")
{
$data = "HiggsBSM:bg2H+-t = ".$_POST["55"]."\n";
fwrite($handle,$data);
}
if($_POST["56"] != "off")
{
$data = "HiggsBSM:allHpair = ".$_POST["56"]."\n";
fwrite($handle,$data);
}
if($_POST["57"] != "off")
{
$data = "HiggsBSM:ffbar2A3H1 = ".$_POST["57"]."\n";
fwrite($handle,$data);
}
if($_POST["58"] != "off")
{
$data = "HiggsBSM:ffbar2A3H2 = ".$_POST["58"]."\n";
fwrite($handle,$data);
}
if($_POST["59"] != "off")
{
$data = "HiggsBSM:ffbar2H+-H1 = ".$_POST["59"]."\n";
fwrite($handle,$data);
}
if($_POST["60"] != "off")
{
$data = "HiggsBSM:ffbar2H+-H2 = ".$_POST["60"]."\n";
fwrite($handle,$data);
}
if($_POST["61"] != "off")
{
$data = "HiggsBSM:ffbar2H+H- = ".$_POST["61"]."\n";
fwrite($handle,$data);
}
if($_POST["62"] != "off")
{
$data = "HiggsBSM:qg2H1q = ".$_POST["62"]."\n";
fwrite($handle,$data);
}
if($_POST["63"] != "off")
{
$data = "HiggsBSM:gg2H1bbbar = ".$_POST["63"]."\n";
fwrite($handle,$data);
}
if($_POST["64"] != "off")
{
$data = "HiggsBSM:qqbar2H1bbbar = ".$_POST["64"]."\n";
fwrite($handle,$data);
}
if($_POST["65"] != "off")
{
$data = "HiggsBSM:gg2H1g(l:t) = ".$_POST["65"]."\n";
fwrite($handle,$data);
}
if($_POST["66"] != "off")
{
$data = "HiggsBSM:qg2H1q(l:t) = ".$_POST["66"]."\n";
fwrite($handle,$data);
}
if($_POST["67"] != "off")
{
$data = "HiggsBSM:qqbar2H1g(l:t) = ".$_POST["67"]."\n";
fwrite($handle,$data);
}
if($_POST["68"] != "off")
{
$data = "HiggsBSM:qg2H2q = ".$_POST["68"]."\n";
fwrite($handle,$data);
}
if($_POST["69"] != "off")
{
$data = "HiggsBSM:gg2H2bbbar = ".$_POST["69"]."\n";
fwrite($handle,$data);
}
if($_POST["70"] != "off")
{
$data = "HiggsBSM:qqbar2H2bbbar = ".$_POST["70"]."\n";
fwrite($handle,$data);
}
if($_POST["71"] != "off")
{
$data = "HiggsBSM:gg2H2g(l:t) = ".$_POST["71"]."\n";
fwrite($handle,$data);
}
if($_POST["72"] != "off")
{
$data = "HiggsBSM:qg2H2q(l:t) = ".$_POST["72"]."\n";
fwrite($handle,$data);
}
if($_POST["73"] != "off")
{
$data = "HiggsBSM:qqbar2H2g(l:t) = ".$_POST["73"]."\n";
fwrite($handle,$data);
}
if($_POST["74"] != "off")
{
$data = "HiggsBSM:qg2A3q = ".$_POST["74"]."\n";
fwrite($handle,$data);
}
if($_POST["75"] != "off")
{
$data = "HiggsBSM:gg2A3bbbar = ".$_POST["75"]."\n";
fwrite($handle,$data);
}
if($_POST["76"] != "off")
{
$data = "HiggsBSM:qqbar2A3bbbar = ".$_POST["76"]."\n";
fwrite($handle,$data);
}
if($_POST["77"] != "off")
{
$data = "HiggsBSM:gg2A3g(l:t) = ".$_POST["77"]."\n";
fwrite($handle,$data);
}
if($_POST["78"] != "off")
{
$data = "HiggsBSM:qg2A3q(l:t) = ".$_POST["78"]."\n";
fwrite($handle,$data);
}
if($_POST["79"] != "off")
{
$data = "HiggsBSM:qqbar2A3g(l:t) = ".$_POST["79"]."\n";
fwrite($handle,$data);
}
if($_POST["80"] != "1.")
{
$data = "HiggsH1:coup2d = ".$_POST["80"]."\n";
fwrite($handle,$data);
}
if($_POST["81"] != "1.")
{
$data = "HiggsH1:coup2u = ".$_POST["81"]."\n";
fwrite($handle,$data);
}
if($_POST["82"] != "1.")
{
$data = "HiggsH1:coup2l = ".$_POST["82"]."\n";
fwrite($handle,$data);
}
if($_POST["83"] != "1.")
{
$data = "HiggsH1:coup2Z = ".$_POST["83"]."\n";
fwrite($handle,$data);
}
if($_POST["84"] != "1.")
{
$data = "HiggsH1:coup2W = ".$_POST["84"]."\n";
fwrite($handle,$data);
}
if($_POST["85"] != "0.")
{
$data = "HiggsH1:coup2Hchg = ".$_POST["85"]."\n";
fwrite($handle,$data);
}
if($_POST["86"] != "1.")
{
$data = "HiggsH2:coup2d = ".$_POST["86"]."\n";
fwrite($handle,$data);
}
if($_POST["87"] != "1.")
{
$data = "HiggsH2:coup2u = ".$_POST["87"]."\n";
fwrite($handle,$data);
}
if($_POST["88"] != "1.")
{
$data = "HiggsH2:coup2l = ".$_POST["88"]."\n";
fwrite($handle,$data);
}
if($_POST["89"] != "1.")
{
$data = "HiggsH2:coup2Z = ".$_POST["89"]."\n";
fwrite($handle,$data);
}
if($_POST["90"] != "1.")
{
$data = "HiggsH2:coup2W = ".$_POST["90"]."\n";
fwrite($handle,$data);
}
if($_POST["91"] != "0.")
{
$data = "HiggsH2:coup2Hchg = ".$_POST["91"]."\n";
fwrite($handle,$data);
}
if($_POST["92"] != "1.")
{
$data = "HiggsH2:coup2H1H1 = ".$_POST["92"]."\n";
fwrite($handle,$data);
}
if($_POST["93"] != "1.")
{
$data = "HiggsH2:coup2A3A3 = ".$_POST["93"]."\n";
fwrite($handle,$data);
}
if($_POST["94"] != "0.")
{
$data = "HiggsH2:coup2H1Z = ".$_POST["94"]."\n";
fwrite($handle,$data);
}
if($_POST["95"] != "0.")
{
$data = "HiggsH2:coup2A3H1 = ".$_POST["95"]."\n";
fwrite($handle,$data);
}
if($_POST["96"] != "0.")
{
$data = "HiggsH2:coup2HchgW = ".$_POST["96"]."\n";
fwrite($handle,$data);
}
if($_POST["97"] != "1.")
{
$data = "HiggsA3:coup2d = ".$_POST["97"]."\n";
fwrite($handle,$data);
}
if($_POST["98"] != "1.")
{
$data = "HiggsA3:coup2u = ".$_POST["98"]."\n";
fwrite($handle,$data);
}
if($_POST["99"] != "1.")
{
$data = "HiggsA3:coup2l = ".$_POST["99"]."\n";
fwrite($handle,$data);
}
if($_POST["100"] != "1.")
{
$data = "HiggsA3:coup2H1Z = ".$_POST["100"]."\n";
fwrite($handle,$data);
}
if($_POST["101"] != "1.")
{
$data = "HiggsA3:coup2H2Z = ".$_POST["101"]."\n";
fwrite($handle,$data);
}
if($_POST["102"] != "0.")
{
$data = "HiggsA3:coup2Z = ".$_POST["102"]."\n";
fwrite($handle,$data);
}
if($_POST["103"] != "0.")
{
$data = "HiggsA3:coup2W = ".$_POST["103"]."\n";
fwrite($handle,$data);
}
if($_POST["104"] != "0.")
{
$data = "HiggsA3:coup2H1H1 = ".$_POST["104"]."\n";
fwrite($handle,$data);
}
if($_POST["105"] != "0.")
{
$data = "HiggsA3:coup2Hchg = ".$_POST["105"]."\n";
fwrite($handle,$data);
}
if($_POST["106"] != "1.")
{
$data = "HiggsA3:coup2HchgW = ".$_POST["106"]."\n";
fwrite($handle,$data);
}
if($_POST["107"] != "5.")
{
$data = "HiggsHchg:tanBeta = ".$_POST["107"]."\n";
fwrite($handle,$data);
}
if($_POST["108"] != "1.")
{
$data = "HiggsHchg:coup2H1W = ".$_POST["108"]."\n";
fwrite($handle,$data);
}
if($_POST["109"] != "0.")
{
$data = "HiggsHchg:coup2H2W = ".$_POST["109"]."\n";
fwrite($handle,$data);
}
if($_POST["110"] != "1")
{
$data = "HiggsH1:parity = ".$_POST["110"]."\n";
fwrite($handle,$data);
}
if($_POST["111"] != "0.")
{
$data = "HiggsH1:etaParity = ".$_POST["111"]."\n";
fwrite($handle,$data);
}
if($_POST["112"] != "1")
{
$data = "HiggsH2:parity = ".$_POST["112"]."\n";
fwrite($handle,$data);
}
if($_POST["113"] != "0.")
{
$data = "HiggsH2:etaParity = ".$_POST["113"]."\n";
fwrite($handle,$data);
}
if($_POST["114"] != "2")
{
$data = "HiggsA3:parity = ".$_POST["114"]."\n";
fwrite($handle,$data);
}
if($_POST["115"] != "0.")
{
$data = "HiggsA3:etaParity = ".$_POST["115"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->

