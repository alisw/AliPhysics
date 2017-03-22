<html>
<head>
<title>Hidden Valley Processes</title>
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

<form method='post' action='HiddenValleyProcesses.php'>
 
<h2>Hidden Valley Processes</h2> 
 
This Hidden Valley (HV) scenario has been developed specifically 
to allow the study of visible consequences of radiation in a 
hidden sector, by recoil effect. A key aspect is therefore that 
the normal timelike showering machinery has been expanded with a 
third kind of radiation, in addition to the QCD and QED ones. 
These three kinds of radiation are fully interleaved, i.e. 
evolution occurs in a common <i>pT</i>-ordered sequence. 
The scenario is described in [<a href="Bibliography.php" target="page">Car10</a>]. Furthermore 
hadronization in the hidden sector has been implemented. 
Three main scenarios for production into and decay out of the 
hidden sector can be compared, in each case either for an 
Abelian or a non-Abelian gauge group in the HV. For further details 
see [<a href="Bibliography.php" target="page">Car11</a>]. 
 
<h3>Particle content and properties</h3> 
 
For simplicity we assume that the HV contains an unbroken <b>SU(N)</b> 
gauge symmetry. This is used in the calculation of production cross 
sections. These could be rescaled by hand for other gauge groups. 
 
<br/><br/><table><tr><td><strong>HiddenValley:Ngauge  </td><td></td><td> <input type="text" name="1" value="3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 1</code>)</td></tr></table>
is <b>U(1)</b> for <code>Ngauge = 1</code>, is <b>SU(N)</b> if 
<code>Ngauge &gt; 1</code>. Note that pair production cross sections 
contains a factor of <code>Ngauge</code> for new particles 
in the fundamental representation of this group. 
   
 
<p/> 
A minimal HV particle content has been introduced. Firstly, there is 
a set of 12 particles that mirrors the Standard Model flavour 
structure, and is charged under both the SM and the HV symmetry groups. 
Each new particle couples flavour-diagonally to a corresponding SM 
state, and has the same SM charge and colour, but in addition is in 
the fundamental representation of the HV colour, as follows: 
<br/><code>Dv</code>, identity 4900001, partner to the normal 
<code>d</code> quark; 
<br/><code>Uv</code>, identity 4900002, partner to the normal 
<code>u</code> quark; 
<br/><code>Sv</code>, identity 4900003, partner to the normal 
<code>s</code> quark; 
<br/><code>Cv</code>, identity 4900004, partner to the normal 
<code>c</code> quark; 
<br/><code>Bv</code>, identity 4900005, partner to the normal 
<code>b</code> quark; 
<br/><code>Tv</code>, identity 4900006, partner to the normal 
<code>t</code> quark; 
<br/><code>Ev</code>, identity 4900011, partner to the normal 
<code>e</code> lepton; 
<br/><code>nuEv</code>, identity 4900012, partner to the normal 
<code>nue</code> neutrino; 
<br/><code>MUv</code>, identity 4900013, partner to the normal 
<code>mu</code> lepton; 
<br/><code>nuMUv</code>, identity 4900014, partner to the normal 
<code>numu</code> neutrino; 
<br/><code>TAUv</code>, identity 4900015, partner to the normal 
<code>tau</code> lepton; 
<br/><code>nuTAUv</code>, identity 4900016, partner to the normal 
<code>nutau</code> neutrino. 
<br/>Collectively we will refer to these states as <code>Fv</code>; 
note, however, that they need not be fermions themselves. 
 
<p/> 
In addition the model contains the HV gauge particle, either 
a HV-gluon or a HV-photon, but not both; see <code>Ngauge</code> 
above: 
<br/><code>gv</code>, identity 4900021, is the massless 
gauge boson of the HV <b>SU(N)</b> group; 
<br/><code>gammav</code>, identity 4900022,  is the massless 
gauge boson of the HV <b>U(1)</b> group. 
 
<p/> 
Finally, for the basic HV scenario, there is a new massive particle 
with only HV charge sitting in the fundamental representation of the 
HV gauge group: 
<br/><code>qv</code>, identity 4900101. 
 
<p/>The typical scenario would be for pair production of one of the 
states presented first above, e.g. <i>g g &rarr; Dv Dvbar</i>. 
Such a <i>Dv</i> can radiate gluons and photons like an SM quark, 
but in addition HV-gluons or HV-photons in a similar fashion. 
Eventually the <i>Dv</i> will decay like <i>Dv &rarr; d + qv</i>. 
The strength of this decay is not set as such, but is implicit in 
your choice of width for the <i>Dv</i> state. Thereafter the 
<i>d</i> and <i>qv</i> can radiate further within their 
respective sectors. The <i>qv</i>, <i>gv</i> or <i>gammav</i> 
are invisible, so their fate need not be considered further. 
 
<p/> 
While not part of the standard scenario, as an alternative there is 
also a kind of <i>Z'</i> resonance: 
<br/><code>Zv</code>, identity 4900023, a boson that can couple 
both to pairs of Standard Model fermions and to <i>qv qvbar</i> 
pairs. Mass, total width and branching ratios can be set as convenient. 
<br/>This opens up for alternative processes 
<i>l^+l^-, q qbar &rarr; Zv &rarr; qv qvbar</i>. 
 
<p/> 
The possibility of a leakage back from the hidden sector will be 
considered in the Hadronization section below. For the <b>U(1)</b> 
case the <i>gammav</i> acquires a mass and can decay back to a 
Standard-Model fermion pair, while the <i>qv</i> remains invisible. 
The <b>SU(N)</b> alternative remains unbroken, so confinement holds 
and the <i>gv</i> is massless. A string like 
<i>qv - gv - ... - gv - qvbar</i> can break by the production of 
new <i>qv - qvbar</i> pairs, which will produce <i>qv-qvbar</i> 
mesons. It would be possible to build a rather sophisticated hidden 
sector by trivial extensions of the HV flavour content. For now, 
however, the <i>qv</i> can be duplicated in up to eight copies 
with the same properties except for the flavour charge. These are 
assigned codes 4900101 - 4900108. This gives a total of 64 possible 
lowest-lying mesons. We also include a duplication of that, into two 
multiplets, corresponding to the pseudoscalar and vector mesons of 
QCD. For now, again, these are assumed to have the same mass and 
other properties. Only the flavour-diagonal ones can decay back into 
the Standard-Model sector, however, while the rest remains in the 
hidden sector. It is therefore only necessary to distinguish a few 
states: 
<br/><code>pivDiag</code>, identity 4900111, a flavour-diagonal 
HV-meson with spin 0 that can decay back into the Standard-Model sector; 
<br/><code>rhovDiag</code>, identity 4900113, a flavour-diagonal 
HV-meson with spin 1 that can decay back into the Standard-Model sector; 
<br/><code>pivUp</code>, identity 4900211, an off-diagonal 
HV-meson with spin 0 that is stable and invisible, with an antiparticle 
<code>pivDn</code> with identity -4900211; the particle is 
the one where the code of the flavour is larger than that of the 
antiflavour; 
<br/><code>rhovUp</code>, identity 4900213, an off-diagonal 
HV-meson with spin 1 that is stable and invisible, with an antiparticle 
<code>rhovDn</code> with identity -4900213; again the particle is 
the one where the code of the flavour is larger than that of the 
antiflavour; 
<br/><code>ggv</code>, identity 4900991, is only rarely used, 
to handle cases where it is kinematically impossible to produce an 
HV-meson on shell, and it therefore is assumed to de-excite by the 
emission of invisible <i>gv-gv </i> v-glueball bound states. 
 
<p/> 
Only the spin of the HV-gluon or HV-photon is determined unambiguously 
to be unity, for the others you can make your choice. 
 
<br/><br/><table><tr><td><strong>HiddenValley:spinFv  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
The spin of the HV partners of the SM fermions, e.g. 
<ei>Dv</ei>, <ei>Uv</ei>, <ei>Ev</ei> and <ei>nuEv</ei>. 
<br/>
<input type="radio" name="2" value="0"><strong>0 </strong>: spin 0.<br/>
<input type="radio" name="2" value="1" checked="checked"><strong>1 </strong>: spin 1/2.<br/>
<input type="radio" name="2" value="2"><strong>2 </strong>: spin 1.<br/>
 
<br/><br/><table><tr><td><strong>HiddenValley:spinqv  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
The spin of <ei>qv</ei> when the <ei>Fv</ei> (the HV partners of 
the SM fermions) have spin 1/2. (While, if they have spin 0 or 1, 
the <ei>qv</ei> spin is fixed at 1/2.) 
<br/>
<input type="radio" name="3" value="0" checked="checked"><strong>0 </strong>: spin 0.<br/>
<input type="radio" name="3" value="1"><strong>1 </strong>: spin 1.<br/>
 
<br/><br/><table><tr><td><strong>HiddenValley:kappa </td><td></td><td> <input type="text" name="4" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
If the <i>Fv</i> have spin 1 then their production 
cross section depends on the presence of anomalous magnetic dipole 
moment, i.e. of a <i>kappa</i> different from unity. For other spins 
this parameter is not used. 
   
 
<br/><br/><strong>HiddenValley:doKinMix</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
allow kinematic mixing or not. 
   
 
<br/><br/><table><tr><td><strong>HiddenValley:kinMix </td><td></td><td> <input type="text" name="6" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
strength of kinetic mixing. 
   
 
<p/> 
You should set the <i>Fv</i> and <i>qv</i> masses appropriately, 
with the latter smaller than the former two to allow decays. 
When <b>U(1)</b> hadronization is switched on, you need to set the 
<i>gammav</i> mass and decay modes. For <b>SU(N)</b> hadronization 
the HV-meson masses should be set to match the <i>qv</i> ones. 
The simplest is to assume that <i>m_qv</i> defines a constituent 
mass, so that  <i>m_HVmeson = 2 m_qv</i>. The <i>hvMesonDiag</i> 
decay modes also need to be set. 
 
<h3>Production processes</h3> 
 
<br/><br/><strong>HiddenValley:all</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of all hard Hidden Valley processes, 
as listed separately in the following. 
   
 
<br/><br/><strong>HiddenValley:gg2DvDvbar</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>g g &rarr; Dv Dvbar</i>. 
Code 4901. 
   
 
<br/><br/><strong>HiddenValley:gg2UvUvbar</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>g g &rarr; Uv Uvbar</i>. 
Code 4902. 
   
 
<br/><br/><strong>HiddenValley:gg2SvSvbar</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>g g &rarr; Sv Svbar</i>. 
Code 4903. 
   
 
<br/><br/><strong>HiddenValley:gg2CvCvbar</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>g g &rarr; Cv Cvbar</i>. 
Code 4904. 
   
 
<br/><br/><strong>HiddenValley:gg2BvBvbar</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>g g &rarr; Bv Bvbar</i>. 
Code 4905. 
   
 
<br/><br/><strong>HiddenValley:gg2TvTvbar</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>g g &rarr; Tv Tvbar</i>. 
Code 4906. 
   
 
<br/><br/><strong>HiddenValley:qqbar2DvDvbar</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>q qbar &rarr; Dv Dvbar</i> 
via intermediate gluon. 
Code 4911. 
   
 
<br/><br/><strong>HiddenValley:qqbar2UvUvbar</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>q qbar &rarr; Uv Uvbar</i> 
via intermediate gluon. 
Code 4912. 
   
 
<br/><br/><strong>HiddenValley:qqbar2SvSvbar</strong>  <input type="radio" name="16" value="on"><strong>On</strong>
<input type="radio" name="16" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>q qbar &rarr; Sv Svbar</i> 
via intermediate gluon. 
Code 4913. 
   
 
<br/><br/><strong>HiddenValley:qqbar2CvCvbar</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>q qbar &rarr; Cv Cvbar</i> 
via intermediate gluon. 
Code 4914. 
   
 
<br/><br/><strong>HiddenValley:qqbar2BvBvbar</strong>  <input type="radio" name="18" value="on"><strong>On</strong>
<input type="radio" name="18" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>q qbar &rarr; Bv Bvbar</i> 
via intermediate gluon. 
Code 4915. 
   
 
<br/><br/><strong>HiddenValley:qqbar2TvTvbar</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>q qbar &rarr; Tv Tvbar</i> 
via intermediate gluon. 
Code 4916. 
   
 
<br/><br/><strong>HiddenValley:ffbar2DvDvbar</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; Dv Dvbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4921. 
   
 
<br/><br/><strong>HiddenValley:ffbar2UvUvbar</strong>  <input type="radio" name="21" value="on"><strong>On</strong>
<input type="radio" name="21" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; Uv Uvbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4922. 
   
 
<br/><br/><strong>HiddenValley:ffbar2SvSvbar</strong>  <input type="radio" name="22" value="on"><strong>On</strong>
<input type="radio" name="22" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; Sv Svbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4923. 
   
 
<br/><br/><strong>HiddenValley:ffbar2CvCvbar</strong>  <input type="radio" name="23" value="on"><strong>On</strong>
<input type="radio" name="23" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; Cv Cvbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4924. 
   
 
<br/><br/><strong>HiddenValley:ffbar2BvBvbar</strong>  <input type="radio" name="24" value="on"><strong>On</strong>
<input type="radio" name="24" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; Bv Bvbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4925. 
   
 
<br/><br/><strong>HiddenValley:ffbar2TvTvbar</strong>  <input type="radio" name="25" value="on"><strong>On</strong>
<input type="radio" name="25" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; Tv Tvbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4926. 
   
 
<br/><br/><strong>HiddenValley:ffbar2EvEvbar</strong>  <input type="radio" name="26" value="on"><strong>On</strong>
<input type="radio" name="26" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; Ev Evbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4931. 
   
 
<br/><br/><strong>HiddenValley:ffbar2nuEvnuEvbar</strong>  <input type="radio" name="27" value="on"><strong>On</strong>
<input type="radio" name="27" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; nuEv nuEvbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4932. 
   
 
<br/><br/><strong>HiddenValley:ffbar2MUvMUvbar</strong>  <input type="radio" name="28" value="on"><strong>On</strong>
<input type="radio" name="28" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; MUv MUvbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4933. 
   
 
<br/><br/><strong>HiddenValley:ffbar2nuMUvnuMUvbar</strong>  <input type="radio" name="29" value="on"><strong>On</strong>
<input type="radio" name="29" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; nuMUv nuMUvbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4934. 
   
 
<br/><br/><strong>HiddenValley:ffbar2TAUvTAUvbar</strong>  <input type="radio" name="30" value="on"><strong>On</strong>
<input type="radio" name="30" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; TAUv TAUvbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4935. 
   
 
<br/><br/><strong>HiddenValley:ffbar2nuTAUvnuTAUvbar</strong>  <input type="radio" name="31" value="on"><strong>On</strong>
<input type="radio" name="31" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production <i>f fbar &rarr; nuTAUv nuTAUvbar</i> 
via intermediate <i>gamma*/Z^*</i>. 
Code 4936. 
   
 
<br/><br/><strong>HiddenValley:ffbar2Zv</strong>  <input type="radio" name="32" value="on"><strong>On</strong>
<input type="radio" name="32" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Production <i>f fbar &rarr; Zv</i> where <i>Zv</i> is a generic 
resonance that couples both SM fermion pairs and a <i>qv qvbar</i> 
pair. Not part of the framework of the above processes, but as an 
alternative. Code 4941. 
   
 
<h3>Timelike showers</h3> 
 
One key point of this HV scenario is that radiation off the 
HV-charged particles is allowed. This is done by the standard 
final-state showering machinery. (HV particles are not produced 
in initial-state radiation.) All the (anti)particles <i>Fv</i> 
and <i>qv</i> have one (negative) unit of HV charge. That is, 
radiation closely mimics the one in QCD. Both QCD, QED and HV 
radiation are interleaved in one common sequence of decreasing 
emission <i>pT</i> scales. Each radiation kind defines a set of 
dipoles, usually spanned between a radiating parton and its recoil 
partner, such that the invariant mass of the pair is not changed 
when a radiation occurs. This need not follow from trivial colour 
assignments, but is often obvious. For instance,  in a decay 
<i>Qv &rarr; q + qv</i> the QCD dipole is between the <i>q</i> and 
the hole after <i>Qv</i>, but <i>qv</i> becomes the recoiler 
should a radiation occur, while the role of <i>q</i> and <i>qv</i> 
is reversed for HV radiation. 
 
<p/>This also includes matrix-element corrections for a number 
of decay processes, with colour, spin and mass effects included 
[<a href="Bibliography.php" target="page">Nor01</a>]. They were calculated within the context of the 
particle content of the MSSM, however, which does not include spin 1 
particles with unit colour charge. In such cases spin 0 is assumed 
instead. By experience, the main effects come from mass and colour 
flow anyway, so this is not a bad approximation. (Furthermore the 
MSSM formulae allow for <i>gamma_5</i> factors from wave 
functions or vertices; these are even less important.) 
 
<p/>An emitted <i>gv</i> can branch in its turn, 
<i>gv &rarr; gv + gv</i>. This radiation may affect momenta 
in the visible sector by recoil effect, but this is a minor 
effect relative to the primary emission of the <i>gv</i>. 
 
<br/><br/><strong>HiddenValley:FSR</strong>  <input type="radio" name="33" value="on"><strong>On</strong>
<input type="radio" name="33" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
switch on final-state shower of <i>gv</i> or <i>gammav</i> 
in a HV production process. 
   
 
<br/><br/><table><tr><td><strong>HiddenValley:alphaFSR </td><td></td><td> <input type="text" name="34" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
fixed alpha scale of <i>gv/gammav</i> emission; corresponds to 
<i>alpha_strong</i> of QCD or <i>alpha_em</i> of QED. For 
shower branchings such as <i>Dv &rarr; Dv + gv</i> the coupling is 
multiplied by <i>C_F = (N^2 - 1) / (2 * N)</i> for an 
<b>SU(N)</b> group and for <i>gv &rarr; gv + gv</i> by <i>N</i>. 
   
 
<br/><br/><table><tr><td><strong>HiddenValley:pTminFSR </td><td></td><td> <input type="text" name="35" value="0.4" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.4</strong></code>; <code>minimum = 0.1</code>)</td></tr></table>
lowest allowed <i>pT</i> of emission. Chosen with same default 
as in normal QCD showers. 
   
 
<h3>Hadronization</h3> 
 
By default the HV particles with no Standard Model couplings 
are not visible. Their presence can only be deduced by the 
observation of missing (transverse) momentum in the event as a 
whole. In the current implementation it is possible to simulate 
two different scenarios where activity can leak back from the 
hidden sector. 
 
<p/> 
The first possibility is relevant for the <b>U(1)</b> scenario. 
The <b>U(1)</b> group may be broken, so that the <i>gammav</i> 
acquires a mass. Furthermore, the <i>gammav</i> may have a 
small mixing angle with the normal photon, or with some <i>Z'</i> 
state or other mediator, and may thus decay back into Standard 
Model particles. The <i>qv</i> still escapes undetected; 
recall that there is no confinement in the <b>U(1)</b> option. 
 
<p/> 
In order to enable this machinery two commands are necessary, 
<code>4900022:m0 = ...</code> to set the  <i>gammav</i> mass 
to the desired value, and <code>4900022:onMode = on</code> to enable 
<i>gammav</i> decays. The default <i>gammav</i> decay 
table contains all Standard Model fermion-antifermion pairs, 
except top, with branching ratios in proportion to their coupling 
to the photon, whenever the production channel is allowed by 
kinematics. This table could easily be tailored to more specific 
models and needs. For instance, for a mass below 1 - 2 GeV, it 
would make sense to construct a table of exclusive hadronic decay 
channels rather than go the way via a hadronizing quark pair. 
 
<p/> 
The <i>gammav</i> are expected to decay so rapidly that no 
secondary vertex will be detectable. However, it is possible to 
set <code>4900022:tau0</code> to a finite lifetime (in mm) that 
will be used to create separated secondary vertices. 
 
<p/> 
The second, more interesting, possibility is relevant for the 
<b>SU(N)</b> scenarios. Here the gauge group remains unbroken, i.e. 
<i>gv</i> is massless, and the partons are confined. Like in 
QCD, the HV-partons can therefore be arranged in one single 
HV-colour-ordered chain, with a <i>qv</i> in one end, a 
<i>qvbar</i> in the other, and a varying number of 
<i>gv</i> in between. Each event will only contain (at most) 
one such string, (i) since perturbative branchings 
<i>gv &rarr; qv qvbar</i> have been neglected, as is a reasonable 
approximation for QCD, and (ii) since HV-colours are assigned in the 
<i>N_C &rarr; infinity</i> limit, just like in the handling of 
string fragmentation in QCD. The HV-string can then fragment by the 
nonperturbative creation of <i>qv qvbar</i> pairs, leading to 
the formation of HV-mesons along the string, each with its 
<i>qv</i> from one vertex and its <i>qvbar</i> from 
the neighbouring one. 
 
<p/> 
Since, so far, we have only assumed there to be one <i>qv</i> 
species, all produced <i>qv qvbar</i> HV-mesons are of the 
same flavour-diagonal species. Such an HV-meson can decay back to 
the normal sector, typically by whatever mediator particle allowed 
production in the first place. In this framework the full energy put 
into the HV sector will leak back to the normal one. To allow more 
flexibility, an ad hoc possibility of <i>n_Flav</i> different 
<i>qv</i> species is introduced. For now they are all assumed 
to have the same mass and other properties, but distinguished by 
some flavour-like property. Only the flavour-diagonal ones can decay, 
meaning that only a fraction (approximately) <i>1/n_Flav</i> of the 
HV-energy leaks back, while the rest remains in the hidden sector. 
 
<p/> 
This scenario contains more parameters than the first one, for the 
<b>U(1)</b> group. They can be subdivided into two sets. One is 
related to particle properties, both for <i>qv</i> and for the 
two different kinds of HV-mesons, here labeled 4900111 and 4900113 
for the diagonal ones, and +-4900211 and +-4900213 for the 
off-diagonal ones. It makes sense to set the HV-meson masses to be 
twice the <i>qv</i> one, as in a simple constituent mass context. 
Furthermore the <i>hvMesonDiag</i> decay modes need to be set up. 
Like with the 
<i>gammav</i> in the <b>U(1)</b> option, the default decay table 
is based on the branching ratios of an off-shell photon. 
 
<p/> 
The second set are fragmentation parameters that extend or replace 
the ones used in normal string fragmentation. Some of them are not 
encoded in the same way as normally, however, but rather scale as 
the <i>qv</i> mass is changed, so as to keep a sensible default 
behaviour. This does not mean that deviations from this set should 
not be explored, or that other scaling rules could be preferred 
within alternative scenarios. These parameters are as follows. 
 
<br/><br/><strong>HiddenValley:fragment</strong>  <input type="radio" name="36" value="on"><strong>On</strong>
<input type="radio" name="36" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
switch on string fragmentation of the HV partonic system. 
Only relevant for <b>SU(N)</b> scenarios. 
   
 
<br/><br/><table><tr><td><strong>HiddenValley:nFlav  </td><td></td><td> <input type="text" name="37" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 8</code>)</td></tr></table>
number of different flavours assumed to exist in the hadronization 
description, leading to approximately <i>1/n_Flav</i> of the 
produced HV-mesons being flavour-diagonal and capable to decay back 
to Standard Model particles. 
   
 
<br/><br/><table><tr><td><strong>HiddenValley:probVector </td><td></td><td> <input type="text" name="38" value="0.75" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.75</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
fraction of HV-mesons that are assigned spin 1 (vector), with the 
remainder spin 0 (pseudoscalar). Assuming the <i>qv</i> have 
spin <i>1/2</i> and the mass splitting is small, spin counting 
predicts that <i>3/4</i> of the mesons should have spin 1. 
   
 
<br/><br/><table><tr><td><strong>HiddenValley:aLund </td><td></td><td> <input type="text" name="39" value="0.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.3</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
The <i>a</i> parameter of the Lund symmetric fragmentation function. 
See the normal <?php $filepath = $_GET["filepath"];
echo "<a href='Fragmentation.php?filepath=".$filepath."' target='page'>";?>fragmentation 
function</a> description for the shape of this function.   
 
<br/><br/><table><tr><td><strong>HiddenValley:bmqv2 </td><td></td><td> <input type="text" name="40" value="0.8" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.8</strong></code>; <code>minimum = 0.2</code>; <code>maximum = 2.0</code>)</td></tr></table>
The <i>b</i> parameter of the Lund symmetric fragmentation function, 
multiplied by the square of the <i>qv</i> mass. This scaling ensures 
that the fragmentation function keeps the same shape when the 
<i>qv</i> mass is changed (neglecting transverse momenta). 
   
 
<br/><br/><table><tr><td><strong>HiddenValley:rFactqv </td><td></td><td> <input type="text" name="41" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
<i>r_qv</i>, i.e. the Bowler correction factor to the Lund symmetric 
fragmentation function, which could be made weaker or stronger than 
its natural value. 
   
 
<br/><br/><table><tr><td><strong>HiddenValley:sigmamqv </td><td></td><td> <input type="text" name="42" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
the width <i>sigma</i> of transverse momenta in the HV fragmentation 
process, normalized to the <i>qv</i> mass. This ensures that 
<i>sigma</i> scales proportionately to <i>m_qv</i>. 
See the normal <?php $filepath = $_GET["filepath"];
echo "<a href='Fragmentation.php?filepath=".$filepath."' target='page'>";?>fragmentation 
<i>pT</i></a> description for conventions for factors of 2. 
   
 
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

if($_POST["1"] != "3")
{
$data = "HiddenValley:Ngauge = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "1")
{
$data = "HiddenValley:spinFv = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "0")
{
$data = "HiddenValley:spinqv = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "1.")
{
$data = "HiddenValley:kappa = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "HiddenValley:doKinMix = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "1.")
{
$data = "HiddenValley:kinMix = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "HiddenValley:all = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "HiddenValley:gg2DvDvbar = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "HiddenValley:gg2UvUvbar = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "HiddenValley:gg2SvSvbar = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "HiddenValley:gg2CvCvbar = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "HiddenValley:gg2BvBvbar = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "HiddenValley:gg2TvTvbar = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "HiddenValley:qqbar2DvDvbar = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "HiddenValley:qqbar2UvUvbar = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off")
{
$data = "HiddenValley:qqbar2SvSvbar = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "HiddenValley:qqbar2CvCvbar = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "off")
{
$data = "HiddenValley:qqbar2BvBvbar = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "HiddenValley:qqbar2TvTvbar = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "HiddenValley:ffbar2DvDvbar = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "off")
{
$data = "HiddenValley:ffbar2UvUvbar = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "off")
{
$data = "HiddenValley:ffbar2SvSvbar = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "off")
{
$data = "HiddenValley:ffbar2CvCvbar = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "off")
{
$data = "HiddenValley:ffbar2BvBvbar = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "off")
{
$data = "HiddenValley:ffbar2TvTvbar = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "off")
{
$data = "HiddenValley:ffbar2EvEvbar = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "off")
{
$data = "HiddenValley:ffbar2nuEvnuEvbar = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "off")
{
$data = "HiddenValley:ffbar2MUvMUvbar = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "off")
{
$data = "HiddenValley:ffbar2nuMUvnuMUvbar = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "off")
{
$data = "HiddenValley:ffbar2TAUvTAUvbar = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "off")
{
$data = "HiddenValley:ffbar2nuTAUvnuTAUvbar = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "off")
{
$data = "HiddenValley:ffbar2Zv = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "off")
{
$data = "HiddenValley:FSR = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "0.1")
{
$data = "HiddenValley:alphaFSR = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "0.4")
{
$data = "HiddenValley:pTminFSR = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "off")
{
$data = "HiddenValley:fragment = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "1")
{
$data = "HiddenValley:nFlav = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "0.75")
{
$data = "HiddenValley:probVector = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "0.3")
{
$data = "HiddenValley:aLund = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
if($_POST["40"] != "0.8")
{
$data = "HiddenValley:bmqv2 = ".$_POST["40"]."\n";
fwrite($handle,$data);
}
if($_POST["41"] != "1.0")
{
$data = "HiddenValley:rFactqv = ".$_POST["41"]."\n";
fwrite($handle,$data);
}
if($_POST["42"] != "0.5")
{
$data = "HiddenValley:sigmamqv = ".$_POST["42"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
