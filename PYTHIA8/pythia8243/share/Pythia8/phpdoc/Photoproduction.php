<html>
<head>
<title>Photoproduction</title>
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

<form method='post' action='Photoproduction.php'>
 
<h2>Photoproduction</h2> 
<ol id="toc">
  <li><a href="#section0">Types of photon processes</a></li>
  <li><a href="#section1">Resolved photon</a></li>
  <li><a href="#section2">Photons from lepton beams</a></li>
</ol>

 
<p> 
Interactions involving one or two photons, either in photon-photon or 
photon-hadron collision or photons emitted from lepton beams. 
Includes both direct and resolved contributions and also soft QCD and MPIs 
for events with resolved photons. Only (quasi-)real photons are considered 
so virtuality of the photons is restricted. The PDF set for resolved photons 
is selected in the <?php $filepath = $_GET["filepath"];
echo "<a href='PDFSelection.php?filepath=".$filepath."' target='page'>";?>PDF selection</a>. 
This page describes some of the special features related to these collisions 
and introduces the relevant parameters. 
</p> 
 
<a name="section0"></a> 
<h3>Types of photon processes</h3> 
 
<p> 
Photons can be either resolved or act as point-like particles (direct). 
Therefore for a photon-photon interaction there are four different 
contributions, resolved-resolved, resolved-direct, direct-resolved and 
direct-direct. In case of photon-hadron collisions there are two 
contributions. With the default value of the parameter below, a mix of 
relevant contributions is generated but each process type can also be 
generated individually. Note that for photon-hadron collisions the code 
for direct contribution depends on which of the beams is photon. 
The sample main program <code>main69.cc</code> demonstrates different 
possibilities. 
</p> 
 
<br/><br/><table><tr><td><strong>Photon:ProcessType  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 4</code>)</td></tr></table>
Sets desired contribution for interactions with one or two photons. 
<br/>
<input type="radio" name="1" value="0" checked="checked"><strong>0 </strong>:  Mix of relevant contributions below. <br/>
<input type="radio" name="1" value="1"><strong>1 </strong>:  Resolved-Resolved: Both colliding photons are  resolved and the partonic content is given by the PDFs. Hard processes  and non-diffractive events can be generated. <br/>
<input type="radio" name="1" value="2"><strong>2 </strong>:  Resolved-Direct: Photon A is resolved and photon B  unresolved, i.e. act as an initiator for the hard process. Hard processes  with a parton and a photon in the initial state can be generated.  In case of photon-hadron collision this provides the direct contribution  when hadron is beam A and photon beam B.<br/>
<input type="radio" name="1" value="3"><strong>3 </strong>:  Direct-Resolved: As above but now photon A is unresolved  and photon B resolved. Direct contribution of photon-hadron when photon  beam A.<br/>
<input type="radio" name="1" value="4"><strong>4 </strong>:  Direct-Direct: Both photons are unresolved. Hard  processes with two photon initiators can be generated.<br/>
 
<p> 
The type of the generated process can be obtained from 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>Info class</a> with method 
<code>int Info::photonMode()</code> which follows the conventions above. 
</p> 
 
<a name="section1"></a> 
<h3>Resolved photon</h3> 
 
<p> 
Photons can either interact directly as an unresolved particle or as a 
hadronic state ("Vector Meson Dominance"). In the latter case the hard 
process can be simulated using PDFs to describe the partonic structure 
of the resolved photon. The evolution equations for photons include an 
additional term that corresponds to <i>gamma &rarr; q qbar</i> splittings. 
Due to this, the PDFs are somewhat different for photons than for hadrons 
and some parts of event generation need special attention. 
</p> 
 
<h4>Process-level generation</h4> 
 
<p> 
Due to the additional term in the evolution equations the quarks in a 
resolved photon may carry a very large fraction <i>(x~1)</i> of the photon 
momentum. In these cases it may happen that, after the hard process, there is 
no energy left to construct the beam remnants. This is true especially if 
a heavy quark is taken out from the beam and a corresponding massive 
antiquark needs to be added to the remnant system to conserve flavour. Even 
though these events are allowed based on the PDFs alone, they are not physical 
and should be rejected. Therefore some amount of errors can be expected when 
generating events close to the edge of phase space, e.g. when collision 
energy is low. 
</p> 
 
<h4>Spacelike showers</h4> 
 
<p> 
The parton showers are generated according to the DGLAP evolution equations. 
Due to the <i>gamma &rarr; q qbar</i> splitting in the photon evolution, 
a few modifications are needed for the ISR algorithm. 
<ul> 
<li> 
The additional term corresponds to a possibility to find the original beam 
photon during the backwards evolution, which is added to the QED part of the 
spacelike shower evolution. If this splitting happens there is no need to 
construct the beam remnants for the given beam. 
</li> 
<li> 
The heavy quark production threshold with photon beams is handled in a 
similar manner as for hadrons, but now the splitting that is forced 
to happen is <i>gamma &rarr; Q Qbar</i>. 
</li> 
<li> 
As the splittings in backwards evolution increases the <i>x</i> of the 
parton taken from the beam, the splittings can lead to a situation where 
there is no room left for massive beam remnants. To make sure that the 
required  remnants can be constructed, splittings that would not leave 
room for the beam remnants are not allowed. 
</li> 
</ul> 
</p> 
 
<h4>MPIs with photon beams</h4> 
 
<p> 
Multiparton interactions with resolved photon beams are generated as with 
hadron beams. The only difference follows again from the additional 
<i>gamma &rarr; q qbar</i> splittings where the beam photon becomes 
unresolved. If this splitting happens during the interleaved evolution 
for either of the photon beams no further MPIs below the branching scale 
<i>pT</i> are allowed since the photon is not resolved any. 
</p> 
 
<p> 
If there have been multiple interactions and a <i>gamma &rarr; q qbar 
</i> splitting occur, the kinematics of this branching are not constructed 
in the spacelike shower. Instead the <i>pT</i> scale of the branching is 
stored and the relevant momenta are then fixed in the beam remnant handling. 
Therefore the status codes for the partons related to this splitting 
actually refer to beam remnants. 
</p> 
 
<p> 
If there are no MPIs before the <i>gamma &rarr; q qbar</i> splitting, 
this splitting is constructed in the spacelike shower in the usual way, 
but the mother beam photon is not added to the event record, since a copy 
of it already exists at the top of the event record. This is unlike the 
documentation of other ISR splittings, where the mother of the branching 
is shown, but consistent with the photon not being added (a second time) 
for events that contain several MPIs. Optionally the photon can be shown, 
using the following flag. 
 
<br/><br/><strong>Photon:showUnres</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Show the evolution steps of the beam photon in the event record, if on. 
   
</p> 
 
<p> 
Based on comparisons with charged hadron production in photon-photon 
collision data from LEP, the default MPI parametrization tuned to 
proton-(anti)proton collisions produces too much hadrons from the 
additional interactions. Such differences are not surprising, given 
that the photon is less hadron-like than the proton, e.g. with less 
well developed PDFs, leaving less room for MPIs. Therefore a different 
parametrization for <i>pT0(eCM)</i> is used in case of photon-photon 
collisions, where the default values are tuned to the LEP data 
(a reference to this study will be added later). By default, 
a logarithmic dependence on <i>eCM</i> is used. 
<br/><b>Note:</b> These parameters override the choices made in 
<?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?>Multiparton Interactions</a> 
when photon-photon collisions are generated. 
</p> 
 
<br/><br/><table><tr><td><strong>PhotonPhoton:pT0parametrization  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
Choice of <ei>pT0</ei> parametrization. See 
<aloc href="MultipartonInteractions">Multiparton Interactions</aloc> for 
further details. 
<br/>
<input type="radio" name="3" value="0"><strong>0 </strong>: Power law in <ei>eCM</ei>.<br/>
<input type="radio" name="3" value="1" checked="checked"><strong>1 </strong>: Logarithmic in <ei>eCM</ei>.<br/>
 
<br/><br/><table><tr><td><strong>PhotonPhoton:ecmRef </td><td></td><td> <input type="text" name="4" value="100.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100.0</strong></code>; <code>minimum = 1.</code>)</td></tr></table>
The <i>ecmRef</i> reference energy scale. 
   
 
<br/><br/><table><tr><td><strong>PhotonPhoton:pT0Ref </td><td></td><td> <input type="text" name="5" value="1.52" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.52</strong></code>; <code>minimum = 0.5</code>; <code>maximum = 10.0</code>)</td></tr></table>
The value of <i>pT0</i> at the reference energy scale. 
   
 
<br/><br/><table><tr><td><strong>PhotonPhoton:ecmPow </td><td></td><td> <input type="text" name="6" value="0.413" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.413</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 0.5</code>)</td></tr></table>
The <i>ecmPow</i> energy rescaling pace. 
   
 
<p/> 
Alternatively, or in combination, a sharp cut can be used. 
<br/><br/><table><tr><td><strong>PhotonPhoton:pTmin </td><td></td><td> <input type="text" name="7" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.0</code>)</td></tr></table>
More details in 
<?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?>Multiparton Interactions</a>. 
   
 
<p> 
A similar study for photon-hadron collisions will follow, current 
recommendation is to use value <i>pT0Ref = 3.0 GeV</i> set in 
<?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?>Multiparton Interactions</a> page. 
</p> 
 
<p> 
The total cross section for photon-photon collisions is paramerized as in 
[<a href="Bibliography.php#refSch97" target="page">Sch97</a>]. Approximate diffractive cross sections have been defined 
according to the assumed VMD contribution. 
</p> 
 
<h4>Beam Remnants</h4> 
 
<p> 
To construct the beam remnants, one should know whether the parton 
taken from the beam is a valence parton or not. The valence partons of a 
photon includes the partons that originate from <i>gamma &rarr; q qbar</i> 
splittings of the original beam photon and the valence partons from the 
hadron-like part of the PDF. In either case, the flavour of the valence 
quarks can fluctuate. Unfortunately the decomposition to the different 
components are typically not provided in the PDF sets and some further 
assumptions are needed to decide the valence content. 
</p> 
 
<p> 
When ISR is applied for photon beams it is possible to end up to the original 
beam photon during the evolution. Therefore there are three possibilities for 
the remnants: 
<ul> 
<li> 
Remnants need to be constructed for both beams. 
</li> 
<li> 
Remnants are constructed only for one side. 
</li> 
<li> 
No need for remnants on either side. 
</li> 
</ul> 
The last case is the simplest as all the partons in the event are already 
generated by the parton showers. In the first case the remnants and 
primordial <i>kT</i> are constructed similarly as for normal hadronic 
interactions [<a href="Bibliography.php#refSjo04" target="page">Sjo04</a>]. For the second case the momenta of the 
remnant partons can not be balanced between the two beams as the kinematics 
of the other side are already fixed. In these cases the momenta are balanced 
between the scattered system and the remnants. 
</p> 
 
<p> 
Since the primordial <i>kT</i> increases the invariant mass of the remnants 
and the scattered system, it may again happen that there is no room for the 
remnant partons after <i>kT</i> is added, so the kinematics can not be 
constructed. In this case new values for <i>kT</i> are sampled. If this 
does not work, a new shower is generated and in some rare cases the 
parton-level generation fails and the hard process is rejected. The inclusion 
of additional MPIs increases the invariant mass of the remnants and takes 
more momentum from the beam particles. Even though the MPIs that would 
not leave enough room for the remnants are rejected, these can still lead 
to a situation where the kinematics cannot be constructed due to the added 
primordial <i>kT</i>. This may cause some amount of errors especially when 
the invariant mass of <i>gamma-gamma</i> system is small. 
</p> 
 
<a name="section2"></a> 
<h3>Photons from lepton beams</h3> 
 
<p> 
Interaction of photons from leptons including photon-photon interactions in 
lepton-lepton collisions and photon-hadron lepton-hadron collisions can be 
set up as described in 
<?php $filepath = $_GET["filepath"];
echo "<a href='PDFSelection.php?filepath=".$filepath."' target='page'>";?>PDF selection</a>. Since the current 
framework can handle only (quasi-)real photons, a upper limit for the 
photon virtuality needs to be set. This can be done with the parameter 
<code>Photon:Q2max</code>. The upper limit for virtuality will set also 
the upper limit for the <i>k_T</i> of the photon, which in turn will 
be the same as the <i>k_T</i> of the scattered lepton. Also some other 
cuts can be imposed. 
 
<br/><br/><table><tr><td><strong>Photon:Q2max </td><td></td><td> <input type="text" name="8" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>maximum = 2.0</code>)</td></tr></table>
Upper limit for (quasi-)real photon virtuality in <i>GeV^2</i>. Too low 
value might cause problems, e.g. if the lower <i>Q^2</i> limit derived from 
kinematics becomes larger than upper limit. 
   
 
<br/><br/><table><tr><td><strong>Photon:Wmin </td><td></td><td> <input type="text" name="9" value="10.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10.0</strong></code>; <code>minimum = 5.0</code>)</td></tr></table>
Lower limit for invariant mass of <i>gamma-gamma</i> system in <i>GeV</i>. 
In lepton-hadron collisions <i>W</i> corresponds to invariant mass of 
photon-hadron system. 
   
 
<br/><br/><table><tr><td><strong>Photon:Wmax </td><td></td><td> <input type="text" name="10" value="-1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.0</strong></code>)</td></tr></table>
Upper limit for invariant mass of <i>gamma-gamma</i> 
(<i>gamma-hadron</i>) system in <i>GeV</i>. 
A value below <code>Photon:Wmin</code> means that the invariant mass of 
the original <i>l+l-</i> (<i>lepton-hadron</i>) system is used as an 
upper limit. 
   
 
<br/><br/><table><tr><td><strong>Photon:thetaAMax </td><td></td><td> <input type="text" name="11" value="-1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.0</strong></code>; <code>maximum = 3.141593</code>)</td></tr></table>
Upper limit for scattering angle of lepton A in <i>rad</i>. A negative 
value means that no cut is applied. Since <i>k_T</i> depends on virtuality 
of the emitted photon, the <code>Photon:Q2max</code> cut is usually  
restrictive unless a very small angle is used. This cut is only applied 
when the colliding beams are defined in their CM frame 
(<code>Beams:frameType=1</code>). Further, in case of <i>2 &rarr; 1</i> 
processes with direct photons the scattered lepton kinematics is modified 
later in the event generation, so accurate rejection can be obtained only 
based on the final lepton momenta in the event record. 
   
 
<br/><br/><table><tr><td><strong>Photon:thetaBMax </td><td></td><td> <input type="text" name="12" value="-1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.0</strong></code>; <code>maximum = 3.141593</code>)</td></tr></table>
As above but for lepton B. 
   
 
<br/><br/><strong>Photon:sampleQ2</strong>  <input type="radio" name="13" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="13" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Determines whether the sampling for the photon virtuality is done. This can 
be used only with <code>PDF:lepton2gammaSet = 2</code> option where 
<i>Q^2</i>-integrated photon is provided by the user with 
<code>Pythia::setPhotonFluxPtr(PDF*, PDF*)</code> method. In this case the 
virtuality and the transverse momentum of the photon (and recoiling particle) 
is set to zero which strictly speaking is kinematically not possible. The 
error here is, however, very small for the cases where the virtualities 
are negligible (e.g. photons from heavy ions).   
</p> 
 
<h4>MPIs with lepton beams</h4> 
 
<p> 
The invariant mass of <i>gamma-gamma</i> or <i>gamma-hadron</i> system 
from lepton beams will vary. Therefore, to generate MPIs and non-diffractive 
events in <i>gamma-gamma</i> and <i>gamma-hadron</i> 
collisions from lepton beams, the MPI framework is initialized with five 
values of <i>W</i> from <code>Photon:Wmin</code> to 
<code>Photon:Wmax</code>. The parameter values are then interpolated 
for the sampled <i>W</i>. 
</p> 
 
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

if($_POST["1"] != "0")
{
$data = "Photon:ProcessType = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "Photon:showUnres = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "1")
{
$data = "PhotonPhoton:pT0parametrization = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "100.0")
{
$data = "PhotonPhoton:ecmRef = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1.52")
{
$data = "PhotonPhoton:pT0Ref = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.413")
{
$data = "PhotonPhoton:ecmPow = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "0.2")
{
$data = "PhotonPhoton:pTmin = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "1.0")
{
$data = "Photon:Q2max = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "10.0")
{
$data = "Photon:Wmin = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "-1.0")
{
$data = "Photon:Wmax = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "-1.0")
{
$data = "Photon:thetaAMax = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "-1.0")
{
$data = "Photon:thetaBMax = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "on")
{
$data = "Photon:sampleQ2 = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
