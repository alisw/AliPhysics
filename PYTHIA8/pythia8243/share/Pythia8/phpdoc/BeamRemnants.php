<html>
<head>
<title>Beam Remnants</title>
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

<form method='post' action='BeamRemnants.php'>
 
<h2>Beam Remnants</h2> 
<ol id="toc">
  <li><a href="#section0">Introduction</a></li>
  <li><a href="#section1">Primordial <ei>kT</ei></a></li>
  <li><a href="#section2">Colour flow</a></li>
  <li><a href="#section3">Further variables</a></li>
</ol>

 
<a name="section0"></a> 
<h3>Introduction</h3> 
 
The <code>BeamParticle</code> class contains information on all partons 
extracted from a beam (so far). As each consecutive multiparton interaction 
defines its respective incoming parton to the hard scattering a 
new slot is added to the list. This information is modified when 
the backwards evolution of the spacelike shower defines a new 
initiator parton. It is used, both for the multiparton interactions 
and the spacelike showers, to define rescaled parton densities based 
on the <i>x</i> and flavours already extracted, and to distinguish 
between valence, sea and companion quarks. Once the perturbative 
evolution is finished, further beam remnants are added to obtain a 
consistent set of flavours. The current physics framework is further 
described in [<a href="Bibliography.php#refSjo04" target="page">Sjo04</a>]. 
 
<p/> 
The introduction of <?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?>rescattering</a> 
in the multiparton interactions framework further complicates the 
processing of events. Specifically, when combined with showers, 
the momentum of an individual parton is no longer uniquely associated 
with one single subcollision. Nevertheless the parton is classified 
with one system, owing to the technical and administrative complications 
of more complete classifications. Therefore the addition of primordial 
<i>kT</i> to the subsystem initiator partons does not automatically 
guarantee overall <i>pT</i> conservation. Various tricks are used to 
minimize the mismatch, with a brute force shift of all parton 
<i>pT</i>'s as a final step. 
 
<p/> 
Much of the above information is stored in a vector of 
<code>ResolvedParton</code> objects, which each contains flavour and 
momentum information, as well as valence/companion information and more. 
The <code>BeamParticle</code> method <code>list()</code> shows the 
contents of this vector, mainly for debug purposes. 
 
<p/> 
The <code>BeamRemnants</code> class takes over for the final step 
of adding primordial <i>kT</i> to the initiators and remnants, 
assigning the relative longitudinal momentum sharing among the 
remnants, and constructing the overall kinematics and colour flow. 
This step couples the two sides of an event, and could therefore 
not be covered in the <code>BeamParticle</code> class, which only 
considers one beam at a time. 
 
<p/> 
The methods of these classes are not intended for general use, 
and so are not described here. 
 
<p/> 
In addition to the parameters described on this page, note that the 
choice of <?php $filepath = $_GET["filepath"];
echo "<a href='PDFSelection.php?filepath=".$filepath."' target='page'>";?>parton densities</a> is made 
in the <code>Pythia</code> class. Then pointers to the pdf's are handed 
on to <code>BeamParticle</code> at initialization, for all subsequent 
usage. 
 
<a name="section1"></a> 
<h3>Primordial <i>kT</i></h3> 
 
The primordial <i>kT</i> of initiators of hard-scattering subsystems 
are selected according to Gaussian distributions in <i>p_x</i> and 
<i>p_y</i> separately. The widths of these distributions are chosen 
to be dependent on the hard scale of the central process and on the mass 
of the whole subsystem defined by the two initiators: 
<br/><i> 
sigma = (sigma_soft * Q_half + sigma_hard * Q) / (Q_half + Q) 
  * m / (m + m_half * y_damp) 
</i><br/> 
Here <i>Q</i> is the hard-process renormalization scale for the 
hardest process and the <i>pT</i> scale for subsequent multiparton 
interactions, <i>m</i> the mass of the system, and 
<i>sigma_soft</i>, <i>sigma_hard</i>, <i>Q_half</i>, 
<i>m_half</i> and <i>y_damp</i> parameters defined below. 
Furthermore each separately defined beam remnant has a distribution 
of width <i>sigma_remn</i>, independently of kinematical variables. 
 
<br/><br/><strong>BeamRemnants:primordialKT</strong>  <input type="radio" name="1" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="1" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow or not selection of primordial <i>kT</i> according to the 
parameter values below. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:primordialKTsoft </td><td></td><td> <input type="text" name="2" value="0.9" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.9</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width <i>sigma_soft</i> in the above equation, assigned as a 
primordial <i>kT</i> to initiators in the soft-interaction limit. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:primordialKThard </td><td></td><td> <input type="text" name="3" value="1.8" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.8</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width <i>sigma_hard</i> in the above equation, assigned as a 
primordial <i>kT</i> to initiators in the hard-interaction limit. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:halfScaleForKT </td><td></td><td> <input type="text" name="4" value="1.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.5</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The scale <i>Q_half</i> in the equation above, defining the 
half-way point between hard and soft interactions. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:halfMassForKT </td><td></td><td> <input type="text" name="5" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The scale <i>m_half</i> in the equation above, defining the 
half-way point between low-mass and high-mass subsystems. 
(Kinematics construction can easily fail if a system is assigned 
a primordial <i>kT</i> value higher than its mass, so the 
mass-dampening is intended to reduce some troubles later on.) 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:reducedKTatHighY </td><td></td><td> <input type="text" name="6" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
For a system of mass <i>m</i> and energy <i>E</i> the 
dampening factor <i>y_damp</i> above is defined as 
<i>y_damp = pow( E/m, r_red)</i>, where <i>r_red</i> is the 
current parameter. The effect is to reduce the primordial <i>kT</i> 
of low-mass systems extra much if they are at large rapidities (recall 
that <i>E/m = cosh(y)</i> before <i>kT</i> is added). The reason 
for this dampening is purely technical, and for reasonable values 
should not have dramatic consequences overall. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:primordialKTremnant </td><td></td><td> <input type="text" name="7" value="0.4" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.4</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width <i>sigma_remn</i>, assigned as a primordial <i>kT</i> 
to beam-remnant partons. 
   
 
<p/> 
A net <i>kT</i> imbalance is obtained from the vector sum of the 
primordial <i>kT</i> values of all initiators and all beam remnants. 
This quantity is compensated by a shift shared equally between 
all partons, except that the dampening factor <i>m / (m_half + m)</i> 
is again used to suppress the role of small-mass systems. 
 
<p/> 
Note that the current <i>sigma</i> definition implies that 
<i>&lt;pT^2&gt; = &lt;p_x^2&gt;+ &lt;p_y^2&gt; = 2 sigma^2</i>. 
It thus cannot be compared directly with the <i>sigma</i> 
of nonperturbative hadronization, where each quark-antiquark 
breakup corresponds to <i>&lt;pT^2&gt; = sigma^2</i> and only 
for hadrons it holds that <i>&lt;pT^2&gt; = 2 sigma^2</i>. 
The comparison is further complicated by the reduction of 
primordial <i>kT</i> values by the overall compensation mechanism. 
 
<br/><br/><strong>BeamRemnants:rescatterRestoreY</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Is only relevant when <?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?>rescattering</a> 
is switched on in the multiparton interactions scenario. For a normal 
interaction the rapidity and mass of a system is preserved when 
primordial <i>kT</i> is introduced, by appropriate modification of the 
incoming parton momenta. Kinematics construction is more complicated for 
a rescattering, and two options are offered. Differences between these 
can be used to explore systematic uncertainties in the rescattering 
framework.<br/> 
The default behaviour is to keep the incoming rescattered parton as is, 
but to modify the unrescattered incoming parton so as to preserve the 
invariant mass of the system. Thereby the rapidity of the rescattering 
is modified.<br/> 
The alternative is to retain the rapidity (and mass) of the rescattered 
system when primordial <i>kT</i> is introduced. This is made at the 
expense of a modified longitudinal momentum of the incoming rescattered 
parton, so that it does not agree with the momentum it ought to have had 
by the kinematics of the previous interaction.<br/> 
For a double rescattering, when both incoming partons have already scattered, 
there is no obvious way to retain the invariant mass of the system in the 
first approach, so the second is always used. 
   
 
<a name="section2"></a> 
<h3>Colour flow</h3> 
 
The colour in the separate subproccsses are tied together via the assignment 
of colour flow in the beam remnants. The assignment of colour flow is not 
known from first principles and therefore it is not an unambiguous procedure. 
Thus two different models have been implemented in <code>Pythia</code>. These 
will be referred to as new and old, based on the time of the implementation. 
 
<p/> 
The old model tries to reconstruct the colour flow in a way that a LO PS would 
produce the beam remnants. The starting point is the junction structure of the 
beam particle (if it is a baryon). The gluons are attached to a quark line and 
quark-antiquark pairs are added as if coming from a gluon splittings. Thus 
this model captures the qualitative behaviour that is expected from leading 
colour QCD. The model is described in  detail in [<a href="Bibliography.php#refSjo04" target="page">Sjo04</a>]. 
 
<p/> 
The new model is built on the full SU(3) colour structure of QCD. The 
starting point is the scattered partons from the MPI. Each of these are 
initially assumed uncorrelated in colour space, allowing the total outgoing 
colour configuration to be calculated as an SU(3) product. Since the beam 
particle is a colour singlet, the beam remnant colour configuration has to be 
the inverse of the outgoing colour configuration. The minimum amount of gluons 
are added to the beam remnant in order to obtain this colour configuration. 
 
<p/> 
The above assumption of uncorrelated MPIs in colour space is a good 
assumption for a few well separated hard MPIs. However if the number of MPIs 
become large and ISR is included, such that the energy scale becomes lower 
(and thus distances becomes larger), the assumption loses its validity. This 
is due to saturation effects. The modelling of saturation is done in crude 
manner, as an exponential suppresion of high multiplet states. 
 
<p/> 
None of the models above can provide a full description of the colour 
flow in an event, however. Therefore additional colour reconfiguration 
is needed. This is referred to as colour reconnection. Several different 
models for colour reconnection are implemented, see 
<?php $filepath = $_GET["filepath"];
echo "<a href='ColourReconnection.php?filepath=".$filepath."' target='page'>";?>Colour Reconection</a>. 
 
<br/><br/><table><tr><td><strong>BeamRemnants:remnantMode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
Switch to choose between the two different colour models for the beam remnant. 
<br/>
<input type="radio" name="9" value="0" checked="checked"><strong>0 </strong>:  The old beam remnant model. <br/>
<input type="radio" name="9" value="1"><strong>1 </strong>:  The new beam remnant model. <br/>
 
<br/><br/><table><tr><td><strong>BeamRemnants:saturation </td><td></td><td> <input type="text" name="10" value="5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 100000</code>)</td></tr></table>
Controls the suppresion due to saturation in the new model. The exact formula 
used is <i>exp(-M / k)</i>, where M is the multiplet size and k is this 
parameter. Thus a small number will result in a large saturation. 
   
 
<a name="section3"></a> 
<h3>Further variables</h3> 
 
<br/><br/><table><tr><td><strong>BeamRemnants:maxValQuark  </td><td></td><td> <input type="text" name="11" value="3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
The maximum valence quark kind allowed in acceptable incoming beams, 
for which multiparton interactions are simulated. Default is that hadrons 
may contain <i>u</i>, <i>d</i> and <i>s</i> quarks, 
but not <i>c</i> and <i>b</i> ones, since sensible 
kinematics has not really been worked out for the latter. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:companionPower  </td><td></td><td> <input type="text" name="12" value="4" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4</strong></code>; <code>minimum = 0</code>; <code>maximum = 4</code>)</td></tr></table>
When a sea quark has been found, a companion antisea quark ought to be 
nearby in <i>x</i>. The shape of this distribution can be derived 
from the gluon mother distribution convoluted with the 
<i>g &rarr; q qbar</i> splitting kernel. In practice, simple solutions 
are only feasible if the gluon shape is assumed to be of the form 
<i>g(x) ~ (1 - x)^p / x</i>, where <i>p</i> is an integer power, 
the parameter above. Allowed values correspond to the cases programmed. 
<br/> 
Since the whole framework is approximate anyway, this should be good 
enough. Note that companions typically are found at small <i>Q^2</i>, 
if at all, so the form is supposed to represent <i>g(x)</i> at small 
<i>Q^2</i> scales, close to the lower cutoff for multiparton interactions. 
   
 
<p/> 
When assigning relative momentum fractions to beam-remnant partons, 
valence quarks are chosen according to a distribution like 
<i>(1 - x)^power / sqrt(x)</i>. This <i>power</i> is given below 
for quarks in mesons, and separately for <i>u</i> and <i>d</i> 
quarks in the proton, based on the approximate shape of low-<i>Q^2</i> 
parton densities. The power for other baryons is derived from the 
proton ones, by an appropriate mixing. The <i>x</i> of a diquark 
is chosen as the sum of its two constituent <i>x</i> values, and can 
thus be above unity. (A common rescaling of all remnant partons and 
particles will fix that.) An additional enhancement of the diquark 
momentum is obtained by its <i>x</i> value being rescaled by the 
<code>valenceDiqEnhance</code> factor. 
 
<br/><br/><table><tr><td><strong>BeamRemnants:valencePowerMeson </td><td></td><td> <input type="text" name="13" value="0.8" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.8</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The abovementioned power for valence quarks in mesons. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:valencePowerUinP </td><td></td><td> <input type="text" name="14" value="3.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3.5</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The abovementioned power for valence <i>u</i> quarks in protons. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:valencePowerDinP </td><td></td><td> <input type="text" name="15" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The abovementioned power for valence <i>d</i> quarks in protons. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:valenceDiqEnhance </td><td></td><td> <input type="text" name="16" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.5</code>; <code>maximum = 10.</code>)</td></tr></table>
Enhancement factor for valence diquarks in baryons, relative to the 
simple sum of the two constituent quarks. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:gluonPower </td><td></td><td> <input type="text" name="17" value="4.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4.0</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The abovementioned power for gluons. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:xGluonCutoff </td><td></td><td> <input type="text" name="18" value="1E-7" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1E-7</strong></code>; <code>minimum = 1E-10</code>; <code>maximum = 1</code>)</td></tr></table>
The gluon PDF is approximated with <i>g(x) ~ (1 - x)^p / x</i>, which 
integrates to infinity when integrated from 0 to 1. This cut-off is 
introduced as a minimum to avoid the problems with infinities. 
   
 
<br/><br/><strong>BeamRemnants:allowJunction</strong>  <input type="radio" name="19" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="19" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
The <code>off</code> option is intended for debug purposes only, as 
follows. When more than one valence quark is kicked out of a baryon 
beam, as part of the multiparton interactions scenario, the subsequent 
hadronization is described in terms of a junction string topology. 
This description involves a number of technical complications that 
may make the program more unstable. As an alternative, by switching 
this option off, junction configurations are rejected (which gives 
an error message that the remnant flavour setup failed), and the 
multiparton interactions and showers are redone until a 
junction-free topology is found. 
   
 
<br/><br/><strong>BeamRemnants:beamJunction</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
This parameter is only relevant if the new colour reconnection scheme is used. 
(see  <?php $filepath = $_GET["filepath"];
echo "<a href='ColourReconnection.php?filepath=".$filepath."' target='page'>";?>colour reconnection</a>) 
This parameter tells whether to form a junction or a di-quark if more 
than two valence quarks are found in the beam remnants. If off a di-quark is 
formed and if on a junction will be formed. 
   
 
<br/><br/><strong>BeamRemnants:allowBeamJunction</strong>  <input type="radio" name="21" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="21" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
This parameter is only relevant if the new Beam remnant model is used. 
This parameter tells whether to allow the formation of junction structures 
in the colour configuration of the scattered partons. 
   
 
<br/><br/><table><tr><td><strong>BeamRemnants:unresolvedHadron  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
Switch to to force either or both of the beam remnants to collapse to a 
single hadron, namely the original incoming one. Must only be used when this 
is physically meaningful, e.g. when a photon can be viewed as emitted from 
a proton that does not break up in the process. 
<br/>
<input type="radio" name="22" value="0" checked="checked"><strong>0 </strong>:  Both hadronic beams are resolved. <br/>
<input type="radio" name="22" value="1"><strong>1 </strong>:  Beam A is unresolved, beam B resolved. <br/>
<input type="radio" name="22" value="2"><strong>2 </strong>:  Beam A is resolved, beam B unresolved. <br/>
<input type="radio" name="22" value="3"><strong>3 </strong>:  Both hadronic beams are unresolved. <br/>
 
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
$data = "BeamRemnants:primordialKT = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0.9")
{
$data = "BeamRemnants:primordialKTsoft = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "1.8")
{
$data = "BeamRemnants:primordialKThard = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "1.5")
{
$data = "BeamRemnants:halfScaleForKT = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1.")
{
$data = "BeamRemnants:halfMassForKT = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.5")
{
$data = "BeamRemnants:reducedKTatHighY = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "0.4")
{
$data = "BeamRemnants:primordialKTremnant = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "BeamRemnants:rescatterRestoreY = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "0")
{
$data = "BeamRemnants:remnantMode = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "5")
{
$data = "BeamRemnants:saturation = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "3")
{
$data = "BeamRemnants:maxValQuark = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "4")
{
$data = "BeamRemnants:companionPower = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0.8")
{
$data = "BeamRemnants:valencePowerMeson = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "3.5")
{
$data = "BeamRemnants:valencePowerUinP = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "2.0")
{
$data = "BeamRemnants:valencePowerDinP = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "2.0")
{
$data = "BeamRemnants:valenceDiqEnhance = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "4.0")
{
$data = "BeamRemnants:gluonPower = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "1E-7")
{
$data = "BeamRemnants:xGluonCutoff = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "on")
{
$data = "BeamRemnants:allowJunction = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "BeamRemnants:beamJunction = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "on")
{
$data = "BeamRemnants:allowBeamJunction = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "0")
{
$data = "BeamRemnants:unresolvedHadron = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
