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
described in [<a href="Bibliography.php" target="page">Sjo04</a>]. 

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

<h3>Primordial <i>kT</i></h3>

The primordial <i>kT</i> of initiators of hard-scattering subsystems 
are selected according to Gaussian distributions in <i>p_x</i> and 
<i>p_y</i> separately. The widths of these distributions are chosen
to be dependent on the hard scale of the central process and on the mass 
of the whole subsystem defined by the two initiators:
<br/><i>
sigma = (sigma_soft * Q_half + sigma_hard * Q) / (Q_half + Q) 
  * m / (m_half + m)  
</i><br/>
Here <i>Q</i> is the hard-process renormalization scale for the 
hardest process and the <i>pT</i> scale for subsequent multiparton
interactions, <i>m</i> the mass of the system, and 
<i>sigma_soft</i>, <i>sigma_hard</i>, <i>Q_half</i> and
<i>m_half</i> parameters defined below. Furthermore each separately
defined beam remnant has a distribution of width <i>sigma_remn</i>, 
independently of kinematical variables.

<br/><br/><strong>BeamRemnants:primordialKT</strong>  <input type="radio" name="1" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="1" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow or not selection of primordial <i>kT</i> according to the
parameter values below.
  

<br/><br/><table><tr><td><strong>BeamRemnants:primordialKTsoft </td><td></td><td> <input type="text" name="2" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width <i>sigma_soft</i> in the above equation, assigned as a 
primordial <i>kT</i> to initiators in the soft-interaction limit.
  

<br/><br/><table><tr><td><strong>BeamRemnants:primordialKThard </td><td></td><td> <input type="text" name="3" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The width <i>sigma_hard</i> in the above equation, assigned as a 
primordial <i>kT</i> to initiators in the hard-interaction limit.
  

<br/><br/><table><tr><td><strong>BeamRemnants:halfScaleForKT </td><td></td><td> <input type="text" name="4" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The scale <i>Q_half</i> in the equation above, defining the 
half-way point between hard and soft interactions. 
  

<br/><br/><table><tr><td><strong>BeamRemnants:halfMassForKT </td><td></td><td> <input type="text" name="5" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The scale <i>m_half</i> in the equation above, defining the 
half-way point between low-mass and high-mass subsystems.
(Kinematics construction can easily fail if a system is assigned 
a primordial <i>kT</i> value higher than its mass, so the 
mass-dampening is intended to reduce some troubles later on.)
  

<br/><br/><table><tr><td><strong>BeamRemnants:primordialKTremnant </td><td></td><td> <input type="text" name="6" value="0.4" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.4</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
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

<br/><br/><strong>BeamRemnants:rescatterRestoreY</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
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
  

<h3>Colour flow</h3>

The colour flows in the separate subprocesses defined in the 
multiparton-interactions scenario are tied together via the assignment
of colour flow in the beam remnant. This is not an unambiguous 
procedure, but currently no parameters are directly associated with it.
However, a simple "minimal" procedure of colour flow only via the beam 
remnants does not result in a scenario in
agreement with data, notably not a sufficiently steep rise of  
<i>&lt;pT&gt;(n_ch)</i>. The true origin of this behaviour and the
correct mechanism to reproduce it remains one of the big unsolved issues 
at the borderline between perturbative and nonperturbative QCD. 
As a simple attempt, an additional step is introduced, wherein the gluons 
of a lower-<i>pT</i> system are merged with the ones in a higher-pT one. 

<br/><br/><strong>BeamRemnants:reconnectColours</strong>  <input type="radio" name="8" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="8" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow or not a system to be merged with another one.
  

<br/><br/><table><tr><td><strong>BeamRemnants:reconnectRange </td><td></td><td> <input type="text" name="9" value="10.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
A system with a hard scale <i>pT</i> can be merged with one of a 
harder scale with a probability that is 
<i>pT0_Rec^2 / (pT0_Rec^2 + pT^2)</i>, where
<i>pT0_Rec</i> is <code>reconnectRange</code> times <i>pT0</i>, 
the latter being the same energy-dependent dampening parameter as 
used for multiparton interactions. 
Thus it is easy to merge a low-<i>pT</i> system with any other,
but difficult to merge two high-<i>pT</i> ones with each other. 
  

<p/>
The procedure is used iteratively. Thus first the reconnection probability
<i>P = pT0_Rec^2 / (pT0_Rec^2 + pT^2)</i> of the lowest-<i>pT</i> 
system is found, and gives the probability for merger with the 
second-lowest one. If not merged, it is tested with the third-lowest one, 
and so on. For the <i>m</i>'th higher system the reconnection
probability thus becomes <i>(1 - P)^(m-1) P</i>. That is, there is 
no explicit dependence on the higher <i>pT</i> scale, but implicitly 
there is via the survival probability of not already having been merged
with a lower-<i>pT</i> system. Also note that the total reconnection
probability for the lowest-<i>pT</i> system in an event with <i>n</i> 
systems becomes <i>1 - (1 - P)^(n-1)</i>. Once the fate of the 
lowest-<i>pT</i> system has been decided, the second-lowest is considered
with respect to the ones above it, then the third-lowest, and so on.  

<p/>
Once it has been decided which systems should be joined, the actual merging
is carried out in the opposite direction. That is, first the hardest
system is studied, and all colour dipoles in it are found (including to 
the beam remnants, as defined by the holes of the incoming partons).
Next each softer system to be merged is studied in turn. Its gluons are,
in decreasing <i>pT</i> order, inserted on the colour dipole <i>i,j</i>
that gives the smallest <i>(p_g p_i)(p_g p_j)/(p_i p_j)</i>, i.e. 
minimizes the "disturbance" on the existing dipole, in terms of 
<i>pT^2</i> or <i>Lambda</i> measure (string length). The insertion
of the gluon means that the old dipole is replaced by two new ones. 
Also the (rather few) quark-antiquark pairs that can be traced back to 
a gluon splitting are treated in close analogy with the gluon case. 
Quark lines that attach directly to the beam remnants cannot be merged 
but are left behind. 

<p/>
The joining procedure can be viewed as a more sophisticated variant of 
the one introduced already in [<a href="Bibliography.php" target="page">Sjo87</a>]. Clearly it is ad hoc. 
It hopefully captures some elements of truth. The lower <i>pT</i> scale 
a system has the larger its spatial extent and therefore the larger its 
overlap with other systems. It could be argued that one should classify 
individual initial-state partons by <i>pT</i> rather than the system 
as a whole. However, for final-state radiation, a soft gluon radiated off 
a hard parton is actually produced at late times and therefore probably 
less likely to reconnect. In the balance, a classification by system 
<i>pT</i> scale appears sensible as a first try. 

<p/>
Note that the reconnection is carried out before resonance decays are
considered. Colour inside a resonance therefore is not reconnected.
This is a deliberate choice, but certainly open to discussion and 
extensions at a later stage, as is the rest of this procedure.

<h3>Further variables</h3>

<br/><br/><table><tr><td><strong>BeamRemnants:maxValQuark  </td><td></td><td> <input type="text" name="10" value="3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
The maximum valence quark kind allowed in acceptable incoming beams,
for which multiparton interactions are simulated. Default is that hadrons
may contain <i>u</i>, <i>d</i> and <i>s</i> quarks, 
but not <i>c</i> and <i>b</i> ones, since sensible
kinematics has not really been worked out for the latter.
  

<br/><br/><table><tr><td><strong>BeamRemnants:companionPower  </td><td></td><td> <input type="text" name="11" value="4" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4</strong></code>; <code>minimum = 0</code>; <code>maximum = 4</code>)</td></tr></table>
When a sea quark has been found, a companion antisea quark ought to be
nearby in <i>x</i>. The shape of this distribution can be derived 
from the gluon mother distribution convoluted with the 
<i>g -> q qbar</i> splitting kernel. In practice, simple solutions 
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

<br/><br/><table><tr><td><strong>BeamRemnants:valencePowerMeson </td><td></td><td> <input type="text" name="12" value="0.8" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.8</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The abovementioned power for valence quarks in mesons.
  

<br/><br/><table><tr><td><strong>BeamRemnants:valencePowerUinP </td><td></td><td> <input type="text" name="13" value="3.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3.5</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The abovementioned power for valence <i>u</i> quarks in protons.
  

<br/><br/><table><tr><td><strong>BeamRemnants:valencePowerDinP </td><td></td><td> <input type="text" name="14" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The abovementioned power for valence <i>d</i> quarks in protons.
  

<br/><br/><table><tr><td><strong>BeamRemnants:valenceDiqEnhance </td><td></td><td> <input type="text" name="15" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.5</code>; <code>maximum = 10.</code>)</td></tr></table>
Enhancement factor for valence diquarks in baryons, relative to the 
simple sum of the two constituent quarks.
  

<br/><br/><strong>BeamRemnants:allowJunction</strong>  <input type="radio" name="16" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="16" value="off"><strong>Off</strong>
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
if($_POST["2"] != "0.5")
{
$data = "BeamRemnants:primordialKTsoft = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "2.0")
{
$data = "BeamRemnants:primordialKThard = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "1.")
{
$data = "BeamRemnants:halfScaleForKT = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1.")
{
$data = "BeamRemnants:halfMassForKT = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.4")
{
$data = "BeamRemnants:primordialKTremnant = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "BeamRemnants:rescatterRestoreY = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "on")
{
$data = "BeamRemnants:reconnectColours = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "10.0")
{
$data = "BeamRemnants:reconnectRange = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "3")
{
$data = "BeamRemnants:maxValQuark = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "4")
{
$data = "BeamRemnants:companionPower = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "0.8")
{
$data = "BeamRemnants:valencePowerMeson = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "3.5")
{
$data = "BeamRemnants:valencePowerUinP = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "2.0")
{
$data = "BeamRemnants:valencePowerDinP = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "2.0")
{
$data = "BeamRemnants:valenceDiqEnhance = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "on")
{
$data = "BeamRemnants:allowJunction = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
