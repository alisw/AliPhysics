<html>
<head>
<title>Total Cross Sections</title>
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

<form method='post' action='TotalCrossSections.php'>
 
<h2>Total Cross Sections</h2> 
<ol id="toc">
  <li><a href="#section0">Master switches</a></li>
  <li><a href="#section1">Set your own cross sections</a></li>
  <li><a href="#section2">Modify the SaS/DL cross sections</a></li>
  <li><a href="#section3">Modify the MBR cross sections</a></li>
  <li><a href="#section4">Modify the ABMST cross sections</a></li>
  <li><a href="#section5">Modify the RPP cross sections</a></li>
  <li><a href="#section6">Coulomb corrections to elastic scattering</a></li>
</ol>

 
The <code>SigmaTotal</code> class returns the total, elastic, diffractive 
and nondiffractive cross sections in hadronic collisions. By implication 
it also has to provide differential elastic and diffractive cross sections, 
since many models start out from the differential expressions and then 
integrate to obtain more inclusive rates. In principle it would have been 
possible to decouple the overall normalization from the differential shape, 
however. 
 
<p/> 
The current page describes the options available for integrated and 
differential cross sections alike. The number of options is especially 
large for diffraction, reflecting the lack of a well-understood theory. 
Conversely, the wide spectrum of options should allow for detailed 
comparisons that eventually will improve our understanding. The 
<?php $filepath = $_GET["filepath"];
echo "<a href='Diffraction.php?filepath=".$filepath."' target='page'>";?>Diffraction</a> page contains those further 
parameters needed to describe the hadronization of a diffractive system, 
or at least those that set diffraction apart from the nondiffractive 
topologies. There are borderline cases, that could have been described 
in either place, such as the ones related to the pomeron-proton cross 
section, which mainly are relevant for the description of MPIs in 
diffractive systems, and therefore have been put on the Diffraction page. 
That page also contains the "hard diffraction" framework, i.e. the 
modelling of diffractive events that contain a hard process. 
 
<p/> 
Several different parametrization options are available for <i>p p</i> 
and  <i>pbar p</i> collisions, of special interest for hadron colliders, 
while the selection for other processes is considerably more limited. 
As a simple generalization, a neutron is assumed to have the same hadronic 
cross section as a proton. 
 
<p/> 
Historically most of the parametrizations used are from 
[<a href="Bibliography.php#refSch94" target="page">Sch94</a>, <a href="Bibliography.php#refSch97" target="page">Sch97</a>] which borrows some of the total cross 
sections from [<a href="Bibliography.php#refDon92" target="page">Don92</a>]. A few parameters allow some possibility 
to vary the basic setup. The allowed combinations of incoming particles 
are <i>p + p</i>, <i>pbar + p</i>, <i>pi+ + p</i>, <i>pi- + p</i>, 
<i>pi0/rho0 + p</i>, <i>phi + p</i>, <i>J/psi + p</i>, 
<i>rho + rho</i>, <i>rho + phi</i>, <i>rho + J/psi</i>, 
<i>phi + phi</i>, <i>phi + J/psi</i>, <i>J/psi + J/psi</i>, 
<i>Pomeron + p</i>, <i>gamma + gamma</i> and <i>gamma + p</i>. 
The strong emphasis on vector mesons is related to the description 
of <i>gamma + p</i> and <i>gamma + gamma</i> interactions in a 
Vector Dominance Model framework (which is not explicitly used in the 
current implementation of photoproduction, but is retained for potential 
future applications). 
 
<p/> 
The other options available for total, elastic and diffractive cross 
sections are: 
<ul> 
<li>A do-it-yourself selection of the main parameters.</li> 
<li>The MBR (Minimum Bias Rockefeller) model [<a href="Bibliography.php#refCie12" target="page">Cie12</a>], which 
is mainly intended for diffractive physics, but also parametrizes the 
total and elastic cross sections.</li> 
<li>The ABMST model [<a href="Bibliography.php#refApp16" target="page">App16</a>], which is based on a quite 
sophisticated Pomeron-inspired framework, and addresses total, elastic and 
single diffractive cross sections. The tuning to single diffractive 
data has mainly been performed at lower energies, so we also include 
variants that (hopefully) improves agreement with LHC data, and 
also introduce simple extensions to double and central diffraction.</li> 
<li>The RPP 2016 parametrization [<a href="Bibliography.php#refPat16" target="page">Pat16</a>], which is also 
Pomeron-inspired. It does not address diffractive cross sections.</li> 
</ul> 
 
<p/> 
The elastic cross section is differential in the squared momentum 
transfer <i>t</i>. The single diffractive additionally is differential 
in the mass of the diffractive system, or in <i>xi = x_Pom</i>, where 
<i>M^2_diff = xi * s</i>. For double diffraction the two masses can 
accordingly be related to <i>xi_1</i> and <i>xi_2</i> values. 
For central diffraction <i>M^2_diff = xi_1 * xi_2 * s</i>, and 
additionally the cross section is differential in <i>t_1</i> and 
<i>t_2</i>. 
 
<a name="section0"></a> 
<h3>Master switches</h3> 
 
The total and elastic cross sections are intimately connected via the 
optical theorem. Therefore the two should be calculated within a 
common setup. The diffractive cross sections are not as easily 
related, and can therefore be chosen separately, hence the two switches 
below. This allows different combinations to be tried out. 
 
<br/><br/><table><tr><td><strong>SigmaTotal:mode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 4</code>)</td></tr></table>
Choice of parametrization of the total and elastic cross sections. 
<br/>
<input type="radio" name="1" value="0"><strong>0 </strong>: Make your own choices (the "own model"), set as  fixed values.  <br/>
<input type="radio" name="1" value="1" checked="checked"><strong>1 </strong>: The DL model for total cross sections, extended to  more processes and to elastic cross sections according to SaS  ("SaS/DL").  <br/>
<input type="radio" name="1" value="2"><strong>2 </strong>: The MBR model for <ei>p p</ei> and <ei>p pbar</ei>,  else as option 1.  <br/>
<input type="radio" name="1" value="3"><strong>3 </strong>: The ABMST parametrizations for <ei>p p</ei> and  <ei>p pbar</ei>, else as option 1.  <br/>
<input type="radio" name="1" value="4"><strong>4 </strong>: The RPP2016 parametrizations for <ei>p p</ei> and  <ei>p pbar</ei>, else as option 1.  <br/>
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:mode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
Choice of parametrization of diffractive cross sections: single, 
double and central ditto. Note that there is no option 4. 
<br/>
<input type="radio" name="2" value="0"><strong>0 </strong>: Make your own choices, set as fixed values.  <br/>
<input type="radio" name="2" value="1" checked="checked"><strong>1 </strong>: The SaS parametrizations, available for a larger  set of incoming hadron combinations.  <br/>
<input type="radio" name="2" value="2"><strong>2 </strong>: The MBR model for <ei>p p</ei> and <ei>p pbar</ei>,  else as option 1.  <br/>
<input type="radio" name="2" value="3"><strong>3 </strong>: The ABMST parametrizations for <ei>p p</ei> and  <ei>p pbar</ei>, else as option 1.  <br/>
 
<p/> 
Note that the total cross section subtracted by the elastic and various 
diffractive ones gives the inelastic nondiffractive cross section, 
which therefore is not set separately. However, since the nondiffractive 
inelastic cross section is what makes up the minimum-bias event class, 
and plays a major role in the description of multiparton interactions, 
it is important that a consistent set is used. 
 
<p/> 
In the following subsections all the parameters available for the 
various values of the master switches are described. A final subsection 
covers the possibility to include Coulomb corrections in elastic scattering, 
and is relevant for all scenarios. 
 
<a name="section1"></a> 
<h3>Set your own cross sections</h3> 
 
The following four parameters can be set for the 
<code>SigmaTotal:mode = 0</code> option. The default values 
are in the right ballpark for LHC physics, but precise numbers 
depend on the energy used. 
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaTot </td><td></td><td> <input type="text" name="3" value="100." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The assumed total cross section in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaEl </td><td></td><td> <input type="text" name="4" value="25." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>25.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The assumed elastic cross section in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaElastic:bSlope </td><td></td><td> <input type="text" name="5" value="18." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>18.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The assumed slope <i>b</i> of the strong-interaction term 
<i>exp(bt)</i>, in units of GeV^-2. 
   
 
<br/><br/><table><tr><td><strong>SigmaElastic:rho </td><td></td><td> <input type="text" name="6" value="0.13" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.13</strong></code>; <code>minimum = -1.</code>; <code>maximum = 1.</code>)</td></tr></table>
The assumed ratio of the real to the imaginary parts of the nuclear 
scattering amplitude. This value is also used in the SaS/DL option. 
   
 
<p/> 
The following four parameters can be set for the 
<code>SigmaDiffractive:mode = 0</code> option. Again the default 
values are in the right ballpark for LHC physics, but with a 
considerable measure of uncertainty. 
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaXB </td><td></td><td> <input type="text" name="7" value="8." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>8.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Single Diffractive cross section <i>A + B &rarr; X + B</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaAX </td><td></td><td> <input type="text" name="8" value="8." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>8.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Single Diffractive cross section <i>A + B &rarr; A + X</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaXX </td><td></td><td> <input type="text" name="9" value="4." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Double Diffractive cross section <i>A + B &rarr; X_1 + X_2</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaAXB </td><td></td><td> <input type="text" name="10" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Central Diffractive cross section <i>A + B &rarr; A + X + B</i> in mb. 
   
 
<p/> 
The key parameter to set the differential shape of single diffraction 
is the <code>SigmaDiffractive:PomFlux</code> switch below. Seven different 
options are included, that provide the differential shape in diffractive 
mass and <i>t</i> of the scattered proton, based on the assumed Pomeron 
flux parametrizations. Only the SaS option contains a (published) 
extension to double diffraction, but the other alternatives have been 
extended in a minimal manner consistent with Pomeron phenomenology. 
These basic shapes can be further modified by the other settings below. 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:PomFlux  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 7</code>)</td></tr></table>
Parametrization of the Pomeron flux <ei>f_Pom/p( x_Pom, t)</ei>. 
<br/>
<input type="radio" name="11" value="1" checked="checked"><strong>1 </strong>:  Schuler and Sj&ouml;strand <ref>Sch94</ref>: based on a  critical Pomeron, giving a mass spectrum roughly like <ei>dm^2/m^2</ei>;  a mass-dependent exponential <ei>t</ei> slope that reduces the rate  of low-mass states.  <br/>
<input type="radio" name="11" value="2"><strong>2 </strong>:  Bruni and Ingelman <ref>Bru93</ref>: also a critical  Pomeron giving close to <ei>dm^2/m^2</ei>,  with a <ei>t</ei> distribution  the sum of two exponentials.  <br/>
<input type="radio" name="11" value="3"><strong>3 </strong>:  a conventional Pomeron description, in the RapGap  manual <ref>Jun95</ref> attributed to Berger et al. and Streng  <ref>Ber87a</ref>, but there (and here) with values updated to a  supercritical Pomeron with <ei>epsilon &gt; 0</ei> (see below),  which gives a stronger peaking towards low-mass diffractive states,  and with a mass-dependent (the <ei>alpha'</ei> below) exponential  <ei>t</ei> slope.  <br/>
<input type="radio" name="11" value="4"><strong>4 </strong>:  a conventional Pomeron description, attributed to  Donnachie and Landshoff <ref>Don84</ref>, again with supercritical Pomeron,  with the same two parameters as option 3 above, but this time with a  power-law <ei>t</ei> distribution.  <br/>
<input type="radio" name="11" value="5"><strong>5 </strong>:  the MBR simulation of (anti)proton-proton interactions  <ref>Cie12</ref>. The mass distribution follows a renormalized-Regge-theory  model, successfully tested using CDF data.  <br/>
<input type="radio" name="11" value="6"><strong>6 </strong>:  The H1 Fit A parametrisation of the Pomeron flux  <ref>H1P06,H1P06a</ref>. The flux factors are motivated by Regge theory,  assuming a Regge trajectory as in options 3 and 4. The flux has been  normalised to 1 at <ei>x_Pomeron = 0.003</ei> and slope parameter and  Pomeron intercept has been fitted to H1 data.  <br/>
<input type="radio" name="11" value="7"><strong>7 </strong>:  The H1 Fit B parametrisation of the Pomeron flux  <ref>H1P06,H1P06a</ref>.  <br/>
 
<p/> 
In options 3, 4, 6, and 7 above, the Pomeron Regge trajectory is 
parametrized as 
<br/><i> 
alpha(t) = 1 + epsilon + alpha' t 
</i><br/> 
The <i>epsilon</i> and <i>alpha'</i> parameters can be set 
separately in options 3 and 4, and additionally <i>alpha'</i> 
is set in option 1, while values are fixed in options 6 and 7: 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:PomFluxEpsilon </td><td></td><td> <input type="text" name="12" value="0.085" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.085</strong></code>; <code>minimum = 0.02</code>; <code>maximum = 0.15</code>)</td></tr></table>
The Pomeron trajectory intercept <i>epsilon</i> above for the 3 and 4 
flux options. For technical reasons <i>epsilon &gt; 0</i> is necessary 
in the current implementation. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:PomFluxAlphaPrime </td><td></td><td> <input type="text" name="13" value="0.25" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.25</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 0.4</code>)</td></tr></table>
The Pomeron trajectory slope <i>alpha'</i> above for the 1, 3 and 4 
flux options. 
   
 
<p/> 
The options above might give vanishing (or even negative) <i>b</i> 
slope values, and also do not enforce the presence of a rapidity gap. 
Furthermore the lowest allowed central diffractive mass is not well-defined; 
it would not be meaningful to go all the way down to the <i>pi pi</i> 
kinematical limit, since exclusive states are not modelled. Therefore 
the following parameters have been introduced to address such issues. 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:OwnbMinDD </td><td></td><td> <input type="text" name="14" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.5</code>; <code>maximum = 5.</code>)</td></tr></table>
In the options with a simple <i>exp(b * t)</i> falloff for the <i>t</i> 
spectrum, ensure that <i>b</i> is at least this large. (Recall that the 
<i>b</i> formula typically contains one term for each incoming hadron 
that does not break up, and for double diffraction such terms are absent. 
This leaves only the pomeron propagator part, which often vanishes in 
the limit of vanishing rapidity gap.) 
   
 
<br/><br/><strong>SigmaDiffractive:OwndampenGap</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Switch on damping of small rapidity gaps in single, double and central 
diffraction. The reason for this option is that the separation between 
diffraction and nondiffraction is blurred for events with small gaps. 
Therefore a damping factor for small gaps is imposed with this option, 
of the form 
<br/><i> 
  1 / (1 + exp( -p * (y - y_gap))) = 1 / (1 + exp(p * y_gap) * (exp(-y))^p), 
</i><br/> 
where <i>y</i> is the rapidity gap(s) in the current event, and 
<i>p</i> and <i>y_gap</i> are two parameters. Thus the damping 
kicks in for <i>y &lt; y_gap</i>, and the transition region from small to 
large damping is of order <i>1/p</i> in <i>y</i>. The <i>exp(-y)</i> 
values are <i>xi</i> for SD, <i>xi_1 * xi_2 * s / m_p^2</i> for DD, 
and <i>xi_1</i> and <i>xi_2</i> for CD. The two parameters of the 
damping are described below. 
<br/><b>Note:</b> if the integrated diffractive cross sections are kept 
fixed, switching on this option will increase the rate of diffractive 
events with large rapidity gaps, so do consistent changes. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:Ownygap </td><td></td><td> <input type="text" name="16" value="2." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.</strong></code>; <code>minimum = 0.1</code>)</td></tr></table>
Assume a damping of small rapidity gaps, as described above, to set in 
around the value <i>y_gap</i> given by this parameter. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:Ownypow </td><td></td><td> <input type="text" name="17" value="5." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5.</strong></code>; <code>minimum = 0.5</code>)</td></tr></table>
Assume a damping of small rapidity gaps, as described above, to set in 
over a rapidity region of width <i>1/p</i>, with <i>p</i> given by 
this parameter. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:OwnmMinCD </td><td></td><td> <input type="text" name="18" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.5</code>)</td></tr></table>
The smallest allowed central diffractive mass, with a sharp cut at 
this value. 
   
 
<a name="section2"></a> 
<h3>Modify the SaS/DL cross sections</h3> 
 
The default description of total, elastic and diffractive interactions was 
parameterized and fit in [<a href="Bibliography.php#refSch94" target="page">Sch94</a>, <a href="Bibliography.php#refSch97" target="page">Sch97</a>]. There is no freedom for 
total and elastic cross sections, except that the <i>rho</i> parameter 
is not modelled but taken from the <code>SigmaElastic:rho</code> 
parameter above. 
 
<p/> 
The following three parameters allow for some modification of the mass 
distribution of the diffractive system, relative to the default setup. 
The parametrized cross sections explicitly depend on them, so that 
integrated diffractive cross section are changed acordingly. 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:mMin </td><td></td><td> <input type="text" name="19" value="0.28" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.28</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Lowest mass of a single or double diffractive system is set to be 
<i>mHadron + mMin</i>. 
   
 
<br/><br/><table><tr><td><strong> </td><td></td><td> <input type="text" name="20" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Normalization factor for the contribution of low-mass resonances 
to the diffractive cross section (<i>cRes</i> in eq. (22) of 
[<a href="Bibliography.php#refSch94" target="page">Sch94</a>]). 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:mResMax </td><td></td><td> <input type="text" name="21" value="1.062" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.062</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
The contribution of low-mass resonances is dampened at around the 
scale <i>mHadron + mResMax</i> (the sum is <i>Mres</i> in eq. (22) 
of [<a href="Bibliography.php#refSch94" target="page">Sch94</a>]). To make sense, we should have 
<i>mResMax > mMin</i>. 
   
 
<p/> 
Central diffraction (CD) was not part of the framework in [<a href="Bibliography.php#refSch94" target="page">Sch94</a>]. 
It has now been added for <i>p p</i> or <i>pbar p</i>, but only for 
multiparticle states, i.e. excluding the low-mass resonance region below 
roughly 1 GeV, as well as other exclusive states. It uses the same 
proton-Pomeron vertex as in single diffraction, twice, to describe 
<i>x_Pomeron</i> and <i>t</i> spectra. This fixes the energy 
dependence, which has been integrated and parametrized. The absolute 
normalization has been left open, however. Furthermore, since CD has not 
been included in previous tunes to data, a special flag is available to 
reproduce the old behaviour (with due complications when one does not want 
to do this). 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:mMinCD </td><td></td><td> <input type="text" name="22" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.5</code>)</td></tr></table>
The smallest allowed central diffractive mass, with a sharp cut at 
this value. 
   
 
<br/><br/><table><tr><td><strong>SigmaTotal:sigmaAXB2TeV </td><td></td><td> <input type="text" name="23" value="1.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.5</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The CD cross section for <i>p p</i> and <i>pbar p</i> collisions, 
normalized to its value at 2 TeV CM energy, expressed in mb. The energy 
dependence is then parametrized, and behaves roughly like 
<i>ln^1.5(s)</i>. 
   
 
<br/><br/><strong>SigmaTotal:zeroAXB</strong>  <input type="radio" name="24" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="24" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
several existing <?php $filepath = $_GET["filepath"];
echo "<a href='Tunes.php?filepath=".$filepath."' target='page'>";?>tunes</a> do not include CD. 
An inclusion of a nonvanishing CD cross section directly affects 
the nondiffractive phenomenology, even if not dramatically, and so 
this flag is used to forcibly set the CD cross section to vanish 
in such tunes. You can switch CD back on <i>after</i> the selection of 
a tune, if you so wish, by resetting <code>SigmaTotal:zeroAXB = off</code>. 
   
 
<p/> 
LHC data have suggested that diffractive cross sections rise slower than 
predicted in the original studies. A likely reason is that unitarization 
effects may dampen the rise of diffractive cross sections relative to 
the default parametrizations. The settings here allows one way to 
introduce a dampening, which is used in some of the existing 
<?php $filepath = $_GET["filepath"];
echo "<a href='Tunes.php?filepath=".$filepath."' target='page'>";?>tunes</a>. 
 
<br/><br/><strong>SigmaDiffractive:dampen</strong>  <input type="radio" name="25" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="25" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow a user to dampen diffractive cross sections; on/off = true/false. 
   
 
<p/> 
When <code>SigmaDiffractive:dampen = on</code>, the three diffractive 
cross sections are damped so that they never can exceed the respective 
values below. Specifically, if the standard parametrization gives 
the cross section <i>sigma_old(s)</i> and a fixed <i>sigma_max</i> 
is set, the actual cross section becomes 
<br/><i> 
sigma_new(s) = sigma_old(s) * sigma_max / (sigma_old(s) + sigma_max). 
</i><br/> 
This reduces to <i>sigma_old(s)</i> at low energies and to 
<i>sigma_max</i> at high ones. Note that the asymptotic value 
is approached quite slowly, however. 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:maxXB </td><td></td><td> <input type="text" name="26" value="65." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>65.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The above <i>sigma_max</i> for <i>A + B &rarr; X + B</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:maxAX </td><td></td><td> <input type="text" name="27" value="65." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>65.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The above <i>sigma_max</i> for <i>A + B &rarr; A + X</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:maxXX </td><td></td><td> <input type="text" name="28" value="65." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>65.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The above <i>sigma_max</i> for <i>A + B &rarr; X_1 + X_2</i> in mb. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:maxAXB </td><td></td><td> <input type="text" name="29" value="3." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The above <i>sigma_max</i> for <i>A + B &rarr; A + X + B</i> in mb. 
   
 
<p/> 
As above, a reduced diffractive cross section automatically translates 
into an increased nondiffractive one, such that the total (and elastic) 
cross section remains fixed. 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:SaSepsilon </td><td></td><td> <input type="text" name="30" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = -0.2</code>; <code>maximum = 0.2</code>)</td></tr></table>
The SaS ansatz starts out from a <i>dM^2/M^2</i> shape of diffractive 
spectra, a shape that then is modified by <i>t</i>-spectra integration 
and small-mass enhancement. For exploratory purposes it is possible to 
modify the base ansatz to be <i>dM^2/M^(2 * (1 + epsilon))</i>. In 
principle the integrated diffractive cross sections ought to be 
recalculated accordingly, but for simplicity they are not modified. 
   
 
<a name="section3"></a> 
<h3>Modify the MBR cross sections</h3> 
 
The MBR differential cross section also comes with a selection of 
parameters that can be changed from their default values, to modify 
diffractive event rates and shapes, while the total and elastic cross 
sections remain unaffected. These parameters are described in the 
following. 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:MBRepsilon </td><td></td><td> <input type="text" name="31" value="0.104" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.104</strong></code>; <code>minimum = 0.02</code>; <code>maximum = 0.15</code>)</td></tr></table>
<a name="anchor1"></a>
<code>parm&nbsp; </code><strong> SigmaDiffractive:MBRalpha &nbsp;</strong> 
 (<code>default = <strong>0.25</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 0.4</code>)<br/>
the parameters of the Pomeron trajectory. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:MBRbeta0 </td><td></td><td> <input type="text" name="32" value="6.566" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>6.566</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 10.0</code>)</td></tr></table>
<a name="anchor2"></a>
<code>parm&nbsp; </code><strong> SigmaDiffractive:MBRsigma0 &nbsp;</strong> 
 (<code>default = <strong>2.82</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)<br/>
the Pomeron-proton coupling, and the total Pomeron-proton cross section. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:MBRm2Min </td><td></td><td> <input type="text" name="33" value="1.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.5</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 3.0</code>)</td></tr></table>
the lowest value of the mass squared of the dissociated system, including 
central diffraction. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:MBRdyminSDflux </td><td></td><td> <input type="text" name="34" value="2.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.3</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)</td></tr></table>
<a name="anchor3"></a>
<code>parm&nbsp; </code><strong> SigmaDiffractive:MBRdyminDDflux &nbsp;</strong> 
 (<code>default = <strong>2.3</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)<br/>
<a name="anchor4"></a>
<code>parm&nbsp; </code><strong> SigmaDiffractive:MBRdyminCDflux &nbsp;</strong> 
 (<code>default = <strong>2.3</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)<br/>
the minimum width of the rapidity gap used in the calculation of 
<i>Ngap(s)</i> (flux renormalization). 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:MBRdyminSD </td><td></td><td> <input type="text" name="35" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)</td></tr></table>
<a name="anchor5"></a>
<code>parm&nbsp; </code><strong> SigmaDiffractive:MBRdyminDD &nbsp;</strong> 
 (<code>default = <strong>2.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)<br/>
<a name="anchor6"></a>
<code>parm&nbsp; </code><strong> SigmaDiffractive:MBRdyminCD &nbsp;</strong> 
 (<code>default = <strong>2.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 5.0</code>)<br/>
the minimum width of the rapidity gap used in the calculation of cross 
sections, i.e. the parameter <i>dy_S</i>, which suppresses the cross 
section at low <i>dy</i> (non-diffractive region). The cross section 
is damped smoothly, such that it is suppressed by a factor of a half 
at around this scale. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:MBRdyminSigSD </td><td></td><td> <input type="text" name="36" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.001</code>; <code>maximum = 5.0</code>)</td></tr></table>
<a name="anchor7"></a>
<code>parm&nbsp; </code><strong> SigmaDiffractive:MBRdyminSigDD &nbsp;</strong> 
 (<code>default = <strong>0.5</strong></code>; <code>minimum = 0.001</code>; <code>maximum = 5.0</code>)<br/>
<a name="anchor8"></a>
<code>parm&nbsp; </code><strong> SigmaDiffractive:MBRdyminSigCD &nbsp;</strong> 
 (<code>default = <strong>0.5</strong></code>; <code>minimum = 0.001</code>; <code>maximum = 5.0</code>)<br/>
the parameter <i>sigma_S</i>, used for the cross section suppression at 
low <i>dy</i> (non-diffractive region). The smaller this value, the  
narrow the rapidity region over which the suppression sets in. 
   
 
<a name="section4"></a> 
<h3>Modify the ABMST cross sections</h3> 
 
The ABMST model provides a detailed description of the total, elastic 
and single diffractive cross sections. The former two components are 
accepted as is, while we have allowed alternative shapes for single 
diffraction, notably to enforce a rapidity gap. The ABMST model does 
not address double and central diffraction, so we have extended it on 
our own, as described below. 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTmodeSD  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
Setup of single diffraction in the ABMST scenario. 
<br/>
<input type="radio" name="37" value="0"><strong>0 </strong>: Keep the pure ABMST ansatz, which notably vanishes above  <ei>|t| = 4 GeV^2</ei>, and has a constant term up to that scale.  <br/>
<input type="radio" name="37" value="1" checked="checked"><strong>1 </strong>: Use a slightly modified ansatz without an upper  <ei>|t|</ei> cut, but instead an exponential fall-off that gives the same  integrated diffractive rate and average <ei>|t|</ei> value. In addition  the low-mass background term is modified as a combination of a linear and  a quadratic term, instead of a qudratic only.  <br/>
<input type="radio" name="37" value="2"><strong>2 </strong>:  Option 0, with a scaling factor of  <ei>k * (s / m_p^2)^q</ei>,  where <ei>k</ei> is <code>SigmaDiffractive:multSD</code> and  <ei>q</ei> is <code>SigmaDiffractive:powSD</code>  <br/>
<input type="radio" name="37" value="3"><strong>3 </strong>:  Option 1, with a scaling factor of  <ei>k * (s / m_p^2)^q</ei>,  where <ei>k</ei> is <code>SigmaDiffractive:multSD</code> and  <ei>q</ei> is <code>SigmaDiffractive:powSD</code>  <br/>
<br/><b>Note:</b> also the <code>SigmaDiffractive:ABMSTdampenGap</code> 
and <code>SigmaDiffractive:ABMSTuseBMin</code> flags below very much affect 
the behaviour; you have to switch them off and use option 0 above to recover 
the pure ABMST model. 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTmultSD </td><td></td><td> <input type="text" name="38" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.01</code>)</td></tr></table>
possibility to rescale the double diffractive cross section by a factor 
<i>k</i> as described above. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTpowSD </td><td></td><td> <input type="text" name="39" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 0.25</code>)</td></tr></table>
possibility to rescale the double diffractive cross section by a factor 
<i>(s / m_p^2)^q</i>, as described above, with <i>q</i> set here. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTmodeDD  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
Setup of double diffraction in the ABMST scenario. Note that ABMST does 
not provide any answer here, so the single-diffractive framework is 
extended by a simple factorized ansatz 
<eq> 
  dsigma_DD( xi_1, xi_2, t) / (dxi_1 dxi_2 dt) 
    = dsigma_SD (xi_1, t) / (dxi_1 dt) * dsigma_SD (xi_2, t) / (dxi_2 dt) 
    / (dsigma_El( t) / dt) . 
</eq> 
The above ansatz is marred by the dip in <ei>dsigma_El / dt</ei> 
by destructive interference, however, so in this extension we only allow 
for Pomerons in the elastic cross section, which is intended to represent 
the bulk of the cross section. As such, the equation gives a parameter-free 
prediction for the double diffractive cross section. For flexibility we 
introduce a (default) option where the absolute normalization can be 
modified, while retaining the shape of the ansatz. 
<br/>
<input type="radio" name="40" value="0"><strong>0 </strong>: Describe the double diffractive cross section by the  simple factorized ansatz introduced above, within the allowed phase-space  limits. Note that the single diffractive cross section is affected by the  choice made for <code>SigmaDiffractive:ABMSTmodeSD</code>.  <br/>
<input type="radio" name="40" value="1" checked="checked"><strong>1 </strong>: The double diffractive  cross section can be rescaled by a factor <ei>k * (s / m_p^2)^q</ei>,  where <ei>k</ei> is <code>SigmaDiffractive:multDD</code> and  <ei>q</ei> is <code>SigmaDiffractive:powDD</code>.  <br/>
<br/><b>Note:</b> also the <code>SigmaDiffractive:ABMSTdampenGap</code> 
and <code>SigmaDiffractive:ABMSTuseBMin</code> flags below very much affect 
the behaviour. 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTmultDD </td><td></td><td> <input type="text" name="41" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.01</code>)</td></tr></table>
possibility to rescale the double diffractive cross section by a factor 
<i>k</i> as described above. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTpowDD </td><td></td><td> <input type="text" name="42" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 0.25</code>)</td></tr></table>
possibility to rescale the double diffractive cross section by a factor 
<i>(s / m_p^2)^q</i>, as described above, with <i>q</i> set here. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTmodeCD  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
Setup of central diffraction in the ABMST scenario. Note that ABMST does 
not provide any answer here, so the single-diffractive framework is 
extended by a simple factorized ansatz 
<eq> 
  dsigma_CD( xi_1, xi_2, t_1, t_2) / (dxi_1 dxi_2 dt_1 dt2_) 
    = dsigma_SD (xi_1, t_1) / (dxi_1 dt_1) 
    * dsigma_SD (xi_2, t_2) / (dxi_2 dt_2) / sigma_total(s) , 
</eq> 
and again a variant is introduced below. 
<br/>
<input type="radio" name="43" value="0" checked="checked"><strong>0 </strong>: Describe the central diffractive cross section by the  simple factorized ansatz introduced above, within the allowed phase-space  limits. Also here, we only allow for Pomerons in the total cross section.  Note that the single diffractive cross section is affected by the  choice made for <code>SigmaDiffractive:ABMSTmodeSD</code>.  <br/>
<input type="radio" name="43" value="1"><strong>1 </strong>: In addition to option 0, the central diffractive  cross section can be rescaled by a factor <ei>k * (s / m_p^2)^q</ei>,  where <ei>k</ei> is <code>SigmaDiffractive:multCD</code> and  <ei>q</ei> is <code>SigmaDiffractive:powCD</code>.  <br/>
<br/><b>Note:</b> also the <code>SigmaDiffractive:ABMSTdampenGap</code> 
and <code>SigmaDiffractive:ABMSTuseBMin</code> flags below very much affect 
the behaviour. 
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTmultCD </td><td></td><td> <input type="text" name="44" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.01</code>)</td></tr></table>
possibility to rescale the central diffractive cross section by a factor 
<i>k</i> as described above. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTpowCD </td><td></td><td> <input type="text" name="45" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 0.25</code>)</td></tr></table>
possibility to rescale the central diffractive cross section by a factor 
<i>(s / m_p^2)^q</i>, as described above, with <i>q</i> set here. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTmMinCD </td><td></td><td> <input type="text" name="46" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.5</code>)</td></tr></table>
The smallest allowed central diffractive mass, with a sharp cut at 
this value. 
   
 
<br/><br/><strong>SigmaDiffractive:ABMSTdampenGap</strong>  <input type="radio" name="47" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="47" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Switch on damping of small rapidity gaps in single, double and central 
diffraction. The reason for this option, on by default, is that the 
the ABMST SD ansats contains terms that peak near <i>xi = 1</i>. 
This leads to very large integrated SD cross sections at higher energies, 
such that the diffractive cross section is larger than the nondiffractive 
one. It then becomes a challenge e.g. how to implement and interpret PDFs, 
which by definition are inclusive, but would have to be split consistently 
between the different contributions. (For the hard-jet subsample it can be 
done e.g. as in [<a href="Bibliography.php#refRas16" target="page">Ras16</a>], but it would be  complicated for 
softer jets in the MPI context.) Furthermore the separation between 
diffraction and nondiffraction is blurred for events with small gaps. 
Therefore a damping factor for small gaps is imposed with this option, 
of the form 
<br/><i> 
  1 / (1 + exp( -p * (y - y_gap))) = 1 / (1 + exp(p * y_gap) * (exp(-y))^p), 
</i><br/> 
where <i>y</i> is the rapidity gap(s) in the current event, and 
<i>p</i> and <i>y_gap</i> are two parameters. Thus the damping 
kicks in for <i>y &lt; y_gap</i>, and the transition region from small to 
large damping is of order <i>1/p</i> in <i>y</i>. The <i>exp(-y)</i> 
values are <i>xi</i> for SD, <i>xi_1 * xi_2 * s / m_p^2</i> for DD, 
and <i>xi_1</i> and <i>xi_2</i> for CD. The two parameters of the 
damping are described below. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTygap </td><td></td><td> <input type="text" name="48" value="2." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.</strong></code>; <code>minimum = 0.1</code>)</td></tr></table>
Assume a damping of small rapidity gaps in the ABMST model, as described 
above, to set in around the value <i>y_gap</i> given by this parameter. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTypow </td><td></td><td> <input type="text" name="49" value="5." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5.</strong></code>; <code>minimum = 0.5</code>)</td></tr></table>
Assume a damping of small rapidity gaps in the ABMST model, as described 
above, to set in over a rapidity region of width <i>1/p</i>, with 
<i>p</i> given by this parameter. 
   
 
<br/><br/><strong>SigmaDiffractive:ABMSTuseBMin</strong>  <input type="radio" name="50" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="50" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
The slope <i>b</i> of an approximate <i>exp(b * t)</i> fall-off 
is <i>xi</i>-dependent in the ABMST model for single diffraction. 
In particular it can become close to zero for large <i>xi</i>, which 
means that the <i>t</i>-integrated cross section becomes very large. 
While the general trend is reasonable, the behaviour in the 
<i>xi &rarr; 1</i> limit is questionable. Therefore it makes sense to 
impose some minimal <i>b</i> slope. For double diffraction such 
issues become even more pressing, since the division by the elastic 
cross section could even lead to a negative <i>b</i> slope, which 
would not be physical. The central diffractive cross section is more 
well-behaved, but for consistency it is meaningful to ensure a minimal 
fall-off also here. Therefore, when this flag is on, a minimal fall-off 
<i>exp(b_min * t)</i> is assumed for each of the three components, 
with the respective <i>b_min</i> value stored in the three parameters 
below. The fall-off is defined relative to the value at <i>t = 0</i>, 
a point that is outside the physical region, but the parametrization 
of the diffractive cross sections can still be used there meaningfully. 
Only positive <i>b_min</i> values are acted on, so the SD/DD/CD 
components can be switched off individually even when this flag is on. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTbMinSD </td><td></td><td> <input type="text" name="51" value="2." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.</strong></code>)</td></tr></table>
Assume a minimal fall-off <i>exp(b_min * t)</i> in the ABMST model 
for single diffraction, as described above. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTbMinDD </td><td></td><td> <input type="text" name="52" value="2." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.</strong></code>)</td></tr></table>
Assume a minimal fall-off <i>exp(b_min * t)</i> in the extension of 
the ABMST model to double diffraction, as described above. 
   
 
<br/><br/><table><tr><td><strong>SigmaDiffractive:ABMSTbMinCD </td><td></td><td> <input type="text" name="53" value="2." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.</strong></code>)</td></tr></table>
Assume a minimal fall-off <i>exp(b_min * (t_1 + t_2))</i> in the extension 
of the ABMST model to central diffraction, as described above. 
   
 
<a name="section5"></a> 
<h3>Modify the RPP cross sections</h3> 
 
The RPP approach only addresses total and (differential) elastic 
cross sections, and there are no free parameters that can be changed. 
 
<a name="section6"></a> 
<h3>Coulomb corrections to elastic scattering</h3> 
 
By default there is no Coulomb-term contribution to the elastic 
(or total) cross section, which of course becomes infinite if this 
contribution is included in the collision between charged particles, 
owing to the <i>1/t^2</i> singularity of <i>t</i>-channel photon 
exchange. You can switch on Coulomb corrections below, however, including 
interference with the conventional strong-interaction term. 
The own, SaS/DL and MBR models share a common machinery to evaluate the 
interference [<a href="Bibliography.php#refBer87" target="page">Ber87</a>], while ABMST and RPP use a slighly 
different expression for the (poorly known) difference in phases 
between the hadronic and the electromagnetic amplitudes. 
 
<br/><br/><strong>SigmaElastic:Coulomb</strong>  <input type="radio" name="54" value="on"><strong>On</strong>
<input type="radio" name="54" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Include Coulomb corrections to the elastic and total cross sections. 
   
 
<br/><br/><table><tr><td><strong>SigmaElastic:tAbsMin </td><td></td><td> <input type="text" name="55" value="5e-5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5e-5</strong></code>; <code>minimum = 1e-10</code>; <code>maximum = 1e-3</code>)</td></tr></table>
since the Coulomb contribution is infinite a lower limit on <i>|t|</i> 
must be set to regularize the divergence, in units of GeV^2. 
This means that the elastic and total cross sections are reduced by 
the amount of the ordinary cross section in the cut-out region, 
but increased by the Coulomb contribution itself and the interference 
term (of either sign). This variable has no effect if Coulomb corrections 
are not switched on or not relevant (e.g. for neutral particles), i.e. 
then <i>t = 0</i> sets the limit. 
   
 
<br/><br/><table><tr><td><strong>SigmaElastic:lambda </td><td></td><td> <input type="text" name="56" value="0.71" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.71</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 2.</code>)</td></tr></table>
the main parameter of the electric form factor 
<i>G(t) = lambda^2 / (lambda + |t|)^2</i>, in units of GeV^2, 
as used in the own, SaS/DL and MBR models. 
   
 
<br/><br/><table><tr><td><strong>SigmaElastic:phaseConst </td><td></td><td> <input type="text" name="57" value="0.577" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.577</strong></code>)</td></tr></table>
The Coulomb term is taken to contain a phase factor 
<i>exp(+- i alpha phi(t))</i>, with + for <i>p p</i> and - for 
<i>pbar p</i>, where <i>phi(t) = - phaseConst - ln(-B t/2)</i>. 
This constant is model dependent [<a href="Bibliography.php#refCah82" target="page">Cah82</a>]. This expression 
is used in the own, SaS/DL and MBR models, where the hadronic cross 
section is modelled as a simple <i>exp(B t)</i>. 
   
 
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

if($_POST["1"] != "1")
{
$data = "SigmaTotal:mode = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "1")
{
$data = "SigmaDiffractive:mode = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "100.")
{
$data = "SigmaTotal:sigmaTot = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "25.")
{
$data = "SigmaTotal:sigmaEl = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "18.")
{
$data = "SigmaElastic:bSlope = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.13")
{
$data = "SigmaElastic:rho = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "8.")
{
$data = "SigmaTotal:sigmaXB = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "8.")
{
$data = "SigmaTotal:sigmaAX = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "4.")
{
$data = "SigmaTotal:sigmaXX = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "1.")
{
$data = "SigmaTotal:sigmaAXB = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "1")
{
$data = "SigmaDiffractive:PomFlux = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "0.085")
{
$data = "SigmaDiffractive:PomFluxEpsilon = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0.25")
{
$data = "SigmaDiffractive:PomFluxAlphaPrime = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "1.")
{
$data = "SigmaDiffractive:OwnbMinDD = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "SigmaDiffractive:OwndampenGap = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "2.")
{
$data = "SigmaDiffractive:Ownygap = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "5.")
{
$data = "SigmaDiffractive:Ownypow = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "1.")
{
$data = "SigmaDiffractive:OwnmMinCD = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "0.28")
{
$data = "SigmaDiffractive:mMin = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "2.0")
{
$data = " = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "1.062")
{
$data = "SigmaDiffractive:mResMax = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "1.")
{
$data = "SigmaDiffractive:mMinCD = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "1.5")
{
$data = "SigmaTotal:sigmaAXB2TeV = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "on")
{
$data = "SigmaTotal:zeroAXB = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "on")
{
$data = "SigmaDiffractive:dampen = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "65.")
{
$data = "SigmaDiffractive:maxXB = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "65.")
{
$data = "SigmaDiffractive:maxAX = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "65.")
{
$data = "SigmaDiffractive:maxXX = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "3.")
{
$data = "SigmaDiffractive:maxAXB = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "0.0")
{
$data = "SigmaDiffractive:SaSepsilon = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "0.104")
{
$data = "SigmaDiffractive:MBRepsilon = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "6.566")
{
$data = "SigmaDiffractive:MBRbeta0 = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "1.5")
{
$data = "SigmaDiffractive:MBRm2Min = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "2.3")
{
$data = "SigmaDiffractive:MBRdyminSDflux = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "2.0")
{
$data = "SigmaDiffractive:MBRdyminSD = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "0.5")
{
$data = "SigmaDiffractive:MBRdyminSigSD = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "1")
{
$data = "SigmaDiffractive:ABMSTmodeSD = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "1.")
{
$data = "SigmaDiffractive:ABMSTmultSD = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "0.0")
{
$data = "SigmaDiffractive:ABMSTpowSD = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
if($_POST["40"] != "1")
{
$data = "SigmaDiffractive:ABMSTmodeDD = ".$_POST["40"]."\n";
fwrite($handle,$data);
}
if($_POST["41"] != "1.")
{
$data = "SigmaDiffractive:ABMSTmultDD = ".$_POST["41"]."\n";
fwrite($handle,$data);
}
if($_POST["42"] != "0.1")
{
$data = "SigmaDiffractive:ABMSTpowDD = ".$_POST["42"]."\n";
fwrite($handle,$data);
}
if($_POST["43"] != "0")
{
$data = "SigmaDiffractive:ABMSTmodeCD = ".$_POST["43"]."\n";
fwrite($handle,$data);
}
if($_POST["44"] != "1.")
{
$data = "SigmaDiffractive:ABMSTmultCD = ".$_POST["44"]."\n";
fwrite($handle,$data);
}
if($_POST["45"] != "0.1")
{
$data = "SigmaDiffractive:ABMSTpowCD = ".$_POST["45"]."\n";
fwrite($handle,$data);
}
if($_POST["46"] != "1.")
{
$data = "SigmaDiffractive:ABMSTmMinCD = ".$_POST["46"]."\n";
fwrite($handle,$data);
}
if($_POST["47"] != "on")
{
$data = "SigmaDiffractive:ABMSTdampenGap = ".$_POST["47"]."\n";
fwrite($handle,$data);
}
if($_POST["48"] != "2.")
{
$data = "SigmaDiffractive:ABMSTygap = ".$_POST["48"]."\n";
fwrite($handle,$data);
}
if($_POST["49"] != "5.")
{
$data = "SigmaDiffractive:ABMSTypow = ".$_POST["49"]."\n";
fwrite($handle,$data);
}
if($_POST["50"] != "on")
{
$data = "SigmaDiffractive:ABMSTuseBMin = ".$_POST["50"]."\n";
fwrite($handle,$data);
}
if($_POST["51"] != "2.")
{
$data = "SigmaDiffractive:ABMSTbMinSD = ".$_POST["51"]."\n";
fwrite($handle,$data);
}
if($_POST["52"] != "2.")
{
$data = "SigmaDiffractive:ABMSTbMinDD = ".$_POST["52"]."\n";
fwrite($handle,$data);
}
if($_POST["53"] != "2.")
{
$data = "SigmaDiffractive:ABMSTbMinCD = ".$_POST["53"]."\n";
fwrite($handle,$data);
}
if($_POST["54"] != "off")
{
$data = "SigmaElastic:Coulomb = ".$_POST["54"]."\n";
fwrite($handle,$data);
}
if($_POST["55"] != "5e-5")
{
$data = "SigmaElastic:tAbsMin = ".$_POST["55"]."\n";
fwrite($handle,$data);
}
if($_POST["56"] != "0.71")
{
$data = "SigmaElastic:lambda = ".$_POST["56"]."\n";
fwrite($handle,$data);
}
if($_POST["57"] != "0.577")
{
$data = "SigmaElastic:phaseConst = ".$_POST["57"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
