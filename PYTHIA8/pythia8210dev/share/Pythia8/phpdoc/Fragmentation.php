<html>
<head>
<title>Fragmentation</title>
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

<form method='post' action='Fragmentation.php'>
 
<h2>Fragmentation</h2> 
 
Fragmentation in PYTHIA is based on the Lund string model 
[<a href="Bibliography.php" target="page">And83, Sjo84</a>]. Several different aspects are involved in 
the physics description, which  here therefore is split accordingly. 
This also, at least partly, reflect the set of classes involved in 
the fragmentation machinery. 
 
<p/> 
The variables collected here have a very wide span of usefulness. 
Some would be central in any hadronization tuning exercise, others 
should not be touched except by experts. 
 
<p/> 
The fragmentation flavour-choice machinery is also used in a few 
other places of the program, notably particle decays, and is thus 
described on the separate <?php $filepath = $_GET["filepath"];
echo "<a href='FlavourSelection.php?filepath=".$filepath."' target='page'>";?>Flavour 
Selection</a> page. 
 
<h3>Fragmentation functions</h3> 
 
The <code>StringZ</code> class handles the choice of longitudinal 
lightcone fraction <i>z</i> according to one of two possible 
shape sets. 
 
<p/> 
The Lund symmetric fragmentation function [<a href="Bibliography.php" target="page">And83</a>] is the 
only alternative for light quarks. It is of the form 
<br/><i> 
    f(z) = (1/z) * (1-z)^a * exp(-b m_T^2 / z) 
</i><br/> 
with the two main free parameters <i>a</i> and <i>b</i> to be 
tuned to data. They are stored in 
 
<br/><br/><table><tr><td><strong>StringZ:aLund </td><td></td><td> <input type="text" name="1" value="0.68" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.68</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
The <i>a</i> parameter of the Lund symmetric fragmentation function. 
   
 
<br/><br/><table><tr><td><strong>StringZ:bLund </td><td></td><td> <input type="text" name="2" value="0.98" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.98</strong></code>; <code>minimum = 0.2</code>; <code>maximum = 2.0</code>)</td></tr></table>
The <i>b</i> parameter of the Lund symmetric fragmentation function. 
   
 
<p/> 
In principle, each flavour can have a different <i>a</i>. Then, 
for going from an old flavour <i>i</i> to a new <i>j</i> one 
the shape is 
<br/><i> 
    f(z) = (1/z) * z^{a_i} * ((1-z)/z)^{a_j} * exp(-b * m_T^2 / z) 
</i><br/> 
This is only implemented for s quarks and diquarks relative to normal quarks: 
 
<br/><br/><table><tr><td><strong>StringZ:aExtraSQuark </td><td></td><td> <input type="text" name="3" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
allows a larger <i>a</i> for s quarks, with total 
<i>a = aLund + aExtraSQuark</i>. 
   
 
<br/><br/><table><tr><td><strong>StringZ:aExtraDiquark </td><td></td><td> <input type="text" name="4" value="0.97" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.97</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
allows a larger <i>a</i> for diquarks, with total 
<i>a = aLund + aExtraDiquark</i>. 
   
 
<p/> 
Finally, the Bowler modification [<a href="Bibliography.php" target="page">Bow81</a>] introduces an extra 
factor 
<br/><i> 
    1/z^{r_Q * b * m_Q^2} 
</i><br/> 
for heavy quarks. To keep some flexibility, a multiplicative factor 
<i>r_Q</i> is introduced, which ought to be unity (provided that 
quark masses were uniquely defined) but can be set in 
 
<br/><br/><table><tr><td><strong>StringZ:rFactC </td><td></td><td> <input type="text" name="5" value="1.32" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.32</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
<i>r_c</i>, i.e. the above parameter for <i>c</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>StringZ:rFactB </td><td></td><td> <input type="text" name="6" value="0.855" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.855</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
<i>r_b</i>, i.e. the above parameter for <i>b</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>StringZ:rFactH </td><td></td><td> <input type="text" name="7" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
<i>r_h</i>, i.e. the above parameter for heavier hypothetical quarks, 
or in general any new coloured particle long-lived enough to hadronize. 
   
 
<p/> 
Within the string framework, the <i>b</i> parameter is universal, 
i.e. common for all flavours. Nevertheless, for fits to experimental 
data, better agreement can be obtained if both <i>a_Q</i> and 
<i>b_Q</i> can be set freely in a general expression 
<br/><i> 
    f(z) = 1/z^{1 + r_Q * b_Q * m_Q^2} * (1-z)^a_Q * exp(-b_Q m_T^2 / z) 
</i><br/> 
The below switches and values can be used to achieve this. They should 
be used with caution and constitute clear deviations from the Lund 
philosophy. 
 
<br/><br/><strong>StringZ:useNonstandardC</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
use the above nonstandard Lund ansatz for <i>c</i> quarks. 
   
 
<br/><br/><strong>StringZ:useNonstandardB</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
use the above nonstandard Lund ansatz for <i>b</i> quarks. 
   
 
<br/><br/><strong>StringZ:useNonstandardH</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
use the above nonstandard Lund ansatz for hypothetical heavier quarks. 
   
 
<br/><br/><table><tr><td><strong>StringZ:aNonstandardC </td><td></td><td> <input type="text" name="11" value="0.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.3</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
The <i>a</i> parameter in the nonstandard Lund ansatz for 
 <i>c</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>StringZ:aNonstandardB </td><td></td><td> <input type="text" name="12" value="0.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.3</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
The <i>a</i> parameter in the nonstandard Lund ansatz for 
 <i>b</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>StringZ:aNonstandardH </td><td></td><td> <input type="text" name="13" value="0.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.3</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
The <i>a</i> parameter in the nonstandard Lund ansatz for 
hypothetical heavier quarks. 
   
 
<br/><br/><table><tr><td><strong>StringZ:bNonstandardC </td><td></td><td> <input type="text" name="14" value="0.8" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.8</strong></code>; <code>minimum = 0.2</code>; <code>maximum = 2.0</code>)</td></tr></table>
The <i>b</i> parameter in the nonstandard Lund ansatz for 
<i>c</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>StringZ:bNonstandardB </td><td></td><td> <input type="text" name="15" value="0.8" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.8</strong></code>; <code>minimum = 0.2</code>; <code>maximum = 2.0</code>)</td></tr></table>
The <i>b</i> parameter in the nonstandard Lund ansatz for 
<i>b</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>StringZ:bNonstandardH </td><td></td><td> <input type="text" name="16" value="0.8" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.8</strong></code>; <code>minimum = 0.2</code>; <code>maximum = 2.0</code>)</td></tr></table>
The <i>b</i> parameter in the nonstandard Lund ansatz for 
hypothetical heavier quarks. 
   
 
<p/> 
As another nonstandard alternative, it is possible to switch over to the 
Peterson/SLAC formula [<a href="Bibliography.php" target="page">Pet83</a>] 
<br/><i> 
     f(z) = 1 / ( z * (1 - 1/z - epsilon/(1-z))^2 ) 
</i><br/> 
for charm, bottom and heavier (defined as above) by the three flags 
 
<br/><br/><strong>StringZ:usePetersonC</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
use Peterson for <i>c</i> quarks. 
   
 
<br/><br/><strong>StringZ:usePetersonB</strong>  <input type="radio" name="18" value="on"><strong>On</strong>
<input type="radio" name="18" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
use Peterson for <i>b</i> quarks. 
   
 
<br/><br/><strong>StringZ:usePetersonH</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
use Peterson for hypothetical heavier quarks. 
   
 
<p/> 
When switched on, the corresponding epsilon values are chosen to be 
 
<br/><br/><table><tr><td><strong>StringZ:epsilonC </td><td></td><td> <input type="text" name="20" value="0.05" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.05</strong></code>; <code>minimum = 0.01</code>; <code>maximum = 0.25</code>)</td></tr></table>
<i>epsilon_c</i>, i.e. the above parameter for <i>c</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>StringZ:epsilonB </td><td></td><td> <input type="text" name="21" value="0.005" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.005</strong></code>; <code>minimum = 0.001</code>; <code>maximum = 0.025</code>)</td></tr></table>
<i>epsilon_b</i>, i.e. the above parameter for <i>b</i> quarks. 
   
 
<br/><br/><table><tr><td><strong>StringZ:epsilonH </td><td></td><td> <input type="text" name="22" value="0.005" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.005</strong></code>; <code>minimum = 0.0001</code>; <code>maximum = 0.25</code>)</td></tr></table>
<i>epsilon_h</i>, i.e. the above parameter for hypothetical heavier 
quarks, normalized to the case where <i>m_h = m_b</i>. The actually 
used parameter is then <i>epsilon = epsilon_h * (m_b^2 / m_h^2)</i>. 
This allows a sensible scaling to a particle with an unknown higher 
mass without the need for a user intervention. 
   
 
<h3>Fragmentation <i>pT</i></h3> 
 
The <code>StringPT</code> class handles the choice of fragmentation 
<i>pT</i>. At each string breaking the quark and antiquark of the pair are 
supposed to receive opposite and compensating <i>pT</i> kicks according 
to a Gaussian distribution in <i>p_x</i> and <i>p_y</i> separately. 
Call <i>sigma_q</i> the width of the <i>p_x</i> and <i>p_y</i> 
distributions separately, i.e. 
<br/><i> 
    d(Prob) = exp( -(p_x^2 + p_y^2) / 2 sigma_q^2). 
</i><br/> 
Then the total squared width is 
<br/><i> 
    &lt;pT^2> = &lt;p_x^2> +  &lt;p_y^2> = 2 sigma_q^2 = sigma^2. 
</i><br/> 
It is this latter number that is stored in 
 
<br/><br/><table><tr><td><strong>StringPT:sigma </td><td></td><td> <input type="text" name="23" value="0.335" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.335</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
the width <i>sigma</i> in the fragmentation process. 
   
 
<p/> 
Since a normal hadron receives <i>pT</i> contributions for two string 
breakings, it has a <i>&lt;p_x^2>_had = &lt;p_y^2>_had = sigma^2</i>, 
and thus <i>&lt;pT^2>_had = 2 sigma^2</i>. 
 
<p/> 
Some studies on isolated particles at LEP has indicated the need for 
a slightly enhanced rate in the high-<i>pT</i> tail of the above 
distribution. This would have to be reviewed in the context of a 
complete retune of parton showers and hadronization, but for the 
moment we stay with the current recipe, to boost the above <i>pT</i> 
by a factor <i>enhancedWidth</i> for a small fraction 
<i>enhancedFraction</i> of the breakups, where 
 
<br/><br/><table><tr><td><strong>StringPT:enhancedFraction </td><td></td><td> <input type="text" name="24" value="0.01" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.</code>)</td></tr></table>
<i>enhancedFraction</i>,the fraction of string breaks with enhanced 
width. 
   
 
<br/><br/><table><tr><td><strong>StringPT:enhancedWidth </td><td></td><td> <input type="text" name="25" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 1.0</code>; <code>maximum = 10.0</code>)</td></tr></table>
<i>enhancedWidth</i>,the enhancement of the width in this fraction. 
   
 
<h3>Jet joining procedure</h3> 
 
String fragmentation is carried out iteratively from both string ends 
inwards, which means that the two chains of hadrons have to be joined up 
somewhere in the middle of the event. This joining is described by 
parameters that in principle follows from the standard fragmentation 
parameters, but in a way too complicated to parametrize. The dependence 
is rather mild, however, so for a sensible range of variation the 
parameters in this section should not be touched. 
 
<br/><br/><table><tr><td><strong>StringFragmentation:stopMass </td><td></td><td> <input type="text" name="26" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
Is used to define a <i>W_min = m_q1 + m_q2 + stopMass</i>, 
where <i>m_q1</i> and <i>m_q2</i> are the masses of the two 
current endpoint quarks or diquarks. 
   
 
<br/><br/><table><tr><td><strong>StringFragmentation:stopNewFlav </td><td></td><td> <input type="text" name="27" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 2.0</code>)</td></tr></table>
Add to <i>W_min</i> an amount <i>stopNewFlav * m_q_last</i>, 
where <i>q_last</i> is the last <i>q qbar</i> pair produced 
between the final two hadrons. 
   
 
<br/><br/><table><tr><td><strong>StringFragmentation:stopSmear </td><td></td><td> <input type="text" name="28" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 0.5</code>)</td></tr></table>
The <i>W_min</i> above is then smeared uniformly in the range 
<i>W_min_smeared = W_min * [ 1 - stopSmear, 1 + stopSmear ]</i>. 
   
 
<p/> 
This <i>W_min_smeared</i> is then compared with the current remaining 
<i>W_transverse</i> to determine if there is energy left for further 
particle production. If not, i.e. if 
<i>W_transverse &lt; W_min_smeared</i>, the final two particles are 
produced from what is currently left, if possible. (If not, the 
fragmentation process is started over.) 
 
<h3>Simplifying systems</h3> 
 
There are a few situations when it is meaningful to simplify the 
original task, one way or another. 
 
<br/><br/><table><tr><td><strong>HadronLevel:mStringMin </td><td></td><td> <input type="text" name="29" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.5</code>; <code>maximum = 1.5</code>)</td></tr></table>
Decides whether a partonic system should be considered as a normal 
string or a ministring, the latter only producing one or two primary 
hadrons. The system mass should be above <i>mStringMin</i> plus the 
sum of quark/diquark constituent masses for a normal string description, 
else the ministring scenario is used. 
   
 
<br/><br/><table><tr><td><strong>FragmentationSystems:mJoin </td><td></td><td> <input type="text" name="30" value="0.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.3</strong></code>; <code>minimum = 0.2</code>; <code>maximum = 1.</code>)</td></tr></table>
When two colour-connected partons are very nearby, with at least 
one being a gluon, they can be joined into one, to avoid technical 
problems of very small string regions. The requirement for joining is 
that the invariant mass of the pair is below <i>mJoin</i>, where a 
gluon only counts with half its momentum, i.e. with its contribution 
to the string region under consideration. (Note that, for technical 
reasons, the 0.2 GeV lower limit is de facto hardcoded.) 
   
 
<br/><br/><table><tr><td><strong>FragmentationSystems:mJoinJunction </td><td></td><td> <input type="text" name="31" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.5</code>; <code>maximum = 2.</code>)</td></tr></table>
When the invariant mass of two of the quarks in a three-quark junction 
string system becomes too small, the system is simplified to a 
quark-diquark simple string. The requirement for this simplification 
is that the diquark mass, minus the two quark masses, falls below 
<i>mJoinJunction</i>. Gluons on the string between the junction and 
the respective quark, if any, are counted as part of the quark 
four-momentum. Those on the two combined legs are clustered with the 
diquark when it is formed. 
   
 
<h3>Ministrings</h3> 
 
The <code>MiniStringFragmentation</code> machinery is only used when a 
string system has so small invariant mass that normal string fragmentation 
is difficult/impossible. Instead one or two particles are produced, 
in the former case shuffling energy-momentum relative to another 
colour singlet system in the event, while preserving the invariant 
mass of that system. With one exception parameters are the same as 
defined for normal string fragmentation, to the extent that they are 
at all applicable in this case. 
 
A discussion of the relevant physics is found in [<a href="Bibliography.php" target="page">Nor00</a>]. 
The current implementation does not completely abide to the scheme 
presented there, however, but has in part been simplified. (In part 
for greater clarity, in part since the class is not quite finished yet.) 
 
<br/><br/><table><tr><td><strong>MiniStringFragmentation:nTry  </td><td></td><td> <input type="text" name="32" value="2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>; <code>maximum = 10</code>)</td></tr></table>
Whenever the machinery is called, first this many attempts are made 
to pick two hadrons that the system fragments to. If the hadrons are 
too massive the attempt will fail, but a new subsequent try could 
involve other flavour and hadrons and thus still succeed. 
After <i>nTry</i> attempts, instead an attempt is made to produce a 
single hadron from the system. Should also this fail, some further 
attempts at obtaining two hadrons will be made before eventually 
giving up. 
   
 
<h3>Junction treatment</h3> 
 
A junction topology corresponds to an Y arrangement of strings 
i.e. where three string pieces have to be joined up in a junction. 
Such topologies can arise if several valence quarks are kicked out 
from a proton beam, or in baryon-number-violating SUSY decays. 
Special attention is necessary to handle the region just around 
the junction, where the baryon number topologically is located. 
The junction fragmentation scheme is described in [<a href="Bibliography.php" target="page">Sjo03</a>]. 
The parameters in this section should not be touched except by experts. 
 
<br/><br/><table><tr><td><strong>StringFragmentation:eNormJunction </td><td></td><td> <input type="text" name="33" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.5</code>; <code>maximum = 10</code>)</td></tr></table>
Used to find the effective rest frame of the junction, which is 
complicated when the three string legs may contain additional 
gluons between the junction and the endpoint. To this end, 
a pull is defined as a weighed sum of the momenta on each leg, 
where the weight is <i>exp(- eSum / eNormJunction)</i>, with 
<i>eSum</i> the summed energy of all partons closer to the junction 
than the currently considered one (in the junction rest frame). 
Should in principle be (close to) <i>sqrt((1 + a) / b)</i>, with 
<i>a</i> and <i>b</i> the parameters of the Lund symmetric 
fragmentation function. 
   
 
<br/><br/><table><tr><td><strong>StringFragmentation:eBothLeftJunction </td><td></td><td> <input type="text" name="34" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.5</code>)</td></tr></table>
Retry (up to 10 times) when the first two considered strings in to a 
junction both have a remaining energy (in the junction rest frame) 
above this number. 
   
 
<br/><br/><table><tr><td><strong>StringFragmentation:eMaxLeftJunction </td><td></td><td> <input type="text" name="35" value="10.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10.0</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Retry (up to 10 times) when the first two considered strings in to a 
junction has a highest remaining energy (in the junction rest frame) 
above a random energy evenly distributed between 
<i>eBothLeftJunction</i> and 
<i>eBothLeftJunction + eMaxLeftJunction</i> 
(drawn anew for each test). 
   
 
<br/><br/><table><tr><td><strong>StringFragmentation:eMinLeftJunction </td><td></td><td> <input type="text" name="36" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Retry (up to 10 times) when the invariant mass-squared of the final leg 
and the leftover momentum of the first two treated legs falls below 
<i>eMinLeftJunction</i> times the energy of the final leg (in the 
junction rest frame). 
   
 
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

if($_POST["1"] != "0.68")
{
$data = "StringZ:aLund = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0.98")
{
$data = "StringZ:bLund = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "0.0")
{
$data = "StringZ:aExtraSQuark = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0.97")
{
$data = "StringZ:aExtraDiquark = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1.32")
{
$data = "StringZ:rFactC = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.855")
{
$data = "StringZ:rFactB = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "1.0")
{
$data = "StringZ:rFactH = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "StringZ:useNonstandardC = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "StringZ:useNonstandardB = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "StringZ:useNonstandardH = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "0.3")
{
$data = "StringZ:aNonstandardC = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "0.3")
{
$data = "StringZ:aNonstandardB = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0.3")
{
$data = "StringZ:aNonstandardH = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "0.8")
{
$data = "StringZ:bNonstandardC = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "0.8")
{
$data = "StringZ:bNonstandardB = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "0.8")
{
$data = "StringZ:bNonstandardH = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "StringZ:usePetersonC = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "off")
{
$data = "StringZ:usePetersonB = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "StringZ:usePetersonH = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "0.05")
{
$data = "StringZ:epsilonC = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "0.005")
{
$data = "StringZ:epsilonB = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "0.005")
{
$data = "StringZ:epsilonH = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "0.335")
{
$data = "StringPT:sigma = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "0.01")
{
$data = "StringPT:enhancedFraction = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "2.0")
{
$data = "StringPT:enhancedWidth = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "1.0")
{
$data = "StringFragmentation:stopMass = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "2.0")
{
$data = "StringFragmentation:stopNewFlav = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "0.2")
{
$data = "StringFragmentation:stopSmear = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "1.")
{
$data = "HadronLevel:mStringMin = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "0.3")
{
$data = "FragmentationSystems:mJoin = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "1.0")
{
$data = "FragmentationSystems:mJoinJunction = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "2")
{
$data = "MiniStringFragmentation:nTry = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "2.0")
{
$data = "StringFragmentation:eNormJunction = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "1.0")
{
$data = "StringFragmentation:eBothLeftJunction = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "10.0")
{
$data = "StringFragmentation:eMaxLeftJunction = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "0.2")
{
$data = "StringFragmentation:eMinLeftJunction = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
