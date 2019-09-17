<html>
<head>
<title>Colour reconnection</title>
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

<form method='post' action='ColourReconnection.php'>
 
<h2>Colour Reconnection</h2> 
<ol id="toc">
  <li><a href="#section0">The MPI-based scheme</a></li>
  <li><a href="#section1">The newer scheme</a></li>
  <li><a href="#section2">The gluon-move scheme</a></li>
  <li><a href="#section3">The <ei>e^+ e^-</ei> colour reconnection schemes</a></li>
</ol>

 
The colour flows in the separate subprocesses defined in the 
multiparton-interactions scenario are tied together via the assignment 
of colour flow in the beam remnant. This is not an unambiguous 
procedure, and currently two different methods are implemented. In the first 
model the colour flow is reconstructed by how a PS could have 
constructed the configuration. In the second model, the full QCD colour 
calculation is taken into account, however the dynamical effects are modeled 
loosely, only an overall saturation is taken into account. The idea is to 
later account for other dynamical effects through colour reconnections. 
 
<p/> 
A simple "minimal" procedure of colour flow only via the beam remnants 
does not result in a scenario in agreement with data, however, 
notably not a sufficiently steep rise of 
<i>&lt;pT&gt;(n_ch)</i>. The true origin of this behaviour and the 
correct mechanism to reproduce it remains one of the big unsolved issues 
at the borderline between perturbative and nonperturbative QCD. Since no final 
answer is known, several models are implemented. The different models also 
rely on the two different colour assignments in the beam remnant. There are 
two, somewhat motivated, models implemented: the original PYTHIA scheme and 
a new scheme that tries to incorporate more of the colour knowledge from QCD. 
 
<p/> 
The original PYTHIA scheme relies on the PS-like colour configuration of the 
beam remnant. This is combined with an additional step, wherein the gluons 
of a lower-<i>pT</i> MPI system are merged with the ones in a higher-pT MPI. 
A more detailed description of the merging can be found below. 
Relative to the other models it tests fewer reconnection possibilities, 
and therefore tends to be reasonably fast. 
 
<p/> 
The new scheme [<a href="Bibliography.php#refChr14a" target="page">Chr14a</a>]relies on the full QCD colour configuration 
in the beam remnant. This is followed up by a colour reconnection, where the 
potential string energy is minimized (ie. the <i>lambda</i> measure is 
minimized). The QCD colour rules are also incorporated in the colour 
reconnection, and determine the probability that a reconnection is allowed. 
The model also allows the creation of junction structures. 
 
<p/> 
In addition to the two models described above, a simple model is implemented, 
wherein gluons can be moved from one location to another so as to reduce the 
total string length. This is one out of a range of simple models developed 
to study potential colour reconnection effects e.g. on top mass 
[<a href="Bibliography.php#refArg14" target="page">Arg14</a>], not from the point of view of having the most realistic 
description, but in order to probe the potential worst-case spread of 
predictions. All of these models are made available separately in 
<code>include/Pythia8Plugins/ColourReconnectionHooks.h</code>, with the 
setup illustrated in <code>examples/main29.cc</code>, but only the 
gluon-move one is sufficiently general and realistic that it has been 
included among the standard options here. 
 
<p/> 
Finally, the SK I and SK II models [<a href="Bibliography.php#refSjo94" target="page">Sjo94</a>] have a smaller range 
of applicability, originally intended for <i>e^+ e^- &rarr; W^+ W^-</i>, 
but in this context offers a more detailed approach. They are not suitable 
for hadronic collisions, since they would only address CR inside a gauge 
boson pair, and not CR in the rest of the event. 
 
<br/><br/><strong>ColourReconnection:reconnect</strong>  <input type="radio" name="1" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="1" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow or not a system to be merged with another one. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:mode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 4</code>)</td></tr></table>
Determine which model is used for colour reconnection. Beware that 
some settings may need to be changed to match the model selected. 
<br/>
<input type="radio" name="2" value="0" checked="checked"><strong>0 </strong>:  The MPI-based original Pythia 8 scheme.  <br/>
<input type="radio" name="2" value="1"><strong>1 </strong>:  The new more QCD based scheme. Should be combined ewith  <code>BeamRemnants:remnantMode = 1.</code>  <br/>
<input type="radio" name="2" value="2"><strong>2 </strong>:  The new gluon-move model.  <br/>
<input type="radio" name="2" value="3"><strong>3 </strong>:  The SK I <ei>e^+ e^-</ei> CR model. Requires  <code>ColourReconnection:forceResonance = on</code> (and default  <code>PartonLevel:earlyResDec = off</code>) to give any CR.  <br/>
<input type="radio" name="2" value="4"><strong>4 </strong>:  The SK II <ei>e^+ e^-</ei> CR model. Requires  <code>ColourReconnection:forceResonance = on</code> (and default  <code>PartonLevel:earlyResDec = off</code>) to give any CR.  <br/>
 
<br/><br/><strong>ColourReconnection:forceHadronLevelCR</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
This flag switches on colour reconnection in the <code>forceHadronLevel</code> 
function. The function is called when only the hadron level of PYTHIA is 
used (see <?php $filepath = $_GET["filepath"];
echo "<a href='HadronLevelStandalone.php?filepath=".$filepath."' target='page'>";?>Hadron-level 
Standalone</a>). The MPI-based model is not available for this setup and 
any resonance decays not already decayed are not included in the CR. 
   
 
<br/><br/><strong>ColourReconnection:forceResonance</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
This parameter allows an additional CR after late resonance decays. All the 
particles from all resonance decays are allowed to reconnect with each 
other. It is mainly intended for <i>H -> WW -> qqqq </i>, where the 
Higgs decay ensures a separation between the W bosons and the MPI 
systems. Reconnections between the decay products from the two <i>W</i> 
bosons is still a possibility, however. This option is not available for 
colored resonances, and not for the MPI-based model. 
   
 
<a name="section0"></a> 
<h3>The MPI-based scheme</h3> 
 
In this scheme partons are classified by which MPI system they belong to. 
The colour flow of two such systems can be fused, and if so the partons 
of the lower-<i>pT</i> system are added to the strings defined by the 
higher-<i>pT</i> system in such a way as to give the smallest total 
string length. The bulk of these lower-<i>pT</i> partons are gluons, 
and this is what the scheme is optimized to handle. 
 
<p/> 
In  detail, an MPI system with a scale <i>pT</i> of the hard 
interaction (normally <i>2 &rarr; 2</i>) can be merged with one of 
a harder scale with a probability that is 
<i>pT0_Rec^2 / (pT0_Rec^2 + pT^2)</i>, where <i>pT0_Rec</i> is 
<code>range</code> times <i>pT0</i>, the latter being the same 
energy-dependent dampening parameter as used for MPIs. 
Thus it is easy to merge a low-<i>pT</i> system with any other, 
but difficult to merge two high-<i>pT</i> ones with each other. 
 
<br/><br/><table><tr><td><strong>ColourReconnection:range </td><td></td><td> <input type="text" name="5" value="1.8" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.8</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
The <code>range</code> parameter defined above. The higher this number is 
the more reconnections can occur. For values above unity the reconnection 
rate tends to saturate, since then most systems are already connected with 
each other. This is why 10 is set as an effective upper limit, beyond 
which it is not meaningful to let the parameter go. 
   
 
<p/> 
The reconnection procedure is applied iteratively. Thus first the 
reconnection probability <i>P = pT0_Rec^2 / (pT0_Rec^2 + pT^2)</i> of the 
lowest-<i>pT</i> system is found, and gives the probability for merger with 
the second-lowest one. If not merged, it is tested with the third-lowest 
one, and so on. For the <i>m</i>'th higher system the reconnection 
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
the one introduced already in [<a href="Bibliography.php#refSjo87" target="page">Sjo87</a>]. Clearly it is ad hoc. 
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
considered by default. Colour inside a resonance therefore is not 
reconnected. The 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='MasterSwitches.php?filepath=".$filepath."' target='page'>";?>PartonLevel:earlyResDec</a></code> 
can be switched on to perform resonance decays before colour reconnection, 
and then the partons from resonance decays could be affected. 
Ideally the time scales of resonance decays and of colour reconnection 
should be picked dynamically, but this is not yet the case. Notably 
the <i>W</i>, <i>Z</i> and <i>t</i> have intermediate decay time 
scales, somewhat but not much shorter than typical hadronization times, 
whereas the <i>H</i> is much  long-lived. 
 
<a name="section1"></a> 
<h3>The newer scheme</h3> 
 
The newer CR scheme builds on the minimization of the string length as well as 
the colour rules from QCD. A main feature of the new model is the introduction 
of junction structures. These are possible outcomes of the reconnection in 
addition to the more common string-string reconnections. The model works by 
constructing all pair of dipoles that are allowed to reconnect by QCD colour 
rules and switching if the new pair has a lower string length. Junctions are 
also allowed to be directly produced from three, and in some special cases, 
four dipoles. This is done iteratively until no further allowed reconnection 
lowers the total string length. 
 
<p/> 
According to QCD colour rules, an uncorrelated triplet and anti-triplet are 
allowed to form a singlet state <i>1/9</i> times. This is reflected in the 
model by giving each dipole a colour number between 0-8 and only dipoles with 
the same colour number are allowed to reconnect. The junction probability is 
given by the product of two triplets, which provides an anti-triplet 
<i>1/3</i> times. This is achieved in the model by allowing reconnections 
between dipoles where modulo three of the color numbers agree. In addition to 
the colour rules, the dipoles also need to be causally connected in order to 
perform a reconnection. The definition of causally connected dipoles is not 
exact, and several different options are available. All the time dilation 
modes introduce a tuneable parameter, which provides a handle on the overall 
amount of colour reconnection. 
 
<p/> 
When the two strings are allowed to reconnect, they will reconnect if it 
lowers the total string length. The total string length is in the model 
defined by an approximation to the <i>lambda</i>-measure. Several options 
for different approximations are available. The <i>lambda</i>-measure is not 
well understood, especially for junction structures, and a tuneable parameter 
is introduced to vary the behaviour between junctions and ordinary strings. 
 
<p/> 
To avoid problems with very low mass string and junction structures, these are 
excluded from participating in the colour reconnections. This is achieved by 
forming the dipole or junction into a pseudo-particle if the invariant mass 
is too low. Especially the approximations made in the <i>lambda</i>-measure 
provides problems at low invariant masses. 
 
<p/> 
The new CR scheme introduce several tuneable parameters, which 
all are listed below. In addition to these, other parameters in PYTHIA also 
need to retuned to account for the new CR. The default values below, together 
with changing <code> MultipartonInteractions:pT0Ref = 2.15</code> and 
<code>ColourReconnection:allowDoubleJunRem = off</code>, provides a good 
starting point. Additional fragmentation variables were also adjusted in the 
first tune, but these provide a smaller change (see [<a href="Bibliography.php#refChr14a" target="page">Chr14a</a>] for a 
complete list). 
 
 
<br/><br/><table><tr><td><strong>ColourReconnection:m0 </td><td></td><td> <input type="text" name="6" value="0.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.3</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 5.</code>)</td></tr></table>
This is the variable used in the lambda measure for the string length. 
See the different choices of lambda measure for exact formulas. This variable 
is in the new model also used as a cut for forming pseudo particles that are 
not colour reconnected. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:junctionCorrection </td><td></td><td> <input type="text" name="7" value="1.20" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.20</strong></code>; <code>minimum = 0.01</code>; <code>maximum = 10.</code>)</td></tr></table>
This variable allows to use a different m0 for junction strings in the lambda 
measure. It is multiplicative correction to the m0 chosen above. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:nColours  </td><td></td><td> <input type="text" name="8" value="9" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>9</strong></code>; <code>minimum = 1</code>; <code>maximum = 30</code>)</td></tr></table>
The number of reconnection colours, this should not be confused with the 
standard number of QCD colours. Each string is given an integer number between 
0 and <code>nColours - 1</code>. Only strings with the same number are allowed 
to do a normal string reconnection. The default value provides 
the standard QCD probability that a triplet and an anti-triplet is in a 
singlet state. The probability for two strings to form a junction structure is 
in QCD given by the product of two triplets, which gives a probability of 1/3. 
Therefore the number of reconnection colours for junction formation is 
<code>iColours % 3</code>, where iColours refer to the integer of the string. 
The behaviour of junction formation therefore only changes slightly with this 
variable. 
   
 
<br/><br/><strong>ColourReconnection:sameNeighbourColours</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
In the normal colour reconnection two neighbouring strings are not allowed 
to have the same colour. Similar two strings originating from a gluon split is 
not allowed to reconnect. The physics motivation for this is that it would 
require colour singlet gluons, and therefore for ordinary physics studies this 
should be turned off. But for testing of extreme scenarios (i.e. 1 colour), 
this variable needs to be turned on, since it is not possible to have 
different neighbouring colours. 
   
 
<br/><br/><strong>ColourReconnection:allowJunctions</strong>  <input type="radio" name="10" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="10" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
This switch disables the formation of junctions in the colour reconnection. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:lambdaForm  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
This allows to switch between different options for what 
<ei>lambda</ei>-measure to use. 
The formula shown are how much each end of a dipole or junction contribute to 
the total <ei>lambda</ei>-measure. The energies are defined in respectively 
the dipole or junction rest frame. 
<br/>
<input type="radio" name="11" value="0" checked="checked"><strong>0 </strong>:  <ei>lambda = ln (1 + sqrt(2) E/m0)</ei>  <br/>
<input type="radio" name="11" value="1"><strong>1 </strong>:  <ei>lambda = ln (1 + 2 E/m0)</ei>  <br/>
<input type="radio" name="11" value="2"><strong>2 </strong>:  <ei>lambda = ln (2 E/m0)</ei>  <br/>
 
<br/><br/><strong>ColourReconnection:allowDoubleJunRem</strong>  <input type="radio" name="12" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="12" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
This parameter tells whether or not to allow a directly connected 
junction-antijunction pair to split into two strings. The lambda measure of 
the junction system is compared to that of the two possible string 
configurations. If the chosen configuration is the junction system, a q-qbar 
system is inserted between the junctions by removing some energy/momentum from 
the other legs. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:timeDilationMode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
Disallow colour reconnection between strings that are not in causal 
contact; if either string has already decayed before the other string forms, 
there is no space-time region in which the reconnection could physically 
occur. The exact defintion of causal contact is not known, hence several 
possible defintions are included. They all include the boost factor, 
<ei>gamma</ei>, and the majority also rely on the typical hadronization scale, 
<ei>r</ei>, which is kept fixed at 1 fm. A tuneable dimensionless parameter is 
included, which can be used to control the overall amount of colour 
reconnection. 
<br/>
<input type="radio" name="13" value="0"><strong>0 </strong>:  All strings are allowed to reconnect.  <br/>
<input type="radio" name="13" value="1"><strong>1 </strong>:  Strings are allowed to reconnect if <ei>gamma &lt  timeDilationPar </ei> and all strings should be causally connected to allow a  reconnection.  <br/>
<input type="radio" name="13" value="2" checked="checked"><strong>2 </strong>:  Strings are allowed to reconnect if <ei>gamma &lt  timeDilationPar * mDip * r </ei> and all strings should be in causal  contact to allow a reconnection.  <br/>
<input type="radio" name="13" value="3"><strong>3 </strong>:  Strings are allowed to reconnect if <ei>gamma &lt  timeDilationPar * mDip * r </ei> and if a single pair of dipoles are in  causal contact the reconnection is allowed.  <br/>
<input type="radio" name="13" value="4"><strong>4 </strong>:  Strings are allowed to reconnect if <ei>gamma &lt  timeDilationPar * mDip' * r </ei> and all strings should be in causal  contact to allow a reconnection. mDip' is the invariant mass at the  formation of the dipole (ie. the first time the colour tag appear in the  perturbative expansion).  <br/>
<input type="radio" name="13" value="5"><strong>5 </strong>:  Strings are allowed to reconnect if <ei>gamma &lt  timeDilationPar * mDip' * r </ei> and if a single pair of dipoles are in  causal contact the reconnection is allowed. mDip' is the invariant mass at  the formation of the dipole (ie. the first time the colour tag appear in  the perturbative expansion).  <br/>
 
<br/><br/><table><tr><td><strong>ColourReconnection:timeDilationPar </td><td></td><td> <input type="text" name="14" value="0.18" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.18</strong></code>; <code>minimum = 0</code>; <code>maximum = 100</code>)</td></tr></table>
This is a tuneable parameter for the time dilation. The definition 
can be seen above under <code>timeDilationMode</code>. 
   
 
<a name="section2"></a> 
<h3>The gluon-move scheme</h3> 
 
This approach contains two steps, a first "move" one and an optional 
second "flip" one. Both are intended to reduce the total "string length" 
<i>lambda</i> measure of an event. For multiparton topologies the 
correct <i>lambda</i> measure can become quite cumbersome, so here it 
is approximated by the sum of <i>lambda</i> contributions from each pair 
of partons connected by a colour string piece. For two partons <i>i</i> 
and <i>j</i> with invariant mass <i>m_ij</i> this contribution 
is defined as <i>lambda_ij = ln(1 + m^2_ij / m2Lambda)</i>. 
The 1 is added ad hoc to avoid problems in the <i>m_ij &rarr; 0</i> 
limit, problems which mainly comes from the approximate treatment, 
and <i>m2Lambda</i> is a parameter set below. 
 
<p/> 
In the move step all final gluons are identified, alternatively only a 
fraction <i>fracGluon</i> of them, and also all colour-connected 
parton pairs. For each gluon and each pair it is calculated how the total 
<i>lambda</i> would shift if the gluon would be removed from its current 
location and inserted in between the pair. The gluon move that gives the 
largest negative shift, if any, is then carried out. Alternatively, only 
shifts  negative than <i>dLambdaCut</i> are considered. The procedure 
is iterated so long as allowed moves can be found. 
 
<p/> 
There is some fine print. If a colour singlet subsystem consists of two 
gluons only then it is not allowed to move any of them, since that then 
would result in in a colour singlet gluon. Also, at most as many moves 
are made as there are gluons, which normally should be enough. A specific 
gluon may be moved more than once, however. Finally, a gluon directly 
connected to a junction cannot be moved, and also no gluon can be inserted 
between it and the junction. This is entirely for practical reasons, but 
should not be a problem, since junctions are rare in this model. 
 
<p/> 
The gluon-move steps will not break the connection between string endpoints, 
in the sense that a quark and an antiquark that are colour-connected via 
a number of gluons will remain so, only the number and identity of the 
intermediate gluons may change. Such a scenario may be too restrictive. 
Therefore an optional second flip step is introduced. In it all such 
colour chains are identified, omitting closed gluon loops. The lambda 
change is defined by what happens if the two colour lines are crossed 
somewhere, e.g. such that two systems <i>q1 - g1 - qbar1</i> and 
<i>q2 - g2 - qbar2</i> are flipped to <i>q1 - g1 - g2 - qbar2</i> 
and <i>q2 - qbar1</i>. The flip that gives the largest <i>lambda</i> 
reduction is carried out, again with <i>dLambdaCut</i> offering a 
possibility to restrict the options. As with the move step, the procedure 
is repeated so long as it is allowed. An important restriction is 
imposed, however, that a given system is only allowed to flip once, 
and not with itself. The practical reason is that repeated flips could 
split off closed gluon loops quite easily, which tends to result in 
bad agreement with data. 
 
<p/> 
As an option, singlet subsystems containing a junction may or may not 
be allowed to take part in the flip step. Since the number of junction 
systems is limited in this model the differences are not so important. 
 
<br/><br/><table><tr><td><strong>ColourReconnection:m2Lambda </td><td></td><td> <input type="text" name="15" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 16.</code>)</td></tr></table>
The <i>m2Lambda</i> parameter used in the definition of the approximate 
<i>lambda</i> measure above. It represents an approximate hadronic 
mass-square scale, cf. <i>m0</i> in the previous model. Its value is 
uncertain up to factors of 2, but the <i>lambda</i> change induced by 
a potential move or flip is rather insensitive to the precise value, 
owing to large cancellations. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:fracGluon </td><td></td><td> <input type="text" name="16" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
The probability that a given gluon will be considered for being moved. 
It thus gives the average fraction of gluons being considered. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:dLambdaCut </td><td></td><td> <input type="text" name="17" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
Restrict gluon moves and colour flips to those that reduce <i>lambda</i> 
by more than this amount. The larger this number, the fewer moves and flips 
will be performed, but those that remain are the ones most likely to produce 
large effects. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:flipMode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 4</code>)</td></tr></table>
Performing the flip step or not. Also possibility to omit the move step. 
<br/>
<input type="radio" name="18" value="0" checked="checked"><strong>0 </strong>:  No flip handling.  <br/>
<input type="radio" name="18" value="1"><strong>1 </strong>:  Allow flips, but not for strings in junction topologies.  <br/>
<input type="radio" name="18" value="2"><strong>2 </strong>:  Allow flips, including for strings in junction topologies.  <br/>
<input type="radio" name="18" value="3"><strong>3 </strong>:  No move handling. Allow flips, but not for strings  in junction topologies.  <br/>
<input type="radio" name="18" value="4"><strong>4 </strong>:  No move handling. Allow flips, including for strings  in junction topologies.  <br/>
 
<a name="section3"></a> 
<h3>The <i>e^+ e^-</i> colour reconnection schemes</h3> 
 
The SK I and SK II models [<a href="Bibliography.php#refSjo94" target="page">Sjo94</a>] were specifically developed for 
<i>e^+ e^- &rarr; W^+ W^- &rarr; q_1 qbar_2 q_3 qbar_4</i> at LEP 2, 
and equally well works for <i>e^+ e^- &rarr; gamma^*/Z^0 gamma^*/Z^0</i>. 
They are not intended to handle hadronic collisions, except in special 
contexts. The prime of these is Higgs decays of the same character as above, 
<i>H^0 &rarr;  W^+ W^- / Z^0 Z^0</i>, since the Higgs is sufficiently 
long-lived that its decay products can be considered separately from 
the rest of the event. The administrative machinery for this possibility 
is not yet in place, however. 
 
<p/> 
The labels I and II refer to the colour-confinement strings being modelled 
either by analogy with type I or type II superconductors. In the former 
model the strings are viewed as transversely extended "bags". The 
likelihood of reconnection is then related to the integrated space-time 
overlap of string pieces from the <i>W^+</i> with those from the 
<i>W^-</i>. In the latter model instead strings are assumed to be 
analogous with vortex lines, where all the topological information 
is stored in a thin core region. Reconnection therefore only can occur 
when these cores pass through each other. 
 
<p/> 
Both of these models are based on a detailed modelling of the space-time 
separation of the <i>W^+</i> and <i>W^-</i> decay vertices, on the 
subsequent shower evolution, on the continued space-time evolution of all 
the string pieces stretched between the showered partons, and on the 
cutoff provided by the strings disappearing by the hadronization process. 
As such, they are more sophisticated than any other reconnection models. 
Unfortunately they would not easily carry over to hadronic collisions, 
where both the initial and the final states are far more complicated, 
and the space-time details less well controlled. 
 
<p/> 
The SK II model has few free parameters, giving more predictive power. 
Conversely, SK I has a a free overall CR strength parameter, making it 
more convenient for tunes to data. The LEP collaborations have used SK I 
as a common reference to establish the existence of CR in <i>W^+ W^-</i> 
events. 
 
<br/><br/><strong>ColourReconnection:lowerLambdaOnly</strong>  <input type="radio" name="19" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="19" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Only allow overlaps that lowers the total string length. 
   
 
<br/><br/><strong>ColourReconnection:singleReconnection</strong>  <input type="radio" name="20" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="20" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Limit the number of reconnections to a single reconnection. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:kI </td><td></td><td> <input type="text" name="21" value="0.6" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.6</strong></code>; <code>minimum = 0.</code>; <code>maximum = 100.</code>)</td></tr></table>
kI is the main free parameter in the reconnection probability for SK I. 
This probability is given by kI times the space-time overlap volume, 
up to saturation effects. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:fragmentationTime </td><td></td><td> <input type="text" name="22" value="1.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.5</strong></code>; <code>minimum = 1.</code>; <code>maximum = 2.</code>)</td></tr></table>
This parameter specifies the average fragmentation time of the string, 
in fm. This is used as an upper limit on the invariant time where strings 
still exist and thus can collide. The expected fragmentation time is 
related to the <i>a</i> and <i>b</i> parameters of the Lund string 
fragmentation function as well as to the string tension. It is therefore 
not a quite free parameter. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:rHadron </td><td></td><td> <input type="text" name="23" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.3</code>; <code>maximum = 1.</code>)</td></tr></table>
Width of the type I string in reconnection calculations, in fm, giving the 
radius of the Gaussian distribution in <i>x</i> and <i>y</i> separately. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:blowR </td><td></td><td> <input type="text" name="24" value="2.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.5</strong></code>; <code>minimum = 1.</code>; <code>maximum = 4.</code>)</td></tr></table>
Technical parameter used in the Monte Carlo sampling of the spatial 
phase space volume in SK I. There is no real reason to change this number. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:blowT </td><td></td><td> <input type="text" name="25" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 1.</code>; <code>maximum = 4.</code>)</td></tr></table>
Technical parameter used in the Monte Carlo sampling of the temporal 
phase space volume in SK I. There is no real reason to change this number. 
   
 
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
$data = "ColourReconnection:reconnect = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0")
{
$data = "ColourReconnection:mode = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "ColourReconnection:forceHadronLevelCR = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "ColourReconnection:forceResonance = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1.8")
{
$data = "ColourReconnection:range = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.3")
{
$data = "ColourReconnection:m0 = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "1.20")
{
$data = "ColourReconnection:junctionCorrection = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "9")
{
$data = "ColourReconnection:nColours = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "ColourReconnection:sameNeighbourColours = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "on")
{
$data = "ColourReconnection:allowJunctions = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "0")
{
$data = "ColourReconnection:lambdaForm = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "on")
{
$data = "ColourReconnection:allowDoubleJunRem = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "2")
{
$data = "ColourReconnection:timeDilationMode = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "0.18")
{
$data = "ColourReconnection:timeDilationPar = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "1.")
{
$data = "ColourReconnection:m2Lambda = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "1.")
{
$data = "ColourReconnection:fracGluon = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "0.")
{
$data = "ColourReconnection:dLambdaCut = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "0")
{
$data = "ColourReconnection:flipMode = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "on")
{
$data = "ColourReconnection:lowerLambdaOnly = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "on")
{
$data = "ColourReconnection:singleReconnection = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "0.6")
{
$data = "ColourReconnection:kI = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "1.5")
{
$data = "ColourReconnection:fragmentationTime = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "0.5")
{
$data = "ColourReconnection:rHadron = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "2.5")
{
$data = "ColourReconnection:blowR = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "2.0")
{
$data = "ColourReconnection:blowT = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
