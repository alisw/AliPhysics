
<html>
<head>
<title>Jet Matching</title>
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

<form method='post' action='JetMatching.php'>

<h2>Jet Matching</h2>

This manual page describes the parton-jet matching interfaces for 
PYTHIA8. In this approach, usually referred to as MLM matching 
[<a href="Bibliography.php" target="page">Man02, Man07</a>], the final jets after parton-shower evolution
and jet clustering are matched to the original partons. The event is 
accepted if a reasonable match is found, and rejected if not. The 
rejection step in an approximate way introduces a Sudakov form factor 
on to the hard processes. Notably the parton shower should not generate 
an emission that would doublecount hard activity already included in 
the matrix-element description. Within this general ansatz, different 
technical solutions can be adopted. We provide two alternatives, one 
based on the algorithm used in ALPGEN [<a href="Bibliography.php" target="page">Man03</a>], and another on 
the one used in Madgraph [<a href="Bibliography.php" target="page">Alw11</a>], both reimplemented from 
scratch here. The main points of these two algorithms are outlined 
further down on this page. 

<p/>We also allow for two alternative sources of external events, 
one in the ALPGEN native format and one in the Madgraph LHEF-based one.
All four combinations of input format and matching style are
implemented. In the following it is therefore important to keep
the two aspects apart, whenever the ALPGEN and Madgraph labels 
are used.    

<p/>Currently all the files of interest are located in the 
<code>examples/</code> subdirectory:
<ul>
<li><code>JetMatching.h</code> contains the machinery for the 
parton-jet matching, in the two <code>JetMatchingAlpgen</code> 
and <code>JetMatchingMadgraph</code> classes.
</li>
<li><code>GeneratorInput.h</code> contains three classes for the 
reading of ALPGEN event and parameter files, and one for the reading
of Madgraph parameters.
</li>
<li><code>CombineMatchingInput.h</code> contains three classes that 
combine the reading of events with the matching of them.
</li>
<li>
<code>main32.cc, main32.cmnd</code> : a sample main program and card
file showing the usage of the previous files/classes.
</li>
</ul>

<h2>Event input source</h2>

External sources of partons are used in the parton-jet matching 
process. The source of the partons has been separated from the
implementation of the matching algorithm. By default, PYTHIA8 
contains a machinery to process Les Houches Event Files (LHEFs)
as described on the
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches Accord</a>
and <?php $filepath = $_GET["filepath"];
echo "<a href='BeamParameters.php?filepath=".$filepath."' target='page'>";?>Beam Parameters</a>
pages. Madgraph5 adheres to this format, but also contains some further 
non-standardized information that can be used. The parsing of the native 
ALPGEN file format is described on the  
<?php $filepath = $_GET["filepath"];
echo "<a href='AlpgenEventInterface.php?filepath=".$filepath."' target='page'>";?>Alpgen Event Interface</a> page.

<p/>Commonly, the source of external partons also contains information
about how a particular type of matching algorithm should be employed.
This information is handled by the <code>AlpgenPar</code> class for 
ALPGEN files, and <code>MadgraphPar</code> for LHEFs.
The user can choose to set default matching parameters using the
<code><?php $filepath = $_GET["filepath"];
echo "<a href='AlpgenEventInterface.php?filepath=".$filepath."' target='page'>";?>Alpgen:setMLM</a></code>
flag for ALPGEN files. For LHEFs, instead, the setting of default 
parameters is controlled with the <code>JetMatching:setMad</code> flag:

<br/><br/><strong>JetMatching:setMad</strong>  <input type="radio" name="1" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="1" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
When enabled, the merging parameters are set according to the values 
in the LHEF header. Specifically, the header must set the
<code>ickkw</code>, <code>xqcut</code>, <code>maxjetflavor</code>
and <code>alpsfact</code> values, and <code>ickkw</code> must be nonzero. 
Note that these labels are Madgraph-specific. For other programs 
with LHEF output, or for Madgraph files lacking this information,
these parameters should be set by the user (or one can rely on the 
default values). The following parameters (described below) must then
be specified:
<ul>
<li> <code>JetMatching:doMerge = ickkw</code>, </li>
<li> <code>JetMatching:qCut = xqcut</code>, </li>
<li> <code>JetMatching:nQmatch = maxjetflavor</code>, </li>
<li> <code>JetMatching:clFact = alpsfact. </code> </li>
</ul>
With this flag on, the values from the LHEF for these parameters take 
precedence over other values.
  

<h2> Jet Matching parameters </h2>

A class <code>JetMatching</code>, derived from <code>UserHooks</code>, is 
used to define the basic structure of a parton-jet matching algorithm.
Two versions are implemented here, based on the FORTRAN code provided
by the ALPGEN and Madgraph packages, respectively:  
<code>JetMatchingAlpgen</code> and <code>JetMatchingMadgraph</code>.
The matching parameters are defined with the <code>JetMatching:*</code>
keyword.

<h3> Scheme and Usage </h3>

<br/><br/><strong>JetMatching:merge</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Master switch to activate parton-jet matching.
When off, all external events are accepted (unless they
are rejected due to weighting or event processing problems).
  

<br/><br/><table><tr><td><strong>JetMatching:scheme  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
The parton-jet MLM-style matching scheme.
<br/>
<input type="radio" name="3" value="1" checked="checked"><strong>1 </strong>:  The one inspired by the Madgraph matching code, here implemented in the <code>JetMatchingMadgraph</code> class.<br/>
<input type="radio" name="3" value="2"><strong>2 </strong>:  The one inspired by the ALPGEN matching code, here implemented in the <code>JetMatchingAlpgen</code> class.<br/>

<h3>Jet algorithm</h3>

The choice of jet algorithm and associated parameters can be adjusted with
the settings below. The PYTHIA8 internal <code>CellJet</code> and
<code>SlowJet</code> routines are used for jet finding.  See the
<?php $filepath = $_GET["filepath"];
echo "<a href='EventAnalysis.php?filepath=".$filepath."' target='page'>";?>Event Analysis</a> page for more details.

<br/><br/><table><tr><td><strong>JetMatching:jetAlgorithm  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
The choice of jet algorithm to use when merging against hard partons.
Currently, only <code>SlowJet</code> with the k<sub>T</sub> algorithm (and
<code>useStandardR = false</code>) is supported for Madgraph-style 
matching, while there is full freedom for the ALPGEN-style matching.
<br/>
<input type="radio" name="4" value="1" checked="checked"><strong>1 </strong>: The <code>CellJet</code> cone algorithm.<br/>
<input type="radio" name="4" value="2"><strong>2 </strong>: The <code>SlowJet</code> clustering algorithm.<br/>

<br/><br/><table><tr><td><strong>JetMatching:slowJetPower  </td><td>  &nbsp;&nbsp;(<code>default = <strong>-1</strong></code>; <code>minimum = -1</code>; <code>maximum = 1</code>)</td></tr></table>
The power to use in the <code>SlowJet</code> algorithm.
<br/>
<input type="radio" name="5" value="-1" checked="checked"><strong>-1 </strong>: The anti-k<sub>T</sub> algorithm.<br/>
<input type="radio" name="5" value="0"><strong>0 </strong>: The Cambridge/Aachen algorithm.<br/>
<input type="radio" name="5" value="1"><strong>1 </strong>: The k<sub>T</sub> algorithm.<br/>

<br/><br/><table><tr><td><strong>JetMatching:nEta  </td><td>  &nbsp;&nbsp;(<code>default = <strong>100</strong></code>; <code>minimum = 50</code>)</td></tr></table>
Specific to the <code>CellJet</code> algorithm, the number of bins 
in pseudorapidity.
</modepick>

<modepick name="JetMatching:nPhi" default="60" min="30">
Specific to the <code>CellJet</code> algorithm, the number of bins in phi.
</modepick>

<parm name="JetMatching:eTseed" default="1.5" min="0.0">
Specific to the <code>CellJet</code> algorithm, the minimum <ei>eT</ei> 
for a cell to be acceptable as the trial center of a jet.
</parm>

<parm name="JetMatching:eTthreshold" default="0.1">
Specific to the <code>CellJet</code> algorithm, cells with 
<ei>eT &lt; eTthreshold</ei> are completely neglected by the jet algorithm.
</parm>

<h3>Merging parameters</h3>

The following options are the three main parameters for the merging
procedure. Although here they are in principle free parameters, they should
be heavily influenced by the hard process generation cuts. 
These values can be set automatically based on the information in the
ALPGEN file or LHEF.

<parm name="JetMatching:eTjetMin" default="20.0" min="5.0">
For the <code>CellJet</code> algorithm, this gives the minimum transverse 
energy inside a cone for a jet to be accepted. For the <code>SlowJet</code> 
algorithm, this is instead the minimum transverse momentum required for a 
cluster to be accepted as a jet. For Madgraph-style matching, this 
parameter should match the <code>qCut</code> parameter described later.
</parm>

<parm name="JetMatching:coneRadius" default="0.7" min="0.1">
For the <code>CellJet</code> algorithm, this gives the size of the cone 
in <ei>(eta, phi)</ei> space drawn around the geometric center of the jet. 
For the <code>SlowJet</code> algorithm, this gives the <ei>R</ei> parameter.
</parm>

<parm name="JetMatching:etaJetMax" default="2.5" min="0.1">
For both jet algorithms, this defines the maximum pseudorapidity that
the detector is assumed to cover. In this context, however, it is tied
to the phase space region in which partons have been generated.
For the Alpgen-style matching, particles within 
<ei>etaJetMax + coneRadius</ei> are passed to the jet algorithm, 
with only jets within <ei>etaJetMax</ei> retained in the merging.
For the Madgraph-style matching, only particles within <ei>etaJetMax</ei> 
are used. 
</parm>

<h3>Exclusive mode</h3>

The following settings determine whether clustered jets which do not
match an original hard parton are allowed. They are typically permitted
in the highest jet multiplicity sample, where the parton shower may
produce extra hard jets, without risk of double counting. Any
extra jet produced by the shower must be softer than any matched light
jet, or else the event is vetoed.

<modepick name="JetMatching:exclusive" default="2" min="0" max="2">
Exclusive or inclusive merging.
<br/>
<input type="radio" name="6" value="0"><strong>0 </strong>:  The merging is run in inclusive mode. All partons must match jets, but  additional jets are allowed, provided they are not harder than the  matched jets. <br/>
<input type="radio" name="6" value="1"><strong>1 </strong>:  The merging is run in exclusive mode.  All partons must match jets,  and no additional jets are allowed. <br/>
<input type="radio" name="6" value="2"><strong>2 </strong>:  If <ei>nJet &lt; nJetMax</ei>, then the merging is run in exclusive  mode, otherwise it is run in inclusive mode. For Madgraph-style matching,  this is checked on an event-by-event basis, which is useful when an LHEF  contains a "soup" of partonic multiplicities.  If <ei>nJetMax &lt; 0</ei> or <ei>nJet &lt; 0</ei>, then the  algorithm defaults to exclusive mode. <br/>

<br/><br/><table><tr><td><strong>JetMatching:nJet  </td><td>  &nbsp;&nbsp;(<code>default = <strong>-1</strong></code>; <code>minimum = -1</code>)</td></tr></table>
When <code>JetMatching:exclusive = 2</code>, <ei>nJet</ei> indicates 
the minimum number of additional light jets in the incoming process.
This value may be set automatically.
</modepick>

<modepick name="JetMatching:nJetMax" default="-1" min="-1">
When <code>JetMatching:exclusive = 2</code>, <ei>nJetMax</ei> is used to 
indicate the maximum number of jets that will be matched. 
</modepick>

<h3>Jet matching</h3>

The following parameters control the criteria for matching a clustered jet
to a hard parton.

<modepick name="JetMatching:jetAllow" default="1" min="1" max="2">
Controls which particles are clustered by the jet algorithm.
<br/>
<input type="radio" name="7" value="1"><strong>1 </strong>:  This option explicitly disallows top quarks, leptons and photons. All other particle types are passed to the jet algorithm. <br/>
<input type="radio" name="7" value="2"><strong>2 </strong>:  No extra particles are disallowed. <br/>

<h3> Alpgen-specific parameters </h3>

<br/><br/><table><tr><td><strong>JetMatching:jetMatch  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
Criteria for matching a clustered jet to a parton.
<br/>
<input type="radio" name="8" value="1" checked="checked"><strong>1 </strong>:  This option can be used with both the <code>CellJet</code> and  <code>SlowJet</code> algorithms. The <ei>delta R</ei> between each parton  and jet is taken, and the minimal value compared against <ei>coneMatchLight * coneRadius</ei> for light jets or <ei>coneMatchHeavy * coneRadiusHeavy</ei> for heavy jets.  Note that by default <ei>coneRadiusHeavy = coneRadius</ei>, see below.  If below this value, the parton and jet are considered to match.  With <code>CellJet</code>, the <ei>delta R</ei> measure is in  <ei>(eta, phi)</ei>, while with <code>SlowJet</code> it is in  <ei>(y, phi)</ei>. <br/>
<input type="radio" name="8" value="2"><strong>2 </strong>:  This option can only be used with the <code>SlowJet</code> algorithm.  The hard partons are inserted into the parton level event as "ghost"  particles, but at the correct <ei>(y, phi)</ei> position. If this  particle is then clustered into a jet, it is considered a match. <br/>

<br/><br/><table><tr><td><strong>JetMatching:coneMatchLight </td><td></td><td> <input type="text" name="9" value="1.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.5</strong></code>; <code>minimum = 0.1</code>)</td></tr></table>
The <i>coneMatchLight</i> parameter used when
<code>JetMatching:jetMatch = 1</code>.
  

<br/><br/><table><tr><td><strong>JetMatching:coneRadiusHeavy </td><td></td><td> <input type="text" name="10" value="-1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.0</strong></code>)</td></tr></table>
The <i>coneRadiusHeavy</i> parameter used when
<code>JetMatching:jetMatch = 1</code>. When assigned a negative value,
the value of <code>JetMatching:coneRadius</code> is used.
  

<br/><br/><table><tr><td><strong>JetMatching:coneMatchHeavy </td><td></td><td> <input type="text" name="11" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.1</code>)</td></tr></table>
The <i>coneMatchHeavy</i> parameter used when
<code>JetMatching:jetMatch = 1</code>.
  

<h3>Madgraph-specific parameters </h3>

<br/><br/><table><tr><td><strong>JetMatching:qCut </td><td></td><td> <input type="text" name="12" value="10.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
k<sub>T</sub> scale for merging shower products into jets.
  

<br/><br/><table><tr><td><strong>JetMatching:nQmatch  </td><td>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 3</code>; <code>maximum = 6</code>)</td></tr></table>
Controls the treatment of heavy quarks.
<br/>
<input type="radio" name="13" value="5" checked="checked"><strong>5 </strong>:  All quarks (except top) are treated as light quarks for matching. <br/>
<input type="radio" name="13" value="4"><strong>4 </strong>:  Bottom quarks are treated separately.  Currently, they are unmatched. <br/>

<br/><br/><table><tr><td><strong>JetMatching:clFact </td><td></td><td> <input type="text" name="14" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>)</td></tr></table>
The <i>clFact</i> parameter determines how jet-to parton matching 
is done. A match is defined as a squared cluster scale that equals:
<br/><i>|clFact| * qCut</i> for inclusive mode,
<br/><i>|clFact| * max(qCut,min(p<sub>T</sub>(parton)))</i>
for exclusive mode, <i>clFact</i> &ge; 0, or
<br/><i>|clFact| * min(k<sub>T</sub>(parton))</i> for exclusive mode, 
<i>clFact</i> &lt; 0.
  

<h2>Alpgen-style parton-jet matching and merging</h2>

This section describes the Alpgen-style MLM merging algorithm for PYTHIA8. 
The most common reference to the algorithm is [<a href="Bibliography.php" target="page">Man02</a>]. In many 
respects, however, the implementation provided in the ALPGEN package should 
be considered the official description of the MLM merging procedure.
Although designed primarily to work with events generated with ALPGEN, it
can in principle also be used with events from a different source. This
should not be done without thought, however, and it is up to the user to
understand the details of the algorithm and the implications of using a
different hard process generator.

<p/>
First, either the <code>CellJet</code> or <code>SlowJet</code> jet 
algorithm is chosen. Both of these algorithms have an <i>R</i> and 
an <i>etaMax</i> parameter. In addition, <code>CellJet</code> has an 
<i>eTmin</i> and <code>SlowJet</code> has a <i>pTmin</i> parameter. 
These are the primary three parameters of the merging procedure, and in 
practice are set dependent on the cuts applied to the matrix element (ME) 
generation. We stress that the merging procedure is not tied to the 
geometry of a specific physical detector, but only to the match between 
the original partons and the resulting jets, using standard jet algorithms 
in the phase space region where partons have been generated.

<p/>
ME samples with different jet multiplicities are run through the event
generator, and the generation interrupted after parton showers have been
applied, but before resonance decays and beam remnants have been
processed. Note in particular that top quarks will not yet be decayed, 
which may lead to slight differences with the PYTHIA 6 interface included 
with the ALPGEN package. In what follows, the hardness measure of 
jets/partons is taken to be <i>eT</i> when <code>CellJet</code>
is used and <i>pT</i> when <code>SlowJet</code> is used. The hard 
system (ignoring all MPI systems) is then analysed:
<ul>
  <li>
    The particles in the original matrix element process are sorted into
    light partons, heavy partons and other particles. For backwards
    compatibility, a light parton is defined as the set <i>(d, u, s, c,
    b, g)</i> with zero mass. A heavy parton is defined as the set
    <i>(c, b, t)</i> with non-zero mass.
  </li>
  <li>
    All particles not originating from the heavy partons or other
    particles are passed to the jet algorithm and clustered.
  </li>
  <li>
    Clustered jets are matched to the light partons in the original ME
    process. There are two different methods which can be used:
    <ul>
      <li>
        Method 1: The following is done for each parton, in order
        of decreasing hardness. The <i>delta R</i> between the parton
        and all jets is calculated and the smallest value taken. If
        this is less than the jet <i>R</i> parameter, possibly
        multiplied by a constant, the jet and parton are considered to
        match, and the jet is removed from further consideration.
        Note that for <code>CellJet</code> the <i>delta R</i> measure 
        is in <i>(eta, phi)</i>, while for <code>SlowJet</code>, it is
        in <i>(y, phi)</i>.
      </li>
      <li>
        Method 2: This method is only possible when using the 
        <code>SlowJet</code> algorithm. Before the clustering is performed, 
        extremely soft "ghost" particles are added to the event at the
        <i>(y, phi)</i> coordinates of the original matrix element
        partons. If such a particle is clustered into a jet, the parton
        and jet are considered to match. The idea of "ghost" particles
        was originally introduced by FastJet as a way to measure jet
        areas [<a href="Bibliography.php" target="page">Cac06</a>] and should not affect clustering with an
        infrared-safe jet algorithm.
      </li>
    </ul>
  <li>
    If there is a light ME parton remaining which has not been matched
    to a jet, then the event is vetoed. If all ME partons have been
    matched to a jet, but there are still some extra jets remaining,
    then two options are possible:
    <ul>
      <li>
        Exclusive mode: the event is vetoed. This is typically used when
        there are ME samples with higher jet multiplicities, which would
        fill in the extra jets.
      </li>
      <li>
        Inclusive mode: the event is retained if the extra jets are softer
        than the softest matched jet. This is typically used when
        there is no ME sample with higher jet multiplicity, so the parton
        shower should be allowed to give extra jets.
    </ul>
  </li>
  <li>
    All particles originating from the heavy partons are passed to the
    jet algorithm and clustered.
  </li>
  <li>
    The clustered jets are again matched to the original partons, but
    there is no requirement for a match to be present; all matched jets
    are immediately discarded. The matching procedure is much the same
    as for light partons, but with two differences when <i>delta R</i>
    matching is used. First, a different <i>R</i> parameter than that
    used by the jet algorithm may optionally be given. Second, all jets
    that are within the given radius of the parton are matched, not
    just the one with the smallest <i>delta R</i> measure. If there
    are still extra jets remaining then in exclusive mode the event is
    immediately vetoed, while in inclusive mode the event is retained if
    the extra jets are softer than the softest <em>light</em> matched jet.
  </li>
</ul>

<p/>
Some different options are provided, specified further above in the 
parameters section. These are set so that, by default, the algorithm 
closely follows the official MLM interface provided in the ALPGEN package. 

<p/>
All vetoing of events is done through the usual
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>User Hooks</a> machinery, and is
therefore already taken into account in the cross section. In the output 
from <code><?php $filepath = $_GET["filepath"];
echo "<a href='EventStatistics.php?filepath=".$filepath."' target='page'>";?>Pythia::stat()</a></code>,
the difference between the "Selected" and "Accepted" columns gives the
number of events that have not survived the vetoing procedure. It is
still the responsibility of the user to add together the results from
runs with different jet multiplicities. In the simplest case, when
ALPGEN input is used and the hard process parameters are used to guide
the merging procedure, it is enough to set the 
<code>JetMatching:nJetMax</code> parameter.

<h2>Madgraph-style parton-jet Merging and Matching</h2>

<p/>
This section describes the Madgraph-style parton-jet matching algorithm 
for PYTHIA8.

<p/>
First, the k<sub>T</sub> jet algorithm is applied using the PYTHIA8 
<code>SlowJet</code> implementation. The <code>useStandardR = false</code>
is used, ie. the <i>(delta R)^2</i> separation is defined as  
<i>2 (cosh(delta y) - cos(delta phi))</i> rather than the more common
<i>(delta y)^2 + delta phi)^2</i>. The <i>R</i>, <i>etaMax</i>, 
and a <i>pTmin</i> parameters are specified. By default, <i>R = 1</i> 
and <i>pTmin = qCut </i>. It is not recommended to change these.
These should match the algorithm parameters used in the Madgraph 
Matrix Element (ME) generation.

<p/>
ME samples with different jet multiplicities are run through the event 
generator, and the generation is interrupted after parton showers have 
been applied, but before resonance decays and beam remnants have been 
processed. In what follows, the hardness measure of jets/partons is taken 
to be <i>k<sub>T</sub></i> relative to <i>qCut</i>.
The hard system (ignoring all MPI systems) is analyzed:
<ul>
  <li>
    The hard partons in the original matrix element process, provided by
    the LHEF, are sorted into light partons, heavy partons and other 
    particles. A heavy parton is defined by the 
    <code>JetMatching:nQmatch</code> or by the <code>maxjetflavor</code>
    value in the LHEF. <i>nQmatch</i> refers to the absolute value of
    the quark PDG identity code.
  </li>
  <li>
    All partons arising from the parton shower are sorted based on their 
    motherhood. A showered parton arising from a heavy parton or "other" 
    parton classified in the previous step is not passed to the jet 
    algorithm. All other partons are clustered into light jets.
  </li>
  <li> It is checked whether there are "too few" or "too many" light jets.
    If the number of light jets is less than the number of light partons 
    defined by <i>nQmatch</i>, the event is vetoed. If the number is 
    larger, the event is vetoed only in exclusive mode (defined below).
  </li>
  <li> In exclusive mode, the number of jets matches the number of light 
    partons. In inclusive mode, the jets are re-clustered until the number 
    of jets equals the number of light partons. Next, each light hard 
    parton is clustered, one at a time, with the jets until a match is found. 
    A match is defined as a squared cluster scale that equals:
    <ul>
      <li><i>|clFact| * qCut</i> for inclusive mode,</li>
      <li><i>|clFact| * max(qCut,min(p<sub>T</sub>(parton)))</i>  
         for exclusive mode, <i>clFact</i> &ge; 0, or</li>
      <li><i>|clFact| * min(k<sub>T</sub>(parton))</i> for exclusive 
         mode, <i>clFact</i> &lt; 0.</li>
    </ul>
    If no match is found, the event is vetoed. When a parton
    matches a jet, the jet is removed from the collection, and
    the process continues. The process terminates when all partons
    are matched to a jet, or a parton is unmatched.  
  </li>
  <li>
    All particles originating from the heavy partons are not used.
  </li>
</ul>
In exclusive mode, it is expected that ME samples with higher parton 
multiplicity are available to fill the phase space above <i>qCut</i>.
The inclusive mode is when there are no such samples, and the parton 
shower is used to fill the phase space.

<p/>
Some different options are provided, specified further above. These
are set so that, by default, the algorithm closely follows the
FORTRAN interface <code>ME2Pythia</code> provided in the Madgraph 
package. 

<p/>
All vetoing of events is done through the usual
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>User Hooks</a> machinery, and is
therefore already taken into account in the cross section. In the output 
from <code><?php $filepath = $_GET["filepath"];
echo "<a href='EventStatistics.php?filepath=".$filepath."' target='page'>";?>Pythia::stat()</a></code>,
the difference between the "Selected" and "Accepted" columns gives the
number of events that have not survived the vetoing procedure. It is
still the responsibility of the user to add together the results from
runs with different jet multiplicities. In the simplest case, when
the hard process parameters are used to guide the merging procedure, 
events will be matched in the exclusive mode.

<p/>
 
<h3>A note on combining UserHooks</h3>

As have been noted above, the matching is implemented using classes
derived from the <code><?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>UserHooks</a></code>
class, thereby gaining access to the event generation process at the
relevant locations. For native ALPGEN files, which do not adhere to
the Les Houches standards, it is also necessary to intervene with
a <code>UserHooks</code>-derived  <code>AlpgenHooks</code> to handle 
the extraction and setting of relevant extra information.   

<p/>
One must then combine multiple <code>UserHooks</code> classes, 
such that the functionality of both is present. A prerequisite 
is that the different <code>UserHooks</code> classes should be 
declared with virtual inheritance, e.g.
<pre>
  class JetMatching : virtual public UserHooks
</pre>
Without this option, when combining two <code>UserHooks</code>-derived 
classes, two copies of the base <code>UserHooks</code> class would be 
created, leading to ambiguities.

<p/>
The two first classes in <code>CombineMatchingInput.h</code> combine
ALPGEN input with the two different matching schemes, e.g. for the first
<pre>
class JetMatchingAlpgenInputAlpgen : public AlpgenHooks, 
  public JetMatchingAlpgen {
public:
  // Constructor and destructor.
  JetMatchingAlpgenInputAlpgen(Pythia& pythia) : AlpgenHooks(pythia), 
    JetMatchingAlpgen() { }
  ~JetMatchingAlpgenInputAlpgen() {}
  // Initialisation.
  virtual bool initAfterBeams() {
    if (!AlpgenHooks::initAfterBeams()) return false;
    if (!JetMatchingAlpgen::initAfterBeams()) return false;
    return true;
  }
  // Process level vetos.
  virtual bool canVetoProcessLevel() { 
    return JetMatchingAlpgen::canVetoProcessLevel();    
  }
  ....
};
</pre>
This class inherits from both <code>AlpgenHooks</code> and 
<code>JetMatchingAlpgen</code>. Any functions which are present 
in both classes should be overridden with a function that calls 
the different parent methods in the desired order. In the
above example, the only shared methods are the constructor and
<code>initAfterBeams()</code>.

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
$data = "JetMatching:setMad = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "JetMatching:merge = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "1")
{
$data = "JetMatching:scheme = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "1")
{
$data = "JetMatching:jetAlgorithm = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "-1")
{
$data = "JetMatching:slowJetPower = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "100")
{
$data = "JetMatching:nEta = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "-1")
{
$data = "JetMatching:nJet = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "1")
{
$data = "JetMatching:jetMatch = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "1.5")
{
$data = "JetMatching:coneMatchLight = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "-1.0")
{
$data = "JetMatching:coneRadiusHeavy = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "1.0")
{
$data = "JetMatching:coneMatchHeavy = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "10.0")
{
$data = "JetMatching:qCut = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "5")
{
$data = "JetMatching:nQmatch = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "1.0")
{
$data = "JetMatching:clFact = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
