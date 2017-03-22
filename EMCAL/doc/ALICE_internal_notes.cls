%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Name:     cernphprep.cls                                       %
%                                                                 %
%  Purpose:  Set LaTeX parameters for CERN PH preprint            %
%                                                                 %
%  Rewrite of CERN preprint class for LHC results (April 2010)    %
%                                                                 %
% 12 April 2010 (version 1.00)                                    %
%   - Special treatment of title page                             %
%   - Basic point size 11pt                                       %
%   - Based on standard LaTeX article.cls                         %
% 16 April 2010 (version 1.01)                                    %
%   - Add \EXPdate command                                        %
%   - Add CERN option to print CERN logo at left                  %
% 18 April 2010 (version 1.02)                                    %
%   - Extend \@alph enumeration to gz (above 200)                 %
% 19 April 2010 (version 1.03)                                    %
%   - Introduce \LHCEXP as flag to tell we are handling LHC exp.  %
%   - Introduce manyauthors flag for handling many authers        %
%     This sets the flag MANYAUTHORS which (re)defines a few      %
%     commands that make it easy to handling long lists of        %
%     authors and the cross-references to institutes, etc.        %
%  6 June  2010 (version 1.04)                                    %
%   - Introduce \CERNCopyright command to define copyright text   %
% 30 June  2010 (version 1.05)                                    %
%   - Introduce \COLLABORATION command to define coll. name.      %
%     Initialized to ``unnamed'' and used in \CERNCopyright       %
%   - Define new text for  \CERNCopyright command as proposed by  %
%     Jens Vigen (after SIPB discussion).                         %
%  8 July  2010 (version 1.06)                                    %
%   - Add COMPASS option                                          %
%                                                                 %
%  Last mod.:  8 July 2010 MG (Michel Goossens, IT/UDS)           %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\fileversion{1.06}
\def\filedate{8 July 2010}
\def\docdate{8 July 2010}
\typeout{Document Style 'cernphprep'. Version <\fileversion>, <\filedate>}
% -------------------------------- identification ---------------------------
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cernphprep}[2010/07/08]
% ---------------------------  declaration of options --
% Define experiment logo for first page
\providecommand*\EXPLOGO{nonelogo}
\providecommand*\LHCEXP{no}
\providecommand*\COLLABORATION{unnamed}
\providecommand*\MANYAUTHORS{no}
\providecommand*\CLEARPAGE{no}
\DeclareOption{nocleardouble}{\renewcommand\CLEARPAGE{yes}}
\DeclareOption{manyauthors}{\renewcommand\MANYAUTHORS{yes}}
\DeclareOption{ALICE}{\renewcommand{\EXPLOGO}{alicelogo}%
                      \renewcommand\LHCEXP{yes}\renewcommand\COLLABORATION{ALICE}}
\DeclareOption{ATLAS}{\renewcommand{\EXPLOGO}{atlaslogo}%
                      \renewcommand\LHCEXP{yes}\renewcommand\COLLABORATION{ATLAS}}
\DeclareOption{CERN}{\renewcommand{\EXPLOGO}{cernlogo}}
\DeclareOption{CMS}{\renewcommand{\EXPLOGO}{cmslogo}%
                    \renewcommand\LHCEXP{yes}\renewcommand\COLLABORATION{CMS}}
\DeclareOption{COMPASS}{\renewcommand{\EXPLOGO}{compasslogo}%
                    \renewcommand\LHCEXP{no}\renewcommand\COLLABORATION{COMPASS}}
\DeclareOption{LHCb}{\renewcommand{\EXPLOGO}{lhcblogo}%
                     \renewcommand\LHCEXP{yes}\renewcommand\COLLABORATION{LHCb}}
\DeclareOption{LHCB}{\renewcommand{\EXPLOGO}{lhcblogo}%
                     \renewcommand\LHCEXP{yes}\renewcommand\COLLABORATION{LHCb}}
\DeclareOption{LHCf}{%\renewcommand{\EXPLOGO}{lhcflogo}%
                     \renewcommand\LHCEXP{yes}\renewcommand\COLLABORATION{LHCf}}
\DeclareOption{LHCF}{%\renewcommand{\EXPLOGO}{lhcflogo}%
                     \renewcommand\LHCEXP{yes}\renewcommand\COLLABORATION{LHCf}}
\DeclareOption{TOTEM}{%\renewcommand{\EXPLOGO}{totemlogo}%
                      \renewcommand\LHCEXP{yes}\renewcommand\COLLABORATION{TOTEM}}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass[11pt,a4paper,twoside]{article}
\renewcommand{\rmdefault}{ptm}
\usepackage{mathptmx}

\RequirePackage{graphicx,amsmath,amssymb,array,tabularx,url,xspace,subfigure}
\RequirePackage{ifthen}
 
\typeout{Document Class "cernphprep" 1.06 (8 July 2010)}

%% Defince collaboration is it cannot be inferred from option 

\providecommand\DEFCOL[1]{\renewcommand\COLLABORATION{#1}}

%%% Customization of visible document elements

\renewcommand\thesection      {\@arabic\c@section}
\renewcommand\thesubsection   {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection.\@arabic\c@subsubsection}
%\renewcommand\theparagraph{\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\baselinestretch{1.}
\renewcommand{\contentsname}{Contents}
\renewcommand{\listfigurename}{List of Figures}
\renewcommand{\listtablename}{List of Tables}
\renewcommand{\refname}{References}
\renewcommand{\indexname}{Index}
\renewcommand{\figurename}{Fig.}
\renewcommand{\tablename}{Table}
\renewcommand{\partname}{Part}
\renewcommand{\appendixname}{Appendix}
\renewcommand{\abstractname}{Abstract}

% Short and long pages

\newcommand{\longpage}[1][1]{\enlargethispage{#1\baselineskip}}
\newcommand{\shortpage}[1][1]{\enlargethispage{-#1\baselineskip}}

%  ****************************************
%  *            Lists                     *
%  ****************************************
%
% Zero the vertical spacing parameters
\def\Nospacing{\itemsep=2pt\topsep=0pt\partopsep=0pt\parskip=0pt\parsep=0pt}
% Redefine the environments in terms of the original values
\newenvironment{Itemize}{\vspace*{-1mm}\itemize\Nospacing}%
                               {\endlist\vspace*{-1mm}}
\newenvironment{Enumerate}{\vspace*{-1mm}\enumerate\Nospacing}
                               {\endlist\vspace*{-1mm}}
\newenvironment{Description}{\vspace*{-1mm}\description\Nospacing}
                               {\endlist\vspace*{-1mm}}
\renewcommand\labelitemi{\normalfont\textendash}
\renewcommand\labelitemii{\normalfont\textendash}
\renewcommand\labelitemiii{\textasteriskcentered}
\renewcommand\labelitemiv{\textperiodcentered}

%  ****************************************
%  *            PAGE LAYOUT               *
%  ****************************************
%
% All margin dimensions measured from a point one inch from top and side
% of page.
 
% SIDE MARGINS:

\textheight 24cm       % Height of text (including footnotes and figures,
                       % excluding running head and foot).
\textwidth 16cm        % Width of text line.
\topmargin -28pt       % offset of top header from voffset (1 inch)
\headheight 13.6pt

% SIDE MARGINS:  %%%% A4 : pagewidth=160mm 
\if@twoside   % Values for two-sided printing:
      %%% margins : inner     =  30 mm =  86 pt
      %%%           outer     =  20 mm =  56 pt
      %%% \odd/evensidemargin are with respect to hoffset -->
   \oddsidemargin   14pt    %   Left margin on odd-numbered pages.
   \evensidemargin -14pt    %   Left margin on even-numbered pages.
\else         % Values for one-sided printing:
   \oddsidemargin    0pt    %   Note that \oddsidemargin = \evensidemargin
   \evensidemargin   0pt
\fi

%%% Captions

\setlength\abovecaptionskip{4pt plus 2pt minus 2pt}   % default is for caption at bottom.
                                                      % use topcapt.sty: \topcaption for cap at top
\setlength\belowcaptionskip{4pt plus 2pt minus 2pt}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{{\small\textbf{#1:} #2}}%
  \ifdim \wd\@tempboxa >\hsize
    {\small\textbf{#1:} #2}\par
  \else
    \global \@minipagefalse
    \hbox to\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}

%%% make footnote robust (if it appears in argument of, e.g., \section)
\let\oldfootnote\footnote
\def\footnote{\protect\oldfootnote}
\RequirePackage{ifthen}

%%% Set up headings for section commands %%%%%%%%%%%%%%%%%%%%%%

\renewcommand\section{\@startsection
{section}%
{1}%
{\z@}%
{-2.5ex \@plus -1ex \@minus -.2ex}%       Was: {-3.5ex \@plus -1ex \@minus -.2ex}%
{0.3ex \@plus.2ex}%                       Was: {2.3ex \@plus.2ex}
{\large\bfseries\mathversion{bold}}}%              Was: \normalfont ...

\renewcommand\subsection{\@startsection
{subsection}%
{2}%
{\z@}%
{-1.0ex\@plus -1ex \@minus -.2ex}%       Was: {-3.25ex\@plus -1ex \@minus -.2ex}%
{0.05ex \@plus .2ex}%                     Was:   {1.5ex \@plus .2ex}%
{\normalsize\bfseries\mathversion{bold}}}%           Was: \normalfont ...

\renewcommand\subsubsection{\@startsection
{subsubsection}%
{3}%
{\z@}%
{-1.0ex\@plus -1ex \@minus -.2ex}%       Was: {-3.25ex\@plus -1ex \@minus -.2ex}
{0.05ex \@plus .2ex}%                     Was:   {1.5ex \@plus .2ex}%
{\normalsize\itshape\bfseries\mathversion{bold}}}%           Was: \normalfont ...

\renewcommand\paragraph{\@startsection
{paragraph}%
{4}%
{\z@}%
{-1.0ex\@plus -1ex \@minus -.2ex}%       Was: {3.25ex \@plus1ex \@minus.2ex}%             Was:
{-1em}%                             Was:  If negative, means in-line. Space is then the indentation.
{\normalsize\bfseries\mathversion{bold}}}%           Was: \normalfont ...

\renewcommand\subparagraph{\@startsection
{subparagraph}%
{5}%
{\parindent}%
{-1.0ex\@plus -1ex \@minus -.2ex}%
{-1em}%
{\normalsize\itshape\bfseries\mathversion{bold}}}%           Was: \normalfont ...

%
%  ****************************************
%  *       Appendix bookkeeping           *
%  ****************************************
%
 
\renewcommand\appendix{\par
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \setcounter{figure}{0}%
  \setcounter{table}{0}%
  \setcounter{equation}{0}%
  \gdef\thesection{\@Alph\c@section}%
  \gdef\thesubsection{\@Alph\c@section.\@arabic\c@subsection}
  \gdef\theequation{\@Alph\c@section.\@arabic\c@equation}
  \@addtoreset {equation}{section}%
  \gdef\thetable{\@Alph\c@section.\@arabic\c@table}%
  \@addtoreset {table}{section}%
  \gdef\thefigure{\@Alph\c@section.\@arabic\c@figure}%
  \@addtoreset {figure}{section}%
}

% 
%  ************************************
%  *   thebibliography environement   *
%  ************************************
%

\setlength\bibindent{1.5em}
\renewenvironment{thebibliography}[1]
     {\section*{\refname}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \Nospacing
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\renewcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% \begin{titlepage}
%% \PHnumber{2010--xxx}
%% \PHdate{DD Month 2010}
%% \EXPnumber{IT-UDS-HUS-001}
%% \EXPdate{23 June 2010}
%% \DEFCOL{CDS-Library}
%% \title{Preparing a CERN PH preprint}
%% \author{Michel Goossens/IT-UDS\thanks{\texttt{Email:michel.goossens@cern.ch}}}
%%   or
%% \author{LHC Collaboration\thanks{see list of authors in Appendix}}
%% \begin{abstract}
%%   ...
%% \end{abstract}
%% \end{titlepage}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% title page %%%%%%%%%%%%%%%%%%%%%%%%%5
\def\maketitle{% Dummy maketitle
      \typeout{Dummy maketitle --- Use titlepage instead}}
% Definition of environment "titlepage"
\newcommand\EXPdate[1]{\def\@EXPdate{#1}}\EXPdate{\empty}
\newcommand\EXPnumber[1]{\def\@EXPnumber{#1}}\EXPnumber{\empty}
\newcommand\PHnumber[1]{\def\@PHnumber{#1}}\PHnumber{To be specified}
\newcommand\PHdate[1]{\def\@PHdate{#1}}\PHdate{\today}
%\providecommand\CERNCopyright{\empty}%CERN, Geneva, 2010. Copyright and licence text}}
\providecommand\CERNCopyright{\parbox{\textwidth}{%
   Copyright CERN, for the benefit of the \COLLABORATION{} Collaboration.\\
   This article is distributed under the terms of Creative Commence  
   Attribution License (CC-BY-3.0), which permits any use provided the
   original author(s) and source are credited.}}
\def\titlepage{\begingroup\clearpage
  \renewcommand{\thefootnote}{\fnsymbol{footnote}}
  \setcounter{footnote}{0}\setcounter{page}{1}\pagestyle{empty}}
%
% CERN title --> 
% prints CERN header 
% then logos, preprint number(s), date 
% at bottom of page: copyright notice 
%
\renewcommand\title[1]{\setlength{\unitlength}{1mm}%
\begin{picture}(0.001,0.001)%
%\graphpaper(0,-280)(210,290)
%\put(80,-15){\makebox(0,0){\LARGE\itshape Final draft version \today}}
\put(80,10){\makebox(0,0){\LARGE EUROPEAN ORGANIZATION FOR NUCLEAR RESEARCH}}
%----
% Experiment's logo (optional) 
\put(0,-20){\makebox(0,0)[lb]{\includegraphics[height=2cm]{\EXPLOGO}}}
% CERN logos
\put(160,-20){\makebox(0,0)[rb]{\includegraphics[height=2cm]{cernlogo}}}
% Experiment's preprint number (optional) 
\put(0,-25){\makebox(0,0)[lb]{\normalsize\@EXPnumber}}
% Experiment date of approval
\put(0,-30){\makebox(0,0)[lb]{\normalsize\@EXPdate}}
% CERN PH Preprint number 
\put(160,-25){\makebox(0,0)[rb]{\normalsize ALICE-INT-2012-xxx}}
% CERN PH Preprint date of approval
\put(160,-30){\makebox(0,0)[rb]{\normalsize\@PHdate}}
%\put(80,-250){\makebox(0,0){\CERNCopyright}} %%% Copyright notice
\end{picture} 

\vspace*{3cm}               
\begin{center}
\mbox{}\\[25mm]
\Large\bf#1\\[10mm]
\end{center}
\gdef\Title{#1}
}%%  end of title
% Submitted{}
\def\Submitted#1{\begin{center}\it#1\end{center}}
% conference{}
\def\Conference#1{\begin{center}\it#1\end{center}}
% note{}
\def\Note#1{\begin{center}\it#1\end{center}}
% dedication{}
\def\Dedication#1{\begin{center}\it#1\end{center}}

% Author or collaboration

\def\Collaboration#1{%% Treat the collaboration
 \begin{center}#1\end{center}\@thanks
 \gdef\Authors{\begingroup\let\thanks\@gobble\let\@thanks\@empty#1\endgroup}}

\def\author#1{%% Treat the list of authors
% \setcounter{footnote}{0}\def\@currentlabel{}%
% \begingroup\def\thefootnote{\arabic{footnote}}
% \def\@makefnmark{\hbox{$^{\rm\@thefnmark)}$}}
 \global\@topnum\z@ \begin{center}{\lineskip.5em
 \begin{tabular}[t]{c}#1\end{tabular}\par}
 \end{center}\par\vskip1.5em\@thanks%\endgroup
 \gdef\Authors{\begingroup\let\thanks\@gobble\let\@thanks\@empty#1\endgroup}}

\newcommand\ShortAuthor[1]{\gdef\Authors{#1}}\gdef\Authors{\empty}
\newcommand\ShortTitle[1]{\gdef\Title{#1}}\gdef\Title{\empty}

\newenvironment{Authlist}{\center}{\endcenter}

\renewenvironment{quotation}
               {\list{}{\listparindent 0pt%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        1em\@plus\p@}%
                \item\relax}
               {\endlist}

%%%%%%%%%%%% Institute definitions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\@ifundefined{c@instfoot}{\newcounter{instfoot}}{}
\def\Iref#1{$^{\Inst{#1})}$}   \let\Instref=\Iref
\def\Inst#1{\@ifundefined{r@#1}{{\bf ??}}%
{\edef\@tempa{\@nameuse{r@#1}}\expandafter\@car\@tempa \@nil\null}}
%
%%AS If we found `Instfoot' inside `Authlist' env, we change the
%    formatting.
%
\def\Instfoot#1#2{\def\@temp{Authlist}%
\ifx\@temp\@currenvir
    \\ \emph{#2} \\[1em]%
\else
    \refstepcounter{instfoot}\label{#1}%
    \xdef\@thefnmark{\arabic{instfoot}}\@footnotetext{#2}%
\fi
}

%
%%%%%%%%%%%% Anotation definitions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%\@ifundefined{c@anotfoot}{\newcounter{anotfoot}}{}
%\@ifundefined{mathrm}{\newcommand{\mathrm}[1]{{\rm#1}}}{}
\def\Aref#1{$^{\mathrm{#1})}$} \let\Anotref=\Aref
\def\Anotfoot#1#2{\xdef\@thefnmark{#1}\@footnotetext{#2}}
%
%%%%%%%%%%%% Anotation references   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\AAref#1#2{$^{\mathrm{#1,#2)}}$}       \let\AnotAnotref=\AAref
\def\IAref#1#2{$^{\Inst{#1},\mathrm{#2)}}$}\let\InstAnotref=\IAref
\def\IIref#1#2{$^{\Inst{#1},\Inst{#2})}$}  \let\InstInstref=\IIref
\def\IIIref#1#2#3{$^{\Inst{#1},\Inst{#2},\Inst{#3})}$}
\def\IIAref#1#2#3{$^{\Inst{#1},\Inst{#2},\mathrm{#3)}}$}
%%%%%%%%%%%% Annospec environment   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newenvironment{Annospec}[1]{% Give width for offset of text
\vskip2mm\bgroup
%\def\@makefnmark{\hbox{$^{\mathrm{\@thefnmark)}}$}}
\long\def\@makefntext##1{\parindent=\z@
\begin{list}{}%
   {\labelwidth=#1\leftmargin\labelwidth% set global width
    \labelsep3pt%  horizontal separation between term and item
    \itemsep0pt%   vertical separation between two items
    \parsep0pt%    vertical separation two paragraphs in an item
    \topsep0pt%    vertical separation text/list
    \footnotesize}
    \item [\hfill $^{\mathrm{\@thefnmark)}}$] ##1\end{list}}}% ***** end of makefntext
{\egroup}% ***** end of environment{Annospec} *****
 
\def\Deceased{\textsuperscript{,*}}
\ifthenelse{\equal{\MANYAUTHORS}{yes}}{% Commands for long author lists
  \gdef\Adef#1{\label{LA#1}}
  \gdef\Aref#1{\textsuperscript{,\ref{LA#1}}}
  \gdef\Arefs#1{\textsuperscript{\ref{LA#1}}} %% \Aref short (without comma)
  \gdef\Idef#1{\label{LI#1}\gdef\INTlatex{#1}}
  \gdef\Irefn#1{\textsuperscript{\ref{LI#1}}}
  \gdef\Irefnn#1#2{\textsuperscript{\ref{LI#1},\ref{LI#2}}}
  \gdef\Idefs#1{\label{LI\INTlatex#1}\textsuperscript{(#1)}}
  \gdef\Irefs#1#2{\textsuperscript{\ref{LI#1#2}#2}}
  \gdef\Irefss#1#2#3#4{\textsuperscript{\ref{LI#1#2}#2,\ref{LI#3#4}#4}}
  \gdef\And{, }
  \gdef\NNospacing{\labelsep=0pt\itemsep=0pt\topsep=0pt\partopsep=0pt\parskip=0pt\parsep=0pt}
  \renewenvironment{Authlist}{\vspace*{-1.5ex}\enumerate\NNospacing}
                              {\endlist\vspace*{-1.5ex}}
}{}
 
% End of titlepage
\def\endtitlepage{%% Reset counters
\setcounter{footnote}{0}\let\titlepage\relax\vfill
\ifthenelse{\equal{\CLEARPAGE}{yes}}{% Commands for long author lists
\clearpage\setcounter{page}{2}}{\cleardoublepage\setcounter{page}{3}}%
\pagestyle{plain}\pagenumbering{arabic}%
\gdef\@thanks{}\gdef\@author{}\gdef\@title{}\let\thanks\relax\endgroup}

% Local CERN stuff
 
\newcommand{\bra}[1]{\left\langle #1\right|} % <bra|
\newcommand{\ket}[1]{\left| #1\right\rangle} % |ket>
\newcommand{\braket}[2]{\left\langle#1 | #2\right\rangle} % <bra|ket>
 
% Figure placement

\setcounter{topnumber}{5}
\setcounter{totalnumber}{6}
\renewcommand{\topfraction}{0.99}
\renewcommand{\dbltopfraction}{0.99}
\renewcommand{\bottomfraction}{0.99}
\renewcommand{\textfraction}{0.15}

\setlength{\intextsep}{10pt plus 2pt minus 4pt}  % Space between in-text figures and text.
\setlength{\floatsep}{10pt plus 2pt minus 4pt}   % Space between adjacent floats moved to top or bottom of text page.
\setlength{\textfloatsep}{10pt plus 2pt minus 4pt}% Space between main text and floats at top or bottom of page.
\setlength{\dblfloatsep}{10pt plus 2pt minus 4pt} % Same for double column
\setlength{\dbltextfloatsep}{10pt plus 2pt minus 4pt}

\parskip.5\baselineskip
\parindent0mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Equations                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\allowdisplaybreaks

% page headers

%\usepackage[fit]{truncate}
%\usepackage{fancyhdr}
%\pagestyle{fancy}
%\fancyhf{}
%\fancyhead[LE,RO]{\thepage}
%\fancyhead[RE]{\small\leftmark}
%\fancyhead[LO]{\small\rightmark}
%\fancyhead[LO]{\truncate{.85\headwidth}{\Title}}
%\fancyhead[RE]{\truncate{.85\headwidth}{\Authors}}


\def\@alph#1{%
  \ifcase#1\or a\or b\or c\or d\or e\or f\or g\or h\or i\or j\or k\or
  l\or m\or n\or o\or p\or q\or r\or s\or t\or u\or v\or w\or x\or
  y\or z%
  \or aa\or ab\or ac\or ad\or ae\or af\or ag\or ah\or ai\or aj\or
  ak\or al\or am\or an\or ao\or ap\or aq\or ar\or as\or at\or au\or
  av\or aw\or ax\or ay\or az%
  \or ba\or bb\or bc\or bd\or be\or bf\or bg\or bh\or bi\or bj\or
  bk\or bl\or bm\or bn\or bo\or bp\or bq\or br\or bs\or bt\or bu\or
  bv\or bw\or bx\or by\or bz%
  \or ca\or cb\or cc\or cd\or ce\or cf\or cg\or ch\or ci\or cj\or
  ck\or cl\or cm\or cn\or co\or cp\or cq\or cr\or cs\or ct\or cu\or
  cv\or cw\or cx\or cy\or cz%
  \or da\or db\or dc\or dd\or de\or df\or dg\or dh\or di\or dj\or
  dk\or dl\or dm\or dn\or do\or dp\or dq\or dr\or ds\or dt\or du\or
  dv\or dw\or dx\or dy\or dz%
  \or ea\or eb\or ec\or ed\or ee\or ef\or eg\or eh\or ei\or ej\or
  ek\or el\or em\or en\or eo\or ep\or eq\or er\or es\or et\or eu\or
  ev\or ew\or ex\or ey\or ez%
  \or fa\or fb\or fc\or fd\or fe\or ff\or fg\or fh\or fi\or fj\or
  fk\or fl\or fm\or fn\or fo\or fp\or fq\or fr\or fs\or ft\or fu\or
  fv\or fw\or fx\or fy\or fz%
  \or ga\or gb\or gc\or gd\or ge\or gf\or gg\or gh\or gi\or gj\or
  gk\or gl\or gm\or gn\or go\or gp\or gq\or gr\or gs\or gt\or gu\or
  gv\or gw\or gx\or gy\or gz%
  \else\@ctrerr\fi}



\endinput
 
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
