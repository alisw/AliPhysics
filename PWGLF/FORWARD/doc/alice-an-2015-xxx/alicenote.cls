%%
%% This is file `alicenote.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% alice_template.dtx  (with options: `cls')
%% 
%% Copyright (c) 2013 ALICE Collaboration.
%% 
%% See the file alice_template.dtx for further comments and documentation.
%% 

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{alicenote}
\newif\if@alice@public@note\@alice@public@notefalse
\newif\if@alice@many@uthors\@alice@many@uthorsfalse
\newif\if@alice@compat\@alice@compatfalse
\DeclareOption{public}{\@alice@public@notetrue}
\DeclareOption{analysis}{\@alice@public@notefalse}
\DeclareOption{manyauthors}{\@alice@many@uthorstrue}
\DeclareOption{compat}{\@alice@compattrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass{article}
\RequirePackage{graphicx}
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{amstext}
\RequirePackage{alicemisc}
\RequirePackage[%
  hyperindex=true,%     Make links on pages
  pageanchor=true,%     Put anchor on every page
  hyperfootnotes=true,% Link to footnotes
  colorlinks=true,%     Colors links
  linkcolor=aliceyellow!25!black,
  anchorcolor=alicered!85!black,
  citecolor=alicepurple!85!black,
  urlcolor=alicered,
  bookmarks=true,%      Make PDF index
  bookmarksopenlevel=2,% By default this deep
  bookmarksnumbered=true % Number bookmarks
]
{hyperref}
\RequirePackage{mathptmx}
\RequirePackage{truncate}
\RequirePackage[a4paper,%
                textwidth=16cm,%
                textheight=24cm,%
                % top=44.27pt,%
                includehead,%
                headheight=13.6pt]{geometry}
\renewcommand{\rmdefault}{ptm}
\renewcommand{\figurename}{Fig.}
\newcommand{\longpage}[1][1]{\enlargethispage{#1\baselineskip}}
\newcommand{\shortpage}[1][1]{\enlargethispage{-#1\baselineskip}}
\def\alice@no@spacing{%
  \itemsep=2pt\topsep=0pt\partopsep=0pt\parskip=0pt\parsep=0pt}
\newenvironment{Itemize}{%
  \vspace*{-1mm}\itemize\alice@no@spacing}{%
  \endlist\vspace*{-1mm}}
\newenvironment{Enumerate}{%
  \vspace*{-1mm}\enumerate\alice@no@spacing}{%
  \endlist\vspace*{-1mm}}
\newenvironment{Description}{%
  \vspace*{-1mm}\description\alice@no@spacing}{%
  \endlist\vspace*{-1mm}}
\renewcommand\labelitemi{\normalfont\textendash}
\renewcommand\labelitemii{\normalfont\textendash}
\renewcommand\labelitemiii{\textasteriskcentered}
\renewcommand\labelitemiv{\textperiodcentered}
\setlength\abovecaptionskip{4pt \@plus 2pt \@minus 2pt}
\setlength\belowcaptionskip{4pt \@plus 2pt \@minus 2pt}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{{\small\textbf{#1:} #2}}%
  \ifdim \wd\@tempboxa >\hsize
    {\small\textbf{#1:} #2}\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
\renewcommand\section{\@startsection{section}{1}{\z@}%
  {-2.5ex \@plus -1ex \@minus -.2ex}%
  {0.3ex \@plus.2ex}%
  {\large\bfseries\mathversion{bold}}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
  {-1.0ex\@plus -1ex \@minus -.2ex}%
  {0.05ex \@plus .2ex}%
  {\normalsize\bfseries\mathversion{bold}}}%
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
  {-1.0ex\@plus -1ex \@minus -.2ex}%
  {0.05ex \@plus .2ex}%
  {\normalsize\itshape\bfseries\mathversion{bold}}}%
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
  {-1.0ex\@plus -1ex \@minus -.2ex}%
  {-1em}%
  {\normalsize\bfseries\mathversion{bold}}}%
\renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
  {-1.0ex\@plus -1ex \@minus -.2ex}%
  {-1em}%
  {\normalsize\itshape\bfseries\mathversion{bold}}}%
\iffalse
\renewcommand\appendix{\par
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \setcounter{figure}{0}%
  \setcounter{table}{0}%
  \setcounter{equation}{0}%
  \gdef\thesection{\@Alph\c@section}%
  \gdef\thesubsection{\thesection.\@arabic\c@subsection}
  \gdef\theequation{\thesection.\@arabic\c@equation}
  \@addtoreset {equation}{section}%
  \gdef\thetable{\thesection.\@arabic\c@table}%
  \@addtoreset {table}{section}%
  \gdef\thefigure{\thesection.\@arabic\c@figure}%
  \@addtoreset {figure}{section}}
\fi
\renewenvironment{thebibliography}[1]{%
  \section*{\refname}%
  \@mkboth{\refname}{\refname}%
  \list{\@biblabel{\@arabic\c@enumiv}}%
  {\settowidth\labelwidth{\@biblabel{#1}}%
    \leftmargin\labelwidth
    \advance\leftmargin\labelsep
    \@openbib@code
    \usecounter{enumiv}%
    \let\p@enumiv\@empty
    \renewcommand\theenumiv{\@arabic\c@enumiv}}%
  \alice@no@spacing
  \sloppy
  \clubpenalty4000
  \@clubpenalty \clubpenalty
  \widowpenalty4000%
  \sfcode`\.\@m}{%
  \def\@noitemerr
  {\@latex@warning{Empty `thebibliography' environment}}%
  \endlist}
\newcommand\notenumber[1]{\gdef\alice@notenumber{#1}}
\let\alice@notenumber\@undefined
\def\alice@collaboration{The ALICE Collaboration}
\def\alice@short@collaboration{The ALICE Collaboration}
\newcommand\acknowledgementsfile[1]{%
  \gdef\alice@acknowledgemenst@file{#1}}
\newcommand\authorlistfile[1]{%
  \gdef\alice@authorlist@file{#1}}
\acknowledgementsfile{acknowledgements.tex}
\authorlistfile{authorlist.tex}
\if@alice@public@note
  \newenvironment{acknowledgments}{%
    \let\\=\@empty{}%
    \parskip=.5ex%
    \parindent=0pt\relax}{\relax}
  \@alice@many@uthorstrue
  \@alice@compattrue
  \AtEndDocument{%
    \section{Acknowledgments}
    \label{app:acknowledgments}
    \begin{acknowledgments}
      \input{\alice@acknowledgemenst@file}
    \end{acknowledgments}
    \section{\alice@collaboration}\label{app:collab}
    \input{\alice@authorlist@file}}
\fi
\newcommand\submitted[1]{\gdef\alice@submitted{#1}}
\let\alice@submitted\@undefined
\newcommand\collaboration[1]{\gdef\alice@user@collaboration{#1}}
\let\alice@user@collaboration\@undefined
\newcommand\conference[1]{\gdef\alice@conference{#1}}
\let\alice@conference\@undefined
\newcommand\titlenote[1]{\gdef\alice@title@note{#1}}
\let\alice@title@note\@undefined
\newcommand\dedication[1]{\gdef\alice@dedication{#1}}
\let\alice@dedication\@undefined
\newcommand\affiliation[1]{(#1)}
\renewcommand\author[2][]{%
  \ifx|#1|%
  \else%
    \gdef\alice@short@author{#1}%
  \fi%
  \gdef\@author{#2}}
\let\alice@short@author\@undefined
\newif\if@alice@institute@in@list\@alice@institute@in@listfalse
\newif\if@alice@institute@out@list\@alice@institute@out@listfalse
\newcommand\institute[1]{\alice@institute{alice@LI@#1}}
\newcounter{alice@inst}
\def\alice@al@test{authorlist}%
\def\alice@institute#1{%
  % \message{alice@institute arguments #1}%
  \ifx\alice@al@test\@currenvir%
    \@alice@institute@in@listtrue%
    \if@alice@institute@out@list%
      \@latex@error{Institute definitions in lists does not go with %
        definitions given out of list}%
    \fi%
    % \message{Label an enumeration with #1}
    \label{#1}%
    \let\next@a\@gobble%
  \else%
    % \message{Institute definition out of list}
    \@alice@institute@out@listtrue%
    \if@alice@institute@in@list%
      \@latex@error{Institute definitions outside lists does not go with %
        definitions given in lists}%
    \fi%
    \refstepcounter{alice@inst}\label{#1}%
    \xdef\@thefnmark{\arabic{alice@inst}}%
    \let\next@a\alice@inst@foot%
  \fi%
  \next@a{#1}}
\def\alice@inst@foot#1#2{%
  % \message{Arguments #1 #2}
  \if@alice@many@uthors%
    % \message{Typesetting foot note text}%
    \@footnotetext{#2}%
  \else%
    % \message{Defining \csname f@#1\endcsname}
    \expandafter\gdef\csname f@#1\endcsname{#2}%
  \fi}
\def\alice@@institute#1{%
  \\%
  \emph{#1}\\[1ex]}
\gdef\alice@nnospacing{%
  \labelsep=0pt%
  \itemsep=0pt%
  \topsep=0pt%
  \partopsep=0pt%
  \parskip=0pt%
  \parsep=0pt}
\newenvironment{authorlist}{%
  \if@alice@many@uthors
    % \vspace*{-1.5ex}
    \enumerate\alice@nnospacing%
  \else%
    \begin{center}%
  \fi}{%
  \if@alice@many@uthors%
    \endlist%
    %\vspace*{-1.5ex}%
  \else%
    \end{center}%
  \fi}

\newcommand\institutereferenceend[1]{}
\def\alice@inst@endchar{)}
\newcommand\inst[2][]{%
  \if@alice@many@uthors%
    \textsuperscript{\alice@inst[#1]{#2}\alice@inst@endchar}%
  \else%
    \ifx|#1|\else%
      \textsuperscript{\alice@@anot{#1}}%
    \fi% Annotation!
    \ifx|#2|\relax\else%
      \\%
      \emph{\alice@@@inst{#2}}\\[1ex]%
      \def\and{\def\and{, }}%
    \fi%
  \fi%
}
\def\alice@inst[#1]#2{%
  \ifx|#2|\relax%
  \else%
    \def\temp@a{\def\temp@a{,}}% Comma separator
    \@for\@el:=#2\do{\temp@a\alice@@inst{\@el}}\fi%
  \ifx|#1|\relax%
  \else%
    \ifx|#2|\relax%
    \else%
      ,%
    \fi%
    \alice@@anot{#1}%
  \fi}
\def\alice@@inst#1{\@ifundefined{r@alice@LI@#1}{%
    \textbf{??}}{%
    \edef\@tempa{\@nameuse{r@alice@LI@#1}}%
    \expandafter\@car\@tempa\@nil}}
\def\alice@@anot#1{%
  \@ifundefined{r@alice@LA@#1}{%
    % \message{Label r@alice@LA@#1 not defined}%
    #1}{%
    \edef\@tempa{\@nameuse{r@alice@LA@#1}}
    \expandafter\@car\@tempa\@nil}}
\def\alice@@@inst#1{%
  \def\temp@a{\def\temp@a{,}}% Comma separator
  \@for\@el:=#1\do{\temp@a\ \csname f@alice@LI@\@el\endcsname}}%
\def\alice@annot#1{#1}
\newif\if@alice@annot@in@list\@alice@annot@in@listfalse
\newif\if@alice@annot@out@list\@alice@annot@out@listfalse
\let\alice@old@makefnmark\@makefnmark
\let\alice@old@makefntext\@makefntext
\newcommand\annotation[1]{%
  \ifx\alice@al@test\@currenvir%  in a list
    \@alice@annot@in@listtrue%
    \if@alice@annot@out@list%
      \typeout{Annotation definitions in lists does not go with %
        definitions given out of list}%
    \fi%
    % \message{Defining label alice@la@#1}%
    \label{alice@LA@#1}%
    \let\next@a\relax%
  \else
    \@alice@annot@out@listtrue%
    \if@alice@annot@in@list%
      \typeout{Annotation definitions outside lists does not go with %
        definitions given in lists}%
    \fi%
    \refstepcounter{alice@inst}\label{#1}%
    \xdef\@thefnmark{#1}%
    \let\next@a\@footnotetext%
   \fi%
   \next@a}
\newenvironment{annotations}[1]{%
  \vskip2mm\bgroup
  % \def\@makefnmark{\hbox{$^{\mathrm{\@thefnmark)}}$}}
  \long\def\@makefntext##1{\parindent=\z@
    \begin{list}{}{\labelwidth=#1\leftmargin\labelwidth% set global width
        \labelsep3pt%  horizontal separation between term and item
        \itemsep0pt%   vertical separation between two items
        \parsep0pt%    vertical separation two paragraphs in an item
        \topsep0pt%    vertical separation text/list
        \footnotesize}%(
    \item [\hfill \textsuperscript{\@thefnmark)}] ##1\end{list}}
}{\egroup
  \let\@makefntext\alice@old@makefntext}
\if@alice@compat
  \def\Idef#1{\institute{#1}}
  \def\Iref#1{\inst[]{#1}}
  \def\Instfoot#1#2{\institute{#1}{#2}}
  \def\Aref#1{\inst[#1]{}}
  \let\Arefs\Aref
  \let\Anotfoot\annotation
  \def\AAref#1#2{\inst[#1,#2]{}}\let\AnotAnotref\AAref%
  \def\IAref#1#2{\inst[#2]{#1}}\let\InstAnotref\IAref%
  \def\IIref#1#2{\inst{#1,#2}}\let\InstInstref\IIref%
  \def\IIAref#1#2#3{\inst[#3]{#1,#2}}\let\InstInstAnotref\IIAref%
  \def\IIIref#1#2#3{\inst{#1,#2,#3}}\let\InstInstInstref\IIIref%
  \let\Irefn\Iref%
  \let\Irefnn\IIref
  \let\Adef\annotation
  \let\Deceased\deceased
  \let\Submitted\submitted
  \let\Collaboration\collaboration
  \let\Conference\conference
  \let\Dedication\dedication
  \let\Note\titlenote
  \let\And\and
  \if@alice@many@uthors\gdef\And{, }\fi
  \newenvironment{Authlist}{%
    \begin{authorlist}}{\end{authorlist}}
  % \show\Irefn
\fi

\long\def\alice@maketitle{%
  \begingroup%
  \if@alice@public@note
    \typeout{Overriding author for public notes}
    \gdef\@author{\alice@collaboration\thanks{See \appendixname~\ref{app:collab}
        for list of collaboration members}}%
    \global\let\alice@short@author\alice@collaboration%
  \fi
  \@ifundefined{alice@notenumber}{%
    \typeout{No ALICE note number defined}
    \def\alice@notenumber{????--???}}{}
  \edef\@alice@note@id{%
    ALICE--\if@alice@public@note PUBLIC\else ANA\fi--\alice@notenumber}
  \begin{tabular}[t]{p{.5\textwidth}@{}p{.5\textwidth}}
    \multicolumn{2}{c}{%
        {\LARGE EUROPEAN ORGANIZATION FOR NUCLEAR RESEARCH}}\\[4ex]
    % \flushleft%
    {\includegraphics[keepaspectratio,height=2cm]{alice_logo}\hfill}
    &
    % \flushright%
    {\hfill\includegraphics[keepaspectratio,height=2cm]{cernlogo}}
    \\
    &\hfill\@alice@note@id\\
    &\hfill\@date
  \end{tabular}
  \hrulefill
  \mbox{}\\[4ex]
  {\Large\bf\@title}\\[4ex]
  \global\let\alice@title\@title
  % \xdef\alice@title{\@title}%
  %\show\alice@title
}
\def\@maketitle{%
  \newpage
  \null
  % \vskip 2em%
  \begin{center}
    \alice@maketitle
    \vskip1.5em%
    \begingroup%
    \let\instiution\alice@@institute%
    \def\and{,\ }%
    % \begin{tabular}[t]{c}
      \@author
    % \end{tabular}
    \endgroup
    \par
    \vskip1.5em
    \@ifundefined{alice@user@collaboration}{}{%
      \alice@user@collaboration\\[2ex]}
    \@ifundefined{alice@submitted}{}{%
      \textit{\alice@submitted}\\[2ex]}
    \@ifundefined{alice@conference}{}{%
      \textit{\alice@conference}\\[2ex]}
    \@ifundefined{alice@title@note}{}{%
      \textit{\alice@title@note}\\[2ex]}
    \@ifundefined{alice@dedication}{}{%
      \textit{\alice@dedication}\\[2ex]}
 \end{center}
  %\xdef\temp@alice@author{\@author}
  \@ifundefined{alice@short@author}{%
    \let\alice@short@author\@author}{}
  \gdef\alice@author{\begingroup%
    \let\thanks\@gobble\let\@thanks\@empty\alice@short@author\endgroup}%
  \gdef\ps@plain{\let\@oddhead\@empty%
    \let\@oddfoot\@empty%
    \let\@evenfoot\@empty%
    \let\@evenhead\@empty}
  \global\let\@makefnmark\alice@old@makefnmark}
\if@alice@public@note
  \def\alice@inst@endchar{}
\fi
\newcommand{\bra}[1]{\ifmmode\left\langle #1\right|\else%
  \@latex@error{Macro \protect\bra cannot be used outside math mode}\fi}
\newcommand{\ket}[1]{\ifmmode\left| #1\right\rangle\else%
  \@latex@error{Macro \protect\ket cannot be used outside math mode}\fi}
\newcommand{\braket}[2]{\ifmmode\left\langle #1 | #2\right\rangle\else%
  \@latex@error{Macro \protect\braket cannot be used outside math mode}\fi}
\setcounter{topnumber}{5}
\setcounter{totalnumber}{6}
\renewcommand{\topfraction}{0.99}
\renewcommand{\dbltopfraction}{0.99}
\renewcommand{\bottomfraction}{0.99}
\renewcommand{\textfraction}{0.15}
\setlength{\intextsep}{10pt plus 2pt minus 4pt}
\setlength{\floatsep}{10pt plus 2pt minus 4pt}
\setlength{\textfloatsep}{10pt plus 2pt minus 4pt}
\setlength{\dblfloatsep}{10pt plus 2pt minus 4pt}
\setlength{\dbltextfloatsep}{10pt plus 2pt minus 4pt}
\setlength{\parskip}{.5\baselineskip}
\setlength{\parindent}{0pt}
\allowdisplaybreaks
\def\ps@alice@note{%
  \let\@oddfoot\@empty%
  \let\@evenfoot\@empty%
  \def\@evenhead{%
    %\show\alice@title%
    \begin{tabular}[t]{@{}p{.15\textwidth}@{}p{.85\textwidth}@{}}%
      \thepage\hfill&\truncate{.85\textwidth}{\hfill\alice@author}\\%
      \hline%
    \end{tabular}}
  \def\@oddhead{%
    \begin{tabular}[t]{@{}p{.85\textwidth}@{}p{.15\textwidth}@{}}%
      \truncate{.85\textwidth}{\alice@title\hfill}&\hfill\thepage\\%
      \hline%
    \end{tabular}}}
\let\ps@headings\ps@alice@note
\pagestyle{headings}
\endinput
%%
%% End of file `alicenote.cls'.
