





<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">



  <link crossorigin="anonymous" href="https://assets-cdn.github.com/assets/frameworks-4ba00b1aa0227e4b7a7961544c3b7938afb2720757a471735991ec4475c829e0.css" integrity="sha256-S6ALGqAifkt6eWFUTDt5OK+ycgdXpHFzWZHsRHXIKeA=" media="all" rel="stylesheet" />
  <link crossorigin="anonymous" href="https://assets-cdn.github.com/assets/github-81684fd7510dbd0631cf199316b86c30d3741f8600001e7998ec8ec766d6423c.css" integrity="sha256-gWhP11ENvQYxzxmTFrhsMNN0H4YAAB55mOyOx2bWQjw=" media="all" rel="stylesheet" />
  
  
  
  

  <meta name="viewport" content="width=device-width">
  
  <title>AliPhysics/AddTask_dNdPtpp_TPCITS.C at master · alisw/AliPhysics</title>
  <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub">
  <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub">
  <meta property="fb:app_id" content="1401488693436528">

    
    <meta content="https://avatars1.githubusercontent.com/u/12711103?v=3&amp;s=400" property="og:image" /><meta content="GitHub" property="og:site_name" /><meta content="object" property="og:type" /><meta content="alisw/AliPhysics" property="og:title" /><meta content="https://github.com/alisw/AliPhysics" property="og:url" /><meta content="AliPhysics - ALICE Analysis Repository" property="og:description" />

  <link rel="assets" href="https://assets-cdn.github.com/">
  <link rel="web-socket" href="wss://live.github.com/_sockets/VjI6MTgyMTA0ODE0OjQ3NzAzYzc1YTEyMjk1OTI1MjFlZjg2MjJiZGNhNzlhOTEwNjIyMDJmMWU3Y2NlMDJiYWEwMDViMDZhNTc3NzE=--13d070a2a6fab074672bc810f599179138cba1f6">
  <meta name="pjax-timeout" content="1000">
  <link rel="sudo-modal" href="/sessions/sudo_modal">
  <meta name="request-id" content="9208:3F90:39AA1DF:52B97EE:59491288" data-pjax-transient>
  

  <meta name="selected-link" value="repo_source" data-pjax-transient>

  <meta name="google-site-verification" content="KT5gs8h0wvaagLKAVWq8bbeNwnZZK1r1XQysX3xurLU">
<meta name="google-site-verification" content="ZzhVyEFwb7w3e0-uOTltm8Jsck2F5StVihD0exw2fsA">
    <meta name="google-analytics" content="UA-3769691-2">

<meta content="collector.githubapp.com" name="octolytics-host" /><meta content="github" name="octolytics-app-id" /><meta content="https://collector.githubapp.com/github-external/browser_event" name="octolytics-event-url" /><meta content="9208:3F90:39AA1DF:52B97EE:59491288" name="octolytics-dimension-request_id" /><meta content="iad" name="octolytics-dimension-region_edge" /><meta content="iad" name="octolytics-dimension-region_render" /><meta content="26876110" name="octolytics-actor-id" /><meta content="mkruegerGitHub" name="octolytics-actor-login" /><meta content="5591bf415a50af885bcb506cfa23a193f538160460949474c4ea008376a96a52" name="octolytics-actor-hash" />
<meta content="/&lt;user-name&gt;/&lt;repo-name&gt;/blob/show" data-pjax-transient="true" name="analytics-location" />




  <meta class="js-ga-set" name="dimension1" content="Logged In">


  

      <meta name="hostname" content="github.com">
  <meta name="user-login" content="mkruegerGitHub">

      <meta name="expected-hostname" content="github.com">
    <meta name="js-proxy-site-detection-payload" content="MzE3MjgxYTA4N2EyMjk1NWFjYWExZTMyZGQ5ZmMzYmE3NTU0ZDRlNWIyOTlhNGU1MDU3Y2E3NjhkMzg1MmYzZHx7InJlbW90ZV9hZGRyZXNzIjoiMTQxLjIuMjQzLjYxIiwicmVxdWVzdF9pZCI6IjkyMDg6M0Y5MDozOUFBMURGOjUyQjk3RUU6NTk0OTEyODgiLCJ0aW1lc3RhbXAiOjE0OTc5NjEwOTcsImhvc3QiOiJnaXRodWIuY29tIn0=">


  <meta name="html-safe-nonce" content="6d0c0e96f316ea519071f245ad4c8a18bdbbff19">

  <meta http-equiv="x-pjax-version" content="f6270ca1a0cb6c87d67ba379b194ad71">
  

      <link href="https://github.com/alisw/AliPhysics/commits/master.atom" rel="alternate" title="Recent Commits to AliPhysics:master" type="application/atom+xml">

  <meta name="description" content="AliPhysics - ALICE Analysis Repository">
  <meta name="go-import" content="github.com/alisw/AliPhysics git https://github.com/alisw/AliPhysics.git">

  <meta content="12711103" name="octolytics-dimension-user_id" /><meta content="alisw" name="octolytics-dimension-user_login" /><meta content="61661378" name="octolytics-dimension-repository_id" /><meta content="alisw/AliPhysics" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="false" name="octolytics-dimension-repository_is_fork" /><meta content="61661378" name="octolytics-dimension-repository_network_root_id" /><meta content="alisw/AliPhysics" name="octolytics-dimension-repository_network_root_nwo" /><meta content="false" name="octolytics-dimension-repository_explore_github_marketplace_ci_cta_shown" />


    <link rel="canonical" href="https://github.com/alisw/AliPhysics/blob/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C" data-pjax-transient>


  <meta name="browser-stats-url" content="https://api.github.com/_private/browser/stats">

  <meta name="browser-errors-url" content="https://api.github.com/_private/browser/errors">

  <link rel="mask-icon" href="https://assets-cdn.github.com/pinned-octocat.svg" color="#000000">
  <link rel="icon" type="image/x-icon" href="https://assets-cdn.github.com/favicon.ico">

<meta name="theme-color" content="#1e2327">



  </head>

  <body class="logged-in env-production page-blob">
    



  <div class="position-relative js-header-wrapper ">
    <a href="#start-of-content" tabindex="1" class="bg-black text-white p-3 show-on-focus js-skip-to-content">Skip to content</a>
    <div id="js-pjax-loader-bar" class="pjax-loader-bar"><div class="progress"></div></div>

    
    
    



        
<div class="header" role="banner">
  <div class="container clearfix">
    <a class="header-logo-invertocat" href="https://github.com/" data-hotkey="g d" aria-label="Homepage" data-ga-click="Header, go to dashboard, icon:logo">
  <svg aria-hidden="true" class="octicon octicon-mark-github" height="32" version="1.1" viewBox="0 0 16 16" width="32"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"/></svg>
</a>


        <div class="header-search scoped-search site-scoped-search js-site-search" role="search">
  <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="/alisw/AliPhysics/search" class="js-site-search-form" data-scoped-search-url="/alisw/AliPhysics/search" data-unscoped-search-url="/search" method="get"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /></div>
    <label class="form-control header-search-wrapper js-chromeless-input-container">
        <a href="/alisw/AliPhysics/blob/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C" class="header-search-scope no-underline">This repository</a>
      <input type="text"
        class="form-control header-search-input js-site-search-focus js-site-search-field is-clearable"
        data-hotkey="s"
        name="q"
        value=""
        placeholder="Search"
        aria-label="Search this repository"
        data-unscoped-placeholder="Search GitHub"
        data-scoped-placeholder="Search"
        autocapitalize="off">
        <input type="hidden" class="js-site-search-type-field" name="type" >
    </label>
</form></div>


      <ul class="header-nav float-left" role="navigation">
        <li class="header-nav-item">
          <a href="/pulls" aria-label="Pull requests you created" class="js-selected-navigation-item header-nav-link" data-ga-click="Header, click, Nav menu - item:pulls context:user" data-hotkey="g p" data-selected-links="/pulls /pulls/assigned /pulls/mentioned /pulls">
            Pull requests
</a>        </li>
        <li class="header-nav-item">
          <a href="/issues" aria-label="Issues you created" class="js-selected-navigation-item header-nav-link" data-ga-click="Header, click, Nav menu - item:issues context:user" data-hotkey="g i" data-selected-links="/issues /issues/assigned /issues/mentioned /issues">
            Issues
</a>        </li>
            <li class="header-nav-item">
              <a href="/marketplace" class="js-selected-navigation-item header-nav-link" data-ga-click="Header, click, Nav menu - item:marketplace context:user" data-selected-links=" /marketplace">
                Marketplace
</a>            </li>
          <li class="header-nav-item">
            <a class="header-nav-link" href="https://gist.github.com/" data-ga-click="Header, go to gist, text:gist">Gist</a>
          </li>
      </ul>

    
<ul class="header-nav user-nav float-right" id="user-links">
  <li class="header-nav-item">
    

  </li>

  <li class="header-nav-item dropdown js-menu-container">
    <a class="header-nav-link tooltipped tooltipped-s js-menu-target" href="/new"
       aria-label="Create new…"
       aria-expanded="false"
       aria-haspopup="true"
       data-ga-click="Header, create new, icon:add">
      <svg aria-hidden="true" class="octicon octicon-plus float-left" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 9H7v5H5V9H0V7h5V2h2v5h5z"/></svg>
      <span class="dropdown-caret"></span>
    </a>

    <div class="dropdown-menu-content js-menu-content">
      <ul class="dropdown-menu dropdown-menu-sw">
        
<a class="dropdown-item" href="/new" data-ga-click="Header, create new repository">
  New repository
</a>

  <a class="dropdown-item" href="/new/import" data-ga-click="Header, import a repository">
    Import repository
  </a>

<a class="dropdown-item" href="https://gist.github.com/" data-ga-click="Header, create new gist">
  New gist
</a>

  <a class="dropdown-item" href="/organizations/new" data-ga-click="Header, create new organization">
    New organization
  </a>



  <div class="dropdown-divider"></div>
  <div class="dropdown-header">
    <span title="alisw/AliPhysics">This repository</span>
  </div>
    <a class="dropdown-item" href="/alisw/AliPhysics/issues/new" data-ga-click="Header, create new issue">
      New issue
    </a>

      </ul>
    </div>
  </li>

  <li class="header-nav-item dropdown js-menu-container">
    <a class="header-nav-link name tooltipped tooltipped-sw js-menu-target" href="/mkruegerGitHub"
       aria-label="View profile and more"
       aria-expanded="false"
       aria-haspopup="true"
       data-ga-click="Header, show menu, icon:avatar">
      <img alt="@mkruegerGitHub" class="avatar" src="https://avatars0.githubusercontent.com/u/26876110?v=3&amp;s=40" height="20" width="20">
      <span class="dropdown-caret"></span>
    </a>

    <div class="dropdown-menu-content js-menu-content">
      <div class="dropdown-menu dropdown-menu-sw">
        <div class="dropdown-header header-nav-current-user css-truncate">
          Signed in as <strong class="css-truncate-target">mkruegerGitHub</strong>
        </div>

        <div class="dropdown-divider"></div>

        <a class="dropdown-item" href="/mkruegerGitHub" data-ga-click="Header, go to profile, text:your profile">
          Your profile
        </a>
        <a class="dropdown-item" href="/mkruegerGitHub?tab=stars" data-ga-click="Header, go to starred repos, text:your stars">
          Your stars
        </a>
        <a class="dropdown-item" href="/explore" data-ga-click="Header, go to explore, text:explore">
          Explore
        </a>
        <a class="dropdown-item" href="https://help.github.com" data-ga-click="Header, go to help, text:help">
          Help
        </a>

        <div class="dropdown-divider"></div>

        <a class="dropdown-item" href="/settings/profile" data-ga-click="Header, go to settings, icon:settings">
          Settings
        </a>

        <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="/logout" class="logout-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="o4m2yNdT5yVFKlQJ9NW1L0QzVEbx0OIvG6t2cAhdJlK1WtxbNZE5MiFyj3I7i8il4uBmUploOBMZXo42kSx6Vw==" /></div>
          <button type="submit" class="dropdown-item dropdown-signout" data-ga-click="Header, sign out, icon:logout">
            Sign out
          </button>
</form>      </div>
    </div>
  </li>
</ul>


    <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="/logout" class="sr-only right-0" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="o0stuM2r4bVUII2oR2kv850q944IdHabX43ZAEHnTjW1mEcrL2k/ojB4VtOIN1J5O/nFmmDMrKddeCFG2JYSMA==" /></div>
      <button type="submit" class="dropdown-item dropdown-signout" data-ga-click="Header, sign out, icon:logout">
        Sign out
      </button>
</form>  </div>
</div>


      

  </div>

  <div id="start-of-content" class="show-on-focus"></div>

    <div id="js-flash-container">
</div>



  <div role="main">
        <div itemscope itemtype="http://schema.org/SoftwareSourceCode">
    <div id="js-repo-pjax-container" data-pjax-container>
      



  


    <div class="pagehead repohead instapaper_ignore readability-menu experiment-repo-nav">
      <div class="container repohead-details-container">

        <ul class="pagehead-actions">
  <li>
        <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="/notifications/subscribe" class="js-social-container" data-autosubmit="true" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="Fj2uOJGHa/28pJjYYlIZPG1yz6QlmrT3nDQgfqn10GC3/GluhwirPtXtFo8F86VpbVA4Fm5XazT0goKFzGPikg==" /></div>      <input class="form-control" id="repository_id" name="repository_id" type="hidden" value="61661378" />

        <div class="select-menu js-menu-container js-select-menu">
          <a href="/alisw/AliPhysics/subscription"
            class="btn btn-sm btn-with-count select-menu-button js-menu-target"
            role="button"
            aria-haspopup="true"
            aria-expanded="false"
            aria-label="Toggle repository notifications menu"
            data-ga-click="Repository, click Watch settings, action:blob#show">
            <span class="js-select-button">
                <svg aria-hidden="true" class="octicon octicon-eye" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M8.06 2C3 2 0 8 0 8s3 6 8.06 6C13 14 16 8 16 8s-3-6-7.94-6zM8 12c-2.2 0-4-1.78-4-4 0-2.2 1.8-4 4-4 2.22 0 4 1.8 4 4 0 2.22-1.78 4-4 4zm2-4c0 1.11-.89 2-2 2-1.11 0-2-.89-2-2 0-1.11.89-2 2-2 1.11 0 2 .89 2 2z"/></svg>
                Watch
            </span>
          </a>
            <a class="social-count js-social-count"
              href="/alisw/AliPhysics/watchers"
              aria-label="16 users are watching this repository">
              16
            </a>

        <div class="select-menu-modal-holder">
          <div class="select-menu-modal subscription-menu-modal js-menu-content">
            <div class="select-menu-header js-navigation-enable" tabindex="-1">
              <svg aria-label="Close" class="octicon octicon-x js-menu-close" height="16" role="img" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"/></svg>
              <span class="select-menu-title">Notifications</span>
            </div>

              <div class="select-menu-list js-navigation-container" role="menu">

                <div class="select-menu-item js-navigation-item selected" role="menuitem" tabindex="0">
                  <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
                  <div class="select-menu-item-text">
                    <input checked="checked" id="do_included" name="do" type="radio" value="included" />
                    <span class="select-menu-item-heading">Not watching</span>
                    <span class="description">Be notified when participating or @mentioned.</span>
                    <span class="js-select-button-text hidden-select-button-text">
                      <svg aria-hidden="true" class="octicon octicon-eye" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M8.06 2C3 2 0 8 0 8s3 6 8.06 6C13 14 16 8 16 8s-3-6-7.94-6zM8 12c-2.2 0-4-1.78-4-4 0-2.2 1.8-4 4-4 2.22 0 4 1.8 4 4 0 2.22-1.78 4-4 4zm2-4c0 1.11-.89 2-2 2-1.11 0-2-.89-2-2 0-1.11.89-2 2-2 1.11 0 2 .89 2 2z"/></svg>
                      Watch
                    </span>
                  </div>
                </div>

                <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
                  <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
                  <div class="select-menu-item-text">
                    <input id="do_subscribed" name="do" type="radio" value="subscribed" />
                    <span class="select-menu-item-heading">Watching</span>
                    <span class="description">Be notified of all conversations.</span>
                    <span class="js-select-button-text hidden-select-button-text">
                      <svg aria-hidden="true" class="octicon octicon-eye" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M8.06 2C3 2 0 8 0 8s3 6 8.06 6C13 14 16 8 16 8s-3-6-7.94-6zM8 12c-2.2 0-4-1.78-4-4 0-2.2 1.8-4 4-4 2.22 0 4 1.8 4 4 0 2.22-1.78 4-4 4zm2-4c0 1.11-.89 2-2 2-1.11 0-2-.89-2-2 0-1.11.89-2 2-2 1.11 0 2 .89 2 2z"/></svg>
                        Unwatch
                    </span>
                  </div>
                </div>

                <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
                  <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
                  <div class="select-menu-item-text">
                    <input id="do_ignore" name="do" type="radio" value="ignore" />
                    <span class="select-menu-item-heading">Ignoring</span>
                    <span class="description">Never be notified.</span>
                    <span class="js-select-button-text hidden-select-button-text">
                      <svg aria-hidden="true" class="octicon octicon-mute" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M8 2.81v10.38c0 .67-.81 1-1.28.53L3 10H1c-.55 0-1-.45-1-1V7c0-.55.45-1 1-1h2l3.72-3.72C7.19 1.81 8 2.14 8 2.81zm7.53 3.22l-1.06-1.06-1.97 1.97-1.97-1.97-1.06 1.06L11.44 8 9.47 9.97l1.06 1.06 1.97-1.97 1.97 1.97 1.06-1.06L13.56 8l1.97-1.97z"/></svg>
                        Stop ignoring
                    </span>
                  </div>
                </div>

              </div>

            </div>
          </div>
        </div>
</form>
  </li>

  <li>
      <div class="js-toggler-container js-social-container starring-container ">
    <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="/alisw/AliPhysics/unstar" class="starred" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="HyVa3e2/G0/um1u/XZ57rqJuEr6o7quUXCgIJOSCkJ4F28wQM0ZfPbBVTQfQlPfajTtEIAqivCqDzGJnWx5dUA==" /></div>
      <button
        type="submit"
        class="btn btn-sm btn-with-count js-toggler-target"
        aria-label="Unstar this repository" title="Unstar alisw/AliPhysics"
        data-ga-click="Repository, click unstar button, action:blob#show; text:Unstar">
        <svg aria-hidden="true" class="octicon octicon-star" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path fill-rule="evenodd" d="M14 6l-4.9-.64L7 1 4.9 5.36 0 6l3.6 3.26L2.67 14 7 11.67 11.33 14l-.93-4.74z"/></svg>
        Unstar
      </button>
        <a class="social-count js-social-count" href="/alisw/AliPhysics/stargazers"
           aria-label="20 users starred this repository">
          20
        </a>
</form>
    <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="/alisw/AliPhysics/star" class="unstarred" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="CXA8+1Q6LiNE1xB7e4LIw93ZU3tETpM7Y9PydmhJUwK90obouwJGAahne3xukaQqFPg89Ux0VTXQp+EpSZgxWw==" /></div>
      <button
        type="submit"
        class="btn btn-sm btn-with-count js-toggler-target"
        aria-label="Star this repository" title="Star alisw/AliPhysics"
        data-ga-click="Repository, click star button, action:blob#show; text:Star">
        <svg aria-hidden="true" class="octicon octicon-star" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path fill-rule="evenodd" d="M14 6l-4.9-.64L7 1 4.9 5.36 0 6l3.6 3.26L2.67 14 7 11.67 11.33 14l-.93-4.74z"/></svg>
        Star
      </button>
        <a class="social-count js-social-count" href="/alisw/AliPhysics/stargazers"
           aria-label="20 users starred this repository">
          20
        </a>
</form>  </div>

  </li>

  <li>
          <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="/alisw/AliPhysics/fork" class="btn-with-count" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="5ueABYCWzrub26MNhjWzcftI+ErsWg+G6eE2GZd4wYzq+9ZKflpeock6BTsUIohQ2g5b6gp1Nc8R8HyHGP4lSg==" /></div>
            <button
                type="submit"
                class="btn btn-sm btn-with-count"
                data-ga-click="Repository, show fork modal, action:blob#show; text:Fork"
                title="Fork your own copy of alisw/AliPhysics to your account"
                aria-label="Fork your own copy of alisw/AliPhysics to your account">
              <svg aria-hidden="true" class="octicon octicon-repo-forked" height="16" version="1.1" viewBox="0 0 10 16" width="10"><path fill-rule="evenodd" d="M8 1a1.993 1.993 0 0 0-1 3.72V6L5 8 3 6V4.72A1.993 1.993 0 0 0 2 1a1.993 1.993 0 0 0-1 3.72V6.5l3 3v1.78A1.993 1.993 0 0 0 5 15a1.993 1.993 0 0 0 1-3.72V9.5l3-3V4.72A1.993 1.993 0 0 0 8 1zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3 10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3-10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"/></svg>
              Fork
            </button>
</form>
    <a href="/alisw/AliPhysics/network" class="social-count"
       aria-label="267 users forked this repository">
      267
    </a>
  </li>
</ul>

        <h1 class="public ">
  <svg aria-hidden="true" class="octicon octicon-repo" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M4 9H3V8h1v1zm0-3H3v1h1V6zm0-2H3v1h1V4zm0-2H3v1h1V2zm8-1v12c0 .55-.45 1-1 1H6v2l-1.5-1.5L3 16v-2H1c-.55 0-1-.45-1-1V1c0-.55.45-1 1-1h10c.55 0 1 .45 1 1zm-1 10H1v2h2v-1h3v1h5v-2zm0-10H2v9h9V1z"/></svg>
  <span class="author" itemprop="author"><a href="/alisw" class="url fn" rel="author">alisw</a></span><!--
--><span class="path-divider">/</span><!--
--><strong itemprop="name"><a href="/alisw/AliPhysics" data-pjax="#js-repo-pjax-container">AliPhysics</a></strong>

</h1>

      </div>
      <div class="container">
        
<nav class="reponav js-repo-nav js-sidenav-container-pjax"
     itemscope
     itemtype="http://schema.org/BreadcrumbList"
     role="navigation"
     data-pjax="#js-repo-pjax-container">

  <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
    <a href="/alisw/AliPhysics" class="js-selected-navigation-item selected reponav-item" data-hotkey="g c" data-selected-links="repo_source repo_downloads repo_commits repo_releases repo_tags repo_branches /alisw/AliPhysics" itemprop="url">
      <svg aria-hidden="true" class="octicon octicon-code" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path fill-rule="evenodd" d="M9.5 3L8 4.5 11.5 8 8 11.5 9.5 13 14 8 9.5 3zm-5 0L0 8l4.5 5L6 11.5 2.5 8 6 4.5 4.5 3z"/></svg>
      <span itemprop="name">Code</span>
      <meta itemprop="position" content="1">
</a>  </span>

    <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
      <a href="/alisw/AliPhysics/issues" class="js-selected-navigation-item reponav-item" data-hotkey="g i" data-selected-links="repo_issues repo_labels repo_milestones /alisw/AliPhysics/issues" itemprop="url">
        <svg aria-hidden="true" class="octicon octicon-issue-opened" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path fill-rule="evenodd" d="M7 2.3c3.14 0 5.7 2.56 5.7 5.7s-2.56 5.7-5.7 5.7A5.71 5.71 0 0 1 1.3 8c0-3.14 2.56-5.7 5.7-5.7zM7 1C3.14 1 0 4.14 0 8s3.14 7 7 7 7-3.14 7-7-3.14-7-7-7zm1 3H6v5h2V4zm0 6H6v2h2v-2z"/></svg>
        <span itemprop="name">Issues</span>
        <span class="Counter">0</span>
        <meta itemprop="position" content="2">
</a>    </span>

  <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
    <a href="/alisw/AliPhysics/pulls" class="js-selected-navigation-item reponav-item" data-hotkey="g p" data-selected-links="repo_pulls /alisw/AliPhysics/pulls" itemprop="url">
      <svg aria-hidden="true" class="octicon octicon-git-pull-request" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M11 11.28V5c-.03-.78-.34-1.47-.94-2.06C9.46 2.35 8.78 2.03 8 2H7V0L4 3l3 3V4h1c.27.02.48.11.69.31.21.2.3.42.31.69v6.28A1.993 1.993 0 0 0 10 15a1.993 1.993 0 0 0 1-3.72zm-1 2.92c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zM4 3c0-1.11-.89-2-2-2a1.993 1.993 0 0 0-1 3.72v6.56A1.993 1.993 0 0 0 2 15a1.993 1.993 0 0 0 1-3.72V4.72c.59-.34 1-.98 1-1.72zm-.8 10c0 .66-.55 1.2-1.2 1.2-.65 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"/></svg>
      <span itemprop="name">Pull requests</span>
      <span class="Counter">4</span>
      <meta itemprop="position" content="3">
</a>  </span>

    <a href="/alisw/AliPhysics/projects" class="js-selected-navigation-item reponav-item" data-selected-links="repo_projects new_repo_project repo_project /alisw/AliPhysics/projects">
      <svg aria-hidden="true" class="octicon octicon-project" height="16" version="1.1" viewBox="0 0 15 16" width="15"><path fill-rule="evenodd" d="M10 12h3V2h-3v10zm-4-2h3V2H6v8zm-4 4h3V2H2v12zm-1 1h13V1H1v14zM14 0H1a1 1 0 0 0-1 1v14a1 1 0 0 0 1 1h13a1 1 0 0 0 1-1V1a1 1 0 0 0-1-1z"/></svg>
      Projects
      <span class="Counter" >0</span>
</a>
    <a href="/alisw/AliPhysics/wiki" class="js-selected-navigation-item reponav-item" data-hotkey="g w" data-selected-links="repo_wiki /alisw/AliPhysics/wiki">
      <svg aria-hidden="true" class="octicon octicon-book" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M3 5h4v1H3V5zm0 3h4V7H3v1zm0 2h4V9H3v1zm11-5h-4v1h4V5zm0 2h-4v1h4V7zm0 2h-4v1h4V9zm2-6v9c0 .55-.45 1-1 1H9.5l-1 1-1-1H2c-.55 0-1-.45-1-1V3c0-.55.45-1 1-1h5.5l1 1 1-1H15c.55 0 1 .45 1 1zm-8 .5L7.5 3H2v9h6V3.5zm7-.5H9.5l-.5.5V12h6V3z"/></svg>
      Wiki
</a>

    <div class="reponav-dropdown js-menu-container">
      <button type="button" class="btn-link reponav-item reponav-dropdown js-menu-target " data-no-toggle aria-expanded="false" aria-haspopup="true">
        Insights
        <svg aria-hidden="true" class="octicon octicon-triangle-down v-align-middle text-gray" height="11" version="1.1" viewBox="0 0 12 16" width="8"><path fill-rule="evenodd" d="M0 5l6 6 6-6z"/></svg>
      </button>
      <div class="dropdown-menu-content js-menu-content">
        <div class="dropdown-menu dropdown-menu-sw">
          <a class="dropdown-item" href="/alisw/AliPhysics/pulse" data-skip-pjax>
            <svg aria-hidden="true" class="octicon octicon-pulse" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path fill-rule="evenodd" d="M11.5 8L8.8 5.4 6.6 8.5 5.5 1.6 2.38 8H0v2h3.6l.9-1.8.9 5.4L9 8.5l1.6 1.5H14V8z"/></svg>
            Pulse
          </a>
          <a class="dropdown-item" href="/alisw/AliPhysics/graphs" data-skip-pjax>
            <svg aria-hidden="true" class="octicon octicon-graph" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M16 14v1H0V0h1v14h15zM5 13H3V8h2v5zm4 0H7V3h2v10zm4 0h-2V6h2v7z"/></svg>
            Graphs
          </a>
        </div>
      </div>
    </div>
</nav>

      </div>
    </div>

<div class="container new-discussion-timeline experiment-repo-nav">
  <div class="repository-content">

    
  <a href="/alisw/AliPhysics/blob/3f05b9f721ee81991e99357be3ce172368940704/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C" class="d-none js-permalink-shortcut" data-hotkey="y">Permalink</a>

  <!-- blob contrib key: blob_contributors:v21:b324cd31f6200004da9c9fd56b267233 -->

  <div class="file-navigation js-zeroclipboard-container">
    
<div class="select-menu branch-select-menu js-menu-container js-select-menu float-left">
  <button class=" btn btn-sm select-menu-button js-menu-target css-truncate" data-hotkey="w"
    
    type="button" aria-label="Switch branches or tags" aria-expanded="false" aria-haspopup="true">
      <i>Branch:</i>
      <span class="js-select-button css-truncate-target">master</span>
  </button>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax>

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <svg aria-label="Close" class="octicon octicon-x js-menu-close" height="16" role="img" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"/></svg>
        <span class="select-menu-title">Switch branches/tags</span>
      </div>

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Filter branches/tags" id="context-commitish-filter-field" class="form-control js-filterable-field js-navigation-enable" placeholder="Filter branches/tags">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" data-filter-placeholder="Filter branches/tags" class="js-select-menu-tab" role="tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" data-filter-placeholder="Find a tag…" class="js-select-menu-tab" role="tab">Tags</a>
            </li>
          </ul>
        </div>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches" role="menu">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/EVE-dev/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="EVE-dev"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                EVE-dev
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/NanoAODdev/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="NanoAODdev"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                NanoAODdev
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/TPCrun3/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="TPCrun3"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                TPCrun3
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/TRDdev/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="TRDdev"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                TRDdev
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/aod-upgrade/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="aod-upgrade"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                aod-upgrade
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/dHLT/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="dHLT"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                dHLT
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/dmeson-backup/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="dmeson-backup"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                dmeson-backup
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/dmeson/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="dmeson"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                dmeson
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/emcal-embedding/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="emcal-embedding"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                emcal-embedding
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/feature-ALIROOT-5805/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="feature-ALIROOT-5805"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                feature-ALIROOT-5805
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/feature-alibits/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="feature-alibits"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                feature-alibits
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/feature-hltdev/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="feature-hltdev"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                feature-hltdev
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/feature-malturan/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="feature-malturan"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                feature-malturan
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/feature-onlineqa/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="feature-onlineqa"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                feature-onlineqa
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/feature-pwghf-dxhfe/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="feature-pwghf-dxhfe"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                feature-pwghf-dxhfe
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/feature-tpccalibonline/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="feature-tpccalibonline"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                feature-tpccalibonline
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/feature-tpcdev/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="feature-tpcdev"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                feature-tpcdev
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open selected"
               href="/alisw/AliPhysics/blob/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="master"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                master
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/patches/v5-06-28-01/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="patches/v5-06-28-01"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                patches/v5-06-28-01
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/pcm/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="pcm"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                pcm
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/pidcent/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="pidcent"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                pidcent
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/prod-hlt/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="prod-hlt"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                prod-hlt
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/rc/TEST-IGNORE-vAN-20151103/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="rc/TEST-IGNORE-vAN-20151103"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                rc/TEST-IGNORE-vAN-20151103
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/rc/TEST-IGNORE-vAN-20160321/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="rc/TEST-IGNORE-vAN-20160321"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                rc/TEST-IGNORE-vAN-20160321
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/rc/TEST-IGNORE-vAN-20170216/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="rc/TEST-IGNORE-vAN-20170216"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                rc/TEST-IGNORE-vAN-20170216
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/rc/rc2016-09-18/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="rc/rc2016-09-18"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                rc/rc2016-09-18
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/rc/vAN-20151229/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="rc/vAN-20151229"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                rc/vAN-20151229
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/rc/vAN-20160927/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="rc/vAN-20160927"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                rc/vAN-20160927
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/rc/vAN-20161127/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="rc/vAN-20161127"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                rc/vAN-20161127
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/rc/vAN-20161223/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="rc/vAN-20161223"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                rc/vAN-20161223
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/rc/vAN-20170119/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="rc/vAN-20170119"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                rc/vAN-20170119
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/relval-cvmfs-forceupdate/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="relval-cvmfs-forceupdate"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                relval-cvmfs-forceupdate
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/revert-1165-master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="revert-1165-master"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                revert-1165-master
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/v5-06-34-01-patches/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="v5-06-34-01-patches"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                v5-06-34-01-patches
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/v5-07-12-patches/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="v5-07-12-patches"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                v5-07-12-patches
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/v5-08-03-patches/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="v5-08-03-patches"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                v5-08-03-patches
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/v5-08-12-01-patches/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="v5-08-12-01-patches"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                v5-08-12-01-patches
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/v5-08-13-01-mcdev/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="v5-08-13-01-mcdev"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                v5-08-13-01-mcdev
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/alisw/AliPhysics/blob/v5-08-13-01-patches/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
               data-name="v5-08-13-01-patches"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                v5-08-13-01-patches
              </span>
            </a>
        </div>

          <div class="select-menu-no-results">Nothing to show</div>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170619/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170619"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170619">
                vAN-20170619
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170618/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170618"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170618">
                vAN-20170618
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170617/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170617"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170617">
                vAN-20170617
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170616/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170616"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170616">
                vAN-20170616
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170615/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170615"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170615">
                vAN-20170615
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170614/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170614"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170614">
                vAN-20170614
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170613/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170613"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170613">
                vAN-20170613
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170610/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170610"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170610">
                vAN-20170610
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170609/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170609"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170609">
                vAN-20170609
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170608/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170608"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170608">
                vAN-20170608
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170607/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170607"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170607">
                vAN-20170607
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170606/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170606"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170606">
                vAN-20170606
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170605/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170605"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170605">
                vAN-20170605
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170604/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170604"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170604">
                vAN-20170604
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170603/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170603"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170603">
                vAN-20170603
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170602/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170602"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170602">
                vAN-20170602
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170601/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170601"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170601">
                vAN-20170601
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170531/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170531"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170531">
                vAN-20170531
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170530/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170530"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170530">
                vAN-20170530
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170529/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170529"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170529">
                vAN-20170529
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170528/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170528"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170528">
                vAN-20170528
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170527/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170527"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170527">
                vAN-20170527
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170526/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170526"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170526">
                vAN-20170526
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170525/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170525"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170525">
                vAN-20170525
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170524/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170524"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170524">
                vAN-20170524
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170523/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170523"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170523">
                vAN-20170523
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170522/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170522"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170522">
                vAN-20170522
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170521/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170521"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170521">
                vAN-20170521
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170520/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170520"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170520">
                vAN-20170520
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170519/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170519"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170519">
                vAN-20170519
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170518/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170518"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170518">
                vAN-20170518
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170517/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170517"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170517">
                vAN-20170517
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170516/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170516"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170516">
                vAN-20170516
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170515/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170515"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170515">
                vAN-20170515
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170514/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170514"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170514">
                vAN-20170514
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170513/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170513"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170513">
                vAN-20170513
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170512/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170512"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170512">
                vAN-20170512
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170511/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170511"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170511">
                vAN-20170511
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170510/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170510"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170510">
                vAN-20170510
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170509/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170509"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170509">
                vAN-20170509
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170508/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170508"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170508">
                vAN-20170508
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170507/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170507"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170507">
                vAN-20170507
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170506/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170506"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170506">
                vAN-20170506
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170505/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170505"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170505">
                vAN-20170505
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170504/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170504"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170504">
                vAN-20170504
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170503/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170503"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170503">
                vAN-20170503
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170502/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170502"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170502">
                vAN-20170502
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170501/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170501"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170501">
                vAN-20170501
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170430/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170430"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170430">
                vAN-20170430
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170429/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170429"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170429">
                vAN-20170429
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170428/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170428"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170428">
                vAN-20170428
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170427/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170427"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170427">
                vAN-20170427
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170426/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170426"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170426">
                vAN-20170426
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170425/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170425"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170425">
                vAN-20170425
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170424/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170424"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170424">
                vAN-20170424
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170423/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170423"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170423">
                vAN-20170423
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170422/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170422"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170422">
                vAN-20170422
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170421/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170421"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170421">
                vAN-20170421
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170420/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170420"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170420">
                vAN-20170420
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170419/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170419"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170419">
                vAN-20170419
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170418/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170418"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170418">
                vAN-20170418
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170417/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170417"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170417">
                vAN-20170417
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170416/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170416"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170416">
                vAN-20170416
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170415/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170415"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170415">
                vAN-20170415
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170414/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170414"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170414">
                vAN-20170414
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170413/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170413"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170413">
                vAN-20170413
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170412/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170412"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170412">
                vAN-20170412
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170411/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170411"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170411">
                vAN-20170411
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170410/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170410"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170410">
                vAN-20170410
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170409/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170409"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170409">
                vAN-20170409
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170408/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170408"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170408">
                vAN-20170408
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170407/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170407"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170407">
                vAN-20170407
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170406/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170406"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170406">
                vAN-20170406
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170405/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170405"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170405">
                vAN-20170405
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170404/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170404"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170404">
                vAN-20170404
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170403/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170403"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170403">
                vAN-20170403
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170402/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170402"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170402">
                vAN-20170402
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170401/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170401"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170401">
                vAN-20170401
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170331/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170331"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170331">
                vAN-20170331
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170330/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170330"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170330">
                vAN-20170330
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170329/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170329"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170329">
                vAN-20170329
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170328/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170328"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170328">
                vAN-20170328
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170327/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170327"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170327">
                vAN-20170327
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170326/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170326"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170326">
                vAN-20170326
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170325/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170325"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170325">
                vAN-20170325
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170324/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170324"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170324">
                vAN-20170324
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170323/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170323"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170323">
                vAN-20170323
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170322/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170322"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170322">
                vAN-20170322
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170321/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170321"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170321">
                vAN-20170321
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170320/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170320"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170320">
                vAN-20170320
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170319/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170319"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170319">
                vAN-20170319
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170318/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170318"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170318">
                vAN-20170318
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170317/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170317"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170317">
                vAN-20170317
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170316/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170316"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170316">
                vAN-20170316
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170315/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170315"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170315">
                vAN-20170315
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170314/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170314"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170314">
                vAN-20170314
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170313/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170313"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170313">
                vAN-20170313
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170312/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170312"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170312">
                vAN-20170312
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170311/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170311"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170311">
                vAN-20170311
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/alisw/AliPhysics/tree/vAN-20170310/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C"
              data-name="vAN-20170310"
              data-skip-pjax="true"
              rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="vAN-20170310">
                vAN-20170310
              </span>
            </a>
        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div>

    </div>
  </div>
</div>

    <div class="BtnGroup float-right">
      <a href="/alisw/AliPhysics/find/master"
            class="js-pjax-capture-input btn btn-sm BtnGroup-item"
            data-pjax
            data-hotkey="t">
        Find file
      </a>
      <button aria-label="Copy file path to clipboard" class="js-zeroclipboard btn btn-sm BtnGroup-item tooltipped tooltipped-s" data-copied-hint="Copied!" type="button">Copy path</button>
    </div>
    <div class="breadcrumb js-zeroclipboard-target">
      <span class="repo-root js-repo-root"><span class="js-path-segment"><a href="/alisw/AliPhysics"><span>AliPhysics</span></a></span></span><span class="separator">/</span><span class="js-path-segment"><a href="/alisw/AliPhysics/tree/master/PWGLF"><span>PWGLF</span></a></span><span class="separator">/</span><span class="js-path-segment"><a href="/alisw/AliPhysics/tree/master/PWGLF/SPECTRA"><span>SPECTRA</span></a></span><span class="separator">/</span><span class="js-path-segment"><a href="/alisw/AliPhysics/tree/master/PWGLF/SPECTRA/ChargedHadrons"><span>ChargedHadrons</span></a></span><span class="separator">/</span><span class="js-path-segment"><a href="/alisw/AliPhysics/tree/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt"><span>dNdPt</span></a></span><span class="separator">/</span><span class="js-path-segment"><a href="/alisw/AliPhysics/tree/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros"><span>macros</span></a></span><span class="separator">/</span><strong class="final-path">AddTask_dNdPtpp_TPCITS.C</strong>
    </div>
  </div>


  
  <div class="commit-tease">
      <span class="float-right">
        <a class="commit-tease-sha" href="/alisw/AliPhysics/commit/fe5e4f52fca2fa9910672e843da4ba03711c3c91" data-pjax>
          fe5e4f5
        </a>
        <relative-time datetime="2017-06-06T12:33:46Z">Jun 6, 2017</relative-time>
      </span>
      <div>
        <img alt="@eperezle" class="avatar" height="20" src="https://avatars0.githubusercontent.com/u/26792896?v=3&amp;s=40" width="20" />
        <a href="/eperezle" class="user-mention" rel="contributor">eperezle</a>
          <a href="/alisw/AliPhysics/commit/fe5e4f52fca2fa9910672e843da4ba03711c3c91" class="message" data-pjax="true" title="IsPythia flag and BG rejection">IsPythia flag and BG rejection</a>
      </div>

    <div class="commit-tease-contributors">
      <button type="button" class="btn-link muted-link contributors-toggle" data-facebox="#blob_contributors_box">
        <strong>2</strong>
         contributors
      </button>
          <a class="avatar-link tooltipped tooltipped-s" aria-label="eperezle" href="/alisw/AliPhysics/commits/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C?author=eperezle"><img alt="@eperezle" class="avatar" height="20" src="https://avatars0.githubusercontent.com/u/26792896?v=3&amp;s=40" width="20" /> </a>
    <a class="avatar-link tooltipped tooltipped-s" aria-label="jgronefe" href="/alisw/AliPhysics/commits/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C?author=jgronefe"><img alt="@jgronefe" class="avatar" height="20" src="https://avatars2.githubusercontent.com/u/26566034?v=3&amp;s=40" width="20" /> </a>


    </div>

    <div id="blob_contributors_box" style="display:none">
      <h2 class="facebox-header" data-facebox-id="facebox-header">Users who have contributed to this file</h2>
      <ul class="facebox-user-list" data-facebox-id="facebox-description">
          <li class="facebox-user-list-item">
            <img alt="@eperezle" height="24" src="https://avatars2.githubusercontent.com/u/26792896?v=3&amp;s=48" width="24" />
            <a href="/eperezle">eperezle</a>
          </li>
          <li class="facebox-user-list-item">
            <img alt="@jgronefe" height="24" src="https://avatars0.githubusercontent.com/u/26566034?v=3&amp;s=48" width="24" />
            <a href="/jgronefe">jgronefe</a>
          </li>
      </ul>
    </div>
  </div>

  <div class="file">
    <div class="file-header">
  <div class="file-actions">

    <div class="BtnGroup">
      <a href="/alisw/AliPhysics/raw/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C" class="btn btn-sm BtnGroup-item" id="raw-url">Raw</a>
        <a href="/alisw/AliPhysics/blame/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C" class="btn btn-sm js-update-url-with-hash BtnGroup-item" data-hotkey="b">Blame</a>
      <a href="/alisw/AliPhysics/commits/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C" class="btn btn-sm BtnGroup-item" rel="nofollow">History</a>
    </div>


        <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="/alisw/AliPhysics/edit/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C" class="inline-form js-update-url-with-hash" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="gAIAxAcEZ1QZPeyC9+CloG+SwlOSAX/pc++jBr+acZl//jsB2TxzLMu1rBP5rIXbXi8fbLU9eGSzFyoqQ1EnEw==" /></div>
          <button class="btn-octicon tooltipped tooltipped-nw" type="submit"
            aria-label="Edit the file in your fork of this project" data-hotkey="e" data-disable-with>
            <svg aria-hidden="true" class="octicon octicon-pencil" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path fill-rule="evenodd" d="M0 12v3h3l8-8-3-3-8 8zm3 2H1v-2h1v1h1v1zm10.3-9.3L12 6 9 3l1.3-1.3a.996.996 0 0 1 1.41 0l1.59 1.59c.39.39.39 1.02 0 1.41z"/></svg>
          </button>
</form>        <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="/alisw/AliPhysics/delete/master/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/AddTask_dNdPtpp_TPCITS.C" class="inline-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="poHewUu4QaobEvZpgdbLgd8l6p4bxKTj5XCWntQKpbDb86yWsG3N3y5LTu4Rw8zYyRsBqyYuKy7gK1egB9cPkg==" /></div>
          <button class="btn-octicon btn-octicon-danger tooltipped tooltipped-nw" type="submit"
            aria-label="Delete the file in your fork of this project" data-disable-with>
            <svg aria-hidden="true" class="octicon octicon-trashcan" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M11 2H9c0-.55-.45-1-1-1H5c-.55 0-1 .45-1 1H2c-.55 0-1 .45-1 1v1c0 .55.45 1 1 1v9c0 .55.45 1 1 1h7c.55 0 1-.45 1-1V5c.55 0 1-.45 1-1V3c0-.55-.45-1-1-1zm-1 12H3V5h1v8h1V5h1v8h1V5h1v8h1V5h1v9zm1-10H2V3h9v1z"/></svg>
          </button>
</form>  </div>

  <div class="file-info">
      117 lines (90 sloc)
      <span class="file-info-divider"></span>
    5.68 KB
  </div>
</div>

    

  <div itemprop="text" class="blob-wrapper data type-c">
      <table class="highlight tab-size js-file-line-container" data-tab-size="8">
      <tr>
        <td id="L1" class="blob-num js-line-number" data-line-number="1"></td>
        <td id="LC1" class="blob-code blob-code-inner js-file-line"><span class="pl-k">void</span> <span class="pl-en">AddTask_dNdPtpp_TPCITS</span>(Int_t cutMode =<span class="pl-c1">222</span> , <span class="pl-k">char</span> *particleMode =<span class="pl-s"><span class="pl-pds">&quot;</span>default<span class="pl-pds">&quot;</span></span>, <span class="pl-k">char</span>* eventTrigger=<span class="pl-s"><span class="pl-pds">&quot;</span>kINT7<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L2" class="blob-num js-line-number" data-line-number="2"></td>
        <td id="LC2" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3" class="blob-num js-line-number" data-line-number="3"></td>
        <td id="LC3" class="blob-code blob-code-inner js-file-line">  TString <span class="pl-smi">stEventTrigger</span>(eventTrigger);</td>
      </tr>
      <tr>
        <td id="L4" class="blob-num js-line-number" data-line-number="4"></td>
        <td id="LC4" class="blob-code blob-code-inner js-file-line">  TString <span class="pl-smi">stParticleMode</span>(particleMode);</td>
      </tr>
      <tr>
        <td id="L5" class="blob-num js-line-number" data-line-number="5"></td>
        <td id="LC5" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L6" class="blob-num js-line-number" data-line-number="6"></td>
        <td id="LC6" class="blob-code blob-code-inner js-file-line">  AliAnalysisManager *mgr = <span class="pl-c1">AliAnalysisManager::GetAnalysisManager</span>();</td>
      </tr>
      <tr>
        <td id="L7" class="blob-num js-line-number" data-line-number="7"></td>
        <td id="LC7" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span> (!mgr) {<span class="pl-c1">Error</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>AddTask_dNdPtpp_TPCITS<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>No analysis manager found.<span class="pl-pds">&quot;</span></span>);<span class="pl-k">return</span> <span class="pl-c1">0</span>;}</td>
      </tr>
      <tr>
        <td id="L8" class="blob-num js-line-number" data-line-number="8"></td>
        <td id="LC8" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L9" class="blob-num js-line-number" data-line-number="9"></td>
        <td id="LC9" class="blob-code blob-code-inner js-file-line">  Bool_t hasMC=(<span class="pl-c1">AliAnalysisManager::GetAnalysisManager</span>()-&gt;<span class="pl-c1">GetMCtruthEventHandler</span>()!=0x0);</td>
      </tr>
      <tr>
        <td id="L10" class="blob-num js-line-number" data-line-number="10"></td>
        <td id="LC10" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L11" class="blob-num js-line-number" data-line-number="11"></td>
        <td id="LC11" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> Switch off all AliInfo (too much output!!!)</span></td>
      </tr>
      <tr>
        <td id="L12" class="blob-num js-line-number" data-line-number="12"></td>
        <td id="LC12" class="blob-code blob-code-inner js-file-line">  <span class="pl-c1">AliLog::SetGlobalLogLevel</span>(AliLog::<span class="pl-c1">kError</span>);</td>
      </tr>
      <tr>
        <td id="L13" class="blob-num js-line-number" data-line-number="13"></td>
        <td id="LC13" class="blob-code blob-code-inner js-file-line">  mgr-&gt;<span class="pl-c1">SetDebugLevel</span>(<span class="pl-c1">0</span>);</td>
      </tr>
      <tr>
        <td id="L14" class="blob-num js-line-number" data-line-number="14"></td>
        <td id="LC14" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L15" class="blob-num js-line-number" data-line-number="15"></td>
        <td id="LC15" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L16" class="blob-num js-line-number" data-line-number="16"></td>
        <td id="LC16" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span>/ Create event cuts</span></td>
      </tr>
      <tr>
        <td id="L17" class="blob-num js-line-number" data-line-number="17"></td>
        <td id="LC17" class="blob-code blob-code-inner js-file-line">  Float_t zvWindow = <span class="pl-c1">30</span>. ;</td>
      </tr>
      <tr>
        <td id="L18" class="blob-num js-line-number" data-line-number="18"></td>
        <td id="LC18" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L19" class="blob-num js-line-number" data-line-number="19"></td>
        <td id="LC19" class="blob-code blob-code-inner js-file-line">  AlidNdPtEventCuts *evtCuts = new <span class="pl-c1">AlidNdPtEventCuts</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>AlidNdPtEventCuts<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>Event cuts<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L20" class="blob-num js-line-number" data-line-number="20"></td>
        <td id="LC20" class="blob-code blob-code-inner js-file-line">  evtCuts-&gt;<span class="pl-c1">SetZvRange</span>(-zvWindow,zvWindow);</td>
      </tr>
      <tr>
        <td id="L21" class="blob-num js-line-number" data-line-number="21"></td>
        <td id="LC21" class="blob-code blob-code-inner js-file-line">  evtCuts-&gt;<span class="pl-c1">SetMeanXYZv</span>(<span class="pl-c1">0.0</span>,<span class="pl-c1">0.0</span>,<span class="pl-c1">0.0</span>);</td>
      </tr>
      <tr>
        <td id="L22" class="blob-num js-line-number" data-line-number="22"></td>
        <td id="LC22" class="blob-code blob-code-inner js-file-line">  evtCuts-&gt;<span class="pl-c1">SetSigmaMeanXYZv</span>(<span class="pl-c1">1.0</span>,<span class="pl-c1">1.0</span>,<span class="pl-c1">10.0</span>);</td>
      </tr>
      <tr>
        <td id="L23" class="blob-num js-line-number" data-line-number="23"></td>
        <td id="LC23" class="blob-code blob-code-inner js-file-line">  evtCuts-&gt;<span class="pl-c1">SetTriggerRequired</span>(<span class="pl-c1">kTRUE</span>);</td>
      </tr>
      <tr>
        <td id="L24" class="blob-num js-line-number" data-line-number="24"></td>
        <td id="LC24" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L25" class="blob-num js-line-number" data-line-number="25"></td>
        <td id="LC25" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> Create geom. acceptance cuts</span></td>
      </tr>
      <tr>
        <td id="L26" class="blob-num js-line-number" data-line-number="26"></td>
        <td id="LC26" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L27" class="blob-num js-line-number" data-line-number="27"></td>
        <td id="LC27" class="blob-code blob-code-inner js-file-line">  Float_t etaWindow = <span class="pl-c1">1</span>. ;</td>
      </tr>
      <tr>
        <td id="L28" class="blob-num js-line-number" data-line-number="28"></td>
        <td id="LC28" class="blob-code blob-code-inner js-file-line">  Float_t ptMin = <span class="pl-c1">0.10</span>;</td>
      </tr>
      <tr>
        <td id="L29" class="blob-num js-line-number" data-line-number="29"></td>
        <td id="LC29" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L30" class="blob-num js-line-number" data-line-number="30"></td>
        <td id="LC30" class="blob-code blob-code-inner js-file-line">  AlidNdPtAcceptanceCuts *accCuts = new <span class="pl-c1">AlidNdPtAcceptanceCuts</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>AlidNdPtAcceptanceCuts<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>Geom. acceptance cuts<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L31" class="blob-num js-line-number" data-line-number="31"></td>
        <td id="LC31" class="blob-code blob-code-inner js-file-line">  accCuts-&gt;<span class="pl-c1">SetEtaRange</span>(-etaWindow,etaWindow);</td>
      </tr>
      <tr>
        <td id="L32" class="blob-num js-line-number" data-line-number="32"></td>
        <td id="LC32" class="blob-code blob-code-inner js-file-line">  accCuts-&gt;<span class="pl-c1">SetPtRange</span>(ptMin,<span class="pl-c1">1</span>.<span class="pl-smi">e10</span>);</td>
      </tr>
      <tr>
        <td id="L33" class="blob-num js-line-number" data-line-number="33"></td>
        <td id="LC33" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L34" class="blob-num js-line-number" data-line-number="34"></td>
        <td id="LC34" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">gROOT</span>-&gt;<span class="pl-c1">LoadMacro</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>$ALICE_PHYSICS/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L35" class="blob-num js-line-number" data-line-number="35"></td>
        <td id="LC35" class="blob-code blob-code-inner js-file-line">  AliESDtrackCuts* esdTrackCuts = <span class="pl-c1">CreatedNdPtTrackCuts</span>(cutMode,hasMC);</td>
      </tr>
      <tr>
        <td id="L36" class="blob-num js-line-number" data-line-number="36"></td>
        <td id="LC36" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span> (!esdTrackCuts) { <span class="pl-c1">printf</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>ERROR: esdTrackCuts could not be created<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span>); <span class="pl-k">return</span>; }</td>
      </tr>
      <tr>
        <td id="L37" class="blob-num js-line-number" data-line-number="37"></td>
        <td id="LC37" class="blob-code blob-code-inner js-file-line">  esdTrackCuts-&gt;<span class="pl-c1">SetHistogramsOn</span>(<span class="pl-c1">kTRUE</span>);</td>
      </tr>
      <tr>
        <td id="L38" class="blob-num js-line-number" data-line-number="38"></td>
        <td id="LC38" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L39" class="blob-num js-line-number" data-line-number="39"></td>
        <td id="LC39" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> Create task</span></td>
      </tr>
      <tr>
        <td id="L40" class="blob-num js-line-number" data-line-number="40"></td>
        <td id="LC40" class="blob-code blob-code-inner js-file-line">  AlidNdPtTask *task = new <span class="pl-c1">AlidNdPtTask</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>AlidNdPtTask<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L41" class="blob-num js-line-number" data-line-number="41"></td>
        <td id="LC41" class="blob-code blob-code-inner js-file-line">  task-&gt;<span class="pl-c1">SetUseMCInfo</span>(hasMC);</td>
      </tr>
      <tr>
        <td id="L42" class="blob-num js-line-number" data-line-number="42"></td>
        <td id="LC42" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L43" class="blob-num js-line-number" data-line-number="43"></td>
        <td id="LC43" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> trigger selection: MB</span></td>
      </tr>
      <tr>
        <td id="L44" class="blob-num js-line-number" data-line-number="44"></td>
        <td id="LC44" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(stEventTrigger.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>kINT7<span class="pl-pds">&quot;</span></span>)) task-&gt;<span class="pl-c1">SelectCollisionCandidates</span>(AliVEvent::<span class="pl-c1">kINT7</span>);</td>
      </tr>
      <tr>
        <td id="L45" class="blob-num js-line-number" data-line-number="45"></td>
        <td id="LC45" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span>(stEventTrigger.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>kMB<span class="pl-pds">&quot;</span></span>)) task-&gt;<span class="pl-c1">SelectCollisionCandidates</span>(AliVEvent::<span class="pl-c1">kMB</span>);</td>
      </tr>
      <tr>
        <td id="L46" class="blob-num js-line-number" data-line-number="46"></td>
        <td id="LC46" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L47" class="blob-num js-line-number" data-line-number="47"></td>
        <td id="LC47" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L48" class="blob-num js-line-number" data-line-number="48"></td>
        <td id="LC48" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> Create cut analysis object</span></td>
      </tr>
      <tr>
        <td id="L49" class="blob-num js-line-number" data-line-number="49"></td>
        <td id="LC49" class="blob-code blob-code-inner js-file-line">  AlidNdPtAnalysis *fdNdPtAnalysispp = new <span class="pl-c1">AlidNdPtAnalysis</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>dNdPtAnalysis<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>dN/dPt Analysis<span class="pl-pds">&quot;</span></span>);</td>
      </tr>
      <tr>
        <td id="L50" class="blob-num js-line-number" data-line-number="50"></td>
        <td id="LC50" class="blob-code blob-code-inner js-file-line">  fdNdPtAnalysispp-&gt;<span class="pl-c1">SetEventCuts</span>(evtCuts);</td>
      </tr>
      <tr>
        <td id="L51" class="blob-num js-line-number" data-line-number="51"></td>
        <td id="LC51" class="blob-code blob-code-inner js-file-line">  fdNdPtAnalysispp-&gt;<span class="pl-c1">SetAcceptanceCuts</span>(accCuts);</td>
      </tr>
      <tr>
        <td id="L52" class="blob-num js-line-number" data-line-number="52"></td>
        <td id="LC52" class="blob-code blob-code-inner js-file-line">  fdNdPtAnalysispp-&gt;<span class="pl-c1">SetTrackCuts</span>(esdTrackCuts);</td>
      </tr>
      <tr>
        <td id="L53" class="blob-num js-line-number" data-line-number="53"></td>
        <td id="LC53" class="blob-code blob-code-inner js-file-line">  fdNdPtAnalysispp-&gt;<span class="pl-c1">SetAnalysisMode</span>(AlidNdPtHelper::<span class="pl-c1">kTPCITS</span>);</td>
      </tr>
      <tr>
        <td id="L54" class="blob-num js-line-number" data-line-number="54"></td>
        <td id="LC54" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(stEventTrigger.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>kINT7<span class="pl-pds">&quot;</span></span>)) fdNdPtAnalysispp-&gt;<span class="pl-c1">SetTriggerMask</span>(AliVEvent::<span class="pl-c1">kINT7</span>);</td>
      </tr>
      <tr>
        <td id="L55" class="blob-num js-line-number" data-line-number="55"></td>
        <td id="LC55" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span>(stEventTrigger.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>kMB<span class="pl-pds">&quot;</span></span>)) fdNdPtAnalysispp-&gt;<span class="pl-c1">SetTriggerMask</span>(AliVEvent::<span class="pl-c1">kMB</span>);</td>
      </tr>
      <tr>
        <td id="L56" class="blob-num js-line-number" data-line-number="56"></td>
        <td id="LC56" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span>fdNdPtAnalysispp-&gt;SetUseSPDClusterVsTrackletRejection(kTRUE);</span></td>
      </tr>
      <tr>
        <td id="L57" class="blob-num js-line-number" data-line-number="57"></td>
        <td id="LC57" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span>fdNdPtAnalysispp-&gt;SetRequireCompleteDAQ(kTRUE);</span></td>
      </tr>
      <tr>
        <td id="L58" class="blob-num js-line-number" data-line-number="58"></td>
        <td id="LC58" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span>fdNdPtAnalysispp-&gt;SetUseTOFBunchCrossing(kFALSE);</span></td>
      </tr>
      <tr>
        <td id="L59" class="blob-num js-line-number" data-line-number="59"></td>
        <td id="LC59" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(hasMC) </td>
      </tr>
      <tr>
        <td id="L60" class="blob-num js-line-number" data-line-number="60"></td>
        <td id="LC60" class="blob-code blob-code-inner js-file-line">  {</td>
      </tr>
      <tr>
        <td id="L61" class="blob-num js-line-number" data-line-number="61"></td>
        <td id="LC61" class="blob-code blob-code-inner js-file-line">    fdNdPtAnalysispp-&gt;<span class="pl-c1">SetUseMCInfo</span>(<span class="pl-c1">kTRUE</span>);</td>
      </tr>
      <tr>
        <td id="L62" class="blob-num js-line-number" data-line-number="62"></td>
        <td id="LC62" class="blob-code blob-code-inner js-file-line">    fdNdPtAnalysispp-&gt;<span class="pl-c1">SetHistogramsOn</span>(<span class="pl-c1">kTRUE</span>);</td>
      </tr>
      <tr>
        <td id="L63" class="blob-num js-line-number" data-line-number="63"></td>
        <td id="LC63" class="blob-code blob-code-inner js-file-line">    fdNdPtAnalysispp-&gt;<span class="pl-c1">SetIsPythia</span>(<span class="pl-c1">kTRUE</span>);</td>
      </tr>
      <tr>
        <td id="L64" class="blob-num js-line-number" data-line-number="64"></td>
        <td id="LC64" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">//</span>fdNdPtAnalysis-&gt;SetHistogramsOn(kFALSE);</span></td>
      </tr>
      <tr>
        <td id="L65" class="blob-num js-line-number" data-line-number="65"></td>
        <td id="LC65" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L66" class="blob-num js-line-number" data-line-number="66"></td>
        <td id="LC66" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L67" class="blob-num js-line-number" data-line-number="67"></td>
        <td id="LC67" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L68" class="blob-num js-line-number" data-line-number="68"></td>
        <td id="LC68" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> SetParticleMode</span></td>
      </tr>
      <tr>
        <td id="L69" class="blob-num js-line-number" data-line-number="69"></td>
        <td id="LC69" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>Pion<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCPion</span>);}</td>
      </tr>
      <tr>
        <td id="L70" class="blob-num js-line-number" data-line-number="70"></td>
        <td id="LC70" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>Proton<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCProton</span>);}</td>
      </tr>
      <tr>
        <td id="L71" class="blob-num js-line-number" data-line-number="71"></td>
        <td id="LC71" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>Kaon<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCKaon</span>);}</td>
      </tr>
      <tr>
        <td id="L72" class="blob-num js-line-number" data-line-number="72"></td>
        <td id="LC72" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>RemainingRest<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCRemainingRest</span>);}</td>
      </tr>
      <tr>
        <td id="L73" class="blob-num js-line-number" data-line-number="73"></td>
        <td id="LC73" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>Rest<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCRest</span>);}</td>
      </tr>
      <tr>
        <td id="L74" class="blob-num js-line-number" data-line-number="74"></td>
        <td id="LC74" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>Lambda<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCLambda</span>);}</td>
      </tr>
      <tr>
        <td id="L75" class="blob-num js-line-number" data-line-number="75"></td>
        <td id="LC75" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>SigmaPlus<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCSigmaPlus</span>);}</td>
      </tr>
      <tr>
        <td id="L76" class="blob-num js-line-number" data-line-number="76"></td>
        <td id="LC76" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>SigmaMinus<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCSigmaMinus</span>);}</td>
      </tr>
      <tr>
        <td id="L77" class="blob-num js-line-number" data-line-number="77"></td>
        <td id="LC77" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>XiMinus<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCXiMinus</span>);}</td>
      </tr>
      <tr>
        <td id="L78" class="blob-num js-line-number" data-line-number="78"></td>
        <td id="LC78" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>OmegaMinus<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCOmegaMinus</span>);}</td>
      </tr>
      <tr>
        <td id="L79" class="blob-num js-line-number" data-line-number="79"></td>
        <td id="LC79" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>Plus<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kPlus</span>);}</td>
      </tr>
      <tr>
        <td id="L80" class="blob-num js-line-number" data-line-number="80"></td>
        <td id="LC80" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>Minus<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMinus</span>);}</td>
      </tr>
      <tr>
        <td id="L81" class="blob-num js-line-number" data-line-number="81"></td>
        <td id="LC81" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>Electron<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCElectron</span>);}</td>
      </tr>
      <tr>
        <td id="L82" class="blob-num js-line-number" data-line-number="82"></td>
        <td id="LC82" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>Muon<span class="pl-pds">&quot;</span></span>)){fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kMCMuon</span>);}</td>
      </tr>
      <tr>
        <td id="L83" class="blob-num js-line-number" data-line-number="83"></td>
        <td id="LC83" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span> <span class="pl-k">if</span> (stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>InclWoSigma<span class="pl-pds">&quot;</span></span>)) {fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kInclWoSimga</span>);}</td>
      </tr>
      <tr>
        <td id="L84" class="blob-num js-line-number" data-line-number="84"></td>
        <td id="LC84" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">else</span>{fdNdPtAnalysispp-&gt;<span class="pl-c1">SetParticleMode</span>(AlidNdPtHelper::<span class="pl-c1">kAllPart</span>);}</td>
      </tr>
      <tr>
        <td id="L85" class="blob-num js-line-number" data-line-number="85"></td>
        <td id="LC85" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> Change binning</span></td>
      </tr>
      <tr>
        <td id="L86" class="blob-num js-line-number" data-line-number="86"></td>
        <td id="LC86" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">const</span> Int_t ptNbins = <span class="pl-c1">81</span>;</td>
      </tr>
      <tr>
        <td id="L87" class="blob-num js-line-number" data-line-number="87"></td>
        <td id="LC87" class="blob-code blob-code-inner js-file-line">  Double_t bins[<span class="pl-c1">82</span>] = {<span class="pl-c1">0.0</span>, <span class="pl-c1">0.05</span>, <span class="pl-c1">0.1</span>, <span class="pl-c1">0.15</span>, <span class="pl-c1">0.2</span>, <span class="pl-c1">0.25</span>, <span class="pl-c1">0.3</span>, <span class="pl-c1">0.35</span>, <span class="pl-c1">0.4</span>, <span class="pl-c1">0.45</span>, <span class="pl-c1">0.5</span>, <span class="pl-c1">0.55</span>, <span class="pl-c1">0.6</span>, <span class="pl-c1">0.65</span>, <span class="pl-c1">0.7</span>, <span class="pl-c1">0.75</span>, <span class="pl-c1">0.8</span>, <span class="pl-c1">0.85</span>, <span class="pl-c1">0.9</span>, <span class="pl-c1">0.95</span>, <span class="pl-c1">1.0</span>, <span class="pl-c1">1.1</span>, <span class="pl-c1">1.2</span>, <span class="pl-c1">1.3</span>, <span class="pl-c1">1.4</span>, <span class="pl-c1">1.5</span>, <span class="pl-c1">1.6</span>, <span class="pl-c1">1.7</span>, <span class="pl-c1">1.8</span>, <span class="pl-c1">1.9</span>, <span class="pl-c1">2.0</span>, <span class="pl-c1">2.2</span>, <span class="pl-c1">2.4</span>, <span class="pl-c1">2.6</span>, <span class="pl-c1">2.8</span>, <span class="pl-c1">3.0</span>, <span class="pl-c1">3.2</span>, <span class="pl-c1">3.4</span>, <span class="pl-c1">3.6</span>, <span class="pl-c1">3.8</span>, <span class="pl-c1">4.0</span>, <span class="pl-c1">4.5</span>, <span class="pl-c1">5.0</span>, <span class="pl-c1">5.5</span>, <span class="pl-c1">6.0</span>, <span class="pl-c1">6.5</span>, <span class="pl-c1">7.0</span>, <span class="pl-c1">8.0</span>, <span class="pl-c1">9.0</span>, <span class="pl-c1">10.0</span>, <span class="pl-c1">11.0</span>, <span class="pl-c1">12.0</span>, <span class="pl-c1">13.0</span>, <span class="pl-c1">14.0</span>, <span class="pl-c1">15.0</span>, <span class="pl-c1">16.0</span>, <span class="pl-c1">18.0</span>, <span class="pl-c1">20.0</span>, <span class="pl-c1">22.0</span>, <span class="pl-c1">24.0</span>, <span class="pl-c1">26.0</span>, <span class="pl-c1">28.0</span>, <span class="pl-c1">30.0</span>, <span class="pl-c1">32.0</span>, <span class="pl-c1">34.0</span>, <span class="pl-c1">36.0</span>, <span class="pl-c1">40.0</span>, <span class="pl-c1">45.0</span>, <span class="pl-c1">50.0</span>, <span class="pl-c1">60.0</span>, <span class="pl-c1">70.0</span>, <span class="pl-c1">80.0</span>, <span class="pl-c1">90.0</span>, <span class="pl-c1">100.0</span>, <span class="pl-c1">110.0</span>, <span class="pl-c1">120.0</span>, <span class="pl-c1">130.0</span>, <span class="pl-c1">140.0</span>, <span class="pl-c1">150.0</span>, <span class="pl-c1">160.0</span>, <span class="pl-c1">180.0</span>, <span class="pl-c1">200.0</span>};</td>
      </tr>
      <tr>
        <td id="L88" class="blob-num js-line-number" data-line-number="88"></td>
        <td id="LC88" class="blob-code blob-code-inner js-file-line">  Double_t* binsPt = new Double_t[<span class="pl-c1">82</span>];</td>
      </tr>
      <tr>
        <td id="L89" class="blob-num js-line-number" data-line-number="89"></td>
        <td id="LC89" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">for</span> (<span class="pl-k">int</span> i=<span class="pl-c1">0</span>; i&lt;<span class="pl-c1">82</span>; i++) {binsPt[i] = bins[i];}</td>
      </tr>
      <tr>
        <td id="L90" class="blob-num js-line-number" data-line-number="90"></td>
        <td id="LC90" class="blob-code blob-code-inner js-file-line">  fdNdPtAnalysispp-&gt;<span class="pl-c1">SetBinsPt</span>(ptNbins, binsPt);</td>
      </tr>
      <tr>
        <td id="L91" class="blob-num js-line-number" data-line-number="91"></td>
        <td id="LC91" class="blob-code blob-code-inner js-file-line">  fdNdPtAnalysispp-&gt;<span class="pl-c1">SetBinsPtCorr</span>(ptNbins, binsPt);</td>
      </tr>
      <tr>
        <td id="L92" class="blob-num js-line-number" data-line-number="92"></td>
        <td id="LC92" class="blob-code blob-code-inner js-file-line">  fdNdPtAnalysispp-&gt;<span class="pl-c1">SetUseMCInfo</span>(hasMC);</td>
      </tr>
      <tr>
        <td id="L93" class="blob-num js-line-number" data-line-number="93"></td>
        <td id="LC93" class="blob-code blob-code-inner js-file-line">  fdNdPtAnalysispp-&gt;<span class="pl-c1">SetHistogramsOn</span>(hasMC);</td>
      </tr>
      <tr>
        <td id="L94" class="blob-num js-line-number" data-line-number="94"></td>
        <td id="LC94" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L95" class="blob-num js-line-number" data-line-number="95"></td>
        <td id="LC95" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L96" class="blob-num js-line-number" data-line-number="96"></td>
        <td id="LC96" class="blob-code blob-code-inner js-file-line">  task-&gt;<span class="pl-c1">AddAnalysisObject</span>( fdNdPtAnalysispp );</td>
      </tr>
      <tr>
        <td id="L97" class="blob-num js-line-number" data-line-number="97"></td>
        <td id="LC97" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L98" class="blob-num js-line-number" data-line-number="98"></td>
        <td id="LC98" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> Add task</span></td>
      </tr>
      <tr>
        <td id="L99" class="blob-num js-line-number" data-line-number="99"></td>
        <td id="LC99" class="blob-code blob-code-inner js-file-line">  mgr-&gt;<span class="pl-c1">AddTask</span>(task);</td>
      </tr>
      <tr>
        <td id="L100" class="blob-num js-line-number" data-line-number="100"></td>
        <td id="LC100" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L101" class="blob-num js-line-number" data-line-number="101"></td>
        <td id="LC101" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L102" class="blob-num js-line-number" data-line-number="102"></td>
        <td id="LC102" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">//</span> Create containers for input</span></td>
      </tr>
      <tr>
        <td id="L103" class="blob-num js-line-number" data-line-number="103"></td>
        <td id="LC103" class="blob-code blob-code-inner js-file-line">  AliAnalysisDataContainer *cinput = mgr-&gt;<span class="pl-c1">GetCommonInputContainer</span>();</td>
      </tr>
      <tr>
        <td id="L104" class="blob-num js-line-number" data-line-number="104"></td>
        <td id="LC104" class="blob-code blob-code-inner js-file-line">  mgr-&gt;<span class="pl-c1">ConnectInput</span>(task, <span class="pl-c1">0</span>, cinput);</td>
      </tr>
      <tr>
        <td id="L105" class="blob-num js-line-number" data-line-number="105"></td>
        <td id="LC105" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L106" class="blob-num js-line-number" data-line-number="106"></td>
        <td id="LC106" class="blob-code blob-code-inner js-file-line">  TString stContainerName;</td>
      </tr>
      <tr>
        <td id="L107" class="blob-num js-line-number" data-line-number="107"></td>
        <td id="LC107" class="blob-code blob-code-inner js-file-line">  stContainerName = <span class="pl-c1">Form</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>dNdPt_pp_<span class="pl-c1">%d</span><span class="pl-pds">&quot;</span></span>,cutMode);</td>
      </tr>
      <tr>
        <td id="L108" class="blob-num js-line-number" data-line-number="108"></td>
        <td id="LC108" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(!stParticleMode.<span class="pl-c1">Contains</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>default<span class="pl-pds">&quot;</span></span>)) {</td>
      </tr>
      <tr>
        <td id="L109" class="blob-num js-line-number" data-line-number="109"></td>
        <td id="LC109" class="blob-code blob-code-inner js-file-line">    stContainerName = stContainerName + <span class="pl-s"><span class="pl-pds">&quot;</span>_<span class="pl-pds">&quot;</span></span> +stParticleMode;</td>
      </tr>
      <tr>
        <td id="L110" class="blob-num js-line-number" data-line-number="110"></td>
        <td id="LC110" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L111" class="blob-num js-line-number" data-line-number="111"></td>
        <td id="LC111" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L112" class="blob-num js-line-number" data-line-number="112"></td>
        <td id="LC112" class="blob-code blob-code-inner js-file-line">  AliAnalysisDataContainer *coutput = mgr-&gt;<span class="pl-c1">CreateContainer</span>(stContainerName,<span class="pl-c1">TList::Class</span>(),AliAnalysisManager::<span class="pl-c1">kOutputContainer</span>, mgr-&gt;<span class="pl-c1">GetCommonFileName</span>());</td>
      </tr>
      <tr>
        <td id="L113" class="blob-num js-line-number" data-line-number="113"></td>
        <td id="LC113" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L114" class="blob-num js-line-number" data-line-number="114"></td>
        <td id="LC114" class="blob-code blob-code-inner js-file-line">  mgr-&gt;<span class="pl-c1">ConnectOutput</span>(task, <span class="pl-c1">1</span>, coutput);</td>
      </tr>
      <tr>
        <td id="L115" class="blob-num js-line-number" data-line-number="115"></td>
        <td id="LC115" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L116" class="blob-num js-line-number" data-line-number="116"></td>
        <td id="LC116" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
</table>

  </div>

  </div>

  <button type="button" data-facebox="#jump-to-line" data-facebox-class="linejump" data-hotkey="l" class="d-none">Jump to Line</button>
  <div id="jump-to-line" style="display:none">
    <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="" class="js-jump-to-line-form" method="get"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /></div>
      <input class="form-control linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" aria-label="Jump to line" autofocus>
      <button type="submit" class="btn">Go</button>
</form>  </div>

  </div>
  <div class="modal-backdrop js-touch-events"></div>
</div>

    </div>
  </div>

  </div>

      
<div class="container site-footer-container">
  <div class="site-footer " role="contentinfo">
    <ul class="site-footer-links float-right">
        <li><a href="https://github.com/contact" data-ga-click="Footer, go to contact, text:contact">Contact GitHub</a></li>
      <li><a href="https://developer.github.com" data-ga-click="Footer, go to api, text:api">API</a></li>
      <li><a href="https://training.github.com" data-ga-click="Footer, go to training, text:training">Training</a></li>
      <li><a href="https://shop.github.com" data-ga-click="Footer, go to shop, text:shop">Shop</a></li>
        <li><a href="https://github.com/blog" data-ga-click="Footer, go to blog, text:blog">Blog</a></li>
        <li><a href="https://github.com/about" data-ga-click="Footer, go to about, text:about">About</a></li>

    </ul>

    <a href="https://github.com" aria-label="Homepage" class="site-footer-mark" title="GitHub">
      <svg aria-hidden="true" class="octicon octicon-mark-github" height="24" version="1.1" viewBox="0 0 16 16" width="24"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"/></svg>
</a>
    <ul class="site-footer-links">
      <li>&copy; 2017 <span title="0.25885s from github-fe136-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="https://github.com/site/terms" data-ga-click="Footer, go to terms, text:terms">Terms</a></li>
        <li><a href="https://github.com/site/privacy" data-ga-click="Footer, go to privacy, text:privacy">Privacy</a></li>
        <li><a href="https://github.com/security" data-ga-click="Footer, go to security, text:security">Security</a></li>
        <li><a href="https://status.github.com/" data-ga-click="Footer, go to status, text:status">Status</a></li>
        <li><a href="https://help.github.com" data-ga-click="Footer, go to help, text:help">Help</a></li>
    </ul>
  </div>
</div>



  <div id="ajax-error-message" class="ajax-error-message flash flash-error">
    <svg aria-hidden="true" class="octicon octicon-alert" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M8.865 1.52c-.18-.31-.51-.5-.87-.5s-.69.19-.87.5L.275 13.5c-.18.31-.18.69 0 1 .19.31.52.5.87.5h13.7c.36 0 .69-.19.86-.5.17-.31.18-.69.01-1L8.865 1.52zM8.995 13h-2v-2h2v2zm0-3h-2V6h2v4z"/></svg>
    <button type="button" class="flash-close js-flash-close js-ajax-error-dismiss" aria-label="Dismiss error">
      <svg aria-hidden="true" class="octicon octicon-x" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"/></svg>
    </button>
    You can't perform that action at this time.
  </div>


    
    <script crossorigin="anonymous" integrity="sha256-moJr+IVGtcuvm8fbBIStk4Dc4SZ+DnVTud0VEMrcYbY=" src="https://assets-cdn.github.com/assets/frameworks-9a826bf88546b5cbaf9bc7db0484ad9380dce1267e0e7553b9dd1510cadc61b6.js"></script>
    
    <script async="async" crossorigin="anonymous" integrity="sha256-KheWZq74zVU4NPSa3hRENOyR+lqo7VDHe7QK0+eukrs=" src="https://assets-cdn.github.com/assets/github-2a179666aef8cd553834f49ade144434ec91fa5aa8ed50c77bb40ad3e7ae92bb.js"></script>
    
    
    
    
  <div class="js-stale-session-flash stale-session-flash flash flash-warn flash-banner d-none">
    <svg aria-hidden="true" class="octicon octicon-alert" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path fill-rule="evenodd" d="M8.865 1.52c-.18-.31-.51-.5-.87-.5s-.69.19-.87.5L.275 13.5c-.18.31-.18.69 0 1 .19.31.52.5.87.5h13.7c.36 0 .69-.19.86-.5.17-.31.18-.69.01-1L8.865 1.52zM8.995 13h-2v-2h2v2zm0-3h-2V6h2v4z"/></svg>
    <span class="signed-in-tab-flash">You signed in with another tab or window. <a href="">Reload</a> to refresh your session.</span>
    <span class="signed-out-tab-flash">You signed out in another tab or window. <a href="">Reload</a> to refresh your session.</span>
  </div>
  <div class="facebox" id="facebox" style="display:none;">
  <div class="facebox-popup">
    <div class="facebox-content" role="dialog" aria-labelledby="facebox-header" aria-describedby="facebox-description">
    </div>
    <button type="button" class="facebox-close js-facebox-close" aria-label="Close modal">
      <svg aria-hidden="true" class="octicon octicon-x" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"/></svg>
    </button>
  </div>
</div>


  </body>
</html>

