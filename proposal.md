# Motivation

The COVID-19 pandemic has made clear, once again, we need robust packages, built
following best practices, in Epidemiology & Infectious Disease Modelling
([Jombart 2021](https://doi.org/10.1016/S1473-3099(20)30996-8)). But the
existence of good packages is just the first step in ensuring that they are
used. The packages also need to be *discoverable*. Lack of discoverability leads
to overlap and missed opportunities for collaboration. The existence of a task
view summarising the main existing tools (in R) is an important first milestone
to increase discoverability of Epidemiology & Infectious Disease Modelling R
packages.

# Scope

Packages included in this task view were identified through recommendations of
IDE experts as well as an automated CRAN search using `pkgsearch::pkg_search()`,
with the keywords: *epidemiology*, *epidemic*, *epi*, *outbreak* and
*transmission*. The list was manually curated for the final selection to satisfy
the conditions described in the previous paragraph.

Packages are deemed in scope if they provide tools, or data, explicitly targeted
at reporting, modelling, or forecasting infectious diseases.

# Packages

A current draft for the task view, with the proposed packages and their
description is available at: <https://github.com/Bisaloo/InfectiousEpi>

Future updates will focus on adding more text and connections between packages
to increase the added value of the task view relative to a simple list.

# Overlap

We expect relatively little overlap with other task views excepted those
focusing on general statistical and data analyses methods, which represent a
large portion of the daily work of epidemiologists. Such general task views are
for example `Bayesian`, `ReproducibleResearch`, `Spatial` or `TimeSeries`

# Maintainers

The primary maintainer is Thibaut Jombart (@thibautjombart). Current
co-maintainers are Hugo Gruson (@Bisaloo) and Matthieu Rolland (@mjrolland).
