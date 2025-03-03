---
name: Epidemiology
topic: Epidemiology
maintainer: Thibaut Jombart, Matthieu Rolland, Hugo Gruson
email: hugo.gruson+ctv@normalesup.org
version: 2025-03-03
source: https://github.com/cran-task-views/Epidemiology/
---

Contributors (in alphabetic order): Neale Batra, Solène Cadiou, Dylan Dijk,
Christopher Endres, Rich FitzJohn, Hugo Gruson, Andreas Handel, Michael Höhle,
Thibaut Jombart, Joseph Larmarange, Sebastian Lequime, Alex Spina, Tim Taylor,
Sean Wu, Achim Zeileis.

## Overview

R is increasingly becoming a standard in epidemiology, providing a wide array of
tools from study design to epidemiological data exploration, modeling,
forecasting, and simulation. This task view provides an overview of packages
specifically developed for epidemiology, including infectious disease
epidemiology (IDE) and environmental epidemiology. It does not include:

- generic tools which are used in these domains but not specifically developed
  for the epidemiological context,
- '*omics*' approaches and genome-wide association studies (GWAS), which can be
  used in epidemiology but form a largely separate domain.

Packages are grouped in the following categories:

1. **Data visualization:** tools dedicated to handling and
   visualization of epidemiological data, *e.g.* epidemic curves ('*epicurves*'),
   exploration of contact tracing networks, etc.
2. **Infectious disease modeling:** IDE-specific packages for the analysis of
   epidemic curves (including outbreak detection / surveillance), estimation of
   transmissibility, short-term forecasting, compartmental models (*e.g.*  SIR
   models), simulation of outbreaks, and reconstruction of transmission trees
3. **Environmental epidemiology:** tools dedicated to the study of
   environmental factors acting as determinants of diseases
4. **Helpers:** tools implementing miscellaneous tasks useful for practicing as
   well as teaching epidemiology, such as sample size calculation, fitting
   discretized Gamma distributions, or handling linelist data.
5. **Data packages:** these packages provide access to both empirical and
   simulated epidemic data; includes a specific section on COVID-19.

Additional links to non specific but highly useful packages (to create tables,
manipulate dates, etc.) are provided in the task view's footnotes.

## Inclusion criteria

Packages included in this task view were identified through recommendations of
expert epidemiologists as well as an automated CRAN search using
`pkgsearch::pkg_search()` with the keywords: *epidemiology*, *epidemic*, *epi*,
*outbreak*, and *transmission*. The list was manually curated for the final
selection to satisfy the conditions described in the previous paragraph.

Packages are deemed in scope if they provide tools, or data, explicitly targeted
at reporting, modeling, or forecasting infectious diseases.

**Your input is welcome!** Please suggest packages we may have missed by
filing an issue in the GitHub repository or by contacting the maintainer.

## Data cleaning and data management

- `r pkg("cleanepi")`: Provides functions to clean epidemiological
  data. It is designed to work with the `linelist` package and provides
  functions to check for missing data, validate dates, and ensure that variables
  are in the correct format.
- `r pkg("epiCleanr")`: A collection of data cleaning utilities for
  epidemiological data.
- `r pkg("epidm")`: Contains utilities and functions for the cleaning,
  processing and management of patient level public health data for surveillance
  and analysis held by the UK Health Security Agency, UKHSA.
- `r pkg("dataquieR")`: Data quality framework and tools to systematically check 
  health data for issues regarding data integrity, completeness, consistency or
  accuracy.

## Data visualization

This section includes packages providing specific tools for the visualization
and exploration of epidemiological data.

- `r pkg("epicontacts", priority = "core")`: Implements a dedicated class for
  contact data, composed of case line lists and contacts between cases. Also
  includes procedures for data handling, interactive graphics, and
  characterizing contact patterns (*e.g.* mixing patterns, serial
  intervals). [RECON](https://www.repidemicsconsortium.org/) package.
- `r pkg("EpiContactTrace")`: Routines for epidemiological contact tracing and
  visualization of networks of contacts.
- `r pkg("EpiCurve")`: Creates simple or stacked epidemic curves for hourly,
  daily, weekly or monthly outcome data.
- `r pkg("epiDisplay")`: Package for data exploration and result presentation.
- `r pkg("epiflows")`: Provides functions and classes designed to handle and
  visualize epidemiological flows of people between locations. Also contains a
  statistical method for predicting disease spread from flow data initially
  described in [Dorigatti et al.
  (2017)](https://doi.org/10.2807%2F1560-7917.ES.2017.22.28.30572).
  [RECON](https://www.repidemicsconsortium.org/) package.
- `r pkg("EpiReport")`: Drafting an epidemiological report in 'Microsoft Word'
  format for a given disease, similar to the Annual Epidemiological Reports
  published by the European Centre for Disease Prevention and Control.
- `r pkg("incidence")`: Functions and classes to compute, handle and visualize
  incidence from dated events for a defined time interval, using various date
  formats. Also provides wrappers for log-linear models of incidence and
  estimation of daily growth
  rate. [RECON](https://www.repidemicsconsortium.org/) package. This package
  is scheduled for deprecation and is replaced by `r pkg("incidence2")`.

## Infectious disease modeling

This section includes packages for specifically dedicated to IDE modeling. Note
that R offers a wealth of options for general-purpose time series modeling,
many of which are listed in the `r view("TimeSeries")` and `r view("Survival")`
task views.

### Epidemics surveillance

Packages below implement surveillance algorithms, but these approaches can be
usefully complemented by spatial analyses. We recommend looking at the `r view("Spatial")`
task view, which has a dedicated section on
`r view("Spatial", "disease mapping and areal data analysis")`.

- `r pkg("Epi", priority = "core")`: Functions for demographic and
  epidemiological analysis in the Lexis diagram, i.e. register and cohort
  follow-up data, in particular representation, manipulation and simulation of
  multistate data - the Lexis suite of functions, which includes interfaces to
  `r pkg("mstate")`, `r pkg("etm")` and `r pkg("cmprsk")` packages. Also
  contains functions for Age-Period-Cohort and Lee-Carter modeling, interval
  censored data, tabulation, plotting, as well as a number of epidemiological
  data sets.
- `r pkg("episensr")`: Basic sensitivity analysis of the observed relative
  risks adjusting for unmeasured confounding and misclassification of the
  exposure/outcome, or both. It follows the bias analysis methods and examples
  from the book by Lash T.L., Fox M.P., and Fink A.K. "Applying Quantitative
  Bias Analysis to Epidemiologic Data", ('Springer', 2009). This tool is also
  provided as an API via the `r pkg("apisensr")` package.
- `r pkg("mem")`: The Moving Epidemic Method, created by T Vega and JE Lozano
  ([2012](https://doi.org/10.1111/j.1750-2659.2012.00422.x),
  [2015](https://doi.org/10.1111/irv.12330)), allows the weekly assessment of
  the epidemic and intensity status to help in routine respiratory infections
  surveillance in health systems. Allows the comparison of different epidemic
  indicators, timing and shape with past epidemics and across different
  regions or countries with different surveillance systems. Also, it gives a
  measure of the performance of the method in terms of sensitivity and
  specificity of the alert week. This tool is also provided as a shiny app with
  the `r pkg("memapp")` package.
- `r pkg("riskCommunicator")`: Estimates flexible epidemiological effect
  measures including both differences and ratios using the parametric
  G-formula developed as an alternative to inverse probability weighting. It
  is useful for estimating the impact of interventions in the presence of
  treatment-confounder-feedback. G-computation was originally described by
  [Robbins (1986)](https://doi.org/10.1016/0270-0255(86)90088-6) and has been
  described in detail by [Ahern, Hubbard, and Galea
  (2009)](https://doi.org/10.1093/aje/kwp015); [Snowden, Rose, and Mortimer
  (2011)](https://doi.org/10.1093/aje/kwq472); and [Westreich et al.
  (2012)](https://doi.org/10.1002/sim.5316).
- `r pkg("RSurveillance")`: Provides a diverse set of functions useful for the
  design and analysis of disease surveillance activities.
- `r pkg("trendeval")`: Provides a coherent interface for evaluating models
  fit with the trending
  package. [RECON](https://www.repidemicsconsortium.org/) package.
- `r pkg("SpatialEpi")`: Methods and data for cluster detection and disease
  mapping.
- `r pkg("surveillance", priority = "core")`: Statistical methods for the
  modeling and monitoring of time series of counts, proportions and
  categorical data, as well as for the flexible covariate based regression
  modeling of spatio-temporal point processes of epidemic phenomena. The
  monitoring methods focus on aberration detection in count data time series
  from public health surveillance of communicable diseases. A recent overview
  of the available monitoring procedures is given by [Salmon et al.
  (2016)](https://doi.org/10.18637%2Fjss.v070.i10) and a recent overview of
  the implemented space-time modeling frameworks for epidemic phenomena is
  given by [Meyer et al. (2017)](https://doi.org/10.18637%2Fjss.v077.i11).
  Also contains back-projection methods to infer time series of exposure from
  disease onset and correction of observed time series for reporting delays
  (nowcasting).
- `r pkg("trending")`: Provides a coherent interface to multiple modeling
  tools for fitting trends along with a standardized approach for generating
  confidence and prediction
  intervals. [RECON](https://www.repidemicsconsortium.org/) package.
- `r pkg("coarseDataTools")`: Functions to analyze coarse data. Specifically, it
  contains functions to (1) fit parametric accelerated failure time models
  to interval-censored survival time data, and (2) estimate the case-fatality
  ratio in scenarios with under-reporting. This package's development was
  motivated by applications to infectious disease: in particular, problems with
  estimating the incubation period and the case fatality ratio of a given
  disease. Sample data files are included in the package. See
  [Reich et al. (2009)](https://doi.org/10.1002/sim.3659),
  [Reich et al. (2012)](https://doi.org/10.1111/j.1541-0420.2011.01709.x), and
  [Lessler et al. (2009)](https://doi.org/10.1016/S1473-3099(09)70069-6).
- `r pkg("cfr")`: Estimate the severity of a disease and ascertainment of cases,
  as discussed in
  [Nishiura et al. (2009)](https://doi.org/10.1371/journal.pone.0006852).
- `r pkg("EpiSignalDetection")`: Exploring epidemiological time series for
  signal detection via methods described in [Salmon et al.
  (2016)](https://doi.org/10.18637/jss.v070.i10). This package also provides a
  shiny interface and automated report generation.
- `r pkg("inctools")`: Tools for estimating incidence from biomarker data in
  cross-sectional surveys, and for calibrating tests for recent infection.
  Implements and extends the method of
  [Kassanjee et al. (2012)](https://doi.org/10.1097/EDE.0b013e3182576c07).
- `r pkg("epiCo")`: Tools for visualization and characterization of
  vector-borne disease outbreaks, with a specific focus on the Colombian
  context.
- `r pkg("sivirep")`: Automated reports for the analysis of surveillance data
  from the SIVIGILA (Sistema Nacional de Vigilancia en Salud Pública) in
  Colombia. In the process, many functions for data cleaning, manipulation, and
  visualization are used under the hood. These functions are also exported and
  can be used outside of the automated reports.
- `r pkg("ATQ")`: Framework for early detection of epidemics using school
  absenteeism data, based on 
  [Vanderkruk et al. (2023)](https://doi.oirg/10.1186/s12889-023-15747-z).
- `r pkg("EVI")`: Calculation of the Epidemic Volatility Index: a measure of
  volatility of newly reported cases, which can serve as an early warning 
  system when volatility exceeds a given threshold. The method is described
  in 
  [Koustolas et al. (2021)](https://www.nature.com/articles/s41598-021-02622-3).
- `r pkg("aedseo")`: Automated and Early Detection of Seasonal Epidemic Onset
  through estimation of the exponential growth rate and computation of the sum 
  of cases exceeding a certain threshold over the past k units of time.

#### Individual-level data

- `r pkg("modelSSE")`: Comprehensive analytical tools are provided to
  characterize infectious disease superspreading from contact tracing
  surveillance data. The underlying theoretical frameworks of this toolkit
  include branching process with transmission heterogeneity ([Lloyd-Smith et al.
  (2005)](https://doi.org/10.1038/nature04153), case cluster size distribution
  ([Nishiura et al. (2012)](https://doi.org/10.1016/j.jtbi.2011.10.039),
  [Blumberg et al. (2014)](https://doi.org/10.1371/journal.ppat.1004452), and
  [Kucharski and
  Althaus (2015)](https://doi.org/10.2807/1560-7917.ES2015.20.25.21167), and
  decomposition of reproduction number ([Zhao et al.
  (2022)](https://doi.org/10.1371/journal.pcbi.1010281).
- `r pkg("nosoi")`: The aim of `r pkg("nosoi")` (pronounced no.si) is to
  provide a flexible agent-based stochastic transmission chain/epidemic
  simulator ([Lequime et al. 2020](https://doi.org/2020.03.03.973107)). The
  package can take into account the influence of multiple variables on the
  transmission process (e.g.  dual-host systems such as arboviruses,
  within-host viral dynamics, transportation, population structure), alone or
  taken together, to create complex but relatively intuitive epidemiological
  simulations.
- `r pkg("epiworldR")` and its Shiny interface `r pkg("epiworldRShiny")`: A
  framework for agent-based modelling with a C++ backend for high performance.

#### Digital Epidemiology

- `r pkg("argo")`: Augmented Regression with General Online data (ARGO) for
  accurate estimation of influenza epidemics in United States on both national
  level and regional level. It replicates the method introduced in paper
  [Yang, S., Santillana, M. and Kou, S.C.
  (2015)](https://doi.org/10.1073/pnas.1515373112) and [Ning, S., Yang, S. and
  Kou, S.C. (2019)](https://doi.org/10.1038/s41598-019-41559-6).
- `r pkg("epitweetr")`: Early Detection of Public Health Threats from 'Twitter'
  Data. This package allows you to automatically monitor trends of tweets by
  time, place and topic aiming at detecting public health threats early through
  the detection of signals (e.g. an unusual increase in the number of tweets).
  It was designed to focus on infectious diseases, and it can be extended to all
  hazards or other fields of study by modifying the topics and keywords. More
  information is available in the ['epitweetr' peer-review
  publication](https://doi.org/10.2807/1560-7917.ES.2022.27.39.2200177).

### Estimation of transmissibility and forecasting

- `r pkg("earlyR")`: Implements a simple, likelihood-based estimation of the
  reproduction number (R0) using a Poisson branching process. This model
  requires knowledge of the serial interval distribution, and dates of symptom
  onsets. It is a simplified version of the model introduced by [Cori et
  al. (2013)](https://doi.org/10.1093%2Faje%2Fkwt133).
- `r pkg("endtoend")`: Computes the expectation of the number of transmissions
  and receptions considering an End-to-End transport model with limited number
  of retransmissions per packet. It provides theoretical results and also
  estimated values based on Monte Carlo simulations. It is also possible to
  consider random data and ACK probabilities.
- `r pkg("EpiEstim", priority = "core")`: Provides tools for estimating
  time-varying transmissibility using the instantaneous reproduction number
  (Rt) introduced in [Cori et al. (2013)](https://doi.org/10.1093/aje/kwt133).
- `r pkg("epinet")`: A collection of epidemic/network-related tools. Simulates
  transmission of diseases through contact networks. Performs Bayesian
  inference on network and epidemic parameters, given epidemic data.
- `r pkg("EpiNow2", priority = "core")`: Provides tools for estimating the
  time-varying reproduction number, rate of spread, and doubling time of
  epidemics while accounting for various delays using the approach introducted
  in [Abbott et al. (2020)](https://doi.org/10.12688/wellcomeopenres.16006.1),
  and [Gostic et al. (2020)](https://doi.org/10.1371/journal.pcbi.1008409).
- `r pkg("nbTransmission")`: Estimates the relative transmission probabilities
  between cases using naive Bayes as introduced in [Leavitt et
  al. (2020)](https://doi.org/10.1093%2Fije%2Fdyaa031). Includes various
  functions to estimate transmission parameters such as the generation/serial
  interval and reproductive number as well as finding the contribution of
  covariates to transmission probabilities and visualizing results.
- `r pkg("R0", priority = "core")`: Estimation of reproduction numbers for
  disease outbreak, based on incidence data including the basic reproduction
  number (R0) and the instantaneous reproduction number (R(t)), alongside
  corresponding 95% Confidence Interval. Also includes routines for plotting
  outputs and for performing sensitivity analyses.
- `r pkg("tsiR")`: The TSIR modeling framework allows users to fit the time
  series SIR model to cumulative case data, which uses a regression equation
  to estimate transmission parameters based on differences in cumulative cases
  between reporting periods. The package supports inference on TSIR parameters
  using GLMs and profile likelihood techniques, as well as forward simulation
  based on a fitted model, as described in [Becker and Grenfell
  (2017)](https://doi.org/10.1371/journal.pone.0185528).
- `r pkg("Bernadette")`: Implements the Bayesian evidence synthesis approach
  described in [Bouranis et al (2022)](https://arxiv.org/abs/2211.15229) to
  modeling the age-specific transmission dynamics of COVID-19 based on daily
  mortality counts. The functionality of `r pkg("Bernadette")` can be used to
  reconstruct the epidemic drivers from publicly available data, to estimate key
  epidemiological quantities like the time-varying rate of disease transmission,
  the latent counts of infections and the reproduction number for a given
  population over time, and to perform model comparison using information
  criteria.
- `r pkg("epicasting")`: Method and tool for generating time series forecasts
  using an ensemble wavelet-based auto-regressive neural network architecture.
  This method provides additional support of exogenous variables.
- `r pkg("epigrowthfit")`: Maximum likelihood estimation of nonlinear mixed
  effects models of epidemic growth. Provides several auxiliary functions,
  including one for computing basic reproduction numbers from fitted values of
  the initial exponential growth rate.
- `r pkg("BayesianFitForecast")`: Bayesian parameter estimation and forecasting
  in epidemiological models, as described in 
  [Karami et al. (2024)](https://doi.org/10.48550/arXiv.2411.05371) and 
  [Grinsztajn et al. (2021)](https://doi.org/10.1002/sim.9164).
- `r pkg("linelistBayes")`: Bayesian estimation of the reproduction number and 
  related metrics, on individual-level (linelist) or aggregated data, using the
  methods described in 
  [Li and White (2021)](https://doi.org/10.1371%2Fjournal.pcbi.1009210).
- `r pkg("rplanes")`: Prepare data and analyze plausibility of both forecasted
  and reported epidemiological signals.
- `r pkg("EpiInvert")`: Incidence curve decomposition (e.g., effect of 
  seasonality) and parameter estimation (e.g., Rt), as described in
  [Alvarez et al. (2021)](https://doi.org/10.1073/pnas.2105112118) and 
  [Alvarez et al. (2022)](https://doi.org/10.3390/biology11040540)
  and short-term forecasts.

### Epidemic simulation models

- `r pkg("EpiILM")`: Provides tools for simulating from discrete-time
  individual level models for infectious disease data analysis. This epidemic
  model class contains spatial and contact-network based models with two
  disease types: Susceptible-Infectious (SI) and
  Susceptible-Infectious-Removed (SIR).
- `r pkg("EpiILMCT")`: Provides tools for simulating from continuous-time
  individual level models of disease transmission, and carrying out infectious
  disease data analyses with the same models. The epidemic models considered
  are distance-based and/or contact network-based models within
  Susceptible-Infectious-Removed (SIR) or
  Susceptible-Infectious-Notified-Removed (SINR) compartmental frameworks. An
  overview of the implemented continuous-time individual level models for
  epidemics is given by [Almutiry and Deardon
  (2019)](https://doi.org/10.1515/ijb-2017-0092).
- `r pkg("EpiModel", priority = "core")`: Tools for simulating mathematical
  models of infectious disease dynamics. Epidemic model classes include
  deterministic compartmental models, stochastic individual-contact models,
  and stochastic network models.  Network models use the robust statistical
  methods of exponential-family random graph models (ERGMs) from the Statnet
  suite of software packages in R. Standard templates for epidemic modeling
  include SI, SIR, and SIS disease types. EpiModel features an API for
  extending these templates to address novel scientific research aims. Full
  methods for EpiModel are detailed in [Jenness et
  al. (2018)](https://doi.org/10.18637/jss.v084.i08).
- `r pkg("odin")`: Provides a generic, fast and computer-efficient platform
  for implementing any deterministic or stochastic compartmental models
  (e.g. SIR, SEIR, SIRS, ...), and can include age stratification or
  spatialization. It uses a domain specific language (DSL) to specify systems
  of ordinary differential equations (ODE) and integrate them. The DSL uses
  R's syntax, but compiles to C in order to efficiently solve the system,
  using interfaces to the packages `r pkg("deSolve")` and `r pkg("dde")`.
- `r pkg("pomp")` Provides a large set of forward simulation algorithms and
  MLE or Bayesian inference techniques to work with state-space models. Models
  may be specified as either deterministic or stochastic and typically follow
  a compartmental model structure. Time may be either discrete or continuous,
  depending on the simulation algorithm chosen. Additionally models may be
  programmed in C and compiled on-the-fly into a format suitable for use in
  the package to speed up simulation and inference. The R package and some of
  its algorithms are described in [King, Nguyen, and Ionides
  (2016)](https://doi.org/10.18637/jss.v069.i12).
- `r pkg("popEpi")`: Enables computation of epidemiological statistics,
  including those where counts or mortality rates of the reference population
  are used. Currently supported: excess hazard models, rates, mean survival
  times, relative survival, and standardized incidence and mortality ratios
  (SIRs/SMRs), all of which can be easily adjusted for by covariates such as
  age. Fast splitting and aggregation of 'Lexis' objects (from package `r
  pkg("Epi", priority = "core")` and other computations achieved using `r
  pkg("data.table")`.
- `r pkg("SimInf")`: Provides an efficient and very flexible framework to
  conduct data-driven epidemiological modeling in realistic large scale
  disease spread simulations. The framework integrates infection dynamics in
  subpopulations as continuous-time Markov chains using the Gillespie
  stochastic simulation algorithm and incorporates available data such as
  births, deaths and movements as scheduled events at predefined time-points.
  Using C code for the numerical solvers and 'OpenMP' (if available) to divide
  work over multiple processors ensures high performance when simulating a
  sample outcome. The package contains template models and can be extended
  with user-defined models. For more details see the paper by [Widgren, Bauer,
  Eriksson and Engblom (2019)](https://doi.org/10.18637%2Fjss.v091.i12).
- `r pkg("finalsize")`: Calculate the final size of a
  susceptible-infectious-recovered epidemic in a population with demographic
  variation in contact patterns and susceptibility to disease, as discussed in
  [Miller (2012)](https://doi.org/10.1007/s11538-012-9749-6).
- `r pkg("shinySIR")`: A Shiny graphical interface to interactively explore
  simple SIR models. Users can also provide their own ODEs.
- `r pkg("EpiDynamics")`: A collection of SEIR, SIR, and SIS models with
  several variations (demography, age groups, risk groups, etc.)
- `r pkg("epichains")`: A re-imagining of the bpmodels package to
  simulate branching process involved in epidemic transmission chains. This
  can also be used for estimation of epidemiological parameters by plugging the
  simulation output into a maximum likelihood or bayesian estimation procedure.
- `r pkg("genSEIR")`: Generalized Susceptible-Exposed-Infected-Recovered (SEIR)
  modeling, as described in 
  [Peng et al. (2020)](https://doi.org/10.1101/2020.02.16.20023465).

### Transmission tree reconstruction

- `r pkg("adegenet")`: primarily a population genetics package,
  `r pkg("adegenet")` implements seqtrack ([Jombart et al.
  2011](https://doi.org/10.1038%2Fhdy.2010.78)), a maximum-parsimony approach
  for reconstructing transmission trees using the Edmonds/Chu-Liu algorithm
- `r pkg("o2geosocial")` (`r pkg("outbreaker2")` module): Bayesian
  reconstruction of who infected whom during past outbreaks using
  routinely-collected surveillance data. Inference of transmission trees using
  genotype, age specific social contacts, distance between cases and onset
  dates of the reported cases ([Robert A, Kucharski AJ, Gastanaduy PA, Paul P,
  Funk S. 2020](https://doi.org/10.1098%2Frsif.2020.0084)).
- `r github("xavierdidelot/o2mod.transphylo")` is a module of
  `r pkg("outbreaker2")` which uses the `r pkg("TransPhylo")` model of
  within-host evolution.
- `r pkg("outbreaker2", priority = "core")`: a modular platform for Bayesian
  reconstruction of disease outbreaks using epidemiological and genetic
  information as introduced in [Jombart T, Cori A, Didelot X, Cauchemez S,
  Fraser C and Ferguson N.
  2014](https://doi.org/10.1371/journal.pcbi.1003457), [Campbell F, Cori A,
  Ferguson N, Jombart T (2019)](https://doi.org/10.1371/journal.pcbi.1006930).
- `r pkg("TransPhylo")`: Inference of transmission tree from a dated
  phylogeny. Includes methods to simulate and analyze outbreaks. The
  methodology is described in [Didelot et al.
  (2014)](https://doi.org/10.1093%2Fmolbev%2Fmsu121), [Didelot et al.
  (2017)](https://doi.org/10.1093%2Fmolbev%2Fmsw275).

### Vaccination and health outcomes

- `r pkg("vaccineff")`: Tools for estimating vaccine effectiveness and related
  metrics.
- `r pkg("excessmort")`: Implementation of method for estimating excess
  mortality and other health related outcomes from weekly or daily count data
  described in [Acosta and Irizarry
  (2021)](https://doi.org/10.1097/ede.0000000000001445).

### Plant epidemiology

- `r pkg("epiphy")`: A toolbox to analyze plant disease epidemics. It provides
  a common framework for plant disease intensity data recorded over time and/or
  space.
- `r pkg("epifitter")`: Analysis and visualization of plant disease progress
  curve data.
- `r pkg("cercospoRa")`: Epidemiological Model for Cercospora Leaf Spot of Sugar
  Beet.

## Environmental epidemiology

Environmental epidemiology is dedicated to the study of physical, chemical, and
biologic agents in the environment acting as determinants of disease. The aims
of environmental epidemiology are to infer causality, to identify environmental
causes of disease, such as from air and water pollutants, dietary contaminants,
built environments, and others.

R packages dedicated to environmental epidemiology include tools dealing with
limits of detection of pollutants (left-censoring issues), and various modeling
approaches to account for multiple correlations between exposures and infer
causality.

- `r pkg("NADA", priority = "core")`: Nondetects and Data Analysis for
  Environmental Data, package containing all the functions derived from the
  methods in [Helsel
  (2011)](https://www.wiley.com/en-us/Statistics+for+Censored+Environmental+Data+Using+Minitab+and+R%2C+2nd+Edition-p-9780470479889).
- `r pkg("EnvStats", priority = "core")`: Package for Environmental
  Statistics, Including US EPA Guidance, graphical and statistical analyses of
  environmental data, with focus on analyzing chemical concentrations and
  physical parameters, usually in the context of mandated environmental
  monitoring. Major environmental statistical methods found in the literature
  and regulatory guidance documents, with extensive help that explains what
  these methods do, how to use them, and where to find them in the
  literature. Numerous built-in data sets from regulatory guidance documents
  and environmental statistics literature [(Millard
  2013)](https://link.springer.com/book/10.1007/978-1-4614-8456-1).
- `r pkg("bkmr")`: Implements Bayesian Kernel Machine Regression, a statistical
  approach for estimating the joint health effects of multiple concurrent
  exposures, as described in [Bobb *et al.*
  (2015)](https://doi.org/10.1093/biostatistics/kxu058)
- `r pkg("mediation", priority = "core")`: Implements parametric and non
  parametric mediation analysis as discussed in [Imai *et al.*
  (2010)](https://doi.org/10.1214/10-sts321).
- `r pkg("mma")`: Implements multiple mediation analysis as described in Yu
  *et al.*  [(2017)](https://doi.org/10.1016%2Fj.sste.2017.02.001).
- `r pkg("HIMA")`: Allows to estimate and test high-dimensional mediation
  effects based on advanced mediator screening and penalized regression
  techniques ([Zhang *et al.*
  2021](https://doi.org/10.1093%2Fbioinformatics%2Fbtab564)).
- `r pkg("eesim")`: Functions to create simulated time series of environmental
  exposures (e.g., temperature, air pollution) and health outcomes for use in
  power analysis and simulation studies in environmental epidemiology.

## Helpers

This section includes packages providing tools to facilitate epidemiological
analysis as well as for training (e.g. computing sample size, contingency
tables, etc).

- `r pkg("incidence2", priority = "core")`: Provides functions and classes to
  compute, handle and visualize incidence from dated events. Improves the
  original `r pkg("incidence")` package in many ways: full flexibility in time
  intervals used, allows multiple stratifications, and is fully compatible with
  [dplyr](https://cran.r-project.org/packages=dplyr) and other tidyverse tools.
  [RECON](https://www.repidemicsconsortium.org/) package.
- `r pkg("DSAIDE")`: Exploration of simulation models (apps) of various
  infectious disease transmission dynamics scenarios. The purpose of the
  package is to help individuals learn about infectious disease epidemiology
  from a dynamical systems perspective. All apps include explanations of the
  underlying models and instructions on what to do with the models.
- `r pkg("epibasix")`: Contains elementary tools for analyzing common
  epidemiological problems, ranging from sample size estimation, through 2x2
  contingency table analysis and basic measures of agreement (kappa,
  sensitivity/specificity). Appropriate print and summary statements are also
  written to facilitate interpretation wherever possible. The target audience
  includes advanced undergraduate and graduate students in epidemiology or
  biostatistics courses, and clinical researchers.
- `r pkg("epiR", priority = "core")`: Tools for the analysis of
  epidemiological and surveillance data. Contains functions for directly and
  indirectly adjusting measures of disease frequency, quantifying measures of
  association on the basis of single or multiple strata of count data
  presented in a contingency table, computation of confidence intervals around
  incidence risk and incidence rate estimates and sample size calculations for
  cross-sectional, case-control and cohort studies. Surveillance tools include
  functions to calculate an appropriate sample size for 1- and 2-stage
  representative freedom surveys, functions to estimate surveillance system
  sensitivity and functions to support scenario tree modeling analyses.
- `r pkg("epitools", priority = "core")`: Tools for training and practicing
  epidemiologists including methods for two-way and multi-way contingency
  tables.
- `r pkg("tidyrates")`: A wrapper around `r pkg("epitools", priority = "core")`
  to allow tidy computation of age-adjusted epidemiological rates.
- `r pkg("epitrix")`: A collection of small functions useful for epidemics
  analysis and infectious disease modeling. This includes computation of
  basic reproduction numbers (R0) from daily growth rates, generation of
  hashed labels to anonymize data, and fitting discretized Gamma
  distributions.
- `r pkg("linelist")`: Implements the `linelist` class for storing case line
  list data, which extends `data.frame` and `tibble` by adding the ability to
  tag key epidemiological variables, validate them, and providing safeguards
  against accidental deletion or alteration of these data to help make data
  pipelines more straightforward and robust.
- `r pkg("powerSurvEpi")`: Functions to calculate power and sample size for
  testing main effect or interaction effect in the survival analysis of
  epidemiological studies (non-randomized studies), taking into account the
  correlation between the covariate of the interest and other covariates. Some
  calculations also take into account the competing risks and stratified
  analysis. This package also includes a set of functions to calculate power
  and sample size for testing main effects in the survival analysis of
  randomized clinical trials.
- `r pkg("AMR")`: Functions to simplify and standardise antimicrobial
  resistance (AMR) data analysis and to work with microbial and antimicrobial
  properties by using evidence-based methods and reliable reference data such
  as LPSN ([Parte *et al.* 2020](https://doi.org/10.1099/ijsem.0.004332)).
- `r pkg("diyar")`: Links records of individuals across multiple datasets.
- `r pkg("epikit")`: Contains tools for formatting inline code, renaming
  redundant columns, aggregating age categories, adding survey weights, finding
  the earliest date of an event, plotting z-curves, generating population counts
  and calculating proportions with confidence intervals.
- `r pkg("epitab")`: Builds contingency tables that cross-tabulate multiple
  categorical variables and also calculates various summary measures. Exports to
  a variety of formats is supported, including: 'HTML', 'LaTeX', and 'Excel'.
- `r pkg("EpiStats")`: Measures of association and impact for case control
  studies and cohort studies.
- `r pkg("EpiForsk")`: Collection of functions used by the Department of
  Epidemiological Research at Statens Serum Institut, Denmark.

## Data

Here are packages providing different epidemiologic datasets, either simulated
or real, useful for research purposes or field applications with a specific
COVID-19 section.

### Social contact data

- `r pkg("contactdata")`: Data package for the supplementary data in [Prem et
  al. (2017)](https://doi.org/10.1371%2Fjournal.pcbi.1005697). Provides easy
  access to contact data for 152 countries, for use in epidemiological,
  demographic or social sciences research.
- `r pkg("socialmixr")`: Provides methods for sampling contact matrices from
  diary data for use in infectious disease modeling, as discussed in [Mossong
  et al. (2008)](https://doi.org/10.1371%2Fjournal.pmed.0050074).

### Epidemic outbreak data

- `r pkg("outbreaks", priority = "core")`: Empirical or simulated disease
  outbreak data, provided either as RData or as text files.
- `r pkg("cholera")`: Amends errors, augments data and aids analysis of
  John Snow's map of the 1854 London cholera outbreak.
- `r pkg("malariaAtlas")`: A suite of tools to allow you to download all
  publicly available parasite rate survey points, mosquito occurrence points and
  raster surfaces from the ['Malaria Atlas Project'](https://malariaatlas.org/)
  servers as well as utility functions for plotting the downloaded data.
- `r pkg("colmozzie")`: Weekly notified dengue cases and climate variables in
  Colombo district Sri Lanka from 2008/ week-52 to 2014/ week-21.
- `r pkg("denguedatahub")`: Centralized access to dengue data worldwide from
  various sources.
- `r pkg("epidatr")`: API package to access real-time epidemiological
  surveillance data for influenza, 'COVID-19', and other diseases in the USA.

### COVID-19

- `r pkg("bets.covid19")`: Implements likelihood inference for early epidemic
  analysis. BETS is short for the four key epidemiological events being
  modeled: Begin of exposure, End of exposure, time of Transmission, and time
  of Symptom onset. The package contains a dataset of the trajectory of
  confirmed cases during the coronavirus disease (COVID-19) early outbreak.
  More detail of the statistical methods can be found in [Zhao et al.
  (2020)](https://arxiv.org/abs/2004.07743).
- `r pkg("corona")`: Manipulate and view coronavirus data and other societally
  relevant data at a basic level.
- `r pkg("coronavirus")`: Provides a daily summary of the Coronavirus
  (COVID-19) cases by state/province. Data source: [Johns Hopkins University
  Center for Systems Science and Engineering (JHU CCSE)
  Coronavirus](https://systems.jhu.edu/research/public-health/ncov/).
- `r pkg("COVID19")`: Download COVID-19 data across governmental sources at
  national, regional, and city level, as described in [Guidotti and Ardia
  (2020)](https://doi.org/10.21105/joss.02376). Includes the time series of
  vaccines, tests, cases, deaths, recovered, hospitalizations, intensive
  therapy, and policy measures by '[Oxford COVID-19 Government Response
  Tracker](https://www.bsg.ox.ac.uk/research/research-projects/coronavirus-government-response-tracker).
  Provides a seamless integration with '[World Bank Open
  Data](https://data.worldbank.org/)', '[Google Mobility
  Reports](https://www.google.com/covid19/mobility/)', 'Apple Mobility
  Reports](<https://covid19.apple.com/mobility>)'.
- `r pkg("covid19.analytics")`: Load and analyze updated time series worldwide
  data of reported cases for the Novel CoronaVirus Disease (CoViD-19) from the
  [Johns Hopkins University Center for Systems Science and Engineering (JHU
  CSSE) data repository](https://github.com/CSSEGISandData/COVID-19). The
  datasets are available in two main modalities, as a time series sequences
  and aggregated for the last day with greater spatial resolution. Several
  analysis, visualization and modeling functions are available in the package
  that will allow the user to compute and visualize total number of cases,
  total number of changes and growth rate globally or for an specific
  geographical location, while at the same time generating models using these
  trends; generate interactive visualizations and generate
  Susceptible-Infected-Recovered (SIR) model for the disease spread.
- `r pkg("covid19br")`: Set of functions to import COVID-19 pandemic data
  into R. The Brazilian COVID-19 data, obtained from the official Brazilian
  repository at <https://covid.saude.gov.br/>, is available at country,
  region, state, and city-levels. The package also downloads the world-level
  COVID-19 data from the John Hopkins University's repository.
- `r pkg("covid19dbcand")`: Provides different datasets parsed from
  '[Drugbank](https://www.drugbank.ca/covid-19)' database using
  `r pkg("dbparser")` package. It is a smaller version from
  `r github("interstellar-Consultation-Services/dbdataset")` package. It
  contains only information about COVID-19 possible treatment.
- `r pkg("covid19italy")`: Provides a daily summary of the Coronavirus
  (COVID-19) cases in Italy by country, region and province level. Data
  source: [Presidenza del Consiglio dei Ministri - Dipartimento della
  Protezione Civile](https://www.protezionecivile.it/).
- `r pkg("covid19france")`: Imports and cleans
  <https://github.com/opencovid19-fr/data> data on COVID-19 in France.
- `r github("covid19R/covid19mobility")`: Scrapes trends in mobility after the
  Covid-19 outbreak from different sources. Currently, the package scrapes
  data from Google (<https://www.google.com/covid19/mobility/>), Apple
  (<https://www.apple.com/covid19/mobility>), and will add others. The data
  returned uses the tidy [Covid19R project data
  standard](https://covid19r.github.io/documentation/) as well as the
  controlled vocabularies for measurement types.
- `r github("covid19R/covid19nytimes")`: Accesses the NY Times Covid-19
  county-level data for the US, described in
  <https://www.nytimes.com/article/coronavirus-county-data-us.html> and
  available at <https://github.com/nytimes/covid-19-data>. It then returns the
  data in a tidy data format according to the Covid19R Project data
  specification. If you plan to use the data or publicly display the data or
  results, please make sure cite the original NY Times source. Please read and
  follow the terms laid out in the data license at
  <https://github.com/nytimes/covid-19-data/blob/master/LICENSE>.
- `r pkg("covid19sf")`: Provides a verity of summary tables of the Covid19
  cases in San Francisco. Data source: [San Francisco, Department of Public
  Health - Population Health Division](https://datasf.org/opendata/).
- `r pkg("covid19swiss")`: Provides a daily summary of the Coronavirus
  (COVID-19) cases in Switzerland cantons and Principality of Liechtenstein.
  Data source: Specialist Unit for Open Government Data Canton of Zurich
  <https://www.zh.ch/de/politik-staat/opendata.html>.
- `r pkg("covid19us")`: A wrapper around the 'COVID Tracking Project API'
  <https://covidtracking.com/api/> providing data on cases of COVID-19 in the
  US.
- `r pkg("CovidMutations")`: A feasible framework for mutation analysis and
  reverse transcription polymerase chain reaction (RT-PCR) assay evaluation of
  COVID-19, including mutation profile visualization, statistics and mutation
  ratio of each assay. The mutation ratio is conducive to evaluating the
  coverage of RT-PCR assays in large-sized samples. [Mercatelli, D. and
  Giorgi, F. M. (2020)](https://doi.org/10.20944/preprints202004.0529.v1).
- `r pkg("CovidMutations")`: A feasible framework for mutation analysis and
  reverse transcription polymerase chain reaction (RT-PCR) assay evaluation of
  COVID-19, including mutation profile visualization, statistics and mutation
  ratio of each assay. The mutation ratio is conducive to evaluating the
  coverage of RT-PCR assays in large-sized samples. [Mercatelli, D. and
  Giorgi, F. M. (2020)](https://doi.org/10.20944/preprints202004.0529.v1).
- `r github("epiforecasts/covidregionaldata")`: An interface to subnational and national level
  COVID-19 data sourced from both official sources, such as Public Health
  England in the UK, and from other COVID-19 data collections, including
  the World Health Organisation (WHO), European Centre for Disease Prevention and
  Control (ECDC), John Hopkins University (JHU), Google Open Data and others.
  This package is designed to streamline COVID-19 data extraction, cleaning,
  and processing from a range of data sources in an open and transparent way.
  For all countries supported, data includes a daily time-series of cases and,
  wherever available, data on deaths, hospitalisations, and tests.

### Other data packages

- `r pkg("nhanesA")`: provides ready access to the National Health and
  Nutrition Examination Survey (NHANES) data tables.
- `r pkg("epimdr")` and `r pkg("epimdr2")`: companion packages with functions
  and data for "Epidemics: Models and Data in R" by O. Bjornstad, for 1st and
  2nd editions respectively.
- `r pkg("epiparameter")` and `r pkg("epiparameterDB")`: a harmonized database
  and standardized tools to manipulate published epidemiological parameters
  (serial interval, generation time, etc.) for a number of diseases.
- `r pkg("hmsidwR")`: A collection of datasets and supporting functions
  accompanying [Health Metrics and the Spread of Infectious Diseases by Federica
  Gazzelloni (2024)](https://doi.org/10.5281%2Fzenodo.10818338).

## Links

- The R Epidemics Consortium ([RECON](https://www.repidemicsconsortium.org)),
  a non-profit organization dedicated to the development of free, open-source
  outbreak analytics resources.

- [*Epiverse*](https://data.org/initiatives/epiverse/), an initiative created
  by [data.org](https://data.org) for the development of open-source resources
  for epidemic preparedness and
  response. [*Epiverse-TRACE*](https://github.com/epiverse-trace) is dedicated
  to creating an ecosystem of R packages for outbreak analytics.

- [The Epidemiologist R Handbook](https://epirhandbook.com/en/).
