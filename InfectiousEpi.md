---
name: InfectiousEpi
topic: Infectious Disease Epidemiology
maintainer: Thibaut Jombart, Matthieu Rolland, Hugo Gruson
email: thibautjombart@gmail.com
version: 2021-02-25
---

Contributors (in alphabetic order): Neale Batra, Rich FitzJohn, Hugo Gruson,
Andreas Handel, Michael Höhle, Thibaut Jombart, Joseph Larmarange, Sebastian
Lequime, Alex Spina, Tim Taylor, Sean Wu

## Overview

R is increasingly becoming a standard in infectious disease epidemiology (IDE),
providing a wide array of tools covering different aspects of outbreak analytics
from study design to data exploration, modelling, forecasting and simulation.
This task view provides an overview of packages specifically developed for
IDE. It does not encompass other applications of epidemiology (such as
environmental epidemiology) or cross-cutting tools used in, but not specific to,
IDE.

Packages were regrouped in the following categories:

1.  **Data exploration and visualisation:** these tools are dedicated to
    handling and visualising IDE data, producing epidemic curves
    ('*epicurves*'), exploring contact tracing data, etc.
2.  **Infectious disease modelling:** these tools encompass various approaches
    for the analysis of epidemic curves (including outbreak detection /
    surveillance), estimating reproduction numbers (*R*), providing short-term
    forecasting, implementing compartmental models (e.g. SIR models), simulating
    outbreaks and reconstructing transmission trees
3.  **Helpers**: tools implementing miscellaneous tasks useful in IDE
4.  **Data packages**: these packages provide access to both empirical and
    simulated epidemic data; includes a specific section on COVID-19.

Additional links to non specific but highly useful packages (to create tables,
manipulate dates, etc.) are provided in the task view's footnotes.

## Inclusion criteria

Packages included in this task view were identified through recommendations of
IDE experts as well as an automated CRAN search using `pkgsearch::pkg_search()`,
with the keywords: *epidemiology*, *epidemic*, *epi*, *outbreak* and
*transmission*. The list was manually curated for the final selection to satisfy
the conditions described in the previous paragraph.

Packages are deemed in scope if they provide tools, or data, explicitly targeted
at reporting, modelling, or forecasting infectious diseases.

**Your input is welcome!** Please suggest packages we may have missed by posting
an issue at:

https://github.com/bisaloo/InfectiousEpi


## Data visualisation and exploration

Here are packages providing help for visualising and exploring epidemic data
(plotting epidemic curves, incidence curves, contact tracing, etc.)

-   `r pkg("epicontacts")`: A collection of tools for representing
    epidemiological contact data, composed of case line lists and contacts
    between cases. Also contains procedures for data handling, interactive
    graphics, and statistics.
-   `r pkg("EpiContactTrace")`: Routines for epidemiological contact tracing and
    visualisation of network of contacts.
-   `r pkg("EpiCurve")`: Creates simple or stacked epidemic curves for hourly,
    daily, weekly or monthly outcome data.
-   `r pkg("epiDisplay")`: Package for data exploration and result presentation.
-   `r pkg("epiflows")`: Provides functions and classes designed to handle and
    visualise epidemiological flows of people between locations. Also contains a
    statistical method for predicting disease spread from flow data initially
    described in [Dorigatti et al.
    (2017)](https://doi.org/10.2807%2F1560-7917.ES.2017.22.28.30572). This
    package is part of the [RECON](https://www.repidemicsconsortium.org/)
    toolkit for outbreak analysis.
-   `r pkg("EpiReport")`: Drafting an epidemiological report in 'Microsoft Word'
    format for a given disease, similar to the Annual Epidemiological Reports
    published by the European Centre for Disease Prevention and Control.
-   `r pkg("i2extras")`: Provides functions to work with 'incidence2' objects,
    including a simplified interface for trend fitting, estimation of growth
    rates, and peak estimation. This package is part of the
    [RECON](https://www.repidemicsconsortium.org/) toolkit for outbreak
    analysis.
-   `r pkg("incidence")`: Provides functions and classes to compute, handle and
    visualise incidence from dated events for a defined time interval. Dates can
    be provided in various standard formats. The class 'incidence' is used to
    store computed incidence and can be easily manipulated, subsetted, and
    plotted. In addition, log-linear models can be fitted to 'incidence' objects
    using 'fit'. This package is part of the
    [RECON](https://www.repidemicsconsortium.org/) toolkit for outbreak
    analysis, but is now scheduled for deprecation and is replaced by
    `r pkg("incidence2")`.
-   `r pkg("incidence2")`: Provides functions and classes to compute, handle and
    visualise incidence from dated events for a defined time interval. Improves
    the original `r pkg("incidence")` package in many ways: full flexibility in
    time intervals used, allows multiple stratification, and is fully compatible
    with `r pkg("dplyr")` and other tidyverse tools. This package is part of the
    [RECON](https://www.repidemicsconsortium.org/) toolkit for outbreak
    analysis.

## Infectious disease modelling

### Epidemics surveillance

-   `r pkg("argo")`: Augmented Regression with General Online data (ARGO) for
    accurate estimation of influenza epidemics in United States on both national
    level and regional level. It replicates the method introduced in paper
    [Yang, S., Santillana, M. and Kou, S.C.
    (2015)](https://doi.org/10.1073/pnas.1515373112) and [Ning, S., Yang, S. and
    Kou, S.C. (2019)](https://doi.org/10.1038/s41598-019-41559-6).
-   `r pkg("Epi")`: Functions for demographic and epidemiological analysis in
    the Lexis diagram, i.e. register and cohort follow-up data, in particular
    representation, manipulation and simulation of multistate data - the Lexis
    suite of functions, which includes interfaces to `r pkg("mstate")`,
    `r pkg("etm")` and `r pkg("cmprsk")` packages. Also contains functions for
    Age-Period-Cohort and Lee-Carter modeling and a function for interval
    censored data and some useful functions for tabulation and plotting, as well
    as a number of epidemiological data sets.
-   `r pkg("episensr")`: Basic sensitivity analysis of the observed relative
    risks adjusting for unmeasured confounding and misclassification of the
    exposure/outcome, or both. It follows the bias analysis methods and examples
    from the book by Lash T.L, Fox M.P, and Fink A.K. "Applying Quantitative
    Bias Analysis to Epidemiologic Data", ('Springer', 2009).
-   `r pkg("mem")`: The Moving Epidemic Method, created by T Vega and JE Lozano
    ([2012](https://doi.org/10.1111/j.1750-2659.2012.00422.x),
    [2015](https://doi.org/10.1111/irv.12330)), allows the weekly assessment of
    the epidemic and intensity status to help in routine respiratory infections
    surveillance in health systems. Allows the comparison of different epidemic
    indicators, timing and shape with past epidemics and across different
    regions or countries with different surveillance systems. Also, it gives a
    measure of the performance of the method in terms of sensitivity and
    specificity of the alert week.
-   `r pkg("memapp")`: The Moving Epidemic Method, created by T Vega and JE
    Lozano ([2012](https://doi.org/10.1111/j.1750-2659.2012.00422.x),
    [2015](https://doi.org/10.1111/irv.12330)), allows the weekly assessment of
    the epidemic and intensity status to help in routine respiratory infections
    surveillance in health systems. Allows the comparison of different epidemic
    indicators, timing and shape with past epidemics and across different
    regions or countries with different surveillance systems. Also, it gives a
    measure of the performance of the method in terms of sensitivity and
    specificity of the alert week. 'memapp' is a web application created in the
    Shiny framework for the `r pkg("mem")` R package.
-   `r pkg("nosoi")`: The aim of `r pkg("nosoi")` (pronounced no.si) is to
    provide a flexible agent-based stochastic transmission chain/epidemic
    simulator ([Lequime et al. 2020](https://doi.org/2020.03.03.973107)). It is
    named after the daimones of plague, sickness and disease that escaped
    Pandora's jar in the Greek mythology. `r pkg("nosoi")` is able to take into
    account the influence of multiple variable on the transmission process (e.g.
    dual-host systems (such as arboviruses), within-host viral dynamics,
    transportation, population structure), alone or taken together, to create
    complex but relatively intuitive epidemiological simulations.
-   `r pkg("riskCommunicator")`: Estimates flexible epidemiological effect
    measures including both differences and ratios using the parametric
    G-formula developed as an alternative to inverse probability weighting. It
    is useful for estimating the impact of interventions in the presence of
    treatment-confounder-feedback. G-computation was originally described by
    [Robbins (1986)](https://doi.org/10.1016/0270-0255(86)90088-6) and has been
    described in detail by [Ahern, Hubbard, and Galea
    (2009)](https://doi.org/10.1093/aje/kwp015); [Snowden, Rose, and Mortimer
    (2011)](https://doi.org/10.1093/aje/kwq472); and [Westreich et al.
    (2012)](https://doi.org/10.1002/sim.5316).
-   `r pkg("RSurveillance")`: Provides a diverse set of functions useful for the
    design and analysis of disease surveillance activities.
-   `r pkg("trendeval")`: Provides a coherent interface for evaluating models
    fit with the trending package. This package is part of the
    [RECON](https://www.repidemicsconsortium.org/) toolkit for outbreak
    analysis.
-   `r pkg("SpatialEpi")`: Methods and data for cluster detection and disease
    mapping.
-   `r pkg("surveillance")`: Statistical methods for the modeling and monitoring
    of time series of counts, proportions and categorical data, as well as for
    the flexible covariate based regression modeling of spatio-temporal point
    processes of epidemic phenomena. The monitoring methods focus on aberration
    detection in count data time series from public health surveillance of
    communicable diseases. A recent overview of the available monitoring
    procedures is given by [Salmon et al.
    (2016)](https://doi.org/10.18637%2Fjss.v070.i10) and a recent overview of
    the implemented space-time modeling frameworks for epidemic phenomena is
    given by [Meyer et al. (2017)](https://doi.org/10.18637%2Fjss.v077.i11).
    Other package feature contain back-projection methods to infer time series
    of exposure from disease onset and correction of observed time series for
    reporting delays (nowcasting).
-   `r github("reconhub/trendbreaker")` implements tools for detecting changes
    in temporal trends of a single response variable. It implements the
    **A**utomatic **S**election of **M**odels and **O**utlier **De**tection for
    **E**pidemmics (ASMODEE), an algorithm originally designed for detecting
    changes in COVID-19 case incidence.
-   `r pkg("trending")`: Provides a coherent interface to multiple modelling
    tools for fitting trends along with a standardised approach for generating
    confidence and prediction intervals.

### Estimation of transmissibility

-   `r pkg("earlyR")`: Implements a simple, likelihood-based estimation of the
    reproduction number (R0) using a branching process with a Poisson
    likelihood. This model requires knowledge of the serial interval
    distribution, and dates of symptom onsets. Infectiousness is determined by
    weighting R0 by the probability mass function of the serial interval on the
    corresponding day. It is a simplified version of the model introduced by
    [Cori et al. (2013)](https://doi.org/10.1093%2Faje%2Fkwt133).
-   `r pkg("endtoend")`: Computes the expectation of the number of transmissions
    and receptions considering an End-to-End transport model with limited number
    of retransmissions per packet. It provides theoretical results and also
    estimated values based on Monte Carlo simulations. It is also possible to
    consider random data and ACK probabilities.
-   `r pkg("EpiEstim")`: Tools to quantify transmissibility throughout an
    epidemic from the analysis of time series of incidence as described in [Cori
    et al. (2013)](https://doi.org/10.1093/aje/kwt133) and [Wallinga and Teunis
    (2004)](https://doi.org/10.1093/aje/kwh255).
-   `r pkg("epimdr")`: Functions, data sets and shiny apps for "Epidemics:
    Models and Data in R" by Ottar N. Bjornstad ([ISBN
    978-3-319-97487-3](https://www.springer.com/gp/book/9783319974866)). The
    package contains functions to study the S(E)IR model, spatial and
    age-structured SIR models; time-series SIR and chain-binomial stochastic
    models; catalytic disease models; coupled map lattice models of spatial
    transmission and network models for social spread of infection. The package
    is also an advanced quantitative companion to the [coursera Epidemics
    Massive Online Open Course](https://www.coursera.org/learn/epidemics).
-   `r pkg("epinet")`: A collection of epidemic/network-related tools. Simulates
    transmission of diseases through contact networks. Performs Bayesian
    inference on network and epidemic parameters, given epidemic data.
-   `r github("epiforecasts/EpiNow2")`: This package estimates the time-varying
    reproduction number, rate of spread, and doubling time using a range of
    open-source tools ([Abbott et
    al.](https://doi.org/10.12688/wellcomeopenres.16006.1)), and current best
    practices ([Gostic et al.](https://doi.org/10.1371/journal.pcbi.1008409)).
    It aims to help users avoid some of the limitations of naive implementations
    in a framework that is informed by community feedback and is under active
    development.
-   `r pkg("nbTransmission")`: Estimates the relative transmission probabilities
    between cases in an infectious disease outbreak or cluster using naive
    Bayes. Included are various functions to use these probabilities to estimate
    transmission parameters such as the generation/serial interval and
    reproductive number as well as finding the contribution of covariates to the
    probabilities and visualizing results. The ideal use is for an infectious
    disease dataset with metadata on the majority of cases but more informative
    data such as contact tracing or pathogen whole genome sequencing on only a
    subset of cases. For a detailed description of the methods see [Leavitt et
    al. (2020)](https://doi.org/10.1093%2Fije%2Fdyaa031).
-   `r pkg("R0")`: Estimation of reproduction numbers for disease outbreak,
    based on incidence data. The R0 package implements several documented
    methods. It is therefore possible to compare estimations according to the
    methods used. Depending on the methods requested by user, basic reproduction
    number (commonly denoted as R0) or real-time reproduction number (referred
    to as R(t)) is computed, along with a 95% Confidence Interval. Plotting
    outputs will give different graphs depending on the methods requested :
    basic reproductive number estimations will only show the epidemic curve
    (collected data) and an adjusted model, whereas real-time methods will also
    show the R(t) variations throughout the outbreak time period. Sensitivity
    analysis tools are also provided, and allow for investigating effects of
    varying Generation Time distribution or time window on estimates.
-   `r pkg("tsiR")`: The TSIR modeling framework allows users to fit the time
    series SIR model to cumulative case data, which uses a regression equation
    to estimate transmission parameters based on differences in cumulative cases
    between reporting periods. The package supports inference on TSIR parameters
    using GLMs and profile likelihood techniques, as well as forward simulation
    based on a fitted model. The package is described in [Becker and Grenfell
    (2017)](https://doi.org/10.1371/journal.pone.0185528).

### Compartmental models

-   `r pkg("EpiILM")`: Provides tools for simulating from discrete-time
    individual level models for infectious disease data analysis. This epidemic
    model class contains spatial and contact-network based models with two
    disease types: Susceptible-Infectious (SI) and
    Susceptible-Infectious-Removed (SIR).
-   `r pkg("EpiILMCT")`: Provides tools for simulating from continuous-time
    individual level models of disease transmission, and carrying out infectious
    disease data analyses with the same models. The epidemic models considered
    are distance-based and/or contact network-based models within
    Susceptible-Infectious-Removed (SIR) or
    Susceptible-Infectious-Notified-Removed (SINR) compartmental frameworks. An
    overview of the implemented continuous-time individual level models for
    epidemics is given by [Almutiry and Deardon
    (2019)](https://doi.org/10.1515/ijb-2017-0092).
-   `r pkg("EpiModel")`: Tools for simulating mathematical models of infectious
    disease dynamics. Epidemic model classes include deterministic compartmental
    models, stochastic individual-contact models, and stochastic network models.
    Network models use the robust statistical methods of exponential-family
    random graph models (ERGMs) from the Statnet suite of software packages
    in R. Standard templates for epidemic modeling include SI, SIR, and SIS
    disease types. EpiModel features an API for extending these templates to
    address novel scientific research aims. Full methods for EpiModel are
    detailed in [Jenness et al. (2018)](https://doi.org/10.18637/jss.v084.i08).
-   `r pkg("odin")`: Generate systems of ordinary differential equations (ODE)
    and integrate them, using a domain specific language (DSL). The DSL uses R's
    syntax, but compiles to C in order to efficiently solve the system. A solver
    is not provided, but instead interfaces to the packages `r pkg("deSolve")`
    and `r pkg("dde")` are generated. With these, while solving the differential
    equations, no allocations are done and the calculations remain entirely in
    compiled code. Alternatively, a model can be transpiled to R for use in
    contexts where a C compiler is not present. After compilation, models can be
    inspected to return information about parameters and outputs, or
    intermediate values after calculations. `r pkg("odin")` is not targeted at
    any particular domain and is suitable for any system that can be expressed
    primarily as mathematical expressions. Additional support is provided for
    working with delays (delay differential equations, DDE), using interpolated
    functions during interpolation, and for integrating quantities that
    represent arrays.
-   `r pkg("pomp")` Provides a large set of forward simulation algorithms and
    MLE or Bayesian inference techniques to work with state-space models. Models
    may be specified as either deterministic or stochastic and typically follow
    a compartmental model structure. Time may be either discrete or continuous,
    depending on the simulation algorithm chosen. Additionally models may be
    programmed in C and compiled on-the-fly into a format suitable for use in
    the package to speed up simulation and inference. The R package and some of
    its algorithms are described in [King, Nguyen, and Ionides
    (2016)](https://doi.org/10.18637/jss.v069.i12).
-   `r pkg("popEpi")`: Enables computation of epidemiological statistics,
    including those where counts or mortality rates of the reference population
    are used. Currently supported: excess hazard models, rates, mean survival
    times, relative survival, and standardized incidence and mortality ratios
    (SIRs/SMRs), all of which can be easily adjusted for by covariates such as
    age. Fast splitting and aggregation of 'Lexis' objects (from package
    `r pkg("Epi")` and other computations achieved using `r pkg("data.table")`.
-   `r pkg("SimInf")`: Provides an efficient and very flexible framework to
    conduct data-driven epidemiological modeling in realistic large scale
    disease spread simulations. The framework integrates infection dynamics in
    subpopulations as continuous-time Markov chains using the Gillespie
    stochastic simulation algorithm and incorporates available data such as
    births, deaths and movements as scheduled events at predefined time-points.
    Using C code for the numerical solvers and 'OpenMP' (if available) to divide
    work over multiple processors ensures high performance when simulating a
    sample outcome. One of our design goals was to make the package extendable
    and enable usage of the numerical solvers from other R extension packages in
    order to facilitate complex epidemiological research. The package contains
    template models and can be extended with user-defined models. For more
    details see the paper by [Widgren, Bauer, Eriksson and Engblom
    (2019)](https://doi.org/10.18637%2Fjss.v091.i12).
-   `r pkg("socialmixr")`: Provides methods for sampling contact matrices from
    diary data for use in infectious disease modelling, as discussed in [Mossong
    et al. (2008)](https://doi.org/10.1371%2Fjournal.pmed.0050074).

### Transmission tree reconstruction

-   `r pkg("adegenet")`: while primarily a population genetics package,
    `r pkg("adegenet")` also implements seqtrack ([Jombart et al.
    2011](https://doi.org/10.1038%2Fhdy.2010.78)), a maximum-parsimony approach
    for reconstructing transmission trees using the Edmonds/Chu-Liu algorithm
-   `r pkg("o2geosocial")` (`r pkg("outbreaker2")` module): Bayesian
    reconstruction of who infected whom during past outbreaks using
    routinely-collected surveillance data. Inference of transmission trees using
    genotype, age specific social contacts, distance between cases and onset
    dates of the reported cases ([Robert A, Kucharski AJ, Gastanaduy PA, Paul P,
    Funk S. 2020](https://doi.org/10.1098%2Frsif.2020.0084)).
-   `r github("xavierdidelot/o2mod.transphylo")` is a module of
    `r pkg("outbreaker2")` which uses the `r pkg("TransPhylo")` model of
    within-host evolution.
-   `r pkg("outbreaker2")`: Bayesian reconstruction of disease outbreaks using
    epidemiological and genetic information. [Jombart T, Cori A, Didelot X,
    Cauchemez S, Fraser C and Ferguson N.
    2014](https://doi.org/10.1371/journal.pcbi.1003457), [Campbell F, Cori A,
    Ferguson N, Jombart T. 2019](https://doi.org/10.1371/journal.pcbi.1006930).
-   `r pkg("TransPhylo")`: Inference of transmission tree from a dated
    phylogeny. Includes methods to simulate and analyse outbreaks. The
    methodology is described in [Didelot et al.
    (2014)](https://doi.org/10.1093%2Fmolbev%2Fmsu121), [Didelot et al.
    (2017)](https://doi.org/10.1093%2Fmolbev%2Fmsw275).

### Vaccination

-   `r pkg("SCCS")`: Self-controlled case series models used to investigate
    associations between time-varying exposures such as vaccines or other drugs
    or non drug exposures and an adverse event can be fitted. Detailed
    information on the self-controlled case series method and its extensions
    with more examples at <https://sccs-studies.info>.

## Helper functions

Here are packages providing useful helper functions for practicing, teaching and
learning epidemiology (eg computing sample size, contingency tables, etc).

-   `r pkg("DSAIDE")`: Exploration of simulation models (apps) of various
    infectious disease transmission dynamics scenarios. The purpose of the
    package is to help individuals learn about infectious disease epidemiology
    (ecology/evolution) from a dynamical systems perspective. All apps include
    explanations of the underlying models and instructions on what to do with
    the models.
-   `r pkg("epibasix")`: Contains elementary tools for analysis of common
    epidemiological problems, ranging from sample size estimation, through 2x2
    contingency table analysis and basic measures of agreement (kappa,
    sensitivity/specificity). Appropriate print and summary statements are also
    written to facilitate interpretation wherever possible. The target audience
    includes advanced undergraduate and graduate students in epidemiology or
    biostatistics courses, and clinical researchers.
-   `r pkg("epiR")`: Tools for the analysis of epidemiological and surveillance
    data. Contains functions for directly and indirectly adjusting measures of
    disease frequency, quantifying measures of association on the basis of
    single or multiple strata of count data presented in a contingency table,
    computation of confidence intervals around incidence risk and incidence rate
    estimates and sample size calculations for cross-sectional, case-control and
    cohort studies. Surveillance tools include functions to calculate an
    appropriate sample size for 1- and 2-stage representative freedom surveys,
    functions to estimate surveillance system sensitivity and functions to
    support scenario tree modelling analyses.
-   `r pkg("epitools")`: Tools for training and practicing epidemiologists
    including methods for two-way and multi-way contingency tables.
-   `r pkg("epitrix")`: A collection of small functions useful for epidemics
    analysis and infectious disease modelling. This includes computation of
    basic reproduction numbers from growth rates, generation of hashed labels to
    anonymise data, and fitting discretised Gamma distributions.
-   `r pkg("powerSurvEpi")`: Functions to calculate power and sample size for
    testing main effect or interaction effect in the survival analysis of
    epidemiological studies (non-randomized studies), taking into account the
    correlation between the covariate of the interest and other covariates. Some
    calculations also take into account the competing risks and stratified
    analysis. This package also includes a set of functions to calculate power
    and sample size for testing main effects in the survival analysis of
    randomized clinical trials.

## Data

Here are packages providing different epidemiologic datasets, either simulated
or real, useful for research purposes or field applications with a specific
COVID-19 section.

### Epidemic outbreak data

-   `r pkg("contactdata")`: Data package for the supplementary data in [Prem et
    al. (2017)](https://doi.org/10.1371%2Fjournal.pcbi.1005697). Provides easy
    access to contact data for 152 countries, for use in epidemiological,
    demographic or social sciences research.
-   `r pkg("outbreaks")`: Empirical or simulated disease outbreak data, provided
    either as RData or as text files.

### COVID-19

-   `r pkg("bets.covid19")`: Implements likelihood inference for early epidemic
    analysis. BETS is short for the four key epidemiological events being
    modeled: Begin of exposure, End of exposure, time of Transmission, and time
    of Symptom onset. The package contains a dataset of the trajectory of
    confirmed cases during the coronavirus disease (COVID-19) early outbreak.
    More detail of the statistical methods can be found in [Zhao et al.
    (2020)](https://arxiv.org/abs/2004.07743).
-   `r pkg("corona")`: Manipulate and view coronavirus data and other societally
    relevant data at a basic level.
-   `r pkg("coronavirus")`: Provides a daily summary of the Coronavirus
    (COVID-19) cases by state/province. Data source: [Johns Hopkins University
    Center for Systems Science and Engineering (JHU CCSE)
    Coronavirus](https://systems.jhu.edu/research/public-health/ncov/).
-   `r pkg("COVID19")`: Download COVID-19 data across governmental sources at
    national, regional, and city level, as described in [Guidotti and Ardia
    (2020)](https://doi.org/10.21105/joss.02376). Includes the time series of
    vaccines, tests, cases, deaths, recovered, hospitalizations, intensive
    therapy, and policy measures by '[Oxford COVID-19 Government Response
    Tracker](https://www.bsg.ox.ac.uk/research/research-projects/coronavirus-government-response-tracker).
    Provides a seamless integration with '[World Bank Open
    Data](https://data.worldbank.org/)', '[Google Mobility
    Reports](https://www.google.com/covid19/mobility/)', 'Apple Mobility
    Reports](<https://covid19.apple.com/mobility>)'.
-   `r pkg("covid19.analytics")`: Load and analyze updated time series worldwide
    data of reported cases for the Novel CoronaVirus Disease (CoViD-19) from the
    [Johns Hopkins University Center for Systems Science and Engineering (JHU
    CSSE) data repository](https://github.com/CSSEGISandData/COVID-19). The
    datasets are available in two main modalities, as a time series sequences
    and aggregated for the last day with greater spatial resolution. Several
    analysis, visualization and modelling functions are available in the package
    that will allow the user to compute and visualize total number of cases,
    total number of changes and growth rate globally or for an specific
    geographical location, while at the same time generating models using these
    trends; generate interactive visualizations and generate
    Susceptible-Infected-Recovered (SIR) model for the disease spread.
-   `r pkg("covid19br")`: Set of functions to import COVID-19 pandemic data
    into R. The Brazilian COVID-19 data, obtained from the official Brazilian
    repository at <https://covid.saude.gov.br/>, is available at country,
    region, state, and city-levels. The package also downloads the world-level
    COVID-19 data from the John Hopkins University's repository.
-   `r pkg("covid19dbcand")`: Provides different datasets parsed from
    '[Drugbank](https://www.drugbank.ca/covid-19)' database using
    `r pkg("dbparser")` package. It is a smaller version from
    `r github("interstellar-Consultation-Services/dbdataset")` package. It
    contains only information about COVID-19 possible treatment.
-   `r pkg("covid19italy")`: Provides a daily summary of the Coronavirus
    (COVID-19) cases in Italy by country, region and province level. Data
    source: [Presidenza del Consiglio dei Ministri - Dipartimento della
    Protezione Civile](http://www.protezionecivile.it/).
-   `r pkg("covid19france")`: Imports and cleans
    <https://github.com/opencovid19-fr/data> data on COVID-19 in France.
-   `r github("nevrome/covid19germany")`: An R package to load, visualise and
    analyse daily updated data on the COVID-19 outbreak in Germany.
-   `r github("covid19R/covid19mobility")`: Scrapes trends in mobility after the
    Covid-19 outbreak from different sources. Currently, the package scrapes
    data from Google (<https://www.google.com/covid19/mobility/>), Apple
    (<https://www.apple.com/covid19/mobility>), and will add others. The data
    returned uses the tidy [Covid19R project data
    standard](https://covid19r.github.io/documentation/) as well as the
    controlled vocabularies for measurement types.
-   `r github("covid19R/covid19nytimes")`: Accesses the NY Times Covid-19
    county-level data for the US, described in
    <https://www.nytimes.com/article/coronavirus-county-data-us.html> and
    available at <https://github.com/nytimes/covid-19-data>. It then returns the
    data in a tidy data format according to the Covid19R Project data
    specification. If you plan to use the data or publicly display the data or
    results, please make sure cite the original NY Times source. Please read and
    follow the terms laid out in the data license at
    <https://github.com/nytimes/covid-19-data/blob/master/LICENSE>.
-   `r pkg("covid19sf")`: Provides a verity of summary tables of the Covid19
    cases in San Francisco. Data source: [San Francisco, Department of Public
    Health - Population Health Division](https://datasf.org/opendata/).
-   `r pkg("covid19swiss")`: Provides a daily summary of the Coronavirus
    (COVID-19) cases in Switzerland cantons and Principality of Liechtenstein.
    Data source: Specialist Unit for Open Government Data Canton of Zurich
    <https://www.zh.ch/de/politik-staat/opendata.html>.
-   `r pkg("covid19us")`: A wrapper around the 'COVID Tracking Project API'
    <https://covidtracking.com/api/> providing data on cases of COVID-19 in the
    US.
-   `r pkg("CovidMutations")`: A feasible framework for mutation analysis and
    reverse transcription polymerase chain reaction (RT-PCR) assay evaluation of
    COVID-19, including mutation profile visualization, statistics and mutation
    ratio of each assay. The mutation ratio is conducive to evaluating the
    coverage of RT-PCR assays in large-sized samples. [Mercatelli, D. and
    Giorgi, F. M. (2020)](https://doi.org/10.20944/preprints202004.0529.v1).
-   `r pkg("oxcgrt")`: The Oxford COVID-19 Government Response Tracker (OxCGRT)
    tracks and compares worldwide government responses to the COVID-19 pandemic
    rigorously and consistently. OxCGRT makes available systematic information
    in a consistent way, aiding those who require information have access to it
    efficiently for their purposes. This package facilitates access to the
    OxCGRT data via its API <https://covidtracker.bsg.ox.ac.uk/> and includes
    functions to calculate the various OxCGRT indices in R.
-   `r github("como-ph/oxcovid19")`: The [OxCOVID19
    Project](https://covid19.eng.ox.ac.uk) aims to increase our understanding of
    the COVID-19 pandemic and elaborate possible strategies to reduce the impact
    on the society through the combined power of statistical, mathematical
    modelling, and machine learning techniques. The OxCOVID19 Database is a
    large, single-centre, multimodal relational database consisting of
    information (using acknowledged sources) related to COVID-19 pandemic. This
    package provides an R-specific interface to the OxCOVID19 Database based on
    widely-used data handling and manipulation approaches in R.

## Links

-   R Epidemics Consortium (RECON), non-profit organisation developing
    professional tools in for field epidemiologists:
    <https://repidemicsconsortium.org/projects/>

-   RECON's task manager, issues regarding infectious epidemiology packages from
    RECON packages and partners: <https://tasks.repidemicsconsortium.org>

-   Survival analysis task view: `r view("Survival")`

-   Data visualisation and graphics task view: `r view("Graphics")`

-   Meta analysis task view: `r view("MetaAnalysis")`

-   Multivariate statistics task view: `r view("Multivariate")`
