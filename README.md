# Simulation studies using the TALE design

The TALE design is a simple algorithm-based approach for phase-I drug combination trials in oncology. The main R Shiny web application
to conduct simulation studies is accessible at [TALEdesign]( https://6kp5ow-francesco-mariani.shinyapps.io/TALEdesign/.).
This repository contains two additional tools that offer some advantages compared to the TALEdesign application:

- [TALEdesign parallelized](https://github.com/framar1997/TALEdesign/blob/main/TALE_application_parallelized.R) -
  Free shiny subscription imposes some limitations on computational resources. For faster simulations, we here provide the R code of the TALEdesign application
  adapted to execute simulations in parallel on the users local machine.

  
- [Core functions](https://github.com/framar1997/TALEdesign/tree/main/core_functions) - In the repository we made available the core functions of the Shiny app so that they can be used as an R package. By editing these functions, users can customize the algorithm characteristics according to their needs.
  We also provide a [tutorial](https://github.com/framar1997/TALEdesign/blob/main/simulation_tutorial.qmd) containing guidelines to simulate a single or multiple trials using
  the provided core functions.

The availability of these additional tools is also meant to provide a clear understanding of the computational procedures, as well as to ensure a reproducible R implementation.
