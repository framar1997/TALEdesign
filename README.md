# Simulation studies using the TALE design.

The TALE design is a simple algorithm-based approach for phase-I dual-agent cancer trials. The main R Shiny web application
to conduct simulation studies is accessible at [TALEdesign]( https://6kp5ow-francesco-mariani.shinyapps.io/TALEdesign/.).
This repository contains two additional tools that offer some advantages compared to the TALEdesign application:

- [Parallelized application](https://github.com/framar1997/TALEdesign/blob/main/TALE_application_parallelized.R).
  Free shiny subscription imposes some limitations on computational resources. For faster simulations, we here provide the R code of the TALEdesign application
  adapted to execute simulations in parallel;
- [Core functions](https://github.com/framar1997/TALEdesign/tree/main/core_functions).
