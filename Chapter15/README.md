# Content

Contains several OpenACC version of the same code, which mimics the structure
of a simplified atmospheric model having only physical parameterizations
and output. For reference a serial and OpenMP versions are also provided.

The code is structured as follows:

* main.f90                 main driver, calls init, time loop, and output
* m_config.f90             configuration information domain size, number of
* steps
* m_fields.f90             global fields
* m_io.f90                 output routine
* m_parameterization.f90   physical parameterizations doing the actual
* computation
* m_physics.f90            driver for the physical paramtrization
* m_setup.f90              code initialization and clean up
* m_timing                 timing routines

Code versions:

* example_serial : serial code
* example_openmp : parallel openmp code
* example_openacc_step1: only OpenACC parallel region 
* example_openacc_step2: with parallel and data region
* example_openacc_step3: further optimization
* example_openacc_alt1: alternate version using openacc routine
* example_openacc_alt2: alternate version calling a cuda function


# Compilation

make example_<name>

For example, in order to compile the serial version of the example, you
can execute the following command:

make example_serial

This creates an executable example_serial/example_serial.

