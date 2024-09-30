// run.cc

#ifndef RUN_CC
#define RUN_CC

#include "../step_flow_nse.h"

// As usual, we have to call the run method. It handles
// the output stream to the terminal.
// Second, we define some output skip that is necessary
// (and really useful) to avoid to much printing
// of solutions. For large time dependent problems it is
// sufficient to print only each tenth solution.
// Third, we perform the time stepping scheme of
// the solution process.
template <int dim>
void StepFlowNSE<dim>::run()
{
    // Definining test cases
    // test_case = "Channel_NSE";
    test_case = "2D-1_NSE";
    // test_case = "Cavity_IV";
    // test_case = "L_shaped_IV";
    // test_case = "Venturi_IV";

    setup_system();

    std::cout << "\n=============================="
              << "=====================================" << std::endl;
    std::cout << "Parameters\n"
              << "==========\n"
              << "Density flow: " << density_fluid << "\n"
              << "Viscosity:    " << kinematic_viscosity << "\n"
              << std::endl;
    // More output can be printed here if wished

    // Apply initial conditions
    {
        AffineConstraints<double> constraints;
        // ConstraintMatrix constraints;
        constraints.close();

        std::vector<bool> component_mask(dim + dim + 1, true);
        VectorTools::project(dof_handler,
                             constraints,
                             QGauss<dim>(degree + 2),
                             InitialValues<dim>(),
                             solution);

        // TODO: eher keine Loesung schreiben, da exp(-gamma) noch nicht initialisiert.
        // output_results (timestep_number,solution);
    }

    // Increment time after having prescribed the initial solution
    time += timestep;
    ++timestep_number;

    // Define a flag that not each *.vtk is printed
    const unsigned int output_skip = 1;

    // Time loop
    do
    {
        std::cout << "Timestep " << timestep_number
                  << " (" << time_stepping_scheme
                  << ")" << ": " << time
                  << " (" << timestep << ")"
                  << "\n=============================="
                  << "====================================================="
                  << std::endl;

        std::cout << std::endl;

        // Compute next time step solution by solving
        // the nonlinear system
        old_old_timestep_solution = old_timestep_solution;
        old_timestep_solution = solution;

        // Solve nonlinear problem
        newton_iteration();

        // Compute functional values
        std::cout << std::endl;
        compute_functional_values();

        // Write solutions
        if ((timestep_number % output_skip == 0))
            output_results(timestep_number, solution);

        // Update time
        time += timestep;
        ++timestep_number;

    } while (timestep_number <= max_no_timesteps);
}

#endif // RUN_CC
