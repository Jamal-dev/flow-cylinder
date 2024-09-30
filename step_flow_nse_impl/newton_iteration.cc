// newton_iteration.cc

#ifndef NEWTON_ITERATION_CC
#define NEWTON_ITERATION_CC

#include "../step_flow_nse.h"

// This is the Newton iteration
// to solve the non-linear system of equations:
//
// A'(U^k)(\delta U, \Psi) = - A(U^k)(\Psi)
//                 U^{k+1} = U^k + \lambda \delta U
//
// with $\lambda\in (0,1]$.
//
// For Newton's method as implemented here, please see
// for instance Chapter 13 in https://doi.org/10.15488/9248
//
// A'(U^k)(\delta U, \Psi) ^= assemble_system_matrix ()
// A(U^k)(\Psi)            ^= assemble_system_rhs ()
//
// The Jacobian A'(U^k)(\delta U, \Psi) is obtained
// by calculating the directional derivatives.
// The minus sign before the residual A(U^k)(\Psi)
// is contained in local_rhs(i) -= ... where we use
//   -= rather than +=
//
//
template <int dim>
void StepFlowNSE<dim>::newton_iteration()
{
    Timer timer_newton;
    const double lower_bound_newton_residuum = 1.0e-10;
    const unsigned int max_no_newton_steps = 20;

    // Decision whether the system matrix should be build
    // at each Newton step
    const double nonlinear_rho = 1e-10; // 0.1;

    // Line search parameters
    unsigned int line_search_step;
    const unsigned int max_no_line_search_steps = 10;
    const double line_search_damping = 0.6;
    double new_newton_residuum;

    // Application of the initial boundary conditions to the
    // variational equations:
    set_initial_bc();
    assemble_system_rhs();

    double newton_residuum = system_rhs.l2_norm();
    double old_newton_residuum = newton_residuum;
    double initial_newton_residuum = newton_residuum;
    unsigned int newton_step = 1;

    unsigned int stop_when_line_search_two_times_max_number = 0;

    if (newton_residuum < lower_bound_newton_residuum)
    {
        std::cout << '\t'
                  << std::scientific
                  << newton_residuum
                  << std::endl;
    }

    while ((newton_residuum > lower_bound_newton_residuum &&
            (newton_residuum / initial_newton_residuum) > lower_bound_newton_residuum) &&
           newton_step < max_no_newton_steps)
    {
        timer_newton.start();
        old_newton_residuum = newton_residuum;

        assemble_system_rhs();
        newton_residuum = system_rhs.l2_norm();

        if (newton_residuum < lower_bound_newton_residuum)
        {
            std::cout << '\t'
                      << std::scientific
                      << newton_residuum << std::endl;
            break;
        }

        // Check if matrix needs to re-build. If not
        // then we save some time and perform quasi-Newton steps.
        if (newton_residuum / old_newton_residuum > nonlinear_rho)
        {
            assemble_system_matrix();
            // Only factorize when matrix is re-built
            A_direct.factorize(system_matrix);
        }

        // Solve Ax = b
        solve();

        line_search_step = 0;
        for (;
             line_search_step < max_no_line_search_steps;
             ++line_search_step)
        {
            solution += newton_update;

            assemble_system_rhs();
            new_newton_residuum = system_rhs.l2_norm();

            if (new_newton_residuum < newton_residuum)
                break;
            else
                solution -= newton_update;

            newton_update *= line_search_damping;
        }

        if (line_search_step == 10)
            stop_when_line_search_two_times_max_number++;

        if (stop_when_line_search_two_times_max_number == 3)
        {
            std::cout << "Aborting Newton as line search does not help to converge anymore." << std::endl;
            abort();
        }

        timer_newton.stop();

        std::cout << std::setprecision(5) << newton_step << '\t'
                  << std::scientific << newton_residuum << '\t'
                  << std::scientific << newton_residuum / initial_newton_residuum << '\t'
                  << std::scientific << newton_residuum / old_newton_residuum << '\t';
        if (newton_residuum / old_newton_residuum > nonlinear_rho)
            std::cout << "r" << '\t';
        else
            std::cout << " " << '\t';
        std::cout << line_search_step << '\t'
                  << std::scientific << timer_newton.cpu_time()
                  << std::endl;

        // Updates
        timer_newton.reset();
        newton_step++;
    }
}

#endif // NEWTON_ITERATION_CC
