// setup_system.cc

#ifndef SETUP_SYSTEM_CC
#define SETUP_SYSTEM_CC

#include "../step_flow_nse.h"

template <int dim>
void StepFlowNSE<dim>::setup_system()
{
    timer.enter_subsection("Setup system.");

    // We set runtime parameters to drive the problem.
    // These parameters could also be read from a parameter file that
    // can be handled by the ParameterHandler object (see step-19)
    set_runtime_parameters();

    system_matrix.clear();

    dof_handler.distribute_dofs(fe);
    DoFRenumbering::Cuthill_McKee(dof_handler);

    // We are dealing with 5 components for this
    // two-dimensional Biot problem
    // (the first component is actually unnecessary
    // and a relict from the original fluid-structure
    // interaction code. But this component
    // might be used for future extensions with other equations.
    // velocity in x and y:                0
    // solid displacement in x and y:      1
    // scalar pressure field:              2
    std::vector<unsigned int> block_component(5, 0);
    block_component[dim] = 1;
    block_component[dim + 1] = 1;
    block_component[dim + dim] = 2;

    DoFRenumbering::component_wise(dof_handler, block_component);

    {
        constraints.clear();
        set_newton_bc();
        DoFTools::make_hanging_node_constraints(dof_handler,
                                                constraints);
    }
    constraints.close();

    std::vector<unsigned int> dofs_per_block(3);
    dofs_per_block = DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_v = dofs_per_block[0],
                       n_u = dofs_per_block[1],
                       n_p = dofs_per_block[2]; // temp

    std::cout << "Cells:\t"
              << triangulation.n_active_cells()
              << std::endl
              << "DoFs:\t"
              << dof_handler.n_dofs()
              << " (" << n_v << '+' << n_u << '+' << n_p << ')'
              << std::endl;

    {
        BlockDynamicSparsityPattern csp(3, 3);

        csp.block(0, 0).reinit(n_v, n_v);
        csp.block(0, 1).reinit(n_v, n_u);
        csp.block(0, 2).reinit(n_v, n_p);

        csp.block(1, 0).reinit(n_u, n_v);
        csp.block(1, 1).reinit(n_u, n_u);
        csp.block(1, 2).reinit(n_u, n_p);

        csp.block(2, 0).reinit(n_p, n_v);
        csp.block(2, 1).reinit(n_p, n_u);
        csp.block(2, 2).reinit(n_p, n_p);

        csp.collect_sizes();

        DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);

        sparsity_pattern.copy_from(csp);
    }

    system_matrix.reinit(sparsity_pattern);

    // Actual solution at time step n
    solution.reinit(3);
    solution.block(0).reinit(n_v);
    solution.block(1).reinit(n_u);
    solution.block(2).reinit(n_p);

    solution.collect_sizes();

    // Old timestep solution at time step n-1
    old_timestep_solution.reinit(3);
    old_timestep_solution.block(0).reinit(n_v);
    old_timestep_solution.block(1).reinit(n_u);
    old_timestep_solution.block(2).reinit(n_p);

    old_timestep_solution.collect_sizes();

    // Old Old timestep solution at time step n-2
    old_old_timestep_solution.reinit(3);
    old_old_timestep_solution.block(0).reinit(n_v);
    old_old_timestep_solution.block(1).reinit(n_u);
    old_old_timestep_solution.block(2).reinit(n_p);

    old_old_timestep_solution.collect_sizes();

    // Updates for Newton's method
    newton_update.reinit(3);
    newton_update.block(0).reinit(n_v);
    newton_update.block(1).reinit(n_u);
    newton_update.block(2).reinit(n_p);

    newton_update.collect_sizes();

    // Residual for  Newton's method
    system_rhs.reinit(3);
    system_rhs.block(0).reinit(n_v);
    system_rhs.block(1).reinit(n_u);
    system_rhs.block(2).reinit(n_p);

    system_rhs.collect_sizes();

    timer.leave_subsection();
}

#endif // SETUP_SYSTEM_CC
