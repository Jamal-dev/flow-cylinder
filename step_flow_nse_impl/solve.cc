// solve.cc

#ifndef SOLVE_CC
#define SOLVE_CC

#include "../step_flow_nse.h"

// In this function, we solve the linear systems
// inside the nonlinear Newton iteration. We just
// use a direct solver from UMFPACK.
template <int dim>
void StepFlowNSE<dim>::solve()
{
    timer.enter_subsection("Solve linear system.");
    Vector<double> sol, rhs;
    sol = newton_update;
    rhs = system_rhs;

    // SparseDirectUMFPACK A_direct;
    // A_direct.factorize(system_matrix);
    A_direct.vmult(sol, rhs);
    newton_update = sol;

    constraints.distribute(newton_update);
    timer.leave_subsection();
}

#endif // SOLVE_CC
