// set_newton_bc.cc

#ifndef SET_NEWTON_BC_CC
#define SET_NEWTON_BC_CC

#include "../step_flow_nse.h"

template <int dim>
void StepFlowNSE<dim>::set_newton_bc()
{
    // Boundary colors
    //
    //   ___3___
    //   |     |
    //  0|     |1
    //   |_____|
    //      2

    std::vector<bool> component_mask(dim + dim + 1, true);
    component_mask[dim + dim] = false; // temp

    component_mask[0] = true;
    component_mask[1] = true;
    component_mask[dim] = true;
    component_mask[dim + 1] = true;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                             constraints,
                                             component_mask);

    component_mask[0] = false;
    component_mask[1] = false;
    component_mask[dim] = false;
    component_mask[dim + 1] = false;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             1,
                                             dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                             constraints,
                                             component_mask);

    component_mask[0] = true;
    component_mask[1] = true;
    component_mask[dim] = true;
    component_mask[dim + 1] = true;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             2,
                                             dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                             constraints,
                                             component_mask);

    component_mask[0] = true;
    component_mask[1] = true;
    component_mask[dim] = true;
    component_mask[dim + 1] = true;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             3,
                                             dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                             constraints,
                                             component_mask);
    if (test_case == "2D-1_NSE" ||
        test_case == "2D-1_IV")
    {
        component_mask[dim] = true;
        component_mask[dim + 1] = true;
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 80,
                                                 dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                                 constraints,
                                                 component_mask);

        VectorTools::interpolate_boundary_values(dof_handler,
                                                 81,
                                                 dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                                 constraints,
                                                 component_mask);
    }
}

#endif // SET_NEWTON_BC_CC
