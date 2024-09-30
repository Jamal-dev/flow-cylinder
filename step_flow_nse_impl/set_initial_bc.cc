// set_initial_bc.cc

#ifndef SET_INITIAL_BC_CC
#define SET_INITIAL_BC_CC

#include "../step_flow_nse.h"

// Here, we impose boundary conditions
// for the whole system.
// `Initial' means that we impose non-homoegeneous Dirichlet
//  condition on the `initial' guess of Newton's method. After,
// these must be set two homogeneous Dirichlet conditions, which
// is done in set_newton_bc ()

template <int dim>
void StepFlowNSE<dim>::set_initial_bc()
{
    // Boundary colors
    //
    //   ___3___
    //   |     |
    //  0|     |1
    //   |_____|
    //      2

    std::map<unsigned int, double> boundary_values;
    std::vector<bool> component_mask(dim + dim + 1, true);
    // 0 = vx  // test function
    // 1 = vy  // test function
    // 2 = ux  // test function
    // 3 = uy  // test function
    // 4 = pressure

    // (Scalar) pressure
    component_mask[dim + dim] = false; // false

    component_mask[0] = true;
    component_mask[1] = true;
    component_mask[dim] = true;
    component_mask[dim + 1] = true;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             // ZeroFunction<dim>(dim+dim+1),
                                             //  Use this function for nonhomogeneous
                                             //  Dirichlet conditions
                                             NonhomDirichletBoundaryValues<dim>(time, test_case),
                                             boundary_values,
                                             component_mask);

    component_mask[0] = false;
    component_mask[1] = false;
    component_mask[dim] = false;
    component_mask[dim + 1] = false;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             1,
                                             dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                             boundary_values,
                                             component_mask);

    component_mask[0] = true;
    component_mask[1] = true;
    component_mask[dim] = true;
    component_mask[dim + 1] = true;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             2,
                                             dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                             boundary_values,
                                             component_mask);

    component_mask[0] = true;
    component_mask[1] = true;
    component_mask[dim] = true;
    component_mask[dim + 1] = true;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             3,
                                             dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                             boundary_values,
                                             component_mask);
    if (test_case == "2D-1_NSE" ||
        test_case == "2D-1_IV")
    {
        component_mask[dim] = true;
        component_mask[dim + 1] = true;
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 80,
                                                 dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                                 boundary_values,
                                                 component_mask);

        VectorTools::interpolate_boundary_values(dof_handler,
                                                 81,
                                                 dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                                 boundary_values,
                                                 component_mask);
    }

    for (typename std::map<unsigned int, double>::const_iterator
             i = boundary_values.begin();
         i != boundary_values.end();
         ++i)
        solution(i->first) = i->second;
}

#endif // SET_INITIAL_BC_CC
