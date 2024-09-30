// COMPUTE_FUNCTIONAL_VALUES_.cc

#ifndef COMPUTE_FUNCTIONAL_VALUES_CC
#define COMPUTE_FUNCTIONAL_VALUES_CC

#include "../step_flow_nse.h"

template <int dim>
void StepFlowNSE<dim>::compute_functional_values()
{

    double p_front, p_back, p_diff;

    p_front = compute_point_value(Point<dim>(0.15, 0.2), dim + dim); // pressure
    p_back = compute_point_value(Point<dim>(0.25, 0.2), dim + dim);  // pressure

    p_diff = p_front - p_back;

    double L2_error_velo, L2_error_dot_u, L2_error_velo_tot;
    double H1_error_velo, H1_error_dot_u, H1_error_velo_tot;

    Vector<float> error_L2_velo(triangulation.n_active_cells());
    Vector<float> error_L2_dot_u(triangulation.n_active_cells());
    Vector<float> error_L2_velo_tot(triangulation.n_active_cells());

    Vector<float> error_H1_velo(triangulation.n_active_cells());
    Vector<float> error_H1_dot_u(triangulation.n_active_cells());
    Vector<float> error_H1_velo_tot(triangulation.n_active_cells());

    BlockVector<double> solution_dot_u;
    solution_dot_u = solution;

    // dot u
    for (unsigned int j = 0; j < solution.block(0).size(); ++j)
    {
        solution_dot_u.block(0)(j) = (solution.block(1)(j) - old_timestep_solution.block(1)(j)) / timestep;
    }

    // Select first dim components
    ComponentSelectFunction<dim> value_select(std::pair(0, dim), dim + dim + 1);

    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                      error_L2_velo,
                                      QGauss<dim>(degree + 2),
                                      VectorTools::L2_norm,
                                      &value_select);

    L2_error_velo = VectorTools::compute_global_error(triangulation,
                                                      error_L2_velo,
                                                      VectorTools::L2_norm);

    VectorTools::integrate_difference(dof_handler,
                                      solution_dot_u,
                                      dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                      error_L2_dot_u,
                                      QGauss<dim>(degree + 2),
                                      VectorTools::L2_norm,
                                      &value_select);

    L2_error_dot_u = VectorTools::compute_global_error(triangulation,
                                                       error_L2_dot_u,
                                                       VectorTools::L2_norm);

    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                      error_H1_velo,
                                      QGauss<dim>(degree + 2),
                                      VectorTools::H1_norm,
                                      &value_select);

    H1_error_velo = VectorTools::compute_global_error(triangulation,
                                                      error_H1_velo,
                                                      VectorTools::H1_norm);

    VectorTools::integrate_difference(dof_handler,
                                      solution_dot_u,
                                      dealii::Functions::ZeroFunction<dim>(dim + dim + 1),
                                      error_H1_dot_u,
                                      QGauss<dim>(degree + 2),
                                      VectorTools::H1_norm,
                                      &value_select);
    H1_error_dot_u = VectorTools::compute_global_error(triangulation,
                                                       error_H1_dot_u,
                                                       VectorTools::H1_norm);

    std::cout << "------------------" << std::endl;
    std::cout << "P-Diff:  " << "   " << std::setprecision(16) << p_diff << std::endl;
    std::cout << "P-front: " << "   " << std::setprecision(16) << p_front << std::endl;
    std::cout << "P-back:  " << "   " << std::setprecision(16) << p_back << std::endl;

    std::cout << "------------------" << std::endl;
    // Compute drag and lift via line integral
    compute_drag_lift_tensor();

    std::cout << "------------------" << std::endl;
    // std::cout << "L2_error_velo:     " << time << "   " << L2_error_velo << std::endl;
    // std::cout << "L2_error_dot_u:    " << time << "   " << L2_error_dot_u << std::endl;
    // std::cout << "H1_error_velo:     " << time << "   " << H1_error_velo << std::endl;
    // std::cout << "H1_error_dot_u:    " << time << "   " << H1_error_dot_u << std::endl;

    std::cout << "------------------" << std::endl;
    std::cout << std::endl;
}

#endif // COMPUTE_FUNCTIONAL_VALUES_CC
