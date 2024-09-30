// assemble_system_matrix.cc

#ifndef ASSEMBLE_SYSTEM_MATRIX_CC
#define ASSEMBLE_SYSTEM_MATRIX_CC

#include "../step_flow_nse.h"

// In this function, we assemble the Jacobian matrix
// for the Newton iteration.
//
// Assembling of the inner most loop is treated with help of
// the fe.system_to_component_index(j).first function from
// the library.
// Using this function makes the assembling process much faster
// than running over all local degrees of freedom.
template <int dim>
void StepFlowNSE<dim>::assemble_system_matrix()
{
    timer.enter_subsection("Assemble Matrix.");
    system_matrix = 0;

    QGauss<dim> quadrature_formula(degree + 2);
    QGauss<dim - 1> face_quadrature_formula(degree + 2);

    FEValues<dim> fe_values(fe, quadrature_formula,
                            update_values |
                                update_quadrature_points |
                                update_JxW_values |
                                update_gradients);

    FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
                                     update_values | update_quadrature_points |
                                         update_normal_vectors | update_gradients |
                                         update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;

    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices(dofs_per_cell);

    // Now, we are going to use the
    // FEValuesExtractors to determine
    // the four principle variables
    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Vector displacements(dim);  // 2
    const FEValuesExtractors::Scalar pressure(dim + dim); // 4

    // We declare Vectors and Tensors for
    // the solutions at the previous Newton iteration:
    std::vector<Vector<double>> old_solution_values(n_q_points,
                                                    Vector<double>(dim + dim + 1));

    std::vector<std::vector<Tensor<1, dim>>> old_solution_grads(n_q_points,
                                                                std::vector<Tensor<1, dim>>(dim + dim + 1));

    std::vector<Vector<double>> old_solution_face_values(n_face_q_points,
                                                         Vector<double>(dim + dim + dim + 1));

    std::vector<std::vector<Tensor<1, dim>>> old_solution_face_grads(n_face_q_points,
                                                                     std::vector<Tensor<1, dim>>(dim + dim + 1));

    // We declare Vectors and Tensors for
    // the solution at the previous time step:
    std::vector<Vector<double>> old_timestep_solution_values(n_q_points,
                                                             Vector<double>(dim + dim + 1));

    std::vector<std::vector<Tensor<1, dim>>> old_timestep_solution_grads(n_q_points,
                                                                         std::vector<Tensor<1, dim>>(dim + dim + 1));

    std::vector<Vector<double>> old_timestep_solution_face_values(n_face_q_points,
                                                                  Vector<double>(dim + dim + 1));

    std::vector<std::vector<Tensor<1, dim>>> old_timestep_solution_face_grads(n_face_q_points,
                                                                              std::vector<Tensor<1, dim>>(dim + dim + 1));

    // Declaring test functions:
    std::vector<Tensor<1, dim>> phi_i_v(dofs_per_cell);
    std::vector<Tensor<2, dim>> phi_i_grads_v(dofs_per_cell);
    std::vector<double> phi_i_p(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_i_grads_p(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_i_u(dofs_per_cell);
    std::vector<Tensor<2, dim>> phi_i_grads_u(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();

    unsigned int cell_counter = 0;
    for (; cell != endc; ++cell)
    {
        fe_values.reinit(cell);
        local_matrix = 0;

        // Old Newton iteration values
        fe_values.get_function_values(solution, old_solution_values);
        fe_values.get_function_gradients(solution, old_solution_grads);

        // Old_timestep_solution values
        fe_values.get_function_values(old_timestep_solution, old_timestep_solution_values);
        fe_values.get_function_gradients(old_timestep_solution, old_timestep_solution_grads);

        // Material id = 0 for first material; later Insa wants a second material with different coefficient,
        // for which we then need id = 1 as well; here as well as in the *.inp mesh file.
        if (cell->material_id() == 0)
        {
            for (unsigned int q = 0; q < n_q_points; ++q)
            {
                for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                    phi_i_v[k] = fe_values[velocities].value(k, q);
                    phi_i_grads_v[k] = fe_values[velocities].gradient(k, q);
                    phi_i_p[k] = fe_values[pressure].value(k, q);
                    phi_i_grads_p[k] = fe_values[pressure].gradient(k, q);
                    phi_i_u[k] = fe_values[displacements].value(k, q);
                    phi_i_grads_u[k] = fe_values[displacements].gradient(k, q);
                }

                const Tensor<1, dim> v = Tensors ::get_v<dim>(q, old_solution_values);

                const Tensor<1, dim> u = Tensors ::get_u<dim>(q, old_solution_values);

                const Tensor<1, dim> old_timestep_u = Tensors ::get_u<dim>(q, old_timestep_solution_values);

                const Tensor<2, dim> grad_v = Tensors ::get_grad_v<dim>(q, old_solution_grads);

                const Tensor<2, dim> grad_u = Tensors ::get_grad_u<dim>(q, old_solution_grads);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {

                    Tensor<2, 2> pI_LinP;
                    pI_LinP.clear();
                    pI_LinP[0][0] = phi_i_p[i];
                    pI_LinP[1][1] = phi_i_p[i];

                    // Stress tensor flow (Newtonian fluid)
                    const Tensor<2, dim> stress_fluid_LinAll =
                        density_fluid * kinematic_viscosity * (phi_i_grads_v[i] + transpose(phi_i_grads_v[i]));

                    const Tensor<1, dim> convection_fluid_LinAll = density_fluid * (phi_i_grads_v[i] * v + grad_v * phi_i_v[i]);

                    const double incompressibility_LinAll = phi_i_grads_v[i][0][0] + phi_i_grads_v[i][1][1];

                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {

                        const unsigned int comp_j = fe.system_to_component_index(j).first;
                        if (comp_j == 0 || comp_j == 1)
                        {
                            // Compute velocities v of Navier-Stokes flow
                            local_matrix(j, i) += (density_fluid * phi_i_v[i] * phi_i_v[j] 
                                        + timestep * theta * convection_fluid_LinAll * phi_i_v[j]
                                                   // Pressure stress
                                                   + timestep * scalar_product(-pI_LinP, phi_i_grads_v[j])
                                                   // Velocities stress
                                                   + timestep * theta * scalar_product(stress_fluid_LinAll, phi_i_grads_v[j])) *
                                                  fe_values.JxW(q);
                        }
                        else if (comp_j == 2 || comp_j == 3)
                        {
                            // Placeholder for some third equation if necessary at some point.
                            local_matrix(j, i) += (phi_i_u[i] * phi_i_u[j]) * fe_values.JxW(q);
                        }
                        else if (comp_j == 4)
                        {
                            // Pressure equation: Newton linearization
                            local_matrix(j, i) += (incompressibility_LinAll * phi_i_p[j]) * fe_values.JxW(q);
                        }
                        // end j dofs
                    }
                    // end i dofs
                }
                // end n_q_points
            }

            cell->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(local_matrix, local_dof_indices,
                                                   system_matrix);
        }
        // Second material id
        else if (cell->material_id() == 1)
        {
            // Currently not implemented in this code.
            std::cout << "Material id not implemented." << std::endl;
            abort();
        }

        cell_counter++;

    } // end cell

    timer.leave_subsection();
}

#endif // ASSEMBLE_SYSTEM_MATRIX_CC
