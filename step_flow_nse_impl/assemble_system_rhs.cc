// assemble_system_rhs.cc

#ifndef ASSEMBLE_SYSTEM_RHS_CC
#define ASSEMBLE_SYSTEM_RHS_CC

#include "../step_flow_nse.h"

// In this function we assemble the semi-linear
// of the right hand side of Newton's method (its residual).
// The framework is in principal the same as for the
// system matrix.
template <int dim>
void StepFlowNSE<dim>::assemble_system_rhs()
{
    timer.enter_subsection("Assemble Rhs.");
    system_rhs = 0;

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

    Vector<double> local_rhs(dofs_per_cell);

    std::vector<unsigned int> local_dof_indices(dofs_per_cell);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Vector displacements(dim);
    const FEValuesExtractors::Scalar pressure(dim + dim);

    std::vector<Vector<double>>
        old_solution_values(n_q_points, Vector<double>(dim + dim + 1));

    std::vector<std::vector<Tensor<1, dim>>>
        old_solution_grads(n_q_points, std::vector<Tensor<1, dim>>(dim + dim + 1));

    std::vector<Vector<double>>
        old_solution_face_values(n_face_q_points, Vector<double>(dim + dim + 1));

    std::vector<std::vector<Tensor<1, dim>>>
        old_solution_face_grads(n_face_q_points, std::vector<Tensor<1, dim>>(dim + dim + 1));

    std::vector<Vector<double>>
        old_timestep_solution_values(n_q_points, Vector<double>(dim + dim + 1));

    std::vector<std::vector<Tensor<1, dim>>>
        old_timestep_solution_grads(n_q_points, std::vector<Tensor<1, dim>>(dim + dim + 1));

    std::vector<Vector<double>>
        old_timestep_solution_face_values(n_face_q_points, Vector<double>(dim + dim + 1));

    std::vector<std::vector<Tensor<1, dim>>>
        old_timestep_solution_face_grads(n_face_q_points, std::vector<Tensor<1, dim>>(dim + dim + 1));

    std::vector<Vector<double>>
        old_old_timestep_solution_values(n_q_points, Vector<double>(dim + dim + 1));

    std::vector<std::vector<Tensor<1, dim>>>
        old_old_timestep_solution_grads(n_q_points, std::vector<Tensor<1, dim>>(dim + dim + 1));

    typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();

    unsigned int cell_counter = 0;
    for (; cell != endc; ++cell)
    {

        fe_values.reinit(cell);
        local_rhs = 0;

        // old Newton iteration
        fe_values.get_function_values(solution, old_solution_values);
        fe_values.get_function_gradients(solution, old_solution_grads);

        // old timestep iteration
        fe_values.get_function_values(old_timestep_solution, old_timestep_solution_values);
        fe_values.get_function_gradients(old_timestep_solution, old_timestep_solution_grads);

        // old old timestep iteration
        fe_values.get_function_values(old_old_timestep_solution, old_old_timestep_solution_values);
        fe_values.get_function_gradients(old_old_timestep_solution, old_old_timestep_solution_grads);

        // The following equations are described in the subdomain with material id = 0 (currently all the domain;
        // but interesting for later, what Insa wants to do)
        if (cell->material_id() == 0)
        {
            for (unsigned int q = 0; q < n_q_points; ++q)
            {
                const Tensor<2, dim> pI = Tensors ::get_pI<dim>(q, old_solution_values);

                const Tensor<1, dim> v = Tensors ::get_v<dim>(q, old_solution_values);

                const Tensor<2, dim> grad_v = Tensors ::get_grad_v<dim>(q, old_solution_grads);

                const Tensor<1, dim> u = Tensors ::get_u<dim>(q, old_solution_values);

                const Tensor<2, dim> grad_u = Tensors ::get_grad_u<dim>(q, old_solution_grads);

                // Stress tensor fluid
                Tensor<2, dim> sigma_fluid;
                sigma_fluid.clear();
                sigma_fluid = density_fluid * kinematic_viscosity * (grad_v + transpose(grad_v));

                const Tensor<1, dim> old_timestep_u = Tensors ::get_u<dim>(q, old_timestep_solution_values);

                const Tensor<1, dim> old_old_timestep_u = Tensors ::get_u<dim>(q, old_old_timestep_solution_values);

                // Convection term of the fluid
                Tensor<1, dim> convection_fluid;
                convection_fluid.clear();
                convection_fluid = density_fluid * grad_v * v;

                // Divergence of the fluid
                const double incompressiblity_fluid = grad_v[0][0] + grad_v[1][1];

                Tensor<1, dim> fluid_force;
                fluid_force.clear();
                fluid_force[0] = force_fluid_x;
                fluid_force[1] = force_fluid_y;

                // Old time step values
                const Tensor<1, dim> old_timestep_v = Tensors ::get_v<dim>(q, old_timestep_solution_values);

                const Tensor<2, dim> old_timestep_grad_v = Tensors ::get_grad_v<dim>(q, old_timestep_solution_grads);

                const Tensor<2, dim> old_timestep_grad_u = Tensors ::get_grad_u<dim>(q, old_timestep_solution_grads);

                // Info: be careful: when force is time-dependent, then
                // this must be taken from the previous time step!!!
                Tensor<1, dim> old_timestep_fluid_force;
                old_timestep_fluid_force.clear();
                old_timestep_fluid_force[0] = force_fluid_x;
                old_timestep_fluid_force[1] = force_fluid_y;

                // TODO
                Tensor<2, dim> old_timestep_sigma_fluid;
                old_timestep_sigma_fluid.clear();
                old_timestep_sigma_fluid =
                    density_fluid * kinematic_viscosity * (old_timestep_grad_v + transpose(old_timestep_grad_v));

                Tensor<1, dim> old_timestep_convection_fluid;
                old_timestep_convection_fluid.clear();
                old_timestep_convection_fluid = density_fluid * old_timestep_grad_v * old_timestep_v;

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                    const unsigned int comp_i = fe.system_to_component_index(i).first;
                    if (comp_i == 0 || comp_i == 1)
                    {
                        // Compute velocities of Navier-Stokes flow
                        const Tensor<1, dim> phi_i_v = fe_values[velocities].value(i, q);
                        const Tensor<2, dim> phi_i_grads_v = fe_values[velocities].gradient(i, q);

                        local_rhs(i) -= (density_fluid * (v - old_timestep_v) * phi_i_v 
                                        + timestep * theta * convection_fluid * phi_i_v 
                                        + timestep * (1.0 - theta) * old_timestep_convection_fluid * phi_i_v 
                                        + timestep * scalar_product(-pI, phi_i_grads_v) 
                                        + timestep * theta * scalar_product(sigma_fluid, phi_i_grads_v) 
                                        + timestep * (1.0 - theta) * scalar_product(old_timestep_sigma_fluid, phi_i_grads_v)
                                         // Right hand side (e.g. gravitation)
                                         - timestep * theta * fluid_force * phi_i_v 
                                         - timestep * (1.0 - theta) * old_timestep_fluid_force * phi_i_v) *
                                        fe_values.JxW(q);
                    }
                    else if (comp_i == 2 || comp_i == 3)
                    {
                        // Placeholder, if some third variable at some point.
                        const Tensor<1, dim> phi_i_u = fe_values[displacements].value(i, q);
                        const Tensor<2, dim> phi_i_grads_u = fe_values[displacements].gradient(i, q);

                        local_rhs(i) -= (u * phi_i_u) * fe_values.JxW(q);
                    }
                    else if (comp_i == 4)
                    {
                        const double phi_i_p = fe_values[pressure].value(i, q);
                        // const Tensor<1,dim> phi_i_grads_t = fe_values[pressure].gradient (i, q);
                        local_rhs(i) -= (incompressiblity_fluid * phi_i_p) * fe_values.JxW(q);
                    }
                    // end i
                }
                // end n_q_points
            }

            /*
            // Pressure Neumann conditions
            for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
              {
                if (cell->face(face)->at_boundary() &&
                (cell->face(face)->boundary_id() == 0)
                )
              {

                fe_face_values.reinit (cell, face);

                fe_face_values.get_function_values (solution, old_solution_face_values);
                fe_face_values.get_function_values (old_timestep_solution, old_timestep_solution_face_values);


                for (unsigned int q=0; q<n_face_q_points; ++q)
                  {
                    Tensor<1,dim> neumann_value;
                    neumann_value[0] = 1.0;
                    neumann_value[1] = 0.0;

                    double fluid_pressure = old_timestep_solution_face_values[q](dim+dim);

                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                  {
                    const unsigned int comp_i = fe.system_to_component_index(i).first;
                    if (comp_i == 0 || comp_i == 1)
                      {
                        local_rhs(i) +=  (timestep * theta * neumann_value * fe_face_values[velocities].value (i, q)
                              ) * fe_face_values.JxW(q);
                      }
                    // end i
                  }
                    // end face_n_q_points
                  }
              }
              }  // end face integrals
            */

            cell->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(local_rhs, local_dof_indices,
                                                   system_rhs);

            // end if (for material id 0)
        }
        // Second material
        else if (cell->material_id() == 1)
        {
            std::cout << "Material id not implemented." << std::endl;
            abort();
        }

        cell_counter++;

    } // end cell

    timer.leave_subsection();
}

#endif // ASSEMBLE_SYSTEM_RHS_CC
