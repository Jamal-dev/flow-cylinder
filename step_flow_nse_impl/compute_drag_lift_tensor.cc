// COMPUTE_DRAG_LIFT_TENSOR.cc

#ifndef COMPUTE_DRAG_LIFT_TENSOR_CC
#define COMPUTE_DRAG_LIFT_TENSOR_CC

#include "../step_flow_nse.h"

// Now, we arrive at the function that is responsible
// to compute the line integrals for the drag and the lift. Note, that
// by a proper transformation via the Gauss theorem, the both
// quantities could also be achieved by domain integral computation.
// Nevertheless, we choose the line integration because deal.II provides
// all routines for face value evaluation.
template <int dim>
void StepFlowNSE<dim>::compute_drag_lift_tensor()
{

    const QGauss<dim - 1> face_quadrature_formula(3);
    FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
                                     update_values | update_gradients | update_normal_vectors |
                                         update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    std::vector<unsigned int> local_dof_indices(dofs_per_cell);
    std::vector<Vector<double>> face_solution_values(n_face_q_points,
                                                     Vector<double>(dim + dim + 1));

    std::vector<std::vector<Tensor<1, dim>>>
        face_solution_grads(n_face_q_points, std::vector<Tensor<1, dim>>(dim + dim + 1));

    Tensor<1, dim> drag_lift_value;

    typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();

    for (; cell != endc; ++cell)
    {

        // First, we are going to compute the forces that
        // act on the cylinder. We notice that only the fluid
        // equations are defined here.
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->at_boundary() &&
                cell->face(face)->boundary_id() == 80)
            {
                fe_face_values.reinit(cell, face);
                fe_face_values.get_function_values(solution, face_solution_values);
                fe_face_values.get_function_gradients(solution, face_solution_grads);

                for (unsigned int q = 0; q < n_face_q_points; ++q)
                {

                    Tensor<2, 2> pI;
                    pI.clear();
                    pI[0][0] = face_solution_values[q](dim + dim);
                    pI[1][1] = face_solution_values[q](dim + dim);

                    Tensor<2, dim> grad_v;
                    grad_v[0][0] = face_solution_grads[q][0][0];
                    grad_v[0][1] = face_solution_grads[q][0][1];
                    grad_v[1][0] = face_solution_grads[q][1][0];
                    grad_v[1][1] = face_solution_grads[q][1][1];

                    Tensor<2, dim> sigma_fluid;
                    sigma_fluid.clear();
                    sigma_fluid = -pI + density_fluid * kinematic_viscosity * (grad_v + transpose(grad_v));

                    drag_lift_value -= sigma_fluid *
                                       fe_face_values.normal_vector(q) * fe_face_values.JxW(q);
                }
            } // end boundary 80 for fluid

    } // end cell

    // 2D-1: 500; 2D-2 and 2D-3: 20 (see Schaefer/Turek 1996)
    if (test_case == "2D-1")
        drag_lift_value *= 500.0;

    std::cout << "Face drag:   " << "   " << std::setprecision(16) << drag_lift_value[0] << std::endl;
    std::cout << "Face lift:   " << "   " << std::setprecision(16) << drag_lift_value[1] << std::endl;

    double reference_value_drag = 5.5787294556197073e+00; // global ref 4
    double reference_value_lift = 1.0610686398307201e-02; // global ref 4
}

#endif // COMPUTE_DRAG_LIFT_TENSOR_CC
