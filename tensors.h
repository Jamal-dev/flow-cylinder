// tensors.h

#ifndef TENSORS_H
#define TENSORS_H

#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <vector>

namespace Tensors
{
  template <int dim>
  inline Tensor<2, dim>
  get_pI(unsigned int q,
         std::vector<Vector<double>> old_solution_values)
  {
    Tensor<2, dim> tmp;
    tmp[0][0] = old_solution_values[q](dim + dim);
    tmp[1][1] = old_solution_values[q](dim + dim);

    return tmp;
  }

  template <int dim>
  inline Tensor<2, dim>
  get_pI_LinP(const double phi_i_p)
  {
    Tensor<2, dim> tmp;
    tmp.clear();
    tmp[0][0] = phi_i_p;
    tmp[1][1] = phi_i_p;

    return tmp;
  }

  // pressure
  template <int dim>
  inline Tensor<1, dim>
  get_grad_p(unsigned int q,
             std::vector<std::vector<Tensor<1, dim>>> old_solution_grads)
  {
    Tensor<1, dim> grad_p;
    grad_p[0] = old_solution_grads[q][dim + dim][0];
    grad_p[1] = old_solution_grads[q][dim + dim][1];

    return grad_p;
  }

  template <int dim>
  inline Tensor<1, dim>
  get_grad_p_LinP(const Tensor<1, dim> phi_i_grad_p)
  {
    Tensor<1, dim> grad_p;
    grad_p[0] = phi_i_grad_p[0];
    grad_p[1] = phi_i_grad_p[1];

    return grad_p;
  }

  template <int dim>
  inline Tensor<2, dim>
  get_grad_u(unsigned int q,
             std::vector<std::vector<Tensor<1, dim>>> old_solution_grads)
  {
    Tensor<2, dim> structure_continuation;
    structure_continuation[0][0] = old_solution_grads[q][dim][0];
    structure_continuation[0][1] = old_solution_grads[q][dim][1];
    structure_continuation[1][0] = old_solution_grads[q][dim + 1][0];
    structure_continuation[1][1] = old_solution_grads[q][dim + 1][1];

    return structure_continuation;
  }

  template <int dim>
  inline double
  get_divergence_u_LinU(const Tensor<2, dim> phi_i_grads_u)
  {
    return (phi_i_grads_u[0][0] + phi_i_grads_u[1][1]);
  }

  template <int dim>
  inline Tensor<2, dim>
  get_grad_v(unsigned int q,
             std::vector<std::vector<Tensor<1, dim>>> old_solution_grads)
  {
    Tensor<2, dim> grad_v;
    grad_v[0][0] = old_solution_grads[q][0][0];
    grad_v[0][1] = old_solution_grads[q][0][1];
    grad_v[1][0] = old_solution_grads[q][1][0];
    grad_v[1][1] = old_solution_grads[q][1][1];

    return grad_v;
  }

  template <int dim>
  inline Tensor<2, dim>
  get_grad_v_T(const Tensor<2, dim> tensor_grad_v)
  {
    Tensor<2, dim> grad_v_T;
    grad_v_T = transpose(tensor_grad_v);

    return grad_v_T;
  }

  template <int dim>
  inline Tensor<2, dim>
  get_grad_v_LinV(const Tensor<2, dim> phi_i_grads_v)
  {
    Tensor<2, dim> tmp;
    tmp[0][0] = phi_i_grads_v[0][0];
    tmp[0][1] = phi_i_grads_v[0][1];
    tmp[1][0] = phi_i_grads_v[1][0];
    tmp[1][1] = phi_i_grads_v[1][1];

    return tmp;
  }

  template <int dim>
  inline Tensor<2, dim>
  get_Identity()
  {
    Tensor<2, dim> identity;
    identity[0][0] = 1.0;
    identity[0][1] = 0.0;
    identity[1][0] = 0.0;
    identity[1][1] = 1.0;

    return identity;
  }

  template <int dim>
  inline Tensor<1, dim>
  get_v(unsigned int q,
        std::vector<Vector<double>> old_solution_values)
  {
    Tensor<1, dim> v;
    v[0] = old_solution_values[q](0);
    v[1] = old_solution_values[q](1);

    return v;
  }

  template <int dim>
  inline Tensor<1, dim>
  get_v_LinV(const Tensor<1, dim> phi_i_v)
  {
    Tensor<1, dim> tmp;
    tmp[0] = phi_i_v[0];
    tmp[1] = phi_i_v[1];

    return tmp;
  }

  template <int dim>
  inline Tensor<1, dim>
  get_u(unsigned int q,
        std::vector<Vector<double>> old_solution_values)
  {
    Tensor<1, dim> u;
    u[0] = old_solution_values[q](dim);
    u[1] = old_solution_values[q](dim + 1);

    return u;
  }

  template <int dim>
  inline Tensor<1, dim>
  get_u_LinU(const Tensor<1, dim> phi_i_u)
  {
    Tensor<1, dim> tmp;
    tmp[0] = phi_i_u[0];
    tmp[1] = phi_i_u[1];

    return tmp;
  }

} // end namespace tensors

#endif // TENSORS_H
