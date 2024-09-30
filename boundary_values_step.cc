// boundary_values.cc

#ifndef BOUNDARY_VALUES_STEP_CC
#define BOUNDARY_VALUES_STEP_CC

#include "boundary_values_step.h"
#include <cmath>
#include <iostream>

template <int dim>
double NonhomDirichletBoundaryValues<dim>::value(const Point<dim> &p,
                                                 const unsigned int component) const
{
  Assert(component < this->n_components,
         ExcIndexRange(component, 0, this->n_components));

  const long double pi = 3.141592653589793238462643;

  double inflow_velocity = 0.2;
  if (_test_case == "Channel_NSE")
    inflow_velocity = 2.0e-3;
  else if (_test_case == "2D-1_NSE")
    inflow_velocity = 0.2;
  else if (_test_case == "Cavity_IV")
    inflow_velocity = 2.0;
  else if (_test_case == "L_shaped_IV")
    inflow_velocity = 2.0;
  else if (_test_case == "Venturi_IV")
    inflow_velocity = 2.0;
  else
  {
    std::cout << "Aborting. In NonhomDirichletBoundaryValues." << std::endl;
    abort();
  }

  if (_test_case == "Channel_NSE" ||
      _test_case == "Channel_IV")
  {
    if (component == 0)
    {
      if (_time < 2.0)
      {
        return ((p(0) == 0) && (p(1) <= 0.4) ? -inflow_velocity *
                                                   (1.0 - std::cos(pi / 2.0 * _time)) / 2.0 *
                                                   1.0 *
                                                   ((p(1) - 0) * (p(1) - 0.4))
                                             : 0);
      }
      else
      {
        return ((p(0) == 0) && (p(1) <= 0.4) ? -inflow_velocity * 1.0 *
                                                   ((p(1) - 0) * (p(1) - 0.4))
                                             : 0);
      }
    }
  }
  else if (_test_case == "2D-1_NSE" ||
           _test_case == "2D-1_IV")
  {
    if (component == 0)
    {
      if (_time < 2.0)
      {
        return ((p(0) == 0) && (p(1) <= 0.41) ? -1.5 * inflow_velocity *
                                                    (1.0 - std::cos(pi / 2.0 * _time)) / 2.0 *
                                                    (4.0 / 0.1681) *
                                                    (std::pow(p(1), 2) - 0.41 * std::pow(p(1), 1))
                                              : 0);
      }
      else
      {
        return ((p(0) == 0) && (p(1) <= 0.41) ? -1.5 * inflow_velocity *
                                                    (4.0 / 0.1681) *
                                                    (std::pow(p(1), 2) - 0.41 * std::pow(p(1), 1))
                                              : 0);
      }
    }
  }
  else if (_test_case == "Cavity_IV")
  {
    if (component == 0)
    {
      if (_time < 2.0)
      {
        return ((p(0) == 0) && (p(1) <= 0.75) && (p(1) >= 0.5) ? -inflow_velocity *
                                                                     (1.0 - std::cos(pi / 2.0 * _time)) / 2.0 *
                                                                     ((p(1) - 0.5) * (p(1) - 0.75))
                                                               : 0);
      }
      else
      {
        return ((p(0) == 0) && (p(1) <= 0.75) && (p(1) >= 0.5) ? -inflow_velocity *
                                                                     ((p(1) - 0.5) * (p(1) - 0.75))
                                                               : 0);
      }
    }

    if (component == 2)
    {
      if (_time < 2.0)
      {
        return ((p(0) == 0) && (p(1) <= 1.0) && (p(1) >= 0.0) ? -0.0 * inflow_velocity *
                                                                    (1.0 - std::cos(pi / 2.0 * _time)) / 2.0 *
                                                                    ((p(1) - 0.0) * (p(1) - 1.0))
                                                              : 0);
      }
      else
      {
        return ((p(0) == 0) && (p(1) <= 1.0) && (p(1) >= 0.0) ? -0.0 * inflow_velocity *
                                                                    ((p(1) - 0.0) * (p(1) - 1.0))
                                                              : 0);
      }
    }
  }
  else if (_test_case == "L_shaped_IV")
  {
    if (component == 0)
    {
      if (_time < 2.0)
      {
        return ((p(0) == 0) && (p(1) <= 1.0) && (p(1) >= 0.75) ? -inflow_velocity *
                                                                     (1.0 - std::cos(pi / 2.0 * _time)) / 2.0 *
                                                                     ((p(1) - 0.75) * (p(1) - 1.0))
                                                               : 0);
      }
      else
      {
        return ((p(0) == 0) && (p(1) <= 1.0) && (p(1) >= 0.75) ? -inflow_velocity *
                                                                     ((p(1) - 0.75) * (p(1) - 1.0))
                                                               : 0);
      }
    }
  }
  else if (_test_case == "Venturi_IV")
  {
    if (component == 0)
    {
      if (_time < 2.0)
      {
        return ((p(0) == 0) && (p(1) <= 1.0) && (p(1) >= 0.0) ? -inflow_velocity *
                                                                    (1.0 - std::cos(pi / 2.0 * _time)) / 2.0 *
                                                                    ((p(1) - 0.0) * (p(1) - 1.0))
                                                              : 0);
      }
      else
      {
        return ((p(0) == 0) && (p(1) <= 1.0) && (p(1) >= 0.0) ? -inflow_velocity *
                                                                    ((p(1) - 0.0) * (p(1) - 1.0))
                                                              : 0);
      }
    }
  }
  else
  {
    std::cout << "Aborting. In NonhomDirichletBoundaryValues." << std::endl;
    abort();
  }

  return 0;
}

template <int dim>
void NonhomDirichletBoundaryValues<dim>::vector_value(const Point<dim> &p,
                                                      Vector<double> &values) const
{
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = NonhomDirichletBoundaryValues<dim>::value(p, c);
}

#endif // BOUNDARY_VALUES_STEP_CC
// Explicit instantiation
// template class NonhomDirichletBoundaryValues<2>;
