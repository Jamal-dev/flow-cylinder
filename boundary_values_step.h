// boundary_values.h

#ifndef BOUNDARY_VALUES_STEP_H
#define BOUNDARY_VALUES_STEP_H

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <vector>

using namespace dealii;

template <int dim>
class NonhomDirichletBoundaryValues : public Function<dim>
{
public:
  NonhomDirichletBoundaryValues(const double time,
                                const std::string test_case)
      : Function<dim>(dim + dim + 1), _time(time), _test_case(test_case)
  {
  }

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const override;

  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &value) const override;

private:
  double _time;
  std::string _test_case;
};
#include "boundary_values_step.cc"

#endif // BOUNDARY_VALUES_STEP_H
