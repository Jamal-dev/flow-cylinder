// initial_values.h

#ifndef INITIAL_VALUES_H
#define INITIAL_VALUES_H

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <vector>

using namespace dealii;

template <int dim>
class InitialValues : public Function<dim>
{
public:
  InitialValues() : Function<dim>(dim + dim + 1) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const override;

  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &value) const override;
};
#include "initial_values.cc"
#endif // INITIAL_VALUES_H
