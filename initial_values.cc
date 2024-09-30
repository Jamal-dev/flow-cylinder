// initial_values.cc

#ifndef INITIAL_VALUES_CC
#define INITIAL_VALUES_CC

#include "initial_values.h"

template <int dim>
double InitialValues<dim>::value(const Point<dim> & /*p*/,
                                 const unsigned int component) const
{
  // Only pressure
  if (component == 4)
  {
    return 0.0; // Initial pressure if desired
  }

  return 0.0;
}

template <int dim>
void InitialValues<dim>::vector_value(const Point<dim> &p,
                                      Vector<double> &values) const
{
  for (unsigned int comp = 0; comp < this->n_components; ++comp)
    values(comp) = InitialValues<dim>::value(p, comp);
}

// Explicit instantiation
// template class InitialValues<2>;

#endif // INITIAL_VALUES_CC
