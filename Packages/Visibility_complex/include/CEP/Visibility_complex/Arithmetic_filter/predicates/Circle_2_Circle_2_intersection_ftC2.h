// ======================================================================
//
// Copyright (c) 1999,2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : include/CGAL/Arithmetic_filter/
// package       : Interval_arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
// ======================================================================

// This file is automatically generated by
// scripts/filtered_predicates_generator.pl

#ifndef CGAL_ARITHMETIC_FILTER_CIRCLE_2_CIRCLE_2_INTERSECTION_FTC2_H
#define CGAL_ARITHMETIC_FILTER_CIRCLE_2_CIRCLE_2_INTERSECTION_FTC2_H

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_MATCHING_BUG_2
template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
#else
static
#endif
/*  */
bool
circle_2_do_intersectC2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, CGAL_IA_PROTECTED, CGAL_IA_CACHE> &c1x,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, CGAL_IA_PROTECTED, CGAL_IA_CACHE> &c1y,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, CGAL_IA_PROTECTED, CGAL_IA_CACHE> &R1_square,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, CGAL_IA_PROTECTED, CGAL_IA_CACHE> &c2x,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, CGAL_IA_PROTECTED, CGAL_IA_CACHE> &c2y,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, CGAL_IA_PROTECTED, CGAL_IA_CACHE> &R2_square)
{
  try
  {
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return circle_2_do_intersectC2(
		c1x.interval(),
		c1y.interval(),
		R1_square.interval(),
		c2x.interval(),
		c2y.interval(),
		R2_square.interval());
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return circle_2_do_intersectC2(
		c1x.exact(),
		c1y.exact(),
		R1_square.exact(),
		c2x.exact(),
		c2y.exact(),
		R2_square.exact());
  }
}

#ifdef CGAL_IA_NEW_FILTERS

struct Static_Filtered_circle_2_do_intersectC2_6
{
  static double _bound;
  //static double ;
  static unsigned number_of_failures; // ?
  static unsigned number_of_updates;

  static bool update_epsilon(
	const Static_filter_error &c1x,
	const Static_filter_error &c1y,
	const Static_filter_error &R1_square,
	const Static_filter_error &c2x,
	const Static_filter_error &c2y,
	const Static_filter_error &R2_square)
  {
    typedef Static_filter_error FT;
  
      FT square_sum(R1_square + R2_square);
      FT a(c2x - c1x);
      FT b(c2y - c1y);
      FT d_square(a*a + b*b);
  
      if (d_square <= square_sum) return true;
  
      FT x = d_square - square_sum;
      return (x*x <= FT(4) * R1_square * R2_square);
  }

  // Call this function from the outside to update the context.
  static void new_bound (const double b) // , const double error = 0)
  {
    _bound = b;
    number_of_updates++;
    // recompute the epsilons: "just" call it over Static_filter_error.
    // That's the tricky part that might not work for everything.
    (void) update_epsilon(b,b,b,b,b,b);
    // TODO: We should verify that all epsilons have really been updated.
  }

  static bool epsilon_variant(
	const Restricted_double &c1x,
	const Restricted_double &c1y,
	const Restricted_double &R1_square,
	const Restricted_double &c2x,
	const Restricted_double &c2y,
	const Restricted_double &R2_square)
  {
    typedef Restricted_double FT;
  
      FT square_sum(R1_square + R2_square);
      FT a(c2x - c1x);
      FT b(c2y - c1y);
      FT d_square(a*a + b*b);
  
      if (d_square <= square_sum) return true;
  
      FT x = d_square - square_sum;
      return (x*x <= FT(4) * R1_square * R2_square);
  }
};

#ifndef CGAL_CFG_MATCHING_BUG_2
template < class CGAL_IA_CT, class CGAL_IA_ET, class CGAL_IA_CACHE >
#else
static
#endif
/*  */
bool
circle_2_do_intersectC2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, true, CGAL_IA_CACHE> &c1x,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, true, CGAL_IA_CACHE> &c1y,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, true, CGAL_IA_CACHE> &R1_square,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, true, CGAL_IA_CACHE> &c2x,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, true, CGAL_IA_CACHE> &c2y,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, true, CGAL_IA_CACHE> &R2_square)
{
//   bool re_adjusted = false;
  const double SAF_bound = Static_Filtered_circle_2_do_intersectC2_6::_bound;

  // Check the bounds.  All arguments must be <= SAF_bound.
  if (
	fabs(c1x.to_double()) > SAF_bound ||
	fabs(c1y.to_double()) > SAF_bound ||
	fabs(R1_square.to_double()) > SAF_bound ||
	fabs(c2x.to_double()) > SAF_bound ||
	fabs(c2y.to_double()) > SAF_bound ||
	fabs(R2_square.to_double()) > SAF_bound)
  {
// re_adjust:
    // Compute the new bound.
    double NEW_bound = 0.0;
    NEW_bound = max(NEW_bound, fabs(c1x.to_double()));
    NEW_bound = max(NEW_bound, fabs(c1y.to_double()));
    NEW_bound = max(NEW_bound, fabs(R1_square.to_double()));
    NEW_bound = max(NEW_bound, fabs(c2x.to_double()));
    NEW_bound = max(NEW_bound, fabs(c2y.to_double()));
    NEW_bound = max(NEW_bound, fabs(R2_square.to_double()));
    // Re-adjust the context.
    Static_Filtered_circle_2_do_intersectC2_6::new_bound(NEW_bound);
  }

  try
  {
    return Static_Filtered_circle_2_do_intersectC2_6::epsilon_variant(
		c1x.dbl(),
		c1y.dbl(),
		R1_square.dbl(),
		c2x.dbl(),
		c2y.dbl(),
		R2_square.dbl());
  }
  catch (...)
  {
    // if (!re_adjusted) {  // It failed, we re-adjust once.
      // re_adjusted = true;
      // goto re_adjust;
    // }
    Static_Filtered_circle_2_do_intersectC2_6::number_of_failures++;
    return circle_2_do_intersectC2(
		c1x.exact(),
		c1y.exact(),
		R1_square.exact(),
		c2x.exact(),
		c2y.exact(),
		R2_square.exact());
  }
}

#ifndef CGAL_CFG_MATCHING_BUG_2
template < class CGAL_IA_CT, class CGAL_IA_ET, class CGAL_IA_CACHE >
#else
static
#endif
/*  */
bool
circle_2_do_intersectC2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, false, CGAL_IA_CACHE> &c1x,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, false, CGAL_IA_CACHE> &c1y,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, false, CGAL_IA_CACHE> &R1_square,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, false, CGAL_IA_CACHE> &c2x,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, false, CGAL_IA_CACHE> &c2y,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, false, CGAL_IA_CACHE> &R2_square)
{
  CGAL_assertion_code(
    const double SAF_bound = Static_Filtered_circle_2_do_intersectC2_6::_bound; )
  CGAL_assertion(!(
	fabs(c1x.to_double()) > SAF_bound ||
	fabs(c1y.to_double()) > SAF_bound ||
	fabs(R1_square.to_double()) > SAF_bound ||
	fabs(c2x.to_double()) > SAF_bound ||
	fabs(c2y.to_double()) > SAF_bound ||
	fabs(R2_square.to_double()) > SAF_bound));

  try
  {
    return Static_Filtered_circle_2_do_intersectC2_6::epsilon_variant(
		c1x.dbl(),
		c1y.dbl(),
		R1_square.dbl(),
		c2x.dbl(),
		c2y.dbl(),
		R2_square.dbl());
  }
  catch (...)
  {
    Static_Filtered_circle_2_do_intersectC2_6::number_of_failures++;
    return circle_2_do_intersectC2(
		c1x.exact(),
		c1y.exact(),
		R1_square.exact(),
		c2x.exact(),
		c2y.exact(),
		R2_square.exact());
  }
}

#endif // CGAL_IA_NEW_FILTERS

CGAL_END_NAMESPACE

#endif // CGAL_ARITHMETIC_FILTER_CIRCLE_2_CIRCLE_2_INTERSECTION_FTC2_H
