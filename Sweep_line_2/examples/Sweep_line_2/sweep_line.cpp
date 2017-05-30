#define CGAL_SL_VERBOSE

//! \file examples/Arrangement_on_surface_2/sweep_line.cpp
// Computing intersection points among curves using the sweep line.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <CGAL/Arr_polycurve_basic_traits_2.h>
#include <list>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT NT;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;

//~ typedef CGAL::Arr_segment_traits_2<Kernel>              Base_traits_2;
//~ typedef CGAL::Arr_polycurve_basic_traits_2<Base_traits_2> Traits;

typedef Traits_2::Curve_2                               Segment_2;

//~ 4
//~ 4 
//~ 5 0
//~ 4 1
//~ 3 1
//~ 2 1

//~ 2
//~ 0 0
//~ 5 0

//~ 5
//~ 5 0
//~ 4 1
//~ 3 1
//~ 2 1
//~ 0 0

//~ 3
//~ 1 2
//~ 3 1
//~ 4 1


int main()
{
  // Construct the input segments.
  const int nbs=7;
  Segment_2 segments[] = {Segment_2 (Point_2 (1, 5), Point_2 (8, 5)),
                          Segment_2 (Point_2 (1, 1), Point_2 (8, 8)),
                          Segment_2 (Point_2 (3, 1), Point_2 (3, 8)),
                          Segment_2 (Point_2 (8, 5), Point_2 (8, 8)),
                          Segment_2 (Point_2 (3, 3), Point_2 (3, 8)),
                          Segment_2 (Point_2 (6, 5), Point_2 (8, 4)),
                          Segment_2 (Point_2 (2, 4), Point_2 (4, 6)),
    };
  
  // Compute all intersection points.
  std::list<Point_2>     pts;

  CGAL::compute_intersection_points (segments, segments + nbs,
                                     std::back_inserter (pts));
  
  // Print the result.
  std::cout << "Found " << pts.size() << " intersection points: " << std::endl; 
  std::copy (pts.begin(), pts.end(),
             std::ostream_iterator<Point_2>(std::cout, "\n"));

  // Compute the non-intersecting sub-segments induced by the input segments.
  std::list<Segment_2>   sub_segs;

  CGAL::compute_subcurves(segments, segments + nbs, std::back_inserter(sub_segs));

  std::cout << "Found " << sub_segs.size()
            << " interior-disjoint sub-segments." << std::endl;

  CGAL_assertion (CGAL::do_curves_intersect (segments, segments + nbs));

  return 0;
}
