namespace CGAL {

/*!
  \ingroup PkgArrangement2Enums
  The enumeration `Arr_curve_end` is used to indicate one of the two ends
  of an \f$ x\f$-monotone curve. It is used by models of the
  `ArrangementOpenBoundaryTraits_2` concept.

  \sa `ArrangementOpenBoundaryTraits_2`
*/
enum Arr_curve_end { ARR_MIN_END, ARR_MAX_END };

/*!
  \ingroup PkgArrangement2Enums

  The enumeration `Arr_halfedge_direction` is defined by
  `CGAL::Arrangement_2::Halfedge` to specify
  the direction of the halfedge.
  
  \sa `CGAL::Arrangement_2::Halfedge`
*/
enum Arr_halfedge_direction { ARR_LEFT_TO_RIGHT, ARR_RIGHT_TO_LEFT };

/*!
 Parameter space type
  \ingroup PkgArrangement2Enums
*/
typedef Box_parameter_space_2 Arr_parameter_space;

/*!
  A parameter space const value
  \ingroup PkgArrangement2Enums
*/
const Arr_parameter_space ARR_LEFT_BOUNDARY = LEFT_BOUNDARY;

/*!
  A parameter space const value
  \ingroup PkgArrangement2Enums
*/
const Arr_parameter_space ARR_RIGHT_BOUNDARY = RIGHT_BOUNDARY;

/*!
  A parameter space const value
  \ingroup PkgArrangement2Enums
*/
const Arr_parameter_space ARR_BOTTOM_BOUNDARY = BOTTOM_BOUNDARY;

/*!
  A parameter space const value
  \ingroup PkgArrangement2Enums
*/
const Arr_parameter_space ARR_TOP_BOUNDARY = TOP_BOUNDARY;

/*!
  A parameter space const value
  \ingroup PkgArrangement2Enums
*/
  const Arr_parameter_space ARR_INTERIOR = INTERIOR;


} /* end namespace CGAL */
