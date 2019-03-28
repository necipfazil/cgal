// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_CONS_FTC2_H
#define CGAL_STRAIGHT_SKELETON_CONS_FTC2_H 1

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/predicates/Straight_skeleton_pred_ftC2.h>
#include <CGAL/Straight_skeleton_2/nt_utils.h>

namespace CGAL { 

namespace CGAL_SS_i
{

template<class K>
Uncertain<Trisegment_collinearity> certified_trisegment_collinearity ( Segment_2<K> const& e0
                                                                     , Segment_2<K> const& e1
                                                                     , Segment_2<K> const& e2
                                                                     );

// Given an oriented 2D straight line segment 'e', computes the normalized coefficients (a,b,c) of the
// supporting line.
// POSTCONDITION: [a,b] is the leftward normal _unit_ (a�+b�=1) vector.
// POSTCONDITION: In case of overflow, an empty optional<> is returned.
template<class K>
optional< Line_2<K> > compute_normalized_line_ceoffC2( Segment_2<K> const& e )
{
  bool finite = true ;
  
  typedef typename K::FT FT ;
  
  FT a (0.0),b (0.0) ,c(0.0)  ;

  if(e.source().y() == e.target().y())
  {
    a = 0 ;
    if(e.target().x() > e.source().x())
    {
      b = 1;
      c = -e.source().y();
    }
    else if(e.target().x() == e.source().x())
    {
      b = 0;
      c = 0;
    }
    else
    {
      b = -1;
      c = e.source().y();
    }

    CGAL_STSKEL_TRAITS_TRACE("Line coefficients for HORIZONTAL line:\n" 
                            << s2str(e) 
                            << "\na="<< n2str(a) << ", b=" << n2str(b) << ", c=" << n2str(c)
                           ) ;
  }
  else if(e.target().x() == e.source().x())
  {
    b = 0;
    if(e.target().y() > e.source().y())
    {
      a = -1;
      c = e.source().x();
    }
    else if (e.target().y() == e.source().y())
    {
      a = 0;
      c = 0;
    }
    else
    {
      a = 1;
      c = -e.source().x();
    }

    CGAL_STSKEL_TRAITS_TRACE("Line coefficients for VERTICAL line:\n"
                            << s2str(e) 
                            << "\na="<< n2str(a) << ", b=" << n2str(b) << ", c=" << n2str(c)
                           ) ;
  }
  else
  {
    FT sa = e.source().y() - e.target().y();
    FT sb = e.target().x() - e.source().x();
    FT l2 = (sa*sa) + (sb*sb) ;

    if ( CGAL_NTS is_finite(l2) )
    {
      FT l = CGAL_SS_i :: inexact_sqrt(l2);
  
      a = sa / l ;
      b = sb / l ;
  
      c = -e.source().x()*a - e.source().y()*b;
      
      CGAL_STSKEL_TRAITS_TRACE("Line coefficients for line:\n"
                               << s2str(e) 
                               << "\nsa="<< n2str(sa) << "\nsb=" << n2str(sb) << "\nl2=" << n2str(l2) << "\nl=" << n2str(l)
                               << "\na="<< n2str(a) << "\nb=" << n2str(b) << "\nc=" << n2str(c)
                               ) ;
    }
    else finite = false ;
    
  }
  
  if ( finite )
    if ( !CGAL_NTS is_finite(a) || !CGAL_NTS is_finite(b) || !CGAL_NTS is_finite(c) ) 
      finite = false ;

  return cgal_make_optional( finite, K().construct_line_2_object()(a,b,c) ) ;
}

template<class K>
optional< Line_2<K> > compute_line_ceoffC2( Segment_2<K> const& e )
{
  bool finite = true ;

  typedef typename K::FT FT ;

  FT a (0.0),b (0.0) ,c(0.0)  ;

  if(e.source().y() == e.target().y())
  {
    a = 0 ;
    if(e.target().x() > e.source().x())
    {
      b = 1;
      c = -e.source().y();
    }
    else if(e.target().x() == e.source().x())
    {
      b = 0;
      c = 0;
    }
    else
    {
      b = -1;
      c = e.source().y();
    }

    CGAL_STSKEL_TRAITS_TRACE("Line coefficients for HORIZONTAL line:\n"
                            << s2str(e)
                            << "\na="<< n2str(a) << ", b=" << n2str(b) << ", c=" << n2str(c)
                           ) ;
  }
  else if(e.target().x() == e.source().x())
  {
    b = 0;
    if(e.target().y() > e.source().y())
    {
      a = -1;
      c = e.source().x();
    }
    else if (e.target().y() == e.source().y())
    {
      a = 0;
      c = 0;
    }
    else
    {
      a = 1;
      c = -e.source().x();
    }

    CGAL_STSKEL_TRAITS_TRACE("Line coefficients for VERTICAL line:\n"
                            << s2str(e)
                            << "\na="<< n2str(a) << ", b=" << n2str(b) << ", c=" << n2str(c)
                           ) ;
  }
  else
  {
    a = e.source().y() - e.target().y();
    b = e.target().x() - e.source().x();
    c = -e.source().x()*a - e.source().y()*b;
  }

  if ( finite )
    if ( !CGAL_NTS is_finite(a) || !CGAL_NTS is_finite(b) || !CGAL_NTS is_finite(c) )
      finite = false ;

  return cgal_make_optional( finite, K().construct_line_2_object()(a,b,c) ) ;
}

template<class FT>
Rational_time<FT> squared_distance_from_point_to_lineC2( FT const& px, FT const& py, FT const& sx, FT const& sy, FT const& tx, FT const& ty )
{
  FT ldx = tx - sx ;
  FT ldy = ty - sy ;
  FT rdx = sx - px ;
  FT rdy = sy - py ;
  
  FT n = CGAL_NTS square(ldx * rdy - rdx * ldy);
  FT d = CGAL_NTS square(ldx) + CGAL_NTS square(ldy);
  
  return Rational_time<FT>(n,d) ;
}

//
// Constructs a Trisegment_2 which stores 3 oriented straight line segments e0,e1,e2 along with their collinearity.
//
// NOTE: If the collinearity cannot be determined reliably, a null trisegment is returned.
//
template<class K>
intrusive_ptr< Trisegment_2<K> > construct_trisegment ( Segment_2<K> const& e0
                                                      , Segment_2<K> const& e1
                                                      , Segment_2<K> const& e2
                                                      )
{
  typedef Trisegment_2<K>                 Trisegment_2 ;
  typedef typename Trisegment_2::Self_ptr Trisegment_2_ptr ;
  
  Uncertain<Trisegment_collinearity> lCollinearity = certified_trisegment_collinearity(e0,e1,e2);
  
  if (is_certain(lCollinearity) ) 
       return Trisegment_2_ptr( new Trisegment_2(e0, e1, e2, lCollinearity) ) ;
  else return Trisegment_2_ptr();
}

// Given 3 oriented straight line segments: e0, e1, e2 
// returns the OFFSET DISTANCE (n/d) at which the offsetted lines
// intersect at a single point, IFF such intersection exist.
// If the lines intersect to the left, the returned distance is positive.
// If the lines intersect to the right, the returned distance is negative.
// If the lines do not intersect, for example, for collinear edges, or parallel edges but with the same orientation,
// returns 0 (the actual distance is undefined in this case, but 0 is a usefull return)
//
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
//
// PRECONDITION: None of e0, e1 and e2 are collinear (but two of them can be parallel)
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
// NOTE: The segments (e0,e1,e2) are stored in the argument as the trisegment st.event()
//
template<class K>
optional< Rational_time< typename K::FT> > compute_normal_offset_lines_isec_timeC2 ( intrusive_ptr< Trisegment_2<K> > const& tri )
{
  typedef typename K::FT  FT ;
  
  typedef Line_2<K> Line_2 ;
  
  typedef optional<Line_2> Optional_line_2 ;
  
  CGAL_STSKEL_TRAITS_TRACE("Computing normal offset lines isec time for: " << tri ) ;
  
  // DETAILS:
  //
  // An offset line is given by:
  //
  //   a*x(t) + b*y(t) + c - t = 0
  //
  // were 't > 0' being to the left of the line.
  // If 3 such offset lines intersect at the same offset distance, the intersection 't',
  // or 'time', can be computed solving for 't' in the linear system formed by 3 such equations.
  // The result is :
  //
  //  t = a2*b0*c1 - a2*b1*c0 - b2*a0*c1 + b2*a1*c0 + b1*a0*c2 - b0*a1*c2
  //      ---------------------------------------------------------------
  //             -a2*b1 + a2*b0 + b2*a1 - b2*a0 + b1*a0 - b0*a1 ;

  Optional_line_2 l0 = compute_line_ceoffC2(tri->e0()) ;
  Optional_line_2 l1 = compute_line_ceoffC2(tri->e1()) ;
  Optional_line_2 l2 = compute_line_ceoffC2(tri->e2()) ;

  if ( l0 && l1 && l2 )
  {
    FT num = (l2->a()*l0->b()*l1->c())
             -(l2->a()*l1->b()*l0->c())
             -(l2->b()*l0->a()*l1->c())
             +(l2->b()*l1->a()*l0->c())
             +(l1->b()*l0->a()*l2->c())
             -(l0->b()*l1->a()*l2->c());

    FT sum_sq_0 = square(l0->a()) + square(l0->b());
    FT sum_sq_1 = square(l1->a()) + square(l1->b());
    FT sum_sq_2 = square(l2->a()) + square(l2->b());

    return boost::make_optional(
      Rational_time<FT>(num,
                        l2->b()*l1->a() - l2->a()*l1->b(),
                        l2->a()*l0->b() - l2->b()*l0->a(),
                        l1->b()*l0->a() - l0->b()*l1->a(),
                        sum_sq_0, sum_sq_1, sum_sq_2));
  }

  return boost::none;
}

// compute the time of the intersection of the bisector of (e0,e1) with the bisector of (e2,e3) on the bisector of (e0,e1)
// in the case the bisector are not collinear
template<class K>
optional< Rational_time_4< typename K::FT> > compute_normal_offset_lines_isec_timeC2 ( const typename K::Segment_2& e0,
                                                                                       const typename K::Segment_2& e1,
                                                                                       const typename K::Segment_2& e2,
                                                                                       const typename K::Segment_2& e3)
{
  typedef typename K::FT FT;
  typedef Line_2<K> Line_2;
  typedef optional<Line_2> Optional_line_2 ;

  Optional_line_2 l0 = compute_line_ceoffC2(e0) ;
  Optional_line_2 l1 = compute_line_ceoffC2(e1) ;
  Optional_line_2 l2 = compute_line_ceoffC2(e2) ;
  Optional_line_2 l3 = compute_line_ceoffC2(e3) ;

  if ( l0 && l1 && l2 && l3)
  {
    FT sum_sq_0 = square(l0->a()) + square(l0->b());
    FT sum_sq_1 = square(l1->a()) + square(l1->b());
    FT sum_sq_2 = square(l2->a()) + square(l2->b());
    FT sum_sq_3 = square(l3->a()) + square(l3->b());

    const FT a0b1_a1b1 = l0->a()*l1->b() - l1->a()*l0->b(),
             a0b2_a2b0 = l0->a()*l2->b() - l2->a()*l0->b(),
             a2b0_a0b2 = l2->a()*l0->b() - l0->a()*l2->b(),
             a1b2_a2b1 = l1->a()*l2->b() - l2->a()*l1->b(),
             a3b0_a0b3 = l3->a()*l0->b() - l0->a()*l3->b(),
             a1b3_a3b1 = l1->a()*l3->b() - l3->a()*l1->b(),
             a0b3_a3b0 = l0->a()*l3->b() - l3->a()*l0->b(),
             a3b1_a1b3 = l3->a()*l1->b() - l1->a()*l3->b(),
             a2b1_a0b2 = l2->a()*l1->b() - l1->a()*l2->b(),
             a2b3_a3b2 = l2->a()*l3->b() - l3->a()*l2->b(),
             d01_exp1 = (a0b1_a1b1 * l3->c() + a3b0_a0b3 * l1->c() + a1b3_a3b1 * l0->c()) * sum_sq_2,
             d01_exp2 = (a0b1_a1b1 * l2->c() + a2b0_a0b2 * l1->c() + a1b2_a2b1 * l0->c()) * sum_sq_3,
             two(2);

   const FT num = ( (a0b1_a1b1 * (a0b1_a1b1 * l2->c() + two * a2b0_a0b2 * l1->c() + two * a1b2_a2b1 * l0->c()) ) * l2->c()
                +   (a0b2_a2b0 * (a0b2_a2b0 * l1->c() + two * a2b1_a0b2 * l0->c() ) ) * l1->c()
                +    square(a1b2_a2b1 * l0->c()) ) * sum_sq_3
                + (( -a0b1_a1b1 * (a0b1_a1b1 * l3->c() + two * a3b0_a0b3 * l1->c() + two * a1b3_a3b1 * l0->c()) ) * l3->c()
                + (  -a0b3_a3b0 * (a0b3_a3b0 * l1->c() + two * a3b1_a1b3 * l0->c())) * l1->c()
                +   - square(a1b3_a3b1 * l0->c())) * sum_sq_2,
            d0 =  a1b2_a2b1 * d01_exp2 - a1b3_a3b1 * d01_exp1,
            d1 = -a0b2_a2b0 * d01_exp2 + a0b3_a3b0 * d01_exp1,
            d2 = -a0b1_a1b1 * ( a0b2_a2b0 * l3->c() + a3b0_a0b3 * l2->c() + a2b3_a3b2 * l0->c() ),
            d3 =  a0b1_a1b1 * ( a1b2_a2b1 * l3->c() + a3b1_a1b3 * l2->c() + a2b3_a3b2 * l1->c() );

    return boost::make_optional(
      Rational_time_4<FT>(num,
                         d0, d1, d2, d3,
                         sum_sq_0, sum_sq_1, sum_sq_1 * sum_sq_2 * sum_sq_3, sum_sq_0 * sum_sq_2* sum_sq_3) );
  }

  return boost::none;
}

// Given two oriented straight line segments e0 and e1 such that e-next follows e-prev, returns
// the coordinates of the midpoint of the segment between e-prev and e-next.
// NOTE: the edges can be oriented e0->e1 or e1->e0
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
optional< Point_2<K> > compute_oriented_midpoint ( Segment_2<K> const& e0, Segment_2<K> const& e1 )
{
  Point_2<K> mp = CGAL::compare_distance(e0.target(), e1.source(), e1.target(), e0.source())
                ? CGAL::midpoint(e0.target(),e1.source())
                : CGAL::midpoint(e1.target(),e0.source());

  CGAL_STSKEL_TRAITS_TRACE("Computing oriented midpoint between:"
                           << "\ne0: " << s2str(e0)
                           << "\ne1: " << s2str(e1)
                           << "\nmp=" << p2str(mp)
                          );

  bool ok = CGAL_NTS is_finite(mp.x()) && CGAL_NTS is_finite(mp.y());

  return cgal_make_optional(ok,mp);
}

//
// Given 3 oriented straight line segments: e0, e1, e2 and the corresponding offseted segments: e0*, e1* and e2*,
// returns the point of the left or right seed (offset vertex) (e0*,e1*) or (e1*,e2*) 
// 
// If the current event (defined by e0,e1,e2) is a propagated event, that is, it follows from a previous event,
// the seeds are skeleten nodes and are given by non-null trisegments.
// If the current event is an initial event the seeds are contour vertices and are given by null trisegmets.
//
// If a seed is a skeleton node, its point has to be computed from the trisegment that defines it.
// That trisegment is exactly the trisegment tree that defined the previous event which produced the skeleton node
// (so the trisegment tree is basically a lazy representation of the seed point).
//
// If a seed is a contour vertex, its point is then simply the target endoint of e0 or e1 (for the left/right seed).
//
// This method returns the specified seed point (left or right)
//
// NOTE: Split events involve 3 edges but only one seed, the left (that is, only e0*,e1* is connected before the event).
// The trisegment tree for a split event has always a null right child even if the event is not an initial event
// (in which case its left child won't be null).
// If you ask for the right child point for a trisegment tree corresponding to a split event you will just get e1.target()
// which is nonsensical for a non initial split event.
//
// NOTE: There is an abnormal collinearity case which ocurrs when e0 and e2 are collinear.
// In this case, these lines do not correspond to an offset vertex (because e0* and e2* are never consecutive before the event),
// so the degenerate seed is neither the left or the right seed. In this case, the SEED ID for the degenerate pseudo seed is UNKOWN.
// If you request the point of such degenerate pseudo seed the oriented midpoint bettwen e0 and e2 is returned.
//
template<class K>
optional< Point_2<K> > compute_seed_pointC2 ( intrusive_ptr< Trisegment_2<K> > const& tri, typename Trisegment_2<K>::SEED_ID sid )
{
  optional< Point_2<K> > p ;

  typedef Trisegment_2<K> Trisegment_2 ;
    
  switch ( sid )
  {
    case Trisegment_2::LEFT :
         
       p = tri->child_l() ? construct_offset_lines_isecC2(tri->child_l())  // this can recurse
                          : compute_oriented_midpoint(tri->e0(),tri->e1()) ;
       break ;                     
             
    case Trisegment_2::RIGHT : 
    
      p = tri->child_r() ? construct_offset_lines_isecC2(tri->child_r()) // this can recurse
                         : compute_oriented_midpoint(tri->e1(),tri->e2()) ;
      break ;        
             
    case Trisegment_2::UNKNOWN : 
    
      p = compute_oriented_midpoint(tri->e0(),tri->e2());
      
      break ;
  }
  
  return p ;  
}

//
// Given the trisegment tree for an event which is known to have a normal collinearity returns the seed point
// of the degenerate seed.
// A normal collinearity occurs when e0,e1 or e1,e2 are collinear.
template<class K>
optional< Point_2<K> > compute_degenerate_seed_pointC2 ( intrusive_ptr< Trisegment_2<K> > const& tri )
{
  return compute_seed_pointC2( tri, tri->degenerate_seed_id() ) ;
}

// Given 3 oriented straight line segments: e0, e1, e2 
// such that two and only two of these edges are collinear, not neccesarily consecutive but with the same orientaton;
// returns the OFFSET DISTANCE (n/d) at which a line perpendicular to the collinear edge passing through
// the degenerate seed point intersects the offset line of the non collinear edge
//
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
optional< Rational_time< typename K::FT> > compute_degenerate_offset_lines_isec_timeC2 ( intrusive_ptr< Trisegment_2<K> > const& tri )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K> Point_2 ;
  typedef Line_2 <K> Line_2 ;
  
  typedef optional<Point_2> Optional_point_2 ;
  typedef optional<Line_2>  Optional_line_2 ;
  
  CGAL_STSKEL_TRAITS_TRACE("Computing degenerate offset lines isec time for: " << tri ) ;
  
  // DETAILS:
  //
  // For simplicity, assume e0,e1 are the collinear edges.
  //
  //   (1)
  //   The bisecting line of e0 and e1 is a line perpendicular to e0 (and e1)
  //   which passes through 'q': the degenerate offset vertex (e0*,e1*)
  //   This "degenerate" bisecting line is given by:
  //
  //     B0(t) = p + t*[l0.a,l0.b]
  //
  //   where p is the projection of q along l0 and l0.a,l0.b are the _normalized_ line coefficients for e0 (or e1 which is the same)
  //   Since [a,b] is a _unit_ vector pointing perpendicularly to the left of e0 (and e1);
  //   any point B0(k) is at a distance k from the line supporting e0 and e1.
  //
  //   (2)
  //   The bisecting line of e0 and e2 is given by the following SEL
  //
  //    l0.a*x(t) + l0.b*y(t) + l0.c + t = 0
  //    l2.a*x(t) + l2.b*y(t) + l2.c + t = 0
  //
  //   where (l0.a,l0.b,l0.c) and (l2.a,l2.b,l0.c) are the normalized line coefficientes of e0 and e2 resp.
  //
  //     B1(t)=[x(t),y(t)]
  //
  //   (3)
  //   These two bisecting lines B0(t) and B1(t) intersect (if they do) in a single point 'p' whose distance
  //   to the lines supporting the 3 edges is exactly 't' (since those expressions are precisely parametrized in a distance)
  //   Solving the following vectorial equation:
  //
  //     [x(y),y(t)] = q + t*[l0.a,l0.b]
  //
  //   for t gives the result we want.
  //
  //
  Optional_line_2 l0 = compute_line_ceoffC2(tri->collinear_edge()) ;
  Optional_line_2 l2 = compute_line_ceoffC2(tri->non_collinear_edge()) ;

  FT sum_sq_0 = square(l0->a()) + square(l0->b());
  FT sum_sq_2 = square(l2->a()) + square(l2->b());

  Optional_point_2 q = compute_degenerate_seed_pointC2(tri);

  if ( l0 && l2 && q )
  {
    FT px, py ;
    line_project_pointC2(l0->a(),l0->b(),l0->c(),q->x(),q->y(),px,py); 
    
    CGAL_STSKEL_TRAITS_TRACE("Seed point: " << p2str(*q) << ".\nProjected seed point: (" << n2str(px) << "," << n2str(py) << ")" ) ;
    
    if ( ! CGAL_NTS is_zero(l0->b()) ) // Non-vertical
    {
      FT num = (l2->a() * l0->b() - l0->a() * l2->b() ) * px + l0->b() * l2->c() - l2->b() * l0->c();
      
      return boost::make_optional(
        Rational_time<FT>(num,
                          l0->b(), -l2->b(), l0->a() * ( l0->a()*l2->b() - l2->a()*l0->b() ),
                          sum_sq_2, sum_sq_0, FT(1)/sum_sq_0 )
      );
    }
    else
    {
      FT num = (l2->a() * l0->b() - l0->a() * l2->b() ) * py - l0->a() * l2->c() + l2->a() * l0->c() ;

      return boost::make_optional(
        Rational_time<FT>(num,
                          l2->a(), - l0->a(), l0->b()*(l0->a() * l2->b() - l0->b() * l2->a() ),
                          sum_sq_0, sum_sq_2, FT(1)/sum_sq_0)
      );
    }
  }

  return boost::none;
}

// compute the time of the intersection of the bisector of (e0,e1) with the bisector of (e2,e3) on the bisector of (e0,e1)
// in the case the bisector are not collinear and e2 and e3 are collinear
template<class K>
optional< Rational_time_4< typename K::FT> > compute_degenerate_offset_lines_isec_timeC2 ( const typename K::Segment_2& e0,
                                                                                           const typename K::Segment_2& e1,
                                                                                           const typename K::Segment_2& e2,
                                                                                           const typename K::Segment_2& e3)
{
  typedef typename K::FT FT ;

  typedef Point_2<K> Point_2 ;
  typedef Line_2 <K> Line_2 ;

  typedef optional<Line_2>  Optional_line_2 ;

  Optional_line_2 l0 = compute_line_ceoffC2(e0) ;
  Optional_line_2 l1 = compute_line_ceoffC2(e1) ;
  Optional_line_2 l2 = compute_line_ceoffC2(e2) ;

  if ( l0 && l1 && l2)
  {
    // define the bisector of e2 and e3 (free from roots)
    Point_2 bisector_pt = e2.target() == e3.source()
                        ? e2.target()
                        : validate( compute_oriented_midpoint(e2, e3) ) ;


    // (a,b,c) is a line perpedincular to the primary edge through bisector_pt.
    // If e0 and e1 are collinear this line is the actual perpendicular bisector.
    FT a = -l2->b(),
       b = l2->a(),
       c = - a * bisector_pt.x() - b * bisector_pt.y();

    FT sum_sq_0 = square(l0->a()) + square(l0->b());
    FT sum_sq_1 = square(l1->a()) + square(l1->b());

    return boost::make_optional(
      Rational_time_4<FT>(   ( a * l0->b() - l0->a() * b ) * l1->c()
                           + ( l1->a() * b - a * l1->b() ) * l0->c()
                           + ( l0->a() * l1->b() - l1->a() * l0->b() ) * c,
                         ( l1->a() * b - a * l1->b() ), a * l0->b() - l0->a() * b, 0, 0,
                         sum_sq_0, sum_sq_1, 0, 0) );
  }

  return boost::none;
}

//
// Calls the appropriate function depending on the collinearity of the edges.
//
template<class K>
optional< Rational_time< typename K::FT > > compute_offset_lines_isec_timeC2 ( intrusive_ptr< Trisegment_2<K> > const& tri )
{
  CGAL_precondition ( tri->collinearity() != TRISEGMENT_COLLINEARITY_ALL ) ;
 
  return tri->collinearity() == TRISEGMENT_COLLINEARITY_NONE ? compute_normal_offset_lines_isec_timeC2    (tri)
                                                             : compute_degenerate_offset_lines_isec_timeC2(tri);
}

//
// Calls the appropriate function depending on the collinearity of the segments.
//
template<class K>
optional< Rational_time_4< typename K::FT> > compute_offset_lines_isec_timeC2 ( const typename K::Segment_2& e0,
                                                                                const typename K::Segment_2& e1,
                                                                                const typename K::Segment_2& e2,
                                                                                const typename K::Segment_2& e3)
{
  CGAL_assertion( !collinear(e0[0], e0[1], e1[0]) || !collinear(e0[0], e0[1], e1[1]) );

  if ( !collinear(e2[0], e2[1], e3[0]) || !collinear(e2[0], e2[1], e3[1]) )
    return compute_normal_offset_lines_isec_timeC2<K>(e0, e1, e2, e3);
  else
    return compute_degenerate_offset_lines_isec_timeC2<K>(e0, e1, e2, e3);
}


// Given 3 oriented line segments e0, e1 and e2 
// such that their offsets at a certian distance intersect in a single point, 
// returns the coordinates (x,y) of such a point.
//
// PRECONDITIONS:
// None of e0, e1 and e2 are collinear (but two of them can be parallel)
// The line coefficients must be normalized: a�+b�==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
optional< Point_2<K> > construct_normal_offset_lines_isecC2 ( intrusive_ptr< Trisegment_2<K> > const& tri )
{
  typedef typename K::FT  FT ;
  
  typedef Line_2<K>  Line_2 ;
  
  typedef optional<Line_2>  Optional_line_2 ;
  
  CGAL_STSKEL_TRAITS_TRACE("Computing normal offset lines isec point for: " << tri ) ;
  
  FT x(0.0),y(0.0) ;
  
  Optional_line_2 l0 = compute_line_ceoffC2(tri->e0()) ;
  Optional_line_2 l1 = compute_line_ceoffC2(tri->e1()) ;
  Optional_line_2 l2 = compute_line_ceoffC2(tri->e2()) ;

  bool ok = false ;
  
  if ( l0 && l1 && l2 )
  {
    FT sum_sq_0 = square(l0->a()) + square(l0->b());
    FT sum_sq_1 = square(l1->a()) + square(l1->b());
    FT sum_sq_2 = square(l2->a()) + square(l2->b());

    FT r0 = CGAL_SS_i::inexact_sqrt(sum_sq_0);
    FT r1 = CGAL_SS_i::inexact_sqrt(sum_sq_1);
    FT r2 = CGAL_SS_i::inexact_sqrt(sum_sq_2);

    FT den = r0 * ( l2->a()*l1->b() - l1->a()*l2->b()) +
             r1 * ( l0->a()*l2->b() - l0->b()*l2->a()) +
             r2 * ( l0->b()*l1->a() - l0->a()*l1->b());

    CGAL_STSKEL_TRAITS_TRACE("den=" << n2str(den) )
  
    if ( ! CGAL_NTS certified_is_zero(den) ) 
    {
      FT numX = r0 * (l2->b()*l1->c() - l1->b()*l2->c() ) +
                r1 * ( l0->b()*l2->c() - l2->b()*l0->c()) +
                r2 * (  l1->b()*l0->c() - l0->b()*l1->c() );
      FT numY = r0 * (l2->a()*l1->c() - l1->a()*l2->c() ) +
                r1 * ( l0->a()*l2->c() - l2->a()*l0->c() ) +
                r2 * ( l1->a()*l0->c()- l0->a()*l1->c() );

      if ( CGAL_NTS is_finite(den) && CGAL_NTS is_finite(numX) && CGAL_NTS is_finite(numY)  )
      {
        ok = true ;
        
        x =  numX / den ;
        y = -numY / den ;
        
       CGAL_STSKEL_TRAITS_TRACE("numX=" << n2str(numX) << "\nnumY=" << n2str(numY)
                               << "\nx=" << n2str(x) << "\ny=" << n2str(y)
                               )
      }
    }
  }
    
    
  return cgal_make_optional(ok,K().construct_point_2_object()(x,y)) ;
}

// Given 3 oriented line segments e0, e1 and e2
// such that their offsets at a certian distance intersect in a single point, 
// returns the coordinates (x,y) of such a point.
// two and only two of the edges are collinear, not neccesarily consecutive but with the same orientaton
//
// PRECONDITIONS:
// The line coefficients must be normalized: a�+b�==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
optional< Point_2<K> > construct_degenerate_offset_lines_isecC2 ( intrusive_ptr< Trisegment_2<K> > const& tri )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K> Point_2 ;
  typedef Line_2<K>  Line_2 ;
  
  typedef optional<Point_2> Optional_point_2 ;
  typedef optional<Line_2>  Optional_line_2 ;
  
  CGAL_STSKEL_TRAITS_TRACE("Computing degenerate offset lines isec point for: " << tri )  ;
  
  FT x(0.0),y(0.0) ;

  Optional_line_2 l0 = compute_normalized_line_ceoffC2(tri->collinear_edge    ()) ;
  Optional_line_2 l2 = compute_normalized_line_ceoffC2(tri->non_collinear_edge()) ;
  
  Optional_point_2 q = compute_degenerate_seed_pointC2(tri);

  bool ok = false ;
  
  if ( l0 && l2 && q )
  {
    FT num, den ;
    
    FT px, py ;
    line_project_pointC2(l0->a(),l0->b(),l0->c(),q->x(),q->y(),px,py); 
    
    CGAL_STSKEL_TRAITS_TRACE("Seed point: " << p2str(*q) << ". Projected seed point: (" << n2str(px) << "," << n2str(py) << ")" ) ;
    
    if ( ! CGAL_NTS is_zero(l0->b()) ) // Non-vertical
    {
      num = (l2->a() * l0->b() - l0->a() * l2->b() ) * px + l0->b() * l2->c() - l2->b() * l0->c() ;
      den = (l0->a() * l0->a() - 1) * l2->b() + ( 1 - l2->a() * l0->a() ) * l0->b() ;
    }
    else
    {
      num = (l2->a() * l0->b() - l0->a() * l2->b() ) * py - l0->a() * l2->c() + l2->a() * l0->c() ;
      den = l0->a() * l0->b() * l2->b() - l0->b() * l0->b() * l2->a() + l2->a() - l0->a() ;
    }
  
    if ( ! CGAL_NTS certified_is_zero(den) && CGAL_NTS is_finite(den) && CGAL_NTS is_finite(num) )
    {
      x = px + l0->a() * num / den  ;
      y = py + l0->b() * num / den  ;
      
      ok = CGAL_NTS is_finite(x) && CGAL_NTS is_finite(y) ;
    }
  }
  

  CGAL_STSKEL_TRAITS_TRACE("\nDegenerate " << (CGAL_NTS is_zero(l0->b()) ? "(vertical)" : "") << " event point:  x=" << n2str(x) << " y=" << n2str(y) )

  return cgal_make_optional(ok,K().construct_point_2_object()(x,y)) ;
}

//
// Calls the appropiate function depending on the collinearity of the edges.
//
template<class K>
optional< Point_2<K> > construct_offset_lines_isecC2 ( intrusive_ptr< Trisegment_2<K> > const& tri )
{
  CGAL_precondition ( tri->collinearity() != TRISEGMENT_COLLINEARITY_ALL ) ;
  
  return tri->collinearity() == TRISEGMENT_COLLINEARITY_NONE ? construct_normal_offset_lines_isecC2    (tri)
                                                             : construct_degenerate_offset_lines_isecC2(tri) ;
}

} // namespace CGAL_SS_i

} // end namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_CONS_FTC2_H //
// EOF //

