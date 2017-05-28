// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
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
// Author(s) : Tali Zvi <talizvi@post.tau.ac.il>,
//             Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein <wein@post.tau.ac.il>
//             Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_SWEEP_LINE_SUBCURVE_H
#define CGAL_SWEEP_LINE_SUBCURVE_H

#include <CGAL/license/Sweep_line_2.h>


/*! \file
 * Defintion of the Sweep_line_subcurve class.
 */

#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Multiset.h>
#include <boost/container/flat_set.hpp>
#include <CGAL/assertions.h>

#include <boost/foreach.hpp>

namespace CGAL {

/*! \class Sweep_line_subcurve
 *
 * This is a wrapper class to X_monotone_curve_2 in the traits class, that
 * contains data that is used when applying the sweep algorithm on a set of
 * x-monotone curves.
 *
 * The information contained in this class is:
 * - the remaining x-monotone curve that is to the right of the current sweep
 *   line.
 * - two event points which are associated with the left and right end of the
 *   curve.
 * - an iterator that points to the location of the subcurve at the status line.
 * - two pointers to subcurves that are the originating subcurves in case of
 *   an overlap, otherwise thay are both NULL.
 */
template <typename Traits_>
class Sweep_line_subcurve {
public:
  typedef Traits_                                   Traits_2;
  typedef typename Traits_2::Point_2                Point_2;
  typedef typename Traits_2::X_monotone_curve_2     X_monotone_curve_2;

  typedef Sweep_line_subcurve<Traits_2>             Self;
  typedef Curve_comparer<Traits_2, Self>            Compare_curves;
  typedef Multiset<Self*, Compare_curves, CGAL_ALLOCATOR(int)>
                                                    Status_line;
  typedef typename Status_line::iterator            Status_line_iterator;

  typedef Sweep_line_event<Traits_2, Self>          Event;
  typedef boost::container::flat_set<const X_monotone_curve_2*> Input_curve_set;
protected:
  // Data members:
  X_monotone_curve_2 m_xcurve;     // The portion of the curve that lies to
                                    // the right of the last event point
                                    // that occured on the curve.

  Event* m_left_event;              // The event associated with the left end.
  Event* m_right_event;             // The event associated with the right end

  Status_line_iterator m_hint;      // The location of the subcurve in the
                                    // status line (the Y-structure).

  Input_curve_set m_input_curves;   // Contains the address of the input curves
                                    // that are containing this subcurve
                                    // SL_SAYS: in case there is no overlap, this
                                    //          should be a simple pointer

public:

  /*! Construct default.
   */
  Sweep_line_subcurve()
  {}

  /*! Construct from a curve.
   */
  Sweep_line_subcurve(const X_monotone_curve_2& curve)
  {
    m_input_curves.insert(&curve);
  }

  /*! Initialize the subcurves by setting the curve. */
  void init(const X_monotone_curve_2& curve) { m_input_curves.insert(&curve); }

  /*! Check if the given event is the matches the right-end event. */
  template <typename SweepEvent>
  bool is_end_point(const SweepEvent* event) const
  { return (m_right_event == (Event*)event); }

  /*! Get the event that corresponds to the left end of the subcurve. */
  Event* left_event() const { return m_left_event; }

  /*! Get the event that corresponds to the right end of the subcurve. */
  Event* right_event() const { return m_right_event; }

  /*! Set the event that corresponds to the left end of the subcurve. */
  template<class SweepEvent>
  void set_left_event(SweepEvent* event) { m_left_event =(Event*)event; }

  /*! Set the event that corresponds to the right end of the subcurve. */
  template<class SweepEvent>
  void set_right_event(SweepEvent* event) { m_right_event = (Event*)event; }

  /*! Get the location of the subcurve in the status line .*/
  Status_line_iterator hint() const { return m_hint; }

  /*! Set the location of the subcurve in the status line .*/
  void set_hint(Status_line_iterator hint) { m_hint = hint; }

  /*! returns the X-monotone curve associated with the subcurve */
  const X_monotone_curve_2& x_monotone_curve() const
  {
    return m_xcurve;
  }

  /*! set the X-monotone curve associated with the subcurve */
  void set_x_monotone_curve(const X_monotone_curve_2& xmc)
  {
    m_xcurve=xmc;
  }

  /*! set as the X-monotone curve a parent curve*/
  void set_x_monotone_curve()
  {
    CGAL_assertion(!m_input_curves.empty());
    m_xcurve=*(*m_input_curves.begin());
  }

  /*! incidcate whether the subcurve represents an overlap */
  bool is_overlap() const
  {
    return m_input_curves.size() > 1;
  }

  bool equal(const Self& other) const
  {
    if (left_event()!=other->left_event() || right_event()!=other.right_event())
      return false;
    return m_input_curves==other.m_input_curves;
  }

  bool has_common_input_curve(const Self& other) const
  {
    BOOST_FOREACH(const X_monotone_curve_2* xmc1, m_input_curves)
      BOOST_FOREACH(const X_monotone_curve_2* xmc2, other.m_input_curves)
      {
        if (xmc1==xmc2)
          return true;
        if (xmc2>xmc1) break;
      }
    return false;
  }

  void set_input_curves(const Self& other)
  {
    m_input_curves=other.m_input_curves;
  }

  void merge_input_curves(const Self& other)
  {
    m_input_curves.insert(other.m_input_curves.begin(),
                          other.m_input_curves.end());
  }

  void swap(Self& other)
  {
    std::swap(m_xcurve, other.m_xcurve);
    std::swap(m_left_event, other.m_left_event);
    std::swap(m_right_event, other.m_right_event);
    std::swap(m_hint, other.m_hint);
    std::swap(m_input_curves, other.m_input_curves);
  }

  const X_monotone_curve_2& x_monotone_input_curve() const
  {
    return *(*m_input_curves.begin());
  }

  // /*! Get all the leaf-nodes in the hierarchy of overlapping subcurves. */
  // template <typename OutputIterator>
  // OutputIterator all_leaves(OutputIterator oi)
  // {
  //   if (m_orig_subcurve1 == NULL) {
  //     *oi++ = this;
  //     return oi;
  //   }
  // 
  //   oi = m_orig_subcurve1->all_leaves(oi);
  //   oi = m_orig_subcurve2->all_leaves(oi);
  //   return oi;
  // }

  // /*! Check if the given subcurve is a node in the overlapping hierarchy. */
  // bool is_inner_node(Self *s)
  // {
  //   if (this == s) return true;
  //   if (m_orig_subcurve1 == NULL) return false;
  //   return (m_orig_subcurve1->is_inner_node(s) ||
  //           m_orig_subcurve2->is_inner_node(s));
  // }
  // 
  // /*! Check if the given subcurve is a leaf in the overlapping hierarchy. */
  // bool is_leaf(Self* s)
  // {
  //   if (m_orig_subcurve1 == NULL) return (this == s);
  //   return (m_orig_subcurve1->is_leaf(s) || m_orig_subcurve2->is_leaf(s));
  // }
  // 
  // /*! Check if the two hierarchies contain the same leaf nodes. */
  // bool has_same_leaves(Self* s)
  // {
  //   std::list<Self*> my_leaves;
  //   std::list<Self*> other_leaves;
  // 
  //   this->all_leaves(std::back_inserter(my_leaves));
  //   s->all_leaves(std::back_inserter(other_leaves));
  // 
  //   typename std::list<Self*>::iterator iter;
  //   for (iter = my_leaves.begin(); iter != my_leaves.end(); ++iter) {
  //     if (std::find(other_leaves.begin(), other_leaves.end(), *iter) ==
  //         other_leaves.end())
  //       return false;
  //   }
  // 
  //   for (iter = other_leaves.begin(); iter != other_leaves.end(); ++iter) {
  //     if (std::find(my_leaves.begin(), my_leaves.end(), *iter) ==
  //         my_leaves.end())
  //       return false;
  //   }
  // 
  //   return true;
  // }
  // 
  // /*! Check if the two hierarchies contain a common leaf node. */
  // bool has_common_leaf(Self *s)
  // {
  //   std::list<Self*> my_leaves;
  //   std::list<Self*> other_leaves;
  // 
  //   this->all_leaves(std::back_inserter(my_leaves));
  //   s->all_leaves(std::back_inserter(other_leaves));
  // 
  //   typename std::list<Self*>::iterator iter;
  //   for (iter = my_leaves.begin(); iter != my_leaves.end(); ++iter) {
  //     if (std::find(other_leaves.begin(), other_leaves.end(), *iter) !=
  //         other_leaves.end())
  //       return true;
  //   }
  //   return false;
  // }
  // 
  // /*! Get all distinct nodes from the two hierarchies. */
  // template <typename OutputIterator>
  // OutputIterator distinct_nodes(Self* s, OutputIterator oi)
  // {
  //   if (m_orig_subcurve1 == NULL) {
  //     if (s->is_leaf(this)) *oi++ = this;
  //     return oi;
  //   }
  // 
  //   if (! s->is_inner_node (m_orig_subcurve1)) *oi++ = m_orig_subcurve1;
  //   else oi++ = m_orig_subcurve1->distinct_nodes(s, oi);
  // 
  //   if (! s->is_inner_node (m_orig_subcurve2)) *oi++ = m_orig_subcurve2;
  //   else oi++ = m_orig_subcurve2->distinct_nodes(s, oi);
  // 
  //   return oi;
  // }
  // 
  /*! Get the depth of the overlap hierarchy. */
  unsigned int overlap_depth()
  {
    return m_input_curves.size();
  }

#ifdef CGAL_SL_VERBOSE
  void Print() const;
#endif
};

#ifdef CGAL_SL_VERBOSE
  template<class Traits>
  void Sweep_line_subcurve<Traits>::Print() const
  {
    /// SL_SAYS: TODO
    // std::cout << "Curve " << this
    //           << "  (" << m_lastCurve << ") "
    //           << " [sc1: " << m_orig_subcurve1
    //           << ", sc2: " << m_orig_subcurve2 << "]";
  }
#endif

} //namespace CGAL

#endif
