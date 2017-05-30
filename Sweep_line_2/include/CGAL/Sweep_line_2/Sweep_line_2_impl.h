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
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Efi Fogel       <efif@post.tau.ac.il>
//               (based on old version by Tali Zvi)

#ifndef CGAL_SWEEP_LINE_2_IMPL_H
#define CGAL_SWEEP_LINE_2_IMPL_H

#include <CGAL/license/Sweep_line_2.h>

#include <CGAL/utility.h>

/*! \file
 * Member-function definitions of the Sweep_line_2 class-template.
 */

namespace CGAL {

//-----------------------------------------------------------------------------
// Initialize the data structures for the sweep-line algorithm.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::_init_structures()
{
  // Initailize the structures maintained by the base sweep-line class.
  Base::_init_structures();

  // Resize the hash to be O(2*n), where n is the number of input curves.
  // SL_SAYS: TODO check with Efi
  m_curves_pair_set = Curve_pair_set(2 * this->m_num_of_subCurves);
}

//-----------------------------------------------------------------------------
// Complete the sweep (complete the data structures).
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::_complete_sweep()
{
  CGAL_SL_PRINT_START_EOL("completing the sweep");

  // Complete the sweep process using base sweep-line class.
  Base::_complete_sweep();

  // Clean the set of curve pairs for which we have computed intersections.
  m_curves_pair_set.clear();

  CGAL_SL_PRINT_END_EOL("completing the sweep");
}

//-----------------------------------------------------------------------------
// Handle the subcurves to the left of the current event point.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::_handle_left_curves()
{
  for(Status_line_iterator sliter=this->m_statusLine.begin();
                           sliter!=this->m_statusLine.end(); ++sliter)
  {
    CGAL_assertion((*sliter)->hint()==sliter);
  }

  
  
  CGAL_SL_PRINT_START("handling left curves at (");
  CGAL_SL_DEBUG(this->PrintEvent(this->m_currentEvent));
  CGAL_SL_PRINT_TEXT(")");
  CGAL_SL_PRINT_EOL();

  this->m_is_event_on_above = false;

  if (! this->m_currentEvent->has_left_curves()) {
    // In case the current event has no left subcurves incident to it, we have
    // to locate a place for it in the status line.
    CGAL_SL_PRINT_TEXT("Handling case: no left curves");
    CGAL_SL_PRINT_EOL();
    this->_handle_event_without_left_curves();

    Status_line_iterator sl_pos = this->m_status_line_insert_hint;

    if (this->m_is_event_on_above) {
      CGAL_SL_PRINT_TEXT("The event is on a curve in the status line");
      CGAL_SL_PRINT_EOL();

      // The current event point starts at the interior of a subcurve that
      // already exists in the status line (this may also indicate an overlap).
      if (! this->m_currentEvent->has_right_curves()) {
        // The event is an isolated point.
        if (this->m_currentEvent->is_query()) {
          // In case of a query point, just notify the visitor about it.
          this->m_is_event_on_above = true;
          this->m_visitor->before_handle_event(this->m_currentEvent);
          return;
        }

        CGAL_assertion(this->m_currentEvent->is_action());
      }

      // Split the curve that has the point on it. We create the two corresponding
      // subcurves and update the events and update the curve in the status line  
      this->m_currentEvent->set_weak_intersection();

      Status_line_iterator sl_it = this->m_status_line_insert_hint; // position of the curve the point is on.
      Subcurve* sc_ptr = static_cast<Subcurve*>(*sl_it); // the corresponding subcurve that will be replaced
      
      const X_monotone_curve_2& xm_curve_to_split = sc_ptr->x_monotone_curve();
      this->m_traits->split_2_object()(xm_curve_to_split,
                                       this->m_currentEvent->point(),
                                       sub_cv1, sub_cv2);

      // create the first subcurves
      Subcurve left_subcurve;
      left_subcurve.set_left_event(sc_ptr->left_event());
      left_subcurve.set_right_event(this->m_currentEvent);
      left_subcurve.set_x_monotone_curve(sub_cv1);
      left_subcurve.set_input_curves(*sc_ptr);
      // create the second subcurves
      Subcurve right_subcurve;
      right_subcurve.set_left_event(this->m_currentEvent);
      right_subcurve.set_right_event(sc_ptr->right_event());
      right_subcurve.set_x_monotone_curve(sub_cv2);
      right_subcurve.set_input_curves(*sc_ptr);
      // remove the subcurve to be replaced from the former endpoint
      right_subcurve.right_event()->remove_curve_from_left(*(*sl_it));
      // add it to the left curve of the event and update the status line
      *sl_it = this->m_currentEvent->push_back_curve_to_left(left_subcurve);
      //SL_SAYS : TODO Ask Efi what the visitor expects
      this->m_visitor->update_event(this->m_currentEvent, reinterpret_cast<Subcurve*>(*sl_it));

      // If necessary, add the subcurves as a right incident curve as well.
      // We also check for overlaps.
      bool is_overlap = _add_curve_to_right(this->m_currentEvent, right_subcurve);

      ++(this->m_status_line_insert_hint);

      CGAL_SL_PRINT_ERASE(*sl_pos);
      this->m_statusLine.erase(sl_pos);

      if (is_overlap) {
        // Handle overlaps.
        /// SL_SAYS TODO I think that sub_cv1 is no longer setted correctly in that case
        this->m_visitor->add_subcurve(sub_cv1, sc_ptr);
        CGAL_SL_PRINT_END_EOL("handling left curves");
        return;
      }
    }
    else {
      // The event is not located on any subcurve.
      this->m_visitor->before_handle_event(this->m_currentEvent);
      CGAL_SL_PRINT_END_EOL("handling left curves");
      return;
    }
  }

  CGAL_SL_PRINT_TEXT("left curves before sorting:");
  CGAL_SL_PRINT_EOL();
  CGAL_SL_DEBUG(if (this->m_currentEvent->left_curves_begin() !=
                    this->m_currentEvent->left_curves_end())
                { this->print_event_info(this->m_currentEvent); });

  this->_sort_left_curves();
  this->m_visitor->before_handle_event(this->m_currentEvent);

  CGAL_SL_PRINT_TEXT("left curves after sorting:");
  CGAL_SL_PRINT_EOL();
  CGAL_SL_DEBUG(if (this->m_currentEvent->left_curves_begin() !=
                    this->m_currentEvent->left_curves_end() )
                { this->print_event_info(this->m_currentEvent); });

  // Remove all left subcurves from the status line, and inform the visitor
  // that we are done handling these subcurves.
  Event_subcurve_iterator left_iter = this->m_currentEvent->left_curves_begin();
  while (left_iter != this->m_currentEvent->left_curves_end()) {
    Subcurve& left_sc = *left_iter;
    this->m_visitor->add_subcurve(left_sc.x_monotone_curve(), &left_sc);
    ++left_iter;

    this->_remove_curve_from_status_line(left_sc);
  }
  
  // SL_SAYS TODO shall we reemove subcurves from the left of the event???

  CGAL_SL_PRINT_END_EOL("handling left curves");
}

//-----------------------------------------------------------------------------
// Handle the subcurves to the right of the current event point.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::_handle_right_curves()
{
  CGAL_SL_PRINT_START("handling right curves at (");
  CGAL_SL_DEBUG(this->PrintEvent(this->m_currentEvent));
  CGAL_SL_PRINT_TEXT(")");
  CGAL_SL_PRINT_EOL();
  
  if (! this->m_currentEvent->has_right_curves()) {
    CGAL_SL_PRINT_END_EOL("handling right curves");
    return;
  }
  
  // Loop over the curves to the right of the status line and handle them:
  // - If we are at the beginning of the curve, we insert it to the status
  //   line, then we look if it intersects any of its neighbors.
  // - If we are at an intersection point between two curves, we add them
  //   to the status line and attempt to intersect them with their neighbors
  // - We also check to see if the two intersect again to the right of the
  //   point. <----- SL_SAYS since all intersections have already been computed
  //                 the comment seem wrong here  
  Event_subcurve_iterator currentOne =
    this->m_currentEvent->right_curves_begin();
  Event_subcurve_iterator rightCurveEnd =
    this->m_currentEvent->right_curves_end();
  

  Subcurve& sc = *currentOne;
  CGAL_assertion( reinterpret_cast<Event*>(sc.left_event())==this->m_currentEvent );
  CGAL_SL_PRINT_INSERT(&sc);
  Status_line_iterator slIter =
    this->m_statusLine.insert_before(this->m_status_line_insert_hint, &sc);
  sc.set_hint(slIter);
  *slIter=sc.right_event()->push_back_curve_to_left(sc); // status line contains only subcurve that are a left curve of an event
  
  CGAL_SL_PRINT_STATUS_LINE();
  if (slIter != this->m_statusLine.begin()) {
    //  get the previous curve in the y-str
    Status_line_iterator prev = cpp11::prev(slIter);
    _intersect(static_cast<Subcurve*>(*prev), static_cast<Subcurve*>(*slIter));
  }
  
  Subcurve* prev_sc = reinterpret_cast<Subcurve*>(*slIter);
  ++currentOne;
  while (currentOne != rightCurveEnd) {
    Subcurve& sc = *currentOne;
  
    CGAL_SL_PRINT_INSERT(&sc);
    slIter =
      this->m_statusLine.insert_before(this->m_status_line_insert_hint, &sc);
    sc.set_hint(slIter);
    *slIter=sc.right_event()->push_back_curve_to_left(sc); // status line contains only subcurve that are a left curve of an event
  
  // SL_SAYS: we should use the curve from the left container!
  
    CGAL_SL_PRINT_STATUS_LINE();
  
    // If the two curves used to be neighbours before, we do not need to
    // intersect them again.
    // SL_SAYS TODO find a way to reenable the test
    // if (!this->m_currentEvent->are_left_neighbours
    //     (static_cast<Subcurve*>(*currentOne), static_cast<Subcurve*>(*prevOne)))
    {
      _intersect(prev_sc, reinterpret_cast<Subcurve*>(*slIter));
    }
  
    prev_sc = reinterpret_cast<Subcurve*>(*slIter);
    ++currentOne;
  }
  
  CGAL_SL_PRINT_STATUS_LINE();
  
  //the next Subcurve at the status line
  ++slIter;
  if (slIter != this->m_statusLine.end())
    _intersect(static_cast<Subcurve*>(prev_sc),
               static_cast<Subcurve*>(*slIter));

  this->m_currentEvent->clear_curves();
  
  CGAL_SL_PRINT_END_EOL("handling right curves");
}

//-----------------------------------------------------------------------------
// Add a subcurve to the right of an event point.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
bool Sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_add_curve_to_right(Event* event, Subcurve& curve, bool overlap_exist)
{
  CGAL_SL_PRINT_START("adding a Curve to the right of (");
  CGAL_SL_DEBUG(this->PrintEvent(event));
  CGAL_SL_PRINT_TEXT(") ");
  CGAL_SL_PRINT_CURVE(&curve);
  CGAL_SL_PRINT_EOL();

  CGAL_assertion(reinterpret_cast<Event*>(curve.left_event())==event);

  Event_subcurve_iterator iter;
  for (iter = event->right_curves_begin(); iter != event->right_curves_end();
       ++iter)
  {
    CGAL_SL_PRINT_CURVE(&(*iter));
    CGAL_SL_PRINT_EOL();
    
    if (curve.has_common_input_curve(*iter))
    {
      if (curve.right_event()==iter->right_event())
      {
        CGAL_SL_PRINT_END_EOL("adding a Curve to the right (complete overlap)");
        iter->merge_input_curves(curve);
        return false;
      }
      // partial overlap, need to split the curve
      if (this->m_queueEventLess(
            reinterpret_cast<Event*>(curve.right_event()), 
            reinterpret_cast<Event*>(iter->right_event())))
        curve.swap(*iter);

      iter->merge_input_curves(curve);
      curve.set_left_event(iter->right_event());
      _add_curve_to_right(reinterpret_cast<Event*>(iter->right_event()), curve, false); // SL_SAYS : TODO NOT SURE ABOUT the value of overlap_exist here
      return false;
    }
  }

  std::pair<bool, Event_subcurve_iterator> pair_res =
    event->add_curve_to_right(curve, this->m_traits);
  bool rc_overlap = pair_res.first;
  Event_subcurve_iterator rc_it = pair_res.second;
  
  if (! rc_overlap) {
    // No overlap occurs.
    CGAL_SL_PRINT_END_EOL("adding a Curve to the right (no overlap)");
    return false;
  }
  
  _handle_overlap(event, curve, rc_it, overlap_exist);
  
  // Inidicate that an overlap has occured:
  CGAL_SL_PRINT_END_EOL("adding a Curve to the right (overlap)");
  return true;
}

//-----------------------------------------------------------------------------
// Compute intersections between the two given curves.
//
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::_intersect(Subcurve* c1,
                                                           Subcurve* c2)
{
  // SL_SAYS TODO optimization in the case we know the curve overlaps and shares their left endpoints
  
  CGAL_SL_PRINT_START("computing intersection of ");
  CGAL_SL_PRINT_CURVE(c1);
  CGAL_SL_PRINT_TEXT(" and ");
  CGAL_SL_PRINT_CURVE(c2);
  CGAL_SL_PRINT_EOL();

  CGAL_assertion(*(c1->hint())==c1);
  CGAL_assertion(*(c2->hint())==c2);

  typedef typename Tr::Multiplicity Multiplicity;

  CGAL_assertion(c1 != c2);

  // look up for (c1,c2) in the table and insert if doesnt exist
  if (!c1->is_overlap() && !c2->is_overlap())
    if (! m_curves_pair_set.insert(
            make_sorted_pair( &(c1->x_monotone_input_curve()),
                              &(c2->x_monotone_input_curve()) ) ).second )
    {
      CGAL_SL_PRINT_END_EOL("computing intersection (intersection already inserted) ");
      return;  //the curves have already been checked for intersection
    }

  std::vector<Object> x_objects;
  this->m_traits->intersect_2_object()(c1->x_monotone_curve(), 
                                       c2->x_monotone_curve(),
                                       std::back_inserter(x_objects)); 

  if (x_objects.empty()) {
    CGAL_SL_PRINT_END_EOL("Computing intersection (no intersection)");
    return; // no intersection at all
  }

  typename std::vector<Object>::iterator vi = x_objects.begin(),
                                         vi_end=x_objects.end();

  // If the two subcurves have a common left-event, and the first intersection
  // object is a point, we can ignore the first intersection (note that in case of
  // an overlap that ends at the common endpoint, we definately want to keep
  // the intersection object).
  if (c1->left_event() == c2->left_event())                                          
  {
    /// SL_SAYS : TODO we know that this can only happen when called from _add_curve_to_right->_handle_overlap
    if (object_cast<std::pair<Point_2, Multiplicity> >(&(*vi)) != NULL)
    {
      CGAL_SL_PRINT_TEXT("Skipping common left endpoint ...");
      CGAL_SL_PRINT_EOL();
      ++vi;
    }
  }

  // If the two subcurves have a common right-event, and the last intersection
  // object is a point, we can ignore last intersection (note that in case of
  // an overlap that ends at the common endpoint, we definately want to keep
  // the intersection object).
  if (c1->right_event() == c2->right_event())
  {
    if (object_cast<std::pair<Point_2, Multiplicity> >(&(*cpp11::prev(vi_end))) != NULL) {
      CGAL_SL_PRINT_TEXT("Skipping common right endpoint...");
      CGAL_SL_PRINT_EOL();
      --vi_end;
    }
  }

  const std::pair<Point_2,Multiplicity>* xp_point;

  //Note : no need to filter intersection points, x-monotone curves are tight
  bool first_intersection=true;
  Subcurve* original_c1 = NULL, * original_c2=NULL;
  Event* original_re_1 = NULL, * original_re_2=NULL;
  for (; vi != vi_end; ++vi) {
    
    if (first_intersection)
    {
      original_c1=c1;
      original_c2=c2;
      original_re_1=c1->right_event();
      original_re_2=c2->right_event();
    }
    
    const X_monotone_curve_2* icv;
    Point_2 xp;
    unsigned int multiplicity = 0;

    xp_point = object_cast<std::pair<Point_2, Multiplicity> >(&(*vi));
    Event* first_event = NULL;
    if (xp_point != NULL) {
      xp = xp_point->first;
      multiplicity = xp_point->second;
      CGAL_SL_PRINT_TEXT("Found an intersection point");
      CGAL_SL_PRINT_EOL();
      first_event = _create_intersection_point(xp, multiplicity, c1, c2, false);
    }
    else {
      icv = object_cast<X_monotone_curve_2>(&(*vi));
      CGAL_assertion(icv != NULL);
      CGAL_SL_PRINT_TEXT("Found an overlap");
      CGAL_SL_PRINT_EOL();

      // TODO EBEB: This code does not work with overlaps that reach the boundary
      Point_2 left_xp = this->m_traits->construct_min_vertex_2_object()(*icv);
      xp = this->m_traits->construct_max_vertex_2_object()(*icv);

      sub_cv1 = *icv;
      _create_intersection_point(xp, 0 , c1 , c2, false);
      first_event = _create_intersection_point(left_xp, 0 , c1 ,c2, true);
    }
    if (first_intersection)
    {
      std::cout << "first intersection! --- " <<  first_event->point() << "\n";
      
      
      //update the right event of the curves
      c1->set_right_event(first_event);
      CGAL_assertion(c1->right_event()==first_event);
      c2->set_right_event(first_event);
      c1->set_hint(original_c1->hint());
      c2->set_hint(original_c2->hint());
      //replace them in the status line
      *(original_c1->hint()) = first_event->push_back_curve_to_left(*c1);
      *(original_c2->hint()) = first_event->push_back_curve_to_left(*c2);
      original_re_1->remove_curve_from_left(*original_c1);
      original_re_2->remove_curve_from_left(*original_c2);
      first_intersection=false;
    }
  }

  CGAL_SL_PRINT_END_EOL("computing intersection");
}

//-----------------------------------------------------------------------------
// Create an intersection-point event between two curves.
//
// SL_SAYS : TODO why references to pointers????
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
typename Sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::Event*
Sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_create_intersection_point(const Point_2& xp,
                           unsigned int multiplicity,
                           Subcurve* c1, Subcurve* c2,
                           bool is_overlap)
{
  CGAL_SL_PRINT_START_EOL("creating an intersection point between");
  CGAL_SL_PRINT_CURVE(c1);
  CGAL_SL_PRINT_EOL();
  CGAL_SL_PRINT_CURVE(c2);
  CGAL_SL_PRINT_EOL();

  // insert the event and check if an event at this point already exists.
  const std::pair<Event*, bool>& pair_res =
    this->_push_event(xp, Base_event::DEFAULT, ARR_INTERIOR, ARR_INTERIOR);

  Event* e = pair_res.first;
  if (pair_res.second) {
    // a new event is created , which indicates that the intersection point
    // cannot be one of the end-points of two curves
    CGAL_SL_PRINT_TEXT("A new event is created .. (");
    CGAL_SL_PRINT(xp);
    CGAL_SL_PRINT_TEXT(")");
    CGAL_SL_PRINT_EOL();

    e->set_intersection();

    this->m_visitor->update_event(e, c1, c2, true);

    // Act according to the multiplicity:
    if (multiplicity == 0) {
      // The multiplicity of the intersection point is unkown or undefined:
      _add_curve_to_right(e, *c1, is_overlap);
      _add_curve_to_right(e, *c2, is_overlap);
      if (! is_overlap) {
        if (e->is_right_curve_bigger(c1, c2)) std::swap(c1, c2);
      }
    }
    else {
      if ((multiplicity % 2) == 1) {
        // The mutiplicity of the intersection point is odd: Swap their
        // order to the right of this point.
        e->add_curve_pair_to_right(*c2, *c1);
      }
      else {
        // The mutiplicity of the intersection point is even, so they
        // maintain their order to the right of this point.
        CGAL_assertion((multiplicity % 2) == 0);
        e->add_curve_pair_to_right(*c1, *c2);
      }
    }
  }
  else {
    // The event already exists, so we need to update it accordingly
    CGAL_SL_PRINT_TEXT("Event already exists, updating.. (");
    CGAL_SL_PRINT(xp);
    CGAL_SL_PRINT_TEXT(")");
    CGAL_SL_PRINT_EOL();
    if (e == this->m_currentEvent) {
      // This can happen when c1 starts at the interior of c2 (or vice versa).
      return e;
    }

    if (!c1->is_end_point(e) && !c2->is_end_point(e)) {
      _add_curve_to_right(e, *c1, is_overlap);
      _add_curve_to_right(e, *c2, is_overlap);
      e->set_intersection();
      this->m_visitor->update_event(e, c1, c2, false);
    }
    else {
      if (!c1->is_end_point(e) && c2->is_end_point(e)) {
        _add_curve_to_right(e, *c1, is_overlap);
        e->set_weak_intersection();
        this->m_visitor->update_event(e, c1);
      }
      else {
        if (c1->is_end_point(e) && !c2->is_end_point(e)) {
          _add_curve_to_right(e, *c2, is_overlap);
          e->set_weak_intersection();
          this->m_visitor->update_event(e, c2);
        }
      }
    }
    if (! is_overlap) {
      if (e->is_right_curve_bigger(c1, c2)) std::swap(c1, c2);
    }

    CGAL_SL_PRINT_EVENT_INFO(e);
  }

  return e;
  CGAL_SL_PRINT_END_EOL("Creating an intersection point");
}

/*! Handle overlap at right insertion to event.
 * \param event the event where that overlap starts (the left end point of the
 *        overlap).
 * \param curve the subcurve that its insertion to the list of right subcurves of
 *        'event' causes the overlap (with *iter).
 * \param iter the existing subcurve at the right subcurves of 'event'
 * \param[in] overlap_exist indicates whether information about the
 *            overlapping curve has been computed already. If overlap_exist
 *            is true, sub_cv1 is a computed X_monotone_curve_2 curve that
 *            should be used as the overlapping curve.
 */
template <typename Tr, typename Vis, typename Subcv, typename Evnt,
          typename Alloc>
void Sweep_line_2<Tr, Vis, Subcv, Evnt, Alloc>::
_handle_overlap(Event* event, Subcurve& curve, Event_subcurve_iterator iter,
                bool overlap_exist)
{
  // An overlap occurs:
  CGAL_SL_PRINT_START("handling overlap at right insertion of (");
  CGAL_SL_DEBUG(this->PrintEvent(event));
  CGAL_SL_PRINT_TEXT(") ");
  CGAL_SL_PRINT_EOL();
  CGAL_SL_PRINT_CURVE(&curve);
  CGAL_SL_PRINT_EOL();
  CGAL_SL_PRINT_CURVE(&(*iter));
  CGAL_SL_PRINT_EOL();

  X_monotone_curve_2 overlap_cv;
  if (overlap_exist) {
    CGAL_SL_PRINT_TEXT("  overlap already computed.");
    CGAL_SL_PRINT_EOL();
    overlap_cv = sub_cv1;
    // std::cout << "Overlap computed!!!!! " << overlap_cv << std::endl;
  }
  else {
    // compute the overlap.
    // --> this part of the code is only called when 
    
    CGAL_SL_PRINT_TEXT("  computing overlap.");
    CGAL_SL_PRINT_EOL();

    // compute intersection points and overlaps
    _intersect(&curve, &(*iter));
    return;
// SL_SAYS : I replace the following code with a call to _intersect(). To be cross-checked with Efi!
//           BTW if correct, I think the function should be split.
#if 0
    CGAL_assertion(iter != event->right_curves_end());
    std::vector<Object>  obj_vec;
    vector_inserter vit(obj_vec);
    this->m_traits->intersect_2_object()(curve->last_curve(),
                                         (*iter)->last_curve(),
                                         vit);

    if (obj_vec.empty()) {
      CGAL_SL_PRINT_END_EOL("handling overlap");
      return;
    }

    std::size_t obj_vec_size = obj_vec.size();
    if (1 == obj_vec_size) {
      //always a curve since an overlap was detected on the right of `event`
      CGAL_assertion( obj_vec.front().is<X_monotone_curve_2>() );
      overlap_cv = object_cast<X_monotone_curve_2>(obj_vec.front());
    }
    else {
      CGAL_SL_PRINT_TEXT("Overlap consists of more than one curve");
      CGAL_SL_PRINT_EOL();
      // Elements in obj_vec are sorted using less-xy.
      for (std::size_t i = 0; i < obj_vec_size; ++i) {
        const CGAL::Object& obj = obj_vec[i];
        if (const X_monotone_curve_2* xcv_ptr =
            object_cast<X_monotone_curve_2>(&obj))
        {
          // We are only interested in overlapping curves.
          // Look for a curve containing the event point.
          if (this->m_traits->compare_xy_2_object()
              (this->m_traits->construct_max_vertex_2_object()(*xcv_ptr),
               event->point()) == LARGER)
          {
            CGAL_assertion(this->m_traits->compare_xy_2_object()
                           (this->m_traits->construct_min_vertex_2_object()
                            (*xcv_ptr), event->point()) != LARGER);
            overlap_cv = *xcv_ptr;
            break;
          }
          else {
            CGAL_SL_PRINT_TEXT("  Skip a curve");
            CGAL_SL_PRINT_EOL();
          }
        }
        else {
          // We are not interested in intersection points.
          CGAL_SL_PRINT_TEXT("  Skip a point");
          CGAL_SL_PRINT_EOL();
        }
      }
    }
#endif
  }

  // Get the right end of overlap_cv (if it is closed from the right).
  Event* right_end;
  Arr_parameter_space  ps_x_r =
    this->m_traits->parameter_space_in_x_2_object()(overlap_cv, ARR_MAX_END);
  Arr_parameter_space  ps_y_r =
    this->m_traits->parameter_space_in_y_2_object()(overlap_cv, ARR_MAX_END);

  CGAL_assertion(ps_x_r != ARR_LEFT_BOUNDARY);
  if ((ps_x_r != ARR_INTERIOR) || (ps_y_r != ARR_INTERIOR)) {
    // The overlapping subcurve is either open from the right, or
    // touches the boundary of the surface. In either case, the curves that
    // are involved in the overlap must also be open or defined at the
    // boundary, so the event associated with their right ends already exists,
    // and we set it as the overlapping subcurve's right event.
    CGAL_assertion(iter->right_event() == curve.right_event());
    right_end = (Event*)(curve.right_event());
  }
  else {
    // The overlapping subcurve has a valid right endpoint.
    // Find the event associated with this point (or create a new event).
    const Point_2& end_overlap =
      this->m_traits->construct_max_vertex_2_object()(overlap_cv);

    const std::pair<Event*, bool>& pair_res =
      this->_push_event(end_overlap, Base_event::OVERLAP, ps_x_r, ps_y_r);

    right_end = pair_res.first;
  }

#if 0
/// SL_SAYS I think this part is no longer needed since x-monotone curves
///         associated to subcurves are tight now
  // Get the left end of overlap_cv (if it is closed from the left).
  Arr_parameter_space  ps_x_l =
    this->m_traits->parameter_space_in_x_2_object()(overlap_cv, ARR_MIN_END);
  Arr_parameter_space  ps_y_l =
    this->m_traits->parameter_space_in_y_2_object()(overlap_cv, ARR_MIN_END);

  CGAL_assertion(ps_x_l != ARR_RIGHT_BOUNDARY);
  if ((ps_x_l == ARR_INTERIOR) && (ps_y_l == ARR_INTERIOR)) {
    // The left end of the overlapping subcurve is regular point, so in case
    // the event is also associated with a regular point (not incident to the
    // surface boundaries), we make sure that the overlapping subcurve does
    // not start to the left of this event.
    if (! event->is_on_boundary()) {
      // If the left endpoint of the overlapping curve is to the left of the
      // event, split the overlapping subcurve so its left endpoint equals
      // the event point.
      const Point_2& begin_overlap =
        this->m_traits->construct_min_vertex_2_object()(overlap_cv);
      Comparison_result res =
        this->m_traits->compare_xy_2_object()(event->point(), begin_overlap);

      CGAL_assertion(res != SMALLER);
      if (res == LARGER) {
        this->m_traits->split_2_object()(overlap_cv, event->point(),
                                         sub_cv1, sub_cv2);
        overlap_cv = sub_cv2;
      }
    }
  }
  else {
    // The left end of the overlapping subcurve is either open, or
    // incident to the surface boundaries. In case the current event is
    // associated with a regular point, it must lie to the right of this
    // curve-end, so we clip the overlapping subcurve accordingly.
    if (! event->is_on_boundary()) {
      this->m_traits->split_2_object()(overlap_cv, event->point(),
                                       sub_cv1, sub_cv2);
      overlap_cv = sub_cv2;
    }
  }
#endif

  // Create the subcurve
  Subcurve overlap_sc;
  overlap_sc.set_hint(this->m_statusLine.end());
  overlap_sc.set_left_event(event);
  overlap_sc.set_right_event(right_end);
  overlap_sc.set_x_monotone_curve(overlap_cv);
  overlap_sc.set_input_curves(*iter);
  overlap_sc.merge_input_curves(curve);

  // Set the two events' attribute to overlap.
  event->set_overlap();

  CGAL_SL_PRINT_CURVE(&curve);
  CGAL_SL_PRINT_TEXT(" + ");
  CGAL_SL_PRINT_CURVE(&(*iter));
  CGAL_SL_PRINT_TEXT(" => ");
  CGAL_SL_PRINT_EOL();
  CGAL_SL_PRINT_TEXT("  ");
  CGAL_SL_PRINT_CURVE(&overlap_sc);
  CGAL_SL_PRINT_EOL();

  _add_curve_to_right(event, overlap_sc);

  CGAL_SL_PRINT_END_EOL("handling overlap");
}

} //namespace CGAL

#endif
