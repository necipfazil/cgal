// Copyright (c) 2018 GeometryFactory (France).
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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Konstantinos Katrioplas


#ifndef CGAL_ISLAND_TRIANGULATE_HOLE_POLYLINE_H
#define CGAL_ISLAND_TRIANGULATE_HOLE_POLYLINE_H

#include <vector>
#include <tuple>
#include <stack>
#include <CGAL/Combination_enumerator.h>



namespace CGAL {
namespace internal {


// Domain structure //
// ---------------- //

template <typename PointRange>
struct Domain
{
  Domain() {}

  // boundary should always be given, the first and last point being different
  // (boundary is considered to be closed)
  Domain(PointRange& boundary)
    : boundary(boundary.begin(),
               boundary.end())
  {
    CGAL_assertion(!boundary.empty() && boundary.back() != boundary.front());
  }

  // constructor with indices
  Domain(const std::vector<int>& ids) : b_ids(ids) {}

  void add_hole(const std::vector<int>& ids)
  {
    holes_list.push_back(ids);

    for(int i=0; i<ids.size(); ++i)
      h_ids.push_back(ids[i]);
  }

  //SL_comments: shouldn't you test also if the boundary is empty?
  bool is_empty() const
  {
    holes_list.empty() ? true : false;
  }

  std::pair<int, int> get_access_edge() const
  {
    CGAL_assertion( b_ids.size()>1 );
    return std::make_pair(b_ids.back(), b_ids.front());
  }

  void print_boundary() const
  {
    print(boundary);
  }

  void print_size() const
  {
    std::cout << boundary.size() << std::endl;
  }

  //data
  // SL_comments: PointRange shouldn't be a template parameter and you should use a vector directly
  PointRange boundary; // not used in the main algorithm
  std::vector<PointRange> holes; // not used

  std::vector<int> b_ids; // comment?
  std::vector<int> h_ids; // comment?
  std::vector<std::vector<int> > holes_list;

};


template <typename PointRange>
void print(PointRange &v)
{
  for(int i=0; i<v.size(); ++i)
  {
    std::cout << v[i] << " ";//<< std::endl;
  }
  //std::cout << std::endl;
}



// partition permutations //
// ---------------------- //

struct Phi
{
  typedef std::pair<std::vector<int>, std::vector<int> > Sub_domains_pair;
  void put(const std::vector<int>& left, const std::vector<int>& right)
  {
    sub_domains_list.push_back(std::make_pair(left, right)); // preallocate list
  }

  std::size_t size() const
  {
    return sub_domains_list.size();
  }

  bool empty() const
  {
    return size() == 0 ? true : false;
  }

  // SL_comments: do you really want a copy?
  std::vector<int> lsubset(const int i) const
  {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < sub_domains_list.size());
    return sub_domains_list[i].first;
  }

  // SL_comments: do you really want a copy?
  std::vector<int> rsubset(const int i) const
  {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < sub_domains_list.size());
    return sub_domains_list[i].second;
  }

  std::vector<Sub_domains_pair> sub_domains_list;
};


void do_permutations(std::vector<std::vector<int> >& hole_list, Phi& subsets)
{
  if(hole_list.empty())
    return;

  std::vector<int> hs;

  for(int n = 0; n < hole_list.size(); ++n)
    hs.push_back(n);

  const int first = hs.front();
  const int last = hs.back();
  //std::sort(hs.begin(), hs.end()); // already sorted  ???

  for(int s = 0; s <= hs.size(); ++s) // s = number of holes on one (left) side
  {
    std::vector<int> p1(s);
    std::vector<int> p2(hole_list.size() - s);

    if(s == 0)
    {
      subsets.put(p1, hs);

      print(p1); std::cout << "-- "; print(hs); std::cout << std::endl;
      continue;
    }

    CGAL::Combination_enumerator<int> permutations(s, first, last+1);

    int p = 0;
    while(!permutations.finished())
    {
      // SL_comments: can't we use std::copy?
      for(int i=0; i<s; ++i)
      {
        p1[i] = permutations[i];
      }

      ++permutations;

      std::sort(p1.begin(), p1.end());
      std::set_symmetric_difference(p1.begin(), p1.end(), hs.begin(), hs.end(),
                                    p2.begin());

      CGAL_assertion(p1.size() == s);
      CGAL_assertion(p2.size() == hs.size() - s);

      print(p1); std::cout << "-- "; print(p2); std::cout << std::endl;

      subsets.put(p1, p2);
    }

  }
}

// split //
// ----- //

template <typename PointRange>
void split_domain(const Domain<PointRange>& init_domain,
                    Domain<PointRange>& left_dom, Domain<PointRange>& right_dom,
                    const int i, const int pid, const int k) // SL_comments: parameter doc + better name for variables
{
  typedef std::vector<int> Ids;
  Ids ids = init_domain.b_ids;
  Ids left;
  Ids right;
  const int n = ids.size();

  // i, k indices of access edge

  Ids::iterator it;
  it = find(ids.begin(), ids.end(), pid);
  CGAL_assertion(it != ids.end());

  // left subset
  // from pid to i
  left.push_back(*it);

  // assume: i is n , k is 0 at start
  while (*it != i) {

    if(it == ids.end()-1)
      it = ids.begin();
    else
      ++it;

    left.push_back(*it);
  }

  // right subset
  // from k to pid
  it = find(ids.begin(), ids.end(), k);
  CGAL_assertion(it != ids.end());

  right.push_back(*it);

  while (*it != pid) {

    if(it == ids.end()-1)
      it = ids.begin();
    else
      ++it;

    right.push_back(*it);
  }


  CGAL_assertion(left.front() == pid);
  CGAL_assertion(left.back() == i);
  CGAL_assertion(right.front() == k);
  CGAL_assertion(right.back() == pid);

  // SL_comments why not directly using left_dom and right_dom containers?
  left_dom.b_ids = left;
  right_dom.b_ids = right;

}

// SL_comments: doc
void reorder_island(std::vector<int>& h_ids, const int v)
{

  std::vector<int>::iterator it = find(h_ids.begin(), h_ids.end(), v);
  CGAL_assertion(it != h_ids.end());

  // 2) rotate by the third vertex of t
  std::size_t dist = std::distance(h_ids.begin(), it); // std::size_t?
  std::rotate(h_ids.begin(), h_ids.begin() + dist, h_ids.end());

  // 3) add the first removed element
  h_ids.push_back(h_ids[0]);

  // 4) reverse. Todo: Check and do it iff reversal is needed.
  std::reverse(h_ids.begin(), h_ids.end());

}

// SL_comments: doc
void merge_id_sets(std::vector<int>& b_ids,
                   const int i, const int v, const int k,
                   std::vector<int>& hole_ids)
{
  reorder_island(hole_ids, v);

  std::size_t initial_b_size = b_ids.size();

  // insertion point = just before k
  typename std::vector<int>::iterator insertion_point = b_ids.begin() + k;

  b_ids.insert(insertion_point, hole_ids.begin(), hole_ids.end());

  CGAL_assertion(*(b_ids.begin() + i) == b_ids[i] );
  CGAL_assertion(b_ids.size() == initial_b_size + hole_ids.size());

}

// SL_comments: doc
template<typename PointRange>
void join_domain(const Domain<PointRange>& domain, Domain<PointRange>& new_domain,
                  const int i, const int v, const int k)
{
  typedef std::vector<int> Ids;
  Ids id_set = domain.b_ids;
  Ids hole_ids = domain.h_ids; // for now assume just one hole.

  merge_id_sets(id_set, i, v, k, hole_ids);
  new_domain.b_ids = id_set;
}


struct Tracer
{
  template <class LookupTable>
  void
  operator()(const LookupTable& lambda, int v0, int v1)
  {
    CGAL_assertion_code( const int n = lambda.n; )
    std::stack<std::pair<int, int> > ranges;
    ranges.push(std::make_pair(v0, v1));

    while(!ranges.empty()) {
      std::pair<int, int> r = ranges.top();
      ranges.pop();
      CGAL_assertion(r.first >= 0 && r.first < n);
      CGAL_assertion(r.second >= 0 && r.second < n);

      // if on border
      if(r.first + 1 == r.second) { continue; }

      int la = lambda.get(r.first, r.second);
      if(la == -1) {
          std::cerr << "out hole" << std::endl;
          //*out_hole++ = std::make_pair(r.first, r.second);
        continue;
      }

      CGAL_assertion(la >= 0 && la < n);
      CGAL_assertion(r.first < la && r.second > la);
      auto triangle = std::make_tuple(r.first, la, r.second);

      collection.push_back(triangle);

      ranges.push(std::make_pair(r.first, la));
      ranges.push(std::make_pair(la, r.second));
    }
  }


  std::vector<std::tuple<int, int, int>> collection;
};



template<typename Kernel, typename WeightCalculator,
         typename WeightTable, typename LambdaTable>
class Triangulate
{
  typedef typename Kernel::Point_3 Point_3;
  typedef typename WeightCalculator::Weight Weight;
  typedef typename std::vector<Point_3> PointRange;


public:

  Triangulate(Domain<PointRange> domain,
              PointRange allpoints,
              WeightTable& W,
              LambdaTable& l,
              const WeightCalculator & WC) :
              Points(allpoints),
              W(W),
              lambda(l),
              domain(domain),
              WC(WC)
              {}

  std::size_t do_triangulation(const int i, const int k, std::size_t& count)
  {
    processDomain(domain, i, k, count);

    lambda.print("data/lambda-rec.dat");
    W.print("data/weight-rec.dat");
  }

  void collect_triangles(std::vector<std::tuple<int, int, int>>& triplets,
                         const int i, const int k)
  {
    Tracer tracer;
    tracer(lambda, i, k);
    triplets = tracer.collection;
  }


private:

  // main loop //
  // --------- //
  void processDomain(Domain<PointRange> domain, const int i, const int k, std::size_t& count)
  {
    // (i, k) = acccess edge

    // domains consisting of only one edge
    if(domain.b_ids.size() == 2)
      return;

    // base case
    if(domain.b_ids.size() == 3 && domain.holes_list.empty())
    {
      //CGAL_assertion(domain.b_ids[0] == i); // access edge source
      //CGAL_assertion(domain.b_ids[2] == k); // access edge target

      int m = domain.b_ids[1]; //third vertex
      std::cout<<"Evaluating t= ("<<i<<","<<m<<","<<k<<")"<<std::endl;
      calculate_weight(i, m, k);
      count++;

      return;
    }
    CGAL_assertion(domain.b_ids.size() >= 3);

    // pid : third vertex

    // CASE I - if there are islands, join until there are no islands.
    for(int pid : domain.h_ids)
    {
      //std::cout << "i= " << i << " k= " << k << std::endl;
      //std::cout << "pid= " << pid << std::endl;

      Domain<PointRange> D1;
      join_domain(domain, D1, i, pid, k);

      // get a new e_D
      std::pair<int, int> e_D1 = D1.get_access_edge();

      processDomain(D1, e_D1.first, e_D1.second, count);

      // calculate weight of triangle t - after the subdomains left and right have been checked
      int m = pid; //third vertex
      std::cout<<"Evaluating t= ("<<i<<","<<m<<","<<k<<")"<<std::endl;
      calculate_weight(i, m, k);
      count++;

    }

    // CASE II
    for(int pid : domain.b_ids)
    {
      //std::cout << "i= " << i << " k= " << k << std::endl;
      //std::cout << "pid= " << pid << std::endl;

      // avoid source & target of e_D
      if(pid == i || pid == k)
      {
        //std::cout << " point aborted" << std::endl;
        continue;
      }

      // split to two sub-domains
      Domain<PointRange> D1;
      Domain<PointRange> D2;
      // essentially splitting boundaries
      split_domain(domain, D1, D2, i, pid, k);
      // D1, D2 have just new boundaries - no hole information.

      // get new access edges for each
      std::pair<int, int> e_D1 = D1.get_access_edge();
      std::pair<int, int> e_D2 = D2.get_access_edge();

      // assign all combination of holes to subdomains and process each pair
      Phi partition_space;
      do_permutations(domain.holes_list, partition_space);

      if(partition_space.empty())
      {
        // when the domain has been merged so that there is no holes inside
        processDomain(D1, e_D1.first, e_D1.second, count);
        processDomain(D2, e_D2.first, e_D2.second, count);
      }
      else
      {
        // when t is formed with a vertex on the boundary of a domain with holes
        for(std::size_t p = 0; p < partition_space.size(); ++p)
        {
          std::vector<int> lholes = partition_space.lsubset(p);
          std::vector<int> rholes = partition_space.rsubset(p);

          for(int lh : lholes)
            D1.add_hole(domain.holes_list[lh]);

          for(int rh : rholes)
            D2.add_hole(domain.holes_list[rh]);

          processDomain(D1, e_D1.first, e_D1.second, count);
          processDomain(D2, e_D2.first, e_D2.second, count);
        }

      }


      // calculate weight of triangle t - after the subdomains left and right have been checked
      int m = pid; //third vertex
      std::cout<<"Evaluating t= ("<<i<<","<<m<<","<<k<<")"<<std::endl;
      calculate_weight(i, m, k);
      count++;

    }
  }


  void calculate_weight(const int i, const int m, const int k)
  {
    std::vector<Point_3> Q;

    // i, m, k are global indices

    const Weight& w_imk = WC(Points, Q, i,m,k, lambda);

    if(w_imk == Weight::NOT_VALID())
    {
      std::cerr << "non-manifold edge"  << std::endl;
      return;
    }

    auto weight_im = W.get(i,m);
    auto weight_mk = W.get(m,k);
    const Weight& w = weight_im + weight_mk + w_imk;

    if(lambda.get(i, k) == -1 || w < W.get(i, k)) {
      W.put(i,k,w);
      lambda.put(i,k, m);
    }
  }



  // data
  PointRange Points;

  WeightTable W;
  LambdaTable lambda;

  Domain<PointRange> domain;
  const WeightCalculator& WC;

};





} // namespace internal
} // namespace CGAL





#endif // CGAL_ISLAND_TRIANGULATE_HOLE_POLYLINE_H
