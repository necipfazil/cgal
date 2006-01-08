
#ifndef GPS_BFS_XOR_VISITOR_H
#define GPS_BFS_XOR_VISITOR_H

#include <CGAL/Boolean_set_operations_2/Gps_bfs_base_visitor.h>

CGAL_BEGIN_NAMESPACE

template <class Arrangement_>
class Gps_bfs_xor_visitor : public Gps_bfs_base_visitor<Arrangement_>
{
  typedef  Arrangement_                                  Arrangement;
  typedef typename Arrangement::Face_iterator            Face_iterator;
  typedef typename Arrangement::Halfedge_iterator        Halfedge_iterator;
  typedef Gps_bfs_base_visitor<Arrangement>              Base;
  typedef typename Base::Edges_hash                      Edges_hash;
  typedef typename Base::Faces_hash                      Faces_hash;

public:

  Gps_bfs_xor_visitor(Edges_hash* edges_hash, Faces_hash* faces_hash, unsigned int n_pgn): 
    Base(edges_hash, faces_hash, n_pgn)
  {}


  void flip_face(Face_iterator f1, Face_iterator f2, Halfedge_iterator he)
  {
    unsigned int ic_f2;
    ic_f2 = this->compute_ic(f1, f2, he);
    (*(this->m_faces_hash))[f2] = ic_f2;
      
    if(ic_f2%2)
      f2->set_contained(true);
  }

   // mark the unbounded_face (true iff contained)
  void visit_ubf(Face_iterator ubf, unsigned int ubf_ic)
  {
    CGAL_assertion(ubf->is_unbounded());
    if(ubf_ic%2)
      ubf->set_contained(true);
  }

  void after_scan(Arrangement& arr)
  {
    typedef typename Arrangement::Traits_2  Traits;
    typedef typename Traits::Compare_endpoints_xy_2 Compare_endpoints_xy_2;
    typedef typename Traits::Construct_opposite_2   Construct_opposite_2;
    typedef typename Arrangement::Edge_iterator     Edge_iterator;

    Compare_endpoints_xy_2 cmp_endpoints =
      m_traits->compare_endpoints_xy_2_object();
    Construct_opposite_2 ctr_opp = m_traits->construct_opposite_2_object();

    for(Edge_iterator eit = arr.edges_begin();
        eit != arr.edges_end();
        ++eit)
    {
      Halfedge_iterator he = eit;
      X_monotone_curve_2&  cv = he->curve();
      bool is_cont = he->face()->contained();
      bool has_same_dir = (cmp_endpoints(cv) == he->direction());
      if((is_cont && !has_same_dir) ||
         (!is_cont && has_same_dir))
        arr.modify_edge(he, ctr_opp(cv));
    }
  }

};

CGAL_END_NAMESPACE

#endif
