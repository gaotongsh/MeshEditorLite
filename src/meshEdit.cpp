#include <cassert>
#include "meshEdit.h"
#include "mutablePriorityQueue.h"
#include "matrix3x3.h"

namespace CMU462 {

VertexIter HalfedgeMesh::splitEdge(EdgeIter e0) {
  // This method should split the given edge and return an iterator to the
  // newly inserted vertex. The halfedge of this vertex should point along
  // the edge that was split, rather than the new edges.

  auto h = e0->halfedge(), ht = h->twin();
  auto f0 = h->face(), f1 = ht->face();
  if (f0->degree() != 3 || f1->degree() != 3) {
    return {};
  }

  // Collect
  // HalfEdges
  auto h1 = h->next(), h2 = h1->next(), ht1 = ht->next(), ht2 = ht1->next();
  // Vertices
  auto v0 = h->vertex(), v1 = ht->vertex(), v2 = h2->vertex(), v3 = ht2->vertex();
  // Faces - already collected

  // Allocate
  // HalfEdges
  auto h1_next = newHalfedge(), h1_next_twin = newHalfedge(),
       h2_next = newHalfedge(), h2_next_twin = newHalfedge(),
       ht_next = newHalfedge(), ht_next_twin = newHalfedge();
  // Vertex
  auto v_new = newVertex();
  // Edges
  auto e_h1_next = newEdge(), e_h2_next = newEdge(), e_ht_next = newEdge();
  // Faces
  auto f2 = newFace(), f3 = newFace();

  // Reassign
  // HalfEdges
  h1_next->setNeighbors(h, h1_next_twin, v2, e_h1_next, f0);
  h1->next() = h1_next;
  h->vertex() = v_new;

  h1_next_twin->setNeighbors(h2, h1_next, v_new, e_h1_next, f2);
  h2->next() = h2_next;
  h2->face() = f2;
  h2_next->setNeighbors(h1_next_twin, h2_next_twin, v0, e_h2_next, f2);

  h2_next_twin->setNeighbors(ht1, h2_next, v_new, e_h2_next, f3);
  ht1->next() = ht_next_twin;
  ht1->face() = f3;
  ht_next_twin->setNeighbors(h2_next_twin, ht_next, v3, e_ht_next, f3);

  ht_next->setNeighbors(ht2, ht_next_twin, v_new, e_ht_next, f1);
  ht->next() = ht_next;

  // Vertices
  v0->halfedge() = h2_next;
  v_new->halfedge() = h;
  v_new->position = (v0->position + v1->position) / 2;

  // Edges
  e_h1_next->halfedge() = h1_next;
  e_h2_next->halfedge() = h2_next;
  e_ht_next->halfedge() = ht_next;

  // Faces
  f0->halfedge() = h;
  f1->halfedge() = ht;
  f2->halfedge() = h2;
  f3->halfedge() = ht1;

  return v_new;
}

set<EdgeIter> HalfedgeMesh::remove_tetrahedron(EdgeIter e) {

  auto h = e->halfedge(), h_twin = h->twin();
  auto v_keep = h->vertex(), v_delete = h_twin->vertex();

  vector<VertexIter> v_keep_neighbor, v_delete_neighbor, v_intersection;
  set<EdgeIter> e_delete;
  do {
    v_keep_neighbor.clear();
    v_delete_neighbor.clear();
    v_intersection.clear();
    auto hi = v_keep->halfedge();
    do {
      v_keep_neighbor.push_back(hi->twin()->vertex());
      hi = hi->twin()->next();
    } while (hi != v_keep->halfedge());
    hi = v_delete->halfedge();
    do {
      v_delete_neighbor.push_back(hi->twin()->vertex());
      hi = hi->twin()->next();
    } while (hi != v_delete->halfedge());

    sort(v_keep_neighbor.begin(), v_keep_neighbor.end());
    sort(v_delete_neighbor.begin(), v_delete_neighbor.end());
    set_intersection(v_keep_neighbor.begin(), v_keep_neighbor.end(),
                     v_delete_neighbor.begin(), v_delete_neighbor.end(),
                     back_inserter(v_intersection));
    if (v_intersection.size() <= 2)
      break;

    cerr << "Tetrahedron spotted!" << endl;

    auto v_to_erase_begin = v_intersection.front();
    for (auto& vi : v_intersection) {
      if (vi->degree() < v_to_erase_begin->degree())
        v_to_erase_begin = vi;
    }

    vector<VertexIter> v_to_erase;
    v_to_erase.push_back(v_to_erase_begin);

    for (Index i = 0; i < v_to_erase.size(); ++i) {
      auto hj = v_to_erase[i]->halfedge();
      do {
        e_delete.insert(hj->edge());

        auto v_neighbor = hj->twin()->vertex();
        if (v_neighbor != v_delete && v_neighbor != v_keep
            && find(v_intersection.begin(), v_intersection.end(), v_neighbor) == v_intersection.end()
            && find(v_to_erase.begin(), v_to_erase.end(), v_neighbor) == v_to_erase.end()) {
          v_to_erase.push_back(v_neighbor);
        }

        hj = hj->twin()->next();
      } while (hj != v_to_erase[i]->halfedge());
    }

    for (auto& vi : v_to_erase) {
      eraseVertex(vi);
    }

  } while (true);

  return e_delete;
}

VertexIter HalfedgeMesh::collapseEdge(EdgeIter e) {
  // This method should collapse the given edge and return an iterator to
  // the new vertex created by the collapse.

  remove_tetrahedron(e);

  auto h = e->halfedge(), h_twin = h->twin();
  auto v_keep = h->vertex(), v_delete = h_twin->vertex();

  // Collect
  // HalfEdges
  vector<HalfedgeIter> h_to_be_deleted;
  h_to_be_deleted.push_back(h);
  h_to_be_deleted.push_back(h_twin);

  vector<HalfedgeIter> h_delete_vec;
  auto hi = h_twin;
  while (hi->twin()->next() != h_twin) {
    hi = hi->twin()->next();
    h_delete_vec.push_back(hi);
  }
  auto h_0 = h_delete_vec.front(), h_0_twin = h_0->twin(),
       h_n = h_delete_vec.back(), h_n_twin = h_n->twin(),
       h_1 = h_0_twin->next();

  auto h_first = h;
  while (h_first->next() != h) {
    h_first = h_first->next();
  }
  auto h_last = h_twin->next();
  auto h_first_prev = h_1;
  while (h_first_prev->next() != h_0_twin) {
    h_first_prev = h_first_prev->next();
  }
  auto h_last_prev = h_n->next();
  while (h_last_prev->next() != h_n) {
    h_last_prev = h_last_prev->next();
  }
  auto h_last_next = h_n->next();

  // Vertices
  v_keep->position = (v_keep->position + v_delete->position) / 2;
  auto v_0 = h_0_twin->vertex(), v_n = h_n_twin->vertex();

  // Edges
  vector<EdgeIter> e_to_be_deleted;
  e_to_be_deleted.push_back(e);
  auto e_0 = h_0->edge(), e_n = h_n->edge();

  // Faces
  vector<FaceIter> f_to_be_deleted;
  auto f_0 = h->face(),
       f_1 = h_0_twin->face(),
       f_n = h_n->face(),
       f_n_p1 = h_twin->face();

  // Reassign
  if (h_0->next() == h_first) {
    h_to_be_deleted.push_back(h_0);
    h_to_be_deleted.push_back(h_0_twin);
    e_to_be_deleted.push_back(e_0);
    f_to_be_deleted.push_back(f_0);

    // HalfEdges
    // We have a corner case here: when the faces on either side of the edge have a degree of 3,
    // while the vertex to be deleted also has a degree of 3. We have to deal with this carefully.
    if (h_last->next() == h_n_twin && f_1 == f_n) {
      h_first->next() = h_last;
    } else {
      h_first->next() = h_1;
    }
    h_first->face() = f_1;
    h_first_prev->next() = h_first;

    // Vertices
    v_0->halfedge() = h_first;

    // Faces
    f_1->halfedge() = h_first;
  } else {
    h_first->next() = h_0;
    f_0->halfedge() = h_0;
  }
  for (auto& hj : h_delete_vec) {
    hj->vertex() = v_keep;
  }
  if (h_last->next() == h_n_twin) {
    h_to_be_deleted.push_back(h_n);
    h_to_be_deleted.push_back(h_n_twin);
    e_to_be_deleted.push_back(e_n);
    f_to_be_deleted.push_back(f_n_p1);

    // HalfEdges
    h_last->next() = h_last_next;
    h_last->face() = f_n;
    h_last_prev->next() = h_last;

    // Vertices
    v_n->halfedge() = h_last_next;

    // Faces
    f_n->halfedge() = h_last;
  } else {
    h_n_twin->next() = h_last;
    f_n_p1->halfedge() = h_last;
  }
  v_keep->halfedge() = h_last;

  // Delete
  for (auto& hj : h_to_be_deleted) {
    deleteHalfedge(hj);
  }
  deleteVertex(v_delete);
  for (auto& ej : e_to_be_deleted) {
    deleteEdge(ej);
  }
  for (auto& fj : f_to_be_deleted) {
    deleteFace(fj);
  }

  return v_keep;
}

FaceIter HalfedgeMesh::eraseVertex(VertexIter v) {
  // This method should replace the given vertex and all its neighboring
  // edges and faces with a single face, returning the new face.

  // Collect
  // HalfEdges
  auto h = v->halfedge(), ht = h->twin();
  vector<HalfedgeIter> h_to_be_deleted, h_outer_first, h_outer_last, h_outer_all;
  h_to_be_deleted.push_back(h);
  h_to_be_deleted.push_back(ht);
  h_outer_first.push_back(h->next());

  // Vertices
  vector<VertexIter> v_outer_vec;
  v_outer_vec.push_back(ht->vertex());

  // Edges
  vector<EdgeIter> e_to_be_deleted;
  e_to_be_deleted.push_back(h->edge());

  // Faces
  auto f0 = h->face();
  vector<FaceIter> f_to_be_deleted;

  // Iterations through remaining
  auto hp = ht, hi = hp->next();
  while (hi != h) {
    h_to_be_deleted.push_back(hi);
    h_to_be_deleted.push_back(hi->twin());
    h_outer_first.push_back(hi->next());

    v_outer_vec.push_back(hi->twin()->vertex());
    e_to_be_deleted.push_back(hi->edge());
    f_to_be_deleted.push_back(hi->face());

    auto hj = hi;
    while (hj->next() != hp) {
      hj = hj->next();
      h_outer_all.push_back(hj);
    }
    h_outer_last.push_back(hj);

    hp = hi->twin();
    hi = hp->next();
  }
  while (hi->next() != hp) {
    hi = hi->next();
    h_outer_all.push_back(hi);
  }
  h_outer_last.push_back(hi);

  // Reassign
  Size n = h_outer_first.size();
  for (Index i = 0; i < n; ++i) {
    h_outer_last[i]->next() = h_outer_first[i];
    v_outer_vec[i]->halfedge() = h_outer_first[i];
  }
  for (auto& hj : h_outer_all) {
    hj->face() = f0;
  }
  f0->halfedge() = h_outer_first.front();

  // Delete
  for (auto& hj : h_to_be_deleted) {
    deleteHalfedge(hj);
  }
  deleteVertex(v);
  for (auto& ej : e_to_be_deleted) {
    deleteEdge(ej);
  }
  for (auto& fj : f_to_be_deleted) {
    deleteFace(fj);
  }

  return f0;
}

EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0) {
  // This method should flip the given edge and return an iterator to the
  // flipped edge.

  // Collect
  // HalfEdges
  auto h = e0->halfedge(), ht = h->twin(),
       h1 = h->next(), h2 = h1->next(),
       ht1 = ht->next(), ht2 = ht1->next(),
       hn = h2, htn = ht2;
  while (hn->next() != h) {
    hn = hn->next();
  }
  while (htn->next() != ht) {
    htn = htn->next();
  }

  // Vertices
  auto v = h->vertex(), vt = ht->vertex(),
       vt1 = h1->twin()->vertex(), v1 = ht1->twin()->vertex();

  // Faces
  auto f = h->face(), ft = ht->face();

  // Reassign
  // HalfEdges
  h->setNeighbors(h2, h->twin(), v1, h->edge(), h->face());
  h1->setNeighbors(ht, h1->twin(), h1->vertex(), h1->edge(), ft);
  hn->next() = ht1;

  ht->setNeighbors(ht2, ht->twin(), vt1, ht->edge(), ht->face());
  ht1->setNeighbors(h, ht1->twin(), ht1->vertex(), ht1->edge(), f);
  htn->next() = h1;

  // Vertices
  v->halfedge() = ht1;
  vt->halfedge() = h1;

  // Faces
  f->halfedge() = h;
  ft->halfedge() = ht;

  return e0;
}

void HalfedgeMesh::splitPolygon(FaceIter f) {
  // Triangulate a polygonal face

  // If already triangle, return
  if (f->degree() == 3) {
    return;
  }

  // Collect
  // HalfEdges
  auto h = f->halfedge();
  size_t degree = 0;
  vector<HalfedgeIter> h_outer_vec;
  do {
    ++degree;
    h_outer_vec.push_back(h);
    h = h->next();
  } while (h != f->halfedge());

  // Vertices
  vector<VertexIter> v_outer_vec;
  v_outer_vec.reserve(degree);
  for (const auto& hi : h_outer_vec) {
    v_outer_vec.push_back(hi->vertex());
  }

  // Allocate
  // HalfEdges
  vector<HalfedgeIter> h_new_vec, h_new_twin_vec;
  h_new_vec.reserve(degree - 3);
  h_new_twin_vec.reserve(degree - 3);

  // Edges
  vector<EdgeIter> e_new_vec;
  e_new_vec.reserve(degree - 3);

  // Faces
  vector<FaceIter> f_new_vec;
  f_new_vec.reserve(degree - 3);

  for (size_t i = 0; i < degree - 3; ++i) {
    h_new_vec.push_back(newHalfedge());
    h_new_twin_vec.push_back(newHalfedge());

    e_new_vec.push_back(newEdge());

    f_new_vec.push_back(newFace());
  }

  // Reassign
  // New edges and faces
  for (size_t i = 0; i < degree - 3; ++i) {
    e_new_vec[i]->halfedge() = h_new_vec[i];
    f_new_vec[i]->halfedge() = h_new_twin_vec[i];
  }

  // Halfedges of the first face
  h_new_vec[0]->setNeighbors(h_outer_vec[0],
                             h_new_twin_vec[0],
                             v_outer_vec[2],
                             e_new_vec[0],
                             f);
  h_outer_vec[1]->next() = h_new_vec[0];

  // Halfedges of the middle faces
  size_t first_v = 0, second_v = 2, third_v = degree - 1, t;
  for (size_t i = 0; i < degree - 4; ++i) {
    if (i % 2 == 0) {
      h_new_twin_vec[i]->setNeighbors(h_new_vec[i + 1],
                                      h_new_vec[i],
                                      v_outer_vec[first_v],
                                      e_new_vec[i],
                                      f_new_vec[i]);
      h_new_vec[i + 1]->setNeighbors(h_outer_vec[third_v],
                                     h_new_twin_vec[i + 1],
                                     v_outer_vec[second_v],
                                     e_new_vec[i + 1],
                                     f_new_vec[i]);
      h_outer_vec[third_v]->next() = h_new_twin_vec[i];
      h_outer_vec[third_v]->face() = f_new_vec[i];
      t = second_v + 1;
      first_v = second_v;
      second_v = third_v;
      third_v = t;
    } else {
      h_new_vec[i + 1]->setNeighbors(h_new_twin_vec[i],
                                     h_new_twin_vec[i + 1],
                                     v_outer_vec[third_v],
                                     e_new_vec[i + 1],
                                     f_new_vec[i]);
      h_new_twin_vec[i]->setNeighbors(h_outer_vec[first_v],
                                      h_new_vec[i],
                                      v_outer_vec[second_v],
                                      e_new_vec[i],
                                      f_new_vec[i]);
      h_outer_vec[first_v]->next() = h_new_vec[i + 1];
      h_outer_vec[first_v]->face() = f_new_vec[i];
      t = second_v - 1;
      first_v = second_v;
      second_v = third_v;
      third_v = t;
    }
  }

  // Halfedges of the final faces
  if (degree % 2 == 0) {
    h_new_twin_vec.back()->setNeighbors(h_outer_vec[second_v],
                                        h_new_vec.back(),
                                        v_outer_vec[first_v],
                                        e_new_vec.back(),
                                        f_new_vec.back());
    h_outer_vec[second_v]->face() = f_new_vec.back();
    h_outer_vec[third_v]->face() = f_new_vec.back();
    h_outer_vec[third_v]->next() = h_new_twin_vec.back();
  } else {
    h_new_twin_vec.back()->setNeighbors(h_outer_vec[first_v],
                                        h_new_vec.back(),
                                        v_outer_vec[second_v],
                                        e_new_vec.back(),
                                        f_new_vec.back());
    h_outer_vec[first_v]->face() = f_new_vec.back();
    h_outer_vec[third_v]->face() = f_new_vec.back();
    h_outer_vec[third_v]->next() = h_new_twin_vec.back();
  }
}

EdgeRecord::EdgeRecord(EdgeIter& _edge) : edge(_edge) {
  // Compute the combined quadric from the edge endpoints.

  auto h = edge->halfedge();
  auto K = h->vertex()->quadric + h->twin()->vertex()->quadric;

  // -> Build the 3x3 linear system whose solution minimizes the quadric error
  //    associated with these two endpoints.

  double A_data[9] = { K(0,0), K(0,1), K(0,2),
                       K(1,0), K(1,1), K(1,2),
                       K(2,0), K(2,1), K(2,2) };
  auto A = Matrix3x3(A_data);
  auto b = Vector3D(-K(0,3), -K(1,3), -K(2,3));

  // -> Use this system to solve for the optimal position, and store it in
  //    EdgeRecord::optimalPoint.

  if (A.det() > 0.01) {
    optimalPoint = A.inv() * b;
  } else {
//    cerr << "detA is too small: " << A.det() << endl;
//    cerr << "inverse answer: " << A.inv() * b << endl;

    auto p0 = h->vertex()->position, p1 = h->twin()->vertex()->position;
    while ((p0 - p1).norm2() > 0.0001) {
      auto c0 = dot(p0, A * p0) - 2 * dot(b, p0),
           c1 = dot(p1, A * p1) - 2 * dot(b, p1);
      if (c0 > c1) {
        p0 = (p0 + p1) / 2;
      } else {
        p1 = (p0 + p1) / 2;
      }
    }
    optimalPoint = p0;

//    cerr << "actual answer: " << optimalPoint << endl << endl;
  }

  // -> Also store the cost associated with collapsing this edg in
  //    EdgeRecord::Cost.

  score = dot(optimalPoint, A * optimalPoint) - 2 * dot(b, optimalPoint) + K(3,3);
}

void MeshResampler::upsample(HalfedgeMesh& mesh)
// This routine should increase the number of triangles in the mesh using Loop
// subdivision.
{
  // Compute new positions for all the vertices in the input mesh, using
  // the Loop subdivision rule, and store them in Vertex::newPosition.
  for (auto v = mesh.verticesBegin(); v != mesh.verticesEnd(); ++v) {
    Vector3D pos;
    auto h = v->halfedge()->twin(), h_start = h;
    Size degree = 0;
    do {
      ++degree;
      pos += h->vertex()->position;
      h = h->next()->twin();
    } while (h != h_start);
    double u = (degree == 3) ? (3.0 / 16) : (3.0 / (8.0 * degree));
    v->newPosition = u * pos + (1 - degree * u) * v->position;
  // -> At this point, we also want to mark each vertex as being a vertex of the
  //    original mesh.
    v->isNew = false;
  }

  // -> Next, compute the updated vertex positions associated with edges, and
  //    store it in Edge::newPosition.
  for (auto e = mesh.edgesBegin(); e != mesh.edgesEnd(); ++e) {
    Vector3D c;
    auto h = e->halfedge();
    c += (3.0 / 8) * h->vertex()->position;
    c += (3.0 / 8) * h->twin()->vertex()->position;
    c += (1.0 / 8) * h->next()->next()->vertex()->position;
    c += (1.0 / 8) * h->twin()->next()->next()->vertex()->position;
    e->newPosition = c;
  }

  // -> Next, we're going to split every edge in the mesh, in any order.  For
  //    future reference, we're also going to store some information about which
  //    subdivided edges come from splitting an edge in the original mesh, and
  //    which edges are new, by setting the flat Edge::isNew. Note that in this
  //    loop, we only want to iterate over edges of the original mesh.
  //    Otherwise, we'll end up splitting edges that we just split (and the
  //    loop will never end!)
  Size n = mesh.nEdges();
  auto e = mesh.edgesBegin();
  for (Index i = 0; i < n; ++i) {
    // get the next edge NOW!
    auto nextEdge = e;
    nextEdge++;

    // now, even if splitting the edge deletes it...
    auto new_pos = e->newPosition;
    auto new_v = mesh.splitEdge(e);
    new_v->isNew = true;
    new_v->halfedge()->edge()->isNew = false;
    new_v->halfedge()->next()->next()->edge()->isNew = true;
    new_v->halfedge()->twin()->next()->edge()->isNew = true;
    new_v->halfedge()->next()->next()->twin()->next()->next()->edge()->isNew = false;
    new_v->newPosition = new_pos;

    // ...we still have a valid reference to the next edge.
    e = nextEdge;
  }

  // -> Now flip any new edge that connects an old and new vertex.
  while (e != mesh.edgesEnd()) {
    auto nextEdge = e;
    nextEdge++;

    auto h = e->halfedge();
    if (e->isNew && (h->vertex()->isNew != h->twin()->vertex()->isNew)) {
      mesh.flipEdge(e);
    }

    e = nextEdge;
  }

  // -> Finally, copy the new vertex positions into final Vertex::position.
  for (auto v = mesh.verticesBegin(); v != mesh.verticesEnd(); ++v) {
    v->position = v->newPosition;
  }
}

void MeshResampler::downsample(HalfedgeMesh& mesh) {
  // Compute initial quadrics for each face by simply writing the plane equation
  // for the face in homogeneous coordinates. These quadrics should be stored
  // in Face::quadric

  for (auto f = mesh.facesBegin(); f != mesh.facesEnd(); ++f) {
    auto n = f->normal(), p = f->halfedge()->vertex()->position;
    auto v = Vector4D(n.x, n.y, n.z, -dot(n, p));
    f->quadric = outer(v, v);
  }

  // -> Compute an initial quadric for each vertex as the sum of the quadrics
  //    associated with the incident faces, storing it in Vertex::quadric

  for (auto v = mesh.verticesBegin(); v != mesh.verticesEnd(); ++v) {
    auto q = Matrix4x4();
    auto h = v->halfedge();
    do {
      q += h->face()->quadric;
      h = h->twin()->next();
    } while (h != v->halfedge());
    v->quadric = q;
  }

  // -> Build a priority queue of edges according to their quadric error cost,
  //    i.e., by building an EdgeRecord for each edge and sticking it in the
  //    queue.

  MutablePriorityQueue<EdgeRecord> queue;
  for (auto e = mesh.edgesBegin(); e != mesh.edgesEnd(); ++e) {
    e->record = EdgeRecord(e);
    queue.insert(e->record);
  }

  // -> Until we reach the target edge budget, collapse the best edge. Remember
  //    to remove from the queue any edge that touches the collapsing edge
  //    BEFORE it gets collapsed, and add back into the queue any edge touching
  //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
  //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
  //    top of the queue.

  Size target = mesh.nFaces() / 4;
  while (mesh.nFaces() > target) {
    // 1. Get the cheapest edge from the queue.
    auto er = queue.top();
    // 2. Remove the cheapest edge from the queue by calling pop().
    queue.pop();
    // 3. Compute the new quadric by summing the quadrics at its two endpoints.
    auto h = er.edge->halfedge();
    auto v0 = h->vertex(), v1 = h->twin()->vertex();
    auto newq = v0->quadric + v1->quadric;
    // 4. Remove any edge touching either of its endpoints from the queue.
    auto hi = v0->halfedge();
    do {
      queue.remove(hi->edge()->record);
      hi = hi->twin()->next();
    } while (hi != v0->halfedge());
    hi = v1->halfedge();
    do {
      queue.remove(hi->edge()->record);
      hi = hi->twin()->next();
    } while (hi != v1->halfedge());

    // [ADDITIONAL] Remove the tetrahedron edges.
    auto e_deleted = mesh.remove_tetrahedron(er.edge);
    for (auto& e : e_deleted) {
      queue.remove(e->record);
    }

    // 5. Collapse the edge.
    auto newv = mesh.collapseEdge(er.edge);
    newv->position = er.optimalPoint;
    // 6. Set the quadric of the new vertex to the quadric computed in Step 3.
    newv->quadric = newq;
    // 7. Insert any edge touching the new vertex into the queue, creating new edge records for each of them.
    hi = newv->halfedge();
    do {
      auto e = hi->edge();
      e->record = EdgeRecord(e);
      queue.insert(e->record);
      hi = hi->twin()->next();
    } while (hi != newv->halfedge());
  }
}

}  // namespace CMU462
