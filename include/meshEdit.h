#ifndef CMU462_MESHEDIT_H
#define CMU462_MESHEDIT_H

#include "halfEdgeMesh.h"

using namespace std;

namespace CMU462 {

class MeshResampler {
 public:
  static void upsample(HalfedgeMesh& mesh);
  static void downsample(HalfedgeMesh& mesh);
};

}  // namespace CMU462

#endif  // CMU462_MESHEDIT_H
