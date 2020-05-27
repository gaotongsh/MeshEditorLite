//
// Created by Tong Gao on 2020/5/17.
//

#ifndef OBJ_LOADER_HALFEDGE_H
#define OBJ_LOADER_HALFEDGE_H

#include <fstream>
#include <vector>

#include "halfEdgeMesh.h"

namespace CMU462 {

class ObjLoaderHalfedgeMesh {
public:
    static const int NUM_ATTR = 6;
    static const int STRIDE = 3 * NUM_ATTR;

    ObjLoaderHalfedgeMesh() = default;

    ~ObjLoaderHalfedgeMesh() {
        if (isData) {
            delete[] vertices;
            delete mesh;
        }
    }

    HalfedgeMesh* load(const char* filename) {
        std::ifstream input(filename);
        if (!input.is_open())
            return nullptr;

        std::string s;
        int nv;
        input >> s >> s >> nv;
        input >> s >> s >> nf;

        std::vector<Vector3D> vertexPositions;
        GLfloat x, y, z;
        for (size_t i = 0; i < nv; ++i) {
            input >> s >> x >> y >> z;
            Vector3D pos = {x, y, z};
            vertexPositions.push_back(pos);
        }

        std::vector<std::vector<Index>> polygons;
        std::string line;
        getline(input, line); // Remove current line
        int num;
        for (size_t i = 0; i < nf; ++i) {
            std::vector<Index> p;
            getline(input, line);
            istringstream iss(line);
            iss >> s;
            while(iss >> num) {
                p.push_back(num);
            }
            polygons.push_back(p);
        }

        input.close();

        mesh = new HalfedgeMesh();
        mesh->build(polygons, vertexPositions);
        mesh->triangulate();

        return mesh;
    }

    GLfloat* getData() {
        // Generate new vertices.
        nf = mesh->nFaces();
        if (isData) {
            delete[] vertices;
        } else {
            isData = true;
        }
        vertices = new GLfloat[STRIDE * nf];

        size_t i = 0;
        for (auto f = mesh->facesBegin(); f != mesh->facesEnd(); f++, i++) {
            auto pos1 = f->halfedge()->vertex()->position;
            auto pos2 = f->halfedge()->next()->vertex()->position;
            auto pos3 = f->halfedge()->next()->next()->vertex()->position;

            vertices[STRIDE * i    ]  = pos1.x;
            vertices[STRIDE * i + 1]  = pos1.y;
            vertices[STRIDE * i + 2]  = pos1.z;

            vertices[STRIDE * i + NUM_ATTR    ]  = pos2.x;
            vertices[STRIDE * i + NUM_ATTR + 1]  = pos2.y;
            vertices[STRIDE * i + NUM_ATTR + 2]  = pos2.z;

            vertices[STRIDE * i + NUM_ATTR * 2    ] = pos3.x;
            vertices[STRIDE * i + NUM_ATTR * 2 + 1] = pos3.y;
            vertices[STRIDE * i + NUM_ATTR * 2 + 2] = pos3.z;

            GLfloat color1 = static_cast<double>(rand()) / RAND_MAX;
            GLfloat color2 = static_cast<double>(rand()) / RAND_MAX;
            GLfloat color3 = static_cast<double>(rand()) / RAND_MAX;
            for (size_t j = 0; j < 3; ++j) {
                vertices[STRIDE * i + NUM_ATTR * j + 3] = color1;
                vertices[STRIDE * i + NUM_ATTR * j + 4] = color2;
                vertices[STRIDE * i + NUM_ATTR * j + 5] = color3;
            }
        }

        return vertices;
    }

    GLint getNf() const { return nf; }

private:
    bool isData;
    GLfloat *vertices;
    GLint nf;
    HalfedgeMesh *mesh;
};

} // namespace CMU462

#endif //OBJ_LOADER_HALFEDGE_H
