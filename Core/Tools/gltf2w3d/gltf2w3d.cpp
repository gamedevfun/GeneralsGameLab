/*
**      Command & Conquer Generals Zero Hour(tm)
**      Copyright 2025 Electronic Arts Inc.
**
**      This program is free software: you can redistribute it and/or modify
**      it under the terms of the GNU General Public License as published by
**      the Free Software Foundation, either version 3 of the License, or
**      (at your option) any later version.
**
**      This program is distributed in the hope that it will be useful,
**      but WITHOUT ANY WARRANTY; without even the implied warranty of
**      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**      GNU General Public License for more details.
**
**      You should have received a copy of the GNU General Public License
**      along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <cmath>

#include "tiny_gltf.h"
#include "json.hpp"
#include "stb_image.h"
#include "stb_image_write.h"

#include "../WW3D/max2w3d/w3d_file.h"
#include "../../Libraries/Source/WWVegas/WWLib/RAWFILE.H"
#include "../../Libraries/Source/WWVegas/WWLib/chunkio.h"

using namespace std;

static void compute_face(const W3dVectorStruct &a, const W3dVectorStruct &b,
                         const W3dVectorStruct &c, W3dVectorStruct &n, float &d)
{
    float ux = b.X - a.X;
    float uy = b.Y - a.Y;
    float uz = b.Z - a.Z;
    float vx = c.X - a.X;
    float vy = c.Y - a.Y;
    float vz = c.Z - a.Z;
    n.X = uy * vz - uz * vy;
    n.Y = uz * vx - ux * vz;
    n.Z = ux * vy - uy * vx;
    float len = sqrtf(n.X * n.X + n.Y * n.Y + n.Z * n.Z);
    if (len > 0.0f) {
        n.X /= len; n.Y /= len; n.Z /= len;
    }
    d = -(n.X * a.X + n.Y * a.Y + n.Z * a.Z);
}

static bool readAccessorFloats(const tinygltf::Model &model, int accessorIndex,
                               std::vector<float> &out)
{
    if (accessorIndex < 0 || accessorIndex >= (int)model.accessors.size()) return false;
    const tinygltf::Accessor &acc = model.accessors[accessorIndex];
    const tinygltf::BufferView &view = model.bufferViews[acc.bufferView];
    const tinygltf::Buffer &buf = model.buffers[view.buffer];
    const unsigned char *data = buf.data.data() + view.byteOffset + acc.byteOffset;
    size_t elemCount = acc.count * tinygltf::GetNumComponentsInType(acc.type);
    size_t elemSize = tinygltf::GetComponentSizeInBytes(acc.componentType);
    out.resize(elemCount);
    const float *src = reinterpret_cast<const float *>(data);
    for (size_t i = 0; i < elemCount; ++i) {
        out[i] = src[i];
    }
    return true;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        printf("Usage: %s <input.gltf> <output.w3d>\n", argv[0]);
        return 1;
    }

    string input = argv[1];
    string output = argv[2];

    tinygltf::TinyGLTF loader;
    tinygltf::Model model;
    string err, warn;
    bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, input);
    if (!warn.empty()) fprintf(stderr, "Warning: %s\n", warn.c_str());
    if (!err.empty()) fprintf(stderr, "Error: %s\n", err.c_str());
    if (!ret) {
        fprintf(stderr, "Failed to load %s\n", input.c_str());
        return 1;
    }
    if (model.meshes.empty()) {
        fprintf(stderr, "No mesh found in %s\n", input.c_str());
        return 1;
    }

    const tinygltf::Primitive &prim = model.meshes[0].primitives[0];

    vector<float> positions;
    vector<float> normals;
    vector<float> texcoords;
    vector<unsigned short> indices;

    if (!readAccessorFloats(model, prim.attributes.find("POSITION")->second, positions)) {
        fprintf(stderr, "Failed to read POSITION\n");
        return 1;
    }
    auto itN = prim.attributes.find("NORMAL");
    if (itN != prim.attributes.end()) {
        readAccessorFloats(model, itN->second, normals);
    }
    auto itT = prim.attributes.find("TEXCOORD_0");
    if (itT != prim.attributes.end()) {
        readAccessorFloats(model, itT->second, texcoords);
    }
    if (prim.indices >= 0) {
        const tinygltf::Accessor &acc = model.accessors[prim.indices];
        const tinygltf::BufferView &view = model.bufferViews[acc.bufferView];
        const tinygltf::Buffer &buf = model.buffers[view.buffer];
        const unsigned char *data = buf.data.data() + view.byteOffset + acc.byteOffset;
        size_t count = acc.count;
        indices.resize(count);
        if (acc.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
            memcpy(indices.data(), data, count * sizeof(unsigned short));
        } else if (acc.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) {
            const unsigned int *src = reinterpret_cast<const unsigned int *>(data);
            for (size_t i = 0; i < count; ++i) indices[i] = static_cast<unsigned short>(src[i]);
        } else {
            fprintf(stderr, "Unsupported index component type\n");
            return 1;
        }
    } else {
        size_t vertCount = positions.size() / 3;
        indices.resize(vertCount);
        for (size_t i = 0; i < vertCount; ++i) indices[i] = static_cast<unsigned short>(i);
    }

    size_t vertexCount = positions.size() / 3;
    vector<W3dVectorStruct> w3dVerts(vertexCount);
    vector<W3dVectorStruct> w3dNormals(vertexCount);
    vector<W3dTexCoordStruct> w3dTex(vertexCount);

    for (size_t i = 0; i < vertexCount; ++i) {
        w3dVerts[i].X = positions[i * 3 + 0];
        w3dVerts[i].Y = positions[i * 3 + 1];
        w3dVerts[i].Z = positions[i * 3 + 2];
        if (!normals.empty()) {
            w3dNormals[i].X = normals[i * 3 + 0];
            w3dNormals[i].Y = normals[i * 3 + 1];
            w3dNormals[i].Z = normals[i * 3 + 2];
        }
        if (!texcoords.empty()) {
            w3dTex[i].U = texcoords[i * 2 + 0];
            w3dTex[i].V = texcoords[i * 2 + 1];
        }
    }

    size_t triCount = indices.size() / 3;
    vector<W3dTriStruct> tris(triCount);
    for (size_t i = 0; i < triCount; ++i) {
        int i0 = indices[i * 3 + 0];
        int i1 = indices[i * 3 + 1];
        int i2 = indices[i * 3 + 2];
        tris[i].Vindex[0] = i0;
        tris[i].Vindex[1] = i1;
        tris[i].Vindex[2] = i2;
        tris[i].Attributes = 0;
        compute_face(w3dVerts[i0], w3dVerts[i1], w3dVerts[i2], tris[i].Normal, tris[i].Dist);
    }

    W3dVectorStruct minv = w3dVerts[0];
    W3dVectorStruct maxv = w3dVerts[0];
    for (size_t i = 1; i < vertexCount; ++i) {
        minv.X = std::min(minv.X, w3dVerts[i].X);
        minv.Y = std::min(minv.Y, w3dVerts[i].Y);
        minv.Z = std::min(minv.Z, w3dVerts[i].Z);
        maxv.X = std::max(maxv.X, w3dVerts[i].X);
        maxv.Y = std::max(maxv.Y, w3dVerts[i].Y);
        maxv.Z = std::max(maxv.Z, w3dVerts[i].Z);
    }
    W3dVectorStruct center;
    center.X = (minv.X + maxv.X) * 0.5f;
    center.Y = (minv.Y + maxv.Y) * 0.5f;
    center.Z = (minv.Z + maxv.Z) * 0.5f;
    float radius = 0.0f;
    for (size_t i = 0; i < vertexCount; ++i) {
        float dx = w3dVerts[i].X - center.X;
        float dy = w3dVerts[i].Y - center.Y;
        float dz = w3dVerts[i].Z - center.Z;
        float r2 = dx*dx + dy*dy + dz*dz;
        if (r2 > radius) radius = r2;
    }
    radius = sqrtf(radius);

    RawFileClass file(output.c_str());
    if (!file.Create()) {
        fprintf(stderr, "Failed to create %s\n", output.c_str());
        return 1;
    }
    ChunkSaveClass csave(&file);
    csave.Begin_Chunk(W3D_CHUNK_MESH);

    W3dMeshHeader3Struct header;
    memset(&header, 0, sizeof(header));
    header.Version = W3D_CURRENT_MESH_VERSION;
    strncpy(header.MeshName, "Mesh", W3D_NAME_LEN);
    strncpy(header.ContainerName, "Container", W3D_NAME_LEN);
    header.NumTris = (uint32)triCount;
    header.NumVertices = (uint32)vertexCount;
    header.NumMaterials = 0;
    header.NumDamageStages = 0;
    header.SortLevel = SORT_LEVEL_NONE;
    header.PrelitVersion = 0;
    header.VertexChannels = W3D_VERTEX_CHANNEL_LOCATION;
    if (!normals.empty()) header.VertexChannels |= W3D_VERTEX_CHANNEL_NORMAL;
    if (!texcoords.empty()) header.VertexChannels |= W3D_VERTEX_CHANNEL_TEXCOORD;
    header.FaceChannels = W3D_FACE_CHANNEL_FACE;
    header.Min = minv;
    header.Max = maxv;
    header.SphCenter = center;
    header.SphRadius = radius;

    csave.Begin_Chunk(W3D_CHUNK_MESH_HEADER3);
    csave.Write(&header, sizeof(header));
    csave.End_Chunk();

    csave.Begin_Chunk(W3D_CHUNK_VERTICES);
    csave.Write(w3dVerts.data(), (uint32)(w3dVerts.size() * sizeof(W3dVectorStruct)));
    csave.End_Chunk();

    if (!normals.empty()) {
        csave.Begin_Chunk(W3D_CHUNK_VERTEX_NORMALS);
        csave.Write(w3dNormals.data(), (uint32)(w3dNormals.size() * sizeof(W3dVectorStruct)));
        csave.End_Chunk();
    }
    if (!texcoords.empty()) {
        csave.Begin_Chunk(W3D_CHUNK_TEXCOORDS);
        csave.Write(w3dTex.data(), (uint32)(w3dTex.size() * sizeof(W3dTexCoordStruct)));
        csave.End_Chunk();
    }

    csave.Begin_Chunk(W3D_CHUNK_TRIANGLES);
    csave.Write(tris.data(), (uint32)(tris.size() * sizeof(W3dTriStruct)));
    csave.End_Chunk();

    csave.End_Chunk();
    file.Close();

    printf("W3D file written to %s\n", output.c_str());
    return 0;
}

