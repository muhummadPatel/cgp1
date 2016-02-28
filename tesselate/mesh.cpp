//
// Mesh
//

#include "mesh.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <list>
#include <sys/stat.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/intersect.hpp>

#include <unordered_map>
#include <algorithm>
#include <utility> //std::pair
#include <boost/functional/hash.hpp>

using namespace std;
using namespace cgp;

GLfloat stdCol[] = {0.7f, 0.7f, 0.75f, 0.4f};
const int raysamples = 5;

bool Sphere::pointInSphere(cgp::Point pnt)
{
    cgp::Vector delvec;

    delvec.diff(c, pnt);
    if(delvec.sqrdlength() < r*r)
        return true;
    else
        return false;
}

bool Mesh::findVert(cgp::Point pnt, int &idx)
{
    bool found = false;
    int i = 0;

    idx = -1;
    // linear search of vertex list
    while(!found && i < (int) verts.size())
    {
        if(verts[i] == pnt)
        {
            found = true;
            idx = i;
        }
        i++;
    }
    return found;
}

bool Mesh::findEdge(vector<Edge> edges, Edge e, int &idx)
{
    bool found = false;
    int i = 0;

    idx = -1;
    // linear search of edge list
    while(!found && i < (int) edges.size())
    {
        if( (edges[i].v[0] == e.v[0] && edges[i].v[1] == e.v[1]) || (edges[i].v[1] == e.v[0] && edges[i].v[0] == e.v[1]) )
        {
            found = true;
            idx = i;
        }
        i++;
    }
    return found;
}

long Mesh::hashVert(cgp::Point pnt, cgp::BoundBox bbox)
{
    long x, y, z;
    float range = 2500.0f;
    long lrangesq, lrange = 2500;

    lrangesq = lrange * lrange;

    // discretise vertex within bounds of the enclosing bounding box
    x = (long) (((pnt.x - bbox.min.x) * range) / bbox.diagLen()) * lrangesq;
    y = (long) (((pnt.y - bbox.min.y) * range) / bbox.diagLen()) * lrange;
    z = (long) (((pnt.z - bbox.min.z) * range) / bbox.diagLen());
    return x+y+z;
}

void Mesh::mergeVerts()
{
    vector<cgp::Point> cleanverts;
    long key;
    int i, p, hitcount = 0;
    // use hashmap to quickly look up vertices with the same coordinates
    std::unordered_map<long, int> idxlookup; // key is concatenation of vertex position, value is index into the cleanverts vector
    cgp::BoundBox bbox;

    // construct a bounding box enclosing all vertices
    for(i = 0; i < (int) verts.size(); i++)
        bbox.includePnt(verts[i]);

    // remove duplicate vertices
    for(i = 0; i < (int) verts.size(); i++)
    {
        key = hashVert(verts[i], bbox);
        if(idxlookup.find(key) == idxlookup.end()) // key not in map
        {
            idxlookup[key] = (int) cleanverts.size(); // put index in map for quick lookup
            cleanverts.push_back(verts[i]);
        }
        else
        {
            hitcount++;
        }
    }
    cerr << "num duplicate vertices found = " << hitcount << " of " << (int) verts.size() << endl;
    cerr << "clean verts = " << (int) cleanverts.size() << endl;
    cerr << "bbox min = " << bbox.min.x << ", " << bbox.min.y << ", " << bbox.min.z << endl;
    cerr << "bbox max = " << bbox.max.x << ", " << bbox.max.y << ", " << bbox.max.z << endl;
    cerr << "bbox diag = " << bbox.diagLen() << endl;


    // re-index triangles
    for(i = 0; i < (int) tris.size(); i++)
        for(p = 0; p < 3; p++)
        {
            key = hashVert(verts[tris[i].v[p]], bbox);
            if(idxlookup.find(key) != idxlookup.end())
                tris[i].v[p] = idxlookup[key];
            else
                cerr << "Error Mesh::mergeVerts: vertex not found in map" << endl;

        }

    verts.clear();
    verts = cleanverts;
}

void Mesh::deriveVertNorms()
{
    vector<int> vinc; // number of faces incident on vertex
    int p, t;
    cgp::Vector n;

    // init structures
    for(p = 0; p < (int) verts.size(); p++)
    {
        vinc.push_back(0);
        norms.push_back(cgp::Vector(0.0f, 0.0f, 0.0f));
    }

    // accumulate face normals into vertex normals
    for(t = 0; t < (int) tris.size(); t++)
    {
        n = tris[t].n; n.normalize();
        for(p = 0; p < 3; p++)
        {
            norms[tris[t].v[p]].add(n);
            vinc[tris[t].v[p]]++;
        }
    }

    // complete average
    for(p = 0; p < (int) verts.size(); p++)
    {
        norms[p].mult(1.0f/((float) vinc[p]));
        norms[p].normalize();
    }
}

void Mesh::deriveFaceNorms()
{
    int t;
    cgp::Vector evec[2];

    for(t = 0; t < (int) tris.size(); t++)
    {
        // right-hand rule for calculating normals, i.e. counter-clockwise winding from front on vertices
        evec[0].diff(verts[tris[t].v[0]], verts[tris[t].v[1]]);
        evec[1].diff(verts[tris[t].v[0]], verts[tris[t].v[2]]);
        evec[0].normalize();
        evec[1].normalize();
        tris[t].n.cross(evec[0], evec[1]);
        tris[t].n.normalize();
    }

}

void Mesh::buildTransform(glm::mat4x4 &tfm)
{
    glm::mat4x4 idt;

    idt = glm::mat4(1.0f);
    tfm = glm::translate(idt, glm::vec3(trx.i, trx.j, trx.k));
    tfm = glm::rotate(tfm, zrot, glm::vec3(0.0f, 0.0f, 1.0f));
    tfm = glm::rotate(tfm, yrot, glm::vec3(0.0f, 1.0f, 0.0f));
    tfm = glm::rotate(tfm, xrot, glm::vec3(1.0f, 0.0f, 0.0f));
    tfm = glm::scale(tfm, glm::vec3(scale));
}

Mesh::Mesh()
{
    col = stdCol;
    scale = 1.0f;
    xrot = yrot = zrot = 0.0f;
    trx = cgp::Vector(0.0f, 0.0f, 0.0f);
}

Mesh::~Mesh()
{
    clear();
}

void Mesh::clear()
{
    verts.clear();
    tris.clear();
    geom.clear();
    col = stdCol;
    scale = 1.0f;
    xrot = yrot = zrot = 0.0f;
    trx = cgp::Vector(0.0f, 0.0f, 0.0f);
    for(int i = 0; i < (int) boundspheres.size(); i++)
        boundspheres[i].ind.clear();
    boundspheres.clear();

    edges.clear();
}

bool Mesh::genGeometry(View * view, ShapeDrawData &sdd)
{
    vector<int> faces;
    int t, p;
    glm::mat4x4 tfm;

    geom.clear();
    geom.setColour(col);

    // transform mesh data structures into a form suitable for rendering
    // by flattening the triangle list
    for(t = 0; t < (int) tris.size(); t++)
        for(p = 0; p < 3; p++)
            faces.push_back(tris[t].v[p]);

    // construct transformation matrix
    buildTransform(tfm);
    geom.genMesh(&verts, &norms, &faces, tfm);

    /*
    // generate geometry for acceleration spheres for testing
    for(p = 0; p < (int) boundspheres.size(); p++)
    {
        glm::mat4x4 idt;

        idt = glm::mat4(1.0f);
        tfm = glm::translate(idt, glm::vec3(trx.i+boundspheres[p].c.x, trx.j+boundspheres[p].c.y, trx.k+boundspheres[p].c.z));
        tfm = glm::rotate(tfm, zrot, glm::vec3(0.0f, 0.0f, 1.0f));
        tfm = glm::rotate(tfm, yrot, glm::vec3(0.0f, 1.0f, 0.0f));
        tfm = glm::rotate(tfm, xrot, glm::vec3(1.0f, 0.0f, 0.0f));
        tfm = glm::scale(tfm, glm::vec3(scale));
        geom.genSphere(boundspheres[p].r, 20, 20, tfm);
    }*/

    // bind geometry to buffers and return drawing parameters, if possible
    if(geom.bindBuffers(view))
    {
        sdd = geom.getDrawParameters();
        return true;
    }
    else
       return false;
}

void Mesh::boxFit(float sidelen)
{
    cgp::Point pnt;
    cgp::Vector shift, diag, halfdiag;
    float scale;
    int v;
    cgp::BoundBox bbox;

    // calculate current bounding box
    for(v = 0; v < (int) verts.size(); v++)
        bbox.includePnt(verts[v]);

    if((int) verts.size() > 0)
    {
        // calculate translation necessary to move center of bounding box to the origin
        diag = bbox.getDiag();
        shift.pntconvert(bbox.min);
        halfdiag = diag; halfdiag.mult(0.5f);
        shift.add(halfdiag);
        shift.mult(-1.0f);

        // scale so that largest side of bounding box fits sidelen
        scale = max(diag.i, diag.j); scale = max(scale, diag.k);
        scale = sidelen / scale;

        // shift center to origin and scale uniformly
        for(v = 0; v < (int) verts.size(); v++)
        {
            pnt = verts[v];
            shift.pntplusvec(pnt, &pnt);
            pnt.x *= scale; pnt.y *= scale; pnt.z *= scale;
            verts[v] = pnt;
        }
    }
}

bool Mesh::readSTL(string filename)
{
    ifstream infile;
    char * inbuffer;
    struct stat results;
    int insize, inpos, numt, t, i;
    cgp::Point vpos;
    Triangle tri;

    // assumes binary format STL file
    infile.open((char *) filename.c_str(), ios_base::in | ios_base::binary);
    if(infile.is_open())
    {
        clear();

        // get the size of the file
        stat((char *) filename.c_str(), &results);
        insize = results.st_size;

        // put file contents in buffer
        inbuffer = new char[insize];
        infile.read(inbuffer, insize);
        if(!infile) // failed to read from the file for some reason
        {
            cerr << "Error Mesh::readSTL: unable to populate read buffer" << endl;
            return false;
        }

        // interpret buffer as STL file
        if(insize <= 84)
        {
            cerr << "Error Mesh::readSTL: invalid STL binary file, too small" << endl;
            return false;
        }

        inpos = 80; // skip 80 character header
        if(inpos+4 >= insize){ cerr << "Error Mesh::readSTL: malformed header on stl file" << endl; return false; }
        numt = (int) (* ((long *) &inbuffer[inpos]));
        inpos += 4;

        t = 0;

        // triangle vertices have consistent outward facing clockwise winding (right hand rule)
        while(t < numt) // read in triangle data
        {
            // normal
            if(inpos+12 >= insize){ cerr << "Error Mesh::readSTL: malformed stl file" << endl; return false; }
            // IEEE floating point 4-byte binary numerical representation, IEEE754, little endian
            tri.n = cgp::Vector((* ((float *) &inbuffer[inpos])), (* ((float *) &inbuffer[inpos+4])), (* ((float *) &inbuffer[inpos+8])));
            inpos += 12;

            // vertices
            for(i = 0; i < 3; i++)
            {
                if(inpos+12 >= insize){ cerr << "Error Mesh::readSTL: malformed stl file" << endl; return false; }
                vpos = cgp::Point((* ((float *) &inbuffer[inpos])), (* ((float *) &inbuffer[inpos+4])), (* ((float *) &inbuffer[inpos+8])));
                tri.v[i] = (int) verts.size();
                verts.push_back(vpos);
                inpos += 12;
            }
            tris.push_back(tri);
            t++;
            inpos += 2; // handle attribute byte count - which can simply be discarded
        }

        // tidy up
        delete inbuffer;
        infile.close();

        cerr << "num vertices = " << (int) verts.size() << endl;
        cerr << "num triangles = " << (int) tris.size() << endl;

        // STL provides a triangle soup so merge vertices that are coincident
        mergeVerts();
        // normal vectors at vertices are needed for rendering so derive from incident faces
        deriveVertNorms();
        if(basicValidity()){
            cerr << "========== basic validity: Pass ==========\n" << endl;

            if(manifoldValidity()){
                cerr << "========== manifold validity: Pass ==========\n" << endl;
            }else{
                cerr << "========== manifold validity: Fail ==========\n" << endl;
            }
        }else{
            cerr << "========== basic validity: Fail ==========\n" << endl;
        }
    }
    else
    {
        cerr << "Error Mesh::readSTL: unable to open " << filename << endl;
        return false;
    }
    return true;
}

bool Mesh::writeSTL(string filename)
{
    ofstream outfile;
    int t, p, numt;
    unsigned short val;

    outfile.open((char *) filename.c_str(), ios_base::out | ios_base::binary);
    if(outfile.is_open())
    {
        outfile.write("File Generated by Tesselator. Binary STL", 80); // skippable header
        numt = (int) tris.size();
        outfile.write((char *) &numt, 4); // number of triangles

        for(t = 0; t < numt; t++)
        {
            // normal
            outfile.write((char *) &tris[t].n.i, 4);
            outfile.write((char *) &tris[t].n.j, 4);
            outfile.write((char *) &tris[t].n.k, 4);

            // triangle vertices
            for(p = 0; p < 3; p++)
            {
                outfile.write((char *) &verts[tris[t].v[p]].x, 4);
                outfile.write((char *) &verts[tris[t].v[p]].y, 4);
                outfile.write((char *) &verts[tris[t].v[p]].z, 4);
            }

            // attribute byte count - null
            val = 0;
            outfile.write((char *) &val, 2);
        }

        // tidy up
        outfile.close();
    }
    else
    {
        cerr << "Error Mesh::writeSTL: unable to open " << filename << endl;
        return false;
    }
    return true;
}

void Mesh::insertEdge(Edge& edge){
    // Get the vertices in this edge in ascending order
    std::vector<int> edge_verts = {std::begin(edge.v), std::end(edge.v)};
    std::sort(edge_verts.begin(), edge_verts.end());

    // Generate the key pair for this edge
    std::pair<int, int> index(edge_verts[0], edge_verts[1]);

    // Push the edge into the correct bucket in the edge map
    edges[index].push_back(edge);
}

bool Mesh::basicValidity()
{
    edges.clear();
    cerr << "\nbasicValidity:________________________" << endl;

    // Loop through all the triangles
    for(int i = 0; i < tris.size(); i++){
        // Check if all the vertex indices in this triangle are valid vertex indices
        std::vector<int> tri_verts = {std::begin(tris[i].v), std::end(tris[i].v)};
        for(int v = 0; v < tri_verts.size(); v++){
            if(tri_verts[v] >= verts.size() || tri_verts[v] < 0){
                cerr << "All vertex indices within bounds: False" << endl;
                return false;
            }
        }

        // Insert all the edges from this triangle into the edge list
        Edge e1 = {{tri_verts[0], tri_verts[1]}, tris[i]};
        Edge e2 = {{tri_verts[1], tri_verts[2]}, tris[i]};
        Edge e3 = {{tri_verts[2], tri_verts[0]}, tris[i]};

        insertEdge(e1);
        insertEdge(e2);
        insertEdge(e3);
    }
    cerr << "All vertex indices within bounds: True" << endl;

    // Compute Euler's characteristic for this mesh and report it to the user
    uint V = verts.size();
    uint E = edges.size();
    uint F = tris.size();
    int eulers_char = V - E + F;
    cerr << "V = " << V << endl;
    cerr << "E = " << E << endl;
    cerr << "F = " << F << endl;
    cerr << "Eulers characteristic: V - E + F = " << eulers_char << endl;
    cerr << "Should be eq to 2 - 2G ____________T" << endl;
    cerr << "Expected genus of loaded mesh: " << (1 - (eulers_char/2.0)) << endl;

    // Check for dangling vertices
    std::vector<bool> verts_used(verts.size(), false);
    for(uint i = 0; i < tris.size(); i++){
        std::vector<int> tri_verts(std::begin(tris[i].v), std::end(tris[i].v));

        // Mark all the vertices in this triangle as used (not dangling)
        for(uint j = 0; j < tri_verts.size(); j++){
            verts_used[tri_verts[j]] = true;
        }
    }
    // Count all unused (dangling) vertices in the mesh and report them to the user
    int dangling_verts = std::count(verts_used.begin(), verts_used.end(), false);
    if(dangling_verts > 0){
        cerr << "No dangling vertices: False" << endl;
        return false;
    }
    cerr << "No dangling vertices: True" << endl;

    // If we reach this point without returning false, the mesh passes basicValidity
    return true;
}


bool Mesh::areOppositeEdges(Edge& edge1, Edge& edge2){
    return (edge1.v[0] == edge2.v[1] && edge1.v[1] == edge2.v[0]);
}

std::pair<int, int> Mesh::getNextEdge(std::pair<int, int> curr, std::vector<std::pair<int, int>>& adjacent_edges, std::vector<Triangle>& visited){
    // Get the triangle we are about to walk across and mark it as visited
    Triangle tri = edges[curr][0].adjacent_tri;
    if(std::find(visited.begin(), visited.end(), tri) != visited.end()){
        // Make sure we arent walking back where we came from
        tri = edges[curr][1].adjacent_tri;
    }
    visited.push_back(tri);

    // Generate the edges in this triangle. The next edge is whichever of these is
    // also in the adjacent_edges list.
    std::vector<std::pair<int, int>> tri_edges;
    tri_edges.push_back(std::make_pair(min(tri.v[0], tri.v[1]), max(tri.v[0], tri.v[1])));
    tri_edges.push_back(std::make_pair(min(tri.v[1], tri.v[2]), max(tri.v[1], tri.v[2])));
    tri_edges.push_back(std::make_pair(min(tri.v[2], tri.v[0]), max(tri.v[2], tri.v[0])));

    // Get rid of the current edge because we've already visited it.
    tri_edges.erase(std::remove(tri_edges.begin(), tri_edges.end(), curr), tri_edges.end());

    // Find and return whichever of the current triangle's edges we can find in the
    // adjacent_edges list. This is the next edge in our march around the vertex.
    if(std::find(adjacent_edges.begin(), adjacent_edges.end(), tri_edges[0]) != adjacent_edges.end()){
        return tri_edges[0];
    }else if(std::find(adjacent_edges.begin(), adjacent_edges.end(), tri_edges[1]) != adjacent_edges.end()){
        return tri_edges[1];
    }

    // We couldn't find the next edge, so return an invalid edge index.
    return std::make_pair(-1, -1);
}

bool Mesh::manifoldValidity()
{
    cerr << "\nmanifoldValidity:________________________" << endl;

    // To keep track of which edges are adjacent to each vertex
    std::vector<std::vector<std::pair<int, int>>> edges_at_vertex(verts.size());

    // Loop over the edges
    for(auto it = edges.begin(); it != edges.end(); it++){
        std::pair<int, int> index = it->first;
        std::vector<Edge> evec = it->second;

        // Check that every edge has 2 triangles associated with it (Closed mesh)
        if(evec.size() != 2){
            cerr << "Every edge has only 2 adjacent triangles: False" << endl;
            return false;
        }

        // Check orientation/consistent winding. Every pair of edges runs opposite to each other.
        if(!areOppositeEdges(evec[0], evec[1])){
            cerr << "Consistent orientation: False" << endl;
            return false;
        }

        // Put this edge into the edges_at_vertex vector
        edges_at_vertex[index.first].push_back(index);
        edges_at_vertex[index.second].push_back(index);
    }
    cerr << "Every edge has only 2 adjacent triangles: True" << endl;
    cerr << "Consistent orientation: True" << endl;

    // For each vertex
    for(int i = 0; i < verts.size(); i++){
        // Get the edges adjacent to that vertex
        std::vector<std::pair<int, int>> adjacent_edges = edges_at_vertex[i];

        std::vector<Triangle> visited; // Keep track of the triangles we visited

        // Try to march around the vertex from edge to edge
        std::pair<int, int> start = adjacent_edges[0];
        std::pair<int, int> curr = getNextEdge(start, adjacent_edges, visited);
        while(curr != start && !adjacent_edges.empty()){
            if(curr.first == -1 && curr.second == -1){
                // We couldn't find a valid next edge to march to
                cerr << "Every vertex has closed ring of triangles around it: False" << endl;
                return false;
            }

            adjacent_edges.erase(std::remove(adjacent_edges.begin(), adjacent_edges.end(), curr), adjacent_edges.end());
            curr = getNextEdge(curr, adjacent_edges, visited);
        }

        // If we didn't end up where we started, or we didn't visit all the triangles around that vertex,
        // then there was no closed ring of triangles around that vertex. (non 2-manifold)
        if(curr != start || visited.size() != edges_at_vertex[i].size()){
            cerr << "Every vertex has closed ring of triangles around it: False" << endl;
            return false;
        }
    }
    cerr << "Every vertex has closed ring of triangles around it: True" << endl;

    // If we reach this point without returning false, then the mesh passes manifoldValidity
    return true;
}
