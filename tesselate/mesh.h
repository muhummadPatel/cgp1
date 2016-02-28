/**
 * @file
 *
 * Data structure representing a triangle mesh in 3D space.
 */

#ifndef _MESH
#define _MESH

#include <vector>
#include <stdio.h>
#include <iostream>
#include "renderer.h"

#include <utility> //std::pair
#include <boost/functional/hash.hpp>

using namespace std;

const int sphperdim = 20;

/**
 * A triangle in 3D space, with 3 indices into a vertex list and an outward facing normal. Triangle winding is counterclockwise.
 */
struct Triangle
{
    int v[3];   ///< index into the vertex list for triangle vertices
    cgp::Vector n;   ///< outward facing unit normal to the triangle
};

// To allow for checking equality of triangles.
inline bool operator==(const Triangle& lhs, const Triangle& rhs)
{
    bool vertsEqual = lhs.v[0] == rhs.v[0] && lhs.v[1] == rhs.v[1] &&lhs.v[2] == rhs.v[2];
    bool normalEqual = true;
    return vertsEqual && normalEqual;
}

/**
 * An edge in 3D space, with 2 indices into a vertex list. The order of the vertices has no special significance.
 */
struct Edge
{
    int v[2];   ///< indices into the vertex list for edge endpoints
    Triangle adjacent_tri;
};

/**
 * A sphere in 3D space, consisting of a center and radius. Used for bounding sphere hierarchy acceleration.
 */
class Sphere
{
public:
    cgp::Point c;  ///< sphere center
    float r;    ///< sphere radius
    std::vector<int> ind; ///< triangle indices included in the bounding sphere. Used for acceleration struct, otherwise ignored

    /**
     * Test whether a point falls inside the sphere
     * @param pnt   point to test for containment
     * @retval true if the point falls within the sphere,
     * @retval false otherwise
     */
    bool pointInSphere(cgp::Point pnt);
};

/**
 * A triangle mesh in 3D space. Ideally this should represent a closed 2-manifold but there are validity tests to ensure this.
 */
class Mesh
{
private:
    std::vector<cgp::Point> verts; ///< vertices of the tesselation structure
    std::vector<cgp::Vector> norms;  ///< per vertex normals
    std::vector<Triangle> tris; ///< edges connecting vertices
    GLfloat * col;              ///< (r,g,b,a) colour
    float scale;                ///< scaling factor
    cgp::Vector trx;                 ///< translation
    float xrot, yrot, zrot;     ///< rotation angles about x, y, and z axes
    std::vector<Sphere> boundspheres; ///< bounding sphere accel structure

    std::unordered_map<std::pair<int, int>, std::vector<Edge>, boost::hash<std::pair<int, int>>> edges;

    /**
     * Search list of vertices to find matching point
     * @param pnt       point to search for in vertex list
     * @param[out] idx  index of point in list if found, otherwise -1
     * @retval true  if the point is found in the vertex list,
     * @retval false otherwise
     */
    bool findVert(cgp::Point pnt, int &idx);

    /**
     * Search a list of edges to find matching edge
     * @param edges     edge list
     * @param e         edge to find
     * @param[out] idx  index of edge in list if found, otherwise -1
     * @retval true  if the point is found in the vertex list,
     * @retval false otherwise
     */
    bool findEdge(vector<Edge> edges, Edge e, int &idx);

    /**
     * Construct a hash key based on a 3D point
     * @param pnt   point to convert to key
     * @param bbox  bounding box enclosing all mesh vertices
     * @retval hash key
     */
    long hashVert(cgp::Point pnt, cgp::BoundBox bbox);

    /// Connect triangles together by merging duplicate vertices
    void mergeVerts();

    /// Generate vertex normals by averaging normals of the surrounding faces
    void deriveVertNorms();

    /// Generate face normals from triangle vertex positions
    void deriveFaceNorms();

    /**
     * Composite rotations, translation and scaling into a single transformation matrix
     * @param tfm   composited transformation matrix
     */
    void buildTransform(glm::mat4x4 &tfm);

    /**
     * Inserts an edge into the edge list. The edge list maps from a pair<int, int>
     * to a vector<Edge>. The pair is made up of the two vertices that make up that
     * edge sorted in ascending order. The vector stores all the instances of that edge
     * in the mesh (e.g. there would be two entries in the vector if the edge is
     * adjacent to 2 triangles).
     * @param edge  The edge to be inserted into the edge list.
     */
    void insertEdge(Edge& edge);

    /**
     * Checks if the given 2 edges are opposites/the same edge. Two edges are opposite
     * if they consist of the same vertices but in reverse order. For example, this would
     * be the case when two triangles share the same edge.
     * @param edge1 The first edge to check
     * @param edge2 The second edge to check
     * @retval true if edge1 and edge2 are the same edge in opposite directions
     * @retval false otherwise
     */
    bool areOppositeEdges(Edge& edge1, Edge& edge2);

    /**
     * Gets the next edge in the adjacent_edges vector. This is the method used to
     * march along the faces around a vertex when checking if the mesh obeys the 2-manifold
     * property. It looks at the current edge and finds the next edge around the vertex
     * to allow us to move around the vertex. It also updates the visited vector to keep track
     * of which faces we have already visited.
     * @param curr the current edge we are looking at
     * @param adjacent_edges the list of edges around the vertex. The next edge is somewhere in this list.
     * @param visited This is the list of triangles/faces we have already visited. We need to update this and make sure we aren't moving backwards.
     * @retval The next edge on our march around the vertex if we can find one.
     * @retval pair(-1, -1) if we cannot find a valid edge to move to next.
     */
    std::pair<int, int> getNextEdge(std::pair<int, int> curr, std::vector<std::pair<int, int>>& adjacent_edges, std::vector<Triangle>& visited);

public:

    ShapeGeometry geom;         ///< renderable version of mesh

    Mesh();

    ~Mesh();

    /// Remove all vertices and triangles, resetting the structure
    void clear();

    /// Test whether mesh is empty of any geometry (true if empty, false otherwise)
    bool empty(){ return verts.empty(); }

    /// Setter for scale
    void setScale(float scf){ scale = scf; }

    /// Getter for scale
    float getScale(){ return scale; }

    /// Setter for translation
    void setTranslation(cgp::Vector tvec){ trx = tvec; }

    /// Getter for translation
    cgp::Vector getTranslation(){ return trx; }

    /// Setter for rotation angles
    void setRotations(float ax, float ay, float az){ xrot = ax; yrot = ay; zrot = az; }

    /// Getter for rotation angles
    void getRotations(float &ax, float &ay, float &az){ ax = xrot; ay = yrot; az = zrot; }

    /// Setter for colour
    void setColour(GLfloat * setcol){ col = setcol; }

    /**
     * Generate triangle mesh geometry for OpenGL rendering
     * @param view      current view parameters
     * @param[out] sdd  openGL parameters required to draw this geometry
     * @retval true  if buffers are bound successfully, in which case sdd is valid,
     * @retval false otherwise
     */
    bool genGeometry(View * view, ShapeDrawData &sdd);

    /**
     * Scale geometry to fit bounding cube centered at origin
     * @param sidelen   length of one side of the bounding cube
     */
    void boxFit(float sidelen);

    /**
     * Read in triangle mesh from STL format binary file
     * @param filename  name of file to load (STL format)
     * @retval true  if load succeeds,
     * @retval false otherwise.
     */
    bool readSTL(string filename);

    /**
     * Write triangle mesh to STL format binary file
     * @param filename  name of file to save (STL format)
     * @retval true  if save succeeds,
     * @retval false otherwise.
     */
    bool writeSTL(string filename);

    /**
     * Basic mesh validity tests - report euler's characteristic, no dangling vertices, edge indices within bounds of the vertex list
     * @retval true if basic validity tests are passed,
     * @retval false otherwise
     * @todo basicValidity requires completing for CGP Prac1
     */
    bool basicValidity();

    /**
     * Check that the mesh is a closed two-manifold - every edge has two incident triangles, every vertex has
     *                                                a closed ring of triangles around it
     * This test does not include self-intersection of individual triangles as this is outside the scope.
     * @retval true if the mesh is two-manifold,
     * @retval false otherwise
     * @todo manifoldValidity requires completing for CGP Prac1
     */
    bool manifoldValidity();
};

#endif
