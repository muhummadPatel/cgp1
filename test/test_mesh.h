#ifndef TILER_TEST_MESH_H
#define TILER_TEST_MESH_H


#include <string>
#include <cppunit/extensions/HelperMacros.h>

#define private public
#include "tesselate/mesh.h"
#undef private

/// Test code for @ref Mesh
class TestMesh : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestMesh);
    CPPUNIT_TEST(testMeshing);
    CPPUNIT_TEST(testLoadsValidMesh);
    CPPUNIT_TEST(testFailsOpenMesh);
    CPPUNIT_TEST(testFailsBadOrientationMesh);
    CPPUNIT_TEST(testFailsNon2ManifoldMesh);
    CPPUNIT_TEST(testFailsDanglingVertexMesh);
    CPPUNIT_TEST(testFailsInvalidVertexMesh);
    CPPUNIT_TEST(testGenus0Mesh);
    CPPUNIT_TEST(testGenus1Mesh);
    CPPUNIT_TEST(testGenus2Mesh);
    CPPUNIT_TEST_SUITE_END();

private:
    Mesh * mesh;

public:

    /// Initialization before unit tests
    void setUp();

    /// Tidying up after unit tests
    void tearDown();

    /**
     * Run standard validity tests on bunny mesh
     * @pre bunny.stl must be located in the <project_root>/test/meshes directory
     */
    void testMeshing();

    /**
     * Test that a perfectly valid mesh can be loaded without errors and passes both
     * validity checks.
     * @pre cube_valid.stl must be located in the <project_root>/test/meshes directory
     */
    void testLoadsValidMesh();

    /**
     * Test that a mesh with a hole in it fails as it it supposed to.
     * @pre cube_open.stl must be located in the <project_root>/test/meshes directory
     */
    void testFailsOpenMesh();

    /**
     * Test that a mesh with inconsistent winding/flipped normals fails as expected.
     * @pre cube_inconsistent_orientation.stl must be located in the <project_root>/test/meshes directory
     */
    void testFailsBadOrientationMesh();

    /**
     * Test that a non 2-manifold mesh fails as expected.
     * @pre cube_point_non_2_manifold.stl must be located in the <project_root>/test/meshes directory
     */
    void testFailsNon2ManifoldMesh();

    /**
     * Test that a mesh with a dangling vertex fails as expected.
     * The dangling vertex will be added to the data structures by the test itself
     * and we then check that the validity checks fail as expected and report the
     * dangling vertex.
     * @pre cube_valid.stl must be located in the <project_root>/test/meshes directory
     */
    void testFailsDanglingVertexMesh();

    /**
     * Test that a mesh with a reference to an invalid vertex index fails the
     * basic validity check. We load a valid mesh and then edit the data structure
     * to manufacture this situation and check that the test reports the fault correctly.
     * @pre cube_valid.stl must be located in the <project_root>/test/meshes directory
     */
    void testFailsInvalidVertexMesh();

    /**
     * Test that the genus of a mesh with genus 0 is correctly reorted to the user.
     * @pre cube_valid.stl must be located in the <project_root>/test/meshes directory
     */
    void testGenus0Mesh();

    /**
     * Test that the genus of a mesh with genus 1 is correctly reported to the user.
     * @pre torus_genus1.stl must be located in the <project_root>/test/meshes directory
     */
    void testGenus1Mesh();

    /**
     * Test that the genus of a mesh with genus 2 is correctly reported to the user.
     * @pre double_torus_genus2.stl must be located in the <project_root>/test/meshes directory
     */
    void testGenus2Mesh();
};

#endif /* !TILER_TEST_MESH_H */
