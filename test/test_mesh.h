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

    void testLoadsValidMesh();

    void testFailsOpenMesh();

    void testFailsBadOrientationMesh();

    void testFailsNon2ManifoldMesh();

    void testFailsDanglingVertexMesh();

    void testFailsInvalidVertexMesh();

    void testGenus0Mesh();

    void testGenus1Mesh();

    void testGenus2Mesh();
};

#endif /* !TILER_TEST_MESH_H */
