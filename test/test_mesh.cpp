#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <test/testutil.h>
#include "test_mesh.h"
#include <stdio.h>
#include <cstdint>
#include <sstream>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>

void TestMesh::setUp()
{
    mesh = new Mesh();
}

void TestMesh::tearDown()
{
    delete mesh;
}

void TestMesh::testMeshing()
{
    mesh->readSTL("../../test/meshes/bunny.stl");
    CPPUNIT_ASSERT(mesh->basicValidity());
    CPPUNIT_ASSERT(!mesh->manifoldValidity()); // bunny has known holes in the bottom
}

void TestMesh::testLoadsValidMesh(){
    mesh->readSTL("../../test/meshes/cube_valid.stl");
    CPPUNIT_ASSERT(mesh->basicValidity());
    CPPUNIT_ASSERT(mesh->manifoldValidity());
}

void TestMesh::testFailsOpenMesh(){
    mesh->readSTL("../../test/meshes/cube_open.stl");
    CPPUNIT_ASSERT(mesh->basicValidity());
    CPPUNIT_ASSERT(!mesh->manifoldValidity());
}

void TestMesh::testFailsBadOrientationMesh(){
    mesh->readSTL("../../test/meshes/cube_inconsistent_orientation.stl");
    CPPUNIT_ASSERT(mesh->basicValidity());
    CPPUNIT_ASSERT(!mesh->manifoldValidity());
}

void TestMesh::testFailsNon2ManifoldMesh(){
    mesh->readSTL("../../test/meshes/cube_point_non_2_manifold.stl");
    CPPUNIT_ASSERT(mesh->basicValidity());
    CPPUNIT_ASSERT(!mesh->manifoldValidity());
}


//#if 0 /* Disabled since it crashes the whole test suite */
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TestMesh, TestSet::perCommit());
//#endif
