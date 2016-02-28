#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <test/testutil.h>
#include "test_mesh.h"
#include <stdio.h>
#include <cstdint>
#include <regex>
#include <sstream>
#include <iostream>
#include <string>
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

void TestMesh::testFailsDanglingVertexMesh(){
    // Load a known valid mesh
    mesh->readSTL("../../test/meshes/cube_valid.stl");

    // Add a dangling vertex to the verts vector
    cgp::Point dangling_vert = cgp::Point();
    mesh->verts.push_back(dangling_vert);

    // Should now fail basic validity because of the dangling vert
    CPPUNIT_ASSERT(!mesh->basicValidity());
}

void TestMesh::testFailsInvalidVertexMesh(){
    // Load a known valid mesh
    mesh->readSTL("../../test/meshes/cube_valid.stl");

    // Change one of the vertices in a triangle to be invalid(> verts.size())
    mesh->tris[0].v[0] = mesh->verts.size() + 5;
    //Should now fail basci validity because of the invalid vertex index
    CPPUNIT_ASSERT(!mesh->basicValidity());

    // Same as above, but with -ve invalid vertex index
    mesh->tris[0].v[0] = -5;
    CPPUNIT_ASSERT(!mesh->basicValidity());
}

void TestMesh::testGenus0Mesh(){
    mesh->readSTL("../../test/meshes/cube_valid.stl");

    // Redirect cerr to stringstream so we can inspect it
    std::ostringstream oss;
    std::streambuf* orig_buf(cerr.rdbuf(oss.rdbuf()));
    mesh->basicValidity(); // Run the validity check so we can capture the output
    cerr.rdbuf(orig_buf); // Reset cerr to print to the screen

    const string output = oss.str();
    cerr << output; // Show the user the captured output

    // Regex magic to extract the genus reported to the user (I LOOOOVEE REGEX XD)
    std::regex regex("Expected genus of loaded mesh: (\\d+)\\s");
    std::smatch match;
    if (std::regex_search(output.begin(), output.end(), match, regex)){
        int genus = std::stoi(match.str(1));
        CPPUNIT_ASSERT(genus == 0);
    }else{
        CPPUNIT_FAIL("Cannot find genus of the mesh in the output.");
    }
}

void TestMesh::testGenus1Mesh(){
    mesh->readSTL("../../test/meshes/torus_genus1.stl");

    // Redirect cerr to stringstream so we can inspect it
    std::ostringstream oss;
    std::streambuf* orig_buf(cerr.rdbuf(oss.rdbuf()));
    mesh->basicValidity(); // Run the validity check so we can capture the output
    cerr.rdbuf(orig_buf); // Reset cerr to print to the screen

    const string output = oss.str();
    cerr << output; // Show the user the captured output

    // Extract the reported genus and check that it matches the expected value
    std::regex regex("Expected genus of loaded mesh: (\\d+)\\s");
    std::smatch match;
    if (std::regex_search(output.begin(), output.end(), match, regex)){
        int genus = std::stoi(match.str(1));
        CPPUNIT_ASSERT(genus == 1);
    }else{
        CPPUNIT_FAIL("Cannot find genus of the mesh in the output.");
    }
}

void TestMesh::testGenus2Mesh(){
    mesh->readSTL("../../test/meshes/double_torus_genus2.stl");

    // Redirect cerr to stringstream so we can inspect it
    std::ostringstream oss;
    std::streambuf* orig_buf(cerr.rdbuf(oss.rdbuf()));
    mesh->basicValidity(); // Run the validity check so we can capture the output
    cerr.rdbuf(orig_buf); // Reset cerr to print to the screen

    const string output = oss.str();
    cerr << output; // Show the user the captured output

    // Extract the reported genus and check that it matches the expected value
    std::regex regex("Expected genus of loaded mesh: (\\d+)\\s");
    std::smatch match;
    if (std::regex_search(output.begin(), output.end(), match, regex)){
        int genus = std::stoi(match.str(1));
        CPPUNIT_ASSERT(genus == 2);
    }else{
        CPPUNIT_FAIL("Cannot find genus of the mesh in the output.");
    }
}

//#if 0 /* Disabled since it crashes the whole test suite */
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TestMesh, TestSet::perCommit());
//#endif
