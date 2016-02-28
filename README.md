##Assignment 1: Closed Orientable 2-Manifold Validity of Meshes

**Author** : Muhummad Yunus Patel  
**Student#** : PTLMUH006  
**Date** : 28-February-2016

**NOTE:** Should you require any clarification, please feel free to contact me
 at muhummad.patel@gmail.com. Thank you. :)

###Description:
This practical checks the validity of a loaded triangle mesh. It checks that the mesh is a closed, orientable, 2-manifold mesh suitable for 3d printing. It reports on the validity of the mesh (by printing results in the terminal/console) when the mesh is loaded.

The solution was implemented mainly in two methods in the mesh.cpp file, viz. basicValidity(), and manifoldValidity(). These methods check the validity of the mesh as requested and report on the various validity aspects. These methods have also been thoroughly tested. Unit tests have been added to check each aspect of these two validity checking methods. The unit tests make use of a number of test meshes located in test/meshes.

###Compiling and Running the Solution:
* See instructions in InstallReadme.txt for instructions on setting up and compiling the project.
* To run the the Tesselate program:
    * First set up and build the project (see previous point).
    * Navigate to build/tesselate
    * Run ./tessviewer
    * Load a mesh by using the options in the menu bar. (Sample meshes provided in /meshes)
    * Note that the result of the validity checks on the loaded mesh will be reported in the console/terminal.
* To run the unit tests:
    * First set up and build the project (see first point).
    * Navigate to build/test
    * Run ./tilertest
    * All the tests should pass and you should see 'OK (10)'.
