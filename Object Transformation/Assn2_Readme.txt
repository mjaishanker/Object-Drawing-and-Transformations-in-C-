Mani Jaishanker
Computer Graphics Assignment 2

#Description
The program will draw a model after reading a object file and parsing it.
It will utilize the midpoint algorithm to draw lines between vertices
and use matrix multiplication to transform the model by translating,
rotating and scaling them. The view can be altered from orthogonal view
to perspective view.

#Run the code in Linux Machine using 
* g++ assn1.cpp -lglut -lGL
* with ./a.out use arg coomand to give location of the .obj file.

#OpenGl
Utilizes the <GL/glut.h>, <GL/gl.h> and <GL/glu.h> libraries to run the
openGL functions

#Users use
* The program will display the initial drawing of an object from the 
parsed file. The view will be in orthogonal at first. 

* The transform the model, the user will press the following keys:
- t - translate
- e - scale
- r - rotate
- v - change view
