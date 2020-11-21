
#include <GL/gl.h>
#include <GL/glu.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <limits>
#include <cassert>
#include <vector>
#include <cstring>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <GL/glut.h>
using namespace std;

#include <GL/glut.h>   // The GL Utility Toolkit (Glut) Header

#define WIDTH 500
#define HEIGHT 500

// Scale the fit the model in window
float scale = sqrt(WIDTH * WIDTH + HEIGHT * HEIGHT);

int x_last, y_last, D = -400;
int translateSwitch = 0, scaleSwitch = 0, rotSwitch = 0;
float centerX, centerY, centerZ, centerW, maxValue = 0.0, minValue = 0.0;
bool viewSwitch = false;
// Perspective view matrix
float presMat[16]
= {
	1.0, 0, 0, 0,
	0, 1.0, 0, 0,
	0, 0, 1.0, 0,
	0, 0, 1.0 / D, 0
};

// Vertices structure to hold main vertices of the model
struct verticesStruct {
	// w is 1 if orthogonal view and z/d if perspective vew
	float x, y, z, w;
	verticesStruct() {};
	verticesStruct(float X, float Y, float Z) { x = X; y = Y; z = Z; };
};

// Face struct to hold the indices of the faces vertices
struct faceStruct {
	int v1, v2, v3;
};

// Structure for holding texture vertices
struct vertTextStruct {
	float x, y;
};

//Structure to hold the normal vertices
struct normalStruct {
	float x, y, z;
};

// Vectors to hold all the vertices and indices
std::vector< verticesStruct > vertices;
std::vector< verticesStruct > MainVertices;
std::vector< vertTextStruct > uvs;
std::vector< normalStruct > normals;
std::vector< faceStruct > faces;

// function to center the model when trasnformation is done
void centerVertices() {
	float xSum = 0.0, ySum = 0.0, zSum = 0.0, wSum = 0.0;
	int vertexSize = MainVertices.size();

	for (size_t i = 0; i < MainVertices.size(); ++i) {
		xSum += MainVertices[i].x;
		ySum += MainVertices[i].y;
		zSum += MainVertices[i].z;
		wSum += MainVertices[i].w;
	}
	centerX = xSum / vertexSize;
	centerY = ySum / vertexSize;
	centerZ = zSum / vertexSize;
	centerW = wSum / vertexSize;

	for (size_t i = 0; i < MainVertices.size(); i++) {
		MainVertices[i].x = MainVertices[i].x - centerX;
		MainVertices[i].y = MainVertices[i].y - centerY;
		MainVertices[i].z = MainVertices[i].z - centerZ - maxValue - minValue;
		MainVertices[i].w = MainVertices[i].w - centerW;
	}
}

// Function to center and normalize the vertices at the begining of the display
void normalizeVertices() {
	float xSum = 0.0, ySum = 0.0, zSum = 0.0, wSum = 0.0;
	float centerX, centerY, centerZ, centerW;
	int vertexSize = MainVertices.size();

	vector<verticesStruct> returnVert;

	for (size_t i = 0; i < MainVertices.size(); ++i) {
		xSum += MainVertices[i].x;
		ySum += MainVertices[i].y;
		zSum += MainVertices[i].z;
		wSum += MainVertices[i].w;
	}
	centerX = xSum / vertexSize;
	centerY = ySum / vertexSize;
	centerZ = zSum / vertexSize;
	centerW = wSum / vertexSize;

	for (size_t i = 0; i < MainVertices.size(); i++) {
		maxValue = max(maxValue, abs(MainVertices[i].x));
		maxValue = max(maxValue, abs(MainVertices[i].y));
		maxValue = max(maxValue, abs(MainVertices[i].z));
		minValue = min(maxValue, abs(MainVertices[i].x));
		minValue = min(maxValue, abs(MainVertices[i].y));
		minValue = min(maxValue, abs(MainVertices[i].z));
	}

	for (size_t i = 0; i < MainVertices.size(); i++) {
		MainVertices[i].x = MainVertices[i].x - centerX;
		MainVertices[i].y = MainVertices[i].y - centerY;
		MainVertices[i].z = MainVertices[i].z - centerZ;
		MainVertices[i].w = MainVertices[i].w - centerW;
	}
	// Scaling the model
	scale = scale / (maxValue);

	for (size_t i = 0; i < MainVertices.size(); i++) {
		MainVertices[i].x = MainVertices[i].x * scale;
		MainVertices[i].y = MainVertices[i].y * scale;
		MainVertices[i].z = MainVertices[i].z * scale;
	}
}

// Parser to read the object file to take vertices.
void readObjFile(const char* path, std::vector < verticesStruct >& main_vertices,
	std::vector < vertTextStruct >& main_uvs, std::vector < normalStruct >& main_normals) {

	std::vector<unsigned int> vertexIndices, uvIndices, normalIndices;
	std::vector<unsigned int> smoothGroup;
	std::vector<verticesStruct> temp_vertices;
	std::vector<vertTextStruct> temp_textures;
	std::vector<normalStruct> temp_normals;

	fstream file;
	file.open(path);
	string line;
	if (!file) {
		cout << "Unable to open file";
		exit(1); // terminate with error
	}
	fstream fin(path, fstream::in);
	while (std::getline(file, line)) {
		std::istringstream lineSS(line);
		std::string lineType;
		lineSS >> lineType;

		if (lineType == "v") {
			verticesStruct vertex;
			lineSS >> vertex.x >> vertex.y >> vertex.z;
			temp_vertices.push_back(vertex);
		}

		// texture
		if (lineType == "vt")
		{
			vertTextStruct uv;
			lineSS >> uv.x >> uv.y;
			temp_textures.push_back(uv);
		}

		// normal
		if (lineType == "vn")
		{
			normalStruct normal;
			lineSS >> normal.x >> normal.y >> normal.z;
			temp_normals.push_back(normal);
		}

		// face
		if (lineType == "f")
		{
			std::string tempStr;
			while (lineSS >> tempStr)
			{
				std::istringstream ref(tempStr);
				std::string vLine, vtLine, vnLine;
				std::getline(ref, vLine, '/');
				std::getline(ref, vtLine, '/');
				std::getline(ref, vnLine, '/');
				int v = atoi(vLine.c_str());
				int vt = atoi(vtLine.c_str());
				int vn = atoi(vnLine.c_str());
				vertexIndices.push_back(v);
			}
		}
	}
	MainVertices = temp_vertices;
	faceStruct tempFace;

	// Store the vertices indices in the vector of faces structure
	for (int i = 0; i < vertexIndices.size() - 2; i = i + 3) {
		tempFace.v1 = vertexIndices[i];
		tempFace.v2 = vertexIndices[i + 1];
		tempFace.v3 = vertexIndices[i + 2];
		faces.push_back(tempFace);
	}
}

/***************************************************************************/

void init_window()
/* Clear the image area, and set up the coordinate system */
{

	/* Clear the window */
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);
	glOrtho(-600, 600, -600, 600, -600, 600);
}

/***************************************************************************/

void write_pixel(int x, int y, double r, double g, double b)
/* Turn on the pixel found at x,y */
{

	glColor3f(r, g, b);
	glBegin(GL_POINTS);
	glVertex3i(x, y, 0);
	glEnd();
}

//***************************************************************************/
// Midpoint algorithm to draw lines
void midPoint(int x0, int y0, int x1, int y1) {
	float m = (y1 - y0);
	m = m / (x1 - x0);
	if ((y1 == y0) || (x1 == x0)) {
		m = 0;
	}

	if (m == 0) {
		if (x1 == x0 && y0 < y1) {
			int x = x0;
			int y = y0;
			for (y = y0; y < y1; y++) {
				write_pixel(x, y, 1.0, 1.0, 1.0);
			}
		}

		else if (x1 == x0 && y1 < y0) {
			int x = x0;
			int y = y1;
			for (y = y1; y < y0; y++) {
				write_pixel(x, y, 1.0, 1.0, 1.0);
			}
		}

		else if (y1 == y0 && x0 < x1) {
			int x = x0;
			int y = y0;
			for (x = x0; x < x1; x++) {
				write_pixel(x, y, 1.0, 1.0, 1.0);
			}
		}

		else if (y1 == y0 && x1 < x0) {
			int x = x1;
			int y = y0;
			for (x = x1; x < x0; x++) {
				write_pixel(x, y, 1.0, 1.0, 1.0);
			}
		}
	}

	if (x0 < x1 && y0 < y1) {
		if (m <= 1 && m > 0) {
			int X = x0;
			int Y = y0;
			int a = y1 - y0;
			int b = -(x1 - x0);
			float d = a + b / 2;
			for (X = x0; X < x1; X++) {
				write_pixel(X, Y, 1.0, 1.0, 1.0);
				if (d < 0) {
					d = d + a;
				}
				else {
					//Y++;
					d = d + a + b;
					Y++;
				}
			}
		}

		else if (m >= 1) {
			int X = x0;
			int Y = y0;
			int a = y1 - y0;
			int b = x1 - x0;
			float d = b - (a / 2);
			for (Y = y0; Y < y1; Y++) {
				write_pixel(X, Y, 1.0, 1.0, 1.0);
				if (d < 0) {
					d = d + b;
				}
				else {
					X++;
					d = d + (b - a);
				}
			}
		}
	}

	if (x1 < x0 && y1 < y0) {
		if (m <= 1 && m > 0) {
			int X = x1;
			int Y = y1;
			int a = y1 - y0;
			int b = -(x1 - x0);
			float d = -a - b / 2;
			for (X = x1; X < x0; X++) {
				write_pixel(X, Y, 1.0, 1.0, 1.0);
				if (d < 0) {
					d = d - a;
				}
				else {
					//Y++;
					d = d - a - b;
					Y++;
				}
			}
		}

		else if (m >= 1) {
			int X = x0;
			int Y = y0;
			int a = y1 - y0;
			int b = -(x1 - x0);
			//printf("a, b is (%d,%d)\n", a, b);
			float d = -b - (a / 2);
			//printf("d is (%f)\n", d);
			for (Y = y0; Y >= y1; --Y) {
				write_pixel(X, Y, 1.0, 1.0, 1.0);
				if (d < 0) {
					d = d + b;
				}
				else {
					d = d + (b + a);
					--X;
				}
			}
		}
	}

	else if (y0 > y1&& x0 < x1) {
		if (m < 0 && m >= -1) {
			int X = x1; // mod
			int Y = y1; // mod
			int a = y1 - y0; //dy
			int b = -(x1 - x0); // -dx
			float d = -a + b / 2;
			for (X = x1 - 1; X >= x0; --X) {
				write_pixel(X, Y, 1.0, 1.0, 1.0);
				if (d < 0) {
					d = d - a;
				}
				else {
					//Y++;
					d = d - a + b;
					Y++;
				}
			}

		}

		else if (m <= -1) {
			int X = x0; // mod
			int Y = y0; // mod
			int a = y1 - y0;
			int b = -(x1 - x0);
			float d = -b + (a / 2);
			for (Y = y0; Y >= y1; Y--) {

				write_pixel(X, Y, 1.0, 1.0, 1.0);
				if (d < 0) {
					d = d - b;
				}
				else {
					X++;
					d = d - (b - a);
				}
			}
		}
	}

	else if (x0 > x1&& y0 < y1) {
		if (m < 0 && m >= -1) {
			int X = x0; // mod
			int Y = y0; // mod
			int a = y1 - y0; //dy
			int b = -(x1 - x0); // -dx
			float d = a - b / 2;
			for (X = x0; X >= x1; X--) {
				write_pixel(X, Y, 1.0, 1.0, 1.0);
				if (d < 0) {
					d = d + a;
				}
				else {
					//Y++;
					d = d + a - b;
					Y++;
				}
			}

		}

		else if (m <= -1) {
			int X = x1; // mod
			int Y = y1; // mod
			int a = y1 - y0;
			int b = x1 - x0;
			float d = b + (a / 2);
			for (Y = y1 - 1; Y >= y0; --Y) {
				write_pixel(X, Y, 1.0, 1.0, 1.0);
				if (d < 0) {
					d = d - b;
				}
				else {
					X++;
					d = d - (b + a);
				}
			}
		}
	}
}

// function to return to center points of the model to center it when using transformation matrix
verticesStruct getCenterVertices() {
	float xSum = 0.0, ySum = 0.0, zSum = 0.0, wSum = 0.0, maxValue = 0.0, minValue = 0.0;
	float centerX, centerY, centerZ, centerW;
	int vertexSize = MainVertices.size();

	verticesStruct CVertices;

	for (size_t i = 0; i < MainVertices.size(); ++i) {
		xSum += MainVertices[i].x;
		ySum += MainVertices[i].y;
		zSum += MainVertices[i].z;
		wSum += MainVertices[i].w;
	}
	centerX = xSum / vertexSize;
	centerY = ySum / vertexSize;
	centerZ = zSum / vertexSize;
	centerW = wSum / vertexSize;

	CVertices.x = centerX;
	CVertices.y = centerY;
	CVertices.z = centerY;

	return CVertices;
}

// Function to multiple 4 by 4 matrix with a 4 by 1 matrix
void matrixMult(float iMat[4][4], int indx) {
	float inputvert[4][1] = { {MainVertices[indx].x}, {MainVertices[indx].y}, {MainVertices[indx].z}, {1.0} };
	float outputvert[4][1] = { { 0.0}, {0.0}, {0.0}, {0.0} };
	int k = 0;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 1; ++j) {
			for (int k = 0; k < 4; ++k)
			{
				outputvert[i][j] += iMat[i][k] * inputvert[k][j];
			}
		}
	}
	MainVertices[indx].x = outputvert[0][0];
	MainVertices[indx].y = outputvert[1][0];
	MainVertices[indx].z = outputvert[2][0];
}

// function to determine the type of transformation and performing the,
void detectTransform(int transformType, char direc) {
	// Main matrices to do composite transformation multiplication
	float identityMat[16] =
	{
		1.0, 0, 0, 0,
		0, 1.0, 0, 0,
		0, 0, 1.0, 0,
		0, 0, 0, 1.0
	};
	float CompositeMatrix[4][4]
	= {
		{1.0, 0, 0, 0},
		{0, 1.0, 0, 0},
		{0, 0, 1.0, 0},
		{0, 0, 0, 1.0}
	};
	float InvCompositeMatrix[4][4]
	= {
		{1.0, 0, 0, 0},
		{0, 1.0, 0, 0 },
		{0, 0, 1.0, 0},
		{0, 0, 0, 1.0}
	};
	float GlobalCompositeMatrix[4][4]
	= {
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0}
	};
	float iMatrx[4][4] =
	{
		{1.0, 0, 0, 0},
		{0, 1.0, 0, 0},
		{0, 0, 1.0, 0},
		{0, 0, 0, 1.0}
	};

	float tempCompMat[4][4] =
	{
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0}
	};
	if (transformType == 1) { // transform is translation
		if (direc == 'a')
			iMatrx[0][3] = -10;
		else if (direc == 'd')
			iMatrx[0][3] = 10;
		else if (direc == 'w')
			iMatrx[1][3] = 10;
		else if (direc == 's')
			iMatrx[1][3] = -10;
		// Translate the model
		for (size_t i = 0; i < MainVertices.size(); i++)
			matrixMult(iMatrx, i);
	}
	else if (transformType == 2) { // transform is scaling
		verticesStruct tempVertices;
		if (direc == 'a')
			iMatrx[0][0] = 1.05, iMatrx[1][1] = 1.05, iMatrx[2][2] = 1.05;
		else if (direc == 'd')
			iMatrx[0][0] = 1.0 / 1.05, iMatrx[1][1] = 1.0 / 1.05, iMatrx[2][2] = 1.0 / 1.05;
		else if (direc == 'w')
			iMatrx[0][0] = 1.05, iMatrx[1][1] = 1.05, iMatrx[2][2] = 1.05;
		else if (direc == 's')
			iMatrx[0][0] = 1.0 / 1.05, iMatrx[1][1] = 1.0 / 1.05, iMatrx[2][2] = 1.0 / 1.05;
		// center vertices
		tempVertices = getCenterVertices();
		// composite multiplication of center, translation matrix, scaling, and inverse translation matrices
		CompositeMatrix[0][3] = tempVertices.x, CompositeMatrix[1][3] = tempVertices.y, CompositeMatrix[2][3] = tempVertices.z;
		InvCompositeMatrix[0][3] = -1.0 * (tempVertices.x), InvCompositeMatrix[1][3] = -1.0 * (tempVertices.y), InvCompositeMatrix[2][3] = -1.0 * (tempVertices.z);
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				for (int k = 0; k < 4; ++k)
				{
					tempCompMat[i][j] += iMatrx[i][k] * CompositeMatrix[k][j];
				}
			}
		}
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				for (int k = 0; k < 4; ++k)
				{
					GlobalCompositeMatrix[i][j] += InvCompositeMatrix[i][k] * tempCompMat[k][j];
				}
			}
		}
		// Multiply final composite matrix with the vertex points
		for (size_t i = 0; i < MainVertices.size(); i++)
			matrixMult(GlobalCompositeMatrix, i);
		// Finally center the model
		centerVertices();
	}

	else if (transformType == 3) { // transform is rotating
		verticesStruct tempVertices;
		if (direc == 'a')
			iMatrx[0][0] = cos(7 * 3.14 / 180), iMatrx[0][2] = sin(7 * 3.14 / 180), iMatrx[2][0] = -sin(7 * 3.14 / 180), iMatrx[2][2] = cos(7 * 3.14 / 180);
		else if (direc == 'd')
			iMatrx[0][0] = cos(-7 * 3.14 / 180), iMatrx[0][2] = sin(-7 * 3.14 / 180), iMatrx[2][0] = -sin(-7 * 3.14 / 180), iMatrx[2][2] = cos(-7 * 3.14 / 180);
		else if (direc == 'w')
			iMatrx[1][1] = cos(7 * 3.14 / 180), iMatrx[1][2] = -sin(7 * 3.14 / 180), iMatrx[2][1] = sin(7 * 3.14 / 180), iMatrx[2][2] = cos(7 * 3.14 / 180);
		else if (direc == 's')
			iMatrx[1][1] = cos(-7 * 3.14 / 180), iMatrx[1][2] = -sin(-7 * 3.14 / 180), iMatrx[2][1] = sin(-7 * 3.14 / 180), iMatrx[2][2] = cos(-7 * 3.14 / 180);
		tempVertices = getCenterVertices();
		
		// composite multiplication of center, translation matrix, rotation, and inverse translation matrices
		CompositeMatrix[0][3] = tempVertices.x, CompositeMatrix[1][3] = tempVertices.y, CompositeMatrix[2][3] = tempVertices.z;
		InvCompositeMatrix[0][3] = -1.0 * (tempVertices.x), InvCompositeMatrix[1][3] = -1.0 * (tempVertices.y), InvCompositeMatrix[2][3] = -1.0 * (tempVertices.z);
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				for (int k = 0; k < 4; ++k)
				{
					tempCompMat[i][j] += iMatrx[i][k] * CompositeMatrix[k][j];
				}
			}
		}
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				for (int k = 0; k < 4; ++k)
				{
					GlobalCompositeMatrix[i][j] += InvCompositeMatrix[i][k] * tempCompMat[k][j];
				}
			}
		}
		// Final composite matrix multiplication of vertex points
		for (size_t i = 0; i < MainVertices.size(); i++) {
			matrixMult(GlobalCompositeMatrix, i);
		}
		// center the model
		centerVertices();
	}
}

// Function to assign the w vertex of the model, gets called when perspective view is called
// performs and produces matrix with z/d
void view(int indx) {
	float inputvert[4] = { MainVertices[indx].x, MainVertices[indx].y, MainVertices[indx].z, 1 };
	float outputvert[4] = { 0, 0, 0, 0 };
	int k = 0;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			outputvert[i] = outputvert[i] + presMat[k++] * inputvert[j];
		}
	}
	//std::cout << "W: " << outputvert[3] << std::endl;
	MainVertices[indx].w = outputvert[3];
}

// If view is orthogonal then set w to 1
void viewOrth(int indx) {
	MainVertices[indx].w = 1;
}

void display(void)   // Create The Display Function
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	      // Clear Screen 

	write_pixel(x_last, y_last, 1.0, 1.0, 1.0);//<-you can get rid of this call if you like
	// CALL YOUR CODE HERE

	//const verticesStruct vS;
	//const vertTextStruct vTS;
	//const normalStruct nS;
	//const faceStruct fS;

	if (viewSwitch == false) {
		for (size_t i = 0; i < MainVertices.size(); i++)
			viewOrth(i);
		for (size_t i = 0; i < faces.size(); i++) {
			// Goes through the faces vector to get the specific face vertices index and looks for that 
			// specific vertex in the mainv vertices vector to draw points between face vertices
			// This performs when view is orthogonal
			midPoint(MainVertices[faces[i].v1 - 1].x, MainVertices[faces[i].v1 - 1].y, MainVertices[faces[i].v2 - 1].x, MainVertices[faces[i].v2 - 1].y);
			midPoint(MainVertices[faces[i].v2 - 1].x, MainVertices[faces[i].v2 - 1].y, MainVertices[faces[i].v3 - 1].x, MainVertices[faces[i].v3 - 1].y);
			midPoint(MainVertices[faces[i].v1 - 1].x, MainVertices[faces[i].v1 - 1].y, MainVertices[faces[i].v3 - 1].x, MainVertices[faces[i].v3 - 1].y);
		}
	}

	else if (viewSwitch == true) {
		for (size_t i = 0; i < MainVertices.size(); i++)
			view(i);
		for (size_t i = 0; i < faces.size(); i++) {
			// Same function as orthogonal view but divides the x and y by z/d to get the perspective view
			float x0 = (MainVertices[faces[i].v1 - 1].x / MainVertices[faces[i].v1 - 1].w);
			float y0 = (MainVertices[faces[i].v1 - 1].y / MainVertices[faces[i].v1 - 1].w);
			float x1 = (MainVertices[faces[i].v2 - 1].x / MainVertices[faces[i].v2 - 1].w);
			float y1 = (MainVertices[faces[i].v2 - 1].y / MainVertices[faces[i].v2 - 1].w);
			float x2 = (MainVertices[faces[i].v3 - 1].x / MainVertices[faces[i].v3 - 1].w);
			float y2 = (MainVertices[faces[i].v3 - 1].y / MainVertices[faces[i].v3 - 1].w);

			// Print the points
			midPoint(x0, y0, x1, y1);
			midPoint(x1, y1, x2, y2);
			midPoint(x0, y0, x2, y2);
		}
	}

	glutSwapBuffers();                                      // Draw Frame Buffer 
}

/***************************************************************************/
void mouse(int button, int state, int x, int y)
{
	/* This function I finessed a bit, the value of the printed x,y should
	   match the screen, also it remembers where the old value was to avoid multiple
	   readings from the same mouse click.  This can cause problems when trying to
	   start a line or curve where the last one ended */
	static int oldx = 0;
	static int oldy = 0;
	int mag;

	y *= -1;  //align y with mouse
	y += 500; //ignore 
	mag = (oldx - x) * (oldx - x) + (oldy - y) * (oldy - y);
	if (mag > 20) {
		printf(" x,y is (%d,%d)\n", x, y);
	}
	oldx = x;
	oldy = y;
	x_last = x;
	y_last = y;
}

/***************************************************************************/
void keyboard(unsigned char key, int x, int y)  // Create Keyboard Function
{

	switch (key) {
	case 27:              // When Escape Is Pressed...
		exit(0);   // Exit The Program
		break;
	case '1':             // stub for new screen
		printf("New screen\n");
		break;
	case 'v':
		std::cout << "view switch: " << viewSwitch << std::endl;
		if (viewSwitch == false) {
			// When view is perspective then center the model for the initial view
			for (size_t i = 0; i < MainVertices.size(); i++) {
				MainVertices[i].x = MainVertices[i].x - centerX;
				MainVertices[i].y = MainVertices[i].y - centerY;
				MainVertices[i].z = MainVertices[i].z - centerZ + maxValue + minValue;
				MainVertices[i].w = MainVertices[i].w - centerW;
			}
			viewSwitch = true;
		}
		break;
	case 't':             // translate left
		translateSwitch = 1;
		scaleSwitch = 0;
		rotSwitch = 0;
		printf("Translating\n");
		break;
	case 'e':             // translate left
		scaleSwitch = 1;
		translateSwitch = 0;
		rotSwitch = 0;
		printf("Scaling\n");
		break;
	case 'r':             // translate left
		rotSwitch = 1;
		translateSwitch = 0;
		scaleSwitch = 0;
		printf("Rotating\n");
		break;
	// wasd direction transformation
	case 'a':
		if (translateSwitch == 1)
			detectTransform(1, 'a');
		else if (scaleSwitch == 1)
			detectTransform(2, 'a');
		else if (rotSwitch == 1)
			detectTransform(3, 'a');
		break;
	case 'd':
		if (translateSwitch == 1)
			detectTransform(1, 'd');
		else if (scaleSwitch == 1)
			detectTransform(2, 'd');
		else if (rotSwitch == 1)
			detectTransform(3, 'd');
		break;
	case 'w':
		if (translateSwitch == 1)
			detectTransform(1, 'w');
		else if (scaleSwitch == 1)
			detectTransform(2, 'w');
		else if (rotSwitch == 1)
			detectTransform(3, 'w');
		break;
	case 's':
		if (translateSwitch == 1)
			detectTransform(1, 's');
		else if (scaleSwitch == 1)
			detectTransform(2, 's');
		else if (rotSwitch == 1)
			detectTransform(3, 's');
		break;
	default:
		break;
	}
}
/***************************************************************************/

int main(int argc, char* argv[])
{
	/* This main function sets up the main loop of the program and continues the
	   loop until the end of the data is reached.  Then the window can be closed
	   using the escape key.						  */

	// read the file
	//for (int i = 0; i < argc; ++i)
	cout << argv[1] << "\n";
	readObjFile(argv[1], vertices, uvs, normals);
	// normalize the model
	normalizeVertices();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Computer Graphics");
	glutDisplayFunc(display);
	glutIdleFunc(display);
	glutMouseFunc(mouse);
	glutKeyboardFunc(keyboard);

	init_window();				             //create_window

	glutMainLoop();                 // Initialize The Main Loop
}
