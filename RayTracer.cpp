/*==================================================================================
* COSC 363  Computer Graphics
* Department of Computer Science and Software Engineering, University of Canterbury.
*
* Artemis Hingston
* irh20
* Assignment 2
* See Report   for details.
*===================================================================================
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include "Sphere.h"
#include "SceneObject.h"
#include "Ray.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Cone.h"
#include "TextureBMP.h"
#include <GL/freeglut.h>
#include <glm/gtc/random.hpp>
using namespace std;

// Parameters for testing
const bool yesRecursion = false;
const float numRecursions = 3.0;
const float yesFog = false;


const float EDIST = 40.0; 	// Dist of image plane from camera
const int NUMDIV = 500;		// Subdivisions along x & y
const int MAX_STEPS = 5;	// levels of recursion
const float XMIN = -10.0;	// Boundary values of the image place 
const float XMAX = 10.0;	// defined so the view axis passes thru centre
const float YMIN = -10.0;
const float YMAX = 10.0;

const float tableY = -15.0;
const float tableMinX = -20.0;
const float tableMaxX = 20.0;
const float tableMinZ = -200.0;
const float tableMaxZ = -100.0;
const float legHeight = 30.0;
const float legInset = 5.0;
const float legY = tableY - legHeight;


// Box Coords
glm::vec3 NETop = glm::vec3(-150.0, 100.0, 200.0);
glm::vec3 NEBase = glm::vec3(-150.0, -50.0, 200.0);
glm::vec3 NWTop = glm::vec3(150.0, 100.0, 200.0);
glm::vec3 NWBase = glm::vec3(150.0, -50.0, 200.0);
glm::vec3 SWTop = glm::vec3(150.0, 100.0, -1000.0);
glm::vec3 SWBase = glm::vec3(150.0, -50.0, -1000.0);
glm::vec3 SETop = glm::vec3(-150.0, 100.0, -1000.0);
glm::vec3 SEBase = glm::vec3(-150.0, -50.0, -1000.0);

TextureBMP texture;
TextureBMP texId[3];

vector<SceneObject*> sceneObjects;


//---The most important function in a ray tracer! ---------------------------------- 
//   Computes the colour value obtained by tracing a ray and finding its 
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step)
{
	glm::vec3 backgroundCol(0.5);						//Background colour = (0,0,0)
	glm::vec3 lightPos(10, 40, -3);					//Light's position
	glm::vec3 color(0);
	SceneObject* obj;

    ray.closestPt(sceneObjects);					//Compare the ray with all objects in the scene
    if(ray.index == -1) return backgroundCol;		//no intersection
	obj = sceneObjects[ray.index];					//object on which the closest point of intersection is found

	if (ray.index == 11) {
		// Chessboard pattern
        int squareSize = 2;
        int ix = floor(ray.hit.x / squareSize);
        int iz = floor(ray.hit.z / squareSize);
        
        if ((ix + iz) % 2 == 0) {
            color = glm::vec3(0, 0, 0); // Black
        } else {
            color = glm::vec3(1, 1, 1); // White
        }
        obj->setColor(color);
	}

	if (ray.index == 6) {
		// Texture Mapping for Table top
		float texcoords = (ray.hit.x - -20.0)/(20.0 - -20.0);	// x1=-15, x2=5
		float texcoordt = (ray.hit.z - -200.0)/(-100.0 - -200.0);	// z1=-60, z2=-90
		if (texcoords > 0 && texcoords < 1 && texcoordt > 0 && texcoordt < 1) {
			color = texId[0].getColorAt(texcoords, texcoordt);
			obj->setColor(color);
		}
	}

	if (ray.index == 0) {
		// Texture Mapping for roof
		float texcoords = (ray.hit.x - -150.0)/(150.0 - -150.0);	// x1=-15, x2=5
		float texcoordt = (ray.hit.z - -1000.0)/(40.0 - -1000.0);	// z1=-60, z2=-90
		if (texcoords > 0 && texcoords < 1 && texcoordt > 0 && texcoordt < 1) {
			color = texId[1].getColorAt(texcoords, texcoordt);
			obj->setColor(color);
		}
	}

	if (ray.index == 16) {
		//Texture mapping for globe
		glm::vec3 normalisedVec = obj->normal(ray.hit);
		float texcoords = 0.5 + (std::atan2(normalisedVec.z, normalisedVec.x))/(2*M_PI);
		float texcoordt = 0.5 + (std::asin(normalisedVec.y)/(M_PI));
		color = texId[2].getColorAt(texcoords, texcoordt);
			obj->setColor(color);
	}

	if (ray.index == 17) {
		//Texture mapping for cylinder
		glm::vec3 normalisedVec = obj->normal(ray.hit);
		float texcoords = 0.5 + (std::atan2(normalisedVec.z, normalisedVec.x))/(2*M_PI);
		float texcoordt = (ray.hit.y - -15.0)/(-7.0 - -15.0);
		color = texId[1].getColorAt(texcoords, texcoordt);
			obj->setColor(color);
	}


	// if (ray.index == 1) {
	// 	// Texture Mapping for ground
	// 	float texcoords = (ray.hit.x - -150.0)/(150.0 - -150.0);	// x1=-15, x2=5
	// 	float texcoordt = (ray.hit.z - -1000.0)/(40.0 - -1000.0);	// z1=-60, z2=-90
	// 	if (texcoords > 0 && texcoords < 1 && texcoordt > 0 && texcoordt < 1) {
	// 		color = texId[2].getColorAt(texcoords, texcoordt);
	// 		obj->setColor(color);
	// 	}
	// }


	color = obj->lighting(lightPos ,-ray.dir, ray.hit);		//Object's colour
	glm::vec3 lightVec = lightPos - ray.hit;
	Ray shadowRay(ray.hit, lightVec);

	float lightDist = glm::length(lightVec);

	// Shadows
	//TODO Loop to get mult random lights & then average

	// loop thru 20ish times 

	shadowRay.closestPt(sceneObjects);

	if ((shadowRay.index > -1) && (shadowRay.dist < lightDist)) {
        SceneObject* shadowObj = sceneObjects[shadowRay.index];
        if (shadowObj->isTransparent() || shadowObj->isRefractive()) {
			color = 0.4f * obj->getColor();	//Lighter for transparent and refractive
        } else {
            color = 0.2f * obj->getColor();
        }
    }

	// Handle reflection
	if (obj->isReflective() && step < MAX_STEPS) {
		float rho = obj->getReflectionCoeff();
		glm::vec3 normalVec = obj->normal(ray.hit);
		glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVec);
		Ray reflectedRay(ray.hit, reflectedDir);
		glm::vec3 reflectedColor = trace(reflectedRay, step + 1);
		color = color + (rho * reflectedColor);
	}

	// Handle refraction
    if (obj->isRefractive() && step < MAX_STEPS) {
        float eta = obj->getRefractiveIndex();
        glm::vec3 normalVec = obj->normal(ray.hit);
        glm::vec3 g = glm::refract(ray.dir, normalVec, 1.0f / eta);
        Ray refrRay(ray.hit, g);
        refrRay.closestPt(sceneObjects);
        if (refrRay.index != -1) {
            glm::vec3 m = sceneObjects[refrRay.index]->normal(refrRay.hit);
            glm::vec3 h = glm::refract(g, -m, eta);
            Ray exitRay(refrRay.hit, h);
            glm::vec3 refractedColor = trace(exitRay, step + 1);
            float refrCoeff = obj->getRefractionCoeff();
            color = (1 - refrCoeff) * color + (refrCoeff * refractedColor);
        }
    }

    // Handle transparency
    if (obj->isTransparent() && step < MAX_STEPS) {
        float tranCoeff = obj->getTransparencyCoeff();
        Ray transparentRay(ray.hit, ray.dir);
        glm::vec3 transparentColor = trace(transparentRay, step + 1);
        color = (1 - tranCoeff) * color + (tranCoeff * transparentColor);
    }

	if (yesFog == true) {
		float z1 = -40;
		float z2 = -900;
		float gamma = ((ray.hit.z)-z1)/(z2-z1);
		color = (1-gamma)*color + gamma +  glm::vec3(0, 0, 0);
	}

	return color;
}


// All the objects and things

void loadTextures() {
	texId[0] = TextureBMP("../oak.bmp");
	texId[1] = TextureBMP("../roof.bmp");
	texId[2] = TextureBMP("../Earth.bmp");
}


void loadWalls() {

	Plane *roof = new Plane(SETop, SWTop, NWTop, NETop);
	roof->setColor(glm::vec3(0, 1, 1));
	roof->setSpecularity(false);
	sceneObjects.push_back(roof);

	Plane *ground = new Plane(NEBase, NWBase, SWBase, SEBase);
	ground->setColor(glm::vec3(0.5, 1, 0.5));
	ground->setSpecularity(false);
	sceneObjects.push_back(ground);

	Plane *northWall = new Plane(NETop, NWTop, NWBase, NEBase);
	northWall->setColor(glm::vec3(1, 0, 0));
	northWall->setSpecularity(false);
	sceneObjects.push_back(northWall);

	Plane *southWall = new Plane(SEBase, SWBase, SWTop, SETop);
	southWall->setColor(glm::vec3(0.89, 1, 0.16));
	southWall->setSpecularity(false);
	sceneObjects.push_back(southWall);

	Plane *westWall = new Plane(SWBase, NWBase, NWTop, SWTop);
	westWall->setColor(glm::vec3(0, 0, 1));
	westWall->setSpecularity(false);
	sceneObjects.push_back(westWall);

	Plane *eastWall = new Plane(NEBase, SEBase, SETop, NETop);
	eastWall->setColor(glm::vec3(1, 0, 1));
	eastWall->setSpecularity(false);
	sceneObjects.push_back(eastWall);
}



void loadTable() {
	glm::vec3 tableColour = glm::vec3(0.4, 0.26, 0.09);

	Plane *tableTop = new Plane(glm::vec3(tableMinX, tableY, tableMaxZ), glm::vec3(tableMaxX, tableY, tableMaxZ), glm::vec3(tableMaxX, tableY, tableMinZ), glm::vec3(tableMinX, tableY, tableMinZ));
	tableTop->setColor(glm::vec3(0.8, 0.8, 0));
	tableTop->setSpecularity(false);
	sceneObjects.push_back(tableTop);

	Cylinder *legNE = new Cylinder(glm::vec3(tableMinX+legInset, legY, tableMaxZ-legInset), 1.5, legHeight);
	legNE->setColor(tableColour);
	sceneObjects.push_back(legNE);

	Cylinder *legNW = new Cylinder(glm::vec3(tableMaxX-legInset, legY, tableMaxZ-legInset), 1.5, legHeight);
	legNW->setColor(tableColour);
	sceneObjects.push_back(legNW);

	Cylinder *legSW = new Cylinder(glm::vec3(tableMaxX-legInset, legY, tableMinZ+legInset), 1.5, legHeight);
	legSW->setColor(tableColour);
	sceneObjects.push_back(legSW);

	Cylinder *legSE = new Cylinder(glm::vec3(tableMinX+legInset, legY, tableMinZ+legInset), 1.5, legHeight);
	legSE->setColor(tableColour);
	sceneObjects.push_back(legSE);
}

void loadChessBoard() {
	Plane *chessBoard = new Plane(glm::vec3(2.0, -14.0, -102.0), glm::vec3(18.0, -14.0, -102.0), glm::vec3(18.0, -14.0, -118.0), glm::vec3(2.0, -14.0, -118.0));
	chessBoard->setColor(glm::vec3(1));
	chessBoard->setSpecularity(false);
	sceneObjects.push_back(chessBoard);

	Plane *chessBoardSide1 = new Plane(glm::vec3(2.0, -15.0, -102.0), glm::vec3(18.0, -15.0, -102.0), glm::vec3(18.0, -14.0, -102.0), glm::vec3(2.0, -14.0, -102.0));
	chessBoardSide1->setColor(glm::vec3(0.24, 0.17, 0.12));
	chessBoardSide1->setSpecularity(false);
	sceneObjects.push_back(chessBoardSide1);

	Plane *chessBoardSide2 = new Plane(glm::vec3(2.0, -15.0, -118.0), glm::vec3(18.0, -15.0, -118.0), glm::vec3(18.0, -14.0, -118.0), glm::vec3(2.0, -14.0, -118.0));
	chessBoardSide2->setColor(glm::vec3(0.24, 0.17, 0.12));
	chessBoardSide2->setSpecularity(false);
	sceneObjects.push_back(chessBoardSide2);

	Plane *chessBoardSide3 = new Plane(glm::vec3(2.0, -15.0, -102.0), glm::vec3(2.0, -15.0, -118.0), glm::vec3(2.0, -14.0, -118.0), glm::vec3(2.0, -14.0, -102.0));
	chessBoardSide3->setColor(glm::vec3(0.24, 0.17, 0.12));
	chessBoardSide3->setSpecularity(false);
	sceneObjects.push_back(chessBoardSide3);

	Plane *chessBoardSide4 = new Plane(glm::vec3(18.0, -15.0, -118.0), glm::vec3(18.0, -15.0, -102.0), glm::vec3(18.0, -14.0, -102.0), glm::vec3(18.0, -14.0, -118.0));
	chessBoardSide4->setColor(glm::vec3(0.24, 0.17, 0.12));
	chessBoardSide4->setSpecularity(false);
	sceneObjects.push_back(chessBoardSide4);

}

void loadMirror() {
	Plane *mirror = new Plane(glm::vec3(-60.0, -20.0, -250.0), glm::vec3(30.0, -20.0, -300.0), glm::vec3(30.0, 25.0, -290.0), glm::vec3(-60.0, 25.0, -240.0));
	mirror->setColor(glm::vec3(0.1, 0.1, 0.1));
	mirror->setReflectivity(true, 1.0);
	// mirror->setSpecularity(false);
	sceneObjects.push_back(mirror);
	
}

void loadCupAndSaucer() {
	Cylinder *cylinder1 = new Cylinder(glm::vec3(-5.0, -15.0, -115.0), 2.0, 5.0);
	cylinder1->setColor(glm::vec3(0.9, 0.9, 1));
	cylinder1->setTransparency(true, 0.5);
	sceneObjects.push_back(cylinder1);

	Cylinder *checker1 = new Cylinder(glm::vec3(15.0, -14.0, -104.0), 1.0, 1.0);
	checker1->setColor(glm::vec3(0.6, 0.2, 6));
	// cylinder1->setTransparency(true, 0.5);
	sceneObjects.push_back(checker1);

	Cylinder *checker2 = new Cylinder(glm::vec3(11.0, -14.0, -108.0), 1.0, 1.0);
	checker2->setColor(glm::vec3(0.9, 0.3, 0.6));
	// cylinder1->setTransparency(true, 0.5);
	sceneObjects.push_back(checker2);

	Cylinder *checker3 = new Cylinder(glm::vec3(9.0, -14.0, -112.0), 1.0, 1.0);
	checker3->setColor(glm::vec3(0.9, 0.3, 0.6));
	// cylinder1->setTransparency(true, 0.5);
	sceneObjects.push_back(checker3);

	Cylinder *checker4 = new Cylinder(glm::vec3(5.0, -14.0, -106.0), 1.0, 1.0);
	checker4->setColor(glm::vec3(0.6, 0.2, 6));
	// cylinder1->setTransparency(true, 0.5);
	sceneObjects.push_back(checker4);

	// Cone *coneBase = new Cone(glm::vec3(0.0, -10.0, -115.0), 5.0, 10.0);
	// coneBase->setColor(glm::vec3(1.0, 1.0, 1.0));
	// coneBase->setTransparency(true, 0.8);
	// sceneObjects.push_back(coneBase);

	// Sphere *ball1 = new Sphere(glm::vec3(1.0, -5.0, -115.0), 1.0);
	// ball1->setColor(glm::vec3(1.0, 0.0, 1.0));
	// sceneObjects.push_back(ball1);

	// Sphere *ball2 = new Sphere(glm::vec3(-2.0, -3.0, -117.0), 2.0);
	// ball2->setColor(glm::vec3(0.0, 0.0, 1.0));
	// sceneObjects.push_back(ball2);
}


void loadHourGlass() {
	Cone *cone1 = new Cone(glm::vec3(-10.0, -15.0, -150.0), 5.0, 10.0);
	cone1->setColor(glm::vec3(1, 0, 0));
	// cone1->setTransparency(true, 0.5);
	sceneObjects.push_back(cone1);
}

void loadCrystalBall() {
	Sphere *glassSphere = new Sphere(glm::vec3(-15.0, -2.0, -150.0), 5.0);
	glassSphere->setColor(glm::vec3(0.0, 0.0, 0.0));
    glassSphere->setReflectivity(true, 0.1f);         // Low reflectivity
    glassSphere->setRefractivity(true, 1.0, 1.5);   // High refractivity, with refractive index 1.5
    // glassSphere->setTransparency(true, 0.8);         // High transparency
	sceneObjects.push_back(glassSphere);

	Cone *coneBase = new Cone(glm::vec3(-15.0, -15.0, -150.0), 5.0, 10.0);
	coneBase->setColor(glm::vec3(1.0, 0.0, 1.0));
	// coneBase->setTransparency(true, 0.9);
	sceneObjects.push_back(coneBase);
}



void loadGlobe() {
	Sphere *sphere1 = new Sphere(glm::vec3(10.0, -5.0, -190.0), 5.0);
	sphere1->setColor(glm::vec3(0.25, 0.88, 0.82));
	// sphere1->setTransparency(true, 0.5);
	sceneObjects.push_back(sphere1);

	Cylinder *cylinder = new Cylinder(glm::vec3(10.0, -15.0, -190.0), 5.0, 5.0);
	cylinder->setColor(glm::vec3(0, 0, 0));
	sceneObjects.push_back(cylinder);
}

void loadOtherThings() {
	Sphere *bigShinySphere = new Sphere(glm::vec3(40.0, -15.0, -240.0), 20.0);
	bigShinySphere->setColor(glm::vec3(0, 0, 0));
	// sphere1->setSpecularity(false);
	bigShinySphere->setReflectivity(true, 0.8);
	sceneObjects.push_back(bigShinySphere);

	// Cylinder *tallPole = new Cylinder(glm::vec3(0.0, -50.0, -250.0), 10.0, 100.0);
	// tallPole->setColor(glm::vec3(0.5, 0.5, 0.5));
	// sceneObjects.push_back(tallPole);

	Cylinder *smallCylinder = new Cylinder(glm::vec3(-15.0, -7.0, -150.0), 5.0, 2.0);
	smallCylinder->setColor(glm::vec3(1.0, 0.5, 1.0));
	sceneObjects.push_back(smallCylinder);
}


// Function to compute the color difference between two colors
float colorDifference(const glm::vec3& col1, const glm::vec3& col2) {
    return glm::length(col1 - col2);
}

// Recursive adaptive sampling function
glm::vec3 adaptiveSample(float xp, float yp, float cellX, float cellY, int depth, float threshold, glm::vec3 eye) {
    if (depth == 0) {
        glm::vec3 dir(xp + 0.5f * cellX, yp + 0.5f * cellY, -EDIST);
        Ray ray(eye, dir);
        return trace(ray, 1);
    }

    // Sample four corners of the pixel
    glm::vec3 colors[4];
    colors[0] = adaptiveSample(xp, yp, cellX / 2, cellY / 2, depth - 1, threshold, eye);
    colors[1] = adaptiveSample(xp + cellX / 2, yp, cellX / 2, cellY / 2, depth - 1, threshold, eye);
    colors[2] = adaptiveSample(xp, yp + cellY / 2, cellX / 2, cellY / 2, depth - 1, threshold, eye);
    colors[3] = adaptiveSample(xp + cellX / 2, yp + cellY / 2, cellX / 2, cellY / 2, depth - 1, threshold, eye);

    // Compute the average color
    glm::vec3 avgColor = (colors[0] + colors[1] + colors[2] + colors[3]) / 4.0f;

    // Check if the color difference exceeds the threshold
    if (colorDifference(colors[0], avgColor) > threshold ||
        colorDifference(colors[1], avgColor) > threshold ||
        colorDifference(colors[2], avgColor) > threshold ||
        colorDifference(colors[3], avgColor) > threshold) {
        // Subdivide further
        avgColor = (adaptiveSample(xp, yp, cellX / 2, cellY / 2, depth - 1, threshold, eye) +
                    adaptiveSample(xp + cellX / 2, yp, cellX / 2, cellY / 2, depth - 1, threshold, eye) +
                    adaptiveSample(xp, yp + cellY / 2, cellX / 2, cellY / 2, depth - 1, threshold, eye) +
                    adaptiveSample(xp + cellX / 2, yp + cellY / 2, cellX / 2, cellY / 2, depth - 1, threshold, eye)) / 4.0f;
    }

    return avgColor;
}


//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display()
{
	float xp, yp;  //grid point
	float cellX = (XMAX - XMIN) / NUMDIV;  //cell width
	float cellY = (YMAX - YMIN) / NUMDIV;  //cell height
	glm::vec3 eye(0., 0., 0.);

	glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	glBegin(GL_QUADS);  //Each cell is a tiny quad.

	for (int i = 0; i < NUMDIV; i++)	//Scan every cell of the image plane
	{
		xp = XMIN + i * cellX;
		for (int j = 0; j < NUMDIV; j++)
		{
			yp = YMIN + j * cellY;
			glm::vec3 col(0);

			if (!yesRecursion) {
				// No anti-aliasing
				glm::vec3 dir(xp + 0.5 * cellX, yp + 0.5 * cellY, -EDIST);	//direction of the primary ray
				Ray ray = Ray(eye, dir);
				col += trace(ray, 1); //Trace the primary ray and get the colour value
			} else {
				// Adaptive anti-aliasing
            	col += adaptiveSample(xp, yp, cellX, cellY, numRecursions, 0.1f, eye);
			}

			glColor3f(col.r, col.g, col.b);
			glVertex2f(xp, yp);			//Draw each cell with its color value
			glVertex2f(xp + cellX, yp);
			glVertex2f(xp + cellX, yp + cellY);
			glVertex2f(xp, yp + cellY);
		}
	}

    glEnd();
    glFlush();
}



//---This function initializes the scene ------------------------------------------- 
//   Specifically, it creates scene objects (spheres, planes, cones, cylinders etc)
//     and add them to the list of scene objects.
//   It also initializes the OpenGL 2D orthographc projection matrix for drawing the
//     the ray traced image.
//----------------------------------------------------------------------------------
void initialize()
{
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);

    glClearColor(0, 0, 0, 1);

	loadTextures();

	loadWalls();
	loadTable();
	loadChessBoard();
	loadGlobe();

	loadCupAndSaucer();
	// loadHourGlass();
	loadMirror();
	loadCrystalBall();
	loadOtherThings();
	// loadSpheres();
	
}


int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(800, 800);
    glutInitWindowPosition(20, 20);
    glutCreateWindow("Raytracing");

    glutDisplayFunc(display);
    initialize();

    glutMainLoop();
    return 0;
}
