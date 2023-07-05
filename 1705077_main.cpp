#include <stdio.h>
#include <stdlib.h>

#include <bits/stdc++.h>

#include <ostream>
#include "1705077_rtx_methods.hpp"
using namespace std;

string cameraPosFile = "camera_pos_mirror.txt";
string sceneFile = "scene_mirror.txt";

bool showIntersectionColor = true;

Camera * cam;
Scene * scene;
Screen * screen;
Floor * sceneFloor = nullptr;
SpotLight * spl;
bool printAllZBuffer = true;
int RTX_MAX_RECURSION_LIMIT;
int windowSize;
double windowWidth = 100, windowHeight = 100;

void capture() {
    captureCamPos = cam->position;
    bitmap_image image(screen->width,screen->height);

    Color bgColor = Color(50,50,50);
    for(int i=0;i<screen->height;i++){
        for(int j=0;j<screen->width;j++){
            image.set_pixel(j,i,bgColor.r,bgColor.g,bgColor.b);
        }
    }
    double planeDistance = (windowHeight/2.0) /tan(deg2rad(cam->viewAngle/2.0));
    Vec3 topleft = cam->position
                   + cam->forward*planeDistance
                   - cam->right*windowWidth/2
                   + cam->up*windowHeight/2;
    double du = ((double)windowWidth)/screen->width;
    double dv = ((double)windowHeight)/screen->height;
//    cout<<windowWidth <<" "<< windowHeight << " du " << du << " dv "<< dv<<"  "<<topleft<< cam->right <<endl;
//    cout<<(image.width() * du) << " " << (image.height() * dv) << endl;
    viewPort.clear();
    viewPort.push_back(topleft);
    viewPort.push_back(topleft +  cam->right * image.width() * du);
    viewPort.push_back(topleft +  cam->right * image.width() * du - cam->up * image.height() * dv);
    viewPort.push_back(topleft - cam->up * image.height() * dv);

// Choose middle of the grid cell
    topleft = topleft + cam->right*(0.5*du) - cam->up*0.5*dv;

    double tNear = 0;

    capturePixelEndpoints.clear();
    capturePixelEndpoints.push_back(topleft);
    capturePixelEndpoints.push_back(topleft - cam->up * (image.height()-1) * dv);
    capturePixelEndpoints.push_back(topleft + cam->right * (image.width() - 1) * du);
    capturePixelEndpoints.push_back(topleft + cam->right * (image.width() - 1) * du - cam->up * (image.height()-1) * dv);

    for (int x = 0; x < image.width(); ++x) {
        for (int y = 0; y < image.height(); ++y) {
            Vec3 curPixel = topleft + cam->right * x * du - cam->up * y * dv;
            Ray ray(cam->position, curPixel - cam->position);
            double tMin;
            int nearestObjectIndex = -1;
            int i = 0;
            for (auto & object:scene->getObjects()) {
                double tCur;
                if(object->tryIntersect(ray, tCur)){
                    if(tCur > tNear && (tCur < tMin || nearestObjectIndex < 0)){
                        tMin = tCur;
                        nearestObjectIndex = i;
                    }
//                    intersectionPoints.push_back(ray.getPointAtT(tCur));
                }
                i++;
            }
            if(nearestObjectIndex >= 0){
                Color c;
                scene->getObjects()[nearestObjectIndex]->actuallyIntersect(ray, tMin, c, scene);
                c.clamp();
                image.set_pixel(x, y, c.rInt(), c.gInt(), c.bInt());
            }
        }
    }
    image.save_image("output_1"+ std::to_string(fileNo) +".bmp" );
    //too lazy to flip to new image so just gonna generate two
    image.save_image("output.bmp" );
    fileNo++;
    std::cout<<"!!!!!!!!!!!!!!!!!!!captured!!!!!!!!!!!!!!!!!!!"<<std::endl;
}

int rotateDir = 0;
double rotationSpeed = 5;
double camMoveSpeed = 8;

int debugCaptureLevels = 1;
//int debugCaptureLevel = debugCaptureLevels-1;
bool drawConstantLine = true;
bool drawAllIntersections = false;
bool drawCamBox = true;


Vec3 captureCamPos;
vector<Vec3> capturePixelEndpoints;
vector<Vec3> intersectionPoints;
vector<Vec3> viewPort;
Ray testRay(Vec3(0,0,0), Vec3(0,0,0));
vector<pair<Vec3, Object *>> testRayIntersections;
int testClosestIntersectionIndex=-1;

int drawgrid =1;

void rotateAndSetNormal(double degrees, Vec3 &l,  Vec3 &u, Vec3 &r){
    double angle = (degrees * M_PI) / 180.0;
    l = l * cos(angle) + u * sin(angle);
    u = r.cross(l);

    u = u *  1/u.magnitude();
}

void drawGrid()
{
    int i;
    if(drawgrid==1)
    {
        glColor3f(0.6, 0.6, 0.6);	//grey
        glBegin(GL_LINES);{
            for(i=-8;i<=8;i++){

                if(i==0)
                    continue;	//SKIP the MAIN axes

                //lines parallel to Y-axis
                glVertex3f(i*999, -90, 0);
                glVertex3f(i*999,  90, 0);

                //lines parallel to X-axis
                glVertex3f(-90, i*999, 0);
                glVertex3f( 90, i*999, 0);
            }
        }glEnd();
    }
}


void KeyboardUpHandler(unsigned char key, int x, int y)
{
    switch (key)
    {

        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
            rotateDir = 0;
            break;
    }
}

void keyboardListener(unsigned char key, int x,int y){
    switch(key){
        case '0':
            capture();
            break;
        case '9':
            drawgrid=1-drawgrid;
            break;
            break;
            break;
//
        case '1':
            rotateAndSetNormal(rotationSpeed, cam->right, cam->forward, cam->up);
            break;
        case '2':
            rotateAndSetNormal(-rotationSpeed, cam->right, cam->forward, cam->up);
            break;

        case '3':
            rotateAndSetNormal(rotationSpeed, cam->forward, cam->up, cam->right);
            break;
        case '4':
            rotateAndSetNormal(-rotationSpeed, cam->forward, cam->up, cam->right);
            break;

        case '5':
            rotateAndSetNormal(rotationSpeed, cam->up,cam->right,  cam->forward);
            break;
        case '6':
            rotateAndSetNormal(-rotationSpeed,cam->up, cam->right,  cam->forward);
            break;


        default:
            break;
    }

}

void specialKeyListener(int key, int x,int y){
    switch(key){
        case GLUT_KEY_DOWN:
            cam->position = cam->position + cam->forward * -camMoveSpeed;//down arrow key
            break;
        case GLUT_KEY_UP: {
            cam->position = cam->position + cam->forward * camMoveSpeed;//down arrow key
            break;
        }

        case GLUT_KEY_RIGHT:
            cam->position = cam->position + cam->right * camMoveSpeed;//down arrow key
            break;
        case GLUT_KEY_LEFT:
            cam->position = cam->position + cam->right * -camMoveSpeed;//down arrow key
            break;

        case GLUT_KEY_PAGE_UP:
            cam->position = cam->position + cam->up * camMoveSpeed;//down arrow key
            break;
        case GLUT_KEY_PAGE_DOWN:
            cam->position = cam->position + cam->up * -camMoveSpeed;//down arrow key
            break;

        case GLUT_KEY_INSERT:
            sceneFloor->position.z += .5;
            break;

        case GLUT_KEY_HOME:
            showIntersectionColor = !showIntersectionColor;
            break;
        case GLUT_KEY_F1:
            drawConstantLine = !drawConstantLine;
            break;
        case GLUT_KEY_F2:
            cout<<"curpos, look at, ucamp:\n"<<cam->position <<endl<<(cam->position + cam->forward * 10)<< endl << cam->up.normalized()<<endl;
            break;
        case GLUT_KEY_F3: {
            testRayIntersections.clear();
            testRay.setStart(cam->position);
            testRay.setDir(cam->forward);
            double tMin;
            testClosestIntersectionIndex = -1;
            int  i = 0;

            for (auto & object:scene->getObjects()) {
                double t;
                if(object->tryIntersect(testRay, t)){
                    if(t < tMin || testClosestIntersectionIndex == -1){
                        tMin = t;
                        testClosestIntersectionIndex = i;
                    }
                    testRayIntersections.push_back(make_pair(testRay.getPointAtT(t),object));
                }

                i++;
            }

            if(testClosestIntersectionIndex != -1) {
                auto * o = scene->getObjects()[testClosestIntersectionIndex];
                auto intersectionPOint = testRay.getPointAtT(tMin);
                auto n = o->normalAtPoint(intersectionPOint);
                auto reflected = Vec3::Reflect(testRay.getDir(), n);
                cout << n.angleBetweenInDegrees(reflected) << " og " << n.angleBetweenInDegrees(testRay.getDir())<<endl;
            }
            cout<<"test intersections: "<< testRayIntersections.size()<< " index " << testClosestIntersectionIndex  <<endl;
        }break;
        case GLUT_KEY_F9:
            drawAllIntersections = !drawAllIntersections;
            break;
        case GLUT_KEY_F8:
            drawCamBox = !drawCamBox;
            break;
        case GLUT_KEY_F6:
            camMoveSpeed *= 2;
            break;
        case GLUT_KEY_F5:
            camMoveSpeed /= 2;
            break;
        default:
            break;
    }
}




void mouseListener(int button, int state, int x, int y) {    //x, y is the x-y of the screen (2D)
    switch (button) {
        case GLUT_LEFT_BUTTON:
            break;

        case GLUT_RIGHT_BUTTON:
            break;

        case GLUT_MIDDLE_BUTTON:
            //........
            break;

        default:
            break;

    }
}
void init(){
    glClearColor(0,0,0,0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(cam->viewAngle,	1,	1,	1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

bool printed = false;
void display(){

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.5,0.5,0.5,0.5);	//color black
//    glClearColor(1,1,((x++) % 255) / 255.0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?
    glCamUpdate(cam);
    //gluLookAt(100,100,100,	0,0,0,	0,0,1);
    //gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    //	gluLookAt(0,0,200,	0,0,0,	0,1,0);

    drawGrid();

    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);



    for(auto& object: scene->getObjects()){
//        object->printShape(cout);
        object->draw();
    }
    for (auto & light:scene->getLights()) {
        light->debugDraw();
    }

    glPointSize(10);


    if(testClosestIntersectionIndex >=0){
        glLineWidth(8);
        {
            glColor(Color(1, 0, 1));
            glBegin(GL_LINES);
            glVertex(testRay.getStart());
            glVertex(testRay.getPointAtT(500));
            glEnd();
        }

    }

    if(drawAllIntersections) {
        glBegin(GL_POINTS);
        {
            for (auto p: intersectionPoints) {
//            glColor3d(0,1,0);
                glColor3d(1, 0, 1);
                glVertex(p);
            }
        }
        glEnd();
        for (auto p: intersectionPoints) {
                glVertex(captureCamPos);
                glVertex(p);
        }
        for (auto p: capturePixelEndpoints) {
                glVertex(captureCamPos);
                glVertex(captureCamPos + (p - captureCamPos)*500);
        }
    }
    if(drawCamBox) {
        glBegin(GL_LINES);
        {
            glColor3d(1, 1, 0);
            for (int i = 0; i < viewPort.size(); i++) {
                glVertex(viewPort[i]);
                glVertex(viewPort[(i + 1) % viewPort.size()]);
            }
            glColor3d(1, 0, 1);


        }
        glEnd();
    }
    if(drawConstantLine){
        Ray r(cam->position,cam->forward);
        double tMin;
        double closestIndex = -1;
        int i = 0;
        for (auto & object:scene->getObjects()) {
            double t;
            if(object->tryIntersect(r, t)){
                if(t < tMin || closestIndex == -1){
                    tMin = t;
                    closestIndex = i;
                }
            }

            i++;
        }
//        cout<<closestIndex<<endl;
        if(closestIndex != -1) {
            Vec3 intersectionPoint = r.getPointAtT(tMin);
            Object * o = scene->getObjects()[closestIndex];
            Vec3 n = o->normalAtPoint(intersectionPoint);
            Vec3 reflected = Vec3::Reflect(r.getDir(), n);
//            cout<<r.getDir()<<" refl: "<<reflectedFrom<<endl;
            {
                glColor(Color(1, 0, 1) * .5f);
                glBegin(GL_LINES);
//                glVertex(cam->position);
//                glVertex(intersectionPoint);
                glVertex(intersectionPoint);
                glVertex(intersectionPoint + reflected * 8);
                glEnd();
            }
            {
                glColor(o->colorAtPoint(intersectionPoint).inverted());
                glBegin(GL_LINES);
                glVertex(intersectionPoint);
                glVertex(intersectionPoint + n * 80);
                glEnd();

                glBegin(GL_POINTS);
                glVertex(intersectionPoint);
                glEnd();
            }

//            spl->plsDeleteMe = spl->dotprodcheck( intersectionPoint) ? Color::Green():Color::Red() ;
        }
    }
    glLineWidth(.1f);

    glutSwapBuffers();
}
void animate(){
    //codes for any changes in Models, Camera
    if(rotateDir != 0){
        rotateAndSetNormal(rotateDir * rotationSpeed, cam->forward,cam->up,cam->right);
    }

    glutPostRedisplay();
}


bool tryLoadData() {

    ifstream fin(sceneFile);
    if(!fin.is_open()){
        cout<<"No scene file"<<endl;
        return false;
    }
    fin >> RTX_MAX_RECURSION_LIMIT;
    fin >> windowSize;
    cout
    << RTX_MAX_RECURSION_LIMIT << endl
    << windowSize<<endl;
    windowHeight = windowWidth = windowSize;
    screen = new Screen(windowSize,windowSize);

    Vec3 CamPos(-1.9230, -59.7890, 55.3205);
    Vec3 camLookPos(-1.9230, -51.7572, 49.3632);
    Vec3 up=Vec3::Up();
    {
        ifstream customFin(cameraPosFile);
        if(customFin.is_open()) {
            customFin >> CamPos;
            customFin >> camLookPos;
            customFin >> up;
            cout<<"CamPos   "<<CamPos<<endl;
            cout<<"camLookPos   "<<camLookPos<<endl;
            cout<<"up   "<<up<<endl;
        }
        customFin.close();
    }
    cam = new Camera(CamPos,camLookPos,up);

    cam->print(cout);
    scene = new Scene();
    sceneFloor = new Floor(Vec3(0,0));
    scene->addObject(sceneFloor);

    string command;
    int numObjects = 0;
    fin>> numObjects;
    cout<<numObjects<<" x objects\n";
    while(numObjects--){
        fin >> command;
        cout<<"############################"<<endl;
        cout<<command<<endl;
        cout<<"############################"<<endl;
        if(command == "triangle"){
            Vec3 pos = Vec3::Zero();
            Vec3 a,b,c;
            fin>> a >> b>>c;
            Material mat;
            fin>>mat;
            Triangle *  t = new Triangle (pos, a,b,c, mat);
            cout <<(*t);
            scene->addObject(t);
        }
        else if(command == "sphere"){
            Vec3 pos;
            fin>> pos;
            double r;
            fin>>r;
            Material mat;
            fin>>mat;
            Sphere * s = new Sphere(pos, r, mat);
            scene->addObject(s);
        }else if(command =="quad"){
            Vec3 pos;
            fin>>pos;
            Vec3 normal,right, up;
            fin>>normal>>right>> up;
            double width, height;
            fin>>width>>height;
            Material mat;
            fin >> mat;
//            mat.alpha = .6f;
            Quad * quad = new Quad(pos,
                                   normal, right, up, width, height, mat);
            cout<<(*quad)<<endl;
            cout<<mat<<endl;
            scene->addObject(quad);
        }
        else if(command =="general"){
            double a,b,c,d,e,f,g,h,i,j;
            fin >> a>> b>> c>> d>> e>> f>> g>> h>> i>> j;
            Vec3 pos;
            fin>>pos;
            double clipCubeLength,clipCubeWidth, clipCubeHeight;
            fin>>clipCubeLength>>clipCubeWidth>> clipCubeHeight;
            Material mat;
            fin >> mat;
            GeneralQuadraticSurface * qs = new GeneralQuadraticSurface(pos,mat,
                                                                     a,b,c,d,e,f,g,h,i,j,
                                                                     clipCubeLength,clipCubeWidth, clipCubeHeight);
            cout<< (*qs)<<endl;
            cout<<mat<<endl;
            scene->addObject(qs);
        }
    }
    int numPointLights;
    fin >> numPointLights;
    while (numPointLights--){
        Vec3 pos;
        fin>>pos;
        Color c;
        fin>>c;
        PointLight * pLight = new PointLight(pos, c);
        scene->addLight(pLight);
    }
    int numSpotLights;
    fin >> numSpotLights;
    while (numSpotLights--){
        Vec3 pos;
        fin>>pos;
        Color c;
        fin>>c;
        Vec3 dir;
        fin>> dir;
        double cutoffAngle;
        fin>>cutoffAngle;
        SpotLight * spotLight = new SpotLight(pos,c,dir,cutoffAngle);
        spl = spotLight;
        scene->addLight(spotLight);
    }
    fin.close();

    return true;
}

int main(int argc, char **argv){
    if(!tryLoadData()){
        cout<<"Data loading failed"<<endl;
        return -1;
    }

    glutInit(&argc,argv);

    glutInitWindowSize(windowSize, windowSize);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("RayTracing");


    init();


    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)


    glutKeyboardFunc(keyboardListener);
    glutKeyboardUpFunc(KeyboardUpHandler);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();

    delete cam;
    delete scene;
    return 0;
}






