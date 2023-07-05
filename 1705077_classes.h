//
// Created by PC on 8/10/2022.
//

#ifndef DEF_1705077_CLASSES_H
#define DEF_1705077_CLASSES_H
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include <windows.h>
#include <GL/glut.h>
#include <ostream>
#define pi (2*acos(0.0))
extern int RTX_MAX_RECURSION_LIMIT;
double RTX_EPSILON = sqrt(1.19209290 * pow(10,-7));
bool LIGHT_OPPOSITE_FACE_CHECK =false;
bool ALLOW_BACKFACE_LIGHTING =true;
double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}
class Vec3;
class Scene;
class Object;
class Camera;
class Screen;
class Floor;
class SpotLight;
class Ray;


class Color{
public:
    double r;
    double g;
    double b;

    friend std::ostream &operator<<(std::ostream &os, const Color &color) {
        os << color.r << " " << color.g << " " << color.b;
        return os;
    }

    friend Color operator * (const Color c, const double factor) {
        return Color(c.r * factor, c.g * factor, c.b * factor);
    }

    friend Color operator * (const Color c1, const Color c2) {
        return Color(c1.r * c2.r, c1.g * c2.g, c1.b * c2.b);
    }

    friend Color operator + (const Color c1,const Color c2) {
        return Color(c1.r + c2.r, c1.g + c2.g, c1.b + c2.b);
    }
    friend Color operator - (const Color c1,const Color c2) {
        return Color(c1.r - c2.r, c1.g - c2.g, c1.b - c2.b);
    }

    friend std::istream &operator >>(std::istream &is, Color &c) {
        double r,g,b;
        is >> c.r >> c.g >> c.b;

        return is;
    }


    Color(double r= 0, double g= 0, double b= 0) : r(r), g(g), b(b) {}

    static Color White() {return Color(1,1,1);}
    static Color Black() {return Color(0,0,0);}
    static Color Violet() {return Color(148/255.0, 0/255.0, 211/255.0);}
    static Color Indigo() {return Color(75/255.0, 0/255.0, 130/255.0);}
    static Color Yellow() {return Color(255/255.0, 255/255.0, 0/255.0);}
    static Color Orange() {return Color(255/255.0, 127/255.0, 0/255.0);}
    static Color Blue()   {return Color(0,0,1);}
    static Color Green()  {return Color(0,1,0);}
    static Color Red() {return Color(1,0,0);}

    int rInt() {return (int)(r *255);}
    int gInt() {return (int)(g *255);}
    int bInt() {return (int)(b *255);}

    bool allZero(){return r == 0 && g == 0 && b == 0;}

    Color lerpTo(Color target, double t){
        Color result(r,g,b);
        return result + (target - result) * t;
    }
    Color inverted()const {return Color(1.0-r, 1.0-g, 1.0-b); }

    void clamp(){
        r =  r < 0
             ? 0
             : r > 1
               ? 1
               : r;

        g =  g < 0
             ? 0
             : g > 1
               ? 1
               : g;

        b =  b < 0
             ? 0
             : b > 1
               ? 1
               : b;

    }
};
/**
 * Base light class
 */
class Light{
public:
    Color color;

    Light(const Color &color) : color(color) {}
    virtual Vec3 getOrigin() =0;

    /// Draw the light's representation in the sceneview
    virtual void debugDraw() =0;
    /// returns whether a ray from this light source to a given point hits a given object
    /// \param scene The scene we raycast in
    /// \param viewRay the ray we raycast in
    /// \param normalAtPoint used to determine whether we hit backface
    /// \param object the object we want to verify that we hit
    /// \param point the point to which we raycast
    /// \return
    virtual bool hitsObjectAtPointInScene(Scene *scene, const Ray &viewRay,
                                          const Vec3 &normalAtPoint,
                                          Object *object,
                                          const Vec3 &point) =0;

    /// Returns the dir this lightsource is casting towards a point
    /// \param point the point
    /// \return the dir this lightsource is casting towards a point
    virtual Vec3 getLightDirForPoint(const Vec3 &point) = 0;
};
/**
 * Scene object
 * Contains the objects and lights in scene
 */
class Scene{
    std::vector<Object *> objects;
    std::vector<Light *> lights;

public:
    Color ambientLighting = Color::White();
    const std::vector<Object *> & getObjects() const;

    const std::vector<Light *> &getLights() const;


    void addObject(Object* object);

    void addLight(Light* light);

    ~Scene(){
        for (auto& s:objects) {
            delete(s);
        }
    }
};
/**
 * Custom generic matric class
 * @tparam Rows
 * @tparam Cols
 */
template<int Rows,int Cols>
class Matrix{
protected:
    double array[Rows][Cols];
public:
    double& operator()(int i, int j)
    {
        return array[i][j];
    }

    friend std::ostream &operator <<(std::ostream &os, const Matrix<Rows,Cols> &mat) {
//        auto oldPrecision = cout.precision(12);
        for(int r =0; r < Rows;r++){
            for(int c =0; c < Cols;c++){
                os << std::setfill(' ') << std::setw(12)
                   << std::fixed << std::setprecision(4)
                   << mat.array[r][c] << "\t";
            }
            os<< std::endl;
        }
//        cout.precision(oldPrecision);
        return os;
    }



    template<int FinalCols>
    inline friend Matrix<Rows,FinalCols> operator * (Matrix<Rows,Cols> &m1, Matrix<Cols,FinalCols> &m2) {
        Matrix<Rows,FinalCols> result;
        for(int i = 0; i < Rows; i++){
            for(int k = 0; k < FinalCols; k++){
                result(i,k) = 0;

                for (int j = 0; j < Cols; j++) {
                    result(i,k) += m1.array[i][j] * m2(j,k);
                }
            }
        }
        return result;
    }

    inline friend Matrix<Rows,Cols> operator * (Matrix<Rows,Cols> &m1, double factor ) {
        Matrix<Rows,Cols> result;
        for(int r =0; r < Rows;r++){
            for(int c =0; c < Cols;c++){
                result.array[r][c] = m1.array[r][c] * factor;
            }
        }
        return result;
    }


    inline friend Matrix<Rows,Cols> operator + (Matrix<Rows,Cols> &m1, Matrix<Rows,Cols> &m2) {
        Matrix<Rows,Cols> result;
        for(int r =0; r < Rows;r++){
            for(int c =0; c < Cols;c++){
                result.array[r][c] = m1.array[r][c] + m2.array[r][c];
            }
        }
        return result;
    }

    Matrix<Rows,Cols> & operator+=(const Matrix<Rows,Cols> & rhs){
        for(int r =0; r < Rows;r++){
            for(int c =0; c < Cols;c++){
                array[r][c] = array[r][c] + rhs.array[r][c];
            }
        }
        return *this;
    }

};
typedef Matrix<4,4> Mat4x4;
//we only need 3x3 determinant
double determinant(Matrix<3, 3> &m){
    double det = 0;
    for (int c = 0; c < 3; ++c) {
        int c1 = (c + 1 +3)% 3;
        int c2 = (c - 1 +3)% 3;
        if(c1 > c2){
            std::swap(c1,c2);
        }
        det += m(0,c)
                * (m(1,c1) * m(2,c2) - m(2,c1) * m(1,c2))
                * ((c%2) != 0? -1 : 1);
    }
    return det;
}


/**
 * Custom Vector3 class
 */
class Vec3{
public:
    double x, y,z;
    Vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    Vec3(const Vec3& v): x(v.x), y(v.y),z(v.z){  }
    Vec3(Matrix<4,1> &m): x(m(0,0)), y(m(1,0)),z(m(2,0)){}

    friend std::ostream &operator<<(std::ostream &os, const Vec3 &v) {
        os << v.x << " "<< v.y<< " " << v.z;
        return os;
    }

    std::ostream& prettyPrint(std::ostream &os, const Vec3 &v) {
        os << '('<< v.x << ", "<< v.y<< ", " << v.z << ')';
        return os;
    }

    friend std::istream &operator>>(std::istream &is, Vec3 &v) {
        return is >> v.x >> v.y >> v.z;
    }

    friend Vec3 operator + (const Vec3 v1, const Vec3 v2) {
        return Vec3(v1.x + v2.x,v1.y + v2.y,v1.z + v2.z);
    }

    friend Vec3 operator - (const Vec3 v1, const Vec3 v2) {
        return Vec3(v1.x - v2.x,v1.y - v2.y,v1.z - v2.z);
    }

    friend Vec3 operator * (const Vec3 v1, const double factor) {
        return Vec3(v1.x * factor,v1.y * factor,v1.z* factor);
    }

    friend Vec3 operator * (const double factor,const Vec3 v1) {
        return Vec3(v1.x * factor,v1.y * factor,v1.z* factor);
    }


    friend Vec3 operator / (const Vec3 v1, double factor) {
        return Vec3(v1.x / factor,v1.y / factor,v1.z / factor);
    }

    double magnitude() const{
        return pow(x * x + y * y + z * z, .5);
    }

    double sqrMagnitude() const{
        return x * x + y * y + z * z;
    }

    void normalize() {
        auto m = magnitude();
        if(m == 0)
            return;
        x /= m;
        y /= m;
        z /= m;
    }

    Vec3 normalized() const{
        Vec3 temp(*this);
        temp.normalize();
        return temp;
    }


    double dot(const Vec3& other) const {
        return x * other.x + y * other.y +z * other.z;
    }

    Vec3 cross(const Vec3& other) const {
        Vec3 result;

        result.x = y * other.z  - z * other.y;
        result.y = z * other.x  - x * other.z;
        result.z = x * other.y  - y * other.x;

        return result;
    }

    bool magnitudeExceeds(double d){
        return (x * x + y * y + z * z)  > (d * d);
    }


    double angleBetweenInDegrees(const Vec3& other) const{
        double m1 = magnitude();
        double m2 = other.magnitude();

        return (acos(this->dot(other) / (m1 * m2)) * 180.0 ) /M_PI;
    }


    /// Rotate via Rodriguez, does not modify original
    /// \param rotateVector
    /// \param angle
    /// \return rotated vector
    Vec3 rotated(const Vec3& rotateVector, const double angle){
        auto cosFactor = cos(angle);
        auto sinFactor = sin(angle);

        return cosFactor * (*this)
               + (1-cosFactor) * (rotateVector.dot(*this)) * (rotateVector)
               + sinFactor * (rotateVector.cross(*this));
    }

    static Vec3 Right() {return Vec3(1,0,0);}
    static Vec3 Up() {return Vec3(0,1,0);}
    static Vec3 Forward() {return Vec3(0,0,1);}
    static Vec3 Zero() {return Vec3(0,0,0);}

    Matrix<4,1> toColMatrix(double wValue = 1){
        Matrix<4,1> result;
        result(0,0) = x;
        result(1,0) = y;
        result(2,0) = z;
        result(3,0) = wValue;
        return result;
    }

    friend Matrix<4,1> operator * (Matrix<4,4>& m, Vec3& v);

    Vec3 getTransformedResult(Mat4x4& transformationMatrix){
        auto t= transformationMatrix * (*this);
        Vec3 result = Vec3(t) / t(0,3);
        return result;
    }

    Vec3 static Reflect(const Vec3 &incidentVec, const Vec3 &normal)
    {
        return incidentVec - 2 * incidentVec.dot(normal) * normal;
    }
    Vec3 static refract(const Vec3 &I, const Vec3 &N, const float &ior)
    {
        float cosi = I.dot(N);
        float etai = 1, etat = ior;
        Vec3 n = N;
        if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -1 * N; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }
};
/// Set a given Vec3 as a glvertex
/// \param v given Vec3
void glVertex(const Vec3 v){
    glVertex3d(v.x, v.y, v.z);
}
/// GLTranslate along a given Vec3
/// \param v given Vec3
void glTranslate(const Vec3 v){
    glTranslated(v.x,v.y,v.z);
}



void setColsOfRow(Matrix<4,4>& mat, int rowIndex, double x, double y , double z, double w = 0) {
    mat(rowIndex,0) = x;
    mat(rowIndex,1) = y;
    mat(rowIndex,2) = z;
    mat(rowIndex,3) = w;
}

void setColsOfRow(Matrix<3,3>& mat, int rowIndex, double x, double y , double z) {
    mat(rowIndex,0) = x;
    mat(rowIndex,1) = y;
    mat(rowIndex,2) = z;
}

void swapCol(Matrix<3,3>& mat, int colIndex, double &r0,double &r1,double &r2){
    std::swap(mat(0, colIndex) , r0);
    std::swap(mat(1, colIndex) , r1);
    std::swap(mat(2, colIndex) , r2);
}

Mat4x4 TranslationMatrix(Vec3 translationAmount){
    Mat4x4 mat;
    for(int r =0; r < 4;r++){
        for(int c =0; c < 4;c++){
            mat(r,c) =(r == c)? 1 : 0;
        }
    }
    mat(0,3) = translationAmount.x;
    mat(1,3) = translationAmount.y;
    mat(2,3) = translationAmount.z;
    mat(3,3) = 1;
    return mat;
}

Mat4x4 ScalingMatrix(Vec3 scaleAmount){
    Mat4x4 mat;
    for(int r =0; r < 4;r++){
        for(int c =0; c < 4;c++){
            double scaleFactor = c == 0
                                 ? scaleAmount.x
                                 : c == 1
                                   ? scaleAmount.y
                                   : scaleAmount.z;
            mat(r,c) =(r == c)? 1 *  scaleFactor: 0;
        }
    }
    mat(3,3) = 1;
    return mat;
}


Mat4x4 IdentityMatrix(){
    Mat4x4 mat;
    for(int r =0; r < 4;r++){
        for(int c =0; c < 4;c++){
            mat(r,c) = (r == c)? 1 : 0;
        }
    }
    return mat;
}

Mat4x4 ZeroMatrix(){
    Mat4x4 mat;
    for(int r =0; r < 4;r++){
        for(int c =0; c < 4;c++){
            mat(r,c) =  0;
        }
    }
    return mat;
}

Mat4x4 RotationMatrix(Vec3 a, double angle,
                      Vec3 i = Vec3::Right(),
                      Vec3 j = Vec3::Up(),
                      Vec3 k = Vec3::Forward()){
    a.normalize();
    auto angleInRadians = deg2rad(angle);

//    offline pdf approach flips signs for some reason. Using slides version
    auto c1 = i.rotated(a, angleInRadians);
    auto c2 = j.rotated(a, angleInRadians);
    auto c3 = k.rotated(a, angleInRadians);

    Mat4x4 result;
//    cout<<result<<endl;
    setColsOfRow(result,0, c1.x, c2.x, c3.x);
    setColsOfRow(result,1, c1.y, c2.y, c3.y);
    setColsOfRow(result,2, c1.z, c2.z, c3.z);
    setColsOfRow(result,3, 0,0,0,1);

    return result;
}


void glColor(Color c){
    glColor3d(c.r,c.g,c.b);
}
void drawRect(double w , double h)
{
    w/=2;
    h/=2;
    //glColor3f(1.0,0.0,0.0);
    glBegin(GL_QUADS);{

    }glEnd();
}

void drawWireCube(double xScale, double yScale, double zScale )
{
//actually no lol
//we need half offsets
//    xScale/=2;
//    yScale/=2;
//    zScale/=2;

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_POLYGON);
    //front
    glVertex3f(xScale, yScale, zScale);
    glVertex3f(xScale, -yScale, zScale);
    glVertex3f(-xScale, -yScale, zScale);
    glVertex3f(-xScale, yScale, zScale);
    glEnd();
    glBegin(GL_POLYGON);
    //back
    glVertex3f(-xScale, yScale, -zScale);
    glVertex3f(-xScale, -yScale, -zScale);
    glVertex3f(xScale, -yScale, -zScale);
    glVertex3f(xScale, yScale, -zScale);
    glEnd();
    glBegin(GL_POLYGON);
    //right
    glVertex3f(xScale, yScale, zScale);
    glVertex3f(xScale, -yScale, zScale);
    glVertex3f(xScale, -yScale, -zScale);
    glVertex3f(xScale, yScale, -zScale);
    glEnd();
    glBegin(GL_POLYGON);
    //left
    glVertex3f(-xScale, yScale, -zScale);
    glVertex3f(-xScale, -yScale, -zScale);
    glVertex3f(-xScale, -yScale, zScale);
    glVertex3f(-xScale, yScale, zScale);
    glEnd();
    glBegin(GL_POLYGON);
    //up
    glVertex3f(-xScale, yScale, zScale);
    glVertex3f(xScale, yScale, zScale);
    glVertex3f(xScale, yScale, -zScale);
    glVertex3f(-xScale, yScale, -zScale);
    glEnd();
    glBegin(GL_POLYGON);
    //down
    glVertex3f(-xScale, -yScale, -zScale);
    glVertex3f(xScale, -yScale, -zScale);
    glVertex3f(xScale, -yScale, zScale);
    glVertex3f(-xScale, -yScale, zScale);
    glEnd();
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

/*
 * Ray calss, Contains ray origin and normalized dir
 * */
class Ray{
    Vec3 start;
    Vec3 dir;
public:
    Ray(const Vec3 &start, const Vec3 &dir) : start(start), dir(dir.normalized()) {

    }

    const Vec3 &getStart() const {
        return start;
    }

    const Vec3 &getDir() const {
        return dir;
    }

    void setStart(const Vec3 &start) {
        Ray::start = start;
    }

    void setDir(const Vec3 &dir) {
        Ray::dir = dir.normalized();
    }
    /// Get the ray endpoint for given parametric value t
    /// basically origin + dir * t
    /// \param t parametric value t
    /// \return
    Vec3 getPointAtT(double t) const{
        return start + dir * t;
    }
    //write appropriate constructor
};
Matrix<4,1> operator * (Matrix<4,4>& m, Vec3& v){
    auto mv = v.toColMatrix();
    return  m * mv;
}
/**
 * Class for desciribing a screen
 */
class Screen{
public:
    int width, height;

    Screen(const int width = 500, const int height= 500) : width(width), height(height) {}
};
static Color RefractionColors[]={
        Color::Violet(),
        Color::Indigo(),
        Color::Blue(),
        Color::Green(),
        Color::Yellow(),
        Color::Orange(),
        Color::Red()
};
/**
 * Class for Material
 */
class Material{
public:
    Color color;
    double alpha = 1;
    double ambientCoeff;
    double diffuseCoeff;
    double specularCoeff;
    double reflectionCoeff;
    int shininess;

    bool wireframe = false;
    bool refracts = false;
    //we have two mediums, we can just keep one eta for each band of color and invert when exit prism
    double entryEta[7]={
            1.3,
            1.25,
            1.2,
            1.175,
            1.150,
            1.125,
            1.1

//            2.42,
//            2.3,
//            1.25,
//            1.2,
//            1.9,
//            1.3,
//            1.15
//            1.7,
//            1.6,
//            1.5,
//            1.4,
//            1.3,
//            1.2,
//            1.1
    };


    Material(const Color &color = Color::White(),
             double ambientCoeff = .5f, double diffuseCoeff = .25f,
             double specularCoeff = .25f, double reflectionCoeff = .25f,
             int shininess = 1) : color(color), ambientCoeff(ambientCoeff),
                                         diffuseCoeff(diffuseCoeff),
                                         specularCoeff(specularCoeff),
                                         reflectionCoeff(reflectionCoeff),
                                         shininess(shininess) {}

    friend std::ostream &operator<<(std::ostream &os, const Material &material) {
        os << "color: " << material.color << " ambientCoeff: " << material.ambientCoeff << " diffuseCoeff: "
           << material.diffuseCoeff << " specularCoeff: " << material.specularCoeff << " reflectionCoeff: "
           << material.reflectionCoeff << " shininess: " << material.shininess;
        return os;
    }

    friend std::istream &operator>>(std::istream &is, Material &material) {
        is  >>  material.color
            >> material.ambientCoeff
            >> material.diffuseCoeff
            >> material.specularCoeff
            >> material.reflectionCoeff
            >> material.shininess;
        return is;
    }

};
void  glMat(const Material& mat){
    if(mat.wireframe){
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }else{
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    }

    if(mat.alpha < 1) {
        glEnable(GL_BLEND); //Enable blending.
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set blending function.
        glColor4d(mat.color.r, mat.color.g, mat.color.b, mat.alpha);
    }else{
        glDisable(GL_BLEND);
        glColor(mat.color);
    }
}
float fresnel(const Vec3 &I, const Vec3 &N, const float &ior)
{
    float cosi =  I.dot(N);
    float etai = 1, etat = ior;
    float kr;
    if (cosi > 0) { std::swap(etai, etat); }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1) {
        kr = 1;
    }
    else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        kr = (Rs * Rs + Rp * Rp) / 2;
    }
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
    return kr;
}
/**
 * Base class for drawable things
 */
class Object{
public:
    std::string name;
    Vec3 position;
    Material mat;
    /**
     * Draws the Object
     */
    virtual void draw()= 0;
    /***
     * Returns the normal vector at point for the object
     * @param v the point in world space
     * @return the normal vector at that point
     */
    virtual Vec3 normalAtPoint(Vec3 v) = 0;
    /***
     * Do insertion check with ray, return true if intersects and fills ref params
     * @param r The ray to check against
     * @param outT The distance along ray to intersection point. Valid only if intersects.
     * @return return true if intersects
     */
    virtual bool tryIntersect(const Ray &r, double &outT) = 0;

    /**
     * Returns the color at a given point
     * Expected to handle points in and out of shape bounds
     * @param point the point in world space
     * @return the color at that point
     */
    virtual Color colorAtPoint(Vec3 point){
        return mat.color;
    }
    /**
     * Prints the shape
     * @param os the output stream
     * @return
     */
    virtual std::ostream& printShape(std::ostream &os){
        os<< "position "<< position<< " color "<<mat.color;
        return os;
    }

    /// Do a recursive raytrace to this object
    /// You are expected to call this after getting the tVal from a tryIntersect call
    /// \param r rayDir
    /// \param tIntersection
    /// \param color The raytrace color result
    /// \param scene the scene we raytrace in
    /// \param reflectionLevel the recursion level/number of reflects
    virtual void
    actuallyIntersect(const Ray &r,
                      double tIntersection,
                      Color &color,
                      Scene *scene,
                      int reflectionLevel = 1) {
        Vec3 intersectionPoint = r.getPointAtT(tIntersection);
        Vec3 normalAtIntersectionPoint = normalAtPoint(intersectionPoint);

        //shadow acne problem happens due to floating point imprecision
        //the intersection point we get may be somewhat inside the object
        //if we cast ray to light from that point then it will be incorrect
        //object normal always points outside
        //so we add that by a small amount to ensure that we check for the line

        color = calcPhongLighting(r, intersectionPoint, normalAtIntersectionPoint, scene);

        if(reflectionLevel < RTX_MAX_RECURSION_LIMIT){
        {
            Color reflectionColorContribution = Color::Black();
                Ray reflectedRay(intersectionPoint  + RTX_EPSILON * normalAtIntersectionPoint, Vec3::Reflect(r.getDir(), normalAtIntersectionPoint));
                double tMin;
                int nearestObjectIndex = -1;
                int i = 0;
                for (auto &object: scene->getObjects()) {
                    double tCur;
                    if (object->tryIntersect(reflectedRay, tCur)) {
                        if (tCur > 0 && (tCur < tMin || nearestObjectIndex < 0)) {
                            tMin = tCur;
                            nearestObjectIndex = i;
                        }
                    }
                    i++;
                }

                if (nearestObjectIndex >= 0) {
                    scene->getObjects()[nearestObjectIndex]->actuallyIntersect(reflectedRay,
                                                                               tMin,
                                                                               reflectionColorContribution,
                                                                               scene,
                                                                               reflectionLevel + 1);
                    reflectionColorContribution = reflectionColorContribution * mat.reflectionCoeff;
                    reflectionColorContribution.clamp();
                }
                Color refractionColorContribution=Color::Black();
                if(mat.refracts){
                    Color refelctionBuildUp = Color::Black();
//                    std::cout<<color<<" refract start" << reflectionLevel<<" : "<<this->name<< std::endl;
//                std::cout<<" refract ?"<< std::endl;
                    //#TODO
                    float krAvg = 0;
                    double refractDot = r.getDir().dot(normalAtIntersectionPoint);
                    auto correctNormal = (refractDot > 0 ? 1: -1) *  normalAtIntersectionPoint;

                    for (int i = 0;i < 7; ++i) {
                        Color refractColor = RefractionColors[i] ;

                        double entryEta = mat.entryEta[i];

                        Vec3 refractedDir;
                        if(refractDot > 0 || refractDot < 0) {
                            double kr;
                            refractedDir = Vec3::refract(r.getDir(), normalAtIntersectionPoint, entryEta);
                            kr = fresnel(r.getDir(), normalAtIntersectionPoint, entryEta);
                            refelctionBuildUp = refelctionBuildUp + kr * reflectionColorContribution * refractColor;
                            refelctionBuildUp.clamp();

                            double kt = 1 - kr;
//                            krAvg += kr;
                            if(kr >= 1) {//no refraction total internal reflection
                                continue;
                            }
                            refractedDir.normalize();
                            Ray refractedRay(intersectionPoint +  correctNormal* RTX_EPSILON, refractedDir);
                            double tMin;
                            int nearestObjectIndex = -1;
                            int j = 0;
                            for (auto &object: scene->getObjects()) {
                                double tCur;
                                if (object->tryIntersect(refractedRay, tCur)) {
                                    if (tCur > 0 && (tCur < tMin || nearestObjectIndex < 0)) {
                                        tMin = tCur;
                                        nearestObjectIndex = j;
                                    }
                                }
                                j++;
                            }
                            if (nearestObjectIndex >= 0) {
                                Color refractionResult = Color::Black();
                                //self intersection is not an issue since faces are diff objects
                                //and we also offset by -normal to ensure that point is inside
                                scene->getObjects()[nearestObjectIndex]->actuallyIntersect(refractedRay,
                                                                                           tMin,
                                                                                           refractionResult,
                                                                                           scene,
                                                                                           reflectionLevel + 1);


                                refractionResult = refractionResult  * refractColor;
                                refractionResult = refractionResult  * kt;

                                refractionResult.clamp();
                                refractionColorContribution =  refractionColorContribution + refractionResult;
                                refractionColorContribution.clamp();
                            }
                        }
                    }
                    reflectionColorContribution = refelctionBuildUp;
                }
                reflectionColorContribution.clamp();
                color = color + reflectionColorContribution + refractionColorContribution;
                color.clamp();
            }
        }else{
//            std::cout<<" reflection level max"<< color<<" : "<<this->name<<std::endl;
        }
    }
    /// Calc phong lighting for a given raytrace
    /// \param r the ray
    /// \param point the point we want to calc the phong lighting for
    /// \param normal The normal at that point
    /// \param scene
    /// \return the phong lighting color result for the point
    Color calcPhongLighting(const Ray &r,const Vec3 &point, const Vec3 normal, Scene *scene){
        using namespace std;
        Color color = Color::Black();
        Color intersectionPointColor = colorAtPoint(point);
        color = color +  scene->ambientLighting *  mat.ambientCoeff * intersectionPointColor;
        color.clamp();
        if(!ALLOW_BACKFACE_LIGHTING && r.getDir().dot(normal) > 0) {//hitting backface
            return color;
        }
        for (auto & light:scene->getLights()) {
            if(light->hitsObjectAtPointInScene(scene, r, normal, this, point)){
                Vec3 lightDir = light->getLightDirForPoint(point);
                Vec3 viewerDir = r.getDir() * -1;//dir from point to ray origin aka cam eye is reverese of passed ray
                Vec3 reflectedLightDir = Vec3::Reflect(lightDir, normal);

                //reflectedFrom light dir . normal instead of lightDir because lightdir is incident and would produce -ve
                //they produce same angle anyways
                double diffuseReflectionFactor = reflectedLightDir.dot(normal) * mat.diffuseCoeff;
                Color diffuseReflection = (light->color) * diffuseReflectionFactor * intersectionPointColor;
                diffuseReflection.clamp();//this has -ve values when dot product < 0. These negative values should not impact ccolr from other lights
                color = color + diffuseReflection;


//                double specularReflectionFactor = pow(reflectedLightDir.dot(viewerDir), mat.shininess) * mat.specularCoeff;
                //if shininess is even then the -ve parts of dot product are removed
                double specularReflectionFactor = reflectedLightDir.dot(viewerDir);
                int specularDotSign = (specularReflectionFactor < 0 ? -1: 1);
                specularReflectionFactor = abs(specularReflectionFactor);
                //preserve sign and reapply
                specularReflectionFactor = specularDotSign
                                         * pow(specularReflectionFactor, mat.shininess)
                                         * mat.specularCoeff;
                Color specularReflection = (light->color) * specularReflectionFactor;
//                Color specularReflection = (light->color) * specularReflectionFactor * intersectionPointColor;
                specularReflection.clamp();
                color = color + specularReflection;
            }
        }
        color.clamp();
        return color;
    }
    Object(const Vec3 &position, const Material &mat) : position(position), mat(mat) {}
};

/***
 * Class for a general quadratic surface drawable object
 */
class GeneralQuadraticSurface:public Object{
public:
    /***
     * The quadratic params
     */
    double A,B,C,D,E,F,G,H,I,J;
    /**
     * The boundary cube dimension we use to clip the Surface when drawing
     */
    double clipCubeLength, clipCubeWidth, clipCubeHeight;

    GeneralQuadraticSurface(const Vec3 &position, const Material &mat, double a, double b, double c, double d, double e,
                            double f, double g, double h, double i, double j, double clipCubeLength,
                            double clipCubeWidth, double clipCubeHeight) : Object(position, mat), A(a), B(b), C(c),
                                                                           D(d), E(e), F(f), G(g), H(h), I(i), J(j),
                                                                           clipCubeLength(clipCubeLength),
                                                                           clipCubeWidth(clipCubeWidth),
                                                                           clipCubeHeight(clipCubeHeight) {}

    friend std::ostream &operator<<(std::ostream &os, const GeneralQuadraticSurface &surface) {
        os << " A: " << surface.A << " B: " << surface.B << " C: " << surface.C
           << " D: " << surface.D << " E: " << surface.E << " F: " << surface.F << " G: " << surface.G << " H: "
           << surface.H << " I: " << surface.I << " J: " << surface.J << " clipCubeLength: " << surface.clipCubeLength
           << " clipCubeWidth: " << surface.clipCubeWidth << " clipCubeHeight: " << surface.clipCubeHeight;
        return os;
    }

public:
    void draw() override {
        {
            glPushMatrix();
            glMat(mat);
            glTranslate(position);
            drawWireCube(
                         clipCubeLength == 0 ? 5: clipCubeLength,
                         clipCubeWidth == 0 ? 5: clipCubeWidth,
                         clipCubeHeight == 0 ? 5: clipCubeHeight
                         );
            glBegin(GL_LINES);
            if(clipCubeLength == 0){
                glVertex(Vec3::Right() * 50);
                glVertex(Vec3::Right() * -50);
            }
            if(clipCubeWidth == 0){
                glVertex(Vec3::Up() * 50);
                glVertex(Vec3::Up() * -50);
            }
            if(clipCubeHeight == 0){
                glVertex(Vec3::Right() * 50);
                glVertex(Vec3::Right() * -50);
            }
            glEnd();
            glPopMatrix();
        }
    }

    bool insideBoundingBox(Vec3 p){
        p = p - position;//relative coords, just need to check x y and z vals directly now
        if(clipCubeHeight != 0 && abs(p.z) > (clipCubeHeight))
            return false;
        if(clipCubeWidth != 0 && abs(p.y) > (clipCubeWidth))
            return false;
        if(clipCubeLength != 0 && abs(p.x) > (clipCubeLength))
            return false;
        return true;
    }
    bool tryIntersect(const Ray &r, double &outT) override {
        using namespace std;
        Vec3 r0 = r.getStart();
        Vec3 rd = r.getDir();
        //place ray eqn in general quadratic surface eqn , x = p.x, y = p.y, z= p.y, p = intersection point
        //this gets us an eqn of t: a * t^2 + b * t + c = 0
        double a = A * rd.x * rd.x + B * rd.y * rd.y + C * rd.z * rd.z
                 + D * rd.x * rd.y + E * rd.y * rd.z + F * rd.x * rd.z
                 ;
        double b = 2 * A * r0.x * rd.x + 2 * B * r0.y * rd.y + 2 * C * r0.z * rd.z
                 + D * (r0.x * rd.y + r0.y * rd.x)
                 + E * (r0.y * rd.z + r0.z * rd.y)
                 + F * (r0.x * rd.z + r0.z * rd.x)
                 + G * rd.x + H * rd.y + I * rd.z
                 ;
        double c = A * r0.x * r0.x + B * r0.y * r0.y + C * r0.z * r0.z
                 + D * r0.x * r0.y + E * r0.y * r0.z + F * r0.x * r0.z
                 + G * r0.x + H * r0.y + I * r0.z
                 + J
                 ;
        if (a == 0) {// t = -c / b
            if(b==0){
                if(c== 0) {//any value of t will work
                    outT = .001;//small +ve value
                    return true;
                }else{
                    return false;
                }
            }else{
                outT = - c / b;
                return true;
            }
        }
        double discriminant = b *b  - 4 * a *c;
//        cout<< "discriminant " << discriminant<<endl;
        if(discriminant < 0)
            return false; //imaginary
        discriminant = sqrt(discriminant);
//        cout<< "sqrt discriminant " << discriminant<<endl;

        double t1 = (-b + discriminant) / (2 * a);
        double t2 = (-b - discriminant) / (2 * a);
//        cout<< "t1 "<< t1<< "t2 "<< pot2<<endl;
        //take smaller +ve t
        bool t1InBoundingBox = insideBoundingBox(r.getPointAtT(t1));
        bool t2InsideBoundingBox = insideBoundingBox(r.getPointAtT(t2));
        if(t1 >= 0 && (t1 <= t2 || t2 < 0 || !t2InsideBoundingBox) && t1InBoundingBox){
            outT = t1;
            return true;
        }else {
            if(t2 >= 0 && (t2 <= t1 || t1 < 0 || !t1InBoundingBox) && t2InsideBoundingBox){
                outT = t2;
                return true;
            }
        }
        return false;
    }
    /// calc normal using partial derivatives at that point
    /// \param v point
    /// \return normal vector at point
    Vec3 normalAtPoint(Vec3 v) override {
        double ddx = 2 * A * v.x +  D * v.y + F * v.z + G;
        double ddy = 2 * B * v.y +  D * v.x + E * v.z + H;
        double ddz = 2 * C * v.z +  E * v.y + F * v.x + I;
        return Vec3(ddx, ddy, ddz).normalized();
    }


};

/**
 * Class for Pointlights
 */
class PointLight: public Light{
public:
    /// the pointlight source position
    Vec3 position;

    Vec3 getOrigin() override {
        return position;
    }

    virtual bool
    hitsObjectAtPointInScene(Scene *scene, const Ray &r, const Vec3 &normalAtPoint, Object *object, const Vec3 &point) {
        Ray rayToPoint(position, point - position);
        if(LIGHT_OPPOSITE_FACE_CHECK){
//            Vec3 normalAtPoint = object->normalAtPoint(point);
            double rDotNormal = r.getDir().dot(normalAtPoint);
            double lDotNormal = rayToPoint.getDir().dot(normalAtPoint);
            //if viewer and light on diff side
            if((rDotNormal> 0 && lDotNormal < 0) || (rDotNormal < 0 && lDotNormal > 0)){
                return false;
            }
        }


        double tVal;
        double tIntersectionOnObject;
        if(object->tryIntersect(rayToPoint, tIntersectionOnObject)){
            for (auto & other: scene->getObjects()) {
                if(other->tryIntersect(rayToPoint, tVal) && other != object && tVal < tIntersectionOnObject)
                    return false;
            }
            return true;
        }
        return false;
    }

    virtual Vec3 getLightDirForPoint(const Vec3 &point) {
        return (point - position).normalized();
    }

    PointLight(const Vec3 &position,const Color &color ) : Light(color), position(position) {}

    void debugDraw() override{
        glColor(color);
        {
            glBegin(GL_POINTS);
            glPointSize(25);
            glVertex(position);
            glEnd();
        }
        {
            glBegin(GL_LINES);
            glVertex(position + Vec3::Up() *50);
            glVertex(position + Vec3::Up() *-50);
            glVertex(position + Vec3::Right() *50);
            glVertex(position + Vec3::Right() *-50);
            glVertex(position + Vec3::Forward() *50);
            glVertex(position + Vec3::Forward() *-50);
            glEnd();
        }

    }
};

/***
 * Class for a spotlight
 */
class SpotLight: public Light{
    double cutoff_angle;
    //debugging usage
    double cutoffMaxDotProduct;
    Vec3 j;
    Vec3 k;
public:
    const Ray ray;

    Vec3 getOrigin() override {
        return ray.getStart();
    }

    virtual Vec3 getLightDirForPoint(const Vec3 &point) {
        return ray.getDir();
    }

    virtual bool
    hitsObjectAtPointInScene(Scene *scene, const Ray &r, const Vec3 &normalAtPoint, Object *object, const Vec3 &point) {

        double tIntersectionOnObject;
        Ray rayToPoint(ray.getStart(), point - ray.getStart());

        if(LIGHT_OPPOSITE_FACE_CHECK){
//            Vec3 normalAtPoint = object->normalAtPoint(point);
            double rDotNormal = r.getDir().dot(normalAtPoint);
            double lDotNormal = rayToPoint.getDir().dot(normalAtPoint);
            //if viewer and light on diff side
            if((rDotNormal> 0 && lDotNormal < 0) || (rDotNormal < 0 && lDotNormal > 0)){
                return false;
            }
        }
        if(!dotprodcheck(rayToPoint, point))
            return false;
        double tVal;

        for (auto & other: scene->getObjects()) {
            if(other->tryIntersect(rayToPoint, tVal) && other != object && tVal < tIntersectionOnObject)
                return false;
        }
        return true;

    }
    bool dotprodcheck(const Vec3 &point){
        Ray rayToPoint(ray.getStart(), point - ray.getStart());
        double dotProduct = rayToPoint.getDir().dot(ray.getDir());

        return dotProduct >= cutoffMaxDotProduct;
    }
    bool dotprodcheck(Ray &rayToPoint,const Vec3 &point){
        double dotProduct = rayToPoint.getDir().dot(ray.getDir());

        return dotProduct >= cutoffMaxDotProduct;
    }



    SpotLight(const Vec3 &position,const Color &color,const Vec3 &lightDirection, double cutoffAngle)
    :Light(color), ray(position, lightDirection)
    {
         this-> cutoff_angle = cutoffAngle;
        cutoffMaxDotProduct = cos(deg2rad(cutoffAngle));

        Vec3 crosser = Vec3::Right();
        auto i = ray.getDir();
        double dot = i.dot(crosser);
        if(dot >=  1){
            crosser = Vec3::Up();
            dot = i.dot(crosser);
            if(dot >=  1){
                crosser = Vec3::Forward();
                dot = i.dot(crosser);
            }
        }
        //ray = i
        j = i.cross(crosser);
        k = i.cross(j);
    }

    void debugDraw() override {
        glColor(color);

        {

            glBegin(GL_POINTS);
            glPointSize(10);
            glVertex(ray.getStart());
            glEnd();
        }
        {

            double dist = 200;
            double radiusAtDist  = dist * tan(deg2rad(cutoff_angle));
            glBegin(GL_LINES);
            Vec3 endpoint = ray.getStart() + ray.getDir() * dist;
            glVertex(ray.getStart());
            glVertex(endpoint);
            glVertex(ray.getStart());
            glVertex(endpoint + j * radiusAtDist);
            glVertex(ray.getStart());
            glVertex(endpoint + j * -radiusAtDist);
            glVertex(ray.getStart());
            glVertex(endpoint + k * radiusAtDist);
            glVertex(ray.getStart());
            glVertex(endpoint + k * -radiusAtDist);
            glEnd();
        }
    }
};

class Sphere: public Object{
private:
    static const int drawStacks = 32;
    static const int drawSlices = 32;

public:
    float radius;

    Sphere(const Vec3 &position, float radius, const Material &mat) : Object(position, mat), radius(radius) {}

    Vec3 normalAtPoint(Vec3 v) override{
        return (v - position).normalized();
    }

    std::ostream &printShape(std::ostream &os) override {
        os<<position << "\tradius: "<<radius;
        Object::printShape(os);
        return os;
    }

    bool tryIntersect(const Ray &r, double &outT) override {
        using namespace std;
        Vec3 r0 =  position - r.getStart() ;
        bool outsideSphere = r0.magnitudeExceeds(radius);
        double tp = r0.dot(r.getDir());//t for closes point that may be inside circle
//        cout<<position<<" r.getStart() "<<r.getStart()<<" r0 "<< r0<<" r.getDir() "<<r.getDir()<< " outsideSphere "<< outsideSphere<< " tp"<< tp <<endl;
        if(outsideSphere && tp < 0)
            return false;
        Vec3 dv = r.getPointAtT(tp) - position;
        if(dv.magnitudeExceeds(radius)){
            return false;
        }
        double t_prime = pow(radius * radius - dv.sqrMagnitude(), .5);
        double t = outsideSphere? (tp - t_prime): (tp + t_prime);
        outT = t;
//        cout<<"t push "<<t <<endl;
        return true ;
    }


    void draw() override {
        static Vec3 drawPoints[drawStacks][drawSlices];
        double h,r;
        //generate points
        for(int i=0; i <= drawStacks; i++)
        {
            h=radius*sin(((double)i/(double)drawStacks) * (pi / 2));
            r=radius*cos(((double)i/(double)drawStacks) * (pi / 2));
            for(int j=0; j <= drawSlices; j++)
            {
                drawPoints[i][j].x=r*cos(((double)j/(double)drawSlices) * 2 * pi);
                drawPoints[i][j].y=r*sin(((double)j/(double)drawSlices) * 2 * pi);
                drawPoints[i][j].z=h;
            }
        }
        {
            glPushMatrix();
            glTranslate(position);
            glMat(mat);

            //draw quads using generated points
            for (int i = 0; i < drawStacks; i++) {
                for (int j = 0; j < drawSlices; j++) {
                    glBegin(GL_QUADS);
                    {
                        //upper hemisphere
                        glVertex3f(drawPoints[i][j].x, drawPoints[i][j].y, drawPoints[i][j].z);
                        glVertex3f(drawPoints[i][j + 1].x, drawPoints[i][j + 1].y, drawPoints[i][j + 1].z);
                        glVertex3f(drawPoints[i + 1][j + 1].x, drawPoints[i + 1][j + 1].y, drawPoints[i + 1][j + 1].z);
                        glVertex3f(drawPoints[i + 1][j].x, drawPoints[i + 1][j].y, drawPoints[i + 1][j].z);
                        //lower hemisphere
                        glVertex3f(drawPoints[i][j].x, drawPoints[i][j].y, -drawPoints[i][j].z);
                        glVertex3f(drawPoints[i][j + 1].x, drawPoints[i][j + 1].y, -drawPoints[i][j + 1].z);
                        glVertex3f(drawPoints[i + 1][j + 1].x, drawPoints[i + 1][j + 1].y, -drawPoints[i + 1][j + 1].z);
                        glVertex3f(drawPoints[i + 1][j].x, drawPoints[i + 1][j].y, -drawPoints[i + 1][j].z);
                    }
                    glEnd();
                }
            }
            glPopMatrix();
        }
    }
};


class Triangle: public Object{
public:
    Vec3 points[3];

    Triangle(const Vec3 &position,const Vec3 &a, const Vec3 &b, const Vec3 &c,const Material &mat)
    : Object(position, mat) {
        points[0] = a;
        points[1] = b;
        points[2] = c;
    }



    std::ostream& printShape(std::ostream &os) override {
        os << points[0] <<std::endl
           << points[1] <<std::endl
           << points[2] <<std::endl;
        Object::printShape(os);
        return os;
    }
    void set(Vec3 a,Vec3 b, Vec3 c){
        points[0] = a;
        points[1] = b;
        points[2] = c;
    }

    friend std::ostream &operator<<(std::ostream &os, const Triangle &triangle) {
        os
                << " points: "<< std::endl
                << triangle.points[0] <<std::endl
                << triangle.points[1] <<std::endl
                << triangle.points[2] <<std::endl;
        return os;
    }

    std::ostream& operator << (std::ostream &os) {
        os << points[0] <<std::endl
           << points[1] <<std::endl
           << points[2] <<std::endl;
        return os;
    }

    /// Intersect with ray using cramers
    /// \param r
    /// \param outT
    /// \return
    bool tryIntersect(const Ray &r, double &outT) override {
        static Matrix<3,3> cramerArray;
        using namespace  std;

        Vec3 r0 = r.getStart();
        Vec3 rd = r.getDir();
        r0 = r0 - position;//local pos now

        setColsOfRow(cramerArray, 0, points[0].x - points[1].x, points[0].x - points[2].x, rd.x);
        setColsOfRow(cramerArray, 1, points[0].y - points[1].y, points[0].y - points[2].y, rd.y);
        setColsOfRow(cramerArray, 2, points[0].z - points[1].z, points[0].z - points[2].z, rd.z);
        double det = determinant(cramerArray);
//        cout<<"det"<< det<< endl;
//        cout<< cramerArray<<endl;

        if(det==0)
            return false;
        //column composed of values  on other side of = sign in cramer
        double swap1= points[0].x - r0.x;
        double swap2= points[0].y - r0.y;
        double swap3= points[0].z - r0.z;
//        cout<<"swap1 "<< swap1 << " swap2 "<< swap2 << " swap3 "<<swap3<<endl;
        swapCol(cramerArray, 0, swap1, swap2, swap3);


//        cout<< cramerArray<<endl;

        double detA = determinant(cramerArray);
        swapCol(cramerArray, 0, swap1, swap2, swap3);

        swapCol(cramerArray, 1, swap1, swap2, swap3);
        double detB = determinant(cramerArray);
        swapCol(cramerArray, 1, swap1, swap2, swap3);

        swapCol(cramerArray, 2, swap1, swap2, swap3);
        double detC = determinant(cramerArray);
        //swapCol(cramerArray, 2, swap1, swap2, swap3);//final swap back not necessary

        //point in triangle = A + beta * (B-A) + gamma * (c-A) = r0 + t * rd
        double beta = detA / det;
        double gamma = detB / det;
        double t = detC / det;

//        cout<< "beta "<< beta << " gamma "<< gamma << " t  "<< t
//        << " point "<< ((1 - beta - gamma) * points[0] + beta * points[1] + gamma * points[2])<<endl;
        if(beta > 0 && gamma > 0 && t > 0 && (beta + gamma) < 1){
            outT = t;
            return true;
        }
        return false;
    }


    Vec3 normalAtPoint(Vec3 v) override {
        return ((points[1]-points[0]).cross(points[2]-points[0]).normalized());
    }
    Vec3 getNormal()  {
        return ((points[1]-points[0]).cross(points[2]-points[0]).normalized());
    }
    bool SameSide(Vec3 p1,Vec3 p2,Vec3  a,Vec3 b) {
        auto cp1 = (b - a).cross(p1 - a);
        auto cp2 = (b - a).cross(p2-a);
        if (cp1.dot(cp2) >= 0)
            return true;
        return false;
    }

    bool contains(Vec3 p){
        static double threshold = .000001;

        p = p - position;
//        if (SameSide(p,points[0], points[1],points[2]) &&
//            SameSide(p,points[1], points[0],points[2]) &&
//            SameSide(p,points[2], points[0],points[1]))
//            return true;
//        return false;

        auto det = .5 * ((points[0] - points[1]).cross(points[0] - points[2])).magnitude();
        auto detA = .5 * ((p - points[1]).cross(p - points[2])).magnitude();
        auto detB = .5 * ((p - points[0]).cross(p - points[2])).magnitude();
        auto detC = .5 * ((p - points[0]).cross(p - points[1])).magnitude();
        if(det ==0)
            return false;
        auto alpha = detA / det;
        auto beta = detB / det;
        auto gamma = detC / det;

        return alpha > 0 && beta > 0 && gamma > 0 && alpha < 1 && abs(alpha + beta + gamma - 1)<threshold;
    }


    ~Triangle() {}

    void draw() override {
        {
            glPushMatrix();
            glTranslated(position.x, position.y, position.z);
            glMat(mat);
            glBegin(GL_TRIANGLES);
            {
                glVertex3f(points[0].x, points[0].y, points[0].z);
                glVertex3f(points[1].x, points[1].y, points[1].z);
                glVertex3f(points[2].x, points[2].y, points[2].z);
            }
            glEnd();
            glPopMatrix();
        }
    };
};


class Quad: public Object{
public:
    Vec3 normal;
    ///dim of quad along own x axis
    double width;
    ///dim of quad along own y axis
    double height;
    ///The x axis of the quad
    Vec3 right;
    ///The y axis of the quad
    Vec3 up;

    Vec3 normalAtPoint(Vec3 v) override {
        return normal;
    }



    Quad(const Vec3 &position, const Vec3 &normal, const Vec3 &right, const Vec3 &up, double width, double height, const Material mat)
            : normal(normal), right(right), up(up), width(width), height(height),  Object(position, mat)  {}

    void draw(){
        {
            glBegin(GL_LINES);
            glColor(mat.color.inverted() + .1f);
            auto origin = position + width *.5f * right + height *.5f * up;
            glVertex(origin);
            glVertex(origin + (normal * 25));
            glEnd();
        }
        glMat(mat);
        glBegin(GL_POINTS);
        glVertex(position);
        glEnd();


        glBegin(GL_QUADS);
        glVertex(position);
        glVertex(position + right * width);
        glVertex(position + right * width + up * height);
        glVertex(position + up * height);
        glEnd();

        glColor(mat.color.inverted());
        glBegin(GL_LINES);
        glVertex(position);
        glVertex(position + normal * 10);
        glEnd();

    }
    bool tryIntersect(const Ray &r, double &outT){
        using namespace std;


        double nDotRayDir = normal.dot(r.getDir());

        if(nDotRayDir == 0) {
            return false;
        }

        double t =(normal.dot(position) - normal.dot(r.getStart())) / nDotRayDir;
        Vec3 point = r.getPointAtT(t);

        if(!containsPointOnPlane(point))
            return false;
        outT = t;
        return true;
    }

    bool containsPointOnPlane(Vec3 point){
        point = point - position;
        double rightDot = right.dot(point);
        double upDot = up.dot(point);

        if( rightDot >= 0 && rightDot <= width && upDot >= 0 && upDot <= height)
            return true;
        return false;
    }

    friend std::ostream& operator << (std::ostream &os, const Quad &quad) {
        os << "position " << quad.position <<std::endl
        << "normal " << quad.normal <<std::endl
        << "right " << quad.right <<std::endl
        << "up " << quad.up <<std::endl
        << "width " << quad.width <<std::endl
        << "height " << quad.height <<std::endl
        ;
        return os;
    }

};
/// get a vec normal to passed vec using 3 axis dot hack no guarantee about which normal
/// \return normal
Vec3 getNormalTo(const Vec3& i){
    Vec3 crosser = Vec3::Right();
    double dot = i.dot(crosser);
    if(dot >=  1){
        crosser = Vec3::Up();
        dot = i.dot(crosser);
        if(dot >=  1){
            crosser = Vec3::Forward();
            dot = i.dot(crosser);
        }
    }
    //ray = i
    return i.cross(crosser);
}

class Prism :public Object{
private:
    Triangle top;
    Triangle bottom;
    Quad ab, bc, ca;
    Vec3 up;
public:


    Prism(const Vec3 position, const Material &mat,
          const Triangle top, const Triangle bottom,
          const Quad ab,
          const Quad bc, const Quad ca, const Vec3 up) : Object(position, mat), top(top), bottom(bottom), ab(ab),
                                                            bc(bc), ca(ca), up(up) {}

    void draw() override {
        top.draw();
        bottom.draw();
//
        ab.draw();
        bc.draw();
        ca.draw();
    }

    Vec3 normalAtPoint(Vec3 v) override {
        if(top.contains(v))
            return  top.getNormal();
        if(bottom.contains(v))
            return  bottom.getNormal();
        if(ab.containsPointOnPlane(v))
            return  ab.normal;
        if(bc.containsPointOnPlane(v))
            return  bc.normal;
        if(ca.containsPointOnPlane(v))
            return  ca.normal;

        return Vec3();
    }

    bool tryIntersect(const Ray &r, double &outT) override {
        bool any  = false;
        double t;
        if(top.tryIntersect(r, t)) {
            if(t > 0 && (t < outT  || !any)){
                any = true;
                outT =t;
            }
        }
        if(bottom.tryIntersect(r, t)) {
            if(t > 0 && (t < outT  || !any)){
                any = true;
                outT =t;
            }
        }
        if(ab.tryIntersect(r, t)) {
            if(t > 0 && (t < outT  || !any)){
                any = true;
                outT =t;
            }
        }
        if(bc.tryIntersect(r, t)) {
            if(t > 0 && (t < outT  || !any)){
                any = true;
                outT =t;
            }
        }
        if(ca.tryIntersect(r, t)) {
            if((t < outT && t > 0 )|| !any){
                any = true;
                outT =t;
            }
        }

        return any;
    }
    Color colorAtPoint(Vec3 v) override{
        if(ab.containsPointOnPlane(v))
            return  ab.mat.color;
        if(bc.containsPointOnPlane(v))
            return  bc.mat.color;
        if(ca.containsPointOnPlane(v))
            return  ca.mat.color;
        if(top.contains(v))
            return  top.colorAtPoint(v);
        if(bottom.contains(v))
            return  bottom.colorAtPoint(v);
        return Color(1,0,1);
    }
    //because c++ is stupid and doesn't allwo creating custom classes inside the constructor
    static Prism * create(const Vec3 position,Vec3 a, Vec3 b, Vec3 c, const Material &mat) {
        Vec3 centroid =(a +b +c) / 3;
        auto up = (centroid - position);
        double halfHeight = up.magnitude();
        double height = halfHeight - 2;
        up.normalize();
        //reverse order to ensure normal
        auto top = Triangle(Vec3::Zero(),c ,b ,a,  Material(Color(0,1,1)));
        auto bottom = Triangle(Vec3::Zero(), a - up * height,b- up * height ,c- up * height,  Material(Color(1,1,0)));

        Vec3 abEdge= b -a;
        double abDistFromCentroid = ((a + b) / 2 - centroid).magnitude();
        Vec3 abNormal = up.cross(abEdge ).normalized();
        auto ab = Quad(bottom.points[0],abNormal,abEdge.normalized(),up,abEdge.magnitude(), height,Material(Color::Red()));

        Vec3 bcEdge= c -b;
        double bcDistFromCentroid = ((b + c) / 2 - centroid).magnitude();
        Vec3 bcNormal = up.cross(bcEdge).normalized();
        auto bc = Quad(bottom.points[1], bcNormal,bcEdge.normalized(), up,bcEdge.magnitude(), height,Material(Color::Blue()));

        Vec3 caEdge= a -c;
        double caDistFromCentroid = ((c + a) / 2 - centroid).magnitude();
        Vec3 caNormal = up.cross(caEdge).normalized();

        auto ca = Quad(bottom.points[2], caNormal,caEdge.normalized(), up,caEdge.magnitude(), height,Material(Color::Green()));

        return new Prism(position, mat, top, bottom, ab , bc , ca, up);
    }



};
class Floor: public Object{
    double width = 1000;
    double height = 1000;
    double checkerWidth = 20;
    double checkerHeight = 20;

public:
    Floor(const Vec3 &position = Vec3::Zero()) : Object(position, Material(Color::White())) {}
    Vec3 getNormal(){return Vec3::Forward();}

    Vec3 normalAtPoint(Vec3 v) override {
        return getNormal();
    }

    bool tryIntersect(const Ray &r, double &outT) override {
        using namespace std;

        Vec3 normal = getNormal();
        double nDotRayDir = normal.dot(r.getDir());
//        cout<<nDotRayDir<<endl;

        if(nDotRayDir == 0) {
            return false;
        }

        double t =(normal.dot(position) - normal.dot(r.getStart())) / nDotRayDir;
        Vec3 point = r.getPointAtT(t);
        point = point - (position - Vec3(width/2, height/2));
        if(point.x < 0 || point.x > width || point.y < 0 || point.y > height)
            return false;
        outT = t;
        return true;
    }



    Color colorAtPoint(Vec3 vec) override{
        auto localPos = vec - (position - Vec3(width/2, height/2));
        return colorAtLocalPoint(localPos);
    }
    Color colorAtLocalPoint(Vec3 localPos){
        int x = localPos.x /  checkerWidth;
        int y = localPos.y / checkerHeight;

//        return ((x + y) % 2 == 0 ? Color::White() : Color::White()*.5f)
        return ((x + y) % 2 == 0 ? Color::White() : Color::Black())
        ;
//               * ((y % 2 + 1 ) * 1.0f / 2.0f);
    }



    void draw() override {
        {
            glPushMatrix();
            glMat(mat);
            glColor(Color::Blue());
            glBegin(GL_POINTS);
            glVertex(position);
            glEnd();
            glBegin(GL_LINES);
            glVertex(Vec3::Zero());
            glVertex(position+getNormal()*80);
            glEnd();
            glTranslate(position +Vec3(-width/2, -height/2));



            glBegin(GL_QUADS);
            {
                for (double x = 0; x < width; x+= checkerWidth){
                    for (double y = 0 ; y < height; y+= checkerHeight){
                        auto c = colorAtLocalPoint(Vec3(x,y));
                        glColor(c) ;
                        glVertex(Vec3(x,y));
                        glVertex(Vec3(x , y + checkerHeight));
                        glVertex(Vec3(x + checkerWidth,y + checkerHeight));
                        glVertex(Vec3(x + checkerWidth,y));
                    }
                }


            }
            glEnd();
            glPopMatrix();
        }
    }


};



class Camera{
public:
    Vec3 position;
    Vec3 forward;
    Vec3 right;
    Vec3 up;
    double viewAngle= 80;

    Camera(): forward(Vec3(0,0,-1)), right(Vec3(1,0,0)), up(Vec3(0,1,0)), position(Vec3::Zero())
    {

    }

    Camera(const Vec3 &position, const Vec3 &lookAtPoint, const Vec3 &up) : position(position), forward(forward) {
        this->position = position;
        this->forward = lookAtPoint - position;
        this->forward.normalize();
        this->right  =  this->forward.cross(up);
        this->right.normalize();
        this->up =  this->right .cross( this->forward);

    }


    void set(const Vec3 &pos,
             const Vec3 forward,
             const Vec3 right,
             const Vec3 up){
        this->position = pos;
        this->forward = forward;
        this->right = right;
        this->up = up;
    }



    void applyTransformation(Mat4x4 transformationMatrix) {
        position.getTransformedResult(transformationMatrix);
    }

    std::ostream & print(std::ostream &ost) {
        ost << "position:" << position << std::endl;
        ost<< "forward:"<< forward<<std::endl;
        ost<< "right:"<< right<<std::endl;
        ost<< "up:"<< up<<std::endl;
        return ost;
    }
};
void glCamUpdate(Camera * cam){
    gluLookAt(cam->position.x,cam->position.y,cam->position.z,
              cam->position.x + cam->forward.x, cam->position.y + cam->forward.y, cam->position.z + cam->forward.z,
              cam->up.x, cam->up.y, cam->up.z);

}



const std::vector<Object *>& Scene::getObjects() const {
    return objects;
}

const std::vector<Light *>& Scene::getLights() const {
    return lights;
}

void Scene::addObject(Object* object){
    objects.push_back(object);
}

void Scene::addLight(Light* light){
    lights.push_back(light);
}
class RefractResult{
public:
    Vec3 intersectionPoint;
    Vec3 refractedDir;
    Vec3 reflectedDir;
    Light * finalLight = nullptr;
    Vec3 correctNormal = Vec3::Zero();
    Color incomingColor = Color::White();
    RefractResult * refractedFrom = nullptr;
    RefractResult * reflectedFrom = nullptr;
    Object * o;

    RefractResult(const Vec3 &intersectionPoint, Object *o) : intersectionPoint(intersectionPoint), o(o) {}
};
RefractResult * traceRay(Object * o,const Ray &r,double tIntersection,Color &color,Scene *scene, std::vector<RefractResult*> &points,int reflectionLevel = 1){
    using namespace std;
    auto intersectionPoint = r.getPointAtT(tIntersection);
    Vec3 normalAtIntersectionPoint = o->normalAtPoint(intersectionPoint);
    RefractResult * rr = new RefractResult(intersectionPoint, o);
    //shadow acne problem happens due to floating point imprecision
    //the intersection point we get may be somewhat inside the object
    //if we cast ray to light from that point then it will be incorrect
    //object normal always points outside
    //so we add that by a small amount to ensure the we check for the line

//        color = calcPhongLighting(r, intersectionPoint, normalAtIntersectionPoint, scene);

    if(reflectionLevel < RTX_MAX_RECURSION_LIMIT){
        {
            Color reflectionColorContribution = Color::Black();
//            if(!mat.refracts) {

            Ray reflectedRay(intersectionPoint  + RTX_EPSILON * normalAtIntersectionPoint, Vec3::Reflect(r.getDir(), normalAtIntersectionPoint));
            rr->reflectedDir =reflectedRay.getDir();
//            cout << normalAtIntersectionPoint << " <- "<< r.getDir()  << " -> "<<rr->reflectedDir <<endl;
            double tMin;
            int nearestObjectIndex = -1;
            int i = 0;
            for (auto &object: scene->getObjects()) {
                double tCur;
                if (object->tryIntersect(reflectedRay, tCur)) {
                    if (tCur > 0 && (tCur < tMin || nearestObjectIndex < 0)) {
                        tMin = tCur;
                        nearestObjectIndex = i;
                    }
                }
                i++;
            }
            for (auto  &light:scene->getLights()) {
                if(light->hitsObjectAtPointInScene(scene,reflectedRay, normalAtIntersectionPoint,o, intersectionPoint)){
                    rr->finalLight = light;
                    break;
                }
            }

            if (nearestObjectIndex >= 0) {
//                scene->getObjects()[nearestObjectIndex]->actuallyIntersect(reflectedRay,
//                                                                           tMin,
//                                                                           reflectionColorContribution,
//                                                                           scene,
//                                                                           reflectionLevel + 1);
                auto rfl = traceRay(scene->getObjects()[nearestObjectIndex], reflectedRay,
                         tMin,
                         reflectionColorContribution,
                         scene,
                         points,reflectionLevel + 1);
                rfl->reflectedFrom = rr;

//                reflectionColorContribution = reflectionColorContribution * o->mat.reflectionCoeff;
//                reflectionColorContribution.clamp();
            }
//            }
            Color refractionColorContribution=Color::Black();
            if(o->mat.refracts){
                double refractDot = r.getDir().dot(normalAtIntersectionPoint);
                auto correctNormal = (refractDot > 0 ? 1: -1) *  normalAtIntersectionPoint;
                rr->correctNormal = correctNormal;
                cout<<o->name << normalAtIntersectionPoint << " normal and corrected :"<< rr->correctNormal << " ray "<< r.getDir()<<" dot "<< refractDot<<endl;
//                    std::cout<<color<<" refract start" << reflectionLevel<<" : "<<this->name<< std::endl;
//                std::cout<<" refract ?"<< std::endl;
                //#TODO
                float krAvg = 0;
//                int  i = ccc;

                for (int i = 0;i < 7; ++i) {
                    Color refractColor = RefractionColors[i] ;
                    double entryEta = o->mat.entryEta[i];
//                    cout<<entryEta<<endl;
                    Vec3 refractedDir;
                    if(refractDot > 0 || refractDot < 0) {
                        double kr;

                        refractedDir = Vec3::refract(r.getDir(), normalAtIntersectionPoint, entryEta);
                        cout<<" dot  ray dir should be > 0 : "<<refractedDir.dot(correctNormal)<<" kr "<< kr<<endl;
                        kr = fresnel(r.getDir(), normalAtIntersectionPoint, entryEta);

                        //#TODO check if this should be done
                        krAvg += kr;
                        if(kr >= 1) {
                            cout<<"skip"<<endl;
                            continue;
                        }

                        refractedDir.normalize();

//                        Ray refractedRay(intersectionPoint, refractedDir);

                        Ray refractedRay(intersectionPoint +  correctNormal* RTX_EPSILON, refractedDir);
                        rr->refractedDir = refractedDir;
                        //                        Ray refractedRay(intersectionPoint + refractedDir * .01, refractedDir);
                        double tMin;
                        int nearestObjectIndex = -1;
                        int j = 0;
                        for (auto &object: scene->getObjects()) {
                            double tCur;
                            if (object->tryIntersect(refractedRay, tCur)) {
                                if (tCur > 0 && (tCur < tMin || nearestObjectIndex < 0)) {
                                    tMin = tCur;
                                    nearestObjectIndex = j;
                                }
                            }
                            j++;
                        }
                        if (nearestObjectIndex >= 0) {
                            Color refractionResult = Color::Black();
                            //self intersection is not an issue since faces are diff objects
                            //and we also offset by -normal to ensure that point is inside
                            auto * rrr = traceRay(scene->getObjects()[nearestObjectIndex], refractedRay,
                                                                                       tMin,
                                                                                       refractionResult,
                                                                                       scene,
                                     points,reflectionLevel + 1);
                            rrr->refractedFrom = rr;
                            rrr->incomingColor = refractColor;

//                                color =  kr * color + (1-kr) * refractionColorContribution * incomingColor;
                            //wait kr isdifferent every time #TODO may be we shoud factor tat in reflection color COntribution
                            // or dot product with color? idk
//                                refractionResult = (1-kr) * refractionResult  * incomingColor;
                            refractionResult = refractionResult  * refractColor;

                            double kt = 1 - kr;


                            refractionResult = refractionResult  * kt;
//                                std::cout
//                                <<refractionResult<<" refract base color "
//                                <<incomingColor<<" refract with "
//                                <<refractionColorContribution <<" at lv "<< reflectionLevel<<" kr "<< kr
//                                <<" : "<<scene->getObjects()[nearestObjectIndex]->name<< std::endl;

                            refractionResult.clamp();
                            refractionColorContribution =  refractionColorContribution + refractionResult;
                            refractionColorContribution.clamp();

//                                std::cout<<refractionColorContribution<<  "  and "<< refractionResult <<std::endl;
                        }
                    }
                }
//                    std::cout<<refractionColorContribution<<" refract finish" << reflectionLevel<<" : "<<this->name<< std::endl;
                reflectionColorContribution = reflectionColorContribution * (krAvg / 7);
            }

            color = color + reflectionColorContribution + refractionColorContribution;
            color.clamp();
        }
    }
    points.push_back(rr);
    return rr;
}
#endif //DEF_1705077_CLASSES_H
