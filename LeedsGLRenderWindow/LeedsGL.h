#ifndef LEEDSGL_H
#define LEEDSGL_H

#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <memory>
#include <cstddef>
#include "Matrix4.h"
#include "RGBAValue.h"
#include "RGBAValueF.h"
#include "RGBAImage.h"
#include "ThreadPool.h"

namespace LeedsGLUtils {
    Matrix4 calculateViewportMatrix(float cx, float cy, float width, float height);

    Matrix4 calculateProjectionFrustum(float left, float right, float bottom, float top, float near, float far);

    Matrix4 calculateProjectionOrtho(float left, float right, float bottom, float top, float near, float far);

    float distancePointLine(Cartesian3 r, Cartesian3 n, Cartesian3 p);
}

// class with vertex attributes
struct InputVertex {
    // position
    Homogeneous4 position;
    // normal
    Cartesian3 normal;
    // color
    RGBAValueF color;
    // texture coords
    Cartesian3 tex_coord;
};

struct TransformedVertex {
    Homogeneous4 position;
    // coords in vsc
    Homogeneous4 v_vcs;
    // normal
    Cartesian3 normal;
    // color
    RGBAValueF color;
    // texture coords
    Cartesian3 tex_coord;
    // transformed w, used for interpolation
    float w;

};

struct Primitive {
    std::vector<TransformedVertex> transformedVertices;
};

struct Fragment {
    int row;
    int col;
    int width;
    int height;
    std::vector<RGBAValueF> colors;
    std::vector<float> depths;

    void resize() {
        int size = width * height;
        colors.clear();
        colors.resize(size, {0, 0, 0, 0});
        depths.clear();
        depths.resize(size);
    }
};


class LeedsGL {
public:
    LeedsGL();

    ~LeedsGL();

    //RENDERING PARAMETERS:
    void setUniform(const std::string &name, const bool value);

    void setUniform(const std::string &name, const Matrix4 &mat);

    void setUniform(const std::string &name, const RGBAValueF &col);

    void setUniform(const std::string &name, const Homogeneous4 &pos);

    void setUniform(const std::string &name, const float val);

    //PIPELINE CONTROL:
    void clearColor(const RGBAValueF &col);

    void clear(std::byte mask);

    void enable(const std::byte function);

    void disable(const std::byte function);

    void texImage2D(RGBAImage const *textureImage);

    void resizeBuffers(unsigned const int width, unsigned const int height);

    void lineWidth(const float width);

    void pointSize(const float size);

    //MAIN PIPELINE IMPLEMENTATION
    void drawArrays(const std::vector<Homogeneous4> &vertices,
                    const std::vector<Homogeneous4> &normals,
                    const std::vector<Cartesian3> &textureCoordinates,
                    const std::vector<RGBAValueF> &colors, std::byte mode);

    void inputAssembly(const std::vector<Homogeneous4> &vertices,
                       const std::vector<Homogeneous4> &normals,
                       const std::vector<Cartesian3> &textureCoordinates,
                       const std::vector<RGBAValueF> &colors,
                       std::vector<InputVertex> &result);

    void transformVertices(std::vector<InputVertex> &vertices,
                           std::vector<TransformedVertex> &result);

    void primitiveAssembly(std::vector<TransformedVertex> &vertices,
                           std::byte mode,
                           std::vector<Primitive> &result);

    void clipAndCull(std::vector<Primitive> &primitives,
                     std::byte mode,
                     std::vector<Primitive> &result);

    void rasterisePrimitives(std::vector<Primitive> &primitives,
                             std::byte mode,
                             std::vector<Fragment> &result);

    void rasterisePoint(int index, const Primitive &point,
                        std::vector<Fragment> &output);

    void rasteriseLine(int index, const Primitive &line,
                       std::vector<Fragment> &output);

    void rasteriseTriangle(int index, const Primitive &triangle,
                           std::vector<Fragment> &output);

    void processFragments(std::vector<Fragment> &fragments);

    //SHADING.
    RGBAValueF CalculateLighting(const Homogeneous4 &n_vcs,
                                 const Homogeneous4 &v_vcs,
                                 const RGBAValueF &em,
                                 const RGBAValueF &am,
                                 const RGBAValueF &diff,
                                 const RGBAValueF &spec,
                                 float shin);

    // texture color
    RGBAValueF textureSampler(const Cartesian3 &uv);

    // generate a color depend on the state
    RGBAValueF calculateColor(const RGBAValueF &color,
                              const Cartesian3 &uv,
                              const Homogeneous4 &n_vcs,
                              const Homogeneous4 &v_vcs);

    //BUFFERS
    RGBAImage frameBuffer;
    RGBAImage swapBuffer;
    std::vector<float> depthBuffer;

    //Masks
    static const std::byte UNKNOWN_MASK{0};
    static const std::byte COLORMASK{1};
    static const std::byte DEPTHMASK{2};
    //Function constants
    static const std::byte DEPTHTEST{3};
    static const std::byte PERSPECTIVE{4};
    //Drawing modes

    static const std::byte POINTS{5};
    static const std::byte LINES{6};
    static const std::byte TRIANGLES{7};

private:

//uniform variables
    RGBAImage const *enabledTexture;
    bool texturingEnabled;
    bool textureModulationEnabled;
    bool UVColourDebug;
    bool lightingEnabled;

    Matrix4 lightMatrix;
    Matrix4 viewPortMatrix;
    Matrix4 projectionMatrix;
    Matrix4 modelviewMatrix;

    Homogeneous4 lightPosition;
    RGBAValueF lightColour;

    RGBAValueF emissiveMaterial;
    RGBAValueF ambientMaterial;
    RGBAValueF diffuseMaterial;
    RGBAValueF specularMaterial;
    float shininessMaterial;

//global states
    float rasterizedLineWidth;
    float rasterizedPointSize;
    RGBAValueF bufferClearColor;
    bool depthTestEnabled;
    bool perspective;

//Queues
    std::vector<InputVertex> inputQueue;
    std::vector<TransformedVertex> transformedQueue;
    std::vector<Primitive> primitivesQueue;
    std::vector<Primitive> clippedPrimitivesQueue;
    std::vector<Fragment> fragmentQueue;

    // texture
    const RGBAImage *texture;
    // multiple threads
    TP::ThreadPool *pool;
};

#endif // LEEDSGL_H
