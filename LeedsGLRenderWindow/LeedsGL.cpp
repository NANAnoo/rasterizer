#include "LeedsGL.h"
#include <math.h>

using namespace TP;
using std::byte;
using std::string;

LeedsGL::LeedsGL() {
    //TODO: Initialize all variables that should be initialized!
    this->pool = new ThreadPool(std::thread::hardware_concurrency());
}

LeedsGL::~LeedsGL() {
    //TODO: Free any resources you allocate.
    delete this->pool;
}

void LeedsGL::clear(byte mask) {
    //TODO: Clear the buffer requested using the provided mask.
    //Look at the usage of this function on LeedsGLRenderWidget to
    //know what to expect of the "mask" parameter. Each bit should
    //tell which buffer to clear.

    if ((mask & COLORMASK) != UNKNOWN_MASK) {
        // clear swapBuffer
        for (int i = 0; i < swapBuffer.height; i ++) {
            for (int j = 0; j < swapBuffer.width; j ++) {
                swapBuffer[i][j].red = std::byte(bufferClearColor.red * 255.f);
                swapBuffer[i][j].green = std::byte(bufferClearColor.green * 255.f);
                swapBuffer[i][j].blue = std::byte(bufferClearColor.blue * 255.f);
                swapBuffer[i][j].alpha = std::byte(bufferClearColor.alpha * 255.f);
            }
        }
    }
    if ((mask & DEPTHMASK) != UNKNOWN_MASK) {
        // TODO: clear depth buffer
    }
}

void LeedsGL::setUniform(const string &name, const bool value) {
    //TODO: set an uniform bool with given name value.
    //uniform variables should decide how your shading happens
}

void LeedsGL::setUniform(const string &name, const Matrix4 &mat) {
    //TODO: set an uniform matrix with given name.
    //uniform variables should decide how your shading happens
}

void LeedsGL::setUniform(const string &name, const RGBAValueF &col) {
    //TODO: set an uniform colour with given name.
    //uniform variables should decide how your shading happens
}

void LeedsGL::setUniform(const std::string &name, const float val) {
    //TODO: set an uniform float with given name.
    //uniform variables should decide how your shading happens
}

void LeedsGL::setUniform(const std::string &name, const Homogeneous4 &pos) {
    //TODO: set an uniform position with given name.
    //uniform variables should decide how your shading happens
}

void LeedsGL::clearColor(const RGBAValueF &col) {
    //TODO: set a specific pipeline value, the color used to clear the buffer.
    bufferClearColor = col;
}

void LeedsGL::resizeBuffers(unsigned const int width, unsigned const int height) {
#ifdef __APPLE__
    // pixel scale on macos
    frameBuffer.Resize(2 * width, 2 * height);
#elif
    frameBuffer.Resize(width, height);
#endif
    swapBuffer.Resize(width, height);
    //TODO: Implement how you resize your depth buffer
}

Matrix4 LeedsGLUtils::calculateViewportMatrix(float cx, float cy, float width, float height) {
    //TODO: return a matrix that given parameters cx,cy (center of your viewport) width height
    //performs the viewport transformation. NDCS -> DCS.
    return Matrix4();
}

Matrix4
LeedsGLUtils::calculateProjectionOrtho(float left, float right, float bottom, float top, float near, float far) {
    //TODO: return a Ortographic projection matrix, with the parameters above.
    //right or left handedness may have effects on other parts of your code,
    //such as shading, clipping, culling, etc.
    return Matrix4();
}

Matrix4
LeedsGLUtils::calculateProjectionFrustum(float left, float right, float bottom, float top, float near, float far) {
    //TODO: return a Perspective projection matrix, with the parameters above.
    //right or left handedness may have effects on other parts of your code,
    //such as shading, clipping, culling, etc.
    return Matrix4();
}


void LeedsGL::texImage2D(RGBAImage const *textureImage) {
    //TODO: set in your pipeline which texture should be used to render.
    //Parameter is a pointer to the texture, be aware of how it is stored, and ownership of the resources.
}

void LeedsGL::enable(const std::byte function) {
    //TODO: enables a pipeline function described by the byte parameter.

}

void LeedsGL::disable(const std::byte function) {
    //TODO: disables a pipeline function described by the byte parameter.

}

void LeedsGL::lineWidth(const float width) {
    //TODO: Set a variable that describes what is the width in pixels of the line to be rasterized.
}

void LeedsGL::pointSize(const float size) {
    //TODO: Set a variable that describes what is the size in pixels of the point to be rasterized.
}

void LeedsGL::drawArrays(const std::vector<Homogeneous4> &vertices, const std::vector<Homogeneous4> &normals,
                         const std::vector<Cartesian3> &textureCoordinates, const std::vector<RGBAValueF> &colors,
                         std::byte mode) {
    //Calls the whole pipeline, step by step.
    inputAssembly(vertices, normals, textureCoordinates, colors, inputQueue);
    transformVertices(inputQueue, transformedQueue);
    primitiveAssembly(transformedQueue, mode, primitivesQueue);
    clipAndCull(primitivesQueue, mode, clippedPrimitivesQueue);
    rasterisePrimitives(clippedPrimitivesQueue, mode, fragmentQueue);
    processFragments(fragmentQueue);
    // swap buffer
    for (int i = 0; i < swapBuffer.height; i ++) {
        for (int j = 0; j < swapBuffer.width; j ++) {
#ifdef __APPLE__
            // pixel resize on mac
            frameBuffer[2 * i][2 * j] = swapBuffer[i][j];
            frameBuffer[2 * i + 1][2 * j] = swapBuffer[i][j];
            frameBuffer[2 * i][2 * j + 1] = swapBuffer[i][j];
            frameBuffer[2 * i + 1][2 * j + 1] = swapBuffer[i][j];
#else
            frameBuffer[i][j] = swapBuffer[i][j];
#endif
        }
    }
}

void LeedsGL::inputAssembly(const std::vector<Homogeneous4> &vertices, const std::vector<Homogeneous4> &normals,
                            const std::vector<Cartesian3> &textureCoordinates, const std::vector<RGBAValueF> &colors,
                            std::vector<InputVertex> &result) {
    //TODO: Check how the input is passed to drawArrays.
    //This function should combine this disjoint information into a series of InputVertex to be processed
    //by the next step.
    // build up  input vertices
    result.resize(vertices.size());
    TaskGroup tasks;
    for (unsigned int i = 0; i < vertices.size(); i++) {
        tasks.emplace_back([i, &vertices, &normals,
                                   &textureCoordinates,
                                   &colors, &result]() {
            result[i].position = vec4(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w);
            if (i < normals.size())
                result[i].normal = vec3(normals[i].x, normals[i].y, normals[i].z);
            if (i < textureCoordinates.size())
                result[i].tex_coord = vec2(textureCoordinates[i].x, textureCoordinates[i].y);
            if (i < colors.size())
                result[i].color = vec4(colors[i].red, colors[i].green, colors[i].blue, colors[i].alpha);
        });
    }
    this->pool->syncGroup(tasks, (int)(tasks.size()) / 100);
}

void LeedsGL::transformVertices(std::vector<InputVertex> &vertices, std::vector<TransformedVertex> &result) {
    //TODO: Transform the input vertices using the matrices set in LeedsGLRenderWidget.
    //Also pass all the necessary information to the next steps of the pipeline.
    //You should check the slides to decide which is the appropriate coordinate system to transform them to.
}

void
LeedsGL::primitiveAssembly(std::vector<TransformedVertex> &vertices, std::byte mode, std::vector<Primitive> &result) {
    //TODO: Assemble the vertices into a primitive according to the selected mode.
}

void LeedsGL::clipAndCull(std::vector<Primitive> &primitives, std::byte mode, std::vector<Primitive> &result) {
    //TODO: Implement clipping and culling. Should have a different behavior for each type of primitive.
    //Pay attention to what type of projection you are using, as your clipping planes will be different.
    //If you choose to skip this step as it is one of your last tasks, just return all the same primitives.
}

void LeedsGL::rasterisePrimitives(std::vector<Primitive> &primitives, std::byte mode, std::vector<Fragment> &results) {
    //TODO: Generate a list of fragments according to what mode is chosen. Should call the "rasterise X" functions
}

void LeedsGL::rasterisePoint(const Primitive &point, std::vector<Fragment> &output) {
    //TODO: Rasterise a point, according to the pointSize.
}

void LeedsGL::rasteriseLine(const Primitive &line, std::vector<Fragment> &output) {
    //TODO: Rasterise a line, according to the linewidth
}


float LeedsGLUtils::distancePointLine(Cartesian3 r, Cartesian3 n, Cartesian3 p) {
    //assumes n is normalized
    return n.dot(r) - n.dot(p);
}

void LeedsGL::rasteriseTriangle(const Primitive &triangle, std::vector<Fragment> &output) {
    TransformedVertex vertex0 = triangle.transformedVertices[0];
    TransformedVertex vertex1 = triangle.transformedVertices[1];
    TransformedVertex vertex2 = triangle.transformedVertices[2];

    // compute a bounding box that starts inverted to frame size
    // clipping will happen in the raster loop proper
    float minX = swapBuffer.width, maxX = 0.0;
    float minY = swapBuffer.height, maxY = 0.0;

    // test against all vertices
    if (vertex0.position.x < minX) minX = vertex0.position.x;
    if (vertex0.position.x > maxX) maxX = vertex0.position.x;
    if (vertex0.position.y < minY) minY = vertex0.position.y;
    if (vertex0.position.y > maxY) maxY = vertex0.position.y;

    if (vertex1.position.x < minX) minX = vertex1.position.x;
    if (vertex1.position.x > maxX) maxX = vertex1.position.x;
    if (vertex1.position.y < minY) minY = vertex1.position.y;
    if (vertex1.position.y > maxY) maxY = vertex1.position.y;

    if (vertex2.position.x < minX) minX = vertex2.position.x;
    if (vertex2.position.x > maxX) maxX = vertex2.position.x;
    if (vertex2.position.y < minY) minY = vertex2.position.y;
    if (vertex2.position.y > maxY) maxY = vertex2.position.y;

    Cartesian3 v0 = Cartesian3(vertex0.position.x, vertex0.position.y, 0);
    Cartesian3 v1 = Cartesian3(vertex1.position.x, vertex1.position.y, 0);
    Cartesian3 v2 = Cartesian3(vertex2.position.x, vertex2.position.y, 0);

    Cartesian3 v0v1 = v1 - v0;
    Cartesian3 n_v0v1 = Cartesian3(-v0v1.y, v0v1.x, 0);
    Cartesian3 v1v2 = v2 - v1;
    Cartesian3 n_v1v2 = Cartesian3(-v1v2.y, v1v2.x, 0);
    Cartesian3 v2v0 = v0 - v2;
    Cartesian3 n_v2v0 = Cartesian3(-v2v0.y, v2v0.x, 0);

    float dAlpha = LeedsGLUtils::distancePointLine(v0, n_v1v2, v1);
    float dBeta = LeedsGLUtils::distancePointLine(v1, n_v2v0, v2);
    float dGamma = LeedsGLUtils::distancePointLine(v2, n_v0v1, v0);

    if (abs(dAlpha - 0) < std::numeric_limits<float>::epsilon() ||
        abs(dBeta - 0) < std::numeric_limits<float>::epsilon() ||
        abs(dGamma - 0) < std::numeric_limits<float>::epsilon())
        return;

    // create a fragment for reuse
    Fragment rasterFragment;

    // loop through the pixels in the bounding box
    for (rasterFragment.row = int(minY); rasterFragment.row <= maxY; rasterFragment.row++) { // per row
        // this is here so that clipping works correctly
        if (rasterFragment.row < 0) continue;
        if (rasterFragment.row >= int(swapBuffer.height)) continue;
        for (rasterFragment.col = int(minX); rasterFragment.col <= maxX; rasterFragment.col++) { // per pixel
            // this is also for correct clipping
            if (rasterFragment.col < 0) continue;
            if (rasterFragment.col >= int(swapBuffer.width)) continue;

            // the pixel in cartesian format
            Cartesian3 pixel(rasterFragment.col + 0.5f, rasterFragment.row + 0.5f, 0.0f);

            // right - we have a pixel inside the frame buffer AND the bounding box
            // note we *COULD* compute gamma = 1.0 - alpha - beta instead
            float alpha = LeedsGLUtils::distancePointLine(pixel, n_v1v2, v1) / dAlpha;
            float beta = LeedsGLUtils::distancePointLine(pixel, n_v2v0, v2) / dBeta;
            float gamma = LeedsGLUtils::distancePointLine(pixel, n_v0v1, v0) / dGamma;

            // now perform the half-plane test
            if ((alpha < 0.0f) || (beta < 0.0f) || (gamma < 0.0f))
                continue;

            //TODO: use the computed baricentric coordinates to interpolate all of the necessary properties for rendering.
            //Be aware of making sure they are perspective correct.

        } // per pixel
    } // per row
}

void LeedsGL::processFragments(std::vector<Fragment> &fragments) {
    //TODO: Process all of the fragments, shading according to the uniform properties.
    //Depth test should go here. We don't explicitly have a pre or post fragment stage.
    //Consider the "shading" as the fragment stage. Decide if the depth test should go before or after, and justify.
}

RGBAValueF LeedsGL::CalculateLighting(const Homogeneous4 &n_vcs, const Homogeneous4 &v_vcs, const RGBAValueF &em,
                                      const RGBAValueF &am, const RGBAValueF &diff, const RGBAValueF &spec,
                                      float shin) {
    if (!(n_vcs.x == 0.0f && n_vcs.y == 0.0f && n_vcs.z == 0.0f)) // we shouldn't try shading if there are no normals
        return RGBAValueF();

    Cartesian3 lightVector;
    Cartesian3 unitNormal = n_vcs.Vector().unit();
    //Directional Light

    Homogeneous4 lp = lightMatrix * lightPosition;

    if (abs(lp.w - 0) < std::numeric_limits<float>::epsilon())
        lightVector = lp.Vector().unit();
    else //point light
        lightVector = (lp - v_vcs).Vector().unit();
    Cartesian3 eyeVector = perspective ? -1 * v_vcs.Point() : Cartesian3(0, 0, 1);
    Cartesian3 bisector = (lightVector + eyeVector).unit();

    RGBAValueF emissive = em;
    RGBAValueF ambient = am.modulate(ambientMaterial);

    float dDot = unitNormal.dot(lightVector);
    dDot = dDot < 0 ? 0 : dDot;

    RGBAValueF diffuse = dDot * diff.modulate(lightColour);

    float sDot = unitNormal.dot(bisector);
    sDot = sDot < 0 ? 0 : sDot;
    sDot = pow(sDot, shin);
    //sDot = ((f.shininess+2)/8.0f)*sDot*dDot;
    sDot = dDot > 0 ? sDot : 0;
    sDot = sDot * dDot * (shin + 2) / 2 * float(M_PI);

    Cartesian3 fs = Cartesian3(spec.red, spec.green, spec.blue);
    Cartesian3 air = Cartesian3(1, 1, 1);
    Cartesian3 a = (air - fs);
    Cartesian3 b = (air + fs);
    Cartesian3 r0 = Cartesian3(a.x / b.x, a.y / b.y, a.z / b.z);
    r0 = Cartesian3(r0.x * r0.x, r0.y * r0.y, r0.z * r0.z);
    Cartesian3 rschlick = r0 + (air - r0) * powf((1 - bisector.dot(lightVector)), 5);
    RGBAValueF updatedSpecular = RGBAValueF(rschlick.x, rschlick.y, rschlick.z, 1);
    RGBAValueF specular = sDot * updatedSpecular.modulate(lightColour);

    return emissive + ambient + diffuse + specular;
}
