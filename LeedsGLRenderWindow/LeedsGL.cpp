#include "LeedsGL.h"
#include <math.h>

struct TicTok {
    const char *label;
    std::chrono::steady_clock::time_point start;

    explicit TicTok(const char *l) : label(l) {
        start = std::chrono::steady_clock::now();
    }

    ~TicTok() {
        std::cout << "[" << label << "] cost : "
                  << (std::chrono::steady_clock::now() - start).count() * 1e-6
                  << std::endl;

    }
};

using namespace TP;
using std::byte;
using std::string;

LeedsGL::LeedsGL() {
#ifdef DEBUG
    this->pool = new ThreadPool(1);
#else
    this->pool = new ThreadPool(std::thread::hardware_concurrency());
#endif
    this->texture = nullptr;
}

LeedsGL::~LeedsGL() {
    // free my thread pool
    delete this->pool;
    // texture is hold by Texture object no need to free here
}

void LeedsGL::clear(byte mask) {//Look at the usage of this function on LeedsGLRenderWidget to
    //know what to expect of the "mask" parameter. Each bit should
    //tell which buffer to clear.

    if ((mask & COLORMASK) != UNKNOWN_MASK) {
        // clear swapBuffer
        for (int i = 0; i < swapBuffer.height; i++) {
            for (int j = 0; j < swapBuffer.width; j++) {
                swapBuffer[i][j] = bufferClearColor;
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
    if ((mask & DEPTHMASK) != UNKNOWN_MASK) {
        // clear depth buffer
        for (auto &depth: depthBuffer) {
            depth = std::numeric_limits<float>::infinity();
        }
    }
}

void LeedsGL::setUniform(const string &name, const bool value) {
    //uniform variables should decide how your shading happens
    if (name == "lightingEnabled") {
        this->lightingEnabled = value;
    } else if (name == "texturingEnabled") {
        this->texturingEnabled = value;
    } else if (name == "textureModulationEnabled") {
        this->textureModulationEnabled = value;
    } else if (name == "UVColourDebug") {
        this->UVColourDebug = value;
    }
}

void LeedsGL::setUniform(const string &name, const Matrix4 &mat) {
    //uniform variables should decide how your shading happens
    if (name == "viewportMatrix") {
        this->viewPortMatrix = mat;
    } else if (name == "projectionMatrix") {
        this->projectionMatrix = mat;
    } else if (name == "modelviewMatrix") {
        this->modelviewMatrix = mat;
    } else if (name == "lightMatrix") {
        this->lightMatrix = mat;
    }
}

void LeedsGL::setUniform(const string &name, const RGBAValueF &col) {
    //uniform variables should decide how your shading happens
    if (name == "lightColour") {
        this->lightColour = col;
    } else if (name == "emissiveMaterial") {
        this->emissiveMaterial = col;
    } else if (name == "ambientMaterial") {
        this->ambientMaterial = col;
    } else if (name == "diffuseMaterial") {
        this->diffuseMaterial = col;
    } else if (name == "specularMaterial") {
        this->specularMaterial = col;
    }
}

void LeedsGL::setUniform(const std::string &name, const float val) {
    //uniform variables should decide how your shading happens
    if (name == "shininessMaterial") {
        this->shininessMaterial = val;
    }
}

void LeedsGL::setUniform(const std::string &name, const Homogeneous4 &pos) {
    //uniform variables should decide how your shading happens
    if (name == "lightPosition") {
        this->lightPosition = pos;
    }
}

void LeedsGL::clearColor(const RGBAValueF &col) {
    // update clear color
    bufferClearColor = col;
}

void LeedsGL::resizeBuffers(unsigned const int width, unsigned const int height) {
#ifdef __APPLE__
    // pixel scale on macos
    swapBuffer.Resize(width / 2, height / 2);
#else
    swapBuffer.Resize(width, height);
#endif
    frameBuffer.Resize(width, height);
    // resize depth buffer
    // set every value as infinity
    depthBuffer.resize(width * height, std::numeric_limits<float>::infinity());
}

Matrix4 LeedsGLUtils::calculateViewportMatrix(float cx, float cy, float width, float height) {
    //performs the viewport transformation. NDCS -> DCS.
#ifdef __APPLE__
    width = width / 2;
    height = height / 2;
#endif
    Matrix4 vp;
    vp[0][0] = width / 2, vp[0][1] = 0, vp[0][2] = 0, vp[0][3] = (cx + 1) * width / 2;
    vp[1][0] = 0, vp[1][1] = height / 2, vp[1][2] = 0, vp[1][3] = (cy + 1) * height / 2;
    vp[2][0] = 0, vp[2][1] = 0, vp[2][2] = 1, vp[2][3] = 0;
    vp[3][0] = 0, vp[3][1] = 0, vp[3][2] = 0, vp[3][3] = 1;
    return vp;
}

Matrix4
LeedsGLUtils::calculateProjectionOrtho(float left, float right, float bottom, float top, float near, float far) {
    //right or left-handedness may have effects on other parts of your code,
    //such as shading, clipping, culling, etc.
    Matrix4 po;
    po[0][0] = 2 / (right - left), po[0][1] = 0, po[0][2] = 0, po[0][3] = -(right + left) / (right - left);
    po[1][0] = 0, po[1][1] = 2 / (top - bottom), po[1][2] = 0, po[1][3] = -(top + bottom) / (top - bottom);
    po[2][0] = 0, po[2][1] = 0, po[2][2] = 1 / (far - near), po[2][3] = -near / (far - near);
    po[3][0] = 0, po[3][1] = 0, po[3][2] = 0, po[3][3] = 1;
    return po;
}

Matrix4
LeedsGLUtils::calculateProjectionFrustum(float left, float right, float bottom, float top, float near, float far) {
    //right or left-handedness may have effects on other parts of your code,
    //such as shading, clipping, culling, etc.
    Matrix4 fr;
    fr[0][0] = 2 * near / (right - left), fr[0][1] = 0, fr[0][2] = -(right + left) / (right - left), fr[0][3] = 0;
    fr[1][0] = 0, fr[1][1] = 2 * near / (top - bottom), fr[1][2] = -(top + bottom) / (top - bottom), fr[1][3] = 0;
    fr[2][0] = 0, fr[2][1] = 0, fr[2][2] = far / (far - near), fr[2][3] = -far * near / (far - near);
    fr[3][0] = 0, fr[3][1] = 0, fr[3][2] = 1, fr[3][3] = 0;
    return fr;
}

void LeedsGL::texImage2D(RGBAImage const *textureImage) {
    //Parameter is a pointer to the texture, be aware of how it is stored, and ownership of the resources.
    texture = textureImage;
}

void LeedsGL::enable(const std::byte function) {
    if ((function & PERSPECTIVE) != UNKNOWN_MASK) {
        perspective = true;
    } else if ((function & DEPTHTEST) != UNKNOWN_MASK) {
        depthTestEnabled = true;
    }
}

void LeedsGL::disable(const std::byte function) {
    if ((function & PERSPECTIVE) != UNKNOWN_MASK) {
        perspective = false;
    } else if ((function & DEPTHTEST) != UNKNOWN_MASK) {
        depthTestEnabled = false;
    }
}

void LeedsGL::lineWidth(const float width) {
    rasterizedLineWidth = width;
}

void LeedsGL::pointSize(const float size) {
    rasterizedPointSize = size;
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
    for (int i = 0; i < swapBuffer.height; i++) {
        for (int j = 0; j < swapBuffer.width; j++) {
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
    //This function should combine this disjoint information into a series of InputVertex to be processed
    //by the next step.
    // build up  input vertices
    TicTok t("inputAssembly");
    result.resize(vertices.size());
    TaskGroup tasks;
    for (unsigned int i = 0; i < vertices.size(); i++) {
        tasks.emplace_back([i, &vertices, &normals,
                                   &textureCoordinates, &colors, &result]() {
            result[i].position = vertices[i];
            if (i < normals.size()) {
                result[i].normal = normals[i];
                result[i].normal.w = 0;
            }
            if (i < textureCoordinates.size())
                result[i].tex_coord = textureCoordinates[i];
            if (i < colors.size())
                result[i].color = colors[i];
        });
    }
    this->pool->syncGroup(tasks);
}

void LeedsGL::transformVertices(std::vector<InputVertex> &vertices, std::vector<TransformedVertex> &result) {
    //Also pass all the necessary information to the next steps of the pipeline.
    //You should check the slides to decide which is the appropriate coordinate system to transform them to.
    // transform all position into NDCS, since camera is always on center
    TicTok t("transformVertices");
    Matrix4 trans;
    // translate the CS into right hand
    trans.SetScale(1, 1, -1);
    Matrix4 mvp = projectionMatrix * trans * modelviewMatrix;
    Matrix4 mv = trans * modelviewMatrix;
    TaskGroup tasks;
    result.resize(vertices.size());
    for (unsigned int i = 0; i < vertices.size(); i++) {
        tasks.emplace_back([i, &vertices, &result, &mvp, &mv]() {
            // Transform position and normal
            result[i].position = mvp * vertices[i].position;
            result[i].w = result[i].position.w;
            result[i].v_vcs = mv * vertices[i].position;
            result[i].normal = (mv * vertices[i].normal).Vector();
            result[i].color = vertices[i].color;
            result[i].tex_coord = vertices[i].tex_coord;
        });
    }
    this->pool->syncGroup(tasks);
}

void
LeedsGL::primitiveAssembly(std::vector<TransformedVertex> &vertices, std::byte mode, std::vector<Primitive> &result) {
    TicTok t("primitiveAssembly");
    int strip_step = int(mode) - 4;
    TaskGroup tasks;
    result.clear();
    result.resize(vertices.size() / strip_step);
    for (int i = 0; i < vertices.size() / strip_step; i++) {
        tasks.emplace_back([i, strip_step, &result, &vertices]() {
            result[i].transformedVertices.clear();
            for (int j = 0; j < strip_step; j++) {
                result[i].transformedVertices.push_back(vertices[i * strip_step + j]);
            }
        });
    }
    this->pool->syncGroup(tasks);
}

void LeedsGL::clipAndCull(std::vector<Primitive> &primitives, std::byte mode, std::vector<Primitive> &result) {
    //Pay attention to what type of projection you are using, as your clipping planes will be different.
    //If you choose to skip this step as it is one of your last tasks, just return all the same primitives.
    const byte INSIDE{0}, LEFT{1}, RIGHT{2}, BOTTOM{4}, TOP{8}, FRONT{16}, BACK{32};
    TicTok t("clipAndCull");
    // clip utils:
    // F1. get intersection info
    const auto is_in_NDCS = [](Homogeneous4 &p) {
        std::byte res = INSIDE;
        if (p.w < 0) {
            return FRONT;
        }
        p = p / p.w;
        if (p.z < 0) res |= FRONT; else if (p.z > 1) res |= BACK;
        if (p.x < -1) res |= LEFT; else if (p.x > 1) res |= RIGHT;
        if (p.y < -1) res |= BOTTOM; else if (p.y > 1) res |= TOP;
        return res;
    };
    // convert position to DCS
    const auto transform = [this](Homogeneous4 &position) {
        position = viewPortMatrix * position;
        position = position / position.w;
    };
    // cut a line with a given surface
    // return an interpaloted vertex
    const auto cutLine =
            [](const TransformedVertex &v0, const TransformedVertex &v1, byte mode) {
                float alpha = 0;
                if (mode == FRONT) {
                    alpha = (0 - v1.position.z) / (v0.position.z - v1.position.z);
                } else if (mode == BACK) {
                    alpha = (1 - v1.position.z) / (v0.position.z - v1.position.z);
                } else if (mode == LEFT) {
                    alpha = (-1 - v1.position.x) / (v0.position.x - v1.position.x);
                } else if (mode == RIGHT) {
                    alpha = (1 - v1.position.x) / (v0.position.x - v1.position.x);
                } else if (mode == BOTTOM) {
                    alpha = (-1 - v1.position.y) / (v0.position.y - v1.position.y);
                } else if (mode == TOP) {
                    alpha = (1 - v1.position.y) / (v0.position.y - v1.position.y);
                }
                float beta = 1 - alpha;
                TransformedVertex res;
                res.position = alpha * v0.position + beta * v1.position;
                res.v_vcs = alpha * v0.v_vcs + beta * v1.v_vcs;
                res.color = alpha * v0.color + beta * v1.color;
                res.normal = alpha * v0.normal + beta * v1.normal;
                res.tex_coord = alpha * v0.tex_coord + beta * v1.tex_coord;
                res.w = res.position.w;
                return res;
            };
    // clear previous data;
    result.clear();
    if (mode == POINTS) {
        for (auto &p: primitives) {
            // check if the point is on the back of the camera
            if (is_in_NDCS(p.transformedVertices[0].position) == INSIDE) {
                // transform to DCS
                transform(p.transformedVertices[0].position);
                result.push_back(p);
            }
        }
    } else if (mode == LINES) {
        for (auto &line: primitives) {
            Homogeneous4 &p1 = line.transformedVertices[0].position;
            Homogeneous4 &p2 = line.transformedVertices[1].position;
            auto res1 = is_in_NDCS(p1);
            auto res2 = is_in_NDCS(p2);
            if (res1 == INSIDE && res2 == INSIDE) {
                // do nothing, just draw it
            } else if ((res1 & res2) != INSIDE) {
                // the line is completely outside, throw it
                continue;
            } else {
                // part inside
                // for each side
                auto v1 =  line.transformedVertices[0];
                auto v2 =  line.transformedVertices[1];
                for (int i = 0; i < 6; i ++) {
                    // calculate side
                    byte which_side = byte(1) << i;
                    if ((res1 & which_side) == INSIDE && (res2 & which_side) == INSIDE) {
                        // both are on same side of this line
                        // no need to cut
                        continue;
                    }
                    if ((res1 & which_side) != INSIDE) {
                        // v2 is in side, swap the two point
                        std::swap(v1, v2);
                    }
                    // v1 inside, v2 outside
                    // get the cut point, update v2
                    v2 = cutLine(v1, v2, which_side);
                    // update side info
                    res1 = is_in_NDCS(v1.position);
                    res2 = is_in_NDCS(v2.position);
                }
                // udpate vertex on line
                line.transformedVertices[0] = v1;
                line.transformedVertices[1] = v2;
            }
            // draw the clipped line
            // transform it into DCS
            transform(line.transformedVertices[0].position);
            transform(line.transformedVertices[1].position);
            result.push_back(line);

        }
    } else if (mode == TRIANGLES) {
        TaskGroup tasks;
        // avoid lock
        std::vector<std::vector<Primitive>> temp_res(primitives.size());
        for (int pri = 0; pri < primitives.size(); pri ++) {

            tasks.emplace_back([pri, &temp_res, is_in_NDCS, cutLine, &primitives, &transform](){
                // clip a triangle:
                // get the position info of vertices
                auto tri = primitives[pri];
                std::vector<TransformedVertex> clipped_vertices;
                // current_vertices contains a polygon from previous iteration
                std::vector<TransformedVertex> current_vertices = tri.transformedVertices;
                for (int i = 0; i < 6; i ++) {
                    std::vector<byte> side_infos(current_vertices.size());
                    // calculate side
                    byte which_side = byte(1) << i;
                    // get location information on this side
                    for (int j = 0; j < side_infos.size(); j ++) {
                        side_infos[j] = is_in_NDCS(current_vertices[j].position);
                    }
                    // analise each line
                    for (int j = 0; j < side_infos.size();j ++) {
                        int begin = j;
                        int end = (j + 1) % side_infos.size();
                        auto begin_v = current_vertices[begin];
                        auto end_v = current_vertices[end];
                        bool begin_inside = (side_infos[begin] & which_side) == INSIDE;
                        bool end_inside = (side_infos[end] & which_side) == INSIDE;
                        if (begin_inside && end_inside) {
                            // both are inside, add to clipped_vertices
                            clipped_vertices.push_back(begin_v);
                            // push the first one, the next one will be added in next iteration
                            // clipped_vertices.push_back(end_v);
                        } else if (!begin_inside && !end_inside) {
                            // both are outside, throw them
                            continue;
                        } else {
                            // one inside, one out side
                            auto intersection = cutLine(begin_v, end_v, which_side);
                            // add points in order
                            if (begin_inside) {
                                clipped_vertices.push_back(begin_v);
                                clipped_vertices.push_back(intersection);
                            } else {
                                clipped_vertices.push_back(intersection);
                                // the next one will be added in next iteration
                                // clipped_vertices.push_back(end_v);
                            }
                        }
                    }
                    // cut finished on this side
                    // update current infos
                    current_vertices = clipped_vertices;
                    clipped_vertices.clear();
                }
                for (auto &v : current_vertices) {
                    // transform position into DCS
                    transform(v.position);
                }
                // generate triangles from current_vertices;
                // 0, 1, 2
                // 0, 2, 3
                // ...
                // 0, n-1, n
                // lock result
                for (int j = 2; j < current_vertices.size(); j ++) {
                    Primitive new_tri;
                    new_tri.transformedVertices = {current_vertices[0],current_vertices[j-1], current_vertices[j]};
                    temp_res[pri].push_back(new_tri);
                }
            });
        }
        pool->syncGroup(tasks);
        for (auto &res:temp_res) {
            for (auto &tri : res) {
                result.push_back(tri);
            }
        }
    }
}

void LeedsGL::rasterisePrimitives(std::vector<Primitive> &primitives, std::byte mode, std::vector<Fragment> &results) {
    TicTok t("rasterisePrimitives");
    TaskGroup tasks;
    // reset result buffers;
    results.clear();
    // resize it in to primitives size
    results.resize(primitives.size());
    for (int i = 0; i < primitives.size(); i++) {
        if (mode == POINTS) {
            rasterisePoint(i, primitives[i], results);
        } else if (mode == LINES) {
            rasteriseLine(i, primitives[i], results);
        } else if (mode == TRIANGLES) {
            tasks.emplace_back([i, this, &primitives, &results]() {
                rasteriseTriangle(i, primitives[i], results);
            });
        }
    }
    pool->syncGroup(tasks);
}

void LeedsGL::rasterisePoint(int index, const Primitive &point, std::vector<Fragment> &output) {
    // get vertex
    auto vertex = point.transformedVertices[0];
    // a function to clip coords
    const auto clip = [](int in, int min, int max) {
        return std::max(min, std::min(in, max));
    };
    // calculate the fragment size according to the point size * position
    int left = clip(int(vertex.position.x - rasterizedPointSize / 2), 0, int(swapBuffer.width) - 1);
    int right = clip(int(vertex.position.x + rasterizedPointSize / 2), 0, int(swapBuffer.width) - 1);
    int top = clip(int(vertex.position.y + rasterizedPointSize / 2), 0, int(swapBuffer.height) - 1);
    int bottom = clip(int(vertex.position.y - rasterizedPointSize / 2), 0, int(swapBuffer.height) - 1);
    output[index].row = bottom;
    output[index].col = left;
    output[index].width = right - left + 1;
    output[index].height = top - bottom + 1;
    output[index].resize();
    RGBAValueF color = calculateColor(vertex.color, Cartesian3(), vertex.normal, vertex.v_vcs);
    for (int x = 0; x < output[index].width; x++) {
        for (int y = 0; y < output[index].height; y++) {
            output[index].colors[y * output[index].width + x] = color;
            output[index].depths[y * output[index].width + x] = vertex.position.z;
        }
    }
}

void LeedsGL::rasteriseLine(int index, const Primitive &line, std::vector<Fragment> &output) {
    // treat the line as a quad and change it into triangles
    // then rasterise Triangles
    TransformedVertex v0 = line.transformedVertices[0];
    TransformedVertex v1 = line.transformedVertices[1];
    Homogeneous4 direction = v0.position - v1.position;
    Cartesian3 normal(direction.y, -direction.x, 0);
    // use the padding to get a rect
    Cartesian3 padding = rasterizedLineWidth * normal.unit() / 2;
    // calculate the four position of the rect
    TransformedVertex start_0, start_1, end_0, end_1;
    start_0.position = v0.position + padding;
    start_0.normal = v0.normal;
    start_0.w = v0.w;
    start_0.color = v0.color;
    start_0.v_vcs = v0.v_vcs;
    start_0.tex_coord = v0.tex_coord;

    start_1.position = v0.position - padding;
    start_1.normal = v0.normal;
    start_1.w = v0.w;
    start_1.color = v0.color;
    start_1.v_vcs = v0.v_vcs;
    start_1.tex_coord = v0.tex_coord;

    end_0.position = v1.position + padding;
    end_0.normal = v1.normal;
    end_0.w = v1.w;
    end_0.color = v1.color;
    end_0.v_vcs = v1.v_vcs;
    end_0.tex_coord = v1.tex_coord;

    end_1.position = v1.position - padding;
    end_1.normal = v1.normal;
    end_1.w = v1.w;
    end_1.color = v1.color;
    end_1.v_vcs = v1.v_vcs;
    end_1.tex_coord = v1.tex_coord;

    // get result of two triangles
    int size = output.size();
    // append one fragment to draw an extra triangle
    output.resize(size + 1);
    // draw one at index
    rasteriseTriangle(index, {{start_1, start_0, end_0}}, output);
    // draw another at the end of output
    rasteriseTriangle(size, {{end_0, end_1, start_1}}, output);
}


float LeedsGLUtils::distancePointLine(Cartesian3 r, Cartesian3 n, Cartesian3 p) {
    //assumes n is normalized
    return n.dot(r) - n.dot(p);
}

void LeedsGL::rasteriseTriangle(int index, const Primitive &triangle, std::vector<Fragment> &output) {
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
    // a function to clip coords
    const auto clip = [](int in, int min, int max) {
        return std::max(min, std::min(in, max));
    };
    // create a fragment for reuse
    Fragment rasterFragment;
    int left = clip(int(minX), 0, int(swapBuffer.width) - 1);
    int right = clip(int(maxX), 0, int(swapBuffer.width) - 1);
    int top = clip(int(maxY), 0, int(swapBuffer.height) - 1);
    int bottom = clip(int(minY), 0, int(swapBuffer.height) - 1);
    // update size of fragment
    rasterFragment.row = bottom;
    rasterFragment.col = left;
    rasterFragment.width = right - left + 1;
    rasterFragment.height = top - bottom + 1;
    // init the fragment buffer
    rasterFragment.resize();
    TP::TaskGroup tasks;
    // loop through the pixels in the bounding box
    for (int y = 0; y < rasterFragment.height; y++) { // per row
        for (int x = 0; x < rasterFragment.width; x++) { // per pixel
            // the pixel in cartesian format
            Cartesian3 pixel(float(rasterFragment.col + x) + 0.5f,
                             float(rasterFragment.row + y) + 0.5f, 0.0f);

            // right - we have a pixel inside the frame buffer AND the bounding box
            // note we *COULD* compute gamma = 1.0 - alpha - beta instead
            float alpha = LeedsGLUtils::distancePointLine(pixel, n_v1v2, v1) / dAlpha;
            float beta = LeedsGLUtils::distancePointLine(pixel, n_v2v0, v2) / dBeta;
            float gamma = LeedsGLUtils::distancePointLine(pixel, n_v0v1, v0) / dGamma;

            // now perform the half-plane test
            if ((alpha < 0.0f) || (beta < 0.0f) || (gamma < 0.0f))
                continue;

            //Be aware of making sure they are perspective correct.
            // interpolation normal_vcs, v_vcs, tex_coords, depth and color (if there is no texture)
            float W = alpha / vertex0.w + beta / vertex1.w + gamma / vertex2.w;
            float depth = alpha * vertex0.position.z + beta * vertex1.position.z + gamma * vertex2.position.z;
            // update alpha beta gamma, Hyperbolic interpolation
            alpha = alpha / vertex0.w, beta = beta / vertex1.w, gamma = gamma / vertex2.w;
            Homogeneous4 n_vcs = (vertex0.normal * alpha + vertex1.normal * beta + vertex2.normal * gamma) / W;
            Homogeneous4 v_vcs = (vertex0.v_vcs * alpha + vertex1.v_vcs * beta + vertex2.v_vcs * gamma) / W;
            Cartesian3 tex_coords =
                    (vertex0.tex_coord * alpha + vertex1.tex_coord * beta + vertex2.tex_coord * gamma) / W;
            RGBAValueF color = (1.f / W) * (alpha * vertex0.color + beta * vertex1.color + gamma * vertex2.color);
            // get color at this pixel
            rasterFragment.colors[y * rasterFragment.width + x] = calculateColor(color, tex_coords, n_vcs, v_vcs);
            rasterFragment.depths[y * rasterFragment.width + x] = depth;
        } // per pixel
    } // per row
    output[index] = rasterFragment;
}

RGBAValueF LeedsGL::textureSampler(const Cartesian3 &uv) {
    if (this->texture == nullptr) {
        return RGBAValueF();
    }
    // transform to the image coords
    float u = std::min(1.f, std::max(0.f, uv.x));
    float v = std::min(1.f, std::max(0.f, uv.y));
    int x = int(u * float(this->texture->width));
    int y = int(v * float(this->texture->height));

    RGBAValue color = (*this->texture)[y][x];
    return {float(color.red) / 255.f,
            float(color.green) / 255.f,
            float(color.blue) / 255.f,
            1.f};
}

void LeedsGL::processFragments(std::vector<Fragment> &fragments) {
    TicTok t("processFragments");
    //Depth test should go here. We don't explicitly have a pre or post fragment stage.
    //Consider the "shading" as the fragment stage. Decide if the depth test should go before or after, and justify.
    for (auto &frag: fragments) {
        for (int x = 0; x < frag.width; x++) {
            for (int y = 0; y < frag.height; y++) {
                // update color from fragment
                // if alpha is 0, there should be no color in frag
                if (frag.colors[y * frag.width + x].alpha > 0) {
                    int px = frag.col + x, py = frag.row + y;
                    int depth_index = py * int(swapBuffer.width) + px;
                    int frag_pixel_index = y * frag.width + x;
                    // depth test
                    if (depthTestEnabled) {
                        if (depthBuffer[depth_index] > frag.depths[frag_pixel_index]) {
                            // this pixel is more closed to the camera
                            // update depth
                            depthBuffer[depth_index] = frag.depths[frag_pixel_index];
                            // update color
                            swapBuffer[py][px] = frag.colors[frag_pixel_index];
                        }
                    } else {
                        // just update color
                        swapBuffer[py][px] = frag.colors[frag_pixel_index];
                    }

                } // set pixel color from frag
            }
        }
    }
}

RGBAValueF LeedsGL::calculateColor(const RGBAValueF &color,
                                   const Cartesian3 &uv,
                                   const Homogeneous4 &n_vcs,
                                   const Homogeneous4 &v_vcs) {
    RGBAValueF res;
    if (texturingEnabled) {
        // get texture color
        res = textureSampler(uv);
        if (textureModulationEnabled && lightingEnabled) {
            // calculate light color
            RGBAValueF light_color = CalculateLighting(n_vcs, v_vcs, emissiveMaterial,
                                                       ambientMaterial,
                                                       diffuseMaterial,
                                                       specularMaterial,
                                                       shininessMaterial);
            // modulate tex_color and light color
            res = res.modulate(light_color);
        }
    } else {
        if (lightingEnabled && !(n_vcs.x == 0.0f && n_vcs.y == 0.0f && n_vcs.z == 0.0f)) {
            // only use light color, if texture is disabled
            res = CalculateLighting(n_vcs, v_vcs, emissiveMaterial,
                                    ambientMaterial,
                                    diffuseMaterial,
                                    specularMaterial,
                                    shininessMaterial);
        } else {
            // use the given color
            res = color;
        }
    }
    return res;
}

RGBAValueF LeedsGL::CalculateLighting(const Homogeneous4 &n_vcs, const Homogeneous4 &v_vcs, const RGBAValueF &em,
                                      const RGBAValueF &am, const RGBAValueF &diff, const RGBAValueF &spec,
                                      float shin) {
    if (n_vcs.x == 0.0f && n_vcs.y == 0.0f && n_vcs.z == 0.0f) // we shouldn't try shading if there are no normals
        return RGBAValueF();

    Cartesian3 lightVector;
    Cartesian3 unitNormal = n_vcs.Vector().unit();
    //Directional Light

    Matrix4 trans;
    trans.SetScale(-1, -1, 1); // transform it into right hand
    Homogeneous4 lp = (trans * lightMatrix) * lightPosition;

    if (abs(lp.w - 0) < std::numeric_limits<float>::epsilon())
        lightVector = lp.Vector().unit();
    else //point light
        lightVector = (lp - v_vcs).Vector().unit();
    Cartesian3 eyeVector = perspective ? -1 * v_vcs.Point() : Cartesian3(0, 0, -1);
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
