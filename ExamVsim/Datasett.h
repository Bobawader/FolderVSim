#pragma once

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <glm/glm.hpp>


class Datasett {
public:
    
    std::vector<glm::vec3> points;       // Original 3D points
    std::vector<size_t> indices;         // Triangle indices
    std::vector<double> delaunayCoords;  // Flattened 2D points for Delaunator
    std::vector<glm::vec3> normals;      // Per-vertex normals
    GLuint VBO, VAO, EBO, normalVBO;     // OpenGL buffers

    // Member functions
    void loadFile(const std::string& filename);                 
    void triangulate();                                         
    void calculateNormals();                                    
    
    //Kalkulasjon
    glm::vec3 computeBarycentric(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) const;
    glm::vec3 getSurfaceNormalAt(const glm::vec3& position) const;
    glm::vec3 calculateClosestPointOnTriangle(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) const;
    float calculateHeightAtPoint(const glm::vec3& point, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) const;


    
    float getHeightAt(const glm::vec3& position) const;
    void setupBuffers(const std::vector<glm::vec3>& points,const std::vector<size_t>& triangles);

    
    //friksjon
    bool isInFrictionZone(size_t triangleIndex) const;
    void assignFrictionCoefficients();
    float getFrictionCoefficientAt(const glm::vec3& position) const;

private:
    std::vector<float> frictionCoefficients;

};






//#pragma once
//#include <glad/glad.h>
//#include <GLFW/glfw3.h>
//#include <iostream>
//#include <vector>
//#include <fstream>
//#include <sstream>
//#include <glm/glm.hpp>
//#include <algorithm>
//#include <limits>
//#include <delaunator.hpp>
//
//
//
//
//
//
//class Datasett {
//public:
//    std::vector<glm::vec3> points;  // Original 3D points
//    std::vector<float> gridVertices;
//    std::vector<size_t> indices;
//    std::vector<glm::vec3> normals;       
//    std::vector<double> delaunayCoords;  // Flattened 2D points for Delaunator
//    GLuint VBO, VAO, EBO;
//
//
//    void loadFile(const std::string& filename) {
//        std::ifstream file(filename);
//        if (!file) {
//            std::cerr << "Failed to open file: " << filename << std::endl;
//            return;
//        }
//
//        std::string line;
//        int numPoints;
//        std::getline(file, line);
//        std::istringstream(line) >> numPoints;
//
//        while (std::getline(file, line)) {
//            std::istringstream ss(line);
//            float x, y, z;
//            ss >> x >> y >> z;
//
//            points.emplace_back(x, z, y); // Swapping z and y if necessary
//
//            // Flatten the 2D points for Delaunator (x, z only)
//            delaunayCoords.push_back(x);
//            delaunayCoords.push_back(z);
//        }
//        file.close();
//
//        std::cout << "Loaded " << points.size() << " points successfully." << std::endl;
//    }
//
//
//    void triangulate() {
//        if (delaunayCoords.empty()) {
//            std::cerr << "No coordinates to triangulate!" << std::endl;
//            return;
//        }
//
//        // Perform Delaunay triangulation
//        delaunator::Delaunator d(delaunayCoords);
//
//        // Store the triangles in the indices member variable
//        indices.clear();
//        indices.insert(indices.end(), d.triangles.begin(), d.triangles.end()); // Directly copy size_t indices
//
//        // Debug output
//        for (std::size_t i = 0; i < indices.size(); i += 3) {
//            std::cout << "Triangle points: [["
//                << d.coords[2 * indices[i]] << ", " << d.coords[2 * indices[i] + 1] << "], ["
//                << d.coords[2 * indices[i + 1]] << ", " << d.coords[2 * indices[i + 1] + 1] << "], ["
//                << d.coords[2 * indices[i + 2]] << ", " << d.coords[2 * indices[i + 2] + 1] << "]]"
//                << std::endl;
//        }
//    }
//
//
//
//
//
//    
//
//    void setupBuffers(const std::vector<glm::vec3>& points, const std::vector<std::size_t>& triangles) {
//        // Create and bind VAO
//        glGenVertexArrays(1, &VAO);
//        glBindVertexArray(VAO);
//
//        // Create VBO for vertex positions
//        glGenBuffers(1, &VBO);
//        glBindBuffer(GL_ARRAY_BUFFER, VBO);
//        glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(glm::vec3), points.data(), GL_STATIC_DRAW);
//
//        // Set vertex attributes for position
//        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
//        glEnableVertexAttribArray(0);
//
//        // Create EBO for triangle indices
//        glGenBuffers(1, &EBO);
//        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
//        glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles.size() * sizeof(std::size_t), triangles.data(), GL_STATIC_DRAW);
//
//        glBindVertexArray(0); // Unbind VAO
//    }
//    //the height data from høydedata downloaded in xyzfile is using Z as height and Y as depth
//    //which is why emplace_back is seemingly doing it in the wrong order now making Z the new Y.
//    //translating into openGL which uses Y as height.
//
//
//
//    //void loadFile(const std::string& filename) {
//    //    std::ifstream file(filename);
//    //    if (!file) {
//    //        std::cerr << "Failed to open file: " << filename << std::endl;
//    //        return;
//    //    }
//
//    //    std::string line;
//    //    int numPoints;
//    //    std::getline(file, line);
//    //    std::istringstream(line) >> numPoints;
//
//    //    // Temporary variables to compute the bounding box
//    //    float minX = std::numeric_limits<float>::max();
//    //    float minY = std::numeric_limits<float>::max();
//    //    float minZ = std::numeric_limits<float>::max();
//    //    float maxX = std::numeric_limits<float>::lowest();
//    //    float maxY = std::numeric_limits<float>::lowest();
//    //    float maxZ = std::numeric_limits<float>::lowest();
//
//    //    // Read points and compute bounding box
//    //    while (std::getline(file, line)) {
//    //        std::istringstream ss(line);
//    //        float x, y, z;
//    //        ss >> x >> y >> z;
//
//    //        // Add the point to the container
//    //        points.emplace_back(x, z, y); // Swapping y and z if necessary
//
//    //        // Update bounding box
//    //        minX = std::min(minX, x);
//    //        minY = std::min(minY, z); // z becomes the new Y (height)
//    //        minZ = std::min(minZ, y);
//    //        maxX = std::max(maxX, x);
//    //        maxY = std::max(maxY, z);
//    //        maxZ = std::max(maxZ, y);
//    //    }
//    //    file.close();
//
//    //    // Calculate center
//    //    glm::vec3 center((minX + maxX) / 2.0f, (minY + maxY) / 2.0f, (minZ + maxZ) / 2.0f);
//
//    //    // Normalize points and center them at the origin
//    //    for (auto& p : points) {
//    //        p -= center;
//    //    }
//
//    //    // Adjusting height so the lowest point is at y = -1
//    //    /*float newMinY = std::numeric_limits<float>::max();
//    //    for (const auto& p : points) {
//    //        newMinY = std::min(newMinY, p.y);
//    //    }
//
//    //    float offsetY = -1.0f - newMinY;
//    //    for (auto& p : points) {
//    //        p.y += offsetY;
//    //    }*/
//
//    //    std::cout << "Loaded " << points.size() << " points successfully." << std::endl;
//    //}
//
//    
//
//
//
//
//    
//
//    void calculateNormals(int numCols, int numRows) {
//        normals.resize(gridVertices.size() / 3, glm::vec3(0.0f, 0.0f, 0.0f));
//
//        // Calculate face normals and accumulate them for each vertex
//        for (size_t i = 0; i < indices.size(); i += 3) {
//            unsigned int idx0 = indices[i];
//            unsigned int idx1 = indices[i + 1];
//            unsigned int idx2 = indices[i + 2];
//
//            glm::vec3 v0 = glm::vec3(gridVertices[idx0 * 3], gridVertices[idx0 * 3 + 1], gridVertices[idx0 * 3 + 2]);
//            glm::vec3 v1 = glm::vec3(gridVertices[idx1 * 3], gridVertices[idx1 * 3 + 1], gridVertices[idx1 * 3 + 2]);
//            glm::vec3 v2 = glm::vec3(gridVertices[idx2 * 3], gridVertices[idx2 * 3 + 1], gridVertices[idx2 * 3 + 2]);
//
//            glm::vec3 edge1 = v1 - v0;
//            glm::vec3 edge2 = v2 - v0;
//            glm::vec3 faceNormal = glm::normalize(glm::cross(edge1, edge2));
//
//            normals[idx0] += faceNormal;
//            normals[idx1] += faceNormal;
//            normals[idx2] += faceNormal;
//        }
//
//        // Normalize the accumulated normals
//        for (auto& normal : normals) {
//            normal = glm::normalize(normal);
//        }
//
//        std::cout << "Calculated normals for " << normals.size() << " vertices." << std::endl;
//    }
//
//    void printNormals() const {
//        for (size_t i = 0; i < normals.size(); ++i) {
//            std::cout << "Normal at vertex " << i << ": ("
//                << normals[i].x << ", "
//                << normals[i].y << ", "
//                << normals[i].z << ")" << std::endl;
//        }
//    }
//
//
//
//};
