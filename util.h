#pragma once

#include <glm/glm.hpp>
#include <json.hpp>

// Returns a random float between 0 and 1, inclusive
float randFloat();

// Returns a random float between min and max, inclusive
float randFloat(float min, float max);

// Returns a random unit vector
glm::vec3 randUnitVec();

// Returns a vector sampled uniformly from a hemisphere
glm::vec3 randUniformHemisphere();

// Returns a vector sampled uniformly from a sphere
glm::vec3 randUniformSphere();

// Returns a vector sampled from a hemisphere using cosine-weighting
glm::vec3 randCosineHemisphere();

// Returns the transformation which transforms the input vector v to point
// in positive z-axis.  v should be normalized.
glm::vec4 orientToZAxis(const glm::vec3& v);

// Returns the transformation which transforms the positive z-axis to
// the input vector v.  v should be normalized
glm::vec4 orientFromZAxis(const glm::vec3& v);

// Rotates the input vector v by quaternion rotation q;
glm::vec3 rotateVector(const glm::vec4& q, const::glm::vec3& v);

glm::vec3 mat4TimesVec3(const glm::mat4& m, const glm::vec3 v, float w);

// Fresenel reflectance, taken from PBRT
float fresnel(float cos_theta, float eta_i, float eta_t);

// Converts degrees to radians
float deg2rad(float deg);

// Takes the sqrt, while ensuring x is non-negative.
inline float safe_sqrt(float x) {
  return std::sqrt(std::max(0.0f, x));
}

// Computes the luminance of an RGB color.
inline float luminance(const glm::vec3& c) {
  //return 0.212f * c[0] + 0.7152f * c[1] + 0.0722f * c[2];
  return 0.212671f * c[0] + 0.715160f * c[1] + 0.072169f * c[2];
}

// Threshold for 0
extern float EPS;

bool nearZero(const glm::vec3& v);


void print_vec3(const glm::vec3& v);

/* Attempts to parse a vec3 from the given data.
 * Quits program if unable to parse. */
glm::vec3 parse_vec3(const nlohmann::json& data);