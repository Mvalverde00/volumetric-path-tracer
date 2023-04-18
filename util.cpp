#include "util.h"


#if defined (_MSC_VER)  // Visual studio
#define thread_local __declspec( thread )
#elif defined (__GCC__) // GCC
#define thread_local __thread
#endif

#include <random>
#include <iostream>

const float PI = 3.141592f;

std::random_device rd;
std::uniform_real_distribution<float> distribution(0.0, 1.0);
std::mt19937 generator = std::mt19937(rd());


#include <thread>
#include <time.h>

float randFloat() {
  static thread_local std::mt19937* generators = nullptr;
  // TODO: DELETE DEBUG
  //if (!generators) generators = new std::mt19937(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
  if (!generators) generators = new std::mt19937(0);
  std::uniform_real_distribution<float> distribution(0.0, 1.0);
  return distribution(*generators);


  return distribution(generator);
}

float randFloat(float min, float max) {
  return min + (randFloat() * (max - min));
}

glm::vec3 randUnitVec() {
  float lambda = acosf(2.0f * randFloat() - 1.0f) - 3.141592653f / 2.0f;
  float phi = 3.141592653f * 2.0f * randFloat();
  float cosl = cos(lambda);
  return glm::vec3(cosl * cosf(phi), cosl * sinf(phi), sinf(lambda));
};


glm::vec3 randUniformHemisphere() {
  float u = randFloat();
  float v = randFloat();

  float z = u;
  float r = std::sqrt(std::max(0.0, 1.0 - z * z));
  float phi = 2.0f * 3.141592653f * v;

  return glm::vec3(r * std::cos(phi), r * std::sin(phi), z);
}

glm::vec3 randUniformSphere() {
  float u = randFloat();
  float v = randFloat();

  float z = 1.0 - 2.0 * u;
  float r = std::sqrt(std::max(0.0, 1.0 - z * z));
  float phi = 2.0f * 3.141592653f * v;

  return glm::vec3(r * std::cos(phi), r * std::sin(phi), z);
}

glm::vec2 randConcentricDisk() {
  float u = 2.0 * randFloat() - 1.0;
  float v = 2.0 * randFloat() - 1.0;

  if (u == 0 && v == 0) return glm::vec2(0.0, 0.0);

  float theta, r;
  if (std::abs(u) > std::abs(v)) {
    r = u;
    theta = (3.141592653f / 4.0) * (v / u);
  } else {
    r = v;
    theta = (3.141592653f / 2.0) - (3.141592653f / 4.0) * (u / v);
  }

  return r * glm::vec2(std::cos(theta), std::sin(theta));
}

glm::vec3 randCosineHemisphere() {
  glm::vec2 d = randConcentricDisk();

  float z = std::sqrt(std::max(0.0, 1.0 - d.x * d.x - d.y * d.y));
  return glm::vec3(d, z);
}

glm::vec4 orientToZAxis(const glm::vec3& v) {
  if (v.z < -0.999999f) return glm::vec4(1.0, 0.0, 0.0, 0.0);

  return glm::normalize(glm::vec4(v.y, -v.x, 0.0, 1.0 + v.z));
}

glm::vec4 orientFromZAxis(const glm::vec3& v) {
  if (v.z < -0.999999f) return glm::vec4(1.0, 0.0, 0.0, 0.0);

  return glm::normalize(glm::vec4(-v.y, v.x, 0.0, 1.0 + v.z));
}

glm::vec3 rotateVector(const glm::vec4& q, const::glm::vec3& v) {
  const glm::vec3 q_axis = glm::vec3(q);
  return 2.0f * glm::dot(q_axis, v) * q_axis
         + (q.w * q.w - glm::dot(q_axis, q_axis)) * v
         + 2.0f * q.w * glm::cross(q_axis, v);
}

glm::vec3 mat4TimesVec3(const glm::mat4& m, const glm::vec3 v, float w) {
  return glm::vec3(m * glm::vec4(v, w));
}

float fresnel(float cos_theta, float eta_i, float eta_t) {
  cos_theta = glm::clamp(cos_theta, -1.0f, 1.0f);

  if (cos_theta <= 0.0f) {
    std::swap(eta_i, eta_t);
    cos_theta = -cos_theta;
  }

  float sin_theta_i = std::sqrt(std::max(0.0f, 1 - cos_theta * cos_theta));
  float sin_theta_t = eta_i / eta_t * sin_theta_i;
  // Total internal reflection
  if (sin_theta_t >= 1.0f) return 1.0f;

  float cos_theta_t = std::sqrt(std::max(0.0f, 1 - sin_theta_t * sin_theta_t));

  float parl = ((eta_t * cos_theta) - (eta_i * cos_theta_t)) /
               ((eta_t * cos_theta) + (eta_i * cos_theta_t));
  float perp = ((eta_i * cos_theta) - (eta_t * cos_theta_t)) /
               ((eta_i * cos_theta) + (eta_t * cos_theta_t));

  return (parl * parl + perp * perp) / 2.0f;
}

float deg2rad(float deg) {
  return deg * PI / 180.0f;
}


float EPS = 0.00001f;


bool nearZero(const glm::vec3& v) {
  return fabs(v.x) < EPS && fabs(v.y) < EPS && fabs(v.z) < EPS;
}

void print_vec3(const glm::vec3& v) {
  std::cout << v.x << ", " << v.y << ", " << v.z << "\n";
}


glm::vec3 parse_vec3(const nlohmann::json& data) {
  if (data.type() == nlohmann::json::value_t::array) {
    return glm::vec3(data[0], data[1], data[2]);
  }
  else if (data.type() == nlohmann::json::value_t::number_float ||
    data.type() == nlohmann::json::value_t::number_integer ||
    data.type() == nlohmann::json::value_t::number_unsigned) {
    return glm::vec3(float(data));
  }

  std::cout << "Error parsing vec3: '" << data << "'.\n";
  exit(1);
}